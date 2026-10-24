#include "DriftCalibrationFitter.h"

#include "AnalysisConfig.h"

#include <TDirectory.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TNamed.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include <RtypesCore.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

bool prepareDataLine(std::string& line)
{
    const std::size_t comment = line.find('#');
    if (comment != std::string::npos) line.erase(comment);
    return line.find_first_not_of(" \t\r\n") != std::string::npos;
}

std::runtime_error parseError(const std::string& fileName,
                              std::size_t lineNumber,
                              const std::string& message)
{
    return std::runtime_error(fileName + ":" + std::to_string(lineNumber) +
                              ": " + message);
}

std::set<unsigned int> parseDetectorList(const std::string& text,
                                         const std::string& fileName,
                                         std::size_t lineNumber)
{
    std::set<unsigned int> result;
    if (text == "-") return result;
    std::istringstream input(text);
    std::string field;
    while (std::getline(input, field, ',')) {
        if (field.empty()) {
            throw parseError(fileName, lineNumber,
                             "empty detector ID in exclusion list");
        }
        std::size_t used = 0;
        unsigned long id = 0;
        try {
            id = std::stoul(field, &used);
        } catch (const std::exception&) {
            throw parseError(fileName, lineNumber,
                             "invalid detector ID '" + field + "'");
        }
        if (used != field.size() || id >= config::kGermaniumIdBins) {
            throw parseError(fileName, lineNumber,
                             "invalid detector ID '" + field + "'");
        }
        result.insert(static_cast<unsigned int>(id));
    }
    return result;
}

std::string groupName(unsigned int index)
{
    std::ostringstream name;
    name << "Group" << std::setfill('0') << std::setw(6) << index;
    return name.str();
}

std::string spectrumPath(unsigned int groupIndex, unsigned int detectorID)
{
    std::ostringstream path;
    path << groupName(groupIndex) << "/Individual/h1_Ge_E_ID"
         << std::setfill('0') << std::setw(2) << detectorID;
    return path.str();
}

std::string graphName(const std::string& prefix, unsigned int detectorID)
{
    std::ostringstream name;
    name << prefix << "_ID" << std::setfill('0') << std::setw(2)
         << detectorID;
    return name.str();
}

struct FileRecord {
    unsigned int groupIndex = 0;
    std::string series;
    unsigned int run = 0;
    std::string fileName;
};

struct PeakMeasurement {
    double measuredEnergy = 0.0;
    double measuredError = 0.0;
    double referenceEnergy = 0.0;
    bool used = true;
};

struct AcceptedPeakRecord {
    unsigned int groupIndex = 0;
    unsigned int detectorID = 0;
    double referenceEnergy = 0.0;
    double measuredEnergy = 0.0;
};

struct ShapePeakWindowDefinition {
    double energyKeV = 0.0;
    double searchHalfWidthKeV = 10.0;
    double fitHalfWidthKeV = 4.0;
};

struct WindowedShapePeakRecord {
    AcceptedPeakRecord peak;
    double centroidErrorKeV = 0.0;
    double sigmaKeV = 0.0;
    double significance = 0.0;
    unsigned int neighbourCount = 0;
    bool accepted = false;
};

constexpr unsigned int kNonlinearDriftDetectorID = 4;
constexpr unsigned int kQuadraticDriftDetectorID = 16;
constexpr double kShapeLowerEnergyKeV = 301.570;
constexpr double kShapeUpperEnergyKeV = 1000.0;
constexpr double kFirstShapeAnchorKeV = 535.561;
constexpr double kSecondShapeAnchorKeV = 627.540;

bool useForDriftFit(unsigned int detectorID, double referenceEnergy)
{
    // Detector 4's main source calibration has a repeatable interpolation
    // feature between the source anchors.  The 536- and 628-keV peaks must
    // therefore measure that feature, not pull the time-dependent affine fit
    // differently from one group to the next.
    return detectorID != kNonlinearDriftDetectorID ||
           referenceEnergy < 450.0 || referenceEnergy > 750.0;
}

double median(std::vector<double> values)
{
    if (values.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const std::size_t middle = values.size()/2U;
    std::nth_element(values.begin(), values.begin()+middle, values.end());
    const double upper = values[middle];
    if (values.size()%2U != 0U) return upper;
    std::nth_element(values.begin(), values.begin()+middle-1U,
                     values.begin()+middle);
    return 0.5*(values[middle-1U]+upper);
}

struct ShapeCorrection {
    bool valid = false;
    double firstAnchorResidual = 0.0;
    double secondAnchorResidual = 0.0;
    std::array<double, 6> normalizedCoefficients{};
};

enum class CorrectionSource : int {
    Measured = 0,
    Interpolated = 1,
    Nearest = 2,
    Missing = 3
};

struct GroupCorrection {
    bool valid = false;
    CorrectionSource source = CorrectionSource::Missing;
    double offset = 0.0;
    double slope = 1.0;
    double offsetError = 0.0;
    double slopeError = 0.0;
    double rms = 0.0;
    unsigned int peakCount = 0;
};

struct PolynomialCorrection {
    bool valid = false;
    CorrectionSource source = CorrectionSource::Missing;
    std::array<double, 3> coefficients{0.0, 1.0, 0.0};
    double rms = 0.0;
    unsigned int peakCount = 0;
    bool nonlinear = false;
};

struct PiecewiseLinearCorrection {
    bool valid = false;
    std::vector<std::pair<double, double>> anchors;
};

double applyPiecewiseLinear(const PiecewiseLinearCorrection& correction,
                            double energy)
{
    if (!correction.valid || correction.anchors.size() < 2U) return energy;
    std::size_t low = 0;
    while (low+2U < correction.anchors.size() &&
           energy > correction.anchors[low+1U].first) {
        ++low;
    }
    const auto& first = correction.anchors[low];
    const auto& second = correction.anchors[low+1U];
    const double slope = (second.second-first.second)/
                         (second.first-first.first);
    return first.second+slope*(energy-first.first);
}

std::vector<double> robustTimeTrajectory(
    const std::vector<double>& raw,
    const std::vector<std::size_t>& chronologicalGroups)
{
    std::vector<double> result(
        raw.size(), std::numeric_limits<double>::quiet_NaN());
    std::vector<double> all;
    for (double value : raw) if (std::isfinite(value)) all.push_back(value);
    if (all.empty()) return result;
    std::vector<std::size_t> validPositions;
    for (std::size_t position = 0; position < chronologicalGroups.size();
         ++position) {
        if (std::isfinite(raw[chronologicalGroups[position]])) {
            validPositions.push_back(position);
        }
    }
    for (std::size_t position = 0; position < chronologicalGroups.size();
         ++position) {
        std::vector<double> window;
        const std::size_t begin = position > 3U ? position-3U : 0U;
        const std::size_t end = std::min(
            chronologicalGroups.size(), position+4U);
        for (std::size_t neighbour = begin; neighbour < end; ++neighbour) {
            const double value = raw[chronologicalGroups[neighbour]];
            if (std::isfinite(value)) window.push_back(value);
        }
        if (window.empty()) {
            const auto after = std::lower_bound(
                validPositions.begin(), validPositions.end(), position);
            if (after == validPositions.begin()) {
                result[chronologicalGroups[position]] =
                    raw[chronologicalGroups[*after]];
            } else if (after == validPositions.end()) {
                result[chronologicalGroups[position]] =
                    raw[chronologicalGroups[validPositions.back()]];
            } else {
                const std::size_t highPosition = *after;
                const std::size_t lowPosition = *(after-1);
                const double fraction =
                    static_cast<double>(position-lowPosition)/
                    static_cast<double>(highPosition-lowPosition);
                const double low = raw[chronologicalGroups[lowPosition]];
                const double high = raw[chronologicalGroups[highPosition]];
                result[chronologicalGroups[position]] =
                    low+fraction*(high-low);
            }
            continue;
        }
        const double centre = median(window);
        std::vector<double> deviations;
        for (double value : window) {
            deviations.push_back(std::abs(value-centre));
        }
        const double limit = std::max(0.6, 4.0*1.4826*median(deviations));
        std::vector<double> retained;
        for (double value : window) {
            if (std::abs(value-centre) <= limit) retained.push_back(value);
        }
        const double localCentre = retained.empty() ? centre : median(retained);
        const double current = raw[chronologicalGroups[position]];
        result[chronologicalGroups[position]] =
            std::isfinite(current) && std::abs(current-centre) <= limit
                ? 0.65*current+0.35*localCentre : localCentre;
    }
    return result;
}

std::vector<WindowedShapePeakRecord> fitWindowedDetector4ShapePeaks(
    TFile& inputFile,
    const std::array<std::vector<GroupCorrection>,
                     config::kGermaniumIdBins>& corrections,
    const std::vector<FileRecord>& files, std::size_t groupCount,
    const std::vector<ShapePeakWindowDefinition>& definitions)
{
    std::vector<std::string> groupSeries(groupCount);
    std::vector<unsigned int> groupFirstRun(
        groupCount, std::numeric_limits<unsigned int>::max());
    for (const FileRecord& file : files) {
        groupSeries[file.groupIndex] = file.series;
        groupFirstRun[file.groupIndex] = std::min(
            groupFirstRun[file.groupIndex], file.run);
    }
    std::map<std::string, std::vector<std::size_t>> seriesGroups;
    for (std::size_t group = 0; group < groupCount; ++group) {
        seriesGroups[groupSeries[group]].push_back(group);
    }
    for (auto& entry : seriesGroups) {
        std::stable_sort(entry.second.begin(), entry.second.end(),
            [&groupFirstRun](std::size_t left, std::size_t right) {
                return groupFirstRun[left] < groupFirstRun[right];
            });
    }

    std::vector<WindowedShapePeakRecord> records;
    for (const auto& seriesEntry : seriesGroups) {
        const std::vector<std::size_t>& groups = seriesEntry.second;
        for (std::size_t position = 0; position < groups.size(); ++position) {
            const std::size_t targetGroup = groups[position];
            const GroupCorrection& targetAffine =
                corrections[kNonlinearDriftDetectorID][targetGroup];
            if (!targetAffine.valid) continue;

            TH1D aggregate(
                "detector4WindowedShape", "", config::kDriftEnergyBins,
                config::kDriftEnergyMin, config::kDriftEnergyMax);
            aggregate.SetDirectory(nullptr);
            const std::size_t begin = position > 2U ? position-2U : 0U;
            const std::size_t end = std::min(groups.size(), position+3U);
            unsigned int neighbourCount = 0;
            for (std::size_t neighbour = begin; neighbour < end; ++neighbour) {
                const std::size_t sourceGroup = groups[neighbour];
                const GroupCorrection& affine =
                    corrections[kNonlinearDriftDetectorID][sourceGroup];
                if (!affine.valid) continue;
                TH1D* source = dynamic_cast<TH1D*>(inputFile.Get(
                    spectrumPath(static_cast<unsigned int>(sourceGroup),
                                 kNonlinearDriftDetectorID).c_str()));
                if (source == nullptr) continue;
                ++neighbourCount;
                for (int bin = 1; bin <= source->GetNbinsX(); ++bin) {
                    const double content = source->GetBinContent(bin);
                    if (content == 0.0) continue;
                    const double correctedEnergy = affine.offset+
                        affine.slope*source->GetBinCenter(bin);
                    const int destination = aggregate.GetXaxis()->FindBin(
                        correctedEnergy);
                    if (destination < 1 ||
                        destination > aggregate.GetNbinsX()) continue;
                    aggregate.SetBinContent(
                        destination,
                        aggregate.GetBinContent(destination)+content);
                }
            }
            if (neighbourCount == 0U) continue;

            for (const ShapePeakWindowDefinition& definition : definitions) {
                WindowedShapePeakRecord record;
                record.peak.groupIndex = static_cast<unsigned int>(targetGroup);
                record.peak.detectorID = kNonlinearDriftDetectorID;
                record.peak.referenceEnergy = definition.energyKeV;
                record.neighbourCount = neighbourCount;
                const double searchMinimum = definition.energyKeV-
                    definition.searchHalfWidthKeV;
                const double searchMaximum = definition.energyKeV+
                    definition.searchHalfWidthKeV;
                const int firstBin = aggregate.GetXaxis()->FindBin(
                    searchMinimum);
                const int lastBin = aggregate.GetXaxis()->FindBin(
                    searchMaximum);
                int maximumBin = firstBin;
                for (int bin = firstBin+1; bin <= lastBin; ++bin) {
                    if (aggregate.GetBinContent(bin) >
                        aggregate.GetBinContent(maximumBin)) maximumBin = bin;
                }
                if (maximumBin <= firstBin || maximumBin >= lastBin) {
                    records.push_back(record);
                    continue;
                }
                const double maximumEnergy = aggregate.GetBinCenter(maximumBin);
                const double fitMinimum = std::max(
                    searchMinimum,
                    maximumEnergy-definition.fitHalfWidthKeV);
                const double fitMaximum = std::min(
                    searchMaximum,
                    maximumEnergy+definition.fitHalfWidthKeV);
                const int lowBin = aggregate.GetXaxis()->FindBin(fitMinimum);
                const int highBin = aggregate.GetXaxis()->FindBin(fitMaximum);
                const double background = 0.5*(
                    aggregate.GetBinContent(lowBin)+
                    aggregate.GetBinContent(highBin));
                const double height = aggregate.GetBinContent(maximumBin);
                TF1 function(
                    "detector4WindowedPeak",
                    "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]+[4]*(x-[1])",
                    fitMinimum, fitMaximum);
                function.SetParameters(
                    std::max(1.0, height-background), maximumEnergy, 1.2,
                    std::max(0.0, background), 0.0);
                function.SetParLimits(0, 0.0, std::max(1.0, height*10.0));
                function.SetParLimits(1, fitMinimum, fitMaximum);
                function.SetParLimits(2, 0.25, 5.0);
                function.SetParLimits(3, 0.0, std::max(1.0, height*10.0));
                const TFitResultPtr fit = aggregate.Fit(
                    &function, "QRSN", "", fitMinimum, fitMaximum);
                const int fitStatus = static_cast<int>(fit);
                const double centroid = function.GetParameter(1);
                record.centroidErrorKeV = function.GetParError(1);
                record.sigmaKeV = function.GetParameter(2);
                const double amplitude = function.GetParameter(0);
                const double amplitudeError = function.GetParError(0);
                record.significance = amplitudeError > 0.0
                    ? amplitude/amplitudeError : 0.0;
                record.accepted = fitStatus == 0 &&
                    std::isfinite(centroid) &&
                    std::isfinite(record.centroidErrorKeV) &&
                    record.centroidErrorKeV > 0.0 &&
                    record.centroidErrorKeV <= 0.4 &&
                    std::abs(centroid-definition.energyKeV) <=
                        definition.searchHalfWidthKeV &&
                    record.sigmaKeV >= 0.4 && record.sigmaKeV <= 4.0 &&
                    record.significance >= 4.0;
                if (record.accepted) {
                    record.peak.measuredEnergy =
                        (centroid-targetAffine.offset)/targetAffine.slope;
                }
                records.push_back(record);
            }
        }
    }
    return records;
}

std::vector<PiecewiseLinearCorrection> determineDetector4Corrections(
    const std::vector<AcceptedPeakRecord>& records,
    const std::array<std::vector<GroupCorrection>,
                     config::kGermaniumIdBins>& corrections,
    const std::vector<std::size_t>& chronologicalGroups)
{
    const std::size_t groupCount = chronologicalGroups.size();
    const double missing = std::numeric_limits<double>::quiet_NaN();
    std::map<double, std::vector<std::vector<double>>> samples;
    for (const AcceptedPeakRecord& record : records) {
        if (record.detectorID != kNonlinearDriftDetectorID ||
            record.referenceEnergy < 250.0 ||
            record.groupIndex >= groupCount) continue;
        const GroupCorrection& affine =
            corrections[kNonlinearDriftDetectorID][record.groupIndex];
        if (!affine.valid) continue;
        const double residual = affine.offset+
            affine.slope*record.measuredEnergy-record.referenceEnergy;
        if (!std::isfinite(residual) || std::abs(residual) > 12.0) continue;
        auto& byGroup = samples[record.referenceEnergy];
        if (byGroup.empty()) byGroup.resize(groupCount);
        byGroup[record.groupIndex].push_back(residual);
    }

    std::map<double, std::vector<double>> trajectories;
    for (const auto& entry : samples) {
        std::vector<double> raw(groupCount, missing);
        std::size_t measuredGroups = 0;
        for (std::size_t group = 0; group < groupCount; ++group) {
            if (entry.second[group].empty()) continue;
            raw[group] = median(entry.second[group]);
            ++measuredGroups;
        }
        if (measuredGroups >= 20U) {
            trajectories.emplace(
                entry.first, robustTimeTrajectory(raw, chronologicalGroups));
        }
    }

    std::vector<PiecewiseLinearCorrection> result(groupCount);
    for (std::size_t group = 0; group < groupCount; ++group) {
        std::vector<std::pair<double, double>> anchors;
        for (const auto& trajectory : trajectories) {
            const double residual = trajectory.second[group];
            if (!std::isfinite(residual)) continue;
            anchors.emplace_back(trajectory.first+residual, trajectory.first);
        }
        if (anchors.size() < 4U ||
            anchors.back().second-anchors.front().second < 700.0) continue;
        bool valid = true;
        for (std::size_t index = 0; index+1U < anchors.size(); ++index) {
            const double inputSpan = anchors[index+1U].first-
                                     anchors[index].first;
            const double outputSpan = anchors[index+1U].second-
                                      anchors[index].second;
            const double slope = outputSpan/inputSpan;
            if (!(inputSpan > 0.0 && slope > 0.9 && slope < 1.1)) {
                valid = false;
                break;
            }
        }
        if (valid) {
            result[group].valid = true;
            result[group].anchors = std::move(anchors);
        }
    }
    return result;
}

[[maybe_unused]] bool solveFourByFour(
                     std::array<std::array<long double, 5>, 4>& matrix,
                     std::array<long double, 4>& solution)
{
    for (std::size_t column = 0; column < 4; ++column) {
        std::size_t pivot = column;
        for (std::size_t row = column+1; row < 4; ++row) {
            if (std::abs(matrix[row][column]) >
                std::abs(matrix[pivot][column])) {
                pivot = row;
            }
        }
        if (std::abs(matrix[pivot][column]) < 1.0e-18L) return false;
        if (pivot != column) std::swap(matrix[pivot], matrix[column]);
        const long double divisor = matrix[column][column];
        for (std::size_t item = column; item < 5; ++item) {
            matrix[column][item] /= divisor;
        }
        for (std::size_t row = 0; row < 4; ++row) {
            if (row == column) continue;
            const long double factor = matrix[row][column];
            for (std::size_t item = column; item < 5; ++item) {
                matrix[row][item] -= factor*matrix[column][item];
            }
        }
    }
    for (std::size_t row = 0; row < 4; ++row) {
        solution[row] = matrix[row][4];
    }
    return true;
}

[[maybe_unused]] ShapeCorrection makeShapeCorrection(double firstResidual,
                                    double secondResidual)
{
    ShapeCorrection result;
    if (!std::isfinite(firstResidual) || !std::isfinite(secondResidual)) {
        return result;
    }
    result.firstAnchorResidual = firstResidual;
    result.secondAnchorResidual = secondResidual;
    const long double span = kShapeUpperEnergyKeV-kShapeLowerEnergyKeV;
    const long double firstInput =
        kFirstShapeAnchorKeV+result.firstAnchorResidual;
    const long double secondInput =
        kSecondShapeAnchorKeV+result.secondAnchorResidual;
    const long double firstT =
        (firstInput-kShapeLowerEnergyKeV)/span;
    const long double secondT =
        (secondInput-kShapeLowerEnergyKeV)/span;
    if (!(firstT > 0.0L && firstT < secondT && secondT < 1.0L)) {
        return result;
    }

    // delta(t)=c2*t^2+...+c5*t^5.  Value and slope are zero at
    // t=0 and t=1, while the two in-beam anchors set the interior shape.
    std::array<std::array<long double, 5>, 4> equations{};
    const auto fillValueRow = [&equations](std::size_t row, long double t,
                                           long double value) {
        long double power = t*t;
        for (std::size_t column = 0; column < 4; ++column) {
            equations[row][column] = power;
            power *= t;
        }
        equations[row][4] = value;
    };
    fillValueRow(0, firstT, -result.firstAnchorResidual);
    fillValueRow(1, secondT, -result.secondAnchorResidual);
    equations[2] = {1.0L, 1.0L, 1.0L, 1.0L, 0.0L};
    equations[3] = {2.0L, 3.0L, 4.0L, 5.0L, 0.0L};
    std::array<long double, 4> solved{};
    if (!solveFourByFour(equations, solved)) return result;
    for (std::size_t index = 0; index < solved.size(); ++index) {
        result.normalizedCoefficients[index+2U] =
            static_cast<double>(solved[index]);
    }

    // Guard against an accidental non-monotonic energy calibration.
    for (unsigned int step = 0; step <= 1000; ++step) {
        const double t = static_cast<double>(step)/1000.0;
        double derivative = 1.0;
        double power = t;
        for (std::size_t degree = 2; degree < 6; ++degree) {
            derivative += degree*result.normalizedCoefficients[degree]*power/
                          static_cast<double>(span);
            power *= t;
        }
        if (!std::isfinite(derivative) || derivative <= 0.9 ||
            derivative >= 1.1) {
            return ShapeCorrection{};
        }
    }
    result.valid = true;
    return result;
}

[[maybe_unused]] std::vector<ShapeCorrection> determineShapeCorrections(
    const std::vector<AcceptedPeakRecord>& records,
    const std::array<std::vector<GroupCorrection>,
                     config::kGermaniumIdBins>& corrections,
    const std::vector<std::size_t>& chronologicalGroups)
{
    const std::size_t groupCount = chronologicalGroups.size();
    const double missing = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> first(groupCount, missing);
    std::vector<double> second(groupCount, missing);
    std::vector<std::vector<double>> firstValues(groupCount);
    std::vector<std::vector<double>> secondValues(groupCount);
    for (const AcceptedPeakRecord& record : records) {
        if (record.detectorID != kNonlinearDriftDetectorID ||
            record.groupIndex >= groupCount) {
            continue;
        }
        const GroupCorrection& correction =
            corrections[record.detectorID][record.groupIndex];
        if (!correction.valid) continue;
        const double residual = correction.offset+
            correction.slope*record.measuredEnergy-record.referenceEnergy;
        if (!std::isfinite(residual) || std::abs(residual) > 8.0) continue;
        if (std::abs(record.referenceEnergy-kFirstShapeAnchorKeV) < 0.1) {
            firstValues[record.groupIndex].push_back(residual);
        } else if (std::abs(record.referenceEnergy-
                            kSecondShapeAnchorKeV) < 0.1) {
            secondValues[record.groupIndex].push_back(residual);
        }
    }
    for (std::size_t group = 0; group < groupCount; ++group) {
        if (!firstValues[group].empty()) first[group] = median(firstValues[group]);
        if (!secondValues[group].empty()) second[group] = median(secondValues[group]);
    }

    const auto smooth = [&chronologicalGroups](const std::vector<double>& raw) {
        std::vector<double> result(raw.size(),
            std::numeric_limits<double>::quiet_NaN());
        std::vector<double> all;
        for (double value : raw) if (std::isfinite(value)) all.push_back(value);
        const double fallback = median(all);
        for (std::size_t position = 0; position < chronologicalGroups.size();
             ++position) {
            std::vector<double> window;
            const std::size_t begin = position > 3U ? position-3U : 0U;
            const std::size_t end = std::min(
                chronologicalGroups.size(), position+4U);
            for (std::size_t neighbour = begin; neighbour < end; ++neighbour) {
                const double value = raw[chronologicalGroups[neighbour]];
                if (std::isfinite(value)) window.push_back(value);
            }
            if (window.empty()) {
                result[chronologicalGroups[position]] = fallback;
                continue;
            }
            const double centre = median(window);
            std::vector<double> deviations;
            deviations.reserve(window.size());
            for (double value : window) {
                deviations.push_back(std::abs(value-centre));
            }
            const double limit = std::max(0.6, 4.0*1.4826*median(deviations));
            std::vector<double> retained;
            for (double value : window) {
                if (std::abs(value-centre) <= limit) retained.push_back(value);
            }
            const double localCentre = retained.empty() ? centre
                                                        : median(retained);
            const double current = raw[chronologicalGroups[position]];
            if (std::isfinite(current) &&
                std::abs(current-centre) <= limit) {
                // Preserve genuine group-to-group curvature while damping
                // centroid noise. Outliers and missing anchors use the robust
                // local trajectory directly.
                result[chronologicalGroups[position]] =
                    0.65*current+0.35*localCentre;
            } else {
                result[chronologicalGroups[position]] = localCentre;
            }
        }
        return result;
    };

    const std::vector<double> smoothedFirst = smooth(first);
    const std::vector<double> smoothedSecond = smooth(second);
    std::vector<ShapeCorrection> result(groupCount);
    for (std::size_t group = 0; group < groupCount; ++group) {
        result[group] = makeShapeCorrection(
            smoothedFirst[group], smoothedSecond[group]);
    }
    return result;
}

[[maybe_unused]] std::array<double, 6> composeShapePolynomial(
    const ShapeCorrection& shape, const GroupCorrection& correction)
{
    const long double span = kShapeUpperEnergyKeV-kShapeLowerEnergyKeV;
    const long double constant =
        (correction.offset-kShapeLowerEnergyKeV)/span;
    const long double linear = correction.slope/span;
    std::array<long double, 6> result{};
    result[0] = correction.offset;
    result[1] = correction.slope;
    std::array<long double, 6> power{};
    power[0] = 1.0L;
    for (std::size_t degree = 0; degree < 6; ++degree) {
        if (degree >= 2U) {
            const long double coefficient =
                shape.normalizedCoefficients[degree];
            for (std::size_t item = 0; item < 6; ++item) {
                result[item] += coefficient*power[item];
            }
        }
        std::array<long double, 6> next{};
        for (std::size_t item = 0; item < 6; ++item) {
            next[item] += constant*power[item];
            if (item+1U < 6U) next[item+1U] += linear*power[item];
        }
        power = next;
    }
    std::array<double, 6> converted{};
    for (std::size_t index = 0; index < 6; ++index) {
        converted[index] = static_cast<double>(result[index]);
    }
    return converted;
}

[[maybe_unused]] double applyShapeCorrection(
    const ShapeCorrection& shape, double energy)
{
    if (!shape.valid || energy <= kShapeLowerEnergyKeV ||
        energy >= kShapeUpperEnergyKeV) {
        return energy;
    }
    const double t = (energy-kShapeLowerEnergyKeV)/
                     (kShapeUpperEnergyKeV-kShapeLowerEnergyKeV);
    double delta = 0.0;
    for (auto coefficient = shape.normalizedCoefficients.rbegin();
         coefficient != shape.normalizedCoefficients.rend(); ++coefficient) {
        delta = delta*t+*coefficient;
    }
    return energy+delta;
}

struct SpectrumAlignment {
    double scale = 1.0;
    double offset = 0.0;
    double score = 0.0;
};

double peakExcessScore(TH1D& histogram, double energy)
{
    const int centre = histogram.GetXaxis()->FindBin(energy);
    // All configured anchors are well inside the 0--2048 keV spectrum, even
    // at the limits of the alignment search, so the sidebands remain valid.

    double peak = 0.0;
    for (int offset = -4; offset <= 4; ++offset) {
        peak = std::max(peak, histogram.GetBinContent(centre+offset));
    }
    double background = 0.0;
    int backgroundBins = 0;
    for (int offset = 16; offset <= 32; ++offset) {
        background += histogram.GetBinContent(centre-offset);
        background += histogram.GetBinContent(centre+offset);
        backgroundBins += 2;
    }
    background /= static_cast<double>(backgroundBins);
    const double significance = std::max(
        0.0, (peak-background)/std::sqrt(std::max(1.0, background)));
    return std::log1p(std::min(1000.0, significance));
}

double alignmentScore(TH1D& histogram,
                      const std::vector<double>& referenceEnergies,
                      double scale,
                      double offset)
{
    double score = 0.0;
    for (const double reference : referenceEnergies) {
        score += peakExcessScore(histogram, offset+scale*reference);
    }
    return score;
}

SpectrumAlignment findSpectrumAlignment(
    TH1D& histogram, const std::vector<double>& referenceEnergies)
{
    SpectrumAlignment best;
    if (referenceEnergies.size() < 2) return best;
    best.score = -1.0;

    // actual spectrum energy = offset + scale*reference energy.  The broad
    // first pass is deliberately capable of finding the approximately 4%
    // detector-2 jump present before 9um run 213.
    for (int scaleStep = 0; scaleStep <= 60; ++scaleStep) {
        const double scale = 0.94+0.002*scaleStep;
        for (int offsetStep = -15; offsetStep <= 15; ++offsetStep) {
            const double offset = static_cast<double>(offsetStep);
            const double score = alignmentScore(
                histogram, referenceEnergies, scale, offset);
            if (score > best.score) best = {scale, offset, score};
        }
    }

    const SpectrumAlignment coarse = best;
    for (int scaleStep = -10; scaleStep <= 10; ++scaleStep) {
        const double scale = coarse.scale+0.0002*scaleStep;
        for (int offsetStep = -10; offsetStep <= 10; ++offsetStep) {
            const double offset = coarse.offset+0.1*offsetStep;
            const double score = alignmentScore(
                histogram, referenceEnergies, scale, offset);
            if (score > best.score) best = {scale, offset, score};
        }
    }
    return best;
}

struct LinearFit {
    bool valid = false;
    double offset = 0.0;
    double slope = 1.0;
    double offsetError = 0.0;
    double slopeError = 0.0;
    double rms = 0.0;
    unsigned int usedPoints = 0;
};

struct QuadraticFit {
    bool valid = false;
    std::array<double, 3> coefficients{0.0, 1.0, 0.0};
    double rms = 0.0;
    unsigned int usedPoints = 0;
};

bool solveThreeByThree(std::array<std::array<long double, 4>, 3>& matrix,
                       std::array<long double, 3>& solution)
{
    for (std::size_t column = 0; column < 3; ++column) {
        std::size_t pivot = column;
        for (std::size_t row = column+1; row < 3; ++row) {
            if (std::abs(matrix[row][column]) >
                std::abs(matrix[pivot][column])) pivot = row;
        }
        if (std::abs(matrix[pivot][column]) < 1.0e-18L) return false;
        if (pivot != column) std::swap(matrix[pivot], matrix[column]);
        const long double divisor = matrix[column][column];
        for (std::size_t item = column; item < 4; ++item) {
            matrix[column][item] /= divisor;
        }
        for (std::size_t row = 0; row < 3; ++row) {
            if (row == column) continue;
            const long double factor = matrix[row][column];
            for (std::size_t item = column; item < 4; ++item) {
                matrix[row][item] -= factor*matrix[column][item];
            }
        }
    }
    for (std::size_t row = 0; row < 3; ++row) solution[row] = matrix[row][3];
    return true;
}

double evaluatePolynomial(const std::array<double, 3>& coefficients,
                          double energy)
{
    return coefficients[0]+energy*(coefficients[1]+energy*coefficients[2]);
}

QuadraticFit weightedQuadraticFit(std::vector<PeakMeasurement> points,
                                  unsigned int minimumPoints,
                                  double minimumSpan)
{
    QuadraticFit result;
    for (unsigned int iteration = 0; iteration < 4; ++iteration) {
        std::array<long double, 5> sums{};
        std::array<long double, 3> right{};
        double minimumEnergy = std::numeric_limits<double>::infinity();
        double maximumEnergy = -std::numeric_limits<double>::infinity();
        unsigned int count = 0;
        for (const PeakMeasurement& point : points) {
            if (!point.used) continue;
            const long double x = point.measuredEnergy;
            const long double y = point.referenceEnergy;
            const long double uncertainty = std::max(0.03, point.measuredError);
            const long double weight = 1.0L/(uncertainty*uncertainty);
            long double power = 1.0L;
            for (std::size_t degree = 0; degree < sums.size(); ++degree) {
                sums[degree] += weight*power;
                power *= x;
            }
            power = 1.0L;
            for (std::size_t degree = 0; degree < right.size(); ++degree) {
                right[degree] += weight*y*power;
                power *= x;
            }
            minimumEnergy = std::min(minimumEnergy, point.referenceEnergy);
            maximumEnergy = std::max(maximumEnergy, point.referenceEnergy);
            ++count;
        }
        if (count < minimumPoints || maximumEnergy-minimumEnergy < minimumSpan) {
            return result;
        }
        std::array<std::array<long double, 4>, 3> equations{};
        for (std::size_t row = 0; row < 3; ++row) {
            for (std::size_t column = 0; column < 3; ++column) {
                equations[row][column] = sums[row+column];
            }
            equations[row][3] = right[row];
        }
        std::array<long double, 3> solved{};
        if (!solveThreeByThree(equations, solved)) return result;
        for (std::size_t index = 0; index < 3; ++index) {
            result.coefficients[index] = static_cast<double>(solved[index]);
        }
        result.usedPoints = count;

        double sumResidualSquared = 0.0;
        double largestResidual = 0.0;
        std::size_t largestIndex = points.size();
        for (std::size_t index = 0; index < points.size(); ++index) {
            if (!points[index].used) continue;
            const double residual = points[index].referenceEnergy-
                evaluatePolynomial(result.coefficients,
                                   points[index].measuredEnergy);
            sumResidualSquared += residual*residual;
            if (std::abs(residual) > largestResidual) {
                largestResidual = std::abs(residual);
                largestIndex = index;
            }
        }
        result.rms = std::sqrt(sumResidualSquared/static_cast<double>(count));
        if (largestResidual <= 1.5 || count <= minimumPoints ||
            largestIndex == points.size()) {
            const double derivativeLow = result.coefficients[1];
            const double derivativeHigh = result.coefficients[1]+
                2.0*result.coefficients[2]*config::kDriftEnergyMax;
            result.valid = std::isfinite(result.coefficients[0]) &&
                std::isfinite(result.coefficients[1]) &&
                std::isfinite(result.coefficients[2]) &&
                derivativeLow > 0.94 && derivativeLow < 1.06 &&
                derivativeHigh > 0.94 && derivativeHigh < 1.06 &&
                std::abs(result.coefficients[0]) < 30.0 && result.rms <= 1.5;
            return result;
        }
        points[largestIndex].used = false;
    }
    return result;
}

LinearFit weightedLinearFit(std::vector<PeakMeasurement>& points,
                            unsigned int minimumPoints,
                            double minimumSpan)
{
    LinearFit result;
    for (unsigned int iteration = 0; iteration < 3; ++iteration) {
        double sumW = 0.0;
        double sumX = 0.0;
        double sumY = 0.0;
        double sumXX = 0.0;
        double sumXY = 0.0;
        double minimumEnergy = std::numeric_limits<double>::infinity();
        double maximumEnergy = -std::numeric_limits<double>::infinity();
        unsigned int count = 0;
        for (const PeakMeasurement& point : points) {
            if (!point.used) continue;
            const double uncertainty = std::max(0.03, point.measuredError);
            const double weight = 1.0/(uncertainty*uncertainty);
            sumW += weight;
            sumX += weight*point.measuredEnergy;
            sumY += weight*point.referenceEnergy;
            sumXX += weight*point.measuredEnergy*point.measuredEnergy;
            sumXY += weight*point.measuredEnergy*point.referenceEnergy;
            minimumEnergy = std::min(minimumEnergy, point.referenceEnergy);
            maximumEnergy = std::max(maximumEnergy, point.referenceEnergy);
            ++count;
        }
        if (count < minimumPoints || maximumEnergy-minimumEnergy < minimumSpan) {
            return result;
        }
        const double determinant = sumW*sumXX-sumX*sumX;
        if (!std::isfinite(determinant) || determinant <= 0.0) return result;
        result.offset = (sumY*sumXX-sumX*sumXY)/determinant;
        result.slope = (sumW*sumXY-sumX*sumY)/determinant;
        result.offsetError = std::sqrt(sumXX/determinant);
        result.slopeError = std::sqrt(sumW/determinant);
        result.usedPoints = count;

        double sumResidualSquared = 0.0;
        double largestResidual = 0.0;
        std::size_t largestIndex = points.size();
        for (std::size_t index = 0; index < points.size(); ++index) {
            const PeakMeasurement& point = points[index];
            if (!point.used) continue;
            const double residual = point.referenceEnergy-
                (result.offset+result.slope*point.measuredEnergy);
            sumResidualSquared += residual*residual;
            if (std::abs(residual) > largestResidual) {
                largestResidual = std::abs(residual);
                largestIndex = index;
            }
        }
        result.rms = std::sqrt(sumResidualSquared/static_cast<double>(count));
        if (largestResidual <= 1.5 || count <= minimumPoints ||
            largestIndex == points.size()) {
            result.valid = std::isfinite(result.offset) &&
                std::isfinite(result.slope) && result.slope > 0.95 &&
                result.slope < 1.05 && std::abs(result.offset) < 25.0 &&
                result.rms <= 1.5;
            return result;
        }
        points[largestIndex].used = false;
    }
    return result;
}

std::vector<std::size_t> chronologicalGroupOrder(
    const std::vector<FileRecord>& files,
    std::size_t groupCount,
    const std::vector<std::string>& requestedSeriesOrder)
{
    std::vector<std::string> encountered;
    std::map<std::string, std::size_t> firstSeen;
    std::vector<std::string> groupSeries(groupCount);
    std::vector<unsigned int> groupRun(
        groupCount, std::numeric_limits<unsigned int>::max());
    for (const FileRecord& file : files) {
        if (firstSeen.emplace(file.series, firstSeen.size()).second) {
            encountered.push_back(file.series);
        }
        groupSeries[file.groupIndex] = file.series;
        groupRun[file.groupIndex] = std::min(groupRun[file.groupIndex], file.run);
    }
    std::map<std::string, std::size_t> rank;
    for (const std::string& series : requestedSeriesOrder) {
        if (!rank.emplace(series, rank.size()).second) {
            throw std::runtime_error("duplicate series '"+series+
                                     "' in --series-order");
        }
    }
    for (const std::string& series : encountered) {
        if (rank.count(series) == 0U) rank.emplace(series, rank.size());
    }
    for (const std::string& series : requestedSeriesOrder) {
        if (firstSeen.count(series) == 0U) {
            throw std::runtime_error("series '"+series+
                                     "' in --series-order is not in input");
        }
    }
    std::vector<std::size_t> order(groupCount);
    for (std::size_t group = 0; group < groupCount; ++group) order[group] = group;
    std::stable_sort(order.begin(), order.end(),
        [&rank, &groupSeries, &groupRun](std::size_t left, std::size_t right) {
            const std::size_t leftRank = rank.at(groupSeries[left]);
            const std::size_t rightRank = rank.at(groupSeries[right]);
            if (leftRank != rightRank) return leftRank < rightRank;
            if (groupRun[left] != groupRun[right]) {
                return groupRun[left] < groupRun[right];
            }
            return left < right;
        });
    return order;
}

void fillMissingCorrections(std::vector<GroupCorrection>& corrections,
                            const std::vector<std::size_t>& order)
{
    std::vector<std::size_t> validPositions;
    for (std::size_t position = 0; position < order.size(); ++position) {
        if (corrections[order[position]].valid) validPositions.push_back(position);
    }
    if (validPositions.empty()) return;

    for (std::size_t position = 0; position < order.size(); ++position) {
        const std::size_t index = order[position];
        if (corrections[index].valid) continue;
        const auto after = std::lower_bound(
            validPositions.begin(), validPositions.end(), position);
        if (after == validPositions.begin()) {
            corrections[index] = corrections[order[*after]];
            corrections[index].source = CorrectionSource::Nearest;
        } else if (after == validPositions.end()) {
            corrections[index] = corrections[order[validPositions.back()]];
            corrections[index].source = CorrectionSource::Nearest;
        } else {
            const std::size_t highPosition = *after;
            const std::size_t lowPosition = *(after-1);
            const std::size_t high = order[highPosition];
            const std::size_t low = order[lowPosition];
            const double fraction = static_cast<double>(position-lowPosition)/
                                    static_cast<double>(highPosition-lowPosition);
            GroupCorrection interpolated;
            interpolated.valid = true;
            interpolated.source = CorrectionSource::Interpolated;
            interpolated.offset = corrections[low].offset + fraction*
                (corrections[high].offset-corrections[low].offset);
            interpolated.slope = corrections[low].slope + fraction*
                (corrections[high].slope-corrections[low].slope);
            interpolated.offsetError = std::max(
                corrections[low].offsetError,
                corrections[high].offsetError);
            interpolated.slopeError = std::max(
                corrections[low].slopeError,
                corrections[high].slopeError);
            corrections[index] = interpolated;
        }
    }
}

void fillMissingPolynomialCorrections(
    std::vector<PolynomialCorrection>& corrections,
    const std::vector<std::size_t>& order)
{
    std::vector<std::size_t> validPositions;
    for (std::size_t position = 0; position < order.size(); ++position) {
        if (corrections[order[position]].valid) validPositions.push_back(position);
    }
    if (validPositions.empty()) return;
    for (std::size_t position = 0; position < order.size(); ++position) {
        const std::size_t index = order[position];
        if (corrections[index].valid) continue;
        const auto after = std::lower_bound(
            validPositions.begin(), validPositions.end(), position);
        if (after == validPositions.begin() || after == validPositions.end()) {
            const std::size_t sourcePosition = after == validPositions.begin()
                ? *after : validPositions.back();
            corrections[index] = corrections[order[sourcePosition]];
            corrections[index].source = CorrectionSource::Nearest;
            continue;
        }
        const std::size_t highPosition = *after;
        const std::size_t lowPosition = *(after-1);
        const double fraction = static_cast<double>(position-lowPosition)/
                                static_cast<double>(highPosition-lowPosition);
        const PolynomialCorrection& low = corrections[order[lowPosition]];
        const PolynomialCorrection& high = corrections[order[highPosition]];
        PolynomialCorrection interpolated;
        interpolated.valid = true;
        interpolated.source = CorrectionSource::Interpolated;
        for (std::size_t coefficient = 0; coefficient < 3; ++coefficient) {
            interpolated.coefficients[coefficient] = low.coefficients[coefficient]+
                fraction*(high.coefficients[coefficient]-
                          low.coefficients[coefficient]);
        }
        interpolated.nonlinear = low.nonlinear || high.nonlinear;
        corrections[index] = interpolated;
    }
}

std::string correctionFileName(const std::string& prefix,
                               const std::string& series)
{
    return prefix + "_" + series + ".cal";
}

std::string nonlinearCorrectionFileName(const std::string& prefix,
                                        const std::string& series)
{
    return prefix + "_" + series + "_det4.mcal";
}

} // namespace

void DriftCalibrationFitter::loadPeakDefinitions(const std::string& fileName)
{
    std::ifstream file(fileName);
    if (!file) {
        throw std::runtime_error("Could not open drift-peak file '" +
                                 fileName + "'");
    }
    peaks_.clear();
    std::string line;
    std::size_t lineNumber = 0;
    while (std::getline(file, line)) {
        ++lineNumber;
        if (!prepareDataLine(line)) continue;
        PeakDefinition peak;
        std::istringstream input(line);
        if (!(input >> peak.energyKeV >> peak.label >>
              peak.searchHalfWidthKeV >> peak.fitHalfWidthKeV)) {
            throw parseError(fileName, lineNumber,
                "expected ENERGY_KEV LABEL SEARCH_HALF_WIDTH_KEV "
                "FIT_HALF_WIDTH_KEV [EXCLUDED_IDS]");
        }
        std::string exclusions;
        if (input >> exclusions) {
            peak.excludedDetectorIDs =
                parseDetectorList(exclusions, fileName, lineNumber);
        }
        std::string extra;
        if (input >> extra) {
            throw parseError(fileName, lineNumber,
                             "unexpected extra value '" + extra + "'");
        }
        if (!std::isfinite(peak.energyKeV) || peak.energyKeV <= 0.0 ||
            !std::isfinite(peak.searchHalfWidthKeV) ||
            peak.searchHalfWidthKeV <= 0.0 ||
            !std::isfinite(peak.fitHalfWidthKeV) ||
            peak.fitHalfWidthKeV <= 0.0 ||
            peak.fitHalfWidthKeV > peak.searchHalfWidthKeV) {
            throw parseError(fileName, lineNumber,
                             "invalid energy or fit/search half-width");
        }
        peaks_.push_back(std::move(peak));
    }
    if (peaks_.empty()) {
        throw std::runtime_error("Drift-peak file contains no peaks");
    }
}

void DriftCalibrationFitter::setMinimumPeakCount(unsigned int count)
{
    if (count < 2) {
        throw std::runtime_error("minimum peak count must be at least two");
    }
    minimumPeakCount_ = count;
}

void DriftCalibrationFitter::setMinimumEnergySpan(double energySpanKeV)
{
    if (!std::isfinite(energySpanKeV) || energySpanKeV <= 0.0) {
        throw std::runtime_error("minimum energy span must be positive");
    }
    minimumEnergySpanKeV_ = energySpanKeV;
}

void DriftCalibrationFitter::setSeriesOrder(
    const std::vector<std::string>& seriesOrder)
{
    if (seriesOrder.empty()) {
        throw std::runtime_error("series order must not be empty");
    }
    seriesOrder_ = seriesOrder;
}

int DriftCalibrationFitter::run(const std::string& inputFileName,
                                const std::string& outputPrefix,
                                const std::string& diagnosticFileName) const
{
    if (peaks_.empty()) {
        std::cerr << "Error: no drift peaks were loaded.\n";
        return 2;
    }
    TFile inputFile(inputFileName.c_str(), "READ");
    if (inputFile.IsZombie()) {
        std::cerr << "Error: could not open '" << inputFileName << "'.\n";
        return 3;
    }
    TTree* provenance = dynamic_cast<TTree*>(inputFile.Get("DriftInputFiles"));
    if (provenance == nullptr) {
        std::cerr << "Error: input contains no DriftInputFiles tree.\n";
        return 4;
    }

    std::vector<FileRecord> files;
    unsigned int maximumGroup = 0;
    {
        TTreeReader reader(provenance);
        TTreeReaderValue<unsigned int> groupIndex(reader, "groupIndex");
        TTreeReaderValue<std::string> series(reader, "series");
        TTreeReaderValue<unsigned int> run(reader, "run");
        TTreeReaderValue<std::string> file(reader, "file");
        while (reader.Next()) {
            files.push_back({*groupIndex, *series, *run, *file});
            maximumGroup = std::max(maximumGroup,
                                    static_cast<unsigned int>(*groupIndex));
        }
    }
    if (files.empty()) {
        std::cerr << "Error: drift input contains no file provenance.\n";
        return 5;
    }
    const std::size_t groupCount = static_cast<std::size_t>(maximumGroup)+1U;
    const std::vector<std::size_t> chronologicalGroups =
        chronologicalGroupOrder(files, groupCount, seriesOrder_);
    if (!seriesOrder_.empty()) {
        std::cout << "Using requested chronological series order: ";
        for (std::size_t index = 0; index < seriesOrder_.size(); ++index) {
            if (index != 0U) std::cout << ", ";
            std::cout << seriesOrder_[index];
        }
        std::cout << ".\n";
    }

    TFile diagnostics(diagnosticFileName.c_str(), "RECREATE");
    if (diagnostics.IsZombie()) {
        std::cerr << "Error: could not create '" << diagnosticFileName
                  << "'.\n";
        return 6;
    }

    unsigned int fitGroupIndex = 0;
    UShort_t fitDetectorID = 0;
    std::string fitLabel;
    double fitReferenceEnergy = 0.0;
    double fitCentroid = 0.0;
    double fitCentroidError = 0.0;
    double fitSigma = 0.0;
    double fitAmplitude = 0.0;
    double fitAmplitudeError = 0.0;
    double fitChiSquare = 0.0;
    int fitNdf = 0;
    int fitStatus = -1;
    bool fitAccepted = false;
    TTree peakTree("DriftPeakFits", "Drift-reference peak fits");
    peakTree.Branch("groupIndex", &fitGroupIndex);
    peakTree.Branch("detectorID", &fitDetectorID);
    peakTree.Branch("label", &fitLabel);
    peakTree.Branch("referenceEnergyKeV", &fitReferenceEnergy);
    peakTree.Branch("centroidKeV", &fitCentroid);
    peakTree.Branch("centroidErrorKeV", &fitCentroidError);
    peakTree.Branch("sigmaKeV", &fitSigma);
    peakTree.Branch("amplitude", &fitAmplitude);
    peakTree.Branch("amplitudeError", &fitAmplitudeError);
    peakTree.Branch("chiSquare", &fitChiSquare);
    peakTree.Branch("ndf", &fitNdf);
    peakTree.Branch("fitStatus", &fitStatus);
    peakTree.Branch("accepted", &fitAccepted);

    unsigned int alignmentGroupIndex = 0;
    UShort_t alignmentDetectorID = 0;
    double alignmentScale = 1.0;
    double alignmentOffset = 0.0;
    double alignmentFitScore = 0.0;
    TTree alignmentTree(
        "DriftSeedAlignment",
        "Detector/group affine alignment used to locate drift peaks");
    alignmentTree.Branch("groupIndex", &alignmentGroupIndex);
    alignmentTree.Branch("detectorID", &alignmentDetectorID);
    alignmentTree.Branch("scale", &alignmentScale);
    alignmentTree.Branch("offsetKeV", &alignmentOffset);
    alignmentTree.Branch("score", &alignmentFitScore);

    std::array<std::vector<GroupCorrection>, config::kGermaniumIdBins>
        corrections;
    for (auto& detectorCorrections : corrections) {
        detectorCorrections.resize(groupCount);
    }
    std::vector<PolynomialCorrection> detector16Corrections(groupCount);
    std::vector<AcceptedPeakRecord> acceptedPeakRecords;

    for (unsigned int group = 0; group < groupCount; ++group) {
        for (unsigned int detector = 0;
             detector < config::kGermaniumIdBins; ++detector) {
            TH1D* histogram = dynamic_cast<TH1D*>(
                inputFile.Get(spectrumPath(group, detector).c_str()));
            std::vector<PeakMeasurement> measurements;

            std::vector<double> alignmentEnergies;
            if (histogram != nullptr && histogram->GetEntries() > 0.0) {
                for (const PeakDefinition& peak : peaks_) {
                    if (peak.excludedDetectorIDs.count(detector) == 0U) {
                        alignmentEnergies.push_back(peak.energyKeV);
                    }
                }
            }
            const SpectrumAlignment alignment =
                histogram == nullptr
                    ? SpectrumAlignment{}
                    : findSpectrumAlignment(*histogram, alignmentEnergies);
            alignmentGroupIndex = group;
            alignmentDetectorID = static_cast<UShort_t>(detector);
            alignmentScale = alignment.scale;
            alignmentOffset = alignment.offset;
            alignmentFitScore = alignment.score;
            alignmentTree.Fill();

            for (const PeakDefinition& peak : peaks_) {
                fitGroupIndex = group;
                fitDetectorID = static_cast<UShort_t>(detector);
                fitLabel = peak.label;
                fitReferenceEnergy = peak.energyKeV;
                fitCentroid = fitCentroidError = fitSigma = 0.0;
                fitAmplitude = fitAmplitudeError = fitChiSquare = 0.0;
                fitNdf = 0;
                fitStatus = -1;
                fitAccepted = false;

                if (peak.excludedDetectorIDs.count(detector) != 0U) {
                    fitStatus = -5;
                    peakTree.Fill();
                    continue;
                }
                if (histogram == nullptr || histogram->GetEntries() == 0.0) {
                    fitStatus = -2;
                    peakTree.Fill();
                    continue;
                }

                const double predictedEnergy =
                    alignment.offset+alignment.scale*peak.energyKeV;
                const double searchMinimum =
                    predictedEnergy-peak.searchHalfWidthKeV;
                const double searchMaximum =
                    predictedEnergy+peak.searchHalfWidthKeV;
                const int firstBin = histogram->GetXaxis()->FindBin(searchMinimum);
                const int lastBin = histogram->GetXaxis()->FindBin(searchMaximum);
                int maximumBin = firstBin;
                for (int bin = firstBin+1; bin <= lastBin; ++bin) {
                    if (histogram->GetBinContent(bin) >
                        histogram->GetBinContent(maximumBin)) {
                        maximumBin = bin;
                    }
                }
                const double maximumEnergy = histogram->GetBinCenter(maximumBin);
                if (maximumBin <= firstBin || maximumBin >= lastBin) {
                    fitStatus = -4;
                    peakTree.Fill();
                    continue;
                }
                const double fitMinimum = std::max(
                    searchMinimum, maximumEnergy-peak.fitHalfWidthKeV);
                const double fitMaximum = std::min(
                    searchMaximum, maximumEnergy+peak.fitHalfWidthKeV);
                const int lowBin = histogram->GetXaxis()->FindBin(fitMinimum);
                const int highBin = histogram->GetXaxis()->FindBin(fitMaximum);
                const double edgeBackground = 0.5*(
                    histogram->GetBinContent(lowBin)+
                    histogram->GetBinContent(highBin));
                const double peakHeight = histogram->GetBinContent(maximumBin);

                TF1 function(
                    "driftPeakFit",
                    "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]+[4]*(x-[1])",
                    fitMinimum, fitMaximum);
                function.SetParameters(
                    std::max(1.0, peakHeight-edgeBackground), maximumEnergy,
                    0.8, std::max(0.0, edgeBackground), 0.0);
                function.SetParNames(
                    "Amplitude", "Centroid", "Sigma", "Background", "Slope");
                function.SetParLimits(0, 0.0, std::max(1.0, peakHeight*10.0));
                function.SetParLimits(1, fitMinimum, fitMaximum);
                function.SetParLimits(2, 0.15, 5.0);
                function.SetParLimits(3, 0.0, std::max(1.0, peakHeight*10.0));
                const TFitResultPtr result = histogram->Fit(
                    &function, "QRSN", "", fitMinimum, fitMaximum);
                fitStatus = static_cast<int>(result);
                fitCentroid = function.GetParameter(1);
                fitCentroidError = function.GetParError(1);
                fitSigma = function.GetParameter(2);
                fitAmplitude = function.GetParameter(0);
                fitAmplitudeError = function.GetParError(0);
                fitChiSquare = function.GetChisquare();
                fitNdf = function.GetNDF();
                const double amplitudeSignificance = fitAmplitudeError > 0.0
                    ? fitAmplitude/fitAmplitudeError : 0.0;
                fitAccepted = fitStatus == 0 &&
                    std::isfinite(fitCentroid) &&
                    std::isfinite(fitCentroidError) &&
                    fitCentroidError > 0.0 && fitCentroidError <= 1.0 &&
                    std::abs(fitCentroid-predictedEnergy) <=
                        peak.searchHalfWidthKeV &&
                    fitSigma >= 0.2 && fitSigma <= 4.0 &&
                    amplitudeSignificance >= 5.0;
                if (fitAccepted) {
                    acceptedPeakRecords.push_back({
                        group, detector, peak.energyKeV, fitCentroid});
                    if (useForDriftFit(detector, peak.energyKeV)) {
                        measurements.push_back({fitCentroid, fitCentroidError,
                                                peak.energyKeV, true});
                    }
                }
                peakTree.Fill();
            }

            const unsigned int minimumPoints =
                detector == kNonlinearDriftDetectorID
                    ? std::max(3U, minimumPeakCount_) : minimumPeakCount_;
            const double minimumSpan =
                detector == kNonlinearDriftDetectorID
                    ? std::max(700.0, minimumEnergySpanKeV_)
                    : minimumEnergySpanKeV_;
            const std::vector<PeakMeasurement> originalMeasurements = measurements;
            LinearFit fit = weightedLinearFit(measurements, minimumPoints,
                                              minimumSpan);
            if (fit.valid) {
                GroupCorrection& correction = corrections[detector][group];
                correction.valid = true;
                correction.source = CorrectionSource::Measured;
                correction.offset = fit.offset;
                correction.slope = fit.slope;
                correction.offsetError = fit.offsetError;
                correction.slopeError = fit.slopeError;
                correction.rms = fit.rms;
                correction.peakCount = fit.usedPoints;

                if (detector == kQuadraticDriftDetectorID) {
                    PolynomialCorrection& polynomial =
                        detector16Corrections[group];
                    polynomial.valid = true;
                    polynomial.source = CorrectionSource::Measured;
                    polynomial.coefficients = {fit.offset, fit.slope, 0.0};
                    polynomial.rms = fit.rms;
                    polynomial.peakCount = fit.usedPoints;
                    const QuadraticFit quadratic = weightedQuadraticFit(
                        originalMeasurements, std::max(4U, minimumPeakCount_),
                        std::max(700.0, minimumEnergySpanKeV_));
                    const double curvatureKeV = quadratic.valid
                        ? std::abs(quadratic.coefficients[2])*0.25*
                          config::kDriftEnergyMax*config::kDriftEnergyMax
                        : 0.0;
                    double linearAllSquared = 0.0;
                    double quadraticAllSquared = 0.0;
                    for (const PeakMeasurement& point : originalMeasurements) {
                        const double linearResidual = point.referenceEnergy-
                            (fit.offset+fit.slope*point.measuredEnergy);
                        const double quadraticResidual = point.referenceEnergy-
                            evaluatePolynomial(quadratic.coefficients,
                                               point.measuredEnergy);
                        linearAllSquared += linearResidual*linearResidual;
                        quadraticAllSquared += quadraticResidual*
                                               quadraticResidual;
                    }
                    const double allCount = static_cast<double>(
                        originalMeasurements.size());
                    const double linearAllRms = std::sqrt(
                        linearAllSquared/std::max(1.0, allCount));
                    const double quadraticAllRms = std::sqrt(
                        quadraticAllSquared/std::max(1.0, allCount));
                    if (quadratic.valid &&
                        curvatureKeV >= 0.4 &&
                        quadraticAllRms+0.15 < linearAllRms &&
                        quadraticAllRms <= 0.8*linearAllRms) {
                        polynomial.coefficients = quadratic.coefficients;
                        polynomial.rms = quadratic.rms;
                        polynomial.peakCount = quadratic.usedPoints;
                        polynomial.nonlinear = true;
                    }
                }
            }
        }
        std::cout << "Fitted drift group " << group+1 << " / "
                  << groupCount << ".\r" << std::flush;
    }
    std::cout << "Fitted " << groupCount << " drift groups.          \n";

    // DriftCorrections is the authoritative record of the robust affine fit;
    // DriftPeakFits retains every independently accepted or rejected peak.
    for (auto& detectorCorrections : corrections) {
        fillMissingCorrections(detectorCorrections, chronologicalGroups);
    }
    fillMissingPolynomialCorrections(detector16Corrections,
                                     chronologicalGroups);

    std::vector<ShapePeakWindowDefinition> detector4WindowDefinitions;
    for (const PeakDefinition& peak : peaks_) {
        if (peak.energyKeV < 450.0 || peak.energyKeV > 750.0) continue;
        detector4WindowDefinitions.push_back({
            peak.energyKeV, std::max(10.0, peak.searchHalfWidthKeV),
            std::max(4.0, peak.fitHalfWidthKeV)});
    }
    const std::vector<WindowedShapePeakRecord> detector4WindowedPeaks =
        fitWindowedDetector4ShapePeaks(
            inputFile, corrections, files, groupCount,
            detector4WindowDefinitions);
    const std::size_t acceptedDetector4WindowedPeaks =
        static_cast<std::size_t>(std::count_if(
            detector4WindowedPeaks.begin(), detector4WindowedPeaks.end(),
            [](const WindowedShapePeakRecord& record) {
                return record.accepted;
            }));
    std::cout << "Detector 4: accepted " << acceptedDetector4WindowedPeaks
              << " / " << detector4WindowedPeaks.size()
              << " sliding-window 536/628-keV shape-anchor fits.\n";
    std::vector<AcceptedPeakRecord> detector4CalibrationRecords;
    detector4CalibrationRecords.reserve(
        acceptedPeakRecords.size()+detector4WindowedPeaks.size());
    for (const AcceptedPeakRecord& record : acceptedPeakRecords) {
        const bool replacedShapeAnchor =
            record.detectorID == kNonlinearDriftDetectorID &&
            record.referenceEnergy >= 450.0 &&
            record.referenceEnergy <= 750.0;
        if (!replacedShapeAnchor) detector4CalibrationRecords.push_back(record);
    }
    for (const WindowedShapePeakRecord& record : detector4WindowedPeaks) {
        if (record.accepted) {
            detector4CalibrationRecords.push_back(record.peak);
        }
    }
    std::vector<PiecewiseLinearCorrection> detector16Piecewise(groupCount);
    for (std::size_t group = 0; group < groupCount; ++group) {
        std::vector<std::pair<double, double>> anchors;
        double affineSquared = 0.0;
        for (const AcceptedPeakRecord& record : acceptedPeakRecords) {
            if (record.groupIndex != group ||
                record.detectorID != kQuadraticDriftDetectorID) continue;
            anchors.emplace_back(record.measuredEnergy,
                                 record.referenceEnergy);
            const GroupCorrection& affine =
                corrections[kQuadraticDriftDetectorID][group];
            const double residual = affine.offset+
                affine.slope*record.measuredEnergy-record.referenceEnergy;
            affineSquared += residual*residual;
        }
        std::sort(anchors.begin(), anchors.end());
        if (anchors.size() < 4U ||
            anchors.back().second-anchors.front().second < 700.0) continue;
        const double affineRms = std::sqrt(
            affineSquared/static_cast<double>(anchors.size()));
        if (affineRms < 0.8) continue;
        bool monotonic = true;
        for (std::size_t index = 0; index+1U < anchors.size(); ++index) {
            const double inputSpan = anchors[index+1U].first-
                                     anchors[index].first;
            const double outputSpan = anchors[index+1U].second-
                                      anchors[index].second;
            const double slope = outputSpan/inputSpan;
            if (!(inputSpan > 0.0 && slope > 0.9 && slope < 1.1)) {
                monotonic = false;
                break;
            }
        }
        if (monotonic) {
            detector16Piecewise[group].valid = true;
            detector16Piecewise[group].anchors = std::move(anchors);
            detector16Corrections[group].nonlinear = true;
        }
    }
    const std::vector<PiecewiseLinearCorrection> detector4Piecewise =
        determineDetector4Corrections(detector4CalibrationRecords, corrections,
                                      chronologicalGroups);
    const std::size_t validDetector4Corrections = static_cast<std::size_t>(
        std::count_if(detector4Piecewise.begin(), detector4Piecewise.end(),
            [](const PiecewiseLinearCorrection& correction) {
                return correction.valid;
            }));
    if (validDetector4Corrections != 0U) {
        std::cout << "Detector 4: made " << validDetector4Corrections
                  << " robust full-range piecewise corrections.\n";
    } else {
        std::cerr << "Warning: detector 4 nonlinear correction could not be "
                     "determined; writing its ordinary affine correction.\n";
    }
    const std::size_t nonlinearDetector16Groups = static_cast<std::size_t>(
        std::count_if(detector16Piecewise.begin(), detector16Piecewise.end(),
            [](const PiecewiseLinearCorrection& correction) {
                return correction.valid;
            }));
    std::cout << "Detector 16: selected a guarded piecewise-linear correction in "
              << nonlinearDetector16Groups << " / " << groupCount
              << " groups.\n";


    unsigned int residualGroupIndex = 0;
    UShort_t residualDetectorID = 0;
    double residualReferenceEnergy = 0.0;
    double residualMeasuredEnergy = 0.0;
    double residualAffineEnergy = 0.0;
    double residualAffineValue = 0.0;
    double residualFinalEnergy = 0.0;
    double residualFinalValue = 0.0;
    bool residualShapeApplied = false;
    bool residualDetector4PiecewiseApplied = false;
    double residualDetector4TrajectoryDifference = 0.0;
    bool residualDetector4TrajectoryOutlier = false;
    bool residualQuadraticApplied = false;
    bool residualPiecewiseApplied = false;
    TTree residualTree(
        "DriftCorrectedPeakResiduals",
        "Accepted peak residuals after affine and final drift corrections");
    residualTree.Branch("groupIndex", &residualGroupIndex);
    residualTree.Branch("detectorID", &residualDetectorID);
    residualTree.Branch("referenceEnergyKeV", &residualReferenceEnergy);
    residualTree.Branch("measuredEnergyKeV", &residualMeasuredEnergy);
    residualTree.Branch("affineEnergyKeV", &residualAffineEnergy);
    residualTree.Branch("affineResidualKeV", &residualAffineValue);
    residualTree.Branch("finalEnergyKeV", &residualFinalEnergy);
    residualTree.Branch("finalResidualKeV", &residualFinalValue);
    residualTree.Branch("shapeApplied", &residualShapeApplied);
    residualTree.Branch("detector4PiecewiseApplied",
                        &residualDetector4PiecewiseApplied);
    residualTree.Branch("detector4TrajectoryDifferenceKeV",
                        &residualDetector4TrajectoryDifference);
    residualTree.Branch("detector4TrajectoryOutlier",
                        &residualDetector4TrajectoryOutlier);
    residualTree.Branch("quadraticApplied", &residualQuadraticApplied);
    residualTree.Branch("piecewiseApplied", &residualPiecewiseApplied);
    for (const AcceptedPeakRecord& record : acceptedPeakRecords) {
        const GroupCorrection& correction =
            corrections[record.detectorID][record.groupIndex];
        if (!correction.valid) continue;
        residualGroupIndex = record.groupIndex;
        residualDetectorID = static_cast<UShort_t>(record.detectorID);
        residualReferenceEnergy = record.referenceEnergy;
        residualMeasuredEnergy = record.measuredEnergy;
        residualAffineEnergy = correction.offset+
            correction.slope*record.measuredEnergy;
        residualAffineValue = residualAffineEnergy-record.referenceEnergy;
        residualShapeApplied = false;
        residualDetector4PiecewiseApplied =
            record.detectorID == kNonlinearDriftDetectorID &&
            detector4Piecewise[record.groupIndex].valid;
        residualDetector4TrajectoryDifference = 0.0;
        residualDetector4TrajectoryOutlier = false;
        if (residualDetector4PiecewiseApplied) {
            for (const auto& anchor :
                 detector4Piecewise[record.groupIndex].anchors) {
                if (std::abs(anchor.second-record.referenceEnergy) >= 0.01) {
                    continue;
                }
                const double trajectoryResidual =
                    anchor.first-anchor.second;
                residualDetector4TrajectoryDifference =
                    residualAffineValue-trajectoryResidual;
                residualDetector4TrajectoryOutlier =
                    std::abs(residualDetector4TrajectoryDifference) > 1.5;
                break;
            }
        }
        residualPiecewiseApplied = record.detectorID ==
                kQuadraticDriftDetectorID &&
            detector16Piecewise[record.groupIndex].valid;
        residualQuadraticApplied = !residualPiecewiseApplied &&
            record.detectorID ==
                kQuadraticDriftDetectorID &&
            detector16Corrections[record.groupIndex].valid &&
            detector16Corrections[record.groupIndex].nonlinear;
        if (residualDetector4PiecewiseApplied) {
            residualFinalEnergy = applyPiecewiseLinear(
                detector4Piecewise[record.groupIndex], residualAffineEnergy);
        } else if (residualPiecewiseApplied) {
            residualFinalEnergy = applyPiecewiseLinear(
                detector16Piecewise[record.groupIndex],
                record.measuredEnergy);
        } else if (residualQuadraticApplied) {
            residualFinalEnergy = evaluatePolynomial(
                detector16Corrections[record.groupIndex].coefficients,
                record.measuredEnergy);
        } else {
            residualFinalEnergy = residualAffineEnergy;
        }
        residualFinalValue = residualFinalEnergy-record.referenceEnergy;
        residualTree.Fill();
    }

    unsigned int correctionGroupIndex = 0;
    UShort_t correctionDetectorID = 0;
    int correctionSource = 0;
    bool correctionValid = false;
    double correctionOffset = 0.0;
    double correctionSlope = 1.0;
    double correctionOffsetError = 0.0;
    double correctionSlopeError = 0.0;
    double correctionRms = 0.0;
    unsigned int correctionPeakCount = 0;
    TTree correctionTree("DriftCorrections",
                         "Measured and interpolated affine drift corrections");
    correctionTree.Branch("groupIndex", &correctionGroupIndex);
    unsigned int correctionChronologicalPosition = 0;
    correctionTree.Branch("chronologicalPosition",
                          &correctionChronologicalPosition);
    correctionTree.Branch("detectorID", &correctionDetectorID);
    correctionTree.Branch("source", &correctionSource);
    correctionTree.Branch("valid", &correctionValid);
    correctionTree.Branch("offsetKeV", &correctionOffset);
    correctionTree.Branch("slope", &correctionSlope);
    correctionTree.Branch("offsetErrorKeV", &correctionOffsetError);
    correctionTree.Branch("slopeError", &correctionSlopeError);
    correctionTree.Branch("rmsKeV", &correctionRms);
    correctionTree.Branch("peakCount", &correctionPeakCount);

    std::array<std::unique_ptr<TGraphErrors>, config::kGermaniumIdBins>
        offsetGraphs;
    std::array<std::unique_ptr<TGraphErrors>, config::kGermaniumIdBins>
        slopeGraphs;
    std::vector<std::size_t> chronologicalPosition(groupCount);
    for (std::size_t position = 0; position < chronologicalGroups.size();
         ++position) {
        chronologicalPosition[chronologicalGroups[position]] = position;
    }
    for (unsigned int detector = 0;
         detector < config::kGermaniumIdBins; ++detector) {
        offsetGraphs[detector] = std::make_unique<TGraphErrors>();
        slopeGraphs[detector] = std::make_unique<TGraphErrors>();
        offsetGraphs[detector]->SetName(
            graphName("g_DriftOffset", detector).c_str());
        slopeGraphs[detector]->SetName(
            graphName("g_DriftSlope", detector).c_str());
        offsetGraphs[detector]->SetTitle(
            ("Germanium ID " + std::to_string(detector) +
             " drift offset;Chronological group;Offset (keV)").c_str());
        slopeGraphs[detector]->SetTitle(
            ("Germanium ID " + std::to_string(detector) +
             " drift slope;Chronological group;Slope").c_str());
        int graphPoint = 0;
        for (unsigned int group = 0; group < groupCount; ++group) {
            const GroupCorrection& correction = corrections[detector][group];
            correctionGroupIndex = group;
            correctionChronologicalPosition = static_cast<unsigned int>(
                chronologicalPosition[group]);
            correctionDetectorID = static_cast<UShort_t>(detector);
            correctionSource = static_cast<int>(correction.source);
            correctionValid = correction.valid;
            correctionOffset = correction.offset;
            correctionSlope = correction.slope;
            correctionOffsetError = correction.offsetError;
            correctionSlopeError = correction.slopeError;
            correctionRms = correction.rms;
            correctionPeakCount = correction.peakCount;
            correctionTree.Fill();
            if (correction.valid) {
                offsetGraphs[detector]->SetPoint(
                    graphPoint, chronologicalPosition[group], correction.offset);
                offsetGraphs[detector]->SetPointError(
                    graphPoint, 0.0, correction.offsetError);
                slopeGraphs[detector]->SetPoint(
                    graphPoint, chronologicalPosition[group], correction.slope);
                slopeGraphs[detector]->SetPointError(
                    graphPoint, 0.0, correction.slopeError);
                ++graphPoint;
            }
        }
    }

    unsigned int nonlinearGroupIndex = 0;
    bool detector4PiecewiseValid = false;
    unsigned int detector4AnchorCount = 0;
    double detector4MinimumAnchorEnergy = 0.0;
    double detector4MaximumAnchorEnergy = 0.0;
    TTree detector4PiecewiseTree(
        "Detector4PiecewiseCorrections",
        "Robust full-range detector-4 piecewise corrections");
    detector4PiecewiseTree.Branch("groupIndex", &nonlinearGroupIndex);
    detector4PiecewiseTree.Branch("valid", &detector4PiecewiseValid);
    detector4PiecewiseTree.Branch("anchorCount", &detector4AnchorCount);
    detector4PiecewiseTree.Branch("minimumAnchorEnergyKeV",
                                  &detector4MinimumAnchorEnergy);
    detector4PiecewiseTree.Branch("maximumAnchorEnergyKeV",
                                  &detector4MaximumAnchorEnergy);

    int detector16Source = 0;
    bool detector16Valid = false;
    bool detector16Nonlinear = false;
    bool detector16PiecewiseApplied = false;
    unsigned int detector16AnchorCount = 0;
    double detector16C0 = 0.0;
    double detector16C1 = 1.0;
    double detector16C2 = 0.0;
    double detector16Rms = 0.0;
    unsigned int detector16PeakCount = 0;
    TTree detector16Tree(
        "Detector16PolynomialCorrections",
        "Guarded detector-16 correction selection");
    detector16Tree.Branch("groupIndex", &nonlinearGroupIndex);
    detector16Tree.Branch("source", &detector16Source);
    detector16Tree.Branch("valid", &detector16Valid);
    detector16Tree.Branch("nonlinear", &detector16Nonlinear);
    detector16Tree.Branch("piecewise", &detector16PiecewiseApplied);
    detector16Tree.Branch("piecewiseAnchorCount", &detector16AnchorCount);
    detector16Tree.Branch("c0", &detector16C0);
    detector16Tree.Branch("c1", &detector16C1);
    detector16Tree.Branch("c2", &detector16C2);
    detector16Tree.Branch("rmsKeV", &detector16Rms);
    detector16Tree.Branch("peakCount", &detector16PeakCount);
    for (std::size_t group = 0; group < groupCount; ++group) {
        nonlinearGroupIndex = static_cast<unsigned int>(group);
        const PiecewiseLinearCorrection& detector4 =
            detector4Piecewise[group];
        detector4PiecewiseValid = detector4.valid;
        detector4AnchorCount = static_cast<unsigned int>(
            detector4.anchors.size());
        detector4MinimumAnchorEnergy = detector4.valid
            ? detector4.anchors.front().second : 0.0;
        detector4MaximumAnchorEnergy = detector4.valid
            ? detector4.anchors.back().second : 0.0;
        detector4PiecewiseTree.Fill();

        const PolynomialCorrection& polynomial = detector16Corrections[group];
        detector16Source = static_cast<int>(polynomial.source);
        detector16Valid = polynomial.valid;
        detector16Nonlinear = polynomial.nonlinear;
        detector16PiecewiseApplied = detector16Piecewise[group].valid;
        detector16AnchorCount = static_cast<unsigned int>(
            detector16Piecewise[group].anchors.size());
        detector16C0 = polynomial.coefficients[0];
        detector16C1 = polynomial.coefficients[1];
        detector16C2 = polynomial.coefficients[2];
        detector16Rms = polynomial.rms;
        detector16PeakCount = polynomial.peakCount;
        detector16Tree.Fill();
    }

    std::map<std::string, std::vector<FileRecord>> filesBySeries;
    for (const FileRecord& file : files) filesBySeries[file.series].push_back(file);
    std::vector<std::string> writtenFiles;
    for (auto& seriesEntry : filesBySeries) {
        std::sort(seriesEntry.second.begin(), seriesEntry.second.end(),
                  [](const FileRecord& left, const FileRecord& right) {
                      return left.run < right.run;
                  });
        const std::string fileName =
            correctionFileName(outputPrefix, seriesEntry.first);
        const std::filesystem::path parent =
            std::filesystem::path(fileName).parent_path();
        if (!parent.empty()) std::filesystem::create_directories(parent);
        std::ofstream output(fileName);
        if (!output) {
            throw std::runtime_error("could not create correction file '" +
                                     fileName + "'");
        }
        output << "# Fine drift correction for series " << seriesEntry.first
               << "\n# Applied after the main calibration. Detector 16 may "
                  "use a guarded quadratic; other entries are affine.\n"
               << "# RUN ID NCOEFF COEFFICIENTS...\n";
        output << std::setprecision(12);
        for (const FileRecord& file : seriesEntry.second) {
            for (unsigned int detector = 0;
                 detector < config::kGermaniumIdBins; ++detector) {
                if (detector4Piecewise[file.groupIndex].valid &&
                    detector == kNonlinearDriftDetectorID) {
                    continue;
                }
                if (detector == kQuadraticDriftDetectorID &&
                    detector16Piecewise[file.groupIndex].valid) {
                    continue;
                }
                if (detector == kQuadraticDriftDetectorID &&
                    detector16Corrections[file.groupIndex].valid) {
                    const PolynomialCorrection& polynomial =
                        detector16Corrections[file.groupIndex];
                    output << file.run << ' ' << detector << ' '
                           << (polynomial.nonlinear ? 3 : 2) << ' '
                           << polynomial.coefficients[0] << ' '
                           << polynomial.coefficients[1];
                    if (polynomial.nonlinear) {
                        output << ' ' << polynomial.coefficients[2];
                    }
                    output << '\n';
                    continue;
                }
                const GroupCorrection& correction =
                    corrections[detector][file.groupIndex];
                if (!correction.valid) continue;
                output << file.run << ' ' << detector << " 2 "
                       << correction.offset << ' ' << correction.slope << '\n';
            }
        }
        writtenFiles.push_back(fileName);
        std::cout << "Wrote " << fileName << ".\n";

        if (validDetector4Corrections != 0U ||
            nonlinearDetector16Groups != 0U) {
            const std::string nonlinearFileName =
                nonlinearCorrectionFileName(outputPrefix, seriesEntry.first);
            std::ofstream nonlinearOutput(nonlinearFileName);
            if (!nonlinearOutput) {
                throw std::runtime_error(
                    "could not create correction file '" +
                    nonlinearFileName + "'");
            }
            nonlinearOutput
                << "# Nonlinear fine drift corrections (detector 4 and, when "
                   "needed, detector 16) for "
                << "series " << seriesEntry.first << '\n'
                << "# Apply after the main calibration and the matching "
                   "run-by-run .cal file.\n"
                << "# RUN ID NPIECES "
                   "[NCOEFF COEFFICIENTS... INPUT_END_KEV]...\n";
            nonlinearOutput << std::setprecision(16);
            for (const FileRecord& file : seriesEntry.second) {
                const PiecewiseLinearCorrection& detector4 =
                    detector4Piecewise[file.groupIndex];
                if (detector4.valid) {
                    const GroupCorrection& correction =
                        corrections[kNonlinearDriftDetectorID][file.groupIndex];
                    if (correction.valid) {
                        nonlinearOutput << file.run << ' '
                            << kNonlinearDriftDetectorID << ' '
                            << detector4.anchors.size()-1U;
                        for (std::size_t index = 0;
                             index+1U < detector4.anchors.size(); ++index) {
                            const auto& first = detector4.anchors[index];
                            const auto& second = detector4.anchors[index+1U];
                            const double pieceSlope =
                                (second.second-first.second)/
                                (second.first-first.first);
                            const double pieceOffset =
                                first.second-pieceSlope*first.first;
                            const double composedOffset = pieceOffset+
                                pieceSlope*correction.offset;
                            const double composedSlope =
                                pieceSlope*correction.slope;
                            const double end =
                                index+2U < detector4.anchors.size()
                                    ? (second.first-correction.offset)/
                                      correction.slope
                                    : config::kDriftEnergyMax;
                            if (!(end > 0.0 &&
                                  end <= config::kDriftEnergyMax)) {
                                throw std::runtime_error(
                                    "detector 4 piecewise correction has invalid "
                                    "boundary for run "+
                                    std::to_string(file.run));
                            }
                            nonlinearOutput << " 2 " << composedOffset << ' '
                                            << composedSlope << ' ' << end;
                        }
                        nonlinearOutput << '\n';
                    }
                }

                const PiecewiseLinearCorrection& piecewise =
                    detector16Piecewise[file.groupIndex];
                if (!piecewise.valid) continue;
                nonlinearOutput << file.run << ' '
                    << kQuadraticDriftDetectorID << ' '
                    << piecewise.anchors.size()-1U;
                for (std::size_t index = 0;
                     index+1U < piecewise.anchors.size(); ++index) {
                    const auto& first = piecewise.anchors[index];
                    const auto& second = piecewise.anchors[index+1U];
                    const double slope = (second.second-first.second)/
                                         (second.first-first.first);
                    const double offset = first.second-slope*first.first;
                    const double end = index+2U < piecewise.anchors.size()
                        ? second.first : config::kDriftEnergyMax;
                    nonlinearOutput << " 2 " << offset << ' ' << slope
                                    << ' ' << end;
                }
                nonlinearOutput << '\n';
            }
            writtenFiles.push_back(nonlinearFileName);
            std::cout << "Wrote " << nonlinearFileName << ".\n";
        }
    }

    diagnostics.cd();
    unsigned int windowGroupIndex = 0;
    double windowReferenceEnergy = 0.0;
    double windowCentroidEnergy = 0.0;
    double windowCentroidError = 0.0;
    double windowSigma = 0.0;
    double windowSignificance = 0.0;
    unsigned int windowNeighbourCount = 0;
    bool windowAccepted = false;
    TTree windowedShapeTree(
        "Detector4WindowedPeakFits",
        "Higher-statistics detector-4 shape-anchor fits");
    windowedShapeTree.Branch("groupIndex", &windowGroupIndex);
    windowedShapeTree.Branch("referenceEnergyKeV", &windowReferenceEnergy);
    windowedShapeTree.Branch("centroidKeV", &windowCentroidEnergy);
    windowedShapeTree.Branch("centroidErrorKeV", &windowCentroidError);
    windowedShapeTree.Branch("sigmaKeV", &windowSigma);
    windowedShapeTree.Branch("significance", &windowSignificance);
    windowedShapeTree.Branch("neighbourCount", &windowNeighbourCount);
    windowedShapeTree.Branch("accepted", &windowAccepted);
    for (const WindowedShapePeakRecord& record : detector4WindowedPeaks) {
        windowGroupIndex = record.peak.groupIndex;
        windowReferenceEnergy = record.peak.referenceEnergy;
        const GroupCorrection& affine =
            corrections[kNonlinearDriftDetectorID][record.peak.groupIndex];
        windowCentroidEnergy = record.accepted
            ? affine.offset+affine.slope*record.peak.measuredEnergy : 0.0;
        windowCentroidError = record.centroidErrorKeV;
        windowSigma = record.sigmaKeV;
        windowSignificance = record.significance;
        windowNeighbourCount = record.neighbourCount;
        windowAccepted = record.accepted;
        windowedShapeTree.Fill();
    }

    unsigned int metadataGroupIndex = 0;
    unsigned int metadataChronologicalPosition = 0;
    std::string metadataSeries;
    unsigned int metadataFirstRun = 0;
    unsigned int metadataLastRun = 0;
    TTree groupMetadataTree(
        "DriftGroupMetadata", "Series and run range for each drift group");
    groupMetadataTree.Branch("groupIndex", &metadataGroupIndex);
    groupMetadataTree.Branch("chronologicalPosition",
                             &metadataChronologicalPosition);
    groupMetadataTree.Branch("series", &metadataSeries);
    groupMetadataTree.Branch("firstRun", &metadataFirstRun);
    groupMetadataTree.Branch("lastRun", &metadataLastRun);
    for (std::size_t group = 0; group < groupCount; ++group) {
        metadataGroupIndex = static_cast<unsigned int>(group);
        metadataChronologicalPosition = static_cast<unsigned int>(
            chronologicalPosition[group]);
        metadataSeries.clear();
        metadataFirstRun = std::numeric_limits<unsigned int>::max();
        metadataLastRun = 0;
        for (const FileRecord& file : files) {
            if (file.groupIndex != group) continue;
            metadataSeries = file.series;
            metadataFirstRun = std::min(metadataFirstRun, file.run);
            metadataLastRun = std::max(metadataLastRun, file.run);
        }
        if (metadataFirstRun == std::numeric_limits<unsigned int>::max()) {
            metadataFirstRun = 0;
        }
        groupMetadataTree.Fill();
    }

    peakTree.Write();
    alignmentTree.Write();
    correctionTree.Write();
    residualTree.Write();
    windowedShapeTree.Write();
    groupMetadataTree.Write();
    detector4PiecewiseTree.Write();
    detector16Tree.Write();
    TDirectory* graphDirectory = diagnostics.mkdir("CorrectionGraphs");
    graphDirectory->cd();
    for (unsigned int detector = 0;
         detector < config::kGermaniumIdBins; ++detector) {
        offsetGraphs[detector]->Write();
        slopeGraphs[detector]->Write();
    }
    diagnostics.cd();
    TNamed inputMetadata("DriftCalibrationInput", inputFileName.c_str());
    inputMetadata.Write();
    TNamed sourceKey(
        "DriftCorrectionSourceCodes",
        "0=measured, 1=interpolated between measured groups, "
        "2=nearest measured group, 3=missing");
    sourceKey.Write();
    if (validDetector4Corrections != 0U) {
        std::ostringstream description;
        description << std::setprecision(12)
            << "Detector 4 uses stable affine anchors outside 450--750 keV. "
               "The 536- and 628-keV shape anchors are fitted after affine "
               "alignment in sliding windows of up to five groups from the "
               "same series. Reference residuals from approximately 301 to "
               "1779 keV are blended with robust seven-group trajectories. "
               "Missing and outlying centroids use the local trajectory. "
               "The resulting monotonic piecewise-linear mappings are stored "
               "in Detector4PiecewiseCorrections.";
        TNamed detector4Metadata(
            "Detector4NonlinearDriftCorrection", description.str().c_str());
        detector4Metadata.Write();
    }
    std::ostringstream orderDescription;
    for (std::size_t index = 0; index < seriesOrder_.size(); ++index) {
        if (index != 0U) orderDescription << ',';
        orderDescription << seriesOrder_[index];
    }
    TNamed orderMetadata(
        "DriftChronologicalSeriesOrder",
        (seriesOrder_.empty() ? "input provenance order"
                              : orderDescription.str()).c_str());
    orderMetadata.Write();
    diagnostics.Close();
    std::cout << "Wrote diagnostics to " << diagnosticFileName << ".\n";
    return 0;
}
