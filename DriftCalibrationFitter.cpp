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

constexpr unsigned int kNonlinearDriftDetectorID = 4;
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

bool solveFourByFour(std::array<std::array<long double, 5>, 4>& matrix,
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

ShapeCorrection determineShapeCorrection(
    const std::vector<AcceptedPeakRecord>& records,
    const std::array<std::vector<GroupCorrection>,
                     config::kGermaniumIdBins>& corrections)
{
    std::vector<double> firstResiduals;
    std::vector<double> secondResiduals;
    for (const AcceptedPeakRecord& record : records) {
        if (record.detectorID != kNonlinearDriftDetectorID ||
            record.groupIndex >= corrections[record.detectorID].size()) {
            continue;
        }
        const GroupCorrection& correction =
            corrections[record.detectorID][record.groupIndex];
        if (!correction.valid) continue;
        const double residual = correction.offset +
            correction.slope*record.measuredEnergy-record.referenceEnergy;
        if (std::abs(record.referenceEnergy-kFirstShapeAnchorKeV) < 0.1) {
            firstResiduals.push_back(residual);
        } else if (std::abs(record.referenceEnergy-
                            kSecondShapeAnchorKeV) < 0.1) {
            secondResiduals.push_back(residual);
        }
    }

    ShapeCorrection result;
    if (firstResiduals.size() < 20U || secondResiduals.size() < 20U) {
        return result;
    }
    result.firstAnchorResidual = median(std::move(firstResiduals));
    result.secondAnchorResidual = median(std::move(secondResiduals));
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

std::array<double, 6> composeShapePolynomial(
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

double applyShapeCorrection(const ShapeCorrection& shape, double energy)
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

void fillMissingCorrections(std::vector<GroupCorrection>& corrections)
{
    std::vector<std::size_t> validIndices;
    for (std::size_t index = 0; index < corrections.size(); ++index) {
        if (corrections[index].valid) validIndices.push_back(index);
    }
    if (validIndices.empty()) return;

    for (std::size_t index = 0; index < corrections.size(); ++index) {
        if (corrections[index].valid) continue;
        const auto after = std::lower_bound(
            validIndices.begin(), validIndices.end(), index);
        if (after == validIndices.begin()) {
            corrections[index] = corrections[*after];
            corrections[index].source = CorrectionSource::Nearest;
        } else if (after == validIndices.end()) {
            corrections[index] = corrections[validIndices.back()];
            corrections[index].source = CorrectionSource::Nearest;
        } else {
            const std::size_t high = *after;
            const std::size_t low = *(after-1);
            const double fraction = static_cast<double>(index-low)/
                                    static_cast<double>(high-low);
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
                fitAccepted = fitStatus == 0 && std::isfinite(fitCentroid) &&
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
            LinearFit fit = weightedLinearFit(
                measurements, minimumPoints, minimumSpan);
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

            }
        }
        std::cout << "Fitted drift group " << group+1 << " / "
                  << groupCount << ".\r" << std::flush;
    }
    std::cout << "Fitted " << groupCount << " drift groups.          \n";

    // DriftCorrections is the authoritative record of the robust affine fit;
    // DriftPeakFits retains every independently accepted or rejected peak.
    for (auto& detectorCorrections : corrections) {
        fillMissingCorrections(detectorCorrections);
    }
    const ShapeCorrection detector4Shape = determineShapeCorrection(
        acceptedPeakRecords, corrections);
    if (detector4Shape.valid) {
        std::cout << "Detector 4: separated stable affine drift from the "
                     "repeatable 536/628-keV response shape (median "
                     "pre-shape residuals "
                  << detector4Shape.firstAnchorResidual << " and "
                  << detector4Shape.secondAnchorResidual << " keV).\n";
    } else {
        std::cerr << "Warning: detector 4 nonlinear correction could not be "
                     "determined; writing its ordinary affine correction.\n";
    }


    unsigned int residualGroupIndex = 0;
    UShort_t residualDetectorID = 0;
    double residualReferenceEnergy = 0.0;
    double residualMeasuredEnergy = 0.0;
    double residualAffineEnergy = 0.0;
    double residualAffineValue = 0.0;
    double residualFinalEnergy = 0.0;
    double residualFinalValue = 0.0;
    bool residualShapeApplied = false;
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
        residualShapeApplied = detector4Shape.valid &&
            record.detectorID == kNonlinearDriftDetectorID &&
            residualAffineEnergy > kShapeLowerEnergyKeV &&
            residualAffineEnergy < kShapeUpperEnergyKeV;
        residualFinalEnergy = residualShapeApplied
            ? applyShapeCorrection(detector4Shape, residualAffineEnergy)
            : residualAffineEnergy;
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
                    graphPoint, group, correction.offset);
                offsetGraphs[detector]->SetPointError(
                    graphPoint, 0.0, correction.offsetError);
                slopeGraphs[detector]->SetPoint(
                    graphPoint, group, correction.slope);
                slopeGraphs[detector]->SetPointError(
                    graphPoint, 0.0, correction.slopeError);
                ++graphPoint;
            }
        }
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
               << "\n# Applied after the main calibration: "
                  "E_corrected = C0 + C1*E_main\n"
               << "# RUN ID NCOEFF C0 C1\n";
        output << std::setprecision(12);
        for (const FileRecord& file : seriesEntry.second) {
            for (unsigned int detector = 0;
                 detector < config::kGermaniumIdBins; ++detector) {
                if (detector4Shape.valid &&
                    detector == kNonlinearDriftDetectorID) {
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

        if (detector4Shape.valid) {
            const std::string nonlinearFileName =
                nonlinearCorrectionFileName(outputPrefix, seriesEntry.first);
            std::ofstream nonlinearOutput(nonlinearFileName);
            if (!nonlinearOutput) {
                throw std::runtime_error(
                    "could not create correction file '" +
                    nonlinearFileName + "'");
            }
            nonlinearOutput
                << "# Detector 4 fine drift and static shape correction for "
                << "series " << seriesEntry.first << '\n'
                << "# Apply after the main calibration and the matching "
                   "run-by-run .cal file.\n"
                << "# RUN ID NPIECES "
                   "[NCOEFF COEFFICIENTS... INPUT_END_KEV]...\n";
            nonlinearOutput << std::setprecision(16);
            for (const FileRecord& file : seriesEntry.second) {
                const GroupCorrection& correction =
                    corrections[kNonlinearDriftDetectorID][file.groupIndex];
                if (!correction.valid) continue;
                const double lowerInputBoundary =
                    (kShapeLowerEnergyKeV-correction.offset)/correction.slope;
                const double upperInputBoundary =
                    (kShapeUpperEnergyKeV-correction.offset)/correction.slope;
                if (!(lowerInputBoundary > 0.0 &&
                      upperInputBoundary > lowerInputBoundary &&
                      upperInputBoundary < config::kDriftEnergyMax)) {
                    throw std::runtime_error(
                        "detector 4 nonlinear correction has invalid piece "
                        "boundaries for run " + std::to_string(file.run));
                }
                const std::array<double, 6> middle =
                    composeShapePolynomial(detector4Shape, correction);
                nonlinearOutput << file.run << ' '
                    << kNonlinearDriftDetectorID << " 3 2 "
                    << correction.offset << ' ' << correction.slope << ' '
                    << lowerInputBoundary << " 6";
                for (const double coefficient : middle) {
                    nonlinearOutput << ' ' << coefficient;
                }
                nonlinearOutput << ' ' << upperInputBoundary << " 2 "
                    << correction.offset << ' ' << correction.slope << ' '
                    << config::kDriftEnergyMax << '\n';
            }
            writtenFiles.push_back(nonlinearFileName);
            std::cout << "Wrote " << nonlinearFileName << ".\n";
        }
    }

    diagnostics.cd();
    peakTree.Write();
    alignmentTree.Write();
    correctionTree.Write();
    residualTree.Write();
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
    if (detector4Shape.valid) {
        std::ostringstream description;
        description << std::setprecision(12)
            << "Detector 4 uses stable affine anchors outside 450--750 keV. "
               "A static smooth correction is zero in value and slope at "
            << kShapeLowerEnergyKeV << " and " << kShapeUpperEnergyKeV
            << " keV; median affine residuals at " << kFirstShapeAnchorKeV
            << " and " << kSecondShapeAnchorKeV << " keV were "
            << detector4Shape.firstAnchorResidual << " and "
            << detector4Shape.secondAnchorResidual
            << " keV. Normalized delta coefficients c0..c5: ";
        for (const double coefficient :
             detector4Shape.normalizedCoefficients) {
            description << coefficient << ' ';
        }
        TNamed detector4Metadata(
            "Detector4NonlinearDriftCorrection", description.str().c_str());
        detector4Metadata.Write();
    }
    diagnostics.Close();
    std::cout << "Wrote diagnostics to " << diagnosticFileName << ".\n";
    return 0;
}
