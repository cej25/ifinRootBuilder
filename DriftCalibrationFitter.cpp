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
                    measurements.push_back({fitCentroid, fitCentroidError,
                                            peak.energyKeV, true});
                }
                peakTree.Fill();
            }

            LinearFit fit = weightedLinearFit(
                measurements, minimumPeakCount_, minimumEnergySpanKeV_);
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
                const GroupCorrection& correction =
                    corrections[detector][file.groupIndex];
                if (!correction.valid) continue;
                output << file.run << ' ' << detector << " 2 "
                       << correction.offset << ' ' << correction.slope << '\n';
            }
        }
        writtenFiles.push_back(fileName);
        std::cout << "Wrote " << fileName << ".\n";
    }

    diagnostics.cd();
    peakTree.Write();
    alignmentTree.Write();
    correctionTree.Write();
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
    diagnostics.Close();
    std::cout << "Wrote diagnostics to " << diagnosticFileName << ".\n";
    return 0;
}
