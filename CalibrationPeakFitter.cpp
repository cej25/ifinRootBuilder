#include "CalibrationPeakFitter.h"

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

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
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

std::string spectrumPath(const std::string& group, unsigned int detectorID)
{
    std::ostringstream path;
    path << group << "/Raw/Individual/h1_Ge_ADC_ID"
         << std::setfill('0') << std::setw(2) << detectorID;
    return path.str();
}

std::string graphName(unsigned int detectorID)
{
    std::ostringstream name;
    name << "g_Ge_CalibrationPoints_ID" << std::setfill('0')
         << std::setw(2) << detectorID;
    return name.str();
}

struct RawInterval {
    double closestChannel = 0.0;
    double closestDifferenceKeV = std::numeric_limits<double>::infinity();
    double minimumChannel = 0.0;
    double maximumChannel = 0.0;
    bool found = false;
};

RawInterval locateEnergy(const GermaniumCalibration& calibration,
                         unsigned int detectorID,
                         double energyKeV,
                         double halfWidthKeV)
{
    RawInterval interval;
    const double lowEnergy = energyKeV - halfWidthKeV;
    const double highEnergy = energyKeV + halfWidthKeV;

    for (unsigned int channel = 0;
         channel < config::kIndividualGermaniumRawBins; ++channel) {
        const auto result = calibration.calibrate(
            static_cast<UShort_t>(detectorID), static_cast<double>(channel));
        if (!result.valid()) continue;

        const double difference = std::abs(result.energy - energyKeV);
        if (difference < interval.closestDifferenceKeV) {
            interval.closestDifferenceKeV = difference;
            interval.closestChannel = static_cast<double>(channel);
        }
        if (result.energy >= lowEnergy && result.energy <= highEnergy) {
            if (!interval.found) {
                interval.minimumChannel = static_cast<double>(channel);
                interval.maximumChannel = static_cast<double>(channel);
                interval.found = true;
            } else {
                interval.minimumChannel = std::min(
                    interval.minimumChannel, static_cast<double>(channel));
                interval.maximumChannel = std::max(
                    interval.maximumChannel, static_cast<double>(channel));
            }
        }
    }
    return interval;
}

double localGain(const GermaniumCalibration& calibration,
                 unsigned int detectorID,
                 double channel)
{
    const double lowChannel = std::max(0.0, channel - 1.0);
    const double highChannel = std::min(
        static_cast<double>(config::kIndividualGermaniumRawBins - 1),
        channel + 1.0);
    const auto low = calibration.calibrate(
        static_cast<UShort_t>(detectorID), lowChannel);
    const auto high = calibration.calibrate(
        static_cast<UShort_t>(detectorID), highChannel);
    if (!low.valid() || !high.valid() || highChannel == lowChannel) return 0.0;
    return (high.energy - low.energy) / (highChannel - lowChannel);
}

} // namespace

void CalibrationPeakFitter::addSeedCalFile(const std::string& fileName)
{
    seedCalibration_.addCalFile(fileName);
}

void CalibrationPeakFitter::addSeedMcalFile(const std::string& fileName)
{
    seedCalibration_.addMcalFile(fileName);
}

void CalibrationPeakFitter::loadPeakDefinitions(const std::string& fileName)
{
    std::ifstream file(fileName);
    if (!file) {
        throw std::runtime_error("Could not open peak-definition file '" +
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
        if (!(input >> peak.group >> peak.energyKeV >> peak.label)) {
            throw parseError(fileName, lineNumber,
                             "expected GROUP ENERGY_KEV LABEL");
        }
        if (input >> peak.searchHalfWidthKeV) {
            if (!(input >> peak.fitHalfWidthKeV)) {
                throw parseError(fileName, lineNumber,
                                 "SEARCH_KEV requires FIT_KEV");
            }
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
        throw std::runtime_error("Peak-definition file contains no peaks");
    }
}

int CalibrationPeakFitter::run(const std::string& inputFileName,
                               const std::string& outputFileName) const
{
    if (seedCalibration_.empty()) {
        std::cerr << "Error: at least one seed calibration stage is required.\n";
        return 2;
    }
    if (peaks_.empty()) {
        std::cerr << "Error: no peak definitions were loaded.\n";
        return 3;
    }

    TFile inputFile(inputFileName.c_str(), "READ");
    if (inputFile.IsZombie()) {
        std::cerr << "Error: could not open '" << inputFileName << "'.\n";
        return 4;
    }
    TFile outputFile(outputFileName.c_str(), "RECREATE");
    if (outputFile.IsZombie()) {
        std::cerr << "Error: could not create '" << outputFileName << "'.\n";
        return 5;
    }

    UShort_t detectorID = 0;
    std::string group;
    std::string label;
    double referenceEnergyKeV = 0.0;
    double predictedAdc = 0.0;
    double centroidAdc = 0.0;
    double centroidErrorAdc = 0.0;
    double sigmaAdc = 0.0;
    double sigmaErrorAdc = 0.0;
    double amplitude = 0.0;
    double amplitudeError = 0.0;
    double background = 0.0;
    double chiSquare = 0.0;
    int ndf = 0;
    int fitStatus = -1;
    bool accepted = false;

    TTree fitTree("PeakFits", "Germanium source-peak fit results");
    fitTree.Branch("detectorID", &detectorID);
    fitTree.Branch("group", &group);
    fitTree.Branch("label", &label);
    fitTree.Branch("referenceEnergyKeV", &referenceEnergyKeV);
    fitTree.Branch("predictedAdc", &predictedAdc);
    fitTree.Branch("centroidAdc", &centroidAdc);
    fitTree.Branch("centroidErrorAdc", &centroidErrorAdc);
    fitTree.Branch("sigmaAdc", &sigmaAdc);
    fitTree.Branch("sigmaErrorAdc", &sigmaErrorAdc);
    fitTree.Branch("amplitude", &amplitude);
    fitTree.Branch("amplitudeError", &amplitudeError);
    fitTree.Branch("background", &background);
    fitTree.Branch("chiSquare", &chiSquare);
    fitTree.Branch("ndf", &ndf);
    fitTree.Branch("fitStatus", &fitStatus);
    fitTree.Branch("accepted", &accepted);

    std::array<std::unique_ptr<TGraphErrors>, config::kGermaniumIdBins> graphs;
    std::array<unsigned int, config::kGermaniumIdBins> graphPoints{};
    for (unsigned int id = 0; id < graphs.size(); ++id) {
        const std::string name = graphName(id);
        graphs[id] = std::make_unique<TGraphErrors>();
        graphs[id]->SetName(name.c_str());
        graphs[id]->SetTitle(("Germanium ID " + std::to_string(id) +
                              " fitted calibration points;ADC channel;"
                              "Reference energy (keV)").c_str());
    }

    std::size_t attemptedFits = 0;
    std::size_t acceptedFits = 0;
    for (unsigned int id = 0; id < config::kGermaniumIdBins; ++id) {
        detectorID = static_cast<UShort_t>(id);
        for (const PeakDefinition& peak : peaks_) {
            ++attemptedFits;
            group = peak.group;
            label = peak.label;
            referenceEnergyKeV = peak.energyKeV;
            predictedAdc = centroidAdc = centroidErrorAdc = 0.0;
            sigmaAdc = sigmaErrorAdc = amplitude = amplitudeError = 0.0;
            background = chiSquare = 0.0;
            ndf = 0;
            fitStatus = -1;
            accepted = false;

            const std::string path = spectrumPath(group, id);
            TH1D* histogram = dynamic_cast<TH1D*>(inputFile.Get(path.c_str()));
            if (histogram == nullptr || histogram->GetEntries() == 0.0) {
                fitStatus = -2;
                fitTree.Fill();
                continue;
            }

            const RawInterval search = locateEnergy(
                seedCalibration_, id, peak.energyKeV,
                peak.searchHalfWidthKeV);
            const RawInterval fitRange = locateEnergy(
                seedCalibration_, id, peak.energyKeV,
                peak.fitHalfWidthKeV);
            predictedAdc = search.closestChannel;
            if (!search.found || !fitRange.found ||
                search.closestDifferenceKeV > peak.searchHalfWidthKeV) {
                fitStatus = -3;
                fitTree.Fill();
                continue;
            }

            const int firstBin = histogram->GetXaxis()->FindBin(
                search.minimumChannel);
            const int lastBin = histogram->GetXaxis()->FindBin(
                search.maximumChannel);
            int maximumBin = firstBin;
            for (int bin = firstBin + 1; bin <= lastBin; ++bin) {
                if (histogram->GetBinContent(bin) >
                    histogram->GetBinContent(maximumBin)) {
                    maximumBin = bin;
                }
            }
            const double peakChannel = histogram->GetBinCenter(maximumBin);
            const double fitMinimum = fitRange.minimumChannel;
            const double fitMaximum = fitRange.maximumChannel;
            if (peakChannel < fitMinimum || peakChannel > fitMaximum ||
                fitMaximum - fitMinimum < 4.0) {
                fitStatus = -4;
                fitTree.Fill();
                continue;
            }

            const int lowBin = histogram->GetXaxis()->FindBin(fitMinimum);
            const int highBin = histogram->GetXaxis()->FindBin(fitMaximum);
            const double edgeBackground = 0.5 *
                (histogram->GetBinContent(lowBin) +
                 histogram->GetBinContent(highBin));
            const double peakHeight = histogram->GetBinContent(maximumBin);
            const double gain = std::abs(localGain(
                seedCalibration_, id, peakChannel));
            const double initialSigma = std::clamp(
                gain > 0.0 ? 0.9 / gain : 2.0, 0.5, 15.0);

            TF1 fitFunction(
                "peakFit",
                "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]+[4]*(x-[1])",
                fitMinimum, fitMaximum);
            fitFunction.SetParameters(
                std::max(1.0, peakHeight - edgeBackground), peakChannel,
                initialSigma, std::max(0.0, edgeBackground), 0.0);
            fitFunction.SetParNames(
                "Amplitude", "Centroid", "Sigma", "Background", "Slope");
            fitFunction.SetParLimits(0, 0.0, std::max(1.0, peakHeight*10.0));
            fitFunction.SetParLimits(1, fitMinimum, fitMaximum);
            fitFunction.SetParLimits(
                2, 0.25, std::max(2.0, fitMaximum-fitMinimum));
            fitFunction.SetParLimits(
                3, 0.0, std::max(1.0, peakHeight*10.0));

            const TFitResultPtr result = histogram->Fit(
                &fitFunction, "QRSN", "", fitMinimum, fitMaximum);
            fitStatus = static_cast<int>(result);
            centroidAdc = fitFunction.GetParameter(1);
            centroidErrorAdc = fitFunction.GetParError(1);
            sigmaAdc = fitFunction.GetParameter(2);
            sigmaErrorAdc = fitFunction.GetParError(2);
            amplitude = fitFunction.GetParameter(0);
            amplitudeError = fitFunction.GetParError(0);
            background = fitFunction.GetParameter(3);
            chiSquare = fitFunction.GetChisquare();
            ndf = fitFunction.GetNDF();

            const double significance = amplitude /
                std::sqrt(std::max(1.0, background));
            accepted = fitStatus == 0 && std::isfinite(centroidAdc) &&
                std::isfinite(centroidErrorAdc) && centroidErrorAdc > 0.0 &&
                centroidAdc >= fitMinimum && centroidAdc <= fitMaximum &&
                sigmaAdc > 0.25 && sigmaAdc < fitMaximum-fitMinimum &&
                significance >= 5.0;
            if (accepted) {
                const int point = static_cast<int>(graphPoints[id]++);
                graphs[id]->SetPoint(point, centroidAdc, referenceEnergyKeV);
                graphs[id]->SetPointError(point, centroidErrorAdc, 0.0);
                ++acceptedFits;
            }
            fitTree.Fill();
        }
        std::cout << "Detector " << std::setw(2) << id << ": "
                  << graphPoints[id] << " accepted peak fits.\n";
    }

    outputFile.cd();
    fitTree.Write();
    TDirectory* graphDirectory = outputFile.mkdir("CalibrationPoints");
    graphDirectory->cd();
    for (const auto& graph : graphs) graph->Write();
    outputFile.cd();
    TNamed inputMetadata("CalibrationPeakFitInput", inputFileName.c_str());
    inputMetadata.Write();
    TNamed seedNote(
        "CalibrationPeakFitSeedNote",
        "Seed calibration was used only to locate fit windows; fitted "
        "centroids are measured in raw ADC channels");
    seedNote.Write();
    outputFile.Close();

    std::cout << "Accepted " << acceptedFits << " / " << attemptedFits
              << " peak fits and wrote " << outputFileName << ".\n";
    return 0;
}
