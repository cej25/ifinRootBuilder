#include "GermaniumDetectorHistograms.h"

#include "AnalysisConfig.h"

#include <TDirectory.h>
#include <TH1D.h>
#include <TH2I.h>

#include <algorithm>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>

GermaniumDetectorHistograms::GermaniumDetectorHistograms() = default;

GermaniumDetectorHistograms::~GermaniumDetectorHistograms() = default;

void GermaniumDetectorHistograms::setCalibratedEnergyAxes()
{
    calibrated_ = true;
    for (auto& item : calibratedSpectra_) {
        item.second->GetXaxis()->SetTitle("Energy [keV]");
    }
    for (auto& item : timeSpectra_) {
        item.second->GetYaxis()->SetTitle("Energy [keV]");
    }
}

void GermaniumDetectorHistograms::setRunningTimeRange(
    double totalRunningTimeSeconds)
{
    if (!timeSpectra_.empty()) {
        throw std::logic_error(
            "Cannot change the running-time range after filling Ge matrices");
    }
    runningTimeMaxSeconds_ = std::max(1.0, totalRunningTimeSeconds + 1.0);
}

TH1D& GermaniumDetectorHistograms::calibratedSpectrumFor(
    UShort_t detectorID)
{
    const auto existing = calibratedSpectra_.find(detectorID);
    if (existing != calibratedSpectra_.end()) {
        return *existing->second;
    }

    std::ostringstream suffix;
    suffix << std::setw(2) << std::setfill('0') << detectorID;
    const std::string name = "h1_Ge_noVeto_E_ID" + suffix.str();
    const std::string title = "Germanium ID " +
        std::to_string(detectorID) + ";" +
        (calibrated_ ? "Energy [keV]" : "Energy [raw units]") +
        ";Counts";

    auto histogram = std::make_unique<TH1D>(
        name.c_str(), title.c_str(),
        config::kIndividualGermaniumCalibratedBins,
        config::kIndividualGermaniumCalibratedMin,
        config::kIndividualGermaniumCalibratedMax);
    histogram->SetDirectory(nullptr);
    histogram->SetOption("HIST");

    TH1D& result = *histogram;
    calibratedSpectra_.emplace(detectorID, std::move(histogram));
    return result;
}

TH1D& GermaniumDetectorHistograms::rawSpectrumFor(UShort_t detectorID)
{
    const auto existing = rawSpectra_.find(detectorID);
    if (existing != rawSpectra_.end()) return *existing->second;

    std::ostringstream suffix;
    suffix << std::setw(2) << std::setfill('0') << detectorID;
    const std::string name = "h1_Ge_noVeto_E_raw_ID" + suffix.str();
    const std::string title = "Germanium ID " +
        std::to_string(detectorID) +
        " raw energy;Energy [ADC channel];Counts";
    auto histogram = std::make_unique<TH1D>(
        name.c_str(), title.c_str(), config::kIndividualGermaniumRawBins,
        config::kIndividualGermaniumRawMin,
        config::kIndividualGermaniumRawMax);
    histogram->SetDirectory(nullptr);
    histogram->SetOption("HIST");
    TH1D& result = *histogram;
    rawSpectra_.emplace(detectorID, std::move(histogram));
    return result;
}

TH2I& GermaniumDetectorHistograms::timeSpectrumFor(UShort_t detectorID)
{
    const auto existing = timeSpectra_.find(detectorID);
    if (existing != timeSpectra_.end()) return *existing->second;

    std::ostringstream suffix;
    suffix << std::setw(2) << std::setfill('0') << detectorID;
    const std::string name = "h2_Ge_noVeto_EvTime_ID" + suffix.str();
    const std::string title = "Germanium ID " +
        std::to_string(detectorID) +
        " energy versus total running time;Total running time [s];" +
        (calibrated_ ? "Energy [keV]" : "Energy [raw units]");
    auto histogram = std::make_unique<TH2I>(
        name.c_str(), title.c_str(), config::kDriftTimeBins, 0.0,
        runningTimeMaxSeconds_, config::kDriftEnergyBins,
        config::kDriftEnergyMin, config::kDriftEnergyMax);
    histogram->SetDirectory(nullptr);
    histogram->SetOption("COLZ");
    TH2I& result = *histogram;
    timeSpectra_.emplace(detectorID, std::move(histogram));
    return result;
}

void GermaniumDetectorHistograms::fill(
    UShort_t detectorID, double calibratedEnergy,
    double rawEnergy, double runningTimeSeconds, bool runningTimeValid)
{
    calibratedSpectrumFor(detectorID).Fill(calibratedEnergy);
    rawSpectrumFor(detectorID).Fill(rawEnergy);
    if (runningTimeValid) {
        timeSpectrumFor(detectorID).Fill(runningTimeSeconds, calibratedEnergy);
    }
}

void GermaniumDetectorHistograms::merge(
    const GermaniumDetectorHistograms& other)
{
    for (const auto& item : other.calibratedSpectra_) {
        calibratedSpectrumFor(item.first).Add(item.second.get());
    }
    for (const auto& item : other.rawSpectra_) {
        rawSpectrumFor(item.first).Add(item.second.get());
    }
    for (const auto& item : other.timeSpectra_) {
        timeSpectrumFor(item.first).Add(item.second.get());
    }
}

void GermaniumDetectorHistograms::write(
    TDirectory& germaniumEnergyDirectory,
    TDirectory& germaniumTimeDirectory) const
{
    TDirectory* individualEnergy =
        germaniumEnergyDirectory.mkdir("Individual");
    TDirectory* calibratedDirectory = individualEnergy != nullptr
        ? individualEnergy->mkdir("Calibrated") : nullptr;
    TDirectory* rawDirectory = individualEnergy != nullptr
        ? individualEnergy->mkdir("Raw") : nullptr;
    TDirectory* individualTime = germaniumTimeDirectory.mkdir("Individual");
    if (individualEnergy == nullptr || calibratedDirectory == nullptr ||
        rawDirectory == nullptr || individualTime == nullptr) {
        throw std::runtime_error(
            "Could not create germanium individual Energy/Time directories");
    }

    calibratedDirectory->cd();
    for (const auto& item : calibratedSpectra_) {
        item.second->Write();
    }
    rawDirectory->cd();
    for (const auto& item : rawSpectra_) {
        item.second->Write();
    }
    individualTime->cd();
    for (const auto& item : timeSpectra_) {
        item.second->Write();
    }
    germaniumEnergyDirectory.cd();
}
