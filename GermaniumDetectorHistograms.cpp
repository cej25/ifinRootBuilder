#include "GermaniumDetectorHistograms.h"

#include "AnalysisConfig.h"

#include <TDirectory.h>
#include <TH1D.h>

#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>

GermaniumDetectorHistograms::GermaniumDetectorHistograms() = default;

GermaniumDetectorHistograms::~GermaniumDetectorHistograms() = default;

void GermaniumDetectorHistograms::setCalibratedEnergyAxes()
{
    calibrated_ = true;
    for (auto& item : spectra_) {
        item.second->GetXaxis()->SetTitle("Energy [keV]");
    }
}

TH1D& GermaniumDetectorHistograms::spectrumFor(UShort_t detectorLUT)
{
    const auto existing = spectra_.find(detectorLUT);
    if (existing != spectra_.end()) {
        return *existing->second;
    }

    std::ostringstream suffix;
    suffix << std::setw(2) << std::setfill('0') << detectorLUT;
    const std::string name = "h1_Ge_noVeto_E_LUT" + suffix.str();
    const std::string title = "Germanium LUT " +
        std::to_string(detectorLUT) + ";" +
        (calibrated_ ? "Energy [keV]" : "Energy [raw units]") +
        ";Counts";

    auto histogram = std::make_unique<TH1D>(
        name.c_str(), title.c_str(),
        config::kGermaniumSpectrumBins,
        config::kGermaniumSpectrumMin,
        config::kGermaniumSpectrumMax);
    histogram->SetDirectory(nullptr);
    histogram->SetOption("HIST");

    TH1D& result = *histogram;
    spectra_.emplace(detectorLUT, std::move(histogram));
    return result;
}

void GermaniumDetectorHistograms::fill(UShort_t detectorLUT, double energy)
{
    spectrumFor(detectorLUT).Fill(energy);
}

void GermaniumDetectorHistograms::write(
    TDirectory& germaniumEnergyDirectory) const
{
    TDirectory* directory =
        germaniumEnergyDirectory.mkdir("Individual");
    if (directory == nullptr) {
        throw std::runtime_error(
            "Could not create germanium/individual_detectors directory");
    }

    directory->cd();
    for (const auto& item : spectra_) {
        item.second->Write();
    }
    germaniumEnergyDirectory.cd();
}
