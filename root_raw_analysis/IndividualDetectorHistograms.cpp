#include "IndividualDetectorHistograms.h"

#include <TDirectory.h>
#include <TH1D.h>

#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>

IndividualDetectorHistograms::IndividualDetectorHistograms(
    const config::DetectorDefinition& definition)
    : definition_(&definition)
{
}

IndividualDetectorHistograms::~IndividualDetectorHistograms() = default;

TH1D& IndividualDetectorHistograms::spectrumFor(unsigned short detectorID)
{
    const auto existing = spectra_.find(detectorID);
    if (existing != spectra_.end()) {
        return *existing->second;
    }

    std::ostringstream suffix;
    suffix << std::setw(2) << std::setfill('0') << detectorID;
    const std::string name = "h1_" +
        std::string(definition_->histogramTag) + "_E_ID" + suffix.str();
    const std::string title = std::string(definition_->title) + " ID " +
        std::to_string(detectorID) + ";Energy [raw units];Counts";

    auto histogram = std::make_unique<TH1D>(
        name.c_str(), title.c_str(),
        definition_->energyBins,
        definition_->energyMin,
        definition_->energyMax);
    histogram->SetDirectory(nullptr);
    histogram->SetOption("HIST");

    TH1D& result = *histogram;
    spectra_.emplace(detectorID, std::move(histogram));
    return result;
}

void IndividualDetectorHistograms::fill(
    unsigned short detectorID, double energy)
{
    spectrumFor(detectorID).Fill(energy);
}

void IndividualDetectorHistograms::merge(
    const IndividualDetectorHistograms& other)
{
    if (definition_->type != other.definition_->type) {
        throw std::logic_error("Cannot merge different detector types");
    }
    for (const auto& item : other.spectra_) {
        spectrumFor(item.first).Add(item.second.get());
    }
}

void IndividualDetectorHistograms::write(
    TDirectory& energyDirectory) const
{
    TDirectory* directory = energyDirectory.mkdir("Individual");
    if (directory == nullptr) {
        throw std::runtime_error(
            "Could not create individual detector directory for " +
            std::string(definition_->name));
    }

    directory->cd();
    for (const auto& item : spectra_) {
        item.second->Write();
    }
    energyDirectory.cd();
}
