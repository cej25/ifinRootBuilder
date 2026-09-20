#include "DetectorHistograms.h"

#include <TDirectory.h>
#include <TH1D.h>

#include <stdexcept>
#include <string>

DetectorHistograms::DetectorHistograms(
    const config::DetectorDefinition& definition)
    : definition_(&definition)
{
    const std::string prefix = "h1_" + std::string(definition.histogramTag) +
        (definition.type == config::kGermaniumType ? "_noVeto" : "");
    const std::string energyName = prefix + "_E";
    const std::string energyTitle = std::string(definition.title)
        + " energy;Energy [raw units];Counts";
    energy_ = std::make_unique<TH1D>(
        energyName.c_str(), energyTitle.c_str(), definition.energyBins,
        definition.energyMin, definition.energyMax);

    const std::string timeName = prefix + "_Time";
    const std::string timeTitle = std::string(definition.title)
        + " relative time;Relative time [ns];Counts";
    time_ = std::make_unique<TH1D>(
        timeName.c_str(), timeTitle.c_str(), config::kTimeBins,
        config::kTimeMinNs, config::kTimeMaxNs);

    const std::string multiplicityName =
        prefix + "_Multiplicity";
    const std::string multiplicityTitle = std::string(definition.title)
        + " hit multiplicity;Hits per event;Events";
    multiplicity_ = std::make_unique<TH1D>(
        multiplicityName.c_str(), multiplicityTitle.c_str(),
        config::kMultiplicityBins, config::kMultiplicityMin,
        config::kMultiplicityMax);

    const std::string lutName = prefix + "_LUT";
    const std::string lutTitle = std::string(definition.title)
        + " LUT occupancy;detectorLUT;Hits";
    lutOccupancy_ = std::make_unique<TH1D>(
        lutName.c_str(), lutTitle.c_str(),
        config::kLutBins, config::kLutMin, config::kLutMax);

    energy_->SetDirectory(nullptr);
    time_->SetDirectory(nullptr);
    multiplicity_->SetDirectory(nullptr);
    lutOccupancy_->SetDirectory(nullptr);
    energy_->SetOption("HIST");
    time_->SetOption("HIST");
    multiplicity_->SetOption("HIST");
    lutOccupancy_->SetOption("HIST");
}

DetectorHistograms::~DetectorHistograms() = default;
DetectorHistograms::DetectorHistograms(DetectorHistograms&&) noexcept = default;
DetectorHistograms& DetectorHistograms::operator=(DetectorHistograms&&) noexcept = default;

unsigned short DetectorHistograms::detectorType() const
{
    return definition_->type;
}

const char* DetectorHistograms::name() const
{
    return definition_->name;
}

const char* DetectorHistograms::directoryName() const
{
    return definition_->directory;
}

void DetectorHistograms::fillHit(unsigned short detectorLUT,
                                 double energy,
                                 double relativeTimeNs)
{
    lutOccupancy_->Fill(detectorLUT);
    energy_->Fill(energy);
    time_->Fill(relativeTimeNs);
}

void DetectorHistograms::fillMultiplicity(unsigned int multiplicity)
{
    multiplicity_->Fill(multiplicity);
}

void DetectorHistograms::setCalibratedEnergyAxis()
{
    energy_->GetXaxis()->SetTitle("Energy [keV]");
}

TDirectory* DetectorHistograms::write(TDirectory& parentDirectory) const
{
    TDirectory* detectorDirectory = parentDirectory.mkdir(directoryName());
    if (detectorDirectory == nullptr) {
        throw std::runtime_error(
            "Could not create output directory for " + std::string(name()));
    }

    detectorDirectory->cd();
    multiplicity_->Write();
    lutOccupancy_->Write();

    TDirectory* energyDirectory = detectorDirectory->mkdir("Energy");
    TDirectory* timeDirectory = detectorDirectory->mkdir("Time");
    if (energyDirectory == nullptr || timeDirectory == nullptr) {
        throw std::runtime_error(
            "Could not create Energy/Time directories for " +
            std::string(directoryName()));
    }
    energyDirectory->cd();
    energy_->Write();
    timeDirectory->cd();
    time_->Write();
    parentDirectory.cd();
    return detectorDirectory;
}
