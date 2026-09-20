#include "AngularCoincidenceHistograms.h"

#include "AnalysisConfig.h"

#include <TDirectory.h>
#include <TH1D.h>

#include <stdexcept>
#include <string>

namespace {

bool inIdRange(unsigned short id, unsigned short minimum,
                unsigned short maximum)
{
    return id >= minimum && id <= maximum;
}

std::unique_ptr<TH1D> makeGatedSpectrum(
    const std::string& name, const std::string& title)
{
    auto histogram = std::make_unique<TH1D>(
        name.c_str(), title.c_str(),
        config::kGammaBins, config::kGammaMin, config::kGammaMax);
    histogram->SetDirectory(nullptr);
    histogram->SetOption("HIST");
    return histogram;
}

} // namespace

AngularCoincidenceHistograms::AngularCoincidenceHistograms()
{
    gates_.reserve(config::kRddsGates.size());
    for (const auto& gate : config::kRddsGates) {
        const std::string baseName =
            "h1_Ge_Si_gg_" + std::string(gate.name);
        gates_.push_back({
            gate.name,
            gate.minimum,
            gate.maximumExclusive,
            makeGatedSpectrum(
                baseName + "_AllvFW_proj",
                std::string("Forward spectrum gated on ") + gate.name +
                    " in all detectors;Energy [raw units];Counts"),
            makeGatedSpectrum(
                baseName + "_AllvBW_proj",
                std::string("Backward spectrum gated on ") + gate.name +
                    " in all detectors;Energy [raw units];Counts")
        });
    }
}

AngularCoincidenceHistograms::~AngularCoincidenceHistograms() = default;

void AngularCoincidenceHistograms::fillEvent(
    const std::vector<double>& gammaEnergies,
    const std::vector<unsigned short>& gammaIDs)
{
    if (gammaEnergies.size() != gammaIDs.size()) {
        throw std::logic_error(
            "Angular coincidence energies and IDs have different sizes");
    }

    for (std::size_t first = 0; first < gammaEnergies.size(); ++first) {
        for (std::size_t second = first + 1;
             second < gammaEnergies.size(); ++second) {
            const double energies[2] = {
                gammaEnergies[first], gammaEnergies[second]};
            const unsigned short ids[2] = {
                gammaIDs[first], gammaIDs[second]};
            for (int orientation = 0; orientation < 2; ++orientation) {
                const double allEnergy = energies[orientation];
                const double angularEnergy = energies[1 - orientation];
                const unsigned short angularID = ids[1 - orientation];
                const bool isForward = inIdRange(
                    angularID, config::kForwardIDMin,
                    config::kForwardIDMax);
                const bool isBackward = inIdRange(
                    angularID, config::kBackwardIDMin,
                    config::kBackwardIDMax);
                for (GateHistograms& gate : gates_) {
                    if (allEnergy < gate.minimum ||
                        allEnergy >= gate.maximumExclusive) {
                        continue;
                    }
                    if (isForward) {
                        gate.forward->Fill(angularEnergy);
                    }
                    if (isBackward) {
                        gate.backward->Fill(angularEnergy);
                    }
                }
            }
        }
    }
}

void AngularCoincidenceHistograms::setCalibratedEnergyAxes()
{
    for (GateHistograms& gate : gates_) {
        gate.forward->GetXaxis()->SetTitle("Energy [keV]");
        gate.backward->GetXaxis()->SetTitle("Energy [keV]");
    }
}

void AngularCoincidenceHistograms::merge(
    const AngularCoincidenceHistograms& other)
{
    if (gates_.size() != other.gates_.size()) {
        throw std::logic_error("Cannot merge different angular gate sets");
    }
    for (std::size_t index = 0; index < gates_.size(); ++index) {
        gates_[index].forward->Add(other.gates_[index].forward.get());
        gates_[index].backward->Add(other.gates_[index].backward.get());
    }
}

void AngularCoincidenceHistograms::write(
    TDirectory& gammaCoincidenceDirectory) const
{
    TDirectory* gatedDirectory =
        gammaCoincidenceDirectory.GetDirectory("Gated");
    if (gatedDirectory == nullptr) {
        throw std::runtime_error("Coincidences/Gated directory is missing");
    }
    TDirectory* forwardDirectory = gatedDirectory->mkdir("AllvFW");
    TDirectory* backwardDirectory = gatedDirectory->mkdir("AllvBW");
    if (forwardDirectory == nullptr || backwardDirectory == nullptr) {
        throw std::runtime_error("Could not create angular gated directories");
    }
    for (const GateHistograms& gate : gates_) {
        forwardDirectory->cd();
        gate.forward->Write();
        backwardDirectory->cd();
        gate.backward->Write();
    }
    gammaCoincidenceDirectory.cd();
}
