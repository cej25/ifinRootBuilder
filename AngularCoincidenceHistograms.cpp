#include "AngularCoincidenceHistograms.h"

#include "AnalysisConfig.h"

#include <TDirectory.h>
#include <TH1D.h>
#include <TH2I.h>

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

std::unique_ptr<TH2I> makeMatrix(
    const char* name, const char* title,
    const char* allAxisTitle, const char* angularAxisTitle)
{
    const std::string fullTitle = std::string(title) + ";" +
        allAxisTitle + ";" + angularAxisTitle;
    auto histogram = std::make_unique<TH2I>(
        name, fullTitle.c_str(), config::kGammaBins,
        config::kGammaMin, config::kGammaMax, config::kGammaBins,
        config::kGammaMin, config::kGammaMax);
    histogram->SetDirectory(nullptr);
    histogram->SetOption("COLZ");
    return histogram;
}

bool insideHalfOpen(double energy, double minimum, double maximumExclusive)
{
    return energy >= minimum && energy < maximumExclusive;
}

double gateWeight(const CoincidenceGateDefinition& gate, double energy)
{
    if (insideHalfOpen(energy, gate.promptMinimum,
                       gate.promptMaximumExclusive)) {
        return 1.0;
    }
    if (insideHalfOpen(energy, gate.lowerMinimum,
                       gate.lowerMaximumExclusive) ||
        insideHalfOpen(energy, gate.upperMinimum,
                       gate.upperMaximumExclusive)) {
        return -gate.sidebandScale();
    }
    return 0.0;
}

} // namespace

AngularCoincidenceHistograms::AngularCoincidenceHistograms()
{
    allVsForward_ = makeMatrix(
        "h2_Ge_Si_gg_AllvFW",
        "All-detector versus forward-detector gamma-gamma matrix with BGO veto and silicon condition",
        "E_{all} [raw units]", "E_{FW} [raw units]");
    allVsBackward_ = makeMatrix(
        "h2_Ge_Si_gg_AllvBW",
        "All-detector versus backward-detector gamma-gamma matrix with BGO veto and silicon condition",
        "E_{all} [raw units]", "E_{BW} [raw units]");
    allProjectionForward_ = makeGatedSpectrum(
        "h1_Ge_Si_gg_AllvFW_proj",
        "ALL-axis projection of all-versus-forward matrix;E_{all} [raw units];Counts");
    allProjectionBackward_ = makeGatedSpectrum(
        "h1_Ge_Si_gg_AllvBW_proj",
        "ALL-axis projection of all-versus-backward matrix;E_{all} [raw units];Counts");
}

AngularCoincidenceHistograms::~AngularCoincidenceHistograms() = default;

void AngularCoincidenceHistograms::configureGateSet(
    const std::vector<CoincidenceGateDefinition>& definitions,
    const char* suffix, const char* angleTitle, bool calibrated,
    std::vector<GateHistograms>& destination)
{
    destination.clear();
    destination.reserve(definitions.size());
    for (const CoincidenceGateDefinition& gate : definitions) {
        const std::string name = "h1_Ge_Si_gg_" + gate.name + "_" +
            suffix + "_proj";
        const std::string title = std::string(angleTitle) +
            " spectrum gated on " + gate.name +
            " in ALL detectors with sideband subtraction;Energy [raw units];Counts";
        auto spectrum = makeGatedSpectrum(name, title);
        spectrum->Sumw2();
        if (calibrated) spectrum->GetXaxis()->SetTitle("Energy [keV]");
        destination.push_back({gate, std::move(spectrum)});
    }
}

void AngularCoincidenceHistograms::configureGates(
    const std::vector<CoincidenceGateDefinition>& forwardGates,
    const std::vector<CoincidenceGateDefinition>& backwardGates)
{
    configureGateSet(forwardGates, "AllvFW", "Forward", calibrated_,
                     forwardGates_);
    configureGateSet(backwardGates, "AllvBW", "Backward", calibrated_,
                     backwardGates_);
}

void AngularCoincidenceHistograms::fillOrientation(
    double allEnergy, double angularEnergy, unsigned short angularID)
{
    const bool isForward = inIdRange(
        angularID, config::kForwardIDMin, config::kForwardIDMax);
    const bool isBackward = inIdRange(
        angularID, config::kBackwardIDMin, config::kBackwardIDMax);
    if (isForward) {
        allVsForward_->Fill(allEnergy, angularEnergy);
        allProjectionForward_->Fill(allEnergy);
        for (GateHistograms& gate : forwardGates_) {
            const double weight = gateWeight(gate.definition, allEnergy);
            if (weight != 0.0) gate.spectrum->Fill(angularEnergy, weight);
        }
    }
    if (isBackward) {
        allVsBackward_->Fill(allEnergy, angularEnergy);
        allProjectionBackward_->Fill(allEnergy);
        for (GateHistograms& gate : backwardGates_) {
            const double weight = gateWeight(gate.definition, allEnergy);
            if (weight != 0.0) gate.spectrum->Fill(angularEnergy, weight);
        }
    }
}

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
            fillOrientation(gammaEnergies[first], gammaEnergies[second],
                            gammaIDs[second]);
            fillOrientation(gammaEnergies[second], gammaEnergies[first],
                            gammaIDs[first]);
        }
    }
}

void AngularCoincidenceHistograms::setCalibratedEnergyAxes()
{
    calibrated_ = true;
    allVsForward_->GetXaxis()->SetTitle("E_{all} [keV]");
    allVsForward_->GetYaxis()->SetTitle("E_{FW} [keV]");
    allVsBackward_->GetXaxis()->SetTitle("E_{all} [keV]");
    allVsBackward_->GetYaxis()->SetTitle("E_{BW} [keV]");
    allProjectionForward_->GetXaxis()->SetTitle("E_{all} [keV]");
    allProjectionBackward_->GetXaxis()->SetTitle("E_{all} [keV]");
    for (GateHistograms& gate : forwardGates_) {
        gate.spectrum->GetXaxis()->SetTitle("Energy [keV]");
    }
    for (GateHistograms& gate : backwardGates_) {
        gate.spectrum->GetXaxis()->SetTitle("Energy [keV]");
    }
}

void AngularCoincidenceHistograms::merge(
    const AngularCoincidenceHistograms& other)
{
    allVsForward_->Add(other.allVsForward_.get());
    allVsBackward_->Add(other.allVsBackward_.get());
    allProjectionForward_->Add(other.allProjectionForward_.get());
    allProjectionBackward_->Add(other.allProjectionBackward_.get());
    const auto mergeGates = [](std::vector<GateHistograms>& destination,
                               const std::vector<GateHistograms>& source) {
        if (destination.size() != source.size()) {
            throw std::logic_error("Cannot merge different angular gate sets");
        }
        for (std::size_t index = 0; index < destination.size(); ++index) {
            if (!(destination[index].definition == source[index].definition)) {
                throw std::logic_error("Cannot merge different angular gates");
            }
            destination[index].spectrum->Add(source[index].spectrum.get());
        }
    };
    mergeGates(forwardGates_, other.forwardGates_);
    mergeGates(backwardGates_, other.backwardGates_);
}

void AngularCoincidenceHistograms::write(
    TDirectory& gammaCoincidenceDirectory) const
{
    gammaCoincidenceDirectory.cd();
    allVsForward_->Write();
    allVsBackward_->Write();
    allProjectionForward_->Write();
    allProjectionBackward_->Write();

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
    for (const GateHistograms& gate : forwardGates_) {
        forwardDirectory->cd();
        gate.spectrum->Write();
    }
    for (const GateHistograms& gate : backwardGates_) {
        backwardDirectory->cd();
        gate.spectrum->Write();
    }
    gammaCoincidenceDirectory.cd();
}
