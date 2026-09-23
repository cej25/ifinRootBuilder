#include "GammaCoincidenceHistograms.h"

#include "AnalysisConfig.h"

#include <TDirectory.h>
#include <TH1D.h>
#include <TH2I.h>

#include <stdexcept>
#include <string>

namespace {

std::unique_ptr<TH1D> makeProjection(const char* name, const char* title)
{
    auto histogram = std::make_unique<TH1D>(
        name, title,
        config::kGammaBins, config::kGammaMin, config::kGammaMax);
    histogram->SetDirectory(nullptr);
    histogram->SetOption("HIST");
    return histogram;
}

std::unique_ptr<TH2I> makeMatrix(const char* name, const char* title)
{
    auto histogram = std::make_unique<TH2I>(
        name, title, config::kGammaBins, config::kGammaMin,
        config::kGammaMax, config::kGammaBins, config::kGammaMin,
        config::kGammaMax);
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

GammaCoincidenceHistograms::GammaCoincidenceHistograms()
{
    gammaGamma_ = std::make_unique<TH2I>(
        "h2_Ge_gg",
        "Germanium gamma-gamma matrix;E_{#gamma 1} [raw units];E_{#gamma 2} [raw units]",
        config::kGammaBins, config::kGammaMin, config::kGammaMax,
        config::kGammaBins, config::kGammaMin, config::kGammaMax);

    gammaGammaSiliconCoincident_ = std::make_unique<TH2I>(
        "h2_Ge_Si_gg",
        "Germanium gamma-gamma matrix with silicon condition;E_{#gamma 1} [raw units];E_{#gamma 2} [raw units]",
        config::kGammaBins, config::kGammaMin, config::kGammaMax,
        config::kGammaBins, config::kGammaMin, config::kGammaMax);

    gammaGamma_->SetDirectory(nullptr);
    gammaGammaSiliconCoincident_->SetDirectory(nullptr);
    gammaGamma_->SetOption("COLZ");
    gammaGammaSiliconCoincident_->SetOption("COLZ");

}

GammaCoincidenceHistograms::~GammaCoincidenceHistograms() = default;

void GammaCoincidenceHistograms::configureGates(
    const std::vector<CoincidenceGateDefinition>& gates)
{
    gates_.clear();
    gates_.reserve(gates.size());
    for (const CoincidenceGateDefinition& gate : gates) {
        const std::string name =
            "h1_Ge_Si_gg_" + gate.name + "_proj";
        const std::string title =
            "Background-subtracted symmetric gamma-gamma gate " +
            gate.name +
            " with BGO veto and silicon condition;Energy [raw units];Counts";
        auto spectrum = makeProjection(name.c_str(), title.c_str());
        spectrum->Sumw2();
        if (calibrated_) spectrum->GetXaxis()->SetTitle("Energy [keV]");
        gates_.push_back({gate, std::move(spectrum)});
    }
}

void GammaCoincidenceHistograms::configureDoubleGates(
    const std::vector<DoubleCoincidenceGateDefinition>& gates)
{
    doubleGates_.clear();
    doubleGates_.reserve(gates.size());
    for (const DoubleCoincidenceGateDefinition& gate : gates) {
        const std::string matrixName =
            "h2_Ge_Si_gg_DG_" + gate.name;
        const std::string matrixTitle =
            "Symmetrised gamma-gamma matrix after required gamma gate " +
            gate.name +
            " with BGO veto and silicon condition;E_{#gamma 1} [raw units];E_{#gamma 2} [raw units]";
        const std::string spectrumName =
            "h1_Ge_Si_gg_DG_" + gate.name + "_proj";
        const std::string spectrumTitle =
            "Double-gated spectrum " + gate.name +
            " with sideband subtraction;Energy [raw units];Counts";
        auto matrix = makeMatrix(matrixName.c_str(), matrixTitle.c_str());
        auto spectrum = makeProjection(
            spectrumName.c_str(), spectrumTitle.c_str());
        spectrum->Sumw2();
        if (calibrated_) {
            matrix->GetXaxis()->SetTitle("E_{#gamma 1} [keV]");
            matrix->GetYaxis()->SetTitle("E_{#gamma 2} [keV]");
            spectrum->GetXaxis()->SetTitle("Energy [keV]");
        }
        doubleGates_.push_back(
            {gate, std::move(matrix), std::move(spectrum)});
    }
}

void GammaCoincidenceHistograms::fillDoubleGates(
    const std::vector<double>& gammaEnergies)
{
    if (gammaEnergies.size() < 3) return;
    for (DoubleGateHistograms& gate : doubleGates_) {
        std::size_t required = gammaEnergies.size();
        for (std::size_t index = 0; index < gammaEnergies.size(); ++index) {
            if (insideHalfOpen(
                    gammaEnergies[index], gate.definition.requiredMinimum,
                    gate.definition.requiredMaximumExclusive)) {
                required = index;
                break;
            }
        }
        if (required == gammaEnergies.size()) continue;

        for (std::size_t first = 0; first < gammaEnergies.size(); ++first) {
            if (first == required) continue;
            for (std::size_t second = first + 1;
                 second < gammaEnergies.size(); ++second) {
                if (second == required) continue;
                const double firstEnergy = gammaEnergies[first];
                const double secondEnergy = gammaEnergies[second];
                gate.matrix->Fill(firstEnergy, secondEnergy);
                gate.matrix->Fill(secondEnergy, firstEnergy);
                const double firstWeight = gateWeight(
                    gate.definition.secondGate, firstEnergy);
                if (firstWeight != 0.0) {
                    gate.spectrum->Fill(secondEnergy, firstWeight);
                }
                const double secondWeight = gateWeight(
                    gate.definition.secondGate, secondEnergy);
                if (secondWeight != 0.0) {
                    gate.spectrum->Fill(firstEnergy, secondWeight);
                }
            }
        }
    }
}

void GammaCoincidenceHistograms::fillProjection(
    double gateEnergy, double projectedEnergy, bool siliconCoincident)
{
    if (!siliconCoincident) return;
    for (GateHistogram& gate : gates_) {
        const double weight = gateWeight(gate.definition, gateEnergy);
        if (weight != 0.0) gate.spectrum->Fill(projectedEnergy, weight);
    }
}

void GammaCoincidenceHistograms::fillEvent(
    const std::vector<double>& gammaEnergies,
    bool siliconCoincident)
{
    // Symmetrise the matrix. Each distinct hit pair contributes in both axis
    // orders. The projections below use exactly the same ordered-pair logic.
    for (std::size_t first = 0; first < gammaEnergies.size(); ++first) {
        for (std::size_t second = first + 1;
             second < gammaEnergies.size(); ++second) {
            const double firstEnergy = gammaEnergies[first];
            const double secondEnergy = gammaEnergies[second];

            gammaGamma_->Fill(firstEnergy, secondEnergy);
            gammaGamma_->Fill(secondEnergy, firstEnergy);
            if (siliconCoincident) {
                gammaGammaSiliconCoincident_->Fill(firstEnergy, secondEnergy);
                gammaGammaSiliconCoincident_->Fill(secondEnergy, firstEnergy);
            }

            fillProjection(firstEnergy, secondEnergy, siliconCoincident);
            fillProjection(secondEnergy, firstEnergy, siliconCoincident);
        }
    }
    if (siliconCoincident) fillDoubleGates(gammaEnergies);
}

void GammaCoincidenceHistograms::setCalibratedEnergyAxes()
{
    calibrated_ = true;
    gammaGamma_->GetXaxis()->SetTitle("E_{#gamma 1} [keV]");
    gammaGamma_->GetYaxis()->SetTitle("E_{#gamma 2} [keV]");
    gammaGammaSiliconCoincident_->GetXaxis()->SetTitle("E_{#gamma 1} [keV]");
    gammaGammaSiliconCoincident_->GetYaxis()->SetTitle("E_{#gamma 2} [keV]");

    for (GateHistogram& gate : gates_) {
        gate.spectrum->GetXaxis()->SetTitle("Energy [keV]");
    }
    for (DoubleGateHistograms& gate : doubleGates_) {
        gate.matrix->GetXaxis()->SetTitle("E_{#gamma 1} [keV]");
        gate.matrix->GetYaxis()->SetTitle("E_{#gamma 2} [keV]");
        gate.spectrum->GetXaxis()->SetTitle("Energy [keV]");
    }
}

void GammaCoincidenceHistograms::merge(
    const GammaCoincidenceHistograms& other)
{
    gammaGamma_->Add(other.gammaGamma_.get());
    gammaGammaSiliconCoincident_->Add(
        other.gammaGammaSiliconCoincident_.get());
    if (gates_.size() != other.gates_.size()) {
        throw std::logic_error("Cannot merge different symmetric gate sets");
    }
    for (std::size_t index = 0; index < gates_.size(); ++index) {
        if (!(gates_[index].definition == other.gates_[index].definition)) {
            throw std::logic_error("Cannot merge different symmetric gates");
        }
        gates_[index].spectrum->Add(other.gates_[index].spectrum.get());
    }
    if (doubleGates_.size() != other.doubleGates_.size()) {
        throw std::logic_error("Cannot merge different double-gate sets");
    }
    for (std::size_t index = 0; index < doubleGates_.size(); ++index) {
        if (!(doubleGates_[index].definition ==
              other.doubleGates_[index].definition)) {
            throw std::logic_error("Cannot merge different double gates");
        }
        doubleGates_[index].matrix->Add(
            other.doubleGates_[index].matrix.get());
        doubleGates_[index].spectrum->Add(
            other.doubleGates_[index].spectrum.get());
    }
}

void GammaCoincidenceHistograms::write(TDirectory& parentDirectory) const
{
    TDirectory* directory = parentDirectory.mkdir("Coincidences");
    if (directory == nullptr) {
        throw std::runtime_error(
            "Could not create Coincidences output directory");
    }

    directory->cd();
    gammaGamma_->Write();
    gammaGammaSiliconCoincident_->Write();

    TDirectory* gatedDirectory = directory->mkdir("Gated");
    if (gatedDirectory == nullptr) {
        throw std::runtime_error(
            "Could not create Coincidences/Gated directory");
    }
    gatedDirectory->cd();
    for (const GateHistogram& gate : gates_) gate.spectrum->Write();

    TDirectory* doubleDirectory = directory->mkdir("DoubleGated");
    if (doubleDirectory == nullptr) {
        throw std::runtime_error(
            "Could not create Coincidences/DoubleGated directory");
    }
    doubleDirectory->cd();
    for (const DoubleGateHistograms& gate : doubleGates_) {
        gate.matrix->Write();
        gate.spectrum->Write();
    }
    parentDirectory.cd();
}
