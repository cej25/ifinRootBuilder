#include "GammaCoincidenceHistograms.h"

#include "AnalysisConfig.h"

#include <TDirectory.h>
#include <TH1D.h>
#include <TH2I.h>

#include <stdexcept>

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

bool insideHalfOpen(double energy, double minimum, double maximumExclusive)
{
    return energy >= minimum && energy < maximumExclusive;
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

    backgroundSubtractedProjectionSiliconCoincident_ = makeProjection(
        "h1_Ge_Si_gg_Gate292to299_proj",
        "Background-subtracted 292-299 gated spectrum with BGO veto and silicon condition;Energy [raw units];Counts");

    gammaGamma_->SetDirectory(nullptr);
    gammaGammaSiliconCoincident_->SetDirectory(nullptr);
    gammaGamma_->SetOption("COLZ");
    gammaGammaSiliconCoincident_->SetOption("COLZ");

    // Weighted positive/negative fills require explicit sum-of-squares arrays.
    backgroundSubtractedProjectionSiliconCoincident_->Sumw2();
}

GammaCoincidenceHistograms::~GammaCoincidenceHistograms() = default;

GammaCoincidenceHistograms::ProjectionWindow
GammaCoincidenceHistograms::projectionWindow(double energy) const
{
    if (insideHalfOpen(energy, config::kGammaGateMin,
                       config::kGammaGateMaxExclusive)) {
        return ProjectionWindow::Prompt;
    }
    if (insideHalfOpen(energy, config::kGammaLowerSidebandMin,
                       config::kGammaLowerSidebandMaxExclusive)) {
        return ProjectionWindow::LowerSideband;
    }
    if (insideHalfOpen(energy, config::kGammaUpperSidebandMin,
                       config::kGammaUpperSidebandMaxExclusive)) {
        return ProjectionWindow::UpperSideband;
    }
    return ProjectionWindow::None;
}

void GammaCoincidenceHistograms::fillProjection(
    double gateEnergy, double projectedEnergy, bool siliconCoincident)
{
    const ProjectionWindow window = projectionWindow(gateEnergy);
    if (window == ProjectionWindow::None) {
        return;
    }

    double subtractionWeight = 1.0;
    if (window != ProjectionWindow::Prompt) {
        subtractionWeight = -config::kGammaSidebandScale;
    }

    if (siliconCoincident) {
        backgroundSubtractedProjectionSiliconCoincident_->Fill(
            projectedEnergy, subtractionWeight);
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
}

void GammaCoincidenceHistograms::setCalibratedEnergyAxes()
{
    gammaGamma_->GetXaxis()->SetTitle("E_{#gamma 1} [keV]");
    gammaGamma_->GetYaxis()->SetTitle("E_{#gamma 2} [keV]");
    gammaGammaSiliconCoincident_->GetXaxis()->SetTitle("E_{#gamma 1} [keV]");
    gammaGammaSiliconCoincident_->GetYaxis()->SetTitle("E_{#gamma 2} [keV]");

    TH1D* projections[] = {
        backgroundSubtractedProjectionSiliconCoincident_.get()};
    for (TH1D* projection : projections) {
        projection->GetXaxis()->SetTitle("Energy [keV]");
    }
}

void GammaCoincidenceHistograms::merge(
    const GammaCoincidenceHistograms& other)
{
    gammaGamma_->Add(other.gammaGamma_.get());
    gammaGammaSiliconCoincident_->Add(
        other.gammaGammaSiliconCoincident_.get());
    backgroundSubtractedProjectionSiliconCoincident_->Add(
        other.backgroundSubtractedProjectionSiliconCoincident_.get());
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
    backgroundSubtractedProjectionSiliconCoincident_->Write();
    parentDirectory.cd();
}
