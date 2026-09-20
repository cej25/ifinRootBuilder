#include "GermaniumConditionHistograms.h"

#include "AnalysisConfig.h"

#include <TCanvas.h>
#include <TColor.h>
#include <TDirectory.h>
#include <TH1D.h>
#include <TH2I.h>
#include <TLegend.h>

#include <algorithm>
#include <stdexcept>

namespace {

std::unique_ptr<TH1D> makeSpectrum(const char* name, const char* title)
{
    auto histogram = std::make_unique<TH1D>(
        name, title, config::kGermaniumSpectrumBins,
        config::kGermaniumSpectrumMin, config::kGermaniumSpectrumMax);
    histogram->SetDirectory(nullptr);
    histogram->SetOption("HIST");
    return histogram;
}

std::unique_ptr<TH1D> makeMultiplicity(const char* name, const char* title)
{
    auto histogram = std::make_unique<TH1D>(
        name, title, config::kMultiplicityBins,
        config::kMultiplicityMin, config::kMultiplicityMax);
    histogram->SetDirectory(nullptr);
    histogram->SetOption("HIST");
    return histogram;
}

std::unique_ptr<TH1D> makeGammaProjection(
    const char* name, const char* title)
{
    auto histogram = std::make_unique<TH1D>(
        name, title, config::kGammaBins,
        config::kGammaMin, config::kGammaMax);
    histogram->SetDirectory(nullptr);
    histogram->SetOption("HIST");
    return histogram;
}

void fillSymmetrisedProjection(
    TH1D& projection, const std::vector<double>& germaniumEnergies)
{
    for (std::size_t first = 0; first < germaniumEnergies.size(); ++first) {
        for (std::size_t second = first + 1;
             second < germaniumEnergies.size(); ++second) {
            const double firstEnergy = germaniumEnergies[first];
            const double secondEnergy = germaniumEnergies[second];
            projection.Fill(firstEnergy);
            projection.Fill(secondEnergy);
        }
    }
}

} // namespace

GermaniumConditionHistograms::GermaniumConditionHistograms()
    : bgoVetoed_(makeSpectrum(
          "h1_Ge_E",
          "Germanium energy after BGO veto;Energy [raw units];Counts")),
      siliconBgoVetoed_(makeSpectrum(
          "h1_Ge_Si_E",
          "Germanium energy after BGO veto with silicon condition;Energy [raw units];Counts")),
      multiplicityAfterBgoVeto_(makeMultiplicity(
          "h1_Ge_Multiplicity",
          "Germanium multiplicity after BGO veto;Germanium hits per event;Events")),
      foldValidMultiplicity_(makeMultiplicity(
          "h1_Ge_Si_FV_Multiplicity",
          "Germanium multiplicity after BGO veto for FoldValid silicon-conditioned events;Germanium hits per event;Events")),
      foldValidEnergy_(makeSpectrum(
          "h1_Ge_Si_FV_E",
          "Germanium energy for FoldValid events after BGO veto and silicon condition;Energy [raw units];Counts")),
      foldInclusiveEnergy_(makeSpectrum(
          "h1_Ge_Si_E",
          "Germanium energy for all FoldValid states after BGO veto and silicon condition;Energy [raw units];Counts"))
{
    unconditionedGammaGammaProjection_ = makeGammaProjection(
        "h1_Ge_noVeto_gg_proj",
        "Projection of symmetrised Ge-Ge matrix without BGO or silicon conditions;Energy [raw units];Counts");

    bgoVetoedSiliconGammaGammaProjection_ = makeGammaProjection(
        "h1_Ge_Si_gg_proj",
        "Projection of symmetrised Ge-Ge matrix after BGO veto with silicon condition;Energy [raw units];Counts");

    energyVsTime_ = std::make_unique<TH2I>(
        "h2_Ge_noVeto_EvTime",
        "Germanium energy versus absolute time;absoluteTime [s];Energy [raw units]",
        config::kDriftTimeBins, config::kDriftTimeMinSeconds,
        config::kDriftTimeMaxSeconds, config::kDriftEnergyBins,
        config::kDriftEnergyMin, config::kDriftEnergyMax);
    energyVsTime_->SetDirectory(nullptr);
    energyVsTime_->SetOption("COLZ");
}

GermaniumConditionHistograms::~GermaniumConditionHistograms() = default;

void GermaniumConditionHistograms::setCalibratedEnergyAxes()
{
    TH1D* spectra[] = {unconditionedGammaGammaProjection_.get(),
                       bgoVetoedSiliconGammaGammaProjection_.get(),
                       bgoVetoed_.get(), siliconBgoVetoed_.get(),
                       foldValidEnergy_.get(), foldInclusiveEnergy_.get()};
    for (TH1D* spectrum : spectra) {
        spectrum->GetXaxis()->SetTitle("Energy [keV]");
    }
    energyVsTime_->GetYaxis()->SetTitle("Energy [keV]");
}

void GermaniumConditionHistograms::fillHit(
    double energy, double absoluteTimeSeconds,
    bool survivesBgoVeto, bool siliconCoincident, bool foldValid)
{
    energyVsTime_->Fill(absoluteTimeSeconds, energy);
    if (!survivesBgoVeto) {
        return;
    }
    bgoVetoed_->Fill(energy);
    if (!siliconCoincident) {
        return;
    }
    siliconBgoVetoed_->Fill(energy);
    foldInclusiveEnergy_->Fill(energy);
    if (foldValid) {
        foldValidEnergy_->Fill(energy);
    }
}

void GermaniumConditionHistograms::fillUnconditionedCoincidences(
    const std::vector<double>& germaniumEnergies)
{
    fillSymmetrisedProjection(
        *unconditionedGammaGammaProjection_, germaniumEnergies);
}

void GermaniumConditionHistograms::fillBgoVetoedSiliconCoincidences(
    const std::vector<double>& germaniumEnergies)
{
    fillSymmetrisedProjection(
        *bgoVetoedSiliconGammaGammaProjection_, germaniumEnergies);
}

void GermaniumConditionHistograms::fillEvent(
    unsigned int germaniumMultiplicityAfterBgoVeto,
    bool siliconCoincident, bool foldValid)
{
    multiplicityAfterBgoVeto_->Fill(germaniumMultiplicityAfterBgoVeto);
    if (siliconCoincident && foldValid) {
        foldValidMultiplicity_->Fill(germaniumMultiplicityAfterBgoVeto);
    }
}

void GermaniumConditionHistograms::merge(
    const GermaniumConditionHistograms& other)
{
    unconditionedGammaGammaProjection_->Add(
        other.unconditionedGammaGammaProjection_.get());
    bgoVetoedSiliconGammaGammaProjection_->Add(
        other.bgoVetoedSiliconGammaGammaProjection_.get());
    bgoVetoed_->Add(other.bgoVetoed_.get());
    siliconBgoVetoed_->Add(other.siliconBgoVetoed_.get());
    multiplicityAfterBgoVeto_->Add(other.multiplicityAfterBgoVeto_.get());
    foldValidMultiplicity_->Add(other.foldValidMultiplicity_.get());
    foldValidEnergy_->Add(other.foldValidEnergy_.get());
    foldInclusiveEnergy_->Add(other.foldInclusiveEnergy_.get());
    energyVsTime_->Add(other.energyVsTime_.get());
}

void GermaniumConditionHistograms::write(
    TDirectory& germaniumDirectory, TDirectory& energyDirectory,
    TDirectory& timeDirectory) const
{
    germaniumDirectory.cd();
    multiplicityAfterBgoVeto_->Write();

    energyDirectory.cd();
    unconditionedGammaGammaProjection_->Write();
    bgoVetoedSiliconGammaGammaProjection_->Write();
    bgoVetoed_->Write();
    siliconBgoVetoed_->Write();

    auto normalizedBaseline = std::make_unique<TH1D>(*bgoVetoed_);
    normalizedBaseline->SetName(
        "h1_Ge_E_NormalizedToSi");
    auto siliconOverlay = std::make_unique<TH1D>(*siliconBgoVetoed_);
    siliconOverlay->SetName(
        "h1_Ge_Si_E_Overlay");
    const double baselineCounts = normalizedBaseline->Integral(
        1, normalizedBaseline->GetNbinsX());
    const double siliconCounts = siliconOverlay->Integral(
        1, siliconOverlay->GetNbinsX());
    if (baselineCounts > 0.0) {
        normalizedBaseline->Scale(siliconCounts / baselineCounts);
    }
    normalizedBaseline->SetLineColor(kAzure + 2);
    normalizedBaseline->SetLineWidth(2);
    normalizedBaseline->SetStats(false);
    normalizedBaseline->SetOption("HIST");
    siliconOverlay->SetLineColor(kRed + 1);
    siliconOverlay->SetLineWidth(2);
    siliconOverlay->SetStats(false);
    siliconOverlay->SetOption("HIST");
    normalizedBaseline->SetMaximum(
        1.10 * std::max(normalizedBaseline->GetMaximum(),
                        siliconOverlay->GetMaximum()));

    TCanvas comparisonCanvas(
        "c_Ge_Si_E_Comparison",
        "BGO-vetoed germanium with and without silicon condition", 1100, 700);
    normalizedBaseline->Draw("HIST");
    siliconOverlay->Draw("HIST SAME");
    TLegend legend(0.56, 0.76, 0.88, 0.88);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.AddEntry(normalizedBaseline.get(),
                    "BGO vetoed (area normalized)", "l");
    legend.AddEntry(siliconOverlay.get(),
                    "BGO vetoed + silicon condition", "l");
    legend.Draw();
    comparisonCanvas.Write();

    TDirectory* foldDirectory = energyDirectory.mkdir("FoldValid");
    if (foldDirectory == nullptr) {
        throw std::runtime_error(
            "Could not create Germanium/Energy/FoldValid directory");
    }
    foldDirectory->cd();
    foldValidMultiplicity_->Write();
    foldValidEnergy_->Write();
    foldInclusiveEnergy_->Write();

    timeDirectory.cd();
    energyVsTime_->Write();
    germaniumDirectory.cd();
}
