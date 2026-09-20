#ifndef GERMANIUM_CONDITION_HISTOGRAMS_H
#define GERMANIUM_CONDITION_HISTOGRAMS_H

#include <memory>
#include <vector>

class TDirectory;
class TH1D;
class TH2I;

class GermaniumConditionHistograms {
public:
    GermaniumConditionHistograms();
    ~GermaniumConditionHistograms();

    GermaniumConditionHistograms(const GermaniumConditionHistograms&) = delete;
    GermaniumConditionHistograms& operator=(
        const GermaniumConditionHistograms&) = delete;

    void setCalibratedEnergyAxes();
    void fillHit(double energy, double absoluteTimeSeconds,
                 bool survivesBgoVeto, bool siliconCoincident,
                 bool foldValid);
    void fillUnconditionedCoincidences(
        const std::vector<double>& germaniumEnergies);
    void fillBgoVetoedSiliconCoincidences(
        const std::vector<double>& germaniumEnergies);
    void fillEvent(unsigned int germaniumMultiplicityAfterBgoVeto,
                   bool siliconCoincident, bool foldValid);
    void merge(const GermaniumConditionHistograms& other);
    void write(TDirectory& germaniumDirectory,
               TDirectory& energyDirectory,
               TDirectory& timeDirectory) const;

private:
    std::unique_ptr<TH1D> unconditionedGammaGammaProjection_;
    std::unique_ptr<TH1D> bgoVetoedSiliconGammaGammaProjection_;
    std::unique_ptr<TH1D> bgoVetoed_;
    std::unique_ptr<TH1D> siliconBgoVetoed_;
    std::unique_ptr<TH1D> multiplicityAfterBgoVeto_;
    std::unique_ptr<TH1D> foldValidMultiplicity_;
    std::unique_ptr<TH1D> foldValidEnergy_;
    std::unique_ptr<TH1D> foldInclusiveEnergy_;
    std::unique_ptr<TH2I> energyVsTime_;
};

#endif
