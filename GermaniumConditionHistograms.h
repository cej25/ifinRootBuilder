#ifndef GERMANIUM_CONDITION_HISTOGRAMS_H
#define GERMANIUM_CONDITION_HISTOGRAMS_H

#include <array>
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
    void setRunningTimeRange(double totalRunningTimeSeconds);
    void fillHit(double energy, double runningTimeSeconds,
                 bool runningTimeValid,
                 bool survivesBgoVeto, bool siliconCoincident,
                 bool foldValid);
    void fillUnconditionedCoincidences(
        const std::vector<double>& germaniumEnergies);
    void fillBgoVetoedSiliconCoincidences(
        const std::vector<double>& germaniumEnergies);
    void fillBgoVetoedSiliconEvent(
        const std::vector<double>& germaniumEnergies,
        const std::vector<unsigned short>& germaniumIDs,
        unsigned int germaniumMultiplicityBeforeBgoVeto,
        unsigned int bgoMultiplicity,
        const std::vector<double>& acceptedSiliconEnergies);
    void fillEvent(unsigned int germaniumMultiplicityAfterBgoVeto,
                   bool siliconCoincident, bool foldValid);
    void merge(const GermaniumConditionHistograms& other);
    void write(TDirectory& germaniumDirectory,
               TDirectory& energyDirectory,
               TDirectory& timeDirectory) const;
    void writeParticleGamma(TDirectory& coincidenceDirectory) const;

private:
    std::unique_ptr<TH1D> unconditionedGammaGammaProjection_;
    std::unique_ptr<TH1D> bgoVetoedSiliconGammaGammaProjection_;
    std::unique_ptr<TH1D> bgoVetoed_;
    std::unique_ptr<TH1D> siliconBgoVetoed_;
    std::unique_ptr<TH1D> siliconBgoVetoedMultiplicityOne_;
    std::array<std::unique_ptr<TH1D>, 5> siliconBgoVetoedByRing_;
    std::unique_ptr<TH2I> siliconEnergyVsGermaniumEnergy_;
    std::unique_ptr<TH1D> multiplicityAfterBgoVeto_;
    std::unique_ptr<TH1D> foldValidMultiplicity_;
    std::unique_ptr<TH1D> foldValidEnergy_;
    std::unique_ptr<TH1D> foldInclusiveEnergy_;
    std::unique_ptr<TH2I> energyVsTime_;
    bool calibrated_ = false;
    double runningTimeMaxSeconds_ = 1.0;
};

#endif
