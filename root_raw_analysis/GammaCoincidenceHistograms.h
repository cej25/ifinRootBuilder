#ifndef GAMMA_COINCIDENCE_HISTOGRAMS_H
#define GAMMA_COINCIDENCE_HISTOGRAMS_H

#include "CoincidenceGateConfig.h"

#include <memory>
#include <vector>

class TDirectory;
class TH1D;
class TH2I;

class GammaCoincidenceHistograms {
public:
    GammaCoincidenceHistograms();
    ~GammaCoincidenceHistograms();

    GammaCoincidenceHistograms(const GammaCoincidenceHistograms&) = delete;
    GammaCoincidenceHistograms& operator=(
        const GammaCoincidenceHistograms&) = delete;

    void fillEvent(const std::vector<double>& gammaEnergies,
                   bool siliconCoincident);
    void configureGates(
        const std::vector<CoincidenceGateDefinition>& gates);
    void configureDoubleGates(
        const std::vector<DoubleCoincidenceGateDefinition>& gates);
    void setCalibratedEnergyAxes();
    void merge(const GammaCoincidenceHistograms& other);
    void write(TDirectory& parentDirectory) const;

private:
    void fillProjection(double gateEnergy, double projectedEnergy,
                        bool siliconCoincident);

    struct GateHistogram {
        CoincidenceGateDefinition definition;
        std::unique_ptr<TH1D> spectrum;
    };

    struct DoubleGateHistograms {
        DoubleCoincidenceGateDefinition definition;
        std::unique_ptr<TH2I> matrix;
        std::unique_ptr<TH1D> spectrum;
    };

    void fillDoubleGates(const std::vector<double>& gammaEnergies);

    std::unique_ptr<TH2I> gammaGamma_;
    std::unique_ptr<TH2I> gammaGammaSiliconCoincident_;
    std::vector<GateHistogram> gates_;
    std::vector<DoubleGateHistograms> doubleGates_;
    bool calibrated_ = false;
};

#endif
