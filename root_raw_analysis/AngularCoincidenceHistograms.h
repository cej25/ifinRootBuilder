#ifndef ANGULAR_COINCIDENCE_HISTOGRAMS_H
#define ANGULAR_COINCIDENCE_HISTOGRAMS_H

#include "CoincidenceGateConfig.h"

#include <memory>
#include <vector>

class TDirectory;
class TH1D;
class TH2I;

class AngularCoincidenceHistograms {
public:
    AngularCoincidenceHistograms();
    ~AngularCoincidenceHistograms();

    AngularCoincidenceHistograms(const AngularCoincidenceHistograms&) = delete;
    AngularCoincidenceHistograms& operator=(
        const AngularCoincidenceHistograms&) = delete;

    void fillEvent(const std::vector<double>& gammaEnergies,
                   const std::vector<unsigned short>& gammaIDs);
    void configureGates(
        const std::vector<CoincidenceGateDefinition>& forwardGates,
        const std::vector<CoincidenceGateDefinition>& backwardGates);
    void setCalibratedEnergyAxes();
    void merge(const AngularCoincidenceHistograms& other);
    void write(TDirectory& gammaCoincidenceDirectory) const;

private:
    struct GateHistograms {
        CoincidenceGateDefinition definition;
        std::unique_ptr<TH1D> spectrum;
    };

    void fillOrientation(double allEnergy, double angularEnergy,
                         unsigned short angularID);
    static void configureGateSet(
        const std::vector<CoincidenceGateDefinition>& definitions,
        const char* suffix, const char* angleTitle,
        bool calibrated, std::vector<GateHistograms>& destination);

    std::unique_ptr<TH2I> allVsForward_;
    std::unique_ptr<TH2I> allVsBackward_;
    std::unique_ptr<TH1D> allProjectionForward_;
    std::unique_ptr<TH1D> allProjectionBackward_;
    std::vector<GateHistograms> forwardGates_;
    std::vector<GateHistograms> backwardGates_;
    bool calibrated_ = false;
};

#endif
