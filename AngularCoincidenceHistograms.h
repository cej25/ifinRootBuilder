#ifndef ANGULAR_COINCIDENCE_HISTOGRAMS_H
#define ANGULAR_COINCIDENCE_HISTOGRAMS_H

#include <memory>
#include <vector>

class TDirectory;
class TH1D;

class AngularCoincidenceHistograms {
public:
    AngularCoincidenceHistograms();
    ~AngularCoincidenceHistograms();

    AngularCoincidenceHistograms(const AngularCoincidenceHistograms&) = delete;
    AngularCoincidenceHistograms& operator=(
        const AngularCoincidenceHistograms&) = delete;

    void fillEvent(const std::vector<double>& gammaEnergies,
                   const std::vector<unsigned short>& gammaIDs);
    void setCalibratedEnergyAxes();
    void merge(const AngularCoincidenceHistograms& other);
    void write(TDirectory& gammaCoincidenceDirectory) const;

private:
    struct GateHistograms {
        const char* name;
        double minimum;
        double maximumExclusive;
        std::unique_ptr<TH1D> forward;
        std::unique_ptr<TH1D> backward;
    };

    std::vector<GateHistograms> gates_;
};

#endif
