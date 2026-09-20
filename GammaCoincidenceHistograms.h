#ifndef GAMMA_COINCIDENCE_HISTOGRAMS_H
#define GAMMA_COINCIDENCE_HISTOGRAMS_H

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
    void setCalibratedEnergyAxes();
    void write(TDirectory& parentDirectory) const;

private:
    enum class ProjectionWindow {
        None,
        Prompt,
        LowerSideband,
        UpperSideband
    };

    ProjectionWindow projectionWindow(double energy) const;
    void fillProjection(double gateEnergy, double projectedEnergy,
                        bool siliconCoincident);

    std::unique_ptr<TH2I> gammaGamma_;
    std::unique_ptr<TH2I> gammaGammaSiliconCoincident_;
    std::unique_ptr<TH1D> backgroundSubtractedProjectionSiliconCoincident_;
};

#endif
