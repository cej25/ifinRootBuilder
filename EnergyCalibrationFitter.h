#ifndef ENERGY_CALIBRATION_FITTER_H
#define ENERGY_CALIBRATION_FITTER_H

#include <string>

class EnergyCalibrationFitter {
public:
    void loadModelOverrides(const std::string& fileName);
    void setMaximumEnergyKeV(double energyKeV);

    int run(const std::string& peakFitFileName,
            const std::string& outputCalibrationFileName,
            const std::string& diagnosticFileName) const;

private:
    struct PiecewiseOverride {
        unsigned int detectorID = 0;
        double boundaryEnergyKeV = 0.0;
        unsigned int lowerDegree = 1;
        unsigned int upperDegree = 1;
    };

    double maximumEnergyKeV_ = 2048.0;
    std::string overrideFileName_;
};

#endif
