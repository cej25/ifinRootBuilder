#ifndef CALIBRATION_PEAK_FITTER_H
#define CALIBRATION_PEAK_FITTER_H

#include "GermaniumCalibration.h"

#include <string>
#include <vector>

class CalibrationPeakFitter {
public:
    void addSeedCalFile(const std::string& fileName);
    void addSeedMcalFile(const std::string& fileName);
    void loadPeakDefinitions(const std::string& fileName);
    void loadAdcOverrides(const std::string& fileName);
    int run(const std::string& inputFileName,
            const std::string& outputFileName) const;

private:
    struct PeakDefinition {
        std::string group;
        double energyKeV = 0.0;
        std::string label;
        double searchHalfWidthKeV = 8.0;
        double fitHalfWidthKeV = 5.0;
    };

    struct AdcOverride {
        std::string group;
        unsigned int detectorID = 0;
        std::string label;
        double centreAdc = 0.0;
        double searchHalfWidthAdc = 0.0;
        double fitHalfWidthAdc = 0.0;
        bool skip = false;
    };

    GermaniumCalibration seedCalibration_;
    std::vector<PeakDefinition> peaks_;
    std::vector<AdcOverride> adcOverrides_;
};

#endif
