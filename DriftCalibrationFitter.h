#ifndef DRIFT_CALIBRATION_FITTER_H
#define DRIFT_CALIBRATION_FITTER_H

#include <set>
#include <string>
#include <vector>

class DriftCalibrationFitter {
public:
    void loadPeakDefinitions(const std::string& fileName);
    void setMinimumPeakCount(unsigned int count);
    void setMinimumEnergySpan(double energySpanKeV);
    int run(const std::string& inputFileName,
            const std::string& outputPrefix,
            const std::string& diagnosticFileName) const;

private:
    struct PeakDefinition {
        double energyKeV = 0.0;
        std::string label;
        double searchHalfWidthKeV = 5.0;
        double fitHalfWidthKeV = 3.0;
        std::set<unsigned int> excludedDetectorIDs;
    };

    unsigned int minimumPeakCount_ = 2;
    double minimumEnergySpanKeV_ = 400.0;
    std::vector<PeakDefinition> peaks_;
};

#endif
