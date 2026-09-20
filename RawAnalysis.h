#ifndef RAW_ANALYSIS_H
#define RAW_ANALYSIS_H

#include "AngularCoincidenceHistograms.h"
#include "DetectorHistograms.h"
#include "GammaCoincidenceHistograms.h"
#include "GermaniumCalibration.h"
#include "GermaniumConditionHistograms.h"
#include "GermaniumDetectorHistograms.h"
#include "IndividualDetectorHistograms.h"
#include "RunCalibrationManager.h"

#include <cstddef>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

class RawAnalysis {
public:
    RawAnalysis();

    void addCalFile(const std::string& fileName);
    void addMcalFile(const std::string& fileName);
    void addRunCalFile(unsigned int firstRun, unsigned int lastRun,
                       const std::string& fileName);
    void addRunMcalFile(unsigned int firstRun, unsigned int lastRun,
                        const std::string& fileName);
    void setDiagnosticsEnabled(bool enabled);
    void excludeGermaniumLUT(unsigned int detectorLUT);

    int run(const std::vector<std::string>& inputPatterns,
            const std::string& outputFileName);

private:
    bool diagnosticsEnabled_ = true;
    std::unordered_set<unsigned short> excludedGermaniumLUTs_;
    std::vector<DetectorHistograms> detectorHistograms_;
    std::unordered_map<unsigned short, std::size_t> detectorIndex_;
    GammaCoincidenceHistograms gammaCoincidences_;
    AngularCoincidenceHistograms angularCoincidences_;
    RunCalibrationManager calibrationManager_;
    GermaniumDetectorHistograms individualGermaniumHistograms_;
    GermaniumConditionHistograms germaniumConditionHistograms_;
    IndividualDetectorHistograms individualSiliconHistograms_;
    IndividualDetectorHistograms individualBgoHistograms_;
    IndividualDetectorHistograms individualLabrHistograms_;
};

#endif
