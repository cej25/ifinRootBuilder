#ifndef RAW_ANALYSIS_H
#define RAW_ANALYSIS_H

#include "AngularCoincidenceHistograms.h"
#include "DetectorHistograms.h"
#include "GammaCoincidenceHistograms.h"
#include "GermaniumCalibration.h"
#include "GermaniumConditionHistograms.h"
#include "GermaniumDetectorHistograms.h"
#include "IndividualDetectorHistograms.h"
#include "DetectorMatchingStatistics.h"
#include "FoldStatistics.h"
#include "GateStatistics.h"
#include "MultiplicityStatistics.h"
#include "RunCalibrationManager.h"

#include <RtypesCore.h>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

class RawAnalysis {
public:
    RawAnalysis();
    ~RawAnalysis();

    void addCalFile(const std::string& fileName);
    void addMcalFile(const std::string& fileName);
    void addRunCalFile(unsigned int firstRun, unsigned int lastRun,
                       const std::string& fileName);
    void addRunMcalFile(unsigned int firstRun, unsigned int lastRun,
                        const std::string& fileName);
    void setDiagnosticsEnabled(bool enabled);
    void setThreadCount(unsigned int threadCount);
    void excludeGermaniumID(unsigned int detectorID);

    int run(const std::vector<std::string>& inputPatterns,
            const std::string& outputFileName);

private:
    void configureCalibratedAxes();
    void processReader(class TTreeReader& reader, bool showProgress = true);
    void merge(const RawAnalysis& other);

    bool diagnosticsEnabled_ = true;
    unsigned int threadCount_ = 1;
    std::unordered_set<unsigned short> excludedGermaniumIDs_;
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
    std::unique_ptr<class TH1D> eventMultiplicity_;
    std::unique_ptr<class TH1D> unknownDetectorTypes_;
    MultiplicityStatistics combinedStatistics_;
    FoldStatistics foldStatistics_;
    GateStatistics gateStatistics_;
    DetectorMatchingStatistics detectorMatchingStatistics_;
    ULong64_t processedEvents_ = 0;
    ULong64_t malformedEvents_ = 0;
    ULong64_t unknownHits_ = 0;
    ULong64_t nonzeroPsdHits_ = 0;
    ULong64_t missingCalibrationEvents_ = 0;
    ULong64_t outOfRangeCalibrationEvents_ = 0;
    ULong64_t nonFiniteCalibrationEvents_ = 0;
    ULong_t firstAbsoluteTime_ = static_cast<ULong_t>(-1);
    ULong_t lastAbsoluteTime_ = 0;
};

#endif
