#ifndef ANALYSIS_TREE_ANALYSIS_H
#define ANALYSIS_TREE_ANALYSIS_H

#include "AngularCoincidenceHistograms.h"
#include "DetectorHistograms.h"
#include "DetectorMatchingStatistics.h"
#include "FoldStatistics.h"
#include "GammaCoincidenceHistograms.h"
#include "GateStatistics.h"
#include "GermaniumConditionHistograms.h"
#include "GermaniumDetectorHistograms.h"
#include "IndividualDetectorHistograms.h"
#include "MultiplicityStatistics.h"
#include "RunningTimeMap.h"

#include <RtypesCore.h>

#include <cstddef>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

class AnalysisTreeAnalysis {
public:
    AnalysisTreeAnalysis();
    ~AnalysisTreeAnalysis();

    void setDiagnosticsEnabled(bool enabled);
    void setProgressEnabled(bool enabled);
    void setThreadCount(unsigned int threadCount);
    int run(const std::vector<std::string>& inputPatterns,
            const std::string& outputFileName);

private:
    void processReader(class TTreeReader& reader,
                       class ProgressReporter* progressReporter);
    void merge(const AnalysisTreeAnalysis& other);
    void configureCalibratedAxes();
    int writeOutput(const std::string& outputFileName);
    void printDiagnostics() const;

    bool diagnosticsEnabled_ = true;
    bool progressEnabled_ = true;
    bool calibratedEnergySeen_ = false;
    bool uncalibratedEnergySeen_ = false;
    unsigned int threadCount_ = 1;
    std::vector<DetectorHistograms> detectorHistograms_;
    std::unordered_map<unsigned short, std::size_t> detectorIndex_;
    GammaCoincidenceHistograms gammaCoincidences_;
    AngularCoincidenceHistograms angularCoincidences_;
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
    RunningTimeMap runningTimeMap_;
    ULong64_t processedEvents_ = 0;
    ULong64_t malformedEvents_ = 0;
    ULong64_t unknownHits_ = 0;
    ULong64_t nonzeroPsdHits_ = 0;
    ULong64_t firstAbsoluteTime_ = static_cast<ULong64_t>(-1);
    ULong64_t lastAbsoluteTime_ = 0;
};

#endif
