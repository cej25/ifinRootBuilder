#ifndef ANALYSIS_TREE_BUILDER_H
#define ANALYSIS_TREE_BUILDER_H

#include "AnalysisConfig.h"
#include "RunCalibrationManager.h"

#include <string>
#include <unordered_set>
#include <vector>

class AnalysisTreeBuilder {
public:
    void addCalFile(const std::string& fileName);
    void addMcalFile(const std::string& fileName);
    void addRunCalFile(unsigned int firstRun, unsigned int lastRun,
                       const std::string& fileName);
    void addRunMcalFile(unsigned int firstRun, unsigned int lastRun,
                        const std::string& fileName);
    void excludeGermaniumID(unsigned int detectorID);
    void setDiagnosticsEnabled(bool enabled);
    void setThreadCount(unsigned int threadCount);

    int run(const std::vector<std::string>& inputPatterns,
            const std::string& outputDirectory) const;

private:
    struct FileResult {
        unsigned long long inputEvents = 0;
        unsigned long long writtenEvents = 0;
        unsigned long long malformedEvents = 0;
        unsigned long long calibrationRejectedEvents = 0;
    };

    FileResult processFile(const std::string& inputFileName,
                           const std::string& outputFileName) const;

    bool diagnosticsEnabled_ = true;
    unsigned int threadCount_ = config::kDefaultThreadCount;
    std::unordered_set<unsigned short> excludedGermaniumIDs_;
    RunCalibrationManager calibrationManager_;
};

#endif
