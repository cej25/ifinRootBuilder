#ifndef DRIFT_SPECTRUM_BUILDER_H
#define DRIFT_SPECTRUM_BUILDER_H

#include "GermaniumCalibration.h"

#include <string>
#include <vector>

class DriftSpectrumBuilder {
public:
    void addBlock(const std::string& series,
                  unsigned int firstRun,
                  unsigned int lastRun,
                  const std::string& inputPattern);
    void addCalFile(const std::string& fileName);
    void addMcalFile(const std::string& fileName);
    void addBreakBefore(const std::string& series, unsigned int run);
    void setGroupSize(unsigned int groupSize);
    void setThreadCount(unsigned int threadCount);
    void setGermaniumTimeGateEnabled(bool enabled);
    int run(const std::string& outputFileName) const;

private:
    struct Block {
        std::string series;
        unsigned int firstRun = 0;
        unsigned int lastRun = 0;
        std::string inputPattern;
    };

    struct BreakPoint {
        std::string series;
        unsigned int run = 0;
    };

    unsigned int groupSize_ = 10;
    unsigned int threadCount_ = 4;
    bool germaniumTimeGateEnabled_ = true;
    GermaniumCalibration calibration_;
    std::vector<Block> blocks_;
    std::vector<BreakPoint> breakPoints_;
};

#endif
