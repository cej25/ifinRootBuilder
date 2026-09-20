#ifndef RUN_CALIBRATION_MANAGER_H
#define RUN_CALIBRATION_MANAGER_H

#include "GermaniumCalibration.h"

#include <string>
#include <vector>

class RunCalibrationManager {
public:
    void addGlobalCalFile(const std::string& fileName);
    void addGlobalMcalFile(const std::string& fileName);
    void addRunCalFile(unsigned int firstRun, unsigned int lastRun,
                       const std::string& fileName);
    void addRunMcalFile(unsigned int firstRun, unsigned int lastRun,
                        const std::string& fileName);

    bool empty() const;
    bool usesRunRanges() const;
    const GermaniumCalibration* calibrationForRun(unsigned int run) const;

    static unsigned int runNumberFromFileName(const std::string& fileName);

private:
    struct Range {
        unsigned int firstRun;
        unsigned int lastRun;
        GermaniumCalibration calibration;
    };

    Range& rangeFor(unsigned int firstRun, unsigned int lastRun);
    void requireCompatibleMode(bool addingRange) const;

    GermaniumCalibration globalCalibration_;
    std::vector<Range> ranges_;
};

#endif
