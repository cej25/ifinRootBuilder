#ifndef RUN_CALIBRATION_MANAGER_H
#define RUN_CALIBRATION_MANAGER_H

#include "GermaniumCalibration.h"

#include <string>
#include <unordered_map>
#include <vector>

class RunCalibrationManager {
public:
    void addGlobalCalFile(const std::string& fileName);
    void addGlobalMcalFile(const std::string& fileName);
    void addRunCalFile(unsigned int firstRun, unsigned int lastRun,
                       const std::string& fileName);
    void addRunMcalFile(unsigned int firstRun, unsigned int lastRun,
                        const std::string& fileName);
    void addRunByRunCalFile(const std::string& fileName);
    void addRunByRunMcalFile(const std::string& fileName);

    bool empty() const;
    bool usesRunDependentCalibration() const;
    const GermaniumCalibration* calibrationForRun(unsigned int run) const;

    static unsigned int runNumberFromFileName(const std::string& fileName);

private:
    struct Range {
        unsigned int firstRun;
        unsigned int lastRun;
        GermaniumCalibration calibration;
    };

    Range& rangeFor(unsigned int firstRun, unsigned int lastRun);
    enum class BaseMode { Global, Ranges };
    void requireCompatibleBaseMode(BaseMode addingMode) const;
    void addRunByRunFile(const std::string& fileName, bool piecewise);
    const GermaniumCalibration* baseCalibrationForRun(unsigned int run) const;
    void invalidateCombinedCalibrations();

    GermaniumCalibration globalCalibration_;
    std::vector<Range> ranges_;
    std::unordered_map<unsigned int, GermaniumCalibration> runByRun_;
    mutable std::unordered_map<unsigned int, GermaniumCalibration>
        combinedByRun_;
};

#endif
