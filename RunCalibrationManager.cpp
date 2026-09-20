#include "RunCalibrationManager.h"

#include <filesystem>
#include <limits>
#include <regex>
#include <stdexcept>

void RunCalibrationManager::requireCompatibleMode(bool addingRange) const
{
    if (addingRange && !globalCalibration_.empty()) {
        throw std::runtime_error(
            "Global --cal/--mcal options cannot be mixed with run-ranged calibrations");
    }
    if (!addingRange && !ranges_.empty()) {
        throw std::runtime_error(
            "Run-ranged calibrations cannot be mixed with global --cal/--mcal options");
    }
}

void RunCalibrationManager::addGlobalCalFile(const std::string& fileName)
{
    requireCompatibleMode(false);
    globalCalibration_.addCalFile(fileName);
}

void RunCalibrationManager::addGlobalMcalFile(const std::string& fileName)
{
    requireCompatibleMode(false);
    globalCalibration_.addMcalFile(fileName);
}

RunCalibrationManager::Range& RunCalibrationManager::rangeFor(
    unsigned int firstRun, unsigned int lastRun)
{
    requireCompatibleMode(true);
    if (firstRun > lastRun) {
        throw std::runtime_error(
            "Calibration range start must not exceed its end");
    }

    for (Range& range : ranges_) {
        if (range.firstRun == firstRun && range.lastRun == lastRun) {
            return range;
        }
        if (firstRun <= range.lastRun && lastRun >= range.firstRun) {
            throw std::runtime_error(
                "Calibration run ranges overlap without being identical");
        }
    }

    ranges_.push_back({firstRun, lastRun, GermaniumCalibration{}});
    return ranges_.back();
}

void RunCalibrationManager::addRunCalFile(
    unsigned int firstRun, unsigned int lastRun, const std::string& fileName)
{
    rangeFor(firstRun, lastRun).calibration.addCalFile(fileName);
}

void RunCalibrationManager::addRunMcalFile(
    unsigned int firstRun, unsigned int lastRun, const std::string& fileName)
{
    rangeFor(firstRun, lastRun).calibration.addMcalFile(fileName);
}

bool RunCalibrationManager::empty() const
{
    return globalCalibration_.empty() && ranges_.empty();
}

bool RunCalibrationManager::usesRunRanges() const
{
    return !ranges_.empty();
}

const GermaniumCalibration* RunCalibrationManager::calibrationForRun(
    unsigned int run) const
{
    if (!ranges_.empty()) {
        for (const Range& range : ranges_) {
            if (run >= range.firstRun && run <= range.lastRun) {
                return &range.calibration;
            }
        }
        throw std::runtime_error(
            "No calibration was configured for run " + std::to_string(run));
    }
    return globalCalibration_.empty() ? nullptr : &globalCalibration_;
}

unsigned int RunCalibrationManager::runNumberFromFileName(
    const std::string& fileName)
{
    const std::string baseName =
        std::filesystem::path(fileName).filename().string();
    const std::regex pattern(R"((?:^|_)([0-9]+)\.root$)");
    std::smatch match;
    if (!std::regex_search(baseName, match, pattern)) {
        throw std::runtime_error(
            "Could not extract the run number from input file '" +
            fileName + "' (expected a name ending in _XXXXXX.root)");
    }

    const unsigned long parsed = std::stoul(match[1].str());
    if (parsed > std::numeric_limits<unsigned int>::max()) {
        throw std::runtime_error(
            "Run number is too large in input file '" + fileName + "'");
    }
    return static_cast<unsigned int>(parsed);
}
