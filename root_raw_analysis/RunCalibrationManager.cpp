#include "RunCalibrationManager.h"

#include <filesystem>
#include <limits>
#include <regex>
#include <stdexcept>
#include <utility>

void RunCalibrationManager::requireCompatibleBaseMode(
    BaseMode addingMode) const
{
    const bool hasGlobal = !globalCalibration_.empty();
    const bool hasRanges = !ranges_.empty();
    if ((addingMode == BaseMode::Global && hasRanges) ||
        (addingMode == BaseMode::Ranges && hasGlobal)) {
        throw std::runtime_error(
            "Global and run-ranged main calibrations cannot be mixed");
    }
}

void RunCalibrationManager::invalidateCombinedCalibrations()
{
    combinedByRun_.clear();
}

void RunCalibrationManager::addGlobalCalFile(const std::string& fileName)
{
    requireCompatibleBaseMode(BaseMode::Global);
    globalCalibration_.addCalFile(fileName);
    invalidateCombinedCalibrations();
}

void RunCalibrationManager::addGlobalMcalFile(const std::string& fileName)
{
    requireCompatibleBaseMode(BaseMode::Global);
    globalCalibration_.addMcalFile(fileName);
    invalidateCombinedCalibrations();
}

RunCalibrationManager::Range& RunCalibrationManager::rangeFor(
    unsigned int firstRun, unsigned int lastRun)
{
    requireCompatibleBaseMode(BaseMode::Ranges);
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
    invalidateCombinedCalibrations();
}

void RunCalibrationManager::addRunMcalFile(
    unsigned int firstRun, unsigned int lastRun, const std::string& fileName)
{
    rangeFor(firstRun, lastRun).calibration.addMcalFile(fileName);
    invalidateCombinedCalibrations();
}

void RunCalibrationManager::addRunByRunCalFile(
    const std::string& fileName)
{
    addRunByRunFile(fileName, false);
}

void RunCalibrationManager::addRunByRunMcalFile(
    const std::string& fileName)
{
    addRunByRunFile(fileName, true);
}

void RunCalibrationManager::addRunByRunFile(
    const std::string& fileName, bool piecewise)
{
    auto loaded = piecewise
        ? GermaniumCalibration::loadRunByRunMcalFile(fileName)
        : GermaniumCalibration::loadRunByRunCalFile(fileName);

    for (auto& entry : loaded) {
        auto existing = runByRun_.find(entry.first);
        if (existing == runByRun_.end()) {
            runByRun_.emplace(entry.first, std::move(entry.second));
        } else {
            existing->second.appendStages(std::move(entry.second));
        }
    }
    invalidateCombinedCalibrations();
}

bool RunCalibrationManager::empty() const
{
    return globalCalibration_.empty() && ranges_.empty() && runByRun_.empty();
}

bool RunCalibrationManager::usesRunDependentCalibration() const
{
    return !ranges_.empty() || !runByRun_.empty();
}

const GermaniumCalibration* RunCalibrationManager::calibrationForRun(
    unsigned int run) const
{
    const GermaniumCalibration* base = baseCalibrationForRun(run);
    if (!runByRun_.empty() && base == nullptr) {
        throw std::runtime_error(
            "Run-by-run fine adjustments require a main calibration "
            "for run " + std::to_string(run));
    }
    const auto adjustment = runByRun_.find(run);
    if (adjustment == runByRun_.end()) {
        if (base != nullptr) {
            return base;
        }
        return nullptr;
    }

    const auto existing = combinedByRun_.find(run);
    if (existing != combinedByRun_.end()) {
        return &existing->second;
    }

    GermaniumCalibration combined;
    if (base != nullptr) {
        combined.appendStages(*base);
    }
    combined.appendStages(adjustment->second);
    return &combinedByRun_.emplace(run, std::move(combined)).first->second;
}

const GermaniumCalibration* RunCalibrationManager::baseCalibrationForRun(
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
