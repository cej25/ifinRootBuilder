#include "RunningTimeMap.h"

#include <iomanip>
#include <ostream>
#include <stdexcept>

void RunningTimeMap::addFile(const std::string& fileName,
                             std::uint64_t firstAbsoluteTime,
                             std::uint64_t lastAbsoluteTime)
{
    if (lastAbsoluteTime < firstAbsoluteTime) {
        throw std::runtime_error(
            "absoluteTime decreases within file '" + fileName + "'");
    }
    if (fileIndex_.count(fileName) != 0) {
        throw std::runtime_error(
            "duplicate input filename in running-time map: '" +
            fileName + "'");
    }

    const std::size_t index = files_.size();
    files_.push_back(FileRange{fileName, firstAbsoluteTime,
                               lastAbsoluteTime, totalSeconds_});
    fileIndex_.emplace(fileName, index);
    totalSeconds_ += static_cast<double>(lastAbsoluteTime - firstAbsoluteTime);
}

double RunningTimeMap::runningTimeSeconds(
    const std::string& fileName, std::uint64_t absoluteTime) const
{
    return runningTimeSeconds(rangeForFile(fileName), absoluteTime);
}

const RunningTimeMap::FileRange& RunningTimeMap::rangeForFile(
    const std::string& fileName) const
{
    const auto found = fileIndex_.find(fileName);
    if (found == fileIndex_.end()) {
        throw std::runtime_error(
            "no running-time information for file '" + fileName + "'");
    }
    return files_.at(found->second);
}

double RunningTimeMap::runningTimeSeconds(
    const FileRange& range, std::uint64_t absoluteTime)
{
    if (absoluteTime < range.firstAbsoluteTime ||
        absoluteTime > range.lastAbsoluteTime) {
        throw std::runtime_error(
            "absoluteTime lies outside the recorded range for file '" +
            range.fileName + "'");
    }
    return range.offsetSeconds +
        static_cast<double>(absoluteTime - range.firstAbsoluteTime);
}

bool RunningTimeMap::tryRunningTimeSeconds(
    const FileRange& range, std::uint64_t absoluteTime,
    double& runningTimeSeconds) noexcept
{
    if (absoluteTime < range.firstAbsoluteTime ||
        absoluteTime > range.lastAbsoluteTime) {
        return false;
    }
    runningTimeSeconds = range.offsetSeconds +
        static_cast<double>(absoluteTime - range.firstAbsoluteTime);
    return true;
}

double RunningTimeMap::totalSeconds() const
{
    return totalSeconds_;
}

std::size_t RunningTimeMap::fileCount() const
{
    return files_.size();
}

void RunningTimeMap::print(std::ostream& output) const
{
    output << "\n=== Input-file absolute-time ranges ===\n"
           << "file  first [s]  last [s]  duration [s]  cumulative end [s]\n";
    for (std::size_t index = 0; index < files_.size(); ++index) {
        const FileRange& file = files_[index];
        const double duration = static_cast<double>(
            file.lastAbsoluteTime - file.firstAbsoluteTime);
        output << std::setw(4) << index + 1 << "  "
               << file.firstAbsoluteTime << "  "
               << file.lastAbsoluteTime << "  "
               << duration << "  " << file.offsetSeconds + duration
               << "  " << file.fileName << '\n';
    }
    output << "Total running time: " << totalSeconds_ << " s\n";
}
