#ifndef RUNNING_TIME_MAP_H
#define RUNNING_TIME_MAP_H

#include <cstdint>
#include <iosfwd>
#include <string>
#include <unordered_map>
#include <vector>

class RunningTimeMap {
public:
    struct FileRange {
        std::string fileName;
        std::uint64_t firstAbsoluteTime = 0;
        std::uint64_t lastAbsoluteTime = 0;
        double offsetSeconds = 0.0;
    };

    void addFile(const std::string& fileName,
                 std::uint64_t firstAbsoluteTime,
                 std::uint64_t lastAbsoluteTime);
    double runningTimeSeconds(const std::string& fileName,
                              std::uint64_t absoluteTime) const;
    const FileRange& rangeForFile(const std::string& fileName) const;
    static double runningTimeSeconds(const FileRange& range,
                                     std::uint64_t absoluteTime);
    static bool tryRunningTimeSeconds(const FileRange& range,
                                      std::uint64_t absoluteTime,
                                      double& runningTimeSeconds) noexcept;
    double totalSeconds() const;
    std::size_t fileCount() const;
    void print(std::ostream& output) const;

private:
    std::vector<FileRange> files_;
    std::unordered_map<std::string, std::size_t> fileIndex_;
    double totalSeconds_ = 0.0;
};

#endif
