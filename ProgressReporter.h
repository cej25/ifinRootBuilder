#ifndef PROGRESS_REPORTER_H
#define PROGRESS_REPORTER_H

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <mutex>
#include <thread>

class ProgressReporter {
public:
    explicit ProgressReporter(std::uint64_t totalEvents,
                              bool enabled = true);
    ~ProgressReporter();

    ProgressReporter(const ProgressReporter&) = delete;
    ProgressReporter& operator=(const ProgressReporter&) = delete;

    void add(std::uint64_t events);
    void finish();

private:
    void run();
    void print(bool final);

    const std::uint64_t totalEvents_;
    const bool enabled_;
    const std::chrono::steady_clock::time_point start_;
    std::atomic<std::uint64_t> processedEvents_{0};
    std::mutex mutex_;
    std::condition_variable wakeup_;
    bool stopping_ = false;
    bool finished_ = false;
    std::thread reporterThread_;
};

#endif
