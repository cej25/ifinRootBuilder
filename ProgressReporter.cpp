#include "ProgressReporter.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

namespace {

std::string formatCount(std::uint64_t value)
{
    std::string text = std::to_string(value);
    for (std::ptrdiff_t position =
             static_cast<std::ptrdiff_t>(text.size()) - 3;
         position > 0; position -= 3) {
        text.insert(static_cast<std::size_t>(position), 1, ',');
    }
    return text;
}

std::string formatDuration(double seconds)
{
    if (!std::isfinite(seconds) || seconds < 0.0) return "--:--:--";
    const auto totalSeconds = static_cast<std::uint64_t>(seconds);
    const auto days = totalSeconds / 86400ULL;
    const auto hours = (totalSeconds / 3600ULL) % 24ULL;
    const auto minutes = (totalSeconds / 60ULL) % 60ULL;
    const auto remainder = totalSeconds % 60ULL;

    std::ostringstream text;
    if (days != 0) text << days << "d ";
    text << std::setfill('0') << std::setw(2) << hours << ':'
         << std::setw(2) << minutes << ':' << std::setw(2) << remainder;
    return text.str();
}

std::string formatRate(double eventsPerSecond)
{
    std::ostringstream text;
    if (eventsPerSecond >= 1.0e6) {
        text << std::fixed << std::setprecision(2)
             << eventsPerSecond / 1.0e6 << " M events/s";
    } else if (eventsPerSecond >= 1.0e3) {
        text << std::fixed << std::setprecision(1)
             << eventsPerSecond / 1.0e3 << " k events/s";
    } else {
        text << std::fixed << std::setprecision(0)
             << eventsPerSecond << " events/s";
    }
    return text.str();
}

} // namespace

ProgressReporter::ProgressReporter(std::uint64_t totalEvents, bool enabled)
    : totalEvents_(totalEvents),
      enabled_(enabled && totalEvents != 0),
      start_(std::chrono::steady_clock::now())
{
    if (enabled_) reporterThread_ = std::thread(&ProgressReporter::run, this);
}

ProgressReporter::~ProgressReporter()
{
    finish();
}

void ProgressReporter::add(std::uint64_t events)
{
    if (enabled_) processedEvents_.fetch_add(events, std::memory_order_relaxed);
}

void ProgressReporter::finish()
{
    if (!enabled_) return;
    {
        std::lock_guard<std::mutex> lock(mutex_);
        if (finished_) return;
        stopping_ = true;
    }
    wakeup_.notify_one();
    if (reporterThread_.joinable()) reporterThread_.join();
    print(true);
    finished_ = true;
}

void ProgressReporter::run()
{
    std::unique_lock<std::mutex> lock(mutex_);
    while (!wakeup_.wait_for(lock, std::chrono::seconds(5),
                             [this] { return stopping_; })) {
        lock.unlock();
        print(false);
        lock.lock();
    }
}

void ProgressReporter::print(bool final)
{
    const std::uint64_t processed = std::min(
        processedEvents_.load(std::memory_order_relaxed), totalEvents_);
    const double elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start_).count();
    const double rate = elapsed > 0.0 ? processed / elapsed : 0.0;
    const double percentage = totalEvents_ != 0
        ? 100.0 * static_cast<double>(processed) / totalEvents_
        : 100.0;

    std::ostringstream line;
    line << "Progress: " << formatCount(processed) << " / "
         << formatCount(totalEvents_) << " (" << std::fixed
         << std::setprecision(2) << percentage << "%) | "
         << formatRate(rate) << " | elapsed " << formatDuration(elapsed)
         << " | ETA ";
    if (processed >= totalEvents_) {
        line << "00:00:00";
    } else if (rate > 0.0) {
        line << formatDuration((totalEvents_ - processed) / rate);
    } else {
        line << "--:--:--";
    }

    std::string output = line.str();
    if (output.size() < 120) output.append(120 - output.size(), ' ');
    std::cout << '\r' << output;
    if (final) std::cout << '\n';
    std::cout << std::flush;
}
