#include "DetectorMatchingStatistics.h"

#include <iomanip>
#include <ostream>
#include <unordered_map>

namespace {

double percentage(std::uint64_t count, std::uint64_t total)
{
    return total == 0 ? 0.0
                      : 100.0 * static_cast<double>(count) /
                            static_cast<double>(total);
}

} // namespace

DetectorMatchingStatistics::EventResult DetectorMatchingStatistics::record(
    const std::vector<unsigned short>& germaniumIDs,
    const std::vector<unsigned short>& bgoIDs)
{
    ++events_;

    std::unordered_map<unsigned short, unsigned int> unmatchedGermanium;
    bool repeatedGermaniumID = false;
    for (const unsigned short id : germaniumIDs) {
        unsigned int& count = unmatchedGermanium[id];
        ++count;
        if (count == 2) {
            repeatedGermaniumID = true;
        }
    }

    unsigned int extraBgoMultiplicity = 0;
    for (const unsigned short id : bgoIDs) {
        auto germanium = unmatchedGermanium.find(id);
        if (germanium != unmatchedGermanium.end() && germanium->second > 0) {
            --germanium->second;
        } else {
            ++extraBgoMultiplicity;
        }
    }

    if (repeatedGermaniumID) {
        ++eventsWithRepeatedGermaniumID_;
    }

    return {extraBgoMultiplicity, repeatedGermaniumID};
}

void DetectorMatchingStatistics::merge(
    const DetectorMatchingStatistics& other)
{
    events_ += other.events_;
    eventsWithRepeatedGermaniumID_ +=
        other.eventsWithRepeatedGermaniumID_;
}

void DetectorMatchingStatistics::print(std::ostream& output) const
{
    output << "\n=== Repeated-germanium diagnostic ===\n"
           << std::fixed << std::setprecision(2)
           << "Events with multiple accepted Ge hits in the same ID: "
           << eventsWithRepeatedGermaniumID_ << " / " << events_ << " ("
           << percentage(eventsWithRepeatedGermaniumID_, events_)
           << "%)\n";
}
