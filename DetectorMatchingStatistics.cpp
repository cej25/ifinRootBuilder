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
    const std::vector<unsigned short>& germaniumLUTs,
    const std::vector<unsigned short>& bgoLUTs)
{
    ++events_;

    std::unordered_map<unsigned short, unsigned int> unmatchedGermanium;
    bool repeatedGermaniumLUT = false;
    for (const unsigned short lut : germaniumLUTs) {
        unsigned int& count = unmatchedGermanium[lut];
        ++count;
        if (count == 2) {
            repeatedGermaniumLUT = true;
        }
    }

    unsigned int extraBgoMultiplicity = 0;
    for (const unsigned short lut : bgoLUTs) {
        auto germanium = unmatchedGermanium.find(lut);
        if (germanium != unmatchedGermanium.end() && germanium->second > 0) {
            --germanium->second;
        } else {
            ++extraBgoMultiplicity;
        }
    }

    if (repeatedGermaniumLUT) {
        ++eventsWithRepeatedGermaniumLUT_;
    }

    return {extraBgoMultiplicity, repeatedGermaniumLUT};
}

void DetectorMatchingStatistics::print(std::ostream& output) const
{
    output << "\n=== Repeated-germanium diagnostic ===\n"
           << std::fixed << std::setprecision(2)
           << "Events with multiple accepted Ge hits in the same LUT: "
           << eventsWithRepeatedGermaniumLUT_ << " / " << events_ << " ("
           << percentage(eventsWithRepeatedGermaniumLUT_, events_)
           << "%)\n";
}
