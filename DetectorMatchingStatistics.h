#ifndef DETECTOR_MATCHING_STATISTICS_H
#define DETECTOR_MATCHING_STATISTICS_H

#include <cstdint>
#include <iosfwd>
#include <vector>

class DetectorMatchingStatistics {
public:
    struct EventResult {
        unsigned int extraBgoMultiplicity;
        bool repeatedGermaniumLUT;
    };

    EventResult record(
        const std::vector<unsigned short>& germaniumLUTs,
        const std::vector<unsigned short>& bgoLUTs);
    void print(std::ostream& output) const;

private:
    std::uint64_t events_ = 0;
    std::uint64_t eventsWithRepeatedGermaniumLUT_ = 0;
};

#endif
