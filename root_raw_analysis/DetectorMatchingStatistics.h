#ifndef DETECTOR_MATCHING_STATISTICS_H
#define DETECTOR_MATCHING_STATISTICS_H

#include <cstdint>
#include <iosfwd>
#include <vector>

class DetectorMatchingStatistics {
public:
    struct EventResult {
        unsigned int extraBgoMultiplicity;
        bool repeatedGermaniumID;
    };

    EventResult record(
        const std::vector<unsigned short>& germaniumIDs,
        const std::vector<unsigned short>& bgoIDs);
    void merge(const DetectorMatchingStatistics& other);
    void print(std::ostream& output) const;

private:
    std::uint64_t events_ = 0;
    std::uint64_t eventsWithRepeatedGermaniumID_ = 0;
};

#endif
