#ifndef MULTIPLICITY_STATISTICS_H
#define MULTIPLICITY_STATISTICS_H

#include "AnalysisConfig.h"

#include <array>
#include <cstdint>
#include <iosfwd>
#include <string>
#include <vector>

class MultiplicityStatistics {
public:
    static constexpr unsigned int kMaximumReportedEventFold = 8;

    MultiplicityStatistics();

    void record(const std::vector<unsigned int>& detectorMultiplicities);
    void print(std::ostream& output, const std::string& label) const;

private:
    static std::size_t indexForType(unsigned short detectorType);
    void printDetectorBreakdown(std::ostream& output,
                                unsigned short detectorType) const;

    std::array<std::vector<std::uint64_t>, config::kDetectors.size()>
        detectorFoldCounts_;
    std::array<std::uint64_t, kMaximumReportedEventFold + 1>
        eventCounts_{};
    std::array<
        std::array<std::uint64_t, config::kDetectors.size()>,
        kMaximumReportedEventFold + 1> detectorHitTotals_{};
};

#endif
