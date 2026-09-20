#ifndef FOLD_STATISTICS_H
#define FOLD_STATISTICS_H

#include <array>
#include <cstdint>
#include <iosfwd>

class FoldStatistics {
public:
    bool record(unsigned int germaniumMultiplicity,
                unsigned int bgoMultiplicity,
                bool allBgoHitsMatched);
    void print(std::ostream& output) const;

private:
    static constexpr unsigned int kMinimumReportedGeMultiplicity = 1;
    static constexpr unsigned int kMaximumReportedGeMultiplicity = 4;
    static constexpr unsigned int kMaximumExplicitBgoMultiplicity = 10;
    static constexpr unsigned int kOverflowBgoCategory = 11;

    std::uint64_t totalEvents_ = 0;
    std::uint64_t validFoldEvents_ = 0;
    std::uint64_t invalidFoldEvents_ = 0;
    std::uint64_t validAndWellMatchedEvents_ = 0;
    std::array<std::array<std::uint64_t, 12>, 4> bgoByGermanium_{};
    std::array<std::uint64_t, 4> foldValidByGermanium_{};
};

#endif
