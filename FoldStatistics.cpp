#include "FoldStatistics.h"

#include <algorithm>
#include <iomanip>
#include <ostream>

namespace {

double percentage(std::uint64_t count, std::uint64_t total)
{
    return total == 0 ? 0.0
                      : 100.0 * static_cast<double>(count) /
                            static_cast<double>(total);
}

} // namespace

bool FoldStatistics::record(unsigned int germaniumMultiplicity,
                            unsigned int bgoMultiplicity,
                            bool allBgoHitsMatched)
{
    ++totalEvents_;
    const bool foldValid = bgoMultiplicity <= germaniumMultiplicity;
    if (foldValid) {
        ++validFoldEvents_;
        if (allBgoHitsMatched) {
            ++validAndWellMatchedEvents_;
        }
    } else {
        ++invalidFoldEvents_;
    }

    if (germaniumMultiplicity >= kMinimumReportedGeMultiplicity &&
        germaniumMultiplicity <= kMaximumReportedGeMultiplicity) {
        const unsigned int bgoCategory = std::min(
            bgoMultiplicity, kOverflowBgoCategory);
        ++bgoByGermanium_[germaniumMultiplicity - 1][bgoCategory];
        if (foldValid) {
            ++foldValidByGermanium_[germaniumMultiplicity - 1];
        }
    }
    return foldValid;
}

void FoldStatistics::merge(const FoldStatistics& other)
{
    totalEvents_ += other.totalEvents_;
    validFoldEvents_ += other.validFoldEvents_;
    invalidFoldEvents_ += other.invalidFoldEvents_;
    validAndWellMatchedEvents_ += other.validAndWellMatchedEvents_;
    for (std::size_t ge = 0; ge < bgoByGermanium_.size(); ++ge) {
        foldValidByGermanium_[ge] += other.foldValidByGermanium_[ge];
        for (std::size_t bgo = 0; bgo < bgoByGermanium_[ge].size(); ++bgo) {
            bgoByGermanium_[ge][bgo] += other.bgoByGermanium_[ge][bgo];
        }
    }
}

void FoldStatistics::print(std::ostream& output) const
{
    output << "\n=== FoldValid statistics ===\n"
           << "FoldValid is true when BGO multiplicity <= Ge multiplicity.\n"
           << "Fold-valid events: " << validFoldEvents_ << " / "
           << totalEvents_ << " (" << std::fixed << std::setprecision(2)
           << percentage(validFoldEvents_, totalEvents_) << "%)\n"
           << "Fold-invalid events (BGO multiplicity > Ge multiplicity): "
           << invalidFoldEvents_
           << " / " << totalEvents_ << " ("
           << percentage(invalidFoldEvents_, totalEvents_) << "%)\n"
           << "FoldValid events with every BGO hit uniquely matched to a "
           << "same-ID Ge hit: " << validAndWellMatchedEvents_ << " / "
           << validFoldEvents_ << " ("
           << percentage(validAndWellMatchedEvents_, validFoldEvents_)
           << "%)\n";

    for (unsigned int germaniumMultiplicity =
             kMinimumReportedGeMultiplicity;
         germaniumMultiplicity <= kMaximumReportedGeMultiplicity;
         ++germaniumMultiplicity) {
        const auto& counts = bgoByGermanium_[germaniumMultiplicity - 1];
        std::uint64_t total = 0;
        for (unsigned int category = 0;
             category <= kOverflowBgoCategory; ++category) {
            total += counts[category];
        }
        const std::uint64_t retained =
            foldValidByGermanium_[germaniumMultiplicity - 1];

        output << "\nGe multiplicity " << germaniumMultiplicity
               << ": " << total << " events; FoldValid " << retained
               << " (" << percentage(retained, total) << "%)\n";
        for (unsigned int bgoMultiplicity = 0;
             bgoMultiplicity <= kMaximumExplicitBgoMultiplicity;
             ++bgoMultiplicity) {
            output << "  BGO multiplicity " << std::setw(2)
                   << bgoMultiplicity << ": "
                   << counts[bgoMultiplicity] << " ("
                   << percentage(counts[bgoMultiplicity], total)
                   << "%)\n";
        }
        output << "  BGO multiplicity 11+: "
               << counts[kOverflowBgoCategory] << " ("
               << percentage(counts[kOverflowBgoCategory], total)
               << "%)\n";
    }
}
