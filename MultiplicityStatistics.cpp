#include "MultiplicityStatistics.h"

#include <algorithm>
#include <iomanip>
#include <numeric>
#include <ostream>
#include <sstream>
#include <stdexcept>

namespace {

std::string formatHitCount(std::uint64_t count, std::uint64_t totalHits)
{
    const double percentage = totalHits > 0
        ? 100.0 * static_cast<double>(count) /
              static_cast<double>(totalHits)
        : 0.0;

    std::ostringstream text;
    text << count << " (" << std::fixed << std::setprecision(1)
         << percentage << "%)";
    return text.str();
}

} // namespace

MultiplicityStatistics::MultiplicityStatistics()
{
    for (auto& counts : detectorFoldCounts_) {
        counts.resize(1, 0);
    }
}

std::size_t MultiplicityStatistics::indexForType(
    unsigned short detectorType)
{
    for (std::size_t index = 0; index < config::kDetectors.size(); ++index) {
        if (config::kDetectors[index].type == detectorType) {
            return index;
        }
    }
    throw std::logic_error("Requested detector type is not configured");
}

void MultiplicityStatistics::record(
    const std::vector<unsigned int>& detectorMultiplicities)
{
    if (detectorMultiplicities.size() != config::kDetectors.size()) {
        throw std::logic_error(
            "MultiplicityStatistics received an unexpected detector count");
    }

    unsigned int eventFold = 0;
    for (std::size_t index = 0;
         index < detectorMultiplicities.size(); ++index) {
        const unsigned int fold = detectorMultiplicities[index];
        auto& counts = detectorFoldCounts_[index];
        if (counts.size() <= fold) {
            counts.resize(static_cast<std::size_t>(fold) + 1, 0);
        }
        ++counts[fold];
        eventFold += fold;
    }

    if (eventFold >= 1 && eventFold <= kMaximumReportedEventFold) {
        ++eventCounts_[eventFold];
        for (std::size_t index = 0;
             index < detectorMultiplicities.size(); ++index) {
            detectorHitTotals_[eventFold][index] +=
                detectorMultiplicities[index];
        }
    }
}

void MultiplicityStatistics::printDetectorBreakdown(
    std::ostream& output,
    unsigned short detectorType) const
{
    const std::size_t index = indexForType(detectorType);
    output << config::kDetectors[index].title << " multiplicity:\n";
    const auto& counts = detectorFoldCounts_[index];
    for (std::size_t fold = 0; fold < counts.size(); ++fold) {
        output << "  multiplicity " << fold << ": " << counts[fold]
               << " events\n";
    }
}

void MultiplicityStatistics::print(std::ostream& output,
                                   const std::string& label) const
{
    output << "\n=== Multiplicity statistics: " << label << " ===\n";
    printDetectorBreakdown(output, config::kGermaniumType);
    printDetectorBreakdown(output, config::kBgoType);
    printDetectorBreakdown(output, config::kSiliconType);

    const std::array<unsigned short, 4> displayOrder{{
        config::kGermaniumType,
        config::kBgoType,
        config::kSiliconType,
        config::kLabrType
    }};

    output << "Event-multiplicity composition (summed detector hits):\n"
           << std::setw(14) << "multiplicity"
           << std::setw(14) << "events"
           << std::setw(24) << "Ge hits (%)"
           << std::setw(24) << "BGO hits (%)"
           << std::setw(24) << "Si hits (%)"
           << std::setw(24) << "LaBr3 hits (%)" << '\n';

    for (unsigned int fold = 1;
         fold <= kMaximumReportedEventFold; ++fold) {
        output << std::setw(14) << fold
               << std::setw(14) << eventCounts_[fold];
        const std::uint64_t totalHits =
            static_cast<std::uint64_t>(fold) * eventCounts_[fold];
        for (const unsigned short type : displayOrder) {
            const std::uint64_t count =
                detectorHitTotals_[fold][indexForType(type)];
            output << std::setw(24) << formatHitCount(count, totalHits);
        }
        output << '\n';
    }
}
