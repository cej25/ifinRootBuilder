#include "GateStatistics.h"

#include "AnalysisConfig.h"

#include <iomanip>
#include <ostream>
#include <sstream>

namespace {

double percentage(std::uint64_t count, std::uint64_t total)
{
    return total == 0 ? 0.0
                      : 100.0 * static_cast<double>(count) /
                            static_cast<double>(total);
}

void printPass(std::ostream& output, const char* label,
               std::uint64_t count, std::uint64_t total)
{
    output << "  " << label << ": " << count << " / " << total
           << " (" << std::fixed << std::setprecision(2)
           << percentage(count, total) << "%)\n";
}

} // namespace

void GateStatistics::recordGermanium(bool passesTiming)
{
    ++germaniumHits_;
    germaniumTimingPass_ += passesTiming;
}

void GateStatistics::recordBgo(bool passesTiming)
{
    ++bgoHits_;
    bgoTimingPass_ += passesTiming;
}

void GateStatistics::recordSilicon(bool passesTiming, bool passesEnergy)
{
    ++siliconHits_;
    siliconTimingPass_ += passesTiming;
    siliconEnergyPass_ += passesEnergy;
    siliconCombinedPass_ += passesTiming && passesEnergy;
}

void GateStatistics::recordEvent(bool hasSiliconCoincidence)
{
    ++events_;
    siliconCoincidenceEvents_ += hasSiliconCoincidence;
}

void GateStatistics::recordBgoVetoDecision(bool survivesVeto)
{
    ++germaniumVetoCandidates_;
    germaniumVetoSurvivors_ += survivesVeto;
}

void GateStatistics::recordGermaniumCounts(
    std::uint64_t total, std::uint64_t timingPass)
{
    germaniumHits_ += total;
    germaniumTimingPass_ += timingPass;
}

void GateStatistics::recordBgoCounts(
    std::uint64_t total, std::uint64_t timingPass)
{
    bgoHits_ += total;
    bgoTimingPass_ += timingPass;
}

void GateStatistics::recordSiliconCounts(
    std::uint64_t total, std::uint64_t timingPass,
    std::uint64_t energyPass, std::uint64_t combinedPass)
{
    siliconHits_ += total;
    siliconTimingPass_ += timingPass;
    siliconEnergyPass_ += energyPass;
    siliconCombinedPass_ += combinedPass;
}

void GateStatistics::merge(const GateStatistics& other)
{
    germaniumHits_ += other.germaniumHits_;
    germaniumTimingPass_ += other.germaniumTimingPass_;
    bgoHits_ += other.bgoHits_;
    bgoTimingPass_ += other.bgoTimingPass_;
    siliconHits_ += other.siliconHits_;
    siliconTimingPass_ += other.siliconTimingPass_;
    siliconEnergyPass_ += other.siliconEnergyPass_;
    siliconCombinedPass_ += other.siliconCombinedPass_;
    events_ += other.events_;
    siliconCoincidenceEvents_ += other.siliconCoincidenceEvents_;
    germaniumVetoCandidates_ += other.germaniumVetoCandidates_;
    germaniumVetoSurvivors_ += other.germaniumVetoSurvivors_;
}

void GateStatistics::print(std::ostream& output) const
{
    const auto rangeLabel = [](const char* prefix, double minimum,
                               double maximum, const char* units) {
        std::ostringstream label;
        label << prefix << ' ' << minimum << "--" << maximum << ' ' << units;
        return label.str();
    };

    output << "\n=== Timing, energy, and veto statistics ===\n"
           << "Raw well-formed hits (before calibration-event rejection):\n";
    printPass(output, rangeLabel("Ge timing", config::kGermaniumTimeMinNs,
                                config::kGermaniumTimeMaxNs, "ns").c_str(),
              germaniumTimingPass_,
              germaniumHits_);
    printPass(output, rangeLabel("BGO timing", config::kBgoVetoTimeMinNs,
                                config::kBgoVetoTimeMaxNs, "ns").c_str(),
              bgoTimingPass_, bgoHits_);
    printPass(output, rangeLabel("Si timing", config::kSiliconTimeMinNs,
                                config::kSiliconTimeMaxNs, "ns").c_str(),
              siliconTimingPass_,
              siliconHits_);
    printPass(output, rangeLabel("Si energy", config::kSiliconEnergyMin,
                                config::kSiliconEnergyMax,
                                "raw units").c_str(), siliconEnergyPass_,
              siliconHits_);
    printPass(output, "Si timing AND energy", siliconCombinedPass_,
              siliconHits_);
    printPass(output, "Events with a qualifying Si hit",
              siliconCoincidenceEvents_, events_);

    output << "Accepted, time-gated Ge hits (after calibration-event rejection):\n";
    printPass(output, "Surviving matching-ID BGO veto",
              germaniumVetoSurvivors_, germaniumVetoCandidates_);
    printPass(output, "Rejected by matching-ID BGO veto",
              germaniumVetoCandidates_ - germaniumVetoSurvivors_,
              germaniumVetoCandidates_);
}
