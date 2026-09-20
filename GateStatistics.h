#ifndef GATE_STATISTICS_H
#define GATE_STATISTICS_H

#include <cstdint>
#include <iosfwd>

class GateStatistics {
public:
    void recordGermanium(bool passesTiming);
    void recordBgo(bool passesTiming);
    void recordSilicon(bool passesTiming, bool passesEnergy);
    void recordEvent(bool hasSiliconCoincidence);
    void recordBgoVetoDecision(bool survivesVeto);
    void merge(const GateStatistics& other);

    void print(std::ostream& output) const;

private:
    std::uint64_t germaniumHits_ = 0;
    std::uint64_t germaniumTimingPass_ = 0;
    std::uint64_t bgoHits_ = 0;
    std::uint64_t bgoTimingPass_ = 0;
    std::uint64_t siliconHits_ = 0;
    std::uint64_t siliconTimingPass_ = 0;
    std::uint64_t siliconEnergyPass_ = 0;
    std::uint64_t siliconCombinedPass_ = 0;
    std::uint64_t events_ = 0;
    std::uint64_t siliconCoincidenceEvents_ = 0;
    std::uint64_t germaniumVetoCandidates_ = 0;
    std::uint64_t germaniumVetoSurvivors_ = 0;
};

#endif
