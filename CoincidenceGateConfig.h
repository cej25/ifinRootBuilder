#ifndef COINCIDENCE_GATE_CONFIG_H
#define COINCIDENCE_GATE_CONFIG_H

#include <string>
#include <vector>

struct CoincidenceGateDefinition {
    std::string name;
    double promptMinimum = 0.0;
    double promptMaximumExclusive = 0.0;
    double lowerMinimum = 0.0;
    double lowerMaximumExclusive = 0.0;
    double upperMinimum = 0.0;
    double upperMaximumExclusive = 0.0;

    double sidebandScale() const;
    bool operator==(const CoincidenceGateDefinition& other) const;
};

class CoincidenceGateConfig {
public:
    static CoincidenceGateConfig load(const std::string& fileName);

    const std::vector<CoincidenceGateDefinition>& symmetric() const;
    const std::vector<CoincidenceGateDefinition>& allVsForward() const;
    const std::vector<CoincidenceGateDefinition>& allVsBackward() const;

private:
    std::vector<CoincidenceGateDefinition> symmetric_;
    std::vector<CoincidenceGateDefinition> allVsForward_;
    std::vector<CoincidenceGateDefinition> allVsBackward_;
};

#endif
