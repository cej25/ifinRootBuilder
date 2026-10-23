#ifndef CALIBRATION_SPECTRUM_BUILDER_H
#define CALIBRATION_SPECTRUM_BUILDER_H

#include <string>
#include <vector>

class CalibrationSpectrumBuilder {
public:
    void addGroup(const std::string& name, const std::string& inputPattern);
    void setGermaniumTimeGateEnabled(bool enabled);
    int run(const std::string& outputFileName) const;

private:
    struct Group {
        std::string name;
        std::vector<std::string> inputPatterns;
    };

    Group& findOrAddGroup(const std::string& name);

    bool germaniumTimeGateEnabled_ = true;
    std::vector<Group> groups_;
};

#endif
