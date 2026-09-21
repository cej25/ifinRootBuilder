#ifndef ANALYSIS_CONFIG_H
#define ANALYSIS_CONFIG_H

#include <array>

namespace config {

inline constexpr const char* kTreeName = "events";
inline constexpr const char* kAnalysisTreeName = "analysis";
inline constexpr unsigned int kAnalysisTreeSchemaVersion = 1;

// The stored branch retains its historical name; analysis terminology uses ID.
inline constexpr const char* kDetectorIDBranch   = "detectorLUT";
inline constexpr const char* kDetectorTypeBranch  = "detectorType";
inline constexpr const char* kEnergyBranch        = "energy";
inline constexpr const char* kPsdBranch           = "psd";
inline constexpr const char* kRelativeNsBranch    = "relativeNsTime";
inline constexpr const char* kRelativePsBranch    = "relativePsTime";
inline constexpr const char* kAbsoluteTimeBranch  = "absoluteTime";

// relativeNsTime is unsigned, so the provisional event-time range is 0--2.2 us.
// Time is displayed in ns after reconstructing it from the coarse and fine parts.
inline constexpr int    kTimeBins  = 2200;
inline constexpr double kTimeMinNs = 0.0;
inline constexpr double kTimeMaxNs = 2200.0;

inline constexpr int    kMultiplicityBins = 101;
inline constexpr double kMultiplicityMin  = -0.5;
inline constexpr double kMultiplicityMax  = 100.5;

// All detector types currently have at most 25 physical IDs: 0--24.
inline constexpr int    kIdBins = 25;
inline constexpr double kIdMin  = -0.5;
inline constexpr double kIdMax  = 24.5;

inline constexpr unsigned short kGermaniumType = 0;
inline constexpr unsigned short kSiliconType   = 1;
inline constexpr unsigned short kBgoType       = 2;
inline constexpr unsigned short kLabrType      = 3;

// Inclusive hit-selection gates. Relative-time values are reconstructed in ns.
inline constexpr double kGermaniumTimeMinNs = 985.0;
inline constexpr double kGermaniumTimeMaxNs = 1100.0;
inline constexpr double kBgoVetoTimeMinNs    = 850.0;
inline constexpr double kBgoVetoTimeMaxNs    = 1150.0;
inline constexpr double kSiliconTimeMinNs    = 940.0;
inline constexpr double kSiliconTimeMaxNs    = 1050.0;
inline constexpr double kSiliconEnergyMin    = 1000.0;
inline constexpr double kSiliconEnergyMax    = 4000.0;

inline constexpr int    kGermaniumSpectrumBins = 4096;
inline constexpr double kGermaniumSpectrumMin  = 0.0;
inline constexpr double kGermaniumSpectrumMax  = 4096.0;
inline constexpr int    kGermaniumIdBins      = 25;

inline constexpr int    kGammaBins = 2048;
inline constexpr double kGammaMin  = 0.0;
inline constexpr double kGammaMax  = 2048.0;
// Half-open energy windows corresponding to the requested inclusive 1-keV
// channels: prompt 292--299, lower 283--289, and upper 310--317.
inline constexpr double kGammaGateMin = 292.0;
inline constexpr double kGammaGateMaxExclusive = 300.0;
inline constexpr double kGammaLowerSidebandMin = 283.0;
inline constexpr double kGammaLowerSidebandMaxExclusive = 290.0;
inline constexpr double kGammaUpperSidebandMin = 310.0;
inline constexpr double kGammaUpperSidebandMaxExclusive = 318.0;
inline constexpr double kGammaSidebandScale =
    (kGammaGateMaxExclusive - kGammaGateMin) /
    ((kGammaLowerSidebandMaxExclusive - kGammaLowerSidebandMin) +
     (kGammaUpperSidebandMaxExclusive - kGammaUpperSidebandMin));

// absoluteTime is currently interpreted as seconds. Values outside this range
// remain visible in the underflow/overflow bins and can be adjusted here.
inline constexpr int    kDriftTimeBins = 600;
inline constexpr double kDriftTimeMinSeconds = 0.0;
inline constexpr double kDriftTimeMaxSeconds = 600.0;
inline constexpr int    kDriftEnergyBins = 2048;
inline constexpr double kDriftEnergyMin = 0.0;
inline constexpr double kDriftEnergyMax = 2048.0;

// Passive configuration record; analysis state and behaviour live in classes.
struct DetectorDefinition {
    unsigned short type;
    const char* name;
    const char* directory;
    const char* title;
    const char* histogramTag;
    int energyBins;
    double energyMin;
    double energyMax;
};

// Provisional mapping supplied from inspection of one tree.
inline constexpr std::array<DetectorDefinition, 4> kDetectors{{
    {kGermaniumType, "germanium", "Germanium", "Germanium", "Ge",
     kGermaniumSpectrumBins, kGermaniumSpectrumMin, kGermaniumSpectrumMax},
    {kSiliconType, "silicon", "Silicon", "Silicon", "Si",
     20000, 0.0, 20000.0},
    {kBgoType,     "bgo",     "BGO",     "BGO", "BGO",
     4000, 0.0, 4000.0},
    {kLabrType,    "labr",    "LaBr",    "LaBr3", "LaBr",
     4000, 0.0, 4000.0}
}};

// Germanium angular groups used for RDDS correlations.
struct GermaniumAngleGroup {
    const char* name;
    const char* title;
    unsigned short firstID;
    unsigned short lastID;
};

inline constexpr std::array<GermaniumAngleGroup, 5> kGermaniumAngleGroups{{
    {"forward",         "Forward (37 deg)",        0,  4},
    {"middle_forward",  "Middle-forward (70 deg)", 5,  9},
    {"middle",          "Middle (90 deg)",         10, 14},
    {"middle_backward", "Middle-backward (110 deg)",15, 19},
    {"backward",        "Backward (143 deg)",      20, 24}
}};

inline constexpr unsigned short kForwardIDMin = 0;
inline constexpr unsigned short kForwardIDMax = 4;
inline constexpr unsigned short kBackwardIDMin = 20;
inline constexpr unsigned short kBackwardIDMax = 24;

// Optional RDDS gates on the all-detector (x) axis. Change the array size and
// add entries here later; forward/backward projected spectra are then created
// automatically. Limits are half-open [minimum, maximumExclusive).
struct RddsGateDefinition {
    const char* name;
    double minimum;
    double maximumExclusive;
};

inline constexpr std::array<RddsGateDefinition, 0> kRddsGates{};

} // namespace config

#endif
