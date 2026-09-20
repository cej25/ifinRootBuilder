#ifndef INDIVIDUAL_DETECTOR_HISTOGRAMS_H
#define INDIVIDUAL_DETECTOR_HISTOGRAMS_H

#include "AnalysisConfig.h"

#include <map>
#include <memory>

class TDirectory;
class TH1D;

class IndividualDetectorHistograms {
public:
    explicit IndividualDetectorHistograms(
        const config::DetectorDefinition& definition);
    ~IndividualDetectorHistograms();

    IndividualDetectorHistograms(const IndividualDetectorHistograms&) = delete;
    IndividualDetectorHistograms& operator=(
        const IndividualDetectorHistograms&) = delete;

    void fill(unsigned short detectorLUT, double energy);
    void write(TDirectory& energyDirectory) const;

private:
    TH1D& spectrumFor(unsigned short detectorLUT);

    const config::DetectorDefinition* definition_;
    std::map<unsigned short, std::unique_ptr<TH1D>> spectra_;
};

#endif
