#ifndef DETECTOR_HISTOGRAMS_H
#define DETECTOR_HISTOGRAMS_H

#include "AnalysisConfig.h"

#include <memory>

class TDirectory;
class TH1D;

class DetectorHistograms {
public:
    explicit DetectorHistograms(const config::DetectorDefinition& definition);
    ~DetectorHistograms();

    DetectorHistograms(const DetectorHistograms&) = delete;
    DetectorHistograms& operator=(const DetectorHistograms&) = delete;
    DetectorHistograms(DetectorHistograms&&) noexcept;
    DetectorHistograms& operator=(DetectorHistograms&&) noexcept;

    unsigned short detectorType() const;
    const char* name() const;
    const char* directoryName() const;

    void fillHit(unsigned short detectorID,
                 double energy,
                 double relativeTimeNs);
    void fillMultiplicity(unsigned int multiplicity);
    void merge(const DetectorHistograms& other);
    void setCalibratedEnergyAxis();
    TDirectory* write(TDirectory& parentDirectory) const;

private:
    const config::DetectorDefinition* definition_;
    std::unique_ptr<TH1D> energy_;
    std::unique_ptr<TH1D> time_;
    std::unique_ptr<TH1D> multiplicity_;
    std::unique_ptr<TH1D> idOccupancy_;
};

#endif
