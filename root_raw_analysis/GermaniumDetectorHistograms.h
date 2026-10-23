#ifndef GERMANIUM_DETECTOR_HISTOGRAMS_H
#define GERMANIUM_DETECTOR_HISTOGRAMS_H

#include <RtypesCore.h>

#include <map>
#include <memory>

class TDirectory;
class TH1D;
class TH2I;

class GermaniumDetectorHistograms {
public:
    GermaniumDetectorHistograms();
    ~GermaniumDetectorHistograms();

    GermaniumDetectorHistograms(const GermaniumDetectorHistograms&) = delete;
    GermaniumDetectorHistograms& operator=(
        const GermaniumDetectorHistograms&) = delete;

    void setCalibratedEnergyAxes();
    void setRunningTimeRange(double totalRunningTimeSeconds);
    void fill(UShort_t detectorID, double calibratedEnergy,
              double rawEnergy, double runningTimeSeconds,
              bool runningTimeValid);
    void merge(const GermaniumDetectorHistograms& other);
    void write(TDirectory& germaniumEnergyDirectory,
               TDirectory& germaniumTimeDirectory) const;

private:
    TH1D& calibratedSpectrumFor(UShort_t detectorID);
    TH1D& rawSpectrumFor(UShort_t detectorID);
    TH2I& timeSpectrumFor(UShort_t detectorID);

    bool calibrated_ = false;
    double runningTimeMaxSeconds_ = 1.0;
    std::map<UShort_t, std::unique_ptr<TH1D>> calibratedSpectra_;
    std::map<UShort_t, std::unique_ptr<TH1D>> rawSpectra_;
    std::map<UShort_t, std::unique_ptr<TH2I>> timeSpectra_;
};

#endif
