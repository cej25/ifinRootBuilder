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
    void fill(UShort_t detectorID, double energy);
    void merge(const GermaniumDetectorHistograms& other);
    void write(TDirectory& germaniumEnergyDirectory) const;

private:
    TH1D& spectrumFor(UShort_t detectorID);

    bool calibrated_ = false;
    std::map<UShort_t, std::unique_ptr<TH1D>> spectra_;
};

#endif
