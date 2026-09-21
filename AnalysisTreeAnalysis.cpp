#include "AnalysisTreeAnalysis.h"

#include "AnalysisConfig.h"
#include "ProgressReporter.h"

#include <ROOT/TTreeProcessorMT.hxx>
#include <TChain.h>
#include <TChainElement.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TNamed.h>
#include <TObjArray.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <unordered_map>
#include <vector>

namespace {

const config::DetectorDefinition& definitionFor(unsigned short type)
{
    const auto found = std::find_if(
        config::kDetectors.begin(), config::kDetectors.end(),
        [type](const config::DetectorDefinition& definition) {
            return definition.type == type;
        });
    if (found == config::kDetectors.end()) {
        throw std::logic_error("Requested detector type is not configured");
    }
    return *found;
}

} // namespace

AnalysisTreeAnalysis::AnalysisTreeAnalysis()
    : individualSiliconHistograms_(definitionFor(config::kSiliconType)),
      individualBgoHistograms_(definitionFor(config::kBgoType)),
      individualLabrHistograms_(definitionFor(config::kLabrType)),
      eventMultiplicity_(std::make_unique<TH1D>(
          "h1_EventMultiplicity", "Event multiplicity;All hits in event;Events",
          config::kMultiplicityBins, config::kMultiplicityMin,
          config::kMultiplicityMax)),
      unknownDetectorTypes_(std::make_unique<TH1D>(
          "h1_UnknownDetectorType",
          "Unrecognised detector type;detectorType;Hits",
          65536, -0.5, 65535.5))
{
    detectorHistograms_.reserve(config::kDetectors.size());
    for (const auto& definition : config::kDetectors) {
        detectorIndex_.emplace(definition.type, detectorHistograms_.size());
        detectorHistograms_.emplace_back(definition);
    }
    eventMultiplicity_->SetDirectory(nullptr);
    unknownDetectorTypes_->SetDirectory(nullptr);
    eventMultiplicity_->SetOption("HIST");
    unknownDetectorTypes_->SetOption("HIST");
}

AnalysisTreeAnalysis::~AnalysisTreeAnalysis() = default;

void AnalysisTreeAnalysis::setDiagnosticsEnabled(bool enabled)
{
    diagnosticsEnabled_ = enabled;
}

void AnalysisTreeAnalysis::setProgressEnabled(bool enabled)
{
    progressEnabled_ = enabled;
}

void AnalysisTreeAnalysis::setThreadCount(unsigned int threadCount)
{
    if (threadCount == 0) {
        throw std::runtime_error("thread count must be at least 1");
    }
    threadCount_ = threadCount;
}

void AnalysisTreeAnalysis::configureCalibratedAxes()
{
    detectorHistograms_[detectorIndex_.at(config::kGermaniumType)]
        .setCalibratedEnergyAxis();
    gammaCoincidences_.setCalibratedEnergyAxes();
    angularCoincidences_.setCalibratedEnergyAxes();
    individualGermaniumHistograms_.setCalibratedEnergyAxes();
    germaniumConditionHistograms_.setCalibratedEnergyAxes();
}

void AnalysisTreeAnalysis::processReader(
    TTreeReader& reader, ProgressReporter* progressReporter)
{
    TTreeReaderValue<ULong64_t> absoluteTime(reader, "absoluteTime");
    TTreeReaderValue<UInt_t> eventMultiplicity(reader, "eventMultiplicity");
    TTreeReaderValue<UInt_t> geMultiplicity(reader, "geMultiplicity");
    TTreeReaderValue<UInt_t> geMultiplicityVetoed(
        reader, "geMultiplicityVetoed");
    TTreeReaderValue<UInt_t> bgoMultiplicity(reader, "bgoMultiplicity");
    TTreeReaderValue<UInt_t> siMultiplicity(reader, "siMultiplicity");
    TTreeReaderValue<UInt_t> labrMultiplicity(reader, "labrMultiplicity");
    TTreeReaderValue<UInt_t> extraBgoMultiplicity(
        reader, "extraBgoMultiplicity");
    TTreeReaderValue<Bool_t> storedFoldValid(reader, "foldValid");
    TTreeReaderValue<Bool_t> siliconCoincident(reader, "siliconCoincident");
    TTreeReaderValue<Bool_t> geEnergyCalibrated(
        reader, "geEnergyCalibrated");
    TTreeReaderValue<UInt_t> geTimingCandidates(
        reader, "geTimingCandidates");
    TTreeReaderValue<UInt_t> geTimingPass(reader, "geTimingPass");
    TTreeReaderValue<UInt_t> bgoTimingPass(reader, "bgoTimingPass");
    TTreeReaderValue<UInt_t> siTimingPass(reader, "siTimingPass");
    TTreeReaderValue<UInt_t> siEnergyPass(reader, "siEnergyPass");
    TTreeReaderValue<UInt_t> siCombinedPass(reader, "siCombinedPass");
    TTreeReaderValue<UInt_t> nonzeroPsdHits(reader, "nonzeroPsdHits");
    TTreeReaderValue<std::vector<UShort_t>> geID(reader, "geID");
    TTreeReaderValue<std::vector<Double_t>> geEnergy(reader, "geEnergy");
    TTreeReaderValue<std::vector<Float_t>> geTime(reader, "geTime");
    TTreeReaderValue<std::vector<UChar_t>> geSurvivesBgoVeto(
        reader, "geSurvivesBgoVeto");
    TTreeReaderValue<std::vector<UShort_t>> bgoID(reader, "bgoID");
    TTreeReaderValue<std::vector<Float_t>> bgoEnergy(reader, "bgoEnergy");
    TTreeReaderValue<std::vector<Float_t>> bgoTime(reader, "bgoTime");
    TTreeReaderValue<std::vector<UShort_t>> siID(reader, "siID");
    TTreeReaderValue<std::vector<Float_t>> siEnergy(reader, "siEnergy");
    TTreeReaderValue<std::vector<Float_t>> siTime(reader, "siTime");
    TTreeReaderValue<std::vector<UShort_t>> labrID(reader, "labrID");
    TTreeReaderValue<std::vector<Float_t>> labrEnergy(reader, "labrEnergy");
    TTreeReaderValue<std::vector<Float_t>> labrTime(reader, "labrTime");
    TTreeReaderValue<std::vector<UShort_t>> unknownDetectorType(
        reader, "unknownDetectorType");

    constexpr std::uint64_t kProgressBatchSize = 10000;
    std::uint64_t pendingProgress = 0;
    while (reader.Next()) {
        ++processedEvents_;
        ++pendingProgress;
        if (progressReporter != nullptr &&
            pendingProgress >= kProgressBatchSize) {
            progressReporter->add(pendingProgress);
            pendingProgress = 0;
        }
        const bool sizesAgree =
            geID->size() == geEnergy->size() &&
            geID->size() == geTime->size() &&
            geID->size() == geSurvivesBgoVeto->size() &&
            bgoID->size() == bgoEnergy->size() &&
            bgoID->size() == bgoTime->size() &&
            siID->size() == siEnergy->size() &&
            siID->size() == siTime->size() &&
            labrID->size() == labrEnergy->size() &&
            labrID->size() == labrTime->size();
        const UInt_t countedVetoSurvivors = static_cast<UInt_t>(std::count_if(
            geSurvivesBgoVeto->begin(), geSurvivesBgoVeto->end(),
            [](UChar_t value) { return value != 0; }));
        if (!sizesAgree || *geMultiplicity != geID->size() ||
            *bgoMultiplicity != bgoID->size() ||
            *siMultiplicity != siID->size() ||
            *labrMultiplicity != labrID->size() ||
            *geMultiplicityVetoed != countedVetoSurvivors) {
            ++malformedEvents_;
            continue;
        }

        calibratedEnergySeen_ = calibratedEnergySeen_ || *geEnergyCalibrated;
        uncalibratedEnergySeen_ =
            uncalibratedEnergySeen_ || !*geEnergyCalibrated;
        firstAbsoluteTime_ = std::min(firstAbsoluteTime_, *absoluteTime);
        lastAbsoluteTime_ = std::max(lastAbsoluteTime_, *absoluteTime);
        nonzeroPsdHits_ += *nonzeroPsdHits;
        for (const UShort_t type : *unknownDetectorType) {
            ++unknownHits_;
            unknownDetectorTypes_->Fill(type);
        }

        gateStatistics_.recordGermaniumCounts(
            *geTimingCandidates, *geTimingPass);
        gateStatistics_.recordBgoCounts(*bgoMultiplicity, *bgoTimingPass);
        gateStatistics_.recordSiliconCounts(
            *siMultiplicity, *siTimingPass, *siEnergyPass, *siCombinedPass);
        gateStatistics_.recordEvent(*siliconCoincident);

        const auto matching = detectorMatchingStatistics_.record(*geID, *bgoID);
        const bool foldValid = foldStatistics_.record(
            *geMultiplicity, *bgoMultiplicity,
            *extraBgoMultiplicity == 0);
        if (foldValid != static_cast<bool>(*storedFoldValid) ||
            matching.extraBgoMultiplicity != *extraBgoMultiplicity) {
            ++malformedEvents_;
            continue;
        }

        std::vector<unsigned int> multiplicities(
            detectorHistograms_.size(), 0U);
        multiplicities[detectorIndex_.at(config::kGermaniumType)] =
            *geMultiplicity;
        multiplicities[detectorIndex_.at(config::kBgoType)] =
            *bgoMultiplicity;
        multiplicities[detectorIndex_.at(config::kSiliconType)] =
            *siMultiplicity;
        multiplicities[detectorIndex_.at(config::kLabrType)] =
            *labrMultiplicity;

        std::vector<double> allGammaEnergies;
        std::vector<double> vetoedGammaEnergies;
        std::vector<unsigned short> vetoedGammaIDs;
        allGammaEnergies.reserve(geEnergy->size());
        for (std::size_t hit = 0; hit < geID->size(); ++hit) {
            const double energy = geEnergy->at(hit);
            allGammaEnergies.push_back(energy);
            detectorHistograms_[detectorIndex_.at(config::kGermaniumType)]
                .fillHit(geID->at(hit), energy, geTime->at(hit));
            individualGermaniumHistograms_.fill(geID->at(hit), energy);
            const bool survives = geSurvivesBgoVeto->at(hit) != 0;
            gateStatistics_.recordBgoVetoDecision(survives);
            germaniumConditionHistograms_.fillHit(
                energy, static_cast<double>(*absoluteTime), survives,
                *siliconCoincident, foldValid);
            if (survives) {
                vetoedGammaEnergies.push_back(energy);
                vetoedGammaIDs.push_back(geID->at(hit));
            }
        }
        germaniumConditionHistograms_.fillUnconditionedCoincidences(
            allGammaEnergies);
        germaniumConditionHistograms_.fillEvent(
            *geMultiplicityVetoed, *siliconCoincident, foldValid);
        gammaCoincidences_.fillEvent(
            vetoedGammaEnergies, *siliconCoincident);
        if (*siliconCoincident) {
            germaniumConditionHistograms_.fillBgoVetoedSiliconCoincidences(
                vetoedGammaEnergies);
            angularCoincidences_.fillEvent(
                vetoedGammaEnergies, vetoedGammaIDs);
        }

        const auto fillDetector = [&](unsigned short type,
                                      const auto& ids,
                                      const auto& energies,
                                      const auto& times,
                                      IndividualDetectorHistograms& individual) {
            DetectorHistograms& histograms =
                detectorHistograms_[detectorIndex_.at(type)];
            for (std::size_t hit = 0; hit < ids.size(); ++hit) {
                histograms.fillHit(ids[hit], energies[hit], times[hit]);
                individual.fill(ids[hit], energies[hit]);
            }
        };
        fillDetector(config::kBgoType, *bgoID, *bgoEnergy, *bgoTime,
                     individualBgoHistograms_);
        fillDetector(config::kSiliconType, *siID, *siEnergy, *siTime,
                     individualSiliconHistograms_);
        fillDetector(config::kLabrType, *labrID, *labrEnergy, *labrTime,
                     individualLabrHistograms_);

        for (std::size_t index = 0; index < detectorHistograms_.size(); ++index) {
            detectorHistograms_[index].fillMultiplicity(multiplicities[index]);
        }
        eventMultiplicity_->Fill(*eventMultiplicity);
        combinedStatistics_.record(multiplicities);

    }
    if (progressReporter != nullptr) progressReporter->add(pendingProgress);
}

void AnalysisTreeAnalysis::merge(const AnalysisTreeAnalysis& other)
{
    for (std::size_t index = 0; index < detectorHistograms_.size(); ++index) {
        detectorHistograms_[index].merge(other.detectorHistograms_[index]);
    }
    gammaCoincidences_.merge(other.gammaCoincidences_);
    angularCoincidences_.merge(other.angularCoincidences_);
    individualGermaniumHistograms_.merge(other.individualGermaniumHistograms_);
    germaniumConditionHistograms_.merge(other.germaniumConditionHistograms_);
    individualSiliconHistograms_.merge(other.individualSiliconHistograms_);
    individualBgoHistograms_.merge(other.individualBgoHistograms_);
    individualLabrHistograms_.merge(other.individualLabrHistograms_);
    eventMultiplicity_->Add(other.eventMultiplicity_.get());
    unknownDetectorTypes_->Add(other.unknownDetectorTypes_.get());
    combinedStatistics_.merge(other.combinedStatistics_);
    foldStatistics_.merge(other.foldStatistics_);
    gateStatistics_.merge(other.gateStatistics_);
    detectorMatchingStatistics_.merge(other.detectorMatchingStatistics_);
    calibratedEnergySeen_ = calibratedEnergySeen_ || other.calibratedEnergySeen_;
    uncalibratedEnergySeen_ =
        uncalibratedEnergySeen_ || other.uncalibratedEnergySeen_;
    processedEvents_ += other.processedEvents_;
    malformedEvents_ += other.malformedEvents_;
    unknownHits_ += other.unknownHits_;
    nonzeroPsdHits_ += other.nonzeroPsdHits_;
    firstAbsoluteTime_ = std::min(firstAbsoluteTime_, other.firstAbsoluteTime_);
    lastAbsoluteTime_ = std::max(lastAbsoluteTime_, other.lastAbsoluteTime_);
}

int AnalysisTreeAnalysis::writeOutput(const std::string& outputFileName)
{
    if (calibratedEnergySeen_) {
        configureCalibratedAxes();
    }
    TFile outputFile(outputFileName.c_str(), "RECREATE");
    if (outputFile.IsZombie()) {
        std::cerr << "Error: could not create '" << outputFileName << "'.\n";
        return 3;
    }
    outputFile.cd();
    eventMultiplicity_->Write();
    unknownDetectorTypes_->Write();
    std::unordered_map<unsigned short, TDirectory*> directories;
    for (const auto& histograms : detectorHistograms_) {
        directories.emplace(histograms.detectorType(), histograms.write(outputFile));
    }
    TDirectory* ge = directories.at(config::kGermaniumType);
    TDirectory* si = directories.at(config::kSiliconType);
    TDirectory* bgo = directories.at(config::kBgoType);
    TDirectory* labr = directories.at(config::kLabrType);
    TDirectory* geEnergy = ge->GetDirectory("Energy");
    TDirectory* geTime = ge->GetDirectory("Time");
    TDirectory* siEnergy = si->GetDirectory("Energy");
    TDirectory* bgoEnergy = bgo->GetDirectory("Energy");
    TDirectory* labrEnergy = labr->GetDirectory("Energy");
    if (!geEnergy || !geTime || !siEnergy || !bgoEnergy || !labrEnergy) {
        return 4;
    }
    germaniumConditionHistograms_.write(*ge, *geEnergy, *geTime);
    individualGermaniumHistograms_.write(*geEnergy);
    individualSiliconHistograms_.write(*siEnergy);
    individualBgoHistograms_.write(*bgoEnergy);
    individualLabrHistograms_.write(*labrEnergy);
    gammaCoincidences_.write(outputFile);
    TDirectory* coincidences = outputFile.GetDirectory("Coincidences");
    if (!coincidences) return 4;
    angularCoincidences_.write(*coincidences);
    outputFile.Close();
    return 0;
}

void AnalysisTreeAnalysis::printDiagnostics() const
{
    std::cout << "Malformed analysis events skipped: " << malformedEvents_ << "\n"
              << "Hits with unknown detectorType: " << unknownHits_ << "\n"
              << "Hits with non-zero psd: " << nonzeroPsdHits_ << "\n";
    if (processedEvents_ > malformedEvents_) {
        std::cout << "Raw absoluteTime range: " << firstAbsoluteTime_
                  << " to " << lastAbsoluteTime_ << "\n";
    }
    combinedStatistics_.print(std::cout, "combined analysis-tree events");
    foldStatistics_.print(std::cout);
    detectorMatchingStatistics_.print(std::cout);
    gateStatistics_.print(std::cout);
}

int AnalysisTreeAnalysis::run(
    const std::vector<std::string>& inputPatterns,
    const std::string& outputFileName)
{
    const auto totalStart = std::chrono::steady_clock::now();
    TChain chain(config::kAnalysisTreeName);
    for (const std::string& pattern : inputPatterns) {
        if (chain.Add(pattern.c_str()) == 0) {
            std::cerr << "Warning: no file matched '" << pattern << "'.\n";
        }
    }
    if (chain.GetNtrees() == 0) {
        std::cerr << "Error: no analysis trees were added.\n";
        return 2;
    }
    std::string expectedConfiguration;
    TObjArray* fileElements = chain.GetListOfFiles();
    for (int index = 0; index < fileElements->GetEntries(); ++index) {
        const auto* element =
            dynamic_cast<const TChainElement*>(fileElements->At(index));
        if (element == nullptr) continue;
        TFile inputFile(element->GetTitle(), "READ");
        TNamed* configuration = nullptr;
        inputFile.GetObject("AnalysisConfiguration", configuration);
        if (configuration == nullptr) {
            std::cerr << "Error: '" << element->GetTitle()
                      << "' has no AnalysisConfiguration metadata.\n";
            return 2;
        }
        const std::string value = configuration->GetTitle();
        if (expectedConfiguration.empty()) {
            expectedConfiguration = value;
        } else if (value != expectedConfiguration) {
            std::cerr << "Error: analysis-tree files were built with "
                      << "different gate/exclusion configurations.\n";
            return 2;
        }
    }
    std::cout << "Reading " << chain.GetNtrees() << " analysis file(s), "
              << chain.GetEntries() << " events.\n";

    const auto processingStart = std::chrono::steady_clock::now();
    ProgressReporter progressReporter(
        static_cast<std::uint64_t>(chain.GetEntries()), progressEnabled_);
    try {
        if (threadCount_ == 1) {
            TTreeReader reader(&chain);
            processReader(reader, &progressReporter);
        } else {
            std::cout << "Processing with " << threadCount_
                      << " worker threads.\n";
            ROOT::TTreeProcessorMT processor(chain, threadCount_);
            std::mutex mutex;
            std::unordered_map<std::thread::id, AnalysisTreeAnalysis*> byThread;
            std::vector<std::unique_ptr<AnalysisTreeAnalysis>> workers;
            processor.Process([&](TTreeReader& reader) {
                AnalysisTreeAnalysis* worker = nullptr;
                {
                    std::lock_guard<std::mutex> lock(mutex);
                    const auto id = std::this_thread::get_id();
                    const auto found = byThread.find(id);
                    if (found != byThread.end()) {
                        worker = found->second;
                    } else {
                        auto state = std::make_unique<AnalysisTreeAnalysis>();
                        state->diagnosticsEnabled_ = false;
                        worker = state.get();
                        workers.push_back(std::move(state));
                        byThread.emplace(id, worker);
                    }
                }
                worker->processReader(reader, &progressReporter);
            });
            for (const auto& worker : workers) merge(*worker);
        }
    } catch (const std::exception& error) {
        std::cerr << "Analysis-tree processing error: " << error.what() << "\n";
        return 5;
    }
    progressReporter.finish();
    const double processingSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - processingStart).count();
    if (calibratedEnergySeen_ && uncalibratedEnergySeen_) {
        std::cerr << "Error: calibrated and uncalibrated analysis trees "
                  << "cannot be combined.\n";
        return 5;
    }
    std::cout << "Processed " << processedEvents_ << " events.\n";
    const int writeStatus = writeOutput(outputFileName);
    if (writeStatus != 0) return writeStatus;
    std::cout << "Wrote histograms to " << outputFileName << ".\n";
    if (diagnosticsEnabled_) printDiagnostics();
    const double totalSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - totalStart).count();
    std::cout << "\n=== Analysis timing ===\n"
              << std::fixed << std::setprecision(3)
              << "Event processing: " << processingSeconds << " s ("
              << std::setprecision(0)
              << (processingSeconds > 0.0 ? processedEvents_ / processingSeconds : 0.0)
              << " events/s)\n"
              << std::setprecision(3) << "Total analysis:   " << totalSeconds
              << " s (" << std::setprecision(0)
              << (totalSeconds > 0.0 ? processedEvents_ / totalSeconds : 0.0)
              << " events/s)\n";
    return 0;
}
