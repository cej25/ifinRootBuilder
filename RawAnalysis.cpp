#include "RawAnalysis.h"

#include "AnalysisConfig.h"
#include "ProgressReporter.h"

#include <ROOT/TTreeProcessorMT.hxx>
#include <RtypesCore.h>
#include <TChain.h>
#include <TChainElement.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TObjArray.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <sstream>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {

bool insideInclusive(double value, double minimum, double maximum)
{
    return value >= minimum && value <= maximum;
}

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

std::pair<std::uint64_t, std::uint64_t> rawAbsoluteTimeRange(
    const std::string& fileName)
{
    TFile inputFile(fileName.c_str(), "READ");
    if (inputFile.IsZombie()) {
        throw std::runtime_error("could not open input file '" + fileName +
                                 "' while reading absoluteTime");
    }
    TTree* tree = nullptr;
    inputFile.GetObject(config::kTreeName, tree);
    if (tree == nullptr) {
        throw std::runtime_error("input file '" + fileName +
                                 "' has no tree called '" +
                                 config::kTreeName + "'");
    }
    const Long64_t entries = tree->GetEntries();
    if (entries == 0) return {0, 0};

    TTreeReader firstReader(tree);
    TTreeReaderValue<ULong_t> firstValue(
        firstReader, config::kAbsoluteTimeBranch);
    if (!firstReader.Next()) {
        throw std::runtime_error("could not read first absoluteTime from '" +
                                 fileName + "'");
    }
    const std::uint64_t first = *firstValue;

    TTreeReader lastReader(tree);
    TTreeReaderValue<ULong_t> lastValue(
        lastReader, config::kAbsoluteTimeBranch);
    if (lastReader.SetEntry(entries - 1) != TTreeReader::kEntryValid) {
        throw std::runtime_error("could not read last absoluteTime from '" +
                                 fileName + "'");
    }
    return {first, static_cast<std::uint64_t>(*lastValue)};
}

std::string formatElapsed(double seconds)
{
    const auto totalMilliseconds = static_cast<unsigned long long>(
        std::llround(seconds * 1000.0));
    const auto hours = totalMilliseconds / 3600000ULL;
    const auto minutes = (totalMilliseconds / 60000ULL) % 60ULL;
    const auto wholeSeconds = (totalMilliseconds / 1000ULL) % 60ULL;
    const auto milliseconds = totalMilliseconds % 1000ULL;

    std::ostringstream text;
    text << std::setfill('0') << std::setw(2) << hours << ':'
         << std::setw(2) << minutes << ':'
         << std::setw(2) << wholeSeconds << '.'
         << std::setw(3) << milliseconds;
    return text.str();
}

} // namespace

RawAnalysis::RawAnalysis()
    : individualSiliconHistograms_(definitionFor(config::kSiliconType)),
      individualBgoHistograms_(definitionFor(config::kBgoType)),
      individualLabrHistograms_(definitionFor(config::kLabrType)),
      eventMultiplicity_(std::make_unique<TH1D>(
          "h1_EventMultiplicity",
          "Event multiplicity;All hits in event;Events",
          config::kMultiplicityBins,
          config::kMultiplicityMin,
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

RawAnalysis::~RawAnalysis() = default;

void RawAnalysis::addCalFile(const std::string& fileName)
{
    calibrationManager_.addGlobalCalFile(fileName);
}

void RawAnalysis::addMcalFile(const std::string& fileName)
{
    calibrationManager_.addGlobalMcalFile(fileName);
}

void RawAnalysis::addRunCalFile(unsigned int firstRun, unsigned int lastRun,
                                const std::string& fileName)
{
    calibrationManager_.addRunCalFile(firstRun, lastRun, fileName);
}

void RawAnalysis::addRunMcalFile(unsigned int firstRun, unsigned int lastRun,
                                 const std::string& fileName)
{
    calibrationManager_.addRunMcalFile(firstRun, lastRun, fileName);
}

void RawAnalysis::setDiagnosticsEnabled(bool enabled)
{
    diagnosticsEnabled_ = enabled;
}

void RawAnalysis::setProgressEnabled(bool enabled)
{
    progressEnabled_ = enabled;
}

void RawAnalysis::setThreadCount(unsigned int threadCount)
{
    if (threadCount == 0) {
        throw std::runtime_error("thread count must be at least 1");
    }
    threadCount_ = threadCount;
}

void RawAnalysis::excludeGermaniumID(unsigned int detectorID)
{
    if (detectorID >= config::kGermaniumIdBins) {
        throw std::runtime_error(
            "germanium detector ID must be between 0 and " +
            std::to_string(config::kGermaniumIdBins - 1));
    }
    excludedGermaniumIDs_.insert(static_cast<unsigned short>(detectorID));
}

void RawAnalysis::configureCalibratedAxes()
{
    const auto germanium = detectorIndex_.find(config::kGermaniumType);
    if (germanium != detectorIndex_.end()) {
        detectorHistograms_[germanium->second].setCalibratedEnergyAxis();
    }
    gammaCoincidences_.setCalibratedEnergyAxes();
    angularCoincidences_.setCalibratedEnergyAxes();
    individualGermaniumHistograms_.setCalibratedEnergyAxes();
    germaniumConditionHistograms_.setCalibratedEnergyAxes();
}

void RawAnalysis::processReader(TTreeReader& reader,
                                ProgressReporter* progressReporter)
{
    TTreeReaderValue<std::vector<UShort_t>> detectorID(
        reader, config::kDetectorIDBranch);
    TTreeReaderValue<std::vector<UShort_t>> detectorType(
        reader, config::kDetectorTypeBranch);
    TTreeReaderValue<std::vector<UShort_t>> energy(
        reader, config::kEnergyBranch);
    TTreeReaderValue<std::vector<UShort_t>> psd(
        reader, config::kPsdBranch);
    TTreeReaderValue<std::vector<UShort_t>> relativeNsTime(
        reader, config::kRelativeNsBranch);
    TTreeReaderValue<std::vector<UShort_t>> relativePsTime(
        reader, config::kRelativePsBranch);
    TTreeReaderValue<ULong_t> absoluteTime(
        reader, config::kAbsoluteTimeBranch);

    const GermaniumCalibration* activeCalibration = nullptr;
    std::string activeFileName;
    TFile* activeInputFile = nullptr;
    const RunningTimeMap::FileRange* activeTimeRange = nullptr;
    if (!calibrationManager_.usesRunRanges()) {
        activeCalibration = calibrationManager_.calibrationForRun(0);
    }

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
        TTree* tree = reader.GetTree();
        TFile* inputFile = tree != nullptr ? tree->GetCurrentFile() : nullptr;
        if (inputFile == nullptr) {
            throw std::runtime_error(
                "TTreeReader did not provide the current input filename");
        }
        if (inputFile != activeInputFile) {
            activeInputFile = inputFile;
            activeFileName = inputFile->GetName();
            activeTimeRange = &runningTimeMap_.rangeForFile(activeFileName);
            if (calibrationManager_.usesRunRanges()) {
                const unsigned int run =
                    RunCalibrationManager::runNumberFromFileName(
                        activeFileName);
                activeCalibration = calibrationManager_.calibrationForRun(run);
            }
        }
        const double runningTimeSeconds = RunningTimeMap::runningTimeSeconds(
            *activeTimeRange, *absoluteTime);

        const std::size_t hitCount = detectorType->size();
        const bool vectorSizesAgree =
            detectorID->size()     == hitCount &&
            energy->size()         == hitCount &&
            psd->size()            == hitCount &&
            relativeNsTime->size() == hitCount &&
            relativePsTime->size() == hitCount;

        if (!vectorSizesAgree) {
            ++malformedEvents_;
            if (diagnosticsEnabled_) {
                std::cerr << "Warning: branch-vector size mismatch in event "
                          << reader.GetCurrentEntry()
                          << "; event skipped.\n";
            }
            continue;
        }

        firstAbsoluteTime_ = std::min(firstAbsoluteTime_, *absoluteTime);
        lastAbsoluteTime_ = std::max(lastAbsoluteTime_, *absoluteTime);

        std::vector<double> relativeTimesNs(hitCount, 0.0);
        std::vector<bool> acceptedGermanium(hitCount, false);
        std::unordered_set<UShort_t> inTimeBgoIDs;
        bool siliconCoincident = false;

        for (std::size_t hit = 0; hit < hitCount; ++hit) {
            const double relativeTimePs =
                1000.0 * relativeNsTime->at(hit) + relativePsTime->at(hit);
            const double relativeTimeNs = relativeTimePs / 1000.0;
            relativeTimesNs[hit] = relativeTimeNs;

            const UShort_t type = detectorType->at(hit);
            if (type == config::kGermaniumType) {
                if (excludedGermaniumIDs_.count(detectorID->at(hit)) != 0) {
                    continue;
                }
                gateStatistics_.recordGermanium(insideInclusive(
                    relativeTimeNs, config::kGermaniumTimeMinNs,
                    config::kGermaniumTimeMaxNs));
            } else if (type == config::kBgoType) {
                const bool passesTiming = insideInclusive(
                    relativeTimeNs, config::kBgoVetoTimeMinNs,
                    config::kBgoVetoTimeMaxNs);
                gateStatistics_.recordBgo(passesTiming);
                if (passesTiming) {
                    inTimeBgoIDs.insert(detectorID->at(hit));
                }
            } else if (type == config::kSiliconType) {
                const bool passesTiming = insideInclusive(
                    relativeTimeNs, config::kSiliconTimeMinNs,
                    config::kSiliconTimeMaxNs);
                const bool passesEnergy = insideInclusive(
                    energy->at(hit), config::kSiliconEnergyMin,
                    config::kSiliconEnergyMax);
                gateStatistics_.recordSilicon(passesTiming, passesEnergy);
                siliconCoincident = siliconCoincident ||
                    (passesTiming && passesEnergy);
            }
        }
        gateStatistics_.recordEvent(siliconCoincident);

        std::vector<double> analysedEnergies(energy->begin(), energy->end());
        std::vector<double> gammaEnergies;
        std::vector<unsigned short> gammaIDs;
        gammaEnergies.reserve(hitCount);
        gammaIDs.reserve(hitCount);

        bool rejectEvent = false;
        for (std::size_t hit = 0; hit < hitCount; ++hit) {
            if (detectorType->at(hit) != config::kGermaniumType ||
                excludedGermaniumIDs_.count(detectorID->at(hit)) != 0 ||
                !insideInclusive(relativeTimesNs[hit],
                                 config::kGermaniumTimeMinNs,
                                 config::kGermaniumTimeMaxNs)) {
                continue;
            }

            const GermaniumCalibration::Result calibrated =
                activeCalibration != nullptr
                    ? activeCalibration->calibrate(
                          detectorID->at(hit), energy->at(hit))
                    : GermaniumCalibration::Result{
                          static_cast<double>(energy->at(hit)),
                          GermaniumCalibration::Failure::None};

            if (!calibrated.valid()) {
                rejectEvent = true;
                if (calibrated.failure ==
                    GermaniumCalibration::Failure::MissingID) {
                    ++missingCalibrationEvents_;
                } else if (calibrated.failure ==
                           GermaniumCalibration::Failure::OutOfRange) {
                    ++outOfRangeCalibrationEvents_;
                } else {
                    ++nonFiniteCalibrationEvents_;
                }
                break;
            }

            analysedEnergies[hit] = calibrated.energy;
            acceptedGermanium[hit] = true;
            gammaEnergies.push_back(calibrated.energy);
            gammaIDs.push_back(detectorID->at(hit));
        }

        if (rejectEvent) {
            continue;
        }

        germaniumConditionHistograms_.fillUnconditionedCoincidences(
            gammaEnergies);

        std::vector<unsigned int> multiplicities(
            detectorHistograms_.size(), 0U);
        std::vector<unsigned short> bgoIDs;
        bgoIDs.reserve(hitCount);
        for (std::size_t hit = 0; hit < hitCount; ++hit) {
            if (detectorType->at(hit) == config::kBgoType) {
                bgoIDs.push_back(detectorID->at(hit));
            }
        }
        const DetectorMatchingStatistics::EventResult matchingResult =
            detectorMatchingStatistics_.record(gammaIDs, bgoIDs);
        const unsigned int germaniumMultiplicity =
            static_cast<unsigned int>(gammaIDs.size());
        const unsigned int bgoMultiplicity =
            static_cast<unsigned int>(bgoIDs.size());
        const bool foldValid = foldStatistics_.record(
            germaniumMultiplicity, bgoMultiplicity,
            matchingResult.extraBgoMultiplicity == 0);

        std::vector<double> bgoVetoedGammaEnergies;
        std::vector<unsigned short> bgoVetoedGammaIDs;
        bgoVetoedGammaEnergies.reserve(gammaIDs.size());
        bgoVetoedGammaIDs.reserve(gammaIDs.size());
        unsigned int germaniumMultiplicityAfterBgoVeto = 0;

        for (std::size_t hit = 0; hit < hitCount; ++hit) {
            if (psd->at(hit) != 0) {
                ++nonzeroPsdHits_;
            }

            const auto found = detectorIndex_.find(detectorType->at(hit));
            if (found == detectorIndex_.end()) {
                ++unknownHits_;
                unknownDetectorTypes_->Fill(detectorType->at(hit));
                continue;
            }

            const bool isGermanium =
                detectorType->at(hit) == config::kGermaniumType;
            if (isGermanium && !acceptedGermanium[hit]) {
                continue;
            }

            const std::size_t detectorIndex = found->second;
            detectorHistograms_[detectorIndex].fillHit(
                detectorID->at(hit), analysedEnergies[hit],
                relativeTimesNs[hit]);
            if (isGermanium) {
                individualGermaniumHistograms_.fill(
                    detectorID->at(hit), analysedEnergies[hit],
                    energy->at(hit), runningTimeSeconds);
                const bool survivesBgoVeto =
                    inTimeBgoIDs.count(detectorID->at(hit)) == 0;
                gateStatistics_.recordBgoVetoDecision(survivesBgoVeto);
                germaniumConditionHistograms_.fillHit(
                    analysedEnergies[hit], runningTimeSeconds,
                    survivesBgoVeto, siliconCoincident, foldValid);
                if (survivesBgoVeto) {
                    ++germaniumMultiplicityAfterBgoVeto;
                    bgoVetoedGammaEnergies.push_back(analysedEnergies[hit]);
                    bgoVetoedGammaIDs.push_back(detectorID->at(hit));
                }
            } else if (detectorType->at(hit) == config::kSiliconType) {
                individualSiliconHistograms_.fill(
                    detectorID->at(hit), analysedEnergies[hit]);
            } else if (detectorType->at(hit) == config::kBgoType) {
                individualBgoHistograms_.fill(
                    detectorID->at(hit), analysedEnergies[hit]);
            } else if (detectorType->at(hit) == config::kLabrType) {
                individualLabrHistograms_.fill(
                    detectorID->at(hit), analysedEnergies[hit]);
            }
            ++multiplicities[detectorIndex];
        }

        germaniumConditionHistograms_.fillEvent(
            germaniumMultiplicityAfterBgoVeto, siliconCoincident, foldValid);
        gammaCoincidences_.fillEvent(
            bgoVetoedGammaEnergies, siliconCoincident);
        if (siliconCoincident) {
            germaniumConditionHistograms_.fillBgoVetoedSiliconCoincidences(
                bgoVetoedGammaEnergies);
            angularCoincidences_.fillEvent(
                bgoVetoedGammaEnergies, bgoVetoedGammaIDs);
        }

        for (std::size_t index = 0;
             index < detectorHistograms_.size(); ++index) {
            detectorHistograms_[index].fillMultiplicity(
                multiplicities[index]);
        }
        eventMultiplicity_->Fill(hitCount);
        combinedStatistics_.record(multiplicities);

    }
    if (progressReporter != nullptr) progressReporter->add(pendingProgress);
}

void RawAnalysis::merge(const RawAnalysis& other)
{
    if (detectorHistograms_.size() != other.detectorHistograms_.size()) {
        throw std::logic_error("Cannot merge incompatible analysis states");
    }
    for (std::size_t index = 0;
         index < detectorHistograms_.size(); ++index) {
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
    processedEvents_ += other.processedEvents_;
    malformedEvents_ += other.malformedEvents_;
    unknownHits_ += other.unknownHits_;
    nonzeroPsdHits_ += other.nonzeroPsdHits_;
    missingCalibrationEvents_ += other.missingCalibrationEvents_;
    outOfRangeCalibrationEvents_ += other.outOfRangeCalibrationEvents_;
    nonFiniteCalibrationEvents_ += other.nonFiniteCalibrationEvents_;
    firstAbsoluteTime_ = std::min(firstAbsoluteTime_, other.firstAbsoluteTime_);
    lastAbsoluteTime_ = std::max(lastAbsoluteTime_, other.lastAbsoluteTime_);
}

int RawAnalysis::run(const std::vector<std::string>& inputPatterns,
                     const std::string& outputFileName)
{
    const auto totalStart = std::chrono::steady_clock::now();

    TChain chain(config::kTreeName);
    for (const auto& pattern : inputPatterns) {
        if (chain.Add(pattern.c_str()) == 0) {
            std::cerr << "Warning: no file matched '" << pattern << "'.\n";
        }
    }
    if (chain.GetNtrees() == 0) {
        std::cerr << "Error: no input trees called '" << config::kTreeName
                  << "' were added.\n";
        return 2;
    }

    std::vector<std::string> inputFiles;
    TObjArray* fileElements = chain.GetListOfFiles();
    inputFiles.reserve(static_cast<std::size_t>(fileElements->GetEntries()));
    for (int index = 0; index < fileElements->GetEntries(); ++index) {
        const auto* element =
            dynamic_cast<const TChainElement*>(fileElements->At(index));
        if (element != nullptr) {
            inputFiles.emplace_back(element->GetTitle());
        }
    }

    std::cout << "Reading " << chain.GetNtrees() << " file(s), "
              << chain.GetEntries() << " events.\n";

    try {
        for (const std::string& fileName : inputFiles) {
            const auto range = rawAbsoluteTimeRange(fileName);
            runningTimeMap_.addFile(fileName, range.first, range.second);
        }
    } catch (const std::exception& error) {
        std::cerr << "Absolute-time scan error: " << error.what() << "\n";
        return 5;
    }
    germaniumConditionHistograms_.setRunningTimeRange(
        runningTimeMap_.totalSeconds());
    individualGermaniumHistograms_.setRunningTimeRange(
        runningTimeMap_.totalSeconds());
    std::cout << "Total running time from " << runningTimeMap_.fileCount()
              << " file(s): " << runningTimeMap_.totalSeconds() << " s.\n";

    if (!calibrationManager_.empty()) {
        configureCalibratedAxes();
        std::cout << "Germanium calibration is enabled"
                  << (calibrationManager_.usesRunRanges()
                          ? " with run-dependent ranges.\n"
                          : ".\n");
        if (calibrationManager_.usesRunRanges()) {
            try {
                for (const std::string& fileName : inputFiles) {
                    const unsigned int run =
                        RunCalibrationManager::runNumberFromFileName(fileName);
                    calibrationManager_.calibrationForRun(run);
                }
            } catch (const std::exception& error) {
                std::cerr << "Calibration selection error: "
                          << error.what() << "\n";
                return 5;
            }
        }
    }

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
            std::mutex workerMutex;
            std::unordered_map<std::thread::id, RawAnalysis*> workerByThread;
            std::vector<std::unique_ptr<RawAnalysis>> workers;

            processor.Process([&](TTreeReader& reader) {
                RawAnalysis* worker = nullptr;
                {
                    std::lock_guard<std::mutex> lock(workerMutex);
                    const std::thread::id id = std::this_thread::get_id();
                    const auto existing = workerByThread.find(id);
                    if (existing != workerByThread.end()) {
                        worker = existing->second;
                    } else {
                        auto state = std::make_unique<RawAnalysis>();
                        state->diagnosticsEnabled_ = false;
                        state->excludedGermaniumIDs_ = excludedGermaniumIDs_;
                        state->calibrationManager_ = calibrationManager_;
                        state->runningTimeMap_ = runningTimeMap_;
                        state->germaniumConditionHistograms_
                            .setRunningTimeRange(runningTimeMap_.totalSeconds());
                        state->individualGermaniumHistograms_
                            .setRunningTimeRange(runningTimeMap_.totalSeconds());
                        if (!calibrationManager_.empty()) {
                            state->configureCalibratedAxes();
                        }
                        worker = state.get();
                        workers.push_back(std::move(state));
                        workerByThread.emplace(id, worker);
                    }
                }
                worker->processReader(reader, &progressReporter);
            });

            for (const auto& worker : workers) {
                merge(*worker);
            }
        }
    } catch (const std::exception& error) {
        std::cerr << "Analysis error: " << error.what() << "\n";
        return 5;
    }
    progressReporter.finish();
    const double processingSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - processingStart).count();

    std::cout << "Processed " << processedEvents_
              << " events.             \n";

    TFile outputFile(outputFileName.c_str(), "RECREATE");
    if (outputFile.IsZombie()) {
        std::cerr << "Error: could not create '" << outputFileName << "'.\n";
        return 3;
    }

    outputFile.cd();
    eventMultiplicity_->Write();
    unknownDetectorTypes_->Write();
    std::unordered_map<unsigned short, TDirectory*> detectorDirectories;
    for (const auto& histograms : detectorHistograms_) {
        detectorDirectories.emplace(
            histograms.detectorType(), histograms.write(outputFile));
    }

    TDirectory* germaniumDirectory =
        detectorDirectories.at(config::kGermaniumType);
    TDirectory* siliconDirectory = detectorDirectories.at(config::kSiliconType);
    TDirectory* bgoDirectory = detectorDirectories.at(config::kBgoType);
    TDirectory* labrDirectory = detectorDirectories.at(config::kLabrType);
    TDirectory* germaniumEnergyDirectory =
        germaniumDirectory->GetDirectory("Energy");
    TDirectory* germaniumTimeDirectory =
        germaniumDirectory->GetDirectory("Time");
    TDirectory* siliconEnergyDirectory =
        siliconDirectory->GetDirectory("Energy");
    TDirectory* bgoEnergyDirectory = bgoDirectory->GetDirectory("Energy");
    TDirectory* labrEnergyDirectory = labrDirectory->GetDirectory("Energy");
    if (germaniumEnergyDirectory == nullptr ||
        germaniumTimeDirectory == nullptr ||
        siliconEnergyDirectory == nullptr || bgoEnergyDirectory == nullptr ||
        labrEnergyDirectory == nullptr) {
        std::cerr << "Error: detector Energy/Time directory is missing.\n";
        return 4;
    }

    germaniumConditionHistograms_.write(
        *germaniumDirectory, *germaniumEnergyDirectory,
        *germaniumTimeDirectory);
    individualGermaniumHistograms_.write(
        *germaniumEnergyDirectory, *germaniumTimeDirectory);
    individualSiliconHistograms_.write(*siliconEnergyDirectory);
    individualBgoHistograms_.write(*bgoEnergyDirectory);
    individualLabrHistograms_.write(*labrEnergyDirectory);
    gammaCoincidences_.write(outputFile);
    TDirectory* gammaCoincidenceDirectory =
        outputFile.GetDirectory("Coincidences");
    if (gammaCoincidenceDirectory == nullptr) {
        std::cerr << "Error: gamma-coincidence directory was not created.\n";
        return 4;
    }
    angularCoincidences_.write(*gammaCoincidenceDirectory);
    outputFile.Close();

    std::cout << "Wrote histograms to " << outputFileName << ".\n";
    if (diagnosticsEnabled_) {
        std::cout << "Malformed events skipped: " << malformedEvents_ << "\n"
                  << "Hits with unknown detectorType: " << unknownHits_ << "\n"
                  << "Hits with non-zero psd: " << nonzeroPsdHits_ << "\n"
                  << "Events rejected for missing Ge calibration: "
                  << missingCalibrationEvents_ << "\n"
                  << "Events rejected outside Ge calibration range: "
                  << outOfRangeCalibrationEvents_ << "\n"
                  << "Events rejected for non-finite Ge calibration: "
                  << nonFiniteCalibrationEvents_ << "\n";
        if (processedEvents_ > malformedEvents_) {
            std::cout << "Raw absoluteTime range: " << firstAbsoluteTime_
                      << " to " << lastAbsoluteTime_ << "\n";
        }
        runningTimeMap_.print(std::cout);
        combinedStatistics_.print(std::cout, "combined accepted events");
        foldStatistics_.print(std::cout);
        detectorMatchingStatistics_.print(std::cout);
        gateStatistics_.print(std::cout);
    }

    const double totalSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - totalStart).count();
    const double processingRate = processingSeconds > 0.0
        ? static_cast<double>(processedEvents_) / processingSeconds : 0.0;
    const double totalRate = totalSeconds > 0.0
        ? static_cast<double>(processedEvents_) / totalSeconds : 0.0;
    std::cout << "\n=== Analysis timing ===\n"
              << "Event processing: " << formatElapsed(processingSeconds)
              << " (" << std::fixed << std::setprecision(0)
              << processingRate << " events/s)\n"
              << "Total analysis:   " << formatElapsed(totalSeconds)
              << " (" << totalRate << " events/s including setup, merge, "
              << "output, and diagnostics)\n";
    return 0;
}
