#include "RawAnalysis.h"

#include "AnalysisConfig.h"
#include "DetectorMatchingStatistics.h"
#include "FoldStatistics.h"
#include "GateStatistics.h"
#include "MultiplicityStatistics.h"

#include <RtypesCore.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include <algorithm>
#include <iostream>
#include <limits>
#include <stdexcept>
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

} // namespace

RawAnalysis::RawAnalysis()
    : individualSiliconHistograms_(definitionFor(config::kSiliconType)),
      individualBgoHistograms_(definitionFor(config::kBgoType)),
      individualLabrHistograms_(definitionFor(config::kLabrType))
{
    detectorHistograms_.reserve(config::kDetectors.size());
    for (const auto& definition : config::kDetectors) {
        detectorIndex_.emplace(definition.type, detectorHistograms_.size());
        detectorHistograms_.emplace_back(definition);
    }
}

void RawAnalysis::addCalFile(const std::string& fileName)
{
    calibrationManager_.addGlobalCalFile(fileName);
}

void RawAnalysis::addMcalFile(const std::string& fileName)
{
    calibrationManager_.addGlobalMcalFile(fileName);
}

void RawAnalysis::addRunCalFile(
    unsigned int firstRun, unsigned int lastRun,
    const std::string& fileName)
{
    calibrationManager_.addRunCalFile(firstRun, lastRun, fileName);
}

void RawAnalysis::addRunMcalFile(
    unsigned int firstRun, unsigned int lastRun,
    const std::string& fileName)
{
    calibrationManager_.addRunMcalFile(firstRun, lastRun, fileName);
}

void RawAnalysis::setDiagnosticsEnabled(bool enabled)
{
    diagnosticsEnabled_ = enabled;
}

void RawAnalysis::excludeGermaniumLUT(unsigned int detectorLUT)
{
    if (detectorLUT >= config::kGermaniumLutBins) {
        throw std::runtime_error(
            "germanium detector LUT must be between 0 and " +
            std::to_string(config::kGermaniumLutBins - 1));
    }
    excludedGermaniumLUTs_.insert(
        static_cast<unsigned short>(detectorLUT));
}

int RawAnalysis::run(const std::vector<std::string>& inputPatterns,
                     const std::string& outputFileName)
{
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

    std::cout << "Reading " << chain.GetNtrees() << " file(s), "
              << chain.GetEntries() << " events.\n";

    if (!calibrationManager_.empty()) {
        const auto germanium = detectorIndex_.find(config::kGermaniumType);
        if (germanium != detectorIndex_.end()) {
            detectorHistograms_[germanium->second].setCalibratedEnergyAxis();
        }
        gammaCoincidences_.setCalibratedEnergyAxes();
        angularCoincidences_.setCalibratedEnergyAxes();
        individualGermaniumHistograms_.setCalibratedEnergyAxes();
        germaniumConditionHistograms_.setCalibratedEnergyAxes();
        std::cout << "Germanium calibration is enabled"
                  << (calibrationManager_.usesRunRanges()
                          ? " with run-dependent ranges.\n"
                          : ".\n");
    }

    TTreeReader reader(&chain);
    TTreeReaderValue<std::vector<UShort_t>> detectorLUT(
        reader, config::kDetectorLUTBranch);
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

    // events->Print() reports absoluteTime/g: a scalar ULong_t leaf.
    TTreeReaderValue<ULong_t> absoluteTime(reader, config::kAbsoluteTimeBranch);

    TH1D eventMultiplicity(
        "h1_EventMultiplicity",
        "Event multiplicity;All hits in event;Events",
        config::kMultiplicityBins,
        config::kMultiplicityMin,
        config::kMultiplicityMax);
    TH1D unknownDetectorTypes(
        "h1_UnknownDetectorType",
        "Unrecognised detector type;detectorType;Hits",
        65536, -0.5, 65535.5);
    eventMultiplicity.SetDirectory(nullptr);
    unknownDetectorTypes.SetDirectory(nullptr);
    eventMultiplicity.SetOption("HIST");
    unknownDetectorTypes.SetOption("HIST");

    ULong64_t processedEvents = 0;
    ULong64_t malformedEvents = 0;
    ULong64_t unknownHits = 0;
    ULong64_t nonzeroPsdHits = 0;
    ULong64_t missingCalibrationEvents = 0;
    ULong64_t outOfRangeCalibrationEvents = 0;
    ULong64_t nonFiniteCalibrationEvents = 0;
    ULong_t firstAbsoluteTime = std::numeric_limits<ULong_t>::max();
    ULong_t lastAbsoluteTime = 0;

    MultiplicityStatistics combinedStatistics;
    FoldStatistics foldStatistics;
    GateStatistics gateStatistics;
    DetectorMatchingStatistics detectorMatchingStatistics;
    const GermaniumCalibration* activeCalibration = nullptr;
    int activeTreeNumber = -1;
    if (!calibrationManager_.usesRunRanges()) {
        activeCalibration = calibrationManager_.calibrationForRun(0);
    }

    while (reader.Next()) {
        ++processedEvents;

        if (calibrationManager_.usesRunRanges() &&
            chain.GetTreeNumber() != activeTreeNumber) {
            activeTreeNumber = chain.GetTreeNumber();
            try {
                const TFile* inputFile = chain.GetCurrentFile();
                if (inputFile == nullptr) {
                    throw std::runtime_error(
                        "TChain did not provide the current input filename");
                }
                const unsigned int run =
                    RunCalibrationManager::runNumberFromFileName(
                        inputFile->GetName());
                activeCalibration =
                    calibrationManager_.calibrationForRun(run);
                std::cout << "Using "
                          << activeCalibration->numberOfStages()
                          << " calibration stage(s) for run " << run
                          << ".\n";
            } catch (const std::exception& error) {
                std::cerr << "Calibration selection error: "
                          << error.what() << "\n";
                return 5;
            }
        }

        const std::size_t hitCount = detectorType->size();
        const bool vectorSizesAgree =
            detectorLUT->size()    == hitCount &&
            energy->size()         == hitCount &&
            psd->size()            == hitCount &&
            relativeNsTime->size() == hitCount &&
            relativePsTime->size() == hitCount;

        if (!vectorSizesAgree) {
            ++malformedEvents;
            if (diagnosticsEnabled_) {
                std::cerr << "Warning: branch-vector size mismatch in event "
                          << reader.GetCurrentEntry()
                          << "; event skipped.\n";
            }
            continue;
        }

        firstAbsoluteTime = std::min(firstAbsoluteTime, *absoluteTime);
        lastAbsoluteTime = std::max(lastAbsoluteTime, *absoluteTime);

        std::vector<double> relativeTimesNs(hitCount, 0.0);
        std::vector<bool> acceptedGermanium(hitCount, false);
        std::unordered_set<UShort_t> inTimeBgoLUTs;
        bool siliconCoincident = false;

        for (std::size_t hit = 0; hit < hitCount; ++hit) {
            const double relativeTimePs =
                1000.0 * relativeNsTime->at(hit) + relativePsTime->at(hit);
            const double relativeTimeNs = relativeTimePs / 1000.0;
            relativeTimesNs[hit] = relativeTimeNs;

            const UShort_t type = detectorType->at(hit);
            if (type == config::kGermaniumType) {
                if (excludedGermaniumLUTs_.count(detectorLUT->at(hit)) != 0) {
                    continue;
                }
                gateStatistics.recordGermanium(insideInclusive(
                    relativeTimeNs,
                    config::kGermaniumTimeMinNs,
                    config::kGermaniumTimeMaxNs));
            } else if (type == config::kBgoType) {
                const bool passesTiming = insideInclusive(
                    relativeTimeNs,
                    config::kBgoVetoTimeMinNs,
                    config::kBgoVetoTimeMaxNs);
                gateStatistics.recordBgo(passesTiming);
                if (passesTiming) {
                    inTimeBgoLUTs.insert(detectorLUT->at(hit));
                }
            } else if (type == config::kSiliconType) {
                const bool passesTiming = insideInclusive(
                    relativeTimeNs,
                    config::kSiliconTimeMinNs,
                    config::kSiliconTimeMaxNs);
                const bool passesEnergy = insideInclusive(
                    energy->at(hit),
                    config::kSiliconEnergyMin,
                    config::kSiliconEnergyMax);
                gateStatistics.recordSilicon(passesTiming, passesEnergy);
                siliconCoincident = siliconCoincident ||
                    (passesTiming && passesEnergy);
            }
        }
        gateStatistics.recordEvent(siliconCoincident);

        std::vector<double> analysedEnergies(energy->begin(), energy->end());
        std::vector<double> gammaEnergies;
        std::vector<unsigned short> gammaLUTs;
        gammaEnergies.reserve(hitCount);
        gammaLUTs.reserve(hitCount);

        bool rejectEvent = false;
        for (std::size_t hit = 0; hit < hitCount; ++hit) {
            if (detectorType->at(hit) != config::kGermaniumType) {
                continue;
            }
            if (excludedGermaniumLUTs_.count(detectorLUT->at(hit)) != 0) {
                continue;
            }
            if (!insideInclusive(relativeTimesNs[hit],
                                 config::kGermaniumTimeMinNs,
                                 config::kGermaniumTimeMaxNs)) {
                continue;
            }

            const GermaniumCalibration::Result calibrated =
                activeCalibration != nullptr
                    ? activeCalibration->calibrate(
                          detectorLUT->at(hit), energy->at(hit))
                    : GermaniumCalibration::Result{
                          static_cast<double>(energy->at(hit)),
                          GermaniumCalibration::Failure::None};

            if (!calibrated.valid()) {
                rejectEvent = true;
                if (calibrated.failure ==
                    GermaniumCalibration::Failure::MissingLUT) {
                    ++missingCalibrationEvents;
                } else if (calibrated.failure ==
                           GermaniumCalibration::Failure::OutOfRange) {
                    ++outOfRangeCalibrationEvents;
                } else {
                    ++nonFiniteCalibrationEvents;
                }
                break;
            }

            analysedEnergies[hit] = calibrated.energy;
            acceptedGermanium[hit] = true;
            gammaEnergies.push_back(calibrated.energy);
            gammaLUTs.push_back(detectorLUT->at(hit));
        }

        if (rejectEvent) {
            continue;
        }

        germaniumConditionHistograms_.fillUnconditionedCoincidences(
            gammaEnergies);

        std::vector<unsigned int> multiplicities(
            detectorHistograms_.size(), 0U);
        std::vector<unsigned short> bgoLUTs;
        bgoLUTs.reserve(hitCount);
        for (std::size_t hit = 0; hit < hitCount; ++hit) {
            if (detectorType->at(hit) == config::kBgoType) {
                bgoLUTs.push_back(detectorLUT->at(hit));
            }
        }
        const DetectorMatchingStatistics::EventResult matchingResult =
            detectorMatchingStatistics.record(gammaLUTs, bgoLUTs);
        const unsigned int germaniumMultiplicity =
            static_cast<unsigned int>(gammaLUTs.size());
        const unsigned int bgoMultiplicity =
            static_cast<unsigned int>(bgoLUTs.size());
        const bool allBgoHitsMatched =
            matchingResult.extraBgoMultiplicity == 0;
        const bool foldValid = foldStatistics.record(
            germaniumMultiplicity, bgoMultiplicity, allBgoHitsMatched);

        std::vector<double> bgoVetoedGammaEnergies;
        std::vector<unsigned short> bgoVetoedGammaLUTs;
        bgoVetoedGammaEnergies.reserve(gammaLUTs.size());
        bgoVetoedGammaLUTs.reserve(gammaLUTs.size());
        unsigned int germaniumMultiplicityAfterBgoVeto = 0;

        for (std::size_t hit = 0; hit < hitCount; ++hit) {
            if (psd->at(hit) != 0) {
                ++nonzeroPsdHits;
            }

            const auto found = detectorIndex_.find(detectorType->at(hit));
            if (found == detectorIndex_.end()) {
                ++unknownHits;
                unknownDetectorTypes.Fill(detectorType->at(hit));
                continue;
            }

            const bool isGermanium =
                detectorType->at(hit) == config::kGermaniumType;
            if (isGermanium && !acceptedGermanium[hit]) {
                continue;
            }

            const std::size_t detectorIndex = found->second;
            detectorHistograms_[detectorIndex].fillHit(
                detectorLUT->at(hit), analysedEnergies[hit],
                relativeTimesNs[hit]);
            if (isGermanium) {
                individualGermaniumHistograms_.fill(
                    detectorLUT->at(hit), analysedEnergies[hit]);
                const bool survivesBgoVeto =
                    inTimeBgoLUTs.count(detectorLUT->at(hit)) == 0;
                gateStatistics.recordBgoVetoDecision(survivesBgoVeto);
                germaniumConditionHistograms_.fillHit(
                    analysedEnergies[hit],
                    static_cast<double>(*absoluteTime),
                    survivesBgoVeto,
                    siliconCoincident,
                    foldValid);
                if (survivesBgoVeto) {
                    ++germaniumMultiplicityAfterBgoVeto;
                    bgoVetoedGammaEnergies.push_back(analysedEnergies[hit]);
                    bgoVetoedGammaLUTs.push_back(detectorLUT->at(hit));
                }
            } else if (detectorType->at(hit) == config::kSiliconType) {
                individualSiliconHistograms_.fill(
                    detectorLUT->at(hit), analysedEnergies[hit]);
            } else if (detectorType->at(hit) == config::kBgoType) {
                individualBgoHistograms_.fill(
                    detectorLUT->at(hit), analysedEnergies[hit]);
            } else if (detectorType->at(hit) == config::kLabrType) {
                individualLabrHistograms_.fill(
                    detectorLUT->at(hit), analysedEnergies[hit]);
            }
            ++multiplicities[detectorIndex];
        }

        germaniumConditionHistograms_.fillEvent(
            germaniumMultiplicityAfterBgoVeto,
            siliconCoincident, foldValid);
        gammaCoincidences_.fillEvent(
            bgoVetoedGammaEnergies, siliconCoincident);
        if (siliconCoincident) {
            germaniumConditionHistograms_.fillBgoVetoedSiliconCoincidences(
                bgoVetoedGammaEnergies);
            angularCoincidences_.fillEvent(
                bgoVetoedGammaEnergies, bgoVetoedGammaLUTs);
        }

        for (std::size_t index = 0; index < detectorHistograms_.size(); ++index) {
            detectorHistograms_[index].fillMultiplicity(multiplicities[index]);
        }
        eventMultiplicity.Fill(hitCount);
        combinedStatistics.record(multiplicities);

        if (processedEvents % 100000 == 0) {
            std::cout << "Processed " << processedEvents << " events.\r"
                      << std::flush;
        }
    }

    std::cout << "Processed " << processedEvents << " events.             \n";

    TFile outputFile(outputFileName.c_str(), "RECREATE");
    if (outputFile.IsZombie()) {
        std::cerr << "Error: could not create '" << outputFileName << "'.\n";
        return 3;
    }

    outputFile.cd();
    eventMultiplicity.Write();
    unknownDetectorTypes.Write();
    std::unordered_map<unsigned short, TDirectory*> detectorDirectories;
    for (const auto& histograms : detectorHistograms_) {
        detectorDirectories.emplace(
            histograms.detectorType(), histograms.write(outputFile));
    }

    TDirectory* germaniumDirectory =
        detectorDirectories.at(config::kGermaniumType);
    TDirectory* siliconDirectory =
        detectorDirectories.at(config::kSiliconType);
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
    individualGermaniumHistograms_.write(*germaniumEnergyDirectory);
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
    if (!diagnosticsEnabled_) {
        return 0;
    }

    std::cout << "Malformed events skipped: " << malformedEvents << "\n"
              << "Hits with unknown detectorType: " << unknownHits << "\n"
              << "Hits with non-zero psd: " << nonzeroPsdHits << "\n"
              << "Events rejected for missing Ge calibration: "
              << missingCalibrationEvents << "\n"
              << "Events rejected outside Ge calibration range: "
              << outOfRangeCalibrationEvents << "\n"
              << "Events rejected for non-finite Ge calibration: "
              << nonFiniteCalibrationEvents << "\n";
    if (processedEvents > malformedEvents) {
        std::cout << "Raw absoluteTime range: " << firstAbsoluteTime
                  << " to " << lastAbsoluteTime << "\n";
    }

    combinedStatistics.print(std::cout, "combined accepted events");
    foldStatistics.print(std::cout);
    detectorMatchingStatistics.print(std::cout);
    gateStatistics.print(std::cout);

    return 0;
}
