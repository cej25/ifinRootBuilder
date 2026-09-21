#include "AnalysisTreeBuilder.h"

#include "AnalysisConfig.h"

#include <RtypesCore.h>
#include <TChain.h>
#include <TChainElement.h>
#include <TFile.h>
#include <TNamed.h>
#include <TObjArray.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <filesystem>
#include <iostream>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {

bool insideInclusive(double value, double minimum, double maximum)
{
    return value >= minimum && value <= maximum;
}

unsigned int extraBgoMultiplicity(
    const std::vector<UShort_t>& germaniumIDs,
    const std::vector<UShort_t>& bgoIDs)
{
    std::unordered_map<UShort_t, unsigned int> unmatchedGermanium;
    for (const UShort_t id : germaniumIDs) {
        ++unmatchedGermanium[id];
    }
    unsigned int extra = 0;
    for (const UShort_t id : bgoIDs) {
        auto ge = unmatchedGermanium.find(id);
        if (ge != unmatchedGermanium.end() && ge->second > 0) {
            --ge->second;
        } else {
            ++extra;
        }
    }
    return extra;
}

std::string outputNameFor(const std::string& inputFileName,
                          const std::string& outputDirectory)
{
    const std::filesystem::path input(inputFileName);
    return (std::filesystem::path(outputDirectory) /
            (input.stem().string() + "_analysis.root")).string();
}

std::string configurationDescription(
    const std::unordered_set<unsigned short>& excludedIDs)
{
    std::vector<unsigned short> sortedIDs(
        excludedIDs.begin(), excludedIDs.end());
    std::sort(sortedIDs.begin(), sortedIDs.end());
    std::ostringstream text;
    text << "schema=" << config::kAnalysisTreeSchemaVersion
         << ";geTime=" << config::kGermaniumTimeMinNs << ':'
         << config::kGermaniumTimeMaxNs
         << ";bgoTime=" << config::kBgoVetoTimeMinNs << ':'
         << config::kBgoVetoTimeMaxNs
         << ";siTime=" << config::kSiliconTimeMinNs << ':'
         << config::kSiliconTimeMaxNs
         << ";siEnergy=" << config::kSiliconEnergyMin << ':'
         << config::kSiliconEnergyMax << ";excludedGe=";
    for (std::size_t index = 0; index < sortedIDs.size(); ++index) {
        if (index != 0) text << ',';
        text << sortedIDs[index];
    }
    return text.str();
}

} // namespace

void AnalysisTreeBuilder::addCalFile(const std::string& fileName)
{
    calibrationManager_.addGlobalCalFile(fileName);
}

void AnalysisTreeBuilder::addMcalFile(const std::string& fileName)
{
    calibrationManager_.addGlobalMcalFile(fileName);
}

void AnalysisTreeBuilder::addRunCalFile(
    unsigned int firstRun, unsigned int lastRun,
    const std::string& fileName)
{
    calibrationManager_.addRunCalFile(firstRun, lastRun, fileName);
}

void AnalysisTreeBuilder::addRunMcalFile(
    unsigned int firstRun, unsigned int lastRun,
    const std::string& fileName)
{
    calibrationManager_.addRunMcalFile(firstRun, lastRun, fileName);
}

void AnalysisTreeBuilder::excludeGermaniumID(unsigned int detectorID)
{
    if (detectorID >= config::kGermaniumIdBins) {
        throw std::runtime_error(
            "germanium detector ID must be between 0 and " +
            std::to_string(config::kGermaniumIdBins - 1));
    }
    excludedGermaniumIDs_.insert(static_cast<unsigned short>(detectorID));
}

void AnalysisTreeBuilder::setDiagnosticsEnabled(bool enabled)
{
    diagnosticsEnabled_ = enabled;
}

void AnalysisTreeBuilder::setThreadCount(unsigned int threadCount)
{
    if (threadCount == 0) {
        throw std::runtime_error("thread count must be at least 1");
    }
    threadCount_ = threadCount;
}

AnalysisTreeBuilder::FileResult AnalysisTreeBuilder::processFile(
    const std::string& inputFileName,
    const std::string& outputFileName) const
{
    const unsigned int run = calibrationManager_.usesRunRanges()
        ? RunCalibrationManager::runNumberFromFileName(inputFileName) : 0;
    const GermaniumCalibration* calibration =
        calibrationManager_.calibrationForRun(run);

    TFile inputFile(inputFileName.c_str(), "READ");
    if (inputFile.IsZombie()) {
        throw std::runtime_error("could not open input file '" +
                                 inputFileName + "'");
    }
    TTree* inputTree = nullptr;
    inputFile.GetObject(config::kTreeName, inputTree);
    if (inputTree == nullptr) {
        throw std::runtime_error("input file '" + inputFileName +
                                 "' has no tree called '" +
                                 config::kTreeName + "'");
    }

    TTreeReader reader(inputTree);
    TTreeReaderValue<std::vector<UShort_t>> detectorID(
        reader, config::kDetectorIDBranch);
    TTreeReaderValue<std::vector<UShort_t>> detectorType(
        reader, config::kDetectorTypeBranch);
    TTreeReaderValue<std::vector<UShort_t>> rawEnergy(
        reader, config::kEnergyBranch);
    TTreeReaderValue<std::vector<UShort_t>> psd(reader, config::kPsdBranch);
    TTreeReaderValue<std::vector<UShort_t>> relativeNsTime(
        reader, config::kRelativeNsBranch);
    TTreeReaderValue<std::vector<UShort_t>> relativePsTime(
        reader, config::kRelativePsBranch);
    TTreeReaderValue<ULong_t> rawAbsoluteTime(
        reader, config::kAbsoluteTimeBranch);

    TFile outputFile(outputFileName.c_str(), "RECREATE");
    if (outputFile.IsZombie()) {
        throw std::runtime_error("could not create output file '" +
                                 outputFileName + "'");
    }
    outputFile.cd();
    TTree outputTree(config::kAnalysisTreeName,
                     "Calibrated and condition-tagged detector events");

    ULong64_t absoluteTime = 0;
    UInt_t eventMultiplicity = 0;
    UInt_t geMultiplicity = 0;
    UInt_t geMultiplicityVetoed = 0;
    UInt_t bgoMultiplicity = 0;
    UInt_t siMultiplicity = 0;
    UInt_t labrMultiplicity = 0;
    UInt_t extraBgo = 0;
    Bool_t foldValid = false;
    Bool_t siliconCoincident = false;
    Bool_t geEnergyCalibrated = calibration != nullptr;
    UInt_t geTimingCandidates = 0;
    UInt_t geTimingPass = 0;
    UInt_t bgoTimingPass = 0;
    UInt_t siTimingPass = 0;
    UInt_t siEnergyPass = 0;
    UInt_t siCombinedPass = 0;
    UInt_t nonzeroPsdHits = 0;

    std::vector<UShort_t> geID, bgoID, siID, labrID;
    std::vector<Double_t> geEnergy;
    std::vector<Float_t> geTime, bgoEnergy, bgoTime;
    std::vector<Float_t> siEnergy, siTime, labrEnergy, labrTime;
    std::vector<UChar_t> geSurvivesBgoVeto;
    std::vector<UShort_t> unknownDetectorType;

    outputTree.Branch("absoluteTime", &absoluteTime);
    outputTree.Branch("eventMultiplicity", &eventMultiplicity);
    outputTree.Branch("geMultiplicity", &geMultiplicity);
    outputTree.Branch("geMultiplicityVetoed", &geMultiplicityVetoed);
    outputTree.Branch("bgoMultiplicity", &bgoMultiplicity);
    outputTree.Branch("siMultiplicity", &siMultiplicity);
    outputTree.Branch("labrMultiplicity", &labrMultiplicity);
    outputTree.Branch("extraBgoMultiplicity", &extraBgo);
    outputTree.Branch("foldValid", &foldValid);
    outputTree.Branch("siliconCoincident", &siliconCoincident);
    outputTree.Branch("geEnergyCalibrated", &geEnergyCalibrated);
    outputTree.Branch("geTimingCandidates", &geTimingCandidates);
    outputTree.Branch("geTimingPass", &geTimingPass);
    outputTree.Branch("bgoTimingPass", &bgoTimingPass);
    outputTree.Branch("siTimingPass", &siTimingPass);
    outputTree.Branch("siEnergyPass", &siEnergyPass);
    outputTree.Branch("siCombinedPass", &siCombinedPass);
    outputTree.Branch("nonzeroPsdHits", &nonzeroPsdHits);
    outputTree.Branch("geID", &geID);
    outputTree.Branch("geEnergy", &geEnergy);
    outputTree.Branch("geTime", &geTime);
    outputTree.Branch("geSurvivesBgoVeto", &geSurvivesBgoVeto);
    outputTree.Branch("bgoID", &bgoID);
    outputTree.Branch("bgoEnergy", &bgoEnergy);
    outputTree.Branch("bgoTime", &bgoTime);
    outputTree.Branch("siID", &siID);
    outputTree.Branch("siEnergy", &siEnergy);
    outputTree.Branch("siTime", &siTime);
    outputTree.Branch("labrID", &labrID);
    outputTree.Branch("labrEnergy", &labrEnergy);
    outputTree.Branch("labrTime", &labrTime);
    outputTree.Branch("unknownDetectorType", &unknownDetectorType);
    outputTree.SetAutoFlush(-100000000LL);

    FileResult result;
    while (reader.Next()) {
        ++result.inputEvents;
        const std::size_t hitCount = detectorType->size();
        if (detectorID->size() != hitCount || rawEnergy->size() != hitCount ||
            psd->size() != hitCount || relativeNsTime->size() != hitCount ||
            relativePsTime->size() != hitCount) {
            ++result.malformedEvents;
            continue;
        }

        absoluteTime = static_cast<ULong64_t>(*rawAbsoluteTime);
        eventMultiplicity = static_cast<UInt_t>(hitCount);
        geMultiplicity = geMultiplicityVetoed = 0;
        bgoMultiplicity = siMultiplicity = labrMultiplicity = 0;
        extraBgo = 0;
        foldValid = false;
        siliconCoincident = false;
        geTimingCandidates = geTimingPass = bgoTimingPass = 0;
        siTimingPass = siEnergyPass = siCombinedPass = 0;
        nonzeroPsdHits = 0;
        geID.clear(); geEnergy.clear(); geTime.clear();
        geSurvivesBgoVeto.clear();
        bgoID.clear(); bgoEnergy.clear(); bgoTime.clear();
        siID.clear(); siEnergy.clear(); siTime.clear();
        labrID.clear(); labrEnergy.clear(); labrTime.clear();
        unknownDetectorType.clear();

        std::unordered_set<UShort_t> inTimeBgoIDs;
        bool rejectEvent = false;
        for (std::size_t hit = 0; hit < hitCount; ++hit) {
            const double time = static_cast<double>(relativeNsTime->at(hit)) +
                static_cast<double>(relativePsTime->at(hit)) / 1000.0;
            const UShort_t type = detectorType->at(hit);
            const UShort_t id = detectorID->at(hit);
            nonzeroPsdHits += psd->at(hit) != 0;

            if (type == config::kGermaniumType) {
                if (excludedGermaniumIDs_.count(id) != 0) {
                    continue;
                }
                ++geTimingCandidates;
                if (!insideInclusive(time, config::kGermaniumTimeMinNs,
                                     config::kGermaniumTimeMaxNs)) {
                    continue;
                }
                ++geTimingPass;
                const auto calibrated = calibration != nullptr
                    ? calibration->calibrate(id, rawEnergy->at(hit))
                    : GermaniumCalibration::Result{
                          static_cast<double>(rawEnergy->at(hit)),
                          GermaniumCalibration::Failure::None};
                if (!calibrated.valid()) {
                    rejectEvent = true;
                    break;
                }
                geID.push_back(id);
                geEnergy.push_back(calibrated.energy);
                geTime.push_back(static_cast<Float_t>(time));
            } else if (type == config::kBgoType) {
                bgoID.push_back(id);
                bgoEnergy.push_back(static_cast<Float_t>(rawEnergy->at(hit)));
                bgoTime.push_back(static_cast<Float_t>(time));
                if (insideInclusive(time, config::kBgoVetoTimeMinNs,
                                    config::kBgoVetoTimeMaxNs)) {
                    ++bgoTimingPass;
                    inTimeBgoIDs.insert(id);
                }
            } else if (type == config::kSiliconType) {
                siID.push_back(id);
                siEnergy.push_back(static_cast<Float_t>(rawEnergy->at(hit)));
                siTime.push_back(static_cast<Float_t>(time));
                const bool timing = insideInclusive(
                    time, config::kSiliconTimeMinNs,
                    config::kSiliconTimeMaxNs);
                const bool energyPass = insideInclusive(
                    rawEnergy->at(hit), config::kSiliconEnergyMin,
                    config::kSiliconEnergyMax);
                siTimingPass += timing;
                siEnergyPass += energyPass;
                siCombinedPass += timing && energyPass;
                siliconCoincident = siliconCoincident ||
                    (timing && energyPass);
            } else if (type == config::kLabrType) {
                labrID.push_back(id);
                labrEnergy.push_back(static_cast<Float_t>(rawEnergy->at(hit)));
                labrTime.push_back(static_cast<Float_t>(time));
            } else {
                unknownDetectorType.push_back(type);
            }
        }

        if (rejectEvent) {
            ++result.calibrationRejectedEvents;
            continue;
        }

        geMultiplicity = static_cast<UInt_t>(geID.size());
        bgoMultiplicity = static_cast<UInt_t>(bgoID.size());
        siMultiplicity = static_cast<UInt_t>(siID.size());
        labrMultiplicity = static_cast<UInt_t>(labrID.size());
        extraBgo = extraBgoMultiplicity(geID, bgoID);
        foldValid = bgoMultiplicity <= geMultiplicity;
        for (const UShort_t id : geID) {
            const bool survives = inTimeBgoIDs.count(id) == 0;
            geSurvivesBgoVeto.push_back(static_cast<UChar_t>(survives));
            geMultiplicityVetoed += survives;
        }

        outputTree.Fill();
        ++result.writtenEvents;
    }

    outputFile.cd();
    TNamed schemaVersion("AnalysisTreeSchemaVersion",
        std::to_string(config::kAnalysisTreeSchemaVersion).c_str());
    TNamed sourceFile("SourceRawFile", inputFileName.c_str());
    const std::string configuration =
        configurationDescription(excludedGermaniumIDs_);
    TNamed configurationRecord("AnalysisConfiguration",
                               configuration.c_str());
    schemaVersion.Write();
    sourceFile.Write();
    configurationRecord.Write();
    outputTree.Write();
    outputFile.Close();
    inputFile.Close();
    return result;
}

int AnalysisTreeBuilder::run(
    const std::vector<std::string>& inputPatterns,
    const std::string& outputDirectory) const
{
    TChain chain(config::kTreeName);
    for (const std::string& pattern : inputPatterns) {
        if (chain.Add(pattern.c_str()) == 0) {
            std::cerr << "Warning: no file matched '" << pattern << "'.\n";
        }
    }
    if (chain.GetNtrees() == 0) {
        std::cerr << "Error: no input trees called '" << config::kTreeName
                  << "' were found.\n";
        return 2;
    }

    std::filesystem::create_directories(outputDirectory);
    std::vector<std::string> inputFiles;
    TObjArray* elements = chain.GetListOfFiles();
    for (int index = 0; index < elements->GetEntries(); ++index) {
        const auto* element =
            dynamic_cast<const TChainElement*>(elements->At(index));
        if (element != nullptr) {
            inputFiles.emplace_back(element->GetTitle());
        }
    }
    if (inputFiles.empty()) {
        std::cerr << "Error: ROOT did not provide any resolved input files.\n";
        return 2;
    }

    std::unordered_set<std::string> outputNames;
    for (const std::string& input : inputFiles) {
        if (!outputNames.insert(outputNameFor(input, outputDirectory)).second) {
            std::cerr << "Error: input files have duplicate base names.\n";
            return 2;
        }
        try {
            const unsigned int run = calibrationManager_.usesRunRanges()
                ? RunCalibrationManager::runNumberFromFileName(input) : 0;
            calibrationManager_.calibrationForRun(run);
        } catch (const std::exception& error) {
            std::cerr << "Calibration selection error: " << error.what()
                      << "\n";
            return 5;
        }
    }

    const auto start = std::chrono::steady_clock::now();
    const unsigned int workers = std::min<unsigned int>(
        threadCount_, static_cast<unsigned int>(inputFiles.size()));
    std::cout << "Building " << inputFiles.size()
              << " analysis-tree file(s) with " << workers
              << " file worker(s).\n";

    std::atomic<std::size_t> nextFile{0};
    std::atomic<bool> failed{false};
    std::mutex outputMutex;
    std::vector<std::string> errors;
    std::vector<FileResult> results(inputFiles.size());
    std::vector<std::thread> threads;
    threads.reserve(workers);
    for (unsigned int worker = 0; worker < workers; ++worker) {
        threads.emplace_back([&]() {
            while (!failed.load()) {
                const std::size_t index = nextFile.fetch_add(1);
                if (index >= inputFiles.size()) {
                    break;
                }
                const std::string output =
                    outputNameFor(inputFiles[index], outputDirectory);
                try {
                    results[index] = processFile(inputFiles[index], output);
                    std::lock_guard<std::mutex> lock(outputMutex);
                    std::cout << "Built " << output << " ("
                              << results[index].writtenEvents << " events).\n";
                } catch (const std::exception& error) {
                    failed = true;
                    std::lock_guard<std::mutex> lock(outputMutex);
                    errors.push_back(inputFiles[index] + ": " + error.what());
                }
            }
        });
    }
    for (std::thread& thread : threads) {
        thread.join();
    }
    if (!errors.empty()) {
        for (const std::string& error : errors) {
            std::cerr << "Error: " << error << "\n";
        }
        return 5;
    }

    unsigned long long inputEvents = 0;
    unsigned long long writtenEvents = 0;
    unsigned long long malformedEvents = 0;
    unsigned long long rejectedEvents = 0;
    for (const FileResult& result : results) {
        inputEvents += result.inputEvents;
        writtenEvents += result.writtenEvents;
        malformedEvents += result.malformedEvents;
        rejectedEvents += result.calibrationRejectedEvents;
    }
    const double seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start).count();
    std::cout << "Processed " << inputEvents << " raw events and wrote "
              << writtenEvents << " analysis events in " << seconds
              << " s (" << (seconds > 0.0 ? inputEvents / seconds : 0.0)
              << " raw events/s).\n";
    if (diagnosticsEnabled_) {
        std::cout << "Malformed events skipped: " << malformedEvents << "\n"
                  << "Calibration-rejected events skipped: "
                  << rejectedEvents << "\n";
    }
    return 0;
}
