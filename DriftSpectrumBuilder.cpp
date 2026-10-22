#include "DriftSpectrumBuilder.h"

#include "AnalysisConfig.h"
#include "RunCalibrationManager.h"

#include <TChain.h>
#include <TChainElement.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TNamed.h>
#include <TObjArray.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include <RtypesCore.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cctype>
#include <condition_variable>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <mutex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <tuple>
#include <vector>

namespace {

constexpr int kDriftBins = 8192;
constexpr double kDriftMinimumKeV = 0.0;
constexpr double kDriftMaximumKeV = 2048.0;

bool validSeriesName(const std::string& name)
{
    if (name.empty()) return false;
    for (const unsigned char character : name) {
        if (!std::isalnum(character) && character != '_' && character != '-') {
            return false;
        }
    }
    return true;
}

bool insideInclusive(double value, double minimum, double maximum)
{
    return value >= minimum && value <= maximum;
}

std::string groupName(unsigned int index)
{
    std::ostringstream name;
    name << "Group" << std::setfill('0') << std::setw(6) << index;
    return name.str();
}

std::string detectorHistogramName(unsigned int detectorID)
{
    std::ostringstream name;
    name << "h1_Ge_E_ID" << std::setfill('0') << std::setw(2) << detectorID;
    return name.str();
}

struct InputFile {
    std::string series;
    unsigned int run = 0;
    std::string fileName;
};

struct SpectrumGroup {
    unsigned int index = 0;
    std::vector<InputFile> files;
};

struct GroupResult {
    std::array<std::unique_ptr<TH1D>, config::kGermaniumIdBins> spectra;
    ULong64_t inputEvents = 0;
    ULong64_t acceptedHits = 0;
};

GroupResult processGroup(const SpectrumGroup& group,
                         const GermaniumCalibration& calibration,
                         bool germaniumTimeGateEnabled)
{
    GroupResult result;
    for (unsigned int detectorID = 0;
         detectorID < result.spectra.size(); ++detectorID) {
        const std::string name = detectorHistogramName(detectorID);
        const std::string title =
            "Germanium ID " + std::to_string(detectorID) +
            " calibrated drift spectrum;Energy (keV);Counts / 0.25 keV";
        result.spectra[detectorID] = std::make_unique<TH1D>(
            name.c_str(), title.c_str(), kDriftBins,
            kDriftMinimumKeV, kDriftMaximumKeV);
        result.spectra[detectorID]->SetDirectory(nullptr);
        result.spectra[detectorID]->SetOption("HIST");
    }

    for (const InputFile& input : group.files) {
        TFile inputFile(input.fileName.c_str(), "READ");
        if (inputFile.IsZombie()) {
            throw std::runtime_error(
                "could not open input file '" + input.fileName + "'");
        }
        TTree* tree = dynamic_cast<TTree*>(inputFile.Get(config::kTreeName));
        if (tree == nullptr) {
            throw std::runtime_error(
                "input file has no tree called '" +
                std::string(config::kTreeName) + "': " + input.fileName);
        }

        TTreeReader reader(tree);
        TTreeReaderValue<std::vector<UShort_t>> detectorID(
            reader, config::kDetectorIDBranch);
        TTreeReaderValue<std::vector<UShort_t>> detectorType(
            reader, config::kDetectorTypeBranch);
        TTreeReaderValue<std::vector<UShort_t>> energy(
            reader, config::kEnergyBranch);
        TTreeReaderValue<std::vector<UShort_t>> relativeNsTime(
            reader, config::kRelativeNsBranch);
        TTreeReaderValue<std::vector<UShort_t>> relativePsTime(
            reader, config::kRelativePsBranch);

        while (reader.Next()) {
            ++result.inputEvents;
            const std::size_t hitCount = detectorType->size();
            if (detectorID->size() != hitCount ||
                energy->size() != hitCount ||
                relativeNsTime->size() != hitCount ||
                relativePsTime->size() != hitCount) {
                continue;
            }
            for (std::size_t hit = 0; hit < hitCount; ++hit) {
                if (detectorType->at(hit) != config::kGermaniumType) continue;
                const unsigned int id = detectorID->at(hit);
                if (id >= result.spectra.size()) continue;
                const double relativeTimeNs =
                    static_cast<double>(relativeNsTime->at(hit)) +
                    static_cast<double>(relativePsTime->at(hit))/1000.0;
                if (germaniumTimeGateEnabled &&
                    !insideInclusive(relativeTimeNs,
                        config::kGermaniumTimeMinNs,
                        config::kGermaniumTimeMaxNs)) {
                    continue;
                }
                const auto calibrated = calibration.calibrate(
                    static_cast<UShort_t>(id), energy->at(hit));
                if (!calibrated.valid()) continue;
                result.spectra[id]->Fill(calibrated.energy);
                ++result.acceptedHits;
            }
        }
    }
    return result;
}

} // namespace

void DriftSpectrumBuilder::addBlock(const std::string& series,
                                    unsigned int firstRun,
                                    unsigned int lastRun,
                                    const std::string& inputPattern)
{
    if (!validSeriesName(series)) {
        throw std::runtime_error(
            "series names may contain only letters, digits, '_' and '-'");
    }
    if (firstRun > lastRun) {
        throw std::runtime_error("block start run must not exceed its end run");
    }
    if (inputPattern.empty()) {
        throw std::runtime_error("block input pattern must not be empty");
    }
    blocks_.push_back({series, firstRun, lastRun, inputPattern});
}

void DriftSpectrumBuilder::addCalFile(const std::string& fileName)
{
    calibration_.addCalFile(fileName);
}

void DriftSpectrumBuilder::addMcalFile(const std::string& fileName)
{
    calibration_.addMcalFile(fileName);
}

void DriftSpectrumBuilder::addBreakBefore(const std::string& series,
                                          unsigned int run)
{
    if (!validSeriesName(series)) {
        throw std::runtime_error("invalid series name in --break-before");
    }
    breakPoints_.push_back({series, run});
}

void DriftSpectrumBuilder::setGroupSize(unsigned int groupSize)
{
    if (groupSize == 0) {
        throw std::runtime_error("group size must be a positive integer");
    }
    groupSize_ = groupSize;
}

void DriftSpectrumBuilder::setThreadCount(unsigned int threadCount)
{
    if (threadCount == 0) {
        throw std::runtime_error("thread count must be a positive integer");
    }
    threadCount_ = threadCount;
}

void DriftSpectrumBuilder::setGermaniumTimeGateEnabled(bool enabled)
{
    germaniumTimeGateEnabled_ = enabled;
}

int DriftSpectrumBuilder::run(const std::string& outputFileName) const
{
    if (blocks_.empty()) {
        std::cerr << "Error: no chronological input blocks were supplied.\n";
        return 2;
    }
    if (calibration_.empty()) {
        std::cerr << "Error: a main calibration is required to build drift "
                     "spectra.\n";
        return 3;
    }

    std::vector<InputFile> inputFiles;
    std::set<std::pair<std::string, unsigned int>> seenRuns;
    std::set<std::pair<std::string, unsigned int>> blockStarts;
    for (const Block& block : blocks_) {
        TChain resolver(config::kTreeName);
        if (resolver.Add(block.inputPattern.c_str()) == 0) {
            std::cerr << "Warning: block " << block.series
                      << " matched no files for '" << block.inputPattern
                      << "'.\n";
            continue;
        }

        std::vector<InputFile> resolved;
        TObjArray* files = resolver.GetListOfFiles();
        for (int index = 0; index < files->GetEntries(); ++index) {
            const auto* element =
                dynamic_cast<const TChainElement*>(files->At(index));
            if (element == nullptr) continue;
            const std::string fileName = element->GetTitle();
            const unsigned int run =
                RunCalibrationManager::runNumberFromFileName(fileName);
            if (run < block.firstRun || run > block.lastRun) continue;
            resolved.push_back({block.series, run, fileName});
        }
        std::sort(resolved.begin(), resolved.end(),
                  [](const InputFile& left, const InputFile& right) {
                      return std::tie(left.run, left.fileName) <
                             std::tie(right.run, right.fileName);
                  });
        if (!resolved.empty()) {
            blockStarts.emplace(resolved.front().series, resolved.front().run);
        }
        for (InputFile& file : resolved) {
            if (!seenRuns.emplace(file.series, file.run).second) {
                throw std::runtime_error(
                    "duplicate series/run in chronological blocks: " +
                    file.series + " run " + std::to_string(file.run));
            }
            inputFiles.push_back(std::move(file));
        }
    }
    if (inputFiles.empty()) {
        std::cerr << "Error: the configured blocks contain no input files.\n";
        return 4;
    }

    const auto isBreakBefore = [this](const InputFile& file) {
        return std::any_of(
            breakPoints_.begin(), breakPoints_.end(),
            [&file](const BreakPoint& point) {
                return point.series == file.series && point.run == file.run;
            });
    };

    std::vector<SpectrumGroup> groups;
    for (const InputFile& file : inputFiles) {
        const bool startNew = groups.empty() ||
            groups.back().files.size() >= groupSize_ || isBreakBefore(file) ||
            blockStarts.count({file.series, file.run}) != 0U;
        if (startNew) {
            groups.push_back(
                {static_cast<unsigned int>(groups.size()), {}});
        }
        groups.back().files.push_back(file);
    }

    TFile outputFile(outputFileName.c_str(), "RECREATE");
    if (outputFile.IsZombie()) {
        std::cerr << "Error: could not create '" << outputFileName << "'.\n";
        return 5;
    }

    unsigned int metadataGroupIndex = 0;
    unsigned int metadataFileCount = 0;
    ULong64_t metadataInputEvents = 0;
    ULong64_t metadataAcceptedHits = 0;
    TTree groupTree("DriftGroups", "Chronological drift-spectrum groups");
    groupTree.Branch("groupIndex", &metadataGroupIndex);
    groupTree.Branch("fileCount", &metadataFileCount);
    groupTree.Branch("inputEvents", &metadataInputEvents);
    groupTree.Branch("acceptedGermaniumHits", &metadataAcceptedHits);

    unsigned int provenanceGroupIndex = 0;
    std::string provenanceSeries;
    unsigned int provenanceRun = 0;
    std::string provenanceFile;
    TTree provenance("DriftInputFiles", "Files in chronological drift groups");
    provenance.Branch("groupIndex", &provenanceGroupIndex);
    provenance.Branch("series", &provenanceSeries);
    provenance.Branch("run", &provenanceRun);
    provenance.Branch("file", &provenanceFile);

    for (const SpectrumGroup& group : groups) {
        for (const InputFile& input : group.files) {
            provenanceGroupIndex = group.index;
            provenanceSeries = input.series;
            provenanceRun = input.run;
            provenanceFile = input.fileName;
            provenance.Fill();
        }
    }

    ULong64_t totalEvents = 0;
    ULong64_t totalAcceptedHits = 0;
    const unsigned int workerCount = std::min<unsigned int>(
        threadCount_, static_cast<unsigned int>(groups.size()));
    std::cout << "Building " << groups.size() << " drift groups from "
              << inputFiles.size() << " files with " << workerCount
              << " worker thread(s).\n";
    std::atomic<std::size_t> nextGroup{0};
    std::atomic<unsigned int> completedGroups{0};
    std::atomic<bool> failed{false};
    std::mutex resultMutex;
    std::condition_variable resultReady;
    std::vector<std::unique_ptr<GroupResult>> results(groups.size());
    std::vector<std::string> errors;
    std::vector<std::thread> workers;
    workers.reserve(workerCount);
    for (unsigned int worker = 0; worker < workerCount; ++worker) {
        workers.emplace_back([&]() {
            while (!failed.load()) {
                const std::size_t index = nextGroup.fetch_add(1);
                if (index >= groups.size()) break;
                try {
                    auto result = std::make_unique<GroupResult>(
                        processGroup(groups[index], calibration_,
                                     germaniumTimeGateEnabled_));
                    {
                        std::lock_guard<std::mutex> lock(resultMutex);
                        results[index] = std::move(result);
                        const unsigned int completed = ++completedGroups;
                        std::cout << "Completed drift group " << index << " ("
                                  << completed << " / " << groups.size()
                                  << ").\n";
                    }
                    resultReady.notify_all();
                } catch (const std::exception& error) {
                    failed = true;
                    {
                        std::lock_guard<std::mutex> lock(resultMutex);
                        errors.push_back(
                            "group " + std::to_string(index) + ": " +
                            error.what());
                    }
                    resultReady.notify_all();
                }
            }
        });
    }

    for (std::size_t index = 0; index < groups.size(); ++index) {
        std::unique_ptr<GroupResult> result;
        {
            std::unique_lock<std::mutex> lock(resultMutex);
            resultReady.wait(lock, [&]() {
                return results[index] != nullptr || failed.load();
            });
            if (results[index] == nullptr) break;
            result = std::move(results[index]);
        }
        const SpectrumGroup& group = groups[index];
        outputFile.cd();
        TDirectory* directory = outputFile.mkdir(groupName(group.index).c_str());
        TDirectory* individual = directory->mkdir("Individual");
        individual->cd();
        for (const auto& spectrum : result->spectra) spectrum->Write();

        metadataGroupIndex = group.index;
        metadataFileCount = static_cast<unsigned int>(group.files.size());
        metadataInputEvents = result->inputEvents;
        metadataAcceptedHits = result->acceptedHits;
        groupTree.Fill();
        totalEvents += result->inputEvents;
        totalAcceptedHits += result->acceptedHits;
    }
    for (std::thread& worker : workers) worker.join();
    if (!errors.empty()) {
        for (const std::string& error : errors) {
            std::cerr << "Error: " << error << "\n";
        }
        return 6;
    }

    outputFile.cd();
    groupTree.Write();
    provenance.Write();
    std::ostringstream settingsText;
    settingsText << "groupSize=" << groupSize_
                 << "; energyBins=" << kDriftBins
                 << "; energyRangeKeV=[" << kDriftMinimumKeV << ","
                 << kDriftMaximumKeV << "]"
                 << "; geTimeGate="
                 << (germaniumTimeGateEnabled_ ? "enabled" : "disabled");
    TNamed settings("DriftSpectrumSettings", settingsText.str().c_str());
    settings.Write();
    outputFile.Close();

    std::cout << "Built " << groups.size() << " drift groups from "
              << inputFiles.size() << " files and " << totalEvents
              << " events; accepted " << totalAcceptedHits
              << " calibrated Ge hits.\nWrote " << outputFileName << ".\n";
    return 0;
}
