#include "CalibrationSpectrumBuilder.h"

#include "AnalysisConfig.h"

#include <TChain.h>
#include <TChainElement.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TNamed.h>
#include <TObjArray.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include <RtypesCore.h>

#include <array>
#include <cctype>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

bool validGroupName(const std::string& name)
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

std::string detectorHistogramName(unsigned int detectorID)
{
    std::ostringstream name;
    name << "h1_Ge_ADC_ID" << std::setfill('0') << std::setw(2)
         << detectorID;
    return name.str();
}

std::string detectorHistogramTitle(unsigned int detectorID)
{
    return "Germanium ID " + std::to_string(detectorID) +
           " raw calibration spectrum;ADC channel;Counts";
}

} // namespace

void CalibrationSpectrumBuilder::addGroup(
    const std::string& name, const std::string& inputPattern)
{
    if (!validGroupName(name)) {
        throw std::runtime_error(
            "group names may contain only letters, digits, '_' and '-'");
    }
    if (inputPattern.empty()) {
        throw std::runtime_error("calibration input pattern must not be empty");
    }
    findOrAddGroup(name).inputPatterns.push_back(inputPattern);
}

void CalibrationSpectrumBuilder::setGermaniumTimeGateEnabled(bool enabled)
{
    germaniumTimeGateEnabled_ = enabled;
}

CalibrationSpectrumBuilder::Group&
CalibrationSpectrumBuilder::findOrAddGroup(const std::string& name)
{
    for (Group& group : groups_) {
        if (group.name == name) return group;
    }
    groups_.push_back({name, {}});
    return groups_.back();
}

int CalibrationSpectrumBuilder::run(const std::string& outputFileName) const
{
    if (groups_.empty()) {
        std::cerr << "Error: no calibration source groups were supplied.\n";
        return 2;
    }

    TFile outputFile(outputFileName.c_str(), "RECREATE");
    if (outputFile.IsZombie()) {
        std::cerr << "Error: could not create '" << outputFileName << "'.\n";
        return 3;
    }

    std::string metadataGroup;
    ULong64_t metadataInputFiles = 0;
    ULong64_t metadataInputEvents = 0;
    ULong64_t metadataMalformedEvents = 0;
    ULong64_t metadataGermaniumHits = 0;
    ULong64_t metadataAcceptedGermaniumHits = 0;
    ULong64_t metadataInvalidGermaniumIDs = 0;

    TTree groupMetadata("CalibrationGroups",
                        "Calibration source-group processing summary");
    groupMetadata.Branch("group", &metadataGroup);
    groupMetadata.Branch("inputFiles", &metadataInputFiles);
    groupMetadata.Branch("inputEvents", &metadataInputEvents);
    groupMetadata.Branch("malformedEvents", &metadataMalformedEvents);
    groupMetadata.Branch("germaniumHits", &metadataGermaniumHits);
    groupMetadata.Branch("acceptedGermaniumHits",
                         &metadataAcceptedGermaniumHits);
    groupMetadata.Branch("invalidGermaniumIDs",
                         &metadataInvalidGermaniumIDs);

    std::string provenanceGroup;
    std::string provenanceFile;
    TTree provenance("CalibrationInputFiles",
                     "Resolved ROOT files used for each calibration group");
    provenance.Branch("group", &provenanceGroup);
    provenance.Branch("file", &provenanceFile);

    for (const Group& group : groups_) {
        TChain chain(config::kTreeName);
        for (const std::string& pattern : group.inputPatterns) {
            if (chain.Add(pattern.c_str()) == 0) {
                std::cerr << "Warning: group '" << group.name
                          << "' matched no files for '" << pattern << "'.\n";
            }
        }
        if (chain.GetNtrees() == 0) {
            std::cerr << "Error: group '" << group.name
                      << "' contains no trees called '" << config::kTreeName
                      << "'.\n";
            return 4;
        }

        provenanceGroup = group.name;
        TObjArray* files = chain.GetListOfFiles();
        for (int index = 0; index < files->GetEntries(); ++index) {
            const auto* element =
                dynamic_cast<const TChainElement*>(files->At(index));
            if (element != nullptr) {
                provenanceFile = element->GetTitle();
                provenance.Fill();
            }
        }

        std::array<std::unique_ptr<TH1D>, config::kGermaniumIdBins>
            individual;
        auto totalAdc = std::make_unique<TH1D>(
            "h1_Ge_ADC",
            "All germanium raw calibration spectrum;ADC channel;Counts",
            config::kIndividualGermaniumRawBins,
            -0.5,
            static_cast<double>(config::kIndividualGermaniumRawBins) - 0.5);
        totalAdc->SetDirectory(nullptr);
        totalAdc->SetOption("HIST");
        for (unsigned int detectorID = 0;
             detectorID < individual.size(); ++detectorID) {
            const std::string name = detectorHistogramName(detectorID);
            const std::string title = detectorHistogramTitle(detectorID);
            individual[detectorID] = std::make_unique<TH1D>(
                name.c_str(), title.c_str(),
                config::kIndividualGermaniumRawBins,
                -0.5,
                static_cast<double>(config::kIndividualGermaniumRawBins) - 0.5);
            individual[detectorID]->SetDirectory(nullptr);
            individual[detectorID]->SetOption("HIST");
        }

        auto idVsAdc = std::make_unique<TH2D>(
            "h2_Ge_IDvADC",
            "Germanium detector ID versus raw ADC channel;ADC channel;Detector ID",
            config::kIndividualGermaniumRawBins,
            -0.5,
            static_cast<double>(config::kIndividualGermaniumRawBins) - 0.5,
            config::kGermaniumIdBins,
            -0.5,
            static_cast<double>(config::kGermaniumIdBins) - 0.5);
        idVsAdc->SetDirectory(nullptr);
        idVsAdc->SetOption("COLZ");

        TTreeReader reader(&chain);
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

        ULong64_t inputEvents = 0;
        ULong64_t malformedEvents = 0;
        ULong64_t germaniumHits = 0;
        ULong64_t acceptedGermaniumHits = 0;
        ULong64_t invalidGermaniumIDs = 0;

        std::cout << "Building calibration group " << group.name << " from "
                  << chain.GetNtrees() << " file(s), " << chain.GetEntries()
                  << " events.\n";

        while (reader.Next()) {
            ++inputEvents;
            const std::size_t hitCount = detectorType->size();
            if (detectorID->size() != hitCount ||
                energy->size() != hitCount ||
                relativeNsTime->size() != hitCount ||
                relativePsTime->size() != hitCount) {
                ++malformedEvents;
                continue;
            }

            for (std::size_t hit = 0; hit < hitCount; ++hit) {
                if (detectorType->at(hit) != config::kGermaniumType) continue;
                ++germaniumHits;

                const double relativeTimeNs =
                    static_cast<double>(relativeNsTime->at(hit)) +
                    static_cast<double>(relativePsTime->at(hit)) / 1000.0;
                if (germaniumTimeGateEnabled_ &&
                    !insideInclusive(relativeTimeNs,
                                     config::kGermaniumTimeMinNs,
                                     config::kGermaniumTimeMaxNs)) {
                    continue;
                }

                const unsigned int id = detectorID->at(hit);
                if (id >= individual.size()) {
                    ++invalidGermaniumIDs;
                    continue;
                }

                const double adc = energy->at(hit);
                totalAdc->Fill(adc);
                individual[id]->Fill(adc);
                idVsAdc->Fill(adc, id);
                ++acceptedGermaniumHits;
            }

            if (inputEvents % 1000000ULL == 0ULL) {
                std::cout << "  " << inputEvents << " / "
                          << chain.GetEntries() << " events\r" << std::flush;
            }
        }
        std::cout << "  processed " << inputEvents << " events; accepted "
                  << acceptedGermaniumHits << " in-time Ge hits.          \n";

        outputFile.cd();
        TDirectory* groupDirectory = outputFile.mkdir(group.name.c_str());
        TDirectory* rawDirectory = groupDirectory->mkdir("Raw");
        rawDirectory->cd();
        totalAdc->Write();
        idVsAdc->Write();
        TDirectory* individualDirectory = rawDirectory->mkdir("Individual");
        individualDirectory->cd();
        for (const auto& histogram : individual) histogram->Write();

        metadataGroup = group.name;
        metadataInputFiles = static_cast<ULong64_t>(chain.GetNtrees());
        metadataInputEvents = inputEvents;
        metadataMalformedEvents = malformedEvents;
        metadataGermaniumHits = germaniumHits;
        metadataAcceptedGermaniumHits = acceptedGermaniumHits;
        metadataInvalidGermaniumIDs = invalidGermaniumIDs;
        groupMetadata.Fill();
    }

    outputFile.cd();
    groupMetadata.Write();
    provenance.Write();
    const std::string settingsText = germaniumTimeGateEnabled_
        ? "Ge timing gate enabled: [" +
              std::to_string(config::kGermaniumTimeMinNs) + ", " +
              std::to_string(config::kGermaniumTimeMaxNs) + "] ns"
        : "Ge timing gate disabled";
    TNamed settings("CalibrationSpectrumSettings", settingsText.c_str());
    settings.Write();
    outputFile.Close();

    std::cout << "Wrote calibration spectra to " << outputFileName << ".\n";
    return 0;
}
