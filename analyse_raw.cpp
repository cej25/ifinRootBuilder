#include "RawAnalysis.h"

#include <TROOT.h>

#include <exception>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void printUsage(const char* program)
{
    std::cerr
        << "Usage: " << program
        << " output.root [--cal file.cal] [--mcal file.mcal] "
        << "[--run-cal first last file.cal] "
        << "[--run-mcal first last file.mcal] "
        << "[--diagnostics|--no-diagnostics] "
        << "[--exclude-ge LUT[,LUT...]] "
        << "input1.root [input2.root ...]\n\n"
        << "Examples:\n"
        << "  " << program << " analysis.root Run_30cm_000123.root\n"
        << "  " << program
        << " analysis.root --cal coarse.cal --mcal fine.mcal "
        << "Run_30cm_*.root\n"
        << "  " << program
        << " analysis.root --run-cal 1 10 run1-10.cal "
        << "--run-mcal 1 10 run1-10.mcal "
        << "--run-cal 11 20 run11-20.cal Run_30cm_*.root\n";
}

unsigned int parseRunNumber(const std::string& text)
{
    std::size_t consumed = 0;
    const unsigned long value = std::stoul(text, &consumed);
    if (consumed != text.size() ||
        value > std::numeric_limits<unsigned int>::max()) {
        throw std::runtime_error("invalid run number '" + text + "'");
    }
    return static_cast<unsigned int>(value);
}

void addExcludedGermaniumLUTs(
    const std::string& text, RawAnalysis& analysis)
{
    if (text.empty()) {
        throw std::runtime_error("empty germanium detector exclusion list");
    }

    std::stringstream values(text);
    std::string value;
    while (std::getline(values, value, ',')) {
        if (value.empty()) {
            throw std::runtime_error(
                "invalid germanium detector exclusion list '" + text + "'");
        }
        std::size_t consumed = 0;
        const unsigned long detectorLUT = std::stoul(value, &consumed);
        if (consumed != value.size() ||
            detectorLUT > std::numeric_limits<unsigned int>::max()) {
            throw std::runtime_error(
                "invalid germanium detector LUT '" + value + "'");
        }
        analysis.excludeGermaniumLUT(
            static_cast<unsigned int>(detectorLUT));
    }
}

} // namespace

int main(int argc, char** argv)
{
    gROOT->SetBatch(true);

    if (argc < 3) {
        printUsage(argv[0]);
        return 1;
    }

    const std::string outputFileName = argv[1];
    RawAnalysis analysis;
    std::vector<std::string> inputPatterns;
    inputPatterns.reserve(static_cast<std::size_t>(argc - 2));

    try {
        for (int argument = 2; argument < argc; ++argument) {
            const std::string value = argv[argument];
            if (value == "--diagnostics" || value == "--no-diagnostics") {
                analysis.setDiagnosticsEnabled(value == "--diagnostics");
            } else if (value == "--exclude-ge") {
                if (argument + 1 >= argc) {
                    std::cerr << "Error: --exclude-ge requires a comma-separated "
                              << "list of detector LUTs.\n";
                    return 1;
                }
                addExcludedGermaniumLUTs(argv[++argument], analysis);
            } else if (value == "--cal" || value == "--mcal") {
                if (argument + 1 >= argc) {
                    std::cerr << "Error: " << value
                              << " requires a filename.\n";
                    return 1;
                }
                const std::string calibrationFile = argv[++argument];
                if (value == "--cal") {
                    analysis.addCalFile(calibrationFile);
                } else {
                    analysis.addMcalFile(calibrationFile);
                }
            } else if (value == "--run-cal" || value == "--run-mcal") {
                if (argument + 3 >= argc) {
                    std::cerr << "Error: " << value
                              << " requires FIRST_RUN LAST_RUN FILE.\n";
                    return 1;
                }
                const unsigned int firstRun = parseRunNumber(argv[++argument]);
                const unsigned int lastRun = parseRunNumber(argv[++argument]);
                const std::string calibrationFile = argv[++argument];
                if (value == "--run-cal") {
                    analysis.addRunCalFile(
                        firstRun, lastRun, calibrationFile);
                } else {
                    analysis.addRunMcalFile(
                        firstRun, lastRun, calibrationFile);
                }
            } else if (!value.empty() && value.front() == '-') {
                std::cerr << "Error: unknown option '" << value << "'.\n";
                return 1;
            } else {
                inputPatterns.push_back(value);
            }
        }
    } catch (const std::exception& error) {
        std::cerr << "Argument error: " << error.what() << "\n";
        return 1;
    }

    if (inputPatterns.empty()) {
        printUsage(argv[0]);
        return 1;
    }

    return analysis.run(inputPatterns, outputFileName);
}
