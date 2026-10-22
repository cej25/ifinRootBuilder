#include "CalibrationPeakFitter.h"

#include <TH1.h>
#include <TROOT.h>

#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void printUsage(const char* program)
{
    std::cerr
        << "Usage: " << program
        << " INPUT.root OUTPUT.root PEAKS.txt --seed-cal FILE "
        << "[--seed-mcal FILE] [--adc-overrides FILE] "
           "[--seed-cal FILE ...]\n\n"
        << "Seed stages are applied from left to right and are used only to "
        << "locate raw-ADC fit windows.\n";
}

} // namespace

int main(int argc, char** argv)
{
    gROOT->SetBatch(true);
    TH1::AddDirectory(false);

    if (argc < 6) {
        printUsage(argv[0]);
        return 1;
    }

    const std::string inputFileName = argv[1];
    const std::string outputFileName = argv[2];
    const std::string peakFileName = argv[3];
    CalibrationPeakFitter fitter;

    try {
        fitter.loadPeakDefinitions(peakFileName);
        for (int argument = 4; argument < argc; ++argument) {
            const std::string option = argv[argument];
            if (option == "--seed-cal") {
                if (++argument >= argc) {
                    throw std::runtime_error("--seed-cal requires FILE");
                }
                fitter.addSeedCalFile(argv[argument]);
            } else if (option == "--seed-mcal") {
                if (++argument >= argc) {
                    throw std::runtime_error("--seed-mcal requires FILE");
                }
                fitter.addSeedMcalFile(argv[argument]);
            } else if (option == "--adc-overrides") {
                if (++argument >= argc) {
                    throw std::runtime_error("--adc-overrides requires FILE");
                }
                fitter.loadAdcOverrides(argv[argument]);
            } else {
                throw std::runtime_error("unknown option '" + option + "'");
            }
        }
    } catch (const std::exception& error) {
        std::cerr << "Argument/configuration error: " << error.what() << "\n";
        printUsage(argv[0]);
        return 1;
    }

    return fitter.run(inputFileName, outputFileName);
}
