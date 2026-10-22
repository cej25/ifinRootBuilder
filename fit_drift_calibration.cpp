#include "DriftCalibrationFitter.h"

#include <TH1.h>
#include <TROOT.h>

#include <cmath>
#include <exception>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

unsigned int parseUnsigned(const std::string& text,
                           const std::string& description)
{
    std::size_t used = 0;
    const unsigned long value = std::stoul(text, &used);
    if (used != text.size() || value > std::numeric_limits<unsigned int>::max()) {
        throw std::runtime_error("invalid " + description + " '" + text + "'");
    }
    return static_cast<unsigned int>(value);
}

double parseDouble(const std::string& text, const std::string& description)
{
    std::size_t used = 0;
    const double value = std::stod(text, &used);
    if (used != text.size() || !std::isfinite(value)) {
        throw std::runtime_error("invalid " + description + " '" + text + "'");
    }
    return value;
}

void printUsage(const char* program)
{
    std::cerr
        << "Usage: " << program
        << " DRIFT_SPECTRA.root OUTPUT_PREFIX DIAGNOSTICS.root "
           "DRIFT_PEAKS.txt [OPTIONS]\n\n"
        << "Writes OUTPUT_PREFIX_SERIES.cal for every input series and, when "
           "detector 4 anchors are available, "
           "OUTPUT_PREFIX_SERIES_det4.mcal. Options:\n"
        << "  --minimum-peaks N       Minimum accepted anchors (default 2)\n"
        << "  --minimum-span KEV      Minimum anchor span (default 400)\n";
}

} // namespace

int main(int argc, char** argv)
{
    gROOT->SetBatch(true);
    TH1::AddDirectory(false);
    if (argc < 5) {
        printUsage(argv[0]);
        return 1;
    }
    DriftCalibrationFitter fitter;
    try {
        fitter.loadPeakDefinitions(argv[4]);
        for (int argument = 5; argument < argc; ++argument) {
            const std::string option = argv[argument];
            if (option == "--minimum-peaks") {
                if (++argument >= argc) {
                    throw std::runtime_error("--minimum-peaks requires N");
                }
                fitter.setMinimumPeakCount(
                    parseUnsigned(argv[argument], "minimum peak count"));
            } else if (option == "--minimum-span") {
                if (++argument >= argc) {
                    throw std::runtime_error("--minimum-span requires KEV");
                }
                fitter.setMinimumEnergySpan(
                    parseDouble(argv[argument], "minimum energy span"));
            } else {
                throw std::runtime_error("unknown option '" + option + "'");
            }
        }
    } catch (const std::exception& error) {
        std::cerr << "Argument/configuration error: " << error.what() << "\n";
        printUsage(argv[0]);
        return 1;
    }
    return fitter.run(argv[1], argv[2], argv[3]);
}
