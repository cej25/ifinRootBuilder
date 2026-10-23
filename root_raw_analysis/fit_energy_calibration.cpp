#include "EnergyCalibrationFitter.h"

#include <cstdlib>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void printUsage(const char* program)
{
    std::cerr
        << "Usage: " << program
        << " PEAK_FITS.root OUTPUT.mcal DIAGNOSTICS.root"
           " [--config calibration_models.txt] [--max-energy KEV]\n";
}

double parsePositiveDouble(const std::string& text,
                           const std::string& option)
{
    std::size_t used = 0;
    const double value = std::stod(text, &used);
    if (used != text.size() || value <= 0.0) {
        throw std::runtime_error(option + " requires a positive number");
    }
    return value;
}

} // namespace

int main(int argc, char** argv)
{
    if (argc < 4) {
        printUsage(argv[0]);
        return 1;
    }

    try {
        EnergyCalibrationFitter fitter;
        for (int argument = 4; argument < argc; ++argument) {
            const std::string option = argv[argument];
            if (option == "--config" && argument + 1 < argc) {
                fitter.loadModelOverrides(argv[++argument]);
            } else if (option == "--max-energy" && argument + 1 < argc) {
                fitter.setMaximumEnergyKeV(parsePositiveDouble(
                    argv[++argument], option));
            } else {
                throw std::runtime_error("Unknown or incomplete option '" +
                                         option + "'");
            }
        }
        return fitter.run(argv[1], argv[2], argv[3]);
    } catch (const std::exception& error) {
        std::cerr << "Calibration-fit error: " << error.what() << '\n';
        return 2;
    }
}
