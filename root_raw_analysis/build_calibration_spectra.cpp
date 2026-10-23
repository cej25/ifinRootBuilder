#include "CalibrationSpectrumBuilder.h"

#include <TH1.h>
#include <TROOT.h>

#include <exception>
#include <iostream>
#include <string>

namespace {

void printUsage(const char* program)
{
    std::cerr
        << "Usage: " << program
        << " OUTPUT.root --group NAME 'INPUT_PATTERN' "
        << "[--group NAME 'INPUT_PATTERN' ...] [--no-ge-time-gate]\n\n"
        << "Quote wildcard patterns so ROOT expands them inside each named "
        << "group. Repeating a group name adds another pattern to that group.\n";
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

    const std::string outputFileName = argv[1];
    CalibrationSpectrumBuilder builder;

    try {
        for (int argument = 2; argument < argc; ++argument) {
            const std::string option = argv[argument];
            if (option == "--group") {
                if (argument + 2 >= argc) {
                    throw std::runtime_error(
                        "--group requires NAME and INPUT_PATTERN");
                }
                const std::string name = argv[++argument];
                const std::string pattern = argv[++argument];
                builder.addGroup(name, pattern);
            } else if (option == "--no-ge-time-gate") {
                builder.setGermaniumTimeGateEnabled(false);
            } else {
                throw std::runtime_error("unknown option '" + option + "'");
            }
        }
    } catch (const std::exception& error) {
        std::cerr << "Argument error: " << error.what() << "\n";
        printUsage(argv[0]);
        return 1;
    }

    return builder.run(outputFileName);
}
