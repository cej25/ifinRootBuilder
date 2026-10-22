#include "DriftSpectrumBuilder.h"

#include <TH1.h>
#include <TROOT.h>

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

void printUsage(const char* program)
{
    std::cerr
        << "Usage: " << program << " OUTPUT.root --mcal MAIN.mcal "
           "--block SERIES FIRST_RUN LAST_RUN 'PATTERN' [OPTIONS]\n\n"
        << "Blocks are processed in command-line order. Options:\n"
        << "  --cal FILE / --mcal FILE       Main calibration stages\n"
        << "  --block SERIES FIRST LAST PAT  Chronological input block\n"
        << "  --group-size N                 Files per group (default 10)\n"
        << "  --threads N                    Group workers (default 4)\n"
        << "  --break-before SERIES RUN      Force a group boundary\n"
        << "  --no-ge-time-gate              Disable the normal Ge time gate\n";
}

} // namespace

int main(int argc, char** argv)
{
    gROOT->SetBatch(true);
    ROOT::EnableThreadSafety();
    TH1::AddDirectory(false);
    if (argc < 8) {
        printUsage(argv[0]);
        return 1;
    }

    DriftSpectrumBuilder builder;
    const std::string outputFileName = argv[1];
    try {
        for (int argument = 2; argument < argc; ++argument) {
            const std::string option = argv[argument];
            if (option == "--cal" || option == "--mcal") {
                if (++argument >= argc) {
                    throw std::runtime_error(option + " requires FILE");
                }
                if (option == "--cal") builder.addCalFile(argv[argument]);
                else builder.addMcalFile(argv[argument]);
            } else if (option == "--block") {
                if (argument + 4 >= argc) {
                    throw std::runtime_error(
                        "--block requires SERIES FIRST_RUN LAST_RUN PATTERN");
                }
                const std::string series = argv[++argument];
                const unsigned int first =
                    parseUnsigned(argv[++argument], "first run");
                const unsigned int last =
                    parseUnsigned(argv[++argument], "last run");
                builder.addBlock(series, first, last, argv[++argument]);
            } else if (option == "--group-size") {
                if (++argument >= argc) {
                    throw std::runtime_error("--group-size requires N");
                }
                builder.setGroupSize(
                    parseUnsigned(argv[argument], "group size"));
            } else if (option == "--threads") {
                if (++argument >= argc) {
                    throw std::runtime_error("--threads requires N");
                }
                builder.setThreadCount(
                    parseUnsigned(argv[argument], "thread count"));
            } else if (option == "--break-before") {
                if (argument + 2 >= argc) {
                    throw std::runtime_error(
                        "--break-before requires SERIES RUN");
                }
                const std::string series = argv[++argument];
                builder.addBreakBefore(
                    series, parseUnsigned(argv[++argument], "break run"));
            } else if (option == "--no-ge-time-gate") {
                builder.setGermaniumTimeGateEnabled(false);
            } else {
                throw std::runtime_error("unknown option '" + option + "'");
            }
        }
    } catch (const std::exception& error) {
        std::cerr << "Argument/configuration error: " << error.what() << "\n";
        printUsage(argv[0]);
        return 1;
    }
    return builder.run(outputFileName);
}
