#include "AnalysisTreeBuilder.h"

#include <TROOT.h>

#include <exception>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

unsigned int parseUnsigned(const std::string& text, const char* label)
{
    std::size_t consumed = 0;
    const unsigned long value = std::stoul(text, &consumed);
    if (consumed != text.size() ||
        value > std::numeric_limits<unsigned int>::max()) {
        throw std::runtime_error("invalid " + std::string(label) +
                                 " '" + text + "'");
    }
    return static_cast<unsigned int>(value);
}

void addExcluded(const std::string& text, AnalysisTreeBuilder& builder)
{
    std::stringstream values(text);
    std::string value;
    while (std::getline(values, value, ',')) {
        if (value.empty()) {
            throw std::runtime_error("invalid empty Ge exclusion");
        }
        builder.excludeGermaniumID(parseUnsigned(value, "Ge ID"));
    }
}

void usage(const char* program)
{
    std::cerr << "Usage: " << program
              << " OUTPUT_DIRECTORY [OPTIONS] INPUT.root [...]\n"
              << "Options: --cal FILE --mcal FILE "
              << "--run-cal FIRST LAST FILE --run-mcal FIRST LAST FILE "
              << "--threads N --exclude-ge ID[,ID...] "
              << "--diagnostics|--no-diagnostics\n";
}

} // namespace

int main(int argc, char** argv)
{
    gROOT->SetBatch(true);
    ROOT::EnableThreadSafety();
    if (argc < 3) {
        usage(argv[0]);
        return 1;
    }

    AnalysisTreeBuilder builder;
    const std::string outputDirectory = argv[1];
    std::vector<std::string> inputs;
    try {
        for (int argument = 2; argument < argc; ++argument) {
            const std::string value = argv[argument];
            if (value == "--threads") {
                if (++argument >= argc) throw std::runtime_error(
                    "--threads requires N");
                builder.setThreadCount(parseUnsigned(argv[argument], "thread count"));
            } else if (value == "--diagnostics" || value == "--no-diagnostics") {
                builder.setDiagnosticsEnabled(value == "--diagnostics");
            } else if (value == "--exclude-ge") {
                if (++argument >= argc) throw std::runtime_error(
                    "--exclude-ge requires IDs");
                addExcluded(argv[argument], builder);
            } else if (value == "--cal" || value == "--mcal") {
                if (++argument >= argc) throw std::runtime_error(
                    value + " requires a file");
                if (value == "--cal") builder.addCalFile(argv[argument]);
                else builder.addMcalFile(argv[argument]);
            } else if (value == "--run-cal" || value == "--run-mcal") {
                if (argument + 3 >= argc) throw std::runtime_error(
                    value + " requires FIRST LAST FILE");
                const unsigned int first = parseUnsigned(argv[++argument], "run");
                const unsigned int last = parseUnsigned(argv[++argument], "run");
                const std::string file = argv[++argument];
                if (value == "--run-cal") builder.addRunCalFile(first, last, file);
                else builder.addRunMcalFile(first, last, file);
            } else if (!value.empty() && value.front() == '-') {
                throw std::runtime_error("unknown option '" + value + "'");
            } else {
                inputs.push_back(value);
            }
        }
    } catch (const std::exception& error) {
        std::cerr << "Argument error: " << error.what() << "\n";
        return 1;
    }
    if (inputs.empty()) {
        usage(argv[0]);
        return 1;
    }
    return builder.run(inputs, outputDirectory);
}
