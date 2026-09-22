#include "AnalysisTreeAnalysis.h"

#include <TH1.h>
#include <TROOT.h>

#include <exception>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

unsigned int parseThreads(const std::string& text)
{
    std::size_t consumed = 0;
    const unsigned long value = std::stoul(text, &consumed);
    if (consumed != text.size() || value == 0 ||
        value > std::numeric_limits<unsigned int>::max()) {
        throw std::runtime_error("invalid thread count '" + text + "'");
    }
    return static_cast<unsigned int>(value);
}

void usage(const char* program)
{
    std::cerr << "Usage: " << program
              << " OUTPUT.root [--threads N] "
              << "[--diagnostics|--no-diagnostics] "
              << "[--progress|--no-progress] "
              << "[--gates coincidence_gates.txt] ANALYSIS.root [...]\n";
}

} // namespace

int main(int argc, char** argv)
{
    gROOT->SetBatch(true);
    TH1::AddDirectory(false);
    if (argc < 3) {
        usage(argv[0]);
        return 1;
    }
    AnalysisTreeAnalysis analysis;
    const std::string output = argv[1];
    std::vector<std::string> inputs;
    try {
        for (int argument = 2; argument < argc; ++argument) {
            const std::string value = argv[argument];
            if (value == "--threads") {
                if (++argument >= argc) throw std::runtime_error(
                    "--threads requires N");
                analysis.setThreadCount(parseThreads(argv[argument]));
            } else if (value == "--diagnostics" || value == "--no-diagnostics") {
                analysis.setDiagnosticsEnabled(value == "--diagnostics");
            } else if (value == "--progress" || value == "--no-progress") {
                analysis.setProgressEnabled(value == "--progress");
            } else if (value == "--gates") {
                if (++argument >= argc) throw std::runtime_error(
                    "--gates requires a filename");
                analysis.loadCoincidenceGates(argv[argument]);
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
    return analysis.run(inputs, output);
}
