#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TKey.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <set>
#include <string>

namespace {

bool closeEnough(double first, double second, double tolerance)
{
    return std::abs(first - second) <=
        tolerance * std::max({1.0, std::abs(first), std::abs(second)});
}

void compareDirectory(TDirectory& first, TDirectory& second,
                      const std::string& path, double tolerance,
                      unsigned long long& compared,
                      unsigned long long& differences)
{
    std::set<std::string> firstNames;
    TIter next(first.GetListOfKeys());
    while (TKey* key = static_cast<TKey*>(next())) {
        const std::string name = key->GetName();
        firstNames.insert(name);
        TObject* firstObject = first.Get(name.c_str());
        TObject* secondObject = second.Get(name.c_str());
        const std::string objectPath = path + "/" + name;
        if (secondObject == nullptr) {
            std::cout << "Missing from second file: " << objectPath << '\n';
            ++differences;
            continue;
        }

        auto* firstDirectory = dynamic_cast<TDirectory*>(firstObject);
        auto* secondDirectory = dynamic_cast<TDirectory*>(secondObject);
        if (firstDirectory != nullptr || secondDirectory != nullptr) {
            if (firstDirectory == nullptr || secondDirectory == nullptr) {
                std::cout << "Object type differs: " << objectPath << '\n';
                ++differences;
            } else {
                compareDirectory(*firstDirectory, *secondDirectory,
                                 objectPath, tolerance,
                                 compared, differences);
            }
            continue;
        }

        auto* firstHistogram = dynamic_cast<TH1*>(firstObject);
        auto* secondHistogram = dynamic_cast<TH1*>(secondObject);
        if (firstHistogram == nullptr && secondHistogram == nullptr) {
            continue; // Canvases and other presentation objects are skipped.
        }
        if (firstHistogram == nullptr || secondHistogram == nullptr ||
            firstHistogram->GetNcells() != secondHistogram->GetNcells()) {
            std::cout << "Histogram shape/type differs: "
                      << objectPath << '\n';
            ++differences;
            continue;
        }

        ++compared;
        bool different = false;
        for (int bin = 0; bin < firstHistogram->GetNcells(); ++bin) {
            if (!closeEnough(firstHistogram->GetBinContent(bin),
                             secondHistogram->GetBinContent(bin), tolerance) ||
                !closeEnough(firstHistogram->GetBinError(bin),
                             secondHistogram->GetBinError(bin), tolerance)) {
                different = true;
                break;
            }
        }
        if (different) {
            std::cout << "Histogram differs: " << objectPath << '\n';
            ++differences;
        }
    }

    TIter nextSecond(second.GetListOfKeys());
    while (TKey* key = static_cast<TKey*>(nextSecond())) {
        if (firstNames.count(key->GetName()) == 0) {
            std::cout << "Only in second file: " << path << "/"
                      << key->GetName() << '\n';
            ++differences;
        }
    }
}

} // namespace

int compare_histograms(const char* firstFileName,
                       const char* secondFileName,
                       double relativeTolerance = 1.0e-9)
{
    TFile first(firstFileName, "READ");
    TFile second(secondFileName, "READ");
    if (first.IsZombie() || second.IsZombie()) {
        std::cerr << "Could not open both input files.\n";
        return 2;
    }
    unsigned long long compared = 0;
    unsigned long long differences = 0;
    compareDirectory(first, second, "", relativeTolerance,
                     compared, differences);
    std::cout << "Compared " << compared << " histograms: "
              << differences << " difference(s).\n";
    return differences == 0 ? 0 : 1;
}
