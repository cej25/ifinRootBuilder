#include "CoincidenceGateConfig.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <unordered_set>

namespace {

enum class Section {
    None,
    Symmetric,
    AllVsForward,
    AllVsBackward,
    Double
};

std::string trim(const std::string& value)
{
    const auto first = std::find_if_not(
        value.begin(), value.end(),
        [](unsigned char character) { return std::isspace(character); });
    const auto last = std::find_if_not(
        value.rbegin(), value.rend(),
        [](unsigned char character) { return std::isspace(character); }).base();
    return first < last ? std::string(first, last) : std::string();
}

bool validName(const std::string& name)
{
    return !name.empty() && std::all_of(
        name.begin(), name.end(), [](unsigned char character) {
            return std::isalnum(character) || character == '_' ||
                   character == '-';
        });
}

void validate(const CoincidenceGateDefinition& gate,
              const std::string& fileName, std::size_t lineNumber)
{
    const double values[] = {
        gate.promptMinimum, gate.promptMaximumExclusive,
        gate.lowerMinimum, gate.lowerMaximumExclusive,
        gate.upperMinimum, gate.upperMaximumExclusive};
    if (!std::all_of(std::begin(values), std::end(values),
                     [](double value) { return std::isfinite(value); }) ||
        gate.promptMinimum >= gate.promptMaximumExclusive ||
        gate.lowerMinimum >= gate.lowerMaximumExclusive ||
        gate.upperMinimum >= gate.upperMaximumExclusive) {
        throw std::runtime_error(
            fileName + ":" + std::to_string(lineNumber) +
            ": every gate and sideband requires MIN < MAX");
    }
    if (!validName(gate.name)) {
        throw std::runtime_error(
            fileName + ":" + std::to_string(lineNumber) +
            ": gate names may contain only letters, numbers, '_' and '-'");
    }
}

} // namespace

double CoincidenceGateDefinition::sidebandScale() const
{
    const double promptWidth = promptMaximumExclusive - promptMinimum;
    const double sidebandWidth =
        (lowerMaximumExclusive - lowerMinimum) +
        (upperMaximumExclusive - upperMinimum);
    return promptWidth / sidebandWidth;
}

bool CoincidenceGateDefinition::operator==(
    const CoincidenceGateDefinition& other) const
{
    return name == other.name &&
        promptMinimum == other.promptMinimum &&
        promptMaximumExclusive == other.promptMaximumExclusive &&
        lowerMinimum == other.lowerMinimum &&
        lowerMaximumExclusive == other.lowerMaximumExclusive &&
        upperMinimum == other.upperMinimum &&
        upperMaximumExclusive == other.upperMaximumExclusive;
}

bool DoubleCoincidenceGateDefinition::operator==(
    const DoubleCoincidenceGateDefinition& other) const
{
    return name == other.name &&
        requiredMinimum == other.requiredMinimum &&
        requiredMaximumExclusive == other.requiredMaximumExclusive &&
        secondGate == other.secondGate;
}

CoincidenceGateConfig CoincidenceGateConfig::load(
    const std::string& fileName)
{
    std::ifstream input(fileName);
    if (!input) {
        throw std::runtime_error(
            "Could not open coincidence gate file '" + fileName + "'");
    }

    CoincidenceGateConfig result;
    Section section = Section::None;
    std::unordered_set<std::string> namesBySection[4];
    std::string rawLine;
    std::size_t lineNumber = 0;
    while (std::getline(input, rawLine)) {
        ++lineNumber;
        const std::size_t comment = rawLine.find('#');
        const std::string line = trim(rawLine.substr(0, comment));
        if (line.empty()) continue;

        if (line.front() == '[' && line.back() == ']') {
            const std::string heading = trim(line.substr(1, line.size() - 2));
            if (heading == "Symmetric") {
                section = Section::Symmetric;
            } else if (heading == "AllvFW") {
                section = Section::AllVsForward;
            } else if (heading == "AllvBW") {
                section = Section::AllVsBackward;
            } else if (heading == "Double") {
                section = Section::Double;
            } else {
                throw std::runtime_error(
                    fileName + ":" + std::to_string(lineNumber) +
                    ": unknown section '[" + heading + "]'");
            }
            continue;
        }
        if (section == Section::None) {
            throw std::runtime_error(
                fileName + ":" + std::to_string(lineNumber) +
                ": gate appears before a section heading");
        }

        std::istringstream fields(line);
        CoincidenceGateDefinition gate;
        DoubleCoincidenceGateDefinition doubleGate;
        if (section == Section::Double) {
            if (!(fields >> doubleGate.name >> doubleGate.requiredMinimum
                  >> doubleGate.requiredMaximumExclusive
                  >> doubleGate.secondGate.promptMinimum
                  >> doubleGate.secondGate.promptMaximumExclusive
                  >> doubleGate.secondGate.lowerMinimum
                  >> doubleGate.secondGate.lowerMaximumExclusive
                  >> doubleGate.secondGate.upperMinimum
                  >> doubleGate.secondGate.upperMaximumExclusive)) {
                throw std::runtime_error(
                    fileName + ":" + std::to_string(lineNumber) +
                    ": expected NAME REQUIRED_MIN REQUIRED_MAX PROMPT_MIN "
                    "PROMPT_MAX LOWER_MIN LOWER_MAX UPPER_MIN UPPER_MAX");
            }
            doubleGate.secondGate.name = doubleGate.name;
            gate = doubleGate.secondGate;
        } else {
            if (!(fields >> gate.name >> gate.promptMinimum
                  >> gate.promptMaximumExclusive >> gate.lowerMinimum
                  >> gate.lowerMaximumExclusive >> gate.upperMinimum
                  >> gate.upperMaximumExclusive)) {
                throw std::runtime_error(
                    fileName + ":" + std::to_string(lineNumber) +
                    ": expected NAME PROMPT_MIN PROMPT_MAX LOWER_MIN "
                    "LOWER_MAX UPPER_MIN UPPER_MAX");
            }
        }
        std::string extra;
        if (fields >> extra) {
            throw std::runtime_error(
                fileName + ":" + std::to_string(lineNumber) +
                ": unexpected extra field '" + extra + "'");
        }
        validate(gate, fileName, lineNumber);
        if (section == Section::Double &&
            (!std::isfinite(doubleGate.requiredMinimum) ||
             !std::isfinite(doubleGate.requiredMaximumExclusive) ||
             doubleGate.requiredMinimum >=
                 doubleGate.requiredMaximumExclusive)) {
            throw std::runtime_error(
                fileName + ":" + std::to_string(lineNumber) +
                ": required double-gate window needs MIN < MAX");
        }

        const std::size_t sectionIndex =
            section == Section::Symmetric ? 0U :
            section == Section::AllVsForward ? 1U :
            section == Section::AllVsBackward ? 2U : 3U;
        if (!namesBySection[sectionIndex].insert(gate.name).second) {
            throw std::runtime_error(
                fileName + ":" + std::to_string(lineNumber) +
                ": duplicate gate name '" + gate.name + "' in section");
        }
        if (section == Section::Symmetric) {
            result.symmetric_.push_back(gate);
        } else if (section == Section::AllVsForward) {
            result.allVsForward_.push_back(gate);
        } else if (section == Section::AllVsBackward) {
            result.allVsBackward_.push_back(gate);
        } else {
            result.doubleGates_.push_back(doubleGate);
        }
    }
    return result;
}

const std::vector<CoincidenceGateDefinition>&
CoincidenceGateConfig::symmetric() const
{
    return symmetric_;
}

const std::vector<CoincidenceGateDefinition>&
CoincidenceGateConfig::allVsForward() const
{
    return allVsForward_;
}

const std::vector<CoincidenceGateDefinition>&
CoincidenceGateConfig::allVsBackward() const
{
    return allVsBackward_;
}

const std::vector<DoubleCoincidenceGateDefinition>&
CoincidenceGateConfig::doubleGates() const
{
    return doubleGates_;
}
