#include "GermaniumCalibration.h"

#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace {

std::runtime_error parseError(const std::string& fileName,
                              std::size_t lineNumber,
                              const std::string& message)
{
    return std::runtime_error(
        fileName + ":" + std::to_string(lineNumber) + ": " + message);
}

bool prepareDataLine(std::string& line)
{
    const std::size_t comment = line.find('#');
    if (comment != std::string::npos) {
        line.erase(comment);
    }

    return line.find_first_not_of(" \t\r\n") != std::string::npos;
}

UShort_t readDetectorID(std::istringstream& input,
                         const std::string& fileName,
                         std::size_t lineNumber)
{
    long long detectorID = -1;
    if (!(input >> detectorID) || detectorID < 0 ||
        detectorID > std::numeric_limits<UShort_t>::max()) {
        throw parseError(fileName, lineNumber, "invalid detector ID");
    }
    return static_cast<UShort_t>(detectorID);
}

std::size_t readPositiveCount(std::istringstream& input,
                              const std::string& fileName,
                              std::size_t lineNumber,
                              const std::string& description)
{
    long long count = 0;
    if (!(input >> count) || count <= 0) {
        throw parseError(fileName, lineNumber,
                         description + " must be a positive integer");
    }
    return static_cast<std::size_t>(count);
}

std::vector<double> readCoefficients(std::istringstream& input,
                                     std::size_t count,
                                     const std::string& fileName,
                                     std::size_t lineNumber)
{
    std::vector<double> coefficients;
    coefficients.reserve(count);
    for (std::size_t index = 0; index < count; ++index) {
        double coefficient = 0.0;
        if (!(input >> coefficient) || !std::isfinite(coefficient)) {
            throw parseError(fileName, lineNumber,
                             "missing or invalid polynomial coefficient");
        }
        coefficients.push_back(coefficient);
    }
    return coefficients;
}

void requireEndOfLine(std::istringstream& input,
                      const std::string& fileName,
                      std::size_t lineNumber)
{
    std::string extra;
    if (input >> extra) {
        throw parseError(fileName, lineNumber,
                         "unexpected extra value '" + extra + "'");
    }
}

bool lineMatchesLayout(const std::string& line,
                       bool piecewise,
                       bool hasDetectorGroup)
{
    std::istringstream input(line);
    long long value = 0;
    if (hasDetectorGroup && !(input >> value)) {
        return false;
    }

    long long detectorID = -1;
    long long itemCount = 0;
    if (!(input >> detectorID >> itemCount) || detectorID < 0 ||
        detectorID > std::numeric_limits<UShort_t>::max() ||
        itemCount <= 0) {
        return false;
    }

    const std::size_t repetitions = piecewise
        ? static_cast<std::size_t>(itemCount) : 1U;
    for (std::size_t item = 0; item < repetitions; ++item) {
        long long coefficientCount = itemCount;
        if (piecewise && (!(input >> coefficientCount) ||
                          coefficientCount <= 0)) {
            return false;
        }
        for (long long coefficient = 0;
             coefficient < coefficientCount; ++coefficient) {
            double number = 0.0;
            if (!(input >> number) || !std::isfinite(number)) {
                return false;
            }
        }
        if (piecewise) {
            double endValue = 0.0;
            if (!(input >> endValue) || !std::isfinite(endValue)) {
                return false;
            }
        }
    }

    std::string extra;
    return !(input >> extra);
}

} // namespace

void GermaniumCalibration::addCalFile(const std::string& fileName)
{
    addFile(fileName, false);
}

void GermaniumCalibration::addMcalFile(const std::string& fileName)
{
    addFile(fileName, true);
}

bool GermaniumCalibration::empty() const
{
    return stages_.empty();
}

std::size_t GermaniumCalibration::numberOfStages() const
{
    return stages_.size();
}

double GermaniumCalibration::evaluatePolynomial(
    const std::vector<double>& coefficients,
    double input)
{
    // Horner evaluation for a0 + a1*x + a2*x^2 + ...
    double result = 0.0;
    for (auto coefficient = coefficients.rbegin();
         coefficient != coefficients.rend(); ++coefficient) {
        result = result * input + *coefficient;
    }
    return result;
}

GermaniumCalibration::Result GermaniumCalibration::calibrate(
    UShort_t detectorID,
    double inputEnergy) const
{
    double energy = inputEnergy;

    for (const Stage& stage : stages_) {
        const auto detector = stage.detectors.find(detectorID);
        if (detector == stage.detectors.end()) {
            return {energy, Failure::MissingID};
        }

        const DetectorCalibration& calibration = detector->second;
        const Piece* selectedPiece = nullptr;

        if (!calibration.piecewise) {
            selectedPiece = &calibration.pieces.front();
        } else {
            if (energy < 0.0) {
                return {energy, Failure::OutOfRange};
            }
            for (const Piece& piece : calibration.pieces) {
                if (energy <= piece.endValue) {
                    selectedPiece = &piece;
                    break;
                }
            }
            if (selectedPiece == nullptr) {
                return {energy, Failure::OutOfRange};
            }
        }

        energy = evaluatePolynomial(selectedPiece->coefficients, energy);
        if (!std::isfinite(energy)) {
            return {energy, Failure::NonFiniteResult};
        }
    }

    return {energy, Failure::None};
}

void GermaniumCalibration::addFile(const std::string& fileName,
                                    bool piecewise)
{
    std::ifstream file(fileName);
    if (!file) {
        throw std::runtime_error(
            "Could not open calibration file '" + fileName + "'");
    }

    Stage stage;
    stage.sourceFile = fileName;

    std::string line;
    std::size_t lineNumber = 0;
    bool layoutKnown = false;
    bool hasDetectorGroup = false;
    while (std::getline(file, line)) {
        ++lineNumber;
        if (!prepareDataLine(line)) {
            continue;
        }

        if (!layoutKnown) {
            const bool idFirst = lineMatchesLayout(line, piecewise, false);
            const bool groupThenId = lineMatchesLayout(line, piecewise, true);
            if (idFirst == groupThenId) {
                throw parseError(
                    fileName, lineNumber,
                    idFirst ? "ambiguous calibration-column layout"
                             : "invalid calibration-column layout");
            }
            hasDetectorGroup = groupThenId;
            layoutKnown = true;
        }

        if (!lineMatchesLayout(line, piecewise, hasDetectorGroup)) {
            throw parseError(fileName, lineNumber,
                             "line does not match the file's column layout");
        }

        std::istringstream input(line);
        if (hasDetectorGroup) {
            long long detectorGroup = 0;
            input >> detectorGroup;
        }
        const UShort_t detectorID =
            readDetectorID(input, fileName, lineNumber);

        if (stage.detectors.find(detectorID) != stage.detectors.end()) {
            throw parseError(fileName, lineNumber,
                             "duplicate detector ID " +
                             std::to_string(detectorID));
        }

        DetectorCalibration detectorCalibration;
        detectorCalibration.piecewise = piecewise;

        if (!piecewise) {
            const std::size_t coefficientCount = readPositiveCount(
                input, fileName, lineNumber, "coefficient count");
            detectorCalibration.pieces.push_back({
                readCoefficients(input, coefficientCount,
                                 fileName, lineNumber),
                std::numeric_limits<double>::infinity()
            });
        } else {
            const std::size_t pieceCount = readPositiveCount(
                input, fileName, lineNumber, "piece count");
            double previousEnd = 0.0;

            for (std::size_t pieceIndex = 0;
                 pieceIndex < pieceCount; ++pieceIndex) {
                const std::size_t coefficientCount = readPositiveCount(
                    input, fileName, lineNumber, "coefficient count");
                std::vector<double> coefficients = readCoefficients(
                    input, coefficientCount, fileName, lineNumber);

                double endValue = 0.0;
                if (!(input >> endValue) || !std::isfinite(endValue) ||
                    endValue <= previousEnd) {
                    throw parseError(
                        fileName, lineNumber,
                        "piece end values must be finite and strictly increasing");
                }

                detectorCalibration.pieces.push_back(
                    {std::move(coefficients), endValue});
                previousEnd = endValue;
            }
        }

        requireEndOfLine(input, fileName, lineNumber);
        stage.detectors.emplace(detectorID,
                                std::move(detectorCalibration));
    }

    if (stage.detectors.empty()) {
        throw std::runtime_error(
            "Calibration file '" + fileName + "' contains no detector data");
    }

    stages_.push_back(std::move(stage));
}
