#include "EnergyCalibrationFitter.h"

#include "AnalysisConfig.h"

#include <TDirectory.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TNamed.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

struct CalibrationPoint {
    unsigned int detectorID = 0;
    std::string group;
    std::string label;
    double adc = 0.0;
    double adcError = 0.0;
    double energyKeV = 0.0;
};

struct Override {
    double boundaryEnergyKeV = 1000.0;
    unsigned int lowerDegree = 4;
    unsigned int upperDegree = 4;
};

struct Piece {
    std::vector<double> coefficients;
    double endAdc = 0.0;
};

struct Model {
    unsigned int detectorID = 0;
    std::vector<Piece> pieces;
    std::vector<bool> used;
    std::vector<double> residuals;
    double rmsKeV = 0.0;
    double maximumResidualKeV = 0.0;
    double bic = 0.0;
    std::string description;
};

bool prepareDataLine(std::string& line)
{
    const std::size_t comment = line.find('#');
    if (comment != std::string::npos) line.erase(comment);
    return line.find_first_not_of(" \t\r\n") != std::string::npos;
}

double median(std::vector<double> values)
{
    if (values.empty()) return 0.0;
    const std::size_t middle = values.size()/2;
    std::nth_element(values.begin(), values.begin()+middle, values.end());
    const double upper = values[middle];
    if (values.size()%2 != 0) return upper;
    std::nth_element(values.begin(), values.begin()+middle-1, values.end());
    return 0.5*(upper + values[middle-1]);
}

std::vector<double> solveLinearSystem(std::vector<std::vector<double>> matrix,
                                      std::vector<double> values)
{
    const std::size_t size = values.size();
    for (std::size_t column = 0; column < size; ++column) {
        std::size_t pivot = column;
        for (std::size_t row = column+1; row < size; ++row) {
            if (std::abs(matrix[row][column]) >
                std::abs(matrix[pivot][column])) {
                pivot = row;
            }
        }
        if (std::abs(matrix[pivot][column]) < 1.0e-14) {
            throw std::runtime_error("singular calibration fit");
        }
        std::swap(matrix[column], matrix[pivot]);
        std::swap(values[column], values[pivot]);

        const double diagonal = matrix[column][column];
        for (std::size_t item = column; item < size; ++item) {
            matrix[column][item] /= diagonal;
        }
        values[column] /= diagonal;

        for (std::size_t row = 0; row < size; ++row) {
            if (row == column) continue;
            const double factor = matrix[row][column];
            for (std::size_t item = column; item < size; ++item) {
                matrix[row][item] -= factor*matrix[column][item];
            }
            values[row] -= factor*values[column];
        }
    }
    return values;
}

double binomial(unsigned int n, unsigned int k)
{
    if (k > n) return 0.0;
    double result = 1.0;
    for (unsigned int index = 1; index <= k; ++index) {
        result *= static_cast<double>(n + 1 - index) /
                  static_cast<double>(index);
    }
    return result;
}

std::vector<double> normalizedToRaw(const std::vector<double>& normalized,
                                    double centre,
                                    double scale)
{
    std::vector<double> raw(normalized.size(), 0.0);
    for (unsigned int power = 0; power < normalized.size(); ++power) {
        const double factor = normalized[power] /
            std::pow(scale, static_cast<int>(power));
        for (unsigned int rawPower = 0; rawPower <= power; ++rawPower) {
            raw[rawPower] += factor*binomial(power, rawPower)*
                std::pow(-centre, static_cast<int>(power-rawPower));
        }
    }
    return raw;
}

double evaluatePolynomial(const std::vector<double>& coefficients, double adc)
{
    double value = 0.0;
    for (auto coefficient = coefficients.rbegin();
         coefficient != coefficients.rend(); ++coefficient) {
        value = value*adc + *coefficient;
    }
    return value;
}

std::vector<double> fitPolynomial(const std::vector<CalibrationPoint>& points,
                                  const std::vector<bool>& used,
                                  unsigned int degree)
{
    double minimum = std::numeric_limits<double>::infinity();
    double maximum = -std::numeric_limits<double>::infinity();
    for (std::size_t index = 0; index < points.size(); ++index) {
        if (!used[index]) continue;
        minimum = std::min(minimum, points[index].adc);
        maximum = std::max(maximum, points[index].adc);
    }
    const double centre = 0.5*(minimum + maximum);
    const double scale = 0.5*(maximum - minimum);
    if (!std::isfinite(scale) || scale <= 0.0) {
        throw std::runtime_error("calibration points have no ADC range");
    }

    const std::size_t parameters = degree + 1;
    std::vector<std::vector<double>> normal(
        parameters, std::vector<double>(parameters, 0.0));
    std::vector<double> right(parameters, 0.0);
    for (std::size_t index = 0; index < points.size(); ++index) {
        if (!used[index]) continue;
        const double coordinate = (points[index].adc-centre)/scale;
        std::vector<double> basis(parameters, 1.0);
        for (std::size_t power = 1; power < parameters; ++power) {
            basis[power] = basis[power-1]*coordinate;
        }
        for (std::size_t row = 0; row < parameters; ++row) {
            right[row] += basis[row]*points[index].energyKeV;
            for (std::size_t column = 0; column < parameters; ++column) {
                normal[row][column] += basis[row]*basis[column];
            }
        }
    }
    return normalizedToRaw(
        solveLinearSystem(std::move(normal), std::move(right)), centre, scale);
}

std::vector<double> residualsForPolynomial(
    const std::vector<CalibrationPoint>& points,
    const std::vector<double>& coefficients)
{
    std::vector<double> residuals;
    residuals.reserve(points.size());
    for (const CalibrationPoint& point : points) {
        residuals.push_back(
            evaluatePolynomial(coefficients, point.adc)-point.energyKeV);
    }
    return residuals;
}

bool rejectWorstOutlier(const std::vector<double>& residuals,
                        std::vector<bool>& used,
                        std::size_t parameterCount)
{
    std::vector<double> activeResiduals;
    for (std::size_t index = 0; index < residuals.size(); ++index) {
        if (used[index]) activeResiduals.push_back(residuals[index]);
    }
    if (activeResiduals.size() <= parameterCount + 3) return false;

    const double centre = median(activeResiduals);
    std::vector<double> deviations;
    deviations.reserve(activeResiduals.size());
    for (const double residual : activeResiduals) {
        deviations.push_back(std::abs(residual-centre));
    }
    const double sigma = 1.4826*median(std::move(deviations));
    const double threshold = std::max(1.0, 3.0*sigma);

    std::size_t worst = residuals.size();
    double worstDeviation = threshold;
    for (std::size_t index = 0; index < residuals.size(); ++index) {
        if (!used[index]) continue;
        const double deviation = std::abs(residuals[index]-centre);
        if (deviation > worstDeviation) {
            worstDeviation = deviation;
            worst = index;
        }
    }
    if (worst == residuals.size()) return false;
    used[worst] = false;
    return true;
}

void finishStatistics(Model& model)
{
    double sumSquares = 0.0;
    std::size_t count = 0;
    model.maximumResidualKeV = 0.0;
    for (std::size_t index = 0; index < model.residuals.size(); ++index) {
        if (!model.used[index]) continue;
        sumSquares += model.residuals[index]*model.residuals[index];
        model.maximumResidualKeV = std::max(
            model.maximumResidualKeV, std::abs(model.residuals[index]));
        ++count;
    }
    model.rmsKeV = count == 0 ? 0.0 :
        std::sqrt(sumSquares/static_cast<double>(count));
}

Model robustPolynomialModel(const std::vector<CalibrationPoint>& points,
                            unsigned int degree)
{
    Model model;
    model.detectorID = points.front().detectorID;
    model.used.assign(points.size(), true);
    std::vector<double> coefficients;
    for (unsigned int iteration = 0; iteration < 8; ++iteration) {
        coefficients = fitPolynomial(points, model.used, degree);
        model.residuals = residualsForPolynomial(points, coefficients);
        if (!rejectWorstOutlier(
                model.residuals, model.used, degree+1)) break;
    }
    coefficients = fitPolynomial(points, model.used, degree);
    model.residuals = residualsForPolynomial(points, coefficients);
    model.pieces.push_back({std::move(coefficients), 0.0});
    finishStatistics(model);

    const std::size_t usedCount = static_cast<std::size_t>(
        std::count(model.used.begin(), model.used.end(), true));
    const std::size_t rejectedCount = points.size()-usedCount;
    const double rss = std::max(
        1.0e-12, model.rmsKeV*model.rmsKeV*usedCount);
    model.bic = usedCount*std::log(rss/usedCount) +
        (degree+1)*std::log(static_cast<double>(usedCount)) +
        2.0*rejectedCount*std::log(static_cast<double>(points.size()));
    model.description = "polynomial degree " + std::to_string(degree);
    return model;
}

double interpolateBoundaryAdc(const std::vector<CalibrationPoint>& points,
                              double boundaryEnergyKeV)
{
    std::vector<const CalibrationPoint*> ordered;
    ordered.reserve(points.size());
    for (const CalibrationPoint& point : points) ordered.push_back(&point);
    std::sort(ordered.begin(), ordered.end(),
              [](const CalibrationPoint* left, const CalibrationPoint* right) {
                  return left->energyKeV < right->energyKeV;
              });
    for (std::size_t index = 1; index < ordered.size(); ++index) {
        const CalibrationPoint& low = *ordered[index-1];
        const CalibrationPoint& high = *ordered[index];
        if (low.energyKeV <= boundaryEnergyKeV &&
            high.energyKeV >= boundaryEnergyKeV &&
            high.energyKeV > low.energyKeV) {
            const double fraction = (boundaryEnergyKeV-low.energyKeV) /
                (high.energyKeV-low.energyKeV);
            return low.adc + fraction*(high.adc-low.adc);
        }
    }
    throw std::runtime_error("piecewise boundary is outside fitted energies");
}

std::pair<std::vector<double>, std::vector<double>> fitContinuousPiecewise(
    const std::vector<CalibrationPoint>& points,
    const std::vector<bool>& used,
    double boundaryAdc,
    unsigned int lowerDegree,
    unsigned int upperDegree)
{
    double scale = 0.0;
    for (std::size_t index = 0; index < points.size(); ++index) {
        if (used[index]) {
            scale = std::max(scale, std::abs(points[index].adc-boundaryAdc));
        }
    }
    if (scale <= 0.0) throw std::runtime_error("invalid piecewise ADC range");

    const std::size_t parameters = 1 + lowerDegree + upperDegree;
    std::vector<std::vector<double>> normal(
        parameters, std::vector<double>(parameters, 0.0));
    std::vector<double> right(parameters, 0.0);
    for (std::size_t index = 0; index < points.size(); ++index) {
        if (!used[index]) continue;
        const double coordinate = (points[index].adc-boundaryAdc)/scale;
        std::vector<double> basis(parameters, 0.0);
        basis[0] = 1.0;
        double powerValue = 1.0;
        if (points[index].adc <= boundaryAdc) {
            for (unsigned int power = 1; power <= lowerDegree; ++power) {
                powerValue *= coordinate;
                basis[power] = powerValue;
            }
        } else {
            for (unsigned int power = 1; power <= upperDegree; ++power) {
                powerValue *= coordinate;
                basis[lowerDegree+power] = powerValue;
            }
        }
        for (std::size_t row = 0; row < parameters; ++row) {
            right[row] += basis[row]*points[index].energyKeV;
            for (std::size_t column = 0; column < parameters; ++column) {
                normal[row][column] += basis[row]*basis[column];
            }
        }
    }
    const std::vector<double> solution = solveLinearSystem(
        std::move(normal), std::move(right));
    std::vector<double> lower(lowerDegree+1, solution[0]);
    std::vector<double> upper(upperDegree+1, solution[0]);
    for (unsigned int power = 1; power <= lowerDegree; ++power) {
        lower[power] = solution[power];
    }
    for (unsigned int power = 1; power <= upperDegree; ++power) {
        upper[power] = solution[lowerDegree+power];
    }
    return {
        normalizedToRaw(lower, boundaryAdc, scale),
        normalizedToRaw(upper, boundaryAdc, scale)
    };
}

double evaluateModel(const Model& model, double adc)
{
    for (const Piece& piece : model.pieces) {
        if (piece.endAdc == 0.0 || adc <= piece.endAdc) {
            return evaluatePolynomial(piece.coefficients, adc);
        }
    }
    return evaluatePolynomial(model.pieces.back().coefficients, adc);
}

std::vector<double> residualsForModel(
    const std::vector<CalibrationPoint>& points,
    const Model& model)
{
    std::vector<double> residuals;
    residuals.reserve(points.size());
    for (const CalibrationPoint& point : points) {
        residuals.push_back(evaluateModel(model, point.adc)-point.energyKeV);
    }
    return residuals;
}

Model robustPiecewiseModel(const std::vector<CalibrationPoint>& points,
                           const Override& specification)
{
    Model model;
    model.detectorID = points.front().detectorID;
    model.used.assign(points.size(), true);
    const double boundaryAdc = interpolateBoundaryAdc(
        points, specification.boundaryEnergyKeV);
    const std::size_t parameterCount =
        1 + specification.lowerDegree + specification.upperDegree;

    for (unsigned int iteration = 0; iteration < 8; ++iteration) {
        auto coefficients = fitContinuousPiecewise(
            points, model.used, boundaryAdc,
            specification.lowerDegree, specification.upperDegree);
        model.pieces = {
            {std::move(coefficients.first), boundaryAdc},
            {std::move(coefficients.second), 0.0}
        };
        model.residuals = residualsForModel(points, model);
        if (!rejectWorstOutlier(
                model.residuals, model.used, parameterCount)) break;
    }
    auto coefficients = fitContinuousPiecewise(
        points, model.used, boundaryAdc,
        specification.lowerDegree, specification.upperDegree);
    model.pieces = {
        {std::move(coefficients.first), boundaryAdc},
        {std::move(coefficients.second), 0.0}
    };
    model.residuals = residualsForModel(points, model);
    finishStatistics(model);
    model.description = "continuous piecewise degrees " +
        std::to_string(specification.lowerDegree) + "/" +
        std::to_string(specification.upperDegree) + " at " +
        std::to_string(specification.boundaryEnergyKeV) + " keV";
    return model;
}

double findMaximumAdc(const Model& model, double maximumEnergyKeV)
{
    double previousEnergy = evaluateModel(model, 0.0);
    for (unsigned int channel = 1;
         channel < config::kIndividualGermaniumRawBins; ++channel) {
        const double energy = evaluateModel(model, static_cast<double>(channel));
        if (energy >= maximumEnergyKeV && previousEnergy < maximumEnergyKeV) {
            double low = static_cast<double>(channel-1);
            double high = static_cast<double>(channel);
            for (unsigned int iteration = 0; iteration < 40; ++iteration) {
                const double middle = 0.5*(low+high);
                if (evaluateModel(model, middle) < maximumEnergyKeV) {
                    low = middle;
                } else {
                    high = middle;
                }
            }
            return 0.5*(low+high);
        }
        previousEnergy = energy;
    }
    throw std::runtime_error("calibration does not reach maximum energy");
}

double highEnergySlope(const std::vector<CalibrationPoint>& points,
                       const std::vector<bool>& used)
{
    std::vector<std::size_t> ordered;
    for (std::size_t index = 0; index < points.size(); ++index) {
        if (used[index]) ordered.push_back(index);
    }
    std::sort(ordered.begin(), ordered.end(),
              [&points](std::size_t left, std::size_t right) {
                  return points[left].energyKeV < points[right].energyKeV;
              });
    if (ordered.size() > 5) ordered.erase(ordered.begin(), ordered.end()-5);

    std::vector<double> slopes;
    for (std::size_t first = 0; first < ordered.size(); ++first) {
        for (std::size_t second = first+1; second < ordered.size(); ++second) {
            const CalibrationPoint& low = points[ordered[first]];
            const CalibrationPoint& high = points[ordered[second]];
            const double adcDifference = high.adc-low.adc;
            if (std::abs(adcDifference) < 10.0) continue;
            const double slope =
                (high.energyKeV-low.energyKeV)/adcDifference;
            if (std::isfinite(slope) && slope > 0.0) slopes.push_back(slope);
        }
    }
    if (slopes.empty()) {
        throw std::runtime_error("could not estimate a high-energy slope");
    }
    return median(std::move(slopes));
}

void setMaximumRange(Model& model,
                     const std::vector<CalibrationPoint>& points,
                     double maximumEnergyKeV)
{
    std::size_t highestUsed = points.size();
    for (std::size_t index = 0; index < points.size(); ++index) {
        if (!model.used[index]) continue;
        if (highestUsed == points.size() ||
            points[index].energyKeV > points[highestUsed].energyKeV) {
            highestUsed = index;
        }
    }
    if (highestUsed == points.size()) {
        throw std::runtime_error("calibration has no retained points");
    }

    if (points[highestUsed].energyKeV >= maximumEnergyKeV) {
        model.pieces.back().endAdc = findMaximumAdc(
            model, maximumEnergyKeV);
        return;
    }

    // Never extrapolate a high-order polynomial beyond the final measured
    // calibration point.  Such a polynomial can turn over even when all of
    // its in-range residuals are good.  Continue instead with a short linear
    // tail whose positive slope is the median of the highest-energy secants.
    const double joinAdc = points[highestUsed].adc;
    const double joinEnergy = evaluateModel(model, joinAdc);
    const double slope = highEnergySlope(points, model.used);
    if (joinEnergy >= maximumEnergyKeV) {
        model.pieces.back().endAdc = findMaximumAdc(
            model, maximumEnergyKeV);
        return;
    }
    const double finalAdc = joinAdc +
        (maximumEnergyKeV-joinEnergy)/slope;
    if (!std::isfinite(finalAdc) || finalAdc <= joinAdc ||
        finalAdc >= config::kIndividualGermaniumRawBins) {
        throw std::runtime_error("invalid high-energy linear extension");
    }
    model.pieces.back().endAdc = joinAdc;
    model.pieces.push_back({
        {joinEnergy-slope*joinAdc, slope}, finalAdc
    });
    model.description += " plus linear high-energy extension";
}

std::map<unsigned int, Override> loadOverrides(const std::string& fileName)
{
    std::map<unsigned int, Override> overrides;
    if (fileName.empty()) {
        overrides.emplace(4U, Override{});
        return overrides;
    }
    std::ifstream file(fileName);
    if (!file) throw std::runtime_error(
        "Could not open model configuration '" + fileName + "'");
    std::string line;
    std::size_t lineNumber = 0;
    while (std::getline(file, line)) {
        ++lineNumber;
        if (!prepareDataLine(line)) continue;
        unsigned int detectorID = 0;
        std::string modelType;
        Override specification;
        std::istringstream input(line);
        if (!(input >> detectorID >> modelType >>
              specification.boundaryEnergyKeV >>
              specification.lowerDegree >> specification.upperDegree) ||
            modelType != "piecewise" ||
            specification.boundaryEnergyKeV <= 0.0 ||
            specification.lowerDegree < 1 || specification.lowerDegree > 6 ||
            specification.upperDegree < 1 || specification.upperDegree > 6) {
            throw std::runtime_error(fileName + ":" +
                std::to_string(lineNumber) + ": invalid model override");
        }
        std::string extra;
        if (input >> extra || !overrides.emplace(
                detectorID, specification).second) {
            throw std::runtime_error(fileName + ":" +
                std::to_string(lineNumber) +
                ": extra value or duplicate detector ID");
        }
    }
    return overrides;
}

std::string detectorDirectoryName(unsigned int detectorID)
{
    std::ostringstream name;
    name << "Detector_" << std::setfill('0') << std::setw(2) << detectorID;
    return name.str();
}

} // namespace

void EnergyCalibrationFitter::loadModelOverrides(const std::string& fileName)
{
    overrideFileName_ = fileName;
}

void EnergyCalibrationFitter::setMaximumEnergyKeV(double energyKeV)
{
    if (!std::isfinite(energyKeV) || energyKeV <= 0.0) {
        throw std::runtime_error("maximum calibration energy must be positive");
    }
    maximumEnergyKeV_ = energyKeV;
}

int EnergyCalibrationFitter::run(
    const std::string& peakFitFileName,
    const std::string& outputCalibrationFileName,
    const std::string& diagnosticFileName) const
{
    TFile inputFile(peakFitFileName.c_str(), "READ");
    if (inputFile.IsZombie()) {
        throw std::runtime_error("Could not open '" + peakFitFileName + "'");
    }
    TTree* peakTree = dynamic_cast<TTree*>(inputFile.Get("PeakFits"));
    if (peakTree == nullptr) {
        throw std::runtime_error("Input file has no PeakFits tree");
    }

    std::array<std::vector<CalibrationPoint>, config::kGermaniumIdBins> points;
    TTreeReader reader(peakTree);
    TTreeReaderValue<UShort_t> detectorID(reader, "detectorID");
    TTreeReaderValue<std::string> group(reader, "group");
    TTreeReaderValue<std::string> label(reader, "label");
    TTreeReaderValue<double> referenceEnergy(reader, "referenceEnergyKeV");
    TTreeReaderValue<double> centroid(reader, "centroidAdc");
    TTreeReaderValue<double> centroidError(reader, "centroidErrorAdc");
    TTreeReaderValue<bool> accepted(reader, "accepted");
    while (reader.Next()) {
        if (!*accepted || *detectorID >= points.size()) continue;
        points[*detectorID].push_back({
            *detectorID, *group, *label, *centroid, *centroidError,
            *referenceEnergy
        });
    }

    const auto overrides = loadOverrides(overrideFileName_);
    std::vector<Model> models;
    for (unsigned int id = 0; id < points.size(); ++id) {
        if (points[id].empty()) continue;
        if (points[id].size() < 6) {
            throw std::runtime_error("Detector " + std::to_string(id) +
                                     " has too few calibration points");
        }

        Model selected;
        const auto override = overrides.find(id);
        if (override != overrides.end()) {
            selected = robustPiecewiseModel(points[id], override->second);
        } else {
            selected = robustPolynomialModel(points[id], 1);
            for (unsigned int degree = 2; degree <= 3; ++degree) {
                Model candidate = robustPolynomialModel(points[id], degree);
                if (candidate.bic < selected.bic) selected = std::move(candidate);
            }
        }
        setMaximumRange(selected, points[id], maximumEnergyKeV_);
        std::cout << "Detector " << std::setw(2) << id << ": "
                  << selected.description << ", RMS "
                  << std::fixed << std::setprecision(3)
                  << selected.rmsKeV << " keV, used "
                  << std::count(selected.used.begin(), selected.used.end(), true)
                  << "/" << points[id].size() << " points.\n";
        models.push_back(std::move(selected));
    }

    std::ofstream calibrationFile(outputCalibrationFileName);
    if (!calibrationFile) throw std::runtime_error(
        "Could not create '" + outputCalibrationFileName + "'");
    calibrationFile << std::setprecision(17);
    calibrationFile << "# Raw ADC to keV calibration generated from "
                    << peakFitFileName << '\n';
    calibrationFile << "# Valid up to " << maximumEnergyKeV_ << " keV.\n";
    for (const Model& model : models) {
        calibrationFile << model.detectorID << ' ' << model.pieces.size();
        for (const Piece& piece : model.pieces) {
            calibrationFile << ' ' << piece.coefficients.size();
            for (const double coefficient : piece.coefficients) {
                calibrationFile << ' ' << coefficient;
            }
            calibrationFile << ' ' << piece.endAdc;
        }
        calibrationFile << '\n';
    }
    calibrationFile.close();

    TFile diagnosticFile(diagnosticFileName.c_str(), "RECREATE");
    if (diagnosticFile.IsZombie()) throw std::runtime_error(
        "Could not create '" + diagnosticFileName + "'");

    UShort_t outputDetectorID = 0;
    unsigned int pieceIndex = 0;
    unsigned int degree = 0;
    double endAdc = 0.0;
    double rmsKeV = 0.0;
    double maximumResidualKeV = 0.0;
    unsigned int pointsUsed = 0;
    unsigned int pointsRejected = 0;
    std::string modelDescription;
    std::vector<double> coefficients;
    TTree modelTree("CalibrationModels", "Selected calibration models");
    modelTree.Branch("detectorID", &outputDetectorID);
    modelTree.Branch("pieceIndex", &pieceIndex);
    modelTree.Branch("degree", &degree);
    modelTree.Branch("endAdc", &endAdc);
    modelTree.Branch("rmsKeV", &rmsKeV);
    modelTree.Branch("maximumResidualKeV", &maximumResidualKeV);
    modelTree.Branch("pointsUsed", &pointsUsed);
    modelTree.Branch("pointsRejected", &pointsRejected);
    modelTree.Branch("description", &modelDescription);
    modelTree.Branch("coefficients", &coefficients);

    std::string pointGroup;
    std::string pointLabel;
    double pointAdc = 0.0;
    double pointEnergyKeV = 0.0;
    double fittedEnergyKeV = 0.0;
    double residualKeV = 0.0;
    bool pointUsed = false;
    TTree residualTree("CalibrationResiduals", "Calibration point residuals");
    residualTree.Branch("detectorID", &outputDetectorID);
    residualTree.Branch("group", &pointGroup);
    residualTree.Branch("label", &pointLabel);
    residualTree.Branch("adc", &pointAdc);
    residualTree.Branch("referenceEnergyKeV", &pointEnergyKeV);
    residualTree.Branch("fittedEnergyKeV", &fittedEnergyKeV);
    residualTree.Branch("residualKeV", &residualKeV);
    residualTree.Branch("used", &pointUsed);

    for (const Model& model : models) {
        outputDetectorID = static_cast<UShort_t>(model.detectorID);
        rmsKeV = model.rmsKeV;
        maximumResidualKeV = model.maximumResidualKeV;
        pointsUsed = static_cast<unsigned int>(
            std::count(model.used.begin(), model.used.end(), true));
        pointsRejected = static_cast<unsigned int>(
            model.used.size()-pointsUsed);
        modelDescription = model.description;
        for (pieceIndex = 0; pieceIndex < model.pieces.size(); ++pieceIndex) {
            coefficients = model.pieces[pieceIndex].coefficients;
            degree = static_cast<unsigned int>(coefficients.size()-1);
            endAdc = model.pieces[pieceIndex].endAdc;
            modelTree.Fill();
        }

        const std::string directoryName = detectorDirectoryName(model.detectorID);
        TDirectory* directory = diagnosticFile.mkdir(directoryName.c_str());
        directory->cd();
        TGraphErrors calibrationGraph;
        calibrationGraph.SetName("g_CalibrationPoints");
        calibrationGraph.SetTitle(
            ("Germanium ID " + std::to_string(model.detectorID) +
             ";Raw ADC channel;Energy (keV)").c_str());
        TGraphErrors residualGraph;
        residualGraph.SetName("g_CalibrationResiduals");
        residualGraph.SetTitle(
            ("Germanium ID " + std::to_string(model.detectorID) +
             ";Energy (keV);Fit - reference (keV)").c_str());
        TGraphErrors curveGraph;
        curveGraph.SetName("g_CalibrationCurve");
        curveGraph.SetTitle(
            ("Germanium ID " + std::to_string(model.detectorID) +
             ";Raw ADC channel;Fitted energy (keV)").c_str());

        for (std::size_t index = 0; index < points[model.detectorID].size();
             ++index) {
            const CalibrationPoint& point = points[model.detectorID][index];
            calibrationGraph.SetPoint(
                static_cast<int>(index), point.adc, point.energyKeV);
            calibrationGraph.SetPointError(
                static_cast<int>(index), point.adcError, 0.0);
            residualGraph.SetPoint(
                static_cast<int>(index), point.energyKeV,
                model.residuals[index]);
            residualGraph.SetPointError(
                static_cast<int>(index), 0.0, 0.0);

            pointGroup = point.group;
            pointLabel = point.label;
            pointAdc = point.adc;
            pointEnergyKeV = point.energyKeV;
            fittedEnergyKeV = evaluateModel(model, point.adc);
            residualKeV = fittedEnergyKeV-pointEnergyKeV;
            pointUsed = model.used[index];
            residualTree.Fill();
        }
        const double finalAdc = model.pieces.back().endAdc;
        for (int index = 0; index <= 1000; ++index) {
            const double adc = finalAdc*index/1000.0;
            curveGraph.SetPoint(index, adc, evaluateModel(model, adc));
            curveGraph.SetPointError(index, 0.0, 0.0);
        }
        calibrationGraph.Write();
        residualGraph.Write();
        curveGraph.Write();
        diagnosticFile.cd();
    }

    modelTree.Write();
    residualTree.Write();
    TNamed sourceFile("CalibrationPeakFitInput", peakFitFileName.c_str());
    sourceFile.Write();
    TNamed outputRange(
        "CalibrationMaximumEnergyKeV",
        std::to_string(maximumEnergyKeV_).c_str());
    outputRange.Write();
    diagnosticFile.Close();

    std::cout << "Wrote calibration to " << outputCalibrationFileName
              << " and diagnostics to " << diagnosticFileName << ".\n";
    return 0;
}
