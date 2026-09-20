#ifndef GERMANIUM_CALIBRATION_H
#define GERMANIUM_CALIBRATION_H

#include <RtypesCore.h>

#include <string>
#include <unordered_map>
#include <vector>

class GermaniumCalibration {
public:
    enum class Failure {
        None,
        MissingID,
        OutOfRange,
        NonFiniteResult
    };

    struct Result {
        double energy;
        Failure failure;

        bool valid() const { return failure == Failure::None; }
    };

    void addCalFile(const std::string& fileName);
    void addMcalFile(const std::string& fileName);

    bool empty() const;
    std::size_t numberOfStages() const;
    Result calibrate(UShort_t detectorID, double inputEnergy) const;

private:
    struct Piece {
        std::vector<double> coefficients;
        double endValue;
    };

    struct DetectorCalibration {
        bool piecewise;
        std::vector<Piece> pieces;
    };

    struct Stage {
        std::string sourceFile;
        std::unordered_map<UShort_t, DetectorCalibration> detectors;
    };

    void addFile(const std::string& fileName, bool piecewise);
    static double evaluatePolynomial(
        const std::vector<double>& coefficients,
        double input);

    std::vector<Stage> stages_;
};

#endif
