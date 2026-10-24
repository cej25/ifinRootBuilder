#include "AnalysisConfig.h"

#include <TAxis.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMultiGraph.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include <RtypesCore.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct CalibrationPoint {
    double adc = 0.0;
    double energy = 0.0;
    double fitted = 0.0;
    double residual = 0.0;
    bool used = false;
};

struct CalibrationPiece {
    unsigned int index = 0;
    double endAdc = 0.0;
    std::vector<double> coefficients;
};

struct ResolutionPoint {
    double energy = 0.0;
    double fwhm = 0.0;
    double error = 0.0;
};

struct DriftPoint {
    unsigned int group = 0;
    unsigned int chronologicalPosition = 0;
    int source = 3;
    double offset = 0.0;
    double offsetError = 0.0;
    double slope = 1.0;
    double slopeError = 0.0;
};

std::string detectorName(unsigned int id)
{
    std::ostringstream name;
    name << "Detector_" << std::setfill('0') << std::setw(2) << id;
    return name.str();
}

void openPdf(TCanvas& canvas, const std::string& path)
{
    canvas.Print((path+"[").c_str());
}

void closePdf(TCanvas& canvas, const std::string& path)
{
    canvas.Print((path+"]").c_str());
}

double derivative(const std::vector<CalibrationPiece>& pieces, double adc)
{
    if (pieces.empty()) return std::numeric_limits<double>::quiet_NaN();
    const CalibrationPiece* selected = &pieces.back();
    for (const CalibrationPiece& piece : pieces) {
        if (adc <= piece.endAdc) {
            selected = &piece;
            break;
        }
    }
    double result = 0.0;
    double power = 1.0;
    for (std::size_t degree = 1;
         degree < selected->coefficients.size(); ++degree) {
        result += static_cast<double>(degree)*
            selected->coefficients[degree]*power;
        power *= adc;
    }
    return result;
}

std::map<unsigned int, std::vector<CalibrationPoint>> loadCalibrationPoints(
    TFile& file)
{
    TTree* tree = dynamic_cast<TTree*>(file.Get("CalibrationResiduals"));
    if (tree == nullptr) {
        throw std::runtime_error("CalibrationResiduals tree is missing");
    }
    std::map<unsigned int, std::vector<CalibrationPoint>> result;
    TTreeReader reader(tree);
    TTreeReaderValue<UShort_t> detector(reader, "detectorID");
    TTreeReaderValue<double> adc(reader, "adc");
    TTreeReaderValue<double> energy(reader, "referenceEnergyKeV");
    TTreeReaderValue<double> fitted(reader, "fittedEnergyKeV");
    TTreeReaderValue<double> residual(reader, "residualKeV");
    TTreeReaderValue<Bool_t> used(reader, "used");
    while (reader.Next()) {
        result[*detector].push_back({
            *adc, *energy, *fitted, *residual, static_cast<bool>(*used)});
    }
    return result;
}

std::map<unsigned int, std::vector<CalibrationPiece>> loadCalibrationModels(
    TFile& file)
{
    TTree* tree = dynamic_cast<TTree*>(file.Get("CalibrationModels"));
    if (tree == nullptr) {
        throw std::runtime_error("CalibrationModels tree is missing");
    }
    std::map<unsigned int, std::vector<CalibrationPiece>> result;
    TTreeReader reader(tree);
    TTreeReaderValue<UShort_t> detector(reader, "detectorID");
    TTreeReaderValue<unsigned int> piece(reader, "pieceIndex");
    TTreeReaderValue<double> endAdc(reader, "endAdc");
    TTreeReaderValue<std::vector<double>> coefficients(reader, "coefficients");
    while (reader.Next()) {
        result[*detector].push_back({*piece, *endAdc, *coefficients});
    }
    for (auto& entry : result) {
        std::sort(entry.second.begin(), entry.second.end(),
            [](const CalibrationPiece& left, const CalibrationPiece& right) {
                return left.index < right.index;
            });
    }
    return result;
}

std::map<unsigned int, std::vector<ResolutionPoint>> loadResolutionPoints(
    TFile& peakFile,
    const std::map<unsigned int, std::vector<CalibrationPiece>>& models)
{
    TTree* tree = dynamic_cast<TTree*>(peakFile.Get("PeakFits"));
    if (tree == nullptr) throw std::runtime_error("PeakFits tree is missing");
    std::map<unsigned int, std::vector<ResolutionPoint>> result;
    TTreeReader reader(tree);
    TTreeReaderValue<UShort_t> detector(reader, "detectorID");
    TTreeReaderValue<double> energy(reader, "referenceEnergyKeV");
    TTreeReaderValue<double> centroid(reader, "centroidAdc");
    TTreeReaderValue<double> sigma(reader, "sigmaAdc");
    TTreeReaderValue<double> sigmaError(reader, "sigmaErrorAdc");
    TTreeReaderValue<Bool_t> accepted(reader, "accepted");
    constexpr double kSigmaToFwhm = 2.354820045;
    while (reader.Next()) {
        if (!static_cast<bool>(*accepted)) continue;
        const auto model = models.find(*detector);
        if (model == models.end()) continue;
        const double gain = std::abs(derivative(model->second, *centroid));
        if (!std::isfinite(gain) || gain <= 0.0) continue;
        result[*detector].push_back({
            *energy, kSigmaToFwhm*(*sigma)*gain,
            kSigmaToFwhm*(*sigmaError)*gain});
    }
    for (auto& entry : result) {
        std::sort(entry.second.begin(), entry.second.end(),
            [](const ResolutionPoint& left, const ResolutionPoint& right) {
                return left.energy < right.energy;
            });
    }
    return result;
}

std::map<unsigned int, std::vector<DriftPoint>> loadDriftPoints(TFile& file)
{
    TTree* tree = dynamic_cast<TTree*>(file.Get("DriftCorrections"));
    if (tree == nullptr) {
        throw std::runtime_error("DriftCorrections tree is missing");
    }
    std::map<unsigned int, std::vector<DriftPoint>> result;
    TTreeReader reader(tree);
    TTreeReaderValue<unsigned int> group(reader, "groupIndex");
    TTreeReaderValue<unsigned int> chronologicalPosition(
        reader, "chronologicalPosition");
    TTreeReaderValue<UShort_t> detector(reader, "detectorID");
    TTreeReaderValue<int> source(reader, "source");
    TTreeReaderValue<Bool_t> valid(reader, "valid");
    TTreeReaderValue<double> offset(reader, "offsetKeV");
    TTreeReaderValue<double> slope(reader, "slope");
    TTreeReaderValue<double> offsetError(reader, "offsetErrorKeV");
    TTreeReaderValue<double> slopeError(reader, "slopeError");
    while (reader.Next()) {
        if (!static_cast<bool>(*valid)) continue;
        result[*detector].push_back({
            *group, *chronologicalPosition, *source, *offset, *offsetError,
            *slope, *slopeError});
    }
    for (auto& entry : result) {
        std::sort(entry.second.begin(), entry.second.end(),
            [](const DriftPoint& left, const DriftPoint& right) {
                return left.chronologicalPosition < right.chronologicalPosition;
            });
    }
    return result;
}

void styleGraph(TGraphErrors& graph, Color_t colour, Style_t marker)
{
    graph.SetLineColor(colour);
    graph.SetMarkerColor(colour);
    graph.SetMarkerStyle(marker);
    graph.SetMarkerSize(0.8);
    graph.SetLineWidth(2);
}

void plotCalibrationPages(
    TFile& diagnosticFile,
    const std::map<unsigned int, std::vector<CalibrationPoint>>& points,
    TFile& output, const std::string& pdf)
{
    TCanvas canvas("cCalibration", "Energy calibration", 1100, 850);
    openPdf(canvas, pdf);
    for (const auto& entry : points) {
        const unsigned int detector = entry.first;
        const std::string path = detectorName(detector)+"/g_CalibrationCurve";
        TGraphErrors* storedCurve = dynamic_cast<TGraphErrors*>(
            diagnosticFile.Get(path.c_str()));
        if (storedCurve == nullptr) continue;
        TGraphErrors curve(*storedCurve), used, rejected, usedResidual,
            rejectedResidual;
        for (const CalibrationPoint& point : entry.second) {
            TGraphErrors& calibration = point.used ? used : rejected;
            int index = calibration.GetN();
            calibration.SetPoint(index, point.adc, point.energy);
            TGraphErrors& residual = point.used ? usedResidual : rejectedResidual;
            index = residual.GetN();
            residual.SetPoint(index, point.energy, point.residual);
        }
        styleGraph(curve, kAzure+2, 1);
        styleGraph(used, kBlack, 20);
        styleGraph(rejected, kRed+1, 24);
        styleGraph(usedResidual, kBlack, 20);
        styleGraph(rejectedResidual, kRed+1, 24);
        curve.SetTitle(("Germanium ID "+std::to_string(detector)+
            ";Raw ADC channel;Energy [keV]").c_str());
        usedResidual.SetTitle(
            ";Reference energy [keV];Fit - reference [keV]");
        canvas.Clear();
        canvas.Divide(1, 2);
        canvas.cd(1); gPad->SetGrid();
        curve.Draw("AL"); used.Draw("P SAME"); rejected.Draw("P SAME");
        TLegend legend(0.14, 0.70, 0.38, 0.88);
        legend.SetBorderSize(0); legend.SetFillStyle(0);
        legend.AddEntry(&curve, "Selected calibration", "l");
        legend.AddEntry(&used, "Used peaks", "p");
        if (rejected.GetN() != 0) legend.AddEntry(&rejected, "Rejected peaks", "p");
        legend.Draw();
        canvas.cd(2); gPad->SetGrid();
        usedResidual.Draw("AP"); rejectedResidual.Draw("P SAME");
        TLine zero(usedResidual.GetXaxis()->GetXmin(), 0.0,
                   usedResidual.GetXaxis()->GetXmax(), 0.0);
        zero.SetLineStyle(2); zero.Draw();
        canvas.Print(pdf.c_str());
        output.cd();
        canvas.Write(("c_Calibration_ID"+
                      (detector < 10 ? std::string("0") : std::string())+
                      std::to_string(detector)).c_str());
    }
    closePdf(canvas, pdf);
}

void plotResolutionPages(
    const std::map<unsigned int, std::vector<ResolutionPoint>>& points,
    TFile& output, const std::string& pdf)
{
    TCanvas canvas("cResolution", "Energy resolution", 1000, 750);
    openPdf(canvas, pdf);
    for (const auto& entry : points) {
        TGraphErrors graph;
        for (const ResolutionPoint& point : entry.second) {
            const int index = graph.GetN();
            graph.SetPoint(index, point.energy, point.fwhm);
            graph.SetPointError(index, 0.0, point.error);
        }
        styleGraph(graph, kBlack, 20);
        graph.SetTitle(("Germanium ID "+std::to_string(entry.first)+
            ";Reference energy [keV];FWHM [keV]").c_str());
        canvas.Clear(); canvas.cd(); gPad->SetGrid();
        graph.Draw("AP");
        std::unique_ptr<TF1> trend;
        if (graph.GetN() >= 4) {
            trend = std::make_unique<TF1>(
                "resolutionTrend",
                "sqrt([0]*[0]+[1]*[1]*x+[2]*[2]*x*x)",
                0.0, config::kDriftEnergyMax);
            trend->SetParameters(0.8, 0.03, 0.0005);
            trend->SetLineColor(kAzure+2);
            trend->SetLineWidth(2);
            graph.Fit(trend.get(), "Q");
            trend->Draw("SAME");
        }
        canvas.Print(pdf.c_str());
        output.cd();
        canvas.Write(("c_Resolution_ID"+
                      (entry.first < 10 ? std::string("0") : std::string())+
                      std::to_string(entry.first)).c_str());
    }
    closePdf(canvas, pdf);
}

void plotDriftPages(
    const std::map<unsigned int, std::vector<DriftPoint>>& points,
    TFile& output, const std::string& pdf)
{
    TCanvas canvas("cDrift", "Run-by-run drift", 1100, 850);
    openPdf(canvas, pdf);
    for (const auto& entry : points) {
        TGraphErrors offset, gain, measuredOffset, measuredGain;
        for (const DriftPoint& point : entry.second) {
            int index = offset.GetN();
            offset.SetPoint(index, point.chronologicalPosition, point.offset);
            offset.SetPointError(index, 0.0, point.offsetError);
            gain.SetPoint(index, point.chronologicalPosition,
                          1000.0*(point.slope-1.0));
            gain.SetPointError(index, 0.0, 1000.0*point.slopeError);
            if (point.source == 0) {
                index = measuredOffset.GetN();
                measuredOffset.SetPoint(
                    index, point.chronologicalPosition, point.offset);
                measuredGain.SetPoint(
                    index, point.chronologicalPosition,
                    1000.0*(point.slope-1.0));
            }
        }
        styleGraph(offset, kAzure+2, 1);
        styleGraph(gain, kAzure+2, 1);
        styleGraph(measuredOffset, kBlack, 20);
        styleGraph(measuredGain, kBlack, 20);
        offset.SetTitle(("Germanium ID "+std::to_string(entry.first)+
            ";Chronological drift group;Offset correction [keV]").c_str());
        gain.SetTitle(
            ";Chronological drift group;Gain correction [per mille]");
        canvas.Clear(); canvas.Divide(1, 2);
        canvas.cd(1); gPad->SetGrid();
        offset.Draw("AL"); measuredOffset.Draw("P SAME");
        TLegend legend(0.14, 0.72, 0.38, 0.88);
        legend.SetBorderSize(0); legend.SetFillStyle(0);
        legend.AddEntry(&offset, "Applied trajectory", "l");
        legend.AddEntry(&measuredOffset, "Directly fitted groups", "p");
        legend.Draw();
        canvas.cd(2); gPad->SetGrid();
        gain.Draw("AL"); measuredGain.Draw("P SAME");
        canvas.Print(pdf.c_str());
        output.cd();
        canvas.Write(("c_Drift_ID"+
                      (entry.first < 10 ? std::string("0") : std::string())+
                      std::to_string(entry.first)).c_str());
    }
    closePdf(canvas, pdf);
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 5) {
        std::cerr
            << "Usage: " << argv[0]
            << " PEAK_FITS.root CALIBRATION_DIAGNOSTICS.root "
               "DRIFT_DIAGNOSTICS.root OUTPUT_DIRECTORY\n";
        return 1;
    }
    try {
        gROOT->SetBatch(true);
        gStyle->SetOptStat(0);
        gStyle->SetOptFit(0);
        gStyle->SetTitleBorderSize(0);
        gStyle->SetLegendBorderSize(0);
        const std::filesystem::path outputDirectory(argv[4]);
        std::filesystem::create_directories(outputDirectory);
        TFile peakFile(argv[1], "READ");
        TFile calibrationFile(argv[2], "READ");
        TFile driftFile(argv[3], "READ");
        if (peakFile.IsZombie() || calibrationFile.IsZombie() ||
            driftFile.IsZombie()) {
            throw std::runtime_error("could not open one or more input files");
        }
        const auto calibrationPoints = loadCalibrationPoints(calibrationFile);
        const auto models = loadCalibrationModels(calibrationFile);
        const auto resolutionPoints = loadResolutionPoints(peakFile, models);
        const auto driftPoints = loadDriftPoints(driftFile);
        const std::filesystem::path rootPath =
            outputDirectory/"calibration_diagnostic_plots.root";
        TFile output(rootPath.string().c_str(), "RECREATE");
        if (output.IsZombie()) {
            throw std::runtime_error("could not create output ROOT file");
        }
        plotCalibrationPages(
            calibrationFile, calibrationPoints, output,
            (outputDirectory/"energy_calibrations.pdf").string());
        plotResolutionPages(
            resolutionPoints, output,
            (outputDirectory/"energy_resolution.pdf").string());
        plotDriftPages(
            driftPoints, output,
            (outputDirectory/"energy_drift.pdf").string());
        output.Close();
        std::cout << "Wrote calibration plots to "
                  << outputDirectory.string() << ".\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Plotting error: " << error.what() << '\n';
        return 2;
    }
}
