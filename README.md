# Raw ROOT-tree analysis

## Build

The program requires:

- CMake 3.16 or newer;
- a C++17 compiler;
- a ROOT installation containing Core, RIO, Tree, TreePlayer, Hist, Graf, and
  Gpad, with implicit multithreading support (Imt).

Make sure the ROOT environment is active before configuring. For example,
`root-config --version` should run successfully.

From the project directory:

```bash
cmake -S . -B build
cmake --build build -j
```

The executables are:

```text
build/build_analysis_tree
build/analyse_tree
build/analyse_raw
build/build_calibration_spectra
build/fit_calibration_peaks
build/fit_energy_calibration
build/build_drift_spectra
build/fit_drift_calibration
build/plot_calibration_diagnostics
```

To perform a completely clean rebuild:

```bash
rm -rf build
cmake -S . -B build
cmake --build build -j
```

## Usage

### Build raw calibration spectra

`build_calibration_spectra` combines named groups of source trees without
modifying or physically merging the input files. It writes a 16,384-channel
raw ADC spectrum for each Ge detector, an ID-versus-ADC matrix, input-file
provenance and processing counters. The normal valid-Ge timing gate is applied;
BGO vetoing, silicon conditions and energy calibration are not applied.

Quote every wildcard so ROOT, rather than the shell, expands it as one named
group. For the current experiment's mid-run sources:

```bash
build/build_calibration_spectra mid_source_spectra.root \
    --group Mid_Eu152 \
    '/data4/cj0041/IFIN-HH/102Mo/rootbuild/Calibrations/152Eu_midexp_*.root' \
    --group Mid_Co56 \
    '/data4/cj0041/IFIN-HH/102Mo/rootbuild/Calibrations/56Co_midexp_*.root' \
    --group Mid_Co60 \
    '/data4/cj0041/IFIN-HH/102Mo/rootbuild/Calibrations/60Co_midexp_*.root'
```

The output structure is:

```text
Mid_Eu152/Raw/h2_Ge_IDvADC
Mid_Eu152/Raw/h1_Ge_ADC
Mid_Eu152/Raw/Individual/h1_Ge_ADC_ID00 ... ID24
Mid_Co56/Raw/...
Mid_Co60/Raw/...
CalibrationGroups
CalibrationInputFiles
CalibrationSpectrumSettings
```

Repeating a group name adds another input pattern to that group. The
`--no-ge-time-gate` option is available for diagnostics but should not normally
be used for the final calibration spectra.

### Fit calibration-source peaks

`fit_calibration_peaks` measures source-peak centroids in the raw per-detector
spectra. It writes every attempted fit to the `PeakFits` tree, including the
raw-ADC centroid and uncertainty, width, amplitude, background, chi-square,
degrees of freedom, fit status, and acceptance flag. Accepted measurements are
also written as one `TGraphErrors` per detector under `CalibrationPoints`.

An existing calibration is supplied only to locate each fit window. Its
coefficients neither constrain the fitted centroid nor enter the new
calibration as measurements:

```bash
build/fit_calibration_peaks \
    mid_source_spectra.root \
    mid_source_peak_fits.root \
    calibration_peaks.txt \
    --seed-cal calibrations/152Eu.cal \
    --seed-mcal calibrations/pint_run_1-10.mcal \
    --adc-overrides calibration_peak_adc_overrides.txt
```

The editable `calibration_peaks.txt` columns are:

```text
GROUP  ENERGY_KEV  LABEL  SEARCH_HALF_WIDTH_KEV  FIT_HALF_WIDTH_KEV
```

Comment out a line to omit it or adjust its windows when a nearby contaminant
requires a different fit. The `SeedAlignment` tree records the detector-wide
affine correction used to move the old seed windows onto the new raw spectra.
This accommodates a stale seed gain without using the old coefficients in the
final calibration.

`calibration_peak_adc_overrides.txt` provides detector- and peak-specific raw
ADC windows when an old high-order seed turns over before reaching a real
source line. The supplied entries recover detector 4's useful 2015.2-,
2034.8- and 2598.5-keV Co-56 lines. A `skip` entry suppresses a known
unreliable detector/peak combination. Window centres only locate peaks; they
are not themselves calibration measurements.

### Fit the main energy calibration

`fit_energy_calibration` converts the accepted raw centroids into the main
raw-ADC-to-keV calibration:

```bash
build/fit_energy_calibration \
    mid_source_peak_fits.root \
    calibrations/midexp_main.mcal \
    midexp_calibration_diagnostics.root \
    --config calibration_models.txt \
    --max-energy 2048
```

For detectors without an override, it robustly fits polynomial degrees 1--3
and selects the simplest statistically supported model using BIC. Isolated
outliers are excluded from the coefficient fit but retained in the diagnostic
`CalibrationResiduals` tree with `used=false`.

`calibration_models.txt` supplies explicit piecewise models for unusually
non-linear detectors. The distributed configuration uses two continuous
fourth-order pieces for detector 4 with a boundary near 1000 keV. The output
is a single `.mcal` stage whose input is raw ADC; use it as `--mcal` without a
preceding main `.cal` stage. Its final raw-ADC boundary corresponds to the
requested maximum energy, so hits above the calibrated range are rejected.

The diagnostic ROOT file contains the selected model and coefficients, every
point and residual, excluded-point flags, and per-detector calibration,
residual and fitted-curve graphs. The high-energy Co-56 points participate in
the fit even when the output is limited to 2048 keV. If a detector's final
accepted source line lies just below that limit, the program adds a short,
continuous linear extension based on the robust high-energy slope instead of
extrapolating a high-order polynomial beyond its measured range.

### Fit experiment-wide gain drift

After making the main source calibration, `build_drift_spectra` applies it to
the experimental trees and accumulates high-resolution calibrated spectra in
successive chronological groups. Input blocks are processed in command-line
order, not alphabetically across the whole experiment. Each block starts a new
group, and the default group size is ten files.

For the present experiment:

```bash
build/build_drift_spectra drift_spectra.root \
    --mcal calibrations/midexp_main.mcal \
    --group-size 10 \
    --threads 24 \
    --break-before 9um 213 \
    --block 9um 1 642 \
      '/data4/cj0041/IFIN-HH/102Mo/rootbuild/9um/Run_9um_*.root' \
    --block 12um 1 634 \
      '/data4/cj0041/IFIN-HH/102Mo/rootbuild/12um/Run_12um_*.root' \
    --block 15um 1 363 \
      '/data4/cj0041/IFIN-HH/102Mo/rootbuild/15um/Run_15um_*.root' \
    --block 15um 364 575 \
      '/data4/cj0041/IFIN-HH/102Mo/rootbuild/15um/Run_15um_*.root' \
    --block 25um 1 700 \
      '/data4/cj0041/IFIN-HH/102Mo/rootbuild/25um/Run_25um_*.root' \
    --block 20um 1 1275 \
      '/data4/cj0041/IFIN-HH/102Mo/rootbuild/20um/Run_20um_*.root'
```

The explicit boundary before `9um` run 213 prevents one ten-file spectrum
from averaging across the known detector-2 gain jump. If inspection shows
that the jump occurs immediately after run 213, change it to
`--break-before 9um 214` and rebuild the drift spectra.

Fit the unshifted reference peaks and create the fine corrections with:

```bash
build/fit_drift_calibration \
    drift_spectra.root \
    calibrations/drift \
    drift_calibration_diagnostics.root \
    drift_peaks.txt \
    --series-order 9um,12um,15um,25um,20um
```

`--series-order` controls interpolation and smoothing order independently of
the order in which the existing `drift_spectra.root` was built. Therefore, an
older drift file made with 20um before 25um can be refitted correctly without
rereading the raw trees. New drift-spectrum files should still be built in the
true acquisition order shown above.

The supplied `drift_peaks.txt` contains the 136.269-, 165.400- and
301.570-keV Ta-181 lines, the 535.561-keV Mo-100 line, the provisional
627.540-keV Sn-114 line, and the 1014.520-, 1460.820- and 1778.987-keV
activation/background references. Each line has:

```text
ENERGY_KEV LABEL SEARCH_HALF_WIDTH_KEV FIT_HALF_WIDTH_KEV [EXCLUDED_IDS]
```

Detector IDs in the final comma-separated field are not offered that peak.
The distributed file excludes detector 4 from the 136- and 165-keV anchors.

For every detector and ten-file group, the fitter measures all usable peaks
after first performing a broad detector/group-wide affine alignment of the
reference pattern. The alignment searches gain scales from approximately
0.94 to 1.06 and offsets of approximately -15 to +15 keV, so a real gain jump
cannot place all of the true peaks outside their individual local windows.
The locally fitted centroids then robustly determine the affine fine adjustment
`E_corrected = C0 + C1*E_main`. A fit requires at least two accepted peaks
spanning at least 400 keV by default. Isolated residuals above 1.5 keV are
removed while enough anchors remain. Missing group fits are interpolated
between the nearest successful groups, or use the nearest successful group at
an endpoint. These choices and all peak-fit results are recorded in
`drift_calibration_diagnostics.root`.

Detector 4 is treated as a controlled special case because its main
piecewise source calibration has a repeatable response feature in the gap
between source anchors. Its time-dependent offset and gain are fitted only
with the stable 301-, 1015-, 1461-, and 1779-keV in-beam anchors (at least
three spanning at least 700 keV). Residuals at the 301-, 536-, 628-, 1015-,
1461- and 1779-keV references then define a full-range piecewise-linear
correction for each group. To obtain enough statistics without broadening a
peak through unresolved gain drift, the 536- and 628-keV spectra from up to
five neighbouring groups in the same target-distance series are first aligned
with their individual affine corrections and then summed. The strong 536-keV
anchor retains a local time trajectory. Because the 628-keV line is much
weaker, its accepted measurements constrain only a broad, robust curvature
relative to the 536--1015-keV interpolation; an individual weak fit can no
longer move one ten-file group directly and broaden the accumulated spectrum.
Series without a measurable 628-keV line inherit the nearest supported
curvature. These shape anchors never enter the affine gain fit. Segment slopes
are checked for monotonicity and reasonable scale. This follows detector 4's
changing response over the full analysis range without allowing false
centroids such as the group-33 1015-keV fit to enter the calibration. The
individual window fits and their neighbour counts are stored in
`Detector4WindowedPeakFits`; the applied 628-keV curvature is stored in
`Detector4PiecewiseCorrections`.

Detector 16 is also checked group by group for nonlinear drift. Mild curvature
may use a guarded quadratic correction. For the sharp response change seen in
the affected late-9um groups, an affine residual RMS above 0.8 keV selects a
continuous piecewise-linear mapping through the accepted anchors, provided at
least four anchors span 700 keV and every segment remains monotonic with a
reasonable slope. This avoids forcing the localised change into one unstable
high-order polynomial. Ordinary groups retain their affine correction.

The diagnostic `DriftSeedAlignment` tree stores the broad alignment for each
detector and group. Re-running `fit_drift_calibration` never rereads the raw
experimental trees; it operates only on `drift_spectra.root`.
`DriftCorrectedPeakResiduals` stores both the affine-only and final residual
for every accepted reference peak. Its `detector4PiecewiseApplied`,
`quadraticApplied` and `piecewiseApplied` branches identify the nonlinear
corrections. `detector4TrajectoryDifferenceKeV` and
`detector4TrajectoryOutlier` distinguish false fitted centroids from points
that were allowed to follow the detector-4 time trajectory.
`Detector4PiecewiseCorrections` and
`Detector16PolynomialCorrections` record the group-by-group choices.

### Plot calibration diagnostics

After the main and drift calibrations have been fitted, make readable
per-detector calibration, resolution, and drift plots with:

```bash
build/plot_calibration_diagnostics \
    mid_source_peak_fits.root \
    midexp_calibration_diagnostics.root \
    drift_calibration_diagnostics.root \
    calibration_plots
```

This creates three multipage PDFs:

- `calibration_plots/energy_calibrations.pdf`: raw ADC versus reference energy,
  the selected calibration curve, used/rejected points, and residuals;
- `calibration_plots/energy_resolution.pdf`: FWHM in keV versus reference
  energy, calculated from each accepted source-peak width using the local
  derivative of that detector's calibration, with a smooth resolution trend;
- `calibration_plots/energy_drift.pdf`: applied offset and gain corrections
  versus the true chronological drift-group position, with directly fitted
  groups distinguished from interpolated groups.

The same canvases are saved in
`calibration_plots/calibration_diagnostic_plots.root` for interactive editing
in ROOT. This plotting step reads only the three diagnostic files and does not
reprocess any raw events.

Because run numbers repeat at different target distances, the command writes
one existing-format run-by-run file per series:

```text
calibrations/drift_9um.cal
calibrations/drift_12um.cal
calibrations/drift_15um.cal
calibrations/drift_20um.cal
calibrations/drift_25um.cal
```

For nonlinear corrections, the same command also writes one run-by-run
piecewise file per series (the historical `_det4` filename is retained for
command-line compatibility):

```text
calibrations/drift_9um_det4.mcal
calibrations/drift_12um_det4.mcal
calibrations/drift_15um_det4.mcal
calibrations/drift_20um_det4.mcal
calibrations/drift_25um_det4.mcal
```

Detector 4 is omitted from the corresponding `.cal` file because its affine
drift and time-dependent full-range adjustment are both contained in
`_det4.mcal`.
Detector 16 is also omitted for the individual runs that require its
piecewise-linear correction; those entries are stored in the same `.mcal`.
Other detector/run combinations pass through unchanged.

Apply the appropriate file after the main calibration when building each
series of analysis trees. For example:

```bash
build/build_analysis_tree analysed/9um \
    --mcal calibrations/midexp_main.mcal \
    --run-by-run-cal calibrations/drift_9um.cal \
    --run-by-run-mcal calibrations/drift_9um_det4.mcal \
    --threads 24 \
    Run_9um_*.root
```

The main `.mcal` converts raw ADC to energy. The generated `.cal` and detector
4 `.mcal` are second-stage corrections whose input is that already-calibrated
energy. Supply the `.cal` first and `_det4.mcal` second as shown. Runs or
detectors omitted from either fine-adjustment file retain the energy from the
preceding stage.

### Recommended two-stage analysis

Stage one reads each raw file, applies germanium calibration and the configured
timing, BGO-veto, silicon, exclusion, multiplicity, and FoldValid logic, then
writes one compact analysis tree per input run:

```bash
mkdir -p analysed
build/build_analysis_tree analysed \
    --run-cal 294 294 calibrations/152Eu.cal \
    --run-mcal 294 294 calibrations/pint_run_1-10.mcal \
    --threads 8 \
    Run_9um_000294.root
```

This creates:

```text
analysed/Run_9um_000294_analysis.root
```

`build_analysis_tree --threads N` processes up to `N` different input files at
once. A single input file therefore uses one file worker; the parallelism is
intended for a multi-run batch. Each output retains its event structure and
stores calibrated and raw Ge energies, separated detector-hit vectors,
relative times, per-file absolute-time metadata, multiplicities, BGO-veto
flags, the silicon condition, and FoldValid.

Stage two reads any number of analysis trees and constructs the usual
histograms:

```bash
build/analyse_tree spectra.root \
    --threads 8 \
    --no-diagnostics \
    analysed/Run_9um_*_analysis.root
```

Stage two uses within-tree multithreading. Calibration and Ge exclusions belong
to stage one and are therefore not accepted by `analyse_tree`. Analysis-tree
files made with different timing/gate/exclusion configurations, or a mixture of
calibrated and uncalibrated files, are rejected rather than silently combined.
Analysis trees made by an older schema must be rebuilt before use.

To process a larger run range:

```bash
build/build_analysis_tree analysed \
    --run-cal 294 350 calibrations/coarse.cal \
    --run-mcal 294 350 calibrations/fine.mcal \
    --threads 24 \
    Run_9um_*.root

build/analyse_tree spectra.root --threads 24 \
    analysed/Run_9um_*_analysis.root
```

For a single calibration file containing separate coefficients for every run,
stage one can instead use the run-by-run options described below:

```bash
build/build_analysis_tree analysed \
    --cal calibrations/main.cal \
    --mcal calibrations/main_fine.mcal \
    --run-by-run-cal calibrations/drift_by_run.cal \
    --threads 24 \
    Run_9um_*.root
```

### Direct raw analysis

The original single-stage analyser remains available for validation and for
small one-off jobs:

The general command is:

```bash
build/analyse_raw OUTPUT.root [OPTIONS] INPUT.root [INPUT2.root ...]
```

Shell wildcards can be used for multiple input files:

```bash
build/analyse_raw analysis.root Run_30um_*.root
```

To validate the two-stage result against a direct result made with the same
calibration and conditions:

```bash
root -l -b -q \
    'compare_histograms.C("direct.root","two_stage.root")'
```

The comparison walks all histogram directories and reports any bin-content or
bin-error differences. Presentation canvases are intentionally skipped.

### Multiple CPU cores

Use `--threads N` to process independent ROOT entry ranges in parallel:

```bash
build/analyse_raw analysis.root --threads 8 Run_30um_*.root
```

Each worker fills its own histograms and diagnostics; the program merges them
before writing the single output file. The default is 4 threads when
`--threads` is omitted. Any positive thread count can be selected explicitly;
in particular, `--threads 1` forces serial processing. Start with the number
of physical CPU cores available on the analysis machine. Storage throughput
and the memory needed for one histogram set per active worker can limit useful
scaling. The same default and override apply to `build_analysis_tree` and
`analyse_tree`.

Every successful run prints both the event-processing time and the total
analysis time, together with the corresponding event rates. This makes, for
example, `--threads 1` and `--threads 8` runs directly comparable.

During event processing, `analyse_raw` and `analyse_tree` display the number
of processed and total events, percentage complete, average event rate,
elapsed time, and estimated time remaining. Progress output is enabled by
default for both serial and threaded runs. Disable it when desired with:

```bash
build/analyse_raw analysis.root --no-progress Run_30um_*.root
```

Use `--progress` to explicitly re-enable it. If both switches are given, the
last one on the command line takes effect.

### Germanium detector and running-time plots

Individual germanium energy spectra are written under
`Germanium/Energy/Individual`. `Raw` contains 16,384-channel ADC spectra and
`Calibrated` contains the corresponding calibrated spectra from 0 to 2048 keV.
Both use the same accepted, in-time Ge hits.

Additional silicon-conditioned products are written during histogramming:

```text
Germanium/Energy/h1_Ge_Si_M1_E
Germanium/Energy/Rings/h1_Ge_Si_E_ID00-04
Germanium/Energy/Rings/h1_Ge_Si_E_ID05-09
Germanium/Energy/Rings/h1_Ge_Si_E_ID10-14
Germanium/Energy/Rings/h1_Ge_Si_E_ID15-19
Germanium/Energy/Rings/h1_Ge_Si_E_ID20-24
Coincidences/ParticleGamma/h2_Ge_Si_EvE
```

These products use germanium hits that pass the Ge timing gate and BGO veto.
The silicon hits pass both the configured silicon timing and energy gates. In
the two-dimensional particle-gamma matrix, every accepted silicon hit is
paired with every accepted germanium hit in the event: two silicon hits and
four germanium hits therefore contribute eight entries. The multiplicity-one
spectrum requires the pre-veto valid Ge multiplicity to be exactly one; its
single Ge hit must then survive the BGO veto. The five ring spectra correspond
to detector IDs 0--4, 5--9, 10--14, 15--19 and 20--24 respectively.

Before event processing, the program reads the first and last `absoluteTime`
from each input file. It concatenates the elapsed time within each file to form
a continuous total-running-time coordinate. The all-detector drift matrix is
stored in `Germanium/Time`, and individual-detector energy-versus-running-time
matrices are stored in `Germanium/Time/Individual`. The total running time is
always printed; detailed per-file absolute-time ranges are included when
diagnostics are enabled.

If an event has an `absoluteTime` outside its file's recorded first-to-last
range, analysis continues. The event still contributes to all energy,
multiplicity and coincidence products; it is omitted only from the total-
running-time matrices. The number of affected events is printed at the end.

### Calibration

One or more global calibration stages can be supplied. They are applied from
left to right in command-line order:

```bash
build/analyse_raw analysis.root \
    --cal calibrations/152Eu.cal \
    --mcal calibrations/pint_run_1-10.mcal \
    Run_30um_000001.root
```

Different calibration chains can be assigned to inclusive run ranges:

```bash
build/analyse_raw analysis.root \
    --run-cal 1 10 coarse_1-10.cal \
    --run-mcal 1 10 fine_1-10.mcal \
    --run-cal 11 20 coarse_11-20.cal \
    --run-mcal 11 20 fine_11-20.mcal \
    Run_30um_*.root
```

Stages with the same run range are applied in command-line order. Run ranges
must not overlap, and every input run must be covered when run-dependent
calibration is used.

A large run-by-run fine-adjustment file is selected with
`--run-by-run-cal FILE` or `--run-by-run-mcal FILE`. It is applied after the
main calibration, so its polynomial input is the already-calibrated energy.
Each non-comment line starts with its run number and then uses the existing
`.cal` or `.mcal` detector format. For example:

```text
# RUN  ID  NCOEFF  C0       C1
294    0   2       0.284    0.198186
294    1   2       0.507    0.179789
295    0   2       0.301    0.198021
295    1   2       0.492    0.179910
```

and a run-by-run piecewise file uses:

```text
# RUN ID NPIECES  NCOEFF COEFFICIENTS... END  [NCOEFF COEFFICIENTS... END ...]
294   0  2        2      0.1 0.2        1000  2      0.2 0.199      2048
```

The optional legacy detector-group value is also accepted between `RUN` and
`ID`. Multiple run-by-run files form successive fine-adjustment stages in
command-line order. A detector ID omitted from a fine-adjustment stage keeps
the energy produced by the preceding stage, and a run omitted from the file
uses only its main calibration. Missing detector IDs in the main calibration
still reject events containing those Ge hits.

The main calibration may use global `--cal`/`--mcal` options or the existing
`--run-cal`/`--run-mcal` ranges. Global and run-ranged *main* calibrations
cannot be mixed with each other, but either form can be followed by run-by-run
fine adjustments. With no calibration option, raw energy values are used.

The run number is read from the final padded numeric field in an input filename
ending in `_XXXXXX.root`.

### Coincidence gates and angular matrices

Both analysis commands accept an editable coincidence-gate file:

```bash
build/analyse_raw spectra.root \
    --gates coincidence_gates.txt \
    --threads 24 \
    Run_30um_*.root
```

or, for the two-stage path:

```bash
build/analyse_tree spectra.root \
    --gates coincidence_gates.txt \
    --threads 24 \
    analysed/Run_30um_*_analysis.root
```

The supplied `coincidence_gates.txt` has four sections:

```text
[Symmetric]
Gate296 294 299 288 293 306 311
Gate159 157 162 151 156 163 168

[AllvFW]
Gate296 294 299 288 293 306 311
Gate159 157 162 151 156 163 168

[AllvBW]
Gate296 294 299 288 293 306 311
Gate159 157 162 151 156 163 168

[Double]
# Add triple-coincidence/double-gate definitions here
```

Each gate line is:

```text
NAME PROMPT_MIN PROMPT_MAX LOWER_MIN LOWER_MAX UPPER_MIN UPPER_MAX
```

All windows are half-open: `MIN <= energy < MAX`. Thus `294 299`
contains the 1-keV channels 294 through 298. The supplied starter gates are
296 keV (prompt 294--298, sidebands 288--292 and 306--310) and 159 keV
(prompt 157--161, sidebands 151--155 and 163--167). The program subtracts both
sidebands and automatically scales them by
`prompt width / combined sideband width`. Any number of uniquely named gates
can be included in each section. Lines beginning with `#` and blank lines are
ignored.

The supplied 296- and 159-keV gates are present in `[Symmetric]`, `[AllvFW]`
and `[AllvBW]`. The symmetric entries project the other gamma from the regular
Ge-Ge matrix. The angular entries gate the ALL axis and project the forward or
backward detector axis respectively.

All configurable gated spectra are calibrated when calibration is supplied,
BGO-vetoed, and silicon-conditioned. `[Symmetric]` gates the regular
symmetrised gamma-gamma matrix. For `[AllvFW]` and `[AllvBW]`, the gate is
placed on the `ALL` x-axis and the resulting spectrum shows the forward or
backward y-axis respectively.

The conditioned angular matrices and their ungated `ALL`-axis projections are
written directly in `Coincidences`:

```text
h2_Ge_Si_gg_AllvFW
h1_Ge_Si_gg_AllvFW_proj
h2_Ge_Si_gg_AllvBW
h1_Ge_Si_gg_AllvBW_proj
```

Configured gated spectra are written beneath `Coincidences/Gated`, with
`AllvFW` and `AllvBW` subdirectories for the angular spectra.

The optional `[Double]` section creates triple-coincidence products without a
full gamma cube. Its lines contain:

```text
NAME REQUIRED_MIN REQUIRED_MAX PROMPT_MIN PROMPT_MAX LOWER_MIN LOWER_MAX UPPER_MIN UPPER_MAX
```

An event must contain at least three BGO-vetoed Ge hits and satisfy the silicon
condition. One gamma in the half-open `REQUIRED_MIN--REQUIRED_MAX` window is
removed. The remaining gammas fill a regular symmetrised gamma-gamma matrix;
the prompt and scaled-sideband windows are then applied to that matrix to make
the double-gated spectrum. If several gammas satisfy the required window, the
event is entered once and the first matching hit is used. Matrices and spectra
are stored together in `Coincidences/DoubleGated` as:

```text
h2_Ge_Si_gg_DG_NAME
h1_Ge_Si_gg_DG_NAME_proj
```

For example:

```text
[Double]
TripleGate 390 410 540 560 520 535 565 580
```

requires a 390--410 keV gamma, forms the matrix from the other gamma rays, and
then gates 540--560 keV with the stated lower and upper sidebands.

### Excluding germanium detectors

Use `--exclude-ge` with a comma-separated list of zero-based Ge ID numbers:

```bash
build/analyse_raw analysis.root \
    --exclude-ge 0,7,12 \
    Run_30um_*.root
```

The option can be repeated:

```bash
build/analyse_raw analysis.root \
    --exclude-ge 0,7 \
    --exclude-ge 12 \
    Run_30um_*.root
```

Valid Ge ID numbers are 0--24.

### Diagnostic output

Detailed terminal diagnostics are enabled by default. Disable them with:

```bash
build/analyse_raw analysis.root --no-diagnostics Run_30um_*.root
```

They can be explicitly enabled with `--diagnostics`. If both switches are
given, the last one on the command line takes effect.

Options can be combined in one command:

```bash
build/analyse_raw analysis.root \
    --no-diagnostics \
    --exclude-ge 7 \
    --cal calibrations/152Eu.cal \
    --mcal calibrations/pint_run_1-10.mcal \
    Run_30um_000001.root
```
