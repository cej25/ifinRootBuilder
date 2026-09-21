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
```

To perform a completely clean rebuild:

```bash
rm -rf build
cmake -S . -B build
cmake --build build -j
```

## Usage

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
before writing the single output file. Omitting the option, or using
`--threads 1`, retains serial processing. Start with the number of physical CPU
cores available on the analysis machine. Storage throughput and the memory
needed for one histogram set per active worker can limit useful scaling.

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

Before event processing, the program reads the first and last `absoluteTime`
from each input file. It concatenates the elapsed time within each file to form
a continuous total-running-time coordinate. The all-detector drift matrix is
stored in `Germanium/Time`, and individual-detector energy-versus-running-time
matrices are stored in `Germanium/Time/Individual`. The total running time is
always printed; detailed per-file absolute-time ranges are included when
diagnostics are enabled.

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
calibration is used. Global and run-dependent calibration options cannot be
mixed in one invocation. With no calibration option, raw energy values are
used.

The run number is read from the final padded numeric field in an input filename
ending in `_XXXXXX.root`.

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
