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

The executable is then:

```text
build/analyse_raw
```

To perform a completely clean rebuild:

```bash
rm -rf build
cmake -S . -B build
cmake --build build -j
```

## Usage

The general command is:

```bash
build/analyse_raw OUTPUT.root [OPTIONS] INPUT.root [INPUT2.root ...]
```

Shell wildcards can be used for multiple input files:

```bash
build/analyse_raw analysis.root Run_30um_*.root
```

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
