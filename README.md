# Raw ROOT-tree analysis

This program analyses the `events` tree. Each tree entry is one built event,
and equal vector indices describe the fields of one detector hit.

## Input

| Branch | ROOT/C++ type |
| --- | --- |
| `detectorLUT` | `std::vector<UShort_t>` |
| `detectorType` | `std::vector<UShort_t>` |
| `energy` | `std::vector<UShort_t>` |
| `psd` | `std::vector<UShort_t>` |
| `relativeNsTime` | `std::vector<UShort_t>` |
| `relativePsTime` | `std::vector<UShort_t>` |
| `absoluteTime` | scalar `ULong_t` (`/g`) |

Relative time is reconstructed in ns as

```text
(1000 * relativeNsTime + relativePsTime) / 1000
```

The detector-type mapping in `AnalysisConfig.h` is:

- `0`: germanium
- `1`: silicon
- `2`: BGO
- `3`: LaBr3

## Definitions and selections

- `EventMultiplicity` is the number of hits stored in the event vectors.
- `GermaniumMultiplicity` is the number of Ge hits passing the inclusive
  985--1100 ns Ge timing gate and, when calibration is enabled, having a valid
  calibration result. It is evaluated before BGO vetoing.
- `GermaniumMultiplicityAfterBGOVeto` counts those valid Ge hits that also
  survive the per-LUT BGO veto.
- BGO, silicon, and LaBr3 multiplicities count all recognised hits of their
  respective type. Their timing and energy conditions do not define hit
  validity.
- `FoldValid` is true exactly when
  `BGOMultiplicity <= GermaniumMultiplicity`.

An in-time BGO hit (850--1150 ns) vetoes a Ge hit with the same LUT. The
silicon condition is event-wide and requires at least one silicon hit with
relative time 940--1050 ns and raw energy 1000--4000. The non-Ge hits remain in
their ordinary spectra and multiplicities whether or not they pass these
condition windows.

For diagnostics, BGO/Ge LUT matching is also calculated one-to-one. The final
report states how many `FoldValid` events have every BGO hit uniquely matched
to a same-LUT Ge hit. This matching result does not alter `FoldValid`. A second
diagnostic reports events containing multiple valid Ge hits in one LUT.

All boundaries are inclusive unless a configured interval is explicitly
described as half-open.

## Build and run

In a ROOT-enabled shell:

```bash
mkdir build
cd build
cmake ..
cmake --build . -j
```

`TreePlayer` is explicitly linked because it supplies `TTreeReader`.

Examples:

```bash
./analyse_raw output.root Run_30cm_000123.root
./analyse_raw output.root Run_30cm_*.root
```

Detailed terminal diagnostics are enabled by default. To retain the normal
progress and output-file messages but suppress the diagnostic report, use:

```bash
./analyse_raw output.root --no-diagnostics Run_30cm_*.root
```

`--diagnostics` explicitly enables them again; if both switches are supplied,
the last one on the command line takes effect.

One or more zero-based germanium detector LUTs can be excluded with a
comma-separated list:

```bash
./analyse_raw output.root --exclude-ge 0,7,12 Run_30cm_*.root
```

The option may be repeated. Excluded LUTs (valid range 0--24) are removed from
all Ge calibration checks, timing statistics, multiplicities, spectra,
matrices, projections, BGO-veto decisions, and `FoldValid`. Other detector
types are unchanged, and `EventMultiplicity` remains the number of hits stored
in the original event vectors.

## Germanium calibration

Calibration stages are applied from left to right:

```bash
./analyse_raw output.root \
    --cal ../calibrations/152Eu.cal \
    --mcal ../calibrations/pint_run_1-10.mcal \
    Run_30cm_000123.root
```

Run-dependent inclusive ranges are supported:

```bash
./analyse_raw output.root \
    --run-cal 1 10 coarse_1-10.cal \
    --run-mcal 1 10 fine_1-10.mcal \
    --run-cal 11 20 coarse_11-20.cal \
    --run-mcal 11 20 fine_11-20.mcal \
    Run_30cm_*.root
```

Distinct ranges must not overlap and every input run must be covered. Global
and run-ranged options cannot be mixed. With no calibration arguments, the
analysis intentionally uses raw energy.

A `.cal` line contains

```text
[optional detectorGroup] detectorLUT coefficientCount a0 a1 a2 ...
```

and evaluates `a0 + a1*x + a2*x^2 + ...`. An `.mcal` line contains

```text
[optional detectorGroup] detectorLUT pieceCount \
    coefficientCount coefficients... pieceEnd [next piece ...]
```

The first piece begins at zero. Piece boundaries are evaluated using the value
entering that calibration stage. The optional leading detector-group value is
accepted and ignored.

If any timing-valid Ge hit has no LUT calibration, lies outside an `.mcal`
range, or produces a non-finite result, the entire event is ignored. Rejection
counts are printed.

## ROOT output structure

Only the following analysis objects are written:

```text
/
  h1_EventMultiplicity
  h1_UnknownDetectorType
  Germanium/
    h1_Ge_noVeto_Multiplicity
    h1_Ge_Multiplicity
    h1_Ge_noVeto_LUT
    Energy/
      h1_Ge_noVeto_E
      h1_Ge_noVeto_gg_proj
      h1_Ge_Si_gg_proj
      h1_Ge_E
      h1_Ge_Si_E
      c_Ge_Si_E_Comparison
      FoldValid/
        h1_Ge_Si_FV_Multiplicity
        h1_Ge_Si_FV_E
        h1_Ge_Si_E
      Individual/
        h1_Ge_noVeto_E_LUTxx ...
    Time/
      h1_Ge_noVeto_Time
      h2_Ge_noVeto_EvTime
  Silicon/
    h1_Si_Multiplicity
    h1_Si_LUT
    Energy/
      h1_Si_E
      Individual/
        h1_Si_E_LUTxx ...
    Time/
      h1_Si_Time
  BGO/
    h1_BGO_Multiplicity
    h1_BGO_LUT
    Energy/
      h1_BGO_E
      Individual/
        h1_BGO_E_LUTxx ...
    Time/
      h1_BGO_Time
  LaBr/
    h1_LaBr_Multiplicity
    h1_LaBr_LUT
    Energy/
      h1_LaBr_E
      Individual/
        h1_LaBr_E_LUTxx ...
    Time/
      h1_LaBr_Time
  Coincidences/
    h2_Ge_gg
    h2_Ge_Si_gg
    Gated/
      h1_Ge_Si_gg_Gate292to299_proj
      AllvFW/
      AllvBW/
```

The ordinary Ge spectrum and individual Ge spectra apply only timing and any
configured calibration. `h1_Ge_noVeto_gg_proj` is the projection that would be
obtained from a symmetrised matrix of all distinct valid Ge-hit pairs, with no
BGO-veto, silicon, or `FoldValid` condition. Consequently, in a
Ge-multiplicity-M event each hit contributes M-1 times; events with fewer than
two valid Ge hits make no contribution. `h1_Ge_E` adds the BGO veto, and
`h1_Ge_Si_E` also adds the silicon condition.
The comparison canvas area-normalises the former spectrum to the latter and
draws both with `HIST`; the stored analysis spectra are not rescaled.

`h1_Ge_Si_gg_proj` repeats the symmetrised-pair projection using only Ge hits
surviving the per-LUT BGO veto and only events passing the silicon
timing-and-energy condition. It does not require `FoldValid`. The two temporary
2D matrices formerly stored beside these projections in `Germanium/Energy`
are no longer created or written.

Everything under `Germanium/Energy/FoldValid` has the BGO veto and silicon
condition. The multiplicity there is the post-veto Ge multiplicity.
`h1_Ge_Si_FV_E` requires `FoldValid`; the `h1_Ge_Si_E` stored in the same
subdirectory accepts either boolean value.

The energy-versus-time matrix uses calibrated Ge energy when calibration is
enabled and interprets `absoluteTime` as seconds. Its current axes are editable
in `AnalysisConfig.h`.

The gamma-gamma matrices use distinct pairs of BGO-veto-surviving Ge hits and
are symmetrised by filling both axis orders. The second matrix additionally
requires the event-wide silicon condition. No `FoldValid` requirement is
applied to either matrix.

The one current gated spectrum requires BGO-vetoed Ge hits and the silicon
condition. It gates on 292--299 keV and subtracts the 283--289 and 310--317 keV
sidebands with width normalisation:

```text
prompt - (8/15) * (lower + upper)
```

It has `Sumw2` enabled, so negative bins are legitimate. Add or edit gates in
`AnalysisConfig.h`.

The RDDS `AllvFW` and `AllvBW` directories are deliberately empty while
`kRddsGates` is empty. Adding gates there creates only the corresponding
forward/backward 1D gated spectra; the all-vs-angle 2D matrices are not stored.
These future spectra use BGO-veto-surviving Ge hits from silicon-conditioned
events. Forward is LUT 0--4 (37 degrees) and backward is LUT 20--24
(143 degrees).

All 1D histograms store ROOT's `HIST` draw option. All retained 2D histograms
store `COLZ`. Energy histograms use one stored-energy unit per bin; this is one
keV per bin for calibrated germanium. LUT occupancy axes cover IDs 0--24.

## Terminal report

The program prints one combined report for all input files:

- detector multiplicity distributions and event-multiplicity composition;
- BGO multiplicity distributions for Ge multiplicities 1--4;
- the `FoldValid` and `FoldValid`-plus-well-matched fractions;
- repeated-Ge-LUT events;
- raw Ge/BGO/silicon gate-pass statistics and BGO-veto survival;
- malformed events, unknown detector types, non-zero PSD values, calibration
  failures, and the observed `absoluteTime` range.

The Ge gate-pass diagnostic is evaluated before whole-event calibration
rejection. Ordinary output histograms contain accepted events only.
