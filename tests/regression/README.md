# Regression tests

The regression tests run complete simulations of the benchmark networks in
`Networks/` with the adaptive CGGA solver and compare the results with
reference files stored in `tests/regression/reference/`. Their purpose is to
detect any change in the computed results caused by a change in the code.

## What is checked

For each network, `regression_test` runs the simulation from start to end and
checks three things against the reference file:

1. **Head series.** The head at the monitored nodes is compared at every
   reference time (every 0.1 s, or every step when the step is longer). The
   largest absolute difference must not exceed the tolerance (default
   0.01 m).
2. **Extreme heads.** The peak and minimum head at each monitored node must
   agree with the reference within the same tolerance.
3. **Mode schedule.** The sequence of quasi-steady (QS), rigid water column
   (RWC) and water hammer (WH) intervals must be the same, and each mode
   transition must occur within 1 s of the reference time. The mode of each
   interval is identified from the time step size (WH: dt < 0.5 s,
   RWC: 0.5 s ≤ dt < 5 s, QS: dt ≥ 5 s).

| Test | Input file | Monitored nodes |
|---|---|---|
| `regression_onizuka_d10` | `Onizuka1986-EPA3-D10.inp` (tank diameter 10 m) | J1, J2, T1 |
| `regression_onizuka_d1` | `Onizuka1986-EPA3.inp` (tank diameter 1 m) | J1, J2, T1 |
| `regression_nault2016` | `EPA3-Nault2016.inp` | J16, J24 |
| `regression_hadimkoy` | `EPA3-hk-burst.inp` (24 h) | 612, 3234 |

For the 24 h Hadımköy case the comparison ends at the end of the water hammer
interval that follows the burst (`--window first-wh`). After this interval the
solver passes through the RWC mode and returns to the QS mode when the
indicators fall below the inertial thresholds. In this network the indicators
stay close to the inertial thresholds for some time, and the instant of the
return to the QS mode differed between otherwise identical runs on different
`windows-latest` machines (13 880 s in one run and 14 052 s in another, with
identical results up to 13 878 s). Such a difference is caused by round-off,
not by a change in the code, so it is excluded from the comparison.

All runs use the threshold values compiled into `CGGASolver`
(`src/Solvers/cggasolver.h`). The tests therefore check that the code gives
the same results for fixed inputs and settings.

## Running the tests

```bash
cmake -B build -A x64
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure                 # all tests
ctest --test-dir build -C Release -L regression --output-on-failure  # regression tests only
```

All four regression tests run on every commit through the continuous
integration workflow. On the `windows-latest` runner, the four simulations
together take about three minutes.

The result of each run is written to `build/regression_results/<test>.csv`
in the same format as the reference file, so a failed comparison can be
inspected directly.

## Reference files

Each reference file is a CSV file with a short header:

```
# EPANET-CGGA regression reference
# input: EPA3-Nault2016.inp
# steps: <number of solver steps>
# modes: QS@0.000 WH@30.000 RWC@374.000 QS@380.000
# peak: J16=...@... J24=...@...
# min:  J16=...@... J24=...@...
time,J16,J24
0.0000,195.871000,196.750000
...
```

The reference files were generated with the Windows (MSVC, x64, Release)
build used by the continuous integration workflow. Results obtained with a
different compiler or platform may differ by more than the tolerance, so the
reference files should be compared only against the same type of build.

## Updating the reference files

When a change in the code is intended to change the results, the reference
files must be regenerated and committed together with the change:

```bash
cmake --build build --config Release --target regression_update
```

The same can be done on GitHub by running the **regression-references**
workflow manually (Actions tab). It builds the code on `windows-latest`,
regenerates all four reference files, and uploads them as an artifact.
