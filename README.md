# EPANET-CGGA

**An adaptive hybrid hydraulic solver for unsteady pipe network analysis, implemented as an object-oriented extension of EPANET 3.**

[![build-and-test](https://github.com/MehmetMelihKosucu/EPANET-CGGA/actions/workflows/ci.yml/badge.svg)](https://github.com/MehmetMelihKosucu/EPANET-CGGA/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Language: C++](https://img.shields.io/badge/Language-C%2B%2B-blue.svg)](https://isocpp.org/)
[![Build: CMake](https://img.shields.io/badge/Build-CMake-064F8C.svg)](https://cmake.org/)
[![Platform: Windows](https://img.shields.io/badge/Platform-Windows-0078D6.svg)](#building)
<!-- [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20199101.svg)](https://doi.org/10.5281/zenodo.20199101) -->

---

## Overview

EPANET-CGGA is a C++ extension of the EPANET 3 hydraulic engine that implements the **Comprehensive Global Gradient Algorithm (CGGA)** of Nault & Karney (2020) within a unified, adaptive hybrid framework. The solver automatically switches between three computational modes on a per-timestep basis:

- **Quasi-Steady (QS)** — steady-state hydraulics for slowly varying demand patterns.
- **Rigid Water Column (RWC)** — inertial unsteady flow neglecting fluid compressibility.
- **Water Hammer (WH)** — fully compressible unsteady flow via the Algebraic Water Hammer (AWH) method.

Mode selection is driven by physical unsteadiness indicators (φA, φR) evaluated against dynamic and inertial thresholds, with a **redo-based escalation pattern** that probes higher accuracy modes when a hydraulic event is detected and a **deferred de-escalation with event-driven lockout** to prevent oscillatory mode switching.

EPANET-CGGA is intended for researchers and practitioners who need transient analysis capability inside an extended-period simulation (EPS) workflow without maintaining separate tools for steady-state and water-hammer analyses.

## Key Features

- **Adaptive hybrid solver** with three coexisting flow modes in a single simulation run.
- **Algebraic Water Hammer (AWH)** method treats each pipe as an undiscretized unit, avoiding the segmentation overhead of traditional Method of Characteristics (MOC) implementations.
- **Object-oriented integration** in the EPANET 3 (`epanet-dev`) codebase, preserving full compatibility with the standard EPANET 3 `.inp` input format.
- **History-management subsystem** with binary-search interpolation across mode transitions and resolution-aware gap handling.
- **MIT-licensed** open-source release with C++14 compatibility.

## Repository Structure

```
EPANET-CGGA/
├── .github/
│   └── workflows/
│       └── ci.yml          # Continuous integration: build and test
├── Networks/               # EPANET .inp benchmark networks
├── src/                    # C++14 source (EPANET 3 base + CGGA extension)
│   ├── Core/               # Engine, hydraulic balance, flow history
│   ├── Elements/           # Network elements
│   ├── Input/              # .inp parsers
│   ├── Models/             # Demand, headloss, leakage models
│   ├── Output/             # Report and binary output writers
│   ├── Solvers/            # GGA, RWC-GGA and CGGA solvers
│   └── Utilities/          # Support utilities
├── tests/                  # Unit tests for the flow history subsystem
│   ├── CMakeLists.txt
│   └── test_flowhistory.cpp
├── CITATION.cff
├── CMakeLists.txt
├── LICENSE
└── README.md
```

## Building

EPANET-CGGA targets **C++14** and is currently developed on Windows with Visual Studio. Build with CMake:

```bash
git clone https://github.com/MehmetMelihKosucu/EPANET-CGGA.git
cd EPANET-CGGA
mkdir build && cd build
cmake -G "Visual Studio 17 2022" ..
cmake --build . --config Release
```

Replace `"Visual Studio 17 2022"` with the generator string matching your installed Visual Studio version (run `cmake --help` to list available generators).

The build produces the executable in `build/bin/Release/`.

## Tests

The flow history subsystem is verified in isolation by a unit test suite
under `tests/`. Build with testing enabled and run:

```bash
cd build
ctest -C Release
```

The suite covers circular-buffer construction and capacity, the linear
interpolation used on the reach-back path, per-pipe history registration,
and the requirement that a reach-back query which cannot be satisfied from
stored data reports failure rather than returning extrapolated values.
Tests run automatically on every commit via GitHub Actions.

## Quick Start

EPANET-CGGA accepts standard EPANET `.inp` files. To enable adaptive hydraulic analysis, add the following option to the `[OPTIONS]` section of the input file:

```
[OPTIONS]
Hyd_Solver        CGGA
```

Run the solver from the command line:

```bash
run-epanet3.exe Networks/Onizuka1986-EPA3.inp output.rpt output.out
```

If `Hyd_Solver` is omitted or set to `GGA`, the solver falls back to the standard EPANET 3 steady-state global gradient algorithm.

### Configuration

The four unsteadiness thresholds governing mode selection are public members
of `CGGASolver` in `src/Solvers/cggasolver.h`:

| Parameter | Symbol | Default | Units |
|---|---|---|---|
| `phiA_D` | φ_A,D | 0.1 | m |
| `phiR_D` | φ_R,D | 0.08 | – |
| `phiA_I` | φ_A,I | 0.05 | m |
| `phiR_I` | φ_R,I | 0.04 | – |

The dynamic thresholds must exceed the inertial ones. Per-mode time steps are
the compile-time constants `Delta_WH` (0.01 s), `Delta_RWC` (1.0 s) and
`Delta_QS` (10.0 s) in the same header.

### Output

Alongside the standard EPANET report, the solver writes an unsteadiness
indicator log recording, at each integer second, the active mode and the
maximum absolute and relative indicators across the network.

### Wave Speed Specification

Pipe wave celerities (used in the Water Hammer mode) are computed by the `getWaveSpeed() const` member function defined in `src/Elements/pipe.cpp`. The function may be modified to adopt user-defined wave-speed models or to read wave celerities directly from the `.inp` file as a custom property.

## Benchmark Networks

The `Networks/` directory contains the input files for the case studies
reported in the accompanying paper, together with a classical reference case:

| Network | Size | Transient event | Source |
|---|---|---|---|
| `Tnet1-EPA3.inp` | 7 junctions, 9 pipes, 1 reservoir, 1 valve | Valve closure at a network boundary (CCV on a 400 mm line) | Streeter & Wylie (1967) |
| `Onizuka1986-EPA3.inp` | 6 junctions, 9 pipes, 3 reservoirs, 1 tank, 2 valves | Linear valve closure over 100 s, analysed at two tank diameters | Onizuka (1986) |
| `EPA3-Nault2016.inp` | 34 junctions, 44 pipes, 3 reservoirs, 8 pumps, 1 valve | Pump shutoff and restart in a looped network | Nault & Karney (2016) |
| `EPA3-hk-burst.inp` | 212 junctions, 231 pipes, 2 reservoirs, 2 valves | Pipe burst at 03:30 within a 24 h extended-period simulation | Hadımköy WDS, Istanbul |

The Hadımköy model includes pressure-dependent leakage. Pipes shorter than the
Courant-limited minimum length were extended so that a common water hammer time
step could be applied across the network; this modification is disclosed in the
accompanying paper.

## Comparison with Existing Tools

| Feature | EPANET-CGGA | TSNet | PTSNet | ALLIEVI | Bentley HAMMER |
|---|---|---|---|---|---|
| License | MIT | MIT | MIT | Proprietary | Commercial |
| Language | C++14 | Python | Python | Windows GUI | Windows GUI |
| Algorithm family | Adaptive CGGA | MOC | MOC | MOC | MOC |
| Adaptive QS/RWC/WH switching | ✅ | ❌ | ❌ | ❌ | ❌ |
| EPANET 3 native integration | ✅ | wraps EPANET externally | wraps EPANET externally | manual | import path |
| Cost to user | Free | Free | Free | License fee | License fee |

## Citation

**Software**:

```bibtex
@software{Kosucu2026EPANETCGGAcode,
  author    = {Kosucu, Mehmet Melih},
  title     = {EPANET-CGGA: Adaptive Hybrid Hydraulic Solver},
  year      = {2026},
  publisher = {Zenodo},
  version   = {v1.0.0},
  doi       = {10.5281/zenodo.20199101},
  url       = {https://github.com/MehmetMelihKosucu/EPANET-CGGA}
}
```

## References

The CGGA framework and its three-mode reduction are described in:

- Nault, J.D., Karney, B.W., 2020. Comprehensive adaptive modelling of 1-D unsteady pipe network hydraulics. Journal of Hydraulic Research 59, 263–279. https://doi.org/10.1080/00221686.2020.1770878.
- Nault, J.D., Karney, B.W., Jung, B., 2016. Algebraic Water Hammer: Global Formulation for Simulating Transient Pipe Network Hydraulics, in: World Environmental and Water Resources Congress 2016. EWRI, West Palm Beach FL, pp. 191–201.
- Streeter, V. L., & Wylie, E. B. (1967). *Hydraulic Transients*. McGraw-Hill.
- Wylie, E. B., & Streeter, V. L. (1993). *Fluid Transients in Systems*. Prentice Hall.
- Onizuka, K., 1986. System Dynamics Approach to Pipe Network Analysis. Journal of Hydraulic Engineering 112, 728–749. https://doi.org/10.1061/(asce)0733-9429(1986)112:8(728).

## Acknowledgments

EPANET-CGGA builds on the object-oriented EPANET 3 codebase developed by the [Open Water Analytics](https://github.com/OpenWaterAnalytics/epanet-dev) community and the original [EPANET](https://github.com/USEPA/EPANET2.2) hydraulic engine from the U.S. Environmental Protection Agency. The CGGA algorithm follows the formulation of Nault & Karney (2020). The author gratefully acknowledges these contributions.

## Contributors

- **Mehmet Melih Koşucu** — Istanbul Technical University, Department of Civil Engineering

Contributions, bug reports, and feature requests are welcome via [GitHub Issues](https://github.com/MehmetMelihKosucu/EPANET-CGGA/issues) and pull requests.

## License

EPANET-CGGA is released under the **MIT License**. See [LICENSE](LICENSE) for the full text. The choice of MIT preserves compatibility with the upstream EPANET (US EPA) and `epanet-dev` (OpenWaterAnalytics) license terms and permits both academic and commercial use.

## Contact

For questions or collaboration enquiries, please open a GitHub Issue or contact the author through Istanbul Technical University.

