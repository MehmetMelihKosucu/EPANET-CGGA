# Changelog

## v1.1.0 (2026-10-23)

This version produced the results of the revised manuscript.

### Changed
- The unsteady friction term is switched off by default. Friction is represented
  by the quasi-steady Darcy-Weisbach term, as in the quasi-steady friction method
  of Bentley HAMMER used for the reference solutions. The term can be switched on
  with the environment variable `CGGA_UNSTEADY_FRICTION=1`.
- When the unsteady friction term is on, the decay coefficient is the Vardy-Brown
  value k = √C*/2 (previously 3√C*).
- The solver source was updated to the version used for the paper
  (`cggasolver`, `hydbalance`, `dualflowhistory`, `epanet3`, `pipe`, `junction`).

### Added
- Regression tests that run complete simulations of the Onizuka (tank diameters
  of 1 m and 10 m), Nault and Karney, and Hadımköy networks and compare the head
  series, the peak and minimum heads, and the mode schedule with stored reference
  files (`tests/regression/`).
- The unit and regression tests run in the continuous integration workflow on
  every push to `main` and on every pull request.
- A manually triggered workflow that regenerates the regression reference files.
- The environment variable `CGGA_THRESHOLDS` sets the four unsteadiness
  thresholds at run time without recompiling.
- Input file `Networks/Onizuka1986-EPA3-D10.inp` (tank diameter of 10 m).

### Fixed
- The forced-step counters of `CGGASolver` are initialized in the constructor.

## v1.0.0 (2026-05-15)

Initial release.
