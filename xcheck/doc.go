// Package xcheck cross-validates the layers of Tetsuhiro-OZ against each
// other and against independent computations:
//
//   - E1–E5 (plan §6.4): the geometric layer (geo, lens) against the
//     wave-optics kernel (optics, via the wave bridge);
//   - E6 (plan §6.4): the wave kernel's diffractive elements (zone plate,
//     blazed grating) against their design formulas;
//   - G7 joint check (plan §4.2): an AR-coated achromat doublet evaluated by
//     the lens joint path against a separate computation built from the geo
//     and coating packages alone.
//
// The package consists of test files only; every check is an automatic
// assertion with numerical thresholds.
package xcheck
