// Package coating implements the thin-film design core of Tetsuhiro-OZ
// (Optics-Zero): the characteristic-matrix (transfer-matrix) optics of
// isotropic multilayer coatings, a dispersion library of common coating and
// substrate materials, and the autonomous design machinery (damped least
// squares with thickness bounds, random restarts and needle-like layer
// growth) that turns a set of spectral targets into a layer stack.
//
// # Conventions
//
//   - Wavelengths and layer thicknesses are in nanometres (nm), incidence
//     angles in degrees, measured from the surface normal.
//
//   - The complex refractive index follows the physics convention
//
//     n = n' + i·k,   k >= 0 absorbing,
//
//     with the time dependence exp(-i·omega·t). A plane wave therefore carries
//     the factor exp(i·k0·n·z) with k0 = 2*pi/lambda, so an absorbing medium
//     attenuates it.
//
//   - s-polarization means E perpendicular to the plane of incidence (the
//     tangential field is Ey), p-polarization means E inside that plane (the
//     tangential field is Ex). "Unpolarized" results are the arithmetic mean of
//     the two.
//
//   - Reflectance, transmittance and absorptance are power (energy) ratios of
//     the semi-infinite substrate system: R + T + A = 1, with
//     A = 1 - R - T (the residue is clamped when it is pure rounding noise).
//
// # Physics
//
// Each layer j with complex index n_j, thickness d_j and free-space wavenumber
// k0 contributes the characteristic matrix (Macleod, "Thin-Film Optical
// Filters", 3rd ed., ch. 2)
//
//	M_j = [ cos(delta_j)        i·sin(delta_j)/eta_j ]
//	      [ i·eta_j·sin(delta_j)  cos(delta_j)       ]
//
// with the complex phase thickness delta_j = k0·d_j·n_j·cos(theta_j) and the
// tilted optical admittance
//
//	eta_j = n_j·cos(theta_j)      (s-polarization)
//	eta_j = n_j / cos(theta_j)    (p-polarization)
//
// where cos(theta_j) = sqrt(n_j^2 - alpha^2)/n_j and alpha = n_0·sin(theta_0)
// is the (complex) Snell invariant of the incidence medium. The stack matrix is
// the ordered product M = M_1·M_2···M_N, the input admittance is
//
//	Y = (M21 + M22·eta_s) / (M11 + M12·eta_s),
//
// and the interface coefficients are
//
//	r = (eta_0 - Y)/(eta_0 + Y),     t = 2·eta_0 / ((eta_0 + Y)·(M11 + M12·eta_s))
//
// with R = |r|^2 and T = Re(eta_s)/Re(eta_0)·|t|^2. For a single layer this
// reduces to the classic Airy summation r = (r1 + r2·exp(-2i·delta))/(1 +
// r1·r2·exp(-2i·delta)); the test suite checks that identity for random
// indices, thicknesses and angles.
//
// # Scope and limitations
//
//   - The substrate is semi-infinite and the whole stack is treated coherently
//     (one coherent pile). TODO: incoherent (intensity-summed) back-surface
//     reflection of a finite substrate, and incoherent thick-layer splitting,
//     are not implemented; a real plane-parallel plate shows an additional
//     reflection halo that this package ignores.
//   - Layers are isotropic. Anisotropic (birefringent) media are out of scope
//     here; Tetsuhiro-WOS implements the 4x4 Berreman method in
//     optics/berreman.go for that case.
//   - Only specular, plane-wave illumination of laterally uniform coatings is
//     modelled (no scattering, roughness, grading, or finite beam effects).
//   - The dispersive materials are literature values for bulk or typical thin
//     films; real evaporated layers shift with deposition conditions (see the
//     doc comment of every library entry for the source and validity range).
package coating
