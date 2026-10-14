package geo

import (
	"fmt"
	"math"
)

// ---------------------------------------------------------------------------
// Chromatic aberrations (paraxial)
// ---------------------------------------------------------------------------
//
// Colour aberrations are the change of the first-order image with wavelength.
// They are computed from the paraxial engine (System.Paraxial), which evaluates
// every medium at the wavelength through the Medium interface — so the
// dispersion enters simply by placing Sellmeier glasses (N_BK7, FusedSilica,
// N_SF5, …) in the gaps of the system.
//
// Definitions used here (both in mm):
//
//   - axial (longitudinal) chromatic aberration: the shift of the paraxial
//     axial focus between two wavelengths, FocusAt(λA) − FocusAt(λB). With the
//     classical Fraunhofer lines the F−C difference is negative for a positive
//     singlet (blue focuses closer than red);
//   - lateral (magnification) chromatic aberration: the change of the
//     paraxial image height of a field between two wavelengths, so the chief
//     ray of the two colours lands at different heights in the image plane.
//
// The classical achromat condition follows from these definitions: a cemented
// doublet satisfies Δf = 0 between F and C when φ₁/V₁ + φ₂/V₂ = 0, and its
// remaining chromatic focus shift is the secondary spectrum.

// FocusAt returns the paraxial axial image position (mm in the system frame)
// of an object at axial position objectZ (mm, math.Inf(-1) for infinity) at
// wavelength lambdaNm. The system must be refractive (no mirrors).
func (s *System) FocusAt(objectZ, lambdaNm float64) (float64, error) {
	p, err := s.Paraxial(lambdaNm, objectZ)
	if err != nil {
		return 0, err
	}
	return p.ImageZ, nil
}

// AxialChromatic returns the axial (longitudinal) chromatic aberration
// between the wavelengths lambdaB and lambdaA, in mm:
//
//	Δz = FocusAt(λA) − FocusAt(λB).
//
// With lambdaA = F (486 nm, blue) and lambdaB = C (656 nm, red) the result is
// the classical F−C focus shift; it is negative when the blue focus lies in
// front of the red one (positive singlet).
func (s *System) AxialChromatic(objectZ, lambdaA, lambdaB float64) (float64, error) {
	zA, err := s.FocusAt(objectZ, lambdaA)
	if err != nil {
		return 0, err
	}
	zB, err := s.FocusAt(objectZ, lambdaB)
	if err != nil {
		return 0, err
	}
	return zA - zB, nil
}

// ChromaticFocalShift returns the axial focus position at each of the given
// wavelengths, expressed relative to the focus at lambdas[0] (mm). It is the
// data of the classical chromatic focal-shift curve; lambdas[0] is the
// reference wavelength (e.g. the d line) and appears as 0.
func (s *System) ChromaticFocalShift(objectZ float64, lambdas []float64) ([]float64, error) {
	if len(lambdas) == 0 {
		return nil, fmt.Errorf("geo: chromatic focal shift: no wavelengths")
	}
	out := make([]float64, len(lambdas))
	z0, err := s.FocusAt(objectZ, lambdas[0])
	if err != nil {
		return nil, err
	}
	for i, l := range lambdas {
		z, err := s.FocusAt(objectZ, l)
		if err != nil {
			return nil, err
		}
		out[i] = z - z0
	}
	return out, nil
}

// LateralChromatic returns the lateral (magnification) chromatic aberration of
// the field f between the wavelengths lambdaB and lambdaA, in mm: the
// difference of the paraxial chief-ray image height at the image plane
// z = imageZ,
//
//	Δh = h(λA) − h(λB),
//
// evaluated in the meridional (y) direction. A positive value means the λA
// chief ray lands farther from the axis than the λB one.
func (s *System) LateralChromatic(f Field, imageZ, lambdaA, lambdaB float64) (float64, error) {
	hA, err := s.ParaxialImageHeight(f, imageZ, lambdaA)
	if err != nil {
		return 0, err
	}
	hB, err := s.ParaxialImageHeight(f, imageZ, lambdaB)
	if err != nil {
		return 0, err
	}
	return hA - hB, nil
}

// DispersiveSystem reports whether any medium of the system varies with
// wavelength. It is a quick diagnostic used by the chromatic routines to tell
// an achromatic (all media constant) system from a dispersive one.
func (s *System) DispersiveSystem() bool {
	for _, m := range s.Mediums {
		if m == nil {
			continue
		}
		if m.N(WavelengthF) != m.N(WavelengthC) {
			return true
		}
	}
	return false
}

// RefocusZ returns the axial position of the wavelength-dependent best focus
// of an infinite object measured from the paraxial image plane of the
// reference wavelength lambdaRef, in mm (a convenience for secondary-spectrum
// reports: RefocusZ(objectZ, λ) = FocusAt(objectZ, λ) − FocusAt(objectZ,
// λRef)).
func (s *System) RefocusZ(objectZ, lambdaRef, lambdaNm float64) (float64, error) {
	zr, err := s.FocusAt(objectZ, lambdaRef)
	if err != nil {
		return 0, err
	}
	z, err := s.FocusAt(objectZ, lambdaNm)
	if err != nil {
		return 0, err
	}
	return z - zr, nil
}

// glassIndexPair documents the dispersion inputs of a glass at the F and C
// lines (used by tests and reports).
func glassIndexPair(m Medium) (nF, nC float64) {
	return m.N(WavelengthF), m.N(WavelengthC)
}

// abbeOf returns the Abbe number Vd = (nd − 1)/(nF − nC) of a medium.
func abbeOf(m Medium) float64 {
	nd := m.N(Wavelengthd)
	den := m.N(WavelengthF) - m.N(WavelengthC)
	if den == 0 {
		return math.Inf(1)
	}
	return (nd - 1) / den
}
