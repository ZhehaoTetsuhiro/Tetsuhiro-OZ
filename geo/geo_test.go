package geo

import (
	"math"
	"testing"
)

// ---------------------------------------------------------------------------
// shared test helpers
// ---------------------------------------------------------------------------

// checkClose fails when |got-want| > tol (absolute).
func checkClose(t *testing.T, name string, got, want, tol float64) {
	t.Helper()
	if math.IsNaN(got) {
		t.Errorf("%s: got NaN, want %.15g", name, want)
		return
	}
	if d := math.Abs(got - want); d > tol {
		t.Errorf("%s: got %.15g, want %.15g (|diff| = %.3g > tol %.3g)", name, got, want, d, tol)
	}
}

// checkRel fails when |got-want| > tol*max(1,|want|).
func checkRel(t *testing.T, name string, got, want, tol float64) {
	t.Helper()
	scale := math.Max(1, math.Abs(want))
	checkClose(t, name, got, want, tol*scale)
}

// checkLess fails when !(got < limit).
func checkLess(t *testing.T, name string, got, limit float64) {
	t.Helper()
	if !(got < limit) {
		t.Errorf("%s: got %.15g, want < %.15g", name, got, limit)
	}
}

// checkStatus fails when got != want.
func checkStatus(t *testing.T, name string, got, want Status) {
	t.Helper()
	if got != want {
		t.Errorf("%s: got status %v, want %v", name, got, want)
	}
}

// mustSystem builds a system and fails the test on validation errors.
func mustSystem(t *testing.T, m Medium, surfaces []Surface, thicknesses []float64, mediums []Medium) *System {
	t.Helper()
	sys, err := NewSystem(m, surfaces, thicknesses, mediums)
	if err != nil {
		t.Fatalf("NewSystem: %v", err)
	}
	return sys
}

// plateSystem returns a plane-parallel plate of thickness T and index n in air,
// entered at z = 0 (surface 0) and left at z = T (surface 1).
func plateSystem(t *testing.T, n, T float64) *System {
	t.Helper()
	return mustSystem(t, Air{},
		[]Surface{
			{Profile: Plane{}, Aperture: CircleAperture(50), Name: "plate front"},
			{Profile: Plane{}, Aperture: CircleAperture(50), Name: "plate rear"},
		},
		[]float64{T},
		[]Medium{Air{}, ConstantIndex(n), Air{}},
	)
}
