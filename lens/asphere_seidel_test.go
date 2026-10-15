package lens

import (
	"math"
	"testing"

	"tetsuhiro-oz/geo"
)

// seidelOf returns the third-order sums of a prescription at the d line for
// the on-axis field and for an off-axis field of 0.05 tangent, at the
// entrance pupil of the prescription.
func seidelOf(t *testing.T, p *Prescription, field geo.Field) *geo.Seidel {
	t.Helper()
	sys, err := p.System()
	if err != nil {
		t.Fatal(err)
	}
	par, err := p.Paraxial(geo.Wavelengthd)
	if err != nil {
		t.Fatal(err)
	}
	sd, err := sys.SeidelCoefficients(geo.Wavelengthd, field, par.EntrancePupilRadius)
	if err != nil {
		t.Fatal(err)
	}
	return sd
}

// checkSeidelSame fails when the five sums of two solutions differ by more
// than tol relative to the largest sum involved.
func checkSeidelSame(t *testing.T, name string, a, b *geo.Seidel, tol float64) {
	t.Helper()
	scale := 0.0
	for _, v := range []float64{a.S1, a.S2, a.S3, a.S4, a.S5, b.S1, b.S2, b.S3, b.S4, b.S5} {
		scale = math.Max(scale, math.Abs(v))
	}
	if scale == 0 {
		t.Fatalf("%s: both solutions vanish — vacuous comparison", name)
	}
	sums := [5][2]float64{{a.S1, b.S1}, {a.S2, b.S2}, {a.S3, b.S3}, {a.S4, b.S4}, {a.S5, b.S5}}
	for i, pair := range sums {
		if d := math.Abs(pair[0] - pair[1]); d > tol*scale {
			t.Errorf("%s: S%d = %g vs %g (|diff| = %.3g > %.3g·%.3g)", name, i+1, pair[0], pair[1], d, tol, scale)
		}
	}
}

// TestPrescriptionAsphereSeidelConsistency ties the prescription layer to the
// Seidel layer across the whole chain Card → geo.Profile → SeidelCoefficients:
// the same front surface expressed as poly_terms, as even-asphere coeffs and
// as a q_asphere with the equivalent single coefficient must give identical
// third-order sums, because all three describe the same third-order surface
// (K = 0, the same r⁴ coefficient). The off-axis field exercises the field
// sums on top of S1.
func TestPrescriptionAsphereSeidelConsistency(t *testing.T) {
	const (
		a4 = -3e-7
		s0 = a4 * 12 * 12 * 12 * 12 // qcon coefficient of the same r⁴ term
	)
	poly := polySinglet(a4)
	even := polySinglet(0)
	even.Cards[0].PolyTerms = nil
	even.Cards[0].Coeffs = []float64{a4}
	q := qSinglet(s0, 0)
	// The three prescriptions describe the same surfaces: check the radii and
	// the built profiles first, so a wiring mistake fails here rather than as
	// a mysterious Seidel difference.
	ps, err := poly.System()
	if err != nil {
		t.Fatal(err)
	}
	es, err := even.System()
	if err != nil {
		t.Fatal(err)
	}
	qs, err := q.System()
	if err != nil {
		t.Fatal(err)
	}
	if _, ok := ps.Surfaces[0].Profile.(geo.PolyAsphere); !ok {
		t.Fatalf("poly card built %T, want geo.PolyAsphere", ps.Surfaces[0].Profile)
	}
	if _, ok := es.Surfaces[0].Profile.(geo.EvenAsphere); !ok {
		t.Fatalf("even card built %T, want geo.EvenAsphere", es.Surfaces[0].Profile)
	}
	if _, ok := qs.Surfaces[0].Profile.(geo.QAsphere); !ok {
		t.Fatalf("q card built %T, want geo.QAsphere", qs.Surfaces[0].Profile)
	}
	for _, field := range []geo.Field{geo.FieldTan(0, 0), geo.FieldTan(0, 0.05)} {
		sdPoly, sdEven, sdQ := seidelOf(t, poly, field), seidelOf(t, even, field), seidelOf(t, q, field)
		if sdPoly.S1 == 0 {
			t.Fatal("the aspheric singlet must have spherical aberration")
		}
		checkSeidelSame(t, "poly_terms vs even coeffs", sdPoly, sdEven, 1e-12)
		checkSeidelSame(t, "q_asphere vs poly_terms", sdQ, sdPoly, 1e-12)
		if sdQ.Wavefront().W040 == 0 {
			t.Fatal("the q asphere must have a ρ⁴ term")
		}
	}
	// Teeth: dropping the r⁴ coefficients (the surfaces back to spheres) moves
	// S1 by far more than the tolerance.
	sphere := polySinglet(0)
	sphere.Cards[0].PolyTerms = nil
	sdA := seidelOf(t, poly, geo.FieldTan(0, 0))
	sdS := seidelOf(t, sphere, geo.FieldTan(0, 0))
	if math.Abs(sdA.S1-sdS.S1) < 1e-3*math.Abs(sdS.S1) {
		t.Fatalf("the r⁴ term is not observable in the prescription chain: S1 %g vs sphere %g", sdA.S1, sdS.S1)
	}
}
