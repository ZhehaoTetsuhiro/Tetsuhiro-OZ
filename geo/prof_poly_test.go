package geo

import (
	"fmt"
	"math"
	"math/rand"
	"testing"
)

// ---------------------------------------------------------------------------
// Shared helpers for the general-asphere tests (prof_poly_test.go and
// prof_qasphere_test.go).
// ---------------------------------------------------------------------------

// aspRadii is the radial sample grid used by the pointwise comparisons. All
// points lie inside the domain of the test profiles (the conic rim and the
// normalisation radii used there).
var aspRadii = []float64{0.01, 0.2, 0.7, 1.5, 3, 5.5, 8}

// aspCentralSlope returns the central-difference derivative of Sag at r.
func aspCentralSlope(p Profile, r, h float64) float64 {
	return (p.Sag(r+h) - p.Sag(r-h)) / (2 * h)
}

// checkSlopeCD fails when the analytic Slope disagrees with the central
// difference of Sag by more than 1e-7·max(1, |central|) at any sample radius.
// The mixed relative/absolute bound keeps the test meaningful near zero slope.
func checkSlopeCD(t *testing.T, name string, p Profile, radii []float64) {
	t.Helper()
	const h = 1e-6
	for _, r := range radii {
		num := aspCentralSlope(p, r, h)
		got := p.Slope(r)
		if math.IsNaN(num) || math.IsNaN(got) {
			t.Errorf("%s: NaN at r=%g (slope %g, central %g)", name, r, got, num)
			continue
		}
		tol := 1e-7 * math.Max(1, math.Abs(num))
		if math.Abs(got-num) > tol {
			t.Errorf("%s: Slope(%g) = %.12g, central difference %.12g (|diff| %.3g > %.3g)",
				name, r, got, num, math.Abs(got-num), tol)
		}
	}
}

// aspRandomRays returns n deterministic rays (fixed seed) that start on the
// negative-z side and point roughly along +z, so they hit a surface near the
// vertex. Origins lie within ±halfSpan mm transversely.
func aspRandomRays(n int, halfSpan float64, seed int64) []struct{ o, d Vec3 } {
	rng := rand.New(rand.NewSource(seed))
	out := make([]struct{ o, d Vec3 }, 0, n)
	for i := 0; i < n; i++ {
		o := Vec3{
			X: (2*rng.Float64() - 1) * halfSpan,
			Y: (2*rng.Float64() - 1) * halfSpan,
			Z: -20,
		}
		// small random tilt (≤ 0.35 rad) of the +z direction
		th := 0.35 * rng.Float64()
		ph := 2 * math.Pi * rng.Float64()
		d := Vec3{X: math.Sin(th) * math.Cos(ph), Y: math.Sin(th) * math.Sin(ph), Z: math.Cos(th)}
		out = append(out, struct{ o, d Vec3 }{o, d})
	}
	return out
}

// checkIntersectResidual checks that every ray hitting the profile does so with
// a sag residual |z − Sag(r)| ≤ 1e-12 mm and that the hit parameter respects tMin.
// It returns the number of hits.
func checkIntersectResidual(t *testing.T, name string, p Profile, rays []struct{ o, d Vec3 }) int {
	t.Helper()
	hits := 0
	worst := 0.0
	for _, rv := range rays {
		tt, ok := p.Intersect(rv.o, rv.d, 0, math.Inf(1))
		if !ok {
			continue
		}
		hits++
		if tt < 0 {
			t.Errorf("%s: intersection t=%g < 0", name, tt)
		}
		q := rv.o.Add(rv.d.Mul(tt))
		r := math.Hypot(q.X, q.Y)
		res := math.Abs(q.Z - p.Sag(r))
		if math.IsNaN(res) {
			t.Errorf("%s: NaN residual at r=%g", name, r)
			continue
		}
		if res > worst {
			worst = res
		}
		if res > 1e-12 {
			t.Errorf("%s: Newton residual %.3g mm > 1e-12 at r=%g", name, res, r)
		}
	}
	if hits == 0 {
		t.Errorf("%s: no ray hit the surface (test geometry broken)", name)
	}
	t.Logf("%s: %d hits, worst residual %.3g mm", name, hits, worst)
	return hits
}

// ---------------------------------------------------------------------------
// PolyAsphere
// ---------------------------------------------------------------------------

// TestPolyAsphereEvenMatchesEvenAsphere is acceptance (c): with only even powers
// the general polynomial asphere is the existing EvenAsphere, pointwise to
// 1e-13 in sag and slope.
func TestPolyAsphereEvenMatchesEvenAsphere(t *testing.T) {
	cases := []struct {
		R, K   float64
		coeffs []float64 // A4, A6, A8, ...
	}{
		{R: 25, K: -0.6, coeffs: []float64{1e-6, 2e-9, -3e-12, 1e-15}},
		{R: -50, K: 0.4, coeffs: []float64{-4e-6, 3e-10, 7e-14}},
		{R: 20, K: -1, coeffs: []float64{1.5e-6, -2e-9}},
		{R: 0, K: 0, coeffs: []float64{4e-6, 0, 1e-13}}, // flat base, pure polynomial
		{R: 33, K: -2.2, coeffs: nil},                   // no aspheric terms
	}
	for ci, c := range cases {
		ev := EvenAsphere{R: c.R, K: c.K, Coeffs: c.coeffs}
		terms := make([]PolyTerm, len(c.coeffs))
		for i, a := range c.coeffs {
			terms[i] = PolyTerm{Power: 2*i + 4, Coeff: a}
		}
		pa := PolyAsphere{R: c.R, K: c.K, Terms: terms}
		for _, r := range aspRadii {
			if d := math.Abs(pa.Sag(r) - ev.Sag(r)); !(d <= 1e-13) {
				t.Errorf("case %d: sag at r=%g: poly %.17g vs even %.17g (|diff| %.3g > 1e-13)",
					ci, r, pa.Sag(r), ev.Sag(r), d)
			}
			if d := math.Abs(pa.Slope(r) - ev.Slope(r)); !(d <= 1e-13) {
				t.Errorf("case %d: slope at r=%g: poly %.17g vs even %.17g (|diff| %.3g > 1e-13)",
					ci, r, pa.Slope(r), ev.Slope(r), d)
			}
		}
		// EvenAsphere reports 1/0 = +Inf for R = 0 (existing behaviour); the
		// polynomial asphere treats R = 0 as a flat base with curvature 0.
		if c.R == 0 {
			checkClose(t, "flat-base curvature", pa.Curvature(), 0, 0)
		} else if d := math.Abs(pa.Curvature() - ev.Curvature()); d != 0 {
			t.Errorf("case %d: curvature differs by %g", ci, d)
		}
		// Intersections: same rays, same hit points. The polynomial hit must meet
		// the acceptance (b) residual of 1e-12 mm on its own. The pointwise
		// difference to the EvenAsphere hit is bounded by 3e-12 mm, because the
		// existing EvenAsphere Newton stops at residual 1e-13·(1+|z|+r) (about
		// 1.8e-12 mm here), which is not changed by this work.
		maxDiff := 0.0
		for ri, rv := range aspRandomRays(300, 6, 11) {
			te, oke := ev.Intersect(rv.o, rv.d, 0, math.Inf(1))
			tp, okp := pa.Intersect(rv.o, rv.d, 0, math.Inf(1))
			if oke != okp {
				t.Errorf("case %d ray %d: hit flags differ (even %v, poly %v)", ci, ri, oke, okp)
				continue
			}
			if oke {
				pe := rv.o.Add(rv.d.Mul(te))
				pp := rv.o.Add(rv.d.Mul(tp))
				if res := math.Abs(pp.Z - pa.Sag(math.Hypot(pp.X, pp.Y))); res > 1e-12 {
					t.Errorf("case %d ray %d: poly residual %.3g mm > 1e-12", ci, ri, res)
				}
				if d := pe.Sub(pp).Norm(); d > maxDiff {
					maxDiff = d
				}
			}
		}
		if maxDiff > 3e-12 {
			t.Errorf("case %d: intersection hit points differ by %.3g mm (> 3e-12)", ci, maxDiff)
		}
	}
}

// TestPolyAsphereOddPowersClosedForm checks the sag and slope of a flat-base
// polynomial with odd powers against the elementary closed forms:
//
//	z = a₁ r + a₃ r³ + a₅ r⁵ ,   z' = a₁ + 3a₃ r² + 5a₅ r⁴ .
func TestPolyAsphereOddPowersClosedForm(t *testing.T) {
	const a1, a3, a5 = 3e-3, -4e-4, 2.5e-6
	p := PolyAsphere{R: 0, Terms: []PolyTerm{{1, a1}, {3, a3}, {5, a5}}}
	for _, r := range aspRadii {
		wantZ := a1*r + a3*r*r*r + a5*math.Pow(r, 5)
		wantS := a1 + 3*a3*r*r + 5*a5*math.Pow(r, 4)
		checkRel(t, "odd sag", p.Sag(r), wantZ, 1e-14)
		checkRel(t, "odd slope", p.Slope(r), wantS, 1e-14)
	}
	// The r¹ term alone is a cone: slope a₁ everywhere, sag linear in r.
	cone := PolyAsphere{Terms: []PolyTerm{{1, 0.05}}}
	checkRel(t, "cone slope at axis", cone.Slope(0), 0.05, 1e-15)
	checkRel(t, "cone sag", cone.Sag(4), 0.2, 1e-15)
	// Odd-only asphere has zero paraxial curvature (no r² term).
	checkClose(t, "odd-only curvature", p.Curvature(), 0, 0)
}

// TestPolyAsphereSlopeCentralDifference is acceptance (a) for the polynomial
// asphere, including odd and sparse powers: the analytic slope agrees with the
// central difference of the sag to 1e-7 (mixed relative/absolute).
func TestPolyAsphereSlopeCentralDifference(t *testing.T) {
	profiles := []struct {
		name string
		p    PolyAsphere
	}{
		{"conic+even+odd", PolyAsphere{R: 30, K: -0.5, Terms: []PolyTerm{
			{2, 1e-4}, {3, 2e-5}, {4, -3e-6}, {5, 4e-7}, {7, -1e-9}}}},
		{"flat cone+odd", PolyAsphere{R: 0, Terms: []PolyTerm{{1, 1e-3}, {3, 2e-4}, {4, 1e-6}}}},
		{"sphere sparse", PolyAsphere{R: -40, K: 0.3, Terms: []PolyTerm{{6, 2e-8}, {9, -1e-12}}}},
	}
	for _, c := range profiles {
		checkSlopeCD(t, "poly "+c.name, c.p, []float64{0.3, 1, 2.5, 4, 6, 8})
	}
}

// TestPolyAsphereIntersectResidual is acceptance (b): the Newton intersection
// of rays with a polynomial asphere (conic seeds and the scan fallback for a
// flat base) converges to a sag residual ≤ 1e-12 mm.
func TestPolyAsphereIntersectResidual(t *testing.T) {
	rays := aspRandomRays(4000, 7, 2026)
	cases := []struct {
		name string
		p    PolyAsphere
	}{
		{"conic base", PolyAsphere{R: 30, K: -0.5, Terms: []PolyTerm{
			{2, 1e-4}, {3, 2e-5}, {4, -3e-6}, {5, 4e-7}, {7, -1e-9}}}},
		{"flat base", PolyAsphere{R: 0, Terms: []PolyTerm{{1, 1e-3}, {3, 2e-4}, {4, 1e-6}}}},
	}
	for _, c := range cases {
		n := checkIntersectResidual(t, "poly intersect "+c.name, c.p, rays)
		if n < len(rays)/4 {
			t.Errorf("poly intersect %s: only %d of %d rays hit (geometry too sparse)", c.name, n, len(rays))
		}
	}
}

// TestPolyAsphereValidate checks the structural error paths.
func TestPolyAsphereValidate(t *testing.T) {
	good := PolyAsphere{R: 20, Terms: []PolyTerm{{2, 1e-4}, {3, 1e-5}}}
	if err := good.Validate(); err != nil {
		t.Fatalf("valid asphere rejected: %v", err)
	}
	bad := []PolyAsphere{
		{R: 20, Terms: []PolyTerm{{0, 1e-4}}},         // power 0
		{R: 20, Terms: []PolyTerm{{-2, 1e-4}}},        // negative power
		{R: 20, Terms: []PolyTerm{{4, math.NaN()}}},   // NaN coefficient
		{R: 20, Terms: []PolyTerm{{4, math.Inf(1)}}},  // infinite coefficient
		{R: math.NaN(), Terms: []PolyTerm{{4, 1e-4}}}, // NaN radius
	}
	for i, p := range bad {
		if err := p.Validate(); err == nil {
			t.Errorf("bad asphere %d accepted: %v", i, p)
		}
	}
}

// TestPolyAsphereParaxialDerivative is acceptance (e) for the polynomial
// asphere. The first-order (paraxial) EFL of a singlet whose front surface is
// the polynomial asphere depends on the coefficients only through the vertex
// curvature c_eff = c + 2·a₂. The analytic derivative
//
//	∂EFL/∂a_i = (∂EFL/∂c_eff) · (∂c_eff/∂a_i) ,   ∂c_eff/∂a₂ = 2, else 0,
//
// with ∂EFL/∂c from the exact ray-transfer matrices (see paraxialEFLdCurv) is
// compared with the central difference of the traced EFL, relative 1e-6. For
// i ≠ 2 the analytic value is exactly zero and the finite difference must be
// zero too: this guards against an accidental dependence on higher-order terms.
func TestPolyAsphereParaxialDerivative(t *testing.T) {
	const lambda = 550.0
	const R = 28.0
	// Representative coefficient magnitudes per power (mm^(1−i)).
	a0 := map[int]float64{1: 1e-3, 2: 1e-4, 3: 1e-5, 4: 1e-6, 5: 1e-7, 6: 1e-8}
	for _, i := range []int{1, 2, 3, 4, 5, 6} {
		i := i
		profile := func(a float64) Profile {
			if i == 2 {
				return PolyAsphere{R: R, K: -0.3, Terms: []PolyTerm{{2, a}}}
			}
			// a₂ = 1e-4 stays on for i ≠ 2 so the surface keeps a paraxial power;
			// the coefficient under test is a_i.
			return PolyAsphere{R: R, K: -0.3, Terms: []PolyTerm{{2, 1e-4}, {i, a}}}
		}
		factor := 0.0
		if i == 2 {
			factor = 2
		}
		checkParaxialDerivative(t, fmt.Sprintf("poly a%d", i), lambda, a0[i], factor, profile)
	}
}

// paraxialEFLdCurv returns the analytic derivative ∂EFL/∂c_j of the paraxial
// effective focal length with respect to the vertex curvature of surface j,
// using the exact ray-transfer matrices of the system (object at infinity):
//
//	EFL = −n_out / C ,   M = A·R_j·B ,   ∂M/∂c_j = A·(∂R_j/∂c_j)·B ,
//	∂R_j/∂c_j = [[0,0],[−(n_out,j − n_in,j), 0]] ,   ∂EFL/∂c_j = n_out·(∂M/∂c_j)_21 / C² .
func paraxialEFLdCurv(t *testing.T, sys *System, lambda float64, j int) float64 {
	t.Helper()
	p, err := sys.Paraxial(lambda, math.Inf(-1))
	if err != nil {
		t.Fatalf("paraxial: %v", err)
	}
	N := len(p.els)
	B := Identity()
	if j > 0 {
		b, err := p.SubMatrix(0, j-1)
		if err != nil {
			t.Fatalf("submatrix: %v", err)
		}
		B = b
	}
	A := p.els[j].transferMatrix()
	if j+1 <= N-1 {
		a, err := p.SubMatrix(j+1, N-1)
		if err != nil {
			t.Fatalf("submatrix: %v", err)
		}
		A = a.Mul(A)
	}
	dR := Mat2{C: -(p.els[j].nA - p.els[j].nB)}
	dM := A.Mul(dR).Mul(B)
	C := p.Matrix().C
	nOut := p.Index[len(p.Index)-1]
	return nOut * dM.C / (C * C)
}

// checkParaxialDerivative compares dEFL/da (analytic chain rule with the
// supplied factor ∂c_eff/∂a) against the central difference of the paraxial EFL
// of the singlet built from profile(a), at the point a0. The analytic value
// ∂EFL/∂c is cross-checked against a finite difference in the curvature
// itself, which tests the matrix derivative independently of the profile.
func checkParaxialDerivative(t *testing.T, name string, lambda, a0, factor float64, profile func(a float64) Profile) {
	t.Helper()
	efl := func(a float64) float64 {
		sys := aspEFLSinglet(t, profile(a))
		p, err := sys.Paraxial(lambda, math.Inf(-1))
		if err != nil {
			t.Fatalf("%s: paraxial: %v", name, err)
		}
		return p.EFL
	}
	sys := aspEFLSinglet(t, profile(a0))
	dEdc := paraxialEFLdCurv(t, sys, lambda, 0)
	analytic := dEdc * factor

	h := 1e-3 * math.Abs(a0)
	if h == 0 {
		h = 1e-9
	}
	fd := (efl(a0+h) - efl(a0-h)) / (2 * h)

	// Sanity: the matrix derivative w.r.t. curvature agrees with a finite
	// difference in the front-surface radius (a Sphere has c = 1/R).
	const dR = 1e-3
	eflR := func(R float64) float64 {
		s := aspEFLSinglet(t, Sphere{R: R})
		p, err := s.Paraxial(lambda, math.Inf(-1))
		if err != nil {
			t.Fatalf("paraxial: %v", err)
		}
		return p.EFL
	}
	R0 := 28.0
	fdCurv := (eflR(R0+dR) - eflR(R0-dR)) / (2 * dR) // dEFL/dR
	anCurv := paraxialEFLdCurv(t, aspEFLSinglet(t, Sphere{R: R0}), lambda, 0)
	// dEFL/dR = dEFL/dc · dc/dR = dEFL/dc · (−1/R²)
	anR := anCurv * (-1 / (R0 * R0))
	if math.Abs(anR-fdCurv) > 1e-6*math.Max(1, math.Abs(fdCurv)) {
		t.Errorf("%s: matrix dEFL/dc disagrees with FD in R: analytic dEFL/dR %.10g, FD %.10g",
			name, anR, fdCurv)
	}

	tol := 1e-6 * math.Max(math.Abs(analytic), math.Abs(fd))
	if tol < 1e-12 {
		tol = 1e-12 // both zero (no dependence): FD must be zero as well
	}
	if math.Abs(analytic-fd) > tol {
		t.Errorf("%s: dEFL/da analytic %.12g, central difference %.12g (|diff| %.3g > %.3g)",
			name, analytic, fd, math.Abs(analytic-fd), tol)
	}
	t.Logf("%s: dEFL/da analytic %.10g, FD %.10g", name, analytic, fd)
}

// aspEFLSinglet builds a thin-to-moderate singlet with the given front profile
// and a sphere at the back, n = 1.5 glass, 4 mm thick, object at infinity. The
// front surface is surface 0, so its curvature index is j = 0.
func aspEFLSinglet(t *testing.T, front Profile) *System {
	t.Helper()
	return mustSystem(t, Air{},
		[]Surface{
			{Profile: front, Aperture: CircleAperture(8), Name: "front"},
			{Profile: Sphere{R: -60}, Aperture: CircleAperture(8), Name: "back"},
		},
		[]float64{4},
		[]Medium{Air{}, ConstantIndex(1.5), Air{}},
	)
}
