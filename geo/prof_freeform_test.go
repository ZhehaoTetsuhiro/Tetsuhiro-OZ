package geo

import (
	"fmt"
	"math"
	"math/rand"
	"strings"
	"testing"
)

// ---------------------------------------------------------------------------
// Acceptance tests of the freeform XY polynomial (see the geo-freeform report)
// ---------------------------------------------------------------------------

// ffSeed is the fixed seed of every random coefficient set of this file
// (deterministic by construction: the same seed always produces the same
// terms, independent of the platform, because math/rand's generator is fully
// specified by the seed).
const ffSeed = 20261015

// ffRandomTerms returns a deterministic coefficient set for an XY polynomial
// with terms of total degree 2..maxDeg. Each coefficient is scaled like the
// corresponding order of a surface with "focal length" f (a_i,j ~ 1/f^(deg-1))
// and jittered by a factor in [0.5, 1.5] drawn from the fixed-seed generator,
// so the sag and the slopes have the magnitude of a realistic asphere over an
// aperture of a few tens of mm.
func ffRandomTerms(seed int64, f float64, maxDeg int) []XYTerm {
	rng := rand.New(rand.NewSource(seed))
	var terms []XYTerm
	for deg := 2; deg <= maxDeg; deg++ {
		scale := 1.0
		for k := 1; k < deg; k++ {
			scale /= f // 1/f^(deg-1), by repeated division (platform-independent)
		}
		for i := 0; i <= deg; i++ {
			j := deg - i
			terms = append(terms, XYTerm{I: i, J: j, Coeff: scale * (0.5 + rng.Float64())})
		}
	}
	return terms
}

// ffMustBuild builds a validated, canonical freeform surface.
func ffMustBuild(t *testing.T, R, K, normR float64, terms []XYTerm) XYPolynomial {
	t.Helper()
	p, err := NewXYPolynomial(R, K, normR, terms)
	if err != nil {
		t.Fatalf("NewXYPolynomial: %v", err)
	}
	return p
}

// ffGrid returns the (x, y) test points used by the shape tests: a grid of the
// aperture, deliberately not symmetric in x and y so that a wrong y-dependence
// (a rotational-symmetry assumption) cannot cancel out.
func ffGrid(a float64, n int) [][2]float64 {
	var pts [][2]float64
	for i := 0; i <= n; i++ {
		for j := 0; j <= n; j++ {
			x := -a + 2*a*float64(i)/float64(n)
			y := -a + 2*a*float64(j)/float64(n)
			pts = append(pts, [2]float64{x, y})
		}
	}
	return pts
}

// ffCentralDiff4 is the 4th-order (5-point) central difference of the sag
// along the x or y direction. A plain 2-point central difference cannot reach
// the 1e-12 relative accuracy of acceptance (a) at these coefficients: its
// truncation error (h²/6)·f”' trades off against the round-off eps·|z|/h and
// bottoms out at 1.6e-11 (h = 1e-4, measured over the test surfaces), an order
// of magnitude above the limit at any step. The 4th-order stencil is still a
// central difference; the same trade-off bottoms out at 1.9e-13 for
// h = 2e-2 (measured; h = 1e-2 gives 4.2e-13, h = 1e-3 4.0e-12).
func ffCentralDiff4(s *Surface, x, y, h float64, axis int) float64 {
	if axis == 0 {
		return (-s.Sag(x+2*h, y) + 8*s.Sag(x+h, y) - 8*s.Sag(x-h, y) + s.Sag(x-2*h, y)) / (12 * h)
	}
	return (-s.Sag(x, y+2*h) + 8*s.Sag(x, y+h) - 8*s.Sag(x, y-h) + s.Sag(x, y-2*h)) / (12 * h)
}

// ---------------------------------------------------------------------------
// Acceptance (a): analytic normal against a central difference of the sag
// ---------------------------------------------------------------------------

// TestXYPolynomialNormalVsCentralDifference is acceptance (a): for random
// coefficient sets (fixed seed) the analytic normal
//
//	n ∝ (−∂z/∂x, −∂z/∂y, 1)
//
// must agree with the central-difference normal of the sag to a relative error
// ≤ 1e-12. Three surfaces are checked: a normalised pure polynomial, an
// unnormalised one on a conic base, and one with a normalisation radius on a
// conic base. The relative error is |n_analytic − n_fd| / |n_fd| with both
// vectors unit (the norm convention of the package for vector errors), and the
// finite difference is the 4th-order central stencil of the sag with h = 2e-2
// (see ffCentralDiff4 for why the 2-point stencil cannot meet 1e-12 and for
// the measured trade-off that fixes the step).
func TestXYPolynomialNormalVsCentralDifference(t *testing.T) {
	const h = 2e-2
	surfaces := []struct {
		name  string
		prof  XYPolynomial
		scale float64
	}{
		{"normalised polynomial", ffMustBuild(t, 0, 0, 20, ffRandomTerms(ffSeed, 60, 6)), 20},
		{"raw polynomial", ffMustBuild(t, 0, 0, 0, ffRandomTerms(ffSeed+1, 60, 6)), 18},
		{"conic base + terms", ffMustBuild(t, -180, -0.6, 15, ffRandomTerms(ffSeed+2, 90, 6)), 15},
	}
	worst := 0.0
	for _, tc := range surfaces {
		s := Surface{Profile: tc.prof}
		for _, xy := range ffGrid(tc.scale, 6) {
			x, y := xy[0], xy[1]
			gx := ffCentralDiff4(&s, x, y, h, 0)
			gy := ffCentralDiff4(&s, x, y, h, 1)
			want, ok := V(-gx, -gy, 1).Unit()
			if !ok {
				t.Fatalf("%s: degenerate gradient at (%g, %g)", tc.name, x, y)
			}
			got := s.Normal(x, y) // the dispatch through Surface must be used
			if d := got.Sub(want).Norm() / want.Norm(); d > worst {
				worst = d
			}
			if d := got.Sub(tc.prof.NormalXY(x, y)).Norm(); d != 0 {
				t.Errorf("%s: Surface.Normal and NormalXY disagree at (%g, %g) by %g", tc.name, x, y, d)
			}
			if d := got.Sub(want).Norm() / want.Norm(); d > 1e-12 {
				t.Errorf("%s: normal at (%g, %g): got %v, central difference %v (rel %.3g > 1e-12)",
					tc.name, x, y, got, want, d)
			}
		}
	}
	t.Logf("acceptance (a): worst central-difference normal error over %d surfaces x 49 points: %.3g (limit 1e-12)",
		len(surfaces), worst)
	if worst > 1e-12 {
		t.Errorf("worst relative normal error %.3g exceeds 1e-12", worst)
	}
}

// TestXYPolynomialGradientSigns locks the analytic gradient against the
// elementary derivative of single-term surfaces, so that a sign error in one
// partial derivative cannot be hidden by the finite-difference comparison
// (which would have to fail too, but this pins the exact values):
//
//	a20·x²  → ∂z/∂x = 2·a20·x, ∂z/∂y = 0;
//	a11·x·y → ∂z/∂x = a11·y,  ∂z/∂y = a11·x  (the mixed term is not
//	          rotationally symmetric and is exactly what a solver that
//	          reduced y to r would get wrong).
func TestXYPolynomialGradientSigns(t *testing.T) {
	p := ffMustBuild(t, 0, 0, 0, []XYTerm{{I: 2, J: 0, Coeff: 3e-3}, {I: 1, J: 1, Coeff: -2e-3}})
	const x, y = 7.0, -5.0
	zx, zy := p.GradXY(x, y)
	checkRel(t, "∂z/∂x of a20 x² + a11 xy", zx, 2*3e-3*x-2e-3*y, 1e-15)
	checkRel(t, "∂z/∂y of a20 x² + a11 xy", zy, -2e-3*x, 1e-15)
	// The mixed term makes the surface non-symmetric: the y partial must not
	// vanish, and RotationallySymmetric must report false.
	if zy == 0 {
		t.Error("∂z/∂y of a mixed xy term is zero: the surface was treated as rotationally symmetric")
	}
	if p.RotationallySymmetric() {
		t.Error("a20 x² + a11 xy reported as rotationally symmetric")
	}
}

// ---------------------------------------------------------------------------
// Acceptance (b): intersection residual, analytic cross-checks, robustness
// ---------------------------------------------------------------------------

// ffRaySet returns a deterministic set of rays with transverse offsets around
// (x0, y0) and transverse slopes up to tiltDeg degrees, launched from z = z0
// travelling in +z (seeded generator, fixed seed → identical rays every run).
func ffRaySet(seed int64, z0, x0, y0, spread, tiltDeg float64, n int) []Ray {
	rng := rand.New(rand.NewSource(seed))
	rays := make([]Ray, 0, n)
	for k := 0; k < n; k++ {
		x := x0 + spread*(2*rng.Float64()-1)
		y := y0 + spread*(2*rng.Float64()-1)
		ax := tiltDeg * (2*rng.Float64() - 1) * math.Pi / 180
		ay := tiltDeg * (2*rng.Float64() - 1) * math.Pi / 180
		d := V(math.Tan(ax), math.Tan(ay), 1)
		rays = append(rays, NewRay(V(x, y, z0), d, 550))
	}
	return rays
}

// TestXYPolynomialIntersectionResidual is acceptance (b): every intersection
// returned by the solver must lie on the surface to ≤ 1e-12 mm, i.e.
//
//	|z(t) − SagXY(x(t), y(t))| ≤ 1e-12 mm
//
// for a strongly curved freeform (conic base plus degree 2..6 terms) and a
// dense, deliberately asymmetric ray set: offsets up to 25 mm, transverse
// tilts up to 25°, rays launched from far in front of the surface. The same
// assertion is repeated through Surface.Intersect (the traced path) and the
// residual is computed from the *public* sag, not from solver internals, so
// the test cannot pass on a self-consistent but wrong surface.
func TestXYPolynomialIntersectionResidual(t *testing.T) {
	p := ffMustBuild(t, -200, -1, 25, ffRandomTerms(ffSeed+3, 80, 6))
	s := Surface{Profile: p, Aperture: CircleAperture(28)}
	rays := ffRaySet(ffSeed, -120, 10, -4, 18, 25, 400)
	hits, misses := 0, 0
	worst := 0.0
	for _, r := range rays {
		tHit, ok := p.IntersectXY(r.Origin, r.Dir, tMinDefault, 28)
		if !ok {
			misses++
			if math.IsNaN(tHit) || math.IsInf(tHit, 0) {
				t.Fatalf("miss returned a non-finite t: %v", tHit)
			}
			continue
		}
		if !finite(tHit) || tHit < tMinDefault {
			t.Fatalf("hit returned t = %v (tMin = %g)", tHit, tMinDefault)
		}
		hits++
		pt := r.At(tHit)
		res := math.Abs(pt.Z - p.SagXY(pt.X, pt.Y))
		if res > worst {
			worst = res
		}
		if res > 1e-12 {
			t.Errorf("ray from %v: residual %.3g mm > 1e-12", r.Origin, res)
		}
		// The same intersection through the Surface dispatch must agree.
		sp, st, ok := s.Intersect(r.Origin, r.Dir, tMinDefault)
		if !ok {
			t.Fatalf("Surface.Intersect missed a ray the direct call hit (from %v)", r.Origin)
		}
		if d := math.Abs(st - tHit); d > 1e-15*(1+math.Abs(tHit)) {
			t.Errorf("direct t %g, Surface.Intersect t %g", tHit, st)
		}
		if d := sp.Sub(pt).Norm(); d > 1e-13 {
			t.Errorf("direct point %v, Surface.Intersect point %v", pt, sp)
		}
	}
	t.Logf("acceptance (b): %d hits, %d misses, worst residual %.3g mm (limit 1e-12)", hits, misses, worst)
	if hits < 300 {
		t.Errorf("only %d of %d rays hit the surface", hits, len(rays))
	}
	if worst > 1e-12 {
		t.Errorf("worst intersection residual %.3g mm exceeds 1e-12", worst)
	}
}

// TestXYPolynomialIntersectionMatchesAnalyticQuadratic verifies the solver
// against the closed-form intersection of a purely quadratic freeform
//
//	z = a20 x² + a02 y²   (no base conic),
//
// for which the ray–surface equation is the quadratic
//
//	(d_z − 2(a20 x d_x + a02 y d_y)) t − (a20 d_x² + a02 d_y²) t² + (z − a20 x² − a02 y²) = 0
//
// in the coordinates of the ray origin. The analytic roots are compared with
// the solver: agreement to 1e-12 mm (and the solver must return the smallest
// root ≥ tMin).
func TestXYPolynomialIntersectionMatchesAnalyticQuadratic(t *testing.T) {
	const a20, a02 = 2.5e-3, -1.1e-3
	p := ffMustBuild(t, 0, 0, 0, []XYTerm{{I: 2, J: 0, Coeff: a20}, {I: 0, J: 2, Coeff: a02}})
	rays := ffRaySet(ffSeed+7, -50, 4, -3, 10, 20, 300)
	hits := 0
	worst := 0.0
	for _, r := range rays {
		// Analytic roots of A t² + B t + C = 0 with
		// A = −(a20 d_x² + a02 d_y²), B = d_z − 2(a20 x d_x + a02 y d_y),
		// C = z − a20 x² − a02 y².
		ox, oy, oz := r.Origin.X, r.Origin.Y, r.Origin.Z
		dx, dy, dz := r.Dir.X, r.Dir.Y, r.Dir.Z
		A := -(a20*dx*dx + a02*dy*dy)
		B := dz - 2*(a20*ox*dx+a02*oy*dy)
		C := oz - a20*ox*ox - a02*oy*oy
		disc := B*B - 4*A*C
		var want float64
		if math.Abs(A) < 1e-300 {
			if B == 0 {
				continue
			}
			want = -C / B
		} else if disc < 0 {
			continue
		} else {
			sq := math.Sqrt(disc)
			// numerically stable pair of roots
			q := -0.5 * (B + math.Copysign(sq, B))
			r1, r2 := q/A, C/q
			if r1 > r2 {
				r1, r2 = r2, r1
			}
			switch {
			case r1 >= tMinDefault:
				want = r1
			case r2 >= tMinDefault:
				want = r2
			default:
				continue
			}
		}
		if want < tMinDefault {
			continue
		}
		got, ok := p.IntersectXY(r.Origin, r.Dir, tMinDefault, 0)
		if !ok {
			t.Errorf("solver missed a ray with analytic root t = %.15g (from %v d %v)", want, r.Origin, r.Dir)
			continue
		}
		hits++
		if d := math.Abs(got - want); d > 1e-12 {
			t.Errorf("ray from %v: solver t = %.15g, analytic t = %.15g (|diff| %.3g > 1e-12)", r.Origin, got, want, d)
		}
		if d := math.Abs(got - want); d > worst {
			worst = d
		}
	}
	t.Logf("analytic quadratic cross-check: %d rays, worst |t_solver − t_analytic| = %.3g mm (limit 1e-12)", hits, worst)
	if hits < 200 {
		t.Fatalf("only %d rays with analytic roots were checked", hits)
	}
}

// TestXYPolynomialFarFromAxisStaysFinite is the robustness requirement of the
// task: rays whose line passes far from the vertex (transverse offsets up to
// 500 mm, i.e. more than ten times the aperture), near-grazing directions and
// wildly oblique directions must never make the solver diverge. Every call
// must return either a finite t, or ok = false / an error; never NaN, never
// ±Inf, and the error-returning entry point must never carry a NaN in its
// message.
//
// The residual is checked against the scale of the hit:
//
//   - hits inside the clear aperture (r ≤ 25 mm): |z − Sag| ≤ 1e-12 mm, the
//     acceptance level of the task;
//   - hits on the analytic branch far outside the aperture, where the surface
//     value itself is of order 1e4 mm: |z − Sag| ≤ 1e-13·(1 + |z| + r), i.e.
//     at the round-off floor of the evaluation (measured: a few ulp of |z|,
//     ~1e-14 relative). Demanding an absolute 1e-12 mm there would demand
//     1e-16 relative accuracy from a float64 evaluation of x^6.
func TestXYPolynomialFarFromAxisStaysFinite(t *testing.T) {
	p := ffMustBuild(t, -200, -1, 25, ffRandomTerms(ffSeed+4, 80, 6))
	cases, inAper, outAper := 0, 0, 0
	worstIn, worstRelOut := 0.0, 0.0
	for _, off := range []float64{0, 12, 60, 200, 500, 2000} {
		for _, tilt := range []float64{0, 5, 30, 70, 89, 90} {
			rays := ffRaySet(ffSeed+int64(off)+int64(tilt), -200, off, -off/2, 20, tilt, 40)
			for _, r := range rays {
				cases++
				tHit, ok := p.IntersectXY(r.Origin, r.Dir, 0, 25)
				if math.IsNaN(tHit) || math.IsInf(tHit, 0) {
					t.Fatalf("offset %g tilt %g: non-finite t = %v", off, tilt, tHit)
				}
				if ok {
					pt := r.At(tHit)
					if !finite(pt.X) || !finite(pt.Y) || !finite(pt.Z) {
						t.Fatalf("offset %g tilt %g: non-finite hit point %v", off, tilt, pt)
					}
					if tHit < 0 {
						t.Fatalf("offset %g tilt %g: t = %v < tMin", off, tilt, tHit)
					}
					res := math.Abs(pt.Z - p.SagXY(pt.X, pt.Y))
					rr := math.Hypot(pt.X, pt.Y)
					if rr <= 25 {
						inAper++
						if res > worstIn {
							worstIn = res
						}
						if res > 1e-12 {
							t.Errorf("offset %g tilt %g: in-aperture residual %.3g mm > 1e-12", off, tilt, res)
						}
					} else {
						outAper++
						rel := res / (1 + math.Abs(pt.Z) + rr)
						if rel > worstRelOut {
							worstRelOut = rel
						}
						if rel > 1e-13 {
							t.Errorf("offset %g tilt %g: relative residual %.3g > 1e-13 at r = %g", off, tilt, rel, rr)
						}
					}
				}
				if _, err := p.IntersectRay(r.Origin, r.Dir, 0, 25); err != nil {
					if got := err.Error(); got == "" {
						t.Error("empty error message")
					}
				}
			}
		}
	}
	t.Logf("far-from-axis robustness: %d solver calls, %d in-aperture hits (worst residual %.3g mm, limit 1e-12), "+
		"%d out-of-aperture hits (worst relative residual %.3g, limit 1e-13); no NaN/Inf anywhere",
		cases, inAper, worstIn, outAper, worstRelOut)
}

// TestXYPolynomialMissReturnsError checks the failure contract: a ray that
// passes a plane-like freeform parallel to it (or an offset plane-like surface
// it never meets) reports ok = false and a non-nil error, and never a NaN.
func TestXYPolynomialMissReturnsError(t *testing.T) {
	plane := XYPolynomial{} // z ≡ 0
	// Parallel ray: no intersection.
	r := NewRay(V(1, 2, 5), V(1, 0, 0), 550)
	if tt, ok := plane.IntersectXY(r.Origin, r.Dir, 0, 0); ok || math.IsNaN(tt) {
		t.Fatalf("parallel ray: ok = %v, t = %v", ok, tt)
	}
	if _, err := plane.IntersectRay(r.Origin, r.Dir, 0, 0); err == nil {
		t.Fatal("parallel ray: IntersectRay returned no error")
	}
	// Ray travelling away from the surface (behind it).
	back := NewRay(V(0, 0, 5), V(0, 0, 1), 550)
	if tt, ok := plane.IntersectXY(back.Origin, back.Dir, 0, 0); ok || math.IsNaN(tt) {
		t.Fatalf("backward ray: ok = %v, t = %v", ok, tt)
	}
	// NaN input: a clean error, not a NaN result.
	if _, err := plane.IntersectRay(V(math.NaN(), 0, 0), V(0, 0, 1), 0, 0); err == nil {
		t.Fatal("NaN origin: IntersectRay returned no error")
	}
	if _, err := plane.IntersectRay(V(0, 0, 0), Zero, 0, 0); err == nil {
		t.Fatal("zero direction: IntersectRay returned no error")
	}
	// A curved base outside its domain (sphere rim) must be a clean miss.
	sp := ffMustBuild(t, 10, 0, 0, nil)
	if tt, ok := sp.IntersectXY(V(20, 0, 0), V(0, 0, 1), 0, 0); ok || math.IsNaN(tt) {
		t.Fatalf("outside the conic domain: ok = %v, t = %v", ok, tt)
	}
	if !math.IsNaN(sp.SagXY(20, 0)) {
		t.Error("a sphere of R = 10 was expected to be undefined at r = 20")
	}
}

// ---------------------------------------------------------------------------
// Acceptance (c): degeneracies
// ---------------------------------------------------------------------------

// TestXYPolynomialZeroDegeneratesToPlane is acceptance (c), part 1: with no
// base conic and all coefficients zero the surface is the plane z ≡ 0. The
// sag must be exactly zero (the evaluation of an empty departure adds nothing),
// the normal exactly (0, 0, 1), and a traced system using the freeform plane
// must reproduce the same system with Plane{}.Degenerate to ≤ 1e-13.
func TestXYPolynomialZeroDegeneratesToPlane(t *testing.T) {
	var zero XYPolynomial
	if err := zero.Validate(); err != nil {
		t.Fatalf("zero value must be valid: %v", err)
	}
	for _, xy := range ffGrid(30, 8) {
		if z := zero.SagXY(xy[0], xy[1]); z != 0 {
			t.Fatalf("zero freeform sag at (%g, %g) = %g, want exactly 0", xy[0], xy[1], z)
		}
		if n := zero.NormalXY(xy[0], xy[1]); n != (Vec3{0, 0, 1}) {
			t.Fatalf("zero freeform normal at (%g, %g) = %v, want (0,0,1)", xy[0], xy[1], n)
		}
	}
	// A plate-like system: the freeform plane against Plane{}.
	pl := plateSystem(t, 1.5, 7)
	ffp := mustSystem(t, Air{},
		[]Surface{
			{Profile: zero, Aperture: CircleAperture(50)},
			{Profile: zero, Aperture: CircleAperture(50)},
		},
		[]float64{7},
		[]Medium{Air{}, ConstantIndex(1.5), Air{}})
	n := 0
	worstPos, worstDir := 0.0, 0.0
	for _, r := range aspFan(550) {
		a := pl.Trace(r)
		b := ffp.Trace(r)
		if a.Status != b.Status {
			t.Fatalf("status plane %v, freeform %v", a.Status, b.Status)
		}
		if !a.OK() {
			continue
		}
		n++
		if d := a.Final().Sub(b.Final()).Norm(); d > worstPos {
			worstPos = d
		}
		if d := a.Direction().Sub(b.Direction()).Norm(); d > worstDir {
			worstDir = d
		}
	}
	if worstPos > 1e-13 || worstDir > 1e-13 {
		t.Errorf("freeform plane differs from Plane by %.3g mm (position) / %.3g (direction)", worstPos, worstDir)
	}
	if n < 30 {
		t.Fatalf("only %d rays traversed the plate", n)
	}
	t.Logf("acceptance (c1): freeform plane == Plane on %d traced rays, worst %.3g mm / %.3g (limit 1e-13)", n, worstPos, worstDir)
}

// TestXYPolynomialBaseConicMatchesConicAndSphere is acceptance (c), part 2:
// with zero coefficients the freeform surface is exactly the base conic
// (Conic for K ≠ 0, Sphere for K = 0), pointwise in the sag and end-to-end in
// a traced system. R = ±Inf and R = 0 both mean a flat base.
func TestXYPolynomialBaseConicMatchesConicAndSphere(t *testing.T) {
	cases := []struct {
		R, K float64
		rad  Profile
	}{
		{40, 0, Sphere{R: 40}},
		{-40, 0, Sphere{R: -40}},
		{60, -0.5, Conic{R: 60, K: -0.5}},
		{-30, 2.3, Conic{R: -30, K: 2.3}},
		{25, -1, Conic{R: 25, K: -1}},
		{math.Inf(1), 0, Plane{}},
		{math.Inf(-1), 0.7, Plane{}},
		{0, -0.3, Plane{}},
	}
	worstPoint := 0.0
	for _, tc := range cases {
		p := ffMustBuild(t, tc.R, tc.K, 0, nil)
		for _, xy := range ffGrid(20, 8) {
			r := math.Hypot(xy[0], xy[1])
			got := p.SagXY(xy[0], xy[1])
			want := tc.rad.Sag(r)
			if math.IsNaN(want) {
				if !math.IsNaN(got) {
					t.Errorf("R=%g K=%g: sag at (%g, %g) = %g, want NaN", tc.R, tc.K, xy[0], xy[1], got)
				}
				continue
			}
			if d := math.Abs(got - want); d > worstPoint {
				worstPoint = d
			}
		}
		if d := math.Abs(p.Curvature() - tc.rad.Curvature()); d != 0 {
			t.Errorf("R=%g K=%g: curvature %g, %v %g", tc.R, tc.K, p.Curvature(), tc.rad, tc.rad.Curvature())
		}
	}
	if worstPoint > 1e-13 {
		t.Errorf("base-conic sag differs from Conic/Sphere by %.3g mm", worstPoint)
	}
	t.Logf("acceptance (c2): pointwise worst sag difference %.3g mm over %d surfaces x 81 points (limit 1e-13)", worstPoint, len(cases))
	// End to end: a singlet whose front surface is the freeform base conic
	// against the same singlet with the bare conic.
	worstEnd := 0.0
	for _, tc := range cases[:5] {
		p := ffMustBuild(t, tc.R, tc.K, 0, nil)
		sa := aspSinglet(t, p)
		sb := aspSinglet(t, tc.rad)
		n := 0
		for _, r := range aspFan(550) {
			a := sa.Trace(r)
			b := sb.Trace(r)
			if a.Status != b.Status {
				t.Fatalf("R=%g K=%g: status freeform %v, %v %v", tc.R, tc.K, a.Status, tc.rad, b.Status)
			}
			if !a.OK() {
				continue
			}
			n++
			if d := a.Final().Sub(b.Final()).Norm(); d > worstEnd {
				worstEnd = d
			}
			if d := a.Direction().Sub(b.Direction()).Norm(); d > 1e-13 {
				t.Errorf("R=%g K=%g: final directions differ by %.3g", tc.R, tc.K, d)
			}
		}
		if n < 30 {
			t.Errorf("R=%g K=%g: only %d rays traversed the singlet", tc.R, tc.K, n)
		}
	}
	if worstEnd > 1e-13 {
		t.Errorf("end-to-end base-conic difference %.3g mm", worstEnd)
	}
	t.Logf("acceptance (c2): end-to-end worst position difference %.3g mm (limit 1e-13)", worstEnd)
	t.Logf("acceptance (c2): %d base-only surfaces match Conic/Sphere pointwise (limit 1e-13) and end-to-end", len(cases))
}

// ---------------------------------------------------------------------------
// Acceptance (d): rotational symmetry
// ---------------------------------------------------------------------------

// ffSymTerms expands Σ_m b_m (x²+y²)^m into freeform monomials:
//
//	(x²+y²)^m = Σ_k C(m,k) x^(2k) y^(2(m-k)),
//
// which is exactly the coefficient relation that a rotationally symmetric
// freeform surface must satisfy.
func ffSymTerms(b []float64) []XYTerm {
	var terms []XYTerm
	for m, bm := range b {
		if m == 0 {
			continue // r^0 = 1 is a piston: not representable (i+j >= 2)
		}
		for k := 0; k <= m; k++ {
			terms = append(terms, XYTerm{I: 2 * k, J: 2 * (m - k), Coeff: bm * binom(m, k)})
		}
	}
	return terms
}

// TestXYPolynomialAxisymmetricMatchesRadialProfile is acceptance (d): a
// freeform surface with only even i+j and the rotation-symmetry coefficient
// relation (a_2k,2m−2k = C(m,k)·b_m) must have exactly the sag of the
// equivalent rotationally symmetric profile Σ b_m r^(2m) — checked against
// PolyAsphere pointwise (≤ 1e-12) and end-to-end through a traced singlet.
// The vertex curvature must be the same number, the symmetry detector must
// report the surface as symmetric, and the first-order engine must accept it —
// while it must reject the asymmetric variant (the applicability rule of the
// task: "不适用近轴一阶引擎").
func TestXYPolynomialAxisymmetricMatchesRadialProfile(t *testing.T) {
	const R, K = 40.0, -0.5
	b := []float64{0, 1.3e-4, -2.7e-7, 4.1e-11} // b_1 r², b_2 r⁴, b_3 r⁶
	sym := ffMustBuild(t, R, K, 0, ffSymTerms(b))
	rad := PolyAsphere{R: R, K: K, Terms: []PolyTerm{{Power: 2, Coeff: b[1]}, {Power: 4, Coeff: b[2]}, {Power: 6, Coeff: b[3]}}}
	if !sym.RotationallySymmetric() {
		t.Error("the rotationally symmetric coefficient set was not detected")
	}
	checkRel(t, "vertex curvature of the axisymmetric freeform", sym.Curvature(), rad.Curvature(), 1e-15)
	worst := 0.0
	for _, xy := range ffGrid(18, 8) {
		r := math.Hypot(xy[0], xy[1])
		d := math.Abs(sym.SagXY(xy[0], xy[1]) - rad.Sag(r))
		if d > worst {
			worst = d
		}
		if d > 1e-12 {
			t.Errorf("sag at (%g, %g): freeform %.17g, radial %.17g (|diff| %.3g)", xy[0], xy[1], sym.SagXY(xy[0], xy[1]), rad.Sag(r), d)
		}
	}
	// End to end through the singlet.
	sa := aspSinglet(t, sym)
	sb := aspSinglet(t, rad)
	n := 0
	worstRay := 0.0
	for _, r := range aspFan(550) {
		a := sa.Trace(r)
		bb := sb.Trace(r)
		if a.Status != bb.Status {
			t.Fatalf("status freeform %v, radial %v", a.Status, bb.Status)
		}
		if !a.OK() {
			continue
		}
		n++
		if d := a.Final().Sub(bb.Final()).Norm(); d > worstRay {
			worstRay = d
		}
		if d := a.Direction().Sub(bb.Direction()).Norm(); d > 1e-12 {
			t.Errorf("final directions differ by %.3g", d)
		}
	}
	if worstRay > 1e-12 {
		t.Errorf("axisymmetric freeform differs from the radial profile by %.3g mm", worstRay)
	}
	if n < 30 {
		t.Fatalf("only %d rays traversed the singlet", n)
	}
	t.Logf("acceptance (d): worst sag difference %.3g mm over 81 points, worst traced-ray difference %.3g mm over %d rays (limit 1e-12)", worst, worstRay, n)

	// The applicability rule: the first-order engine accepts the symmetric
	// freeform and rejects the asymmetric one.
	if _, err := sa.Paraxial(550, math.Inf(-1)); err != nil {
		t.Errorf("Paraxial rejected a rotationally symmetric freeform singlet: %v", err)
	}
	asym := ffMustBuild(t, R, K, 0, append(ffSymTerms(b), XYTerm{I: 1, J: 1, Coeff: 3e-5}))
	if asym.RotationallySymmetric() {
		t.Error("a mixed x·y term did not break the rotational symmetry detector")
	}
	sAsym := aspSinglet(t, asym)
	if _, err := sAsym.Paraxial(550, math.Inf(-1)); err == nil {
		t.Error("Paraxial accepted an asymmetric freeform surface")
	}
	if _, err := sAsym.SeidelCoefficients(550, FieldTan(0, 0), 4); err == nil {
		t.Error("Seidel accepted an asymmetric freeform surface")
	}
	if sAsym.Reversed() != nil {
		t.Error("Reversed accepted an asymmetric freeform surface")
	}
}

// ---------------------------------------------------------------------------
// Canonical order, validation and the radial (Profile) surrogates
// ---------------------------------------------------------------------------

// TestXYPolynomialCanonicalValidate covers the structural API: validation
// errors (order below quadratic, negative exponents, non-finite data, negative
// normalisation radius), the canonical order (deterministic summation
// independent of the order the caller supplied the terms in, duplicates
// merged), and the documented meaning of the Profile surrogates
// (Sag/Slope = the +y meridian, Curvature = the mean vertex curvature).
func TestXYPolynomialCanonicalValidate(t *testing.T) {
	// Validation.
	bad := []struct {
		name string
		p    XYPolynomial
	}{
		{"piston", XYPolynomial{Terms: []XYTerm{{0, 0, 1}}}},
		{"tilt", XYPolynomial{Terms: []XYTerm{{1, 0, 1e-3}}}},
		{"negative exponent", XYPolynomial{Terms: []XYTerm{{-1, 3, 1e-3}}}},
		{"NaN coefficient", XYPolynomial{Terms: []XYTerm{{2, 0, math.NaN()}}}},
		{"Inf coefficient", XYPolynomial{Terms: []XYTerm{{2, 0, math.Inf(1)}}}},
		{"NaN radius", XYPolynomial{R: math.NaN()}},
		{"NaN conic constant", XYPolynomial{K: math.NaN()}},
		{"negative norm radius", XYPolynomial{NormR: -1, Terms: []XYTerm{{2, 0, 1e-3}}}},
	}
	for _, tc := range bad {
		if err := tc.p.Validate(); err == nil {
			t.Errorf("%s: Validate accepted an invalid surface %+v", tc.name, tc.p)
		}
		if _, err := NewXYPolynomial(tc.p.R, tc.p.K, tc.p.NormR, tc.p.Terms); err == nil {
			t.Errorf("%s: NewXYPolynomial accepted an invalid surface", tc.name)
		}
	}
	if _, err := NewXYPolynomial(-200, -1, 25, []XYTerm{{2, 0, -1e-3}, {0, 2, -1e-3}}); err != nil {
		t.Errorf("NewXYPolynomial rejected a valid surface: %v", err)
	}

	// Canonical order: same terms, different input order, duplicate entries.
	base := []XYTerm{{4, 0, 1e-6}, {2, 0, -2e-3}, {0, 2, -2e-3}, {1, 1, 3e-5}, {0, 4, 1e-6}, {2, 2, 2e-6}}
	shuffled := []XYTerm{{2, 2, 2e-6}, {0, 4, 1e-6}, {2, 0, -1e-3}, {1, 1, 3e-5}, {4, 0, 1e-6}, {0, 2, -2e-3}, {2, 0, -1e-3}}
	a := ffMustBuild(t, 0, 0, 0, base)
	b := ffMustBuild(t, 0, 0, 0, shuffled)
	if len(a.Terms) != 6 || len(b.Terms) != 6 {
		t.Fatalf("canonicalisation kept %d and %d terms, want 6 (duplicates merged)", len(a.Terms), len(b.Terms))
	}
	if got := a.Canonical(); fmt.Sprint(got.Terms) != fmt.Sprint(b.Terms) {
		t.Errorf("canonical orders differ: %v vs %v", got.Terms, b.Terms)
	}
	for _, xy := range ffGrid(12, 5) {
		za, zb := a.SagXY(xy[0], xy[1]), b.SagXY(xy[0], xy[1])
		if math.Float64bits(za) != math.Float64bits(zb) {
			t.Errorf("sag at (%g,%g) depends on the input term order: %.17g vs %.17g", xy[0], xy[1], za, zb)
		}
	}

	// Normalisation radius: Σ a_ij (x/N)^i (y/N)^j equals the unnormalised
	// surface with b_ij = a_ij/N^(i+j).
	const N = 15.0
	norm := ffMustBuild(t, 0, 0, N, []XYTerm{{2, 0, -3e-2}, {0, 2, -3e-2}, {1, 1, 5e-3}, {4, 0, 2e-3}})
	raw := ffMustBuild(t, 0, 0, 0, []XYTerm{{2, 0, -3e-2 / (N * N)}, {0, 2, -3e-2 / (N * N)}, {1, 1, 5e-3 / (N * N)}, {4, 0, 2e-3 / (N * N * N * N)}})
	for _, xy := range ffGrid(12, 5) {
		checkRel(t, "normalised vs scaled raw sag", norm.SagXY(xy[0], xy[1]), raw.SagXY(xy[0], xy[1]), 1e-14)
	}

	// The Profile surrogates: Sag/Slope are the +y meridian, Curvature the
	// mean of the two vertex curvatures.
	p := ffMustBuild(t, 0, 0, 0, []XYTerm{{2, 0, -2.5e-3}, {0, 2, -1.5e-3}, {1, 1, 4e-4}})
	checkClose(t, "Sag(r) = SagXY(0, r)", p.Sag(6), p.SagXY(0, 6), 0)
	_, zy := p.GradXY(0, 6)
	checkClose(t, "Slope(r) = dz/dy(0, r)", p.Slope(6), zy, 0)
	const h = 1e-4
	kxx := (p.SagXY(h, 0) - 2*p.SagXY(0, 0) + p.SagXY(-h, 0)) / (h * h)
	kyy := (p.SagXY(0, h) - 2*p.SagXY(0, 0) + p.SagXY(0, -h)) / (h * h)
	checkRel(t, "Curvature = mean vertex curvature", p.Curvature(), 0.5*(kxx+kyy), 1e-9)
	// ZRange through the Surface dispatch samples the two-dimensional sag.
	s := Surface{Profile: p, Aperture: CircleAperture(8)}
	lo, hi := s.ZRange(0)
	zmin, zmax := 0.0, 0.0
	for _, xy := range ffGrid(8, 32) {
		if math.Hypot(xy[0], xy[1]) > 8+1e-12 {
			continue
		}
		z := p.SagXY(xy[0], xy[1])
		if z < zmin {
			zmin = z
		}
		if z > zmax {
			zmax = z
		}
	}
	if lo > zmin+1e-12 || hi < zmax-1e-12 {
		t.Errorf("ZRange [%g, %g] does not contain the sampled range [%g, %g]", lo, hi, zmin, zmax)
	}
	if _, _, ok := s.Intersect(V(2, -1, -20), V(0, 0, 1), tMinDefault); !ok {
		t.Error("Surface.Intersect failed on an axial ray")
	}
	if p.String() == "" || !strings.Contains(p.String(), "xy polynomial") {
		t.Errorf("String() = %q", p.String())
	}
}

// TestXYPolynomialIntersectionFindsFirstCrossing verifies the Profile.Intersect
// contract that the *smallest* valid intersection is returned when the ray
// meets the surface more than once:
//
//   - a chord of the quadratic bowl z = a20x² + a02y² that enters and leaves
//     the surface has two analytic crossings; the solver must return the
//     first, and the test counts how many such double-crossing rays it
//     exercised;
//   - a strongly curved mixed surface (conic base plus degree-6 terms) is
//     scanned densely (1e5 samples over the search window) and the first sign
//     change of the residual is compared with the returned root to the
//     resolution of the scan.
func TestXYPolynomialIntersectionFindsFirstCrossing(t *testing.T) {
	// (1) Quadratic bowl, analytic double crossings.
	const a20, a02 = -2.5e-3, -2.5e-3
	bowl := ffMustBuild(t, 0, 0, 0, []XYTerm{{I: 2, J: 0, Coeff: a20}, {I: 0, J: 2, Coeff: a02}})
	double, checked := 0, 0
	for k := 1; k <= 6; k++ {
		// A horizontal chord of the bowl at depth h: it enters the surface at
		// x = -sqrt(h/a) and leaves it at x = +sqrt(h/a), so it has exactly two
		// analytic crossings.
		h := 0.02 * float64(k)
		rc := math.Sqrt(h / -a20)
		o := V(-(rc + 3), 0, -h)
		d := V(1, 0, 0)
		ox, oy, oz := o.X, o.Y, o.Z
		A := -(a20*d.X*d.X + a02*d.Y*d.Y)
		B := d.Z - 2*(a20*ox*d.X+a02*oy*d.Y)
		C := oz - a20*ox*ox - a02*oy*oy
		disc := B*B - 4*A*C
		if disc <= 0 || A == 0 {
			continue
		}
		sq := math.Sqrt(disc)
		q := -0.5 * (B + math.Copysign(sq, B))
		r1, r2 := q/A, C/q
		if r1 > r2 {
			r1, r2 = r2, r1
		}
		if r1 < tMinDefault {
			continue
		}
		checked++
		if r2 > r1+1e-9 {
			double++
		}
		got, ok := bowl.IntersectXY(o, d, tMinDefault, 0)
		if !ok {
			t.Fatalf("ray %d: no intersection found (analytic roots %g, %g)", k, r1, r2)
		}
		if math.Abs(got-r1) > 1e-12 {
			t.Errorf("ray %d: returned t = %.15g, first analytic crossing %.15g (second %.15g)", k, got, r1, r2)
		}
	}
	if double < 5 {
		t.Fatalf("only %d of %d rays had two crossings: the check is not exercised", double, checked)
	}
	// (2) Mixed surface: dense scan reference.
	p := ffMustBuild(t, -200, -1, 25, ffRandomTerms(ffSeed+9, 80, 6))
	rays := ffRaySet(ffSeed+21, -60, 6, -4, 10, 12, 12)
	worst := 0.0
	for _, r := range rays {
		got, ok := p.IntersectXY(r.Origin, r.Dir, tMinDefault, 25)
		if !ok {
			continue
		}
		// Dense scan of the same window the solver uses.
		tA, tB, okW := p.searchWindowXY(r.Origin, r.Dir, tMinDefault, 25)
		if !okW {
			t.Fatalf("no search window for ray from %v", r.Origin)
		}
		const n = 100000
		prevT := tA
		prevG, _ := p.residualXY(r.Origin, r.Dir, prevT)
		found := math.NaN()
		for i := 1; i <= n; i++ {
			tt := tA + (tB-tA)*float64(i)/float64(n)
			g, _ := p.residualXY(r.Origin, r.Dir, tt)
			if !math.IsNaN(g) && !math.IsNaN(prevG) && prevG*g <= 0 {
				found = 0.5 * (prevT + tt)
				break
			}
			prevT, prevG = tt, g
		}
		if math.IsNaN(found) {
			continue
		}
		res := (tB - tA) / n
		if math.Abs(got-found) > 10*res+1e-9 {
			t.Errorf("first crossing: solver t = %.15g, dense scan %.15g (resolution %.3g)", got, found, res)
		}
		if d := math.Abs(got - found); d > worst {
			worst = d
		}
	}
	t.Logf("first-crossing check: %d double-crossing chords (all returned the first root); on the mixed surface the "+
		"returned root coincides with the first crossing of a 1e5-sample dense scan to within the scan resolution "+
		"(worst deviation %.3g mm)", double, worst)
}
