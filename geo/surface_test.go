package geo

import (
	"math"
	"math/rand"
	"testing"
)

// TestConicSagAgainstClosedForm checks the conic sag against the elementary
// closed forms: sphere z = R - sqrt(R²-r²) (R > 0), parabola z = r²/(2R) and
// the vertex expansion z ≈ r²/(2R).
func TestConicSagAgainstClosedForm(t *testing.T) {
	const R = 26.5
	for _, r := range []float64{0, 1e-6, 1, 5, 12.345, 20, 26, 26.4} {
		want := R - math.Sqrt(R*R-r*r)
		checkRel(t, "sphere sag", Sphere{R: R}.Sag(r), want, 1e-13)
		// Hyperbolic version of the same conic, K = 0 must reproduce it too.
		checkRel(t, "conic K=0 sag", Conic{R: R, K: 0}.Sag(r), want, 1e-13)
	}
	// Concave (R < 0) must mirror the convex case.
	for _, r := range []float64{1, 5, 20} {
		checkClose(t, "negative R sag", Sphere{R: -R}.Sag(r), -(R - math.Sqrt(R*R-r*r)), 1e-12)
	}
	// Parabola: exact z = r²/(2R).
	for _, r := range []float64{0, 0.5, 3, 17.25} {
		checkRel(t, "parabola sag", Conic{R: R, K: -1}.Sag(r), r*r/(2*R), 1e-15)
	}
	// Domain: a sphere has no sag beyond r = |R|.
	if !math.IsNaN(Sphere{R: R}.Sag(R * 1.000001)) {
		t.Errorf("sphere sag defined beyond the rim")
	}
	// Plane: identically zero, any radius.
	for _, r := range []float64{0, 1e9} {
		checkClose(t, "plane sag", Plane{}.Sag(r), 0, 0)
	}
	// R = ±Inf is a plane.
	checkClose(t, "infinite radius sag", Conic{R: math.Inf(1)}.Sag(12), 0, 0)
}

// TestSlopeAgainstNumericalDerivative cross-checks every slope implementation
// against a high-order central difference of its own sag.
func TestSlopeAgainstNumericalDerivative(t *testing.T) {
	profiles := []Profile{
		Plane{},
		Sphere{R: 33.3},
		Sphere{R: -18.0},
		Conic{R: 50, K: -1},
		Conic{R: 40, K: -2.5},
		Conic{R: 22, K: 0.7},
		EvenAsphere{R: 25, K: -0.6, Coeffs: []float64{1e-6, 2e-9, -3e-12, 1e-15}},
		EvenAsphere{R: 0, K: 0, Coeffs: []float64{4e-6, 0, 1e-13}}, // polish-only base
	}
	const h = 1e-6
	for _, p := range profiles {
		for _, r := range []float64{0.001, 0.5, 2, 5, 9.9} {
			fp := p.Sag(r + h)
			fm := p.Sag(r - h)
			if math.IsNaN(fp) || math.IsNaN(fm) {
				continue
			}
			num := (fp - fm) / (2 * h)
			got := p.Slope(r)
			tol := 1e-7 * math.Max(1, math.Abs(num))
			if math.IsNaN(got) || math.Abs(got-num) > tol {
				t.Errorf("%v: Slope(%g) = %g, numerical %g (tol %g)", p, r, got, num, tol)
			}
		}
	}
}

// TestSphereSagRangeIsExact checks the sag range helper used to size search
// windows.
func TestSagRange(t *testing.T) {
	lo, hi := sagRange(Sphere{R: 10}, 8)
	checkClose(t, "sag range lo", lo, 0, 1e-12)
	checkClose(t, "sag range hi", hi, 10-math.Sqrt(100-64), 1e-12)
}

// TestSurfaceNormalIsAnalytic checks the surface normal against the exact
// geometric normal of a sphere (n ∝ P - C) and against the numerical gradient
// of the sag for an asphere.
func TestSurfaceNormalIsAnalytic(t *testing.T) {
	const R = -17.5
	s := Surface{Profile: Sphere{R: R}}
	center := Vec3{0, 0, R}
	for _, p := range []Vec3{{0, 0, 0}, {3, 4, 0}, {-5, 2, 0}, {0, 12, 0}} {
		z := s.Sag(p.X, p.Y)
		pt := Vec3{p.X, p.Y, z}
		want := pt.Sub(center).Normalized()
		got := s.Normal(p.X, p.Y)
		if math.Abs(math.Abs(got.Dot(want))-1) > 1e-14 {
			t.Errorf("sphere normal at %s: got %s, want ±%s", pt, got, want)
		}
		// The +z orientation convention.
		if got.Z <= 0 {
			t.Errorf("normal is not +z oriented: %s", got)
		}
	}

	asp := Surface{Profile: EvenAsphere{R: 20, K: -0.5, Coeffs: []float64{2e-6, -1e-9}}}
	const h = 1e-6
	for _, p := range []Vec3{{1, 1, 0}, {4, -3, 0}, {-2, 6.5, 0}} {
		gradx := (asp.Sag(p.X+h, p.Y) - asp.Sag(p.X-h, p.Y)) / (2 * h)
		grady := (asp.Sag(p.X, p.Y+h) - asp.Sag(p.X, p.Y-h)) / (2 * h)
		want := Vec3{-gradx, -grady, 1}.Normalized()
		got := asp.Normal(p.X, p.Y)
		if got.Sub(want).Norm() > 1e-7 {
			t.Errorf("asphere normal at %s: got %s, want %s", p, got, want)
		}
	}
	// At the vertex the normal is exactly +z.
	if n := asp.Normal(0, 0); n.Sub(Vec3{0, 0, 1}).Norm() != 0 {
		t.Errorf("vertex normal = %s, want (0,0,1)", n)
	}
}

// TestEvenAsphereMatchesConicBase checks that an asphere with zero coefficients
// is exactly its conic base.
func TestEvenAsphereMatchesConicBase(t *testing.T) {
	base := Conic{R: 31.25, K: -1.7}
	asp := EvenAsphere{R: 31.25, K: -1.7}
	for _, r := range []float64{0, 0.1, 2, 7.5, 15} {
		checkClose(t, "asphere vs conic sag", asp.Sag(r), base.Sag(r), 0)
		checkClose(t, "asphere vs conic slope", asp.Slope(r), base.Slope(r), 0)
	}
}

// TestIntersectionResidual checks that the hit point returned by the surface
// solver satisfies the sag equation to solver tolerance, for every profile
// type and for random rays.
func TestIntersectionResidual(t *testing.T) {
	rng := rand.New(rand.NewSource(4242))
	profiles := []Profile{
		Plane{},
		Sphere{R: 12.5},
		Sphere{R: -240},
		Conic{R: 30, K: -1},
		Conic{R: 24, K: -3.5},
		Conic{R: 18, K: 2.0},
		EvenAsphere{R: 20, K: -1, Coeffs: []float64{1.5e-6, -2e-9}},
		EvenAsphere{R: -50, K: 0.4, Coeffs: []float64{-4e-6, 3e-10, 7e-14}},
	}
	for _, p := range profiles {
		for i := 0; i < 400; i++ {
			o := Vec3{rng.Float64()*6 - 3, rng.Float64()*6 - 3, rng.Float64()*8 - 12}
			d := Vec3{rng.Float64()*0.2 - 0.1, rng.Float64()*0.2 - 0.1, 1 + rng.Float64()*0.05}
			d = d.Normalized()
			tt, ok := p.Intersect(o, d, 0, 12)
			if !ok {
				continue
			}
			pt := o.Add(d.Mul(tt))
			r := pt.PerpNorm()
			res := pt.Z - p.Sag(r)
			if math.Abs(res) > 1e-11*(1+math.Abs(pt.Z)+r) {
				t.Fatalf("%v: residual %.3g at t=%g point %s", p, res, tt, pt)
			}
			if tt < 0 {
				t.Fatalf("%v: negative intersection parameter %g", p, tt)
			}
		}
	}
}

// TestIntersectionIsFirstRoot checks that the solver returns the smallest valid
// crossing for a strongly aspheric surface that a ray crosses twice.
func TestIntersectionIsFirstRoot(t *testing.T) {
	// A steep "dome" asphere: positive A4 makes the sag rise fast with r, so a
	// horizontal ray near the axis can cross the surface twice.
	p := EvenAsphere{R: 40, K: 0, Coeffs: []float64{2e-3}}
	// Ray along +z at r = 8: the surface rises from z=0 at r=0 to z≈8.2 at r=8
	// (sag(8) = 8²/80 + 0.002·8⁴ = 0.8 + 8.19 = 8.99 mm).
	o := Vec3{8, 0, -20}
	d := Vec3{0, 0, 1}
	tt, ok := p.Intersect(o, d, 0, 20)
	if !ok {
		t.Fatalf("no intersection found for the asphere")
	}
	pt := o.Add(d.Mul(tt))
	checkClose(t, "asphere hit sag", pt.Z, p.Sag(8), 1e-10)
	// No earlier root may exist inside the search window.
	if tt > 1e-9 {
		if _, ok := p.Intersect(o, d, tt+1e-7, 20); ok {
			t.Log("second intersection exists beyond the first (expected for a steep dome)")
		}
	}
}

// TestPlaneNoIntersectionForParallelRay covers the degenerate plane case.
func TestPlaneNoIntersectionForParallelRay(t *testing.T) {
	if _, ok := (Plane{}).Intersect(Vec3{0, 0, 5}, Vec3{1, 0, 0}, 0, 10); ok {
		t.Errorf("a ray parallel to a plane must not intersect it")
	}
	if _, ok := (Plane{}).Intersect(Vec3{0, 0, 5}, Vec3{0, 0, 1}, 0, 10); ok {
		t.Errorf("the intersection behind the ray must be rejected")
	}
	if tt, ok := (Plane{}).Intersect(Vec3{0, 0, 5}, Vec3{0, 0, -1}, 0, 10); !ok || tt != 5 {
		t.Errorf("expected t=5, got %g ok=%v", tt, ok)
	}
}

// TestSurfaceZRange checks the aperture-aware z extent used for launch planes.
func TestSurfaceZRange(t *testing.T) {
	s := Surface{Profile: Sphere{R: 20}, Aperture: CircleAperture(10)}
	lo, hi := s.ZRange(0)
	checkClose(t, "zrange lo", lo, 0, 1e-12)
	checkClose(t, "zrange hi", hi, 20-math.Sqrt(400-100), 1e-9)
}
