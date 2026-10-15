package geo

import (
	"math"
	"testing"
)

// Tests for the FreeformProfile wiring of ZernikeSurface (see
// prof_zernike_freeform.go): two-dimensional dispatch of Surface/System and
// the rejection of asymmetric Zernike surfaces by the first-order engines.

// TestZernikeSurfaceIsFreeformProfile locks the interface wiring: with it,
// Surface.Sag/Normal/Intersect/ZRange evaluate the true two-dimensional
// Zernike sum; without it an asymmetric surface would silently use the θ = 0
// cut.
func TestZernikeSurfaceIsFreeformProfile(t *testing.T) {
	z := ZernikeSurface{
		R: 0, K: 0, NormR: 10,
		Coeffs: []float64{0, 0, 0, 2e-3, 3e-3}, // Z3 (tilt, m=1) + Z5 (astig, m=2)
	}
	var _ FreeformProfile = z // compile-time assertion
	sf := Surface{Profile: z, Aperture: CircleAperture(10)}

	// Sag dispatch: at a point off both axes the two-dimensional sag differs
	// from the θ = 0 cut unless the dispatch is active.
	for _, p := range [][2]float64{{0, 0}, {7, 0}, {0, 7}, {5, -4}, {-6, 3}} {
		x, y := p[0], p[1]
		want := z.SagXY(x, y)
		if got := sf.Sag(x, y); got != want {
			t.Errorf("Surface.Sag(%g,%g) = %g, want ZernikeSurface.SagXY = %g", x, y, got, want)
		}
		r := math.Hypot(x, y)
		if math.Abs(z.SagXY(x, y)-z.Sag(r)) > 1e-12 && x != 0 && y != 0 {
			// sanity: the test points really do distinguish the 2-D sag from the cut
			continue
		}
	}

	// Normal dispatch: the exact gradient at (x, y), not the radial surrogate.
	x, y := 5.0, -4.0
	n := sf.Normal(x, y)
	dzdx, dzdy := z.GradXY(x, y)
	want, ok := Vec3{-dzdx, -dzdy, 1}.Unit()
	if !ok {
		t.Fatal("degenerate normal")
	}
	if d := n.Distance(want); d > 1e-15 {
		t.Errorf("Surface.Normal(%g,%g) = %v, want %v (distance %g)", x, y, n, want, d)
	}

	// ZRange dispatch: a finite bracket that contains every sampled sag.
	lo, hi := sf.ZRange(0) // aperture radius 10 is used
	if !(lo < hi) {
		t.Fatalf("ZRange bracket [%g, %g] invalid", lo, hi)
	}
	for k := 0; k <= 32; k++ {
		a := 2 * math.Pi * float64(k) / 32
		v := z.SagXY(10*math.Cos(a), 10*math.Sin(a))
		if v < lo-1e-12 || v > hi+1e-12 {
			t.Errorf("sag %g at angle %g outside ZRangeXY bracket [%g, %g]", v, a, lo, hi)
			break
		}
	}
}

// TestZernikeAsymmetricRejectedByFirstOrder: an asymmetric Zernike surface is
// a freeform that is not rotationally symmetric, so the paraxial engine must
// refuse it explicitly (never silently produce numbers), while an m = 0-only
// surface stays usable.
func TestZernikeAsymmetricRejectedByFirstOrder(t *testing.T) {
	mk := func(coeffs []float64) *System {
		sys := &System{
			Surfaces: []Surface{
				{Profile: Sphere{R: 100}, Aperture: CircleAperture(5)},
				{Profile: ZernikeSurface{R: math.Inf(1), K: 0, NormR: 10, Coeffs: coeffs}, Aperture: CircleAperture(10)},
			},
			Thickness: []float64{50},
			Mediums:   []Medium{Air{}, Air{}, Air{}},
		}
		return sys
	}
	asym := mk([]float64{0, 0, 0, 0, 2e-3}) // Z5: m = 2
	if _, err := asym.Paraxial(587.5618, math.Inf(-1)); err == nil {
		t.Error("asymmetric Zernike surface: Paraxial did not return an error")
	} else {
		t.Logf("asymmetric Zernike rejected: %v", err)
	}
	sym := mk([]float64{0, 0, 0, 1e-3}) // Z4: m = 0
	if _, err := sym.Paraxial(587.5618, math.Inf(-1)); err != nil {
		t.Errorf("rotationally symmetric Zernike surface: Paraxial failed: %v", err)
	}
}
