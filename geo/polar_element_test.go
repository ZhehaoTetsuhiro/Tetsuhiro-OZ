package geo

import (
	"math"
	"math/cmplx"
	"testing"
)

// polUnitaryResidual returns max |(J†J − I)_ik| for a 2×2 Jones matrix.
func polUnitaryResidual(j JonesMat) float64 {
	worst := 0.0
	for i := 0; i < 2; i++ {
		for k := 0; k < 2; k++ {
			var s complex128
			for m := 0; m < 2; m++ {
				s += cmplx.Conj(j[m][i]) * j[m][k]
			}
			if i == k {
				s -= 1
			}
			if d := cmplx.Abs(s); d > worst {
				worst = d
			}
		}
	}
	return worst
}

// TestPolElementsUnitaryAndProjector checks the wave plates (unitary) and the
// polarizer (orthogonal projector) on a grid of parameters.
func TestPolElementsUnitaryAndProjector(t *testing.T) {
	worstU, worstP := 0.0, 0.0
	for _, delta := range []float64{0.3, math.Pi / 2, math.Pi, 2.7} {
		for _, phi := range []float64{0, 0.2, 0.9, 1.6} {
			if d := polUnitaryResidual(RetarderJones(delta, phi)); d > worstU {
				worstU = d
			}
		}
	}
	for _, phi := range []float64{0, 0.4, 1.1} {
		// An ideal polarizer is an orthogonal projector: P² = P and P† = P.
		p := PolarizerJones(phi)
		p2 := p.Mul(p)
		for i := 0; i < 2; i++ {
			for k := 0; k < 2; k++ {
				if d := cmplx.Abs(p2[i][k] - p[i][k]); d > worstP {
					worstP = d
				}
				if d := cmplx.Abs(p[i][k] - cmplx.Conj(p[k][i])); d > worstP {
					worstP = d
				}
			}
		}
	}
	if worstU > 1e-12 {
		t.Fatalf("retarder not unitary: residual %.3g", worstU)
	}
	if worstP > 1e-12 {
		t.Fatalf("polarizer not idempotent: residual %.3g", worstP)
	}
}

// polOnePlaneTrace returns a system with one plane interface (air to air, so the
// only effect is the element) and a ray along +z, for element tests.
func polOnePlaneTrace(t *testing.T) (*System, Ray) {
	t.Helper()
	sys, err := NewSystem(Air{}, []Surface{NewSurface(Plane{}, 0)}, nil, nil)
	if err != nil {
		t.Fatal(err)
	}
	return sys, NewRay(V(0, 0, -3), AxisZ, 550)
}

// TestPolQWPCircularRegression is the G5 circular-polarization regression: a
// 45° linear state through a quarter-wave plate with the fast axis along ex is
// circular to 1e-12 (|S3| = S0 with S1 = S2 = 0), left-handed in the optics
// naming for +45° and right-handed for −45°.
func TestPolQWPCircularRegression(t *testing.T) {
	sys, ray := polOnePlaneTrace(t)
	for _, c := range []struct {
		theta float64
		hand  string
		sign  float64
	}{
		{math.Pi / 4, "左旋", 1},
		{-math.Pi / 4, "右旋", -1},
	} {
		opt := &PolOptions{Elements: []PolElementAt{{AfterSurface: 0, J: QuarterWaveJones(0)}}}
		res, err := sys.TracePolarized(ray, JonesLinear(c.theta), opt)
		if err != nil {
			t.Fatal(err)
		}
		el := res.Ellipse
		polCheckReal(t, "S3/S0", el.S3/el.S0, c.sign, 1e-12)
		polCheckReal(t, "S1/S0", el.S1/el.S0, 0, 1e-12)
		polCheckReal(t, "S2/S0", el.S2/el.S0, 0, 1e-12)
		polCheckReal(t, "degree", el.Degree, 1, 1e-12)
		if el.Handedness != c.hand {
			t.Errorf("handedness %q, want %q", el.Handedness, c.hand)
		}
	}
	// A QWP with the fast axis at 45° leaves a linear state along ±45° unchanged up to phase.
	opt := &PolOptions{Elements: []PolElementAt{{AfterSurface: 0, J: QuarterWaveJones(math.Pi / 4)}}}
	res, err := sys.TracePolarized(ray, JonesLinear(math.Pi/4), opt)
	if err != nil {
		t.Fatal(err)
	}
	polCheckReal(t, "QWP(45°) on 45° linear: linear", res.Ellipse.AxisRatio, 0, 1e-12)
}

// TestPolHWPRotatesLinearPolarization: a half-wave plate with fast axis φ turns
// linear polarization at α into linear polarization at 2φ − α (mod π), to 1e-12.
func TestPolHWPRotatesLinearPolarization(t *testing.T) {
	sys, ray := polOnePlaneTrace(t)
	for _, c := range [][2]float64{{0.3, 0.5}, {-0.7, 1.2}, {1.1, -0.35}, {0, 0.25}} {
		alpha, phi := c[0], c[1]
		opt := &PolOptions{Elements: []PolElementAt{{AfterSurface: 0, J: HalfWaveJones(phi)}}}
		res, err := sys.TracePolarized(ray, JonesLinear(alpha), opt)
		if err != nil {
			t.Fatal(err)
		}
		want := 2*phi - alpha
		d := math.Remainder(res.Ellipse.Azimuth-want, math.Pi)
		polCheckReal(t, "HWP azimuth (mod π)", d, 0, 1e-12)
		polCheckReal(t, "HWP keeps linear state", res.Ellipse.AxisRatio, 0, 1e-12)
		polCheckReal(t, "HWP keeps intensity", res.Transmittance, 1, 1e-12)
	}
}

// TestPolCrossedPolarizersExtinction is the extinction regression: a polarizer at
// 0° followed by one at 90° transmits nothing (to 1e-12 in power), while the
// single 0° polarizer transmits the whole aligned state.
func TestPolCrossedPolarizersExtinction(t *testing.T) {
	sys, ray := polOnePlaneTrace(t)
	crossed := &PolOptions{Elements: []PolElementAt{
		{AfterSurface: 0, J: PolarizerJones(0)},
		{AfterSurface: 0, J: PolarizerJones(math.Pi / 2)},
	}}
	res, err := sys.TracePolarized(ray, JonesLinear(0.4), crossed)
	if err != nil {
		t.Fatal(err)
	}
	polCheckReal(t, "crossed-polarizer transmittance", res.Transmittance, 0, 1e-12)
	polCheckReal(t, "crossed-polarizer output amplitude", cmplxAbs2(res.Jones), 0, 1e-12)

	aligned := &PolOptions{Elements: []PolElementAt{{AfterSurface: 0, J: PolarizerJones(0)}}}
	res, err = sys.TracePolarized(ray, JonesLinear(0.0), aligned)
	if err != nil {
		t.Fatal(err)
	}
	polCheckReal(t, "aligned-polarizer transmittance", res.Transmittance, 1, 1e-12)
}

// cmplxAbs2 returns |v|² of a two-component complex vector.
func cmplxAbs2(v [2]complex128) float64 {
	return real(v[0])*real(v[0]) + imag(v[0])*imag(v[0]) + real(v[1])*real(v[1]) + imag(v[1])*imag(v[1])
}
