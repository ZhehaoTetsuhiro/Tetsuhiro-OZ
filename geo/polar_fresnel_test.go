package geo

import (
	"math"
	"math/cmplx"
	"math/rand"
	"testing"
)

// polCheck fails when |got − want| > tol for complex values.
func polCheck(t *testing.T, name string, got, want complex128, tol float64) {
	t.Helper()
	if d := cmplx.Abs(got - want); !(d <= tol) {
		t.Errorf("%s: got %.17g%+.17gi, want %.17g%+.17gi (|diff| = %.3g > %.3g)",
			name, real(got), imag(got), real(want), imag(want), d, tol)
	}
}

// polCheckReal fails when |got − want| > tol (NaN fails).
func polCheckReal(t *testing.T, name string, got, want, tol float64) {
	t.Helper()
	if d := math.Abs(got - want); !(d <= tol) {
		t.Errorf("%s: got %.17g, want %.17g (|diff| = %.3g > %.3g)", name, got, want, d, tol)
	}
}

// TestPolNormalIncidenceAnalytic is acceptance (a): at normal incidence
// t = 2n1/(n1+n2) and r = (n1−n2)/(n1+n2), identical for s and p, to 1e-12.
func TestPolNormalIncidenceAnalytic(t *testing.T) {
	pairs := [][2]float64{
		{1, 1.5}, {1.5, 1}, {1, 2.4}, {1.33, 1.5168}, {1.7, 1.0},
		{1, 1}, {1.0003, 1.4585}, {1, 3.5}, {2.0, 1.2},
	}
	for _, p := range pairs {
		n1, n2 := p[0], p[1]
		a, cosT, err := FresnelAmplitudes(complex(n1, 0), complex(n2, 0), 1)
		if err != nil {
			t.Fatalf("n1=%g n2=%g: %v", n1, n2, err)
		}
		wantT := complex(2*n1/(n1+n2), 0)
		wantR := complex((n1-n2)/(n1+n2), 0)
		polCheck(t, "t_s", a.TS, wantT, 1e-12)
		polCheck(t, "t_p", a.TP, wantT, 1e-12)
		polCheck(t, "r_s", a.RS, wantR, 1e-12)
		polCheck(t, "r_p", a.RP, wantR, 1e-12)
		polCheck(t, "cosT", cosT, 1, 1e-12)
	}
}

// TestPolBrewsterAngle is acceptance (b): at θB = atan(n2/n1), r_p = 0 to 1e-12.
// The s coefficient there is the closed form |r_s| = |n1²−n2²|/(n1²+n2²), which
// is checked as well so that the test cannot pass with trivial coefficients.
func TestPolBrewsterAngle(t *testing.T) {
	pairs := [][2]float64{{1, 1.5}, {1, 1.33}, {1.5, 1}, {1, 2.4}, {1.7, 1.2}}
	for _, p := range pairs {
		n1, n2 := p[0], p[1]
		thetaB := math.Atan2(n2, n1)
		a, _, err := FresnelAmplitudes(complex(n1, 0), complex(n2, 0), math.Cos(thetaB))
		if err != nil {
			t.Fatalf("n1=%g n2=%g: %v", n1, n2, err)
		}
		polCheck(t, "r_p at Brewster", a.RP, 0, 1e-12)
		wantRS := math.Abs(n1*n1-n2*n2) / (n1*n1 + n2*n2)
		polCheckReal(t, "|r_s| at Brewster", cmplx.Abs(a.RS), wantRS, 1e-12)
	}
}

// TestPolEnergyConservationLossless is acceptance (c): for a lossless interface
// R + (n2·cosθ2)/(n1·cosθ1)·T = 1 to 1e-12 for both polarizations, with T = |t|².
func TestPolEnergyConservationLossless(t *testing.T) {
	rng := rand.New(rand.NewSource(20261015))
	worst := 0.0
	checked := 0
	for i := 0; i < 20000; i++ {
		n1 := 1 + 1.5*rng.Float64()
		n2 := 1 + 2.5*rng.Float64()
		cosI := 0.01 + 0.99*rng.Float64()
		if (n1/n2)*(n1/n2)*(1-cosI*cosI) >= 1 {
			continue // total internal reflection: no transmitted power
		}
		a, cosT, err := FresnelAmplitudes(complex(n1, 0), complex(n2, 0), cosI)
		if err != nil {
			t.Fatalf("case %d: %v", i, err)
		}
		if imag(cosT) != 0 {
			t.Fatalf("case %d: cosT is not real below the critical angle: %v", i, cosT)
		}
		fac := n2 * real(cosT) / (n1 * cosI)
		for _, c := range [][2]complex128{{a.RS, a.TS}, {a.RP, a.TP}} {
			R := real(c[0] * cmplx.Conj(c[0]))
			T := real(c[1] * cmplx.Conj(c[1]))
			if e := math.Abs(R + fac*T - 1); e > worst {
				worst = e
			}
		}
		checked++
	}
	if checked < 10000 {
		t.Fatalf("only %d lossless cases checked", checked)
	}
	if worst > 1e-12 {
		t.Fatalf("worst energy residual %.3g over %d cases exceeds 1e-12", worst, checked)
	}
	t.Logf("worst energy residual %.3g over %d random lossless interfaces", worst, checked)
}

// TestPolFresnelTextbookForms checks the coefficients against the textbook
// trigonometric forms (Hecht, Optics, §4.4; Born & Wolf §1.5), independently of
// the admittance formulation used in the code:
//
//	r_s = −sin(θ1−θ2)/sin(θ1+θ2)     t_s = 2 sinθ2 cosθ1 / sin(θ1+θ2)
//	r_p = −tan(θ1−θ2)/tan(θ1+θ2)    t_p = 2 sinθ2 cosθ1 / (sin(θ1+θ2) cos(θ1−θ2))
//
// where the minus sign of r_p is the unified convention (Hecht's r_p is the
// Born & Wolf one).
func TestPolFresnelTextbookForms(t *testing.T) {
	rng := rand.New(rand.NewSource(7))
	worst := 0.0
	for i := 0; i < 5000; i++ {
		n1 := 1 + 1.2*rng.Float64()
		n2 := 1 + 1.2*rng.Float64()
		th1 := 1.5 * rng.Float64() // 0 … 86°
		if n1*math.Sin(th1)/n2 >= 1 {
			continue
		}
		th2 := math.Asin(n1 * math.Sin(th1) / n2)
		a, _, err := FresnelAmplitudes(complex(n1, 0), complex(n2, 0), math.Cos(th1))
		if err != nil {
			t.Fatalf("case %d: %v", i, err)
		}
		rsT := -math.Sin(th1-th2) / math.Sin(th1+th2)
		tsT := 2 * math.Sin(th2) * math.Cos(th1) / math.Sin(th1+th2)
		rpT := -math.Tan(th1-th2) / math.Tan(th1+th2)
		tpT := 2 * math.Sin(th2) * math.Cos(th1) / (math.Sin(th1+th2) * math.Cos(th1-th2))
		for _, d := range []float64{
			math.Abs(real(a.RS) - rsT), math.Abs(real(a.TS) - tsT),
			math.Abs(real(a.RP) - rpT), math.Abs(real(a.TP) - tpT),
		} {
			if d > worst {
				worst = d
			}
		}
		if imag(a.RS) != 0 || imag(a.TP) != 0 {
			t.Fatalf("case %d: coefficients are not real for real indices", i)
		}
	}
	if worst > 1e-12 {
		t.Fatalf("worst deviation from textbook forms %.3g exceeds 1e-12", worst)
	}
	t.Logf("worst deviation from the textbook trigonometric forms: %.3g", worst)
}

// TestPolTotalInternalReflection checks the evanescent branch: beyond the critical
// angle |r_s| = |r_p| = 1 and cosθ2 is purely imaginary with Im > 0 (the forward
// wave decays), the convention of package coating.
func TestPolTotalInternalReflection(t *testing.T) {
	n1, n2 := 1.5, 1.0
	for _, cosI := range []float64{0.1, 0.3, 0.6, 0.7} {
		sinI := math.Sqrt(1 - cosI*cosI)
		if n1*sinI/n2 <= 1 {
			t.Fatalf("test angle is not beyond the critical angle")
		}
		a, cosT, err := FresnelAmplitudes(complex(n1, 0), complex(n2, 0), cosI)
		if err != nil {
			t.Fatal(err)
		}
		polCheckReal(t, "|r_s| TIR", cmplx.Abs(a.RS), 1, 1e-12)
		polCheckReal(t, "|r_p| TIR", cmplx.Abs(a.RP), 1, 1e-12)
		if real(cosT) != 0 || !(imag(cosT) > 0) {
			t.Errorf("cosθ2 = %v beyond the critical angle, want i·κ with κ > 0", cosT)
		}
	}
}

// TestPolMetalMirrorLimit checks that MetalMirrorAmplitudes is the n2 → ∞ limit
// of the dielectric coefficients. The deviation is O(2/(n2·cosθ)), so with
// n2 = 1e15 and cosθ ≥ 0.2 it is below 1e-13 and the 1e-12 bound is not tight.
func TestPolMetalMirrorLimit(t *testing.T) {
	for _, cosI := range []float64{1, 0.8, 0.5, 0.2} {
		a, _, err := FresnelAmplitudes(complex(1, 0), complex(1e15, 0), cosI)
		if err != nil {
			t.Fatal(err)
		}
		m := MetalMirrorAmplitudes()
		polCheck(t, "r_s → −1", a.RS, m.RS, 1e-12)
		polCheck(t, "r_p → −1", a.RP, m.RP, 1e-12)
	}
}

// TestPolFresnelErrors covers the rejected inputs.
func TestPolFresnelErrors(t *testing.T) {
	if _, _, err := FresnelAmplitudes(1, 1.5, 0); err == nil {
		t.Error("grazing incidence (cosI = 0) must be rejected")
	}
	if _, _, err := FresnelAmplitudes(1, 1.5, 1.2); err == nil {
		t.Error("cosI > 1 must be rejected")
	}
	if _, _, err := FresnelAmplitudes(0, 1.5, 0.5); err == nil {
		t.Error("zero index must be rejected")
	}
}
