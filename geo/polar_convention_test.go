package geo_test

import (
	"math"
	"math/cmplx"
	"math/rand"
	"testing"

	"tetsuhiro-oz/geo"
	"tetsuhiro-oz/optics"
)

// Convention checks against optics/polar.go, the wave kernel's polarization
// code. This file is an external test package so that it may import optics.

// polAzimuthDiffMod returns the difference of two azimuths modulo π, mapped to
// (−π/2, π/2]. Azimuths are defined only modulo π, and the branch can flip when
// S2 is at round-off level.
func polAzimuthDiffMod(a, b float64) float64 {
	return math.Remainder(a-b, math.Pi)
}

// TestPolEllipseMatchesOptics is the Jones-convention comparison: for random and
// special Jones vectors geo.EllipseOf must reproduce optics.StokesAt field by
// field (S0..S3, azimuth, ellipticity, axis ratio, degree and handedness).
func TestPolEllipseMatchesOptics(t *testing.T) {
	rng := rand.New(rand.NewSource(20261015 + 2))
	cases := [][2]complex128{
		{1, 0}, {0, 1}, {1, 1}, {1, -1}, {1, 1i}, {1, -1i},
		{complex(0.6, 0), complex(0, 0.8)}, {0, 0}, {2i, 1},
	}
	for i := 0; i < 5000; i++ {
		cases = append(cases, [2]complex128{
			complex(rng.NormFloat64(), rng.NormFloat64()),
			complex(rng.NormFloat64(), rng.NormFloat64()),
		})
	}
	worst := 0.0
	for i, c := range cases {
		st := optics.StokesAt(c[0], c[1])
		el := geo.EllipseOf(c[0], c[1])
		for _, d := range []float64{
			math.Abs(el.S0 - st.S0), math.Abs(el.S1 - st.S1),
			math.Abs(el.S2 - st.S2), math.Abs(el.S3 - st.S3),
			math.Abs(el.Ellipticity - st.Ellipt), math.Abs(el.AxisRatio - st.AxRatio),
			math.Abs(el.Degree - st.Degree),
		} {
			if d > worst {
				worst = d
			}
		}
		// The azimuth of a circular state (S1 = S2 = 0) is undefined: round-off in
		// S2 decides the branch of atan2, so it is compared only for elliptical states.
		if st.S0 > 0 && math.Hypot(st.S1, st.S2) > 1e-9*st.S0 {
			if d := math.Abs(polAzimuthDiffMod(el.Azimuth, st.Azimuth)); d > worst {
				worst = d
			}
		}
		if el.Handedness != st.Handedness {
			t.Fatalf("case %d (%v, %v): handedness %q, optics says %q", i, c[0], c[1], el.Handedness, st.Handedness)
		}
	}
	if worst > 1e-12 {
		t.Fatalf("EllipseOf differs from optics.StokesAt by %.3g > 1e-12", worst)
	}
	t.Logf("EllipseOf vs optics.StokesAt over %d Jones vectors: worst |Δ| = %.3g", len(cases), worst)
}

// TestPolHandednessTimeDomain checks the sign convention against the physical
// field. With the e^{-iωt} convention the real field is E(θ) = Re(J·e^{-iθ}),
// θ = ωt, and the signed area it sweeps in the (ex, ey) plane over one period is
//
//	∮ (Ex dEy − Ey dEx) = π·S3 .
//
// So S3 > 0 means counterclockwise rotation for a viewer facing the beam (ex to
// the right, ey up, ex × ey = d), the left-handed case of optics/polar.go. The
// integrand contains only the harmonics 0 and ±2, so the N-point trapezoidal sum
// is exact for any N ≥ 3; N = 64 is used, which keeps the round-off of the sum at
// the 1e-15 level.
func TestPolHandednessTimeDomain(t *testing.T) {
	rng := rand.New(rand.NewSource(77))
	const N = 64
	worst := 0.0
	for i := 0; i < 200; i++ {
		jx := complex(rng.NormFloat64(), rng.NormFloat64())
		jy := complex(rng.NormFloat64(), rng.NormFloat64())
		area := 0.0
		for k := 0; k < N; k++ {
			th := 2 * math.Pi * float64(k) / N
			ph := cmplx.Exp(complex(0, -th))
			x := real(jx * ph)
			y := real(jy * ph)
			dx := real(-1i * jx * ph) // dx/dθ
			dy := real(-1i * jy * ph) // dy/dθ
			area += (x*dy - y*dx) * 2 * math.Pi / N
		}
		el := geo.EllipseOf(jx, jy)
		// The round-off of the sum scales with the intensity S0 of the state.
		if d := math.Abs(area-math.Pi*el.S3) / el.S0; d > worst {
			worst = d
		}
	}
	if worst > 1e-12 {
		t.Fatalf("time-domain swept area differs from π·S3 by %.3g·S0 > 1e-12·S0", worst)
	}
	t.Logf("swept area vs π·S3 over 200 random states: worst |Δ|/S0 = %.3g", worst)
}

// TestPolCircularStatesMatchOptics checks the circular states used in the
// regressions: JonesCircular(true) is left-handed in optics.StokesAt.
func TestPolCircularStatesMatchOptics(t *testing.T) {
	l := geo.JonesCircular(true)
	st := optics.StokesAt(l[0], l[1])
	if st.Handedness != "左旋" || math.Abs(st.S3-1) > 1e-12 {
		t.Fatalf("JonesCircular(true): optics handedness %q, S3 = %g", st.Handedness, st.S3)
	}
	r := geo.JonesCircular(false)
	st = optics.StokesAt(r[0], r[1])
	if st.Handedness != "右旋" || math.Abs(st.S3+1) > 1e-12 {
		t.Fatalf("JonesCircular(false): optics handedness %q, S3 = %g", st.Handedness, st.S3)
	}
}
