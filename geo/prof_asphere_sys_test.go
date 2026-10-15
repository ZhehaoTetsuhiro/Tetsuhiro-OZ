package geo

import (
	"math"
	"testing"
)

// aspSinglet returns a singlet (object at infinity, n = 1.5, 4 mm thick) whose
// front surface is the given profile and whose back surface is a sphere.
// The clear apertures are 8 mm.
func aspSinglet(t *testing.T, front Profile) *System {
	t.Helper()
	return aspEFLSinglet(t, front)
}

// aspFan returns a deterministic fan of collimated rays entering from z = −40
// (object at infinity), spread over the aperture of the singlet.
func aspFan(lambdaNm float64) []Ray {
	var rays []Ray
	for i := -4; i <= 4; i++ {
		for j := -3; j <= 3; j++ {
			x := 0.9 * float64(i)
			y := 0.9 * float64(j)
			rays = append(rays, NewRay(V(x, y, -40), V(0, 0, 1), lambdaNm))
		}
	}
	return rays
}

// TestAsphereSystemTraceMatchesEvenAsphere checks the full sequential trace
// (intersection, normal, Snell refraction, propagation) of a singlet whose front
// surface is the even-only PolyAsphere against the same singlet with the
// existing EvenAsphere: the final ray must agree to 1e-11 in position and
// direction, and the two traces must have the same status.
func TestAsphereSystemTraceMatchesEvenAsphere(t *testing.T) {
	ev := EvenAsphere{R: 25, K: -0.6, Coeffs: []float64{1e-5, -2e-8, 3e-11}}
	pa := PolyAsphere{R: 25, K: -0.6, Terms: []PolyTerm{{4, 1e-5}, {6, -2e-8}, {8, 3e-11}}}
	sa := aspSinglet(t, ev)
	sb := aspSinglet(t, pa)
	ok := 0
	for _, r := range aspFan(550) {
		ta := sa.Trace(r)
		tb := sb.Trace(r)
		if ta.Status != tb.Status {
			t.Fatalf("ray %v: status even %v, poly %v", r.Origin, ta.Status, tb.Status)
		}
		if !ta.OK() {
			continue
		}
		ok++
		if d := ta.Final().Sub(tb.Final()).Norm(); d > 1e-11 {
			t.Errorf("ray from %v: final positions differ by %.3g mm", r.Origin, d)
		}
		if d := ta.Direction().Sub(tb.Direction()).Norm(); d > 1e-11 {
			t.Errorf("ray from %v: final directions differ by %.3g", r.Origin, d)
		}
	}
	if ok < 30 {
		t.Errorf("only %d rays traversed the singlet", ok)
	}
}

// TestAsphereSystemTraceMatchesConicWhenZero checks that a Qcon asphere with all
// coefficients zero traces exactly like the conic base (1e-11 in position and
// direction).
func TestAsphereSystemTraceMatchesConicWhenZero(t *testing.T) {
	const R, K = 25.0, -0.6
	sc := aspSinglet(t, Conic{R: R, K: K})
	sq := aspSinglet(t, QAsphere{Form: QCon, R: R, K: K, NormR: 8, Coeffs: []float64{0, 0, 0}})
	for _, r := range aspFan(550) {
		tc := sc.Trace(r)
		tq := sq.Trace(r)
		if tc.Status != tq.Status {
			t.Fatalf("status conic %v, qcon %v", tc.Status, tq.Status)
		}
		if !tc.OK() {
			continue
		}
		if d := tc.Final().Sub(tq.Final()).Norm(); d > 1e-11 {
			t.Errorf("final positions differ by %.3g mm", d)
		}
		if d := tc.Direction().Sub(tq.Direction()).Norm(); d > 1e-11 {
			t.Errorf("final directions differ by %.3g", d)
		}
	}
}

// TestAsphereTraceStepsLieOnSurface traces rays through an asphere with odd and
// even terms and checks that every recorded hit point lies on the surface to
// 1e-12 mm, i.e. the intersection used by the tracer is the surface itself.
func TestAsphereTraceStepsLieOnSurface(t *testing.T) {
	front := PolyAsphere{R: 30, K: -0.4, Terms: []PolyTerm{{2, 1e-4}, {3, 2e-5}, {5, -3e-7}}}
	sys := aspSinglet(t, front)
	zv := sys.VertexZ(0)
	hits := 0
	for _, r := range aspFan(550) {
		tr := sys.Trace(r)
		if len(tr.Steps) == 0 {
			continue
		}
		st := tr.Steps[0]
		hits++
		local := V(st.Point.X, st.Point.Y, st.Point.Z-zv)
		if res := math.Abs(local.Z - front.Sag(math.Hypot(local.X, local.Y))); res > 1e-12 {
			t.Errorf("hit at (%g,%g): residual %.3g mm > 1e-12", st.Point.X, st.Point.Y, res)
		}
	}
	if hits < 30 {
		t.Errorf("only %d rays reached the asphere", hits)
	}
}

// TestAsphereSurfaceNormalMatchesGradient checks that the analytic surface
// normal used by the tracer (Surface.Normal = unit(−∂z/∂x, −∂z/∂y, 1), with the
// slope supplied by the new profiles) agrees with the central-difference gradient
// of the sag, to 1e-7 componentwise.
func TestAsphereSurfaceNormalMatchesGradient(t *testing.T) {
	profiles := []Profile{
		PolyAsphere{R: 30, K: -0.4, Terms: []PolyTerm{{2, 1e-4}, {3, 2e-5}, {5, -3e-7}}},
		QAsphere{Form: QCon, R: 30, K: -0.4, NormR: 9, Coeffs: []float64{1e-4, -2e-5, 3e-6}},
		QAsphere{Form: QBfs, R: 45, NormR: 9, Coeffs: []float64{2e-4, -3e-5, 1e-5}},
	}
	pts := [][2]float64{{0.6, 0.2}, {-1.4, 2.1}, {3.2, -2.5}, {5.1, 4.4}, {0.0, 0.0}}
	const h = 1e-6
	for _, p := range profiles {
		s := Surface{Profile: p}
		for _, xy := range pts {
			x, y := xy[0], xy[1]
			gx := (s.Sag(x+h, y) - s.Sag(x-h, y)) / (2 * h)
			gy := (s.Sag(x, y+h) - s.Sag(x, y-h)) / (2 * h)
			want, okw := V(-gx, -gy, 1).Unit()
			if !okw {
				t.Fatalf("%v: degenerate gradient at (%g,%g)", p, x, y)
			}
			got := s.Normal(x, y)
			if d := got.Sub(want).Norm(); d > 1e-7 {
				t.Errorf("%v: normal at (%g,%g) = %v, gradient %v (|diff| %.3g)", p, x, y, got, want, d)
			}
		}
	}
}
