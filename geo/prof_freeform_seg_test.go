package geo

import (
	"fmt"
	"math"
	"sort"
	"testing"
)

// ---------------------------------------------------------------------------
// A real freeform element: the off-axis segment of a parabola
// ---------------------------------------------------------------------------
//
// An "off-axis parabola" (OAP) is manufactured as a segment cut out of a
// larger parent paraboloid and used in a frame whose z axis is the surface
// normal at the segment centre. In such a frame the surface is neither
// rotationally symmetric nor expressible as a conicoid: its exact sag solves
// the implicit equation
//
//	(d + x c + z s)² + y² = 4f (z c − x s) + d² ,  c = cosθ, s = sinθ,
//	tanθ = d/(2f),
//
// where f is the parent focal length and d the decentre of the segment centre
// from the parent vertex. Its Taylor series about the segment centre is an
// ordinary XY polynomial with i+j ≥ 2, which is exactly what the freeform
// profile is for. This file checks that a truncated series built by the
// package's XYPolynomial reproduces the exact surface (sag and normal) and
// focuses a collimated beam onto the parent focus with a spot frozen at the
// level of the truncation.
//
// The exact reference used here is an independent implementation (the
// quadratic root of the implicit equation, with the normal taken from the
// parent paraboloid and rotated into the local frame); it is not the package's
// series machinery, and its own consistency is checked against the parent
// constraint below.

// ---------------------------------------------------------------------------
// Test-local bivariate power series (canonical order, deterministic sums)
// ---------------------------------------------------------------------------

// ffMon is one monomial c·x^i·y^j of a truncated bivariate series.
type ffMon struct {
	i, j int
	c    float64
}

// ffSerCanon sorts by (i+j, i, j), drops zero coefficients and merges equal
// monomials, so every series built from the same input has the same order of
// summation.
func ffSerCanon(ms []ffMon) []ffMon {
	sort.Slice(ms, func(a, b int) bool {
		if ms[a].i+ms[a].j != ms[b].i+ms[b].j {
			return ms[a].i+ms[a].j < ms[b].i+ms[b].j
		}
		if ms[a].i != ms[b].i {
			return ms[a].i < ms[b].i
		}
		return ms[a].j < ms[b].j
	})
	w := 0
	for _, m := range ms {
		if m.c == 0 {
			continue
		}
		if w > 0 && ms[w-1].i == m.i && ms[w-1].j == m.j {
			ms[w-1].c += m.c
			continue
		}
		ms[w] = m
		w++
	}
	return ms[:w]
}

// ffSerTrunc drops every monomial of total degree above maxDeg.
func ffSerTrunc(a []ffMon, maxDeg int) []ffMon {
	out := make([]ffMon, 0, len(a))
	for _, m := range a {
		if m.i+m.j <= maxDeg {
			out = append(out, m)
		}
	}
	return out
}

// ffSerAdd returns a+b truncated to maxDeg.
func ffSerAdd(a, b []ffMon, maxDeg int) []ffMon {
	return ffSerCanon(ffSerTrunc(append(append([]ffMon(nil), a...), b...), maxDeg))
}

// ffSerScale returns s·a.
func ffSerScale(a []ffMon, s float64) []ffMon {
	out := make([]ffMon, len(a))
	for k, m := range a {
		out[k] = ffMon{m.i, m.j, m.c * s}
	}
	return out
}

// ffSerMul returns a·b truncated to maxDeg.
func ffSerMul(a, b []ffMon, maxDeg int) []ffMon {
	var out []ffMon
	for _, ma := range a {
		for _, mb := range b {
			if ma.i+mb.i+ma.j+mb.j > maxDeg {
				continue
			}
			out = append(out, ffMon{ma.i + mb.i, ma.j + mb.j, ma.c * mb.c})
		}
	}
	return ffSerCanon(out)
}

// ffSerEval evaluates a series at (x, y).
func ffSerEval(a []ffMon, x, y float64) float64 {
	var s float64
	for _, m := range a {
		s += m.c * powInt(x, m.i) * powInt(y, m.j)
	}
	return s
}

// ffSerRecip returns 1/a for a series whose constant term is non-zero, by the
// fixed-point iteration u ← (1/a₀)(1 − (a − a₀)u).
func ffSerRecip(a []ffMon, maxDeg int) []ffMon {
	c0 := 0.0
	var rest []ffMon
	for _, m := range a {
		if m.i == 0 && m.j == 0 {
			c0 = m.c
		} else {
			rest = append(rest, m)
		}
	}
	inv := 1 / c0
	u := []ffMon{{0, 0, inv}}
	one := []ffMon{{0, 0, 1}}
	for k := 0; k < 8; k++ {
		ru := ffSerScale(ffSerMul(rest, u, maxDeg), -1)
		u = ffSerScale(ffSerAdd(one, ru, maxDeg), inv)
	}
	return u
}

// ---------------------------------------------------------------------------
// The exact off-axis parabola segment (test-local reference)
// ---------------------------------------------------------------------------

// ffOAPSegment is the exact off-axis segment of the parent parabola
// z_p = −(x_p²+y_p²)/(4f) in the local frame centred at (d, 0, −d²/(4f)) with
// local z along the surface normal there (a rotation Ry(θ) with
// tanθ = d/(2f) — the same convention as CoordBreak, whose Pose
// {Dx: d, Dz: −d²/(4f), TiltY: θ} places such a frame).
type ffOAPSegment struct {
	f, d, theta float64
}

// cos, sin and the local-to-parent map.
func (s ffOAPSegment) cs() (c, sn float64) { return math.Cos(s.theta), math.Sin(s.theta) }

// toParent maps a local point into the parent frame (p_P = P0 + Ry(θ)·p).
func (s ffOAPSegment) toParent(x, y, z float64) Vec3 {
	c, sn := s.cs()
	return V(s.d+x*c+z*sn, y, -s.d*s.d/(4*s.f)-x*sn+z*c)
}

// toLocal maps a parent-frame direction into the local frame (Qᵀ·v).
func (s ffOAPSegment) toLocal(v Vec3) Vec3 {
	c, sn := s.cs()
	return V(v.X*c-v.Z*sn, v.Y, v.X*sn+v.Z*c)
}

// sag returns the exact local sag: the root of
//
//	s²z² + z(2(d+xc)s + 4fc) + (d+xc)² + y² − d² − 4fxs = 0
//
// that vanishes at the origin, computed in the cancellation-free form.
func (s ffOAPSegment) sag(x, y float64) float64 {
	c, sn := s.cs()
	A2 := sn * sn
	A1 := 2*(s.d+x*c)*sn + 4*s.f*c
	A0 := (s.d+x*c)*(s.d+x*c) + y*y - s.d*s.d - 4*s.f*x*sn
	disc := A1*A1 - 4*A2*A0
	if disc < 0 {
		return math.NaN()
	}
	return -2 * A0 / (A1 + math.Sqrt(disc))
}

// normal returns the exact unit normal in the local frame: the parent normal
// (x_p/(2f), y_p/(2f), 1) rotated into the local frame.
func (s ffOAPSegment) normal(x, y float64) Vec3 {
	z := s.sag(x, y)
	if math.IsNaN(z) {
		return Vec3{}
	}
	pp := s.toParent(x, y, z)
	np, ok := V(pp.X/(2*s.f), pp.Y/(2*s.f), 1).Unit()
	if !ok {
		return Vec3{}
	}
	n, ok := s.toLocal(np).Unit()
	if !ok {
		return Vec3{}
	}
	return n
}

// series returns the Taylor series of the sag to maxDeg, by iterating the
// fixed point z = −(A0 + s²z²)/A1 (each pass doubles the correct order).
func (s ffOAPSegment) series(maxDeg int) []ffMon {
	c, sn := s.cs()
	A2 := []ffMon{{0, 0, sn * sn}}
	A1 := []ffMon{{0, 0, 2*s.d*sn + 4*s.f*c}, {1, 0, 2 * c * sn}}
	A0 := []ffMon{{2, 0, c * c}, {0, 2, 1}}
	R := ffSerRecip(A1, maxDeg)
	z := []ffMon{}
	for k := 0; k < 8; k++ {
		inner := ffSerMul(A2, ffSerMul(z, z, maxDeg), maxDeg)
		z = ffSerScale(ffSerMul(ffSerAdd(A0, inner, maxDeg), R, maxDeg), -1)
	}
	return z
}

// ffOAPProfile adapts the exact segment to the Profile interface so that it can
// be traced. The intersection delegates to the parent conic through the rigid
// frame transform (a rotation preserves the ray parameter of a unit
// direction), so it is exact.
type ffOAPProfile struct {
	seg ffOAPSegment
}

// SagXY implements FreeformProfile.
func (p ffOAPProfile) SagXY(x, y float64) float64 { return p.seg.sag(x, y) }

// NormalXY implements FreeformProfile.
func (p ffOAPProfile) NormalXY(x, y float64) Vec3 { return p.seg.normal(x, y) }

// IntersectXY implements FreeformProfile via the parent parabola.
func (p ffOAPProfile) IntersectXY(o, d Vec3, tMin, rMax float64) (float64, bool) {
	c, sn := p.seg.cs()
	// local → parent: rotation Ry(θ) applied to the ray, origin moved to P0.
	op := V(p.seg.d+o.X*c+o.Z*sn, o.Y, -p.seg.d*p.seg.d/(4*p.seg.f)-o.X*sn+o.Z*c)
	dp := V(d.X*c+d.Z*sn, d.Y, -d.X*sn+d.Z*c)
	return Conic{R: -2 * p.seg.f, K: -1}.Intersect(op, dp, tMin, rMax)
}

// ZRangeXY samples the sag over the aperture disc.
func (p ffOAPProfile) ZRangeXY(rMax float64) (float64, float64) {
	zmin, zmax := math.Inf(1), math.Inf(-1)
	for k := 0; k <= 16; k++ {
		r := rMax * float64(k) / 16
		for j := 0; j < 24; j++ {
			a := 2 * math.Pi * float64(j) / 24
			z := p.seg.sag(r*math.Cos(a), r*math.Sin(a))
			if math.IsNaN(z) {
				continue
			}
			if z < zmin {
				zmin = z
			}
			if z > zmax {
				zmax = z
			}
		}
	}
	return zmin, zmax
}

// Sag implements Profile (the +y meridian).
func (p ffOAPProfile) Sag(r float64) float64 { return p.seg.sag(0, r) }

// Slope implements Profile (dz/dy on the +y meridian, by a central difference
// of the exact sag with a 4th-order stencil — this reference is only used for
// shape and system checks, never inside the analysed solver).
func (p ffOAPProfile) Slope(r float64) float64 {
	const h = 1e-4
	return (-p.seg.sag(0, r+2*h) + 8*p.seg.sag(0, r+h) - 8*p.seg.sag(0, r-h) + p.seg.sag(0, r-2*h)) / (12 * h)
}

// Curvature implements Profile: the mean vertex curvature
// −(1 + cos²θ)/(2 d sinθ + 4 f cosθ).
func (p ffOAPProfile) Curvature() float64 {
	c, sn := p.seg.cs()
	return -(1 + c*c) / (2*p.seg.d*sn + 4*p.seg.f*c)
}

// Intersect implements Profile.
func (p ffOAPProfile) Intersect(o, d Vec3, tMin, rMax float64) (float64, bool) {
	return p.IntersectXY(o, d, tMin, rMax)
}

// RotationallySymmetric implements FreeformProfile: the segment is symmetric
// about its own axis only for θ = 0 (the undecentred parent).
func (p ffOAPProfile) RotationallySymmetric() bool { return p.seg.theta == 0 }

// String implements Profile.
func (p ffOAPProfile) String() string {
	return fmt.Sprintf("exact OAP segment f=%g d=%g", p.seg.f, p.seg.d)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

// ffSegmentCase is the configuration used by the tests of this file: a parent
// paraboloid of 100 mm focal length with the segment decentred by 20 mm
// (θ = atan(0.1) = 5.71°) and an 8 mm beam semi-aperture.
const (
	ffSegF = 100.0
	ffSegD = 20.0
	ffSegA = 8.0
)

// TestOffAxisParabolaSegmentSeries checks the freeform representation of the
// off-axis segment against the exact reference:
//
//   - the reference itself must satisfy the parent paraboloid constraint,
//     z_p + (x_p²+y_p²)/(4f) = 0, to round-off (measured 4e-16 mm);
//   - the truncated series must converge to the exact sag as the degree grows
//     (orders 4/6/8 measured 4.9e-8 / 2.6e-12 / 3.3e-16 mm over the aperture);
//   - the order-8 XYPolynomial built from the series must match the exact
//     sag to ≤ 1e-13 mm and the exact normal to ≤ 1e-9 (the normal is
//     dominated by the derivative of the residual series term).
func TestOffAxisParabolaSegmentSeries(t *testing.T) {
	theta := math.Atan(ffSegD / (2 * ffSegF))
	seg := ffOAPSegment{f: ffSegF, d: ffSegD, theta: theta}
	pts := ffGrid(ffSegA*math.Sqrt2, 8) // 81 points covering the aperture disc

	// (1) The exact reference against the parent constraint.
	worstSelf := 0.0
	for _, xy := range pts {
		z := seg.sag(xy[0], xy[1])
		p := seg.toParent(xy[0], xy[1], z)
		res := math.Abs(p.Z + (p.X*p.X+p.Y*p.Y)/(4*ffSegF))
		if res > worstSelf {
			worstSelf = res
		}
	}
	if worstSelf > 1e-14 {
		t.Errorf("the exact segment sag violates the parent constraint by %.3g mm", worstSelf)
	}

	// (2) Series convergence.
	bounds := map[int]float64{4: 1e-6, 6: 1e-10, 8: 1e-14}
	for _, deg := range []int{4, 6, 8} {
		ms := seg.series(deg)
		worst := 0.0
		for _, xy := range pts {
			d := math.Abs(ffSerEval(ms, xy[0], xy[1]) - seg.sag(xy[0], xy[1]))
			if d > worst {
				worst = d
			}
		}
		t.Logf("  series order %d (%d terms): worst sag error %.3g mm (bound %.0e)", deg, len(ms), worst, bounds[deg])
		if worst > bounds[deg] {
			t.Errorf("order %d: worst sag error %.3g mm exceeds %.0e", deg, worst, bounds[deg])
		}
	}

	// (3) The order-8 XYPolynomial against the exact sag and normal.
	ms := seg.series(8)
	terms := make([]XYTerm, len(ms))
	for k, m := range ms {
		terms[k] = XYTerm{I: m.i, J: m.j, Coeff: m.c}
	}
	ff := ffMustBuild(t, 0, 0, 0, terms)
	if ff.RotationallySymmetric() {
		t.Error("the tilted segment was reported as rotationally symmetric")
	}
	worstSag, worstNormal := 0.0, 0.0
	for _, xy := range pts {
		x, y := xy[0], xy[1]
		sg := ff.SagXY(x, y)
		if d := math.Abs(sg - seg.sag(x, y)); d > worstSag {
			worstSag = d
		}
		na, nb := ff.NormalXY(x, y), seg.normal(x, y)
		if d := na.Sub(nb).Norm(); d > worstNormal {
			worstNormal = d
		}
	}
	t.Logf("  order-8 XYPolynomial: worst sag error %.3g mm (bound 1e-13), "+
		"worst normal error %.3g (bound 1e-9); reference self-check %.3g mm", worstSag, worstNormal, worstSelf)
	if worstSag > 1e-13 {
		t.Errorf("order-8 freeform sag differs from the exact segment by %.3g mm", worstSag)
	}
	if worstNormal > 1e-9 {
		t.Errorf("order-8 freeform normal differs from the exact segment by %.3g", worstNormal)
	}
}

// TestOffAxisParabolaSegmentSpotRMSFrozen traces the collimated beam through
// the freeform off-axis segment (order 8) and through the exact reference,
// both placed with the coordinate break Pose {Dx: d, Dz: −d²/(4f), TiltY: θ}.
// The reflected beam must converge onto the parent focus (0, 0, −f):
//
//   - the two spellings must agree ray by ray (≤ 1e-9 mm; measured 1.8e-14 mm,
//     the reflection of the 4e-15 mm series residual at the focus);
//   - the spot RMS at the focus is frozen for both (the exact one at solver
//     round-off, the freeform one at the series truncation level);
//   - the order-6 series, truncated more coarsely, shows a larger residual —
//     the check that the frozen numbers follow the truncation and not the
//     sampling.
func TestOffAxisParabolaSegmentSpotRMSFrozen(t *testing.T) {
	// Frozen with the Hexapolar(8) sampling at the analytic parent focus
	// (0, 0, −f); both are at solver round-off because the order-8 series
	// already reproduces the exact segment to 4e-15 mm.
	const (
		frozenExactRMS    = 5.16558252666134e-15
		frozenFreeformRMS = 6.14181564025635e-15
	)
	theta := math.Atan(ffSegD / (2 * ffSegF))
	seg := ffOAPSegment{f: ffSegF, d: ffSegD, theta: theta}
	pose := func() *CoordBreak {
		return posePtr(CoordBreak{Dx: ffSegD, Dz: -ffSegD * ffSegD / (4 * ffSegF), TiltY: theta})
	}
	build := func(prof Profile) *System {
		return mustSystem(t, Air{}, []Surface{{
			Profile: prof, Mirror: true,
			Aperture: CircleAperture(ffSegA + 2), Pose: pose(), Name: "OAP segment",
		}}, nil, nil)
	}
	terms := func(deg int) []XYTerm {
		ms := seg.series(deg)
		out := make([]XYTerm, len(ms))
		for k, m := range ms {
			out[k] = XYTerm{I: m.i, J: m.j, Coeff: m.c}
		}
		return out
	}
	exactSys := build(ffOAPProfile{seg: seg})
	ffSys := build(ffMustBuild(t, 0, 0, 0, terms(8)))
	ffSys6 := build(ffMustBuild(t, 0, 0, 0, terms(6)))
	// The beam is centred on the segment, i.e. 20 mm off the global axis; the
	// parent focus is on the global axis at (0, 0, −f).
	pup := Pupil{Z: -40, Radius: ffSegA, Offset: V(ffSegD, 0, 0)}
	focus := V(0, 0, -ffSegF)

	sdE, err := exactSys.SpotDiagram(FieldTan(0, 0), pup, Hexapolar(8), -ffSegF, testLambda)
	if err != nil {
		t.Fatal(err)
	}
	stE := sdE.Stats(ExplicitReference(focus))
	checkClose(t, "exact OAP segment spot RMS", stE.RMS, frozenExactRMS, 1e-16)
	checkLess(t, "exact OAP segment RMS bound", stE.RMS, 1e-12)

	sdF, err := ffSys.SpotDiagram(FieldTan(0, 0), pup, Hexapolar(8), -ffSegF, testLambda)
	if err != nil {
		t.Fatal(err)
	}
	stF := sdF.Stats(ExplicitReference(focus))
	checkClose(t, "freeform OAP segment spot RMS", stF.RMS, frozenFreeformRMS, 1e-15)
	checkLess(t, "freeform OAP segment RMS bound", stF.RMS, 1e-12)

	worstRay := 0.0
	sagErr := 0.0
	for i := range sdE.Points {
		pe, pf := sdE.Points[i], sdF.Points[i]
		if pe.Status != pf.Status {
			t.Fatalf("ray %d: status %v vs %v", i, pe.Status, pf.Status)
		}
		if !pe.Status.OK() {
			continue
		}
		if d := pe.Hit.Sub(pf.Hit).Norm(); d > worstRay {
			worstRay = d
		}
	}
	if worstRay > 1e-9 {
		t.Errorf("the freeform segment focuses %.3g mm away from the exact segment", worstRay)
	}
	// Density independence: the focused spot is at the truncation level for
	// every sampling.
	worstFF := 0.0
	for _, n := range []int{4, 8, 16, 32} {
		sd, err := ffSys.SpotDiagram(FieldTan(0, 0), pup, Hexapolar(n), -ffSegF, testLambda)
		if err != nil {
			t.Fatal(err)
		}
		r := sd.Stats(ExplicitReference(focus)).RMS
		if r > worstFF {
			worstFF = r
		}
		if r > 1e-12 {
			t.Errorf("freeform segment hexapolar(%d): RMS %.3g mm > 1e-12", n, r)
		}
		t.Logf("  freeform segment hexapolar(%2d): RMS=%.6g mm", n, r)
	}
	// The coarser order-6 series has an orders-of-magnitude larger truncation
	// residual: the frozen values track the truncation, not the sampling.
	sd6, err := ffSys6.SpotDiagram(FieldTan(0, 0), pup, Hexapolar(8), -ffSegF, testLambda)
	if err != nil {
		t.Fatal(err)
	}
	st6 := sd6.Stats(ExplicitReference(focus))
	if !(st6.RMS > 10*stF.RMS) {
		t.Errorf("order-6 truncation RMS %.3g is not larger than the order-8 value %.3g", st6.RMS, stF.RMS)
	}
	for _, xy := range ffGrid(ffSegA, 4) {
		if d := math.Abs(ffMustBuild(t, 0, 0, 0, terms(6)).SagXY(xy[0], xy[1]) - seg.sag(xy[0], xy[1])); d > sagErr {
			sagErr = d
		}
	}
	t.Logf("freeform off-axis segment: exact RMS=%.6g, order-8 RMS=%.6g, order-6 RMS=%.6g mm; "+
		"per-ray exact/freeform difference %.3g mm; order-6 sag residual %.3g mm",
		stE.RMS, stF.RMS, st6.RMS, worstRay, sagErr)
}
