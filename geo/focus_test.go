package geo

import (
	"math"
	"testing"
)

// axisCrossingZ returns the z at which the traced ray (assumed to lie in the
// y-z plane) crosses the optical axis.
func axisCrossingZ(res TraceResult) float64 {
	p := res.Ray.Origin
	d := res.Ray.Dir
	if d.Y == 0 {
		return math.Inf(signOf(-p.Y))
	}
	t := -p.Y / d.Y
	return p.Z + t*d.Z
}

// distancePointToRay returns the perpendicular distance from the point q to the
// ray line through the trace result.
func distancePointToRay(res TraceResult, q Vec3) float64 {
	v := q.Sub(res.Ray.Origin)
	d := res.Ray.Dir
	return v.Sub(d.Mul(v.Dot(d))).Norm()
}

// TestSingleSphericalSurfaceFocus checks the paraxial focus of a single
// refracting spherical surface, s' = n' R / (n' - n), against real ray tracing
// at a vanishing aperture.
func TestSingleSphericalSurfaceFocus(t *testing.T) {
	for _, tc := range []struct {
		R, n1, n2, h float64
	}{
		{20, 1, 1.5168, 1e-4},
		{20, 1, 1.5168, 1e-6},
		{-30, 1, 1.62, 1e-4}, // concave, seen from the object side
		{25, 1.5, 1.0, 1e-4}, // exit into air from an immersion medium
		{100, 1.33, 1.5168, 1e-4},
	} {
		sys := mustSystem(t, ConstantIndex(tc.n1),
			[]Surface{{Profile: Sphere{R: tc.R}, Aperture: CircleAperture(30)}},
			nil, []Medium{ConstantIndex(tc.n1), ConstantIndex(tc.n2)})
		res := sys.Trace(NewRay(Vec3{0, tc.h, -50}, Vec3{0, 0, 1}, Wavelengthd))
		if !res.OK() {
			t.Fatalf("trace failed: %v", res.Error())
		}
		want := tc.n2 * tc.R / (tc.n2 - tc.n1)
		checkRel(t, "spherical surface focus", axisCrossingZ(res), want, 1e-6)
		// The paraxial engine must agree exactly.
		p, err := sys.Paraxial(Wavelengthd, math.Inf(-1))
		if err != nil {
			t.Fatalf("Paraxial: %v", err)
		}
		checkRel(t, "paraxial back focal distance", p.BFL, want, 1e-12)
		checkRel(t, "paraxial back focal plane", p.BackFocalZ, want, 1e-12)
	}
}

// TestPlanoConvexFocusConvergence verifies that the real marginal-ray focus of
// a plano-convex lens converges to the paraxial focus R/(n-1), and that the
// error falls by a factor four when the aperture is halved (longitudinal
// spherical aberration is quadratic in the aperture).
func TestPlanoConvexFocusConvergence(t *testing.T) {
	n := N_BK7.N(Wavelengthd)
	const f = 100.0
	R := f * (n - 1) // thin-lens maker's formula: 1/f = (n-1)(1/R1 - 1/R2), R2 = Inf

	surfaces := []Surface{
		{Profile: Sphere{R: R}, Aperture: CircleAperture(20), Name: "convex"},
		{Profile: Plane{}, Aperture: CircleAperture(20), Name: "flat"},
	}
	sys := mustSystem(t, Air{}, surfaces, []float64{4}, []Medium{Air{}, N_BK7, Air{}})
	p, err := sys.Paraxial(Wavelengthd, math.Inf(-1))
	if err != nil {
		t.Fatalf("Paraxial: %v", err)
	}
	checkRel(t, "plano-convex EFL = R/(n-1)", p.EFL, R/(n-1), 1e-12)
	checkClose(t, "plano-convex EFL = f", p.EFL, f, 1e-9)

	focusAt := func(h float64) float64 {
		res := sys.Trace(NewRay(Vec3{0, h, -100}, Vec3{0, 0, 1}, Wavelengthd))
		if !res.OK() {
			t.Fatalf("h=%g: trace failed: %v", h, res.Error())
		}
		return axisCrossingZ(res)
	}
	h1, h2, h3 := 4.0, 2.0, 1.0
	e1 := focusAt(h1) - p.BackFocalZ
	e2 := focusAt(h2) - p.BackFocalZ
	e3 := focusAt(h3) - p.BackFocalZ
	// Undercorrected positive singlet: the marginal focus is short of the
	// paraxial focus.
	if !(e1 < 0) {
		t.Fatalf("marginal focus %g is not short of the paraxial focus %g", focusAt(h1), p.BackFocalZ)
	}
	checkClose(t, "halving the aperture quarters the error", e2/e1, 0.25, 0.01)
	checkClose(t, "halving again quarters the error", e3/e2, 0.25, 0.01)
	t.Logf("longitudinal SA: h=%g: %.6g mm, h=%g: %.6g mm, h=%g: %.6g mm", h1, e1, h2, e2, h3, e3)
}

// TestParabolicMirrorPerfectFocus verifies the textbook property of a
// parabolic mirror (K = -1): parallel rays at any height are reflected exactly
// through the focus at R/2, with no spherical aberration at all.
func TestParabolicMirrorPerfectFocus(t *testing.T) {
	const R = 200.0
	focus := Vec3{0, 0, R / 2}
	sys := mustSystem(t, Air{},
		[]Surface{{Profile: Conic{R: R, K: -1}, Mirror: true, Aperture: CircleAperture(60)}},
		nil, nil)
	worst := 0.0
	for _, h := range []float64{1, 5, 10, 25, 40, 55} {
		res := sys.Trace(NewRay(Vec3{0, h, -50}, Vec3{0, 0, 1}, Wavelengthd))
		if !res.OK() {
			t.Fatalf("h=%g: trace failed: %v", h, res.Error())
		}
		d := distancePointToRay(res, focus)
		if d > worst {
			worst = d
		}
		// The reflected ray must travel back towards the object side.
		if res.Ray.Dir.Z >= 0 {
			t.Errorf("h=%g: reflected ray does not travel backwards (d=%s)", h, res.Ray.Dir)
		}
		// The hit point must satisfy the parabola z = r²/(2R) exactly.
		checkClose(t, "parabola hit sag", res.Steps[0].Point.Z, h*h/(2*R), 1e-11)
	}
	if worst > 1e-9 {
		t.Fatalf("parabolic mirror focus error %.3g exceeds 1e-9 mm", worst)
	}
	t.Logf("worst parabolic mirror focus error over 6 heights: %.3g mm", worst)
}

// TestSphericalMirrorHasSphericalAberration is the companion of
// TestParabolicMirrorPerfectFocus: the sphere with the same vertex radius must
// show clear spherical aberration, with the marginal focus short of R/2.
func TestSphericalMirrorHasSphericalAberration(t *testing.T) {
	const R = 200.0
	sys := mustSystem(t, Air{},
		[]Surface{{Profile: Sphere{R: R}, Mirror: true, Aperture: CircleAperture(60)}},
		nil, nil)
	focusAt := func(h float64) float64 {
		res := sys.Trace(NewRay(Vec3{0, h, -50}, Vec3{0, 0, 1}, Wavelengthd))
		if !res.OK() {
			t.Fatalf("h=%g: trace failed: %v", h, res.Error())
		}
		return axisCrossingZ(res)
	}
	// Paraxial limit: the sphere must focus at R/2.
	checkClose(t, "spherical mirror paraxial focus", focusAt(1e-4), R/2, 1e-6)
	// With aperture: the rays cross the axis clearly in front of R/2 and the
	// spread grows with the height (undercorrected mirror).
	f10, f20, f40 := focusAt(10), focusAt(20), focusAt(40)
	if !(f40 < f20 && f20 < f10 && f10 < R/2) {
		t.Fatalf("expected monotone shortening of the focus with height, got %g, %g, %g (R/2=%g)", f10, f20, f40, R/2)
	}
	// The longitudinal aberration is quadratic in the aperture.
	l10 := R/2 - f10
	l20 := R/2 - f20
	checkClose(t, "mirror SA ∝ h²", l20/l10, 4, 0.2)
	t.Logf("spherical mirror axis crossings: h=10: %.6g, h=20: %.6g, h=40: %.6g (paraxial %.1f)", f10, f20, f40, R/2)
}

// TestMirrorDoubleReflectionIsIdentity checks that reflecting the same ray
// twice at the same surface point restores the original direction (the
// reflection map is an involution) — the ray must return along its own path.
func TestMirrorDoubleReflectionIsIdentity(t *testing.T) {
	sys := mustSystem(t, Air{},
		[]Surface{{Profile: Conic{R: -500, K: -0.8}, Mirror: true, Aperture: CircleAperture(50)}},
		nil, nil)
	for _, h := range []float64{0, 3, 12} {
		res := sys.Trace(NewRay(Vec3{0, h, -80}, Vec3{0, 0, 1}, Wavelengthd))
		if !res.OK() {
			t.Fatalf("h=%g: trace failed: %v", h, res.Error())
		}
		back := sys.Trace(res.Ray)
		if !back.OK() {
			t.Fatalf("h=%g: reverse trace failed: %v", h, back.Error())
		}
		if d := back.Ray.Dir.Sub(Vec3{0, 0, 1}).Norm(); d > 1e-15 {
			t.Errorf("h=%g: double reflection direction error %.3g", h, d)
		}
		checkClose(t, "double reflection returns to the hit point", back.Ray.Origin.Sub(res.Ray.Origin).Norm(), 0, 1e-12)
	}
}

// TestTracingThroughReversedSystem checks the Reversed helper: tracing a ray
// through the reversed system with the reversed outgoing ray must reproduce the
// incoming ray exactly.
func TestTracingThroughReversedSystem(t *testing.T) {
	sys := mustSystem(t, Air{},
		[]Surface{
			{Profile: Sphere{R: 40}, Aperture: CircleAperture(20), Name: "L1 front"},
			{Profile: Conic{R: -35, K: -0.5}, Aperture: CircleAperture(20), Name: "L1 rear"},
			{Profile: Plane{}, Aperture: CircleAperture(20), Name: "stop"},
			{Profile: Sphere{R: 60}, Aperture: CircleAperture(20), Name: "L2 front"},
		},
		[]float64{6, 3, 25},
		[]Medium{Air{}, N_BK7, Air{}, N_SK16, Air{}})
	sys.StopIndex = 2
	rev := sys.Reversed()
	start := Vec3{1.5, 3.5, -30}
	dir := Vec3{0.01, -0.02, 1}.Normalized()
	res := sys.Trace(NewRay(start, dir, Wavelengthd))
	if !res.OK() {
		t.Fatalf("forward trace failed: %v", res.Error())
	}
	back := rev.Trace(NewRay(res.Ray.Origin, res.Ray.Dir.Neg(), Wavelengthd))
	if !back.OK() {
		t.Fatalf("reverse trace failed: %v", back.Error())
	}
	// The reversed ray retraces the forward path: at the far end it must be
	// travelling along -dir (the reverse of the original input direction) and
	// its position must lie exactly on the original input line.
	if d := back.Ray.Dir.Add(dir).Norm(); d > 1e-12 {
		t.Errorf("reversed system direction error %.3g", d)
	}
	per := start.Sub(back.Ray.Origin)
	per = per.Sub(dir.Mul(per.Dot(dir)))
	if per.Norm() > 1e-9 {
		t.Errorf("reversed ray is off the original input line by %.3g", per.Norm())
	}
	// Optical and geometrical path lengths match up to the free air gap
	// between the forward start plane and surface 0.
	gap := start.Distance(back.Ray.Origin)
	checkRel(t, "reversed OPL + gap", back.OPL()+gap, res.OPL(), 1e-12)
	checkRel(t, "reversed path + gap", back.Ray.Path+gap, res.Ray.Path, 1e-12)
}

// TestHelmholtzLagrangeInvariant checks that H = n(ȳu - yū) computed from the
// paraxial marginal and chief rays is constant across every surface of a
// multi-surface system (tolerance 1e-10 relative).
func TestHelmholtzLagrangeInvariant(t *testing.T) {
	n1 := N_BK7.N(Wavelengthd)
	n2 := N_SF5.N(Wavelengthd)
	sys := mustSystem(t, Air{},
		[]Surface{
			{Profile: Sphere{R: 61.0}, Aperture: CircleAperture(12)},
			{Profile: Sphere{R: -44.0}, Aperture: CircleAperture(12)},
			{Profile: Sphere{R: -130.0}, Aperture: CircleAperture(12)},
			{Profile: Sphere{R: 90.0}, Aperture: CircleAperture(25)},
		},
		[]float64{6, 2, 40},
		[]Medium{Air{}, N_BK7, N_SF5, Air{}, Air{}})
	sys.StopIndex = 1
	_ = n1
	_ = n2

	p, err := sys.Paraxial(Wavelengthd, math.Inf(-1))
	if err != nil {
		t.Fatalf("Paraxial: %v", err)
	}
	const (
		radius = 6.0
		tanY   = 0.087488663525924 // tan(5°)
	)
	h0 := p.LagrangeAt(0, radius, tanY)
	if h0 == 0 {
		t.Fatalf("degenerate invariant")
	}
	// The invariant is preserved by the refraction and the transfer of the
	// paraxial trace: first verify only at the surface planes.
	for i := 0; i < sys.NumSurfaces(); i++ {
		hi := p.LagrangeAt(i, radius, tanY)
		checkRel(t, "H at surface", hi, h0, 1e-10)
	}
	// Also verify the general evaluation helper.
	checkRel(t, "H via Lagrange()", p.Lagrange(radius, tanY), h0, 1e-12)
}

// TestSkewInvariantRealRays checks the exact rotational invariant of real
// (non-paraxial) rays: L = n(x·d_y - y·d_x) is conserved through refraction and
// reflection at any surface of revolution, because the surface normal satisfies
// x·N_y - y·N_x = 0.
func TestSkewInvariantRealRays(t *testing.T) {
	sys := mustSystem(t, Air{},
		[]Surface{
			{Profile: Sphere{R: 40}, Aperture: CircleAperture(20)},
			{Profile: EvenAsphere{R: -35, K: -0.8, Coeffs: []float64{1e-6}}, Aperture: CircleAperture(20)},
			{Profile: Plane{}, Aperture: CircleAperture(20)},
		},
		[]float64{6, 3},
		[]Medium{Air{}, N_BK7, N_SF5, Air{}})
	worst := 0.0
	for _, start := range []Vec3{
		{3.1, 7.7, -30}, {-4.4, 2.2, -30}, {0.5, -9.3, -30}, {11.0, 0.0, -30},
	} {
		dir := Vec3{0.02, -0.11, 1}.Normalized()
		res := sys.Trace(NewRay(start, dir, Wavelengthd))
		if !res.OK() {
			continue
		}
		// Value before the first surface.
		nBefore := SIndex(sys, 0, Wavelengthd)
		l0 := nBefore * (start.X*dir.Y - start.Y*dir.X)
		for _, st := range res.Steps {
			l := st.NAfter * (st.Point.X*st.DirOut.Y - st.Point.Y*st.DirOut.X)
			if d := math.Abs(l - l0); d > worst {
				worst = d
			}
		}
	}
	if worst > 1e-13 {
		t.Fatalf("skew invariant drift %.3g exceeds 1e-13", worst)
	}
	t.Logf("worst skew invariant drift: %.3g", worst)
}

// SIndex is a tiny helper for tests: the index in front of surface i.
func SIndex(s *System, i int, lambda float64) float64 { return s.IndexBefore(i, lambda) }
