package geo

import (
	"math"
	"testing"
)

// ---------------------------------------------------------------------------
// Real-ray replay
// ---------------------------------------------------------------------------

// TestGhostReplayMatchesForwardTrace checks that the replay of a ghost path
// reproduces the forward trace exactly on the part the two share: the
// refractions in front of the first reflection site. The replayer walks the
// same Surface.Intersect / Surface.Normal / Refract machinery (through the
// public Surface API and the frame chain of the system), so the hit points and
// directions must agree to round-off — the test uses 1e-15 mm, the values are
// bit-identical in practice.
func TestGhostReplayMatchesForwardTrace(t *testing.T) {
	sys := ghostDoublet(t)
	r := LaunchRay(FieldTan(0.01, 0.005), Pupil{Z: -5, Radius: 0.2}, PupilPointAt(0, 0), Wavelengthd)
	res := sys.Trace(r)
	if !res.OK() {
		t.Fatalf("forward trace failed: %v", res.Status)
	}
	path, err := sys.GhostPathOf(GhostPair{I: 2, J: 1})
	if err != nil {
		t.Fatal(err)
	}
	hit := sys.ReplayGhost(path, r, 190.0)
	if len(hit.Steps) != len(path.Nodes) {
		t.Fatalf("replay stopped after %d of %d nodes (status %v)", len(hit.Steps), len(path.Nodes), hit.Status)
	}
	for i := 0; i < path.Pair.I; i++ {
		a, b := res.Steps[i], hit.Steps[i]
		if d := a.Point.Sub(b.Point).Norm(); d > 1e-15 {
			t.Errorf("step %d: hit point differs by %.3g mm (%s vs %s)", i, d, a.Point, b.Point)
		}
		if d := a.DirOut.Sub(b.DirOut).Norm(); d > 1e-15 {
			t.Errorf("step %d: direction differs by %.3g (%s vs %s)", i, d, a.DirOut, b.DirOut)
		}
	}
	// The first reflection must be the exact mirror law at the recorded
	// surface normal, and it must send the ray backwards.
	first := hit.Steps[path.Pair.I]
	want := Reflect(first.DirIn, first.Normal)
	if d := first.DirOut.Sub(want).Norm(); d > 1e-15 {
		t.Errorf("reflection at surface %d: d' differs from d − 2(d·n)n by %.3g", first.Surface, d)
	}
	if !(first.DirOut.Z < 0) {
		t.Errorf("the first reflection must send the ray towards −z, got %s", first.DirOut)
	}
	// The second reflection returns the ray to +z, happens in the medium the
	// ray arrives from (the flint, gap 2 — not the incident medium of the
	// forward pass), and is followed by the second crossing of surface 2.
	var second, last GhostStep
	nRefl := 0
	for _, st := range hit.Steps {
		if st.Action == GhostReflect {
			nRefl++
			if nRefl == 2 {
				second = st
			}
		}
		last = st
	}
	if nRefl != 2 {
		t.Fatalf("got %d reflections, want 2", nRefl)
	}
	if second.Surface != path.Pair.J || !(second.DirOut.Z > 0) || !second.Backward {
		t.Errorf("second reflection = %s, want the backwards hit of surface %d leaving towards +z",
			second, path.Pair.J)
	}
	if second.NBefore != sys.Mediums[2].N(Wavelengthd) || second.NAfter != second.NBefore {
		t.Errorf("second reflection indices: %g → %g, want the flint %g on both sides",
			second.NBefore, second.NAfter, sys.Mediums[2].N(Wavelengthd))
	}
	if last.Surface != 2 || last.Action != GhostRefract {
		t.Errorf("last node = %s, want a refraction at surface 2", last)
	}
	if last.NBefore != sys.Mediums[2].N(Wavelengthd) || last.NAfter != sys.Mediums[3].N(Wavelengthd) {
		t.Errorf("the last refraction must cross from the flint %g into the air %g, got %g → %g",
			sys.Mediums[2].N(Wavelengthd), sys.Mediums[3].N(Wavelengthd), last.NBefore, last.NAfter)
	}
}

// TestGhostReplayStatuses checks the failure statuses of the replay: a ray
// vignetted at a node stops there with the forward trace's StatusVignetted,
// and a malformed node list (a node whose gap does not follow the medium
// bookkeeping) is reported as degenerate instead of being traced.
func TestGhostReplayStatuses(t *testing.T) {
	// A plate whose rear surface has a small aperture: the ghost reflection at
	// the rear happens far off axis, so it vignettes there.
	sys := mustSystem(t, Air{},
		[]Surface{
			{Profile: Plane{}, Aperture: CircleAperture(50)},
			{Profile: Plane{}, Aperture: CircleAperture(1)},
		},
		[]float64{5}, []Medium{Air{}, ConstantIndex(1.5), Air{}})
	path, err := sys.GhostPathOf(GhostPair{I: 1, J: 0})
	if err != nil {
		t.Fatal(err)
	}
	r := NewRay(V(0, 2, -1), V(0, 0, 1), Wavelengthd)
	hit := sys.ReplayGhost(path, r, 20)
	checkStatus(t, "vignetted ghost replay", hit.Status, StatusVignetted)
	if hit.FailedAt != 1 {
		t.Errorf("vignetting must happen at node 1 (the rear reflection), got %d", hit.FailedAt)
	}
	if len(hit.Steps) != 2 {
		t.Errorf("got %d steps, want 2 (the front refraction and the failing rear hit)", len(hit.Steps))
	}

	// A hand-built node list that does not follow the medium bookkeeping is a
	// programming error, not a physical outcome.
	bad := GhostPath{Pair: GhostPair{I: 1, J: 0}, Nodes: []GhostNode{
		{Surface: 0, Action: GhostRefract, Medium: 1}, // wrong: the ray starts in gap 0
	}}
	hit = sys.ReplayGhost(bad, r, 20)
	checkStatus(t, "malformed ghost path", hit.Status, StatusDegenerate)

	// A degenerate ray is rejected before any interaction.
	degenerate := NewRay(V(0, 0, -1), Zero, Wavelengthd)
	hit = sys.ReplayGhost(path, degenerate, 20)
	checkStatus(t, "degenerate launch ray", hit.Status, StatusDegenerate)
	if len(hit.Steps) != 0 {
		t.Errorf("a degenerate ray must not reach a surface, got %d steps", len(hit.Steps))
	}
}

// ---------------------------------------------------------------------------
// Paraxial ghost trace
// ---------------------------------------------------------------------------

// ghostParaxialRecurse is the independent transcription of the paraxial ghost
// path used as the oracle of TestGhostParaxialRepeatability: it propagates the
// paraxial ray (y, u) through the ghost path p with the y–u recurrence
//
//	transfer    y' = y + (z' − z) u
//	refraction  u' = η u + (η − 1) c y ,  η = n_in/n_out
//	reflection  u' = −u − 2 c y
//
// and returns the height and the slope (dy/dz) at the plane z = imageZ. The
// zs slice holds the vertex positions of the surfaces; an interval crossed
// backwards gives a negative transfer, which is exactly the reverse
// propagation convention (see the GhostPath documentation).
func ghostParaxialRecurse(els []paraxialElement, zs []float64, p GhostPair, y, u, imageZ float64) (float64, float64) {
	n := len(els)
	cur := 0
	moveTo := func(k int) {
		y += (zs[k] - zs[cur]) * u
		cur = k
	}
	for k := 0; k < p.I; k++ {
		if k > 0 {
			moveTo(k)
		}
		e := els[k]
		eta := e.nB / e.nA
		u = eta*u + (eta-1)*e.c*y
	}
	moveTo(p.I)
	u = -u - 2*els[p.I].c*y
	for k := p.I - 1; k > p.J; k-- {
		moveTo(k)
		e := els[k]
		eta := e.nA / e.nB // reversed pass: the incident medium is the one behind
		u = eta*u + (eta-1)*e.c*y
	}
	moveTo(p.J)
	u = -u - 2*els[p.J].c*y
	for k := p.J + 1; k < n; k++ {
		moveTo(k)
		e := els[k]
		eta := e.nB / e.nA
		u = eta*u + (eta-1)*e.c*y
	}
	y += (imageZ - zs[n-1]) * u
	return y, u
}

// TestGhostThinLensParaxialPower verifies the paraxial ghost chain against a
// closed form derived from scratch: for a thin element in air of index n and
// vertex curvatures c₁, c₂ the four interactions at one plane (refraction,
// reflection, reflection, refraction) compose to
//
//	Δū = ū' − ū = −(3n−1)(c₁ − c₂) y   ⇒   Φ_ghost = (3n−1)(c₁ − c₂) ,
//
// so a collimated beam at height y₀ leaves with ū' = −Φ_ghost y₀ and crosses
// the axis at f_ghost = 1/Φ_ghost = R/(2(3n−1)) for R₁ = R, R₂ = −R — one
// seventh of the focal length R/(2(n−1)) at n = 1.5. The thin element is the
// first-order idealisation (zero thickness): a real ray cannot traverse a
// zero-thickness biconvex lens (the two caps leave no positive glass
// thickness), which is why this test is paraxial-only; the real-ray replay is
// validated on thick systems below and by TestGhostPlateAnalytic.
func TestGhostThinLensParaxialPower(t *testing.T) {
	const n = 1.5
	const R = 100.0
	thin := mustSystem(t, Air{},
		[]Surface{
			{Profile: Sphere{R: R}, Aperture: CircleAperture(10)},
			{Profile: Sphere{R: -R}, Aperture: CircleAperture(10)},
		},
		[]float64{0},
		[]Medium{Air{}, ConstantIndex(n), Air{}})
	phi := (3*n - 1) * (1/R - (-1 / R)) // 2(3n−1)/R
	fGhost := 1 / phi

	// The ghost focus and the 1/7 ratio of the thin lens of the same radii.
	p, err := thin.Paraxial(Wavelengthd, math.Inf(-1))
	if err != nil {
		t.Fatal(err)
	}
	fLens := p.EFL
	if math.Abs(fLens-R/(2*(n-1))) > 1e-9 {
		t.Fatalf("lens focal length %g, want %g", fLens, R/(2*(n-1)))
	}
	// EFL = −1/ū' of the marginal ray, so the ghost focus is at Φ_ghost⁻¹.
	checkRel(t, "thin-lens ghost focal length", fGhost, R/(2*(3*n-1)), 1e-15)
	checkRel(t, "f_ghost / f_lens", fGhost/fLens, (n-1)/(3*n-1), 1e-12)

	// The landing height at an arbitrary image plane, y(z) = y₀(1 − Φ z); at
	// z = f_ghost the collimated beam crosses the axis.
	const y0 = 2.5
	for _, z := range []float64{1, 5, 7, fGhost} {
		got, err := thin.ParaxialGhostLanding(GhostPair{I: 1, J: 0}, GhostParaxial{Y: y0}, z, Wavelengthd)
		if err != nil {
			t.Fatal(err)
		}
		want := y0 * (1 - phi*z)
		checkClose(t, "collimated ghost landing", got.Y, want, 1e-12*(1+math.Abs(want)))
		checkClose(t, "ghost landing is meridional", got.X, 0, 1e-15)
	}
	// The power does not depend on the launch slope, so a tilted collimated
	// beam keeps its slope and picks up the same power: y(z) = y₀ + (u₀ − Φ y₀) z.
	for _, u := range []float64{1e-3, 1e-2} {
		z := 7.0
		got, err := thin.ParaxialGhostLanding(GhostPair{I: 1, J: 0}, GhostParaxial{Y: y0, U: u}, z, Wavelengthd)
		if err != nil {
			t.Fatal(err)
		}
		checkClose(t, "tilted ghost landing", got.Y, y0+(u-phi*y0)*z, 1e-12*(1+math.Abs(y0)))
	}
}

// TestGhostPlateAnalytic verifies the reversed pass of the paraxial chain
// against an exact closed form. A plane-parallel plate of index n and
// thickness T has no power, so the ghost path (reflect at the rear, reflect at
// the front, out through the rear again) is pure geometry:
//
//	y(z) = y₀ + 3 T u₀/n + u₀ (z − T) ,
//
// which pins the three conventions of the reversed pass at once: the transfer
// of the backwards leg (z' − z is negative), the medium swap of the reversed
// refraction (η = n_after/n_before = n instead of 1/n — using 1/n would put
// u₀/n² into the last leg and miss the closed form by a factor n²), and the
// reflection sign flip. The real-ray replay converges to the same value as
// u₀³ (the first neglected order).
func TestGhostPlateAnalytic(t *testing.T) {
	const n = 1.7
	const T = 5.0
	const y0 = 0.0
	sys := plateSystem(t, n, T)
	path, err := sys.GhostPathOf(GhostPair{I: 1, J: 0})
	if err != nil {
		t.Fatal(err)
	}
	const zImg = 30.0
	analytic := func(u0 float64) float64 { return y0 + 3*T*u0/n + u0*(zImg-T) }

	for _, c := range []struct {
		u0, tol float64
	}{
		{1e-4, 1e-9}, // measured residual 2.9e-12 mm
		{1e-3, 1e-7}, // measured residual 2.9e-9 mm
	} {
		u0, tol := c.u0, c.tol
		got, err := sys.ParaxialGhostLanding(GhostPair{I: 1, J: 0}, GhostParaxial{Y: y0, U: u0}, zImg, Wavelengthd)
		if err != nil {
			t.Fatal(err)
		}
		checkClose(t, "plate paraxial ghost", got.Y, analytic(u0), 1e-15)
		r := NewRay(V(0, y0, 0), V(0, u0, 1), Wavelengthd)
		hit := sys.ReplayGhost(path, r, zImg)
		if !hit.OK() {
			t.Fatalf("plate ghost replay failed: %v", hit.Status)
		}
		checkStatus(t, "plate ghost status", hit.Status, StatusOK)
		if d := math.Abs(hit.Landing.Y - analytic(u0)); d > tol {
			t.Errorf("plate real ghost at u₀=%g: got %.12g, analytic %.12g (|diff| = %.3g > %.0e mm)",
				u0, hit.Landing.Y, analytic(u0), d, tol)
		}
	}
	// The residual is the first neglected order of the paraxial expansion, so
	// it must grow like u₀³.
	real := func(u0 float64) float64 {
		hit := sys.ReplayGhost(path, NewRay(V(0, y0, 0), V(0, u0, 1), Wavelengthd), zImg)
		return math.Abs(hit.Landing.Y - analytic(u0))
	}
	d4, d2 := real(1e-4), real(1e-2)
	if ratio := d2 / d4; ratio < 1e5 || ratio > 1e7 {
		t.Errorf("plate ghost residual ratio %.3g, want about 1e6 (cubic in u₀)", ratio)
	}
}

// TestGhostParaxialRepeatability checks the two independent transcriptions of
// the paraxial ghost chain — the y–u recurrence and the (y, ū) 2×2 matrix
// product — against each other for every pair of the doublet and for a range
// of launches. The acceptance threshold of the task is 1e-6; the two chains
// agree to a few 1e-15 mm, so the test asserts 1e-12 and the measured values
// are reported by the report of this task.
func TestGhostParaxialRepeatability(t *testing.T) {
	sys := ghostDoublet(t)
	p, err := sys.Paraxial(Wavelengthd, math.Inf(-1))
	if err != nil {
		t.Fatal(err)
	}
	els := sys.paraxialElements(Wavelengthd)
	zs := make([]float64, len(els))
	for i := range els {
		zs[i] = sys.VertexZ(i)
	}
	n0 := sys.Mediums[0].N(Wavelengthd)
	launches := [][2]float64{{1, 0}, {0.5, 1e-3}, {2, 0.01}, {-1.5, -2e-3}}
	worstY, worstU := 0.0, 0.0
	for _, pr := range GhostPairs(sys.NumSurfaces()) {
		m, err := sys.GhostABCD(pr, p.BackFocalZ, Wavelengthd)
		if err != nil {
			t.Fatalf("%s: GhostABCD: %v", pr, err)
		}
		// Every elementary matrix of the chain has determinant 1 (transfers
		// [[1,L/n],[0,1]], refraction and reflection [[1,0],[c,1]]), so the
		// symplectic invariant det M = 1 checks the composition.
		if d := math.Abs(m.Det() - 1); d > 1e-12 {
			t.Errorf("%s: det(reversed ABCD chain) = %.15g, want 1 (|diff| = %.3g)", pr, m.Det(), d)
		}
		// The landing reported by the analysis must be the chain applied to the
		// launch, and the independent y–u recurrence must agree with it.
		for _, lu := range launches {
			land, err := sys.ParaxialGhostLanding(pr, GhostParaxial{Y: lu[0], U: lu[1]}, p.BackFocalZ, Wavelengthd)
			if err != nil {
				t.Fatal(err)
			}
			want, _ := m.Apply(lu[0], n0*lu[1])
			checkClose(t, "landing from the ABCD chain", land.Y, want, 1e-12)

			y1, u1 := ghostParaxialRecurse(els, zs, pr, lu[0], lu[1], p.BackFocalZ)
			y2, u2, _ := ghostParaxialChain(els, zs, n0, pr, lu[0], lu[1], p.BackFocalZ)
			if d := math.Abs(y1 - y2); d > worstY {
				worstY = d
			}
			if d := math.Abs(u1 - u2); d > worstU {
				worstU = d
			}
			if d := math.Abs(y1 - y2); d > 1e-12 {
				t.Errorf("%s launch (%g, %g): recurrence %.15g vs matrices %.15g (|diff| = %.3g > 1e-12)",
					pr, lu[0], lu[1], y1, y2, d)
			}
		}
	}
	t.Logf("paraxial recurrence vs ABCD chain: max |Δy| = %.3g mm, max |Δu| = %.3g", worstY, worstU)
}

// TestGhostRealConvergesToParaxial is the high-precision version of the
// acceptance criterion (a): with a vanishing sampling aperture and a vanishing
// field the real-ray centroid of every ghost path of the doublet must converge
// to the paraxial landing point of its reference ray. The residual is the
// third-order aberration of the path at the reference ray, so the thresholds
// are absolute lengths: 1e-6 mm for the 0.05° field of this test.
func TestGhostRealConvergesToParaxial(t *testing.T) {
	sys := ghostDoublet(t)
	p, err := sys.Paraxial(Wavelengthd, math.Inf(-1))
	if err != nil {
		t.Fatal(err)
	}
	field := FieldTan(0.05*math.Pi/180, 0)
	pup := Pupil{Z: -5, Radius: 1e-3}
	rep, err := sys.GhostAnalysis(GhostOptions{
		Field: field, Pupil: pup, ImageZ: p.BackFocalZ, LambdaNm: Wavelengthd,
	})
	if err != nil {
		t.Fatal(err)
	}
	if len(rep.Paths) != 3 {
		t.Fatalf("got %d paths, want 3", len(rep.Paths))
	}
	for _, r := range rep.Paths {
		if r.Skipped() || r.Spot.N != r.Spot.NTotal {
			t.Fatalf("%s: %d of %d rays reached the image plane (%s)", r.Pair, r.Spot.N, r.Spot.NTotal, r.Note)
		}
		d := r.Spot.Centroid.Sub(r.Paraxial).Norm()
		if d > 1e-6 {
			t.Errorf("%s: real centroid %s vs paraxial %s differ by %.3g mm > 1e-6 mm",
				r.Pair, r.Spot.Centroid, r.Paraxial, d)
		}
		t.Logf("%s: |centroid − paraxial| = %.3g mm at a 1e-3 mm pupil, 0.05° field", r.Pair, d)
	}
}

// TestGhostOnAxisSymmetry checks that the ghost of an on-axis field keeps the
// rotational symmetry of the system: the paraxial landing point is the axis
// and the real-ray centroid is on the axis to machine precision.
func TestGhostOnAxisSymmetry(t *testing.T) {
	sys := ghostDoublet(t)
	p, err := sys.Paraxial(Wavelengthd, math.Inf(-1))
	if err != nil {
		t.Fatal(err)
	}
	rep, err := sys.GhostAnalysis(GhostOptions{
		Field: FieldTan(0, 0), Pupil: Pupil{Z: -5, Radius: 0.2}, ImageZ: p.BackFocalZ, LambdaNm: Wavelengthd,
		Points: Hexapolar(2),
	})
	if err != nil {
		t.Fatal(err)
	}
	for _, r := range rep.Paths {
		if r.Skipped() {
			t.Fatalf("%s: %s", r.Pair, r.Note)
		}
		if d := r.Paraxial.PerpNorm(); d > 1e-12 {
			t.Errorf("%s: paraxial landing %s is %.3g mm off the axis", r.Pair, r.Paraxial, d)
		}
		if d := r.Spot.Centroid.PerpNorm(); d > 1e-9 {
			t.Errorf("%s: on-axis centroid %s is %.3g mm off the axis", r.Pair, r.Spot.Centroid, d)
		}
	}
}
