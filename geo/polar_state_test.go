package geo

import (
	"errors"
	"math"
	"math/cmplx"
	"math/rand"
	"testing"
)

// ---------------------------------------------------------------------------
// Independent reference: a 3D boundary-value solver for one plane interface.
// It uses neither the Fresnel formulas nor the s/p basis of the package. The
// unknown reflected and transmitted fields are expanded in two orthonormal
// vectors each and fixed by the continuity of the tangential components of E
// and H = n·k̂×E (η0 = 1, e^{i(k·r−ωt)}), which is a 4×4 complex linear system.
// ---------------------------------------------------------------------------

// polRandUnit returns a random unit vector.
func polRandUnit(rng *rand.Rand) Vec3 {
	for {
		v := V(rng.NormFloat64(), rng.NormFloat64(), rng.NormFloat64())
		if u, ok := v.Unit(); ok {
			return u
		}
	}
}

// polRandTransverse returns a random complex field transverse to the unit d.
func polRandTransverse(rng *rand.Rand, d Vec3) cvec {
	var e cvec
	for k := range e {
		e[k] = complex(rng.NormFloat64(), rng.NormFloat64())
	}
	c := e.dot(d)
	return e.plus(polReal(d, -c))
}

// polCrossReal returns v×e for a real vector v and a complex field e.
func polCrossReal(v Vec3, e cvec) cvec {
	return cvec{
		complex(v.Y, 0)*e[2] - complex(v.Z, 0)*e[1],
		complex(v.Z, 0)*e[0] - complex(v.X, 0)*e[2],
		complex(v.X, 0)*e[1] - complex(v.Y, 0)*e[0],
	}
}

// polPerp returns two unit vectors completing the unit vector d to a triad. The
// construction is independent of TransverseFrame.
func polPerp(d Vec3) (Vec3, Vec3) {
	a := AxisX
	if math.Abs(d.X) > 0.9 {
		a = AxisY
	}
	u1 := d.Cross(a).Normalized()
	return u1, d.Cross(u1)
}

// polSolve4 solves A·x = b (complex 4×4) by Gaussian elimination with partial
// pivoting.
func polSolve4(A [4][4]complex128, b [4]complex128) ([4]complex128, bool) {
	for col := 0; col < 4; col++ {
		piv := col
		for r := col + 1; r < 4; r++ {
			if cmplx.Abs(A[r][col]) > cmplx.Abs(A[piv][col]) {
				piv = r
			}
		}
		if cmplx.Abs(A[piv][col]) < 1e-300 {
			return [4]complex128{}, false
		}
		A[col], A[piv] = A[piv], A[col]
		b[col], b[piv] = b[piv], b[col]
		for r := col + 1; r < 4; r++ {
			f := A[r][col] / A[col][col]
			for c := col; c < 4; c++ {
				A[r][c] -= f * A[col][c]
			}
			b[r] -= f * b[col]
		}
	}
	var x [4]complex128
	for r := 3; r >= 0; r-- {
		s := b[r]
		for c := r + 1; c < 4; c++ {
			s -= A[r][c] * x[c]
		}
		x[r] = s / A[r][r]
	}
	return x, true
}

// polRefSolve returns the reflected and transmitted electric fields (and the
// reflected and transmitted directions) of the plane wave e (transverse to din)
// at a plane interface with unit normal nrm, from a medium of real index n1 into
// a medium of real index n2. Total internal reflection is reported as an error.
func polRefSolve(din, nrm Vec3, n1, n2 float64, e cvec) (er, et cvec, dr, dt Vec3, err error) {
	dt, ok := Refract(din, nrm, n1, n2)
	if !ok {
		return er, et, dr, dt, errors.New("total internal reflection")
	}
	dr = din.Sub(nrm.Mul(2 * din.Dot(nrm)))
	u1, u2 := polPerp(dr)
	w1, w2 := polPerp(dt)
	t1, t2 := polPerp(nrm)
	us := [2]Vec3{u1, u2}
	ws := [2]Vec3{w1, w2}
	hi := polCrossReal(din, e) // H_inc = n1·(din×e)
	var A [4][4]complex128
	var b [4]complex128
	for k, tk := range [2]Vec3{t1, t2} {
		for j := 0; j < 2; j++ {
			A[k][j] = complex(us[j].Dot(tk), 0)
			A[k][2+j] = -complex(ws[j].Dot(tk), 0)
			A[2+k][j] = complex(n1*dr.Cross(us[j]).Dot(tk), 0)
			A[2+k][2+j] = -complex(n2*dt.Cross(ws[j]).Dot(tk), 0)
		}
		b[k] = -e.dot(tk)
		b[2+k] = -complex(n1, 0) * hi.dot(tk)
	}
	x, ok := polSolve4(A, b)
	if !ok {
		return er, et, dr, dt, errors.New("singular boundary-value system")
	}
	er = polReal(us[0], x[0]).plus(polReal(us[1], x[1]))
	et = polReal(ws[0], x[2]).plus(polReal(ws[1], x[3]))
	return er, et, dr, dt, nil
}

// polMaxDiff returns max_k |a_k − b_k|.
func polMaxDiff(a, b cvec) float64 {
	m := 0.0
	for k := range a {
		if d := cmplx.Abs(a[k] - b[k]); d > m {
			m = d
		}
	}
	return m
}

// TestPolTransverseFrame checks the transverse basis of TransverseFrame: it is
// orthonormal, transverse to d and right-handed (ex×ey = d), and on axis it is
// the global (x, y) frame.
func TestPolTransverseFrame(t *testing.T) {
	ex, ey := TransverseFrame(AxisZ)
	if ex != AxisX || ey != AxisY {
		t.Fatalf("on-axis frame = (%v, %v), want the global (x, y)", ex, ey)
	}
	rng := rand.New(rand.NewSource(3))
	for i := 0; i < 2000; i++ {
		d := polRandUnit(rng)
		ex, ey := TransverseFrame(d)
		polCheckReal(t, "|ex|", ex.Norm(), 1, 1e-12)
		polCheckReal(t, "|ey|", ey.Norm(), 1, 1e-12)
		polCheckReal(t, "ex·d", ex.Dot(d), 0, 1e-12)
		polCheckReal(t, "ey·d", ey.Dot(d), 0, 1e-12)
		polCheckReal(t, "ex·ey", ex.Dot(ey), 0, 1e-12)
		if r := ex.Cross(ey).Sub(d).Norm(); r > 1e-12 {
			t.Fatalf("ex×ey ≠ d: residual %.3g for d=%v", r, d)
		}
	}
	// Singular direction d ∥ ŷ falls back to ex = x̂, which is transverse to ŷ.
	ex, _ = TransverseFrame(AxisY)
	if ex != AxisX {
		t.Fatalf("fallback frame for d = ŷ: ex = %v, want x̂", ex)
	}
}

// TestPolInteractionMatchesReference compares the vector algebra of one interface
// (polInteract with the Fresnel coefficients) with the independent 3D boundary
// solver, for random directions and normals (both orientations, arbitrary skew
// geometry), random indices and random transverse fields. The residual is
// required to be below 1e-12, which is about 1e4 times the round-off level.
func TestPolInteractionMatchesReference(t *testing.T) {
	rng := rand.New(rand.NewSource(20261015 + 1))
	worst := 0.0
	cases := 0
	for i := 0; i < 6000; i++ {
		n1 := 1 + 0.8*rng.Float64()
		n2 := 1 + 1.2*rng.Float64()
		din := polRandUnit(rng)
		nrm := polRandUnit(rng)
		cosI := math.Abs(din.Dot(nrm))
		if cosI < 1e-3 {
			continue
		}
		e := polRandTransverse(rng, din)
		er, et, dr, dt, err := polRefSolve(din, nrm, n1, n2, e)
		if err != nil {
			continue // total internal reflection
		}
		a, _, err := FresnelAmplitudes(complex(n1, 0), complex(n2, 0), cosI)
		if err != nil {
			t.Fatal(err)
		}
		gotT, _ := polInteract(din, nrm, dt, false, a.TS, a.TP, e)
		gotR, _ := polInteract(din, nrm, dr, true, a.RS, a.RP, e)
		if d := polMaxDiff(gotT, et); d > worst {
			worst = d
		}
		if d := polMaxDiff(gotR, er); d > worst {
			worst = d
		}
		// The reference fields must themselves be transverse to their directions.
		if cmplx.Abs(et.dot(dt)) > 1e-12 || cmplx.Abs(er.dot(dr)) > 1e-12 {
			t.Fatalf("case %d: reference field not transverse to its direction", i)
		}
		cases++
	}
	if cases < 4000 {
		t.Fatalf("only %d valid cases", cases)
	}
	if worst > 1e-12 {
		t.Fatalf("worst deviation from the 3D reference %.3g exceeds 1e-12", worst)
	}
	t.Logf("worst deviation from the independent 3D boundary solver: %.3g over %d cases", worst, cases)
}

// TestPolPerfectMirrorTangentialField checks the metal-mirror model at random
// geometries: the total tangential field vanishes at the surface, E_inc + E_ref
// has no component along the two tangents of the surface, and the reflected field
// is transverse to the mirror direction.
func TestPolPerfectMirrorTangentialField(t *testing.T) {
	rng := rand.New(rand.NewSource(99))
	m := MetalMirrorAmplitudes()
	worst := 0.0
	for i := 0; i < 2000; i++ {
		din := polRandUnit(rng)
		nrm := polRandUnit(rng)
		if math.Abs(din.Dot(nrm)) < 1e-3 {
			continue
		}
		e := polRandTransverse(rng, din)
		dr := din.Sub(nrm.Mul(2 * din.Dot(nrm)))
		er, _ := polInteract(din, nrm, dr, true, m.RS, m.RP, e)
		t1, t2 := polPerp(nrm)
		for _, tk := range []Vec3{t1, t2} {
			if d := cmplx.Abs(e.dot(tk) + er.dot(tk)); d > worst {
				worst = d
			}
		}
		if d := cmplx.Abs(er.dot(dr)); d > worst {
			worst = d
		}
	}
	if worst > 1e-12 {
		t.Fatalf("tangential field at a perfect mirror: residual %.3g > 1e-12", worst)
	}
	t.Logf("perfect mirror: worst tangential residual %.3g", worst)
}

// TestPolVectorEnergyConservation checks power balance of the vector algebra at a
// lossless interface for arbitrary (also elliptical) polarization: the reflected
// and transmitted powers n·cosθ·|E|² add up to the incident power to 1e-12.
func TestPolVectorEnergyConservation(t *testing.T) {
	rng := rand.New(rand.NewSource(555))
	worst := 0.0
	cases := 0
	for i := 0; i < 4000; i++ {
		n1 := 1 + 0.8*rng.Float64()
		n2 := 1 + 1.2*rng.Float64()
		din := polRandUnit(rng)
		nrm := polRandUnit(rng)
		cosI := math.Abs(din.Dot(nrm))
		if cosI < 1e-3 || (n1/n2)*(n1/n2)*(1-cosI*cosI) >= 1 {
			continue
		}
		e := polRandTransverse(rng, din)
		a, _, err := FresnelAmplitudes(complex(n1, 0), complex(n2, 0), cosI)
		if err != nil {
			t.Fatal(err)
		}
		dr := din.Sub(nrm.Mul(2 * din.Dot(nrm)))
		dt, ok := Refract(din, nrm, n1, n2)
		if !ok {
			continue
		}
		er, _ := polInteract(din, nrm, dr, true, a.RS, a.RP, e)
		et, _ := polInteract(din, nrm, dt, false, a.TS, a.TP, e)
		pin := n1 * cosI * e.intensity()
		pr := n1 * cosI * er.intensity()
		pt := n2 * math.Abs(dt.Dot(nrm)) * et.intensity()
		if d := math.Abs(pr+pt-pin) / pin; d > worst {
			worst = d
		}
		cases++
	}
	if cases < 3000 {
		t.Fatalf("only %d cases", cases)
	}
	if worst > 1e-12 {
		t.Fatalf("relative power imbalance %.3g > 1e-12", worst)
	}
	t.Logf("worst relative power imbalance %.3g over %d random elliptical states", worst, cases)
}

// TestPolNormalIncidenceTrace checks the normal-incidence limit through the
// trace: on axis, the output Jones vector is t·(Jx, Jy) with t = 2n1/(n1+n2),
// and the transmitted power fraction is (n2/n1)·t² = 1 − R.
func TestPolNormalIncidenceTrace(t *testing.T) {
	const n1, n2 = 1.0, 1.5
	sys, err := NewSystem(Air{}, []Surface{NewSurface(Plane{}, 0)}, nil, []Medium{ConstantIndex(n1), ConstantIndex(n2)})
	if err != nil {
		t.Fatal(err)
	}
	ray := NewRay(V(0.1, -0.2, -3), AxisZ, 550)
	tt := 2 * n1 / (n1 + n2)
	R := (n1 - n2) * (n1 - n2) / ((n1 + n2) * (n1 + n2))
	for _, in := range [][2]complex128{{1, 0}, {0, 1}, {0.6, complex(0, 0.8)}} {
		res, err := sys.TracePolarized(ray, in, nil)
		if err != nil {
			t.Fatal(err)
		}
		polCheck(t, "Jones x", res.Jones[0], complex(tt, 0)*in[0], 1e-12)
		polCheck(t, "Jones y", res.Jones[1], complex(tt, 0)*in[1], 1e-12)
		polCheckReal(t, "transmittance", res.Transmittance, 1-R, 1e-12)
	}
}

// TestPolMirrorReversesHandedness checks a perfect mirror at normal incidence: the
// field is reversed, E_out = −E_in, and in the output frame of the reflected ray
// (ex × ey = d_out = −ẑ) the Jones vector is (Jx, −Jy). A left circular state
// therefore becomes right-handed, which is the physical helicity flip of a mirror
// seen by a viewer facing the reflected beam.
func TestPolMirrorReversesHandedness(t *testing.T) {
	sys, err := NewSystem(Air{}, []Surface{NewSurface(Plane{}, 0).AsMirror()}, nil, nil)
	if err != nil {
		t.Fatal(err)
	}
	ray := NewRay(V(0.1, 0.2, -5), AxisZ, 550)
	in := JonesCircular(true)
	res, err := sys.TracePolarized(ray, in, nil)
	if err != nil {
		t.Fatal(err)
	}
	polCheck(t, "reflected Jones x", res.Jones[0], in[0], 1e-12)
	polCheck(t, "reflected Jones y", res.Jones[1], -in[1], 1e-12)
	polCheck(t, "field = −E_in (x)", res.Field[0], -in[0], 1e-12)
	polCheck(t, "field = −E_in (y)", res.Field[1], -in[1], 1e-12)
	if res.Ellipse.Handedness != "右旋" {
		t.Errorf("reflected circular state handedness %q, want 右旋", res.Ellipse.Handedness)
	}
	polCheckReal(t, "mirror transmittance", res.Transmittance, 1, 1e-12)
}

// TestPolParallelPlateOblique is acceptance (d): a tilted ray through a
// plane-parallel plate (two refractions, the plate's faces parallel to each other
// so the plane of incidence does not change). The output is compared with two
// independent references: the closed form E_out = (1−r_s²)E_s + (1−r_p²)E_p, and
// the numerical 3D boundary solver applied to both faces. Required agreement
// 1e-9 (the achieved value is at the round-off level).
func TestPolParallelPlateOblique(t *testing.T) {
	const nGlass = 1.5168
	const thick = 2.5
	const wl = 587.5618
	sys, err := NewSystem(Air{}, []Surface{NewSurface(Plane{}, 0), NewSurface(Plane{}, 0)},
		[]float64{thick}, []Medium{Air{}, ConstantIndex(nGlass), Air{}})
	if err != nil {
		t.Fatal(err)
	}
	inputs := [][2]complex128{
		{1, 0},
		{0, 1},
		{complex(1/math.Sqrt2, 0), complex(1/math.Sqrt2, 0)},
		{0.6, complex(0, 0.8)},
		{complex(0.3, -0.5), complex(-0.2, 0.7)},
	}
	worstAn, worstNum, worstT := 0.0, 0.0, 0.0
	for _, theta := range []float64{0.05, 0.3, 0.6, 0.9, 1.2, 1.45} {
		d := V(0, math.Sin(theta), math.Cos(theta))
		ray := NewRay(V(0.25, -1.5, -4), d, wl)
		ex, ey := TransverseFrame(d)
		a, _, err := FresnelAmplitudes(1, complex(nGlass, 0), math.Cos(theta))
		if err != nil {
			t.Fatal(err)
		}
		As := 1 - a.RS*a.RS // t_s(12)·t_s(21) = 1 − r_s²
		Ap := 1 - a.RP*a.RP // t_p(12)·t_p(21) = 1 − r_p²
		s := d.Cross(AxisZ).Normalized()
		p := s.Cross(d)
		for _, in := range inputs {
			res, err := sys.TracePolarized(ray, in, nil)
			if err != nil {
				t.Fatal(err)
			}
			ein := polReal(ex, in[0]).plus(polReal(ey, in[1]))
			ean := polReal(s, As*ein.dot(s)).plus(polReal(p, Ap*ein.dot(p)))
			_, e1, _, dt1, err := polRefSolve(d, AxisZ, 1, nGlass, ein)
			if err != nil {
				t.Fatal(err)
			}
			_, e2, _, _, err := polRefSolve(dt1, AxisZ, nGlass, 1, e1)
			if err != nil {
				t.Fatal(err)
			}
			if x := polMaxDiff(cvec(res.Field), ean); x > worstAn {
				worstAn = x
			}
			if x := polMaxDiff(cvec(res.Field), e2); x > worstNum {
				worstNum = x
			}
			jAn := [2]complex128{ean.dot(ex), ean.dot(ey)}
			polCheck(t, "plate Jones x", res.Jones[0], jAn[0], 1e-9)
			polCheck(t, "plate Jones y", res.Jones[1], jAn[1], 1e-9)
			// OPL to the exit face: 4/cosθ in air (origin at z = −4) plus n·T/cosθ_g in glass.
			cosG := math.Sqrt(1 - math.Sin(theta)*math.Sin(theta)/(nGlass*nGlass))
			opl := 4/math.Cos(theta) + nGlass*thick/cosG
			wantPhase := 2 * math.Pi * opl / (wl * 1e-6)
			if d := math.Abs(res.OPLPhase - wantPhase); d > 1e-12*wantPhase {
				t.Fatalf("OPL phase %.17g, want %.17g (2π·OPL/λ with OPL = %.17g)", res.OPLPhase, wantPhase, opl)
			}
			tAn := ean.intensity() / ein.intensity()
			if x := math.Abs(res.Transmittance - tAn); x > worstT {
				worstT = x
			}
		}
	}
	if worstAn > 1e-9 || worstNum > 1e-9 || worstT > 1e-9 {
		t.Fatalf("parallel plate: closed form %.3g, numerical reference %.3g, transmittance %.3g (limit 1e-9)",
			worstAn, worstNum, worstT)
	}
	t.Logf("parallel plate: |field − closed form| ≤ %.3g, |field − 3D solver| ≤ %.3g, |ΔT| ≤ %.3g",
		worstAn, worstNum, worstT)
}

// skewLens is a thick biconvex lens with spherical surfaces.
func skewLens(t *testing.T) *System {
	t.Helper()
	sys, err := NewSystem(Air{},
		[]Surface{NewSurface(Sphere{R: 60}, 14), NewSurface(Sphere{R: -60}, 14)},
		[]float64{6}, []Medium{Air{}, ConstantIndex(1.6), Air{}})
	if err != nil {
		t.Fatal(err)
	}
	return sys
}

// TestPolSkewRayThroughLensChain checks a skew ray (off the meridional plane, so
// the plane of incidence turns from surface to surface) through a spherical lens.
// The polarized trace is compared, surface by surface, with the 3D boundary
// solver driven by the same geometry (directions, normals and indices of the
// trace). A fixed-basis approximation is also evaluated to show that the basis
// rotation changes the answer for this geometry.
func TestPolSkewRayThroughLensChain(t *testing.T) {
	sys := skewLens(t)
	ray := NewRay(V(3, -4, -30), V(-0.08, 0.1, 1), 550)
	ex, ey := TransverseFrame(ray.Dir)
	worstChain, worstNaive, maxTurn := 0.0, math.Inf(1), 0.0
	for _, in := range [][2]complex128{{1, 0}, {0, 1}, {0.6, complex(0, 0.8)}} {
		res, err := sys.TracePolarized(ray, in, nil)
		if err != nil {
			t.Fatalf("trace: %v", err)
		}
		if len(res.Steps) != 2 {
			t.Fatalf("expected 2 interactions, got %d", len(res.Steps))
		}
		for _, st := range res.Steps {
			// s, p_in, p_out form orthonormal triads with the directions of their interaction.
			polCheckReal(t, "|s|", st.S.Norm(), 1, 1e-12)
			polCheckReal(t, "s·n", st.S.Dot(st.Normal), 0, 1e-12)
			polCheckReal(t, "s·d_in", st.S.Dot(st.DirIn), 0, 1e-12)
			polCheckReal(t, "p_in·d_in", st.PIn.Dot(st.DirIn), 0, 1e-12)
			polCheckReal(t, "p_out·d_out", st.POut.Dot(st.DirOut), 0, 1e-12)
			if r := st.POut.Sub(st.S.Cross(st.DirOut)).Norm(); r > 1e-12 {
				t.Fatalf("p_out ≠ s×d_out (residual %.3g)", r)
			}
		}
		// The plane of incidence turns between the two surfaces.
		turn := math.Acos(math.Min(1, math.Abs(res.Steps[0].S.Dot(res.Steps[1].S))))
		if turn < 1e-3 {
			t.Fatalf("planes of incidence of the two surfaces coincide (%.3g rad)", turn)
		}
		if turn > maxTurn {
			maxTurn = turn
		}
		e := polReal(ex, in[0]).plus(polReal(ey, in[1]))
		for _, st := range res.Steps {
			_, et, _, _, err := polRefSolve(st.DirIn, st.Normal, st.NBefore, st.NAfter, e)
			if err != nil {
				t.Fatalf("reference at surface %d: %v", st.Surface, err)
			}
			e = et
		}
		if d := polMaxDiff(cvec(res.Field), e); d > worstChain {
			worstChain = d
		}
		// Fixed-basis approximation: the s/p basis of the first surface is reused
		// at the second surface.
		s0, p0 := res.Steps[0].S, res.Steps[0].PIn
		f := polReal(ex, in[0]).plus(polReal(ey, in[1]))
		for _, st := range res.Steps {
			f = polReal(s0, st.Amp.TS*f.dot(s0)).plus(polReal(p0, st.Amp.TP*f.dot(p0)))
		}
		if d := polMaxDiff(cvec(res.Field), f); d < worstNaive {
			worstNaive = d
		}
	}
	if worstChain > 1e-12 {
		t.Fatalf("skew-ray chain: deviation from the 3D reference %.3g > 1e-12", worstChain)
	}
	if worstNaive < 1e-3 {
		t.Fatalf("fixed-basis approximation differs by only %.3g; the geometry does not probe the basis rotation", worstNaive)
	}
	t.Logf("skew ray through a lens: |field − 3D chain| ≤ %.3g; plane of incidence turns by %.3g rad; fixed-basis approximation differs by ≥ %.3g",
		worstChain, maxTurn, worstNaive)
}

// TestPolSurfaceMatrixAccumulation checks the 2×2 accumulation: the product of the
// per-interaction matrices (and of the element matrices placed on the path) must
// reproduce the 3D field, for a skew ray through the lens and for a tilted plate
// with a quarter-wave plate placed behind its second face.
func TestPolSurfaceMatrixAccumulation(t *testing.T) {
	lens := skewLens(t)
	ray := NewRay(V(3, -4, -30), V(-0.08, 0.1, 1), 550)
	plate, err := NewSystem(Air{}, []Surface{NewSurface(Plane{}, 0), NewSurface(Plane{}, 0)},
		[]float64{2.5}, []Medium{Air{}, ConstantIndex(1.5168), Air{}})
	if err != nil {
		t.Fatal(err)
	}
	platRay := NewRay(V(0.25, -1.5, -4), V(0, math.Sin(0.8), math.Cos(0.8)), 587.5618)
	cases := []struct {
		sys *System
		ray Ray
		opt *PolOptions
	}{
		{lens, ray, nil},
		{plate, platRay, &PolOptions{Elements: []PolElementAt{{AfterSurface: 1, J: QuarterWaveJones(0.3)}}}},
	}
	worst := 0.0
	for _, c := range cases {
		m, err := c.sys.JonesMatrix(c.ray, c.opt)
		if err != nil {
			t.Fatal(err)
		}
		for _, in := range [][2]complex128{{1, 0}, {0, 1}, {0.6, complex(0, 0.8)}} {
			res, err := c.sys.TracePolarized(c.ray, in, c.opt)
			if err != nil {
				t.Fatal(err)
			}
			w := res.Matrix.Apply(in)
			for k := 0; k < 2; k++ {
				if d := cmplx.Abs(w[k] - res.Jones[k]); d > worst {
					worst = d
				}
			}
			for i := 0; i < 2; i++ {
				for k := 0; k < 2; k++ {
					if d := cmplx.Abs(res.Matrix[i][k] - m[i][k]); d > worst {
						worst = d
					}
				}
			}
		}
	}
	if worst > 1e-12 {
		t.Fatalf("surface-by-surface matrix accumulation differs from the 3D field by %.3g > 1e-12", worst)
	}
	t.Logf("2×2 accumulation vs 3D field and vs JonesMatrix: worst |Δ| = %.3g", worst)
}

// TestPolJonesMatrixLinearity checks that the system Jones matrix reproduces the
// output of an arbitrary input vector, i.e. the polarized trace is linear in the
// input, and that the columns are the outputs of the basis vectors.
func TestPolJonesMatrixLinearity(t *testing.T) {
	sys := skewLens(t)
	ray := NewRay(V(3, -4, -30), V(-0.08, 0.1, 1), 550)
	m, err := sys.JonesMatrix(ray, nil)
	if err != nil {
		t.Fatal(err)
	}
	rng := rand.New(rand.NewSource(11))
	for i := 0; i < 20; i++ {
		v := [2]complex128{complex(rng.NormFloat64(), rng.NormFloat64()), complex(rng.NormFloat64(), rng.NormFloat64())}
		res, err := sys.TracePolarized(ray, v, nil)
		if err != nil {
			t.Fatal(err)
		}
		w := m.Apply(v)
		polCheck(t, "M·v (x)", res.Jones[0], w[0], 1e-12)
		polCheck(t, "M·v (y)", res.Jones[1], w[1], 1e-12)
	}
}

// TestPolTraceErrors covers the failure modes: total internal reflection, a
// vignetted ray and an invalid ray. The partial result is kept with the error.
func TestPolTraceErrors(t *testing.T) {
	// Glass to air at 60° (beyond the 41.8° critical angle of n = 1.5): TIR.
	// A plate entered from air never gives TIR at its back face, so the ray
	// starts inside the glass.
	tir, err := NewSystem(ConstantIndex(1.5), []Surface{NewSurface(Plane{}, 0)}, nil,
		[]Medium{ConstantIndex(1.5), Air{}})
	if err != nil {
		t.Fatal(err)
	}
	res, err := tir.TracePolarized(NewRay(V(0, -2, -2), V(0, math.Sin(math.Pi/3), math.Cos(math.Pi/3)), 550),
		[2]complex128{1, 0}, nil)
	if err == nil {
		t.Errorf("expected a TIR error, got a complete trace (%v)", res.Trace)
	}
	// Vignetted at the second surface.
	small, err := NewSystem(Air{}, []Surface{NewSurface(Plane{}, 0.5), NewSurface(Plane{}, 0.5)},
		[]float64{1}, []Medium{Air{}, ConstantIndex(1.5), Air{}})
	if err != nil {
		t.Fatal(err)
	}
	vign, err := small.TracePolarized(NewRay(V(0.9, 0, -2), AxisZ, 550), [2]complex128{1, 0}, nil)
	if err == nil {
		t.Errorf("expected a vignetting error")
	}
	if len(vign.Steps) != 0 {
		t.Errorf("a ray vignetted at the first surface must report no completed step")
	}
	// Invalid ray (zero direction).
	if _, err := small.TracePolarized(Ray{Origin: V(0, 0, -2), WavelengthNm: 550}, [2]complex128{1, 0}, nil); err == nil {
		t.Errorf("expected an error for an invalid ray")
	}
}
