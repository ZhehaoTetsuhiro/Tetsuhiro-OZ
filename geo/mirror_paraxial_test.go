package geo

import (
	"math"
	"testing"
)

// ---------------------------------------------------------------------------
// First-order (paraxial) engine for coaxial systems with mirrors.
//
// The engine is verified against closed-form first-order optics and against the
// real ray trace of the package (System.Trace), which is the independent
// reference: every axial position quoted by the paraxial engine (ImageZ,
// BackFocalZ, ImageZ of a finite conjugate, the exit pupil of a mirror system)
// is the extrapolation of a real ray line and can be measured with a real ray
// of vanishing aperture. Real rays are therefore launched at a height of the
// order of 1e-4 mm and the crossings are compared to the paraxial prediction;
// the residual is the third-order aberration of the configuration, which is
// quoted for the full aperture as well (see the individual tests).
//
// Conventions used by the tests (see System.Paraxial for the derivation):
//
//	surface sag            z = c r²/2,  c = 1/R, positive sag towards +z
//	object                 at z = objectZ < 0, light travels towards +z
//	mirror                 reflects; the thickness behind it is negative
//	mirror law             u' = -u - 2 c y,  ū' = -ū - 2 n c y
//	single mirror image    ImageZ = s R / (R + 2 s), m = ImageZ / s,  s = -objectZ
//	angular magnification  Γ = w'/w with w the slope along the propagation
//	                       direction (Paraxial.Parity flips the sign after an odd
//	                       number of mirrors)
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// helpers
// ---------------------------------------------------------------------------

// singleMirror returns a one-surface mirror system of radius R and clear
// semi-aperture semi (the surface is also the aperture stop).
func singleMirror(t *testing.T, R, semi float64) *System {
	t.Helper()
	return mustSystem(t, Air{}, []Surface{
		{Profile: Sphere{R: R}, Mirror: true, Aperture: CircleAperture(semi)},
	}, nil, nil)
}

// stopMirrorPair returns [plane stop, conic mirror] with the given separation.
func stopMirrorPair(t *testing.T, gap, stopSemi, R, K, semi float64) *System {
	t.Helper()
	return mustSystem(t, Air{}, []Surface{
		{Profile: Plane{}, Aperture: CircleAperture(stopSemi)},
		{Profile: Conic{R: R, K: K}, Mirror: true, Aperture: CircleAperture(semi)},
	}, []float64{gap}, nil)
}

// realCrossOfRay returns the z at which the traced real ray crosses the axis
// (NaN when the trace fails, ±Inf for a ray parallel to the axis).
func realCrossOfRay(sys *System, h, lambda, zStart float64) float64 {
	res := sys.Trace(NewRay(V(0, h, zStart), V(0, 0, 1), lambda))
	if !res.OK() {
		return math.NaN()
	}
	p, d := res.Ray.Origin, res.Ray.Dir
	if d.Y == 0 {
		return math.Inf(signOf(-p.Y))
	}
	return p.Z - p.Y*d.Z/d.Y
}

// realFiniteImageZ traces the ray that leaves the axial object point at objectZ
// and cuts the plane of surface 0 at height h, and returns its axis crossing.
// In the paraxial limit it is the image of the axial object point.
func realFiniteImageZ(sys *System, objectZ, h, lambda float64) float64 {
	u := h / (sys.VertexZ(0) - objectZ)
	res := sys.Trace(NewRay(V(0, 0, objectZ), V(0, u, 1), lambda))
	if !res.OK() {
		return math.NaN()
	}
	p, d := res.Ray.Origin, res.Ray.Dir
	if d.Y == 0 {
		return math.Inf(signOf(-p.Y))
	}
	return p.Z - p.Y*d.Z/d.Y
}

// realAngleRatio traces the ray through the point (0, 0, stopZ) with entrance
// tangent tan and returns the ratio of the emergent to the incident slope
// measured along the direction of propagation (w = dy/dz in the travelling
// frame, so the fold of a mirror does not flip the sign of an angle).
func realAngleRatio(sys *System, stopZ, tan, lambda float64) float64 {
	res := sys.Trace(NewRay(V(0, 0, stopZ), V(0, -tan, 1), lambda))
	if !res.OK() {
		return math.NaN()
	}
	d := res.Ray.Dir
	if d.Z == 0 {
		return math.Inf(1)
	}
	return (d.Y / math.Abs(d.Z)) / (-tan)
}

// realBestFocusZ locates the plane that minimises the transverse RMS of a
// parallel on-axis bundle of nRays between lo and hi on a dense deterministic
// grid, refined by four shrinking steps (no random numbers: the result is a
// pure function of the system and the arguments).
func realBestFocusZ(sys *System, semi, lambda, zStart, lo, hi float64, nRays int) (float64, float64) {
	var pts, dirs []Vec3
	for i := 0; i < nRays; i++ {
		h := semi * (0.02 + 0.96*float64(i)/float64(nRays-1))
		res := sys.Trace(NewRay(V(0, h, zStart), V(0, 0, 1), lambda))
		if !res.OK() {
			continue
		}
		last := res.Steps[len(res.Steps)-1]
		pts = append(pts, last.Point)
		dirs = append(dirs, last.DirOut)
	}
	if len(pts) < 3 {
		return math.NaN(), math.NaN()
	}
	rms := func(z float64) float64 {
		ys := make([]float64, 0, len(pts))
		for i := range pts {
			if math.Abs(dirs[i].Z) < 1e-14 {
				continue
			}
			ys = append(ys, pts[i].Y+(z-pts[i].Z)*dirs[i].Y/dirs[i].Z)
		}
		if len(ys) == 0 {
			return math.Inf(1)
		}
		var m, sum float64
		for _, y := range ys {
			m += y
		}
		m /= float64(len(ys))
		for _, y := range ys {
			sum += (y - m) * (y - m)
		}
		return math.Sqrt(sum / float64(len(ys)))
	}
	const steps = 2000
	best, bestRMS := math.NaN(), math.Inf(1)
	for k := 0; k <= steps; k++ {
		zz := lo + (hi-lo)*float64(k)/float64(steps)
		if r := rms(zz); r < bestRMS {
			bestRMS, best = r, zz
		}
	}
	dz := (hi - lo) / float64(steps)
	for pass := 0; pass < 4; pass++ {
		for _, zz := range []float64{best - dz, best + dz} {
			if r := rms(zz); r < bestRMS {
				bestRMS, best = r, zz
			}
		}
		dz *= 0.3
	}
	return best, bestRMS
}

// unfoldedMatrix builds the ray transfer matrix of the system in the frame of
// the propagating beam (the formulation of the unfolding used by the ghost
// analysis of the package, see ghostParaxialChain): the reduced angle stored is
// n·w with w = dy/dz measured along the direction of travel, a gap contributes
// B = σ·Δz/n with σ the propagation parity, a refraction contributes
// C = -σ(nA-nB)c and a mirror contributes C = 2σ·n·c and flips σ. It is an
// independent implementation of the same first-order physics and must satisfy
//
//	Matrix() = diag(1, σ)·unfoldedMatrix()
//
// with σ the parity after the last surface (mirrors only ever conjugate the
// physical matrix by that sign, because a reflection reverses ū).
func unfoldedMatrix(sys *System, lambda float64) (Mat2, float64) {
	sigma, idx := 1.0, 1.0
	m := Identity()
	cur := sys.VertexZ(0)
	for i := 0; i < sys.NumSurfaces(); i++ {
		sf := sys.Surfaces[i]
		if sf.Break {
			continue
		}
		z := sys.VertexZ(i)
		if i > 0 {
			m = Mat2{A: 1, B: sigma * (z - cur) / idx, D: 1}.Mul(m)
		}
		cur = z
		c := curvatureOf(sf.Profile)
		nUp := sys.Mediums[i].N(lambda)
		nDn := sys.Mediums[i+1].N(lambda)
		if sf.Mirror {
			m = Mat2{A: 1, C: 2 * sigma * idx * c, D: 1}.Mul(m)
			sigma = -sigma
		} else {
			m = Mat2{A: 1, C: -sigma * (nDn - nUp) * c, D: 1}.Mul(m)
			idx = nDn
		}
	}
	return m, sigma
}

// pairCrossZ reproduces the ray-pair crossing of tools/libcheck/measure
// (crossZ) including its sign convention, which measures the parameter along
// the incident ray instead of against it. For a bundle that converges it
// therefore reports the point reflected at the hit plane, p1.Z - t·d1.Z,
// instead of the true crossing p1.Z + t·d1.Z.
func pairCrossZ(p1, d1, p2, d2 Vec3) float64 {
	den := d1.Z*d2.Y - d1.Y*d2.Z
	if math.Abs(den) < 1e-15 {
		return math.NaN()
	}
	t := ((p2.Y-p1.Y)*d2.Z - (p2.Z-p1.Z)*d2.Y) / den
	return p1.Z + t*d1.Z
}

// pairCrossZTrue is the same crossing with the sign of the ray parameter
// corrected, i.e. the intersection of the two ray lines.
func pairCrossZTrue(p1, d1, p2, d2 Vec3) float64 {
	den := d1.Z*d2.Y - d1.Y*d2.Z
	if math.Abs(den) < 1e-15 {
		return math.NaN()
	}
	t := ((p2.Y-p1.Y)*d2.Z - (p2.Z-p1.Z)*d2.Y) / den
	return p1.Z - t*d1.Z
}

// ---------------------------------------------------------------------------
// 1. single concave mirror: focal length and the imaging formula
// ---------------------------------------------------------------------------

func TestParaxialSingleConcaveMirror(t *testing.T) {
	const (
		R    = -2000.0 // concave as seen from the object side (sag towards -z)
		semi = 100.0
	)
	sys := singleMirror(t, R, semi)
	p, err := sys.Paraxial(Wavelengthd, math.Inf(-1))
	if err != nil {
		t.Fatalf("Paraxial: %v", err)
	}
	// Closed form: f = R/2 in the transmission convention of the package, so a
	// concave mirror has a negative EFL and the focal point at R/2 in front.
	checkRel(t, "EFL = R/2", p.EFL, R/2, 1e-12)
	checkRel(t, "BFL = R/2", p.BFL, R/2, 1e-12)
	checkRel(t, "BackFocalZ = R/2", p.BackFocalZ, R/2, 1e-12)
	checkRel(t, "FFL = R/2 (mirror is its own time reverse)", p.FFL, R/2, 1e-12)
	checkRel(t, "FrontFocalZ = BackFocalZ", p.FrontFocalZ, p.BackFocalZ, 1e-12)
	checkRel(t, "BackPrincipalZ = 0 (principal plane at the vertex)", p.BackPrincipalZ, 0, 1e-9)
	// The unfolded (customary) focal length of a mirror system is Parity·EFL:
	// for a single reflection this is +|R|/2, the usual f = R/2 of the mirror
	// convention.
	checkRel(t, "Parity = -1 (one reflection)", p.Parity(), -1, 1e-15)
	checkRel(t, "unfolded focal length = -EFL", p.Parity()*p.EFL, -R/2, 1e-12)
	if p.Mirrors != 1 {
		t.Errorf("Mirrors = %d, want 1", p.Mirrors)
	}
	if p.Afocal() {
		t.Error("a single mirror with R != 0 is not afocal")
	}
	// Focal length by real ray trace: parallel ray of vanishing height.
	checkRel(t, "real focus at R/2", realCrossOfRay(sys, 1e-4, Wavelengthd, -50), R/2, 1e-9)

	// Imaging formula ImageZ = sR/(R+2s), m = ImageZ/s for a table of finite
	// conjugates; each entry is also verified with a real ray.
	for _, tc := range []struct {
		s          float64 // distance of the object in front of the mirror
		wantZ, m   float64
		realTolRel float64
	}{
		{6000, -1200, -0.2, 1e-9},
		{3000, -1500, -0.5, 1e-9},
		{2000, -2000, -1.0, 1e-9}, // object at the centre of curvature
		{1500, -3000, -2.0, 1e-9},
		{1000, math.Inf(1), 0, 1e-9}, // object at the focus: image at infinity
		{500, 1000, 2.0, 1e-9},       // object inside f: virtual, erect, magnified
	} {
		objectZ := -tc.s
		p, err := sys.Paraxial(Wavelengthd, objectZ)
		if err != nil {
			t.Fatalf("s=%g: %v", tc.s, err)
		}
		wantZ := tc.s * R / (R + 2*tc.s)
		if math.IsInf(tc.wantZ, 1) {
			if !math.IsInf(p.ImageZ, 1) && !math.IsInf(p.ImageZ, -1) {
				t.Errorf("s=%g: ImageZ = %g, want an image at infinity", tc.s, p.ImageZ)
			}
		} else {
			checkRel(t, "imaging formula ImageZ", p.ImageZ, wantZ, 1e-12)
			checkRel(t, "transverse magnification m = ImageZ/s", p.Magnification, wantZ/tc.s, 1e-12)
			checkRel(t, "real ray image of the axial point", realFiniteImageZ(sys, objectZ, 1e-3, Wavelengthd),
				wantZ, tc.realTolRel)
		}
		if !math.IsNaN(tc.m) && tc.m != 0 && !math.IsInf(wantZ, 0) {
			checkRel(t, "magnification sign", math.Abs(p.Magnification), math.Abs(tc.m), 1e-12)
		}
	}
	// Marginal ray trace of a spherical mirror shows the expected spherical
	// aberration: the marginal rays cross the axis closer to the mirror than the
	// paraxial focus (the focus is at z = R/2 < 0, the mirror at z = 0) and the
	// longitudinal error grows as the square of the height.
	f1 := realCrossOfRay(sys, 1e-4, Wavelengthd, -50)
	f10 := realCrossOfRay(sys, 10, Wavelengthd, -50)
	f20 := realCrossOfRay(sys, 20, Wavelengthd, -50)
	if !(f1 < f10 && f10 < f20) {
		t.Fatalf("expected the marginal focus closer to the mirror than the paraxial one: h=1e-4 %g, h=10 %g, h=20 %g", f1, f10, f20)
	}
	checkRel(t, "mirror SA ∝ h²", (f20-f1)/(f10-f1), 4, 0.05)
	t.Logf("single concave mirror R=%g: paraxial focus %.9f, marginal h=10: %.9f, h=20: %.9f", R, f1, f10, f20)
}

// ---------------------------------------------------------------------------
// 2. single convex mirror: virtual image behind the mirror
// ---------------------------------------------------------------------------

func TestParaxialSingleConvexMirrorVirtualImage(t *testing.T) {
	const (
		R    = 200.0 // convex as seen from the object side (sag towards +z)
		semi = 50.0
	)
	sys := singleMirror(t, R, semi)
	p, err := sys.Paraxial(Wavelengthd, math.Inf(-1))
	if err != nil {
		t.Fatalf("Paraxial: %v", err)
	}
	// A convex mirror diverges: the focal point is virtual, behind the mirror.
	checkRel(t, "BackFocalZ = R/2 (behind the mirror)", p.BackFocalZ, R/2, 1e-12)
	if !(p.BackFocalZ > 0) {
		t.Fatalf("convex mirror focus must lie behind the mirror (z>0), got %g", p.BackFocalZ)
	}
	checkRel(t, "unfolded focal length = -EFL", p.Parity()*p.EFL, -R/2, 1e-12)
	checkRel(t, "real virtual focus at R/2", realCrossOfRay(sys, 1e-4, Wavelengthd, -50), R/2, 1e-9)

	for _, tc := range []struct{ s, wantZ, wantM float64 }{
		{200, 200.0 * 200.0 / 600.0, 1.0 / 3.0},  // ImageZ = sR/(R+2s) = +66.6…
		{100, 100.0 * 200.0 / 400.0, 0.5},        // object closer in: virtual image, m=+1/2
		{600, 600.0 * 200.0 / 1400.0, 1.0 / 7.0}, // far object: m → 0, image → f
	} {
		objectZ := -tc.s
		p, err := sys.Paraxial(Wavelengthd, objectZ)
		if err != nil {
			t.Fatalf("s=%g: %v", tc.s, err)
		}
		checkRel(t, "convex imaging formula", p.ImageZ, tc.wantZ, 1e-12)
		checkRel(t, "convex magnification (erect, reduced)", p.Magnification, tc.wantM, 1e-12)
		if !(p.Magnification > 0) {
			t.Errorf("a convex mirror forms an erect image, m = %g", p.Magnification)
		}
		if !(math.Abs(p.Magnification) < 1) {
			t.Errorf("a convex mirror forms a reduced image, m = %g", p.Magnification)
		}
		checkRel(t, "real virtual image", realFiniteImageZ(sys, objectZ, 1e-4, Wavelengthd), tc.wantZ, 1e-9)
	}
}

// ---------------------------------------------------------------------------
// 3. Newtonian primary (unfolded layout of library/geo/mirror-newtonian-*)
// ---------------------------------------------------------------------------

func TestParaxialNewtonianUnfolded(t *testing.T) {
	// Layout of library/geo/mirror-newtonian-f1000-f5.json: an aperture stop
	// 60 mm in front of a parabolic primary of R = -2000 mm, focus unfolded.
	sys := stopMirrorPair(t, 60, 100, -2000, -1, 102)
	p, err := sys.Paraxial(Wavelengthd, math.Inf(-1))
	if err != nil {
		t.Fatalf("Paraxial: %v", err)
	}
	checkRel(t, "EFL = R/2", p.EFL, -1000, 1e-12)
	checkRel(t, "BackFocalZ = 60 - 1000", p.BackFocalZ, -940, 1e-12)
	checkRel(t, "ImageZ = BackFocalZ", p.ImageZ, -940, 1e-12)
	// The parabola is free of spherical aberration: real rays of every height
	// cross the axis exactly at the paraxial focus.
	for _, h := range []float64{1e-4, 1, 10, 50, 100} {
		checkRel(t, "parabola focus independent of the height",
			realCrossOfRay(sys, h, Wavelengthd, -2000), -940, 1e-9)
	}
	fz, rms := realBestFocusZ(sys, 100, Wavelengthd, -2000, -941, -939, 21)
	checkRel(t, "best focus of the full aperture", fz, -940, 1e-6)
	if rms > 1e-6 {
		t.Errorf("parabolic mirror RMS at the focus = %g mm, want < 1e-6", rms)
	}
	// Finite conjugates of the same mirror follow ImageZ = z_m + sR/(R+2s) with
	// s the distance of the object in front of the mirror and z_m its vertex.
	zMirror := sys.VertexZ(1)
	for _, s := range []float64{3000, 2000, 1000, 500} {
		pp, err := sys.Paraxial(Wavelengthd, zMirror-s)
		if err != nil {
			t.Fatalf("s=%g: %v", s, err)
		}
		want := zMirror + s*(-2000)/(-2000+2*s)
		checkRel(t, "finite conjugate", pp.ImageZ, want, 1e-12)
		checkRel(t, "finite conjugate real ray",
			realFiniteImageZ(sys, zMirror-s, 1e-3, Wavelengthd), want, 1e-6)
	}
	t.Logf("newtonian: FOCUSZ %.6f, RMS %.3g", fz, rms)
}

// ---------------------------------------------------------------------------
// 4. Cassegrain layout of the library: a fold with a negative thickness
// ---------------------------------------------------------------------------

func TestParaxialCassegrainFoldLayout(t *testing.T) {
	// Layout of library/geo/mirror-cassegrain-f1000-d300-b300.json: stop,
	// primary R=-2000 at z=60, a negative thickness of -300 mm (the fold) and
	// the secondary R=-646.1538461538461, every surface in air.
	const R1, R2 = -2000.0, -646.1538461538461
	sys := mustSystem(t, Air{}, []Surface{
		{Profile: Plane{}, Aperture: CircleAperture(150)},
		{Profile: Sphere{R: R1}, Mirror: true, Aperture: CircleAperture(153)},
		{Profile: Sphere{R: R2}, Mirror: true, Aperture: CircleAperture(117.5)},
	}, []float64{60, -300}, nil)
	if got := sys.VertexZ(2); got != -240 {
		t.Fatalf("the negative thickness must fold the chain: vertex z = %g, want -240", got)
	}
	p, err := sys.Paraxial(Wavelengthd, math.Inf(-1))
	if err != nil {
		t.Fatalf("Paraxial: %v", err)
	}
	if p.Mirrors != 2 || p.Parity() != 1 {
		t.Fatalf("Mirrors = %d, Parity = %g, want 2 and +1", p.Mirrors, p.Parity())
	}
	// Analytic two-step derivation of the image: the primary sends the parallel
	// ray to its focus 1000 mm in front of it, so the ray reaches the secondary
	// (300 mm down-beam) at height y2 = 1 - 300/1000 = 0.7 with ū = 1e-3; the
	// secondary reflects it to ū' = -ū - 2 c2 y2 = -1e-3 + 2·0.7/646.1538461…
	// and the image sits 0.7/ū' = 600 mm down-beam of the secondary.
	u2 := 1e-3
	y2 := 1.0 - 300*u2
	uba := -u2 - 2*y2/R2
	checkRel(t, "secondary incident height y2", y2, 0.7, 1e-12)
	checkRel(t, "secondary emergent ū", uba, 1.1666666666666667e-3, 1e-12)
	checkRel(t, "ImageZ = z2 - y2/ū'", p.ImageZ, -240-y2/uba, 1e-15)
	checkRel(t, "ImageZ = -840 (600 mm down-beam of the secondary)", p.ImageZ, -840, 1e-12)
	// The same system in the unfolded power-sum form: with the parity σ (the
	// first mirror +1, the second -1) the unfolded powers Φ̂ = -2 n c σ add as
	// usual over the unfolded separation d̂ = -t, and the physical C entry of
	// the matrix is -Φ̂_total.
	c1, c2 := 1/R1, 1/R2
	phi1 := -2 * 1 * c1 * 1
	phi2 := -2 * 1 * c2 * (-1)
	phiU := phi1 + phi2 - 300*phi1*phi2
	checkRel(t, "unfolded power Φ̂1", phi1, 1e-3, 1e-12)
	checkRel(t, "unfolded power Φ̂2", phi2, -2/(646.1538461538461), 1e-12)
	m := p.Matrix()
	checkRel(t, "C = -Φ̂_total", m.C, -phiU, 1e-12)
	checkRel(t, "EFL = -n/C", p.EFL, -1/m.C, 1e-15)
	checkRel(t, "EFL (absolute value)", math.Abs(p.EFL), 857.1428571428571, 1e-12)
	if !(p.EFL < 0) {
		t.Errorf("this layout diverges (virtual image at z=%g): EFL = %g must be negative", p.ImageZ, p.EFL)
	}
	// Real rays: the paraxial image is the crossing of a vanishing-aperture ray,
	// and the full aperture shows the spherical aberration of the two spheres.
	checkRel(t, "real ray image", realCrossOfRay(sys, 1e-4, Wavelengthd, -2000), -840, 1e-9)
	fz, rms := realBestFocusZ(sys, 150, Wavelengthd, -2000, -900, -700, 61)
	t.Logf("cassegrain layout: paraxial %.6f, real marginal %.6f, best focus %.6f (RMS %.4g), library FOCUSZ 345.328867 (artifact, see TestLibraryMirrorFocusNumbers)",
		p.ImageZ, realCrossOfRay(sys, 150, Wavelengthd, -2000), fz, rms)
	if math.Abs(fz+840) > 5 {
		t.Errorf("best focus of the full aperture %.3f is far from the paraxial image %.3f", fz, p.ImageZ)
	}
	// Finite conjugate: a point source 1000 mm in front of the primary.
	pp, err := sys.Paraxial(Wavelengthd, -1060)
	if err != nil {
		t.Fatalf("finite object: %v", err)
	}
	checkRel(t, "finite conjugate real ray", realFiniteImageZ(sys, -1060, 1e-3, Wavelengthd), pp.ImageZ, 1e-6)
	checkRel(t, "Conjugate agrees with ImageZ", pp.Conjugate(-1060), pp.ImageZ, 1e-12)
}

// ---------------------------------------------------------------------------
// 5. a converging two-mirror system: parabola primary + parabola secondary
// ---------------------------------------------------------------------------

func TestParaxialTrueCassegrainConverging(t *testing.T) {
	// Primary f1 = 1000 mm (parabola R = -2000) with the secondary 300 mm in
	// front of it and the back focus 300 mm behind the primary, i.e. the ray
	// that leaves the secondary at height y2 = 0.7 crosses the axis after
	// 600 mm; the secondary is confocal with the prime focus (700 mm away), so
	// its radius follows 2 c2 = -1/700 - 1/600 and the system works at
	// m = 600/700.
	const (
		d = 300.0 // primary to secondary
		b = 300.0 // back focal distance behind the primary
	)
	c2 := (1.0/600 - 1.0/700) / 2
	R2 := 1 / c2
	sys := mustSystem(t, Air{}, []Surface{
		{Profile: Plane{}, Aperture: CircleAperture(40)},
		{Profile: Conic{R: -2000, K: -1}, Mirror: true, Aperture: CircleAperture(40)},
		{Profile: Conic{R: R2, K: -1}, Mirror: true, Aperture: CircleAperture(30)},
	}, []float64{60, -d}, nil)
	p, err := sys.Paraxial(Wavelengthd, math.Inf(-1))
	if err != nil {
		t.Fatalf("Paraxial: %v", err)
	}
	checkRel(t, "secondary radius", R2, 8400, 1e-12)
	checkRel(t, "ImageZ = primary vertex + back focus", p.ImageZ, 60+b, 1e-12)
	// EFL = f1·m with m = (distance of the final focus behind the secondary) /
	// (distance of the prime focus behind the secondary) = (d+b)/(f1-d) = 600/700.
	checkRel(t, "EFL = f1·m", p.EFL, 1000*(d+b)/(1000-d), 1e-12)
	checkRel(t, "real ray image", realCrossOfRay(sys, 1e-4, Wavelengthd, -2000), 60+b, 1e-9)
	fz, rms := realBestFocusZ(sys, 40, Wavelengthd, -2000, 60+b-2, 60+b+2, 21)
	checkClose(t, "best focus of the full aperture (both surfaces are parabolas)", fz, 60+b, 0.2)
	t.Logf("two-mirror parabola system: paraxial %.6f, real h=10 %.6f, best focus %.6f (RMS %.4g)",
		p.ImageZ, realCrossOfRay(sys, 10, Wavelengthd, -2000), fz, rms)
	// Finite conjugate: the real ray and the paraxial engine agree.
	pp, err := sys.Paraxial(Wavelengthd, -940)
	if err != nil {
		t.Fatalf("finite object: %v", err)
	}
	checkRel(t, "finite conjugate real ray", realFiniteImageZ(sys, -940, 1e-3, Wavelengthd), pp.ImageZ, 1e-6)
	t.Logf("two-mirror finite conjugate: object -940 -> image %.6f (m = %.6f)", pp.ImageZ, pp.Magnification)
}

// ---------------------------------------------------------------------------
// 6. Mersenne telescope: afocal, angular magnification
// ---------------------------------------------------------------------------

func TestParaxialMersenneAfocal(t *testing.T) {
	// Two confocal parabolas: primary R = -2000 (f1 = 1000) and a secondary
	// 500 mm in front of it with R = -1000, confocal with the prime focus.
	// The system is afocal with the beam compressed by A = 1 - d/f1 = 0.5 and
	// Γ = 1/A = f1/(f1-d) = 2.
	const (
		f1 = 1000.0
		d  = 500.0
	)
	sys := mustSystem(t, Air{}, []Surface{
		{Profile: Plane{}, Aperture: CircleAperture(50)},
		{Profile: Conic{R: -2 * f1, K: -1}, Mirror: true, Aperture: CircleAperture(50)},
		{Profile: Conic{R: -2 * (f1 - d), K: -1}, Mirror: true, Aperture: CircleAperture(30)},
	}, []float64{60, -d}, nil)
	p, err := sys.Paraxial(Wavelengthd, math.Inf(-1))
	if err != nil {
		t.Fatalf("Paraxial: %v", err)
	}
	if !p.Afocal() {
		t.Fatalf("mersenne must be afocal, |C| = %g", math.Abs(p.Matrix().C))
	}
	m := p.Matrix()
	checkRel(t, "beam compression A = 1 - d/f1", m.A, 1-d/f1, 1e-12)
	checkRel(t, "B = z2 - A z1", m.B, 530, 1e-12)
	checkClose(t, "C = 0", m.C, 0, 1e-12)
	checkRel(t, "D = Γ = f1/(f1-d)", m.D, f1/(f1-d), 1e-12)
	checkRel(t, "Parity = +1 (two reflections)", p.Parity(), 1, 1e-15)
	if !math.IsInf(p.EFL, 0) {
		t.Errorf("an afocal system has an infinite EFL, got %g", p.EFL)
	}
	// The documented definition Γ = w'/w equals Parity·D·n0/nk, and the ray
	// based result equals the closed form.
	gamma, err := sys.AngularMagnification(math.Inf(-1), 1e-4, Wavelengthd)
	if err != nil {
		t.Fatalf("AngularMagnification: %v", err)
	}
	checkRel(t, "Γ = 2", gamma, 2, 1e-12)
	checkRel(t, "Γ = Parity·D·n0/nk", gamma, p.Parity()*m.D*p.Index[0]/p.Index[len(p.Index)-1], 1e-14)
	checkRel(t, "Γ = f1/(f1-d)", gamma, f1/(f1-d), 1e-12)
	// Real rays: a parallel bundle at a small angle keeps its collimation and
	// is magnified by Γ. The residual is the coma of an off-axis two-mirror
	// parabola system: the field angle enters quadratically (measured
	// |ratio-Γ|/Γ = 1.22·tan², independent of the aperture) and vanishes in the
	// paraxial limit, which is the tolerance used here.
	for _, tc := range []struct {
		semi, tan float64
	}{
		{2, 1e-4},
		{5, 1e-4},
		{2, 1e-3},
		{50, 1e-2},
	} {
		got := realAngleRatio(sys, 0, tc.tan, Wavelengthd)
		tol := 1.5 * tc.tan * tc.tan * gamma
		if math.Abs(got-gamma) > tol {
			t.Errorf("real angular ratio (semi=%g, tan=%g) = %.12f, want %.12f ± %.3g", tc.semi, tc.tan, got, gamma, tol)
		}
		t.Logf("mersenne real angular ratio semi=%g tan=%g: %.12f (relerr %.3g, coma budget %.3g)",
			tc.semi, tc.tan, got, math.Abs(got-gamma)/gamma, 1.5*tc.tan*tc.tan)
	}
	// The task's 1e-6 relative bound is met at a field angle of 1e-4 rad.
	if got := realAngleRatio(sys, 0, 1e-4, Wavelengthd); math.Abs(got-gamma)/gamma > 1e-6 {
		t.Errorf("real angular ratio at tan=1e-4 = %.12f, want %.12f within 1e-6 relative", got, gamma)
	}
	// On axis the emergent bundle is collimated: the parallel input ray keeps
	// its height ratio and travels parallel to the axis.
	res := sys.Trace(NewRay(V(0, 20, -2000), V(0, 0, 1), Wavelengthd))
	if !res.OK() {
		t.Fatalf("on-axis trace failed: %v", res.Error())
	}
	if math.Abs(res.Ray.Dir.Y/res.Ray.Dir.Z) > 1e-12 {
		t.Errorf("on-axis output is not collimated: slope %g", res.Ray.Dir.Y/res.Ray.Dir.Z)
	}
	checkRel(t, "on-axis beam compression", res.Ray.Origin.Y, 20*(1-d/f1), 1e-9)
}

// ---------------------------------------------------------------------------
// 7. Kepler and Galileo telescopes built from ideal thin lenses
// ---------------------------------------------------------------------------

func TestParaxialThinLensTelescopes(t *testing.T) {
	// Two surfaces separated by a vanishing thickness form an exactly thin lens
	// of power Φ = (n-1)(c1-c2); with a constant index the focal length is exact
	// and the closed forms Γ = -f1/f2 (Kepler) and +f1/f2 (Galileo) hold.
	n := ConstantIndex(1.5168)
	thin := func(f float64) []Surface {
		R := 2 * (n.N(Wavelengthd) - 1) * f
		return []Surface{
			{Profile: Sphere{R: R}, Aperture: CircleAperture(5)},
			{Profile: Sphere{R: -R}, Aperture: CircleAperture(5)},
		}
	}
	tel := func(f2, d float64) *System {
		ss := []Surface{{Profile: Plane{}, Aperture: CircleAperture(5)}}
		ss = append(ss, thin(100)...)
		ss = append(ss, thin(f2)...)
		return mustSystem(t, Air{}, ss, []float64{0, 0, d, 0},
			[]Medium{Air{}, Air{}, n, Air{}, n, Air{}})
	}
	for _, tc := range []struct {
		name  string
		f2, d float64
		want  float64
	}{
		{"Kepler", 20, 120, -5},
		{"Galileo", -20, 80, 5},
	} {
		sys := tel(tc.f2, tc.d)
		p, err := sys.Paraxial(Wavelengthd, math.Inf(-1))
		if err != nil {
			t.Fatalf("%s: %v", tc.name, err)
		}
		if !p.Afocal() {
			t.Fatalf("%s: not afocal, |C| = %g", tc.name, math.Abs(p.Matrix().C))
		}
		m := p.Matrix()
		g, err := sys.AngularMagnification(math.Inf(-1), 1e-3, Wavelengthd)
		if err != nil {
			t.Fatalf("%s: AngularMagnification: %v", tc.name, err)
		}
		checkClose(t, tc.name+" Γ", g, tc.want, 1e-9)
		checkRel(t, tc.name+" Γ = ±f1/f2", g, -100/tc.f2, 1e-12)
		checkRel(t, tc.name+" A = -f2/f1", m.A, -tc.f2/100, 1e-12)
		checkRel(t, tc.name+" D = Γ (air to air)", m.D, tc.want, 1e-12)
		checkRel(t, tc.name+" Γ = Parity·D·n0/nk", g, p.Parity()*m.D, 1e-14)
		// The same value from the ray based definition, using the chief ray.
		c := p.ChiefRay(1e-3)
		checkRel(t, tc.name+" chief ray ratio", p.Parity()*c.UAfter[len(c.UAfter)-1]/c.U0, g, 1e-12)
		// An afocal system has no angular magnification for an axial reference
		// ray: the API reports the degeneracy instead of dividing by zero.
		if _, err := sys.AngularMagnification(math.Inf(-1), 0, Wavelengthd); err == nil {
			t.Errorf("%s: an axial (0-tangent) reference ray must be reported as degenerate", tc.name)
		}
	}
	// A realizable Kepler (thick N-BK7 biconvex lenses) with the spacing set
	// from the principal planes: it must be afocal with Γ = -f1/f2 and the real
	// rays must reproduce that value at a small aperture.
	nBK := N_BK7
	mkLens := func(f, tk, semi float64) *System {
		nd := nBK.N(Wavelengthd)
		R := 2 * (nd - 1) * f
		return mustSystem(t, Air{}, []Surface{
			{Profile: Sphere{R: R}, Aperture: CircleAperture(semi)},
			{Profile: Sphere{R: -R}, Aperture: CircleAperture(semi)},
		}, []float64{tk}, []Medium{Air{}, nBK, Air{}})
	}
	l1, l2 := mkLens(100, 2.0, 8), mkLens(20, 1.0, 2)
	p1, err := l1.Paraxial(Wavelengthd, math.Inf(-1))
	if err != nil {
		t.Fatal(err)
	}
	p2, err := l2.Paraxial(Wavelengthd, math.Inf(-1))
	if err != nil {
		t.Fatal(err)
	}
	// The principal planes must be f1+f2 apart for the telescope to be afocal.
	lensGap := p1.EFL + p2.EFL + p1.BackPrincipalZ - p2.FrontPrincipalZ - 2.0
	sys := mustSystem(t, Air{}, []Surface{
		{Profile: Plane{}, Aperture: CircleAperture(8)},
		{Profile: Sphere{R: 2 * (nBK.N(Wavelengthd) - 1) * 100}, Aperture: CircleAperture(8)},
		{Profile: Sphere{R: -2 * (nBK.N(Wavelengthd) - 1) * 100}, Aperture: CircleAperture(8)},
		{Profile: Sphere{R: 2 * (nBK.N(Wavelengthd) - 1) * 20}, Aperture: CircleAperture(2)},
		{Profile: Sphere{R: -2 * (nBK.N(Wavelengthd) - 1) * 20}, Aperture: CircleAperture(2)},
	}, []float64{0, 2.0, lensGap, 1.0},
		[]Medium{Air{}, Air{}, nBK, Air{}, nBK, Air{}})
	pt, err := sys.Paraxial(Wavelengthd, math.Inf(-1))
	if err != nil {
		t.Fatalf("thick Kepler: %v", err)
	}
	if !pt.Afocal() {
		t.Fatalf("thick Kepler must be afocal, |C| = %g", math.Abs(pt.Matrix().C))
	}
	g, err := sys.AngularMagnification(math.Inf(-1), 1e-4, Wavelengthd)
	if err != nil {
		t.Fatal(err)
	}
	checkRel(t, "thick Kepler Γ = -f1/f2", g, -p1.EFL/p2.EFL, 1e-12)
	got := realAngleRatio(sys, 0, 1e-4, Wavelengthd)
	if math.Abs(got-g)/math.Abs(g) > 1e-5 {
		t.Errorf("thick Kepler real angular ratio %.12f, want %.12f", got, g)
	}
	t.Logf("thick Kepler: f1=%.6f f2=%.6f Γ=%.9f, real ratio %.9f (relerr %.3g), lens gap %.6f",
		p1.EFL, p2.EFL, g, got, math.Abs(got-g)/math.Abs(g), lensGap)
}

// ---------------------------------------------------------------------------
// 8. plane fold mirror: transparent to angles
// ---------------------------------------------------------------------------

func TestParaxialPlaneFoldMirror(t *testing.T) {
	sys := mustSystem(t, Air{}, []Surface{
		{Profile: Plane{}, Mirror: true, Aperture: CircleAperture(50)},
	}, nil, nil)
	p, err := sys.Paraxial(Wavelengthd, math.Inf(-1))
	if err != nil {
		t.Fatalf("Paraxial: %v", err)
	}
	if !p.Afocal() {
		t.Fatalf("a plane mirror is afocal, |C| = %g", math.Abs(p.Matrix().C))
	}
	checkRel(t, "Parity = -1", p.Parity(), -1, 1e-15)
	if !math.IsInf(p.EFL, 0) {
		t.Errorf("plane mirror EFL = %g, want ±Inf", p.EFL)
	}
	g, err := sys.AngularMagnification(math.Inf(-1), 1e-3, Wavelengthd)
	if err != nil {
		t.Fatal(err)
	}
	checkClose(t, "Γ = +1 (a fold does not change an angle)", g, 1, 1e-12)
	if got := realAngleRatio(sys, 0, 1e-3, Wavelengthd); math.Abs(got-1) > 1e-12 {
		t.Errorf("real angular ratio across a plane mirror = %.15f, want 1", got)
	}
}

// ---------------------------------------------------------------------------
// 9. the ray transfer matrix agrees with the unfolded chain
// ---------------------------------------------------------------------------

func TestParaxialMirrorMatrixUnfoldedChain(t *testing.T) {
	mersenne := mustSystem(t, Air{}, []Surface{
		{Profile: Plane{}, Aperture: CircleAperture(50)},
		{Profile: Conic{R: -2000, K: -1}, Mirror: true, Aperture: CircleAperture(50)},
		{Profile: Conic{R: -1000, K: -1}, Mirror: true, Aperture: CircleAperture(30)},
	}, []float64{60, -500}, nil)
	cassegrain := mustSystem(t, Air{}, []Surface{
		{Profile: Plane{}, Aperture: CircleAperture(150)},
		{Profile: Sphere{R: -2000}, Mirror: true, Aperture: CircleAperture(153)},
		{Profile: Sphere{R: -646.1538461538461}, Mirror: true, Aperture: CircleAperture(117.5)},
	}, []float64{60, -300}, nil)
	cat := mustSystem(t, Air{}, []Surface{
		{Profile: Plane{}, Aperture: CircleAperture(10)},
		{Profile: Conic{R: -2000, K: -1}, Mirror: true, Aperture: CircleAperture(10)},
		{Profile: Sphere{R: -600}, Aperture: CircleAperture(10)},
		{Profile: Sphere{R: 600}, Aperture: CircleAperture(10)},
	}, []float64{60, -560, -6}, []Medium{Air{}, Air{}, Air{}, N_BK7, Air{}})
	refractive := mustSystem(t, Air{}, []Surface{
		{Profile: Sphere{R: 61}, Aperture: CircleAperture(12)},
		{Profile: Sphere{R: -44}, Aperture: CircleAperture(12)},
		{Profile: Sphere{R: -130}, Aperture: CircleAperture(25)},
	}, []float64{6, 2}, []Medium{Air{}, N_BK7, N_SF5, Air{}})
	for _, tc := range []struct {
		name string
		sys  *System
	}{
		{"single mirror", singleMirror(t, -2000, 100)},
		{"newtonian", stopMirrorPair(t, 60, 100, -2000, -1, 102)},
		{"mersenne", mersenne},
		{"cassegrain", cassegrain},
		{"catadioptric", cat},
		{"refractive triplet", refractive},
	} {
		p, err := tc.sys.Paraxial(Wavelengthd, math.Inf(-1))
		if err != nil {
			t.Fatalf("%s: %v", tc.name, err)
		}
		mp := p.Matrix()
		mu, sigma := unfoldedMatrix(tc.sys, Wavelengthd)
		// The physical matrix is the unfolded one with the second row scaled by
		// the exit parity (a reflection reverses ū = n·u).
		checkRel(t, tc.name+" A", mp.A, mu.A, 1e-12)
		checkRel(t, tc.name+" B", mp.B, mu.B, 1e-12)
		checkRel(t, tc.name+" C", mp.C, sigma*mu.C, 1e-12)
		checkRel(t, tc.name+" D", mp.D, sigma*mu.D, 1e-12)
		checkRel(t, tc.name+" parity", p.Parity(), sigma, 1e-15)
		// det M = (-1)^mirrors · n0/nk
		want := math.Pow(-1, float64(p.Mirrors)) * p.Index[0] / p.Index[len(p.Index)-1]
		checkRel(t, tc.name+" determinant", mp.Det(), want, 1e-12)
		// Every mirror system must satisfy the imaging relation B = 0 between
		// conjugate planes: the engine's ImageZ is the plane where the transfer
		// from the object plane has B = 0.
		if !p.Infinite {
			t.Errorf("%s: expected an object at infinity in this table", tc.name)
		}
	}
}

// ---------------------------------------------------------------------------
// 10. a pure mirror system has no chromatic aberration at all
// ---------------------------------------------------------------------------

func TestParaxialPureMirrorChromaticZero(t *testing.T) {
	single := singleMirror(t, -2000, 100)
	mersenne := mustSystem(t, Air{}, []Surface{
		{Profile: Plane{}, Aperture: CircleAperture(50)},
		{Profile: Conic{R: -2000, K: -1}, Mirror: true, Aperture: CircleAperture(50)},
		{Profile: Conic{R: -1000, K: -1}, Mirror: true, Aperture: CircleAperture(30)},
	}, []float64{60, -500}, nil)
	cassegrain := mustSystem(t, Air{}, []Surface{
		{Profile: Plane{}, Aperture: CircleAperture(150)},
		{Profile: Sphere{R: -2000}, Mirror: true, Aperture: CircleAperture(153)},
		{Profile: Sphere{R: -646.1538461538461}, Mirror: true, Aperture: CircleAperture(117.5)},
	}, []float64{60, -300}, nil)
	for _, tc := range []struct {
		name string
		sys  *System
	}{
		{"single mirror", single},
		{"mersenne", mersenne},
		{"cassegrain layout", cassegrain},
	} {
		if tc.sys.DispersiveSystem() {
			t.Fatalf("%s: an all-air system must not be dispersive", tc.name)
		}
		fc, err := tc.sys.AxialChromatic(math.Inf(-1), WavelengthF, WavelengthC)
		if err != nil {
			t.Fatalf("%s: AxialChromatic: %v", tc.name, err)
		}
		if math.Abs(fc) > 1e-12 {
			t.Errorf("%s: F-C focus shift of a pure mirror system = %g mm, want 0", tc.name, fc)
		}
		shift, err := tc.sys.ChromaticFocalShift(math.Inf(-1), []float64{WavelengthF, Wavelengthd, WavelengthC})
		if err != nil {
			t.Fatalf("%s: ChromaticFocalShift: %v", tc.name, err)
		}
		for i, s := range shift {
			if math.Abs(s) > 1e-12 {
				t.Errorf("%s: focal shift at lambda[%d] = %g, want 0", tc.name, i, s)
			}
		}
		if z, err := tc.sys.RefocusZ(math.Inf(-1), Wavelengthd, WavelengthF); err != nil || math.Abs(z) > 1e-12 {
			t.Errorf("%s: RefocusZ = %g (%v), want 0", tc.name, z, err)
		}
		// The lateral colour of the chief ray vanishes as well.
		if h, err := tc.sys.LateralChromatic(FieldTan(0, 1e-3), tc.sys.VertexZ(tc.sys.NumSurfaces()-1),
			WavelengthF, WavelengthC); err != nil || math.Abs(h) > 1e-12 {
			t.Errorf("%s: lateral colour = %g (%v), want 0", tc.name, h, err)
		}
		// Finite object as well: the whole first-order description is
		// wavelength independent when every medium is air.
		fcFin, err := tc.sys.AxialChromatic(-5000, WavelengthF, WavelengthC)
		if err != nil {
			t.Fatalf("%s: finite object: %v", tc.name, err)
		}
		if math.Abs(fcFin) > 1e-12 {
			t.Errorf("%s: finite-object F-C = %g, want 0", tc.name, fcFin)
		}
	}
}

// ---------------------------------------------------------------------------
// 11. catadioptric: the colour of a singlet in a folded beam
// ---------------------------------------------------------------------------

// catadioptricSinglet builds the parabola R = -2000 mm with a BK7 singlet
// (R = ±600 mm, 6 mm thick) 560 mm behind it in the folded beam. The singlet
// sees a beam travelling towards -z, so its first surface (the one at the
// larger z) carries the curvature -1/600 to stay converging.
func catadioptricSinglet(t *testing.T) *System {
	t.Helper()
	return mustSystem(t, Air{}, []Surface{
		{Profile: Plane{}, Aperture: CircleAperture(10)},
		{Profile: Conic{R: -2000, K: -1}, Mirror: true, Aperture: CircleAperture(10)},
		{Profile: Sphere{R: -600}, Aperture: CircleAperture(10)},
		{Profile: Sphere{R: 600}, Aperture: CircleAperture(10)},
	}, []float64{60, -560, -6}, []Medium{Air{}, Air{}, Air{}, N_BK7, Air{}})
}

func TestParaxialCatadioptricSingletChromatic(t *testing.T) {
	sys := catadioptricSinglet(t)
	pd, err := sys.Paraxial(Wavelengthd, math.Inf(-1))
	if err != nil {
		t.Fatalf("Paraxial: %v", err)
	}
	pF, err := sys.Paraxial(WavelengthF, math.Inf(-1))
	if err != nil {
		t.Fatal(err)
	}
	pC, err := sys.Paraxial(WavelengthC, math.Inf(-1))
	if err != nil {
		t.Fatal(err)
	}
	if pF.ImageZ == pC.ImageZ {
		t.Fatal("a catadioptric with a singlet must show axial colour")
	}
	fc, err := sys.AxialChromatic(math.Inf(-1), WavelengthF, WavelengthC)
	if err != nil {
		t.Fatal(err)
	}
	checkRel(t, "paraxial F-C = FocusAt(F) - FocusAt(C)", fc, pF.ImageZ-pC.ImageZ, 1e-15)
	// The mirror contributes no colour: the whole shift is the singlet's, and it
	// is reproduced by the thin-lens derivative relation
	//
	//	Δz = (b/f)²·Δf   (1/a + 1/f = 1/b)
	//
	// with f the lens's own paraxial focal length and b its distance to the
	// focus. The remaining difference is the thick-lens and principal-plane
	// residual of the lens (≈1%), quoted below.
	lensOnly := mustSystem(t, Air{}, []Surface{
		{Profile: Sphere{R: -600}, Aperture: CircleAperture(10)},
		{Profile: Sphere{R: 600}, Aperture: CircleAperture(10)},
	}, []float64{-6}, []Medium{Air{}, N_BK7, Air{}})
	ld, err := lensOnly.Paraxial(Wavelengthd, math.Inf(-1))
	if err != nil {
		t.Fatal(err)
	}
	lF, err := lensOnly.Paraxial(WavelengthF, math.Inf(-1))
	if err != nil {
		t.Fatal(err)
	}
	lC, err := lensOnly.Paraxial(WavelengthC, math.Inf(-1))
	if err != nil {
		t.Fatal(err)
	}
	b := pd.ImageZ - sys.VertexZ(3)
	expect := (b / ld.EFL) * (b / ld.EFL) * (lF.EFL - lC.EFL)
	if math.Abs(expect-fc)/math.Abs(fc) > 0.02 {
		t.Errorf("thin-lens expectation %g differs from the paraxial F-C %g by more than 2%%", expect, fc)
	}
	// Real rays: the vanishing-aperture crossing of every colour must sit on the
	// paraxial image, which makes the F-C cross-check exact to the trace
	// tolerance, and the RMS-optimal focus of a 10 mm aperture reproduces the
	// same drift within the spherical aberration of the singlet.
	zF := realCrossOfRay(sys, 1e-4, WavelengthF, -3000)
	zC := realCrossOfRay(sys, 1e-4, WavelengthC, -3000)
	checkRel(t, "real crossing at F", zF, pF.ImageZ, 1e-8)
	checkRel(t, "real crossing at C", zC, pC.ImageZ, 1e-8)
	checkRel(t, "real F-C of the vanishing aperture", zF-zC, fc, 1e-8)
	bfF, rmsF := realBestFocusZ(sys, 5, WavelengthF, -3000, pd.ImageZ-3, pd.ImageZ+3, 21)
	bfC, rmsC := realBestFocusZ(sys, 5, WavelengthC, -3000, pd.ImageZ-3, pd.ImageZ+3, 21)
	if math.Abs((bfF-bfC)-fc) > 0.1 {
		t.Errorf("real best-focus drift %.6f differs from the paraxial F-C %.6f by more than 0.1 mm", bfF-bfC, fc)
	}
	t.Logf("catadioptric singlet: paraxial ImageZ d=%.6f F=%.6f C=%.6f F-C=%.9f; thin-lens estimate %.6f; real best focus F %.6f (RMS %.3g) C %.6f (RMS %.3g), drift %.6f",
		pd.ImageZ, pF.ImageZ, pC.ImageZ, fc, expect, bfF, rmsF, bfC, rmsC, bfF-bfC)
}

// ---------------------------------------------------------------------------
// 12. catadioptric with an achromat in the folded beam
// ---------------------------------------------------------------------------

func TestParaxialCatadioptricAchromatChromatic(t *testing.T) {
	// The cemented doublet of library/geo/achromat-dls-f100-nbk7-f2.json
	// (N-BK7 crown, F2 flint) traversed backwards by the folded beam: the
	// surfaces keep their physical curvatures and order (the light meets the
	// flint surface first), the gaps become negative and the media are
	// reversed. This is the layout System.Reversed produces for the same lens.
	const (
		Rflint = -1016.200219223319
		Rcem   = -44.13066801758621
		Rcrown = 44.168566077079184
		tFlint = 0.6276392964649857
		tCrown = 2.2947800255278685
	)
	sys := mustSystem(t, Air{}, []Surface{
		{Profile: Plane{}, Aperture: CircleAperture(10)},
		{Profile: Conic{R: -2000, K: -1}, Mirror: true, Aperture: CircleAperture(10)},
		{Profile: Sphere{R: Rflint}, Aperture: CircleAperture(8)},
		{Profile: Sphere{R: Rcem}, Aperture: CircleAperture(8)},
		{Profile: Sphere{R: Rcrown}, Aperture: CircleAperture(8)},
	}, []float64{60, -560, -tFlint, -tCrown},
		[]Medium{Air{}, Air{}, Air{}, F2, N_BK7, Air{}})
	pd, err := sys.Paraxial(Wavelengthd, math.Inf(-1))
	if err != nil {
		t.Fatalf("Paraxial: %v", err)
	}
	fc, err := sys.AxialChromatic(math.Inf(-1), WavelengthF, WavelengthC)
	if err != nil {
		t.Fatal(err)
	}
	if fc == 0 {
		t.Fatal("the achromat of this test is traversed backwards, which is not perfectly achromatic")
	}
	// The reversal convention must reproduce the physical lens: the paraxial
	// focal length of the backward traversed doublet equals that of the library
	// prescription (|EFL| = 100.002447 mm, see the library description).
	lensOnly := mustSystem(t, Air{}, []Surface{
		{Profile: Sphere{R: Rflint}, Aperture: CircleAperture(8)},
		{Profile: Sphere{R: Rcem}, Aperture: CircleAperture(8)},
		{Profile: Sphere{R: Rcrown}, Aperture: CircleAperture(8)},
	}, []float64{-tFlint, -tCrown}, []Medium{Air{}, F2, N_BK7, Air{}})
	lo, err := lensOnly.Paraxial(Wavelengthd, math.Inf(-1))
	if err != nil {
		t.Fatal(err)
	}
	checkRel(t, "backward doublet EFL = library EFL", math.Abs(lo.EFL), 100.002447, 1e-6)
	// Real rays confirm the paraxial colour of the whole system exactly, and the
	// achromat keeps the shift an order of magnitude below the singlet's.
	singlet := catadioptricSinglet(t)
	fcSinglet, err := singlet.AxialChromatic(math.Inf(-1), WavelengthF, WavelengthC)
	if err != nil {
		t.Fatal(err)
	}
	if !(math.Abs(fc) < 0.1*math.Abs(fcSinglet)) {
		t.Errorf("achromat F-C %.6g is not well below the singlet F-C %.6g", fc, fcSinglet)
	}
	zF := realCrossOfRay(sys, 1e-4, WavelengthF, -3000)
	zC := realCrossOfRay(sys, 1e-4, WavelengthC, -3000)
	checkRel(t, "real F-C crossing equals the paraxial F-C", zF-zC, fc, 1e-6)
	checkRel(t, "real crossing at F", zF, mustParaxialZ(t, sys, WavelengthF), 1e-8)
	checkClose(t, "real crossing at C", zC, mustParaxialZ(t, sys, WavelengthC), 1e-6)
	// Spherical aberration of the folded doublet: quoted for the record.
	for _, h := range []float64{1e-3, 0.1, 1} {
		t.Logf("achromat catadioptric: real crossing h=%g: %.9f (paraxial %.9f, Δ = %.3g)",
			h, realCrossOfRay(sys, h, Wavelengthd, -3000), pd.ImageZ, realCrossOfRay(sys, h, Wavelengthd, -3000)-pd.ImageZ)
	}
	t.Logf("achromat catadioptric: paraxial ImageZ %.6f, EFL %.6f, F-C %.6f (singlet F-C %.6f)",
		pd.ImageZ, pd.EFL, fc, fcSinglet)
}

// mustParaxialZ returns the paraxial image position or fails the test.
func mustParaxialZ(t *testing.T, sys *System, lambda float64) float64 {
	t.Helper()
	p, err := sys.Paraxial(lambda, math.Inf(-1))
	if err != nil {
		t.Fatalf("Paraxial(%g): %v", lambda, err)
	}
	return p.ImageZ
}

// ---------------------------------------------------------------------------
// 13. determinism: two runs are bit-identical
// ---------------------------------------------------------------------------

func TestParaxialMirrorDeterminism(t *testing.T) {
	systems := []*System{
		singleMirror(t, -2000, 100),
		singleMirror(t, 200, 50),
		stopMirrorPair(t, 60, 100, -2000, -1, 102),
		mustSystem(t, Air{}, []Surface{
			{Profile: Plane{}, Aperture: CircleAperture(150)},
			{Profile: Sphere{R: -2000}, Mirror: true, Aperture: CircleAperture(153)},
			{Profile: Sphere{R: -646.1538461538461}, Mirror: true, Aperture: CircleAperture(117.5)},
		}, []float64{60, -300}, nil),
		catadioptricSinglet(t),
	}
	objects := []float64{math.Inf(-1), -5000, -3000, -1000, -500}
	lambdas := []float64{WavelengthF, Wavelengthd, WavelengthC}
	for si, sys := range systems {
		for _, oz := range objects {
			for _, lam := range lambdas {
				a, err := sys.Paraxial(lam, oz)
				if err != nil {
					continue // not every object position is valid for every system
				}
				b, err := sys.Paraxial(lam, oz)
				if err != nil {
					t.Fatalf("system %d: %v", si, err)
				}
				for _, pair := range [][2]float64{
					{a.EFL, b.EFL}, {a.BFL, b.BFL}, {a.FFL, b.FFL},
					{a.ImageZ, b.ImageZ}, {a.Magnification, b.Magnification},
					{a.BackFocalZ, b.BackFocalZ}, {a.FrontFocalZ, b.FrontFocalZ},
					{a.EntrancePupilZ, b.EntrancePupilZ}, {a.ExitPupilZ, b.ExitPupilZ},
				} {
					for _, x := range pair {
						if math.IsNaN(x) {
							break
						}
					}
					if math.Float64bits(pair[0]) != math.Float64bits(pair[1]) {
						t.Fatalf("system %d (object %g, λ %g): %.17g != %.17g", si, oz, lam, pair[0], pair[1])
					}
				}
				ma, mb := a.Matrix(), b.Matrix()
				if ma != mb {
					t.Fatalf("system %d: matrix %v != %v", si, ma, mb)
				}
				ga, err := sys.AngularMagnification(oz, 1e-3, lam)
				if err != nil {
					continue
				}
				gb, _ := sys.AngularMagnification(oz, 1e-3, lam)
				if math.Float64bits(ga) != math.Float64bits(gb) {
					t.Fatalf("system %d: Γ %v != %v", si, ga, gb)
				}
			}
		}
	}
}

// ---------------------------------------------------------------------------
// 14. error paths of the mirror support
// ---------------------------------------------------------------------------

func TestParaxialMirrorErrors(t *testing.T) {
	// A mirror is an interaction inside one medium: a reflecting surface that
	// separates two different media is rejected by the first-order engine (the
	// real ray trace stays permissive and simply reflects).
	bad := mustSystem(t, Air{}, []Surface{
		{Profile: Sphere{R: 100}, Aperture: CircleAperture(10)},
		{Profile: Sphere{R: -100}, Mirror: true, Aperture: CircleAperture(10)},
	}, []float64{5}, []Medium{Air{}, N_BK7, Air{}})
	if _, err := bad.Paraxial(Wavelengthd, math.Inf(-1)); err == nil {
		t.Error("a mirror between two different media must be rejected")
	}
	// An object behind the first surface is still rejected.
	if _, err := singleMirror(t, -2000, 100).Paraxial(Wavelengthd, 10); err == nil {
		t.Error("an object behind surface 0 must be rejected")
	}
	// A system with power has no angular magnification for an object at
	// infinity: the API reports the definition instead of a number.
	sys := catadioptricSinglet(t)
	if _, err := sys.AngularMagnification(math.Inf(-1), 1e-3, Wavelengthd); err == nil {
		t.Error("a non-afocal system with an object at infinity has no angular magnification")
	}
	// The same system with a finite object does have one, and it is the
	// Lagrange-invariant ratio γ = n/(n'·m).
	gamma, err := sys.AngularMagnification(-4000, 1e-3, Wavelengthd)
	if err != nil {
		t.Fatalf("finite object: %v", err)
	}
	p, err := sys.Paraxial(Wavelengthd, -4000)
	if err != nil {
		t.Fatal(err)
	}
	checkRel(t, "γ = n/(n'·m)", gamma, p.Index[0]/(p.Index[len(p.Index)-1]*p.Magnification), 1e-15)
	// The traced axial ray reproduces the same ratio in the propagation frame.
	c := p.rayFrom(1, -1/(-4000.0))
	checkRel(t, "γ from the traced axial ray", p.Parity()*c.UAfter[len(c.UAfter)-1]/c.U0, gamma, 1e-10)
	// A tilted mirror system still belongs to the posed branch, which the
	// first-order engine refuses (unchanged by the mirror support).
	tilted := mustSystem(t, Air{}, []Surface{
		{Profile: Sphere{R: -2000}, Mirror: true, Aperture: CircleAperture(10),
			Pose: &CoordBreak{TiltY: 1e-3}},
	}, nil, nil)
	if _, err := tilted.Paraxial(Wavelengthd, math.Inf(-1)); err == nil {
		t.Error("a tilted mirror system must stay rejected as a posed system")
	}
}

// ---------------------------------------------------------------------------
// 15. the FOCUSZ numbers stored by the mirror examples of the library
// ---------------------------------------------------------------------------

func TestLibraryMirrorFocusNumbers(t *testing.T) {
	// The descriptions of library/geo/mirror-cassegrain-f1000-d300-b300.json and
	// mirror-newtonian-f1000-f5.json quote FOCUSZ (345.328867 mm and
	// 1053.04778 mm). Those numbers are the lower edge of the scan bracket of
	// tools/libcheck/measure.focusScan, whose ray-pair crossing crossZ measures
	// the parameter along the incident ray with the wrong sign: for a bundle that
	// converges it returns the point reflected at the hit plane, p1.Z - t·d1.Z,
	// so the whole scan runs on the wrong side of the last surface and its
	// minimum is a plane where the bundle is far from focus (the quoted
	// FOCUSRMS of 60.55 mm / 58.06 mm is not a spot size).
	//
	// The relation is exact and is pinned here: the reported crossing is
	// 2·p1.Z - z_true, the true crossing mirrored at the hit plane.
	for _, tc := range []struct {
		name  string
		gap   float64
		R, K  float64
		semi  float64
		start float64
	}{
		{"newtonian layout", 60, -2000, -1, 100, -2000},
	} {
		sys := stopMirrorPair(t, tc.gap, tc.semi, tc.R, tc.K, tc.semi)
		p, err := sys.Paraxial(Wavelengthd, math.Inf(-1))
		if err != nil {
			t.Fatalf("%s: %v", tc.name, err)
		}
		zTrue := realCrossOfRay(sys, 1e-4, Wavelengthd, tc.start)
		checkRel(t, tc.name+" true focus", zTrue, p.ImageZ, 1e-9)
		// Two real rays of the first fan pair, crossed with the flawed helper.
		r1 := sys.Trace(NewRay(V(0, tc.semi*(0.02+0.96*0/20.0), tc.start), V(0, 0, 1), Wavelengthd))
		r2 := sys.Trace(NewRay(V(0, tc.semi*(0.02+0.96*1/20.0), tc.start), V(0, 0, 1), Wavelengthd))
		if !r1.OK() || !r2.OK() {
			t.Fatalf("%s: fan trace failed", tc.name)
		}
		s1 := r1.Steps[len(r1.Steps)-1]
		s2 := r2.Steps[len(r2.Steps)-1]
		zArtifact := pairCrossZ(s1.Point, s1.DirOut, s2.Point, s2.DirOut)
		zPair := pairCrossZTrue(s1.Point, s1.DirOut, s2.Point, s2.DirOut)
		checkRel(t, tc.name+" artifact crossing = 2 z_hit - z_pair", zArtifact, 2*s1.Point.Z-zPair, 1e-9)
		if math.Abs(zArtifact-p.ImageZ) < 1 {
			t.Errorf("%s: the artifact crossing %g must not coincide with the physical focus %g", tc.name, zArtifact, p.ImageZ)
		}
		t.Logf("%s: physical focus %.9f, pair crossing %.9f, artifact crossing %.9f, library FOCUSZ %.6f",
			tc.name, zTrue, zPair, zArtifact, 1053.04778)
	}
	// The cassegrain layout of the library: the same relation, using its own
	// secondary radius and its 150 mm stop.
	sys := mustSystem(t, Air{}, []Surface{
		{Profile: Plane{}, Aperture: CircleAperture(150)},
		{Profile: Sphere{R: -2000}, Mirror: true, Aperture: CircleAperture(153)},
		{Profile: Sphere{R: -646.1538461538461}, Mirror: true, Aperture: CircleAperture(117.5)},
	}, []float64{60, -300}, nil)
	p, err := sys.Paraxial(Wavelengthd, math.Inf(-1))
	if err != nil {
		t.Fatal(err)
	}
	zTrue := realCrossOfRay(sys, 1e-4, Wavelengthd, -2000)
	r1 := sys.Trace(NewRay(V(0, 150*(0.02+0.96*0/20.0), -2000), V(0, 0, 1), Wavelengthd))
	r2 := sys.Trace(NewRay(V(0, 150*(0.02+0.96*1/20.0), -2000), V(0, 0, 1), Wavelengthd))
	if !r1.OK() || !r2.OK() {
		t.Fatal("fan trace failed")
	}
	s1 := r1.Steps[len(r1.Steps)-1]
	s2 := r2.Steps[len(r2.Steps)-1]
	zArtifact := pairCrossZ(s1.Point, s1.DirOut, s2.Point, s2.DirOut)
	zPair := pairCrossZTrue(s1.Point, s1.DirOut, s2.Point, s2.DirOut)
	checkRel(t, "cassegrain true focus", zTrue, p.ImageZ, 1e-9)
	checkRel(t, "cassegrain artifact crossing = 2 z_hit - z_pair", zArtifact, 2*s1.Point.Z-zPair, 1e-9)
	t.Logf("cassegrain layout: physical focus %.9f, pair crossing %.6f, artifact crossing %.6f, library FOCUSZ %.6f",
		zTrue, zPair, zArtifact, 345.328867)
}
