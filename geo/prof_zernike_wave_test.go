// Acceptance item (e) of the Zernike surface task, in an external test package
// because it crosses the geometric -> wave bridge: a Z4 coefficient on a mirror
// must produce exactly the wavefront that the wave package calls "W020 of
// defocus" (W = W020·ρ² in the pupil), checked three independent ways:
//
//  1. the ray-traced OPD map of the Zernike mirror against the analytic
//     W020·ρ²;
//  2. the diffraction PSF of that map against the PSF of the wave package's
//     own synthetic W020·ρ² grid, and against the Maréchal value of the wave
//     package's defocus test (wave/psf_test.go);
//  3. against a *real* spherical defocus of the same amount: the identical
//     mirror without the Z4 term, with the detector shifted by the paraxial
//     Δz = 8·λ·F#²·W020 that produces W020.
package geo_test

import (
	"math"
	"testing"

	"tetsuhiro-oz/geo"
	"tetsuhiro-oz/wave"
)

// zsMirror builds a single-mirror system on the given profile.
func zsMirror(t *testing.T, prof geo.Profile) *geo.System {
	t.Helper()
	sys, err := geo.NewSystem(geo.Air{}, []geo.Surface{{
		Profile:  prof,
		Mirror:   true,
		Aperture: geo.CircleAperture(41), // pupil radius is 40
	}}, nil, nil)
	if err != nil {
		t.Fatalf("NewSystem: %v", err)
	}
	return sys
}

// TestZernikeZ4MirrorWaveBridge checks acceptance (e): the defocus of
// W020 = 0.25 λ added to a flat-based mirror surface as a Z4 coefficient is
// the same wavefront -- and gives the same diffraction peak -- as the wave
// package's defocus definition.
func TestZernikeZ4MirrorWaveBridge(t *testing.T) {
	const (
		lambdaNm = 550.0
		w020     = 0.25 // waves of defocus, W(ρ) = W020·ρ², the wave-bridge definition
		normR    = 40.0 // unit-disk normalisation radius of the Zernike part = pupil radius
		mirrorR  = -2000.0
		pupilZ   = -100.0
		np       = 64
		pad      = 4
	)
	focusZ := mirrorR / 2 // -1000: the Newton focus of the parabola
	lamMM := lambdaNm * 1e-6
	// A Z4 coefficient c4 on a mirror adds the round-trip wavefront
	// 2·c4·Z4(ρ) = 4c4·ρ² − 2c4 (a reflected wavefront picks up twice the
	// surface sag; the piston is fixed by the chief reference), so the
	// equivalent defocus of the wave bridge is W020_equiv = 4·c4. For
	// W020 = 0.25 λ the coefficient is therefore c4 = W020·λ/4.
	c4 := w020 * lamMM / 4

	pts := wave.SquareGrid(np)
	pup := geo.Pupil{Z: pupilZ, Radius: normR}
	field := geo.FieldTan(0, 0)

	z4 := zsMirror(t, geo.ZernikeSurface{R: mirrorR, K: -1, NormR: normR,
		Coeffs: []float64{0, 0, 0, c4}})
	plain := zsMirror(t, geo.Conic{R: mirrorR, K: -1})

	gZ4, err := z4.OPDGrid(field, pup, pts, focusZ, lambdaNm, geo.OPDOptions{})
	if err != nil {
		t.Fatalf("Z4 mirror OPDGrid: %v", err)
	}
	gPlain, err := plain.OPDGrid(field, pup, pts, focusZ, lambdaNm, geo.OPDOptions{})
	if err != nil {
		t.Fatalf("plain mirror OPDGrid: %v", err)
	}
	if gPlain.RMSWaves() > 1e-9 {
		t.Fatalf("the bare parabolic mirror is not stigmatic: RMS = %.3g waves", gPlain.RMSWaves())
	}

	// (1) Traced wavefront against the analytic W020·ρ², for every pupil
	// sample. The chief is at ρ = 0 where the map is zero by construction.
	worstRel, worstRho := 0.0, 0.0
	for i, p := range gZ4.Pupils {
		rho := math.Hypot(p.PX, p.PY)
		if rho > 1 {
			continue // outside the unit disc: blocked in the wave bridge
		}
		got := gZ4.Waves(gZ4.Raw[i])
		want := w020 * rho * rho
		if d := math.Abs(got - want); d > worstRel {
			worstRel, worstRho = d, rho
		}
	}
	if worstRel > 5e-3*w020 {
		t.Errorf("traced wavefront deviates from W020·ρ² by %.3e waves (%.2f%% of W020, at ρ = %.3f)",
			worstRel, 100*worstRel/w020, worstRho)
	}
	t.Logf("(1) traced OPD vs W020·ρ²: worst |dev| = %.3e waves = %.3f%% of W020 (at ρ = %.3f); traced σ = %.5f λ, target %.5f λ",
		worstRel, 100*worstRel/w020, worstRho, gZ4.RMSWaves(), w020/math.Sqrt(12))

	// (2) Diffraction PSF through the wave bridge. The traced map must give the
	// same peak as the synthetic W020·ρ² grid of the wave package and as the
	// Maréchal prediction exp(−(2πσ)²), σ = W020/√12 (the reference of
	// wave/psf_test.go TestPSFDefocusMarechal).
	psfZ4, err := wave.PSFFromGrid(gZ4, -focusZ, np, pad)
	if err != nil {
		t.Fatalf("PSF of the Zernike mirror: %v", err)
	}
	synthetic := &geo.OPDGrid{
		Lambda: lambdaNm,
		Pupil:  pup,
		Pupils: pts,
		Raw:    make([]float64, len(pts)),
		OPD:    make([]float64, len(pts)),
	}
	for i, p := range pts {
		synthetic.Raw[i] = w020 * (p.PX*p.PX + p.PY*p.PY) * lamMM
		synthetic.OPD[i] = synthetic.Raw[i]
	}
	psfSynth, err := wave.PSFFromGrid(synthetic, -focusZ, np, pad)
	if err != nil {
		t.Fatalf("PSF of the synthetic defocus: %v", err)
	}
	sigma := w020 / math.Sqrt(12)
	marechal := math.Exp(-(2 * math.Pi * sigma) * (2 * math.Pi * sigma))
	if d := math.Abs(psfZ4.Strehl - psfSynth.Strehl); d > 0.01 {
		t.Errorf("Strehl of the Z4 mirror %.5f differs from the wave-bridge defocus grid %.5f by %.4f",
			psfZ4.Strehl, psfSynth.Strehl, d)
	}
	if d := math.Abs(psfZ4.Strehl - marechal); d > 0.015 {
		t.Errorf("Strehl of the Z4 mirror %.5f differs from Maréchal %.5f by %.4f",
			psfZ4.Strehl, marechal, d)
	}
	// The peak must be at the focus (chief ray), i.e. at the grid centre.
	if _, r, c := psfZ4.Peak(); r != np*pad/2 || c != np*pad/2 {
		t.Errorf("PSF peak at (%d, %d), want the centre (%d, %d)", r, c, np*pad/2, np*pad/2)
	}
	t.Logf("(2) Strehl: Z4 mirror %.5f, synthetic W020·ρ² %.5f, Maréchal %.5f",
		psfZ4.Strehl, psfSynth.Strehl, marechal)

	// (3) The same defocus as a *real* spherical wavefront change: shift the
	// detector of the bare parabola by the paraxial Δz = 8·λ·F#²·W020.
	fNumber := (-mirrorR / 2) / (2 * normR)
	dz := 8 * lamMM * fNumber * fNumber * w020
	gShift, err := plain.OPDGrid(field, pup, pts, focusZ+dz, lambdaNm, geo.OPDOptions{})
	if err != nil {
		t.Fatalf("shifted-focus OPDGrid: %v", err)
	}
	psfShift, err := wave.PSFFromGrid(gShift, -focusZ, np, pad)
	if err != nil {
		t.Fatalf("PSF of the shifted focus: %v", err)
	}
	if d := math.Abs(psfZ4.Strehl - psfShift.Strehl); d > 0.02 {
		t.Errorf("Strehl of the Z4 mirror %.5f differs from the equivalent spherical defocus %.5f by %.4f",
			psfZ4.Strehl, psfShift.Strehl, d)
	}
	worstShift := 0.0
	for i, p := range gShift.Pupils {
		rho := math.Hypot(p.PX, p.PY)
		if rho > 1 {
			continue
		}
		worstShift = math.Max(worstShift, math.Abs(gShift.Waves(gShift.Raw[i])-w020*rho*rho))
	}
	t.Logf("(3) detector shift Δz = %.5f mm (F# = %.2f): Strehl %.5f, wavefront vs W020·ρ² worst |dev| = %.3e waves",
		dz, fNumber, psfShift.Strehl, worstShift)
	if worstShift > 5e-3*w020 {
		t.Errorf("the equivalent spherical defocus differs from W020·ρ² by %.3e waves", worstShift)
	}

	// (4) The literal "plane mirror" of the task: the same Z4 coefficient on a
	// *flat-based* mirror, measured differentially against the identical bare
	// plane mirror (the huge reference-sphere defocus of a collimated beam
	// cancels in the difference). The reflected wavefront doubles the surface
	// sag, so the difference must be 2c4(Z4(ρ) - Z4(0)) = W020·λ·ρ².
	planeZ4 := zsMirror(t, geo.ZernikeSurface{R: math.Inf(1), K: 0, NormR: normR,
		Coeffs: []float64{0, 0, 0, c4}})
	planeFlat := zsMirror(t, geo.Plane{})
	gPZ4, err := planeZ4.OPDGrid(field, pup, pts, -1000, lambdaNm, geo.OPDOptions{})
	if err != nil {
		t.Fatalf("flat Z4 mirror OPDGrid: %v", err)
	}
	gPF, err := planeFlat.OPDGrid(field, pup, pts, -1000, lambdaNm, geo.OPDOptions{})
	if err != nil {
		t.Fatalf("flat mirror OPDGrid: %v", err)
	}
	worstFlat := 0.0
	for i, p := range gPZ4.Pupils {
		rho := math.Hypot(p.PX, p.PY)
		if rho > 1 {
			continue
		}
		diff := gPZ4.Waves(gPZ4.Raw[i]) - gPF.Waves(gPF.Raw[i])
		worstFlat = math.Max(worstFlat, math.Abs(diff-w020*rho*rho))
	}
	if worstFlat > 5e-3*w020 {
		t.Errorf("flat-mirror Z4 wavefront differs from W020·ρ² by %.3e waves", worstFlat)
	}
	t.Logf("(4) flat mirror + Z4 (differential): worst |dev| from W020·ρ² = %.3e waves = %.3f%% of W020",
		worstFlat, 100*worstFlat/w020)
}

// TestZernikeZ11MirrorSphericalAberration is the higher-order companion of the
// Z4 test: a primary spherical term Z11 = 6ρ⁴-6ρ²+1 on the same parabolic
// mirror must appear in the wavefront at the Newton focus as 2·c11·Z11(ρ),
// i.e. the surface sag is doubled by the reflection and no spurious lower
// order appears. This also exercises the full sequential trace with a
// non-quadratic Zernike term (the rotationally symmetric r-only path).
func TestZernikeZ11MirrorSphericalAberration(t *testing.T) {
	const (
		lambdaNm = 550.0
		normR    = 40.0
		mirrorR  = -2000.0
		pupilZ   = -100.0
		np       = 64
	)
	focusZ := mirrorR / 2
	lamMM := lambdaNm * 1e-6
	c11 := 0.05 * lamMM // 0.05 λ of primary spherical aberration
	mirror := zsMirror(t, geo.ZernikeSurface{R: mirrorR, K: -1, NormR: normR,
		Coeffs: []float64{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, c11}})
	pup := geo.Pupil{Z: pupilZ, Radius: normR}
	pts := wave.SquareGrid(np)
	g, err := mirror.OPDGrid(geo.FieldTan(0, 0), pup, pts, focusZ, lambdaNm, geo.OPDOptions{})
	if err != nil {
		t.Fatalf("Z11 mirror OPDGrid: %v", err)
	}
	worst, worstRho := 0.0, 0.0
	for i, p := range g.Pupils {
		rho := math.Hypot(p.PX, p.PY)
		if rho > 1 {
			continue
		}
		// Chief-referenced map: 2c11(Z11(ρ) - Z11(0)) *is not* the form here,
		// the chief sits at ρ = 0 where Z11(0) = 1, so the raw map must be
		// 2c11(Z11(ρ) - 1) = 2c11(6ρ⁴ - 6ρ²) in waves.
		want := 2 * c11 * (6*rho*rho*rho*rho - 6*rho*rho) / lamMM
		got := g.Waves(g.Raw[i])
		if d := math.Abs(got - want); d > worst {
			worst, worstRho = d, rho
		}
	}
	peak := 2 * c11 * 6 / lamMM // |2c11·Z11| at the rim in waves
	if worst > 0.01*peak {
		t.Errorf("traced Z11 wavefront deviates from 2c11(Z11(ρ)-1) by %.3e waves (%.2f%% of the %.3f λ rim value, at ρ = %.3f)",
			worst, 100*worst/peak, peak, worstRho)
	}
	t.Logf("Z11 on the parabola: traced map vs 2c11(Z11(ρ)-1), worst |dev| = %.3e waves = %.2f%% of the %.4f λ rim value",
		worst, 100*worst/peak, peak)
}
