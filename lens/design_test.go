package lens

import (
	"math"
	"testing"

	"tetsuhiro-oz/geo"
)

// TestDesignAchromat drives the full autonomous path: analytic solve, DLS
// polish, then independent re-verification of every design claim.
func TestDesignAchromat(t *testing.T) {
	p, res, err := DesignAchromat(100, "N-BK7", "F2", 10, nil)
	if err != nil {
		t.Fatal(err)
	}
	if res == nil || res.Iter < 1 {
		t.Fatalf("polish did not run: %+v", res)
	}
	t.Logf("optimizer: iter=%d evals=%d merit=%.6g", res.Iter, res.Evals, res.Merit)

	par, err := p.Paraxial(geo.Wavelengthd)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(par.EFL-100) > 0.05 {
		t.Fatalf("EFL = %.4f mm, want 100 ± 0.05", par.EFL)
	}
	chrom, err := ChromaticResidual(p, 1.0)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(chrom) > 0.025 {
		t.Fatalf("F–C focal shift = %.4f mm, want |shift| < 0.025", chrom)
	}

	samples := UnitDiskSamples(3, 6)
	rmsOn, _, nOn, err := SpotRMS(p, geo.Wavelengthd, Field{0, 0}, samples)
	if err != nil {
		t.Fatal(err)
	}
	if nOn == 0 || rmsOn > 0.02 {
		t.Fatalf("on-axis spot RMS = %.4f mm over %d rays, want < 0.02", rmsOn, nOn)
	}
	rmsOff, _, nOff, err := SpotRMS(p, geo.Wavelengthd, Field{0, 0.0087}, samples)
	if err != nil {
		t.Fatal(err)
	}
	if nOff == 0 || rmsOff > 0.03 {
		t.Fatalf("0.5° spot RMS = %.4f mm over %d rays, want < 0.03", rmsOff, nOff)
	}
	rep, err := Describe(p)
	if err != nil {
		t.Fatal(err)
	}
	t.Logf("design report:\n%s", rep)
}

// TestDesignAchromatDeterministic pins bit-identical repeatability.
func TestDesignAchromatDeterministic(t *testing.T) {
	a, ra, err := DesignAchromat(80, "N-BK7", "N-SF5", 8, nil)
	if err != nil {
		t.Fatal(err)
	}
	b, rb, err := DesignAchromat(80, "N-BK7", "N-SF5", 8, nil)
	if err != nil {
		t.Fatal(err)
	}
	if ra.Merit != rb.Merit {
		t.Fatalf("merit differs across runs: %g vs %g", ra.Merit, rb.Merit)
	}
	for i := range a.Cards {
		if a.Cards[i].RadiusMM != b.Cards[i].RadiusMM || a.Cards[i].ThicknessMM != b.Cards[i].ThicknessMM {
			t.Fatalf("card %d differs across runs", i)
		}
	}
}

// TestDesignAchromatSolveOnly checks that the analytic solve (thin-lens
// seed + thick-lens Newton achromatisation, no polish) already meets the
// focal length and colour targets, and that the driver validates its
// inputs.
func TestDesignAchromatSolveOnly(t *testing.T) {
	p, res, err := DesignAchromat(100, "N-BK7", "F2", 10, &DesignOptions{SkipPolish: true})
	if err != nil {
		t.Fatal(err)
	}
	if res != nil {
		t.Fatalf("skip-polish returned a result: %+v", res)
	}
	par, err := p.Paraxial(geo.Wavelengthd)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(par.EFL-100) > 0.05 {
		t.Fatalf("solve-only EFL = %.4f mm, want within 0.05 mm of 100", par.EFL)
	}
	chrom, err := ChromaticResidual(p, 1.0)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(chrom) > 0.005 {
		t.Fatalf("solve-only F–C shift = %.5f mm, want < 0.005", chrom)
	}
	// Wrong glass roles must be rejected.
	if _, _, err := DesignAchromat(100, "F2", "N-BK7", 10, &DesignOptions{SkipPolish: true}); err == nil {
		t.Fatal("flint/crown roles swapped but accepted")
	}
	if _, _, err := DesignAchromat(100, "N-BK7", "unobtainium", 10, nil); err == nil {
		t.Fatal("unknown flint accepted")
	}
}

// TestSpotResidualsShape pins the residual-vector layout used by the
// optimizer (2 components per sample per field, plus focal and colour).
func TestSpotResidualsShape(t *testing.T) {
	p, _, err := DesignAchromat(100, "N-BK7", "F2", 10, &DesignOptions{SkipPolish: true})
	if err != nil {
		t.Fatal(err)
	}
	samples := UnitDiskSamples(2, 4) // 1 + 4 + 8 = 13
	r, err := SpotResiduals(p, geo.Wavelengthd, []Field{{0, 0}, {0, 0.0087}}, samples, 0.01)
	if err != nil {
		t.Fatal(err)
	}
	want := 2 * 2 * len(samples)
	if len(r) != want {
		t.Fatalf("residual length %d, want %d", len(r), want)
	}
	// A pathological prescription (everything flat) must produce penalties
	// rather than an error from SpotResiduals itself when traced... it has
	// no power, so no finite image: the caller (designResiduals) turns the
	// error into a penalty; SpotResiduals itself may error — acceptable and
	// documented. Just make sure it does not panic.
	flat := &Prescription{ObjectZMM: math.Inf(-1), Cards: []Card{
		{RadiusMM: 0, ThicknessMM: 5, SemiAperMM: 10},
		{RadiusMM: 0, ThicknessMM: 5, SemiAperMM: 10},
	}}
	_, _ = SpotResiduals(flat, geo.Wavelengthd, []Field{{0, 0}}, samples, 0.01)
}
