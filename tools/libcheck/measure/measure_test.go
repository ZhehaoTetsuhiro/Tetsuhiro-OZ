package measure

import (
	"math"
	"strings"
	"testing"

	"tetsuhiro-oz/geo"
	"tetsuhiro-oz/lens"
)

// refractive and focus examples used by the clause tests
func sampleRefract() *Metrics {
	return &Metrics{Mode: ModeRefract, EFLMM: 100.123456789, BFLMM: 95.5, FNumber: 5, ColorFCMM: 0.0123,
		Field: 0.5, FieldUnit: "deg", RMSAxisMM: 0.00123456, RMSFieldMM: 0.0045, Rays: 37, FieldRays: 37}
}

// TestClauseRoundTrip: a clause written from a measurement reads back to the
// same numbers, within the library tolerance (nine significant digits are
// far inside 1e-6 relative).
func TestClauseRoundTrip(t *testing.T) {
	for _, m := range []*Metrics{
		sampleRefract(),
		{Mode: ModeFocus, FocusZMM: -12.3456789, FocusRMSMM: 0.0009, Rays: 21},
		{Mode: ModeAfocal, SlopeRMSUrad: 412.5, SlopeMeanUrad: -3.25, Mag: -2.0251, Rays: 21},
		{Mode: ModeCollimate, SlopeRMSUrad: 1.5, SlopeMeanUrad: 0.2, Rays: 21},
	} {
		desc := "示例。" + ExpectClause(m)
		mode, vals, err := ParseClause(desc)
		if err != nil {
			t.Fatalf("%s: parse: %v", m.Mode, err)
		}
		if probs := Compare(m, mode, vals); len(probs) != 0 {
			t.Errorf("%s: round trip mismatch: %v", m.Mode, probs)
		}
	}
}

// TestCompareTolerance: a 1e-3 relative error in a description is caught; a
// 1e-9 relative error is inside the stated tolerance 1e-6.
func TestCompareTolerance(t *testing.T) {
	m := sampleRefract()
	mode, vals, err := ParseClause("x。" + ExpectClause(m))
	if err != nil {
		t.Fatal(err)
	}
	vals["EFL"] *= 1 + 1e-9
	if probs := Compare(m, mode, vals); len(probs) != 0 {
		t.Errorf("1e-9 relative change rejected: %v", probs)
	}
	vals["EFL"] = m.EFLMM * (1 + 1e-3)
	if probs := Compare(m, mode, vals); len(probs) != 1 || !strings.HasPrefix(probs[0], "EFL:") {
		t.Errorf("1e-3 relative change not reported as EFL: %v", probs)
	}
}

// TestCompareReportsMissingAndExtraKeys: a clause that drops or adds a number
// is a failure, not a silent pass.
func TestCompareReportsMissingAndExtraKeys(t *testing.T) {
	m := sampleRefract()
	mode, vals, err := ParseClause("x。" + ExpectClause(m))
	if err != nil {
		t.Fatal(err)
	}
	delete(vals, "BFL")
	vals["EXTRA"] = 1
	probs := Compare(m, mode, vals)
	joined := strings.Join(probs, " | ")
	if !strings.Contains(joined, "clause misses BFL") || !strings.Contains(joined, "unexpected key EXTRA") {
		t.Errorf("missing/extra keys not reported: %v", probs)
	}
}

// TestParseClauseRequiresMode: a description without the clause is rejected.
func TestParseClauseRequiresMode(t *testing.T) {
	if _, _, err := ParseClause("没有预期指标的描述"); err == nil {
		t.Error("description without a clause was accepted")
	}
}

// TestFingerprintIdentity: names do not enter the fingerprint, parameters do,
// and a flat surface written as +Inf or as 0 is the same design.
func TestFingerprintIdentity(t *testing.T) {
	base := func(name string, r float64) *lens.Prescription {
		return &lens.Prescription{Name: name, ObjectZMM: math.Inf(-1), Stop: 0, Cards: []lens.Card{
			{Name: "stop", RadiusMM: math.Inf(1), ThicknessMM: 5, SemiAperMM: 5},
			{Name: "front", RadiusMM: r, ThicknessMM: 3, Glass: "N-BK7", SemiAperMM: 5},
			{Name: "rear", RadiusMM: -r, ThicknessMM: 10, SemiAperMM: 5},
		}}
	}
	a, err := Fingerprint(base("alpha", 50))
	if err != nil {
		t.Fatal(err)
	}
	b, _ := Fingerprint(base("beta", 50))
	if a != b {
		t.Error("renaming changed the fingerprint")
	}
	c, _ := Fingerprint(base("alpha", 50.000001))
	if a == c {
		t.Error("a 1e-6 mm radius change did not change the fingerprint")
	}
	p := base("flat", 0)
	p.Cards[0].RadiusMM = 0
	d1, _ := Fingerprint(p)
	q := base("flat", 0)
	q.Cards[0].RadiusMM = math.Inf(1)
	d2, _ := Fingerprint(q)
	if d1 != d2 {
		t.Error("flat surface written as 0 and as +Inf gives different fingerprints")
	}
}

// ---------- posed systems (decentred / tilted) ----------

// posedLens builds a plano-convex singlet with an optional decentre of its
// two surfaces; posedLens(0, 0) is the coaxial reference.
func posedLens(decenterY float64) *lens.Prescription {
	p := &lens.Prescription{ObjectZMM: math.Inf(-1), Stop: 0, Cards: []lens.Card{
		{Name: "aperture stop", RadiusMM: math.Inf(1), ThicknessMM: 2.5, SemiAperMM: 12.5},
		{Name: "L1", RadiusMM: 51.68, ThicknessMM: 5, Glass: "N-BK7", SemiAperMM: 16, DecenterYMM: decenterY},
		{Name: "rear", RadiusMM: math.Inf(1), ThicknessMM: 30, SemiAperMM: 15, DecenterYMM: decenterY},
	}}
	return p
}

// TestPosedMeasurementUsesRealRays: a decentred prescription has no
// first-order solution (geo refuses it), so Measure must classify it as
// ModePosed and locate a focus with the real tracer; the focus must stay in
// the neighbourhood of the coaxial image and the RMS must pass Check.
func TestPosedMeasurementUsesRealRays(t *testing.T) {
	ref, err := Measure(posedLens(0))
	if err != nil {
		t.Fatal(err)
	}
	if _, perr := posedLens(0).Paraxial(500); perr != nil {
		t.Fatalf("the coaxial reference must have a paraxial solution: %v", perr)
	}
	if ref.Mode != ModeRefract {
		t.Fatalf("the coaxial reference must measure as %q, got %q", ModeRefract, ref.Mode)
	}
	p := posedLens(0.5)
	if _, perr := p.Paraxial(500); perr == nil {
		t.Fatal("a decentred prescription must be refused by the paraxial engine")
	}
	m, err := Measure(p)
	if err != nil {
		t.Fatal(err)
	}
	if m.Mode != ModePosed {
		t.Fatalf("mode %q, want %q", m.Mode, ModePosed)
	}
	if m.Rays < 16 {
		t.Errorf("only %d of 21 rays reach the last surface", m.Rays)
	}
	if !(m.FocusRMSMM > 0) || m.FocusRMSMM > 0.05 {
		t.Errorf("focus RMS %gmm is not a sane real-ray focus", m.FocusRMSMM)
	}
	if m.FocusZMM < 90 || m.FocusZMM > 115 {
		t.Errorf("focus z %gmm is far from the ~100 mm image of the design", m.FocusZMM)
	}
	if probs := Check(p, m); len(probs) > 0 {
		t.Errorf("posed measurement failed the library assertions: %v", probs)
	}
}

// TestPosedScanFindsTiltedMirrorFocus: a spherical mirror of radius R tilted
// by θ about x focuses the parallel bundle at R/2 along the reflected axis
// (0, sin2θ, −cos2θ) from its vertex — the closed form of the mirror focus.
// The posed scan must find that plane to a few percent.
func TestPosedScanFindsTiltedMirrorFocus(t *testing.T) {
	const R, tiltDeg, zMirror = -400.0, 10.0, 60.0
	p := &lens.Prescription{ObjectZMM: math.Inf(-1), Stop: 0, Cards: []lens.Card{
		{Name: "aperture stop", RadiusMM: math.Inf(1), ThicknessMM: zMirror, SemiAperMM: 25},
		{Name: "primary", RadiusMM: R, Mirror: true, ThicknessMM: 200, SemiAperMM: 25.5, TiltXDeg: tiltDeg},
	}}
	m, err := Measure(p)
	if err != nil {
		t.Fatal(err)
	}
	if m.Mode != ModePosed {
		t.Fatalf("mode %q, want %q", m.Mode, ModePosed)
	}
	th := 2 * tiltDeg * math.Pi / 180
	wantZ := zMirror - math.Abs(R)/2*math.Cos(th)
	if rel := math.Abs(m.FocusZMM-wantZ) / math.Abs(wantZ); !(rel <= 0.02) {
		t.Errorf("focus z %.6g, want %.6g (relative %.3g)", m.FocusZMM, wantZ, rel)
	}
	// the axial focus of the untilted mirror is at z = R/2; for the tilted one
	// the spot is dominated by its astigmatism, but it stays small against the
	// aperture
	if !(m.FocusRMSMM > 0) || m.FocusRMSMM > 0.5 {
		t.Errorf("focus RMS %gmm is not a sane tilted-mirror focus", m.FocusRMSMM)
	}
}

// TestPosedParallelOutputHasNoFocus: a 45° plane mirror folds the parallel
// bundle but does not focus it; the posed scan must report that instead of a
// fabricated focus.
func TestPosedParallelOutputHasNoFocus(t *testing.T) {
	p := &lens.Prescription{ObjectZMM: math.Inf(-1), Stop: 0, Cards: []lens.Card{
		{Name: "aperture stop", RadiusMM: math.Inf(1), ThicknessMM: 60, SemiAperMM: 12},
		{Name: "fold", RadiusMM: 0, Mirror: true, ThicknessMM: 50, SemiAperMM: 15, TiltXDeg: 45},
	}}
	if _, err := Measure(p); err == nil {
		t.Fatal("a plane-mirror fold has no focus; Measure must say so")
	}
}

// TestFingerprintCoversPoseAndAspheres: the pose fields and the two aspheric
// descriptions are part of the design identity, and a prescription without
// them keeps the historical fingerprint byte for byte.
func TestFingerprintCoversPoseAndAspheres(t *testing.T) {
	base := func() *lens.Prescription {
		return &lens.Prescription{ObjectZMM: math.Inf(-1), Stop: 0, Cards: []lens.Card{
			{Name: "stop", RadiusMM: math.Inf(1), ThicknessMM: 5, SemiAperMM: 5},
			{Name: "front", RadiusMM: 50, ThicknessMM: 3, Glass: "N-BK7", SemiAperMM: 5},
			{Name: "rear", RadiusMM: -50, ThicknessMM: 10, SemiAperMM: 5},
		}}
	}
	ref, err := Fingerprint(base())
	if err != nil {
		t.Fatal(err)
	}
	// the golden value documents the fingerprint format of the library as it
	// was before the pose and asphere fields existed; the new fields are
	// omitempty, so a legacy prescription keeps this exact hash
	const golden = "79149bce07bfe0411a7466eb0435ec35b6e0ebce3c882b57aa15efcb7a4e8743"
	if ref != golden {
		t.Errorf("legacy fingerprint changed: %s, want %s", ref, golden)
	}
	kinds := []struct {
		name string
		mut  func(p *lens.Prescription)
	}{
		{"decentre", func(p *lens.Prescription) { p.Cards[1].DecenterYMM = 0.1 }},
		{"tilt", func(p *lens.Prescription) { p.Cards[2].TiltXDeg = 0.05 }},
		{"poly terms", func(p *lens.Prescription) {
			p.Cards[1].PolyTerms = []lens.PolyTermCard{{Power: 4, Coeff: -1e-6}}
		}},
		{"q asphere", func(p *lens.Prescription) {
			p.Cards[1].QAsphere = &lens.QCard{Form: "qcon", NormRadiusMM: 5, Coeffs: []float64{1e-3}}
		}},
	}
	for _, k := range kinds {
		q := base()
		k.mut(q)
		fp, err := Fingerprint(q)
		if err != nil {
			t.Fatalf("%s: %v", k.name, err)
		}
		if fp == ref {
			t.Errorf("%s does not enter the fingerprint", k.name)
		}
	}
}

// ---------- focus scan and the ray-pair crossing ----------

// TestCrossZIntersection: crossZ must return the intersection of the two ray
// lines. A regression: the parameter t used to carry the wrong sign, which
// mirrored every crossing at the hit plane (2·z_hit − z) and put the whole
// focusScan bracket of every focus-mode example on the wrong side of the last
// surface (the historical FOCUSZ values of the mirror library were artifacts).
func TestCrossZIntersection(t *testing.T) {
	// y = z and y = 1 − z cross at (y, z) = (0.5, 0.5).
	got := crossZ(geo.V(0, 0, 0), geo.V(0, 1, 1), geo.V(0, 1, 0), geo.V(0, -1, 1))
	if math.Abs(got-0.5) > 1e-12 {
		t.Errorf("crossZ = %g, want 0.5", got)
	}
	// The mirrored value (2·z_hit − z = −0.5) must not be returned.
	if math.Abs(got+0.5) < 1e-12 {
		t.Errorf("crossZ returned the mirrored crossing %g", got)
	}
}

// TestFocusScanParabolicMirror: a parabolic mirror of R = −2000 focuses a
// parallel bundle at R/2, i.e. 1000 mm in front of the mirror vertex — for a
// mirror 60 mm behind the stop the focus is at z = −940 with a vanishing RMS.
// Before the crossZ fix the scan bracket sat near +1053 (the crossing
// mirrored at the mirror) and the measured FOCUSZ was that artifact.
func TestFocusScanParabolicMirror(t *testing.T) {
	p := &lens.Prescription{ObjectZMM: math.Inf(-1), Stop: 0, Cards: []lens.Card{
		{Name: "aperture stop", RadiusMM: math.Inf(1), ThicknessMM: 60, SemiAperMM: 100},
		{Name: "primary", RadiusMM: -2000, ConicK: -1, Mirror: true, ThicknessMM: 200, SemiAperMM: 102},
	}}
	m, err := Measure(p)
	if err != nil {
		t.Fatal(err)
	}
	if m.Mode != ModeFocus {
		t.Fatalf("mode %q, want %q", m.Mode, ModeFocus)
	}
	if math.Abs(m.FocusZMM+940) > 0.01 {
		t.Errorf("FOCUSZ = %.9g mm, want −940 (R/2 of the parabola)", m.FocusZMM)
	}
	if !(m.FocusRMSMM < 1e-3) {
		t.Errorf("FOCUSRMS = %.9g mm, want ~0 for a perfect parabola", m.FocusRMSMM)
	}
}
