package measure

import (
	"math"
	"strings"
	"testing"

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
