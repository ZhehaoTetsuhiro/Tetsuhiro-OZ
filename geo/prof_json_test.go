package geo_test

import (
	"encoding/json"
	"math"
	"strings"
	"testing"

	"tetsuhiro-oz/geo"
	"tetsuhiro-oz/lens"
)

// asphereCardPrescription returns a prescription that uses the three new card
// descriptions (polynomial with odd powers, Qcon, Qbfs) next to a plain conic.
func asphereCardPrescription() *lens.Prescription {
	return &lens.Prescription{
		Name:      "asphere-json",
		ObjectZMM: math.Inf(-1),
		Cards: []lens.Card{
			{
				Name: "poly", RadiusMM: 40, ConicK: -0.3, ThicknessMM: 4,
				Glass: "N-BK7", SemiAperMM: 10,
				PolyTerms: []lens.PolyTermCard{
					{Power: 2, Coeff: 1e-4}, {Power: 3, Coeff: 2.0 / 3 * 1e-5},
					{Power: 5, Coeff: -3e-7}, {Power: 8, Coeff: 1.0 / 7 * 1e-11},
				},
			},
			{
				Name: "qcon", RadiusMM: -60, ConicK: 0.2, ThicknessMM: 3, SemiAperMM: 9,
				QAsphere: &lens.QCard{Form: "qcon", NormRadiusMM: 9.0 / 3, Coeffs: []float64{5e-5, -1e-6, 2e-8}},
			},
			{
				Name: "qbfs", RadiusMM: 80, ThicknessMM: 0, SemiAperMM: 9,
				QAsphere: &lens.QCard{Form: "qbfs", NormRadiusMM: 9, Coeffs: []float64{2e-4, -7.5e-6, 1.25e-7}},
			},
		},
	}
}

// TestAsphereCardJSONRoundTrip is acceptance (f): a prescription with the new
// asphere descriptions survives a JSON round trip with every coefficient
// bit-identical, re-encodes to the same bytes, and builds geo profiles whose
// sags are bit-identical to those of the original.
func TestAsphereCardJSONRoundTrip(t *testing.T) {
	orig := asphereCardPrescription()
	buf, err := json.MarshalIndent(orig, "", "  ")
	if err != nil {
		t.Fatalf("marshal: %v", err)
	}
	var back lens.Prescription
	if err := json.Unmarshal(buf, &back); err != nil {
		t.Fatalf("unmarshal: %v", err)
	}
	bitsEq := func(name string, a, b float64) {
		t.Helper()
		if math.Float64bits(a) != math.Float64bits(b) {
			t.Errorf("%s: %.17g became %.17g after JSON round trip", name, a, b)
		}
	}
	// Poly terms.
	pa, pb := orig.Cards[0].PolyTerms, back.Cards[0].PolyTerms
	if len(pa) != len(pb) {
		t.Fatalf("poly terms: %d became %d", len(pa), len(pb))
	}
	for i := range pa {
		if pa[i].Power != pb[i].Power {
			t.Errorf("poly term %d power %d became %d", i, pa[i].Power, pb[i].Power)
		}
		bitsEq("poly coeff", pa[i].Coeff, pb[i].Coeff)
	}
	// Q-type cards.
	for _, idx := range []int{1, 2} {
		qa, qb := orig.Cards[idx].QAsphere, back.Cards[idx].QAsphere
		if qa == nil || qb == nil {
			t.Fatalf("card %d: q_asphere lost in round trip", idx)
		}
		if qa.Form != qb.Form {
			t.Errorf("card %d form %q became %q", idx, qa.Form, qb.Form)
		}
		bitsEq("norm radius", qa.NormRadiusMM, qb.NormRadiusMM)
		if len(qa.Coeffs) != len(qb.Coeffs) {
			t.Fatalf("card %d: %d coeffs became %d", idx, len(qa.Coeffs), len(qb.Coeffs))
		}
		for i := range qa.Coeffs {
			bitsEq("q coeff", qa.Coeffs[i], qb.Coeffs[i])
		}
	}
	// Second encoding is byte-identical.
	buf2, err := json.MarshalIndent(&back, "", "  ")
	if err != nil {
		t.Fatalf("re-marshal: %v", err)
	}
	if string(buf) != string(buf2) {
		t.Errorf("re-encoded JSON differs from the first encoding")
	}
	// The geo profiles built from both prescriptions are the same surfaces.
	sa, err := orig.System()
	if err != nil {
		t.Fatalf("original System(): %v", err)
	}
	sb, err := back.System()
	if err != nil {
		t.Fatalf("round-trip System(): %v", err)
	}
	for i := range sa.Surfaces {
		for _, r := range []float64{0.05, 1.3, 4.7, 8.2} {
			za := sa.Surfaces[i].Profile.Sag(r)
			zb := sb.Surfaces[i].Profile.Sag(r)
			if math.Float64bits(za) != math.Float64bits(zb) {
				t.Errorf("surface %d sag at r=%g: %.17g vs %.17g", i, r, za, zb)
			}
		}
	}
}

// TestAsphereCardProfileTypes checks which geo profile each card description
// produces and that the parameters are passed through.
func TestAsphereCardProfileTypes(t *testing.T) {
	sys, err := asphereCardPrescription().System()
	if err != nil {
		t.Fatalf("System(): %v", err)
	}
	poly, ok := sys.Surfaces[0].Profile.(geo.PolyAsphere)
	if !ok {
		t.Fatalf("poly card built %T, want geo.PolyAsphere", sys.Surfaces[0].Profile)
	}
	if poly.R != 40 || poly.K != -0.3 || len(poly.Terms) != 4 || poly.Terms[2] != (geo.PolyTerm{Power: 5, Coeff: -3e-7}) {
		t.Errorf("poly asphere parameters not passed through: %+v", poly)
	}
	qc, ok := sys.Surfaces[1].Profile.(geo.QAsphere)
	if !ok || qc.Form != geo.QCon || qc.NormR != 3 || qc.K != 0.2 || len(qc.Coeffs) != 3 {
		t.Errorf("qcon card built %#v, want geo.QAsphere{QCon, NormR 3, K 0.2, 3 coeffs}", sys.Surfaces[1].Profile)
	}
	qb, ok := sys.Surfaces[2].Profile.(geo.QAsphere)
	if !ok || qb.Form != geo.QBfs || qb.NormR != 9 || qb.R != 80 || len(qb.Coeffs) != 3 {
		t.Errorf("qbfs card built %#v, want geo.QAsphere{QBfs, R 80, NormR 9, 3 coeffs}", sys.Surfaces[2].Profile)
	}
}

// TestAsphereCardLegacyJSONUnchanged is the backward-compatibility test: the
// wire form of legacy prescriptions (no poly_terms / q_asphere) is byte-identical
// to the form written before these descriptions existed, and legacy documents
// still load to the same cards.
func TestAsphereCardLegacyJSONUnchanged(t *testing.T) {
	const wantEncoded = `{
  "name": "legacy-golden",
  "object_z_mm": "inf",
  "cards": [
    {
      "name": "s1",
      "radius_mm": 50,
      "conic_k": -0.5,
      "coeffs": [
        0.0001,
        -0.000002
      ],
      "thickness_mm": 3,
      "glass": "N-BK7",
      "semi_aperture_mm": 12
    },
    {
      "name": "s2",
      "radius_mm": -50,
      "thickness_mm": 0,
      "mirror": true
    }
  ],
  "stop": 1
}`
	p := lens.Prescription{
		Name:      "legacy-golden",
		ObjectZMM: math.Inf(-1),
		Cards: []lens.Card{
			{Name: "s1", RadiusMM: 50, ConicK: -0.5, Coeffs: []float64{1e-4, -2e-6}, ThicknessMM: 3, Glass: "N-BK7", SemiAperMM: 12},
			{Name: "s2", RadiusMM: -50, ThicknessMM: 0, Mirror: true},
		},
		Stop: 1,
	}
	got, err := json.MarshalIndent(&p, "", "  ")
	if err != nil {
		t.Fatalf("marshal: %v", err)
	}
	if string(got) != wantEncoded {
		t.Errorf("legacy encoding changed:\n%s\nwant:\n%s", got, wantEncoded)
	}
	const legacy = `{"name":"old","object_z_mm":"inf","cards":[{"radius_mm":30,"conic_k":-1,"thickness_mm":5,"glass":"N-SF5","semi_aperture_mm":9},{"radius_mm":-30,"thickness_mm":0}],"stop":0}`
	var q lens.Prescription
	if err := json.Unmarshal([]byte(legacy), &q); err != nil {
		t.Fatalf("legacy document rejected: %v", err)
	}
	if len(q.Cards) != 2 || q.Cards[0].RadiusMM != 30 || q.Cards[0].ConicK != -1 || q.Cards[0].PolyTerms != nil || q.Cards[0].QAsphere != nil {
		t.Errorf("legacy cards loaded wrongly: %+v", q.Cards)
	}
	sys, err := q.System()
	if err != nil {
		t.Fatalf("legacy System(): %v", err)
	}
	if _, ok := sys.Surfaces[0].Profile.(geo.Conic); !ok {
		t.Errorf("legacy conic card built %T, want geo.Conic", sys.Surfaces[0].Profile)
	}
	for _, key := range []string{"poly_terms", "q_asphere"} {
		if strings.Contains(string(got), key) {
			t.Errorf("legacy encoding contains new key %q", key)
		}
	}
}

// TestAsphereCardValidation checks the error paths of the new descriptions and
// the mutual exclusivity of the aspheric fields.
func TestAsphereCardValidation(t *testing.T) {
	mk := func(c lens.Card) *lens.Prescription {
		return &lens.Prescription{Cards: []lens.Card{c, {RadiusMM: -30, ThicknessMM: 0}}}
	}
	bad := map[string]lens.Card{
		"coeffs and poly_terms": {RadiusMM: 40, Coeffs: []float64{1e-5},
			PolyTerms: []lens.PolyTermCard{{Power: 4, Coeff: 1e-5}}},
		"poly and q_asphere": {RadiusMM: 40,
			PolyTerms: []lens.PolyTermCard{{Power: 4, Coeff: 1e-5}},
			QAsphere:  &lens.QCard{Form: "qcon", NormRadiusMM: 5, Coeffs: []float64{1e-5}}},
		"unknown q form": {RadiusMM: 40,
			QAsphere: &lens.QCard{Form: "qx", NormRadiusMM: 5, Coeffs: []float64{1e-5}}},
		"missing q form": {RadiusMM: 40,
			QAsphere: &lens.QCard{NormRadiusMM: 5, Coeffs: []float64{1e-5}}},
		"qbfs with conic K": {RadiusMM: 40, ConicK: -0.5,
			QAsphere: &lens.QCard{Form: "qbfs", NormRadiusMM: 5, Coeffs: []float64{1e-5}}},
		"q coefficients without norm radius": {RadiusMM: 40,
			QAsphere: &lens.QCard{Form: "qcon", NormRadiusMM: 0, Coeffs: []float64{1e-5}}},
		"poly power zero": {RadiusMM: 40,
			PolyTerms: []lens.PolyTermCard{{Power: 0, Coeff: 1e-5}}},
	}
	for name, c := range bad {
		p := mk(c)
		if err := p.Validate(); err == nil {
			t.Errorf("%s: Validate accepted an invalid card", name)
		}
		if _, err := p.System(); err == nil {
			t.Errorf("%s: System() accepted an invalid card", name)
		}
	}
	// Zero Q coefficients need no normalisation radius: the bare base is valid.
	ok := mk(lens.Card{RadiusMM: 40, QAsphere: &lens.QCard{Form: "qcon", Coeffs: []float64{0, 0}}})
	if err := ok.Validate(); err != nil {
		t.Errorf("zero Q coefficients rejected: %v", err)
	}
}

// TestAsphereCardCloneIsDeep checks that Clone duplicates the new slices and the
// QAsphere pointer, so editing a clone never changes the original.
func TestAsphereCardCloneIsDeep(t *testing.T) {
	orig := asphereCardPrescription()
	cl := orig.Clone()
	cl.Cards[0].PolyTerms[0].Coeff = 42
	cl.Cards[1].QAsphere.Coeffs[0] = 42
	cl.Cards[1].QAsphere.Form = "qbfs"
	if orig.Cards[0].PolyTerms[0].Coeff == 42 {
		t.Errorf("Clone shares the PolyTerms slice")
	}
	if orig.Cards[1].QAsphere.Coeffs[0] == 42 {
		t.Errorf("Clone shares the QAsphere coefficient slice")
	}
	if orig.Cards[1].QAsphere.Form != "qcon" {
		t.Errorf("Clone shares the QAsphere struct")
	}
}
