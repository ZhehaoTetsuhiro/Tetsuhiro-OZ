package lens

import (
	"encoding/json"
	"math"
	"os"
	"strings"
	"testing"

	"tetsuhiro-oz/geo"
)

func biconvexBK7() *Prescription {
	return &Prescription{
		Name:      "biconvex BK7",
		ObjectZMM: math.Inf(-1),
		Cards: []Card{
			{Name: "front", RadiusMM: 100, ThicknessMM: 6, Glass: "N-BK7", SemiAperMM: 12.7},
			{Name: "rear", RadiusMM: -100, ThicknessMM: 95, Glass: "", SemiAperMM: 12.7},
		},
		Stop: 0,
	}
}

// TestPrescriptionEFL checks the prescription → geo.System → paraxial chain
// against the closed-form thick-lens focal length
//
//	1/f = (n-1)(1/R1 - 1/R2 + (n-1)d/(n R1 R2)).
func TestPrescriptionEFL(t *testing.T) {
	p := biconvexBK7()
	par, err := p.Paraxial(geo.Wavelengthd)
	if err != nil {
		t.Fatal(err)
	}
	n := geo.N_BK7.N(geo.Wavelengthd)
	const R1, R2, d = 100.0, -100.0, 6.0
	invF := (n - 1) * (1/R1 - 1/R2 + (n-1)*d/(n*R1*R2))
	fThick := 1 / invF
	if math.Abs(par.EFL-fThick) > 0.05 {
		t.Fatalf("EFL = %.6f mm, want %.6f (thick-lens formula)", par.EFL, fThick)
	}
	// The thin-lens value is a sanity bracket (slightly longer focal length).
	if math.Abs(par.EFL-fThick) > math.Abs(fThick-n/((n-1)*0.02)) {
		t.Fatalf("EFL differs from both thin and thick formula")
	}
	if !(par.BFL > 0) || !(par.FNumber > 0) {
		t.Fatalf("degenerate first-order numbers: BFL=%g F#=%g", par.BFL, par.FNumber)
	}
}

func TestPrescriptionJSONRoundTrip(t *testing.T) {
	p := biconvexBK7()
	dir := t.TempDir()
	path := dir + "/doublet.json"
	if err := p.SaveJSON(path); err != nil {
		t.Fatal(err)
	}
	q, err := LoadJSON(path)
	if err != nil {
		t.Fatal(err)
	}
	a, err := p.Paraxial(geo.Wavelengthd)
	if err != nil {
		t.Fatal(err)
	}
	b, err := q.Paraxial(geo.Wavelengthd)
	if err != nil {
		t.Fatal(err)
	}
	if a.EFL != b.EFL || a.BFL != b.BFL {
		t.Fatalf("JSON round trip changed first order: %g/%g vs %g/%g", a.EFL, a.BFL, b.EFL, b.BFL)
	}
}

func TestPrescriptionMetadataRoundTrip(t *testing.T) {
	p := biconvexBK7()
	p.Description = "双凸单透镜，N-BK7，f≈190mm；轴上球差为主"
	p.Category = "singlet"
	p.Difficulty = 1
	p.Tags = []string{"N-BK7", "F/7.5", "visible"}

	dir := t.TempDir()
	path := dir + "/meta.json"
	if err := p.SaveJSON(path); err != nil {
		t.Fatal(err)
	}
	raw, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	for _, key := range []string{`"description"`, `"category"`, `"difficulty"`, `"tags"`} {
		if !strings.Contains(string(raw), key) {
			t.Fatalf("saved JSON is missing %s:\n%s", key, raw)
		}
	}
	q, err := LoadJSON(path)
	if err != nil {
		t.Fatal(err)
	}
	if q.Description != p.Description || q.Category != p.Category || q.Difficulty != p.Difficulty {
		t.Fatalf("metadata round trip: got %q/%q/%d, want %q/%q/%d",
			q.Description, q.Category, q.Difficulty, p.Description, p.Category, p.Difficulty)
	}
	if len(q.Tags) != len(p.Tags) {
		t.Fatalf("tags round trip: got %v, want %v", q.Tags, p.Tags)
	}
	for i := range p.Tags {
		if q.Tags[i] != p.Tags[i] {
			t.Fatalf("tags round trip: got %v, want %v", q.Tags, p.Tags)
		}
	}
	// The metadata is inert: the first order must be untouched.
	a, err := p.Paraxial(geo.Wavelengthd)
	if err != nil {
		t.Fatal(err)
	}
	b, err := q.Paraxial(geo.Wavelengthd)
	if err != nil {
		t.Fatal(err)
	}
	if a.EFL != b.EFL || a.BFL != b.BFL {
		t.Fatalf("metadata changed the first order: %g/%g vs %g/%g", a.EFL, a.BFL, b.EFL, b.BFL)
	}
	// Omitted when empty.
	plain := biconvexBK7()
	buf, err := json.Marshal(plain)
	if err != nil {
		t.Fatal(err)
	}
	for _, key := range []string{`"description"`, `"category"`, `"difficulty"`, `"tags"`} {
		if strings.Contains(string(buf), key) {
			t.Fatalf("empty metadata was not omitted (%s):\n%s", key, buf)
		}
	}
	// Clone duplicates the tag slice: editing the copy must not touch the
	// original.
	clone := p.Clone()
	clone.Tags[0] = "mutated"
	if p.Tags[0] != "N-BK7" {
		t.Fatalf("Clone shares the tags backing array: %v", p.Tags)
	}
}

func TestPrescriptionGlassErrors(t *testing.T) {
	p := biconvexBK7()
	p.Cards[0].Glass = "unobtainium"
	if err := p.Validate(); err == nil {
		t.Fatal("unknown glass accepted")
	}
	if _, err := p.System(); err == nil {
		t.Fatal("System() accepted unknown glass")
	}
}

func TestPrescriptionPlaneShorthand(t *testing.T) {
	for _, r := range []float64{0, math.Inf(1), math.Inf(-1)} {
		prof, err := (Card{RadiusMM: r}).profile()
		if err != nil {
			t.Fatal(err)
		}
		if _, ok := prof.(geo.Plane); !ok {
			t.Fatalf("RadiusMM=%v produced %T, want geo.Plane", r, prof)
		}
	}
	// Pure polynomial asphere on a flat base.
	prof, err := (Card{RadiusMM: 0, Coeffs: []float64{1e-8, 1e-12}}).profile()
	if err != nil {
		t.Fatal(err)
	}
	if _, ok := prof.(geo.EvenAsphere); !ok {
		t.Fatalf("flat asphere produced %T, want geo.EvenAsphere", prof)
	}
}

func TestPrescriptionVarsPlumbing(t *testing.T) {
	p := biconvexBK7()
	specs := []VarSpec{{0, VarRadius, 0}, {0, VarThickness, 0}, {1, VarRadius, 0}}
	x, err := p.Vars(specs)
	if err != nil {
		t.Fatal(err)
	}
	if x[0] != 100 || x[1] != 6 || x[2] != -100 {
		t.Fatalf("Vars: %v", x)
	}
	x[0], x[1], x[2] = 90, 5, -110
	if err := p.SetVars(specs, x); err != nil {
		t.Fatal(err)
	}
	if p.Cards[0].RadiusMM != 90 || p.Cards[0].ThicknessMM != 5 || p.Cards[1].RadiusMM != -110 {
		t.Fatalf("SetVars did not land: %+v", p.Cards)
	}
	if _, err := p.Vars([]VarSpec{{9, VarRadius, 0}}); err == nil {
		t.Fatal("out-of-range spec accepted")
	}
}

func TestPrescriptionMirrorRejectedByParaxial(t *testing.T) {
	p := biconvexBK7()
	p.Cards[1].Mirror = true
	if _, err := p.Paraxial(geo.Wavelengthd); err == nil {
		t.Fatal("paraxial analysis accepted a mirror system")
	}
	// The system itself must still build for real-ray tracing.
	if _, err := p.System(); err != nil {
		t.Fatalf("System() must build mirror systems: %v", err)
	}
}

// TestPrescriptionCloneIndependence makes sure Clone deep-copies Coeffs.
func TestPrescriptionCloneIndependence(t *testing.T) {
	p := biconvexBK7()
	p.Cards[0].Coeffs = []float64{1e-8}
	q := p.Clone()
	q.Cards[0].Coeffs[0] = 9e-9
	if p.Cards[0].Coeffs[0] != 1e-8 {
		t.Fatal("Clone shares the Coeffs slice")
	}
}
