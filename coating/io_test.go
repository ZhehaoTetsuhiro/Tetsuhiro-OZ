package coating

import (
	"path/filepath"
	"strings"
	"testing"
)

// TestStackJSONRoundTrip pins the persistence format: a designed stack
// survives save/load with identical spectra.
func TestStackJSONRoundTrip(t *testing.T) {
	s, _, err := BroadbandAR(Air, BK7, 420, 680, Options{})
	if err != nil {
		t.Fatal(err)
	}
	dir := t.TempDir()
	path := filepath.Join(dir, "ar.json")
	if err := s.SaveJSON(path); err != nil {
		t.Fatal(err)
	}
	q, err := LoadStackJSON(path)
	if err != nil {
		t.Fatal(err)
	}
	if len(q.Layers) != len(s.Layers) {
		t.Fatalf("layer count %d vs %d", len(q.Layers), len(s.Layers))
	}
	for _, lam := range []float64{420, 550, 680} {
		a, err := s.Evaluate(lam, 0, PolUnpolarized)
		if err != nil {
			t.Fatal(err)
		}
		b, err := q.Evaluate(lam, 0, PolUnpolarized)
		if err != nil {
			t.Fatal(err)
		}
		if a.R != b.R || a.T != b.T {
			t.Fatalf("spectra differ at %g nm: %g/%g vs %g/%g", lam, a.R, a.T, b.R, b.T)
		}
	}
	// The saved file references materials by name.
	raw, err := LoadStackJSON(path)
	if err != nil {
		t.Fatal(err)
	}
	if raw.Incident.Name != "Air" || raw.Substrate.Name != "BK7" {
		t.Fatalf("media not preserved: %s / %s", raw.Incident.Name, raw.Substrate.Name)
	}
}

func TestStackJSONErrors(t *testing.T) {
	bad := `{"incident":"Air","substrate":"Unobtainium","layers":[]}`
	var s Stack
	if err := s.UnmarshalJSON([]byte(bad)); err == nil || !strings.Contains(err.Error(), "unknown substrate") {
		t.Fatalf("unknown substrate accepted: %v", err)
	}
	bad2 := `{"incident":"Air","substrate":"BK7","layers":[{"mat":"Nope","d_nm":100}]}`
	if err := s.UnmarshalJSON([]byte(bad2)); err == nil {
		t.Fatal("unknown layer material accepted")
	}
	bad3 := `{"incident":"Air","substrate":"BK7","layers":[{"mat":"MgF2","d_nm":-5}]}`
	if err := s.UnmarshalJSON([]byte(bad3)); err == nil {
		t.Fatal("negative thickness accepted")
	}
	// Custom (non-library) materials cannot be serialised.
	custom := NewStack(Air, BK7, Layer{Material: NewMaterial("mysubstance", ConstIndex{N: 1.9}), ThicknessNm: 50})
	if _, err := custom.MarshalJSON(); err == nil {
		t.Fatal("custom material serialised instead of failing")
	}
}
