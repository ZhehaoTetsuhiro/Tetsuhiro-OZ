package optics

import (
	"math"
	"strings"
	"testing"
)

// TestExtendedSourcePresetRuns pins the built-in preset: an incoherent
// extended emitter, a rotatable slit and a screen must simulate, record one
// plane, and keep the slit's Roll at 0.
func TestExtendedSourcePresetRuns(t *testing.T) {
	var ex *Example
	for i := range Examples() {
		if strings.Contains(Examples()[i].Name, "扩展光源") {
			ex = &Examples()[i]
			break
		}
	}
	if ex == nil {
		t.Fatal("the catalog has no 扩展光源 preset")
	}
	cfg := ex.Config
	if !cfg.IsScene() {
		t.Fatalf("preset %q is not a positioned scene", ex.Name)
	}
	if cfg.Sources[0].Type != "extended" {
		t.Errorf("preset source type = %q, want extended", cfg.Sources[0].Type)
	}
	res, err := Simulate(cfg)
	if err != nil {
		t.Fatal(err)
	}
	if len(res.Planes) != 1 {
		t.Fatalf("preset recorded %d planes, want 1: %v", len(res.Planes), res.Warnings)
	}
	if !(res.Planes[0].Stats.Power > 0) {
		t.Error("the screen recorded no power")
	}
	// 6x3 = 18 point sources -> 18 incoherent units, capped at maxRetainedParts.
	if n := len(res.Planes[0].Parts); n != 18 {
		t.Errorf("preset produced %d parts, want 18 (6x3 samples)", n)
	}
	// The slit is unrotated in the preset, so the far-field streak is narrow
	// in x and wide in y.
	if res.Planes[0].Stats.RMSX >= res.Planes[0].Stats.RMSY {
		t.Errorf("unrotated slit: rms_x %.4g should be below rms_y %.4g",
			res.Planes[0].Stats.RMSX, res.Planes[0].Stats.RMSY)
	}
}

// TestExtendedSourcePresetRollSweep is the physics check on the preset: rotating
// the slit turns the far-field streak, which rms_x/rms_y measures directly.
func TestExtendedSourcePresetRollSweep(t *testing.T) {
	var cfg Config
	for _, ex := range Examples() {
		if strings.Contains(ex.Name, "扩展光源") {
			cfg = ex.Config
			break
		}
	}
	rep, err := RunScan(cfg, ScanSpec{
		Param:  "scene.components[0].roll",
		Values: []float64{0, math.Pi / 2},
		Metric: "rms_y",
	})
	if err != nil {
		t.Fatal(err)
	}
	if len(rep.Points) != 2 {
		t.Fatalf("got %d points, want 2", len(rep.Points))
	}
	for i, p := range rep.Points {
		if p.Error != "" {
			t.Fatalf("point %d: %s", i, p.Error)
		}
	}
	if !(rep.Points[1].Metric < rep.Points[0].Metric) {
		t.Errorf("rotating the slit by 90° must shrink rms_y: %.6g -> %.6g",
			rep.Points[0].Metric, rep.Points[1].Metric)
	}
}
