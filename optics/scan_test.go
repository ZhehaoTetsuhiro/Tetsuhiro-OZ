package optics

import (
	"math"
	"testing"
)

// scanCfg is the 1016-style layout: an incoherent emitter, one rotatable
// slit, and a screen.
func scanCfg() Config {
	w := 0.03
	return sceneCfg(128, w, []SourceSpec{{
		ID: "src", Label: "光源", Type: "extended", Pos: &Vec3{Z: -0.30}, Dir: &Vec3{Z: 1},
		Params: map[string]any{"power": 1e-3, "samples_x": 1, "samples_y": 1, "standoff": 0.14},
	}},
		ComponentSpec{Type: "aperture", Label: "单缝", Pos: v3(0, 0, -0.15), Shape: RectOutline(4e-4, 0.02)},
		bigDet("屏", v3(0, 0, 0.15), math.Pi, w))
}

// TestScanSubstitutesValue pins the core contract: the requested parameter is
// actually substituted, one run per value, and the metric is reported in the
// same order.
func TestScanSubstitutesValue(t *testing.T) {
	rep, err := RunScan(scanCfg(), ScanSpec{
		Param:  "scene.components[0].roll",
		Values: []float64{0, 0.5, 1.0},
		Metric: "power",
	})
	if err != nil {
		t.Fatal(err)
	}
	if len(rep.Points) != 3 {
		t.Fatalf("got %d points, want 3", len(rep.Points))
	}
	for i, want := range []float64{0, 0.5, 1.0} {
		if rep.Points[i].Value != want {
			t.Errorf("point %d: value %g, want %g", i, rep.Points[i].Value, want)
		}
		if rep.Points[i].Error != "" {
			t.Errorf("point %d: unexpected error %q", i, rep.Points[i].Error)
		}
		if !(rep.Points[i].Metric > 0) {
			t.Errorf("point %d: metric %g, want > 0", i, rep.Points[i].Metric)
		}
	}
}

// TestScanActuallyChangesTheConfig is the guard against a silent no-op: the
// reported metric must depend on the swept parameter. Rotating the slit moves
// the far-field streak from horizontal to vertical, so rms_y must drop.
func TestScanActuallyChangesTheConfig(t *testing.T) {
	rep, err := RunScan(scanCfg(), ScanSpec{
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
	if rep.Points[0].Metric <= rep.Points[1].Metric {
		t.Errorf("rotating the slit by 90° must shrink rms_y: %.6g -> %.6g",
			rep.Points[0].Metric, rep.Points[1].Metric)
	}
}

// TestScanSweepsLensFocal: a lens parameter sweep must move the focus, which
// the centroid/rms of the focal plane shows.
func TestScanSweepsLensFocal(t *testing.T) {
	w := 0.02
	cfg := sceneCfg(128, w, []SourceSpec{{
		Type: "plane", Pos: &Vec3{Z: -0.1}, Dir: &Vec3{Z: 1},
		Params: map[string]any{"power": 1e-3},
	}},
		ComponentSpec{Type: "lens", Label: "透镜", Pos: v3(0, 0, 0), Params: map[string]any{"f": 0.1},
			Shape: CircleOutline(2e-3)},
		bigDet("屏", v3(0, 0, 0.1), math.Pi, w))
	rep, err := RunScan(cfg, ScanSpec{
		Param:  "scene.components[0].params.f",
		Values: []float64{0.05, 0.1, 0.2},
		Metric: "rms_x",
	})
	if err != nil {
		t.Fatal(err)
	}
	if len(rep.Points) != 3 {
		t.Fatalf("got %d points, want 3", len(rep.Points))
	}
	// At f = 0.1 the detector sits exactly at the focus (smallest spot); the
	// 0.05 m lens focuses in front of it and the 0.2 m lens behind it, so both
	// give a larger spot than the matched one.
	if !(rep.Points[1].Metric < rep.Points[0].Metric && rep.Points[1].Metric < rep.Points[2].Metric) {
		t.Errorf("f=0.1 should be the tightest focus: %v", rep.Points)
	}
}

// TestScanSupportsEveryDocumentedPath covers the path grammar the API
// advertises, including source parameters and shape parameters.
func TestScanSupportsEveryDocumentedPath(t *testing.T) {
	cfg := scanCfg()
	cfg.Sources[0].Params["tilt_x"] = 0
	for _, path := range []string{
		"wavelength",
		"sources[0].params.tilt_x",
		"scene.components[0].roll",
		"scene.components[0].yaw",
		"scene.components[0].pitch",
		"scene.components[0].pos.z",
		"scene.components[0].shape.params.width",
		"scene.components[0].shape.params.rotation",
	} {
		rep, err := RunScan(cfg, ScanSpec{Param: path, Values: []float64{1e-3}, Metric: "power"})
		if err != nil {
			t.Errorf("path %q: %v", path, err)
			continue
		}
		if len(rep.Points) != 1 || rep.Points[0].Error != "" {
			t.Errorf("path %q: %+v", path, rep.Points)
		}
	}
}

// TestScanRejectsBadRequests covers the guard rails: a missing param, an
// empty value list, an oversized sweep, an unknown metric, an out-of-range
// index and an unsupported path must all be refused (or reported per point)
// instead of producing a plausible-looking report.
func TestScanRejectsBadRequests(t *testing.T) {
	cfg := scanCfg()
	if _, err := RunScan(cfg, ScanSpec{Values: []float64{1}}); err == nil {
		t.Error("a scan without a param must be refused")
	}
	if _, err := RunScan(cfg, ScanSpec{Param: "scene.components[0].roll"}); err == nil {
		t.Error("a scan without values must be refused")
	}
	vals := make([]float64, maxScanPoints+1)
	if _, err := RunScan(cfg, ScanSpec{Param: "scene.components[0].roll", Values: vals}); err == nil {
		t.Error("an oversized sweep must be refused")
	}
	if _, err := RunScan(cfg, ScanSpec{Param: "scene.components[0].roll", Values: []float64{0}, Metric: "nope"}); err == nil {
		t.Error("an unknown metric must be refused")
	}
	// An out-of-range index is reported per point (the run continues), not
	// silently ignored.
	rep, err := RunScan(cfg, ScanSpec{Param: "scene.components[9].roll", Values: []float64{0}})
	if err != nil {
		t.Fatal(err)
	}
	if len(rep.Points) != 1 || rep.Points[0].Error == "" {
		t.Errorf("out-of-range index: %+v", rep.Points)
	}
	// A path that is not sweepable at all (or an unknown field) is reported
	// per point, not silently ignored.
	for _, p := range []string{"grid.size", "scene.components[0].nope", "bogus"} {
		rep, err := RunScan(cfg, ScanSpec{Param: p, Values: []float64{64}, Metric: "power"})
		if err != nil {
			t.Fatalf("path %q: %v", p, err)
		}
		if len(rep.Points) != 1 || rep.Points[0].Error == "" {
			t.Errorf("path %q must be refused, got %+v", p, rep.Points)
		}
	}
}

// TestScanReportsInvalidValue: a value that makes the configuration invalid
// is reported as an error for that point, and the sweep continues.
func TestScanReportsInvalidValue(t *testing.T) {
	rep, err := RunScan(scanCfg(), ScanSpec{
		Param:  "scene.components[0].shape.params.width",
		Values: []float64{4e-4, -1e-4, 8e-4},
		Metric: "power",
	})
	if err != nil {
		t.Fatal(err)
	}
	if len(rep.Points) != 3 {
		t.Fatalf("got %d points, want 3", len(rep.Points))
	}
	if rep.Points[0].Error != "" || rep.Points[2].Error != "" {
		t.Errorf("valid values must succeed: %+v", rep.Points)
	}
	if rep.Points[1].Error == "" {
		t.Error("a negative slit width must be reported as an error for that point")
	}
}

// TestScanIsDeterministic: identical requests give identical reports.
func TestScanIsDeterministic(t *testing.T) {
	spec := ScanSpec{Param: "scene.components[0].roll", Values: []float64{0, 0.3, 0.6}, Metric: "power"}
	a, err := RunScan(scanCfg(), spec)
	if err != nil {
		t.Fatal(err)
	}
	b, err := RunScan(scanCfg(), spec)
	if err != nil {
		t.Fatal(err)
	}
	if len(a.Points) != len(b.Points) {
		t.Fatalf("reports differ in length: %d vs %d", len(a.Points), len(b.Points))
	}
	for i := range a.Points {
		if a.Points[i].Value != b.Points[i].Value || a.Points[i].Metric != b.Points[i].Metric ||
			a.Points[i].Error != b.Points[i].Error {
			t.Errorf("point %d differs: %+v vs %+v", i, a.Points[i], b.Points[i])
		}
	}
}

// TestScanDoesNotMutateTheInput: the caller's configuration must come back
// unchanged (the sweep works on copies).
func TestScanDoesNotMutateTheInput(t *testing.T) {
	cfg := scanCfg()
	before := cfg.Scene.Components[0].Roll
	beforeW := cfg.Scene.Components[0].Shape.Params["width"]
	if _, err := RunScan(cfg, ScanSpec{
		Param:  "scene.components[0].roll",
		Values: []float64{1.5},
		Metric: "power",
	}); err != nil {
		t.Fatal(err)
	}
	if _, err := RunScan(cfg, ScanSpec{
		Param:  "scene.components[0].shape.params.width",
		Values: []float64{9e-4},
		Metric: "power",
	}); err != nil {
		t.Fatal(err)
	}
	if cfg.Scene.Components[0].Roll != before {
		t.Errorf("scan mutated roll: %g -> %g", before, cfg.Scene.Components[0].Roll)
	}
	if cfg.Scene.Components[0].Shape.Params["width"] != beforeW {
		t.Errorf("scan mutated the shape width: %v -> %v", beforeW, cfg.Scene.Components[0].Shape.Params["width"])
	}
}
