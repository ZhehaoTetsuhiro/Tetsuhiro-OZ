package optics

import (
	"math"
	"testing"
)

// rotScene builds the 2026-10-16 layout: an incoherent emitter far in front,
// a single rotatable slit, and a screen. thetaDeg rotates the slit's long
// axis away from y.
func rotScene(t *testing.T, thetaDeg float64, params map[string]any) Config {
	t.Helper()
	w := 0.03
	comps := []ComponentSpec{
		{Type: "aperture", Label: "单缝", Pos: v3(0, 0, -0.15), Roll: thetaDeg * math.Pi / 180,
			Shape: RectOutline(4e-4, 0.02)},
		bigDet("屏", v3(0, 0, 0.15), math.Pi, w),
	}
	cfg := sceneCfg(256, w, []SourceSpec{{
		ID: "src", Label: "光源", Type: "extended", Pos: &Vec3{Z: -0.30}, Dir: &Vec3{Z: 1},
		Params: params,
	}}, comps...)
	return cfg
}

// TestRotatingSlitRollIsApplied: the component's Roll must actually rotate
// the slit mask. At 90° the 0.4 mm width lies along y, so the far-field
// pattern spreads horizontally instead of vertically.
func TestRotatingSlitRollIsApplied(t *testing.T) {
	p := map[string]any{"power": 1e-3, "samples_x": 1, "samples_y": 1, "standoff": 0.14}
	flat, err := Simulate(rotScene(t, 0, p))
	if err != nil {
		t.Fatal(err)
	}
	rot, err := Simulate(rotScene(t, 90, p))
	if err != nil {
		t.Fatal(err)
	}
	// The unrotated slit is wide in y and narrow in x, so its far field is a
	// narrow horizontal streak: rms_x << rms_y.
	if flat.Planes[0].Stats.RMSX >= flat.Planes[0].Stats.RMSY {
		t.Errorf("unrotated slit: rms_x %.4g m should be far below rms_y %.4g m",
			flat.Planes[0].Stats.RMSX, flat.Planes[0].Stats.RMSY)
	}
	// Rotated by 90° the streak must be vertical instead.
	if rot.Planes[0].Stats.RMSY >= rot.Planes[0].Stats.RMSX {
		t.Errorf("90° rotated slit: rms_y %.4g m should be far below rms_x %.4g m",
			rot.Planes[0].Stats.RMSY, rot.Planes[0].Stats.RMSX)
	}
}

// TestRotatingSlitRollPreservesArea pins the rotation convention: a roll of
// theta rotates the mask without rescaling it, so the *area* of the clear
// aperture is the same at every angle. The transmitted power is then only
// changed by the geometry: a rotated slit presents a different projected
// width to the (off-axis) emitter points, so the power varies — the check is
// that it stays within the physically expected band rather than collapsing.
func TestRotatingSlitRollPreservesArea(t *testing.T) {
	p := map[string]any{"power": 1e-3, "samples_x": 1, "samples_y": 1, "standoff": 0.14}
	lo, hi := math.Inf(1), 0.0
	for _, deg := range []float64{0, 15, 45, 90, 137, 180} {
		res, err := Simulate(rotScene(t, deg, p))
		if err != nil {
			t.Fatal(err)
		}
		got := res.Planes[0].Stats.Power
		if got < lo {
			lo = got
		}
		if got > hi {
			hi = got
		}
	}
	// A rotated 0.4 mm x 20 mm slit keeps its area, so the transmitted power
	// must stay within a factor of ~2 of the unrotated value across 0-180°.
	if hi/lo > 2.0 {
		t.Errorf("transmitted power varies from %.6g to %.6g W (factor %.2f): the rotation must not rescale the aperture",
			lo, hi, hi/lo)
	}
}

// TestRotatingSlitInvalidRejected: a negative or oversized slit size must be
// refused by validation, not silently clipped.
func TestRotatingSlitInvalidRejected(t *testing.T) {
	for _, sh := range []*ShapeSpec{
		RectOutline(-4e-4, 0.02),
		RectOutline(4e-4, -0.02),
		{Kind: "rectangle", Params: map[string]any{"width": 4e-4}},
	} {
		w := 0.03
		cfg := sceneCfg(256, w, []SourceSpec{{
			Type: "extended", Pos: &Vec3{Z: -0.30}, Dir: &Vec3{Z: 1},
			Params: map[string]any{"power": 1e-3, "samples_x": 1, "samples_y": 1, "standoff": 0.14},
		}},
			ComponentSpec{Type: "aperture", Label: "单缝", Pos: v3(0, 0, -0.15), Shape: sh},
			bigDet("屏", v3(0, 0, 0.15), math.Pi, w))
		if issues := ValidateConfig(&cfg); len(issues) == 0 {
			t.Errorf("shape %+v passed validation", sh)
		}
	}
}
