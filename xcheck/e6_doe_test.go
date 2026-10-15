package xcheck

import (
	"fmt"
	"math"
	"testing"

	"tetsuhiro-oz/optics"
)

// E6 (plan §6.4): the wave kernel's diffractive elements against their design
// formulas — the "geometric ↔ wave" cross-check for DOEs.
//
//   - Zone plate: the WOS binary phase zone plate (`zone_plate`, kind
//     "phase") must focus at its design focal length f within 1 % of f.
//   - Grating: a fully blazed grating (`grating`, kind "blazed", 2π blaze
//     depth) must send the first order to the design direction
//     sin θ = λ/Λ: after propagating L the transverse centroid of the beam
//     must sit at L·tan θ within 1 %.

// ptrBool is a local helper for the *bool Polarized field of optics.Config.
func ptrBool(b bool) *bool { return &b }

// TestE6ZonePlateFocusVsDesign measures the on-axis focus of a WOS phase zone
// plate (f = 50 mm, radius 2 mm, λ = 632.8 nm) by scanning planes around the
// design focal length and asserts |z_peak − f| ≤ 1 % f. The finer pixel-scale
// evidence (sub-pixel offset, depth of focus) is logged, not asserted.
func TestE6ZonePlateFocusVsDesign(t *testing.T) {
	const (
		f      = 0.05
		radius = 0.002
		lam    = 632.8e-9
		size   = 1024
		width  = 0.004
	)
	dx := width / size
	// Fine planes at f ± k·dx (k = 0, 1, 2) bracket the sub-pixel focus; the
	// coarse planes at ±1 mm catch a focus that is grossly off design.
	zs := []float64{f - 1e-3, f - 2*dx, f - dx, f, f + dx, f + 2*dx, f + 1e-3}

	cfg := optics.Config{
		Grid:       optics.GridSpec{Size: size, Width: width},
		Wavelength: lam,
		Polarized:  ptrBool(false),
		Method:     "asm",
		Evanescent: "decay",
		Source:     optics.SourceSpec{Type: "plane", Params: map[string]any{"power": 1e-3}},
		Elements:   []optics.ElementSpec{{Type: "zone_plate", Params: map[string]any{"f": f, "radius": radius, "kind": "phase"}}},
	}
	cur := 0.0
	for i, z := range zs {
		cfg.Elements = append(cfg.Elements,
			optics.ElementSpec{Type: "propagate", Params: map[string]any{"distance": z - cur}},
			optics.ElementSpec{Type: "sensor", Params: map[string]any{"label": fmt.Sprintf("e6zp%02d", i)}})
		cur = z
	}
	res, err := optics.Simulate(cfg)
	if err != nil {
		t.Fatal(err)
	}
	var planes []*optics.Plane
	for _, pl := range res.Planes {
		if len(pl.Label) >= 4 && pl.Label[:4] == "e6zp" {
			planes = append(planes, pl)
		}
	}
	if len(planes) != len(zs) {
		t.Fatalf("got %d sensor planes, want %d", len(planes), len(zs))
	}
	c := size / 2
	best, bestI := 0, -1.0
	onAxis := make([]float64, len(zs))
	for k, pl := range planes {
		v := pl.Ex[c*size+c]
		onAxis[k] = real(v)*real(v) + imag(v)*imag(v)
		if onAxis[k] > bestI {
			best, bestI = k, onAxis[k]
		}
	}
	zPeak := zs[best]
	if d := math.Abs(zPeak - f); d > 0.01*f {
		t.Errorf("zone plate: on-axis peak at %.6f m, design f = %.6f m, |diff| %.3g m > 1%% (%.3g m)",
			zPeak, f, d, 0.01*f)
	}
	fNum := f / (2 * radius)
	dof := lam * fNum * fNum
	t.Logf("zone plate (f=%.0f mm, R=%.0f mm, F#=%.1f): on-axis peak %+.3f µm = %+.2f px (dx=%.2f µm); "+
		"|diff|/f = %.4f%%; depth of focus λF#² = %.3g µm",
		f*1e3, radius*1e3, fNum, (zPeak-f)*1e6, (zPeak-f)/dx, dx*1e6, 100*math.Abs(zPeak-f)/f, dof*1e6)
}

// TestE6BlazedGratingFirstOrderVsDesign propagates a Gaussian beam
// (w0 = 1 mm) through a fully blazed grating (period Λ = 100 µm, blaze depth
// 2π) and asserts that the transverse centroid of the first order after
// L = 50 mm equals L·tan θ within 1 %, θ = asin(λ/Λ).
func TestE6BlazedGratingFirstOrderVsDesign(t *testing.T) {
	const (
		lam    = 632.8e-9
		period = 100e-6
		L      = 0.05
		size   = 1024
		width  = 0.008
		waist  = 1e-3
	)
	sinT := lam / period
	tanT := sinT / math.Sqrt(1-sinT*sinT)
	want := L * tanT

	cfg := optics.Config{
		Grid:       optics.GridSpec{Size: size, Width: width},
		Wavelength: lam,
		Polarized:  ptrBool(false),
		Method:     "asm",
		Evanescent: "decay",
		Source:     optics.SourceSpec{Type: "gaussian", Params: map[string]any{"power": 1e-3, "waist": waist}},
		Elements: []optics.ElementSpec{
			{Type: "grating", Params: map[string]any{"kind": "blazed", "period": period, "blaze_depth": 2 * math.Pi}},
			{Type: "propagate", Params: map[string]any{"distance": L}},
			{Type: "sensor", Params: map[string]any{"label": "e6g"}},
		},
	}
	res, err := optics.Simulate(cfg)
	if err != nil {
		t.Fatal(err)
	}
	var pl *optics.Plane
	for _, p := range res.Planes {
		if p.Label == "e6g" {
			pl = p
		}
	}
	if pl == nil {
		t.Fatal("sensor plane e6g not found")
	}
	cx := pl.Stats.CentroidX
	if d := math.Abs(cx - want); d > 0.01*want {
		t.Errorf("grating: first-order centroid at %+.6f mm, design %+.6f mm (L·tan θ, sin θ = λ/Λ = %.6f), |diff| %.3g mm > 1%% (%.3g mm)",
			cx*1e3, want*1e3, sinT, d*1e3, 0.01*want*1e3)
	}
	// The first order must also stay on the design side (sign) of the axis.
	if cx <= 0 {
		t.Errorf("grating: first-order centroid %+.6f mm is on the wrong side of the axis", cx*1e3)
	}
	t.Logf("blazed grating (Λ=%.0f µm, λ=%.1f nm): first-order centroid at %+.6f mm vs design %+.6f mm; "+
		"relative deviation %.4f%%; power conservation dev=%.3g",
		period*1e6, lam*1e9, cx*1e3, want*1e3, 100*math.Abs(cx-want)/want, math.Abs(pl.Stats.Power-1e-3)/1e-3)
}
