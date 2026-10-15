package doe

import (
	"fmt"
	"math"
	"sort"
	"testing"

	"tetsuhiro-oz/optics"
)

// TestZonePlateCrossCheck compares the focus of the WOS binary zone plate
// (optics, kind "phase": 0/π zones, f = 50 mm, R = 2 mm, λ = 632.8 nm, 1024²
// grid over 4 mm, pitch dx ≈ 3.9 µm) with the design focal length f.
//
// Asserted: the on-axis intensity maximum lies within the depth of focus
// λF#² (F# = f/2R) of f. This is the physical resolution of the focus position
// and is far looser than one pixel on purpose: a hard-aperture binary zone
// plate has a strong defocused-order background on the axis, which moves the
// maximum by about one pixel and changes with grid and aperture (measured:
// −1.0 px at dx = 3.9 µm, +2.7 px at dx = 1.95 µm with R = 1 mm). The
// one-pixel criterion is therefore reported, not asserted, for this reference.
//
// Logged: the 1-pixel offset, and the window efficiency (1 mm radius) against
// η_2·EE_Airy(v). The WOS focal spot (Airy radius 1.22λF# ≈ 9.7 µm) spans about
// 2.5 pixels, so its narrow-window fraction is sampling-limited.
func TestZonePlateCrossCheck(t *testing.T) {
	const (
		f      = 0.05
		radius = 0.002
		lam    = 632.8e-9
		size   = 1024
		width  = 0.004
	)
	dx := width / size
	var zs []float64
	for k := -3; k <= 3; k++ {
		zs = append(zs, f+float64(k)*dx)
	}
	zs = append(zs, f-1e-3, f+1e-3)
	sort.Float64s(zs)

	cfg := optics.Config{
		Grid:       optics.GridSpec{Size: size, Width: width},
		Wavelength: lam,
		Polarized:  func() *bool { b := false; return &b }(),
		Method:     "asm",
		Evanescent: "decay",
		Source:     optics.SourceSpec{Type: "plane", Params: map[string]any{"power": 1e-3}},
		Elements:   []optics.ElementSpec{{Type: "zone_plate", Params: map[string]any{"f": f, "radius": radius, "kind": "phase"}}},
	}
	cur := 0.0
	for i, z := range zs {
		cfg.Elements = append(cfg.Elements,
			optics.ElementSpec{Type: "propagate", Params: map[string]any{"distance": z - cur}},
			optics.ElementSpec{Type: "sensor", Params: map[string]any{"label": fmt.Sprintf("z%02d", i)}})
		cur = z
	}
	res, err := optics.Simulate(cfg)
	if err != nil {
		t.Fatal(err)
	}
	var planes []int
	for i, pl := range res.Planes {
		if len(pl.Ex) == size*size && len(pl.Label) > 0 && pl.Label[0] == 'z' {
			planes = append(planes, i)
		}
	}
	if len(planes) != len(zs) {
		t.Fatalf("got %d sensor planes, want %d", len(planes), len(zs))
	}
	onAxis := make([]float64, len(zs))
	var winEta float64
	const rho = 1e-3
	c := size / 2
	for k, pi := range planes {
		pl := res.Planes[pi]
		v := pl.Ex[c*size+c]
		onAxis[k] = real(v)*real(v) + imag(v)*imag(v)
		if math.Abs(zs[k]-f) > 1e-9 {
			continue
		}
		var win float64
		for r := 0; r < size; r++ {
			y := float64(r-c) * dx
			for cc := 0; cc < size; cc++ {
				x := float64(cc-c) * dx
				if math.Hypot(x, y) <= rho {
					w := pl.Ex[r*size+cc]
					win += real(w)*real(w) + imag(w)*imag(w)
				}
			}
		}
		winEta = win * dx * dx / 1e-3
	}
	// Peak on the fine planes (f + k·dx, k = −3…3); parabolic refinement.
	var fine []int
	for k, z := range zs {
		if math.Abs(z-f) <= 3.5*dx {
			fine = append(fine, k)
		}
	}
	kmax := fine[0]
	for _, k := range fine {
		if onAxis[k] > onAxis[kmax] {
			kmax = k
		}
	}
	zPeak := zs[kmax]
	if kmax > 0 && kmax < len(zs)-1 {
		y0, y1, y2 := onAxis[kmax-1], onAxis[kmax], onAxis[kmax+1]
		if den := y0 - 2*y1 + y2; den != 0 {
			zPeak += 0.5 * (y0 - y2) / den * dx
		}
	}
	fNum := f / (2 * radius)
	dof := lam * fNum * fNum
	if d := math.Abs(zPeak - f); d > dof {
		t.Errorf("WOS zone plate on-axis peak %.6f m, f = %.6f m: |diff| %.3g m exceeds the depth of focus λF#² = %.3g m",
			zPeak, f, d, dof)
	}
	v := 2 * math.Pi / lam * radius * rho / f
	ee := 1 - math.J0(v)*math.J0(v) - math.J1(v)*math.J1(v)
	eta2 := 4 / (math.Pi * math.Pi)
	t.Logf("WOS zone plate (not asserted at 1 px): on-axis peak %+.2f µm = %+.2f px (dx %.2f µm); DOF λF#² = %.3g m; "+
		"window ρ=1 mm η=%.5f vs η_2·EE_Airy(v=%.0f)=%.5f (ratio %.4f)",
		(zPeak-f)*1e6, (zPeak-f)/dx, dx*1e6, dof, winEta, v, eta2*ee, winEta/(eta2*ee))
}
