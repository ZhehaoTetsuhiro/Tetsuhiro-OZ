package xcheck

import (
	"math"
	"testing"

	"tetsuhiro-oz/geo"
	"tetsuhiro-oz/optics"
)

// E5 (plan §6.4): focal length of an ideal thin lens, the optics lens element
// against the geometric paraxial EFL.
//
// Geometric side: a zero-thickness biconvex thin lens (n = 1.5, R1 = +100 mm,
// R2 = −100 mm, air on both sides) traced by geo's paraxial engine.
// Its radii follow from the design focal length f0 = 100 mm by
// R = 2(n−1)·f0, and the analytic thin-lens power 1/f = (n−1)(1/R1 − 1/R2)
// is checked against the engine to 1e-12.
//
// Wave side: the optics "lens" element built with f = f0 (SI, metres) is applied
// to a unit plane wave on a 256 × 256 grid (pixel 2 µm). Its wavefront
// φ(r) = −k·r²/(2f) is fitted by least squares through the origin on r ≤ 0.2 mm
// (|φ| ≤ 2.2 rad, so no phase wrapping); the focal length is f_w = −k/(2c).
//
// Criterion: |f_w − EFL_geo| / EFL_geo ≤ 1e-6.
func TestE5ThinLensEFL(t *testing.T) {
	const (
		f0 = 100.0 // mm, design focal length
		n  = 1.5
	)
	R := 2 * (n - 1) * f0
	sys := &geo.System{
		Surfaces: []geo.Surface{
			geo.NewSurface(geo.Sphere{R: R}, 10),
			geo.NewSurface(geo.Sphere{R: -R}, 10),
		},
		Thickness: []float64{0},
		Mediums:   []geo.Medium{geo.Air{}, geo.ConstantIndex(n), geo.Air{}},
	}
	if err := sys.Validate(); err != nil {
		t.Fatal(err)
	}
	par, err := sys.Paraxial(geo.Wavelengthd, math.Inf(-1))
	if err != nil {
		t.Fatal(err)
	}
	eflGeo := par.EFL
	analytic := 1 / ((n - 1) * (1/R - 1/(-R)))
	if e := relErr(eflGeo, analytic); e > 1e-12 {
		t.Errorf("E5: geo EFL %.15f vs analytic thin-lens %.15f: relative difference %.3e", eflGeo, analytic, e)
	}

	el, err := optics.NewElement(optics.ElementSpec{Type: "lens", Params: map[string]any{"f": f0 * 1e-3}})
	if err != nil {
		t.Fatal(err)
	}
	const N = 256
	dxM := 2e-6
	field := optics.NewField(N, dxM, false)
	for i := range field.Ex {
		field.Ex[i] = 1
	}
	ctx := &optics.Context{Wavelength: dLineM}
	if err := el.Apply(field, ctx); err != nil {
		t.Fatal(err)
	}
	k := 2 * math.Pi / dLineM
	var sRP, sRR float64 // Σ r²·φ and Σ r⁴ for φ = c·r²
	nFit := 0
	for j := 0; j < N; j++ {
		for i := 0; i < N; i++ {
			x, y := field.X(i), field.Y(j)
			r2 := x*x + y*y
			if r2 == 0 || r2 > 0.2e-3*0.2e-3 {
				continue
			}
			v := field.Ex[j*N+i]
			ph := math.Atan2(imag(v), real(v))
			sRP += r2 * ph
			sRR += r2 * r2
			nFit++
		}
	}
	c := sRP / sRR
	fWave := -k / (2 * c) * 1e3 // mm
	dev := relErr(fWave, eflGeo)
	t.Logf("E5 geo thin lens: EFL = %.12f mm (analytic %.12f mm)", eflGeo, analytic)
	t.Logf("E5 optics lens element (f = %.1f mm): wavefront-fitted f = %.12f mm over %d samples; |Δf|/f = %.3e",
		f0, fWave, nFit, dev)
	if dev > 1e-6 {
		t.Errorf("E5: wave focal length %.9f mm vs geo EFL %.9f mm: relative difference %.3e > 1e-6", fWave, eflGeo, dev)
	}
}
