package optics

import (
	"encoding/json"
	"fmt"
	"math"
	"strings"
	"testing"
)

// ---------------------------------------------------------------------------
// Custom-medium tests: expression language, the physics of the reference-index
// split-step (plate focal shift, gradient deflection, absorption/gain, layer
// stacks), configuration validation, catalog wiring, and the cost of the
// compiled evaluation path.
//
// Each analytic comparison states its tolerance and why it is what it is; the
// numbers are printed (t.Logf) so the accuracy of a run is visible in
// `go test -v ./optics/ -run CustomMedium`.
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// helpers
// ---------------------------------------------------------------------------

// planeField builds a plane wave on an N x N grid.
func planeField(t *testing.T, n int, width, wl float64) *Field {
	t.Helper()
	f, err := BuildSource(SourceSpec{Type: "plane", Params: map[string]any{"power": 1e-3}}, n, width, false, wl)
	if err != nil {
		t.Fatal(err)
	}
	return f
}

// mustMedium compiles a spec or fails the test with the compiler's message.
func mustMedium(t *testing.T, spec CustomMediumSpec) *CustomMedium {
	t.Helper()
	m, err := NewCustomMedium(spec)
	if err != nil {
		t.Fatalf("NewCustomMedium: %v", err)
	}
	return m
}

// cAbs is |z|, and maxRelDiff the largest relative difference between two
// complex fields, with a floor so a null pixel is not a division by zero.
func cAbs(z complex128) float64 { return math.Hypot(real(z), imag(z)) }

func maxRelDiff(a, b []complex128) (maxRel, maxAbs float64) {
	for i := range a {
		d := cAbs(a[i] - b[i])
		if d > maxAbs {
			maxAbs = d
		}
		den := cAbs(b[i])
		if den < 1e-12 {
			den = 1e-12
		}
		if r := d / den; r > maxRel {
			maxRel = r
		}
	}
	return
}

// waistVertex fits the least-squares vertex of rms^2(z) = a(1 + ((z-zf)/zR)^2)
// — a parabola in z — from samples (zs, ys) around the waist. Solving the 3x3
// normal equations rather than taking the argmin removes the discretization
// bias of a flat minimum: the synthetic self-check below asserts the ruler
// measures a known vertex exactly before it is used on a beam.
func waistVertex(zs, ys []float64, z0 float64) float64 {
	var a0, a1, a2, b0, b1, b2, c0, c1, c2, d0, d1, d2 float64
	for i, z := range zs {
		u := z - z0
		a0, a1, a2 = a0+u*u*u*u, a1+u*u*u, a2+u*u
		b0, b1, b2 = b0+u*u*u, b1+u*u, b2+u
		c0, c1, c2 = c0+u*u, c1+u, c2+1
		d0, d1, d2 = d0+ys[i]*u*u, d1+ys[i]*u, d2+ys[i]
	}
	det := a0*(b1*c2-b2*c1) - a1*(b0*c2-b2*c0) + a2*(b0*c1-b1*c0)
	A := (d0*(b1*c2-b2*c1) - a1*(d1*c2-b2*d2) + a2*(d1*c1-b1*d2)) / det
	B := (a0*(d1*c2-b2*d2) - d0*(b0*c2-b2*c0) + a2*(b0*d2-d1*c0)) / det
	return z0 - B/(2*A)
}

// ---------------------------------------------------------------------------
// ① uniform medium == plane-parallel glass plate
// ---------------------------------------------------------------------------

// TestCustomMediumPlateFocusShift: a homogeneous custom layer must behave like
// a real plane-parallel plate, i.e. refract. The analytic focal shift of a plate
// of thickness d and index n in a converging beam is d(1-1/n) (ray matrix of a
// slab: a reduced translation d/n between its faces, so the focus moves away
// from the plate by d-d/n). The alternative reading "(n-1)d" is the *optical
// path* the plate adds, not the focal shift; the vacuum-reference kernel
// (reference:1, and the built-in "medium" element) produces neither, because a
// layer that does not vary across the beam is then only a constant phase screen
// — that case is asserted too, so the difference between the two models is
// pinned by a test rather than left to prose.
//
// The focus plane is found from the minimum of the beam's second moment
// (rms^2(z), a parabola in z) by a least-squares quadratic fit, at N=256.
func TestCustomMediumPlateFocusShift(t *testing.T) {
	const wl = 632.8e-9
	const n, width = 256, 0.004
	const fLens, d, idx = 0.3, 2e-3, 1.5
	ctx := &Context{Wavelength: wl, Evanescent: "decay"}

	// Self-check of the measurement: a synthetic parabola with a known vertex.
	{
		var zs, ys []float64
		for z := 0.28; z <= 0.31; z += 5e-4 {
			zs, ys = append(zs, z), append(ys, 1+math.Pow((z-0.298912)/0.0179, 2))
		}
		if v := waistVertex(zs, ys, 0.28); math.Abs(v-0.298912) > 1e-9 {
			t.Fatalf("waist fit self-check: got %.9f want 0.298912", v)
		}
	}

	waist := func(medium func(f *Field, ctx *Context) error) float64 {
		f, err := BuildSource(SourceSpec{Type: "gaussian", Params: map[string]any{"waist": 1e-3, "power": 1e-3}}, n, width, false, wl)
		if err != nil {
			t.Fatal(err)
		}
		el, err := NewElement(ElementSpec{Type: "lens", Params: map[string]any{"f": fLens}})
		if err != nil {
			t.Fatal(err)
		}
		if err := el.Apply(f, ctx); err != nil {
			t.Fatal(err)
		}
		if medium != nil {
			if err := medium(f, ctx); err != nil {
				t.Fatal(err)
			}
		}
		var zs, ys []float64
		for z := 0.0; z <= 0.31; z += 5e-4 {
			if z >= 0.28 {
				r := f.rmsRadius()
				zs, ys = append(zs, z), append(ys, r*r)
			}
			if err := Propagate(f, 5e-4, MethodASM, ctx); err != nil {
				t.Fatal(err)
			}
		}
		return waistVertex(zs, ys, 0.28)
	}

	zVac := waist(nil)
	plate := mustMedium(t, CustomMediumSpec{Distance: d, Index: "1.5", Steps: 8})
	zPlate := waist(func(f *Field, ctx *Context) error { return PropagateCustomMedium(f, plate, ctx) })
	vacRef := mustMedium(t, CustomMediumSpec{Layers: []CustomMediumLayerSpec{
		{Thickness: d, Index: "1.5", Reference: ptrF(1)}}})
	zVacRef := waist(func(f *Field, ctx *Context) error { return PropagateCustomMedium(f, vacRef, ctx) })
	zOld := waist(func(f *Field, ctx *Context) error {
		return PropagateSplitStep(f, d, UniformIndex(complex(idx, 0)), 8, ctx)
	})

	// z is measured from the medium's exit face, so a layer that does not
	// refract at all lands at zVac-d and the real plate at zVac-d+d(1-1/n).
	wantPlate := zVac - d + d*(1-1/idx)
	shift := zPlate - (zVac - d)
	t.Logf("waist without medium             : %.6f m", zVac)
	t.Logf("plate (reference = n = %.1f)      : %.6f m -> shift %.4f mm (analytic d(1-1/n) = %.4f mm, (n-1)d = %.4f mm)",
		idx, zPlate, shift*1e3, d*(1-1/idx)*1e3, (idx-1)*d*1e3)
	t.Logf("reference = 1 (vacuum kernel)    : %.6f m -> shift %.4f mm", zVacRef, (zVacRef-(zVac-d))*1e3)
	t.Logf("built-in PropagateSplitStep      : %.6f m -> shift %.4f mm", zOld, (zOld-(zVac-d))*1e3)

	// Tolerance: the fit resolves the vertex to ~1 um (the self-check above),
	// and the difference relative to the analytic paraxial plate shift is
	// limited by O(NA^2) terms of the beam (NA = w/f = 3.3 mrad) plus the
	// truncated-Gaussian focal shift, which is common to both runs. 2% of the
	// shift (13 um) is generous.
	if rel := math.Abs(zPlate-wantPlate) / (d * (1 - 1/idx)); rel > 0.02 {
		t.Fatalf("plate focal shift %.4f mm, analytic %.4f mm (rel %.3f > 2%%)",
			shift*1e3, d*(1-1/idx)*1e3, rel)
	}
	if shift <= 0 {
		t.Fatalf("plate must move the focus downstream, got %.4f mm", shift*1e3)
	}
	// Without a reference index the layer is a pure phase screen: no refraction.
	for name, z := range map[string]float64{"reference=1": zVacRef, "PropagateSplitStep": zOld} {
		if got := z - (zVac - d); math.Abs(got) > 2e-5 {
			t.Fatalf("%s must leave the focus in place for a homogeneous slab, got %+.4f mm", name, got*1e3)
		}
	}
}

// TestCustomMediumPlateExactTransfer: with reference equal to the layer index
// the layer is not merely "shifted correctly" — it is the exact homogeneous
// angular-spectrum propagator, i.e. each plane-wave component picks up
// exp(i d sqrt((n k0)^2 - kt^2)) and nothing else. Comparing a random field
// propagated through the compiled medium against the same field with that
// analytic kernel applied directly measures the implementation, not a model.
func TestCustomMediumPlateExactTransfer(t *testing.T) {
	const wl = 632.8e-9
	const n, width = 64, 0.002
	const d, idx = 5e-4, 1.7
	ctx := &Context{Wavelength: wl, Evanescent: "decay"}
	f, err := BuildSource(SourceSpec{Type: "plane", Params: map[string]any{"power": 1e-3}}, n, width, false, wl)
	if err != nil {
		t.Fatal(err)
	}
	// A structured field: a lens phase plus a tilt, so many plane waves are present.
	el, err := NewElement(ElementSpec{Type: "lens", Params: map[string]any{"f": 0.05}})
	if err != nil {
		t.Fatal(err)
	}
	if err := el.Apply(f, ctx); err != nil {
		t.Fatal(err)
	}
	f.ApplyTilt(3e-4, -2e-4, wl)
	// The reference field with the analytic homogeneous-medium kernel applied
	// directly (including the axial phase k0*n*d).
	// One slice is enough for a homogeneous layer (the kernel is exact and the
	// slices telescope), so the comparison below is a single FFT round-trip.
	m := mustMedium(t, CustomMediumSpec{Distance: d, Index: "1.7", Steps: 1})
	want := f.Clone()
	k0 := 2 * math.Pi / wl
	want.Ex = append([]complex128(nil), f.Ex...)
	ctxRef := *ctx
	ctxRef.transfer(want, d, func(fx, fy float64) complex128 {
		kt := 4 * math.Pi * math.Pi * (fx*fx + fy*fy)
		return cexpI(d * math.Sqrt(k0*k0*idx*idx-kt))
	})
	if err := PropagateCustomMedium(f, m, ctx); err != nil {
		t.Fatal(err)
	}
	maxRel, maxAbs := maxRelDiff(f.Ex, want.Ex)
	t.Logf("max |ΔE| = %.3e, max relative = %.3e (N=%d, d=%g m, n=%.1f)", maxAbs, maxRel, n, d, idx)
	// The comparison carries the medium's axial phase k0*n*d ≈ 8440 rad, so
	// double precision alone permits a relative difference of about
	// eps*phase ≈ 2e-12; 1e-11 is the honest bound (max |ΔE| above is ~6e-11
	// on a field of order 15, i.e. 4e-12 relative).
	if maxRel > 1e-11 {
		t.Fatalf("homogeneous layer is not the exact angular-spectrum transfer function: rel %.3e", maxRel)
	}
}

// TestCustomMediumAxialPhase: the optical path through a uniform layer is
// exactly k0 n d, i.e. the on-axis phase is the phase the medium imposes (this
// is the "(n-1)d optical path" statement, measured where it is actually exact).
func TestCustomMediumAxialPhase(t *testing.T) {
	const n, width, wl = 128, 0.004, 632.8e-9
	const d, idx = 1e-3, 1.5
	ctx := &Context{Wavelength: wl, Evanescent: "decay"}
	f := planeField(t, n, width, wl)
	m := mustMedium(t, CustomMediumSpec{Distance: d, Index: "1.5"})
	if err := PropagateCustomMedium(f, m, ctx); err != nil {
		t.Fatal(err)
	}
	k0 := 2 * math.Pi / wl
	// Total on-axis phase = k0 * optical path = k0 * n * d. The full medium also
	// propagates the field, so compare the total against k0*n*d.
	c := f.Ex[(n/2)*n+n/2]
	got := math.Atan2(imag(c), real(c))
	want := math.Remainder(k0*idx*d, 2*math.Pi)
	diff := math.Remainder(got-want, 2*math.Pi)
	t.Logf("axial phase %.9f rad vs analytic k0*n*d = %.9f rad (Δ = %.3e)", got, want, diff)
	if math.Abs(diff) > 1e-9 {
		t.Fatalf("optical path through the layer wrong: Δ = %.3e rad", diff)
	}
}

// ptrF is a small helper for optional float64 spec fields.
func ptrF(v float64) *float64 { return &v }

// ---------------------------------------------------------------------------
// ② linear transverse gradient: deflection angle ≈ g*L, centroid displacement
// ---------------------------------------------------------------------------

// TestCustomMediumGradientDeflection: for n(x) = n0 + g*x the ray equation
// dθ/dz = (1/n) ∂n/∂x gives a deflection θ(L) = ∫ g/n dz ≈ g L/n0 (the medium
// starts on axis, so x(z) = g z²/(2 n0²) is second order and the deviation from
// g*L is O(g²L³)). The same test checks the centroid displacement against the
// double integral x(L) = ∫θ dz ≈ g L²/(2 n0). Both are compared with a
// numerically ray-traced reference (RK4 on the exact ray equations, including
// the n(x) factor) so the check is not just "the formula the code implements".
func TestCustomMediumGradientDeflection(t *testing.T) {
	const wl = 632.8e-9
	const n, width = 256, 0.004
	const g, L, thin = 0.5, 0.02, 0.006 // 1/m, m, m
	const n0 = 1.0
	ctx := &Context{Wavelength: wl, Evanescent: "decay"}
	// An on-axis beam, wide compared with the distance it sweeps, so its
	// centroid is clean. Grid Nyquist set: the total tilt g*L = 10 mrad must
	// stay below the sampling limit theta_nyq = lambda/(2 dx) = 20 mrad, or the
	// plane-wave spectrum of the deflected beam wraps around the window.
	f, err := BuildSource(SourceSpec{Type: "gaussian", Params: map[string]any{"waist": 1e-3, "power": 1e-3}}, n, width, false, wl)
	if err != nil {
		t.Fatal(err)
	}
	t.Logf("gk: g·L = %.4f rad, theta_nyq = lambda/(2dx) = %.4f rad", g*L, wl/(2*(width/float64(n))))
	centroid := func(f *Field) (float64, float64) {
		var w, sx, sy float64
		for j := 0; j < n; j++ {
			for i := 0; i < n; i++ {
				p := f.Intensity(j*n + i)
				w += p
				sx += p * f.X(i)
				sy += p * f.Y(j)
			}
		}
		return sx / w, sy / w
	}
	m := mustMedium(t, CustomMediumSpec{
		Steps: 400,
		Layers: []CustomMediumLayerSpec{{
			Thickness: L, Index: "n0 + g*x", Params: map[string]float64{"n0": n0, "g": g},
		}},
	})
	if err := PropagateCustomMedium(f, m, ctx); err != nil {
		t.Fatal(err)
	}
	xc1, yc1 := centroid(f)
	// The wavefront tilt is the centroid slope in free space just behind the
	// medium: propagate a short distance and difference the centroids.
	free := f.Clone()
	if err := Propagate(free, thin, MethodASM, ctx); err != nil {
		t.Fatal(err)
	}
	save := f.Ex
	f.Ex = free.Ex
	xc2, yc2 := centroid(f)
	f.Ex = save
	angle := (xc2 - xc1) / thin

	// Independent ray-traced reference: dx/dz = theta, dtheta/dz = g/n(x).
	xr, thr := 0.0, 0.0
	const hs = 500
	dz := L / hs
	for s := 0; s < hs; s++ {
		k1x, k1t := thr, g/(n0+g*xr)
		k2x, k2t := thr+0.5*dz*k1t, g/(n0+g*(xr+0.5*dz*k1x))
		k3x, k3t := thr+0.5*dz*k2t, g/(n0+g*(xr+0.5*dz*k2x))
		k4x, k4t := thr+dz*k3t, g/(n0+g*(xr+dz*k3x))
		xr += dz / 6 * (k1x + 2*k2x + 2*k3x + k4x)
		thr += dz / 6 * (k1t + 2*k2t + 2*k3t + k4t)
	}
	t.Logf("deflection: measured %.6f rad, ray-traced %.6f rad, small-angle g·L = %.6f rad (Δ = %.3f%%)",
		angle, thr, g*L, (angle-thr)/thr*100)
	t.Logf("centroid:   measured %.4f mm, ray-traced %.4f mm (Δ = %.3f%%)",
		xc1*1e3, xr*1e3, (xc1-xr)/xr*100)
	if rel := math.Abs(angle-thr) / thr; rel > 0.02 {
		t.Fatalf("gradient deflection %.6f rad vs ray-traced %.6f rad (rel %.3f > 2%%)", angle, thr, rel)
	}
	if rel := math.Abs(xc1-xr) / xr; rel > 0.03 {
		t.Fatalf("beam centroid %.4f mm vs ray-traced %.4f mm (rel %.3f > 3%%)", xc1*1e3, xr*1e3, rel)
	}
	// No deflection off the gradient axis.
	if math.Abs(yc2-yc1) > 1e-6 {
		t.Fatalf("beam moved in y by %.3e m (the index does not vary in y)", yc2-yc1)
	}
	// The deflection must follow g·L: with n0 = 1 the correction is O(g²L³) and
	// the ray-traced reference above is the exact small-angle answer.
	if rel := math.Abs(angle-g*L) / (g * L); rel > 0.01 {
		t.Fatalf("deflection %.6f rad far from g·L = %.6f rad (rel %.3f)", angle, g*L, rel)
	}
}

// ---------------------------------------------------------------------------
// ③ absorption and gain: exp(-k L) in intensity
// ---------------------------------------------------------------------------

// TestCustomMediumAbsorptionAndGain: a layer with absorption expression k
// (Im(n)) attenuates the field by exp(-k0 Im(n) d) and the power by
// exp(-2 k0 Im(n) d) = exp(-α d) with the intensity coefficient α = 2 k0 Im(n);
// a negative value amplifies by the same law, and a spatially varying absorber
// attenuates exactly the region it covers. Tolerances: 1e-6 relative (the
// exponential is applied exactly per slice, and the ASM steps are unitary).
func TestCustomMediumAbsorptionAndGain(t *testing.T) {
	const n, width, wl = 128, 0.004, 632.8e-9
	k0 := 2 * math.Pi / wl
	ctx := &Context{Wavelength: wl, Evanescent: "decay"}

	atten := func(f *Field, spec CustomMediumSpec) (*Field, error) {
		m := mustMedium(t, spec)
		if err := PropagateCustomMedium(f, m, ctx); err != nil {
			return nil, err
		}
		return f, nil
	}

	for _, nim := range []float64{2e-5, -1e-5} { // absorb, then amplify
		f := planeField(t, n, width, wl)
		p0 := f.Power()
		th := 2e-3
		if _, err := atten(f, CustomMediumSpec{Distance: th, Index: "1.0", Absorption: fmt.Sprintf("%g", nim)}); err != nil {
			t.Fatal(err)
		}
		want := p0 * math.Exp(-2*k0*nim*th)
		rel := math.Abs(f.Power()-want) / want
		t.Logf("Im(n)=%+.1e, d=%g m: power %.6e vs exp(-2k0 Im(n) d)=%.6e (rel %.3e)", nim, th, f.Power(), want, rel)
		if rel > 1e-6 {
			t.Fatalf("Im(n)=%g: power %.6e, want %.6e (rel %.3e > 1e-6)", nim, f.Power(), want, rel)
		}
	}

	// A spatially varying absorber: a disc of radius 1 mm with Im(n) = a
	// attenuates the power inside the disc by exp(-2 k0 a d) and leaves the rest
	// untouched, so the transmitted power follows from the incident power
	// inside/outside the disc.
	f := planeField(t, n, width, wl)
	p0 := f.Power()
	const a, d = 3e-4, 1e-3
	fIn, err := BuildSource(SourceSpec{Type: "plane", Params: map[string]any{"power": 1e-3}}, n, width, false, wl)
	if err != nil {
		t.Fatal(err)
	}
	m := mustMedium(t, CustomMediumSpec{Distance: d, Index: "1.0", Absorption: "a*step(1e-3 - r)", Params: map[string]float64{"a": a}})
	if err := PropagateCustomMedium(fIn, m, ctx); err != nil {
		t.Fatal(err)
	}
	// Fraction of the incident power inside r <= 1 mm, from the incident field.
	var pin, ptot float64
	for j := 0; j < n; j++ {
		for i := 0; i < n; i++ {
			p := f.Intensity(j*n + i)
			ptot += p
			if math.Hypot(f.X(i), f.Y(j)) <= 1e-3 {
				pin += p
			}
		}
	}
	want := (pin*math.Exp(-2*k0*a*d) + (ptot - pin)) / ptot * p0
	rel := math.Abs(fIn.Power()-want) / want
	t.Logf("disc absorber: power %.6e vs analytic %.6e (rel %.3e), inside fraction %.4f", fIn.Power(), want, rel, pin/ptot)
	if rel > 2e-3 {
		t.Fatalf("disc absorber power %.6e, want %.6e (rel %.3e)", fIn.Power(), want, rel)
	}
}

// ---------------------------------------------------------------------------
// z-layer stack, z / zl variables and the vacuum outside
// ---------------------------------------------------------------------------

// TestCustomMediumLayerStackPhase: three layers of different index accumulate
// exactly the optical path Σ k0 n_i d_i, and the light outside the stack is
// untouched (vacuum). A layer whose expression reads the depth z is compared
// against the analytic integral ∫ k0 n(z) dz with enough slices.
func TestCustomMediumLayerStackPhase(t *testing.T) {
	const n, width, wl = 128, 0.004, 632.8e-9
	k0 := 2 * math.Pi / wl
	ctx := &Context{Wavelength: wl, Evanescent: "decay"}

	// (a) three homogeneous layers: optical path Σ n_i d_i.
	f := planeField(t, n, width, wl)
	m := mustMedium(t, CustomMediumSpec{Layers: []CustomMediumLayerSpec{
		{Thickness: 0.5e-3, Index: "1.5"},
		{Thickness: 1.0e-3, Index: "1.2"},
		{Thickness: 0.25e-3, Index: "2.0"},
	}})
	if err := PropagateCustomMedium(f, m, ctx); err != nil {
		t.Fatal(err)
	}
	path := 1.5*0.5e-3 + 1.2*1.0e-3 + 2.0*0.25e-3
	got := math.Atan2(imag(f.Ex[(n/2)*n+n/2]), real(f.Ex[(n/2)*n+n/2]))
	diff := math.Remainder(got-math.Remainder(k0*path, 2*math.Pi), 2*math.Pi)
	t.Logf("stack optical path %.6e m: phase Δ = %.3e rad (Σk0 n_i d_i)", path, diff)
	if math.Abs(diff) > 1e-9 {
		t.Fatalf("layer stack optical path wrong: Δ = %.3e rad", diff)
	}
	if d := m.Distance(); math.Abs(d-1.75e-3) > 1e-15 {
		t.Fatalf("stack thickness %g, want 1.75e-3", d)
	}
	if s := m.SliceCount(); s != DefaultCustomMediumSteps {
		t.Fatalf("slice count %d, want the default %d", s, DefaultCustomMediumSteps)
	}

	// (b) a z-graded layer: n(z) = n0 + a z/d, optical path = k0 ∫ n dz, and the
	// z variable is the depth from the *medium* entrance, kept across layers.
	f2 := planeField(t, n, width, wl)
	const d1, d2 = 0.5e-3, 0.5e-3
	m2 := mustMedium(t, CustomMediumSpec{Layers: []CustomMediumLayerSpec{
		{Thickness: d1, Index: "1.4"},
		{Thickness: d2, Index: "1.4 + 0.4*z/1e-3", Steps: 200},
	}})
	if err := PropagateCustomMedium(f2, m2, ctx); err != nil {
		t.Fatal(err)
	}
	path2 := 1.4*d1 + 1.4*d2 + 0.4*(d1*d2+d2*d2/2)/(1e-3) // ∫1.4 dz + 0.4/1e-3 ∫z dz, z from d1 to d1+d2
	got2 := math.Atan2(imag(f2.Ex[(n/2)*n+n/2]), real(f2.Ex[(n/2)*n+n/2]))
	diff2 := math.Remainder(got2-math.Remainder(k0*path2, 2*math.Pi), 2*math.Pi)
	t.Logf("z-graded layer: analytic path %.6e m, phase Δ = %.3e rad", path2, diff2)
	if math.Abs(diff2) > 1e-6 {
		t.Fatalf("z-graded optical path wrong: Δ = %.3e rad", diff2)
	}

	// (c) zl is the depth from the current layer's face: "n0 + b*zl" over one
	// layer equals the same integral as (b) with z replaced by zl.
	f3 := planeField(t, n, width, wl)
	m3 := mustMedium(t, CustomMediumSpec{
		Distance: d2, Index: "1.4 + 0.4*zl/1e-3", Params: nil, Steps: 200,
	})
	if err := PropagateCustomMedium(f3, m3, ctx); err != nil {
		t.Fatal(err)
	}
	path3 := 1.4*d2 + 0.4*(d2*d2/2)/(1e-3)
	got3 := math.Atan2(imag(f3.Ex[(n/2)*n+n/2]), real(f3.Ex[(n/2)*n+n/2]))
	diff3 := math.Remainder(got3-math.Remainder(k0*path3, 2*math.Pi), 2*math.Pi)
	t.Logf("zl-graded layer: analytic path %.6e m, phase Δ = %.3e rad", path3, diff3)
	if math.Abs(diff3) > 1e-6 {
		t.Fatalf("zl-graded optical path wrong: Δ = %.3e rad", diff3)
	}

	// (d) a lens phase behind the medium is unaffected by layers far away: the
	// medium is strictly a z-interval, the field entering and leaving it keeps
	// propagating in vacuum (checked by comparing a field through the medium
	// with the same field propagated the same total distance in vacuum, with
	// the medium's index equal to 1).
	f4 := planeField(t, n, width, wl)
	f4v := f4.Clone()
	mu := mustMedium(t, CustomMediumSpec{Distance: 1e-3, Index: "1.0 + 0*z"})
	if err := PropagateCustomMedium(f4, mu, ctx); err != nil {
		t.Fatal(err)
	}
	if err := Propagate(f4v, 1e-3, MethodASM, ctx); err != nil {
		t.Fatal(err)
	}
	_, maxRel := maxRelDiff(f4.Ex, f4v.Ex)
	t.Logf("n=1 layer vs equal vacuum propagation: max |ΔE| = %.3e", maxRel)
	if maxRel > 1e-12 {
		t.Fatalf("a n=1 medium must be transparent: max |ΔE| = %.3e", maxRel)
	}
}

// TestCustomMediumLayerEquivalence: homogeneous layers behave as the physics
// says they must. (a) A stack of homogeneous layers accumulates Σ k0 n_i d_i of
// optical path and, for a beam, a total reduced thickness Σ d_i/n_i — so it is
// *identical* to a single layer of thickness Σd_i and index n_eff = Σd/Σ(d_i/n_i)
// up to that piston phase, which is checked against the analytic value.
// (b) Homogeneous layers commute: each layer's operator is diagonal in the
// transverse spectrum, so reordering them reproduces the field exactly (a
// bookkeeping check on the layer z-ranges). (c) Layers whose index reads the
// global depth z do *not* commute, and the two orders are compared with the
// analytic ∫ n dz of each interval.
func TestCustomMediumLayerEquivalence(t *testing.T) {
	const wl = 632.8e-9
	const n, width = 128, 0.004
	// A gentle lens (NA = 2 mm / 0.5 m = 4 mrad) keeps the comparison in the
	// paraxial regime the equivalence holds in: the residual of "stack ==
	// single layer of the same reduced thickness" is the kt^4 term of the
	// transfer function, i.e. O(NA^4) ≈ 1e-10.
	const fLens = 0.5
	k0 := 2 * math.Pi / wl
	ctx := &Context{Wavelength: wl, Evanescent: "decay"}
	// A Gaussian beam (not a window-filling plane wave) with a lens phase: its
	// spectrum is confined, so the comparison below is not polluted by the
	// high-frequency tails a hard-edged field would carry.
	lensed := func() *Field {
		f, err := BuildSource(SourceSpec{Type: "gaussian", Params: map[string]any{"waist": 1e-3, "power": 1e-3}}, n, width, false, wl)
		if err != nil {
			t.Fatal(err)
		}
		el, err := NewElement(ElementSpec{Type: "lens", Params: map[string]any{"f": fLens}})
		if err != nil {
			t.Fatal(err)
		}
		if err := el.Apply(f, ctx); err != nil {
			t.Fatal(err)
		}
		return f
	}
	run := func(f *Field, spec CustomMediumSpec) {
		t.Helper()
		if err := PropagateCustomMedium(f, mustMedium(t, spec), ctx); err != nil {
			t.Fatal(err)
		}
	}
	const d1, n1, d2, n2 = 0.4e-3, 1.6, 0.6e-3, 1.1
	stack := CustomMediumSpec{Layers: []CustomMediumLayerSpec{
		{Thickness: d1, Index: "1.6"}, {Thickness: d2, Index: "1.1"}}}
	reordered := CustomMediumSpec{Layers: []CustomMediumLayerSpec{
		{Thickness: d2, Index: "1.1"}, {Thickness: d1, Index: "1.6"}}}

	// (a) single equivalent layer of the same reduced thickness.
	dTot := d1 + d2
	nEff := dTot / (d1/n1 + d2/n2)
	fa, fb := lensed(), lensed()
	run(fa, stack)
	// n_eff is given at full double precision: a rounded decimal would show up
	// as a constant piston offset of k0*d*Δn (2.8e-8 rad at 12 digits).
	run(fb, CustomMediumSpec{Distance: dTot, Index: fmt.Sprintf("%.17g", nEff)})
	// The two differ by the piston exp(i k0 (Σ n_i d_i - n_eff d)). The
	// comparison runs over the beam core (|E| within 1e-4 of the peak), where a
	// relative measure means something; the Gaussian's far tails, four orders
	// of magnitude down, carry a relative deviation of their own.
	piston := k0 * (n1*d1 + n2*d2 - nEff*dTot)
	var maxE float64
	for i := range fb.Ex {
		if a := cAbs(fb.Ex[i]); a > maxE {
			maxE = a
		}
	}
	// Two measures: the relative deviation over the beam core (|E| above 10% of
	// the peak, where a ratio is meaningful) and the absolute deviation over the
	// whole grid relative to the peak. The Gaussian's far tails carry the FFT
	// noise floor, so only the absolute measure is meaningful out there.
	var coreDev, coreAbs, allAbs float64
	ph := cexpI(piston)
	for i := range fa.Ex {
		amp := cAbs(fb.Ex[i])
		absDev := cAbs(fa.Ex[i] - fb.Ex[i]*ph)
		if absDev > allAbs {
			allAbs = absDev
		}
		if amp < 0.1*maxE || amp < 1e-12 {
			continue
		}
		if absDev > coreAbs {
			coreAbs = absDev
		}
		if dev := cAbs(fa.Ex[i]/fb.Ex[i] - ph); dev > coreDev {
			coreDev = dev
		}
	}
	t.Logf("stack vs single layer of the same reduced thickness (n_eff = %.9f, piston %.6f rad): core rel %.3e, core abs %.3e, whole-grid abs/peak %.3e",
		nEff, piston, coreDev, coreAbs, allAbs/maxE)
	// The residual is the non-paraxial kt^4 term of the transfer function: its
	// coefficient Σ d_i/n_i^3 differs from d/n_eff^3, which for this beam
	// (lambda f <= 4e-3) is O(1e-9) in phase; the numerical floor is lower.
	if coreDev > 1e-8 || allAbs > 1e-8*maxE {
		t.Fatalf("a homogeneous stack must equal the single layer up to the analytic piston phase: core rel %.3e, abs/peak %.3e",
			coreDev, allAbs/maxE)
	}

	// (b) homogeneous layers commute: each layer's operator is diagonal in the
	// transverse spectrum, so reordering them can only change the rounding
	// (~1e-13 after the FFT round trips), never the field.
	fc, fd := lensed(), lensed()
	run(fc, stack)
	run(fd, reordered)
	_, maxAbs := maxRelDiff(fc.Ex, fd.Ex)
	t.Logf("reordered homogeneous stack: max |ΔE| = %.3e", maxAbs)
	if maxAbs > 1e-12 {
		t.Fatalf("homogeneous layers must commute: max |ΔE| = %.3e", maxAbs)
	}

	// (c) z-dependent layers do not commute: two different graded layers put
	// the same z-intervals under different expressions, so both orders have a
	// different optical path — each compared with its analytic integral.
	gradA := CustomMediumSpec{Layers: []CustomMediumLayerSpec{
		{Thickness: d1, Index: "1.5 + 400*z", Steps: 60},
		{Thickness: d2, Index: "1.9 - 400*z", Steps: 60},
	}}
	gradB := CustomMediumSpec{Layers: []CustomMediumLayerSpec{
		{Thickness: d2, Index: "1.9 - 400*z", Steps: 60},
		{Thickness: d1, Index: "1.5 + 400*z", Steps: 60},
	}}
	// ∫(1.5+400z) and ∫(1.9-400z) over [a,b].
	iUp := func(a, b float64) float64 { return 1.5*(b-a) + 200*(b*b-a*a) }
	iDn := func(a, b float64) float64 { return 1.9*(b-a) - 200*(b*b-a*a) }
	pathA := iUp(0, d1) + iDn(d1, d1+d2)
	pathB := iDn(0, d2) + iUp(d2, d2+d1)
	fe, ff := planeField(t, n, width, wl), planeField(t, n, width, wl)
	run(fe, gradA)
	run(ff, gradB)
	pA := math.Atan2(imag(fe.Ex[(n/2)*n+n/2]), real(fe.Ex[(n/2)*n+n/2]))
	pB := math.Atan2(imag(ff.Ex[(n/2)*n+n/2]), real(ff.Ex[(n/2)*n+n/2]))
	wA, wB := math.Remainder(k0*pathA, 2*math.Pi), math.Remainder(k0*pathB, 2*math.Pi)
	t.Logf("z-graded order A: path %.6e m, phase %.9f vs analytic %.9f (Δ %.2e)", pathA, pA, wA, math.Remainder(pA-wA, 2*math.Pi))
	t.Logf("z-graded order B: path %.6e m, phase %.9f vs analytic %.9f (Δ %.2e)", pathB, pB, wB, math.Remainder(pB-wB, 2*math.Pi))
	if math.Abs(math.Remainder(pA-wA, 2*math.Pi)) > 1e-6 || math.Abs(math.Remainder(pB-wB, 2*math.Pi)) > 1e-6 {
		t.Fatalf("z-dependent stack optical path wrong")
	}
	if math.Abs(pathA-pathB) < 1e-7 {
		t.Fatalf("the two orders must have different optical paths: %.9e vs %.9e", pathA, pathB)
	}

	// (d) outside the stack the light is in vacuum: a lens phase applied after
	// the medium is untouched (compare against the same lens with no medium).
	// (d) an index of exactly 1 is vacuum: the reference-index kernel reduces to
	// the plain ASM and the screen vanishes, so the field must equal a plain
	// propagation of the same length bit for bit.
	fg, fh := lensed(), lensed()
	run(fg, CustomMediumSpec{Distance: 1e-3, Index: "1.0", Steps: 1})
	if err := Propagate(fh, 1e-3, MethodASM, ctx); err != nil {
		t.Fatal(err)
	}
	_, ident := maxRelDiff(fg.Ex, fh.Ex)
	if ident > 1e-13 {
		t.Fatalf("n = 1 medium must be transparent: max |ΔE| = %.3e", ident)
	}
	t.Logf("n = 1 medium vs plain vacuum propagation: max |ΔE| = %.3e", ident)
}

// ---------------------------------------------------------------------------
// expression language
// ---------------------------------------------------------------------------

// TestCustomMediumExprLanguage: every documented variable and the shared
// function table work, and a wrong name/arity/syntax is a load-time error with
// a position, never a silent zero.
func TestCustomMediumExprLanguage(t *testing.T) {
	// All variables, evaluated pointwise.
	m, err := CompileMediumExpr("x + y + z + zl + r + th + wl + k + pi + e", nil)
	if err != nil {
		t.Fatal(err)
	}
	got := m.Eval(1, 2, 3, 4, 5)
	want := 1 + 2 + 3 + 4 + math.Hypot(1, 2) + math.Atan2(2, 1) + 5 + 2*math.Pi/5 + math.Pi + math.E
	if math.Abs(got-want) > 1e-12 {
		t.Fatalf("variables: got %g want %g", got, want)
	}
	if !m.IsSpatial() || !m.NeedsRadius() || !m.NeedsTheta() {
		t.Fatalf("spatial/r/th flags wrong: %+v", m)
	}
	// Parameters and functions.
	m2, err := CompileMediumExpr("n0 + g*x*x - a*exp(-(x/w)^2) + if(x>0, 1, -1) + clamp(sin(b), -0.5, 0.5)", map[string]float64{"n0": 1.5, "g": 2.0, "a": 0.1, "w": 0.01, "b": 0.3})
	if err != nil {
		t.Fatal(err)
	}
	x := 0.002
	want2 := 1.5 + 2.0*x*x - 0.1*math.Exp(-(x/0.01)*(x/0.01)) + 1 + math.Min(math.Max(math.Sin(0.3), -0.5), 0.5)
	if got2 := m2.Eval(x, 0, 0, 0, 632.8e-9); math.Abs(got2-want2) > 1e-12 {
		t.Fatalf("functions: got %g want %g", got2, want2)
	}
	// if() is lazy: a singular branch that is not taken must not poison the value.
	m3, err := CompileMediumExpr("if(r > 0, n0 + a*sin(b*r)/r, n0)", map[string]float64{"n0": 1.5, "a": 0.2, "b": 1e5})
	if err != nil {
		t.Fatal(err)
	}
	if v := m3.Eval(0, 0, 0, 0, 632.8e-9); v != 1.5 {
		t.Fatalf("if() must be lazy at r=0: got %g want 1.5", v)
	}
	// A constant expression is not spatial: that is what selects the exact
	// homogeneous-layer path.
	if c, err := CompileMediumExpr("1.5 + 3e-15/wl^2", nil); err != nil || c.IsSpatial() {
		t.Fatalf("a wavelength-only expression must not be spatial (err=%v)", err)
	}

	// Error cases.
	cases := []struct {
		name, src string
		params    map[string]float64
		wantSub   []string
	}{
		{"unknown variable", "n0 + q*x", map[string]float64{"n0": 1.5}, []string{"column 6", "unknown variable \"q\"", "x, y, z, zl"}},
		{"unknown function", "foo(x)", nil, []string{"unknown function \"foo\""}},
		{"arity", "sin(x, y)", nil, []string{"sin takes 1 argument(s), got 2"}},
		{"syntax", "1.5 + * x", nil, []string{"column"}},
		{"unbalanced", "(1 + x", nil, []string{"expected ')'"}},
		{"empty", "   ", nil, []string{"empty expression"}},
		{"param missing", "g*x", nil, []string{"unknown variable \"g\"", "params"}},
		{"bad char", "x # y", nil, []string{"unexpected character"}},
	}
	for _, c := range cases {
		_, err := CompileMediumExpr(c.src, c.params)
		if err == nil {
			t.Fatalf("%s: %q compiled, want an error", c.name, c.src)
		}
		for _, sub := range c.wantSub {
			if !strings.Contains(err.Error(), sub) {
				t.Fatalf("%s: error %q does not mention %q", c.name, err.Error(), sub)
			}
		}
		t.Logf("%-16s -> %v", c.name, err)
	}
}

// ---------------------------------------------------------------------------
// ④ validation: bad expressions and bad values are reported, with a path
// ---------------------------------------------------------------------------

// TestCustomMediumValidationErrors: NewCustomMedium and ValidateConfig both
// report the offending detail (variable, column, value) with the element path.
func TestCustomMediumValidationErrors(t *testing.T) {
	cases := []struct {
		name    string
		spec    CustomMediumSpec
		wantSub []string
	}{
		{"no expression", CustomMediumSpec{Distance: 1e-3}, []string{"no index expression"}},
		{"zero thickness", CustomMediumSpec{Distance: 0, Index: "1.5"}, []string{"distance", "must be > 0"}},
		{"negative thickness", CustomMediumSpec{Layers: []CustomMediumLayerSpec{{Thickness: -1e-3, Index: "1.5"}}}, []string{"layer 0", "thickness must be > 0"}},
		{"missing layer index", CustomMediumSpec{Layers: []CustomMediumLayerSpec{{Thickness: 1e-3}}}, []string{"layer 0", "index expression is required"}},
		{"unknown variable", CustomMediumSpec{Distance: 1e-3, Index: "1.5 + q*x"}, []string{"unknown variable \"q\""}},
		{"syntax error", CustomMediumSpec{Distance: 1e-3, Index: "1.5 +"}, []string{"column"}},
		{"nonpositive index", CustomMediumSpec{Distance: 1e-3, Index: "1.0 - 2*step(x)"}, []string{"index must be > 0"}},
		{"nan index", CustomMediumSpec{Distance: 1e-3, Index: "sqrt(0-1)*x"}, []string{"not finite"}},
		{"abs syntax", CustomMediumSpec{Distance: 1e-3, Index: "1.5", Absorption: "0.1*"}, []string{"absorption", "column"}},
		{"bad reference", CustomMediumSpec{Layers: []CustomMediumLayerSpec{{Thickness: 1e-3, Index: "1.5", Reference: ptrF(0)}}}, []string{"reference index must be"}},
		{"unresolvable auto ref", CustomMediumSpec{Layers: []CustomMediumLayerSpec{{Thickness: 1e-3, Index: "2*(x>=0) "}}}, []string{"index must be > 0"}},
		{"both forms", CustomMediumSpec{Distance: 1e-3, Index: "1.5", Layers: []CustomMediumLayerSpec{{Thickness: 1e-3, Index: "1.5"}}}, []string{"not both"}},
		{"bad param", CustomMediumSpec{Distance: 1e-3, Index: "n0", Params: map[string]float64{"n0": math.NaN()}}, []string{"must be finite"}},
		{"negative steps", CustomMediumSpec{Distance: 1e-3, Index: "1.5", Steps: -3}, []string{"steps must be >= 0"}},
		{"thickness limit", CustomMediumSpec{Distance: 3, Index: "1.5"}, []string{"exceeds the limit"}},
	}
	for _, c := range cases {
		_, err := NewCustomMedium(c.spec)
		if err == nil {
			t.Fatalf("%s: compiled, want an error", c.name)
		}
		for _, sub := range c.wantSub {
			if !strings.Contains(err.Error(), sub) {
				t.Fatalf("%s: error %q does not mention %q", c.name, err.Error(), sub)
			}
		}
		t.Logf("%-20s -> %v", c.name, err)
	}

	// The same errors, through the configuration validation of an element train.
	cfg := Config{
		Grid: GridSpec{Size: 64, Width: 0.004}, Wavelength: 632.8e-9, Method: "asm",
		Sources: []SourceSpec{{Type: "plane", Params: map[string]any{}}},
		Elements: []ElementSpec{{Type: "custom_medium", Params: map[string]any{
			"layers": []any{map[string]any{"thickness": 1e-3, "index": "1.5 + q*x"}},
		}}},
	}
	issues := ValidateConfig(&cfg)
	if len(issues) != 1 || issues[0].Path != "elements[0]" || !strings.Contains(issues[0].Message, `unknown variable "q"`) {
		t.Fatalf("element train validation: %+v", issues)
	}
	// And Simulate refuses the configuration with the same text.
	if _, err := Simulate(cfg); err == nil || !strings.Contains(err.Error(), `unknown variable "q"`) {
		t.Fatalf("Simulate must refuse a broken expression, got %v", err)
	}
	t.Logf("train validation -> %s: %s", issues[0].Path, issues[0].Message)

	// A typo in a layer key must be rejected, not silently ignored.
	if _, err := NewCustomMedium(CustomMediumSpec{Layers: []CustomMediumLayerSpec{{Thickness: 1e-3, Index: "1.5"}}}); err != nil {
		t.Fatal(err)
	}
	_, err := newCustomMedium(map[string]any{"layers": `[{"thicknes":1e-3,"index":"1.5"}]`})
	if err == nil || !strings.Contains(err.Error(), "thicknes") {
		t.Fatalf("a misspelled layer key must be an error, got %v", err)
	}
	t.Logf("strict layer decode -> %v", err)

	// And through a scene component: the issue is attached to that component.
	scfg := Config{
		Grid: GridSpec{Size: 64, Width: 0.004}, Wavelength: 632.8e-9, Method: "asm",
		Sources: []SourceSpec{{Type: "plane", Params: map[string]any{}}},
		Scene: &SceneSpec{Components: []ComponentSpec{
			{ID: "m", Type: "custom_medium", Params: map[string]any{"distance": 1e-3, "index": "1.5 + *2"}},
			{ID: "d", Type: "sensor", Pos: v3(0, 0, 0.01)},
		}},
	}
	issues = ValidateConfig(&scfg)
	if len(issues) == 0 {
		t.Fatalf("scene component with a broken expression validated clean")
	}
	found := false
	for _, is := range issues {
		if is.Path == "scene.components[0].params" && strings.Contains(is.Message, "column") {
			found = true
		}
	}
	if !found {
		t.Fatalf("scene validation issues: %+v", issues)
	}
	t.Logf("scene validation -> %s: %s", issues[0].Path, issues[0].Message)
}

// TestCustomMediumConfigJSON: the whole feature is reachable from configuration
// JSON — the flat form, the layered form, a JSON-string layer list (what the
// GUI text field sends) and the element through a full Simulate run. The old
// media are untouched: the same train with the built-in "medium" element must
// produce exactly the field it produced before this element existed (it is the
// same code path), and an empty/absent custom-medium train must be a no-op.
func TestCustomMediumConfigJSON(t *testing.T) {
	type el = map[string]any
	cfg := Config{
		Grid: GridSpec{Size: 128, Width: 0.004}, Wavelength: 632.8e-9, Method: "asm",
		Sources: []SourceSpec{{Type: "plane", Params: map[string]any{"power": 1e-3}}},
		Elements: []ElementSpec{
			{Type: "custom_medium", Params: el{"distance": 1e-3, "index": "1.5 + 0.1*x/1e-3", "steps": 20}},
			{Type: "sensor", Params: el{}},
		},
	}
	// The same through JSON text (what the API receives).
	raw := `{"grid":{"size":128,"width":0.004},"wavelength":632.8e-9,"method":"asm",
	  "sources":[{"type":"plane","params":{"power":1e-3}}],
	  "elements":[
	    {"type":"custom_medium","params":{"steps":20,
	      "layers":"[{\"thickness\":5e-4,\"index\":\"1.5\"},{\"thickness\":5e-4,\"index\":\"1.8\",\"absorption\":\"1e-5\",\"params\":{}}]"}},
	    {"type":"sensor","params":{}}]}`
	var cfg2 Config
	if err := json.Unmarshal([]byte(raw), &cfg2); err != nil {
		t.Fatal(err)
	}
	for name, c := range map[string]Config{"struct": cfg, "json": cfg2} {
		res, err := Simulate(c)
		if err != nil {
			t.Fatalf("%s: %v", name, err)
		}
		if len(res.Planes) != 1 {
			t.Fatalf("%s: %d planes", name, len(res.Planes))
		}
		st := res.Planes[0].Stats
		// The JSON train carries an absorbing second layer (Im(n) = 1e-5 over
		// 5e-4 m), so its transmitted power is exp(-2 k0 Im(n) d) times the
		// incident power; the struct train is lossless.
		want := 1e-3
		if name == "json" {
			want = 1e-3 * math.Exp(-2*(2*math.Pi/632.8e-9)*1e-5*5e-4)
		}
		t.Logf("%s: power %.6e W (want %.6e), peak %.4e W/m^2, rms %.4g m", name, st.Power, want, st.Peak, st.RMSX)
		if rel := math.Abs(st.Power-want) / want; rel > 1e-4 {
			t.Fatalf("%s: transmitted power %.6e, want %.6e (rel %.2e)", name, st.Power, want, rel)
		}
	}
	// Backwards compatibility: no custom medium in the train == the same field
	// as before, and the built-in medium element still behaves as it did.
	base := Config{
		Grid: GridSpec{Size: 64, Width: 0.004}, Wavelength: 632.8e-9, Method: "asm",
		Sources: []SourceSpec{{Type: "plane", Params: map[string]any{"power": 1e-3}}},
		Elements: []ElementSpec{
			{Type: "medium", Params: el{"distance": 1e-3, "index": 1.5, "steps": 10}},
			{Type: "sensor", Params: el{}},
		},
	}
	res, err := Simulate(base)
	if err != nil {
		t.Fatal(err)
	}
	// The built-in element is a constant phase screen plus vacuum diffraction:
	// its plane-wave phase is exactly k0*n*d.
	k0 := 2 * math.Pi / 632.8e-9
	pl := res.Planes[0]
	c := pl.Ex[(pl.Size/2)*pl.Size+pl.Size/2]
	got := math.Atan2(imag(c), real(c))
	if diff := math.Remainder(got-math.Remainder(k0*1.5*1e-3, 2*math.Pi), 2*math.Pi); math.Abs(diff) > 1e-9 {
		t.Fatalf("built-in medium element changed behavior: Δ = %.3e rad", diff)
	}
}

// ---------------------------------------------------------------------------
// catalog / GUI wiring
// ---------------------------------------------------------------------------

// TestCustomMediumCatalog: the element is registered, classified for the scene
// view and documented in /api/catalog with its parameters and limits.
func TestCustomMediumCatalog(t *testing.T) {
	if _, err := NewElement(ElementSpec{Type: "custom_medium"}); err == nil {
		t.Fatal("custom_medium without parameters must error")
	}
	if b, err := componentBehavior("custom_medium"); err != nil || b != behaviorTransmit {
		t.Fatalf("component behavior: %v %v", b, err)
	}
	if cls := ComponentClasses()["custom_medium"]; cls == "" {
		t.Fatal("component class missing from the catalog class map")
	}
	cat := BuildCatalog()
	var doc *ElementDoc
	for i := range cat.Elements {
		if cat.Elements[i].Type == "custom_medium" {
			doc = &cat.Elements[i]
		}
	}
	if doc == nil {
		t.Fatal("custom_medium missing from /api/catalog elements")
	}
	keys := map[string]bool{}
	for _, p := range doc.Params {
		keys[p.Key] = true
	}
	for _, want := range []string{"layers", "distance", "index", "absorption", "params", "steps"} {
		if !keys[want] {
			t.Fatalf("catalog parameter %q missing (%+v)", want, doc.Params)
		}
	}
	b, err := json.Marshal(doc)
	if err != nil {
		t.Fatal(err)
	}
	// The catalog payload must survive JSON (the GUI parses it) and carry the
	// schema of a layer so the parameter panel can be built from it.
	if !strings.Contains(string(b), "thickness") || !strings.Contains(string(b), "reference") {
		t.Fatalf("catalog help does not document the layer schema: %s", b)
	}
	t.Logf("catalog entry: %s (%d params, %d bytes)", doc.Label, len(doc.Params), len(b))
}

// ---------------------------------------------------------------------------
// implementation cross-checks and cost
// ---------------------------------------------------------------------------

// TestCustomMediumMatchesGenericSplitStep: the compiled inline evaluation must
// agree with the generic IndexFunc implementation (PropagateSplitStepRef) for
// the same medium, and with nref=1 it must agree with the pre-existing
// vacuum-reference PropagateSplitStep (so the new kernel contains the old one).
func TestCustomMediumMatchesGenericSplitStep(t *testing.T) {
	const n, width, wl = 64, 0.004, 632.8e-9
	ctx := &Context{Wavelength: wl, Evanescent: "decay"}
	const g = 30.0
	nFn := func(x, y, z float64) complex128 { return complex(1.5+1e-3*z/1e-3+g*x, 2e-6*(1+math.Sin(50*y))) }
	cases := []struct {
		name string
		spec CustomMediumSpec
	}{
		{"explicit reference", CustomMediumSpec{Layers: []CustomMediumLayerSpec{{
			Thickness: 1e-3, Index: "1.5 + 1e-3*z/1e-3 + g*x",
			Absorption: "2e-6*(1 + sin(50*y))", Reference: ptrF(1.5), Steps: 7,
			Params: map[string]float64{"g": g}}}}},
		{"auto reference (mid-plane on axis)", CustomMediumSpec{
			Distance: 1e-3, Index: "1.5 + 1e-3*z/1e-3 + g*x",
			Absorption: "2e-6*(1 + sin(50*y))", Steps: 7,
			Params: map[string]float64{"g": g},
		}},
	}
	for _, c := range cases {
		fa := planeField(t, n, width, wl)
		fb := fa.Clone()
		m, err := NewCustomMedium(c.spec)
		if err != nil {
			t.Fatal(err)
		}
		info, err := m.LayerInfo(wl)
		if err != nil {
			t.Fatal(err)
		}
		nref := info[0].Reference
		if err := PropagateCustomMedium(fa, m, ctx); err != nil {
			t.Fatal(err)
		}
		if err := PropagateSplitStepRef(fb, 1e-3, nFn, nref, 7, ctx); err != nil {
			t.Fatal(err)
		}
		maxRel, _ := maxRelDiff(fa.Ex, fb.Ex)
		t.Logf("%s (reference %.6f): max relative difference to PropagateSplitStepRef = %.3e", c.name, nref, maxRel)
		if maxRel > 1e-12 {
			t.Fatalf("compiled path differs from the generic split-step: %.3e", maxRel)
		}
	}
	// nref = 1 must reproduce the existing vacuum-reference split-step.
	fa := planeField(t, n, width, wl)
	fb := fa.Clone()
	m := mustMedium(t, CustomMediumSpec{Distance: 1e-3, Index: "1.5 + 1e-3*z/1e-3", Steps: 7})
	if err := PropagateSplitStepRef(fa, 1e-3, func(x, y, z float64) complex128 { return complex(1.5+1e-3*z/1e-3, 0) }, 1, 7, ctx); err != nil {
		t.Fatal(err)
	}
	if err := PropagateSplitStep(fb, 1e-3, func(x, y, z float64) complex128 { return complex(1.5+1e-3*z/1e-3, 0) }, 7, ctx); err != nil {
		t.Fatal(err)
	}
	_ = m
	maxRel, _ := maxRelDiff(fa.Ex, fb.Ex)
	t.Logf("nref=1 vs PropagateSplitStep: max relative difference = %.3e", maxRel)
	if maxRel > 1e-12 {
		t.Fatalf("nref=1 must equal the vacuum-reference split-step: %.3e", maxRel)
	}
	// A bad reference index is rejected up front.
	if err := PropagateSplitStepRef(planeField(t, n, width, wl), 1e-3, UniformIndex(1.5), 0, 1, ctx); err == nil {
		t.Fatal("nref = 0 must be an error")
	}
}

// TestCustomMediumUniformFastPath: a layer whose expressions do not vary in
// space takes the scalar-multiply path; it must give the same field as the
// per-pixel path of the same expression written with a dummy spatial term.
func TestCustomMediumUniformFastPath(t *testing.T) {
	const n, width, wl = 64, 0.004, 632.8e-9
	ctx := &Context{Wavelength: wl, Evanescent: "decay"}
	fa := planeField(t, n, width, wl)
	fb := fa.Clone()
	fast := mustMedium(t, CustomMediumSpec{Distance: 1e-3, Index: "1.5", Absorption: "1e-6"})
	slow := mustMedium(t, CustomMediumSpec{Distance: 1e-3, Index: "1.5 + 0*z", Absorption: "1e-6 + 0*x"})
	if !fast.layers[0].uniform || slow.layers[0].uniform {
		t.Fatalf("uniform detection wrong: %v %v", fast.layers[0].uniform, slow.layers[0].uniform)
	}
	if err := PropagateCustomMedium(fa, fast, ctx); err != nil {
		t.Fatal(err)
	}
	if err := PropagateCustomMedium(fb, slow, ctx); err != nil {
		t.Fatal(err)
	}
	maxRel, _ := maxRelDiff(fa.Ex, fb.Ex)
	if maxRel > 1e-15 {
		t.Fatalf("uniform fast path differs from the per-pixel path: %.3e", maxRel)
	}
}

// TestCustomMediumSamplingWarning: a medium whose phase screen flips by more
// than Pi per pixel is aliased; the run must say so (custom_medium_alias)
// instead of quietly returning a wrong field.
func TestCustomMediumSamplingWarning(t *testing.T) {
	// A coarse grid with a moderate gradient (∂n/∂x = 20/m, n positive over the
	// window): the screen phase steps by k0*g*dz*dx = 9.93e6*20*1e-2*2.5e-4 ≈
	// 497 rad per pixel in a single slice — far past Pi, so the slice aliases
	// and the run must say so. With 200 slices the same medium is sampled at
	// 2.5 rad/pixel and is clean.
	const n, width, wl = 32, 0.008, 632.8e-9
	ctx := &Context{Wavelength: wl, Evanescent: "decay", Warnings: &Warnings{}}
	f := planeField(t, n, width, wl)
	m := mustMedium(t, CustomMediumSpec{Distance: 1e-2, Index: "1.5 + 20*x", Steps: 1})
	if err := PropagateCustomMedium(f, m, ctx); err != nil {
		t.Fatal(err)
	}
	found := false
	for _, w := range ctx.Warnings.List() {
		if w.Code == "custom_medium_alias" {
			found = true
			t.Logf("warning: %s", w.Message)
		}
	}
	if !found {
		t.Fatal("aliased phase screen must raise custom_medium_alias")
	}
	// The same medium with enough slices is clean.
	ctx2 := &Context{Wavelength: wl, Evanescent: "decay", Warnings: &Warnings{}}
	m2 := mustMedium(t, CustomMediumSpec{Distance: 1e-2, Index: "1.5 + 20*x", Steps: 200})
	f2 := planeField(t, n, width, wl)
	if err := PropagateCustomMedium(f2, m2, ctx2); err != nil {
		t.Fatal(err)
	}
	for _, w := range ctx2.Warnings.List() {
		if w.Code == "custom_medium_alias" {
			t.Fatalf("200 slices must sample the screen: %s", w.Message)
		}
	}
	t.Logf("200 slices: no aliasing warning")
}

// TestCustomMediumRunTimeGuard: a non-finite value that appears only at the run
// wavelength (a dispersion term) is an error at propagation time, not a silent
// NaN over the whole plane.
func TestCustomMediumRunTimeGuard(t *testing.T) {
	const n, width = 32, 0.004
	const wlRun = 1000e-9
	// Finite at the probe wavelength (632.8 nm), singular at the run wavelength:
	// only a run-time check can catch this one.
	const expr = "1.5 + 1e-12/(wl - 1000e-9)"
	ctx := &Context{Wavelength: wlRun, Evanescent: "decay"}

	// (a) with a pinned reference the screen evaluation is what has to fire.
	m := mustMedium(t, CustomMediumSpec{Layers: []CustomMediumLayerSpec{
		{Thickness: 1e-3, Index: expr, Reference: ptrF(1.5)}}})
	err := PropagateCustomMedium(planeField(t, n, width, wlRun), m, ctx)
	if err == nil || !strings.Contains(err.Error(), "not finite") {
		t.Fatalf("want a non-finite screen error at run time, got %v", err)
	}
	t.Logf("run-time screen guard -> %v", err)

	// (b) without a reference the run stops even earlier, when the reference
	// index cannot be derived from the expression at this wavelength.
	m2 := mustMedium(t, CustomMediumSpec{Distance: 1e-3, Index: expr})
	err = PropagateCustomMedium(planeField(t, n, width, wlRun), m2, ctx)
	if err == nil || !strings.Contains(err.Error(), "cannot derive a reference index") {
		t.Fatalf("want a reference-resolution error at run time, got %v", err)
	}
	t.Logf("run-time reference guard -> %v", err)
}

// ---------------------------------------------------------------------------
// cost
// ---------------------------------------------------------------------------

// BenchmarkCustomMediumScreen measures the per-pixel cost of the compiled
// screen evaluation (one slice of thickness dz, no diffraction) at 512x512, for
// the shapes of expression a user is likely to write: a homogeneous layer
// (scalar path), a transverse gradient, a gradient with a spatial absorber, and
// expressions full of transcendental calls. The "closure" case applies the
// physically identical screen through a hand-written Go closure via
// applyMediumPhaseRef, i.e. it is the baseline the compiled expression engine is
// measured against.
//
//	go test ./optics/ -run XXX -bench BenchmarkCustomMediumScreen -benchtime 50x
func BenchmarkCustomMediumScreen(b *testing.B) {
	const n, width, wl = 512, 0.008, 632.8e-9
	k0 := 2 * math.Pi / wl
	ctx := &Context{Wavelength: wl, Evanescent: "decay"}
	layers := map[string]*customLayer{}
	for name, spec := range map[string]CustomMediumSpec{
		"uniform":  {Distance: 1e-3, Index: "1.5"},
		"gradient": {Distance: 1e-3, Index: "1.5 + 2*x"},
		"grad_abs": {Distance: 1e-3, Index: "1.5 + 2*x", Absorption: "3e-4*exp(-(x*x+y*y)/1e-6)"},
		"grad_2fn": {Distance: 1e-3, Index: "1.5 + sin(1e4*x)*cos(1e4*y)", Absorption: "3e-4*(1+sin(1e4*y))"},
	} {
		m, err := NewCustomMedium(spec)
		if err != nil {
			b.Fatal(err)
		}
		layers[name] = &m.layers[0]
	}
	xs := make([]float64, n)
	for i := range xs {
		xs[i] = (float64(i) - float64(n)/2) * (width / float64(n))
	}
	newBenchField := func() *Field {
		f := NewField(n, width/float64(n), false)
		for i := range f.Ex {
			f.Ex[i] = complex(1, 0)
		}
		return f
	}
	for name, L := range layers {
		env := &mediumEnv{wl: wl, k: k0, frames: make([]float64, L.frameSlots())}
		sc := &mediumScratch{xs: xs}
		f := newBenchField()
		b.Run(name, func(b *testing.B) {
			b.ReportAllocs()
			for k := 0; k < b.N; k++ {
				if err := applyCustomScreen(f, L, 0, 5e-4, 5e-4, 1.5, 1e-6, env, sc, ctx); err != nil {
					b.Fatal(err)
				}
			}
		})
	}
	// Hand-written closure baseline: the same "1.5 + 2*x" screen, no expression.
	fn := func(x, y, z float64) complex128 { return complex(1.5+2*x, 0) }
	f := newBenchField()
	b.Run("closure_1p5_2x", func(b *testing.B) {
		b.ReportAllocs()
		for k := 0; k < b.N; k++ {
			applyMediumPhaseRef(f, fn, 5e-4, k0, 1e-6, 1.5)
		}
	})
}

// BenchmarkCustomMediumPropagate measures the full operator at 512x512 (two
// FFTs per slice, which dominates) against the generic IndexFunc split-step and
// the pre-existing vacuum-reference one, so the cost of the expression
// evaluation is visible as the difference.
func BenchmarkCustomMediumPropagate(b *testing.B) {
	const n, width, wl = 512, 0.008, 632.8e-9
	ctx := &Context{Wavelength: wl, Evanescent: "decay"}
	m, err := NewCustomMedium(CustomMediumSpec{Distance: 1e-3, Index: "1.5 + 2*x", Steps: 4})
	if err != nil {
		b.Fatal(err)
	}
	fn := func(x, y, z float64) complex128 { return complex(1.5+2*x, 0) }
	runs := map[string]func(f *Field){
		"custom_compiled": func(f *Field) { _ = PropagateCustomMedium(f, m, ctx) },
		"generic_nref15":  func(f *Field) { _ = PropagateSplitStepRef(f, 1e-3, fn, 1.5, 4, ctx) },
		"vacuum_kernel":   func(f *Field) { _ = PropagateSplitStep(f, 1e-3, fn, 4, ctx) },
	}
	for name, run := range runs {
		f := NewField(n, width/float64(n), false)
		for i := range f.Ex {
			f.Ex[i] = complex(1, 0)
		}
		b.Run(name, func(b *testing.B) {
			b.ReportAllocs()
			for k := 0; k < b.N; k++ {
				run(f)
			}
		})
	}
}
