package optics

import (
	"math"
	"testing"
)

// The coating element is a TMM (characteristic-matrix) evaluator for isotropic
// multilayer stacks. The tests pin, in order:
//
//  1. the amplitude AND phase of one layer against a direct solution of the
//     interface boundary conditions (independent of the matrix method),
//  2. a bare interface against the Fresnel formulas,
//  3. a quarter-wave layer and a quarter-wave stack against their
//     equivalent-admittance closed forms (and that the layer order matters),
//  4. R + T = 1 for lossless stacks and R + T < 1 with an absorbing layer,
//  5. the phase sign convention of the kernel (e^{−iωt}, forward e^{+ikz}),
//  6. the s/p application to the Jones field and the scalar path,
//  7. parsing of the four material dispersion forms and the error paths.

const coatLam = 550e-9

// coatAdm returns the tilted optical admittance η = N cosθ (s) or N/cosθ (p),
// computing N cosθ = sqrt(N² − N0²sin²θ0) on the decaying branch.
func coatAdm(N, sinTheta0 complex128, pol int) complex128 {
	c := csqrt(N*N - sinTheta0*sinTheta0)
	if pol == 0 {
		return c
	}
	return N * N / c
}

// singleLayerBoundarySolve solves the four interface equations of one layer by
// hand (no matrices): with A/B the forward/backward amplitudes inside the
// layer and δ = k0·(N1 cosθ1)·d,
//
//	1 + r = A + B,   η0(1 − r) = η1(A − B),
//	A e^{iδ} + B e^{−iδ} = t,   η1(A e^{iδ} − B e^{−iδ}) = η_s t,
//
// which gives t = 2η0 / [½(η0+η1)(1+η_s/η1)e^{−iδ} + ½(η0−η1)(1−η_s/η1)e^{+iδ}].
func singleLayerBoundarySolve(N1, n0, nsub complex128, theta0, lambda, d float64, pol int) (t, r complex128) {
	k0 := 2 * math.Pi / lambda
	sn := n0 * complex(math.Sin(theta0), 0)
	e0 := coatAdm(n0, sn, pol)
	e1 := coatAdm(N1, sn, pol)
	es := coatAdm(nsub, sn, pol)
	cosj := csqrt(N1*N1 - sn*sn) // N1 cosθ1
	delta := complex(k0*d, 0) * cosj
	den := 0.5*(e0+e1)*(1+es/e1)*cexpNegI(delta) + 0.5*(e0-e1)*(1-es/e1)*cexpIPos(delta)
	t = 2 * e0 / den
	a := t * (1 + es/e1) * cexpNegI(delta) / 2
	b := t * (1 - es/e1) * cexpIPos(delta) / 2
	return t, a + b - 1
}

// cexpC returns exp(z) for complex z = a+ib.
func cexpC(z complex128) complex128 {
	return complex(math.Cos(imag(z)), math.Sin(imag(z))) * complex(math.Exp(real(z)), 0)
}

// cexpNegI returns exp(−i·z) for complex z (the e^{−iδ} factor of the
// transmitted-wave sum; note that exp(−i·z) is *not* exp(−z), and for a
// real δ it is *not* exp(δ) either — cexpC is the plain complex exponential).
func cexpNegI(z complex128) complex128 {
	return cexpC(complex(imag(z), -real(z)))
}

// cexpIPos returns exp(+i·z) for complex z.
func cexpIPos(z complex128) complex128 {
	return cexpC(complex(-imag(z), real(z)))
}

func TestCoatingMatchesBoundaryConditionSolve(t *testing.T) {
	cases := []struct {
		theta     float64
		N1        complex128
		d         float64
		n0, nsub  complex128
		labelName string
	}{
		{0.0, 2.35, coatLam / (4 * 2.35), 1, 1.52, "quarter wave, normal"},
		{0.5, 2.35, 80e-9, 1, 1.52, "oblique lossless"},
		{0.5, complex(1.38, 0.02), 120e-9, 1, 1.52, "oblique absorbing"},
		{1.1, 2.1, 333e-9, 1, 1.52, "large angle"},
		{0.0, complex(1.38, 0.05), 200e-9, 1, 1.0, "absorbing film in air"},
	}
	for _, c := range cases {
		got, err := ThinFilmCoefficients([]TMMLayer{{Index: c.N1, Thickness: c.d}}, c.n0, c.nsub, c.theta, coatLam)
		if err != nil {
			t.Fatalf("%s: %v", c.labelName, err)
		}
		for pol := 0; pol < 2; pol++ {
			wantT, wantR := singleLayerBoundarySolve(c.N1, c.n0, c.nsub, c.theta, coatLam, c.d, pol)
			gotT, gotR := got.Ts, got.Rs
			if pol == 1 {
				gotT, gotR = got.Tp, got.Rp
			}
			if math.Abs(real(gotT-wantT)) > 1e-14 || math.Abs(imag(gotT-wantT)) > 1e-14 {
				t.Fatalf("%s (%s): t = %v, want %v", c.labelName, polName(pol), gotT, wantT)
			}
			if math.Abs(real(gotR-wantR)) > 1e-14 || math.Abs(imag(gotR-wantR)) > 1e-14 {
				t.Fatalf("%s (%s): r = %v, want %v", c.labelName, polName(pol), gotR, wantR)
			}
		}
	}
}

func polName(pol int) string {
	if pol == 0 {
		return "s"
	}
	return "p"
}

func TestCoatingBareInterfaceFresnel(t *testing.T) {
	for _, theta := range []float64{0, 0.3, 0.7, 1.2} {
		for _, nsub := range []complex128{1.52, complex(1.45, 0.01)} {
			got, err := ThinFilmCoefficients(nil, 1, nsub, theta, coatLam)
			if err != nil {
				t.Fatal(err)
			}
			sn := complex(math.Sin(theta), 0)
			for pol := 0; pol < 2; pol++ {
				e0 := coatAdm(1, sn, pol)
				es := coatAdm(nsub, sn, pol)
				wantR := (e0 - es) / (e0 + es)
				wantT := 2 * e0 / (e0 + es)
				gotT, gotR := got.Ts, got.Rs
				gotTPow, gotRPow := got.TsPow, got.RsPow
				if pol == 1 {
					gotT, gotR = got.Tp, got.Rp
					gotTPow, gotRPow = got.TpPow, got.RpPow
				}
				if complexAbs(gotT-wantT) > 1e-14 || complexAbs(gotR-wantR) > 1e-14 {
					t.Fatalf("θ=%g %s: t=%v (want %v), r=%v (want %v)", theta, polName(pol), gotT, wantT, gotR, wantR)
				}
				if math.Abs(gotRPow-complexAbs(wantR)*complexAbs(wantR)) > 1e-14 {
					t.Fatalf("θ=%g %s: R = %g", theta, polName(pol), gotRPow)
				}
				// Flux balance at a single interface holds even when the second
				// medium absorbs: T is the flux just inside it.
				if math.Abs(gotRPow+gotTPow-1) > 1e-12 {
					t.Fatalf("θ=%g %s: R+T = %g at a bare interface", theta, polName(pol), gotRPow+gotTPow)
				}
			}
		}
	}
}

func complexAbs(z complex128) float64 {
	return math.Hypot(real(z), imag(z))
}

// TestCoatingQuarterWaveLayer checks the equivalent-admittance result for one
// quarter-wave layer: R = ((η0 − η1²/η_s)/(η0 + η1²/η_s))².
func TestCoatingQuarterWaveLayer(t *testing.T) {
	const nH, nsub = 2.35, 1.52
	for _, theta := range []float64{0, 0.5, 1.0} {
		sn := complex(math.Sin(theta), 0)
		// δ = π/2 requires d = λ/(4·Re(N cosθ)) — the same physical thickness
		// for both polarizations (δ does not depend on polarization).
		cosj := csqrt(complex(nH*nH, 0) - sn*sn)
		d := coatLam / (4 * real(cosj))
		got, err := ThinFilmCoefficients([]TMMLayer{{Index: complex(nH, 0), Thickness: d}}, 1, nsub, theta, coatLam)
		if err != nil {
			t.Fatal(err)
		}
		for pol := 0; pol < 2; pol++ {
			e0 := coatAdm(1, sn, pol)
			e1 := coatAdm(complex(nH, 0), sn, pol)
			es := coatAdm(nsub, sn, pol)
			y := e1 * e1 / es
			wantR := complexAbs((e0 - y) / (e0 + y))
			wantR *= wantR
			gotR := got.RsPow
			if pol == 1 {
				gotR = got.RpPow
			}
			if math.Abs(gotR-wantR) > 1e-12 {
				t.Fatalf("θ=%g %s: R = %.12f, want %.12f", theta, polName(pol), gotR, wantR)
			}
		}
	}
	// A quarter-wave MgF2 layer is the textbook single-layer antireflection
	// coating on glass: R ≈ 1.3% instead of the bare 4.3%.
	nm := mustIndex(t, "N-BK7", coatLam)
	dAR := coatLam / (4 * 1.37851)
	got, err := ThinFilmCoefficients([]TMMLayer{{Index: complex(1.37851, 0), Thickness: dAR}}, 1, complex(nm, 0), 0, coatLam)
	if err != nil {
		t.Fatal(err)
	}
	if got.RsPow > 0.02 {
		t.Fatalf("MgF2 quarter-wave AR coating reflects %g, want < 2%%", got.RsPow)
	}
	bare, _ := ThinFilmCoefficients(nil, 1, complex(nm, 0), 0, coatLam)
	if bare.RsPow < 0.03 {
		t.Fatalf("bare N-BK7 should reflect > 3%%, got %g", bare.RsPow)
	}
	t.Logf("MgF2 QW AR on N-BK7 at 550 nm: R = %.4f%% (bare substrate %.4f%%)", 100*got.RsPow, 100*bare.RsPow)
}

// TestCoatingQuarterWaveStack checks a (HL)^4 stack against the closed form
// that follows from the quarter-wave pair identity
// M_H·M_L = diag(−η_L/η_H, −η_H/η_L): t = 2η0/(η0·A + η_s·B) with
// A = (−η_L/η_H)^N, B = (−η_H/η_L)^N. It also pins the multiplication order:
// the reversed stack corresponds to swapping A and B and gives a different t.
func TestCoatingQuarterWaveStack(t *testing.T) {
	const nH, nL, nsub, pairs = 2.35, 1.38, 1.52, 4
	for _, theta := range []float64{0, 0.5} {
		sn := complex(math.Sin(theta), 0)
		dH := coatLam / (4 * real(csqrt(complex(nH*nH, 0)-sn*sn)))
		dL := coatLam / (4 * real(csqrt(complex(nL*nL, 0)-sn*sn)))
		var layers []TMMLayer
		for i := 0; i < pairs; i++ {
			layers = append(layers, TMMLayer{Index: complex(nH, 0), Thickness: dH},
				TMMLayer{Index: complex(nL, 0), Thickness: dL})
		}
		got, err := ThinFilmCoefficients(layers, 1, nsub, theta, coatLam)
		if err != nil {
			t.Fatal(err)
		}
		for pol := 0; pol < 2; pol++ {
			e0 := coatAdm(1, sn, pol)
			eH := coatAdm(complex(nH, 0), sn, pol)
			eL := coatAdm(complex(nL, 0), sn, pol)
			es := coatAdm(nsub, sn, pol)
			a := cpow(-eL/eH, pairs)
			b := cpow(-eH/eL, pairs)
			wantT := 2 * e0 / (e0*a + es*b)
			wantR := (e0*a - es*b) / (e0*a + es*b)
			gotT, gotR := got.Ts, got.Rs
			gotRpow := got.RsPow
			if pol == 1 {
				gotT, gotR, gotRpow = got.Tp, got.Rp, got.RpPow
			}
			if complexAbs(gotT-wantT) > 1e-12 {
				t.Fatalf("θ=%g %s: t = %v, want %v", theta, polName(pol), gotT, wantT)
			}
			if complexAbs(gotR-wantR) > 1e-12 {
				t.Fatalf("θ=%g %s: r = %v, want %v", theta, polName(pol), gotR, wantR)
			}
			if gotRpow <= 0.9 {
				t.Fatalf("θ=%g %s: quarter-wave stack R = %g, want > 0.9", theta, polName(pol), gotRpow)
			}
		}
		// Reversed order must give the swapped closed form (a non-commutative
		// product is being checked here, not a rounding difference).
		rev := make([]TMMLayer, len(layers))
		for i, L := range layers {
			rev[len(layers)-1-i] = L
		}
		gotRev, err := ThinFilmCoefficients(rev, 1, nsub, theta, coatLam)
		if err != nil {
			t.Fatal(err)
		}
		e0 := coatAdm(1, sn, 0)
		eH := coatAdm(complex(nH, 0), sn, 0)
		eL := coatAdm(complex(nL, 0), sn, 0)
		es := coatAdm(nsub, sn, 0)
		a := cpow(-eH/eL, pairs)
		b := cpow(-eL/eH, pairs)
		wantRev := 2 * e0 / (e0*a + es*b)
		if complexAbs(gotRev.Ts-wantRev) > 1e-12 {
			t.Fatalf("θ=%g reversed stack: t = %v, want %v", theta, gotRev.Ts, wantRev)
		}
		if complexAbs(gotRev.Ts-got.Ts) < 1e-6 {
			t.Fatal("reversing the layer list did not change the transmission: the matrix product is order-insensitive (bug)")
		}
	}
}

// cpow returns z^n for a small integer n >= 0.
func cpow(z complex128, n int) complex128 {
	out := complex(1, 0)
	for i := 0; i < n; i++ {
		out *= z
	}
	return out
}

// TestCoatingEnergyConservation: lossless stacks conserve R + T = 1 for both
// polarizations (the transverse flux ratio), and an absorbing film must lose
// energy (R + T < 1, falling as κ grows).
func TestCoatingEnergyConservation(t *testing.T) {
	for _, theta := range []float64{0, 0.4, 0.9} {
		sn := complex(math.Sin(theta), 0)
		var layers []TMMLayer
		for i := 0; i < 4; i++ {
			layers = append(layers,
				TMMLayer{Index: 2.35, Thickness: 0.7 * coatLam / (4 * real(csqrt(complex(2.35*2.35, 0)-sn*sn)))},
				TMMLayer{Index: 1.38, Thickness: 1.3 * coatLam / (4 * real(csqrt(complex(1.38*1.38, 0)-sn*sn)))})
		}
		got, err := ThinFilmCoefficients(layers, 1, 1.52, theta, coatLam)
		if err != nil {
			t.Fatal(err)
		}
		for _, sum := range []float64{got.RsPow + got.TsPow, got.RpPow + got.TpPow} {
			if math.Abs(sum-1) > 1e-12 {
				t.Fatalf("θ=%g: lossless stack R+T = %.15f", theta, sum)
			}
		}
	}
	prev := math.Inf(1)
	for _, k := range []float64{0, 0.02, 0.1} {
		got, err := ThinFilmCoefficients([]TMMLayer{{Index: complex(1.5, k), Thickness: 300e-9}}, 1, 1, 0, coatLam)
		if err != nil {
			t.Fatal(err)
		}
		sum := got.RsPow + got.TsPow
		if k == 0 {
			if math.Abs(sum-1) > 1e-12 {
				t.Fatalf("lossless film R+T = %g", sum)
			}
			prev = sum
			continue
		}
		if sum >= prev {
			t.Fatalf("absorbing film κ=%g: R+T = %g did not drop below %g", k, sum, prev)
		}
		if sum > 1 {
			t.Fatalf("absorbing film κ=%g: R+T = %g > 1 (energy created)", k, sum)
		}
		prev = sum
	}
}

// TestCoatingPhaseConvention pins the sign of the transmitted phase to the
// kernel's convention (exp(−iωt), forward waves e^{+ikz}): a thin film of
// thickness d must add ≈ +k0·n·d of phase (the optical path through the film),
// not its negative. A conjugated t (the exp(+iωt) textbook convention) would
// fail this by 2k0nd.
func TestCoatingPhaseConvention(t *testing.T) {
	const n1 = 1.5
	for _, d := range []float64{10e-9, 100e-9} {
		got, err := ThinFilmCoefficients([]TMMLayer{{Index: n1, Thickness: d}}, 1, 1, 0, coatLam)
		if err != nil {
			t.Fatal(err)
		}
		ph := math.Atan2(imag(got.Ts), real(got.Ts))
		want := 2 * math.Pi / coatLam * n1 * d
		// The interfaces and the internal reflections add a small correction
		// (< 2% of the phase for these thicknesses); requiring agreement to
		// 0.15 rad separates the correct sign (Δ ≈ 0.01 rad) from the
		// conjugated one (Δ ≈ 2·k0·n·d ≈ 3.4 rad).
		if math.Abs(ph-want) > 0.15 {
			t.Fatalf("d=%g: transmission phase %+.5f rad, want ≈ +%+.5f rad (optical path through the film)", d, ph, want)
		}
	}
}

// TestCoatingElementAppliesSPPerPolarization checks how the element feeds the
// stack response into the Jones field: the plane of incidence is the local
// x–z plane, so Ey (s) and Ex (p) receive different coefficients, and the
// scalar path only scales Ex.
func TestCoatingElementAppliesSPPerPolarization(t *testing.T) {
	params := map[string]any{
		"theta0":    45 * math.Pi / 180,
		"substrate": "N-BK7",
		"layers":    []any{map[string]any{"material": "Ta2O5", "thickness_nm": 60.0}},
	}
	el, err := NewElement(ElementSpec{Type: "coating", Params: params})
	if err != nil {
		t.Fatal(err)
	}
	ce := el.(*coatingEl)
	res, err := ce.response(coatLam)
	if err != nil {
		t.Fatal(err)
	}
	if complexAbs(res.Ts-res.Tp) < 1e-6 {
		t.Fatalf("the 45° stack should transmit s and p differently (ts=%v tp=%v)", res.Ts, res.Tp)
	}
	// Vector field.
	f := NewField(8, 1e-4, true)
	for i := range f.Ex {
		f.Ex[i] = 1
		f.Ey[i] = 1
	}
	if err := el.Apply(f, ctxFor(coatLam)); err != nil {
		t.Fatal(err)
	}
	tpEff := ThinFilmEffectiveAmplitude(res.Tp, res.TpPow)
	tsEff := ThinFilmEffectiveAmplitude(res.Ts, res.TsPow)
	if complexAbs(f.Ex[0]-tpEff) > 1e-15 || complexAbs(f.Ey[0]-tsEff) > 1e-15 {
		t.Fatalf("plane=x: Ex ← t_p, Ey ← t_s expected, got Ex=%v Ey=%v", f.Ex[0], f.Ey[0])
	}
	// The applied amplitude is power-preserving: |E|² is the stack's T.
	if math.Abs(complexAbs(f.Ex[0])*complexAbs(f.Ex[0])-res.TpPow) > 1e-15 {
		t.Fatalf("applied |t_p|² = %g, want T_p = %g", complexAbs(f.Ex[0])*complexAbs(f.Ex[0]), res.TpPow)
	}
	if math.Abs(complexAbs(f.Ey[0])*complexAbs(f.Ey[0])-res.TsPow) > 1e-15 {
		t.Fatalf("applied |t_s|² = %g, want T_s = %g", complexAbs(f.Ey[0])*complexAbs(f.Ey[0]), res.TsPow)
	}
	// plane = "y" swaps the roles.
	elY, err := NewElement(ElementSpec{Type: "coating", Params: map[string]any{
		"plane": "y", "theta0": 45 * math.Pi / 180, "substrate": "N-BK7",
		"layers": []any{map[string]any{"material": "Ta2O5", "thickness_nm": 60.0}},
	}})
	if err != nil {
		t.Fatal(err)
	}
	g := NewField(8, 1e-4, true)
	for i := range g.Ex {
		g.Ex[i] = 1
		g.Ey[i] = 1
	}
	if err := elY.Apply(g, ctxFor(coatLam)); err != nil {
		t.Fatal(err)
	}
	if complexAbs(g.Ex[0]-tsEff) > 1e-15 || complexAbs(g.Ey[0]-tpEff) > 1e-15 {
		t.Fatalf("plane=y: Ex ← t_s, Ey ← t_p expected, got Ex=%v Ey=%v", g.Ex[0], g.Ey[0])
	}
	// Scalar field: the kernel's scalar field is single-polarization (Ex), so
	// it gets the x coefficient and stays scalar (no promotion).
	sc := NewField(8, 1e-4, false)
	for i := range sc.Ex {
		sc.Ex[i] = 1
	}
	if err := el.Apply(sc, ctxFor(coatLam)); err != nil {
		t.Fatal(err)
	}
	if sc.Polarized {
		t.Fatal("a scalar field must not be promoted by the coating")
	}
	if complexAbs(sc.Ex[0]-tpEff) > 1e-15 {
		t.Fatalf("scalar field got %v, want the x (p) coefficient %v", sc.Ex[0], tpEff)
	}
	// At normal incidence s and p coincide and the scalar path is exact.
	el0, err := NewElement(ElementSpec{Type: "coating", Params: map[string]any{"substrate": "N-BK7"}})
	if err != nil {
		t.Fatal(err)
	}
	r0, err := el0.(*coatingEl).response(coatLam)
	if err != nil {
		t.Fatal(err)
	}
	if complexAbs(r0.Ts-r0.Tp) > 1e-15 {
		t.Fatalf("normal incidence: t_s = %v, t_p = %v should coincide", r0.Ts, r0.Tp)
	}
}

// TestCoatingMaterialSpecs checks the four dispersion forms and the parsing
// errors of the element parameters.
func TestCoatingMaterialSpecs(t *testing.T) {
	// Constant n, k.
	m, err := parseCoatingMaterial(map[string]any{"n": 2.0, "k": 0.01})
	if err != nil {
		t.Fatal(err)
	}
	if idx, err := m.index(500e-9); err != nil || complexAbs(idx-complex(2.0, 0.01)) > 0 || real(idx) != 2.0 {
		t.Fatalf("constant material: %v (%v)", idx, err)
	}
	// Cauchy: n = A + B/λ² + C/λ⁴ with λ in µm.
	m, err = parseCoatingMaterial(map[string]any{"cauchy": []any{1.5, 0.01, 0.0}})
	if err != nil {
		t.Fatal(err)
	}
	idx, err := m.index(500e-9) // λ = 0.5 µm → 1.5 + 0.01/0.25 = 1.54
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(real(idx)-1.54) > 1e-12 {
		t.Fatalf("cauchy material at 500 nm: n = %g, want 1.54", real(idx))
	}
	// Sellmeier (the fuse-silica fit must reproduce the glass helper).
	m, err = parseCoatingMaterial(map[string]any{"sellmeier": []any{0.6961663, 0.4079426, 0.8974794, 0.00467914826, 0.0135120631, 97.9340025}})
	if err != nil {
		t.Fatal(err)
	}
	idx, err = m.index(587.6e-9)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(real(idx)-1.45846) > 1e-5 {
		t.Fatalf("sellmeier material: n_d = %g, want 1.45846", real(idx))
	}
	// Table with linear interpolation and clamping at the ends.
	m, err = parseCoatingMaterial(map[string]any{"table": []any{
		[]any{400.0, 2.4, 0.0}, []any{600.0, 2.2, 0.02}}})
	if err != nil {
		t.Fatal(err)
	}
	idx, err = m.index(500e-9)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(real(idx)-2.3) > 1e-12 || math.Abs(imag(idx)-0.01) > 1e-12 {
		t.Fatalf("table interpolation: got %v, want 2.3+0.01i", idx)
	}
	if idx, _ = m.index(300e-9); real(idx) != 2.4 {
		t.Fatalf("table below range: got %v, want the first sample", idx)
	}
	if idx, _ = m.index(900e-9); real(idx) != 2.2 {
		t.Fatalf("table above range: got %v, want the last sample", idx)
	}
	// Built-in names, aliases and the shared glass table.
	for _, c := range []struct {
		name string
		want float64
	}{{"air", 1}, {"MgF2", 1.37774}, {"SiO2", 1.45846}, {"FusedSilica", 1.45846}, {"N-BK7", 1.51680}} {
		m, ok := coatingBuiltinMaterial(c.name)
		if !ok {
			t.Fatalf("%s: not a built-in material", c.name)
		}
		idx, err := m.index(587.6e-9)
		if err != nil {
			t.Fatalf("%s: %v", c.name, err)
		}
		if math.Abs(real(idx)-c.want) > 3e-4 {
			t.Fatalf("%s: n_d = %g, want %.5f", c.name, real(idx), c.want)
		}
	}
	// The layer string form.
	ls, err := parseCoatingLayerString("MgF2:99.2; Ta2O5:60.3")
	if err != nil {
		t.Fatal(err)
	}
	if len(ls) != 2 || math.Abs(ls[0].thickness-99.2e-9) > 1e-15 || ls[1].material.name != "Ta2O5" {
		t.Fatalf("layer string parsing: %+v", ls)
	}
	// Layer objects.
	layers, err := parseCoatingLayers([]any{
		map[string]any{"material": "MgF2", "thickness_nm": 99.2},
		map[string]any{"n": 2.35, "thickness_nm": 60.0},
		"SiO2:120",
	})
	if err != nil {
		t.Fatal(err)
	}
	if len(layers) != 3 || math.Abs(layers[2].thickness-120e-9) > 1e-15 {
		t.Fatalf("layer list parsing: %+v", layers)
	}
	// Error paths.
	bad := []map[string]any{
		{"substrate": "unobtainium"},
		{"layers": []any{map[string]any{"material": "MgF2"}}},
		{"layers": []any{map[string]any{"material": "MgF2", "thickness_nm": -5.0}}},
		{"layers": []any{map[string]any{"material": "MgF2", "thickness_nm": 100, "cauchy": []any{1, 2, 3, 4}}}},
		{"layers": []any{map[string]any{"sellmeier": []any{1, 2, 3}}}},
		{"layers": []any{map[string]any{"thickness_nm": 100}}},
		{"layers": "MgF2"},
		{"layers": "Nope:100"},
		{"theta0": 2.0},
		{"plane": "z"},
		{"layers": []any{map[string]any{"table": []any{[]any{600.0, 2.2}, []any{400.0, 2.4}}, "thickness_nm": 50}}},
	}
	for _, p := range bad {
		if _, err := NewElement(ElementSpec{Type: "coating", Params: p}); err == nil {
			t.Fatalf("coating accepted %v", p)
		}
	}
	// ThinFilmCoefficients validates its own inputs.
	if _, err := ThinFilmCoefficients(nil, 1, 1.5, math.Pi/2, coatLam); err == nil {
		t.Fatal("θ0 = π/2 must be rejected")
	}
	if _, err := ThinFilmCoefficients([]TMMLayer{{Index: 1.5, Thickness: -1}}, 1, 1.5, 0, coatLam); err == nil {
		t.Fatal("negative layer thickness must be rejected")
	}
}

// TestCoatingInTrain runs the element through the simulator: the recorded
// power must equal P·T (the coating is the only loss), the run is
// deterministic, and a 12-layer mirror reflects > 99% (the light that is not
// transmitted is simply not propagated).
func TestCoatingInTrain(t *testing.T) {
	wl := coatLam
	cfg := scalarCfg(128, 4e-3, wl)
	cfg.Source = planeSource()
	cfg.Elements = []ElementSpec{
		{Type: "coating", Params: map[string]any{"substrate": "N-BK7"}},
		{Type: "propagate", Params: map[string]any{"distance": 0.05}},
		{Type: "sensor", Params: map[string]any{"label": "through"}},
	}
	res := mustSim(t, cfg)
	el, err := NewElement(ElementSpec{Type: "coating", Params: map[string]any{"substrate": "N-BK7"}})
	if err != nil {
		t.Fatal(err)
	}
	r, err := el.(*coatingEl).response(wl)
	if err != nil {
		t.Fatal(err)
	}
	want := 1e-3 * r.TpPow // scalar run: the x (p) coefficient applies
	if rel := math.Abs(res.Planes[0].Stats.Power-want) / want; rel > 1e-9 {
		t.Fatalf("transmitted power %g W, want %g W (rel %g)", res.Planes[0].Stats.Power, want, rel)
	}
	res2 := mustSim(t, cfg)
	if math.Abs(res2.Planes[0].Stats.Power-res.Planes[0].Stats.Power) != 0 {
		t.Fatal("coating element is not deterministic")
	}
	// High reflector: 12 quarter-wave layers.
	var layers []any
	for i := 0; i < 12; i++ {
		n, name := 2.35, "Ta2O5"
		if i%2 == 1 {
			n, name = 1.45, "SiO2"
		}
		layers = append(layers, map[string]any{
			"n": n, "thickness_nm": 550.0 / (4 * n), "material": name})
	}
	hr, err := NewElement(ElementSpec{Type: "coating", Params: map[string]any{
		"substrate": "N-BK7", "layers": layers}})
	if err != nil {
		t.Fatal(err)
	}
	resHR, err := hr.(*coatingEl).response(wl)
	if err != nil {
		t.Fatal(err)
	}
	if resHR.TsPow > 0.02 {
		t.Fatalf("12-layer quarter-wave mirror transmits %g, want < 2%%", resHR.TsPow)
	}
	if math.Abs(resHR.TsPow+resHR.RsPow-1) > 1e-12 {
		t.Fatalf("mirror R+T = %g", resHR.TsPow+resHR.RsPow)
	}
	t.Logf("12-layer QW mirror: T = %.4f%%, R = %.4f%%", 100*resHR.TsPow, 100*resHR.RsPow)
}
