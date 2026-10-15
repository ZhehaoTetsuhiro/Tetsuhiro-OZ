package optics

import (
	"encoding/json"
	"fmt"
	"math"
	"strings"
)

// ---------------------------------------------------------------------------
// Thin-film coating (coating) — isotropic multilayer stacks evaluated with the
// characteristic-matrix (transfer-matrix, TMM) method.
//
// The element is a thin element whose transmittance is *not* a position
// dependent mask but the response of a stack of sub-wavelength layers: the
// film thicknesses (tens of nm) cannot be resolved as propagation on the wave
// grid, so the stack is evaluated analytically and collapsed into one plane.
//
// Stack model (MacLeod, "Thin-Film Optical Filters", ch. 2, in the kernel's
// sign convention). Layer j has the complex index Nj = nj + i·κj (κj > 0
// absorbs, matching the kernel's exp(−iωt) / e^{+ikz} convention,
// docs/PHYSICS.md §1) and thickness dj. With the tilted optical admittance
//
//	ηj = Nj·cosθj      (s polarization, E ⊥ plane of incidence)
//	ηj = Nj/cosθj      (p polarization, E ∥ plane of incidence)
//
// and the layer phase δj = 2π·Nj·cosθj·dj/λ = k0·(Nj·cosθj)·dj, each layer is
//
//	Mj = [[ cos δj,      −i·sin δj / ηj ],
//	      [ −i·ηj·sin δj,   cos δj       ]]
//
// (the −i follows from exp(−iωt) with forward waves e^{+ikz}, see
// ThinFilmCoefficients; textbooks using exp(+iωt) write +i, which conjugates
// the phase of t). The stack matrix is M = M1·M2·…·Mk with the first layer
// next to the incident medium, and with [B; C] = M·[1; η_sub] the input
// admittance is Y = C/B. The amplitude coefficients of the stack (referred to
// the plane in front of the stack and to the exit plane behind it) are
//
//	t = 2·η_inc/(η_inc·B + C),
//	r = (η_inc·B − C)/(η_inc·B + C),
//
// with power quantities R = |r|² and T = Re(η_sub)/Re(η_inc)·|t|², so that
// R + T = 1 exactly for a lossless stack (tested).
//
// Oblique incidence uses Snell's law in complex form,
//
//	Nj·cosθj = sqrt(Nj² − N0²·sin²θ0),
//
// on the branch Re ≥ 0, Im ≥ 0 (csqrt below): propagating waves advance in +z
// and evanescent/absorbed ones decay, so absorbing layers and total internal
// reflection are handled without special cases. θ0 is the incidence angle in
// the *incident medium*.
//
// Relation to optics/berreman.go: that file is the 4×4 Berreman propagator for
// an *anisotropic bulk* medium — it is the right tool when birefringence
// couples the polarizations over a propagation distance. An isotropic stack
// keeps s and p decoupled, so the 2×2 characteristic matrix used here is both
// exact and much cheaper; nothing is duplicated (berreman.go owns the
// eigenmode machinery for ε tensors, this file owns the interface response of
// isotropic layers). Use a `medium`/`biaxial` element for bulk propagation
// and this element for the coating on its surface.
//
// Polarization in WOS. The kernel is a Jones 2-component (Ex, Ey) system, and
// the element takes the plane of incidence to be the local x–z plane unless
// the parameter `plane` says otherwise:
//
//	plane = "x" (default): s is the local y axis (Ey), p the local x axis (Ex)
//	                       Ex ← t_p·Ex,  Ey ← t_s·Ey
//	plane = "y":           s is the local x axis (Ex), p is y
//	                       Ex ← t_s·Ex,  Ey ← t_p·Ey
//
// A scalar run (Polarized = false) carries the x component only (the kernel's
// scalar field is single-polarization, docs/PHYSICS.md §7), so it receives the
// x coefficient — the p one by default, the s one with plane = "y".
//
// Amplitude applied to the field. The kernel's field does not carry the
// admittance change of the medium it travels in (a beam leaving a component
// continues on the same grid and |E|² is the power density everywhere), so the
// raw stack coefficient t would apply a |t|² loss instead of the physical
// T = Re(η_sub)/Re(η_inc)·|t|²: for a bare air→glass interface 37% instead of
// the 4.3% Fresnel reflection. The element therefore applies
// √(Re(η_sub)/Re(η_inc))·t (ThinFilmEffectiveAmplitude), which keeps the
// energy books exact — the beam loses exactly the reflected and absorbed power
// — and leaves the phase, i.e. the coating's wavefront contribution, untouched.
// ThinFilmCoefficients still returns the raw coefficients of thin-film theory.
//
// Limit (documented, not hidden): θ0 is one scalar design angle for the whole
// element. A thin element is a local operator and carries no per-point
// wavevector, so the angle-of-incidence spread of a converging/diverging beam
// (or of a tilted illumination) is not modelled: every point gets the θ0
// response. For a beam spanning Δθ the result is the θ0 answer, i.e. an
// approximation whose quality improves as the beam becomes more collimated.
// ---------------------------------------------------------------------------

// TMMLayer is one isotropic layer of a stack: complex refractive index
// (N = n + i·κ, κ > 0 absorbs) and physical thickness in metres.
type TMMLayer struct {
	Index     complex128
	Thickness float64
}

// ThinFilmResponse is the complete characteristic-matrix answer for one
// wavelength and one incidence angle. Ts/Tp and Rs/Rp are complex amplitude
// coefficients (t and r of the stack); TsPow/TpPow and RsPow/RpPow are the
// corresponding power quantities (T = Re(η_sub)/Re(η_inc)·|t|², R = |r|²).
type ThinFilmResponse struct {
	Ts, Tp       complex128
	Rs, Rp       complex128
	TsPow, TpPow float64
	RsPow, RpPow float64
}

// csqrt returns the square root of z on the branch Re ≥ 0 with the sign of the
// imaginary part following Im(z): the branch that makes
// exp(i·k0·(Nj·cosθj)·z) decay towards +z for absorbing layers and keeps
// propagating waves advancing. The standard math/cmplx.Sqrt branch is the same
// one, but spelling it out documents the convention the TMM relies on.
func csqrt(z complex128) complex128 {
	r := math.Hypot(real(z), imag(z))
	re := math.Sqrt(math.Max(0, (r+real(z))/2))
	im := math.Sqrt(math.Max(0, (r-real(z))/2))
	if imag(z) < 0 {
		im = -im
	}
	return complex(re, im)
}

// cosSin returns cos δ and sin δ for a complex δ = a + i·b:
// cos = cos a·cosh b − i·sin a·sinh b, sin = sin a·cosh b + i·cos a·sinh b.
func complexCosSin(d complex128) (complex128, complex128) {
	a, b := real(d), imag(d)
	ca, sa := math.Cos(a), math.Sin(a)
	cb, sb := math.Cosh(b), math.Sinh(b)
	return complex(ca*cb, -sa*sb), complex(sa*cb, ca*sb)
}

// ThinFilmCoefficients evaluates an isotropic multilayer stack between two
// homogeneous media with the characteristic-matrix method described in the
// file header of elem_coating.go. incident and substrate are the complex
// indices of the media around the stack, theta0 is the incidence angle (rad)
// in the incident medium, lambda the vacuum wavelength (m).
//
// The returned coefficients are referenced to the plane in front of the stack
// (incident wave) and to the plane behind it (transmitted wave), the usual
// convention of thin-film theory: the sub-micron physical thickness is inside
// the matrix, and the few-hundred-nanometre offset of the two reference planes
// is part of the model, not a missing term. Concretely a bare slab limit
// reproduces the propagation phase, t → exp(i·k0·N·d) (tested).
//
// Sign convention: the kernel's time factor is exp(−iωt) with forward waves
// e^{+ikz} (docs/PHYSICS.md §1), which makes the off-diagonal entries of the
// layer matrix −i·sin δ/η (and −i·η·sin δ). Textbooks that write +i (MacLeod,
// "Thin-Film Optical Filters", ch. 2) use the conjugate convention exp(+iωt);
// |t|, R and T are the same either way, the *phase* of t is not, and the phase
// is what this element contributes to the wavefront — so the kernel's
// convention is used, and elem_coating_test.go pins the phase against a direct
// solution of the interface boundary conditions.
func ThinFilmCoefficients(layers []TMMLayer, incident, substrate complex128, theta0, lambda float64) (ThinFilmResponse, error) {
	var out ThinFilmResponse
	if !(lambda > 0) || math.IsInf(lambda, 0) {
		return out, fmt.Errorf("coating: wavelength must be positive and finite")
	}
	if math.IsNaN(theta0) || theta0 <= -math.Pi/2 || theta0 >= math.Pi/2 {
		return out, fmt.Errorf("coating: incidence angle must be in (−π/2, π/2) rad")
	}
	if incident == 0 || substrate == 0 {
		return out, fmt.Errorf("coating: incident and substrate media need a non-zero refractive index")
	}
	for i, L := range layers {
		if L.Thickness < 0 || math.IsNaN(L.Thickness) || math.IsInf(L.Thickness, 0) {
			return out, fmt.Errorf("coating: layer %d thickness must be finite and >= 0", i)
		}
		if L.Index == 0 {
			return out, fmt.Errorf("coating: layer %d has a zero refractive index", i)
		}
	}
	k0 := 2 * math.Pi / lambda
	// Snell invariant N0·sinθ0 (the same for every medium).
	sn := incident * complex(math.Sin(theta0), 0)

	// Tilted admittances. η = N·cosθ = sqrt(N² − (N0 sinθ0)²) for s, and
	// η = N²/(N·cosθ) for p; the incident and substrate values also give the
	// power normalisation Re(η)/Re(η0).
	etaIncS := csqrt(incident*incident - sn*sn)
	etaSubS := csqrt(substrate*substrate - sn*sn)
	etaIncP := incident * incident / etaIncS
	etaSubP := substrate * substrate / etaSubS

	var resS, resP ThinFilmResponse
	for _, pol := range [2]int{0, 1} { // 0 = s, 1 = p
		m11 := complex(1, 0)
		m12 := complex(0, 0)
		m21 := complex(0, 0)
		m22 := complex(1, 0)
		for _, L := range layers {
			if L.Thickness == 0 {
				continue // a zero-thickness layer is a no-op
			}
			cosj := csqrt(L.Index*L.Index - sn*sn) // Nj·cosθj
			var eta complex128
			if pol == 0 {
				eta = cosj
			} else {
				eta = L.Index * L.Index / cosj
			}
			delta := complex(k0*L.Thickness, 0) * cosj
			cd, sd := complexCosSin(delta)
			// Mj = [[cd, −i·sd/η], [−i·η·sd, cd]] (kernel sign convention, see
			// ThinFilmCoefficients) and M ← M·Mj: the characteristic matrix of
			// the stack is M1·M2·…·Mk with the first layer of the list
			// adjacent to the incident medium. The order is not commutative
			// for a real design (the reversed product describes the mirror
			// stack) — pinned by a test.
			negI := complex(0, -1)
			mj11, mj12, mj21, mj22 := cd, negI*sd/eta, negI*eta*sd, cd
			n11 := m11*mj11 + m12*mj21
			n12 := m11*mj12 + m12*mj22
			n21 := m21*mj11 + m22*mj21
			n22 := m21*mj12 + m22*mj22
			m11, m12, m21, m22 = n11, n12, n21, n22
		}
		var etaInc, etaSub complex128
		if pol == 0 {
			etaInc, etaSub = etaIncS, etaSubS
		} else {
			etaInc, etaSub = etaIncP, etaSubP
		}
		b := m11 + m12*etaSub
		c := m21 + m22*etaSub
		// 2·η0 = η0·E(0) + H(0) = t·(η0·B + C) from the interface boundary
		// conditions, so t = 2η0/(η0·B + C) and r = (η0·B − C)/(η0·B + C).
		den := etaInc*b + c
		t := 2 * etaInc / den
		r := (etaInc*b - c) / den
		rPow := real(r)*real(r) + imag(r)*imag(r)
		tPow := real(t)*real(t) + imag(t)*imag(t)
		if re := real(etaInc); re != 0 {
			tPow *= real(etaSub) / re
		}
		if pol == 0 {
			resS.Ts, resS.Rs, resS.TsPow, resS.RsPow = t, r, tPow, rPow
		} else {
			resP.Tp, resP.Rp, resP.TpPow, resP.RpPow = t, r, tPow, rPow
		}
	}
	out = ThinFilmResponse{
		Ts: resS.Ts, Tp: resP.Tp, Rs: resS.Rs, Rp: resP.Rp,
		TsPow: resS.TsPow, TpPow: resP.TpPow, RsPow: resS.RsPow, RpPow: resP.RpPow,
	}
	return out, nil
}

// ThinFilmEffectiveAmplitude converts a raw stack amplitude coefficient t (the
// field ratio between the incident plane in front of the stack and the exit
// plane inside the substrate) into the multiplier a thin element of THIS kernel
// applies to the field:
//
//	√(Re(η_sub)/Re(η_inc))·t,   with |result|² = T = Re(η_sub)/Re(η_inc)·|t|²
//
// The kernel's fields do not carry the impedance/admittance change of the
// medium they travel in (a beam after a component continues on the same grid,
// |E|² is the power density everywhere), so multiplying by the raw t would
// silently apply a |t|² loss instead of the physical T: for a bare air→glass
// interface that is 37% instead of the 4.3% Fresnel reflection. Applying the
// power-preserving amplitude keeps the energy books right (the reflected part
// is simply not propagated, exactly like a mirror's reflectance), and the phase
// of the coating — the part that matters for wavefronts and interferometers —
// is unchanged. Documented in docs/PHYSICS.md §5.2.
func ThinFilmEffectiveAmplitude(t complex128, T float64) complex128 {
	mag2 := real(t)*real(t) + imag(t)*imag(t)
	if mag2 == 0 || T <= 0 {
		return 0
	}
	return t * complex(math.Sqrt(T/mag2), 0)
}

// ---- materials -------------------------------------------------------------

// coatingMaterialKind tags how a material's index is evaluated.
type coatingMaterialKind int

const (
	coatingConstant  coatingMaterialKind = iota // n + i·k, no dispersion
	coatingCauchy                               // n(λ) = A + B/λ² + C/λ⁴, λ in µm
	coatingSellmeier                            // n² = 1 + Σ B_i·λ²/(λ² − C_i), λ in µm
	coatingTable                                // (λ[nm], n, k) samples, linear in λ
	coatingGlass                                // built-in Sellmeier glass (thick lens table)
)

// coatingMaterial is a parsed layer/substrate material.
type coatingMaterial struct {
	kind   coatingMaterialKind
	name   string
	n, k   float64      // constant
	cauchy [3]float64   // A, B, C (λ in µm)
	sell   [6]float64   // B1..B3, C1..C3 (λ in µm, C in µm²)
	tab    [][3]float64 // λ (m), n, k — ascending in λ
	glass  string       // built-in glass name
}

// coatingFilmMaterials are the built-in layer materials. The list is small on
// purpose: coatings in practice use a handful of standard films, and any other
// material can be given inline with the four dispersion forms
// (n/k, cauchy, sellmeier, table). Provenance:
//
//   - MgF2: Sellmeier fit of the ordinary ray, M. J. Dodge, "Refractive
//     properties of magnesium fluoride", Appl. Opt. 23, 1980–1985 (1984),
//     0.2–7 µm. C_i are the squared oscillator wavelengths of that fit.
//   - SiO2: fused-silica Sellmeier fit, I. H. Malitson, J. Opt. Soc. Am. 55,
//     1205 (1965) — the same entry the thick lens uses (glass "FusedSilica").
//   - Ta2O5: constant-index approximation n = 2.146, κ = 4.2e-6, the value of
//     the ion-assisted e-beam film of L. Gao, F. Lemarchand, M. Lequime,
//     Opt. Express 20, 15734 (2012) as published in the refractiveindex.info
//     database at 587.6 nm. It carries no dispersion; use an inline Cauchy,
//     Sellmeier or table model when the spectral shape matters.
var coatingFilmMaterials = []coatingMaterial{
	{kind: coatingSellmeier, name: "MgF2",
		sell: [6]float64{0.48755108, 0.39875031, 2.3120353,
			0.001882178397, 0.008951888472, 566.1355913}},
	{kind: coatingGlass, name: "SiO2", glass: "FusedSilica"},
	{kind: coatingConstant, name: "Ta2O5", n: 2.146, k: 4.2e-6},
}

// coatingMaterialNames lists the built-in names usable for layers (the
// substrates additionally accept the glasses of thickLensGlasses).
func coatingMaterialNames() []string {
	out := []string{"air", "vacuum"}
	for _, m := range coatingFilmMaterials {
		out = append(out, m.name)
	}
	out = append(out, ThickLensGlassNames()...)
	return out
}

// coatingBuiltinMaterial resolves a built-in material name, case-insensitively.
func coatingBuiltinMaterial(name string) (coatingMaterial, bool) {
	key := glassKey(name)
	if key == "air" || key == "vacuum" || key == "1" {
		return coatingMaterial{kind: coatingConstant, name: "air", n: 1}, true
	}
	for _, m := range coatingFilmMaterials {
		if glassKey(m.name) == key {
			if m.kind == coatingGlass {
				m.glass = resolveGlassName(m.glass)
			}
			return m, true
		}
	}
	for _, g := range thickLensGlasses {
		if glassKey(g.name) == key {
			return coatingMaterial{kind: coatingGlass, name: g.name, glass: g.name}, true
		}
	}
	return coatingMaterial{}, false
}

// resolveGlassName returns the canonical spelling of a built-in glass name.
func resolveGlassName(name string) string {
	for _, g := range thickLensGlasses {
		if glassKey(g.name) == glassKey(name) {
			return g.name
		}
	}
	return name
}

// index evaluates the material's complex index n + i·k at a wavelength in
// metres (κ > 0 absorbs).
func (m coatingMaterial) index(lambda float64) (complex128, error) {
	if !(lambda > 0) || math.IsInf(lambda, 0) {
		return 0, fmt.Errorf("wavelength must be positive and finite, got %g m", lambda)
	}
	switch m.kind {
	case coatingConstant:
		return complex(m.n, m.k), nil
	case coatingCauchy:
		um := lambda * 1e6
		l2 := um * um
		return complex(m.cauchy[0]+m.cauchy[1]/l2+m.cauchy[2]/(l2*l2), 0), nil
	case coatingSellmeier:
		g := thickLensGlass{name: m.name, b: [3]float64{m.sell[0], m.sell[1], m.sell[2]},
			c: [3]float64{m.sell[3], m.sell[4], m.sell[5]}}
		n, err := sellmeierIndex(g, lambda)
		if err != nil {
			return 0, err
		}
		return complex(n, 0), nil
	case coatingTable:
		if len(m.tab) == 0 {
			return 0, fmt.Errorf("material %q has an empty index table", m.name)
		}
		if lambda <= m.tab[0][0] {
			return complex(m.tab[0][1], m.tab[0][2]), nil
		}
		last := m.tab[len(m.tab)-1]
		if lambda >= last[0] {
			return complex(last[1], last[2]), nil
		}
		for i := 1; i < len(m.tab); i++ {
			if lambda <= m.tab[i][0] {
				a, bb := m.tab[i-1], m.tab[i]
				f := (lambda - a[0]) / (bb[0] - a[0])
				return complex(a[1]+f*(bb[1]-a[1]), a[2]+f*(bb[2]-a[2])), nil
			}
		}
		return complex(last[1], last[2]), nil
	case coatingGlass:
		n, _, _, err := ThickLensGlassIndex(m.glass, lambda)
		if err != nil {
			return 0, err
		}
		return complex(n, 0), nil
	}
	return 0, fmt.Errorf("unknown material kind")
}

// coatingMaterialJSON is the JSON form of an inline material:
//
//	{"material": "N-BK7"}                       built-in by name
//	{"n": 2.35}              {"n": 2.35, "k": 0.001}      constant
//	{"cauchy": [2.2, 0.01, 0]}                    n = A + B/λ² + C/λ⁴, λ in µm
//	{"sellmeier": [B1,B2,B3,C1,C2,C3]}            n² = 1 + Σ B_i λ²/(λ² − C_i)
//	{"table": [[400, 2.4, 0], [700, 2.3, 0]]}     λ in nm, n, k (k optional)
type coatingMaterialJSON struct {
	Material  string      `json:"material"`
	Name      string      `json:"name"`
	N         *float64    `json:"n"`
	K         *float64    `json:"k"`
	Cauchy    []float64   `json:"cauchy"`
	Sellmeier []float64   `json:"sellmeier"`
	Table     [][]float64 `json:"table"`
}

// parseCoatingMaterial turns one JSON value (a built-in name or an inline
// object) into a coatingMaterial.
func parseCoatingMaterial(v any) (coatingMaterial, error) {
	switch x := v.(type) {
	case nil:
		return coatingMaterial{}, nil // caller supplies the default
	case string:
		m, ok := coatingBuiltinMaterial(x)
		if !ok {
			return coatingMaterial{}, fmt.Errorf("unknown material %q (built-in: %s)", x, strings.Join(coatingMaterialNames(), ", "))
		}
		return m, nil
	}
	b, err := json.Marshal(v)
	if err != nil {
		return coatingMaterial{}, fmt.Errorf("material must be a name or an object")
	}
	var spec coatingMaterialJSON
	if err := json.Unmarshal(b, &spec); err != nil {
		return coatingMaterial{}, fmt.Errorf("material must be a name or an object with n/k, cauchy, sellmeier or table")
	}
	name := spec.Material
	if name == "" {
		name = spec.Name
	}
	if len(spec.Table) > 0 {
		mat := coatingMaterial{kind: coatingTable, name: orElse(name, "table")}
		for _, row := range spec.Table {
			if len(row) < 2 {
				return coatingMaterial{}, fmt.Errorf("material %q: table rows are [λ_nm, n, k?]", mat.name)
			}
			k := 0.0
			if len(row) > 2 {
				k = row[2]
			}
			mat.tab = append(mat.tab, [3]float64{row[0] * 1e-9, row[1], k})
		}
		for i := 1; i < len(mat.tab); i++ {
			if mat.tab[i][0] <= mat.tab[i-1][0] {
				return coatingMaterial{}, fmt.Errorf("material %q: table wavelengths must be strictly increasing", mat.name)
			}
		}
		return mat, nil
	}
	if len(spec.Sellmeier) > 0 {
		if len(spec.Sellmeier) != 6 {
			return coatingMaterial{}, fmt.Errorf("material %q: sellmeier needs [B1,B2,B3,C1,C2,C3]", orElse(name, "sellmeier"))
		}
		var s [6]float64
		copy(s[:], spec.Sellmeier)
		return coatingMaterial{kind: coatingSellmeier, name: orElse(name, "sellmeier"), sell: s}, nil
	}
	if len(spec.Cauchy) > 0 {
		if len(spec.Cauchy) < 1 || len(spec.Cauchy) > 3 {
			return coatingMaterial{}, fmt.Errorf("material %q: cauchy needs [A] or [A,B] or [A,B,C] (λ in µm)", orElse(name, "cauchy"))
		}
		var c [3]float64
		copy(c[:], spec.Cauchy)
		return coatingMaterial{kind: coatingCauchy, name: orElse(name, "cauchy"), cauchy: c}, nil
	}
	if spec.N != nil {
		m := coatingMaterial{kind: coatingConstant, name: orElse(name, "constant"), n: *spec.N}
		if spec.K != nil {
			m.k = *spec.K
		}
		if m.n <= 0 {
			return coatingMaterial{}, fmt.Errorf("material %q: n must be > 0", m.name)
		}
		return m, nil
	}
	if name != "" {
		m, ok := coatingBuiltinMaterial(name)
		if !ok {
			return coatingMaterial{}, fmt.Errorf("unknown material %q (built-in: %s)", name, strings.Join(coatingMaterialNames(), ", "))
		}
		return m, nil
	}
	return coatingMaterial{}, fmt.Errorf("material object needs one of material, n, cauchy, sellmeier, table")
}

func orElse(s, def string) string {
	if s != "" {
		return s
	}
	return def
}

// ---- the element -----------------------------------------------------------

// coatingLayer is one parsed layer of the stack.
type coatingLayer struct {
	material  coatingMaterial
	thickness float64 // m
}

// coatingLayerJSON is the JSON form of one layer: a "material:nm" string or an
// object with a material spec plus thickness_nm.
type coatingLayerJSON struct {
	Material    string      `json:"material"`
	Name        string      `json:"name"`
	ThicknessNM *float64    `json:"thickness_nm"`
	Thickness   *float64    `json:"thickness"` // metres, for programmatic callers
	N           *float64    `json:"n"`
	K           *float64    `json:"k"`
	Cauchy      []float64   `json:"cauchy"`
	Sellmeier   []float64   `json:"sellmeier"`
	Table       [][]float64 `json:"table"`
}

// coatingEl implements the "coating" element.
type coatingEl struct {
	incident  coatingMaterial
	substrate coatingMaterial
	layers    []coatingLayer
	theta0    float64
	planeX    bool // true: s ∥ y (default); false: s ∥ x
}

// parseCoatingLayerString parses "MgF2:99.2; TiO2:60.3" (separators ';' or ',',
// blank entries skipped) into layers.
func parseCoatingLayerString(s string) ([]coatingLayer, error) {
	var out []coatingLayer
	for _, item := range strings.FieldsFunc(s, func(r rune) bool { return r == ';' || r == ',' || r == '\n' }) {
		item = strings.TrimSpace(item)
		if item == "" {
			continue
		}
		name, num := item, ""
		if i := strings.LastIndex(item, ":"); i >= 0 {
			name, num = strings.TrimSpace(item[:i]), strings.TrimSpace(item[i+1:])
		}
		if name == "" || num == "" {
			return nil, fmt.Errorf("coating: layer %q must look like \"material:thickness_nm\"", item)
		}
		var nm float64
		if _, err := fmt.Sscanf(num, "%g", &nm); err != nil {
			return nil, fmt.Errorf("coating: layer %q has a non-numeric thickness", item)
		}
		mat, ok := coatingBuiltinMaterial(name)
		if !ok {
			return nil, fmt.Errorf("coating: unknown material %q in %q", name, item)
		}
		out = append(out, coatingLayer{material: mat, thickness: nm * 1e-9})
	}
	return out, nil
}

// parseCoatingLayers accepts a "mat:nm; mat:nm" string or a JSON list whose
// entries are strings or objects (coatingLayerJSON).
func parseCoatingLayers(v any) ([]coatingLayer, error) {
	if v == nil {
		return nil, nil
	}
	if s, ok := v.(string); ok {
		return parseCoatingLayerString(s)
	}
	b, err := json.Marshal(v)
	if err != nil {
		return nil, fmt.Errorf("layers must be a list of layers or a \"material:nm\" string")
	}
	var items []json.RawMessage
	if err := json.Unmarshal(b, &items); err != nil {
		return nil, fmt.Errorf("layers must be a list of layers or a \"material:nm; ...\" string")
	}
	out := make([]coatingLayer, 0, len(items))
	for i, raw := range items {
		var asString string
		if err := json.Unmarshal(raw, &asString); err == nil {
			one, err := parseCoatingLayerString(asString)
			if err != nil {
				return nil, err
			}
			out = append(out, one...)
			continue
		}
		var spec coatingLayerJSON
		if err := json.Unmarshal(raw, &spec); err != nil {
			return nil, fmt.Errorf("coating: layer %d must be a \"material:nm\" string or an object", i)
		}
		nm := 0.0
		switch {
		case spec.ThicknessNM != nil:
			nm = *spec.ThicknessNM
		case spec.Thickness != nil:
			nm = *spec.Thickness * 1e9
		default:
			return nil, fmt.Errorf("coating: layer %d needs thickness_nm (nm)", i)
		}
		if !(nm >= 0) || math.IsInf(nm, 0) {
			return nil, fmt.Errorf("coating: layer %d thickness must be finite and >= 0", i)
		}
		// Reuse the material parser: forward the whole object (it ignores
		// thickness_nm) so both spellings share one code path.
		matAny := map[string]any{}
		if spec.N != nil {
			matAny["n"] = *spec.N
		}
		if spec.K != nil {
			matAny["k"] = *spec.K
		}
		if spec.Cauchy != nil {
			matAny["cauchy"] = spec.Cauchy
		}
		if spec.Sellmeier != nil {
			matAny["sellmeier"] = spec.Sellmeier
		}
		if spec.Table != nil {
			matAny["table"] = spec.Table
		}
		if n := orElse(spec.Material, spec.Name); n != "" {
			matAny["material"] = n
		}
		if len(matAny) == 0 {
			return nil, fmt.Errorf("coating: layer %d has no material", i)
		}
		mat, err := parseCoatingMaterial(matAny)
		if err != nil {
			return nil, fmt.Errorf("coating: layer %d: %v", i, err)
		}
		out = append(out, coatingLayer{material: mat, thickness: nm * 1e-9})
	}
	return out, nil
}

func newCoating(p map[string]any) (Element, error) {
	e := &coatingEl{theta0: pfd(p, "theta0", 0), planeX: true}
	switch pl := strings.ToLower(strings.TrimSpace(ps(p, "plane", "x"))); pl {
	case "x", "":
		e.planeX = true
	case "y":
		e.planeX = false
	default:
		return nil, fmt.Errorf("coating: unknown plane %q (use \"x\" or \"y\")", pl)
	}
	if e.theta0 <= -math.Pi/2 || e.theta0 >= math.Pi/2 || math.IsNaN(e.theta0) {
		return nil, fmt.Errorf("coating: theta0 must be in (−π/2, π/2) rad")
	}
	// Incident medium: air by default.
	if v, ok := p["incident"]; ok && v != nil {
		m, err := parseCoatingMaterial(v)
		if err != nil {
			return nil, fmt.Errorf("coating: incident medium: %v", err)
		}
		e.incident = m
	} else {
		e.incident = coatingMaterial{kind: coatingConstant, name: "air", n: 1}
	}
	// Substrate: N-BK7 by default (a coating needs something to sit on).
	if v, ok := p["substrate"]; ok && v != nil {
		m, err := parseCoatingMaterial(v)
		if err != nil {
			return nil, fmt.Errorf("coating: substrate: %v", err)
		}
		e.substrate = m
	} else {
		e.substrate = coatingMaterial{kind: coatingGlass, name: "N-BK7", glass: "N-BK7"}
	}
	// Default stack: a single quarter-wave MgF2 antireflection layer designed
	// at 550 nm (MgF2 n = 1.3785 → d = λ/(4n) = 99.75 nm), the textbook AR
	// coating and a useful starting point.
	if v, ok := p["layers"]; ok && v != nil {
		layers, err := parseCoatingLayers(v)
		if err != nil {
			return nil, err
		}
		e.layers = layers
	} else {
		// Thickness from the material's own index, not a copied constant.
		mat := coatingFilmMaterials[0]
		n, err := mat.index(550e-9)
		if err != nil {
			return nil, fmt.Errorf("coating: default MgF2 layer: %v", err)
		}
		e.layers = []coatingLayer{{material: mat, thickness: 550e-9 / (4 * real(n))}}
	}
	return e, nil
}

// response evaluates the whole stack at one wavelength.
func (e *coatingEl) response(lambda float64) (ThinFilmResponse, error) {
	inc, err := e.incident.index(lambda)
	if err != nil {
		return ThinFilmResponse{}, fmt.Errorf("coating: incident medium: %v", err)
	}
	sub, err := e.substrate.index(lambda)
	if err != nil {
		return ThinFilmResponse{}, fmt.Errorf("coating: substrate: %v", err)
	}
	layers := make([]TMMLayer, len(e.layers))
	for i, L := range e.layers {
		N, err := L.material.index(lambda)
		if err != nil {
			return ThinFilmResponse{}, fmt.Errorf("coating: layer %d (%s): %v", i, L.material.name, err)
		}
		layers[i] = TMMLayer{Index: N, Thickness: L.thickness}
	}
	return ThinFilmCoefficients(layers, inc, sub, e.theta0, lambda)
}

// Apply multiplies the field by the stack's complex transmission, per
// polarization (see the file header for the s/p → Ex/Ey mapping). The applied
// amplitude is power-preserving: the raw TMM coefficient is rescaled by
// ThinFilmEffectiveAmplitude, so the element removes exactly the reflected (and
// absorbed) power and nothing else.
func (e *coatingEl) Apply(f *Field, ctx *Context) error {
	res, err := e.response(ctx.Wavelength)
	if err != nil {
		return err
	}
	tx := ThinFilmEffectiveAmplitude(res.Tp, res.TpPow)
	ty := ThinFilmEffectiveAmplitude(res.Ts, res.TsPow)
	if !e.planeX {
		tx, ty = ThinFilmEffectiveAmplitude(res.Ts, res.TsPow), ThinFilmEffectiveAmplitude(res.Tp, res.TpPow)
	}
	if !f.Polarized {
		// Scalar run: only Ex is carried (the kernel's scalar field is
		// single-polarization), so the x coefficient applies as-is.
		for i := range f.Ex {
			f.Ex[i] *= tx
		}
		return nil
	}
	for i := range f.Ex {
		f.Ex[i] *= tx
		f.Ey[i] *= ty
	}
	return nil
}

func init() {
	RegisterElement("coating", newCoating)
}
