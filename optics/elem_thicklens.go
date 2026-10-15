package optics

import (
	"fmt"
	"math"
	"strings"
)

// ---------------------------------------------------------------------------
// Thick lens (thick_lens) — a two-surface lens with a real centre thickness,
// modelled as one phase mask.
//
// Geometry and sign convention (the Chinese write-up is docs/PHYSICS.md §5.2):
//
//	z is the propagation direction (the component's surface normal) and r is
//	the distance from the lens axis. Surface i has vertex radius of curvature
//	Ri (metres) and conic constant ki; its sag measured from its own vertex
//	along +z is the exact conic sag
//
//	    sag_i(r) = r² / (Ri·(1 + sqrt(1 − (1+ki)·r²/Ri²)))
//	             = r²/(2Ri) + O(r⁴)                       (paraxial)
//
//	Ri = 0 or ±Inf is a plane surface (sag 0). Ri > 0 puts the centre of
//	curvature on the +z side of the vertex: for a biconvex lens the first
//	surface is R1 > 0 and the second R2 < 0, the same convention as the
//	classical lensmaker formula
//
//	    1/f = (n−1)·(1/R1 − 1/R2 + (n−1)·d/(n·R1·R2)).
//
//	Surface 2 is d behind surface 1 (d = centre thickness), so the glass
//	thickness at radius r is
//
//	    thick(r) = d + sag2(r) − sag1(r)
//	             ≈ d − (r²/2)·(1/R1 − 1/R2).
//
// Wave model (thin-element style: the whole lens is collapsed onto one plane).
// The optical path from a plane in front of the lens to a plane behind it is
//
//	    OPL(r) = (A + B − d) + n·d + (n−1)·(sag2(r) − sag1(r))
//	           = (A + B − d) + (n−1)·thick(r),
//
// because a ray at radius r crosses air down to sag1(r), glass through
// thick(r) and then air again up to B, while the on-axis ray crosses d of
// glass. The transmitted field is therefore multiplied by
//
//	    t(r) = exp( i·k0·(n(λ) − 1)·thick(r) ),      k0 = 2π/λ
//
// which is exactly the optical-path difference and gives the transverse phase
//
//	    φ(r) ≈ −k0·(n−1)·(r²/2)·(1/R1 − 1/R2) = −k0·r²/(2·f_mask),
//	    f_mask = 1/((n−1)·(1/R1 − 1/R2)).
//
// That is the kernel's converging-lens convention exp(−i·k·r²/(2f)) (see
// docs/PHYSICS.md §1), i.e. the wavefront leaving the element is a sphere
// converging towards the focus at f_mask on the +z side for a biconvex lens.
//
// Relation to the classical thick-lens formula (the "model" parameter):
//
//   - model = "sag" (default): the exact mask above. Its focal distance from
//     the element plane is f_mask = 1/((n−1)(1/R1 − 1/R2)), i.e. the d → 0
//     limit of the classical formula.
//   - model = "equivalent": the paraxial equivalent thin lens carrying the
//     full thick-lens power,
//
//	    t(r) = exp( i·k0·(n−1)·d ) · exp( −i·k0·r²/(2·f_geo) ),
//	    f_geo = 1/((n−1)·(1/R1 − 1/R2 + (n−1)·d/(n·R1·R2)))
//
//     which focuses exactly at f_geo measured from the element plane, so the
//     catalogue focal length of a real lens is reproduced. Conic constants
//     are ignored here (they belong to the sag model).
//
// The two differ by the power-stacking term (n−1)²·d/(n·R1·R2): it comes from
// the second surface sitting d/n behind the first, which a single phase plane
// cannot express — the sag model has no propagation inside the glass. For a
// 12 mm thick N-BK7 lens with R1 = +60 mm, R2 = −60 mm the focal lengths
// differ by ≈ 3.5% (f_mask = 58.2 mm, f_geo = 60.3 mm at 632.8 nm), both
// verified numerically in elem_thicklens_test.go. A physically complete thick
// lens would be the two surface phases plus a `medium(index)` propagation
// inside the glass (d/n); that is the ray-trace route this element replaces.
//
// Sampling note: a hard semi-aperture truncates a large quadratic phase, so
// its edge must be sampled (see docs/PHYSICS.md §2): the local phase step per
// pixel at the rim is k0·(n−1)·(1/R1 − 1/R2)·r·dx ≈ 2π·r·dx/(λ·f), which must
// stay below π. Below 1 pixel the edge is unresolvable anyway, so the sag is
// band-limited by the element parameter edge_sigma (an erf edge, like the
// aperture element), and the closed form of the mask is used for the phase.
// ---------------------------------------------------------------------------

// thickLensGlass is one built-in Sellmeier glass entry:
//
//	n²(λ) = 1 + Σ_i B_i·λ²/(λ² − C_i),   λ in µm, C_i in µm².
type thickLensGlass struct {
	name string
	b    [3]float64
	c    [3]float64
	lo   float64 // documented validity range, µm
	hi   float64
	src  string // provenance of the coefficients
}

// thickLensGlasses holds the built-in glasses. The kernel has no other
// dispersion facility (optics/medium.go, optics/birefringent.go and
// optics/berreman.go all take a fixed index), so the small Sellmeier table
// lives here and is shared with the coating element (elem_coating.go) for its
// substrates.
//
// Sources: the Schott data sheets publish the N-glass Sellmeier fits with
// C_i already in µm² (N-BK7, N-SF5, N-SK16, F2); fused silica is the
// Malitson 1965 fit (JOSA 55, 1205), the standard 0.21–3.7 µm dispersion of
// fused quartz.
var thickLensGlasses = []thickLensGlass{
	{"N-BK7", [3]float64{1.03961212, 0.231792344, 1.01046945},
		[3]float64{0.00600069867, 0.0200179144, 103.560653}, 0.30, 2.50,
		"Schott N-BK7 data sheet (Sellmeier fit)"},
	{"FusedSilica", [3]float64{0.6961663, 0.4079426, 0.8974794},
		[3]float64{0.00467914826, 0.0135120631, 97.9340025}, 0.21, 3.70,
		"Malitson, J. Opt. Soc. Am. 55, 1205 (1965)"},
	{"N-SF5", [3]float64{1.52481889, 0.187085527, 1.42729015},
		[3]float64{0.011254756, 0.0588995392, 129.141675}, 0.35, 2.50,
		"Schott N-SF5 data sheet (Sellmeier fit)"},
	{"N-SK16", [3]float64{1.34317774, 0.241144399, 0.994317969},
		[3]float64{0.00704687339, 0.0229005, 92.7508526}, 0.35, 2.50,
		"Schott N-SK16 data sheet (Sellmeier fit)"},
	{"F2", [3]float64{1.34533359, 0.209073176, 0.937357162},
		[3]float64{0.00997743871, 0.0470450767, 111.886764}, 0.35, 2.50,
		"Schott F2 data sheet (Sellmeier fit)"},
}

// glassKey normalises a user-supplied glass name so "N-BK7", "n bk7" and
// "nbk7" all select the same entry.
func glassKey(name string) string {
	var b strings.Builder
	for _, r := range strings.ToLower(strings.TrimSpace(name)) {
		if (r >= 'a' && r <= 'z') || (r >= '0' && r <= '9') {
			b.WriteRune(r)
		}
	}
	key := b.String()
	// Common aliases.
	switch key {
	case "sio2", "silica", "fusedsilica", "quartz":
		return "fusedsilica"
	case "nsf5":
		return "nsf5"
	case "nsk16":
		return "nsk16"
	}
	return key
}

// ThickLensGlassNames returns the built-in glass names, in table order.
func ThickLensGlassNames() []string {
	out := make([]string, len(thickLensGlasses))
	for i, g := range thickLensGlasses {
		out[i] = g.name
	}
	return out
}

// ThickLensGlassIndex returns the refractive index of a built-in glass at a
// vacuum wavelength in metres (Sellmeier fits, sources in thickLensGlasses).
// The wavelength range each fit is documented for is returned as the second
// and third value (metres) for callers that want to warn; the value is still
// evaluated outside it (extrapolation is smooth, the fit just loses accuracy).
func ThickLensGlassIndex(name string, wavelength float64) (float64, float64, float64, error) {
	key := glassKey(name)
	for _, g := range thickLensGlasses {
		if glassKey(g.name) != key {
			continue
		}
		n, err := sellmeierIndex(g, wavelength)
		return n, g.lo * 1e-6, g.hi * 1e-6, err
	}
	return 0, 0, 0, fmt.Errorf("unknown glass %q (built-in: %s)", name, strings.Join(ThickLensGlassNames(), ", "))
}

// sellmeierIndex evaluates one Sellmeier glass entry at a wavelength in metres.
func sellmeierIndex(g thickLensGlass, wavelength float64) (float64, error) {
	if !(wavelength > 0) || math.IsInf(wavelength, 0) {
		return 0, fmt.Errorf("wavelength must be positive and finite, got %g m", wavelength)
	}
	um := wavelength * 1e6
	l2 := um * um
	n2 := 1.0
	for i := 0; i < 3; i++ {
		den := l2 - g.c[i]
		if den == 0 {
			return 0, fmt.Errorf("wavelength %g m is exactly on a Sellmeier resonance of %s", wavelength, g.name)
		}
		n2 += g.b[i] * l2 / den
	}
	if !(n2 > 0) || math.IsInf(n2, 0) {
		return 0, fmt.Errorf("glass %s has no real index at %g m (inside an absorption band)", g.name, wavelength)
	}
	return math.Sqrt(n2), nil
}

// ConicSag returns the sag of a conic surface at squared radial distance r2
// (metres²): sag = r2 / (R·(1 + sqrt(1 − (1+k)·r2/R²))), the exact conic
// solution of the standard sag equation (k = 0 sphere, k = −1 paraboloid,
// k < −1 hyperboloid, k > 0 oblate ellipsoid). The argument is r² because the
// element evaluates the phase in a pixel loop where r² is already available.
//
// R = 0 or ±Inf is a plane surface (sag 0). Past the rim of the conic
// ((1+k)·r2/R² > 1) the square root would be imaginary; the sag is clamped to
// its rim value there instead of returning NaN (a NaN in a thin element
// spreads over the whole field after the next propagation — the aperture is
// expected to cut the beam before the rim, and the clamp keeps the element
// finite even when it does not).
func ConicSag(r2, R, k float64) float64 {
	if R == 0 || math.IsInf(R, 0) {
		return 0
	}
	u := (1 + k) * r2 / (R * R)
	if u > 1 {
		u = 1
	}
	return r2 / (R * (1 + math.Sqrt(1-u)))
}

// surfacePower returns 1/R with a plane surface (R = 0 or ±Inf) counting as
// zero power.
func surfacePower(R float64) float64 {
	if R == 0 || math.IsInf(R, 0) {
		return 0
	}
	return 1 / R
}

// thinLensPower is the paraxial power a two-surface lens carries when the
// glass thickness is ignored, 1/f_mask = (n−1)·(1/R1 − 1/R2).
func thinLensPower(R1, R2, n float64) float64 {
	return (n - 1) * (surfacePower(R1) - surfacePower(R2))
}

// thickLensPower is the classical thick-lens power
// 1/f_geo = (n−1)·(1/R1 − 1/R2 + (n−1)·d/(n·R1·R2)), the paraxial power of
// the two surfaces separated by d/n of glass. The last term is dropped as
// soon as one surface is plane (it is the product of the two surface powers).
func thickLensPower(R1, R2, d, n float64) float64 {
	p := thinLensPower(R1, R2, n)
	if p1, p2 := surfacePower(R1), surfacePower(R2); p1 != 0 && p2 != 0 {
		p += (n - 1) * (n - 1) * d * p1 * p2 / n
	}
	return p
}

// ThickLensMaskFocalLength returns the focal length the sag phase mask of
// thick_lens focuses to, measured from the element plane:
// 1/f = (n−1)·(1/R1 − 1/R2). It is infinite for a zero-power lens
// (both surfaces plane, or R1 = R2).
func ThickLensMaskFocalLength(R1, R2, n float64) float64 {
	p := thinLensPower(R1, R2, n)
	if p == 0 {
		return math.Inf(1)
	}
	return 1 / p
}

// ThickLensGeometricFocalLength returns the classical thick-lens focal length
// (the lens-catalogue value, measured from the principal planes)
// 1/f = (n−1)·(1/R1 − 1/R2 + (n−1)·d/(n·R1·R2)). The element applies this
// power only in model = "equivalent"; in the default sag model the focus sits
// at ThickLensMaskFocalLength, smaller by the power-stacking term.
func ThickLensGeometricFocalLength(R1, R2, d, n float64) float64 {
	p := thickLensPower(R1, R2, d, n)
	if p == 0 {
		return math.Inf(1)
	}
	return 1 / p
}

// thickLensEl implements the "thick_lens" element.
type thickLensEl struct {
	material  string  // built-in glass name, or "constant"
	nConst    float64 // index used when material = "constant"
	r1, r2    float64 // vertex radii of curvature, m (0 or ±Inf = plane)
	k1, k2    float64 // conic constants (0 = spherical)
	centre    float64 // centre thickness d, m
	aperture  float64 // semi-aperture radius, m (0 = unlimited)
	block     bool    // true (default): opaque outside the semi-aperture
	edgeSigma float64 // erf softening of the aperture edge, m (0 = hard edge)
	equiv     bool    // model = "equivalent": apply the geometric power
	x0, y0    float64
}

func newThickLens(p map[string]any) (Element, error) {
	e := &thickLensEl{
		material:  ps(p, "material", "N-BK7"),
		nConst:    pfd(p, "n", 1.5),
		r1:        pfd(p, "R1", 0.05),
		r2:        pfd(p, "R2", -0.05),
		k1:        pfd(p, "k1", 0),
		k2:        pfd(p, "k2", 0),
		centre:    pfd(p, "d", 0.004),
		aperture:  pfd(p, "aperture", 0),
		block:     pbool(p, "block", true),
		edgeSigma: pfd(p, "edge_sigma", 0),
		x0:        pfd(p, "x", 0),
		y0:        pfd(p, "y", 0),
	}
	switch m := strings.ToLower(strings.TrimSpace(ps(p, "model", "sag"))); m {
	case "sag", "wave":
		e.equiv = false
	case "equivalent", "thin", "principal":
		e.equiv = true
	default:
		return nil, fmt.Errorf("thick_lens: unknown model %q (use \"sag\" or \"equivalent\")", m)
	}
	for _, R := range []float64{e.r1, e.r2} {
		if math.IsNaN(R) {
			return nil, fmt.Errorf("thick_lens: radii of curvature must be finite or ±Inf (0 = plane)")
		}
	}
	if math.IsNaN(e.k1) || math.IsNaN(e.k2) {
		return nil, fmt.Errorf("thick_lens: conic constants must be finite")
	}
	if e.centre < 0 {
		return nil, fmt.Errorf("thick_lens: centre thickness d must be >= 0")
	}
	if e.aperture < 0 {
		return nil, fmt.Errorf("thick_lens: aperture must be >= 0")
	}
	if e.edgeSigma < 0 {
		return nil, fmt.Errorf("thick_lens: edge_sigma must be >= 0")
	}
	if e.constantGlass() {
		if !(e.nConst > 0) {
			return nil, fmt.Errorf("thick_lens: n must be > 0")
		}
	} else if _, _, _, err := ThickLensGlassIndex(e.material, 632.8e-9); err != nil {
		return nil, fmt.Errorf("thick_lens: %v", err)
	}
	// The conic sag exists only inside its own rim |R|; an aperture beyond it
	// would silently clip the sag (ConicSag clamps) and put the wrong phase on
	// the outer part of the beam, so it is rejected instead.
	for _, R := range []float64{e.r1, e.r2} {
		if e.aperture > 0 && R != 0 && !math.IsInf(R, 0) && e.aperture > math.Abs(R) {
			return nil, fmt.Errorf("thick_lens: aperture %g m exceeds the surface rim |R| = %g m", e.aperture, math.Abs(R))
		}
	}
	return e, nil
}

// constantGlass reports whether the element uses the flat `n` parameter
// instead of a Sellmeier glass.
func (e *thickLensEl) constantGlass() bool {
	key := glassKey(e.material)
	return key == "constant" || key == "custom" || key == "n"
}

// refractiveIndex evaluates n(λ) for the configured material.
func (e *thickLensEl) refractiveIndex(wavelength float64) (float64, float64, float64, error) {
	if e.constantGlass() {
		return e.nConst, 0, 0, nil
	}
	return ThickLensGlassIndex(e.material, wavelength)
}

// Apply multiplies the field by exp(i·k0·(n−1)·thick(r)) (model "sag") or by
// the paraxial equivalent thin lens (model "equivalent"); see the file header
// for the derivation and the sign convention.
func (e *thickLensEl) Apply(f *Field, ctx *Context) error {
	wl := ctx.Wavelength
	n, lo, hi, err := e.refractiveIndex(wl)
	if err != nil {
		return fmt.Errorf("thick_lens: %v", err)
	}
	if lo > 0 && (wl < lo || wl > hi) {
		ctx.Warnings.Add("thick_lens_glass_range",
			fmt.Sprintf("thick_lens: 波长 %g nm 超出玻璃 %s 的标定范围 %.0f–%.0f nm，Selmeier 拟合外推",
				wl*1e9, e.material, lo*1e9, hi*1e9), wl)
	}
	k0 := 2 * math.Pi / wl
	dn := n - 1

	// The equivalent model needs the paraxial thick-lens power of this
	// wavelength; the sag model evaluates the two surface sags per pixel.
	power := 0.0
	if e.equiv {
		power = thickLensPower(e.r1, e.r2, e.centre, n)
	}
	// The exact sag is only defined inside the surface rim |R| (ConicSag clamps
	// beyond it). If light actually reaches that region the phase is wrong, so
	// the run is told about it instead of silently clamping.
	rim := math.Inf(1)
	for _, R := range []float64{e.r1, e.r2} {
		if R != 0 && !math.IsInf(R, 0) {
			if a := math.Abs(R); a < rim {
				rim = a
			}
		}
	}
	maxActiveR2 := 0.0
	ap, sig := e.aperture, e.edgeSigma
	nx := f.N
	for j := 0; j < nx; j++ {
		dy := f.Y(j) - e.y0
		for i := 0; i < nx; i++ {
			dx := f.X(i) - e.x0
			r2 := dx*dx + dy*dy
			idx := j*nx + i
			// Aperture weight: 1 inside the semi-aperture, 0 outside, with
			// the erf edge of `edge_sigma` in between.
			w := 1.0
			if ap > 0 {
				d := ap - math.Sqrt(r2)
				if sig > 0 {
					if d < -8*sig {
						w = 0
					} else {
						w = real(smoothStep(d, sig))
					}
				} else if d < 0 {
					w = 0
				}
			}
			if w == 0 {
				if e.block {
					idx := j*nx + i
					f.Ex[idx] = 0
					if f.Polarized {
						f.Ey[idx] = 0
					}
				}
				// block = false: the light outside the clear aperture is
				// simply not refracted (t = 1), it is not absorbed.
				continue
			}
			// Track where the phase is actually applied to light: only that
			// region can hit the surface rim of the sag model.
			if r2 > maxActiveR2 && (f.Ex[idx] != 0 || (f.Polarized && f.Ey[idx] != 0)) {
				maxActiveR2 = r2
			}
			var phase float64
			if e.equiv {
				phase = k0 * (dn*e.centre - 0.5*r2*power)
			} else {
				phase = k0 * dn * (e.centre + ConicSag(r2, e.r2, e.k2) - ConicSag(r2, e.r1, e.k1))
			}
			t := cexpI(phase)
			if w != 1 {
				if e.block {
					t *= complex(w, 0)
				} else {
					// Interpolate towards the unobstructed field (t = 1).
					t = complex(1-w, 0) + complex(w, 0)*t
				}
			}
			f.Ex[idx] *= t
			if f.Polarized {
				f.Ey[idx] *= t
			}
		}
	}
	if !math.IsInf(rim, 0) && maxActiveR2 > rim*rim {
		ctx.Warnings.Add("thick_lens_sag_clip",
			fmt.Sprintf("thick_lens: 光斑延伸到 |R| = %g m 的镜面边界之外（r = %g m），该处矢高被截断为边界值",
				rim, math.Sqrt(maxActiveR2)), rim)
	}
	return nil
}

func init() {
	RegisterElement("thick_lens", newThickLens)
}
