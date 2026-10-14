package geo

import (
	"fmt"
	"math"
)

// ---------------------------------------------------------------------------
// Third-order (Seidel) aberration sums
// ---------------------------------------------------------------------------

// This file implements the classical third-order aberration theory of
// rotationally symmetric systems (Seidel 1857; the form used here follows
// Welford, "Aberrations of Optical Systems", chapter 8).
//
// Notation for one surface j: c = 1/R is the vertex curvature, n and n' the
// refractive indices before and after it (n' = −n for a mirror, which is the
// standard way of folding a reflecting surface into the refraction
// recurrence), y the height and u the slope dy/dz of the paraxial marginal
// ray, ȳ and ū the same for the paraxial chief ray of the field, and
//
//	A  = n·(c·y  + u)   the marginal-ray refraction invariant  (n·i),
//	Ā  = n·(c·ȳ + ū)   the chief-ray refraction invariant,
//	H  = n·(ȳ·u − y·ū) the (system-constant) Lagrange invariant,
//	Δ(u/n) = u'/n' − u/n,   Δ(1/n) = 1/n' − 1/n.
//
// The five sums are
//
//	S1 = −Σ A²·y·Δ(u/n)              spherical aberration
//	S2 = −Σ A·Ā·y·Δ(u/n)             coma
//	S3 = −Σ Ā²·y·Δ(u/n)              astigmatism
//	S4 = −Σ H²·c·Δ(1/n)              field curvature (Petzval term)
//	S5 = −Σ [ (Ā³/A)·y·Δ(u/n) + (Ā/A)·H²·c·Δ(1/n) ]   distortion
//
// with the conic/aspheric contribution of a surface added to S1 (the r⁴ sag
// term of a conic of constant K contributes a wavefront term
// −(n'−n)·(K·c³/8 + A4)·r⁴, cf. asphereParams):
//
//	S1 += (n'−n)·(K·c³ + 8·A4)·y⁴ .
//
// Every sum is in mm when the marginal ray is normalized to the pupil radius
// passed to SeidelCoefficients (the sums scale as y₀⁴). The aberration
// function they generate is, with the pupil coordinate ρ = r/y₀ and the
// relative field H̃ = 1 at the field of the solution,
//
//	W(ρ, θ) = W040·ρ⁴ + W131·H̃·ρ³cosθ + W222·H̃²·ρ²cos²θ
//	        + W220·H̃²·ρ² + W311·H̃³·ρ·cosθ                (mm)
//
// with the wavefront coefficients of WavefrontCoefficients. The sign
// convention is the physical one of this package: W is the optical path of a
// ray to the reference sphere centred at the image point, measured relative to
// the chief ray, so W040 = −S1/8 (some texts print W040 = +S1/8 — the sign of
// the wavefront is then defined the other way round). The relations were
// verified against real ray tracing in the test suite (small-aperture limit),
// see seidel_test.go.

// seidelSurface is the per-surface data of the third-order sums.
type seidelSurface struct {
	c      float64
	nB, nA float64 // indices before/after (nA = −nB for a mirror)
	t      float64 // distance to the next surface (mm)
	mirror bool
	K, A4  float64 // conic constant and r⁴ aspheric coefficient (mm⁻³)
	label  string
}

// asphereParams returns the conic constant K and the r⁴ coefficient A4 of a
// profile. Profiles that do not expose them (a custom Profile implementation)
// are treated as spheres (K = 0, A4 = 0).
func asphereParams(p Profile) (K, A4 float64) {
	switch v := p.(type) {
	case Conic:
		return v.K, 0
	case *Conic:
		if v != nil {
			return v.K, 0
		}
	case Sphere, *Sphere, Plane, *Plane, nil:
		return 0, 0
	case EvenAsphere:
		if len(v.Coeffs) > 0 {
			A4 = v.Coeffs[0]
		}
		return v.K, A4
	case *EvenAsphere:
		if v != nil {
			if len(v.Coeffs) > 0 {
				A4 = v.Coeffs[0]
			}
			return v.K, A4
		}
	}
	return 0, 0
}

// seidelSurfaces extracts the paraxial description of every surface at the
// given wavelength, folding mirrors into the refraction model with n' = −n.
func (s *System) seidelSurfaces(lambdaNm float64) []seidelSurface {
	out := make([]seidelSurface, len(s.Surfaces))
	for i, sf := range s.Surfaces {
		t := 0.0
		if i < len(s.Thickness) {
			t = s.Thickness[i]
		}
		nB := s.Mediums[i].N(lambdaNm)
		nA := s.Mediums[i+1].N(lambdaNm)
		if sf.Mirror {
			nA = -nB
		}
		K, A4 := asphereParams(sf.Profile)
		out[i] = seidelSurface{
			c: curvatureOf(sf.Profile), nB: nB, nA: nA, t: t,
			mirror: sf.Mirror, K: K, A4: A4, label: sf.Label(),
		}
	}
	return out
}

// traceThirdOrderRay propagates the paraxial ray (y, u) through the elements
// with the standard recurrence
//
//	refraction: n'·u' = n·u − y·Φ,  Φ = (n' − n)·c
//	transfer:   y' = y + t·u'
//
// and returns the height/slope at the incidence side of every surface plus the
// slope after refraction.
func traceThirdOrderRay(els []seidelSurface, y, u float64) (ys, us, usAfter []float64) {
	n := len(els)
	ys, us, usAfter = make([]float64, n), make([]float64, n), make([]float64, n)
	for i := range els {
		ys[i], us[i] = y, u
		phi := (els[i].nA - els[i].nB) * els[i].c
		u2 := (els[i].nB*u - y*phi) / els[i].nA
		usAfter[i] = u2
		y += els[i].t * u2
		u = u2
	}
	return ys, us, usAfter
}

// SeidelTerm is the third-order contribution of one surface to the sums.
type SeidelTerm struct {
	Surface            int
	Label              string
	Mirror             bool
	S1, S2, S3, S4, S5 float64
	Petzval            float64 // c·Δ(1/n) of the surface
}

// Seidel holds the third-order aberration sums of a system for one field and
// wavelength, together with the data that interprets them. See the file
// comment of seidel.go for the formulas and conventions.
type Seidel struct {
	Lambda      float64
	Field       Field
	PupilRadius float64 // marginal-ray height y₀ at surface 0 (mm)
	H           float64 // Lagrange invariant of the marginal/chief ray pair
	// The five sums, in mm.
	S1, S2, S3, S4, S5 float64
	// Petzval is the Petzval sum Σ c·Δ(1/n) (1/mm); its reciprocal is the
	// Petzval curvature radius.
	Petzval float64
	// Terms holds the per-surface contributions (their sums are S1…S5 and
	// Petzval).
	Terms []SeidelTerm
	// MarginalHeight and MarginalSlope are the exit data of the marginal ray
	// at the last surface (its height yₖ in mm and slope uₖ' = dy/dz).
	MarginalHeight, MarginalSlope float64
	// IndexImageSide is the refractive index after the last surface.
	IndexImageSide float64
	// ImageHeight is the paraxial image height of the field at the plane in
	// which the exit rays of the marginal ray would cross the axis — provided
	// for reporting; it is NaN when the system has no real image (uₖ' = 0).
	ImageHeight float64
}

// SeidelCoefficients computes the third-order aberration sums of the system at
// wavelength lambdaNm for the field f. pupilRadius is the height of the
// paraxial marginal ray at surface 0 (normally the semi-diameter of the
// entrance pupil; the sums scale with its fourth power — pass the same radius
// that normalizes the pupil coordinate ρ of the samplings).
//
// The field enters only through its meridional magnitude |h| = hypot of the
// field tangents (or object heights); the system is rotationally symmetric, so
// a field of a given magnitude has the aberrations computed here up to a
// rotation about the axis.
//
// Mirrors are supported through the n' = −n folding rule of the paraxial
// recurrence, but the object must still lie in front of the first surface in
// the (possibly folded) layout, and the thicknesses must describe that layout.
// The sums are third-order: they assume small apertures and fields and
// spherical/aspheric (not general freeform) surfaces.
func (s *System) SeidelCoefficients(lambdaNm float64, f Field, pupilRadius float64) (*Seidel, error) {
	if err := s.Validate(); err != nil {
		return nil, err
	}
	if !finitePositive(pupilRadius) {
		return nil, fmt.Errorf("geo: Seidel: pupil radius %g is not positive", pupilRadius)
	}
	if f.Finite && !(f.ObjectZ < s.VertexZ(0)) {
		return nil, fmt.Errorf("geo: Seidel: object at z=%g is not in front of surface 0 (z=%g)", f.ObjectZ, s.VertexZ(0))
	}
	els := s.seidelSurfaces(lambdaNm)
	n := len(els)
	for i := range els {
		// A mirror is folded in as n' = −n, so only the magnitude of the two
		// indices may be checked here.
		if !finitePositive(els[i].nB) || !finitePositive(math.Abs(els[i].nA)) {
			return nil, fmt.Errorf("geo: Seidel: zero or non-finite refractive index at surface %d", i)
		}
	}

	// Marginal ray: height pupilRadius at the plane of surface 0.
	var y0, u0 float64
	y0 = pupilRadius
	if f.Finite {
		u0 = -y0 / f.ObjectZ // ray from the axial object point (z = ObjectZ) to y₀ at z = 0
	}
	ysM, usM, uaM := traceThirdOrderRay(els, y0, u0)

	// Chief ray of the field: the ray that crosses the stop vertex. Its
	// starting data at the plane of surface 0 follow from the linear
	// condition "the traced ray passes through the stop centre", evaluated
	// with the ray transfer matrix from surface 0 to the stop plane (this is
	// the true paraxial chief ray; for a stop at surface 0 it reduces to the
	// straight line through the stop vertex).
	z0 := s.VertexZ(0)
	tx, ty := f.Transverse()
	hMag := math.Hypot(tx, ty)
	// (y, u) matrix from the plane of surface 0 to the plane of the stop.
	ma, mb := 1.0, 0.0
	mc, md := 0.0, 1.0
	for j := 0; j < s.StopIndex; j++ {
		phi := (els[j].nA - els[j].nB) * els[j].c
		// refraction: u' = −(φ/n')·y + (n/n')·u
		na, nb := ma, mb
		nc := -phi/els[j].nA*ma + els[j].nB/els[j].nA*mc
		nd := -phi/els[j].nA*mb + els[j].nB/els[j].nA*md
		// transfer over the following gap: y' = y + t·u
		ma, mb = na+els[j].t*nc, nb+els[j].t*nd
		mc, md = nc, nd
	}
	var uC0, yC0 float64
	if f.Finite {
		// The ray leaves the object point (hMag, ObjectZ): y0 = hMag + u0·D
		// with D the distance to the plane of surface 0; y_stop = a·y0+b·u0.
		D := z0 - f.ObjectZ
		den := ma*D + mb
		if den == 0 {
			return nil, fmt.Errorf("geo: Seidel: chief ray undefined (the stop images to infinity)")
		}
		uC0 = -ma * hMag / den
		yC0 = hMag + uC0*D
	} else {
		uC0 = -hMag
		if ma == 0 {
			return nil, fmt.Errorf("geo: Seidel: chief ray undefined (the stop images to infinity)")
		}
		yC0 = -mb * uC0 / ma
	}
	ysC, usC, uaC := traceThirdOrderRay(els, yC0, uC0)

	// Lagrange invariant (constant across the surfaces; evaluated at surface 0).
	H := els[0].nB * (ysC[0]*usM[0] - ysM[0]*usC[0])

	sd := &Seidel{Lambda: lambdaNm, Field: f, PupilRadius: pupilRadius, H: H}
	sd.Terms = make([]SeidelTerm, n)
	for j := 0; j < n; j++ {
		e := els[j]
		y := ysM[j]
		u := usM[j]
		u2 := uaM[j]
		yb := ysC[j]
		ub := usC[j]
		A := e.nB * (e.c*y + u)
		Ab := e.nB * (e.c*yb + ub)
		duN := u2/e.nA - u/e.nB // Δ(u/n)
		d1N := 1/e.nA - 1/e.nB  // Δ(1/n)
		t := SeidelTerm{Surface: j, Label: e.label, Mirror: e.mirror}
		// S1..S3 with the aspheric (conic + A4) contribution folded into S1.
		t.S1 = -A * A * y * duN
		if e.mirror || e.K != 0 || e.A4 != 0 {
			kc := e.K * e.c * e.c * e.c
			t.S1 += (e.nA - e.nB) * (kc + 8*e.A4) * y * y * y * y
		}
		t.S2 = -A * Ab * y * duN
		t.S3 = -Ab * Ab * y * duN
		t.Petzval = e.c * d1N
		t.S4 = -H * H * e.c * d1N
		// S5: the classical (Ā/A) form; the division is singular for a
		// marginal ray that reaches the surface at exactly normal incidence
		// (A = 0). For such a surface the contribution is taken as zero: this
		// is exact for a plane or non-deviating surface and an approximation
		// otherwise (the classical formula has no limit there).
		scale := math.Abs(e.nB) * (math.Abs(e.c*y) + math.Abs(u))
		if scale > 0 && math.Abs(A) > 1e-13*scale {
			t.S5 = -((Ab*Ab*Ab/A)*y*duN + (Ab/A)*H*H*e.c*d1N)
		}
		sd.Terms[j] = t
		sd.S1 += t.S1
		sd.S2 += t.S2
		sd.S3 += t.S3
		sd.S4 += t.S4
		sd.S5 += t.S5
		sd.Petzval += t.Petzval
	}
	sd.MarginalHeight = ysM[n-1]
	sd.MarginalSlope = uaM[n-1]
	sd.IndexImageSide = els[n-1].nA
	if uaM[n-1] != 0 {
		// The axial marginal ray crosses the axis at the paraxial image plane;
		// the chief ray height there is the paraxial image height of the field.
		zk := s.VertexZ(n - 1)
		zImg := zk - ysM[n-1]/uaM[n-1]
		sd.ImageHeight = ysC[n-1] + (zImg-zk)*uaC[n-1]
	} else {
		sd.ImageHeight = math.Inf(signOf(-ysC[n-1]))
	}
	return sd, nil
}

// PetzvalRadius returns the radius of the Petzval surface, 1/Petzval (mm). It
// is +Inf for a system with a vanishing Petzval sum and NaN when the sum is
// not finite. The Petzval surface curves towards the lens for a positive
// (focusing) component, so a thin lens in air of focal length f and index n
// has PetzvalRadius = −n·f.
func (sd *Seidel) PetzvalRadius() float64 {
	if sd.Petzval == 0 {
		return math.Inf(1)
	}
	return 1 / sd.Petzval
}

// WavefrontCoeffs are the third-order wavefront aberration coefficients of
// Welford's polynomial (mm), with the pupil coordinate normalized to
// PupilRadius and the field normalized to the field of the Seidel solution:
//
//	W(ρx, ρy) = W040·ρ⁴ + W131·H·ρ³cosθ + W222·H²·ρ²cos²θ
//	          + W220·H²·ρ² + W311·H³·ρ·cosθ ,  H = 1 at the solution's field.
//
// They are derived from the sums by W040 = −S1/8, W131 = −S2/2,
// W222 = −S3/2, W220 = −(S3+S4)/4 and W311 = −S5/2 (see the file comment for
// the sign convention; the relations are validated against real-ray OPD maps
// in seidel_test.go).
type WavefrontCoeffs struct {
	W040 float64 // spherical aberration (mm)
	W131 float64 // coma (mm)
	W222 float64 // astigmatism (mm)
	W220 float64 // field curvature / defocus (mm)
	W311 float64 // distortion (mm)
}

// Wavefront returns the third-order wavefront coefficients implied by the
// sums (see WavefrontCoeffs).
func (sd *Seidel) Wavefront() WavefrontCoeffs {
	return WavefrontCoeffs{
		W040: -sd.S1 / 8,
		W131: -sd.S2 / 2,
		W222: -sd.S3 / 2,
		W220: -(sd.S3 + sd.S4) / 4,
		W311: -sd.S5 / 2,
	}
}

// TransverseDistortion returns the third-order transverse image-height error
// of the field in mm, i.e. how far the image point of the chief ray falls from
// the paraxial image point (along the meridional direction of the field):
//
//	δh' = −W311 / (n'·u'ₖ) = S5 / (4·n'·u'ₖ) ,
//
// with u'ₖ = MarginalSlope and n' = IndexImageSide. The relation is the
// transverse aberration generated by the ρ·cosθ wavefront term of W311
// (the transverse displacement is the wavefront gradient divided by
// n'·u'ₖ, cf. the validated conversions in seidel_test.go). NaN when the
// marginal ray is afocal (u'ₖ = 0).
func (sd *Seidel) TransverseDistortion() float64 {
	if sd.MarginalSlope == 0 || sd.IndexImageSide == 0 {
		return math.NaN()
	}
	return -sd.Wavefront().W311 / (sd.IndexImageSide * sd.MarginalSlope)
}

// RelativeDistortion returns the relative distortion δh'/h' of the field
// (a fraction: 0.01 = 1% pincushion for a positive value, since positive h'
// and a positive image-height error both push away from the axis). NaN when
// the paraxial image height is not finite or zero (an on-axis field).
func (sd *Seidel) RelativeDistortion() float64 {
	if !finitePositive(math.Abs(sd.ImageHeight)) {
		return math.NaN()
	}
	return sd.TransverseDistortion() / sd.ImageHeight
}

// String implements fmt.Stringer.
func (sd *Seidel) String() string {
	w := sd.Wavefront()
	return fmt.Sprintf("seidel[λ=%gnm S1=%.6g S2=%.6g S3=%.6g S4=%.6g S5=%.6g P=%.6g W040=%.6g W131=%.6g]",
		sd.Lambda, sd.S1, sd.S2, sd.S3, sd.S4, sd.S5, sd.Petzval, w.W040, w.W131)
}
