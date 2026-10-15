package geo

import (
	"fmt"
	"math"
)

// Mat2 is a 2x2 ray transfer (ABCD) matrix. It acts on the column vector
// (y, ū) where y is the ray height (mm) and ū = n·u the reduced angle (rad,
// u = dy/dz the real paraxial slope and n the refractive index of the local
// medium).
//
// In these coordinates the elementary matrices are
//
//	transfer over the axial distance t in index n:  [[1, t/n], [0, 1]]
//	refraction at a surface of power Φ = (n'-n)c:   [[1, 0], [-Φ, 1]]
//	reflection at a mirror of curvature c:          [[1, 0], [-2·n·c, -1]]
//
// and the imaging condition between two conjugate planes is B = 0, with
// transverse magnification A. The effective focal length of a system is
// EFL = -n'/C (n' = image-side index), and the matrix of a thin lens in air of
// focal length f is [[1,0],[-1/f,1]].
//
// The mirror matrix has determinant -1: the reflection reverses the direction
// of propagation, so ū → -ū. The unfolding of a mirror system (System.Paraxial)
// is the equivalent formulation in which the index of every medium behind the
// mirror is negated (n → -n) and the reflection keeps the refraction form
// [[1,0],[-Φ,1]]; both descriptions produce the same heights and the same
// physical slopes.
type Mat2 struct {
	A, B, C, D float64
}

// Identity is the unit ray transfer matrix.
func Identity() Mat2 { return Mat2{A: 1, D: 1} }

// Mul returns the product m·o, i.e. the transfer obtained by applying o first
// and m afterwards.
func (m Mat2) Mul(o Mat2) Mat2 {
	return Mat2{
		A: m.A*o.A + m.B*o.C,
		B: m.A*o.B + m.B*o.D,
		C: m.C*o.A + m.D*o.C,
		D: m.C*o.B + m.D*o.D,
	}
}

// Det returns the determinant A·D - B·C (equals n/n' for a system between media
// of index n and n').
func (m Mat2) Det() float64 { return m.A*m.D - m.B*m.C }

// Apply propagates the ray (y, ū) through the matrix.
func (m Mat2) Apply(y, uBar float64) (float64, float64) {
	return m.A*y + m.B*uBar, m.C*y + m.D*uBar
}

// String implements fmt.Stringer.
func (m Mat2) String() string {
	return fmt.Sprintf("[[%g %g] [%g %g]]", m.A, m.B, m.C, m.D)
}

// paraxialElement is the first-order description of one surface: vertex
// curvature c = 1/R, indices before/after and the axial distance to the next
// surface. Higher-order shape terms (conic constant, aspheric coefficients) do
// not enter first-order optics.
//
// mirror marks a reflecting surface. Its two media are equal (a mirror is an
// interaction inside one medium, see System.Paraxial) and its thickness may be
// negative: the axial position of every following surface is the running sum of
// the signed thicknesses, so a fold is written exactly as the real trace reads
// it (e.g. library/geo/mirror-cassegrain-*.json has a thickness of -300 mm
// behind the primary).
type paraxialElement struct {
	c      float64 // curvature 1/mm
	nB, nA float64 // refractive index before/after
	t      float64 // distance to the next surface, mm (signed, negative after a fold)
	mirror bool    // true: reflect instead of refract
}

// surfaceMatrix returns the matrix of the interaction at the surface:
//
//	refraction: [[1,0],[-Φ,1]],  Φ = (nA - nB)·c
//	mirror:     [[1,0],[-2·nB·c,-1]]
//
// The mirror entry follows from the paraxial reflection law in the y-u frame of
// the package: at a surface of sag z = c·r²/2 the reflected slope of the ray
// line is u' = -u - 2·c·y, so with ū = n·u and nA = nB = n,
//
//	ū' = n·(-u - 2cy) = -ū - 2·n·c·y.
//
// The determinant is -1 because the reflected ray travels back along -z;
// unfolding (negating the index of the medium behind the mirror, n → -n) turns
// the same map into an ordinary refraction with Φ = -2·n·c and determinant +1,
// which is how the systems of this package are laid out (see System.Paraxial).
func (e paraxialElement) surfaceMatrix() Mat2 {
	if e.mirror {
		return Mat2{A: 1, C: -2 * e.nB * e.c, D: -1}
	}
	phi := (e.nA - e.nB) * e.c
	return Mat2{A: 1, C: -phi, D: 1}
}

// transferMatrix returns the transfer matrix of the following gap.
func (e paraxialElement) transferMatrix() Mat2 {
	if e.nA == 0 {
		return Mat2{A: 1, D: 1}
	}
	return Mat2{A: 1, B: e.t / e.nA, D: 1}
}

// curvatureOf returns the paraxial vertex curvature of a profile.
func curvatureOf(p Profile) float64 {
	if p == nil {
		return 0
	}
	return p.Curvature()
}

// paraxialElements derives the first-order element list of the system at the
// given wavelength.
func (s *System) paraxialElements(lambdaNm float64) []paraxialElement {
	els := make([]paraxialElement, len(s.Surfaces))
	for i := range s.Surfaces {
		t := 0.0
		if i < len(s.Thickness) {
			t = s.Thickness[i]
		}
		els[i] = paraxialElement{
			c:      curvatureOf(s.Surfaces[i].Profile),
			nB:     s.Mediums[i].N(lambdaNm),
			nA:     s.Mediums[i+1].N(lambdaNm),
			t:      t,
			mirror: s.Surfaces[i].Mirror,
		}
	}
	return els
}

// reverseElements returns the element list of the time-reversed system. The
// transformation is the reflection z → -z combined with the reversal of the
// direction of travel: the surfaces are listed backwards, the sign of every
// curvature flips (a surface seen from the other side has -c), the media swap
// and the mirror flags are kept. A mirror is its own time reverse
// (u' = -u - 2cy is an involution), so the construction covers reflecting
// surfaces without further changes; a gap keeps its physical value, which is
// what makes the reversed chain of a folded layout come out with the signs of
// the reversed travel.
func reverseElements(els []paraxialElement) []paraxialElement {
	n := len(els)
	out := make([]paraxialElement, n)
	for j := 0; j < n; j++ {
		e := els[n-1-j]
		out[j] = paraxialElement{c: -e.c, nB: e.nA, nA: e.nB, mirror: e.mirror}
	}
	for j := 0; j < n-1; j++ {
		out[j].t = els[n-2-j].t
	}
	return out
}

// traceParaxial propagates the paraxial ray through the elements using the
// standard y–nu (reduced angle) recurrence:
//
//	refraction: n' u' = n u - y Φ ,   Φ = (n'-n) c   (the refraction invariant
//	            A = n(yc + u) = n·i is preserved),
//	reflection: u' = -u - 2 c y      (the mirror is a refraction into an
//	            equivalent medium of index -n, see paraxialElement.surfaceMatrix),
//	transfer:   y' = y + t u' .       (t is the signed thickness)
//
// It returns the heights and slopes at the incidence side of every surface
// (ys, us, with us the real slope u = dy/dz in the physical z frame, so after an
// odd number of mirrors the ray travels towards -z and its slope keeps the sign
// of the ray line, not of the direction of travel) and the real slopes after
// each surface (usAfter).
func traceParaxial(els []paraxialElement, y, u float64) (ys, us, usAfter []float64) {
	n := len(els)
	ys = make([]float64, n)
	us = make([]float64, n)
	usAfter = make([]float64, n)
	for i := range els {
		ys[i], us[i] = y, u
		_, ubar := els[i].surfaceMatrix().Apply(y, els[i].nB*u)
		u2 := ubar / els[i].nA
		usAfter[i] = u2
		y = y + els[i].t*u2
		u = u2
	}
	return ys, us, usAfter
}

// ParaxialRay is the paraxial trace of a single ray: heights and slopes at
// every surface, in the system's global frame (surface i sits at
// Paraxial.VertexZ[i]).
type ParaxialRay struct {
	Y       []float64 // height y at the incidence side of surface i (mm)
	U       []float64 // real slope u = dy/dz at the incidence side
	UAfter  []float64 // real slope after refraction/reflection at surface i
	N       []float64 // refractive index at the incidence side
	VertexZ []float64 // vertex position of surface i (mm)
	// Y0,U0 are the starting values at surface 0.
	Y0, U0 float64
	// UBar0 is the reduced starting angle n0·u0.
	UBar0 float64
}

// HeightAt returns the ray height at axial position z (mm) by linear transfer
// from the last surface.
func (p *ParaxialRay) HeightAt(z float64) float64 {
	k := len(p.Y) - 1
	return p.Y[k] + (z-p.VertexZ[k])*p.UAfter[k]
}

// SlopeAt returns the real slope dy/dz after the last surface.
func (p *ParaxialRay) SlopeAt() float64 {
	if len(p.UAfter) == 0 {
		return 0
	}
	return p.UAfter[len(p.UAfter)-1]
}

// Paraxial holds the complete first-order (Gaussian) description of a system at
// one wavelength and one object position. All axial positions are absolute
// (mm) in the system frame described in package geo's documentation; the
// distances BFL/FFL are signed distances from the last/first surface vertex.
type Paraxial struct {
	Lambda   float64   // wavelength (nm)
	ObjectZ  float64   // object plane (mm); -Inf for an object at infinity
	Infinite bool      // object at infinity
	Index    []float64 // index in front of surface i (len = NumSurfaces+1,
	// the last entry is the image-side index)
	VertexZ []float64 // vertex position of every surface (mm)

	EFL float64 // effective focal length (mm), EFL = -y0/u'_k for a parallel input
	FFL float64 // front focal length: signed distance from surface 0 to the front focus
	BFL float64 // back focal length: signed distance from the last surface to the back focus

	FrontFocalZ     float64 // front focal point (mm)
	BackFocalZ      float64 // back focal point (mm)
	FrontPrincipalZ float64 // front principal plane (mm)
	BackPrincipalZ  float64 // back principal plane (mm)

	ImageZ        float64 // paraxial image of the axial object point (mm)
	Magnification float64 // transverse magnification at ImageZ (NaN for an infinite object)

	StopIndex int // aperture stop surface index

	// Mirrors is the number of reflecting surfaces in the chain. The
	// propagation direction of the emergent beam is Paraxial.Parity() times the
	// incident one, and the determinant of the ray transfer matrix is
	// (-1)^Mirrors·n0/nk.
	Mirrors int

	// EntrancePupilZ / ExitPupilZ are the axial positions of the images of the
	// aperture stop in the object and image space (mm). The pupil planes can
	// be virtual (inside the system). EntrancePupilRadius/ExitPupilRadius are
	// the corresponding semi-diameters (mm).
	EntrancePupilZ      float64
	ExitPupilZ          float64
	EntrancePupilRadius float64
	ExitPupilRadius     float64
	StopRadius          float64 // semi-diameter of the aperture stop (mm)

	FNumber float64 // EFL / (2·EntrancePupilRadius)
	ImageNA float64 // n'·sin(u') of the marginal ray at the entrance pupil edge

	// EntrancePupilMag is the A entry of the transfer from the entrance pupil
	// plane to the stop plane (y_stop = A·y_entrance pupil), so the entrance
	// pupil radius is StopRadius/|EntrancePupilMag|. ExitPupilMag is the A
	// entry of the transfer from the stop plane to the exit pupil plane
	// (y_exit pupil = A·y_stop), so the exit pupil radius is
	// |ExitPupilMag|·StopRadius.
	EntrancePupilMag float64
	ExitPupilMag     float64

	els []paraxialElement
}

// Paraxial derives the complete first-order description of the system at
// wavelength lambdaNm for an object at axial position objectZ (mm, use
// math.Inf(-1) for an object at infinity; objectZ must be < the first vertex
// otherwise).
//
// The paraxial engine traces the standard y–u (y–n·u) paraxial rays:
//
//	marginal ray: y = 1 mm at surface 0, u = 0 (parallel input),
//	chief ray:    passes through the centre of the aperture stop.
//
// # Mirrors and unfolding
//
// Coaxial systems with reflecting surfaces are supported. Values are reported in
// the physical frame of the package: the vertex of surface i sits at
// VertexZ(i) = OriginZ + Σ Thickness[0..i-1] and the thicknesses keep their
// signs, so a fold is written exactly as the real ray trace reads it (a
// negative thickness behind a mirror, e.g. the −300 mm of
// library/geo/mirror-cassegrain-*.json). ImageZ, VertexZ and Magnification are
// therefore positions and ratios in one single z axis and can be compared
// directly with a real ray trace.
//
// The reflecting surface is handled by unfolding. With the surface sag written
// as z = c·r²/2 (c = 1/R, positive sag towards +z) the paraxial reflection law
// for the slope of the ray line is
//
//	u' = -u - 2·c·y ,
//
// so with ū = n·u  the interaction reads ū' = -ū - 2·n·c·y (the matrix
// [[1,0],[-2nc,-1]] of Mat2), while the gap behind the mirror is simply the
// signed thickness. The classical unfolding is the same statement in the frame
// of the propagating beam: negate the index of every medium behind the mirror
// (n → -n) and travel along -z; the reflection then becomes an ordinary
// refraction of power Φ = (n' - n)·c = -2·n·c with ū conserved (determinant +1),
// the curvature seen by the beam is ĉ = -c and the gap is t̂ = -t, all of which
// the signed-thickness layout and the medium list already carry. Every derived
// quantity (EFL, focal points, ImageZ, Magnification, pupils) is the standard
// linear extrapolation of the traced rays, so the first-order formulas of a
// refracting system hold verbatim for the unfolded chain — with two caveats for
// mirror systems, both of which follow from the direction reversal:
//
//   - the effective focal length keeps the sign of the surface radius
//     (a concave mirror seen from the object side has R < 0 and EFL = R/2 < 0,
//     its back focal point lies at -f in front of the vertex);
//   - FrontPrincipalZ and BackPrincipalZ are derived from EFL with the
//     refracting relations and are therefore not consistent with each other for
//     a single mirror (BackPrincipalZ = 0 is the physical one, the principal
//     plane of a mirror is its vertex); the physically verified quantities of a
//     mirror system are the focal points, ImageZ, Magnification and the angular
//     magnification;
//   - the Helmholtz–Lagrange invariant n(ȳu − yū) of the physical slopes
//     changes sign at every mirror; the invariant that is conserved, and the one
//     behind the relation m·γ = n/n' used by AngularMagnification, is the one
//     built from the slopes measured along the direction of propagation.
//
// A mirror is an interaction inside one medium: the media on its two sides must
// have the same index at the wavelength (otherwise the system is rejected; the
// real ray trace is more permissive and simply reflects). Everything else is
// unchanged: the first-order engine is defined for coaxial systems only, so a
// surface or a coordinate break with a decentre or a tilt about x or y still
// returns an error (a roll about z does not change the axis and is accepted);
// trace the real rays of such a system instead (SpotDiagram, RayFan and OPDGrid
// take the real chief ray as reference when the system has poses).
//
// Paraxial.Parity() converts the signed EFL into the customary focal length of
// the mirror convention (positive for a converging system): for a single
// concave mirror of R < 0 the engine reports EFL = R/2 < 0 while
// Parity()·EFL = |R|/2 is the usual f = R/2.
func (s *System) Paraxial(lambdaNm float64, objectZ float64) (*Paraxial, error) {
	if err := s.Validate(); err != nil {
		return nil, err
	}
	if err := s.requireCoaxial("paraxial analysis"); err != nil {
		return nil, err
	}
	mirrors := 0
	for i, sf := range s.Surfaces {
		nB := s.Mediums[i].N(lambdaNm)
		if sf.Mirror {
			nA := s.Mediums[i+1].N(lambdaNm)
			if nA != nB {
				return nil, fmt.Errorf("geo: paraxial analysis: surface %d is a mirror between "+
					"different media (%s, %s); a mirror reflects inside one medium, so the "+
					"media on both of its sides must be equal", i, s.Mediums[i].Name(), s.Mediums[i+1].Name())
			}
			mirrors++
		}
	}
	infinite := math.IsInf(objectZ, -1) || objectZ < -1e15
	if !infinite && objectZ >= s.VertexZ(0) {
		return nil, fmt.Errorf("geo: object at z=%g is not in front of surface 0 (z=%g)", objectZ, s.VertexZ(0))
	}
	els := s.paraxialElements(lambdaNm)
	p := &Paraxial{Lambda: lambdaNm, ObjectZ: objectZ, Infinite: infinite, els: els, Mirrors: mirrors}
	n := len(els)
	p.Index = make([]float64, n+1)
	p.VertexZ = make([]float64, n)
	for i := 0; i < n; i++ {
		p.Index[i] = els[i].nB
		p.VertexZ[i] = s.VertexZ(i)
	}
	p.Index[n] = els[n-1].nA
	p.StopIndex = s.StopIndex

	// Marginal ray of the axial field: y = 1 at surface 0.
	ys, _, uAfter := traceParaxial(els, 1, 0)
	yk := ys[n-1]
	uk := uAfter[n-1]
	if uk != 0 {
		p.EFL = -1 / uk
		p.BFL = -yk / uk
		p.BackFocalZ = p.VertexZ[n-1] + p.BFL
	} else {
		p.EFL = math.Inf(signOf(-1))
		p.BFL = math.Inf(signOf(-ys[n-1]))
		p.BackFocalZ = math.Inf(signOf(p.BFL))
	}
	// Front focal point from the time-reversed system.
	rels := reverseElements(els)
	rys, _, ruAfter := traceParaxial(rels, 1, 0)
	ryk := rys[n-1]
	ruk := ruAfter[n-1]
	if ruk != 0 {
		p.FFL = ryk / ruk // signed distance from surface 0 (z = 0) to the front focus
	} else {
		p.FFL = math.Inf(signOf(ryk))
	}
	p.FrontFocalZ = p.VertexZ[0] + p.FFL
	if !math.IsInf(p.EFL, 0) {
		p.BackPrincipalZ = p.BackFocalZ - p.EFL
		p.FrontPrincipalZ = p.FrontFocalZ + p.EFL
	} else {
		p.BackPrincipalZ = math.NaN()
		p.FrontPrincipalZ = math.NaN()
	}

	// Image plane and magnification.
	if infinite {
		p.ImageZ = p.BackFocalZ
		p.Magnification = 0
	} else {
		// Ray from the axial object point (y = 0 at objectZ) to y = 1 at surface 0.
		u0 := -1 / objectZ
		iys, _, iuAfter := traceParaxial(els, 1, u0)
		iyk := iys[n-1]
		iuk := iuAfter[n-1]
		if iuk != 0 {
			p.ImageZ = p.VertexZ[n-1] - iyk/iuk
		} else {
			p.ImageZ = math.Inf(signOf(-iyk))
		}
		// Parallel ray from the object point of height 1: its height at the
		// image plane is the transverse magnification.
		p.Magnification = yk + (p.ImageZ-p.VertexZ[n-1])*uk
	}

	// Pupils: images of the aperture stop in the object and image space.
	stop := s.StopIndex
	sr := s.Surfaces[stop].MaxRadius()
	p.StopRadius = sr
	p.EntrancePupilZ, p.EntrancePupilRadius, p.EntrancePupilMag = math.NaN(), math.NaN(), math.NaN()
	p.ExitPupilZ, p.ExitPupilRadius, p.ExitPupilMag = math.NaN(), math.NaN(), math.NaN()

	// Entrance pupil: conjugate of the stop plane through the elements before
	// the stop, evaluated in the object space.
	mIn := Identity()
	for i := 0; i < stop; i++ {
		mIn = els[i].surfaceMatrix().Mul(mIn)
		mIn = els[i].transferMatrix().Mul(mIn)
	}
	n0 := els[0].nB
	if mIn.A != 0 {
		L := -mIn.B * n0 / mIn.A // distance in front of surface 0
		p.EntrancePupilZ = s.VertexZ(0) - L
		p.EntrancePupilMag = mIn.A
		if sr > 0 && !math.IsInf(sr, 1) {
			p.EntrancePupilRadius = sr / math.Abs(mIn.A)
		}
	} else {
		// The stop is imaged at infinity (telecentric object space).
		p.EntrancePupilZ = math.Inf(-1)
	}
	// Exit pupil: conjugate of the stop plane through the elements after it.
	// Appending a gap L after the last surface left-multiplies the transfer
	// from the stop plane by [[1, L/nk],[0,1]], so the conjugate condition
	// B = 0 reads B + L·D/nk = 0: an appended gap is governed by the D entry
	// (the entrance pupil, whose gap is prepended, by the A entry). The stop
	// image sits at z = VertexZ(n−1) + L with semi-diameter
	// |A + L·C/nk|·StopRadius — the A entry of the conjugated transfer, which
	// equals |1/D| because the subsystem determinant is ±1.
	mOut := Identity()
	for i := stop; i < n; i++ {
		mOut = els[i].surfaceMatrix().Mul(mOut)
		mOut = els[i].transferMatrix().Mul(mOut)
	}
	nk := p.Index[n]
	if mOut.D != 0 {
		L := -mOut.B * nk / mOut.D
		p.ExitPupilZ = p.VertexZ[n-1] + L
		mag := mOut.A + L*mOut.C/nk
		p.ExitPupilMag = mag
		if sr > 0 && !math.IsInf(sr, 1) {
			p.ExitPupilRadius = sr * math.Abs(mag)
		}
	} else {
		// The stop images at infinity (telecentric image space).
		p.ExitPupilZ = math.Inf(1)
	}

	// F-number and image-side numerical aperture from the entrance pupil edge.
	if p.EntrancePupilRadius > 0 && !math.IsInf(p.EntrancePupilRadius, 1) {
		if uk != 0 {
			p.FNumber = p.EFL / (2 * p.EntrancePupilRadius)
		}
		uEdge := uk * p.EntrancePupilRadius
		p.ImageNA = math.Abs(nk * uEdge)
	} else {
		p.FNumber = math.NaN()
		p.ImageNA = math.NaN()
	}
	return p, nil
}

// MarginalRay returns the paraxial ray that starts at height y0 = radius at
// surface 0 and travels parallel to the axis (the marginal ray of the axial
// field; radius is normally the entrance pupil semi-diameter in mm).
func (p *Paraxial) MarginalRay(radius float64) *ParaxialRay {
	return p.rayFrom(radius, 0)
}

// ChiefRay returns the paraxial chief ray of the field whose tangent in the y
// direction is tanY: it passes through the centre of the aperture stop. For an
// object at infinity tanY = tan(field angle); for a finite object it is the
// height of the object point (tanY is then the object height in mm).
func (p *Paraxial) ChiefRay(tanY float64) *ParaxialRay {
	zStop := p.VertexZ[p.StopIndex]
	var y0, u0 float64
	if p.Infinite {
		u0 = -tanY
		y0 = zStop * tanY
	} else {
		u0 = (0 - tanY) / (zStop - p.ObjectZ)
		y0 = -u0 * zStop
	}
	return p.rayFrom(y0, u0)
}

// rayFrom traces a paraxial ray with the given starting values at surface 0.
func (p *Paraxial) rayFrom(y0, u0 float64) *ParaxialRay {
	ys, us, ua := traceParaxial(p.els, y0, u0)
	n := len(p.els)
	idx := make([]float64, n)
	for i := 0; i < n; i++ {
		idx[i] = p.els[i].nB
	}
	return &ParaxialRay{
		Y: ys, U: us, UAfter: ua, N: idx,
		VertexZ: append([]float64(nil), p.VertexZ...),
		Y0:      y0, U0: u0, UBar0: p.Index[0] * u0,
	}
}

// Lagrange returns the optical (Helmholtz–Lagrange) invariant
//
//	H = n (ȳ u - y ū)
//
// of the ray pair (marginal ray at height y0 = radius at surface 0, chief ray
// of field tangent tanY), with u the real slope dy/dz of the meridional ray.
// H is constant across all surfaces of an optical system; the invariant is used
// by the aberration routines and is checked by the test suite.
//
// For a system with mirrors the evaluation uses the physical slopes, and H then
// changes sign at every reflection (the ray comes back along -z). The conserved
// quantity of a folded chain is the one built from the slopes measured along the
// direction of propagation, w = Parity()·u, and it is that invariant that gives
// the imaging relation m·γ = n/n' behind AngularMagnification.
func (p *Paraxial) Lagrange(radius, tanY float64) float64 {
	m := p.MarginalRay(radius)
	c := p.ChiefRay(tanY)
	n := len(m.Y)
	if n == 0 {
		return 0
	}
	nn := p.els[0].nB
	return nn * (c.Y[0]*m.U[0] - m.Y[0]*c.U[0])
}

// LagrangeAt returns the invariant evaluated at surface i (incidence side); a
// constant value across i verifies the first-order invariants.
func (p *Paraxial) LagrangeAt(i int, radius, tanY float64) float64 {
	m := p.MarginalRay(radius)
	c := p.ChiefRay(tanY)
	return p.els[i].nB * (c.Y[i]*m.U[i] - m.Y[i]*c.U[i])
}

// Matrix returns the ray transfer matrix of the whole system in (y, ū) form,
// from the plane of surface 0 to the plane of the last surface. Its determinant
// is (-1)^(number of mirrors) · n0/nk.
func (p *Paraxial) Matrix() Mat2 {
	m := Identity()
	for i := range p.els {
		m = p.els[i].surfaceMatrix().Mul(m)
		m = p.els[i].transferMatrix().Mul(m)
	}
	return m
}

// SubMatrix returns the ray transfer matrix from surface 'from' to surface
// 'to' (both boundaries included as planes, from <= to).
func (p *Paraxial) SubMatrix(from, to int) (Mat2, error) {
	if from < 0 || to >= len(p.els) || from > to {
		return Mat2{}, fmt.Errorf("geo: invalid matrix range [%d,%d]", from, to)
	}
	m := Identity()
	for i := from; i <= to; i++ {
		m = p.els[i].surfaceMatrix().Mul(m)
		m = p.els[i].transferMatrix().Mul(m)
	}
	return m, nil
}

// Parity returns the sign of the propagation direction of the emergent beam
// relative to the incident one: +1 when the light leaves the system travelling
// the way it entered (no mirror or an even number of them), -1 after an odd
// number of reflections. It is the factor that converts a real slope u = dy/dz
// of the emergent ray into a slope measured along the direction of propagation,
// w = Parity()·u (see AngularMagnification); profiles and step positions are
// unaffected by it.
func (p *Paraxial) Parity() float64 {
	if p.Mirrors%2 != 0 {
		return -1
	}
	return 1
}

// afocalTol is the residual power (1/mm) below which a system counts as
// afocal: |C| <= afocalTol, i.e. an effective focal length beyond 1e12 mm.
// The composition of the elementary matrices of a telescope leaves a rounding
// residue of the order of 1e-17/mm, so the test must not demand an exact zero.
const afocalTol = 1e-12

// Afocal reports whether the system has no first-order power: a beam entering
// parallel to the axis leaves parallel to the axis (C of the ray transfer
// matrix vanishes, |C| <= 1e-12 /mm). An afocal system has an infinite EFL; its
// finite first-order quantity is the angular magnification
// (System.AngularMagnification).
func (p *Paraxial) Afocal() bool {
	return math.Abs(p.Matrix().C) <= afocalTol
}

// AngularMagnification returns the first-order angular magnification of the
// paraxial system:
//
//   - afocal system (C = 0): Γ = w'/w, the ratio of the slope of the emergent
//     ray to the slope of the incident ray of the field whose tangent is
//     fieldTan, both slopes measured along the direction of propagation
//     (w = dy/dz in the frame that travels with the beam, see Paraxial.Parity).
//     The chief ray is used: for an object at infinity the parallel ray at angle
//     fieldTan that passes through the centre of the aperture stop, for a finite
//     object the ray from the object point at height fieldTan through the stop
//     centre. When C = 0 the system maps ū ↦ D·ū for every ray, so
//     Γ = Parity()·D·n0/nk does not depend on the chosen ray; a telescope is the
//     standard example (Kepler f1 = 100, f2 = 20, d = 120 gives Γ = -5).
//   - imaging system (C ≠ 0): γ = n0/(nk·m), the ratio w'/w of the axial ray,
//     where m is the transverse magnification (System.Paraxial returns it in
//     Magnification) and n0, nk the object- and image-side indices. This is the
//     Lagrange-invariant relation m·γ = n0/nk written for the ray that leaves
//     the axial object point; it is defined for a finite object only and is
//     independent of fieldTan.
//
// The sign convention of the propagation-frame slopes makes a plane fold mirror
// transparent to angles (Γ = +1) and gives the classical negative Γ for an
// inverting Kepler telescope, so the sign of the magnification is the sign of
// the image.
func (p *Paraxial) AngularMagnification(fieldTan float64) (float64, error) {
	n := len(p.els)
	if n == 0 {
		return 0, fmt.Errorf("geo: angular magnification: empty system")
	}
	if p.Afocal() {
		c := p.ChiefRay(fieldTan)
		if c.U0 == 0 {
			return 0, fmt.Errorf("geo: angular magnification: the reference ray of the field "+
				"sits on the axis (field tangent %g) and has no angle to compare", fieldTan)
		}
		return p.Parity() * c.UAfter[n-1] / c.U0, nil
	}
	if p.Infinite {
		return 0, fmt.Errorf("geo: angular magnification: a system with power (EFL=%g) has no "+
			"angular magnification for an object at infinity; use a finite object or an afocal system", p.EFL)
	}
	if p.Magnification == 0 {
		return 0, fmt.Errorf("geo: angular magnification: transverse magnification is zero")
	}
	n0 := p.Index[0]
	nk := p.Index[len(p.Index)-1]
	return n0 / (nk * p.Magnification), nil
}

// AngularMagnification returns the first-order angular magnification of the
// system at wavelength lambdaNm for an object at axial position objectZ, with
// the reference field tangent fieldTan (mm of object height for a finite object,
// tangent of the field angle for an object at infinity — the same convention as
// ChiefRay). See Paraxial.AngularMagnification for the definition: it is the
// ratio of the emergent to the incident slope of the chief ray for an afocal
// system, and the axial-ray ratio γ = n/(n'·m) for an imaging system.
//
// Examples: a Kepler telescope of f1 = 100 mm, f2 = 20 mm and d = f1 + f2
// returns -5 (inverted), the Galilean counterpart +5 (erect), and a single
// plane fold mirror +1.
func (s *System) AngularMagnification(objectZ, fieldTan, lambdaNm float64) (float64, error) {
	p, err := s.Paraxial(lambdaNm, objectZ)
	if err != nil {
		return 0, err
	}
	return p.AngularMagnification(fieldTan)
}

// ReciprocalFocus returns the axial position at which a point source placed at
// axial position zObject (mm) is imaged, for the same wavelength. The system may
// contain mirrors (see System.Paraxial); a folded chain simply extrapolates the
// emergent ray line back to the axis in the physical z frame.
func (p *Paraxial) Conjugate(zObject float64) float64 {
	n := len(p.els)
	if n == 0 {
		return math.NaN()
	}
	u0 := -1 / zObject
	ys, _, ua := traceParaxial(p.els, 1, u0)
	if ua[n-1] == 0 {
		return math.Inf(signOf(-ys[n-1]))
	}
	return p.VertexZ[n-1] - ys[n-1]/ua[n-1]
}

// signOf returns -1 for x < 0 and +1 otherwise (0 gives +1).
func signOf(x float64) int {
	if x < 0 {
		return -1
	}
	return 1
}

// String implements fmt.Stringer.
func (p *Paraxial) String() string {
	return fmt.Sprintf("paraxial[λ=%gnm EFL=%g BFL=%g FFL=%g epZ=%g epR=%g F/%g NA=%g]",
		p.Lambda, p.EFL, p.BFL, p.FFL, p.EntrancePupilZ, p.EntrancePupilRadius, p.FNumber, p.ImageNA)
}
