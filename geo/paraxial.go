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
//
// and the imaging condition between two conjugate planes is B = 0, with
// transverse magnification A. The effective focal length of a system is
// EFL = -n'/C (n' = image-side index), and the matrix of a thin lens in air of
// focal length f is [[1,0],[-1/f,1]].
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
type paraxialElement struct {
	c      float64 // curvature 1/mm
	nB, nA float64 // refractive index before/after
	t      float64 // distance to the next surface, mm
}

// refractMatrix returns the refraction matrix [[1,0],[-Φ,1]] in (y, ū).
func (e paraxialElement) refractMatrix() Mat2 {
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
			c:  curvatureOf(s.Surfaces[i].Profile),
			nB: s.Mediums[i].N(lambdaNm),
			nA: s.Mediums[i+1].N(lambdaNm),
			t:  t,
		}
	}
	return els
}

// reverseElements returns the element list of the time-reversed system.
func reverseElements(els []paraxialElement) []paraxialElement {
	n := len(els)
	out := make([]paraxialElement, n)
	for j := 0; j < n; j++ {
		e := els[n-1-j]
		out[j] = paraxialElement{c: -e.c, nB: e.nA, nA: e.nB}
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
//	transfer:   y' = y + t u' .
//
// It returns the heights and slopes at the incidence side of every surface
// (ys, us, with us the real slope u = dy/dz) and the real slopes after each
// surface (usAfter).
func traceParaxial(els []paraxialElement, y, u float64) (ys, us, usAfter []float64) {
	n := len(els)
	ys = make([]float64, n)
	us = make([]float64, n)
	usAfter = make([]float64, n)
	for i := range els {
		ys[i], us[i] = y, u
		phi := (els[i].nA - els[i].nB) * els[i].c
		u2 := (els[i].nB*u - y*phi) / els[i].nA
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

	// EntrancePupilZ / ExitPupilZ are the axial positions of the images of the
	// aperture stop in the object and image space (mm). The pupil planeses can
	// be virtual (inside the system). EntrancePupilRadius/ExitPupilRadius are
	// the corresponding semi-diameters (mm).
	EntrancePupilZ      float64
	ExitPupilZ          float64
	EntrancePupilRadius float64
	ExitPupilRadius     float64
	StopRadius          float64 // semi-diameter of the aperture stop (mm)

	FNumber float64 // EFL / (2·EntrancePupilRadius)
	ImageNA float64 // n'·sin(u') of the marginal ray at the entrance pupil edge

	// PrincipalMag is A of the transfer from the object space to the stop
	// plane (the entrance pupil magnification).
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
// Mirrors are rejected: the first-order treatment of a reflecting surface
// requires unfolding the system, which is not part of this revision. Use real
// ray tracing for mirror systems.
func (s *System) Paraxial(lambdaNm float64, objectZ float64) (*Paraxial, error) {
	if err := s.Validate(); err != nil {
		return nil, err
	}
	for i, sf := range s.Surfaces {
		if sf.Mirror {
			return nil, fmt.Errorf("geo: paraxial analysis: surface %d is a mirror; "+
				"the paraxial engine only handles refraction (trace real rays instead)", i)
		}
	}
	infinite := math.IsInf(objectZ, -1) || objectZ < -1e15
	if !infinite && objectZ >= s.VertexZ(0) {
		return nil, fmt.Errorf("geo: object at z=%g is not in front of surface 0 (z=%g)", objectZ, s.VertexZ(0))
	}
	els := s.paraxialElements(lambdaNm)
	p := &Paraxial{Lambda: lambdaNm, ObjectZ: objectZ, Infinite: infinite, els: els}
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
		mIn = els[i].refractMatrix().Mul(mIn)
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
	mOut := Identity()
	for i := stop; i < n; i++ {
		mOut = els[i].refractMatrix().Mul(mOut)
		mOut = els[i].transferMatrix().Mul(mOut)
	}
	nk := p.Index[n]
	if mOut.A != 0 {
		L := -mOut.B * nk / mOut.A
		p.ExitPupilZ = p.VertexZ[n-1] + L
		p.ExitPupilMag = mOut.A
		if sr > 0 && !math.IsInf(sr, 1) {
			p.ExitPupilRadius = sr / math.Abs(mOut.A)
		}
	} else {
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
// of field tangent tanY). H is constant across all surfaces of an optical
// system; the invariant is used by the aberration routines and is checked by
// the test suite.
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
// from the plane of surface 0 to the plane of the last surface.
func (p *Paraxial) Matrix() Mat2 {
	m := Identity()
	for i := range p.els {
		m = p.els[i].refractMatrix().Mul(m)
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
		m = p.els[i].refractMatrix().Mul(m)
		m = p.els[i].transferMatrix().Mul(m)
	}
	return m, nil
}

// ReciprocalFocus returns the axial position at which a point source placed at
// axial position zObject (mm) is imaged, for the same wavelength.
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
