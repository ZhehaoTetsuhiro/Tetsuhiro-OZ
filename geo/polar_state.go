package geo

import (
	"errors"
	"fmt"
	"math"
)

// cvec is a complex electric-field amplitude in the global frame of package geo.
type cvec [3]complex128

// polReal returns the real vector v scaled by the complex factor a.
func polReal(v Vec3, a complex128) cvec {
	return cvec{a * complex(v.X, 0), a * complex(v.Y, 0), a * complex(v.Z, 0)}
}

// plus returns e + f.
func (e cvec) plus(f cvec) cvec {
	return cvec{e[0] + f[0], e[1] + f[1], e[2] + f[2]}
}

// dot returns the bilinear product e·v with a real vector v (no conjugation).
func (e cvec) dot(v Vec3) complex128 {
	return e[0]*complex(v.X, 0) + e[1]*complex(v.Y, 0) + e[2]*complex(v.Z, 0)
}

// intensity returns |e|² = Σ|e_k|².
func (e cvec) intensity() float64 {
	s := 0.0
	for _, c := range e {
		s += real(c)*real(c) + imag(c)*imag(c)
	}
	return s
}

// TransverseFrame returns the right-handed transverse basis (ex, ey) of a ray
// travelling along the direction d:
//
//	ex = ŷ×d / |ŷ×d| ,   ey = d×ex ,   so that ex×ey = d .
//
// On axis (d = +ẑ) this is the global (x, y) frame of optics/polar.go: with the
// beam travelling towards the viewer, ex points right and ey up. The rule is
// singular for d ∥ ±ŷ, where ex = x̂ is used instead. Jones vectors of this
// package are always expressed in this frame, at the input and at the output.
func TransverseFrame(d Vec3) (ex, ey Vec3) {
	u, ok := d.Unit()
	if !ok {
		u = AxisZ
	}
	ex = AxisY.Cross(u)
	if n := ex.Norm(); n > 1e-12 {
		ex = ex.Div(n)
	} else {
		ex = AxisX
	}
	ey = u.Cross(ex)
	return ex, ey
}

// polBasis is the local s/p basis of one interaction.
type polBasis struct {
	s, pIn, pOut Vec3
}

// polInteract propagates the transverse field e through one interaction at a
// surface with unit normal nrm (either orientation). din and dout are the unit
// incident and outgoing directions. reflected selects a reflected wave (dout is
// the mirror direction); cs and cp are the s and p coefficients. The returned
// field is transverse to dout.
//
// The s/p basis is rebuilt from din and nrm at every call, so the rotation of the
// plane of incidence from one surface to the next is carried exactly by the
// vector algebra; no basis is stored between surfaces.
func polInteract(din, nrm, dout Vec3, reflected bool, cs, cp complex128, e cvec) (cvec, polBasis) {
	s := din.Cross(nrm)
	if n := s.Norm(); n > 1e-15 {
		s = s.Div(n)
	} else {
		// Normal incidence: the coefficients of s and p coincide there, so any
		// direction transverse to din serves as s.
		s, _ = TransverseFrame(din)
	}
	pIn := s.Cross(din)
	var pOut Vec3
	if reflected {
		pOut = dout.Cross(s) // = −s×dout, the unified convention for reflection
	} else {
		pOut = s.Cross(dout)
	}
	out := polReal(s, cs*e.dot(s)).plus(polReal(pOut, cp*e.dot(pIn)))
	return out, polBasis{s: s, pIn: pIn, pOut: pOut}
}

// interactionJones returns the 2×2 matrix of one interaction, from the transverse
// frame of the incident ray (TransverseFrame(din)) to that of the outgoing ray
// (TransverseFrame(dout)): J[i][j] = eout_i·M·ein_j with M = cs·s sᵀ + cp·pOut pInᵀ.
func interactionJones(din, dout Vec3, b polBasis, cs, cp complex128) JonesMat {
	var ein, eout [2]Vec3
	ein[0], ein[1] = TransverseFrame(din)
	eout[0], eout[1] = TransverseFrame(dout)
	var j JonesMat
	for i := 0; i < 2; i++ {
		for k := 0; k < 2; k++ {
			j[i][k] = cs*complex(eout[i].Dot(b.s)*b.s.Dot(ein[k]), 0) +
				cp*complex(eout[i].Dot(b.pOut)*b.pIn.Dot(ein[k]), 0)
		}
	}
	return j
}

// applyElement applies the Jones matrix j to the field e, expressed in the
// transverse frame of the ray direction d.
func applyElement(e cvec, d Vec3, j JonesMat) cvec {
	ex, ey := TransverseFrame(d)
	a, b := e.dot(ex), e.dot(ey)
	return polReal(ex, j[0][0]*a+j[0][1]*b).plus(polReal(ey, j[1][0]*a+j[1][1]*b))
}

// PolStep is one interaction of a polarized trace.
type PolStep struct {
	Surface int     // index of the surface in System.Surfaces
	Label   string  // surface label
	Mirror  bool    // reflective interaction
	NBefore float64 // refractive index in front of the surface
	NAfter  float64 // refractive index behind the surface
	DirIn   Vec3    // unit incident direction
	Normal  Vec3    // unit surface normal at the hit point, +z oriented
	DirOut  Vec3    // unit outgoing direction (reflected or refracted)
	CosI    float64 // |cos| of the angle of incidence
	CosT    float64 // |cos| of the angle of the outgoing wave
	S       Vec3    // s unit vector, normal to the plane of incidence
	PIn     Vec3    // p vector of the incident wave (S×DirIn)
	POut    Vec3    // p vector of the outgoing wave
	Amp     Amplitudes
	// J is the 2×2 matrix of this interaction from the transverse frame of the
	// incident ray to that of the outgoing ray (elements not included).
	J JonesMat
	// Field is the electric field (global frame) just after this interaction,
	// including any element placed after this surface.
	Field [3]complex128
}

// PolElementAt places a 2×2 Jones element after the interaction with a surface.
type PolElementAt struct {
	AfterSurface int
	J            JonesMat
}

// PolOptions configures a polarized trace. A nil pointer or the zero value uses
// FresnelStep at every surface and no elements.
type PolOptions struct {
	// Amplitude[i], when non-nil, replaces the amplitude model of surface i (for
	// example CoatedAmplitudes). Missing or nil entries use FresnelStep.
	Amplitude []AmplitudeFunc
	// Elements act in the transverse frame of the ray after the interaction with
	// surface AfterSurface, in slice order.
	Elements []PolElementAt
}

// PolTrace is the outcome of TracePolarized.
type PolTrace struct {
	Trace    TraceResult   // the geometrical trace (System.Trace, unmodified)
	Steps    []PolStep     // one entry per completed interaction
	InFrame  [2]Vec3       // input transverse basis (ex, ey)
	OutFrame [2]Vec3       // output transverse basis (ex, ey)
	Field    [3]complex128 // output electric field, global frame
	Jones    [2]complex128 // output Jones vector in OutFrame
	// Matrix is the system Jones matrix accumulated surface by surface (product
	// of the interaction and element matrices along the path), so that
	// Jones = Matrix·in. It is checked against the 3D field in the tests.
	Matrix JonesMat
	// Transmittance is P_out/P_in with P = n·cosθ·|E|² evaluated at the entrance
	// and at the exit surface element. It contains the Fresnel and the element
	// losses of the single traced path; it does not model the change of ray-tube
	// area, which a single ray does not have.
	Transmittance float64
	// OPLPhase = 2π·OPL/λ is the propagation phase of the optical path length.
	// It is reported but not applied to Jones (see the limitations in the
	// report and the package documentation).
	OPLPhase float64
	Ellipse  Ellipse // polarization ellipse of Jones
}

// TracePolarized traces r through the system and propagates the polarization
// state of the ray. The input Jones vector in is given in the input frame
// TransverseFrame(r.Dir). At each completed interaction the 3D field is
// multiplied by the amplitude coefficients of that surface in the s/p basis of
// that very interaction; the trace itself is the unmodified System.Trace.
//
// On error (invalid ray, vignetting, TIR, amplitude model failure) the partial
// result is returned together with the error.
func (s *System) TracePolarized(r Ray, in [2]complex128, opt *PolOptions) (PolTrace, error) {
	var out PolTrace
	if opt == nil {
		opt = &PolOptions{}
	}
	if !r.Valid() {
		return out, errors.New("geo: polarization: invalid ray")
	}
	out.InFrame[0], out.InFrame[1] = TransverseFrame(r.Dir)
	e := polReal(out.InFrame[0], in[0]).plus(polReal(out.InFrame[1], in[1]))
	res := s.Trace(r)
	out.Trace = res
	lambda := r.WavelengthNm
	out.OPLPhase = 2 * math.Pi * res.Ray.OPL / (lambda * 1e-6)
	if len(res.Steps) == 0 {
		return out, res.Error()
	}
	st0 := res.Steps[0]
	pin := st0.NBefore * math.Abs(r.Dir.Dot(st0.Normal)) * e.intensity()
	acc := JonesMat{{1, 0}, {0, 1}}

	for _, st := range res.Steps {
		if st.Status != StatusOK {
			break
		}
		amp, err := polAmplitude(opt, st, lambda)
		if err != nil {
			return out, fmt.Errorf("geo: polarization at surface %d (%s): %w", st.Surface, st.Label, err)
		}
		cs, cp := amp.TS, amp.TP
		if st.Mirror {
			cs, cp = amp.RS, amp.RP
		}
		next, basis := polInteract(st.DirIn, st.Normal, st.DirOut, st.Mirror, cs, cp, e)
		e = next
		jStep := interactionJones(st.DirIn, st.DirOut, basis, cs, cp)
		acc = jStep.Mul(acc)
		for _, el := range opt.Elements {
			if el.AfterSurface == st.Surface {
				e = applyElement(e, st.DirOut, el.J)
				acc = el.J.Mul(acc)
			}
		}
		out.Steps = append(out.Steps, PolStep{
			Surface: st.Surface, Label: st.Label, Mirror: st.Mirror,
			NBefore: st.NBefore, NAfter: st.NAfter,
			DirIn: st.DirIn, Normal: st.Normal, DirOut: st.DirOut,
			CosI: math.Abs(st.DirIn.Dot(st.Normal)), CosT: math.Abs(st.DirOut.Dot(st.Normal)),
			S: basis.s, PIn: basis.pIn, POut: basis.pOut,
			Amp:   amp,
			J:     jStep,
			Field: [3]complex128(e),
		})
	}
	if !res.OK() {
		return out, res.Error()
	}

	last := out.Steps[len(out.Steps)-1]
	nOut := last.NAfter
	if last.Mirror {
		nOut = last.NBefore
	}
	cOut := math.Abs(res.Ray.Dir.Dot(last.Normal))
	if pin > 0 {
		out.Transmittance = nOut * cOut * e.intensity() / pin
	}
	out.Field = [3]complex128(e)
	out.OutFrame[0], out.OutFrame[1] = TransverseFrame(res.Ray.Dir)
	out.Jones = [2]complex128{e.dot(out.OutFrame[0]), e.dot(out.OutFrame[1])}
	out.Matrix = acc
	out.Ellipse = EllipseOf(out.Jones[0], out.Jones[1])
	return out, nil
}

// JonesMatrix returns the 2×2 Jones matrix of the path of r. Column k is the
// output Jones vector (in the output frame) produced by the input basis vector
// k (in the input frame). It traces the ray twice.
func (s *System) JonesMatrix(r Ray, opt *PolOptions) (JonesMat, error) {
	a, err := s.TracePolarized(r, [2]complex128{1, 0}, opt)
	if err != nil {
		return JonesMat{}, err
	}
	b, err := s.TracePolarized(r, [2]complex128{0, 1}, opt)
	if err != nil {
		return JonesMat{}, err
	}
	return JonesMat{{a.Jones[0], b.Jones[0]}, {a.Jones[1], b.Jones[1]}}, nil
}

// polAmplitude returns the amplitude coefficients of one step from the options.
func polAmplitude(opt *PolOptions, st Step, lambda float64) (Amplitudes, error) {
	f := AmplitudeFunc(FresnelStep)
	if st.Surface < len(opt.Amplitude) && opt.Amplitude[st.Surface] != nil {
		f = opt.Amplitude[st.Surface]
	}
	return f(st, lambda)
}
