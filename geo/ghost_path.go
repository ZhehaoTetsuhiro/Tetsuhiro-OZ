package geo

import (
	"fmt"
	"math"
)

// ---------------------------------------------------------------------------
// Double-reflection (ghost) paths
// ---------------------------------------------------------------------------
//
// A ghost path is an ordinary ray path with one extra pair of Fresnel
// reflections: the ray crosses the surfaces in front of I, reflects at surface
// I instead of transmitting, propagates backwards (towards −z) through the
// surfaces in front of it, reflects at surface J < I, and then propagates
// forwards again through every surface behind J — surface I a second time — to
// the image plane. Written as a node list the path is
//
//	0 … I−1 (refract, +z) , I (reflect) , I−1 … J+1 (refract, −z) ,
//	J (reflect) , J+1 … n−1 (refract, +z) , image plane
//
// For a cemented doublet (0 = front, 1 = cement, 2 = rear) the pair
// I = 2, J = 1 is the classical two-inner-surface ghost: reflect at the rear
// surface, travel back to the cement, reflect there and leave through the rear
// surface again. Every pair 0 <= J < I <= n−1 is a physical double-reflection
// path; GhostPairs enumerates them in a deterministic order (I ascending, then
// J ascending).
//
// Sign conventions of the reverse pass — the whole content of a ghost trace:
//
//   - The direction vector is never mirrored: after the reflection at I it
//     simply acquires a negative z component. Every axial transfer uses the
//     signed vertex distance z_next − z_current, which is negative while the
//     ray runs backwards; no coordinate flipping is needed anywhere.
//   - The refraction at a surface crossed towards −z takes its media from the
//     side the ray arrives from: crossing surface m backwards, the incident
//     medium is Mediums[m+1] and the medium beyond it is Mediums[m] — the
//     opposite assignment of the forward pass. Snell's law itself is symmetric
//     in the two directions (transmitting a ray backwards through a surface
//     undoes the forward refraction); geo.Refract orients the surface normal
//     against the incoming ray on its own, so this medium swap is the only
//     change the vector tracer needs.
//   - Paraxially (see ghostParaxialRecurse) the same statement in (y, u) form —
//     u = dy/dz and η = n_in/n_out the ratio of the incident to the exit index
//     — is
//
//	refraction:  u' = η u + (η − 1) c y   (either direction of travel),
//	reflection:  u' = −u − 2 c y          (either direction of travel),
//	transfer:    y' = y + (z' − z) u      (either direction of travel),
//
//     with c the vertex curvature of the surface. For the backwards pass
//     η = nA/nB instead of nB/nA; for a plane (c = 0) the transfer and the
//     refraction leave u alone and the reflection simply flips it.
//
// Documented limitations of the model (see also ghost_analysis.go):
//
//   - The reflectance is a real normal-incidence Fresnel value (or a
//     user-supplied R): no coating phase, no angle dependence of the Fresnel
//     coefficients and no polarization split.
//   - Exactly two reflections: higher orders (4, 6, … reflections) and any
//     scattering are out of scope.
//   - The reflection sites must be refracting surfaces: mirror surfaces are
//     rejected (their nominal action is already a reflection).

// GhostPair identifies one double-reflection (ghost) path by its two
// reflecting surfaces. Light travels in +z when it reaches I, reflects there,
// travels in −z to J (< I), reflects there and travels in +z again to the
// image plane.
type GhostPair struct {
	I int // surface of the first reflection (the one further along the axis)
	J int // surface of the second reflection (J < I)
}

// String implements fmt.Stringer.
func (p GhostPair) String() string { return fmt.Sprintf("reflect(I=%d, J=%d)", p.I, p.J) }

// GhostAction is what happens to the ray at one node of a ghost path.
type GhostAction int

const (
	// GhostRefract: the ray crosses the surface following Snell's law with the
	// media on the two sides of the surface.
	GhostRefract GhostAction = iota
	// GhostReflect: the ray reflects (a Fresnel ghost reflection).
	GhostReflect
)

// String implements fmt.Stringer.
func (a GhostAction) String() string {
	switch a {
	case GhostRefract:
		return "refract"
	case GhostReflect:
		return "reflect"
	}
	return fmt.Sprintf("ghost-action(%d)", int(a))
}

// GhostNode is one interaction of a ghost path.
type GhostNode struct {
	Surface  int         // index of the surface
	Action   GhostAction // refraction or (ghost) reflection
	Backward bool        // true when the ray reaches the node travelling towards −z
	// Medium is the index of the gap the ray is in just before the node, i.e.
	// an index into System.Mediums. Surface belongs to Mediums[Surface] and
	// Mediums[Surface+1]; the incident medium is Medium (the side the ray comes
	// from), the exit medium is the other one. Reflections leave it unchanged.
	Medium int
}

// String implements fmt.Stringer.
func (n GhostNode) String() string {
	dir := "+z"
	if n.Backward {
		dir = "-z"
	}
	return fmt.Sprintf("s%d:%s(%s,gap%d)", n.Surface, n.Action, dir, n.Medium)
}

// GhostPath is the complete node sequence of one double-reflection path.
type GhostPath struct {
	Pair  GhostPair
	Nodes []GhostNode
}

// String implements fmt.Stringer.
func (p GhostPath) String() string {
	s := "ghost[" + p.Pair.String()
	for _, n := range p.Nodes {
		s += " " + n.String()
	}
	return s + "]"
}

// GhostPairs returns every double-reflection pair of a system with n surfaces,
// in a deterministic order: I ascending, then J ascending. The count is
// n(n−1)/2 and the result is nil for n < 2.
func GhostPairs(n int) []GhostPair {
	if n < 2 {
		return nil
	}
	out := make([]GhostPair, 0, n*(n-1)/2)
	for i := 1; i < n; i++ {
		for j := 0; j < i; j++ {
			out = append(out, GhostPair{I: i, J: j})
		}
	}
	return out
}

// GhostNodes returns the node sequence of the pair p of a system with n
// surfaces. It only depends on the geometry of the sequence (the surfaces are
// visited in the order of the path); the media of the nodes are the gap
// indices described in GhostNode.
func GhostNodes(n int, p GhostPair) ([]GhostNode, error) {
	if p.I < 0 || p.I >= n || p.J < 0 || p.J >= p.I {
		return nil, fmt.Errorf("geo: ghost pair %s is out of range for a system with %d surfaces", p, n)
	}
	nodes := make([]GhostNode, 0, 2*p.I+n-2*p.J)
	// Forward pass: the ray crosses every surface in front of the first
	// reflection, which happens at I.
	for m := 0; m < p.I; m++ {
		nodes = append(nodes, GhostNode{Surface: m, Action: GhostRefract, Medium: m})
	}
	nodes = append(nodes, GhostNode{Surface: p.I, Action: GhostReflect, Medium: p.I})
	// Reverse pass: the surfaces strictly between I and J are crossed towards
	// −z, so the ray is in the gap on the +z side of the surface (m+1).
	for m := p.I - 1; m > p.J; m-- {
		nodes = append(nodes, GhostNode{Surface: m, Action: GhostRefract, Backward: true, Medium: m + 1})
	}
	nodes = append(nodes, GhostNode{Surface: p.J, Action: GhostReflect, Backward: true, Medium: p.J + 1})
	// Second forward pass: every surface behind J, the first reflection site
	// among them.
	for m := p.J + 1; m < n; m++ {
		nodes = append(nodes, GhostNode{Surface: m, Action: GhostRefract, Medium: m})
	}
	return nodes, nil
}

// GhostPathOf returns the node sequence of the pair p for this system,
// rejecting pairs that involve elements the double-reflection model cannot
// describe: coordinate-break elements (no profile) and mirrors (their nominal
// action is already a reflection, so they cannot be a refracting node nor a
// Fresnel reflection site).
func (s *System) GhostPathOf(p GhostPair) (GhostPath, error) {
	if err := s.Validate(); err != nil {
		return GhostPath{}, err
	}
	nodes, err := GhostNodes(len(s.Surfaces), p)
	if err != nil {
		return GhostPath{}, err
	}
	for _, nd := range nodes {
		sf := &s.Surfaces[nd.Surface]
		if sf.Break {
			return GhostPath{}, fmt.Errorf("geo: %s: element %d is a coordinate break", p, nd.Surface)
		}
		if sf.Mirror {
			return GhostPath{}, fmt.Errorf("geo: %s: surface %d is a mirror; the ghost model covers "+
				"Fresnel reflections at refracting surfaces only", p, nd.Surface)
		}
	}
	return GhostPath{Pair: p, Nodes: nodes}, nil
}

// ---------------------------------------------------------------------------
// Relative energy of a ghost path
// ---------------------------------------------------------------------------

// FresnelAmplitude returns the normal-incidence amplitude reflectance
//
//	r = (n1 − n2)/(n1 + n2)
//
// of the interface between two media of refractive index n1 and n2. The sign
// of r depends on the direction of travel through the interface; the intensity
// reflectance (FresnelReflectance) does not.
func FresnelAmplitude(n1, n2 float64) float64 { return (n1 - n2) / (n1 + n2) }

// FresnelReflectance returns the normal-incidence intensity (power)
// reflectance
//
//	R = r² = ((n1 − n2)/(n1 + n2))² ,  r = FresnelAmplitude(n1, n2) .
//
// R is symmetric in n1 and n2, so it does not matter from which side a ray
// arrives at a surface: a ghost reflection at a surface between air and glass
// contributes the same reflectance whichever pass of the path hits it. NaN is
// returned for n1 = −n2 = 0, the only degenerate case.
func FresnelReflectance(n1, n2 float64) float64 {
	if n1+n2 == 0 {
		return math.NaN()
	}
	r := (n1 - n2) / (n1 + n2)
	return r * r
}

// GhostReflectance returns the intensity reflectance R of surface i at the
// given wavelength. When override[i] exists and is not NaN it is used as it is
// (a supplied per-surface R, e.g. from a coating run); otherwise R is the
// normal-incidence Fresnel reflectance of the two media of the surface,
// FresnelReflectance(Mediums[i], Mediums[i+1]).
//
// The reflectance is a property of the interface alone: it is the same for the
// forwards and the backwards pass of a ghost path, so a path contributes
// exactly one R per reflection site.
func (s *System) GhostReflectance(i int, lambdaNm float64, override []float64) (float64, error) {
	if i < 0 || i >= len(s.Surfaces) {
		return 0, fmt.Errorf("geo: ghost reflectance: surface %d out of range [0,%d]", i, len(s.Surfaces)-1)
	}
	if i < len(override) && !math.IsNaN(override[i]) {
		r := override[i]
		if math.IsInf(r, 0) || r < 0 || r > 1 {
			return 0, fmt.Errorf("geo: ghost reflectance: surface %d: supplied R = %g is outside [0,1]", i, r)
		}
		return r, nil
	}
	n1 := s.Mediums[i].N(lambdaNm)
	n2 := s.Mediums[i+1].N(lambdaNm)
	if !finitePositive(n1) || !finitePositive(n2) {
		return 0, fmt.Errorf("geo: ghost reflectance: surface %d: refractive indices %g, %g are not positive and finite", i, n1, n2)
	}
	return FresnelReflectance(n1, n2), nil
}

// GhostEnergy returns the relative energy of the ghost path p at the given
// wavelength together with the two reflectances it is made of:
//
//	Energy = R_I · R_J
//
// with R_I the reflectance of the first reflection site and R_J that of the
// second (see GhostReflectance for how each R is obtained). It is a purely
// algebraic estimate — the transmission of the surfaces the ray crosses on the
// way is not included, which is the standard relative ordering of ghost paths.
func (s *System) GhostEnergy(p GhostPair, lambdaNm float64, override []float64) (energy, reflectI, reflectJ float64, err error) {
	if err := s.Validate(); err != nil {
		return 0, 0, 0, err
	}
	if p.I >= len(s.Surfaces) || p.J < 0 || p.J >= p.I {
		return 0, 0, 0, fmt.Errorf("geo: ghost energy: pair %s out of range", p)
	}
	reflectI, err = s.GhostReflectance(p.I, lambdaNm, override)
	if err != nil {
		return 0, 0, 0, err
	}
	reflectJ, err = s.GhostReflectance(p.J, lambdaNm, override)
	if err != nil {
		return 0, 0, 0, err
	}
	return reflectI * reflectJ, reflectI, reflectJ, nil
}
