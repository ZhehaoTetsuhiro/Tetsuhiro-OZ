package geo

import (
	"fmt"
	"math"
)

// ---------------------------------------------------------------------------
// Replaying a ghost path with real rays
// ---------------------------------------------------------------------------
//
// ReplayGhost is a minimal sequential replayer: it walks the node list of a
// GhostPath and applies the same intersection (Surface.Intersect), normal
// (Surface.Normal), Snell (Refract) and mirror (Reflect) machinery the forward
// Trace uses, but in the path order and with the media of the side the ray
// actually arrives from. It is deliberately a separate walker: the forward
// Trace is left untouched.
//
// Semantic differences with System.Trace:
//
//   - The order of the interactions is the node list, not 0…n−1, and a surface
//     may be visited twice (the double pass of the ghost).
//   - At a refraction node the incident index is Mediums[Medium] and the exit
//     index the medium on the other side; Trace always uses Mediums[i] before
//     and Mediums[i+1] behind, i.e. the forward assignment.
//   - A GhostReflect node reflects even though Surface.Mirror is false; the
//     reflectance of that interface is not part of the trace (it enters the
//     relative energy estimate of the path, see System.GhostEnergy).
//   - Mirror and coordinate-break elements are rejected (System.GhostPathOf):
//     in Trace a mirror reflects nominally and a break has no interaction.
//   - The trace ends on the plane z = imageZ, like traceToPlane, and the
//     object-side geometry in front of surface 0 is not part of the replay:
//     the launch ray is taken as it is given.
//   - Path length is accumulated as in Trace, but no optical path length is
//     reported: ghost energy does not use it.
//
// Documented limitation: the replay is strictly sequential. The path order is
// taken as given, so occlusion by a surface that is not the next node (a
// reflected ray re-hitting the same concave surface, a baffle) is not detected
// — the same sequential assumption Trace makes. A zero-thickness element with
// curved surfaces is not geometrically realisable for a real ray (the two caps
// leave no glass thickness off the vertex), so such a path escapes; the
// paraxial chain has no such restriction, which is why the idealised thin lens
// is exercised through it.

// GhostStep records one interaction of a ghost replay.
type GhostStep struct {
	GhostNode
	Point    Vec3    // hit point (mm, global frame)
	Normal   Vec3    // unit surface normal at the hit point, +z oriented (global)
	DirIn    Vec3    // unit direction before the interaction
	DirOut   Vec3    // unit direction after the interaction
	T        float64 // path length of the segment ending at this node (mm)
	NBefore  float64 // refractive index in front of the surface, along the ray
	NAfter   float64 // refractive index beyond it (equal to NBefore for a reflection)
	IncAngle float64 // incidence angle (rad, in [0, π/2])
	Status   Status  // outcome of this interaction
}

// GhostHit is the outcome of replaying one ghost path for one ray.
type GhostHit struct {
	Status   Status      // StatusOK when the ray reached the image plane
	FailedAt int         // index in Path.Nodes of the failing interaction, −1 when OK
	Landing  Vec3        // image-plane intersection (valid when Status is StatusOK)
	Steps    []GhostStep // one entry per node reached, in path order
	Path     float64     // geometrical path length from the launch point (mm)
}

// OK reports whether the ray reached the image plane.
func (h GhostHit) OK() bool { return h.Status == StatusOK }

// String implements fmt.Stringer.
func (h GhostHit) String() string {
	return fmt.Sprintf("ghosthit[status=%s steps=%d landing=%s path=%.9gmm]",
		h.Status, len(h.Steps), h.Landing, h.Path)
}

// ReplayGhost propagates the ray r through the ghost path p and to the image
// plane z = imageZ, returning the landing point and one GhostStep per
// interaction. Statuses are the ones of the forward trace: StatusOK when the
// ray reached the image plane, otherwise the status of the interaction that
// stopped it (StatusVignetted, StatusTIR, StatusEscape, StatusDegenerate).
func (s *System) ReplayGhost(p GhostPath, r Ray, imageZ float64) GhostHit {
	out := GhostHit{FailedAt: -1}
	if !r.Valid() || !finite(imageZ) {
		out.Status = StatusDegenerate
		return out
	}
	if err := s.Validate(); err != nil {
		out.Status = StatusDegenerate
		return out
	}
	lambda := r.WavelengthNm
	o, d := r.Origin, r.Dir
	gap := 0 // the gap (index into Mediums) the ray travels in
	for k, nd := range p.Nodes {
		fail := func(st Status, step *GhostStep) GhostHit {
			step.Status = st
			out.Steps = append(out.Steps, *step)
			out.Status, out.FailedAt = st, k
			return out
		}
		if nd.Surface < 0 || nd.Surface+1 >= len(s.Mediums) {
			step := GhostStep{GhostNode: nd, DirIn: d, DirOut: d}
			return fail(StatusDegenerate, &step)
		}
		if nd.Medium != gap {
			// The node list does not follow the medium bookkeeping of the
			// path: a malformed argument, not a physical outcome.
			step := GhostStep{GhostNode: nd, DirIn: d, DirOut: d}
			return fail(StatusDegenerate, &step)
		}
		surf := &s.Surfaces[nd.Surface]
		step := GhostStep{GhostNode: nd, DirIn: d, DirOut: d}
		if surf.Break || surf.Profile == nil {
			return fail(StatusDegenerate, &step)
		}
		beyond := nd.Surface
		if !nd.Backward {
			beyond = nd.Surface + 1
		}
		nB := s.Mediums[gap].N(lambda)
		nA := s.Mediums[beyond].N(lambda)
		if nd.Action == GhostReflect {
			// A reflection does not change the medium: report the incident
			// index on both sides (the gap on the other side is not entered).
			nA = nB
		}
		step.NBefore, step.NAfter = nB, nA
		if !finitePositive(nB) || !finitePositive(s.Mediums[beyond].N(lambda)) {
			return fail(StatusDegenerate, &step)
		}
		e := s.elementFrame(nd.Surface)
		local := e.toLocal(o)
		point, t, hit := surf.Intersect(local, e.dirToLocal(d), tMinDefault)
		if !hit {
			return fail(StatusEscape, &step)
		}
		if t < 0 {
			t = 0
		}
		o = e.toGlobal(point)
		nrm := e.dirToGlobal(surf.Normal(point.X, point.Y))
		if nrm.Norm2() == 0 {
			return fail(StatusDegenerate, &step)
		}
		step.Point, step.Normal, step.T = o, nrm, t
		out.Path += t
		cos := math.Abs(d.Dot(nrm))
		if cos > 1 {
			cos = 1
		}
		step.IncAngle = math.Acos(cos)
		if surf.Vignetted(point.X, point.Y) {
			return fail(StatusVignetted, &step)
		}
		newDir := d
		switch nd.Action {
		case GhostReflect:
			newDir = Reflect(d, nrm)
		default:
			if nB != nA {
				nd2, ok := Refract(d, nrm, nB, nA)
				if !ok {
					return fail(StatusTIR, &step)
				}
				newDir = nd2
			}
		}
		step.Status = StatusOK
		step.DirOut = newDir
		out.Steps = append(out.Steps, step)
		d = newDir
		if nd.Action != GhostReflect {
			gap = beyond // crossing the surface changes the gap the ray is in
		}
	}
	if d.Z == 0 {
		out.Status = StatusEscape
		return out
	}
	t := (imageZ - o.Z) / d.Z
	if !finite(t) || t < 0 {
		out.Status = StatusEscape
		return out
	}
	out.Landing = o.Add(d.Mul(t))
	out.Path += t
	out.Status = StatusOK
	return out
}

// ---------------------------------------------------------------------------
// Paraxial ghost trace
// ---------------------------------------------------------------------------

// GhostParaxial is the paraxial launch state of a ghost analysis: the ray
// height and slope at the plane of surface 0 (z = System.VertexZ(0)), in the
// meridional (Y, U) and the sagittal (X, V) half spaces. The two halves are
// independent because every profile of the package is rotationally symmetric,
// so the same scalar engine is applied twice.
type GhostParaxial struct {
	Y, U float64 // height (mm) and slope dy/dz at the plane of surface 0
	X, V float64 // height (mm) and slope dx/dz at the plane of surface 0
}

// GhostParaxialOfRay extracts the paraxial launch state of the ray r at the
// plane z0: the height where the ray line pierces that plane and the slopes
// dy/dz, dx/dz. ok is false when the ray does not travel towards +z (dz <= 0)
// or is degenerate.
func GhostParaxialOfRay(r Ray, z0 float64) (GhostParaxial, bool) {
	if !r.Valid() || !finite(z0) || !(r.Dir.Z > 0) {
		return GhostParaxial{}, false
	}
	u := r.Dir.Y / r.Dir.Z
	v := r.Dir.X / r.Dir.Z
	dy := z0 - r.Origin.Z
	return GhostParaxial{
		Y: r.Origin.Y + dy*u,
		U: u,
		X: r.Origin.X + dy*v,
		V: v,
	}, true
}

// ghostParaxialChain propagates the same paraxial ray with explicit 2×2 ray
// transfer matrices in (y, ū) coordinates, ū = σ n u the reduced angle of the
// segment the ray is in (σ = +1 travelling towards +z, σ = −1 towards −z, n
// the index of the segment):
//
//	transfer    [[1, L/n], [0, 1]]      L = travelled distance, ≥ 0
//	refraction  [[1, 0], [−σΦ, 1]]      Φ = (n_out − n_in) c
//	reflection  [[1, 0], [2σn c, 1]]    after which σ flips
//
// The reflection matrix is the usual mirror folded into a refraction into an
// index −n; flipping σ is what keeps ū = σ n u consistent across it. The
// returned matrix is the product from the plane of surface 0 to the image
// plane. It is an independent transcription of the same path as
// the independent y–u recurrence oracle of the tests (TestGhostParaxialRepeatability),
// with which it agrees to round-off.
func ghostParaxialChain(els []paraxialElement, zs []float64, n0 float64, p GhostPair, y0, u0, imageZ float64) (float64, float64, Mat2) {
	n := len(els)
	y, ubar := y0, n0*u0
	sigma, idx := 1.0, n0
	cur := 0
	m := Identity()
	apply := func(mm Mat2) {
		y, ubar = mm.Apply(y, ubar)
		m = mm.Mul(m)
	}
	moveTo := func(k int) {
		L := sigma * (zs[k] - zs[cur])
		apply(Mat2{A: 1, B: L / idx, D: 1})
		cur = k
	}
	refract := func(c, nUp, nDown float64) {
		apply(Mat2{A: 1, C: -sigma * (nDown - nUp) * c, D: 1})
		idx = nDown
	}
	reflect := func(c float64) {
		apply(Mat2{A: 1, C: 2 * sigma * idx * c, D: 1})
		sigma = -sigma
	}
	for k := 0; k < p.I; k++ {
		if k > 0 {
			moveTo(k)
		}
		refract(els[k].c, els[k].nB, els[k].nA)
	}
	moveTo(p.I)
	reflect(els[p.I].c)
	for k := p.I - 1; k > p.J; k-- {
		moveTo(k)
		refract(els[k].c, els[k].nA, els[k].nB)
	}
	moveTo(p.J)
	reflect(els[p.J].c)
	for k := p.J + 1; k < n; k++ {
		moveTo(k)
		refract(els[k].c, els[k].nB, els[k].nA)
	}
	// Final transfer to the image plane.
	L := sigma * (imageZ - zs[cur])
	apply(Mat2{A: 1, B: L / idx, D: 1})
	uImg := sigma * ubar / idx
	return y, uImg, m
}

// GhostABCD returns the reversed ray transfer matrix of the ghost path p from
// the plane of surface 0 to the plane z = imageZ, in the (y, ū) coordinates of
// the package (ū = σ n u is the reduced angle of the segment the ray is in,
// σ = +1 towards +z and −1 towards −z, n the index of that segment). The chain
// is built from
//
//	transfer    [[1, L/n], [0, 1]]      L = travelled distance, ≥ 0
//	refraction  [[1, 0], [−σΦ, 1]]      Φ = (n_out − n_in) c
//	reflection  [[1, 0], [2σn c, 1]]    after which σ flips
//
// so the backwards pass enters with positive travelled distances and with the
// refraction power of the reversed direction, and every elementary matrix has
// determinant 1: det M = 1 is the symplectic invariant of the chain, checked
// by the tests. The matrix does not depend on the launch; the reference
// (y₀, ū₀ = n₀u₀) of GhostParaxial is applied to it by ParaxialGhostLanding.
func (s *System) GhostABCD(p GhostPair, imageZ, lambdaNm float64) (Mat2, error) {
	if err := s.Validate(); err != nil {
		return Mat2{}, err
	}
	if err := s.requireCoaxial("paraxial ghost analysis"); err != nil {
		return Mat2{}, err
	}
	if _, err := s.GhostPathOf(p); err != nil {
		return Mat2{}, err
	}
	if !finite(imageZ) {
		return Mat2{}, fmt.Errorf("geo: paraxial ghost: image plane z is not finite")
	}
	if !finitePositive(lambdaNm) {
		return Mat2{}, fmt.Errorf("geo: paraxial ghost: wavelength %g is not positive", lambdaNm)
	}
	els := s.paraxialElements(lambdaNm)
	zs := make([]float64, len(els))
	for i := range els {
		zs[i] = s.VertexZ(i)
	}
	_, _, m := ghostParaxialChain(els, zs, s.Mediums[0].N(lambdaNm), p, 0, 0, imageZ)
	return m, nil
}

// ParaxialGhostLanding traces the paraxial launch state through the ghost path
// p with the reversed ABCD chain (GhostABCD) and returns the landing point on
// the plane z = imageZ (mm, global frame). The meridional and the sagittal
// channel use the same matrix because every profile of the package is
// rotationally symmetric.
//
// The path must be coaxial (no decentre, no tilt — the first-order engine of
// the package is defined for coaxial systems only, see requireCoaxial) and may
// not contain mirrors or coordinate breaks at its nodes (System.GhostPathOf).
func (s *System) ParaxialGhostLanding(p GhostPair, launch GhostParaxial, imageZ, lambdaNm float64) (Vec3, error) {
	m, err := s.GhostABCD(p, imageZ, lambdaNm)
	if err != nil {
		return Vec3{}, err
	}
	n0 := s.Mediums[0].N(lambdaNm)
	yImg, _ := m.Apply(launch.Y, n0*launch.U)
	xImg, _ := m.Apply(launch.X, n0*launch.V)
	return V(xImg, yImg, imageZ), nil
}
