package geo

import (
	"fmt"
	"math"
)

// ---------------------------------------------------------------------------
// Coordinate breaks (decentre and tilt)
// ---------------------------------------------------------------------------

// CoordBreak is a coordinate break in the sense of Zemax's COORDBRK surface: a
// decentre of a frame followed by tilts about its axes. It is the only way a
// surface or a sequence of surfaces is decentred or tilted in package geo.
//
// Convention (the same everywhere in the package; see also the package
// documentation). Let Q be the rotation
//
//	Q = Rx(TiltX) · Ry(TiltY) · Rz(TiltZ) ,
//
// with the right-handed active rotations
//
//	Rx(a) = [[1, 0, 0], [0, c, -s], [0, s, c]]
//	Ry(a) = [[c, 0, s], [0, 1, 0], [-s, 0, c]]
//	Rz(a) = [[c, -s, 0], [s, c, 0], [0, 0, 1]]      (c = cos a, s = sin a).
//
// The break maps local coordinates p_new of the new frame into the coordinates
// p_old of the enclosing frame:
//
//	p_old = d + Q · p_new ,        d = (Dx, Dy, Dz) .
//
// Read as a sequence of operations: first the origin is moved by the
// displacement d measured in the enclosing frame (the decentre (Dx, Dy) and an
// axial shift Dz), then the axes are tilted about the displaced origin: first
// about the x axis by TiltX, then about the tilted y axis by TiltY, then about
// the tilted z axis by TiltZ (intrinsic order). The tilts are the angles of the
// Q matrix above; a positive angle is a right-handed (counter-clockwise, seen
// from the positive axis) rotation. Dz is zero for an ordinary decentre; the
// return break of a tilted break has a non-zero Dz in general (see Inverse).
//
// A surface placed by a break has its vertex at the new origin and its local
// optical axis along the new z axis, Q·(0,0,1) in the enclosing frame. For a
// surface carrying no break the local axis is the global +z.
//
// Return semantics. A break followed by its return break (Inverse) at zero axial
// spacing is the identity transform: the return has displacement d' = -Qᵀ d and
// rotation Q' = Qᵀ, which is again written as Rx(a)·Ry(b)·Rz(c) (the angles a, b,
// c are obtained by a Cardan decomposition of Qᵀ). Placed after a break with
// nothing in between, the return restores the optical axis of the enclosing
// frame. A non-zero axial spacing between a break and its return is a
// translation along the broken axis and must be included in the Dz of the
// return; the identity is tested for zero spacing only.
//
// The angle fields are in radians. |TiltY| = π/2 is a gimbal lock of the
// inverse decomposition; the forward transform is well-defined there.
type CoordBreak struct {
	Dx, Dy       float64 // decentre (mm) in the enclosing frame, applied first
	Dz           float64 // axial displacement (mm) in the enclosing frame, applied with the decentre
	TiltX, TiltY float64 // tilts about x, then about the tilted y axis (rad)
	TiltZ        float64 // roll about the tilted z axis (rad), applied last
}

// IsIdentity reports whether the break leaves the frame unchanged (all
// parameters exactly zero).
func (cb CoordBreak) IsIdentity() bool {
	return cb.Dx == 0 && cb.Dy == 0 && cb.Dz == 0 && cb.TiltX == 0 && cb.TiltY == 0 && cb.TiltZ == 0
}

// decentredOrTilted reports whether the break moves the optical axis, i.e.
// whether it has a decentre, an axial shift or a tilt about x or y. A roll about
// z alone leaves the axis (and the first-order optics of rotationally symmetric
// surfaces) unchanged.
func (cb CoordBreak) decentredOrTilted() bool {
	return cb.Dx != 0 || cb.Dy != 0 || cb.Dz != 0 || cb.TiltX != 0 || cb.TiltY != 0
}

// finite reports whether every parameter is a finite number.
func (cb CoordBreak) finite() bool {
	for _, v := range [...]float64{cb.Dx, cb.Dy, cb.Dz, cb.TiltX, cb.TiltY, cb.TiltZ} {
		if math.IsNaN(v) || math.IsInf(v, 0) {
			return false
		}
	}
	return true
}

// Apply maps a point in the local coordinates of the break into the enclosing
// frame: p_old = d + Q·p_new.
func (cb CoordBreak) Apply(p Vec3) Vec3 {
	return cb.displacement().Add(cb.rotation().mulVec(p))
}

// displacement returns the decentre vector d = (Dx, Dy, Dz).
func (cb CoordBreak) displacement() Vec3 { return Vec3{cb.Dx, cb.Dy, cb.Dz} }

// Inverse returns the return break of cb: the break whose composition with cb
// (cb first, then the return, zero spacing) is the identity. With Q the rotation
// of cb and d its displacement, the return has displacement d' = -Qᵀ d (in
// general with a non-zero Dz, even when cb has none) and the rotation Qᵀ written
// as Rx(a)·Ry(b)·Rz(c). Inverse(Inverse(cb)) equals cb up to round-off.
func (cb CoordBreak) Inverse() CoordBreak {
	qt := cb.rotation().transpose()
	dn := qt.mulVec(cb.displacement()).Neg()
	a, b, c := eulerXYZ(qt)
	return CoordBreak{Dx: dn.X, Dy: dn.Y, Dz: dn.Z, TiltX: a, TiltY: b, TiltZ: c}
}

// rotation returns Q = Rx(TiltX)·Ry(TiltY)·Rz(TiltZ).
func (cb CoordBreak) rotation() mat3 {
	return rotX(cb.TiltX).mul(rotY(cb.TiltY)).mul(rotZ(cb.TiltZ))
}

// ---------------------------------------------------------------------------
// Small linear algebra used by the frame machinery (unexported)
// ---------------------------------------------------------------------------

// mat3 is a 3×3 matrix, row-major: m[row][col].
type mat3 [3][3]float64

// ident3 returns the identity matrix.
func ident3() mat3 { return mat3{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}} }

// mul returns m·n.
func (m mat3) mul(n mat3) mat3 {
	var out mat3
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			out[i][j] = m[i][0]*n[0][j] + m[i][1]*n[1][j] + m[i][2]*n[2][j]
		}
	}
	return out
}

// mulVec returns m·v.
func (m mat3) mulVec(v Vec3) Vec3 {
	return Vec3{
		m[0][0]*v.X + m[0][1]*v.Y + m[0][2]*v.Z,
		m[1][0]*v.X + m[1][1]*v.Y + m[1][2]*v.Z,
		m[2][0]*v.X + m[2][1]*v.Y + m[2][2]*v.Z,
	}
}

// transpose returns mᵀ.
func (m mat3) transpose() mat3 {
	var out mat3
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			out[i][j] = m[j][i]
		}
	}
	return out
}

// col returns column j of m (the image of the j-th unit vector).
func (m mat3) col(j int) Vec3 { return Vec3{m[0][j], m[1][j], m[2][j]} }

// rotX, rotY, rotZ are the right-handed active rotations about the axes.
func rotX(a float64) mat3 {
	c, s := math.Cos(a), math.Sin(a)
	return mat3{{1, 0, 0}, {0, c, -s}, {0, s, c}}
}

func rotY(a float64) mat3 {
	c, s := math.Cos(a), math.Sin(a)
	return mat3{{c, 0, s}, {0, 1, 0}, {-s, 0, c}}
}

func rotZ(a float64) mat3 {
	c, s := math.Cos(a), math.Sin(a)
	return mat3{{c, -s, 0}, {s, c, 0}, {0, 0, 1}}
}

// eulerXYZ returns the angles (a, b, c) with m = Rx(a)·Ry(b)·Rz(c), for a
// rotation matrix m. The decomposition is the Cardan (XYZ) one:
//
//	b = asin(m[0][2]),  a = atan2(-m[1][2], m[2][2]),  c = atan2(-m[0][1], m[0][0]).
//
// At the gimbal lock cos b = 0 the split between a and c is arbitrary; the
// convention c = 0 is used and a = atan2(m[2][1], m[1][1]).
func eulerXYZ(m mat3) (a, b, c float64) {
	b = math.Asin(math.Max(-1, math.Min(1, m[0][2])))
	if math.Abs(math.Cos(b)) > 1e-9 {
		a = math.Atan2(-m[1][2], m[2][2])
		c = math.Atan2(-m[0][1], m[0][0])
		return a, b, c
	}
	return math.Atan2(m[2][1], m[1][1]), b, 0
}

// frame is a rigid transform from a local coordinate system into the
// enclosing (global) frame: P = o + r·p, where the columns of r are the local
// axes expressed in the global frame.
type frame struct {
	o     Vec3
	r     mat3
	plain bool // r is exactly the identity: the untilted fast path (exact arithmetic)
}

// toGlobal maps a local point to the global frame.
func (f frame) toGlobal(p Vec3) Vec3 {
	if f.plain {
		return f.o.Add(p) // exact: the arithmetic of the untilted layout
	}
	return f.o.Add(f.r.mulVec(p))
}

// toLocal maps a global point to the local frame: p = rᵀ (P − o).
func (f frame) toLocal(P Vec3) Vec3 {
	if f.plain {
		return P.Sub(f.o) // exact: the arithmetic of the untilted layout
	}
	return f.r.transpose().mulVec(P.Sub(f.o))
}

// dirToGlobal maps a direction (or a unit normal) to the global frame.
func (f frame) dirToGlobal(d Vec3) Vec3 {
	if f.plain {
		return d
	}
	return f.r.mulVec(d)
}

// dirToLocal maps a global direction into the local frame.
func (f frame) dirToLocal(d Vec3) Vec3 {
	if f.plain {
		return d
	}
	return f.r.transpose().mulVec(d)
}

// apply returns the frame obtained by the coordinate break cb applied to f:
// origin o + r·d, axes r·Q (d = (Dx, Dy, Dz)).
func (f frame) apply(cb CoordBreak) frame {
	r := f.r.mul(cb.rotation())
	return frame{
		o:     f.o.Add(f.r.mulVec(cb.displacement())),
		r:     r,
		plain: f.plain && r == ident3(),
	}
}

// ---------------------------------------------------------------------------
// Chain of frames of a system
// ---------------------------------------------------------------------------

// chainStart returns the frame of the system origin: the vertex of element 0
// sits at z = OriginZ and the axes are the global axes.
func (s *System) chainStart() frame {
	return frame{o: Vec3{0, 0, s.OriginZ}, r: ident3(), plain: true}
}

// place returns the frame e in which the local coordinates of element i are
// defined, given the chain frame c that precedes it, and the reference frame
// ref from which the gap Thickness[i] to the next element is measured.
//
//   - A coordinate break element (Surface.Break) has e = c·Pose and changes the
//     chain: its own frame is the reference of the following gap.
//   - A surface with a pose (Surface.Pose) has e = c·Pose, but the chain is not
//     changed: the gap is measured from c (auto-return to the optical axis).
//   - A plain surface has e = c.
func (s *System) place(i int, c frame) (e, ref frame) {
	sf := &s.Surfaces[i]
	if sf.Break {
		if sf.Pose == nil {
			return c, c
		}
		e = c.apply(*sf.Pose)
		return e, e
	}
	if sf.Pose == nil {
		return c, c
	}
	return c.apply(*sf.Pose), c
}

// advance returns the chain frame that precedes element i+1: the origin is
// moved by Thickness[i] along the axis of ref. After the last element the chain
// is returned unchanged.
func (s *System) advance(ref frame, i int) frame {
	if i < 0 || i >= len(s.Thickness) {
		return ref
	}
	if ref.plain {
		return frame{o: Vec3{ref.o.X, ref.o.Y, ref.o.Z + s.Thickness[i]}, r: ref.r, plain: true} // exact: untilted chain
	}
	return frame{o: ref.o.Add(ref.r.col(2).Mul(s.Thickness[i])), r: ref.r}
}

// elementFrame returns the frame of element i (see place).
func (s *System) elementFrame(i int) frame {
	c := s.chainStart()
	for k := 0; k < i; k++ {
		_, ref := s.place(k, c)
		c = s.advance(ref, k)
	}
	e, _ := s.place(i, c)
	return e
}

// hasPose reports whether the system contains a coordinate-break element or a
// surface with a non-identity pose (decentre, tilt or roll).
func (s *System) hasPose() bool {
	for _, sf := range s.Surfaces {
		if sf.Break || (sf.Pose != nil && !sf.Pose.IsIdentity()) {
			return true
		}
	}
	return false
}

// hasAsymmetricFreeform reports whether the system contains a freeform surface
// that is not rotationally symmetric (see FreeformProfile). Such a surface has
// no rotationally symmetric surrogate, so the operations that rely on one
// (the first- and third-order engines, Reversed) refuse it.
func (s *System) hasAsymmetricFreeform() bool {
	for _, sf := range s.Surfaces {
		if sf.Break || sf.Profile == nil {
			continue
		}
		if fp, ok := sf.Profile.(FreeformProfile); ok && !fp.RotationallySymmetric() {
			return true
		}
	}
	return false
}

// requireCoaxial returns an error naming the first element whose pose decentres
// or tilts the optical axis, or whose profile is an asymmetric freeform surface.
// The first-order (paraxial) and third-order (Seidel) engines are defined for
// rotationally symmetric coaxial systems only and must not silently produce
// numbers for a decentred, tilted or freeform layout. A rotationally symmetric
// XYPolynomial is not freeform in this sense (it is an ordinary rotationally
// symmetric profile) and passes.
func (s *System) requireCoaxial(what string) error {
	for i, sf := range s.Surfaces {
		if sf.Pose != nil && sf.Pose.decentredOrTilted() {
			kind := "surface"
			if sf.Break {
				kind = "coordinate break"
			}
			return errNotCoaxial(what, i, kind, sf.Label(), *sf.Pose)
		}
		if sf.Break || sf.Profile == nil {
			continue
		}
		if fp, ok := sf.Profile.(FreeformProfile); ok && !fp.RotationallySymmetric() {
			return errFreeform(what, i, sf.Label())
		}
	}
	return nil
}

// errFreeform builds the error of the coaxial-only engines for a freeform
// surface that is not rotationally symmetric.
func errFreeform(what string, i int, label string) error {
	return fmt.Errorf("geo: %s: surface %d (%s) is a freeform surface that is not rotationally symmetric; "+
		"first- and third-order optics are defined for rotationally symmetric surfaces only: "+
		"use real ray tracing (Trace, SpotDiagram, RayFan, OPDGrid with an explicit Pupil)",
		what, i, label)
}

// errNotCoaxial builds the error of the coaxial-only engines.
func errNotCoaxial(what string, i int, kind, label string, p CoordBreak) error {
	return fmt.Errorf("geo: %s: %s %d (%s) is decentred or tilted (dx=%g dy=%g tiltX=%g tiltY=%g rad); "+
		"first- and third-order optics are defined for coaxial systems only: "+
		"use real ray tracing (Trace, SpotDiagram, RayFan, OPDGrid with an explicit Pupil)",
		what, kind, i, label, p.Dx, p.Dy, p.TiltX, p.TiltY)
}
