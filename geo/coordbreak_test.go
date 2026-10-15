package geo

import (
	"math"
	"math/rand"
	"testing"
)

// These tests pin the coordinate-break convention analytically and check the
// algebra the rest of the package relies on: the return break, the Cardan
// decomposition and the frame chain.

// randomBreak draws a coordinate break with moderate parameters from rng.
func randomBreak(rng *rand.Rand) CoordBreak {
	return CoordBreak{
		Dx:    rng.Float64()*4 - 2,
		Dy:    rng.Float64()*4 - 2,
		Dz:    rng.Float64()*4 - 2,
		TiltX: rng.Float64()*1.2 - 0.6,
		TiltY: rng.Float64()*1.2 - 0.6,
		TiltZ: rng.Float64()*2 - 1,
	}
}

// TestCoordBreakConventionAnalytic checks the explicit formulas of the
// convention p_old = d + Rx(TiltX)·Ry(TiltY)·Rz(TiltZ)·p_new at angles where
// the answer is known in closed form.
func TestCoordBreakConventionAnalytic(t *testing.T) {
	// Pure decentre: the origin moves, the axes do not.
	p := CoordBreak{Dx: 1.5, Dy: -2}.Apply(V(0.25, 0.5, 3))
	checkClose(t, "decentre x", p.X, 1.75, 1e-15)
	checkClose(t, "decentre y", p.Y, -1.5, 1e-15)
	checkClose(t, "decentre z", p.Z, 3, 1e-15)

	// Tilt about x by +90°: the local y axis goes to the global +z axis,
	// Rx(π/2)·(0,1,0) = (0, cos, sin) = (0,0,1).
	q := CoordBreak{TiltX: math.Pi / 2}.Apply(V(0, 1, 0))
	checkClose(t, "Rx(90) y->z x", q.X, 0, 1e-15)
	checkClose(t, "Rx(90) y->z y", q.Y, 0, 1e-15)
	checkClose(t, "Rx(90) y->z z", q.Z, 1, 1e-15)

	// Tilt about y by +90°: Ry(π/2)·(1,0,0) = (cos, 0, -sin) = (0,0,-1).
	r := CoordBreak{TiltY: math.Pi / 2}.Apply(V(1, 0, 0))
	checkClose(t, "Ry(90) x->-z x", r.X, 0, 1e-15)
	checkClose(t, "Ry(90) x->-z z", r.Z, -1, 1e-15)

	// Roll about z by +90°: Rz(π/2)·(1,0,0) = (0,1,0).
	s := CoordBreak{TiltZ: math.Pi / 2}.Apply(V(1, 0, 0))
	checkClose(t, "Rz(90) x->y", s.Y, 1, 1e-15)
	checkClose(t, "Rz(90) x->x", s.X, 0, 1e-15)

	// Decentre first, then tilt about the decentred origin: the rotated point
	// is measured from the new origin (1, 0, 0).
	u := CoordBreak{Dx: 1, TiltX: math.Pi / 2}.Apply(V(0, 1, 0))
	checkClose(t, "decentre then tilt x", u.X, 1, 1e-15)
	checkClose(t, "decentre then tilt z", u.Z, 1, 1e-15)
}

// TestCoordBreakInverseIsReturn checks that the return break undoes the break
// in both orders, for random breaks: B∘B⁻¹ = B⁻¹∘B = identity to 1e-13 mm.
func TestCoordBreakInverseIsReturn(t *testing.T) {
	rng := rand.New(rand.NewSource(20261015))
	worst := 0.0
	for k := 0; k < 200; k++ {
		cb := randomBreak(rng)
		inv := cb.Inverse()
		p := V(rng.Float64()*20-10, rng.Float64()*20-10, rng.Float64()*20-10)
		a := cb.Apply(inv.Apply(p)).Sub(p).Norm()
		b := inv.Apply(cb.Apply(p)).Sub(p).Norm()
		worst = math.Max(worst, math.Max(a, b))
	}
	if worst > 1e-13 {
		t.Errorf("B∘B⁻¹ / B⁻¹∘B not identity: worst deviation %.3g mm (want <= 1e-13)", worst)
	}
	t.Logf("worst |B∘B⁻¹ − id| over 200 random breaks = %.3g mm", worst)
}

// TestCoordBreakInverseOfInverse checks Inverse(Inverse(B)) = B as a map.
func TestCoordBreakInverseOfInverse(t *testing.T) {
	rng := rand.New(rand.NewSource(7))
	worst := 0.0
	for k := 0; k < 100; k++ {
		cb := randomBreak(rng)
		back := cb.Inverse().Inverse()
		p := V(rng.Float64(), rng.Float64(), rng.Float64())
		worst = math.Max(worst, cb.Apply(p).Sub(back.Apply(p)).Norm())
	}
	if worst > 1e-12 {
		t.Errorf("Inverse(Inverse(B)) differs from B: %.3g mm", worst)
	}
}

// maxRotationError returns the largest entry error of the reconstruction
// Rx(a)·Ry(b)·Rz(c) of the Cardan angles of m.
func maxRotationError(m mat3) float64 {
	a, b, c := eulerXYZ(m)
	r := rotX(a).mul(rotY(b)).mul(rotZ(c))
	worst := 0.0
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			worst = math.Max(worst, math.Abs(r[i][j]-m[i][j]))
		}
	}
	return worst
}

// TestEulerXYZDecomposition checks that eulerXYZ reproduces every rotation
// Rx(a)·Ry(b)·Rz(c) for |b| <= 1.4 rad (1e-12, the rotation entries are
// well-conditioned there) and the gimbal case b = ±π/2 (1e-12 as well: the
// split is arbitrary, the matrix is not). Near the gimbal lock (|cos b| = 1e-7)
// the entries of the matrix carry only ~1e-16/1e-7 relative information, so the
// reconstruction is checked to 1e-8 there.
func TestEulerXYZDecomposition(t *testing.T) {
	rng := rand.New(rand.NewSource(11))
	worst := 0.0
	for k := 0; k < 300; k++ {
		m := rotX(rng.Float64()*6 - 3).mul(rotY(rng.Float64()*2.8 - 1.4)).mul(rotZ(rng.Float64()*6 - 3))
		worst = math.Max(worst, maxRotationError(m))
	}
	for _, m := range []mat3{
		rotX(0.3).mul(rotY(math.Pi / 2)).mul(rotZ(1.1)),
		rotX(0.8).mul(rotY(-math.Pi / 2)).mul(rotZ(-0.2)),
		rotX(0).mul(rotY(0)).mul(rotZ(0)),
	} {
		worst = math.Max(worst, maxRotationError(m))
	}
	if worst > 1e-12 {
		t.Errorf("Cardan decomposition does not reproduce the rotation: max entry error %.3g", worst)
	}
	near := rotX(0.5).mul(rotY(math.Pi/2 - 1e-7)).mul(rotZ(-0.4))
	if e := maxRotationError(near); e > 1e-8 {
		t.Errorf("near-gimbal reconstruction error %.3g, want <= 1e-8", e)
	}
	t.Logf("Cardan reconstruction: worst entry error %.3g (|b|<=1.4 and gimbal), near-gimbal %.3g",
		worst, maxRotationError(near))
}

// TestFrameRoundTrip checks toLocal∘toGlobal and dir transforms of a tilted
// frame, and that a frame built from a break maps the break's origin correctly.
func TestFrameRoundTrip(t *testing.T) {
	rng := rand.New(rand.NewSource(3))
	f := frame{o: V(0.4, -0.3, 12), r: ident3()}.apply(randomBreak(rng))
	for k := 0; k < 50; k++ {
		p := V(rng.Float64(), rng.Float64(), rng.Float64())
		back := f.toLocal(f.toGlobal(p))
		if d := back.Sub(p).Norm(); d > 1e-13 {
			t.Fatalf("toLocal(toGlobal(p)) != p: %.3g", d)
		}
		d := V(rng.Float64()-0.5, rng.Float64()-0.5, 1).Normalized()
		if e := f.dirToLocal(f.dirToGlobal(d)).Sub(d).Norm(); e > 1e-13 {
			t.Fatalf("dirToLocal(dirToGlobal(d)) != d: %.3g", e)
		}
	}
	// The origin of the frame of a break is the decentre seen from the parent.
	cb := CoordBreak{Dx: 0.7, Dy: -0.2, TiltX: 0.4}
	child := frame{o: V(0, 0, 5), r: ident3()}.apply(cb)
	if d := child.o.Sub(V(0.7, -0.2, 5)).Norm(); d > 1e-15 {
		t.Errorf("break origin %v, want (0.7, -0.2, 5)", child.o)
	}
}

// TestChainTrivialForCoaxialSystem checks the chain machinery on a coaxial
// system: VertexZ and the element frames are the plain thickness sums.
func TestChainTrivialForCoaxialSystem(t *testing.T) {
	sys := mustSystem(t, Air{},
		[]Surface{NewSurface(Sphere{R: 30}, 5), NewSurface(Plane{}, 5), NewSurface(Sphere{R: -20}, 5)},
		[]float64{3, 2.5},
		[]Medium{Air{}, ConstantIndex(1.5), ConstantIndex(1.6), Air{}})
	sys.OriginZ = -1.25
	checkClose(t, "vertex 0", sys.VertexZ(0), -1.25, 0)
	checkClose(t, "vertex 1", sys.VertexZ(1), 1.75, 0)
	checkClose(t, "vertex 2", sys.VertexZ(2), 4.25, 0)
	if sys.hasPose() {
		t.Fatal("a coaxial system must not report poses")
	}
	if e := sys.elementFrame(2); e.r != ident3() || e.o != V(0, 0, 4.25) {
		t.Errorf("element 2 frame %+v, want identity at z=4.25", e)
	}
}
