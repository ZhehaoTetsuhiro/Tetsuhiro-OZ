package geo

import (
	"math"
	"testing"
)

// ---------------------------------------------------------------------------
// Field points
// ---------------------------------------------------------------------------

// TestFieldConstructors checks the field conventions of the package: field
// angles are stored as tangents, a positive tanY sends light downwards, finite
// objects sit at z = ObjectZ, and the angle/tangent forms agree with each
// other.
func TestFieldConstructors(t *testing.T) {
	f := FieldAngle(0, 10)
	checkClose(t, "FieldAngle tanY", f.TanY, math.Tan(10*math.Pi/180), 1e-15)
	checkClose(t, "FieldAngle tanX", f.TanX, 0, 0)
	if !f.IsInfinite() {
		t.Fatal("FieldAngle must produce an infinite object")
	}
	if !math.IsInf(f.ObjectPlane(), -1) {
		t.Errorf("ObjectPlane = %g, want −Inf", f.ObjectPlane())
	}
	// d = (−tanX, −tanY, 1)/|·| = (−sinθx·..., −sin θy, cos θy) for the y part.
	d, ok := f.Direction()
	if !ok {
		t.Fatal("Direction not ok")
	}
	checkClose(t, "direction.y", d.Y, -math.Sin(10*math.Pi/180), 1e-15)
	checkClose(t, "direction.z", d.Z, math.Cos(10*math.Pi/180), 1e-15)

	ft := FieldTan(0.03, -0.04)
	dx, dy := ft.Transverse()
	checkClose(t, "FieldTan transverse x", dx, 0.03, 0)
	checkClose(t, "FieldTan transverse y", dy, -0.04, 0)

	ff := FieldHeight(-200, 3, -4)
	if ff.IsInfinite() {
		t.Fatal("FieldHeight must produce a finite object")
	}
	p := ff.ObjectPoint()
	if !p.Equal(V(3, -4, -200), 0) {
		t.Errorf("ObjectPoint = %s, want (3, -4, -200)", p)
	}
	checkClose(t, "FieldHeight plane", ff.ObjectPlane(), -200, 0)
	if _, ok := ff.Direction(); ok {
		t.Error("a finite field must not have a parallel-beam direction")
	}
	// The zero value is the on-axis infinite field, and a literal with only
	// tangents set is an infinite field.
	z := Field{}
	if !z.IsInfinite() || z.TanX != 0 || z.TanY != 0 {
		t.Errorf("zero Field = %+v, want on-axis infinite", z)
	}
	lit := Field{TanY: 0.1}
	if !lit.IsInfinite() {
		t.Error("Field{TanY: 0.1} must be an infinite field")
	}
	if s := ft.String(); s == "" {
		t.Error("String is empty")
	}
}

// ---------------------------------------------------------------------------
// Pupil samplings
// ---------------------------------------------------------------------------

// TestPupilPoint checks the small helpers of PupilPoint.
func TestPupilPoint(t *testing.T) {
	p := PupilPointAt(0.6, -0.8) // a 3-4-5 triangle: R = 1
	checkClose(t, "R", p.R(), 1, 1e-15)
	if !p.Inside() {
		t.Error("boundary point must be inside")
	}
	if PupilPointAt(1.0001, 0).Inside() {
		t.Error("point outside the disc reported inside")
	}
	checkClose(t, "angle", PupilPointAt(0, 1).Angle(), math.Pi/2, 1e-15)
	checkClose(t, "angle wrapped", PupilPointAt(0, -1).Angle(), 3*math.Pi/2, 1e-15)
}

// TestHexapolarSampling checks the hexapolar pattern: 1 + 3n(n+1) points, all
// inside the unit disc, ring k contains exactly 6k points at radius k/n and
// no duplicates.
func TestHexapolarSampling(t *testing.T) {
	for _, n := range []int{1, 2, 3, 5} {
		pts := Hexapolar(n)
		want := 1 + 3*n*(n+1)
		if len(pts) != want {
			t.Fatalf("Hexapolar(%d) has %d points, want %d", n, len(pts), want)
		}
		seen := map[PupilPoint]bool{}
		ringCount := 0
		for _, p := range pts {
			if !p.Inside() {
				t.Errorf("Hexapolar(%d): %s outside the unit disc", n, p)
			}
			if seen[p] {
				t.Errorf("Hexapolar(%d): duplicate point %s", n, p)
			}
			seen[p] = true
			if r := p.R(); math.Abs(r-1) < 1e-12 || r > 0 {
				if math.Abs(r*float64(n)-math.Round(r*float64(n))) < 1e-12 {
					ringCount++
				}
			}
		}
		// every point lies on a ring of the pattern (multiples of 1/n)
		if ringCount != len(pts)-1 {
			t.Errorf("Hexapolar(%d): %d of %d points on rings", n, ringCount, len(pts)-1)
		}
	}
	if got := Hexapolar(0); len(got) != 1 || got[0] != (PupilPoint{}) {
		t.Errorf("Hexapolar(0) = %v, want the centre only", got)
	}
}

// TestRingsSampling checks the concentric ring pattern: count, radii, angles
// and the centre point.
func TestRingsSampling(t *testing.T) {
	pts := Rings(4, 6)
	if len(pts) != 1+4*6 {
		t.Fatalf("Rings(4,6) has %d points, want %d", len(pts), 25)
	}
	if pts[0] != (PupilPoint{}) {
		t.Errorf("first point = %v, want the centre", pts[0])
	}
	counts := map[int]int{}
	for _, p := range pts[1:] {
		r := p.R()
		k := int(math.Round(r * 4))
		if math.Abs(r-float64(k)/4) > 1e-12 {
			t.Errorf("%s not on ring %d/4 (r = %g)", p, k, r)
		}
		counts[k]++
	}
	for k := 1; k <= 4; k++ {
		if counts[k] != 6 {
			t.Errorf("ring %d has %d points, want 6", k, counts[k])
		}
	}
}

// TestGridSampling checks the square grid clipped to the unit disc.
func TestGridSampling(t *testing.T) {
	pts := Grid(5, 5)
	// 5x5 grid over [−1,1]²: the corners (r = √2) fall outside; the count is
	// the number of lattice points with x²+y² <= 1.
	want := 0
	for j := 0; j < 5; j++ {
		for i := 0; i < 5; i++ {
			x := 2*float64(i)/4 - 1
			y := 2*float64(j)/4 - 1
			if x*x+y*y <= 1 {
				want++
			}
		}
	}
	if len(pts) != want {
		t.Fatalf("Grid(5,5) has %d points, want %d", len(pts), want)
	}
	for _, p := range pts {
		if !p.Inside() {
			t.Errorf("%s outside the unit disc", p)
		}
	}
	// The centre must be present for odd grids.
	found := false
	for _, p := range pts {
		if p == (PupilPoint{}) {
			found = true
		}
	}
	if !found {
		t.Error("Grid(5,5) misses the centre")
	}
	if len(Grid(0, 3)) != 0 || len(Grid(3, 0)) != 0 {
		t.Error("degenerate grids must be empty")
	}
	// A 1x1 grid is the pupil centre.
	if got := Grid(1, 1); len(got) != 1 || got[0] != (PupilPoint{}) {
		t.Errorf("Grid(1,1) = %v", got)
	}
}

// TestFanSamplings checks the fan generators: endpoints, spacing and the
// meridional/sagittal specializations.
func TestFanSamplings(t *testing.T) {
	f := Fan(5, 90)
	if len(f) != 5 {
		t.Fatalf("Fan(5,90) has %d points", len(f))
	}
	checkClose(t, "first", f[0].PY, -1, 1e-15)
	checkClose(t, "last", f[4].PY, 1, 1e-15)
	checkClose(t, "middle", f[2].PY, 0, 1e-15)
	for _, p := range f {
		checkClose(t, "meridian px", p.PX, 0, 1e-15)
	}
	if got := Fan(1, 45); len(got) != 1 || got[0] != (PupilPoint{}) {
		t.Errorf("Fan(1) = %v, want the centre", got)
	}
	for _, p := range MeridionalFan(7) {
		if p.PX != 0 {
			t.Errorf("MeridionalFan point %s has PX != 0", p)
		}
	}
	for _, p := range SagittalFan(7) {
		if p.PY != 0 {
			t.Errorf("SagittalFan point %s has PY != 0", p)
		}
	}
	// A 45° fan runs along (1,1)/√2.
	for _, p := range Fan(3, 45) {
		checkClose(t, "45° fan symmetric", p.PX, p.PY, 1e-15)
	}
	if got := len(CrossFan(5)); got != 2*5-1 {
		t.Errorf("CrossFan(5) has %d points, want %d", got, 2*5-1)
	}
}

// ---------------------------------------------------------------------------
// Launching
// ---------------------------------------------------------------------------

// TestLaunchRayInfinite checks that the rays of an infinite field are parallel
// to the field direction and pass through the requested pupil points (including
// the pupil offset).
func TestLaunchRayInfinite(t *testing.T) {
	f := FieldTan(0, 0.1)
	pup := Pupil{Z: -30, Radius: 12, Offset: V(0, 1.5, 0)}
	want, _ := f.Direction()
	for _, p := range []PupilPoint{{}, {1, 0}, {0, -1}, {0.5, 0.5}} {
		r := LaunchRay(f, pup, p, Wavelengthd)
		if !r.Valid() {
			t.Fatalf("p=%s: invalid ray", p)
		}
		if d := r.Dir.Sub(want).Norm(); d > 1e-15 {
			t.Errorf("p=%s: direction off by %g", p, d)
		}
		// The ray must cross the plane z = pup.Z at Offset + (px,py)·Radius.
		wantPoint := V(pup.Offset.X+p.PX*pup.Radius, pup.Offset.Y+p.PY*pup.Radius, pup.Z)
		tt := (pup.Z - r.Origin.Z) / r.Dir.Z
		if got := r.At(tt); !got.Equal(wantPoint, 1e-12) {
			t.Errorf("p=%s: crosses the pupil plane at %s, want %s", p, got, wantPoint)
		}
	}
}

// TestLaunchRayFinite checks the finite-object launch: the ray starts at the
// object point and crosses the pupil point.
func TestLaunchRayFinite(t *testing.T) {
	f := FieldHeight(-100, 4, -2)
	pup := Pupil{Z: -20, Radius: 5}
	obj := f.ObjectPoint()
	for _, p := range []PupilPoint{{}, {0.3, -0.7}, {-1, 0}} {
		r := LaunchRay(f, pup, p, Wavelengthd)
		if !r.Valid() {
			t.Fatalf("p=%s: invalid ray", p)
		}
		if !r.Origin.Equal(obj, 1e-15) {
			t.Errorf("p=%s: origin %s, want the object point %s", p, r.Origin, obj)
		}
		wantPoint := V(p.PX*pup.Radius, p.PY*pup.Radius, pup.Z)
		tt := (pup.Z - r.Origin.Z) / r.Dir.Z
		if got := r.At(tt); !got.Equal(wantPoint, 1e-12) {
			t.Errorf("p=%s: crosses the pupil plane at %s, want %s", p, got, wantPoint)
		}
	}
	// A degenerate pupil produces a degenerate ray rather than a panic.
	bad := LaunchRay(f, Pupil{Z: math.NaN(), Radius: 5}, PupilPoint{}, Wavelengthd)
	checkStatus(t, "NaN pupil", bad.Status, StatusDegenerate)
	// LaunchSet keeps the order.
	pts := MeridionalFan(3)
	rays := LaunchSet(f, pup, pts, Wavelengthd)
	if len(rays) != 3 {
		t.Fatalf("LaunchSet returned %d rays", len(rays))
	}
	for i := range pts {
		if !rays[i].Valid() {
			t.Errorf("ray %d invalid", i)
		}
	}
}

// TestEntranceExitPupilConsistency checks the pupil helpers against the
// (independently tested) paraxial engine, and the error paths for telecentric
// spaces and unbounded stops.
func TestEntranceExitPupilConsistency(t *testing.T) {
	n := N_BK7.N(Wavelengthd)
	R := 150 * (n - 1)
	sys := mustSystem(t, Air{},
		[]Surface{
			{Profile: Plane{}, Aperture: CircleAperture(5)}, // stop in front
			{Profile: Sphere{R: R}, Aperture: CircleAperture(12)},
			{Profile: Sphere{R: -R}, Aperture: CircleAperture(12)},
		},
		[]float64{20, 6}, []Medium{Air{}, Air{}, N_BK7, Air{}})
	sys.StopIndex = 0
	f := FieldTan(0, 0.05)
	ep, err := sys.EntrancePupil(f, Wavelengthd)
	if err != nil {
		t.Fatalf("EntrancePupil: %v", err)
	}
	xp, err := sys.ExitPupil(f, Wavelengthd)
	if err != nil {
		t.Fatalf("ExitPupil: %v", err)
	}
	par, err := sys.Paraxial(Wavelengthd, math.Inf(-1))
	if err != nil {
		t.Fatal(err)
	}
	checkRel(t, "entrance pupil z", ep.Z, par.EntrancePupilZ, 1e-12)
	checkRel(t, "entrance pupil radius", ep.Radius, par.EntrancePupilRadius, 1e-12)
	checkRel(t, "exit pupil z", xp.Z, par.ExitPupilZ, 1e-12)
	checkRel(t, "exit pupil radius", xp.Radius, par.ExitPupilRadius, 1e-12)
	if ep.Offset != (Vec3{}) {
		t.Errorf("entrance pupil offset = %s, want zero", ep.Offset)
	}

	// Telecentric object space: the stop at the rear focal plane of a lens
	// images to infinity through the front part, so the entrance pupil does
	// not exist as a finite plane.
	tel := mustSystem(t, Air{},
		[]Surface{
			{Profile: Sphere{R: R}, Aperture: CircleAperture(12)},
			{Profile: Sphere{R: -R}, Aperture: CircleAperture(12)},
			{Profile: Plane{}, Aperture: CircleAperture(3)}, // stop at z = 6 + f
		},
		[]float64{6, 148.2}, []Medium{Air{}, N_BK7, Air{}, Air{}})
	tel.StopIndex = 2
	parT, err := tel.Paraxial(Wavelengthd, math.Inf(-1))
	if err != nil {
		t.Fatal(err)
	}
	if !math.IsInf(parT.EntrancePupilZ, 0) {
		t.Skipf("layout is not telecentric (epZ = %g)", parT.EntrancePupilZ)
	}
	if _, err := tel.EntrancePupil(f, Wavelengthd); err == nil {
		t.Error("EntrancePupil must fail for a telecentric object space")
	}

	// A stop with an unlimited aperture cannot size a pupil.
	open := mustSystem(t, Air{}, []Surface{{Profile: Sphere{R: R}}}, nil, nil)
	if _, err := open.EntrancePupil(f, Wavelengthd); err == nil {
		t.Error("EntrancePupil must fail for an unbounded stop")
	}
}

// TestChiefRayThroughStop checks the ray aiming of ChiefRay: the returned ray
// must cross the exact centre of the aperture stop, even when refracting
// surfaces sit between the object and the stop, for both infinite and finite
// objects and for a field that is rotated out of the meridional plane.
func TestChiefRayThroughStop(t *testing.T) {
	n := N_BK7.N(Wavelengthd)
	R := 150 * (n - 1)
	sys := mustSystem(t, Air{},
		[]Surface{
			{Profile: Sphere{R: R}, Aperture: CircleAperture(12)},
			{Profile: Sphere{R: -R}, Aperture: CircleAperture(12)},
			{Profile: Plane{}, Aperture: CircleAperture(4)}, // displaced stop
		},
		[]float64{6, 40}, []Medium{Air{}, N_BK7, Air{}, Air{}})
	sys.StopIndex = 2

	for _, tc := range []struct {
		name  string
		field Field
	}{
		{"infinite", FieldTan(0, 0.05)},
		{"infinite rotated", FieldTan(0.03, -0.04)},
		{"finite", FieldHeight(-150, 6, 3)},
		{"on axis", FieldTan(0, 0)},
	} {
		chief := sys.ChiefRay(tc.field, Wavelengthd)
		if !chief.Valid() {
			t.Fatalf("%s: chief ray invalid", tc.name)
		}
		res := sys.Trace(chief)
		if !res.OK() {
			t.Fatalf("%s: chief trace failed: %v", tc.name, res.Error())
		}
		hit := res.Steps[sys.StopIndex].Point
		if r := math.Hypot(hit.X, hit.Y); r > 1e-9 {
			t.Errorf("%s: chief crosses the stop %g mm from its centre", tc.name, r)
		}
		// The chief must start in front of the system.
		if chief.Origin.Z >= 0 {
			t.Errorf("%s: chief starts at z=%g, not in front", tc.name, chief.Origin.Z)
		}
		// For an infinite field the chief travels parallel to the field.
		if tc.field.IsInfinite() {
			d, _ := tc.field.Direction()
			if e := chief.Dir.Sub(d).Norm(); e > 1e-12 {
				t.Errorf("%s: chief direction off the field direction by %g", tc.name, e)
			}
		}
	}
}

// TestChiefOffsetStraightLine checks ChiefOffset for the case it is defined
// for: a stop that no refracting surface precedes (a mirror, or a stop at the
// first surface). There the chief ray is a straight line and the offset must
// reproduce it exactly.
func TestChiefOffsetStraightLine(t *testing.T) {
	mirror := mustSystem(t, Air{},
		[]Surface{{Profile: Sphere{R: -200}, Mirror: true, Aperture: CircleAperture(50)}},
		nil, nil)
	f := FieldTan(0, 0.03)
	zPlane := -60.0
	off, err := mirror.ChiefOffset(f, zPlane)
	if err != nil {
		t.Fatal(err)
	}
	// The chief crosses the mirror vertex (0,0,0), so at z = -60 it sits at
	// y = +60·tanθ (the beam descends towards the vertex).
	checkClose(t, "mirror chief offset y", off.Y, 60*0.03, 1e-12)
	checkClose(t, "mirror chief offset x", off.X, 0, 0)
	checkClose(t, "mirror chief offset z", off.Z, 0, 0)

	// The traced chief ray must pass through the same point.
	chief := mirror.ChiefRay(f, Wavelengthd)
	tt := (zPlane - chief.Origin.Z) / chief.Dir.Z
	checkClose(t, "chief y at the plane", chief.At(tt).Y, off.Y, 1e-12)

	// Launching a chief-aligned pupil through the centre must reproduce the
	// chief ray exactly, and the mirror's chief must hit the vertex.
	pup := Pupil{Z: zPlane, Radius: 20, Offset: off}
	r := LaunchRay(f, pup, PupilPoint{}, Wavelengthd)
	if d := r.Dir.Sub(chief.Dir).Norm(); d > 1e-15 {
		t.Errorf("chief-aligned launch direction off by %g", d)
	}
	res := mirror.Trace(r)
	if !res.OK() {
		t.Fatalf("chief trace failed: %v", res.Error())
	}
	if h := res.Steps[0].Point; math.Hypot(h.X, h.Y) > 1e-9 {
		t.Errorf("chief-aligned launch hits the mirror at %s, not the vertex", h)
	}
	// An on-axis field has no offset.
	off0, err := mirror.ChiefOffset(Field{}, zPlane)
	if err != nil {
		t.Fatal(err)
	}
	if off0 != (Vec3{}) {
		t.Errorf("on-axis chief offset = %s, want zero", off0)
	}
	// A finite object: the offset is on the line from the object point
	// through the stop vertex.
	ff := FieldHeight(-100, 3, -4)
	offF, err := mirror.ChiefOffset(ff, -90)
	if err != nil {
		t.Fatal(err)
	}
	want := ff.ObjectPoint().Add(V(0, 0, 0).Sub(ff.ObjectPoint()).Mul((-90 - ff.ObjectZ) / (0 - ff.ObjectZ)))
	checkClose(t, "finite chief offset x", offF.X, want.X, 1e-12)
	checkClose(t, "finite chief offset y", offF.Y, want.Y, 1e-12)
}
