package geo

import (
	"math"
	"math/rand"
	"testing"
)

// TestSnellVectorToMachinePrecision verifies Snell's law n1 sin i = n2 sin t
// numerically for a large set of random incidence geometries. The sines are
// computed from vector cross products (|d × n| = sin θ, no cancellation), so
// the residual must be at the 1e-15..1e-16 level; the test requires < 1e-12.
// It also checks that the transmitted ray stays in the plane of incidence and
// is a unit vector.
func TestSnellVectorToMachinePrecision(t *testing.T) {
	rng := rand.New(rand.NewSource(20261015))
	const n1, n2 = 1.0, 1.5234
	worst := 0.0
	checked := 0
	for i := 0; i < 20000; i++ {
		d := randomUnit(rng)
		n := randomUnit(rng)
		cosi := math.Abs(d.Dot(n))
		if 1-cosi*cosi > (n2/n1)*(n2/n1) {
			continue // would be total internal reflection
		}
		out, ok := Refract(d, n, n1, n2)
		if !ok {
			t.Fatalf("Refract reported TIR for a geometry inside the critical angle")
		}
		if math.Abs(out.Norm()-1) > 1e-15 {
			t.Fatalf("transmitted direction is not a unit vector: |d'| = %.17g", out.Norm())
		}
		sinI := d.Cross(n).Norm()
		sinT := out.Cross(n).Norm()
		scale := math.Max(n1*sinI, n2*sinT)
		if scale < 1e-8 {
			continue // too close to normal incidence to say anything
		}
		checked++
		if err := math.Abs(n1*sinI-n2*sinT) / scale; err > worst {
			worst = err
		}
		// Coplanarity: d, n and d' must be coplanar, i.e. d' ⟂ (d × n).
		if axis := d.Cross(n); axis.Norm() > 1e-7 {
			if math.Abs(out.Dot(axis.Normalized())) > 1e-14 {
				t.Fatalf("transmitted ray left the plane of incidence: %.3g", math.Abs(out.Dot(axis.Normalized())))
			}
		}
	}
	if checked < 10000 {
		t.Fatalf("only %d geometries checked", checked)
	}
	if worst > 1e-12 {
		t.Fatalf("worst relative Snell residual %.3g exceeds 1e-12", worst)
	}
	t.Logf("worst relative Snell residual over %d random cases: %.3g", checked, worst)
}

// TestSnellReversibility checks optical reversibility: refracting the reversed
// transmitted ray must reproduce the incident ray. Geometries within 10° of the
// critical angle are excluded because the transmitted angle's sensitivity
// (dθ'/dθ → ∞) amplifies round-off; inside that range the check is still
// required to hold to 1e-10.
func TestSnellReversibility(t *testing.T) {
	rng := rand.New(rand.NewSource(7))
	worst := 0.0
	for i := 0; i < 5000; i++ {
		d := randomUnit(rng)
		n := randomUnit(rng)
		out, ok := Refract(d, n, 1.0, 1.72)
		if !ok {
			continue
		}
		back, ok := Refract(out.Neg(), n, 1.72, 1.0)
		if !ok {
			t.Fatalf("reverse refraction reported TIR")
		}
		e := back.Neg().Sub(d).Norm()
		if e > worst {
			worst = e
		}
	}
	if worst > 3e-11 {
		t.Fatalf("reversibility error %.3g too large", worst)
	}
	t.Logf("worst reversibility error: %.3g", worst)
}

// TestTotalInternalReflection checks the critical angle exactly (±1e-9 rad) and
// that a ray inside a glass plate beyond the critical angle stops with
// StatusTIR at the glass/air interface.
func TestTotalInternalReflection(t *testing.T) {
	const n1, n2 = 1.5, 1.0
	crit := math.Asin(n2 / n1)
	normal := Vec3{0, 0, 1}

	below := Vec3{math.Sin(crit - 1e-9), 0, math.Cos(crit - 1e-9)}
	if _, ok := Refract(below, normal, n1, n2); !ok {
		t.Fatalf("refraction refused just below the critical angle")
	}
	above := Vec3{math.Sin(crit + 1e-9), 0, math.Cos(crit + 1e-9)}
	if _, ok := Refract(above, normal, n1, n2); ok {
		t.Fatalf("refraction allowed just beyond the critical angle")
	}
	// At the critical angle the transmitted ray is exactly tangent.
	atCrit := Vec3{math.Sin(crit), 0, math.Cos(crit)}
	tang, ok := Refract(atCrit, normal, n1, n2)
	if !ok {
		t.Fatalf("critical angle should still transmit (grazing)")
	}
	checkClose(t, "grazing cos", math.Abs(tang.Dot(normal)), 0, 1e-9)

	// A ray inside the plate travelling beyond the critical angle must TIR at
	// the glass/air interface.
	ang := crit + 0.05
	dir := Vec3{0, -math.Sin(ang), math.Cos(ang)} // travels down and forward
	plate := plateSystem(t, n1, 10)
	sub, err := plate.Sub(1, 2) // only the glass/air interface
	if err != nil {
		t.Fatalf("Sub: %v", err)
	}
	res := sub.Trace(NewRay(Vec3{0, 5, -8}, dir, Wavelengthd))
	checkStatus(t, "TIR trace", res.Status, StatusTIR)
	if res.FailedAt != 0 {
		t.Errorf("TIR reported at surface %d, want 0", res.FailedAt)
	}
	if len(res.Steps) != 1 {
		t.Fatalf("expected 1 step, got %d", len(res.Steps))
	}
	// The sub-system keeps the absolute frame: the rear face of the 10 mm plate
	// sits at z = 10.
	checkClose(t, "TIR hit z", res.Steps[0].Point.Z, 10, 1e-12)
	checkClose(t, "TIR incidence", res.Steps[0].IncAngle, ang, 1e-12)
	// ... and its direction is unchanged by the failed refraction.
	checkClose(t, "TIR direction kept", res.Ray.Dir.Sub(dir).Norm(), 0, 1e-15)
}

// TestReflectionLaw verifies the law of reflection d' = d - 2(d·n)n and the
// equality of the incidence and reflection angles.
func TestReflectionLaw(t *testing.T) {
	rng := rand.New(rand.NewSource(99))
	worst := 0.0
	for i := 0; i < 5000; i++ {
		d := randomUnit(rng)
		n := randomUnit(rng)
		r := Reflect(d, n)
		if math.Abs(r.Norm()-1) > 1e-14 {
			t.Fatalf("reflected direction is not unit: %.17g", r.Norm())
		}
		ai := IncidenceAngle(d, n)
		ar := IncidenceAngle(r, n)
		if e := math.Abs(ai - ar); e > worst {
			worst = e
		}
		// Reflection must be independent of the normal's orientation.
		r2 := Reflect(d, n.Neg())
		if r.Sub(r2).Norm() > 1e-14 {
			t.Fatalf("reflection depends on the normal orientation")
		}
	}
	if worst > 1e-13 {
		t.Fatalf("worst reflection angle mismatch %.3g", worst)
	}
	// 45° fold mirror: a ray along +z hits a plane whose normal is tilted 45°
	// in the y-z plane and must leave along +y.
	dir := Vec3{0, 0, 1}
	nn := Vec3{0, -math.Sin(math.Pi / 4), math.Cos(math.Pi / 4)}
	out := Reflect(dir, nn)
	checkClose(t, "45° mirror out.x", out.X, 0, 1e-15)
	checkClose(t, "45° mirror out.y", out.Y, 1, 1e-15)
	checkClose(t, "45° mirror out.z", out.Z, 0, 1e-15)
}

// TestPlaneParallelPlateShift checks the exact lateral displacement of a ray
// that passes through a tilted plane-parallel plate: the exit point is shifted
// in y by T·tan t' with sin t' = sin i / n; the direction is preserved exactly
// and the optical path length is exact, too.
func TestPlaneParallelPlateShift(t *testing.T) {
	const (
		n = 1.5168
		T = 12.5
	)
	for _, deg := range []float64{0, 5, 15, 30, 45, 60, 75} {
		ang := deg * math.Pi / 180
		sinI := math.Sin(ang)
		cosI := math.Cos(ang)
		sinT := sinI / n
		cosT := math.Sqrt(1 - sinT*sinT)
		tant := sinT / cosT

		sys := plateSystem(t, n, T)
		start := Vec3{0, 0, -3}
		dir := Vec3{0, sinI, cosI}
		res := sys.Trace(NewRay(start, dir, Wavelengthd))
		if !res.OK() {
			t.Fatalf("deg=%g: trace failed: %v", deg, res.Error())
		}
		entry := res.Steps[0].Point
		exit := res.Steps[1].Point
		checkClose(t, "entry y", entry.Y, 3*sinI/cosI, 1e-9)
		checkClose(t, "entry z", entry.Z, 0, 1e-12)
		checkClose(t, "exit y", exit.Y, entry.Y+T*tant, 1e-9)
		checkClose(t, "exit z", exit.Z, T, 1e-12)
		if d := res.Ray.Dir.Sub(dir).Norm(); d > 1e-15 {
			t.Errorf("deg=%g: exit direction changed by %.3g", deg, d)
		}
		want := 3/cosI + n*T/cosT
		checkRel(t, "OPL", res.OPL(), want, 1e-12)
		checkRel(t, "air segment", res.Steps[0].T, 3/cosI, 1e-12)
		checkRel(t, "glass segment", res.Steps[1].T, T/cosT, 1e-12)
	}
}

// TestTraceStatuses exercises the status model: escape (missed geometry),
// vignetting (outside the clear aperture), backwards rays and a clean hit.
func TestTraceStatuses(t *testing.T) {
	sys := mustSystem(t, Air{},
		[]Surface{
			{Profile: Plane{}, Aperture: CircleAperture(50)},
			{Profile: Plane{}, Aperture: RectangleAperture{HalfWidth: 4, HalfHeight: 8}},
		},
		[]float64{20}, nil)

	res := sys.Trace(NewRay(Vec3{0, 0, -10}, Vec3{0, 0, 1}, Wavelengthd))
	checkStatus(t, "axial ray", res.Status, StatusOK)
	if len(res.Steps) != 2 {
		t.Fatalf("expected 2 steps, got %d", len(res.Steps))
	}
	checkClose(t, "axial hit y", res.Steps[0].Point.Y, 0, 1e-12)

	// Rectangle aperture: (3.9, 7.9) passes, (4.5, 0) is vignetted.
	res = sys.Trace(NewRay(Vec3{3.9, 7.9, -10}, Vec3{0, 0, 1}, Wavelengthd))
	if !res.OK() {
		t.Errorf("ray inside the rectangle aperture failed: %v", res.Error())
	}
	res = sys.Trace(NewRay(Vec3{4.5, 0, -10}, Vec3{0, 0, 1}, Wavelengthd))
	checkStatus(t, "rectangle vignetting", res.Status, StatusVignetted)
	if res.FailedAt != 1 {
		t.Errorf("rectangle vignetting reported at surface %d, want 1", res.FailedAt)
	}

	// Circle aperture on surface 0.
	circle := mustSystem(t, Air{},
		[]Surface{{Profile: Sphere{R: -40}, Aperture: CircleAperture(5)}}, nil, nil)
	res = circle.Trace(NewRay(Vec3{0, 6, -10}, Vec3{0, 0, 1}, Wavelengthd))
	checkStatus(t, "circle vignetting", res.Status, StatusVignetted)
	if res.FailedAt != 0 {
		t.Errorf("circle vignetting reported at surface %d, want 0", res.FailedAt)
	}

	// A ray travelling backwards never reaches the surfaces ahead.
	res = sys.Trace(NewRay(Vec3{0, 0, 30}, Vec3{0, 0, 1}, Wavelengthd))
	checkStatus(t, "backwards ray", res.Status, StatusEscape)

	// Degenerate directions are caught up front.
	res = sys.Trace(NewRay(Vec3{0, 0, -10}, Vec3{0, 0, 0}, Wavelengthd))
	checkStatus(t, "degenerate ray", res.Status, StatusDegenerate)
}

// TestAnnulusAperture checks the ring aperture.
func TestAnnulusAperture(t *testing.T) {
	a := AnnulusAperture{Inner: 2, Outer: 5}
	for _, tc := range []struct {
		x, y float64
		want bool
	}{{0, 0, false}, {1.9, 0, false}, {2, 0, true}, {3, 4, true}, {5.0001, 0, false}, {0, -5, true}} {
		if got := a.Contains(tc.x, tc.y); got != tc.want {
			t.Errorf("annulus.Contains(%g,%g) = %v, want %v", tc.x, tc.y, got, tc.want)
		}
	}
	checkClose(t, "annulus max radius", a.MaxRadius(), 5, 0)
	sys := mustSystem(t, Air{}, []Surface{{Profile: Plane{}, Aperture: a}}, nil, nil)
	if res := sys.Trace(NewRay(Vec3{0, 3, -5}, Vec3{0, 0, 1}, Wavelengthd)); !res.OK() {
		t.Errorf("ray in the annulus failed: %v", res.Error())
	}
	res := sys.Trace(NewRay(Vec3{0, 0.5, -5}, Vec3{0, 0, 1}, Wavelengthd))
	checkStatus(t, "ray in the hole", res.Status, StatusVignetted)
}

// randomUnit returns a random unit vector, sampled uniformly on the sphere.
func randomUnit(rng *rand.Rand) Vec3 {
	for {
		v := Vec3{rng.NormFloat64(), rng.NormFloat64(), rng.NormFloat64()}
		if n := v.Norm(); n > 1e-3 && n < 4 {
			return v.Div(n)
		}
	}
}
