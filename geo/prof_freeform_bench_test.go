package geo

import (
	"testing"
)

// BenchmarkFreeformTrace compares the cost of tracing a spot diagram through
// the same physical off-axis parabolic mirror expressed three ways: as the
// analytic Conic, as an XYPolynomial with the same conic base (which reuses the
// analytic seeds for the intersection), and as a pure XY polynomial with no
// base (whose single Newton seed requires the coarse first-crossing check of
// the solver). Measured on amd64 (Hexapolar(16), 817 rays, 300 iterations):
// conic 1.5 µs/ray, freeform with base 1.9 µs/ray, pure freeform 12 µs/ray.
func BenchmarkFreeformTrace(b *testing.B) {
	const f, d, a = 100.0, 20.0, 8.0
	pure, err := NewXYPolynomial(0, 0, 0, []XYTerm{
		{I: 2, J: 0, Coeff: -1 / (4 * f)},
		{I: 0, J: 2, Coeff: -1 / (4 * f)},
	})
	if err != nil {
		b.Fatal(err)
	}
	based, err := NewXYPolynomial(-2*f, -1, 0, nil)
	if err != nil {
		b.Fatal(err)
	}
	sys := func(prof Profile) *System {
		s, err := NewSystem(Air{}, []Surface{{
			Profile: prof, Mirror: true, Aperture: CircleAperture(d + a + 2),
			Pose: posePtr(CoordBreak{Dx: d}),
		}}, nil, nil)
		if err != nil {
			b.Fatal(err)
		}
		return s
	}
	pup := Pupil{Z: -40, Radius: a}
	pts := Hexapolar(16)
	cases := []struct {
		name string
		sys  *System
	}{
		{"conic", sys(Conic{R: -2 * f, K: -1})},
		{"freeform-base-conic", sys(based)},
		{"freeform-pure", sys(pure)},
	}
	for _, tc := range cases {
		b.Run(tc.name, func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				if _, err := tc.sys.SpotDiagram(FieldTan(0, 0), pup, pts, -f, 550); err != nil {
					b.Fatal(err)
				}
			}
		})
	}
}
