package geo

import (
	"math"
	"testing"
)

func TestScratchSym(t *testing.T) {
	n := N_BK7.N(Wavelengthd)
	R := 150 * (n - 1)
	build := func(stop int) *System {
		sys := mustSystem(t, Air{},
			[]Surface{
				{Profile: Sphere{R: R}, Aperture: CircleAperture(12)},
				{Profile: Plane{}, Aperture: CircleAperture(4)},
				{Profile: Sphere{R: -R}, Aperture: CircleAperture(12)},
			},
			[]float64{3, 3},
			[]Medium{Air{}, N_BK7, N_BK7, Air{}})
		sys.StopIndex = stop
		return sys
	}
	sym := build(1)
	objZ, _ := magnificationMinusOne(sym, Wavelengthd)
	t.Logf("objZ=%.6f", objZ)
	for _, z := range []float64{objZ, objZ * 0.7, objZ * 0.5, -100, math.Inf(-1)} {
		p, err := sym.Paraxial(Wavelengthd, z)
		if err != nil {
			t.Fatal(err)
		}
		h := 0.5
		f := FieldHeight(z, 0, h)
		if math.IsInf(z, -1) {
			f = FieldTan(0, 0.005)
		}
		sd, err := sym.SeidelCoefficients(Wavelengthd, f, 4)
		if err != nil {
			t.Fatal(err)
		}
		t.Logf("objZ=%.3f mag=%.6g S1=%.6g S2=%.6g S5=%.6g  S2/S1=%.4g S5/S1=%.4g",
			z, p.Magnification, sd.S1, sd.S2, sd.S5, sd.S2/sd.S1, sd.S5/sd.S1)
	}
	// displaced stop at 1:1
	front := build(0)
	for _, stop := range []int{0, 2} {
		fs := build(stop)
		o2, _ := magnificationMinusOne(fs, Wavelengthd)
		sd, err := fs.SeidelCoefficients(Wavelengthd, FieldHeight(o2, 0, 0.5), 4)
		if err != nil {
			t.Fatal(err)
		}
		t.Logf("stop=%d objZ=%.4f S1=%.6g S2=%.6g S5=%.6g S2/S1=%.4g S5/S1=%.4g", stop, o2, sd.S1, sd.S2, sd.S5, sd.S2/sd.S1, sd.S5/sd.S1)
	}
	_ = front
}
