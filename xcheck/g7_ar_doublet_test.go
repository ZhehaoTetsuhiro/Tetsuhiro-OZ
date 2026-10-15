package xcheck

import (
	"math"
	"testing"

	"tetsuhiro-oz/coating"
	"tetsuhiro-oz/geo"
	"tetsuhiro-oz/lens"
)

// G7 joint check (plan §4.2): an AR-coated achromat doublet at its design
// wavelength and F/#. The joint evaluation (lens.Prescription.Transmittance,
// which traces every ray and evaluates each interface's coating with the TMM
// at the ray's own incidence angle and polarization) must agree with a
// separate, independent computation to ≤ 1e-9.
//
// The independent computation re-derives everything from the geo and coating
// packages alone, without going through the lens joint path:
//
//  1. launch the sample rays from first principles (parallel rays arriving at
//     the entrance pupil for an object at infinity);
//  2. trace each ray with geo.System.Trace and read each interface's
//     incidence angle from the trace steps;
//  3. build each interface's thin-film stack directly (coating.NewMaterial +
//     coating.Stack) and evaluate it with coating.Stack.Evaluate;
//  4. multiply the interface transmittances along the ray and average.
//
// Physical sanity is logged: the AR coating must lift the system
// transmittance, and the marginal rays must reach several degrees of
// incidence (the F/5 pupil), so the check really exercises off-normal TMM.

// geoDisp adapts a geo.Medium to the coating package's Dispersion interface.
type geoDisp struct{ m geo.Medium }

func (d geoDisp) IndexNm(lambdaNm float64) complex128 { return complex(d.m.N(lambdaNm), 0) }
func (d geoDisp) Lossless() bool                      { return true }
func (d geoDisp) Model() string                       { return "xcheck geo medium " + d.m.Name() }

func TestG7ARCoatedDoubletTransmittance(t *testing.T) {
	p, _, err := lens.DesignAchromat(100, "N-BK7", "F2", 10, nil)
	if err != nil {
		t.Fatal(err)
	}
	sys, err := p.System()
	if err != nil {
		t.Fatal(err)
	}
	par, err := p.Paraxial(geo.Wavelengthd)
	if err != nil {
		t.Fatal(err)
	}

	// MgF2 quarter-wave AR at the d-line on the two air–glass interfaces.
	qw, err := coating.QWOTNm(coating.MgF2, geo.Wavelengthd)
	if err != nil {
		t.Fatal(err)
	}
	ar := lens.Coating{Name: "MgF2 QWOT", Layers: []coating.Layer{{Material: coating.MgF2, ThicknessNm: qw}}}
	n := len(sys.Surfaces)
	coatings := make([]lens.Coating, n)
	coatings[0] = ar
	coatings[n-1] = ar

	fields := []lens.Field{{TanX: 0, TanY: 0}, {TanX: 0, TanY: 0.0087}} // on-axis and 0.5°
	samples := lens.UnitDiskSamples(3, 6)
	pols := []coating.Polarization{coating.PolS, coating.PolP, coating.PolUnpolarized}

	// Joint path.
	tab, err := p.Transmittance(lens.TransmittanceSpec{
		Wavelengths: []float64{geo.Wavelengthd},
		Fields:      fields,
		Samples:     samples,
		Pols:        pols,
		Coatings:    coatings,
	})
	if err != nil {
		t.Fatal(err)
	}
	// Uncoated reference (same geometry, empty coatings).
	tabBare, err := p.Transmittance(lens.TransmittanceSpec{
		Wavelengths: []float64{geo.Wavelengthd},
		Fields:      fields,
		Samples:     samples,
		Pols:        pols,
	})
	if err != nil {
		t.Fatal(err)
	}

	// Independent computation from the geo + coating packages alone.
	zp := par.EntrancePupilZ
	rp := par.EntrancePupilRadius
	zStop := par.VertexZ[par.StopIndex]
	var maxIncDeg float64
	for fi, fld := range fields {
		sumS, sumP, nPass := 0.0, 0.0, 0
		for _, smp := range samples {
			// Rays from first principles: parallel to the field direction,
			// arriving at the entrance pupil at height (Px, Py)·rp.
			dir, ok := geo.V(-fld.TanX, -fld.TanY, 1).Unit()
			if !ok {
				t.Fatal("degenerate field direction")
			}
			dz := zp - zStop
			ch := geo.V(dir.X/dir.Z*dz, dir.Y/dir.Z*dz, 0)
			origin := geo.V(ch.X+smp.Px*rp, ch.Y+smp.Py*rp, zp)
			res := sys.Trace(geo.NewRay(origin, dir, geo.Wavelengthd))
			if !res.OK() {
				continue
			}
			ts, tp := 1.0, 1.0
			for _, st := range res.Steps {
				if st.Status != geo.StatusOK {
					continue
				}
				before := sys.Mediums[st.Surface]
				after := sys.Mediums[st.Surface+1]
				stack := coating.Stack{
					Incident:  coating.NewMaterial(before.Name(), geoDisp{before}),
					Substrate: coating.NewMaterial(after.Name(), geoDisp{after}),
					Layers:    coatings[st.Surface].Layers,
				}
				angleDeg := st.IncAngle * 180 / math.Pi
				if angleDeg > maxIncDeg {
					maxIncDeg = angleDeg
				}
				sS, err := stack.Evaluate(geo.Wavelengthd, angleDeg, coating.PolS)
				if err != nil {
					t.Fatal(err)
				}
				sP, err := stack.Evaluate(geo.Wavelengthd, angleDeg, coating.PolP)
				if err != nil {
					t.Fatal(err)
				}
				ts *= sS.T
				tp *= sP.T
			}
			nPass++
			sumS += ts
			sumP += tp
		}
		if nPass == 0 {
			t.Fatalf("field %d: no ray traversed the system", fi)
		}
		manual := []float64{sumS / float64(nPass), sumP / float64(nPass), 0.5 * (sumS + sumP) / float64(nPass)}
		for k, pol := range pols {
			joint := tab.T[0][fi][k]
			if d := math.Abs(joint - manual[k]); d > 1e-9 {
				t.Errorf("field %d pol %v: joint T = %.15f, independent computation = %.15f (|diff| %.3g > 1e-9)",
					fi, pol, joint, manual[k], d)
			}
			t.Logf("field %.2f° pol %v: joint T = %.10f, independent = %.10f, |diff| = %.3g", math.Atan(fields[fi].TanY)*180/math.Pi, pol, joint, manual[k], math.Abs(joint-manual[k]))
		}
		if tabBare.T[0][fi][0] >= tab.T[0][fi][0] {
			t.Errorf("field %d: AR coating did not raise the s-pol transmittance (bare %.6f, coated %.6f)",
				fi, tabBare.T[0][fi][0], tab.T[0][fi][0])
		}
		t.Logf("field %.2f°: s-pol transmittance bare %.6f -> AR-coated %.6f; throughput coated %.6f",
			math.Atan(fields[fi].TanY)*180/math.Pi, tabBare.T[0][fi][0], tab.T[0][fi][0], tab.Throughput[0][fi][0])
	}
	if maxIncDeg < 4 {
		t.Errorf("max incidence angle on a coated interface %.2f° < 4°: the F/5 aperture was not exercised", maxIncDeg)
	}
	t.Logf("max incidence angle over all traced interface hits: %.3f° (F/5 pupil)", maxIncDeg)
}
