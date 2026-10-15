package geo

import (
	"math"
	"testing"

	"tetsuhiro-oz/coating"
)

// polInterfaceStep returns a trace-like step for a plane interface at z = 0 with
// the incidence angle theta (rad) in the y–z plane, from index n1 into n2.
func polInterfaceStep(t *testing.T, theta, n1, n2 float64) Step {
	t.Helper()
	din := V(0, math.Sin(theta), math.Cos(theta))
	dout, ok := Refract(din, AxisZ, n1, n2)
	if !ok {
		t.Fatalf("total internal reflection at θ = %g", theta)
	}
	return Step{DirIn: din, Normal: AxisZ, DirOut: dout, NBefore: n1, NAfter: n2,
		IncAngle: theta, CosInc: math.Cos(theta), Status: StatusOK}
}

// polStackConst builds a coating stack from constant (lossless or absorbing)
// indices, for the coated-surface tests.
func polStackConst(nIn, nSub complex128, layers ...coating.Layer) coating.Stack {
	in := coating.NewMaterial("incident", coating.ConstIndex{N: real(nIn), K: imag(nIn)})
	sub := coating.NewMaterial("substrate", coating.ConstIndex{N: real(nSub), K: imag(nSub)})
	return coating.NewStack(in, sub, layers...)
}

// TestPolCoatingZeroThicknessMatchesFresnel validates the coated-surface model in
// the convention of this package: a stack whose only layer has zero thickness is
// the bare interface, so CoatedAmplitudes must reproduce FresnelStep for s and p
// at every angle. This is also the check of the sign of r_p against the
// characteristic-matrix code.
func TestPolCoatingZeroThicknessMatchesFresnel(t *testing.T) {
	stk := polStackConst(1, 1.52, coating.Layer{Material: coating.NewMaterial("film", coating.ConstIndex{N: 1.38}), ThicknessNm: 0})
	f := CoatedAmplitudes(stk)
	worst := 0.0
	for k := 0; k <= 16; k++ {
		theta := 1.39 * float64(k) / 16
		st := polInterfaceStep(t, theta, 1, 1.52)
		a, err := f(st, 550)
		if err != nil {
			t.Fatalf("θ=%g: %v", theta, err)
		}
		b, err := FresnelStep(st, 550)
		if err != nil {
			t.Fatal(err)
		}
		for _, d := range []complex128{a.RS - b.RS, a.TS - b.TS, a.RP - b.RP, a.TP - b.TP} {
			if x := cmplxAbsC(d); x > worst {
				worst = x
			}
		}
	}
	if worst > 1e-12 {
		t.Fatalf("coated model with a zero-thickness layer differs from Fresnel by %.3g", worst)
	}
	t.Logf("zero-thickness stack vs bare Fresnel interface: worst |Δ| = %.3g", worst)
}

// cmplxAbsC returns |z|.
func cmplxAbsC(z complex128) float64 {
	return math.Hypot(real(z), imag(z))
}

// TestPolCoatingLosslessEnergy checks a lossless quarter-wave antireflection
// layer (n = 1.38 on glass 1.52, λ0 = 550 nm) at oblique incidence:
// |r|² + (n_s·cosθ_s)/(n0·cosθ0)·|t|² = 1 for both polarizations, with the
// electric-field t_p produced by CoatedAmplitudes. The transmittance is also
// compared with the one reported by package coating.
func TestPolCoatingLosslessEnergy(t *testing.T) {
	const n0, ns, nf, lambda = 1.0, 1.52, 1.38, 550.0
	layer := coating.Layer{Material: coating.NewMaterial("MgF2-like", coating.ConstIndex{N: nf}), ThicknessNm: lambda / (4 * nf)}
	stk := polStackConst(n0, ns, layer)
	f := CoatedAmplitudes(stk)
	worstE, worstTMM := 0.0, 0.0
	for k := 0; k <= 14; k++ {
		theta := 1.2 * float64(k) / 14
		st := polInterfaceStep(t, theta, n0, ns)
		a, err := f(st, lambda)
		if err != nil {
			t.Fatal(err)
		}
		cosI := math.Cos(theta)
		cosS := math.Sqrt(1 - (n0*math.Sin(theta)/ns)*(n0*math.Sin(theta)/ns))
		fac := ns * cosS / (n0 * cosI)
		for _, pr := range []struct {
			r, tt complex128
			pol   coating.Polarization
		}{{a.RS, a.TS, coating.PolS}, {a.RP, a.TP, coating.PolP}} {
			R := real(pr.r * cmplxConj(pr.r))
			T := fac * real(pr.tt*cmplxConj(pr.tt))
			if d := math.Abs(R + T - 1); d > worstE {
				worstE = d
			}
			smp, err := stk.Evaluate(lambda, theta*180/math.Pi, pr.pol)
			if err != nil {
				t.Fatal(err)
			}
			if d := math.Abs(R - smp.R); d > worstTMM {
				worstTMM = d
			}
			if d := math.Abs(T - smp.T); d > worstTMM {
				worstTMM = d
			}
		}
	}
	if worstE > 1e-12 {
		t.Fatalf("lossless coated interface: energy residual %.3g > 1e-12", worstE)
	}
	if worstTMM > 1e-12 {
		t.Fatalf("coated interface: power differs from coating.Evaluate by %.3g > 1e-12", worstTMM)
	}
	t.Logf("lossless AR-coated interface: energy residual %.3g; |ΔR|,|ΔT| vs coating ≤ %.3g", worstE, worstTMM)
}

// cmplxConj returns the complex conjugate.
func cmplxConj(z complex128) complex128 {
	return complex(real(z), -imag(z))
}

// TestPolCoatedPlateTraceTransmittance runs a coated face through the trace: a
// plane-parallel plate whose front face carries the AR stack (air to glass) and
// whose back face is bare. For a pure s or pure p input the path transmittance
// of the single traversal must equal the product of the two faces, the front one
// taken from coating.Evaluate and the back one from the Fresnel formula.
func TestPolCoatedPlateTraceTransmittance(t *testing.T) {
	const ns, lambda = 1.52, 550.0
	layer := coating.Layer{Material: coating.NewMaterial("film", coating.ConstIndex{N: 1.38}), ThicknessNm: lambda / (4 * 1.38)}
	stk := polStackConst(1, ns, layer)
	sys, err := NewSystem(Air{}, []Surface{NewSurface(Plane{}, 0), NewSurface(Plane{}, 0)},
		[]float64{4}, []Medium{Air{}, ConstantIndex(ns), Air{}})
	if err != nil {
		t.Fatal(err)
	}
	worst := 0.0
	for _, theta := range []float64{0.2, 0.7, 1.1} {
		d := V(0, math.Sin(theta), math.Cos(theta))
		ray := NewRay(V(0.1, -1, -3), d, lambda)
		opt := &PolOptions{Amplitude: []AmplitudeFunc{CoatedAmplitudes(stk), nil}}
		for _, pol := range []struct {
			in  [2]complex128
			tmm coating.Polarization
		}{{[2]complex128{1, 0}, coating.PolS}, {[2]complex128{0, 1}, coating.PolP}} {
			res, err := sys.TracePolarized(ray, pol.in, opt)
			if err != nil {
				t.Fatal(err)
			}
			smp, err := stk.Evaluate(lambda, theta*180/math.Pi, pol.tmm)
			if err != nil {
				t.Fatal(err)
			}
			// Back face: glass to air, Fresnel. The angle inside the glass follows from Snell.
			sinG := math.Sin(theta) / ns
			cosG := math.Sqrt(1 - sinG*sinG)
			a, cosAir, err := FresnelAmplitudes(ns, 1, cosG)
			if err != nil {
				t.Fatal(err)
			}
			tBack := 0.0
			if pol.tmm == coating.PolS {
				tBack = real(a.TS * cmplxConj(a.TS))
			} else {
				tBack = real(a.TP * cmplxConj(a.TP))
			}
			backFac := real(cosAir) / (ns * cosG)
			want := smp.T * tBack * backFac
			if x := math.Abs(res.Transmittance - want); x > worst {
				worst = x
			}
		}
	}
	if worst > 1e-12 {
		t.Fatalf("coated plate: path transmittance differs from the face product by %.3g > 1e-12", worst)
	}
	t.Logf("coated plate traced transmittance vs face product: worst |Δ| = %.3g", worst)
}

// TestPolCoatingMismatchRejected: the stack must sit on the media of the trace.
func TestPolCoatingMismatchRejected(t *testing.T) {
	stk := polStackConst(1.2, 1.52)
	st := polInterfaceStep(t, 0.3, 1, 1.52)
	if _, err := CoatedAmplitudes(stk)(st, 550); err == nil {
		t.Error("an incidence index that differs from the trace must be rejected")
	}
	stk2 := polStackConst(1, 1.7)
	if _, err := CoatedAmplitudes(stk2)(st, 550); err == nil {
		t.Error("a substrate index that differs from the trace must be rejected")
	}
}
