package geo

import (
	"fmt"
	"math"

	"tetsuhiro-oz/coating"
)

// CoatedAmplitudes returns an AmplitudeFunc for a surface carrying the thin-film
// stack stk. The complex amplitudes come from the characteristic-matrix
// evaluation of package coating at the actual angle of incidence and wavelength
// of the step:
//
//	r_s = Rc(s),  t_s = Tc(s),  r_p = Rc(p),  t_p = Tc(p)·cosθ0/cosθs .
//
// coating reports the tangential amplitude ratio Tc of the transmitted wave; the
// factor cosθ0/cosθs converts it to the electric-field ratio used by this
// package. The reflection coefficients already follow the unified p convention,
// which the tests check against FresnelAmplitudes on a zero-thickness layer.
//
// The incidence medium of the stack must have the index in front of the surface
// and, for a refracting surface, the substrate the index behind it; a mismatch
// is reported as an error. The coating is evaluated for the lossless or lossy
// stack exactly as given, with the e^{-iωt} sign convention of package coating.
func CoatedAmplitudes(stk coating.Stack) AmplitudeFunc {
	return func(st Step, lambdaNm float64) (Amplitudes, error) {
		if err := stk.Validate(); err != nil {
			return Amplitudes{}, err
		}
		n0 := stk.Incident.IndexNm(lambdaNm)
		ns := stk.Substrate.IndexNm(lambdaNm)
		if math.Abs(real(n0)-st.NBefore) > 1e-9 {
			return Amplitudes{}, fmt.Errorf("geo: coating incidence index %.12g does not match the medium in front of the surface (%.12g)", real(n0), st.NBefore)
		}
		if !st.Mirror && math.Abs(real(ns)-st.NAfter) > 1e-9 {
			return Amplitudes{}, fmt.Errorf("geo: coating substrate index %.12g does not match the medium behind the surface (%.12g)", real(ns), st.NAfter)
		}
		angleDeg := st.IncAngle * 180 / math.Pi
		sS, err := stk.Evaluate(lambdaNm, angleDeg, coating.PolS)
		if err != nil {
			return Amplitudes{}, err
		}
		sP, err := stk.Evaluate(lambdaNm, angleDeg, coating.PolP)
		if err != nil {
			return Amplitudes{}, err
		}
		cosI := math.Abs(st.DirIn.Dot(st.Normal))
		sinI := math.Sin(st.IncAngle)
		// cosθs from Snell's law in the substrate, on the same branch as coating.
		cosS := polSqrt(ns*ns-n0*n0*complex(sinI*sinI, 0)) / ns
		if cosS == 0 {
			return Amplitudes{}, fmt.Errorf("geo: coating substrate cosine vanishes")
		}
		tp := sP.Tc * complex(cosI, 0) / cosS
		return Amplitudes{RS: sS.Rc, TS: sS.Tc, RP: sP.Rc, TP: tp}, nil
	}
}
