package geo

import (
	"math"
	"testing"
)

// thinLensBK7 is the thin singlet of the field-curvature test: a biconvex BK7
// lens in air, R1 = +100 mm, R2 = −100 mm, centre thickness 0.05 mm (thin for
// the Petzval comparison), aperture stop on the front surface, semi-aperture
// 12.5 mm. Its paraxial focal length is about 96.76 mm.
func thinLensBK7(t *testing.T) *System {
	t.Helper()
	return mustSystem(t, Air{},
		[]Surface{
			{Profile: Sphere{R: 100}, Aperture: CircleAperture(12.5), Name: "front"},
			{Profile: Sphere{R: -100}, Aperture: CircleAperture(12.5), Name: "back"},
		},
		[]float64{0.05}, []Medium{Air{}, N_BK7, Air{}})
}

// fitSagPoly fits z = z0 + A·y² + B·y⁴ to the points (y_k, z_k) by least
// squares (normal equations, 3×3 Gauss–Jordan elimination) and returns z0, A, B.
func fitSagPoly(y, z []float64) (z0, a, b float64) {
	var m [3][4]float64
	for k := range y {
		basis := [3]float64{1, y[k] * y[k], y[k] * y[k] * y[k] * y[k]}
		for i := 0; i < 3; i++ {
			for j := 0; j < 3; j++ {
				m[i][j] += basis[i] * basis[j]
			}
			m[i][3] += basis[i] * z[k]
		}
	}
	for c := 0; c < 3; c++ {
		p := c
		for r := c + 1; r < 3; r++ {
			if math.Abs(m[r][c]) > math.Abs(m[p][c]) {
				p = r
			}
		}
		m[c], m[p] = m[p], m[c]
		for r := 0; r < 3; r++ {
			if r == c {
				continue
			}
			f := m[r][c] / m[c][c]
			for j := c; j < 4; j++ {
				m[r][j] -= f * m[c][j]
			}
		}
	}
	return m[0][3] / m[0][0], m[1][3] / m[1][1], m[2][3] / m[2][2]
}

// petzvalSumRadius returns the Petzval radius of a thin-or-thick refracting
// stack from its surface data, R_P = 1/Σ c_i·(1/n'_i − 1/n_i) (the exact
// third-order Petzval sum, independent of the thicknesses and of the stop).
// For a thin lens in air it reduces to R_P = −n·f.
func petzvalSumRadius(radii []float64, indices []float64) float64 {
	// indices[i] is the index before surface i; indices[i+1] after it.
	sum := 0.0
	for i, r := range radii {
		sum += (1 / r) * (1/indices[i+1] - 1/indices[i])
	}
	return 1 / sum
}

// TestFieldCurvatureThinLensPetzval checks the field curvature of a thin lens
// measured from real tangential and sagittal foci against its analytic Petzval
// curvature. Quantities (sags z(y) = z0 + A·y² + B·y⁴ fitted over the field,
// A = y² coefficient, curvature = 2A):
//
//	P    Petzval sag coefficient, A_P = 1/(2·R_P) with R_P the analytic thin-lens
//	     radius −n·f (n = N_BK7 at 587.5618 nm, f the paraxial focal length); the
//	     exact surface-sum radius of the same singlet agrees to 0.01 %.
//	The classical third-order structure of the two focal surfaces (the tangential
//	surface lies three times as far from the Petzval surface as the sagittal one,
//	z_T − z_P = 3(z_S − z_P), z_S − z_P = astigmatic sag) fixes the Petzval
//	surface from the two measured surfaces: A_P = (3·A_S − A_T)/2.
//
// Checks, all ≤ 1 % at the small fields used (θ = 0.02°…0.16°, image height
// ≈ 0.03…0.27 mm, pupil pair at ±1e-3 of the entrance-pupil radius):
//  1. the measured Petzval radius (3A_S − A_T)/2 against the analytic R_P;
//  2. the tangential sag predicted from the sagittal sag and the analytic
//     Petzval, A_T,pred = 3·A_S − 2·A_P, against the measured A_T;
//  3. the Seidel sums of geo (third order) at the smallest field: the Petzval
//     term −S4/(2·NA²·y²) against the analytic R_P, and the astigmatic term
//     −S3/(2·NA²·y²) against the measured A_S − A_P.
func TestFieldCurvatureThinLensPetzval(t *testing.T) {
	const (
		tol      = 1e-2
		rimFrac  = 1e-3
		fieldMax = 0.16
		nFields  = 8
	)
	sys := thinLensBK7(t)
	p, err := sys.Paraxial(Wavelengthd, math.Inf(-1))
	if err != nil {
		t.Fatalf("Paraxial: %v", err)
	}
	n := N_BK7.N(Wavelengthd)
	rThin := -n * p.EFL
	rExact := petzvalSumRadius([]float64{100, -100}, []float64{1, n, 1})
	aPAnalytic := 1 / (2 * rExact) // A_P of the exact surface-sum Petzval
	theta := make([]float64, nFields)
	for k := range theta {
		theta[k] = fieldMax * float64(k+1) / nFields
	}
	foci, err := sys.FieldFoci(theta, Wavelengthd, rimFrac)
	if err != nil {
		t.Fatalf("FieldFoci: %v", err)
	}
	yT, zT := make([]float64, nFields), make([]float64, nFields)
	yS, zS := make([]float64, nFields), make([]float64, nFields)
	for k, fp := range foci {
		yT[k], zT[k] = fp.YT, fp.ZT
		yS[k], zS[k] = fp.YS, fp.ZS
	}
	z0T, aT, bT := fitSagPoly(yT, zT)
	z0S, aS, bS := fitSagPoly(yS, zS)
	aPMeas := (3*aS - aT) / 2
	rPMeas := 1 / (2 * aPMeas)
	aTPred := 3*aS - 2*aPAnalytic

	// Third-order (Seidel) terms of geo at the smallest field.
	sd, err := sys.SeidelCoefficients(Wavelengthd, FieldAngle(0, theta[0]), p.EntrancePupilRadius)
	if err != nil {
		t.Fatalf("SeidelCoefficients: %v", err)
	}
	na2 := p.ImageNA * p.ImageNA
	y1 := foci[0].YS
	aPSeidel := -sd.S4 / (2 * na2 * y1 * y1)
	aAstigSeidel := -sd.S3 / (2 * na2 * y1 * y1)

	t.Logf("paraxial focus z0 (tangential %.9f, sagittal %.9f) vs Gaussian BFL focus %.9f",
		z0T, z0S, p.BackFocalZ)
	t.Logf("fit: tangential A_T=%.9g B_T=%.6g; sagittal A_S=%.9g B_S=%.6g", aT, bT, aS, bS)
	t.Logf("analytic: R_P(thin, −n·f)=%.9g mm; R_P(exact surface sum)=%.9g mm; f=%.9g mm, n=%.9g",
		rThin, rExact, p.EFL, n)
	t.Logf("1. Petzval radius from the measured surfaces 1/(2·(3A_S−A_T)/2)=%.9g mm; rel err vs exact %.3e, vs thin-lens %.3e",
		rPMeas, math.Abs(rPMeas-rExact)/math.Abs(rExact), math.Abs(rPMeas-rThin)/math.Abs(rThin))
	t.Logf("2. A_T measured %.9g, predicted 3A_S−2A_P %.9g; rel err %.3e", aT, aTPred, math.Abs(aT-aTPred)/math.Abs(aT))
	t.Logf("3. Seidel Petzval R_P=%.9g mm (rel err %.3e vs exact); Seidel astig −S3/(2NA²y²)=%.9g vs measured A_S−A_P=%.9g (rel err %.3e)",
		1/(2*aPSeidel), math.Abs(1/(2*aPSeidel)-rExact)/math.Abs(rExact),
		aAstigSeidel, aS-aPMeas, math.Abs(aAstigSeidel-(aS-aPMeas))/math.Abs(aAstigSeidel))

	if d := math.Abs(rPMeas-rExact) / math.Abs(rExact); d > tol {
		t.Errorf("measured Petzval radius %.9g mm vs analytic %.9g mm: rel %.3e > %.0e", rPMeas, rExact, d, tol)
	}
	if d := math.Abs(rPMeas-rThin) / math.Abs(rThin); d > tol {
		t.Errorf("measured Petzval radius %.9g mm vs thin-lens −n·f %.9g mm: rel %.3e > %.0e", rPMeas, rThin, d, tol)
	}
	if d := math.Abs(aT-aTPred) / math.Abs(aT); d > tol {
		t.Errorf("tangential sag coefficient %.9g vs prediction 3A_S−2A_P %.9g: rel %.3e > %.0e", aT, aTPred, d, tol)
	}
	if d := math.Abs(1/(2*aPSeidel)-rExact) / math.Abs(rExact); d > tol {
		t.Errorf("Seidel Petzval radius %.9g mm vs analytic %.9g mm: rel %.3e > %.0e", 1/(2*aPSeidel), rExact, d, tol)
	}
	if d := math.Abs(aAstigSeidel-(aS-aPMeas)) / math.Abs(aAstigSeidel); d > tol {
		t.Errorf("Seidel astigmatic term %.9g vs measured %.9g: rel %.3e > %.0e", aAstigSeidel, aS-aPMeas, d, tol)
	}
	if d := math.Abs(z0T - p.BackFocalZ); d > 1e-5 {
		t.Errorf("tangential focus on axis %.9f differs from the Gaussian focus %.9f by %.3e mm", z0T, p.BackFocalZ, d)
	}
}
