package geo

import (
	"fmt"
	"math"
)

// QForm selects the Forbes Q-type departure.
type QForm int

// Forbes Q-type forms.
const (
	// QCon is the conic-based form: the base is the conic (R, K) and the
	// departure is a sag term, Forbes (2010) Optics Express 18(13) 13851,
	// Eq. (5.1)–(5.3):
	//
	//	z(ρ) = cρ²/(1+φ) + u⁴ Σ_m s_m Q_m^con(u²),  φ = sqrt(1−(1+K)c²ρ²),
	//	u = ρ/ρ_max,  Q_m^con(x) = P_m^(0,4)(2x−1)  (Jacobi, P_m(1) = 1).
	QCon QForm = iota
	// QBfs is the best-fit-sphere form: the base is a sphere of radius R
	// (K = 0) and the departure is measured along the sphere normal, Forbes
	// (2007) Optics Express 15(8) 5218, Eq. (2.1), with the orthonormal basis of
	// Forbes (2010) Optics Express 18(19) 19700:
	//
	//	z(ρ) = cρ²/(1+sqrt(1−c²ρ²)) + u²(1−u²)/sqrt(1−c²ρ²) Σ_m a_m Q_m(u²).
	QBfs
)

// String implements fmt.Stringer ("qcon" / "qbfs").
func (f QForm) String() string {
	switch f {
	case QCon:
		return "qcon"
	case QBfs:
		return "qbfs"
	}
	return fmt.Sprintf("qform(%d)", int(f))
}

// QAsphere is an axially symmetric Forbes Q-type asphere (see QForm). The
// coefficients s_m (QCon) or a_m (QBfs), m = 0, 1, …, multiply the orthogonal
// basis members in order; a zero coefficient list (or all-zero coefficients)
// gives the bare base surface, which is Conic{R, K} for QCon and Sphere{R} for
// QBfs with K = 0.
//
// NormR is the normalisation radius ρ_max (mm). Forbes uses the semi-aperture
// CA/2. It must be positive whenever a coefficient is non-zero.
//
// Verification status of the formulas: the Qcon and Qbfs definitions above are
// taken from Forbes (2010, OE 18(13) 13851 and OE 18(19) 19700); the Qbfs
// recurrence constants are checked against the printed matrix of Forbes (2010)
// Eq. (A.17) and the slope orthonormality of Eq. (2.2) in prof_qasphere_test.go.
// The Qbfs form with a conic (non-sphere) base, which Forbes (2011, OE 19(10)
// 9923) introduced as a generalisation, is NOT implemented here (待核对: its
// exact formula was not verified); QBfs requires K = 0.
type QAsphere struct {
	Form   QForm
	R      float64   // base radius (mm): conic vertex radius (QCon) or best-fit sphere radius (QBfs)
	K      float64   // conic constant of the base; must be 0 for QBfs
	NormR  float64   // normalisation radius ρ_max (mm)
	Coeffs []float64 // s_m (QCon) or a_m (QBfs), m = 0, 1, …
}

// Validate reports structural errors: an unknown form, NaN parameters, QBfs
// with a non-zero K, non-finite coefficients, and a non-positive normalisation
// radius while coefficients are active.
func (p QAsphere) Validate() error {
	if p.Form != QCon && p.Form != QBfs {
		return fmt.Errorf("geo: q asphere: unknown form %v", p.Form)
	}
	if math.IsNaN(p.R) || math.IsNaN(p.K) || math.IsNaN(p.NormR) {
		return fmt.Errorf("geo: q asphere %v: NaN parameter (R=%v K=%v norm=%v)", p.Form, p.R, p.K, p.NormR)
	}
	if p.Form == QBfs && p.K != 0 {
		return fmt.Errorf("geo: q asphere qbfs: base is the best-fit sphere, K must be 0 (got %g)", p.K)
	}
	for i, a := range p.Coeffs {
		if math.IsNaN(a) || math.IsInf(a, 0) {
			return fmt.Errorf("geo: q asphere %v: coefficient %d is not finite (%v)", p.Form, i, a)
		}
	}
	if p.active() && !(p.NormR > 0) {
		return fmt.Errorf("geo: q asphere %v: normalisation radius must be > 0 (got %g)", p.Form, p.NormR)
	}
	return nil
}

// active reports whether any departure coefficient is non-zero.
func (p QAsphere) active() bool {
	for _, a := range p.Coeffs {
		if a != 0 {
			return true
		}
	}
	return false
}

// baseK is the conic constant of the base used by the sag (0 for QBfs).
func (p QAsphere) baseK() float64 {
	if p.Form == QBfs {
		return 0
	}
	return p.K
}

// Sag implements Profile.
func (p QAsphere) Sag(r float64) float64 {
	if p.Form == QBfs && p.K != 0 {
		return math.NaN()
	}
	z := conicBaseSag(r, p.R, p.baseK())
	if math.IsNaN(z) || !p.active() {
		return z
	}
	if !(p.NormR > 0) {
		return math.NaN()
	}
	u := r / p.NormR
	x := u * u
	switch p.Form {
	case QCon:
		s, _ := qconSum(p.Coeffs, x)
		return z + x*x*s
	case QBfs:
		psi2 := 1 - baseCurvature(p.R)*baseCurvature(p.R)*r*r
		if psi2 < 0 {
			return math.NaN()
		}
		s, _ := qbfsSum(p.Coeffs, x)
		return z + x*(1-x)*s/math.Sqrt(psi2)
	}
	return math.NaN()
}

// Slope implements Profile: the analytic dz/dr of the sag above.
//
//	QCon:  dz/dr = conic' + (2u³/ρ_max)·[2S(x) + x·S'(x)],          x = u²,
//	QBfs:  dz/dr = sphere' + (2u/ρ_max)·F'(x)/ψ + c²r/ψ³ · F(x),
//	       F(x) = x(1−x)S(x), ψ = sqrt(1 − c²r²), c = 1/R.
func (p QAsphere) Slope(r float64) float64 {
	if p.Form == QBfs && p.K != 0 {
		return math.NaN()
	}
	s0 := conicBaseSlope(r, p.R, p.baseK())
	if math.IsNaN(s0) || !p.active() {
		return s0
	}
	if !(p.NormR > 0) {
		return math.NaN()
	}
	u := r / p.NormR
	x := u * u
	switch p.Form {
	case QCon:
		s, ds := qconSum(p.Coeffs, x)
		return s0 + (2*u*x/p.NormR)*(2*s+x*ds)
	case QBfs:
		c := baseCurvature(p.R)
		psi2 := 1 - c*c*r*r
		if psi2 <= 0 {
			return math.NaN()
		}
		psi := math.Sqrt(psi2)
		s, ds := qbfsSum(p.Coeffs, x)
		fx := x * (1 - x) * s
		fpx := (1-2*x)*s + x*(1-x)*ds
		return s0 + (2*u/p.NormR)*fpx/psi + c*c*r*fx/(psi*psi2)
	}
	return math.NaN()
}

// Curvature implements Profile: the paraxial vertex curvature.
//
// QCon: the departure u⁴(…) has no r² term, so the vertex curvature is the base
// curvature 1/R.
//
// QBfs: the departure u²(1−u²)/ψ·S(u²) contributes S(0)·r²/ρ_max², so
//
//	c_eff = c + 2·S(0)/ρ_max² = c + 2·Σ_m a_m Q_m(0)/ρ_max² .
func (p QAsphere) Curvature() float64 {
	c := baseCurvature(p.R)
	if p.Form == QBfs && p.active() && p.NormR > 0 {
		s0, _ := qbfsSum(p.Coeffs, 0)
		c += 2 * s0 / (p.NormR * p.NormR)
	}
	return c
}

// Intersect implements Profile: conic-base analytic seeds, Newton on the exact
// sag equation, polish to machine precision.
func (p QAsphere) Intersect(o, d Vec3, tMin, rMax float64) (float64, bool) {
	if p.Form == QBfs && p.K != 0 {
		return 0, false
	}
	hasBase := !math.IsInf(p.R, 0) && p.R != 0
	t, ok := solveIntersection(p, p.R, p.baseK(), hasBase, o, d, tMin, rMax)
	if !ok {
		return 0, false
	}
	return polishIntersection(p, o, d, t, tMin), true
}

// String implements Profile.
func (p QAsphere) String() string {
	return fmt.Sprintf("q asphere %v R=%g K=%g norm=%g (%d coeffs)", p.Form, p.R, p.K, p.NormR, len(p.Coeffs))
}

// ---------------------------------------------------------------------------
// Basis polynomials
// ---------------------------------------------------------------------------

// qconABC returns the coefficients of the Jacobi (α=0, β=4) three-term
// recurrence in t = 2x−1 for the step to order n+1 (Forbes 2010 Eq. (5.2),
// written in x = u²; here in t, the standard form):
//
//	P_{n+1}(t) = (A_n t + B_n) P_n(t) − C_n P_{n−1}(t),
//	A_n = (2n+5)(n+3) / ((n+1)(n+5)),
//	B_n = −4(2n+5) / ((n+1)(n+2)(n+5)),
//	C_n = n(n+3)(n+4) / ((n+1)(n+2)(n+5)).
func qconABC(n int) (A, B, C float64) {
	fn := float64(n)
	A = (2*fn + 5) * (fn + 3) / ((fn + 1) * (fn + 5))
	B = -4 * (2*fn + 5) / ((fn + 1) * (fn + 2) * (fn + 5))
	C = fn * (fn + 3) * (fn + 4) / ((fn + 1) * (fn + 2) * (fn + 5))
	return A, B, C
}

// qconBasis fills z[m] = Q_m^con(x) = P_m^(0,4)(2x−1) and dz[m] = dQ_m^con/dx
// for m = 0 … len(z)−1, by the three-term recurrence and its t-derivative.
func qconBasis(x float64, z, dz []float64) {
	n := len(z)
	if n == 0 {
		return
	}
	t := 2*x - 1
	z[0], dz[0] = 1, 0
	if n == 1 {
		return
	}
	// P_1 = 3t − 2, dP_1/dt = 3.
	pm2, dpm2 := 1.0, 0.0
	pm1, dpm1 := 3*t-2, 3.0
	z[1], dz[1] = pm1, 2*dpm1
	for m := 2; m < n; m++ {
		A, B, C := qconABC(m - 1)
		lin := A*t + B
		pm := lin*pm1 - C*pm2
		dpm := A*pm1 + lin*dpm1 - C*dpm2
		z[m], dz[m] = pm, 2*dpm
		pm2, dpm2 = pm1, dpm1
		pm1, dpm1 = pm, dpm
	}
}

// qconSum returns S(x) = Σ_m s_m Q_m^con(x) and dS/dx at x = u².
func qconSum(s []float64, x float64) (sum, dsum float64) {
	n := len(s)
	if n == 0 {
		return 0, 0
	}
	z := make([]float64, n)
	dz := make([]float64, n)
	qconBasis(x, z, dz)
	for m, c := range s {
		sum += c * z[m]
		dsum += c * dz[m]
	}
	return sum, dsum
}

// qbfsConstants returns the recurrence constants f_m (m < n), g_m (m < n−1),
// h_m (m < n−2) of Forbes (2010) Appendix A:
//
//	f_0 = 2,  f_1 = √19/2,  g_0 = −1/2,
//	h_{m−2} = −m(m−1) / (2 f_{m−2}),                           (A.14)
//	g_{m−1} = −(1 + g_{m−2} h_{m−2}) / f_{m−1},                (A.15)
//	f_m     = sqrt(m(m+1) + 3 − g_{m−1}² − h_{m−2}²),          (A.16)
//
// evaluated in this order for m = 2, 3, …
func qbfsConstants(n int) (f, g, h []float64) {
	if n <= 0 {
		return nil, nil, nil
	}
	f = make([]float64, n)
	f[0] = 2
	if n > 1 {
		f[1] = math.Sqrt(19) / 2
		g = make([]float64, n-1)
		g[0] = -0.5
	}
	if n > 2 {
		h = make([]float64, n-2)
	}
	for m := 2; m < n; m++ {
		hm2 := -float64(m*(m-1)) / (2 * f[m-2])
		h[m-2] = hm2
		gm1 := -(1 + g[m-2]*hm2) / f[m-1]
		g[m-1] = gm1
		f[m] = math.Sqrt(float64(m*(m+1)+3) - gm1*gm1 - hm2*hm2)
	}
	return f, g, h
}

// qbfsBasis fills q[m] = Q_m^bfs(x) and dq[m] = dQ_m^bfs/dx for m < len(q).
//
// The auxiliary Jacobi-type polynomials P_m (Forbes 2010 Eq. (2.6)) obey
// P_0 = 2, P_1 = 6 − 8x, P_{m+1} = (2 − 4x) P_m − P_{m−1}, and the orthogonal
// basis is generated by Eq. (2.7):
//
//	Q_m = [P_m − g_{m−1} Q_{m−1} − h_{m−2} Q_{m−2}] / f_m ,
//
// so that Q_0 = 1 and Q_1 = (13 − 16x)/√19.
func qbfsBasis(x float64, q, dq []float64) {
	n := len(q)
	if n == 0 {
		return
	}
	f, g, h := qbfsConstants(n)
	c := 2 - 4*x
	var pm2, pm1, dpm2, dpm1 float64 // P_{m−2}, P_{m−1} and their x-derivatives
	for m := 0; m < n; m++ {
		var p, dp float64
		switch m {
		case 0:
			p, dp = 2, 0
		case 1:
			p, dp = 6-8*x, -8
		default:
			p = c*pm1 - pm2
			dp = -4*pm1 + c*dpm1 - dpm2
		}
		qm, dqm := p, dp
		if m >= 1 {
			qm -= g[m-1] * q[m-1]
			dqm -= g[m-1] * dq[m-1]
		}
		if m >= 2 {
			qm -= h[m-2] * q[m-2]
			dqm -= h[m-2] * dq[m-2]
		}
		q[m] = qm / f[m]
		dq[m] = dqm / f[m]
		pm2, dpm2 = pm1, dpm1
		pm1, dpm1 = p, dp
	}
}

// qbfsSum returns S(x) = Σ_m a_m Q_m^bfs(x) and dS/dx at x = u².
func qbfsSum(a []float64, x float64) (sum, dsum float64) {
	n := len(a)
	if n == 0 {
		return 0, 0
	}
	q := make([]float64, n)
	dq := make([]float64, n)
	qbfsBasis(x, q, dq)
	for m, c := range a {
		sum += c * q[m]
		dsum += c * dq[m]
	}
	return sum, dsum
}
