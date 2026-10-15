package geo

import (
	"fmt"
	"math"
	"testing"
)

// ---------------------------------------------------------------------------
// Independent references for the Q-polynomial tests.
// ---------------------------------------------------------------------------

// aspBinom returns the binomial coefficient C(n, k) as a float64 (n, k ≥ 0).
func aspBinom(n, k int) float64 {
	if k < 0 || k > n {
		return 0
	}
	r := 1.0
	for i := 1; i <= k; i++ {
		r = r * float64(n-k+i) / float64(i)
	}
	return r
}

// aspJacobi04 evaluates P_n^(0,4)(t) from the explicit binomial sum of the
// Jacobi polynomial (Abramowitz & Stegun ch. 22, the z-form of the explicit
// sum):
//
//	P_n^(α,β)(t) = Σ_s C(n+α, n−s) C(n+β, s) ((t−1)/2)^s ((t+1)/2)^(n−s),
//
// with α = 0, β = 4. It is independent of the three-term recurrence used by
// qconBasis.
func aspJacobi04(n int, t float64) float64 {
	sum := 0.0
	for s := 0; s <= n; s++ {
		sum += aspBinom(n, s) * aspBinom(n+4, s) * math.Pow((t-1)/2, float64(s)) * math.Pow((t+1)/2, float64(n-s))
	}
	return sum
}

// aspGaussLegendre returns the n-point Gauss–Legendre nodes and weights on
// [−1, 1] (Newton iteration on P_n; exact for polynomials of degree ≤ 2n−1).
func aspGaussLegendre(n int) (x, w []float64) {
	x = make([]float64, n)
	w = make([]float64, n)
	for i := 0; i < n; i++ {
		xi := math.Cos(math.Pi * (float64(i) + 0.75) / (float64(n) + 0.5))
		var dp float64
		for it := 0; it < 100; it++ {
			p0, p1 := 1.0, xi // P_0, P_1
			for k := 2; k <= n; k++ {
				fk := float64(k)
				p0, p1 = p1, ((2*fk-1)*xi*p1-(fk-1)*p0)/fk
			}
			// p1 = P_n(xi), p0 = P_{n−1}(xi)
			dp = float64(n) * (xi*p1 - p0) / (xi*xi - 1)
			dx := p1 / dp
			xi -= dx
			if math.Abs(dx) < 1e-15 {
				break
			}
		}
		x[i] = xi
		w[i] = 2 / ((1 - xi*xi) * dp * dp)
	}
	return x, w
}

// aspQPolyP returns the Forbes auxiliary polynomials P_0..P_{n−1} at x by the
// recurrence of Forbes (2010) Eq. (2.6): P_0 = 2, P_1 = 6 − 8x,
// P_{m+1} = (2 − 4x) P_m − P_{m−1}. It is written here, independently of
// qbfsBasis, so that Eq. (2.5) can be checked against it.
func aspQPolyP(n int, x float64) []float64 {
	p := make([]float64, n)
	if n > 0 {
		p[0] = 2
	}
	if n > 1 {
		p[1] = 6 - 8*x
	}
	for m := 2; m < n; m++ {
		p[m] = (2-4*x)*p[m-1] - p[m-2]
	}
	return p
}

// ---------------------------------------------------------------------------
// Qcon basis
// ---------------------------------------------------------------------------

// TestQconBasisAgainstExplicitJacobi cross-validates the Qcon basis Q_m^con(x) =
// P_m^(0,4)(2x−1) (Forbes 2010 Eq. (5.1)–(5.2)): the three-term recurrence and
// its derivative against the explicit binomial sum and its finite difference,
// plus P_m(1) = 1 and the explicit low orders P_1 = 3t−2, P_2 = 7t²−7t+1.
func TestQconBasisAgainstExplicitJacobi(t *testing.T) {
	xs := []float64{0, 0.05, 0.13, 0.3, 0.5, 0.77, 0.95, 1, 1.2}
	const N = 12
	z := make([]float64, N)
	dz := make([]float64, N)
	for _, x := range xs {
		qconBasis(x, z, dz)
		tt := 2*x - 1
		for m := 0; m < N; m++ {
			want := aspJacobi04(m, tt)
			if d := math.Abs(z[m] - want); d > 1e-12*math.Max(1, math.Abs(want)) {
				t.Errorf("Q_%d^con(%g): recurrence %.15g, explicit sum %.15g", m, x, z[m], want)
			}
			const h = 1e-6
			num := (aspJacobi04(m, 2*(x+h)-1) - aspJacobi04(m, 2*(x-h)-1)) / (2 * h)
			if d := math.Abs(dz[m] - num); d > 1e-7*math.Max(1, math.Abs(num)) {
				t.Errorf("dQ_%d^con/dx(%g): analytic %.12g, central difference %.12g", m, x, dz[m], num)
			}
		}
	}
	// Normalisation at the edge of the aperture: P_m^(0,4)(1) = 1 (x = 1).
	qconBasis(1, z, dz)
	for m := 0; m < N; m++ {
		checkClose(t, fmt.Sprintf("Q_%d^con(1)", m), z[m], 1, 1e-13)
	}
	// Explicit low orders in x (t = 2x−1): Q_1 = 6x − 5, Q_2 = 7(2x−1)² − 7(2x−1) + 1.
	qconBasis(0.3, z, dz)
	checkClose(t, "Q_1^con(0.3)", z[1], 6*0.3-5, 1e-14)
	tt := 2*0.3 - 1
	checkClose(t, "Q_2^con(0.3)", z[2], 7*tt*tt-7*tt+1, 1e-13)
}

// TestQconJacobiOrthogonality checks the orthogonality of the Qcon basis against
// its Jacobi weight, ∫_{−1}^{1} P_m P_n (1+t)^4 dt = δ_mn · 32/(2n+5), by exact
// Gauss–Legendre quadrature (the integrand is a polynomial of degree ≤ 2N−1).
func TestQconJacobiOrthogonality(t *testing.T) {
	const M = 9
	nodes, weights := aspGaussLegendre(24)
	// vals[k][m] = Q_m^con at node k, through x = (t+1)/2.
	vals := make([][]float64, len(nodes))
	dummy := make([]float64, M)
	for k, tt := range nodes {
		vals[k] = make([]float64, M)
		qconBasis((tt+1)/2, vals[k], dummy)
	}
	for m := 0; m < M; m++ {
		for n := 0; n < M; n++ {
			s := 0.0
			for k, tt := range nodes {
				s += weights[k] * vals[k][m] * vals[k][n] * math.Pow(1+tt, 4)
			}
			want := 0.0
			if m == n {
				want = 32 / float64(2*n+5)
			}
			if d := math.Abs(s - want); d > 1e-12*math.Max(1, math.Abs(want)) {
				t.Errorf("∫P_%d P_%d (1+t)^4 = %.15g, want %.15g", m, n, s, want)
			}
		}
	}
}

// ---------------------------------------------------------------------------
// Qbfs basis
// ---------------------------------------------------------------------------

// TestQbfsForbesIdentities checks the Qbfs basis against the identities printed
// in Forbes (2010, OE 18(19) 19700):
//   - Q_0 = 1 and Q_1 = 19^(−1/2)(13 − 16x) (Eq. 2.7 text);
//   - Eq. (2.5): P_m = f_m Q_m + g_{m−1} Q_{m−1} + h_{m−2} Q_{m−2}, with P from the
//     independent recurrence aspQPolyP;
//   - the constant matrix of Eq. (A.17) (f_2, f_3, g_1, g_2, h_0, h_1).
func TestQbfsForbesIdentities(t *testing.T) {
	const M = 14
	f, g, h := qbfsConstants(M)
	// (A.17): printed entries of the change-of-basis matrix.
	checkRel(t, "A.17 f_2", f[2], 4*math.Sqrt(10.0/19), 1e-14)
	checkRel(t, "A.17 f_3", f[3], 0.5*math.Sqrt(509.0/10), 1e-14)
	checkRel(t, "A.17 g_1", g[1], -5/(2*math.Sqrt(19)), 1e-14)
	checkRel(t, "A.17 g_2", g[2], -17/(2*math.Sqrt(190)), 1e-14)
	checkRel(t, "A.17 h_0", h[0], -0.5, 1e-14)
	checkRel(t, "A.17 h_1", h[1], -6/math.Sqrt(19), 1e-14)
	checkRel(t, "A.17 f_0", f[0], 2, 0)
	checkRel(t, "A.17 f_1", f[1], math.Sqrt(19)/2, 1e-14)

	q := make([]float64, M)
	dq := make([]float64, M)
	for _, x := range []float64{0, 0.11, 0.3, 0.62, 0.9, 1.0} {
		qbfsBasis(x, q, dq)
		checkRel(t, fmt.Sprintf("Q_0(%g)", x), q[0], 1, 1e-15)
		checkRel(t, fmt.Sprintf("Q_1(%g)", x), q[1], (13-16*x)/math.Sqrt(19), 1e-14)
		p := aspQPolyP(M, x)
		for m := 0; m < M; m++ {
			rhs := f[m] * q[m]
			if m >= 1 {
				rhs += g[m-1] * q[m-1]
			}
			if m >= 2 {
				rhs += h[m-2] * q[m-2]
			}
			if d := math.Abs(rhs - p[m]); d > 1e-11*math.Max(1, math.Abs(p[m])) {
				t.Errorf("Eq.(2.5) m=%d x=%g: f Q + g Q + h Q = %.15g, P_m = %.15g", m, x, rhs, p[m])
			}
		}
	}
}

// aspSbfs returns the normal-slope basis S_m(u) = d/du[u²(1−u²)Q_m(u²)]
// (Forbes 2010, OE 18(19), Eq. (2.3)) at u, written as 2u·[(1−2x)Q_m + x(1−x)Q_m'] with
// x = u² (differentiate the product u²(1−u²) and the chain rule dx/du = 2u).
func aspSbfs(u float64, m int, q, dq []float64) float64 {
	x := u * u
	return 2 * u * ((1-2*x)*q[m] + x*(1-x)*dq[m])
}

// TestQbfsSlopeOrthonormal is the defining property of the Qbfs basis (Forbes
// Forbes 2010, OE 18(19), Eq. (2.2)):
//
//	(2/π) ∫_0^1 S_m(u) S_n(u) / sqrt(1 − u²) du = δ_mn .
//
// The integrand is even in u and polynomial, so the Gauss–Chebyshev rule with
// N nodes, (1/N) Σ_k f(cos θ_k), θ_k = (2k−1)π/(2N), is exact for it.
func TestQbfsSlopeOrthonormal(t *testing.T) {
	const M = 11
	const N = 64
	q := make([]float64, M)
	dq := make([]float64, M)
	// S_m at each Chebyshev node.
	s := make([][]float64, M)
	for m := range s {
		s[m] = make([]float64, N)
	}
	for k := 0; k < N; k++ {
		u := math.Cos(float64(2*k+1) * math.Pi / float64(2*N))
		qbfsBasis(u*u, q, dq)
		for m := 0; m < M; m++ {
			s[m][k] = aspSbfs(u, m, q, dq)
		}
	}
	for m := 0; m < M; m++ {
		for n := 0; n < M; n++ {
			sum := 0.0
			for k := 0; k < N; k++ {
				sum += s[m][k] * s[n][k]
			}
			got := sum / N
			want := 0.0
			if m == n {
				want = 1
			}
			if d := math.Abs(got - want); d > 1e-11 {
				t.Errorf("(2/π)∫S_%d S_%d/√(1−u²) = %.15g, want δ = %g (|diff| %.3g)", m, n, got, want, d)
			}
		}
	}
}

// TestQbfsSlopeFormulaVsDerivative checks Forbes' S_m(u) = d/du[u²(1−u²)Q_m(u²)]
// (the closed form aspSbfs) against a central difference of the product.
func TestQbfsSlopeFormulaVsDerivative(t *testing.T) {
	const M = 8
	q := make([]float64, M)
	dq := make([]float64, M)
	qq := make([]float64, M)
	dqq := make([]float64, M)
	const h = 1e-6
	for _, u := range []float64{0.05, 0.3, 0.55, 0.8, 0.95} {
		for m := 0; m < M; m++ {
			qbfsBasis(u*u, q, dq)
			F := func(v float64) float64 {
				qbfsBasis(v*v, qq, dqq)
				return v * v * (1 - v*v) * qq[m]
			}
			num := (F(u+h) - F(u-h)) / (2 * h)
			got := aspSbfs(u, m, q, dq)
			if d := math.Abs(got - num); d > 1e-7*math.Max(1, math.Abs(num)) {
				t.Errorf("S_%d(%g): formula %.12g, central difference %.12g", m, u, got, num)
			}
		}
	}
}

// ---------------------------------------------------------------------------
// QAsphere: sag, slope, curvature, intersection and degeneracies.
// ---------------------------------------------------------------------------

// aspQconProfile and aspQbfsProfile return the representative Q-type asphere
// used throughout the tests.
func aspQconProfile() QAsphere {
	return QAsphere{Form: QCon, R: 30, K: -0.4, NormR: 9,
		Coeffs: []float64{1e-4, -2e-5, 3e-6, 1e-7, -2e-8, 5e-9}}
}

func aspQbfsProfile() QAsphere {
	return QAsphere{Form: QBfs, R: 45, K: 0, NormR: 9,
		Coeffs: []float64{2e-4, -3e-5, 1e-5, -4e-6, 7e-7, 2e-7, -5e-8}}
}

// TestQAsphereSlopeCentralDifference is acceptance (a) for both Forbes forms:
// the analytic slope agrees with the central difference of the sag to 1e-7
// (mixed relative/absolute) at every sample radius.
func TestQAsphereSlopeCentralDifference(t *testing.T) {
	checkSlopeCD(t, "qcon", aspQconProfile(), aspRadii)
	checkSlopeCD(t, "qbfs", aspQbfsProfile(), aspRadii)
	// Qcon with only the first and a high-order member, and a sag that has a
	// pure r⁴ shape from s_0 alone.
	checkSlopeCD(t, "qcon sparse", QAsphere{Form: QCon, R: -22, K: 0.5, NormR: 8,
		Coeffs: []float64{0, 0, 4e-6, 0, -1e-9}}, aspRadii)
}

// TestQAsphereLeadingTermClosedForm checks the single-member closed forms.
// Q_0^con = 1 gives z = conic + s_0 (r/ρ_max)⁴; Q_0^bfs = 1 gives
// z = sphere + a_0 u²(1−u²)/sqrt(1 − c²r²). Derivatives are also closed form.
func TestQAsphereLeadingTermClosedForm(t *testing.T) {
	const R, K, normR, s0 = 31.25, -1.7, 10.0, 2.5e-6
	qc := QAsphere{Form: QCon, R: R, K: K, NormR: normR, Coeffs: []float64{s0}}
	base := Conic{R: R, K: K}
	for _, r := range aspRadii {
		u := r / normR
		checkRel(t, "qcon s0 sag", qc.Sag(r), base.Sag(r)+s0*math.Pow(u, 4), 1e-14)
		checkRel(t, "qcon s0 slope", qc.Slope(r), base.Slope(r)+4*s0*math.Pow(r, 3)/math.Pow(normR, 4), 1e-13)
	}
	const Rb, a0 = 60.0, 3e-5
	qb := QAsphere{Form: QBfs, R: Rb, NormR: normR, Coeffs: []float64{a0}}
	c := 1 / Rb
	for _, r := range aspRadii {
		u := r / normR
		psi := math.Sqrt(1 - c*c*r*r)
		want := Sphere{R: Rb}.Sag(r) + a0*u*u*(1-u*u)/psi
		checkRel(t, "qbfs a0 sag", qb.Sag(r), want, 1e-13)
	}
}

// TestQAsphereZeroCoefficientsIsConic is acceptance (d): with all Q coefficients
// zero, the Qcon asphere is the conic base Conic{R, K} to 1e-13 in sag and slope,
// and the Qbfs asphere (K = 0) is Sphere{R}. Intersections agree to 3e-12 mm;
// the bound is set by the existing Conic Newton tolerance (see the EvenAsphere
// note in prof_poly_test.go).
func TestQAsphereZeroCoefficientsIsConic(t *testing.T) {
	const R, K = 31.25, -1.7
	conic := Conic{R: R, K: K}
	for _, coeffs := range [][]float64{nil, {0}, {0, 0, 0, 0}} {
		q := QAsphere{Form: QCon, R: R, K: K, NormR: 10, Coeffs: coeffs}
		for _, r := range aspRadii {
			if d := math.Abs(q.Sag(r) - conic.Sag(r)); !(d <= 1e-13) {
				t.Errorf("qcon zero coeffs %v: sag at r=%g differs by %.3g", coeffs, r, d)
			}
			if d := math.Abs(q.Slope(r) - conic.Slope(r)); !(d <= 1e-13) {
				t.Errorf("qcon zero coeffs %v: slope at r=%g differs by %.3g", coeffs, r, d)
			}
		}
		if q.Curvature() != conic.Curvature() {
			t.Errorf("qcon zero coeffs: curvature %g, conic %g", q.Curvature(), conic.Curvature())
		}
		maxDiff := 0.0
		for _, rv := range aspRandomRays(300, 6, 77) {
			tq, okq := q.Intersect(rv.o, rv.d, 0, math.Inf(1))
			tc, okc := conic.Intersect(rv.o, rv.d, 0, math.Inf(1))
			if okq != okc {
				t.Errorf("qcon zero coeffs: hit flags differ (qcon %v, conic %v)", okq, okc)
				continue
			}
			if okq {
				d := rv.o.Add(rv.d.Mul(tq)).Sub(rv.o.Add(rv.d.Mul(tc))).Norm()
				if d > maxDiff {
					maxDiff = d
				}
			}
		}
		if maxDiff > 3e-12 {
			t.Errorf("qcon zero coeffs: hit points differ by %.3g mm", maxDiff)
		}
	}
	// Qbfs with K = 0 and zero coefficients is the sphere.
	const Rs = 45.0
	qb := QAsphere{Form: QBfs, R: Rs, NormR: 9, Coeffs: []float64{0, 0, 0}}
	sph := Sphere{R: Rs}
	for _, r := range aspRadii {
		checkClose(t, "qbfs zero coeffs sag", qb.Sag(r), sph.Sag(r), 1e-13)
		checkClose(t, "qbfs zero coeffs slope", qb.Slope(r), sph.Slope(r), 1e-13)
	}
}

// TestQAsphereIntersectResidual is acceptance (b) for both forms: the ray
// intersection converges to a sag residual ≤ 1e-12 mm.
func TestQAsphereIntersectResidual(t *testing.T) {
	rays := aspRandomRays(4000, 7, 4242)
	n1 := checkIntersectResidual(t, "qcon intersect", aspQconProfile(), rays)
	n2 := checkIntersectResidual(t, "qbfs intersect", aspQbfsProfile(), rays)
	if n1 < len(rays)/4 || n2 < len(rays)/4 {
		t.Errorf("too few hits: qcon %d, qbfs %d of %d", n1, n2, len(rays))
	}
}

// TestQAsphereCurvature checks the paraxial vertex curvature. For Qcon the
// departure has no r² term (curvature 1/R). For Qbfs the vertex curvature is
// c + 2 S(0)/ρ_max² with S(0) = Σ a_m Q_m(0); with Q_0(0) = 1 and
// Q_1(0) = 13/√19 this gives a closed form for two members. The numerical
// slope/r at r → 0 must agree with the curvature.
func TestQAsphereCurvature(t *testing.T) {
	qc := aspQconProfile()
	checkClose(t, "qcon curvature", qc.Curvature(), 1/qc.R, 1e-16)
	// Numerical vertex curvature: z'(r)/r → c_eff.
	for _, p := range []QAsphere{qc, aspQbfsProfile()} {
		const r = 1e-4
		num := p.Slope(r) / r
		if d := math.Abs(num - p.Curvature()); d > 1e-6*math.Max(1, math.Abs(p.Curvature())) {
			t.Errorf("%v: slope/r at r=%g is %.12g, Curvature() %.12g", p.Form, r, num, p.Curvature())
		}
	}
	// Closed form with two Qbfs members.
	const R, normR, a0, a1 = 45.0, 9.0, 2e-4, -3e-5
	q := QAsphere{Form: QBfs, R: R, NormR: normR, Coeffs: []float64{a0, a1}}
	want := 1/R + 2*(a0+a1*13/math.Sqrt(19))/(normR*normR)
	checkRel(t, "qbfs two-member curvature", q.Curvature(), want, 1e-13)
}

// TestQAsphereValidate checks the structural error paths of QAsphere.
func TestQAsphereValidate(t *testing.T) {
	ok := []QAsphere{
		aspQconProfile(), aspQbfsProfile(),
		{Form: QCon, R: 20, K: -1}, // no coefficients: NormR not needed
		{Form: QBfs, R: 20, Coeffs: []float64{0, 0}},
	}
	for i, q := range ok {
		if err := q.Validate(); err != nil {
			t.Errorf("valid q asphere %d rejected: %v", i, err)
		}
	}
	bad := []QAsphere{
		{Form: QForm(7), R: 20, NormR: 5, Coeffs: []float64{1e-4}},
		{Form: QBfs, R: 20, K: -0.5, NormR: 5, Coeffs: []float64{1e-4}},
		{Form: QCon, R: math.NaN(), NormR: 5, Coeffs: []float64{1e-4}},
		{Form: QCon, R: 20, NormR: 5, Coeffs: []float64{math.Inf(1)}},
		{Form: QCon, R: 20, NormR: 0, Coeffs: []float64{1e-4}},
		{Form: QBfs, R: 20, NormR: -3, Coeffs: []float64{1e-4}},
	}
	for i, q := range bad {
		if err := q.Validate(); err == nil {
			t.Errorf("bad q asphere %d accepted: %v", i, q)
		}
	}
	if s := QCon.String() + " " + QBfs.String(); s != "qcon qbfs" {
		t.Errorf("form names: %q", s)
	}
}

// TestQAsphereParaxialDerivative is acceptance (e) for the Q-type asphere. The
// first-order EFL depends on the Q coefficients through the vertex curvature:
// zero for Qcon (no r² term) and 2·Q_m(0)/ρ_max² for the Qbfs member a_m. The
// analytic chain-rule derivative (exact ray-transfer matrices, see
// paraxialEFLdCurv) agrees with the central difference of the traced EFL to 1e-6
// relative. For Qbfs the factor is non-trivial for every m.
func TestQAsphereParaxialDerivative(t *testing.T) {
	const lambda = 550.0
	const R, normR = 28.0, 8.0
	// Qbfs: analytic factor 2 Q_m(0)/ρ_max² for member m.
	qz := make([]float64, 6)
	dqz := make([]float64, 6)
	qbfsBasis(0, qz, dqz)
	a0 := []float64{1e-4, 1e-5, -2e-5, 3e-6, -1e-6}
	for m := 0; m < len(a0); m++ {
		m := m
		profile := func(a float64) Profile {
			c := make([]float64, m+1)
			c[m] = a
			return QAsphere{Form: QBfs, R: R, NormR: normR, Coeffs: c}
		}
		factor := 2 * qz[m] / (normR * normR)
		checkParaxialDerivative(t, fmt.Sprintf("qbfs a%d", m), lambda, a0[m], factor, profile)
	}
	// Qcon: every member leaves the paraxial model unchanged.
	for m := 0; m < 4; m++ {
		m := m
		profile := func(a float64) Profile {
			c := make([]float64, m+1)
			c[m] = a
			return QAsphere{Form: QCon, R: R, K: -0.2, NormR: normR, Coeffs: c}
		}
		checkParaxialDerivative(t, fmt.Sprintf("qcon s%d", m), lambda, 1e-5, 0, profile)
	}
}
