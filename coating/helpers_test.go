package coating

import (
	"math"
	"math/cmplx"
	"math/rand"
	"testing"
)

// This file holds the helpers shared by the package tests.

// sqrtPrincipal returns the square root of q with Im >= 0, written out locally
// so the reference implementations below do not depend on package internals
// (the package uses the exp(-i·omega·t) convention, for which the forward wave
// exp(i·k_z·z) must decay inside a lossy medium).
func sqrtPrincipal(q complex128) complex128 {
	r := cmplx.Sqrt(q)
	if imag(r) < 0 || (imag(r) == 0 && real(r) < 0) {
		r = -r
	}
	return r
}

// refAdmittance returns the tilted admittance of a medium for the reference
// implementations (etaS = n·cos(theta), etaP = n/cos(theta)).
func refAdmittance(n, alpha complex128, pol Polarization) complex128 {
	c := sqrtPrincipal(n*n - alpha*alpha)
	if pol == PolP {
		return n * n / c
	}
	return c
}

// airySingleLayer evaluates the closed-form Airy summation of a single layer
// between two semi-infinite media:
//
//	r = (r1 + r2·e^{+2i·delta}) / (1 + r1·r2·e^{+2i·delta})
//	t = t1·t2·e^{+i·delta} / (1 + r1·r2·e^{+2i·delta})
//
// with delta = k0·d·n1·cos(theta1) and the Fresnel coefficients of the two
// interfaces. It is an independent formulation of the same physics as the
// characteristic-matrix product.
func airySingleLayer(lambdaNm float64, n0, n1, ns complex128, dNm, thetaDeg float64, pol Polarization) (r, tt complex128) {
	k0 := 2 * math.Pi / lambdaNm
	alpha := n0 * complex(math.Sin(thetaDeg*math.Pi/180), 0)
	e0 := refAdmittance(n0, alpha, pol)
	e1 := refAdmittance(n1, alpha, pol)
	es := refAdmittance(ns, alpha, pol)
	delta := complex(k0*dNm, 0) * sqrtPrincipal(n1*n1-alpha*alpha)

	r1 := (e0 - e1) / (e0 + e1)
	r2 := (e1 - es) / (e1 + es)
	t1 := 2 * e0 / (e0 + e1)
	t2 := 2 * e1 / (e1 + es)

	phase := cmplx.Exp(complex(0, 2) * delta)
	den := 1 + r1*r2*phase
	return (r1 + r2*phase) / den, t1 * t2 * cmplx.Exp(complex(0, 1)*delta) / den
}

// recursionReflectance evaluates the complex amplitude reflectance of a stack
// through the classic admittance recurrence, here written for the
// exp(-i·omega·t) convention of the package:
//
//	Y_front = eta_j·(-i·eta_j·sin(delta_j) + Y_back·cos(delta_j)) /
//	                (eta_j·cos(delta_j) - i·Y_back·sin(delta_j)),
//
// starting at the substrate (Y = eta_substrate) and working towards the
// incidence medium. It is a second, independent implementation used to
// cross-check the matrix product (and it exercises complex indices, where the
// closed-form Airy expression does not apply).
func recursionReflectance(lambdaNm float64, n0 complex128, layers []Layer, ns complex128, thetaDeg float64, pol Polarization) complex128 {
	k0 := 2 * math.Pi / lambdaNm
	alpha := n0 * complex(math.Sin(thetaDeg*math.Pi/180), 0)
	eta0 := refAdmittance(n0, alpha, pol)
	Y := refAdmittance(ns, alpha, pol)
	for i := len(layers) - 1; i >= 0; i-- {
		n := layers[i].Material.IndexNm(lambdaNm)
		e := refAdmittance(n, alpha, pol)
		delta := complex(k0*layers[i].ThicknessNm, 0) * sqrtPrincipal(n*n-alpha*alpha)
		cd, sd := cmplx.Cos(delta), cmplx.Sin(delta)
		num := -complex(0, 1)*e*sd + Y*cd
		den := e*cd - complex(0, 1)*Y*sd
		Y = e * num / den
	}
	return (eta0 - Y) / (eta0 + Y)
}

// closeComplex fails the test when two complex numbers differ by more than tol
// (relative to the magnitude of the reference value).
func closeComplex(t *testing.T, name string, got, want complex128, tol float64) {
	t.Helper()
	if d := cmplx.Abs(got - want); d > tol*(1+cmplx.Abs(want)) {
		t.Errorf("%s: got %v, want %v (|diff| = %g)", name, got, want, d)
	}
}

// closeFloat fails the test when two real numbers differ by more than tol
// (relative to the magnitude of the reference value).
func closeFloat(t *testing.T, name string, got, want, tol float64) {
	t.Helper()
	if d := math.Abs(got - want); d > tol*(1+math.Abs(want)) {
		t.Errorf("%s: got %.12g, want %.12g (|diff| = %g)", name, got, want, d)
	}
}

// randomRealMaterial builds a lossless material with a random index in [nMin, nMax].
func randomRealMaterial(rng *rand.Rand, name string, nMin, nMax float64) Material {
	return NewMaterial(name, ConstIndex{N: nMin + rng.Float64()*(nMax-nMin)})
}

// randomStack builds a random stack of n layers with real or complex indices.
func randomStack(rng *rand.Rand, n int, lossy bool, nMin, nMax, dMin, dMax float64) Stack {
	layers := make([]Layer, n)
	for i := range layers {
		nv := nMin + rng.Float64()*(nMax-nMin)
		kv := 0.0
		if lossy {
			kv = rng.Float64() * 0.05
		}
		layers[i] = Layer{
			Material:    NewMaterial("random", ConstIndex{N: nv, K: kv}),
			ThicknessNm: dMin + rng.Float64()*(dMax-dMin),
		}
	}
	return NewStack(Air, NewMaterial("random substrate", ConstIndex{N: 1.4 + rng.Float64() * 0.3}), layers...)
}
