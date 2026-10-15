package geo

import (
	"math"
	"math/cmplx"
)

// JonesMat is a 2×2 complex Jones matrix acting on the (Ex, Ey) components of a
// ray in its transverse frame (see TransverseFrame).
type JonesMat [2][2]complex128

// Mul returns the product a·b (b acts first).
func (a JonesMat) Mul(b JonesMat) JonesMat {
	var m JonesMat
	for i := 0; i < 2; i++ {
		for j := 0; j < 2; j++ {
			m[i][j] = a[i][0]*b[0][j] + a[i][1]*b[1][j]
		}
	}
	return m
}

// Apply returns a·v for a Jones vector v.
func (a JonesMat) Apply(v [2]complex128) [2]complex128 {
	return [2]complex128{
		a[0][0]*v[0] + a[0][1]*v[1],
		a[1][0]*v[0] + a[1][1]*v[1],
	}
}

// RetarderJones returns the Jones matrix of a linear retarder (wave plate) with
// retardance delta (rad) whose fast axis is at angle phi (rad) from ex. With the
// e^{-iωt} convention the slow component accumulates the extra phase delta:
// phase_slow − phase_fast = delta. In the fast/slow basis the matrix is
// diag(e^{-iδ/2}, e^{+iδ/2}) up to a global phase, rotated by phi:
//
//	J = U·diag(e^{-iδ/2}, e^{+iδ/2})·Uᵀ ,  U = [[cosφ, −sinφ], [sinφ, cosφ]] .
func RetarderJones(delta, phi float64) JonesMat {
	c, s := math.Cos(phi), math.Sin(phi)
	f := cmplx.Exp(complex(0, -delta/2)) // fast axis
	w := cmplx.Exp(complex(0, delta/2))  // slow axis
	return JonesMat{
		{complex(c*c, 0)*f + complex(s*s, 0)*w, complex(c*s, 0) * (f - w)},
		{complex(c*s, 0) * (f - w), complex(s*s, 0)*f + complex(c*c, 0)*w},
	}
}

// QuarterWaveJones returns a quarter-wave plate (δ = π/2) with fast axis at phi.
func QuarterWaveJones(phi float64) JonesMat { return RetarderJones(math.Pi/2, phi) }

// HalfWaveJones returns a half-wave plate (δ = π) with fast axis at phi.
func HalfWaveJones(phi float64) JonesMat { return RetarderJones(math.Pi, phi) }

// PolarizerJones returns an ideal linear polarizer whose transmission axis is at
// angle phi (rad) from ex: the projector onto (cosφ, sinφ).
func PolarizerJones(phi float64) JonesMat {
	c, s := math.Cos(phi), math.Sin(phi)
	return JonesMat{
		{complex(c*c, 0), complex(c*s, 0)},
		{complex(c*s, 0), complex(s*s, 0)},
	}
}

// JonesLinear returns the Jones vector of linear polarization at angle theta
// (rad) from ex, normalized to unit intensity.
func JonesLinear(theta float64) [2]complex128 {
	return [2]complex128{complex(math.Cos(theta), 0), complex(math.Sin(theta), 0)}
}

// JonesCircular returns the normalized circular Jones vector (1, ±i)/√2. With
// left true it is left-handed in the optics/polar.go naming (S3 > 0).
func JonesCircular(left bool) [2]complex128 {
	r := complex(1/math.Sqrt2, 0)
	if left {
		return [2]complex128{r, complex(0, 1/math.Sqrt2)}
	}
	return [2]complex128{r, complex(0, -1/math.Sqrt2)}
}
