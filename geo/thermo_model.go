package geo

import (
	"fmt"
	"math"
)

// Thermal model of a glass medium and of the geometry of a surface (G9).
//
// Index. For a glass at temperature change ΔT = T − 20 °C the absolute index
// follows the TIE-19 integral (see ThermalGlass). A lens in air needs the
// index relative to air, which is what the catalogue indices and Glass()
// provide; the air index itself depends on temperature (TIE-19 eq. 8–10,
// n_air − 1 ∝ 1/T), so
//
//	Δn_rel(λ, ΔT) = [ Δn_abs(λ, ΔT) − n_rel0 · (n_air(T) − n_air(20 °C)) ] / n_air(T),
//
// with n_rel0 the catalogue index and T = 20 °C + ΔT. The air term is not
// small: for N-BK7 at the e-line the relative coefficient is about 3.0e-6 /K
// while the absolute one is 1.6e-6 /K. Its 1/T nonlinearity (about 6 % over
// 20…40 °C) is needed to reproduce the datasheet relative values.
// The geometric paraxial engine works with air = 1 exactly, so the relative
// index is the one that makes the EFL of a lens in air correct.
//
// Geometry. An isotropic change of the element by the factor s = 1 + α ΔT
// scales every length of the surface (radius, sag, aperture) by s; an even
// asphere A_k r^(2k+4) therefore carries A_k·s^-(2k+3) (see ScaleProfile).

// airIndexAt returns the air index relative to vacuum at temperature tC (°C),
// 1 atm and wavelength lambdaNm, and its exact temperature derivative (1/K).
// It is TIE-19 eq. (8)–(9) with P = P0: n_air − 1 ∝ 1/(1 + 3.4785e-3 (T − 15)),
// i.e. the ideal-gas 1/T law. The 1/T nonlinearity is not negligible here:
// over 20…40 °C it changes the air contribution to the relative coefficient
// by about 6 % (≈0.1e-6 /K for n ≈ 1.7), which the datasheets' +20/+40 °C
// averages contain.
func airIndexAt(lambdaNm, tC float64) (n, dn float64) {
	l := lambdaNm / 1000 // µm
	l2 := l * l
	n15 := 1 + 1e-8*(6432.8+2949810*l2/(146*l2-1)+25540*l2/(41*l2-1))
	den := 1 + 3.4785e-3*(tC-15)
	n = 1 + (n15-1)/den
	dn = -3.4785e-3 * (n - 1) / den
	return n, dn
}

// deltaNAbs is the TIE-19 eq. (3) absolute index change at ΔT (K).
func (g ThermalGlass) deltaNAbs(lambdaNm, dT, n0 float64) float64 {
	l := lambdaNm / 1000
	den := l*l - g.LambdaTK*g.LambdaTK
	c := (n0*n0 - 1) / (2 * n0)
	return c * (g.D0*dT + g.D1*dT*dT + g.D2*dT*dT*dT + (g.E0*dT+g.E1*dT*dT)/den)
}

// dnAbsDT is the TIE-19 eq. (2) absolute coefficient at ΔT (K), in 1/K.
func (g ThermalGlass) dnAbsDT(lambdaNm, dT, n0 float64) float64 {
	l := lambdaNm / 1000
	den := l*l - g.LambdaTK*g.LambdaTK
	c := (n0*n0 - 1) / (2 * n0)
	return c * (g.D0 + 2*g.D1*dT + 3*g.D2*dT*dT + (g.E0+2*g.E1*dT)/den)
}

// DeltaNRel returns the change of the index relative to air, Δn_rel(λ, ΔT),
// for a glass whose catalogue index at 20 °C is nRel0. It is the exact
// relation n_rel(T) = n_abs(T)/n_air(T) with n_abs(T) = n_rel0·n_air(20 °C) +
// Δn_abs(ΔT):
//
//	Δn_rel = [ Δn_abs(ΔT) − n_rel0 · (n_air(T) − n_air(20 °C)) ] / n_air(T),
//
// T = 20 °C + ΔT. dT is in kelvin, lambdaNm in nm.
func (g ThermalGlass) DeltaNRel(lambdaNm, dT, nRel0 float64) float64 {
	nAir20, _ := airIndexAt(lambdaNm, 20)
	nAir, _ := airIndexAt(lambdaNm, 20+dT)
	return (g.deltaNAbs(lambdaNm, dT, nRel0) - nRel0*(nAir-nAir20)) / nAir
}

// DnRelDT returns the derivative d n_rel / dT (1/K) at temperature change dT,
// the exact derivative of DeltaNRel:
//
//	dn_rel/dT = [ dn_abs/dT − n_rel(T) · dn_air/dT(T) ] / n_air(T).
func (g ThermalGlass) DnRelDT(lambdaNm, dT, nRel0 float64) float64 {
	nAir, dnAir := airIndexAt(lambdaNm, 20+dT)
	nRel := nRel0 + g.DeltaNRel(lambdaNm, dT, nRel0)
	return (g.dnAbsDT(lambdaNm, dT, nRel0) - nRel*dnAir) / nAir
}

// DnAbsDT returns the absolute (vacuum-referenced) coefficient dn_abs/dT
// (1/K) of TIE-19 eq. (2). It is exposed for cross-checks against datasheets;
// the lens model uses DnRelDT.
func (g ThermalGlass) DnAbsDT(lambdaNm, dT, nRel0 float64) float64 {
	return g.dnAbsDT(lambdaNm, dT, nRel0)
}

// ThermalMedium is a catalogue glass at temperature change DT (K) from 20 °C:
// its index is the base (catalogue) index plus the relative index change of
// the glass. With DT == 0 the result is the base index bit for bit, so a
// system built with ThermalMedium at ΔT = 0 is identical to the original.
type ThermalMedium struct {
	Base  Medium
	Glass ThermalGlass
	DT    float64 // temperature change from 20 °C (K)
}

// N implements Medium.
func (m ThermalMedium) N(wavelengthNm float64) float64 {
	n0 := m.Base.N(wavelengthNm)
	if m.DT == 0 {
		return n0
	}
	return n0 + m.Glass.DeltaNRel(wavelengthNm, m.DT, n0)
}

// Name implements Medium.
func (m ThermalMedium) Name() string {
	if m.DT == 0 {
		return m.Base.Name()
	}
	return fmt.Sprintf("%s@%+.4gK", m.Base.Name(), m.DT)
}

// ScaleProfile returns the profile scaled isotropically by s > 0 about its
// vertex: z'(r) = s·z(r/s). Sphere, Conic and EvenAsphere are supported;
// their radii scale by s and the aspheric coefficients by s^-(2k+3) for the
// A_(2k+4) term. A plane is returned unchanged. s == 1 returns an equal
// profile, bit for bit.
func ScaleProfile(p Profile, s float64) (Profile, error) {
	if !(s > 0) || math.IsInf(s, 0) {
		return nil, fmt.Errorf("geo: ScaleProfile: scale factor %v must be positive and finite", s)
	}
	switch q := p.(type) {
	case Plane:
		return q, nil
	case Sphere:
		return Sphere{R: q.R * s}, nil
	case Conic:
		return Conic{R: q.R * s, K: q.K}, nil
	case EvenAsphere:
		out := EvenAsphere{R: q.R * s, K: q.K}
		if q.Coeffs != nil {
			out.Coeffs = make([]float64, len(q.Coeffs))
			for k, a := range q.Coeffs {
				if s == 1 {
					out.Coeffs[k] = a
				} else {
					out.Coeffs[k] = a * math.Pow(s, -float64(2*k+3))
				}
			}
		}
		return out, nil
	}
	return nil, fmt.Errorf("geo: ScaleProfile: unsupported profile %T", p)
}

// ScaleAperture returns the aperture scaled by s. A nil aperture stays nil.
func ScaleAperture(a Aperture, s float64) Aperture {
	switch q := a.(type) {
	case nil:
		return nil
	case CircleAperture:
		return CircleAperture(float64(q) * s)
	case RectangleAperture:
		return RectangleAperture{HalfWidth: q.HalfWidth * s, HalfHeight: q.HalfHeight * s}
	case AnnulusAperture:
		return AnnulusAperture{Inner: q.Inner * s, Outer: q.Outer * s}
	}
	return a
}

// ScaleSurface returns the surface with its profile and clear aperture scaled
// isotropically by s. The interaction type and the label are kept.
func ScaleSurface(sf Surface, s float64) (Surface, error) {
	p, err := ScaleProfile(sf.Profile, s)
	if err != nil {
		return Surface{}, err
	}
	sf.Profile = p
	sf.Aperture = ScaleAperture(sf.Aperture, s)
	return sf, nil
}

// ThinLensEFLDrift is the analytic temperature derivative dEFL/dT (mm/K) of a
// thin lens in air under the model above. With s = 1 + α ΔT the radii are
// R·s, so the curvatures are c_i/s and their derivative is −(α/s) c_i:
//
//	Φ = (n − 1)(c1 − c2),  dΦ/dT = (c1 − c2)(n' − (α/s)(n − 1)),
//	dEFL/dT = −Φ'/Φ² = EFL · ( α/s − n'/(n − 1) ).
//
// The thickness of a thin lens does not enter. nRel and dnRelDT are the
// relative index and its derivative at dT; cte is α and dT the temperature
// change from 20 °C (K).
func ThinLensEFLDrift(efl, nRel, dnRelDT, cte, dT float64) float64 {
	s := 1 + cte*dT
	return efl * (cte/s - dnRelDT/(nRel-1))
}
