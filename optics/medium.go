package optics

import (
	"fmt"
	"math"
)

// IndexFunc returns the complex refractive index n(x, y, z) at a point. The
// real part is the phase index and the imaginary part is the absorption
// (Im > 0) or gain (Im < 0) coefficient.
type IndexFunc func(x, y, z float64) complex128

// UniformIndex returns an IndexFunc for a constant (possibly complex) index.
func UniformIndex(n complex128) IndexFunc {
	return func(x, y, z float64) complex128 { return n }
}

// PropagateSplitStep advances the field by distance z through an inhomogeneous
// medium with refractive index n(x,y,z), using symmetric split-step (Strang)
// propagation: each of 'steps' slices applies [ASM dz/2, medium phase dz,
// ASM dz/2]. Because n may vary in x, y and z, stratified n(z), gradient
// n(x,y) and absorption/gain media are supported. This is a paraxial BPM
// approximation: the medium contributes the phase k0*(n-1)*dz in real space
// while the ASM kernel keeps the vacuum dispersion sqrt(1-(lambda0 f)^2).
func PropagateSplitStep(f *Field, z float64, n IndexFunc, steps int, ctx *Context) error {
	if z == 0 {
		return nil
	}
	if steps <= 0 {
		steps = 1
	}
	if ctx == nil {
		ctx = &Context{Evanescent: "decay"}
	}
	dz := z / float64(steps)
	k := 2 * math.Pi / ctx.Wavelength
	for s := 0; s < steps; s++ {
		zc := z*float64(s)/float64(steps) + dz/2
		if err := Propagate(f, dz/2, MethodASM, ctx); err != nil {
			return err
		}
		applyMediumPhase(f, n, zc, k, dz)
		if err := Propagate(f, dz/2, MethodASM, ctx); err != nil {
			return err
		}
	}
	return nil
}

// applyMediumPhase multiplies the field by exp(i k0 (n-1) dz) at every pixel
// for the slice centered at zc. Complex n: real part -> phase, imaginary part
// -> absorption (Im>0) or gain (Im<0).
func applyMediumPhase(f *Field, n IndexFunc, zc, k, dz float64) {
	nn := f.N
	for j := 0; j < nn; j++ {
		y := f.Y(j)
		for i := 0; i < nn; i++ {
			nv := n(f.X(i), y, zc)
			ph := k * (real(nv) - 1) * dz
			t := complex(math.Exp(-k*imag(nv)*dz), 0) * cexpI(ph)
			idx := j*nn + i
			f.Ex[idx] *= t
			if f.Polarized {
				f.Ey[idx] *= t
			}
		}
	}
}

// ---------------------------------------------------------------------------
// Reference-index split-step (reference-index BPM)
//
// PropagateSplitStep above removes the carrier with reference index 1: its
// diffraction half-steps use the *vacuum* transfer function
// exp(i z sqrt(k0^2 - kt^2)) and the medium contributes the whole phase
// k0 (n-1) dz. That is the classic "phase screen + vacuum diffraction" BPM.
// Its limit is well known and easy to see: a homogeneous layer of index n
// whose n does not vary across the beam contributes only a constant phase, so
// a slab does not refract — a plane-parallel plate leaves the focus exactly
// where it was, while a real plate moves it by d(1-1/n).
//
// PropagateSplitStepRef removes the carrier with an arbitrary reference index
// nref instead. Inside the layer the diffraction half-steps then use the exact
// homogeneous-medium transfer function of index nref,
//
//	H(fx,fy,z) = exp[ i z sqrt((nref k0)^2 - (2 pi fx)^2 - (2 pi fy)^2) ]
//	             (transverse wavevector below nref k0: propagating)
//	           = exp[ -|z| sqrt((2 pi f)^2 - (nref k0)^2) ]
//	             (beyond nref k0: evanescent inside the medium)
//
// and the screen contributes only the *difference* to that reference,
// exp(i k0 (n - nref) dz) * exp(-k0 Im(n) dz). For nref = 1 this is exactly
// the existing vacuum-reference split-step (tested to 1e-12); for
// nref = the layer's index it makes a homogeneous layer exact in the
// angular-spectrum sense — a plane-parallel plate then shows its true
// refraction: reduced optical depth d/nref inside the plate, the exact axial
// phase k0 n d, and the analytic focal shift d(1-1/nref) (see
// docs/PHYSICS.md §5.3 and TestCustomMediumPlateFocusShift).
//
// The operator per slice is the same symmetric (Strang) composition
// [ASM_nref dz/2, screen dz, ASM_nref dz/2]; the split-step error terms are
// those of the reference-index BPM, i.e. O(dz^3) local / O(dz^2) global in the
// commutator of the diffraction operator and the transverse index variation.
// The reference index must therefore be close to the local index for the
// expansion to be accurate: |n - nref| << nref. See docs/PHYSICS.md §5.3.
// ---------------------------------------------------------------------------

// PropagateSplitStepRef advances the field by distance z through an
// inhomogeneous medium with index n(x,y,z) using a symmetric split-step with a
// *reference index* nref for the diffraction kernel (see the block comment
// above). nref must be > 0; the medium's imaginary part is absorption (Im > 0)
// or gain (Im < 0) applied without reference, i.e. exp(-k0 Im(n) dz).
//
// It is the generic (index-function) form of the reference-index split-step;
// a compiled medium evaluates its own expressions inline instead (see
// propagateCustomLayer), and the two are cross-checked to 1e-12 in the tests.
func PropagateSplitStepRef(f *Field, z float64, n IndexFunc, nref float64, steps int, ctx *Context) error {
	if z == 0 {
		return nil
	}
	if steps <= 0 {
		steps = 1
	}
	if !(nref > 0) || math.IsInf(nref, 0) || math.IsNaN(nref) {
		return fmt.Errorf("split-step: reference index must be > 0 (got %g)", nref)
	}
	if ctx == nil {
		ctx = &Context{Evanescent: "decay"}
	}
	dz := z / float64(steps)
	k := 2 * math.Pi / ctx.Wavelength
	zero, reg, _, _ := ctx.evanescentPolicy(dz / 2)
	for s := 0; s < steps; s++ {
		zc := z*float64(s)/float64(steps) + dz/2
		pIn := f.Power()
		propagateMediumHalf(f, dz/2, nref, ctx)
		applyMediumPhaseRef(f, n, zc, k, dz, nref)
		propagateMediumHalf(f, dz/2, nref, ctx)
		asmWarn(f, dz/2, ctx, pIn, f.Power(), zero, reg)
	}
	return nil
}

// propagateMediumHalf applies half an angular-spectrum step inside a medium of
// reference index nref. The evanescent policy (Evanescent / EvanescentLimit /
// BackwardRegularize / TikhonovAlpha) is the same one the plain ASM step uses.
func propagateMediumHalf(f *Field, z, nref float64, ctx *Context) {
	if z == 0 {
		return
	}
	zero, reg, alpha, limit := ctx.evanescentPolicy(z)
	asmWrapWarn(f, z, ctx, float64(f.N)*f.DX)
	ctx.transfer(f, z, mediumTF(ctx.Wavelength, z, nref, zero, reg, alpha, limit))
}

// mediumTF returns the angular-spectrum transfer function of a homogeneous
// medium of index nref: H = exp(i z sqrt((nref k0)^2 - (2 pi f)^2)). The
// propagating branch is evaluated in the factored form
// sqrt((1-rho)(1+rho)) with rho = lambda f / nref — the same cancellation-safe
// arrangement asmTF uses — and wavevectors beyond nref k0 are evanescent inside
// the medium (decay forward, zero or Tikhonov-regularized backward). nref = 1
// gives back exactly asmTF.
func mediumTF(wl, z, nref float64, zeroEv, reg bool, alpha, limit float64) func(fx, fy float64) complex128 {
	k := 2 * math.Pi / wl * nref // wavenumber inside the reference medium
	az := math.Abs(z)
	backward := z < 0
	return func(fx, fy float64) complex128 {
		rho := math.Hypot(wl*fx, wl*fy) / nref
		if rho <= 1 {
			return cexpI(k * z * math.Sqrt((1-rho)*(1+rho)))
		}
		decay := k * az * math.Sqrt((rho-1)*(rho+1))
		if zeroEv || (limit > 0 && decay > limit) {
			return 0
		}
		if backward {
			if reg {
				a := math.Exp(decay)
				return complex(a/(1+alpha*alpha*a*a), 0)
			}
			return complex(math.Exp(decay), 0)
		}
		return complex(math.Exp(-decay), 0)
	}
}

// applyMediumPhaseRef multiplies the field by the reference-relative screen
// exp(i k0 (n - nref) dz) * exp(-k0 Im(n) dz) at every pixel of the slice
// centered at zc. With nref = 1 this is applyMediumPhase above.
func applyMediumPhaseRef(f *Field, n IndexFunc, zc, k, dz, nref float64) {
	nn := f.N
	for j := 0; j < nn; j++ {
		y := f.Y(j)
		for i := 0; i < nn; i++ {
			nv := n(f.X(i), y, zc)
			ph := k * (real(nv) - nref) * dz
			t := complex(math.Exp(-k*imag(nv)*dz), 0) * cexpI(ph)
			idx := j*nn + i
			f.Ex[idx] *= t
			if f.Polarized {
				f.Ey[idx] *= t
			}
		}
	}
}
