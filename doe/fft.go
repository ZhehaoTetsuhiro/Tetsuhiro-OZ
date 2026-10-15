package doe

import (
	"fmt"
	"math"
	"math/bits"
	"runtime"
	"sync"
)

// fftPlan holds the bit-reversal permutation and twiddle factors of one
// power-of-two length. Plans are immutable once built.
type fftPlan struct {
	n   int
	rev []int
	tw  []complex128 // tw[k] = exp(−2πi k/n), k < n/2
}

func newPlan(n int) *fftPlan {
	p := &fftPlan{n: n, rev: make([]int, n), tw: make([]complex128, n/2)}
	lg := bits.TrailingZeros(uint(n))
	for i := 0; i < n; i++ {
		p.rev[i] = int(bits.Reverse(uint(i)) >> (bits.UintSize - lg))
	}
	for k := range p.tw {
		ang := -2 * math.Pi * float64(k) / float64(n)
		p.tw[k] = complex(math.Cos(ang), math.Sin(ang))
	}
	return p
}

// transform performs an in-place radix-2 FFT of a (len(a) == p.n). The
// inverse transform uses conjugated twiddles and divides by n.
func (p *fftPlan) transform(a []complex128, inverse bool) {
	n := p.n
	for i := 0; i < n; i++ {
		if j := p.rev[i]; i < j {
			a[i], a[j] = a[j], a[i]
		}
	}
	for size := 2; size <= n; size <<= 1 {
		half := size >> 1
		step := n / size
		for start := 0; start < n; start += size {
			for k := 0; k < half; k++ {
				w := p.tw[k*step]
				if inverse {
					w = complex(real(w), -imag(w))
				}
				u := a[start+k]
				v := a[start+k+half] * w
				a[start+k] = u + v
				a[start+k+half] = u - v
			}
		}
	}
	if inverse {
		s := complex(1/float64(n), 0)
		for i := range a {
			a[i] *= s
		}
	}
}

// FFT2 transforms an n×n row-major array in place (n a power of two). The
// forward transform is unnormalised, U(k) = Σ_x u(x) exp(−2πi k·x/n); the
// inverse divides by n². Rows and columns are processed in parallel over the
// available CPUs; the result does not depend on the scheduling.
func FFT2(a []complex128, n int, inverse bool) error {
	if n < 1 || n&(n-1) != 0 || len(a) != n*n {
		return fmt.Errorf("doe: FFT2 needs an n×n array with n a power of two (n=%d, len=%d)", n, len(a))
	}
	p := newPlan(n)
	workers := runtime.GOMAXPROCS(0)
	if workers > n {
		workers = n
	}
	parallel := func(fn func(lo, hi int)) {
		var wg sync.WaitGroup
		for w := 0; w < workers; w++ {
			lo := w * n / workers
			hi := (w + 1) * n / workers
			wg.Add(1)
			go func(lo, hi int) {
				defer wg.Done()
				fn(lo, hi)
			}(lo, hi)
		}
		wg.Wait()
	}
	// Rows (contiguous).
	parallel(func(lo, hi int) {
		for r := lo; r < hi; r++ {
			p.transform(a[r*n:(r+1)*n], inverse)
		}
	})
	// Columns: gather into a per-worker buffer, transform, scatter back.
	parallel(func(lo, hi int) {
		col := make([]complex128, n)
		for c := lo; c < hi; c++ {
			for r := 0; r < n; r++ {
				col[r] = a[r*n+c]
			}
			p.transform(col, inverse)
			for r := 0; r < n; r++ {
				a[r*n+c] = col[r]
			}
		}
	})
	return nil
}

// fft1 transforms a 1-D slice in place (length a power of two).
func fft1(a []complex128, inverse bool) error {
	n := len(a)
	if n < 1 || n&(n-1) != 0 {
		return fmt.Errorf("doe: FFT length %d is not a power of two", n)
	}
	newPlan(n).transform(a, inverse)
	return nil
}
