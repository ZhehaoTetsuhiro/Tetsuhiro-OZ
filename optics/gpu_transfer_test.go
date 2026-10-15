//go:build cuda

package optics

// Transfer-bandwidth benchmarks for the CUDA backend: a full two-way
// (H2D + D2H) copy of an n*n complex128 buffer, pageable vs pinned.
//
// These live in their own file because they exercise the in-process
// staging internals (gpuAlloc/gpuHostAlloc/gpuCopyH2D/gpuCopyD2H) which
// only exist in the CUDA backend; gpu_test.go now also builds for the
// MUSA helper backend (`-tags musa`), where the device is owned by a
// separate process and no in-process staging exists.

import (
	"testing"
	"unsafe"
)

// benchTransfer measures a full two-way (H2D + D2H) copy of an n*n complex128
// buffer, pageable vs pinned.
func benchTransfer(b *testing.B, n int, pinned bool) {
	nn := n * n
	nb := nn * 16
	a := make([]complex128, nn)
	d, err := gpuAlloc(nb)
	if err != nil {
		b.Fatal(err)
	}
	defer gpuRelease(d, nb)
	var hp unsafe.Pointer
	if pinned {
		hp, err = gpuHostAlloc(nb)
		if err != nil {
			b.Fatal(err)
		}
		defer gpuHostFree(hp)
	}
	src := nb
	b.SetBytes(int64(nb) * 2)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if pinned {
			gpuCopyH2D(d, hp, src)
			gpuCopyD2H(hp, d, src)
			copy(a, unsafe.Slice((*complex128)(hp), nn))
		} else {
			gpuCopyH2D(d, unsafe.Pointer(&a[0]), src)
			gpuCopyD2H(unsafe.Pointer(&a[0]), d, src)
		}
	}
}

func BenchmarkTransferPageable4096(b *testing.B) { benchTransfer(b, 4096, false) }
func BenchmarkTransferPinned4096(b *testing.B)   { benchTransfer(b, 4096, true) }
