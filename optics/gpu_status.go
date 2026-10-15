package optics

import (
	"sync"
	"sync/atomic"
	"time"
)

// GPUStatus is a snapshot of the optional GPU FFT backend, intended for
// health endpoints and dashboards: it answers "is the GPU enabled, and is it
// actually being used?".
//
// Attempts counts transforms that the dispatch points offered to the GPU
// (i.e. the backend was enabled and the grid was large enough); Failures
// counts those that errored and fell back to the pure-Go CPU path. A
// non-zero, recent LastError together with Failures > 0 while Enabled is
// the signature of a silently degraded backend (for example a helper that
// refuses oversized requests), which is exactly what health checks should
// flag instead of letting a big run grind along on the CPU.
type GPUStatus struct {
	Available       bool   `json:"available"` // a GPU device / helper is present
	Enabled         bool   `json:"enabled"`   // GPU dispatch is switched on
	Backend         string `json:"backend"`   // human-readable backend description
	Attempts        uint64 `json:"attempts"`  // transforms offered to the GPU
	Failures        uint64 `json:"failures"`  // transforms that fell back to the CPU
	LastError       string `json:"last_error,omitempty"`
	LastErrorUnixMs int64  `json:"last_error_unix_ms,omitempty"`
}

var (
	gpuStatAttempts  atomic.Uint64
	gpuStatFailures  atomic.Uint64
	gpuStatMu        sync.Mutex
	gpuStatLastErr   string
	gpuStatLastErrMs int64
)

// GPUState returns a snapshot of the backend state. It is cheap (two atomic
// loads and a short mutex) and safe for concurrent use from any goroutine.
func GPUState() GPUStatus {
	st := GPUStatus{
		Available: GPUAvailable(),
		Enabled:   GPUEnabled(),
		Backend:   GPUBackendInfo(),
		Attempts:  gpuStatAttempts.Load(),
		Failures:  gpuStatFailures.Load(),
	}
	gpuStatMu.Lock()
	st.LastError = gpuStatLastErr
	st.LastErrorUnixMs = gpuStatLastErrMs
	gpuStatMu.Unlock()
	return st
}

// gpuNoteAttempt records that a dispatch point is about to attempt a GPU
// transform (called only when the backend is enabled and the grid is big
// enough for it).
func gpuNoteAttempt() { gpuStatAttempts.Add(1) }

// gpuNoteFailure records one failed GPU transform plus the error that made
// the caller fall back to the CPU path.
func gpuNoteFailure(err error) {
	gpuStatFailures.Add(1)
	if err != nil {
		gpuStatMu.Lock()
		gpuStatLastErr = err.Error()
		gpuStatLastErrMs = time.Now().UnixMilli()
		gpuStatMu.Unlock()
	}
}

// resetGPUStats clears the counters (used by tests).
func resetGPUStats() {
	gpuStatAttempts.Store(0)
	gpuStatFailures.Store(0)
	gpuStatMu.Lock()
	gpuStatLastErr = ""
	gpuStatLastErrMs = 0
	gpuStatMu.Unlock()
}
