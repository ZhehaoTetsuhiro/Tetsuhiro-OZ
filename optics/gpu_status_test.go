package optics

import (
	"fmt"
	"testing"
)

// TestGPUStateCounters checks the degradation counters that back the health
// endpoint and the server's fallback warning: offered transforms increment
// Attempts, failed ones increment Failures and record the error text and
// timestamp.
func TestGPUStateCounters(t *testing.T) {
	resetGPUStats()
	defer resetGPUStats()

	gpuNoteAttempt()
	gpuNoteAttempt()
	gpuNoteFailure(fmt.Errorf("boom: oversized request"))

	st := GPUState()
	if st.Attempts != 2 || st.Failures != 1 {
		t.Fatalf("counters: attempts=%d failures=%d, want 2/1", st.Attempts, st.Failures)
	}
	if st.LastError != "boom: oversized request" {
		t.Fatalf("last error %q, want %q", st.LastError, "boom: oversized request")
	}
	if st.LastErrorUnixMs <= 0 {
		t.Fatalf("last error timestamp not set")
	}
	if st.Backend == "" {
		t.Fatalf("backend description is empty")
	}
}

// TestGPUStateCleanAfterReset pins the reset helper used by the dispatch
// test (and by operators clearing state after a restart).
func TestGPUStateCleanAfterReset(t *testing.T) {
	gpuNoteAttempt()
	gpuNoteFailure(fmt.Errorf("x"))
	resetGPUStats()
	st := GPUState()
	if st.Attempts != 0 || st.Failures != 0 || st.LastError != "" || st.LastErrorUnixMs != 0 {
		t.Fatalf("reset left state behind: %+v", st)
	}
}
