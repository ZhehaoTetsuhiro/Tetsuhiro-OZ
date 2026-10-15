package main

import (
	"os"
	"strings"
	"testing"
)

// TestVersionMatchesVERSIONFile keeps the compiled-in version (served by
// /api/version and printed by `oz version`) equal to the repository VERSION
// file. The integrator owns VERSION; this test makes a mismatch a failure
// instead of a silent drift.
func TestVersionMatchesVERSIONFile(t *testing.T) {
	buf, err := os.ReadFile("../../VERSION")
	if err != nil {
		t.Fatal(err)
	}
	want := strings.TrimSpace(string(buf))
	if version != want {
		t.Fatalf("version constant %q does not match VERSION file %q", version, want)
	}
}
