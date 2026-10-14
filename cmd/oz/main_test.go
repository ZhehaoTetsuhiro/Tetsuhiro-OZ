package main

import (
	"path/filepath"
	"testing"

	"tetsuhiro-oz/lens"
)

// TestDesignEvalRoundTrip drives the CLI wiring: design writes a
// prescription that eval can load, describe and re-read.
func TestDesignEvalRoundTrip(t *testing.T) {
	dir := t.TempDir()
	out := filepath.Join(dir, "design.json")
	if err := cmdDesign([]string{"achromat", "--efl", "50", "--semi", "6", "--out", out}); err != nil {
		t.Fatal(err)
	}
	p, err := lens.LoadJSON(out)
	if err != nil {
		t.Fatal(err)
	}
	if len(p.Cards) != 3 {
		t.Fatalf("designed prescription has %d cards, want 3", len(p.Cards))
	}
	if err := cmdEval([]string{out}); err != nil {
		t.Fatal(err)
	}
}

// TestTolSmoke exercises the tolerance subcommand end to end with a small
// run count.
func TestTolSmoke(t *testing.T) {
	dir := t.TempDir()
	out := filepath.Join(dir, "design.json")
	if err := cmdDesign([]string{"achromat", "--efl", "50", "--semi", "6", "--no-polish", "--out", out}); err != nil {
		t.Fatal(err)
	}
	rep := filepath.Join(dir, "rep.json")
	err := cmdTol([]string{
		"--radius", "0.05", "--index", "0.0005",
		"--runs", "60", "--seed", "5", "--target-efl", "50", "--max-focal", "0.1",
		"--json", rep,
		out,
	})
	if err != nil {
		t.Fatal(err)
	}
	if _, err := lens.LoadJSON(out); err != nil { // prescription still intact
		t.Fatal(err)
	}
}

// TestUsagePaths checks the command dispatch guards.
func TestUsagePaths(t *testing.T) {
	if err := cmdDesign([]string{"achromat"}); err == nil {
		t.Fatal("design without --efl accepted")
	}
	if err := cmdDesign([]string{"triplet", "--efl", "50"}); err == nil {
		t.Fatal("unknown lens type accepted")
	}
}
