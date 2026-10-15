//go:build musa && !cuda

// GPU acceleration backend for the wave-optics kernel on Moore Threads
// (MUSA) devices.
//
// Unlike gpu.go (CUDA/cuFFT, in-process via cgo), this backend talks to a
// small helper process — tools/musa-fft/fft_helper.c, binary
// `wos-musa-fft` — over a UNIX socket. Reason: the MUSA userspace driver
// reserves its device address window [4GB, ~873GB) inside the process at
// init time and aborts with musaErrorInvalidDevice ("invalid device
// ordinal") if any foreign VMA lives in that range. The Go runtime's heap
// arena sits at 0xc000000000 (768 GB) — permanently inside that window —
// so no Go process can initialize the device, while a plain C process can.
// The helper (a C process) owns the MUSA context and serves FFTs; the Go
// side never touches the device. Full analysis with evidence:
// docs/notes-musa-rootcause.md.
//
// Everything else is identical to the CUDA backend: same dispatch points
// (fft2D / fft1DAny), same fallback semantics (any error -> pure-Go CPU
// path), same thresholds and the same SetGPU/GPUEnabled/GPUBackendInfo
// API, so every propagation method (asm, asm_pad, fresnel_*, fraunhofer),
// element application (thin lenses, Berreman, birefringent), coherence and
// field metrics that touch the FFT automatically use the GPU when enabled.
//
// This file is pure Go (no cgo): `go build -tags musa` needs no MUSA SDK
// at build time; only running the GPU path needs the helper binary and a
// MUSA device. Helper discovery (first hit wins):
//
//	1. $WOS_MUSA_HELPER            — explicit binary path
//	2. "wos-musa-fft" next to the wos executable
//	3. "wos-musa-fft" in $PATH
//
// Socket path: $WOS_MUSA_SOCKET, default /tmp/wos-musa-fft.sock. An
// already-running daemon on that socket is reused; otherwise one is
// spawned lazily on the first transform request.
//
// Protocol (little-endian; both sides are x86-64 on this platform):
//
//	request  : u32 magic("WUFF"), u32 op(0=info,1=fft2d,2=fft1d),
//	           u32 dir(0=fwd,1=inv), u32 n, u64 bytes, then payload
//	reply    : u32 status(0=ok), u32 rsvd, u64 bytes, then payload
//
// Inverse transforms come back unnormalized; this file applies the 1/n
// (1-D) / 1/n^2 (2-D) scaling, matching the pure-Go fft conventions.

package optics

import (
	"encoding/binary"
	"fmt"
	"io"
	"net"
	"os"
	"os/exec"
	"path/filepath"
	"sync"
	"sync/atomic"
	"time"
	"unsafe"
)

const (
	helperMagic = 0x46465557 // "WUFF"
	hOpInfo     = 0
	hOpFFT2D    = 1
	hOpFFT1D    = 2
	hDirFwd     = 0
	hDirInv     = 1

	helperSockDefault = "/tmp/wos-musa-fft.sock"
	helperConnectTO   = 300 * time.Millisecond // probe an already-running daemon
	helperSpawnTO     = 10 * time.Second       // wait for a freshly spawned daemon
	helperIOTimeout   = 5 * time.Minute        // per-request guard
	helperIdleExitSec = "900"                  // daemon reaps itself when unused
)

var (
	gpuOnce    sync.Once
	gpuPresent bool
	gpuEnabled atomic.Bool

	// helperMu serializes all helper traffic and guards helperConn,
	// helperCmd, helperPath, helperSock and gpuName. The kernel serializes
	// GPU calls anyway (same rationale as gpu.go's gpuExecMu).
	helperMu   sync.Mutex
	helperConn net.Conn
	helperCmd  *exec.Cmd
	helperPath string
	helperSock string
	gpuName    string
)

func gpuInit() {
	helperMu.Lock()
	defer helperMu.Unlock()

	helperSock = os.Getenv("WOS_MUSA_SOCKET")
	if helperSock == "" {
		helperSock = helperSockDefault
	}
	// Passive probe: reuse a running daemon if one is listening. Otherwise
	// report available when the helper binary exists and let the first
	// transform request spawn it (spawning during init would launch a
	// process even for callers that merely ask GPUAvailable).
	helperPath = helperFindBinary()
	if conn, ok := helperDial(); ok {
		helperConn = conn
		gpuName = helperQueryName(conn)
		gpuPresent = true
	} else if helperPath != "" {
		gpuName = "helper " + filepath.Base(helperPath)
		gpuPresent = true
	}
	if os.Getenv("WOS_GPU") == "1" {
		gpuEnabled.Store(true)
	}
}

// GPUAvailable reports whether a MUSA device is reachable through the
// helper (build tag `musa` compiled in and either a daemon or the helper
// binary present).
func GPUAvailable() bool { gpuOnce.Do(gpuInit); return gpuPresent }

// GPUEnabled reports whether GPU acceleration is currently active.
func GPUEnabled() bool { gpuOnce.Do(gpuInit); return gpuPresent && gpuEnabled.Load() }

// SetGPU turns GPU acceleration on or off. It returns false if the request
// cannot be honoured (no helper / no MUSA device / backend not compiled in).
func SetGPU(on bool) bool {
	gpuOnce.Do(gpuInit)
	if on && !gpuPresent {
		return false
	}
	gpuEnabled.Store(on)
	return true
}

// GPUBackendInfo returns a human-readable description of the active backend.
func GPUBackendInfo() string {
	gpuOnce.Do(gpuInit)
	helperMu.Lock()
	name := gpuName
	helperMu.Unlock()
	if !gpuPresent {
		return "CPU (pure Go; no MUSA helper reachable)"
	}
	if gpuEnabled.Load() {
		return "GPU: MUSA/muFFT via wos-musa-fft helper on " + name
	}
	return "GPU available (MUSA/muFFT helper on " + name + ") but disabled; CPU in use"
}

func gpuReady() bool { gpuOnce.Do(gpuInit); return gpuPresent && gpuEnabled.Load() }

// helperFindBinary resolves the helper executable: $WOS_MUSA_HELPER,
// then next to the running executable, then $PATH.
func helperFindBinary() string {
	if p := os.Getenv("WOS_MUSA_HELPER"); p != "" {
		if fi, err := os.Stat(p); err == nil && !fi.IsDir() {
			return p
		}
	}
	if exe, err := os.Executable(); err == nil {
		p := filepath.Join(filepath.Dir(exe), "wos-musa-fft")
		if fi, err := os.Stat(p); err == nil && !fi.IsDir() {
			return p
		}
	}
	if p, err := exec.LookPath("wos-musa-fft"); err == nil {
		return p
	}
	return ""
}

// helperDial connects to the daemon socket. Callers must hold helperMu.
func helperDial() (net.Conn, bool) {
	c, err := net.DialTimeout("unix", helperSock, helperConnectTO)
	if err != nil {
		return nil, false
	}
	return c, true
}

// helperQueryName asks the daemon for the device name (op=info). A daemon
// that does not answer is still usable, so failures are non-fatal.
func helperQueryName(c net.Conn) string {
	var hdr [24]byte
	binary.LittleEndian.PutUint32(hdr[0:], helperMagic)
	binary.LittleEndian.PutUint32(hdr[4:], hOpInfo)
	_ = c.SetDeadline(time.Now().Add(helperConnectTO))
	if _, err := c.Write(hdr[:]); err != nil {
		_ = c.SetDeadline(time.Time{})
		return "wos-musa-fft helper"
	}
	var rep [16]byte
	if _, err := io.ReadFull(c, rep[:]); err != nil {
		_ = c.SetDeadline(time.Time{})
		return "wos-musa-fft helper"
	}
	nb := binary.LittleEndian.Uint64(rep[8:])
	if binary.LittleEndian.Uint32(rep[0:]) != 0 || nb == 0 || nb > 4096 {
		_ = c.SetDeadline(time.Time{})
		return "wos-musa-fft helper"
	}
	buf := make([]byte, nb)
	if _, err := io.ReadFull(c, buf); err != nil {
		_ = c.SetDeadline(time.Time{})
		return "wos-musa-fft helper"
	}
	_ = c.SetDeadline(time.Time{})
	return string(buf)
}

// helperAcquireLocked returns a live daemon connection, spawning the
// helper if none is running. Callers must hold helperMu.
func helperAcquireLocked() (net.Conn, error) {
	if helperConn != nil {
		return helperConn, nil
	}
	if conn, ok := helperDial(); ok {
		helperConn = conn
		gpuName = helperQueryName(conn)
		return conn, nil
	}
	if helperPath == "" {
		return nil, fmt.Errorf("musa helper binary not found")
	}
	cmd := exec.Command(helperPath, "--socket", helperSock, "--idle-exit", helperIdleExitSec)
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	if err := cmd.Start(); err != nil {
		return nil, fmt.Errorf("spawn musa helper: %w", err)
	}
	helperCmd = cmd
	deadline := time.Now().Add(helperSpawnTO)
	for time.Now().Before(deadline) {
		if conn, ok := helperDial(); ok {
			helperConn = conn
			gpuName = helperQueryName(conn)
			return conn, nil
		}
		time.Sleep(100 * time.Millisecond)
	}
	_ = cmd.Process.Kill()
	helperCmd = nil
	return nil, fmt.Errorf("musa helper did not come up within %s", helperSpawnTO)
}

// helperDropLocked closes a broken connection so the next request
// re-acquires (reusing or respawning the daemon). Callers hold helperMu.
func helperDropLocked() {
	if helperConn != nil {
		_ = helperConn.Close()
		helperConn = nil
	}
}

// helperRequest runs one transform through the daemon and copies the
// result back into a. op/dir are protocol values; n is the 1-D length or
// the 2-D side. Returns an error on any protocol/IO failure — the caller
// then falls back to the pure-Go path.
func helperRequest(op, dir, n int, a []complex128) error {
	helperMu.Lock()
	defer helperMu.Unlock()
	c, err := helperAcquireLocked()
	if err != nil {
		return err
	}
	if len(a) == 0 {
		return nil
	}
	nb := len(a) * int(unsafe.Sizeof(complex128(0)))

	var hdr [24]byte
	binary.LittleEndian.PutUint32(hdr[0:], helperMagic)
	binary.LittleEndian.PutUint32(hdr[4:], uint32(op))
	binary.LittleEndian.PutUint32(hdr[8:], uint32(dir))
	binary.LittleEndian.PutUint32(hdr[12:], uint32(n))
	binary.LittleEndian.PutUint64(hdr[16:], uint64(nb))

	// Zero-copy view of the complex buffer (x86-64: little-endian doubles,
	// no padding between real/imaginary of complex128).
	payload := unsafe.Slice((*byte)(unsafe.Pointer(&a[0])), nb)

	_ = c.SetDeadline(time.Now().Add(helperIOTimeout))
	if _, err := c.Write(hdr[:]); err != nil {
		helperDropLocked()
		return fmt.Errorf("musa helper header: %w", err)
	}
	if _, err := c.Write(payload); err != nil {
		helperDropLocked()
		return fmt.Errorf("musa helper payload: %w", err)
	}
	var rep [16]byte
	if _, err := io.ReadFull(c, rep[:]); err != nil {
		helperDropLocked()
		return fmt.Errorf("musa helper reply: %w", err)
	}
	if status := binary.LittleEndian.Uint32(rep[0:]); status != 0 {
		helperDropLocked()
		return fmt.Errorf("musa helper status %d", status)
	}
	if got := binary.LittleEndian.Uint64(rep[8:]); got != uint64(nb) {
		helperDropLocked()
		return fmt.Errorf("musa helper reply size %d != %d", got, nb)
	}
	if _, err := io.ReadFull(c, payload); err != nil {
		helperDropLocked()
		return fmt.Errorf("musa helper result: %w", err)
	}
	_ = c.SetDeadline(time.Time{})
	return nil
}

// gpuFFT2D performs an in-place 2-D transform of a length n*n (row-major)
// via the helper, matching the sign/normalization convention of the
// pure-Go fft2D.
func gpuFFT2D(a []complex128, n int, inverse bool) error {
	if !gpuReady() {
		return errGPUUnavailable
	}
	if len(a) != n*n {
		return fmt.Errorf("gpuFFT2D: length %d != %d*%d", len(a), n, n)
	}
	dir := hDirFwd
	if inverse {
		dir = hDirInv
	}
	if err := helperRequest(hOpFFT2D, dir, n, a); err != nil {
		return err
	}
	if inverse {
		s := complex(1/float64(n*n), 0)
		for i := range a {
			a[i] *= s
		}
	}
	return nil
}

// gpuFFT1D performs an in-place 1-D transform of arbitrary length via the
// helper, matching the sign/normalization convention of the pure-Go
// fft1DAny.
func gpuFFT1D(a []complex128, inverse bool) error {
	if !gpuReady() {
		return errGPUUnavailable
	}
	n := len(a)
	if n == 0 {
		return nil
	}
	dir := hDirFwd
	if inverse {
		dir = hDirInv
	}
	if err := helperRequest(hOpFFT1D, dir, n, a); err != nil {
		return err
	}
	if inverse {
		s := complex(1/float64(n), 0)
		for i := range a {
			a[i] *= s
		}
	}
	return nil
}
