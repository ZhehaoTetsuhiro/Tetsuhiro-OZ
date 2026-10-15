#!/usr/bin/env python3
"""Protocol test for the wos-musa-fft daemon.

Connects to the UNIX socket, runs 1-D and 2-D forward transforms plus an
inverse round trip, and verifies results against a reference DFT computed
with cmath. Exit code 0 = all good.

Usage: python3 test_client.py [socket_path]
"""
import cmath
import random
import socket
import struct
import sys
import time

MAGIC = 0x46465557
OP_FFT2D, OP_FFT1D = 1, 2
DIR_FWD, DIR_INV = 0, 1
HDR_REQ = "<IIIIQ"  # magic, op, dir, n, bytes
HDR_REP = "<IIQ"    # status, rsvd, bytes


def connect_with_retry(path, tries=50, delay=0.1):
    last = None
    for _ in range(tries):
        try:
            s = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
            s.connect(path)
            return s
        except OSError as e:  # daemon not up yet
            last = e
            time.sleep(delay)
    raise SystemExit(f"could not connect to {path}: {last}")


def recv_exact(s, n):
    buf = b""
    while len(buf) < n:
        chunk = s.recv(n - len(buf))
        if not chunk:
            raise SystemExit("connection closed mid-reply")
        buf += chunk
    return buf


def request(s, op, direction, n, data):
    payload = b"".join(struct.pack("<dd", c.real, c.imag) for c in data)
    s.sendall(struct.pack(HDR_REQ, MAGIC, op, direction, n, len(payload)) + payload)
    status, _rsvd, nbytes = struct.unpack(HDR_REP, recv_exact(s, 16))
    if status != 0:
        raise SystemExit(f"daemon error status={status}")
    buf = recv_exact(s, nbytes)
    return [complex(*struct.unpack_from("<dd", buf, 16 * i)) for i in range(nbytes // 16)]


def dft(a, inv=False):
    n = len(a)
    sign = 1.0 if inv else -1.0
    return [
        sum(a[j] * cmath.exp(sign * 2j * cmath.pi * k * j / n) for j in range(n))
        for k in range(n)
    ]


def dft2(a, n):
    rows = [dft(a[r * n:(r + 1) * n]) for r in range(n)]
    out = [0j] * (n * n)
    for c in range(n):
        col = dft([rows[r][c] for r in range(n)])
        for r in range(n):
            out[r * n + c] = col[r]
    return out


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "/tmp/wos-musa-fft.sock"
    s = connect_with_retry(path)
    random.seed(20261015)

    n = 16
    a = [complex(random.random(), random.random()) for _ in range(n)]
    got = request(s, OP_FFT1D, DIR_FWD, n, a)
    ref = dft(a)
    err = max(abs(got[i] - ref[i]) for i in range(n))
    print(f"1D forward  n={n:5d}  max |daemon - DFT| = {err:.3e}")
    assert err < 1e-10, "1D forward mismatch"

    back = request(s, OP_FFT1D, DIR_INV, n, got)
    rterr = max(abs(back[i] / n - a[i]) for i in range(n))
    print(f"1D round-trip          max |inv/n - in| = {rterr:.3e}")
    assert rterr < 1e-12, "1D round trip mismatch"

    n2 = 8
    m = [complex(random.random(), random.random()) for _ in range(n2 * n2)]
    got2 = request(s, OP_FFT2D, DIR_FWD, n2, m)
    ref2 = dft2(m, n2)
    err2 = max(abs(got2[i] - ref2[i]) for i in range(n2 * n2))
    print(f"2D forward  n={n2}x{n2}   max |daemon - DFT| = {err2:.3e}")
    assert err2 < 1e-10, "2D forward mismatch"

    print("ALL OK")
    s.close()


if __name__ == "__main__":
    main()
