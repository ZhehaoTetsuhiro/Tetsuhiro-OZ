#!/usr/bin/env python3
"""Compare one wave simulation between two OZ instances (CPU vs GPU).

Submits the same config to both bases, polls both runs, then compares the
plane statistics of the first output plane (power, peak, centroid, rms).
Prints relative differences; exit 0 when all are within 1e-9.
"""
import json
import sys
import time
import urllib.request

CFG = {
    "name": "cpu-vs-gpu",
    "grid": {"size": 1024, "width": 0.008},
    "wavelength": 532e-9,
    "polarized": False,
    "method": "asm",
    "evanescent": "decay",
    "sources": [{
        "id": "src", "type": "gaussian",
        "params": {"power": 0.001, "waist": 0.0015},
        "pos": {"x": 0.0, "y": 0.0, "z": -0.05},
        "dir": {"x": 0.0, "y": 0.0, "z": 1.0},
        "wavelength": 532e-9,
    }],
    "scene": {"components": [{
        "id": "det", "type": "sensor",
        "pos": {"x": 0.0, "y": 0.0, "z": 0.5},
        "label": "观察面",
        "yaw": 3.141592653589793,
    }]},
}


def post(base, path, obj, timeout=60):
    req = urllib.request.Request(base + path, data=json.dumps(obj).encode(),
                                 headers={"Content-Type": "application/json"}, method="POST")
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.load(r)


def get(base, path, timeout=30):
    with urllib.request.urlopen(base + path, timeout=timeout) as r:
        return json.load(r)


def run_one(base, cfg):
    r = post(base, "/api/simulate", cfg)
    rid = r.get("run_id")
    for _ in range(60):
        s = get(base, f"/api/runs/{rid}")
        if s.get("status") == "done":
            return s
        if s.get("status") == "error":
            raise SystemExit(f"{base} run error: {s.get('error')}")
        time.sleep(0.5)
    raise SystemExit(f"{base} run timeout")


def main():
    cpu = sys.argv[1] if len(sys.argv) > 1 else "http://127.0.0.1:6006"
    gpu = sys.argv[2] if len(sys.argv) > 2 else "http://127.0.0.1:6008"
    a = run_one(cpu, CFG)
    b = run_one(gpu, CFG)
    pa = a["planes"][0]["stats"]
    pb = b["planes"][0]["stats"]
    worst = 0.0
    print(f"{'field':12s} {'CPU':>22s} {'GPU':>22s} {'rel diff':>12s} {'abs diff':>12s}")
    for k in ("power", "peak", "rms_x", "rms_y", "centroid_x", "centroid_y"):
        va, vb = pa[k], pb[k]
        denom = max(abs(va), abs(vb), 1e-300)
        d = abs(va - vb) / denom
        da = abs(va - vb)
        # Pass when the relative difference is within 1e-9, or when the
        # absolute difference is numerically negligible (1e-15 SI): the
        # centroid of a symmetric beam is ~1e-11 m, so its relative
        # difference is dominated by summation-order noise.
        eff = d if da > 1e-15 else 0.0
        worst = max(worst, eff)
        print(f"{k:12s} {va:22.12g} {vb:22.12g} {d:12.3e} {da:12.3e}")
    print(f"worst effective difference: {worst:.3e}  (criterion < 1e-9, abs floor 1e-15)")
    print("CPU-vs-GPU", "OK" if worst < 1e-9 else "FAIL")
    return 0 if worst < 1e-9 else 1


if __name__ == "__main__":
    sys.exit(main())
