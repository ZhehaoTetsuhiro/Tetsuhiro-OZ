#!/usr/bin/env python3
"""Cross-check the Tetsuhiro-OZ coating TMM against the Python `tmm` package.

Feed the same complex refractive indices (exported verbatim from the Go
side) and the same stacks into Steve Byrnes' `tmm` (an independently
written, widely used implementation) and compare every R/T point. This
isolates the transfer-matrix algorithm from the material models.

Usage: python3 crosscheck.py [oz_tmm.json]
Exit code 0 = agreement within 1e-9 in R and T.
"""
import json
import math
import sys

from tmm import coh_tmm


def nk(tables, name, lam):
    tab = tables[name]
    if lam <= tab[0][0]:
        return complex(tab[0][1], tab[0][2])
    if lam >= tab[-1][0]:
        return complex(tab[-1][1], tab[-1][2])
    lo, hi = 0, len(tab) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if tab[mid][0] <= lam:
            lo = mid
        else:
            hi = mid
    l0, n0, k0 = tab[lo]
    l1, n1, k1 = tab[hi]
    w = (lam - l0) / (l1 - l0)
    return complex(n0 + w * (n1 - n0), k0 + w * (k1 - k0))


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "/tmp/oz_tmm.json"
    data = json.load(open(path))
    mats = data["materials"]
    worst_r, worst_t = 0.0, 0.0
    total = 0
    for case in data["cases"]:
        case_worst = 0.0
        for pt in case["points"]:
            lam, ang, pol = pt["lam_nm"], pt["angle_deg"], pt["pol"]
            n_list = [nk(mats, case["incident"], lam)]
            n_list += [nk(mats, l["mat"], lam) for l in case["layers"]]
            n_list += [nk(mats, case["substrate"], lam)]
            d_list = [float("inf")] + [l["d_nm"] for l in case["layers"]] + [float("inf")]  # semi-infinite media
            res = coh_tmm(pol, n_list, d_list, math.radians(ang), lam)
            dr = abs(res["R"] - pt["R"])
            dt = abs(res["T"] - pt["T"])
            worst_r = max(worst_r, dr)
            worst_t = max(worst_t, dt)
            case_worst = max(case_worst, dr, dt)
            total += 1
        print(f"{case['name']}: {len(case['points'])} points, worst |d| = {case_worst:.3e}")
    print(f"checked {total} points; worst |dR| = {worst_r:.3e}, worst |dT| = {worst_t:.3e}")
    ok = worst_r < 1e-9 and worst_t < 1e-9
    print("CROSSCHECK", "OK" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
