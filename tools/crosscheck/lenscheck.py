#!/usr/bin/env python3
"""Independent cross-check of Tetsuhiro-OZ first-order data and marginal ray
against the rayoptics package (0.9.5).

Input: the JSON written by
    go run ./tools/crosscheck --lens <prescription.json>
(object at infinity). The rayoptics model is rebuilt from exactly the same
numbers (radius, thickness, nd, vd, semi-aperture, stop index), so both sides
describe the same geometry and only the computation differs.

Quantities (all infinite-conjugate):
  EFL, BFL     parallel-input focal lengths; rayoptics fod.efl / fod.bfl.
  F/#          EFL / (2 EPR).  EPR_ref is the physical stop semi-aperture when
               the stop is surface 0 (the entrance pupil is then the stop
               itself); otherwise it is the OZ paraxial EPR, which is fed to
               rayoptics as its EPD.
  EPR          stop surface 0: OZ EPR against the physical stop semi-aperture.
               Interior stop: rayoptics' own paraxial marginal height at the
               physical stop, obtained with EPD = 2*EPR_OZ, must equal the
               stop semi-aperture (non-circular check of the paraxial EP).
  principal planes  front: distance from the first vertex; back: distance from
               the last vertex (rayoptics pp1 / ppk).
  marginal ray a parallel ray at height frac*EPR (frac = 0.999999 on both sides)
               traced by real ray tracing; its transverse intercept with the
               paraxial focal plane is compared.

Thresholds (fixed by the cross-check spec, never relaxed):
  EFL and F/# relative 1e-6; BFL absolute 1e-6 mm; principal planes and EPR
  absolute 1e-6 mm; marginal-ray intercept absolute 1e-3 mm.
  The marginal intercept carries a residual from the finite rayoptics object
  gap (OBJ_DISTANCE), see above; it is ~1e-9 mm at the default gap.

Output: a table, a machine-readable line "RESULT {json}" and the verdict
"LENS-CROSSCHECK OK|FAIL". Exit code 0 only when every gated check passes.
"""
import json
import math
import sys

from rayoptics.environment import OpticalModel, PupilSpec
from rayoptics.raytr.trace import trace_base

# rayoptics needs a finite object gap. A finite gap D tilts the marginal ray by
# h/D, which moves its focal-plane intercept by ~BFL*h/D (measured: 1.3e-7 mm at
# D = 1e10 mm, 1.3e-9 mm at 1e12 mm). EFL, BFL and EPR are unaffected.
OBJ_DISTANCE = 1e12
WAVELENGTH_D = 587.5618  # geo.Wavelengthd
FLAT_RADIUS = 1e18  # rayoptics radius of a plane surface

THRESH = {
    "efl_rel": 1e-6,
    "fno_rel": 1e-6,
    "bfl_abs_mm": 1e-6,
    "pp_abs_mm": 1e-6,
    "epr_abs_mm": 1e-6,
    "marginal_abs_mm": 1e-3,
}


def build(data, epr_ref):
    """Rebuild the prescription in rayoptics (radius mode, EPD = 2*epr_ref)."""
    opm = OpticalModel(radius_mode=True)  # rayoptics default is curvature
    sm = opm['seq_model']
    osp = opm['optical_spec']
    sm.gaps[0].thi = OBJ_DISTANCE
    for s in data["surfaces"]:
        r = s["radius_mm"]
        if r == 0 or math.isinf(r):
            r = FLAT_RADIUS
        sm.add_surface([r, s["thickness_mm"], s["nd"], s["vd"], s["semi_mm"]])
    sm.stop_surface = data["stop_index"] + 1  # rayoptics: object is surface 0
    osp['pupil'] = PupilSpec(osp, key=('object', 'epd'), value=2.0 * epr_ref)
    opm.update_model()
    return opm


def main(argv):
    path = argv[1] if len(argv) > 1 else "/tmp/oz_lens.json"
    with open(path) as fh:
        data = json.load(fh)
    n = len(data["surfaces"])
    stop = data["stop_index"]
    stop_first = stop == 0
    epr_oz = data["epr_mm"]
    stop_semi = data["stop_semi_mm"]
    epr_ref = stop_semi if stop_first else epr_oz

    opm = build(data, epr_ref)
    fod = opm['analysis_results']['parax_data'].fod
    efl, bfl = float(fod.efl), float(fod.bfl)
    epr_ro = float(fod.enp_radius)
    fno_ro = efl / (2.0 * epr_ref)
    pp1, ppk = float(fod.pp1), float(fod.ppk)
    # rayoptics' paraxial marginal height at the physical stop (ax_ray[k][0],
    # k = rayoptics surface index, 0 = object).
    stop_h_ro = float(opm['analysis_results']['parax_data'].ax_ray[stop + 1][0])

    # Marginal ray: rayoptics' own ray tracing, relative pupil (0, frac), on
    # the on-axis infinite field. ray[k] is the list of interactions in the
    # local frame of optical surface k (ray[0] = object, ray[1..n] = lens,
    # ray[n+1] = image); the last optical surface has its vertex at z = 0.
    frac = data["marginal_pupil_fraction"]
    fld = opm['optical_spec']['fov'].fields[0]
    pkg = trace_base(opm, [0.0, frac], fld, WAVELENGTH_D, apply_vignetting=False)
    ray = pkg.ray if hasattr(pkg, 'ray') else pkg[0]
    p_last, d_last = ray[n][0], ray[n][1]
    y_ro = float(p_last[1] + d_last[1] * (bfl - p_last[2]) / d_last[2])

    vertex_last = float(data["vertex_last_z_mm"])
    oz = {
        "efl": data["efl_mm"], "bfl": data["bfl_mm"], "fno": data["fno"],
        "epr": data["epr_mm"], "fpp": data["front_principal_z_mm"],
        "bpp": data["back_principal_z_mm"] - vertex_last,
        "marg": data["marginal_y_at_focus_mm"],
    }
    ro = {"efl": efl, "bfl": bfl, "fno": fno_ro, "epr": epr_ro,
          "fpp": pp1, "bpp": ppk, "marg": y_ro}

    checks = {}
    checks["efl_rel"] = abs(oz["efl"] - ro["efl"]) / abs(ro["efl"])
    checks["fno_rel"] = abs(oz["fno"] - ro["fno"]) / abs(ro["fno"])
    checks["bfl_abs_mm"] = abs(oz["bfl"] - ro["bfl"])
    checks["pp_abs_mm"] = max(abs(oz["fpp"] - ro["fpp"]), abs(oz["bpp"] - ro["bpp"]))
    if stop_first:
        checks["epr_abs_mm"] = abs(epr_oz - stop_semi)
    else:
        checks["epr_abs_mm"] = abs(stop_h_ro - stop_semi)
    checks["marginal_abs_mm"] = abs(oz["marg"] - ro["marg"])

    print(f"{'quantity':12s} {'OZ':>22s} {'rayoptics':>22s} {'|diff|':>12s}")
    rows = [
        ("EFL", oz["efl"], ro["efl"]), ("BFL", oz["bfl"], ro["bfl"]),
        ("F/#", oz["fno"], ro["fno"]), ("EPR", oz["epr"], ro["epr"]),
        ("front PP z", oz["fpp"], ro["fpp"]), ("back PP", oz["bpp"], ro["bpp"]),
        ("marginal y", oz["marg"], ro["marg"]),
    ]
    for name, a, b in rows:
        print(f"{name:12s} {a:22.15g} {b:22.15g} {abs(a - b):12.3e}")
    print(f"stop index {stop} (surface {'first' if stop_first else 'interior'}), "
          f"stop semi {stop_semi:.12g} mm, EPR_OZ {epr_oz:.12g} mm, "
          f"rayoptics paraxial stop height {stop_h_ro:.12g} mm, "
          f"marginal frac {frac!r}, OZ marginal status {data['marginal_status']}")

    ok = data["marginal_status"] == "ok"
    for k, v in checks.items():
        passed = v <= THRESH[k]
        ok = ok and passed
        print(f"  {k:14s} {v:.3e}  <= {THRESH[k]:.0e}  {'pass' if passed else 'FAIL'}")
    result = {
        "checks": {k: float(v) for k, v in checks.items()},
        "oz": {k: float(v) for k, v in oz.items()},
        "ro": {k: float(v) for k, v in ro.items()},
        "stop_height_ro": stop_h_ro,
        "epr_ref": float(epr_ref),
        "ok": bool(ok),
    }
    print("RESULT " + json.dumps(result, sort_keys=True))
    print("LENS-CROSSCHECK", "OK" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
