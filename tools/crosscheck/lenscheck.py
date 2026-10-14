#!/usr/bin/env python3
"""Cross-check Tetsuhiro-OZ lens first-order data + marginal ray against the
independent rayoptics package (v0.9.5).

Input: the JSON emitted by
    go run ./tools/crosscheck --lens <prescription.json>

The comparison uses the SAME nd/Vd values (our Sellmeier evaluated at the
d-line), so it isolates the paraxial engine and the real-ray tracer:
  * EFL(d), BFL(d), F/# must agree to ~1e-6 relative;
  * the marginal ray's transverse intercept at the paraxial focus plane must
    agree (spherical-aberration residual included) within a small absolute
    band.
"""
import json
import math
import sys

from rayoptics.environment import *  # noqa: F401,F403 (OpticalModel, PupilSpec, trace)


def build(data):
    opm = OpticalModel()
    sm = opm['seq_model']
    osp = opm['optical_spec']
    # The fresh model already carries [Obj, Img]; the object distance is the
    # first gap (1e10 ~ infinity). Do NOT add another object surface: a
    # duplicate would leave the object gap at zero and update_model() would
    # divide by zero (obj2enp_dist = 0) inside rayoptics.
    sm.gaps[0].thi = data.get("object_distance_mm", 1e10)
    for s in data["surfaces"]:
        r = s["radius_mm"]
        if r == 0 or math.isinf(r):
            r = 1e18
        sm.add_surface([r, s["thickness_mm"], s["nd"], s["vd"], s["semi_mm"]])
    sm.stop_surface = 1
    osp['aperture'] = PupilSpec(osp, key=('object', 'epd'), value=2 * data["pupil_radius"])
    opm.update_model()
    return opm, sm


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "/tmp/oz_lens.json"
    data = json.load(open(path))
    opm, sm = build(data)
    fod = opm['analysis_results']['parax_data'].fod
    efl = fod.efl
    bfd = getattr(fod, 'bfl', None)
    if bfd is None:
        bfd = getattr(fod, 'bfd', None)
    fno = getattr(fod, 'fno', None)

    # Marginal ray: the real ray from the on-axis object point through the
    # rim of the entrance pupil — rayoptics' own definition (trace_base with a
    # relative-pupil coordinate), traced with the same pupil fraction the Go
    # side used.
    from rayoptics.raytr.trace import trace_base
    if data.get("marginal_status") != "ok":
        print("marginal trace status:", data.get("marginal_status"))
    frac = data.get("marginal_pupil_fraction", 1.0)
    fld = osp['fov'].fields[0]  # rayoptics wants a Field instance, not a tuple
    pkg = trace_base(opm, [0.0, frac], fld, data["wavelength_nm"], apply_vignetting=False)
    ray = pkg.ray if hasattr(pkg, 'ray') else pkg[0]
    # Interfaces: [0]=obj, [1..n]=lens, [n+1]=image.
    last_optical = ray[len(data["surfaces"])]
    p, d = last_optical[0], last_optical[1]
    z_focus = p[2] + bfd
    t = (z_focus - p[2]) / d[2]
    y_focus = p[1] + d[1] * t

    ours = {
        "efl": data["efl_d_mm"], "bfl": data["bfl_d_mm"], "fno": data["fno_d"],
        "marginal_y": data["marginal_y_at_focus_mm"],
    }
    theirs = {"efl": efl, "bfl": bfd, "fno": fno, "marginal_y": y_focus}

    print(f"{'quantity':10s} {'OZ':>16s} {'rayoptics':>16s} {'diff':>12s}")
    worst = 0.0
    for k in ours:
        a, b = ours[k], theirs[k]
        if a is None or b is None or (isinstance(a, float) and math.isnan(a)):
            print(f"{k:10s} {'n/a':>16s} {b!r:>16s}")
            continue
        d_ = abs(a - b)
        worst = max(worst, d_)
        print(f"{k:10s} {a:16.8g} {b:16.8g} {d_:12.3e}")

    efl_rel = abs(ours["efl"] - theirs["efl"]) / max(abs(theirs["efl"]), 1e-300)
    ok_efl = efl_rel < 5e-6
    ok_bfl = abs(ours["bfl"] - theirs["bfl"]) < 5e-4
    ok_fno = abs(ours["fno"] - theirs["fno"]) / max(theirs["fno"], 1e-300) < 5e-6
    ok_marg = abs(ours["marginal_y"] - theirs["marginal_y"]) < 0.01  # mm
    print(f"checks: EFL rel {efl_rel:.2e}; BFL abs {abs(ours['bfl']-theirs['bfl']):.2e};"
          f" F# rel {abs(ours['fno']-theirs['fno'])/theirs['fno']:.2e};"
          f" marginal-y abs {abs(ours['marginal_y']-theirs['marginal_y']):.2e} mm")
    ok = ok_efl and ok_bfl and ok_fno and ok_marg
    print("LENS-CROSSCHECK", "OK" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
