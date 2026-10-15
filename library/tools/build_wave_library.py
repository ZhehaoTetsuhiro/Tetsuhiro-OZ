#!/usr/bin/env python3
"""生成 library/wave/ 下的全部波动光学示例（仅标准库）。

    python3 library/tools/build_wave_library.py                       # 全量生成
    python3 library/tools/build_wave_library.py --only diff,beam      # 只生成若干族
    python3 library/tools/build_wave_library.py --out /tmp/lib        # 输出到别处

生成的每个 JSON 都是完整、可直接 POST /api/simulate 的 Config，并内嵌
meta{name, zod{category, difficulty, description, expected, tags, theory}, sample?}。
同时写出 index.json（文件名 → 元数据）。
"""

import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from wavegen.common import Registry, write_all, wave_configs_only, quantum_configs  # noqa: E402

import importlib  # noqa: E402


def _fam(modname):
    """按需导入族模块，缺失时抛 ImportError（便于分阶段生成）。"""
    return importlib.import_module("wavegen." + modname)


FAMILIES = [
    ("diff", ["fam_diff.single_slit", "fam_diff.single_slit_fresnel",
              "fam_diff.circular_aperture", "fam_diff.double_slit_diffraction",
              "fam_diff.multi_slit", "fam_diff.shapes"]),
    ("diff2", ["fam_diff2.gratings_sin_amp", "fam_diff2.gratings_binary_amp",
               "fam_diff2.gratings_phase_sin", "fam_diff2.gratings_phase_binary",
               "fam_diff2.gratings_blazed", "fam_diff2.gratings_rotation",
               "fam_diff2.gratings_missing", "fam_diff2.grating_spectrum",
               "fam_diff2.zone_plates", "fam_diff2.talbot", "fam_diff2.knife_edge",
               "fam_diff2.vortex_farfield", "fam_diff2.zernike_psf",
               "fam_diff2.fourf_filter", "fam_diff2.diffraction_extras"]),
    ("interf", ["fam_interf.young_double_slit", "fam_interf.two_hole_young",
                "fam_interf.multi_beam_interf", "fam_interf.newton_rings",
                "fam_interf.equal_thickness", "fam_interf.michelson_fringes",
                "fam_interf.mach_zehnder_phase", "fam_interf.two_source_fringes",
                "fam_interf.parallel_plate", "fam_interf.visibility_study"]),
    ("beam", ["fam_beam.gauss_waist_sweep", "fam_beam.gauss_z_sweep",
              "fam_beam.gouy_phase", "fam_beam.lg_modes", "fam_beam.hg_modes",
              "fam_beam.bessel_beams", "fam_beam.spherical_waves", "fam_beam.tilted_plane",
              "fam_beam.gauss_lens", "fam_beam.vortex_prop", "fam_beam.wavelength_beam"]),
    ("pol", ["fam_pol.malus_law", "fam_pol.polarizer_pairs", "fam_pol.qwp_sweep",
             "fam_pol.hwp_sweep", "fam_pol.ellipsometry", "fam_pol.rotator_sweep",
             "fam_pol.custom_jones_examples", "fam_pol.polarization_interference",
             "fam_pol.vectorial_focus", "fam_pol.stokes_scenes"]),
    ("med", ["fam_med.uniform_index", "fam_med.absorption_sweep", "fam_med.gain_medium",
             "fam_med.splitstep_convergence", "fam_med.uniaxial_sweep",
             "fam_med.uniaxial_thickness", "fam_med.biaxial_sweep", "fam_med.kerr_effects",
             "fam_med.saturable_sweep", "fam_med.self_focusing"]),
    ("inst", ["fam_inst.michelson_family", "fam_inst.mach_zehnder_family",
              "fam_inst.fabry_perot", "fam_inst.sagnac_ring", "fam_inst.fourf_processing",
              "fam_inst.vortex_coronagraph", "fam_inst.grating_spectrometer",
              "fam_inst.zoneplate_microscope", "fam_inst.telescopes",
              "fam_inst.beam_expander", "fam_inst.spatial_filter",
              "fam_inst.phase_contrast", "fam_inst.confocal", "fam_inst.dic_microscope",
              "fam_inst.metalens_instruments", "fam_inst.imaging_systems"]),
    ("sweep", ["fam_sweep.sweep_double_slit_d", "fam_sweep.sweep_slit_width",
               "fam_sweep.sweep_aperture_diameter", "fam_sweep.sweep_wavelength",
               "fam_sweep.sweep_grating_period", "fam_sweep.sweep_focal_length",
               "fam_sweep.sweep_propagation_distance", "fam_sweep.sweep_slit_count",
               "fam_sweep.sweep_index", "fam_sweep.sweep_polarizer_angle",
               "fam_sweep.sweep_retardance", "fam_sweep.sweep_spiral_charge",
               "fam_sweep.sweep_axicon_alpha", "fam_sweep.sweep_waist",
               "fam_sweep.sweep_diffuser_sigma", "fam_sweep.sweep_zernike",
               "fam_sweep.sweep_absorption"]),
    ("quantum", ["fam_quantum.quantum_all"]),
]


def build(outdir, only=None):
    reg = Registry()
    ran = []
    for name, fns in FAMILIES:
        if only and name not in only:
            continue
        for dotted in fns:
            modname, funcname = dotted.split(".")
            fn = getattr(_fam(modname), funcname)
            fn(reg)
            ran.append(dotted)
    index = write_all(reg, outdir)
    return reg, index, ran


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(os.path.dirname(HERE), "wave"))
    ap.add_argument("--only", default=None, help="逗号分隔的族名（diff,diff2,interf,beam,pol,med,inst,sweep,quantum）")
    ap.add_argument("--json", action="store_true", help="只打印统计")
    a = ap.parse_args()
    only = set(a.only.split(",")) if a.only else None
    reg, index, ran = build(a.out, only)
    counts = reg.counts()
    wave = wave_configs_only(reg)
    quant = quantum_configs(reg)
    sizes = {}
    for f, _m, c in reg.items:
        key = c["grid"]["size"] if "grid" in c else "quantum"
        sizes[key] = sizes.get(key, 0) + 1
    out = {
        "total": len(reg.items),
        "wave_top_level": len(wave),
        "quantum_subdir": len(quant),
        "by_category": counts,
        "grid_sizes": sizes,
        "families_run": len(ran),
    }
    print(json.dumps(out, ensure_ascii=False, indent=2))
    return out


if __name__ == "__main__":
    main()
