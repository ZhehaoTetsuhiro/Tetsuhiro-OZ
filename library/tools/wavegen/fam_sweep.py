"""参数扫掠家族（sweep_*）—— 目标 100+。

每个家族固定几何、只扫一个物理量，逐条给出「该参数下应有的解析结果」，便于教学时
一条条对照；也是回归测试的天然素材（同一公式、不同取值）。
"""

import math

from wavegen.common import *   # noqa: F401,F403
from wavegen.fam_interf import mz_comps, michelson_comps, d_port_z, d_port_nx
from wavegen.fam_pol import _pol, _ret

W20 = 0.02
W10 = 0.01
W8 = 0.008
PI = math.pi
CAT = "参数扫掠"


def _add(reg, fname, cfg, sample=False):
    z = cfg["meta"]["zod"]
    reg.add(fname, CAT, z["difficulty"], z["description"], z["expected"], cfg,
            tags=z["tags"], theory=z.get("theory"), sample=sample)


def _mk(name, comps, srcs, size, width, wl, desc, expect, tags, theory,
        difficulty="进阶", polarized=False, method="asm", sample=False):
    return build_config(name=name, size=size, width=width, wl=wl, components=comps,
                        sources=srcs, polarized=polarized, method=method,
                        category=CAT, difficulty=difficulty, description=desc,
                        expected=expect, tags=tags, theory=theory, sample=sample)


# ==========================================================================
# 1. 双缝缝距扫掠：Δx = λL/d
# ==========================================================================

def sweep_double_slit_d(reg):
    z, a, wl, W = 1.0, 0.15e-3, WL["633"], W20
    N = pick_grid(a, W, minimum=4.0)
    for i, d in enumerate([0.20e-3, 0.30e-3, 0.40e-3, 0.50e-3, 0.60e-3, 0.80e-3,
                           1.00e-3, 1.25e-3, 1.50e-3, 2.00e-3, 2.50e-3, 3.00e-3]):
        fring = fringe_spacing(wl, z, d)
        check_sampling(W / N, {"条纹间距": fring}, 4.0)
        cfg = _mk("扫掠·双缝缝距 d=%.2f mm → 条纹间距 %.3f mm" % (d * 1e3, fring * 1e3),
                  [aperture("ds", (0, 0, 0), double_slit(d, a, W * 0.95),
                            label="双缝 d=%.2f mm a=0.15 mm" % (d * 1e3)),
                   det_rect("det", z, W * 0.95, W * 0.95, label="屏（1 m）")],
                  [plane_source(wl=wl)], N, W, wl,
                  "固定缝宽与屏距，把双缝间距从 0.20 mm 扫到 3.00 mm（本条第 %d 条）。" % (i + 1),
                  "条纹间距 Δx=λL/d=%.4f mm，窗口内约占 %.1f 条/mm；"
                  "d 加倍则条纹密度加倍（严格反比）。" % (fring * 1e3, 1.0 / (fring * 1e3)),
                  ["缝距扫掠", "条纹密度"], "Δx = λL/d", "入门", sample=(i < 2))
        cfg["meta"]["zod"]["sweep"] = {"param": "d", "value_mm": d * 1e3, "index": i + 1, "of": 12}
        _add(reg, "sweep_double_slit_d_%02d.json" % (i + 1), cfg, sample=(i < 2))


# ==========================================================================
# 2. 单缝宽度扫掠：x1 = λL/a
# ==========================================================================

def sweep_slit_width(reg):
    z, wl, W = 1.5, WL["633"], W20
    for i, a in enumerate([0.10e-3, 0.15e-3, 0.20e-3, 0.30e-3, 0.40e-3,
                           0.50e-3, 0.60e-3, 0.80e-3, 1.00e-3, 1.20e-3]):
        N = pick_grid(a, W, minimum=3.0)
        x1 = sinc_first_zero(wl, z, a)
        check_sampling(W / N, {"缝宽": a}, 3.0)
        cfg = _mk("扫掠·单缝宽度 a=%.2f mm → 第一暗点 %.2f mm" % (a * 1e3, x1 * 1e3),
                  [aperture("slit", (0, 0, 0), rectangle(a, W * 0.95),
                            label="单缝 a=%.2f mm" % (a * 1e3)),
                   det_rect("det", z, W * 0.95, W * 0.95, label="远场屏")],
                  [plane_source(wl=wl)], N, W, wl,
                  "固定波长与屏距 1.5 m，把单缝宽度从 0.10 mm 扫到 1.20 mm（第 %d 条）。" % (i + 1),
                  "中央亮纹宽度 2x₁=2λL/a=%.3f mm（反比于缝宽）；缝宽加倍 → 图样减半。",
                  ["缝宽扫掠", "衍射展宽"], "x₁ = λL/a", "入门", sample=(i < 2))
        cfg["meta"]["zod"]["sweep"] = {"param": "a", "value_mm": a * 1e3, "index": i + 1, "of": 10}
        _add(reg, "sweep_slit_width_%02d.json" % (i + 1), cfg, sample=(i < 2))


# ==========================================================================
# 3. 圆孔直径扫掠：艾里半径 1.22λf/D
# ==========================================================================

def sweep_aperture_diameter(reg):
    f, wl, W = 0.3, WL["633"], W8
    for i, D in enumerate([0.4e-3, 0.6e-3, 0.8e-3, 1.0e-3, 1.5e-3,
                           2.0e-3, 2.5e-3, 3.0e-3, 4.0e-3, 5.0e-3]):
        r1 = airy_radius(wl, f, D)
        N = pick_grid(r1, W, minimum=4.0)
        check_sampling(W / N, {"口径 D": D}, 4.0)
        cfg = _mk("扫掠·圆孔直径 D=%.1f mm → 艾里半径 %.1f µm" % (D * 1e3, r1 * 1e6),
                  [aperture("ap", (0, 0, 0), circle(D / 2), label="圆孔 D=%.1f mm" % (D * 1e3)),
                   lens("fl", (0, 0, 0.02), f=f, aperture=W * 0.49, label="傅里叶透镜 f=0.3 m"),
                   det_circle("det", 0.32, W * 0.49, label="焦面")],
                  [plane_source(wl=wl)], N, W, wl,
                  "同焦距 0.3 m，把圆孔直径从 0.4 mm 扫到 5.0 mm（第 %d 条）。" % (i + 1),
                  "艾里第一暗环 r₁=1.22λf/D=%.1f µm，反比于口径；直径加倍 → 焦斑减半、"
                  "峰值强度变为 4 倍（能量集中）。" % (r1 * 1e6),
                  ["口径扫掠", "艾里斑", "分辨率"], "r₁ = 1.22λf/D", "入门", sample=(i < 2))
        cfg["meta"]["zod"]["sweep"] = {"param": "D", "value_mm": D * 1e3, "index": i + 1, "of": 10}
        _add(reg, "sweep_aperture_d_%02d.json" % (i + 1), cfg, sample=(i < 2))


# ==========================================================================
# 4. 波长扫掠（单缝）
# ==========================================================================

def sweep_wavelength(reg):
    z, a, W = 1.5, 0.25e-3, W20
    N = pick_grid(a, W, minimum=3.0)
    keys = ["405", "450", "473", "500", "532", "550", "589", "633", "700", "780"]
    for i, key in enumerate(keys):
        wl = WL[key]
        x1 = sinc_first_zero(wl, z, a)
        cfg = _mk("扫掠·波长 λ=%s nm → 第一暗点 %.2f mm" % (key, x1 * 1e3),
                  [aperture("slit", (0, 0, 0), rectangle(a, W * 0.95), label="单缝 a=0.25 mm"),
                   det_rect("det", z, W * 0.95, W * 0.95, label="远场屏")],
                  [plane_source(wl=wl, label="λ=%s nm" % key)], N, W, wl,
                  "固定缝宽与屏距，把波长从 405 nm 扫到 780 nm（第 %d 条）。" % (i + 1),
                  "x₁=λL/a=%.3f mm 与波长严格成正比；本条的中央峰宽度相当于 405 nm 的 %.2f 倍。"
                  % (x1 * 1e3, wl / WL["405"]),
                  ["波长扫掠", "色散"], "x₁ = λL/a ∝ λ", "入门", sample=(i < 2))
        cfg["meta"]["zod"]["sweep"] = {"param": "lambda_nm", "value": key, "index": i + 1, "of": len(keys)}
        _add(reg, "sweep_wavelength_%02d.json" % (i + 1), cfg, sample=(i < 2))


# ==========================================================================
# 5. 光栅周期扫掠：x_q = qfλ/Λ
# ==========================================================================

def sweep_grating_period(reg):
    f, wl, W, N = 0.3, WL["633"], W8, 1024
    for i, Lam in enumerate([60e-6, 80e-6, 100e-6, 120e-6, 150e-6,
                             180e-6, 200e-6, 250e-6, 300e-6, 400e-6]):
        x1 = grating_order_pos(1, f, wl, Lam)
        check_sampling(W / N, {"周期": Lam}, 5.0)
        check_sampling(W / N, {"±1 级": x1}, 15.0)
        n_ord = int(W * 0.45 / x1)
        cfg = _mk("扫掠·光栅周期 Λ=%.0f µm → ±1 级在 %.2f mm" % (Lam * 1e6, x1 * 1e3),
                  [comp("g", "grating", (0, 0, 0),
                        params={"kind": "phase_sin", "period": Lam, "modulation": 1.6},
                        label="相位光栅 Λ=%.0f µm" % (Lam * 1e6)),
                   lens("fl", (0, 0, 0.02), f=f, aperture=W * 0.49, label="傅里叶透镜 f=0.3 m"),
                   det_rect("det", 0.32, W * 0.95, W * 0.95, label="焦面")],
                  [plane_source(wl=wl)], N, W, wl,
                  "同焦距 0.3 m，把光栅周期从 60 µm 扫到 400 µm（第 %d 条）。" % (i + 1),
                  "级次位置 x₁=fλ/Λ=%.3f mm，反比于周期；窗口内可见约 ±%d 级。"
                  "周期越小色散越大（但各级越稀疏）。" % (x1 * 1e3, n_ord),
                  ["光栅周期", "级次位置", "色散"], "x_q = q f λ/Λ", "进阶", sample=(i < 2))
        cfg["meta"]["zod"]["sweep"] = {"param": "period_um", "value": Lam * 1e6, "index": i + 1, "of": 10}
        _add(reg, "sweep_grating_period_%02d.json" % (i + 1), cfg, sample=(i < 2))


# ==========================================================================
# 6. 焦距扫掠：焦斑与艾里半径
# ==========================================================================

def sweep_focal_length(reg):
    wl, W, R = WL["633"], W10, 2.0e-3
    for i, f in enumerate([0.05, 0.08, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.75, 1.0]):
        r1 = airy_radius(wl, f, 2 * R)
        N = W = None
        for Wtry in (0.004, 0.006, 0.01):
            if Wtry / 2.0 * 0.95 < R:
                continue
            try:
                N = pick_grid(r1, Wtry, minimum=3.0)
                W = Wtry
                break
            except AssertionError:
                continue
        assert N, "no grid fits f=%.3f" % f
        check_sampling(W / N, {"光瞳半径": R}, 20.0)
        cfg = _mk("扫掠·焦距 f=%.0f mm → 艾里半径 %.1f µm" % (f * 1e3, r1 * 1e6),
                  [aperture("ap", (0, 0, 0), circle(R), label="光瞳 D=4 mm"),
                   lens("l", (0, 0, 0.01), f=f, aperture=W * 0.49, label="透镜 f=%.0f mm" % (f * 1e3)),
                   detector("det", 0.01 + f, label="焦面",
                            params={"strehl_aperture": R, "strehl_distance": f})],
                  [plane_source(wl=wl)], N, W, wl,
                  "固定光瞳直径 4 mm（F 数随焦距变化），把焦距从 50 mm 扫到 1 m（第 %d 条）。" % (i + 1),
                  "艾里半径 r₁=1.22λf/D=%.1f µm ∝ f；F 数 F#=f/D=%.1f，焦深 ≈±2λF#²=%.0f µm。"
                  "长焦系统的焦斑大但焦深长——这是「分辨率 vs 焦深」的基本取舍。"
                  % (r1 * 1e6, f / (2 * R), 2 * wl * (f / (2 * R)) ** 2 * 1e6),
                  ["焦距扫掠", "艾里斑", "焦深", "F数"], "r₁ = 1.22λF#；焦深 ∝ F#²",
                  "入门", sample=(i < 2))
        cfg["meta"]["zod"]["sweep"] = {"param": "f_mm", "value": f * 1e3, "index": i + 1, "of": 10}
        _add(reg, "sweep_focal_length_%02d.json" % (i + 1), cfg, sample=(i < 2))


# ==========================================================================
# 7. 传播距离扫掠：菲涅耳数
# ==========================================================================

def sweep_propagation_distance(reg):
    wl, a, W = WL["633"], 0.3e-3, W10
    N = pick_grid(a, W, minimum=3.0)
    for i, z in enumerate([5e-3, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 4.0]):
        NF = fresnel_number(a, wl, z)
        scale = math.sqrt(wl * z) if z > 0 else 1e-9
        check_sampling(W / N, {"缝宽": a}, 3.0)
        cfg = _mk("扫掠·传播距离 z=%.0f mm → 菲涅耳数 N_F=%.1f" % (z * 1e3, NF),
                  [aperture("slit", (0, 0, 0), rectangle(a, W * 0.95), label="单缝 a=0.3 mm"),
                   det_rect("det", z, W * 0.95, W * 0.95, label="屏（z=%.0f mm）" % (z * 1e3))],
                  [plane_source(wl=wl)], N, W, wl,
                  "固定 0.3 mm 单缝，把观察距离从 5 mm 扫到 4 m（第 %d 条）。" % (i + 1),
                  "菲涅耳数 N_F=a²/(λz)=%.2f：N_F≫1 是近场（几何投影 + 边缘振荡，特征尺度 √(λz)=%.0f µm），"
                  "N_F≪1 是远场（sinc² 图样，第一暗点 λz/a=%.2f mm）。本条位于%s。"
                  % (NF, scale * 1e6, sinc_first_zero(wl, z, a) * 1e3,
                     "近场" if NF > 1 else "远场"),
                  ["距离扫掠", "菲涅耳数", "近场到远场"], "N_F = a²/(λz)",
                  "进阶", sample=(i < 3))
        cfg["meta"]["zod"]["sweep"] = {"param": "z_mm", "value": z * 1e3, "index": i + 1, "of": 10}
        _add(reg, "sweep_distance_%02d.json" % (i + 1), cfg, sample=(i < 3))


# ==========================================================================
# 8. 缝数扫掠：主极大锐度 ∝ 1/N
# ==========================================================================

def sweep_slit_count(reg):
    d, a, wl, z, W = 1.0e-3, 0.15e-3, WL["633"], 1.0, W20
    N = pick_grid(a, W, minimum=3.0)
    for i, n in enumerate([2, 3, 4, 5, 6, 8, 10, 16]):
        fring = fringe_spacing(wl, z, d)
        comps = []
        for m in range(n):
            comps.append(aperture("s%d" % m, (0, 0, 0), rectangle(a, W * 0.95), parallel=True,
                                  x=(m - (n - 1) / 2.0) * d))
        comps.append(det_rect("det", z, W * 0.95, W * 0.95, label="屏（1 m）"))
        cfg = _mk("扫掠·缝数 N=%d → 主极大宽度 ∝ 1/N" % n, comps,
                  [plane_source(wl=wl)], N, W, wl,
                  "固定缝距 1.0 mm、缝宽 0.15 mm，把缝数从 2 增到 16（第 %d 条）。" % (i + 1),
                  "主极大位置仍是 Δx=λL/d=%.2f mm，但角宽度压缩到 ∝1/N=1/%d（%s），"
                  "次极大数 N−2=%d、次极大强度上限约 1/N²=%.4f。缝数越多「谱线」越锐利——"
                  "这就是光栅分辨率 ∝ 刻线数 N 的来源。" % (fring * 1e3, n,
                                                            "比双缝锐 %.0f 倍" % (n / 2.0) if n > 2 else "与双缝同",
                                                            max(0, n - 2), 1.0 / n ** 2),
                  ["缝数扫掠", "主极大锐度", "光栅分辨率"], "主极大宽度 ∝ 1/(N d)",
                  "进阶", sample=(i < 2))
        cfg["meta"]["zod"]["sweep"] = {"param": "N", "value": n, "index": i + 1, "of": 8}
        _add(reg, "sweep_slit_count_%02d.json" % (i + 1), cfg, sample=(i < 2))


# ==========================================================================
# 9. 折射率扫掠（干涉读出）
# ==========================================================================

def sweep_index(reg):
    wl, W, d = WL["633"], W10, 0.01
    for i, dn in enumerate([0.0, 5e-5, 1e-4, 1.5e-4, 2e-4, 2.5e-4, 3e-4, 4e-4]):
        n = 1.0 + dn
        ph = k_of(wl) * dn * d
        cfg = _mk("扫掠·折射率 n=%.6f → 相位 %.2f×2π" % (n, ph / (2 * PI)),
                  mz_comps(W, A=0.08,
                           arm_extra=[comp("plate", "medium", (0, 0, 0.04), shape=circle(W * 0.4),
                                           params={"distance": d, "index": n, "absorption": 0.0,
                                                   "steps": 4},
                                           label="介质板 d=10 mm、n=%.6f" % n)],
                           d1_label="端口 1", d2_label="端口 2"),
                  [plane_source(wl=wl)], 512, W, wl,
                  "固定 10 mm 介质厚度，把折射率从 1.00000 扫到 1.00040（第 %d 条）。" % (i + 1),
                  "相位 φ=k(n−1)d=%.3f rad=%.3f×2π，端口功率 = cos²(φ/2)=%.3f。"
                  "Δn 的周期为 λ/d=%.0f ppm=%.4f —— 用干涉法测折射率的分辨率可达 1e-7 量级。"
                  % (ph, ph / (2 * PI), math.cos(ph / 2) ** 2, wl / d * 1e6, wl / d),
                  ["折射率扫掠", "干涉计量", "灵敏度"], "φ = k(n−1)d，周期 Δn = λ/d",
                  "高级", sample=(i < 2))
        cfg["meta"]["zod"]["sweep"] = {"param": "n", "value": n, "index": i + 1, "of": 8}
        _add(reg, "sweep_index_%02d.json" % (i + 1), cfg, sample=(i < 2))


# ==========================================================================
# 10. 偏振片角度扫掠（马吕斯）
# ==========================================================================

def sweep_polarizer_angle(reg):
    wl, W = WL["633"], W10
    for i, deg in enumerate([0, 5, 10, 15, 20, 30, 40, 45, 50, 60, 70, 80, 90]):
        clip = W * 0.6
        comps = [aperture("clip", (0, 0, 0), rectangle(clip, clip), label="限光阑"),
                 _pol("p1", 0.02, 0.0, "起偏器 0°"),
                 _pol("p2", 0.04, math.radians(deg), "检偏器 %.0f°" % deg),
                 det_rect("det", 0.06, clip, clip, label="读数")]
        cfg = _mk("扫掠·检偏角 %.0f° → I/I₀=%.4f" % (deg, math.cos(math.radians(deg)) ** 2),
                  comps, [source("src", "plane", {"power": 1e-3, "polarization": "x"},
                                 (0.0, 0.0, SRC_Z), wavelength=wl, label="x 偏振平面波")],
                  512, W, wl,
                  "固定入射偏振方向，把检偏器从 0° 转到 90°（第 %d 条）。" % (i + 1),
                  "马吕斯定律 I/I₀=cos²θ=%.4f；θ=45° 时恰好一半，90° 完全消光。"
                  "扫掠给出完整的一条 cos² 曲线（每 5° 一点）。" % (math.cos(math.radians(deg)) ** 2),
                  ["偏振片角度扫掠", "马吕斯定律", "cos²"], "I = I₀cos²θ",
                  "入门", polarized=True, sample=(i < 3))
        cfg["meta"]["zod"]["sweep"] = {"param": "angle_deg", "value": deg, "index": i + 1, "of": 13}
        _add(reg, "sweep_polarizer_angle_%02d.json" % (i + 1), cfg, sample=(i < 3))


# ==========================================================================
# 11. 波片延迟扫掠
# ==========================================================================

def sweep_retardance(reg):
    wl, W = WL["633"], W10
    for i, frac in enumerate([0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 1.0]):
        ret = frac * 2 * PI
        clip = W * 0.6
        comps = [aperture("clip", (0, 0, 0), rectangle(clip, clip), label="限光阑"),
                 _pol("p0", 0.02, 0.0, "起偏器 0°"),
                 _ret("r", 0.04, ret, PI / 4, "波片（快轴 45°、延迟 %.3f×2π）" % frac),
                 _pol("a", 0.06, PI / 2, "检偏器 90°（正交）"),
                 det_rect("det", 0.08, clip, clip, label="读数")]
        I = math.sin(ret / 2) ** 2
        cfg = _mk("扫掠·波片延迟 %.3f×2π → 正交透过率 %.3f" % (frac, I),
                  comps, [source("src", "plane", {"power": 1e-3, "polarization": "x"},
                                 (0.0, 0.0, SRC_Z), wavelength=wl, label="x 偏振平面波")],
                  512, W, wl,
                  "固定起偏 0°、快轴 45°、检偏 90°（正交），把波片相位延迟从 0 扫到 2π（第 %d 条）。" % (i + 1),
                  "正交检偏下透过率 I=sin²(δ/2)=%.3f：δ=0 与 2π 时全暗（等价于无波片）、"
                  "δ=π 时全亮（半波片把偏振转 90°）、δ=π/2 时 50%%（四分之一波片产生圆偏振）。"
                  "这条曲线就是判断波片规格的标准测量法（旋转正交偏振片法）。" % I,
                  ["波片延迟扫掠", "正交偏振", "波片定标"], "I = sin²(δ/2)",
                  "进阶", polarized=True, sample=(i < 3))
        cfg["meta"]["zod"]["sweep"] = {"param": "retardance_2pi", "value": frac, "index": i + 1, "of": 8}
        _add(reg, "sweep_retardance_%02d.json" % (i + 1), cfg, sample=(i < 3))


# ==========================================================================
# 12. 螺旋相位拓扑荷扫掠
# ==========================================================================

def sweep_spiral_charge(reg):
    w0, wl, W = 1.2e-3, WL["633"], W8
    N = 512
    for i, l in enumerate([-5, -3, -1, 0, 1, 2, 3, 5]):
        cfg = _mk("扫掠·拓扑荷 l=%d" % l,
                  [comp("spp", "spiral_phase", (0, 0, 0), shape=circle(W * 0.45),
                        params={"charge": l}, label="螺旋相位板 l=%d" % l),
                   detector("det", 0.3, label="传播 0.3 m 后")],
                  [gauss_source(w0, wl=wl, label="高斯 w₀=1.2 mm")], N, W, wl,
                  "把螺旋相位板的拓扑荷从 −5 扫到 +5（第 %d 条）：l 是相位绕轴一周的 2π 圈数。" % (i + 1),
                  ("l=0：普通高斯（轴上最亮）。" if l == 0 else
                   "l=%d：轴上强度为零（相位奇点），环状光斑半径随 |l| 增大；"
                   "相位图绕轴 %d 圈 2π。l 的符号只改变螺旋旋向（±OAM），强度分布相同。"
                   % (l, abs(l))),
                  ["拓扑荷扫掠", "光学涡旋", "OAM"], "相位 exp(ilθ)；轴上 ∝ r^{2|l|}",
                  "进阶", sample=(i < 3))
        cfg["meta"]["zod"]["sweep"] = {"param": "l", "value": l, "index": i + 1, "of": 8}
        _add(reg, "sweep_spiral_charge_%02d.json" % (i + 1), cfg, sample=(i < 3))


# ==========================================================================
# 13. 轴锥镜锥角扫掠
# ==========================================================================

def sweep_axicon_alpha(reg):
    w0, wl = 1.5e-3, WL["633"]
    for i, alpha in enumerate([0.005, 0.008, 0.01, 0.015, 0.02, 0.03, 0.05, 0.08]):
        n_ax = 1.5
        core = bessel_core_radius(wl, (n_ax - 1) * alpha)
        if core / (W8 / 1024.0) >= 3.0:
            N, W = 1024, W8
        elif core / (0.004 / 1024.0) >= 3.0:
            N, W = 1024, 0.004
        else:
            N, W = 2048, 0.004
        zmax = w0 / math.tan((n_ax - 1) * alpha)
        cfg = _mk("扫掠·轴锥镜 α=%.1f mrad → 中心斑 %.1f µm" % (alpha * 1e3, core * 1e6),
                  [comp("ax", "axicon", (0, 0, 0), shape=circle(W * 0.45),
                        params={"alpha": alpha, "index": n_ax}, label="轴锥镜 α=%.1f mrad" % (alpha * 1e3)),
                   detector("det", min(0.15, 0.8 * zmax), label="贝塞尔区")],
                  [gauss_source(w0, wl=wl, label="高斯 w₀=1.5 mm")], N, W, wl,
                  "把轴锥镜的锥角从 5 mrad 扫到 80 mrad（第 %d 条）：锥角决定贝塞尔光束的「锥面波矢」。"
                  % (i + 1),
                  "中心亮斑半径 2.405/(k sinβ)=%.1f µm（反比于锥角 β≈(n−1)α）；"
                  "无衍射区长度 z_max=w₀/tanβ=%.0f mm（也反比于锥角）。"
                  "锥角越大 → 光斑越小、无衍射区越短——这就是「亚波长光斑」的另一条路线。"
                  % (core * 1e6, zmax * 1e3),
                  ["轴锥镜扫掠", "贝塞尔光束", "无衍射区"], "核心半径 ≈2.405/(k sinβ)，z_max=w₀/tanβ",
                  "进阶", sample=(i < 3))
        cfg["meta"]["zod"]["sweep"] = {"param": "alpha_mrad", "value": alpha * 1e3, "index": i + 1, "of": 8}
        _add(reg, "sweep_axicon_alpha_%02d.json" % (i + 1), cfg, sample=(i < 3))


# ==========================================================================
# 14. 高斯束腰扫掠（同一距离的展宽）
# ==========================================================================

def sweep_waist(reg):
    z, wl = 0.5, WL["633"]
    for i, w0 in enumerate([0.05e-3, 0.1e-3, 0.15e-3, 0.25e-3, 0.4e-3, 0.6e-3, 1.0e-3, 1.5e-3]):
        wz = waist_at(wl, w0, z)
        W = pick_window(2.2 * wz)
        N = pick_grid(w0, W, minimum=4.0)
        check_sampling(W / N, {"束腰": w0}, 4.0)
        cfg = _mk("扫掠·束腰 w₀=%.2f mm → 0.5 m 后 %.2f mm" % (w0 * 1e3, wz * 1e3),
                  [detector("det", z, label="z=0.5 m 处")],
                  [gauss_source(w0, wl=wl, label="高斯 w₀=%.2f mm" % (w0 * 1e3))], N, W, wl,
                  "把高斯光束的束腰从 0.05 mm 扫到 1.5 mm（第 %d 条），都在 0.5 m 处观察。" % (i + 1),
                  "w(0.5)=%.3f mm；瑞利长度 z_R=πw₀²/λ=%.3f m，发散半角 θ=λ/(πw₀)=%.2f mrad。"
                  "束腰小 → 展宽剧烈（本条的 w/w₀=%.2f）。"
                  % (wz * 1e3, rayleigh_range(wl, w0), gauss_divergence(wl, w0) * 1e3, wz / w0),
                  ["束腰扫掠", "发散角", "瑞利长度"], "z_R=πw₀²/λ，θ=λ/(πw₀)",
                  "入门", sample=(i < 3))
        cfg["meta"]["zod"]["sweep"] = {"param": "w0_um", "value": w0 * 1e6, "index": i + 1, "of": 8}
        _add(reg, "sweep_waist_%02d.json" % (i + 1), cfg, sample=(i < 3))


# ==========================================================================
# 15. 漫射体强度扫掠（散斑对比度）
# ==========================================================================

def sweep_diffuser_sigma(reg):
    wl, W, z = WL["633"], W8, 0.2
    N = 512
    for i, sigma in enumerate([0.2, 0.5, 1.0, 2.0, 3.1416, 4.5]):
        cfg = _mk("扫掠·漫射体相位标准差 σ=%.2f rad" % sigma,
                  [comp("d", "diffuser", (0, 0, 0), shape=circle(W * 0.45),
                        params={"sigma": sigma, "correlation": 20e-6, "seed": 7},
                        label="漫射体 σ=%.2f rad、相关长度 20 µm" % sigma),
                   det_circle("det", z, W * 0.45, label="散斑场")],
                  [plane_source(wl=wl)], N, W, wl,
                  "把漫射体的相位标准差从 0.2 rad 扫到 4.5 rad（第 %d 条）：σ 是相位屏的「粗糙程度」。"
                  % (i + 1),
                  "散斑对比度 V=σ_I/⟨I⟩ 随 σ 增大从接近 0 单调升到 1（σ≳π 时饱和为完全散斑）；"
                  "散斑颗粒尺度 ≈λz/D=%.0f µm（D 为照明口径 %.0f mm）。"
                  "σ=π 是本内核默认的「强漫射」值（测试实测对比度≈1）。"
                  % (speckle_size(wl, z, W * 0.9) * 1e6, W * 0.9 * 1e3),
                  ["散斑扫掠", "漫射体", "对比度"], "σ_I/⟨I⟩ → 1（强漫射）；散斑尺度 λz/D",
                  "进阶", sample=(i < 2))
        cfg["meta"]["zod"]["sweep"] = {"param": "sigma_rad", "value": sigma, "index": i + 1, "of": 6}
        _add(reg, "sweep_diffuser_%02d.json" % (i + 1), cfg, sample=(i < 2))


# ==========================================================================
# 16. 泽尼克系数扫掠（各类像差）
# ==========================================================================

def sweep_zernike(reg):
    wl, f, W, R = WL["633"], 0.3, 0.006, 2.0e-3
    r1 = airy_radius(wl, f, 2 * R)
    N = pick_grid(r1, W, minimum=4.0)
    specs = [("c4", "离焦 Z4"), ("c5", "像散 Z5"), ("c7", "彗差 Z7"),
             ("c11", "球差 Z11"), ("c10", "三叶 Z10"), ("c3", "倾斜 Z3"),
             ("c9", "二级像散 Z9"), ("c21", "高阶球差 Z21")]
    k_idx = 0
    for i, (key, lbl) in enumerate(specs):
        for j, val in enumerate([0.1, 0.25]):
            k_idx += 1
            params = {"radius": R, key: val}
            cfg = _mk("扫掠·%s = %.2f λ" % (lbl, val),
                      [lens("l", (0, 0, 0), f=f, aperture=R, label="理想透镜 f=0.3 m"),
                       comp("ab", "zernike", (0, 0, 0.02), shape=circle(R), params=params,
                            label="%s = %.2f λ" % (lbl, val)),
                       detector("det", f, label="焦面",
                                params={"strehl_aperture": R, "strehl_distance": f})],
                      [plane_source(wl=wl)], N, W, wl,
                      "在理想透镜后加入 %s=%.2f λ 的泽尼克像差（Noll 序、RMS 归一化）。" % (lbl, val),
                      "波前 RMS = %.2f λ → 理论 Strehl≈exp(−(2πσ)²)=%.3f；"
                      "焦面 PSF 按该像差的特征变形（离焦→同心环、像散→十字/焦线、彗差→彗尾、"
                      "球差→中心亮环）。无像差时艾里半径 r₁=%.1f µm。"
                      % (val, math.exp(-(2 * PI * val) ** 2), r1 * 1e6),
                      ["泽尼克扫掠", "像差", "Strehl"], "Strehl ≈ exp(−(2πσ)²)",
                      "高级", sample=(i == 0 and j == 0))
            cfg["meta"]["zod"]["sweep"] = {"param": key, "value_lambda": val,
                                           "index": len(reg.items) + 1, "of": len(specs) * 2}
            _add(reg, "sweep_zernike_%02d_%d.json" % (i + 1, j + 1), cfg,
                 sample=(i == 0 and j == 0))


# ==========================================================================
# 17. 吸收系数扫掠
# ==========================================================================

def sweep_absorption(reg):
    wl, W, d = WL["633"], W10, 0.01
    for i, T in enumerate([0.99, 0.95, 0.9, 0.8, 0.7, 0.5, 0.3, 0.1]):
        np_ = _abs_n_local(T, d, wl)
        alpha = 2 * k_of(wl) * np_
        cfg = _mk("扫掠·吸收介质透过率 %.0f%% → n''=%.3e" % (T * 100, np_),
                  [comp("m", "medium", (0, 0, 0), shape=circle(W * 0.49),
                        params={"distance": d, "index": 1.5, "absorption": np_, "steps": 10},
                        label="吸收介质 d=10 mm"),
                   detector("det", 0.1, shape=circle(W * 0.49), label="透过功率")],
                  [plane_source(wl=wl)], 512, W, wl,
                  "固定 10 mm 厚度，把目标透过率从 99%% 扫到 10%%（第 %d 条）。" % (i + 1),
                  "复折射率虚部 n''=%.4e（α=2kn''=%.0f /m）；强度透过 T=exp(−αd)=%.3f。"
                  "强度与场的关系是「2k n''」而不是「k n''」——这是复折射率的定义 "
                  "n=n′+iκ、E∝exp(iknd) 的直接结果。" % (np_, alpha, T),
                  ["吸收扫掠", "复折射率", "比尔-朗伯"], "T = exp(−2 k n'' d)",
                  "入门", sample=(i < 2))
        cfg["meta"]["zod"]["sweep"] = {"param": "T", "value": T, "index": i + 1, "of": 8}
        _add(reg, "sweep_absorption_%02d.json" % (i + 1), cfg, sample=(i < 2))


def _abs_n_local(T, d, wl):
    return -math.log(T) / (2 * k_of(wl) * d)
