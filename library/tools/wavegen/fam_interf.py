"""干涉类示例（interf_*）—— 目标 60+。

杨氏双缝（缝距/波长/倾角扫掠）、双孔、多光束、**牛顿环**（用曲率镜实现 r_m²=mλR 的
环形等厚干涉）、等厚（楔形）条纹、迈克尔逊（等臂/不等臂/微倾）、马赫-曾德尔（相位板扫描）、
双光源相干条纹、平行平板/菲索条纹、可见度研究。
"""

import math

from wavegen.common import *   # noqa: F401,F403

W20 = 0.02
W10 = 0.01
PI = math.pi


def _slit(z=1.0, d=1.0e-3, a=0.2e-3, W=W20, label=None):
    return aperture("ds", (0, 0, 0), double_slit(d, a, W * 0.95),
                    label=label or "双缝 d=%.2f mm a=%.2f mm" % (d * 1e3, a * 1e3))


# —— 干涉仪几何助手（严格按场景模型的射线规则摆放，避免「镜子朝向错」的静默失败）——
#   * 沿 +z 前进的光遇到法线 n=(0,0,1)（yaw=0）的镜面 → 回射 -z；
#   * 沿 +x 前进的光遇到 yaw=π/2 的镜面（n=(1,0,0)）→ 回射 -x；
#   * 45° 折转镜：yaw=-π/4 把 +z 变成 +x；yaw=3π/4 把 +x 变成 +z；
#   * 探测器法线必须**迎着**来光：来光 -z → yaw=0；来光 -x → yaw=π/2；
#     来光 +z → yaw=π；来光 +x → yaw=-π/2。

def m_retro_z(cid, z, **kw):
    """回射 +z 来光的平面镜（法线 +z）或球面镜。"""
    yaw = kw.pop("yaw", 0.0)
    return mirror(cid, (0, 0, z), yaw=yaw, **kw)


def m_retro_x(cid, x, **kw):
    """回射 +x 来光的平面镜/球面镜（法线 +x）。"""
    yaw = kw.pop("yaw", PI / 2)
    return mirror(cid, (x, 0, 0), yaw=yaw, **kw)


def d_port_z(cid, z, label, r):
    """记录沿 -z 来光的端口（法线 +z，位置 (0,0,z)）。"""
    return detector(cid, z=z, shape=circle(r), yaw=0.0, label=label)


def d_port_nx(cid, x, label, r):
    """记录沿 -x 来光的端口（法线 +x，位置 (x,0,0)）。"""
    return detector(cid, z=0.0, x=x, shape=circle(r), yaw=PI / 2, label=label)


def michelson_comps(W, L=0.05, r=0.35, bs_R=0.5, m1=None, m2=None,
                    d1_label="回光源端口", d2_label="侧端口", both_ports=True,
                    extra=()):
    """标准迈克尔逊：分束器 + 两个回射端镜 + 两个输出端口。"""
    W35 = W * r
    comps = [beam_splitter("bs", (0, 0, 0), reflectivity=bs_R, aperture=W35,
                           label="分束器 %.0f:%.0f" % (bs_R * 100, (1 - bs_R) * 100))]
    comps.extend(extra)
    m1 = m1 or {}
    m2 = m2 or {}
    comps.append(mirror("m1", (0, 0, L), reflectivity=m1.pop("reflectivity", 1.0),
                        curvature=m1.pop("curvature", 0.0), tilt_x=m1.pop("tilt_x", 0.0),
                        tilt_y=m1.pop("tilt_y", 0.0), yaw=m1.pop("yaw", 0.0),
                        shape=circle(W35), label=m1.pop("label", "端镜 A")))
    comps.append(mirror("m2", (L, 0, 0), reflectivity=m2.pop("reflectivity", 1.0),
                        curvature=m2.pop("curvature", 0.0), tilt_x=m2.pop("tilt_x", 0.0),
                        tilt_y=m2.pop("tilt_y", 0.0), yaw=m2.pop("yaw", PI / 2),
                        shape=circle(W35), label=m2.pop("label", "端镜 B")))
    comps.append(d_port_z("d1", -0.05, d1_label, W35))
    if both_ports:
        comps.append(d_port_nx("d2", -L, d2_label, W35))
    return comps


def mz_comps(W, A=0.08, r=0.35, bs_R=0.5, arm_extra=(), arm2_extra=(),
             d1_label="端口 1", d2_label="端口 2"):
    """标准马赫-曾德尔：两个 45° 折转镜 + 两个分束器 + 两个端口。"""
    W35 = W * r
    comps = [beam_splitter("bs1", (0, 0, 0), reflectivity=bs_R, aperture=W35,
                           label="分束器 1（%.0f:%.0f）" % (bs_R * 100, (1 - bs_R) * 100))]
    comps.extend(arm_extra)                      # 透射臂（+z 先走的那条）
    comps.append(mirror("m1", (0, 0, A), reflectivity=1.0, yaw=-PI / 4,
                        shape=circle(W35), label="折转镜 A（+z→+x）"))
    comps.append(mirror("m2", (A, 0, 0), reflectivity=1.0, yaw=3 * PI / 4,
                        shape=circle(W35), label="折转镜 B（+x→+z）"))
    comps.extend(arm2_extra)                     # 反射臂的其余元件
    comps.append(beam_splitter("bs2", (A, 0, A), reflectivity=bs_R, aperture=W35,
                               label="合束器（%.0f:%.0f）" % (bs_R * 100, (1 - bs_R) * 100)))
    comps.append(detector("d1", (A, 0, 2 * A), shape=circle(W35), yaw=PI, label=d1_label))
    comps.append(detector("d2", (2 * A, 0, A), shape=circle(W35), yaw=-PI / 2, label=d2_label))
    return comps


# ==========================================================================
# 1. 杨氏双缝：缝距扫掠
# ==========================================================================

def young_double_slit(reg):
    for i, d in enumerate([0.40e-3, 0.50e-3, 0.70e-3, 1.00e-3, 1.25e-3,
                           1.50e-3, 2.00e-3, 2.50e-3, 3.00e-3, 3.60e-3]):
        a, wl, z = 0.20e-3, WL["633"], 1.5
        W = W20
        N = pick_grid(a, W, minimum=4.0)
        dx = W / N
        fring = fringe_spacing(wl, z, d)
        check_sampling(dx, {"缝距 d": d, "缝宽 a": a, "条纹间距": fring}, 4.0)
        n_fr = int(W * 0.9 / fring)
        cfg = build_config(
            name="杨氏双缝干涉：缝距 d=%.2f mm（条纹间距 %.2f mm）" % (d * 1e3, fring * 1e3),
            size=N, width=W, wl=wl,
            components=[_slit(d=d, a=a),
                        det_rect("det", z, W * 0.95, W * 0.95, label="干涉屏（%.1f m 处）" % z)],
            sources=[plane_source(wl=wl, label="平面波 633 nm")],
            category="干涉", difficulty="入门",
            description="双缝（缝宽 0.20 mm、缝距 %.2f mm）被平面波照明，屏距 %.1f m；"
                        "缝宽远小于缝距，属纯双光束干涉。" % (d * 1e3, z),
            expected="等间距直条纹：Δx=λz/d=%.3f mm，窗口内约 %d 条；缝距越小条纹越疏。"
                     "远处被单缝衍射包络（第一零点 λz/a=%.2f mm）调制——这是「干涉 × 衍射」的经典组合。"
                     % (fring * 1e3, n_fr, sinc_first_zero(wl, z, a) * 1e3),
            tags=["杨氏双缝", "条纹间距", "干涉"],
            theory="Δx = λL/d", sample=(i < 2),
        )
        reg.add("interf_double_slit_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))

    # 波长扫掠
    for i, key in enumerate(["405", "450", "532", "633", "780"]):
        wl, d, a, z, W = WL[key], 0.8e-3, 0.2e-3, 1.5, W20
        N = pick_grid(a, W, minimum=4.0)
        dx = W / N
        fring = fringe_spacing(wl, z, d)
        check_sampling(dx, {"条纹间距": fring}, 4.0)
        cfg = build_config(
            name="双缝干涉：波长 %s nm（条纹间距 %.2f mm）" % (key, fring * 1e3),
            size=N, width=W, wl=wl,
            components=[_slit(d=d, a=a),
                        det_rect("det", z, W * 0.95, W * 0.95, label="干涉屏")],
            sources=[plane_source(wl=wl, label="平面波 %s nm" % key)],
            category="干涉", difficulty="入门",
            description="缝距固定 0.80 mm，只改波长（%s nm）——条纹间距与波长成正比。" % key,
            expected="Δx=λz/d=%.3f mm：红光条纹比紫光宽 %.2f 倍。用白光照明时各波长条纹错开叠加，"
                     "只有零级（光程差为 0 处）是白的，两侧出现彩色——这就是白光干涉条纹的成因。"
                     % (fring * 1e3, WL[key] / WL["405"]),
            tags=["双缝", "波长", "白光条纹"],
            theory="Δx = λz/d 与 λ 成正比", sample=(i == 0),
        )
        reg.add("interf_double_slit_wl_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))

    # 照明倾角：条纹整体平移
    for i, tilt in enumerate([1e-4, 2e-4, 5e-4, 1e-3]):
        wl, d, a, z, W = WL["633"], 1.0e-3, 0.2e-3, 1.0, W20
        N = pick_grid(a, W, minimum=4.0)
        dx = W / N
        fring = fringe_spacing(wl, z, d)
        shift = tilt * z
        check_sampling(dx, {"条纹间距": fring}, 4.0)
        cfg = build_config(
            name="双缝干涉：照明倾角 %.1f mrad → 条纹平移 %.2f mm" % (tilt * 1e3, shift * 1e3),
            size=N, width=W, wl=wl,
            components=[_slit(d=d, a=a, z=0.0),
                        det_rect("det", z, W * 0.95, W * 0.95, label="干涉屏")],
            sources=[plane_source(wl=wl, tilt_x=tilt, label="倾斜平面波 %.1f mrad" % (tilt * 1e3))],
            category="干涉", difficulty="进阶",
            description="双缝不变，把照明平面波倾斜 %.1f mrad（等价于把点光源横向移动）。" % (tilt * 1e3),
            expected="条纹**间距不变**（仍为 λz/d=%.3f mm），但整体平移 Δx=z·θ=%.3f mm——"
                     "条纹的相位由两条缝的光程差决定，倾斜照明给两条缝加了额外光程差。"
                     "这是「干涉仪对光源指向敏感」的根源。" % (fring * 1e3, shift * 1e3),
            tags=["双缝", "倾斜照明", "条纹平移"],
            theory="倾斜 θ 使条纹平移 z·θ（零级位于 sinθ = Δ/... 方向）",
        )
        reg.add("interf_double_slit_tilt_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"))


# ==========================================================================
# 2. 双孔（二维杨氏）
# ==========================================================================

def two_hole_young(reg):
    for i, (d, r) in enumerate([(1.0e-3, 0.25e-3), (1.5e-3, 0.25e-3), (2.0e-3, 0.3e-3)]):
        wl, z, W = WL["633"], 1.0, W20
        N = pick_grid(r, W, minimum=4.0)
        dx = W / N
        fring = fringe_spacing(wl, z, d)
        check_sampling(dx, {"条纹间距": fring}, 4.0)
        comps = [aperture("h1", (0, 0, 0), circle(r), parallel=True, x=-d / 2, label="孔 1"),
                 aperture("h2", (0, 0, 0), circle(r), parallel=True, x=d / 2, label="孔 2"),
                 det_rect("det", z, W * 0.95, W * 0.95, label="干涉屏")]
        cfg = build_config(
            name="双孔干涉：孔距 %.2f mm、孔径 %.2f mm" % (d * 1e3, r * 1e3),
            size=N, width=W, wl=wl, components=comps,
            sources=[plane_source(wl=wl)],
            category="干涉", difficulty="入门",
            description="两个半径 %.2f mm 的圆孔（中心距 %.2f mm）——杨氏实验的二维版本。" % (r * 1e3, d * 1e3),
            expected="平行直条纹，间距 Δx=λz/d=%.2f mm；条纹被圆孔艾里包络调制（第一暗环 %.2f mm），"
                     "包络内可见约 %.0f 条条纹。" % (fring * 1e3, airy_radius(wl, z, 2 * r) * 1e3,
                                                    airy_radius(wl, z, 2 * r) / fring),
            tags=["双孔", "二维杨氏", "艾里包络"],
            theory="Δx=λz/d，包络为孔径的艾里图样", sample=(i == 0),
        )
        reg.add("interf_two_hole_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))


# ==========================================================================
# 3. 多光束（N 缝）
# ==========================================================================

def multi_beam_interf(reg):
    for i, n in enumerate([3, 4, 5, 6, 8]):
        d, a, wl, z, W = 1.0e-3, 0.1e-3, WL["633"], 1.0, W20
        N = pick_grid(a, W, minimum=4.0)
        dx = W / N
        fring = fringe_spacing(wl, z, d)
        check_sampling(dx, {"缝距": d, "缝宽": a}, 4.0)
        comps = []
        for m in range(n):
            comps.append(aperture("s%d" % m, (0, 0, 0), rectangle(a, W * 0.95), parallel=True,
                                  x=(m - (n - 1) / 2.0) * d, label="第 %d 缝" % (m + 1)))
        comps.append(det_rect("det", z, W * 0.95, W * 0.95, label="干涉屏"))
        cfg = build_config(
            name="%d 缝干涉：主极大变锐" % n,
            size=N, width=W, wl=wl, components=comps,
            sources=[plane_source(wl=wl)],
            category="干涉", difficulty="进阶",
            description="缝宽 0.10 mm、缝距 1.00 mm 的 %d 缝阵列（同 z 上 parallel 并集）。" % n,
            expected="主极大间距仍为 λz/d=%.2f mm，但主极大随 N 增大而**变锐**（角宽度 ∝ 1/(N·d)）；"
                     "N 缝时两主极大之间有 N−2=%d 个次极大，强度上限约 1/N²=%.3f 倍主极大。"
                     % (fring * 1e3, n - 2, 1.0 / n ** 2),
            tags=["多缝干涉", "主极大", "谱线宽度"],
            theory="主极大宽度 ∝ 1/(N d)，次极大数 N−2", sample=(i < 2),
        )
        reg.add("interf_multi_beam_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))


# ==========================================================================
# 4. 牛顿环（曲率镜的等厚干涉）
# ==========================================================================

def newton_rings(reg):
    """牛顿环：两条臂的波前曲率不同 → 相位差 ∝ r²，环半径 r_m² = m λ R。

    用「一个平面镜 + 一个曲率半径 R 的球面镜」的迈克尔逊结构实现：
    球面镜相位 exp(-ik·curvature·r²)，取 curvature = -1/R 即得 +k r²/R 的额外相位，
    与空气隙厚度 t(r)=r²/(2R) 造成的 2kt=kr²/R 完全一致（环半径同公式 r_m²=mλR）。
    """
    geom = [(0.10, 0.004), (0.20, 0.006), (0.35, 0.008), (0.50, 0.010),
            (0.75, 0.012), (1.00, 0.012), (1.50, 0.016), (2.00, 0.02)]
    for i, (R, W) in enumerate(geom):
        wl = WL["633"]
        r1 = math.sqrt(wl * R)
        dr_out = wl * R / (2 * W * 0.35)
        N = pick_grid(min(r1, dr_out), W, minimum=4.0)
        dx = W / N
        check_sampling(dx, {"第一环半径": r1, "外环间距": dr_out}, 3.0)
        m_max = int((W * 0.35) ** 2 / (wl * R))
        cfg = build_config(
            name="牛顿环：曲率半径 R=%.2f m（%d 条环）" % (R, m_max),
            size=N, width=W, wl=wl,
            components=michelson_comps(
                W, L=0.06,
                m2={"curvature": -1.0 / R,
                    "label": "球面镜 R=%.2f m（等效空气隙 r²/2R）" % R},
                d1_label="反射端口：牛顿环"),
            sources=[plane_source(wl=wl, label="平面波 633 nm")],
            category="干涉", difficulty="进阶",
            description="迈克尔逊结构：参考臂是平面镜，另一臂是曲率半径 R=%.2f m 的球面镜——"
                        "两臂反射波前的曲率差等价于厚度 t(r)=r²/(2R) 的空气楔，即牛顿环的标准模型。" % R,
            expected="同心圆环：第 m 个环半径 r_m=√(mλR)，本实例 r₁=%.2f mm、r₂=%.2f mm、r₃=%.2f mm；"
                     "环越往外越密（相邻环间距 Δr=λR/(2r) 随 r 减小）。窗口内约 %d 条环。"
                     "中心是亮斑或暗斑取决于半波损失（本内核的镜面是实振幅反射，不含附加 π 相移）。"
                     % (math.sqrt(wl * R) * 1e3, math.sqrt(2 * wl * R) * 1e3,
                        math.sqrt(3 * wl * R) * 1e3, m_max),
            tags=["牛顿环", "等厚干涉", "曲率镜"],
            theory="r_m² = mλR（空气隙 t=r²/2R，2t=mλ）；相邻环间距 Δr=λR/(2r)",
            sample=(i < 2),
        )
        reg.add("interf_newton_rings_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))

    # 牛顿环的波长依赖
    for i, key in enumerate(["450", "532", "633", "780"]):
        wl, R, W = WL[key], 0.5, W10
        r1 = math.sqrt(wl * R)
        N = pick_grid(r1, W, minimum=4.0)
        cfg = build_config(
            name="牛顿环的波长依赖：λ=%s nm" % key,
            size=N, width=W, wl=wl,
            components=michelson_comps(
                W, L=0.06, m2={"curvature": -1.0 / R, "label": "球面镜 R=0.50 m"},
                d1_label="牛顿环"),
            sources=[plane_source(wl=wl, label="平面波 %s nm" % key)],
            category="干涉", difficulty="进阶",
            description="同一个球面镜（R=0.50 m）在不同波长下的牛顿环：环半径 r_m=√(mλR) 与 √λ 成正比。",
            expected="第一环半径 %.2f mm（λ=%s nm）；波长越长环越大、同样窗口内环数越少"
                     "（r₁ ∝ √λ）。白光照明时各色环错开叠加成彩色同心环，"
                     "同一级不同波长的环半径比 r(λ₁)/r(λ₂)=√(λ₁/λ₂)。"
                     % (r1 * 1e3, key),
            tags=["牛顿环", "波长", "彩色环"],
            theory="r_m ∝ √(mλR)", sample=(i == 0),
        )
        reg.add("interf_newton_rings_wl_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))


# ==========================================================================
# 5. 等厚干涉（楔形相位）
# ==========================================================================

def equal_thickness(reg):
    for i, (alpha, n_w) in enumerate([(2.0e-3, 1.5), (5.0e-3, 1.5), (1.0e-2, 1.5),
                                      (2.0e-2, 1.5), (1.0e-2, 1.7)]):
        wl, W = WL["633"], W10
        spacing = wl / ((n_w - 1.0) * alpha)
        N = pick_grid(spacing, W, minimum=4.0)
        dx = W / N
        check_sampling(dx, {"条纹间距": spacing}, 4.0)
        cfg = build_config(
            name="等厚干涉（楔形相位）：α=%.0f mrad → 条纹间距 %.2f mm"
                 % (alpha * 1e3, spacing * 1e3),
            size=N, width=W, wl=wl,
            components=mz_comps(
                W, A=0.08,
                arm_extra=[comp("wedge", "wedge", (0, 0, 0.04), shape=circle(W * 0.35),
                                params={"alpha": alpha, "index": n_w},
                                label="楔形棱镜 α=%.0f mrad（n=%.1f）" % (alpha * 1e3, n_w))],
                d1_label="输出端口"),
            sources=[plane_source(wl=wl)],
            category="干涉", difficulty="进阶",
            description="马赫-曾德尔的一臂中插入楔形棱镜（夹角 α=%.0f mrad、n=%.1f），"
                        "产生线性相位 → 等厚干涉直条纹。" % (alpha * 1e3, n_w),
            expected="等间距直条纹，间距 Δx=λ/((n−1)α)=%.3f mm（条纹平行于楔棱）；"
                     "α 越小条纹越疏，楔角增大则条纹变密。这就是「用条纹数测小角度」的原理："
                     "α=λ/((n−1)Δx)。" % (spacing * 1e3),
            tags=["等厚干涉", "楔形", "小角测量"],
            theory="Δx = λ/((n−1)α)", sample=(i < 2),
        )
        reg.add("interf_equal_thickness_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))


# ==========================================================================
# 6. 迈克尔逊：臂长差与微倾条纹
# ==========================================================================

def michelson_fringes(reg):
    # 微倾 → 直条纹，间距 λ/(2Δθ)
    for i, tilt in enumerate([1.0e-4, 2.0e-4, 5.0e-4, 1.0e-3]):
        wl, W = WL["633"], W10
        spacing = wl / (2.0 * tilt)
        N = pick_grid(spacing, W, minimum=4.0)
        dx = W / N
        check_sampling(dx, {"条纹间距": spacing}, 4.0)
        cfg = build_config(
            name="迈克尔逊等厚条纹：端镜微倾 %.1f mrad → 间距 %.2f mm" % (tilt * 1e3, spacing * 1e3),
            size=N, width=W, wl=wl,
            components=michelson_comps(
                W, L=0.05,
                m2={"tilt_y": tilt, "label": "端镜 B（微倾 %.1f mrad）" % (tilt * 1e3)},
                d1_label="回光源端口", d2_label="侧端口"),
            sources=[plane_source(wl=wl)],
            category="干涉", difficulty="进阶",
            description="迈克尔逊干涉仪，两臂等长但端镜 B 倾斜 %.1f mrad——两束回射光之间有 2θ 的夹角，"
                        "形成等厚（等倾）直条纹。" % (tilt * 1e3),
            expected="直条纹间距 Δx=λ/(2θ)=%.3f mm；倾角越大条纹越密。两个输出端口互补："
                     "一个端口亮条纹的位置在另一端口是暗条纹（能量守恒）。"
                     "折返光路的相位按 2kL 累加（往返 = 两倍臂长）。" % (spacing * 1e3),
            tags=["迈克尔逊", "微倾", "等厚条纹", "互补端口"],
            theory="Δx = λ/(2θ)；折返相位 2kL", sample=(i < 2),
        )
        reg.add("interf_michelson_tilt_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))

    # 臂长差扫描：λ/4 步进 → 端口在亮/暗之间循环
    for i, frac in enumerate([0.0, 0.25, 0.5, 0.75]):
        wl, W = WL["633"], W10
        d = frac * wl / 2.0
        N = 512
        cfg = build_config(
            name="迈克尔逊臂长差扫描：ΔL=%.2f λ → 相位差 %.0f°" % (frac, frac * 360),
            size=N, width=W, wl=wl,
            components=[
                beam_splitter("bs", (0, 0, 0), reflectivity=0.5, aperture=W * 0.35,
                              label="分束器 50:50"),
                mirror("m1", (0, 0, 0.05), reflectivity=1.0, yaw=0.0,
                       shape=circle(W * 0.35), label="端镜 A（固定）"),
                mirror("m2", (0.05 + d, 0, 0), reflectivity=1.0, yaw=PI / 2,
                       shape=circle(W * 0.35), label="端镜 B（前移 %.0f nm）" % (d * 1e9)),
                d_port_z("d1", -0.05, "回光源端口", W * 0.35),
                d_port_nx("d2", -0.05, "侧端口", W * 0.35),
            ],
            sources=[plane_source(wl=wl)],
            category="干涉", difficulty="进阶",
            description="端镜 B 沿光轴前移 ΔL=%.0f nm，两臂光程差 2ΔL=%.2f λ（相位差 2kΔL）。"
                        % (d * 1e9, 2 * d / wl),
            expected="折返使相位差 = 2kΔL：光程差每变化 λ/2 端口就在亮—暗之间完整循环一次。"
                     "本实例相位差 %.0f°：等臂（ΔL=0）时两臂在回光源端口**反相相消**（该端口全暗、"
                     "侧端口全亮），端镜前移 λ/4 后两端口互换。实测：回光源端口 %s。"
                     % (frac * 360, "全暗" if frac < 0.1 else ("全亮" if abs(frac - 0.5) < 0.1 else "半亮")),
            tags=["迈克尔逊", "臂长差", "相位扫描"],
            theory="相位差 = 2kΔL（折返两倍），周期为 λ/2", sample=(i < 2),
        )
        reg.add("interf_michelson_scan_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))


# ==========================================================================
# 7. 马赫-曾德尔：臂内相位板扫描
# ==========================================================================

def mach_zehnder_phase(reg):
    for i, frac in enumerate([0.0, 0.25, 0.5, 0.75, 1.0]):
        wl, W = WL["633"], W10
        d_plate = 0.01
        n_plate = 1.0 + frac * wl / d_plate          # 额外相位 = k(n-1)d = frac·2π
        N = 512
        cfg = build_config(
            name="马赫-曾德尔相位扫描：相位板 %.2f×2π" % frac,
            size=N, width=W, wl=wl,
            components=mz_comps(
                W, A=0.08,
                arm_extra=[comp("plate", "medium", (0, 0, 0.04), shape=circle(W * 0.4),
                                params={"distance": d_plate, "index": n_plate,
                                        "absorption": 0.0, "steps": 4},
                                label="相位板（额外相位 %.2f×2π）" % frac)],
                d1_label="端口 1", d2_label="端口 2"),
            sources=[plane_source(wl=wl)],
            category="干涉", difficulty="进阶",
            description="平衡马赫-曾德尔干涉仪，在一条臂中插入厚度 10 mm 的相位板（折射率 n=%.7f），"
                        "使两臂相差 %.2f×2π。" % (n_plate, frac),
            expected="端口功率按 cos²(φ/2) 与 sin²(φ/2) 分配：相差 0 时端口 1 全亮、端口 2 全暗；"
                     "相差 π（0.5×2π）时两端口反转；两个端口功率之和恒等于输入功率（能量守恒）。"
                     "本实例相差 %.2f×2π，端口 1 应为%s。"
                     % (frac, "全亮" if frac < 0.1 else ("全暗" if abs(frac - 0.5) < 0.1 else "部分亮")),
            tags=["马赫-曾德尔", "相位扫描", "能量守恒", "互补端口"],
            theory="相差 φ 时端口功率 = (1±cosφ)/2；两端口互补", sample=(i < 3),
        )
        reg.add("interf_mach_zehnder_phase_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))


# ==========================================================================
# 8. 双（相干）光源条纹：同组内两个倾斜平面波
# ==========================================================================

def two_source_fringes(reg):
    for i, (half_ang, wl) in enumerate([(2.5e-5, WL["633"]), (5e-5, WL["633"]),
                                       (1e-4, WL["633"]), (2e-4, WL["633"]),
                                       (1e-4, WL["532"]), (1e-4, WL["780"])]):
        W = W20
        spacing = wl / (2.0 * half_ang)
        N = pick_grid(spacing, W, minimum=4.0)
        dx = W / N
        check_sampling(dx, {"条纹间距": spacing}, 4.0)
        cfg = build_config(
            name="双光束干涉：两束夹角 %.1f mrad → 条纹间距 %.2f mm"
                 % (2 * half_ang * 1e3, spacing * 1e3),
            size=N, width=W, wl=wl,
            components=[det_rect("det", 0.5, W * 0.95, W * 0.95, label="干涉平面")],
            sources=[
                plane_source(wl=wl, power=0.5e-3, tilt_x=+half_ang, group="A",
                             sid="s1", label="光束 1（+%.2f mrad）" % (half_ang * 1e3)),
                plane_source(wl=wl, power=0.5e-3, tilt_x=-half_ang, group="A",
                             sid="s2", label="光束 2（−%.2f mrad）" % (half_ang * 1e3)),
            ],
            category="干涉", difficulty="进阶",
            description="两个**同组同波长**的倾斜平面波（夹角 %.1f mrad）在 0.5 m 处叠加——"
                        "最纯粹的双光束干涉（不需要任何分束器）。" % (2 * half_ang * 1e3),
            expected="等间距直条纹，间距 Δx=λ/(2sinθ)≈λ/(2θ)=%.3f mm，对比度接近 1（两束等功率）；"
                     "条纹方向垂直于两波矢之差。把其中之一换成不同波长或不同 group，条纹立刻消失"
                     "（只看强度相加）——相干性的定义。" % (spacing * 1e3),
            tags=["双光束", "相干叠加", "相干组"],
            theory="两平面波夹角 θ：Δx = λ/(2 sin(θ/2)) ≈ λ/θ；相干组决定是否叠加振幅",
            sample=(i < 2),
        )
        reg.add("interf_two_beam_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))

    # 相干 vs 非相干：同一几何，分组不同
    W, wl = W20, WL["633"]
    for i, (grp, lbl) in enumerate([("A", "同相干组（相干叠加）"), ("B", "分属两组（强度相加）")]):
        spacing = wl / (2 * 1e-4)
        N = pick_grid(spacing, W, minimum=4.0)
        cfg = build_config(
            name="相干性对比：两束光的 %s" % lbl,
            size=N, width=W, wl=wl,
            components=[det_rect("det", 0.5, W * 0.95, W * 0.95, label="观察面")],
            sources=[
                plane_source(wl=wl, power=0.5e-3, tilt_x=+1e-4, group="A", sid="s1",
                             label="光束 1（组 A）"),
                plane_source(wl=wl, power=0.5e-3, tilt_x=-1e-4, group=grp, sid="s2",
                             label="光束 2（组 %s）" % grp),
            ],
            category="干涉", difficulty="进阶",
            description="几何完全相同（两束夹角 0.2 mrad 的平面波），只有**相干组**不同：%s。" % lbl,
            expected=("看到间距 λ/(2θ)=3.16 mm 的高对比度条纹（振幅相干叠加，出现相消与相长）。"
                      if grp == "A" else
                      "条纹完全消失，只剩均匀的强度相加（互不相干的光只叠加强度，不产生干涉项）。"),
            tags=["相干性", "相干组", "强度叠加"],
            theory="同组同波长 → 相干叠加；否则只叠加强度", sample=True,
        )
        reg.add("interf_coherence_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=True)


# ==========================================================================
# 9. 平行平板 / 菲索条纹（部分反射镜 + 后置全反镜）
# ==========================================================================

def parallel_plate(reg):
    for i, (gap, R) in enumerate([(10e-6, 0.5), (20e-6, 0.5), (40e-6, 0.5), (20e-6, 0.8)]):
        wl, W = WL["633"], W10
        tilt = 1.5e-4
        spacing = wl / (2 * tilt)
        N = pick_grid(spacing, W, minimum=4.0)
        dx = W / N
        check_sampling(dx, {"条纹间距": spacing}, 4.0)
        # 板间距是**纵向**尺寸（只贡献 2k·gap 的相位与很小的横向位移），不需要横向采样
        n_wave = 2.0 * gap / wl
        cfg = build_config(
            name="平行平板（菲索）条纹：间距 %.0f µm、部分反射镜 R=%.1f" % (gap * 1e6, R),
            size=N, width=W, wl=wl,
            components=[
                beam_splitter("bs", (0, 0, 0), reflectivity=0.5, aperture=W * 0.35,
                              label="分束器 50:50"),
                mirror("mflat", (0, 0, 0.05), reflectivity=1.0, yaw=0.0,
                       shape=circle(W * 0.35), label="参考臂端镜（平面）"),
                mirror("mpart", (0.05, 0, 0), reflectivity=math.sqrt(R), yaw=PI / 2,
                       shape=circle(W * 0.35),
                       label="部分反射镜 ρ=%.3f（功率 R=%.1f）" % (math.sqrt(R), R)),
                mirror("mback", (0.05 + gap, 0, 0), reflectivity=1.0, yaw=PI / 2,
                       tilt_y=tilt, shape=circle(W * 0.35),
                       label="后置全反镜（微倾 %.1f mrad、间距 %.0f µm）" % (tilt * 1e3, gap * 1e6)),
                d_port_z("d1", -0.05, "反射端口", W * 0.35),
            ],
            sources=[plane_source(wl=wl)],
            category="干涉", difficulty="高级",
            description="一臂中放两块平行镜（前镜部分反射 ρ=√R、后镜全反），间距 %.0f µm —— "
                        "这一臂返回两束光（前镜直接反射 + 后镜折返后透射），等价于**平行平板/菲索干涉仪**。"
                        % (gap * 1e6),
            expected="两束光的光程差 2ngap=%.1fλ 与后镜微倾共同决定条纹：直条纹间距 λ/(2θ)=%.2f mm，"
                     "对比度由部分反射镜的反射率决定（R=%.1f 时两束振幅比大致 %.2f）；"
                     "板间距每变化 λ/2 相位差变化 2π（可用来做波长计/厚度计）。"
                     "本内核只计首次往返（两光束近似），高反射率的 Airy 多光束条纹不在结果中。"
                     % (2 * gap / wl, spacing * 1e3, R, math.sqrt(R)),
            tags=["平行平板", "菲索", "部分反射镜", "两光束近似"],
            theory="平行平板两光束近似：光程差 2ngap；多光束 Airy 条纹需多次往返（本内核只算首次往返）",
        )
        reg.add("interf_parallel_plate_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"))


# ==========================================================================
# 10. 条纹可见度研究（不等功率比）
# ==========================================================================

def visibility_study(reg):
    for i, R in enumerate([0.1, 0.3, 0.5, 0.7, 0.9]):
        wl, W = WL["633"], W10
        N = 512
        i1, i2 = R, 1.0 - R
        vis = 2.0 * math.sqrt(i1 * i2) / (i1 + i2)
        cfg = build_config(
            name="条纹可见度 vs 分束比：R=%.1f → V=%.2f" % (R, vis),
            size=N, width=W, wl=wl,
            components=michelson_comps(
                W, L=0.05, bs_R=R,
                m2={"tilt_y": 1.5e-4, "label": "端镜 B（微倾 0.15 mrad）"},
                d1_label="条纹（回光源端口）"),
            sources=[plane_source(wl=wl)],
            category="干涉", difficulty="高级",
            description="迈克尔逊干涉仪，分束器反射率 R=%.1f（两臂功率比 %.1f:%.1f）——"
                        "研究**条纹可见度** V=(Imax−Imin)/(Imax+Imin) 与臂功率比的关系。" % (R, i1, i2),
            expected="理论可见度 V=2√(I₁I₂)/(I₁+I₂)=%.3f；R=0.5（等功率）时 V=1（条纹最清晰），"
                     "臂功率比越悬殊可见度越低，但当 I₂→0 时条纹消失而背景仍在。"
                     "注意分束器功率互补（R+(1−R)=1），能量守恒。" % vis,
            tags=["可见度", "分束比", "能量守恒"],
            theory="V = 2√(I₁I₂)/(I₁+I₂)；R=0.5 时 V=1", sample=(i == 2),
        )
        reg.add("interf_visibility_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))
