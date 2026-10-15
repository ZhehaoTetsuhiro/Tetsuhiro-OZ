"""衍射类示例（第二批）：光栅族、波带片、泰伯自成像、直边衍射、涡旋远场、像差 PSF、4f 滤波。"""

import math

from wavegen.common import *   # noqa: F401,F403

W20 = 0.02
W8 = 0.008
W4 = 0.004
W10 = 0.01

# 傅里叶透镜 + 焦面探测器：光栅级次位置 x_q = q f λ / Λ
FF_F = 0.3
FF_Z = 0.32


def _ff_comps(comps, W=W8):
    """在元件列表后追加傅里叶透镜 + 焦面探测器。"""
    out = list(comps)
    out.append(lens("fl", (0, 0, 0.02), f=FF_F, aperture=W * 0.49,
                    label="傅里叶透镜 f=%.2f m" % FF_F))
    out.append(det_rect("det", FF_Z, W * 0.95, W * 0.95, label="焦面（衍射级）"))
    return out


def _mk_ff(name, comps, desc, expect, tags, theory, difficulty="进阶", N=1024, W=W8,
           wl=None, sample=False, label="焦面（衍射级）"):
    """傅里叶透镜焦面几何的 Config（光栅族统一入口）。"""
    wl = wl or WL["633"]
    out = list(comps)
    out.append(lens("fl", (0, 0, 0.02), f=FF_F, aperture=W * 0.49,
                    label="傅里叶透镜 f=%.2f m" % FF_F))
    out.append(det_rect("det", FF_Z, W * 0.95, W * 0.95, label=label))
    return build_config(name=name, size=N, width=W, wl=wl, components=out,
                        sources=[plane_source(wl=wl)],
                        category="衍射", difficulty=difficulty, description=desc,
                        expected=expect, tags=tags, theory=theory, sample=sample)


def _reg(reg, filename, cfg, sample=False):
    """把 _mk_ff 产出的 Config 登记进 Registry（元数据已内嵌）。"""
    z = cfg["meta"]["zod"]
    reg.add(filename, z["category"], z["difficulty"], z["description"], z["expected"],
            cfg, tags=z["tags"], theory=z.get("theory"), sample=sample)


# ==========================================================================
# 光栅 1：正弦振幅光栅 —— 级次权重 1 : (m/2)² : (m/2)²
# ==========================================================================

def gratings_sin_amp(reg):
    for i, (Lam, m) in enumerate([(60e-6, 0.8), (80e-6, 0.8), (100e-6, 0.8),
                                  (150e-6, 0.8), (200e-6, 0.8)]):
        N, W = 1024, W8
        dx = W / N
        check_sampling(dx, {"光栅周期 Λ": Lam}, 5.0)
        x1 = grating_order_pos(1, FF_F, WL["633"], Lam)
        window_ok(W, {"±1 级": x1})
        check_sampling(dx, {"±1 级位置": x1}, 20.0)
        eff1 = (m / 2.0) ** 2 / (1 + 2 * (m / 2.0) ** 2)
        nm = "正弦振幅光栅：周期 Λ=%.0f µm" % (Lam * 1e6)
        cfg = _mk_ff(
            name=nm, difficulty="入门",
            comps=[comp("g", "grating", (0, 0, 0),
                        params={"kind": "amplitude_sin", "period": Lam, "modulation": m},
                        label="正弦振幅光栅 Λ=%.0f µm" % (Lam * 1e6))],
            desc="透过率 t(x)=[1+m·cos(2πx/Λ)]/(1+m) 的连续振幅光栅（Λ=%.0f µm，对比度 m=%.2f），"
                 "经焦距 %.2f m 的傅里叶透镜到焦面读数。" % (Lam * 1e6, m, FF_F),
            expect="只有 0 级与 ±1 级：±1 级位于 x=±fλ/Λ=±%.2f mm，级次间强度比 I₁:I₀ = (m/2)²:1，"
                   "归一化后各级份额 0 级 ≈%.1f%%、±1 级各 ≈%.1f%%；±2 级及以上严格为零（正弦光栅只有三束衍射光）。"
                   % (x1 * 1e3, 100 * (1 - 2 * eff1), 100 * eff1),
            tags=["正弦振幅光栅", "三束衍射", "傅里叶光学"],
            theory="正弦振幅光栅 → 仅 0、±1 三束；I±1/I0 = (m/2)²", sample=(i == 0))
        _reg(reg, "diff_grating_sin_amp_%02d.json" % (i + 1), cfg, sample=(i == 0))


# ==========================================================================
# 光栅 2：二元振幅光栅 + 占空比（缺级）
# ==========================================================================

def gratings_binary_amp(reg):
    for i, duty in enumerate([0.25, 0.4, 0.5, 0.6, 0.75]):
        Lam, N, W = 100e-6, 1024, W8
        dx = W / N
        check_sampling(dx, {"周期": Lam}, 5.0)
        x1 = grating_order_pos(1, FF_F, WL["633"], Lam)
        missing = "偶数级（±2、±4…）全部缺失" if abs(duty - 0.5) < 1e-9 else \
                  "缺口条件不满足，各级次都存在"
        cfg = _mk_ff(
            name="二元振幅光栅：占空比 %.2f" % duty,
            comps=[comp("g", "grating", (0, 0, 0),
                        params={"kind": "amplitude_binary", "period": Lam, "duty": duty},
                        label="二元振幅光栅 占空比 %.2f" % duty)],
            desc="占空比 %.2f 的方波振幅光栅（周期 100 µm、缝宽 %0.1f µm），傅里叶透镜焦面观察级次。"
                 % (duty, duty * 100),
            expect="第 q 级强度 ∝ sinc²(q·duty)：%s。±1 级在 x=±%.2f mm，零级最亮（占 %.0f%%）。"
                   "占空比即「缝宽/周期」，因此光栅的缺级规律可以直接由它控制。"
                   % (missing, x1 * 1e3, 100 * duty ** 2 if duty != 1 else 100),
            tags=["二元光栅", "占空比", "缺级"],
            theory="I_q ∝ sinc²(q·duty)，duty=1/2 时偶数级缺失")
        _reg(reg, "diff_grating_bin_amp_%02d.json" % (i + 1), cfg)


# ==========================================================================
# 光栅 3：正弦相位光栅（Raman-Nath）—— 贝塞尔级次
# ==========================================================================

def gratings_phase_sin(reg):
    for i, m in enumerate([0.4, 0.8, 1.2, 1.6, 2.0]):
        Lam, N, W = 100e-6, 1024, W8
        dx = W / N
        check_sampling(dx, {"周期": Lam}, 5.0)
        x1 = grating_order_pos(1, FF_F, WL["633"], Lam)
        # Raman-Nath: 第 q 级振幅 = J_q(m/2)
        try:
            from math import factorial
            def J(n, x):
                s = 0.0
                for kk in range(0, 14):
                    s += (-1) ** kk / (factorial(kk) * factorial(kk + n)) * (x / 2.0) ** (2 * kk + n)
                return s
            j0, j1, j2 = J(0, m / 2.0), J(1, m / 2.0), J(2, m / 2.0)
            r10 = (j1 / j0) ** 2
            r21 = (j2 / j1) ** 2
        except Exception:
            r10 = r21 = float("nan")
        cfg = _mk_ff(
            name="正弦相位光栅（Raman-Nath）：峰峰相位 %.1f rad" % m,
            comps=[comp("g", "grating", (0, 0, 0),
                        params={"kind": "phase_sin", "period": Lam, "modulation": m},
                        label="正弦相位光栅 m=%.1f rad" % m)],
            desc="薄相位光栅 t(x)=exp[i(m/2)sin(2πx/Λ)]，峰峰相位 m=%.1f rad、周期 100 µm。" % m,
            expect="与振幅光栅不同，相位光栅有**多级**衍射：第 q 级振幅为贝塞尔函数 J_q(m/2)。"
                   "理论比 I₁/I₀=(J₁/J₀)²=%.3f、I₂/I₁=(J₂/J₁)²=%.3f。m=%.1f rad 时 0 级仍最强；"
                   "m→2.405 rad 时 J₀=0，零级完全消失（光栅「点亮」）。"
                   % (r10, r21, m),
            tags=["相位光栅", "Raman-Nath", "贝塞尔级次"],
            theory="第 q 级振幅 = J_q(m/2)，能量守恒 Σ J_q² = 1", sample=(i == 0))
        _reg(reg, "diff_grating_phase_sin_%02d.json" % (i + 1), cfg, sample=(i == 0))


# ==========================================================================
# 光栅 4：二元相位光栅（π 相位）
# ==========================================================================

def gratings_phase_binary(reg):
    for i, duty in enumerate([0.3, 0.5, 0.7, 0.25]):
        Lam, N, W = 100e-6, 1024, W8
        dx = W / N
        check_sampling(dx, {"周期": Lam}, 5.0)
        x1 = grating_order_pos(1, FF_F, WL["633"], Lam)
        cfg = _mk_ff(
            name="二元相位光栅：占空比 %.2f（π 相位台阶）" % duty,
            comps=[comp("g", "grating", (0, 0, 0),
                        params={"kind": "phase_binary", "period": Lam, "duty": duty,
                                "modulation": 2.0},
                        label="二元相位光栅 占空比 %.2f" % duty)],
            desc="占空比 %.2f 的 0/π 二元相位光栅（周期 100 µm）——表面浮雕像的经典模型，理论衍射效率"
                 "在 ±1 级各约 40.5%%。" % duty,
            expect="等占空比（0.5）时 ±1 级最强、零级被抵消（0/π 台阶的傅里叶系数在偶数级为零）；"
                   "±1 级位置 x=±%.2f mm。偏离 0.5 后零级重新变亮。这就是「相位光栅比振幅光栅亮」"
                   "（无吸收损失）的原因。" % (x1 * 1e3),
            tags=["二元相位光栅", "浮雕", "衍射效率"],
            theory="0/π 方波相位：奇数级 ∝ 1/q²，偶数级为零")
        _reg(reg, "diff_grating_bin_phase_%02d.json" % (i + 1), cfg)


# ==========================================================================
# 光栅 5：闪耀光栅（把能量搬进 +1 级）
# ==========================================================================

def gratings_blazed(reg):
    for i, dep in enumerate([math.pi / 4, math.pi / 2, math.pi, 3 * math.pi / 2, 2 * math.pi]):
        Lam, N, W = 100e-6, 1024, W8
        dx = W / N
        check_sampling(dx, {"周期": Lam}, 5.0)
        x1 = grating_order_pos(1, FF_F, WL["633"], Lam)
        frac = dep / (2 * math.pi)
        cfg = _mk_ff(
            name="闪耀光栅：锯齿相位深度 %.2f rad" % dep,
            comps=[comp("g", "grating", (0, 0, 0),
                        params={"kind": "blazed", "period": Lam, "blaze_depth": dep},
                        label="闪耀光栅 深度 %.2f rad" % dep)],
            desc="锯齿（blazed）相位光栅，周期 100 µm、峰峰相位 %.2f rad = %.2f×2π；闪耀角把 0 级能量整体"
                 "搬向 +1 级。" % (dep, frac),
            expect="峰峰相位 2π 时 +1 级拿到 ~100%% 的功率（其余级次为零）——「闪耀」的本质是让每个周期的"
                   "相位斜坡恰好补偿 +1 级在该方向的相位差；深度偏离 2π 时能量按 sinc² 泄漏到零级与其它级。"
                   "本实例深度 %.2f rad（%.0f%% of 2π），+1 级位于 x=+%.2f mm。"
                   % (dep, 100 * frac, x1 * 1e3),
            tags=["闪耀光栅", "blaze", "光栅光谱仪核心"],
            theory="闪耀条件：相位斜坡 2πp 时全部能量进入 +1 级")
        _reg(reg, "diff_grating_blazed_%02d.json" % (i + 1), cfg)


# ==========================================================================
# 光栅 6：方位旋转 + 相位偏移
# ==========================================================================

def gratings_rotation(reg):
    for i, rot in enumerate([0.0, math.pi / 6, math.pi / 4, math.pi / 2]):
        Lam, N, W = 100e-6, 1024, W8
        dx = W / N
        check_sampling(dx, {"周期": Lam}, 5.0)
        x1 = grating_order_pos(1, FF_F, WL["633"], Lam)
        cfg = _mk_ff(
            name="光栅方位旋转 %.0f°：级次方向随光栅转动" % math.degrees(rot),
            difficulty="入门",
            comps=[comp("g", "grating", (0, 0, 0),
                        params={"kind": "phase_sin", "period": Lam, "modulation": 1.6,
                                "rotation": rot},
                        label="光栅旋转 %.0f°" % math.degrees(rot))],
            desc="相位光栅绕光轴旋转 %.0f°（沟槽方向随之转动），傅里叶透镜焦面观察。" % math.degrees(rot),
            expect="衍射级次始终沿**垂直于沟槽**的方向展开：光栅转 %.0f°，级次连线也转 %.0f°；"
                   "±1 级到中心的距离不变（%.2f mm）——这决定了光谱仪中「色散方向」与「沟槽方向」正交的几何。"
                   % (math.degrees(rot), math.degrees(rot), x1 * 1e3),
            tags=["光栅旋转", "色散方向", "沟槽方向"],
            theory="色散方向 ⟂ 沟槽方向；级次位置只取决于周期与方向", sample=(i == 0))
        _reg(reg, "diff_grating_rot_%02d.json" % (i + 1), cfg, sample=(i == 0))


# ==========================================================================
# 光栅 7：缺级（双缝 = 光栅的 a/d 比值）
# ==========================================================================

def gratings_missing(reg):
    for i, (a, d) in enumerate([(0.3e-3, 0.6e-3), (0.3e-3, 0.9e-3), (0.2e-3, 0.8e-3)]):
        wl, z, N, W = WL["633"], 1.5, 1024, W20
        dx = W / N
        check_sampling(dx, {"缝宽 a": a, "缝距 d": d}, 2.0)
        fring = fringe_spacing(wl, z, d)
        order_missing = 1 if abs(d / a - 2) < 1e-6 else (3 if abs(d / a - 3) < 1e-6 else None)
        cfg = build_config(
            name="光栅缺级：a/d=1/%.0f 时第 %.0f 级缺失" % (round(d / a), round(d / a)),
            size=N, width=W, wl=wl,
            components=[aperture("ds", (0, 0, 0), double_slit(d, a, W * 0.95),
                                 label="双缝 a/d=1/%.0f" % round(d / a)),
                        det_rect("det", z, W * 0.95, W * 0.95, label="远场")],
            sources=[plane_source(wl=wl)],
            category="衍射", difficulty="高级",
            description="双缝的缝距 d、缝宽 a 满足 a/d=1/%.0f，这是最简的「缺级」演示（N 缝光栅缺级的最简情形）。"
                        % round(d / a),
            expected="干涉级次 x_m=m·λz/d（间距 %.3f mm）中，**第 %.0f 级恰好落在衍射包络的零点**"
                     "（包络零点 λz/a=%.2f mm）而被抹掉 → 该级「缺失」；其余级次正常。"
                     "这是光栅设计里用占空比控制缺级的原理。"
                     % (fring * 1e3, d / a, sinc_first_zero(wl, z, a) * 1e3),
            tags=["缺级", "包络零点", "占空比设计"],
            theory="缺级条件：d/a = 整数 m 时第 m 级消失",
            sample=(i == 0),
        )
        reg.add("diff_grating_missing_%02d.json" % (i + 1), "衍射", "高级",
                cfg["meta"]["zod"]["description"], cfg["meta"]["zod"]["expected"], cfg,
                tags=cfg["meta"]["zod"]["tags"], theory=cfg["meta"]["zod"].get("theory"),
                sample=cfg["meta"].get("sample", False))


# ==========================================================================
# 光栅 8：柱面傅里叶透镜做成的光栅光谱仪（色散方向聚焦、沟槽方向不变）
# ==========================================================================

def grating_spectrum(reg):
    combos = [(1.0e-4, "633"), (1.0e-4, "532"), (1.5e-4, "633"), (0.6e-4, "633"),
              (2.0e-4, "633"), (1.0e-4, "405")]
    for i, (Lam, key) in enumerate(combos):
        wl, N, W = WL[key], 1024, W8
        dx = W / N
        check_sampling(dx, {"周期 Λ": Lam}, 5.0)
        x1 = grating_order_pos(1, FF_F, wl, Lam)
        check_sampling(dx, {"±1 级": x1}, 20.0)
        window_ok(W, {"±1 级": x1})
        cfg = build_config(
            name="光栅光谱仪（柱面傅里叶透镜）：Λ=%.0f µm、λ=%s nm" % (Lam * 1e6, key),
            size=N, width=W, wl=wl,
            components=[
                comp("g", "grating", (0, 0, 0),
                     params={"kind": "phase_sin", "period": Lam, "modulation": 2.0},
                     label="光栅 Λ=%.0f µm" % (Lam * 1e6)),
                cyl_lens("cl", (0, 0, 0.02), f=FF_F, rotation=0.0,
                         label="柱面傅里叶透镜 f=%.2f m（只沿 x 聚焦）" % FF_F),
                det_rect("det", FF_Z, W * 0.95, W * 0.95, label="焦面：衍射级成谱线"),
            ],
            sources=[plane_source(wl=wl, label="平面波 %s nm" % key)],
            category="衍射", difficulty="进阶",
            description="光栅后接柱面透镜：只沿色散方向（x）聚焦、沟槽方向（y）完全不变，"
                        "因此每个级次是一条垂直于色散轴的谱线，而不是一个点。",
            expected="级次是**谱线**：I(y)/I(0)=1 在整行上成立（y 方向严格不变），x 方向位置 x_q=q·fλ/Λ；"
                     "±1 级位于 ±%.3f mm（λ=%s nm）。用球面透镜会把 y 也聚焦，级次退化成 λf/W_y≈%.0f µm 的点。"
                     % (x1 * 1e3, key, wl * FF_F / (W * 0.95) * 1e6),
            tags=["柱面透镜", "光谱仪", "谱线", "一维几何"],
            theory="x_q = q f λ / Λ；柱面透镜对 y 不做任何事", sample=(i == 0),
        )
        reg.add("diff_grating_spectrum_%02d.json" % (i + 1), "衍射", "进阶",
                cfg["meta"]["zod"]["description"], cfg["meta"]["zod"]["expected"], cfg,
                tags=cfg["meta"]["zod"]["tags"], theory=cfg["meta"]["zod"].get("theory"),
                sample=cfg["meta"].get("sample", False))


# ==========================================================================
# 波带片
# ==========================================================================

def zone_plates(reg):
    specs = [(0.02, "phase"), (0.03, "phase"), (0.05, "phase"), (0.08, "phase"),
             (0.10, "phase"), (0.05, "amplitude"), (0.03, "amplitude"), (0.15, "phase")]
    for i, (f, kind) in enumerate(specs):
        N, W, wl = 1024, W4, WL["633"]
        dx = W / N
        m_target = 8
        r_out = zone_radius(m_target, wl, f)
        if r_out > W * 0.45:
            r_out = W * 0.45
        m_act = zone_count(r_out, wl, f)
        dr_out = wl * f / (2 * r_out)
        check_sampling(dx, {"外环宽度": dr_out}, 4.0)
        spot = wl * f / (2 * r_out)
        cfg = build_config(
            name="菲涅耳波带片（%s）：f=%.0f mm，%d 个波带" % ("相位型 0/π" if kind == "phase" else "振幅型", f * 1e3, m_act),
            size=N, width=W, wl=wl,
            components=[
                comp("zp", "zone_plate", (0, 0, 0),
                     shape=circle(r_out),
                     params={"f": f, "radius": r_out, "kind": kind},
                     label="波带片 f=%.0f mm" % (f * 1e3)),
                det_circle("det", f, max(2 * spot, 20 * dx), label="一级焦点"),
            ],
            sources=[plane_source(wl=wl)],
            category="衍射", difficulty="进阶",
            description="环半径 r_m²=mλf+(mλ/2)² 的%s波带片，共 %d 个波带（外半径 %.2f mm，最外环宽 %.0f µm）。"
                        % ("相位" if kind == "phase" else "振幅", m_act, r_out * 1e3, dr_out * 1e6),
            expected="光严格汇聚到 **f=%.0f mm** 处（一级焦点，不是 2f）：焦斑尺度 λf/(2r_out)≈%.1f µm。"
                     "相位型（0/π）效率理论值 ≈40.5%%，振幅型（明暗交替）≈10%%——两者焦点位置相同但亮度差 4 倍。"
                     % (f * 1e3, spot * 1e6),
            tags=["波带片", "菲涅耳环", "X射线光学"],
            theory="r_m² = mλf + (mλ/2)²；相位型一级效率 4/π²≈40.5%",
            sample=(i < 2),
        )
        reg.add("diff_zone_plate_%02d.json" % (i + 1), "衍射", "进阶",
                cfg["meta"]["zod"]["description"], cfg["meta"]["zod"]["expected"], cfg,
                tags=cfg["meta"]["zod"]["tags"], theory=cfg["meta"]["zod"].get("theory"),
                sample=cfg["meta"].get("sample", False))

    # 波长色散：波带片是强色散元件
    for i, key in enumerate(["450", "532", "633", "780"]):
        N, W, f = 1024, W4, 0.05
        wl = WL[key]
        dx = W / N
        r_out = min(zone_radius(8, wl, f), W * 0.45)
        cfg = build_config(
            name="波带片色散：λ=%s nm 时一级焦点移到 f·λ₀/λ" % key,
            size=N, width=W, wl=wl,
            components=[
                comp("zp", "zone_plate", (0, 0, 0), shape=circle(r_out),
                     params={"f": 0.05, "radius": r_out, "kind": "phase"},
                     label="波带片（设计 f=50 mm @633 nm）"),
                det_circle("det", 0.05, 0.2e-3, label="设计焦点位置"),
            ],
            sources=[plane_source(wl=wl)],
            category="衍射", difficulty="高级",
            description="波带片的环半径按 λ₀=633 nm、f=50 mm 写死，此处用 λ=%s nm 照明——衍射元件的"
                        "强色散：焦点位置随波长移动。" % key,
            expected="焦点 z ≈ f·λ₀/λ = %.1f mm（λ=%s nm）；波长偏离设计值时环半径不再满足 r²=mλf，"
                     "焦斑弥散、轴上强度下降。这是衍射光学元件（DOE）与折射透镜（无焦移）最本质的区别。"
                     % (0.05 * 633e-9 / wl * 1e3, key),
            tags=["波带片", "色散", "DOE"],
            theory="衍射元件焦点随波长移动 f' = f λ₀/λ",
        )
        reg.add("diff_zone_plate_disp_%02d.json" % (i + 1), "衍射", "高级",
                cfg["meta"]["zod"]["description"], cfg["meta"]["zod"]["expected"], cfg,
                tags=cfg["meta"]["zod"]["tags"], theory=cfg["meta"]["zod"].get("theory"))


# ==========================================================================
# 泰伯自成像（Talbot）
# ==========================================================================

def talbot(reg):
    Lam, wl = 100e-6, WL["633"]
    zT = talbot_distance(wl, Lam)
    for i, frac in enumerate([0.0, 0.25, 0.5, 0.75, 1.0]):
        N, W = 1024, 0.002
        dx = W / N
        check_sampling(dx, {"周期 Λ": Lam}, 10.0)
        z = frac * zT
        zname = {0.0: "光栅面", 0.25: "z_T/4（分数泰伯）", 0.5: "z_T/2（半周期平移像）",
                 0.75: "3z_T/4", 1.0: "z_T（自成像）"}[frac]
        cfg = build_config(
            name="泰伯效应：z=%.2f·z_T（%s）" % (frac, zname),
            size=N, width=W, wl=wl,
            components=[
                comp("g", "grating", (0, 0, 0),
                     params={"kind": "amplitude_binary", "period": Lam, "duty": 0.5},
                     label="振幅光栅 Λ=100 µm"),
                det_rect("det", z, W * 0.95, W * 0.95, label=zname),
            ],
            sources=[plane_source(wl=wl)],
            category="衍射", difficulty="高级",
            description="周期 100 µm 的二元振幅光栅后自由传播 z=%.2f·z_T（泰伯距离 z_T=2Λ²/λ=%.2f mm）——"
                        "光栅的「自成像」现象（无透镜成像）。" % (frac, zT * 1e3),
            expected=("z=z_T 处光栅**原样复现**（自成像，周期仍是 100 µm）" if frac == 1.0 else
                      "z=z_T/2 处出现**半周期平移**的像（亮暗互换、整体平移 Λ/2）" if frac == 0.5 else
                      "z=z_T/4 处是 2 倍空间频率的中间态" if frac == 0.25 else
                      "分数泰伯位置的非平凡自成像" if frac == 0.75 else "光栅面上的原始方波结构"),
            tags=["泰伯效应", "自成像", "分数泰伯"],
            theory="z_T = 2Λ²/λ；z_T/2 处半周期平移像，z_T/p·q 处出现 p 倍频结构",
            sample=(i < 3),
        )
        reg.add("diff_talbot_%02d.json" % (i + 1), "衍射", "高级",
                cfg["meta"]["zod"]["description"], cfg["meta"]["zod"]["expected"], cfg,
                tags=cfg["meta"]["zod"]["tags"], theory=cfg["meta"]["zod"].get("theory"),
                sample=cfg["meta"].get("sample", False))

    # 不同周期的泰伯距离对比
    for i, Lam in enumerate([50e-6, 80e-6, 150e-6]):
        N, W, wl = 1024, 0.002, WL["633"]
        dx = W / N
        check_sampling(dx, {"周期 Λ": Lam}, 8.0)
        zT = talbot_distance(wl, Lam)
        cfg = build_config(
            name="泰伯距离 ∝ Λ²：Λ=%.0f µm → z_T=%.2f mm" % (Lam * 1e6, zT * 1e3),
            size=N, width=W, wl=wl,
            components=[
                comp("g", "grating", (0, 0, 0),
                     params={"kind": "phase_binary", "period": Lam, "duty": 0.5, "modulation": 2.0},
                     label="相位光栅 Λ=%.0f µm" % (Lam * 1e6)),
                det_rect("det", zT, W * 0.95, W * 0.95, label="z_T 处"),
            ],
            sources=[plane_source(wl=wl)],
            category="衍射", difficulty="高级",
            description="相位光栅（0/π）的泰伯自成像：周期 %.0f µm，z_T=2Λ²/λ=%.2f mm。" % (Lam * 1e6, zT * 1e3),
            expected="相位光栅在 z_T 处形成振幅型的自成像（相位结构转成强度可见），周期仍为 %.0f µm；"
                     "z_T 与 Λ² 成正比——周期缩小一半，泰伯距离缩到 1/4（%.2f mm）。"
                     % (Lam * 1e6, talbot_distance(wl, Lam / 2) * 1e3),
            tags=["泰伯效应", "相位光栅", "Λ² 标度"],
            theory="z_T = 2Λ²/λ，强烈依赖周期",
        )
        reg.add("diff_talbot_scale_%02d.json" % (i + 1), "衍射", "高级",
                cfg["meta"]["zod"]["description"], cfg["meta"]["zod"]["expected"], cfg,
                tags=cfg["meta"]["zod"]["tags"], theory=cfg["meta"]["zod"].get("theory"))


# ==========================================================================
# 直边（刀口）衍射
# ==========================================================================

def knife_edge(reg):
    for i, z in enumerate([0.05, 0.1, 0.2, 0.4, 0.8]):
        N, W, wl = 512, W10, WL["633"]
        dx = W / N
        scale = edge_fringe_scale(wl, z)
        check_sampling(dx, {"条纹尺度": scale}, 3.0)
        NF = fresnel_number(2.0e-3, wl, z)
        cfg = build_config(
            name="直边（刀口）菲涅耳衍射：z=%.0f mm，N_F=%.1f" % (z * 1e3, NF),
            size=N, width=W, wl=wl,
            components=[
                aperture("knife", (0, 0, 0), rectangle(W * 0.98, W * 0.98),
                         x=-W * 0.49, label="半平面挡板（边缘过光轴）"),
                det_rect("det", z, W * 0.95, W * 0.95, label="近场屏"),
            ],
            sources=[plane_source(wl=wl)],
            category="衍射", difficulty="进阶",
            description="用矩形挡板遮住半个视场（边缘位于 x=0）→ 直边衍射（刀口实验）；传播到 z=%.0f mm。" % (z * 1e3),
            expected="几何亮区一侧出现**明暗相间的条纹**，条纹尺度 ~√(λz/2)=%.0f µm；紧贴几何阴影边界处强度为"
                     "入射的 1/4（不是 1/2），阴影区指数衰减但很快趋零。z 越大条纹越粗、越多（菲涅耳数 N_F=%.1f）。"
                     % (scale * 1e6, NF),
            tags=["直边衍射", "刀口", "菲涅耳积分"],
            theory="菲涅耳积分：I(0)=I₀/4 在几何边界上；条纹尺度 √(λz/2)",
            sample=(i < 2),
        )
        reg.add("diff_knife_edge_%02d.json" % (i + 1), "衍射", "进阶",
                cfg["meta"]["zod"]["description"], cfg["meta"]["zod"]["expected"], cfg,
                tags=cfg["meta"]["zod"]["tags"], theory=cfg["meta"]["zod"].get("theory"),
                sample=cfg["meta"].get("sample", False))


# ==========================================================================
# 螺旋相位（光学涡旋）的远场甜甜圈
# ==========================================================================

def vortex_farfield(reg):
    for i, l in enumerate([1, 2, 3, 4]):
        w0, N, W, wl = 1.2e-3, 512, W8, WL["633"]
        dx = W / N
        check_sampling(dx, {"束腰": w0}, 10.0)
        cfg = build_config(
            name="光学涡旋（l=%d）的传播与甜甜圈" % l,
            size=N, width=W, wl=wl,
            components=[
                comp("spp", "spiral_phase", (0, 0, 0), shape=circle(W * 0.45),
                     params={"charge": l}, label="螺旋相位板 l=%d" % l),
                det_circle("det", 0.3, W * 0.45, label="传播 0.3 m 后"),
            ],
            sources=[gauss_source(w0, wl=wl, label="高斯光束 w₀=%.1f mm" % (w0 * 1e3))],
            category="衍射", difficulty="进阶",
            description="拓扑荷 l=%d 的螺旋相位板作用在高斯光束上，传播 0.3 m 观察强度与相位结构。" % l,
            expected="轴上强度**严格为零**（相位奇点），光斑呈甜甜圈，环半径随 l 增大而增大；"
                     "相位图绕轴一圈共 %d 次 2π 跳变（拓扑荷的定义）。l=1 的环最紧凑。" % l,
            tags=["光学涡旋", "OAM", "相位奇点", "甜甜圈"],
            theory="涡旋相位 exp(ilθ)；轴上强度 ∝ r^{2|l|} → 光轴处为零", sample=(i == 0),
        )
        reg.add("diff_vortex_farfield_%02d.json" % (i + 1), "衍射", "进阶",
                cfg["meta"]["zod"]["description"], cfg["meta"]["zod"]["expected"], cfg,
                tags=cfg["meta"]["zod"]["tags"], theory=cfg["meta"]["zod"].get("theory"),
                sample=cfg["meta"].get("sample", False))


# ==========================================================================
# 泽尼克像差 → 衍射受限 PSF
# ==========================================================================

def zernike_psf(reg):
    specs = [("离焦 Z4", {"c4": 0.25}), ("离焦 Z4 半波", {"c4": 0.5}),
             ("像散 Z5", {"c5": 0.3}), ("像散 Z6", {"c6": 0.3}),
             ("彗差 Z7", {"c7": 0.3}), ("彗差 Z8", {"c8": 0.3}),
             ("球差 Z11", {"c11": 0.3}), ("三叶 Z10", {"c10": 0.3})]
    for i, (lbl, coefs) in enumerate(specs):
        wl, f, N, W = WL["633"], 0.3, 512, 0.006
        dx = W / N
        R = 2.0e-3
        check_sampling(dx, {"光瞳半径 R": R}, 20.0)
        r_airy = airy_radius(wl, f, 2 * R)
        check_sampling(dx, {"艾里半径": r_airy}, 4.0)
        params = dict(coefs)
        params["radius"] = R
        cfg = build_config(
            name="像差 PSF：%s（%.2f λ）" % (lbl, list(coefs.values())[0]),
            size=N, width=W, wl=wl,
            components=[
                lens("l1", (0, 0, 0), f=f, aperture=R, label="理想透镜 f=0.3 m"),
                comp("ab", "zernike", (0, 0, 0.02), shape=circle(R), params=params, label=lbl),
                det_circle("det", f, W * 0.45, label="焦面",
                           params={"strehl_aperture": R, "strehl_distance": f}),
            ],
            sources=[plane_source(wl=wl)],
            category="衍射", difficulty="高级",
            description="在理想透镜后加入 %s 像差（Noll 序，系数 %.2f λ，RMS 归一化），观察焦面 PSF 的退化。"
                        % (lbl, list(coefs.values())[0]),
            expected="无像差时是艾里图样（Strehl≈1，第一暗环 %.1f µm）；加入 %s 后峰值下降、图样按像差类型"
                     "变形：离焦→同心环能量外流、像散→焦点分裂成两条正交焦线、彗差→单向彗尾不对称、"
                     "球差→中心亮斑变小但光环变强。Strehl 指标在结果里给出（本实例为 %.2f λ 的像差，"
                     "Strehl 明显低于 1）。"
                     % (r_airy * 1e6, lbl, list(coefs.values())[0]),
            tags=["泽尼克", "像差", "PSF", "Strehl"],
            theory="Strehl ≈ exp(-(2πσ)²)，σ 为波前 RMS（单位 λ）", sample=(i < 2),
        )
        reg.add("diff_zernike_psf_%02d.json" % (i + 1), "衍射", "高级",
                cfg["meta"]["zod"]["description"], cfg["meta"]["zod"]["expected"], cfg,
                tags=cfg["meta"]["zod"]["tags"], theory=cfg["meta"]["zod"].get("theory"),
                sample=cfg["meta"].get("sample", False))


# ==========================================================================
# 4f 傅里叶滤波
# ==========================================================================

def fourf_filter(reg):
    f1 = 0.075
    N, W, wl = 512, W8, WL["633"]
    dx = W / N
    Lam = 200e-6
    x1 = grating_order_pos(1, f1, wl, Lam)
    check_sampling(dx, {"周期 Λ": Lam, "级次位置": x1}, 8.0)
    filters = [
        ("低通（只过零级）", circle(x1 * 0.5), "只让零级通过 → 输出均匀背景，条纹被抹平"),
        ("高通（挡掉零级）", ring(x1 * 0.5, W * 0.4), "挡掉零级 → 只剩 ±1 级干涉，输出是纯正弦条纹（边缘增强）"),
        ("带通（只过一个 ±1 级）", rectangle(x1 * 0.9, W * 0.4), "只留 +1 级 → 输出是倾斜的单一空间频率条纹"),
    ]
    for i, (lbl, fshape, eff) in enumerate(filters):
        comps = [
            comp("g", "grating", (0, 0, 0),
                 params={"kind": "amplitude_sin", "period": Lam, "modulation": 0.8},
                 label="物：正弦振幅光栅 Λ=100 µm"),
            lens("l1", (0, 0, f1), f=f1, aperture=W * 0.45, label="傅里叶透镜 f₁=75 mm"),
            comp("flt", "aperture", (0, 0, 2 * f1), shape=fshape, label="频谱面滤波：%s" % lbl),
            lens("l2", (0, 0, 3 * f1), f=f1, aperture=W * 0.45, label="逆变换透镜 f₂=75 mm"),
            det_rect("det", 4 * f1, W * 0.95, W * 0.95, label="输出面（4f 像面）"),
        ]
        cfg = build_config(
            name="4f 傅里叶滤波系统：%s" % lbl,
            size=N, width=W, wl=wl, components=comps,
            sources=[plane_source(wl=wl)],
            category="衍射", difficulty="高级",
            description="经典 4f 系统：物（正弦光栅）→ 傅里叶透镜 → **频谱面插入滤波器** → 逆变换透镜 → 像面；"
                        "滤波器规格 %s（±1 级位于 x=±%.2f mm）。" % (lbl, x1 * 1e3),
            expected="频谱面上的三个光点（0、±1 级）被滤波器裁剪后，像面重建出对应的空间频率成分：%s。"
                     "这正是「4f 相干光学图像处理」的物理内核——滤波即空间频率选择。" % eff,
            tags=["4f 系统", "傅里叶光学", "空间滤波"],
            theory="4f：物面→频谱面→像面，滤波即乘传递函数",
            sample=(i < 2),
        )
        reg.add("diff_4f_filter_%02d.json" % (i + 1), "衍射", "高级",
                cfg["meta"]["zod"]["description"], cfg["meta"]["zod"]["expected"], cfg,
                tags=cfg["meta"]["zod"]["tags"], theory=cfg["meta"]["zod"].get("theory"),
                sample=cfg["meta"].get("sample", False))


# ==========================================================================
# 互补/双孔等补充结构
# ==========================================================================

def diffraction_extras(reg):
    # 双孔（杨氏干涉的二维版）：两个圆孔
    for i, (d, r) in enumerate([(1.5e-3, 0.2e-3), (2.5e-3, 0.25e-3)]):
        wl, z, N, W = WL["633"], 1.0, 1024, W20
        dx = W / N
        check_sampling(dx, {"孔距 d": d, "孔半径 r": r}, 5.0)
        fring = fringe_spacing(wl, z, d)
        check_sampling(dx, {"条纹间距": fring}, 5.0)
        comps = [
            aperture("h1", (0, 0, 0), circle(r), parallel=True, x=-d / 2, label="孔 1"),
            aperture("h2", (0, 0, 0), circle(r), parallel=True, x=d / 2, label="孔 2"),
            det_rect("det", z, W * 0.95, W * 0.95, label="干涉屏"),
        ]
        cfg = build_config(
            name="双孔干涉（圆孔版杨氏）：孔径 %.2f mm、孔距 %.1f mm" % (r * 1e3, d * 1e3),
            size=N, width=W, wl=wl, components=comps,
            sources=[plane_source(wl=wl)],
            category="衍射", difficulty="进阶",
            description="两个半径 %.2f mm 的圆孔（中心距 %.1f mm）替代狭缝，观察二维干涉图样（用同 z 的"
                        "parallel 并集实现）。" % (r * 1e3, d * 1e3),
            expected="干涉条纹间距 Δx=λz/d=%.2f mm，条纹方向垂直于两孔连线；整体被圆孔的艾里包络调制"
                     "（艾里第一暗环 %.1f mm）。孔越大条纹越少、孔距越大条纹越密。" % (fring * 1e3, airy_radius(wl, z, 2 * r) * 1e3),
            tags=["双孔", "杨氏实验", "二维干涉"],
            theory="Δx = λz/d，包络为圆孔艾里图样", sample=(i == 0),
        )
        reg.add("diff_two_hole_%02d.json" % (i + 1), "衍射", "进阶",
                cfg["meta"]["zod"]["description"], cfg["meta"]["zod"]["expected"], cfg,
                tags=cfg["meta"]["zod"]["tags"], theory=cfg["meta"]["zod"].get("theory"),
                sample=cfg["meta"].get("sample", False))

    # 方孔阵列（二维光栅）
    for i, n in enumerate([3, 4]):
        wl, z, N, W = WL["633"], 1.0, 1024, W20
        d, r = 2.5e-3, 0.4e-3
        dx = W / N
        check_sampling(dx, {"孔半径": r, "孔距": d}, 10.0)
        comps = []
        k = 0
        for ix in range(n):
            for iy in range(n):
                comps.append(aperture("h%d" % k, (0, 0, 0), circle(r), parallel=True,
                                      x=(ix - (n - 1) / 2.0) * d,
                                      y=(iy - (n - 1) / 2.0) * d))
                k += 1
        comps.append(det_rect("det", z, W * 0.95, W * 0.95, label="二维衍射图"))
        cfg = build_config(
            name="%d×%d 圆孔阵列的二维衍射（二维光栅）" % (n, n),
            size=N, width=W, wl=wl, components=comps,
            sources=[plane_source(wl=wl)],
            category="衍射", difficulty="高级",
            description="%d×%d 个半径 0.4 mm 的圆孔按 2.5 mm 的二维周期排列（parallel 并集），"
                        "远场 %.0f m 观察二维衍射。"
                        % (n, n, z),
            expected="二维点阵图样：主极大出现在 (m·λz/d, n·λz/d) 的二维格点上（间距 %.2f mm），"
                     "每个主极大的包络是圆孔艾里图样的二维复制；阵列越大主极大越锐。" % (fringe_spacing(wl, z, d) * 1e3),
            tags=["二维光栅", "孔阵列", "点阵衍射"],
            theory="二维周期 a 的阵列 → 远场为倒格点阵，间隔 λz/a", sample=(i == 0),
        )
        reg.add("diff_hole_array_%02d.json" % (i + 1), "衍射", "高级",
                cfg["meta"]["zod"]["description"], cfg["meta"]["zod"]["expected"], cfg,
                tags=cfg["meta"]["zod"]["tags"], theory=cfg["meta"]["zod"].get("theory"),
                sample=cfg["meta"].get("sample", False))

    # 边缘软化（edge_sigma）消除硬边混叠
    for i, (es, lbl) in enumerate([(0.0, "理想硬边 edge_sigma=0"),
                                   (2e-6, "软化 2 µm"),
                                   (5e-6, "软化 5 µm")]):
        wl, z, N, W = WL["633"], 1.0, 1024, W20
        a = 0.3e-3
        dx = W / N
        check_sampling(dx, {"缝宽": a}, 8.0)
        cfg = build_config(
            name="硬边 vs 软化边缘：%s" % lbl,
            size=N, width=W, wl=wl,
            components=[aperture("slit", (0, 0, 0), rectangle(a, W * 0.95),
                                 edge_sigma=(es if es else None), label="单缝 a=0.3 mm（%s）" % lbl),
                        det_rect("det", z, W * 0.95, W * 0.95, label="远场")],
            sources=[plane_source(wl=wl)],
            category="衍射", difficulty="高级",
            description="缝宽 0.3 mm 的单缝，边缘从理想硬边改为 erf 软化（edge_sigma=%s）：这是**数值正则化**"
                        "手段，用来抑制硬边全带宽频谱带来的混叠伪影。" % lbl,
            expected="主极大与第一暗点（λz/a=%.2f mm）在两种情况下位置相同；硬边版本的高频副极大里混有"
                     "网格混叠造成的「毛刺」与不对称，软化后副极大更干净、动态范围更平滑。"
                     "软化尺度必须远小于特征尺寸（这里 %.1f µm ≪ 300 µm）才不改变物理。"
                     % (sinc_first_zero(wl, z, a) * 1e3, es * 1e6),
            tags=["edge_sigma", "混叠抑制", "数值正则化"],
            theory="erf 软化边缘等效于乘一个高斯切趾，抑制高频拖尾", sample=(i == 1),
        )
        reg.add("diff_edge_soft_%02d.json" % (i + 1), "衍射", "高级",
                cfg["meta"]["zod"]["description"], cfg["meta"]["zod"]["expected"], cfg,
                tags=cfg["meta"]["zod"]["tags"], theory=cfg["meta"]["zod"].get("theory"),
                sample=cfg["meta"].get("sample", False))

    # 环形孔阵列（光子筛风格的近似）
    for i, (nring, rout) in enumerate([(4, 1.0e-3), (6, 1.2e-3)]):
        wl, z, N, W = WL["633"], 1.0, 1024, W20
        dx = W / N
        dr = rout / nring
        check_sampling(dx, {"环宽": dr / 2}, 4.0)
        cfg = build_config(
            name="同心环阵列（光子筛风格）：%d 环" % nring,
            size=N, width=W, wl=wl,
            components=[
                aperture("ap", (0, 0, 0), ring(rout - dr, rout), label="外环"),
                det_rect("det", z, W * 0.95, W * 0.95, label="远场"),
            ],
            sources=[plane_source(wl=wl)],
            category="衍射", difficulty="高级",
            description="半径端到端的单环孔（外半径 %.1f mm、环宽 %.0f µm）的远场；串接多个同心环即是"
                        "「光子筛/环阵列」的简化模型。" % (rout * 1e3, dr * 1e6),
            expected="环孔远场是 J₀ 型中心亮斑 + 同心环；环越窄中心斑越大（类似圆孔但缺少低频），"
                     "第一暗环尺度 ≈1.22λz/(2r_out)=%.2f mm。与实心圆孔相比中心峰的旁瓣更强。"
                     % (airy_radius(wl, z, 2 * rout) * 1e3),
            tags=["环形孔", "光子筛", "J₀"],
            theory="环孔远场 ∝ 组合贝塞尔函数，环宽是新的特征尺度",
        )
        reg.add("diff_ring_array_%02d.json" % (i + 1), "衍射", "高级",
                cfg["meta"]["zod"]["description"], cfg["meta"]["zod"]["expected"], cfg,
                tags=cfg["meta"]["zod"]["tags"], theory=cfg["meta"]["zod"].get("theory"))
