"""仪器类示例（inst_*）—— 目标 60+。

迈克尔逊族、马赫-曾德尔族、法布里-珀罗/平行平板、萨格纳克环、4f 图像处理、
涡旋日冕仪、光栅光谱仪、波带片显微镜、望远镜、扩束器、空间滤波器、相衬显微镜、
共聚焦针孔、偏振 DIC、超表面透镜仪器、成像系统。
"""

import math

from wavegen.common import *   # noqa: F401,F403
from wavegen.fam_interf import mz_comps, michelson_comps, d_port_z, d_port_nx
from wavegen.fam_pol import _pol, _ret, pol_stage

W10 = 0.01
W8 = 0.008
W20 = 0.02
PI = math.pi


# ==========================================================================
# 1. 迈克尔逊族
# ==========================================================================

def michelson_family(reg):
    specs = [
        ("等臂（输出端口全暗）", 0.0, 0.0, "回光源端口应为全暗、侧端口全亮——等臂时两臂反相相消"),
        ("臂长差 λ/4（相位差 π/2）", 0.25, 0.0, "两端口各得一半功率"),
        ("臂长差 λ/2（相位差 π）", 0.5, 0.0, "两端口反转：回光源端口变全亮"),
        ("臂长差 λ（相位差 2π）", 1.0, 0.0, "回到等臂状态，端口恢复原始分配"),
        ("臂长差 10λ（大量级）", 10.0, 0.0, "相位差 20π 等价于 0——干涉仪只对 2kΔL mod 2π 敏感"),
        ("等臂 + 微倾 0.1 mrad", 0.0, 1e-4, "等臂直条纹：间距 λ/(2θ)=3.16 mm"),
        ("微倾 + λ/4 臂长差", 0.25, 1e-4, "条纹对比度不变、整体相位平移（条纹移动 1/4 个周期）"),
    ]
    for i, (lbl, frac, tilt, eff) in enumerate(specs):
        wl, W = WL["633"], W10
        d = frac * wl / 2.0
        N = 512
        cfg = build_config(
            name="迈克尔逊干涉仪：%s" % lbl,
            size=N, width=W, wl=wl,
            components=[
                beam_splitter("bs", (0, 0, 0), reflectivity=0.5, aperture=W * 0.35,
                              label="分束器 50:50"),
                mirror("m1", (0, 0, 0.05), reflectivity=1.0, yaw=0.0,
                       shape=circle(W * 0.35), label="端镜 A（固定，回射）"),
                mirror("m2", (0.05 + d, 0, 0), reflectivity=1.0, yaw=PI / 2,
                       tilt_y=tilt, shape=circle(W * 0.35),
                       label="端镜 B（%s）" % ("微倾 %.1f mrad" % (tilt * 1e3) if tilt else "平面回射")),
                d_port_z("d1", -0.05, "回光源端口", W * 0.35),
                d_port_nx("d2", -0.05, "侧端口", W * 0.35),
            ],
            sources=[plane_source(wl=wl, label="HeNe 632.8 nm 平面波")],
            category="仪器", difficulty="进阶",
            description="迈克尔逊干涉仪（分束器 + 两个回射端镜 + 两个输出端口）：%s。"
                        "端镜 B 沿光轴移动 ΔL=%.0f nm。" % (lbl, d * 1e9),
            expected="%s。折返光路的相位按 2kL 累加（往返 = 两倍臂长），相位差 = 2kΔL = %.0f°。"
                     "两个端口严格互补，功率之和 = 输入功率（无吸收元件）。" % (eff, 4 * PI * frac * 0.5 / PI * 180 / 2 * 2),
            tags=["迈克尔逊", "分束器", "互补端口", "折返相位"],
            theory="相位差 2kΔL；BS 约定 t=√(1−R)、r=i√R", sample=(i < 4),
        )
        reg.add("inst_michelson_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))

    # 多波长「白光」干涉
    for i, nwl in enumerate([2, 3, 5]):
        wl0, W = WL["633"], W10
        d = wl0 / 4.0
        cols = [WL[k] for k in (["633", "532"] if nwl == 2 else
                               (["633", "532", "450"] if nwl == 3 else
                                ["633", "550", "500", "450", "405"]))]
        srcs = [plane_source(wl=c, power=1e-3 / nwl, sid="s%d" % k, group="g%d" % k,
                             label="λ=%.0f nm" % (c * 1e9)) for k, c in enumerate(cols)]
        cfg = build_config(
            name="迈克尔逊多波长干涉：%d 个波长同时入射" % nwl,
            size=512, width=W, wl=wl0, components=[
                beam_splitter("bs", (0, 0, 0), reflectivity=0.5, aperture=W * 0.35,
                              label="分束器 50:50"),
                mirror("m1", (0, 0, 0.05), reflectivity=1.0, yaw=0.0,
                       shape=circle(W * 0.35), label="端镜 A"),
                mirror("m2", (0.05 + d, 0, 0), reflectivity=1.0, yaw=PI / 2,
                       tilt_y=1.5e-4, shape=circle(W * 0.35), label="端镜 B（微倾）"),
                d_port_z("d1", -0.05, "白光条纹端口", W * 0.35),
            ],
            sources=srcs,
            category="仪器", difficulty="高级",
            description="迈克尔逊干涉仪同时输入 %d 个不同波长（各自独立相干组）：波长不同 → "
                        "同一臂长差对应的相位不同 → 各色条纹错开叠加。" % nwl,
            expected="不同波长的条纹间距与相位都不同（相位差 2kΔL ∝ 1/λ），叠加后只有在光程差接近 0 处"
                     "各色同时接近极大 → 出现「白/亮」的中心条纹，两侧被彩色条纹占据（相干长度决定"
                     "能看清的条纹数）。这正是白光干涉测零光程差、以及傅里叶变换光谱仪的原理。",
            tags=["白光干涉", "多波长", "相干长度", "FTIR"],
            theory="各波长独立传播与叠加（不同相干组）；相位差 2kΔL 随 λ 变化", sample=(i == 0),
        )
        reg.add("inst_michelson_multi_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))


# ==========================================================================
# 2. 马赫-曾德尔族（不平衡臂 + 相位板扫描）
# ==========================================================================

def mach_zehnder_family(reg):
    for i, (lbl, R, nph) in enumerate([
            ("平衡 50:50，无相位板", 0.5, 0.0),
            ("平衡 50:50，相位板 π/2", 0.5, 0.25),
            ("平衡 50:50，相位板 π（暗口）", 0.5, 0.5),
            ("不对称分束 70:30", 0.7, 0.0),
            ("不对称分束 90:10", 0.9, 0.0),
            ("平衡 + 一臂插入楔形（直条纹）", 0.5, -1.0)]):
        wl, W = WL["633"], W10
        A = 0.08
        arm_extra = []
        if nph >= 0 and nph > 0:
            d_plate = 0.01
            arm_extra = [comp("plate", "medium", (0, 0, 0.04), shape=circle(W * 0.4),
                              params={"distance": d_plate,
                                      "index": 1.0 + nph * wl / d_plate,
                                      "absorption": 0.0, "steps": 4},
                              label="相位板（%.2f×2π）" % nph)]
        elif nph < 0:
            arm_extra = [comp("wedge", "wedge", (0, 0, 0.04), shape=circle(W * 0.4),
                              params={"alpha": 0.01, "index": 1.5},
                              label="楔形棱镜 α=10 mrad（等厚直条纹）")]
        comps = mz_comps(W, A=A, bs_R=R, arm_extra=arm_extra, d1_label="端口 1", d2_label="端口 2")
        cfg = build_config(
            name="马赫-曾德尔干涉仪：%s" % lbl,
            size=512, width=W, wl=wl, components=comps,
            sources=[plane_source(wl=wl, label="平面波 633 nm")],
            category="仪器", difficulty="进阶",
            description="马赫-曾德尔干涉仪（两个 45° 折转镜 + 两个分束器）：%s。" % lbl,
            expected=("平衡时端口 1 全亮、端口 2 全暗（两臂相消）；两端口之和 = 输入功率。" if nph == 0.0 and R == 0.5 else
                      "端口分配随相位 %.2f×2π 按 (1±cosφ)/2 变化，两端口互补。" % nph if nph >= 0 else
                      "一臂插入楔形棱镜 → 线性相位 → 输出为等厚直条纹（间距 λ/((n−1)α)=126.6 µm），"
                      "这是「空间相位斜坡」的最直观演示。"),
            tags=["马赫-曾德尔", "分束器", "相位扫描", "能量守恒"],
            theory="BS 约定 t=√(1−R)、r=i√R；两端口互补", sample=(i < 4),
        )
        reg.add("inst_mach_zehnder_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))


# ==========================================================================
# 3. 法布里-珀罗 / 平行平板（两光束近似）
# ==========================================================================

def fabry_perot(reg):
    for i, (gap, R, tilt) in enumerate([(20e-6, 0.3, 1.0e-4), (20e-6, 0.5, 1.0e-4),
                                        (20e-6, 0.8, 1.0e-4), (40e-6, 0.5, 1.0e-4),
                                        (10e-6, 0.5, 2.0e-4), (20e-6, 0.5, 2.0e-4)]):
        wl, W = WL["633"], W10
        spacing = wl / (2 * tilt)
        N = pick_grid(spacing, W, minimum=4.0)
        cfg = build_config(
            name="法布里-珀罗/平行平板：间隙 %.0f µm、R=%.1f" % (gap * 1e6, R),
            size=N, width=W, wl=wl,
            components=[
                beam_splitter("bs", (0, 0, 0), reflectivity=0.5, aperture=W * 0.35,
                              label="分束器 50:50"),
                mirror("mflat", (0, 0, 0.05), reflectivity=1.0, yaw=0.0,
                       shape=circle(W * 0.35), label="参考臂端镜"),
                mirror("mpart", (0.05, 0, 0), reflectivity=math.sqrt(R), yaw=PI / 2,
                       shape=circle(W * 0.35),
                       label="前镜（功率反射率 R=%.1f）" % R),
                mirror("mback", (0.05 + gap, 0, 0), reflectivity=1.0, yaw=PI / 2,
                       tilt_y=tilt, shape=circle(W * 0.35),
                       label="后镜（间隙 %.0f µm、微倾 %.1f mrad）" % (gap * 1e6, tilt * 1e3)),
                d_port_z("d1", -0.05, "反射端口（干涉条纹）", W * 0.35),
            ],
            sources=[plane_source(wl=wl)],
            category="仪器", difficulty="高级",
            description="由「前镜部分反射 + 后镜全反」构成的平行平板（法布里-珀罗腔的两光束模型）："
                        "间隙 %.0f µm、前镜功率反射率 R=%.1f、后镜微倾 %.1f mrad。" % (gap * 1e6, R, tilt * 1e3),
            expected="返回两束光的光程差 2ngap=%.1fλ（间隙每变 λ/2 相位走一圈 2π），叠加后是间距 "
                     "λ/(2θ)=%.2f mm 的直条纹。理想腔的精细度 F=π√R/(1−R)=%.1f、自由光谱范围 "
                     "Δλ=λ²/(2nL)=%.3f nm、峰值透射=1。**注意本内核只计首次往返**（两光束近似），"
                     "高反射率时 Airy 多光束的窄峰不在结果中（会给出 scene_cycle_dropped 告警）。"
                     % (2 * gap / wl, spacing * 1e3, fp_finesse(R), fp_fsr(wl, gap) * 1e9),
            tags=["法布里-珀罗", "平行平板", "精细度", "自由光谱范围"],
            theory="Δλ_FSR=λ²/(2nL)，F=π√R/(1−R)，T_max=1（无吸收）；两光束近似", sample=(i < 2),
        )
        reg.add("inst_fabry_perot_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))


# ==========================================================================
# 4. 萨格纳克环
# ==========================================================================

def _sagnac_loop(L=0.1, r=0.005):
    """方形萨格纳克环：BS(0,0,0) → m1(0,0,L) → m2(L,0,L) → m3(L,0,0) → 回到 BS。"""
    return [
        beam_splitter("bs", (0, 0, 0), reflectivity=0.5, aperture=r, label="分束器 50:50"),
        mirror("m1", (0, 0, L), reflectivity=1.0, yaw=-PI / 4, shape=circle(r),
               label="环镜 1（+z→+x）"),
        mirror("m2", (L, 0, L), reflectivity=1.0, yaw=-3 * PI / 4, shape=circle(r),
               label="环镜 2（+x→−z）"),
        mirror("m3", (L, 0, 0), reflectivity=1.0, yaw=-PI / 4, shape=circle(r),
               label="环镜 3（−z→−x，回到分束器）"),
    ]


def sagnac_ring(reg):
    specs = [("三角形/方形环：互易输出为零", 0.5, 0.0),
             ("环内插入双折射元件（非互易）", 0.5, 1.0),
             ("分束比 70:30 的环", 0.7, 0.0),
             ("分束比 30:70 的环", 0.3, 0.0),
             ("带楔形相位的环", 0.5, 2.0)]
    for i, (lbl, R, kind) in enumerate(specs):
        wl, W = WL["633"], 0.04
        L = 0.1
        comps = _sagnac_loop(L)
        comps[0] = beam_splitter("bs", (0, 0, 0), reflectivity=R, aperture=W * 0.15,
                                 label="分束器 %.0f:%.0f" % (R * 100, (1 - R) * 100))
        if kind == 1.0:
            comps.append(comp("xtal", "uniaxial", (0, 0, L / 2), shape=circle(W * 0.15),
                              params={"distance": 0.005, "n_o": 1.5, "n_e": 1.6},
                              label="环内单轴晶体（打破互易对称）"))
        elif kind == 2.0:
            comps.append(comp("wedge", "wedge", (0, 0, L / 2), shape=circle(W * 0.15),
                              params={"alpha": 0.002, "index": 1.5},
                              label="环内楔形（相位斜坡）"))
        comps.append(detector("dout", (-0.06, 0, 0), yaw=PI / 2, shape=circle(W * 0.15),
                              label="输出端口"))
        comps.append(detector("dback", (0, 0, -0.05), yaw=0.0, shape=circle(W * 0.15),
                              label="回光源端口"))
        cfg = build_config(
            name="萨格纳克环干涉仪：%s" % lbl,
            size=512, width=W, wl=wl, components=comps,
            sources=[plane_source(wl=wl, power=1e-3, label="平面波 633 nm")],
            polarized=(kind == 1.0),
            category="仪器", difficulty="高级",
            description="方形萨格纳克环（%.0f×%.0f mm，3 个环镜 + 1 个分束器）：光从分束器分成顺时针与"
                        "逆时针两束，沿同一环路反向传播后在分束器处重新合束。%s" % (L * 1e3, L * 1e3, lbl),
            expected=("互易光路的萨格纳克输出端口**严格为零**：两束反向光走过的光程完全相同 → "
                      "相位差恒为 0 → 在输出端口相消、功率全部回到入射端口（这是萨格纳克环抗击振动/"
                      "温度扰动的根本原因，也是它被用作反射镜的原因）。" if kind == 0.0 else
                      "环内的非互易元件（双折射/楔形）使顺逆两束的偏振或相位演化不同 → "
                      "输出端口重新出现信号（非零功率）。真实萨格纳克陀螺正是用旋转引入的光程差 "
                      "Δφ=8πNA/(λc)·Ω，本内核没有旋转参考系，因此只能用「非互易元件」模拟这一效应。"),
            tags=["萨格纳克", "环形干涉仪", "互易性", "非互易"],
            theory="互易环输出为零；旋转引入 Δφ = 8πNAΩ/(λc)", sample=(i < 2),
        )
        reg.add("inst_sagnac_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))


# ==========================================================================
# 5. 4f 图像处理
# ==========================================================================

def fourf_processing(reg):
    f1 = 0.075
    specs = [("低通滤波（挡高频）", circle(0.0003), "去掉细结构、保留大尺度轮廓"),
             ("高通滤波（挡零级）", ring(0.0003, W8 * 0.4), "边缘增强（只留变化剧烈的细节）"),
             ("方向滤波（竖直狭缝）", rectangle(0.0002, W8 * 0.4), "只保留水平方向的空间频率（条纹变竖直）"),
             ("方向滤波（水平狭缝）", rectangle(W8 * 0.4, 0.0002), "只保留竖直方向的空间频率"),
             ("相位对比（λ/4 相位环）", None, "零级加 π/2 相位 → 相位物体变成强度可见")]
    for i, (lbl, fshape, eff) in enumerate(specs):
        wl, W = WL["633"], W8
        Lam = 200e-6
        x1 = grating_order_pos(1, f1, wl, Lam)
        N = 512
        base = [aperture("obj", (0, 0, 0), rectangle(W * 0.5, W * 0.5),
                         label="物：方孔（含各种空间频率）"),
                lens("l1", (0, 0, f1), f=f1, aperture=W * 0.49, label="傅里叶透镜 f=75 mm")]
        if fshape is not None:
            base.append(comp("flt", "aperture", (0, 0, 2 * f1), shape=fshape, label="频谱面滤波：%s" % lbl))
        else:
            base.append(comp("pc", "medium", (0, 0, 2 * f1), shape=circle(0.0012),
                             params={"distance": 0.005, "index": 1.0 + wl / (4 * 0.005),
                                     "absorption": 0.0, "steps": 2},
                             label="相位环（λ/4，只作用零级）"))
            base.append(comp("stop", "aperture", (0, 0, 2 * f1), shape=ring(0.0012, W8 * 0.49),
                             label="环形挡板（挡掉零级以外的区域）"))
        base.append(lens("l2", (0, 0, 3 * f1), f=f1, aperture=W * 0.49, label="逆变换透镜 f=75 mm"))
        base.append(det_rect("det", 4 * f1, W * 0.95, W * 0.95, label="像面"))
        cfg = build_config(
            name="4f 相干图像处理：%s" % lbl,
            size=N, width=W, wl=wl, components=base,
            sources=[plane_source(wl=wl)],
            category="仪器", difficulty="高级",
            description="4f 系统（物 → 傅里叶透镜 → **频谱面处理** → 逆变换透镜 → 像）：%s。" % lbl,
            expected="%s。物是 8 mm 方孔（频谱为二维 sinc，零级最大、高频向四周展开），"
                     "在频谱面上做空间频率选择后，像面重建出对应的滤波结果——"
                     "这就是 4f 相干光学图像处理（阿贝成像理论）的全部内容。" % eff,
            tags=["4f系统", "傅里叶光学", "图像处理", "阿贝成像"],
            theory="4f：物面 → 频谱面（傅里叶变换）→ 像面（逆变换）；滤波=乘传递函数", sample=(i < 2),
        )
        reg.add("inst_4f_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))


# ==========================================================================
# 6. 涡旋日冕仪
# ==========================================================================

def vortex_coronagraph(reg):
    for i, (l, Lyot) in enumerate([(1, 0.6), (2, 0.6), (1, 0.9), (3, 0.9), (1, 0.3)]):
        wl, W, f = WL["633"], W8, 0.075
        N = 512
        comps = [
            comp("spp", "spiral_phase", (0, 0, 0), shape=circle(W * 0.4),
                 params={"charge": l}, label="涡旋相位板（拓扑荷 %d）" % l),
            lens("l1", (0, 0, f), f=f, aperture=W * 0.49, label="傅里叶透镜 f=75 mm"),
            comp("lyot", "aperture", (0, 0, 2 * f), shape=circle(W * 0.49 * Lyot),
                 label="Lyot 光阑（透过率 %.0f%% 直径）" % (Lyot * 100)),
            lens("l2", (0, 0, 3 * f), f=f, aperture=W * 0.49, label="逆变换透镜"),
            det_rect("det", 4 * f, W * 0.95, W * 0.95, label="最终像面"),
        ]
        cfg = build_config(
            name="涡旋日冕仪：拓扑荷 l=%d、Lyot 光阑 %.0f%%" % (l, Lyot * 100),
            size=N, width=W, wl=wl, components=comps,
            sources=[plane_source(wl=wl, label="轴上恒星（平面波）")],
            category="仪器", difficulty="高级",
            description="涡旋日冕仪：轴上星光经拓扑荷 %d 的涡旋相位板后在傅里叶面上形成**甜甜圈**"
                        "（轴上严格为零），再用 Lyot 光阑抑制外圈，最后成像。" % l,
            expected="轴上恒星的像在中心被**抑制**（涡旋把中心能量搬到环上，环被 Lyot 光阑挡住）；"
                     "对照偏移的视场点（另见 vortex_coronagraph 的偏心实例）会保留。"
                     "这是高对比度成像（寻找系外行星）的核心器件：理想涡旋日冕仪在 l=2 时"
                     "对轴上光的抑制度最高，且对低阶像差比传统光阑日冕仪宽容。",
            tags=["涡旋日冕仪", "高对比度成像", "Lyot光阑", "系外行星"],
            theory="涡旋相位 exp(ilθ) 使轴上光在焦面形成暗核；Lyot 光阑挡掉环状能量", sample=(i < 2),
        )
        reg.add("inst_vortex_coronagraph_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))

    # 偏心视场点：轴上被抑制、轴外保留
    for i, (dx, l) in enumerate([(0.0, 2), (2e-4, 2), (5e-4, 2)]):
        wl, W, f = WL["633"], W8, 0.075
        objs = [
            comp("spp", "spiral_phase", (0, 0, 0), shape=circle(W * 0.4),
                 params={"charge": l}, label="涡旋相位板 l=%d" % l),
            lens("l1", (0, 0, f), f=f, aperture=W * 0.49, label="傅里叶透镜"),
            comp("lyot", "aperture", (0, 0, 2 * f), shape=circle(W * 0.49 * 0.6),
                 label="Lyot 光阑 60%"),
            lens("l2", (0, 0, 3 * f), f=f, aperture=W * 0.49, label="逆变换透镜"),
            det_rect("det", 4 * f, W * 0.95, W * 0.95, label="像面"),
        ]
        srcs = [plane_source(wl=wl, tilt_x=dx, label="视场点（倾角 %.2f mrad）" % (dx * 1e3))]
        cfg = build_config(
            name="涡旋日冕仪视场点：倾角 %.2f mrad" % (dx * 1e3),
            size=512, width=W, wl=wl, components=objs, sources=srcs,
            category="仪器", difficulty="高级",
            description="涡旋日冕仪的离轴视场点（倾角 %.2f mrad）：星光的焦点偏离涡旋奇点，"
                        "抑制效果随之下降。" % (dx * 1e3),
            expected=("轴上（0 mrad）：像被强烈抑制（中心近零）。" if dx == 0 else
                      "倾角 %.2f mrad 时焦点离开奇点，涡旋不再把光全部搬进环里 → "
                      "残余光穿透 Lyot 光阑、像面出现明显的残差（对比度下降）。"
                      "因此「日冕仪的有效视场」由奇点对准精度决定。" % (dx * 1e3)),
            tags=["涡旋日冕仪", "视场", "对准误差"],
            theory="轴外点偏离奇点 → 抑制比下降，有效视场受限", sample=(dx == 0),
        )
        reg.add("inst_vortex_coronagraph_fov_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))


# ==========================================================================
# 7. 光栅光谱仪
# ==========================================================================

def grating_spectrometer(reg):
    wls = ["405", "450", "500", "550", "600", "650", "700"]
    f = 0.3
    Lam = 100e-6
    for i, key in enumerate(wls):
        wl, N, W = WL[key], 1024, W8
        x1 = grating_order_pos(1, f, wl, Lam)
        check_sampling(W / N, {"±1 级位置": x1}, 20.0)
        cfg = build_config(
            name="光栅光谱仪：λ=%s nm 的级次位置" % key,
            size=N, width=W, wl=wl,
            components=[
                aperture("es", (0, 0, 0), rectangle(0.5e-3, W8 * 0.9), label="入射狭缝 0.5 mm"),
                comp("g", "grating", (0, 0, 0), params={"kind": "blazed", "period": Lam,
                                                        "blaze_depth": 2 * PI},
                     label="闪耀光栅 Λ=100 µm（能量全进 +1 级）"),
                cyl_lens("cl", (0, 0, 0.02), f=f, rotation=0.0,
                         label="柱面傅里叶透镜 f=0.3 m（只沿色散方向聚焦）"),
                det_rect("det", 0.32, W8 * 0.95, W8 * 0.95, label="焦面（谱线）"),
            ],
            sources=[plane_source(wl=wl, label="λ=%s nm" % key)],
            category="仪器", difficulty="进阶",
            description="单色仪/光谱仪的核心几何：入射狭缝 → 闪耀光栅（周期 100 µm，闪耀深度 2π 把能量"
                        "全部送进 +1 级）→ 柱面傅里叶透镜 → 焦面 CCD。" % (),
            expected="谱线位置 x = f·λ/Λ = %.3f mm（λ=%s nm）；色散率 dλ/dx = Λ/f = %.2f nm/mm。"
                     "闪耀光栅使 +1 级拿到几乎全部功率（±2 级与 0 级都很暗）；柱面透镜让谱线在"
                     "垂直方向保持狭缝的高度而不是收缩成点。改变波长谱线位置线性移动 → 这就是摄谱仪。"
                     % (x1 * 1e3, key, Lam / f * 1e9),
            tags=["光谱仪", "闪耀光栅", "色散", "谱线"],
            theory="x = fλ/Λ；色散率 dλ/dx = Λ/f", sample=(i == 0),
        )
        reg.add("inst_spectrometer_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))

    # 多波长同时（彩色谱）
    for i, (n, wl_list) in enumerate([(3, ["450", "550", "650"]),
                                      (5, ["450", "500", "550", "600", "650"]),
                                      (7, wls)]):
        N, W = 1024, W8
        f = 0.3
        srcs = [plane_source(wl=WL[k], power=1e-3 / n, sid="s%s" % k, group="g%s" % k,
                             label="λ=%s nm" % k) for k in wl_list]
        cfg = build_config(
            name="多波长光谱仪：%d 条谱线同时成像" % n,
            size=N, width=W, wl=WL["550"],
            components=[
                aperture("es", (0, 0, 0), rectangle(0.5e-3, W8 * 0.9), label="入射狭缝"),
                comp("g", "grating", (0, 0, 0), params={"kind": "blazed", "period": 100e-6,
                                                        "blaze_depth": 2 * PI},
                     label="闪耀光栅 Λ=100 µm"),
                cyl_lens("cl", (0, 0, 0.02), f=f, rotation=0.0, label="柱面傅里叶透镜 f=0.3 m"),
                det_rect("det", 0.32, W8 * 0.95, W8 * 0.95, label="焦面（多条谱线）"),
            ],
            sources=srcs,
            category="仪器", difficulty="高级",
            description="%d 个波长同时进入光栅光谱仪（各自独立相干组）——焦面上同时出现 %d 条谱线，"
                        "位置 x=fλ/Λ 与波长成正比。" % (n, n),
            expected="谱线按 x=fλ/Λ 线性排开：从 450 nm 的 %.2f mm 到 650 nm 的 %.2f mm，"
                     "相邻谱线间隔 ≈f·Δλ/Λ。用 GUI 的「图像视图」(field=color) 可以直接看到"
                     "按波长着色的彩色谱线——这正是光谱仪把波长变成位置的过程。"
                     % (WL["450"] * f / 100e-6 * 1e3, WL["650"] * f / 100e-6 * 1e3),
            tags=["光谱仪", "多波长", "彩色谱", "色散"],
            theory="x = fλ/Λ 线性色散；多波长独立叠加", sample=(i == 0),
        )
        reg.add("inst_spectrometer_multi_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))


# ==========================================================================
# 8. 波带片显微镜 / 光子筛
# ==========================================================================

def zoneplate_microscope(reg):
    for i, (f, kind, key) in enumerate([(0.05, "phase", "633"), (0.03, "phase", "450"),
                                        (0.08, "phase", "633"), (0.05, "amplitude", "633"),
                                        (0.05, "phase", "780")]):
        wl, N, W = WL[key], 1024, 0.004
        dx = W / N
        r_out = min(zone_radius(8, wl, f), W * 0.45)
        m_act = zone_count(r_out, wl, f)
        dr = wl * f / (2 * r_out)
        check_sampling(dx, {"外环宽度": dr}, 3.0)
        spot = wl * f / (2 * r_out)
        cfg = build_config(
            name="波带片显微镜：%s f=%.0f mm、λ=%s nm" % ("相位型" if kind == "phase" else "振幅型", f * 1e3, key),
            size=N, width=W, wl=wl,
            components=[
                comp("zp", "zone_plate", (0, 0, 0), shape=circle(r_out),
                     params={"f": f, "radius": r_out, "kind": kind},
                     label="波带片（%d 波带、外半径 %.2f mm）" % (m_act, r_out * 1e3)),
                det_circle("det", f, max(2 * spot, 12 * dx), label="焦点（成像面）"),
            ],
            sources=[plane_source(wl=wl)],
            category="仪器", difficulty="进阶",
            description="把波带片当作「透镜」用于成像/聚焦：设计波长 %s nm、一级焦距 %.0f mm、"
                        "%d 个波带、最外环宽 %.1f µm（这是加工精度的关键指标）。"
                        % (key, f * 1e3, m_act, dr * 1e6),
            expected="光汇聚到 %.0f mm 的一级焦点，焦斑 ≈λf/(2r_out)=%.1f µm（接近衍射极限 λf/D）；"
                     "相位型（0/π）效率理论 40.5%%，振幅型（明暗交替）≈10%%。"
                     "波带片是 X 射线/极紫外（无合适透镜材料）与软 X 射线显微术的关键元件，"
                     "也是「光子筛」（把环换成环形排列的小孔）的前身。"
                     % (f * 1e3, spot * 1e6),
            tags=["波带片", "显微术", "X射线光学", "衍射效率"],
            theory="r_m²=mλf+(mλ/2)²；效率：相位型 4/π²≈40.5%、振幅型 1/π²≈10%", sample=(i < 2),
        )
        reg.add("inst_zoneplate_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))


# ==========================================================================
# 9. 望远镜（开普勒 / 伽利略 / 扩束器）
# ==========================================================================

def telescopes(reg):
    specs = [("开普勒望远镜 f1=300 f2=50（放大 6×）", 0.3, 0.05, False),
             ("开普勒望远镜 f1=300 f2=25（放大 12×）", 0.3, 0.025, False),
             ("伽利略望远镜 f1=300 f2=−50（放大 6×）", 0.3, -0.05, False),
             ("伽利略望远镜 f1=200 f2=−25（放大 8×）", 0.2, -0.025, False)]
    for i, (lbl, f1, f2, _x) in enumerate(specs):
        wl, N, W = WL["633"], 512, W10
        cfg = build_config(
            name="望远镜：%s" % lbl,
            size=N, width=W, wl=wl,
            components=[
                lens("obj", (0, 0, 0), f=f1, aperture=W * 0.4, label="物镜 f=%.0f mm" % (f1 * 1e3)),
                lens("eye", (0, 0, f1 + f2), f=f2, aperture=W * 0.4, label="目镜 f=%.0f mm" % (f2 * 1e3)),
                det_rect("det", f1 + 2 * f2 + 0.05, W * 0.9, W * 0.9, label="出瞳后"),
            ],
            sources=[plane_source(wl=wl, tilt_x=2e-4, label="轴上点（倾角 0.2 mrad）")],
            category="仪器", difficulty="进阶",
            description="%s：物镜与目镜间距 = f₁+f₂ = %.0f mm（两透镜共焦）。物镜口径 8 mm。" % (lbl, (f1 + f2) * 1e3),
            expected="角放大率 M=−f₁/f₂=%.1f：入射 0.2 mrad 的平行光出射为 %.2f mrad。"
                     "开普勒式在两透镜间有实焦点（可放十字丝/视场光阑），伽利略式（f₂<0）"
                     "为短筒正像系统。出射光束直径 = D·|f₂/f₁| = %.2f mm（这就是「出瞳」）。"
                     % (-f1 / f2, 2e-4 * abs(f1 / f2) * 1e3, 8.0 * abs(f2 / f1)),
            tags=["望远镜", "开普勒", "伽利略", "角放大率"],
            theory="M = −f₁/f₂；两透镜共焦 → 平行光进、平行光出", sample=(i < 2),
        )
        reg.add("inst_telescope_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))


def beam_expander(reg):
    for i, (f1, f2, d0) in enumerate([(0.02, 0.2, 1.0e-3), (0.02, 0.4, 1.0e-3),
                                      (0.05, 0.2, 1.5e-3)]):
        wl, N, W = WL["633"], 512, W10
        M = f2 / f1
        w_out = d0 * M
        cfg = build_config(
            name="激光扩束器：%d× 扩束（f1=%.0f → f2=%.0f mm）" % (M, f1 * 1e3, f2 * 1e3),
            size=N, width=W, wl=wl,
            components=[
                lens("l1", (0, 0, 0), f=f1, aperture=W * 0.4, label="短焦透镜 f=%.0f mm" % (f1 * 1e3)),
                lens("l2", (0, 0, f1 + f2), f=f2, aperture=W * 0.4, label="长焦透镜 f=%.0f mm" % (f2 * 1e3)),
                det_rect("det", f1 + f2 + 0.3, W * 0.9, W * 0.9, label="扩束后 0.3 m 处"),
            ],
            sources=[gauss_source(d0, wl=wl, label="高斯 w₀=%.1f mm" % (d0 * 1e3))],
            category="仪器", difficulty="进阶",
            description="伽利略/开普勒型扩束器：束腰 %.1f mm 的高斯光束经 f₁=%.0f mm、f₂=%.0f mm "
                        "的共焦透镜组，束腰放大 M=f₂/f₁=%.0f×。" % (d0 * 1e3, f1 * 1e3, f2 * 1e3, M),
            expected="出射束腰 %.2f mm（%.0f×放大），发散角相应缩小到 1/%.0f：θ_out=λ/(πw_out)=%.2f mrad。"
                     "扩束器同时压缩发散角——这是激光雷达、干涉测量里「准直」的标准手段。"
                     "两透镜共焦（间距 f₁+f₂）是准直条件。" % (w_out * 1e3, M, M,
                                                            gauss_divergence(wl, w_out) * 1e3),
            tags=["扩束器", "准直", "束腰变换", "望远镜"],
            theory="M = f₂/f₁；w_out = M·w_in，θ_out = θ_in/M", sample=(i == 0),
        )
        reg.add("inst_expander_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))


# ==========================================================================
# 10. 空间滤波器 / 共聚焦针孔
# ==========================================================================

def spatial_filter(reg):
    for i, (d_pin, lbl) in enumerate([(0.5e-3, "针孔 0.5 mm（宽松）"), (0.2e-3, "针孔 0.2 mm"),
                                      (0.05e-3, "针孔 0.05 mm（严格）"), (0.02e-3, "针孔 0.02 mm（近单模）")]):
        wl, N, W, f = WL["633"], 512, W10, 0.1
        N_pix = pick_grid(d_pin, W, minimum=2.0)
        cfg = build_config(
            name="空间滤波器：%s" % lbl,
            size=512, width=W, wl=wl,
            components=[
                aperture("obj", (0, 0, 0), circle(1.5e-3), label="待净化光束（带高频噪声）"),
                lens("l1", (0, 0, 0.0 + f), f=f, aperture=W * 0.4, label="聚焦透镜 f=100 mm"),
                aperture("pin", (0, 0, 2 * f), circle(d_pin / 2), label="针孔 %s" % lbl),
                lens("l2", (0, 0, 3 * f), f=f, aperture=W * 0.4, label="再准直透镜 f=100 mm"),
                det_rect("det", 4 * f + 0.2, W * 0.9, W * 0.9, label="净化后的光束"),
            ],
            sources=[plane_source(wl=wl)],
            category="仪器", difficulty="进阶",
            description="空间滤波器（4f 针孔）：透镜把光束聚焦到针孔，针孔只让近轴的低空间频率通过，"
                        "再由第二个透镜恢复平行光；针孔直径 %s。" % lbl,
            expected="针孔直径 %.2f mm 对应焦面上允许的空间频率上限 ρ_max=%.2f /mm（截止角 θ=%.2f mrad）。"
                     "针孔越小输出越干净（接近理想高斯），但透过率也越低；"
                     "针孔小于聚焦艾里斑时会因衍射而展宽输出。这是激光实验里「净化波前」的标准做法。"
                     % (d_pin * 1e3, d_pin / 2 / (wl * f) * 1e-3, d_pin / 2 / f * 1e3),
            tags=["空间滤波", "针孔", "4f", "波前净化"],
            theory="焦面针孔 = 低通滤波器；截止频率 ρ = r_p/(λf)", sample=(i < 2),
        )
        reg.add("inst_spatial_filter_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))


def confocal(reg):
    for i, (d_pin, dz) in enumerate([(0.05e-3, 0.0), (0.02e-3, 0.0), (0.1e-3, 0.0),
                                     (0.02e-3, 0.5e-3), (0.02e-3, 2.0e-3)]):
        wl, N, W, f = WL["633"], 512, W10, 0.1
        cfg = build_config(
            name="共聚焦针孔：针孔 %.0f µm、离焦 %.1f mm" % (d_pin * 1e3, dz * 1e3),
            size=N, width=W, wl=wl,
            components=[
                lens("l1", (0, 0, 0), f=f, aperture=W * 0.4, label="物镜 f=100 mm"),
                aperture("pin", (0, 0, f + dz), circle(d_pin / 2),
                         label="共聚焦针孔 %.0f µm" % (d_pin * 1e3)),
                det_circle("det", f + dz + 0.01, 3e-3, label="针孔后功率"),
            ],
            sources=[plane_source(wl=wl)],
            category="仪器", difficulty="高级",
            description="共聚焦显微镜的核心：焦面处放针孔（%.0f µm），只收集来自焦平面的光。"
                        "本例把针孔沿轴移动 %.1f mm 模拟**离焦**。" % (d_pin * 1e3, dz * 1e3),
            expected="焦斑（艾里）半径 λf/D=%.1f µm：针孔 %.0f µm（%.1f 个艾里半径）时收集到约 %.0f%% 的焦面光；"
                     "离焦 %.1f mm 后光斑急剧膨胀（z 处光斑 ∝ 离焦量），针孔迅速把光挡掉 → "
                     "这就是共聚焦显微镜的**光学切片**能力（轴向分辨率由针孔尺寸与 NA 共同决定）。"
                     % (1.22 * wl * f / (2 * W * 0.4) * 1e6, d_pin * 1e3,
                        (d_pin / 2) / (1.22 * wl * f / (2 * W * 0.4)), 50.0, dz * 1e3),
            tags=["共聚焦", "针孔", "光学切片", "轴向分辨率"],
            theory="针孔只收集焦面光；离焦光斑放大 → 被针孔拒绝", sample=(i < 2),
        )
        reg.add("inst_confocal_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))


# ==========================================================================
# 11. 相衬显微镜（泽尼克）
# ==========================================================================

def phase_contrast(reg):
    f1 = 0.075
    for i, (phase_q, lbl) in enumerate([(0.5, "π/2（正相衬）"), (1.0, "π（暗场式）"),
                                        (0.25, "π/4"), (3.0, "3π/2（负相衬）")]):
        wl, N, W = WL["633"], 512, W8
        ph = phase_q * PI / 2
        d_plate = 0.005
        n_plate = 1.0 + ph / (k_of(wl) * d_plate)
        comps = [
            comp("spec", "medium", (0, 0, 0), shape=circle(W * 0.3),
                 params={"distance": 0.005, "index": 1.0 + WL["633"] / 5.0, "absorption": 0.0,
                         "steps": 2},
                 label="相位物体（0.2λ 相位延迟）"),
            lens("l1", (0, 0, f1), f=f1, aperture=W * 0.45, label="物镜（傅里叶透镜）f=75 mm"),
            comp("phase_ring", "medium", (0, 0, 2 * f1), shape=circle(0.0005),
                 params={"distance": d_plate, "index": n_plate, "absorption": 0.0, "steps": 2},
                 label="相衬环（在零级上加 %.2f×2π 相位）" % (phase_q / 4)),
            lens("l2", (0, 0, 3 * f1), f=f1, aperture=W * 0.45, label="成像透镜"),
            det_rect("det", 4 * f1, W * 0.9, W * 0.9, label="像面"),
        ]
        cfg = build_config(
            name="相衬显微镜（泽尼克）：相移环 %.2f×2π" % (phase_q / 4),
            size=N, width=W, wl=wl, components=comps,
            sources=[plane_source(wl=wl)],
            category="仪器", difficulty="高级",
            description="泽尼克相衬：样品是纯相位物体（%.2fλ 相位、吸收为零，原本在明场下**看不见**），"
                        "在频谱面的零级位置上加一个相位环（%.2f×2π），使散射光与直射光之间产生"
                        "固定的相位偏置 → 相位差转成强度差。" % (0.2, phase_q / 4),
            expected="无相位环时像面强度均匀（相位物体不可见）；加上 %.2f×2π 的相衬环后，"
                     "样品区与背景区出现明显的明暗对比（相衬成像）。π/2 相移给出「正相衬」"
                     "（光程长的区域更暗），β=3π/2 给出「负相衬」（明暗反转）——"
                     "这正是 1953 年泽尼克诺贝尔奖的核心思想。" % (phase_q / 4),
            tags=["相衬显微镜", "泽尼克", "相位→强度", "生物成像"],
            theory="零级加 π/2 相位 → 散射光与直射光干涉，I ≈ 1 ± 2φ（φ 为样品相位）", sample=(i < 2),
        )
        reg.add("inst_phase_contrast_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))


# ==========================================================================
# 12. 偏振 DIC（微分干涉对比）
# ==========================================================================

def dic_microscope(reg):
    for i, shear in enumerate([0.2e-3, 0.4e-3, 0.8e-3, 1.2e-3]):
        wl, N, W, f1 = WL["633"], 512, W10, 0.075
        comps = [
            _pol("pol1", 0.0, 0.0, "起偏器 45°", r=W * 0.45),
            comp("obj", "medium", (0, 0, 0.02), shape=rectangle(1.0e-3, 1.0e-3),
                 params={"distance": 0.005, "index": 1.0 + WL["633"] / 5.0, "absorption": 0.0,
                         "steps": 2},
                 label="相位物体（方孔，0.2λ）"),
            comp("prism", "uniaxial", (0, 0, 0.04), shape=circle(W * 0.45),
                 params={"distance": shear, "n_o": 1.5, "n_e": 1.5 + WL["633"] * 0.5 / max(shear, 1e-9)},
                 label="DIC 棱镜（剪切 %.2f mm）" % (shear * 1e3)),
            _pol("pol2", 0.06, PI / 2 + PI / 4, "检偏器 135°（与起偏正交）", r=W * 0.45),
            det_rect("det", 0.09, W * 0.9, W * 0.9, label="像面"),
        ]
        cfg = build_config(
            name="偏振 DIC 显微镜：剪切量 %.2f mm" % (shear * 1e3),
            size=N, width=W, wl=wl, components=comps,
            sources=[source("src", "plane", {"power": 1e-3, "polarization": "d"},
                            (0.0, 0.0, SRC_Z), wavelength=wl, label="45° 线偏振")],
            polarized=True,
            category="仪器", difficulty="高级",
            description="微分干涉对比（DIC）：45° 线偏振 → 样品 → **双折射剪切棱镜**（把一束光剪成"
                        "两束偏振正交、横向错开 %.2f mm 的光）→ 正交检偏器。" % (shear * 1e3),
            expected="样品被「剪开」的两点相减后成像，得到**浮雕状**的强度图（强度 ∝ 相位的一阶导数）"
                     "——无染色即可看到透明样品的边缘与厚度梯度，这是 DIC 相对相衬的优势（无光晕、"
                     "可做光学切片）。剪切量 %.2f mm 决定浮雕的感知厚度：过大会出现两个分离的像。"
                     % (shear * 1e3),
            tags=["DIC", "微分干涉", "双折射棱镜", "偏振显微"],
            theory="DIC：正交偏振的剪切光束之差 → 强度 ∝ dφ/dx（相位梯度成像）", sample=(i == 0),
        )
        reg.add("inst_dic_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))


# ==========================================================================
# 13. 超表面透镜与成像系统
# ==========================================================================

def metalens_instruments(reg):
    for i, (wl0, key, f) in enumerate([(633e-9, "633", 0.3), (633e-9, "450", 0.3),
                                       (633e-9, "780", 0.3), (450e-9, "633", 0.3),
                                       (633e-9, "633", 0.5)]):
        wl, N, W = WL[key], 512, 0.006
        R = W * 0.45
        D = 2 * R
        r1 = airy_radius(wl, f, D)
        N2 = pick_grid(r1, W, minimum=4.0)
        zf = f * wl0 / wl
        cfg = build_config(
            name="超表面透镜成像（色散）：设计 λ₀=%d nm，使用 λ=%s nm" % (wl0 * 1e9, key),
            size=N2, width=W, wl=wl,
            components=[
                comp("ml", "metalens_native", (0, 0, 0), shape=circle(R),
                     params={"f": f, "wl0": wl0, "sign": 1},
                     label="超表面透镜（设计 λ₀=%d nm，f=%.2f m）" % (wl0 * 1e9, f)),
                detector("det", zf, label="实际焦点（%.0f mm）" % (zf * 1e3),
                         params={"strehl_aperture": R, "strehl_distance": zf}),
                detector("det0", f, label="设计焦点位置（%.0f mm）" % (f * 1e3)),
            ],
            sources=[plane_source(wl=wl, label="平面波 %s nm" % key)],
            category="仪器", difficulty="高级",
            description="超表面（metasurface）透镜：把设计波长 λ₀=%d nm 下的双曲相位掩模 "
                        "φ(r)=−(2π/λ₀)(√(r²+f²)−f) 直接写进元件；此处用 %s nm 照明，"
                        "焦点按衍射规律移到 z=f·λ₀/λ。" % (wl0 * 1e9, key),
            expected="实际焦点 z=f·λ₀/λ=%.1f mm（设计点 %.0f mm）；偏差 %.1f mm。"
                     "设计点处的光斑明显弥散（离焦 + 球差混合）。超表面透镜最致命的短板就是"
                     "**强色散**（焦点随波长移动 ∝1/λ），这也是「消色差超表面」研究的动机。"
                     % (zf * 1e3, f * 1e3, abs(zf - f) * 1e3),
            tags=["超表面", "metasurface", "色散", "平面光学"],
            theory="φ = −(2π/λ₀)(√(r²+f²)−f)；z_f = f λ₀/λ（衍射元件的强色散）", sample=(i < 2),
        )
        reg.add("inst_metalens_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))

    # 超表面 vs 折射透镜：同焦距对比
    for i, (kind, lbl) in enumerate([("metalens_native", "超表面透镜（双曲相位）"),
                                     ("lens", "理想折射透镜（薄透镜相位）")]):
        wl, N, W, f = WL["633"], 512, 0.006, 0.3
        R = W * 0.45
        r1 = airy_radius(wl, f, 2 * R)
        N2 = pick_grid(r1, W, minimum=4.0)
        if kind == "lens":
            el = lens("e", (0, 0, 0), f=f, aperture=R, label=lbl + " f=0.3 m")
        else:
            el = comp("e", "metalens_native", (0, 0, 0), shape=circle(R),
                      params={"f": f, "wl0": wl, "sign": 1}, label=lbl + " f=0.3 m")
        cfg = build_config(
            name="超表面 vs 折射透镜：%s" % lbl,
            size=N2, width=W, wl=wl,
            components=[el, detector("det", f, label="焦面",
                                     params={"strehl_aperture": R, "strehl_distance": f})],
            sources=[plane_source(wl=wl)],
            category="仪器", difficulty="高级",
            description="在**设计波长**（633 nm）下比较 %s 的聚焦性能。%s" % (lbl,
                        "超表面在 λ₀ 处严格聚焦到 f（相位掩模写死），因此 Strehl 应接近 1。" if kind == "metalens_native"
                        else "理想薄透镜的相位 exp(−ikr²/2f) 是傍轴近似；口径边缘的精确球面波相位略大，"
                             "因此 Strehl 会略低于 1。"),
            expected="两者在设计波长下焦斑都接近艾里极限（r₁=%.1f µm），但走的路不同：超表面用"
                     "写死的分布式相位（可做到亚波长厚度、任意波前），折射透镜用本体厚度积累相位。"
                     "区别在**带宽**：超表面随波长迅速失焦（见 inst_metalens_02/03），"
                     "折射透镜的焦点只在色散修正下缓慢变化。" % (r1 * 1e6),
            tags=["超表面", "对比", "聚焦效率", "色差"],
            theory="超表面相位写死 → 单色完美、宽带失色；折射透镜相位 ∝ k → 色散由材料决定", sample=(i == 0),
        )
        reg.add("inst_metalens_vs_lens_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))


def imaging_systems(reg):
    # 单透镜成像（物在有限距离）—— 用光栅作为「物」
    for i, (object_z, f) in enumerate([(-0.2, 0.1), (-0.15, 0.1), (-0.3, 0.15)]):
        wl, N, W = WL["633"], 512, W10
        a = object_z
        b = 1.0 / (1.0 / f - 1.0 / a)   # 1/v = 1/f - 1/u，u=-a>0
        mag = -b / a
        cfg = build_config(
            name="单透镜成像：物距 %.0f mm、f=%.0f mm → 像距 %.0f mm、放大 %.2f×"
                 % (-a * 1e3, f * 1e3, b * 1e3, abs(mag)),
            size=N, width=W, wl=wl,
            components=[
                comp("g", "grating", (-a, 0, 0), params={"kind": "amplitude_binary",
                                                         "period": 0.5e-3, "duty": 0.5},
                     label="物：周期 0.5 mm 的振幅光栅"),
                lens("l", (0, 0, 0), f=f, aperture=W * 0.45, label="成像透镜 f=%.0f mm" % (f * 1e3)),
                det_rect("det", b, W * 0.9, W * 0.9, label="像面"),
            ],
            sources=[plane_source(wl=wl)],
            category="仪器", difficulty="进阶",
            description="阿贝成像：光栅（周期 0.5 mm）放在物距 %.0f mm，用 f=%.0f mm 透镜成像，"
                        "像距按高斯公式 b=f·a/(a−f)=%.0f mm。" % (-a * 1e3, f * 1e3, b * 1e3),
            expected="像面出现光栅的像，放大率 M=b/a=%.2f（倒立）；由于孔径有限（8 mm 口径）"
                     "只有有限个衍射级参与成像，因此像的对比度取决于「通过了几级」——"
                     "这正是阿贝成像理论：**成像 = 两次衍射**，分辨率由能收集的衍射级数决定。"
                     % abs(mag),
            tags=["阿贝成像", "透镜公式", "分辨率"],
            theory="1/f = 1/a + 1/b；M = b/a；分辨率由收集的衍射级数决定", sample=(i < 2),
        )
        reg.add("inst_imaging_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))

    # 暗场显微镜（中心遮挡）
    for i, (ob, lbl) in enumerate([(0.0, "无遮挡（明场）"), (0.3, "中心遮挡 30%（暗场）"),
                                   (0.6, "中心遮挡 60%（强暗场）")]):
        wl, N, W, f = WL["633"], 512, W10, 0.1
        comps = [lens("l", (0, 0, 0), f=f, aperture=W * 0.45, label="物镜 f=100 mm")]
        if ob == 0.0:
            comps.append(aperture("stop", (0, 0, 0.001), circle(W * 0.45), label="无遮挡"))
        else:
            comps.append(aperture("stop", (0, 0, 0.001),
                                  ring(W * 0.45 * ob, W * 0.45), label="环形光阑（遮挡比 %.0f%%）" % (ob * 100)))
        comps.append(det_rect("det", f, W * 0.9, W * 0.9, label="焦面"))
        cfg = build_config(
            name="暗场/环形光阑成像：%s" % lbl,
            size=N, width=W, wl=wl, components=comps,
            sources=[plane_source(wl=wl)],
            category="仪器", difficulty="进阶",
            description="在物镜入瞳上加环形光阑（%s）：挡掉低频、只让大角度散射光成像。" % lbl,
            expected=("中央亮斑变窄（等效数值孔径提高）、旁瓣增强——这就是「环形光阑超分辨」的代价与收益："
                      "分辨率提高但成像对比度下降（旁瓣是把双刃剑）。暗场显微镜用中心遮挡只收集散射光，"
                      "背景被压暗、透明小颗粒变得可见。" if ob > 0 else
                      "无遮挡的标准艾里图样（r₁=%.1f µm）。" % (1.22 * wl * f / (2 * W * 0.45) * 1e6)),
            tags=["暗场", "环形光阑", "超分辨", "艾里斑"],
            theory="中心遮挡 → 截止频率提高、旁瓣能量增加（暗场/环形照明）",
        )
        reg.add("inst_darkfield_%02d.json" % (i + 1), "仪器", "进阶",
                cfg["meta"]["zod"]["description"], cfg["meta"]["zod"]["expected"], cfg,
                tags=cfg["meta"]["zod"]["tags"], theory=cfg["meta"]["zod"].get("theory"),
                sample=(i < 2))
