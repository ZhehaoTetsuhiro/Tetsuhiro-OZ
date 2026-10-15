"""衍射类示例（diff_*）—— 目标 100+。

覆盖：单缝/圆孔/方孔/椭圆/多边/星形/十字/环形/超椭圆/自定义顶点孔径的夫琅禾费与
菲涅耳衍射；振幅/相位/二元/闪耀光栅与缺级；波带片；泰伯自成像；直边衍射；涡旋远场；
泽尼克像差的衍射受限 PSF；4f 滤波。
"""

import math

from wavegen.common import *   # noqa: F401,F403

# 常用网格
W20 = 0.02      # 20 mm 窗口
W8 = 0.008      # 8 mm
W4 = 0.004      # 4 mm
W10 = 0.01


def _det_end(cid="det", z=1.5, w=W20, h=W20, label="观察屏"):
    return det_rect(cid, z, w * 0.9, h * 0.9, label=label)


# ==========================================================================
# 1. 单缝夫琅禾费（缝宽扫掠 / 波长扫掠）
# ==========================================================================

def single_slit(reg):
    for i, a in enumerate([0.10e-3, 0.15e-3, 0.20e-3, 0.30e-3, 0.40e-3, 0.60e-3, 0.80e-3]):
        wl, z = WL["633"], 1.5
        N, W = 512, W20
        dx = W / N
        check_sampling(dx, {"缝宽 a": a}, 2.5)
        x1 = sinc_first_zero(wl, z, a)
        n_zeros = int((W * 0.45) / x1)
        cfg = build_config(
            name="单缝夫琅禾费衍射：缝宽 a=%.2f mm" % (a * 1e3),
            size=N, width=W, wl=wl,
            components=[
                aperture("slit", (0, 0, 0), rectangle(a, W * 0.95), label="单缝 a=%.2f mm" % (a * 1e3)),
                _det_end(z=z),
            ],
            sources=[plane_source(wl=wl, label="平面波 633 nm")],
            category="衍射", difficulty="入门",
            description="平面波照明宽 %.2f mm 的单缝，传播 %.1f m 后记录远场分布（菲涅耳数 N_F=a²/(λz)=%.3g≪1，属夫琅禾费区）。"
                        % (a * 1e3, z, fresnel_number(a, wl, z)),
            expected="沿 x 为 sinc² 图样：第一暗点 x₁=λz/a=%.2f mm，其后暗点等间距 %.2f mm，"
                     "窗口内可见约 %d 个亮纹；缝越宽中央峰越窄（反比）。缝方向的 y 几乎不变。" % (x1 * 1e3, x1 * 1e3, n_zeros),
            tags=["单缝", "夫琅禾费", "sinc²"],
            theory="I(θ) ∝ sinc²(π a sinθ/λ)，第一暗点 x = λz/a",
            sample=(i < 2),
        )
        reg.add("diff_single_slit_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))

    # 波长扫掠
    for i, key in enumerate(["405", "450", "532", "633", "780"]):
        wl, z, a = WL[key], 1.5, 0.25e-3
        x1 = sinc_first_zero(wl, z, a)
        cfg = build_config(
            name="单缝衍射：波长扫掠 %s nm" % key,
            size=512, width=W20, wl=wl,
            components=[aperture("slit", (0, 0, 0), rectangle(a, W20 * 0.95), label="单缝 a=0.25 mm"),
                        _det_end(z=z)],
            sources=[plane_source(wl=wl, label="平面波 %s nm" % key)],
            category="衍射", difficulty="入门",
            description="同一缝宽 a=0.25 mm、同一距离 z=%.1f m，只换波长（%s nm）。" % (z, key),
            expected="第一暗点 x₁=λz/a=%.2f mm：波长越长衍射角越大、图样按 λ 线性展宽（红光的中央峰比紫光宽 %.2f 倍）。"
                     % (x1 * 1e3, WL[key] / WL["405"]),
            tags=["单缝", "波长", "色散"],
            theory="x₁ = λz/a 与 λ 成正比",
            sample=(i == 0),
        )
        reg.add("diff_single_slit_wl_%02d.json" % (i + 1), "衍射", "入门",
                cfg["meta"]["zod"]["description"], cfg["meta"]["zod"]["expected"], cfg,
                tags=cfg["meta"]["zod"]["tags"], theory=cfg["meta"]["zod"].get("theory"),
                sample=cfg["meta"].get("sample", False))


# ==========================================================================
# 2. 单缝菲涅耳（近场）—— 菲涅耳数扫掠
# ==========================================================================

def single_slit_fresnel(reg):
    for i, z in enumerate([0.02, 0.05, 0.10, 0.20, 0.40]):
        a = 0.40e-3
        wl = WL["633"]
        NF = fresnel_number(a, wl, z)
        dx = W10 / 512
        check_sampling(dx, {"缝宽 a": a}, 2.0)
        cfg = build_config(
            name="单缝菲涅耳（近场）衍射：z=%.0f mm，N_F=%.1f" % (z * 1e3, NF),
            size=512, width=W10, wl=wl,
            components=[aperture("slit", (0, 0, 0), rectangle(a, W10 * 0.95), label="单缝 a=0.40 mm"),
                        _det_end(z=z, w=W10, h=W10, label="近场屏")],
            sources=[plane_source(wl=wl)],
            category="衍射", difficulty="进阶",
            description="缝宽 0.40 mm，观察距离缩到 %.0f mm（菲涅耳数 N_F=a²/(λz)=%.2f）。"
                        "距离不够远，边缘子波的相位差不可忽略。" % (z * 1e3, NF),
            expected="近场不满足远场条件：图样是菲涅耳区结构——亮纹不是等间距，中央强度随 z 起伏"
                     "（N_F≈%s 时轴上%s）。距离越近越接近缝的几何投影（边缘仅 %s）。"
                     % ("%.2f" % NF, "为极大" if NF % 2 > 1 else "为极小",
                        "轻微波纹" if NF > 10 else "明显振荡"),
            tags=["单缝", "菲涅耳衍射", "近场"],
            theory="N_F = a²/(λz)；N_F≪1 才进入夫琅禾费区",
        )
        reg.add("diff_single_slit_fresnel_%02d.json" % (i + 1), "衍射", "进阶",
                cfg["meta"]["zod"]["description"], cfg["meta"]["zod"]["expected"], cfg,
                tags=cfg["meta"]["zod"]["tags"], theory=cfg["meta"]["zod"].get("theory"))


# ==========================================================================
# 3. 圆孔衍射：艾里斑（傅里叶透镜焦面）
# ==========================================================================

def circular_aperture(reg):
    for i, D in enumerate([0.5e-3, 0.8e-3, 1.0e-3, 1.5e-3, 2.0e-3, 3.0e-3]):
        wl, f = WL["633"], 0.3
        W = W8
        r1 = airy_radius(wl, f, D)
        N = 512 if r1 / (W / 512.0) >= 6.0 else 1024
        dx = W / N
        check_sampling(dx, {"口径 D": D}, 8.0)
        check_sampling(dx, {"第一暗环": r1}, 6.0)
        window_ok(W, {"第三暗环": 3.24 * r1})
        cfg = build_config(
            name="圆孔衍射（艾里斑）：口径 D=%.1f mm" % (D * 1e3),
            size=N, width=W, wl=wl,
            components=[
                aperture("ap", (0, 0, 0), circle(D / 2), label="圆孔 D=%.1f mm" % (D * 1e3)),
                lens("fl", (0, 0, 0.02), f=f, aperture=W8 / 2, label="傅里叶透镜 f=0.3 m"),
                det_circle("det", 0.32, W8 / 2 * 0.95, label="焦面"),
            ],
            sources=[plane_source(wl=wl)],
            category="衍射", difficulty="入门",
            description="直径 %.1f mm 圆孔 + 焦距 0.3 m 傅里叶透镜，焦面即孔径的二维傅里叶变换（艾里图样）。" % (D * 1e3),
            expected="中央艾里斑 + 同心暗环：第一暗环半径 r₁=1.22λf/D=%.3f mm，第二暗环 2.23r₁、第三 3.24r₁；"
                     "口径越大光斑越小（r₁∝1/D），中央峰约占总功率 84%%。" % (r1 * 1e3),
            tags=["艾里斑", "圆孔", "衍射极限"],
            theory="r₁ = 1.22 λ f / D",
            sample=(i < 2),
        )
        reg.add("diff_circular_airy_%02d.json" % (i + 1), "衍射", "入门",
                cfg["meta"]["zod"]["description"], cfg["meta"]["zod"]["expected"], cfg,
                tags=cfg["meta"]["zod"]["tags"], theory=cfg["meta"]["zod"].get("theory"),
                sample=cfg["meta"].get("sample", False))

    # 直接传播到远场的艾里斑（不用透镜）
    for i, D in enumerate([1.0e-3, 2.0e-3, 3.0e-3]):
        wl, z = WL["633"], 4.0
        N, W = 1024, W20
        dx = W / N
        check_sampling(dx, {"口径 D": D}, 16.0)
        r1 = airy_radius(wl, z, D)
        check_sampling(dx, {"第一暗环": r1}, 6.0)
        cfg = build_config(
            name="圆孔远场艾里斑（直接传播 z=4 m）：D=%.0f mm" % (D * 1e3),
            size=N, width=W, wl=wl,
            components=[aperture("ap", (0, 0, 0), circle(D / 2), label="圆孔"),
                        _det_end(z=z)],
            sources=[plane_source(wl=wl)],
            category="衍射", difficulty="入门",
            description="不用透镜，直径 %.0f mm 的圆孔直接传播 %.0f m 到屏上（夫琅禾费近似成立）。" % (D * 1e3, z),
            expected="艾里图样，第一暗环 r₁=1.22λz/D=%.2f mm（约 %d 个像素）；与用 f=z 的透镜等价——"
                     "自由空间传播本身就是一次傅里叶变换。" % (r1 * 1e3, r1 / dx),
            tags=["艾里斑", "远场", "直接传播"],
            theory="r₁ = 1.22 λ z / D",
        )
        reg.add("diff_circular_farfield_%02d.json" % (i + 1), "衍射", "入门",
                cfg["meta"]["zod"]["description"], cfg["meta"]["zod"]["expected"], cfg,
                tags=cfg["meta"]["zod"]["tags"], theory=cfg["meta"]["zod"].get("theory"))


# ==========================================================================
# 4. 双缝：衍射包络 × 干涉条纹
# ==========================================================================

def double_slit_diffraction(reg):
    for i, d in enumerate([0.4e-3, 0.6e-3, 1.0e-3, 1.5e-3]):
        a, wl, z, N, W = 0.20e-3, WL["633"], 1.5, 1024, W20
        dx = W / N
        check_sampling(dx, {"缝距 d": d, "缝宽 a": a}, 4.0)
        fring = fringe_spacing(wl, z, d)
        check_sampling(dx, {"条纹间距": fring}, 5.0)
        env = sinc_first_zero(wl, z, a)
        cfg = build_config(
            name="双缝衍射：缝距 d=%.1f mm（包络零点 %.1f mm）" % (d * 1e3, env * 1e3),
            size=N, width=W, wl=wl,
            components=[aperture("ds", (0, 0, 0),
                                 double_slit(d, a, W * 0.95), label="双缝 d=%.1f mm a=0.2 mm" % (d * 1e3)),
                        _det_end(z=z)],
            sources=[plane_source(wl=wl)],
            category="衍射", difficulty="入门",
            description="缝宽固定 0.20 mm、缝距 %0.1f mm 的双缝，%.1f m 处观察。这是「干涉 × 衍射」的复合图样。"
                        % (d * 1e3, z),
            expected="细密干涉条纹（间距 Δx=λz/d=%.3f mm，约 %d 条/mm）被单缝衍射包络调制："
                     "包络第一零点在 λz/a=%.2f mm 处；超过零点条纹重新变亮（π 相位跳变）。"
                     % (fring * 1e3, 1.0 / (fring * 1e3), env * 1e3),
            tags=["双缝", "干涉+衍射", "包络"],
            theory="I ∝ cos²(π d x/(λz)) · sinc²(π a x/(λz))",
            sample=(i < 2),
        )
        reg.add("diff_double_slit_env_%02d.json" % (i + 1), "衍射", "入门",
                cfg["meta"]["zod"]["description"], cfg["meta"]["zod"]["expected"], cfg,
                tags=cfg["meta"]["zod"]["tags"], theory=cfg["meta"]["zod"].get("theory"),
                sample=cfg["meta"].get("sample", False))

    # 缝宽扫掠：包络零点移动
    for i, a in enumerate([0.10e-3, 0.15e-3, 0.30e-3]):
        d, wl, z, N, W = 1.0e-3, WL["633"], 1.5, 1024, W20
        dx = W / N
        check_sampling(dx, {"缝宽 a": a, "缝距 d": d}, 3.0)
        env = sinc_first_zero(wl, z, a)
        cfg = build_config(
            name="双缝衍射：缝宽 a=%.2f mm 的包络调制" % (a * 1e3),
            size=N, width=W, wl=wl,
            components=[aperture("ds", (0, 0, 0),
                                 double_slit(d, a, W * 0.95), label="双缝 a=%.2f mm" % (a * 1e3)),
                        _det_end(z=z)],
            sources=[plane_source(wl=wl)],
            category="衍射", difficulty="进阶",
            description="缝距固定 1.0 mm，缝宽取 %.2f mm——只有衍射包络随缝宽改变，条纹间距不变。" % (a * 1e3),
            expected="条纹间距仍为 λz/d=%.3f mm；包络零点移到 λz/a=%.2f mm，包络内干涉条纹数约为 2a/d 的倒数×2："
                     "缝越窄包络越宽、可见条纹越多。" % (fringe_spacing(wl, z, d) * 1e3, env * 1e3),
            tags=["双缝", "缝宽", "包络宽度"],
            theory="包络第一零点由缝宽 a 决定，条纹间距由缝距 d 决定",
        )
        reg.add("diff_double_slit_width_%02d.json" % (i + 1), "衍射", "进阶",
                cfg["meta"]["zod"]["description"], cfg["meta"]["zod"]["expected"], cfg,
                tags=cfg["meta"]["zod"]["tags"], theory=cfg["meta"]["zod"].get("theory"))


# ==========================================================================
# 5. 多缝（N 缝阵列）：主极大/次极大
# ==========================================================================

def multi_slit(reg):
    for i, n in enumerate([3, 4, 5, 6, 8, 12]):
        d, a, wl, z, N, W = 0.5e-3, 0.10e-3, WL["633"], 1.5, 1024, W20
        dx = W / N
        check_sampling(dx, {"缝距 d": d, "缝宽 a": a}, 2.0)
        fring = fringe_spacing(wl, z, d)
        # 多缝用「多个矩形孔径并排」表达：aperture 元件在内核只支持双缝，
        # 这里改用 N 个 aperture 元件的 parallel 语义（并集）实现 N 缝阵列。
        comps = [aperture("slit0", (0, 0, 0), rectangle(a, W * 0.95),
                          label="第 1 缝", parallel=True, x=-0.5 * (n - 1) * d)]
        for m in range(1, n):
            comps.append(aperture("slit%d" % m, (0, 0, 0), rectangle(a, W * 0.95),
                                  label="第 %d 缝" % (m + 1), parallel=True,
                                  x=(-0.5 * (n - 1) + m) * d))
        comps.append(_det_end(z=z))
        cfg = build_config(
            name="%d 缝衍射：主极大与次极大" % n,
            size=N, width=W, wl=wl, components=comps,
            sources=[plane_source(wl=wl)],
            category="衍射", difficulty="进阶",
            description="缝宽 0.10 mm、缝距 0.50 mm 的 %d 缝阵列（用同 z 上 parallel 并集的孔径描述）。" % n,
            expected="主极大位置 x_m=m·λz/d（间距 %.3f mm）固定不变；两条主极大之间出现 %d 个次极大，"
                     "次极大强度 < 主极大的 1/8 且随 N 增大迅速变暗 → 主极大变锐（谱线宽度 ∝ 1/N）。"
                     % (fring * 1e3, n - 2),
            tags=["多缝", "主极大", "光栅雏形"],
            theory="主极大 x_m = m λ z / d；N 缝时主极大宽度 ∝ 1/N",
            sample=(i < 2),
        )
        reg.add("diff_multi_slit_%02d.json" % (i + 1), "衍射", "进阶",
                cfg["meta"]["zod"]["description"], cfg["meta"]["zod"]["expected"], cfg,
                tags=cfg["meta"]["zod"]["tags"], theory=cfg["meta"]["zod"].get("theory"),
                sample=cfg["meta"].get("sample", False))


# ==========================================================================
# 6. 各种孔径形状（正方形/椭圆/三角/多边形/星形/十字/环形/超椭圆/自定义）
# ==========================================================================

def shapes(reg):
    # 方孔 / 矩孔
    for i, (lbl, sh, dim) in enumerate([
        ("方孔 1.0 mm", rectangle(1.0e-3, 1.0e-3), 1.0e-3),
        ("矩孔 2.0×0.5 mm", rectangle(2.0e-3, 0.5e-3), 0.5e-3),
        ("方孔 2.0 mm", square(2.0e-3), 2.0e-3),
    ]):
        wl, z, N, W = WL["633"], 2.0, 1024, W20
        dx = W / N
        check_sampling(dx, {"最短边": dim}, 8.0)
        x1 = sinc_first_zero(wl, z, dim)
        cfg = build_config(
            name="方/矩孔夫琅禾费衍射：%s" % lbl,
            size=N, width=W, wl=wl,
            components=[aperture("ap", (0, 0, 0), sh, label=lbl), _det_end(z=z)],
            sources=[plane_source(wl=wl)],
            category="衍射", difficulty="入门",
            description="%s 的矩形孔径，传播 %.0f m 观察远场。" % (lbl, z),
            expected="二维 sinc² 图样，沿两方向的零点由各自边长决定：最短边方向第一零点 λz/a=%.2f mm，"
                     "垂直方向 %.2f mm；图样为十字形亮带 + 副极大网格。" % (x1 * 1e3, x1 * 1e3),
            tags=["方孔", "矩形孔", "二维sinc"],
            theory="I ∝ sinc²(π a x/(λz))·sinc²(π b y/(λz))",
        )
        reg.add("diff_shape_rect_%02d.json" % (i + 1), "衍射", "入门",
                cfg["meta"]["zod"]["description"], cfg["meta"]["zod"]["expected"], cfg,
                tags=cfg["meta"]["zod"]["tags"], theory=cfg["meta"]["zod"].get("theory"),
                sample=(i == 0))

    # 椭圆
    for i, (a, b) in enumerate([(0.6e-3, 0.3e-3), (1.0e-3, 0.4e-3), (0.4e-3, 0.4e-3)]):
        wl, z, N, W = WL["633"], 2.0, 1024, W20
        dx = W / N
        check_sampling(dx, {"半轴 b": b}, 8.0)
        cfg = build_config(
            name="椭圆孔衍射：a=%.1f mm, b=%.1f mm" % (a * 1e3, b * 1e3),
            size=N, width=W, wl=wl,
            components=[aperture("ap", (0, 0, 0), ellipse(a, b),
                                 label="椭圆孔 a=%.1f b=%.1f mm" % (a * 1e3, b * 1e3)),
                        _det_end(z=z)],
            sources=[plane_source(wl=wl)],
            category="衍射", difficulty="进阶",
            description="半轴 %.1f×%.1f mm 的椭圆孔，远场 %.0f m。" % (a * 1e3, b * 1e3, z),
            expected="中央亮斑沿两方向的角宽度反比于对应半轴：x 方向零点 λz/(2a)=%.2f mm、y 方向 λz/(2b)=%.2f mm；"
                     "圆孔极限（a=b）退化为艾里图样。" % (sinc_first_zero(wl, z, 2 * a) * 1e3,
                                                          sinc_first_zero(wl, z, 2 * b) * 1e3),
            tags=["椭圆孔", "各向异性"],
            theory="椭圆孔的远场是贝塞尔-傅里叶变换，尺度 ∝ 1/半轴",
        )
        reg.add("diff_shape_ellipse_%02d.json" % (i + 1), "衍射", "进阶",
                cfg["meta"]["zod"]["description"], cfg["meta"]["zod"]["expected"], cfg,
                tags=cfg["meta"]["zod"]["tags"], theory=cfg["meta"]["zod"].get("theory"))

    # 三角孔
    for i, r in enumerate([0.8e-3, 1.5e-3]):
        wl, z, N, W = WL["633"], 2.0, 1024, W20
        dx = W / N
        check_sampling(dx, {"外接圆半径": r}, 16.0)
        cfg = build_config(
            name="三角孔衍射：外接圆半径 %.1f mm" % (r * 1e3),
            size=N, width=W, wl=wl,
            components=[aperture("ap", (0, 0, 0), triangle(r),
                                 label="三角孔 R=%.1f mm" % (r * 1e3), rotation=0.0),
                        _det_end(z=z)],
            sources=[plane_source(wl=wl)],
            category="衍射", difficulty="进阶",
            description="等边三角孔径（外接圆半径 %.1f mm），远场 %.0f m。" % (r * 1e3, z),
            expected="六角形对称的衍射星芒（三角孔的三重对称在强度图样上表现为六重对称），"
                     "第一零点尺度 λz/(2R)=%.2f mm；旋转孔径会让图样整体旋转。" % (sinc_first_zero(wl, z, 2 * r) * 1e3),
            tags=["三角孔", "对称性"],
            theory="孔径的夫琅禾费图样带孔径的旋转对称性（三次对称 → 图样六重对称）",
        )
        reg.add("diff_shape_triangle_%02d.json" % (i + 1), "衍射", "进阶",
                cfg["meta"]["zod"]["description"], cfg["meta"]["zod"]["expected"], cfg,
                tags=cfg["meta"]["zod"]["tags"], theory=cfg["meta"]["zod"].get("theory"))

    # 多边形
    for i, (r, sides) in enumerate([(1.0e-3, 5), (1.0e-3, 6), (1.0e-3, 8), (0.8e-3, 12)]):
        wl, z, N, W = WL["633"], 2.0, 1024, W20
        dx = W / N
        check_sampling(dx, {"外接圆半径": r}, 16.0)
        cfg = build_config(
            name="正 %d 边形孔衍射：R=%.1f mm" % (sides, r * 1e3),
            size=N, width=W, wl=wl,
            components=[aperture("ap", (0, 0, 0), polygon(r, sides),
                                 label="正 %d 边形孔" % sides),
                        _det_end(z=z)],
            sources=[plane_source(wl=wl)],
            category="衍射", difficulty="进阶",
            description="外接圆半径 %.1f mm 的正 %d 边形孔，远场 %.0f m。" % (r * 1e3, sides, z),
            expected="%d 重对称的星芒图样：边数越多越接近圆孔的艾里斑（%d 边形已有明显同心环结构，"
                     "第一暗环尺度 ≈1.22λz/(2R)=%.2f mm）。" % (sides, sides, airy_radius(wl, z, 2 * r) * 1e3),
            tags=["多边形孔", "对称性", "趋近圆孔"],
            theory="正多边形孔 → 边数 →∞ 时退化为圆孔的 2J₁(x)/x",
        )
        reg.add("diff_shape_polygon_%02d.json" % (i + 1), "衍射", "进阶",
                cfg["meta"]["zod"]["description"], cfg["meta"]["zod"]["expected"], cfg,
                tags=cfg["meta"]["zod"]["tags"], theory=cfg["meta"]["zod"].get("theory"))

    # 星形
    for i, (r, inner, points) in enumerate([(1.0e-3, 0.4e-3, 5), (1.0e-3, 0.3e-3, 6), (0.8e-3, 0.5e-3, 8)]):
        wl, z, N, W = WL["633"], 2.0, 1024, W20
        dx = W / N
        check_sampling(dx, {"外接圆半径": r, "内半径": inner}, 12.0)
        cfg = build_config(
            name="星形孔衍射：%d 角星（R=%.1f mm）" % (points, r * 1e3),
            size=N, width=W, wl=wl,
            components=[aperture("ap", (0, 0, 0), star(r, inner, points),
                                 label="%d 角星孔" % points),
                        _det_end(z=z)],
            sources=[plane_source(wl=wl)],
            category="衍射", difficulty="高级",
            description="%d 角星形孔径（外接圆半径 %.1f mm、谷底 %.1f mm），远场 %.0f m。" % (points, r * 1e3, inner * 1e3, z),
            expected="%d 重对称的尖刺状星芒：尖角越锐（inner/radius 越小）星芒越长越细；"
                     "星芒方向对应孔径尖角方向的进一步细化 → 用孔径整形可以得到任意「衍射星冠」。" % points,
            tags=["星形孔", "衍射星芒", "孔径整形"],
            theory="尖角处的高频空间频率分量 → 远场沿尖角方向的长条亮斑",
        )
        reg.add("diff_shape_star_%02d.json" % (i + 1), "衍射", "高级",
                cfg["meta"]["zod"]["description"], cfg["meta"]["zod"]["expected"], cfg,
                tags=cfg["meta"]["zod"]["tags"], theory=cfg["meta"]["zod"].get("theory"))

    # 十字
    for i, (ln, wd) in enumerate([(3.0e-3, 0.3e-3), (2.0e-3, 0.5e-3)]):
        wl, z, N, W = WL["633"], 2.0, 1024, W20
        dx = W / N
        check_sampling(dx, {"臂宽": wd}, 8.0)
        cfg = build_config(
            name="十字孔衍射：臂长 %.1f mm 臂宽 %.2f mm" % (ln * 1e3, wd * 1e3),
            size=N, width=W, wl=wl,
            components=[aperture("ap", (0, 0, 0), cross(ln, wd), label="十字孔"),
                        _det_end(z=z)],
            sources=[plane_source(wl=wl)],
            category="衍射", difficulty="进阶",
            description="十字形孔径（两臂端到端 %.1f mm、臂宽 %.2f mm），远场 %.0f m。" % (ln * 1e3, wd * 1e3, z),
            expected="「十字叉」衍射图：亮纹沿两臂方向展开（纵向狭缝 → 横向展开、横向狭缝 → 纵向展开），"
                     "两臂各自的 sinc 乘积 = 十字形亮带；臂宽决定亮带宽度 λz/w=%.2f mm。"
                     % (sinc_first_zero(wl, z, wd) * 1e3),
            tags=["十字孔", "可分离孔径"],
            theory="叉形 = 两个正交矩形的并集，sinc 结构沿两轴分开",
        )
        reg.add("diff_shape_cross_%02d.json" % (i + 1), "衍射", "进阶",
                cfg["meta"]["zod"]["description"], cfg["meta"]["zod"]["expected"], cfg,
                tags=cfg["meta"]["zod"]["tags"], theory=cfg["meta"]["zod"].get("theory"))

    # 环形孔
    for i, (rin, rout) in enumerate([(0.5e-3, 1.0e-3), (0.8e-3, 1.0e-3), (0.2e-3, 1.0e-3), (0.95e-3, 1.0e-3)]):
        wl, f, N, W = WL["633"], 0.3, 1024, W8
        dx = W / N
        check_sampling(dx, {"外半径": rout, "内半径": rin}, 8.0)
        r1 = airy_radius(wl, f, 2 * rout)
        cfg = build_config(
            name="环孔（环形光阑）衍射：内径/外径 = %.1f/%.1f" % (rin / rout, 1.0),
            size=N, width=W, wl=wl,
            components=[aperture("ap", (0, 0, 0), ring(rin, rout),
                                 label="环形孔 rin=%.2f mm" % (rin * 1e3)),
                        lens("fl", (0, 0, 0.02), f=f, aperture=W8 / 2, label="傅里叶透镜"),
                        det_circle("det", 0.32, W / 2 * 0.95, label="焦面")],
            sources=[plane_source(wl=wl)],
            category="衍射", difficulty="进阶",
            description="外径 1.0 mm 的环形孔径，内径 %.2f mm（遮挡比 %.1f），经焦距 0.3 m 透镜到焦面。"
                        % (rin * 1e3, rin / rout),
            expected="遮挡比越大中央亮斑越窄（环形孔径等效口径变大），但旁瓣能量升高："
                     "遮挡比 0.5 时中央峰宽约缩小 %.0f%%，第一暗环位置 r≈1.22λf/(2·rout)=%.3f mm。" 
                     % (100 * (1 - rin / rout), r1 * 1e3),
            tags=["环形孔", "环形光阑", "超分辨"],
            theory="环形瞳 → 中央峰变窄、斯特列尔比下降（旁瓣增强）",
            sample=(i == 0),
        )
        reg.add("diff_shape_ring_%02d.json" % (i + 1), "衍射", "进阶",
                cfg["meta"]["zod"]["description"], cfg["meta"]["zod"]["expected"], cfg,
                tags=cfg["meta"]["zod"]["tags"], theory=cfg["meta"]["zod"].get("theory"),
                sample=cfg["meta"].get("sample", False))

    # 超椭圆
    for i, (a, b, order) in enumerate([(1.0e-3, 1.0e-3, 2.5), (1.0e-3, 1.0e-3, 5.0),
                                       (1.0e-3, 0.6e-3, 8.0)]):
        wl, z, N, W = WL["633"], 2.0, 1024, W20
        dx = W / N
        check_sampling(dx, {"最短半轴": min(a, b)}, 8.0)
        cfg = build_config(
            name="超椭圆孔衍射：阶数 m=%.1f" % order,
            size=N, width=W, wl=wl,
            components=[aperture("ap", (0, 0, 0), superellipse(a, b, order),
                                 label="超椭圆 a=%.1f b=%.1f m=%.1f" % (a * 1e3, b * 1e3, order)),
                        _det_end(z=z)],
            sources=[plane_source(wl=wl)],
            category="衍射", difficulty="高级",
            description="超椭圆孔径 (|u|/a)^m+(|v|/b)^m≤1，m=%.1f：m→2 是椭圆、m→∞ 趋近矩形。" % order,
            expected="介于椭圆与矩形之间的图样：m=2.5 仍接近艾里型同心环，m≥5 已出现矩形的十字亮带与"
                     "方形副极大网格（角越尖 → 副极大越强）。",
            tags=["超椭圆孔", "形状内插"],
            theory="超椭圆是椭圆→矩形之间的连续内插，图样随之从 J₁ 型过渡到 sinc 型",
        )
        reg.add("diff_shape_superellipse_%02d.json" % (i + 1), "衍射", "高级",
                cfg["meta"]["zod"]["description"], cfg["meta"]["zod"]["expected"], cfg,
                tags=cfg["meta"]["zod"]["tags"], theory=cfg["meta"]["zod"].get("theory"))

    # 自定义顶点
    verts = [
        ("0.8e-3,0;0.4e-3,0.7e-3;-0.4e-3,0.7e-3;-0.8e-3,0;-0.4e-3,-0.7e-3;0.4e-3,-0.7e-3",
         "六边形", 0.8e-3),
        ("0,1e-3;0.3e-3,0.3e-3;1e-3,0.3e-3;0.5e-3,-0.2e-3;0.7e-3,-1e-3;0,-0.5e-3;-0.7e-3,-1e-3;-0.5e-3,-0.2e-3;-1e-3,0.3e-3;-0.3e-3,0.3e-3",
         "十角星", 1.0e-3),
        ("1e-3,0.4e-3;0.6e-3,0.4e-3;0.6e-3,1e-3;-0.6e-3,1e-3;-0.6e-3,0.4e-3;-1e-3,0.4e-3;-1e-3,-0.4e-3;-0.6e-3,-0.4e-3;-0.6e-3,-1e-3;0.6e-3,-1e-3;0.6e-3,-0.4e-3;1e-3,-0.4e-3",
         "加号（粗十字）", 0.8e-3),
    ]
    for i, (v, lbl, scale) in enumerate(verts):
        wl, z, N, W = WL["633"], 2.0, 1024, W20
        dx = W / N
        check_sampling(dx, {"顶点尺度": scale}, 12.0)
        cfg = build_config(
            name="自定义顶点孔径衍射：%s" % lbl,
            size=N, width=W, wl=wl,
            components=[aperture("ap", (0, 0, 0), custom(v), label="自定义孔径：%s" % lbl),
                        _det_end(z=z)],
            sources=[plane_source(wl=wl)],
            category="衍射", difficulty="高级",
            description="用顶点列表手工定义的%s孔径（自动闭合），远场 %.0f m——演示「任意形状孔径」的衍射计算。"
                        % (lbl, z),
            expected="衍射图样是第一类贝塞尔型线积分的结果，其对称性由顶点表决定；%s 产生对应旋向的"
                     "星芒。改变顶点即改变图样，可用于逆向设计孔径。" % lbl,
            tags=["自定义孔径", "顶点列表", "自由曲面"],
            theory="孔径函数 A(x,y) 的二维傅里叶变换，任意多边形顶点可用格林公式解析给出",
        )
        reg.add("diff_shape_custom_%02d.json" % (i + 1), "衍射", "高级",
                cfg["meta"]["zod"]["description"], cfg["meta"]["zod"]["expected"], cfg,
                tags=cfg["meta"]["zod"]["tags"], theory=cfg["meta"]["zod"].get("theory"))
