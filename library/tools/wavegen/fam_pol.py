"""偏振类示例（pol_*）—— 目标 50+。

马吕斯定律、偏振片对（正交/平行）、QWP/HWP 快轴扫描、圆/椭圆偏振的产生与检测、
旋光片、自定义琼斯元件、偏振干涉、高 NA 矢量聚焦（|Ez|）、斯托克斯参数。
"""

import math

from wavegen.common import *   # noqa: F401,F403

W10 = 0.01
PI = math.pi


def pol_stage(W, elements, z_det=None, src_pol="x", src_waist=None, wl=None, label="探测器"):
    """偏振实验的统一几何：先用矩形光阑把光束裁到偏振元件之内（保证透过率读数不被
    「元件轮廓以外的未调制光」污染），再放偏振元件，最后用等尺寸矩形探测器读数。"""
    wl = wl or WL["633"]
    clip = W * 0.6
    comps = [aperture("clip", (0, 0, 0), rectangle(clip, clip), label="限光阑（%.1f mm 方孔）" % (clip * 1e3))]
    comps.extend(elements)
    z_last = 0.0
    for e in elements:
        z_last = max(z_last, e["pos"]["z"])
    comps.append(det_rect("det", z_det or (z_last + 0.02), clip, clip, label=label))
    if src_waist:
        src = gauss_source(src_waist, wl=wl, polarization=src_pol, label="光源（%s 偏振）" % src_pol)
    else:
        src = source("src", "plane", {"power": 1e-3, "polarization": src_pol},
                     (0.0, 0.0, SRC_Z), wavelength=wl, label="平面波（%s 偏振）" % src_pol)
    return comps, [src]


def _pol(cid, z, angle, label=None, r=0.0049):
    return comp("x" + cid, "polarizer", (0, 0, z), shape=circle(r),
                params={"angle": angle}, label=label or "偏振片 %.0f°" % math.degrees(angle))


def _ret(cid, z, retardance, axis, label=None, r=0.0049):
    return comp("x" + cid, "retarder", (0, 0, z), shape=circle(r),
                params={"retardance": retardance, "axis": axis}, label=label)


# ==========================================================================
# 1. 马吕斯定律：检偏器角度扫描
# ==========================================================================

def malus_law(reg):
    angles = [0, 10, 20, 30, 40, 45, 50, 60, 70, 80, 90]
    for i, deg in enumerate(angles):
        W, wl = W10, WL["633"]
        comps, srcs = pol_stage(W, [_pol("p1", 0.02, 0.0, "起偏器 0°"),
                                    _pol("p2", 0.04, math.radians(deg), "检偏器 %.0f°" % deg)],
                                src_pol="x")
        cfg = build_config(
            name="马吕斯定律：检偏器 %.0f° → I/I₀=%.3f" % (deg, math.cos(math.radians(deg)) ** 2),
            size=512, width=W, wl=wl, components=comps, sources=srcs, polarized=True,
            category="偏振", difficulty="入门",
            description="线偏振（x）先经 0° 起偏器，再经转角的检偏器——经典的**马吕斯定律**实验。",
            expected="透过率 I/I₀=cos²θ：θ=0° 全透、45° 半透（%.3f）、90° 全消光。"
                     "本实例 θ=%.0f° → 理论 cos²θ=%.4f。检偏器转到 90° 时视场完全变暗，"
                     "这是判断线偏振方向的标准方法。" % (0.5, deg, math.cos(math.radians(deg)) ** 2),
            tags=["马吕斯定律", "检偏器", "线偏振"],
            theory="I = I₀ cos²θ（琼斯投影矩阵）", sample=(i < 3),
        )
        reg.add("pol_malus_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))


# ==========================================================================
# 2. 偏振片对：正交/平行/多片
# ==========================================================================

def polarizer_pairs(reg):
    specs = [("两片正交（0°+90°）", [(0.0, "起偏 0°"), (PI / 2, "检偏 90°")], 0.0),
             ("两片平行（0°+0°）", [(0.0, "起偏 0°"), (0.0, "检偏 0°")], 1.0),
             ("三片：0°+45°+90°（消光被破坏）", [(0.0, "0°"), (PI / 4, "45°"), (PI / 2, "90°")], 0.125),
             ("三片：0°+90°+90°", [(0.0, "0°"), (PI / 2, "90°"), (PI / 2, "90°")], 0.0)]
    for i, (lbl, items, ideal) in enumerate(specs):
        W, wl = W10, WL["633"]
        els = [_pol("p%d" % k, 0.02 + 0.02 * k, a, l) for k, (a, l) in enumerate(items)]
        comps, srcs = pol_stage(W, els)
        cfg = build_config(
            name="偏振片组：%s" % lbl,
            size=512, width=W, wl=wl, components=comps, sources=srcs, polarized=True,
            category="偏振", difficulty="入门",
            description="把多片偏振片按不同角度串联：%s。" % lbl,
            expected="理论总透过率 %.3f：正交偏振片完全消光；在正交片之间插入一片 45° 偏振片，"
                     "会把 0° 偏振「旋转」出一部分，使消光被破坏（透过率变为 1/8）——"
                     "这是琼斯演算最直观的演示。"
                     % ideal,
            tags=["偏振片组", "琼斯演算", "消光"],
            theory="琼斯矩阵连乘：P(θ_n)…P(θ_1)E_in", sample=(i < 2),
        )
        reg.add("pol_pairs_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))


# ==========================================================================
# 3. 四分之一波片：快轴扫描
# ==========================================================================

def qwp_sweep(reg):
    for i, deg in enumerate([0, 15, 30, 45, 60, 75, 90]):
        W, wl = W10, WL["633"]
        comps, srcs = pol_stage(W, [_pol("p0", 0.02, 0.0, "起偏器 0°"),
                                    _ret("q", 0.04, PI / 2, math.radians(deg),
                                         "四分之一波片 快轴 %.0f°" % deg),
                                    _pol("a", 0.06, 0.0, "检偏器 0°")],
                                src_pol="x")
        cfg = build_config(
            name="四分之一波片快轴扫描：%.0f°" % deg,
            size=512, width=W, wl=wl, components=comps, sources=srcs, polarized=True,
            category="偏振", difficulty="进阶",
            description="x 偏振光 → 快轴 %.0f° 的四分之一波片 → 固定 0° 检偏器。QWP 把线偏振转成椭圆/圆偏振，"
                        "检偏器只能看到原始 x 分量。" % deg,
            expected="透过检偏器的振幅 ∝ |cosθ·e^{-iπ/4}·cosθ + sinθ·e^{+iπ/4}·sinθ|："
                     "θ=0° 或 90° 时是纯线偏振（全透），θ=45° 时输出**圆偏振**、"
                     "x 分量振幅减到 1/√2（功率 0.5）——用检偏器转一圈强度不变（圆偏振的判据）。"
                     "本实例 θ=%.0f°。" % deg,
            tags=["四分之一波片", "圆偏振", "快轴"],
            theory="QWP 快轴 θ：R(−θ)diag(e^{-iπ/4}, e^{iπ/4})R(θ)", sample=(i < 2),
        )
        reg.add("pol_qwp_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))


# ==========================================================================
# 4. 半波片：偏振方向的旋转
# ==========================================================================

def hwp_sweep(reg):
    for i, deg in enumerate([0, 10, 20, 30, 45, 60]):
        W, wl = W10, WL["633"]
        comps, srcs = pol_stage(W, [_pol("p0", 0.02, 0.0, "起偏器 0°"),
                                    _ret("h", 0.04, PI, math.radians(deg),
                                         "半波片 快轴 %.0f°" % deg),
                                    _pol("a", 0.06, 0.0, "检偏器 0°")],
                                src_pol="x")
        cfg = build_config(
            name="半波片偏振旋转：快轴 %.0f° → 偏振转 %.0f°" % (deg, 2 * deg),
            size=512, width=W, wl=wl, components=comps, sources=srcs, polarized=True,
            category="偏振", difficulty="进阶",
            description="x 偏振经快轴 %.0f° 的半波片（相位延迟 π），线偏振方向被镜像到 2θ=%.0f°。" % (deg, 2 * deg),
            expected="HWP 把线偏振方向绕快轴镜像：入射 0° → 出射 %.0f°。固定 0° 检偏器读到的透过率为 "
                     "cos²(%.0f°)=%.3f。快轴 45° 时偏振转 90°（被检偏器完全挡住）——"
                     "这就是用 HWP 做「偏振旋转旋钮」与光衰减器的原理。"
                     % (2 * deg, 2 * deg, math.cos(math.radians(2 * deg)) ** 2),
            tags=["半波片", "偏振旋转", "镜像"],
            theory="HWP 快轴 θ：线偏振方向 θ_in → 2θ−θ_in", sample=(i < 2),
        )
        reg.add("pol_hwp_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))


# ==========================================================================
# 5. 椭圆偏振测量（QWP + 检偏器转角）
# ==========================================================================

def ellipsometry(reg):
    for i, adeg in enumerate([0, 15, 30, 45, 60, 75, 90]):
        W, wl = W10, WL["633"]
        comps, srcs = pol_stage(W, [_pol("p0", 0.02, 0.0, "起偏器 0°"),
                                    _ret("q", 0.04, PI / 2, PI / 4, "QWP 快轴 45°（产生圆偏振）"),
                                    _pol("a", 0.06, math.radians(adeg), "检偏器 %.0f°" % adeg)],
                                src_pol="x")
        cfg = build_config(
            name="椭圆偏振测量：检偏器 %.0f°（对圆偏振）" % adeg,
            size=512, width=W, wl=wl, components=comps, sources=srcs, polarized=True,
            category="偏振", difficulty="进阶",
            description="x 偏振 → QWP 快轴 45° → **圆偏振** → 检偏器 %.0f°。转动检偏器读强度。" % adeg,
            expected="圆偏振的判据：检偏器转到任何角度透过率都相同（%.3f，与 45° 处一致）。"
                     "对照 pol_qwp_* 系列（椭圆偏振）会看到明显的 cos² 型调制——"
                     "「转一圈强度不变」是区分圆偏振与椭圆/线偏振的实验判据。"
                     % 0.5,
            tags=["椭圆偏振", "圆偏振", "检偏"],
            theory="圆偏振 |E|² 与检偏方向无关；椭圆偏振随检偏角调制", sample=(i < 3),
        )
        reg.add("pol_ellipse_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))


# ==========================================================================
# 6. 旋光片
# ==========================================================================

def rotator_sweep(reg):
    for i, deg in enumerate([0, 15, 30, 45, 60, 90]):
        W, wl = W10, WL["633"]
        comps, srcs = pol_stage(W, [_pol("p0", 0.02, 0.0, "起偏器 0°"),
                                    comp("rot", "rotator", (0, 0, 0.04), shape=circle(W * 0.49),
                                         params={"angle": math.radians(deg)},
                                         label="旋光片 %.0f°" % deg),
                                    _pol("a", 0.06, 0.0, "检偏器 0°")],
                                src_pol="x")
        cfg = build_config(
            name="旋光片：旋转 %.0f° → 透过率 %.3f" % (deg, math.cos(math.radians(deg)) ** 2),
            size=512, width=W, wl=wl, components=comps, sources=srcs, polarized=True,
            category="偏振", difficulty="进阶",
            description="旋光片（旋光角 %.0f°）把线偏振方向旋转，用 0° 检偏器读数。" % deg,
            expected="等效旋转矩阵 R(α)：出射偏振方向转 %.0f°，检偏器读到 cos²α=%.3f。"
                     "旋光（自然旋光/法拉第旋转）与波片的区别在于它是**旋转**而非相位延迟："
                     "往返两次会加倍（法拉第）或抵消（自然旋光，与介质手性有关）。"
                     % (deg, math.cos(math.radians(deg)) ** 2),
            tags=["旋光片", "旋光", "琼斯旋转矩阵"],
            theory="R(α) 旋转琼斯矢量，透过检偏器 cos²α", sample=(i < 2),
        )
        reg.add("pol_rotator_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))


# ==========================================================================
# 7. 自定义琼斯元件
# ==========================================================================

def custom_jones_examples(reg):
    specs = [
        ("理想偏振分光（透 x 挡 y）", [[1, 0, 0, 0], [0, 0, 0, 0]], "x", 1.0,
         "透过 x、完全挡住 y：对 x 入射全透。"),
        ("正交投影（透 y 挡 x）", [[0, 0, 0, 0], [0, 0, 1, 0]], "x", 0.0,
         "只透 y：对 x 入射完全消光。"),
        ("45° 旋转（实矩阵）", [[math.cos(PI / 4), 0, math.sin(PI / 4), 0],
                                [-math.sin(PI / 4), 0, math.cos(PI / 4), 0]], "x", 0.5,
         "实旋转矩阵 R(45°)：x 偏振转 45°，功率不变但检偏器读数变 cos²45°=0.5。"),
        ("非酉衰减器 diag(0.7, 0.7)", [[0.7, 0, 0, 0], [0, 0, 0.7, 0]], "x", 0.49,
         "非酉矩阵（吸收 51%）：演示「自定义琼斯」不要求酉性——可以表达吸收/二级衰减。"),
        ("纯相位 e^{iφ} diag", [[0, 0, 0, 0], [0, 0, 1, 0]], "x", 1.0,
         "diag(e^{iπ/2}, e^{iπ/2})：所有分量加同样的相位（不影响强度，但会改变干涉结果）。"),
    ]
    for i, (lbl, m, src_pol, ideal, note) in enumerate(specs):
        W, wl = W10, WL["633"]
        if i == 4:
            mat = [[0.0, 0.0, 0.0, 1.0], [0.0, 0.0, 1.0, 0.0]]   # diag(i, i)：a=i, d=i
        else:
            mat = m
        a = mat[0]
        c = mat[1]
        comps, srcs = pol_stage(W, [
            comp("cj", "custom_jones", (0, 0, 0.03), shape=circle(W * 0.49),
                 params={"a_re": a[0], "a_im": a[1], "b_re": a[2], "b_im": a[3],
                         "c_re": c[0], "c_im": c[1], "d_re": c[2], "d_im": c[3]},
                 label="自定义琼斯元件：%s" % lbl)], src_pol=src_pol)
        cfg = build_config(
            name="自定义琼斯元件：%s" % lbl,
            size=512, width=W, wl=wl, components=comps, sources=srcs, polarized=True,
            category="偏振", difficulty="高级",
            description="直接给出 2×2 复琼斯矩阵 [[a b],[c d]] = [[%s],[%s]]，作用于 [Ex,Ey]ᵀ。"
                        % (a, c),
            expected="%s 理论透过率 %.2f。" % (note, ideal),
            tags=["自定义琼斯", "2×2 矩阵", "偏振光学元件"],
            theory="E_out = J·E_in，J 可为任意复矩阵（不必酉）",
        )
        reg.add("pol_custom_jones_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"))


# ==========================================================================
# 8. 偏振干涉（正交偏振 + 检偏器）
# ==========================================================================

def polarization_interference(reg):
    for i, adeg in enumerate([0, 30, 45, 60, 90]):
        W, wl = W10, WL["633"]
        A = 0.08
        comps = [
            beam_splitter("bs1", (0, 0, 0), reflectivity=0.5, aperture=W * 0.35, label="分束器 1"),
            _pol("pA", 0.04, 0.0, "臂 A：偏振片 0°", r=W * 0.4),
            mirror("m1", (0, 0, A), reflectivity=1.0, yaw=-PI / 4, shape=circle(W * 0.35), label="折转镜 A"),
            _pol("pB", 0.04, PI / 2, "臂 B：偏振片 90°", r=W * 0.4),
            mirror("m2", (A, 0, 0), reflectivity=1.0, yaw=3 * PI / 4, shape=circle(W * 0.35), label="折转镜 B"),
            beam_splitter("bs2", (A, 0, A), reflectivity=0.5, aperture=W * 0.35, label="合束器"),
            _pol("ana", A + 0.02, math.radians(adeg), "检偏器 %.0f°" % adeg, r=W * 0.4),
            detector("d1", (A, 0, 2 * A + 0.02), shape=circle(W * 0.4), yaw=PI, label="端口 1"),
        ]
        cfg = build_config(
            name="偏振干涉：两臂正交偏振 + 检偏器 %.0f°" % adeg,
            size=512, width=W, wl=wl, components=comps,
            sources=[plane_source(wl=wl, label="x 偏振平面波")], polarized=True,
            category="偏振", difficulty="高级",
            description="马赫-曾德尔的两条臂分别放 0° 与 90° 偏振片，使两臂偏振**正交**——"
                        "正交偏振不能干涉（经典意义下）。在合束后加一个 %.0f° 检偏器把两者投影到同一方向，"
                        "干涉条纹重新出现。" % adeg,
            expected="检偏器角度决定可见度：0° 或 90° 时只看到一臂（无干涉、均匀场）；45° 时两臂等权，"
                     "干涉对比度最大。这说明「正交偏振光不干涉」是因为没有把它们投影到同一偏振态上——"
                     "用检偏器投影后干涉恢复到 (1+cosφ)/2 的规律。" if adeg == 45 else
                     "检偏器 %.0f° 使两臂投影振幅分别为 cos%.0f° 与 sin%.0f°，可见度随之变化；"
                     "45° 时两臂等权、对比度最大。" % (adeg, adeg, adeg),
            tags=["偏振干涉", "正交偏振", "检偏器", "可见度"],
            theory="正交偏振态的相干叠加需先投影到同一基；可见度 ∝ |cosθ·sinθ|", sample=(i == 2),
        )
        reg.add("pol_interference_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))


# ==========================================================================
# 9. 高 NA 矢量聚焦（|Ez| 纵场）
# ==========================================================================

def vectorial_focus(reg):
    for i, NA in enumerate([0.3, 0.5, 0.7, 0.9]):
        W, wl, f = W10, WL["633"], 0.3
        D = 2 * W * 0.45
        R = D / 2
        r1 = airy_radius(wl, f, D)
        N = pick_grid(r1, W, minimum=3.0)
        cfg = build_config(
            name="高 NA 矢量聚焦：NA=%.1f（纵场 Ez）" % NA,
            size=N, width=W, wl=wl,
            components=[
                lens("l", (0, 0, 0), f=f, aperture=R, label="高 NA 透镜 f=0.3 m"),
                detector("det", f, label="焦点（可看 |Ez|²）",
                         params={"strehl_aperture": R, "strehl_distance": f}),
            ],
            sources=[gauss_source(wl * 0.6 / (2 * NA) if NA > 0.2 else 1e-3, wl=wl,
                                  polarization="x", label="x 偏振高斯光")],
            method="vectorial",
            category="偏振", difficulty="高级",
            description="用**矢量角谱法**（method=vectorial）模拟 x 偏振光被高 NA 透镜（NA=%.1f）聚焦："
                        "矢量模型会由散度条件重构纵分量 Ez 并把三分量一起传播。" % NA,
            expected="焦斑沿 x 方向被压缩（偏振方向的聚焦更紧）、出现**纵向场 Ez**（在 field=ez 视图可见）；"
                     "NA 越大纵场越强（NA→1 时 |Ez|² 的峰值可与横向相当）。艾里第一暗环 %.1f µm。"
                     "标量模型完全看不到 Ez——这正是高 NA 聚焦必须用矢量理论的原因（超分辨、光镊、"
                     "激光加工都要考虑它）。" % (r1 * 1e6),
            tags=["矢量衍射", "高NA", "纵场Ez", "角谱法"],
            theory="矢量角谱：由 ∇·E=0 解出 Ez，三分量独立传播", sample=(i == 0),
        )
        reg.add("pol_vectorial_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))


# ==========================================================================
# 10. 斯托克斯参数场景（多种偏振态并排）
# ==========================================================================

def stokes_scenes(reg):
    specs = [("线偏振 45°", "d", []),
             ("圆偏振（右旋）", "r", []),
             ("椭圆偏振（QWP 30°）", "x", [("q", PI / 2, math.radians(30))]),
             ("经旋光片 30° 后的线偏振", "x", [("rot", 0.5235987755982988, math.radians(30))])]
    for i, (lbl, spol, extra) in enumerate(specs):
        W, wl = W10, WL["633"]
        els = []
        for k, (kind, val, arg) in enumerate(extra):
            if kind == "q":
                els.append(_ret("q", 0.03, val, arg, "QWP 快轴 %.0f°" % math.degrees(arg), r=W * 0.49))
            else:
                els.append(comp("xrot", "rotator", (0, 0, 0.03), shape=circle(W * 0.49),
                                params={"angle": arg}, label="旋光片 %.0f°" % math.degrees(arg)))
        comps, srcs = pol_stage(W, els, label="偏振读数面", src_pol=spol)
        cfg = build_config(
            name="斯托克斯参数读数：%s" % lbl,
            size=512, width=W, wl=wl, components=comps, sources=srcs, polarized=True,
            category="偏振", difficulty="高级",
            description="生成 %s，在探测面上用 /api/runs/{id}/inspect/{pid} 读取该像素的琼斯矢量、"
                        "斯托克斯参数 (s0,s1,s2,s3)、偏振方位角、椭率与偏振度。" % lbl,
            expected="各偏振态的斯托克斯参数：线偏振 45° → (1,0,±1,0)（s3=0）；圆偏振 → (1,0,0,±1)；"
                     "椭圆偏振 → s3 介于 0 与 ±1 之间；自然光/部分偏振 → s0²>s1²+s2²+s3²。"
                     "方位角 ψ=½atan(s2/s1)，椭率 χ=½asin(s3/s0)。GUI 在偏振视图叠加椭圆阵列。",
            tags=["斯托克斯参数", "偏振椭圆", "偏振度"],
            theory="S=(s0,s1,s2,s3) 由琼斯矢量算得；纯态 s0²=s1²+s2²+s3²", sample=(i < 2),
        )
        reg.add("pol_stokes_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))

    # 部分偏振：两个不相干光源（不同波长）合成
    for i, (n_src, lbl) in enumerate([(1, "单一偏振态（完全偏振）"), (2, "两个不同波长的正交偏振（部分偏振）")]):
        W, wl = W10, WL["633"]
        comps = [aperture("clip", (0, 0, 0), rectangle(W * 0.6, W * 0.6), label="限光阑"),
                 detector("det", 0.02, shape=rectangle(W * 0.6, W * 0.6), label="偏振读数面")]
        srcs = [plane_source(wl=wl, power=0.5e-3 if n_src == 2 else 1e-3, sid="s1",
                             label="633 nm（x 偏振）")]
        if n_src == 2:
            srcs.append(source("s2", "plane", {"power": 0.5e-3, "polarization": "y"},
                               (0.0, 0.0, SRC_Z), wavelength=WL["532"],
                               group="B", label="532 nm（y 偏振，不相干）"))
        cfg = build_config(
            name="偏振度：%s" % lbl,
            size=512, width=W, wl=wl, components=comps, sources=srcs, polarized=True,
            category="偏振", difficulty="高级",
            description="%s：两个不同波长的光源属于互不相干的相干单元，合成后是**部分偏振**光。"
                        "用 inspect 读偏振度 DOP=s1²+s2²+s3² 开方/s0。" % lbl,
            expected=("完全偏振：DOP=1。" if n_src == 1 else
                      "部分偏振：DOP<1（两个正交偏振、互不相干的分量叠加，"
                      "s1、s2 部分抵消）。这是「偏振度」的实际意义，也是天空散射光、"
                      "多波长激光合成的真实情形。"),
            tags=["偏振度", "部分偏振", "相干分组"],
            theory="多相干单元的斯托克斯参数按强度相加 → DOP 可 <1", sample=True,
        )
        reg.add("pol_dop_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=True)
