"""介质类示例（med_*）—— 目标 40+。

均匀介质的折射率（相位）与吸收/增益、split-step 子步收敛、单轴晶体双折射、
双轴晶体（Berreman 4×4）、克尔非线性、可饱和吸收体、自聚焦。
"""

import math

from wavegen.common import *   # noqa: F401,F403
from wavegen.fam_interf import mz_comps, michelson_comps, d_port_z, d_port_nx

W10 = 0.01
W8 = 0.008
PI = math.pi


def _medium(cid, z, d, index=1.5, absorption=0.0, steps=20, r=0.004, yaw=0.0, label=None):
    p = {"distance": float(d), "index": float(index), "absorption": float(absorption),
         "steps": int(steps)}
    return comp(cid, "medium", (0, 0, z), yaw=yaw, shape=circle(r), params=p, label=label)


def _abs_n(n_w, d, wl=None):
    """把吸收系数换算成复折射率虚部：I/I0 = exp(-2 k n'' d)。"""
    wl = wl or WL["633"]
    return -math.log(max(1e-12, n_w)) / (2 * k_of(wl) * d)


def _phase_index(phase, d, wl=None):
    """要让介质产生 phase 的额外相位（k(n−1)d），所需的折射率。"""
    wl = wl or WL["633"]
    return 1.0 + phase / (k_of(wl) * d)


# ==========================================================================
# 1. 均匀介质：折射率（相位）
# ==========================================================================

def uniform_index(reg):
    d, wl = 0.01, WL["633"]
    for i, frac in enumerate([0.0, 0.25, 0.5, 0.75, 1.0, 1.25]):
        n = _phase_index(frac * 2 * PI, d, wl)
        cfg = build_config(
            name="介质相位板：折射率 n=%.7f（额外相位 %.2f×2π）" % (n, frac),
            size=512, width=W10, wl=wl,
            components=mz_comps(W10, A=0.08,
                                arm_extra=[_medium("plate", 0.04, d, index=n, absorption=0.0,
                                                   steps=8, r=W10 * 0.4,
                                                   label="介质板 d=10 mm、n=%.7f" % n)],
                                d1_label="端口 1（干涉读出相位）"),
            sources=[plane_source(wl=wl)],
            category="介质", difficulty="进阶",
            description="马赫-曾德尔的一臂插入厚度 10 mm 的均匀介质（复折射率 n=%.7f + i·0），"
                        "用干涉读出它引入的相位 k(n−1)d。" % n,
            expected="介质引起的相位 = k(n−1)d = %.2f×2π：端口功率 = cos²(φ/2)，"
                     "因此折射率每变化 Δn=λ/d=%.0f ppm，端口就在亮暗之间走一格——"
                     "这就是干涉式折射率计（灵敏度可达 1e-7 量级）的原理。"
                     % (frac, WL["633"] / d * 1e6),
            tags=["折射率", "相位", "干涉式测量"],
            theory="φ = k(n−1)d；周期 Δn = λ/d", sample=(i < 2),
        )
        reg.add("med_index_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))


# ==========================================================================
# 2. 吸收扫描（直通）：I/I₀ = exp(-2 k n'' d)
# ==========================================================================

def absorption_sweep(reg):
    d, wl = 0.01, WL["633"]
    for i, T in enumerate([0.9, 0.7, 0.5, 0.3, 0.1, 0.01]):
        np_ = _abs_n(T, d, wl)
        alpha = 2 * k_of(wl) * np_
        cfg = build_config(
            name="均匀介质吸收：透过率 %.0f%%（α=%.0f /m）" % (T * 100, alpha),
            size=512, width=W10, wl=wl,
            components=[
                _medium("m", 0.02, d, index=1.5, absorption=np_, steps=10, r=W10 * 0.49,
                        label="吸收介质 d=10 mm、n''=%.3e" % np_),
                detector("det", 0.1, shape=circle(W10 * 0.49), label="透过功率"),
            ],
            sources=[plane_source(wl=wl)],
            category="介质", difficulty="入门",
            description="平面波穿过 10 mm 厚的吸收介质（复折射率实部 1.5、虚部 %.3e）。" % np_,
            expected="比尔-朗伯定律：I/I₀=exp(−2k n''d)=exp(−αd)=%.3f（α=2kn''=%.0f /m）。"
                     "虚部与强度的关系是「2k n''」，因为强度 ∝ |E|²。"
                     "本内核用 split-step 实现：吸收在子步里逐段施加，与解析结果一致。" % (T, alpha),
            tags=["吸收", "比尔-朗伯", "复折射率"],
            theory="I(z) = I₀ exp(−2 k n'' z)", sample=(i < 3),
        )
        reg.add("med_absorption_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))

    # 厚度扫描（同一吸收系数）
    for i, dd in enumerate([0.002, 0.005, 0.01, 0.02, 0.05]):
        np_ = _abs_n(0.5, 0.01, WL["633"])   # 10 mm 透过 50%
        T = math.exp(-2 * k_of(WL["633"]) * np_ * dd)
        cfg = build_config(
            name="吸收介质厚度扫描：d=%.0f mm → T=%.1f%%" % (dd * 1e3, T * 100),
            size=512, width=W10, wl=WL["633"],
            components=[
                _medium("m", 0.02, dd, index=1.5, absorption=np_, steps=10, r=W10 * 0.49,
                        label="吸收介质 d=%.0f mm" % (dd * 1e3)),
                detector("det", 0.1, shape=circle(W10 * 0.49), label="透过功率"),
            ],
            sources=[plane_source(wl=WL["633"])],
            category="介质", difficulty="入门",
            description="吸收系数固定（n''=%.3e，即 α=%.0f /m），只改厚度。" % (np_, 2 * k_of(WL["633"]) * np_),
            expected="透过率按指数随厚度下降：T=exp(−αd)。10 mm 时 50%%，20 mm 时 25%%——"
                     "指数衰减律（朗伯定律）是浓度/厚度测量（分光光度计）的基础。" % (),
            tags=["吸收", "厚度", "朗伯定律"],
            theory="T = exp(−αd)", sample=False,
        )
        reg.add("med_thickness_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"))


# ==========================================================================
# 3. 增益介质（负吸收）
# ==========================================================================

def gain_medium(reg):
    d = 0.01
    for i, G in enumerate([1.2, 1.5, 2.0, 4.0]):
        np_ = -_abs_n(G, d)      # 负虚部 = 增益
        cfg = build_config(
            name="增益介质：小信号增益 %.1f×" % G,
            size=512, width=W10, wl=WL["633"],
            components=[
                _medium("m", 0.02, d, index=1.5, absorption=np_, steps=10, r=W10 * 0.49,
                        label="增益介质 d=10 mm、n''=%.3e（负虚部）" % np_),
                detector("det", 0.1, shape=circle(W10 * 0.49), label="输出功率"),
            ],
            sources=[plane_source(wl=WL["633"], power=1e-4)],
            category="介质", difficulty="进阶",
            description="复折射率虚部取负（n''=%.3e）→ 介质提供**增益**而不是吸收；"
                        "输入功率取 1e-4 W 以避免非物理的无限放大。" % np_,
            expected="输出功率按 exp(+2k|n''|d)=exp(α_g d)=%.1f× 放大（小信号增益）。"
                     "这是激光增益介质（掺铒光纤、半导体光放大器）的最简模型；"
                     "真实器件还有增益饱和（见可饱和吸收体/克尔元件）。" % G,
            tags=["增益", "负吸收", "光放大"],
            theory="负虚部 → exp(+2k|n''|d) 放大（本模型不含饱和）",
        )
        reg.add("med_gain_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"))


# ==========================================================================
# 4. split-step 子步收敛
# ==========================================================================

def splitstep_convergence(reg):
    d, wl = 0.01, WL["633"]
    for i, steps in enumerate([1, 2, 5, 10, 40]):
        cfg = build_config(
            name="split-step 子步收敛：steps=%d" % steps,
            size=512, width=W10, wl=wl,
            components=[
                lens("l", (0, 0, 0), f=0.15, aperture=W10 * 0.4, label="聚焦透镜 f=150 mm"),
                _medium("m", 0.05, d, index=1.5, absorption=0.0, steps=steps, r=W10 * 0.49,
                        label="介质（split-step 子步 %d）" % steps),
                detector("det", 0.15, label="焦点"),

            ],
            sources=[gauss_source(0.6e-3, wl=wl, label="高斯 w₀=0.6 mm")],
            category="介质", difficulty="高级",
            description="聚焦光束穿过均匀介质，改变 split-step 的子步数 steps=%d：子步越少，"
                        "「先衍射后加相位」的算子分裂误差越大。" % steps,
            expected="均匀介质里没有非线性/不均匀折射率，理论上任意 steps 都给出同一结果"
                     "（介质只是把 k·n·d 的相位与吸收施加到场上）——本实例用来验证算子分裂的收敛："
                     "steps=%d 的焦斑峰值与 steps=40 的差异应在数值精度内。"
                     "当介质换成强度相关（克尔）或分层折射率时，子步数才真正决定精度。" % steps,
            tags=["split-step", "BPM", "算子分裂", "收敛"],
            theory="对称分步傅里叶：exp(−iH₁h/2)exp(−iH₂h)exp(−iH₁h/2)", sample=(i < 2),
        )
        reg.add("med_splitstep_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))


# ==========================================================================
# 5. 单轴晶体（双折射）
# ==========================================================================

def uniaxial_sweep(reg):
    wl, d = WL["633"], 0.002
    for i, (no, ne) in enumerate([(1.5, 1.5), (1.5, 1.55), (1.5, 1.6), (1.5, 1.65),
                                  (1.5, 1.7), (1.5, 1.9)]):
        dn = ne - no
        dphi = k_of(wl) * dn * d
        N = 512
        cfg = build_config(
            name="单轴晶体双折射：n_e−n_o=%.3f（相位差 %.2f×2π）" % (dn, dphi / (2 * PI)),
            size=N, width=W10, wl=wl,
            components=[
                comp("p45", "polarizer", (0, 0, 0), shape=circle(W10 * 0.49),
                     params={"angle": PI / 4}, label="起偏器 45°"),
                comp("xtal", "uniaxial", (0, 0, 0.03), shape=circle(W10 * 0.49),
                     params={"distance": d, "n_o": no, "n_e": ne},
                     label="单轴晶体（光轴∥x）d=2 mm"),
                comp("ana", "polarizer", (0, 0, 0.06), shape=circle(W10 * 0.49),
                     params={"angle": PI / 2 + PI / 4}, label="检偏器 135°（正交）"),
                detector("det", 0.08, shape=circle(W10 * 0.49), label="干涉色读数"),
            ],
            sources=[source("src", "plane", {"power": 1e-3, "polarization": "d"},
                            (0.0, 0.0, SRC_Z), wavelength=wl, label="45° 线偏振平面波")],
            polarized=True,
            category="介质", difficulty="进阶",
            description="45° 线偏振光穿过 2 mm 单轴晶体（光轴沿 x）：o 光（Ey）与 e 光（Ex）感受到不同折射率，"
                        "出射后由正交的 135° 检偏器读出两者干涉。",
            expected="两分量相位差 δ=k(n_e−n_o)d=%.2f rad=%.2f×2π；正交检偏器下透过率 = sin²(δ/2)。"
                     "δ 随 Δn 线性变化：Δn 每变化 λ/d=%.0f ppm，干涉色走一整圈——"
                     "这就是偏光显微镜、相位延迟片与电光调制器的工作基础。"
                     % (dphi, dphi / (2 * PI), wl / d * 1e6),
            tags=["双折射", "单轴晶体", "o光e光", "相位延迟"],
            theory="δ = k(n_e−n_o)d；正交偏光下 I = I₀sin²(δ/2)", sample=(i < 3),
        )
        reg.add("med_uniaxial_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))


def uniaxial_thickness(reg):
    wl, no, ne = WL["633"], 1.5, 1.62
    for i, d in enumerate([0.1e-3, 0.25e-3, 0.5e-3, 0.75e-3, 1.0e-3]):
        dphi = k_of(wl) * (ne - no) * d
        cfg = build_config(
            name="双折射厚度扫描：d=%.2f mm → δ=%.2f×2π" % (d * 1e3, dphi / (2 * PI)),
            size=512, width=W10, wl=wl,
            components=[
                comp("p45", "polarizer", (0, 0, 0), shape=circle(W10 * 0.49),
                     params={"angle": PI / 4}, label="起偏器 45°"),
                comp("xtal", "uniaxial", (0, 0, 0.03), shape=circle(W10 * 0.49),
                     params={"distance": d, "n_o": no, "n_e": ne},
                     label="单轴晶体 d=%.2f mm" % (d * 1e3)),
                comp("ana", "polarizer", (0, 0, 0.06), shape=circle(W10 * 0.49),
                     params={"angle": 3 * PI / 4}, label="检偏器 135°（正交）"),
                detector("det", 0.08, shape=circle(W10 * 0.49), label="读数"),
            ],
            sources=[source("src", "plane", {"power": 1e-3, "polarization": "d"},
                            (0.0, 0.0, SRC_Z), wavelength=wl, label="45° 线偏振")],
            polarized=True,
            category="介质", difficulty="进阶",
            description="Δn 固定 0.12，扫描晶体厚度——等价于「一块晶体切多厚就能当几分之一波片」。",
            expected="相位差 δ=2πΔn·d/λ：厚度 %.2f mm 时 δ=%.3f×2π。"
                     "δ=0.5×2π 时相当于半波片（正交偏光下全透），δ=0.25×2π 时相当于四分之一波片"
                     "（正交偏光下透过 50%%）。石英/云母波片就是按厚度分级挑选的。" % (d * 1e3, dphi / (2 * PI)),
            tags=["波片", "厚度分级", "双折射"],
            theory="δ ∝ d；HWP: δ=π，QWP: δ=π/2", sample=(i < 2),
        )
        reg.add("med_uniaxial_thick_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))


# ==========================================================================
# 6. 双轴晶体（Berreman 4×4）
# ==========================================================================

def biaxial_sweep(reg):
    specs = [(1.60, 1.58, 1.56), (1.70, 1.60, 1.50), (1.80, 1.55, 1.45),
             (1.55, 1.60, 1.65), (2.00, 1.60, 1.40)]
    for i, (nx, ny, nz) in enumerate(specs):
        wl, d = WL["633"], 0.001
        cfg = build_config(
            name="双轴晶体（Berreman 4×4）：n=(%.2f, %.2f, %.2f)" % (nx, ny, nz),
            size=512, width=W10, wl=wl,
            components=[
                comp("p45", "polarizer", (0, 0, 0), shape=circle(W10 * 0.49),
                     params={"angle": PI / 4}, label="起偏器 45°"),
                comp("xtal", "biaxial", (0, 0, 0.03), shape=circle(W10 * 0.49),
                     params={"distance": d, "n_x": nx, "n_y": ny, "n_z": nz},
                     label="双轴晶体 d=1 mm"),
                comp("ana", "polarizer", (0, 0, 0.06), shape=circle(W10 * 0.49),
                     params={"angle": 3 * PI / 4}, label="检偏器 135°"),
                detector("det", 0.08, shape=circle(W10 * 0.49), label="读数"),
            ],
            sources=[source("src", "plane", {"power": 1e-3, "polarization": "d"},
                            (0.0, 0.0, SRC_Z), wavelength=wl, label="45° 线偏振")],
            polarized=True,
            category="介质", difficulty="高级",
            description="三主轴折射率互不相等的**双轴晶体**（n_x=%.2f, n_y=%.2f, n_z=%.2f），"
                        "用完整 Berreman 4×4 各向异性传播求解（含 o/e 耦合，比单轴模型慢）。" % (nx, ny, nz),
            expected="双轴晶体的两个本征偏振态都不是纯 o/e 光，偏振态随传播连续演化（含 o↔e 能量交换），"
                     "正交检偏器下的透过率不再是简单的 sin²(kΔnd/2) 而是 4×4 特征值问题的结果。"
                     "本实例 Δn 约 %.3f，δ≈%.2f rad——与同 Δn 的单轴晶体相比，透过率会略有差别。"
                     % (max(nx, ny, nz) - min(nx, ny, nz), k_of(wl) * (max(nx, ny, nz) - min(nx, ny, nz)) * d),
            tags=["双轴晶体", "Berreman", "4×4", "各向异性"],
            theory="Berreman 4×4 传递矩阵：求解 [Ex,Ey,Hx,Hy] 的本征模",
        )
        reg.add("med_biaxial_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"))


# ==========================================================================
# 7. 克尔介质（强度相关折射率）
# ==========================================================================

def kerr_effects(reg):
    wl = WL["633"]
    for i, (n2, L) in enumerate([(1e-6, 1e-3), (1e-6, 5e-3), (5e-6, 1e-3),
                                 (1e-6, 1e-2), (2e-6, 2e-3), (-1e-6, 1e-3)]):
        cfg = build_config(
            name="克尔介质：n₂=%.0e m²/W、L=%.0f mm" % (n2, L * 1e3),
            size=512, width=W8, wl=wl,
            components=[
                comp("kerr", "kerr", (0, 0, 0.02), shape=circle(W8 * 0.49),
                     params={"n2": n2, "length": L},
                     label="克尔片 n₂=%.0e m²/W L=%.0f mm" % (n2, L * 1e3)),
                detector("det", 0.3, label="远场"),
            ],
            sources=[gauss_source(0.6e-3, wl=wl, power=1e-3, label="高斯 w₀=0.6 mm、1 mW")],
            category="介质", difficulty="高级",
            description="薄克尔片：折射率 n=n₀+n₂·I 随局部强度变化，因此中心（强度高）与边缘（强度低）"
                        "获得不同相位 φ=k·n₂·I·L —— 光束自己给自己写了一块透镜。",
            expected="相位剖面 φ(r)=k·n₂·I(r)·L，中心峰值相位 %.3f rad（I_peak≈%.3g W/m²，"
                     "n₂ 为「有效值」：真实介质 n₂~1e-14–1e-20 m²/W，需要 GW/cm² 级峰值功率）。"
                     "n₂>0 时中心相位滞后最多 → 等效**自聚焦**透镜（远场光斑收缩）；"
                     "n₂<0 则自散焦。这是 Kerr 透镜模式锁定、自聚焦成丝与光限幅的物理内核。"
                     % (k_of(wl) * n2 * (2 * 1e-3 / (PI * (0.6e-3) ** 2)) * L,
                        2 * 1e-3 / (PI * (0.6e-3) ** 2)),
            tags=["克尔效应", "自聚焦", "非线性", "n₂"],
            theory="φ = k·n₂·I·L；B 积分 = k∫n₂I dz", sample=(i < 2),
        )
        reg.add("med_kerr_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))

    # 克尔介质 + 双光子吸收
    for i, tpa in enumerate([0.0, 1e-9, 1e-8]):
        cfg = build_config(
            name="克尔介质 + 双光子吸收：β=%.0e m/W" % tpa,
            size=512, width=W8, wl=WL["633"],
            components=[
                comp("kerr", "kerr", (0, 0, 0.02), shape=circle(W8 * 0.49),
                     params={"n2": 2e-6, "length": 2e-3, "tpa": tpa},
                     label="克尔片（含双光子吸收 β=%.0e m/W）" % tpa),
                detector("det", 0.3, label="远场"),
            ],
            sources=[gauss_source(0.6e-3, wl=WL["633"], power=1e-3, label="高斯 w₀=0.6 mm")],
            category="介质", difficulty="高级",
            description="克尔片再加**双光子吸收**（β=%.0e m/W）：高强度区被额外吸收，"
                        "与自聚焦的「光强越集中越聚焦」形成竞争。" % tpa,
            expected="β=0 时只有自聚焦；β>0 时中心衰减、光斑剖面变平甚至出现「反饱和」形变"
                     "（强度越高吸收越强 → 峰值被削平）。双光子吸收是自相关仪、光限幅器与"
                     "双光子显微镜的基础；它与克尔自聚焦的竞争决定成丝阈值。",
            tags=["双光子吸收", "克尔", "光限幅"],
            theory="dI/dz = −βI²（双光子）+ n₂I 相位", sample=(i == 0),
        )
        reg.add("med_kerr_tpa_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))


# ==========================================================================
# 8. 可饱和吸收体
# ==========================================================================

def saturable_sweep(reg):
    for i, (a0, isat, L) in enumerate([(500, 1e6, 1e-3), (2000, 1e6, 1e-3), (2000, 1e7, 1e-3),
                                       (5000, 1e6, 1e-3), (2000, 1e5, 1e-3)]):
        cfg = build_config(
            name="可饱和吸收体：α₀=%.0f /m、I_sat=%.0e W/m²" % (a0, isat),
            size=512, width=W8, wl=WL["633"],
            components=[
                comp("sa", "saturable_absorber", (0, 0, 0.02), shape=circle(W8 * 0.49),
                     params={"alpha0": a0, "isat": isat, "length": L},
                     label="可饱和吸收体 α₀=%.0f /m" % a0),
                detector("det", 0.3, label="输出"),
            ],
            sources=[gauss_source(0.6e-3, wl=WL["633"], power=2e-3, label="高斯 w₀=0.6 mm、2 mW")],
            category="介质", difficulty="高级",
            description="可饱和吸收体 α(I)=α₀/(1+I/I_sat)：弱光被强吸收、强光几乎透明。"
                        "输入为高斯光束，中心强、边缘弱。" % (),
            expected="输出剖面的中心透过率高于边缘（中心 exp(−α₀L/(1+I₀/I_sat))）："
                     "弱光极限透过率 exp(−α₀L)=%.3f，中心（I₀≈%.2g W/m²）处透过率提高到 %.3f。"
                     "结果是一个**中心鼓起的透过率剖面**——这就是被动调 Q / 可饱和吸收被动锁模"
                     "（激光脉冲的自动整形与短脉冲产生）的核心机制。"
                     % (math.exp(-a0 * L), 2e-3 / (PI * 0.6e-3 ** 2),
                        math.exp(-a0 * L / (1 + (2e-3 / (PI * 0.6e-3 ** 2)) / isat))),
            tags=["可饱和吸收", "被动调Q", "锁模"],
            theory="α(I) = α₀/(1+I/I_sat)：强度越高吸收越低", sample=(i < 2),
        )
        reg.add("med_saturable_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))


# ==========================================================================
# 9. 自聚焦（克尔 + 传播）
# ==========================================================================

def self_focusing(reg):
    for i, (n2, L) in enumerate([(2e-6, 5e-3), (5e-6, 10e-3)]):
        cfg = build_config(
            name="克尔自聚焦：n₂=%.0e、L=%.0f mm 后的焦点" % (n2, L * 1e3),
            size=512, width=W8, wl=WL["633"],
            components=[
                comp("kerr", "kerr", (0, 0, 0.0), shape=circle(W8 * 0.49),
                     params={"n2": n2, "length": L},
                     label="克尔介质 n₂=%.0e m²/W、等效长度 %.0f mm" % (n2, L * 1e3)),
                detector("det", 0.2, label="远场（自聚焦后）"),
            ],
            sources=[gauss_source(0.8e-3, wl=WL["633"], power=2e-3, label="高斯 w₀=0.8 mm")],
            category="介质", difficulty="高级",
            description="高斯光束穿过克尔介质：中心相位滞后（n₂>0）形成等效自聚焦透镜，"
                        "传播 0.2 m 后观察光斑收缩。" % (),
            expected="远场光斑小于同等条件下的线性高斯（衍射）光斑：非线性焦点距离 z_f≈w₀/√(2n₂I₀/n₀)（B 积分模型）。"
                     "峰值强度 I₀≈%.2g W/m²、非线性相位 B=%.2f rad —— B≈1–3 时自聚焦最明显；"
                     "B 过大时会出现自聚焦-衍射的周期振荡（成丝）。"
                     % (2e-3 / (PI * 0.8e-3 ** 2), k_of(WL["633"]) * n2 * (2e-3 / (PI * 0.8e-3 ** 2)) * L),
            tags=["自聚焦", "克尔", "B积分", "成丝"],
            theory="自聚焦焦距 z_f ≈ w₀/sqrt(2 n₂ I / n₀)；临界功率 P_cr=3.77λ²/(8πn₂)",
        )
        reg.add("med_self_focusing_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"))

    # 介质中的衍射：菲涅耳数与介质折射率（波长在介质中变短）
    for i, n in enumerate([1.0, 1.33, 1.5, 1.8]):
        wl = WL["633"]
        d = 0.02
        cfg = build_config(
            name="介质中的衍射：n=%.2f（介质内波长 %.0f nm）" % (n, wl / n * 1e9),
            size=512, width=W10, wl=wl,
            components=[
                _medium("m", 0.0, d, index=n, absorption=0.0, steps=20, r=W10 * 0.49,
                        label="介质 d=20 mm、n=%.2f" % n),
                aperture("slit", (0, 0, 0.001), rectangle(0.3e-3, W10 * 0.95), label="单缝 0.3 mm"),
                detector("det", 1.5, label="远场（介质之后）"),
            ],
            sources=[plane_source(wl=wl)],
            category="介质", difficulty="高级",
            description="先把平面波送进 20 mm 厚的介质（n=%.2f）再加单缝衍射——介质内的波长为 λ/n=%.0f nm，"
                        "出射后回到空气中的 λ。" % (n, wl / n * 1e9),
            expected="介质内传播只改变相位（k·n·d）与吸收，不改变场的横向分布，因此出射后的**衍射图样与"
                     "真空完全一致**（第一暗点仍为 λz/a=%.2f mm，λ 用真空波长）。"
                     "这与「介质中波长缩短、衍射角变小」的常见误解形成对照：衍射角由光在**外部介质**"
                     "（出射侧）中的波长决定，只要出射到空气就是 λ。" % (sinc_first_zero(wl, 1.5, 0.3e-3) * 1e3),
            tags=["介质中的衍射", "波长", "相位"],
            theory="均匀介质只乘 exp(iknd) 与吸收；横向分布不变",
        )
        reg.add("med_diffraction_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=(i == 0))
