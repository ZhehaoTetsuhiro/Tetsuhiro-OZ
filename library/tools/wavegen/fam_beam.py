"""传播/光束类示例（beam_*）—— 目标 70+。

高斯光束的束腰演化/发散/瑞利长度/古伊相位、LG 与 HG 高阶模、贝塞尔（无衍射）光束、
球面波、倾斜平面波位移、涡旋光束、透镜对高斯光束的变换（束腰变换公式）。
"""

import math

from wavegen.common import *   # noqa: F401,F403

W8 = 0.008
W20 = 0.02
W10 = 0.01


def _det_plane(cid="det", z=0.5):
    return detector(cid, z, label="观察面")


# ==========================================================================
# 高斯光束：束腰扫掠 —— w(z) = w0√(1+(z/z_R)²)
# ==========================================================================

def gauss_waist_sweep(reg):
    z = 0.5
    for i, w0 in enumerate([0.1e-3, 0.2e-3, 0.3e-3, 0.5e-3, 0.8e-3, 1.2e-3, 2.0e-3]):
        wz = waist_at(WL["633"], w0, z)
        zR = rayleigh_range(WL["633"], w0)
        W = pick_window(2.2 * wz)
        N = 512 if W <= 0.01 else 1024
        dx = W / N
        check_sampling(dx, {"束腰 w0": w0}, 4.0)
        window_ok(W, {"z 处半径": 2.2 * wz})
        cfg = build_config(
            name="高斯光束传播：w₀=%.2f mm（0.5 m 后展宽 %.2f mm）" % (w0 * 1e3, wz * 1e3),
            size=N, width=W, wl=WL["633"],
            components=[_det_plane(z=z)],
            sources=[gauss_source(w0, wl=WL["633"], label="高斯 w₀=%.2f mm" % (w0 * 1e3))],
            category="传播/光束", difficulty="入门",
            description="基模高斯光束（束腰 %.2f mm）自由传播 %.2f m。瑞利长度 z_R=πw₀²/λ=%.3f m。"
                        % (w0 * 1e3, z, zR),
            expected="束腰处最窄，之后按 w(z)=w₀√(1+(z/z_R)²) 展宽到 %.2f mm（为 w₀ 的 %.2f 倍）；"
                     "束腰越小发散越快（θ=λ/πw₀=%.2f mrad）——这是「小光斑 ⇔ 大发散角」的测不准关系式表达。"
                     % (wz * 1e3, wz / w0, gauss_divergence(WL["633"], w0) * 1e3),
            tags=["高斯光束", "瑞利长度", "束腰演化"],
            theory="w(z)=w₀√(1+(z/z_R)²)，z_R=πw₀²/λ", sample=(i < 2),
        )
        reg.add("beam_gauss_waist_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))


# ==========================================================================
# 高斯光束：距离扫掠
# ==========================================================================

def gauss_z_sweep(reg):
    w0 = 0.3e-3
    zR = rayleigh_range(WL["633"], w0)
    for i, z in enumerate([0.05, 0.1, 0.2, 0.45, 0.8, 1.5]):
        wz = waist_at(WL["633"], w0, z)
        W = pick_window(2.2 * wz)
        N = 512 if W <= 0.01 else 1024
        dx = W / N
        check_sampling(dx, {"束腰": w0}, 8.0)
        window_ok(W, {"z 处半径": 2.2 * wz})
        cfg = build_config(
            name="高斯光束距离扫掠：z=%.2f m（z/z_R=%.2f）" % (z, z / zR),
            size=N, width=W, wl=WL["633"],
            components=[_det_plane(z=z)],
            sources=[gauss_source(w0, wl=WL["633"], label="高斯 w₀=0.30 mm")],
            category="传播/光束", difficulty="入门",
            description="束腰 0.30 mm 的高斯光束，观察面从 %.0f mm 移到 %.0f mm（z_R=%.3f m）。"
                        % (0.05 * 1e3, z * 1e3, zR),
            expected="z≪z_R 时几乎不展宽（准直区），z≫z_R 后按 z·θ=z·λ/(πw₀) 线性发散："
                     "本处 z/z_R=%.2f，光斑半径 %.2f mm。远场半角 θ=%.2f mrad。" 
                     % (z / zR, wz * 1e3, gauss_divergence(WL["633"], w0) * 1e3),
            tags=["高斯光束", "准直区", "远场发散"],
            theory="近场准直 / 远场线性发散，转折点即 z_R", sample=(i < 2),
        )
        reg.add("beam_gauss_z_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))


# ==========================================================================
# 古伊相位
# ==========================================================================

def gouy_phase(reg):
    w0 = 0.3e-3
    zR = rayleigh_range(WL["633"], w0)
    for i, z in enumerate([0.1, 0.45, 1.0, 2.0]):
        N, W = 512, W10
        dx = W / N
        g = gouy(WL["633"], w0, z)
        cfg = build_config(
            name="高斯光束的古伊相位：z/z_R=%.2f → ζ=%.1f°" % (z / zR, math.degrees(g)),
            size=N, width=W, wl=WL["633"],
            components=[_det_plane(z=z)],
            sources=[gauss_source(w0, wl=WL["633"], label="高斯 w₀=0.30 mm")],
            category="传播/光束", difficulty="高级",
            description="束腰 0.30 mm 的高斯光束传播到 z=%.2f m；用 phase_u（二维去包裹）视图读波前，"
                        "观察**古伊相位**——高斯光束相对平面波多出的轴向相位。" % z,
            expected="弯折的等相位面在轴上累积 ζ(z)=atan(z/z_R)=%.1f°（本处）。z 从 0 到 ∞ 时 ζ 从 0 单调增到 90°；"
                     "这是高斯光束与平面波的全部相位差（焦点附近波前是平的、远场是球面的）。" % math.degrees(g),
            tags=["古伊相位", "波前", "相位分析"],
            theory="ζ(z)=atan(z/z_R)，总古伊相移 π（各半面 π/2）",
        )
        reg.add("beam_gouy_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"))


# ==========================================================================
# LG 高阶模
# ==========================================================================

def lg_modes(reg):
    combos = [(0, 1), (0, 2), (0, 3), (0, 4), (0, 5),
              (1, 0), (1, 1), (1, 2), (1, 3),
              (2, 0), (2, 1), (2, 2),
              (0, -1), (0, -3), (1, -2)]
    for i, (p, l) in enumerate(combos):
        w0, N, W = 0.35e-3, 512, W8
        dx = W / N
        check_sampling(dx, {"束腰": w0}, 6.0)
        n_rings = p + 1
        cfg = build_config(
            name="拉盖尔-高斯模 LG(%d,%d)：%d 个同心环 + 螺旋相位 %d" % (p, l, n_rings, l),
            size=N, width=W, wl=WL["633"],
            components=[_det_plane(z=0.3)],
            sources=[lg_source(w0, p, l, wl=WL["633"], label="LG(%d,%d)" % (p, l))],
            category="传播/光束", difficulty="进阶",
            description="拉盖尔-高斯模 LG_p^l（p=%d 径向阶、l=%d 方位角阶），束腰 0.35 mm，传播 0.3 m。"
                        % (p, l),
            expected="强度呈 %d 个同心环；l≠0 时轴上强度为零、相位绕轴 %d 圈 2π（轨道角动量 lħ/光子）；"
                     "l 的符号决定相位螺旋旋向（l>0 与 l<0 的图样在**强度**上相同、在**相位**上相反）。"
                     % (n_rings, abs(l)),
            tags=["拉盖尔-高斯", "OAM", "高阶模"],
            theory="LG_p^l 有 p+1 个径向节环，相位因子 e^{ilθ}", sample=(i < 2),
        )
        reg.add("beam_lg_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))


# ==========================================================================
# HG 高阶模
# ==========================================================================

def hg_modes(reg):
    combos = [(1, 0), (2, 0), (3, 0), (4, 0), (5, 0),
              (0, 1), (0, 2), (0, 3),
              (1, 1), (1, 2), (2, 1), (2, 2), (3, 3), (4, 4)]
    for i, (m, n) in enumerate(combos):
        w0, N, W = 0.35e-3, 512, W8
        dx = W / N
        check_sampling(dx, {"束腰": w0}, 6.0)
        cfg = build_config(
            name="厄米-高斯模 HG(%d,%d)：%d×%d 个亮瓣" % (m, n, m + 1, n + 1),
            size=N, width=W, wl=WL["633"],
            components=[_det_plane(z=0.3)],
            sources=[hg_source(w0, m, n, wl=WL["633"], label="HG(%d,%d)" % (m, n))],
            category="传播/光束", difficulty="进阶",
            description="厄米-高斯模 HG_mn（m=%d、n=%d），束腰 0.35 mm，传播 0.3 m。" % (m, n),
            expected="强度呈 %d×%d 的矩形亮瓣阵列（沿 x 有 %d 个节点、沿 y 有 %d 个节点），"
                     "整体被高斯包络限制；阶数越高光束越宽（M²=2m+1 与 2n+1）。"
                     % (m + 1, n + 1, m, n),
            tags=["厄米-高斯", "高阶模", "M²"],
            theory="HG_mn 由厄米多项式 H_m(x)H_n(y) 调制高斯，节点数 = (m,n)",
            sample=(i < 2),
        )
        reg.add("beam_hg_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))


# ==========================================================================
# 贝塞尔光束（轴锥镜）：无衍射
# ==========================================================================

def bessel_beams(reg):
    for i, alpha in enumerate([0.005, 0.008, 0.01, 0.015, 0.02, 0.03, 0.05]):
        w0, n_ax = 1.5e-3, 1.5
        core = bessel_core_radius(WL["633"], (n_ax - 1) * alpha)
        # 自适应网格：中心斑至少要占 3 个像素
        if core / (0.008 / 1024.0) >= 3.0:
            N, W = 1024, 0.008
        elif core / (0.004 / 1024.0) >= 3.0:
            N, W = 1024, 0.004
        else:
            N, W = 2048, 0.004
        dx = W / N
        check_sampling(dx, {"束腰": w0}, 6.0)
        check_sampling(dx, {"贝塞尔中心斑": core}, 3.0)
        zmax = w0 / math.tan((n_ax - 1) * alpha) if alpha > 0 else 1.0
        cfg = build_config(
            name="贝塞尔光束（轴锥镜 α=%.1f mrad）：中心斑 %.1f µm" % (alpha * 1e3, core * 1e6),
            size=N, width=W, wl=WL["633"],
            components=[
                comp("ax", "axicon", (0, 0, 0), shape=circle(W * 0.45),
                     params={"alpha": alpha, "index": n_ax},
                     label="轴锥镜 α=%.1f mrad" % (alpha * 1e3)),
                _det_plane(z=min(0.15, 0.9 * zmax)),
            ],
            sources=[gauss_source(w0, wl=WL["633"], label="高斯 w₀=1.5 mm")],
            category="传播/光束", difficulty="进阶",
            description="轴锥镜（底面夹角 α=%.1f mrad、n=1.5）把高斯光束转成无衍射（贝塞尔）光束；"
                        "理论无衍射区长度 z_max=w₀/tan((n−1)α)≈%.0f mm。"
                        % (alpha * 1e3, zmax * 1e3),
            expected="中心亮斑直径 ≈2×2.405/(k·sin((n−1)α))=%.1f µm（比同尺寸高斯聚焦更小），"
                     "外侧为同心环；在 z_max 之前横向轮廓**不随距离变化**（无衍射），之后退化为普通高斯。"
                     % (core * 2 * 1e6),
            tags=["贝塞尔光束", "轴锥镜", "无衍射光束"],
            theory="中心斑 ~2.405/(k sinβ)，β=(n−1)α；无衍射区 z_max=w₀/tanβ", sample=(i < 2),
        )
        reg.add("beam_bessel_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))

    # 贝塞尔 vs 高斯：同一距离两个面
    for i, z in enumerate([0.05, 0.15, 0.3]):
        N, W, w0 = 1024, W8, 1.5e-3
        dx = W / N
        zmax = w0 / math.tan(0.5 * 0.03)
        cfg = build_config(
            name="贝塞尔 vs 高斯：z=%.0f mm 处的无衍射对比" % (z * 1e3),
            size=N, width=W, wl=WL["633"],
            components=[
                comp("ax", "axicon", (0, 0, 0), shape=circle(W * 0.45),
                     params={"alpha": 0.03, "index": 1.5}, label="轴锥镜 α=30 mrad"),
                det_circle("d1", z, W * 0.45, label="贝塞尔光束（轴锥镜后 %.0f mm）" % (z * 1e3)),
            ],
            sources=[gauss_source(w0, wl=WL["633"], label="高斯 w₀=1.5 mm")],
            category="传播/光束", difficulty="进阶",
            description="轴锥镜产生的贝塞尔光束在 z=%.0f mm 处的横向轮廓（无衍射区约 %.0f mm）。" % (z * 1e3, zmax * 1e3),
            expected="轮廓在 z_max=%.0f mm 内几乎与 z 无关（中心斑 ~%.1f µm 恒定）；"
                     "对照同束腰高斯光束在同一距离已展宽 %.1f 倍——这展示了贝塞尔光束「无衍射」的价值"
                     "（激光加工、光学镊子）。" % (zmax * 1e3, bessel_core_radius(WL["633"], 0.015) * 1e6,
                                                 waist_at(WL["633"], w0, z) / w0),
            tags=["无衍射光束", "贝塞尔", "对比"],
            theory="贝塞尔光束由锥面波矢叠加而成，横向轮廓在 z_max 内自再现",
        )
        reg.add("beam_bessel_vs_gauss_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"))


# ==========================================================================
# 球面波
# ==========================================================================

def spherical_waves(reg):
    for i, (R, conv) in enumerate([(0.1, False), (0.2, False), (0.5, False),
                                   (0.5, True), (1.0, True), (0.2, True)]):
        N, W = 512, W10
        dx = W / N
        scale = (WL["633"] * R / dx) ** 0.5
        check_sampling(dx, {"菲涅耳带尺度": scale}, 2.0)
        cfg = build_config(
            name="球面波（%s）：曲率半径 R=%.2f m" % ("会聚" if conv else "发散", R),
            size=N, width=W, wl=WL["633"],
            components=[_det_plane(z=0.2)],
            sources=[spherical_source(R, converging=conv, wl=WL["633"],
                                      label="球面波 R=%.2f m" % R)],
            category="传播/光束", difficulty="入门",
            description="%s球面波，等价于点光源位于光源面后方（或前方）%.2f m 处。" % ("会聚" if conv else "发散", R),
            expected=("强度沿半径按 1/r² 衰减（越靠边缘越暗），等相位面为同心球面（phase_u 上是同心环）；"
                      "发散波的高斯相位导数为正、曲率半径向外扩大。" if not conv else
                      "会聚球面波，波前向轴心收缩；若继续传播会在 R 处（轴上）形成焦点（本面取在 z=0.2 m，"
                      "焦点在光源前方 %.2f m 处）。" % R),
            tags=["球面波", "点光源", "波前曲率"],
            theory="U ∝ exp(ikr)/r；会聚/发散由二次相位符号区分",
            sample=(i < 2),
        )
        reg.add("beam_spherical_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))


# ==========================================================================
# 倾斜平面波的横向位移：Δx = z·tanθ
# ==========================================================================

def tilted_plane(reg):
    for i, tilt in enumerate([0.5e-3, 1e-3, 2e-3, 5e-3, 1e-2]):
        z, N, W = 0.2, 512, W10
        dx = W / N
        shift = tilt_shift(z, tilt)
        window_ok(W, {"横向位移": shift})
        cfg = build_config(
            name="倾斜平面波的横移：θ=%.1f mrad → Δx=%.2f mm" % (tilt * 1e3, shift * 1e3),
            size=N, width=W, wl=WL["633"],
            components=[_det_plane(z=z)],
            sources=[plane_source(wl=WL["633"], tilt_x=tilt,
                                  label="平面波（倾角 %.1f mrad）" % (tilt * 1e3))],
            category="传播/光束", difficulty="入门",
            description="平面波以 %.1f mrad 的倾角入射（绕 y 轴；内核用**精确方向余弦**而非傍轴近似），"
                        "在 z=%.0f mm 处读数。" % (tilt * 1e3, z * 1e3),
            expected="光斑整体横移 Δx=z·tanθ=%.3f mm（= %.1f 个像素）；强度仍是均匀的（平面波在自由空间中"
                     "形状不变）。相位图沿 x 是均匀的斜坡（phase_x 呈条纹），倾斜量越大条纹越密——"
                     "这正是干涉仪测角/测波前的原理。" % (shift * 1e3, shift / dx),
            tags=["倾斜入射", "方向余弦", "相位斜坡"],
            theory="Δx = z·tanθ，用精确方向余弦而非傍轴近似", sample=(i < 2),
        )
        reg.add("beam_tilt_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))


# ==========================================================================
# 透镜对高斯光束的变换（束腰变换公式）
# ==========================================================================

def gauss_lens(reg):
    for i, (w0, f, z_lens) in enumerate([(1.0e-3, 0.3, 0.0), (1.5e-3, 0.3, 0.0),
                                         (1.0e-3, 0.5, 0.0), (2.0e-3, 0.2, 0.0),
                                         (1.0e-3, 0.3, 0.15), (0.5e-3, 0.25, 0.0)]):
        W = W8
        w_out = WL["633"] * f / (math.pi * w0)
        N = pick_grid(w_out, W, minimum=4.0)
        dx = W / N
        check_sampling(dx, {"焦斑半径": w_out}, 4.0)
        check_sampling(dx, {"入射束腰": w0}, 8.0)
        zr = rayleigh_range(WL["633"], w0)
        cfg = build_config(
            name="高斯光束经薄透镜的束腰变换：w₀=%.1f mm, f=%.0f mm → w₁=%.1f µm"
                 % (w0 * 1e3, f * 1e3, w_out * 1e6),
            size=N, width=W, wl=WL["633"],
            components=[
                lens("l", (0, 0, 0), f=f, aperture=W * 0.45, label="薄透镜 f=%.0f mm" % (f * 1e3)),
                detector("det", z_lens + f, label="焦点附近", params={"strehl_aperture": 2.5 * w0,
                                                                     "strehl_distance": f}),
            ],
            sources=[gauss_source(w0, wl=WL["633"], label="高斯 w₀=%.1f mm" % (w0 * 1e3))],
            category="传播/光束", difficulty="进阶",
            description="束腰 %.1f mm 的高斯光束经焦距 %.0f mm 的薄透镜聚焦%s。"
                        % (w0 * 1e3, f * 1e3, "（束腰在透镜上）" if z_lens == 0 else
                           "（束腰距透镜 %.0f mm）" % (z_lens * 1e3)),
            expected="新束腰 w₁≈λf/(πw₀)=%.1f µm，位于透镜后 ≈%.0f mm；入射束腰越大焦斑越小"
                     "（反过来成像时 w₁·w₀ = λf/π 恒定）。入射瑞利长度 z_R=%.2f m，"
                     "透镜后焦深 ≈2z_R₁=2πw₁²/λ=%.2f mm。"
                     % (w_out * 1e6, f * 1e3, zr, 2 * math.pi * w_out ** 2 / WL["633"] * 1e3),
            tags=["束腰变换", "聚焦", "高斯光学"],
            theory="w₁ = λf/(πw₀)（束腰在透镜上时），共轭关系 w₁w₀ = λf/π", sample=(i < 2),
        )
        reg.add("beam_gauss_lens_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))


# ==========================================================================
# 涡旋光束的传播
# ==========================================================================

def vortex_prop(reg):
    for i, l in enumerate([1, 3, 5, 8, 10]):
        w0, N, W = 1.5e-3, 512, W8
        dx = W / N
        check_sampling(dx, {"束腰": w0}, 6.0)
        cfg = build_config(
            name="涡旋光束传播：拓扑荷 l=%d 的环随距离长大" % l,
            size=N, width=W, wl=WL["633"],
            components=[
                comp("spp", "spiral_phase", (0, 0, 0), shape=circle(W * 0.45),
                     params={"charge": l}, label="螺旋相位板 l=%d" % l),
                detector("d1", 0.1, label="z=100 mm"),
                detector("d2", 0.3, label="z=300 mm", params=None),
            ],
            sources=[gauss_source(w0, wl=WL["633"], label="高斯 w₀=1.5 mm")],
            category="传播/光束", difficulty="进阶",
            description="拓扑荷 l=%d 的光学涡旋在两个距离上读数：环随传播自然放大（衍射）。" % l,
            expected="两个面上的环半径比 ≈ 距离比（0.3/0.1=3，若已在远场区）；环宽 ∝ 距离；"
                     "轴上始终为零（相位奇点不会消失，拓扑荷是守恒量）。l=%d 的环比 l=1 粗、半径更大。"
                     % l,
            tags=["光学涡旋", "OAM", "传播"],
            theory="涡旋光束的环半径按远场衍射规律增长，拓扑荷守恒", sample=(i == 0),
        )
        reg.add("beam_vortex_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))

    # 轴锥镜叠加螺旋相位：高阶贝塞尔
    for i, l in enumerate([1, 2, 3]):
        w0, N, W = 1.5e-3, 1024, W8
        cfg = build_config(
            name="高阶贝塞尔光束：螺旋相位 + 轴锥镜（l=%d）" % l,
            size=N, width=W, wl=WL["633"],
            components=[
                comp("spp", "spiral_phase", (0, 0, 0), shape=circle(W * 0.45),
                     params={"charge": l}, label="螺旋相位板 l=%d" % l),
                comp("ax", "axicon", (0, 0, 0.01), shape=circle(W * 0.45),
                     params={"alpha": 0.03, "index": 1.5}, label="轴锥镜 α=30 mrad"),
                detector("det", 0.15, label="贝塞尔区"),
            ],
            sources=[gauss_source(w0, wl=WL["633"], label="高斯 w₀=1.5 mm")],
            category="传播/光束", difficulty="高级",
            description="螺旋相位板（l=%d）与轴锥镜串联：产生**高阶贝塞尔光束**（涡旋 + 无衍射）。" % l,
            expected="中心是**暗核**（拓扑荷 l=%d），外围是无衍射的同心环；"
                     "高阶贝塞尔光束的中心暗核尺寸比 l=1 更大，环结构在无衍射区内不变形。"
                     "可用于「无衍射光镊」与超分辨光片照明。" % l,
            tags=["高阶贝塞尔", "涡旋", "无衍射", "光镊"],
            theory="exp(ilθ)·J_l(k_r r)：径向为 l 阶贝塞尔函数",
        )
        reg.add("beam_high_order_bessel_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"))


# ==========================================================================
# 波长对光束的影响
# ==========================================================================

def wavelength_beam(reg):
    for i, key in enumerate(["405", "532", "633", "780", "1064"]):
        w0, z, N, W = 0.3e-3, 0.5, 512, W8
        wl = WL[key]
        dx = W / N
        wz = waist_at(wl, w0, z)
        cfg = build_config(
            name="波长对高斯光束发散的影响：λ=%s nm" % key,
            size=N, width=W, wl=wl,
            components=[_det_plane(z=z)],
            sources=[gauss_source(w0, wl=wl, label="高斯 w₀=0.30 mm @ %s nm" % key)],
            category="传播/光束", difficulty="入门",
            description="同一束腰 0.30 mm、同一距离 %.1f m，只改波长。" % z,
            expected="z 处半径 w(z)=%.3f mm（z_R=πw₀²/λ=%.2f m 随 λ 反比变化）；"
                     "波长越长瑞利长度越短、同距离展宽越大（θ=λ/(πw₀)=%.2f mrad）。"
                     % (wz * 1e3, rayleigh_range(wl, w0), gauss_divergence(wl, w0) * 1e3),
            tags=["波长", "发散角", "瑞利长度"],
            theory="z_R ∝ 1/λ，θ ∝ λ", sample=(i == 0),
        )
        reg.add("beam_wavelength_%02d.json" % (i + 1), cfg["meta"]["zod"]["category"],
                cfg["meta"]["zod"]["difficulty"], cfg["meta"]["zod"]["description"],
                cfg["meta"]["zod"]["expected"], cfg, tags=cfg["meta"]["zod"]["tags"],
                theory=cfg["meta"]["zod"].get("theory"), sample=cfg["meta"].get("sample", False))
