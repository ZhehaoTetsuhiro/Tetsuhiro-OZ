"""量子光学示例（library/wave/quantum/*.json）—— 目标 40+。

注意：这些文件使用 **QuantumConfig** 协议（POST /api/quantum，同步返回 QuantumResult），
与波动光学的 Config 协议不同，因此放在 library/wave/ 的 quantum/ 子目录下，
校验脚本对它们调用 /api/quantum 而不是 /api/validate。
"""

import math

from wavegen.common import Registry

STATE_TYPES = ("vacuum", "fock", "coherent", "squeezed_vacuum", "two_mode_squeezed", "thermal")
GATE_TYPES = ("phase_shift", "beam_splitter", "displacement", "squeeze", "loss")


def qcfg(name, modes, cutoff, state, gates, description, expected, tags, theory,
         difficulty="进阶", postselect=None, sample=False):
    cfg = {
        "name": name,
        "modes": int(modes),
        "cutoff": int(cutoff),
        "state": state,
        "gates": gates,
    }
    if postselect is not None:
        cfg["postselect"] = postselect
    meta = {
        "name": name,
        "zod": {
            "category": "量子",
            "difficulty": difficulty,
            "description": description,
            "expected": expected,
            "tags": list(tags),
            "theory": theory,
            "protocol": "quantum",   # POST /api/quantum
        },
    }
    if sample:
        meta["sample"] = True
    cfg["meta"] = meta
    return cfg


def _add(reg, fname, cfg, reg_q):
    z = cfg["meta"]["zod"]
    reg.add(fname, z["category"], z["difficulty"], z["description"], z["expected"], cfg,
            tags=z["tags"], theory=z.get("theory"), sample=cfg["meta"].get("sample", False))


def _fs(occ):
    return {"type": "fock", "params": {"occupation": list(occ)}}


def _bs(m0, m1, R):
    return {"type": "beam_splitter", "params": {"mode0": m0, "mode1": m1, "reflectivity": R}}


def _ps(m, ph):
    return {"type": "phase_shift", "params": {"mode": m, "phase": ph}}


def _loss(m, T):
    return {"type": "loss", "params": {"mode": m, "transmittance": T}}


# ==========================================================================
# 1. Hong-Ou-Mandel 干涉
# ==========================================================================

def q_hom(reg):
    for i, R in enumerate([0.5, 0.3, 0.4, 0.6, 0.7]):
        cfg = qcfg("Hong-Ou-Mandel 干涉：分束比 R=%.1f" % R, 2, 4, _fs([1, 1]), [_bs(0, 1, R)],
                   "两个**全同光子**从不同端口同时入射 50:50 分束器（|1,1⟩ → BS）："
                   "不可分辨光子的量子干涉。",
                   ("R=0.5 时完美聚束：P(1,1)=0（两光子从不分开走），P(2,0)=P(0,2)=1/2。"
                    "这是 HOM 干涉的经典结果，也是量子光学里「光子不可分辨性」的最直接证据。"
                    if abs(R - 0.5) < 1e-9 else
                    "R=%.1f 时两光子不再全同地分配到两个输出端口，P(1,1)=2R(1−R)=%.3f 不再为零"
                    "（只有 50:50 分束器才有完美相消）。" % (R, 2 * R * (1 - R))),
                   ["HOM", "双光子干涉", "聚束", "全同性"],
                   "P(1,1) = |R − (1−R)|² /... 对 50:50 为 0；P(2,0)=P(0,2)=1/2", sample=(i == 0))
        _add(reg, "quantum/q_hom_%02d.json" % (i + 1), cfg, True)

    # HOM 的相位扫描（在一条臂上加相位，等效时间延迟）
    for i, ph in enumerate([0.0, math.pi / 4, math.pi / 2, 3 * math.pi / 4, math.pi]):
        cfg = qcfg("HOM 干涉臂相位扫描：φ=%.2f rad" % ph, 2, 4, _fs([1, 1]),
                   [_ps(0, ph), _bs(0, 1, 0.5)],
                   "在 HOM 干涉仪的入射臂上加 φ=%.2f rad 的相位（模拟两光子的延迟/相位不匹配）。"
                   % ph,
                   "相位差在**两光子之间**只影响不可分辨性；理想分束器上的 HOM 相消由分束器本身"
                   "（对称性）保证，因此即使加相位 P(0,2)=P(2,0)=1/2、P(1,1)=0 仍然成立——"
                   "这是 HOM 干涉与经典干涉最本质的区别（相消来自交换对称性，不来自相位）。",
                   ["HOM", "相位", "交换对称性"], "HOM 相消与单光子相位无关（对称 BS 约定）",
                   "高级", sample=(i == 0))
        _add(reg, "quantum/q_hom_phase_%02d.json" % (i + 1), cfg, True)


# ==========================================================================
# 2. 单光子马赫-曾德尔（单光子干涉）
# ==========================================================================

def q_single_photon_mz(reg):
    for i, ff in enumerate([0.0, 0.125, 0.25, 0.375, 0.5, 0.75, 1.0]):
        ph = ff * 2 * math.pi
        cfg = qcfg("单光子马赫-曾德尔：φ=%.2f×2π" % ff, 2, 4, _fs([1, 0]),
                   [_bs(0, 1, 0.5), _ps(1, ph), _bs(0, 1, 0.5)],
                   "单个光子进入马赫-曾德尔干涉仪（|1,0⟩ → BS → 相位 φ → BS）："
                   "每个端口探测到的是概率，但**概率本身**随相位呈现干涉。",
                   "P(1,0)=(1−cosφ)/2=%.3f、P(0,1)=(1+cosφ)/2=%.3f；"
                   "φ=0 时光子确定地回到输入端（相消）、φ=π 时确定地出现在另一个端口。"
                   "单光子的干涉是「一个光子与自己干涉」：实验中降低到单光子级仍能看到完整干涉条纹，"
                   "这排除了「两个光子相互干涉」的经典解释。" % ((1 - math.cos(ph)) / 2,
                                                                 (1 + math.cos(ph)) / 2),
                   ["单光子干涉", "马赫-曾德尔", "概率", "波粒二象性"],
                   "P(1,0)=(1−cosφ)/2，P(0,1)=(1+cosφ)/2", sample=(i < 2))
        _add(reg, "quantum/q_mz_%02d.json" % (i + 1), cfg, True)


# ==========================================================================
# 3. 相干态：泊松统计
# ==========================================================================

def q_coherent(reg):
    for i, a in enumerate([0.5, 1.0, 1.5, 2.0, 3.0]):
        cfg = qcfg("相干态 |α⟩：α=%.1f（泊松统计）" % a, 1, 25,
                   {"type": "coherent", "params": {"mode": 0, "alpha_re": a, "alpha_im": 0.0}},
                   [],
                   "单模相干态 |α⟩（α=%.1f）：激光的量子态，光子数服从泊松分布。" % a,
                   "⟨n⟩=|α|²=%.2f，方差 = ⟨n⟩（泊松），g²(0)=1（相干光）；"
                   "正交分量 var_x=var_p=1/4（散粒噪声基线，|α| 越大越接近经典极限）。"
                   "这就是「激光是最接近经典的光」的量子描述。" % a ** 2,
                   ["相干态", "泊松分布", "g2", "散粒噪声"],
                   "P(n)=e^{−|α|²}|α|^{2n}/n!；⟨n⟩=Var(n)=|α|²；Δx=1/2", sample=(i < 2))
        _add(reg, "quantum/q_coherent_%02d.json" % (i + 1), cfg, True)

    # 位移后的真空
    for i, a in enumerate([0.5, 1.0]):
        cfg = qcfg("位移真空态 D(α)|0⟩：α=%.1f" % a, 1, 15,
                   {"type": "vacuum", "params": {}},
                   [{"type": "displacement", "params": {"mode": 0, "alpha_re": a, "alpha_im": 0.0}}],
                   "对真空施加位移算符 D(α)=exp(αa†−α*a)：D(α)|0⟩ 就是相干态 |α⟩。",
                   "⟨n⟩=|α|²=%.2f、g²(0)=1、正交分量方差仍为 1/4 —— 位移只是把真空的"
                   "量子噪声球平移到相空间原点之外。" % a ** 2,
                   ["位移算符", "真空态", "相干态生成"],
                   "D(α)|0⟩ = |α⟩（相干态的标准生成方式）", sample=(i == 0))
        _add(reg, "quantum/q_displace_%02d.json" % (i + 1), cfg, True)


# ==========================================================================
# 4. 压缩态
# ==========================================================================

def q_squeezed(reg):
    for i, r in enumerate([0.25, 0.5, 0.75, 1.0, 1.5]):
        cfg = qcfg("压缩真空态：r=%.2f（%.1f dB）" % (r, 20 * r / math.log(10)),
                   1, 25, {"type": "squeezed_vacuum", "params": {"mode": 0, "r": r, "phase": 0.0}},
                   [],
                   "单模压缩真空 S(r)|0⟩：把一个正交分量的噪声压到散粒噪声以下，另一个分量相应放大。",
                   "var_x=e^{−2r}/4=%.4f（低于散粒噪声 1/4=0.25，即 %.1f dB 压缩）、"
                   "var_p=e^{2r}/4=%.4f；乘积 var_x·var_p=1/16（最小不确定态，纯态取等号）。"
                   "⟨n⟩=sinh²r=%.3f（压缩会产生光子对）。压缩光是 LIGO 提高灵敏度的手段。"
                   % (math.exp(-2 * r) / 4, -10 * math.log10(math.exp(-2 * r)) ,
                      math.exp(2 * r) / 4, math.sinh(r) ** 2),
                   ["压缩态", "正交分量", "海森堡极限", "LIGO"],
                   "var_x=e^{−2r}/4，var_p=e^{2r}/4，乘积 = 1/16", sample=(i < 2))
        _add(reg, "quantum/q_squeezed_%02d.json" % (i + 1), cfg, True)

    # 先压缩再位移：明亮压缩光
    for i, (r, a) in enumerate([(0.5, 1.0), (1.0, 2.0)]):
        cfg = qcfg("明亮压缩光：S(%.1f) 后 D(%.1f)" % (r, a), 1, 30,
                   {"type": "vacuum", "params": {}},
                   [{"type": "squeeze", "params": {"mode": 0, "r": r, "phase": 0.0}},
                    {"type": "displacement", "params": {"mode": 0, "alpha_re": a, "alpha_im": 0.0}}],
                   "先压缩真空再位移 → 「明亮压缩光」（有经典幅度但某分量噪声低于散粒噪声）。",
                   "⟨n⟩≈|α|²+sinh²r=%.2f，var_x 仍为 e^{−2r}/4=%.4f —— 位移不改变压缩量，"
                   "这就是实验上产生压缩光的标准方案（光学参量放大 + 强本振）。"
                   % (a ** 2 + math.sinh(r) ** 2, math.exp(-2 * r) / 4),
                   ["明亮压缩光", "压缩+位移", "连续变量"],
                   "D(α)S(r)|0⟩：位移不改变压缩方差", sample=(i == 0))
        _add(reg, "quantum/q_bright_squeezed_%02d.json" % (i + 1), cfg, True)


# ==========================================================================
# 5. 双模压缩（EPR 纠缠）
# ==========================================================================

def q_tms(reg):
    for i, r in enumerate([0.5, 1.0, 1.5]):
        cfg = qcfg("双模压缩真空（EPR）：r=%.1f" % r, 2, 20,
                   {"type": "two_mode_squeezed", "params": {"r": r}}, [],
                   "双模压缩真空 sech(r)Σtanh(r)ⁿ|n,n⟩：模式 0 与 1 的光子数**完全关联**（EPR 纠缠源）。",
                   "⟨n₀⟩=⟨n₁⟩=sinh²r=%.3f；P(n,m) 只在 n=m 上非零（完美关联）；"
                   "单模约化态是热态（g²=2）。这是连续变量量子信息（量子隐形传态、EPR 判据）"
                   "的标准资源。" % math.sinh(r) ** 2,
                   ["双模压缩", "EPR", "纠缠", "光子数关联"],
                   "|ψ⟩=sech r Σ tanhⁿr |n,n⟩", sample=(i < 2))
        _add(reg, "quantum/q_tms_%02d.json" % (i + 1), cfg, True)

    # 双模压缩 + 后选择：herald 单光子态
    for i, r in enumerate([0.5, 1.0, 1.5]):
        cfg = qcfg("双模压缩后选择（heralded 单光子）：r=%.1f" % r, 2, 20,
                   {"type": "two_mode_squeezed", "params": {"r": r}}, [],
                   "对双模压缩态的模式 1 做**后选择**（只保留恰好 1 个光子），模式 0 就被预报为单光子态。",
                   "后选择概率 p=Σ_n |c_n|²·n·... ≈ sech²r·tanh²r/(1−tanh⁴r)（该分支的概率），"
                   "条件态接近 |1⟩；r 越小纯度越高但概率越低。这就是「预报式单光子源」"
                   "（heralded single photon）的量子光学原型，也是 KLM 线性光学量子计算的基础资源。",
                   ["后选择", "预报单光子", "herald", "KLM"],
                   "Postselect 模式 1 = 1 光子 → 模式 0 的条件态 ≈ |1⟩",
                   postselect={"modes": [1], "counts": [1]}, sample=(i == 0))
        _add(reg, "quantum/q_tms_postselect_%02d.json" % (i + 1), cfg, True)


# ==========================================================================
# 6. 热态
# ==========================================================================

def q_thermal(reg):
    for i, n in enumerate([0.5, 1.0, 2.0, 3.0]):
        cfg = qcfg("热态：平均光子数 %.1f（g²=2）" % n, 1, 25,
                   {"type": "thermal", "params": {"mean_n": [n]}}, [],
                   "单模**热态**（平均光子数 %.1f）：混合态，走密度矩阵后端，光子数服从几何分布。" % n,
                   "P(n)=n̄ⁿ/(1+n̄)^{n+1}（几何/玻色-爱因斯坦分布）、g²(0)=2（聚束，"
                   "比相干光的 1 高出一倍）、Var(n)=n̄+n̄²>n̄（超泊松统计）。"
                   "热态是黑体辐射、白炽灯与散斑场的量子描述。" % (),
                   ["热态", "几何分布", "g2=2", "混合态", "密度矩阵"],
                   "P(n)=n̄ⁿ/(1+n̄)^{n+1}；g²(0)=2；Var(n)=n̄(n̄+1)", sample=(i < 2))
        _add(reg, "quantum/q_thermal_%02d.json" % (i + 1), cfg, True)

    # 热态 + 损耗：仍为热态（温度降低）
    for i, (n, T) in enumerate([(2.0, 0.5), (2.0, 0.25)]):
        cfg = qcfg("热态经过损耗：n̄=2.0、T=%.2f" % T, 1, 25,
                   {"type": "thermal", "params": {"mean_n": [n]}}, [_loss(0, T)],
                   "热态经过透射率 %.2f 的损耗信道：热态在分束器模型下**仍然是热态**"
                   "（只是平均光子数乘以 T）。" % T,
                   "⟨n⟩=T·n̄=%.2f，g²(0) 仍为 2（热态在损耗下保持热态的统计特征）；"
                   "迹守恒（Tr ρ=1）。这与 Fock 态完全不同（见 loss 系列），"
                   "是「经典噪声无法用衰减消除」的量子表述。" % (T * n),
                   ["热态", "损耗信道", "经典噪声"], "热态经损耗仍为热态：⟨n⟩→T⟨n⟩", sample=(i == 0))
        _add(reg, "quantum/q_thermal_loss_%02d.json" % (i + 1), cfg, True)


# ==========================================================================
# 7. 损耗信道
# ==========================================================================

def q_loss(reg):
    for i, T in enumerate([0.1, 0.3, 0.5, 0.7, 0.9]):
        cfg = qcfg("单光子损耗信道：T=%.1f" % T, 1, 5, _fs([1]), [_loss(0, T)],
                   "单光子态 |1⟩ 经过透射率 %.1f 的损耗信道（振幅阻尼 Kraus 算符）。" % T,
                   "|1⟩ → 二项分布：P(0)=1−T=%.2f、P(1)=T=%.2f，⟨n⟩=T=%.2f。"
                   "损耗把纯态变成混合态（需要密度矩阵描述）——光子「丢失」的概率就是 1−T。"
                   % (1 - T, T, T),
                   ["损耗信道", "振幅阻尼", "二项分布", "混合态"],
                   "E_l=√((1−T)^l/l!)T^{n/2}a^l；|1⟩→T|1⟩⟨1|+(1−T)|0⟩⟨0|", sample=(i < 2))
        _add(reg, "quantum/q_loss_%02d.json" % (i + 1), cfg, True)

    # 相干态经损耗：仍相干
    for i, T in enumerate([0.5, 0.25]):
        cfg = qcfg("相干态经损耗：α=2、T=%.2f" % T, 1, 25,
                   {"type": "coherent", "params": {"mode": 0, "alpha_re": 2.0, "alpha_im": 0.0}},
                   [_loss(0, T)],
                   "相干态 |α=2⟩ 经过损耗：经典光衰减时**保持相干**（只是幅度变小）。",
                   "⟨n⟩=T|α|²=%.2f、g²(0)=1（仍是泊松/相干光）。与热态、Fock 态对比："
                   "损耗对不同量子态的「伤害」完全不同——这是量子光学中「经典态最鲁棒」的原因。"
                   % (T * 4),
                   ["相干态", "损耗", "鲁棒性"], "相干态经损耗仍是相干态 α→√T·α", sample=(i == 0))
        _add(reg, "quantum/q_coherent_loss_%02d.json" % (i + 1), cfg, True)


# ==========================================================================
# 8. 玻色采样与多模干涉
# ==========================================================================

def q_boson_sampling(reg):
    specs = [("3 模 3 光子均分网络", 3, 4, [1, 1, 1], [_bs(0, 1, 0.5), _bs(1, 2, 0.5)]),
             ("3 模 2 光子（HOM 扩展）", 3, 4, [1, 1, 0], [_bs(0, 1, 0.5), _bs(1, 2, 0.5),
                                                          _ps(2, 1.0)]),
             ("4 模 2 光子（T 型网络）", 4, 3, [1, 1, 0, 0],
              [_bs(0, 1, 0.5), _bs(1, 2, 0.5), _bs(2, 3, 0.5), _ps(1, 1.5707963267948966)]),
             ("4 模 4 光子（玻色采样雏形）", 4, 3, [1, 1, 1, 1],
              [_bs(0, 1, 0.5), _bs(1, 2, 0.5), _bs(2, 3, 0.5), _bs(0, 3, 0.5), _ps(1, 0.7853981633974483)]),
             ("5 模 3 光子干涉", 5, 3, [1, 1, 1, 0, 0],
              [_bs(0, 1, 0.5), _bs(1, 2, 0.5), _bs(2, 3, 0.5), _bs(3, 4, 0.5)])]
    for i, (lbl, M, C, occ, gates) in enumerate(specs):
        cfg = qcfg("线性光学网络：%s" % lbl, M, C, _fs(occ), gates,
                   "%s：%d 个不可分辨光子进入 %d 模线性光学网络（分束器 + 相位）。" % (lbl, sum(occ), M),
                   "输出必须从**全模式联合分布 joint_full**（长度 (cutoff+1)^M，little-endian 下标）读取："
                   "两两联合分布只是它的边缘。玻色子的聚束（boson bunching）使输出分布与"
                   "「把光子当成经典小球独立随机走」的预言显著不同——这正是玻色采样的计算能力来源"
                   "（也是「量子优越性」实验的核心）。",
                   ["玻色采样", "多模干涉", "联合分布", "聚束"],
                   "输出概率 = |Perm(U_S)|²（子矩阵的积和式）", sample=(i < 2))
        _add(reg, "quantum/q_boson_%02d.json" % (i + 1), cfg, True)


# ==========================================================================
# 9. KLM / 后选择门
# ==========================================================================

def q_klm(reg):
    specs = [("后选择「模式 1 为空」（herald 单光子）", 2, 4, [1, 1], [_bs(0, 1, 0.5)],
              {"modes": [1], "counts": [0]}),
             ("后选择「模式 1 有 1 光子」", 2, 4, [1, 1], [_bs(0, 1, 0.5)],
              {"modes": [1], "counts": [1]}),
             ("后选择两模都为 0（真空投影）", 2, 4, [1, 1], [_bs(0, 1, 0.5)],
              {"modes": [0, 1], "counts": [0, 0]}),
             ("非线性符号门雏形：后选择 + 相位", 3, 3, [1, 1, 0],
              [_bs(0, 1, 0.5), _ps(1, 1.5707963267948966), _bs(0, 2, 0.5)],
              {"modes": [2], "counts": [0]})]
    for i, (lbl, M, C, occ, gates, post) in enumerate(specs):
        cfg = qcfg("KLM 后选择门：%s" % lbl, M, C, _fs(occ), gates,
                   "线性光学量子计算的核心工具——**测量 + 后选择**：%s。"
                   "线性光学本身没有非线性，但「测量到特定光子数」的事件会在剩余模式上"
                   "诱导出非线性演化（KLM 方案）。" % lbl,
                   "返回结果里除常规统计外还有 postselect_probability（该分支的概率）；"
                   "剩余态已被重新归一化。本内核用 Postselect{modes, counts} 接口实现投影，"
                   "这是 heralded 门、纠缠制备与量子纠错后选择的标准形式。",
                   ["KLM", "后选择", "heralded", "投影测量"],
                   "投影到指定模式光子数子空间并归一化；给出分支概率",
                   difficulty="高级", postselect=post, sample=(i < 2))
        _add(reg, "quantum/q_klm_%02d.json" % (i + 1), cfg, True)

    # 单光子马赫-曾德尔 + 后选择（对比）
    cfg = qcfg("单光子 MZ + 后选择：φ=π/2 时的条件态", 3, 4, _fs([1, 0, 0]),
               [_bs(0, 1, 0.5), _ps(1, math.pi / 2), _bs(0, 2, 0.5)],
               "单光子 MZ 的输出再与第三个模式干涉，并后选择模式 2 为空。",
               "后选择把连续变量型的干涉转换为离散的条件态：只保留「模式 2 为真空」的分支后，"
               "模式 0、1 的条件分布与无条件分布在定义上完全不同（这是量子态制备/过滤的常用技巧）。",
               ["后选择", "态制备", "条件分布"], "Postselect 之后的态为条件态",
               difficulty="高级", postselect={"modes": [2], "counts": [0]})
    _add(reg, "quantum/q_mz_postselect.json", cfg, True)


# ==========================================================================
# 10. 联合分布 / 符合计数 / 参数扫掠
# ==========================================================================

def q_joint(reg):
    for i, ph in enumerate([0.0, math.pi / 2, math.pi, 3 * math.pi / 2]):
        cfg = qcfg("双光子符合计数 vs 相位：φ=%.2f rad" % ph, 2, 5, _fs([1, 1]),
                   [_bs(0, 1, 0.5), _ps(0, ph), _bs(0, 1, 0.5)],
                   "两次 50:50 分束器构成的双光子干涉仪（中间加相位 φ=%.2f）——"
                   "两级 HOM 的复合结构，输出符合计数对相位敏感。" % ph,
                   "两级对称分束器的组合使两光子重新合束：输出 P(1,1) 随相位在 0 与 1 之间振荡"
                   "（第二级分束器把已经聚束的光子对重新分开或者再次聚束）。"
                   "从 joint_distributions['0,1'] 直接读 P(n₀,n₁)。",
                   ["符合计数", "联合分布", "双光子干涉仪"], "joint 分布随相位振荡", sample=(i == 0))
        _add(reg, "quantum/q_joint_%02d.json" % (i + 1), cfg, True)

    # 截断扫描：验证数值收敛
    for i, C in enumerate([2, 4, 8, 16]):
        cfg = qcfg("截断扫描：cutoff=%d（压缩态 r=1.5 的收敛）" % C, 1, C,
                   {"type": "squeezed_vacuum", "params": {"mode": 0, "r": 1.5, "phase": 0.0}}, [],
                   "压缩真空态是 Fock 空间的无穷叠加，数值上必须截断到 cutoff=%d。" % C,
                   "压缩态的系数按 tanhⁿr 衰减（tanh1.5=0.905），需要 cutoff≳%d 才能把范数收敛到 "
                   "1e-3 以内；范数 norm 会随着 cutoff 增大趋近 1（这是判断截断是否足够的标准做法）。"
                   "⟨n⟩=sinh²1.5=%.3f 应在 cutoff≥%d 时稳定。"
                   % (int(abs(math.log(1e-4) / math.log(1 / math.tanh(1.5)))), math.sinh(1.5) ** 2,
                      int(math.sinh(1.5) ** 2 * 3)),
                   ["截断", "收敛性", "数值误差"], "平面波截断误差 ~tanh^{cutoff}(r)",
                   difficulty="高级", sample=(i < 2))
        _add(reg, "quantum/q_cutoff_%02d.json" % (i + 1), cfg, True)

    # 真空态（基线）
    cfg = qcfg("真空态基线：⟨n⟩=0、var=1/4", 1, 4, {"type": "vacuum", "params": {}}, [],
               "全真空态 |0⟩：量子光学的基线（没有任何光子，但存在零点涨落）。",
               "⟨n⟩=0，光子数分布 P(0)=1；正交分量 var_x=var_p=1/4（零点/散粒噪声基线，"
               "不是零！）。真空涨落是压缩、Casimir 效应与自发辐射的根源。",
               ["真空态", "零点涨落", "散粒噪声"], "var(x_θ)=1/4", difficulty="入门", sample=True)
    _add(reg, "quantum/q_vacuum.json", cfg, True)


def quantum_all(reg):
    q_hom(reg)
    q_single_photon_mz(reg)
    q_coherent(reg)
    q_squeezed(reg)
    q_tms(reg)
    q_thermal(reg)
    q_loss(reg)
    q_boson_sampling(reg)
    q_klm(reg)
    q_joint(reg)
