# Tetsuhiro-OZ (Optics-Zero)

**一台完整的光学模拟器**：完整几何光学（geo / lens / coating）与 Tetsuhiro-WOS
全部波动/量子光学能力（`optics/` 内核）在同一仓库、同一二进制中交付。

当前版本：**0.1.0-Preview**（2026-10-15）。

## 形态

- **单仓单二进制**：`oz`（`cmd/oz`）。`go.mod` 零第三方依赖（纯 Go 标准库）。
- 波动/量子内核已从 Tetsuhiro-WOS **内化**到本仓 `optics/`；WOS 仓库冻结为存档
  （tag `archive-pre-oz-absorb`）。
- 单一 HTTP 服务：`oz serve` 同时提供几何/膜系端点与全部波动/量子端点，
  以及内嵌 GUI（含几何/膜系模式页）。

## 模块划分

    geo/       几何光学核心：面型（平面/球面/圆锥/偶次与多项式非球面/Forbes Q 非球面/
               自由曲面 XY 多项式/泽尼克面）、矢量折射/反射追迹、近轴光学、光线扇/点列图、
               波前/OPL、赛德尔像差、色差、采样与渐晕、坐标断点（偏心/倾斜）、偏振追迹、
               几何 MTF/场曲、杂散光鬼像、热光学
    lens/      镜头设计：序列处方模型（JSON）、DLS 优化、消色差自主设计、
               Monte-Carlo 公差、灵敏度/逆灵敏度/补偿器、镀膜-镜头联合评价、热漂移
    coating/   镀膜设计：TMM 多层膜光谱（R/T/A、角度、偏振、相位）、材料色散库、
               膜系优化（DLS/层增广）、AR/高反/带通等模板
    wave/      几何→衍射桥：OPD→PSF/Strehl/MTF（FFT 由 optics 内核提供）
    optics/    波动/量子内核（内化自 WOS）：自研并行 FFT（含 Bluestein）、ASM 族传播、
               菲涅耳/夫琅禾费、非均匀/复折射率介质、部分相干、各向异性（Berreman）、
               28 类参数化元件与脚本元件、量子光学（Fock 基）、场景光路、MUSA/CUDA GPU
    doe/       衍射光学元件（相位型 DOE 设计、量化效率）
    xcheck/    几何↔波动互验（E1–E6 与 G7 联合示例，自动断言）
    cmd/oz/    CLI 与单一 HTTP 服务（API + GUI）
    library/   wave/（760 例）+ geo/（624 例）示例库与索引
    docs/      物理模型、内核、API、GUI、验证文档

## 快速开始

    go build ./cmd/oz && ./oz help         # CLI: design / eval / tol / coating / lens / geo / doe / serve
    ./oz design achromat --efl 100 --out lens.json
    ./oz coating design --template broadband-ar --out ar.json
    ./oz lens thermal lens.json --delta-t 20           # 热漂移（EFL/BFL 随温度）
    ./oz lens sens lens.json                           # 灵敏度表（--inverse/--compensate）
    ./oz geo mtf lens.json --fields 0,0.5              # 几何 MTF 与场曲
    ./oz geo ghost lens.json                           # 杂散光/鬼像路径
    ./oz doe focus --focal-mm 10 --levels 4            # DOE 一级效率与焦点偏差
    ./oz serve --addr :1130                # 本地 HTTP 服务（JSON API + GUI）
    go test ./...                          # 全量精度测试（全部确定性）

**GUI 里的 AI 助手**：顶栏「✦ AI 助手」打开右侧栏，填上自己的 OpenAI 兼容端点
（含本地 llama.cpp / Ollama）与 API Key 即可。助手会自动读取本二进制的元件规范
与当前光路状态，能解释元件、检查配置、设计光路（`oz-config` 一键应用）或写新的
脚本元件定义（`oz-element` 一键落盘、即刻可用）。密钥只存在浏览器里，服务端不记录。

GPU（MUSA，本机 MTT S4000）：`go build -tags musa -o oz-musa ./cmd/oz`，
`WOS_GPU=1 WOS_MUSA_HELPER=./wos-musa-fft ./oz-musa serve --addr :6008 --gpu`。
helper 源码见 `tools/musa-fft/`。

## 状态（2026-10-15）

- 测试：全仓 `go test ./...` 全绿；GPU/CPU 一致性 < 1e-9（MUSA 实测）。
- 示例库：`library/wave` 760 例（整库校验通过、抽样实跑通过）、`library/geo` 624 例
  （`tools/libcheck` 逐例复核 624/624 通过、两跑 sha256 逐位一致；含 singlet 140 /
  achromat 102 / 偏心倾斜 90 / 非球面 61 / Forbes Q 44 / 镜面 45；镜面家族为会聚设计、
  FOCUSZ 为实测最佳焦面）。
- 证据链：`docs/VERIFICATION.md`（TMM 对拍 1e-14、rayoptics 对拍、互验 E1–E6 等）。
- 图形与交互：内嵌 GUI（键盘+鼠标、移动端适配、几何/膜系模式页、高级分析折叠区、AI 助手右侧栏）。

详细运维（服务端口、外链、日志）：见 `/root/autodl-tmp/README-services.md`。
