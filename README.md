# Tetsuhiro-OZ (Optics-Zero)

基于 Tetsuhiro-WOS 的扩展光学模拟器：在波动光学内核之上，增加**几何光学仿真、
镜头（透镜）自主设计、镀膜自主设计**等模块，目标是一套精确、功能全面的光学设计工具。

当前版本：**0.1.0-Preview**（首版预览，2026-10-15）。

模块划分（开发中）：

    geo/       几何光学核心：面型（平面/球面/圆锥/非球面）、矢量折射/反射光线追迹、
               近轴光学（y-u 追迹、ABCD、焦距/主面/入瞳出瞳）、光线扇/点列图、
               波前/OPL、赛德尔像差、色差、采样与渐晕
    lens/      镜头设计：序列面型处方模型、一阶求解、像差分析、阻尼最小二乘优化、
               处方进出（JSON）、经典设计（双胶合消色差、三片式 Cooke）
    coating/   镀膜设计：特征矩阵法（TMM）多层膜光谱（R/T/A、角度、偏振）、
               材料色散库、膜系优化（DLS/层增广）、AR/高反/带通等设计目标
    cmd/oz/    CLI 与本地 Web 服务（HTTP API + 可视化）
    docs/      物理模型、内核开发、API 与使用文档

与 Tetsuhiro-WOS 的关系：本模块通过 `replace twos => ../Tetsuhiro-WOS` 直接复用
WOS 的 optics 内核（含 MUSA/cuFFT GPU 加速的 FFT），几何光学与波动光学互相验证。

## 快速开始

    go build ./cmd/oz && ./oz help     # CLI: design / eval / tol / coating / serve
    ./oz design achromat --efl 100 --out lens.json   # 自主设计消色差双胶合
    ./oz coating design --template broadband-ar --out ar.json
    ./oz serve --addr :1130            # 本地 HTTP 服务（JSON API + 单页前端）
    go test ./lens/ ./coating/ ./geo/  # 精度测试（全部确定性）

## 状态

- geo：追迹/近轴/采样/点列图/光线扇/OPD/赛德尔/色差/畸变 全部完成并验真（60 项测试）
- coating：TMM/材料/相位物理层已验真（与独立实现交叉对拍 1e-14）；优化器与模板完成；JSON I/O 完成
- lens：处方模型 + DLS 优化器 + 消色差自主设计（EFL 偏差 2.5µm、色差 11µm、点列图 4~7µm）+ Monte-Carlo 公差（确定性）已完成
- cmd/oz：CLI（design/eval/tol/coating/serve）+ HTTP 服务与单页前端可用
- wave：几何→衍射桥（OPD→PSF/Strehl/MTF），对照 Airy/解析 MTF/Maréchal 验证，真实双胶合 Strehl=0.985
- 内置示例库：波动光学（WOS）/几何光学（OZ）各 ≥520 例编写中
- docs：VERIFICATION.md（精确性证据）已建；其余待补
