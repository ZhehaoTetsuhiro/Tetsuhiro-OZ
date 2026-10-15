# HTTP API（统一服务 `oz serve`）

Tetsuhiro-OZ 只有一个 HTTP 服务：`oz serve`。它在同一个 mux 上提供几何/膜系 API（`/api/version`、`/api/lens/*`、`/api/coating/*`）、波动/量子 API（`/api/catalog`、`/api/simulate`、`/api/runs/*` 等）和 GUI（`/`）。本文件是该服务的接口参考；物理约定见 `docs/kernel/PHYSICS.md`，内核结构见 `docs/kernel/KERNEL.md`，GPU 见 `docs/kernel/GPU.md`，GUI 见 `docs/GUI.md`。

基址：`http://localhost:1130`（默认端口，见下）。全部 JSON 响应带 `Access-Control-Allow-Origin: *`；二进制与 PNG/SVG 端点除外。

## 1. 启动与参数

    oz serve [选项]

| 选项 | 默认 | 说明 |
|---|---|---|
| `-addr` | `:1130` | 监听地址。旧 WOS 端口用 `-addr :1120` |
| `-max-run-mb` | `512` | 已完成仿真结果在内存中的保留上限（MB，LRU 驱逐）。WOS 原参数名，单位 MB |
| `-elements` | 空 | 额外的元件定义目录（逗号分隔，优先级高于默认目录） |
| `-no-elements-watch` | 关 | 关闭定义文件自动重载（默认每 2 秒检查一次） |
| `-gpu` | 关 | 启用 GPU 后端（需以 `-tags cuda` 或 `-tags musa` 构建；等价于 `WOS_GPU=1`） |
| `-check-elements` | 关 | 校验元件定义文件，打印结果后退出（有错误则退出码 1） |
| `-gen-go <文件>` | 空 | 把元件定义导出为 Go 原生元件源码（stdout），`-gen-name` 指定类型名 |

环境变量（由 optics 内核读取，无需参数）：

| 变量 | 作用 |
|---|---|
| `WOS_GPU=1` | 启动时启用 GPU 后端（不可用时降级为 CPU，并在 `/api/health` 中可见） |
| `WOS_MUSA_HELPER` | MUSA helper 可执行文件路径（默认在可执行文件同目录或 PATH 中查找 `wos-musa-fft`） |
| `WOS_MUSA_SOCKET` | MUSA helper 的 UNIX socket 路径（默认 `/tmp/wos-musa-fft.sock`） |

从旧 `wos` 迁移：`wos -addr :1120 -gpu` → `oz serve -addr :1120 -gpu`；环境变量与元件定义目录（`elements/`、`~/.wos/elements`）不变。

## 2. 路径归属（单一 mux）

| 路径 | 归属 | 说明 |
|---|---|---|
| `/api/version`、`/api/lens/design`、`/api/lens/eval`、`/api/lens/tol`、`/api/lens/thermal`、`/api/lens/sens`、`/api/lens/parax`、`/api/lens/group/layout`、`/api/lens/group/image`、`/api/lens/group/afocal`、`/api/lens/group/dispersion`、`/api/lens/group/set`、`/api/lens/group/prescription`、`/api/geo/mtf`、`/api/geo/ghost`、`/api/doe/focus`、`/api/coating/design`、`/api/coating/spectrum` | OZ | 精确路径，优先于 `/api/` 前缀 |
| `/api/...`（其余） | WOS | 仿真 API，由 `server.Server.Handler()` 处理；未知路径 404 |
| `/oz/` | 重定向 | 301 → `/#oz`（合并前的几何/膜系页面地址） |
| `/` 及静态文件 | GUI | 内嵌在二进制中：`index.html`、`app.js`、`oz.js`、`style.css`、`favicon.svg`、字体 |

## 3. 约定

- 请求与响应：JSON（`Content-Type: application/json; charset=utf-8`），除非端点另有说明。
- 错误：`{"error": "…"}`，状态码 400（请求无效）、404（资源不存在）、405（方法不允许）、409（运行未完成）、429（排队已满）、500（内部错误或响应含 NaN/Inf）。
- 方法：计算端点为 POST，查询端点为 GET；错误方法返回 405。`GET /api/health` 不限制方法（既有行为）。
- 工作量上限（超出即 400，计算前拒绝）：光谱点数 ≤ 10000；公差运行次数 ≤ 20000；设计迭代 ≤ 10000；热漂移 |ΔT| ≤ 500 K；灵敏度参数 ≤ 128 个、补偿器迭代 ≤ 10000；MTF 视场 ≤ 64、每视场光线 ≤ 200000、LSF 分箱 ≤ 4096、频率点 ≤ 512；鬼像双反射路径 ≤ 512 条、采样环数 ≤ 8；DOE 网格 ≤ 2048（2 的幂）、级数 ≤ 4096；透镜组元件数 ≤ 512、Hexapolar 采样环数 ≤ 8、焦面搜索粗采样 ≤ 4096、单次编辑 ≤ 512 项。仿真网格与元件数由 `/api/validate` 检查。
- 请求体上限：`/api/quantum` 1 MB；`/api/validate`、`/api/simulate`、`/api/convert` 2 MB；`/api/lens/design`、`/api/lens/eval`、`/api/lens/thermal`、`/api/lens/sens`、`/api/geo/mtf`、`/api/geo/ghost`、`/api/doe/focus`、`/api/coating/design` 1 MB；`/api/lens/tol`、`/api/lens/parax`、`/api/lens/group/*` 4 MB；`/api/coating/spectrum` 8 MB。
- 确定性：全部计算端点对相同请求给出逐字节相同的响应（无时钟、无 map 迭代顺序、无随机数）；`cmd/oz/analysis_test.go` 的 `TestAnalysisEndpointDeterminism` 与 `cmd/oz/group_api_test.go` 的 `TestGroupEndpointDeterminism` 逐个端点断言。
- `GET /api/version` 返回的 `version` 与仓库根目录 `VERSION` 文件一致（由 `cmd/oz` 的测试保证）。

## 4. 几何与膜系（OZ）

本节每个端点与对应的 CLI 子命令共用同一个实现（`cmd/oz/analysis.go`、`cmd/oz/group_api.go` 的 `run*` 函数，`cmd/oz/analysis_ghost.go` 为鬼像）：`oz lens thermal`、`oz lens sens`、`oz lens group ...`、`oz lens parax`、`oz geo mtf`、`oz geo ghost`、`oz doe focus` 与 HTTP 端点给出同样的数字（CLI 用 `--json` 可导出与本文件字段一致的响应）。

### GET /api/version

    {"name":"tetsuhiro-oz","version":"0.1.0-Preview","go":"go1.25.10","time":"2026-10-15T18:43:17+08:00"}

### POST /api/lens/design

消色差双胶合自主设计（与 `oz design achromat` 同一调用）。

请求（字段均可省略，省略时取默认）：

| 字段 | 类型 | 默认 | 说明 |
|---|---|---|---|
| `efl_mm` | number | 必填，>0 | 有效焦距（mm） |
| `crown` / `flint` | string | `N-BK7` / `F2` | 玻璃目录名 |
| `semi_aperture_mm` | number | `efl_mm/10` | 半口径（mm） |
| `field_deg` | number | 未指定时取 0.5° 的默认视场 | 离轴视场半角（度）；显式给出时生效 |
| `skip_polish` | bool | false | 只做解析解与厚透镜消色差，不做 DLS 抛光 |
| `max_iter` | int | 0（默认上限） | DLS 迭代上限，最大 10000（超出返回 400） |

响应：`{"prescription": {...}, "report": "…", "optimizer": {"iterations", "evals", "merit", "converged"}}`。`prescription` 可直接作为 `/api/lens/eval`、`/api/lens/tol` 的请求体（或其 `prescription` 字段）。

    curl -s -X POST localhost:1130/api/lens/design -d '{"efl_mm":100,"semi_aperture_mm":10}'

### POST /api/lens/eval

请求体为一个处方（`prescription` JSON）。响应：`{"report": "…", "metrics": {"efl_mm", "bfl_mm", "f_number", "color_shift_mm"}}`。`color_shift_mm` 为 F–C 轴向色差。

### POST /api/lens/tol

确定性 Monte-Carlo 公差分析（同一 `seed` 同一结果）。

    {
      "prescription": { ... },
      "spec":     {"radius_mm":0.05, "thickness_mm":0.05, "conic_k":0, "coef_abs":0, "index_delta":0.0005},
      "criteria": {"target_efl_mm":100, "max_focal_err_mm":0.1, "max_color_shift_mm":0,
                   "max_rms_on_axis_mm":0, "max_rms_off_axis_mm":0},
      "runs": 2000,
      "seed": 1
    }

`criteria` 至少给出一个判据时才计算良率（`yield`）。响应为 `ToleranceReport` JSON：`runs`、`seed`、`yield`（满足全部判据的比例）、`nominal`、`on_axis_rms`、`off_axis_rms`、`color_shift`、`focal_err`（分布统计）、`fail_reasons`、`failed_traces`。`runs` 缺省为 2000，最大 20000（超出返回 400），耗时随次数线性增长。

### POST /api/lens/thermal

一阶热漂移（G9）：把处方从 20 °C 变到 20 °C + ΔT，报告 EFL/BFL 的变化与每开尔文割线斜率（同一调用为 `oz lens thermal`）。温度模型：每个玻璃面按 CTE 缩放半径/非球面系数/半口径（胶合面取两侧玻璃 CTE 的均值，空气面取镜筒 CTE），玻璃间隙按玻璃 CTE、空气间隙按镜筒 CTE 缩放，玻璃折射率叠加温致 Δn（见 `lens/thermal_model.go` 与 `geo/thermo_model.go`）。

请求：

| 字段 | 类型 | 默认 | 说明 |
|---|---|---|---|
| `prescription` | object | 必填 | 处方 JSON（`/api/lens/design` 的 `prescription` 可直接复用） |
| `delta_t_k` | number | 必填，≠0 | 温度变化（K，相对 20 °C），\|ΔT\| ≤ 500（超出返回 400） |
| `lambda_nm` | number | `587.5618` | 波长 nm |
| `mount_cte_per_k` | number | `0`（刚性镜筒） | 隔圈/镜筒的线膨胀系数 1/K，\|α\| ≤ 0.01 |

响应：

| 字段 | 说明 |
|---|---|
| `delta_t_k`、`lambda_nm`、`mount_cte_per_k` | 回显请求参数 |
| `efl0_mm`、`bfl0_mm` | 20 °C 的 EFL/BFL（mm） |
| `efl1_mm`、`bfl1_mm` | 20 °C+ΔT 的 EFL/BFL（mm） |
| `efl_per_k_mm`、`bfl_per_k_mm` | 割线斜率 (X1−X0)/ΔT（mm/K） |

    curl -s -X POST localhost:1130/api/lens/thermal -d '{"prescription":{…},"delta_t_k":10,"mount_cte_per_k":0}'

    {"delta_t_k":10,"lambda_nm":587.5618,"mount_cte_per_k":0,
     "efl0_mm":60.000000,"efl1_mm":60.001695,"bfl0_mm":58.147130,"bfl1_mm":58.148709,
     "efl_per_k_mm":0.0001695347,"bfl_per_k_mm":0.00015790}

### POST /api/lens/sens

灵敏度表（∂f/∂x），可选逆灵敏度公差分配与补偿器验算（G6；同一调用为 `oz lens sens`）。参数集合按卡序固定生成：每面曲率 c = 1/R（平面也用曲率）、每面到下一面的厚度、非零圆锥常数、非零非球面系数、每个玻璃面的折射率 Δn；标签为 `s<卡>.<c|thk|K|A4|…|dn>`。

请求（除 `prescription` 外全部可选）：

| 字段 | 类型 | 默认 | 说明 |
|---|---|---|---|
| `prescription` | object | 必填 | 处方 JSON |
| `lambda_nm` | number | `587.5618` | 波长 nm |
| `metric` | string | `efl` | `efl`：∂EFL/∂x（mm/参数单位）；`merit`：固定像面点列 merit 的 ∂M/∂x（G6 口径） |
| `rel_step` | number | `1e-3` | 最粗相对步长 |
| `fixed_step` | bool | false | 只用固定步长（快，不给收敛判据） |
| `inverse` | bool | false | 逆灵敏度分配（自动要求 `metric=merit`） |
| `budget_mm` | number | — | 逆灵敏度的 merit 预算 ΔM（mm，>0，`inverse` 时必填） |
| `mode` | string | `rss` | `rss`（统计）或 `linear`（最坏情况） |
| `max_curvature_per_mm` | number | `2e-5` | 曲率公差上限（≈R=50 mm 处 0.05 mm 半径公差） |
| `max_thickness_mm` | number | `0.05` | 厚度公差上限 mm |
| `max_index` | number | `5e-4` | 折射率公差上限 Δn |
| `max_conic` | number | `0.01` | 圆锥常数上限 |
| `max_coeff` | number | `1e-7` | 非球面系数上限 |
| `compensate` | bool | false | 补偿器验算（自动要求 `metric=merit`） |
| `perturb_radius_mm` / `perturb_thickness_mm` / `perturb_index` | number | 0 | 扰动量 dR/dt/Δn（加到每个有限半径/间隔/玻璃上） |
| `compensator_vars` | string[] | 空 | 补偿器变量 `s<卡>.<R\|thk\|K\|A4…>`；留空时只用像面补偿 |
| `compensator_ranges_mm` | number[] | 空 | 各变量半行程 mm（可选，与变量一一对应） |
| `image_shift` | bool | false | 在给出变量时额外把像面加入补偿器 |
| `shift_range_mm` | number | `1` | 像面补偿半行程 mm |
| `max_iter` | number | `100` | 补偿器 DLS 迭代上限 |

响应：`metric`、`lambda_nm`、`nominal`（f 在名义处的值）、`evals`（求值次数）、`rows`（每行 `param`、`deriv`、`err`（不确定度）、`naive`（普通中心差分）、`rel_step`、`step`、`converged`）；`inverse` 时另有 `allocation`（`mode`、`budget_mm`、`predicted_mm`、`rel_err`、`rows`：`param`、`sens`、`tol`、`share`、`capped`）；`compensate` 时另有 `compensation`（`nominal_mm`、`compensated_nominal_mm`、`perturbed_mm`、`compensated_mm`、`recovery`（可为 `null`：扰动未使 merit 恶化）、`image_shift_mm`、`iterations`、`converged`）。

`inverse` 的分配使用 `lens.AllocateTolerances`：每个未封顶参数贡献相同份额（linear：|s·t| = ΔM/N；rss：(s·t)² = ΔM²/N），触到上限的参数固定在 `max_*` 并把剩余预算分给其余参数；预算无法放置（全部参数先到上限）返回 400，提示放宽 `max_*` 或降低 `budget_mm`。参数 > 128 返回 400。

    curl -s -X POST localhost:1130/api/lens/sens -d '{"prescription":{…},"metric":"merit","inverse":true,"budget_mm":0.005}'

### POST /api/lens/group/layout

透镜组（`lens.Group`，元件级可编辑模型：透镜/反射镜/光阑、前后表面面型、玻璃、厚度、可负的 gap、偏心与倾斜、fold 折转）的逐面全局布局。请求体是 group 的 JSON wire form：既可以直接平铺在顶层（`name`/`elements`/`stop`/`origin_z_mm`），也可以嵌在 `"group"` 字段里（两种写法等价，响应逐字节相同）。同一模型与 `oz lens group layout` 同源。

group 的 wire form（`lens.Group` 的 JSON）：

| 字段 | 类型 | 默认 | 说明 |
|---|---|---|---|
| `name` | string | `""` | 组名 |
| `object_z_mm` | number 或 `"inf"` | `"inf"` | 物面轴向位置（mm）；`"inf"` = 无穷远物 |
| `elements` | array | 必填 | 元件表（见下） |
| `stop` | int | `0` | 光阑元件下标（该元件的前表面是光阑面） |
| `origin_z_mm` | number | `0` | 第一面顶点的全局 z（mm） |

元件（`elements[]`）：`kind`（`"lens"`/`"mirror"`/`"stop"`）、`front`/`rear`（面型：`radius_mm`、`conic_k`、`coeffs`（A4,A6,…）、`poly_terms`、`q_asphere`、`xy_polynomial`、`zernike`、`semi_aperture_mm`；`rear` 只属于 `lens`）、`glass`（元件内介质）、`gap_glass`（到下一元件的间隔介质，默认空气——胶合面/双胶合靠它表达）、`thickness_mm`（中心厚度，单面元件必须为 0）、`gap_mm`（到下一元件的轴向距离，**可负** = 轴向折返，如卡塞格林）、`decenter_x_mm`/`decenter_y_mm`、`tilt_x_deg`/`tilt_y_deg`/`tilt_z_deg`（先偏心后倾斜，`geo.CoordBreak` 语义）、`fold`（`0`=auto：倾斜的反射镜把链轴转到反射光束；`1`=off）。

链的语义（`lens/group.go` 的文档与 `Group.Layout`）：面按链序排布，元件的前表面在上一元件末面之后 `gap_mm` 处（沿链轴量）；链轴初始为全局 +z，被倾斜反射镜折转后跟随反射光束；元件的厚度沿元件自身轴量（元件刚性），**posed 元件不会移动后续元件**（自动复位）。反射镜自身不计基板厚度。

响应：

| 字段 | 说明 |
|---|---|
| `surfaces[]` | 逐面（面顺序 = 链序）：`element`（元件下标）、`face`（`"front"`/`"rear"`）、`surface`（该面在 `Group.System()` 中的下标，含坐标断点）、`vertex`/`axis`/`normal`（全局 `{x,y,z}` mm）、`folds`（该面之前的折转次数） |
| `exit_vertex` / `exit_axis` | 末面顶点与出射光束方向（含折转、反射；`Group.ExitFrame`） |
| `folds` | 总折转次数 |
| `num_surfaces` | 光学面数（不含坐标断点） |

    curl -s -X POST localhost:1130/api/lens/group/layout -d '{"elements":[
      {"kind":"stop","front":{"semi_aperture_mm":150},"thickness_mm":0,"gap_mm":60},
      {"kind":"mirror","front":{"radius_mm":-2000,"semi_aperture_mm":153},"thickness_mm":0,"gap_mm":-300},
      {"kind":"mirror","front":{"radius_mm":-646.1538461538461,"semi_aperture_mm":117.5},"thickness_mm":0,"gap_mm":600}]}'

    {"surfaces":[
      {"element":0,"face":"front","surface":0,"vertex":{"x":0,"y":0,"z":0},"axis":{"x":0,"y":0,"z":1},"normal":{"x":0,"y":0,"z":1},"folds":0},
      {"element":1,"face":"front","surface":1,"vertex":{"x":0,"y":0,"z":60}, ...},
      {"element":2,"face":"front","surface":2,"vertex":{"x":0,"y":0,"z":-240}, ...}],
     "exit_vertex":{"x":0,"y":0,"z":-240},"exit_axis":{"x":0,"y":0,"z":1},"folds":0,"num_surfaces":3}

（库例 `library/geo/mirror-cassegrain-f1000-d300-b300.json` 的元件形式：负 gap −300 使副镜落在 z = −240；主镜与副镜都未倾斜，所以链不折转。）

**折转链的适用范围**（与 `lens` 包的实现一致，示例与测试都在这两个范围内构造）：①折转镜处在沿链轴传播的光束里（典型：45° 平面镜折转准直/会聚光束，后续元件落在反射链轴上）；②未倾斜的反射镜用**负 gap** 表达轴向折返（卡塞格林）。另有两点实现边界，`layout` 与成像质心不受影响、但读数要按下面的说明理解：③折转镜后面紧跟“返回型”反射镜时，`Group.ExitFrame()` 的出射方向会残留前一元件坐标系的量，`image`/`dispersion` 的焦面会落在错误平面上（此类组合请改用等价的面序，例如把折转镜放在链尾）；④“返回型”反射镜（如抛物面主镜）之后接折转镜时，链轴与光束反平行，`exit_axis` 与 `focus_mm` 的符号按链轴给出（焦面与像点本身正确，如上面的牛顿式示例：像 (0,−100,−840) 而 `focus_mm` 读作 −100）。这两处都在 `lens/group_build.go` 的 `plan()`/`foldBreak` 内，本轮接线未改 `lens` 包。

### POST /api/lens/group/image

实光线成像位置：在组的光瞳上采样（默认 Hexapolar 3 环 = 37 条），追迹后沿**出射轴**搜索使 RMS 最小的焦面（粗扫描 + 黄金分割，确定性；`lens.Group.ImagePosition`）。对折转、偏心、倾斜的含镜系统同样有效（一阶引擎不适用时它仍可用）。

请求 = group 的 wire form（平铺或 `"group"`）+ 束参数：

| 字段 | 类型 | 默认 | 说明 |
|---|---|---|---|
| `object_z_mm` | number 或 `"inf"` | group 自身的物距 | 物面 z；数字 = 有限物 |
| `field_deg` | number | `0` | 无穷远物的视场角（度，子午 y 方向；\\|θ\\|<90） |
| `field_tan` | number | 由 `field_deg` 换算 | 无穷远物 = 场角切线；有限物 = 物高 mm（给定时优先于 `field_deg`） |
| `pupil_radius_mm` | number | `0`（自动） | 采样光瞳半径 mm；0 = 用光阑元件的半口径，在最前面之前 1 mm、以真实主光线为中心，采样按 0.1 % 内缩 |
| `pupil_z_mm` | number | `0` | 显式光瞳平面的 z（配合 `pupil_radius_mm` 使用） |
| `samples` | int | `3` | Hexapolar 环数（1…8） |
| `lambda_nm` | number | `587.5618` | 波长 nm |
| `range_mm` | number | 0（自动） | 焦面搜索窗半宽 mm（0 = 0.5·\\|估计焦距\\|+10） |
| `scan` | int | 0（65） | 粗扫描点数（≤ 4096） |
| `tol_mm` | number | 0（1e-9） | 黄金分割的轴向容差 mm |

响应（纯数值，可确定性序列化）：

| 字段 | 说明 |
|---|---|
| `point` | 最佳焦面上追迹光线的质心（全局 mm） |
| `rms_mm` / `geo_radius_mm` | 关于质心的 RMS 半径 / 最大半径 mm |
| `axis` / `plane_o` | 焦面法向（= 出射轴）与平面参考点（末面顶点） |
| `focus_mm` | 从 `plane_o` 沿 `axis` 到焦面的有符号距离 mm |
| `rays` / `failed` | 进入统计的光线数 / 未到达末面的采样数 |
| `chief_hit` / `chief_dir` | 主光线（采样中心）的落点与出射方向 |
| `lambda_nm` / `field_tan_y` / `object_z_mm` / `samples` | 回显束参数（`object_z_mm` 为 `null` = 无穷远） |

    curl -s -X POST localhost:1130/api/lens/group/image -d '{"elements":[
      {"kind":"mirror","front":{"radius_mm":-400,"conic_k":-1,"semi_aperture_mm":24}}]}'

    {"point":{"x":0,"y":0,"z":-200.00000000000617},"rms_mm":5.94e-13,"geo_radius_mm":7.38e-13,
     "axis":{"x":0,"y":0,"z":-1},"plane_o":{"x":0,"y":0,"z":0},"focus_mm":200.00000000000608,
     "rays":37,"failed":0, ...}

（抛物面 K=−1 对轴上平行光无球差：像点恰在 f = \\|R\\|/2 = 200 mm，RMS 在搜索容差量级；负号是因为反射光束沿 −z 出射，焦面到顶点的距离沿轴读数为 +200。）

### POST /api/lens/group/afocal

无焦组的角放大率 Γ（实光线，`lens.Group.AngularMagnification`）：以两个倾斜 ±δ = 1e-6 的平行束在弧矢（x）方向测量出射斜率之比；另有 D_in/D_out（束径比，\\|Γ\\| 的独立几何交叉验证）与出射准直残差。**出射束准直不好（> 1e-2 rad）的系统返回 400**（有光焦度的系统没有角放大率）。

请求 = group 的 wire form + `object_z_mm`/`field_tan`/`pupil_radius_mm`/`samples`/`lambda_nm`（同上）；对象必须是无穷远物（`field_tan` 对 Γ 无影响，Γ 由弧矢方向定义）。

响应：`gamma`、`beam_ratio`、`spread_rad`、`rays`、`entry_dir`/`exit_dir`（探测束入射/出射方向）、`axis_in`/`axis_out`、`lambda_nm`、`samples`。

    curl -s -X POST localhost:1130/api/lens/group/afocal -d '{"group":{...开普勒望远镜...}}'
    {"gamma":-1.0000000000201266,"beam_ratio":0.9940377442497541,"spread_rad":3.497e-4,
     "rays":37,"axis_out":{"x":0,"y":0,"z":1}, ...}

（两片 f = 100 mm 的薄透镜组成开普勒望远镜：Γ = −1；在准直空间插入 45° 平面折转镜后 Γ 逐位不变，出射轴转为 +y，`TestGroupAngularMagnificationFoldInvariant`。）

### POST /api/lens/group/dispersion

两个波长的最佳焦面轴向差（`lens.Group.Dispersion`，同一束参数与搜索选项）。

请求 = group 的 wire form + 束参数 + `lambda_a_nm`（默认 F 线 486.1327）/`lambda_b_nm`（默认 C 线 656.2725）。

响应：`shift_mm`（`a.focus_mm − b.focus_mm`，沿出射轴）、`lambda_a_nm`/`lambda_b_nm`、`a`/`b`（两次 `image` 的完整响应）。**纯反射镜系统（介质恒定）的位移恰好为 0**：每个波长的追迹逐位相同。

    curl -s -X POST localhost:1130/api/lens/group/dispersion -d '{"group":{...双胶合...}}'
    {"shift_mm":-0.4065509316004352,"lambda_a_nm":486.1327,"lambda_b_nm":656.2725,
     "a":{"point":{"x":0,"y":0,"z":73.6214...},"focus_mm":64.6214...},"b":{...,"focus_mm":65.0279...}}

### POST /api/lens/group/set

把一串参数编辑应用到 group 上（按顺序，`lens.Group.SetParam`），返回更新后的 group。参数语法与 CLI 的 `-p` 相同：`lens.GroupParam.String()` 的写法（`e0.front.R`、`e0.rear.c`、`e0.front.K`、`e0.front.A4`、`e0.front.poly0`/`q0`/`xy0`/`z1`、`e0.front.semi`、`e0.thk`、`e0.gap`、`e0.z`、`e0.dx`/`dy`、`e0.tx`/`ty`/`tz`）以及长别名（`radius_mm`、`curvature`、`conic_k`、`semi_aperture_mm`、`thickness_mm`、`gap_mm`、`position_mm`、`decenter_x_mm`、`tilt_x_deg`，…）；面向 `f`/`r`/`1`/`2` 与 `front`/`rear` 等价。`R = 0`（或 ±Inf）表示平面，`gap_mm` 可为负。

请求：`{"group": <group wire form>, "sets": [{"param": "e1.front.R", "value": -1600}, …]}`（`sets` 最多 512 项）。响应：`{"group": <更新后的 group>, "applied": n}`。

编辑失败返回 400，**错误文本就是 `lens` 包的文本**（含元件下标与原因），例如：

    {"error":"lens: parameter e0.rear.R: lens: element 0 (stop) has no rear surface"}

    curl -s -X POST localhost:1130/api/lens/group/set -d '{"group":{...},"sets":[{"param":"e1.front.R","value":-1600}]}'

### POST /api/lens/group/prescription

把元件级的 group 转成「卡片链」处方（`lens.Group.Prescription`）：纯轴向或轴向折返（负厚度，如卡塞格林）的布局可以精确表达，非末元件的 3D 折转、含中心厚度的倾斜两面元件等不可表达的布局返回 400 并给出原因。响应：`{"prescription": <处方 JSON>}`。

    curl -s -X POST localhost:1130/api/lens/group/prescription -d '{"elements":[...卡塞格林...]}'
    {"prescription":{"name":"…","object_z_mm":"inf","cards":[
      {"radius_mm":0,"thickness_mm":60,"semi_aperture_mm":150},
      {"radius_mm":-2000,"thickness_mm":-300,"semi_aperture_mm":153,"mirror":true},
      {"radius_mm":-646.1538461538461,"thickness_mm":600,"semi_aperture_mm":117.5,"mirror":true}],"stop":0}}

（输出与 `library/geo/mirror-cassegrain-f1000-d300-b300.json` 的卡片逐位一致——`cmd/oz` 的 `TestLensGroupToPrescriptionCLI` 就是这么断言的。反向的 处方→group 在 `lens` 包内（`Prescription.Group()`），暂未开放为端点。）

### POST /api/lens/parax

处方的一阶（近轴）报告，**含反射镜**（`geo.System.Paraxial` 的展开法：镜面按 ū' = −ū − 2ncy 处理，负厚度写折返；`Mirrors`/`Parity`/`Afocal` 与角放大率为本次接线新增的公开量）。与 `oz lens parax` 同源。偏心/倾斜（posed）系统仍返回 400（一阶引擎只接受共轴系统，含 z 轴 roll）。

请求：

| 字段 | 类型 | 默认 | 说明 |
|---|---|---|---|
| `prescription` | object | 必填 | 处方 JSON |
| `object_z_mm` | number 或 `"inf"` | 处方自身的物距 | 物面位置 mm（`"inf"` = 无穷远） |
| `field_tan` | number | — | 参考光线的场切线（afocal 的 Γ 需要 ≠0；有限物的 Γ 与它无关） |
| `lambda_nm` | number | `587.5618` | 波长 nm |

响应（不存在或非有限的量用 `null`，并在对应的 `*_note` 里说明原因——JSON 没有 NaN/Inf）：

| 字段 | 说明 |
|---|---|
| `lambda_nm`、`object_z_mm` | 回显（`null` = 无穷远物） |
| `mirrors`、`parity`、`afocal` | 反射面数、出射传播方向奇偶性 `(-1)^mirrors`、\\|C\\| ≤ 1e-12 /mm |
| `efl_mm`、`focal_length_mm` | 带符号 EFL（凹镜 R<0 时 EFL = R/2 < 0）与习惯焦距 `Parity·EFL` |
| `ffl_mm`、`bfl_mm`、`image_z_mm` | 前/后焦距与轴向像面位置（同一 z 轴，可直接与实光线比较） |
| `magnification` | 横向放大率（无穷远物为 0，引擎约定） |
| `entrance_pupil_z_mm`/`entrance_pupil_radius_mm`、`exit_pupil_z_mm`/`exit_pupil_radius_mm` | 入/出瞳位置与半口径 |
| `gamma` | 角放大率：afocal 用主光线比值；有光焦度的系统在有限物距下给 γ = n/(n'·m)。不可用时 `null` |
| `gamma_note` | Γ 未给出的原因（afocal 缺 `field_tan`、无穷远物 + 有光焦度、posed 系统…） |
| `fc_shift_mm` | F–C 焦移（`FocusAt(F) − FocusAt(C)`，mm）；**纯镜系统恰好 0** |
| `fc_note` | 焦移未给出的原因（afocal 系统没有像面等） |

    curl -s -X POST localhost:1130/api/lens/parax -d '{"prescription":{...卡塞格林...}}'
    {"lambda_nm":587.5618,"object_z_mm":null,"mirrors":2,"parity":1,"afocal":false,
     "efl_mm":-857.1428571428569,"focal_length_mm":-857.1428571428569,"ffl_mm":1713.0612244897957,
     "bfl_mm":-599.9999999999998,"image_z_mm":-839.9999999999998,"magnification":0,
     "entrance_pupil_z_mm":0,"entrance_pupil_radius_mm":150,
     "exit_pupil_z_mm":-728.5714285714287,"exit_pupil_radius_mm":214.2857142857143,
     "gamma":null,"gamma_note":"无穷远物 + 有光焦度的系统没有角放大率（Γ）：给 --object-z 有限物距或换 afocal 系统",
     "fc_shift_mm":0}

单镜（R = −2000、空气）的闭式：`EFL = R/2 = −1000`；物在 z = −3000 时 `image_z_mm = −1500`、`magnification = −0.5`、`gamma = −2 = 1/m`（`TestLensParaxCLI`）。

### POST /api/geo/mtf

几何 MTF 与场曲摘要（G11；`oz geo mtf`）。每个视场用等面积圆盘采样（`geo.EqualAreaDisc`）在近轴入瞳上追迹真实光线，把命中像面的坐标分成 LSFT（y，子午）与 LSFS（x，弧矢），|Σ p_k·exp(−2πi f x_k)| 即该频率的几何 MTF（`geo.BinLSF`/`geo.LSF.MTF`）；每个视场按自己的光斑范围独立分箱，响应给出该视场的 Δ 与 Nyquist（`bin_width_*`、`nyquist_*`），请求频率超过 Nyquist 时数值仍给出但已混叠。场曲用 `geo.System.FieldFoci`（入瞳 ±rim_frac 的子午/弧矢光线对）给出 T/S 焦点，`sag_*` 为相对像面的离焦，Petzval 面取 (3S−T)/2。

请求：

| 字段 | 类型 | 默认 | 说明 |
|---|---|---|---|
| `prescription` | object | 必填 | 处方 JSON |
| `lambda_nm` | number | `587.5618` | 波长 nm |
| `fields_deg` | number[] | `[0]` | 子午面视场角（度，\|θ\|<90）；有限物距时按入瞳处的场角换算物高 |
| `cycles_per_mm` | number[] | `[10,20,30,40,50]` | MTF 频率（cycles/mm，非负） |
| `rays` | int | `2000` | 每视场光线数（4…200000） |
| `bins` | int | `201` | LSF 分箱数（2…4096） |
| `image_z_mm` | number | 近轴像面 | 像面轴向位置 mm（0 = 近轴像面） |
| `rim_frac` | number | `1e-3` | 场曲光线对的入瞳半径比例 |

响应：`lambda_nm`、`image_z_mm`、`rays`、`bins`、`cycles_per_mm`、`fields`（每项：`theta_deg`、`tan_y`、`rays_hit`/`rays_total`、`centroid_x_mm`/`centroid_y_mm`、`rms_radius_mm`（关于主光线）、`bin_width_tangential_mm`/`bin_width_sagittal_mm`、`nyquist_tangential_per_mm`/`nyquist_sagittal_per_mm`、`mtf_tangential[]`、`mtf_sagittal[]`）、`foci`（每项：`theta_deg`、`z_tangential_mm`/`z_sagittal_mm`、`y_tangential_mm`/`y_sagittal_mm`、`sag_tangential_mm`/`sag_sagittal_mm`/`sag_petzval_mm`）、`foci_note`（场曲不可用时说明原因：有限物距、偏心/倾斜或非近轴系统）。

### POST /api/geo/ghost

双反射鬼像路径表（G10；`oz geo ghost`）。`geo.System.GhostAnalysis` 枚举每对 0 ≤ J < I ≤ n−1 的反射面组合，相对能量取法向入射 Fresnel 乘积 R_I·R_J（可被 `reflectance` 覆盖），按能量降序（同能量按面索引升序）排序，并给出近轴落点（反向 ABCD 链）与小孔径采样（Hexapolar）的真实杂散光斑。

请求：

| 字段 | 类型 | 默认 | 说明 |
|---|---|---|---|
| `prescription` | object | 必填 | 处方 JSON（面数 n 满足 n(n−1)/2 ≤ 512） |
| `lambda_nm` | number | `587.5618` | 波长 nm |
| `field_deg` | number | `0` | 子午视场角（度） |
| `image_z_mm` | number | 近轴像面 | 杂散光斑评估的像面 |
| `pupil_radius_mm` | number | 入瞳半径×`1e-3` | 采样光瞳半径 mm |
| `rings` | int | `1` | Hexapolar 采样环数（1…8） |
| `min_relative` | number | `0` | 相对能量阈值，低于它的路径不列出 |
| `max_paths` | int | `0`（全部） | 最多列出的路径数 |
| `reflectance` | (number\|null)[] | 空（Fresnel） | 每面强度反射率（长度 = 面数；`null` 项回落 Fresnel 值） |

响应：`lambda_nm`、`image_z_mm`、`field_deg`、`pupil_z_mm`、`pupil_radius_mm`、`rings`、`min_relative`、`paths_enumerated`、`paths`（每项：`pair`（`reflect(I=2, J=1)`）、`i`/`j`、`relative_energy`、`reflect_i`/`reflect_j`、`paraxial_x_mm`/`paraxial_y_mm`（非同轴/无法一阶分析时为 `null`）、`chief_x_mm`/`chief_y_mm`、`chief_status`、`spot_rays`/`spot_rays_total`、`spot_centroid_x_mm`/`spot_centroid_y_mm`、`spot_rms_mm`、`spot_geo_mm`、`note`（非空：镜面/坐标断点等无法分析的路径））。

    curl -s -X POST localhost:1130/api/geo/ghost -d '{"prescription":{…},"min_relative":1e-3,"max_paths":8}'

### POST /api/doe/focus

N 级量化相位聚焦 DOE 的一级效率与焦点偏差（G8；`oz doe focus`）。调用 `doe.SimulateFocus`：M×M 网格上的量化透镜透过率 + 带限 ASM 传播到设计焦面求一级效率（窗口内功率/入射功率），并在轴上强度谱上扫描出轴向峰值（抛物线细化）。单位：接口外为 mm/nm，内部按 doe 包的 µm 工作。

请求：

| 字段 | 类型 | 默认 | 说明 |
|---|---|---|---|
| `focal_mm` | number | 必填，>0 | 设计焦距 mm |
| `levels` | int | 必填 | 相位级数 N（≥2，≤4096） |
| `lambda_nm` | number | `632.8` | 波长 nm |
| `waist_um` | number | 由 `fresnel_number` 给出 | DOE 上高斯照明 1/e 振幅半径 µm |
| `fresnel_number` | number | `32` | N_F = w0²/(λf)，用于换算 w0 = sqrt(N_F·λ·f) |
| `grid_size` | int | `1024` | 网格 M（2 的幂，16…2048） |
| `pixel_um` | number | `8·w0/M` | DOE 与焦面的像素间距 µm（默认与 doe 包测试同配置） |
| `window_radii` | number | `2.5` | 效率窗口半径（焦斑半径的倍数） |

响应：`focal_mm`/`focal_um`、`lambda_nm`、`levels`、`grid_size`、`pixel_um`、`waist_um`、`fresnel_number`、`spot_radius_um`（w_f = λF/(πw0)）、`eta_first_order`（数值一级效率）、`eta_analytic`（η_N = [sin(π/N)/(π/N)]²）、`eta_rel_err`、`axial_peak_um`、`focus_deviation_um`/`focus_deviation_px`（轴向峰值 − 设计焦面）、`peak_offset_px`（焦面强度峰值的横向偏移）、`window_fraction`、`incident_power`。默认网格下 |η − η_N| ≈ 3.5e-4（M=2048 时 ≈ 9e-5）。

### POST /api/coating/design

内置膜系模板设计。

    {"template":"broadband-ar", "lambda_nm":550, "from_nm":420, "to_nm":680,
     "pairs":12, "cutoff_nm":0, "incident":"Air", "substrate":"BK7"}

模板：`single-ar`、`vcoat`、`broadband-ar`、`hr`（高反，`pairs` 为 H/L 对数）、`lwp`、`swp`（长/短波通，需要 `cutoff_nm > 0`）、`splitter`（50:50）。未知模板、材料或截止波长无效返回 400。响应：`{"stack": {...}, "summary": "…", "merit": number}`。

### POST /api/coating/spectrum

请求：`{"stack": <stack>, "from_nm":400, "to_nm":700, "step_nm":5, "angle_deg":0, "pol":"u"}`。`pol` 取 `s`、`p`、`u`（缺省为非偏振 `u`）。点数 `(to-from)/step+1` 最多 10000，超出返回 400。响应：`{"rows": [{"lambda_nm", "R", "T", "A"}, …]}`（端点含两端，共 `(to-from)/step+1` 行）。

## 5. 波动与量子（WOS，同一服务）

以下端点的完整语义见 `docs/kernel/API.md`（字段与示例）与 `docs/kernel/INTEGRATION.md`（接入）；本节列出路径与要点。

| 方法 | 路径 | 说明 |
|---|---|---|
| GET | `/api/catalog` | 元件/光源/算法/偏振/量子门与态/示例目录（GUI 由它驱动）。脚本元件带 `custom:true` 与 `source` |
| GET | `/api/health` | 存活与 GPU 状态：`{"status":"ok","time":…,"gpu":{available,enabled,backend,attempts,failures,degraded,last_error?,last_error_time?}}` |
| GET | `/api/elements` | 已加载的脚本元件（名字、标签、行为、类别、参数、表达式、来源）与 `version` |
| POST | `/api/elements/reload` | 重新扫描定义目录，返回 `{dirs, loaded, skipped, errors, version}` |
| GET | `/api/elements/{name}/preview.png` | 元件掩膜预览（PNG）。`kind=amp|phase`、`size`、`width`、`wl`、`mask`，其余查询参数作为元件参数 |
| POST | `/api/validate` | 配置校验：`{"ok":bool,"issues":[{"path","message"}]}` |
| POST | `/api/convert` | 旧的 `elements` 序列转为定位场景：`{"scene":{…},"sources":[…]}`，不计算 |
| POST | `/api/simulate` | 提交仿真，202 返回 `{"run_id","status":"running"}`；串行执行，排队上限 8 |
| GET | `/api/runs/{id}` | 状态（`running`/`done`/`error`）与元数据（网格、波长、平面列表、警告、场景） |
| GET | `/api/runs/{id}/planes/{pid}` | 平面场：`fmt=bin`（float32 小端 N×N 行主序，`X-Grid-Size` 头）或 `fmt=png` |
| GET | `/api/runs/{id}/profiles/{pid}` | 一维剖面 JSON：`axis=x|y`、`field`、`coord`、`part` |
| GET | `/api/runs/{id}/inspect/{pid}` | 像素读数：琼斯矢量、斯托克斯参数、偏振椭圆、相位与波前 PV/RMS。`x`、`y`（像素）缺省为强度质心 |
| GET | `/api/runs/{id}/scene` | 路由出的光路几何与光源颜色；仅当元数据含 `scene` 时存在，否则 404 |
| POST | `/api/quantum` | Fock 基线性光学同步计算（`fmt=png|svg` 返回图表） |

仿真配置（`/api/validate`、`/api/simulate`、`/api/convert`）的字段见 `optics.Config`（`docs/kernel/INTEGRATION.md` §1.1）。`/api/simulate` 校验失败返回 400 并给出第一个错误；完整错误列表用 `/api/validate`。

平面视图参数（`/api/runs/{id}/planes/{pid}`）：`field`（`total`、`amplitude`、`ex`、`ey`、`ez`、`phase_x|y|z|u`、`pol_azimuth|ellip|s1|s2|s3|degree|intensity`、`color`，默认 `total`；`color` 仅 `fmt=png`）、`part`（相干单元序号，默认合成）、`scale`（`lin|log`）、`cmap`（`inferno|phase|gray|diverging`）、`mask`、`exposure`、`gamma`、`pmin`、`pmax`。

量子请求：

    {"modes":2, "cutoff":4,
     "state":{"type":"fock","params":{"occupation":[1,1]}},
     "gates":[{"type":"beam_splitter","params":{"mode0":0,"mode1":1,"reflectivity":0.5}}]}

响应为 `QuantumResult`（`mean_photons`、`g2`、`photon_distributions`、`quadratures`、`joint_distributions`、`joint_full`（全模式联合分布，little-endian 下标）、`norm`）。限制：`(cutoff+1)^modes ≤ 2^20`。

## 6. 最小示例

    # 服务
    oz serve -addr :1130 &
    curl -s localhost:1130/api/version
    # 几何：设计 → 评估
    curl -s -X POST localhost:1130/api/lens/design -d '{"efl_mm":100}' > d.json
    curl -s -X POST localhost:1130/api/lens/eval -d "$(python3 -c 'import json;print(json.dumps(json.load(open("d.json"))["prescription"]))')"
    # 分析：热漂移 / 灵敏度（可加 "inverse":true,"budget_mm":0.005）/ 几何 MTF / 鬼像路径
    P="$(python3 -c 'import json;print(json.dumps(json.load(open("d.json"))["prescription"]))')"
    curl -s -X POST localhost:1130/api/lens/thermal -d "{\"prescription\":$P,\"delta_t_k\":10}"
    curl -s -X POST localhost:1130/api/lens/sens -d "{\"prescription\":$P}"
    curl -s -X POST localhost:1130/api/lens/parax -d "{\"prescription\":$P}"          # 一阶（含镜）
    curl -s -X POST localhost:1130/api/geo/mtf -d "{\"prescription\":$P,\"fields_deg\":[0,0.5],\"cycles_per_mm\":[10,30]}"
    curl -s -X POST localhost:1130/api/geo/ghost -d "{\"prescription\":$P}"
    # 透镜组（3D，含镜）：布局 / 成像位置 / 角放大率 / 色散 / 编辑 / 转处方
    G='{"name":"newtonian-45","elements":[{"kind":"stop","front":{"semi_aperture_mm":20},"gap_mm":60},
        {"kind":"mirror","front":{"radius_mm":-2000,"conic_k":-1,"semi_aperture_mm":40},"gap_mm":-900},
        {"kind":"mirror","front":{"radius_mm":0,"semi_aperture_mm":14},"gap_mm":120,"tilt_x_deg":45}]}'
    curl -s -X POST localhost:1130/api/lens/group/layout -d "$G"          # 顶点 (0,0,0)/(0,0,60)/(0,0,-840)
    curl -s -X POST localhost:1130/api/lens/group/image -d "$G"           # 像 (0,-100,-840)，RMS ~6e-14
    curl -s -X POST localhost:1130/api/lens/group/set -d "{\"group\":$G,\"sets\":[{\"param\":\"e0.gap_mm\",\"value\":50}]}"
    curl -s -X POST localhost:1130/api/lens/group/dispersion -d "$G"      # 纯镜：shift_mm 恰为 0
    curl -s -X POST localhost:1130/api/lens/group/afocal -d '{"elements":[
        {"kind":"stop","front":{"semi_aperture_mm":4},"gap_mm":0},
        {"kind":"lens","front":{"radius_mm":51.68000345005885,"semi_aperture_mm":20},
         "rear":{"radius_mm":0,"semi_aperture_mm":20},"glass":"N-BK7","thickness_mm":2,"gap_mm":198.68},
        {"kind":"lens","front":{"radius_mm":51.68000345005885,"semi_aperture_mm":20},
         "rear":{"radius_mm":0,"semi_aperture_mm":20},"glass":"N-BK7","thickness_mm":2}]}'   # Γ = −1
    curl -s -X POST localhost:1130/api/lens/group/prescription -d "$G"
    # DOE（G8）：一级效率与解析 η_N 的对比
    curl -s -X POST localhost:1130/api/doe/focus -d '{"focal_mm":10,"levels":4}'
    # 膜系光谱
    curl -s -X POST localhost:1130/api/coating/design -d '{"template":"single-ar","lambda_nm":550}' > c.json
    curl -s -X POST localhost:1130/api/coating/spectrum -d "{\"stack\":$(python3 -c 'import json;print(json.dumps(json.load(open("c.json"))["stack"]))'),\"from_nm\":500,\"to_nm\":600,\"step_nm\":10}"
    # 波动：量子 HOM
    curl -s -X POST localhost:1130/api/quantum -d '{"modes":2,"cutoff":4,"state":{"type":"fock","params":{"occupation":[1,1]}},"gates":[{"type":"beam_splitter","params":{"mode0":0,"mode1":1,"reflectivity":0.5}}]}'

端到端的端点检查见 `cmd/oz/service_test.go` 与 `server/endpoints_test.go`（httptest），以及 `docs/VERIFICATION.md`。
