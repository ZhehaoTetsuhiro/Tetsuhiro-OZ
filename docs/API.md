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
| `/api/version`、`/api/lens/design`、`/api/lens/eval`、`/api/lens/tol`、`/api/coating/design`、`/api/coating/spectrum` | OZ | 精确路径，优先于 `/api/` 前缀 |
| `/api/...`（其余） | WOS | 仿真 API，由 `server.Server.Handler()` 处理；未知路径 404 |
| `/oz/` | 重定向 | 301 → `/#oz`（合并前的几何/膜系页面地址） |
| `/` 及静态文件 | GUI | 内嵌在二进制中：`index.html`、`app.js`、`oz.js`、`style.css`、`favicon.svg`、字体 |

## 3. 约定

- 请求与响应：JSON（`Content-Type: application/json; charset=utf-8`），除非端点另有说明。
- 错误：`{"error": "…"}`，状态码 400（请求无效）、404（资源不存在）、405（方法不允许）、409（运行未完成）、429（排队已满）、500（内部错误或响应含 NaN/Inf）。
- 方法：计算端点为 POST，查询端点为 GET；错误方法返回 405。`GET /api/health` 不限制方法（既有行为）。
- 工作量上限（超出即 400，计算前拒绝）：光谱点数 ≤ 10000；公差运行次数 ≤ 20000；设计迭代 ≤ 10000。仿真网格与元件数由 `/api/validate` 检查。
- 请求体上限：`/api/quantum` 1 MB；`/api/validate`、`/api/simulate`、`/api/convert` 2 MB；`/api/lens/design` 与 `/api/lens/eval` 1 MB；`/api/lens/tol` 4 MB；`/api/coating/design` 1 MB；`/api/coating/spectrum` 8 MB。
- `GET /api/version` 返回的 `version` 与仓库根目录 `VERSION` 文件一致（由 `cmd/oz` 的测试保证）。

## 4. 几何与膜系（OZ）

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
    # 膜系光谱
    curl -s -X POST localhost:1130/api/coating/design -d '{"template":"single-ar","lambda_nm":550}' > c.json
    curl -s -X POST localhost:1130/api/coating/spectrum -d "{\"stack\":$(python3 -c 'import json;print(json.dumps(json.load(open("c.json"))["stack"]))'),\"from_nm\":500,\"to_nm\":600,\"step_nm\":10}"
    # 波动：量子 HOM
    curl -s -X POST localhost:1130/api/quantum -d '{"modes":2,"cutoff":4,"state":{"type":"fock","params":{"occupation":[1,1]}},"gates":[{"type":"beam_splitter","params":{"mode0":0,"mode1":1,"reflectivity":0.5}}]}'

端到端的端点检查见 `cmd/oz/service_test.go` 与 `server/endpoints_test.go`（httptest），以及 `docs/VERIFICATION.md`。
