# 图形界面（GUI）

`oz serve` 内嵌一个单页界面，地址 `/`。波动光学、量子光学、几何/膜系三个模式共用同一页面，由顶部的模式按钮（或键盘 `m`）切换。界面文件位于 `server/web/`（`index.html`、`app.js`、`oz.js`、`style.css`、`favicon.svg`、`fonts/`），编译进二进制，运行时不需要外部文件。

## 1. 为什么是标签，而不是独立的 `/oz/` 页面

采用同一页面的模式标签：

- 方案 §6.3 与 §8.2 要求单套 GUI（「GUI 合一」）；几何与波动共享一次启动、一个地址、一套快捷键习惯。
- 几何/膜系前端 `oz.js` 是独立 IIFE：不向全局作用域写入任何名字，与波动前端 `app.js`（全局作用域）互不干扰。
- 旧地址 `/oz/` 以 301 跳转到 `/#oz`，直接打开几何/膜系标签；旧页面的 HTML 已删除。

## 2. 模式

| 按钮 | 模式 | 左侧 | 右侧 |
|---|---|---|---|
| 波动光学 | `wave` | 光路场景（元件、位置）、光源、全局参数 | 图像/相位/偏振/剖面/立体视图 |
| 量子光学 | `quantum` | 态、门序列 | 光子数分布、g²、正交分量、联合分布 |
| 几何/膜系 | `oz` | 消色差设计参数、公差参数、透镜组（3D）面板、膜系模板 | 设计结果与处方 JSON、公差报告、透镜组布局画布与测量、膜系与光谱 |

`#oz` 哈希让页面直接进入几何/膜系模式（`/oz/` 的跳转目标）。

几何/膜系模式的操作（全部调用服务端 API，见 `docs/API.md` §4）：

1. **设计**：`POST /api/lens/design`，显示报告、优化器统计与处方 JSON；设计成功后启用「评估」「运行公差」与「高级分析」里的按钮。
2. **评估当前设计**：`POST /api/lens/eval`，显示 EFL、BFL、F/#、F–C 色差。
3. **运行公差**：`POST /api/lens/tol`，参数为半径/厚度/Δn 与运行次数；显示良率、各分布统计与失败原因。焦距判据取自设计时的目标 EFL。
4. **高级分析（折叠区）**：面板由 `oz.js` 在初始化时构建（`index.html` 不感知），复用页面既有的 `ozGrid/ozField/ozActions/ozDetails/ozOut` 类，因此手机断点与触控目标规则同样生效。「热漂移」调用 `POST /api/lens/thermal`（ΔT 与镜筒 α 取自面板输入，λ 用 d 线），显示 EFL/BFL 变化与 mm/K；「灵敏度与公差分配」调用 `POST /api/lens/sens`（`metric=merit`、`inverse=true`、预算取自面板），显示逐参数 ∂M/∂x、分配的 ±公差与封顶标记。两个按钮在首次设计成功前禁用。
5. **透镜组（3D，含镜）（折叠区）**：见 §2.1。
6. **设计膜系**：`POST /api/coating/design`，显示膜层摘要。
7. **设计并画光谱**：先设计膜系，再 `POST /api/coating/spectrum`（400–750 nm，2 nm 步长，非偏振），在画布上绘制 R(λ)。画布按容器宽度与设备像素比重绘。

### 2.1 透镜组（3D，含镜）

元件级透镜组（`lens.Group`：透镜/反射镜/光阑、前后表面面型与圆锥/非球面参数、玻璃与间隔介质、中心厚度、可负的 gap、偏心与倾斜、fold 折转）的交互面板，与 `docs/API.md` §4 的 `/api/lens/group/*` 一一对应。面板与表格同样由 `oz.js` 在初始化时构建，`index.html` 只保留左侧的空容器。

- **载入示例**：下拉里有三个构造好的 group（JSON 直接内嵌在 `oz.js` 里，`cmd/oz` 的 `TestGroupExamplesInGUI` 按 `/*oz:example*/` 标记提取并逐个解码 + `Validate()` + `System()` 校验）：卡塞格林（两镜，库例 F1000 d300 b300）、BK7/F2 双胶合 f=100 mm（`oz design achromat` 的设计值）、折叠式开普勒（45° 折转，Γ=−1）。载入后自动调用一次「布局」并绘出 3D 光轴。
- **元件表格**：每个元件一块，含类型（lens/mirror/stop）、前后半径与圆锥常数、半口径、玻璃、间隔介质 `gap_glass`、中心厚度、到下一元件的 gap、偏心 dx/dy、倾斜 tx/ty/tz、折转 fold。空的输入框表示「不修改」；在空的后表面 R 里填一个数会给单面元件补上后表面。**类型/玻璃/间隔介质/fold 由前端写进 group JSON**（`frames` 之类不能由数值参数表达的字段），数值输入作为 `sets` 发到 `POST /api/lens/group/set` 由服务端逐项校验 —— 编辑失败时**服务端错误文本原样显示**在结果区（如 `lens: parameter e0.rear.R: lens: element 0 (stop) has no rear surface`），返回的更新后 group 会重建整张表格。
- **布局**：`POST /api/lens/group/layout`，结果表给出逐面「顶点/轴向/折转」，并在画布上绘制：顶点连线的链轴（灰）、每个面的局部轴（蓝实线）与顶点法向（橙虚线）、出射方向（绿虚线）、原点处的 x/y/z 三向标。**画布支持拖拽旋转**（Pointer Events，鼠标与触摸同样有效；`pointerdown` 按下、`pointermove` 转身、`pointerup/leave` 结束），另有「等轴测 / XZ 面 / XY 面」三个预设视角按钮；缩放自动适配包围盒，`devicePixelRatio` 后备像素，`oz:mode`/`resize`/`orientationchange` 时重绘。
- **成像位置 / 角放大率 / 色散 / 转处方**：分别调用 `POST /api/lens/group/image`（像点、RMS、焦面距离、主光线）、`/afocal`（Γ、束径比、准直残差）、`/dispersion`（F–C 焦移与两个焦面）、`/prescription`（卡片链处方 JSON）。束参数用面板上的视场角、波长与光瞳半径（0 = 自动用人眼口径/光阑口径）。非无焦系统的「角放大率」、不可表达为卡片的「转处方」都会把服务端的原因显示在结果区。
- **分派到元素表格的量**：kind/glass/gap_glass/fold 走 JSON，其余走 `/set` 的 `sets`（参数名与 CLI 的 `-p` 同语法，见 `docs/API.md` §4）。手机断点（`.ozGrid` 单列、按钮 44 px、画布 `clamp(180px, 44vw, 340px)` + `touch-action: none`）与前两个折叠区完全一致。

### 2.2 波动模式的「参数扫描」面板

左栏「光源 / 全局参数」之下多出一节 **参数扫描**（`server/web/scan.js`，独立 IIFE；`index.html` 只多一行 `<script>`）。它把当前配置里所有**可扫的数值字段**列成下拉（全局波长、每个光源的参数、每个元件的 roll/yaw/pitch/位置、形状参数与元件参数），用户填一组取值、选一个指标，点「运行扫描」即逐点重复仿真：

- **取值**框接受两种写法：`0,15,30,45`（逗号列表）与 `0-90-5`（起-止-步长；`5e-7-8e-7-1.5e-7` 这类科学计数法同样支持）。
- **指标**下拉：功率、峰值、质心 x/y、RMS x/y、斯特列尔比、最大强度。
- 结果区给逐点表格（含逐点错误，某取值让配置非法时只该点报错，其余继续），下面是把指标随参数变化的折线画布；「导出 CSV」保存 `value,metric,error` 三列。
- 面板每秒钟刷新一次参数列表（配置异步载入、用户增删元件/光源后自动跟上），已选项会保留。
- 典型用法：模板「扩展光源 + 旋转单缝（部分相干）」→ 参数选 `scene.components[0].roll` → 取值 `0-90-15` → 指标 `rms_y`，即可定量看到远场条纹随缝旋转由横转竖。

## 3. 键盘

| 键 | 作用 |
|---|---|
| `m` | 循环模式：波动 → 量子 → 几何/膜系 → 波动 |
| `?` / `F1` | 快捷键帮助（三种模式通用） |
| `Esc` | 关闭对话框 |
| `Tab` | 在控件间移动焦点 |

几何/膜系模式下，波动的快捷键（空格运行、数字键切视图、`i` 插入元件等）不生效；空格不会触发波动仿真。输入框内的按键不触发快捷键。

## 3b. AI 助手（右侧栏，默认隐藏）

顶栏「✦ AI 助手」按钮开关右侧栏（也可用 `?ai=1` / `#ai` 直接打开，便于存书签）。三种模式都能用：波动模式解释元件与检查配置、量子模式读门序列、几何/膜系模式读设计与光谱结果文本。

面板结构（`index.html` 的 `#aiPanel`，逻辑在独立的 `ai.js` IIFE 里）：

1. **提供方设置**（折叠区）：提供方预设（OpenAI / DeepSeek / 通义千问 / Moonshot / 智谱 / 本地 llama.cpp / Ollama / 自定义）、Base URL、模型、API Key、温度，加一个「测试连接」按钮。设置存在浏览器 `localStorage`（键 `oz.ai.settings`），**不随配置 JSON 导出、不进服务端日志**；服务端只在每次提问时把它转发给用户自己配置的端点。
2. **规范摘要**：一行「规范摘要 N 字符 · M 节」，点击展开全文。摘要由服务端 `GET /api/ai/context` 现算，来源就是本二进制的元件目录，因此回答「AI 到底读了什么」不需要猜。
3. **对话区**：`#aiLog` 自己滚动（`flex:1` + `overflow-y:auto`），不牵动左右两栏；流式回答带光标，`aria-live="polite"`。
4. **快捷提问**：解释当前元件 / 检查当前配置 / 设计光路 / 新元件定义——把问题填进输入框（不自动发送，可先改）。
5. **输入框**：Enter 发送、Shift+Enter 换行；发送中可「停止」（AbortController 取消请求），「清空对话」重开一轮。

助手的回复里出现 ```` ```oz-config ```` 代码块时，GUI 在代码块上放一个「**应用到场景**」按钮：先 `POST /api/validate` 校验，通过后走与「打开 JSON 文件」完全相同的 `adoptConfig()` 路径（规范化、补默认光源、渲染、运行）。出现 ```` ```oz-element ```` 时代之以「**写入元件定义**」：调 `POST /api/ai/element`，内核校验通过才落盘到元件定义目录，随后重载——元件立刻出现在插入列表（带「· 脚本」标记），无需重启。校验失败的元件定义会被 400 拒绝并显示内核给出的原因，磁盘上不会留下坏文件。

助手读到的「当前状态」由 `ai.js` 只读收集：波动模式的 `S.config`、当前选中元件/光源及其文档、`S.meta` 的平面与警告、`#warnings` 文本；量子模式的 `S.qconfig`；几何/膜系模式直接读 `oz.js` 渲染好的结果文本（`#ozLensOut`、`#ozCoatOut`）。`app.js` 与 `oz.js` 都不需要为 AI 改一行代码。

按键隔离：焦点落在 `#aiPanel` 内时，`app.js` 的全局键盘处理器直接返回——在助手面板里按空格不会触发波动仿真、按 `m` 不会切模式、按 `i` 不会弹出插入元件。

## 4. 移动端

页面面向手机与平板：

- `<meta name="viewport" content="width=device-width, initial-scale=1">`。
- 宽度 ≤ 700 px：页头按钮换行；左右两栏堆叠为单列，页面整体纵向滚动，不出现横向滚动条。助手栏同样占满整行（`#aiPanel { width:100% }`），对话区限高 60vh。
- 触控设备（`pointer: coarse`）与几何/膜系标签内的按钮、输入框最小高度 44 px。
- 光谱画布高度 `clamp(150px, 38vw, 300px)`，宽度跟随容器，按 `devicePixelRatio` 设置后备像素；透镜组 3D 画布同样处理（`clamp(180px, 44vw, 340px)`，`touch-action: none` 让拖拽旋转在触屏上可用）。

验证方法（浏览器开发者工具，或 Chrome DevTools Protocol 的 `Emulation.setDeviceMetricsOverride`）：视口 390×844，检查 `document.documentElement.scrollWidth <= innerWidth`，检查几何/膜系标签内按钮的 `getBoundingClientRect().height >= 44`。2026-10-15 的实测结果见 `/root/autodl-tmp/wt-notes/server-merge.md`；高级分析折叠区的同一检查（390×844 + 触控模拟：`scrollWidth` 390 = `innerWidth`、面板内所有元素无横向溢出、3 个输入框与 2 个按钮均 44 px）见 `/root/autodl-tmp/wt-notes/wiring.md`。

**透镜组（3D）面板的检查状态**（如实记录）：本次接线在无浏览器的构建环境里完成，GUI 只做了**静态审查 + 端点自测**：`node --check server/web/oz.js`（`cmd/oz` 的 `TestGUIScriptsParse` 在没有 node 时会 skip，本机把 `OZ_NODE` 指向内置的 node 26.7 可跑）、`TestGroupExamplesInGUI`（示例 JSON 解码）、`TestServeEndpoints`/`TestGUIServedAtRoot`（页面/脚本/样式的结构断言与 MIME）、以及用 `curl` 对新端点的逐个自测（数字见 `/root/autodl-tmp/wt-notes/wiring2.md`）。尚未做的：真机/`Emulation.setDeviceMetricsOverride` 的 390×844 实测与拖拽手势实测。

## 5. 缓存

GUI 资源返回 `Cache-Control: no-cache`：浏览器每次加载都向服务端确认，未变化时收到 304，不会把旧的 `index.html` 与新脚本混用。

## 6. 开发与检查

- 修改前端后需要重新编译二进制（资源是内嵌的）：`go build -o oz ./cmd/oz`。
- 语法检查：`node --check server/web/app.js && node --check server/web/oz.js && node --check server/web/ai.js`。`cmd/oz` 的测试 `TestGUIScriptsParse` 在有 node 时执行同样的检查（`OZ_NODE` 指定 node 路径）。
- 页面与资源的结构性断言（视口标签、分栏断点、触控断点、`clamp`、API 调用、AI 侧栏三件套 `#aiBtn`/`#aiPanel`/`ai.js`）见 `cmd/oz/serve_test.go` 与 `cmd/oz/service_test.go`。
- AI 端点的行为测试（提供方转发、SSE 中继、错误映射、规范摘要、元件落盘校验）见 `server/ai_test.go`。
