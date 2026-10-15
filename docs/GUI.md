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
| 几何/膜系 | `oz` | 消色差设计参数、公差参数、膜系模板 | 设计结果与处方 JSON、公差报告、膜系与光谱 |

`#oz` 哈希让页面直接进入几何/膜系模式（`/oz/` 的跳转目标）。

几何/膜系模式的操作（全部调用服务端 API，见 `docs/API.md` §4）：

1. **设计**：`POST /api/lens/design`，显示报告、优化器统计与处方 JSON；设计成功后启用「评估」「运行公差」与「高级分析」里的按钮。
2. **评估当前设计**：`POST /api/lens/eval`，显示 EFL、BFL、F/#、F–C 色差。
3. **运行公差**：`POST /api/lens/tol`，参数为半径/厚度/Δn 与运行次数；显示良率、各分布统计与失败原因。焦距判据取自设计时的目标 EFL。
4. **高级分析（折叠区）**：面板由 `oz.js` 在初始化时构建（`index.html` 不感知），复用页面既有的 `ozGrid/ozField/ozActions/ozDetails/ozOut` 类，因此手机断点与触控目标规则同样生效。「热漂移」调用 `POST /api/lens/thermal`（ΔT 与镜筒 α 取自面板输入，λ 用 d 线），显示 EFL/BFL 变化与 mm/K；「灵敏度与公差分配」调用 `POST /api/lens/sens`（`metric=merit`、`inverse=true`、预算取自面板），显示逐参数 ∂M/∂x、分配的 ±公差与封顶标记。两个按钮在首次设计成功前禁用。
5. **设计膜系**：`POST /api/coating/design`，显示膜层摘要。
6. **设计并画光谱**：先设计膜系，再 `POST /api/coating/spectrum`（400–750 nm，2 nm 步长，非偏振），在画布上绘制 R(λ)。画布按容器宽度与设备像素比重绘。

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
- 光谱画布高度 `clamp(150px, 38vw, 300px)`，宽度跟随容器，按 `devicePixelRatio` 设置后备像素。

验证方法（浏览器开发者工具，或 Chrome DevTools Protocol 的 `Emulation.setDeviceMetricsOverride`）：视口 390×844，检查 `document.documentElement.scrollWidth <= innerWidth`，检查几何/膜系标签内按钮的 `getBoundingClientRect().height >= 44`。2026-10-15 的实测结果见 `/root/autodl-tmp/wt-notes/server-merge.md`；高级分析折叠区的同一检查（390×844 + 触控模拟：`scrollWidth` 390 = `innerWidth`、面板内所有元素无横向溢出、3 个输入框与 2 个按钮均 44 px）见 `/root/autodl-tmp/wt-notes/wiring.md`。

## 5. 缓存

GUI 资源返回 `Cache-Control: no-cache`：浏览器每次加载都向服务端确认，未变化时收到 304，不会把旧的 `index.html` 与新脚本混用。

## 6. 开发与检查

- 修改前端后需要重新编译二进制（资源是内嵌的）：`go build -o oz ./cmd/oz`。
- 语法检查：`node --check server/web/app.js && node --check server/web/oz.js && node --check server/web/ai.js`。`cmd/oz` 的测试 `TestGUIScriptsParse` 在有 node 时执行同样的检查（`OZ_NODE` 指定 node 路径）。
- 页面与资源的结构性断言（视口标签、分栏断点、触控断点、`clamp`、API 调用、AI 侧栏三件套 `#aiBtn`/`#aiPanel`/`ai.js`）见 `cmd/oz/serve_test.go` 与 `cmd/oz/service_test.go`。
- AI 端点的行为测试（提供方转发、SSE 中继、错误映射、规范摘要、元件落盘校验）见 `server/ai_test.go`。
