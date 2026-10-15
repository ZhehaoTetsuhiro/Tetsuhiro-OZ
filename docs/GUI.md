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

## 4. 移动端

页面面向手机与平板：

- `<meta name="viewport" content="width=device-width, initial-scale=1">`。
- 宽度 ≤ 700 px：页头按钮换行；左右两栏堆叠为单列，页面整体纵向滚动，不出现横向滚动条。
- 触控设备（`pointer: coarse`）与几何/膜系标签内的按钮、输入框最小高度 44 px。
- 光谱画布高度 `clamp(150px, 38vw, 300px)`，宽度跟随容器，按 `devicePixelRatio` 设置后备像素。

验证方法（浏览器开发者工具，或 Chrome DevTools Protocol 的 `Emulation.setDeviceMetricsOverride`）：视口 390×844，检查 `document.documentElement.scrollWidth <= innerWidth`，检查几何/膜系标签内按钮的 `getBoundingClientRect().height >= 44`。2026-10-15 的实测结果见 `/root/autodl-tmp/wt-notes/server-merge.md`。

## 5. 缓存

GUI 资源返回 `Cache-Control: no-cache`：浏览器每次加载都向服务端确认，未变化时收到 304，不会把旧的 `index.html` 与新脚本混用。

## 6. 开发与检查

- 修改前端后需要重新编译二进制（资源是内嵌的）：`go build -o oz ./cmd/oz`。
- 语法检查：`node --check server/web/app.js && node --check server/web/oz.js`。`cmd/oz` 的测试 `TestGUIScriptsParse` 在有 node 时执行同样的检查（`OZ_NODE` 指定 node 路径）。
- 页面与资源的结构性断言（视口标签、分栏断点、触控断点、`clamp`、API 调用）见 `cmd/oz/serve_test.go` 与 `cmd/oz/service_test.go`。
