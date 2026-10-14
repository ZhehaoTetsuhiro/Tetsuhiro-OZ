# 本地 Qwen3.8-27B 在 OZ 开发中的使用记录

服务：`127.0.0.1:8080`（llama.cpp MUSA 版，OpenAI 兼容；幂等启动脚本
`/root/autodl-tmp/models/serve-qwen-musa.sh`；实测生成 ~5.4 t/s）。

## 使用原则（用户指令的落地方式）

- **核心 / 难度极高的模块**（geo 面型与求交、lens 优化器、coating TMM、wave 传播与 PSF）
  由人工设计实现，不外包给本地模型；
- 本地模型只承担**非核心、规格清晰、结果可被编译器与测试判定**的模块；
- **每一行产出都逐行审查**，并以 `gofmt` / `go vet` / `go test` 作为最终判据；
- 记录请求与响应原始文件，便于复现与追责。

## 用例 1：`plot` 包（SVG 图表工具，非核心）

- 请求脚本：`/root/autodl-tmp/scripts/qwen_gen_plot.py`（中文规格：纯标准库、确定性输出、
  `viewBox` 响应式、nice-number 刻度 5~8 个、图例、XML 转义、空数据/单点不 panic、
  `WriteFile` 与 `SVG()` 一致）
- 原始响应：`/root/autodl-tmp/qwen-codegen-20261015_013812.json`
  （completion 4969 tokens，耗时 1041.8 s，`truncated=0`）
- 产出：`plot/svg.go`（475 行）+ `plot/svg_test.go`（119 行，6 个测试）

### 逐行审查结论

| # | 问题 | 严重度 | 处理 |
|---|------|--------|------|
| 1 | `func niceNum(range float64, round bool)`——用 Go 关键字 `range` 当参数名，**无法编译** | 致命 | 参数改名 `span` |
| 2 | 用到 `math.Log10/Pow/Floor/Ceil` 却未导入 `math` | 致命 | 补 `"math"` 导入 |
| 3 | 导入了 `"fmt"` 但全文件未使用 | 致命 | 删除该导入 |
| 4 | 刻度值可能略微超出数据区间，网格线/刻度文字会画到坐标框外 | 轻微 | `xToPixel/yToPixel` 增加区间夹取 |
| 5 | 图例 X 固定为 `plotWidth-120`，窄图会越出左边界 | 轻微 | 增加下界保护 |
| 6 | 转义（单遍、不双重转义）、确定性（无 map 迭代）、`viewBox` + `width:100%`、单点散点、空图表、`WriteFile` 一致性 | — | 审查通过，保留原实现 |

- 实证：沙盒 `/tmp/qwenplot` 与仓库内均 `gofmt -l` 为空、`go vet` 干净、
  `go test ./plot/` 6/6 PASS。
- 落地：新包 `plot/`，并由 CLI 接线使用：`oz coating spectrum <stack.json> --svg out.svg`
  输出 R/T 光谱图（可嵌入报告或网页）。

## 结论

本地 Qwen 对"规格清晰 + 可机器判定"的非核心模块可用，但**当前水平必须逐行审查**：
首轮产出含 3 个编译级错误（关键字做参数名、缺导入、多导入），直接采用会立刻破坏构建。
"每段代码都要检查"不是形式要求，而是必需步骤。
