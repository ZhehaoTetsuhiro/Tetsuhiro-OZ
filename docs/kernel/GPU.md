# GPU 加速设计（docs/GPU.md，可选/远期）

本文给出内核的 GPU 加速接口设计。当前实现为纯 Go 零依赖（仅标准库）；GPU 后端是可选加速，
必须保持 CPU 路径不变、通过构建标签（build tag）按需启用，并以 CPU 结果为基准验证一致性。

## 1. 加速目标

按热点排序（见 README 性能表，1024² 单次 ASM ~150ms）：

1. 2-D/1-D FFT（fft2D / fft1DAny，占 ASM 绝大部分耗时）。
2. 传递函数逐像素相乘（ctx.transfer、propASMPadCore 的频域循环）。
3. split-step 的介质相位屏（applyMediumPhase，逐像素复指数）。

## 2. 后端选择

- CUDA（cuFFT）经 cgo：吞吐最高，NVIDIA 硬件；用 //go:build cuda 标签隔离。
- OpenCL：跨厂商（NVIDIA/AMD/Intel/Apple），可动态发现设备。
- 混合精度：FFT 用 complex64 可再降一半带宽；相位敏感路径保留 complex128（或做误差回归测试）。

## 3. 接口抽象

在 fft.go 引入后端分派：

    // 内部：fft1DAny/fft2D 在 gpu 标签下改走 fftBackend.Do(...)
    type fftBackend interface {
        Plan1D(n int) FFTPlan
        Plan2D(n int) FFTPlan
        Exec(p FFTPlan, a []complex128, inverse bool)
    }

传播/介质层不感知后端：它们只调用 fft2D/fft1DAny 与逐像素循环；GPU 版把逐像素循环也下推为内核
（transfer、applyMediumPhase 的核函数化）。

## 4. 内存与数据搬运

- 1024² 复 128 场 ~16 MB/分量；PCIe 传输 ~ms 级，与单次 FFT 相当，需尽量“驻留显存”。
- 建议把整条光路的平面数组整体搬入显存，多平面传播在 GPU 内完成，仅最终传感器回传。
- server 的 LRU 预算（-max-run-mb）需计入显存池。

## 5. 正确性保障

- 复用现有精度测试作为“CPU 参照”：同一输入，GPU 与 CPU 相对 L2 误差须 < 1e-6（单精度则 < 1e-4）。
- Parseval / 冲激谱 / 往返三项 FFT 测试在 GPU 后端上重放。
- 非 2 的幂（Bluestein）在 GPU 上应退化为 CPU 或由 cuFFT 的任意长度支持，需单独验证。

## 6. 状态

当前仅落地本设计文档；未引入 cgo/CUDA 依赖，以维持零依赖单二进制分发。

## 7. 实现状态（已落地，2026-10）

第 1、3 节的设计已实现于 `optics/gpu.go`（`//go:build cuda`）与 `optics/gpu_stub.go`（默认构建）：

- 分派点：`fft2D` / `fft1DAny` 在启用时调用 `gpuFFT2D` / `gpuFFT1D`；纯 Go 实现保留为
  `fft2DCPU` / `fft1DAnyCPU`，任何 CUDA 调用失败即回退，所以无 GPU 机器结果不变。
- cuFFT Z2Z（double）单次 2-D 变换替代逐行/逐列；任意长度（含非 2 的幂）由 cuFFT 原生支持。
- 显存缓冲按字节大小做自由表复用（`gpuAlloc`/`gpuRelease`），避免每次 `cudaMalloc/cudaFree`。
- 线程安全：所有显存拷贝与 cuFFT 执行经单一互斥量串行化（句柄并发执行无保证）。
- 启用：`-gpu` 命令行开关或 `WOS_GPU=1`；`optics.SetGPU/GPUEnabled/GPUAvailable/GPUBackendInfo`。
- 阈值：`gpuMin2D=64`、`gpuMin1D=1024`（小于此的网格传输开销大于收益）。

实测（Tesla T4，`-bench FFT2D`，单次含往返）：1024² CPU 14 ms / GPU 12 ms（1.1×）；
2048² 57 ms / 45 ms（1.3×）；4096² 331 ms / 179 ms（1.85×）。

未做（后续可加）：真正“驻留显存”的多平面光路（把 `ctx.transfer`/`applyMediumPhase` 的逐像素循环
下推为 CUDA 内核，目标是每次传播只做一次 H2D/D2H）；本条需引入 nvcc 编译的 `.cu` 内核，
而当前后端只依赖运行库（cudart/cuFFT），保持 `go build -tags cuda` 无需 nvcc。

实测补充：页面锁定（pinned）暂存内存在本机反而更慢（4.1 GB/s vs 可分页 6.8 GB/s），
故未启用 pinned 路径；`gpuHostAlloc/gpuCopyH2D/gpuCopyD2H` 保留给基准与将来设备使用。

## 8. MUSA（摩尔线程）后端：helper 进程架构（已落地，2026-10）

`go build -tags musa` 启用摩尔线程 GPU（MUSA/muFFT）FFT 加速。与 CUDA 版
（进程内 cgo）不同，MUSA 版把 FFT 交给一个独立的 C helper 进程
（`tools/musa-fft/`，产物 `wos-musa-fft`），Go 侧经 UNIX socket 调用。

**为什么必须用独立进程**：MUSA 用户态驱动在初始化时向进程地址空间预留自己的
设备地址窗口 [4GB, ~873GB)，窗口内出现任何外来 VMA 都会让
`musaGetDeviceCount` 返回 101（"invalid device ordinal"）。Go 运行时的堆
arena 固定在 0xc000000000（768GB）——正好落在窗口内，因此**任何 Go 进程
（任意 Go 版本、PIE 与否、cgo 与否、ulimit/环境如何）都无法初始化该驱动**；
纯 C 进程不受影响。完整证据链（strace/maps/参数化实验）见
[notes-musa-rootcause.md](notes-musa-rootcause.md)。

架构与语义（与 CUDA 后端逐项对齐）：

- 分派点与回退完全一致：`fft2D`/`fft1DAny` 出错即回退纯 Go；
  `-gpu` 开关 / `WOS_GPU=1`、`SetGPU/GPUEnabled/GPUBackendInfo` API 相同；
  阈值沿用 `gpuMin2D=64`/`gpuMin1D=1024`。
- helper 发现顺序：`$WOS_MUSA_HELPER` → 可执行文件同目录 → `$PATH`；
  socket 默认 `/tmp/wos-musa-fft.sock`（`$WOS_MUSA_SOCKET` 可改），已有
  守护进程直接复用，否则惰性拉起（`--idle-exit 900` 自动回收）。
- Go 侧为**纯 Go（无 cgo）**：`-tags musa` 构建不需要 MUSA SDK；只有运行
  GPU 路径才需要 helper 二进制与设备。
- GPU 与 CPU 的一致性由 `go test -tags musa ./optics/ -run TestGPU` 验证
  （2-D/1-D 往返、asm/asm_pad 整步传播，判据 < 1e-9，实测通过）。

构建与运行：

    bash tools/musa-fft/build.sh           # gcc 编译 helper 并 --check 自检
    go build -tags musa -o wos ./cmd/wos   # 不需要 cgo / MUSA SDK
    WOS_GPU=1 ./wos -addr :1120            # 启用 GPU（自动发现/拉起 helper）

实测（MTT S4000 48GB、驱动 2.7.0、MUSA SDK 3.1.0，同机 15 核 CPU 对照；
helper 架构含 socket + PCIe 双向搬运，为端到端值）：

| 网格 | CPU（纯 Go，15 核配额） | GPU（helper 端到端） |
|---|---|---|
| 1024² | ~5 ms | ~16 ms |
| 2048² | ~26–30 ms | ~56–80 ms |
| 4096² | ~220–470 ms* | ~230–320 ms* |

\* 并发负载下波动较大，只作量级参考。结论（诚实版）：该驱动/该库下 helper 的
搬运开销在小网格上大于 GPU 收益，4096² 与 15 核 CPU 打平；GPU 路径的价值在于
正确可用、接口/语义与 CUDA 版一致、以及 CPU 紧张或超大网格场景。muFFT 本体
很快（512² 每变换 daemon 侧 ~2–3 ms 含搬运）。

已知坑（MUSA 3.1.0）：muFFT 对**主机指针静默空操作**——`mufftExecZ2Z` 传主机
地址会返回 `MUFFT_SUCCESS` 但不动缓冲区，必须使用设备内存（helper 的正式路径
与 `--check` 都已如此）。

已知坑二：helper 对单次请求设有载荷上限（默认 **16 GiB**，`--max-bytes` 可调）。
超过上限的请求会被拒绝，该次变换按设计静默回退 CPU——例如 16384² 网格每个请求
就是 4 GiB。默认值曾为 3 GiB，导致大网格全程静默走 CPU、GPU 使用率读数恒为 0
（排查口诀：helper 关闭 verbose 时看 /proc/&lt;pid&gt;/io 的 rchar 是否增长；verbose
时直接看 `bad request ... bytes=` 行）。

降级可视性（已内置）：`GET /api/health` 的 `gpu` 字段报告 available/enabled/
attempts/failures/degraded/last_error（degraded = 已启用且最近 10 分钟内有失败）；
服务器日志在"启用但持续回退 CPU"时每分钟最多打一条警告。内核计数由
`optics.GPUState()` 提供（`gpu_status.go`），CUDA 与 MUSA 两条后端共用。
