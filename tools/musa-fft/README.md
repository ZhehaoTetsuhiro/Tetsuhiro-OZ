# tools/musa-fft — WOS 的 MUSA（摩尔线程）FFT helper

## 为什么需要独立进程

MUSA 用户态驱动（libmusart）在设备初始化时会向进程地址空间预留自己的
"设备地址窗口"：[4GB, ~873GB)（两段 PROT_NONE 匿名映射）。窗口内出现任何
外来 VMA 都会导致初始化失败（musaGetDeviceCount → 101 "invalid device
ordinal"）。Go 运行时的堆 arena 固定在 0xc000000000（768GB），正好在窗口
内 —— 因此任何 Go 进程（任意 Go 版本、PIE 与否、cgo 与否）都无法初始化该
驱动；纯 C 进程没有这个问题。完整证据链见
[`../../docs/notes-musa-rootcause.md`](../../docs/notes-musa-rootcause.md)。

于是 GPU FFT 由本 helper（C 进程）承接：它持有 MUSA 上下文、muFFT 计划与
显存缓冲，通过 UNIX socket 为 Go 侧（`optics/gpu_musa.go`，`-tags musa`）
提供 2D/1D 复数 FFT。Go 侧为纯标准库、无需 cgo；helper 不存在或任何调用
失败时自动回退纯 Go CPU 路径，结果永远正确。

## 构建

    bash build.sh          # 产出 ../../../bin/wos-musa-fft 并跑一次 --check 自检

## 运行

    wos-musa-fft --socket /tmp/wos-musa-fft.sock --idle-exit 900

- Go 侧按 `$WOS_MUSA_HELPER` → 可执行文件同目录 → `$PATH` 的顺序查找，
  并自动拉起（socket 上已有守护进程则直接复用）。
- `--check`：一次性设备自检（直流分量 + 往返精度），退出码 0/1。
- `--idle-exit SECONDS`：无连接空闲 N 秒后退出（0 = 不退）。
- `--max-bytes N`：单次请求载荷上限（默认 16 GiB）。超限请求会被拒绝、该次
  变换由调用方回退 CPU（16384² 单请求即 4 GiB，32768² 为 16 GiB）。

## 协议（小端；两侧均为 x86-64）

    请求 24 字节头：magic "WUFF"(u32) | op(u32) | dir(u32) | n(u32) | bytes(u64)
         op: 0=info（回设备名，无载荷）/ 1=fft2d / 2=fft1d
         dir: 0=正变换 / 1=逆变换（未归一化，调用方补 1/n 或 1/n²）
         bytes: 载荷字节数 = n*n*16（2D）或 n*16（1D），随后跟载荷
    应答 16 字节头：status(u32, 0=ok) | rsvd(u32) | bytes(u64)，随后跟载荷

## 测试

    python3 test_client.py [socket_path]

对照 cmath 参考 DFT 验证 1D/2D 正变换与逆往返（当前实测 1e-14 / 1e-16 级）。

## 注意（MUSA 3.1.0 的坑）

muFFT 对**主机指针静默空操作**：`mufftExecZ2Z` 传主机地址会返回
`MUFFT_SUCCESS` 但不动缓冲区。必须使用设备内存（`musaMalloc` + H2D/D2H）；
本 helper 的正式路径与 `--check` 都如此。
