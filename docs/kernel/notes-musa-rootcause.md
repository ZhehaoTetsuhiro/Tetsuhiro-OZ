# MUSA + Go 兼容性问题：根因与证据（2026-10-15 调试记录）

## 结论

摩尔线程 MUSA 运行时（libmusart/libsrv_um）在设备初始化时会向进程地址空间
**预留自己的"设备地址窗口"**：两段 PROT_NONE 匿名映射
[0x100000000, 0xba00000000) ≈ [4GB, 744GB) 与 [0xba00200000, 0xd9c0200000) ≈ [746GB, 873GB)。
**该窗口内出现任何外来 VMA 都会导致设备初始化失败**（musaGetDeviceCount → 101
musaErrorInvalidDevice / "invalid device ordinal"）。

Go 运行时（linux/amd64）把堆 arena 固定在 0xc000000000（768GB）——**正好落在驱动窗口内部**，
因此任何 Go 进程（cgo 或纯 Go、所有 Go 版本 1.22–1.25、PIE 与否）都无法在该驱动上初始化 MUSA。
C / Python（ctypes）进程不受影响（它们的映射在 0x55.. / 0x7f.. 高位，窗口为空）。

## 关键证据

1. C 探针 `/tmp/m2`：musaGetDeviceCount → count=1 err=0，30/30 次成功。
   Go 探针（cgo，-lmusart）：err=101，30/30 次失败。Python ctypes 成功。
2. strace 对比：两边系统调用几乎一致；区别在 mmap：
   - C: mmap(0x100000000, 794568949760, PROT_NONE, MAP_FIXED_NOREPLACE) = 成功
        mmap(0xba00200000, 136365146112, PROT_NONE, ...) = 成功
   - Go: 第一条同样成功；第二条 = -1 EEXIST（被 Go 的 0xc000000000 堆块占了）
3. 进程 maps：
   - C-musa 进程：`100000000-ba00000000 ---p` + `ba00200000-d9c01f0000 ---p` 两段窗口都在。
   - Go-musa 进程：第一段在；第二段缺失（失败）；Go 堆在 `c000000000-c000400000`。
   - 纯 Go 进程（无 musa）：无窗口段，只有 0xc000000000 堆块。
4. 参数化实验（纯 C）：在 [4GB,873GB) 内任意位置放 1GB 预留 → 立刻 101；
   4GB 以下（256MB/1GB/2GB）无影响；1TB 以上无影响。窗口边界精确吻合。
5. 排除项：GOMEMLIMIT、ulimit -v、GOGC、GOMAXPROCS、asyncpreemptoff、PIE/非PIE、
   MUSA_VISIBLE_DEVICES、LD_LIBRARY_PATH、ENV（env -i）均不影响。

## 采用的方案（WOS MUSA 后端 v2）

独立 helper 进程架构：
- `tools/musa-fft/fft_helper.c`（C、gcc、链接 libmusart+libmufft）常驻进程，
  在 C 进程地址空间里持有 MUSA 上下文 + 持久化计划/缓冲，通过 AF_UNIX socket 提供
  fft2d/fft1d(fwd/inv) RPC。C 进程默认布局里驱动窗口为空 → 工作正常。
- Go 侧（optics/gpu_musa.go）：优先尝试进程内直连（未来驱动修了就能用）；
  否则自动 fork helper 并用 socket 调用；失败一律回退纯 Go CPU 路径（语义不变）。
- 构建：tools/musa-fft/build.sh；启动脚本经 WOS_MUSA_HELPER 指定 helper 路径。

## 复现脚本（证据保留）

- /tmp/m2 (C 探针源 /tmp/m2.c)、/tmp/mres (+/tmp/mres.c 参数化预留)、
  /root/autodl-tmp/tmp-probe/probe2 (cgo 探针)、probe5 (pthread 变体)、
  /tmp/minit (+minit.c) 等。strace/maps 痕迹在 /tmp/tr-*.txt /tmp/g2-*.txt /tmp/full-*.txt。

## 附加发现（MUSA 3.1.0 / muFFT）

- muFFT 对**主机指针静默空操作**：mufftExecZ2Z 传主机地址返回
  MUFFT_SUCCESS 但不动缓冲区（输出保持原值）。必须使用设备内存
  （musaMalloc + H2D/D2H）。调试初期曾因此出现"执行成功但结果为零"的假象。
- MUSA clang 不接受 `-Od3` 编译标志（llama.cpp 上游 MUSA 编译选项里有，
  本机需改为 `-O3`）——与 WOS 无关，但同属该 SDK 的兼容性坑，记录备查。
