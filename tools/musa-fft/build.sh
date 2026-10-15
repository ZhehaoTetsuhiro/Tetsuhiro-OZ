#!/bin/bash
# ============================================================================
# build.sh — 编译 WOS 的 MUSA FFT helper（纯 C 进程，承接 GPU FFT）
#
# 用法：
#   bash build.sh [输出路径]        # 默认输出到 ../../../bin/wos-musa-fft
#
# 背景：MUSA 驱动要求设备地址窗口 [4GB,~873GB) 无外来映射，而 Go 运行时
# 的堆恰好占用其中 768GB 处 —— 任何 Go 进程都无法初始化该驱动；纯 C 进程
# 正常。因此 GPU FFT 由本 helper（C 进程）经 UNIX socket 提供给 Go 侧。
# 详情见仓库 docs/notes-musa-rootcause.md。
# ============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
OUT="${1:-$HERE/../../../bin/wos-musa-fft}"
MUSA_DIR="${MUSA_DIR:-/usr/local/musa}"

gcc -O2 -Wall -o "$OUT" "$HERE/fft_helper.c" \
  -I"$MUSA_DIR/include" -L"$MUSA_DIR/lib" \
  -lmufft -lmusart -lmusa -lpthread \
  -Wl,-rpath,"$MUSA_DIR/lib"

echo "built: $OUT"
"$OUT" --check
