# Tetsuhiro-WOS 文档存档

本目录是 Tetsuhiro-WOS 的原始文档，随内核一并内化进 Tetsuhiro-OZ
（2026-10-15，`archive-pre-oz-absorb`）。注意：

- 文中出现的 `cmd/wos`、端口 `:1120`、`import "twos/optics"` 等均为
  **WOS 时代**的写法；在 OZ 中对应 `cmd/oz`、默认端口 `:1130`
  （部署实例用 `:6006`/`:6008`，见 `/root/autodl-tmp/README-services.md`）
  与 `import "tetsuhiro-oz/optics"`。
- 现行 API 与 GUI 说明已更新在上一级目录：`../API.md`、`../GUI.md`。
- 内核代码（`optics/`）本身未重写，`go test ./optics/` 与
  `go test -tags musa ./optics/` 全绿。
