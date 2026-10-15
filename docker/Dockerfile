# Tetsuhiro-OZ 容器镜像
#
# OZ 是纯 Go（CGO_ENABLED=0）静态单文件，故镜像只需要二进制 + 元件定义文件。
# 本地手工构建：
#
#   export CGO_ENABLED=0 GOEXPERIMENT=nodwarf5 GOFLAGS= GOTOOLCHAIN=go1.25.10
#   go build -trimpath -ldflags="-s -w" -o dist/oz ./cmd/oz
#   docker build -f docker/Dockerfile --build-arg VERSION=0.1.0-Preview \
#     --build-arg REVISION=$(git rev-parse HEAD) -t tetsuhiro-oz:0.1.0-Preview .
#
# 运行（默认监听 :1130；如需兼容旧脚本可 -addr :1120）：
#
#   docker run -p 1130:1130 tetsuhiro-oz:0.1.0-Preview
#
FROM alpine:3.20

ARG VERSION=dev
ARG REVISION=unknown

LABEL org.opencontainers.image.title="Tetsuhiro-OZ" \
      org.opencontainers.image.description="完整光学模拟器：几何光学 + 波动/量子光学（纯 Go 静态单文件）+ Web GUI" \
      org.opencontainers.image.source="https://github.com/ZhehaoTetsuhiro/Tetsuhiro-OZ" \
      org.opencontainers.image.url="https://github.com/ZhehaoTetsuhiro/Tetsuhiro-OZ/releases" \
      org.opencontainers.image.licenses="AGPL-3.0-or-later" \
      org.opencontainers.image.version="${VERSION}" \
      org.opencontainers.image.revision="${REVISION}"

RUN addgroup -S oz && adduser -S -G oz -h /home/oz oz

COPY dist/oz /usr/local/bin/oz
COPY elements /usr/share/oz/elements

USER oz
WORKDIR /home/oz
EXPOSE 1130

ENTRYPOINT ["/usr/local/bin/oz", "serve", "-addr", ":1130", "-elements", "/usr/share/oz/elements"]
