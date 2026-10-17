"""JS 语法检查（node --check 的替代品）：用 QuickJS 解析 GUI 脚本。

本机（以及任何没有 node 的构建环境）跑不了 `node --check`，而 `cmd/oz` 的
TestGUIScriptsParse 在没有 node 时会 skip —— 也就是说改完前端可能完全没有
语法检查。这个脚本用 pip 的 quickjs 包做同一件事：只编译、不执行，所以
document 之类不存在也不会报错。

    python3 scripts/check-js.py                       # 默认四个 GUI 脚本
    python3 scripts/check-js.py server/web/oz.js      # 指定文件

退出码：0 = 全部通过，1 = 有语法错误，2 = 没有 quickjs（改用 node --check）。
"""
import sys

try:
    import quickjs
except ImportError:  # 没有依赖不算失败：交给 node
    print("quickjs 未安装（pip install quickjs）；请改用 node --check", file=sys.stderr)
    sys.exit(2)

DEFAULT_FILES = [
    "server/web/app.js",
    "server/web/oz.js",
    "server/web/ai.js",
    "server/web/scan.js",
]


def check(path: str) -> bool:
    src = open(path, encoding="utf-8").read()
    ctx = quickjs.Context()
    try:
        # 包在函数体里编译：语法错误照样抛出，顶层 await 之类也不会误判。
        ctx.eval("(function(){\n" + src + "\n})")
        return True
    except Exception as e:  # quickjs 的 JS 异常类型随版本而变，统一兜住
        print("FAIL %s: %s" % (path, str(e)[:400]))
        return False


def main(argv):
    files = argv[1:] or DEFAULT_FILES
    ok = True
    for f in files:
        if check(f):
            print("OK   %s" % f)
        else:
            ok = False
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
