#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""波动光学示例库校验器（仅标准库）。

    python3 scripts/check_wave_library.py                    # 校验全部示例
    python3 scripts/check_wave_library.py --base http://127.0.0.1:6008
    python3 scripts/check_wave_library.py --simulate 48      # 额外抽样实跑并写 REPORT.md
    python3 scripts/check_wave_library.py --json             # 机读输出

做三件事：

1. **校验**：遍历 ``library/wave/*.json``，逐个 ``POST /api/validate``，要求 ``ok=true``；
   子目录 ``library/wave/quantum/*.json`` 用的是量子协议，改为 ``POST /api/quantum``
   并要求 200 + 含 ``mean_photons``。任何失败都会列出文件与 issues。
2. **统计**：按 ``meta.zod.category`` 统计数量，并检查各类别的最低数量要求。
3. **抽样实跑**（``--simulate N``）：从**网格 ≤512** 的示例中按类别均匀抽取 N 个
   （优先 ``meta.sample=true``），真实 ``POST /api/simulate`` 并轮询到 ``done``，
   记录耗时/告警/平面统计与解析预期对比，写入 ``library/wave/REPORT.md``。

退出码：0 = 全部通过（抽样运行无 error），1 = 有失败。
"""

import argparse
import json
import math
import os
import random
import sys
import time
import urllib.error
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
WAVE_DIR = os.path.join(ROOT, "library", "wave")
INDEX = os.path.join(WAVE_DIR, "index.json")
REPORT = os.path.join(WAVE_DIR, "REPORT.md")

DEFAULT_BASE = "http://127.0.0.1:6010"      # WOS CPU 实例（备用：6008 MUSA 实例）
FALLBACK_BASES = ("http://127.0.0.1:6008", "http://127.0.0.1:1120")

# 各类别的最低数量要求（硬指标）
MIN_BY_CATEGORY = {
    "干涉": 60,
    "衍射": 100,
    "传播/光束": 70,
    "偏振": 50,
    "介质": 40,
    "量子": 40,
    "仪器": 60,
    "参数扫掠": 100,
}
MIN_TOTAL_WAVE = 520          # library/wave/*.json（不含子目录）的数量下限
SIM_GRID_MAX = 512            # 抽样实跑的最大网格（避免占满服务）


# --------------------------------------------------------------------------
# HTTP 小工具（标准库）
# --------------------------------------------------------------------------

def post_json(base, path, obj, timeout=120):
    data = json.dumps(obj).encode("utf-8")
    req = urllib.request.Request(base + path, data=data,
                                 headers={"Content-Type": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return r.status, json.loads(r.read().decode("utf-8"))
    except urllib.error.HTTPError as e:
        body = e.read().decode("utf-8", "replace")
        try:
            return e.code, json.loads(body)
        except Exception:
            return e.code, {"error": body[:400]}


def get_json(base, path, timeout=120):
    try:
        with urllib.request.urlopen(base + path, timeout=timeout) as r:
            return r.status, json.loads(r.read().decode("utf-8"))
    except urllib.error.HTTPError as e:
        return e.code, {"error": e.read().decode("utf-8", "replace")[:200]}


def reachable(base):
    try:
        st, j = get_json(base, "/api/health", timeout=8)
        return st == 200 and j.get("status") == "ok", j
    except Exception as e:
        return False, {"error": str(e)}


def pick_base(explicit):
    """确定可用的服务地址：显式指定 > 默认 > 备用。"""
    cands = [explicit] if explicit else [DEFAULT_BASE] + list(FALLBACK_BASES)
    for b in cands:
        if not b:
            continue
        ok, info = reachable(b)
        if ok:
            return b, info
    return None, None


# --------------------------------------------------------------------------
# 收集示例
# --------------------------------------------------------------------------

def collect(wave_dir=WAVE_DIR):
    wave, quantum = [], []
    for name in sorted(os.listdir(wave_dir)):
        p = os.path.join(wave_dir, name)
        if not os.path.isfile(p) or not name.endswith(".json"):
            continue
        if name == "index.json":
            continue
        wave.append((name, p))
    qdir = os.path.join(wave_dir, "quantum")
    if os.path.isdir(qdir):
        for name in sorted(os.listdir(qdir)):
            if name.endswith(".json"):
                quantum.append(("quantum/" + name, os.path.join(qdir, name)))
    return wave, quantum


def load(p):
    with open(p, "r", encoding="utf-8") as fh:
        return json.load(fh)


def meta_of(cfg, fname):
    m = cfg.get("meta") or {}
    zod = m.get("zod") or {}
    return {
        "name": m.get("name") or cfg.get("name") or fname,
        "category": zod.get("category") or "（未分类）",
        "difficulty": zod.get("difficulty") or "-",
        "description": zod.get("description") or "",
        "expected": zod.get("expected") or "",
        "tags": zod.get("tags") or [],
        "sample": bool(m.get("sample")),
    }


# --------------------------------------------------------------------------
# 1) 校验
# --------------------------------------------------------------------------

def validate_all(base, wave, quantum, quiet=False):
    failures = []
    for name, path in wave:
        cfg = load(path)
        st, res = post_json(base, "/api/validate", cfg, timeout=60)
        ok = bool(res.get("ok")) if isinstance(res, dict) else False
        if st != 200 or not ok:
            failures.append((name, st, res.get("issues") if isinstance(res, dict) else res))
    qfail = []
    for name, path in quantum:
        cfg = load(path)
        st, res = post_json(base, "/api/quantum", cfg, timeout=60)
        if st != 200 or not isinstance(res, dict) or "mean_photons" not in res:
            qfail.append((name, st, str(res)[:200]))
    return failures, qfail


def selftest(base):
    """反向自检：故意构造一个坏配置，确认服务确实会判它不合法。"""
    bad = {"grid": {"size": 255, "width": 0.01}, "wavelength": 0.0,
           "sources": [{"id": "s", "type": "nonexistent", "params": {}}],
           "scene": {"components": []}}
    st, res = post_json(base, "/api/validate", bad, timeout=30)
    return (st == 200 and isinstance(res, dict) and res.get("ok") is False
            and len(res.get("issues") or []) > 0)


# --------------------------------------------------------------------------
# 2) 统计
# --------------------------------------------------------------------------

def counts_of(wave, quantum):
    c = {}
    for name, path in wave + quantum:
        m = meta_of(load(path), name)
        c[m["category"]] = c.get(m["category"], 0) + 1
    return c


# --------------------------------------------------------------------------
# 3) 抽样实跑
# --------------------------------------------------------------------------

def pick_samples(wave, n, seed=20261015):
    """按类别均匀抽样：优先 meta.sample，网格必须 ≤512。"""
    by_cat = {}
    for name, path in wave:
        cfg = load(path)
        size = (cfg.get("grid") or {}).get("size", 1 << 30)
        if size > SIM_GRID_MAX:
            continue
        m = meta_of(cfg, name)
        by_cat.setdefault(m["category"], []).append((name, path, m, size))
    rnd = random.Random(seed)
    chosen = []
    cats = sorted(by_cat, key=lambda k: -len(by_cat[k]))
    # 先每类取优先项，再轮流补足
    for cat in cats:
        items = [it for it in by_cat[cat] if it[2]["sample"]] or by_cat[cat]
        rnd.shuffle(items)
        chosen.append(items[0])
    pool = [it for cat in cats for it in by_cat[cat]
            if it not in chosen]
    rnd.shuffle(pool)
    i = 0
    while len(chosen) < n and i < len(pool):
        chosen.append(pool[i])
        i += 1
    return chosen[:n]


def run_one(base, path, timeout=300):
    cfg = load(path)
    t0 = time.time()
    st, res = post_json(base, "/api/simulate", cfg, timeout=60)
    if st != 202 or not isinstance(res, dict) or "run_id" not in res:
        return {"status": "submit_failed", "detail": str(res)[:200], "elapsed": 0.0}
    rid = res["run_id"]
    while time.time() - t0 < timeout:
        st2, j = get_json(base, "/api/runs/" + rid, timeout=60)
        if st2 != 200 or not isinstance(j, dict):
            return {"status": "poll_failed", "detail": str(j)[:200],
                    "elapsed": time.time() - t0}
        if j.get("status") != "running":
            j["elapsed"] = time.time() - t0
            if j.get("status") == "done" and not (j.get("planes") or []):
                j["status"] = "no_planes"      # 探测器没收到光（几何没对上）
            return j
        time.sleep(0.05)
    return {"status": "timeout", "elapsed": time.time() - t0}


def run_samples(base, wave, n, timeout=300):
    rows = []
    for name, path, m, size in pick_samples(wave, n):
        cfg = load(path)
        j = run_one(base, path, timeout=timeout)
        planes = []
        for p in (j.get("planes") or []):
            st = p.get("stats") or {}
            planes.append({"id": p.get("id"), "power": st.get("power"),
                           "peak": st.get("peak"), "strehl": st.get("strehl"),
                           "rms_x": st.get("rms_x"), "label": p.get("label")})
        rows.append({
            "file": name, "name": m["name"], "category": m["category"],
            "difficulty": m["difficulty"], "grid": size,
            "status": j.get("status"), "elapsed_ms": j.get("elapsed_ms"),
            "wall_s": round(j.get("elapsed", 0.0), 2),
            "warnings": [w.get("code") for w in (j.get("warnings") or [])],
            "planes": planes, "detail": j.get("detail"),
            "expected": m["expected"],
        })
        print("  %-38s %-8s %6.2fs  %s" % (name[:38], rows[-1]["status"],
                                           rows[-1]["wall_s"], rows[-1]["warnings"]))
    return rows


def run_quantum_samples(base, quantum, n=10):
    rows = []
    step = max(1, len(quantum) // n)
    for name, path in quantum[::step][:n]:
        cfg = load(path)
        st, res = post_json(base, "/api/quantum", cfg, timeout=60)
        ok = st == 200 and isinstance(res, dict) and "mean_photons" in res
        rows.append({"file": name, "status": "done" if ok else "error",
                     "mean_photons": res.get("mean_photons") if ok else None,
                     "g2": res.get("g2") if ok else None,
                     "norm": res.get("norm") if ok else None,
                     "postselect_probability": res.get("postselect_probability") if ok else None})
    return rows


# --------------------------------------------------------------------------
# 报告
# --------------------------------------------------------------------------

def write_report(path, base, counts, wave_n, quantum_n, failures, qfail,
                 sim_rows, q_rows, min_checks):
    L = []
    L.append("# 波动光学示例库实跑报告（library/wave/REPORT.md）\n")
    L.append("由 `scripts/check_wave_library.py --simulate N` 自动生成。\n")
    L.append("- 服务地址：`%s`\n" % base)
    L.append("- 生成时间：%s\n" % time.strftime("%Y-%m-%d %H:%M:%S"))
    L.append("- 规模：`library/wave/*.json` **%d** 个（可直接 POST /api/simulate），"
             "`library/wave/quantum/*.json` **%d** 个（POST /api/quantum）\n" % (wave_n, quantum_n))
    L.append("- 校验：波动示例 `/api/validate` 全部 `ok=true` 失败 %d 个；"
             "量子示例 `/api/quantum` 返回正常失败 %d 个\n" % (len(failures), len(qfail)))

    L.append("\n## 1. 分类计数与硬指标\n\n")
    L.append("| 分类 | 数量 | 要求 | 结论 |\n|---|---|---|---|\n")
    for cat, req in min_checks:
        have = counts.get(cat, 0)
        L.append("| %s | %d | ≥%d | %s |\n" % (cat, have, req, "✅" if have >= req else "❌"))
    L.append("| **合计** | **%d** | — | — |\n" % sum(counts.values()))

    L.append("\n## 2. 抽样实跑（真实 /api/simulate，网格 ≤%d）\n\n" % SIM_GRID_MAX)
    ok = sum(1 for r in sim_rows if r["status"] == "done")
    L.append("共实跑 **%d** 个示例，`done` **%d** 个。耗时含排队与轮询。\n\n" % (len(sim_rows), ok))
    L.append("| 文件 | 分类 | 难度 | 网格 | 状态 | 服务端耗时 | 墙钟 | 告警 | 输出平面（功率 W / 峰值 W·m⁻²） |\n")
    L.append("|---|---|---|---|---|---|---|---|---|\n")
    for r in sim_rows:
        pl = "；".join("%s: P=%.3g peak=%.3g" % (p["id"], p["power"] or 0.0, p["peak"] or 0.0)
                       for p in r["planes"]) or "—"
        L.append("| `%s` | %s | %s | %d | %s | %s | %.2fs | %s | %s |\n"
                 % (r["file"], r["category"], r["difficulty"], r["grid"], r["status"],
                    ("%.0f ms" % r["elapsed_ms"]) if r["elapsed_ms"] else "—",
                    r["wall_s"], ",".join(r["warnings"]) or "—", pl))

    L.append("\n### 2.1 各示例的解析预期（抽样）\n\n")
    for r in sim_rows:
        L.append("**`%s`** — %s（%s，%s）\n\n> %s\n\n"
                 % (r["file"], r["name"], r["category"], r["difficulty"], r["expected"]))

    L.append("\n## 3. 量子示例抽样（/api/quantum，同步返回）\n\n")
    L.append("| 文件 | 状态 | ⟨n⟩ | g²(0) | norm | 后选择概率 |\n|---|---|---|---|---|---|\n")
    for r in q_rows:
        mn = ", ".join("%.3f" % x for x in (r["mean_photons"] or [])) or "—"
        g2 = ", ".join("%.3f" % x for x in (r["g2"] or [])) or "—"
        pp = r["postselect_probability"]
        L.append("| `%s` | %s | %s | %s | %s | %s |\n"
                 % (r["file"], r["status"], mn, g2,
                    ("%.4f" % r["norm"]) if r["norm"] is not None else "—",
                    ("%.4f" % pp) if pp is not None else "—"))

    L.append("\n## 4. 已知限制（本内核对这些现象是近似表示）\n\n")
    L.append("""1. **多光束（Airy）条纹**：定位场景只计入**首次往返**，高反射率平行平板/法布里-珀罗腔的
   多光束窄峰不在结果里（会给出 `scene_cycle_dropped` 告警，见 `inst_fabry_perot_*`、
   `interf_parallel_plate_*`）。文中给出的精细度 F、自由光谱范围是解析对照值。
2. **牛顿环**用「平面镜 + 曲率镜」的两光束模型实现，等价于厚度 t=r²/2R 的空气楔相位 kr²/R；
   真实的玻璃-空气-玻璃三层反射（含半波损失与多反射）不在模型中，中心的亮/暗取决于附加 π 相移。
3. **萨格纳克环**：scene 模型没有旋转参考系，因此**不能用旋转产生非互易相位**；
   示例演示的是「互易环输出恒为零」（这是真实物理），旋转灵敏度只能靠环内非互易元件间接展示。
4. **量子示例**走 `/api/quantum`（Fock 基线性光学），与波动光学内核相互独立；
   其 "config" 不含 grid/wavelength，不能用 `/api/simulate` 运行。
5. **标量薄元件**：元件模型是逐点复透过率（+可选琼斯矩阵），因此**忽略多次反射、
   体散射与真正的 3D 曲面折射**；单轴/双轴晶体走近似传播（Berreman 4×4 用于双轴），
   大角度斜入射与界面反射（菲涅耳系数）不建模。
6. **数值上限**：网格分辨率 dx=W/N 决定能表示的最小结构（本库生成时强制「特征尺寸 ≥4 像素」），
   亚像素结构会被静默平均掉；远场大角度衍射需要 `asm_pad`/`asm_shift`（见 PHYSICS.md §3）。
7. **超表面/波带片**：`metalens_native` 与 `zone_plate` 是**单波长写死的相位掩模**，
   因此色散（焦点移到 f·λ₀/λ）是物理正确的，但**没有**真实的亚波长结构散射与偏振转换效率。
8. **介质元件**为薄片近似（相位 k·n·d + 吸收 exp(−2k n'' d)），配合 split-step 表达
   增益/克尔/可饱和吸收；不包含 Δn 的横向扩散、热效应与增益饱和。
""")
    L.append("\n## 5. 复现方式\n\n")
    L.append("```bash\npython3 library/tools/build_wave_library.py --out library/wave   # 重新生成全部示例\n"
             "python3 scripts/check_wave_library.py                              # 校验全部示例\n"
             "python3 scripts/check_wave_library.py --simulate 48                # 抽样实跑 + 刷新本报告\n"
             "```\n")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("".join(L))


# --------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default=None, help="服务地址（默认自动探测 %s）" % DEFAULT_BASE)
    ap.add_argument("--simulate", type=int, default=0,
                    help="抽样实跑的示例数（≥40 才满足交付要求；0 = 不实跑）")
    ap.add_argument("--quantum-samples", type=int, default=10)
    ap.add_argument("--json", action="store_true", help="只输出机读 JSON")
    ap.add_argument("--timeout", type=int, default=300, help="单个示例的最长等待秒数")
    ap.add_argument("--seed", type=int, default=20261015)
    a = ap.parse_args()

    base, health = pick_base(a.base)
    if base is None:
        print("× 找不到可用的 WOS 服务（试过 %s）" % (a.base or DEFAULT_BASE), file=sys.stderr)
        return 2

    wave, quantum = collect()
    if not wave and not quantum:
        print("× 在 %s 下没有找到示例" % WAVE_DIR, file=sys.stderr)
        return 2

    ok_self = selftest(base)
    failures, qfail = validate_all(base, wave, quantum)
    counts = counts_of(wave, quantum)

    min_checks = sorted(MIN_BY_CATEGORY.items())
    cat_ok = all(counts.get(c, 0) >= r for c, r in min_checks)

    sim_rows, q_rows = [], []
    if a.simulate:
        sim_rows = run_samples(base, wave, a.simulate, timeout=a.timeout)
        q_rows = run_quantum_samples(base, quantum, a.quantum_samples)

    summary = {
        "base": base,
        "service": health,
        "validator_selftest_detects_bad_config": ok_self,
        "wave_examples": len(wave),
        "quantum_examples": len(quantum),
        "validate_failures": len(failures),
        "quantum_failures": len(qfail),
        "by_category": counts,
        "category_minimums_ok": cat_ok,
        "simulated": len(sim_rows),
        "simulated_done": sum(1 for r in sim_rows if r["status"] == "done"),
        "index_exists": os.path.isfile(INDEX),
    }
    if a.json:
        print(json.dumps(summary, ensure_ascii=False, indent=2))
    else:
        print("服务地址        : %s（%s）" % (base, (health or {}).get("gpu", {}).get("backend", "?")))
        print("校验器自检      : %s" % ("通过（能识别非法配置）" if ok_self else "失败！"))
        print("波动示例（顶层）: %d 个；量子示例：%d 个" % (len(wave), len(quantum)))
        print("validate 失败   : %d 个" % len(failures))
        for n, st, iss in failures[:20]:
            print("   - %s (HTTP %s) %s" % (n, st, iss))
        print("quantum 失败    : %d 个" % len(qfail))
        for n, st, iss in qfail[:20]:
            print("   - %s (HTTP %s) %s" % (n, st, iss))
        print("分类计数        :")
        for c, r in min_checks:
            have = counts.get(c, 0)
            print("   %-10s %4d / ≥%d  %s" % (c, have, r, "OK" if have >= r else "不足"))
        extra = {k: v for k, v in counts.items() if k not in dict(min_checks)}
        if extra:
            print("   其它类别:", extra)
        print("合计            : %d" % sum(counts.values()))
        if sim_rows:
            print("抽样实跑        : %d 个，done %d 个"
                  % (len(sim_rows), sum(1 for r in sim_rows if r["status"] == "done")))
        if a.simulate:
            write_report(REPORT, base, counts, len(wave), len(quantum), failures, qfail,
                         sim_rows, q_rows, min_checks)
            print("报告已写入      : %s" % REPORT)

    bad = len(failures) + len(qfail) + (0 if ok_self else 1)
    if not cat_ok or len(wave) < MIN_TOTAL_WAVE:
        bad += 1
    if sim_rows and any(r["status"] != "done" for r in sim_rows):
        bad += 1
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
