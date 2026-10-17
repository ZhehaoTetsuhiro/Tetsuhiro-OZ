/* Tetsuhiro-OZ 波动模式的「参数扫描」面板（左栏，折叠区）。
 *
 * 独立 IIFE：不向全局作用域写任何名字（app.js 的 S、oz.js 的私有 $ 互不干扰）。
 * 只做三件事：
 *   1. 从当前配置里挑出可扫描的参数（元件 roll/yaw/pitch/位置、形状与元件参数、
 *      光源参数、波长），生成一份参数下拉列表；
 *   2. 让用户填一组取值（逗号列表或 a-b-step 范围），选一个指标，POST /api/scan；
 *   3. 把返回的点画成曲线，并列出每个点的数值（含逐点错误）。
 *
 * 全部计算在服务端完成（optics.RunScan），这里只做表单、绘图与结果展示。
 */
(function () {
  "use strict";

  // $ 同时接受 "#id" 与裸 id（与 oz.js/ai.js 的写法一致，避免两种风格混用出错）。
  const $ = (id) => document.getElementById(String(id).replace(/^#/, ""));
  const DEG = 180 / Math.PI;

  // 指标：与内核 optics.planeMetric 的取值一一对应。
  const METRICS = [
    ["power", "功率 (W)"],
    ["peak", "峰值 (W/m²)"],
    ["centroid_x", "质心 x (m)"],
    ["centroid_y", "质心 y (m)"],
    ["rms_x", "RMS x (m)"],
    ["rms_y", "RMS y (m)"],
    ["strehl", "斯特列尔比"],
    ["intensity_max", "最大强度 (W/m²)"],
  ];

  // scanPaths 列出当前配置里可扫描的数值字段：路径 → 中文标签。
  // 与内核 setConfigValue 支持的路径保持一致（grid.* 不可扫：它是采样而非物理）。
  function scanPaths() {
    const out = [["wavelength", "波长 (m)"]];
    (S.config.sources || []).forEach((s, i) => {
      const name = s.label || s.id || ("光源 " + (i + 1));
      Object.keys(s.params || {}).forEach((k) => {
        if (typeof s.params[k] === "number") out.push(["sources[" + i + "].params." + k, name + " · " + k]);
      });
    });
    ((S.config.scene && S.config.scene.components) || []).forEach((c, i) => {
      const name = c.label || c.id || ("元件 " + (i + 1));
      out.push(["scene.components[" + i + "].roll", name + " · 自转 roll"]);
      out.push(["scene.components[" + i + "].yaw", name + " · 偏航 yaw"]);
      out.push(["scene.components[" + i + "].pitch", name + " · 俯仰 pitch"]);
      out.push(["scene.components[" + i + "].pos.x", name + " · 位置 x"]);
      out.push(["scene.components[" + i + "].pos.y", name + " · 位置 y"]);
      out.push(["scene.components[" + i + "].pos.z", name + " · 位置 z"]);
      Object.keys((c.shape && c.shape.params) || {}).forEach((k) => {
        if (typeof c.shape.params[k] === "number") out.push(["scene.components[" + i + "].shape.params." + k, name + " · 形状 " + k]);
      });
      Object.keys(c.params || {}).forEach((k) => {
        if (typeof c.params[k] === "number") out.push(["scene.components[" + i + "].params." + k, name + " · " + k]);
      });
    });
    return out;
  }

  // buildPanel 在左栏「全局参数」之后插入扫描面板（只建一次）。
  function buildPanel() {
    if ($("#scanPanel")) return;
    const host = $("waveGlobals");
    if (!host) return;
    const sec = document.createElement("section");
    sec.id = "scanPanel";
    sec.innerHTML =
      '<h2 style="margin-top:10px">参数扫描</h2>' +
      '<div class="prow"><label for="scanParam" title="选择一个数值参数，逐点代入后重复仿真">扫描参数</label>' +
      '<select id="scanParam"></select></div>' +
      '<div class="prow"><label for="scanValues" title="逗号列表（0,15,30）或范围（0-90-5 = 起-止-步长）">取值</label>' +
      '<input type="text" id="scanValues" value="0,15,30,45" spellcheck="false" autocapitalize="off" ' +
      'placeholder="0,15,30 或 0-90-5"></div>' +
      '<div class="prow"><label for="scanMetric">指标</label><select id="scanMetric"></select></div>' +
      '<div class="ozActions">' +
      '<button id="scanRunBtn" type="button">运行扫描</button>' +
      '<button id="scanCsvBtn" type="button" class="secondary" disabled>导出 CSV</button>' +
      '</div>' +
      '<pre class="ozOut" id="scanOut" aria-live="polite">选择参数与取值后点「运行扫描」。</pre>' +
      '<canvas id="scanPlot" aria-label="扫描曲线" height="150"></canvas>';
    host.appendChild(sec);

    $("scanRunBtn").addEventListener("click", () => runScan());
    $("scanCsvBtn").addEventListener("click", () => exportCsv());
    // 参数列表依赖当前配置：每次渲染左栏时刷新（保留已选项）。
    fillParams();
    fillMetrics();
  }

  // fillParams 刷新参数下拉（保留已选项）。配置可能还没载入（app.js 异步读
  // catalog），所以这里对空配置必须静默返回，不能抛异常把面板留半成品。
  function fillParams() {
    const sel = $("#scanParam");
    if (!sel || !S.config) return;
    let paths;
    try {
      paths = scanPaths();
    } catch (e) {
      return; // 配置形状异常：等下一次刷新
    }
    if (!paths.length) return;
    const prev = sel.value;
    sel.innerHTML = "";
    paths.forEach(([p, label]) => {
      const o = document.createElement("option");
      o.value = p;
      o.textContent = label;
      sel.appendChild(o);
    });
    if (prev && paths.some(([p]) => p === prev)) sel.value = prev;
  }

  function fillMetrics() {
    const sel = $("#scanMetric");
    if (!sel) return;
    sel.innerHTML = "";
    METRICS.forEach(([k, label]) => {
      const o = document.createElement("option");
      o.value = k;
      o.textContent = label;
      if (k === "power") o.selected = true;
      sel.appendChild(o);
    });
  }

  function fmt(v) {
    if (v === null || v === undefined || Number.isNaN(v)) return "–";
    const a = Math.abs(v);
    if (a !== 0 && (a >= 1e4 || a < 1e-3)) return v.toExponential(3);
    return v.toPrecision(6);
  }

  // parseValues 与服务端 CLI 同一套语法：逗号列表 + a-b / a-b-step 范围。
  function parseValues(s) {
    const out = [];
    for (const part0 of String(s).split(",")) {
      const part = part0.trim();
      if (!part) continue;
      const bits = [];
      let start = 0;
      for (let i = 1; i < part.length; i++) {
        if (part[i] !== "-") continue;
        const prev = part[i - 1];
        if ((prev >= "0" && prev <= "9") || prev === ".") { bits.push(part.slice(start, i)); start = i + 1; }
      }
      bits.push(part.slice(start));
      if (bits.length === 1) {
        const v = Number(bits[0]);
        if (!Number.isFinite(v)) throw new Error("无法解析取值 " + part);
        out.push(v);
        continue;
      }
      if (bits.length > 3) throw new Error("范围 " + part + " 最多是 a-b-step 三段");
      const lo = Number(bits[0]), hi = Number(bits[1]);
      const step = bits.length === 3 ? Number(bits[2]) : 1;
      if (![lo, hi, step].every(Number.isFinite)) throw new Error("无法解析范围 " + part);
      if (step === 0) throw new Error("范围 " + part + " 的步长不能为 0");
      if ((hi > lo) !== (step > 0)) throw new Error("范围 " + part + " 的步长方向与起止相反");
      const n = Math.floor((hi - lo) / step + 1e-9);
      if (n + 1 > 512) throw new Error("范围 " + part + " 超过 512 个点");
      for (let k = 0; k <= n; k++) out.push(lo + step * k);
    }
    if (!out.length) throw new Error("没有给出任何取值");
    return out;
  }

  let lastReport = null;

  async function runScan() {
    const param = $("#scanParam").value;
    const metric = $("#scanMetric").value;
    let values;
    try {
      values = parseValues($("scanValues").value);
    } catch (e) {
      $("scanOut").textContent = "取值列表有问题：" + e.message;
      return;
    }
    $("scanRunBtn").disabled = true;
    $("scanOut").textContent = "计算中…（" + values.length + " 个点）";
    try {
      const r = await fetch("/api/scan", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ config: S.config, scan: { param: param, values: values, metric: metric } }),
      });
      const j = await r.json().catch(() => ({}));
      if (!r.ok) { $("scanOut").textContent = "请求被拒绝：" + (j.error || r.status); return; }
      const rid = j.run_id;
      for (;;) {
        const m = await fetch("/api/runs/" + rid).then((x) => x.json());
        if (m.status === "done") { lastReport = m.scan; showReport(m.scan); return; }
        if (m.status === "error") { $("scanOut").textContent = "计算错误：" + (m.error || "?"); return; }
        await new Promise((res) => setTimeout(res, 250));
      }
    } catch (e) {
      $("scanOut").textContent = "请求失败：" + e.message;
    } finally {
      $("scanRunBtn").disabled = false;
      $("scanCsvBtn").disabled = !lastReport;
    }
  }

  function showReport(rep) {
    if (!rep) { $("scanOut").textContent = "空报告"; return; }
    const lines = ["参数 " + rep.param + "（指标 " + rep.metric + "）"];
    (rep.points || []).forEach((p, i) => {
      if (p.error) { lines.push("  [" + i + "] " + fmt(p.value) + "  错误: " + String(p.error).split("\n")[0]); return; }
      lines.push("  [" + i + "] " + fmt(p.value) + "  " + fmt(p.metric));
    });
    if (rep.warnings && rep.warnings.length) {
      lines.push("警告（去重）:");
      rep.warnings.forEach((w) => lines.push("  " + w));
    }
    $("scanOut").textContent = lines.join("\n");
    drawCurve(rep);
  }

  // drawCurve 把指标随参数的变化画成折线（线性坐标，自动范围）。
  function drawCurve(rep) {
    const cv = $("#scanPlot");
    if (!cv) return;
    const W = cv.clientWidth || 320, H = 150;
    cv.width = W; cv.height = H;
    const ctx = cv.getContext("2d");
    ctx.clearRect(0, 0, W, H);
    const pts = (rep.points || []).filter((p) => !p.error && Number.isFinite(p.metric));
    if (pts.length < 2) {
      ctx.fillStyle = "#33475c";
      ctx.font = "12px system-ui, sans-serif";
      ctx.fillText("点数不足（<2），无法画曲线", 12, 20);
      return;
    }
    let vmin = Infinity, vmax = -Infinity;
    for (const p of pts) { if (p.metric < vmin) vmin = p.metric; if (p.metric > vmax) vmax = p.metric; }
    if (vmax - vmin < 1e-30) { vmax = vmin + 1; }
    const pad = 34;
    const px = (i) => pad + (i / (pts.length - 1)) * (W - pad * 2);
    const py = (v) => H - pad - ((v - vmin) / (vmax - vmin)) * (H - pad * 2);
    ctx.strokeStyle = "rgba(20,40,60,.25)";
    ctx.beginPath(); ctx.moveTo(pad, H - pad); ctx.lineTo(W - pad, H - pad); ctx.stroke();
    ctx.strokeStyle = "#1460c8";
    ctx.lineWidth = 1.6;
    ctx.beginPath();
    pts.forEach((p, i) => { const x = px(i), y = py(p.metric); if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y); });
    ctx.stroke();
    ctx.fillStyle = "#33475c";
    ctx.font = "11px system-ui, sans-serif";
    ctx.fillText(fmt(vmin), 2, H - pad + 4);
    ctx.fillText(fmt(vmax), 2, pad - 2);
    ctx.fillText(fmt(pts[0].value), pad - 6, H - pad + 14);
    ctx.fillText(fmt(pts[pts.length - 1].value), W - pad - 24, H - pad + 14);
  }

  function exportCsv() {
    const rep = lastReport;
    if (!rep) return;
    const rows = [["value", rep.metric || "metric", "error"].join(",")];
    (rep.points || []).forEach((p) => rows.push([p.value, p.error ? "" : p.metric, JSON.stringify(p.error || "")].join(",")));
    const blob = new Blob([rows.join("\n") + "\n"], { type: "text/csv" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = "scan.csv";
    document.body.appendChild(a);
    a.click();
    a.remove();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
    if (typeof setStatus === "function") setStatus("已导出 scan.csv");
  }

  // init：面板构建 + 定期刷新参数列表（app.js 异步载入配置与目录，
  // 之后用户还会增删元件/光源，所以用一个低频轮询把下拉列表跟上）。
  function init() {
    buildPanel();
    fillMetrics();
    const tick = () => { if (S.config) fillParams(); };
    tick();
    setInterval(tick, 1000);
    // 切回波动模式时立即刷新（oz/quantum 模式下左栏是隐藏的）。
    document.addEventListener("oz:mode", (ev) => { if (ev.detail === "wave") tick(); });
  }

  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", init);
  else init();
})();
