/* Tetsuhiro-OZ 几何 / 膜系标签。
 *
 * 与波动、量子标签共用同一个页面（index.html），由 app.js 的 setMode("oz") 切换可见性。
 * 本文件是独立 IIFE：不向全局作用域写入任何名字（app.js 的 $、clone 等互不干扰），
 * 只读写本标签自己的 DOM（id 以 oz 开头）并调用 /api/version、/api/lens/*、/api/coating/*。
 * 计算全部在服务端完成；这里只做表单读取、结果展示与光谱绘制。 */
(function () {
  "use strict";

  const $ = (id) => document.getElementById(id);
  const LAMBDA_FROM = 400, LAMBDA_TO = 750;

  let prescription = null; // 最近一次设计得到的处方（评估与公差分析使用）
  let designEFL = 0;       // 该处方的目标焦距（mm），公差分析的焦距判据
  let stack = null;        // 最近一次膜系设计得到的膜层结构
  let rows = null;         // 最近一次光谱 [{lambda_nm, R, T, A}]

  // api POSTs JSON and returns the decoded answer; a non-2xx answer becomes an
  // Error carrying the server's {"error": ...} text.
  async function api(path, body) {
    const r = await fetch(path, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(body),
    });
    const text = await r.text();
    let d;
    try { d = JSON.parse(text); } catch (e) { d = { raw: text }; }
    if (!r.ok) throw new Error(d && d.error ? d.error : "HTTP " + r.status);
    return d;
  }

  function show(id, value, isErr) {
    const el = $(id);
    el.classList.toggle("err", !!isErr);
    el.textContent = typeof value === "string" ? value : JSON.stringify(value, null, 2);
  }

  function num(id) {
    const v = parseFloat($(id).value);
    return isNaN(v) ? 0 : v;
  }

  function txt(id) {
    return $(id).value.trim();
  }

  // busy disables a button while its request runs and always restores it.
  async function busy(btn, fn) {
    btn.disabled = true;
    try { return await fn(); } finally { btn.disabled = false; }
  }

  function setStatus(msg, isErr) {
    const el = $("ozLensStatus");
    el.textContent = msg;
    el.classList.toggle("err", !!isErr);
  }

  async function design() {
    const d = await api("/api/lens/design", {
      efl_mm: num("ozEfl"),
      semi_aperture_mm: num("ozSemi"),
      field_deg: num("ozField"),
      crown: txt("ozCrown") || "N-BK7",
      flint: txt("ozFlint") || "F2",
    });
    prescription = d.prescription || null;
    designEFL = num("ozEfl");
    $("ozEvalBtn").disabled = !prescription;
    $("ozTolBtn").disabled = !prescription;
    const opt = d.optimizer || {};
    setStatus("设计完成 · 迭代 " + (opt.iterations ?? "—") + " 次 · merit " +
      (typeof opt.merit === "number" ? opt.merit.toPrecision(4) : "—"), false);
    show("ozLensOut", d.report || d, false);
    show("ozPrescOut", prescription || {}, false);
  }

  async function evaluate() {
    if (!prescription) { show("ozLensOut", "请先完成一次设计。", true); return; }
    const d = await api("/api/lens/eval", prescription);
    const m = d.metrics || {};
    const lines = [d.report || ""];
    if (typeof m.efl_mm === "number") lines.push("EFL = " + m.efl_mm.toFixed(4) + " mm");
    if (typeof m.bfl_mm === "number") lines.push("BFL = " + m.bfl_mm.toFixed(4) + " mm");
    if (typeof m.f_number === "number") lines.push("F/# = " + m.f_number.toFixed(3));
    if (typeof m.color_shift_mm === "number") lines.push("F–C 色差 = " + m.color_shift_mm.toFixed(5) + " mm");
    show("ozLensOut", lines.join("\n"), false);
  }

  async function tolerance() {
    if (!prescription) { show("ozTolOut", "请先完成一次设计。", true); return; }
    const d = await api("/api/lens/tol", {
      prescription: prescription,
      spec: { radius_mm: num("ozTolR"), thickness_mm: num("ozTolT"), index_delta: num("ozTolN") },
      criteria: { target_efl_mm: designEFL, max_focal_err_mm: 0.1 },
      runs: Math.max(1, Math.round(num("ozTolRuns"))),
      seed: 1,
    });
    show("ozTolOut", {
      runs: d.runs,
      yield: d.yield,
      nominal: d.nominal,
      on_axis_rms: d.on_axis_rms,
      off_axis_rms: d.off_axis_rms,
      color_shift: d.color_shift,
      focal_err: d.focal_err,
      fail_reasons: d.fail_reasons,
    }, false);
  }

  async function designCoating() {
    const d = await api("/api/coating/design", {
      template: $("ozTemplate").value,
      lambda_nm: num("ozLambda"),
      from_nm: num("ozFrom"),
      to_nm: num("ozTo"),
      pairs: Math.round(num("ozPairs")),
      cutoff_nm: num("ozCutoff"),
      incident: txt("ozInc") || "Air",
      substrate: txt("ozSub") || "BK7",
    });
    stack = d.stack || null;
    show("ozCoatOut", d.summary || d, false);
    return d;
  }

  async function spectrum() {
    const d = await designCoating();
    if (!stack) return;
    const s = await api("/api/coating/spectrum", {
      stack: stack, from_nm: LAMBDA_FROM, to_nm: LAMBDA_TO, step_nm: 2, angle_deg: 0, pol: "u",
    });
    rows = s.rows || [];
    $("ozCoatOut").textContent = (d.summary || "") + "\n\n光谱：" + rows.length + " 点（" +
      LAMBDA_FROM + "–" + LAMBDA_TO + " nm，0°，非偏振）";
    drawPlot();
  }

  /* --- 光谱图：宽度跟随容器，按 devicePixelRatio 设置后备像素 --- */
  function fitCanvas() {
    const cv = $("ozPlot");
    const dpr = window.devicePixelRatio || 1;
    const w = Math.max(1, cv.clientWidth), h = Math.max(1, cv.clientHeight);
    cv.width = Math.round(w * dpr);
    cv.height = Math.round(h * dpr);
    const c = cv.getContext("2d");
    c.setTransform(dpr, 0, 0, dpr, 0, 0);
    return { c: c, w: w, h: h };
  }

  function drawPlot() {
    if ($("ozOutput").hidden) return;
    const f = fitCanvas(), c = f.c, w = f.w, h = f.h;
    c.clearRect(0, 0, w, h);
    const x0 = Math.max(34, Math.round(w * 0.085)), x1 = w - 10, y0 = 22, y1 = h - 22;
    c.font = "11px system-ui,sans-serif";
    c.strokeStyle = "#7d97ad"; c.lineWidth = 1;
    c.beginPath(); c.moveTo(x0, y0); c.lineTo(x0, y1); c.lineTo(x1, y1); c.stroke();
    c.textAlign = "right"; c.fillStyle = "#5a7d9e";
    for (const v of [0, 0.5, 1]) {
      const y = y1 - v * (y1 - y0);
      c.strokeStyle = v === 0 ? "#7d97ad" : "rgba(125,151,173,0.35)";
      c.beginPath(); c.moveTo(x0, y); c.lineTo(x1, y); c.stroke();
      c.fillText(Math.round(v * 100) + "%", x0 - 5, y + 4);
    }
    c.textAlign = "center";
    for (let lam = LAMBDA_FROM; lam <= LAMBDA_TO; lam += 50) {
      const x = x0 + (lam - LAMBDA_FROM) / (LAMBDA_TO - LAMBDA_FROM) * (x1 - x0);
      c.strokeStyle = "#7d97ad";
      c.beginPath(); c.moveTo(x, y1); c.lineTo(x, y1 + 4); c.stroke();
      c.fillStyle = "#5a7d9e"; c.fillText(String(lam), x, y1 + 15);
    }
    c.textAlign = "left"; c.fillStyle = "#5a7d9e";
    c.fillText("R(λ) %", x0, y0 - 9);
    c.textAlign = "right";
    c.fillText("λ (nm) · 400–750 nm · 0°", x1, y0 - 9);
    if (!rows || !rows.length) {
      c.textAlign = "center"; c.fillStyle = "#5a7d9e";
      c.fillText("点「设计并画光谱」绘制反射率曲线", (x0 + x1) / 2, (y0 + y1) / 2);
      return;
    }
    c.strokeStyle = "#1976d2"; c.lineWidth = 2; c.beginPath();
    let first = true;
    for (const p of rows) {
      const x = x0 + (p.lambda_nm - LAMBDA_FROM) / (LAMBDA_TO - LAMBDA_FROM) * (x1 - x0);
      const y = y1 - p.R * (y1 - y0);
      if (first) { c.moveTo(x, y); first = false; } else { c.lineTo(x, y); }
    }
    c.stroke();
  }

  // guard wraps a click handler so a failure is shown in the result pane
  // instead of an unhandled promise rejection.
  function guard(btnId, outId, fn) {
    $(btnId).addEventListener("click", async function () {
      const btn = $(btnId);
      try {
        await busy(btn, fn);
      } catch (e) {
        show(outId, "错误：" + e.message, true);
      }
    });
  }

  async function loadVersion() {
    try {
      const r = await fetch("/api/version");
      const d = await r.json();
      $("ozVersion").textContent = "v" + d.version;
    } catch (e) {
      $("ozVersion").textContent = "";
    }
  }

  function init() {
    guard("ozDesignBtn", "ozLensOut", design);
    guard("ozEvalBtn", "ozLensOut", evaluate);
    guard("ozTolBtn", "ozTolOut", tolerance);
    guard("ozCoatBtn", "ozCoatOut", designCoating);
    guard("ozSpecBtn", "ozCoatOut", spectrum);
    window.addEventListener("resize", drawPlot);
    window.addEventListener("orientationchange", drawPlot);
    if (window.ResizeObserver) new ResizeObserver(drawPlot).observe($("ozPlot"));
    // app.js announces mode changes; the canvas needs a layout pass to be sized.
    document.addEventListener("oz:mode", function (ev) {
      if (ev.detail === "oz") drawPlot();
    });
    loadVersion();
    drawPlot();
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init);
  } else {
    init();
  }
})();
