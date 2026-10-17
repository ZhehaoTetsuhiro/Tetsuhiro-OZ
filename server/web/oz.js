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

  // setGroupStatus 写透镜组卡片标题旁的一行摘要（与立体视图的 #stats 同位置感）。
  function setGroupStatus(msg, isErr) {
    const el = $("ozGroupStatus");
    if (!el) return;
    el.textContent = msg || "";
    el.classList.toggle("err", !!isErr);
  }

  /* --- 高级分析（热漂移 / 灵敏度）：折叠区由本文件构建，index.html 不感知 ---
   * 面板复用页面已有的 ozGrid/ozField/ozActions/ozDetails/ozOut 类，所以手机断点与
   * 触控目标规则（style.css 的 @media 规则）对它同样生效。 */

  function buildAdvanced() {
    // 高级分析是消色差设计的延伸（温度与灵敏度都针对那条处方），所以挂在
    // #ozPanelAchromat 里：跟随该子模式一起显隐，而不是常驻在左栏底部。
    const host = $("ozPanelAchromat");
    if (!host || $("ozAdv")) return;
    const det = document.createElement("details");
    det.className = "ozDetails";
    det.id = "ozAdv";
    det.innerHTML =
      '<summary>高级分析（热漂移 · 灵敏度 · 逆灵敏度）</summary>' +
      '<div class="ozGrid">' +
      '<label class="ozField"><span>温度变化 ΔT（K）</span>' +
      '<input type="number" id="ozAdvDT" value="10" step="1" inputmode="decimal"></label>' +
      '<label class="ozField"><span>镜筒线膨胀 α（1/K）</span>' +
      '<input type="number" id="ozAdvCTE" value="0" step="0.000001" min="0" inputmode="decimal"></label>' +
      '<label class="ozField"><span>merit 预算 ΔM（mm）</span>' +
      '<input type="number" id="ozAdvBudget" value="0.005" step="0.001" min="0" inputmode="decimal"></label>' +
      '</div>' +
      '<div class="ozActions">' +
      '<button id="ozAdvThermalBtn" type="button" class="secondary" disabled>热漂移</button>' +
      '<button id="ozAdvSensBtn" type="button" class="secondary" disabled>灵敏度与公差分配</button>' +
      '</div>' +
      '<pre class="ozOut" id="ozAdvOut" aria-live="polite">完成一次设计后可用。</pre>';
    host.appendChild(det);
  }

  function setAdvanced(enabled) {
    for (const id of ["ozAdvThermalBtn", "ozAdvSensBtn"]) {
      if ($(id)) $(id).disabled = !enabled;
    }
  }

  async function thermal() {
    if (!prescription) { show("ozAdvOut", "请先完成一次设计。", true); return; }
    const d = await api("/api/lens/thermal", {
      prescription: prescription,
      delta_t_k: num("ozAdvDT"),
      mount_cte_per_k: num("ozAdvCTE"),
    });
    show("ozAdvOut", [
      "λ = " + d.lambda_nm.toFixed(4) + " nm · ΔT = " + d.delta_t_k + " K · α = " + d.mount_cte_per_k + " 1/K",
      "EFL  " + d.efl0_mm.toFixed(6) + " → " + d.efl1_mm.toFixed(6) + " mm  (" + d.efl_per_k_mm.toExponential(3) + " mm/K)",
      "BFL  " + d.bfl0_mm.toFixed(6) + " → " + d.bfl1_mm.toFixed(6) + " mm  (" + d.bfl_per_k_mm.toExponential(3) + " mm/K)",
    ].join("\n"), false);
  }

  async function sensitivity() {
    if (!prescription) { show("ozAdvOut", "请先完成一次设计。", true); return; }
    const d = await api("/api/lens/sens", {
      prescription: prescription,
      metric: "merit",
      inverse: true,
      budget_mm: num("ozAdvBudget"),
    });
    const lines = ["标量 = " + d.metric + " · 名义 " + d.nominal.toPrecision(6) + " mm · " + d.evals + " 次求值"];
    for (const r of d.rows) {
      lines.push(r.param + "  ∂ = " + r.deriv.toPrecision(6) + (r.converged ? "" : "（未收敛）"));
    }
    if (d.allocation) {
      const a = d.allocation;
      lines.push("", "逆灵敏度分配（" + a.mode + "，预算 " + a.budget_mm + " mm，预测 " +
        a.predicted_mm.toPrecision(6) + " mm，相对误差 " + a.rel_err.toExponential(2) + "）");
      for (const r of a.rows) {
        lines.push(r.param + "  公差 ±" + r.tol.toPrecision(4) + (r.capped ? "（封顶）" : ""));
      }
    }
    show("ozAdvOut", lines.join("\n"), false);
  }

  /* --- 透镜组（3D，含镜）--------------------------------------------------
   * 「透镜组」是几何光学的默认子模式：面板直接挂在 index.html 的
   * #ozPanelGroup 里（没有折叠区，也不需要点「载入示例」之外的手动展开），
   * 按钮全部走 /api/lens/group/*（docs/API.md §4）。表格里的 kind / 玻璃 /
   * 间隔介质 / fold 由前端写进 group JSON，数值参数走 /set 的 sets —— 服务端
   * 校验并把错误文本原样回显。示例组的 JSON 由 cmd/oz 的测试用「oz:example」
   * 标记提取并逐个解码校验（见 TestGroupExamplesInGUI）。 */

  const GROUP_EXAMPLES = [
    {
      label: "卡塞格林（两镜 · 库例 F1000 d300 b300）",
      group: /*oz:example*/{
        "name": "卡塞格林（两镜，球面，库例 F1000 d300 b300）",
        "object_z_mm": "inf",
        "elements": [
          {
            "kind": "stop",
            "front": { "radius_mm": 0, "semi_aperture_mm": 150 },
            "thickness_mm": 0,
            "gap_mm": 60
          },
          {
            "kind": "mirror",
            "front": { "radius_mm": -2000, "semi_aperture_mm": 153 },
            "thickness_mm": 0,
            "gap_mm": -300
          },
          {
            "kind": "mirror",
            "front": { "radius_mm": 8400, "semi_aperture_mm": 117.5 },
            "thickness_mm": 0,
            "gap_mm": 600
          }
        ],
        "stop": 0
      }/*:oz:example*/
    },
    {
      label: "BK7/F2 双胶合 f=100 mm（设计值）",
      group: /*oz:example*/{
        "name": "BK7/F2 双胶合 f=100 mm（设计值）",
        "object_z_mm": "inf",
        "elements": [
          {
            "name": "crown front",
            "kind": "lens",
            "front": { "radius_mm": 44.1685612361813, "semi_aperture_mm": 10 },
            "rear": { "radius_mm": -44.13070097892217, "semi_aperture_mm": 10 },
            "glass": "N-BK7",
            "gap_glass": "F2",
            "thickness_mm": 2.294765132716998,
            "gap_mm": 0.6284096143319601
          },
          {
            "name": "flint rear",
            "kind": "lens",
            "front": { "radius_mm": -1016.2003173810598, "semi_aperture_mm": 10 },
            "thickness_mm": 0,
            "gap_mm": 95
          }
        ],
        "stop": 0
      }/*:oz:example*/
    },
    {
      label: "折叠式开普勒（45° 折转 · Γ=−1）",
      group: /*oz:example*/{
        "name": "折叠式开普勒（45° 折转 + 两片 f=100）",
        "object_z_mm": "inf",
        "elements": [
          {
            "kind": "stop",
            "front": { "radius_mm": 0, "semi_aperture_mm": 6 },
            "thickness_mm": 0,
            "gap_mm": 15
          },
          {
            "kind": "mirror",
            "front": { "radius_mm": 0, "semi_aperture_mm": 20 },
            "thickness_mm": 0,
            "gap_mm": 20,
            "tilt_x_deg": 45
          },
          {
            "kind": "lens",
            "front": { "radius_mm": 51.680003450058855, "semi_aperture_mm": 20 },
            "rear": { "radius_mm": 0, "semi_aperture_mm": 20 },
            "glass": "N-BK7",
            "thickness_mm": 2,
            "gap_mm": 198.68143462914773
          },
          {
            "kind": "lens",
            "front": { "radius_mm": 51.680003450058855, "semi_aperture_mm": 20 },
            "rear": { "radius_mm": 0, "semi_aperture_mm": 20 },
            "glass": "N-BK7",
            "thickness_mm": 2,
            "gap_mm": 0
          }
        ],
        "stop": 0
      }/*:oz:example*/
    },
  ];

  let groupObj = null;                         // 当前透镜组（wire form）
  let groupLayout = null;                      // 最近一次 layout 响应（画布数据）
  let groupView = { yaw: -0.62, pitch: 0.45 }; // 画布视角（弧度）
  let groupSel = 0;                            // 当前选中的元件下标

  const GROUP_KINDS = [["lens", "lens 透镜"], ["mirror", "mirror 反射镜"], ["stop", "stop 光阑"]];

  function buildGroup() {
    const host = $("ozPanelGroup");
    if (!host || $("ozGroup")) return;
    const opts = GROUP_EXAMPLES.map(function (ex, i) {
      return '<option value="' + i + '">' + ex.label + "</option>";
    }).join("");
    const kindOpts = GROUP_KINDS.map(function (k) {
      return '<option value="' + k[0] + '">' + k[1] + "</option>";
    }).join("");
    const det = document.createElement("details");
    det.className = "ozDetails";
    det.id = "ozGroup";
    det.open = true; // 默认子模式：展开着，没有多余的点击
    // 左栏是「波动光学」那套编辑方式：物（位置/大小）+ 元件列表 + 选中元件的
    // 参数（位置/大小/形状），外加 / 复制 / 删除。画布与结果文本在右侧主区
    // （#ozGroupView），一次只显示一个元件的字段，而不是把所有元件摊在左栏。
    det.innerHTML =
      '<summary>透镜组（3D，含镜）</summary>' +
      '<div class="subhead">物（位置 · 大小）</div>' +
      '<div class="prow"><label for="ozObjMode" title="无穷远＝平行光入射（用视场角）；有限距离＝物在 z 处（用物高）">物在</label>' +
      '<select id="ozObjMode"><option value="inf">无穷远</option><option value="near">有限距离</option></select></div>' +
      '<div class="prow" id="ozObjDistRow" hidden><label for="ozObjDist">物距（mm）</label>' +
      '<input type="number" id="ozObjDist" value="500" min="0.001" step="10" inputmode="decimal"></div>' +
      '<div class="prow" id="ozObjDegRow"><label for="ozObjFieldDeg">视场角（°）</label>' +
      '<input type="number" id="ozObjFieldDeg" value="0" step="0.1" inputmode="decimal"></div>' +
      '<div class="prow" id="ozObjHRow" hidden><label for="ozObjH" title="有限距离物的物高（mm，子午 y 方向）">物高（mm）</label>' +
      '<input type="number" id="ozObjH" value="0" step="0.1" inputmode="decimal"></div>' +
      '<div class="prow"><label for="ozGroupLambda">波长（nm）</label>' +
      '<input type="number" id="ozGroupLambda" value="587.5618" step="1" inputmode="decimal"></div>' +
      '<div class="prow"><label for="ozGroupPupil" title="采样光瞳半径；0 = 用光阑半口径">光瞳半径（mm）</label>' +
      '<input type="number" id="ozGroupPupil" value="0" min="0" step="1" inputmode="decimal"></div>' +
      '<div class="subhead">元件（<span id="ozGroupCount">0</span>）</div>' +
      '<div id="ozGroupList" role="listbox" aria-label="元件列表" tabindex="0"></div>' +
      '<div class="prow"><label for="ozGroupAddKind">新增</label><select id="ozGroupAddKind">' + kindOpts + "</select></div>" +
      '<div class="row">' +
      '<button id="ozGroupAddBtn" type="button">+ 元件</button>' +
      '<button id="ozGroupDupBtn" type="button">复制</button>' +
      '<button id="ozGroupDelBtn" type="button">删除</button>' +
      "</div>" +
      '<div id="ozGroupEditor" aria-label="选中元件参数"></div>' +
      '<div class="ozActions"><button id="ozGroupApplyBtn" type="button" disabled>应用编辑</button></div>' +
      '<div class="prow"><label for="ozGroupExample">示例组</label><select id="ozGroupExample">' + opts + "</select></div>" +
      '<div class="ozActions">' +
      '<button id="ozGroupLoadBtn" type="button" class="secondary">载入示例</button>' +
      '<button id="ozGroupLayoutBtn" type="button" class="secondary" disabled>布局</button>' +
      '<button id="ozGroupImageBtn" type="button" class="secondary" disabled>成像位置</button>' +
      '<button id="ozGroupAfocalBtn" type="button" class="secondary" disabled>角放大率</button>' +
      '<button id="ozGroupDispBtn" type="button" class="secondary" disabled>色散</button>' +
      '<button id="ozGroupPrescBtn" type="button" class="secondary" disabled>转处方</button>' +
      "</div>" +
      '<p class="hint">选元件（↑/↓ 也行）→ 改下面的位置/大小/形状 → 「应用编辑」，服务端逐项校验并原样回显错误；+ 元件 / 复制 / 删除 改的是光路本身，校验失败会退回。</p>';
    host.appendChild(det);
  }

  // buildGroupView 在右侧主区搭出画布与结果文本（index.html 只留 #ozGroupView
  // 这个空容器）。元件编辑在左栏，这里不再有表格。
  function buildGroupView() {
    const host = $("ozGroupView");
    if (!host || $("ozGroupPlot")) return;
    host.innerHTML =
      '<div class="ozGroupBar">' +
      '<button id="ozGroupIsoBtn" type="button" class="secondary miniBtn">等轴测</button>' +
      '<button id="ozGroupXZBtn" type="button" class="secondary miniBtn">XZ 面</button>' +
      '<button id="ozGroupXYBtn" type="button" class="secondary miniBtn">XY 面</button>' +
      '<span class="ozHint" id="ozGroupViewHint">拖拽画布旋转视角</span>' +
      "</div>" +
      '<div id="ozGroupCanvasWrap">' +
      '<canvas id="ozGroupPlot" aria-label="透镜组 3D 布局（拖拽旋转）"></canvas>' +
      "</div>" +
      '<pre class="ozOut" id="ozGroupOut" aria-live="polite">点「载入示例」开始。</pre>';
  }

  function setGroupButtons(enabled) {
    const ids = ["ozGroupApplyBtn", "ozGroupLayoutBtn", "ozGroupImageBtn",
      "ozGroupAfocalBtn", "ozGroupDispBtn", "ozGroupPrescBtn"];
    for (const id of ids) {
      if ($(id)) $(id).disabled = !enabled;
    }
  }

  /* --- 选中元件的参数面板：一行一个 .prow（与波动光学的参数面板同款） --- */
  function pRow(parent, label, title) {
    const row = document.createElement("div");
    row.className = "prow";
    const l = document.createElement("label");
    l.textContent = label;
    if (title) l.title = title;
    row.appendChild(l);
    parent.appendChild(row);
    return row;
  }

  function pNum(parent, label, id, value, step, title) {
    const row = pRow(parent, label, title);
    const i = document.createElement("input");
    i.type = "number";
    i.id = id;
    i.step = String(step);
    i.inputMode = "decimal";
    i.value = (value === undefined || value === null) ? "" : String(value);
    row.appendChild(i);
    return i;
  }

  function pText(parent, label, id, value, title) {
    const row = pRow(parent, label, title);
    const i = document.createElement("input");
    i.type = "text";
    i.id = id;
    i.spellcheck = false;
    i.autocapitalize = "off";
    i.value = value === undefined || value === null ? "" : String(value);
    row.appendChild(i);
    return i;
  }

  function pSel(parent, label, id, options, value, title) {
    const row = pRow(parent, label, title);
    const sel = document.createElement("select");
    sel.id = id;
    for (const o of options) {
      const op = document.createElement("option");
      op.value = o[0];
      op.textContent = o[1];
      sel.appendChild(op);
    }
    sel.value = value;
    row.appendChild(sel);
    return sel;
  }

  function pSubhead(parent, text) {
    const h = document.createElement("div");
    h.className = "subhead";
    h.textContent = text;
    parent.appendChild(h);
  }

  // kindDot 给元件列表的小圆点上色：透镜蓝、反射镜橙、光阑灰（与波动的 .sc 一致）。
  function kindDot(kind) {
    if (kind === "mirror") return "#e8590c";
    if (kind === "stop") return "#868e96";
    return "#1976d2";
  }

  // elementTitle 是列表行上的一句话摘要：类型 + 关键尺寸 + 口径。
  function elementTitle(e) {
    const bits = [e.kind];
    if (e.front && e.front.radius_mm) bits.push("R=" + fx(e.front.radius_mm));
    if (e.rear) bits.push("R2=" + fx(e.rear.radius_mm));
    if (e.glass) bits.push(e.glass);
    return bits.join(" · ");
  }

  function renderGroupList() {
    const list = $("ozGroupList");
    if (!list) return;
    list.textContent = "";
    const n = groupObj && groupObj.elements ? groupObj.elements.length : 0;
    const count = $("ozGroupCount");
    if (count) count.textContent = String(n);
    if (!groupObj) return;
    groupObj.elements.forEach(function (e, i) {
      const b = document.createElement("button");
      b.type = "button";
      b.className = "sc" + (i === groupSel ? " sel" : "");
      b.dataset.i = String(i);
      b.setAttribute("role", "option");
      b.setAttribute("aria-selected", i === groupSel ? "true" : "false");
      const dot = document.createElement("span");
      dot.className = "dot";
      dot.style.background = kindDot(e.kind);
      b.appendChild(dot);
      const name = document.createElement("span");
      name.textContent = "[" + i + "] " + elementTitle(e);
      b.appendChild(name);
      const tag = document.createElement("span");
      tag.className = "tag";
      tag.textContent = e.front && e.front.semi_aperture_mm ? "口径 " + e.front.semi_aperture_mm : "";
      b.appendChild(tag);
      b.addEventListener("click", function () { selectGroupEl(i); });
      list.appendChild(b);
    });
  }

  // renderGroupEditor rebuilds the parameter panel of the selected element:
  // 位置（gap、偏心、倾斜）→ 大小（半口径）→ 形状（R、K、后表面、玻璃、间隔介质、折转）。
  function renderGroupEditor() {
    const host = $("ozGroupEditor");
    if (!host) return;
    host.textContent = "";
    if (!groupObj || !groupObj.elements || !groupObj.elements.length) {
      const p = document.createElement("p");
      p.className = "hint";
      p.textContent = "没有元件：用「+ 元件」加一个。";
      host.appendChild(p);
      return;
    }
    if (groupSel >= groupObj.elements.length) groupSel = groupObj.elements.length - 1;
    const e = groupObj.elements[groupSel];
    const box = document.createElement("div");
    box.className = "ozGroupEdit";
    host.appendChild(box);
    pStyleSubhead(box, "选中元件 [" + groupSel + "]");
    pSel(box, "类型", "ozGEKind", GROUP_KINDS, e.kind, "lens 有前后两面与玻璃；mirror/stop 只有前表面");

    pStyleSubhead(box, "位置");
    pNum(box, "gap（mm）", "ozGEGap", e.gap_mm, 1, "到下一元件的轴向距离，可负（轴向折返）");
    pNum(box, "偏心 dx（mm）", "ozGEDx", e.decenter_x_mm, 0.1);
    pNum(box, "偏心 dy（mm）", "ozGEDy", e.decenter_y_mm, 0.1);
    pNum(box, "倾斜 tx（°）", "ozGETx", e.tilt_x_deg, 0.5);
    pNum(box, "倾斜 ty（°）", "ozGETy", e.tilt_y_deg, 0.5);
    pNum(box, "倾斜 tz（°）", "ozGETz", e.tilt_z_deg, 0.5);

    pStyleSubhead(box, "大小");
    pNum(box, "半口径（mm）", "ozGESemi",
      e.front ? e.front.semi_aperture_mm : null, 0.5, "前表面通光半口径");
    if (e.rear) {
      pNum(box, "后半口径（mm）", "ozGESemi2", e.rear.semi_aperture_mm, 0.5);
    }

    pStyleSubhead(box, "形状");
    pNum(box, "前表面 R（mm）", "ozGER1", e.front ? e.front.radius_mm : null, 1, "0 = 平面");
    pNum(box, "前表面 K", "ozGEK1", e.front ? e.front.conic_k : null, 0.1, "圆锥常数，0 = 球面");
    if (e.kind === "lens") {
      pNum(box, "后表面 R（mm）", "ozGER2", e.rear ? e.rear.radius_mm : null, 1);
      pNum(box, "后表面 K", "ozGEK2", e.rear ? e.rear.conic_k : null, 0.1);
      pNum(box, "中心厚度（mm）", "ozGEThk", e.thickness_mm, 0.1);
      pText(box, "玻璃", "ozGEGlass", e.glass, "如 N-BK7 / F2 / Air");
    }
    if (e.kind !== "stop") {
      pText(box, "间隔介质", "ozGEGapGlass", e.gap_glass, "到下一元件的介质，空 = 空气；胶合面写后一片的玻璃");
    }
    pSel(box, "折转 fold", "ozGEFold",
      [["auto", "auto（倾斜镜折转链）"], ["off", "off（不折转）"]],
      (e.fold === 1 || e.fold === "off") ? "off" : "auto",
      "auto：倾斜的反射镜把链轴转到反射光束");
  }

  // pStyleSubhead 是编辑器内部的小标题（比 .subhead 轻一点）。
  function pStyleSubhead(parent, text) {
    const h = document.createElement("h4");
    h.textContent = text;
    parent.appendChild(h);
  }

  function selectGroupEl(i) {
    if (!groupObj || !groupObj.elements) return;
    groupSel = Math.max(0, Math.min(i, groupObj.elements.length - 1));
    renderGroupList();
    renderGroupEditor();
  }

  // numOrNull reads a number input; "" (untouched) means "keep the server value".
  function numOrNull(id) {
    const el = $(id);
    if (!el) return null;
    const s = el.value.trim();
    if (s === "") return null;
    const v = parseFloat(s);
    return isNaN(v) ? null : v;
  }

  // collectSelected applies the editor of the selected element to a copy of
  // groupObj: kind / glass / gap_glass / fold (and a rear surface invented when
  // a lens has none) go into the JSON itself, the numbers become the sets of
  // /api/lens/group/set. Only the selected element is sent.
  function collectSelected() {
    const i = groupSel;
    const g = JSON.parse(JSON.stringify(groupObj));
    const e = g.elements[i];
    const sets = [];
    const push = function (param, value) {
      if (value !== null && value !== undefined) sets.push({ param: param, value: value });
    };
    const kind = $("ozGEKind").value;
    e.kind = kind;
    if (kind === "lens") {
      if (!e.rear) {
        const r2 = numOrNull("ozGER2");
        e.rear = { radius_mm: r2 === null ? 0 : r2 };
        if (e.front && e.front.semi_aperture_mm) e.rear.semi_aperture_mm = e.front.semi_aperture_mm;
      }
      const glass = $("ozGEGlass").value.trim();
      if (glass) e.glass = glass; else delete e.glass;
    } else {
      delete e.glass;
      delete e.rear;
      e.thickness_mm = 0;
    }
    const gapGlass = $("ozGEGapGlass") ? $("ozGEGapGlass").value.trim() : "";
    if (kind !== "stop" && gapGlass) e.gap_glass = gapGlass; else delete e.gap_glass;
    e.fold = ($("ozGEFold").value === "off") ? 1 : 0;
    if (e.fold === 0) delete e.fold;

    push("e" + i + ".front.R", numOrNull("ozGER1"));
    push("e" + i + ".front.K", numOrNull("ozGEK1"));
    push("e" + i + ".front.semi", numOrNull("ozGESemi"));
    if (e.rear) {
      push("e" + i + ".rear.R", numOrNull("ozGER2"));
      push("e" + i + ".rear.K", numOrNull("ozGEK2"));
      push("e" + i + ".rear.semi", numOrNull("ozGESemi2"));
    }
    if (kind === "lens") push("e" + i + ".thk", numOrNull("ozGEThk"));
    push("e" + i + ".gap", numOrNull("ozGEGap"));
    push("e" + i + ".dx", numOrNull("ozGEDx"));
    push("e" + i + ".dy", numOrNull("ozGEDy"));
    push("e" + i + ".tx", numOrNull("ozGETx"));
    push("e" + i + ".ty", numOrNull("ozGETy"));
    push("e" + i + ".tz", numOrNull("ozGETz"));
    return { group: g, sets: sets };
  }

  function fx(v) {
    return (typeof v === "number" && isFinite(v)) ? v.toPrecision(6) : String(v);
  }

  function fv(p) {
    return "(" + fx(p.x) + ", " + fx(p.y) + ", " + fx(p.z) + ")";
  }

  // objectIsFinite reports whether the group images a finite object (its
  // object_z_mm is a number, not "inf").
  function objectIsFinite() {
    return !!groupObj && typeof groupObj.object_z_mm === "number" && isFinite(groupObj.object_z_mm);
  }

  // groupBody is the shared body of the group measurements: the group itself
  // plus the bundle parameters of the panel. An object at infinity carries a
  // field angle (field_deg); a finite object carries its height in mm
  // (field_tan) — the API rejects field_deg for a finite object on purpose.
  function groupBody(extra) {
    const body = {
      group: groupObj,
      lambda_nm: num("ozGroupLambda"),
      pupil_radius_mm: num("ozGroupPupil"),
    };
    if (objectIsFinite()) {
      body.field_tan = num("ozObjH");
    } else {
      body.field_deg = num("ozObjFieldDeg");
    }
    for (const k in extra) {
      body[k] = extra[k];
    }
    return body;
  }

  // renderObjectBlock writes the 物 panel from groupObj (∞ / 有限距离) and
  // toggles which size field applies: 视场角 (°) for an object at infinity,
  // 物高 (mm) for a finite one.
  function renderObjectBlock() {
    const mode = $("ozObjMode");
    if (!mode) return;
    const finite = objectIsFinite();
    mode.value = finite ? "near" : "inf";
    const dist = $("ozObjDist");
    if (finite && dist) dist.value = String(Math.abs(groupObj.object_z_mm));
    $("ozObjDistRow").hidden = !finite;
    $("ozObjHRow").hidden = !finite;
    $("ozObjDegRow").hidden = finite;
  }

  // applyObjectMode pushes the 物 panel into groupObj: "inf" keeps the design
  // object at infinity, "near" stores the object plane as a negative z (in
  // front of the first surface, the convention of lens.Group).
  function applyObjectMode() {
    if (!groupObj) return;
    const mode = $("ozObjMode").value;
    if (mode === "near") {
      const d = Math.abs(num("ozObjDist")) || 500;
      groupObj.object_z_mm = -d;
    } else {
      groupObj.object_z_mm = "inf";
    }
    renderObjectBlock();
  }

  async function groupLoad() {
    const idx = parseInt($("ozGroupExample").value, 10);
    const ex = GROUP_EXAMPLES[isNaN(idx) ? 0 : idx];
    groupObj = JSON.parse(JSON.stringify(ex.group));
    groupLayout = null;
    groupSel = 0;
    renderObjectBlock();
    renderGroupList();
    renderGroupEditor();
    setGroupButtons(true);
    setGroupStatus(ex.label + " · " + groupObj.elements.length + " 元件");
    show("ozGroupOut", "已载入「" + ex.label + "」：" + groupObj.elements.length + " 个元件。", false);
    await groupLayoutCall();
  }

  async function groupApply() {
    if (!groupObj) { show("ozGroupOut", "先「载入示例」。", true); return; }
    const packed = collectSelected();
    const d = await api("/api/lens/group/set", packed);
    groupObj = d.group;
    groupLayout = null;
    renderGroupList();
    renderGroupEditor();
    await groupLayoutCall();
    show("ozGroupOut", "已应用 " + d.applied + " 项编辑（第 " + groupSel + " 个元件）。\n" + groupSummary(), false);
  }

  // fixStop points groupObj.stop at the first stop element (or drops the field
  // when the group has none, which the server then judges).
  function fixStop() {
    let idx = -1;
    groupObj.elements.forEach(function (e, i) {
      if (e.kind === "stop" && idx < 0) idx = i;
    });
    if (idx < 0) delete groupObj.stop;
    else groupObj.stop = idx;
  }

  // structural runs a structural edit of the optical path (insert / duplicate /
  // delete). The group is client-held, so the edit is applied locally and then
  // validated by one layout call; if the server rejects it the previous group
  // comes back untouched and the error text is shown as usual.
  async function structural(label, mutate) {
    if (!groupObj) { show("ozGroupOut", "先「载入示例」。", true); return; }
    const prev = JSON.parse(JSON.stringify(groupObj));
    const prevSel = groupSel;
    mutate();
    fixStop();
    if (!groupObj.elements.length) {
      groupObj = prev;
      groupSel = prevSel;
      show("ozGroupOut", "至少要保留一个元件。", true);
      return;
    }
    groupSel = Math.max(0, Math.min(groupSel, groupObj.elements.length - 1));
    renderGroupList();
    renderGroupEditor();
    try {
      await groupLayoutCall();
    } catch (err) {
      groupObj = prev;
      groupSel = prevSel;
      renderGroupList();
      renderGroupEditor();
      throw new Error(label + "失败：" + err.message);
    }
  }

  // newElement is the default element of each kind: a plano-convex-ish lens,
  // a concave mirror or a flat stop, each with a sane gap.
  function newElement(kind) {
    if (kind === "stop") {
      return { kind: "stop", front: { radius_mm: 0, semi_aperture_mm: 10 }, thickness_mm: 0, gap_mm: 20 };
    }
    if (kind === "mirror") {
      return { kind: "mirror", front: { radius_mm: -400, semi_aperture_mm: 10 }, thickness_mm: 0, gap_mm: 40 };
    }
    return {
      kind: "lens",
      front: { radius_mm: 100, semi_aperture_mm: 10 },
      rear: { radius_mm: -100, semi_aperture_mm: 10 },
      glass: "N-BK7",
      thickness_mm: 4,
      gap_mm: 20,
    };
  }

  function groupInsert() {
    const kind = $("ozGroupAddKind").value;
    return structural("新增元件", function () {
      groupObj.elements.splice(groupSel + 1, 0, newElement(kind));
      groupSel = groupSel + 1;
    });
  }

  function groupDup() {
    return structural("复制元件", function () {
      const src = JSON.parse(JSON.stringify(groupObj.elements[groupSel]));
      if (src.name) src.name = src.name + " copy";
      groupObj.elements.splice(groupSel + 1, 0, src);
      groupSel = groupSel + 1;
    });
  }

  function groupDel() {
    return structural("删除元件", function () {
      groupObj.elements.splice(groupSel, 1);
    });
  }

  function groupSummary() {
    const lines = [];
    groupObj.elements.forEach(function (e, i) {
      let s = "[" + i + "] " + e.kind;
      if (e.front && e.front.radius_mm) s += " R1=" + fx(e.front.radius_mm);
      if (e.rear) s += " R2=" + fx(e.rear.radius_mm);
      if (e.thickness_mm) s += " thk=" + fx(e.thickness_mm);
      s += " gap=" + fx(e.gap_mm);
      if (e.glass) s += " glass=" + e.glass;
      if (e.tilt_x_deg) s += " tx=" + e.tilt_x_deg + "°";
      if (e.tilt_y_deg) s += " ty=" + e.tilt_y_deg + "°";
      lines.push(s);
    });
    return lines.join("\n");
  }

  async function groupLayoutCall() {
    if (!groupObj) { show("ozGroupOut", "先「载入示例」。", true); return; }
    const d = await api("/api/lens/group/layout", groupBody());
    groupLayout = d;
    drawGroupPlot();
    setGroupStatus(groupObj.elements.length + " 元件 · " + d.num_surfaces + " 面 · 折转 " + d.folds + " 次");
    const lines = ["布局：" + d.num_surfaces + " 面，折转 " + d.folds + " 次"];
    for (const s of d.surfaces) {
      lines.push("  [e" + s.element + " " + s.face + "] 面 " + s.surface + "  顶点 " + fv(s.vertex) +
        " mm  轴向 " + fv(s.axis) + "  折转 " + s.folds);
    }
    lines.push("  出射：顶点 " + fv(d.exit_vertex) + " mm，光束方向 " + fv(d.exit_axis));
    show("ozGroupOut", lines.join("\n"), false);
  }

  async function groupImageCall() {
    if (!groupObj) { show("ozGroupOut", "先「载入示例」。", true); return; }
    const d = await api("/api/lens/group/image", groupBody());
    show("ozGroupOut", [
      "成像位置 · λ " + fx(d.lambda_nm) + " nm · 物 " +
      (d.object_z_mm === null ? "∞" : fx(d.object_z_mm) + " mm") +
      " · 场值 " + fx(d.field_tan_y) + " · " + d.samples + " 环采样",
      "像点 " + fv(d.point) + " mm · RMS " + fx(d.rms_mm) + " mm · 最大半径 " + fx(d.geo_radius_mm) + " mm",
      "光线 " + d.rays + "（被拦 " + d.failed + "）",
      "焦面：距离末面顶点 " + fx(d.focus_mm) + " mm，法向 " + fv(d.axis),
      "主光线落点 " + fv(d.chief_hit) + " mm，方向 " + fv(d.chief_dir),
    ].join("\n"), false);
  }

  async function groupAfocalCall() {
    if (!groupObj) { show("ozGroupOut", "先「载入示例」。", true); return; }
    const d = await api("/api/lens/group/afocal", groupBody());
    show("ozGroupOut", [
      "角放大率 · λ " + fx(d.lambda_nm) + " nm · " + d.rays + " 条光线",
      "Γ = " + fx(d.gamma) + " · 入/出射束直径比 " + fx(d.beam_ratio) + " · 准直残差 " + fx(d.spread_rad) + " rad",
      "入射轴 " + fv(d.axis_in) + " → 出射轴 " + fv(d.axis_out),
    ].join("\n"), false);
  }

  async function groupDispersionCall() {
    if (!groupObj) { show("ozGroupOut", "先「载入示例」。", true); return; }
    const d = await api("/api/lens/group/dispersion", groupBody());
    show("ozGroupOut", [
      "色散（F–C）· λa " + fx(d.lambda_a_nm) + " nm → 焦面 " + fx(d.a.focus_mm) + " mm",
      "                λb " + fx(d.lambda_b_nm) + " nm → 焦面 " + fx(d.b.focus_mm) + " mm",
      "焦移 (a − b) = " + fx(d.shift_mm) + " mm（沿出射轴，从末面顶点量起）",
    ].join("\n"), false);
  }

  async function groupPrescriptionCall() {
    if (!groupObj) { show("ozGroupOut", "先「载入示例」。", true); return; }
    const d = await api("/api/lens/group/prescription", groupBody());
    show("ozGroupOut", "处方（卡片链）：\n" + JSON.stringify(d.prescription, null, 1), false);
  }

  /* --- 3D 布局画布：正交投影，拖拽旋转（yaw 绕 y、pitch 绕 x） --- */

  function groupRotate(p) {
    const cy = Math.cos(groupView.yaw), sy = Math.sin(groupView.yaw);
    const cx = Math.cos(groupView.pitch), sx = Math.sin(groupView.pitch);
    const x1 = p.x * cy + p.z * sy;
    const z1 = -p.x * sy + p.z * cy;
    const y1 = p.y * cx - z1 * sx;
    return { x: x1, y: y1 };
  }

  function drawGroupPlot() {
    const cv = $("ozGroupPlot");
    if (!cv) return;
    const dpr = window.devicePixelRatio || 1;
    const w = Math.max(1, cv.clientWidth), h = Math.max(1, cv.clientHeight);
    cv.width = Math.round(w * dpr);
    cv.height = Math.round(h * dpr);
    const c = cv.getContext("2d");
    c.setTransform(dpr, 0, 0, dpr, 0, 0);
    c.clearRect(0, 0, w, h);
    c.font = "11px system-ui,sans-serif";
    const hint = $("ozGroupViewHint");
    if (!groupLayout || !groupLayout.surfaces || !groupLayout.surfaces.length) {
      c.fillStyle = "#5a7d9e";
      c.textAlign = "center";
      c.fillText("点「载入示例」再「布局」绘制 3D 光轴", w / 2, h / 2);
      if (hint) hint.textContent = "拖拽画布旋转视角";
      return;
    }
    const pts = [];
    for (const s of groupLayout.surfaces) pts.push(s.vertex);
    pts.push(groupLayout.exit_vertex);
    // 有限距离的物：按真实物高在物面上画一支箭头（位置 = object_z_mm，
    // 长度 = 物高 mm）。∞ 物没有可画的长度，只在状态栏标视场角。
    const objPts = [];
    if (objectIsFinite()) {
      const z0 = groupObj.object_z_mm;
      const h0 = num("ozObjH");
      objPts.push({ x: 0, y: 0, z: z0 }, { x: 0, y: h0, z: z0 });
      pts.push(objPts[0], objPts[1]);
    }
    // 世界包围盒（决定轴线/法向的线段长度）
    let wlo = { x: Infinity, y: Infinity, z: Infinity };
    let whi = { x: -Infinity, y: -Infinity, z: -Infinity };
    for (const p of pts) {
      wlo.x = Math.min(wlo.x, p.x); whi.x = Math.max(whi.x, p.x);
      wlo.y = Math.min(wlo.y, p.y); whi.y = Math.max(whi.y, p.y);
      wlo.z = Math.min(wlo.z, p.z); whi.z = Math.max(whi.z, p.z);
    }
    const worldSpan = Math.max(whi.x - wlo.x, whi.y - wlo.y, whi.z - wlo.z, 1e-6);
    // 视图包围盒（缩放到画布）
    let lo = { x: Infinity, y: Infinity }, hi = { x: -Infinity, y: -Infinity };
    const view = pts.map(groupRotate);
    for (const v of view) {
      lo.x = Math.min(lo.x, v.x); lo.y = Math.min(lo.y, v.y);
      hi.x = Math.max(hi.x, v.x); hi.y = Math.max(hi.y, v.y);
    }
    const cvx = (lo.x + hi.x) / 2, cvy = (lo.y + hi.y) / 2;
    const spanX = hi.x - lo.x, spanY = hi.y - lo.y;
    let scale = 1;
    if (spanX > 1e-9 || spanY > 1e-9) {
      scale = Math.min((w * 0.78) / Math.max(spanX, 1e-9), (h * 0.72) / Math.max(spanY, 1e-9));
    }
    const map = function (p) {
      const v = groupRotate(p);
      return { x: w / 2 + (v.x - cvx) * scale, y: h / 2 - (v.y - cvy) * scale };
    };
    const seg = function (a, b) {
      const p = map(a), q = map(b);
      c.beginPath(); c.moveTo(p.x, p.y); c.lineTo(q.x, q.y); c.stroke();
    };
    // 1) 链（顶点连线）
    c.strokeStyle = "#7d97ad";
    c.lineWidth = 1;
    for (let i = 1; i < groupLayout.surfaces.length; i++) {
      seg(groupLayout.surfaces[i - 1].vertex, groupLayout.surfaces[i].vertex);
    }
    // 2) 逐面：轴线（蓝，实线）与顶点法向（橙，虚线）
    for (const s of groupLayout.surfaces) {
      const v = s.vertex;
      c.strokeStyle = "#1976d2";
      c.lineWidth = 1.5;
      seg(v, { x: v.x + s.axis.x * worldSpan * 0.06, y: v.y + s.axis.y * worldSpan * 0.06, z: v.z + s.axis.z * worldSpan * 0.06 });
      c.strokeStyle = "#ff922b";
      c.lineWidth = 1;
      if (c.setLineDash) c.setLineDash([3, 3]);
      seg(v, { x: v.x + s.normal.x * worldSpan * 0.05, y: v.y + s.normal.y * worldSpan * 0.05, z: v.z + s.normal.z * worldSpan * 0.05 });
      if (c.setLineDash) c.setLineDash([]);
      const q = map(v);
      c.fillStyle = s.folds > 0 ? "#e8590c" : "#1976d2";
      c.beginPath();
      c.arc(q.x, q.y, s.face === "front" ? 3 : 2, 0, 2 * Math.PI);
      c.fill();
      if (s.face === "front") {
        c.fillStyle = "#5a7d9e";
        c.textAlign = "left";
        c.fillText("e" + s.element, q.x + 5, q.y - 5);
      }
    }
    // 3) 出射光束方向（虚线）
    const ev = groupLayout.exit_vertex, ea = groupLayout.exit_axis;
    c.strokeStyle = "#2f9e44";
    c.lineWidth = 1.5;
    if (c.setLineDash) c.setLineDash([5, 4]);
    seg(ev, { x: ev.x + ea.x * worldSpan * 0.2, y: ev.y + ea.y * worldSpan * 0.2, z: ev.z + ea.z * worldSpan * 0.2 });
    if (c.setLineDash) c.setLineDash([]);
    // 4) 坐标三向标（画在第一个顶点）
    const o = groupLayout.surfaces[0].vertex;
    const triad = [[{ x: 1, y: 0, z: 0 }, "x"], [{ x: 0, y: 1, z: 0 }, "y"], [{ x: 0, y: 0, z: 1 }, "z"]];
    c.strokeStyle = "#adb5bd";
    c.lineWidth = 1;
    for (const t of triad) {
      seg(o, { x: o.x + t[0].x * worldSpan * 0.08, y: o.y + t[0].y * worldSpan * 0.08, z: o.z + t[0].z * worldSpan * 0.08 });
      const q = map({ x: o.x + t[0].x * worldSpan * 0.11, y: o.y + t[0].y * worldSpan * 0.11, z: o.z + t[0].z * worldSpan * 0.11 });
      c.fillStyle = "#868e96";
      c.textAlign = "center";
      c.fillText(t[1], q.x, q.y - 2);
    }
    // 5) 物（有限距离时按真实物高画箭头；物高 0 就是一个物点）
    if (objPts.length) {
      const a = map(objPts[0]), b = map(objPts[1]);
      c.strokeStyle = "#c92a2a";
      c.lineWidth = 2;
      c.beginPath(); c.moveTo(a.x, a.y); c.lineTo(b.x, b.y); c.stroke();
      // 箭头尖（物高非 0 时）
      if (Math.hypot(b.x - a.x, b.y - a.y) > 6) {
        const ang = Math.atan2(b.y - a.y, b.x - a.x);
        const L = 8;
        c.beginPath();
        c.moveTo(b.x, b.y);
        c.lineTo(b.x - L * Math.cos(ang - 0.4), b.y - L * Math.sin(ang - 0.4));
        c.moveTo(b.x, b.y);
        c.lineTo(b.x - L * Math.cos(ang + 0.4), b.y - L * Math.sin(ang + 0.4));
        c.stroke();
      }
      c.fillStyle = "#c92a2a";
      c.beginPath(); c.arc(a.x, a.y, 3, 0, 2 * Math.PI); c.fill();
      c.textAlign = "left";
      c.fillText("物", a.x + 6, a.y + 4);
    }
    c.fillStyle = "#5a7d9e";
    c.textAlign = "right";
    const viewLabel = viewName();
    c.fillText(viewLabel, w - 6, 14);
    if (hint) {
      hint.textContent = viewLabel + " · 拖拽画布旋转 · " +
        (objectIsFinite()
          ? "物 " + fx(Math.abs(groupObj.object_z_mm)) + " mm，物高 " + fx(num("ozObjH")) + " mm"
          : "物 ∞，视场角 " + fx(num("ozObjFieldDeg")) + "°");
    }
  }

  function viewName() {
    const yaw = groupView.yaw, pitch = groupView.pitch;
    if (Math.abs(yaw) < 1e-6 && Math.abs(pitch) < 1e-6) return "XY 面（沿 z 看）";
    if (Math.abs(yaw) < 1e-6 && Math.abs(pitch + Math.PI / 2) < 1e-6) return "XZ 面（沿 y 看）";
    return "等轴测";
  }

  function bindGroupCanvas() {
    const cv = $("ozGroupPlot");
    if (!cv) return;
    let dragging = false, lx = 0, ly = 0;
    cv.addEventListener("pointerdown", function (ev) {
      dragging = true;
      lx = ev.clientX;
      ly = ev.clientY;
      if (cv.setPointerCapture) {
        try { cv.setPointerCapture(ev.pointerId); } catch (e) { /* 触摸事件可能不支持捕获 */ }
      }
      ev.preventDefault();
    });
    cv.addEventListener("pointermove", function (ev) {
      if (!dragging) return;
      groupView.yaw += (ev.clientX - lx) * 0.012;
      groupView.pitch += (ev.clientY - ly) * 0.012;
      groupView.pitch = Math.max(-1.5, Math.min(1.5, groupView.pitch));
      lx = ev.clientX;
      ly = ev.clientY;
      drawGroupPlot();
    });
    const stop = function () { dragging = false; };
    cv.addEventListener("pointerup", stop);
    cv.addEventListener("pointercancel", stop);
    cv.addEventListener("pointerleave", stop);
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
    setAdvanced(!!prescription);
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

  /* --- 几何光学子模式（下拉）------------------------------------------------
   * 「透镜组」是默认子模式：index.html 的 #ozModeSel 里 selected 的选项，
   * 左栏只显示 #ozPanelGroup，右栏只显示 #ozCardGroup。其它三个子模式
   * （消色差设计 / 公差分析 / 膜系设计）从同一个下拉切换。选择记在
   * localStorage（键 oz.submode），与主题、立体视图拖动方向同一套做法：
   * 隐私模式下写入失败就静默退回默认值。 */

  const SUBMODE_KEY = "oz.submode";
  const SUBMODES = {
    group: { panels: ["ozPanelGroup"], cards: ["ozCardGroup"] },
    achromat: { panels: ["ozPanelAchromat"], cards: ["ozCardAchromat"] },
    tol: { panels: ["ozPanelTol"], cards: ["ozCardTol"] },
    coat: { panels: ["ozPanelCoat"], cards: ["ozCardCoat"] },
  };
  // 子模式循环顺序（k 键）：透镜组 → 消色差 → 公差 → 膜系 → 透镜组。
  const SUBMODE_NEXT = { group: "achromat", achromat: "tol", tol: "coat", coat: "group" };

  function currentSubmode() {
    const sel = $("ozModeSel");
    const v = sel ? sel.value : "group";
    return SUBMODES[v] ? v : "group";
  }

  function setSubmode(name) {
    if (!SUBMODES[name]) name = "group";
    const sel = $("ozModeSel");
    if (sel && sel.value !== name) sel.value = name;
    for (const key of Object.keys(SUBMODES)) {
      const on = key === name;
      for (const id of SUBMODES[key].panels) if ($(id)) $(id).hidden = !on;
      for (const id of SUBMODES[key].cards) if ($(id)) $(id).hidden = !on;
    }
    // 画布只有在可见时才有尺寸：切回「透镜组」立刻补一次重绘。
    if (name === "group") drawGroupPlot();
    if (name === "coat") drawPlot();
    try { localStorage.setItem(SUBMODE_KEY, name); } catch (e) { /* 隐私模式 */ }
    // 让 app.js（与 AI 助手）看到当前子模式。
    document.dispatchEvent(new CustomEvent("oz:submode", { detail: name }));
  }

  function loadSubmode() {
    let v = "group";
    try { v = localStorage.getItem(SUBMODE_KEY) || "group"; } catch (e) { v = "group"; }
    setSubmode(SUBMODES[v] ? v : "group");
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
    buildAdvanced();
    buildGroup();
    buildGroupView();
    guard("ozDesignBtn", "ozLensOut", design);
    guard("ozEvalBtn", "ozLensOut", evaluate);
    guard("ozTolBtn", "ozTolOut", tolerance);
    guard("ozAdvThermalBtn", "ozAdvOut", thermal);
    guard("ozAdvSensBtn", "ozAdvOut", sensitivity);
    guard("ozGroupLoadBtn", "ozGroupOut", groupLoad);
    guard("ozGroupApplyBtn", "ozGroupOut", groupApply);
    guard("ozGroupAddBtn", "ozGroupOut", groupInsert);
    guard("ozGroupDupBtn", "ozGroupOut", groupDup);
    guard("ozGroupDelBtn", "ozGroupOut", groupDel);
    guard("ozGroupLayoutBtn", "ozGroupOut", groupLayoutCall);
    guard("ozGroupImageBtn", "ozGroupOut", groupImageCall);
    guard("ozGroupAfocalBtn", "ozGroupOut", groupAfocalCall);
    guard("ozGroupDispBtn", "ozGroupOut", groupDispersionCall);
    guard("ozGroupPrescBtn", "ozGroupOut", groupPrescriptionCall);
    guard("ozCoatBtn", "ozCoatOut", designCoating);
    guard("ozSpecBtn", "ozCoatOut", spectrum);
    // The 3D layout canvas: view presets, drag-to-rotate, redraw on layout.
    const setView = function (yaw, pitch) {
      groupView.yaw = yaw;
      groupView.pitch = pitch;
      drawGroupPlot();
    };
    $("ozGroupIsoBtn").addEventListener("click", function () { setView(-0.62, 0.45); });
    $("ozGroupXZBtn").addEventListener("click", function () { setView(0, -Math.PI / 2); });
    $("ozGroupXYBtn").addEventListener("click", function () { setView(0, 0); });
    bindGroupCanvas();
    // 几何光学子模式：下拉切换（透镜组为默认），k 键循环。
    const sel = $("ozModeSel");
    if (sel) sel.addEventListener("change", function () { setSubmode(sel.value); });
    loadSubmode();
    // 物：∞ ↔ 有限距离的切换只改 group 的 object_z_mm，测量时按物距选
    // field_deg（∞）或 field_tan=物高（有限）。改完补一次画布重绘，让物
    // 箭头的位置/长度跟着变。
    $("ozObjMode").addEventListener("change", function () { applyObjectMode(); drawGroupPlot(); });
    $("ozObjDist").addEventListener("change", function () {
      if (objectIsFinite()) { applyObjectMode(); drawGroupPlot(); }
    });
    $("ozObjH").addEventListener("change", function () { drawGroupPlot(); });
    $("ozObjFieldDeg").addEventListener("change", function () { drawGroupPlot(); });
    // 元件列表的 ↑/↓（与波动光学的元件列表同习惯）。
    $("ozGroupList").addEventListener("keydown", function (ev) {
      if (ev.key === "ArrowDown") { ev.preventDefault(); selectGroupEl(groupSel + 1); }
      else if (ev.key === "ArrowUp") { ev.preventDefault(); selectGroupEl(groupSel - 1); }
    });
    window.addEventListener("resize", drawPlot);
    window.addEventListener("orientationchange", drawPlot);
    window.addEventListener("resize", drawGroupPlot);
    window.addEventListener("orientationchange", drawGroupPlot);
    if (window.ResizeObserver) new ResizeObserver(drawPlot).observe($("ozPlot"));
    if (window.ResizeObserver) new ResizeObserver(drawGroupPlot).observe($("ozGroupPlot"));
    // app.js announces mode changes; the canvases need a layout pass to be sized.
    document.addEventListener("oz:mode", function (ev) {
      if (ev.detail === "oz") {
        drawPlot();
        drawGroupPlot();
      }
    });
    // k 键（app.js）与下拉之外的入口都通过这个事件切换子模式；setSubmode
    // 自己也会发同名事件，currentSubmode() 判等后不会重入。
    document.addEventListener("oz:submode", function (ev) {
      if (ev.detail && ev.detail !== currentSubmode()) setSubmode(ev.detail);
    });
    // 「透镜组」是默认子模式：进页面就有内容可看，自动载入第一个示例
    // （失败不阻塞初始化，结果区会显示服务端的错误文本）。
    groupLoad().catch(function (e) { show("ozGroupOut", "错误：" + e.message, true); });
    loadVersion();
    drawPlot();
    drawGroupPlot();
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init);
  } else {
    init();
  }
})();
