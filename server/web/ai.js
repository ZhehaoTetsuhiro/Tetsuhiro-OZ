/* Tetsuhiro-OZ AI 助手（右侧栏，默认隐藏）。
 *
 * 独立 IIFE：不向全局作用域写入任何名字（与 app.js 的全局 S、oz.js 的私有 $
 * 互不干扰）。本文件只做三件事：
 *   1. 让用户配置自己的 AI 提供方（base_url / model / api_key），存在浏览器
 *      localStorage 里——密钥不落服务端、不进日志、只随每次请求转发给用户
 *      自己配置的提供方；
 *   2. 把「规范 + 当前状态」打包交给服务端 /api/ai/chat（服务端注入规范摘要），
 *      因此助手看到的元件目录永远是本二进制真正实现的那些；
 *   3. 解析助手回复里的 ```oz-config / ```oz-element 代码块，给出「应用到场景」
 *      与「写入元件定义」按钮——配置走 app.js 的 adoptConfig（与打开文件同一条
 *      校验+运行路径），元件定义走 /api/ai/element（内核校验通过才落盘）。
 *
 * 状态读取全部是只读的：S.config / S.meta / curComp() / curSrc() 来自 app.js，
 * 几何与膜系结果取自 oz.js 渲染好的 DOM 文本，两边都不需要为 AI 改一行代码。
 */
(function () {
  "use strict";

  const $ = (id) => document.getElementById(id);
  const SETTINGS_KEY = "oz.ai.settings";

  // 提供方预设：base_url 填到 /v1 根即可，服务端会补 /chat/completions。
  const PRESETS = [
    { id: "custom", label: "自定义（OpenAI 兼容）", base: "", model: "" },
    { id: "openai", label: "OpenAI", base: "https://api.openai.com/v1", model: "gpt-4o-mini" },
    { id: "deepseek", label: "DeepSeek", base: "https://api.deepseek.com/v1", model: "deepseek-chat" },
    { id: "qwen", label: "通义千问（百炼兼容模式）", base: "https://dashscope.aliyuncs.com/compatible-mode/v1", model: "qwen-plus" },
    { id: "moonshot", label: "Moonshot（Kimi）", base: "https://api.moonshot.cn/v1", model: "moonshot-v1-8k" },
    { id: "glm", label: "智谱 GLM", base: "https://open.bigmodel.cn/api/paas/v4", model: "glm-4-flash" },
    { id: "llama", label: "本地 llama.cpp", base: "http://127.0.0.1:8080/v1", model: "qwen3" },
    { id: "ollama", label: "本地 Ollama", base: "http://127.0.0.1:11434/v1", model: "qwen2.5" },
  ];

  const state = {
    open: false,
    settings: loadSettings(),
    messages: [],       // {role, content} —— 只保留用户与助手的对话
    busy: false,
    abort: null,        // 当前请求的 AbortController（「停止」按钮用）
    specChars: 0,
    specSections: [],
  };

  function loadSettings() {
    const def = { preset: "custom", base_url: "", model: "", api_key: "", temperature: 0.3 };
    try {
      const raw = localStorage.getItem(SETTINGS_KEY);
      if (!raw) return def;
      const s = JSON.parse(raw);
      return Object.assign(def, {
        preset: String(s.preset || "custom"),
        base_url: String(s.base_url || ""),
        model: String(s.model || ""),
        api_key: String(s.api_key || ""),
        temperature: Number(s.temperature) > 0 ? Number(s.temperature) : 0.3,
      });
    } catch (e) { return def; }
  }

  function saveSettings() {
    try { localStorage.setItem(SETTINGS_KEY, JSON.stringify(state.settings)); } catch (e) { /* 隐私模式：忽略 */ }
  }

  // ---------- 小工具 ----------
  function esc(s) {
    return String(s).replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;").replace(/"/g, "&quot;");
  }

  function setStatus(msg, isErr) {
    const el = $("aiStatus");
    if (!el) return;
    el.textContent = msg;
    el.classList.toggle("err", !!isErr);
  }

  async function api(path, body) {
    const r = await fetch(path, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(body),
      signal: state.abort ? state.abort.signal : undefined,
    });
    const text = await r.text();
    let d;
    try { d = JSON.parse(text); } catch (e) { d = { raw: text }; }
    if (!r.ok) throw new Error(d && d.error ? d.error : "HTTP " + r.status);
    return d;
  }

  // ---------- 当前状态（只读收集） ----------
  // 助手要「自动读取规范、当前状态」：规范由服务端注入，当前状态在这里收集。
  function collectState() {
    const st = { mode: typeof S !== "undefined" ? S.mode : "wave" };
    try {
      if (st.mode === "wave" && typeof S !== "undefined" && S.config) {
        st.config = S.config;
        const c = typeof curComp === "function" ? curComp() : null;
        const src = typeof curSrc === "function" ? curSrc() : null;
        if (c) {
          const doc = typeof docFor === "function" ? docFor(c.type) : null;
          st.selected = { kind: "component", type: c.type, label: c.label || (doc && doc.label) || c.type,
            help: doc ? doc.help : "", params: c.params || {} };
        } else if (src) {
          const doc = typeof srcDocFor === "function" ? srcDocFor(src.type) : null;
          st.selected = { kind: "source", type: src.type, label: src.label || (doc && doc.label) || src.type,
            help: doc ? doc.help : "", params: src.params || {} };
        }
        if (S.meta) {
          st.run = {
            status: S.meta.status || "done",
            elapsed_ms: S.meta.elapsed_ms,
            grid: S.meta.grid ? S.meta.grid.size : 0,
            planes: (S.meta.planes || []).map((p) => (p.label || p.id) + (p.path ? " [" + p.path + "]" : "")),
            warnings: S.meta.warnings || [],
          };
        }
        const warnBox = $("warnings");
        if (warnBox && warnBox.textContent.trim()) st.warnings = [warnBox.textContent.trim().slice(0, 800)];
      }
      if (st.mode === "quantum" && typeof S !== "undefined" && S.qconfig) {
        st.quantum = JSON.stringify(S.qconfig);
      }
      // 几何/膜系模式的结果在 oz.js 的 DOM 里（只读取文本）。
      if (st.mode === "oz") {
        const lens = $("ozLensOut"), coat = $("ozCoatOut");
        if (lens && lens.textContent.trim()) st.lens = lens.textContent.trim().slice(0, 3000);
        if (coat && coat.textContent.trim()) st.coating = coat.textContent.trim().slice(0, 3000);
      }
    } catch (e) { /* 状态收集失败不阻断对话 */ }
    return st;
  }

  // ---------- 渲染 ----------
  // 助手回复的轻量渲染：先整体转义，再还原代码围栏；```oz-config / ```oz-element
  // 额外挂上动作按钮。不用第三方 markdown，也就没有注入面。
  function renderContent(text) {
    const parts = String(text).split(/```([a-zA-Z-]*)\n?/);
    let html = "";
    for (let i = 0; i < parts.length; i += 2) {
      html += renderText(parts[i]);
      const lang = (parts[i + 1] || "").trim();
      const code = parts[i + 2] === undefined ? "" : parts[i + 2].replace(/\n$/, "");
      if (lang === "") break;
      const kind = lang === "oz-config" ? "config" : lang === "oz-element" ? "element" : "";
      html += '<pre class="aiCode">' +
        (kind ? '<button class="aiApply" data-kind="' + kind + '" data-code="' + esc(code) + '">' +
          (kind === "config" ? "应用到场景" : "写入元件定义") + "</button>" : "") +
        '<code>' + esc(code) + "</code></pre>";
    }
    return html;
  }

  function renderText(t) {
    let s = esc(t);
    s = s.replace(/`([^`\n]+)`/g, "<code>$1</code>");
    s = s.split(/\n{2,}/).map((p) => "<p>" + p.replace(/\n/g, "<br>") + "</p>").join("");
    return s;
  }

  function renderLog() {
    const box = $("aiLog");
    box.innerHTML = "";
    if (!state.messages.length) {
      box.innerHTML = '<p class="hint">还没有对话。可以问「这个元件是做什么的？」「当前配置有什么问题？」「设计一个杨氏双缝干涉光路」，也可以直接说「新建一个轴锥镜元件」。</p>';
    }
    state.messages.forEach((m, i) => {
      const div = document.createElement("div");
      div.className = "aiMsg " + (m.role === "user" ? "me" : "ai");
      if (m.role === "user") {
        div.textContent = m.content;
      } else {
        div.innerHTML = renderContent(m.content);
        if (m.pending) div.innerHTML += '<span class="aiCursor" aria-hidden="true">▍</span>';
      }
      div.dataset.idx = String(i);
      box.appendChild(div);
    });
    box.scrollTop = box.scrollHeight;
    // 重新绑定动作按钮（innerHTML 会丢掉监听器）。
    box.querySelectorAll(".aiApply").forEach((b) => {
      b.addEventListener("click", () => onApply(b.dataset.kind, b.dataset.code, b));
    });
  }

  // ---------- 动作：应用配置 / 写入元件 ----------
  async function onApply(kind, code, btn) {
    btn.disabled = true;
    try {
      const obj = JSON.parse(code);
      if (kind === "config") {
        if (typeof adoptConfig !== "function") throw new Error("页面尚未就绪");
        setStatus("正在校验并应用配置…");
        const v = await fetch("/api/validate", {
          method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(obj),
        }).then((r) => r.json());
        if (v && v.ok === false) {
          const first = (v.issues || [])[0];
          throw new Error("配置校验未通过" + (first ? "：" + first.path + " " + first.message : ""));
        }
        await adoptConfig(obj, "AI 助手：已应用配置");
        setStatus("配置已应用");
      } else {
        setStatus("正在写入元件定义…");
        const d = await api("/api/ai/element", { definition: obj });
        setStatus("已写入 " + d.path);
        if (typeof reloadElements === "function") {
          await reloadElements(null).catch(() => {});
        }
        if (typeof pollElements === "function") pollElements();
      }
    } catch (e) {
      setStatus(e.message, true);
    } finally {
      btn.disabled = false;
    }
  }

  // ---------- 对话 ----------
  function providerPayload() {
    const s = state.settings;
    return {
      base_url: s.base_url.trim(),
      api_key: s.api_key.trim(),
      model: s.model.trim(),
      temperature: s.temperature,
    };
  }

  function providerReady() {
    const s = state.settings;
    return s.base_url.trim() !== "" && s.model.trim() !== "";
  }

  async function send() {
    const input = $("aiInput");
    const text = input.value.trim();
    if (!text || state.busy) return;
    if (!providerReady()) {
      setStatus("请先填写提供方地址与模型（可点「提供方设置」选预设）", true);
      $("aiSettings").open = true;
      $("aiBase").focus();
      return;
    }
    input.value = "";
    state.messages.push({ role: "user", content: text });
    state.messages.push({ role: "assistant", content: "", pending: true });
    renderLog();
    state.busy = true;
    setWire(true);
    setStatus("正在询问 " + state.settings.model + " …");
    state.abort = new AbortController();

    const body = {
      provider: providerPayload(),
      messages: state.messages.filter((m) => !m.pending),
      state: collectState(),
      stream: true,
    };
    try {
      const r = await fetch("/api/ai/chat", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(body),
        signal: state.abort.signal,
      });
      if (!r.ok || !r.body) {
        const t = await r.text();
        let d = {};
        try { d = JSON.parse(t); } catch (e) { d = { error: t }; }
        throw new Error(d.error || "HTTP " + r.status);
      }
      const reader = r.body.getReader();
      const dec = new TextDecoder();
      let buf = "";
      const last = () => state.messages[state.messages.length - 1];
      for (;;) {
        const { value, done } = await reader.read();
        if (done) break;
        buf += dec.decode(value, { stream: true });
        let idx;
        while ((idx = buf.indexOf("\n\n")) >= 0) {
          const frame = buf.slice(0, idx);
          buf = buf.slice(idx + 2);
          const line = frame.split("\n").find((l) => l.startsWith("data:"));
          if (!line) continue;
          let ev;
          try { ev = JSON.parse(line.slice(5).trim()); } catch (e) { continue; }
          if (ev.type === "delta") {
            last().content += ev.text;
            renderLog();
          } else if (ev.type === "error") {
            last().content += (last().content ? "\n\n" : "") + "（出错了）" + ev.error;
            renderLog();
            setStatus(ev.error, true);
          } else if (ev.type === "done") {
            if (ev.content && !last().content) last().content = ev.content;
          }
        }
      }
      const fin = last();
      fin.pending = false;
      renderLog();
      setStatus(fin.content ? "回答完成" : "没有收到内容");
    } catch (e) {
      const fin = state.messages[state.messages.length - 1];
      if (fin) fin.pending = false;
      if (e && e.name === "AbortError") {
        setStatus("已停止");
      } else {
        setStatus(e.message, true);
        if (fin && !fin.content) fin.content = "请求失败：" + e.message;
      }
      renderLog();
    } finally {
      state.busy = false;
      state.abort = null;
      setWire(false);
      $("aiInput").focus();
    }
  }

  function stop() {
    if (state.abort) state.abort.abort();
  }

  function setWire(busy) {
    $("aiSend").disabled = busy;
    $("aiStop").disabled = !busy;
    $("aiClear").disabled = busy;
    $("aiInput").disabled = busy;
  }

  // ---------- 设置面板 ----------
  function fillPresets() {
    const sel = $("aiPreset");
    PRESETS.forEach((p) => {
      const o = document.createElement("option");
      o.value = p.id;
      o.textContent = p.label;
      sel.appendChild(o);
    });
  }

  function applyPreset(id) {
    const p = PRESETS.find((x) => x.id === id);
    if (!p) return;
    state.settings.preset = p.id;
    if (p.base) $("aiBase").value = p.base;
    if (p.model) $("aiModel").value = p.model;
    if (id !== "custom") {
      state.settings.base_url = $("aiBase").value;
      state.settings.model = $("aiModel").value;
    }
    saveSettings();
  }

  function readForm() {
    state.settings.base_url = $("aiBase").value;
    state.settings.model = $("aiModel").value;
    state.settings.api_key = $("aiKey").value;
    const t = Number($("aiTemp").value);
    state.settings.temperature = t > 0 ? t : 0.3;
    saveSettings();
  }

  async function testConnection() {
    readForm();
    if (!providerReady()) { setStatus("请先填写提供方地址与模型", true); return; }
    const btn = $("aiTest");
    btn.disabled = true;
    setStatus("正在测试连接…");
    try {
      const d = await api("/api/ai/chat", {
        provider: providerPayload(),
        messages: [{ role: "user", content: "用一句话确认你能看到下面的规范摘要，并说出其中列出的第一个元件类型。" }],
        state: { mode: "test" },
        stream: false,
      });
      setStatus("连接正常：" + String(d.content || "").slice(0, 120));
    } catch (e) {
      setStatus("连接失败：" + e.message, true);
    } finally {
      btn.disabled = false;
    }
  }

  // ---------- 规范摘要（回答「AI 到底读了什么」） ----------
  async function loadSpec() {
    try {
      const r = await fetch("/api/ai/context");
      const d = await r.json();
      state.specChars = d.chars || 0;
      state.specSections = d.sections || [];
      const el = $("aiSpec");
      if (el) {
        el.textContent = "规范摘要 " + state.specChars.toLocaleString() + " 字符 · " +
          state.specSections.length + " 节（来自本二进制的元件目录，点此查看）";
      }
    } catch (e) { /* 摘要加载失败不影响对话 */ }
  }

  function showSpec() {
    const box = $("aiSpecText");
    if (!box) return;
    if (box.hidden) {
      fetch("/api/ai/context").then((r) => r.json()).then((d) => {
        box.textContent = d.spec || "";
        box.hidden = false;
      }).catch(() => { box.textContent = "规范摘要加载失败"; box.hidden = false; });
    } else {
      box.hidden = true;
    }
  }

  // ---------- 开关 ----------
  function setOpen(open) {
    state.open = open;
    $("aiPanel").hidden = !open;
    const btn = $("aiBtn");
    btn.classList.toggle("active", open);
    btn.setAttribute("aria-pressed", open ? "true" : "false");
    if (open) {
      loadSpec();
      renderLog();
      $("aiInput").focus();
    }
  }

  function wire() {
    $("aiBtn").addEventListener("click", () => setOpen(!state.open));
    $("aiClose").addEventListener("click", () => setOpen(false));
    $("aiSend").addEventListener("click", send);
    $("aiStop").addEventListener("click", stop);
    $("aiClear").addEventListener("click", () => {
      if (state.busy) return;
      state.messages = [];
      renderLog();
      setStatus("对话已清空");
    });
    $("aiTest").addEventListener("click", testConnection);
    $("aiPreset").addEventListener("change", (e) => applyPreset(e.target.value));
    ["aiBase", "aiModel", "aiKey", "aiTemp"].forEach((id) => {
      $(id).addEventListener("change", readForm);
    });
    const specLink = $("aiSpec");
    if (specLink) specLink.addEventListener("click", showSpec);
    $("aiInput").addEventListener("keydown", (e) => {
      // Enter 发送，Shift+Enter 换行；输入框内的按键不触发 app.js 的全局快捷键。
      if (e.key === "Enter" && !e.shiftKey) { e.preventDefault(); send(); }
      e.stopPropagation();
    });
    document.querySelectorAll(".aiQuick button").forEach((b) => {
      b.addEventListener("click", () => {
        const q = b.dataset.q;
        if (!q) return;
        $("aiInput").value = q;
        $("aiInput").focus();
      });
    });
  }

  function init() {
    fillPresets();
    const s = state.settings;
    $("aiBase").value = s.base_url;
    $("aiModel").value = s.model;
    $("aiKey").value = s.api_key;
    $("aiTemp").value = String(s.temperature);
    $("aiPreset").value = PRESETS.some((p) => p.id === s.preset) ? s.preset : "custom";
    wire();
    renderLog();
    loadSpec();
    // 侧栏默认隐藏；?ai=1 可直接打开（例如把带助手的地址存成书签）。
    if (location.hash === "#ai" || new URLSearchParams(location.search).get("ai") === "1") setOpen(true);
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init);
  } else {
    init();
  }
})();
