package server

// AI assistant endpoints (AI 助手):
//
//	POST /api/ai/chat        chat with the user's own provider (SSE when stream:true)
//	GET  /api/ai/context     the specification digest the assistant is given
//	POST /api/ai/element     validate and write a scripted element definition file
//
// The provider (base URL, API key, model) is supplied by the browser with every
// request and is never stored, logged or written to disk by this server: the key
// stays in the user's browser (localStorage) and is forwarded to the provider
// over TLS exactly like a direct call would be. Keeping the call server-side is
// what lets the assistant read the live specification and the current state
// without the browser having to know either, and it sidesteps the CORS rules
// most OpenAI-compatible providers enforce.
//
// The specification digest (see aiSpec) is derived from optics.BuildCatalog(),
// so the assistant documents exactly the elements, sources, methods and shapes
// this binary implements — no hand-kept copy that could drift.

import (
	"bufio"
	"bytes"
	"context"
	"encoding/json"
	"fmt"
	"io"
	"log"
	"net/http"
	"net/url"
	"os"
	"path/filepath"
	"strings"
	"sync"
	"time"

	"tetsuhiro-oz/optics"
)

// Limits that keep the assistant endpoints bounded.
const (
	// aiMaxBody caps the chat request (messages + state). The state may carry
	// the whole configuration, so this is generous but finite.
	aiMaxBody = 4 << 20
	// aiMaxStateChars caps the rendered "current state" section of the system
	// prompt, so one huge configuration cannot crowd out the specification.
	aiMaxStateChars = 24000
	// aiMaxSpecChars caps the specification digest itself.
	aiMaxSpecChars = 48000
	// aiMaxMessages bounds the conversation length a client may replay.
	aiMaxMessages = 60
	// aiMaxTokensDefault / aiMaxTokensLimit bound the completion length.
	aiMaxTokensDefault = 2048
	aiMaxTokensLimit   = 8192
	// aiTimeout is the hard deadline of one provider call; the client-side
	// context cancels it earlier when the browser goes away.
	aiTimeout = 15 * time.Minute
	// aiErrBodyLimit caps how much of a provider error body is relayed.
	aiErrBodyLimit = 2000
	// aiProviderURLLimit / aiKeyLimit / aiModelLimit cap the credential fields.
	aiProviderURLLimit = 512
	aiKeyLimit         = 512
	aiModelLimit       = 128
)

// aiProvider is the user's own model endpoint. Every field comes from the
// browser; nothing here is persisted.
type aiProvider struct {
	BaseURL     string  `json:"base_url"`
	APIKey      string  `json:"api_key"`
	Model       string  `json:"model"`
	Temperature float64 `json:"temperature"`
	MaxTokens   int     `json:"max_tokens"`
}

// aiMessage is one conversation turn.
type aiMessage struct {
	Role    string `json:"role"`
	Content string `json:"content"`
}

// aiSelected describes the element or source the user is looking at, so the
// assistant can answer "what does this do?" about the right thing.
type aiSelected struct {
	Kind   string         `json:"kind"` // "component" | "source" | "element"
	Type   string         `json:"type"`
	Label  string         `json:"label,omitempty"`
	Help   string         `json:"help,omitempty"`
	Params map[string]any `json:"params,omitempty"`
}

// aiRunSummary is the digest of the last simulation the GUI shows.
type aiRunSummary struct {
	Status    string   `json:"status,omitempty"`
	ElapsedMS float64  `json:"elapsed_ms,omitempty"`
	Grid      int      `json:"grid,omitempty"`
	Planes    []string `json:"planes,omitempty"`
	Warnings  []string `json:"warnings,omitempty"`
}

// aiState is the live GUI state the assistant should read before answering.
// Raw blobs are relayed verbatim by the server; only their length is capped.
type aiState struct {
	Mode     string          `json:"mode,omitempty"`
	Config   json.RawMessage `json:"config,omitempty"`
	Selected *aiSelected     `json:"selected,omitempty"`
	Warnings []string        `json:"warnings,omitempty"`
	Run      *aiRunSummary   `json:"run,omitempty"`
	// Quantum and geometry/coating modes keep their working document in the
	// page; the GUI passes it here as text so the assistant sees the same
	// numbers the user does.
	Quantum string `json:"quantum,omitempty"`
	Lens    string `json:"lens,omitempty"`
	Coating string `json:"coating,omitempty"`
}

// aiChatRequest is the POST /api/ai/chat body.
type aiChatRequest struct {
	Provider aiProvider  `json:"provider"`
	Messages []aiMessage `json:"messages"`
	State    *aiState    `json:"state"`
	Stream   bool        `json:"stream"`
}

// aiChatResponse is the non-streaming answer.
type aiChatResponse struct {
	Content string `json:"content"`
	Model   string `json:"model,omitempty"`
	Usage   any    `json:"usage,omitempty"`
}

// handleAI dispatches every /api/ai/ path.
func (s *Server) handleAI(w http.ResponseWriter, r *http.Request) {
	rest := strings.Trim(strings.TrimPrefix(r.URL.Path, "/api/ai"), "/")
	switch rest {
	case "":
		writeErr(w, http.StatusNotFound, "use /api/ai/chat, /api/ai/context or /api/ai/element")
		return
	case "chat":
		s.handleAIChat(w, r)
	case "context":
		s.handleAIContext(w, r)
	case "element":
		s.handleAIElement(w, r)
	default:
		writeErr(w, http.StatusNotFound, "unknown ai sub-resource "+rest)
	}
}

// handleAIContext answers GET /api/ai/context with the specification digest the
// assistant receives. It is the honest answer to "what does the AI read?": the
// GUI shows its size and can display it in full.
func (s *Server) handleAIContext(w http.ResponseWriter, r *http.Request) {
	if r.Method != http.MethodGet {
		writeErr(w, http.StatusMethodNotAllowed, "GET only")
		return
	}
	spec := aiSpec()
	writeJSON(w, http.StatusOK, map[string]any{
		"spec":        spec,
		"chars":       len(spec),
		"sections":    aiSpecSections(spec),
		"generated":   time.Now().Format(time.RFC3339),
		"spec_source": "optics.BuildCatalog()",
	})
}

// handleAIChat forwards one conversation turn to the user's provider.
func (s *Server) handleAIChat(w http.ResponseWriter, r *http.Request) {
	if r.Method != http.MethodPost {
		writeErr(w, http.StatusMethodNotAllowed, "POST only")
		return
	}
	body, err := io.ReadAll(http.MaxBytesReader(w, r.Body, aiMaxBody))
	if err != nil {
		writeErr(w, http.StatusBadRequest, "failed to read request body: "+err.Error())
		return
	}
	var req aiChatRequest
	if err := json.Unmarshal(body, &req); err != nil {
		writeErr(w, http.StatusBadRequest, "invalid JSON: "+err.Error())
		return
	}
	if err := req.Provider.validate(); err != nil {
		writeErr(w, http.StatusBadRequest, err.Error())
		return
	}
	if len(req.Messages) == 0 {
		writeErr(w, http.StatusBadRequest, "messages must not be empty")
		return
	}
	if len(req.Messages) > aiMaxMessages {
		writeErr(w, http.StatusBadRequest, fmt.Sprintf("too many messages (%d, limit %d)", len(req.Messages), aiMaxMessages))
		return
	}
	for i, m := range req.Messages {
		switch m.Role {
		case "user", "assistant":
		default:
			writeErr(w, http.StatusBadRequest, fmt.Sprintf("messages[%d].role must be user or assistant", i))
			return
		}
	}

	upstream := aiUpstreamMessages(aiSystemPrompt(req.State), req.Messages)
	if req.Stream {
		s.streamAIChat(w, r, req.Provider, upstream)
		return
	}
	var full strings.Builder
	content, usage, err := callAIProvider(r.Context(), req.Provider, upstream, false, func(d string) { full.WriteString(d) })
	if err != nil {
		writeErr(w, aiErrorStatus(err), err.Error())
		return
	}
	writeJSON(w, http.StatusOK, aiChatResponse{Content: content, Model: req.Provider.Model, Usage: usage})
}

// aiErrorStatus maps a provider failure to an HTTP status: a rejected
// credential or model is the caller's mistake (400), anything else is an
// upstream problem (502).
func aiErrorStatus(err error) int {
	msg := err.Error()
	if strings.Contains(msg, "provider answered") {
		return http.StatusBadGateway
	}
	return http.StatusBadRequest
}

// streamAIChat relays the provider's token stream to the browser as SSE. Each
// event is a JSON object: {"type":"delta","text":…}, {"type":"done"} or
// {"type":"error","error":…}.
func (s *Server) streamAIChat(w http.ResponseWriter, r *http.Request, p aiProvider, msgs []aiMessage) {
	flusher, ok := w.(http.Flusher)
	if !ok {
		// No streaming possible (HTTP/1.0 test client): fall back to one shot.
		content, _, err := callAIProvider(r.Context(), p, msgs, false, func(string) {})
		if err != nil {
			writeErr(w, aiErrorStatus(err), err.Error())
			return
		}
		writeJSON(w, http.StatusOK, aiChatResponse{Content: content, Model: p.Model})
		return
	}
	w.Header().Set("Content-Type", "text/event-stream; charset=utf-8")
	w.Header().Set("Cache-Control", "no-cache")
	w.Header().Set("Connection", "keep-alive")
	w.Header().Set("Access-Control-Allow-Origin", "*")
	w.WriteHeader(http.StatusOK)
	flusher.Flush()

	send := func(v any) {
		buf, err := json.Marshal(v)
		if err != nil {
			return
		}
		// A payload may not contain a raw newline: it would end the event.
		buf = bytes.ReplaceAll(buf, []byte("\n"), []byte(" "))
		_, _ = w.Write(append([]byte("data: "), buf...))
		_, _ = w.Write([]byte("\n\n"))
		flusher.Flush()
	}
	// A browser that navigates away or clicks 停止 cancels r.Context(); the
	// provider call dies with it and the loop below simply stops relaying.
	content, _, err := callAIProvider(r.Context(), p, msgs, true, func(delta string) {
		if delta == "" {
			return
		}
		send(map[string]any{"type": "delta", "text": delta})
	})
	if err != nil {
		// Whatever already streamed stays visible; the error is appended so the
		// user sees both the partial answer and why it stopped.
		send(map[string]any{"type": "error", "error": err.Error()})
		return
	}
	send(map[string]any{"type": "done", "content": content})
}

// callAIProvider performs one chat completion against the user's provider. It
// is used by both the streaming and the one-shot path: onDelta receives every
// streamed chunk as it arrives (and, in the one-shot path, the whole answer),
// and the returned string is the concatenation of the answer.
func callAIProvider(ctx context.Context, p aiProvider, msgs []aiMessage, stream bool, onDelta func(string)) (string, any, error) {
	payload := map[string]any{
		"model":       p.Model,
		"messages":    msgs,
		"temperature": p.temperature(),
		"max_tokens":  p.maxTokens(),
		"stream":      stream,
	}
	buf, err := json.Marshal(payload)
	if err != nil {
		return "", nil, fmt.Errorf("failed to encode provider request: %v", err)
	}
	ctx, cancel := context.WithTimeout(ctx, aiTimeout)
	defer cancel()
	req, err := http.NewRequestWithContext(ctx, http.MethodPost, p.endpoint(), bytes.NewReader(buf))
	if err != nil {
		return "", nil, fmt.Errorf("invalid provider URL: %v", err)
	}
	req.Header.Set("Content-Type", "application/json")
	req.Header.Set("Accept", "application/json")
	if p.APIKey != "" {
		req.Header.Set("Authorization", "Bearer "+p.APIKey)
	}
	resp, err := http.DefaultClient.Do(req)
	if err != nil {
		return "", nil, fmt.Errorf("provider request failed: %v", err)
	}
	defer resp.Body.Close()
	if resp.StatusCode < 200 || resp.StatusCode > 299 {
		detail, _ := io.ReadAll(io.LimitReader(resp.Body, aiErrBodyLimit))
		msg := strings.TrimSpace(string(detail))
		if msg == "" {
			msg = resp.Status
		}
		return "", nil, fmt.Errorf("provider answered HTTP %d: %s", resp.StatusCode, msg)
	}
	if !stream {
		var out struct {
			Choices []struct {
				Message struct {
					Content string `json:"content"`
				} `json:"message"`
			} `json:"choices"`
			Usage any `json:"usage"`
		}
		if err := json.NewDecoder(io.LimitReader(resp.Body, aiMaxBody)).Decode(&out); err != nil {
			return "", nil, fmt.Errorf("provider answer is not valid JSON: %v", err)
		}
		if len(out.Choices) == 0 {
			return "", nil, fmt.Errorf("provider answer has no choices")
		}
		onDelta(out.Choices[0].Message.Content)
		return out.Choices[0].Message.Content, out.Usage, nil
	}
	// Streaming: relay each delta as it arrives and accumulate the full text.
	var text strings.Builder
	sc := bufio.NewScanner(io.LimitReader(resp.Body, aiMaxBody*16))
	sc.Buffer(make([]byte, 0, 64*1024), 4<<20)
	for sc.Scan() {
		line := strings.TrimSpace(sc.Text())
		if !strings.HasPrefix(line, "data:") {
			continue
		}
		data := strings.TrimSpace(strings.TrimPrefix(line, "data:"))
		if data == "" {
			continue
		}
		if data == "[DONE]" {
			break
		}
		var chunk struct {
			Choices []struct {
				Delta struct {
					Content string `json:"content"`
				} `json:"delta"`
			} `json:"choices"`
		}
		if err := json.Unmarshal([]byte(data), &chunk); err != nil {
			continue // keep-alive comments and provider-specific frames
		}
		for _, c := range chunk.Choices {
			if c.Delta.Content == "" {
				continue
			}
			text.WriteString(c.Delta.Content)
			onDelta(c.Delta.Content)
		}
	}
	if err := sc.Err(); err != nil {
		return text.String(), nil, fmt.Errorf("provider stream ended early: %v", err)
	}
	return text.String(), nil, nil
}

// validate rejects a provider configuration that cannot work, before any
// network call is made.
func (p *aiProvider) validate() error {
	u, err := url.Parse(strings.TrimSpace(p.BaseURL))
	if err != nil || u.Host == "" {
		return fmt.Errorf("provider.base_url must be an absolute http(s) URL")
	}
	switch u.Scheme {
	case "http", "https":
	default:
		return fmt.Errorf("provider.base_url must use http or https (got %q)", u.Scheme)
	}
	if len(p.BaseURL) > aiProviderURLLimit {
		return fmt.Errorf("provider.base_url is too long (limit %d)", aiProviderURLLimit)
	}
	if len(p.APIKey) > aiKeyLimit {
		return fmt.Errorf("provider.api_key is too long (limit %d)", aiKeyLimit)
	}
	if strings.TrimSpace(p.Model) == "" {
		return fmt.Errorf("provider.model is required")
	}
	if len(p.Model) > aiModelLimit {
		return fmt.Errorf("provider.model is too long (limit %d)", aiModelLimit)
	}
	return nil
}

// endpoint returns the chat-completions URL of the provider: the base URL with
// any trailing slash removed, plus /chat/completions unless already present.
func (p *aiProvider) endpoint() string {
	base := strings.TrimRight(strings.TrimSpace(p.BaseURL), "/")
	if strings.HasSuffix(base, "/chat/completions") {
		return base
	}
	return base + "/chat/completions"
}

func (p *aiProvider) temperature() float64 {
	if p.Temperature <= 0 {
		return 0.3
	}
	if p.Temperature > 2 {
		return 2
	}
	return p.Temperature
}

func (p *aiProvider) maxTokens() int {
	if p.MaxTokens <= 0 {
		return aiMaxTokensDefault
	}
	if p.MaxTokens > aiMaxTokensLimit {
		return aiMaxTokensLimit
	}
	return p.MaxTokens
}

// aiUpstreamMessages prepends the system prompt to the user's conversation.
func aiUpstreamMessages(system string, msgs []aiMessage) []aiMessage {
	out := make([]aiMessage, 0, len(msgs)+1)
	out = append(out, aiMessage{Role: "system", Content: system})
	return append(out, msgs...)
}

// ---- scripted element writing -------------------------------------------

// aiElementRequest is the POST /api/ai/element body: one element definition.
type aiElementRequest struct {
	Definition json.RawMessage `json:"definition"`
}

// handleAIElement validates an element definition and, when it is sound,
// writes it into the first scripted-element directory and reloads the
// registry. The definition is compiled by the kernel itself
// (optics.ParseElementDefinition), so a broken expression, a reserved name or
// an unknown field is refused before anything touches the disk.
func (s *Server) handleAIElement(w http.ResponseWriter, r *http.Request) {
	if r.Method != http.MethodPost {
		writeErr(w, http.StatusMethodNotAllowed, "POST only")
		return
	}
	body, err := io.ReadAll(http.MaxBytesReader(w, r.Body, 1<<20))
	if err != nil {
		writeErr(w, http.StatusBadRequest, "failed to read request body: "+err.Error())
		return
	}
	var req aiElementRequest
	if err := json.Unmarshal(body, &req); err != nil {
		writeErr(w, http.StatusBadRequest, "invalid JSON: "+err.Error())
		return
	}
	raw := bytes.TrimSpace(req.Definition)
	if len(raw) == 0 {
		writeErr(w, http.StatusBadRequest, "definition is required")
		return
	}
	// Validate first: the kernel reports the exact reason a definition is bad.
	if _, err := optics.ParseElementDefinition(raw, ""); err != nil {
		writeErr(w, http.StatusBadRequest, err.Error())
		return
	}
	var def optics.ElementDefinition
	if err := json.Unmarshal(raw, &def); err != nil {
		writeErr(w, http.StatusBadRequest, "invalid JSON: "+err.Error())
		return
	}
	dir, err := aiElementDir()
	if err != nil {
		writeErr(w, http.StatusInternalServerError, err.Error())
		return
	}
	path := filepath.Join(dir, def.Name+".json")
	pretty, err := json.MarshalIndent(&def, "", "  ")
	if err != nil {
		writeErr(w, http.StatusInternalServerError, "failed to encode definition: "+err.Error())
		return
	}
	pretty = append(pretty, '\n')
	if err := os.WriteFile(path, pretty, 0o644); err != nil {
		writeErr(w, http.StatusInternalServerError, "failed to write "+path+": "+err.Error())
		return
	}
	report := optics.ReloadElementDefinitions(nil)
	log.Printf("元件定义已由 AI 助手写入 %s（重载后 %d 个脚本元件）", path, len(report.Loaded))
	writeJSON(w, http.StatusOK, map[string]any{
		"ok": true, "name": def.Name, "path": path, "report": report,
	})
}

// aiElementDir returns the directory new definitions are written to: the first
// configured directory that exists, or the first one created on demand.
func aiElementDir() (string, error) {
	dirs := optics.ElementDirs()
	if len(dirs) == 0 {
		dirs = optics.DefaultElementDirs()
	}
	if len(dirs) == 0 {
		return "", fmt.Errorf("no element definition directory is configured")
	}
	for _, d := range dirs {
		if st, err := os.Stat(d); err == nil && st.IsDir() {
			return d, nil
		}
	}
	if err := os.MkdirAll(dirs[0], 0o755); err != nil {
		return "", fmt.Errorf("failed to create %s: %v", dirs[0], err)
	}
	return dirs[0], nil
}

// ---- specification digest ------------------------------------------------

var (
	aiSpecOnce sync.Once
	aiSpecText string
)

// aiSpec returns the specification the assistant reads, built once from the
// live catalog so it can never drift from what this binary implements.
func aiSpec() string {
	aiSpecOnce.Do(func() { aiSpecText = buildAISpec() })
	return aiSpecText
}

// aiSpecHead is the fixed part of the system prompt: conventions, the output
// protocol the GUI understands, and the element-definition file format.
const aiSpecHead = `你是 Tetsuhiro-OZ 的内置光学助手。OZ 是一台波动/几何/量子光学模拟器：一个 Go 内核（零第三方依赖）+ 一个内嵌网页 GUI，GUI 通过 HTTP API 与内核通信。你回答用户关于光学与当前光路的问题，并在需要时给出可直接应用的配置或元件定义。

【单位与约定】
- 长度一律是米（m），功率是瓦（W），角度是弧度（rad），波长是米（632.8e-9 表示 632.8 nm）。
- 波动配置 config 的顶层字段：grid{size,width}、wavelength、polarized、method、evanescent("decay"|"zero")、evanescent_limit、backward_regularize、tikhonov_alpha、bandlimit{fraction,sigma}、sources[]、scene{components[]}。
- 场景 scene 由带位置的元件组成：pos{x,y,z}（z 沿光轴）、yaw/pitch/roll、shape{kind,params}（通光孔形状）、params（元件参数）。光路按元件位置自动路由，与列表顺序无关；场景至少要有一个 sensor 才能记录光。
- 元件类型、参数名与取值范围以下面的元件目录为准，不要编造。脚本元件（custom:true）来自 elements/*.json 定义文件，可随写随用。

【你能做的事】
1. 解释元件、参数、物理含义、单位与典型量级（例如“这个视场角下色差有多大”）。
2. 依据当前配置找出问题并给出具体改法（写出字段路径与新值）。
3. 设计光路：输出一个完整、可直接运行的配置。
4. 定义新的脚本元件：输出元件定义 JSON，服务端会用内核校验后写入元件目录，立即可在 GUI 的插入列表中使用。

【输出协议（GUI 会解析，必须遵守）】
- 完整配置放在 ` + "```oz-config" + ` 代码块里，内容是完整 config JSON（含 grid/wavelength/sources/scene），一个字段都不要省略；GUI 会提供「应用到场景」按钮。
- 元件定义放在 ` + "```oz-element" + ` 代码块里；GUI 会提供「写入元件定义」按钮（服务端校验通过才落盘）。
- 其余解释用普通文字。不要编造不存在的元件类型、参数名或端点。
- 用中文回答，简洁、可执行；给出数值时带上单位。

【脚本元件定义文件格式（elements/*.json）】
- name：小写标识符（字母开头，可含数字与下划线，≤48 字符），不能与内置元件同名。
- label / help：中文标签与说明；behavior："transmit"（默认）或 "mirror"；class：lens|mirror|splitter|detector|stop|other。
- params：[{key,label,unit,kind(float|int|choice|bool|text),min,max,step,default,choices[],help,show_if}]；key 是标识符，避免使用 x,y,r,th,wl,k,pi,e 与函数名，也不能用保留名 kind,size,width,wl,mask,scale,cmap,pmin,pmax,gamma,part。
- phase / amp：在元件局部坐标 (x,y)（米）上的表达式，如 "-k*(sqrt(r*r+f*f) - f)"、"(1 + m*cos(2*pi*x/Lambda))/(1+m)"。可用变量 x,y,r=sqrt(x²+y²),th=atan2(y,x),wl（波长 m）,k=2π/λ,pi,e，以及本元件各参数名。
- 可用函数：sin cos tan asin acos atan atan2 sinh cosh tanh exp log log10 sqrt abs sign floor ceil round min max pow hypot mod clamp if(条件,真,假) step rect smoothstep erf sinc(=sin(πx)/(πx)) rad deg。
- phase_unit："rad"（默认）或 "waves"；jones 可选 2×2 琼斯矩阵 {a:{re,im},b:{re,im},c:{re,im},d:{re,im}}，各项为表达式。
- 示例：{"name":"metalens","label":"超表面透镜","behavior":"transmit","class":"lens","params":[{"key":"f","label":"焦距","unit":"m","kind":"float","min":0.01,"max":10,"step":0.01,"default":0.3}],"phase":"-k*(sqrt(r*r+f*f) - f)"}

【HTTP 端点（GUI 用的就是这些）】
GET /api/catalog（元件/光源/方法/示例文档）、GET /api/health、GET /api/elements、POST /api/elements/reload、POST /api/validate（配置校验）、POST /api/simulate（提交仿真）、GET /api/runs/{id}[/planes|profiles|inspect|scene]、POST /api/convert（旧元件序列转场景）、POST /api/quantum、POST /api/lens/{design,eval,tol,thermal,sens}、POST /api/geo/{mtf,ghost}、POST /api/doe/focus、POST /api/coating/{design,spectrum}。
`

// aiSpecTail introduces the live state section.
const aiSpecTail = "\n【当前状态】\n"

// buildAISpec renders the specification digest from the live catalog.
func buildAISpec() string {
	cat := optics.BuildCatalog()
	var b strings.Builder
	b.WriteString(aiSpecHead)

	writeDocs := func(title string, docs []optics.ElementDoc) {
		if len(docs) == 0 {
			return
		}
		b.WriteString("\n【" + title + "】\n")
		for _, d := range docs {
			entry := map[string]any{"type": d.Type, "label": d.Label}
			if h := aiTrim(d.Help, 160); h != "" {
				entry["help"] = h
			}
			if d.Custom {
				entry["custom"] = true
			}
			if d.Class != "" {
				entry["class"] = d.Class
			}
			params := make([]map[string]any, 0, len(d.Params))
			for _, p := range d.Params {
				if p.Kind == "nested" {
					continue
				}
				e := map[string]any{"key": p.Key, "label": p.Label}
				if p.Unit != "" {
					e["unit"] = p.Unit
				}
				if p.Kind != "" {
					e["kind"] = p.Kind
				}
				if p.Default != nil {
					e["default"] = p.Default
				}
				if len(p.Choices) > 0 {
					e["choices"] = p.Choices
				}
				if p.Min != 0 || p.Max != 0 {
					e["range"] = []float64{p.Min, p.Max}
				}
				if h := aiTrim(p.Help, 120); h != "" {
					e["help"] = h
				}
				if p.ShowIf != "" {
					e["show_if"] = p.ShowIf
				}
				params = append(params, e)
			}
			if len(params) > 0 {
				entry["params"] = params
			}
			line, err := json.Marshal(entry)
			if err != nil {
				continue
			}
			b.Write(line)
			b.WriteByte('\n')
		}
	}
	writeDocs("元件目录", cat.Elements)
	writeDocs("光源目录", cat.Sources)
	if len(cat.Shapes) > 0 {
		b.WriteString("\n【通光孔形状 scene.components[].shape】\n")
		for _, sh := range cat.Shapes {
			keys := make([]string, 0, len(sh.Params))
			for _, p := range sh.Params {
				keys = append(keys, p.Key)
			}
			line, _ := json.Marshal(map[string]any{"kind": sh.Kind, "label": sh.Label, "params": keys})
			b.Write(line)
			b.WriteByte('\n')
		}
	}
	if len(cat.Methods) > 0 {
		b.WriteString("\n【传播方法 config.method】\n")
		for _, m := range cat.Methods {
			mm, ok := m.(map[string]string)
			if !ok {
				continue
			}
			line, _ := json.Marshal(map[string]string{"key": mm["key"], "label": mm["label"], "help": aiTrim(mm["help"], 140)})
			b.Write(line)
			b.WriteByte('\n')
		}
	}
	if len(cat.Examples) > 0 {
		names := make([]string, 0, len(cat.Examples))
		for _, ex := range cat.Examples {
			names = append(names, ex.Name)
		}
		b.WriteString("\n【内置示例（GUI 左侧「模板」下拉，完整配置见 /api/catalog）】\n")
		b.WriteString(strings.Join(names, "、") + "\n")
	}
	if len(cat.Quantum.States) > 0 || len(cat.Quantum.Gates) > 0 {
		b.WriteString("\n【量子模式 /api/quantum：初态与门（modes/cutoff/state/gates）】\n")
		for _, d := range append(append([]optics.ElementDoc(nil), cat.Quantum.States...), cat.Quantum.Gates...) {
			keys := make([]string, 0, len(d.Params))
			for _, p := range d.Params {
				keys = append(keys, p.Key)
			}
			line, _ := json.Marshal(map[string]any{"type": d.Type, "label": d.Label, "params": keys, "help": aiTrim(d.Help, 120)})
			b.Write(line)
			b.WriteByte('\n')
		}
	}
	b.WriteString(aiSpecTail)
	out := b.String()
	if len(out) > aiMaxSpecChars {
		out = out[:aiMaxSpecChars] + "\n…（规范摘要超出长度上限，已截断）\n"
	}
	return out
}

// aiTrim shortens a help string so the digest stays compact.
func aiTrim(s string, n int) string {
	s = strings.Join(strings.Fields(s), " ")
	if len(s) <= n {
		return s
	}
	return s[:n] + "…"
}

// aiSpecSections lists the 【…】 headings of the digest, so the GUI can show
// what the assistant reads without dumping 40 KB of text.
func aiSpecSections(spec string) []string {
	var out []string
	for _, line := range strings.Split(spec, "\n") {
		if strings.HasPrefix(line, "【") && strings.HasSuffix(strings.TrimSpace(line), "】") {
			out = append(out, strings.Trim(line, "【】"))
		}
	}
	return out
}

// aiSystemPrompt assembles the system message: the specification digest plus
// the live GUI state.
func aiSystemPrompt(st *aiState) string {
	var b strings.Builder
	b.WriteString(aiSpec())
	if st == nil {
		b.WriteString("（GUI 未提供当前状态：仅依据规范回答，涉及当前光路时先问用户。）\n")
		return b.String()
	}
	if st.Mode != "" {
		b.WriteString("模式: " + st.Mode + "\n")
	}
	if len(st.Config) > 0 {
		b.WriteString("当前配置 config（JSON，单位 SI）: " + aiCap(string(st.Config), aiMaxStateChars) + "\n")
	}
	if st.Selected != nil {
		sel, err := json.Marshal(st.Selected)
		if err == nil {
			b.WriteString("用户正在查看的" + aiSelectedKind(st.Selected.Kind) + ": " + string(sel) + "\n")
		}
	}
	if len(st.Warnings) > 0 {
		b.WriteString("GUI 显示的警告: " + strings.Join(st.Warnings, " | ") + "\n")
	}
	if st.Run != nil {
		run, err := json.Marshal(st.Run)
		if err == nil {
			b.WriteString("最近一次仿真: " + string(run) + "\n")
		}
	}
	if st.Quantum != "" {
		b.WriteString("量子配置（文本）: " + aiCap(st.Quantum, 4000) + "\n")
	}
	if st.Lens != "" {
		b.WriteString("几何设计结果（文本）: " + aiCap(st.Lens, 4000) + "\n")
	}
	if st.Coating != "" {
		b.WriteString("膜系设计结果（文本）: " + aiCap(st.Coating, 4000) + "\n")
	}
	return b.String()
}

func aiSelectedKind(kind string) string {
	switch kind {
	case "source":
		return "光源"
	case "element":
		return "元件"
	default:
		return "元件"
	}
}

// aiCap truncates a blob to n bytes on a rune boundary.
func aiCap(s string, n int) string {
	if len(s) <= n {
		return s
	}
	cut := s[:n]
	for len(cut) > 0 && !utf8Cont(cut[len(cut)-1]) {
		cut = cut[:len(cut)-1]
	}
	return cut + "…（已截断）"
}

// utf8Cont reports whether b is a UTF-8 continuation byte, i.e. cutting here
// would split a character.
func utf8Cont(b byte) bool { return b&0xC0 == 0x80 }
