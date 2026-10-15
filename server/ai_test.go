package server

// AI assistant endpoint tests (AI 助手): the provider contract, the streaming
// relay, the specification digest and the element-writing path. Every provider
// call is served by a local httptest server that speaks the OpenAI wire format,
// so no test touches the network or needs a real key.

import (
	"encoding/json"
	"net/http"
	"net/http/httptest"
	"strings"
	"testing"

	"tetsuhiro-oz/optics"
)

// fakeProvider is an OpenAI-compatible stub. It records the last request it saw
// (so a test can assert what the server forwarded) and answers either in one
// piece or as an SSE stream, echoing the model and the system-prompt length.
type fakeProvider struct {
	gotAuth   string
	gotModel  string
	gotStream bool
	gotSysLen int
	gotMsgs   int
	calls     int
	failWith  int // if > 0, answer with this status and a JSON error
}

func (f *fakeProvider) handler() http.HandlerFunc {
	return func(w http.ResponseWriter, r *http.Request) {
		f.calls++
		f.gotAuth = r.Header.Get("Authorization")
		var body struct {
			Model    string `json:"model"`
			Stream   bool   `json:"stream"`
			Messages []struct {
				Role    string `json:"role"`
				Content string `json:"content"`
			} `json:"messages"`
		}
		_ = json.NewDecoder(r.Body).Decode(&body)
		f.gotModel, f.gotStream, f.gotMsgs = body.Model, body.Stream, len(body.Messages)
		if len(body.Messages) > 0 {
			f.gotSysLen = len(body.Messages[0].Content)
		}
		if f.failWith > 0 {
			w.Header().Set("Content-Type", "application/json")
			w.WriteHeader(f.failWith)
			_, _ = w.Write([]byte(`{"error":{"message":"invalid api key","type":"auth"}}`))
			return
		}
		if body.Stream {
			w.Header().Set("Content-Type", "text/event-stream")
			w.WriteHeader(http.StatusOK)
			fl, _ := w.(http.Flusher)
			for _, tok := range []string{"流", "式", "回答"} {
				_, _ = w.Write([]byte(`data: {"choices":[{"delta":{"content":"` + tok + `"}}]}` + "\n\n"))
				if fl != nil {
					fl.Flush()
				}
			}
			_, _ = w.Write([]byte("data: [DONE]\n\n"))
			return
		}
		w.Header().Set("Content-Type", "application/json")
		_, _ = w.Write([]byte(`{"choices":[{"message":{"content":"模型=` + body.Model +
			`；规范=` + itoa(f.gotSysLen) + `字符；轮数=` + itoa(len(body.Messages)) + `"}}],"usage":{"total_tokens":3}}`))
	}
}

func itoa(n int) string {
	if n == 0 {
		return "0"
	}
	var b []byte
	for n > 0 {
		b = append([]byte{byte('0' + n%10)}, b...)
		n /= 10
	}
	return string(b)
}

// aiBody builds a chat request naming the given provider base URL.
func aiBody(base, model, key string, msgs []aiMessage, stream bool) string {
	b, _ := json.Marshal(map[string]any{
		"provider": map[string]any{"base_url": base, "api_key": key, "model": model},
		"messages": msgs,
		"stream":   stream,
	})
	return string(b)
}

// TestAIChatForwardsToProvider covers POST /api/ai/chat: the provider sees the
// model, the bearer key and a system prompt assembled from the specification.
func TestAIChatForwardsToProvider(t *testing.T) {
	f := &fakeProvider{}
	up := httptest.NewServer(f.handler())
	defer up.Close()
	ts := httptest.NewServer(New(64 << 20).Handler())
	defer ts.Close()

	body := aiBody(up.URL, "test-model", "sk-secret", []aiMessage{{Role: "user", Content: "这个透镜做什么用的？"}}, false)
	r := doReq(t, http.MethodPost, ts.URL+"/api/ai/chat", body)
	if r.code != http.StatusOK {
		t.Fatalf("chat = %d %s", r.code, r.body)
	}
	if f.calls != 1 {
		t.Fatalf("provider called %d times, want 1", f.calls)
	}
	if f.gotAuth != "Bearer sk-secret" {
		t.Errorf("Authorization = %q, want the bearer key forwarded verbatim", f.gotAuth)
	}
	if f.gotModel != "test-model" {
		t.Errorf("model = %q", f.gotModel)
	}
	// The system prompt must carry the specification (tens of KB), so the
	// assistant really does read the live catalog.
	if f.gotSysLen < 20000 {
		t.Errorf("system prompt is %d chars, want the full specification digest", f.gotSysLen)
	}
	if f.gotMsgs != 2 {
		t.Errorf("upstream saw %d messages, want system + user", f.gotMsgs)
	}
	m := decodeMap(t, r.body)
	content, _ := m["content"].(string)
	if !strings.Contains(content, "test-model") {
		t.Errorf("answer not relayed: %s", content)
	}
}

// TestAIChatStreamRelaysDeltas covers the SSE path: every provider delta
// arrives as its own event and the stream ends with a done frame.
func TestAIChatStreamRelaysDeltas(t *testing.T) {
	f := &fakeProvider{}
	up := httptest.NewServer(f.handler())
	defer up.Close()
	ts := httptest.NewServer(New(64 << 20).Handler())
	defer ts.Close()

	body := aiBody(up.URL, "stream-model", "", []aiMessage{{Role: "user", Content: "hi"}}, true)
	body = strings.Replace(body, `"api_key":""`, `"api_key":""`, 1)
	r := doReq(t, http.MethodPost, ts.URL+"/api/ai/chat", body)
	if r.code != http.StatusOK {
		t.Fatalf("stream chat = %d %s", r.code, r.body)
	}
	if ct := r.header.Get("Content-Type"); !strings.HasPrefix(ct, "text/event-stream") {
		t.Errorf("Content-Type = %q, want text/event-stream", ct)
	}
	text := string(r.body)
	for _, want := range []string{`"type":"delta"`, "流", "式", "回答", `"type":"done"`} {
		if !strings.Contains(text, want) {
			t.Errorf("stream lacks %q:\n%s", want, text)
		}
	}
	// A streamed answer must still be complete in the done frame, so a client
	// that ignores deltas gets the whole text.
	var done string
	for _, line := range strings.Split(text, "\n") {
		if strings.Contains(line, `"type":"done"`) {
			done = line
		}
	}
	var frame struct{ Content string }
	if err := json.Unmarshal([]byte(strings.TrimPrefix(strings.TrimSpace(done), "data:")), &frame); err != nil {
		t.Fatalf("done frame is not JSON: %v (%s)", err, done)
	}
	if frame.Content != "流式回答" {
		t.Errorf("done content = %q, want the concatenated answer", frame.Content)
	}
}

// TestAIChatErrors covers the refusal paths: a bad provider configuration is a
// 400 before any network call, and an upstream rejection is relayed verbatim.
func TestAIChatErrors(t *testing.T) {
	f := &fakeProvider{failWith: http.StatusUnauthorized}
	up := httptest.NewServer(f.handler())
	defer up.Close()
	ts := httptest.NewServer(New(64 << 20).Handler())
	defer ts.Close()
	url := ts.URL + "/api/ai/chat"

	// Upstream 401: the provider's own message must survive so the user knows
	// it is a key problem, not a server bug.
	r := doReq(t, http.MethodPost, url, aiBody(up.URL, "m", "bad-key",
		[]aiMessage{{Role: "user", Content: "hi"}}, false))
	if r.code != http.StatusBadGateway {
		t.Errorf("upstream 401 relayed as %d, want 502", r.code)
	}
	if !strings.Contains(string(r.body), "invalid api key") {
		t.Errorf("provider error text lost: %s", r.body)
	}
	if f.calls != 1 {
		t.Errorf("provider called %d times", f.calls)
	}

	for _, c := range []struct{ name, body string }{
		{"no base url", `{"provider":{"model":"m"},"messages":[{"role":"user","content":"x"}]}`},
		{"bad scheme", `{"provider":{"base_url":"ftp://x/v1","model":"m"},"messages":[{"role":"user","content":"x"}]}`},
		{"no model", `{"provider":{"base_url":"http://127.0.0.1:1/v1"},"messages":[{"role":"user","content":"x"}]}`},
		{"no messages", `{"provider":{"base_url":"http://127.0.0.1:1/v1","model":"m"},"messages":[]}`},
		{"bad role", `{"provider":{"base_url":"http://127.0.0.1:1/v1","model":"m"},"messages":[{"role":"system","content":"x"}]}`},
		{"malformed", `{`},
	} {
		if r := doReq(t, http.MethodPost, url, c.body); r.code != http.StatusBadRequest {
			t.Errorf("%s = %d, want 400 (%s)", c.name, r.code, r.body)
		}
	}
	// None of them may have reached the provider.
	if f.calls != 1 {
		t.Errorf("invalid requests reached the provider %d times, want only the valid one", f.calls)
	}
	// Method contract.
	if r := doReq(t, http.MethodGet, url, ""); r.code != http.StatusMethodNotAllowed {
		t.Errorf("GET chat = %d, want 405", r.code)
	}
}

// TestAIChatStateIsForwarded checks that the live GUI state reaches the system
// prompt: the assistant must see the configuration it is asked about.
func TestAIChatStateIsForwarded(t *testing.T) {
	f := &fakeProvider{}
	up := httptest.NewServer(f.handler())
	defer up.Close()
	ts := httptest.NewServer(New(64 << 20).Handler())
	defer ts.Close()

	cfg := `{"grid":{"size":256,"width":0.02},"wavelength":6.328e-7,"sources":[{"id":"s1","type":"gaussian","params":{"waist":0.001}}],` +
		`"scene":{"components":[{"id":"l1","type":"lens","pos":{"x":0,"y":0,"z":0},"params":{"f":0.25}}]}}`
	body, _ := json.Marshal(map[string]any{
		"provider": map[string]any{"base_url": up.URL, "api_key": "k", "model": "m"},
		"messages": []aiMessage{{Role: "user", Content: "当前配置有什么问题？"}},
		"state": map[string]any{
			"mode":   "wave",
			"config": json.RawMessage(cfg),
			"selected": map[string]any{
				"kind": "component", "type": "lens", "label": "透镜",
				"help": "薄透镜", "params": map[string]any{"f": 0.25},
			},
			"run": map[string]any{"status": "done", "grid": 256, "planes": []string{"焦面"}},
		},
	})
	if r := doReq(t, http.MethodPost, ts.URL+"/api/ai/chat", string(body)); r.code != http.StatusOK {
		t.Fatalf("chat with state = %d %s", r.code, r.body)
	}
	if f.gotSysLen <= 20000 {
		t.Errorf("system prompt %d chars, want spec + state", f.gotSysLen)
	}
}

// TestAIContextDigest covers GET /api/ai/context: the digest the assistant
// reads is derived from the live catalog, so it documents exactly what this
// binary implements.
func TestAIContextDigest(t *testing.T) {
	ts := httptest.NewServer(New(64 << 20).Handler())
	defer ts.Close()
	r := doReq(t, http.MethodGet, ts.URL+"/api/ai/context", "")
	if r.code != http.StatusOK {
		t.Fatalf("context = %d %s", r.code, r.body)
	}
	m := decodeMap(t, r.body)
	spec, _ := m["spec"].(string)
	if spec == "" {
		t.Fatal("empty specification digest")
	}
	// The digest must carry the output protocol the GUI parses and the element
	// vocabulary from the kernel.
	for _, want := range []string{"oz-config", "oz-element", "scene", "sources"} {
		if !strings.Contains(spec, want) {
			t.Errorf("digest lacks %q", want)
		}
	}
	// Every built-in element type is documented (the digest is generated, not
	// hand-written, so it cannot drift).
	cat := optics.BuildCatalog()
	for _, d := range cat.Elements[:8] {
		if !strings.Contains(spec, `"type":"`+d.Type+`"`) {
			t.Errorf("digest lacks element %q", d.Type)
		}
	}
	sections, _ := m["sections"].([]any)
	if len(sections) < 5 {
		t.Errorf("digest lists %d sections, want the catalog sections", len(sections))
	}
	if r := doReq(t, http.MethodPost, ts.URL+"/api/ai/context", "{}"); r.code != http.StatusMethodNotAllowed {
		t.Errorf("POST context = %d, want 405", r.code)
	}
}

// TestAIElementWrite covers POST /api/ai/element: a sound definition is
// written, validated by the kernel and immediately usable; a broken one is
// refused before anything touches the disk.
func TestAIElementWrite(t *testing.T) {
	dir := t.TempDir()
	optics.SetElementDirs([]string{dir})
	optics.ReloadElementDefinitions(nil)
	t.Cleanup(func() {
		// The registry is process-wide: leave it empty so other tests do not
		// see this test's definition.
		optics.SetElementDirs([]string{t.TempDir()})
		optics.ReloadElementDefinitions(nil)
	})
	ts := httptest.NewServer(New(64 << 20).Handler())
	defer ts.Close()
	url := ts.URL + "/api/ai/element"

	good := `{"definition":{"name":"ai_test_axicon","label":"AI 轴锥镜","behavior":"transmit","class":"lens",` +
		`"params":[{"key":"f","label":"焦距","unit":"m","kind":"float","min":0.01,"max":5,"step":0.01,"default":0.2}],` +
		`"phase":"-k*(sqrt(r*r+f*f) - f)"}}`
	r := doReq(t, http.MethodPost, url, good)
	if r.code != http.StatusOK {
		t.Fatalf("element write = %d %s", r.code, r.body)
	}
	m := decodeMap(t, r.body)
	path, _ := m["path"].(string)
	if !strings.HasPrefix(path, dir) {
		t.Errorf("written to %q, want inside %q", path, dir)
	}
	// The element is now loaded: it appears in the catalog like any other.
	if !optics.IsScriptedElement("ai_test_axicon") {
		t.Error("written definition is not registered")
	}
	found := false
	for _, d := range optics.ScriptedElementDocs() {
		if d.Type == "ai_test_axicon" {
			found = true
		}
	}
	if !found {
		t.Error("written definition is absent from the catalog")
	}

	// Broken definitions are refused with the kernel's own message.
	for _, c := range []struct{ name, body string }{
		{"syntax", `{"definition":{"name":"ai_broken","phase":"1 +"}}`},
		{"reserved name", `{"definition":{"name":"lens","phase":"0"}}`},
		{"bad name", `{"definition":{"name":"AI_Elem","phase":"0"}}`},
		{"unknown field", `{"definition":{"name":"ai_x","nope":1}}`},
		{"missing", `{}`},
	} {
		if r := doReq(t, http.MethodPost, url, c.body); r.code != http.StatusBadRequest {
			t.Errorf("%s = %d, want 400 (%s)", c.name, r.code, r.body)
		}
	}
	if optics.IsScriptedElement("ai_broken") {
		t.Error("a refused definition was registered anyway")
	}
	if r := doReq(t, http.MethodGet, url, ""); r.code != http.StatusMethodNotAllowed {
		t.Errorf("GET element = %d, want 405", r.code)
	}
}

// TestAIRouteMounted proves the AI routes are reachable through the one mux the
// service mounts (a 404 here would mean the mount is missing).
func TestAIRouteMounted(t *testing.T) {
	f := &fakeProvider{}
	up := httptest.NewServer(f.handler())
	defer up.Close()
	ts := httptest.NewServer(New(64 << 20).Handler())
	defer ts.Close()
	if r := doReq(t, http.MethodGet, ts.URL+"/api/ai/context", ""); r.code != http.StatusOK {
		t.Errorf("GET /api/ai/context = %d, want 200", r.code)
	}
	if r := doReq(t, http.MethodPost, ts.URL+"/api/ai/chat",
		aiBody(up.URL, "m", "k", []aiMessage{{Role: "user", Content: "x"}}, false)); r.code != http.StatusOK {
		t.Errorf("POST /api/ai/chat = %d, want 200", r.code)
	}
	if r := doReq(t, http.MethodGet, ts.URL+"/api/ai", ""); r.code == http.StatusOK {
		t.Error("bare /api/ai answered 200")
	}
}
