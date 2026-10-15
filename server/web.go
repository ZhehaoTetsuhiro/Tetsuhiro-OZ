package server

import (
	"embed"
	"io/fs"
	"mime"
	"net/http"
)

// webFS is the single-page GUI: index.html, app.js (wave/quantum), oz.js
// (geometry/coating tab), style.css, favicon.svg and the vendored font. It is
// compiled into the binary, so `oz serve` needs no files next to it.
//
//go:embed web
var webFS embed.FS

// The embedded GUI ships a vendored web font; register its MIME type so the
// file server answers with font/woff2 instead of sniffing an octet-stream.
func init() {
	_ = mime.AddExtensionType(".woff2", "font/woff2")
	_ = mime.AddExtensionType(".woff", "font/woff")
}

// WebHandler serves the embedded GUI at the path it is mounted on (normally
// "/"). It answers GET/HEAD only for files that exist in the bundle. Responses
// must be revalidated on every load: index.html, app.js and oz.js change
// together with the binary, and a browser that keeps an old index.html beside
// a new script misbehaves. no-cache still lets the browser reuse a copy when
// the server answers 304 Not Modified.
func WebHandler() http.Handler {
	sub, err := fs.Sub(webFS, "web")
	if err != nil {
		// The embed pattern is fixed at compile time; this cannot fail.
		panic("server: embedded web bundle: " + err.Error())
	}
	files := http.FileServer(http.FS(sub))
	return http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		w.Header().Set("Cache-Control", "no-cache")
		files.ServeHTTP(w, r)
	})
}
