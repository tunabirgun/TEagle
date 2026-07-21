"""TEagle local backend — stdlib HTTP server exposing the real scientific core.
Serves the web UI and JSON APIs. No external web framework (reproducibility + fewer deps).
All request validation and scientific behaviour live in engine.py (shared with the native app);
this module is only the HTTP transport."""
from __future__ import annotations
import json, os, sys
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import engine
from engine import BadRequest
from teagle_core import appdirs

WEB = appdirs.resource("web") or os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "web"))
MIME = {".html": "text/html", ".css": "text/css", ".js": "application/javascript",
        ".svg": "image/svg+xml", ".ico": "image/x-icon"}


class Handler(BaseHTTPRequestHandler):
    timeout = 30                                # fail fast on a lying/short Content-Length instead of pinning a worker

    def log_message(self, *a):  # quiet
        pass

    def _send(self, code, body, ctype="application/json"):
        data = body if isinstance(body, (bytes, bytearray)) else json.dumps(body).encode()
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def _body(self):
        n = int(self.headers.get("Content-Length", 0) or 0)
        if n < 0:
            raise BadRequest("invalid Content-Length")
        try:
            data = json.loads(self.rfile.read(n) or b"{}")
        except (json.JSONDecodeError, ValueError, RecursionError):   # deeply nested JSON raises RecursionError -> 400, not 500
            raise BadRequest("malformed JSON request body")
        if not isinstance(data, dict):
            raise BadRequest("request body must be a JSON object")
        return data

    def do_GET(self):
        path = self.path.split("?")[0]
        if path == "/api/health":
            return self._send(200, engine.run_health())
        if path == "/api/env":
            return self._send(200, engine.run_env())
        if path == "/api/wsl/status":                  # full annotation-stack status (slower)
            return self._send(200, engine.run_wsl_status())
        if path == "/api/wsl/install_log":
            return self._send(200, engine.run_wsl_install_log())
        rel = "index.html" if path in ("/", "") else path.lstrip("/")
        fp = os.path.abspath(os.path.join(WEB, rel))
        if not fp.startswith(WEB) or not os.path.isfile(fp):   # path-traversal guard
            return self._send(404, {"error": "not found"})
        ext = os.path.splitext(fp)[1]
        try:
            with open(fp, "rb") as f:
                data = f.read()
        except OSError as e:                                 # TOCTOU race / permission / disconnect -> 500, not a bare traceback
            return self._send(500, {"error": f"could not read {rel}: {type(e).__name__}"})
        self._send(200, data, MIME.get(ext, "application/octet-stream"))

    # POST path -> engine handler. Each returns a dict or raises BadRequest (-> 400).
    _ROUTES = {
        "/api/fetch": engine.run_fetch,
        "/api/analyze": engine.run_analyze,
        "/api/wsl/install": lambda b: engine.run_wsl_install(),
        "/api/eta": engine.run_eta,
        "/api/annotate": engine.run_annotate,
        "/api/splice": engine.run_splice,
        "/api/miniprot": engine.run_miniprot,
        "/api/primers": engine.run_primers,
        "/api/pcr": engine.run_pcr,
    }

    def do_POST(self):
        try:
            n = int(self.headers.get("Content-Length", 0) or 0)
            if n < 0:                                            # a negative length must not bypass the cap
                return self._send(400, {"error": "invalid Content-Length"})
            if n > 80_000_000:                                   # request-size cap
                return self._send(413, {"error": "request too large (> 80 MB)"})
            handler = self._ROUTES.get(self.path)
            if handler is None:
                return self._send(404, {"error": "unknown endpoint"})
            self._send(200, handler(self._body()))
        except BadRequest as e:
            self._send(400, {"error": str(e)})               # malformed request -> 400, not 500
        except Exception as e:
            import traceback
            traceback.print_exc()                            # full detail to the server log
            self._send(500, {"error": f"{type(e).__name__}: {e}"})


def make_server(host="127.0.0.1", port=8765):
    """Bind the engine socket (raises OSError on a busy/blocked port) so the caller can bind in the
    main thread and surface the failure, instead of losing it in a daemon thread."""
    return ThreadingHTTPServer((host, port), Handler)


def serve(host="127.0.0.1", port=8765):
    srv = make_server(host, port)
    print(f"TEagle backend on http://{host}:{port}  (primer3 {primers.PRIMER3_VERSION})")
    srv.serve_forever()


if __name__ == "__main__":
    serve(port=int(sys.argv[1]) if len(sys.argv) > 1 else 8765)
