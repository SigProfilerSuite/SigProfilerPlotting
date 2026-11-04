"""Security-hardened local HTTP server for the SigProfilerPlotting dashboard.

The server is deliberately minimal (stdlib ``http.server`` only, no new runtime
dependency) and deliberately strict about what it will talk to: it binds only
to loopback addresses, rejects requests with an unrecognized ``Host`` header,
and rejects state-changing requests with an unrecognized ``Origin`` header.
The browser's native file picker is the only way to supply a matrix -- the
server never accepts a client-supplied filesystem path, which removes the
path-traversal attack surface for the one piece of user input that matters
most.
"""

from __future__ import annotations

import json
import mimetypes
import re
import threading
import traceback
from email import policy
from email.parser import BytesParser
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlsplit

from .matrix import MatrixError, load_matrix
from .renderers import render_sample

STATIC_ROOT = (Path(__file__).parent / "static").resolve()

_EXTRA_MIME_TYPES = {".js": "text/javascript; charset=utf-8", ".css": "text/css; charset=utf-8"}

_NAME_RE = re.compile(r'name="([^"]*)"')

# Defense-in-depth: 'unsafe-inline' is only for style-src (Plotly sets inline
# styles on its generated elements); the dashboard's own JS has no inline
# scripts, so script-src is left at the 'self' default.
_CONTENT_SECURITY_POLICY = (
    "default-src 'self'; style-src 'self' 'unsafe-inline'; "
    "img-src 'self' data: blob:; connect-src 'self'"
)


class DashboardHTTPServer(ThreadingHTTPServer):
    """A :class:`ThreadingHTTPServer` that also holds the uploaded matrix.

    Only one matrix is held at a time, matching the spec ("browse for or
    upload *one* matrix"); uploading a new one simply replaces it.
    """

    daemon_threads = True

    def __init__(self, server_address, handler_class):
        super().__init__(server_address, handler_class)
        self.dashboard_matrix = None


def make_server(host: str = "127.0.0.1", port: int = 0) -> DashboardHTTPServer:
    """Build (but do not start) the dashboard server.

    Callers start it with ``srv.serve_forever()`` (typically in a daemon
    thread) and stop it with ``srv.shutdown()`` followed by
    ``srv.server_close()`` -- both are the unmodified ``ThreadingHTTPServer``
    implementations, so that exact sequence works cleanly. The bound port is
    always readable via ``srv.server_address[1]``, including when ``port=0``
    asks the OS for an ephemeral port.
    """

    if host not in ("127.0.0.1", "localhost", "::1"):
        raise ValueError("The dashboard server only binds to loopback addresses.")
    return DashboardHTTPServer((host, port), DashboardRequestHandler)


def _allowed_host(header_value: str | None, port: int) -> bool:
    if not header_value:
        return False
    host = header_value.split(",", 1)[0].strip()
    return host in (f"127.0.0.1:{port}", f"localhost:{port}")


def _allowed_origin(header_value: str | None, port: int) -> bool:
    if not header_value:
        return False
    origin = header_value.strip().rstrip("/")
    return origin in (f"http://127.0.0.1:{port}", f"http://localhost:{port}")


def _extract_matrix_payload(raw: bytes, content_type: str) -> bytes:
    """Pull the uploaded matrix bytes out of a (possibly multipart) POST body."""
    if not content_type.startswith("multipart/form-data"):
        if not raw:
            raise MatrixError("No matrix data was received.")
        return raw
    # Hand the body to the stdlib email parser by prefixing the real
    # Content-Type header -- this reuses a well-tested MIME multipart parser
    # instead of a hand-rolled boundary splitter (the `cgi` module used by
    # older examples was removed in Python 3.13).
    header = f"Content-Type: {content_type}\r\nMIME-Version: 1.0\r\n\r\n".encode("utf-8")
    message = BytesParser(policy=policy.compat32).parsebytes(header + raw)
    if not message.is_multipart():
        raise MatrixError("Malformed multipart upload.")
    for part in message.get_payload():
        disposition = part.get("Content-Disposition", "")
        match = _NAME_RE.search(disposition)
        field_name = match.group(1) if match else None
        if field_name == "matrix" or "filename=" in disposition:
            payload = part.get_payload(decode=True)
            if not payload:
                raise MatrixError("The uploaded matrix field is empty.")
            return payload
    raise MatrixError("The upload did not include a 'matrix' file field.")


class DashboardRequestHandler(BaseHTTPRequestHandler):
    server_version = "SigProfilerPlottingDashboard/1.0"

    def log_message(self, format, *args):  # noqa: A002 - matches base signature
        pass  # the launcher prints its own user-facing status lines

    def _port(self) -> int:
        return self.server.server_address[1]

    def _reject(self, status: HTTPStatus, message: str) -> None:
        self._write_json(status, {"error": message})

    def _write_json(self, status: HTTPStatus, payload: dict) -> None:
        body = json.dumps(payload).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Content-Security-Policy", _CONTENT_SECURITY_POLICY)
        self.end_headers()
        self.wfile.write(body)

    def _write_bytes(self, status: HTTPStatus, body: bytes, content_type: str) -> None:
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Content-Security-Policy", _CONTENT_SECURITY_POLICY)
        self.end_headers()
        self.wfile.write(body)

    def _check_host(self) -> bool:
        if not _allowed_host(self.headers.get("Host"), self._port()):
            self._reject(HTTPStatus.BAD_REQUEST, "Unrecognized Host header.")
            return False
        return True

    # -- GET -----------------------------------------------------------
    def do_GET(self) -> None:
        if not self._check_host():
            return
        parsed = urlsplit(self.path)
        path = parsed.path
        try:
            if path == "/":
                self._serve_index()
            elif path == "/plotly.js":
                self._serve_plotlyjs()
            elif path.startswith("/static/"):
                self._serve_static(path[len("/static/"):])
            elif path == "/api/render":
                self._serve_render(parse_qs(parsed.query))
            else:
                self._reject(HTTPStatus.NOT_FOUND, "Not found.")
        except BrokenPipeError:
            pass
        except Exception as exc:  # noqa: BLE001 - see _handle_unexpected_error
            self._handle_unexpected_error(exc)

    def _serve_index(self) -> None:
        index_path = STATIC_ROOT / "index.html"
        if not index_path.is_file():
            self._reject(HTTPStatus.NOT_FOUND, "The dashboard UI has not been built yet.")
            return
        self._write_bytes(HTTPStatus.OK, index_path.read_bytes(), "text/html; charset=utf-8")

    def _serve_plotlyjs(self) -> None:
        import plotly.offline as pyo

        self._write_bytes(
            HTTPStatus.OK, pyo.get_plotlyjs().encode("utf-8"), "text/javascript; charset=utf-8"
        )

    def _serve_static(self, relative: str) -> None:
        # Resolve and confine to STATIC_ROOT before touching the filesystem,
        # so a crafted "/static/../../..." path can never escape the
        # package's static directory (the dashboard's only disk-read route).
        candidate = (STATIC_ROOT / relative).resolve()
        if candidate != STATIC_ROOT and STATIC_ROOT not in candidate.parents:
            self._reject(HTTPStatus.NOT_FOUND, "Not found.")
            return
        if not candidate.is_file():
            self._reject(HTTPStatus.NOT_FOUND, "Not found.")
            return
        content_type = _EXTRA_MIME_TYPES.get(candidate.suffix) or (
            mimetypes.guess_type(str(candidate))[0] or "application/octet-stream"
        )
        self._write_bytes(HTTPStatus.OK, candidate.read_bytes(), content_type)

    def _serve_render(self, query: dict) -> None:
        matrix = self.server.dashboard_matrix
        if matrix is None:
            self._reject(HTTPStatus.CONFLICT, "Upload a matrix before requesting a render.")
            return
        sample = query.get("sample", [None])[0]
        unit = query.get("unit", ["counts"])[0]
        if not sample:
            self._reject(HTTPStatus.BAD_REQUEST, "The 'sample' query parameter is required.")
            return
        try:
            figure = render_sample(matrix, sample, unit)
        except MatrixError as exc:
            self._reject(HTTPStatus.BAD_REQUEST, str(exc))
            return
        self._write_json(HTTPStatus.OK, figure)

    # -- POST ------------------------------------------------------------
    def do_POST(self) -> None:
        if not self._check_host():
            return
        try:
            self._dispatch_post()
        except BrokenPipeError:
            pass
        except Exception as exc:  # noqa: BLE001 - see _handle_unexpected_error
            self._handle_unexpected_error(exc)

    def _dispatch_post(self) -> None:
        parsed = urlsplit(self.path)
        if parsed.path == "/api/shutdown":
            if not _allowed_origin(self.headers.get("Origin"), self._port()):
                self._reject(HTTPStatus.FORBIDDEN, "Unrecognized Origin header.")
                return
            self._serve_shutdown()
            return
        if parsed.path != "/api/matrix":
            self._reject(HTTPStatus.NOT_FOUND, "Not found.")
            return
        if not _allowed_origin(self.headers.get("Origin"), self._port()):
            self._reject(HTTPStatus.FORBIDDEN, "Unrecognized Origin header.")
            return
        length = int(self.headers.get("Content-Length", 0) or 0)
        raw = self.rfile.read(length) if length else b""
        content_type = self.headers.get("Content-Type", "")
        try:
            payload = _extract_matrix_payload(raw, content_type)
            matrix = load_matrix(payload)
        except MatrixError as exc:
            self._reject(HTTPStatus.BAD_REQUEST, str(exc))
            return
        self.server.dashboard_matrix = matrix
        self._write_json(
            HTTPStatus.OK,
            {
                "context": matrix.context,
                "samples": list(matrix.columns),
                "labels_count": len(matrix.labels),
            },
        )

    def _serve_shutdown(self) -> None:
        # ThreadingHTTPServer runs each request on its own worker thread, so
        # it's safe to signal shutdown from here -- but server.shutdown()
        # blocks the calling thread until the serve_forever() loop (running
        # on a *different* thread, typically the launcher's main thread)
        # notices and exits. Do that handoff on a background thread so this
        # handler can send its response and return immediately rather than
        # stalling the client for the shutdown poll interval.
        self._write_json(HTTPStatus.OK, {"status": "stopping"})
        threading.Thread(target=self.server.shutdown, daemon=True).start()

    def _handle_unexpected_error(self, exc: Exception) -> None:
        # A bug anywhere in parsing/rendering must not silently drop the
        # connection (the client only ever sees an opaque "Failed to fetch"
        # for that) -- log the real traceback server-side for diagnosis and
        # still return a normal JSON error response.
        traceback.print_exc()
        try:
            self._reject(HTTPStatus.INTERNAL_SERVER_ERROR, f"Internal server error: {exc}")
        except Exception:
            pass  # the connection is already broken; nothing more to do
