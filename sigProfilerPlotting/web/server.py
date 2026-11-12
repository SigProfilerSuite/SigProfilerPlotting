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

import io
import json
import mimetypes
import os
import secrets
import socket
import tempfile
import threading
import traceback
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlsplit

from . import examples
from .matrix import MatrixError, load_matrix
from .renderers import render_sample

STATIC_ROOT = (Path(__file__).parent / "static").resolve()

_EXTRA_MIME_TYPES = {
    ".js": "text/javascript; charset=utf-8",
    ".css": "text/css; charset=utf-8",
}

# A real SBS4608 matrix can run into the hundreds of MB of plain text
# (tens of thousands of sample columns); this is a safety cap against a
# truly unreasonable upload, not a tight budget -- the point of streaming
# straight to a temp file (_dispatch_post) is that memory use no longer
# scales with this number at all. Configurable since "unreasonable" is a
# judgment call that depends on the machine running the dashboard.
MAX_UPLOAD_BYTES = int(
    os.environ.get("SIGPROFILERPLOTTING_MAX_UPLOAD_BYTES", 4 * 1024 * 1024 * 1024)
)
_UPLOAD_CHUNK_SIZE = 1024 * 1024

# A rasterized chart PNG is at most a few MB even at a generous export scale;
# this is a safety cap against abuse, not a tight budget.
MAX_EXPORT_PNG_BYTES = 64 * 1024 * 1024
# Matches the static SBS4608 PDF export's own page-width convention (810pt,
# i.e. 11.25in) closely enough that a plot exported either way looks similar
# printed -- not meant to match it exactly, since this path serves a quick
# snapshot of whatever is on screen rather than a release-quality artifact.
_EXPORT_PDF_WIDTH_IN = 11.0

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
        # Initialize our state before HTTPServer binds.  socketserver calls
        # server_close() if binding fails, so the override below must be safe
        # even when super().__init__() never completes.
        self.dashboard_matrix = None
        # The temp file backing dashboard_matrix's lazy, disk-backed reads
        # (see matrix.py) -- owned by the server, not by MatrixData, so it
        # can be cleaned up here when replaced by a new upload or when the
        # server shuts down, rather than ever leaking.
        self._matrix_tempfile: str | None = None
        self._matrix_lock = threading.RLock()
        # A fresh random token per launch -- the launcher's pidfile records
        # it alongside the pid/port, and treats a pid/port whose live
        # /api/identity doesn't echo this exact token as unrelated to this
        # server (stale pidfile, reused pid, or a different process on the
        # same port), rather than ever trusting the pidfile's pid alone.
        self.instance_token = secrets.token_hex(16)
        super().__init__(server_address, handler_class)

    def set_dashboard_matrix(self, matrix, tempfile_path: str | None) -> None:
        """Replace the active matrix, cleaning up whatever temp file the
        previous one was reading from (if any) -- uploading a new matrix
        must not leak the old one's backing file."""
        with self._matrix_lock:
            previous = self._matrix_tempfile
            self.dashboard_matrix = matrix
            self._matrix_tempfile = tempfile_path
            if previous is not None:
                try:
                    os.unlink(previous)
                except OSError:
                    pass

    def server_close(self) -> None:
        with self._matrix_lock:
            if self._matrix_tempfile is not None:
                try:
                    os.unlink(self._matrix_tempfile)
                except OSError:
                    pass
                self._matrix_tempfile = None
            self.dashboard_matrix = None
        super().server_close()


class IPv6DashboardHTTPServer(DashboardHTTPServer):
    """IPv6 loopback variant used when ``make_server`` receives ``::1``."""

    address_family = socket.AF_INET6


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
    server_class = IPv6DashboardHTTPServer if host == "::1" else DashboardHTTPServer
    return server_class((host, port), DashboardRequestHandler)


def _allowed_host(header_value: str | None, port: int) -> bool:
    if not header_value:
        return False
    host = header_value.split(",", 1)[0].strip()
    return host in (f"127.0.0.1:{port}", f"localhost:{port}", f"[::1]:{port}")


def _allowed_origin(header_value: str | None, port: int) -> bool:
    if not header_value:
        return False
    origin = header_value.strip().rstrip("/")
    return origin in (
        f"http://127.0.0.1:{port}",
        f"http://localhost:{port}",
        f"http://[::1]:{port}",
    )


def _stream_upload_to_tempfile(rfile, content_length: int) -> str:
    """Write exactly `content_length` bytes from `rfile` to a new private
    temp file in bounded chunks, never holding the whole upload in memory
    at once (see the module docstring and matrix.py's). The caller owns
    the returned path and must clean it up on any failure path.
    """
    fd, path = tempfile.mkstemp(suffix=".sppmatrix")
    try:
        with os.fdopen(fd, "wb") as handle:
            remaining = content_length
            while remaining > 0:
                chunk = rfile.read(min(remaining, _UPLOAD_CHUNK_SIZE))
                if not chunk:
                    raise MatrixError(
                        "The upload ended early; fewer bytes were received than "
                        "the Content-Length header promised."
                    )
                handle.write(chunk)
                remaining -= len(chunk)
    except Exception:
        try:
            os.unlink(path)
        except OSError:
            pass
        raise
    return path


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
                self._serve_static(path[len("/static/") :])
            elif path == "/api/render":
                self._serve_render(parse_qs(parsed.query))
            elif path == "/api/cosmic/builds":
                self._serve_cosmic_builds()
            elif path == "/api/examples":
                self._serve_examples_list()
            elif path == "/api/identity":
                self._serve_identity()
            else:
                self._reject(HTTPStatus.NOT_FOUND, "Not found.")
        except BrokenPipeError:
            pass
        except Exception as exc:  # noqa: BLE001 - see _handle_unexpected_error
            self._handle_unexpected_error(exc)

    def _serve_index(self) -> None:
        index_path = STATIC_ROOT / "index.html"
        if not index_path.is_file():
            self._reject(
                HTTPStatus.NOT_FOUND, "The dashboard UI has not been built yet."
            )
            return
        self._write_bytes(
            HTTPStatus.OK, index_path.read_bytes(), "text/html; charset=utf-8"
        )

    def _serve_plotlyjs(self) -> None:
        import plotly.offline as pyo

        self._write_bytes(
            HTTPStatus.OK,
            pyo.get_plotlyjs().encode("utf-8"),
            "text/javascript; charset=utf-8",
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
        sample = query.get("sample", [None])[0]
        unit = query.get("unit", ["counts"])[0]
        compare_cosmic = query.get("compare_cosmic", [None])[0] or None
        if not sample:
            self._reject(
                HTTPStatus.BAD_REQUEST, "The 'sample' query parameter is required."
            )
            return
        try:
            # MatrixData reads its selected sample lazily from the uploaded
            # temp file.  Hold the same lock used by replacement/cleanup so
            # another upload cannot unlink that file mid-render.
            with self.server._matrix_lock:
                matrix = self.server.dashboard_matrix
                if matrix is None:
                    self._reject(
                        HTTPStatus.CONFLICT,
                        "Upload a matrix before requesting a render.",
                    )
                    return
                figure = render_sample(
                    matrix, sample, unit, compare_cosmic=compare_cosmic
                )
        except MatrixError as exc:
            self._reject(HTTPStatus.BAD_REQUEST, str(exc))
            return
        self._write_json(HTTPStatus.OK, figure)

    def _serve_cosmic_builds(self) -> None:
        from . import cosmic

        catalogues = [
            {
                "build": info.build,
                "version": info.version,
                "source_path": info.source_path,
                "checksum": info.checksum,
                "signature_count": info.signature_count,
            }
            for info in cosmic.catalogue_metadata()
        ]
        self._write_json(
            HTTPStatus.OK,
            {
                "builds": [c["build"] for c in catalogues],
                "catalogues": catalogues,
                "catalogue_dir": (
                    str(cosmic.catalogue_dir()) if cosmic.catalogue_dir() else None
                ),
            },
        )

    def _serve_examples_list(self) -> None:
        self._write_json(HTTPStatus.OK, {"examples": examples.available()})

    def _serve_identity(self) -> None:
        # Loopback-only (every GET already passes _check_host()) proof of
        # which process is actually listening here, and under which launch's
        # token -- see DashboardHTTPServer.instance_token.
        self._write_json(
            HTTPStatus.OK, {"pid": os.getpid(), "token": self.server.instance_token}
        )

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
        if parsed.path not in (
            "/api/shutdown",
            "/api/matrix",
            "/api/export/pdf",
            "/api/examples/load",
        ):
            self._reject(HTTPStatus.NOT_FOUND, "Not found.")
            return
        # Every POST route is state-changing (or at least does real work), so
        # all of them require a matching Origin, not just the upload/shutdown
        # routes this check originally guarded.
        if not _allowed_origin(self.headers.get("Origin"), self._port()):
            self._reject(HTTPStatus.FORBIDDEN, "Unrecognized Origin header.")
            return
        if parsed.path == "/api/shutdown":
            self._serve_shutdown()
            return
        if parsed.path == "/api/export/pdf":
            self._serve_export_pdf()
            return
        if parsed.path == "/api/examples/load":
            self._serve_example_load(parse_qs(parsed.query))
            return

        header = self.headers.get("Content-Length")
        try:
            length = int(header)
        except (TypeError, ValueError):
            self._reject(
                HTTPStatus.LENGTH_REQUIRED, "A valid Content-Length header is required."
            )
            return
        if length < 0:
            self._reject(HTTPStatus.BAD_REQUEST, "Content-Length must not be negative.")
            return
        if length == 0:
            self._reject(HTTPStatus.BAD_REQUEST, "No matrix data was received.")
            return
        if length > MAX_UPLOAD_BYTES:
            self._reject(
                HTTPStatus.REQUEST_ENTITY_TOO_LARGE,
                f"The upload ({length} bytes) exceeds the {MAX_UPLOAD_BYTES}-byte limit.",
            )
            return

        # Stream the body straight to a temp file in bounded chunks -- this
        # is the fix for the memory scaling plotting_A's review measured:
        # a real SBS4608 upload (tens of thousands of sample columns) must
        # never be held whole in memory, here or in matrix.py, which reads
        # lazily from this same file afterward.
        tmp_path = None
        try:
            tmp_path = _stream_upload_to_tempfile(self.rfile, length)
            matrix = load_matrix(tmp_path)
        except MatrixError as exc:
            if tmp_path is not None:
                try:
                    os.unlink(tmp_path)
                except OSError:
                    pass
            self._reject(HTTPStatus.BAD_REQUEST, str(exc))
            return
        self._respond_with_loaded_matrix(matrix, tmp_path)

    def _respond_with_loaded_matrix(self, matrix, tempfile_path: str | None) -> None:
        self.server.set_dashboard_matrix(matrix, tempfile_path)
        self._write_json(
            HTTPStatus.OK,
            {
                "context": matrix.context,
                "samples": list(matrix.columns),
                "labels_count": len(matrix.labels),
                "warnings": list(matrix.warnings),
            },
        )

    def _serve_example_load(self, query: dict) -> None:
        name = query.get("name", [None])[0]
        if not name:
            self._reject(
                HTTPStatus.BAD_REQUEST, "The 'name' query parameter is required."
            )
            return
        path = examples.example_path(name)
        if path is None:
            self._reject(
                HTTPStatus.NOT_FOUND,
                "Unknown example, or it is not bundled in this install.",
            )
            return
        try:
            # Read straight from the package's own bundled file -- no
            # temp-file copy needed, since load_matrix already reads lazily
            # from whatever path it's given (see matrix.py). Pass None as
            # the owned tempfile_path below: this path belongs to the
            # package, not to us, so it must never be unlinked the way an
            # upload's temp file is when replaced or on shutdown.
            matrix = load_matrix(path)
        except MatrixError as exc:
            self._reject(HTTPStatus.BAD_REQUEST, str(exc))
            return
        self._respond_with_loaded_matrix(matrix, None)

    def _serve_shutdown(self) -> None:
        # Defense-in-depth alongside the Origin check in _dispatch_post: the
        # launcher's own --stop only ever calls this once it has separately
        # verified (via /api/identity) that this exact process holds this
        # exact token, so a stale/foreign pidfile can never reach here with
        # a valid token.
        if self.headers.get("X-Dashboard-Token") != self.server.instance_token:
            self._reject(HTTPStatus.FORBIDDEN, "Invalid or missing dashboard token.")
            return
        # ThreadingHTTPServer runs each request on its own worker thread, so
        # it's safe to signal shutdown from here -- but server.shutdown()
        # blocks the calling thread until the serve_forever() loop (running
        # on a *different* thread, typically the launcher's main thread)
        # notices and exits. Do that handoff on a background thread so this
        # handler can send its response and return immediately rather than
        # stalling the client for the shutdown poll interval.
        self._write_json(HTTPStatus.OK, {"status": "stopping"})
        threading.Thread(target=self.server.shutdown, daemon=True).start()

    def _serve_export_pdf(self) -> None:
        """Wrap a browser-rendered plot PNG in a one-page PDF.

        The PNG itself is rasterized client-side by Plotly.js (the same
        mechanism the toolbar's own "download plot as png" button already
        uses) -- this endpoint only ever wraps an already-decoded raster
        image, so unlike the static SBS4608 PDF export it needs neither
        Kaleido nor a local Chrome/Chromium install.
        """
        header = self.headers.get("Content-Length")
        try:
            length = int(header)
        except (TypeError, ValueError):
            self._reject(
                HTTPStatus.LENGTH_REQUIRED, "A valid Content-Length header is required."
            )
            return
        if length <= 0:
            self._reject(HTTPStatus.BAD_REQUEST, "No image data was received.")
            return
        if length > MAX_EXPORT_PNG_BYTES:
            self._reject(
                HTTPStatus.REQUEST_ENTITY_TOO_LARGE,
                f"The image ({length} bytes) exceeds the {MAX_EXPORT_PNG_BYTES}-byte limit.",
            )
            return
        body = self.rfile.read(length)

        from PIL import Image

        try:
            with Image.open(io.BytesIO(body)) as image:
                image.load()
                if image.mode in ("RGBA", "LA") or (
                    image.mode == "P" and "transparency" in image.info
                ):
                    # Every chart already paints its own opaque white
                    # background, so this only matters defensively -- PIL's
                    # PDF writer requires RGB, and naively discarding alpha
                    # instead of flattening it onto white would blacken
                    # whatever was meant to show through as transparent.
                    image = image.convert("RGBA")
                    flattened = Image.new("RGB", image.size, "white")
                    flattened.paste(image, mask=image.getchannel("A"))
                    image = flattened
                elif image.mode != "RGB":
                    image = image.convert("RGB")
                # A fixed page width keeps the PDF a sensible physical size
                # regardless of how large Plotly happened to rasterize the
                # image; height follows from the image's own aspect ratio.
                dpi = image.width / _EXPORT_PDF_WIDTH_IN
                buffer = io.BytesIO()
                image.save(buffer, "PDF", resolution=dpi, quality=95, subsampling=0)
        except (OSError, ValueError) as exc:
            self._reject(
                HTTPStatus.BAD_REQUEST, f"Could not build a PDF from the image: {exc}"
            )
            return

        pdf_bytes = buffer.getvalue()
        self.send_response(HTTPStatus.OK)
        self.send_header("Content-Type", "application/pdf")
        self.send_header("Content-Length", str(len(pdf_bytes)))
        self.send_header("Content-Security-Policy", _CONTENT_SECURITY_POLICY)
        self.send_header("Content-Disposition", 'attachment; filename="plot.pdf"')
        self.end_headers()
        self.wfile.write(pdf_bytes)

    def _handle_unexpected_error(self, exc: Exception) -> None:
        # A bug anywhere in parsing/rendering must not silently drop the
        # connection (the client only ever sees an opaque "Failed to fetch"
        # for that) -- log the real traceback server-side for diagnosis and
        # still return a normal JSON error response.
        traceback.print_exc()
        try:
            self._reject(
                HTTPStatus.INTERNAL_SERVER_ERROR, f"Internal server error: {exc}"
            )
        except Exception:
            pass  # the connection is already broken; nothing more to do
