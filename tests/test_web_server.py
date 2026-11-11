import http.client
import json
import os
import socket
import threading

import pytest

from sigProfilerPlotting.web.server import make_server

SBS6_LABELS = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]


def _sbs96_tsv():
    from sigProfilerPlotting.web.matrix import REFERENCE_CONTEXTS

    labels = REFERENCE_CONTEXTS["SBS96"]
    lines = ["MutationType\tS1"]
    for index, label in enumerate(labels):
        lines.append(f"{label}\t{index + 1}")
    return "\n".join(lines) + "\n"


@pytest.fixture
def running_server():
    server = make_server("127.0.0.1", 0)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        yield server
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)


def _connection(server):
    port = server.server_address[1]
    return http.client.HTTPConnection("127.0.0.1", port, timeout=5)


def _matrix_body(content):
    # The dashboard's own frontend posts the raw file bytes as the whole
    # request body (app.js's upload()) -- not multipart/form-data -- so the
    # server can stream it straight to a temp file without ever parsing a
    # MIME envelope out of it first.
    return content.encode("utf-8"), "application/octet-stream"


def _sbs6_tsv():
    lines = ["MutationType\tS1\tS2"]
    for index, label in enumerate(SBS6_LABELS):
        lines.append(f"{label}\t{index}\t{index * 2}")
    return "\n".join(lines) + "\n"


def _id89_r_export_tsv():
    from sigProfilerPlotting.web.matrix import REFERENCE_CONTEXTS

    labels = list(REFERENCE_CONTEXTS["ID89"])
    labels[labels.index("Del(C):R(6,)")] = "Del(C):R(6,9)"
    lines = ["ID89Signature"]  # R row-name export: omitted top-left cell
    lines.extend(f"{label}\t{index + 1}" for index, label in enumerate(labels))
    return "\n".join(lines) + "\n"


def _id476_r_export_tsv():
    from sigProfilerPlotting.web.matrix import REFERENCE_CONTEXTS

    # The incomplete columns reproduce only the *shape* of the Liu R export;
    # all values and names here are synthetic.
    lines = ["Lower83\tComplete476\tEmpty"]
    for index, label in enumerate(REFERENCE_CONTEXTS["ID476"]):
        lower = str(index + 10) if index < 83 else "NA"
        lines.append(f"{label}\t{lower}\t{index + 1}\tNA")
    return "\n".join(lines) + "\n"


# SigProfilerPlotting does not bundle or redistribute COSMIC's signature
# catalogues -- configure a small, synthetic, user-supplied-style catalogue
# directory for the tests below instead of relying on real bundled data.
@pytest.fixture(autouse=True)
def _cosmic_catalogue(tmp_path):
    from sigProfilerPlotting.web import cosmic
    from sigProfilerPlotting.web.matrix import REFERENCE_CONTEXTS

    labels = list(REFERENCE_CONTEXTS["SBS96"])
    directory = tmp_path / "cosmic"
    directory.mkdir()
    names = ["SBS1", "SBS2", "SBS3", "SBS4", "SBS5", "SBS13"]
    lines = ["\t".join(["Type", *names])]
    for idx, label in enumerate(labels):
        row = [label] + [
            str(((idx * 7 + col * 13 + 3) % 23) / 23.0) for col in range(len(names))
        ]
        lines.append("\t".join(row))
    (directory / "COSMIC_v3.4_SBS_GRCh37.txt").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )
    cosmic.configure_catalogue_dir(directory)
    yield
    cosmic.configure_catalogue_dir(None)


def test_index_page_is_served(running_server):
    conn = _connection(running_server)
    conn.request("GET", "/")
    resp = conn.getresponse()
    assert resp.status == 200
    assert b"<html" in resp.read().lower()
    conn.close()


def test_forged_host_header_is_rejected(running_server):
    port = running_server.server_address[1]
    conn = http.client.HTTPConnection("127.0.0.1", port, timeout=5)
    conn.putrequest("GET", "/", skip_host=True)
    conn.putheader("Host", "evil.example.com")
    conn.endheaders()
    resp = conn.getresponse()
    assert resp.status == 400
    conn.close()


def test_cross_origin_post_is_rejected(running_server):
    body, content_type = _matrix_body(_sbs6_tsv())
    conn = _connection(running_server)
    headers = {
        "Content-Type": content_type,
        "Content-Length": str(len(body)),
        "Origin": "http://evil.example.com",
    }
    conn.request("POST", "/api/matrix", body=body, headers=headers)
    resp = conn.getresponse()
    assert resp.status == 403
    conn.close()


def test_static_path_traversal_is_rejected(running_server):
    conn = _connection(running_server)
    conn.request("GET", "/static/../../../../etc/passwd")
    resp = conn.getresponse()
    assert resp.status == 404
    conn.close()


def test_matrix_upload_and_render_round_trip(running_server):
    port = running_server.server_address[1]
    origin = f"http://127.0.0.1:{port}"
    body, content_type = _matrix_body(_sbs6_tsv())

    conn = _connection(running_server)
    conn.request(
        "POST",
        "/api/matrix",
        body=body,
        headers={
            "Content-Type": content_type,
            "Content-Length": str(len(body)),
            "Origin": origin,
        },
    )
    resp = conn.getresponse()
    assert resp.status == 200
    manifest = json.loads(resp.read())
    assert manifest["context"] == "SBS6"
    assert manifest["samples"] == ["S1", "S2"]
    conn.close()

    conn = _connection(running_server)
    conn.request("GET", "/api/render?sample=S1&unit=counts")
    resp = conn.getresponse()
    assert resp.status == 200
    figure = json.loads(resp.read())
    assert figure["kind"] == "standard"
    conn.close()


@pytest.mark.parametrize(
    ("matrix_text", "expected_context", "sample", "expected_warnings"),
    [
        (_id89_r_export_tsv, "ID89", "ID89Signature", 0),
        (_id476_r_export_tsv, "ID476", "Complete476", 2),
    ],
)
def test_r_style_id_matrix_upload_and_render_round_trip(
    running_server, matrix_text, expected_context, sample, expected_warnings
):
    port = running_server.server_address[1]
    origin = f"http://127.0.0.1:{port}"
    body, content_type = _matrix_body(matrix_text())

    conn = _connection(running_server)
    conn.request(
        "POST",
        "/api/matrix",
        body=body,
        headers={
            "Content-Type": content_type,
            "Content-Length": str(len(body)),
            "Origin": origin,
        },
    )
    response = conn.getresponse()
    assert response.status == 200
    manifest = json.loads(response.read())
    conn.close()

    assert manifest["context"] == expected_context
    assert manifest["samples"] == [sample]
    assert manifest["labels_count"] == int(expected_context[2:])
    assert len(manifest["warnings"]) == expected_warnings
    assert running_server.dashboard_matrix.values is None

    conn = _connection(running_server)
    conn.request("GET", f"/api/render?sample={sample}&unit=counts")
    response = conn.getresponse()
    assert response.status == 200
    figure = json.loads(response.read())
    conn.close()

    assert figure["kind"] == "standard"
    if expected_context == "ID476":
        assert "83 of 476" in manifest["warnings"][0]
        assert "Empty" in manifest["warnings"][1]


def test_render_without_a_prior_upload_is_rejected(running_server):
    conn = _connection(running_server)
    conn.request("GET", "/api/render?sample=S1")
    resp = conn.getresponse()
    assert resp.status == 409
    conn.close()


def test_unexpected_render_error_returns_500_without_dropping_the_connection(
    running_server, monkeypatch
):
    # A bug anywhere in render_sample (anything that isn't MatrixError) must
    # not silently kill the connection -- that's indistinguishable from a
    # server crash to the browser ("Failed to fetch", no diagnostic info).
    import sigProfilerPlotting.web.server as server_module

    def _boom(*args, **kwargs):
        raise RuntimeError("boom")

    monkeypatch.setattr(server_module, "render_sample", _boom)
    port = running_server.server_address[1]
    origin = f"http://127.0.0.1:{port}"
    body, content_type = _matrix_body(_sbs6_tsv())
    conn = _connection(running_server)
    conn.request(
        "POST",
        "/api/matrix",
        body=body,
        headers={
            "Content-Type": content_type,
            "Content-Length": str(len(body)),
            "Origin": origin,
        },
    )
    resp = conn.getresponse()
    resp.read()
    conn.close()

    conn = _connection(running_server)
    conn.request("GET", "/api/render?sample=S1")
    resp = conn.getresponse()
    assert resp.status == 500
    payload = json.loads(resp.read())
    assert "boom" in payload["error"]
    conn.close()

    # The server itself must still be alive and answering other requests.
    conn = _connection(running_server)
    conn.request("GET", "/")
    resp = conn.getresponse()
    assert resp.status == 200
    resp.read()
    conn.close()


def test_malformed_matrix_upload_returns_actionable_error(running_server):
    port = running_server.server_address[1]
    origin = f"http://127.0.0.1:{port}"
    bad_tsv = "MutationType\tS1\nC>A\tnot-a-number\n"
    body, content_type = _matrix_body(bad_tsv)
    conn = _connection(running_server)
    conn.request(
        "POST",
        "/api/matrix",
        body=body,
        headers={
            "Content-Type": content_type,
            "Content-Length": str(len(body)),
            "Origin": origin,
        },
    )
    resp = conn.getresponse()
    assert resp.status == 400
    payload = json.loads(resp.read())
    assert "error" in payload
    conn.close()


def test_upload_with_a_missing_content_length_is_rejected(running_server):
    port = running_server.server_address[1]
    origin = f"http://127.0.0.1:{port}"
    body, content_type = _matrix_body(_sbs6_tsv())
    conn = _connection(running_server)
    conn.putrequest("POST", "/api/matrix", skip_host=False)
    conn.putheader("Content-Type", content_type)
    conn.putheader("Origin", origin)
    conn.endheaders()
    conn.send(body)
    resp = conn.getresponse()
    assert resp.status == 411
    conn.close()


def test_upload_with_a_negative_content_length_is_rejected(running_server):
    port = running_server.server_address[1]
    origin = f"http://127.0.0.1:{port}"
    body, content_type = _matrix_body(_sbs6_tsv())
    conn = _connection(running_server)
    conn.putrequest("POST", "/api/matrix", skip_host=False)
    conn.putheader("Content-Type", content_type)
    conn.putheader("Content-Length", "-5")
    conn.putheader("Origin", origin)
    conn.endheaders()
    resp = conn.getresponse()
    assert resp.status == 400
    conn.close()


def test_upload_exceeding_the_configured_size_limit_is_rejected(
    running_server, monkeypatch
):
    import sigProfilerPlotting.web.server as server_module

    monkeypatch.setattr(server_module, "MAX_UPLOAD_BYTES", 10)
    port = running_server.server_address[1]
    origin = f"http://127.0.0.1:{port}"
    body, content_type = _matrix_body(_sbs6_tsv())
    assert len(body) > 10
    conn = _connection(running_server)
    conn.request(
        "POST",
        "/api/matrix",
        body=body,
        headers={
            "Content-Type": content_type,
            "Content-Length": str(len(body)),
            "Origin": origin,
        },
    )
    resp = conn.getresponse()
    assert resp.status == 413
    conn.close()

    # The server itself must still be alive and answering other requests.
    conn = _connection(running_server)
    conn.request("GET", "/")
    resp = conn.getresponse()
    assert resp.status == 200
    resp.read()
    conn.close()


def test_a_failed_upload_does_not_leak_its_temp_file(running_server):
    import sigProfilerPlotting.web.server as server_module

    port = running_server.server_address[1]
    origin = f"http://127.0.0.1:{port}"
    bad_tsv = "MutationType\tS1\nC>A\tnot-a-number\n"
    body, content_type = _matrix_body(bad_tsv)

    captured = {}
    real_stream = server_module._stream_upload_to_tempfile

    def spy(rfile, length):
        path = real_stream(rfile, length)
        captured["path"] = path
        return path

    server_module._stream_upload_to_tempfile = spy
    try:
        conn = _connection(running_server)
        conn.request(
            "POST",
            "/api/matrix",
            body=body,
            headers={
                "Content-Type": content_type,
                "Content-Length": str(len(body)),
                "Origin": origin,
            },
        )
        resp = conn.getresponse()
        assert resp.status == 400
        resp.read()
        conn.close()
    finally:
        server_module._stream_upload_to_tempfile = real_stream

    assert "path" in captured
    import os

    assert not os.path.exists(captured["path"])


def test_uploading_a_replacement_matrix_cleans_up_the_previous_temp_file(
    running_server,
):
    port = running_server.server_address[1]
    origin = f"http://127.0.0.1:{port}"

    def upload(content):
        body, content_type = _matrix_body(content)
        conn = _connection(running_server)
        conn.request(
            "POST",
            "/api/matrix",
            body=body,
            headers={
                "Content-Type": content_type,
                "Content-Length": str(len(body)),
                "Origin": origin,
            },
        )
        resp = conn.getresponse()
        assert resp.status == 200
        resp.read()
        conn.close()

    upload(_sbs6_tsv())
    first_tempfile = running_server._matrix_tempfile
    assert first_tempfile is not None and os.path.exists(first_tempfile)

    upload(_sbs96_tsv())
    second_tempfile = running_server._matrix_tempfile
    assert second_tempfile != first_tempfile
    assert not os.path.exists(first_tempfile)
    assert os.path.exists(second_tempfile)


def test_replacement_waits_for_an_in_flight_lazy_render(running_server, monkeypatch):
    import sigProfilerPlotting.web.server as server_module

    port = running_server.server_address[1]
    origin = f"http://127.0.0.1:{port}"

    def upload(content):
        body, content_type = _matrix_body(content)
        conn = _connection(running_server)
        conn.request(
            "POST",
            "/api/matrix",
            body=body,
            headers={
                "Content-Type": content_type,
                "Content-Length": str(len(body)),
                "Origin": origin,
            },
        )
        response = conn.getresponse()
        status = response.status
        response.read()
        conn.close()
        return status

    assert upload(_sbs6_tsv()) == 200
    old_path = running_server._matrix_tempfile
    real_render = server_module.render_sample
    real_set_matrix = running_server.set_dashboard_matrix
    render_started = threading.Event()
    release_render = threading.Event()
    replacement_attempted = threading.Event()

    def blocking_render(*args, **kwargs):
        render_started.set()
        assert release_render.wait(timeout=5)
        return real_render(*args, **kwargs)

    monkeypatch.setattr(server_module, "render_sample", blocking_render)

    def observed_set_matrix(matrix, tempfile_path):
        replacement_attempted.set()
        return real_set_matrix(matrix, tempfile_path)

    monkeypatch.setattr(running_server, "set_dashboard_matrix", observed_set_matrix)
    outcomes = {}

    def render_request():
        conn = _connection(running_server)
        conn.request("GET", "/api/render?sample=S1&unit=counts")
        response = conn.getresponse()
        outcomes["render"] = response.status
        response.read()
        conn.close()

    def replacement_upload():
        outcomes["upload"] = upload(_sbs96_tsv())

    render_thread = threading.Thread(target=render_request)
    upload_thread = threading.Thread(target=replacement_upload)
    render_thread.start()
    assert render_started.wait(timeout=5)
    upload_thread.start()
    assert replacement_attempted.wait(timeout=5)

    # The replacement has reached the swap, but it cannot unlink while the lazy
    # renderer is still using the old backing file.
    upload_thread.join(timeout=0.1)
    assert upload_thread.is_alive()
    assert old_path is not None and os.path.exists(old_path)

    release_render.set()
    render_thread.join(timeout=5)
    upload_thread.join(timeout=5)
    assert outcomes == {"render": 200, "upload": 200}
    assert not os.path.exists(old_path)


def test_html_response_carries_a_content_security_policy(running_server):
    conn = _connection(running_server)
    conn.request("GET", "/")
    resp = conn.getresponse()
    csp = resp.getheader("Content-Security-Policy")
    resp.read()
    conn.close()
    assert csp is not None
    assert "default-src 'self'" in csp


def test_make_server_rejects_non_loopback_host():
    with pytest.raises(ValueError):
        make_server("0.0.0.0", 0)


def test_bind_failure_preserves_the_original_socket_error(monkeypatch):
    import sigProfilerPlotting.web.server as server_module

    def fail_bind(self):
        raise OSError("address already in use")

    monkeypatch.setattr(server_module.DashboardHTTPServer, "server_bind", fail_bind)
    with pytest.raises(OSError, match="address already in use"):
        make_server("127.0.0.1", 0)


def test_ipv6_loopback_server_uses_ipv6_and_accepts_bracketed_host():
    if not socket.has_ipv6:
        pytest.skip("IPv6 is unavailable on this platform")
    try:
        server = make_server("::1", 0)
    except OSError as exc:
        pytest.skip(f"IPv6 loopback cannot be bound on this runner: {exc}")
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        assert server.address_family == socket.AF_INET6
        port = server.server_address[1]
        conn = http.client.HTTPConnection("::1", port, timeout=5)
        conn.putrequest("GET", "/", skip_host=True)
        conn.putheader("Host", f"[::1]:{port}")
        conn.endheaders()
        response = conn.getresponse()
        assert response.status == 200
        response.read()
        conn.close()
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)


def test_shutdown_requires_matching_origin(running_server):
    conn = _connection(running_server)
    conn.request(
        "POST",
        "/api/shutdown",
        body=b"",
        headers={"Content-Length": "0", "Origin": "http://evil.example.com"},
    )
    resp = conn.getresponse()
    assert resp.status == 403
    resp.read()
    conn.close()

    # Rejected shutdown must not actually stop the server.
    conn = _connection(running_server)
    conn.request("GET", "/")
    resp = conn.getresponse()
    assert resp.status == 200
    resp.read()
    conn.close()


def test_identity_endpoint_reports_pid_and_token(running_server):
    import os

    conn = _connection(running_server)
    conn.request("GET", "/api/identity")
    resp = conn.getresponse()
    assert resp.status == 200
    body = json.loads(resp.read())
    assert body["pid"] == os.getpid()
    assert isinstance(body["token"], str) and len(body["token"]) >= 16
    conn.close()


def test_shutdown_requires_matching_token(running_server):
    port = running_server.server_address[1]
    origin = f"http://127.0.0.1:{port}"
    conn = _connection(running_server)
    conn.request(
        "POST",
        "/api/shutdown",
        body=b"",
        headers={
            "Content-Length": "0",
            "Origin": origin,
            "X-Dashboard-Token": "wrong-token",
        },
    )
    resp = conn.getresponse()
    assert resp.status == 403
    resp.read()
    conn.close()

    # Rejected shutdown must not actually stop the server.
    conn = _connection(running_server)
    conn.request("GET", "/")
    resp = conn.getresponse()
    assert resp.status == 200
    resp.read()
    conn.close()


def test_quit_button_shuts_down_the_running_server():
    # Exercises the real lifecycle the "Quit" button drives: a same-origin
    # POST to /api/shutdown must make serve_forever() actually return, the
    # way the launcher's own SIGTERM handler does, rather than merely
    # returning a 200 while the server keeps running.
    server = make_server("127.0.0.1", 0)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        port = server.server_address[1]
        origin = f"http://127.0.0.1:{port}"
        conn = _connection(server)
        conn.request("GET", "/api/identity")
        token = json.loads(conn.getresponse().read())["token"]
        conn.close()

        conn = _connection(server)
        conn.request(
            "POST",
            "/api/shutdown",
            body=b"",
            headers={
                "Content-Length": "0",
                "Origin": origin,
                "X-Dashboard-Token": token,
            },
        )
        resp = conn.getresponse()
        assert resp.status == 200
        payload = json.loads(resp.read())
        assert payload["status"] == "stopping"
        conn.close()

        thread.join(timeout=5)
        assert (
            not thread.is_alive()
        ), "serve_forever() should have returned after /api/shutdown"
    finally:
        server.server_close()


def test_cosmic_builds_endpoint_lists_configured_genome_builds(running_server):
    conn = _connection(running_server)
    conn.request("GET", "/api/cosmic/builds")
    resp = conn.getresponse()
    assert resp.status == 200
    body = json.loads(resp.read())
    assert body["builds"] == ["GRCh37"]
    assert len(body["catalogues"]) == 1
    info = body["catalogues"][0]
    assert info["build"] == "GRCh37"
    assert info["version"] == "3.4"
    assert info["signature_count"] == 6
    assert "checksum" in info and "source_path" in info
    assert body["catalogue_dir"] == info["source_path"].rsplit("/", 1)[0]
    conn.close()


def test_cosmic_builds_endpoint_reports_no_builds_without_a_configured_directory(
    running_server,
):
    from sigProfilerPlotting.web import cosmic

    cosmic.configure_catalogue_dir(None)
    conn = _connection(running_server)
    conn.request("GET", "/api/cosmic/builds")
    resp = conn.getresponse()
    assert resp.status == 200
    body = json.loads(resp.read())
    assert body["builds"] == []
    assert body["catalogues"] == []
    assert body["catalogue_dir"] is None
    conn.close()


def test_render_with_compare_cosmic_returns_the_comparison_payload(running_server):
    port = running_server.server_address[1]
    origin = f"http://127.0.0.1:{port}"
    body, content_type = _matrix_body(_sbs96_tsv())

    conn = _connection(running_server)
    conn.request(
        "POST",
        "/api/matrix",
        body=body,
        headers={
            "Content-Type": content_type,
            "Content-Length": str(len(body)),
            "Origin": origin,
        },
    )
    resp = conn.getresponse()
    assert resp.status == 200
    resp.read()
    conn.close()

    conn = _connection(running_server)
    conn.request("GET", "/api/render?sample=S1&unit=counts&compare_cosmic=GRCh37")
    resp = conn.getresponse()
    assert resp.status == 200
    figure = json.loads(resp.read())
    assert figure["kind"] == "standard"
    compare = figure["cosmic_compare"]
    assert compare["build"] == "GRCh37"
    assert len(compare["matches"]) == 5
    conn.close()


def test_render_with_an_unknown_cosmic_build_is_rejected(running_server):
    port = running_server.server_address[1]
    origin = f"http://127.0.0.1:{port}"
    body, content_type = _matrix_body(_sbs96_tsv())

    conn = _connection(running_server)
    conn.request(
        "POST",
        "/api/matrix",
        body=body,
        headers={
            "Content-Type": content_type,
            "Content-Length": str(len(body)),
            "Origin": origin,
        },
    )
    resp = conn.getresponse()
    resp.read()
    conn.close()

    conn = _connection(running_server)
    conn.request("GET", "/api/render?sample=S1&compare_cosmic=not-a-real-build")
    resp = conn.getresponse()
    assert resp.status == 400
    conn.close()


def _png_bytes(size=(40, 20), mode="RGB", color="white"):
    import io

    from PIL import Image

    buffer = io.BytesIO()
    Image.new(mode, size, color).save(buffer, "PNG")
    return buffer.getvalue()


def test_export_pdf_wraps_a_png_in_a_one_page_pdf(running_server):
    import io

    from pypdf import PdfReader

    port = running_server.server_address[1]
    origin = f"http://127.0.0.1:{port}"
    body = _png_bytes()

    conn = _connection(running_server)
    conn.request(
        "POST",
        "/api/export/pdf",
        body=body,
        headers={
            "Content-Type": "image/png",
            "Content-Length": str(len(body)),
            "Origin": origin,
        },
    )
    resp = conn.getresponse()
    assert resp.status == 200
    assert resp.getheader("Content-Type") == "application/pdf"
    assert resp.getheader("Content-Disposition") == 'attachment; filename="plot.pdf"'
    pdf_bytes = resp.read()
    conn.close()

    assert pdf_bytes.startswith(b"%PDF")
    reader = PdfReader(io.BytesIO(pdf_bytes))
    assert len(reader.pages) == 1


def test_export_pdf_flattens_transparency_onto_white(running_server):
    import io

    from PIL import Image
    from pypdf import PdfReader

    port = running_server.server_address[1]
    origin = f"http://127.0.0.1:{port}"
    body = _png_bytes(mode="RGBA", color=(0, 0, 0, 0))  # fully transparent

    conn = _connection(running_server)
    conn.request(
        "POST",
        "/api/export/pdf",
        body=body,
        headers={
            "Content-Type": "image/png",
            "Content-Length": str(len(body)),
            "Origin": origin,
        },
    )
    resp = conn.getresponse()
    assert resp.status == 200
    pdf_bytes = resp.read()
    conn.close()

    reader = PdfReader(io.BytesIO(pdf_bytes))
    page = reader.pages[0]
    resources = page.get("/Resources", {}).get_object()
    xobject = next(iter(resources.get("/XObject", {}).get_object().values())).get_object()
    # A transparent source must not be embedded as-is (which would blacken
    # whatever showed through in a viewer that defaults to a black canvas);
    # it must have been flattened onto an opaque RGB image first.
    assert xobject.get("/ColorSpace") in ("/DeviceRGB", None, "/DeviceGray")
    image = Image.open(io.BytesIO(xobject.get_data()))
    assert image.mode in ("RGB", "L", "CMYK")


def test_export_pdf_requires_matching_origin(running_server):
    body = _png_bytes()

    conn = _connection(running_server)
    conn.request(
        "POST",
        "/api/export/pdf",
        body=body,
        headers={
            "Content-Type": "image/png",
            "Content-Length": str(len(body)),
            "Origin": "http://evil.example.com",
        },
    )
    resp = conn.getresponse()
    assert resp.status == 403
    conn.close()


def test_export_pdf_rejects_malformed_image_data(running_server):
    port = running_server.server_address[1]
    origin = f"http://127.0.0.1:{port}"
    body = b"not a real png"

    conn = _connection(running_server)
    conn.request(
        "POST",
        "/api/export/pdf",
        body=body,
        headers={
            "Content-Type": "image/png",
            "Content-Length": str(len(body)),
            "Origin": origin,
        },
    )
    resp = conn.getresponse()
    assert resp.status == 400
    assert b"Could not build a PDF" in resp.read()
    conn.close()


def test_export_pdf_rejects_an_oversized_image(running_server, monkeypatch):
    from sigProfilerPlotting.web import server as server_module

    monkeypatch.setattr(server_module, "MAX_EXPORT_PNG_BYTES", 10)
    port = running_server.server_address[1]
    origin = f"http://127.0.0.1:{port}"
    body = _png_bytes()
    assert len(body) > 10

    conn = _connection(running_server)
    conn.request(
        "POST",
        "/api/export/pdf",
        body=body,
        headers={
            "Content-Type": "image/png",
            "Content-Length": str(len(body)),
            "Origin": origin,
        },
    )
    resp = conn.getresponse()
    assert resp.status == 413
    conn.close()
