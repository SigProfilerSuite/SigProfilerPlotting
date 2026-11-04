import http.client
import io
import json
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


def _multipart_body(filename, content):
    boundary = "----testboundary123"
    body = (
        f"--{boundary}\r\n"
        f'Content-Disposition: form-data; name="matrix"; filename="{filename}"\r\n'
        f"Content-Type: text/tab-separated-values\r\n\r\n"
    ).encode("utf-8")
    body += content.encode("utf-8")
    body += f"\r\n--{boundary}--\r\n".encode("utf-8")
    return body, f"multipart/form-data; boundary={boundary}"


def _sbs6_tsv():
    lines = ["MutationType\tS1\tS2"]
    for index, label in enumerate(SBS6_LABELS):
        lines.append(f"{label}\t{index}\t{index * 2}")
    return "\n".join(lines) + "\n"


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
    port = running_server.server_address[1]
    body, content_type = _multipart_body("matrix.tsv", _sbs6_tsv())
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
    body, content_type = _multipart_body("matrix.tsv", _sbs6_tsv())

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


def test_render_without_a_prior_upload_is_rejected(running_server):
    conn = _connection(running_server)
    conn.request("GET", "/api/render?sample=S1")
    resp = conn.getresponse()
    assert resp.status == 409
    conn.close()


def test_unexpected_render_error_returns_500_without_dropping_the_connection(running_server, monkeypatch):
    # A bug anywhere in render_sample (anything that isn't MatrixError) must
    # not silently kill the connection -- that's indistinguishable from a
    # server crash to the browser ("Failed to fetch", no diagnostic info).
    import sigProfilerPlotting.web.server as server_module

    def _boom(*args, **kwargs):
        raise RuntimeError("boom")

    monkeypatch.setattr(server_module, "render_sample", _boom)
    port = running_server.server_address[1]
    origin = f"http://127.0.0.1:{port}"
    body, content_type = _multipart_body("matrix.tsv", _sbs6_tsv())
    conn = _connection(running_server)
    conn.request(
        "POST", "/api/matrix", body=body,
        headers={"Content-Type": content_type, "Content-Length": str(len(body)), "Origin": origin},
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
    body, content_type = _multipart_body("bad.tsv", bad_tsv)
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


def test_shutdown_requires_matching_origin(running_server):
    conn = _connection(running_server)
    conn.request(
        "POST", "/api/shutdown", body=b"",
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
        conn.request(
            "POST", "/api/shutdown", body=b"",
            headers={"Content-Length": "0", "Origin": origin},
        )
        resp = conn.getresponse()
        assert resp.status == 200
        payload = json.loads(resp.read())
        assert payload["status"] == "stopping"
        conn.close()

        thread.join(timeout=5)
        assert not thread.is_alive(), "serve_forever() should have returned after /api/shutdown"
    finally:
        server.server_close()


def test_cosmic_builds_endpoint_lists_bundled_genome_builds(running_server):
    conn = _connection(running_server)
    conn.request("GET", "/api/cosmic/builds")
    resp = conn.getresponse()
    assert resp.status == 200
    body = json.loads(resp.read())
    assert body["builds"] == ["GRCh37", "GRCh38", "mm9", "mm10", "mm39", "rn6", "rn7"]
    conn.close()


def test_render_with_compare_cosmic_returns_the_comparison_payload(running_server):
    port = running_server.server_address[1]
    origin = f"http://127.0.0.1:{port}"
    body, content_type = _multipart_body("matrix.tsv", _sbs96_tsv())

    conn = _connection(running_server)
    conn.request(
        "POST", "/api/matrix", body=body,
        headers={"Content-Type": content_type, "Content-Length": str(len(body)), "Origin": origin},
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
    assert figure["kind"] == "sbs96_compare"
    assert figure["build"] == "GRCh37"
    assert len(figure["matches"]) == 5
    conn.close()


def test_render_with_an_unknown_cosmic_build_is_rejected(running_server):
    port = running_server.server_address[1]
    origin = f"http://127.0.0.1:{port}"
    body, content_type = _multipart_body("matrix.tsv", _sbs96_tsv())

    conn = _connection(running_server)
    conn.request(
        "POST", "/api/matrix", body=body,
        headers={"Content-Type": content_type, "Content-Length": str(len(body)), "Origin": origin},
    )
    resp = conn.getresponse()
    resp.read()
    conn.close()

    conn = _connection(running_server)
    conn.request("GET", "/api/render?sample=S1&compare_cosmic=not-a-real-build")
    resp = conn.getresponse()
    assert resp.status == 400
    conn.close()
