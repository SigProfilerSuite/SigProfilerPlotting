import http.client
import io
import json
import threading

import pytest

from sigProfilerPlotting.web.server import make_server

SBS6_LABELS = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]


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
