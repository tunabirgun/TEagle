"""Bounded tests for the optional local HTTP transport."""
from __future__ import annotations

import json
import threading
import urllib.error
import urllib.request

import pytest

import server
from teagle_core import provenance


def test_held_miniprot_tier_is_not_exposed_by_the_http_api():
    assert "/api/miniprot" not in server.Handler._ROUTES


def test_serve_binds_ephemeral_loopback_socket_and_reports_primer3(monkeypatch, capsys):
    """Exercise the startup line after a real bind without entering the blocking server loop."""
    real_make_server = server.make_server
    state = {}

    def bounded_make_server(host="127.0.0.1", port=8765):
        srv = real_make_server(host, 0)
        state["address"] = srv.server_address
        srv.serve_forever = lambda: state.update(loop_entered=True)
        original_close = srv.server_close

        def close():
            original_close()
            state["closed"] = True

        srv.server_close = close
        return srv

    monkeypatch.setattr(server, "make_server", bounded_make_server)
    server.serve(host="127.0.0.1", port=0)

    assert state["address"][0] == "127.0.0.1"
    assert state["address"][1] > 0
    assert state["loop_entered"] is True
    assert state["closed"] is True
    output = capsys.readouterr().out
    assert f"TEagle backend on http://127.0.0.1:{state['address'][1]}" in output
    assert "http://127.0.0.1:0" not in output
    assert f"primer3 {server.primers.PRIMER3_VERSION}" in output


def test_api_json_is_literal_utf8_and_nonfinite_payload_fails_closed(monkeypatch):
    monkeypatch.setitem(server.Handler._ROUTES, "/api/unicode", lambda _body: {"record": "é–record"})
    monkeypatch.setitem(server.Handler._ROUTES, "/api/nonfinite", lambda _body: {"score": float("nan")})
    srv = server.make_server("127.0.0.1", 0)
    worker = threading.Thread(target=srv.serve_forever, daemon=True)
    worker.start()
    base = f"http://127.0.0.1:{srv.server_port}"
    try:
        request = urllib.request.Request(base + "/api/unicode", data=b"{}",
                                         headers={"Content-Type": "application/json"}, method="POST")
        with urllib.request.urlopen(request, timeout=5) as response:
            body = response.read()
            assert response.headers.get_content_type() == "application/json"
            assert response.headers.get_content_charset() == "utf-8"
        assert "é–record".encode("utf-8") in body and b"\\u00e9" not in body
        assert json.loads(body) == {"record": "é–record"}

        request = urllib.request.Request(base + "/api/nonfinite", data=b"{}",
                                         headers={"Content-Type": "application/json"}, method="POST")
        with pytest.raises(urllib.error.HTTPError) as failure:
            urllib.request.urlopen(request, timeout=5)
        assert failure.value.code == 500
        error_body = failure.value.read()
        assert b"NaN" not in error_body and b"Infinity" not in error_body
        assert json.loads(error_body)["error"].startswith("ValueError:")
    finally:
        srv.shutdown()
        srv.server_close()
        worker.join(timeout=5)


def test_browser_shaped_source_wire_reaches_engine_through_http_transport():
    accession = "M11240.1"
    sequence = "ACGT" * 30
    raw = {
        "accession": accession, "taxid": "7227", "source": "NCBI nuccore",
        "endpoint": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
        "retrievedUtc": "2026-08-21T00:00:00+00:00",
        "sourceUrl": f"https://www.ncbi.nlm.nih.gov/nuccore/{accession}", "length": len(sequence),
    }
    source = provenance.normalize_source_context(raw, [(accession, sequence)])
    payload = json.dumps({"sequence": f">{accession}\n{sequence}", "source": source}).encode()
    srv = server.make_server("127.0.0.1", 0)
    worker = threading.Thread(target=srv.serve_forever, daemon=True)
    worker.start()
    try:
        request = urllib.request.Request(
            f"http://127.0.0.1:{srv.server_port}/api/analyze", data=payload,
            headers={"Content-Type": "application/json"}, method="POST")
        with urllib.request.urlopen(request, timeout=20) as response:
            result = json.loads(response.read())
        assert result["records"][0]["id"] == accession
        assert result["records"][0]["identitySource"] == "validated_fetch"
        assert result["provenance"]["input"]["sourceUrl"] == raw["sourceUrl"]
        assert provenance.verify_manifest(result["provenance"])
    finally:
        srv.shutdown()
        srv.server_close()
        worker.join(timeout=5)
