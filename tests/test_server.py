"""Bounded tests for the optional local HTTP transport."""
from __future__ import annotations

import server


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
