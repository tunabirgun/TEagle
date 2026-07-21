"""WSL2 install-from-app logic (elevated Windows-side path). The feature can't be end-to-end tested
where WSL is already present, so these lock down the decomposable parts: the elevated .bat contents,
the NUL-stripping Windows-log reader, the absent-vs-broken classification, and the dialog routing.
No real `wsl --install` / UAC is ever triggered."""
import os
import pytest
from teagle_core import wsl


def _wsl2(res):
    return [c for c in res["components"] if c["key"] == "wsl2"][0]


# ---------- elevated .bat ----------
def test_wsl2_bat_script_has_install_fallback_and_terminal_marker():
    bat = wsl._wsl2_bat_script(r"C:\tmp\wsl.log")
    assert "wsl.exe --install -d Ubuntu --no-launch" in bat      # primary: register distro, no OOBE
    assert "retrying: wsl --install" in bat                      # fallback for older wsl.exe
    assert "--set-default-version 2" in bat and "|| rem" in bat  # tolerant pre-reboot
    assert "DONE-WSL" in bat                                     # terminal marker the poller stops on


# ---------- Windows-side log reader ----------
def test_wsl2_install_log_strips_utf16_nuls(monkeypatch, tmp_path):
    monkeypatch.setattr(wsl.appdirs, "user_data_dir", lambda: str(tmp_path))
    p = os.path.join(str(tmp_path), wsl._WIN_WSL_LOG)
    with open(p, "wb") as f:                                     # wsl.exe emits UTF-16LE
        f.write(b"[teagle] start\n" + "Ubuntu".encode("utf-16-le") + b"\n[teagle] DONE-WSL 0\n")
    log = wsl.wsl2_install_log(50)
    assert "\x00" not in log and "Ubuntu" in log and "DONE-WSL" in log


def test_wsl2_installing_terminal_markers(monkeypatch, tmp_path):
    monkeypatch.setattr(wsl.appdirs, "user_data_dir", lambda: str(tmp_path))
    p = os.path.join(str(tmp_path), wsl._WIN_WSL_LOG)
    open(p, "w").write("[teagle] running...\n");            assert wsl._wsl2_installing() is True
    open(p, "w").write("...\n[teagle] DONE-WSL 0\n");        assert wsl._wsl2_installing() is False
    open(p, "w").write("...\n[teagle] FAILED - declined\n"); assert wsl._wsl2_installing() is False


# ---------- absent vs registered-but-broken ----------
def test_components_status_absent_offers_install(monkeypatch):
    monkeypatch.setattr(wsl, "available", lambda: {"wsl2": False, "distro": None, "error": "no distro"})
    w = _wsl2(wsl.components_status())
    assert "Install WSL" in w["guide"]
    assert w.get("installable") == (os.name == "nt")            # Windows-only capability


def test_components_status_broken_distro_offers_unregister(monkeypatch):
    monkeypatch.setattr(wsl, "available", lambda: {"wsl2": False, "distro": "Ubuntu-24.04", "error": "launch failed"})
    w = _wsl2(wsl.components_status())
    assert "unregister Ubuntu-24.04" in w["guide"]              # broken ext4.vhdx needs unregister+reinstall


def test_install_wsl2_is_windows_only(monkeypatch):
    if os.name == "nt":
        pytest.skip("would attempt a real elevated install on Windows")
    assert wsl.install_wsl2()["started"] is False               # non-Windows: refused, never launches


# ---------- dialog routing (no real install) ----------
def test_dialog_wsl2_button_and_install_all_routing():
    pytest.importorskip("PySide6")
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    import sys
    native = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "app", "native")
    if native not in sys.path:
        sys.path.insert(0, native)
    from PySide6.QtWidgets import QApplication
    QApplication.instance() or QApplication([])
    from install_dialog import InstallDialog
    dlg = InstallDialog()
    dlg._render_components({"wsl2": False, "installing": False, "components": [
        {"key": "wsl2", "name": "WSL2", "desc": "x", "ok": False, "detail": "not installed",
         "installable": True, "guide": "Click Install WSL..."}]})
    b = dlg._rows["wsl2"]["btn"]
    assert b.isEnabled() and b.text() == "Install WSL"          # live even though WSL is absent
    ops = []
    dlg.engine.submit = lambda op, body=None, key=None: ops.append(op)
    dlg._wsl2_ok = False; dlg._install_all()
    assert ops == ["wsl_install_wsl2"] and dlg._wsl2_installing is True   # absent -> elevated WSL installer
    ops.clear(); dlg._wsl2_installing = False; dlg._wsl2_ok = True; dlg._install_all()
    assert ops == ["wsl_install"]                              # present -> normal in-WSL stack (unchanged path)


def test_dialog_reboot_pending_not_destructive_unregister():
    """A successful install that needs a reboot must NOT advise the destructive `wsl --unregister`."""
    pytest.importorskip("PySide6")
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    import sys
    native = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "app", "native")
    if native not in sys.path:
        sys.path.insert(0, native)
    from PySide6.QtWidgets import QApplication
    QApplication.instance() or QApplication([])
    from install_dialog import InstallDialog
    dlg = InstallDialog()
    dlg.engine.submit = lambda *a, **k: None                   # don't fire the real _refresh() job
    # a completed elevated install reports DONE-WSL 0 -> sticky reboot-pending state (status flickers via _refresh,
    # then the next _render_components restores the reboot message — assert the flag here, the message after render)
    dlg._on_done("wsl2_log", {"log": "[teagle] installing WSL2 + Ubuntu\n[teagle] DONE-WSL 0 (restart Windows)\n"})
    assert dlg._wsl2_reboot_pending is True
    # the very next components probe shows the just-registered distro as won't-start (pre-reboot) with unregister guidance
    dlg._render_components({"wsl2": False, "installing": False, "components": [
        {"key": "wsl2", "name": "WSL2", "desc": "x", "ok": False,
         "detail": "'Ubuntu-24.04' registered but won't start", "installable": True,
         "guide": "In an Administrator PowerShell run:  wsl --unregister Ubuntu-24.04  then click Install WSL"}]})
    det = dlg._rows["wsl2"]["detail"].text().lower()
    assert "restart" in det and "unregister" not in det        # reboot advice, NOT destructive unregister
    assert dlg._rows["wsl2"]["btn"].isVisible() is False        # no button to wipe/reinstall the good install
    assert "restart" in dlg.statusLine.text().lower()          # reboot message persists through _render_components (C2)
    # once WSL actually comes up (post-reboot), the sticky state clears
    dlg._render_components({"wsl2": True, "ready": True, "installing": False, "components": [
        {"key": "wsl2", "name": "WSL2", "desc": "x", "ok": True, "detail": "Ubuntu-24.04"}]})
    assert dlg._wsl2_reboot_pending is False
    # a nonzero DONE-WSL is a soft failure, still reboot-pending
    dlg._on_done("wsl2_log", {"log": "[teagle] DONE-WSL 1\n"})
    assert dlg._wsl2_reboot_pending is True
