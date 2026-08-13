"""WSL2 install-from-app logic (elevated Windows-side path). The feature can't be end-to-end tested
where WSL is already present, so these lock down the decomposable parts: the elevated .bat contents,
the NUL-stripping Windows-log reader, the absent-vs-broken classification, and the dialog routing.
No real `wsl --install` / UAC is ever triggered."""
import os
from pathlib import Path
import pytest
from teagle_core import wsl


def _wsl2(res):
    return [c for c in res["components"] if c["key"] == "wsl2"][0]


def _assert_cleanup_target_contract(iss: str) -> None:
    cleanup = iss[iss.index("function CleanTeagleData(): Boolean;"):iss.index("procedure CurStepChanged")]
    capture = cleanup.index("HasRecordedDistro := TryReadRecordedWslDistro(Distro);")
    assert "WindowsCleanupOk := DelTree(ExpandConstant('{localappdata}\\TEagle')" in cleanup
    delete = cleanup.index("WindowsCleanupOk := DelTree(ExpandConstant('{localappdata}\\TEagle')")
    invoke = cleanup.index("WslCleanupOk := CleanWslBackend(Distro);")
    assert capture < invoke < delete
    assert "if HasRecordedDistro then\n    WslCleanupOk := CleanWslBackend(Distro);" in cleanup
    assert "SaveStringToFile(marker, Distro, False)" in cleanup
    assert "HasRecordedDistro and (not WslCleanupOk)" in cleanup
    assert "if not WindowsCleanupOk then" in cleanup
    assert "if MarkerRestored then" in cleanup
    assert "if MarkerRestored then\n      FailureText := FailureText + 'The optional WSL backend could not be cleaned." in cleanup
    assert "retry marker could not be restored" in cleanup
    assert "Result := WindowsCleanupOk and WslCleanupOk;" in cleanup

    backend = iss[iss.index("function CleanWslBackend"):iss.index("function CleanTeagleData")]
    assert "'-d \"' + Distro + '\" -- bash" in backend
    assert "if not Exec(" in backend and "if rc <> 0 then" in backend
    assert "Result := True;" in backend
    assert "distroArg" not in backend


# ---------- elevated .bat ----------
def test_wsl2_bat_script_has_install_fallback_and_terminal_marker():
    bat = wsl._wsl2_bat_script(r"C:\tmp\wsl.log")
    assert "wsl.exe --install -d Ubuntu --no-launch" in bat      # primary: register distro, no OOBE
    assert "retrying: wsl --install" in bat                      # fallback for older wsl.exe
    assert "--set-default-version 2" in bat and "|| rem" in bat  # tolerant pre-reboot
    assert "DONE-WSL" in bat                                     # terminal marker the poller stops on


def test_windows_clean_install_preserves_shared_micromamba_root():
    """The opt-in cleanup removes TEagle data without deleting other applications' conda envs."""
    root = Path(__file__).resolve().parents[1]
    iss = (root / "installer" / "teagle.iss").read_text(encoding="utf-8")
    clean = next(line for line in iss.splitlines() if "Exec(wsl, '-d" in line)
    assert "~/micromamba/envs/te " in clean
    assert "~/micromamba/envs/teagle-vrna" in clean
    assert "~/teagle_genomes" in clean and "~/teagle_libs" in clean
    assert "~/micromamba " not in clean
    assert "~/bin/micromamba" not in clean


def test_windows_cleanup_captures_recorded_distro_before_deleting_marker():
    """A recorded non-default distro remains the cleanup target after LocalAppData is removed."""
    root = Path(__file__).resolve().parents[1]
    iss = (root / "installer" / "teagle.iss").read_text(encoding="utf-8")
    _assert_cleanup_target_contract(iss)


def test_windows_cleanup_missing_or_invalid_marker_never_uses_default_distro():
    """No marker means no Linux-side deletion; malformed marker text cannot inject wsl.exe arguments."""
    root = Path(__file__).resolve().parents[1]
    iss = (root / "installer" / "teagle.iss").read_text(encoding="utf-8")
    reader = iss[iss.index("function TryReadRecordedWslDistro"):iss.index("function CleanWslBackend")]
    cleanup = iss[iss.index("function CleanTeagleData(): Boolean;"):iss.index("procedure CurStepChanged")]
    validator = iss[iss.index("function IsSafeDistroName"):iss.index("function TryReadRecordedWslDistro")]

    assert "if not FileExists(marker) then" in reader
    assert "Result := False;" in reader
    assert "if HasRecordedDistro then\n    WslCleanupOk := CleanWslBackend(Distro);" in cleanup
    assert "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789._-" in validator
    assert "default distro" not in iss.lower()


def test_windows_cleanup_promises_match_the_safe_skip_behavior():
    root = Path(__file__).resolve().parents[1]
    iss = (root / "installer" / "teagle.iss").read_text(encoding="utf-8")
    assert iss.count("valid recorded distribution") >= 3
    assert "failed cleanup retains" in iss
    assert "retained for a later retry" in iss
    assert "including the WSL annotation backend) will be removed" not in iss


def test_windows_cleanup_gate_rejects_injected_old_defects():
    """Prove the source gate fails if deletion comes first or cleanup loses its explicit distro."""
    root = Path(__file__).resolve().parents[1]
    iss = (root / "installer" / "teagle.iss").read_text(encoding="utf-8")
    capture = "  HasRecordedDistro := TryReadRecordedWslDistro(Distro);\n"
    invoke = "  if HasRecordedDistro then\n    WslCleanupOk := CleanWslBackend(Distro);\n"
    delete = "  WindowsCleanupOk := DelTree(ExpandConstant('{localappdata}\\TEagle'), True, True, True);\n"

    deleted_first = iss.replace(invoke + delete, delete + invoke, 1)
    with pytest.raises(AssertionError):
        _assert_cleanup_target_contract(deleted_first)

    implicit_default = iss.replace("'-d \"' + Distro + '\" -- bash", "'-- bash", 1)
    with pytest.raises(AssertionError):
        _assert_cleanup_target_contract(implicit_default)

    unchecked_exit = iss.replace("  if rc <> 0 then", "  if False then", 1)
    with pytest.raises(AssertionError):
        _assert_cleanup_target_contract(unchecked_exit)

    erased_retry = iss.replace("SaveStringToFile(marker, Distro, False)",
                               "SaveStringToFile(marker, '', False)", 1)
    with pytest.raises(AssertionError):
        _assert_cleanup_target_contract(erased_retry)

    unchecked_windows_delete = iss.replace(
        "WindowsCleanupOk := DelTree(ExpandConstant('{localappdata}\\TEagle'), True, True, True);",
        "DelTree(ExpandConstant('{localappdata}\\TEagle'), True, True, True);", 1)
    with pytest.raises(AssertionError):
        _assert_cleanup_target_contract(unchecked_windows_delete)

    false_restoration_claim = iss.replace("if MarkerRestored then\n      FailureText :=",
                                          "if True then\n      FailureText :=", 1)
    with pytest.raises(AssertionError):
        _assert_cleanup_target_contract(false_restoration_claim)


def test_inno_source_has_no_fallback_version_and_callers_check_cleanup_result():
    root = Path(__file__).resolve().parents[1]
    iss = (root / "installer" / "teagle.iss").read_text(encoding="utf-8")
    assert '#error "MyAppVersion must be injected by installer/build_installer.ps1"' in iss
    preamble = iss[:iss.index("#define MyAppName")]
    assert '#define MyAppVersion "' not in preamble
    assert iss.count("if not CleanTeagleData() then") == 2
    assert "RaiseException('Clean install stopped" in iss


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
    assert b.isEnabled() and b.text().upper() == "INSTALL WSL"  # live even though WSL is absent (labels are uppercased)
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
