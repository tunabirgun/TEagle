"""WSL backend contracts for 4.1.0 (lane A6). Hermetic: WSL, subprocess and the elevated launcher are
mocked; the production shell scripts are executed in a local POSIX bash against fake tools in a temp HOME."""
from __future__ import annotations

import hashlib
import math
import os
import shutil
import subprocess
import sys
import time

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "app", "backend"))
import envcheck
from teagle_core import wsl

_RESOLVE = wsl.resolve_distro                      # conftest stubs the module attribute; keep the real function
GIB = 1073741824


def _bash():
    cands = []
    if os.name == "nt":
        cands += [r"C:\Program Files\Git\bin\bash.exe", r"C:\Program Files (x86)\Git\bin\bash.exe"]
    cands.append(shutil.which("bash"))
    path = next((p for p in cands if p and os.path.isfile(p) and "WindowsApps" not in p
                 and "System32" not in p), None)
    if not path:
        pytest.skip("a POSIX bash is required to execute the production WSL scripts")
    return path


def _run_bash(script, home, timeout=90):
    env = os.environ.copy()
    env["HOME"] = os.fspath(home).replace("\\", "/")
    return subprocess.run([_bash(), "-s"], input=script, text=True, cwd=home, env=env,
                          capture_output=True, timeout=timeout)


def _fake_exe(path, body):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("#!/usr/bin/env bash\n" + body, encoding="utf-8", newline="\n")
    os.chmod(path, 0o755)
    return path


# ---------- wsl_backend-2: a dead session is not a clean zero-hit result ----------
def _annotate_env(monkeypatch, script_result):
    monkeypatch.setattr(wsl, "resolve_species", lambda s: {"ok": True})
    monkeypatch.setattr(wsl, "env_status", lambda: {"wsl2": True, "ready": True, "repeatmasker": "4.2.4",
                                                     "dfam": True, "dfam_uncurated": False,
                                                     "dfam_library": {"version": "4.0", "partitions": ["x"]}})
    monkeypatch.setattr(wsl, "_wsl", lambda *a, **k: (0, "", ""))
    monkeypatch.setattr(wsl, "_wsl_script", script_result)


def test_annotate_dead_session_is_not_a_sealed_zero_hit_result(monkeypatch):
    _annotate_env(monkeypatch, lambda s, timeout=90: (
        4294967295, "", "The Windows Subsystem for Linux instance has terminated."))
    r = wsl.annotate(">x\nACGT", species="Homo sapiens")
    assert r["ok"] is False
    assert "nothing was searched" in r["error"]


def test_annotate_nonzero_session_rc_fails_even_with_a_zero_exit_marker(monkeypatch):
    _annotate_env(monkeypatch, lambda s, timeout=90: (1, "EXIT 0\n", ""))
    assert wsl.annotate(">x\nACGT", species="Homo sapiens")["ok"] is False


def test_annotate_clean_zero_hit_run_stays_ok(monkeypatch):
    _annotate_env(monkeypatch, lambda s, timeout=90: (0, "EXIT 0\n", ""))
    r = wsl.annotate(">x\nACGT", species="Homo sapiens")
    assert r["ok"] is True and r["n_hits"] == 0


# ---------- wsl_backend-12: a failed lineage probe is not "0 families" ----------
@pytest.mark.parametrize("result", [
    (4294967295, "", "The Windows Subsystem for Linux instance has terminated."),
    (4, "TEAGLE_FAMDB_FAILED 2\n", ""),
    (0, "", ""),
])
def test_lineage_family_count_fails_closed(monkeypatch, result):
    monkeypatch.setattr(wsl, "_wsl_script", lambda s, **k: result)
    r = wsl.dfam_lineage_families("Drosophila melanogaster")
    assert r["ok"] is False and "families" not in r


def test_lineage_family_count_reports_a_genuine_count(monkeypatch):
    monkeypatch.setattr(wsl, "_wsl_script", lambda s, **k: (0, "399\n", ""))
    assert wsl.dfam_lineage_families("Drosophila melanogaster") == {
        "ok": True, "species": "Drosophila melanogaster", "families": 399,
        "include_uncurated": False, "uncurated_available": None}


def test_zero_family_count_for_an_ambiguous_name_is_not_reported_as_a_count(monkeypatch):
    def script(s, **k):
        if " lineage -ad " in s:
            return 0, "Ambiguous search term, found 141 results\n", ""
        return 0, "0\n", ""
    monkeypatch.setattr(wsl, "_wsl_script", script)
    r = wsl.dfam_lineage_families("drosophila")
    assert r["ok"] is False and "141" in r["error"] and "families" not in r


def test_lineage_count_script_propagates_famdb_failure(monkeypatch, tmp_path):
    captured = []
    monkeypatch.setattr(wsl, "_wsl_script", lambda s, **k: (captured.append(s) or (0, "1\n", "")))
    wsl.dfam_lineage_families("Drosophila melanogaster")
    _fake_exe(tmp_path / "bin" / "micromamba", "echo 'usage: famdb' >&2\nexit 2\n")
    run = _run_bash(captured[0], tmp_path)
    assert run.returncode != 0 and "TEAGLE_FAMDB_FAILED" in run.stdout


# ---------- wsl_backend-3: distro selection ----------
def _table(rows):
    lines = ["  NAME            STATE           VERSION"]
    for name, state, default in rows:
        lines.append(f"{'*' if default else ' '} {name:<15} {state:<15} 2")
    return ("\r\n".join(lines) + "\r\n").encode("utf-16-le")


def _fake_wsl_run(rows, usable, probes=None):
    class R:
        def __init__(self, rc, out=b""):
            self.returncode, self.stdout, self.stderr = rc, out, b""

    def run(cmd, **kw):
        if cmd[:3] == ["wsl.exe", "-l", "-v"]:
            return R(0, _table(rows))
        name = cmd[cmd.index("-d") + 1]
        if probes is not None:
            probes.append(name)
        return R(0 if usable.get(name) else 127)
    return run


@pytest.fixture
def distro_env(monkeypatch, tmp_path):
    monkeypatch.setattr(wsl, "_DISTRO_FILE", str(tmp_path / "wsl_distro.txt"))
    monkeypatch.setattr(wsl, "_distro_cache", None)
    return tmp_path


def test_container_runtime_distro_is_never_the_backend(monkeypatch, distro_env):
    rows = [("docker-desktop", "Running", True), ("Ubuntu-24.04", "Stopped", False)]
    monkeypatch.setattr(wsl.subprocess, "run", _fake_wsl_run(rows, {"docker-desktop": False, "Ubuntu-24.04": True}))
    assert _RESOLVE() == "Ubuntu-24.04"


def test_docker_only_machine_has_no_usable_distro(monkeypatch, distro_env):
    rows = [("docker-desktop", "Running", True), ("docker-desktop-data", "Stopped", False)]
    probes = []
    monkeypatch.setattr(wsl.subprocess, "run", _fake_wsl_run(rows, {}, probes))
    assert _RESOLVE() is None
    assert probes == []                                    # excluded names are not even started


def test_default_distro_without_bash_or_python_is_skipped(monkeypatch, distro_env):
    rows = [("Alpine", "Running", True), ("Ubuntu", "Stopped", False)]
    monkeypatch.setattr(wsl.subprocess, "run", _fake_wsl_run(rows, {"Alpine": False, "Ubuntu": True}))
    assert _RESOLVE() == "Ubuntu"


def test_qualifying_default_is_preferred_and_alone_is_probed(monkeypatch, distro_env):
    rows = [("Archive", "Stopped", False), ("Ubuntu", "Stopped", True)]
    probes = []
    monkeypatch.setattr(wsl.subprocess, "run", _fake_wsl_run(rows, {"Archive": True, "Ubuntu": True}, probes))
    assert _RESOLVE() == "Ubuntu" and probes == ["Ubuntu"]


def test_remembered_distro_beats_a_changed_default(monkeypatch, distro_env):
    rows = [("Archive", "Stopped", False), ("Ubuntu", "Stopped", True)]
    monkeypatch.setattr(wsl.subprocess, "run", _fake_wsl_run(rows, {"Archive": True, "Ubuntu": True}))
    wsl.set_distro("Archive")
    assert open(wsl._DISTRO_FILE, encoding="utf-8").read() == "Archive"
    assert _RESOLVE() == "Archive"


def test_remembered_distro_that_no_longer_qualifies_falls_back(monkeypatch, distro_env):
    rows = [("Archive", "Stopped", False), ("Ubuntu", "Stopped", True)]
    open(wsl._DISTRO_FILE, "w", encoding="utf-8").write("Archive")
    monkeypatch.setattr(wsl.subprocess, "run", _fake_wsl_run(rows, {"Archive": False, "Ubuntu": True}))
    assert _RESOLVE() == "Ubuntu"


def test_every_wsl_call_names_the_chosen_distro(monkeypatch):
    seen = []
    monkeypatch.setattr(wsl, "resolve_distro", lambda: "Ubuntu-24.04")

    class R:
        returncode, stdout, stderr = 0, b"ok", b""
    monkeypatch.setattr(wsl.subprocess, "run", lambda cmd, **k: (seen.append(cmd) or R()))
    wsl._wsl("echo ok")
    wsl._wsl_script("echo ok")
    assert all(c[:3] == ["wsl.exe", "-d", "Ubuntu-24.04"] for c in seen) and len(seen) == 2


def test_environment_panel_does_not_call_docker_desktop_a_backend(monkeypatch):
    monkeypatch.setattr(envcheck.shutil, "which", lambda n: "wsl.exe")

    class R:
        stdout = ("  NAME            STATE           VERSION\r\n"
                  "* docker-desktop   Running         2\r\n")
    monkeypatch.setattr(envcheck.subprocess, "run", lambda *a, **k: R())
    assert envcheck.probe_backends()["wsl2"] == "no distro"
    R.stdout += "  Ubuntu-24.04    Stopped         2\r\n"
    assert envcheck.probe_backends()["wsl2"] == "available"


# ---------- wsl_backend-4: diagnose a failed liveness probe before advising ----------
def _available(monkeypatch, outcome):
    monkeypatch.setattr(wsl, "resolve_distro", lambda: "Ubuntu-24.04")
    monkeypatch.setattr(wsl, "_distro_cache", "Ubuntu-24.04")
    monkeypatch.setattr(wsl, "_DISTRO_FILE", os.devnull)

    def fake(script, stdin=None, timeout=600):
        if isinstance(outcome, BaseException):
            raise outcome
        return outcome
    monkeypatch.setattr(wsl, "_wsl", fake)
    return wsl.available()


def test_cold_start_timeout_is_reported_as_still_starting(monkeypatch):
    r = _available(monkeypatch, subprocess.TimeoutExpired("wsl.exe", 30))
    assert r["wsl2"] is False and r["reason"] == "timeout"


def test_missing_bash_is_an_unsuitable_distro_not_a_broken_one(monkeypatch):
    r = _available(monkeypatch, (127, "", "/bin/sh: 1: exec: bash: not found"))
    assert r["reason"] == "no_bash"


def test_terminated_instance_is_transient(monkeypatch):
    r = _available(monkeypatch, (4294967295, "", "The Windows Subsystem for Linux instance has terminated."))
    assert r["reason"] == "starting"


def test_registration_failure_is_confirmed_from_the_wsl_error_code(monkeypatch):
    r = _available(monkeypatch, (4294967295, "", "Wsl/Service/CreateInstance/0x80070002"))
    assert r["reason"] == "start_failed"


def _status_for(monkeypatch, reason, tmp_path, log=None, distro="Ubuntu-24.04"):
    monkeypatch.setattr(wsl.appdirs, "user_data_dir", lambda: str(tmp_path))
    if log is not None:
        (tmp_path / wsl._WIN_WSL_LOG).write_text(log, encoding="utf-8")
    monkeypatch.setattr(wsl, "available", lambda: {"wsl2": False, "distro": distro, "reason": reason,
                                                    "error": "x"})
    return wsl.components_status()


def _w(res):
    return next(c for c in res["components"] if c["key"] == "wsl2")


@pytest.mark.parametrize("reason", ["timeout", "starting", "no_bash", None])
def test_unconfirmed_failures_never_advise_unregister(monkeypatch, tmp_path, reason):
    text = (_w(_status_for(monkeypatch, reason, tmp_path))["guide"]).lower()
    assert "unregister" not in text


def test_timeout_guide_says_wait_and_refresh(monkeypatch, tmp_path):
    w = _w(_status_for(monkeypatch, "timeout", tmp_path))
    assert "starting" in w["detail"].lower() and "refresh" in w["guide"].lower()


def test_confirmed_start_failure_warns_that_unregister_deletes_the_distro_data(monkeypatch, tmp_path):
    guide = _w(_status_for(monkeypatch, "start_failed", tmp_path))["guide"].lower()
    assert "unregister ubuntu-24.04" in guide and "permanently" in guide and "home" in guide


def test_pending_reboot_is_decided_by_the_windows_log_not_dialog_state(monkeypatch, tmp_path):
    res = _status_for(monkeypatch, "start_failed", tmp_path, log="x\n[teagle] DONE-WSL 0 (restart)\n")
    w = _w(res)
    assert res["reboot_pending"] is True and "restart windows" in w["detail"].lower()
    assert "unregister" not in w["guide"].lower() and "unregister" not in w["detail"].lower()


def test_install_log_older_than_the_last_boot_is_not_a_pending_reboot(monkeypatch, tmp_path):
    monkeypatch.setattr(wsl, "_boot_time", lambda: time.time() + 60)
    res = _status_for(monkeypatch, "start_failed", tmp_path, log="x\n[teagle] DONE-WSL 0 (restart)\n")
    assert res["reboot_pending"] is False


# ---------- wsl_backend-5: the Windows-side install is reported separately ----------
def test_windows_install_flag_is_separate_from_the_in_wsl_lock(monkeypatch, tmp_path):
    res = _status_for(monkeypatch, None, tmp_path, log="[teagle] launching\n", distro=None)
    assert res["installing"] is False and res["wsl2_installing"] is True


def test_windows_install_log_without_a_marker_goes_stale(monkeypatch, tmp_path):
    res = _status_for(monkeypatch, None, tmp_path, log="[teagle] launching\n", distro=None)
    p = tmp_path / wsl._WIN_WSL_LOG
    old = time.time() - 3 * 3600
    os.utime(p, (old, old))
    res = wsl.components_status()
    assert res["wsl2_installing"] is False and res["wsl2_install_stale"] is True


def test_windows_install_log_from_before_the_last_boot_is_stale(monkeypatch, tmp_path):
    _status_for(monkeypatch, None, tmp_path, log="[teagle] launching\n", distro=None)
    monkeypatch.setattr(wsl, "_boot_time", lambda: time.time() + 60)
    assert wsl.components_status()["wsl2_installing"] is False


def _dialog():
    pytest.importorskip("PySide6")
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    native = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "app", "native")
    if native not in sys.path:
        sys.path.insert(0, native)
    from PySide6.QtWidgets import QApplication
    QApplication.instance() or QApplication([])
    import install_dialog
    d = install_dialog.InstallDialog()
    d._ops = []
    d.engine.submit = lambda op, body=None, key=None: d._ops.append(op)
    return d, install_dialog


_WSL2_ROW = {"key": "wsl2", "name": "WSL2", "desc": "x", "ok": False, "detail": "not installed",
             "installable": True}


def test_reopened_dialog_attaches_to_a_running_windows_install():
    d, _ = _dialog()
    d._render_components({"wsl2": False, "installing": False, "wsl2_installing": True,
                          "components": [dict(_WSL2_ROW)]})
    assert d._wsl2_installing is True and d._poll.isActive() and d._busy
    d._ops.clear()
    d._tick()
    assert d._ops[0] == "wsl2_install_log"                    # the Windows log, not the in-WSL log
    d._poll.stop()


def test_reopened_dialog_reads_the_pending_reboot_from_the_backend():
    d, _ = _dialog()
    d._render_components({"wsl2": False, "installing": False, "reboot_pending": True,
                          "components": [dict(_WSL2_ROW, detail="'Ubuntu' registered but won't start")]})
    assert d._wsl2_reboot_pending is True
    assert "restart" in d._rows["wsl2"]["detail"].text().lower()
    assert "restart" in d.statusLine.text().lower()


def test_stale_windows_install_releases_the_dialog():
    d, _ = _dialog()
    d._wsl2_installing = True
    d._poll.start()
    d._set_busy(True)
    d._render_components({"wsl2": False, "installing": False, "wsl2_installing": False,
                          "wsl2_install_stale": True, "components": [dict(_WSL2_ROW)]})
    assert d._wsl2_installing is False and not d._poll.isActive() and not d._busy
    assert "no longer running" in d.statusLine.text()


# ---------- wsl_backend-15: first-open probe failure schedules a retry ----------
def test_first_probe_failure_retries_then_points_to_refresh(monkeypatch):
    d, mod = _dialog()
    calls = []

    class FakeTimer:
        @staticmethod
        def singleShot(ms, *args):                          # (ms, fn) or (ms, context, fn)
            calls.append((ms, args[-1]))
    monkeypatch.setattr(mod, "QTimer", FakeTimer)
    comps = [dict(_WSL2_ROW, detail="probe failed")]
    d._render_components({"error": "timed out", "components": comps})
    assert calls and calls[0][0] >= 1000 and "retrying" in d.statusLine.text()
    assert "wsl2" in d._rows                                   # rows are drawn, not an empty window
    for _ in range(5):
        d._render_components({"error": "timed out", "components": comps})
    assert len(calls) <= 3 and "refresh" in d.statusLine.text().lower()


def test_successful_probe_resets_the_retry_budget(monkeypatch):
    d, mod = _dialog()
    monkeypatch.setattr(mod, "QTimer", type("T", (), {"singleShot": staticmethod(lambda *a: None)}))
    for _ in range(4):
        d._render_components({"error": "x"})
    d._render_components({"wsl2": True, "installing": False, "ready": False, "components": []})
    calls = []
    monkeypatch.setattr(mod, "QTimer", type("T", (), {"singleShot": staticmethod(lambda ms, *args: calls.append(ms))}))
    d._render_components({"error": "x"})
    assert calls


# ---------- wsl_backend-6: components install independently ----------
def _script_for(monkeypatch, steps, keys):
    monkeypatch.setattr(wsl, "_windows_backing_free_bytes", lambda _d=None: None)
    for k, body in steps.items():
        monkeypatch.setitem(wsl._STEP, k, body)
    return wsl._build_script(keys)


def test_a_failing_step_does_not_block_independent_steps(monkeypatch, tmp_path):
    script = _script_for(monkeypatch, {
        "minimap2": 'echo "[teagle] STEP minimap2 START"\nfail "solver died"\n',
        "dfam_root": 'echo "[teagle] STEP dfam_root START"\necho "[teagle] STEP dfam_root OK"\n',
        "dfam_curated": 'echo "[teagle] STEP dfam_curated START"\necho "[teagle] STEP dfam_curated OK"\n',
    }, ["minimap2", "dfam_root", "dfam_curated"])
    run = _run_bash(script, tmp_path)
    out = run.stdout
    assert "STEP dfam_root OK" in out and "STEP dfam_curated OK" in out
    assert run.returncode != 0 and "[teagle] DONE" not in out
    failed = [l for l in out.splitlines() if l.startswith("[teagle] FAILED")]
    assert len(failed) == 1 and "minimap2" in failed[0]        # one terminal line: mid-run notes must not end polling


def test_dependents_of_a_failed_step_are_skipped_not_run(monkeypatch, tmp_path):
    ran = tmp_path / "ran.txt"
    body = lambda k: f'echo "[teagle] STEP {k} START"\necho {k} >> "{ran.as_posix()}"\n'
    script = _script_for(monkeypatch, {
        "micromamba": 'echo "[teagle] STEP micromamba START"\nfail "no download"\n',
        "repeatmasker": body("repeatmasker"), "famdb_conf": body("famdb_conf"), "dfam_root": body("dfam_root"),
    }, ["micromamba", "repeatmasker", "dfam_root", "famdb_conf"])
    run = _run_bash(script, tmp_path)
    assert ran.read_text().split() == ["dfam_root"]
    assert "STEP repeatmasker SKIPPED" in run.stdout and "STEP famdb_conf SKIPPED" in run.stdout
    last = [l for l in run.stdout.splitlines() if l.startswith("[teagle] FAILED")][-1]
    assert "micromamba" in last


def test_all_steps_succeeding_still_ends_with_done(monkeypatch, tmp_path):
    script = _script_for(monkeypatch, {"minimap2": 'echo "[teagle] STEP minimap2 OK"\n'}, ["minimap2"])
    run = _run_bash(script, tmp_path)
    assert run.returncode == 0 and "[teagle] DONE" in run.stdout and "[teagle] FAILED" not in run.stdout


# ---------- wsl_backend-7: completed vs part-finished annotation ----------
def _annot_fixture(tmp_path, spsig="lib-a"):
    acc = "GCF_000000000.1"
    gdir = tmp_path / "teagle_genomes" / acc
    adir = gdir / "annot"
    adir.mkdir(parents=True)
    fasta = ">chr1\n" + "A" * 200 + "\n"
    sha = hashlib.sha256(fasta.encode()).hexdigest()
    (gdir / "genome.fna").write_bytes(fasta.encode())
    (gdir / "meta.txt").write_text(f"accession={acc}\ntarget=genome.fna\nsha256={sha}\ntarget_sha256={sha}\n"
                                   f"n_seqs=1\nbytes={len(fasta)}\nfasta_bytes={len(fasta)}\n", encoding="utf-8")
    (gdir / ".done").touch()
    (adir / ".split_done").touch()
    out = adir / "chunk_0001.out"
    out.write_text("100 1.0 0.0 0.0 chr1 1 50 (50) + teA LTR/Gypsy\n", encoding="utf-8")
    (adir / "chunk_0001.out.sha256").write_text(hashlib.sha256(out.read_bytes()).hexdigest() + "\n")
    (adir / ".runsig").write_text(f"genome_sha256={sha}|species={spsig}|sens=|chunkbp=1000000")
    return gdir, adir, sha, acc


def _annot_script(acc, spsig):
    return (wsl._ANNOT_SCRIPT.replace("__ACC__", acc).replace("__PA__", "1").replace("__CHUNKBP__", "1000000")
            .replace("__SENS__", "").replace("__SPSIG__", spsig).replace("__SP__", "")
            .replace("__HOST_FREE_B__", "-1"))


def test_a_successful_run_is_marked_complete_and_a_rerun_with_other_settings_says_so(tmp_path):
    gdir, adir, sha, acc = _annot_fixture(tmp_path)
    ok = _run_bash(_annot_script(acc, "lib-a"), tmp_path)
    assert "ANNOT_OK" in ok.stdout, ok.stdout + ok.stderr
    assert (adir / ".complete").exists()
    other = _run_bash(_annot_script(acc, "lib-b"), tmp_path)
    assert other.returncode != 0
    assert "completed annotation" in other.stdout and "different settings" in other.stdout
    assert "part-finished" not in other.stdout


def test_an_unfinished_run_is_still_called_part_finished(tmp_path):
    gdir, adir, sha, acc = _annot_fixture(tmp_path)
    other = _run_bash(_annot_script(acc, "lib-b"), tmp_path)
    assert "part-finished" in other.stdout and "different settings" in other.stdout


def test_genome_annotate_reports_the_state_behind_a_settings_refusal(monkeypatch):
    monkeypatch.setattr(wsl, "env_status", lambda: {"wsl2": True, "ready": True, "repeatmasker": "4.2.4",
                                                     "dfam": True, "dfam_library": {"partitions": ["a"]}})
    monkeypatch.setattr(wsl, "resolve_species", lambda s: {"ok": True})
    monkeypatch.setattr(wsl, "_wsl_script", lambda s, **k: (
        1, "FAIL a completed annotation of this genome exists with different settings (x) than requested\n", ""))
    r = wsl.genome_annotate("GCF_000000000.1", species="Homo sapiens")
    assert r["ok"] is False and r["annotation_state"] == "complete"
    assert "different settings" in r["error"]


# ---------- wsl_backend-8: disk sizing counts the FASTA, not the 2bit ----------
def test_budget_sizes_a_2bit_cache_from_the_fasta_it_must_re_download(monkeypatch):
    monkeypatch.setattr(wsl, "_wsl_script", lambda *a, **k: (0, f"cores=8\nmemkb={16 * 1048576}\navailb={20 * GIB}\n", ""))
    monkeypatch.setattr(wsl, "_windows_backing_free_bytes", lambda _d=None: None)
    fasta = 3 * GIB
    b = wsl.annotate_budget(fasta // 4, fasta_bytes=fasta, fasta_cached=False)
    assert b["disk_needed_gb"] == math.ceil(2 * fasta / GIB) + 2          # zip+FASTA peak, then FASTA+chunk copy
    cached = wsl.annotate_budget(fasta // 4, fasta_bytes=fasta, fasta_cached=True)
    assert cached["disk_needed_gb"] == math.ceil(fasta / GIB) + 2
    assert b["disk_ok"] is True
    big = wsl.annotate_budget(0, fasta_bytes=10 * GIB, fasta_cached=False)
    assert big["disk_needed_gb"] == 22 and big["disk_ok"] is False


def test_prepare_records_the_source_fasta_size_in_the_cache_metadata():
    assert "fasta_bytes=%s" in wsl._PREP_SCRIPT


def test_annotation_checks_free_space_before_re_downloading_the_fasta(tmp_path):
    acc = "GCF_000000000.1"
    gdir = tmp_path / "teagle_genomes" / acc
    gdir.mkdir(parents=True)
    (gdir / "genome.2bit").write_bytes(b"x")
    (gdir / "meta.txt").write_text(f"accession={acc}\ntarget=genome.2bit\nsha256={'a' * 64}\n"
                                   f"target_sha256={'b' * 64}\nn_seqs=1\nbytes=1\nfasta_bytes={10 ** 15}\n")
    (gdir / ".done").touch()
    called = tmp_path / "datasets-called"
    _fake_exe(tmp_path / "micromamba" / "envs" / "te" / "bin" / "datasets", f"touch '{called.as_posix()}'\nexit 0\n")
    run = _run_bash(_annot_script(acc, "lib-a"), tmp_path)
    assert "FAIL insufficient usable disk" in run.stdout and not called.exists()


# ---------- wsl_backend-10: FamDB check is graded on its exit status ----------
def _probe_text(famdb):
    return ("=RM=\nRepeatMasker version 4.2.4\n=MM=\n2.31-r1\n=FAMDB=\n" + famdb +
            "=FILES=\n=OPTFILES=\n=SCAN=\n=FATOTWOBIT=\nabsent\n=VRNA=\nabsent\n")


def _famdb_check(monkeypatch, famdb):
    monkeypatch.setattr(wsl, "available", lambda: {"wsl2": True})
    monkeypatch.setattr(wsl, "_wsl_script", lambda s, timeout=180: (0, _probe_text(famdb), ""))
    return next(c for c in wsl.integrity_check()["checks"] if c["name"] == "FamDB loads")


def test_famdb_error_text_mentioning_version_is_not_a_healthy_library(monkeypatch):
    c = _famdb_check(monkeypatch, "FAMDB_RC=1\nTraceback: unsupported file format version\n")
    assert c["ok"] is False


def test_famdb_nonzero_exit_fails_even_with_a_header(monkeypatch):
    assert _famdb_check(monkeypatch, "FAMDB_RC=1\nDatabase: Dfam\nVersion: 4.0\n")["ok"] is False


def test_famdb_clean_exit_with_header_passes(monkeypatch):
    assert _famdb_check(monkeypatch, "FAMDB_RC=0\nDatabase : Dfam\nVersion : 4.0\n")["ok"] is True


def test_integrity_probe_records_the_famdb_exit_status():
    assert "FAMDB_RC=" in wsl._INTEGRITY_PROBE


# ---------- wsl_backend-11: scan checksum probe is bounded, scaled and caught ----------
def test_scan_checksum_probe_timeout_is_a_plain_error_with_room_for_a_large_genome(monkeypatch):
    monkeypatch.setattr(wsl, "available", lambda: {"wsl2": True})
    calls = []

    def script(s, timeout=90):
        calls.append(timeout)
        raise subprocess.TimeoutExpired("wsl.exe", timeout)
    monkeypatch.setattr(wsl, "_wsl_script", script)
    r = wsl.genome_scan("GCF_000001405.40", "q\tACGT\tACGT\n")
    assert r["ok"] is False and "timed out" in r["error"] and "Traceback" not in r["error"]
    assert calls[-1] >= 600                       # 30 s was too short to hash a multi-GB genome read cold


def test_genome_prepare_tool_probe_timeout_is_a_plain_error(monkeypatch):
    monkeypatch.setattr(wsl, "available", lambda: {"wsl2": True})

    def inline(*a, **k):
        raise subprocess.TimeoutExpired("wsl.exe", 20)
    monkeypatch.setattr(wsl, "_wsl", inline)
    r = wsl.genome_prepare("GCF_000001405.40")
    assert r["ok"] is False and "did not answer" in r["error"] and "Traceback" not in r["error"]


# ---------- wsl_backend-13: the shared 'te' environment is only adopted when TEagle owns it ----------
def _env_fixture(tmp_path, with_pkgs=()):
    env = tmp_path / "micromamba" / "envs" / "te"
    (env / "conda-meta").mkdir(parents=True)
    for p in with_pkgs:
        (env / "conda-meta" / p).write_text("{}")
    calls = tmp_path / "mm_calls.log"
    _fake_exe(tmp_path / "bin" / "micromamba",
              'echo "$@" >> "$HOME/mm_calls.log"\n'
              'case "$1" in --version) echo 2.8.1;; run) exit 1;; '
              'create) mkdir -p "$HOME/micromamba/envs/te/conda-meta";; esac\nexit 0\n')
    return env, calls


def test_a_users_own_te_environment_is_not_adopted_or_overwritten(monkeypatch, tmp_path):
    env, calls = _env_fixture(tmp_path, with_pkgs=["numpy-1.26.0-py312.json"])
    script = _script_for(monkeypatch, {}, ["repeatmasker"])
    run = _run_bash(script, tmp_path)
    log = calls.read_text() if calls.exists() else ""
    assert "install" not in log, log
    assert not (env / ".teagle-owned").exists()
    assert run.returncode != 0 and "not created by TEagle" in run.stdout


def test_a_legacy_teagle_environment_is_adopted_and_marked(monkeypatch, tmp_path):
    env, calls = _env_fixture(tmp_path, with_pkgs=[f"repeatmasker-{wsl._CONDA_PINS['repeatmasker']}-pl5.json"])
    run = _run_bash(_script_for(monkeypatch, {}, ["repeatmasker"]), tmp_path)
    assert (env / ".teagle-owned").exists()
    assert "install" in calls.read_text()


def test_a_new_environment_is_marked_as_teagle_owned(monkeypatch, tmp_path):
    env = tmp_path / "micromamba" / "envs" / "te"
    _fake_exe(tmp_path / "bin" / "micromamba",
              'case "$1" in --version) echo 2.8.1;; run) exit 1;; '
              'create) mkdir -p "$HOME/micromamba/envs/te/conda-meta";; esac\nexit 0\n')
    _run_bash(_script_for(monkeypatch, {}, ["repeatmasker"]), tmp_path)
    assert (env / ".teagle-owned").exists()


def test_a_different_micromamba_already_in_home_bin_is_kept_as_a_backup(monkeypatch, tmp_path):
    _fake_exe(tmp_path / "bin" / "micromamba", 'case "$1" in --version) echo 1.5.0;; esac\nexit 0\n')
    _fake_exe(tmp_path / ".local" / "bin" / "micromamba",
              f'case "$1" in --version) echo {wsl._MICROMAMBA_VERSION};; esac\nexit 0\n')
    run = _run_bash(_script_for(monkeypatch, {}, ["micromamba"]), tmp_path)
    assert run.returncode == 0, run.stdout + run.stderr
    backup = tmp_path / "bin" / "micromamba.pre-teagle"
    assert backup.exists() and "1.5.0" in backup.read_text()
    assert wsl._MICROMAMBA_VERSION in (tmp_path / "bin" / "micromamba").read_text()


# ---------- wsl_backend-14/19: elevated install ----------
def test_elevated_install_has_no_interactive_fallback_and_stages_no_script():
    exe, params = wsl._wsl2_elevated_command(r"C:\Users\x\AppData\Local\TEagle\wsl_win_install.log")
    assert exe.lower().endswith("cmd.exe")
    installs = [part for part in params.split("--install")[1:]]
    assert installs and all(p.lstrip().startswith("-d Ubuntu --no-launch") for p in installs)
    assert "DONE-WSL" in params and ".bat" not in params.lower()


@pytest.mark.skipif(os.name != "nt", reason="runs the generated cmd.exe command line")
def test_elevated_command_logs_the_install_status_and_terminal_marker(tmp_path):
    stub = tmp_path / "wslstub.cmd"
    stub.write_text("@echo off\r\necho stub %*\r\nif \"%1\"==\"--install\" exit /b 3\r\nexit /b 0\r\n")
    log = tmp_path / "wsl log.txt"
    exe, params = wsl._wsl2_elevated_command(str(log), wsl_exe=str(stub))
    done = subprocess.run(f'{exe} {params}', capture_output=True, timeout=60)
    text = log.read_bytes().decode("utf-8", "ignore").replace("\x00", "")
    assert "stub --install -d Ubuntu --no-launch" in text
    assert "stub --set-default-version 2" in text and "stub --status" in text
    assert "DONE-WSL 3" in text and text.index("--status") < text.index("DONE-WSL")


@pytest.mark.skipif(os.name != "nt", reason="Windows-only elevation path")
@pytest.mark.parametrize("is_admin", [True, False])
def test_install_wsl2_launches_without_writing_a_batch_file_and_flags_over_the_shoulder(monkeypatch, tmp_path, is_admin):
    monkeypatch.setattr(wsl.appdirs, "user_data_dir", lambda: str(tmp_path))
    monkeypatch.setattr(wsl, "_wsl2_thread", None)
    monkeypatch.setattr(wsl, "_account_is_administrator", lambda: is_admin)
    launched = []
    monkeypatch.setattr(wsl, "_shell_execute_runas", lambda exe, params: launched.append((exe, params)) or 42)
    r = wsl.install_wsl2()
    wsl._wsl2_thread.join(timeout=5)
    assert r["started"] is True and launched and not (tmp_path / "install_wsl.bat").exists()
    first = (tmp_path / wsl._WIN_WSL_LOG).read_text(encoding="utf-8")
    assert ("different account" in first.lower()) is (not is_admin)
    assert bool(r.get("warning")) is (not is_admin)


# ---------- wsl_backend-16: resumed progress counts from the finished chunks ----------
def test_resumed_run_numbers_its_progress_from_the_chunks_already_done(tmp_path):
    gdir, adir, sha, acc = _annot_fixture(tmp_path)
    pend = adir / "chunk_0002.fa"
    pend.write_bytes(b">chr2\nACGT\n")
    (adir / "chunk_0002.input.sha256").write_text(hashlib.sha256(pend.read_bytes()).hexdigest() + "\n")
    _fake_exe(tmp_path / "bin" / "micromamba", 'chunk="${@: -1}"\nprintf "" > "$chunk.out"\nexit 0\n')
    run = _run_bash(_annot_script(acc, "lib-a"), tmp_path)
    assert "RepeatMasker chunk 2/2" in run.stdout and "CHUNK_DONE 2/2" in run.stdout, run.stdout
    assert "chunk 1/2" not in run.stdout


# ---------- wsl_backend-17 / 18: wording ----------
def test_uncurated_root_row_quotes_the_measured_coverage_table():
    desc = next(c[3] for c in wsl._COMP_META if c[0] == "dfam_unc_root")
    assert "very few families" not in desc
    assert wsl.curated_coverage_sentence() in desc


@pytest.mark.parametrize("st,needle", [
    ({"wsl2": False, "repeatmasker": None, "dfam": False, "dfam_library": None}, "WSL is not installed"),
    ({"wsl2": True, "repeatmasker": None, "dfam": True, "dfam_library": {"v": 1}}, "RepeatMasker is not installed"),
    ({"wsl2": True, "repeatmasker": "4.2.4", "dfam": False, "dfam_library": None}, "Dfam"),
    ({"wsl2": True, "repeatmasker": "4.2.4", "dfam": True, "dfam_library": None}, "FamDB could not open"),
])
def test_not_ready_error_names_the_failing_condition(monkeypatch, st, needle):
    monkeypatch.setattr(wsl, "resolve_species", lambda s: {"ok": True})
    monkeypatch.setattr(wsl, "env_status", lambda: dict(st, ready=False))
    r = wsl.annotate(">x\nACGT", species="Homo sapiens")
    assert r["ok"] is False and needle in r["error"] and "None" not in r["error"]


# ---------- ui_family_splice_prov-2: the searched family set is what RepeatMasker searches ----------
_RM_LOG = ("Search Engine: rmblast [ 2.14.1+ ]\nSpecies/Taxa Search:\n  Homo sapiens [NCBI Taxonomy ID: 9606]\n"
           "  Lineage: root; cellular organisms; Eukaryota\n")


def _run_annotate_through_bash(monkeypatch, tmp_path, species, rm_body):
    _fake_exe(tmp_path / "bin" / "micromamba", rm_body)
    monkeypatch.setattr(wsl, "resolve_species", lambda s: {"ok": True})
    monkeypatch.setattr(wsl, "env_status", lambda: {"wsl2": True, "ready": True, "repeatmasker": "4.2.4",
                                                     "dfam": True, "dfam_uncurated": False,
                                                     "dfam_library": {"version": "4.0", "partitions": ["x"]}})
    monkeypatch.setattr(wsl, "_wsl", lambda *a, **k: (0, "", ""))

    def run(script, timeout=90):
        rid = __import__("re").search(r"/tmp/(teagle_[0-9a-f]+)", script).group(1)
        d = tmp_path / "tmpdir" / rid
        script = script.replace(f"/tmp/{rid}", d.as_posix())
        pre = f'mkdir -p "{d.as_posix()}"; printf ">x\\nACGT\\n" > "{d.as_posix()}/q.fa"\n'
        p = _run_bash(pre + script, tmp_path)
        return p.returncode, p.stdout, p.stderr
    monkeypatch.setattr(wsl, "_wsl_script", run)
    return wsl.annotate(">x\nACGT", species=species)


_RM_OK = ('printf "%s" \'' + _RM_LOG.replace("'", "") + '\'\n'
          'for a in "$@"; do q="$a"; done\n'
          'printf "  100 1.0 0.0 0.0 x 1 4 (0) + Alu SINE/Alu 1 4 (0) 1\\n" > "$q.out"\nexit 0\n')


def test_no_organism_is_labelled_with_the_lineage_repeatmasker_defaults_to(monkeypatch, tmp_path):
    r = _run_annotate_through_bash(monkeypatch, tmp_path, None, _RM_OK)
    assert r["ok"] is True, r
    assert r["species"] == "Homo sapiens (RepeatMasker default lineage)"
    assert "all installed" not in r["species"] and r["species_defaulted"] is True
    assert r["lineage_taxid"] == "9606" and r["n_hits"] == 1


def test_no_organism_falls_back_to_the_documented_default_when_the_log_has_no_taxon_line(monkeypatch, tmp_path):
    r = _run_annotate_through_bash(monkeypatch, tmp_path, None,
                                   'for a in "$@"; do q="$a"; done\nprintf "" > "$q.out"\nexit 0\n')
    assert r["species"] == "Homo sapiens (RepeatMasker default lineage)"


def test_named_organism_keeps_its_name_and_the_resolved_taxon(monkeypatch, tmp_path):
    r = _run_annotate_through_bash(monkeypatch, tmp_path, "human", _RM_OK)
    assert r["species"] == "human" and r["lineage_searched"] == "Homo sapiens"


def test_unknown_species_advice_no_longer_claims_a_search_of_all_families(monkeypatch):
    monkeypatch.setattr(wsl, "_wsl_script", lambda s, timeout=90: (1, "No species found", ""))
    r = wsl.resolve_species("Notasp")
    assert r["ok"] is False and "all installed families" not in r["error"] and "Homo sapiens" in r["error"]


def test_genome_path_labels_the_default_lineage_too():
    assert "SPECIES_SEARCHED" in wsl._ANNOT_SCRIPT


# ---------- ui_family_splice_prov-3: transcript strand comes from minimap2's ts tag ----------
def _sam(flag, ts=None, cigar="50M200N50M"):
    tag = f"\tts:A:{ts}" if ts else ""
    return f"q1\t{flag}\tref\t101\t60\t{cigar}\t*\t0\t0\t*\t*\tNM:i:0{tag}\n"


@pytest.mark.parametrize("flag,ts,gene,align", [
    (0, "+", "+", "+"), (0, "-", "-", "+"), (16, "+", "-", "-"), (16, "-", "+", "-"),
])
def test_gene_strand_combines_alignment_orientation_with_the_ts_tag(flag, ts, gene, align):
    r = wsl._parse_sam_splice(_sam(flag, ts))
    assert r["strand"] == gene and r["alignment_strand"] == align and r["strand_determined"] is True


def test_spliced_alignment_without_a_ts_tag_does_not_claim_a_determined_strand():
    r = wsl._parse_sam_splice(_sam(16, None))
    assert r["strand_determined"] is False and r["alignment_strand"] == "-"


def test_antisense_transcript_reads_splice_sites_on_the_gene_strand(monkeypatch):
    # plus-strand text of a MINUS-strand gene: intron reads CT..AC; read as the gene it is GT..AG
    exon1, intron, exon2 = "A" * 50, "CT" + "T" * 196 + "AC", "C" * 50
    g = exon1 + intron + exon2
    sam = f"q1\t0\tref\t1\t60\t50M{len(intron)}N50M\t*\t0\t0\t*\t*\tNM:i:0\tts:A:-\n"
    monkeypatch.setattr(wsl, "available", lambda: {"wsl2": True})

    def inline(script, stdin=None, timeout=600):
        if "minimap2 -a" in script:
            return 0, sam, ""
        if "yes" in script:
            return 0, "yes", ""
        return 0, "2.31-r1", ""
    monkeypatch.setattr(wsl, "_wsl", inline)
    r = wsl.splice_align(f">g\n{g}\n", ">t\nACGT\n")
    assert r["ok"] and r["strand"] == "-" and r["canonical_introns"] == 1
    assert (r["introns"][0]["donor"], r["introns"][0]["acceptor"]) == ("GT", "AG")


# ---------- packaging_ci-2: a frozen build has nothing to install ----------
def test_frozen_build_does_not_report_install_needed_on_first_run(monkeypatch, tmp_path):
    monkeypatch.setattr(envcheck, "FROZEN", True)
    monkeypatch.setattr(envcheck, "STATE", str(tmp_path / "none.json"))
    reqs = [("psutil", "7.2.2")]
    monkeypatch.setattr(envcheck, "parse_requirements", lambda: reqs)
    monkeypatch.setattr(envcheck, "_installed", lambda n: "7.2.2")
    monkeypatch.setattr(envcheck, "_importable", lambda n: True)
    monkeypatch.setattr(envcheck, "probe_backends", lambda: {})
    rep = envcheck.check()
    assert rep["needs_install"] is False and rep["first_run"] is False and rep["bundled"] is True


def test_frozen_build_still_flags_a_genuinely_missing_dependency(monkeypatch, tmp_path):
    monkeypatch.setattr(envcheck, "FROZEN", True)
    monkeypatch.setattr(envcheck, "STATE", str(tmp_path / "none.json"))
    monkeypatch.setattr(envcheck, "parse_requirements", lambda: [("psutil", "7.2.2")])
    monkeypatch.setattr(envcheck, "_installed", lambda n: None)
    monkeypatch.setattr(envcheck, "probe_backends", lambda: {})
    assert envcheck.check()["needs_install"] is True


def test_source_run_still_reports_first_run(monkeypatch, tmp_path):
    monkeypatch.setattr(envcheck, "FROZEN", False)
    monkeypatch.setattr(envcheck, "STATE", str(tmp_path / "none.json"))
    monkeypatch.setattr(envcheck, "parse_requirements", lambda: [("psutil", "7.2.2")])
    monkeypatch.setattr(envcheck, "_installed", lambda n: "7.2.2")
    monkeypatch.setattr(envcheck, "_importable", lambda n: True)
    monkeypatch.setattr(envcheck, "probe_backends", lambda: {})
    rep = envcheck.check()
    assert rep["first_run"] is True and rep["needs_install"] is True


# ---------- installer clean install removes the 'te' environment only when TEagle owns it ----------
def _clean_install_command(iss_text=None):
    import re
    text = iss_text if iss_text is not None else (
        __import__("pathlib").Path(__file__).resolve().parent.parent / "installer" / "teagle.iss").read_text(encoding="utf-8")
    return re.search(r'-- bash -lc "([^"]*)"\'', text).group(1)


def _run_clean_install(tmp_path, env_files):
    env = tmp_path / "micromamba" / "envs" / "te"
    for rel in env_files:
        (env / rel).parent.mkdir(parents=True, exist_ok=True)
        (env / rel).write_text("{}")
    (tmp_path / "micromamba" / "envs" / "teagle-vrna").mkdir(parents=True)
    (tmp_path / "teagle_genomes").mkdir()
    (tmp_path / "bin").mkdir()
    (tmp_path / "bin" / "micromamba").write_text("x")
    run = _run_bash(_clean_install_command(), tmp_path)
    assert run.returncode == 0, run.stdout + run.stderr
    return env


def test_clean_install_keeps_an_unowned_te_environment_and_the_shared_micromamba(tmp_path):
    env = _run_clean_install(tmp_path, ["conda-meta/numpy-1.26.0-py312.json"])
    assert (env / "conda-meta" / "numpy-1.26.0-py312.json").exists()
    assert (tmp_path / "bin" / "micromamba").exists()
    assert not (tmp_path / "teagle_genomes").exists()
    assert not (tmp_path / "micromamba" / "envs" / "teagle-vrna").exists()


def test_clean_install_removes_a_marked_te_environment(tmp_path):
    env = _run_clean_install(tmp_path, [".teagle-owned", "conda-meta/numpy-1.26.0-py312.json"])
    assert not env.exists() and (tmp_path / "bin" / "micromamba").exists()


def test_clean_install_removes_a_legacy_teagle_environment_by_its_pinned_repeatmasker(tmp_path):
    env = _run_clean_install(tmp_path, [f"conda-meta/repeatmasker-{wsl._CONDA_PINS['repeatmasker']}-pl5321.json"])
    assert not env.exists()


def test_clean_install_ownership_test_matches_the_backend_marker_and_pin():
    cmd = _clean_install_command()
    assert ".teagle-owned" in cmd and f"repeatmasker-{wsl._CONDA_PINS['repeatmasker']}-" in cmd
