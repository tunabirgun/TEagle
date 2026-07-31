"""Visible-progress and honest-label contract for the backend installer.

A user who starts the optional 3.9 GiB Dfam partition must be able to watch it move. Before this,
curl drew its own carriage-return meter: the whole 40-60 minute transfer was ONE line with no newline
in it, so `tail -n` handed the log panel a single 250 kB line and the download looked frozen. These
lock down the three parts that fix: the step reports whole lines, the reader normalises any CR meter
it still meets, and a component that is not installed says Install rather than Repair."""
import os
import sys

import pytest
from teagle_core import wsl


def _step(key="dfam_unc_euk"):
    return wsl._STEP[key]


# ---------- the download step reports progress as whole lines ----------
def test_dfam_step_silences_curl_meter_and_reports_lines():
    s = _step()
    assert "curl -sS -L --fail -C -" in s          # -s kills the CR meter, -S keeps real errors
    assert "sleep 15" in s and "kill -0 \"$dl_pid\"" in s   # a watcher polls the file while curl runs
    assert 'echo "[teagle] ' in s.split("dl_pid=$!", 1)[1]  # and reports on its own newline-terminated lines


def test_dfam_step_percentage_is_derived_from_a_head_not_the_resumed_get():
    """A resumed (-C -) GET answers 206 with Content-Length = REMAINING bytes. Reading the total off
    the transfer would make every percentage wrong, so the total comes from a separate HEAD."""
    s = _step()
    assert "curl -sIL" in s and "content-length" in s.lower()
    head, rest = s.split("total_b=$(curl -sIL", 1)
    assert "-C -" not in rest.split(")", 1)[0]     # the HEAD itself must not be a resumed request
    # a percentage is only printed when the total is known AND consistent with the bytes on disk
    assert '[ "$total_b" -gt 0 ] && [ "$now_b" -le "$total_b" ]' in s


def _bash(script: str):
    """Run a bash snippet and return (rc, stdout). Skips where no bash is available."""
    import shutil, subprocess, tempfile, os as _os
    bash = shutil.which("bash")
    if not bash:
        pytest.skip("no bash on PATH")
    fd, p = tempfile.mkstemp(suffix=".sh")
    with _os.fdopen(fd, "wb") as f:
        f.write(script.replace("\r\n", "\n").encode())
    try:
        r = subprocess.run([bash, p], capture_output=True, text=True, timeout=60)
        return r.returncode, r.stdout.strip()
    finally:
        _os.unlink(p)


@pytest.mark.parametrize("rc,total,have,expect", [
    (0, 4_000_000_000, 4_000_000_000, "complete"),   # curl clean -> complete, whatever the sizes say
    (0, 0, 123, "complete"),                          # curl clean with an unknown total is still complete
    (33, 4_000_000_000, 4_000_000_000, "complete"),   # HTTP 416 on a finished .gz: errored but all there
    (18, 4_000_000_000, 2_000_000_000, "partial"),    # cut short against a known total
    (18, 0, 2_000_000_000, "partial"),                # cut short with NO known total — the case that regressed
    (18, 0, 0, "partial"),
])
def test_transfer_state_calls_an_unknown_total_partial_not_complete(rc, total, have, expect):
    """Deleting a multi-GB partial costs the whole transfer again, so only positive evidence that the
    download finished may license it. The size HEAD can fail or a server can omit Content-Length; when
    the first version of this guard required a known total, that case fell through to the md5 gate and
    the partial was deleted despite the panel promising it resumes."""
    fn = wsl._PRELUDE[wsl._PRELUDE.index("transfer_state(){"):wsl._PRELUDE.index("# A truncated/corrupt")]
    out_rc, out = _bash(f"{fn}\ntransfer_state {rc} {total} {have}\n")
    assert out_rc == 0
    assert out == expect


def test_dfam_step_never_deletes_a_partial_and_still_checks_it():
    """A partial is checked anyway — a complete file whose size was never advertised must not loop
    forever — but the delete-on-mismatch gate is reached only on the complete path."""
    s = _step()
    i_state = s.index("transfer_state ")
    i_del = s.index("rm -f")
    assert i_state < i_del                                  # completeness is decided before any delete
    partial_branch = s[i_state:s.index('if [ "$verified" -eq 0 ]')]
    assert "rm -f" not in partial_branch                    # nothing removes the file on the partial path
    assert "md5sum -c -" in partial_branch                  # but it is still verified
    assert "partial file is KEPT" in partial_branch


def test_dfam_step_announces_the_silent_multi_minute_stages():
    s = _step()
    assert "verifying checksum" in s and "decompressing" in s   # md5 + gunzip on 3.9 GiB are minutes of silence
    # ...but the cost warning is conditional on the actual size: the root partition is 1.8 MB
    assert 'slow=""; [ "$now_b" -gt' in s


def test_env_creation_recovers_a_prefix_that_is_not_an_environment():
    """Installing a Dfam partition before the environment exists creates envs/te/share/... . micromamba
    then refuses that prefix ("Non-conda folder exists at prefix") for good, bricking the backend."""
    p = wsl._PRELUDE
    body = p[p.index("mm_create(){"):p.index("mm_install(){")]   # scope every claim to mm_create itself
    assert "Non-conda folder exists at prefix" in body           # the failure this recovers from, named
    i_stash = body.index('STASH="$ENV.stash')
    i_create = body.index('"$MM" create -y -n te', i_stash)
    assert i_stash < i_create                                    # set aside BEFORE create, or create aborts
    # a multi-GB library is moved back, never re-downloaded and never duplicated on a disk that may not hold two copies
    assert 'mv -f "$f" "$FAMDIR/"' in body and "cp -a" not in body
    # a failed create restores what was there rather than leaving the user with neither
    assert 'mv "$STASH" "$ENV"' in body


def test_integrity_probe_covers_every_dfam_partition_and_grades_only_the_required_ones():
    """A user who has just spent 40 minutes on the optional partition must see it named by the
    integrity check; its absence must not be graded as a fault, because it is optional."""
    assert set(wsl._DFAM_REQUIRED) | set(wsl._DFAM_OPTIONAL) == set(wsl._DFAM_FILES)   # derived, not restated
    assert set(wsl._DFAM_REQUIRED).isdisjoint(wsl._DFAM_OPTIONAL)
    assert all(k in wsl._ALL_STEPS for k in wsl._DFAM_REQUIRED)
    for k in wsl._DFAM_FILES:                                   # every partition file appears in the probe
        assert wsl._DFAM_FILES[k][0] in wsl._INTEGRITY_PROBE, k


def test_integrity_check_reports_optional_partitions_without_failing_on_them(monkeypatch):
    probe = ("=RM=\nRepeatMasker version 4.2.4\n=MM=\n2.31-r1302\n=FAMDB=\nVersion : 4.0\n"
             "=FILES=\npresent dfam40.0.h5 1\npresent dfam40.curated.consensus.0.h5 1\n"
             "=OPTFILES=\npresent dfam40.uncurated.consensus.0.h5\nabsent dfam40.uncurated.consensus.1.h5\n"
             "=SCAN=\nisPcr present datasets present\n")
    monkeypatch.setattr(wsl, "available", lambda: {"wsl2": True, "distro": "Ubuntu", "error": None})
    monkeypatch.setattr(wsl, "_wsl_script", lambda *a, **k: (0, probe, ""))
    r = wsl.integrity_check()
    files = [c for c in r["checks"] if c["name"] == "Dfam library files present"][0]
    assert r["ok"] is True and files["ok"] is True              # one optional partition missing is not a fault
    assert "uncurated.consensus.0.h5 present" in files["detail"]
    assert "uncurated.consensus.1.h5 not installed" in files["detail"]


def test_every_conda_install_step_creates_the_environment_through_the_recovery():
    """mm_create carries the non-conda-prefix recovery; a step that called mm_install without it would
    still abort with no way back."""
    for key in ("repeatmasker", "minimap2", "genomescan"):
        b = wsl._STEP[key]
        assert "mm_install" in b and "mm_create" in b, key
        assert b.index("mm_create") < b.index("mm_install"), key


def test_prelude_kills_background_helpers_on_exit():
    """Closing the app kills the WSL session and this script; an orphaned watcher would keep writing
    progress for a download that is no longer running."""
    assert "cleanup(){" in wsl._PRELUDE and "kill $BG" in wsl._PRELUDE
    assert "trap cleanup EXIT" in wsl._PRELUDE


@pytest.mark.parametrize("key", sorted(wsl._DFAM_FILES))
def test_every_dfam_step_is_valid_bash(key, tmp_path):
    """`bash -n` the real generated script — a quoting slip in the watcher would only surface mid-download."""
    import shutil, subprocess
    bash = shutil.which("bash")
    if not bash:
        pytest.skip("no bash on PATH")
    p = tmp_path / f"{key}.sh"
    p.write_bytes(wsl._build_script([key]).encode())
    r = subprocess.run([bash, "-n", str(p)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr


# ---------- the reader normalises any CR meter that still reaches it ----------
def test_install_log_collapses_a_carriage_return_meter_to_one_line(monkeypatch):
    """micromamba's solver still draws one, and it must not arrive as a single unbounded line."""
    meter = "\r".join(f"  {p} 3938M  {p}%" for p in range(0, 101))
    raw = "[teagle] START\n[teagle] STEP dfam_unc_euk START\n" + meter
    monkeypatch.setattr(wsl, "_wsl", lambda *a, **k: (0, raw, ""))
    out = wsl.install_log(200)
    assert "\r" not in out
    assert out.splitlines()[-1] == "  100 3938M  100%"     # the state the meter ended on, not its history
    assert len(out.splitlines()) == 3                      # two real lines plus the collapsed meter


def test_install_log_read_is_bounded_by_bytes(monkeypatch):
    """`tail -n` alone cannot bound a CR meter — it is one line however long it grows."""
    seen = {}
    monkeypatch.setattr(wsl, "_wsl", lambda s, **k: (seen.setdefault("cmd", s), "x", "")[1:] + ("",))
    wsl.install_log(200)
    assert "tail -c" in seen["cmd"]


def test_install_log_returns_the_tail_it_was_asked_for(monkeypatch):
    monkeypatch.setattr(wsl, "_wsl", lambda *a, **k: (0, "\n".join(str(i) for i in range(500)), ""))
    assert wsl.install_log(200).splitlines() == [str(i) for i in range(300, 500)]


# ---------- an absent component says Install, not Repair ----------
def _dialog():
    pytest.importorskip("PySide6")
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    native = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "app", "native")
    if native not in sys.path:
        sys.path.insert(0, native)
    from PySide6.QtWidgets import QApplication
    QApplication.instance() or QApplication([])
    from install_dialog import InstallDialog
    return InstallDialog()


def test_optional_partition_offers_install_before_it_exists():
    """The two uncurated partitions are never installed by default, so "Repair" was the first word a
    user met for a download they had not yet run."""
    dlg = _dialog()
    dlg._render_components({"wsl2": True, "installing": False, "ready": False, "components": [
        {"key": "dfam_unc_euk", "name": "Dfam uncurated · Eukaryota (optional)", "desc": "3.9 GiB",
         "ok": False, "detail": "not installed (optional)", "repairable": True},
        {"key": "dfam_root", "name": "Dfam 4.0 root library", "desc": "root", "ok": True,
         "detail": "present", "repairable": True}]})
    assert dlg._rows["dfam_unc_euk"]["btn"].text().upper() == "INSTALL"   # labels are uppercased
    # present: the same idempotent step re-runs and fixes what is missing — a repair. It does NOT
    # re-download a file already on disk, so it must not claim "Reinstall".
    assert dlg._rows["dfam_root"]["btn"].text().upper() == "REPAIR"


def test_dialog_attaches_to_an_install_that_was_already_running():
    """Opening the window while a download runs used to leave the log blank for its whole run: polling
    started only when the user clicked something in that session."""
    dlg = _dialog()
    sent = []
    dlg.engine.submit = lambda op, body=None, key=None: sent.append(op)
    core_ready = {"wsl2": True, "installing": True, "ready": True, "components": [
        {"key": "dfam_unc_euk", "name": "Dfam uncurated · Eukaryota (optional)", "desc": "3.9 GiB",
         "ok": False, "detail": "not installed (optional)", "repairable": True}]}
    dlg._render_components(core_ready)
    assert dlg._poll.isActive() and dlg._busy
    # the status must not read "ready" over the top of a running install
    assert "running" in dlg.statusLine.text() and "log" in dlg.statusLine.text()
    dlg._poll.stop()


def test_dialog_recovers_when_an_install_dies_without_a_terminal_marker():
    """A run killed outright — app closed, WSL shut down — writes no DONE/FAILED. Polling for one would
    leave the window disabled for good."""
    import time as _t
    dlg = _dialog()
    dlg.engine.submit = lambda op, body=None, key=None: None
    running = {"wsl2": True, "installing": True, "ready": False, "components": [
        {"key": "dfam_unc_euk", "name": "Dfam uncurated · Eukaryota (optional)", "desc": "3.9 GiB",
         "ok": False, "detail": "not installed (optional)", "repairable": True}]}
    dlg._render_components(running)
    assert dlg._poll.isActive() and dlg._busy
    stopped = dict(running, installing=False)
    dlg._render_components(stopped)                    # inside the grace window: a just-started run, keep polling
    assert dlg._poll.isActive()
    dlg._poll_since = _t.monotonic() - 60               # past it, with no marker in the log
    dlg._render_components(stopped)
    assert dlg._poll.isActive()                         # ONE lockless probe is not proof — see the test below
    dlg._render_components(stopped)
    assert not dlg._poll.isActive() and not dlg._busy
    assert "stopped before it finished" in dlg.statusLine.text()
    assert "resumes" in dlg.statusLine.text()           # and says the download is not lost


def test_a_successful_install_is_never_announced_as_dead():
    """The lock state and the log come from two separate WSL round trips taken at different moments. A
    run that finishes between the log read and the lock check reads as lockless-with-no-marker on that
    one tick — so a 40-minute download that just SUCCEEDED would announce itself as stopped."""
    import time as _t
    dlg = _dialog()
    dlg.engine.submit = lambda op, body=None, key=None: None
    comps = [{"key": "dfam_unc_euk", "name": "Dfam uncurated · Eukaryota (optional)", "desc": "3.9 GiB",
              "ok": True, "detail": "present", "repairable": True}]
    dlg._render_components({"wsl2": True, "installing": True, "ready": True, "components": comps})
    dlg._poll_since = _t.monotonic() - 60
    # the tick where the lock is already gone but this tick's log snapshot predates the DONE line
    dlg._render_components({"wsl2": True, "installing": False, "ready": True, "components": comps})
    assert dlg._poll.isActive() and dlg._busy
    assert "stopped" not in dlg.statusLine.text()
    # the next tick's log carries the marker, and the normal completion path takes over
    dlg._refresh = lambda: None                         # _on_done ends with a re-probe that repaints the line
    dlg._on_done("log", {"log": "[teagle] STEP dfam_unc_euk OK\n[teagle] DONE 2026-07-31T07:23:04Z\n"})
    assert not dlg._poll.isActive() and not dlg._busy
    assert "finished" in dlg.statusLine.text() and "stopped" not in dlg.statusLine.text()


_GIB = 1073741824
_EUK_UNPACKED = 24242422200      # what the test EXPECTS, stated independently of the value under test
_EUK_ARCHIVE = 4130166278


def _space_need(total, archive, have, unpacked):
    fn = wsl._PRELUDE[wsl._PRELUDE.index("space_need(){"):wsl._PRELUDE.index("# A truncated/corrupt")]
    rc, out = _bash(f"{fn}\nspace_need {total} {archive} {have} {unpacked}\n")
    assert rc == 0
    return int(out)


def test_disk_gate_counts_only_what_is_left_to_fetch_plus_what_it_unpacks_to():
    """The gate must not demand bytes already on disk a second time (a resume), and must not drop the
    download half when the size request fails."""
    u, a = _EUK_UNPACKED, _EUK_ARCHIVE
    assert _space_need(a, a, 0, u) == a + u + _GIB                     # fresh: whole archive still to come
    assert _space_need(a, a, a // 2, u) == a - a // 2 + u + _GIB       # resumed: only the remainder
    assert _space_need(a, a, a, u) == u + _GIB                         # archive already complete
    assert _space_need(a, a, a * 2, u) == u + _GIB                     # never negative
    assert _space_need(0, a, 0, u) == a + u + _GIB                     # unknown live total -> measured size


def test_disk_gate_fails_on_a_full_disk_and_abstains_only_when_df_cannot_answer():
    """`avail_b` of 0 is a FULL disk, not an unknown one; only an empty reading means df could not say."""
    s = _step("dfam_unc_euk")
    assert '[ -n "$avail_b" ] && [ "$avail_b" -lt "$need_b" ]' in s
    assert '"$avail_b" -eq 0 ] ||' not in s          # the earlier form waved a full disk through


def test_measured_sizes_cover_every_partition_and_drive_the_panel_text():
    assert set(wsl._DFAM_UNPACKED_B) == set(wsl._DFAM_FILES)
    assert set(wsl._DFAM_ARCHIVE_B) == set(wsl._DFAM_FILES)
    assert wsl._DFAM_UNPACKED_B["dfam_unc_euk"] == _EUK_UNPACKED
    assert wsl._DFAM_ARCHIVE_B["dfam_unc_euk"] == _EUK_ARCHIVE
    # the panel figure is computed from the same measurement, so text and gate cannot disagree
    desc = [c[3] for c in wsl._COMP_META if c[0] == "dfam_unc_euk"][0]
    assert wsl._gib(_EUK_UNPACKED) in desc and "14 GB" not in desc


def test_curated_coverage_is_stated_per_lineage_and_never_as_one_number():
    """The app told users "the curated library holds just 9 families for Drosophila melanogaster, and
    copia, gypsy, hobo and mdg1 are not among them". Measured 2026-07-31 against famdb: curated-only
    D. melanogaster holds 399 models and DOES contain Copia_I, Copia_LTR, Gypsy_I, Gypsy_LTR, hobo and
    MDG1_I/LTR. The 9 was the yeast figure on the wrong organism."""
    cov = wsl._DFAM_CURATED_COVERAGE
    assert cov["Drosophila melanogaster"] == (399, 998)
    assert cov["Saccharomyces cerevisiae"] == (9, 398)
    assert cov["Homo sapiens"][0] == cov["Homo sapiens"][1]      # curated alone is the whole set for human
    for sp, (c, u) in cov.items():
        assert 0 < c <= u, sp
    sent = wsl.curated_coverage_sentence()
    for n in (cov["Homo sapiens"][0], cov["Arabidopsis thaliana"][0], cov["Arabidopsis thaliana"][1]):
        assert str(n) in sent                                    # built from the table, not retyped
    import pathlib
    src = pathlib.Path(wsl.__file__).parent.parent.parent / "native" / "main.py"
    text = src.read_text(encoding="utf-8")
    assert "just 9 families" not in text                          # the false claim is gone from the UI
    assert "copia, gypsy, hobo and mdg1 are not among them" not in text


_MULTI = [
    {"key": "viennarna", "name": "ViennaRNA (primer QC)", "desc": "optional",
     "ok": False, "detail": "not installed (optional)", "repairable": True},
    {"key": "dfam_unc_root", "name": "Dfam uncurated · root (optional)", "desc": "0.3 MiB",
     "ok": False, "detail": "not installed (optional)", "repairable": True},
    {"key": "dfam_unc_euk", "name": "Dfam uncurated · Eukaryota (optional)", "desc": "3.9 GiB",
     "ok": False, "detail": "not installed (optional)", "repairable": True},
]


def _glyphs(dlg):
    return {k: r["icon"].text() for k, r in dlg._rows.items()}


def test_only_the_running_component_is_marked_working_even_in_a_window_that_did_not_start_it():
    """The install outlives the dialog, and a reopened window cannot remember what a previous one
    clicked — the running step is read back from the script's own markers instead. A single-component
    fixture cannot catch this, so this one carries three."""
    dlg = _dialog()
    dlg.engine.submit = lambda op, body=None, key=None: None
    running = {"wsl2": True, "installing": True, "ready": False, "components": _MULTI}
    dlg._render_components(running)                       # fresh window attaching to a run already going
    assert dlg._active_key is None
    # its first log tick names the step actually running
    dlg._on_done("log", {"log": "[teagle] START\n[teagle] STEP dfam_unc_euk START\n"})
    assert dlg._active_key == "dfam_unc_euk"
    dlg._render_components(running)
    g = _glyphs(dlg)
    assert g["dfam_unc_euk"] != g["viennarna"], g          # only the running one carries the working mark
    assert g["viennarna"] == g["dfam_unc_root"], g
    # when that step ends, nothing is left marked as running
    dlg._refresh = lambda: None
    dlg._on_done("log", {"log": "[teagle] STEP dfam_unc_euk START\n[teagle] STEP dfam_unc_euk OK\n"
                                "[teagle] DONE 2026-07-31T07:23:04Z\n"})
    assert dlg._active_key is None


@pytest.mark.parametrize("log,expect", [
    ("[teagle] STEP micromamba START", "micromamba"),
    ("[teagle] STEP micromamba START\n[teagle] STEP micromamba OK", None),
    ("[teagle] STEP micromamba OK\n[teagle] STEP dfam_root START", "dfam_root"),
    ("[teagle] STEP a START\n[teagle] STEP a OK\n[teagle] STEP b START\n[teagle] STEP b OK", None),
    ("", None),
])
def test_running_step_reads_the_scripts_own_markers(log, expect):
    from install_dialog import InstallDialog
    assert InstallDialog._running_step(log) == expect


def test_integrity_check_finishing_mid_install_does_not_hand_the_buttons_back():
    dlg = _dialog()
    dlg.engine.submit = lambda op, body=None, key=None: None
    dlg._render_components({"wsl2": True, "installing": True, "ready": False, "components": _MULTI})
    assert dlg._busy and dlg._poll.isActive()
    dlg._on_done("wsl_integrity", {"ok": True, "checks": []})
    assert dlg._busy, "an install is still running; its buttons must stay disabled"
    dlg._poll.stop()


def test_starting_a_component_says_which_one_and_where_progress_appears():
    dlg = _dialog()
    dlg._render_components({"wsl2": True, "installing": False, "ready": False, "components": [
        {"key": "dfam_unc_euk", "name": "Dfam uncurated · Eukaryota (optional)", "desc": "3.9 GiB",
         "ok": False, "detail": "not installed (optional)", "repairable": True}]})
    sent = []
    dlg.engine.submit = lambda op, body=None, key=None: sent.append((op, body))
    dlg._repair("dfam_unc_euk")
    assert sent == [("wsl_repair", {"component": "dfam_unc_euk"})]
    line = dlg.statusLine.text()
    assert "installing" in line and "Eukaryota" in line and "log" in line


# ---------- the header exposes it at the top level ----------
def test_main_window_header_has_a_backend_button():
    """It used to be reachable only from inside panel 03, where a user who had not opened that card
    never saw it — yet it is what turns on family naming, splice detection and genome scans."""
    pytest.importorskip("PySide6")
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    native = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "app", "native")
    if native not in sys.path:
        sys.path.insert(0, native)
    from PySide6.QtWidgets import QApplication
    QApplication.instance() or QApplication([])
    import main as m
    opened = []
    orig = m.MainWindow._open_installer                          # patched on the CLASS: the buttons bind at build time
    m.MainWindow._open_installer = lambda self: opened.append("installer")
    try:
        w = m.MainWindow()
        assert w.backendBtn.text().upper() == "BACKEND"
        assert w.backendBtn.parent() is not None                 # lives in the header, not inside a card
        w.backendBtn.click()
        w.wslInstallBtn.click()
    finally:
        m.MainWindow._open_installer = orig
    # both entry points route to the same handler, so they can never drift apart
    assert opened == ["installer", "installer"]
