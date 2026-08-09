"""Visible-progress and honest-label contract for the backend installer.

A user who starts the optional 3.9 GiB Dfam partition must be able to watch it move. Before this,
curl drew its own carriage-return meter: the whole 40-60 minute transfer was ONE line with no newline
in it, so `tail -n` handed the log panel a single 250 kB line and the download looked frozen. These
lock down the three parts that fix: the step reports whole lines, the reader normalises any CR meter
it still meets, and a component that is not installed says Install rather than Repair."""
import os
import subprocess
import sys
import threading
from contextlib import contextmanager
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

import pytest
from teagle_core import wsl


def _step(key="dfam_unc_euk"):
    return wsl._STEP[key]


def _dfam_python_downloader() -> str:
    """Return the exact heredoc executed in a python3-only fresh WSL."""
    marker = 'python3 - "$dest" "$url" <<\'PY\'\n'
    _, remainder = _step().split(marker, 1)
    return remainder.split("\nPY\n", 1)[0]


@contextmanager
def _download_server(payload: bytes, behavior: str):
    """Hermetic HTTP server for the stdlib downloader's wire-level contract."""
    state = {"requests": 0, "ranges": []}

    class Handler(BaseHTTPRequestHandler):
        protocol_version = "HTTP/1.1"

        def log_message(self, *_args):
            pass

        def do_GET(self):
            state["requests"] += 1
            range_header = self.headers.get("Range")
            state["ranges"].append(range_header)
            if behavior == "transient" and state["requests"] == 1:
                self.send_response(503)
                self.send_header("Content-Length", "0")
                self.end_headers()
                return
            start = int(range_header.split("=", 1)[1].split("-", 1)[0]) if range_header else 0
            body = payload[start:] if range_header and behavior != "ignore_range" else payload
            status = 206 if range_header and behavior != "ignore_range" else 200
            self.send_response(status)
            if status == 206:
                content_range = ("bytes malformed" if behavior == "malformed_range"
                                 else f"bytes {start}-{len(payload) - 1}/{len(payload)}")
                self.send_header("Content-Range", content_range)
            advertised = len(body) if behavior != "failed_transfer" else len(body) + 5
            self.send_header("Content-Length", str(advertised))
            self.end_headers()
            self.wfile.write(body if behavior != "failed_transfer" else body[:2])
            self.wfile.flush()
            if behavior == "failed_transfer":
                self.close_connection = True

    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    worker = threading.Thread(target=server.serve_forever, daemon=True)
    worker.start()
    try:
        yield f"http://127.0.0.1:{server.server_port}/dfam.gz", state
    finally:
        server.shutdown()
        server.server_close()
        worker.join(timeout=5)


def _run_embedded_downloader(tmp_path, url: str, initial: bytes = b""):
    destination = tmp_path / "dfam.gz"
    if initial:
        destination.write_bytes(initial)
    return subprocess.run([sys.executable, "-", str(destination), url], input=_dfam_python_downloader(),
                          text=True, capture_output=True, timeout=30), destination


def test_python_only_dfam_downloader_fetches_a_fresh_archive(tmp_path):
    payload = b"TEagle-Dfam" * 100
    with _download_server(payload, "normal") as (url, state):
        result, destination = _run_embedded_downloader(tmp_path, url)
    assert result.returncode == 0, result.stderr
    assert destination.read_bytes() == payload
    assert state["ranges"] == [None]


def test_python_only_dfam_downloader_resumes_an_honored_range(tmp_path):
    payload, partial = b"TEagle-Dfam" * 100, b"TEagle-Dfam" * 17
    with _download_server(payload, "normal") as (url, state):
        result, destination = _run_embedded_downloader(tmp_path, url, partial)
    assert result.returncode == 0, result.stderr
    assert destination.read_bytes() == payload
    assert state["ranges"] == [f"bytes={len(partial)}-"]


def test_python_only_dfam_downloader_restarts_safely_when_range_is_ignored(tmp_path):
    payload, partial = b"TEagle-Dfam" * 100, b"obsolete-partial"
    with _download_server(payload, "ignore_range") as (url, state):
        result, destination = _run_embedded_downloader(tmp_path, url, partial)
    assert result.returncode == 0, result.stderr
    assert destination.read_bytes() == payload
    assert state["ranges"] == [f"bytes={len(partial)}-"]
    assert not destination.with_name(destination.name + ".restart").exists()


def test_python_only_dfam_downloader_retries_a_transient_failure(tmp_path):
    payload = b"TEagle-Dfam" * 100
    with _download_server(payload, "transient") as (url, state):
        result, destination = _run_embedded_downloader(tmp_path, url)
    assert result.returncode == 0, result.stderr
    assert destination.read_bytes() == payload
    assert state["requests"] == 2


@pytest.mark.parametrize("behavior", ["malformed_range", "failed_transfer"])
def test_python_only_dfam_downloader_keeps_a_partial_on_malformed_or_failed_transfer(tmp_path, behavior):
    payload, partial = b"TEagle-Dfam" * 100, b"TEagle-Dfam" * 17
    with _download_server(payload, behavior) as (url, _state):
        result, destination = _run_embedded_downloader(tmp_path, url, partial)
    assert result.returncode != 0
    # A malformed Content-Range is rejected before append. A real mid-stream failure can retain the
    # bytes received before the socket closed, but it must never discard the resumable prefix.
    if behavior == "malformed_range":
        assert destination.read_bytes() == partial
    else:
        assert destination.read_bytes().startswith(partial)
    assert not destination.with_name(destination.name + ".restart").exists()


# ---------- the download step reports progress as whole lines ----------
def test_dfam_step_silences_curl_meter_and_reports_lines():
    s = _step()
    assert "curl -sS -L --fail -C -" in s          # -s kills the CR meter, -S keeps real errors
    assert "elif command -v python3" in s and "urllib.request" in s
    assert 'headers = {"Range": f"bytes={have}-"}' in s  # stdlib path resumes an existing partial
    assert "sleep 15" in s and "kill -0 \"$dl_pid\"" in s   # a watcher polls the file while curl runs
    assert 'echo "[teagle] ' in s.split("dl_pid=$!", 1)[1]  # and reports on its own newline-terminated lines


def test_range_ignored_restart_progress_watches_the_active_restart_file():
    watcher = _step().split('while kill -0 "$dl_pid"', 1)[1].split('wait "$dl_pid"', 1)[0]
    assert '[ -f "dfam40.uncurated.consensus.1.h5.gz.restart" ]' in watcher
    assert 'active="dfam40.uncurated.consensus.1.h5.gz.restart"' in watcher
    assert 'stat -c %s "$active"' in watcher


def test_dfam_step_percentage_is_derived_from_a_head_not_the_resumed_get():
    """A resumed (-C -) GET answers 206 with Content-Length = REMAINING bytes. Reading the total off
    the transfer would make every percentage wrong, so the total comes from a separate HEAD."""
    s = _step()
    assert "dfam_head_bytes" in s and "content-length" in s.lower()
    assert "curl -sIL" in s and 'Request(sys.argv[1], method="HEAD")' in s
    head, rest = s.split("total_b=$(dfam_head_bytes", 1)
    assert "-C -" not in rest.split(")", 1)[0]     # the HEAD itself must not be a resumed request
    # a percentage is only printed when the total is known AND consistent with the bytes on disk
    assert '[ "$total_b" -gt 0 ] && [ "$now_b" -le "$total_b" ]' in s


def _bash(script: str):
    """Run a bash snippet and return (rc, stdout). Skips where no bash is available."""
    import shutil, subprocess
    bash = shutil.which("bash")
    if not bash:
        pytest.skip("no bash on PATH")
    # A native temporary path can be mangled when Windows launches bash/WSL; deliver script bytes directly.
    r = subprocess.run([bash, "-s"], input=script.replace("\r\n", "\n").encode(),
                       capture_output=True, timeout=60)
    return r.returncode, r.stdout.decode(errors="replace").strip()


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


def test_integrity_probe_conditions_success_on_real_exit_codes_and_signatures():
    probe = wsl._INTEGRITY_PROBE
    assert 'ispcr_out=$("$ENV/bin/isPcr" 2>&1); ispcr_rc=$?' in probe
    assert '"isPcr - Standalone v 33x2"*' in probe and '[ "$ispcr_rc" -eq 255 ]' in probe
    assert 'datasets_rc=$?' in probe and 'datasets version 18.34.0 OK' in probe
    assert '"faToTwoBit - Convert DNA"*' in probe and '[ "$fa_rc" -eq 255 ]' in probe
    assert "isPcr executable ran" not in probe and "datasets executable ran" not in probe


def test_integrity_check_reports_optional_partitions_without_failing_on_them(monkeypatch):
    probe = ("=RM=\nRepeatMasker version 4.2.4\n=MM=\n2.31-r1302\n=FAMDB=\nVersion : 4.0\n"
             "=FILES=\npresent dfam40.0.h5 1\npresent dfam40.curated.consensus.0.h5 1\n"
             "=OPTFILES=\npresent dfam40.uncurated.consensus.0.h5\nabsent dfam40.uncurated.consensus.1.h5\n"
             "=SCAN=\nisPcr - Standalone v 33x2\nisPcr signature 33x2 OK\n"
             "datasets version: 18.34.0\ndatasets version 18.34.0 OK\n"
             "ispcr package 33\ndatasets package 18.34.0\n"
             "=FATOTWOBIT=\nabsent\n=VRNA=\nabsent\n")
    monkeypatch.setattr(wsl, "available", lambda: {"wsl2": True, "distro": "Ubuntu", "error": None})
    monkeypatch.setattr(wsl, "_wsl_script", lambda *a, **k: (0, probe, ""))
    r = wsl.integrity_check()
    files = [c for c in r["checks"] if c["name"] == "Dfam library files present"][0]
    assert r["ok"] is True and files["ok"] is True              # one optional partition missing is not a fault
    assert "uncurated.consensus.0.h5 present" in files["detail"]
    assert "uncurated.consensus.1.h5 not installed" in files["detail"]


@pytest.mark.parametrize("section, body, expected", [
    ("FATOTWOBIT", "faToTwoBit - Convert DNA\nfaToTwoBit signature OK\npresent 482", True),
    ("FATOTWOBIT", "faToTwoBit - Convert DNA\nfaToTwoBit signature OK\npresent 481", False),
    ("VRNA", "present 2.7.2", True),
    ("VRNA", "present 2.7.1", False),
])
def test_integrity_grades_installed_optional_tools_at_their_pinned_version(monkeypatch, section, body, expected):
    probe = ("=RM=\nRepeatMasker version 4.2.4\n=MM=\n2.31-r1302\n=FAMDB=\nVersion : 4.0\n"
             "=FILES=\npresent dfam40.0.h5 1\npresent dfam40.curated.consensus.0.h5 1\n"
             "=OPTFILES=\nabsent dfam40.uncurated.consensus.0.h5\nabsent dfam40.uncurated.consensus.1.h5\n"
             "=SCAN=\nisPcr - Standalone v 33x2\nisPcr signature 33x2 OK\n"
             "datasets version: 18.34.0\ndatasets version 18.34.0 OK\n"
             "ispcr package 33\ndatasets package 18.34.0\n"
             f"={section}=\n{body}\n")
    other = "VRNA" if section == "FATOTWOBIT" else "FATOTWOBIT"
    probe += f"={other}=\nabsent\n"
    monkeypatch.setattr(wsl, "available", lambda: {"wsl2": True, "distro": "Ubuntu", "error": None})
    monkeypatch.setattr(wsl, "_wsl_script", lambda *a, **k: (0, probe, ""))
    result = wsl.integrity_check()
    optional = [check for check in result["checks"] if section.lower().replace("vrna", "viennarna") in check["name"].lower()][0]
    assert optional["ok"] is expected
    assert result["ok"] is expected


@pytest.mark.parametrize("section, bad_body, check_name", [
    ("SCAN", "/env/bin/isPcr: cannot execute: required file not found\nisPcr BROKEN rc=127 wrong signature\n"
             "datasets version: 18.34.0\ndatasets version 18.34.0 OK\n"
             "ispcr package 33\ndatasets package 18.34.0", "isPcr + NCBI datasets"),
    ("SCAN", "unexpected isPcr banner\nisPcr BROKEN rc=1 wrong signature\n"
             "datasets version: 18.34.0\ndatasets version 18.34.0 OK\n"
             "ispcr package 33\ndatasets package 18.34.0", "isPcr + NCBI datasets"),
    ("SCAN", "isPcr - Standalone v 33x2\nisPcr signature 33x2 OK\n"
             "datasets version: 18.34.1\ndatasets BROKEN rc=0 wrong version/signature\n"
             "ispcr package 33\ndatasets package 18.34.0", "isPcr + NCBI datasets"),
    ("FATOTWOBIT", "/env/bin/faToTwoBit: cannot execute: required file not found\n"
                     "faToTwoBit BROKEN rc=127 wrong signature\npresent 482", "faToTwoBit"),
    ("FATOTWOBIT", "unexpected faToTwoBit banner\nfaToTwoBit BROKEN rc=1 wrong signature\npresent 482", "faToTwoBit"),
])
def test_integrity_rejects_unusable_or_wrong_signature_executables(monkeypatch, section, bad_body, check_name):
    healthy_scan = ("isPcr - Standalone v 33x2\nisPcr signature 33x2 OK\n"
                    "datasets version: 18.34.0\ndatasets version 18.34.0 OK\n"
                    "ispcr package 33\ndatasets package 18.34.0")
    probe = ("=RM=\nRepeatMasker version 4.2.4\n=MM=\n2.31-r1302\n=FAMDB=\nVersion : 4.0\n"
             "=FILES=\npresent dfam40.0.h5 1\npresent dfam40.curated.consensus.0.h5 1\n"
             "=OPTFILES=\nabsent dfam40.uncurated.consensus.0.h5\nabsent dfam40.uncurated.consensus.1.h5\n"
             f"=SCAN=\n{bad_body if section == 'SCAN' else healthy_scan}\n"
             f"=FATOTWOBIT=\n{bad_body if section == 'FATOTWOBIT' else 'absent'}\n=VRNA=\nabsent\n")
    monkeypatch.setattr(wsl, "available", lambda: {"wsl2": True, "distro": "Ubuntu", "error": None})
    monkeypatch.setattr(wsl, "_wsl_script", lambda *a, **k: (0, probe, ""))
    result = wsl.integrity_check()
    failed = next(check for check in result["checks"] if check_name in check["name"])
    assert failed["ok"] is False
    assert result["ok"] is False
    assert "cannot execute" in failed["detail"] or "BROKEN" in failed["detail"]


@pytest.mark.parametrize("fa_state, expected_detail", [
    ("absent", "not installed (optional"),
    ("version drift", "version drift"),
    ("broken (pinned metadata, executable missing)", "broken"),
])
def test_component_status_reports_optional_fatotwobit_without_breaking_required_genome_scan(monkeypatch, fa_state, expected_detail):
    status = ("micromamba=2.8.1\nrepeatmasker=4.2.4\nminimap2=2.31-r1302\ngenomescan=1\n"
              "ispcr_pin=1\ndatasets_pin=1\n"
              f"fatotwobit={fa_state}\n"
              "viennarna=absent\ndfam_root=1\ndfam_curated=1\ndfam_unc_root=0\ndfam_unc_euk=0\n"
              "famdb_conf=1\ninstalling=0\ndisk_free_gb=100\n")
    monkeypatch.setattr(wsl, "available", lambda: {"wsl2": True, "distro": "Ubuntu", "error": None})
    monkeypatch.setattr(wsl, "_wsl_script", lambda *a, **k: (0, status, ""))
    components = {c["key"]: c for c in wsl.components_status()["components"]}
    assert components["genomescan"]["ok"] is True
    assert expected_detail in components["genomescan"]["detail"]


def test_status_probe_distinguishes_absent_from_present_broken_viennarna():
    probe = wsl._STATUS_PROBE
    assert '[ ! -d "$VRNAENV" ]' in probe and 'viennarna=absent' in probe
    assert 'viennarna=broken/unusable' in probe
    assert 'import RNA;print(RNA.__version__)' in probe


@pytest.mark.parametrize("vrna_state, expected_ok, expected_detail", [
    ("absent", False, "not installed (optional)"),
    ("broken/unusable", False, "broken/unusable"),
    ("2.7.2", True, "v2.7.2"),
])
def test_component_status_reports_viennarna_absent_broken_and_exact(monkeypatch, vrna_state, expected_ok, expected_detail):
    status = ("micromamba=2.8.1\nrepeatmasker=4.2.4\nminimap2=2.31-r1302\ngenomescan=1\n"
              "ispcr_pin=1\ndatasets_pin=1\nfatotwobit=absent\n"
              f"viennarna={vrna_state}\n"
              "dfam_root=1\ndfam_curated=1\ndfam_unc_root=0\ndfam_unc_euk=0\n"
              "famdb_conf=1\ninstalling=0\ndisk_free_gb=100\n")
    monkeypatch.setattr(wsl, "available", lambda: {"wsl2": True, "distro": "Ubuntu", "error": None})
    monkeypatch.setattr(wsl, "_wsl_script", lambda *a, **k: (0, status, ""))
    component = next(c for c in wsl.components_status()["components"] if c["key"] == "viennarna")
    assert component["ok"] is expected_ok
    assert expected_detail in component["detail"]


def test_every_conda_install_step_creates_the_environment_through_the_recovery():
    """mm_create carries the non-conda-prefix recovery; a step that called mm_install without it would
    still abort with no way back."""
    for key in ("repeatmasker", "minimap2", "genomescan"):
        b = wsl._STEP[key]
        assert "mm_install" in b and "mm_create" in b, key
        assert b.index("mm_create") < b.index("mm_install"), key


def test_failed_package_probe_forces_file_restoration_even_when_metadata_survives():
    """A missing binary must not be accepted as healthy because conda still records the package."""
    helper = wsl._PRELUDE[wsl._PRELUDE.index("mm_install(){"):]
    first_install = helper.split("mm_reset_cache", 1)[0]
    assert 'install --force-reinstall -y -n te' in first_install
    for key in ("repeatmasker", "minimap2", "miniprot", "genomescan"):
        assert "mm_install" in wsl._STEP[key], key
    assert "install --force-reinstall -y -n teagle-vrna" in wsl._STEP["viennarna"]
    assert "install --force-reinstall -y -n te" in wsl._STEP["genomescan"]  # optional faToTwoBit


def test_backend_package_pins_cover_every_installed_tool():
    assert wsl._CONDA_PINS == {
        "repeatmasker": "4.2.4", "minimap2": "2.31", "miniprot": "0.18",
        "ispcr": "33", "ncbi-datasets-cli": "18.34.0", "ucsc-fatotwobit": "482",
        "viennarna": "2.7.2",
    }
    assert "repeatmasker=4.2.4" in wsl._STEP["repeatmasker"]
    assert "minimap2=2.31" in wsl._STEP["minimap2"]
    assert "miniprot=0.18" in wsl._STEP["miniprot"]
    assert "ispcr=33" in wsl._STEP["genomescan"]
    assert "ncbi-datasets-cli=18.34.0" in wsl._STEP["genomescan"]
    assert "ucsc-fatotwobit=482" in wsl._STEP["genomescan"]
    assert "viennarna=2.7.2" in wsl._STEP["viennarna"]


def test_micromamba_bootstrap_is_versioned_and_checksum_gated():
    s = wsl._STEP["micromamba"]
    assert wsl._MICROMAMBA_VERSION == "2.8.1"
    assert "/2.8.1-0/micromamba-linux-64.tar.bz2" in s
    assert wsl._MICROMAMBA_SHA256 in s
    assert "hashlib.sha256(data).hexdigest()" in s
    assert "sha256sum -c -" in s
    assert "/latest" not in s


def test_prelude_kills_background_helpers_on_exit():
    """Closing the app kills the WSL session and this script; an orphaned watcher would keep writing
    progress for a download that is no longer running."""
    assert "cleanup(){" in wsl._PRELUDE and "kill $BG" in wsl._PRELUDE
    assert "trap cleanup EXIT" in wsl._PRELUDE


@pytest.mark.parametrize("key", sorted(wsl._STEP))
def test_every_install_step_is_valid_bash(key):
    """`bash -n` every real generated script before any package or database is downloaded."""
    import shutil, subprocess
    bash = shutil.which("bash")
    if not bash:
        pytest.skip("no bash on PATH")
    script = wsl._build_script([key]).replace("\r\n", "\n").encode()
    r = subprocess.run([bash, "-n", "-s"], input=script, capture_output=True, timeout=60)
    assert r.returncode == 0, r.stderr.decode(errors="replace")


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
    yeast_counts = wsl.dfam_coverage_counts("Saccharomyces cerevisiae")
    assert yeast_counts == (9, 398, 389)
    assert cov["Homo sapiens"][0] == cov["Homo sapiens"][1]      # curated alone is the whole set for human
    for sp, (c, u) in cov.items():
        assert 0 < c <= u, sp
    sent = wsl.curated_coverage_sentence()
    for n in (cov["Homo sapiens"][0], cov["Arabidopsis thaliana"][0], cov["Arabidopsis thaliana"][1]):
        assert str(n) in sent                                    # built from the table, not retyped
    import pathlib
    src = pathlib.Path(wsl.__file__).parent.parent.parent / "native" / "main.py"
    text = src.read_text(encoding="utf-8")
    backend_text = pathlib.Path(wsl.__file__).read_text(encoding="utf-8")
    assert "just 9 families" not in text                          # the false claim is gone from the UI
    assert "copia, gypsy, hobo and mdg1 are not among them" not in text
    assert "421 more" not in text + backend_text
    assert text.count("{_DFAM_YEAST_ADDITIONAL}") >= 1
    assert backend_text.count("{_DFAM_YEAST_ADDITIONAL}") >= 2
    desc = next(c[3] for c in wsl._COMP_META if c[0] == "dfam_unc_euk")
    assert all(str(n) in desc for n in yeast_counts)


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
