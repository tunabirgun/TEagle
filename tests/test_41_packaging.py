"""Packaging, launcher and documentation contracts for 4.1.0 (lane A7). Behavioural and hermetic."""
from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_READY = {"python": "3.12.0", "python_ok": True, "app_version": "x", "packages_ok": True, "needs_install": False,
          "first_run": False, "installed_now": False, "error": None, "backends": {}}


def _launcher(monkeypatch, argv, calls):
    mod = _load("teagle_launcher_under_test", ROOT / "app" / "teagle.py")
    monkeypatch.setattr(mod.envcheck, "ensure", lambda auto_install=True: calls.append(auto_install) or dict(_READY))
    monkeypatch.setattr(sys, "argv", ["teagle.py", *argv])
    return mod


def test_check_flag_never_installs(monkeypatch, capsys):
    calls = []
    assert _launcher(monkeypatch, ["--check"], calls).main() == 0
    assert calls == [False]


def test_source_launch_still_installs_by_default(monkeypatch):
    calls = []
    mod = _launcher(monkeypatch, [], calls)
    fake_main = type(sys)("main")
    fake_main.main = lambda: 0
    monkeypatch.setitem(sys.modules, "main", fake_main)
    assert mod.main() == 0
    assert calls == [True]


def test_server_flag_is_retired(monkeypatch):
    calls, served = [], []
    fake_server = type(sys)("server")
    fake_server.serve = lambda **k: served.append(k)
    monkeypatch.setitem(sys.modules, "server", fake_server)
    mod = _launcher(monkeypatch, ["--server"], calls)
    with pytest.raises(SystemExit) as exc:
        mod.main()
    assert exc.value.code == 2 and not served and not calls


def test_no_module_or_test_imports_the_retired_web_server():
    tracked = subprocess.run(["git", "ls-files", "*.py"], cwd=ROOT, capture_output=True, text=True, check=True).stdout.split()
    assert tracked
    for rel in tracked:
        path = ROOT / rel
        if rel == "tests/test_41_packaging.py" or not path.exists():
            continue
        for ln in path.read_text(encoding="utf-8").splitlines():
            assert not ln.strip().startswith(("import server", "from server ")), (rel, ln)
    assert not (ROOT / "app" / "web").exists() and not (ROOT / "app" / "backend" / "server.py").exists()


import ast
import os

SPEC = ROOT / "installer" / "teagle_native.spec"


def _spec_nodes():
    return ast.parse(SPEC.read_text(encoding="utf-8"), filename=str(SPEC)).body


def _spec_assigned(name):
    for node in _spec_nodes():
        if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == name for t in node.targets):
            return ast.literal_eval(node.value)
    raise AssertionError(name)


def test_spec_bundles_no_web_server_module():
    assert "server" not in _spec_assigned("hiddenimports")
    assert "server" not in _spec_assigned("_ALLOWED_TOP")


def test_spec_writes_version_resource_to_the_build_workpath(tmp_path):
    head = []
    for node in _spec_nodes():
        if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == "datas" for t in node.targets):
            break
        if not isinstance(node, (ast.Import, ast.ImportFrom)):
            head.append(node)
    tracked = ROOT / "installer" / "_version_info.txt"
    before = tracked.read_bytes() if tracked.exists() else None
    ns = {"os": os, "re": __import__("re"), "SPECPATH": str(ROOT / "installer"), "workpath": str(tmp_path / "work")}
    exec(compile(ast.Module(body=head, type_ignores=[]), str(SPEC), "exec"), ns)
    written = tmp_path / "work" / "_version_info.txt"
    assert written.is_file() and ns["VERSION"] in written.read_text(encoding="utf-8")
    assert ns["_verinfo"] == str(written)
    assert (tracked.read_bytes() if tracked.exists() else None) == before


def test_spec_fails_the_build_when_package_metadata_cannot_be_collected():
    loop = next(n for n in _spec_nodes() if isinstance(n, ast.For) and getattr(n.target, "id", "") == "dist")

    def missing(dist):
        raise RuntimeError(f"No package metadata was found for {dist}")

    ns = {"datas": [], "copy_metadata": missing}
    with pytest.raises(RuntimeError):
        exec(compile(ast.Module(body=[loop], type_ignores=[]), str(SPEC), "exec"), ns)


def _native_launcher(monkeypatch, fake_main):
    mod = _load("teagle_native_under_test", ROOT / "installer" / "teagle_native.py")
    monkeypatch.setitem(sys.modules, "main", fake_main)
    monkeypatch.setattr(mod, "_enable_kill_on_close_job", lambda: None)
    monkeypatch.setattr(mod, "_windows_single_instance", lambda: False)
    boxes = []
    monkeypatch.setattr(mod, "_msgbox", lambda *a, **k: boxes.append(a))
    return mod, boxes


def _close_kept(mod):
    for h in mod._KEEP:
        if hasattr(h, "close"):
            h.close()


def test_windowed_launch_without_std_handles_logs_instead_of_crashing(monkeypatch, tmp_path):
    fake = type(sys)("main")
    fake.main = lambda: sys.stderr.write("engine fault trace\n") and 0
    mod, boxes = _native_launcher(monkeypatch, fake)
    monkeypatch.setenv("LOCALAPPDATA", str(tmp_path))
    monkeypatch.setattr(sys, "stdout", None)
    monkeypatch.setattr(sys, "stderr", None)
    try:
        assert mod.main() == 0
        sys.stderr.flush()
    finally:
        _close_kept(mod)
    assert not boxes
    assert "engine fault trace" in (tmp_path / "TEagle" / "logs" / "teagle.log").read_text(encoding="utf-8")


def test_unwritable_log_location_falls_back_to_devnull(monkeypatch, tmp_path):
    blocker = tmp_path / "blocked"
    blocker.write_text("a file where the LOCALAPPDATA directory should be")
    fake = type(sys)("main")
    fake.main = lambda: sys.stderr.write("trace\n") and 0
    mod, boxes = _native_launcher(monkeypatch, fake)
    monkeypatch.setenv("LOCALAPPDATA", str(blocker))
    monkeypatch.setattr(sys, "stdout", None)
    monkeypatch.setattr(sys, "stderr", None)
    try:
        assert mod.main() == 0
    finally:
        _close_kept(mod)
    assert not boxes


def test_existing_std_streams_are_left_alone(monkeypatch, tmp_path, capsys):
    fake = type(sys)("main")
    fake.main = lambda: print("visible") or 0
    mod, _ = _native_launcher(monkeypatch, fake)
    monkeypatch.setenv("LOCALAPPDATA", str(tmp_path))
    assert mod.main() == 0
    assert "visible" in capsys.readouterr().out and not (tmp_path / "TEagle").exists()


def test_second_launch_raises_the_running_window_or_says_so(monkeypatch):
    fake = type(sys)("main")
    fake.main = lambda: pytest.fail("a second instance must not start")
    mod, boxes = _native_launcher(monkeypatch, fake)
    monkeypatch.setattr(mod, "_windows_single_instance", lambda: True)
    monkeypatch.setattr(sys, "platform", "win32")
    raised = []
    monkeypatch.setattr(mod, "_raise_running_window", lambda *a: raised.append(1) or True)
    assert mod.main() == 0 and raised and not boxes
    monkeypatch.setattr(mod, "_raise_running_window", lambda *a: False)
    assert mod.main() == 0 and len(boxes) == 1 and "already running" in boxes[0][0]


@pytest.mark.skipif(sys.platform != "win32", reason="Win32 window lookup")
def test_raise_running_window_reports_absence():
    mod = _load("teagle_native_window_probe", ROOT / "installer" / "teagle_native.py")
    assert mod._raise_running_window("TEagle no such window 7f3a") is False


def test_job_object_permits_breakaway_so_links_can_outlive_the_app():
    mod = _load("teagle_native_job_flags", ROOT / "installer" / "teagle_native.py")
    assert mod._JOB_LIMIT_FLAGS & 0x2000 and mod._JOB_LIMIT_FLAGS & 0x800


@pytest.mark.skipif(sys.platform != "win32", reason="Win32 job objects")
def test_real_job_object_lets_a_child_break_away():
    code = (
        "import importlib.util,subprocess,sys\n"
        "s=importlib.util.spec_from_file_location('tn',sys.argv[1]);m=importlib.util.module_from_spec(s);s.loader.exec_module(m)\n"
        "m._enable_kill_on_close_job()\n"
        "subprocess.run([sys.executable,'-c','pass'],creationflags=0x01000000,check=True)\n"   # CREATE_BREAKAWAY_FROM_JOB
    )
    r = subprocess.run([sys.executable, "-c", code, str(ROOT / "installer" / "teagle_native.py")], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr


def _ignored(rel):
    r = subprocess.run(["git", "check-ignore", "--no-index", "-q", rel], cwd=ROOT)
    if r.returncode not in (0, 1):
        pytest.skip("git check-ignore unavailable")
    return r.returncode == 0


def test_workspace_ignore_rule_does_not_hide_ci_workflows():
    assert not _ignored(".github/workflows/new-check.yml")
    assert _ignored("workflows/scratch.json")
    assert _ignored("installer/_version_info.txt")


import re
import shutil

ISS = ROOT / "installer" / "teagle.iss"


def _iss_section(name):
    text = ISS.read_text(encoding="utf-8")
    m = re.search(rf"(?ms)^\[{name}\]\s*$(.*?)(?=^\[|\Z)", text)
    assert m, name
    return [ln for ln in m.group(1).splitlines() if ln.strip() and not ln.lstrip().startswith(";")]


def test_installer_upgrade_replaces_the_whole_previous_bundle():
    assert any(re.search(r'Name:\s*"\{app\}\\_internal"', ln) and "filesandordirs" in ln for ln in _iss_section("InstallDelete"))


def test_installer_ships_the_licence_texts_the_binary_needs():
    shipped = set()
    for ln in _iss_section("Files"):
        src = re.search(r'Source:\s*"([^"]+)"', ln).group(1)
        dest = re.search(r'DestDir:\s*"([^"]+)"', ln).group(1)
        if dest.endswith("licenses"):
            found = list((ROOT / "installer").glob(src.replace("\\", "/"))) if "*" in src else [ROOT / "installer" / src.replace("\\", "/")]
            assert found and all(f.resolve().is_file() for f in found), src
            shipped |= {f.name for f in found}
    assert {"LICENSE", "THIRD-PARTY-NOTICES.md", "LGPL-3.0.txt", "GPL-3.0.txt",
            "HMMER-BSD-3-Clause.txt", "Easel-BSD-2-Clause.txt"} <= shipped
    assert "BSD three-clause" in (ROOT / "installer" / "licenses" / "HMMER-BSD-3-Clause.txt").read_text(encoding="utf-8")
    assert "GNU LESSER GENERAL PUBLIC LICENSE" in (ROOT / "installer" / "licenses" / "LGPL-3.0.txt").read_text(encoding="utf-8", errors="replace")


def test_installer_refuses_windows_older_than_qt_supports():
    ver = re.search(r"(?m)^MinVersion=(\d+)\.(\d+)\.(\d+)", ISS.read_text(encoding="utf-8"))
    assert ver and tuple(map(int, ver.groups())) >= (10, 0, 17763)


def test_clean_install_removes_registry_settings_and_tolerates_a_missing_wsl():
    code = ISS.read_text(encoding="utf-8")
    clean = code[code.index("function CleanTeagleData"):code.index("procedure CurStepChanged")]
    assert r"RegDeleteKeyIncludingSubkeys(HKCU, 'Software\TEagle\TEagle')" in clean
    wslfn = code[code.index("function CleanWslBackend"):code.index("function CleanTeagleData")]
    assert wslfn.index("not FileExists(wsl)") < wslfn.index("Exec(")
    assert "WslDistroMayExist(Distro)" in wslfn and wslfn.index("WslDistroMayExist(Distro)") < wslfn.index("Exec(")
    assert r"Lxss\' + keys[I]" in code


def _iscc():
    import os
    for base in (os.environ.get("LOCALAPPDATA"), os.environ.get("ProgramFiles(x86)"), os.environ.get("ProgramFiles")):
        if base:
            c = Path(base) / ("Programs" if base == os.environ.get("LOCALAPPDATA") else "") / "Inno Setup 6" / "ISCC.exe"
            if c.is_file():
                return c
    return None


@pytest.mark.skipif(_iscc() is None, reason="Inno Setup 6 is not installed")
def test_installer_script_compiles(tmp_path):
    inst = tmp_path / "installer"
    shutil.copytree(ROOT / "installer", inst, ignore=shutil.ignore_patterns("__pycache__"))
    for name in ("LICENSE", "THIRD-PARTY-NOTICES.md"):
        shutil.copy2(ROOT / name, tmp_path / name)
    (tmp_path / "dist" / "TEagle").mkdir(parents=True)
    (tmp_path / "dist" / "TEagle" / "TEagle.exe").write_bytes(b"MZ")
    r = subprocess.run([str(_iscc()), "/Q", "/DMyAppVersion=0.0.1", str(inst / "teagle.iss")], capture_output=True, text=True, timeout=300)
    assert r.returncode == 0, r.stdout + r.stderr


_LOCAL_ONLY = ("report/", "manuscript/", "verification/")
# manuscript-guardrail tests pending the owner's decision on keeping them in the public tree
_PENDING_OWNER = {
    "test_version_locked_benchmark_docs_match_their_derived_evidence_version",
    "test_report_benchmark_spec_is_an_exact_mirror_when_maintained",
    "test_mobile_dna_manuscript_submission_guardrails",
    "test_submission_assets_use_reader_facing_evidence_and_current_status",
}


def _read_calls(func):
    out = []
    for n in ast.walk(func):
        if isinstance(n, ast.Call) and isinstance(n.func, ast.Name) and n.func.id == "_read":
            out += [c.value for c in ast.walk(n) if isinstance(c, ast.Constant) and isinstance(c.value, str)]
        if isinstance(n, ast.For) and isinstance(n.iter, ast.Tuple) and any(
                isinstance(c, ast.Call) and getattr(c.func, "id", "") == "_read" for c in ast.walk(n)):
            out += [e.value for e in n.iter.elts if isinstance(e, ast.Constant) and isinstance(e.value, str)]
    return out


def test_no_doc_test_mixes_tracked_files_with_local_artifacts_under_a_whole_test_skip():
    """`_read` skips the whole test when a report/ or manuscript/ file is absent, so a test that also asserts on
    tracked files would silently stop checking them in CI. Such tests must read local files through `_local`."""
    tree = ast.parse((ROOT / "tests" / "test_docs_track_build.py").read_text(encoding="utf-8"))
    mixed = []
    for f in tree.body:
        if isinstance(f, ast.FunctionDef) and f.name.startswith("test_") and f.name not in _PENDING_OWNER:
            paths = _read_calls(f)
            local = [p for p in paths if p.startswith(_LOCAL_ONLY)]
            tracked = [p for p in paths if "/" in p or p.endswith((".md", ".cff"))]
            tracked = [p for p in tracked if not p.startswith(_LOCAL_ONLY)]
            if local and tracked:
                mixed.append(f.name)
    assert not mixed, mixed


def test_the_mixed_test_guard_detects_the_old_pattern():
    old = ast.parse('def test_x():\n    a = _read("README.md")\n    b = _read("report/teagle_manual.tex")\n').body[0]
    paths = _read_calls(old)
    assert any(p.startswith(_LOCAL_ONLY) for p in paths) and "README.md" in paths


def _txt(rel):
    return (ROOT / rel).read_text(encoding="utf-8")


def _flat(rel):
    return " ".join(_txt(rel).split())


def test_source_launch_docs_state_the_exact_python_the_launcher_enforces():
    import envcheck
    ver = ".".join(map(str, envcheck.REQUIRED_PY))
    for rel in ("app/README.md", "README.md"):
        assert f"Python {ver}" in _txt(rel), rel
    assert "Python 3.10" not in _txt("app/README.md")


def test_unclassified_dde_group_is_named_where_the_panel_is_enumerated():
    from teagle_core import classify
    assert "unclassified DDE" in classify.DOMAINS_TESTED
    assert "unclassified DDE" in _flat("app/README.md")


def test_docs_agree_that_a_4_0_0_evaluation_exists_and_is_scoped_to_that_release():
    doi = "10.5281/zenodo.22159977"
    for rel in ("README.md", "app/README.md", "docs/TRY-IT-YOURSELF.md", "CITATION.cff"):
        assert doi in _txt(rel), rel
    for stale in ("Neither evidence set is presented as a 4.0.0 execution",
                  "no 4.0.0 benchmark or prospective wet-lab result is claimed",
                  "do not establish TEagle 4.0.0 performance"):
        assert stale not in "\n".join(_flat(r) for r in ("README.md", "app/README.md", "docs/TRY-IT-YOURSELF.md", "CITATION.cff"))


def test_fetch_panel_documentation_does_not_promise_a_taxid_the_panel_omits():
    assert "taxid" not in _flat("app/README.md").split("Fetch by accession")[1].split("Using it")[0]


def test_third_party_notices_cover_every_pinned_package_and_profile_accession():
    import re as _re
    notices = _txt("THIRD-PARTY-NOTICES.md")
    low = notices.lower()
    for line in _txt("app/backend/requirements.txt").splitlines():
        m = _re.match(r"\s*([A-Za-z0-9_.-]+)\s*==\s*([^\s#;]+)", line)
        if m:
            name = m.group(1).lower().replace("_", "-")
            assert name.replace("-", "") in low.replace("-", "").replace("_", ""), name
            assert m.group(2) in notices, (name, m.group(2))
    hmm = _txt("app/backend/data/te_domains.hmm")
    accs = _re.findall(r"(?m)^ACC\s+(\S+)", hmm)
    assert len(accs) == len(_re.findall(r"(?m)^HMMER3", hmm)) == 30
    for acc in accs:
        assert acc in notices, acc
    assert "37.x" not in notices
    for needed in ("HMMER", "Easel", "psutil", "PSF-2.0", "OpenSSL", "LGPL-3.0"):
        assert needed in notices, needed


def test_tracked_documents_contain_no_mojibake():
    # Frozen 4.0.0 deposit files: recorded evidence is never edited, so these two are exempt permanently.
    recorded_evidence = {"reviewer_package/README.md", "benchmarks/protocol/BENCHMARK_SPEC_v1.md"}
    files = subprocess.run(["git", "ls-files", "*.md", "*.txt", "*.cff", "*.ps1", "*.iss", "*.yml"], cwd=ROOT,
                           capture_output=True, text=True, check=True).stdout.split()
    marker = "\u00e2\u20ac"
    bad = [f for f in files if (ROOT / f).exists() and marker in (ROOT / f).read_text(encoding="utf-8", errors="replace")]
    assert not [f for f in bad if f not in recorded_evidence], bad


def test_docs_img_index_does_not_cite_untracked_harness_paths():
    text = _txt("docs/img/README.md")
    assert "verification/" not in text
    tracked = set(subprocess.run(["git", "ls-files", "docs/img"], cwd=ROOT, capture_output=True, text=True,
                                 check=True).stdout.split())
    for name in re.findall(r"`((?:teagle-[\w-]+|[\w]+)\.(?:png|svg))`", text):
        assert f"docs/img/{name}" in tracked or f"docs/img/screenshots/{name}" in tracked, name


def test_installer_probe_of_a_missing_pyinstaller_is_not_a_terminating_error_under_ps51():
    src = _txt("installer/build_installer.ps1")
    probe = next(ln for ln in src.splitlines() if "m.version('PyInstaller')" in ln)
    assert "$ErrorActionPreference = 'Continue'" in probe


def test_test_fixtures_that_duplicate_bundled_examples_stay_byte_identical():
    ex = ROOT / "app" / "backend" / "teagle_core" / "data" / "examples"
    shared = [p for p in (ROOT / "tests" / "fixtures").glob("*.fasta") if (ex / p.name).exists()]
    assert shared
    for p in shared:
        assert p.read_bytes() == (ex / p.name).read_bytes(), p.name


def test_walkthrough_multicopy_step_gives_a_second_product():
    """TRY-IT-YOURSELF: the specimen pasted twice as background makes one pair prime at both copies."""
    import engine
    from teagle_core import examples, sequtil
    seq = sequtil.parse_fasta(examples.load("M11240"))[0][1]
    rc = lambda s: s[::-1].translate(str.maketrans("ACGT", "TGCA"))
    res = engine.run_pcr({"sequence": ">x\n" + seq + "\n", "fwd": seq[1000:1021], "rev": rc(seq[1300:1321]),
                          "background": ">dup\n" + seq + seq + "\n"})
    keys = {a["product_key"] for a in res["amplicons"]}
    assert len(keys) == 3                      # the specimen plus one product per pasted copy
    bg = [a for a in res["amplicons"] if "custom background" in a["source"]]
    assert len(bg) == 2 and all(a["matches_intended_product"] is None for a in bg)


def test_clean_install_treats_a_missing_wsl_registry_key_as_nothing_to_clean():
    code = ISS.read_text(encoding="utf-8")
    fn = code[code.index("function WslDistroMayExist"):code.index("function CleanWslBackend")]
    assert "RegKeyExists(HKCU" in fn and fn.index("RegKeyExists") < fn.index("RegGetSubkeyNames")
    assert fn.index("Result := False;") < fn.index("RegKeyExists")


def test_second_launch_leaves_the_running_instances_log_alone(monkeypatch, tmp_path):
    log = tmp_path / "TEagle" / "logs" / "teagle.log"
    log.parent.mkdir(parents=True)
    log.write_bytes(b"x" * 1_100_000)
    fake = type(sys)("main")
    fake.main = lambda: pytest.fail("a second instance must not start")
    mod, _ = _native_launcher(monkeypatch, fake)
    monkeypatch.setattr(mod, "_windows_single_instance", lambda: True)
    monkeypatch.setattr(mod, "_raise_running_window", lambda *a: True)
    monkeypatch.setattr(sys, "platform", "win32")
    monkeypatch.setenv("LOCALAPPDATA", str(tmp_path))
    monkeypatch.setattr(sys, "stdout", None)
    monkeypatch.setattr(sys, "stderr", None)
    try:
        assert mod.main() == 0
    finally:
        _close_kept(mod)
    assert log.stat().st_size == 1_100_000
