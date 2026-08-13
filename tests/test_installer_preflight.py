"""Hermetic gates for the Windows packaging preflight; no bundle is built."""
from __future__ import annotations

import ast
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "installer" / "build_installer.ps1"
SPEC = ROOT / "installer" / "teagle_native.spec"
BUILD_LOCK = ROOT / "installer" / "requirements-build.txt"


def _locked_pyinstaller_version(path: Path = BUILD_LOCK) -> str:
    matches = re.findall(
        r"(?m)^\s*PyInstaller==([^\s#]+)\s*$",
        path.read_text(encoding="utf-8"),
    )
    assert len(matches) == 1, "requirements-build.txt must contain exactly one PyInstaller==<version> pin"
    return matches[0]


PYINSTALLER_PIN = _locked_pyinstaller_version()


def _stage_preflight(tmp_path: Path, *, envcheck_exit: int, pyi_pin: str) -> Path:
    root = tmp_path / "fixture"
    installer = root / "installer"
    backend = root / "app" / "backend"
    core = backend / "teagle_core"
    installer.mkdir(parents=True)
    core.mkdir(parents=True)
    shutil.copy2(SOURCE, installer / SOURCE.name)
    (installer / "requirements-build.txt").write_text(
        f"PyInstaller=={pyi_pin}\n", encoding="utf-8"
    )
    (core / "__init__.py").write_text('__version__ = "3.8.0"\n', encoding="utf-8")
    (backend / "envcheck.py").write_text(
        "import sys\n"
        "print('certifi required 2026.7.22; installed 2026.6.17' if "
        f"{envcheck_exit} else 'runtime pins and imports exact')\n"
        f"raise SystemExit({envcheck_exit})\n",
        encoding="utf-8",
    )
    return installer / SOURCE.name


def _run_preflight(script: Path, python_exe: Path | str | None, *, env=None) -> subprocess.CompletedProcess[str]:
    powershell = shutil.which("powershell")
    if not powershell:
        pytest.skip("Windows PowerShell is unavailable")
    command = [powershell, "-NoProfile", "-File", str(script), "-PreflightOnly"]
    if python_exe is not None:
        command += ["-PythonExe", str(python_exe)]
    return subprocess.run(
        command,
        capture_output=True,
        env=env,
        text=True,
        timeout=30,
    )


def test_runtime_pin_or_import_failure_stops_before_build(tmp_path):
    script = _stage_preflight(
        tmp_path,
        envcheck_exit=1,
        pyi_pin=PYINSTALLER_PIN,
    )
    result = _run_preflight(script, sys.executable)
    assert result.returncode != 0
    assert "Runtime dependency preflight failed" in result.stdout + result.stderr
    assert not (script.parents[1] / "dist").exists()


def test_stale_pyinstaller_stops_before_build(tmp_path):
    script = _stage_preflight(tmp_path, envcheck_exit=0, pyi_pin="0.0.0")
    result = _run_preflight(script, sys.executable)
    assert result.returncode != 0
    assert "PyInstaller version mismatch" in result.stdout + result.stderr
    assert not (script.parents[1] / "dist").exists()


def test_absent_python_executable_cannot_pass(tmp_path):
    script = _stage_preflight(
        tmp_path,
        envcheck_exit=0,
        pyi_pin=PYINSTALLER_PIN,
    )
    result = _run_preflight(script, tmp_path / "missing-python.exe")
    assert result.returncode != 0
    output = " ".join((result.stdout + result.stderr).split())
    assert "not a runnable Python 3.12 executable" in output
    assert not (script.parents[1] / "dist").exists()


def test_exact_preflight_can_pass_without_building(tmp_path):
    script = _stage_preflight(
        tmp_path,
        envcheck_exit=0,
        pyi_pin=PYINSTALLER_PIN,
    )
    result = _run_preflight(script, sys.executable)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "Preflight passed; no build requested" in result.stdout
    assert not (script.parents[1] / "dist").exists()


def test_localappdata_python_discovery_works_without_python_on_path(tmp_path):
    local_app_data = os.environ.get("LOCALAPPDATA")
    if not local_app_data or not Path(sys.executable).is_relative_to(Path(local_app_data)):
        pytest.skip("test interpreter is not a per-user LOCALAPPDATA Python")
    script = _stage_preflight(
        tmp_path,
        envcheck_exit=0,
        pyi_pin=PYINSTALLER_PIN,
    )
    clean_env = os.environ.copy()
    clean_env["PATH"] = ""
    result = _run_preflight(script, None, env=clean_env)
    assert result.returncode == 0, result.stdout + result.stderr
    assert str(Path(sys.executable).resolve()) in result.stdout


@pytest.mark.parametrize(
    "contents",
    [
        "",
        "PyInstaller>=6.21.0\n",
        "PyInstaller==6.21.0\nPyInstaller==6.21.1\n",
    ],
)
def test_build_lock_parser_rejects_missing_loose_or_duplicate_pins(tmp_path, contents):
    lock = tmp_path / "requirements-build.txt"
    lock.write_text(contents, encoding="utf-8")
    with pytest.raises(AssertionError, match="exactly one PyInstaller"):
        _locked_pyinstaller_version(lock)


def test_full_release_path_cannot_skip_freeze_and_orders_gates_before_mutation():
    source = SOURCE.read_text(encoding="utf-8")
    assert "SkipFreeze" not in source
    runtime_gate = source.index("& $python $envcheck --no-install")
    build_tool_gate = source.index("PyInstaller version mismatch")
    bundle_cleanup = source.index('$bundle = Join-Path $root "dist\\TEagle"')
    freeze = source.index('& $python -m PyInstaller "installer\\teagle_native.spec" --noconfirm')
    assert runtime_gate < bundle_cleanup < freeze
    assert build_tool_gate < bundle_cleanup


def test_installer_spec_filters_dependency_tests_without_dropping_runtime_payloads():
    """Exercise the spec's collector without running Analysis or building a bundle."""
    tree = ast.parse(SPEC.read_text(encoding="utf-8"), filename=str(SPEC))
    helper_names = {"_is_test_payload", "_has_test_filename", "_collect_runtime_package"}
    nodes = [
        node
        for node in tree.body
        if (
            isinstance(node, ast.Assign)
            and any(
                isinstance(target, ast.Name) and target.id == "_TEST_COMPONENTS"
                for target in node.targets
            )
        )
        or (isinstance(node, ast.FunctionDef) and node.name in helper_names)
    ]

    runtime_datas = [
        ("pyhmmer-runtime.pyd", "pyhmmer"),
        ("cacert.pem", "certifi"),
        ("primer3-runtime.c", "primer3/src/libprimer3"),
        ("latest_runtime.c", "primer3/src/libprimer3"),
    ]
    test_datas = [
        ("test_hmmer.py", "pyhmmer/tests"),
        ("test_certify.py", "certifi/tests"),
        ("khash_test.c", "primer3/src/libprimer3/klib/test"),
        ("long_seq_tm_test_main.c", "primer3/src/libprimer3"),
    ]
    runtime_modules = ["pyhmmer.plan7", "certifi.core", "primer3.bindings", "openpyxl.workbook"]
    test_modules = ["pyhmmer.tests.test_plan7", "certifi.tests.test_certify"]
    namespace = {
        "os": os,
        "re": re,
        "collect_all": lambda package: (
            runtime_datas + test_datas,
            [("runtime.dll", ".")],
            runtime_modules + test_modules,
        ),
        "collect_submodules": lambda package: ["pyhmmer.hmmer", "pyhmmer.tests.test_hmmer"],
    }
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(SPEC), "exec"), namespace)

    datas, binaries, hiddenimports = namespace["_collect_runtime_package"]("fixture")
    assert datas == runtime_datas
    assert binaries == [("runtime.dll", ".")]
    assert hiddenimports == runtime_modules + ["pyhmmer.hmmer"]
    assert "d, b, h = _collect_runtime_package(pkg)" in SPEC.read_text(encoding="utf-8")
