"""Hermetic gates for the Windows packaging preflight; no bundle is built."""
from __future__ import annotations

import ast
import importlib.util
import os
import re
import shutil
import struct
import subprocess
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "installer" / "build_installer.ps1"
SPEC = ROOT / "installer" / "teagle_native.spec"
INNO_SCRIPT = ROOT / "installer" / "teagle.iss"
INNO_METADATA = ROOT / "installer" / "sanitize_inno_metadata.py"
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


def test_inno_metadata_minimisation_is_gated_before_and_after_compile():
    source = SOURCE.read_text(encoding="utf-8")
    inno = INNO_SCRIPT.read_text(encoding="utf-8")
    assert "prepare-templates $innoSource $innoTemp $ver" in source
    assert "verify-setup $setupOutput $ver" in source
    assert source.index("prepare-templates $innoSource $innoTemp $ver") < source.index("& $sanitizedIscc")
    assert source.index("& $sanitizedIscc") < source.index("verify-setup $setupOutput $ver")
    assert "DisablePrecompiledFileVerifications=setup setupcustomstyle setupldr" in inno
    assert INNO_METADATA.is_file()


def test_inno_metadata_helper_has_narrow_fail_closed_contract():
    source = INNO_METADATA.read_text(encoding="utf-8")
    tree = ast.parse(source, filename=str(INNO_METADATA))
    assignments = {
        target.id: ast.literal_eval(node.value)
        for node in tree.body
        if isinstance(node, (ast.Assign, ast.AnnAssign))
        for target in (
            node.targets if isinstance(node, ast.Assign) else [node.target]
        )
        if isinstance(target, ast.Name)
        and target.id in {
            "TEMPLATE_NAMES", "VERSIONED_TEMPLATE_NAMES", "GENERATOR_COMMENT", "OPTIONAL_STRING_KEYS",
        }
    }
    assert assignments["TEMPLATE_NAMES"] == (
        "Setup.e32", "SetupCustomStyle.e32", "SetupLdr.e32", "SetupLdr.e64", "unins000.exe",
    )
    assert assignments["VERSIONED_TEMPLATE_NAMES"] == (
        "Setup.e32", "SetupCustomStyle.e32", "unins000.exe",
    )
    assert assignments["GENERATOR_COMMENT"] == "This installation was built with Inno Setup."
    assert assignments["OPTIONAL_STRING_KEYS"] == (
        "Comments", "InternalName", "LegalCopyright", "LegalTrademarks",
        "OriginalFilename", "PrivateBuild", "SpecialBuild",
    )
    assert "must occur exactly once; found" in source
    assert "changed outside the approved byte span" in source
    assert "installer retains optional PE metadata" in source
    assert "installer retains the Inno Setup generator Comments value" in source
    assert "uninstaller template already uses the release version" in source
    assert "installed uninstaller retains optional PE metadata" in source


def _load_inno_metadata_helper():
    spec = importlib.util.spec_from_file_location("teagle_sanitize_inno_metadata", INNO_METADATA)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_inno_template_cleaner_changes_only_one_exact_same_length_span(tmp_path):
    helper = _load_inno_metadata_helper()
    source_dir = tmp_path / "source"
    copied_dir = tmp_path / "copied"
    source_dir.mkdir()
    copied_dir.mkdir()
    needle = helper.GENERATOR_COMMENT.encode("utf-16le")
    fixed = struct.pack(
        "<13I", 0xFEEF04BD, 0x00010000, (51 << 16) | 1054, 0,
        0, 0, 0x3F, 0, 0x00040004, 1, 0, 0, 0,
    )
    old_file_version = "51.1054.0.0".encode("utf-16le")
    payload = b"prefix-resource-bytes" + fixed + old_file_version + needle + b"suffix-resource-bytes"
    for name in helper.TEMPLATE_NAMES:
        (source_dir / name).write_bytes(payload)
        (copied_dir / name).write_bytes(payload)

    helper.prepare_templates(source_dir, copied_dir, "3.8.1")

    for name in helper.TEMPLATE_NAMES:
        result = (copied_dir / name).read_bytes()
        assert len(result) == len(payload)
        assert result[:len(b"prefix-resource-bytes")] == b"prefix-resource-bytes"
        assert result[-len(b"suffix-resource-bytes"):] == b"suffix-resource-bytes"
        assert needle not in result
        comment_start = payload.index(needle)
        assert result[comment_start:comment_start + len(needle)] == b"\0" * len(needle)
        assert (source_dir / name).read_bytes() == payload
    for name in helper.VERSIONED_TEMPLATE_NAMES:
        result = (copied_dir / name).read_bytes()
        version_offset = result.index(struct.pack("<I", 0xFEEF04BD))
        version_fields = struct.unpack_from("<13I", result, version_offset)
        assert version_fields[2:4] == ((3 << 16) | 8, 1 << 16)
        assert "3.8.1".encode("utf-16le") in result
        assert old_file_version not in result


@pytest.mark.parametrize("copies", [0, 2])
def test_inno_template_cleaner_fails_closed_when_known_comment_count_drifts(tmp_path, copies):
    helper = _load_inno_metadata_helper()
    source_dir = tmp_path / "source"
    copied_dir = tmp_path / "copied"
    source_dir.mkdir()
    copied_dir.mkdir()
    needle = helper.GENERATOR_COMMENT.encode("utf-16le")
    fixed = struct.pack(
        "<13I", 0xFEEF04BD, 0x00010000, (51 << 16) | 1054, 0,
        0, 0, 0x3F, 0, 0x00040004, 1, 0, 0, 0,
    )
    payload = b"prefix" + fixed + "51.1054.0.0".encode("utf-16le") + (needle * copies) + b"suffix"
    for name in helper.TEMPLATE_NAMES:
        (source_dir / name).write_bytes(payload)
        (copied_dir / name).write_bytes(payload)
    with pytest.raises(ValueError, match="known generator Comments value must occur exactly once"):
        helper.prepare_templates(source_dir, copied_dir, "3.8.1")
    assert (copied_dir / helper.TEMPLATE_NAMES[0]).read_bytes() == payload


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
