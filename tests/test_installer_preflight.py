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

from installer_metadata_fixture import manifest as _installer_manifest
from installer_metadata_fixture import synthetic_pe as _synthetic_installer_pe


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "installer" / "build_installer.ps1"
SPEC = ROOT / "installer" / "teagle_native.spec"
INNO_SCRIPT = ROOT / "installer" / "teagle.iss"
INNO_METADATA = ROOT / "installer" / "sanitize_inno_metadata.py"
IMPORT_GATE = ROOT / "installer" / "check_runtime_imports.py"
BUILD_LOCK = ROOT / "installer" / "requirements-build.txt"


def _locked_pyinstaller_version(path: Path = BUILD_LOCK) -> str:
    matches = re.findall(
        r"(?m)^\s*PyInstaller==([^\s#]+)\s*$",
        path.read_text(encoding="utf-8"),
    )
    assert len(matches) == 1, "requirements-build.txt must contain exactly one PyInstaller==<version> pin"
    return matches[0]


PYINSTALLER_PIN = _locked_pyinstaller_version()


def _locked_inno_version(path: Path = BUILD_LOCK) -> str:
    matches = re.findall(
        r"(?m)^\s*InnoSetup==([^\s#]+)\s*$",
        path.read_text(encoding="utf-8"),
    )
    assert len(matches) == 1, "requirements-build.txt must contain exactly one InnoSetup==<version> pin"
    return matches[0]


INNO_PIN = _locked_inno_version()


def _stage_preflight(
    tmp_path: Path, *, envcheck_exit: int, pyi_pin: str, inno_pin: str = INNO_PIN,
) -> Path:
    root = tmp_path / "fixture"
    installer = root / "installer"
    backend = root / "app" / "backend"
    core = backend / "teagle_core"
    installer.mkdir(parents=True)
    core.mkdir(parents=True)
    shutil.copy2(SOURCE, installer / SOURCE.name)
    shutil.copy2(IMPORT_GATE, installer / IMPORT_GATE.name)
    (installer / "requirements-build.txt").write_text(
        f"PyInstaller=={pyi_pin}\nInnoSetup=={inno_pin}\n", encoding="utf-8"
    )
    (core / "__init__.py").write_text('__version__ = "3.8.0"\n', encoding="utf-8")
    (backend / "requirements.txt").write_text(
        "primer3-py==2.3.0\npyhmmer==0.12.1\npsutil==7.2.2\nPySide6==6.11.1\n"
        "openpyxl==3.1.5\net_xmlfile==2.0.0\ncertifi==2026.7.22\npypdf==6.15.0\n",
        encoding="utf-8",
    )
    (backend / "envcheck.py").write_text(
        "import sys\n"
        "print('certifi required 2026.7.22; installed 2026.6.17' if "
        f"{envcheck_exit} else 'runtime pins and imports exact')\n"
        f"raise SystemExit({envcheck_exit})\n",
        encoding="utf-8",
    )
    (backend / "runtime_fixture.py").write_text(
        "import primer3, pyhmmer, psutil, PySide6, openpyxl, et_xmlfile, certifi, pypdf\n",
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


def test_unsupported_inno_setup_version_stops_before_build(tmp_path):
    script = _stage_preflight(
        tmp_path, envcheck_exit=0, pyi_pin=PYINSTALLER_PIN, inno_pin="0.0.0",
    )
    result = _run_preflight(script, sys.executable)
    assert result.returncode != 0
    assert "Inno Setup version mismatch" in result.stdout + result.stderr
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
    import_parity_gate = source.index("& $python $runtimeImportGate --root $root")
    runtime_gate = source.index("& $python $envcheck --no-install")
    build_tool_gate = source.index("PyInstaller version mismatch")
    bundle_cleanup = source.index('$bundle = Join-Path $root "dist\\TEagle"')
    freeze = source.index('& $python -m PyInstaller "installer\\teagle_native.spec" --noconfirm')
    assert import_parity_gate < runtime_gate < bundle_cleanup < freeze
    assert build_tool_gate < bundle_cleanup


def test_runtime_import_parity_gate_observably_fails_on_undeclared_import(tmp_path):
    root = tmp_path / "fixture"
    shutil.copytree(ROOT / "app" / "backend", root / "app" / "backend")
    shutil.copytree(ROOT / "app" / "native", root / "app" / "native")
    shutil.copy2(ROOT / "app" / "teagle.py", root / "app" / "teagle.py")
    (root / "installer").mkdir(parents=True)
    shutil.copy2(ROOT / "installer" / "teagle_native.py", root / "installer" / "teagle_native.py")
    shutil.copy2(IMPORT_GATE, root / "installer" / IMPORT_GATE.name)
    clean = subprocess.run(
        [sys.executable, str(root / "installer" / IMPORT_GATE.name), "--root", str(root)],
        capture_output=True, text=True, timeout=30,
    )
    assert clean.returncode == 0, clean.stdout + clean.stderr
    engine = root / "app" / "backend" / "engine.py"
    engine.write_text("import undeclared_runtime_package\n" + engine.read_text(encoding="utf-8"), encoding="utf-8")
    mutated = subprocess.run(
        [sys.executable, str(root / "installer" / IMPORT_GATE.name), "--root", str(root)],
        capture_output=True, text=True, timeout=30,
    )
    assert mutated.returncode != 0
    assert "undeclared runtime import undeclared_runtime_package" in mutated.stderr


def test_runtime_import_parity_scans_the_source_launcher(tmp_path):
    root = tmp_path / "fixture"
    shutil.copytree(ROOT / "app" / "backend", root / "app" / "backend")
    shutil.copytree(ROOT / "app" / "native", root / "app" / "native")
    shutil.copy2(ROOT / "app" / "teagle.py", root / "app" / "teagle.py")
    (root / "installer").mkdir(parents=True)
    shutil.copy2(ROOT / "installer" / "teagle_native.py", root / "installer" / "teagle_native.py")
    shutil.copy2(IMPORT_GATE, root / "installer" / IMPORT_GATE.name)
    launcher = root / "app" / "teagle.py"
    launcher.write_text("import source_launcher_dependency\n" + launcher.read_text(encoding="utf-8"), encoding="utf-8")
    result = subprocess.run(
        [sys.executable, str(root / "installer" / IMPORT_GATE.name), "--root", str(root)],
        capture_output=True, text=True, timeout=30,
    )
    assert result.returncode != 0
    assert "undeclared runtime import source_launcher_dependency" in result.stderr


def test_runtime_import_parity_rejects_unused_mapped_requirement(tmp_path):
    root = tmp_path / "fixture"
    shutil.copytree(ROOT / "app" / "backend", root / "app" / "backend")
    shutil.copytree(ROOT / "app" / "native", root / "app" / "native")
    shutil.copy2(ROOT / "app" / "teagle.py", root / "app" / "teagle.py")
    (root / "installer").mkdir(parents=True)
    shutil.copy2(ROOT / "installer" / "teagle_native.py", root / "installer" / "teagle_native.py")
    gate = root / "installer" / IMPORT_GATE.name
    shutil.copy2(IMPORT_GATE, gate)
    gate.write_text(
        gate.read_text(encoding="utf-8").replace(
            '    "pypdf": "pypdf",',
            '    "pypdf": "pypdf",\n    "unused-runtime": "unused_runtime",',
        ),
        encoding="utf-8",
    )
    requirements = root / "app" / "backend" / "requirements.txt"
    requirements.write_text(requirements.read_text(encoding="utf-8") + "unused-runtime==1.0\n", encoding="utf-8")
    result = subprocess.run(
        [sys.executable, str(gate), "--root", str(root)], capture_output=True, text=True, timeout=30,
    )
    assert result.returncode != 0
    assert "unused requirement unused-runtime" in result.stderr


@pytest.mark.parametrize("distribution,parent", [("psutil", "pyhmmer"), ("et_xmlfile", "openpyxl")])
def test_runtime_import_parity_requires_pinned_transitive_closure(tmp_path, distribution, parent):
    root = tmp_path / "fixture"
    shutil.copytree(ROOT / "app" / "backend", root / "app" / "backend")
    shutil.copytree(ROOT / "app" / "native", root / "app" / "native")
    shutil.copy2(ROOT / "app" / "teagle.py", root / "app" / "teagle.py")
    (root / "installer").mkdir(parents=True)
    shutil.copy2(ROOT / "installer" / "teagle_native.py", root / "installer" / "teagle_native.py")
    shutil.copy2(IMPORT_GATE, root / "installer" / IMPORT_GATE.name)
    requirements = root / "app" / "backend" / "requirements.txt"
    lines = [line for line in requirements.read_text(encoding="utf-8").splitlines()
             if not line.startswith(distribution + "==")]
    requirements.write_text("\n".join(lines) + "\n", encoding="utf-8")
    result = subprocess.run(
        [sys.executable, str(root / "installer" / IMPORT_GATE.name), "--root", str(root)],
        capture_output=True, text=True, timeout=30,
    )
    assert result.returncode != 0
    assert f"missing transitive requirement {distribution} (required by {parent})" in result.stderr


def test_frozen_spec_collects_and_copies_metadata_for_pinned_transitives():
    requirements = (ROOT / "app" / "backend" / "requirements.txt").read_text(encoding="utf-8")
    spec = SPEC.read_text(encoding="utf-8")
    assert "psutil==7.2.2" in requirements
    assert "et_xmlfile==2.0.0" in requirements
    assert '"psutil", "openpyxl", "et_xmlfile"' in spec
    assert '"pyhmmer", "psutil", "openpyxl", "et_xmlfile"' in spec
    assert '"primer3", "pyhmmer", "psutil"' in spec
    excludes = spec[spec.index("excludes = ["):spec.index("_icon =", spec.index("excludes = ["))]
    allowed = spec[spec.index("_ALLOWED_TOP = {"):spec.index("_counts =", spec.index("_ALLOWED_TOP = {"))]
    assert '"typing_extensions"' in excludes
    assert '"typing_extensions"' not in allowed


def _load_bundle_gate_helpers():
    tree = ast.parse(SPEC.read_text(encoding="utf-8"), filename=str(SPEC))
    helper_names = {
        "_is_forbidden_bundle_path",
        "_without_forbidden_bundle_entries",
        "_assert_bundle_inventory_clean",
        "_bundle_inventory",
    }
    nodes = [
        node for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name in helper_names
    ]
    namespace = {"os": os}
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(SPEC), "exec"), namespace)
    return namespace


def test_frozen_spec_filters_unused_qt_native_chains_before_collect():
    source = SPEC.read_text(encoding="utf-8")
    analysis = source.index("a = Analysis(")
    filter_binaries = source.index("a.binaries = _without_forbidden_bundle_entries(a.binaries)")
    filter_datas = source.index("a.datas = _without_forbidden_bundle_entries(a.datas)")
    collect = source.index("coll = COLLECT(")
    assert analysis < filter_binaries < collect
    assert analysis < filter_datas < collect

    helpers = _load_bundle_gate_helpers()
    entries = [
        ("PySide6/Qt6Core.dll", "core", "BINARY"),
        ("PySide6/Qt6Pdf.dll", "pdf", "BINARY"),
        ("PySide6/Qt6QuickWidgets.dll", "quick", "BINARY"),
        ("PySide6/Qt6QmlModels.dll", "qml", "BINARY"),
        ("PySide6/Qt6VirtualKeyboard.dll", "keyboard", "BINARY"),
        ("PySide6/plugins/imageformats/qpdf.dll", "qpdf", "BINARY"),
        ("PySide6/plugins/platforminputcontexts/qtvirtualkeyboardplugin.dll", "vk", "BINARY"),
        ("PySide6/qml/QtQuick/qtquick2plugin.dll", "qml-plugin", "BINARY"),
    ]
    assert helpers["_without_forbidden_bundle_entries"](entries) == [entries[0]]


@pytest.mark.parametrize(
    "mutation",
    [
        "_internal\\PySide6\\plugins\\imageformats\\qpdf.dll",
        "_internal/PySide6/Qt6Pdf.dll",
        "_internal/PySide6/Qt6Quick.dll",
        "_internal/PySide6/Qt6QmlWorkerScript.dll",
        "_internal/PySide6/Qt6VirtualKeyboard.dll",
        "_internal/PySide6/qml/QtQuick/qtquick2plugin.dll",
        "_internal/PySide6/plugins/platforminputcontexts/qtvirtualkeyboardplugin.dll",
    ],
)
def test_frozen_bundle_inventory_gate_fails_closed_on_unused_qt_mutation(mutation):
    helpers = _load_bundle_gate_helpers()
    with pytest.raises(SystemExit, match="BUNDLE GUARD: forbidden unused Qt payload"):
        helpers["_assert_bundle_inventory_clean"](["_internal/PySide6/Qt6Core.dll", mutation], "fixture")


def test_frozen_spec_checks_collect_table_and_assembled_inventory_after_collect():
    source = SPEC.read_text(encoding="utf-8")
    collect = source.index("coll = COLLECT(")
    toc_gate = source.index(
        '_assert_bundle_inventory_clean((entry[0] for entry in coll.toc), "the COLLECT table")'
    )
    disk_gate = source.index(
        '_assert_bundle_inventory_clean(_bundle_inventory(coll.name), "the assembled bundle")'
    )
    assert collect < toc_gate < disk_gate


def test_frozen_executable_version_info_omits_optional_identity_strings():
    source = SPEC.read_text(encoding="utf-8")
    assert "StringStruct('InternalName'" not in source
    assert "StringStruct('OriginalFilename'" not in source


def test_inno_setup_compiler_version_is_exactly_locked_and_gated_preflight():
    lock = BUILD_LOCK.read_text(encoding="utf-8")
    matches = re.findall(r"(?m)^InnoSetup==([^\s#]+)\s*$", lock)
    assert matches == [INNO_PIN] == ["6.7.3"]
    source = SOURCE.read_text(encoding="utf-8")
    assert "Inno Setup version mismatch" in source
    assert source.index("Inno Setup version mismatch") < source.index("if ($PreflightOnly)")
    assert source.index("if ($PreflightOnly)") < source.index('$bundle = Join-Path $root "dist\\TEagle"')


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
    old_file_version = "51.1054.0.0".encode("utf-16le")
    payload = _synthetic_installer_pe(_installer_manifest())
    for name in helper.TEMPLATE_NAMES:
        (source_dir / name).write_bytes(payload)
        (copied_dir / name).write_bytes(payload)

    helper.prepare_templates(source_dir, copied_dir, "3.8.1")

    for name in helper.TEMPLATE_NAMES:
        result = (copied_dir / name).read_bytes()
        assert len(result) == len(payload)
        assert needle not in result
        comment_start = payload.index(needle)
        assert result[comment_start:comment_start + len(needle)] == b"\0" * len(needle)
        _, _, manifest = helper._manifest_resource_span(result, name)
        helper._validate_manifest(
            manifest, name, helper.NEUTRAL_MANIFEST_IDENTITY,
            helper.NEUTRAL_MANIFEST_DESCRIPTION, 0,
        )
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
    payload = _synthetic_installer_pe(_installer_manifest(), comment_copies=copies)
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
