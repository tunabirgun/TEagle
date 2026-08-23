from __future__ import annotations

import importlib.util
from pathlib import Path
import struct

import pytest

from installer_metadata_fixture import manifest as _manifest
from installer_metadata_fixture import synthetic_pe as _synthetic_pe


ROOT = Path(__file__).resolve().parents[1]
HELPER_PATH = ROOT / "installer" / "sanitize_inno_metadata.py"
BUILD_SCRIPT = ROOT / "installer" / "build_installer.ps1"


def _load_helper():
    spec = importlib.util.spec_from_file_location("teagle_installer_metadata_gate", HELPER_PATH)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_build_sanitizes_frozen_executable_before_selftest_and_packaging():
    source = BUILD_SCRIPT.read_text(encoding="utf-8")
    freeze = source.index('& $python -m PyInstaller "installer\\teagle_native.spec" --noconfirm')
    sanitize = source.index("& $python $metadataHelper sanitize-rich-header $dist")
    selftest = source.index('$env:TEAGLE_SELFTEST = "1"')
    compile_setup = source.index('& $sanitizedIscc "/DMyAppVersion=$ver"')
    assert freeze < sanitize < selftest < compile_setup


def test_manifest_patch_is_same_length_structural_and_confined():
    helper = _load_helper()
    source = _synthetic_pe(_manifest())
    result, spans = helper._patch_manifest(source, "Setup.e32")
    assert len(result) == len(source)
    assert helper._outside_spans(result, spans) == helper._outside_spans(source, spans)
    start, end, payload = helper._manifest_resource_span(result, "Setup.e32")
    assert result[start:end] == payload
    helper._validate_manifest(
        payload,
        "Setup.e32",
        helper.NEUTRAL_MANIFEST_IDENTITY,
        helper.NEUTRAL_MANIFEST_DESCRIPTION,
        0,
    )
    assert b"JR.Inno.Setup" not in payload
    assert b"Inno Setup" not in payload


@pytest.mark.parametrize(
    ("identity", "description", "extra"),
    [
        ("JR.Inno.Setup.Next", "Inno Setup", ""),
        ("JR-Inno-Setup", "Inno Setup", ""),
        ("JR.Inno.Setup", "InnoSetup Pro", ""),
        ("JR.Inno.Setup", "Inno Setup", "<!-- Built by Inno Installer -->"),
    ],
)
def test_manifest_patch_fails_closed_on_generator_variants(identity, description, extra):
    helper = _load_helper()
    source = _synthetic_pe(_manifest(identity, description, extra))
    with pytest.raises(ValueError, match="identity/description is not exact|generator-identifying residue"):
        helper._patch_manifest(source, "Setup.e32")


def test_manifest_validator_rejects_generator_residue_after_neutralisation():
    helper = _load_helper()
    payload = _manifest("TEagle.Client", "TEagle App", "<!-- Inno_Setup derivative -->")
    with pytest.raises(ValueError, match="generator-identifying residue"):
        helper._validate_manifest(payload, "setup", "TEagle.Client", "TEagle App", 0)


def test_manifest_resource_parser_rejects_unknown_code_page():
    helper = _load_helper()
    source = bytearray(_synthetic_pe(_manifest()))
    struct.pack_into("<I", source, 0x400 + 0x60 + 8, 1200)
    with pytest.raises(ValueError, match="unexpected code page"):
        helper._manifest_resource_span(bytes(source), "Setup.e32")


def test_rich_header_rewrite_preserves_layout_and_recomputes_checksum(tmp_path):
    helper = _load_helper()
    path = tmp_path / "TEagle.exe"
    source = _synthetic_pe(_manifest(), rich=True)
    path.write_bytes(source)
    rich_start, rich_end = helper._rich_header_span(source, path.name)
    checksum_offset = helper._pe_layout(source, path.name)["checksum_offset"]

    helper.sanitize_rich_header(path)

    result = path.read_bytes()
    assert len(result) == len(source)
    assert result[rich_start:rich_end] == b"\0" * (rich_end - rich_start)
    approved = [(rich_start, rich_end), (checksum_offset, checksum_offset + 4)]
    assert helper._outside_spans(result, approved) == helper._outside_spans(source, approved)
    stored_checksum = struct.unpack_from("<I", result, checksum_offset)[0]
    assert stored_checksum == helper._pe_checksum(result, checksum_offset)
    helper._pe_layout(result, path.name)
    with pytest.raises(ValueError, match="aligned Rich marker; found 0"):
        helper._rich_header_span(result, path.name)


@pytest.mark.parametrize("mutation", ["missing", "lowercase", "duplicate", "missing_dans"])
def test_rich_header_rewrite_fails_closed_on_mutations(tmp_path, mutation):
    helper = _load_helper()
    source = bytearray(_synthetic_pe(_manifest(), rich=True))
    rich_start, rich_end = helper._rich_header_span(bytes(source), "fixture.exe")
    rich_marker = rich_end - 8
    if mutation == "missing":
        source[rich_marker:rich_end] = b"\0" * 8
    elif mutation == "lowercase":
        source[rich_marker:rich_marker + 4] = b"RICH"
    elif mutation == "duplicate":
        source[0xC0:0xC4] = b"Rich"
        source[0xC4:0xC8] = source[rich_marker + 4:rich_end]
    else:
        source[rich_start:rich_start + 4] = b"\0" * 4
    path = tmp_path / "fixture.exe"
    path.write_bytes(source)
    before = path.read_bytes()
    with pytest.raises(ValueError, match="Rich marker|DanS marker"):
        helper.sanitize_rich_header(path)
    assert path.read_bytes() == before


def test_rich_header_rewrite_refuses_authenticode_table(tmp_path):
    helper = _load_helper()
    path = tmp_path / "signed.exe"
    source = _synthetic_pe(_manifest(), rich=True, certificate=True)
    path.write_bytes(source)
    with pytest.raises(ValueError, match="Authenticode-signed"):
        helper.sanitize_rich_header(path)
    assert path.read_bytes() == source
