"""Minimise optional metadata in copied Inno Setup release templates."""
from __future__ import annotations

import argparse
import ctypes
from pathlib import Path
import struct


TEMPLATE_NAMES = (
    "Setup.e32",
    "SetupCustomStyle.e32",
    "SetupLdr.e32",
    "SetupLdr.e64",
    "unins000.exe",
)
VERSIONED_TEMPLATE_NAMES = (
    "Setup.e32",
    "SetupCustomStyle.e32",
    "unins000.exe",
)
GENERATOR_COMMENT = "This installation was built with Inno Setup."
REQUIRED_STRING_VALUES = {
    "CompanyName": "Tuna Birgun",
    "FileDescription": "TEagle Setup",
    "ProductName": "TEagle",
}
OPTIONAL_STRING_KEYS = (
    "Comments",
    "InternalName",
    "LegalCopyright",
    "LegalTrademarks",
    "OriginalFilename",
    "PrivateBuild",
    "SpecialBuild",
)


class _VSFixedFileInfo(ctypes.Structure):
    _fields_ = [
        ("signature", ctypes.c_uint32),
        ("struct_version", ctypes.c_uint32),
        ("file_version_ms", ctypes.c_uint32),
        ("file_version_ls", ctypes.c_uint32),
        ("product_version_ms", ctypes.c_uint32),
        ("product_version_ls", ctypes.c_uint32),
        ("file_flags_mask", ctypes.c_uint32),
        ("file_flags", ctypes.c_uint32),
        ("file_os", ctypes.c_uint32),
        ("file_type", ctypes.c_uint32),
        ("file_subtype", ctypes.c_uint32),
        ("file_date_ms", ctypes.c_uint32),
        ("file_date_ls", ctypes.c_uint32),
    ]


def _atomic_write(path: Path, data: bytes) -> None:
    temporary = path.with_name(f".{path.name}.teagle-metadata-tmp")
    if temporary.exists():
        raise FileExistsError(f"metadata temporary already exists: {temporary}")
    try:
        temporary.write_bytes(data)
        temporary.replace(path)
    finally:
        if temporary.exists():
            temporary.unlink()


def _version_tuple(version_text: str) -> tuple[int, int, int, int]:
    parts = tuple(int(part) for part in version_text.split("."))
    if not 1 <= len(parts) <= 4 or any(part < 0 or part > 0xFFFF for part in parts):
        raise ValueError(f"invalid release version: {version_text}")
    return (parts + (0, 0, 0, 0))[:4]


def _replace_exact_span(data: bytes, old: bytes, new: bytes, label: str) -> tuple[bytes, tuple[int, int]]:
    if len(old) != len(new):
        raise ValueError(f"{label} replacement must preserve byte length")
    count = data.count(old)
    if count != 1:
        raise ValueError(f"{label} must occur exactly once; found {count}")
    offset = data.index(old)
    return data[:offset] + new + data[offset + len(old):], (offset, offset + len(old))


def _patch_uninstaller_version(data: bytes, version_text: str) -> tuple[bytes, list[tuple[int, int]]]:
    """Replace only the template's generator file version with the TEagle release version."""
    fixed_signature = struct.pack("<I", 0xFEEF04BD)
    signature_count = data.count(fixed_signature)
    if signature_count != 1:
        raise ValueError(f"uninstaller fixed version signature must occur exactly once; found {signature_count}")
    offset = data.index(fixed_signature)
    if len(data) - offset < 52:
        raise ValueError("uninstaller fixed version information is truncated")
    fixed = list(struct.unpack_from("<13I", data, offset))
    if fixed[1] != 0x00010000:
        raise ValueError("uninstaller fixed version structure is unexpected")
    old_version = (
        fixed[2] >> 16,
        fixed[2] & 0xFFFF,
        fixed[3] >> 16,
        fixed[3] & 0xFFFF,
    )
    new_version = _version_tuple(version_text)
    if old_version == new_version:
        raise ValueError("uninstaller template already uses the release version; version gate was not exercised")
    old_text = ".".join(str(part) for part in old_version).encode("utf-16le")
    new_label = ".".join(str(int(part)) for part in version_text.split("."))
    new_text_value = new_label.encode("utf-16le")
    if len(new_text_value) > len(old_text):
        raise ValueError("release version does not fit the uninstaller template FileVersion field")
    padded_new_text = new_text_value + b"\0" * (len(old_text) - len(new_text_value))
    data, text_span = _replace_exact_span(
        data, old_text, padded_new_text, "uninstaller FileVersion string"
    )
    fixed[2] = (new_version[0] << 16) | new_version[1]
    fixed[3] = (new_version[2] << 16) | new_version[3]
    old_fixed = struct.pack("<13I", *struct.unpack_from("<13I", data, offset))
    new_fixed = struct.pack("<13I", *fixed)
    data = data[:offset] + new_fixed + data[offset + len(new_fixed):]
    return data, [(offset, offset + len(old_fixed)), text_span]


def _outside_spans(data: bytes, spans: list[tuple[int, int]]) -> bytes:
    kept = bytearray()
    cursor = 0
    for start, end in sorted(spans):
        if start < cursor or end < start:
            raise ValueError("approved metadata byte spans overlap or are invalid")
        kept.extend(data[cursor:start])
        cursor = end
    kept.extend(data[cursor:])
    return bytes(kept)


def prepare_templates(source_dir: Path, copied_dir: Path, version_text: str) -> None:
    """Patch only known same-length template metadata spans without changing PE layout."""
    _version_tuple(version_text)
    needle = GENERATOR_COMMENT.encode("utf-16le")
    replacement = b"\0" * len(needle)
    for name in TEMPLATE_NAMES:
        source_path = source_dir / name
        copied_path = copied_dir / name
        if not source_path.is_file():
            raise FileNotFoundError(f"required installed Inno Setup template is missing: {source_path}")
        if not copied_path.is_file():
            raise FileNotFoundError(f"required copied Inno Setup template is missing: {copied_path}")
        source = source_path.read_bytes()
        copied = copied_path.read_bytes()
        if copied != source:
            raise ValueError(f"copied Inno template differs before metadata minimisation: {name}")
        expected, comment_span = _replace_exact_span(
            source, needle, replacement, f"{name} known generator Comments value"
        )
        spans = [comment_span]
        if name in VERSIONED_TEMPLATE_NAMES:
            expected, version_spans = _patch_uninstaller_version(expected, version_text)
            spans.extend(version_spans)
        _atomic_write(copied_path, expected)
        after = copied_path.read_bytes()
        if len(after) != len(source) or after != expected:
            raise ValueError(f"{name} differs from the verified metadata-minimised bytes")
        if _outside_spans(after, spans) != _outside_spans(source, spans):
            raise ValueError(f"{name} changed outside the approved byte spans")
        print(f"Minimised copied Inno template metadata: {name}")


def _query_version_strings(path: Path, keys: tuple[str, ...]) -> tuple[dict[str, str], tuple[int, ...]]:
    if not hasattr(ctypes, "windll"):
        raise RuntimeError("installer PE metadata verification requires Windows")
    version = ctypes.windll.version
    size = version.GetFileVersionInfoSizeW(str(path), None)
    if not size:
        raise ValueError(f"installer has no readable version resource: {path}")
    data = ctypes.create_string_buffer(size)
    if not version.GetFileVersionInfoW(str(path), 0, size, data):
        raise ValueError(f"could not read installer version resource: {path}")

    pointer = ctypes.c_void_p()
    length = ctypes.c_uint()
    if not version.VerQueryValueW(data, "\\VarFileInfo\\Translation", ctypes.byref(pointer), ctypes.byref(length)):
        raise ValueError("installer version resource has no translation table")
    if length.value < 4:
        raise ValueError("installer version translation table is truncated")
    raw_translation = ctypes.string_at(pointer, 4)
    language = int.from_bytes(raw_translation[:2], "little")
    code_page = int.from_bytes(raw_translation[2:], "little")
    block = f"{language:04x}{code_page:04x}"

    values: dict[str, str] = {}
    for key in keys:
        query = f"\\StringFileInfo\\{block}\\{key}"
        pointer = ctypes.c_void_p()
        length = ctypes.c_uint()
        if version.VerQueryValueW(data, query, ctypes.byref(pointer), ctypes.byref(length)) and length.value:
            values[key] = ctypes.wstring_at(pointer, length.value).rstrip("\0").strip()
        else:
            values[key] = ""

    pointer = ctypes.c_void_p()
    length = ctypes.c_uint()
    if not version.VerQueryValueW(data, "\\", ctypes.byref(pointer), ctypes.byref(length)):
        raise ValueError("installer version resource has no fixed information")
    if length.value < ctypes.sizeof(_VSFixedFileInfo):
        raise ValueError("installer fixed version information is truncated")
    fixed = ctypes.cast(pointer, ctypes.POINTER(_VSFixedFileInfo)).contents
    fixed_versions = (
        fixed.file_version_ms >> 16,
        fixed.file_version_ms & 0xFFFF,
        fixed.file_version_ls >> 16,
        fixed.file_version_ls & 0xFFFF,
        fixed.product_version_ms >> 16,
        fixed.product_version_ms & 0xFFFF,
        fixed.product_version_ls >> 16,
        fixed.product_version_ls & 0xFFFF,
    )
    return values, fixed_versions


def verify_setup(path: Path, version_text: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"installer is missing: {path}")
    version_keys = tuple(REQUIRED_STRING_VALUES) + ("FileVersion", "ProductVersion") + OPTIONAL_STRING_KEYS
    values, fixed_versions = _query_version_strings(path, version_keys)
    expected = dict(REQUIRED_STRING_VALUES)
    expected.update({"FileVersion": version_text, "ProductVersion": version_text})
    actual = {key: values[key] for key in expected}
    if actual != expected:
        raise ValueError(f"installer required PE metadata is not exact: {actual}")
    populated_optional = {key: values[key] for key in OPTIONAL_STRING_KEYS if values[key]}
    if populated_optional:
        raise ValueError(f"installer retains optional PE metadata: {populated_optional}")
    expected_fixed = _version_tuple(version_text)
    if fixed_versions != expected_fixed + expected_fixed:
        raise ValueError("installer fixed file/product version does not match the release version")
    if GENERATOR_COMMENT.encode("utf-16le") in path.read_bytes():
        raise ValueError("installer retains the Inno Setup generator Comments value")
    print(f"Verified minimal installer PE metadata: {path.name}")


def verify_uninstaller(path: Path, version_text: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"installed uninstaller is missing: {path}")
    version_keys = (
        "Comments", "CompanyName", "FileDescription", "FileVersion", "InternalName",
        "LegalCopyright", "LegalTrademarks", "OriginalFilename", "PrivateBuild",
        "ProductName", "ProductVersion", "SpecialBuild",
    )
    values, fixed_versions = _query_version_strings(path, version_keys)
    expected = {
        "CompanyName": "Tuna Birgun",
        "FileDescription": "Setup/Uninstall",
        "FileVersion": version_text,
        "ProductName": "TEagle",
        "ProductVersion": version_text,
    }
    actual = {key: values[key] for key in expected}
    if actual != expected:
        raise ValueError(f"installed uninstaller required PE metadata is not exact: {actual}")
    optional = set(version_keys) - set(expected)
    populated_optional = {key: values[key] for key in optional if values[key]}
    if populated_optional:
        raise ValueError(f"installed uninstaller retains optional PE metadata: {populated_optional}")
    expected_fixed = _version_tuple(version_text)
    if fixed_versions != expected_fixed + expected_fixed:
        raise ValueError("installed uninstaller fixed file/product version does not match the release version")
    raw = path.read_bytes()
    if GENERATOR_COMMENT.encode("utf-16le") in raw:
        raise ValueError("installed uninstaller retains the Inno Setup generator Comments value")
    print(f"Verified minimal installed-uninstaller PE metadata: {path.name}")


def main() -> int:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    prepare = subparsers.add_parser("prepare-templates")
    prepare.add_argument("source_dir", type=Path)
    prepare.add_argument("copied_dir", type=Path)
    prepare.add_argument("version")
    verify = subparsers.add_parser("verify-setup")
    verify.add_argument("path", type=Path)
    verify.add_argument("version")
    verify_uninstall = subparsers.add_parser("verify-uninstaller")
    verify_uninstall.add_argument("path", type=Path)
    verify_uninstall.add_argument("version")
    args = parser.parse_args()
    if args.command == "prepare-templates":
        prepare_templates(args.source_dir.resolve(), args.copied_dir.resolve(), args.version)
    elif args.command == "verify-setup":
        verify_setup(args.path.resolve(), args.version)
    else:
        verify_uninstaller(args.path.resolve(), args.version)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
