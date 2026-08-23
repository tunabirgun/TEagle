"""Minimise optional metadata in copied Inno Setup release templates."""
from __future__ import annotations

import argparse
import ctypes
from pathlib import Path
import re
import struct
import xml.etree.ElementTree as ET


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
SOURCE_MANIFEST_IDENTITY = "JR.Inno.Setup"
NEUTRAL_MANIFEST_IDENTITY = "TEagle.Client"
SOURCE_MANIFEST_DESCRIPTION = "Inno Setup"
NEUTRAL_MANIFEST_DESCRIPTION = "TEagle App"
_GENERATOR_MANIFEST_TOKEN = re.compile(
    r"(?i)\b(?:jr[.\s_-]*)?inno(?:[.\s_-]*setup|[.\s_-]*installer)?\b"
)
_ASSEMBLY_V1 = "urn:schemas-microsoft-com:asm.v1"
_RT_MANIFEST = 24
_MANIFEST_RESOURCE_ID = 1
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


def _pe_layout(data: bytes, label: str) -> dict[str, object]:
    if len(data) < 0x40 or data[:2] != b"MZ":
        raise ValueError(f"{label} is not a PE file")
    pe_offset = struct.unpack_from("<I", data, 0x3C)[0]
    if pe_offset < 0x40 or pe_offset + 24 > len(data) or data[pe_offset:pe_offset + 4] != b"PE\0\0":
        raise ValueError(f"{label} has an invalid PE header")
    coff_offset = pe_offset + 4
    section_count = struct.unpack_from("<H", data, coff_offset + 2)[0]
    optional_size = struct.unpack_from("<H", data, coff_offset + 16)[0]
    optional_offset = coff_offset + 20
    if optional_offset + optional_size > len(data):
        raise ValueError(f"{label} has a truncated PE optional header")
    magic = struct.unpack_from("<H", data, optional_offset)[0]
    if magic == 0x10B:
        directory_count_offset, directory_offset = 92, 96
    elif magic == 0x20B:
        directory_count_offset, directory_offset = 108, 112
    else:
        raise ValueError(f"{label} has an unsupported PE optional-header magic")
    if optional_size < directory_offset:
        raise ValueError(f"{label} has no PE data-directory table")
    directory_count = struct.unpack_from("<I", data, optional_offset + directory_count_offset)[0]
    available_directories = (optional_size - directory_offset) // 8
    if directory_count > available_directories:
        raise ValueError(f"{label} PE data-directory count exceeds its optional header")
    section_offset = optional_offset + optional_size
    if not section_count or section_offset + section_count * 40 > len(data):
        raise ValueError(f"{label} has a truncated PE section table")
    sections: list[tuple[int, int, int, int]] = []
    for index in range(section_count):
        offset = section_offset + index * 40
        virtual_size, virtual_address, raw_size, raw_offset = struct.unpack_from("<IIII", data, offset + 8)
        if raw_size and raw_offset + raw_size > len(data):
            raise ValueError(f"{label} has a section outside the file")
        sections.append((virtual_address, virtual_size, raw_offset, raw_size))
    directories: list[tuple[int, int]] = []
    for index in range(directory_count):
        directories.append(struct.unpack_from("<II", data, optional_offset + directory_offset + index * 8))
    certificate = directories[4] if len(directories) > 4 else (0, 0)
    if bool(certificate[0]) != bool(certificate[1]):
        raise ValueError(f"{label} has an inconsistent PE certificate-table entry")
    if certificate[0] and certificate[0] + certificate[1] > len(data):
        raise ValueError(f"{label} has a certificate table outside the file")
    return {
        "pe_offset": pe_offset,
        "checksum_offset": optional_offset + 64,
        "directories": directories,
        "sections": sections,
        "certificate": certificate,
    }


def _rva_to_file_offset(
    data: bytes, layout: dict[str, object], rva: int, size: int, label: str
) -> int:
    if size < 0:
        raise ValueError(f"{label} has a negative resource size")
    for virtual_address, virtual_size, raw_offset, raw_size in layout["sections"]:
        mapped_size = max(virtual_size, raw_size)
        if virtual_address <= rva and rva + size <= virtual_address + mapped_size:
            delta = rva - virtual_address
            if delta + size > raw_size:
                raise ValueError(f"{label} resource extends beyond its raw PE section")
            return raw_offset + delta
    raise ValueError(f"{label} resource RVA is not mapped by a PE section")


def _manifest_resource_span(data: bytes, label: str) -> tuple[int, int, bytes]:
    layout = _pe_layout(data, label)
    directories = layout["directories"]
    if len(directories) <= 2 or not directories[2][0] or not directories[2][1]:
        raise ValueError(f"{label} has no PE resource directory")
    resource_rva, resource_size = directories[2]
    resource_base = _rva_to_file_offset(data, layout, resource_rva, resource_size, label)
    resource_end = resource_base + resource_size
    if resource_end > len(data):
        raise ValueError(f"{label} PE resource directory is truncated")

    def entries(relative_offset: int) -> list[tuple[int, int]]:
        offset = resource_base + relative_offset
        if offset < resource_base or offset + 16 > resource_end:
            raise ValueError(f"{label} has an invalid PE resource directory offset")
        named, identified = struct.unpack_from("<HH", data, offset + 12)
        count = named + identified
        if offset + 16 + count * 8 > resource_end:
            raise ValueError(f"{label} has a truncated PE resource directory")
        return [struct.unpack_from("<II", data, offset + 16 + i * 8) for i in range(count)]

    def one_id(level: list[tuple[int, int]], wanted: int, context: str) -> int:
        matches = [target for name, target in level if not name & 0x80000000 and name == wanted]
        if len(matches) != 1:
            raise ValueError(f"{label} must contain exactly one {context}; found {len(matches)}")
        return matches[0]

    type_target = one_id(entries(0), _RT_MANIFEST, "RT_MANIFEST resource type")
    if not type_target & 0x80000000:
        raise ValueError(f"{label} RT_MANIFEST type does not point to a directory")
    name_target = one_id(entries(type_target & 0x7FFFFFFF), _MANIFEST_RESOURCE_ID, "manifest resource ID 1")
    if not name_target & 0x80000000:
        raise ValueError(f"{label} manifest resource ID does not point to a language directory")
    languages = entries(name_target & 0x7FFFFFFF)
    if len(languages) != 1:
        raise ValueError(f"{label} must contain exactly one manifest language; found {len(languages)}")
    language_name, data_target = languages[0]
    if language_name & 0x80000000 or data_target & 0x80000000:
        raise ValueError(f"{label} manifest language entry is structurally invalid")
    if language_name != 1033:
        raise ValueError(f"{label} manifest has an unexpected language ID: {language_name}")
    data_entry = resource_base + data_target
    if data_entry < resource_base or data_entry + 16 > resource_end:
        raise ValueError(f"{label} manifest data entry is outside the resource directory")
    manifest_rva, manifest_size, code_page, reserved = struct.unpack_from("<IIII", data, data_entry)
    if not manifest_size or reserved:
        raise ValueError(f"{label} manifest data entry is invalid")
    if code_page not in (0, 1252):
        raise ValueError(f"{label} manifest has an unexpected code page: {code_page}")
    start = _rva_to_file_offset(data, layout, manifest_rva, manifest_size, label)
    end = start + manifest_size
    return start, end, data[start:end]


def _parse_manifest(payload: bytes, label: str) -> tuple[ET.Element, ET.Element, ET.Element, str]:
    try:
        text = payload.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise ValueError(f"{label} manifest is not strict UTF-8") from exc
    try:
        root = ET.fromstring(text)
    except ET.ParseError as exc:
        raise ValueError(f"{label} manifest is not structurally valid XML") from exc
    assembly_tag = f"{{{_ASSEMBLY_V1}}}assembly"
    identity_tag = f"{{{_ASSEMBLY_V1}}}assemblyIdentity"
    description_tag = f"{{{_ASSEMBLY_V1}}}description"
    if root.tag != assembly_tag or root.attrib != {"manifestVersion": "1.0"}:
        raise ValueError(f"{label} manifest assembly root is unexpected")
    identities = [child for child in root if child.tag == identity_tag]
    descriptions = [child for child in root if child.tag == description_tag]
    if len(identities) != 1 or len(descriptions) != 1:
        raise ValueError(f"{label} manifest must have one assemblyIdentity and one description")
    identity = identities[0]
    description = descriptions[0]
    required_identity = {"name", "processorArchitecture", "version", "type"}
    if set(identity.attrib) != required_identity:
        raise ValueError(f"{label} manifest assemblyIdentity attributes are unexpected")
    if identity.attrib["processorArchitecture"] not in {"x86", "amd64"}:
        raise ValueError(f"{label} manifest processor architecture is unexpected")
    if identity.attrib["version"] != "1.0.0.0" or identity.attrib["type"] != "win32":
        raise ValueError(f"{label} manifest identity semantics are unexpected")
    return root, identity, description, text


def _manifest_shape(element: ET.Element) -> tuple[object, ...]:
    attributes = dict(element.attrib)
    text = element.text
    if element.tag == f"{{{_ASSEMBLY_V1}}}assemblyIdentity":
        attributes["name"] = "<application-identity>"
    elif element.tag == f"{{{_ASSEMBLY_V1}}}description":
        text = "<application-description>"
    return (
        element.tag,
        tuple(sorted(attributes.items())),
        text,
        element.tail,
        tuple(_manifest_shape(child) for child in element),
    )


def _validate_manifest(
    payload: bytes, label: str, expected_identity: str, expected_description: str,
    expected_generator_matches: int,
) -> ET.Element:
    root, identity, description, text = _parse_manifest(payload, label)
    if identity.attrib["name"] != expected_identity or description.text != expected_description:
        raise ValueError(
            f"{label} manifest identity/description is not exact: "
            f"{identity.attrib['name']!r}, {description.text!r}"
        )
    matches = _GENERATOR_MANIFEST_TOKEN.findall(text)
    if len(matches) != expected_generator_matches:
        raise ValueError(
            f"{label} manifest has unexpected generator-identifying residue; "
            f"expected {expected_generator_matches} exact source fields, found {len(matches)}"
        )
    return root


def _patch_manifest(data: bytes, label: str) -> tuple[bytes, list[tuple[int, int]]]:
    start, end, payload = _manifest_resource_span(data, label)
    before = _validate_manifest(
        payload, label, SOURCE_MANIFEST_IDENTITY, SOURCE_MANIFEST_DESCRIPTION, 2
    )
    patched, identity_span = _replace_exact_span(
        payload,
        SOURCE_MANIFEST_IDENTITY.encode("utf-8"),
        NEUTRAL_MANIFEST_IDENTITY.encode("utf-8"),
        f"{label} manifest assembly identity",
    )
    patched, description_span = _replace_exact_span(
        patched,
        SOURCE_MANIFEST_DESCRIPTION.encode("utf-8"),
        NEUTRAL_MANIFEST_DESCRIPTION.encode("utf-8"),
        f"{label} manifest description",
    )
    after = _validate_manifest(
        patched, label, NEUTRAL_MANIFEST_IDENTITY, NEUTRAL_MANIFEST_DESCRIPTION, 0
    )
    if len(patched) != len(payload) or _manifest_shape(before) != _manifest_shape(after):
        raise ValueError(f"{label} manifest structure or non-identifying semantics changed")
    result = data[:start] + patched + data[end:]
    return result, [
        (start + identity_span[0], start + identity_span[1]),
        (start + description_span[0], start + description_span[1]),
    ]


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
        expected, manifest_spans = _patch_manifest(source, name)
        expected, comment_span = _replace_exact_span(
            expected, needle, replacement, f"{name} known generator Comments value"
        )
        spans = manifest_spans + [comment_span]
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
    raw = path.read_bytes()
    _, _, manifest = _manifest_resource_span(raw, path.name)
    _validate_manifest(
        manifest, path.name, NEUTRAL_MANIFEST_IDENTITY, NEUTRAL_MANIFEST_DESCRIPTION, 0
    )
    if GENERATOR_COMMENT.encode("utf-16le") in raw:
        raise ValueError("installer retains the Inno Setup generator Comments value")
    print(f"Verified minimal installer PE metadata and RT_MANIFEST: {path.name}")


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
    _, _, manifest = _manifest_resource_span(raw, path.name)
    _validate_manifest(
        manifest, path.name, NEUTRAL_MANIFEST_IDENTITY, NEUTRAL_MANIFEST_DESCRIPTION, 0
    )
    if GENERATOR_COMMENT.encode("utf-16le") in raw:
        raise ValueError("installed uninstaller retains the Inno Setup generator Comments value")
    print(f"Verified minimal installed-uninstaller PE metadata and RT_MANIFEST: {path.name}")


def _rich_header_span(data: bytes, label: str) -> tuple[int, int]:
    layout = _pe_layout(data, label)
    pe_offset = layout["pe_offset"]
    rich_offsets = [
        offset for offset in range(0x40, pe_offset - 7, 4)
        if data[offset:offset + 4] == b"Rich"
    ]
    if len(rich_offsets) != 1:
        raise ValueError(f"{label} must contain exactly one aligned Rich marker; found {len(rich_offsets)}")
    rich_offset = rich_offsets[0]
    xor_key = struct.unpack_from("<I", data, rich_offset + 4)[0]
    encoded_dans = 0x536E6144 ^ xor_key
    starts = [
        offset for offset in range(0x40, rich_offset, 4)
        if struct.unpack_from("<I", data, offset)[0] == encoded_dans
    ]
    if len(starts) != 1:
        raise ValueError(f"{label} must contain exactly one encoded DanS marker; found {len(starts)}")
    start = starts[0]
    decoded = [
        struct.unpack_from("<I", data, offset)[0] ^ xor_key
        for offset in range(start, rich_offset, 4)
    ]
    if len(decoded) < 6 or decoded[:4] != [0x536E6144, 0, 0, 0] or (len(decoded) - 4) % 2:
        raise ValueError(f"{label} Rich header structure is unexpected")
    return start, rich_offset + 8


def _pe_checksum(data: bytes, checksum_offset: int) -> int:
    if checksum_offset < 0 or checksum_offset + 4 > len(data):
        raise ValueError("PE checksum field is outside the file")
    total = 0
    for offset in range(0, len(data), 2):
        if checksum_offset <= offset < checksum_offset + 4:
            word = 0
        elif offset + 1 < len(data):
            word = data[offset] | (data[offset + 1] << 8)
        else:
            word = data[offset]
        total = (total + word) & 0xFFFFFFFF
        total = (total & 0xFFFF) + (total >> 16)
    total = (total & 0xFFFF) + (total >> 16)
    return (total + len(data)) & 0xFFFFFFFF


def sanitize_rich_header(path: Path) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"frozen executable is missing: {path}")
    source = path.read_bytes()
    layout = _pe_layout(source, path.name)
    if layout["certificate"] != (0, 0):
        raise ValueError("refusing to rewrite an Authenticode-signed executable")
    rich_start, rich_end = _rich_header_span(source, path.name)
    expected = source[:rich_start] + b"\0" * (rich_end - rich_start) + source[rich_end:]
    checksum_offset = layout["checksum_offset"]
    checksum = _pe_checksum(expected, checksum_offset)
    expected = (
        expected[:checksum_offset]
        + struct.pack("<I", checksum)
        + expected[checksum_offset + 4:]
    )
    _pe_layout(expected, path.name)
    if any(expected[rich_start:rich_end]) or b"Rich" in expected[0x40:layout["pe_offset"]]:
        raise ValueError("Rich-header neutralisation did not remove the compiler record")
    approved = [(rich_start, rich_end), (checksum_offset, checksum_offset + 4)]
    if len(expected) != len(source) or _outside_spans(expected, approved) != _outside_spans(source, approved):
        raise ValueError("frozen executable changed outside Rich-header and checksum spans")
    _atomic_write(path, expected)
    after = path.read_bytes()
    if after != expected or _pe_checksum(after, checksum_offset) != checksum:
        raise ValueError("frozen executable Rich-header rewrite did not verify")
    print(
        f"Neutralised Rich header ({rich_end - rich_start} bytes) and recomputed PE checksum: {path.name}"
    )


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
    sanitize_rich = subparsers.add_parser("sanitize-rich-header")
    sanitize_rich.add_argument("path", type=Path)
    args = parser.parse_args()
    if args.command == "prepare-templates":
        prepare_templates(args.source_dir.resolve(), args.copied_dir.resolve(), args.version)
    elif args.command == "verify-setup":
        verify_setup(args.path.resolve(), args.version)
    elif args.command == "verify-uninstaller":
        verify_uninstaller(args.path.resolve(), args.version)
    else:
        sanitize_rich_header(args.path.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
