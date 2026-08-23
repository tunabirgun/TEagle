from __future__ import annotations

import struct


def manifest(identity: str = "JR.Inno.Setup", description: str = "Inno Setup", extra: str = "") -> bytes:
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\r\n'
        '<assembly xmlns="urn:schemas-microsoft-com:asm.v1" manifestVersion="1.0">\r\n'
        '<assemblyIdentity name="'
        + identity
        + '" processorArchitecture="x86" version="1.0.0.0" type="win32"/>\r\n'
        '<description>'
        + description
        + '</description>\r\n'
        + extra
        + '<trustInfo xmlns="urn:schemas-microsoft-com:asm.v3"><security><requestedPrivileges>'
        '<requestedExecutionLevel level="asInvoker" uiAccess="false"/>'
        '</requestedPrivileges></security></trustInfo>\r\n'
        '</assembly>\r\n'
    ).encode("utf-8")


def synthetic_pe(
    manifest_bytes: bytes, *, rich: bool = False, certificate: bool = False,
    comment_copies: int = 1,
) -> bytes:
    pe_offset = 0x100
    optional_size = 0xE0
    section_rva = 0x1000
    section_raw = 0x400
    manifest_raw = 0x600
    manifest_rva = section_rva + manifest_raw - section_raw
    resource_directory_size = manifest_raw - section_raw + len(manifest_bytes)
    comment = "This installation was built with Inno Setup.".encode("utf-16le") * comment_copies
    fixed = struct.pack(
        "<13I", 0xFEEF04BD, 0x00010000, (51 << 16) | 1054, 0,
        0, 0, 0x3F, 0, 0x00040004, 1, 0, 0, 0,
    )
    version_text = "51.1054.0.0".encode("utf-16le")
    suffix = fixed + version_text + comment
    raw_size = ((manifest_raw + len(manifest_bytes) + len(suffix) + 0x1FF) // 0x200) * 0x200
    data = bytearray(section_raw + raw_size)
    data[:2] = b"MZ"
    struct.pack_into("<I", data, 0x3C, pe_offset)
    data[pe_offset:pe_offset + 4] = b"PE\0\0"
    coff = pe_offset + 4
    struct.pack_into("<HHIIIHH", data, coff, 0x14C, 1, 0, 0, 0, optional_size, 0x102)
    optional = coff + 20
    struct.pack_into("<H", data, optional, 0x10B)
    struct.pack_into("<I", data, optional + 56, 0x2000)
    struct.pack_into("<I", data, optional + 60, section_raw)
    struct.pack_into("<I", data, optional + 92, 16)
    struct.pack_into("<II", data, optional + 96 + 2 * 8, section_rva, resource_directory_size)
    if certificate:
        struct.pack_into("<II", data, optional + 96 + 4 * 8, len(data) - 8, 8)
    section = optional + optional_size
    data[section:section + 8] = b".rsrc\0\0\0"
    struct.pack_into("<IIII", data, section + 8, raw_size, section_rva, raw_size, section_raw)
    struct.pack_into("<I", data, section + 36, 0x40000040)

    base = section_raw
    struct.pack_into("<IIHHHH", data, base, 0, 0, 0, 0, 0, 1)
    struct.pack_into("<II", data, base + 16, 24, 0x80000020)
    struct.pack_into("<IIHHHH", data, base + 0x20, 0, 0, 0, 0, 0, 1)
    struct.pack_into("<II", data, base + 0x30, 1, 0x80000040)
    struct.pack_into("<IIHHHH", data, base + 0x40, 0, 0, 0, 0, 0, 1)
    struct.pack_into("<II", data, base + 0x50, 1033, 0x60)
    struct.pack_into("<IIII", data, base + 0x60, manifest_rva, len(manifest_bytes), 1252, 0)
    data[manifest_raw:manifest_raw + len(manifest_bytes)] = manifest_bytes
    suffix_start = manifest_raw + len(manifest_bytes)
    data[suffix_start:suffix_start + len(suffix)] = suffix

    if rich:
        key = 0x13579BDF
        decoded = (0x536E6144, 0, 0, 0, 0x01020304, 7)
        rich_start = 0x80
        for index, value in enumerate(decoded):
            struct.pack_into("<I", data, rich_start + index * 4, value ^ key)
        rich_marker = rich_start + len(decoded) * 4
        data[rich_marker:rich_marker + 4] = b"Rich"
        struct.pack_into("<I", data, rich_marker + 4, key)
    return bytes(data)
