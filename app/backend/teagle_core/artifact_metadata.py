"""Deterministic metadata minimisation for generated scientific artifacts.

TEagle retains one human-attribution value only: ``Tuna Birgun``. Office Open
XML package edits are namespace-aware and atomic so metadata removal cannot
leave dangling package relationships or corrupt document XML.
"""
from __future__ import annotations

import os
import posixpath
import struct
import tempfile
import zlib
from collections import Counter
from pathlib import Path
from urllib.parse import unquote, urlsplit
from xml.dom import Node, minidom
from xml.parsers.expat import ExpatError
from zipfile import ZIP_DEFLATED, ZipFile, ZipInfo


AUTHOR = "Tuna Birgun"
_FIXED_ZIP_TIME = (1980, 1, 1, 0, 0, 0)
_PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"
_PNG_CRITICAL_CHUNKS = {b"IHDR", b"PLTE", b"IDAT", b"IEND"}
# These ancillary chunks affect faithful display or physical sizing rather than authorship or history.
# Everything else is optional metadata and is removed instead of attempting to interpret arbitrary payloads.
_PNG_RENDERING_CHUNKS = {
    b"cHRM", b"gAMA", b"iCCP", b"sBIT", b"sRGB", b"cICP", b"mDCV", b"cLLI",
    b"bKGD", b"hIST", b"pHYs", b"sPLT", b"tRNS",
}
_PNG_CANONICAL_AUTHOR = b"Author\0" + AUTHOR.encode("latin-1")

_CONTENT_TYPES_PART = "[Content_Types].xml"
_ROOT_RELS_PART = "_rels/.rels"
_CORE_PART = "docProps/core.xml"
_CORE_PART_KEY = _CORE_PART.casefold()
_ALWAYS_DROP_PART_KEYS = {
    "docprops/app.xml",
    "docprops/custom.xml",
}

_CONTENT_TYPES_NS = "http://schemas.openxmlformats.org/package/2006/content-types"
_RELATIONSHIPS_NS = "http://schemas.openxmlformats.org/package/2006/relationships"
_CORE_PROPERTIES_NS = (
    "http://schemas.openxmlformats.org/package/2006/metadata/core-properties"
)
_DC_NS = "http://purl.org/dc/elements/1.1/"
_CORE_CONTENT_TYPE = (
    "application/vnd.openxmlformats-package.core-properties+xml"
)
_CORE_REL_TYPE = (
    "http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties"
)
_OPTIONAL_METADATA_REL_TYPES = {
    "http://schemas.openxmlformats.org/officeDocument/2006/relationships/extended-properties",
    "http://schemas.openxmlformats.org/officeDocument/2006/relationships/custom-properties",
    "http://schemas.openxmlformats.org/package/2006/relationships/metadata/thumbnail",
}

_WORD_NS = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
_WORD_2010_NS = "http://schemas.microsoft.com/office/word/2010/wordml"
_WORD_2012_NS = "http://schemas.microsoft.com/office/word/2012/wordml"
_WORD_2016_CID_NS = (
    "http://schemas.microsoft.com/office/word/2016/wordml/cid"
)
_WORD_2018_DU_NS = (
    "http://schemas.microsoft.com/office/word/2018/wordml/cex"
)

_SVG_NS = "http://www.w3.org/2000/svg"
_RDF_NS = "http://www.w3.org/1999/02/22-rdf-syntax-ns#"
_XMLNS_NS = "http://www.w3.org/2000/xmlns/"

_CORE_XML = (
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
    '<cp:coreProperties '
    f'xmlns:cp="{_CORE_PROPERTIES_NS}" '
    f'xmlns:dc="{_DC_NS}">'
    f'<dc:creator>{AUTHOR}</dc:creator>'
    '</cp:coreProperties>'
).encode("utf-8")


def _local_name(node: Node) -> str:
    return node.localName or node.nodeName.rsplit(":", 1)[-1]


def _element_children(node: Node) -> list[Node]:
    return [child for child in node.childNodes if child.nodeType == Node.ELEMENT_NODE]


def _parse_xml(data: bytes, part_name: str) -> minidom.Document:
    try:
        return minidom.parseString(data)
    except (ExpatError, ValueError) as exc:
        raise ValueError(f"invalid XML in OOXML part {part_name}: {exc}") from exc


def _xml_bytes(document: minidom.Document) -> bytes:
    return document.toxml(encoding="utf-8")


def _normalise_part_name(value: str, *, base: str = "") -> str | None:
    """Resolve an OPC part URI to a package-relative, decoded POSIX path."""
    parsed = urlsplit(value)
    if parsed.scheme or parsed.netloc:
        return None
    path = unquote(parsed.path).replace("\\", "/")
    if not path:
        return None
    combined = path.lstrip("/") if path.startswith("/") else posixpath.join(base, path)
    normalised = posixpath.normpath(combined).lstrip("/")
    if normalised in {"", "."} or normalised == ".." or normalised.startswith("../"):
        return None
    return normalised


def _relationship_source_base(rels_part: str) -> str | None:
    if rels_part == _ROOT_RELS_PART:
        return ""
    marker = "/_rels/"
    if marker not in rels_part or not rels_part.endswith(".rels"):
        return None
    prefix, rel_name = rels_part.rsplit(marker, 1)
    source_part = posixpath.join(prefix, rel_name[:-5])
    return posixpath.dirname(source_part)


def _relationship_elements(document: minidom.Document, part_name: str) -> list[Node]:
    root = document.documentElement
    if _local_name(root) != "Relationships" or root.namespaceURI != _RELATIONSHIPS_NS:
        raise ValueError(f"invalid Relationships root in {part_name}")
    elements = _element_children(root)
    if any(
        _local_name(element) != "Relationship"
        or element.namespaceURI != _RELATIONSHIPS_NS
        for element in elements
    ):
        raise ValueError(f"unexpected element in Relationships part {part_name}")
    return elements


def _relationship_target_part(rels_part: str, relationship: Node) -> str | None:
    if relationship.getAttribute("TargetMode").casefold() == "external":
        return None
    base = _relationship_source_base(rels_part)
    if base is None:
        return None
    return _normalise_part_name(relationship.getAttribute("Target"), base=base)


def _is_thumbnail_part(part_name: str) -> bool:
    key = part_name.casefold()
    return key.startswith("docprops/thumbnail.")


def _metadata_parts_to_drop(entries: dict[str, bytes]) -> set[str]:
    drop = {
        name
        for name in entries
        if name.casefold() in _ALWAYS_DROP_PART_KEYS or _is_thumbnail_part(name)
    }
    for name, data in entries.items():
        if not name.endswith(".rels"):
            continue
        document = _parse_xml(data, name)
        for relationship in _relationship_elements(document, name):
            if relationship.getAttribute("Type") not in _OPTIONAL_METADATA_REL_TYPES:
                continue
            target = _relationship_target_part(name, relationship)
            if target is not None:
                drop.add(target)
    return drop


def _qualified_child_name(root: Node, local_name: str) -> str:
    return f"{root.prefix}:{local_name}" if root.prefix else local_name


def _clean_content_types_xml(data: bytes, drop_parts: set[str]) -> bytes:
    document = _parse_xml(data, _CONTENT_TYPES_PART)
    root = document.documentElement
    if _local_name(root) != "Types" or root.namespaceURI != _CONTENT_TYPES_NS:
        raise ValueError("invalid content-types root")
    drop_keys = {part.casefold() for part in drop_parts}
    for element in list(_element_children(root)):
        if _local_name(element) != "Override" or element.namespaceURI != _CONTENT_TYPES_NS:
            continue
        target = _normalise_part_name(element.getAttribute("PartName"))
        target_key = target.casefold() if target else ""
        if target_key in drop_keys or target_key == _CORE_PART_KEY:
            root.removeChild(element)
            element.unlink()
    override = document.createElementNS(
        _CONTENT_TYPES_NS,
        _qualified_child_name(root, "Override"),
    )
    override.setAttribute("PartName", f"/{_CORE_PART}")
    override.setAttribute("ContentType", _CORE_CONTENT_TYPE)
    root.appendChild(override)
    return _xml_bytes(document)


def _next_relationship_id(used: set[str]) -> str:
    number = 1
    while f"rId{number}" in used:
        number += 1
    return f"rId{number}"


def _clean_relationships_xml(
    part_name: str,
    data: bytes,
    drop_parts: set[str],
) -> bytes:
    document = _parse_xml(data, part_name)
    root = document.documentElement
    relationships = _relationship_elements(document, part_name)
    drop_keys = {part.casefold() for part in drop_parts}
    for relationship in list(relationships):
        rel_type = relationship.getAttribute("Type")
        target = _relationship_target_part(part_name, relationship)
        target_key = target.casefold() if target else ""
        should_drop = (
            rel_type in _OPTIONAL_METADATA_REL_TYPES
            or target_key in drop_keys
            or (part_name == _ROOT_RELS_PART and rel_type == _CORE_REL_TYPE)
            or (part_name == _ROOT_RELS_PART and target_key == _CORE_PART_KEY)
        )
        if should_drop:
            root.removeChild(relationship)
            relationship.unlink()

    if part_name == _ROOT_RELS_PART:
        # Root relationship IDs are not referenced from another package part,
        # so malformed duplicate or missing IDs can be repaired without
        # changing any document-level r:id reference.
        used: set[str] = set()
        for relationship in _relationship_elements(document, part_name):
            rel_id = relationship.getAttribute("Id")
            if not rel_id or rel_id in used:
                rel_id = _next_relationship_id(used)
                relationship.setAttribute("Id", rel_id)
            used.add(rel_id)
        core = document.createElementNS(
            _RELATIONSHIPS_NS,
            _qualified_child_name(root, "Relationship"),
        )
        core.setAttribute("Id", _next_relationship_id(used))
        core.setAttribute("Type", _CORE_REL_TYPE)
        core.setAttribute("Target", _CORE_PART)
        root.appendChild(core)
    return _xml_bytes(document)


def _word_metadata_attribute_kind(attribute: Node) -> str | None:
    namespace = attribute.namespaceURI
    local = attribute.localName or attribute.nodeName.rsplit(":", 1)[-1]
    if namespace == _WORD_NS:
        if local.startswith("rsid"):
            return "remove"
        if local in {"date", "initials"}:
            return "remove"
        if local == "author":
            return "author"
    if namespace == _WORD_2010_NS and local in {"paraId", "textId"}:
        return "remove"
    if namespace == _WORD_2012_NS and local == "docId":
        return "remove"
    is_microsoft_word_namespace = bool(
        namespace
        and namespace.startswith("http://schemas.microsoft.com/office/word/")
    )
    if (
        (namespace == _WORD_2016_CID_NS or is_microsoft_word_namespace)
        and local == "durableId"
    ):
        return "remove"
    if (
        (namespace == _WORD_2018_DU_NS or is_microsoft_word_namespace)
        and local == "dateUtc"
    ):
        return "remove"
    return None


def _all_elements(document: minidom.Document) -> list[Node]:
    return list(document.getElementsByTagName("*"))


def _clean_word_xml(part_name: str, data: bytes) -> bytes:
    document = _parse_xml(data, part_name)
    for element in _all_elements(document):
        if (
            element.namespaceURI == _WORD_NS
            and _local_name(element) in {"rsid", "rsids"}
        ):
            if element.parentNode is not None:
                element.parentNode.removeChild(element)
                element.unlink()
            continue
        attributes = [
            element.attributes.item(index)
            for index in range(element.attributes.length)
        ]
        for attribute in attributes:
            kind = _word_metadata_attribute_kind(attribute)
            if kind == "remove":
                element.removeAttributeNode(attribute)
            elif kind == "author":
                element.setAttributeNS(attribute.namespaceURI, attribute.nodeName, AUTHOR)
    return _xml_bytes(document)


def _clean_package_xml(name: str, data: bytes, drop_parts: set[str]) -> bytes:
    key = name.casefold()
    if key == _CORE_PART_KEY:
        return _CORE_XML
    if name == _CONTENT_TYPES_PART:
        return _clean_content_types_xml(data, drop_parts)
    if name.endswith(".rels"):
        return _clean_relationships_xml(name, data, drop_parts)
    if key.startswith("word/") and key.endswith(".xml"):
        return _clean_word_xml(name, data)
    return data


def _copy_zip_info(old: ZipInfo) -> ZipInfo:
    info = ZipInfo(old.filename, date_time=_FIXED_ZIP_TIME)
    info.compress_type = old.compress_type
    info.create_system = 0
    info.external_attr = old.external_attr
    info.internal_attr = old.internal_attr
    info.flag_bits = 0
    info.comment = b""
    info.extra = b""
    return info


def _replace_file(source: str | os.PathLike[str], target: str | os.PathLike[str]) -> None:
    os.replace(source, target)


def clean_ooxml_metadata(path: str | os.PathLike[str]) -> None:
    """Rewrite a generated DOCX/XLSX package with minimal deterministic metadata."""
    target = Path(path)
    with ZipFile(target, "r") as source:
        infos = source.infolist()
        names = [info.filename for info in infos]
        duplicates = sorted(
            name for name, count in Counter(names).items() if count > 1
        )
        if duplicates:
            raise ValueError(f"duplicate OOXML package parts: {', '.join(duplicates)}")
        entries = {info.filename: source.read(info) for info in infos}

    if _CONTENT_TYPES_PART not in entries or _ROOT_RELS_PART not in entries:
        raise ValueError("OOXML package lacks content types or root relationships")
    drop_parts = _metadata_parts_to_drop(entries)
    drop_keys = {part.casefold() for part in drop_parts}

    fd, tmp_name = tempfile.mkstemp(
        prefix=target.name + ".", suffix=".tmp", dir=target.parent
    )
    os.close(fd)
    try:
        with ZipFile(tmp_name, "w", compression=ZIP_DEFLATED) as destination:
            for old in infos:
                key = old.filename.casefold()
                if key in drop_keys or (key == _CORE_PART_KEY and old.filename != _CORE_PART):
                    continue
                destination.writestr(
                    _copy_zip_info(old),
                    _clean_package_xml(old.filename, entries[old.filename], drop_parts),
                )
            if _CORE_PART not in entries:
                info = ZipInfo(_CORE_PART, date_time=_FIXED_ZIP_TIME)
                info.compress_type = ZIP_DEFLATED
                info.create_system = 0
                destination.writestr(info, _CORE_XML)

        problems = ooxml_metadata_violations(tmp_name)
        if problems:
            raise ValueError("cleaned OOXML package failed validation: " + "; ".join(problems))
        _replace_file(tmp_name, target)
    finally:
        if os.path.exists(tmp_name):
            os.unlink(tmp_name)


def _core_metadata_problems(data: bytes) -> list[str]:
    problems: list[str] = []
    try:
        document = _parse_xml(data, _CORE_PART)
    except ValueError as exc:
        return [str(exc)]
    root = document.documentElement
    if _local_name(root) != "coreProperties" or root.namespaceURI != _CORE_PROPERTIES_NS:
        return ["invalid core metadata root"]
    children = _element_children(root)
    creators = [
        child
        for child in children
        if child.namespaceURI == _DC_NS and _local_name(child) == "creator"
    ]
    if len(creators) != 1:
        problems.append("author is not exactly Tuna Birgun")
    else:
        creator = creators[0]
        creator_text = "".join(
            child.data
            for child in creator.childNodes
            if child.nodeType in {Node.TEXT_NODE, Node.CDATA_SECTION_NODE}
        )
        creator_has_elements = bool(_element_children(creator))
        creator_attributes = [
            creator.attributes.item(index)
            for index in range(creator.attributes.length)
            if creator.attributes.item(index).namespaceURI != _XMLNS_NS
        ]
        if creator_text != AUTHOR or creator_has_elements or creator_attributes:
            problems.append("author is not exactly Tuna Birgun")
    if len(children) != 1 or any(child not in creators for child in children):
        problems.append("non-author core metadata retained")
    if any(
        child.nodeType in {Node.TEXT_NODE, Node.CDATA_SECTION_NODE} and child.data.strip()
        for child in root.childNodes
        if child not in creators
    ):
        problems.append("non-author core metadata retained")
    root_attributes = [
        root.attributes.item(index)
        for index in range(root.attributes.length)
        if root.attributes.item(index).namespaceURI != _XMLNS_NS
        and not root.attributes.item(index).nodeName.startswith("xmlns")
    ]
    if root_attributes:
        problems.append("non-author core metadata retained")
    return problems


def _content_types_problems(data: bytes, drop_parts: set[str]) -> list[str]:
    try:
        document = _parse_xml(data, _CONTENT_TYPES_PART)
    except ValueError as exc:
        return [str(exc)]
    root = document.documentElement
    if _local_name(root) != "Types" or root.namespaceURI != _CONTENT_TYPES_NS:
        return ["invalid content-types root"]
    drop_keys = {part.casefold() for part in drop_parts}
    core_overrides: list[Node] = []
    problems: list[str] = []
    for element in _element_children(root):
        if _local_name(element) != "Override" or element.namespaceURI != _CONTENT_TYPES_NS:
            continue
        target = _normalise_part_name(element.getAttribute("PartName"))
        target_key = target.casefold() if target else ""
        if target_key in drop_keys:
            problems.append(f"content-type reference retained: {target}")
        if target_key == _CORE_PART_KEY:
            core_overrides.append(element)
    if len(core_overrides) != 1:
        problems.append("core content-type override is not unique")
    elif core_overrides[0].getAttribute("ContentType") != _CORE_CONTENT_TYPE:
        problems.append("core content type is incorrect")
    return problems


def _word_metadata_problems(part_name: str, data: bytes) -> list[str]:
    try:
        document = _parse_xml(data, part_name)
    except ValueError as exc:
        return [str(exc)]
    for element in _all_elements(document):
        if (
            element.namespaceURI == _WORD_NS
            and _local_name(element) in {"rsid", "rsids"}
        ):
            return [f"Word revision identifier or author metadata retained: {part_name}"]
        for index in range(element.attributes.length):
            attribute = element.attributes.item(index)
            kind = _word_metadata_attribute_kind(attribute)
            if kind == "remove" or (kind == "author" and attribute.value != AUTHOR):
                return [f"Word revision identifier or author metadata retained: {part_name}"]
    return []


def ooxml_metadata_violations(path: str | os.PathLike[str]) -> list[str]:
    """Return semantic metadata-policy and package-integrity violations."""
    problems: list[str] = []
    try:
        with ZipFile(path, "r") as package:
            infos = package.infolist()
            names = [info.filename for info in infos]
            entries = {info.filename: package.read(info) for info in infos}
            if package.comment:
                problems.append("ZIP archive comment retained")
    except Exception as exc:
        return [f"invalid OOXML package: {exc}"]

    duplicate_names = sorted(
        name for name, count in Counter(names).items() if count > 1
    )
    problems.extend(f"duplicate OOXML package part: {name}" for name in duplicate_names)
    name_keys = {name.casefold(): name for name in names}
    drop_parts = {
        name
        for name in names
        if name.casefold() in _ALWAYS_DROP_PART_KEYS or _is_thumbnail_part(name)
    }

    relationship_documents: dict[str, minidom.Document] = {}
    parsed_relationships: dict[str, list[Node]] = {}
    for name, data in entries.items():
        if not name.endswith(".rels"):
            continue
        try:
            document = _parse_xml(data, name)
            relationships = _relationship_elements(document, name)
            relationship_documents[name] = document
            parsed_relationships[name] = relationships
        except ValueError as exc:
            problems.append(str(exc))
            continue
        seen_ids: set[str] = set()
        for relationship in relationships:
            rel_id = relationship.getAttribute("Id")
            if not rel_id or rel_id in seen_ids:
                problems.append(f"missing or duplicate relationship ID: {name}")
            seen_ids.add(rel_id)
            rel_type = relationship.getAttribute("Type")
            target = _relationship_target_part(name, relationship)
            if rel_type in _OPTIONAL_METADATA_REL_TYPES:
                problems.append(f"optional metadata relationship retained: {name}")
                if target is not None:
                    drop_parts.add(target)
            if (
                target is not None
                and relationship.getAttribute("TargetMode").casefold() != "external"
                and target.casefold() not in name_keys
            ):
                problems.append(f"dangling package relationship: {name} -> {target}")

    drop_keys = {part.casefold() for part in drop_parts}
    for name, relationships in parsed_relationships.items():
        for relationship in relationships:
            target = _relationship_target_part(name, relationship)
            if target is not None and target.casefold() in drop_keys:
                problems.append(f"relationship to optional metadata retained: {name}")

    for part in sorted(drop_parts):
        if part.casefold() in name_keys:
            problems.append(f"optional metadata part retained: {name_keys[part.casefold()]}")

    core_name = name_keys.get(_CORE_PART_KEY)
    if core_name is None:
        problems.append("core metadata part missing")
    else:
        problems.extend(_core_metadata_problems(entries[core_name]))
        if core_name != _CORE_PART:
            problems.append("core metadata part uses non-canonical name")

    if _CONTENT_TYPES_PART not in entries:
        problems.append("content-types part missing")
    else:
        problems.extend(_content_types_problems(entries[_CONTENT_TYPES_PART], drop_parts))

    root_relationships = relationship_documents.get(_ROOT_RELS_PART)
    if root_relationships is None:
        problems.append("root relationships part missing or invalid")
    else:
        core_relationships = [
            relationship
            for relationship in _relationship_elements(root_relationships, _ROOT_RELS_PART)
            if relationship.getAttribute("Type") == _CORE_REL_TYPE
        ]
        if len(core_relationships) != 1:
            problems.append("core relationship is not unique")
        else:
            relationship = core_relationships[0]
            target = _relationship_target_part(_ROOT_RELS_PART, relationship)
            if target is None or target.casefold() != _CORE_PART_KEY:
                problems.append("core relationship target is incorrect")

    for info in infos:
        if info.date_time != _FIXED_ZIP_TIME:
            problems.append(f"non-deterministic ZIP timestamp: {info.filename}")
        if info.comment:
            problems.append(f"ZIP entry comment retained: {info.filename}")
        if info.extra:
            problems.append(f"ZIP entry extra metadata retained: {info.filename}")
        key = info.filename.casefold()
        if key.startswith("word/") and key.endswith(".xml"):
            problems.extend(_word_metadata_problems(info.filename, entries[info.filename]))
    return sorted(set(problems))


def clean_pdf_metadata(path: str | os.PathLike[str]) -> None:
    """Retain document navigation and pages, but only the author metadata."""
    from pypdf import PdfReader, PdfWriter
    from pypdf.generic import NameObject

    target = Path(path)
    reader = PdfReader(str(target))
    writer = PdfWriter()
    # Clone the catalog so bookmarks, named destinations, page mode, forms and
    # initial-view settings survive. Remove the metadata reference and then
    # discard unreferenced objects so the old XMP stream is absent physically,
    # not merely hidden from readers.
    writer.clone_document_from_reader(reader)
    writer._info = None
    writer.add_metadata({"/Author": AUTHOR})
    writer.xmp_metadata = None
    writer._ID = None
    writer.root_object.pop(NameObject("/Metadata"), None)
    writer.compress_identical_objects(remove_duplicates=True, remove_unreferenced=True)
    fd, tmp_name = tempfile.mkstemp(prefix=target.name + ".", suffix=".tmp", dir=target.parent)
    os.close(fd)
    try:
        with open(tmp_name, "wb") as stream:
            writer.write(stream)
        os.replace(tmp_name, target)
    finally:
        if os.path.exists(tmp_name):
            os.unlink(tmp_name)


def _png_chunk(chunk_type: bytes, payload: bytes) -> bytes:
    if len(chunk_type) != 4:
        raise ValueError("PNG chunk type must be four bytes")
    checksum = zlib.crc32(chunk_type + payload) & 0xFFFFFFFF
    return struct.pack(">I", len(payload)) + chunk_type + payload + struct.pack(">I", checksum)


def _png_chunks(payload: bytes) -> list[tuple[bytes, bytes]]:
    """Parse and CRC-check a complete PNG without an image-decoding dependency."""
    if not payload.startswith(_PNG_SIGNATURE):
        raise ValueError("missing PNG signature")
    chunks: list[tuple[bytes, bytes]] = []
    offset = len(_PNG_SIGNATURE)
    plte_seen = False
    idat_seen = False
    idat_ended = False
    while offset < len(payload):
        if len(payload) - offset < 12:
            raise ValueError("truncated PNG chunk")
        length = struct.unpack(">I", payload[offset:offset + 4])[0]
        if length > 0x7FFFFFFF:
            raise ValueError("PNG chunk length exceeds the specification limit")
        end = offset + 12 + length
        if end > len(payload):
            raise ValueError("PNG chunk length exceeds the file")
        chunk_type = payload[offset + 4:offset + 8]
        if len(chunk_type) != 4 or any(
            not (65 <= byte <= 90 or 97 <= byte <= 122) for byte in chunk_type
        ):
            raise ValueError("PNG chunk type contains a non-letter byte")
        if 97 <= chunk_type[2] <= 122:
            raise ValueError("PNG chunk type uses the reserved lowercase bit")
        if 65 <= chunk_type[0] <= 90 and chunk_type not in _PNG_CRITICAL_CHUNKS:
            label = chunk_type.decode("ascii")
            raise ValueError(f"unknown critical PNG chunk: {label}")
        data = payload[offset + 8:offset + 8 + length]
        expected_crc = struct.unpack(">I", payload[offset + 8 + length:end])[0]
        actual_crc = zlib.crc32(chunk_type + data) & 0xFFFFFFFF
        if expected_crc != actual_crc:
            label = chunk_type.decode("ascii", errors="replace")
            raise ValueError(f"invalid CRC for {label} chunk")
        if not chunks:
            if chunk_type != b"IHDR" or length != 13:
                raise ValueError("IHDR must be the first PNG chunk")
        elif chunk_type == b"IHDR":
            raise ValueError("PNG contains more than one IHDR chunk")
        if chunk_type == b"PLTE":
            if plte_seen:
                raise ValueError("PNG contains more than one PLTE chunk")
            if idat_seen:
                raise ValueError("PLTE appears after IDAT")
            plte_seen = True
        if chunk_type == b"IDAT":
            if idat_ended:
                raise ValueError("PNG IDAT chunks are not consecutive")
            idat_seen = True
        elif idat_seen:
            idat_ended = True
        chunks.append((chunk_type, data))
        offset = end
        if chunk_type == b"IEND":
            if offset != len(payload):
                raise ValueError("data follows the IEND chunk")
            break
    if not idat_seen:
        raise ValueError("PNG contains no IDAT chunk")
    if chunks[-1] != (b"IEND", b""):
        raise ValueError("IEND must be the final empty PNG chunk")
    return chunks


def _png_metadata_problems(chunks: list[tuple[bytes, bytes]]) -> list[str]:
    problems: list[str] = []
    author_positions = [
        index for index, item in enumerate(chunks)
        if item == (b"tEXt", _PNG_CANONICAL_AUTHOR)
    ]
    if author_positions != [1]:
        problems.append("PNG does not contain one canonical author chunk immediately after IHDR")
    for chunk_type, payload in chunks:
        if chunk_type in _PNG_CRITICAL_CHUNKS or chunk_type in _PNG_RENDERING_CHUNKS:
            continue
        if chunk_type == b"tEXt" and payload == _PNG_CANONICAL_AUTHOR:
            continue
        label = chunk_type.decode("ascii", errors="replace")
        problems.append(f"PNG optional metadata retained: {label}")
    return problems


def clean_png_metadata(path: str | os.PathLike[str]) -> None:
    """Preserve PNG image chunks and replace optional metadata with one author field."""
    target = Path(path)
    chunks = _png_chunks(target.read_bytes())
    cleaned = bytearray(_PNG_SIGNATURE)
    for index, (chunk_type, payload) in enumerate(chunks):
        if chunk_type not in _PNG_CRITICAL_CHUNKS and chunk_type not in _PNG_RENDERING_CHUNKS:
            continue
        cleaned.extend(_png_chunk(chunk_type, payload))
        if index == 0:
            cleaned.extend(_png_chunk(b"tEXt", _PNG_CANONICAL_AUTHOR))
    result = bytes(cleaned)
    result_chunks = _png_chunks(result)
    violations = _png_metadata_problems(result_chunks)
    if violations:
        raise ValueError("invalid cleaned PNG metadata: " + "; ".join(violations))
    _atomic_write_bytes(target, result, suffix=".png")


def _svg_metadata_elements(document: minidom.Document) -> list[Node]:
    root = document.documentElement
    svg_namespace = root.namespaceURI
    return [
        element
        for element in _all_elements(document)
        if _local_name(element) == "metadata" and element.namespaceURI == svg_namespace
    ]


def _svg_author_metadata(document: minidom.Document) -> Node:
    root = document.documentElement
    svg_namespace = root.namespaceURI
    metadata_name = _qualified_child_name(root, "metadata")
    metadata = document.createElementNS(svg_namespace, metadata_name)
    rdf = document.createElementNS(_RDF_NS, "rdf:RDF")
    rdf.setAttributeNS(_XMLNS_NS, "xmlns:rdf", _RDF_NS)
    rdf.setAttributeNS(_XMLNS_NS, "xmlns:dc", _DC_NS)
    description = document.createElementNS(_RDF_NS, "rdf:Description")
    creator = document.createElementNS(_DC_NS, "dc:creator")
    sequence = document.createElementNS(_RDF_NS, "rdf:Seq")
    item = document.createElementNS(_RDF_NS, "rdf:li")
    item.appendChild(document.createTextNode(AUTHOR))
    sequence.appendChild(item)
    creator.appendChild(sequence)
    description.appendChild(creator)
    rdf.appendChild(description)
    metadata.appendChild(rdf)
    return metadata


def _atomic_write_bytes(target: Path, payload: bytes, *, suffix: str) -> None:
    fd, tmp_name = tempfile.mkstemp(prefix=target.name + ".", suffix=suffix, dir=target.parent)
    try:
        with os.fdopen(fd, "wb") as stream:
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        _replace_file(tmp_name, target)
    finally:
        if os.path.exists(tmp_name):
            os.unlink(tmp_name)


def clean_svg_metadata(path: str | os.PathLike[str]) -> None:
    """Atomically replace SVG metadata with a minimal RDF/DC author record."""
    target = Path(path)
    original = target.read_bytes()
    try:
        document = minidom.parseString(original)
    except (ExpatError, ValueError) as exc:
        raise ValueError(f"invalid SVG: {target}: {exc}") from exc
    root = document.documentElement
    if _local_name(root).casefold() != "svg":
        raise ValueError(f"invalid SVG root: {target}")
    for metadata in _svg_metadata_elements(document):
        if metadata.parentNode is not None:
            metadata.parentNode.removeChild(metadata)
            metadata.unlink()
    # Pretty-printers commonly leave indentation-only text nodes directly under the SVG root. They have
    # no rendering meaning, but a one-space blank line makes generated assets fail `git diff --check`.
    for child in list(root.childNodes):
        if child.nodeType == Node.TEXT_NODE and not child.data.strip():
            root.removeChild(child)
            child.unlink()
    root.insertBefore(_svg_author_metadata(document), root.firstChild)
    _atomic_write_bytes(target, document.toxml(encoding="utf-8"), suffix=".svg")


def _non_namespace_attributes(element: Node) -> list[Node]:
    return [
        element.attributes.item(index)
        for index in range(element.attributes.length)
        if element.attributes.item(index).namespaceURI != _XMLNS_NS
        and not element.attributes.item(index).nodeName.startswith("xmlns")
    ]


def _single_child(element: Node, namespace: str, local_name: str) -> Node | None:
    children = _element_children(element)
    if len(children) != 1:
        return None
    child = children[0]
    if child.namespaceURI != namespace or _local_name(child) != local_name:
        return None
    return child


def _svg_metadata_is_author_only(document: minidom.Document) -> bool:
    root = document.documentElement
    metadata_elements = _svg_metadata_elements(document)
    if len(metadata_elements) != 1 or metadata_elements[0].parentNode is not root:
        return False
    metadata = metadata_elements[0]
    if _non_namespace_attributes(metadata):
        return False
    rdf = _single_child(metadata, _RDF_NS, "RDF")
    description = _single_child(rdf, _RDF_NS, "Description") if rdf else None
    creator = _single_child(description, _DC_NS, "creator") if description else None
    sequence = _single_child(creator, _RDF_NS, "Seq") if creator else None
    item = _single_child(sequence, _RDF_NS, "li") if sequence else None
    if not all((rdf, description, creator, sequence, item)):
        return False
    if any(
        _non_namespace_attributes(element)
        for element in (rdf, description, creator, sequence, item)
    ):
        return False
    if _element_children(item):
        return False
    item_text = "".join(
        child.data
        for child in item.childNodes
        if child.nodeType in {Node.TEXT_NODE, Node.CDATA_SECTION_NODE}
    )
    return item_text == AUTHOR


def media_metadata_violations(path: str | os.PathLike[str]) -> list[str]:
    """Return semantic PDF/PNG/SVG metadata-policy violations."""
    target = Path(path)
    suffix = target.suffix.lower()
    problems: list[str] = []
    if suffix == ".pdf":
        from pypdf import PdfReader

        reader = PdfReader(str(target))
        metadata = {
            str(key): str(value)
            for key, value in (reader.metadata or {}).items()
            if value not in (None, "")
        }
        if metadata != {"/Author": AUTHOR}:
            problems.append(f"PDF metadata is not author-only: {metadata}")
        if "/Metadata" in reader.root_object:
            problems.append("PDF XMP metadata retained")
        if reader.trailer.get("/ID") is not None:
            problems.append("PDF document ID retained")
    elif suffix == ".png":
        try:
            chunks = _png_chunks(target.read_bytes())
        except ValueError as exc:
            problems.append(f"invalid PNG: {exc}")
        else:
            problems.extend(_png_metadata_problems(chunks))
    elif suffix == ".svg":
        try:
            document = minidom.parseString(target.read_bytes())
        except (ExpatError, ValueError) as exc:
            problems.append(f"invalid SVG: {exc}")
        else:
            root = document.documentElement
            if _local_name(root).casefold() != "svg":
                problems.append("invalid SVG root")
            elif not _svg_metadata_is_author_only(document):
                problems.append("SVG metadata is not author-only RDF/DC")
    else:
        raise ValueError(f"unsupported media file: {target}")
    return sorted(set(problems))
