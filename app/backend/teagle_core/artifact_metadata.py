"""Deterministic metadata minimisation for generated scientific artifacts.

Project artifacts retain one human-attribution value only: ``Tuna Birgun``. Every cleaner and validator takes an
``author`` keyword (default: that value); ``author=None`` writes and requires no attribution at all, which is what
the app uses for files it writes for end users. Office Open XML package edits are namespace-aware and atomic so
metadata removal cannot leave dangling package relationships or corrupt document XML.
"""
from __future__ import annotations

import io
import os
import posixpath
import re
import struct
import tempfile
import warnings
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


def _png_author_chunk(author: str) -> bytes:
    return b"Author\0" + author.encode("latin-1")


_PNG_MAX_BYTES = 512 * 1024 * 1024
_PNG_MAX_CHUNK_BYTES = 256 * 1024 * 1024
_PNG_MAX_PIXELS = 100_000_000
_PNG_VALID_DEPTHS = {
    0: {1, 2, 4, 8, 16}, 2: {8, 16}, 3: {1, 2, 4, 8},
    4: {8, 16}, 6: {8, 16},
}
_PNG_AUTHENTICITY_CHUNKS = {b"caBX"}

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
_CUSTOM_XML_REL_TYPE = (
    "http://schemas.openxmlformats.org/officeDocument/2006/relationships/customXml"
)
_CUSTOM_XML_PROPS_REL_TYPE = (
    "http://schemas.openxmlformats.org/officeDocument/2006/relationships/customXmlProps"
)
_CUSTOM_XML_PROPS_CONTENT_TYPE = (
    "application/vnd.openxmlformats-officedocument.customXmlProperties+xml"
)
_OOXML_UNSAFE_PART_PATTERNS = (
    (re.compile(r"(^|/)vbaproject\.bin$", re.I), "macro project"),
    (re.compile(r"(^|/)vbadata\.xml$", re.I), "macro data"),
    (re.compile(r"^_xmlsignatures/", re.I), "package signature"),
    (re.compile(r"(^|/)origin\.sigs$", re.I), "package signature origin"),
    (re.compile(r"(^|/)activex/", re.I), "ActiveX part"),
    (re.compile(r"(^|/)controls?/", re.I), "control part"),
    (re.compile(r"(^|/)ctrlprops/", re.I), "control properties"),
    (re.compile(r"(^|/)embeddings/", re.I), "embedded object"),
)
_OOXML_UNSAFE_REL_SUFFIXES = (
    "/vbaproject", "/vbadata", "/digital-signature/origin",
    "/digital-signature/signature", "/activexcontrol", "/control",
    "/oleobject", "/package", "/attachedtemplate",
)
_OOXML_TRACKED_CHANGE_NAMES = {
    "ins", "del", "movefrom", "moveto", "movefromrangestart",
    "movefromrangeend", "movetorangestart", "movetorangeend",
    "cellins", "celldel", "cellmerge",
    "customxmlinsrangestart", "customxmlinsrangeend",
    "customxmldelrangestart", "customxmldelrangeend",
    "customxmlmovefromrangestart", "customxmlmovefromrangeend",
    "customxmlmovetorangestart", "customxmlmovetorangeend",
}
_OOXML_FORM_NAMES = {"sdt", "ffdata", "formfield", "legacydrawing", "legacydrawinghf"}

_WORD_NS = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
_WORD_2010_NS = "http://schemas.microsoft.com/office/word/2010/wordml"
_WORD_2012_NS = "http://schemas.microsoft.com/office/word/2012/wordml"
_WORD_2016_CID_NS = (
    "http://schemas.microsoft.com/office/word/2016/wordml/cid"
)
_WORD_2018_DU_NS = (
    "http://schemas.microsoft.com/office/word/2018/wordml/cex"
)
_BIBLIOGRAPHY_NS = (
    "http://schemas.openxmlformats.org/officeDocument/2006/bibliography"
)
_CUSTOM_XML_NS = (
    "http://schemas.openxmlformats.org/officeDocument/2006/customXml"
)
_UUID_PATTERN = re.compile(
    r"^\{[0-9A-Fa-f]{8}-[0-9A-Fa-f]{4}-[0-9A-Fa-f]{4}-[0-9A-Fa-f]{4}-[0-9A-Fa-f]{12}\}$"
)

_SVG_NS = "http://www.w3.org/2000/svg"
_RDF_NS = "http://www.w3.org/1999/02/22-rdf-syntax-ns#"
_XMLNS_NS = "http://www.w3.org/2000/xmlns/"
_XML_NS = "http://www.w3.org/XML/1998/namespace"
_XLINK_NS = "http://www.w3.org/1999/xlink"

def _core_xml(author: str | None) -> bytes:
    creator = f"<dc:creator>{author}</dc:creator>" if author else ""
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        f'<cp:coreProperties xmlns:cp="{_CORE_PROPERTIES_NS}" xmlns:dc="{_DC_NS}">'
        f'{creator}</cp:coreProperties>'
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


def _is_custom_xml_item(part_name: str) -> bool:
    key = part_name.casefold()
    return (
        key.startswith("customxml/")
        and "/_rels/" not in key
        and key.endswith(".xml")
        and posixpath.basename(key).startswith("item")
        and not posixpath.basename(key).startswith("itemprops")
    )


def _non_namespace_attribute_names(element: Node) -> set[str]:
    return {
        attribute.nodeName
        for index in range(element.attributes.length)
        if (attribute := element.attributes.item(index)).namespaceURI != _XMLNS_NS
        and not attribute.nodeName.startswith("xmlns")
    }


def _has_non_whitespace_content(element: Node) -> bool:
    return any(
        child.nodeType == Node.ELEMENT_NODE
        or (
            child.nodeType in {Node.TEXT_NODE, Node.CDATA_SECTION_NODE}
            and child.data.strip()
        )
        for child in element.childNodes
    )


def _empty_bibliography_residue_parts(
    entries: dict[str, bytes],
) -> tuple[set[str], list[str]]:
    """Identify only the complete, empty Word bibliography custom XML package."""
    names_by_key = {name.casefold(): name for name in entries}
    drop_parts: set[str] = set()
    problems: list[str] = []
    for source_name, source_data in entries.items():
        if not _is_custom_xml_item(source_name):
            continue
        try:
            source_document = _parse_xml(source_data, source_name)
        except ValueError as exc:
            problems.append(f"invalid custom XML part: {exc}")
            continue
        source = source_document.documentElement
        if (
            source.namespaceURI != _BIBLIOGRAPHY_NS
            or _local_name(source) != "Sources"
            or _has_non_whitespace_content(source)
        ):
            continue
        if _non_namespace_attribute_names(source) != {"SelectedStyle", "StyleName"}:
            problems.append(
                f"malformed empty Word bibliography customXml: {source_name}"
            )
            continue

        directory = posixpath.dirname(source_name)
        filename = posixpath.basename(source_name)
        props_name = f"itemProps{filename[4:]}"
        expected_props = posixpath.join(directory, props_name)
        rels_name = posixpath.join(directory, "_rels", f"{filename}.rels")
        rels_actual = names_by_key.get(rels_name.casefold())
        if rels_actual is None:
            problems.append(
                f"malformed empty Word bibliography customXml: missing relationships for {source_name}"
            )
            continue
        try:
            item_relationships = _relationship_elements(
                _parse_xml(entries[rels_actual], rels_actual), rels_actual
            )
        except ValueError as exc:
            problems.append(f"malformed empty Word bibliography customXml: {exc}")
            continue
        if len(item_relationships) != 1:
            problems.append(
                f"malformed empty Word bibliography customXml: unexpected relationships for {source_name}"
            )
            continue
        item_relationship = item_relationships[0]
        props_target = _relationship_target_part(rels_actual, item_relationship)
        if (
            item_relationship.getAttribute("Type") != _CUSTOM_XML_PROPS_REL_TYPE
            or props_target is None
            or props_target.casefold() != expected_props.casefold()
        ):
            problems.append(
                f"malformed empty Word bibliography customXml: invalid property relationship for {source_name}"
            )
            continue
        props_actual = names_by_key.get(expected_props.casefold())
        if props_actual is None:
            problems.append(
                f"malformed empty Word bibliography customXml: missing properties for {source_name}"
            )
            continue
        try:
            props_document = _parse_xml(entries[props_actual], props_actual)
        except ValueError as exc:
            problems.append(f"malformed empty Word bibliography customXml: {exc}")
            continue
        datastore = props_document.documentElement
        datastore_id = datastore.getAttributeNS(_CUSTOM_XML_NS, "itemID")
        schema_refs = _element_children(datastore)
        if (
            datastore.namespaceURI != _CUSTOM_XML_NS
            or _local_name(datastore) != "datastoreItem"
            or _non_namespace_attribute_names(datastore) != {"ds:itemID"}
            or not _UUID_PATTERN.fullmatch(datastore_id)
            or len(schema_refs) != 1
            or schema_refs[0].namespaceURI != _CUSTOM_XML_NS
            or _local_name(schema_refs[0]) != "schemaRefs"
            or _non_namespace_attribute_names(schema_refs[0])
        ):
            problems.append(
                f"malformed empty Word bibliography customXml: invalid datastore properties for {source_name}"
            )
            continue
        schema_references = _element_children(schema_refs[0])
        if (
            len(schema_references) != 1
            or schema_references[0].namespaceURI != _CUSTOM_XML_NS
            or _local_name(schema_references[0]) != "schemaRef"
            or _non_namespace_attribute_names(schema_references[0]) != {"ds:uri"}
            or schema_references[0].getAttributeNS(_CUSTOM_XML_NS, "uri")
            != _BIBLIOGRAPHY_NS
        ):
            problems.append(
                f"malformed empty Word bibliography customXml: invalid schema reference for {source_name}"
            )
            continue

        inbound: list[tuple[str, Node]] = []
        for relationship_part, relationship_data in entries.items():
            if not relationship_part.endswith(".rels") or relationship_part == rels_actual:
                continue
            try:
                relationships = _relationship_elements(
                    _parse_xml(relationship_data, relationship_part), relationship_part
                )
            except ValueError as exc:
                problems.append(f"malformed empty Word bibliography customXml: {exc}")
                inbound = []
                break
            inbound.extend(
                (relationship_part, relationship)
                for relationship in relationships
                if _relationship_target_part(relationship_part, relationship)
                and _relationship_target_part(relationship_part, relationship).casefold()
                == source_name.casefold()
            )
        if len(inbound) != 1 or inbound[0][1].getAttribute("Type") != _CUSTOM_XML_REL_TYPE:
            problems.append(
                f"malformed empty Word bibliography customXml: invalid inbound relationship for {source_name}"
            )
            continue

        try:
            content_types = _parse_xml(entries[_CONTENT_TYPES_PART], _CONTENT_TYPES_PART)
        except (KeyError, ValueError) as exc:
            problems.append(f"malformed empty Word bibliography customXml: {exc}")
            continue
        overrides = [
            element
            for element in _element_children(content_types.documentElement)
            if _local_name(element) == "Override"
            and element.namespaceURI == _CONTENT_TYPES_NS
            and _normalise_part_name(element.getAttribute("PartName"))
            and _normalise_part_name(element.getAttribute("PartName")).casefold()
            == expected_props.casefold()
        ]
        if (
            len(overrides) != 1
            or overrides[0].getAttribute("ContentType") != _CUSTOM_XML_PROPS_CONTENT_TYPE
        ):
            problems.append(
                f"malformed empty Word bibliography customXml: invalid content type for {source_name}"
            )
            continue
        drop_parts.update({source_name, rels_actual, props_actual})
    return drop_parts, problems


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


def _external_relationship_problem(part_name: str, relationship: Node) -> str | None:
    target = relationship.getAttribute("Target").strip()
    mode = relationship.getAttribute("TargetMode").casefold()
    rel_type = relationship.getAttribute("Type").casefold()
    parsed = urlsplit(target)
    looks_machine_local = (
        parsed.scheme.casefold() == "file"
        or target.startswith(("\\\\", "//"))
        or bool(re.match(r"^[A-Za-z]:[\\/]", target))
    )
    if mode == "external":
        safe_hyperlink = (
            rel_type.endswith("/hyperlink")
            and parsed.scheme.casefold() in {"http", "https", "mailto"}
        )
        if not safe_hyperlink:
            return f"external file relationship: {part_name} -> {target or '<empty>'}"
    elif looks_machine_local:
        return f"machine-local relationship target: {part_name} -> {target}"
    return None


def _ooxml_unsafe_features(entries: dict[str, bytes]) -> list[str]:
    """Find package features whose meaning or authenticity a rewrite could invalidate."""
    problems: list[str] = []
    for name, data in entries.items():
        key = name.casefold()
        for pattern, label in _OOXML_UNSAFE_PART_PATTERNS:
            if pattern.search(name):
                problems.append(f"OOXML {label} cannot be rewritten: {name}")
        if key == _CONTENT_TYPES_PART.casefold():
            lowered = data.lower()
            if b"macroenabled" in lowered or b"vbaproject" in lowered:
                problems.append("OOXML macro-enabled content type cannot be rewritten")
            if b"digital-signature" in lowered:
                problems.append("OOXML signature content type cannot be rewritten")
            if b"activex" in lowered or b"controlproperties" in lowered:
                problems.append("OOXML ActiveX/control content type cannot be rewritten")
        if name.endswith(".rels"):
            try:
                relationships = _relationship_elements(_parse_xml(data, name), name)
            except ValueError:
                continue
            for relationship in relationships:
                rel_type = relationship.getAttribute("Type").casefold()
                if any(rel_type.endswith(suffix) for suffix in _OOXML_UNSAFE_REL_SUFFIXES):
                    problems.append(
                        f"OOXML unsafe relationship cannot be rewritten: {name} -> {rel_type}"
                    )
                external_problem = _external_relationship_problem(name, relationship)
                if external_problem:
                    problems.append(external_problem)
        if key.startswith("word/") and key.endswith(".xml"):
            try:
                document = _parse_xml(data, name)
            except ValueError:
                continue
            for element in _all_elements(document):
                local = _local_name(element).casefold()
                if local in _OOXML_TRACKED_CHANGE_NAMES:
                    problems.append(f"OOXML tracked change cannot be rewritten: {name}:{local}")
                if local in _OOXML_FORM_NAMES:
                    problems.append(f"OOXML form/control cannot be rewritten: {name}:{local}")
    return sorted(set(problems))


def _clean_word_xml(part_name: str, data: bytes, author: str | None = AUTHOR) -> bytes:
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
            if kind == "remove" or (kind == "author" and not author):
                element.removeAttributeNode(attribute)
            elif kind == "author":
                element.setAttributeNS(attribute.namespaceURI, attribute.nodeName, author)
    return _xml_bytes(document)


def _clean_package_xml(name: str, data: bytes, drop_parts: set[str], author: str | None = AUTHOR) -> bytes:
    key = name.casefold()
    if key == _CORE_PART_KEY:
        return _core_xml(author)
    if name == _CONTENT_TYPES_PART:
        return _clean_content_types_xml(data, drop_parts)
    if name.endswith(".rels"):
        return _clean_relationships_xml(name, data, drop_parts)
    if key.startswith("word/") and key.endswith(".xml"):
        return _clean_word_xml(name, data, author)
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


def clean_ooxml_metadata(path: str | os.PathLike[str], *, author: str | None = AUTHOR) -> None:
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
    unsafe = _ooxml_unsafe_features(entries)
    if unsafe:
        raise ValueError("OOXML rewrite refused: " + "; ".join(unsafe))
    drop_parts = _metadata_parts_to_drop(entries)
    bibliography_drop_parts, bibliography_problems = _empty_bibliography_residue_parts(entries)
    if bibliography_problems:
        raise ValueError("; ".join(bibliography_problems))
    drop_parts.update(bibliography_drop_parts)
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
                    _clean_package_xml(old.filename, entries[old.filename], drop_parts, author),
                )
            if _CORE_PART not in entries:
                info = ZipInfo(_CORE_PART, date_time=_FIXED_ZIP_TIME)
                info.compress_type = ZIP_DEFLATED
                info.create_system = 0
                destination.writestr(info, _core_xml(author))

        problems = ooxml_metadata_violations(tmp_name, author=author)
        if problems:
            raise ValueError("cleaned OOXML package failed validation: " + "; ".join(problems))
        _replace_file(tmp_name, target)
    finally:
        if os.path.exists(tmp_name):
            os.unlink(tmp_name)


def _core_metadata_problems(data: bytes, author: str | None = AUTHOR) -> list[str]:
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
    if not author:
        if creators:
            problems.append("author metadata retained")
    elif len(creators) != 1:
        problems.append(f"author is not exactly {author}")
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
        if creator_text != author or creator_has_elements or creator_attributes:
            problems.append(f"author is not exactly {author}")
    if len(children) != (1 if author else 0) or any(child not in creators for child in children):
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


def _word_metadata_problems(part_name: str, data: bytes, author: str | None = AUTHOR) -> list[str]:
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
            if kind == "remove" or (kind == "author" and attribute.value != author):
                return [f"Word revision identifier or author metadata retained: {part_name}"]
    return []


def ooxml_metadata_violations(path: str | os.PathLike[str], *, author: str | None = AUTHOR) -> list[str]:
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
    problems.extend(_ooxml_unsafe_features(entries))
    name_keys = {name.casefold(): name for name in names}
    drop_parts = {
        name
        for name in names
        if name.casefold() in _ALWAYS_DROP_PART_KEYS or _is_thumbnail_part(name)
    }
    bibliography_drop_parts, bibliography_problems = _empty_bibliography_residue_parts(entries)
    drop_parts.update(bibliography_drop_parts)
    problems.extend(bibliography_problems)
    problems.extend(
        f"empty Word bibliography customXml residue retained: {part}"
        for part in sorted(bibliography_drop_parts)
        if part.casefold() in name_keys
    )

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
        problems.extend(_core_metadata_problems(entries[core_name], author))
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
            problems.extend(_word_metadata_problems(info.filename, entries[info.filename], author))
    return sorted(set(problems))


def _pdf_unsafe_features(path: Path, reader=None) -> list[str]:
    """Find authenticity- or functionality-bearing structures before PDF reconstruction."""
    from pypdf import PdfReader
    from pypdf.generic import ArrayObject, DictionaryObject, IndirectObject

    problems: list[str] = []
    raw = path.read_bytes()
    if len(re.findall(rb"(?m)^%%EOF[ \t]*(?:\r?\n|$)", raw)) != 1:
        problems.append("PDF incremental or ambiguous EOF history cannot be rewritten")
        return problems
    if reader is None:
        reader = PdfReader(str(path), strict=True)
    if reader.is_encrypted:
        problems.append("encrypted PDF cannot be rewritten")
        return problems
    root = reader.root_object
    if "/Perms" in root:
        problems.append("PDF /Perms or DocMDP signature policy cannot be rewritten")
    if "/AcroForm" in root:
        problems.append("PDF AcroForm cannot be rewritten")
    names = root.get("/Names")
    if names is not None:
        try:
            names = names.get_object()
            if "/EmbeddedFiles" in names:
                problems.append("PDF embedded-file name tree cannot be rewritten")
        except Exception as exc:
            problems.append(f"PDF name tree cannot be validated: {exc}")

    seen: set[tuple[int, int] | int] = set()

    def walk(value) -> None:
        if isinstance(value, IndirectObject):
            marker = (value.idnum, value.generation)
            if marker in seen:
                return
            seen.add(marker)
            try:
                value = value.get_object()
            except Exception as exc:
                problems.append(f"PDF indirect object cannot be validated: {exc}")
                return
        elif isinstance(value, (DictionaryObject, ArrayObject)):
            marker = id(value)
            if marker in seen:
                return
            seen.add(marker)
        if isinstance(value, DictionaryObject):
            field_type = str(value.get("/FT", ""))
            object_type = str(value.get("/Type", ""))
            subtype = str(value.get("/Subtype", ""))
            if field_type == "/Sig" or object_type == "/Sig" or "/ByteRange" in value:
                problems.append("PDF signature dictionary cannot be rewritten")
            if subtype == "/FileAttachment" or any(key in value for key in ("/EmbeddedFiles", "/EF", "/AF")):
                problems.append("PDF attachment or embedded file cannot be rewritten")
            if str(value.get("/TransformMethod", "")) == "/DocMDP":
                problems.append("PDF DocMDP transform cannot be rewritten")
            for child in value.values():
                walk(child)
        elif isinstance(value, ArrayObject):
            for child in value:
                walk(child)

    walk(root)
    return sorted(set(problems))


def clean_pdf_metadata(path: str | os.PathLike[str], *, author: str | None = AUTHOR) -> None:
    """Retain document navigation and pages, but only the author metadata (none when author is None)."""
    from pypdf import PdfReader, PdfWriter
    from pypdf.generic import NameObject

    target = Path(path)
    unsafe = _pdf_unsafe_features(target)
    if unsafe:
        raise ValueError("PDF rewrite refused: " + "; ".join(unsafe))
    reader = PdfReader(str(target), strict=True)
    writer = PdfWriter()
    # Clone the catalog so bookmarks, named destinations, page mode, forms and
    # initial-view settings survive. Remove the metadata reference and then
    # discard unreferenced objects so the old XMP stream is absent physically,
    # not merely hidden from readers.
    writer.clone_document_from_reader(reader)
    writer._info = None
    if author:
        writer.add_metadata({"/Author": author})
    writer.xmp_metadata = None
    writer._ID = None
    writer.root_object.pop(NameObject("/Metadata"), None)
    writer.compress_identical_objects(remove_duplicates=True, remove_unreferenced=True)
    fd, tmp_name = tempfile.mkstemp(prefix=target.name + ".", suffix=".tmp", dir=target.parent)
    os.close(fd)
    try:
        with open(tmp_name, "wb") as stream:
            writer.write(stream)
            stream.flush()
            os.fsync(stream.fileno())
        problems = media_metadata_violations(tmp_name, expected_format=".pdf", author=author)
        if problems:
            raise ValueError("cleaned PDF failed validation: " + "; ".join(problems))
        _replace_file(tmp_name, target)
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
    if len(payload) > _PNG_MAX_BYTES:
        raise ValueError("PNG exceeds the bounded decoder input limit")
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
        if length > _PNG_MAX_CHUNK_BYTES:
            raise ValueError("PNG chunk exceeds the bounded decoder limit")
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
    _validate_png_structure(chunks)
    return chunks


def _validate_png_structure(chunks: list[tuple[bytes, bytes]]) -> None:
    ihdr = chunks[0][1]
    width, height, depth, color, compression, filter_method, interlace = struct.unpack(
        ">IIBBBBB", ihdr
    )
    if not width or not height or width * height > _PNG_MAX_PIXELS:
        raise ValueError("PNG dimensions are zero or exceed the bounded pixel limit")
    if color not in _PNG_VALID_DEPTHS or depth not in _PNG_VALID_DEPTHS[color]:
        raise ValueError("invalid PNG bit-depth/color-type combination")
    if compression != 0 or filter_method != 0 or interlace not in {0, 1}:
        raise ValueError("unsupported PNG compression, filter, or interlace method")
    palettes = [payload for chunk_type, payload in chunks if chunk_type == b"PLTE"]
    if color == 3 and not palettes:
        raise ValueError("indexed PNG has no PLTE chunk")
    if color in {0, 4} and palettes:
        raise ValueError("grayscale PNG must not contain a PLTE chunk")
    if palettes:
        entries = len(palettes[0]) // 3
        if not len(palettes[0]) or len(palettes[0]) % 3 or entries > 256:
            raise ValueError("invalid PNG PLTE length")
        if color == 3 and entries > 2 ** depth:
            raise ValueError("PNG palette exceeds the indexed bit depth")
    for chunk_type, payload in chunks:
        if chunk_type == b"pHYs" and (
            len(payload) != 9 or payload[-1] not in {0, 1}
        ):
            raise ValueError("invalid PNG pHYs chunk")
        if chunk_type == b"iCCP":
            separator = payload.find(b"\0")
            if separator < 1 or separator > 79 or separator + 2 > len(payload):
                raise ValueError("invalid PNG iCCP chunk")
            if payload[separator + 1] != 0:
                raise ValueError("unsupported PNG iCCP compression method")
            try:
                profile = zlib.decompress(payload[separator + 2:])
            except zlib.error as exc:
                raise ValueError(f"invalid compressed PNG ICC profile: {exc}") from exc
            if not profile or len(profile) > 32 * 1024 * 1024:
                raise ValueError("PNG ICC profile is empty or exceeds the bounded limit")


def _strict_png_decode(payload: bytes, chunks: list[tuple[bytes, bytes]]) -> None:
    """Decode with the bundled Qt runtime used by the frozen application."""
    from PySide6.QtCore import QBuffer, QByteArray, QIODevice
    from PySide6.QtGui import QImageReader

    width, height = struct.unpack(">II", chunks[0][1][:8])
    buffer = QBuffer()
    buffer.setData(QByteArray(payload))
    if not buffer.open(QIODevice.OpenModeFlag.ReadOnly):
        raise ValueError("Qt PNG decoder buffer could not be opened")
    reader = QImageReader(buffer, b"png")
    reader.setDecideFormatFromContent(True)
    image = reader.read()
    if image.isNull():
        raise ValueError(f"Qt PNG decoder rejected the image: {reader.errorString()}")
    if image.width() != width or image.height() != height:
        raise ValueError("Qt PNG decoder dimensions disagree with IHDR")


def independent_pillow_png_decode(payload: bytes) -> None:
    """Run the development-only independent decoder check when Pillow is available."""
    from PIL import Image

    chunks = _png_chunks(payload)
    width, height = struct.unpack(">II", chunks[0][1][:8])
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("error", Image.DecompressionBombWarning)
            with Image.open(io.BytesIO(payload)) as candidate:
                if candidate.format != "PNG" or candidate.size != (width, height):
                    raise ValueError("Pillow PNG decoder disagrees with IHDR")
                candidate.load()
    except (OSError, SyntaxError, Image.DecompressionBombWarning) as exc:
        raise ValueError(f"Pillow PNG decoder rejected the image: {exc}") from exc


def png_decoder_evidence() -> dict[str, str]:
    """Return the strict decoder identities used by this runtime."""
    from PySide6.QtCore import qVersion

    evidence = {"primary": f"Qt {qVersion()} QImageReader PNG"}
    try:
        from PIL import __version__ as pillow_version
    except ImportError:
        evidence["independent"] = "Pillow unavailable"
    else:
        evidence["independent"] = f"Pillow {pillow_version}"
    return evidence


def _png_unsafe_features(chunks: list[tuple[bytes, bytes]]) -> list[str]:
    problems: list[str] = []
    chunk_types = {chunk_type for chunk_type, _ in chunks}
    if b"eXIf" in chunk_types:
        problems.append("PNG EXIF/orientation metadata cannot be safely rewritten")
    for chunk_type in sorted(chunk_types & _PNG_AUTHENTICITY_CHUNKS):
        problems.append(
            f"PNG authenticity-bearing {chunk_type.decode('ascii')} chunk cannot be rewritten"
        )
    return problems


def _png_metadata_problems(chunks: list[tuple[bytes, bytes]], author: str | None = AUTHOR) -> list[str]:
    problems: list[str] = []
    canonical = _png_author_chunk(author) if author else None
    author_positions = [
        index for index, item in enumerate(chunks)
        if canonical is not None and item == (b"tEXt", canonical)
    ]
    if author_positions != ([1] if author else []):
        problems.append("PNG does not contain one canonical author chunk immediately after IHDR" if author
                        else "PNG author metadata retained")
    for chunk_type, payload in chunks:
        if chunk_type in _PNG_CRITICAL_CHUNKS or chunk_type in _PNG_RENDERING_CHUNKS:
            continue
        if chunk_type == b"tEXt" and payload == canonical:
            continue
        label = chunk_type.decode("ascii", errors="replace")
        problems.append(f"PNG optional metadata retained: {label}")
    return problems


def clean_png_metadata(path: str | os.PathLike[str], *, author: str | None = AUTHOR) -> None:
    """Preserve PNG image chunks and replace optional metadata with one author field (none when author is None)."""
    target = Path(path)
    original = target.read_bytes()
    chunks = _png_chunks(original)
    unsafe = _png_unsafe_features(chunks)
    if unsafe:
        raise ValueError("PNG rewrite refused: " + "; ".join(unsafe))
    _strict_png_decode(original, chunks)
    cleaned = bytearray(_PNG_SIGNATURE)
    for index, (chunk_type, payload) in enumerate(chunks):
        if chunk_type not in _PNG_CRITICAL_CHUNKS and chunk_type not in _PNG_RENDERING_CHUNKS:
            continue
        cleaned.extend(_png_chunk(chunk_type, payload))
        if index == 0 and author:
            cleaned.extend(_png_chunk(b"tEXt", _png_author_chunk(author)))
    result = bytes(cleaned)
    result_chunks = _png_chunks(result)
    _strict_png_decode(result, result_chunks)
    violations = _png_metadata_problems(result_chunks, author)
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


def _xml_comments(node: Node) -> list[Node]:
    comments: list[Node] = []
    for child in node.childNodes:
        if child.nodeType == Node.COMMENT_NODE:
            comments.append(child)
        comments.extend(_xml_comments(child))
    return comments


def _svg_reference_is_external(value: str) -> bool:
    value = value.strip().strip('"\'')
    if not value or value.startswith("#") or value.casefold().startswith("data:"):
        return False
    return True


def _svg_unsafe_features(document: minidom.Document, raw: bytes) -> list[str]:
    problems: list[str] = []
    if re.search(br"<!\s*(?:DOCTYPE|ENTITY)\b", raw, flags=re.I):
        problems.append("SVG DTD or entity declaration cannot be rewritten")
    for element in _all_elements(document):
        local = _local_name(element).casefold()
        if local == "script":
            problems.append("SVG script cannot be rewritten")
        for index in range(element.attributes.length):
            attribute = element.attributes.item(index)
            attr_local = (attribute.localName or attribute.nodeName.rsplit(":", 1)[-1]).casefold()
            if attribute.namespaceURI == _XML_NS and attr_local == "base":
                problems.append("SVG xml:base cannot be rewritten")
            if attr_local in {"href", "src"} and _svg_reference_is_external(attribute.value):
                problems.append(f"SVG external or file reference cannot be rewritten: {attribute.value}")
            if attr_local == "style":
                for reference in re.findall(r"url\(([^)]+)\)", attribute.value, flags=re.I):
                    if _svg_reference_is_external(reference):
                        problems.append(f"SVG external style reference cannot be rewritten: {reference}")
                if re.search(r"@import\b", attribute.value, flags=re.I):
                    problems.append("SVG external style import cannot be rewritten")
        if local == "style":
            text = "".join(
                child.data for child in element.childNodes
                if child.nodeType in {Node.TEXT_NODE, Node.CDATA_SECTION_NODE}
            )
            for reference in re.findall(r"url\(([^)]+)\)", text, flags=re.I):
                if _svg_reference_is_external(reference):
                    problems.append(f"SVG external stylesheet reference cannot be rewritten: {reference}")
            if re.search(r"@import\b", text, flags=re.I):
                problems.append("SVG external stylesheet import cannot be rewritten")
    return sorted(set(problems))


def _svg_author_metadata(document: minidom.Document, author: str = AUTHOR) -> Node:
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
    item.appendChild(document.createTextNode(author))
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


def clean_svg_metadata(path: str | os.PathLike[str], *, author: str | None = AUTHOR) -> None:
    """Atomically replace SVG metadata with a minimal RDF/DC author record (no metadata when author is None)."""
    target = Path(path)
    original = target.read_bytes()
    try:
        document = minidom.parseString(original)
    except (ExpatError, ValueError) as exc:
        raise ValueError(f"invalid SVG: {target}: {exc}") from exc
    root = document.documentElement
    if _local_name(root).casefold() != "svg" or root.namespaceURI != _SVG_NS:
        raise ValueError(f"invalid SVG root: {target}")
    unsafe = _svg_unsafe_features(document, original)
    if unsafe:
        raise ValueError("SVG rewrite refused: " + "; ".join(unsafe))
    for comment in _xml_comments(document):
        if comment.parentNode is not None:
            comment.parentNode.removeChild(comment)
            comment.unlink()
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
    if author:
        root.insertBefore(_svg_author_metadata(document, author), root.firstChild)
    payload = document.toxml(encoding="utf-8")
    payload = b"\n".join(
        b"" if not line.strip(b" \t") else line
        for line in payload.split(b"\n")
    )
    try:
        candidate = minidom.parseString(payload)
    except (ExpatError, ValueError) as exc:
        raise ValueError(f"cleaned SVG is invalid: {exc}") from exc
    problems = _svg_unsafe_features(candidate, payload)
    if not _svg_metadata_is_author_only(candidate, author):
        problems.append("SVG metadata is not author-only RDF/DC" if author else "SVG metadata retained")
    if _xml_comments(candidate):
        problems.append("SVG XML comment retained")
    if problems:
        raise ValueError("cleaned SVG failed validation: " + "; ".join(sorted(set(problems))))
    _atomic_write_bytes(target, payload, suffix=".svg")


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


def _svg_metadata_is_author_only(document: minidom.Document, author: str | None = AUTHOR) -> bool:
    root = document.documentElement
    metadata_elements = _svg_metadata_elements(document)
    if not author:
        return not metadata_elements
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
    return item_text == author


def media_metadata_violations(
    path: str | os.PathLike[str], *, expected_format: str | None = None, author: str | None = AUTHOR,
) -> list[str]:
    """Return semantic PDF/PNG/SVG metadata-policy violations (author=None: no attribution allowed)."""
    target = Path(path)
    suffix = (expected_format or target.suffix).lower()
    problems: list[str] = []
    if suffix == ".pdf":
        from pypdf import PdfReader

        try:
            problems.extend(_pdf_unsafe_features(target))
            if problems:
                return sorted(set(problems))
            reader = PdfReader(str(target), strict=True)
        except Exception as exc:
            return [f"invalid PDF: {exc}"]
        metadata = {
            str(key): str(value)
            for key, value in (reader.metadata or {}).items()
            if value not in (None, "")
        }
        if metadata != ({"/Author": author} if author else {}):
            problems.append(f"PDF metadata is not author-only: {metadata}" if author
                            else f"PDF metadata retained: {metadata}")
        if "/Metadata" in reader.root_object:
            problems.append("PDF XMP metadata retained")
        if reader.trailer.get("/ID") is not None:
            problems.append("PDF document ID retained")
    elif suffix == ".png":
        try:
            raw = target.read_bytes()
            chunks = _png_chunks(raw)
            problems.extend(_png_unsafe_features(chunks))
            if problems:
                return sorted(set(problems))
            _strict_png_decode(raw, chunks)
        except ValueError as exc:
            problems.append(f"invalid PNG: {exc}")
        else:
            problems.extend(_png_metadata_problems(chunks, author))
    elif suffix == ".svg":
        raw = target.read_bytes()
        try:
            document = minidom.parseString(raw)
        except (ExpatError, ValueError) as exc:
            problems.append(f"invalid SVG: {exc}")
        else:
            root = document.documentElement
            if _local_name(root).casefold() != "svg" or root.namespaceURI != _SVG_NS:
                problems.append("invalid SVG root")
            else:
                problems.extend(_svg_unsafe_features(document, raw))
                if _xml_comments(document):
                    problems.append("SVG XML comment retained")
                if not _svg_metadata_is_author_only(document, author):
                    problems.append("SVG metadata is not author-only RDF/DC" if author else "SVG metadata retained")
    else:
        raise ValueError(f"unsupported media file: {target}")
    return sorted(set(problems))
