"""Run provenance manifest — packs exact tool versions, checksums, params, environment
with every result (see Deliverables/09). Values are REAL (captured at runtime)."""
from __future__ import annotations
import copy, sys, platform, json, hashlib, datetime, re
from urllib.parse import urlsplit
from . import __version__ as TEAGLE_VERSION
from .primers import PRIMER3_VERSION            # guarded: "unavailable" if primer3 fails to load
from .domains import PYHMMER_VERSION            # guarded likewise
from .gff3 import is_public_accession, is_public_record_id, public_record_identity

# wall-clock / machine-specific fields that must NOT enter the content-addressed seal
_SEAL_EXCLUDE_TOP = ("createdUtc", "environment")
# volatile / label-only input fields: retrieval time, and the raw FASTA-header id (differs between the
# NCBI and ENA-fallback fetch paths for the SAME accession+sequence). The sequence sha256 + the resolved
# accession/organism/taxid still seal the scientific identity.
_SEAL_EXCLUDE_INPUT = ("retrievedUtc", "id", "assemblyName", "displayLocus", "chromName",
                       "source", "endpoint", "sourceUrl")   # serving-DB labels: recorded, never sealed (NCBI vs ENA fallback
                                                            # must not change the content-addressed seal for identical bytes)

SCHEMA_VERSION = "0.2"
CANONICALIZATION = "teagle-json-utf8-v1"
_SOURCE_PROVIDERS = {
    "NCBI nuccore": "ncbi",
    "NCBI feature table (efetch rettype=ft)": "ncbi",
    "NCBI E-utilities (efetch seq_start/seq_stop)": "ncbi",
    "ENA (EMBL-EBI)": "ena",
}
_SOURCE_URLS = {
    "ncbi": ("www.ncbi.nlm.nih.gov", "/nuccore/{accession}"),
    "ena": ("www.ebi.ac.uk", "/ena/browser/view/{accession}"),
}
_FETCH_ENDPOINTS = {
    "NCBI nuccore": ("ncbi_nuccore", "ncbi_efetch",
                      "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"),
    "ENA (EMBL-EBI)": ("ena_fasta", "ena_fasta", None),
}
_COORD_SOURCE = "NCBI E-utilities (efetch seq_start/seq_stop)"
_COORD_ENDPOINT = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
_UTC_SECONDS = re.compile(r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(?:Z|\+00:00)$")
SOURCE_CONTEXT_VERSION = "1.0"


def _seal_source(manifest: dict) -> dict:
    seal_src = {k: copy.deepcopy(v) for k, v in manifest.items()
                if k not in _SEAL_EXCLUDE_TOP and k != "manifestSha256"}
    if isinstance(seal_src.get("input"), dict):
        seal_src["input"] = {k: v for k, v in seal_src["input"].items()
                             if k not in _SEAL_EXCLUDE_INPUT}
    return seal_src


def canonical_manifest_bytes(manifest: dict) -> bytes:
    """Return the versioned canonical UTF-8 bytes covered by ``manifestSha256``."""
    if not isinstance(manifest, dict):
        raise TypeError("manifest must be a JSON object")
    if manifest.get("schemaVersion") != SCHEMA_VERSION:
        raise ValueError(f"unsupported manifest schema: {manifest.get('schemaVersion')!r}")
    if manifest.get("canonicalization") != CANONICALIZATION:
        raise ValueError(f"unsupported manifest canonicalization: {manifest.get('canonicalization')!r}")
    return json.dumps(_seal_source(manifest), ensure_ascii=False, allow_nan=False,
                      sort_keys=True, separators=(",", ":")).encode("utf-8")


def verify_manifest(manifest: dict) -> bool:
    """Verify the declared contract and content-addressed seal without mutating the manifest."""
    try:
        claimed = manifest.get("manifestSha256")
        if not isinstance(claimed, str) or len(claimed) != 64:
            return False
        if not _semantic_public_fields_valid(manifest):
            return False
        observed = hashlib.sha256(canonical_manifest_bytes(manifest)).hexdigest()
        return observed == claimed.lower()
    except (AttributeError, TypeError, ValueError):
        return False


def manifest_json_bytes(manifest: dict) -> bytes:
    """Serialize one verified manifest as strict UTF-8 JSON and reverify the round trip."""
    if not verify_manifest(manifest):
        raise ValueError("manifest seal verification failed")
    encoded = json.dumps(
        manifest, ensure_ascii=False, allow_nan=False, sort_keys=True, indent=2,
    ).encode("utf-8")
    if not verify_manifest(json.loads(encoded)):
        raise ValueError("serialized manifest failed seal verification")
    return encoded


def _record_descriptors(input_records: list | None) -> list | None:
    if input_records is None:
        return None
    if not isinstance(input_records, list):
        raise ValueError("input_records must be an ordered list")
    records = []
    for rec in input_records:
        if not isinstance(rec, dict) or "id" not in rec:
            raise ValueError("each input record must contain an id")
        if "sequence" in rec:
            seq = rec["sequence"]
            if not isinstance(seq, str):
                raise ValueError("input record sequence must be text")
            length = len(seq)
            digest = hashlib.sha256(seq.encode("utf-8")).hexdigest()
        else:
            length, digest = rec.get("length"), str(rec.get("sha256", "")).lower()
            if isinstance(length, bool) or not isinstance(length, int) or length < 0:
                raise ValueError("input record length must be a non-negative integer")
            if len(digest) != 64 or any(ch not in "0123456789abcdef" for ch in digest):
                raise ValueError("input record sha256 must be a 64-character hexadecimal SHA-256")
        identity_source = rec.get("identitySource", rec.get("identityTrust", "local_input"))
        if identity_source not in {"local_input", "validated_fetch", "generated"}:
            raise ValueError("input record identity source is invalid")
        if "idSha256" in rec:
            identity = {"id": rec["id"], "idSha256": rec["idSha256"]}
            if identity_source != "local_input" or not _identity_fields_valid(identity):
                raise ValueError("input record public identity binding is invalid")
        else:
            identity = public_record_identity(rec["id"], trust=identity_source)
            if identity_source != "local_input" and identity != {"id": rec["id"]}:
                raise ValueError("input record public identity binding is invalid")
        records.append({**identity, "identitySource": identity_source, "length": length, "sha256": digest})
    return records


def _source_provider(value) -> str | None:
    return _SOURCE_PROVIDERS.get(value) if isinstance(value, str) else None


def _public_source_url(value, accession=None, source_label=None) -> str | None:
    """Return an exact accession-bound NCBI or ENA record URL, otherwise ``None``."""
    acc = str(accession or "")
    if not is_public_accession(acc):
        return None
    if not isinstance(value, str) or not value or value != value.strip() or "\\" in value or "%" in value:
        return None
    if any(ord(ch) < 32 or ord(ch) == 127 for ch in value):
        return None
    try:
        parsed = urlsplit(value)
        if (parsed.scheme != "https" or not parsed.hostname or parsed.query or parsed.fragment or
                parsed.username is not None or parsed.password is not None or parsed.port is not None):
            return None
        host = parsed.hostname.rstrip(".").lower()
    except (TypeError, ValueError):
        return None
    provider = None
    for candidate, (expected_host, path_template) in _SOURCE_URLS.items():
        if host == expected_host and parsed.path == path_template.format(accession=acc):
            provider = candidate
            break
    declared_provider = _source_provider(source_label)
    if provider is None or (declared_provider is not None and declared_provider != provider):
        return None
    expected_host, path_template = _SOURCE_URLS[provider]
    canonical = f"https://{expected_host}{path_template.format(accession=acc)}"
    return canonical if value == canonical else None


def _identity_fields_valid(container: dict, key: str = "id") -> bool:
    value = container.get(key)
    digest = container.get(key + "Sha256")
    if not isinstance(value, str) or not is_public_record_id(value):
        return False
    if digest is None:
        return not value.startswith("record_") or bool(re.fullmatch(r"record_[0-9a-f]{12,64}", value))
    return (isinstance(digest, str) and bool(re.fullmatch(r"[0-9a-f]{64}", digest)) and
            value == f"record_{digest[:12]}")


def _identity_binding_valid(container: dict, *, key: str = "id", source_key: str = "identitySource") -> bool:
    """Validate an identity together with the provenance class that licensed its public form."""
    if not isinstance(container, dict) or not _identity_fields_valid(container, key):
        return False
    source = container.get(source_key)
    value = container.get(key)
    digest = container.get(key + "Sha256")
    if source == "local_input":
        return digest is not None
    if source == "validated_fetch":
        return digest is None and public_record_identity(value, trust="validated_fetch") == {"id": value}
    if source == "generated":
        return digest is None and public_record_identity(value, trust="generated") == {"id": value}
    return False


def _nested_identity_bindings_valid(value) -> bool:
    """Reject resealed compact/full-digest mismatches anywhere in sealed parameters."""
    if isinstance(value, list):
        return all(_nested_identity_bindings_valid(item) for item in value)
    if not isinstance(value, dict):
        return True
    for key in tuple(value):
        if key.endswith("Sha256") and key[:-6] in value and not _identity_fields_valid(value, key[:-6]):
            return False
    if "sourceIdSha256" in value:
        source = str(value.get("source", ""))
        for suffix in (" (input)", " (custom background)"):
            if source.endswith(suffix):
                source = source[:-len(suffix)]
                break
        if not _identity_fields_valid({"id": source, "idSha256": value["sourceIdSha256"]}):
            return False
    return all(_nested_identity_bindings_valid(item) for item in value.values())


def _strict_utc(value) -> str | None:
    if not isinstance(value, str) or not _UTC_SECONDS.fullmatch(value):
        return None
    try:
        parsed = datetime.datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError:
        return None
    return value if parsed.utcoffset() == datetime.timedelta(0) else None


def _numeric_taxid(value) -> str | None:
    if value is None or value == "":
        return None
    text = str(value)
    return text if text.isascii() and text.isdigit() and int(text) > 0 else None


def _source_record_observations(records) -> list | None:
    if not isinstance(records, list):
        return None
    observed = []
    for record in records:
        if isinstance(record, (tuple, list)) and len(record) == 2:
            record_id, sequence = record
            if not isinstance(sequence, str):
                return None
            length = len(sequence)
            digest = hashlib.sha256(sequence.encode("utf-8")).hexdigest()
        elif isinstance(record, dict) and "id" in record:
            record_id = record["id"]
            if "sequence" in record:
                sequence = record["sequence"]
                if not isinstance(sequence, str):
                    return None
                length = len(sequence)
                digest = hashlib.sha256(sequence.encode("utf-8")).hexdigest()
            else:
                length, digest = record.get("length"), str(record.get("sha256", "")).lower()
                if (isinstance(length, bool) or not isinstance(length, int) or length < 0 or
                        not re.fullmatch(r"[0-9a-f]{64}", digest)):
                    return None
        else:
            return None
        observed.append({"id": _header_accession(str(record_id).split(None, 1)[0]),
                         "length": length, "sha256": digest})
    return observed


def _header_accession(token: str) -> str:
    """The accession a FASTA header token names: ENA writes 'ENA|<accession>|<accession.version>'."""
    parts = token.split("|")
    return parts[2] if len(parts) == 3 and parts[0] == "ENA" and parts[1] == parts[2].rsplit(".", 1)[0] else token


def _unversioned(accession: str) -> str:
    return accession.rsplit(".", 1)[0]


def _wire_record_binding_valid(binding) -> bool:
    return (isinstance(binding, dict) and set(binding) == {"id", "length", "sha256"} and
            isinstance(binding.get("id"), str) and bool(binding["id"]) and
            isinstance(binding.get("length"), int) and not isinstance(binding["length"], bool) and
            binding["length"] >= 0 and isinstance(binding.get("sha256"), str) and
            bool(re.fullmatch(r"[0-9a-f]{64}", binding["sha256"])))


def normalize_source_context(source, records) -> dict:
    """Normalize one raw fetch response into the exact serializable operation-wire schema."""
    local = {"sourceContextVersion": SOURCE_CONTEXT_VERSION, "inputType": "local_input"}
    if not isinstance(source, dict) or not isinstance(records, list):
        return local
    observations = _source_record_observations(records)
    if not observations:
        return local
    timestamp = _strict_utc(source.get("retrievedUtc"))
    taxid = _numeric_taxid(source.get("taxid"))
    label = source.get("source")
    if label in _FETCH_ENDPOINTS:
        accession = str(source.get("accession") or "")
        provider, endpoint_enum, endpoint = _FETCH_ENDPOINTS[label]
        # fetch.retrieve records the resolved accession-version but builds its URLs from the text the user
        # typed, so an unversioned accession ('M11240') names the same record as 'M11240.1'. The canonical
        # versioned URL is emitted; a URL naming any other record, version or host is still refused.
        forms = (accession, _unversioned(accession))
        expected_endpoints = {endpoint or f"https://www.ebi.ac.uk/ena/browser/api/fasta/{form}" for form in forms}
        source_url = None
        if any(_public_source_url(source.get("sourceUrl"), form, label) for form in forms):
            host, path = _SOURCE_URLS[_source_provider(label)]
            source_url = _public_source_url(f"https://{host}{path.format(accession=accession)}", accession, label)
        if (not timestamp or not is_public_accession(accession) or source.get("endpoint") not in expected_endpoints or
                source_url is None or len(records) != 1):
            return local
        header, sequence = records[0]
        if _header_accession(str(header).split(None, 1)[0]) != accession:
            return local
        length = source.get("length")
        if length is not None and (isinstance(length, bool) or not isinstance(length, int) or length != len(sequence)):
            return local
        out = {"sourceContextVersion": SOURCE_CONTEXT_VERSION, "inputType": "fetched_accession",
               "source": provider, "endpoint": endpoint_enum, "accession": accession,
               "retrievedUtc": timestamp, "sourceUrl": source_url, "records": observations}
        if taxid:
            out["taxid"] = taxid
        return out
    if (label != _COORD_SOURCE or source.get("endpoint") != _COORD_ENDPOINT or
            source.get("retrievalType") != "coordinate" or source.get("coordSystem") != "1-based-inclusive" or
            not timestamp or not is_public_accession(source.get("assemblyAccession")) or
            not isinstance(source.get("regions"), list) or not source["regions"] or
            len(source["regions"]) != len(records)):
        return local
    normalized_regions = []
    for raw, record in zip(source["regions"], records):
        if not isinstance(raw, dict):
            return local
        accession, start, end, strand = raw.get("chrAccession"), raw.get("start"), raw.get("stop"), raw.get("strand")
        if (not is_public_accession(accession) or any(isinstance(v, bool) or not isinstance(v, int)
                                                     for v in (start, end, strand)) or
                start < 1 or end < start or strand not in {1, 2}):
            return local
        expected_header = f"{accession}:{start}-{end}" + ("(-)" if strand == 2 else "")
        if record[0] != expected_header or len(record[1]) != end - start + 1:
            return local
        normalized_regions.append({"chrAccession": accession, "start": start, "end": end, "strand": strand})
        observations[len(normalized_regions) - 1]["id"] = f"{accession}:{start}-{end}"
    first_accession = normalized_regions[0]["chrAccession"]
    out = {"sourceContextVersion": SOURCE_CONTEXT_VERSION, "inputType": "fetched_coordinate",
           "source": "ncbi_coordinate", "endpoint": "ncbi_efetch_range", "accession": first_accession,
           "assemblyAccession": source["assemblyAccession"], "retrievedUtc": timestamp,
           "coordSystem": "1-based-inclusive", "regions": normalized_regions,
           "sourceUrl": f"https://www.ncbi.nlm.nih.gov/nuccore/{first_accession}", "records": observations}
    if taxid:
        out["taxid"] = taxid
    return out


def validate_source_context(source, records) -> dict:
    """Revalidate an exact serialized source context against the actual operation records."""
    local = {"sourceContextVersion": SOURCE_CONTEXT_VERSION, "inputType": "local_input"}
    observations = _source_record_observations(records)
    if not isinstance(source, dict) or not observations:
        return local
    input_type = source.get("inputType")
    common = {"sourceContextVersion", "inputType", "source", "endpoint", "accession",
              "retrievedUtc", "sourceUrl", "records"}
    allowed = common
    if input_type == "fetched_coordinate":
        allowed = allowed | {"assemblyAccession", "coordSystem", "regions"}
    elif input_type != "fetched_accession":
        return local
    if "taxid" in source:
        allowed = allowed | {"taxid"}
    if set(source) != allowed or source.get("sourceContextVersion") != SOURCE_CONTEXT_VERSION:
        return local
    timestamp = _strict_utc(source.get("retrievedUtc"))
    taxid = _numeric_taxid(source.get("taxid"))
    if not timestamp or ("taxid" in source and taxid != source["taxid"]):
        return local
    bindings = source.get("records")
    if (not isinstance(bindings, list) or len(bindings) != len(observations) or
            any(not _wire_record_binding_valid(binding) for binding in bindings)):
        return local
    if input_type == "fetched_accession":
        accession = source.get("accession")
        pair = (source.get("source"), source.get("endpoint"))
        if (pair not in {("ncbi_nuccore", "ncbi_efetch"), ("ena_fasta", "ena_fasta")} or
                not is_public_accession(accession) or len(bindings) != 1 or bindings[0]["id"] != accession or
                bindings != observations):
            return local
        label = "NCBI nuccore" if pair[0] == "ncbi_nuccore" else "ENA (EMBL-EBI)"
        if _public_source_url(source.get("sourceUrl"), accession, label) != source.get("sourceUrl"):
            return local
    else:
        if (source.get("source") != "ncbi_coordinate" or source.get("endpoint") != "ncbi_efetch_range" or
                source.get("coordSystem") != "1-based-inclusive" or
                not is_public_accession(source.get("assemblyAccession")) or
                not isinstance(source.get("regions"), list) or len(source["regions"]) != len(bindings)):
            return local
        for index, (region, binding) in enumerate(zip(source["regions"], bindings)):
            expected_id = f"{region.get('chrAccession')}:{region.get('start')}-{region.get('end')}"
            observed = observations[index]
            if (not isinstance(region, dict) or set(region) != {"chrAccession", "start", "end", "strand"} or
                    not is_public_accession(region.get("chrAccession")) or
                    any(isinstance(region.get(k), bool) or not isinstance(region.get(k), int)
                        for k in ("start", "end", "strand")) or region["start"] < 1 or
                    region["end"] < region["start"] or region["strand"] not in {1, 2} or
                    binding["id"] != expected_id or binding["length"] != region["end"] - region["start"] + 1 or
                    observed["id"] not in ({expected_id, expected_id + "(-)"} if region["strand"] == 2 else {expected_id}) or
                    observed["length"] != binding["length"] or observed["sha256"] != binding["sha256"]):
                return local
        accession = source["regions"][0]["chrAccession"]
        if source.get("accession") != accession or source.get("sourceUrl") != (
                f"https://www.ncbi.nlm.nih.gov/nuccore/{accession}"):
            return local
    return copy.deepcopy(source)


def _semantic_public_fields_valid(manifest: dict) -> bool:
    inp = manifest.get("input")
    if not isinstance(inp, dict) or not _identity_binding_valid(inp):
        return False
    records = inp.get("records", [])
    if (not isinstance(records, list) or any(
            not isinstance(rec, dict) or
            set(rec) != ({"id", "identitySource", "length", "sha256"} |
                         ({"idSha256"} if "idSha256" in rec else set())) or
            not _identity_binding_valid(rec)
            for rec in records)):
        return False
    if not _nested_identity_bindings_valid(manifest.get("parameters")):
        return False
    # input.id is outside the seal, so it is tied to the sealed per-record bindings instead: a generated
    # label must not be rewritable
    if len(records) == 1 and any(inp.get(k) != records[0].get(k) for k in ("id", "idSha256", "identitySource")):
        return False
    if len(records) > 1 and (inp["id"] != "multiple_records" or inp["identitySource"] != "generated"):
        return False
    input_type = inp.get("inputType", "local_input")
    if (inp.get("sourceContextVersion") != SOURCE_CONTEXT_VERSION or
            input_type not in {"local_input", "fetched_accession", "fetched_coordinate"}):
        return False
    base_keys = {"id", "identitySource", "length", "sha256", "sourceContextVersion", "inputType"}
    if "idSha256" in inp:
        base_keys.add("idSha256")
    if records:
        base_keys.add("records")
    source_keys = {
        "fetched_accession": {"source", "endpoint", "accession", "retrievedUtc", "sourceUrl"},
        "fetched_coordinate": {"source", "endpoint", "accession", "assemblyAccession", "retrievedUtc",
                               "coordSystem", "regions", "sourceUrl"},
    }.get(input_type, set())
    if "taxid" in inp:
        source_keys.add("taxid")
    if set(inp) != base_keys | source_keys:
        return False
    if input_type == "local_input" and inp.get("identitySource") == "validated_fetch":
        return False
    multi_region = (input_type == "fetched_coordinate" and isinstance(inp.get("regions"), list) and
                    len(inp["regions"]) > 1)
    if input_type != "local_input" and inp.get("identitySource") != ("generated" if multi_region else "validated_fetch"):
        return False
    if input_type == "local_input" and any(key in inp for key in
            ("accession", "assemblyAccession", "taxid", "source", "endpoint", "sourceUrl", "regions")):
        return False
    if input_type != "local_input":
        expected_pair = {"fetched_accession": {("ncbi_nuccore", "ncbi_efetch"),
                                                ("ena_fasta", "ena_fasta")},
                         "fetched_coordinate": {("ncbi_coordinate", "ncbi_efetch_range")}}[input_type]
        if (not is_public_accession(inp.get("accession")) or
                (inp.get("source"), inp.get("endpoint")) not in expected_pair or
                not _strict_utc(inp.get("retrievedUtc")) or
                ("taxid" in inp and _numeric_taxid(inp["taxid"]) != inp["taxid"])):
            return False
        if any(rec.get("identitySource") != "validated_fetch" for rec in records):
            return False
    elif any(rec.get("identitySource") == "validated_fetch" for rec in records):
        return False
    if input_type == "fetched_coordinate" and not is_public_accession(inp.get("assemblyAccession")):
        return False
    if input_type == "fetched_coordinate" and not inp.get("regions"):
        return False
    for region in inp.get("regions", []):
        if (not isinstance(region, dict) or set(region) != {"chrAccession", "start", "end", "strand"} or
                not is_public_accession(region.get("chrAccession")) or
                any(isinstance(region.get(k), bool) or not isinstance(region.get(k), int)
                    for k in ("start", "end", "strand")) or region["start"] < 1 or
                region["end"] < region["start"] or region["strand"] not in {1, 2}):
            return False
    if input_type == "fetched_accession":
        expected_ids = [inp["accession"]]
        if inp.get("id") != inp["accession"] or (records and
                [record.get("id") for record in records] != expected_ids):
            return False
    elif input_type == "fetched_coordinate":
        expected_ids = [f"{region['chrAccession']}:{region['start']}-{region['end']}"
                        for region in inp["regions"]]
        if records and [record.get("id") for record in records] != expected_ids:
            return False
        if ((len(expected_ids) == 1 and inp.get("id") != expected_ids[0]) or
                (len(expected_ids) > 1 and
                 (inp.get("id") != "multiple_records" or inp.get("identitySource") != "generated"))):
            return False
    if "sourceUrl" in inp:
        label = {"ncbi_nuccore": "NCBI nuccore", "ena_fasta": "ENA (EMBL-EBI)",
                 "ncbi_coordinate": _COORD_SOURCE}.get(inp.get("source"))
        if _public_source_url(inp["sourceUrl"], inp.get("accession"), label) != inp["sourceUrl"]:
            return False
    return True


def public_source_fields(source, records=None) -> dict:
    """Return manifest source fields only after exact wire-schema and record-binding validation."""
    validated = validate_source_context(source, records)
    validated.pop("records", None)
    return validated


def _software(run_type: str) -> list:
    """Only the tools that actually produced THIS run enter the software list (and thus the seal),
    so an unused tool's version/availability can never change a byte-identical run's checksum."""
    sw = [{"name": "python", "version": platform.python_version()},
          {"name": "TEagle core", "version": TEAGLE_VERSION}]
    # the versions are those of the Python bindings; the wrapped Primer3 core release is read from the library
    if run_type == "primer":
        sw.insert(0, {"name": "primer3-py (Primer3 bindings)", "version": PRIMER3_VERSION})
        core = _libprimer3_release()
        if core:
            sw.insert(1, {"name": "libprimer3 (Primer3 core)", "version": core})
    elif run_type == "analysis":
        sw.insert(0, {"name": "pyhmmer (HMMER 3 bindings)", "version": PYHMMER_VERSION})
    return sw


def _libprimer3_release() -> str | None:
    try:
        from primer3.thermoanalysis import get_libprimer3_version
        return str(get_libprimer3_version()).rsplit(None, 1)[-1]
    except Exception:
        return None


_VOLATILE_SOURCE_FIELDS = ("retrievedUtc", "source", "endpoint", "sourceUrl")


def _integral_floats_as_int(value):
    """Write an integral float (60.0) as an int in sealed parameters. Python keeps the two apart when it
    serializes, so the seal would otherwise depend on how a value was typed (opt_tm 60 or 60.0) and would
    break on a re-serialization by a JSON tool that does not (JavaScript, jq)."""
    if isinstance(value, float) and value.is_integer():
        return int(value)
    if isinstance(value, dict):
        return {k: _integral_floats_as_int(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_integral_floats_as_int(v) for v in value]
    return value


def _without_volatile_bindings(params):
    """Drop the retrieval time and serving-database labels from the input bindings nested in sealed parameters.

    In-silico PCR copies the assay and design input bindings into its parameters. The top-level input already
    keeps these fields out of the seal; a copy left inside `parameters` would make the same primers and the
    same sequence seal differently after a refetch or from a different serving database."""
    def strip(binding):
        return ({k: v for k, v in binding.items() if k not in _VOLATILE_SOURCE_FIELDS}
                if isinstance(binding, dict) else binding)

    if not isinstance(params, dict):
        return params
    out = _integral_floats_as_int(params)
    for key in ("design_input", "assay_input"):
        if key in out:
            out[key] = strip(out[key])
    if isinstance(out.get("search_inputs"), list):
        out["search_inputs"] = [strip(binding) for binding in out["search_inputs"]]
    return out


def build_manifest(run_type: str, input_seq: str, input_id: str, params: dict,
                   not_run: list | None = None, source: dict | None = None,
                   references: list | None = None, databases: list | None = None,
                   input_sha256: str | None = None, input_length: int | None = None,
                   input_records: list | None = None, result_descriptor: dict | None = None,
                   input_identity: dict | None = None):
    """Build a content-addressed manifest.

    ``input_sha256`` and ``input_length`` bind a previously verified, non-materialised input such as a
    multi-gigabyte cached genome. Ordinary sequence workflows omit them and retain the original behavior.
    Whole-genome annotation must supply both; hashing its intentionally empty ``input_seq`` would otherwise
    seal the empty string rather than the assembly that RepeatMasker read.
    """
    if (input_sha256 is None) != (input_length is None):
        raise ValueError("input_sha256 and input_length must be supplied together")
    if input_sha256 is not None:
        candidate = str(input_sha256).lower()
        if len(candidate) != 64 or any(ch not in "0123456789abcdef" for ch in candidate):
            raise ValueError("input_sha256 must be a 64-character hexadecimal SHA-256")
        input_hash = candidate
    else:
        input_hash = hashlib.sha256(input_seq.encode()).hexdigest()
    if input_length is not None:
        if isinstance(input_length, bool) or not isinstance(input_length, int) or input_length < 0:
            raise ValueError("input_length must be a non-negative integer")
        recorded_length = input_length
    else:
        recorded_length = len(input_seq)
    if run_type == "genome_annotate" and input_sha256 is None:
        raise ValueError("genome_annotate requires the verified genome SHA-256 and length")
    if run_type == "genome_annotate" and recorded_length < 1:
        raise ValueError("genome_annotate requires a positive verified genome length")
    if result_descriptor is not None:
        if run_type != "analysis" or set(result_descriptor) != {"resultSchemaVersion", "resultPayloadSha256"}:
            raise ValueError("result_descriptor is valid only for analysis and must use the defined fields")
        digest = result_descriptor.get("resultPayloadSha256")
        if (not isinstance(result_descriptor.get("resultSchemaVersion"), str) or
                not result_descriptor["resultSchemaVersion"] or not isinstance(digest, str) or
                len(digest) != 64 or any(ch not in "0123456789abcdef" for ch in digest.lower())):
            raise ValueError("invalid analysis result descriptor")
    identity = copy.deepcopy(input_identity) if input_identity is not None else public_record_identity(input_id)
    identity_source = identity.pop("identitySource", None)
    if identity_source is None:
        identity_source = ("local_input" if identity.get("idSha256") else "generated")
    identity["identitySource"] = identity_source
    if not _identity_binding_valid(identity):
        raise ValueError("input public identity binding is invalid")
    inp = {**identity, "length": recorded_length, "sha256": input_hash}
    records = _record_descriptors(input_records)
    if records is not None:
        inp["records"] = records
    source_records = (input_records if input_records is not None else
                      [{"id": input_id, "sequence": input_seq}])
    inp.update(public_source_fields(source, source_records))
    m = {
        "teagleVersion": TEAGLE_VERSION,
        "schemaVersion": SCHEMA_VERSION,
        "canonicalization": CANONICALIZATION,
        "runType": run_type,
        "createdUtc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        "input": inp,
        "software": _software(run_type),
        "databases": databases or [],
        "parameters": _without_volatile_bindings(params),
        "references": references or [],
        "environment": {
            "os": platform.platform(),
            "machine": platform.machine(),
            "python": sys.version.split()[0],
        },
        "notRun": not_run or [
            "DB-backed family naming (Dfam / RepeatMasker) — available via the WSL backend; not run in this step",
            "External NCBI Primer-BLAST",
            "De novo family discovery",
        ],
    }
    if result_descriptor is not None:
        m["resultDescriptor"] = copy.deepcopy(result_descriptor)
    # Content-addressed seal over the scientific content only: exclude wall-clock time, the machine,
    # and the volatile retrieval timestamp, so an identical run (same input, versions, params, DBs)
    # yields an identical seal across machines and refetches.
    m["manifestSha256"] = hashlib.sha256(canonical_manifest_bytes(m)).hexdigest()
    return m
