"""Build a compact, sequence-free index of the per-case TEagle manifests.

The locally retained raw manifests are deliberately not copied to the reviewer package:
they contain fetched sequence and detector payloads.  This program extracts only the identities,
hashes, sealed settings, versions, classifications, and bounded LTR/TIR/TSD
measurements needed to audit the benchmark.  It never executes TEagle or a
comparator.

Run from the repository root::

    python benchmarks/make_provenance_index.py
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
from collections import Counter, defaultdict
from pathlib import Path, PurePosixPath
from tempfile import NamedTemporaryFile
from typing import Any, Iterable


def _find_repository_root(script_path: Path) -> Path:
    """Allow the byte-identical reviewer-package mirror to run from either path."""
    for candidate in script_path.resolve().parents:
        if (candidate / "benchmarks" / "raw").is_dir():
            return candidate
    raise RuntimeError("could not locate repository root containing benchmarks/raw")


ROOT = _find_repository_root(Path(__file__))
RAW_GROUPS = (
    ("broader", ROOT / "benchmarks" / "raw" / "teagle"),
    ("holdout", ROOT / "benchmarks" / "raw" / "teagle_holdout"),
)
CORPUS_FILES = {
    "broader": ROOT / "benchmarks" / "corpus.tsv",
    "holdout": ROOT / "benchmarks" / "corpus_holdout.tsv",
}
DEFAULT_OUTPUTS = (
    ROOT / "benchmarks" / "provenance",
    ROOT / "reviewer_package" / "results",
)
SCHEMA_VERSION = "1.1"
MAX_PUBLIC_STRING = 256
MAX_COMPLETENESS_STRING = 512
MAX_PARAMETER_COUNT = 64
MAX_PARAMETER_KEY_LENGTH = 64
MAX_ABS_PARAMETER_VALUE = 1_000_000_000_000
MAX_SOFTWARE_ENTRIES = 16
MAX_DATABASE_ENTRIES = 16
MAX_RAW_STRUCTURAL_ENTRIES = 32
MAX_TERMINAL_REPEAT_ENTRIES = 2
MAX_TSD_ENTRIES = 1
PARAMETER_KEY_RE = re.compile(r"^[a-z][a-z0-9]*(?:_[a-z0-9]+)*$")
TSD_STATES = frozenset(
    {
        "found",
        "not assessable (no flanking sequence in the record)",
        "absent (flanks present, no direct repeat found)",
    }
)
TSD_CONGRUENCE = frozenset(
    {"congruent", "incongruent", "not assessed", "not diagnostic", "not recorded"}
)
SEAL_EXCLUDE_TOP = frozenset({"manifestSha256", "createdUtc", "environment"})
SEAL_EXCLUDE_INPUT = frozenset(
    {
        "retrievedUtc", "id", "assemblyName", "displayLocus", "chromName",
        "source", "endpoint", "sourceUrl",
    }
)

# These are deliberately exhaustive: a new field must be reviewed before it can
# enter the public compact extract.  In particular, no source sequence, peptide,
# primer, coordinate, motif, free-text explanation, or nested raw payload is copied.
CASE_KEYS = frozenset(
    {
        "case_id", "corpus", "accession", "panel", "record_scope",
        "input_length", "input_sha256", "sealed_input_sha256",
        "source_manifest", "source_manifest_sha256", "internal_manifest_sha256",
        "provenance_schema_version", "run_type", "teagle_version", "software",
        "databases", "parameters", "records",
    }
)
RECORD_KEYS = frozenset(
    {"class", "order", "superfamily", "confidence", "completeness", "structural"}
)
STRUCTURAL_KEYS = frozenset({"terminal_repeats", "target_site_duplications"})
TERMINAL_REPEAT_KEYS = frozenset(
    {"kind", "length", "identity_pct", "length_is_lower_bound", "threshold", "tsd_state"}
)
TSD_KEYS = frozenset({"length", "matched_expected", "congruence"})
# Parameter names such as ``orf_min_aa`` are required sealed numeric evidence,
# so key checking is exact rather than a broad substring ban.
FORBIDDEN_KEYS = frozenset(
    {
        "seq", "sequence", "seq_preview", "preview", "orfs", "orf", "primers", "primer",
        "protein", "peptide", "dna", "motif", "coords", "coordinates", "element_span",
        "five_prime", "three_prime", "upstream", "downstream", "pos", "position",
    }
)


class ProvenanceIndexError(ValueError):
    """A raw manifest cannot safely support the compact public extract."""


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _bounded_string(value: Any, context: str, *, maximum: int = MAX_PUBLIC_STRING,
                    optional: bool = False) -> str | None:
    if value is None and optional:
        return None
    if not isinstance(value, str) or not value or len(value) > maximum:
        raise ProvenanceIndexError(f"{context}: expected a nonempty string of at most {maximum} characters")
    if any(ord(character) < 32 or ord(character) == 127 for character in value):
        raise ProvenanceIndexError(f"{context}: control characters are not permitted")
    return value


def _bounded_number(value: Any, context: str, *, minimum: float, maximum: float,
                    integer: bool = False) -> int | float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ProvenanceIndexError(f"{context}: expected a numeric scalar (boolean is not numeric here)")
    if not math.isfinite(value) or value < minimum or value > maximum:
        raise ProvenanceIndexError(f"{context}: numeric value is non-finite or outside [{minimum}, {maximum}]")
    if integer and not isinstance(value, int):
        raise ProvenanceIndexError(f"{context}: expected an integer")
    return value


def _normalise_parameters(value: Any, context: str) -> dict[str, int | float]:
    if not isinstance(value, dict) or not value or len(value) > MAX_PARAMETER_COUNT:
        raise ProvenanceIndexError(
            f"{context}: parameters must be a nonempty flat mapping with at most {MAX_PARAMETER_COUNT} entries"
        )
    compact: dict[str, int | float] = {}
    for key, parameter_value in value.items():
        if (not isinstance(key, str) or len(key) > MAX_PARAMETER_KEY_LENGTH
                or not PARAMETER_KEY_RE.fullmatch(key) or key.startswith("_")):
            raise ProvenanceIndexError(f"{context}: unsafe parameter key {key!r}")
        compact[key] = _bounded_number(
            parameter_value, f"{context}: parameter {key!r}",
            minimum=-MAX_ABS_PARAMETER_VALUE, maximum=MAX_ABS_PARAMETER_VALUE,
        )
    return {key: compact[key] for key in sorted(compact)}


def compute_internal_manifest_sha256(provenance: dict[str, Any]) -> str:
    """Recompute TEagle's content seal independently from the original manifest.

    This intentionally mirrors the documented algorithm, not a stored summary:
    remove the existing digest, wall clock and environment; remove only the
    documented volatile input labels; then hash canonical compact JSON.
    """
    if not isinstance(provenance, dict):
        raise ProvenanceIndexError("provenance seal: manifest is not an object")
    seal_source = {key: value for key, value in provenance.items() if key not in SEAL_EXCLUDE_TOP}
    sealed_input = seal_source.get("input")
    if not isinstance(sealed_input, dict):
        raise ProvenanceIndexError("provenance seal: input is not an object")
    seal_source["input"] = {
        key: value for key, value in sealed_input.items() if key not in SEAL_EXCLUDE_INPUT
    }
    try:
        encoded = json.dumps(
            seal_source, sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode("utf-8")
    except (TypeError, ValueError) as error:
        raise ProvenanceIndexError(f"provenance seal: content is not canonical JSON: {error}") from error
    return sha256_bytes(encoded)


def _require(mapping: dict[str, Any], key: str, context: str) -> Any:
    value = mapping.get(key)
    if value is None or value == "" or value == [] or value == {}:
        raise ProvenanceIndexError(f"{context}: missing required {key!r}")
    return value


def _require_sha256(value: Any, context: str) -> str:
    if not isinstance(value, str) or len(value) != 64 or any(c not in "0123456789abcdef" for c in value):
        raise ProvenanceIndexError(f"{context}: required SHA-256 is absent or malformed")
    return value


def _normalised_text_sha256(path: Path) -> str:
    return sha256_bytes(path.read_text(encoding="utf-8").encode("utf-8"))


def validated_source_paths(
    raw_groups: tuple[tuple[str, Path], ...] = RAW_GROUPS,
    corpus_files: dict[str, Path] | None = CORPUS_FILES,
) -> list[tuple[str, Path]]:
    """Fail closed unless every configured run ledger and successful case is present."""
    sources: list[tuple[str, Path]] = []
    for corpus, directory in raw_groups:
        if not directory.is_dir():
            raise ProvenanceIndexError(f"{corpus}: locally retained per-case directory is missing: {directory}")
        run_path = directory / "_run.json"
        if not run_path.is_file():
            raise ProvenanceIndexError(f"{corpus}: run ledger is missing: {run_path}")
        try:
            run = json.loads(run_path.read_text(encoding="utf-8"))
        except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
            raise ProvenanceIndexError(f"{corpus}: run ledger cannot be read: {error}") from error
        counts = run.get("counts")
        if not isinstance(counts, dict):
            raise ProvenanceIndexError(f"{corpus}: run ledger lacks counts")
        corpus_rows = run.get("corpus_rows")
        analysed = counts.get("analysed")
        skipped = counts.get("skipped_existing")
        failed = counts.get("failed")
        if (any(isinstance(value, bool) or not isinstance(value, int) or value < 0
                for value in (corpus_rows, analysed, skipped, failed)) or corpus_rows == 0):
            raise ProvenanceIndexError(f"{corpus}: run ledger counts are malformed")
        expected_successes = corpus_rows - failed
        if expected_successes <= 0 or analysed + skipped != expected_successes:
            raise ProvenanceIndexError(
                f"{corpus}: run ledger accounting disagrees: rows={corpus_rows}, analysed={analysed}, "
                f"skipped={skipped}, failed={failed}"
            )
        case_paths = sorted(path for path in directory.glob("*.json") if path.name != "_run.json")
        if len(case_paths) != expected_successes:
            raise ProvenanceIndexError(
                f"{corpus}: expected {expected_successes} successful per-case manifests from the run "
                f"ledger, found {len(case_paths)}"
            )
        if corpus_files is not None:
            corpus_path = corpus_files.get(corpus)
            if corpus_path is None or not corpus_path.is_file():
                raise ProvenanceIndexError(f"{corpus}: configured corpus file is missing")
            with corpus_path.open(encoding="utf-8", newline="") as handle:
                row_count = sum(1 for _ in csv.DictReader(handle, delimiter="\t"))
            if row_count != corpus_rows:
                raise ProvenanceIndexError(
                    f"{corpus}: corpus has {row_count} rows but the run ledger records {corpus_rows}"
                )
            if run.get("corpus_sha256") != _normalised_text_sha256(corpus_path):
                raise ProvenanceIndexError(f"{corpus}: corpus hash disagrees with the run ledger")
        sources.extend((corpus, path) for path in case_paths)
    return sources


def _normalise_software(items: Any, context: str) -> list[dict[str, str]]:
    if not isinstance(items, list) or not items or len(items) > MAX_SOFTWARE_ENTRIES:
        raise ProvenanceIndexError(f"{context}: missing software versions")
    compact = []
    seen: dict[str, str] = {}
    for item in items:
        if not isinstance(item, dict) or set(item) != {"name", "version"}:
            raise ProvenanceIndexError(f"{context}: malformed software entry")
        name = _bounded_string(item.get("name"), f"{context}: software name", maximum=128)
        version = _bounded_string(item.get("version"), f"{context}: software version", maximum=64)
        assert name is not None and version is not None
        if name in seen:
            detail = "conflicting versions" if seen[name] != version else "duplicate name"
            raise ProvenanceIndexError(f"{context}: software {name!r} has {detail}")
        seen[name] = version
        compact.append({"name": name, "version": version})
    return sorted(compact, key=lambda x: (x["name"], x["version"]))


def _normalise_databases(items: Any, context: str) -> list[dict[str, str]]:
    if not isinstance(items, list) or not items or len(items) > MAX_DATABASE_ENTRIES:
        raise ProvenanceIndexError(f"{context}: missing database checksums")
    compact = []
    seen: dict[str, tuple[str | None, str]] = {}
    for item in items:
        if not isinstance(item, dict) or not set(item).issubset({"name", "sha256", "version", "file"}):
            raise ProvenanceIndexError(f"{context}: malformed database entry")
        name = _bounded_string(item.get("name"), f"{context}: database name", maximum=128)
        assert name is not None
        digest = _require_sha256(item.get("sha256"), context)
        version = None
        for label in ("version", "file"):
            if label in item:
                bounded = _bounded_string(item[label], f"{context}: database {label}", maximum=128)
                if label == "version":
                    version = bounded
        if name in seen:
            detail = "conflicting versions/checksums" if seen[name] != (version, digest) else "duplicate name"
            raise ProvenanceIndexError(f"{context}: database {name!r} has {detail}")
        seen[name] = (version, digest)
        entry = {"name": name, "sha256": digest}
        if version is not None:
            entry["version"] = version
        compact.append(entry)
    return sorted(compact, key=lambda x: (x["name"], x["sha256"]))


def _terminal_repeat(item: dict[str, Any], *, input_length: int) -> dict[str, Any] | None:
    kind = _bounded_string(item.get("type"), "structural type", maximum=96)
    if kind == "LTR (terminal direct repeat)":
        lower_bound = item.get("length_is_lower_bound", False)
        if not isinstance(lower_bound, bool):
            raise ProvenanceIndexError("LTR structural record: length_is_lower_bound must be boolean")
        tsd_state = _bounded_string(item.get("tsd_state"), "LTR structural record: tsd_state", maximum=80)
        if tsd_state not in TSD_STATES:
            raise ProvenanceIndexError("LTR structural record: unknown tsd_state")
        return {
            "kind": "LTR",
            "length": _bounded_number(item.get("ltr_len"), "LTR structural record: length",
                                      minimum=1, maximum=input_length, integer=True),
            "identity_pct": _bounded_number(item.get("identity"), "LTR structural record: identity",
                                             minimum=0, maximum=100),
            "length_is_lower_bound": lower_bound,
            "tsd_state": tsd_state,
        }
    if kind == "Sub-threshold terminal direct repeat (advisory)":
        return {
            "kind": "LTR_advisory_subthreshold",
            "length": _bounded_number(item.get("ltr_len"), "advisory LTR structural record: length",
                                      minimum=1, maximum=input_length, integer=True),
            "identity_pct": _bounded_number(item.get("identity"), "advisory LTR structural record: identity",
                                             minimum=0, maximum=100),
            "threshold": _bounded_number(item.get("threshold"), "advisory LTR structural record: threshold",
                                         minimum=0, maximum=100),
        }
    if kind == "TIR (terminal inverted repeat)":
        tsd_state = _bounded_string(item.get("tsd_state"), "TIR structural record: tsd_state", maximum=80)
        if tsd_state not in TSD_STATES:
            raise ProvenanceIndexError("TIR structural record: unknown tsd_state")
        return {
            "kind": "TIR",
            "length": _bounded_number(item.get("tir_len"), "TIR structural record: length",
                                      minimum=1, maximum=input_length, integer=True),
            "identity_pct": _bounded_number(item.get("identity"), "TIR structural record: identity",
                                             minimum=0, maximum=100),
            "tsd_state": tsd_state,
        }
    return None


def _compact_structural(items: Any, *, input_length: int,
                        parameters: dict[str, int | float]) -> dict[str, list[dict[str, Any]]]:
    terminal_repeats: list[dict[str, Any]] = []
    target_site_duplications: list[dict[str, Any]] = []
    if not isinstance(items, list) or len(items) > MAX_RAW_STRUCTURAL_ENTRIES:
        raise ProvenanceIndexError("structural payload must be a bounded list")
    tsd_min = _bounded_number(parameters.get("tsd_min_len"), "parameter tsd_min_len",
                              minimum=1, maximum=input_length, integer=True)
    tsd_max = _bounded_number(parameters.get("tsd_max_len"), "parameter tsd_max_len",
                              minimum=tsd_min, maximum=input_length, integer=True)
    for item in items:
        if not isinstance(item, dict):
            raise ProvenanceIndexError("structural entry must be an object")
        terminal_repeat = _terminal_repeat(item, input_length=input_length)
        if terminal_repeat is not None:
            terminal_repeats.append(terminal_repeat)
        if item.get("type") == "TSD (target-site duplication)":
            matched_expected = item.get("matched_expected")
            if not isinstance(matched_expected, bool):
                raise ProvenanceIndexError("TSD structural record: matched_expected must be boolean")
            raw_congruence = item.get("tsd_congruence")
            congruence = "not recorded" if raw_congruence is None else _bounded_string(
                raw_congruence, "TSD structural record: congruence", maximum=32
            )
            if congruence not in TSD_CONGRUENCE:
                raise ProvenanceIndexError("TSD structural record: unknown congruence")
            target_site_duplications.append(
                {
                    "length": _bounded_number(item.get("length"), "TSD structural record: length",
                                              minimum=tsd_min, maximum=tsd_max, integer=True),
                    "matched_expected": matched_expected,
                    "congruence": congruence,
                }
            )
    if len(terminal_repeats) > MAX_TERMINAL_REPEAT_ENTRIES:
        raise ProvenanceIndexError("analysis record exceeds the terminal-repeat result cap")
    terminal_kinds = [entry["kind"] for entry in terminal_repeats]
    if len(terminal_kinds) != len(set(terminal_kinds)):
        raise ProvenanceIndexError("analysis record contains duplicate terminal-repeat kinds")
    if "LTR" in terminal_kinds and len(terminal_kinds) != 1:
        raise ProvenanceIndexError("accepted LTR cannot coexist with another terminal-repeat result")
    if len(target_site_duplications) > MAX_TSD_ENTRIES:
        raise ProvenanceIndexError("analysis record contains more than one TSD result")
    return {"terminal_repeats": terminal_repeats, "target_site_duplications": target_site_duplications}


def _compact_record(record: dict[str, Any], context: str, *, input_length: int,
                    parameters: dict[str, int | float]) -> dict[str, Any]:
    if not isinstance(record, dict):
        raise ProvenanceIndexError(f"{context}: analysis record is not an object")
    classification = record.get("classification")
    if not isinstance(classification, dict):
        raise ProvenanceIndexError(f"{context}: missing classification")
    completeness = classification.get("completeness")
    # A string tier is the only completeness value needed for the benchmark claim;
    # the detailed domain/orf payload is deliberately excluded.
    if completeness is not None and not isinstance(completeness, dict):
        raise ProvenanceIndexError(f"{context}: malformed completeness object")
    tier = completeness.get("tier") if isinstance(completeness, dict) else None
    return {
        "class": _bounded_string(classification.get("class"), f"{context}: classification class"),
        "order": _bounded_string(classification.get("order"), f"{context}: classification order", optional=True),
        "superfamily": _bounded_string(classification.get("superfamily"), f"{context}: classification superfamily"),
        "confidence": _bounded_string(classification.get("confidence"), f"{context}: classification confidence",
                                      maximum=32),
        "completeness": _bounded_string(tier, f"{context}: completeness tier",
                                        maximum=MAX_COMPLETENESS_STRING, optional=True),
        "structural": _compact_structural(
            record.get("structural"), input_length=input_length, parameters=parameters
        ),
    }


def compact_manifest(raw: dict[str, Any], *, corpus: str, source_manifest: str,
                     source_manifest_sha256: str) -> dict[str, Any]:
    """Validate and reduce one raw per-case output to the public allowlist."""
    if not isinstance(raw, dict):
        raise ProvenanceIndexError("source manifest is not an object")
    corpus = _bounded_string(corpus, "corpus", maximum=32)
    if corpus not in {"broader", "holdout"}:
        raise ProvenanceIndexError(f"unsupported corpus {corpus!r}")
    source_manifest = _bounded_string(source_manifest, "source manifest path", maximum=256)
    assert source_manifest is not None
    source_path = PurePosixPath(source_manifest)
    expected_prefix = "teagle" if corpus == "broader" else "teagle_holdout"
    if ("\\" in source_manifest or ":" in source_manifest or source_path.is_absolute()
            or len(source_path.parts) != 2 or source_path.parts[0] != expected_prefix
            or source_path.name in {"", ".", ".."} or source_path.suffix != ".json"):
        raise ProvenanceIndexError(
            f"source manifest must be {expected_prefix}/<basename>.json without traversal"
        )
    context = source_manifest
    corpus_row = raw.get("corpus_row")
    result = raw.get("result")
    if not isinstance(corpus_row, dict) or not isinstance(result, dict):
        raise ProvenanceIndexError(f"{context}: missing corpus row or result")
    provenance = result.get("provenance")
    if not isinstance(provenance, dict):
        raise ProvenanceIndexError(f"{context}: missing provenance")
    accession = _bounded_string(raw.get("accession"), f"{context}: accession", maximum=128)
    corpus_accession = _bounded_string(corpus_row.get("accession"), f"{context}: corpus accession", maximum=128)
    if accession != corpus_accession:
        raise ProvenanceIndexError(f"{context}: top-level and corpus-row accession differ")
    parameters = _normalise_parameters(provenance.get("parameters"), context)
    records = result.get("records")
    if not isinstance(records, list) or len(records) != 1:
        raise ProvenanceIndexError(f"{context}: benchmark requires exactly one analysis record")
    sealed_input = provenance.get("input")
    if not isinstance(sealed_input, dict):
        raise ProvenanceIndexError(f"{context}: missing sealed input")
    input_length = _bounded_number(raw.get("input_length"), f"{context}: input length",
                                   minimum=1, maximum=1_000_000_000, integer=True)
    sealed_input_length = _bounded_number(sealed_input.get("length"), f"{context}: sealed input length",
                                          minimum=1, maximum=1_000_000_000, integer=True)
    if sealed_input_length != input_length:
        raise ProvenanceIndexError(
            f"{context}: provenance.input.length {sealed_input_length} != analysed input_length {input_length}"
        )
    provenance_schema = _bounded_string(
        provenance.get("schemaVersion"), f"{context}: provenance schema version", maximum=32
    )
    run_type = _bounded_string(provenance.get("runType"), f"{context}: run type", maximum=32)
    if run_type != "analysis":
        raise ProvenanceIndexError(f"{context}: expected analysis run, got {run_type!r}")
    teagle_version = _bounded_string(
        provenance.get("teagleVersion"), f"{context}: TEagle version", maximum=64
    )
    software = _normalise_software(provenance.get("software"), context)
    databases = _normalise_databases(provenance.get("databases"), context)
    core_versions = [item["version"] for item in software if item["name"] == "TEagle core"]
    if core_versions != [teagle_version]:
        raise ProvenanceIndexError(
            f"{context}: TEagle core software version must occur once and equal teagleVersion"
        )
    stored_internal_sha256 = _require_sha256(provenance.get("manifestSha256"), context)
    recomputed_internal_sha256 = compute_internal_manifest_sha256(provenance)
    if recomputed_internal_sha256 != stored_internal_sha256:
        raise ProvenanceIndexError(f"{context}: internal manifestSha256 verification failed")
    compact = {
        "case_id": f"{corpus}:{source_path.stem}",
        "corpus": corpus,
        "accession": accession,
        "panel": _bounded_string(corpus_row.get("panel"), f"{context}: panel", maximum=64),
        "record_scope": _bounded_string(corpus_row.get("record_scope"), f"{context}: record scope", maximum=32),
        "input_length": input_length,
        "input_sha256": _require_sha256(raw.get("input_sha256"), context),
        "sealed_input_sha256": _require_sha256(sealed_input.get("sha256"), context),
        "source_manifest": source_manifest,
        "source_manifest_sha256": _require_sha256(source_manifest_sha256, context),
        "internal_manifest_sha256": stored_internal_sha256,
        "provenance_schema_version": provenance_schema,
        "run_type": run_type,
        "teagle_version": teagle_version,
        "software": software,
        "databases": databases,
        "parameters": parameters,
        "records": [
            _compact_record(records[0], context, input_length=input_length, parameters=parameters)
        ],
    }
    if compact["record_scope"] not in {"element", "containing_record"}:
        raise ProvenanceIndexError(f"{context}: unknown record_scope")
    assert_safe_case(compact)
    return compact


def _walk_keys(value: Any) -> Iterable[str]:
    if isinstance(value, dict):
        for key, child in value.items():
            yield key
            yield from _walk_keys(child)
    elif isinstance(value, list):
        for child in value:
            yield from _walk_keys(child)


def assert_safe_case(case: dict[str, Any]) -> None:
    """Fail if the allowlist is bypassed or a sequence-bearing field leaks."""
    if not isinstance(case, dict) or set(case) != CASE_KEYS:
        raise ProvenanceIndexError(f"case allowlist mismatch: {sorted(set(case) ^ CASE_KEYS)}")
    for key, maximum in {
        "case_id": 320, "corpus": 32, "accession": 128, "panel": 64,
        "record_scope": 32, "source_manifest": 256, "provenance_schema_version": 32,
        "run_type": 32, "teagle_version": 64,
    }.items():
        _bounded_string(case[key], f"compact case {key}", maximum=maximum)
    source_manifest = case["source_manifest"]
    source_path = PurePosixPath(source_manifest)
    expected_prefix = "teagle" if case["corpus"] == "broader" else "teagle_holdout"
    if ("\\" in source_manifest or ":" in source_manifest or source_path.is_absolute()
            or len(source_path.parts) != 2 or source_path.parts[0] != expected_prefix
            or source_path.name in {"", ".", ".."} or source_path.suffix != ".json"):
        raise ProvenanceIndexError("compact source_manifest violates the public relative-path policy")
    input_length = _bounded_number(case["input_length"], "compact case input_length",
                                   minimum=1, maximum=1_000_000_000, integer=True)
    for key in ("input_sha256", "sealed_input_sha256", "source_manifest_sha256",
                "internal_manifest_sha256"):
        _require_sha256(case[key], f"compact case {key}")
    parameters = _normalise_parameters(case["parameters"], "compact case")
    if parameters != case["parameters"]:
        raise ProvenanceIndexError("compact parameters are not canonically sorted")
    software = _normalise_software(case["software"], "compact case")
    databases = _normalise_databases(case["databases"], "compact case")
    if software != case["software"] or databases != case["databases"]:
        raise ProvenanceIndexError("compact software/database entries are not canonical")
    core_versions = [item["version"] for item in software if item["name"] == "TEagle core"]
    if core_versions != [case["teagle_version"]]:
        raise ProvenanceIndexError("compact TEagle core version disagrees with teagle_version")
    if not isinstance(case["records"], list) or len(case["records"]) != 1:
        raise ProvenanceIndexError("compact case must contain exactly one analysis record")
    for record in case["records"]:
        if (not isinstance(record, dict) or set(record) != RECORD_KEYS
                or not isinstance(record["structural"], dict)
                or set(record["structural"]) != STRUCTURAL_KEYS):
            raise ProvenanceIndexError("record or structural allowlist mismatch")
        _bounded_string(record["class"], "compact classification class")
        _bounded_string(record["order"], "compact classification order", optional=True)
        _bounded_string(record["superfamily"], "compact classification superfamily")
        _bounded_string(record["confidence"], "compact classification confidence", maximum=32)
        _bounded_string(record["completeness"], "compact completeness", maximum=MAX_COMPLETENESS_STRING,
                        optional=True)
        terminal_repeats = record["structural"]["terminal_repeats"]
        tsds = record["structural"]["target_site_duplications"]
        if not isinstance(terminal_repeats, list) or len(terminal_repeats) > MAX_TERMINAL_REPEAT_ENTRIES:
            raise ProvenanceIndexError("compact terminal-repeat list exceeds its cap")
        if not isinstance(tsds, list) or len(tsds) > MAX_TSD_ENTRIES:
            raise ProvenanceIndexError("compact TSD list exceeds its cap")
        terminal_kinds = []
        for terminal_repeat in record["structural"]["terminal_repeats"]:
            if not isinstance(terminal_repeat, dict) or not set(terminal_repeat).issubset(TERMINAL_REPEAT_KEYS):
                raise ProvenanceIndexError("terminal-repeat allowlist mismatch")
            kind = _bounded_string(terminal_repeat.get("kind"), "compact terminal-repeat kind", maximum=32)
            if kind not in {"LTR", "LTR_advisory_subthreshold", "TIR"}:
                raise ProvenanceIndexError("unknown compact terminal-repeat kind")
            terminal_kinds.append(kind)
            _bounded_number(terminal_repeat.get("length"), "compact terminal-repeat length",
                            minimum=1, maximum=input_length, integer=True)
            _bounded_number(terminal_repeat.get("identity_pct"), "compact terminal-repeat identity",
                            minimum=0, maximum=100)
            if kind == "LTR":
                if not isinstance(terminal_repeat.get("length_is_lower_bound"), bool):
                    raise ProvenanceIndexError("compact LTR lower-bound flag must be boolean")
            if kind == "LTR_advisory_subthreshold":
                _bounded_number(terminal_repeat.get("threshold"), "compact advisory threshold",
                                minimum=0, maximum=100)
            if kind in {"LTR", "TIR"} and terminal_repeat.get("tsd_state") not in TSD_STATES:
                raise ProvenanceIndexError("unknown compact terminal-repeat tsd_state")
        if len(terminal_kinds) != len(set(terminal_kinds)):
            raise ProvenanceIndexError("duplicate compact terminal-repeat kind")
        if "LTR" in terminal_kinds and len(terminal_kinds) != 1:
            raise ProvenanceIndexError("compact accepted LTR cannot coexist with another terminal repeat")
        for tsd in record["structural"]["target_site_duplications"]:
            if not isinstance(tsd, dict) or set(tsd) != TSD_KEYS:
                raise ProvenanceIndexError("TSD allowlist mismatch")
            _bounded_number(tsd["length"], "compact TSD length",
                            minimum=parameters["tsd_min_len"], maximum=parameters["tsd_max_len"], integer=True)
            if not isinstance(tsd["matched_expected"], bool):
                raise ProvenanceIndexError("compact TSD matched_expected must be boolean")
            if tsd["congruence"] not in TSD_CONGRUENCE:
                raise ProvenanceIndexError("unknown compact TSD congruence")
    bad = sorted({key for key in _walk_keys(case) if key.lower() in FORBIDDEN_KEYS})
    if bad:
        raise ProvenanceIndexError(f"forbidden sequence-bearing keys in compact output: {bad}")
    encoded = json.dumps(case, sort_keys=True, ensure_ascii=True)
    # Long uninterrupted nucleotide strings are a second guard for a leak through
    # a value rather than a field name. SHA-256 values are hexadecimal, not DNA.
    if re.search(r"[ACGTUN]{25,}", encoded, flags=re.IGNORECASE):
        raise ProvenanceIndexError("long sequence-like value in compact output")


def _source_snapshot(paths: Iterable[Path]) -> dict[Path, tuple[int, int]]:
    return {path: (path.stat().st_size, path.stat().st_mtime_ns) for path in paths}


def _assert_sources_unchanged(snapshot: dict[Path, tuple[int, int]]) -> None:
    after = _source_snapshot(snapshot)
    changed = [str(path) for path, value in snapshot.items() if after[path] != value]
    if changed:
        raise ProvenanceIndexError(f"source manifests changed while indexing: {changed}")


def reject_duplicate_case_ids(cases: Iterable[dict[str, Any]]) -> None:
    case_ids = [case.get("case_id") for case in cases]
    duplicates = sorted(case_id for case_id, count in Counter(case_ids).items() if count > 1)
    if duplicates:
        raise ProvenanceIndexError(f"duplicate case IDs: {duplicates}")


def _derived_structural_values(cases: list[dict[str, Any]]) -> dict[str, Any]:
    terminal_repeats = [repeat for case in cases for record in case["records"]
                        for repeat in record["structural"]["terminal_repeats"]]
    tsds = [tsd for case in cases for record in case["records"]
            for tsd in record["structural"]["target_site_duplications"]]
    by_kind = Counter(repeat["kind"] for repeat in terminal_repeats)
    lengths = defaultdict(list)
    for repeat in terminal_repeats:
        lengths[repeat["kind"]].append(repeat["length"])
    return {
        "terminal_repeat_counts": dict(sorted(by_kind.items())),
        "terminal_repeat_length_ranges": {
            kind: {"min": min(values), "max": max(values)}
            for kind, values in sorted(lengths.items())
        },
        "ltr_lower_bound_count": sum(
            1 for repeat in terminal_repeats if repeat.get("length_is_lower_bound")
        ),
        "tsd_count": len(tsds),
        "tsd_length_range": None if not tsds else {
            "min": min(tsd["length"] for tsd in tsds),
            "max": max(tsd["length"] for tsd in tsds),
        },
        "tsd_matched_expected_count": sum(1 for tsd in tsds if tsd["matched_expected"]),
    }


def make_index(raw_groups: tuple[tuple[str, Path], ...] = RAW_GROUPS) -> tuple[dict[str, Any], dict[str, Any]]:
    """Read all current per-case manifests and return deterministic index and summary objects."""
    configured_corpora = CORPUS_FILES if raw_groups == RAW_GROUPS else None
    sources = validated_source_paths(raw_groups, configured_corpora)
    paths = [path for _, path in sources]
    snapshot = _source_snapshot(paths)
    cases = []
    for corpus, path in sources:
        raw_bytes = path.read_bytes()
        raw = json.loads(raw_bytes.decode("utf-8"))
        relative = path.relative_to(ROOT / "benchmarks" / "raw").as_posix()
        cases.append(compact_manifest(raw, corpus=corpus, source_manifest=relative,
                                      source_manifest_sha256=sha256_bytes(raw_bytes)))
    _assert_sources_unchanged(snapshot)
    reject_duplicate_case_ids(cases)
    cases.sort(key=lambda case: case["case_id"])
    index = {
        "schema_version": SCHEMA_VERSION,
        "description": "Compact sequence-free per-case provenance evidence extracted from raw TEagle manifests.",
        "hash_semantics": {
            "input_sha256": (
                "SHA-256 of the fetch-level FASTA text supplied by the benchmark runner, including its "
                "header and text representation."
            ),
            "sealed_input_sha256": (
                "SHA-256 of the analysed residue string recorded inside TEagle's content-addressed manifest."
            ),
            "relationship": (
                "These hashes cover different representations and are expected to differ; they must not be equated."
            ),
            "source_manifest_sha256": "SHA-256 of the complete private per-case source-manifest bytes.",
            "internal_manifest_sha256": "Recomputed and verified TEagle content-addressed manifest seal.",
        },
        "internal_seal_verification": {
            "algorithm": "SHA-256 of canonical JSON (sorted keys and compact separators).",
            "excluded_top_level_fields": sorted(SEAL_EXCLUDE_TOP),
            "excluded_input_fields": sorted(SEAL_EXCLUDE_INPUT),
        },
        "excluded_field_classes": [
            "raw nucleotide and amino-acid sequences, previews, motifs, and primers",
            "ORFs, domain hit payloads, and unbounded nested detector output",
            "structural coordinates and record spans",
            "free-text result explanations and source citations",
        ],
        "cases": cases,
    }
    summary = make_summary(index)
    return index, summary


def make_summary(index: dict[str, Any]) -> dict[str, Any]:
    cases = index["cases"]
    corpus_counts = dict(sorted(Counter(case["corpus"] for case in cases).items()))
    software = Counter((item["name"], item["version"]) for case in cases for item in case["software"])
    databases = Counter((item["name"], item["sha256"]) for case in cases for item in case["databases"])
    parameter_sets = Counter(json.dumps(case["parameters"], sort_keys=True, separators=(",", ":")) for case in cases)
    return {
        "schema_version": SCHEMA_VERSION,
        "case_count": len(cases),
        "corpus_counts": corpus_counts,
        "unique_accession_count": len({case["accession"] for case in cases}),
        "unique_input_sha256_count": len({case["input_sha256"] for case in cases}),
        "unique_sealed_input_sha256_count": len({case["sealed_input_sha256"] for case in cases}),
        "source_manifest_sha256_coverage": {
            "present": sum(bool(case["source_manifest_sha256"]) for case in cases), "total": len(cases),
        },
        "internal_manifest_sha256_coverage": {
            "present": sum(bool(case["internal_manifest_sha256"]) for case in cases), "total": len(cases),
        },
        "internal_manifest_sha256_verification": {
            "verified": len(cases), "total": len(cases),
        },
        "version_coverage": {
            "teagle_versions": dict(sorted(Counter(case["teagle_version"] for case in cases).items())),
            "provenance_schema_versions": dict(sorted(Counter(case["provenance_schema_version"] for case in cases).items())),
            "run_types": dict(sorted(Counter(case["run_type"] for case in cases).items())),
            "software": [
                {"name": name, "version": version, "case_count": count}
                for (name, version), count in sorted(software.items())
            ],
            "databases": [
                {"name": name, "sha256": digest, "case_count": count}
                for (name, digest), count in sorted(databases.items())
            ],
        },
        "parameter_coverage": {
            "complete_case_count": len(cases),
            "unique_sealed_parameter_sets": len(parameter_sets),
            "parameter_names": sorted({key for case in cases for key in case["parameters"]}),
        },
        "derived_structural_values": _derived_structural_values(cases),
        "historical_orphan_file_count": {
            "status": "historical_non_reconstructable_from_current_manifests",
            "note": "The compact extract indexes the current 129 broader and 43 holdout per-case manifests; it cannot establish a historical orphan-file count.",
        },
    }


def _atomic_json_write(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = (json.dumps(value, indent=2, sort_keys=True, ensure_ascii=True) + "\n").encode("utf-8")
    with NamedTemporaryFile("wb", delete=False, dir=path.parent, prefix=f".{path.name}.", suffix=".tmp") as handle:
        handle.write(payload)
        handle.flush()
        os.fsync(handle.fileno())
        temporary = Path(handle.name)
    os.replace(temporary, path)


def write_outputs(index: dict[str, Any], summary: dict[str, Any], output_dirs: Iterable[Path]) -> None:
    for directory in output_dirs:
        _atomic_json_write(directory / "provenance_index.json", index)
        _atomic_json_write(directory / "provenance_summary.json", summary)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", action="append", type=Path,
                        help="Write only to this directory; repeat for mirrors.")
    args = parser.parse_args(argv)
    index, summary = make_index()
    write_outputs(index, summary, args.output_dir or DEFAULT_OUTPUTS)
    print(f"wrote {summary['case_count']} compact provenance cases to "
          f"{', '.join(str(path) for path in (args.output_dir or DEFAULT_OUTPUTS))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
