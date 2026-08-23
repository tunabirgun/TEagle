"""Focused contract tests for the compact, sequence-free provenance index."""
from __future__ import annotations

import copy
import hashlib
import importlib.util
import json
import re
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "benchmarks" / "make_provenance_index.py"
SPEC = importlib.util.spec_from_file_location("make_provenance_index", SCRIPT)
assert SPEC and SPEC.loader
provenance_index = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(provenance_index)
def _raw_case():
    """Return a sequence-free synthetic source manifest for contract tests.

    The publication tree intentionally excludes fetched and per-case sequence
    payloads.  This fixture exercises the raw-manifest validator without making
    those private working outputs a test dependency.
    """
    raw = {
        "accession": "SYNTHETIC",
        "input_length": 100,
        "input_sha256": "0123456789abcdef" * 4,
        "corpus_row": {
            "accession": "SYNTHETIC",
            "panel": "synthetic contract fixture",
            "record_scope": "element",
        },
        "result": {
            "provenance": {
                "schemaVersion": "0.1",
                "runType": "analysis",
                "teagleVersion": "3.8.0",
                "software": [{"name": "TEagle core", "version": "3.8.0"}],
                "databases": [{"name": "Synthetic profiles", "sha256": "b" * 64}],
                "parameters": {"tir_min_len": 10, "tsd_min_len": 2, "tsd_max_len": 8},
                "input": {"length": 100, "sha256": "fedcba9876543210" * 4},
            },
            "records": [{
                "classification": {
                    "class": "Class I",
                    "order": "LTR",
                    "superfamily": "Copia",
                    "confidence": "High",
                    "completeness": {"tier": "complete"},
                },
                "structural": [],
            }],
        },
    }
    provenance = raw["result"]["provenance"]
    provenance["manifestSha256"] = provenance_index.compute_internal_manifest_sha256(provenance)
    return raw


def _compact(raw):
    return provenance_index.compact_manifest(
        raw,
        corpus="holdout",
        source_manifest="teagle_holdout/001_SYNTHETIC.json",
        source_manifest_sha256="0" * 64,
    )


def _reseal(raw):
    provenance = raw["result"]["provenance"]
    provenance["manifestSha256"] = provenance_index.compute_internal_manifest_sha256(provenance)
    return raw


def _write_synthetic_run_group(directory, *, cases=1, recorded_cases=None):
    directory.mkdir(parents=True, exist_ok=True)
    recorded = cases if recorded_cases is None else recorded_cases
    run = {
        "corpus_rows": recorded,
        "counts": {"analysed": recorded, "skipped_existing": 0, "failed": 0},
    }
    (directory / "_run.json").write_text(json.dumps(run), encoding="utf-8")
    for number in range(cases):
        (directory / f"{number + 1:03d}_SYNTHETIC.json").write_text("{}", encoding="utf-8")
    return directory


def _all_source_paths():
    return sorted(
        path for _, directory in provenance_index.RAW_GROUPS
        for path in directory.glob("*.json") if path.name != "_run.json"
    )


def _stored_index_and_summary():
    canonical_dir = ROOT / "benchmarks" / "provenance"
    return (
        json.loads((canonical_dir / "provenance_index.json").read_text(encoding="utf-8")),
        json.loads((canonical_dir / "provenance_summary.json").read_text(encoding="utf-8")),
    )


def _private_sources_available():
    return all(directory.exists() and any(directory.glob("[0-9]*.json"))
               for _, directory in provenance_index.RAW_GROUPS)


def test_compact_index_is_safe_and_internally_consistent():
    index, summary = _stored_index_and_summary()
    assert provenance_index.make_summary(index) == summary
    assert index["schema_version"] == "1.1"
    assert summary["case_count"] == 172
    assert summary["corpus_counts"] == {"broader": 129, "holdout": 43}
    assert summary["unique_accession_count"] == 159
    assert summary["unique_input_sha256_count"] == 159
    assert summary["unique_sealed_input_sha256_count"] == 159
    assert summary["internal_manifest_sha256_verification"] == {"verified": 172, "total": 172}
    assert "expected to differ" in index["hash_semantics"]["relationship"]
    assert all(case["input_sha256"] != case["sealed_input_sha256"] for case in index["cases"])

    case_ids = [case["case_id"] for case in index["cases"]]
    assert case_ids == sorted(case_ids)
    assert len(case_ids) == len(set(case_ids))
    for case in index["cases"]:
        provenance_index.assert_safe_case(case)

    rendered = json.dumps(index, sort_keys=True)
    assert not re.search(r"[ACGTUN]{25,}", rendered, flags=re.IGNORECASE)
    forbidden = provenance_index.FORBIDDEN_KEYS
    assert not {key.lower() for key in provenance_index._walk_keys(index)} & forbidden

    # Exact extraction, source-byte hashes and seal recomputation require the
    # ignored maintainer manifests.  Exercise that stronger gate when those
    # private working outputs are present, without making them public inputs.
    if _private_sources_available():
        source_paths = _all_source_paths()
        before = {path: (path.stat().st_size, path.stat().st_mtime_ns) for path in source_paths}
        index_a, summary_a = provenance_index.make_index()
        index_b, summary_b = provenance_index.make_index()
        after = {path: (path.stat().st_size, path.stat().st_mtime_ns) for path in source_paths}
        assert before == after, "indexing must not alter source-manifest bytes or mtimes"
        assert index_a == index_b == index
        assert summary_a == summary_b == summary
        for case in index_a["cases"]:
            source = ROOT / "benchmarks" / "raw" / case["source_manifest"]
            assert hashlib.sha256(source.read_bytes()).hexdigest() == case["source_manifest_sha256"]


def test_summary_covers_versions_parameters_and_bounded_tir_tsd_values():
    index, summary = _stored_index_and_summary()
    assert provenance_index.make_summary(index) == summary
    assert summary["source_manifest_sha256_coverage"] == {"present": 172, "total": 172}
    assert summary["internal_manifest_sha256_coverage"] == {"present": 172, "total": 172}
    assert summary["internal_manifest_sha256_verification"] == {"verified": 172, "total": 172}
    assert summary["version_coverage"]["teagle_versions"] == {"3.8.0": 172}
    assert summary["version_coverage"]["provenance_schema_versions"] == {"0.1": 172}
    assert summary["parameter_coverage"]["complete_case_count"] == 172
    assert summary["parameter_coverage"]["unique_sealed_parameter_sets"] == 1
    derived = summary["derived_structural_values"]
    assert derived["terminal_repeat_counts"] == {
        "LTR": 62, "LTR_advisory_subthreshold": 7, "TIR": 43,
    }
    assert derived["tsd_count"] == 14
    assert derived["tsd_length_range"] == {"min": 2, "max": 8}
    assert derived["tsd_matched_expected_count"] == 0
    assert summary["historical_orphan_file_count"]["status"] == (
        "historical_non_reconstructable_from_current_manifests"
    )


def test_source_group_gate_rejects_missing_empty_and_count_mismatched_inputs(tmp_path):
    with pytest.raises(provenance_index.ProvenanceIndexError, match="directory is missing"):
        provenance_index.validated_source_paths((("broader", tmp_path / "missing"),), None)

    empty = _write_synthetic_run_group(tmp_path / "empty", cases=0, recorded_cases=1)
    with pytest.raises(provenance_index.ProvenanceIndexError, match="expected 1"):
        provenance_index.validated_source_paths((("broader", empty),), None)

    mismatch = _write_synthetic_run_group(tmp_path / "mismatch", cases=1, recorded_cases=2)
    with pytest.raises(provenance_index.ProvenanceIndexError, match="expected 2"):
        provenance_index.validated_source_paths((("broader", mismatch),), None)


@pytest.mark.parametrize("source_manifest", [
    "../001_SYNTHETIC.json",
    "teagle_holdout/../001_SYNTHETIC.json",
    "C:/Users/example/001_SYNTHETIC.json",
    r"teagle_holdout\001_SYNTHETIC.json",
    "teagle/001_SYNTHETIC.json",
    "teagle_holdout/nested/001_SYNTHETIC.json",
])
def test_source_manifest_policy_rejects_traversal_absolute_and_unexpected_paths(source_manifest):
    with pytest.raises(provenance_index.ProvenanceIndexError, match="source manifest"):
        provenance_index.compact_manifest(
            _raw_case(), corpus="holdout", source_manifest=source_manifest,
            source_manifest_sha256="0" * 64,
        )


@pytest.mark.parametrize("source_manifest", [
    "../001_SYNTHETIC.json", "C:/Users/example/001_SYNTHETIC.json",
    r"teagle_holdout\001_SYNTHETIC.json",
])
def test_safe_case_rejects_path_disclosures(source_manifest):
    compact = _compact(_raw_case())
    compact["source_manifest"] = source_manifest
    with pytest.raises(provenance_index.ProvenanceIndexError, match="source_manifest"):
        provenance_index.assert_safe_case(compact)


@pytest.mark.parametrize(
    "parameters",
    [
        {"nested": {"payload": 1}},
        {"_private": 1},
        {"payload": "A" * 1_000_000},
        {"text": "1"},
        {"flag": True},
        {"nan": float("nan")},
        {"infinite": float("inf")},
        {"x" * 65: 1},
        {f"p_{number}": number for number in range(65)},
    ],
)
def test_parameter_contract_rejects_nested_private_unbounded_or_nonnumeric_values(parameters):
    raw = _raw_case()
    raw["result"]["provenance"]["parameters"] = parameters
    with pytest.raises(provenance_index.ProvenanceIndexError):
        _compact(raw)


def test_length_hash_identity_and_internal_seal_are_verified_independently():
    raw = _raw_case()
    provenance = raw["result"]["provenance"]
    assert provenance_index.compute_internal_manifest_sha256(provenance) == provenance["manifestSha256"]
    case = _compact(raw)
    assert case["input_length"] == provenance["input"]["length"]
    assert case["input_sha256"] != case["sealed_input_sha256"]

    mismatch = _raw_case()
    mismatch["result"]["provenance"]["input"]["length"] += 1
    with pytest.raises(provenance_index.ProvenanceIndexError, match="analysed input_length"):
        _compact(mismatch)

    tampered = _raw_case()
    tampered["result"]["provenance"]["parameters"]["tir_min_len"] += 1
    with pytest.raises(provenance_index.ProvenanceIndexError, match="verification failed"):
        _compact(tampered)

    forged = _raw_case()
    forged["result"]["provenance"]["manifestSha256"] = "0" * 64
    with pytest.raises(provenance_index.ProvenanceIndexError, match="verification failed"):
        _compact(forged)

    volatile = _raw_case()
    volatile_provenance = volatile["result"]["provenance"]
    volatile_provenance["createdUtc"] = "2099-01-01T00:00:00+00:00"
    volatile_provenance["environment"] = {"machine": "different"}
    volatile_provenance["input"]["id"] = "different display header"
    assert _compact(volatile)["internal_manifest_sha256"] == provenance["manifestSha256"]


@pytest.mark.parametrize("conflicting", [False, True])
def test_duplicate_or_conflicting_software_and_database_names_fail(conflicting):
    raw = _raw_case()
    provenance = raw["result"]["provenance"]
    duplicate_software = copy.deepcopy(provenance["software"][0])
    if conflicting:
        duplicate_software["version"] = "different"
    provenance["software"].append(duplicate_software)
    with pytest.raises(provenance_index.ProvenanceIndexError, match="software"):
        _compact(raw)

    raw = _raw_case()
    provenance = raw["result"]["provenance"]
    duplicate_database = copy.deepcopy(provenance["databases"][0])
    if conflicting:
        duplicate_database["sha256"] = "f" * 64
    provenance["databases"].append(duplicate_database)
    with pytest.raises(provenance_index.ProvenanceIndexError, match="database"):
        _compact(raw)


def test_teagle_core_version_and_list_bounds_are_enforced():
    raw = _raw_case()
    raw["result"]["provenance"]["software"][-1]["version"] = "0.0.0"
    with pytest.raises(provenance_index.ProvenanceIndexError, match="equal teagleVersion"):
        _compact(raw)

    raw = _raw_case()
    raw["result"]["provenance"]["software"][0]["name"] = "x" * 129
    with pytest.raises(provenance_index.ProvenanceIndexError):
        _compact(raw)

    raw = _raw_case()
    raw["result"]["provenance"]["databases"][0]["version"] = "2026.1"
    case = _compact(_reseal(raw))
    assert case["databases"][0]["version"] == "2026.1"


def test_public_string_fields_are_bounded():
    mutations = [
        lambda raw: raw.update(accession="A" * 129),
        lambda raw: raw["corpus_row"].update(panel="p" * 65),
        lambda raw: raw["result"]["records"][0]["classification"].update(**{"class": "x" * 257}),
        lambda raw: raw["result"]["records"][0]["classification"].update(
            completeness={"tier": "x" * 513}
        ),
    ]
    for mutate in mutations:
        raw = _raw_case()
        mutate(raw)
        with pytest.raises(provenance_index.ProvenanceIndexError):
            _compact(raw)

    raw = _raw_case()
    raw["result"]["provenance"]["databases"] = [
        {"name": f"db-{number}", "sha256": f"{number:064x}"} for number in range(17)
    ]
    with pytest.raises(provenance_index.ProvenanceIndexError):
        _compact(raw)


def test_structural_measurements_are_bounded_unique_and_enum_checked():
    input_length = _raw_case()["input_length"]
    good_state = "absent (flanks present, no direct repeat found)"
    invalid_payloads = [
        [{"type": "TIR (terminal inverted repeat)", "tir_len": input_length + 1,
          "identity": 90, "tsd_state": good_state}],
        [{"type": "TIR (terminal inverted repeat)", "tir_len": 10,
          "identity": 101, "tsd_state": good_state}],
        [{"type": "TIR (terminal inverted repeat)", "tir_len": True,
          "identity": 90, "tsd_state": good_state}],
        [{"type": "TIR (terminal inverted repeat)", "tir_len": 10,
          "identity": 90, "tsd_state": "unknown"}],
        [{"type": "TSD (target-site duplication)", "length": 13,
          "matched_expected": False, "tsd_congruence": "not diagnostic"}],
        [{"type": "TSD (target-site duplication)", "length": 3,
          "matched_expected": False, "tsd_congruence": "unknown"}],
        [{"type": "TSD (target-site duplication)", "length": 3,
          "matched_expected": 0, "tsd_congruence": "not diagnostic"}],
        [{"type": "TIR (terminal inverted repeat)", "tir_len": 10,
          "identity": 90, "tsd_state": good_state}] * 2,
        [{"type": "TSD (target-site duplication)", "length": 3,
          "matched_expected": False, "tsd_congruence": "not diagnostic"}] * 2,
        [{"type": "ignored"}] * 33,
    ]
    for structural in invalid_payloads:
        raw = _raw_case()
        raw["result"]["records"][0]["structural"] = structural
        with pytest.raises(provenance_index.ProvenanceIndexError):
            _compact(raw)


def test_missing_required_identity_version_and_hash_fields_fail():
    paths = [
        ("raw", "accession"), ("raw", "input_length"), ("raw", "input_sha256"),
        ("corpus", "accession"), ("corpus", "panel"), ("corpus", "record_scope"),
        ("provenance", "schemaVersion"), ("provenance", "runType"),
        ("provenance", "teagleVersion"), ("provenance", "manifestSha256"),
        ("input", "length"), ("input", "sha256"),
    ]
    for location, key in paths:
        raw = _raw_case()
        containers = {
            "raw": raw,
            "corpus": raw["corpus_row"],
            "provenance": raw["result"]["provenance"],
            "input": raw["result"]["provenance"]["input"],
        }
        del containers[location][key]
        with pytest.raises(provenance_index.ProvenanceIndexError):
            _compact(raw)


def test_duplicate_case_ids_and_multiple_analysis_records_fail():
    with pytest.raises(provenance_index.ProvenanceIndexError, match="duplicate case IDs"):
        provenance_index.reject_duplicate_case_ids([{"case_id": "x"}, {"case_id": "x"}])
    raw = _raw_case()
    raw["result"]["records"].append(copy.deepcopy(raw["result"]["records"][0]))
    with pytest.raises(provenance_index.ProvenanceIndexError, match="exactly one"):
        _compact(raw)


def test_committed_index_and_reviewer_mirror_are_identical():
    canonical_dir = ROOT / "benchmarks" / "provenance"
    reviewer_dir = ROOT / "reviewer_package" / "results"
    assert (ROOT / "benchmarks" / "make_provenance_index.py").read_bytes() == (
        ROOT / "reviewer_package" / "scripts" / "make_provenance_index.py"
    ).read_bytes()
    for name in ("provenance_index.json", "provenance_summary.json"):
        canonical = canonical_dir / name
        reviewer = reviewer_dir / name
        assert canonical.read_bytes() == reviewer.read_bytes()
    stored_index = json.loads((canonical_dir / "provenance_index.json").read_text(encoding="utf-8"))
    stored_summary = json.loads((canonical_dir / "provenance_summary.json").read_text(encoding="utf-8"))
    assert provenance_index.make_summary(stored_index) == stored_summary
    for case in stored_index["cases"]:
        provenance_index.assert_safe_case(case)

    # Maintainer worktrees retain ignored private per-case outputs and can prove
    # exact regeneration.  Public source checkouts intentionally omit those
    # sequence-bearing payloads and therefore validate the compact evidence
    # intrinsically instead of treating an absent private input as a failure.
    if all(directory.exists() and any(directory.glob("[0-9]*.json"))
           for _, directory in provenance_index.RAW_GROUPS):
        generated_index, generated_summary = provenance_index.make_index()
        assert stored_index == generated_index
        assert stored_summary == generated_summary


def test_docs_use_transparent_organism_label_semantics_and_protocol_is_mirrored():
    for path in (ROOT / "benchmarks" / "CORPUS_SCOPE.md", ROOT / "reviewer_package" / "corpora" / "CORPUS_SCOPE.md",
                 ROOT / "reviewer_package" / "README.md"):
        text = path.read_text(encoding="utf-8")
        assert "50 species" not in text
        assert "56 distinct literal" in text
        assert "not a taxonomic species count" in text
    scope = (ROOT / "benchmarks" / "CORPUS_SCOPE.md").read_text(encoding="utf-8")
    scope_mirror = (ROOT / "reviewer_package" / "corpora" / "CORPUS_SCOPE.md").read_text(encoding="utf-8")
    assert scope == scope_mirror
    assert (ROOT / "benchmarks" / "CORPUS_SCOPE.md").read_bytes() == (
        ROOT / "reviewer_package" / "corpora" / "CORPUS_SCOPE.md"
    ).read_bytes()
    assert "Historical, superseded accounting" in scope
    assert "order accuracy 0.986" in scope
    assert "distinct analysed input hashes as the primary unit" in scope
    assert "66/66" in scope and "67/67" in scope
    assert "do not estimate performance in a target population" in scope

    reviewer = (ROOT / "reviewer_package" / "README.md").read_text(encoding="utf-8")
    assert "pending Phase 2" not in reviewer
    assert "Phase 3 has completed" not in reviewer
    assert "current local reviewer evidence set after independent reconciliation" in reviewer
    assert "Public deposition under a persistent identifier remains pending" not in reviewer
    assert "https://doi.org/10.5281/zenodo.21848814" in reviewer
    assert "proprietary" not in reviewer.lower()
    assert "AGPL-3.0-or-later" in reviewer
    reviewer_flat = " ".join(reviewer.split())
    for phrase in ("single-locus", "curated, non-random", "abstention denominator",
                   "not prospective wet-lab validation"):
        assert phrase in reviewer_flat
    assert "expected to differ and must not be equated" in reviewer

    canonical_path = ROOT / "benchmarks" / "protocol" / "BENCHMARK_SPEC_v1.md"
    canonical = canonical_path.read_text(encoding="utf-8")
    mirror_path = ROOT / "report" / "BENCHMARK_SPEC_v1.md"
    if mirror_path.exists():
        mirror = mirror_path.read_text(encoding="utf-8")
        assert canonical == mirror
        assert canonical_path.read_bytes() == mirror_path.read_bytes()
    assert "Amendment 4" in canonical
    assert "no preregistration timing claim" in canonical
    assert "reruns neither raw TEagle nor" in canonical
    assert "shortfall against the pre-registered plan" not in canonical
    assert "recorded *before* TEagle is run" not in canonical
    assert "cannot establish that timing" in canonical
    assert "This is the paper's contribution" not in canonical
    assert "Historical planned thesis (not supported by the current paired result)" in canonical
    assert "no comparator errors" in canonical
    for nonexistent in ("benchmarks/CORPUS_DEFECT.md", "benchmarks/run_tesorter.sh",
                        "benchmarks/analyse.py", "benchmarks/RESULTS.md"):
        assert nonexistent not in canonical
    for existing in ("benchmarks/CORPUS_SCOPE.md", "benchmarks/run_teagle.py",
                     "benchmarks/run_tesorter.py", "benchmarks/score.py",
                     "benchmarks/head_to_head.py", "benchmarks/figures.py",
                     "benchmarks/make_table1.py", "benchmarks/make_table2.py",
                     "benchmarks/make_table3.py", "benchmarks/make_tables.py",
                     "reviewer_package/README.md"):
        assert (ROOT / existing).exists(), existing
