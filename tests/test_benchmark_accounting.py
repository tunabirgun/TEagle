"""Regression guards for hash-level benchmark accounting; fixtures contain no sequence data."""
from __future__ import annotations

import copy
import glob
import json
import os
import shutil
import subprocess
import sys

import pytest


ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [os.path.join(ROOT, "benchmarks"), os.path.join(ROOT, "manuscript")]
import head_to_head  # noqa: E402
import make_table2  # noqa: E402
import score  # noqa: E402
VALUES_PATH = os.path.join(ROOT, "manuscript", "values.py")
if os.path.isfile(VALUES_PATH):
    import values  # noqa: E402
else:
    values = None


def _score_validators():
    validators = [make_table2.assert_score_integrity]
    if values is not None:
        validators.append(values.assert_score_integrity)
    return validators


def _require_values():
    if values is None:
        pytest.skip("manuscript working sources are not distributed in public source checkouts")


def _require_raw(*tools):
    missing = []
    for tool in tools:
        directory = os.path.join(ROOT, "benchmarks", "raw", tool)
        case_pattern = "*.default.json" if tool == "tesorter" else "[0-9]*.json"
        if (not os.path.isfile(os.path.join(directory, "_run.json"))
                or not glob.glob(os.path.join(directory, case_pattern))):
            missing.append(tool)
    if missing:
        pytest.skip("locally retained, sequence-bearing raw evidence is absent: " + ", ".join(missing))


def _require_head_to_head_raw():
    _require_values()
    _require_raw("teagle", "tesorter")
    if not glob.glob(os.path.join(ROOT, "benchmarks", "raw", "tesorter", "_work", "*.fa")):
        pytest.skip("locally retained comparator FASTA evidence is absent")


def _pair_row(accession="ACC", sha="a" * 64, order="LTR", observed="LTR", abstained=False):
    return {"accession": accession, "input_sha256": sha, "expected_class": "CLASS_I",
            "expected_order": order, "observed_order": observed, "abstained_at_order": abstained,
            "observed_te_class": "LTR/Copia", "confidence": "High", "single_element_input": True}


def _scores():
    with open(os.path.join(ROOT, "benchmarks", "raw", "scores.json"), encoding="utf-8") as fh:
        return json.load(fh)


def _comparison():
    with open(os.path.join(ROOT, "benchmarks", "raw", "head_to_head.json"), encoding="utf-8") as fh:
        return json.load(fh)


def _run_record(tool):
    with open(os.path.join(ROOT, "benchmarks", "raw", tool, "_run.json"), encoding="utf-8") as fh:
        return json.load(fh)


def test_duplicate_consistent_accession_rows_collapse():
    first = _pair_row()
    second = dict(first, panel="duplicate annotation")
    got = head_to_head.collapse_accession_rows([first, second], source="synthetic")
    assert list(got) == ["ACC"]
    assert got["ACC"]["input_sha256"] == first["input_sha256"]


def test_duplicate_conflict_fails_before_pairing():
    with pytest.raises(ValueError, match="ambiguous rows"):
        head_to_head.collapse_accession_rows([_pair_row(), _pair_row(observed="LINE")], source="synthetic")
    with pytest.raises(ValueError, match="conflicting order-level"):
        score.collapse_unique_inputs([
            {"accession": "A", "input_sha256": "b" * 64, "expected_order": "LTR",
             "observed_order": "LTR", "abstained_at_order": False},
            {"accession": "B", "input_sha256": "b" * 64, "expected_order": "LTR",
             "observed_order": "LINE", "abstained_at_order": False},
        ], "order")


def test_accession_universe_is_derived_from_retained_and_stratified_hashed_rows():
    scores = {"cases": [_pair_row("A")], "stratified_cases": [_pair_row("B", "b" * 64)]}
    got = head_to_head.successful_accessions(scores)
    assert set(got) == {"A", "B"}
    with pytest.raises(ValueError, match="input_sha256"):
        head_to_head.successful_accessions({"cases": [_pair_row("A", "")], "stratified_cases": []})


def test_current_unique_input_summary_and_values_bindings():
    _require_values()
    scores = _scores()
    unique = scores["overall_unique_inputs"]
    for level in ("class", "order"):
        assert unique[level]["n_gradable"] == 93
        assert unique[level]["n_answered"] == 66
        assert unique[level]["n_correct"] == 66
    assert scores["independence"]["element_only_order_withholding"] == {
        "n_gradable": 82, "n_withheld": 16, "n_answered": 66, "withholding_rate": 0.1951}
    gradable = [c for c in scores["cases"]
                if not (c["expected_class"] != "NOT_TE" and c["expected_order"] == "OTHER")]
    assert all(c["input_sha256"] for c in gradable)
    rendered = values.build()
    assert "corpus_species" not in rendered
    assert rendered["broad_distinct_inputs"] == len({c["input_sha256"] for c in gradable})
    assert rendered["broad_unique_order_answered"] == 66


def test_missing_comparator_accounting_and_fixed_mcnemar_tables():
    comparison = _comparison()
    assert comparison["historical_vs_current"] == {
        "historical_tesorter_corpus_rows": 136,
        "current_teagle_corpus_rows": 135,
        "corpus_row_count_mismatch": True,
    }
    expected = {"default": (64, 2, 9, 7, 0.0654), "matched": (61, 5, 9, 7, 0.424)}
    for mode, (both, te_only, ts_only, neither, p) in expected.items():
        result = comparison["modes"][mode]
        accounting = result["accounting"]
        assert accounting["current_successful_unique_accessions"] == 116
        assert accounting["comparator_outputs"] == 115
        assert accounting["paired"] == 82
        assert accounting["excluded_with_output"] == 33
        assert accounting["missing_comparator_output"] == {
            "n": 1,
            "cases": [{"accession": "NC_000001.10", "reason": "no TEsorter output for current successful accession"}],
        }
        assert accounting["extra_comparator_outputs"] == {"n": 0, "accessions": []}
        assert result["n_paired"] + accounting["excluded_with_output"] + accounting["missing_comparator_output"]["n"] == 116
        assert (result["both_correct"], result["teagle_only"], result["tesorter_only"], result["neither"]) == (
            both, te_only, ts_only, neither)
        assert result["mcnemar_p_two_sided_exact"] == p
        assert len({pair["input_sha256"] for pair in result["pairs"]}) == result["n_paired"]


def test_consumer_integrity_checks_reject_corrupted_aggregate():
    broken = copy.deepcopy(_scores())
    broken["overall"]["order"]["n_answered"] += 1
    for validator in _score_validators():
        with pytest.raises(ValueError, match="overall order"):
            validator(broken)


@pytest.mark.parametrize("family,path", [
    ("unique aggregate", ("overall_unique_inputs", "order", "n_correct")),
    ("withholding aggregate", ("independence", "element_only_order_withholding", "n_withheld")),
])
def test_consumers_reject_corrupted_unique_input_families(family, path):
    broken = copy.deepcopy(_scores())
    target = broken
    for key in path[:-1]:
        target = target[key]
    target[path[-1]] = 999
    for validator in _score_validators():
        with pytest.raises(ValueError, match="unique|withholding"):
            validator(broken)


def test_consumers_reject_missing_hash_and_duplicate_hash_conflict():
    missing = copy.deepcopy(_scores())
    missing["cases"][0]["input_sha256"] = ""
    conflicting = copy.deepcopy(_scores())
    duplicate_hash = next(sha for sha in {c["input_sha256"] for c in conflicting["cases"]}
                          if sum(c["input_sha256"] == sha for c in conflicting["cases"]) > 1)
    duplicate_rows = [c for c in conflicting["cases"] if c["input_sha256"] == duplicate_hash]
    duplicate_rows[1]["observed_order"] = "SYNTHETIC_CONFLICT"
    for validator in _score_validators():
        with pytest.raises(ValueError, match="input_sha256"):
            validator(missing)
        with pytest.raises(ValueError):
            validator(conflicting)


def test_score_cli_rejects_incomplete_raw_case_universe(tmp_path):
    _require_raw("teagle")
    rawdir = tmp_path / "raw"
    rawdir.mkdir()
    shutil.copy2(os.path.join(ROOT, "benchmarks", "raw", "teagle", "_run.json"), rawdir / "_run.json")
    one = sorted(glob.glob(os.path.join(ROOT, "benchmarks", "raw", "teagle", "[0-9]*.json")))[0]
    shutil.copy2(one, rawdir / "one.json")
    proc = subprocess.run(
        [sys.executable, os.path.join(ROOT, "benchmarks", "score.py"),
         "--corpus", os.path.join(ROOT, "benchmarks", "corpus.tsv"),
         "--rawdir", str(rawdir), "--out", str(tmp_path / "scores.json")],
        text=True, capture_output=True, check=False)
    assert proc.returncode != 0
    assert "raw-case ledger invalid" in proc.stderr
    assert not (tmp_path / "scores.json").exists()


def test_score_cli_rejects_duplicate_raw_case_key(tmp_path):
    _require_raw("teagle")
    rawdir = tmp_path / "raw"
    rawdir.mkdir()
    shutil.copy2(os.path.join(ROOT, "benchmarks", "raw", "teagle", "_run.json"), rawdir / "_run.json")
    one = sorted(glob.glob(os.path.join(ROOT, "benchmarks", "raw", "teagle", "[0-9]*.json")))[0]
    shutil.copy2(one, rawdir / "copy-a.json")
    shutil.copy2(one, rawdir / "copy-b.json")
    proc = subprocess.run(
        [sys.executable, os.path.join(ROOT, "benchmarks", "score.py"),
         "--corpus", os.path.join(ROOT, "benchmarks", "corpus.tsv"),
         "--rawdir", str(rawdir), "--out", str(tmp_path / "scores.json")],
        text=True, capture_output=True, check=False)
    assert proc.returncode != 0
    assert "duplicate raw case key" in proc.stderr
    assert not (tmp_path / "scores.json").exists()


def test_current_raw_case_universe_reconciles_exactly():
    _require_raw("teagle")
    rawdir = os.path.join(ROOT, "benchmarks", "raw", "teagle")
    files = sorted(glob.glob(os.path.join(rawdir, "[0-9]*.json")))
    records = score.validate_raw_case_universe(
        os.path.join(ROOT, "benchmarks", "corpus.tsv"), rawdir, files)
    assert len(records) == 129


@pytest.mark.parametrize("mutation,match", [
    (("corpus_sha256",), "corpus_sha256"),
    (("corpus_path",), "corpus_path"),
    (("corpus_rows",), "corpus_rows"),
    (("counts", "analysed"), "counts.analysed"),
    (("counts", "skipped_existing"), "skipped_existing"),
    (("failures", 0, "accession"), "ambiguous"),
])
def test_raw_case_guard_rejects_corrupt_run_ledger(tmp_path, mutation, match):
    _require_raw("teagle")
    rawdir = tmp_path / "raw"
    rawdir.mkdir()
    run = _run_record("teagle")
    if mutation == ("corpus_sha256",):
        run["corpus_sha256"] = "0" * 64
    elif mutation == ("corpus_path",):
        run["corpus_path"] = "wrong.tsv"
    elif mutation == ("corpus_rows",):
        run["corpus_rows"] = 999
    elif mutation == ("counts", "analysed"):
        run["counts"]["analysed"] = 999
    elif mutation == ("counts", "skipped_existing"):
        run["counts"]["skipped_existing"] = 999
    else:
        run["failures"][0]["accession"] = "AF123535"
    with open(rawdir / "_run.json", "w", encoding="utf-8") as fh:
        json.dump(run, fh)
    files = sorted(glob.glob(os.path.join(ROOT, "benchmarks", "raw", "teagle", "[0-9]*.json")))
    with pytest.raises(ValueError, match=match):
        score.validate_raw_case_universe(
            os.path.join(ROOT, "benchmarks", "corpus.tsv"), str(rawdir), files)


def _set_path(obj, path, value):
    target = obj
    for key in path[:-1]:
        target = target[key]
    target[path[-1]] = value


@pytest.mark.parametrize("path", [
    ("modes", "default", "n_paired"),
    ("modes", "default", "both_correct"),
    ("modes", "default", "mcnemar_discordant"),
    ("modes", "default", "mcnemar_p_two_sided_exact"),
    ("modes", "default", "teagle_accuracy"),
    ("modes", "default", "abstention_cost_teagle_silent_tesorter_right"),
    ("modes", "default", "abstention_vindicated_teagle_silent_tesorter_wrong"),
    ("modes", "default", "tesorter_no_call"),
    ("modes", "default", "order_errors", "with_negative_controls_denominator"),
    ("modes", "default", "order_errors", "tesorter_cases"),
    ("modes", "default", "excluded", "n_excluded"),
    ("modes", "default", "excluded", "by_reason", "negative_control"),
    ("modes", "default", "excluded", "negative_controls", "n"),
    ("modes", "default", "excluded", "cases", 0, "reason"),
    ("modes", "default", "input_identity", "n_checked"),
    ("modes", "default", "input_identity", "per_case", 0, "verdict"),
    ("modes", "default", "accounting", "current_successful_unique_accessions"),
    ("modes", "default", "accounting", "missing_comparator_output", "cases", 0, "reason"),
    ("historical_vs_current", "historical_tesorter_corpus_rows"),
    ("historical_tesorter_run", "corpus_rows"),
])
def test_head_to_head_consumer_rejects_corrupted_bound_family(path):
    _require_head_to_head_raw()
    broken = copy.deepcopy(_comparison())
    _set_path(broken, path, 999)
    with pytest.raises(ValueError, match="head_to_head.json"):
        values.assert_head_to_head_integrity(
            broken, _scores(), _run_record("tesorter"), _run_record("teagle"))


def test_head_to_head_consumer_rejects_corrupted_pair_truth():
    _require_head_to_head_raw()
    broken = copy.deepcopy(_comparison())
    broken["modes"]["default"]["pairs"][0]["teagle_correct"] = 999
    with pytest.raises(ValueError, match="head_to_head.json"):
        values.assert_head_to_head_integrity(
            broken, _scores(), _run_record("tesorter"), _run_record("teagle"))


@pytest.mark.parametrize("field", [
    "input_sha256", "accession", "expected", "teagle",
    "teagle_abstained", "teagle_confidence", "single_element_input",
])
def test_head_to_head_consumer_binds_pair_detail_to_scores(field):
    _require_head_to_head_raw()
    broken = copy.deepcopy(_comparison())
    pair = broken["modes"]["default"]["pairs"][0]
    if field == "input_sha256":
        pair[field] = ("0" if pair[field][0] != "0" else "1") + pair[field][1:]
    elif field == "accession":
        pair[field] = "SYNTHETIC_ACCESSION"
    elif field in ("expected", "teagle"):
        pair[field] = "TIR" if pair[field] != "TIR" else "LTR"
    elif field in ("teagle_abstained", "single_element_input"):
        pair[field] = not pair[field]
    else:
        pair[field] = "Synthetic confidence"
    with pytest.raises(ValueError, match="head_to_head.json"):
        values.assert_head_to_head_integrity(
            broken, _scores(), _run_record("tesorter"), _run_record("teagle"))


@pytest.mark.parametrize("field", [
    "input_sha256", "expected_order", "teagle", "teagle_abstained",
    "teagle_confidence", "single_element_input",
])
def test_head_to_head_consumer_binds_excluded_detail_to_scores(field):
    _require_head_to_head_raw()
    broken = copy.deepcopy(_comparison())
    case = broken["modes"]["default"]["excluded"]["cases"][0]
    if field == "input_sha256":
        case[field] = ("0" if case[field][0] != "0" else "1") + case[field][1:]
    elif field in ("expected_order", "teagle"):
        case[field] = "TIR" if case[field] != "TIR" else "LTR"
    elif field in ("teagle_abstained", "single_element_input"):
        case[field] = not case[field]
    else:
        case[field] = "Synthetic confidence"
    with pytest.raises(ValueError, match="head_to_head.json"):
        values.assert_head_to_head_integrity(
            broken, _scores(), _run_record("tesorter"), _run_record("teagle"))


@pytest.mark.parametrize("field", [
    "teagle_input_residue_sha256", "tesorter_input_residue_sha256", "verdict",
])
def test_head_to_head_consumer_recomputes_identity_hash_relation(field):
    _require_head_to_head_raw()
    broken = copy.deepcopy(_comparison())
    identity = broken["modes"]["default"]["input_identity"]["per_case"][0]
    if field == "verdict":
        identity[field] = "mismatch" if identity[field] != "mismatch" else "verified"
    else:
        identity[field] = ("0" if identity[field][0] != "0" else "1") + identity[field][1:]
    with pytest.raises(ValueError, match="head_to_head.json"):
        values.assert_head_to_head_integrity(
            broken, _scores(), _run_record("tesorter"), _run_record("teagle"))


def test_holdout_scores_match_current_schema_and_regenerate_exactly(tmp_path):
    canonical = os.path.join(ROOT, "benchmarks", "raw", "scores_holdout.json")
    reviewer = os.path.join(ROOT, "reviewer_package", "results", "scores_holdout.json")
    with open(canonical, "rb") as fh:
        canonical_bytes = fh.read()
    with open(reviewer, "rb") as fh:
        assert fh.read() == canonical_bytes
    stored = json.loads(canonical_bytes)
    for required in ("overall_unique_inputs", "independence", "stratified_cases"):
        assert required in stored

    rawdir = os.path.join(ROOT, "benchmarks", "raw", "teagle_holdout")
    run_path = os.path.join(rawdir, "_run.json")
    if not os.path.exists(run_path):
        pytest.skip("locally retained, sequence-bearing holdout outputs are absent")
    with open(run_path, encoding="utf-8") as stream:
        expected_cases = json.load(stream)["counts"]["analysed"]
    if len(glob.glob(os.path.join(rawdir, "[0-9]*.json"))) != expected_cases:
        pytest.skip("complete sequence-bearing holdout outputs are not retained")

    generated = tmp_path / "scores_holdout.json"
    proc = subprocess.run(
        [sys.executable, os.path.join(ROOT, "benchmarks", "score.py"),
         "--corpus", os.path.join(ROOT, "benchmarks", "corpus_holdout.tsv"),
         "--rawdir", rawdir,
         "--out", str(generated)], text=True, capture_output=True, check=False)
    assert proc.returncode == 0, proc.stderr
    assert generated.read_bytes() == canonical_bytes
