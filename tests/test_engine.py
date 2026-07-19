"""Engine-adapter validation tests. These guard the request-coercion logic that used to live
in server.py (_num / _clean_params / _require_nt) and is now the single source of truth shared
by the HTTP server and the native app. Each encodes a hard-won 400-not-500 fix."""
import math
import pytest
import engine
from engine import BadRequest


# --- non-finite / mistyped scalar params are a clean BadRequest, never a 500 (obs 4592/4625/4647).
# Non-scalar values (list/dict) pass through _clean_params and are rejected downstream at the _num call. ---
@pytest.mark.parametrize("bad", [float("nan"), float("inf"), float("-inf"), "abc"])
def test_clean_params_rejects_nonfinite_and_nonnumeric(bad):
    with pytest.raises(BadRequest):
        engine._clean_params({"max_mm": bad})


def test_clean_params_passes_valid_numbers():
    p = engine._clean_params({"max_mm": 2, "tp": "5", "prod_min": 70.0})
    assert p["max_mm"] == 2 and p["tp"] == 5 and p["prod_min"] == 70.0


def test_num_bool_is_rejected():
    with pytest.raises(BadRequest):
        engine._num(True, 0)


# --- non-nucleotide input (e.g. an accession pasted into the sequence box) -> BadRequest ---
def test_require_nt_rejects_accession_like_text():
    with pytest.raises(BadRequest):
        engine._require_nt("NC_003075.7")


def test_require_nt_accepts_real_sequence():
    engine._require_nt("ACGTACGTACGTNNNN")          # no raise


def test_require_nt_rejects_empty():
    with pytest.raises(BadRequest):
        engine._require_nt("")


# --- run_pcr: the float-tp / NaN-max_mm crash class is a 400 through the real entry point ---
def test_run_pcr_nonfinite_param_is_badrequest():
    body = {"sequence": "ACGT" * 40, "fwd": "ACGTACGTACGTACGTAC", "rev": "TTGGTTGGTTGGTTGGTT",
            "params": {"max_mm": float("nan")}}
    with pytest.raises(BadRequest):
        engine.run_pcr(body)


def test_run_pcr_requires_primers():
    with pytest.raises(BadRequest):
        engine.run_pcr({"sequence": "ACGT" * 40})


def test_run_pcr_nonstring_background_is_badrequest():
    body = {"sequence": "ACGT" * 40, "fwd": "ACGTACGTACGTACGTAC", "rev": "TTGGTTGGTTGGTTGGTT",
            "background": 5}
    with pytest.raises(BadRequest):
        engine.run_pcr(body)


def test_run_pcr_malformed_target_span_is_badrequest():
    body = {"sequence": "ACGT" * 40, "fwd": "ACGTACGTACGTACGTAC", "rev": "TTGGTTGGTTGGTTGGTT",
            "target_span": "whoops"}
    with pytest.raises(BadRequest):
        engine.run_pcr(body)


# --- run_annotate: non-string species must not reach wsl (was a 500, obs 4640) ---
def test_run_annotate_nonstring_species_is_badrequest():
    with pytest.raises(BadRequest):
        engine.run_annotate({"sequence": ">x\nACGTACGTACGT", "species": 5})


def test_run_annotate_empty_sequence_is_badrequest():
    with pytest.raises(BadRequest):
        engine.run_annotate({"sequence": ""})


# --- run_primers: included region must be a well-formed pair ---
def test_run_primers_malformed_included_is_badrequest():
    with pytest.raises(BadRequest):
        engine.run_primers({"sequence": "ACGT" * 60, "included": [1, 2, 3]})


# --- run_fetch: a non-string accession is a BadRequest, a bad accession is a soft {ok:False} ---
def test_run_fetch_nonstring_accession_is_badrequest():
    with pytest.raises(BadRequest):
        engine.run_fetch({"accession": 12345})


# --- run_analyze end to end on a trivial sequence returns the expected shape ---
def test_run_analyze_shape():
    r = engine.run_analyze({"sequence": "ACGTACGTACGTACGTACGT"})
    assert "records" in r and "provenance" in r and "references" in r
    assert r["records"][0]["valid"] is True
    assert len(r["provenance"]["manifestSha256"]) == 64
