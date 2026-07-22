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


def test_run_genome_pcr_requires_primers():
    with pytest.raises(BadRequest):
        engine.run_genome_pcr({"organism": "Homo sapiens"})


def test_run_genome_pcr_rejects_non_dna_primers():
    with pytest.raises(BadRequest):
        engine.run_genome_pcr({"fwd": "hello world", "rev": "TTTT" * 5, "organism": "Homo sapiens"})


def test_run_genome_pcr_requires_resolvable_assembly():
    with pytest.raises(BadRequest):                     # no organism and no accession -> cannot resolve a genome
        engine.run_genome_pcr({"fwd": "ACGT" * 5, "rev": "TTTT" * 5})


def test_run_genome_pcr_local_scan_seals_assembly(monkeypatch):
    raw = ">NC_000012.12:100+336 pair 237bp ACGTACGTACGTACGTACGT TTTTGGGGCCCCAAAATTTT\nACGT\n"
    seen = {}
    def fake_scan(acc, query, **k):
        seen["acc"] = acc; seen["query"] = query
        return {"ok": True, "raw": raw, "isPcr_version": "33x2", "target": "genome.2bit",
                "sha256": "abc123", "n_seqs": 24}
    monkeypatch.setattr(engine.wsl, "genome_scan", fake_scan)
    r = engine.run_genome_pcr({"fwd": "ACGTACGTACGTACGTACGT", "rev": "TTTTGGGGCCCCAAAATTTT",
                               "organism": "Homo sapiens"})
    assert r["ok"] and seen["acc"] == "GCF_000001405.40"
    assert r["assemblyAccession"] == "GCF_000001405.40" and r["scope"] == "whole-genome"
    assert r["n_amplicons"] == 1 and r["amplicons"][0]["on_target"] is False
    # a local scan IS sealed + reproducible: the assembly (accession + source-FASTA sha256) and isPcr
    # version are in the content-addressed manifest
    m = r["provenance"]
    assert m["runType"] == "genome-scan" and m["manifestSha256"]
    assert any(d["version"] == "GCF_000001405.40" and d["sha256"] == "abc123" for d in m["databases"])
    assert m["parameters"]["isPcr_version"] == "33x2"
    assert "fwdonly" in seen["query"] and "revonly" in seen["query"]   # all three isPcr query rows sent


def test_run_genome_pcr_seal_ignores_assembly_name(monkeypatch):
    # the human-readable assemblyName is a label, not sealing content: the SAME genome (accession + source
    # sha256) + primers must seal identically whether the name is populated (dropdown) or empty (bare accession)
    raw = ">NC_1:100+300 pair 201bp AAAAAAAAAAAAAAAAAAAA CCCCCCCCCCCCCCCCCCCC\n"
    monkeypatch.setattr(engine.wsl, "genome_scan",
                        lambda acc, query, **k: {"ok": True, "raw": raw, "isPcr_version": "33x2",
                                                 "target": "genome.2bit", "sha256": "SAME", "n_seqs": 1})
    base = {"fwd": "A" * 20, "rev": "C" * 20, "assemblyAccession": "GCF_9.9", "taxid": "1"}
    r1 = engine.run_genome_pcr({**base, "assemblyName": "Some Assembly v1"})
    r2 = engine.run_genome_pcr({**base, "assemblyName": ""})
    assert r1["provenance"]["manifestSha256"] == r2["provenance"]["manifestSha256"]


def test_run_genome_pcr_passes_through_need_prepare(monkeypatch):
    monkeypatch.setattr(engine.wsl, "genome_scan",
                        lambda acc, query, **k: {"ok": False, "error": "genome not prepared", "need_prepare": True})
    r = engine.run_genome_pcr({"fwd": "ACGT" * 5, "rev": "TTTT" * 5, "organism": "Homo sapiens"})
    assert r["ok"] is False and r.get("need_prepare") is True


def test_run_genome_prepare_resolves_accession(monkeypatch):
    seen = {}
    monkeypatch.setattr(engine.wsl, "genome_prepare",
                        lambda acc, name, **k: (seen.update(acc=acc, name=name),
                        {"ok": True, "accession": acc, "target": "genome.2bit", "sha256": "s", "n_seqs": 24, "bytes": 100})[1])
    r = engine.run_genome_prepare({"organism": "Saccharomyces cerevisiae"})
    assert r["ok"] and seen["acc"] == "GCF_000146045.2" and r["assemblyAccession"] == "GCF_000146045.2"


@pytest.mark.parametrize("bad", [float("nan"), float("inf"), float("-inf")])
def test_run_pcr_nonfinite_target_span_is_badrequest(bad):
    # a NaN/Inf target_span would otherwise seal into the reproducibility manifest
    with pytest.raises(BadRequest):
        engine.run_pcr({"sequence": "ACGT" * 60, "fwd": "ACGT", "rev": "ACGT", "target_span": [bad, 100]})


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


# --- run_fetch_coords: malformed transport is a BadRequest; an unresolvable request is a soft error ---
def test_run_fetch_coords_missing_regions_is_badrequest():
    with pytest.raises(BadRequest):
        engine.run_fetch_coords({"organism": "Homo sapiens"})


def test_run_fetch_coords_bad_strand_is_badrequest():
    with pytest.raises(BadRequest):
        engine.run_fetch_coords({"regions": "chr1:1-10", "strand": "sense", "organism": "Homo sapiens"})


def test_run_fetch_coords_nonstring_organism_is_badrequest():
    with pytest.raises(BadRequest):
        engine.run_fetch_coords({"regions": "chr1:1-10", "organism": 123})


def test_run_fetch_coords_no_organism_is_soft_error():
    r = engine.run_fetch_coords({"regions": "chr1:1-10", "organism": "", "customQuery": ""})
    assert r["ok"] is False and "organism" in r["error"].lower()


def test_run_fetch_coords_comma_only_region_is_soft_error():
    # a comma-only coordinate group must be a soft {ok:False}, not an uncaught ValueError/500
    r = engine.run_fetch_coords({"regions": "chr13:,-100", "organism": "Homo sapiens"})
    assert r["ok"] is False


def test_run_fetch_coords_custom_non_json_is_soft_error(monkeypatch):
    # a custom-organism resolve against a non-JSON Datasets body must be a soft {ok:False}, not a 500
    monkeypatch.setattr(engine.fetch, "_get", lambda *a, **k: "<html>rate limited</html>")
    r = engine.run_fetch_coords({"regions": "chr1:1-100", "organism": "", "customQuery": "Foo bar"})
    assert r["ok"] is False and r["error"]


def test_run_fetch_coords_curated_uses_pinned_assembly(monkeypatch):
    # a curated organism routes to its pinned assembly accession without touching resolve_assembly
    seen = {}
    def fake(regions, acc, name, org, taxid="", strand="+", refresh=False):
        seen.update(accession=acc, organism=org, strand=strand, taxid=taxid)
        return {"fasta": ">x\nACGT", "runType": "coordinate", "source": {}}
    monkeypatch.setattr(engine.fetch, "retrieve_coords", fake)
    r = engine.run_fetch_coords({"regions": "chr13:1-100", "strand": "-", "organism": "Homo sapiens"})
    assert r["ok"] is True
    assert seen["accession"] == engine.fetch.COORD_ASSEMBLIES["Homo sapiens"]["assemblyAccession"]
    assert seen["organism"] == "Homo sapiens" and seen["strand"] == "-"
    assert seen["taxid"] == "9606"                            # curated path threads the pinned taxid into the seal


def test_run_fetch_coords_custom_routes_through_resolve(monkeypatch):
    monkeypatch.setattr(engine.fetch, "resolve_assembly",
                        lambda q: {"organism": "Foo bar", "taxid": "9", "assemblyName": "FooAsm",
                                   "assemblyAccession": "GCF_999.1"})
    seen = {}
    monkeypatch.setattr(engine.fetch, "retrieve_coords",
                        lambda regions, acc, name, org, taxid="", strand="+", refresh=False:
                        (seen.update(accession=acc, taxid=taxid), {"runType": "coordinate", "source": {}})[1])
    r = engine.run_fetch_coords({"regions": "chr1:1-10", "organism": "", "customQuery": "GCF_999.1"})
    assert r["ok"] is True and seen["accession"] == "GCF_999.1" and seen["taxid"] == "9"


# --- RNA detection must read the sequence body, not a header that merely contains 'U' ---
def test_analyze_rna_flag_ignores_header_u():
    r = engine.run_analyze({"sequence": ">Human U2 snRNA gene\nACGTACGTACGTACGTACGT"})
    assert not any("RNA" in n for n in r["records"][0].get("notes", []))


def test_analyze_rna_flag_detects_body_u():
    r = engine.run_analyze({"sequence": ">x\nACGUACGUACGUACGU"})
    assert any("RNA" in n for n in r["records"][0].get("notes", []))


# --- a Primer3 environment fault must surface as a 500 (re-raised), not a user BadRequest ---
def test_run_primers_env_fault_is_not_badrequest(monkeypatch):
    def boom(*a, **k):
        raise RuntimeError("primer3 shared library unavailable")
    monkeypatch.setattr(engine.primers, "design_primers", boom)
    with pytest.raises(RuntimeError):
        engine.run_primers({"sequence": "ACGT" * 60})


# --- run_analyze end to end on a trivial sequence returns the expected shape ---
def test_run_analyze_shape():
    r = engine.run_analyze({"sequence": "ACGTACGTACGTACGTACGT"})
    assert "records" in r and "provenance" in r and "references" in r
    assert r["records"][0]["valid"] is True
    assert len(r["provenance"]["manifestSha256"]) == 64
