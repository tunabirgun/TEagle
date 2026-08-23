"""Engine-adapter validation tests. These guard the request-coercion logic that used to live
in server.py (_num / _clean_params / _require_nt) and is now the single source of truth shared
by the HTTP server and the native app. Each encodes a hard-won 400-not-500 fix."""
import copy, hashlib, json, math
import pytest
import engine
from engine import BadRequest
from teagle_core import provenance


def _fetched_source(accession, sequence, **updates):
    raw = {
        "accession": accession, "taxid": "7227", "source": "NCBI nuccore",
        "endpoint": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
        "retrievedUtc": "2026-08-21T00:00:00+00:00",
        "sourceUrl": f"https://www.ncbi.nlm.nih.gov/nuccore/{accession}", "length": len(sequence),
    }
    raw.update(updates)
    return provenance.normalize_source_context(raw, [(accession, sequence)])


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


def test_run_genome_pcr_returns_offtarget_summary(monkeypatch):
    raw = "\n".join([
        ">NC_1:100+300 pair 201bp AAAAAAAAAAAAAAAAAAAA CCCCCCCCCCCCCCCCCCCC",
        ">NC_1:900+1100 pair 201bp AAAAAAAAAAAAAAAAAAAA CCCCCCCCCCCCCCCCCCCC",
        ">NC_2:10+210 pair 201bp AAAAAAAAAAAAAAAAAAAA CCCCCCCCCCCCCCCCCCCC",
        ">NC_1:5+55 fwdonly 51bp AAAAAAAAAAAAAAAAAAAA AAAAAAAAAAAAAAAAAAAA",
    ]) + "\n"
    monkeypatch.setattr(engine.wsl, "genome_scan",
                        lambda acc, query, **k: {"ok": True, "raw": raw, "isPcr_version": "33x2",
                                                 "target": "genome.2bit", "sha256": "s", "n_seqs": 5})
    r = engine.run_genome_pcr({"fwd": "A" * 20, "rev": "C" * 20, "organism": "Homo sapiens"})
    s = r["summary"]
    assert s["n_pair"] == 3 and s["n_single"] == 1 and s["n_sources"] == 2
    assert s["n_pair"] + s["n_single"] == r["n_amplicons"]
    assert dict(s["per_source"])["NC_1"] == 2 and s["tier"] == "low-copy"


def test_run_genome_pcr_summary_not_in_seal(monkeypatch):
    # the off-target summary is DERIVED from the discovered products; folding it into the seal would make the
    # SAME primers+genome+params seal differently as the product set changes. Two scans with different product
    # sets but identical inputs MUST produce an identical manifest seal.
    base = {"fwd": "A" * 20, "rev": "C" * 20, "assemblyAccession": "GCF_9.9", "taxid": "1"}
    one = ">NC_1:100+300 pair 201bp AAAAAAAAAAAAAAAAAAAA CCCCCCCCCCCCCCCCCCCC\n"
    many = "\n".join(f">NC_{i}:100+300 pair 201bp AAAAAAAAAAAAAAAAAAAA CCCCCCCCCCCCCCCCCCCC"
                     for i in range(1, 8)) + "\n"
    def scan_of(raw):
        return lambda acc, query, **k: {"ok": True, "raw": raw, "isPcr_version": "33x2",
                                        "target": "genome.2bit", "sha256": "SAME", "n_seqs": 3}
    monkeypatch.setattr(engine.wsl, "genome_scan", scan_of(one))
    r1 = engine.run_genome_pcr(dict(base))
    monkeypatch.setattr(engine.wsl, "genome_scan", scan_of(many))
    r2 = engine.run_genome_pcr(dict(base))
    assert r1["summary"]["n_pair"] == 1 and r2["summary"]["n_pair"] == 7       # different product sets
    assert r1["provenance"]["manifestSha256"] == r2["provenance"]["manifestSha256"]   # identical seal


def test_run_genome_pcr_marks_on_target_from_design_locus(monkeypatch):
    # a product overlapping the design locus is the on-target; the sibling elsewhere is off-target
    raw = "\n".join([
        ">NC_1:1000+1240 pair 241bp AAAAAAAAAAAAAAAAAAAA CCCCCCCCCCCCCCCCCCCC",   # overlaps locus -> on-target
        ">NC_1:5000+5240 pair 241bp AAAAAAAAAAAAAAAAAAAA CCCCCCCCCCCCCCCCCCCC",   # elsewhere -> off-target
    ]) + "\n"
    monkeypatch.setattr(engine.wsl, "genome_scan",
                        lambda acc, query, **k: {"ok": True, "raw": raw, "isPcr_version": "33x2",
                                                 "target": "genome.2bit", "sha256": "s", "n_seqs": 2})
    r = engine.run_genome_pcr({"fwd": "A" * 20, "rev": "C" * 20, "organism": "Homo sapiens",
                               "design_locus": {"accession": "NC_1", "start": 950, "stop": 1300}})
    assert r["has_locus"] and r["n_on_target"] == 1 and r["summary"]["n_on"] == 1 and r["summary"]["n_off"] == 1
    on = [a for a in r["amplicons"] if a.get("on_target")]
    assert len(on) == 1 and on[0]["start"] == 1000


def test_run_genome_pcr_no_locus_is_neutral(monkeypatch):
    raw = ">NC_1:1000+1240 pair 241bp AAAAAAAAAAAAAAAAAAAA CCCCCCCCCCCCCCCCCCCC\n"
    monkeypatch.setattr(engine.wsl, "genome_scan",
                        lambda acc, query, **k: {"ok": True, "raw": raw, "isPcr_version": "33x2",
                                                 "target": "genome.2bit", "sha256": "s", "n_seqs": 1})
    r = engine.run_genome_pcr({"fwd": "A" * 20, "rev": "C" * 20, "organism": "Homo sapiens"})   # no design_locus
    assert r["has_locus"] is False and r["n_on_target"] == 0 and not any(a.get("on_target") for a in r["amplicons"])


def test_run_genome_pcr_design_locus_not_in_seal(monkeypatch):
    # the on/off labeling from the design locus is DERIVED, not sealed: same primers+genome seal identically
    raw = ">NC_1:1000+1240 pair 241bp AAAAAAAAAAAAAAAAAAAA CCCCCCCCCCCCCCCCCCCC\n"
    monkeypatch.setattr(engine.wsl, "genome_scan",
                        lambda acc, query, **k: {"ok": True, "raw": raw, "isPcr_version": "33x2",
                                                 "target": "genome.2bit", "sha256": "SAME", "n_seqs": 1})
    base = {"fwd": "A" * 20, "rev": "C" * 20, "assemblyAccession": "GCF_9.9", "taxid": "1"}
    r1 = engine.run_genome_pcr(dict(base))                                                     # no locus
    r2 = engine.run_genome_pcr({**base, "design_locus": {"accession": "NC_1", "start": 950, "stop": 1300}})
    assert r2["n_on_target"] == 1 and r1["n_on_target"] == 0                                   # different labeling
    assert r1["provenance"]["manifestSha256"] == r2["provenance"]["manifestSha256"]            # identical seal


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
def test_run_pcr_nonfinite_intended_span_is_badrequest(bad):
    with pytest.raises(BadRequest):
        engine.run_pcr({"sequence": "ACGT" * 60, "fwd": "ACGT", "rev": "ACGT",
                        "intended_product_span": [bad, 100]})


def test_run_pcr_rejects_ambiguous_legacy_target_span():
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
        engine.run_primers({"sequence": "ACGT" * 60, "included_region": [1, 2, 3]})


def _design_bound_pcr(monkeypatch):
    fwd, rev = "GACTGACTGTCAGTCAGGCT", "TTGGCCATTGGCACTGGCAT"
    seq = "C" * 30 + fwd + "A" * 130 + engine.sequtil.reverse_complement(rev) + "C" * 30
    left, right = 30, 200
    target = [80, 20]
    candidate = {
        "id": "P1", "left_seq": fwd, "right_seq": rev,
        "left_pos": [left, left + len(fwd)], "right_pos": [right - len(rev), right],
        "left_tm": 60.0, "right_tm": 60.0, "left_gc": 50.0, "right_gc": 50.0,
        "product_size": right - left, "penalty": 1.0,
        "primer3_target": target, "included_region": None,
        "intended_product_span": [left, right], "intended_endpoint_tolerance_bp": 0,
    }
    monkeypatch.setattr(engine.primers, "design_primers", lambda *a, **k: {
        "candidates": [copy.deepcopy(candidate)], "explain_left": "", "explain_right": "", "explain_pair": ""})
    monkeypatch.setattr(engine, "_prime_vienna", lambda c: None)
    monkeypatch.setattr(engine.oligoqc, "qc_pair", lambda *a: {"ok": False, "error": "not run"})
    design_source = {"accession": "DESIGN1", "source": "design fixture"}
    designed = engine.run_primers({"sequence": ">record\n" + seq, "primer3_target": target,
                                   "source": design_source})
    c = designed["candidates"][0]
    body = {
        "sequence": ">record\n" + seq, "fwd": c["left_seq"], "rev": c["right_seq"],
        "primer3_target": c["primer3_target"], "included_region": c["included_region"],
        "intended_product_span": c["intended_product_span"],
        "intended_endpoint_tolerance_bp": c["intended_endpoint_tolerance_bp"],
        "design_provenance": c["design_provenance"], "design_sequence": ">record\n" + seq,
        "source": design_source,
        "params": {"max_mm": 0, "tp": 5,
        "prod_min": 70, "prod_max": 500},
    }
    return body


def test_design_and_pcr_manifests_bind_distinct_geometry(monkeypatch):
    body = _design_bound_pcr(monkeypatch)
    result = engine.run_pcr(body)
    assert result["status"] == "complete"
    assert result["mechanismCount"] >= result["uniqueProductCount"] >= 1
    assert any(a["on_target"] is True for a in result["amplicons"])
    params = result["provenance"]["parameters"]
    assert params["primer3_target"] == body["primer3_target"]
    assert params["intended_product_span"] == body["intended_product_span"]
    assert params["included_region"] is None
    assert params["intended_endpoint_tolerance_bp"] == 0
    assert params["design_manifest_sha256"] == body["design_provenance"]["manifestSha256"]
    from teagle_core.gff3 import public_record_identity
    expected = public_record_identity("record")["id"]
    assert result["bindings"]["design"]["id"] == expected
    assert result["bindings"]["assay"]["id"] == expected
    assert "accession" not in result["bindings"]["design"]
    assert "accession" not in result["bindings"]["assay"]


@pytest.mark.parametrize("field,value", [
    ("primer3_target", [81, 20]),
    ("included_region", [0, 210]),
    ("intended_product_span", [31, 200]),
    ("intended_endpoint_tolerance_bp", 1),
    ("fwd", "GACTGACTGTCAGTCAGGCA"),
])
def test_pcr_rejects_candidate_mutations_against_design_seal(monkeypatch, field, value):
    body = _design_bound_pcr(monkeypatch)
    body[field] = value
    with pytest.raises(BadRequest):
        engine.run_pcr(body)


def test_pcr_rejects_manifest_parameter_and_seal_mutations(monkeypatch):
    body = _design_bound_pcr(monkeypatch)
    body["design_provenance"] = copy.deepcopy(body["design_provenance"])
    body["design_provenance"]["parameters"]["candidate"]["product_size"] += 1
    with pytest.raises(BadRequest):
        engine.run_pcr(body)


@pytest.mark.parametrize("mutation", [
    ("input", "id", "X01005"),
    ("input", "length", 999),
    ("input", "sha256", "0" * 64),
    ("candidate", "selected_record", {"index": 0, "id": "other", "length": 230, "sha256": "0" * 64}),
    ("candidate", "left_pos", [31, 51]),
    ("candidate", "right_pos", [179, 199]),
    ("candidate", "product_size", 169),
    ("candidate", "intended_product_span", [31, 200]),
    ("candidate", "intended_contributors", {"left": "R", "right": "F", "orientation": "R+F"}),
])
def test_pcr_rejects_resealed_internal_design_contradictions(monkeypatch, mutation):
    body = _design_bound_pcr(monkeypatch)
    manifest = copy.deepcopy(body["design_provenance"])
    scope, field, value = mutation
    target = manifest["input"] if scope == "input" else manifest["parameters"]["candidate"]
    target[field] = value
    manifest["manifestSha256"] = hashlib.sha256(provenance.canonical_manifest_bytes(manifest)).hexdigest()
    if scope == "input" and field == "id":
        assert not provenance.verify_manifest(manifest), "a resealed compact/full identity mismatch must fail semantically"
    else:
        assert provenance.verify_manifest(manifest), "mutation must pass the cryptographic seal before semantic rejection"
    body["design_provenance"] = manifest
    with pytest.raises(BadRequest):
        engine.run_pcr(body)
    body = _design_bound_pcr(monkeypatch)
    body["design_provenance"] = copy.deepcopy(body["design_provenance"])
    body["design_provenance"]["manifestSha256"] = "0" * 64
    with pytest.raises(BadRequest):
        engine.run_pcr(body)


def test_pcr_rejects_false_assessed_claim_for_different_record_identity(monkeypatch):
    body = _design_bound_pcr(monkeypatch)
    body["sequence"] = body["sequence"].replace(">record", ">different-record", 1)
    body["design_assessment"] = "assess"
    with pytest.raises(BadRequest):
        engine.run_pcr(body)


@pytest.mark.parametrize("field,value", [
    ("primer3_target", [250, 20]),
    ("included_region", [100, 50]),
    ("intended_product_span", [200, 100]),
    ("intended_endpoint_tolerance_bp", -1),
])
def test_pcr_rejects_invalid_or_crossed_geometry(field, value):
    body = {"sequence": ">x\n" + "ACGT" * 70, "fwd": "ACGT" * 5, "rev": "TGCA" * 5,
            field: value}
    with pytest.raises(BadRequest):
        engine.run_pcr(body)


def test_pcr_completion_status_distinguishes_zero_and_capped(monkeypatch):
    zero = engine.run_pcr({"sequence": "A" * 200, "fwd": "C" * 18, "rev": "G" * 18,
                           "params": {"max_mm": 0}})
    assert zero["status"] == "complete_zero_products"
    monkeypatch.setattr(engine.primers, "in_silico_pcr", lambda *a, stats=None, **k:
                        (stats.update(amplicons_capped=True, sites_capped=False,
                                      max_amplicons=4000, max_sites=4000) or []))
    capped = engine.run_pcr({"sequence": "A" * 200, "fwd": "A" * 18, "rev": "T" * 18})
    assert capped["status"] == "partial_capped" and capped["truncated"] is True


@pytest.mark.parametrize("params", [
    {"max_mm": -2}, {"max_mm": 19}, {"tp": -1}, {"tp": 2.9}, {"tp": 19},
    {"min_anneal": -7}, {"min_anneal": 17},
    {"prod_min": 500, "prod_max": 100}, {"prod_min": 70.5}, {"max_mm": True},
])
def test_pcr_rejects_silent_numeric_normalization(params):
    with pytest.raises(BadRequest):
        engine.run_pcr({"sequence": "A" * 200, "fwd": "A" * 18, "rev": "T" * 18,
                        "params": params})


def test_pcr_rejects_primers_shorter_than_sealed_annealing_floor():
    with pytest.raises(BadRequest):
        engine.run_pcr({"sequence": "A" * 200, "fwd": "A" * 17, "rev": "T" * 17})


def test_default_primer_manifest_seals_all_effective_primer3_parameters(monkeypatch):
    body = _design_bound_pcr(monkeypatch)
    manifest = body["design_provenance"]
    for name, value in engine.primers.DESIGN_DEFAULTS.items():
        assert manifest["parameters"][name] == value


def test_pcr_reports_withheld_geometry_separately_from_search_completion(monkeypatch):
    # A withheld status is derived from a real sealed design whose record differs from the assay record;
    # a caller cannot assert withholding without preserving that design identity.
    body = _design_bound_pcr(monkeypatch)
    body["sequence"] = body["sequence"].replace(">record", ">assay-record", 1)
    body["design_assessment"] = "withheld_template_mismatch"
    body["source"] = {"accession": "ASSAY1", "source": "local fixture"}
    result = engine.run_pcr(body)
    assert result["criteria"]["intended_product_status"] == "withheld"
    assert result["criteria"]["target_coverage_status"] == "withheld"
    from teagle_core.gff3 import public_record_identity
    assert result["bindings"]["design"]["id"] == public_record_identity("record")["id"]
    assert result["bindings"]["assay"]["id"] == public_record_identity("assay-record")["id"]
    assert "accession" not in result["bindings"]["assay"]
    assert all(a["matches_intended_product"] is None and a["covers_primer3_target"] is None
               and a["on_target"] is None for a in result["amplicons"])


def test_custom_background_preserves_identity_but_never_inherits_intended_call(monkeypatch):
    body = _design_bound_pcr(monkeypatch)
    design_seq = body["design_sequence"].split("\n", 1)[1]
    body["background"] = ">background-record\n" + design_seq
    result = engine.run_pcr(body)
    assay = [a for a in result["amplicons"] if a["source"].endswith("(input)")]
    background = [a for a in result["amplicons"] if a["source"].endswith("(custom background)")]
    assert assay and any(a["on_target"] is True for a in assay)
    assert background and all(a["matches_intended_product"] is None and
                              a["covers_primer3_target"] is None and a["on_target"] is None
                              for a in background)
    from teagle_core.gff3 import public_record_identity
    expected_background = public_record_identity("background-record")
    assert result["bindings"]["backgrounds"][0]["id"] == expected_background["id"]
    assert result["bindings"]["backgrounds"][0]["idSha256"] == expected_background["idSha256"]
    params = result["provenance"]["parameters"]
    assert params["assay_input"]["id"] == public_record_identity("record")["id"]
    assert params["design_input"]["id"] == public_record_identity("record")["id"]
    assert params["search_inputs"][1]["role"] == "background"


# --- run_fetch: a non-string accession is a BadRequest, a bad accession is a soft {ok:False} ---
def test_run_fetch_nonstring_accession_is_badrequest():
    with pytest.raises(BadRequest):
        engine.run_fetch({"accession": 12345})


def test_raw_native_fetch_normalizes_to_serializable_bound_source_context(monkeypatch):
    accession = "M11240.1"
    sequence = "ACGT" * 10
    raw = {
        "accession": accession, "taxid": "7227", "source": "NCBI nuccore",
        "endpoint": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
        "retrievedUtc": "2026-08-21T00:00:00+00:00",
        "sourceUrl": f"https://www.ncbi.nlm.nih.gov/nuccore/{accession}",
        "length": len(sequence), "fasta": f">{accession}\n{sequence}",
    }
    monkeypatch.setattr(engine.fetch, "retrieve", lambda *a, **k: copy.deepcopy(raw))
    result = engine.run_fetch({"accession": accession})
    wire = json.loads(json.dumps(result["sourceContext"]))
    assert provenance.validate_source_context(wire, [(accession, sequence)]) == wire
    assert wire["records"][0]["sha256"] == hashlib.sha256(sequence.encode()).hexdigest()


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
    assert r["resultSchemaVersion"] == "1.0" and r["status"] == "complete"
    assert r["records"][0]["processing"]["structural"]["status"] == "complete"


def test_analysis_manifest_binds_ordered_record_identity_and_boundaries():
    a = engine.run_analyze({"sequence": ">whole\n" + "ACGT" * 6})
    b = engine.run_analyze({"sequence": ">left\n" + "ACGT" * 3 + "\n>right\n" + "ACGT" * 3})
    assert a["provenance"]["input"]["sha256"] == b["provenance"]["input"]["sha256"]
    assert a["provenance"]["manifestSha256"] != b["provenance"]["manifestSha256"]
    from teagle_core.gff3 import public_record_identity
    assert [x["id"] for x in b["provenance"]["input"]["records"]] == [
        public_record_identity("left")["id"], public_record_identity("right")["id"]]


def test_analysis_stage_failure_is_structured_and_withholds_structural_absence(monkeypatch):
    monkeypatch.setattr(engine.structural, "detect_all", lambda seq: (_ for _ in ()).throw(RuntimeError("boom")))
    r = engine.run_analyze({"sequence": ">x\n" + "ACGT" * 30})
    rec = r["records"][0]
    assert r["status"] == "complete_with_stage_failures"
    assert rec["processing"]["structural"] == {"status": "failed", "errorType": "RuntimeError"}
    assert rec["processing"]["classification"]["status"] == "withheld"
    if rec["classification"]:
        assert rec["classification"]["structural_unavailable"] is True
        assert "structural scan did not run" in " ".join(rec["classification"]["evidence"])
    assert engine.validate_analysis_result_document(engine.analysis_result_document(r))


def test_classification_failure_leaves_retroviral_stage_not_run_and_validates(monkeypatch):
    monkeypatch.setattr(engine.classify, "classify", lambda *a, **k: (_ for _ in ()).throw(RuntimeError("boom")))
    result = engine.run_analyze({"sequence": ">x\n" + "ACGT" * 30})
    rec = result["records"][0]
    assert rec["processing"]["classification"] == {"status": "failed", "errorType": "RuntimeError"}
    assert rec["processing"]["retroviral"] == {"status": "not_run"}
    assert rec["retroviral"] is None
    assert engine.validate_analysis_result_document(engine.analysis_result_document(result))


def test_analysis_status_distinguishes_zero_partial_and_not_run(monkeypatch):
    monkeypatch.setattr(engine.structural, "detect_all", lambda seq: [])
    def fake_orfs(seq, min_aa=40):
        return [{"start": i, "end": i + 3} for i in range(engine.domains.MAX_ORFS_SCANNED + 1)]
    monkeypatch.setattr(engine.sequtil, "find_orfs", fake_orfs)
    def fake_domains(seq, evalue=1e-3, max_orfs=engine.domains.MAX_ORFS_SCANNED):
        return []
    monkeypatch.setattr(engine.domains, "scan_domains", fake_domains)
    r = engine.run_analyze({"sequence": ">ok\n" + "ACGT" * 30 + "\n>bad\nACGTZ"})
    ok, bad = r["records"]
    assert ok["processing"]["structural"] == {"status": "complete", "count": 0}
    assert ok["processing"]["domains"]["status"] == "partial"
    assert ok["processing"]["domains"]["orfsNotScanned"] == 1
    assert bad["processing"]["structural"]["status"] == "not_run"
    assert bad["processing"]["classification"]["status"] == "not_run"


def test_analysis_result_document_is_complete_and_rejects_tampering():
    result = engine.run_analyze({"sequence": ">one\n" + "ACGT" * 8 + "\n>two\n" + "TGCA" * 8})
    doc = engine.analysis_result_document(result)
    assert doc["format"] == "TEagle analysis result" and doc["resultSchemaVersion"] == "1.0"
    assert len(doc["records"]) == 2 and doc["status"] == result["status"]
    assert doc["provenance"] == result["provenance"]
    assert engine.validate_analysis_result_document(doc)
    doc["records"].reverse()
    assert not engine.validate_analysis_result_document(doc)
    doc = engine.analysis_result_document(result)
    doc["status"] = "complete_with_stage_failures"
    assert not engine.validate_analysis_result_document(doc)


def test_analysis_result_strict_json_retains_source_url_and_reverifies_seals():
    source_url = "https://www.ncbi.nlm.nih.gov/nuccore/M11240.1"
    sequence = "ACGT" * 30
    result = engine.run_analyze({
        "sequence": ">M11240.1\n" + sequence,
        "source": _fetched_source("M11240.1", sequence),
    })
    encoded = engine.analysis_result_json_bytes(result)
    assert b"NaN" not in encoded and b"Infinity" not in encoded
    document = json.loads(encoded)
    assert document["provenance"]["input"]["sourceUrl"] == source_url
    assert engine.validate_analysis_result_document(document)


def test_analysis_result_strict_json_rejects_nonfinite_mutation_and_preserves_source():
    result = engine.run_analyze({"sequence": ">x\n" + "ACGT" * 30})
    result["records"][0]["composition"]["gc"] = float("nan")
    with pytest.raises(ValueError):
        engine.analysis_result_json_bytes(result)


def test_analysis_result_public_ids_hide_local_headers_without_record_conflation():
    from teagle_core.gff3 import public_record_identity

    raw_ids = [r"C:\Users\person\private\sample.fa", "/home/person/private/sample.fa",
               "../../private/sample.fa", "patient_042", "donor-A1", "sample_tuna_birgun",
               "descriptive specimen", "M11240.1", "chr1", "P1", "record_deadbeefdead",
               "input_sequence"]
    fasta = "".join(f">{record_id}\n{'ACGT' * 30}\n" for record_id in raw_ids)
    result = engine.run_analyze({"sequence": fasta})
    expected = [public_record_identity(record_id) for record_id in raw_ids]
    assert [{k: rec[k] for k in ("id", "idSha256")} for rec in result["records"]] == expected
    assert len({rec["id"] for rec in result["records"]}) == len(raw_ids)
    bound = result["provenance"]["input"]["records"]
    assert [item["idSha256"] for item in bound] == [item["idSha256"] for item in expected]
    encoded = engine.analysis_result_json_bytes(result).decode("utf-8")
    assert all(record_id not in encoded for record_id in raw_ids)
    assert engine.validate_analysis_result_document(json.loads(encoded))


def test_analysis_result_hashes_accession_looking_local_headers():
    from teagle_core.gff3 import public_record_identity
    accessions = ["M11240", "M11240.1", "NC_000013.11", "GCF_009914755.1"]
    fasta = "".join(f">{accession}\n{'ACGT' * 30}\n" for accession in accessions)
    result = engine.run_analyze({"sequence": fasta})
    expected = [public_record_identity(accession) for accession in accessions]
    assert [record["id"] for record in result["records"]] == [item["id"] for item in expected]
    assert [record["idSha256"] for record in result["records"]] == [item["idSha256"] for item in expected]
    assert engine.validate_analysis_result_document(json.loads(engine.analysis_result_json_bytes(result)))


def test_analysis_result_preserves_accession_only_from_validated_fetch_context():
    accession = "M11240.1"
    sequence = "ACGT" * 30
    result = engine.run_analyze({"sequence": f">{accession}\n{sequence}",
                                 "source": _fetched_source(accession, sequence)})
    assert result["records"][0]["id"] == accession
    assert result["records"][0]["identitySource"] == "validated_fetch"
    assert "idSha256" not in result["records"][0]
    assert engine.validate_analysis_result_document(json.loads(engine.analysis_result_json_bytes(result)))


def test_analysis_source_schema_omits_free_text_and_forged_trust_fields():
    accession = "M11240.1"
    sequence = "ACGT" * 30
    raw = {
        "accession": accession, "source": "NCBI nuccore",
        "endpoint": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
        "retrievedUtc": "2026-08-21T00:00:00+00:00",
        "sourceUrl": f"https://www.ncbi.nlm.nih.gov/nuccore/{accession}", "length": len(sequence),
        "organism": "patient_042", "title": "descriptive specimen", "assemblyName": "donor-A1",
        "endpointPath": r"C:\private\endpoint", "identityTrust": "validated_fetch", "unknown": "/private/path",
    }
    source = provenance.normalize_source_context(raw, [(accession, sequence)])
    result = engine.run_analyze({"sequence": f">{accession}\n{sequence}", "source": source})
    public_input = result["provenance"]["input"]
    assert public_input["inputType"] == "fetched_accession"
    for forbidden in ("organism", "title", "assemblyName", "endpointPath", "identityTrust", "unknown"):
        assert forbidden not in public_input
    encoded = engine.analysis_result_json_bytes(result).decode("utf-8")
    assert all(value not in encoded for value in ("patient_042", "descriptive specimen", "donor-A1",
                                                   r"C:\private\endpoint", "/private/path"))


@pytest.mark.parametrize("mutation", [
    {"endpoint": "efetch"},
    {"retrievedUtc": "2026-08-21T00:00:00.123+00:00"},
    {"sourceUrl": "https://www.ncbi.nlm.nih.gov.attacker.example/nuccore/M11240.1"},
    {"identitySource": "validated_fetch"},
])
def test_forged_or_inconsistent_fetch_context_cannot_preserve_accession(mutation):
    accession = "M11240.1"
    sequence = "ACGT" * 30
    source = _fetched_source(accession, sequence)
    source.update(mutation)
    result = engine.run_analyze({"sequence": f">{accession}\n{sequence}", "source": source})
    assert result["records"][0]["id"] == engine.gff3.public_record_identity(accession)["id"]
    assert result["records"][0]["identitySource"] == "local_input"
    assert result["provenance"]["input"]["inputType"] == "local_input"


def test_validated_wire_survives_all_sequence_operation_manifests(monkeypatch):
    accession = "M11240.1"
    sequence = "ACGT" * 60
    fasta = f">{accession}\n{sequence}"
    source = _fetched_source(accession, sequence)

    monkeypatch.setattr(engine.wsl, "annotate", lambda *a, **k: {
        "ok": True, "hits": [], "dfam_library": {"version": "fixture"},
        "repeatmasker_version": "fixture", "library_kind": "curated", "species": None})
    monkeypatch.setattr(engine.wsl, "splice_align", lambda *a, **k: {
        "ok": True, "exons": [], "introns": [], "minimap2_version": "fixture"})
    monkeypatch.setattr(engine.wsl, "protein_align", lambda *a, **k: {
        "ok": True, "features": [], "miniprot_version": "fixture"})
    monkeypatch.setattr(engine.primers, "design_primers", lambda *a, **k: {
        "candidates": [], "explain_left": "", "explain_right": "", "explain_pair": ""})
    monkeypatch.setattr(engine.primers, "in_silico_pcr", lambda *a, **k: [])

    results = [
        engine.run_annotate({"sequence": fasta, "source": source}),
        engine.run_splice({"sequence": fasta, "transcript": "ACGT" * 20, "source": source}),
        engine.run_miniprot({"sequence": fasta, "protein": ">p\nMPEPTIDE", "source": source}),
        engine.run_primers({"sequence": fasta, "source": source}),
        engine.run_pcr({"sequence": fasta, "source": source,
                        "fwd": "ACGTACGTACGTACGTAC", "rev": "TGCATGCATGCATGCATG"}),
    ]
    for result in results:
        manifest_input = result["provenance"]["input"]
        assert manifest_input["id"] == accession
        assert manifest_input["identitySource"] == "validated_fetch"
        assert manifest_input["inputType"] == "fetched_accession"
        assert manifest_input["source"] == "ncbi_nuccore"
        assert manifest_input["sourceUrl"] == f"https://www.ncbi.nlm.nih.gov/nuccore/{accession}"
        assert provenance.verify_manifest(result["provenance"])


def test_forged_wire_downgrades_all_sequence_operation_manifests(monkeypatch):
    accession = "M11240.1"
    sequence = "ACGT" * 60
    source = _fetched_source(accession, sequence)
    source["records"][0]["sha256"] = "0" * 64
    monkeypatch.setattr(engine.wsl, "splice_align", lambda *a, **k: {
        "ok": True, "exons": [], "introns": [], "minimap2_version": "fixture"})
    result = engine.run_splice({"sequence": f">{accession}\n{sequence}",
                                "transcript": "ACGT" * 20, "source": source})
    manifest_input = result["provenance"]["input"]
    assert manifest_input["inputType"] == "local_input"
    assert manifest_input["identitySource"] == "local_input"
    assert "accession" not in manifest_input and "sourceUrl" not in manifest_input


def test_pcr_sources_and_bindings_hide_local_fasta_headers(monkeypatch):
    from teagle_core.gff3 import public_record_identity

    def fake_pcr(_fwd, _rev, _seq, source, stats=None, **_kwargs):
        if stats is not None:
            stats.update(amplicons_capped=False, sites_capped=False)
        return [{"source": source, "start": 0, "end": 80, "length": 80, "product_length": 80,
                 "product_key": source + ":0:80", "mechanism_key": "F+R", "fwd_mm": 0, "rev_mm": 0}]

    monkeypatch.setattr(engine.primers, "in_silico_pcr", fake_pcr)
    assay_id = r"C:\Users\person\private\assay.fa"
    background_id = "/home/person/private/background.fa"
    result = engine.run_pcr({"sequence": f">{assay_id}\n" + "A" * 200,
                             "background": f">{background_id}\n" + "A" * 200,
                             "fwd": "A" * 18, "rev": "T" * 18})
    assay = public_record_identity(assay_id)
    background = public_record_identity(background_id)
    sources = {amplicon["source"] for amplicon in result["amplicons"]}
    assert sources == {assay["id"] + " (input)", background["id"] + " (custom background)"}
    assert {amplicon["sourceIdSha256"] for amplicon in result["amplicons"]} == {
        assay["idSha256"], background["idSha256"]}
    assert all("idsha256=" in amplicon["product_key"] for amplicon in result["amplicons"])
    assert result["bindings"]["assay"]["idSha256"] == assay["idSha256"]
    assert result["bindings"]["backgrounds"][0]["idSha256"] == background["idSha256"]
    encoded = json.dumps(result, ensure_ascii=False)
    assert assay_id not in encoded and background_id not in encoded


@pytest.mark.parametrize("family,mutate", [
    ("composition", lambda d: d["records"][0]["composition"].__setitem__("gc", 99.0)),
    ("structural", lambda d: d["records"][0]["structural"].append({"type": "fabricated"})),
    ("domains", lambda d: d["records"][0]["domains"].append({"domain": "fabricated"})),
    ("orfs", lambda d: d["records"][0]["orfs"].append({"start": 0, "end": 3})),
    ("classification", lambda d: d["records"][0]["classification"].__setitem__("superfamily", "fabricated")),
    ("notes", lambda d: d["records"][0]["notes"].append("fabricated")),
    ("processing", lambda d: d["records"][0]["processing"]["structural"].__setitem__("status", "fabricated")),
    ("invalid", lambda d: d["records"][0]["invalid"].append({"pos": 0, "char": "Z"})),
    ("references", lambda d: d["references"].append({"name": "fabricated"})),
])
def test_analysis_result_payload_seal_rejects_scientific_mutation(family, mutate):
    result = engine.run_analyze({"sequence": ">x\n" + "ACGT" * 30})
    doc = engine.analysis_result_document(result)
    mutate(doc)
    assert not engine.validate_analysis_result_document(doc), family


def _json_roundtrip(value):
    return json.loads(json.dumps(value, ensure_ascii=False, allow_nan=False))


def _reseal_all_unkeyed_hashes(document):
    document["resultPayloadSha256"] = engine.scientific_result_sha256(document)
    document["provenance"]["resultDescriptor"] = {
        "resultSchemaVersion": document["resultSchemaVersion"],
        "resultPayloadSha256": document["resultPayloadSha256"],
    }
    document["provenance"]["manifestSha256"] = hashlib.sha256(
        provenance.canonical_manifest_bytes(document["provenance"])).hexdigest()


@pytest.mark.parametrize("family,mutate", [
    ("composition", lambda d: d["records"][0]["composition"].__setitem__("gc", 99.0)),
    ("structural", lambda d: d["records"][0]["structural"].append({"type": "fabricated"})),
    ("domains", lambda d: d["records"][0]["domains"].append({"domain": "fabricated"})),
    ("orfs", lambda d: d["records"][0]["orfs"].append({"start": 0, "end": 3})),
    ("classification", lambda d: d["records"][0]["classification"].__setitem__("superfamily", "fabricated")),
    ("notes", lambda d: d["records"][0]["notes"].append("fabricated")),
    ("processing", lambda d: d["records"][0]["processing"]["structural"].__setitem__("status", "fabricated")),
    ("invalid", lambda d: d["records"][0]["invalid"].append({"pos": 0, "char": "Z"})),
    ("references", lambda d: d["references"].append({"name": "fabricated"})),
])
def test_manifest_descriptor_rejects_json_roundtrip_mutation_with_recomputed_payload(family, mutate):
    doc = _json_roundtrip(engine.analysis_result_document(
        engine.run_analyze({"sequence": ">x\n" + "ACGT" * 30})))
    sealed_descriptor = copy.deepcopy(doc["provenance"]["resultDescriptor"])
    mutate(doc)
    doc["resultPayloadSha256"] = engine.scientific_result_sha256(doc)
    assert engine.verify_scientific_result_payload(doc)
    assert provenance.verify_manifest(doc["provenance"])
    assert doc["provenance"]["resultDescriptor"] == sealed_descriptor
    assert doc["provenance"]["resultDescriptor"]["resultPayloadSha256"] != doc["resultPayloadSha256"]
    assert not engine.validate_analysis_result_document(doc), family


@pytest.mark.parametrize("case,mutate", [
    ("unknown_status", lambda p: p["structural"].__setitem__("status", "fabricated")),
    ("count_mismatch", lambda p: p["structural"].__setitem__("count", 999)),
    ("stage_insertion", lambda p: p.__setitem__("fabricated", {"status": "not_run"})),
    ("stage_omission", lambda p: p.pop("retroviral")),
    ("stage_duplicate", lambda p: p.__setitem__("structural_duplicate", copy.deepcopy(p["structural"]))),
    ("unknown_stage_field", lambda p: p["structural"].__setitem__("fabricated", True)),
    ("classification_outcome_mismatch", lambda p: p["classification"].__setitem__("status", "withheld")),
    ("partial_without_cap", lambda p: p["domains"].__setitem__("status", "partial")),
])
def test_processing_contract_rejects_mutation_even_after_all_unkeyed_hashes_are_recomputed(case, mutate):
    doc = _json_roundtrip(engine.analysis_result_document(
        engine.run_analyze({"sequence": ">x\n" + "ACGT" * 30})))
    mutate(doc["records"][0]["processing"])
    _reseal_all_unkeyed_hashes(doc)
    assert engine.verify_scientific_result_payload(doc)
    assert provenance.verify_manifest(doc["provenance"])
    assert not engine.validate_analysis_result_document(doc), case


def test_analysis_result_descriptor_and_manifest_seal_mutations_fail_closed():
    result = engine.run_analyze({"sequence": ">x\n" + "ACGT" * 30})
    descriptor_doc = _json_roundtrip(engine.analysis_result_document(result))
    descriptor_doc["provenance"]["resultDescriptor"]["resultPayloadSha256"] = "0" * 64
    descriptor_doc["provenance"]["manifestSha256"] = hashlib.sha256(
        provenance.canonical_manifest_bytes(descriptor_doc["provenance"])).hexdigest()
    assert provenance.verify_manifest(descriptor_doc["provenance"])
    assert not engine.validate_analysis_result_document(descriptor_doc)

    seal_doc = _json_roundtrip(engine.analysis_result_document(result))
    seal_doc["provenance"]["parameters"]["orf_min_aa"] += 1
    assert not provenance.verify_manifest(seal_doc["provenance"])
    assert not engine.validate_analysis_result_document(seal_doc)


def test_analysis_result_digest_and_provenance_seal_have_no_circular_dependency():
    result = engine.run_analyze({"sequence": ">x\n" + "ACGT" * 30})
    digest = result["resultPayloadSha256"]
    without_provenance = copy.deepcopy(result)
    without_provenance.pop("provenance")
    assert engine.scientific_result_sha256(without_provenance) == digest
    assert engine.scientific_result_sha256(result) == digest
    assert result["provenance"]["resultDescriptor"] == {
        "resultSchemaVersion": engine.ANALYSIS_RESULT_SCHEMA,
        "resultPayloadSha256": digest,
    }
    assert provenance.verify_manifest(result["provenance"])


def test_analysis_statuses_are_centralized_and_nonexportable_fail_closed():
    empty = engine.run_analyze({"sequence": ""})
    invalid = engine.run_analyze({"sequence": ">bad\nACGTZ"})
    mixed = engine.run_analyze({"sequence": ">ok\n" + "ACGT" * 8 + "\n>bad\nACGTZ"})
    assert empty["status"] == engine.ANALYSIS_STATUS_NO_INPUT
    assert invalid["status"] == engine.ANALYSIS_STATUS_NO_VALID_RECORDS
    assert mixed["status"] == engine.ANALYSIS_STATUS_COMPLETE_WITH_INVALID_RECORDS
    assert mixed["status"] == engine.derive_analysis_status(mixed["records"])
    with pytest.raises(ValueError):
        engine.analysis_result_document(empty)
    with pytest.raises(ValueError):
        engine.analysis_result_document(invalid)
    assert engine.validate_analysis_result_document(engine.analysis_result_document(mixed))


def test_analysis_result_validator_rejects_unknown_schema_status_and_nonfinite():
    result = engine.run_analyze({"sequence": ">x\n" + "ACGT" * 30})
    for key, value in (("resultSchemaVersion", "9.9"), ("status", "unknown")):
        doc = engine.analysis_result_document(result)
        doc[key] = value
        doc["resultPayloadSha256"] = engine.scientific_result_sha256(doc)
        assert not engine.validate_analysis_result_document(doc)
    doc = engine.analysis_result_document(result)
    doc["records"][0]["composition"]["gc"] = float("nan")
    assert not engine.validate_analysis_result_document(doc)


def test_analysis_status_derivation_cannot_drift_from_constructor_or_validator():
    result = engine.run_analyze({"sequence": ">ok\n" + "ACGT" * 8 + "\n>bad\nACGTZ"})
    assert result["status"] == engine.derive_analysis_status(result["records"])
    doc = engine.analysis_result_document(result)
    doc["status"] = engine.ANALYSIS_STATUS_COMPLETE
    doc["resultPayloadSha256"] = engine.scientific_result_sha256(doc)
    assert not engine.validate_analysis_result_document(doc)
