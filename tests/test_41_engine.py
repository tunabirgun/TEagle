"""4.1 audit regressions for engine state and provenance (lane A4)."""
import copy
import pytest

import engine
from engine import BadRequest

F, R = "A" * 20, "C" * 20


def _pair(src, lo, hi, name="pair"):
    return f">{src}:{lo}+{hi} {name} {hi - lo + 1}bp {F} {R}"


def _scan(monkeypatch, *rows):
    raw = "\n".join(rows) + "\n"
    monkeypatch.setattr(engine.wsl, "genome_scan",
                        lambda acc, query, **k: {"ok": True, "raw": raw, "isPcr_version": "33x2",
                                                 "target": "genome.2bit", "sha256": "s", "n_seqs": 2})


def _genome(locus, **extra):
    return engine.run_genome_pcr({"fwd": F, "rev": R, "assemblyAccession": "GCF_9.9", "taxid": "1",
                                  "design_locus": locus, **extra})


def _on(res):
    return [(a["start"], a["end"]) for a in res["amplicons"] if a.get("on_target")]


# ---------- genomepcr_gel-1 ----------
def test_genome_on_target_is_the_designed_product_not_the_longest(monkeypatch):
    _scan(monkeypatch, _pair("NC_1", 5001, 5200), _pair("NC_1", 12001, 12450))
    r = _genome({"accession": "NC_1", "start": 1, "end": 20000, "strand": 1,
                 "intended_product_span": [5000, 5200]})
    assert _on(r) == [(5001, 5200)] and r["n_on_target"] == 1 and r["has_locus"] is True
    assert r["summary"]["n_on"] == 1 and r["summary"]["n_off"] == 1


def test_genome_on_target_beats_a_spanning_product(monkeypatch):
    _scan(monkeypatch, _pair("NC_1", 1101, 1400), _pair("NC_1", 1101, 4900), _pair("NC_1", 5101, 5400))
    r = _genome({"accession": "NC_1", "start": 1001, "end": 6000, "strand": 1,
                 "intended_product_span": [100, 400]})
    assert _on(r) == [(1101, 1400)]


def test_genome_tie_resolves_to_the_designed_copy_not_the_leftmost(monkeypatch):
    _scan(monkeypatch, _pair("NC_1", 1101, 1400), _pair("NC_1", 5101, 5400))
    r = _genome({"accession": "NC_1", "start": 1001, "end": 6000, "strand": 1,
                 "intended_product_span": [4100, 4400]})
    assert _on(r) == [(5101, 5400)]


@pytest.mark.parametrize("strand", [2, "-"])
def test_genome_span_maps_through_a_minus_strand_fetch(monkeypatch, strand):
    # record position [100,300) of the reverse-complemented fetch NC_1:1000-2999 is genome 2700-2899
    _scan(monkeypatch, _pair("NC_1", 1100, 1299), _pair("NC_1", 2700, 2899))
    r = _genome({"accession": "NC_1", "start": 1000, "end": 2999, "strand": strand,
                 "intended_product_span": [100, 300]})
    assert _on(r) == [(2700, 2899)]


def test_genome_designed_product_absent_promotes_nothing(monkeypatch):
    _scan(monkeypatch, _pair("NC_1", 12001, 12450))
    r = _genome({"accession": "NC_1", "start": 1, "end": 20000, "strand": 1,
                 "intended_product_span": [5000, 5200]})
    assert _on(r) == [] and r["n_on_target"] == 0 and r["has_locus"] is True
    assert r["summary"]["tier"] == "intended-locus-negative"


def test_genome_window_without_span_is_neutral_not_overlap(monkeypatch):
    _scan(monkeypatch, _pair("NC_1", 1000, 1240))
    for locus in ({"accession": "NC_1", "start": 950, "stop": 1300},
                  {"accession": "NC_1", "start": 950, "end": 1300, "strand": 1}):
        r = _genome(locus)
        assert r["has_locus"] is False and r["n_on_target"] == 0 and _on(r) == []


def test_genome_missing_window_end_is_not_coerced_to_zero(monkeypatch):
    _scan(monkeypatch, _pair("NC_1", 500, 700))
    r = _genome({"accession": "NC_1", "start": 1001, "intended_product_span": [0, 100]})
    assert r["has_locus"] is False and _on(r) == []


def test_genome_endpoint_tolerance_is_honoured(monkeypatch):
    _scan(monkeypatch, _pair("NC_1", 5001, 5203))
    base = {"accession": "NC_1", "start": 1, "end": 20000, "strand": 1, "intended_product_span": [5000, 5200]}
    assert _on(_genome(base)) == []
    assert _on(_genome({**base, "intended_endpoint_tolerance_bp": 5})) == [(5001, 5203)]


def test_genome_single_primer_products_never_on_target(monkeypatch):
    _scan(monkeypatch, _pair("NC_1", 5001, 5200, "fwdonly"))
    r = _genome({"accession": "NC_1", "start": 1, "end": 20000, "strand": 1,
                 "intended_product_span": [5000, 5200]})
    assert _on(r) == []


def test_genome_span_outside_window_is_a_bad_request(monkeypatch):
    _scan(monkeypatch, _pair("NC_1", 5001, 5200))
    with pytest.raises(BadRequest):
        _genome({"accession": "NC_1", "start": 1, "end": 1000, "strand": 1, "intended_product_span": [900, 1500]})


def test_genome_span_labelling_does_not_change_the_seal(monkeypatch):
    _scan(monkeypatch, _pair("NC_1", 5001, 5200), _pair("NC_1", 12001, 12450))
    base = {"fwd": F, "rev": R, "assemblyAccession": "GCF_9.9", "taxid": "1"}
    plain = engine.run_genome_pcr(dict(base))
    spanned = _genome({"accession": "NC_1", "start": 1, "end": 20000, "strand": 1,
                       "intended_product_span": [5000, 5200]})
    assert plain["provenance"]["manifestSha256"] == spanned["provenance"]["manifestSha256"]


# ---------- provenance-1 ----------
_NCBI_RANGE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def _coord_fetch(regions):
    """(fasta text, normalized source context) for fetched coordinate regions [(start, stop, strand)]."""
    raw = {"accession": "NC_000013.11", "taxid": "9606", "assemblyAccession": "GCF_000001405.40",
           "regions": [{"chrAccession": "NC_000013.11", "start": s, "stop": e, "strand": st} for s, e, st in regions],
           "coordSystem": "1-based-inclusive", "retrievalType": "coordinate",
           "source": "NCBI E-utilities (efetch seq_start/seq_stop)", "endpoint": _NCBI_RANGE,
           "retrievedUtc": "2026-07-21T10:00:00+00:00"}
    recs = [(f"NC_000013.11:{s}-{e}" + ("(-)" if st == 2 else ""), "ACGT" * ((e - s + 1) // 4) + "A" * ((e - s + 1) % 4))
            for s, e, st in regions]
    from teagle_core import provenance
    return "".join(f">{h}\n{q}\n" for h, q in recs), provenance.normalize_source_context(raw, recs)


def test_multi_region_coordinate_analysis_verifies_and_exports():
    from teagle_core import provenance
    text, src = _coord_fetch([(1, 40, 1), (101, 160, 2)])
    assert src["inputType"] == "fetched_coordinate" and len(src["regions"]) == 2
    res = engine.analyze(text, src)
    m = res["provenance"]
    assert m["input"]["id"] == "multiple_records" and m["input"]["identitySource"] == "generated"
    assert [r["id"] for r in res["records"]] == ["NC_000013.11:1-40", "NC_000013.11:101-160"]
    assert provenance.verify_manifest(m)
    assert engine.analysis_result_document(res)["status"] == "complete"
    assert engine.analysis_result_json_bytes(res)


def _resealed(manifest, mutate):
    from teagle_core import provenance
    m = copy.deepcopy(manifest)
    mutate(m)
    import hashlib
    m["manifestSha256"] = hashlib.sha256(provenance.canonical_manifest_bytes(m)).hexdigest()
    return m


def test_multi_region_identity_rules_still_reject_contradictions():
    from teagle_core import provenance
    text, src = _coord_fetch([(1, 40, 1), (101, 160, 1)])
    multi = engine.analyze(text, src)["provenance"]
    assert not provenance.verify_manifest(_resealed(multi, lambda m: m["input"].update(
        id="NC_000013.11:1-40", identitySource="validated_fetch")))
    assert not provenance.verify_manifest(_resealed(multi, lambda m: m["input"].update(id="record")))
    assert not provenance.verify_manifest(_resealed(multi, lambda m: m["input"]["records"][0].update(
        identitySource="generated")))
    one_text, one_src = _coord_fetch([(1, 40, 1)])
    single = engine.analyze(one_text, one_src)["provenance"]
    assert provenance.verify_manifest(single)
    assert not provenance.verify_manifest(_resealed(single, lambda m: m["input"].update(
        id="multiple_records", identitySource="generated")))


# ---------- provenance-2 ----------
_ENA_API = "https://www.ebi.ac.uk/ena/browser/api/fasta/"
_SEQ = "ACGT" * 60


def _accession_raw(typed, resolved, ena=False, **over):
    """The context fetch.retrieve builds: the resolved accession-version beside URLs from the typed text."""
    raw = {"accession": resolved, "taxid": "7227", "retrievedUtc": "2026-08-21T00:00:00+00:00", "length": len(_SEQ)}
    if ena:
        raw.update(source="ENA (EMBL-EBI)", endpoint=_ENA_API + typed,
                   sourceUrl="https://www.ebi.ac.uk/ena/browser/view/" + typed)
    else:
        raw.update(source="NCBI nuccore", endpoint=_NCBI_RANGE, sourceUrl="https://www.ncbi.nlm.nih.gov/nuccore/" + typed)
    raw.update(over)
    return raw


def _normalize(raw, header):
    from teagle_core import provenance
    return provenance.normalize_source_context(raw, [(header, _SEQ)])


def test_unversioned_typed_accession_keeps_fetched_provenance():
    ctx = _normalize(_accession_raw("M11240", "M11240.1"), "M11240.1 Drosophila melanogaster")
    assert ctx["inputType"] == "fetched_accession" and ctx["accession"] == "M11240.1"
    assert ctx["sourceUrl"] == "https://www.ncbi.nlm.nih.gov/nuccore/M11240.1"
    res = engine.analyze(f">M11240.1 Drosophila\n{_SEQ}\n", ctx)
    assert res["records"][0]["id"] == "M11240.1" and res["provenance"]["input"]["accession"] == "M11240.1"


def test_ena_header_form_and_unversioned_endpoint_keep_fetched_provenance():
    from teagle_core import provenance
    ctx = _normalize(_accession_raw("M11240", "M11240.1", ena=True), "ENA|M11240|M11240.1 Drosophila melanogaster")
    assert ctx["inputType"] == "fetched_accession" and ctx["source"] == "ena_fasta"
    res = engine.analyze(f">ENA|M11240|M11240.1 Drosophila\n{_SEQ}\n", ctx)
    assert res["records"][0]["id"] == "M11240.1" and provenance.verify_manifest(res["provenance"])
    ncbi = _normalize(_accession_raw("M11240", "M11240.1"), "M11240.1 Drosophila melanogaster")
    ncbi_res = engine.analyze(f">M11240.1 Drosophila\n{_SEQ}\n", ncbi)
    assert res["provenance"]["manifestSha256"] == ncbi_res["provenance"]["manifestSha256"]


@pytest.mark.parametrize("raw_over,header", [
    ({"sourceUrl": "https://www.ncbi.nlm.nih.gov/nuccore/M11241"}, "M11240.1 x"),            # URL names another record
    ({"sourceUrl": "https://example.org/nuccore/M11240"}, "M11240.1 x"),                    # foreign host
    ({"sourceUrl": "https://www.ncbi.nlm.nih.gov/nuccore/M11240.2"}, "M11240.1 x"),         # other version
    ({"sourceUrl": "https://www.ncbi.nlm.nih.gov/nuccore/M11240?x=1"}, "M11240.1 x"),
    ({}, "M11241.1 x"),                                                                      # header names another record
    ({}, "ENA|M11240|M11241.1 x"),
])
def test_forged_or_mismatched_accession_context_still_downgrades(raw_over, header):
    assert _normalize(_accession_raw("M11240", "M11240.1", **raw_over), header)["inputType"] == "local_input"


def test_typed_version_that_differs_from_resolved_still_downgrades():
    raw = _accession_raw("M11240.2", "M11240.1")
    assert _normalize(raw, "M11240.1 x")["inputType"] == "local_input"


# ---------- engine_state-1 ----------
def _trunc_notes(seq):
    rec = engine.analyze(">r\n" + seq + "\n")["records"][0]
    return rec, [n for n in rec["notes"] if "runs off" in n or "truncat" in n.lower()]


@pytest.mark.parametrize("acc", ["M80343", "X01005", "X05424", "M11240", "M12927"])
def test_intact_reference_elements_carry_no_truncation_note(acc):
    from helpers import fixture_seq
    rec, notes = _trunc_notes(fixture_seq(acc))
    assert notes == [], (acc, notes)


def test_truncation_note_names_the_three_prime_end_for_a_cut_forward_element():
    from helpers import fixture_seq
    seq = fixture_seq("M11240")[:4000]
    rec, notes = _trunc_notes(seq)
    assert len(notes) == 1 and "3′" in notes[0] and "5′" not in notes[0]
    assert rec["classification"]["completeness"]["coding"]["truncated_frame"] is True


def test_truncation_note_names_the_five_prime_end_for_the_reverse_complement():
    from helpers import fixture_seq
    from teagle_core import sequtil
    seq = sequtil.reverse_complement(fixture_seq("M11240")[:4000])
    rec, notes = _trunc_notes(seq)
    assert len(notes) == 1 and "5′" in notes[0] and "3′" not in notes[0]


# ---------- structural-9 ----------
def _sealed(res):
    return res["provenance"]["parameters"]


def test_seal_carries_the_extension_and_poly_a_signal_thresholds_and_tsd_table():
    p = _sealed(engine.analyze(">r\n" + "ACGT" * 50 + "\n"))
    for key in ("ltr_extend_x_drop", "polya_signal_dse_lo", "polya_signal_dse_hi",
                "polya_signal_dse_min_gt_frac", "polya_signal_dse_min_t", "tsd_expectation_by_superfamily"):
        assert key in p, key
    assert p["tsd_expectation_by_superfamily"]["Tc1/Mariner"]["length"] == 2


def test_seal_follows_the_applied_values_not_a_restated_copy(monkeypatch):
    from teagle_core import structural
    monkeypatch.setattr(structural.find_polya_signal, "__defaults__", (25, 70, 0.7, 5))
    monkeypatch.setattr(structural._extend, "__defaults__", (45,))
    monkeypatch.setitem(structural.TSD_EXPECT, "Tc1/Mariner", (3, "test"))
    p = _sealed(engine.analyze(">r\n" + "ACGT" * 50 + "\n"))
    assert (p["polya_signal_dse_lo"], p["polya_signal_dse_hi"], p["polya_signal_dse_min_gt_frac"],
            p["polya_signal_dse_min_t"], p["ltr_extend_x_drop"]) == (25, 70, 0.7, 5, 45)
    assert p["tsd_expectation_by_superfamily"]["Tc1/Mariner"] == {"length": 3}


@pytest.mark.parametrize("key,attr", sorted(engine._SEALED_STRUCTURAL_CONSTANTS.items()))
def test_sealed_structural_constants_follow_their_module(monkeypatch, key, attr):
    from teagle_core import structural
    if not hasattr(structural, attr):
        pytest.skip(f"structural.{attr} is not defined")
    monkeypatch.setattr(structural, attr, getattr(structural, attr) + 1)
    assert engine._detector_parameters()[key] == getattr(structural, attr)


def test_a_changed_seal_changes_the_manifest_hash(monkeypatch):
    from teagle_core import structural
    text = ">r\n" + "ACGT" * 50 + "\n"
    before = engine.analyze(text)["provenance"]["manifestSha256"]
    monkeypatch.setattr(structural._extend, "__defaults__", (31,))
    assert engine.analyze(text)["provenance"]["manifestSha256"] != before


# ---------- fetch_refs-11 ----------
def _ref_keys(res):
    return {r["key"] for r in res["references"]}


def test_domain_negative_analysis_still_cites_hmmer_and_pfam():
    from helpers import make_ltr_element
    res = engine.analyze(">r\n" + make_ltr_element(seed=1) + "\n")
    rec = res["records"][0]
    assert rec["domains"] == [] and rec["processing"]["domains"]["orfsScanned"] > 0
    assert {"HMMER", "Pfam"} <= _ref_keys(res)


def test_analysis_that_searched_no_orf_does_not_cite_hmmer():
    res = engine.analyze(">r\n" + "ACGT" * 5 + "\n")
    assert res["records"][0]["processing"]["domains"]["orfsScanned"] == 0
    assert not {"HMMER", "Pfam"} & _ref_keys(res)


# ---------- engine_state-6 ----------
def test_empty_record_is_invalid_with_a_reason_and_never_classified():
    import json
    res = engine.analyze(">a\n>b\n" + "ACGT" * 30 + "\n")
    a, b = res["records"]
    assert a["valid"] is False and a["classification"] is None and a["invalidReason"]
    assert all(s == {"status": "not_run"} for s in a["processing"].values())
    assert b["valid"] is True and b["classification"] is not None and "invalidReason" not in b
    assert res["status"] == engine.ANALYSIS_STATUS_COMPLETE_WITH_INVALID_RECORDS
    assert "empty" in res["warning"].lower()
    doc = engine.analysis_result_document(res)
    assert engine.validate_analysis_result_document(json.loads(engine.analysis_result_json_bytes(res)))
    assert doc["records"][0]["invalidReason"] == a["invalidReason"]


@pytest.mark.parametrize("seq", ["N" * 40, "NNRYNNKMNNSWNN"])
def test_record_without_informative_bases_is_not_a_complete_negative(seq):
    res = engine.analyze(f">n\n{seq}\n")
    rec = res["records"][0]
    assert rec["valid"] is False and rec["classification"] is None and "no informative" in rec["invalidReason"]
    assert res["status"] == engine.ANALYSIS_STATUS_NO_VALID_RECORDS
    assert "no informative" in res["warning"]


def test_ordinary_record_is_unchanged_by_the_empty_record_rule():
    from helpers import make_ltr_element
    rec = engine.analyze(">r\n" + make_ltr_element(seed=2) + "\n")["records"][0]
    assert rec["valid"] is True and "invalidReason" not in rec and rec["classification"] is not None


def test_invalid_reason_is_only_accepted_on_an_invalid_record():
    res = engine.analyze(">a\n>b\n" + "ACGT" * 30 + "\n")
    ok = copy.deepcopy(res)
    ok["records"][1]["invalidReason"] = "forged"
    assert engine._analysis_processing_valid(res["records"][0]) and engine._analysis_processing_valid(res["records"][1])
    assert not engine._analysis_processing_valid(ok["records"][1])
    bare = copy.deepcopy(res["records"][0])
    del bare["invalidReason"]
    assert not engine._analysis_processing_valid(bare)


# ---------- engine_state-14 ----------
def test_stage_summary_counts_stage_states_of_valid_records_only(monkeypatch):
    res = engine.analyze(">bad\nACGTX\n>ok\n" + "ACGT" * 30 + "\n")
    assert res["stageSummary"]["structural"] == {"complete": 1}
    monkeypatch.setattr(engine.structural, "detect_all", lambda seq: 1 / 0)
    failed = engine.analyze(">ok\n" + "ACGT" * 30 + "\n")
    assert failed["stageSummary"]["structural"] == {"failed": 1}


def test_partial_domain_stages_are_visible_in_the_stage_summary():
    from helpers import fixture_seq
    partial = 0
    for acc in ("M11240", "M12927", "M80343", "X01005", "X05424"):
        res = engine.analyze(f">{acc}\n{fixture_seq(acc)}\n")
        partial += res["stageSummary"]["domains"].get("partial", 0)
    assert partial > 0


def test_tampered_stage_summary_is_rejected_and_old_documents_without_it_still_verify():
    res = engine.analyze(">ok\n" + "ACGT" * 30 + "\n")
    doc = engine.analysis_result_document(res)
    forged = copy.deepcopy(doc)
    forged["stageSummary"]["domains"] = {"partial": 1}
    forged["resultPayloadSha256"] = engine.scientific_result_sha256(forged)
    forged["provenance"]["resultDescriptor"]["resultPayloadSha256"] = forged["resultPayloadSha256"]
    assert not engine.validate_analysis_result_document(forged)
    legacy = copy.deepcopy(doc)
    del legacy["stageSummary"]
    legacy["resultPayloadSha256"] = engine.scientific_result_sha256(legacy)
    legacy["provenance"] = _resealed(legacy["provenance"], lambda m: m["resultDescriptor"].update(
        resultPayloadSha256=legacy["resultPayloadSha256"]))
    assert engine.validate_analysis_result_document(legacy)


def test_multi_record_warning_names_the_first_valid_record_as_default():
    res = engine.analyze(">bad\nACGTX\n>ok\n" + "ACGT" * 30 + "\n>ok2\n" + "ACGT" * 30 + "\n")
    assert "first valid record by default" in res["warning"]


# ---------- persona-U7 ----------
def test_domain_stage_that_searched_no_orf_is_not_run_not_complete():
    from teagle_core import classify
    seq = "ACGTNNACGT" * 9 + "ACGTAC"
    res = engine.analyze(f">short\n{seq}\n")
    rec = res["records"][0]
    assert rec["orfs"] == []
    assert rec["processing"]["domains"] == {"status": "not_run", "count": 0, "orfsScanned": 0, "orfsNotScanned": 0}
    assert res["stageSummary"]["domains"] == {"not_run": 1}
    assert rec["classification"] == classify.classify(
        rec["structural"], rec["domains"], seq=seq, domains_ok=True, orfs_unscanned=0, orfs=[], structural_ok=True)
    assert engine.analysis_result_document(res)


def test_domain_stage_that_read_an_orf_stays_complete_even_without_hits():
    from helpers import make_ltr_element
    rec = engine.analyze(">r\n" + make_ltr_element(seed=1) + "\n")["records"][0]
    assert rec["domains"] == [] and rec["processing"]["domains"]["status"] == "complete"


def test_not_run_domain_state_is_only_valid_when_nothing_was_scanned():
    rec = engine.analyze(">short\n" + "ACGT" * 20 + "\n")["records"][0]
    assert engine._analysis_processing_valid(rec)
    bad = copy.deepcopy(rec)
    bad["processing"]["domains"]["orfsScanned"] = 1
    assert not engine._analysis_processing_valid(bad)
    bad = copy.deepcopy(rec)
    bad["processing"]["domains"] = {"status": "not_run"}
    assert not engine._analysis_processing_valid(bad)


# ---------- engine_state-13 ----------
def _rna_note(rec):
    return any("RNA input" in n for n in rec["notes"])


def test_rna_note_is_per_record():
    dna, rna = "ACGT" * 30, "ACGU" * 30
    res = engine.analyze(f">d\n{dna}\n>r\n{rna}\n>d2\n{dna}\n")
    assert [_rna_note(r) for r in res["records"]] == [False, True, False]


def test_rna_in_text_before_the_first_header_flags_no_record():
    res = engine.analyze("ACGU stray line\n>d\n" + "ACGT" * 30 + "\n")
    assert [_rna_note(r) for r in res["records"]] == [False]


def test_rna_in_a_multiline_record_and_bare_input_still_flag():
    assert _rna_note(engine.analyze("ACGT" * 10 + "\n" + "ACGU" * 10 + "\n")["records"][0])
    assert _rna_note(engine.analyze(">r\n" + "ACGT" * 10 + "\nACGU\n")["records"][0])


# ---------- engine_state-9 ----------
def _stub_operations(monkeypatch):
    monkeypatch.setattr(engine.wsl, "annotate", lambda *a, **k: {
        "ok": True, "hits": [], "dfam_library": {"version": "fixture"},
        "repeatmasker_version": "fixture", "library_kind": "curated", "species": None})
    monkeypatch.setattr(engine.wsl, "splice_align", lambda *a, **k: {
        "ok": True, "exons": [], "introns": [], "minimap2_version": "fixture"})
    monkeypatch.setattr(engine.wsl, "protein_align", lambda *a, **k: {
        "ok": True, "features": [], "miniprot_version": "fixture"})
    monkeypatch.setattr(engine.primers, "in_silico_pcr", lambda *a, **k: [])


def _operations(text):
    sixty = "ACGTACGTACGTACGTAC"
    return {
        "annotate": engine.run_annotate({"sequence": text}),
        "splice": engine.run_splice({"sequence": text, "transcript": "ACGT" * 20}),
        "miniprot": engine.run_miniprot({"sequence": text, "protein": ">p\nMPEPTIDE"}),
        "primers": engine.run_primers({"sequence": text}),
        "pcr": engine.run_pcr({"sequence": text, "fwd": sixty, "rev": "TGCATGCATGCATGCATG"}),
    }


def test_bare_input_operations_carry_the_analysis_record_identity(monkeypatch):
    from teagle_core import provenance
    _stub_operations(monkeypatch)
    from helpers import make_ltr_element
    seq = make_ltr_element(seed=3)
    analysed = engine.analyze(seq)["records"][0]
    for name, result in _operations(seq).items():
        inp = result["provenance"]["input"]
        assert (inp["id"], inp["identitySource"]) == (analysed["id"], analysed["identitySource"]), name
        assert "idSha256" not in inp and provenance.verify_manifest(result["provenance"]), name


def test_headed_input_operations_keep_the_analysis_record_identity(monkeypatch):
    _stub_operations(monkeypatch)
    from helpers import make_ltr_element
    text = ">Copia-like consensus\n" + make_ltr_element(seed=3) + "\n"
    analysed = engine.analyze(text)["records"][0]
    assert analysed["identitySource"] == "local_input"
    for name, result in _operations(text).items():
        inp = result["provenance"]["input"]
        assert (inp["id"], inp.get("idSha256")) == (analysed["id"], analysed["idSha256"]), name


# ---------- provenance-8 ----------
def _annotation(custom, **over):
    base = {"ok": True, "accession": "GCF_000000000.1", "repeat_names": [], "genome_sha256": "a" * 64,
            "genome_bp": 200, "chunks": 1, "complete": True, "total_hits": 0, "total_bp": 0, "te_bp": 0,
            "repeatmasker_version": "4.2.4", "include_uncurated": False, "sensitivity": "default"}
    if custom:
        base.update({"custom_library": {"name": "c.fa", "sha256": "c" * 64, "sequences": 2, "bytes": 20},
                     "library_kind": "user-supplied FASTA library", "dfam_library": None, "species": None})
    else:
        base.update({"custom_library": None, "library_kind": "installed Dfam partitions, curated only",
                     "dfam_library": {"version": "4.0", "partitions": ["0 [root]: curated"]},
                     "species": "Homo sapiens"})
    base.update(over)
    return base


@pytest.mark.parametrize("custom", [True, False])
def test_genome_annotation_manifest_names_the_assembly(monkeypatch, custom):
    from teagle_core import provenance
    monkeypatch.setattr(engine.wsl, "genome_annotate", lambda *a, **k: _annotation(custom))
    body = {"assemblyAccession": "GCF_000000000.1"}
    if custom:
        body["custom_library"] = "c.fa"
    m = engine.run_genome_annotate(body)["provenance"]
    assert {"name": "RefSeq genome assembly", "version": "GCF_000000000.1", "sha256": "a" * 64} in m["databases"]
    assert m["databases"][0]["name"] == ("User-supplied repeat library" if custom else "Dfam (via famdb)")
    assert provenance.verify_manifest(m)


@pytest.mark.parametrize("backend_accession", ["GCF_999999999.1", "not an accession"])
def test_genome_annotation_refuses_to_seal_a_mismatched_assembly(monkeypatch, backend_accession):
    monkeypatch.setattr(engine.wsl, "genome_annotate",
                        lambda *a, **k: _annotation(False, accession=backend_accession))
    r = engine.run_genome_annotate({"assemblyAccession": "GCF_000000000.1"})
    assert r["ok"] is False and r["partial"] is True and "provenance" not in r


# ---------- provenance-11 ----------
@pytest.mark.parametrize("forged", ["multiple_records", "record", "record_000000000000", "P1"])
def test_top_level_input_id_cannot_be_rewritten_without_detection(forged):
    from teagle_core import provenance
    manifest = engine.analyze(">x\n" + "ACGT" * 30 + "\n")["provenance"]
    assert provenance.verify_manifest(manifest)
    tampered = copy.deepcopy(manifest)
    tampered["input"]["id"] = forged
    assert tampered["input"]["id"] != manifest["input"]["id"]
    assert not provenance.verify_manifest(tampered)


def test_generated_single_record_and_multi_record_ids_still_verify():
    from teagle_core import provenance
    bare = engine.analyze("ACGT" * 30)["provenance"]
    assert bare["input"]["id"] == "input_sequence" and provenance.verify_manifest(bare)
    multi = engine.analyze(">a\n" + "ACGT" * 30 + "\n>b\n" + "TGCA" * 30 + "\n")["provenance"]
    assert multi["input"]["id"] == "multiple_records" and provenance.verify_manifest(multi)
    renamed = copy.deepcopy(multi)
    renamed["input"]["id"] = "record"
    assert not provenance.verify_manifest(renamed)


# ---------- provenance-12 ----------
def test_software_entries_name_the_bindings_not_the_wrapped_tool():
    from teagle_core import domains, primers
    analysis = {s["name"]: s["version"] for s in engine.analyze(">r\n" + "ACGT" * 30 + "\n")["provenance"]["software"]}
    assert "HMMER (pyhmmer)" not in analysis and analysis["pyhmmer (HMMER 3 bindings)"] == domains.PYHMMER_VERSION
    primer = engine.run_primers({"sequence": "ACGT" * 100})["provenance"]
    names = {s["name"]: s["version"] for s in primer["software"]}
    assert "primer3 (primer3-py)" not in names and names["primer3-py (Primer3 bindings)"] == primers.PRIMER3_VERSION


def test_primer_manifest_records_the_libprimer3_core_release():
    from primer3.thermoanalysis import get_libprimer3_version
    names = {s["name"]: s["version"] for s in
             engine.run_primers({"sequence": "ACGT" * 100})["provenance"]["software"]}
    assert get_libprimer3_version().endswith(names["libprimer3 (Primer3 core)"])


# ---------- provenance-4 ----------
def _fetched_contexts():
    from helpers import fixture_seq
    seq = fixture_seq("M11240")[:1500]

    def ctx(raw, header):
        from teagle_core import provenance
        return provenance.normalize_source_context({**raw, "length": len(seq)}, [(header, seq)])
    ncbi1 = ctx(_accession_raw("M11240.1", "M11240.1", retrievedUtc="2026-08-21T00:00:00+00:00"), "M11240.1 x")
    ncbi2 = ctx(_accession_raw("M11240.1", "M11240.1", retrievedUtc="2026-09-30T12:34:56+00:00"), "M11240.1 x")
    ena = ctx(_accession_raw("M11240.1", "M11240.1", ena=True), "ENA|M11240|M11240.1 x")
    assert {c["inputType"] for c in (ncbi1, ncbi2, ena)} == {"fetched_accession"}
    return seq, ncbi1, ncbi2, ena


def test_pcr_seal_is_invariant_to_retrieval_time_and_serving_database():
    seq, *contexts = _fetched_contexts()
    text = f">M11240.1\n{seq}\n"
    fwd, rev = seq[:20], "TGCAACGTGCAAACGTGCAA"
    seals, bindings = set(), []
    for ctx in contexts:
        res = engine.run_pcr({"sequence": text, "source": ctx, "fwd": fwd, "rev": rev})
        seals.add(res["provenance"]["manifestSha256"])
        bindings.append(res["bindings"]["assay"])
    assert len(seals) == 1
    assert bindings[0]["retrievedUtc"] != bindings[1]["retrievedUtc"] and bindings[2]["source"] == "ena_fasta"


def test_design_bound_pcr_seal_is_invariant_to_how_the_design_was_fetched():
    seq, ncbi1, ncbi2, ena = _fetched_contexts()
    text = f">M11240.1\n{seq}\n"
    seals = set()
    for ctx in (ncbi1, ncbi2, ena):
        design = engine.run_primers({"sequence": text, "source": ctx})["candidates"][0]
        prov = design["design_provenance"]
        res = engine.run_pcr({"sequence": text, "source": ctx, "fwd": design["left_seq"], "rev": design["right_seq"],
                              "design_provenance": prov, "design_sequence": text,
                              "primer3_target": design["primer3_target"], "included_region": design["included_region"],
                              "intended_product_span": design["intended_product_span"],
                              "intended_endpoint_tolerance_bp": design["intended_endpoint_tolerance_bp"],
                              "design_assessment": "assess"})
        seals.add(res["provenance"]["manifestSha256"])
    assert len(seals) == 1


def test_pcr_seal_still_binds_the_source_identity():
    seq, ncbi1, *_ = _fetched_contexts()
    text = f">M11240.1\n{seq}\n"
    base = {"fwd": seq[:20], "rev": "TGCAACGTGCAAACGTGCAA"}
    one = engine.run_pcr({**base, "sequence": text, "source": ncbi1})["provenance"]["manifestSha256"]
    local = engine.run_pcr({**base, "sequence": text})["provenance"]["manifestSha256"]
    assert one != local


# ---------- engine_state-8 ----------
def _qt_app():
    import os
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def test_destroying_an_engine_does_not_wait_for_a_running_job():
    import threading, time
    import shiboken6
    _qt_app()
    from engine_worker import Engine
    started, release = threading.Event(), threading.Event()
    Engine.OPS["_blocking"] = lambda b: (started.set(), release.wait(30))
    try:
        eng = Engine()
        eng.submit("_blocking", {})
        assert started.wait(5)
        t0 = time.monotonic()
        shiboken6.delete(eng)
        assert time.monotonic() - t0 < 2
    finally:
        release.set()
        Engine.OPS.pop("_blocking", None)


def test_engine_still_delivers_results_and_waits_for_headless_callers():
    from PySide6.QtCore import QCoreApplication
    import time
    _qt_app()
    from engine_worker import Engine
    eng = Engine(max_threads=2)
    got = {}
    eng.done.connect(lambda k, r: got.__setitem__(k, ("done", r)))
    eng.user_error.connect(lambda k, m: got.__setitem__(k, ("user_error", m)))
    eng.failed.connect(lambda k, m, t: got.__setitem__(k, ("failed", m)))
    eng.submit("analyze", {"sequence": "ACGT" * 30}, key="ok")
    eng.submit("pcr", {"sequence": "ACGT" * 40}, key="bad")
    assert eng.wait(20000)
    deadline = time.monotonic() + 10
    while len(got) < 2 and time.monotonic() < deadline:
        QCoreApplication.processEvents()
        time.sleep(0.01)
    assert got["ok"][0] == "done" and got["bad"][0] == "user_error"


def test_a_running_job_does_not_keep_the_process_alive_after_the_application_exits():
    import os, subprocess, sys, textwrap, time
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    script = textwrap.dedent("""
        import sys, time
        sys.path[:0] = [r"%s", r"%s"]
        from PySide6.QtCore import QCoreApplication
        app = QCoreApplication([])
        from engine_worker import Engine
        Engine.OPS["_long"] = lambda b: time.sleep(60)
        eng = Engine()
        eng.submit("_long", {})
        time.sleep(0.5)
    """ % (os.path.join(root, "app", "native"), os.path.join(root, "app", "backend")))
    t0 = time.monotonic()
    subprocess.run([sys.executable, "-c", script], timeout=40, check=True,
                   env={**os.environ, "QT_QPA_PLATFORM": "offscreen"})
    assert time.monotonic() - t0 < 30


# ---------- provenance-3 ----------
def _like_json_tools(value):
    """What a JSON consumer that does not keep the int/float distinction (JavaScript, jq) writes back."""
    if isinstance(value, float) and value.is_integer():
        return int(value)
    if isinstance(value, dict):
        return {k: _like_json_tools(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_like_json_tools(v) for v in value]
    return value


def test_sealed_parameters_do_not_depend_on_how_a_number_was_typed():
    from helpers import fixture_seq
    seq = fixture_seq("M11240")[:1500]
    default = engine.run_primers({"sequence": seq})["provenance"]
    typed_int = engine.run_primers({"sequence": seq, "params": {"opt_tm": 60, "min_gc": 40}})["provenance"]
    typed_float = engine.run_primers({"sequence": seq, "params": {"opt_tm": 60.0, "min_gc": 40.0}})["provenance"]
    assert default["manifestSha256"] == typed_int["manifestSha256"] == typed_float["manifestSha256"]


@pytest.mark.parametrize("make", [
    lambda: engine.analyze(">r\n" + "ACGT" * 30 + "\n")["provenance"],
    lambda: engine.run_primers({"sequence": __import__("helpers").fixture_seq("M11240")[:1500]})["provenance"],
])
def test_manifests_survive_a_round_trip_through_number_normalising_json_tools(make):
    import json
    from teagle_core import provenance
    manifest = make()
    assert provenance.verify_manifest(json.loads(json.dumps(manifest)))
    assert provenance.verify_manifest(_like_json_tools(json.loads(json.dumps(manifest))))


def test_bare_input_design_and_pcr_stay_bound_to_one_record():
    from helpers import fixture_seq
    seq = fixture_seq("M11240")[:1500]
    design = engine.run_primers({"sequence": seq})["candidates"][0]
    res = engine.run_pcr({"sequence": seq, "fwd": design["left_seq"], "rev": design["right_seq"],
                          "design_provenance": design["design_provenance"], "design_sequence": seq,
                          "primer3_target": design["primer3_target"], "included_region": design["included_region"],
                          "intended_product_span": design["intended_product_span"],
                          "intended_endpoint_tolerance_bp": design["intended_endpoint_tolerance_bp"],
                          "design_assessment": "assess"})
    assert res["bindings"]["designAssessment"] == "assess"
    assert res["provenance"]["input"]["id"] == design["design_provenance"]["input"]["id"] == "input_sequence"


# ---------- F4: schema 1.1, preamble warning, orientation vote, sealed windows, budget size ----------
def _as_schema(doc, version, drop_summary=False):
    d = copy.deepcopy(doc)
    d["resultSchemaVersion"] = version
    if drop_summary:
        d.pop("stageSummary", None)
    d["resultPayloadSha256"] = engine.scientific_result_sha256(d)
    d["provenance"] = _resealed(d["provenance"], lambda m: m["resultDescriptor"].update(
        resultSchemaVersion=version, resultPayloadSha256=d["resultPayloadSha256"]))
    return d


def test_schema_is_1_1_and_both_it_and_a_4_0_0_shaped_1_0_document_validate():
    from helpers import make_ltr_element
    assert engine.ANALYSIS_RESULT_SCHEMA == "1.1"
    res = engine.analyze(">r\n" + make_ltr_element(seed=2) + "\n")
    doc = engine.analysis_result_document(res)
    assert doc["resultSchemaVersion"] == "1.1" and engine.validate_analysis_result_document(doc)
    legacy = _as_schema(doc, "1.0", drop_summary=True)
    assert engine.validate_analysis_result_document(legacy)
    assert not engine.validate_analysis_result_document(_as_schema(doc, "1.2"))
    assert not engine.validate_analysis_result_document(_as_schema(doc, "2.0"))


def test_domains_not_run_validates_under_1_1_and_is_refused_under_1_0_rules():
    res = engine.analyze(">short\n" + "ACGTNNACGT" * 9 + "ACGTAC\n")
    assert res["records"][0]["processing"]["domains"]["status"] == "not_run"
    doc = engine.analysis_result_document(res)
    assert engine.validate_analysis_result_document(doc)
    assert not engine.validate_analysis_result_document(_as_schema(doc, "1.0", drop_summary=True))


def test_invalid_reason_is_a_1_1_field_that_a_1_0_document_may_not_carry():
    res = engine.analyze(">e\n>ok\n" + "ACGT" * 30 + "\n")
    doc = engine.analysis_result_document(res)
    assert doc["records"][0]["invalidReason"] and engine.validate_analysis_result_document(doc)
    assert not engine.validate_analysis_result_document(_as_schema(doc, "1.0", drop_summary=True))


def test_preamble_warning_composes_with_the_multi_record_warning_before_sealing():
    s = "ACGT" * 30
    res = engine.analyze(f"JUNKTEXT12\n>a\n{s}\n>b\n{s}\n")
    assert "records found" in res["warning"]
    assert "10 characters before the first FASTA header belong to no record and were not analysed" in res["warning"]
    assert engine.verify_scientific_result_payload(res)
    assert engine.validate_analysis_result_document(engine.analysis_result_document(res))
    single = engine.analyze(f"X\n>a\n{s}\n")
    assert single["warning"] == "1 character before the first FASTA header belongs to no record and was not analysed."
    assert engine.verify_scientific_result_payload(single)
    assert "warning" not in engine.analyze(f"\n  \n>a\n{s}\n")
    assert "warning" not in engine.analyze(f">a\n{s}\n")
    bad = engine.analyze(f"LEAD\n>a\nACGTZ\n")
    assert "No valid nucleotide record" in bad["warning"] and "4 characters before" in bad["warning"]


def _undetermined_ltr_seq():
    import random
    from teagle_core import structural
    for seed in range(40):
        r = random.Random(seed)
        b = lambda n: "".join(r.choice("ACGT") for _ in range(n))
        ltr = b(400)
        seq = b(100) + ltr + b(1500) + ltr + b(100)
        row = next(e for e in structural.detect_all(seq) if e["type"].startswith("LTR (terminal"))
        if row["cis_orientation"] == "undetermined":
            return seq
    raise AssertionError("no undetermined LTR pair among 40 seeds")


def _hit(domain, strand):
    return {"domain": domain, "label": domain, "class": "x", "hmm": "x", "pfam": "x", "score": 90.0, "evalue": 1e-30,
            "confidence": "high", "orf": 0, "strand": strand, "aa": [1, 50], "nt": [300, 450], "dna": "", "protein": ""}


def _fake_domains(monkeypatch, hits):
    import inspect

    def scan_domains(*args, **kwargs):
        return list(hits)
    scan_domains.__signature__ = inspect.signature(engine.domains.scan_domains)   # the seal reads its defaults
    monkeypatch.setattr(engine.domains, "scan_domains", scan_domains)


def _spy_orientation(monkeypatch):
    calls = []
    real = engine.structural.detect_all
    monkeypatch.setattr(engine.structural, "detect_all", lambda seq, orientation=None: (calls.append(orientation), real(seq, orientation))[1])
    return calls


def _ltr_row(rec):
    return next(e for e in rec["structural"] if e["type"].startswith("LTR (terminal"))


@pytest.mark.parametrize("strand", ["+", "-"])
def test_undetermined_ltr_orientation_is_taken_from_the_domain_strand(monkeypatch, strand):
    seq = _undetermined_ltr_seq()
    _fake_domains(monkeypatch, [_hit("RT", strand), _hit("INT", strand)])
    calls = _spy_orientation(monkeypatch)
    rec = engine.analyze(f">e\n{seq}\n")["records"][0]
    assert calls == [None, strand]
    row = _ltr_row(rec)
    assert row["cis_orientation"] == strand and row["cis_orientation_source"] == "domain-strand"
    assert "cis_note" not in row
    assert rec["processing"]["structural"]["count"] == len(rec["structural"])


@pytest.mark.parametrize("hits", [
    [("RT", "+"), ("INT", "-")],
    [("GAG", "+")],
    [],
])
def test_orientation_is_not_taken_from_a_split_or_irrelevant_domain_vote(monkeypatch, hits):
    seq = _undetermined_ltr_seq()
    _fake_domains(monkeypatch, [_hit(d, st) for d, st in hits])
    calls = _spy_orientation(monkeypatch)
    row = _ltr_row(engine.analyze(f">e\n{seq}\n")["records"][0])
    assert calls == [None] and row["cis_orientation"] == "undetermined"


def test_a_cis_element_orientation_is_not_overridden_by_the_domain_strand(monkeypatch):
    from helpers import fixture_seq
    seq = fixture_seq("M11240")
    _fake_domains(monkeypatch, [_hit("RT", "-")])
    calls = _spy_orientation(monkeypatch)
    row = _ltr_row(engine.analyze(f">e\n{seq}\n")["records"][0])
    assert calls == [None] and row["cis_orientation"] == "+" and row["cis_orientation_source"] == "cis-elements"


def test_a_failed_orientation_rerun_keeps_the_first_structural_result(monkeypatch):
    seq = _undetermined_ltr_seq()
    _fake_domains(monkeypatch, [_hit("RT", "-")])
    real = engine.structural.detect_all

    def flaky(s, orientation=None):
        if orientation:
            raise RuntimeError("boom")
        return real(s)
    monkeypatch.setattr(engine.structural, "detect_all", flaky)
    res = engine.analyze(f">e\n{seq}\n")
    rec = res["records"][0]
    assert _ltr_row(rec)["cis_orientation"] == "undetermined"
    assert rec["processing"]["structural"] == {"status": "complete", "count": len(rec["structural"])}
    assert res["status"] == engine.ANALYSIS_STATUS_COMPLETE


@pytest.mark.parametrize("name", ["_TIR_SEED_WINDOW", "_TIR_JOINT_WINDOW", "_TAIL_WINDOW", "CHANCE_LEVEL_P"])
def test_new_structural_windows_and_chance_level_are_sealed_from_the_module(monkeypatch, name):
    from teagle_core import structural
    key = "structural_" + name.lstrip("_").lower()
    before = engine._detector_parameters()[key]
    assert before == getattr(structural, name)
    monkeypatch.setattr(structural, name, before * 2)
    assert engine._detector_parameters()[key] == before * 2


@pytest.mark.parametrize("key,name", [("classify_chance_tir_max_bp", "CHANCE_TIR_MAX_BP"),
                                      ("classify_chance_tir_sequences", "CHANCE_TIR_SEQUENCES")])
def test_classification_chance_reference_is_sealed_from_the_module(monkeypatch, key, name):
    from teagle_core import classify
    before = engine._detector_parameters()[key]
    assert before == getattr(classify, name)
    monkeypatch.setattr(classify, name, before + 1)
    assert engine._detector_parameters()[key] == before + 1


def test_scan_scope_note_appears_only_for_a_long_record_with_no_terminal_repeat():
    import random
    r = random.Random(3)
    long_seq = "".join(r.choice("ACGT") for _ in range(30000))
    rec = engine.analyze(f">long\n{long_seq}\n")["records"][0]
    assert not any(e["type"].startswith(("LTR (terminal", "TIR (terminal")) for e in rec["structural"])
    note = next(n for n in rec["notes"] if "terminal-repeat" in n)
    assert "30000 bp" in note and "6000" in note and "1200" in note and "600" in note
    short = engine.analyze(">s\n" + "".join(r.choice("ACGT") for _ in range(900)) + "\n")["records"][0]
    assert not any("terminal-repeat" in n for n in short["notes"])


def test_scan_scope_note_is_absent_when_a_terminal_repeat_was_found():
    import random
    r = random.Random(4)
    b = lambda n: "".join(r.choice("ACGT") for _ in range(n))
    ltr = b(400)
    rec = engine.analyze(f">e\n{b(3000)}{ltr}{b(9000)}{ltr}{b(3000)}\n")["records"][0]
    assert any(e["type"].startswith("LTR (terminal") for e in rec["structural"])
    assert not any("terminal-repeat" in n for n in rec["notes"])


def test_annotate_budget_receives_the_source_fasta_size_and_cache_state(monkeypatch):
    seen = []
    monkeypatch.setattr(engine.wsl, "annotate_budget", lambda n, **kw: seen.append((n, kw)) or {"ok": True})
    engine.run_annotate_budget({"genome_bytes": 5, "fasta_bytes": 4 * 2 ** 30, "fasta_cached": True,
                                "library_mode": "custom"})
    engine.run_annotate_budget({"genome_bytes": 5, "fasta_bytes": "123", "library_mode": "custom"})
    engine.run_annotate_budget({"genome_bytes": 5, "library_mode": "custom"})
    assert seen == [(5, {"fasta_bytes": 4 * 2 ** 30, "fasta_cached": True}),
                    (5, {"fasta_bytes": 123, "fasta_cached": None}), (5, {})]
    for bad in ({"fasta_bytes": -1}, {"fasta_bytes": True}, {"fasta_bytes": 1, "fasta_cached": "yes"}):
        with pytest.raises(BadRequest):
            engine.run_annotate_budget({"library_mode": "custom", **bad})


def _flanked_tc1(flank, seed):
    import random
    from helpers import fixture_seq
    r = random.Random(seed)
    b = lambda n: "".join(r.choice("ACGT") for _ in range(n))
    return b(flank) + "TA" + fixture_seq("X01005") + "TA" + b(flank)


def _scope_note(seq):
    rec = engine.analyze(f">r\n{seq}\n")["records"][0]
    return rec, next((n for n in rec["notes"] if "terminal-repeat" in n), None)


def test_scan_scope_note_names_the_tir_window_for_a_record_the_ltr_window_already_covers():
    seq = _flanked_tc1(1300, 2)                                    # ~4.2 kb: the LTR scan reads it all, the TIR scan does not
    rec, note = _scope_note(seq)
    assert not any(e["type"].startswith(("LTR (terminal", "TIR (terminal")) for e in rec["structural"])
    assert note is not None, "a TIR missed beyond its 1,200 bp window left no record note"
    assert "TIR" in note and "1200" in note and "LTR seeds" not in note


def test_scan_scope_note_names_both_windows_for_a_very_long_record_and_neither_for_a_short_one():
    import random
    r = random.Random(3)
    rec, note = _scope_note("".join(r.choice("ACGT") for _ in range(30000)))
    assert "LTR seeds" in note and "TIR seeds" in note and "6000" in note and "1200" in note
    _rec, none = _scope_note("".join(r.choice("ACGT") for _ in range(900)))
    assert none is None
