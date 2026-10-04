"""4.1 audit regressions for the GFF3/BED export, NCBI fetch and genome-scan summaries (lane A5)."""
import io
import json
import os
import random
import sys
import urllib.error

import pytest

_BE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "app", "backend")
if _BE not in sys.path:
    sys.path.insert(0, _BE)
import helpers                                                     # noqa: E402
import engine                                                      # noqa: E402
from teagle_core import classify, examples, fetch, genomepcr, gff3, refs   # noqa: E402


def _rows(text):
    rows = (l.split("\t") for l in text.splitlines() if l and not l.startswith("#"))
    return [r for r in rows if len(r) == 9]


def _attr(row, key):
    for part in row[8].split(";"):
        if part.startswith(key + "="):
            return part[len(key) + 1:]
    return None


def _ft(*lines):
    return fetch.parse_feature_table(">Feature gb|X.1|Y\n" + "\n".join(lines) + "\n")


def _no_pacing(monkeypatch):
    monkeypatch.setattr(fetch.time, "sleep", lambda s: None)


# ------------------------------------------------------------------ genomepcr_gel-3
def _amps(*headers):
    return genomepcr.parse_ispcr("".join(h + "\n" for h in headers))


_BANNED = ("do not amplify anywhere", "no additional priming site", "no genome-wide priming site")


def test_verdict_with_locus_and_single_primer_product_does_not_claim_no_additional_site():
    amps = _amps(">NC_1:100+300 pair 201bp AAA CCC", ">NC_9:100+500 fwdonly 401bp AAA AAA")
    amps[0]["on_target"] = True
    s = genomepcr.summarize(amps, has_locus=True)
    assert s["tier"] == "intended-only"
    assert "single-primer" in s["verdict"] and "no additional priming site" not in s["verdict"]


def test_verdict_with_locus_and_only_single_primer_products_keeps_the_hedge():
    s = genomepcr.summarize(_amps(">NC_9:100+500 fwdonly 401bp AAA AAA"), has_locus=True)
    assert s["tier"] == "none"
    assert "single-primer" in s["verdict"] and "under current settings" in s["verdict"]
    assert not any(b in s["verdict"] for b in _BANNED)


def test_verdict_with_locus_and_nothing_found_keeps_the_hedge_too():
    s = genomepcr.summarize([], has_locus=True)
    assert s["tier"] == "none" and "under current settings" in s["verdict"]
    assert "do not amplify anywhere" not in s["verdict"]


def test_verdict_without_locus_and_only_single_primer_products():
    s = genomepcr.summarize(_amps(">NC_9:100+500 revonly 401bp CCC CCC"), has_locus=False)
    assert s["tier"] == "none"
    assert "single-primer" in s["verdict"] and not any(b in s["verdict"] for b in _BANNED)


@pytest.mark.parametrize("n_off, tier", [(0, "intended-only"), (1, "limited-off-target"), (6, "many-off-target")])
def test_single_primer_clause_is_appended_without_moving_the_tier(n_off, tier):
    hdrs = [">NC_1:100+300 pair 201bp AAA CCC"] + [f">NC_{i}:100+300 pair 201bp AAA CCC" for i in range(2, 2 + n_off)]
    amps = _amps(*hdrs, ">NC_9:100+500 fwdonly 401bp AAA AAA")
    amps[0]["on_target"] = True
    s = genomepcr.summarize(amps, has_locus=True)
    assert s["tier"] == tier and "1 single-primer" in s["verdict"]


# ------------------------------------------------------------------ exports-1 / spec_conformance-16
def _flanked(seed=3, flank=2000):
    r = random.Random(seed)
    rb = lambda n: "".join(r.choice("ACGT") for _ in range(n))
    tsd = rb(5)
    return rb(flank) + tsd + helpers.fixture_seq("M11240") + tsd + rb(flank)


def _ltr(rec):
    return next(e for e in rec["structural"] if e["type"].startswith("LTR"))


def _te_row(rows, term):
    return next(r for r in rows if r[2] == term and _attr(r, "ID") == "TE_1")


def test_flanked_element_parent_spans_the_element_not_the_record():
    rec = engine.analyze(">x\n" + _flanked())["records"][0]
    a, b = _ltr(rec)["element_span"]
    assert (a, b) != (0, rec["composition"]["length"])
    rows = _rows(gff3.to_gff3(rec, seqid="x"))
    te = _te_row(rows, gff3.so_term_for(rec["classification"]))
    assert (int(te[3]), int(te[4])) == (a + 1, b)


def test_flanked_tsds_sit_in_a_wrapper_outside_the_element_and_every_child_is_contained():
    rec = engine.analyze(">x\n" + _flanked())["records"][0]
    rows = _rows(gff3.to_gff3(rec, seqid="x"))
    span = {_attr(r, "ID"): (int(r[3]), int(r[4])) for r in rows}
    te = span["TE_1"]
    tsds = [r for r in rows if r[2] == "target_site_duplication"]
    assert len(tsds) == 2
    for t in tsds:
        assert not (te[0] <= int(t[3]) and int(t[4]) <= te[1])
    for r in rows:
        parent = _attr(r, "Parent")
        if parent:
            lo, hi = span[parent]
            assert lo <= int(r[3]) and int(r[4]) <= hi, (r[2], parent)
    wrapper = next(r for r in rows if _attr(r, "Parent") is None and r[2] == "repeat_region")
    assert _attr(_te_row(rows, "LTR_retrotransposon"), "Parent") == _attr(wrapper, "ID")


def test_flanked_anchored_export_and_bed_follow_the_element_span():
    rec = engine.analyze(">x\n" + _flanked())["records"][0]
    a, b = _ltr(rec)["element_span"]
    off = 1_000_000
    g = _rows(gff3.to_gff3(rec, seqid="NC_004354.4", offset=off, identity_trust="validated_fetch"))
    te = _te_row(g, "LTR_retrotransposon")
    assert (int(te[3]), int(te[4])) == (a + 1 + off, b + off)
    bed = [l.split("\t") for l in gff3.to_bed(rec, "NC_004354.4", offset=off,
                                              identity_trust="validated_fetch").splitlines()]
    tagged = [r for r in bed if "|outcome=" in r[3]]
    assert len(tagged) == 1 and (int(tagged[0][1]), int(tagged[0][2])) == (a + off, b + off)


def test_bare_specimens_keep_the_element_as_parent_and_domain_only_keeps_the_record():
    for acc in ("M11240", "X01005"):
        rec = engine.analyze(examples.load(acc))["records"][0]
        rows = _rows(gff3.to_gff3(rec, seqid=acc))
        assert (int(rows[0][3]), int(rows[0][4])) == (1, rec["composition"]["length"]), acc
    rec = engine.analyze(examples.load("M80343"))["records"][0]
    assert not [e for e in rec["structural"] if e["type"][:3] in ("LTR", "TIR")]
    row = _rows(gff3.to_gff3(rec, seqid="M80343"))[0]
    assert (int(row[3]), int(row[4])) == (1, rec["composition"]["length"])


def test_withheld_ltr_candidate_children_use_neutral_terms_and_labels():
    st = [{"type": "LTR (terminal direct repeat)", "five_prime": [100, 300], "three_prime": [1100, 1300],
           "element_span": [100, 1300]},
          {"type": "PBS (primer-binding site)", "pos": [300, 318]},
          {"type": "PPT (polypurine tract)", "pos": [1080, 1100], "purine_frac": 0.9}]
    cl = classify.classify(st, [])
    assert cl["te_class"] == "repeat/structural-only" and cl["outcome"] == "withheld"
    rec = {"classification": cl, "composition": {"length": 1500}, "structural": st, "domains": []}
    text = gff3.to_gff3(rec, seqid="x")
    terms = [r[2] for r in _rows(text)]
    assert "long_terminal_repeat" not in terms and "RR_tract" not in terms
    assert terms.count("direct_repeat") == 2
    drs = [r for r in _rows(text) if r[2] == "direct_repeat"]
    assert all("SO:0000314" in r[8] and "Name=terminal_direct_repeat" in r[8] for r in drs)
    pbs = next(r for r in _rows(text) if r[2] == "primer_binding_site")
    assert _attr(pbs, "context") == "ltr_origin_unconfirmed"
    names = [l.split("\t")[3] for l in gff3.to_bed(rec, "x").splitlines()]
    assert not any(n.startswith("LTR|") for n in names)


def test_classified_ltr_element_still_exports_ltr_terms():
    rows = _rows(gff3.to_gff3(engine.analyze(examples.load("M11240"))["records"][0], seqid="M11240"))
    assert any(r[2] == "long_terminal_repeat" for r in rows) and any(r[2] == "RR_tract" for r in rows)


def test_withheld_without_a_direct_repeat_candidate_keeps_ltr_terms():
    # an ERV-style withheld superfamily is still an LTR retroelement; only the neutral candidate relabels
    st = [{"type": "LTR (terminal direct repeat)", "five_prime": [0, 200], "three_prime": [1000, 1200]}]
    cl = {"te_class": "LTR/Gypsy", "outcome": "withheld", "order": "withheld",
          "completeness": {"tier": "intact", "present": ["RT", "INT"], "missing": []}}
    text = gff3.to_gff3({"classification": cl, "composition": {"length": 1200}, "structural": st, "domains": []})
    assert "long_terminal_repeat" in text and "direct_repeat\t" not in text


# ------------------------------------------------------------------ exports-3
def _partial_rec(unscanned, **cl_extra):
    return {"classification": {"te_class": "none", "n_domains": 0, "orfs_unscanned": unscanned,
                               "completeness": {}, **cl_extra},
            "composition": {"length": 20000}, "structural": [], "domains": []}


def test_partial_domain_stage_is_not_exported_as_a_complete_negative():
    rec = _partial_rec(63)
    parent = _rows(gff3.to_gff3(rec))[0]
    assert _attr(parent, "domain_scan_scope") == "partial" and _attr(parent, "orfs_not_scanned") == "63"
    assert _attr(parent, "domain_evidence_status") == "tested_empty"
    name = gff3.to_bed(rec).splitlines()[0].split("\t")[3]
    assert "|domain_scan=partial|orfs_not_scanned=63" in name


def test_complete_domain_stage_carries_no_partial_marker():
    rec = _partial_rec(0)
    assert "domain_scan_scope" not in _rows(gff3.to_gff3(rec))[0][8]
    assert "domain_scan=" not in gff3.to_bed(rec).splitlines()[0]


def test_unavailable_domain_scan_is_not_also_labelled_partial():
    rec = _partial_rec(5, domains_unavailable=True)
    assert "domain_scan_scope" not in _rows(gff3.to_gff3(rec))[0][8]


def test_bundled_partial_specimen_exports_the_partial_marker():
    rec = engine.analyze(examples.load("M12927"))["records"][0]
    assert rec["processing"]["domains"]["status"] == "partial"
    assert "orfs_not_scanned=" in _rows(gff3.to_gff3(rec, seqid="M12927"))[0][8]
    assert "|domain_scan=partial|" in gff3.to_bed(rec, "M12927")


# ------------------------------------------------------------------ exports-9
def test_dna_transposon_with_unassessed_tirs_degrades_to_the_generic_term():
    d = classify.classify(structural=[], structural_ok=False,
                          domains=[{"domain": "TPase", "class": "hAT", "score": 200, "nt": [100, 900], "strand": "+"}])
    assert any("TIR" in m for m in d["completeness"].get("unassessed", []))
    assert gff3.so_term_for(d) == "DNA_transposon"
    row = _rows(gff3.to_gff3({"classification": d, "composition": {"length": 1000}, "structural": [],
                              "domains": []}))[0]
    assert row[2] == "DNA_transposon"


def test_dna_transposon_with_recovered_tirs_keeps_the_tir_element_term():
    rec = engine.analyze(examples.load("X01005"))["records"][0]
    assert gff3.so_term_for(rec["classification"]) == "terminal_inverted_repeat_element"


# ------------------------------------------------------------------ gene model: fetch_refs-1, exports-2/11/12
_M80343_FT = (">Feature gb|M80343.1|HUMTNL22\n1\t909\t5'UTR\n911\t1927\tCDS\n\t\t\tprot_desc\tORF1\n"
              "\t\t\tprotein_id\tgb|AAB59367.1|\n1991\t5818\tCDS\n\t\t\tprot_desc\tORF2\n"
              "\t\t\tprotein_id\tgb|AAB59368.1|\n5818\t6049\t3'UTR\n")


def test_separate_unspliced_cds_features_give_no_introns():
    gm = fetch.build_gene_model(fetch.parse_feature_table(_M80343_FT))
    assert gm["introns"] == [] and gm["counts"]["introns"] == 0 and gm["derived_introns"] is False
    bed = gff3.gene_model_to_bed(gm, "M80343.1", sequence_length=6050)
    assert "\tintron_" not in bed


def test_gypsy_style_orf_cds_features_give_no_introns():
    gm = fetch.build_gene_model(_ft("1080\t2435\tCDS", "\t\t\tproduct\ta", "2438\t5470\tCDS", "\t\t\tproduct\tb",
                                    "5561\t7000\tCDS", "\t\t\tproduct\tc"))
    assert gm["introns"] == []


def test_two_genes_on_opposite_strands_give_only_within_gene_introns():
    gm = fetch.build_gene_model(_ft("101\t200\texon", "\t\t\tgene\tgeneA", "301\t400\texon", "\t\t\tgene\tgeneA",
                                    "2000\t1901\texon", "\t\t\tgene\tgeneB", "2400\t2301\texon", "\t\t\tgene\tgeneB"))
    assert {(i["start"], i["end"], i["strand"]) for i in gm["introns"]} == {(200, 300, "+"), (2000, 2300, "-")}
    assert all(i.get("derived") for i in gm["introns"])


def test_two_mrna_joins_for_different_genes_give_no_intergenic_intron():
    gm = fetch.build_gene_model(_ft("100\t200\tmRNA", "300\t400", "\t\t\tgene\tA", "1000\t1100\tmRNA", "1300\t1400",
                                    "\t\t\tgene\tB"))
    assert {(i["start"], i["end"]) for i in gm["introns"]} == {(200, 299), (1100, 1299)}


def test_single_joined_cds_still_yields_its_intron_marked_derived():
    gm = fetch.build_gene_model(_ft("100\t200\tCDS", "300\t400", "\t\t\tcodon_start\t1"))
    assert [(i["start"], i["end"], i.get("derived")) for i in gm["introns"]] == [(200, 299, True)]
    text = gff3.gene_model_to_gff3(gm, "chr1", sequence="A" * 500)
    assert "derived=true" in next(l for l in text.splitlines() if "\tintron\t" in l)


def test_completing_a_two_unit_model_creates_no_cross_unit_intron():
    gm = fetch.build_gene_model(_ft("101\t200\texon", "\t\t\tgene\tA", "501\t600\texon", "\t\t\tgene\tA",
                                    "2001\t2100\texon", "\t\t\tgene\tB", "2301\t2400\texon", "\t\t\tgene\tB",
                                    "251\t300\tCDS", "\t\t\tcodon_start\t1"))
    assert {(i["start"], i["end"]) for i in gm["introns"]} == {(200, 250), (300, 500), (2100, 2300)}


def test_absent_codon_start_means_frame_one_and_gff3_exports():
    gm = fetch.build_gene_model(fetch.parse_feature_table(_M80343_FT))
    assert [c["phase"] for c in gm["cds"]] == [0, 0]
    text = gff3.gene_model_to_gff3(gm, "M80343.1", sequence="A" * 6050)
    assert "\tCDS\t911\t1927\t.\t+\t0\t" in text


def test_two_segment_cds_without_codon_start_carries_the_running_phase():
    gm = fetch.build_gene_model(_ft("2424\t2610\tCDS", "3397\t3542", "\t\t\tproduct\tinsulin"))
    assert [c["phase"] for c in gm["cds"]] == [0, 2]


def test_minus_strand_cds_without_codon_start_starts_in_frame_one():
    gm = fetch.build_gene_model(_ft("400\t301\tCDS", "200\t101", "\t\t\tproduct\tp"))
    assert {c["start"]: c["phase"] for c in gm["cds"]} == {300: 0, 100: 2}


@pytest.mark.parametrize("value", ["0", "4", "x", ""])
def test_present_but_invalid_codon_start_is_still_unavailable(value):
    gm = fetch.build_gene_model(fetch.parse_feature_table(f">Feature gb|X.1|Y\n1\t9\tCDS\n\t\t\tcodon_start\t{value}\n"))
    assert gm["cds"][0]["phase"] is None


def _cds_only_model():
    return fetch.build_gene_model(_ft("101\t200\tCDS", "\t\t\tcodon_start\t1", "301\t400\tCDS", "\t\t\tcodon_start\t1"))


def test_derived_exons_are_marked_in_bed_and_fasta_but_annotated_ones_are_not():
    gm = _cds_only_model()
    seq = "ACGT" * 110
    bed = [l.split("\t") for l in gff3.gene_model_to_bed(gm, "chr1", sequence_length=len(seq)).splitlines()]
    exons = [b for b in bed if b[3].startswith("exon_")]
    assert exons and all("derived=true" in b[3] for b in exons)
    assert not any("derived=true" in b[3] for b in bed if b[3].startswith("CDS_"))
    fa = gff3.gene_model_to_fasta(gm, seq, "chr1", feature_types=("exon", "CDS"))
    heads = [l for l in fa.splitlines() if l.startswith(">")]
    assert all("derived=true" in h for h in heads if ".exon" in h)
    assert not any("derived=true" in h for h in heads if ".cds" in h)
    annotated = fetch.build_gene_model(_ft("101\t200\texon", "\t\t\tgene\tA"))
    assert "derived" not in gff3.gene_model_to_bed(annotated, "chr1", sequence_length=500)


def test_derived_exon_rows_do_not_inherit_cds_only_attributes():
    text = gff3.gene_model_to_gff3(_cds_only_model(), "chr1", sequence="A" * 500)
    for r in _rows(text):
        if r[2] == "exon":
            assert _attr(r, "derived") == "true" and _attr(r, "cds_group") is None and r[7] == "."


def test_segments_of_one_cds_share_an_id_and_name_and_distinct_cds_do_not():
    model = {"exons": [{"start": 0, "end": 30, "strand": "+"}],
             "cds": [{"start": 1, "end": 6, "strand": "+", "phase": 0, "cds_group": "cds_1"},
                     {"start": 11, "end": 18, "strand": "+", "phase": 2, "cds_group": "cds_1"},
                     {"start": 12, "end": 19, "strand": "+", "phase": 0, "cds_group": "cds_2"}]}
    cds = [r for r in _rows(gff3.gene_model_to_gff3(model, "chr1", sequence="A" * 30)) if r[2] == "CDS"]
    one = [r for r in cds if _attr(r, "cds_group") == "cds_1"]
    two = [r for r in cds if _attr(r, "cds_group") == "cds_2"]
    assert len({_attr(r, "ID") for r in one}) == 1 and len({_attr(r, "Name") for r in one}) == 1
    assert _attr(two[0], "ID") != _attr(one[0], "ID")
    assert [r[7] for r in one] == ["0", "2"]
    exon_ids = [_attr(r, "ID") for r in _rows(gff3.gene_model_to_gff3(model, "chr1", sequence="A" * 30))
                if r[2] == "exon"]
    assert not set(exon_ids) & {_attr(r, "ID") for r in cds}


def test_cds_without_a_group_keeps_unique_per_row_ids():
    model = {"cds": [{"start": 1, "end": 6, "strand": "+", "phase": 0},
                     {"start": 11, "end": 18, "strand": "+", "phase": 2}]}
    ids = [_attr(r, "ID") for r in _rows(gff3.gene_model_to_gff3(model, "chr1", sequence="A" * 30))]
    assert len(ids) == len(set(ids)) == 2


def test_complete_gene_model_membership_check_is_not_quadratic():
    import time
    exons = [{"start": i * 1000, "end": i * 1000 + 500, "strand": "+"} for i in range(8000)]
    cds = [{"start": i * 1000 + 10, "end": i * 1000 + 400, "strand": "+"} for i in range(8000)]
    t = time.perf_counter()
    fetch.complete_gene_model({"exons": exons, "introns": [], "cds": cds, "derived_exons": False,
                               "derived_introns": False, "counts": {}})
    assert time.perf_counter() - t < 1


# ------------------------------------------------------------------ fetch_refs-2
def _stub_custom_store(monkeypatch, tmp_path, *resolutions):
    monkeypatch.setattr(fetch, "_CUSTOM_ASM_PATH", str(tmp_path / "c.json"))
    queue = list(resolutions)
    monkeypatch.setattr(fetch, "resolve_assembly", lambda q: queue.pop(0))


_T2T = {"organism": "Homo sapiens", "taxid": "9606", "assemblyName": "T2T-CHM13v2.0",
        "assemblyAccession": "GCF_009914755.1"}


def _sim(acc):
    return {"organism": "Drosophila simulans", "taxid": "7240", "assemblyName": "x", "assemblyAccession": acc}


def test_adding_an_alternative_assembly_of_a_curated_organism_is_refused(monkeypatch, tmp_path):
    _stub_custom_store(monkeypatch, tmp_path, _T2T)
    with pytest.raises(fetch.CoordError, match="already pinned"):
        fetch.add_custom_assembly("GCF_009914755.1")
    assert fetch.load_custom_assemblies() == {}


def test_second_assembly_of_a_custom_organism_is_refused_and_keeps_the_first(monkeypatch, tmp_path):
    _stub_custom_store(monkeypatch, tmp_path, _sim("GCF_016746395.2"), _sim("GCA_000754195.3"))
    fetch.add_custom_assembly("GCF_016746395.2")
    with pytest.raises(fetch.CoordError):
        fetch.add_custom_assembly("GCA_000754195.3")
    assert fetch.load_custom_assemblies()["Drosophila simulans"]["assemblyAccession"] == "GCF_016746395.2"


def test_re_adding_the_same_accession_is_idempotent_and_new_organisms_are_added(monkeypatch, tmp_path):
    _stub_custom_store(monkeypatch, tmp_path, _sim("GCF_016746395.2"), _sim("GCF_016746395.2"))
    fetch.add_custom_assembly("a")
    assert fetch.add_custom_assembly("a")["assemblyAccession"] == "GCF_016746395.2"
    assert fetch.all_assemblies()["Drosophila simulans"]["custom"] is True


def test_re_adding_the_curated_accession_is_a_noop(monkeypatch, tmp_path):
    same = {**_T2T, "assemblyName": "GRCh38.p14", "assemblyAccession": "GCF_000001405.40"}
    _stub_custom_store(monkeypatch, tmp_path, same)
    fetch.add_custom_assembly("GCF_000001405.40")
    assert fetch.load_custom_assemblies() == {}


# ------------------------------------------------------------------ fetch_refs-3 / -14
def _page(chrs, token=None):
    reps = [{"role": "assembled-molecule", "chr_name": c, "ucsc_style_name": u, "refseq_accession": a,
             "length": 100} for c, u, a in chrs]
    return {"reports": reps, **({"next_page_token": token} if token else {})}


def test_live_assembly_map_follows_next_page_token(monkeypatch):
    calls = []
    pages = {None: _page([("1", "chr1", "NC_1.1"), ("2", "chr2", "NC_2.1")], "tok"),
             "tok": _page([("MT", "chrM", "NC_MT.1")])}

    def fake(path, *a, **k):
        calls.append(path)
        return pages["tok" if "page_token=tok" in path else None]
    monkeypatch.setattr(fetch, "_datasets_json", fake)
    doc = fetch._fetch_assembly_map_live("GCF_TEST.1")
    assert [m["refseqAccession"] for m in doc["molecules"]] == ["NC_1.1", "NC_2.1", "NC_MT.1"]
    assert len(calls) == 2 and "page_token=tok" in calls[1]


def test_build_script_sequence_report_follows_next_page_token(monkeypatch):
    import urllib.request
    import teagle_core.build_assembly_maps as bam
    seen = []

    class R(io.BytesIO):
        def __enter__(self): return self
        def __exit__(self, *a): return False

    def fake_urlopen(req, timeout=0):
        seen.append(req.full_url)
        body = _page([("MT", "chrM", "NC_MT.1")]) if "page_token=tok" in req.full_url else \
            _page([("1", "chr1", "NC_1.1")], "tok")
        return R(json.dumps(body).encode())
    monkeypatch.setattr(urllib.request, "urlopen", fake_urlopen)
    reps = bam.sequence_report("GCF_TEST.1")
    assert [r["refseq_accession"] for r in reps] == ["NC_1.1", "NC_MT.1"] and len(seen) == 2


def test_stale_cached_custom_map_without_the_pagination_stamp_is_rebuilt(monkeypatch, tmp_path):
    monkeypatch.setattr(fetch, "CACHE_DIR", str(tmp_path))
    monkeypatch.setattr(fetch, "_assembly_map_dir", lambda: str(tmp_path / "none"))
    fetch._write_cache("asm_GCF_TEST.1", {"assemblyAccession": "GCF_TEST.1",
                                          "molecules": [{"chrName": "1", "refseqAccession": "NC_1.1", "length": 5}]})
    rebuilt = {"assemblyAccession": "GCF_TEST.1", "molecules": [
        {"chrName": "1", "refseqAccession": "NC_1.1", "length": 5},
        {"chrName": "MT", "ucscStyleName": "chrM", "refseqAccession": "NC_MT.1", "length": 5}]}
    monkeypatch.setattr(fetch, "_fetch_assembly_map_live", lambda acc: dict(rebuilt))
    assert fetch.resolve_chrom("GCF_TEST.1", "chrM")["refseqAccession"] == "NC_MT.1"


def test_bundled_large_assembly_maps_include_the_mitochondrion():
    assert fetch.resolve_chrom("GCF_000001215.4", "chrM")["refseqAccession"] == "NC_024511.2"
    assert fetch.resolve_chrom("GCF_000002035.6", "chrM")["refseqAccession"] == "NC_002333.2"


_PLANT = {"assemblyAccession": "GCF_TEST.1", "molecules": [
    {"chrName": "1", "ucscStyleName": "", "refseqAccession": "NC_1.1", "length": 10},
    {"chrName": "Pltd", "ucscStyleName": "", "refseqAccession": "NC_PL.1", "length": 10},
    {"chrName": "MT", "ucscStyleName": "chrM", "refseqAccession": "NC_MT.1", "length": 10}]}


@pytest.mark.parametrize("name", ["Pt", "chrC", "ChrC", "C", "chloroplast", "Pltd"])
def test_plastid_synonyms_resolve(monkeypatch, name):
    monkeypatch.setattr(fetch, "_load_assembly_map", lambda acc, refresh=False: _PLANT)
    assert fetch.resolve_chrom("GCF_TEST.1", name)["refseqAccession"] == "NC_PL.1"


def test_an_exact_chromosome_name_is_never_shadowed_by_an_alias(monkeypatch):
    amap = {"assemblyAccession": "GCF_TEST.1", "molecules": [
        {"chrName": "C", "refseqAccession": "NC_C.1", "length": 10},
        {"chrName": "Pltd", "refseqAccession": "NC_PL.1", "length": 10}]}
    monkeypatch.setattr(fetch, "_load_assembly_map", lambda acc, refresh=False: amap)
    assert fetch.resolve_chrom("GCF_TEST.1", "C")["refseqAccession"] == "NC_C.1"


def test_unknown_chromosome_error_lists_the_available_names():
    with pytest.raises(fetch.CoordError, match="2L"):
        fetch.resolve_chrom("GCF_000001215.4", "bogus99")


# ------------------------------------------------------------------ fetch_refs-4
def _esummary(acc_ver, slen, status="", replacedby=""):
    rec = {"accessionversion": acc_ver, "organism": "Caenorhabditis elegans", "taxid": 6239, "title": "Tc1",
           "slen": slen, "moltype": "dna", "status": status, "replacedby": replacedby}
    return json.dumps({"result": {"uids": ["1"], "1": rec}})


def _ena_world(monkeypatch, tmp_path, header, seq, esummary=None):
    _no_pacing(monkeypatch)
    monkeypatch.setattr(fetch, "CACHE_DIR", str(tmp_path))
    urls = []

    def fake_get(url, *a, **k):
        urls.append(url)
        if "esummary" in url:
            return esummary or _esummary("X01005.1", len(seq))
        if "rettype=ft" in url:
            return ""
        if "ebi.ac.uk" in url:
            return header + "\n" + "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)) + "\n"
        raise fetch.FetchError("HTTP 503 from source")
    monkeypatch.setattr(fetch, "_get", fake_get)
    return urls


@pytest.mark.parametrize("typed", ["X01005", "X01005.1"])
def test_ena_fallback_yields_fetched_accession_provenance(monkeypatch, tmp_path, typed):
    from teagle_core import provenance, sequtil
    seq = helpers.fixture_seq("X01005")
    urls = _ena_world(monkeypatch, tmp_path, ">ENA|X01005|X01005.1 Caenorhabditis elegans Tc1", seq)
    out = engine.run_fetch({"accession": typed, "refresh": True})
    ctx = out["sourceContext"]
    assert ctx["inputType"] == "fetched_accession" and ctx["source"] == "ena_fasta"
    assert ctx["accession"] == "X01005.1"
    assert out["fasta"].startswith(">X01005.1")
    assert any(u.endswith("/X01005.1") for u in urls if "ebi.ac.uk" in u)
    assert provenance.validate_source_context(ctx, sequtil.parse_fasta(out["fasta"]))["inputType"] == "fetched_accession"


def test_ena_record_naming_a_different_accession_stays_local_input(monkeypatch, tmp_path):
    seq = helpers.fixture_seq("X01005")
    _ena_world(monkeypatch, tmp_path, ">ENA|X01005|X01006.1 other", seq)
    out = engine.run_fetch({"accession": "X01005", "refresh": True})
    assert out["sourceContext"]["inputType"] == "local_input"
    assert out["fasta"].startswith(">ENA|X01005|X01006.1")


def test_unversioned_typed_accession_still_seals_as_fetched_from_ncbi(monkeypatch, tmp_path):
    _no_pacing(monkeypatch)
    monkeypatch.setattr(fetch, "CACHE_DIR", str(tmp_path))
    seq = helpers.fixture_seq("X01005")

    def fake_get(url, *a, **k):
        if "esummary" in url:
            return _esummary("X01005.1", len(seq))
        if "rettype=ft" in url:
            return ""
        return ">X01005.1 Tc1\n" + seq + "\n"
    monkeypatch.setattr(fetch, "_get", fake_get)
    out = engine.run_fetch({"accession": "X01005", "refresh": True})
    assert out["sourceContext"]["inputType"] == "fetched_accession"
    assert out["sourceUrl"].endswith("/X01005.1")


# ------------------------------------------------------------------ fetch_refs-5 / -16
def test_resolve_surfaces_replaced_status_and_live_records_stay_clean(monkeypatch):
    monkeypatch.setattr(fetch, "_get", lambda *a, **k: _esummary("M11240.1", 5146, "replaced", "X04456"))
    m = fetch.resolve("M11240")
    assert m["status"] == "replaced" and m["replacedBy"] == "X04456"
    monkeypatch.setattr(fetch, "_get", lambda *a, **k: _esummary("X01005.1", 1610))
    assert not fetch.resolve("X01005").get("status")


def test_replaced_status_stays_out_of_the_seal(monkeypatch, tmp_path):
    from teagle_core import provenance, sequtil
    _no_pacing(monkeypatch)
    monkeypatch.setattr(fetch, "CACHE_DIR", str(tmp_path))
    seq = helpers.fixture_seq("X01005")

    def fake_get(url, *a, **k):
        if "esummary" in url:
            return _esummary("X01005.1", len(seq), "replaced", "X99999")
        if "rettype=ft" in url:
            return ""
        return ">X01005.1 Tc1\n" + seq + "\n"
    monkeypatch.setattr(fetch, "_get", fake_get)
    out = engine.run_fetch({"accession": "X01005", "refresh": True})
    assert out["status"] == "replaced"
    assert "status" not in json.dumps(out["sourceContext"]) and "X99999" not in json.dumps(out["sourceContext"])


def test_copia_and_line_examples_carry_the_replacement_notes():
    cat = {a: t for a, _l, _o, t in examples.CATALOG}
    assert "X04456" in cat["M11240"] and "AH005269" in cat["M80343"]
    assert "ORF1" in cat["M80343"] and "RT, no integrase" not in cat["M80343"]


def test_gypsy_example_does_not_promise_a_retroviral_env_call():
    rec = engine.analyze(examples.load("M12927"))["records"][0]
    cl = rec["classification"]
    assert cl["te_class"].startswith("LTR/Gypsy") and not cl.get("is_erv")
    text = {a: t for a, _l, _o, t in examples.CATALOG}["M12927"]
    assert "env-bearing errantivirus" not in text and "not an ERV call" in text


# ------------------------------------------------------------------ fetch_refs-6 / -7
def test_error_example_is_a_te_sized_record_not_a_chromosome():
    with pytest.raises(fetch.FetchError) as e:
        fetch.validate_accession("not!an!acc")
    assert "NC_003075.7" not in str(e.value)
    for f in os.listdir(fetch._assembly_map_dir()):
        mols = json.load(open(os.path.join(fetch._assembly_map_dir(), f), encoding="utf-8"))["molecules"]
        assert not any(m["refseqAccession"] in str(e.value) for m in mols)


@pytest.mark.parametrize("acc", ["NZ_CP012345.1", "NZ_JAAGXY010000001.1", "m11240"])
def test_nz_and_lowercase_accessions_validate(acc):
    assert fetch.validate_accession(acc) == acc


def test_bare_nz_prefix_is_still_invalid():
    with pytest.raises(fetch.FetchError):
        fetch.validate_accession("NZ_")


@pytest.mark.parametrize("acc, word", [("NP_000509.1", "protein"), ("WP_012345678.1", "protein"),
                                       ("GCF_000001405.40", "coordinate"), ("GCA_000001405.29", "coordinate")])
def test_protein_and_assembly_accessions_get_a_targeted_message(acc, word):
    with pytest.raises(fetch.FetchError, match=word):
        fetch.validate_accession(acc)


# ------------------------------------------------------------------ fetch_refs-8
def _stub_retrieve(monkeypatch, tmp_path):
    monkeypatch.setattr(fetch, "CACHE_DIR", str(tmp_path))
    monkeypatch.setattr(fetch, "resolve", lambda acc: {"accession": "M12927.1", "organism": "x", "taxid": 1,
                                                       "title": "t", "length": 80, "moltype": "DNA",
                                                       "source": "NCBI nuccore"})
    monkeypatch.setattr(fetch, "fetch_fasta", lambda acc, served=None: ">M12927.1\n" + "ACGT" * 20)
    _no_pacing(monkeypatch)


def test_a_failed_feature_fetch_is_retried_on_the_next_cache_hit(monkeypatch, tmp_path):
    _stub_retrieve(monkeypatch, tmp_path)
    good = fetch.build_gene_model(_ft("2\t30\tCDS", "\t\t\tcodon_start\t1"))
    calls = []

    def feats(acc):
        calls.append(acc)
        if len(calls) == 1:
            raise fetch.FetchError("HTTP 429 from source")
        return good
    monkeypatch.setattr(fetch, "fetch_features", feats)
    first = fetch.retrieve("M12927")
    assert "features" not in first and first["featuresStatus"] == "failed"
    second = fetch.retrieve("M12927")
    assert second["fromCache"] is True and "features" in second and len(calls) == 2
    cached = fetch._read_cache("M12927")
    assert cached["features"] and cached.get("fromCache") is not True and "servedUtc" not in cached


def test_a_record_without_annotation_is_not_re_requested(monkeypatch, tmp_path):
    _stub_retrieve(monkeypatch, tmp_path)
    calls = []
    monkeypatch.setattr(fetch, "fetch_features", lambda acc: calls.append(acc) or fetch.build_gene_model([]))
    assert fetch.retrieve("M12927")["featuresStatus"] == "none"
    fetch.retrieve("M12927")
    assert len(calls) == 1


def test_a_legacy_cache_with_an_old_gene_model_is_rebuilt_once(monkeypatch, tmp_path):
    _stub_retrieve(monkeypatch, tmp_path)
    fetch._write_cache("M12927", {"accession": "M12927.1", "fasta": ">M12927.1\nACGT", "retrievedUtc": "2026-01-01T00:00:00+00:00",
                                  "features": {"counts": {"exons": 1, "introns": 1, "cds": 1},
                                               "cds": [{"start": 1, "end": 9, "strand": "+"}]}})
    good = fetch.build_gene_model(_ft("2\t30\tCDS"))
    calls = []
    monkeypatch.setattr(fetch, "fetch_features", lambda acc: calls.append(acc) or good)
    m = fetch.retrieve("M12927")
    assert calls == ["M12927"] and m["features"]["cds"][0]["phase"] == 0
    assert m["retrievedUtc"] == "2026-01-01T00:00:00+00:00" and m["fasta"].startswith(">M12927.1")
    fetch.retrieve("M12927")
    assert len(calls) == 1


def test_stale_gene_model_that_cannot_be_rebuilt_offline_is_dropped_but_the_fasta_is_served(monkeypatch, tmp_path):
    _stub_retrieve(monkeypatch, tmp_path)
    fetch._write_cache("M12927", {"accession": "M12927.1", "fasta": ">M12927.1\nACGT",
                                  "features": {"counts": {"exons": 1, "introns": 0, "cds": 0}}})
    monkeypatch.setattr(fetch, "fetch_features", lambda acc: (_ for _ in ()).throw(fetch.FetchError("offline")))
    m = fetch.retrieve("M12927")
    assert m["fasta"].startswith(">M12927.1") and "features" not in m


# ------------------------------------------------------------------ fetch_refs-9
@pytest.mark.parametrize("key", ["Dfam", "RepeatMasker"])
def test_annotate_references_carry_no_version_or_partition_claim(key):
    import re
    name = refs.REFS[key]["name"]
    assert not re.search(r"curated|\d+\.\d+", name), name
    assert key in {r["key"] for r in refs.for_run("annotate")}


# ------------------------------------------------------------------ fetch_refs-10
class _Resp(io.BytesIO):
    def __enter__(self): return self
    def __exit__(self, *a): return False


def _http_error(code, headers=None):
    return urllib.error.HTTPError("http://x", code, "err", headers or {}, None)


def _scripted_urlopen(monkeypatch, outcomes):
    seen = []

    def fake(req, timeout=0, context=None):
        seen.append(req.full_url)
        o = outcomes.pop(0)
        if isinstance(o, Exception):
            raise o
        return _Resp(o)
    monkeypatch.setattr(fetch.urllib.request, "urlopen", fake)
    sleeps = []
    monkeypatch.setattr(fetch.time, "sleep", lambda s: sleeps.append(s))
    return seen, sleeps


def test_get_retries_rate_limit_and_server_errors_then_succeeds(monkeypatch):
    seen, sleeps = _scripted_urlopen(monkeypatch, [_http_error(429, {"Retry-After": "2"}), _http_error(503), b"ok"])
    assert fetch._get("https://eutils.ncbi.nlm.nih.gov/x") == "ok"
    assert len(seen) == 3 and 2 in sleeps


def test_get_gives_up_after_two_retries_with_a_plain_message(monkeypatch):
    seen, _ = _scripted_urlopen(monkeypatch, [_http_error(429)] * 3)
    with pytest.raises(fetch.FetchError, match="rate-limiting"):
        fetch._get("https://eutils.ncbi.nlm.nih.gov/x")
    assert len(seen) == 3


def test_get_does_not_retry_a_client_error(monkeypatch):
    seen, _ = _scripted_urlopen(monkeypatch, [_http_error(404)])
    with pytest.raises(fetch.FetchError, match="404"):
        fetch._get("https://eutils.ncbi.nlm.nih.gov/x")
    assert len(seen) == 1


def test_get_explains_a_connectivity_failure_without_retrying(monkeypatch):
    _scripted_urlopen(monkeypatch, [urllib.error.URLError("getaddrinfo failed")])
    with pytest.raises(fetch.FetchError, match="internet connection"):
        fetch._get("https://eutils.ncbi.nlm.nih.gov/x")


def test_ncbi_requests_are_paced(monkeypatch):
    seen, sleeps = _scripted_urlopen(monkeypatch, [b"a", b"b"])
    monkeypatch.setattr(fetch, "_last_ncbi", [fetch.time.monotonic()])
    fetch._get("https://eutils.ncbi.nlm.nih.gov/a")
    assert sleeps and 0 < sleeps[0] <= fetch._NCBI_GAP + 0.01


# ------------------------------------------------------------------ fetch_refs-12
def test_coordinate_cache_hit_keeps_the_original_retrieval_time(monkeypatch, tmp_path):
    monkeypatch.setattr(fetch, "CACHE_DIR", str(tmp_path))
    monkeypatch.setattr(fetch, "resolve_chrom", lambda acc, key, refresh=False: {
        "refseqAccession": "NC_000013.11", "length": 114364328, "ucscStyleName": "chr13", "chrName": "13"})
    monkeypatch.setattr(fetch, "fetch_fasta_range",
                        lambda acc, s, e, strand="+": f">{acc}:{s}-{e}\n" + "ACGT" * ((e - s + 4) // 4))
    first = fetch.retrieve_coords("chr13:1-100", "GCF_TEST.1", "GRCh38", "Homo sapiens")
    t0 = first["source"]["retrievedUtc"]
    monkeypatch.setattr(fetch, "_now", lambda: "2099-01-01T00:00:00+00:00")
    hit = fetch.retrieve_coords("chr13:1-100", "GCF_TEST.1", "GRCh38", "Homo sapiens")
    assert hit["fromCache"] is True and hit["source"]["retrievedUtc"] == t0
    assert hit["servedUtc"] == "2099-01-01T00:00:00+00:00"


# ------------------------------------------------------------------ fetch_refs-13
@pytest.mark.parametrize("raw, expected", [
    ("GAAA01000001.1", ["GAAA01000001.1"]),
    ("gnl|WGS:JAABCD|FUN_000123-T1", []),
    ("gb|AB123456.1|", ["AB123456.1"]),
    ("ref|NM_000207.3|", ["NM_000207.3"]),
    ("NM_000207.3", ["NM_000207.3"]),
    ("lcl|local_id_7", []),
])
def test_transcript_picker_accessions_are_anchored(raw, expected):
    gm = fetch.build_gene_model(_ft("1\t100\tmRNA", f"\t\t\ttranscript_id\t{raw}", "\t\t\tgene\tg"))
    assert [t["accession"] for t in gm["transcripts"]] == expected


def test_structural_stage_failure_is_declared_on_the_parent():
    d = classify.classify([], [], structural_ok=False)
    parent = _rows(gff3.to_gff3({"classification": d, "composition": {"length": 500}, "structural": [],
                                 "domains": []}))[0]
    assert _attr(parent, "structural_evidence_status") == "unavailable_failed"


def test_cache_written_from_the_typed_accession_is_repaired_on_read(monkeypatch, tmp_path):
    monkeypatch.setattr(fetch, "CACHE_DIR", str(tmp_path))
    fetch._write_cache("X01005", {"accession": "X01005.1", "source": "ENA (EMBL-EBI)", "fasta": ">X01005.1\nACGT",
                                  "endpoint": fetch.ENA + "X01005",
                                  "sourceUrl": "https://www.ebi.ac.uk/ena/browser/view/X01005",
                                  "featuresStatus": "none", "featuresModel": fetch.GENE_MODEL_VERSION})
    hit = fetch.retrieve("X01005")
    assert hit["endpoint"] == fetch.ENA + "X01005.1" and hit["sourceUrl"].endswith("/X01005.1")


def test_wrapper_covers_only_the_element_and_its_tsds_never_other_flank_features():
    rec = engine.analyze(">x\n" + _flanked())["records"][0]
    tsd = next(e for e in rec["structural"] if e["type"].startswith("TSD"))
    rec["domains"] = list(rec["domains"]) + [{"domain": "RT", "nt": [100, 400], "strand": "+", "score": 50.0}]
    rows = _rows(gff3.to_gff3(rec, seqid="x"))
    wrapper = next(r for r in rows if r[2] == "repeat_region")
    assert (int(wrapper[3]), int(wrapper[4])) == (tsd["upstream"][0] + 1, tsd["downstream"][1])
    flank_domain = next(r for r in rows if r[2] == "polypeptide_domain" and int(r[3]) == 101)
    assert _attr(flank_domain, "Parent") is None


def test_uncredited_tir_does_not_bound_the_element():
    st = [{"type": "TIR (terminal inverted repeat)", "five_prime": [2000, 2030], "three_prime": [2900, 2930],
           "element_span": [2000, 2930]}]
    dm = [{"domain": "TPase", "class": "hAT", "score": 200.0, "nt": [100, 900], "strand": "+"}]
    cl = classify.classify(st, dm)
    assert cl["te_class"].startswith("DNA/") and any("TIR" in m for m in cl["completeness"]["missing"])
    row = _rows(gff3.to_gff3({"classification": cl, "composition": {"length": 3000}, "structural": st,
                              "domains": dm}))[0]
    assert (int(row[3]), int(row[4])) == (1, 3000)


# ------------------------------------------------------------------ F5: export carries the classification fields
def _parent(rec):
    return next(r for r in _rows(gff3.to_gff3(rec, seqid="x")) if _attr(r, "ID") == "TE_1")


def test_line_export_names_wicker_order_and_domain_architecture_separately():
    rec = engine.analyze(examples.load("M80343"))["records"][0]
    cl = rec["classification"]
    row = _parent(rec)
    assert _attr(row, "wicker_order") == "LINE"
    assert _attr(row, "domain_architecture") == cl["domain_architecture"] != "unassigned"
    assert _attr(row, "classification_order") == cl["domain_architecture"]      # deprecated alias


def test_bed_order_is_the_wicker_order():
    rec = engine.analyze(examples.load("M80343"))["records"][0]
    assert "|order=LINE|" in gff3.to_bed(rec, "x").splitlines()[0]


def test_ltr_and_tir_rows_export_tsd_state():
    for kind, five, three in (("LTR (long terminal repeat)", [0, 100], [900, 1000]),
                              ("TIR (terminal inverted repeat)", [0, 30], [970, 1000])):
        rec = {"classification": {"te_class": "repeat/structural-only", "completeness": {"tier": "structural-only"}},
               "composition": {"length": 1000},
               "structural": [{"type": kind, "five_prime": five, "three_prime": three,
                               "tsd_state": "not assessable (no flanking sequence in the record)"}],
               "domains": []}
        rows = [r for r in _rows(gff3.to_gff3(rec, seqid="x")) if _attr(r, "ID") != "TE_1"]
        assert rows and all(_attr(r, "tsd_state") == "not assessable (no flanking sequence in the record)"
                            for r in rows), kind


def test_sub_confident_pbs_is_not_exported_as_a_primer_binding_site():
    def pbs(confident):
        return {"classification": {"te_class": "LTR/Gypsy", "completeness": {"tier": "intact", "present": ["RT"]}},
                "composition": {"length": 5000},
                "structural": [{"type": "PBS (primer-binding site)", "pos": [100, 118], "confident": confident,
                                "priming_trna": "tRNA-Lys3" if confident else "undetermined"}],
                "domains": []}
    weak = [r for r in _rows(gff3.to_gff3(pbs(False), seqid="x")) if _attr(r, "ID") != "TE_1"]
    assert [r[2] for r in weak] == ["sequence_feature"] and _attr(weak[0], "Name") == "PBS_candidate"
    assert _attr(weak[0], "Ontology_term") == gff3.SO["sequence_feature"]
    strong = [r for r in _rows(gff3.to_gff3(pbs(True), seqid="x")) if _attr(r, "ID") != "TE_1"]
    assert [r[2] for r in strong] == ["primer_binding_site"]


def test_pfam_reference_matches_the_bundled_panel():
    n = sum(1 for l in open(os.path.join(_BE, "data", "te_domains.hmm"), encoding="utf-8") if l.startswith("NAME "))
    name = refs.REFS["Pfam"]["name"]
    assert f"{n} accession versions" in name
    assert "retrieved 2026-07-17/23/28" in name
    assert "Deliverables" not in (refs.__doc__ or "")


def test_genomepcr_intended_only_verdict_both_branches():
    one = genomepcr.parse_ispcr(">NC_1:100+300 pair 201bp AAA CCC")
    one[0]["on_target"] = True
    s = genomepcr.summarize(one, has_locus=True)
    assert s["verdict"] == "1 on-target product; no additional priming site detected in this assembly under current settings"
    s = genomepcr.summarize(one + genomepcr.parse_ispcr(">NC_1:500+700 fwdonly 201bp AAA AAA"), has_locus=True)
    assert s["tier"] == "intended-only" and s["verdict"].startswith("1 on-target product; no additional forward+reverse")
    assert s["verdict"].endswith("1 single-primer (F+F/R+R) product(s) also predicted")
