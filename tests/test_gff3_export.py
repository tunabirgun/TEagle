"""GFF3 / BED export (roadmap 4.1).

Column 3 of a GFF3 file is an assertion every downstream tool trusts without re-deriving it. The whole
point of gating it on the completeness tier is that TEagle must not tell a genome browser an element is
an LTR_retrotransposon when its own report says structural-only — the browser has no way to see the
hedge, so the hedge has to be in the term itself.
"""
import os
import sys

import pytest

_BE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "app", "backend")
if _BE not in sys.path:
    sys.path.insert(0, _BE)
from teagle_core import examples, gff3                          # noqa: E402
import engine                                                    # noqa: E402


def _rec(acc):
    return engine.analyze(examples.load(acc))["records"][0]


def _rows(text):
    return [l.split("\t") for l in text.splitlines() if l and not l.startswith("#")]


# ---------------- the completeness gate ----------------
def test_structural_only_never_claims_a_specific_subclass():
    cl = {"te_class": "LTR/structural-only",
          "completeness": {"tier": "structural-only (terminal repeats, no coding domain detected)"}}
    assert gff3.so_term_for(cl) == "repeat_region"


def test_unclassified_exports_as_repeat_region():
    assert gff3.so_term_for({"te_class": "none"}) == "repeat_region"
    assert gff3.so_term_for({}) == "repeat_region"
    assert gff3.so_term_for(None) == "repeat_region"


@pytest.mark.parametrize("te_class,present,missing,expected", [
    # a retrotransposon subclass may claim its term only when RT is in the ledger;
    # a DNA/ TIR-element term only when neither terminal inverted repeat is missing.
    ("LTR/Copia", ["GAG", "PR", "RT", "RNaseH", "INT"], [], "LTR_retrotransposon"),
    ("LTR/Gypsy", ["RT", "INT"], ["GAG"], "LTR_retrotransposon"),
    ("LINE", ["ORF1", "EN", "RT"], [], "LINE_element"),
    ("DIRS", ["RT"], [], "YR_retrotransposon"),
    ("DNA/Tc1", ["TPase", "TIR (5′)", "TIR (3′)"], [], "terminal_inverted_repeat_element"),
])
def test_classes_with_coding_evidence_get_their_term(te_class, present, missing, expected):
    cl = {"te_class": te_class,
          "completeness": {"tier": "intact / autonomous-consistent", "present": present, "missing": missing}}
    assert gff3.so_term_for(cl) == expected


def test_dna_transposon_without_recovered_tirs_degrades_to_the_generic_term():
    """A transposase with no terminal inverted repeat recovered must NOT export the SO term whose defining
    feature is exactly those repeats; it degrades to DNA_transposon (DNA-mediated only), what the evidence
    supports. Regression for the CRITICAL where a 5'/3'-truncated copy asserted its own absent ends."""
    from teagle_core import classify
    d = classify.classify(structural=[],
                          domains=[{"domain": "TPase", "class": "hAT", "score": 200, "nt": [100, 900], "strand": "+"}])
    assert d["te_class"].startswith("DNA/")
    assert any("TIR" in m for m in d["completeness"]["missing"])
    assert gff3.so_term_for(d) == "DNA_transposon"


def test_ltr_fragment_without_rt_degrades_to_repeat_region():
    """A relic that kept a coding domain (env/capsid) but lost pol has no reverse-transcriptase evidence,
    so column 3 cannot assert LTR_retrotransposon. Regression for the second export-overclaim CRITICAL."""
    from teagle_core import classify
    d = classify.classify(
        structural=[{"type": "LTR (terminal direct repeat)", "five_prime": [0, 200], "three_prime": [1000, 1200]}],
        domains=[{"domain": "ENV", "class": "", "score": 100, "nt": [300, 600], "strand": "+"}])
    assert d["te_class"] == "LTR/partial"
    assert "RT" not in d["completeness"]["present"]
    assert gff3.so_term_for(d) == "repeat_region"


def test_every_term_used_has_a_verified_so_accession():
    """Terms and accessions were checked against SO-Ontologies so-simple.obo. A term without an
    accession would export an Ontology_term a consumer cannot resolve."""
    for name, acc in gff3.SO.items():
        assert acc.startswith("SO:") and acc[3:].isdigit(), (name, acc)
    for term in list(gff3._CLASS_TERM.values()) + list(gff3._STRUCT_TERM.values()):
        assert term in gff3.SO, f"{term} is emitted but has no accession"


# ---------------- document structure ----------------
def test_gff3_has_the_required_directives():
    text = gff3.to_gff3(_rec("M11240"), seqid="M11240")
    assert text.splitlines()[0] == "##gff-version 3"
    assert any(l.startswith("##sequence-region M11240 1 ") for l in text.splitlines())


def test_every_feature_row_has_nine_columns_and_sane_coordinates():
    text = gff3.to_gff3(_rec("M80343"), seqid="M80343")
    rows = _rows(text)
    assert rows
    for r in rows:
        assert len(r) == 9, r
        start, end = int(r[3]), int(r[4])
        assert start >= 1, "GFF3 is 1-based; a start below 1 is invalid"
        assert end >= start, r


def test_sub_features_reference_their_parent():
    rows = _rows(gff3.to_gff3(_rec("M11240"), seqid="M11240"))
    parent_id = [a for a in rows[0][8].split(";") if a.startswith("ID=")][0][3:]
    children = [r for r in rows[1:] if "Parent=" in r[8]]
    assert children, "an element with LTRs and domains must emit sub-features"
    for c in children:
        assert f"Parent={parent_id}" in c[8]


def test_wicker_code_travels_in_a_lowercase_attribute():
    """GFF3 reserves capitalised tags for its predefined set; a custom TEID= tag is non-compliant."""
    attrs = _rows(gff3.to_gff3(_rec("M11240"), seqid="M11240"))[0][8]
    assert "wicker_code=" in attrs
    assert "TEID=" not in attrs


def test_reserved_characters_are_percent_encoded():
    attrs = _rows(gff3.to_gff3(_rec("M11240"), seqid="M11240"))[0][8]
    # the scope sentence contains commas, which are the GFF3 attribute list separator
    assert "%2C" in attrs
    body = attrs.split("domains_tested_scope=")[-1]
    assert ";" not in body or body.count(";") == 0


def test_embedded_fasta_is_emitted_when_a_sequence_is_supplied():
    text = gff3.to_gff3(_rec("X01005"), seqid="X01005", sequence="ACGT" * 30)
    assert "##FASTA" in text
    assert ">X01005" in text
    assert text.rstrip().endswith("ACGT") or "ACGT" in text.split("##FASTA")[1]


def test_no_fasta_section_when_no_sequence_given():
    assert "##FASTA" not in gff3.to_gff3(_rec("X01005"), seqid="X01005")


def test_lower_bound_flag_survives_into_the_export():
    """A terminal-repeat length limited by the record is hedged in the UI; the file must carry the hedge
    too, or the browser shows a measurement TEagle never claimed."""
    rec = _rec("M11240")
    for ev in rec.get("structural") or []:
        ev["length_is_lower_bound"] = True
    text = gff3.to_gff3(rec, seqid="M11240")
    assert "length_is_lower_bound=true" in text


# ---------------- BED ----------------
def test_bed_is_zero_based_half_open():
    rec = _rec("M11240")
    gff_rows = _rows(gff3.to_gff3(rec, seqid="M11240"))
    bed_rows = [l.split("\t") for l in gff3.to_bed(rec, "M11240").splitlines() if l]
    assert len(bed_rows) == len(gff_rows)
    for g, b in zip(gff_rows, bed_rows):
        assert int(b[1]) == int(g[3]) - 1, "BED start must be the GFF3 start minus one"
        assert int(b[2]) == int(g[4])


def test_bed_rows_have_six_columns():
    for line in gff3.to_bed(_rec("M80343"), "M80343").splitlines():
        assert len(line.split("\t")) == 6, line


def test_bed_score_is_a_valid_0_to_1000_integer():
    """BED's score column is an integer 0-1000 that drives shading; a raw HMMER bit score is neither
    integer-bounded nor on that scale, and some browsers reject the row outright."""
    for line in gff3.to_bed(_rec("M11240"), "M11240").splitlines():
        if not line:
            continue
        score = line.split("\t")[4]
        assert score.isdigit(), f"BED score must be a non-negative integer, got {score!r}"
        assert 0 <= int(score) <= 1000, f"BED score out of range: {score}"


# ---------------- structural sub-features that were being dropped ----------------
def test_ppt_evidence_survives_into_the_export():
    """M11240 (copia) detects a polypurine tract. build_features once read only five_prime/three_prime/
    start, so PBS/PPT (keyed 'pos') and TSD ('upstream'/'downstream') were silently dropped — a scientific
    export claiming less evidence than the analysis produced."""
    rows = _rows(gff3.to_gff3(_rec("M11240"), seqid="M11240"))
    assert any(r[2] == "RR_tract" for r in rows), "PPT evidence must survive into the GFF3"


def test_multi_record_export_embeds_the_selected_records_own_sequence():
    """The UI once passed analyzed_clean (every record of a multi-FASTA paste concatenated) as the export
    sequence, so a non-first record's GFF3 embedded the wrong bases under its own coordinates. Each record
    now carries its own sequence; the embedded FASTA must be that record's, not the whole input's."""
    res = engine.analyze(">recA\n" + "ACGTACGT" * 20 + "\n>recB\n" + "GGGGCCCC" * 15)
    recs = res["records"]
    assert len(recs) == 2
    rb = recs[1]
    assert rb.get("seq") and len(rb["seq"]) == rb["composition"]["length"]
    text = gff3.to_gff3(rb, seqid="recB", sequence=rb.get("seq"))
    emb = "".join(l for l in text.split("##FASTA", 1)[1].splitlines() if l and not l.startswith(">"))
    assert len(emb) == rb["composition"]["length"], "embedded FASTA length must equal recB's own length"
    assert set(emb) <= set("GC"), "recB is all G/C; recA's A/T must not leak into recB's embedded sequence"


def test_tsd_and_pbs_sub_features_are_not_dropped():
    rec = {"classification": {"te_class": "DNA/Tc1",
                              "completeness": {"tier": "intact / autonomous-consistent",
                                               "present": ["TPase", "TIR (5′)", "TIR (3′)"], "missing": []}},
           "composition": {"length": 2000},
           "structural": [
               {"type": "TSD (target-site duplication)", "upstream": [0, 2], "downstream": [1998, 2000]},
               {"type": "PBS (primer-binding site)", "pos": [10, 28]},
           ],
           "domains": []}
    terms = [r[2] for r in _rows(gff3.to_gff3(rec, seqid="locus"))]
    assert terms.count("target_site_duplication") == 2, "a TSD is two flanking repeats; both must be emitted"
    assert "primer_binding_site" in terms
    bed_names = [l.split("\t")[3] for l in gff3.to_bed(rec, "locus").splitlines() if l]
    assert bed_names.count("TSD") == 2 and "PBS" in bed_names
