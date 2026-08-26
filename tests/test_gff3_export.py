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
    ("DNA/Tc1-Mariner", ["TPase", "TIR (5′)", "TIR (3′)"], [], "terminal_inverted_repeat_element"),
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
    """An isolated ENV hit plus terminal direct repeats does not establish LTR origin, so neither the
    classification nor GFF3 column 3 may assert an LTR retrotransposon."""
    from teagle_core import classify
    d = classify.classify(
        structural=[{"type": "LTR (terminal direct repeat)", "five_prime": [0, 200], "three_prime": [1000, 1200]}],
        domains=[{"domain": "ENV", "class": "", "score": 100, "nt": [300, 600], "strand": "+"}])
    assert d["te_class"] == "repeat/structural-only"
    assert d["outcome"] == "withheld"
    assert "RT" not in d["completeness"]["present"]
    assert gff3.so_term_for(d) == "repeat_region"
    attrs = _rows(gff3.to_gff3({"classification": d, "completeness": d["completeness"],
                                "structural": [{"type": "LTR (terminal direct repeat)",
                                                "five_prime": [0, 200], "three_prime": [1000, 1200]}],
                                "domains": [], "composition": {"length": 1200}}))[0][8]
    assert "classification_outcome=withheld" in attrs
    assert "classification_order=withheld" in attrs


def test_trusted_ltr_fragment_gff3_preserves_hedge_while_using_neutral_so_term():
    from teagle_core import classify
    structural = [{"type": "LTR (terminal direct repeat)", "five_prime": [0, 200],
                   "three_prime": [1000, 1200]}]
    domains = [{"domain": "GAG", "hmm": "Gag_p24", "class": "", "score": 100,
                "nt": [300, 600], "strand": "+"},
               {"domain": "PR", "hmm": "RVP", "class": "", "score": 90,
                "nt": [650, 800], "strand": "+"}]
    cl = classify.classify(structural, domains)
    parent = _rows(gff3.to_gff3({"classification": cl, "completeness": cl["completeness"],
                                 "structural": structural, "domains": domains,
                                 "composition": {"length": 1200}}))[0]
    assert parent[2] == "repeat_region"
    assert "wicker_code=LTR/partial" in parent[8]
    assert "classification_outcome=complete" in parent[8]
    assert "classification_order=unassigned" in parent[8]
    assert "domain_evidence_status=tested_with_hits" in parent[8]
    assert "classification_explanation=" in parent[8]
    assert "RT/pol not detected" in parent[8]


def test_failed_domain_gff3_and_bed_preserve_withheld_unavailable_state():
    from teagle_core import classify
    structural = [{"type": "LTR (terminal direct repeat)", "five_prime": [0, 100],
                   "three_prime": [300, 400]}]
    cl = classify.classify(structural, [], domains_ok=False)
    rec = {"classification": cl, "completeness": cl["completeness"],
           "structural": structural, "domains": [], "composition": {"length": 400}}
    parent = _rows(gff3.to_gff3(rec))[0]
    assert parent[2] == "repeat_region"
    assert "classification_outcome=withheld" in parent[8]
    assert "classification_order=withheld" in parent[8]
    assert "domain_evidence_status=unavailable_failed" in parent[8]
    bed_parent = gff3.to_bed(rec).splitlines()[0].split("\t")
    assert "|outcome=withheld|order=withheld|domains=unavailable_failed" in bed_parent[3]


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
    safe = gff3.public_record_id("M11240")
    assert any(l.startswith(f"##sequence-region {safe} 1 ") for l in text.splitlines())


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
    assert all("=" in item for item in attrs.split(";"))


def test_embedded_fasta_is_emitted_when_a_sequence_is_supplied():
    text = gff3.to_gff3(_rec("X01005"), seqid="X01005", sequence="ACGT" * 30)
    assert "##FASTA" in text
    assert f">{gff3.public_record_id('X01005')}" in text
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


def test_esc_percent_encoded_first():
    """A literal '%' in a value (e.g. a PBS note's '61.1%') must become '%25', and the '%' introduced for
    other reserved chars ('%3B') must not be double-encoded. Regression: _esc omitted '%'."""
    assert gff3._esc("61.1%") == "61.1%25"
    assert gff3._esc("a;b") == "a%3Bb"                       # ';' still encodes to a single, un-doubled %3B


def test_ppt_purine_frac_and_tsd_congruence_reach_the_export():
    """A PPT row must carry its purine_frac quality metric; a TSD row must carry the superfamily-length
    congruence verdict classify() computed (the hedge shown on screen). Both were dropped from the file."""
    from teagle_core import classify
    st = [{"type": "TIR (terminal inverted repeat)", "five_prime": [0, 20], "three_prime": [880, 900]},
          {"type": "TSD (target-site duplication)", "upstream": [0, 3], "downstream": [900, 903], "length": 3},
          {"type": "PPT (polypurine tract)", "pos": [850, 862], "purine_frac": 0.92}]
    dm = [{"domain": "TPase", "class": "dna:hAT", "nt": [100, 800], "strand": "+", "score": 200.0}]
    cl = classify.classify(st, dm)
    rec = {"classification": cl, "composition": {"length": 1000}, "structural": st, "domains": dm}
    g = gff3.to_gff3(rec, seqid="x")
    assert "purine_frac=0.92" in g
    assert "tsd_congruence=incongruent" in g                 # hAT expects 8 bp; a 3 bp TSD is incongruent


def test_pbs_confidence_hedge_survives_into_the_export():
    """A diverged PBS is reported on screen with priming_trna='undetermined' + a note; the GFF3 must carry
    that hedge, not a bare primer_binding_site a browser would read as a confident call. (An export that
    silently drops a hedge is a defect — this locks the invariant so it can't regress.)"""
    rec = {"classification": {"te_class": "LTR/Gypsy", "superfamily": "Gypsy", "confidence": "High",
                              "completeness": {"tier": "intact", "present": ["RT", "INT"], "missing": [], "scope": "x"}},
           "composition": {"length": 5000},
           "structural": [{"type": "PBS (primer-binding site)", "pos": [100, 118],
                           "priming_trna": "undetermined", "best_match": "tRNA-Lys3", "identity": 61.1,
                           "confident": False, "note": "priming tRNA undetermined — often diverged"}],
           "domains": []}
    g = gff3.to_gff3(rec, seqid="x")
    row = next((l for l in g.splitlines() if "primer_binding_site" in l), "")
    assert "priming_trna=undetermined" in row
    assert "confident=false" in row
    assert "identity_pct=61.1" in row
    assert "note=" in row                                    # the explanatory hedge travels too


def test_poly_a_tail_survives_into_the_export():
    """M80343 (LINE-1) has a 3' poly-A tail — a dedicated on-screen genome-viewer track. It was dropped
    from GFF3/BED because its type ('poly-A tail') had no _STRUCT_TERM entry. It now exports as the generic
    sequence_feature (SO:0000110): no genomic-specific SO term exists, and polyA_sequence (SO:0000610) is an
    mRNA feature, so asserting it here would misdescribe the genomic tract."""
    rows = _rows(gff3.to_gff3(_rec("M80343"), seqid="M80343"))
    tail = [r for r in rows if r[2] == "sequence_feature"]
    assert tail, "the poly-A tail must survive into the GFF3 (not be silently dropped)"
    assert all("SO:0000110" in r[8] for r in tail)
    assert any("Name=poly-A" in r[8] or "Name=poly-T" in r[8] for r in tail)


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
    rec = {"classification": {"te_class": "DNA/Tc1-Mariner",
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
    assert sum(name.startswith("TSD|") for name in bed_names) == 2
    assert any(name.startswith("PBS|") for name in bed_names)


# ---------------- fetched gene-model and splice coordinate exports ----------------
def _gene_model():
    return {"exons": [{"start": 12, "end": 20, "strand": "+", "note": "2"},
                       {"start": 0, "end": 6, "strand": "+", "note": "1"}],
            "introns": [{"start": 6, "end": 12, "strand": "+", "donor": "GT",
                         "acceptor": "AG", "canonical": True}],
            "cds": [{"start": 2, "end": 6, "strand": "+", "phase": 0, "cds_group": "cds_1"},
                    {"start": 12, "end": 18, "strand": "+", "phase": 2, "cds_group": "cds_1"}]}


def test_gene_model_gff3_and_bed_convert_0_based_half_open_exactly_and_sort():
    model = _gene_model()
    g_rows = [row for row in _rows(gff3.gene_model_to_gff3(
        model, "chr1", sequence="ACGT" * 5)) if len(row) == 9]
    b_rows = [line.split("\t") for line in gff3.gene_model_to_bed(
        model, "chr1", sequence_length=20).splitlines()]
    assert [(int(r[3]), int(r[4])) for r in g_rows] == [(1, 6), (3, 6), (7, 12), (13, 18), (13, 20)]
    assert [(int(r[1]), int(r[2])) for r in b_rows] == [(0, 6), (2, 6), (6, 12), (12, 18), (12, 20)]
    assert all(int(b[1]) == int(g[3]) - 1 and int(b[2]) == int(g[4])
               for g, b in zip(g_rows, b_rows))
    assert "donor=GT" in g_rows[2][8] and "acceptor=AG" in g_rows[2][8]
    cds = [row for row in g_rows if row[2] == "CDS"]
    assert [row[7] for row in cds] == ["0", "2"]
    assert all("cds_group=cds_1" in row[8] for row in cds)


@pytest.mark.parametrize("phase", [None, -1, 3, True, "0"])
def test_gene_model_gff3_fails_closed_when_cds_phase_is_unavailable_or_invalid(phase):
    model = _gene_model()
    model["cds"][0]["phase"] = phase
    with pytest.raises(ValueError, match="phase"):
        gff3.gene_model_to_gff3(model, "chr1", sequence="ACGT" * 5)
    assert "\tCDS_" in gff3.gene_model_to_bed(model, "chr1", sequence_length=20)
    assert ">chr1.cds" in gff3.gene_model_to_fasta(
        model, "ACGT" * 5, "chr1", feature_types=("CDS",), identity_trust="validated_fetch")


@pytest.mark.parametrize("mutation", [
    lambda m: m["exons"][0].update(start=-1),
    lambda m: m["exons"][0].update(end=12),
    lambda m: m["exons"][0].update(start=20, end=20),
    lambda m: m["exons"][0].update(strand="?"),
    lambda m: m["exons"][0].update(start=True),
])
def test_gene_model_coordinate_and_strand_mutations_fail(mutation):
    import copy
    model = copy.deepcopy(_gene_model())
    mutation(model)
    with pytest.raises(ValueError):
        gff3.gene_model_to_gff3(model, "chr1", sequence="ACGT" * 5)
    with pytest.raises(ValueError):
        gff3.gene_model_to_bed(model, "chr1", sequence_length=20)


def test_gene_model_rows_are_deterministic_under_input_reordering():
    import copy
    a = _gene_model()
    b = copy.deepcopy(a)
    b["exons"].reverse()
    b["cds"].reverse()
    assert gff3.gene_model_to_gff3(a, "chr1", sequence="ACGT" * 5) == \
        gff3.gene_model_to_gff3(b, "chr1", sequence="ACGT" * 5)
    assert gff3.gene_model_to_bed(a, "chr1", sequence_length=20) == \
        gff3.gene_model_to_bed(b, "chr1", sequence_length=20)


def test_path_like_record_id_is_rewritten_without_losing_identity():
    record_id = "../../sample A;unit=1\\x"
    safe = gff3.sanitize_seqid(record_id)
    stem = gff3.safe_filename_stem(record_id)
    assert safe.startswith("record_") and "/" not in safe and "\\" not in safe
    assert "/" not in stem and "\\" not in stem and stem not in {".", ".."}
    text = gff3.gene_model_to_gff3(_gene_model(), record_id, sequence="ACGT" * 5)
    assert f"##sequence-region {safe} 1 20" in text
    digest = gff3._identity_sha256(record_id)
    assert f"original_record_sha256={digest}" in text
    assert record_id not in text
    bed = gff3.gene_model_to_bed(_gene_model(), record_id, sequence_length=20)
    assert f"record_sha256={digest}" in bed and record_id not in bed
    assert gff3.safe_filename_stem("CON") != "CON"


def test_te_annotation_writers_also_sanitize_path_like_record_ids():
    record_id = "../unsafe sample"
    safe = gff3.sanitize_seqid(record_id)
    rec = _rec("M11240")
    text = gff3.to_gff3(rec, record_id, sequence=rec["seq"])
    rows = [row for row in _rows(text) if len(row) == 9]
    assert all(row[0] == safe for row in rows)
    digest = gff3._identity_sha256(record_id)
    assert f"original_record_sha256={digest}" in rows[0][8]
    assert record_id not in text
    assert f">{safe}" in text and f">{record_id}" not in text
    bed = gff3.to_bed(rec, record_id)
    assert all(line.split("\t")[0] == safe for line in bed.splitlines())
    assert f"record_sha256={digest}" in bed and record_id not in bed


@pytest.mark.parametrize("record_id", [
    r"C:\Users\person\private\sample.fa",
    "C:private_sample.fa",
    "/home/person/private/sample.fa",
    "../../private/sample.fa",
])
def test_filesystem_like_ids_never_reappear_in_any_public_coordinate_export(record_id):
    digest = gff3._identity_sha256(record_id)
    assert gff3.safe_filename_stem(record_id) == f"record_{digest[:12]}"
    assert gff3.safe_filename_stem(record_id, default="../unsafe") == f"record_{digest[:12]}"
    outputs = [
        gff3.gene_model_to_gff3(_gene_model(), record_id, sequence="ACGT" * 5),
        gff3.gene_model_to_bed(_gene_model(), record_id, sequence_length=20),
        gff3.gene_model_to_fasta(_gene_model(), "ACGT" * 5, record_id, feature_types=("exon",)),
    ]
    for output in outputs:
        assert record_id not in output
        assert digest in output


def test_gene_model_fasta_uses_exact_bound_sequence_and_strand():
    model = {"exons": [{"start": 0, "end": 4, "strand": "+"},
                       {"start": 4, "end": 8, "strand": "-"}],
             "cds": [{"start": 2, "end": 6, "strand": "+", "phase": 0}]}
    text = gff3.gene_model_to_fasta(model, "AAAACCGT", "record A", feature_types=("exon",))
    records = [block for block in text.split(">") if block]
    assert len(records) == 2
    assert records[0].splitlines()[1] == "AAAA"
    assert records[1].splitlines()[1] == "ACGG"  # reverse complement of CCGT
    assert f"record_sha256={gff3._identity_sha256('record A')}" in text
    assert "record A" not in text and "coordinate_system=0-based-half-open" in text
    mutated = {"exons": [{"start": 4, "end": 9, "strand": "+"}]}
    with pytest.raises(ValueError, match="exceeds"):
        gff3.gene_model_to_fasta(mutated, "AAAACCGT", "record")
    with pytest.raises(ValueError, match="IUPAC"):
        gff3.gene_model_to_fasta(model, "AAAAXCGT", "record", feature_types=("exon",))
    with pytest.raises(ValueError, match="not FASTA"):
        gff3.gene_model_to_fasta(model, ">record\nAAAACCGT", "record", feature_types=("exon",))


def test_gene_model_export_has_no_nan_or_machine_path_leakage():
    text = gff3.gene_model_to_gff3(_gene_model(), "chr1", sequence="ACGT" * 5,
                                    source="NCBI annotation", model_kind="fetched_gene_model")
    assert "NaN" not in text and "C:\\" not in text and "/Users/" not in text
    assert "model_kind=fetched_gene_model" in text


@pytest.mark.parametrize("record_id", [
    r"C:\Users\person\private\sample.fa",
    "/home/person/private/sample.fa",
    "../../private/sample.fa",
    r"\\server\share\sample.fa",
])
def test_public_record_identity_rewrites_filesystem_values_with_full_digest(record_id):
    identity = gff3.public_record_identity(record_id)
    digest = gff3._identity_sha256(record_id)
    assert identity == {"id": f"record_{digest[:12]}", "idSha256": digest}
    assert gff3.public_record_id(record_id) == identity["id"]


@pytest.mark.parametrize("record_id", [
    "M11240", "M11240.1", "AM040263", "NC_000013.11", "GCF_009914755.1",
    "chr1", "chr1:100-200",
])
def test_public_record_identity_preserves_recognized_ids_only_in_validated_fetch_context(record_id):
    assert gff3.public_record_identity(record_id, trust="validated_fetch") == {"id": record_id}
    assert "idSha256" in gff3.public_record_identity(record_id)


@pytest.mark.parametrize("record_id", ["input_sequence", "record_deadbeefdead", "P1"])
def test_public_record_identity_preserves_generated_ids_only_in_generated_context(record_id):
    assert gff3.public_record_identity(record_id, trust="generated") == {"id": record_id}
    assert "idSha256" in gff3.public_record_identity(record_id)


@pytest.mark.parametrize("record_id", [
    "patient_042", "donor-A1", "sample_tuna_birgun", "descriptive specimen header",
    " specimen ", "specimen\tprivate", r"C:\private\sample.fa", "../../sample.fa",
])
def test_public_record_identity_rewrites_sample_like_and_free_text_headers(record_id):
    digest = gff3._identity_sha256(record_id)
    expected = {"id": f"record_{digest[:12]}", "idSha256": digest}
    assert gff3.public_record_identity(record_id) == expected
    assert gff3.sanitize_seqid(record_id) == expected["id"]
    assert gff3.safe_filename_stem(record_id) == expected["id"]
