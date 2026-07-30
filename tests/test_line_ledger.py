"""Regression tests for the LINE completeness ledger and the non-LTR call (2026-07-27).

Before this: the non-LTR branch of _completeness hardcoded expected = present = ["RT"], missing = [],
so the ledger could not be contradicted by its own evidence and a full-length L1 returned a verdict
byte-identical to a dead 5'-truncated fragment. The call itself also read "RT without integrase and
without LTRs -> LINE" as though the absence were positive evidence, while DIRS-group (tyrosine
recombinase) and Penelope-like (GIY-YIG) elements share exactly that absence and neither enzyme is
in the bundled Pfam panel.
"""
from teagle_core import classify


def _rt_domain():
    return {"domain": "RT", "hmm": "RVT_1", "nt": [100, 1200], "strand": "+", "score": 210.0, "class": ""}


_POLYA = [{"type": "poly-A tail", "start": 5900, "end": 5990}]


def test_line_ledger_distinguishes_full_length_from_fragment():
    full = classify.classify(_POLYA, [_rt_domain()])
    frag = classify.classify([], [_rt_domain()])
    assert full["completeness"] != frag["completeness"]


def test_line_ledger_records_a_missing_tail():
    frag = classify.classify([], [_rt_domain()])["completeness"]
    assert frag["missing"], "a LINE with no 3' tail must record something missing"
    assert any("tail" in m for m in frag["missing"])
    assert "RT" in frag["present"]


def _orf1():
    return {"domain": "ORF1", "hmm": "Transposase_22", "nt": [10, 90], "strand": "+", "score": 200.0, "class": ""}


def _en(nt=(95, 140)):
    return {"domain": "EN", "hmm": "Exo_endo_phos", "nt": list(nt), "strand": "+", "score": 70.0, "class": ""}


def test_line_ledger_credits_a_recovered_tail():
    led = classify.classify(_POLYA, [_rt_domain()])["completeness"]
    assert any("tail" in p for p in led["present"])
    # RT + tail alone is NOT complete: ORF1 and the ORF2 endonuclease are modelled and were not found
    assert set(led["missing"]) == {"ORF1", "EN"}


def test_full_line_architecture_reaches_intact():
    led = classify.classify(_POLYA, [_orf1(), _en(), _rt_domain()])["completeness"]
    assert led["missing"] == []
    assert "intact" in led["tier"]


def test_endonuclease_is_only_credited_upstream_of_the_rt():
    """Pfam's endonuclease family also covers host DNase I, so EN counts only in the ORF2p EN-RT order."""
    downstream = classify.classify(_POLYA, [_orf1(), _en(nt=(2000, 2100)), _rt_domain()])["completeness"]
    assert "EN" in downstream["missing"]
    upstream = classify.classify(_POLYA, [_orf1(), _en(), _rt_domain()])["completeness"]
    assert "EN" in upstream["present"]


def test_line_tier_never_claims_intact():
    """ORF1 and the endonuclease are not in the bundled panel, so a LINE cannot be shown autonomous."""
    for structural_ev in ([], _POLYA):
        tier = classify.classify(structural_ev, [_rt_domain()])["completeness"]["tier"]
        assert "intact" not in tier.lower()
        assert "autonomous" not in tier.lower()


def test_line_without_a_tail_is_only_a_candidate():
    """The tail is the sole positive evidence; without it the call rests on an absence."""
    assert classify.classify([], [_rt_domain()])["confidence"] == "Candidate"
    assert classify.classify(_POLYA, [_rt_domain()])["confidence"] == "Moderate"


def test_line_call_names_the_alternative_it_cannot_exclude():
    """DIRS is now tested directly (a YR branch ahead of the RT branch), so it is no longer a hedge.
    Penelope-like elements stay unexcluded because their GIY-YIG endonuclease is not in the panel."""
    ev = " ".join(classify.classify(_POLYA, [_rt_domain()])["evidence"]).lower()
    assert "penelope" in ev
    assert "not excluded" in ev


def test_yr_with_rt_classifies_as_dirs_not_line():
    """A tyrosine recombinase in place of the DDE integrase is a DIRS-group element. Before the branch
    reorder these fell through to 'LINE (non-LTR)' with an affirmative evidence string."""
    yr = {"domain": "YR", "hmm": "Phage_integrase", "nt": [2000, 2400], "strand": "+", "score": 30.0, "class": ""}
    cl = classify.classify([], [_rt_domain(), yr])
    assert cl["te_class"] == "DIRS"
    assert "dirs" in cl["superfamily"].lower()
    ev = " ".join(cl["evidence"]).lower()
    assert "host" in ev and "recombinase" in ev          # the promiscuity caveat must travel with the call


# ---------------- the path that shipped a NameError (fix loop round 1, CRITICAL) ----------------
def _tir(a=100, b=130, c=900, d=930):
    return {"type": "TIR (terminal inverted repeat)", "tir_len": b - a, "identity": 95.0,
            "five_prime": [a, b], "three_prime": [c, d], "element_span": [a, d]}


def _tsd(length=8):
    return {"type": "TSD (target-site duplication)", "length": length, "motif": "A" * length,
            "upstream": [100 - length, 100], "downstream": [930, 930 + length]}


def _tpase(cls="dna:hAT"):
    return {"domain": "TPase", "hmm": "Dimer_Tnp_hAT", "nt": [200, 800], "strand": "+",
            "score": 120.0, "class": cls}


def test_dna_transposon_with_tir_and_tsd_does_not_crash():
    """classify.py called structural_mod.tsd_congruence() while the module was never imported, so the
    TEXTBOOK Class II case — a transposase enclosed by TIRs with a flanking target-site duplication —
    raised NameError. No existing test reached it: the benchmark's Tc1 and Ac are bare elements with no
    flanking sequence, so has_tsd was always False."""
    cl = classify.classify([_tir(), _tsd()], [_tpase()])
    assert cl["te_class"] == "DNA/hAT"


def test_tsd_congruence_reaches_the_evidence_string():
    congruent = " ".join(classify.classify([_tir(), _tsd(8)], [_tpase()])["evidence"])
    assert "target-site duplication" in congruent
    assert "congruent" in congruent, "a congruent TSD length should be stated"
    incongruent = " ".join(classify.classify([_tir(), _tsd(3)], [_tpase()])["evidence"])
    assert "3 bp where" in incongruent or "coincidental" in incongruent, \
        "an incongruent TSD must be reported, not silently accepted as confirming the termini"


def test_tsd_refine_prefers_superfamily_length_over_coincidental_flank():
    """detect_all picks the LONGEST exact flanking repeat before the superfamily is known; once classify
    resolves a Tc1/Mariner call it must re-detect with the diagnostic 2 bp TA preferred, so a coincidental
    longer flank does not flip a genuinely complete element to 'incongruent'. (fix loop round 5, MODERATE)"""
    # element [100, 930] flanked by an exact 6-mer TATATA on both sides: longest-first reports 6 bp, but the
    # immediate 2 bp flank is also an exact TA repeat, so the superfamily-aware re-run must prefer the 2 bp TA.
    chars = list("C" * 940)
    chars[94:100] = list("TATATA"); chars[930:936] = list("TATATA")
    seq = "".join(chars)
    tsd = {"type": "TSD (target-site duplication)", "length": 6, "motif": "TATATA",
           "upstream": [94, 100], "downstream": [930, 936]}
    ev = " ".join(classify.classify([_tir(), tsd], [_tpase(cls="dna:Tc1-Mariner")], seq=seq)["evidence"])
    assert tsd["length"] == 2 and tsd["motif"] == "TA"        # corrected in place -> the GFF3 export carries 2 bp
    assert "congruent" in ev and "incongruent" not in ev      # the corrected length now corroborates the termini


def test_tsd_refine_is_a_noop_without_the_sequence():
    # backward compatible: with no seq (every non-engine caller) the ungated longest-first length stands
    tsd = {"type": "TSD (target-site duplication)", "length": 6, "motif": "GGGGTA",
           "upstream": [94, 100], "downstream": [930, 936]}
    ev = " ".join(classify.classify([_tir(), tsd], [_tpase(cls="dna:Tc1-Mariner")])["evidence"])
    assert tsd["length"] == 6                                 # unchanged
    assert "incongruent" in ev or "6 bp where" in ev


def test_ltr_congruent_tsd_evidence_states_the_observed_length():
    """The LTR congruent-TSD sentence once read 'a target-site duplication of the length Copia duplicates
    (...) flanks the element' — grammatically broken and missing the observed length. (fix loop round 5, LOW)"""
    struct = [{"type": "LTR (terminal direct repeat)", "element_span": [0, 9000], "identity": 98.0},
              {"type": "TSD (target-site duplication)", "length": 5, "motif": "ACGTA",
               "upstream": [95, 100], "downstream": [9000, 9005]}]
    doms = [{"domain": "INT", "nt": [100, 300], "strand": "+", "score": 40.0, "class": "retro"},
            {"domain": "RT", "nt": [500, 1400], "strand": "+", "score": 90.0, "class": "retro"}]
    ev = " ".join(classify.classify(struct, doms)["evidence"])
    assert "5 bp target-site duplication" in ev               # the observed length is stated, not omitted
    assert "congruent" in ev and "incongruent" not in ev
