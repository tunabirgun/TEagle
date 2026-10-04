import copy
import random
from types import SimpleNamespace

import pytest

import engine
from engine import BadRequest
from teagle_core import oligoqc, primers, sequtil
from helpers import fixture_seq

RC = sequtil.reverse_complement


def rnd(n, seed):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


def mutate(s, positions):
    comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
    return "".join(comp[b] if i in positions else b for i, b in enumerate(s))


def design(seq, **params):
    return engine.run_primers({"sequence": ">t\n" + seq, "params": {"prod_min": 150, "prod_max": 500, **params}})


@pytest.fixture(autouse=True)
def _no_wsl_fold(monkeypatch):
    monkeypatch.setattr(engine, "_prime_vienna", lambda c: None)


# ---- spec_conformance-11: the product carries the oligo over its whole footprint ----------------

def test_tailed_product_carries_the_full_oligo_at_default_mismatch_budget():
    core_f, core_r = "GACTGACTGTCAGTCAGGCT", "TTGGCCATTGGCACTGGCAT"
    tail_f, tail_r = "GGATCCGG", "AAGCTTAAA"
    seq = "C" * 50 + core_f + "A" * 200 + RC(core_r) + "C" * 50
    fwd, rev = tail_f + core_f, tail_r + core_r
    pair = [a for a in primers.in_silico_pcr(fwd, rev, seq, "t", prod_min=100, prod_max=600)
            if not a["single_primer"]]
    assert len(pair) == 1
    a = pair[0]
    assert a["seq"].startswith(fwd) and a["seq"].endswith(RC(rev))
    assert len(a["seq"]) == a["product_length"]


def test_untailed_site_with_mismatches_reports_the_primer_bases_in_the_product():
    fwd, rev = rnd(20, 1), rnd(20, 2)
    site = mutate(fwd, {2, 6})
    seq = rnd(60, 3) + site + rnd(200, 4) + RC(rev) + rnd(60, 5)
    pair = [a for a in primers.in_silico_pcr(fwd, rev, seq, "t", prod_min=100, prod_max=600)
            if not a["single_primer"]]
    assert len(pair) == 1
    assert pair[0]["seq"][:20] == fwd
    assert pair[0]["seq"][-20:] == RC(rev)
    assert pair[0]["fwd_mm"] == 2


def test_overlapping_footprints_keep_product_length_and_start_with_the_forward_oligo():
    seq = rnd(240, 6)
    fwd, rev = seq[100:120], RC(seq[110:130])
    pair = [a for a in primers.in_silico_pcr(fwd, rev, seq, "t", max_mm=0, prod_min=20, prod_max=100)
            if not a["single_primer"]]
    assert len(pair) == 1
    a = pair[0]
    assert (a["start"], a["end"]) == (100, 130)
    assert a["seq"] == seq[100:130] and len(a["seq"]) == a["product_length"]


# ---- primers2-1: IUPAC codes in the template must not make design impossible ----------------------

def test_iupac_template_designs_and_primers_avoid_the_ambiguous_bases():
    seq = fixture_seq("M11240")[:800]
    amb = seq[:300] + "RYKM" + seq[304:]
    res = design(amb)
    assert res["candidates"]
    for c in res["candidates"]:
        assert c["left_pos"][1] <= 300 or c["left_pos"][0] >= 304
        assert c["right_pos"][1] <= 300 or c["right_pos"][0] >= 304
    assert res["template_ambiguity"]["masked_to_n"] == 4
    assert "N" in res["template_ambiguity"]["note"]


def test_single_ambiguity_code_in_a_full_element_still_designs():
    seq = fixture_seq("M11240")
    res = design(seq[:4000] + "Y" + seq[4001:])
    assert res["candidates"]


def test_clean_template_reports_no_ambiguity_masking():
    res = design(fixture_seq("M11240"))
    assert res["template_ambiguity"]["masked_to_n"] == 0


def test_template_rejection_is_not_blamed_on_parameters():
    seq = fixture_seq("M11240")[:800]
    with pytest.raises(BadRequest) as e:
        design(seq[:300] + "-" + seq[301:])
    msg = str(e.value)
    assert "parameters rejected" not in msg
    assert "template" in msg.lower()


# ---- primers2-14 / persona-U27: design Tm conditions explicit and sealed --------------------------

@pytest.mark.parametrize("acc", ["M11240", "X05424"])
def test_explicit_tm_conditions_leave_primer3_output_unchanged(acc):
    import primer3
    seq = fixture_seq(acc)
    p = {"prod_min": 150, "prod_max": 500}
    got = primers.design_primers(seq, p)["candidates"]
    r = primers.resolved_design_parameters(p)
    bare = {"PRIMER_OPT_SIZE": r["opt_size"], "PRIMER_MIN_SIZE": r["min_size"], "PRIMER_MAX_SIZE": r["max_size"],
            "PRIMER_OPT_TM": r["opt_tm"], "PRIMER_MIN_TM": r["min_tm"], "PRIMER_MAX_TM": r["max_tm"],
            "PRIMER_MIN_GC": r["min_gc"], "PRIMER_MAX_GC": r["max_gc"],
            "PRIMER_PRODUCT_SIZE_RANGE": [150, 500], "PRIMER_NUM_RETURN": r["num_return"],
            "PRIMER_MAX_POLY_X": r["max_poly_x"], "PRIMER_GC_CLAMP": r["gc_clamp"]}
    ref = primer3.bindings.design_primers({"SEQUENCE_ID": "teagle", "SEQUENCE_TEMPLATE": seq}, bare)
    assert [(c["left_seq"], c["right_seq"], c["left_tm"], c["right_tm"]) for c in got] == [
        (ref[f"PRIMER_LEFT_{i}_SEQUENCE"], ref[f"PRIMER_RIGHT_{i}_SEQUENCE"],
         round(ref[f"PRIMER_LEFT_{i}_TM"], 1), round(ref[f"PRIMER_RIGHT_{i}_TM"], 1))
        for i in range(ref["PRIMER_PAIR_NUM_RETURNED"])]


def test_design_tm_conditions_are_sealed_and_reported_beside_the_qc_conditions():
    res = design(fixture_seq("M11240"))
    sealed = res["provenance"]["parameters"]["tm_conditions"]
    assert sealed == res["design_tm_conditions"] == primers.TM_CONDITIONS
    assert res["candidates"][0]["design_provenance"]["parameters"]["tm_conditions"] == sealed
    assert res["oligoqc_conditions"] == oligoqc.CONDITIONS
    assert res["oligoqc_conditions"] != res["design_tm_conditions"]
    sealed_keys = set(res["candidates"][0]["design_provenance"]["parameters"])
    assert not {"mv_conc", "viennarna", "oligoqc_conditions"} & sealed_keys


def test_qc_tm_is_labelled_with_its_conditions():
    assert oligoqc.qc_oligo("GACTGACTGTCAGTCAGGCT")["tm_basis"] == "qc_conditions"


# ---- primers2-4: one annealing core per oligo across every searched sequence ---------------------

def _unit(f, r, seed):
    return rnd(60, seed) + f + rnd(150, seed + 1) + RC(r) + rnd(60, seed + 2)


def _core_case():
    fwd, rev = rnd(26, 11), rnd(22, 12)
    bad = mutate(fwd, {0, 1, 2, 3, 4, 5, 6, 12})
    return fwd, rev, _unit(fwd, rev, 20), _unit(bad, rev, 30)


def _run(specimen, background, fwd, rev):
    body = {"sequence": ">s\n" + specimen, "fwd": fwd, "rev": rev,
            "params": {"prod_min": 70, "prod_max": 400}}
    if background:
        body["background"] = ">b\n" + background
    return engine.run_pcr(body)


def _call(a):
    return (a["fwd_anneal_len"], a["fwd_tail5"], a["fwd_mm"], a["rev_anneal_len"], a["rev_tail5"], a["rev_mm"])


def test_core_decision_is_the_same_for_a_site_in_specimen_and_background():
    fwd, rev, good, bad = _core_case()
    both = _run(good + rnd(1000, 40) + bad, bad, fwd, rev)
    from_bad = [a for a in both["amplicons"] if a["fwd_mm"] or a["fwd_tail5"]]
    assert from_bad == [], "the imperfect copy was dropped in the specimen but reported in the background"
    assert len(both["amplicons"]) == 1
    same = _run(bad, bad, fwd, rev)
    assert len(same["amplicons"]) == 2
    assert _call(same["amplicons"][0]) == _call(same["amplicons"][1]) and same["amplicons"][0]["fwd_tail5"] == 6


@pytest.mark.parametrize("which", ["good", "bad"])
def test_specimen_result_does_not_depend_on_the_backgrounds(which):
    fwd, rev, good, bad = _core_case()
    specimen = {"good": good, "bad": bad}[which]
    alone = _run(specimen, None, fwd, rev)
    for bg in (good, bad, good + rnd(1000, 41) + bad):
        with_bg = _run(specimen, bg, fwd, rev)
        mine = [a for a in with_bg["amplicons"] if a["source"] == alone["amplicons"][0]["source"]]
        assert mine == alone["amplicons"], "a background changed the specimen's own products"
        assert with_bg["status"] == alone["status"]


def test_background_with_a_longer_site_is_still_reported_at_the_specimen_core():
    fwd, rev, good, bad = _core_case()
    res = _run(bad, good, fwd, rev)
    assert [a["fwd_anneal_len"] for a in res["amplicons"]] == [20, 20]
    assert res["criteria"]["annealing_core"]["F"]["forward_strand"] == 20


def test_core_lengths_are_reported():
    fwd, rev, good, bad = _core_case()
    res = _run(good + rnd(1000, 40) + bad, bad, fwd, rev)
    assert res["criteria"]["annealing_core"]["F"]["forward_strand"] == 26


def test_background_without_a_specimen_site_keeps_its_own_core():
    fwd, rev, good, bad = _core_case()
    res = _run(rnd(300, 90), bad, fwd, rev)
    assert [a["fwd_tail5"] for a in res["amplicons"]] == [6]


def test_a_lone_search_keeps_the_per_sequence_core():
    fwd, rev, _good, bad = _core_case()
    amps = primers.in_silico_pcr(fwd, rev, bad, "b", prod_min=70, prod_max=400)
    assert len(amps) == 1 and amps[0]["fwd_tail5"] == 6 and amps[0]["fwd_mm"] == 2


# ---- primers2-5: design and sealed PCR apply one geometry rule ------------------------------------

def test_product_range_that_forces_overlapping_primers_is_rejected_at_design():
    with pytest.raises(BadRequest) as e:
        design(fixture_seq("M11240"), prod_min=30, prod_max=40)
    assert "overlap" in str(e.value).lower()


def test_every_designed_pair_passes_sealed_pcr():
    seq = fixture_seq("M11240")
    res = design(seq, prod_min=60, prod_max=80)
    assert res["candidates"]
    for c in res["candidates"]:
        out = engine.run_pcr({
            "sequence": ">t\n" + seq, "fwd": c["left_seq"], "rev": c["right_seq"],
            "primer3_target": c["primer3_target"], "included_region": c["included_region"],
            "intended_product_span": c["intended_product_span"],
            "intended_endpoint_tolerance_bp": c["intended_endpoint_tolerance_bp"],
            "design_provenance": c["design_provenance"], "design_sequence": ">t\n" + seq,
            "params": {"prod_min": 40, "prod_max": 120}})
        assert out["criteria"]["intended_product_status"] == "assessed"


# ---- primers2-6: the exact 3' window is honest ---------------------------------------------------

def test_strict_window_beyond_min_anneal_is_reported_and_cautioned():
    fwd, rev = rnd(20, 51), rnd(20, 52)
    seq = rnd(200, 53) + mutate(fwd, {0}) + rnd(200, 54) + RC(rev) + rnd(60, 55)
    body = {"sequence": ">s\n" + seq, "fwd": fwd, "rev": rev,
            "params": {"max_mm": 0, "tp": 20, "prod_min": 70, "prod_max": 600}}
    res = engine.run_pcr(body)
    assert len(res["amplicons"]) == 1
    assert res["criteria"]["three_prime_strict"] == 20
    assert res["criteria"]["three_prime_strict_effective"] == primers.MIN_PRIMER_SIZE
    assert any("exact" in c for c in res["cautions"])
    body["params"]["tp"] = 5
    assert not engine.run_pcr(body)["cautions"]


# ---- primers2-7: no product shorter than the primers ---------------------------------------------

def _crossed():
    seq = rnd(240, 61)
    return seq, seq[100:120], RC(seq[95:115])


def test_crossed_three_prime_ends_make_no_product():
    seq, f, r = _crossed()
    assert primers.in_silico_pcr(f, r, seq, "t", max_mm=0, prod_min=1, prod_max=100) == []


def test_tailed_crossed_pair_is_not_rescued_by_its_tails():
    seq, f, r = _crossed()
    assert primers.in_silico_pcr("G" * 40 + f, "G" * 40 + r, seq, "t", max_mm=0, prod_min=70, prod_max=300) == []


def test_primer_length_self_primed_product_is_not_reported():
    f = "GACTGACTGTACAGTCAGTC"                                   # palindromic 20-mer
    assert RC(f) == f
    seq = rnd(60, 62) + f + rnd(60, 63)
    assert primers.in_silico_pcr(f, f, seq, "t", max_mm=0, prod_min=1, prod_max=100) == []


# ---- primers2-9: primer input is normalised and validated ----------------------------------------

def _pcr_body(fwd, rev):
    seq = rnd(60, 71) + "GCTAATTGACTTCAGCATAT" + rnd(150, 72) + RC("TTGGCCATTGGCACTGGCAT") + rnd(60, 73)
    return {"sequence": ">s\n" + seq, "fwd": fwd, "rev": rev}


@pytest.mark.parametrize("fwd", ["GCTAATTGACTTCAGCATAT ", " gctaattgacttc agcatat\n", "GCUAAUUGACUUCAGCAUAU"])
def test_primer_whitespace_case_and_u_are_normalised(fwd):
    res = engine.run_pcr(_pcr_body(fwd, "TTGGCCATTGGCACTGGCAT"))
    assert res["status"] == "complete"
    assert res["provenance"]["parameters"]["fwd"] == "GCTAATTGACTTCAGCATAT"


@pytest.mark.parametrize("bad", ["GCTAATTGACTT-CAGCATAT", "GCTAATTGACTT1CAGCATAT", "GCTAATTGACTTEQCAGCATAT"])
def test_non_iupac_primer_is_rejected_with_the_offending_characters(bad):
    with pytest.raises(BadRequest) as e:
        engine.run_pcr(_pcr_body(bad, "TTGGCCATTGGCACTGGCAT"))
    assert "IUPAC" in str(e.value) and "forward primer" in str(e.value)


def test_in_silico_pcr_normalises_whitespace_too():
    seq, f, r = _crossed()
    assert primers.in_silico_pcr(f + " ", r, seq, "t", max_mm=0) == primers.in_silico_pcr(f, r, seq, "t", max_mm=0)


# ---- primers2-10: an ambiguous template base is a mismatch whatever the primer base --------------

def test_primer_n_does_not_match_template_n():
    f, r = rnd(20, 81), rnd(20, 82)
    site = f[:9] + "N" + f[10:]
    seq = rnd(60, 83) + site + rnd(150, 84) + RC(r) + rnd(60, 85)
    degen = f[:9] + "N" + f[10:]
    assert primers.in_silico_pcr(degen, r, seq, "t", max_mm=0, prod_min=70, prod_max=400) == []
    amps = primers.in_silico_pcr(degen, r, seq, "t", max_mm=1, prod_min=70, prod_max=400)
    assert [a["fwd_mm"] for a in amps] == [1]


# ---- primers2-11: the caps fire only when something was actually skipped --------------------------

def _tandem(n):
    f, r = rnd(20, 91), rnd(20, 92)
    return f, r, "".join(f + "A" * 30 for _ in range(n)) + "T" * 10 + RC(r)


def test_site_cap_is_set_only_when_a_further_site_exists(monkeypatch):
    monkeypatch.setattr(primers, "MAX_SITES", 5)
    for n, expect in ((5, False), (6, True)):
        f, r, seq = _tandem(n)
        stats = {}
        primers.in_silico_pcr(f, r, seq, "t", max_mm=0, prod_min=1, prod_max=10000, stats=stats)
        assert stats["sites_capped"] is expect


def test_amplicon_cap_is_set_only_when_a_further_product_exists(monkeypatch):
    f, r = rnd(20, 93), rnd(20, 94)
    monkeypatch.setattr(primers, "MAX_AMPS", 4)
    for n, expect in ((4, False), (5, True)):
        seq = "".join(f + "A" * 30 for _ in range(1)) + "".join("A" * 30 + RC(r) for _ in range(n))
        stats = {}
        amps = primers.in_silico_pcr(f, r, seq, "t", max_mm=0, prod_min=1, prod_max=10000, stats=stats)
        assert len(amps) == min(n, 4)
        assert stats["amplicons_capped"] is expect


# ---- primers2-13: API defaults equal the documented Standard preset -------------------------------

def test_design_defaults_equal_the_standard_preset():
    d = primers.DESIGN_DEFAULTS
    assert (d["prod_min"], d["prod_max"]) == (150, 500)


# ---- primers2-16: a non-negative Primer3 dG is not a structure -----------------------------------

def test_positive_primer3_dg_is_no_structure():
    assert oligoqc._p3_kcal(SimpleNamespace(structure_found=True, dg=400.0, tm=10.0)) == (None, None)
    assert oligoqc._p3_kcal(SimpleNamespace(structure_found=True, dg=-3000.0, tm=10.0))[0] == -3.0


def test_real_pair_with_positive_end_stability_reports_none(monkeypatch):
    a, b = "TGCTGCGAGATTTGTGGAGT", "GTCACCACTGTGGCAGAGAA"
    for x, y in ((a, b), (b, a)):
        dg = oligoqc._p3_end_stability(x, y)
        assert dg is None or dg < 0
    end = oligoqc.qc_pair(a, b)["end_stability"]
    assert end["p3"] is None or end["p3"] < 0
    if end["p3"] is None:
        assert end["agree"] == "none"


def test_sealed_design_on_an_ambiguous_template_passes_sealed_pcr():
    seq = fixture_seq("M11240")[:800]
    amb = seq[:300] + "RYKM" + seq[304:]
    c = design(amb)["candidates"][0]
    out = engine.run_pcr({
        "sequence": ">t\n" + amb, "fwd": c["left_seq"], "rev": c["right_seq"],
        "primer3_target": c["primer3_target"], "included_region": c["included_region"],
        "intended_product_span": c["intended_product_span"],
        "intended_endpoint_tolerance_bp": c["intended_endpoint_tolerance_bp"],
        "design_provenance": c["design_provenance"], "design_sequence": ">t\n" + amb,
        "params": {"prod_min": 100, "prod_max": 600}})
    assert out["criteria"]["intended_product_status"] == "assessed"
    assert c["design_provenance"]["parameters"]["template_ambiguity_masked"] == 4


def test_in_silico_pcr_rejects_foreign_cores():
    seq, f, r = _crossed()
    with pytest.raises(ValueError):
        primers.in_silico_pcr(f, r, seq, "t", cores={"F": 20})
