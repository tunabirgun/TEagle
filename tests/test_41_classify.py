"""4.1.0 classification-core regressions: retrovirus/ERV withhold, RVT_3 recode, withheld-confidence ladder,
completeness kind, nested integrases, strand-aware tails, failed-scan withholding, wording. Hermetic."""
import os
import random
import re

import pytest

from teagle_core import classify, domains, retroviral, sequtil, structural
from helpers import fixture_seq

LTR = [{"type": "LTR (terminal direct repeat)", "five_prime": [0, 968], "three_prime": [8504, 9472],
        "element_span": [0, 9472], "identity": 98.0}]
RETRO_UNDETERMINED = "LTR retrotransposon (superfamily undetermined)"


def _h(code, hmm, s, e, strand="+", orf=0, score=50.0, cls="retro"):
    return {"domain": code, "hmm": hmm, "label": code, "class": cls, "nt": [s, e], "strand": strand,
            "orf": orf, "score": score}


def _tir(n=25):
    return [{"type": "TIR (terminal inverted repeat)", "tir_len": n, "five_prime": [0, n],
             "three_prime": [900, 900 + n], "element_span": [0, 900 + n], "identity": 96.0}]


def _herv(env_hmm="HERV-K_env_2", gag_hmm="Gag_p24"):
    return [_h("GAG", gag_hmm, 1141, 1390), _h("PR", "RVP", 3693, 3762),
            _h("RT", "RVT_1", 4130, 4613, score=200.0), _h("RNaseH", "RNase_H", 4700, 4900),
            _h("INT", "rve", 5804, 6095), _h("ENV", env_hmm, 6792, 8499)]


def _text(cl):
    return " ".join([cl["explanation"], *cl["evidence"]])


# ---- classify-1 / owner decision (a) -----------------------------------------------------------------
def test_hervk_layout_withholds_the_superfamily_instead_of_naming_gypsy():
    cl = classify.classify(LTR, _herv())
    assert cl["te_class"] == "LTR/unclassified"
    assert not cl["superfamily"].startswith("Gypsy")
    assert cl["superfamily"] == classify.RETROVIRAL_UNDETERMINED
    assert cl["is_erv"] is True
    assert cl["confidence"] in ("Moderate", "Candidate")
    assert "errantivirus" not in _text(cl)
    assert any("retrovirus-specific" in e for e in cl["evidence"])
    assert not cl["superfamily"].endswith("chromovirus")


@pytest.mark.parametrize("code,hmm", [("GAG", "Gag_p24"), ("GAG", "Gag_p24_C"), ("GAG", "Gag_p10"),
                                      ("ENV", "HERV-K_env_2"), ("ENV", "GP41")])
def test_each_retrovirus_specific_profile_withholds_gypsy(code, hmm):
    doms = [_h("RT", "RVT_1", 4000, 4500, score=150), _h("INT", "rve", 5000, 5300), _h(code, hmm, 6000, 6400)]
    cl = classify.classify(LTR, doms)
    assert cl["te_class"] == "LTR/unclassified" and cl["superfamily"] == classify.RETROVIRAL_UNDETERMINED
    assert cl["confidence"] != "High"


@pytest.mark.parametrize("code,hmm", [("GAG", "PEG10_N-capsid"), ("GAG", "zf-CCHC_5"), ("ENV", "TLV_coat")])
def test_gypsy_type_profiles_keep_the_gypsy_call_and_errantivirus_wording(code, hmm):
    doms = [_h("RT", "RVT_1", 4000, 4500, score=150), _h("INT", "rve", 5000, 5300), _h(code, hmm, 6000, 6400)]
    cl = classify.classify(LTR, doms)
    assert cl["te_class"] == "LTR/Gypsy" and cl["superfamily"].startswith("Gypsy")
    assert cl["confidence"] == "High"
    assert ("errantivirus" in _text(cl)) == (code == "ENV")


def test_copia_order_is_unaffected_by_a_retroviral_profile():
    doms = [_h("GAG", "Gag_p24", 100, 400), _h("INT", "rve", 1000, 1300), _h("RT", "RVT_2", 2000, 2600, score=200)]
    assert classify.classify(LTR, doms)["te_class"] == "LTR/Copia"


def test_withheld_retroviral_label_with_tsd_does_not_print_an_ltr_target_site_claim():
    tsd = {"type": "TSD (target-site duplication)", "length": 5, "upstream": [0, 5], "downstream": [9472, 9477]}
    cl = classify.classify(LTR + [tsd], _herv())
    assert cl["te_class"] == "LTR/unclassified"
    assert "LTR elements duplicate" not in _text(cl)


def test_withheld_retroviral_call_keeps_the_architecture_and_the_erv_transcript_model():
    cl = classify.classify(LTR, _herv())
    assert cl["order"] == "GAG–PR–RT–RNaseH–INT–ENV"
    plain = [{**d, "nt": list(d["nt"])} for d in _herv()]
    arch = retroviral.transcript_architecture(LTR, plain, cl, 9472)
    assert arch is not None and arch["intron"] == [1141, 6792]


@pytest.mark.parametrize("acc,te_class", [("M11240", "LTR/Copia"), ("M12927", "LTR/Gypsy"), ("M80343", "LINE"),
                                          ("X01005", "DNA/Tc1-Mariner"), ("X05424", "DNA/hAT")])
def test_regression_anchors_keep_their_calls(acc, te_class):
    s = fixture_seq(acc)
    cl = classify.classify(structural.detect_all(s), domains.scan_domains(s), seq=s)
    assert cl["te_class"] == te_class


# ---- classify-3 --------------------------------------------------------------------------------------
@pytest.mark.parametrize("extra,second", [
    ("tir12", None), ("tir25", None), (None, _h("EN", "Exo_endo_phos", 3000, 3300, orf=1, cls="retro:LINE?")),
    (None, _h("HEL", "Helitron_like_N", 3000, 3300, orf=1, cls="dna:Helitron"))])
def test_withheld_lone_rt_is_never_moderate(extra, second):
    st = _tir(12) if extra == "tir12" else _tir(25) if extra == "tir25" else []
    doms = [_h("RT", "RVT_1", 100, 1200, score=150)] + ([second] if second else [])
    cl = classify.classify(st, doms)
    assert cl["outcome"] == "withheld"
    assert cl["confidence"] == "Candidate"
    assert "(candidate confidence)" in cl["explanation"]


# ---- classify-4 / classify-11 ------------------------------------------------------------------------
def test_solo_copia_ltr_is_not_filed_under_dna_transposon_in_the_completeness_line():
    s = fixture_seq("M11240")[:276]
    cl = classify.classify(structural.detect_all(s), domains.scan_domains(s), seq=s)
    assert cl["te_class"] == "repeat/structural-only"
    assert "DNA transposon" not in cl["completeness"]["kind"]


def test_long_tir_without_domains_is_still_a_dna_transposon_kind():
    cl = classify.classify(_tir(25), [])
    assert cl["te_class"] == "DNA/non-autonomous" and cl["completeness"]["kind"] == "DNA transposon"


@pytest.mark.parametrize("n,dom", [(25, _h("GAG", "Gag_p24", 300, 600)),
                                   (12, _h("EN", "Exo_endo_phos", 300, 600, cls="retro:LINE?"))])
def test_tir_plus_non_transposase_domain_kind_follows_the_call(n, dom):
    cl = classify.classify(_tir(n), [dom])
    assert "DNA transposon" not in cl["completeness"]["kind"]


def test_helitron_evidence_does_not_claim_tir_absence_when_a_tir_was_detected():
    cl = classify.classify(_tir(25), [_h("HEL", "Helitron_like_N", 300, 700, cls="dna:Helitron")])
    assert "absence of TIR" not in _text(cl)
    assert "terminal inverted repeat was detected" in _text(cl)


def test_long_tir_with_a_host_like_domain_keeps_the_class_ii_candidate_and_mentions_the_tir():
    cl = classify.classify(_tir(31), [_h("CHR", "Chromo", 300, 450)])
    assert cl["te_class"] == "DNA/non-autonomous"
    assert "31 bp" in _text(cl) and "CHR" in _text(cl)


# ---- classify-5 --------------------------------------------------------------------------------------
def test_rt_plus_chromodomain_uses_the_ltr_ledger():
    cl = classify.classify([], [_h("RT", "RVT_1", 100, 700, score=150), _h("CHR", "Chromo", 800, 950)])
    comp = cl["completeness"]
    assert cl["te_class"] == "LTR/unclassified"
    assert "LINE" not in comp["kind"] and "LINE" not in comp["tier"]
    assert not {"ORF1", "EN", "3′ tail (poly-A / poly-T)"} & set(comp["missing"])
    assert "CHR" in comp["present"]


def test_line_ledger_is_unchanged_for_an_orf1_line():
    comp = classify.classify([], [_h("RT", "RVT_1", 500, 1200, score=150), _h("ORF1", "Transposase_22", 10, 90)])["completeness"]
    assert "LINE" in comp["kind"]


# ---- classify-6 --------------------------------------------------------------------------------------
def _nested(int_up, int_down):
    return [_h("RT", "RVT_1", 3000, 3500, orf=0, score=90), int_down, int_up]


def test_integrase_in_the_rts_own_orf_beats_a_stronger_integrase_elsewhere():
    doms = _nested(_h("INT", "rve", 1000, 1300, orf=2, score=34), _h("INT", "rve", 4700, 5000, orf=0, score=23))
    cl = classify.classify(LTR, doms)
    assert cl["te_class"] == "LTR/Gypsy"


def test_bracketing_integrases_in_other_orfs_are_conflicting_not_a_confident_copia_or_gypsy():
    doms = _nested(_h("INT", "rve", 1000, 1300, orf=2, score=34), _h("INT", "rve", 4700, 5000, orf=1, score=23))
    cl = classify.classify(LTR, doms)
    assert cl["te_class"] == "LTR/unclassified"
    assert cl["superfamily"] == RETRO_UNDETERMINED
    assert "conflicting" in _text(cl)
    assert cl["confidence"] != "High"
    assert classify.classify(LTR, [d for d in doms if d["score"] != 23])["te_class"] == "LTR/Copia"


def test_single_integrase_stays_high():
    cl = classify.classify(LTR, [_h("RT", "RVT_1", 3000, 3500, score=90), _h("INT", "rve", 4700, 5000, score=23)])
    assert (cl["te_class"], cl["confidence"]) == ("LTR/Gypsy", "High")


# ---- classify-8 / classify-18 / owner decision (b) ---------------------------------------------------
_CODON = dict(zip("ACDEFGHIKLMNPQRSTVWYX", "GCT TGC GAT GAA TTT GGT CAT ATT AAA CTG ATG AAT CCG CAG CGT AGC ACC GTT TGG TAT GCT".split()))


def _back_translate(name):
    hmm = next(h for h in domains._hmms() if str(h.name) == name)
    return "ATG" + "".join(_CODON[a] for a in hmm.consensus.upper()) + "TAA"


def test_rvt3_is_coded_as_an_rnase_h_like_domain():
    code, label, _cls, acc = domains.DOMAIN_INFO["RVT_3"]
    assert (code, acc) == ("RNaseH", "PF13456") and "RNase H" in label


@pytest.mark.parametrize("name", ["RNase_H", "RVT_3"])
def test_an_rnase_h_fold_alone_is_never_reported_as_reverse_transcriptase_or_class_i(name):
    seq = _back_translate(name)
    hits = domains.scan_domains(seq)
    assert hits and "RT" not in {h["domain"] for h in hits}
    assert not classify.classify(structural.detect_all(seq), hits, seq=seq)["class"].startswith("Class I")


def test_integrase_consensus_does_not_also_name_a_transposase():
    seq = _back_translate("rve")
    hits = domains.scan_domains(seq)
    assert "INT" in {h["domain"] for h in hits}
    assert "TPase" not in {h["domain"] for h in hits}
    assert not classify.classify([], hits, seq=seq)["class"].startswith("Class II")


def _hit(code, s, e, score, strand="+", hmm="x"):
    return {"domain": code, "hmm": hmm, "nt": [s, e], "strand": strand, "score": score}


def test_int_and_tpase_overlap_keeps_the_stronger_only():
    kept = domains._dedup_domains([_hit("INT", 0, 300, 60), _hit("TPase", 30, 300, 30)])
    assert [k["domain"] for k in kept] == ["INT"]
    kept = domains._dedup_domains([_hit("INT", 0, 300, 30), _hit("TPase", 30, 300, 90)])
    assert [k["domain"] for k in kept] == ["TPase"]


def test_int_tpase_shadowing_is_narrow():
    assert len(domains._dedup_domains([_hit("INT", 0, 300, 60), _hit("TPase", 600, 900, 30)])) == 2
    assert len(domains._dedup_domains([_hit("INT", 0, 300, 60), _hit("TPase", 0, 300, 30, strand="-")])) == 2
    assert len(domains._dedup_domains([_hit("RNaseH", 0, 300, 60), _hit("TPase", 0, 300, 30)])) == 2


# ---- classify-9 --------------------------------------------------------------------------------------
@pytest.mark.parametrize("acc,arch,wicker", [("M80343", "ORF1–EN–RT", "LINE"), ("M11240", "INT–RT", "LTR"),
                                             ("X01005", "TPase", "TIR")])
def test_domain_architecture_and_wicker_order_are_separate_fields(acc, arch, wicker):
    s = fixture_seq(acc)
    cl = classify.classify(structural.detect_all(s), domains.scan_domains(s), seq=s)
    assert (cl["domain_architecture"], cl["wicker_order"]) == (arch, wicker)
    assert "tail" not in cl["domain_architecture"]


def test_wicker_order_states_for_withheld_unassigned_and_failed_scans():
    assert classify.classify(LTR, [])["wicker_order"] == "withheld"
    assert classify.classify([], [])["wicker_order"] == "unassigned"
    assert classify.classify([], [_h("GAG", "Gag_p24", 10, 300)])["wicker_order"] == "unassigned"
    cl = classify.classify(LTR, [], domains_ok=False)
    assert cl["domain_architecture"] == "unavailable" and cl["wicker_order"] == "withheld"
    assert classify.classify([], [])["domain_architecture"] == "none detected"
    assert classify.classify([], [_h("TPase", "x", 1, 300, cls="dna:hAT")], structural_ok=False)["domain_architecture"] == "TPase"


# ---- classify-10 -------------------------------------------------------------------------------------
_DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "app", "backend", "data")


def test_data_readme_table_equals_the_bundled_profile_file():
    hmm_text = open(os.path.join(_DATA, "te_domains.hmm"), encoding="utf-8").read()
    file_rows = re.findall(r"^NAME\s+(\S+)\nACC\s+(PF\d+\.\d+)", hmm_text, flags=re.M)
    readme = open(os.path.join(_DATA, "README.md"), encoding="utf-8").read()
    rows = re.findall(r"^\| (\S+) \| (PF\d+\.\d+) \| (\S+) \|", readme, flags=re.M)
    assert [r[:2] for r in rows] == file_rows
    assert [r[0] for r in rows if r[0] in domains.DOMAIN_INFO] == [r[0] for r in rows]
    assert {r[0]: r[2] for r in rows} == {k: v[0] for k, v in domains.DOMAIN_INFO.items()}
    assert f"{len(domains.DOMAIN_INFO)} concatenated" in readme


# ---- classify-12 / classify-17 -----------------------------------------------------------------------
def _line(strand="+"):
    if strand == "+":
        return [_h("ORF1", "Transposase_22", 10, 90), _h("EN", "Exo_endo_phos", 100, 400, cls="retro:LINE?"),
                _h("RT", "RVT_1", 500, 1200, score=150)]
    return [_h("ORF1", "Transposase_22", 1500, 1590, strand="-"),
            _h("EN", "Exo_endo_phos", 700, 1000, strand="-", cls="retro:LINE?"),
            _h("RT", "RVT_1", 100, 600, strand="-", score=150)]


_POLY_A = {"type": "poly-A tail (3')", "length": 20, "pos": [2000, 2020]}
_POLY_T = {"type": "poly-T (5')", "length": 20, "pos": [0, 20]}


@pytest.mark.parametrize("strand,tail,credited", [("+", _POLY_A, True), ("+", _POLY_T, False),
                                                   ("-", _POLY_T, True), ("-", _POLY_A, False)])
def test_line_tail_is_credited_only_on_the_strand_it_belongs_to(strand, tail, credited):
    cl = classify.classify([tail], _line(strand))
    assert cl["te_class"] == "LINE"
    assert ("intact" in cl["completeness"]["tier"]) is credited
    assert (cl["confidence"] == "Moderate") is credited
    assert any("not credited" in e for e in cl["evidence"]) is (not credited)


def test_polya_signal_motif_row_is_not_counted_as_a_tail():
    motif = {"type": "polyA-signal (motif, advisory)", "pos": [2000, 2006]}
    cl = classify.classify([motif], [_h("ORF1", "Transposase_22", 10, 90), _h("RT", "RVT_1", 500, 1200, score=150)])
    assert cl["te_class"] == "LINE" and cl["confidence"] == "Candidate"


# ---- classify-13 -------------------------------------------------------------------------------------
@pytest.mark.parametrize("structure_ok", [True, False])
def test_dirs_needs_rt_and_recombinase_on_one_strand(structure_ok):
    doms = [_h("RT", "RVT_1", 100, 1200, score=150), _h("YR", "Phage_integrase", 9000, 9400, strand="-", cls="retro:YR?")]
    cl = classify.classify([], doms, structural_ok=structure_ok)
    assert cl["te_class"] != "DIRS" and "DIRS" not in (cl["superfamily"] or "")
    same = [doms[0], {**doms[1], "strand": "+"}]
    assert classify.classify([], same, structural_ok=structure_ok)["te_class"] == "DIRS"


# ---- classify-14 -------------------------------------------------------------------------------------
def test_line_coding_evidence_plus_direct_repeats_is_a_conflict_not_an_ltr_call():
    cl = classify.classify(LTR, _line("+"))
    assert cl["te_class"] == "retro/unassigned" and cl["confidence"] == "Candidate"
    text = _text(cl)
    assert "ORF1p" in text and "direct repeats" in text
    assert {"ORF1", "EN"} <= set(cl["completeness"]["present"])


# ---- classify-16 -------------------------------------------------------------------------------------
def test_terminal_direct_repeats_do_not_lift_a_transposase_call():
    tp = _h("TPase", "Dimer_Tnp_hAT", 200, 800, cls="dna:hAT", score=90)
    assert classify.classify(LTR, [tp])["confidence"] == "Candidate"
    assert classify.classify(_tir(30), [tp])["confidence"] == "Moderate"


# ---- engine_state-5 / spec_conformance-4 / -6 --------------------------------------------------------
def test_scan_domains_raises_when_pyhmmer_is_unavailable(monkeypatch):
    monkeypatch.setattr(domains, "pyhmmer", None)
    monkeypatch.setattr(domains, "PYHMMER_ERROR", "ImportError: simulated")
    with pytest.raises(RuntimeError, match="pyhmmer unavailable"):
        domains.scan_domains("ATG" + "GCT" * 200 + "TAA")


def test_engine_reports_the_domain_stage_as_failed_when_pyhmmer_is_unavailable(monkeypatch):
    import engine
    monkeypatch.setattr(domains, "pyhmmer", None)
    rec = engine.run_analyze({"sequence": ">x\n" + fixture_seq("M11240")})["records"][0]
    assert rec["processing"]["domains"]["status"] == "failed"
    cl = rec["classification"]
    assert cl["domains_unavailable"] is True and "none detected in tested panel" not in cl["explanation"]
    assert any("domain scan unavailable" in n for n in rec["notes"])


@pytest.mark.parametrize("st", [_tir(31), _tir(12), [{"type": "poly-A tail (3')", "length": 20, "pos": [0, 20]}]])
def test_failed_domain_scan_withholds_every_structure_only_call(st):
    cl = classify.classify(st, [], domains_ok=False)
    wording = _text(cl).lower()
    assert cl["outcome"] == "withheld" and cl["explanation"].startswith("Classification withheld")
    assert not cl["class"].startswith("Class II")
    assert "no transposase" not in wording and "non-autonomous" not in wording and "no clear te" not in wording
    assert "no coding domain" not in " ".join(str(v) for v in (cl["completeness"] or {}).values()).lower()
    assert cl["confidence"] == "Candidate"


def test_working_domain_scan_still_calls_a_long_tir_non_autonomous():
    assert classify.classify(_tir(31), [])["te_class"] == "DNA/non-autonomous"


def test_transposase_without_structure_keeps_the_advisory_strongest_profile_family():
    tp = _h("TPase", "Dimer_Tnp_hAT", 100, 500, cls="dna:hAT", score=90)
    cl = classify.classify([], [tp], structural_ok=False)
    assert cl["te_class"] == "DNA/hAT" and "hAT" in cl["superfamily"] and "advisory" in cl["superfamily"]
    assert cl["confidence"] == "Candidate"
    assert cl["completeness"]["unassessed"] == ["TIR (5′)", "TIR (3′)"]
    two = [tp, _h("TPase", "DDE_3", 600, 900, cls="dna:Tc1-Mariner", score=40)]
    assert any("multiple transposase classes" in e for e in classify.classify([], two, structural_ok=False)["evidence"])


# ---- spec_conformance-17 -----------------------------------------------------------------------------
def test_chance_tir_sentence_quotes_its_constants_and_the_floor_sits_above_them():
    text = _text(classify.classify(_tir(31), []))
    assert f"{classify.CHANCE_TIR_MAX_BP} bp across {classify.CHANCE_TIR_SEQUENCES:,} random sequences" in text
    assert structural.MIN_TIR_STANDALONE > classify.CHANCE_TIR_MAX_BP
    assert "4,500" not in text


# ---- persona-U16 -------------------------------------------------------------------------------------
def test_coding_axis_counts_orfs_and_recognises_the_canonical_line_layout():
    s = fixture_seq("M80343")
    orfs = sequtil.find_orfs(s)
    cl = classify.classify(structural.detect_all(s), domains.scan_domains(s), seq=s, orfs=orfs)
    note = cl["completeness"]["coding"]["note"]
    assert "separate reading frame" not in note
    assert "2 separate ORFs" in note and "canonical LINE ORF1/ORF2 organisation" in note
    assert "degenerate copy" not in note.split("canonical")[0]


def test_coding_axis_never_claims_more_than_six_frames():
    orfs = [{"start": i * 1000, "end": i * 1000 + 900, "strand": "+", "frame": i % 3} for i in range(8)]
    doms = [_h("RT", "RVT_1", i * 1000 + 10, i * 1000 + 300, orf=i, score=60 + i) for i in range(8)]
    note = classify.classify([], doms, orfs=orfs)["completeness"]["coding"]["note"]
    assert "8 separate ORFs" in note
    assert not re.search(r"(?<!ORFs)\b[7-9] separate reading frame", note)


# ---- persona-U24 -------------------------------------------------------------------------------------
@pytest.mark.parametrize("doms,sf", [
    ([_h("INT", "rve", 100, 400), _h("RT", "RVT_2", 600, 1200, score=150)], "Copia"),
    ([_h("RT", "RVT_1", 100, 700, score=150), _h("INT", "rve", 800, 1100)], "Gypsy")])
def test_superfamily_evidence_does_not_present_a_superfamily_as_an_order(doms, sf):
    text = _text(classify.classify(LTR, doms))
    assert "Ty1 order" not in text and "Ty3 order" not in text
    assert f"{sf}" in text and "superfamily" in text



# ---- F3: wiring of the structural helpers ------------------------------------------------------------
def _flanked(acc, seed, tsd, flank=300):
    rng = random.Random(seed * 7919 + len(acc))
    rnd = lambda n: "".join(rng.choices("ACGT", k=n))
    left, right, el = rnd(flank), rnd(flank), fixture_seq(acc)
    return left + tsd + el + tsd + right, len(left) + len(tsd), len(left) + len(tsd) + len(el)


def _classified(seq):
    st = structural.detect_all(seq)
    return st, classify.classify(st, domains.scan_domains(seq), seq=seq)


def _row(rows, prefix):
    return next((r for r in rows if r["type"].startswith(prefix)), None)


def test_flanked_tc1_reports_the_trimmed_tir_and_its_ta_target_site():
    ok = []
    for seed in range(30):
        seq, s, e = _flanked("X01005", seed, "TA")
        st, cl = _classified(seq)
        tir, tsd = _row(st, "TIR"), _row(st, "TSD")
        ok.append(cl["te_class"] == "DNA/Tc1-Mariner" and tir["tir_len"] == 54 and tir["element_span"] == [s, e]
                  and bool(tsd) and tsd["motif"] == "TA" and tsd["upstream"] == [s - 2, s]
                  and tsd["downstream"] == [e, e + 2] and tsd["tsd_congruence"] == "congruent"
                  and len([r for r in st if r["type"].startswith("TSD")]) == 1)
    print(f"flanked Tc1 through classify: {sum(ok)}/{len(ok)}")
    assert sum(ok) / len(ok) >= 0.9


def test_target_site_trim_is_gated_on_the_tc1_call():
    seq, s, e = _flanked("X05424", 3, "TA")                       # an hAT-called TIR is never trimmed to expose a TA
    st, cl = _classified(seq)
    assert cl["te_class"] == "DNA/hAT"
    assert not any(r.get("boundary_adjustment") for r in st if r["type"].startswith("TIR"))


def test_bare_specimens_keep_their_anchors_through_the_tc1_wiring():
    for acc, te_class, conf in (("X01005", "DNA/Tc1-Mariner", None), ("X05424", "DNA/hAT", None),
                                ("M11240", "LTR/Copia", "High"), ("M12927", "LTR/Gypsy", "High"),
                                ("M80343", "LINE", "Moderate")):
        s = fixture_seq(acc)
        st, cl = _classified(s)
        assert cl["te_class"] == te_class and (conf is None or cl["confidence"] == conf), acc
    st, _ = _classified(fixture_seq("X01005"))
    assert _row(st, "TIR")["tir_len"] == 54 and _row(st, "TSD") is None


def _tsd_row(length, motif, **kw):
    return {"type": "TSD (target-site duplication)", "length": length, "motif": motif,
            "upstream": [100 - length, 100], "downstream": [930, 930 + length], **kw}


def test_the_tsd_motif_reaches_the_congruence_test():
    tc1 = _tpase("dna:Tc1-Mariner")
    ta = _tsd_row(2, "TA", chance_probability=0.06, chance_level=True)
    gc = _tsd_row(2, "GC", chance_probability=0.06, chance_level=True)
    classify.classify([*_tir(), ta], [tc1])
    ev = " ".join(classify.classify([*_tir(), gc], [tc1])["evidence"])
    assert ta["tsd_congruence"] == "congruent" and gc["tsd_congruence"] == "incongruent"
    assert "GC" in ev and "TA" in ev and "not credited" in ev


def _tpase(cls="dna:hAT"):
    return {"domain": "TPase", "hmm": "Dimer_Tnp_hAT", "nt": [200, 800], "strand": "+", "score": 120.0,
            "class": cls}


@pytest.mark.parametrize("cls,length,flag,captured", [
    ("dna:hAT", 8, False, True), ("dna:hAT", 8, True, False),
    ("dna:CACTA", 3, False, True), ("dna:CACTA", 3, True, False),
    ("dna:Tc1-Mariner", 2, True, False)])
def test_insertion_site_captured_requires_a_tsd_that_is_not_chance_level(cls, length, flag, captured):
    motif = "TA" if length == 2 else "A" * length
    cl = classify.classify([*_tir(), _tsd_row(length, motif, chance_probability=0.001, chance_level=flag)], [_tpase(cls)])
    text = _text(cl)
    assert ("insertion site itself is captured" in text) is captured
    assert ("both element termini are present" in text) is captured
    if not captured:
        assert "arise by chance" in text


def test_a_tsd_row_without_a_chance_estimate_is_not_credited_as_the_insertion_site():
    text = _text(classify.classify([*_tir(), _tsd_row(8, "A" * 8)], [_tpase()]))
    assert "insertion site itself is captured" not in text and "not assessed" in text


def test_tsd_expectation_is_quoted_as_the_range_the_congruence_test_uses():
    struct = [{"type": "LTR (terminal direct repeat)", "element_span": [0, 9000], "identity": 98.0},
              _tsd_row(2, "AC")]
    doms = [_h("INT", "rve", 100, 400), _h("RT", "RVT_2", 600, 1200, score=150)]
    text = _text(classify.classify(struct, doms))
    assert "4-6 bp" in text and "duplicate 5 bp" not in text


def _long(n, seed=1):
    return "".join(random.Random(seed).choices("ACGT", k=n))


def test_missing_tir_names_the_part_of_a_long_record_the_scan_did_not_read():
    long_text = _text(classify.classify([], [_tpase()], seq=_long(5000)))
    assert "no terminal inverted-repeat pair recovered" in long_text
    assert "5,000 bp" in long_text and "only within" in long_text
    assert "only within" not in _text(classify.classify([], [_tpase()], seq=_long(700)))
    assert "only within" not in _text(classify.classify([], [_tpase()]))


def test_missing_ltr_names_the_part_of_a_long_record_the_scan_did_not_read():
    doms = [_h("INT", "rve", 100, 400), _h("RT", "RVT_2", 600, 1200, score=150)]
    assert "only within 6,000 bp" in _text(classify.classify([], doms, seq=_long(30000)))
    assert "only within" not in _text(classify.classify([], doms, seq=_long(9000)))
