"""Retroviral transcript architecture (env splice) + LTR cis-elements (PBS / PPT).
Unit tests use mock structural/domain input; the @network test validates against real HERV-K113.
These check CORRECTNESS (strand handling, geometry guards, honest hedging), not just presence."""
import os, sys
import pytest
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "app", "backend"))
from teagle_core import retroviral, structural           # noqa: E402


def _ltr(s5=0, e5=968, s3=8504, e3=9472):
    return {"type": "LTR (terminal direct repeat)", "five_prime": [s5, e5],
            "three_prime": [s3, e3], "element_span": [s5, e3]}


def _dom(code, s, e, strand="+"):
    return {"domain": code, "nt": [s, e], "strand": strand}


def _erv_plus():
    return [_dom("GAG", 1141, 1390), _dom("PR", 3693, 3762), _dom("RT", 4130, 4613),
            _dom("INT", 5804, 6095), _dom("ENV", 6792, 8499)]


def _erv_minus(n=9472):
    # reverse-complement orientation of _erv_plus(): a +[a,b] domain maps to [n-b, n-a] on strand "-"
    return [_dom(d["domain"], n - d["nt"][1], n - d["nt"][0], "-") for d in _erv_plus()]


def test_erv_env_splice_architecture_plus_strand():
    arch = retroviral.transcript_architecture([_ltr()], _erv_plus(),
                                              {"is_erv": True, "superfamily": "Gypsy (Ty3)"}, 9472)
    assert arch is not None and arch["strand"] == "+"
    assert arch["leader_exon"] == [968, 1141]              # 5' LTR end -> gag start
    assert arch["intron"] == [1141, 6792]                  # gag-pro-pol -> env start (the single large intron)
    assert arch["env_exon"][0] == 6792
    assert arch["approximate"] is True
    assert "frameshift" not in arch["note"] or "readthrough" in arch["note"]   # mechanism not asserted as fact
    assert "e.g." in arch["subsplice_note"].lower()        # rec/np9 conditional, not asserted for this locus


def test_minus_strand_erv_builds_a_sane_mirrored_model_not_garbled():
    # a reverse-complement-oriented provirus is a normal coordinate-fetch input; the model must mirror,
    # not emit a zero-length leader / tiny intron with the domains sprawling the wrong way
    arch = retroviral.transcript_architecture([_ltr()], _erv_minus(),
                                              {"is_erv": True, "superfamily": "Gypsy (Ty3)"}, 9472)
    assert arch is not None and arch["strand"] == "-"
    lex, intr, eex = arch["leader_exon"], arch["intron"], arch["env_exon"]
    assert eex[0] < intr[0] < lex[0]                       # env low, intron middle, leader high (mirrored)
    assert all(b - a > 0 for a, b in (lex, intr, eex))     # no zero-length span
    assert intr[1] - intr[0] > 4000                        # the gag-pol intron is large, not a few bp


def test_env_before_gag_or_empty_body_returns_none():
    # + strand but the only body domain sits AFTER env (env-before-gag layout) -> empty body -> no model
    doms = [_dom("ENV", 1000, 2000), _dom("GAG", 5000, 6000), _dom("RT", 6100, 6500)]
    assert retroviral.transcript_architecture([_ltr()], doms, {"is_erv": True, "superfamily": "Gypsy"}, 9472) is None
    # env only, no gag/pol at all -> no model
    assert retroviral.transcript_architecture([_ltr()], [_dom("ENV", 6792, 8499)],
                                              {"is_erv": True, "superfamily": "Gypsy"}, 9472) is None


def test_non_erv_and_missing_ltr_get_no_model():
    assert retroviral.transcript_architecture([_ltr()], _erv_plus(), {"is_erv": False}, 9472) is None
    assert retroviral.transcript_architecture([], _erv_plus(), {"is_erv": True}, 9472) is None


def test_find_ppt_reports_a_clean_purine_tract():
    # leading pyrimidines / N before the purine run must NOT be included in the reported PPT
    seq = "C" * 90 + "TTCT" + "AAGAAAAGGGGGAAA" + "G" * 968     # 3' LTR starts at 109
    ppt = structural.find_ppt(seq, ltr_three_prime_start=109)
    assert ppt is not None and ppt["pos"][1] == 109
    assert ppt["motif"][0] in "AG"                              # tract starts on a purine
    assert "T" not in ppt["motif"][:1] and "N" not in ppt["motif"]
    assert ppt["purine_frac"] >= 0.82


def test_find_pbs_names_trna_only_when_confident():
    lys3 = "TGGCGCCCGAACAGGGAC"
    pbs = structural.find_pbs("N" * 968 + lys3 + "ACGT" * 20, ltr_five_prime_end=968)
    assert pbs["confident"] is True and pbs["priming_trna"] == "tRNA-Lys3" and pbs["identity"] == 100.0
    # a diverged (HERV-K-like) match is detected but NOT hard-named — priming tRNA undetermined
    div = "TGGTGCCCAACGTGGAGG"
    pbs2 = structural.find_pbs("N" * 968 + div + "ACGT" * 20, ltr_five_prime_end=968)
    assert pbs2["confident"] is False and pbs2["priming_trna"] == "undetermined"
    assert pbs2["best_match"] == "tRNA-Lys3" and "undetermined" in pbs2["note"]


@pytest.mark.network
def test_hervk113_live_architecture_and_cis_elements():
    from teagle_core import fetch, domains, classify
    s = fetch.retrieve("AY037928")["fasta"]
    s = (s.split("\n", 1)[1].replace("\n", "") if s.startswith(">") else s).upper()
    st = structural.detect_all(s)
    dm = domains.scan_domains(s)
    cl = classify.classify(st, dm)
    arch = retroviral.transcript_architecture(st, dm, cl, len(s))
    assert arch is not None and arch["strand"] == "+"
    assert arch["intron"][1] - arch["intron"][0] > 4000                        # large gag-pol intron
    pbs = next(e for e in st if e["type"].startswith("PBS"))
    assert pbs["best_match"].startswith("tRNA-Lys") and pbs["confident"] is False  # Lys closest, honestly hedged
    ppt = next(e for e in st if e["type"].startswith("PPT"))
    assert ppt["motif"][0] in "AG" and ppt["purine_frac"] >= 0.82              # clean purine tract
    # the same provirus reverse-complemented must NOT produce a garbled model
    rc = s.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    rst, rdm = structural.detect_all(rc), domains.scan_domains(rc)
    rarch = retroviral.transcript_architecture(rst, rdm, classify.classify(rst, rdm), len(rc))
    assert rarch is None or (rarch["intron"][1] - rarch["intron"][0] > 4000 and rarch["strand"] == "-")
