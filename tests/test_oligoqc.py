"""Primer/oligo secondary-structure QC (oligoqc): dual-engine ΔG (primer3 + ViennaRNA), kcal/mol units,
flag tiers, engine concordance, and graceful degradation when ViennaRNA is absent."""
import os, sys
import pytest

_BE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "app", "backend")
if _BE not in sys.path:
    sys.path.insert(0, _BE)
from teagle_core import oligoqc


def test_engines_available():
    a = oligoqc.available()
    assert a["primer3"] is True                          # primer3 is a hard dependency
    assert a["primer3_version"] != "unavailable"
    # ViennaRNA is optional; if present it must report a version
    if a["viennarna"]:
        assert a["viennarna_version"] != "unavailable"


def test_dg_is_kcal_per_mol_not_cal():
    # a real hairpin ΔG is single/low-double digits kcal/mol, NEVER thousands (the cal/mol -> kcal/mol bug)
    o = oligoqc.qc_oligo("GCGCGCGCATATATGCGCGCGC")
    dg = o["hairpin"]["p3"]
    assert dg is not None and -40.0 < dg < 0.0           # negative, and in kcal/mol range


def test_hairpin_prone_oligo_flags():
    o = oligoqc.qc_oligo("GCGCGCGCATATATGCGCGCGC")        # strong stem-loop
    assert o["hairpin"]["flag"] in ("caution", "warn")
    assert o["hairpin"]["p3"] <= oligoqc._HAIRPIN_CAUTION


def test_benign_oligo_is_ok():
    o = oligoqc.qc_oligo("ATCTGGCGGCGGAGTGGGCG")          # a real GAPDH-ish primer, no strong structure
    assert o["ok"] and o["hairpin"]["flag"] == "ok"      # weak/no hairpin
    assert 40.0 <= o["gc"] <= 80.0 and o["len"] == 20


def test_three_prime_complementary_pair_warns():
    # 3'-ends complementary -> hetero-dimer + end-stability should escalate
    r = oligoqc.qc_pair("AAAAAAAAAAGGGCCCTTTT", "AAAAAAAAAAAAAGGGCCC")
    assert r["ok"]
    assert r["hetero_dimer"]["flag"] in ("caution", "warn")
    assert r["worst"] in ("caution", "warn")


def test_pair_shape_and_worst_rollup():
    r = oligoqc.qc_pair("ATCTGGCGGCGGAGTGGGCG", "GCCGCCTACGCCACCAAGAC")
    for k in ("left", "right", "hetero_dimer", "end_stability", "worst", "conditions", "engines"):
        assert k in r
    assert r["worst"] in ("ok", "caution", "warn")
    assert r["conditions"]["temp_c"] == 37.0             # sealed, IDT-comparable conditions


def test_gc_clamp_and_polyx():
    assert oligoqc.gc_clamp("AAAAAGGCGC") == 5           # last 5 = GGCGC -> 5 G/C
    assert oligoqc.gc_clamp("AAAAATAGCT") == 2           # last 5 = TAGCT -> G,C
    assert oligoqc.gc_clamp("GGGGGAAAAA") == 0           # last 5 = AAAAA
    assert oligoqc.longest_poly_x("ATTTTTGC") == 5
    assert oligoqc.longest_poly_x("ACGT") == 1


def test_tier_thresholds():
    assert oligoqc._tier(None, -9, -5) == "ok"           # no structure
    assert oligoqc._tier(-3.0, -9, -5) == "ok"           # weaker than caution
    assert oligoqc._tier(-6.0, -9, -5) == "caution"
    assert oligoqc._tier(-10.0, -9, -5) == "warn"


def test_engine_agreement_labels():
    assert oligoqc._agree(-8.0, -7.0) == "agree"         # within band
    assert oligoqc._agree(-28.0, -9.0) == "disagree"     # far apart (the palindrome case)
    assert oligoqc._agree(-5.0, None) == "single"        # one engine only
    assert oligoqc._agree(None, None) == "none"


def test_worst_of_engines_drives_flag():
    # if EITHER engine says warn, the shown flag is warn (a concerning call is never hidden)
    m = oligoqc._metric(-10.0, -1.0, oligoqc._DIMER_WARN, oligoqc._DIMER_CAUTION)
    assert m["flag"] == "warn" and m["agree"] == "disagree"


def test_primer_seal_stable_across_qc_engines(monkeypatch):
    # the advisory secondary-structure QC must NEVER change the Primer3 design seal: same primers seal
    # identically whether or not ViennaRNA is present, and the QC method refs are attached UNSEALED
    import os as _os, sys as _sys
    _BEp = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "app", "backend")
    if _BEp not in _sys.path:
        _sys.path.insert(0, _BEp)
    _sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "tests"))
    import engine
    from helpers import fixture_seq
    body = {"sequence": fixture_seq("M11240"), "params": {"prod_min": 150, "prod_max": 500}}
    seal_with = engine.run_primers(body)["provenance"]["manifestSha256"]
    monkeypatch.setattr(oligoqc, "RNA", None)
    r = engine.run_primers(body)
    assert r["provenance"]["manifestSha256"] == seal_with            # seal unaffected by QC engine presence
    assert [x["key"] for x in r["provenance"]["references"]] == ["Primer3"]   # only the design method is sealed
    assert [x["key"] for x in r["oligoqc_references"]][0] == "SantaLucia1998"  # QC refs reported separately


def test_qc_pair_guards_empty_primer():
    assert oligoqc.qc_pair("", "GCCGCCTACGCCACCAAGAC")["ok"] is False       # empty primer -> guarded, no calc_* on ""
    assert oligoqc.qc_pair("ATCTGGCGGCGGAGTGGGCG", "")["ok"] is False


def test_graceful_without_viennarna(monkeypatch):
    # ViennaRNA absent -> cross-check returns None, primer3 primary still works, flags still computed
    monkeypatch.setattr(oligoqc, "RNA", None)
    o = oligoqc.qc_oligo("GCGCGCGCATATATGCGCGCGC")
    assert o["ok"] and o["hairpin"]["p3"] is not None and o["hairpin"]["vrna"] is None
    assert o["hairpin"]["agree"] == "n/a"                     # ViennaRNA not installed -> cross-check couldn't run (not 'single')
    r = oligoqc.qc_pair("ATCTGGCGGCGGAGTGGGCG", "GCCGCCTACGCCACCAAGAC")
    assert r["ok"] and r["hetero_dimer"]["vrna"] is None
