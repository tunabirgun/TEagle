from teagle_core import structural
from helpers import make_ltr_element, make_tir_element


def test_ltr_detected_and_deterministic():
    seq = make_ltr_element(seed=1, ltr_len=160, divergence=0.02)
    a = structural.find_ltr(seq)
    b = structural.find_ltr(seq)
    assert a is not None, "LTR should be detected on a synthetic LTR element"
    assert a == b, "detection must be deterministic"
    assert 140 <= a["ltr_len"] <= 165           # ~160 bp, allow heuristic boundary slack
    assert a["identity"] >= 95
    s0, s1 = a["element_span"]
    assert s0 < s1 and (s1 - s0) > 1500          # spans the whole element


def test_no_ltr_on_random_sequence():
    import random
    r = random.Random(99)
    seq = "".join(r.choice("ACGT") for _ in range(2000))
    assert structural.find_ltr(seq) is None      # no terminal direct repeat -> no false positive


def test_tir_detected():
    seq = make_tir_element(seed=7, tir_len=24)
    tir = structural.find_tir(seq)
    assert tir is not None
    assert 18 <= tir["tir_len"] <= 30
    assert tir["identity"] >= 95


def test_tir_boundary_perfect_vs_imperfect():
    import random
    from teagle_core.sequtil import reverse_complement
    rng = random.Random(11)
    core = "".join(rng.choice("ACGT") for _ in range(54))
    mid = "".join(rng.choice("ACGT") for _ in range(400))
    # a perfect 54 bp TIR must be reported at exactly 54, not extended into the internal sequence
    perfect = core + mid + reverse_complement(core)
    tp = structural.find_tir(perfect)
    assert tp and tp["tir_len"] == 54 and tp["identity"] == 100.0
    # an imperfect ~28 bp TIR (a few 3'-side mismatches) must extend to ~28, not collapse to a short perfect core
    tir28 = "".join(rng.choice("ACGT") for _ in range(28))
    three = list(reverse_complement(tir28))
    for i in (2, 9, 17, 23):
        three[i] = {"A": "C", "C": "A", "G": "T", "T": "G"}[three[i]]
    imperfect = tir28 + mid + "".join(three)
    ti = structural.find_tir(imperfect)
    assert ti and ti["tir_len"] >= 24        # extends the imperfect TIR, not collapsed to a short perfect core


def test_polya():
    seq = "ACGT" * 20 + "A" * 12
    hits = structural.find_polya(seq, min_run=8)
    assert any(h["type"].startswith("poly-A") and h["length"] >= 12 for h in hits)


def test_tsd_direct_repeat_flank():
    # element flanked by an identical 6 bp direct repeat = TSD
    tsd = "ACGTAC"
    elem = "GG" + "T" * 100 + "GG"
    seq = "N" * 3 + tsd + elem + tsd + "N" * 3
    start = 3 + len(tsd)
    end = start + len(elem)
    d = structural.find_tsd(seq, start, end)
    assert d is not None and d["motif"] == tsd and d["length"] == 6


def test_detect_all_shape():
    seq = make_ltr_element(seed=2)
    ev = structural.detect_all(seq)
    assert any(e["type"].startswith("LTR") for e in ev)
    for e in ev:
        assert "type" in e and "method" in e or e["type"].startswith(("TSD", "poly"))
