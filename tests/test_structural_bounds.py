"""Regression tests for the terminal-repeat measurement defects (2026-07-27).

Before this: find_ltr capped its seeding window at 1800 bp and reported whatever fraction of the
repeat fell inside both windows as a measured length (a 2.5 kb LTR read as 700 bp, a 4.4 kb LTR was
not detected at all); _terminal_tir stopped at 60 bp, so 100/200/400 bp TIRs all read as exactly 60;
and the anchor-density floor required ~93% identity, which silently dropped older LTR pairs and made
the 80% identity gate unreachable. None of it carried a flag — the truncated number was presented as
evidence and propagated into element_span, TSD/PBS/PPT, primer targets and the sealed manifest.
"""
import random

from teagle_core import structural
from teagle_core.sequtil import reverse_complement


def _rnd(n, seed):
    return "".join(random.Random(seed).choices("ATGC", k=n))


def _ltr_element(ltr_len, inner=4000, flank=200, seed=1):
    """A synthetic LTR element: flank + LTR + internal + identical LTR + flank."""
    r = random.Random(seed)
    ltr = "".join(r.choices("ATGC", k=ltr_len))
    return ("".join(r.choices("ATGC", k=flank)) + ltr + "".join(r.choices("ATGC", k=inner))
            + ltr + "".join(r.choices("ATGC", k=flank)))


# ---------- LTR length is measured, not window-limited ----------
def test_ltr_longer_than_the_old_window_is_measured():
    for ltr_len in (1800, 2500, 4400):
        ev = structural.find_ltr(_ltr_element(ltr_len, seed=ltr_len))
        assert ev is not None, f"{ltr_len} bp LTR not detected"
        # exact to a few bp: ungapped extension can run a base or two into a chance match past the end
        assert abs(ev["ltr_len"] - ltr_len) <= 5, (ltr_len, ev["ltr_len"])


def test_ltr_element_span_covers_both_copies():
    ltr_len, inner, flank = 2500, 4000, 200
    ev = structural.find_ltr(_ltr_element(ltr_len, inner, flank, seed=9))
    lo, hi = ev["element_span"]
    assert abs(lo - flank) <= 5
    assert abs(hi - (flank + 2 * ltr_len + inner)) <= 5


# ---------- divergent LTR pairs are still found ----------
def test_diverged_ltr_pair_is_detected():
    """LTR-LTR divergence is the basis of insertion dating, so the detector must survive it.
    The old anchor-density floor needed ~93% identity; 85% is an ordinary older insertion."""
    r = random.Random(31)
    ltr = "".join(r.choices("ATGC", k=1200))
    mutated = list(ltr)
    for i in r.sample(range(1200), 180):              # 15% substituted
        mutated[i] = r.choice("ACGT")
    seq = ("".join(r.choices("ATGC", k=200)) + ltr + "".join(r.choices("ATGC", k=3000))
           + "".join(mutated) + "".join(r.choices("ATGC", k=200)))
    ev = structural.find_ltr(seq)
    assert ev is not None, "85%-identity LTR pair not detected"
    assert abs(ev["ltr_len"] - 1200) <= 30
    assert ev["identity"] >= 80


# ---------- TIR length is measured, not clipped to the scan limit ----------
def test_tir_longer_than_the_old_limit_is_measured():
    for tir_len in (100, 200, 400, 900):
        r = random.Random(tir_len)
        tir = "".join(r.choices("ATGC", k=tir_len))
        seq = tir + "".join(r.choices("ATGC", k=2000)) + reverse_complement(tir)
        ev = structural.find_tir(seq)
        assert ev is not None, f"{tir_len} bp TIR not detected"
        assert abs(ev["tir_len"] - tir_len) <= 5, (tir_len, ev["tir_len"])


def test_short_canonical_tirs_still_found():
    for tir_len in (10, 11, 28):
        r = random.Random(500 + tir_len)
        tir = "".join(r.choices("ATGC", k=tir_len))
        ev = structural.find_tir(tir + "".join(r.choices("ATGC", k=1500)) + reverse_complement(tir))
        assert ev is not None and ev["tir_len"] >= tir_len


# ---------- a length limited by the record is declared, never presented as measured ----------
def test_repeat_running_off_the_record_is_flagged_as_a_lower_bound():
    r = random.Random(77)
    ltr = "".join(r.choices("ATGC", k=1000))
    seq = "".join(r.choices("ATGC", k=150)) + ltr + "".join(r.choices("ATGC", k=3000)) + ltr
    ev = structural.find_ltr(seq)                     # record stops exactly at the 3' LTR end
    assert ev is not None
    assert ev.get("length_is_lower_bound") is True
    assert "lower bound" in ev.get("bound_reason", "")


def test_complete_element_is_not_flagged():
    ev = structural.find_ltr(_ltr_element(900, seed=4))
    assert ev is not None
    assert ev.get("length_is_lower_bound") is None


# ---------- low-complexity sequence is not terminal-repeat evidence ----------
def test_simple_repeats_are_not_called_terminal_repeats():
    for seq in ("AT" * 3000, "TA" * 40, "CAG" * 2000, "GATC" * 1500, "A" * 5000):
        assert structural.find_ltr(seq) is None, seq[:12]
        assert structural.find_tir(seq) is None, seq[:12]     # (AT)n is its own reverse complement


# ---------- negative control: the detectors must stay quiet on random sequence ----------
def test_false_positive_rate_on_random_sequence_stays_low():
    ltr_hits = tir_hits = 0
    trials = 60
    for i in range(trials):
        seq = _rnd(6000, seed=90000 + i)
        if structural.find_ltr(seq):
            ltr_hits += 1
        if structural.find_tir(seq):
            tir_hits += 1
    assert ltr_hits == 0, f"{ltr_hits}/{trials} random sequences called as LTR pairs"
    assert tir_hits <= trials * 0.05, f"{tir_hits}/{trials} random sequences called as TIRs"


def test_subthreshold_terminal_repeat_is_reported_but_never_credited():
    """A pair rejected on identity is REPORTED evidence, never CREDITED evidence.

    The advisory near-miss row was first named "LTR candidate below identity threshold". Four consumers
    dispatch on a "LTR" prefix — classify.has_ltr, retroviral.py's ltr lookup, figures.py's band colour and
    main.py's evidence roll-up — so the rejected pair was read as a confirmed terminal repeat and came back
    as LTR/Copia at HIGH confidence. That is strictly worse than the silent discard the row replaced: the
    detector said "not an LTR" and the classifier answered "LTR, high confidence".

    Two independent guards, both asserted here: the type must not begin with a credited prefix, and an
    `advisory` row must be excluded even if it does."""
    from teagle_core import classify

    def dom(code, nt):
        return {"domain": code, "nt": list(nt), "strand": "+", "score": 90.0, "aa": [1, 100], "class": "retro"}

    doms = [dom("INT", (400, 700)), dom("RT", (800, 1400))]
    near = {"type": "Sub-threshold terminal direct repeat (advisory)", "identity": 78.4, "advisory": True,
            "five_prime": [0, 300], "three_prime": [1500, 1800], "element_span": [0, 1800]}

    # guard 1: the shipped name must not collide with any credited prefix
    for prefix in ("LTR", "TIR", "TSD", "poly"):
        assert not near["type"].startswith(prefix), f"advisory type collides with the {prefix!r} prefix"

    cl = classify.classify([near], doms)
    assert cl["te_class"] != "LTR/Copia", cl

    # guard 2: the advisory row must change NOTHING. This is stricter than asserting the absence of an
    # "LTR" substring, and it stays correct now that RT + a DDE integrase name the LTR *order* on domain
    # architecture alone: that call is reached with or without the rejected pair, so the pair is not what
    # produced it. An earlier form of this test asserted `"LTR" not in te_class`, which conflated "the
    # advisory repeat was credited" with "an LTR order was called by any route" and would have blocked a
    # correct domain-based call.
    without = classify.classify([], doms)
    assert cl["te_class"] == without["te_class"], (cl["te_class"], without["te_class"])
    assert cl["superfamily"] == without["superfamily"], (cl["superfamily"], without["superfamily"])
    assert cl["confidence"] == without["confidence"], (cl["confidence"], without["confidence"])

    # guard 3: the advisory flag blocks it even under a colliding name
    legacy = dict(near, type="LTR candidate below identity threshold")
    assert classify.classify([legacy], doms)["te_class"] != "LTR/Copia"
    assert classify.classify([legacy], doms)["te_class"] == without["te_class"]

    # and a genuine accepted LTR is unaffected
    real = {"type": "LTR (terminal direct repeat)", "five_prime": [0, 300],
            "three_prime": [1500, 1800], "element_span": [0, 1800]}
    good = classify.classify([real], doms)
    assert good["te_class"] == "LTR/Copia" and good["confidence"] == "High", good


def test_standalone_tir_floor_clears_chance_and_the_solo_ltr_case():
    """MIN_TIR_STANDALONE must stay above what chance produces, and above the solo-LTR case it guards.

    The floor decides whether a terminal inverted repeat may name a Class II element with NO transposase
    to corroborate it. Both numbers behind it are re-derived here rather than restated, so a change to
    find_tir that lengthens chance hits, or that changes what the copia LTR yields, fails the suite instead
    of silently loosening the gate.

    Kept deliberately cheap: 400 random sequences, enough to catch a floor that has fallen into the bulk of
    the chance distribution (13-15 bp), not enough to resolve the extreme tail. The 4,500-trial measurement
    that set the constant is recorded in structural.py.
    """
    import random
    from teagle_core import structural

    rng = random.Random(31337)
    longest_chance = 0
    for i in range(400):
        seq = "".join(rng.choice("ACGT") for _ in range(4000))
        hit = structural.find_tir(seq)
        if hit:
            longest_chance = max(longest_chance, hit["tir_len"])

    assert longest_chance < structural.MIN_TIR_STANDALONE, (
        f"chance inverted repeats reached {longest_chance} bp, at or above the standalone floor of "
        f"{structural.MIN_TIR_STANDALONE} bp -- a bare repeat that long would name a DNA transposon on noise")

    # The case the conservative branch was written for: a solo copia 5' LTR trips the TIR scan and was once
    # filed as a DNA transposon. It must stay below the floor. Built here as a direct repeat pair whose ends
    # happen to be weakly inverted, which is what a real LTR does.
    rng2 = random.Random(4242)
    ltr = "TGTTGGAATATAC" + "".join(rng2.choice("ACGT") for _ in range(250)) + "GTATATTCCAACA"
    hit = structural.find_tir(ltr)
    if hit:
        assert hit["tir_len"] < structural.MIN_TIR_STANDALONE, (
            f"a solo LTR yields a {hit['tir_len']} bp inverted repeat, at or above the standalone floor -- "
            f"a Class I fragment would be called Class II")
