"""Regression tests for the terminal-repeat measurement defects (2026-07-27).

Before this: find_ltr capped its seeding window at 1800 bp and reported whatever fraction of the
repeat fell inside both windows as a measured length (a 2.5 kb LTR read as 700 bp, a 4.4 kb LTR was
not detected at all); _terminal_tir stopped at 60 bp, so 100/200/400 bp TIRs all read as exactly 60;
and the anchor-density floor required ~93% identity, which silently dropped older LTR pairs and made
the 80% identity gate unreachable. None of it carried a flag — the truncated number was presented as
evidence and propagated into element_span, TSD/PBS/PPT, primer targets and the sealed manifest.
"""
import random
import time

import pytest

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


def _mutated(unit, n, rate, seed):
    r = random.Random(seed)
    s = unit * n
    return "".join(r.choice("ACGT") if r.random() < rate else c for c in s)


def _terminal_scans(seq):
    t0 = time.perf_counter()
    structural.find_ltr(seq), structural.find_tir(seq), structural.detect_all(seq)
    return time.perf_counter() - t0


@pytest.mark.parametrize("name,seq", [
    ("(AT)n", "AT" * 3000), ("(CAG)n", "CAG" * 2000), ("(GATC)n", "GATC" * 1500), ("A5000", "A" * 5000),
    ("(AT)n 5% diverged", _mutated("AT", 3000, 0.05, 1)), ("(CAG)n 5% diverged", _mutated("CAG", 2000, 0.05, 2)),
    ("random 6 kb", _rnd(6000, seed=77))], ids=lambda v: v if len(v) < 30 else "seq")
def test_low_complexity_and_random_records_are_scanned_quickly(name, seq):
    """Seeds from a low-complexity word pair with each other thousands of times; the seed search must mask them
    (4.1 first took 20-45 s on (AT)n and (GATC)n, 4.0.0 milliseconds)."""
    assert _terminal_scans(seq) < 2.0, name


def test_dust_mask_flags_microsatellites_and_leaves_random_sequence_alone():
    assert all(structural._dust("AT" * 100)) and all(structural._dust(_mutated("CAG", 100, 0.05, 3)))
    assert sum(structural._dust(_rnd(6000, seed=5))) == 0
    flank = structural._dust(_rnd(300, seed=6) + "A" * 80 + _rnd(300, seed=7))
    assert all(flank[300:380]) and not any(flank[:230]) and not any(flank[450:])


def _diverged(unit, n, rate, rng):
    return "".join(rng.choice("ACGT".replace(c, "")) if rng.random() < rate else c for c in unit * n)


def _tir_rows(seq):
    return [r for r in [structural.find_tir(seq)] + [e for e in structural.detect_all(seq) if e["type"].startswith("TIR")] if r]


@pytest.mark.parametrize("unit,n", [("AT", 3000), ("CAG", 2000)])
def test_diverged_microsatellite_arrays_yield_no_terminal_inverted_repeat(unit, n):
    """A short-period tandem repeat pairs with its own reverse complement, so a record that is nothing else holds no TIR."""
    for seed in range(8):
        for rate in (0.02, 0.05):
            seq = _diverged(unit, n, rate, random.Random(seed * 1000 + int(rate * 100)))
            assert not _tir_rows(seq), (unit, rate, seed, _tir_rows(seq)[0]["tir_len"])


def test_specimens_between_diverged_microsatellite_flanks_carry_no_arm_inside_a_flank():
    """The specimens' own ends lie 2 kb inside the record, beyond the TIR seed window: a row whose arm is in a flank is spurious."""
    from helpers import fixture_seq
    for acc in ("M11240", "M12927", "M80343", "X01005", "X05424"):
        for seed in range(3):
            for rate in (0.02, 0.05):
                rng = random.Random(seed * 7 + int(rate * 100))
                seq = _diverged("AT", 1000, rate, rng) + fixture_seq(acc) + _diverged("AT", 1000, rate, rng)
                for r in _tir_rows(seq):
                    assert r["five_prime"][1] > 2000 and r["three_prime"][0] < len(seq) - 2000, (acc, rate, seed, r["five_prime"])


@pytest.mark.parametrize("arm", ["ATATATATATATATATATATATATA", "TATATATATATATATATATATATATATATCTAT", "ATATATATATATATATATATATATATATGT"])
def test_arms_that_are_a_microsatellite_with_a_few_substitutions_are_not_a_terminal_repeat(arm):
    """_is_simple passes a 35-bp (TA)n arm that carries one C."""
    assert not structural.find_tir(arm + _rnd(2000, seed=3) + reverse_complement(arm))


def test_a_tsd_inside_a_microsatellite_run_does_not_corroborate_a_short_terminal_repeat(monkeypatch):
    """11-bp arms between 40-bp poly-A runs share a 12-bp 'TSD' that the run alone supplies (chance arms stay uncorroborated)."""
    def record(seed):
        r = random.Random(seed)
        arm = "".join(r.choices("ACGT", k=11))
        return _rnd(100, seed) + "A" * 40 + arm + _rnd(2000, seed + 1) + reverse_complement(arm) + "A" * 40 + _rnd(100, seed + 2)
    seqs = [record(s) for s in range(30)]
    corroborated = lambda: [s for s in seqs if (structural.find_tir(s) or {}).get("boundary_tests")]
    assert not corroborated()
    monkeypatch.setattr(structural, "_tsd_in_microsatellite", lambda *a: False)
    assert corroborated(), "the control no longer exercises the guard"



def test_a_long_random_ltr_is_not_rejected_as_a_simple_repeat():
    """The distinct-4-mer floor scaled with length past the 256 possible 4-mers, rejecting every repeat over 5.1 kb."""
    ltr_len = 5500
    ev = structural.find_ltr(_ltr_element(ltr_len, inner=1500, flank=200, seed=11))
    assert ev is not None and abs(ev["ltr_len"] - ltr_len) <= 5
    assert not structural._is_simple(_rnd(ltr_len, seed=12))
    assert structural._is_simple("ATTA" * 2000) and structural._is_simple("AT" * 3000)


def _tc1_in_flanks(runlen, rate, seed, tsd="TA"):
    """X01005 (Tc1, 54-bp arms) inserted at a TA site: random flank, (TA)n run, TA TSD, element, TA TSD, run, flank."""
    from helpers import fixture_seq
    core = fixture_seq("X01005").upper()
    r = random.Random(seed)
    flank = lambda: "".join(r.choices("ACGT", k=r.randint(100, 400)))
    run = lambda: _diverged("TA", runlen // 2, rate, r)
    left = flank() + run() + tsd
    return left + core + tsd + run() + flank(), len(left), len(left) + len(core)


def _tc1_arm(seq, core_start):
    """(arm length after trimming the palindromic TA TSD, whether that TSD was found, offset of the arm start from the element start)."""
    row = structural.find_tir(seq)
    if not row:
        return None
    hit = structural.palindromic_target_site(seq, row, "TA")
    row2, tsd = hit if hit else (row, None)
    return row2["tir_len"], bool(tsd), row2["five_prime"][0] - core_start


def test_a_tc1_next_to_a_ta_run_keeps_its_own_terminal_repeat():
    """The arms ran through the TA site into the (TA)n run, and the microsatellite guard then discarded the whole candidate
    (27 of 30 records dropped to Candidate for a run of 60 bp or more); the run is trimmed off instead."""
    cells = [(n, d) for n in (20, 40, 60, 100, 200) for d in (0.0, 0.05)]
    seeds = range(6)
    got = {(n, d): [_tc1_arm(*_tc1_in_flanks(n, d, 100 * n + s + int(d * 100))[:2]) for s in seeds] for n, d in cells}
    free = [_tc1_arm(*_tc1_in_flanks(0, 0.0, s)[:2]) for s in seeds]
    rate = lambda xs: sum(bool(x) and x[0] == 54 and x[2] == 0 for x in xs) / len(xs)
    assert rate(free) >= 0.8, free
    for (n, d), xs in got.items():
        assert all(xs), (n, d, xs)
        assert all(abs(x[2]) <= 8 and x[0] <= 62 for x in xs), (n, d, xs)     # no arm through the run (they were 70-250 bp)
        assert rate(xs) >= (1.0 if d == 0 else 0.6), (n, d, xs)               # a substitution at the run's edge can leave a few of its bases in the arm
        assert sum(x[1] for x in xs) >= len(xs) - 1, (n, d, xs)               # the TA TSD is still reported


def _headed_element(head, arm_len, tsd, seed):
    """An element whose own arms begin with a simple run (head + random), between random flanks with a TSD."""
    r = random.Random(seed)
    rnd = lambda n: "".join(r.choices("ACGT", k=n))
    arm = head + rnd(arm_len - len(head))
    el = arm + rnd(1500) + reverse_complement(arm)
    tsd = tsd or rnd(8)
    left = rnd(r.randint(150, 500)) + tsd
    return left + el + tsd + rnd(r.randint(150, 500)), len(left), arm_len


@pytest.mark.parametrize("head", ["A" * 10, "CA" * 8, "CAG" * 4])
def test_a_terminal_repeat_that_begins_with_a_simple_run_keeps_its_full_length(head):
    """The run trim cut a real arm's own head between random flanks (54-bp arms with an A10 or (CA)8 head reported as
    38-46 bp, 0 of 8 at full length). Only a run whose unit is its own reverse complement can pair across an insertion
    that split it; these heads cannot."""
    for tsd in (None, "TA"):
        for arm_len in (30, 54):
            exact = inside = 0
            for seed in range(8):
                seq, start, L = _headed_element(head, arm_len, tsd, 7000 + seed * 37 + arm_len + len(head))
                row = structural.find_tir(seq)
                f = row["five_prime"][0] if row else None
                exact += bool(row) and start - 4 <= f <= start + 1 and row["tir_len"] >= L - 1
                inside += bool(row) and f > start + 1
            assert exact >= 6 and inside == 0, (head, tsd, arm_len, exact, inside)


def test_a_terminal_repeat_made_of_a_unit_that_cannot_pair_with_itself_is_kept():
    """Ginger1 (Dfam) has a (TCTA)n terminal inverted repeat. (TCTA)n cannot pair with its own reverse complement, so an
    array of it yields no chance inverted repeat, and its arms are evidence; the period-4-6 guard rejected them."""
    found = 0
    for seed in range(8):
        seq, start, L = _headed_element("TCTA" * 12, 63, None, 9100 + seed)
        row = structural.find_tir(seq)
        found += bool(row) and abs(row["five_prime"][0] - start) <= 2 and row["tir_len"] >= L - 2
    assert found >= 7, found
    assert not structural._period_run("TCTA" * 10) and structural._period_run("GATC" * 10) and structural._period_run("CAG" * 10)


def test_only_a_run_whose_unit_is_its_own_reverse_complement_is_trimmed():
    for w, p, paired in (("TA" * 8, 2, True), ("GATC" * 4, 4, True), ("AATT" * 4, 4, True), ("AAATTT" * 3, 6, True),
                         ("A" * 12, 1, False), ("CA" * 8, 2, False), ("CAG" * 5, 3, False), ("AAT" * 5, 3, False)):
        assert structural._self_paired(w, p) is paired, w
    assert structural._self_paired("TATATACATATATA", 2)                     # a substitution does not change the unit


@pytest.mark.parametrize("unit", ["AATT", "GATC", "ACGT", "AAATTT"])
def test_palindromic_period_4_to_6_arrays_yield_no_terminal_inverted_repeat(unit):
    """(AATT)n, (GATC)n, (ACGT)n and (AAATTT)n are their own reverse complements; a diverged copy was reported as a TIR whose
    arms ran to the 1-kb limit (the guard tested periods 1-3)."""
    for seed in range(4):
        for rate in (0.02, 0.05, 0.10):
            seq = _diverged(unit, 6000 // len(unit), rate, random.Random(555 + seed * 131 + int(rate * 1000) + len(unit)))
            assert not _tir_rows(seq), (unit, rate, seed, _tir_rows(seq)[0]["tir_len"])


@pytest.mark.parametrize("unit", ["AATT", "GATC"])
def test_specimens_between_diverged_palindromic_flanks_carry_no_arm_inside_a_flank(unit):
    from helpers import fixture_seq
    for acc in ("M80343", "X01005", "X05424"):
        for seed in range(2):
            rng = random.Random(seed * 11 + len(unit))
            seq = _diverged(unit, 500, 0.03, rng) + fixture_seq(acc) + _diverged(unit, 500, 0.03, rng)
            for r in _tir_rows(seq):
                assert r["five_prime"][1] > 2000 and r["three_prime"][0] < len(seq) - 2000, (unit, acc, seed, r["five_prime"])


def test_short_arm_at_the_edge_of_a_record_cut_from_a_microsatellite_array_is_not_a_terminal_repeat():
    """A record cut from a (CAG)n array has no outside sequence at its ends: the 10-bp arm CAGCTGCAGC of a 20%-diverged
    array was never tested for period 3 (the window was under 20 bp); the 16 bp inside each arm are now read with it."""
    rows = []
    for seed in range(10):
        seq = _diverged("CAG", 2000, 0.2, random.Random(555 + seed * 131 + 200 + 3))
        rows += _tir_rows(seq)
    assert not rows, [r["tir_len"] for r in rows]
    seq = _diverged("CAG", 2000, 0.2, random.Random(90210 + 200 + 3))
    assert not _tir_rows(seq)


def test_the_tandem_run_test_rarely_rejects_a_random_window():
    """Random 10-40 bp windows read as a run of period 1-6: 0.04% at 10 bp and 0.08% at 20 bp at GC 0.5 (100,000 windows per cell)."""
    r = random.Random(41)
    for length, ceiling in ((10, 0.005), (20, 0.005), (26, 0.005), (40, 0.001)):
        n = 20000
        hits = sum(structural._period_run("".join(r.choices("ACGT", k=length))) for _ in range(n))
        assert hits / n < ceiling, (length, hits / n)


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
