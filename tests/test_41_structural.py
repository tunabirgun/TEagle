"""4.1 structural detectors on realistic inputs: genomic flanks, LTR copies with indels, tails that are
not at the record end, reverse-oriented elements. Every case builds on the bundled specimens so that a
regression against the 4.0.0 anchors (copia LTR 276, gypsy 482, Tc1 TIR 54, Ac TIR 11, L1.2 tail 27) fails here.

Success rates over seeded random flanks are asserted against floors and printed, because some ambiguity is
intrinsic: when one end of a repeat overshoots by a base, the TSD that still matches there is one base shorter,
and the neighbouring boundary that restores the full length is explained just as well by a chance base."""
import random
import time

import pytest

from helpers import fixture_seq
from teagle_core import sequtil, structural
from teagle_core.sequtil import reverse_complement

SEEDS = range(60)
FLANK = 300
TSD_LEN = {"M11240": 5, "M12927": 5, "X05424": 8}


def _rnd(rng, n):
    return "".join(rng.choices("ACGT", k=n))


def embed(acc, seed, tsd=None, flank=FLANK, element=None):
    """flank + TSD + element + TSD + flank; returns the record and the true element span."""
    rng = random.Random(seed * 7919 + len(acc))
    el = element if element is not None else fixture_seq(acc)
    t = tsd or ""
    left, right = _rnd(rng, flank), _rnd(rng, flank)
    s = len(left) + len(t)
    return left + t + el + t + right, s, s + len(el)


def planted_tsd(acc, seed):
    if acc not in TSD_LEN:
        return "TA"
    return _rnd(random.Random(seed * 31 + 5), TSD_LEN[acc])


def _row(ev, prefix):
    return next((e for e in ev if e["type"].startswith(prefix)), None)


def _flanked_success(acc, kind):
    ok = []
    for seed in SEEDS:
        seq, s, e = embed(acc, seed, planted_tsd(acc, seed))
        ev = structural.detect_all(seq)
        r, t = _row(ev, kind), _row(ev, "TSD")
        good = bool(r and t) and r["element_span"] == [s, e] and t["upstream"][1] == s and t["downstream"][0] == e
        if good:    # the length the superfamily would pick (classify re-runs find_tsd with `expect`) is the planted one
            good = structural.find_tsd(seq, s, e, expect=TSD_LEN[acc])["length"] == TSD_LEN[acc]
        ok.append(good)
    return sum(ok) / len(ok)


# ---------- regression anchors on the bare specimens ----------
@pytest.mark.parametrize("acc,kind,key,length,ident", [
    ("M11240", "LTR", "ltr_len", 276, 100.0), ("M12927", "LTR", "ltr_len", 482, 100.0),
    ("X01005", "TIR", "tir_len", 54, 100.0), ("X05424", "TIR", "tir_len", 11, 90.9)])
def test_bare_specimen_anchors_do_not_move(acc, kind, key, length, ident):
    r = _row(structural.detect_all(fixture_seq(acc)), kind)
    assert r[key] == length and r["identity"] == ident
    assert r["five_prime"][0] == 0 and r["element_span"][1] == len(fixture_seq(acc))


def test_bare_l1_tail_anchor():
    s = fixture_seq("M80343")
    r = _row(structural.detect_all(s), "poly-A")
    assert r["length"] == 27 and r["pos"] == [len(s) - 27, len(s)]


def test_bare_specimens_stay_fast():
    for acc in ("M11240", "M12927", "M80343", "X01005", "X05424"):
        s = fixture_seq(acc)
        t = []
        for _ in range(5):
            t0 = time.perf_counter()
            structural.detect_all(s)
            t.append(time.perf_counter() - t0)
        assert sorted(t)[2] < 0.05, (acc, sorted(t))     # 4.0.0 took 0.3-1.2 ms; a generous ceiling, not a benchmark


# ---------- flanked records: boundaries and TSD ----------
@pytest.mark.parametrize("acc,kind,floor", [("M11240", "LTR", 0.9), ("M12927", "LTR", 0.75), ("X05424", "TIR", 0.85)])
def test_flanked_element_boundaries_and_tsd(acc, kind, floor):
    rate = _flanked_success(acc, kind)
    print(f"flanked {acc}: {rate:.2f} over {len(SEEDS)} seeds (4.0.0 measured 0.57 / 0.42 / 0.00)")
    assert rate >= floor, rate


def test_flanked_tc1_target_site_recovered_for_the_superfamily_that_owns_it():
    """TA is its own reverse complement, so the arms run through the target site (TIR 56-58 rather than 54).
    Only the Tc1/Mariner caller may trim it: for any other TIR, trimming palindromic outer bases would
    manufacture a TSD."""
    ok = []
    for seed in SEEDS:
        seq, s, e = embed("X01005", seed, "TA")
        tir = _row(structural.detect_all(seq), "TIR")
        res = structural.palindromic_target_site(seq, tir)
        ok.append(bool(res) and res[0]["element_span"] == [s, e] and res[0]["tir_len"] == 54
                  and res[1]["motif"] == "TA" and res[1]["upstream"] == [s - 2, s])
    print(f"flanked Tc1 target site: {sum(ok) / len(ok):.2f} over {len(SEEDS)} seeds")
    assert sum(ok) / len(ok) >= 0.9


def test_palindromic_target_site_is_not_invented_for_other_elements():
    """A TIR whose outer pair is not the superfamily's palindromic motif gets no TSD."""
    seq, s, e = embed("X05424", 3, planted_tsd("X05424", 3))
    tir = _row(structural.detect_all(seq), "TIR")
    assert structural.palindromic_target_site(seq, tir, motif="TA") is None


def test_no_tsd_control_does_not_fabricate_long_target_sites():
    """Same flanked elements, nothing planted. Joint boundary search raises the chance of a TSD, so the rate
    of fabricated >= 4 bp TSDs is measured, not assumed: 4.0.0 gave 0/100 (copia), 2/100 (gypsy), 1/100 (Tc1)."""
    for acc, kind, ceiling in (("M11240", "LTR", 0.07), ("M12927", "LTR", 0.08), ("X01005", "TIR", 0.05), ("X05424", "TIR", 0.03)):
        n = 0
        for seed in range(100):
            seq, s, e = embed(acc, seed)
            t = _row(structural.detect_all(seq), "TSD")
            n += bool(t and t["length"] >= 4)
        print(f"no-TSD control {acc}: {n}/100 TSD rows of >= 4 bp")
        assert n / 100 <= ceiling, (acc, n)


# ---------- LTR copies with indels ----------
def _del(seq, pos, k):
    return seq[:pos] + seq[pos + k:]


@pytest.mark.parametrize("k,ident", [(1, 99.6), (2, 99.3), (3, 98.9), (10, 96.4)])
def test_single_indel_in_a_ltr_copy_reads_as_one_indel_not_divergence(k, ident):
    cop = fixture_seq("M11240")
    q = _del(cop, len(cop) - 276 + 138, k)
    ev = structural.find_ltr(q)
    assert ev and ev["ltr_len"] == 276 and ev["identity"] == ident and ev["indels"] == 1
    assert ev["five_prime"] == [0, 276]
    seq, s, e = embed("M11240", 3, "GATTC", element=q)
    ev = structural.find_ltr(seq)
    assert ev["five_prime"] == [s, s + 276] and ev["three_prime"] == [e - (276 - k), e]
    assert "length_is_lower_bound" not in ev            # the repeat ended inside the record, so it was measured


def test_ltr_with_an_indel_is_classified_as_ltr_not_advisory():
    cop = fixture_seq("M11240")
    ev = structural.detect_all(_del(cop, len(cop) - 276 + 138, 1))
    assert _row(ev, "LTR (terminal")["identity"] >= 99
    assert not any(e.get("advisory") for e in ev)
    assert _row(ev, "TIR") is None                      # no spurious TIR row from the out-of-phase half


def test_diverged_ltr_pairs_with_several_indels_keep_their_full_length():
    def mutate(r, s):
        a = list(s)
        for i in range(len(a)):
            if r.random() < 0.04:
                a[i] = r.choice("ACGT")
        for _ in range(3):
            p, k = r.randrange(20, len(a) - 20), r.choice((1, 1, 2, 3, 5))
            if r.random() < 0.5:
                del a[p:p + k]
            else:
                a[p:p] = [r.choice("ACGT") for _ in range(k)]
        return "".join(a)
    full = 0
    for seed in range(40):
        r = random.Random(seed)
        ltr = _rnd(r, 500)
        seq = _rnd(r, 200) + ltr + _rnd(r, 3000) + mutate(r, ltr) + _rnd(r, 200)
        ev = structural.find_ltr(seq)
        full += bool(ev) and abs(ev["ltr_len"] - 500) <= 12 and ev["five_prime"][0] >= 195
    print(f"LTR pairs, 4% substitutions + 3 indels: {full}/40 measured at full length (4.0.0: 6/100)")
    assert full >= 36


def test_lower_bound_names_the_actual_cause():
    cop = fixture_seq("M11240")
    ev = structural.find_ltr(_del(cop, len(cop) - 276 + 138, 1))        # bare record: the 3' copy ends at the record end
    assert ev["length_is_lower_bound"] is True and "record ends" in ev["bound_reason"]


def test_ambiguous_bases_are_neither_matches_nor_mismatches():
    cop = fixture_seq("M11240")
    q = cop[:len(cop) - 276 + 100] + "N" * 40 + cop[len(cop) - 276 + 140:]
    ev = structural.find_ltr(q)
    assert ev["identity"] == 100.0 and ev["ltr_len"] == 276 and ev["ambiguous_columns"] == 40


def test_ambiguity_codes_are_not_a_target_site_duplication():
    seq = "ACGTA" + "RYW" + "C" * 60 + "RYW" + "GATTA"
    assert structural.find_tsd(seq, 8, 68, min_tsd=3) is None
    seq = "ACGTA" + "ACG" + "C" * 60 + "ACG" + "GATTA"
    assert structural.find_tsd(seq, 8, 68, min_tsd=3)["motif"] == "ACG"


# ---------- TIR search in flanked records ----------
def _tir_element(r, tir_len, mismatch_at=None, inner=1500):
    tir = _rnd(r, tir_len)
    three = list(reverse_complement(tir))
    if mismatch_at is not None:
        three[mismatch_at] = {"A": "C", "C": "A", "G": "T", "T": "G"}[three[mismatch_at]]
    return tir + _rnd(r, inner) + "".join(three)


def test_imperfect_tir_with_a_central_mismatch_is_found_when_a_tsd_corroborates_it():
    found = 0
    for seed in range(40):
        r = random.Random(seed)
        el = _tir_element(r, 16, mismatch_at=8)
        tsd = _rnd(r, 8)
        seq = _rnd(r, 300) + tsd + el + tsd + _rnd(r, 300)
        s, e = 308, 308 + len(el)
        tir = structural.find_tir(seq)
        found += bool(tir) and tir["element_span"] == [s, e]
    print(f"16 bp TIR with a central mismatch, 8 bp TSD, 300 bp flanks: {found}/40")
    assert found >= 34


def test_weak_tir_without_a_tsd_is_not_called_in_flanked_records():
    """Without the TSD the same arms are indistinguishable from chance."""
    called = 0
    for seed in range(60):
        r = random.Random(seed)
        el = _tir_element(r, 16, mismatch_at=8)
        called += bool(structural.find_tir(_rnd(r, 300) + el + _rnd(r, 300)))
    assert called <= 6, called


def test_subterminal_inverted_repeat_does_not_outrank_the_corroborated_termini():
    """Bundled Ac carries a 13 bp internal inverted repeat that out-counts its 11 bp terminal one; the TSD
    decides, and the winner is reported as corroborated."""
    seq, s, e = embed("X05424", 11, planted_tsd("X05424", 11))
    tir = structural.find_tir(seq)
    assert tir["element_span"] == [s, e] and "boundary_uncertain" not in tir


def test_uncorroborated_arms_are_reported_as_uncertain_and_the_tsd_state_says_so():
    """Tc1 in random flanks with no target-site duplication planted: the arms are strong enough to stand but
    no TSD fixes the boundary, so the TIR row says so and the TSD state is 'not assessed', not 'absent'."""
    uncertain = 0
    for seed in range(20):
        seq, s, e = embed("X01005", seed)
        ev = structural.detect_all(seq)
        tir = _row(ev, "TIR")
        assert tir is not None
        if tir.get("boundary_uncertain"):
            uncertain += 1
            if _row(ev, "TSD") is None:               # a chance 2-3 bp match at the arms is still reported as found
                assert tir["tsd_state"].startswith("not assessed (element boundary uncertain")
    assert uncertain >= 15, uncertain


def test_tsd_state_distinguishes_missing_partial_and_absent_flanks():
    cop = fixture_seq("M11240")
    assert _row(structural.detect_all(cop), "LTR")["tsd_state"].startswith("not assessable")
    seq = "GAT" + cop + "TTG"                             # 3 bp each side: a 5 bp TSD could not be tested
    st = _row(structural.detect_all(seq), "LTR")["tsd_state"]
    assert st.startswith("partially assessable") and "3 bp" in st
    seq, s, e = embed("M11240", 5)
    assert _row(structural.detect_all(seq), "LTR")["tsd_state"].split(" ")[0] in ("absent", "found")


LIVE_CHANCE_TIR_SEED = 41


def live_chance_tir(seed=LIVE_CHANCE_TIR_SEED, per_cell=300, gcs=(0.5, 0.35), lengths=(1500, 4000, 6000)):
    """The live chance-TIR measurement: find_tir on seeded random sequences. Returns (longest bp, hits, n).
    This, not a recorded benchmark file, is what MIN_TIR_STANDALONE and the classify sentence are checked against."""
    rng = random.Random(seed)
    top = hits = n = 0
    for gc in gcs:
        weights = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
        for L in lengths:
            for _ in range(per_cell):
                tir = structural.find_tir("".join(rng.choices("ACGT", weights, k=L)))
                n += 1
                if tir:
                    hits += 1
                    top = max(top, tir["tir_len"])
    return top, hits, n


def test_chance_tir_ceiling_with_the_flanked_search_stays_below_the_standalone_floor():
    """The flanked-record search must not lift what chance produces above MIN_TIR_STANDALONE. 1,800 random
    sequences, uniform and AT-rich, checked for the ceiling and for the false-positive rate."""
    top, hits, n = live_chance_tir()
    print(f"chance TIRs: {hits}/{n}, longest {top} bp (floor {structural.MIN_TIR_STANDALONE})")
    assert top < structural.MIN_TIR_STANDALONE
    assert hits / n <= 0.05


def test_fast_live_chance_tir_check_at_the_recorded_seed():
    """A few hundred sequences at a recorded seed: the longest chance TIR stays below the standalone floor with
    margin enough to notice a find_tir change that moves the bulk of the distribution."""
    top, hits, n = live_chance_tir(seed=LIVE_CHANCE_TIR_SEED, per_cell=75)
    print(f"seed {LIVE_CHANCE_TIR_SEED}: {hits}/{n} chance TIRs, longest {top} bp")
    assert n == 450 and top < structural.MIN_TIR_STANDALONE


def test_tir_lengths_are_not_extended_by_chance_when_the_arms_have_no_corroboration():
    """Arms placed by repeat strength alone are reported as the exact core that qualified them, so chance
    cannot lengthen them past what a perfect run gives (see structural.find_tir)."""
    rng = random.Random(41)
    uncertain = []
    for _ in range(900):
        tir = structural.find_tir("".join(rng.choices("ACGT", k=4000)))
        if tir and tir.get("boundary_uncertain"):
            uncertain.append(tir)
    assert uncertain, "the construct must produce uncorroborated hits for this test to mean anything"
    assert all(t["identity"] == 100.0 and t["tir_len"] < structural.MIN_TIR_STANDALONE for t in uncertain)


# ---------- tails that are not at the record end ----------
def test_flanked_l1_tail_and_its_target_site_duplication_are_found():
    el = fixture_seq("M80343")
    ok = []
    for seed in SEEDS:
        rng = random.Random(seed)
        tsd, left, right = _rnd(rng, 12), _rnd(rng, 300), _rnd(rng, 300)
        seq = left + tsd + el + tsd + right
        s = 312
        e = s + len(el)
        ev = structural.detect_all(seq)
        a, t = _row(ev, "poly-A"), _row(ev, "TSD")
        ok.append(bool(a and t) and a["pos"] == [e - 27, e] and t["upstream"][1] == s and t["downstream"][0] == e
                  and a.get("tsd_corroborated") is True)
    print(f"flanked L1.2 tail + TSD: {sum(ok)}/{len(ok)}")
    assert sum(ok) / len(ok) >= 0.9


def test_uncorroborated_tract_inside_a_flank_is_not_a_tail():
    el = fixture_seq("M80343")
    for seed in range(40):
        rng = random.Random(seed)
        ev = structural.detect_all(_rnd(rng, 300) + el + _rnd(rng, 300))
        assert _row(ev, "poly-A") is None


def _chance_tail_rate(gc, plant):
    """Records with no TSD anywhere: random flanks, and optionally an internal A14 run (150-450 bp from the
    3' end) that is not at the record end. Returns how many of 400 seeded records get a tail credited."""
    weights = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
    rng = random.Random(4100 + int(gc * 100) + plant)
    credited = 0
    for _ in range(400):
        seq = "".join(rng.choices("ACGT", weights, k=3000))
        if plant:
            at = len(seq) - rng.randint(150, 450)
            seq = seq[:at] + "A" * 14 + seq[at + 14:]
        rows, tsds = structural._tails_both(seq, structural.POLYA_MIN_RUN)
        credited += bool(tsds or any(r.get("tsd_corroborated") for r in rows))
    return credited


@pytest.mark.parametrize("gc,plant", [(0.5, True), (0.35, True), (0.35, False)])
def test_chance_direct_repeat_does_not_corroborate_a_tail(gc, plant):
    """About 3,000 left-copy x end placements are tested per tract, so a bare 8-mer match is expected by chance
    in several percent of records. The tail-TSD acceptance is family-wise: control with nothing planted."""
    n = _chance_tail_rate(gc, plant)
    print(f"chance tail credit gc={gc} plant={plant}: {n}/400")
    assert n <= 4


def test_tail_tsd_row_carries_its_family_wise_chance():
    el = fixture_seq("M80343")
    rng = random.Random(3)
    tsd = _rnd(rng, 12)
    seq = _rnd(rng, 300) + tsd + el + tsd + _rnd(rng, 300)
    t = _row(structural.detect_all(seq), "TSD")
    single = structural._chance_p(seq, t["length"])
    assert t["n_tests"] > 1000
    assert t["chance_probability"] == pytest.approx(1 - (1 - single) ** t["n_tests"], rel=0.02)
    assert t["chance_level"] is (t["chance_probability"] > structural.CHANCE_LEVEL_P) is False


def test_tail_tolerates_an_interruption_and_keeps_its_length_otherwise():
    base = "ACGT" * 25
    assert structural.find_polya(base + "A" * 27)[0]["length"] == 27
    t = structural.find_polya(base + "A" * 13 + "G" + "A" * 13)[0]
    assert t["pos"] == [len(base), len(base) + 27]
    assert structural.find_polya(base + "AAAAAAAAGAC")[0]["length"] == 8      # slack: a short non-A end
    assert structural.find_polya(base) == []


def test_reverse_oriented_element_has_a_poly_t_tract_and_its_tsd():
    el = fixture_seq("M80343")
    rc = reverse_complement(el)
    rng = random.Random(5)
    tsd = _rnd(rng, 12)
    seq = _rnd(rng, 300) + tsd + rc + tsd + _rnd(rng, 300)
    ev = structural.detect_all(seq)
    t5 = _row(ev, "poly-T")
    assert t5["pos"] == [312, 339]
    assert _row(ev, "TSD")["upstream"][1] == 312


# ---------- cis-elements follow the element's orientation ----------
@pytest.mark.parametrize("acc", ["M11240", "M12927"])
def test_reverse_oriented_ltr_element_reports_mapped_cis_elements_on_the_right_strand(acc):
    fwd = fixture_seq(acc)
    n = len(fwd)
    f = structural.detect_all(fwd)
    r = structural.detect_all(reverse_complement(fwd))
    assert _row(f, "LTR")["cis_orientation"] == "+" and _row(r, "LTR")["cis_orientation"] == "-"
    for kind in ("PPT", "polyA-signal"):
        pf, pr = _row(f, kind), _row(r, kind)
        assert pr["pos"] == [n - pf["pos"][1], n - pf["pos"][0]]
        assert pr["motif"] == pf["motif"] and pr["strand"] == "-" and pf["strand"] == "+"
    ltr5 = _row(r, "LTR")["five_prime"]
    assert _row(r, "polyA-signal")["pos"][1] <= ltr5[1]    # the signal sits in the copy that is the element's 3' LTR


def test_caller_orientation_wins():
    cop = fixture_seq("M11240")
    forced = structural.detect_all(cop, orientation="-")
    assert _row(forced, "LTR")["cis_orientation_source"] == "caller"
    assert _row(forced, "PPT") is None                    # no polypurine tract reads on the minus strand of this element


def test_unknown_orientation_withholds_cis_elements_and_says_why():
    """A random LTR pair has no polypurine tract, confident PBS or poly(A) signal on either strand: the
    elements are not reported on a strand the evidence does not support, and the LTR row says so."""
    undetermined = 0
    for seed in range(30):
        r = random.Random(seed)
        ltr = _rnd(r, 400)
        ev = structural.detect_all(_rnd(r, 100) + ltr + _rnd(r, 1500) + ltr + _rnd(r, 100))
        row = _row(ev, "LTR (terminal")
        assert row is not None
        if row["cis_orientation"] == "undetermined":
            undetermined += 1
            assert "orientation of the element is undetermined" in row["cis_note"]
            assert not [e for e in ev if e["type"].startswith(("PBS", "PPT", "polyA"))]
    assert undetermined >= 20, undetermined


def test_pbs_rows_carry_their_chance_probability_and_the_panel_scope():
    r = random.Random(1)
    row = next(x for x in (structural.find_pbs(_rnd(r, 300), 50) for _ in range(200)) if x)
    assert 0 < row["chance_probability"] <= 1 and row["panel"] == ["tRNA-Lys3"]
    ltr = _row(structural.detect_all(fixture_seq("M11240")), "LTR")
    assert "tRNA-Lys3 only" in ltr["pbs_scope"] or "tRNA-Lys3" in ltr["pbs_scope"]


# ---------- TSD congruence and chance accounting ----------
def test_ltr_tsd_congruence_accepts_the_4_to_6_bp_range():
    for sf in ("Copia (Ty1)", "Gypsy (Ty3)"):
        for n, verdict in ((3, "incongruent"), (4, "congruent"), (5, "congruent"), (6, "congruent"), (7, "incongruent")):
            c = structural.tsd_congruence(n, sf)
            assert c["verdict"] == verdict and c["expected"] == 5 and c["expected_range"] == [4, 6]


def test_tc1_tsd_congruence_requires_the_ta_motif_when_a_motif_is_given():
    assert structural.tsd_congruence(2, "Tc1/Mariner")["verdict"] == "congruent"        # length only, as before
    assert structural.tsd_congruence(2, "Tc1/Mariner", motif="TA")["verdict"] == "congruent"
    c = structural.tsd_congruence(2, "Tc1/Mariner", motif="GC")
    assert c["verdict"] == "incongruent" and c["expected_motif"] == "TA"
    assert structural.tsd_congruence(2, "MULE (Mutator)")["verdict"] == "not diagnostic"


def test_tsd_chance_is_adjusted_for_the_lengths_and_boundaries_tested():
    seq = "ACGTAC" + "GATTACACTGCAGTCCTAGG" * 3 + "ACGTAC"
    one = structural.find_tsd(seq, 6, 66)
    many = structural.find_tsd(seq, 6, 66, placements=36)
    p = structural._chance_p(seq, 6)
    assert one["n_tests"] == structural.MAX_TSD - 6 + 1 and many["n_tests"] == 36 * one["n_tests"]
    assert one["chance_probability"] >= p and many["chance_probability"] > 10 * one["chance_probability"]
    assert many["chance_probability"] == pytest.approx(1 - (1 - one["chance_probability"]) ** 36, rel=0.02)
    assert many["chance_level"] is (many["chance_probability"] > structural.CHANCE_LEVEL_P)
    exp = structural.find_tsd(seq, 6, 66, expect=6)                  # one pre-specified length is one test
    assert exp["n_tests"] == 1 and exp["chance_probability"] == pytest.approx(p, rel=0.01)


def test_ltr_boundary_snap_records_the_shifts_searched_and_the_motif_it_used():
    seen = moved = bonus = 0
    for acc in ("M11240", "M12927"):
        for seed in range(100):
            seq, s, e = embed(acc, seed)
            ev = structural.detect_all(seq)
            ltr, tsd = _row(ev, "LTR"), _row(ev, "TSD")
            if not ltr:
                continue
            seen += 1
            if "boundary_adjustment" in ltr:
                moved += 1
                assert ltr["boundary_tests"] >= 36
                assert not tsd or tsd["n_tests"] >= 36
                if ltr["termini"].get("selected_for_motif"):
                    bonus += 1
                    assert ltr["termini"]["canonical"] is True
            else:
                assert "selected_for_motif" not in ltr["termini"]
    print(f"LTR snap: {moved}/{seen} moved, {bonus} used the canonical bonus")
    assert moved and bonus


def test_tsd_corroborated_tir_row_counts_its_boundary_placements():
    n = 0
    for seed in SEEDS:
        seq, s, e = embed("X01005", seed, planted_tsd("X01005", seed))
        ev = structural.detect_all(seq)
        tir, tsd = _row(ev, "TIR"), _row(ev, "TSD")
        if tir and tsd and not tir.get("boundary_uncertain") and tsd["length"] >= structural._TIR_WEAK_TSD:
            n += 1
            assert tsd["n_tests"] > 1 and tir["boundary_tests"] > 1
    assert n


def test_tsd_rows_carry_composition_adjusted_chance():
    seq = "ACGTAC" + "GATTACACTGCAGTCCTAGG" * 3 + "ACGTAC"
    t = structural.find_tsd(seq, 6, 66)
    assert t["length"] == 6 and t["chance_level"] is False and 0 < t["chance_probability"] < 1e-3
    short = structural.find_tsd("GC" + "ACGTTGCA" * 8 + "GC", 2, 66)
    assert short["length"] == 2 and short["chance_level"] is True and short["chance_probability"] > 1e-3
    at_rich = "AT" * 20 + "G" + "AT" * 20
    p_uniform = structural._chance_p("ACGT" * 50, 5)
    assert structural._chance_p(at_rich, 5) > p_uniform                                    # AT-rich flanks match more often


# ---------- named scan windows ----------
def test_scan_windows_are_named_and_reported():
    sc = structural.scan_scope("ACGT" * 5000)
    assert sc["ltr_seed_window_bp"] == structural._LTR_SEED_WINDOW
    assert sc["tir_seed_window_bp"] == structural._TIR_SEED_WINDOW
    assert sc["tail_window_bp"] == structural._TAIL_WINDOW and sc["covers_whole_record"] is False
    assert structural.scan_scope("ACGT" * 20)["covers_whole_record"] is True


# ---------- sequtil: text before the first header ----------
def test_text_before_the_first_header_is_reported_not_silently_dropped():
    assert sequtil.fasta_preamble("ACGTACGTAAAA\n>b\nCCCCGGGG") == "ACGTACGTAAAA"
    assert sequtil.parse_fasta("ACGTACGTAAAA\n>b\nCCCCGGGG") == [("b", "CCCCGGGG")]     # parsing itself is unchanged
    for blank in (">a\nACGT", "\n  \n>a\nACGT", "ACGT", "", None):
        assert sequtil.fasta_preamble(blank) == ""
    assert sequtil.fasta_preamble("AC GT\r\nGG\r\n>x\r\nAA") == "ACGTGG"
