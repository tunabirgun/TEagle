from teagle_core import primers, sequtil
from helpers import fixture_seq

RC = sequtil.reverse_complement


def test_design_primers_real_and_deterministic():
    seq = fixture_seq("M11240")                       # copia, 5146 bp
    p = {"prod_min": 150, "prod_max": 500, "min_tm": 57, "max_tm": 63}
    d1 = primers.design_primers(seq, p)
    d2 = primers.design_primers(seq, p)
    assert d1["candidates"], "Primer3 should return candidates for copia"
    assert d1 == d2, "Primer3 design must be deterministic for identical input"
    c = d1["candidates"][0]
    assert 55 <= c["left_tm"] <= 65 and 55 <= c["right_tm"] <= 65
    assert 150 <= c["product_size"] <= 500
    assert c["left_seq"] == seq[c["left_pos"][0]:c["left_pos"][1]]     # left coords map to the template
    assert c["right_seq"] == RC(seq[c["right_pos"][0]:c["right_pos"][1]])


def _amplicon_template(fwd, rev, mid_len=200, flank=50):
    mid = "A" * mid_len
    seq = "C" * flank + fwd + mid + RC(rev) + "C" * flank
    left = flank
    right = flank + len(fwd) + mid_len + len(rev)
    return seq, left, right


def test_in_silico_pcr_on_target():
    fwd, rev = "GACTGACTGTCAGTCAGGCT", "TTGGCCATTGGCACTGGCAT"
    seq, left, right = _amplicon_template(fwd, rev)
    amps = primers.in_silico_pcr(fwd, rev, seq, "t", max_mm=1, tp=5,
                                 prod_min=100, prod_max=600, target_span=[left, right])
    assert len(amps) == 1
    a = amps[0]
    assert a["start"] == left and a["end"] == right
    assert a["fwd_mm"] == 0 and a["rev_mm"] == 0 and a["on_target"] is True


def test_in_silico_pcr_three_prime_mismatch_rejected():
    fwd, rev = "GACTGACTGTCAGTCAGGCT", "TTGGCCATTGGCACTGGCAT"
    seq, left, right = _amplicon_template(fwd, rev)
    # put a 3'-terminal mismatch into the forward binding site on the template
    bad_last = "A" if fwd[-1] != "A" else "C"
    fwd_site = fwd[:-1] + bad_last
    seq2 = seq[:left] + fwd_site + seq[left + len(fwd):]
    amps = primers.in_silico_pcr(fwd, rev, seq2, "t", max_mm=2, tp=5,
                                 prod_min=100, prod_max=600)
    assert amps == [], "a 3'-terminal mismatch must abolish the amplicon under the strict 3' rule"


def test_in_silico_pcr_requires_inward_pair():
    # two same-orientation forward hits, no reverse-capable site -> no amplicon
    fwd = "GACTGACTGTCAGTCAGGCT"
    rev = "TTGGCCATTGGCACTGGCAT"          # its RC does not appear in the template
    seq = "C" * 40 + fwd + "A" * 150 + fwd + "C" * 40
    amps = primers.in_silico_pcr(fwd, rev, seq, "t", max_mm=1, tp=5, prod_min=50, prod_max=500)
    assert amps == [], "independent same-strand hits must not form an amplicon"


def test_base_ok_iupac_matching():
    assert primers._base_ok("R", "A") and primers._base_ok("R", "G")       # R = A/G
    assert not primers._base_ok("R", "C") and not primers._base_ok("R", "T")
    assert primers._base_ok("N", "A") and primers._base_ok("Y", "T") and primers._base_ok("A", "A")
    assert not primers._base_ok("A", "G")                                   # a concrete primer base stays strict
    assert not primers._base_ok("R", "N")                                   # an ambiguous TEMPLATE base is a mismatch (conservative)


def test_in_silico_pcr_matches_iupac_degenerate_primer():
    # a degenerate consensus primer (IUPAC at the 3' end) must still bind — the genome-scan path (isPcr) is
    # ambiguity-aware, so local in-silico PCR must agree rather than silently report the pair as non-binding
    concrete, rev = "GACTGACTGTCAGTCAGGCA", "TTGGCCATTGGCACTGGCAT"          # concrete fwd ends in A
    seq, left, right = _amplicon_template(concrete, rev)                    # template carries the literal A at the fwd 3' end
    degen = concrete[:-1] + "R"                                             # R (A/G) covers the template's A
    amps = primers.in_silico_pcr(degen, rev, seq, "t", max_mm=0, tp=5, prod_min=100, prod_max=600)
    pair = [a for a in amps if not a.get("single_primer")]
    assert len(pair) == 1 and pair[0]["start"] == left                      # degenerate 3' base binds -> amplicon found


def test_in_silico_pcr_degenerate_primer_stays_specific():
    # R (A/G) must NOT match a template carrying C at that 3' position -> the strict 3' rule still abolishes it
    concrete, rev = "GACTGACTGTCAGTCAGGCC", "TTGGCCATTGGCACTGGCAT"          # concrete fwd ends in C
    seq, _l, _r = _amplicon_template(concrete, rev)
    degen = concrete[:-1] + "R"                                             # R does not cover C
    amps = primers.in_silico_pcr(degen, rev, seq, "t", max_mm=0, tp=5, prod_min=100, prod_max=600)
    assert [a for a in amps if not a.get("single_primer")] == []


# ---- non-templated 5' tails ------------------------------------------------------------------------
# A primer carrying a restriction site, adapter, barcode or promoter at its 5' end does not anneal over its
# full length on the first cycle, yet the tail is copied into the product. Requiring a whole-primer match
# rejects that whole class of assay; anchoring at the 3' end accepts it and accounts for the tail.

def test_tailed_primer_binds_and_product_includes_the_tail():
    core_f, core_r = "GACTGACTGTCAGTCAGGCT", "TTGGCCATTGGCACTGGCAT"
    seq, left, right = _amplicon_template(core_f, core_r)          # template carries the cores only
    # BamHI- and HindIII-like 5' additions. Neither may end in a base that pairs with the flank, or the
    # junction extends the match by a base and the core is genuinely one longer than the tail implies:
    # the flank here is poly-C, so a forward tail must not end in C and a reverse tail must not end in G.
    tail_f, tail_r = "GGATCCGG", "AAGCTTAAA"
    amps = primers.in_silico_pcr(tail_f + core_f, tail_r + core_r, seq, "t",
                                 max_mm=0, tp=5, prod_min=100, prod_max=600)
    pair = [a for a in amps if not a["single_primer"]]
    assert len(pair) == 1
    a = pair[0]
    assert a["start"] == left and a["end"] == right                # coordinates stay on the template
    assert a["length"] == right - left                             # templated span unchanged by the tails
    assert a["fwd_tail5"] == len(tail_f) and a["rev_tail5"] == len(tail_r)
    assert a["product_length"] == a["length"] + len(tail_f) + len(tail_r)
    assert a["fwd_anneal_len"] == len(core_f) and a["rev_anneal_len"] == len(core_r)


def test_untailed_primer_reports_no_tail():
    fwd, rev = "GACTGACTGTCAGTCAGGCT", "TTGGCCATTGGCACTGGCAT"
    seq, _l, _r = _amplicon_template(fwd, rev)
    a = [x for x in primers.in_silico_pcr(fwd, rev, seq, "t", max_mm=0, tp=5,
                                          prod_min=100, prod_max=600) if not x["single_primer"]][0]
    assert a["fwd_tail5"] == 0 and a["rev_tail5"] == 0
    assert a["product_length"] == a["length"]                      # the two lengths coincide when nothing dangles


def test_tail_is_not_invented_when_the_whole_primer_anneals():
    # The search tries the full primer first and stops at the first length that binds, so a primer that
    # matches end to end can never be reported as tailed even though a shorter core would also match.
    fwd, rev = "GACTGACTGTCAGTCAGGCT", "TTGGCCATTGGCACTGGCAT"
    seq, _l, _r = _amplicon_template(fwd, rev)
    for a in primers.in_silico_pcr(fwd, rev, seq, "t", max_mm=1, tp=5, prod_min=100, prod_max=600):
        assert a["fwd_tail5"] == 0 and a["rev_tail5"] == 0


def test_min_anneal_bounds_how_far_a_primer_may_be_trimmed():
    core_f, core_r = "GACTGACTGTCAGTCAGGCT", "TTGGCCATTGGCACTGGCAT"      # 20 nt each
    seq, _l, _r = _amplicon_template(core_f, core_r)
    tail = "GGATCCGGAT"                                                   # 10 nt, so the full primer is 30 nt
    kw = dict(max_mm=0, tp=5, prod_min=100, prod_max=600)
    ok = primers.in_silico_pcr(tail + core_f, tail + core_r, seq, "t", min_anneal=20, **kw)
    assert [a for a in ok if not a["single_primer"]], "a 20 nt core must bind when min_anneal is 20"
    strict = primers.in_silico_pcr(tail + core_f, tail + core_r, seq, "t", min_anneal=21, **kw)
    assert strict == [], "a floor above the available core must abolish the site, not shorten it further"


def test_size_window_applies_to_the_product_not_the_templated_span():
    # The window is about the band. A tailed pair whose templated span sits just below the ceiling but whose
    # product exceeds it must be excluded, otherwise prod_max would not mean what a user takes it to mean.
    core_f, core_r = "GACTGACTGTCAGTCAGGCT", "TTGGCCATTGGCACTGGCAT"
    seq, left, right = _amplicon_template(core_f, core_r)
    span = right - left
    tail = "GGATCCGGAT"                                                   # 10 nt on each primer -> +20 bp
    kw = dict(max_mm=0, tp=5, prod_min=100)
    assert primers.in_silico_pcr(tail + core_f, tail + core_r, seq, "t", prod_max=span + 20, **kw)
    assert primers.in_silico_pcr(tail + core_f, tail + core_r, seq, "t", prod_max=span + 19, **kw) == []
    # the same template with untailed primers is unaffected by the distinction
    assert primers.in_silico_pcr(core_f, core_r, seq, "t", prod_max=span, **kw)


def test_min_anneal_default_is_the_design_minimum():
    import inspect
    assert (inspect.signature(primers.in_silico_pcr).parameters["min_anneal"].default
            == primers.MIN_PRIMER_SIZE), "the annealing floor must be the shortest oligo the tool designs"


def test_pcr_seal_records_every_matcher_threshold_at_defaults():
    # A run left entirely at defaults is the case where an under-specified seal is invisible: nothing was
    # typed, so nothing would be recorded unless the resolved values are sealed explicitly.
    import engine
    core_f, core_r = "GACTGACTGTCAGTCAGGCT", "TTGGCCATTGGCACTGGCAT"
    seq, _l, _r = _amplicon_template(core_f, core_r)
    sealed = engine.run_pcr({"sequence": ">t\n" + seq, "fwd": core_f, "rev": core_r})["provenance"]["parameters"]
    for name in ("max_mm", "tp", "min_anneal", "prod_min", "prod_max"):
        assert name in sealed, f"{name} decides which products are reported but is absent from the seal"
    assert sealed["min_anneal"] == primers.MIN_PRIMER_SIZE
