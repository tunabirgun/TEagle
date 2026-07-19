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
