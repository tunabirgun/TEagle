"""LTR cis-element detectors added in the sub-structure batch: the advisory poly(A)-signal motif and the
non-canonical terminal-motif panel. Both are BADGES — neither may move a boundary or assert a cleavage site."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "app", "backend"))
from teagle_core import structural


def _ltr(l5s, l5e, l3s, l3e):
    return {"five_prime": [l5s, l5e], "three_prime": [l3s, l3e], "element_span": [l5s, l3e]}


def _seq_with_pas(hexamer="AATAAA", dse="TTTTGTTTGTTTTGTTGTTTTGTTTGTTTTGTTGTTTTGTTTGTTG", tail=""):
    """Build a record whose 3' LTR carries `hexamer` followed by the real DSE window.

    The detector scores 20-60 nt PAST the hexamer (where the downstream element actually sits, given the
    hexamer is 10-30 nt upstream of the cleavage site and the DSE begins ~15 nt after it), so the spacer
    below is sized to put `dse` inside that window rather than in the hexamer-to-cleavage gap."""
    pre = "C" * 100                                  # 5' LTR copy region placeholder
    ltr3 = "A" * 20 + hexamer + "C" * 20 + dse
    return pre + "G" * 100 + ltr3 + tail, len(pre) + 100


def test_polya_signal_requires_a_downstream_element():
    """A hexamer alone is a chance word (~1 per 4 kb); without the GU/U-rich DSE nothing may be reported."""
    seq, s3 = _seq_with_pas(dse="C" * 50)            # window present but not GU/U-rich -> gate fails
    assert structural.find_polya_signal(seq, _ltr(0, 100, s3, s3 + 56)) is None


def test_polya_signal_reports_a_gated_hit_with_its_hedge():
    seq, s3 = _seq_with_pas()
    pas = structural.find_polya_signal(seq, _ltr(0, 100, s3, s3 + 56))
    assert pas is not None and pas["motif"] == "AATAAA" and pas["variant"] == "canonical"
    assert pas["dse_state"] == "found"
    assert pas["type"].startswith("polyA-signal")    # never "poly-A", which is the LINE TAIL feature
    low = pas["note"].lower()
    assert "advisory" in low and "does not locate" in low
    assert "cleavage site" in low                    # the limit that matters must travel with the call


def test_polya_signal_not_assessable_when_the_record_ends():
    """Three outcomes, not two: no downstream window means the question could not be asked."""
    seq, s3 = _seq_with_pas(dse="")                  # record stops right after the hexamer
    pas = structural.find_polya_signal(seq, _ltr(0, 100, s3, len(seq)))
    assert pas is not None and pas["dse_state"].startswith("not assessable")
    assert pas["confident"] is False


def test_noncanonical_motif_panel_is_exact_match_only():
    """A mismatch-tolerant search over eight motifs would match nearly any terminus and become noise."""
    # TG..GA is the non-canonical motif TGGA from LTR_retriever's default panel
    seq = "TG" + "C" * 96 + "GA" + "T" * 50 + "TG" + "C" * 96 + "GA"
    tm = structural._termini_motif(seq, 0, 100, 150, 250)
    assert tm["canonical"] is False
    assert tm["noncanonical_motif"] == "TGGA" and tm["motif_tier"] == "non-canonical"
    # a terminus matching nothing stays "none" rather than being snapped onto the nearest motif
    seq2 = "AG" + "C" * 96 + "TT" + "T" * 50 + "AG" + "C" * 96 + "TT"
    tm2 = structural._termini_motif(seq2, 0, 100, 150, 250)
    assert tm2["noncanonical_motif"] is None and tm2["motif_tier"] == "none"
    assert "not evidence against" in tm2["note"]


def test_noncanonical_panel_matches_ltr_retriever_default():
    """Verified against the LTR_retriever 2.9.0 manual (Table 1): -motif [TCCA TGCT TACA TACT TGGA TATA TGTA TGCA]."""
    assert structural._CANONICAL_MOTIF == "TGCA"
    assert set(structural._NONCANONICAL_MOTIFS) == {"TCCA", "TGCT", "TACA", "TACT", "TGGA", "TATA", "TGTA"}
