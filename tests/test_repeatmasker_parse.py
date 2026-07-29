"""Regression tests for the RepeatMasker .out parser (2026-07-28, roadmap 3.3).

Two defects motivated these. The consensus-side columns were discarded entirely, so a hit's coverage of
its family consensus — the thing that separates a full-length element from a 300 bp scrap — was not
available. And RepeatMasker reverses those columns by strand: '+' rows read `begin end (left)` while 'C'
rows read `(left) end begin`. Reading them positionally produced negative consensus lengths on roughly
half of all real hits; the first draft of this parser did exactly that and returned a coverage of -2292%.

The fragment-ID column is the other half: RepeatMasker gives one ID to the several alignment blocks of a
single interrupted element, so a fragmented L1 was rendering as N apparently independent copies.
"""
from teagle_core import wsl

# Real .out column layout. L1PA5 is a plain '+' hit; AluSx is a 'C' hit (the reversed columns);
# the two L1HS rows share fragment ID 7 — one element interrupted, not two copies.
SAMPLE = """
   SW   perc perc perc  query      position in query     matching  repeat      position in repeat
score   div. del. ins.  sequence   begin end   (left)    repeat  class/family  begin end  (left)  ID

 1200   12.5  1.0  0.5  chr1        100   600  (9400) +  L1PA5   LINE/L1          1   501 (5500)   1
 1500   18.2  2.0  1.0  chr1        800  1300  (8700) C  AluSx   SINE/Alu      (12)   300     1    2
  900   10.0  0.5  0.0  chr1       2000  2400  (7600) +  L1HS    LINE/L1          1   401 (5600)   7
  850   10.4  0.6  0.0  chr1       3000  3500  (6500) +  L1HS    LINE/L1        402   902 (5099)   7
"""


def _by_family(hits, name):
    return next(h for h in hits if h["family"] == name)


def test_indel_columns_are_kept():
    h = _by_family(wsl.parse_out(SAMPLE), "L1PA5")
    assert h["pct_del"] == 1.0 and h["pct_ins"] == 0.5


def test_plus_strand_consensus_coordinates():
    h = _by_family(wsl.parse_out(SAMPLE), "L1PA5")
    assert (h["cons_start"], h["cons_end"], h["cons_left"]) == (1, 501, 5500)
    assert h["cons_length"] == 6001


def test_minus_strand_columns_are_not_read_positionally():
    """The 'C'-row trap: (left) end begin. Read in '+' order this yields a negative consensus length."""
    h = _by_family(wsl.parse_out(SAMPLE), "AluSx")
    assert (h["cons_start"], h["cons_end"], h["cons_left"]) == (1, 300, 12)
    assert h["cons_length"] == 312
    assert 90 <= h["cons_coverage_pct"] <= 100          # a near-full-length Alu


def test_every_coverage_is_physically_possible():
    for h in wsl.parse_out(SAMPLE):
        if h["cons_length"] is not None:
            assert h["cons_length"] > 0, h["family"]
        if h["cons_coverage_pct"] is not None:
            assert 0 <= h["cons_coverage_pct"] <= 100, (h["family"], h["cons_coverage_pct"])


def test_fragments_sharing_an_id_merge_into_one_hit():
    hits = wsl.parse_out(SAMPLE)
    l1hs = [h for h in hits if h["family"] == "L1HS"]
    assert len(l1hs) == 1, "an interrupted element must not read as several independent copies"
    h = l1hs[0]
    assert h["n_fragments"] == 2
    assert (h["q_start"], h["q_end"]) == (1999, 3500)   # outermost query span, 0-based start
    assert (h["cons_start"], h["cons_end"]) == (1, 902)
    assert "not 2 separate copies" in h["fragment_note"]


def test_distinct_families_never_merge():
    hits = wsl.parse_out(SAMPLE)
    assert {h["family"] for h in hits} == {"L1PA5", "AluSx", "L1HS"}


def test_malformed_row_drops_only_itself():
    bad = SAMPLE + "\n  xx   bad  row  that  cannot  parse\n"
    assert len(wsl.parse_out(bad)) == len(wsl.parse_out(SAMPLE))


def test_divergence_is_labelled_raw_not_kimura():
    """The parser must not imply a Kimura correction it does not apply."""
    doc = (wsl.parse_out.__doc__ or "").lower()
    assert "raw" in doc and "kimura" in doc
