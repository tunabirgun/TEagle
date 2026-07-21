"""In-silico PCR single-primer products + gel imaging (co-migration, intensity, single-primer,
no-on-target cue). Engine tests need only the backend; the figure tests run headless via offscreen Qt."""
import os, sys
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
import pytest
from teagle_core.primers import in_silico_pcr, reverse_complement

_NATIVE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "app", "native")
if _NATIVE not in sys.path:
    sys.path.insert(0, _NATIVE)


# ---------- engine: single-primer (self-priming) products ----------
def test_single_primer_product_across_inverted_repeat():
    F, R = "GGATCCAAGCTTGAATTCCG", "TGACTGACTGACTGACTGAC"      # R absent from the template
    tmpl = ("T" * 60) + F + ("A" * 160) + reverse_complement(F) + ("C" * 60)
    amps = in_silico_pcr(F, R, tmpl, max_mm=0, tp=5, prod_min=70, prod_max=1000)
    sp = [a for a in amps if a["single_primer"]]
    assert sp, "expected a single-primer (F+F) product across the inverted repeat"
    assert sp[0]["start"] == 60 and sp[0]["end"] == 260 and sp[0]["length"] == 200
    assert sp[0]["fwd_primer"] == sp[0]["rev_primer"] == "F"


def test_two_primer_product_is_not_flagged_single():
    F = "ACGTACGTAAGGCCTTACGT"
    rcR = "TTGGCCAATTGGCCAATTGG"
    R = reverse_complement(rcR)
    tmpl = ("GAGAGAGAGA" * 6) + F + ("CACACACACA" * 16) + rcR + ("GAGAGAGAGA" * 6)
    amps = in_silico_pcr(F, R, tmpl, max_mm=0, tp=5, prod_min=70, prod_max=1000)
    pair = [a for a in amps if not a["single_primer"]]
    assert pair and pair[0]["fwd_primer"] != pair[0]["rev_primer"]


# ---------- gel imaging ----------
@pytest.fixture(scope="module")
def figures():
    pytest.importorskip("PySide6")
    import figures as f
    return f


def test_comigration_collapses_to_one_band_on_target_wins(figures):
    P = figures.GELPAL["dark"]
    amps = [{"length": 200, "on_target": True, "single_primer": False, "fwd_mm": 0, "rev_mm": 0},
            {"length": 200, "on_target": False, "single_primer": False, "fwd_mm": 1, "rev_mm": 0}]
    bands = figures._lane_bands(amps, P)
    assert len(bands) == 1                                   # equal sizes co-migrate into one band
    assert bands[0]["color"] == P["on"] and bands[0]["on"]  # on-target colour is never painted over
    assert bands[0]["count"] == 2 and "1 on-target" in bands[0]["t"] and "1 off-target" in bands[0]["t"]


def test_single_primer_band_colour_and_disjoint_count(figures):
    P = figures.GELPAL["dark"]
    bands = figures._lane_bands([{"length": 350, "on_target": False, "single_primer": True,
                                  "fwd_mm": 0, "rev_mm": 0}], P)
    assert bands[0]["color"] == P["single"] and bands[0]["single"]
    assert bands[0]["t"].count("single-primer") == 1 and "off-target" not in bands[0]["t"]   # counted once


def test_band_intensity_tracks_mismatches(figures):
    assert figures._band_opacity(0) == 1.0
    assert figures._band_opacity(1) < 1.0 and figures._band_opacity(2) < figures._band_opacity(1)
    assert figures._band_opacity(99) == 0.4                 # floored so a faint band stays visible


def test_gel_legend_and_no_on_target_caption(figures):
    lanes = [{"label": "P1", "amplicons": [{"length": 300, "on_target": False, "single_primer": True,
                                            "fwd_mm": 0, "rev_mm": 0}]}]
    svg = figures.svg_gel({"lanes": lanes}, "dark")
    assert ">single-primer<" in svg                         # legend gains the swatch only when present
    assert ">no on-target<" in svg                          # lane has a band but none intended


def test_gel_wraps_after_ten_lanes(figures):
    import re
    def W_H(n):
        L = [{"label": f"P{i}", "amplicons": [{"length": 200, "on_target": True, "fwd_mm": 0, "rev_mm": 0}]}
             for i in range(n)]
        svg = figures.svg_gel({"lanes": L}, "dark")
        m = re.search(r'width="(\d+)" height="(\d+)"', svg)
        return int(m.group(1)), int(m.group(2)), svg
    w2, h2, _ = W_H(2)
    w11, h11, svg11 = W_H(11)
    w25, h25, svg25 = W_H(25)
    assert w11 == w25 and w11 >= w2                 # width caps at a full 10-lane row once wrapping
    assert h11 > h2 and h25 > h11                   # height grows with each added row
    assert svg11.count(">L</text>") == 2            # 11 lanes -> 2 rows -> 2 ladders
    assert svg25.count(">L</text>") == 3            # 25 lanes -> 3 rows -> 3 ladders


def test_gel_regions_follow_wrapped_rows(figures):
    L = [{"label": f"P{i}", "amplicons": [{"length": 200, "on_target": True, "fwd_mm": 0, "rev_mm": 0,
          "start": 0, "end": 200}]} for i in range(15)]
    regs = figures.gel_regions({"lanes": L})
    assert len(regs) == 15
    ys = sorted({round(r["y0"]) for r in regs})
    assert max(ys) - min(ys) > 100                  # second-row boxes sit well below the first row


def test_gel_regions_group_by_size_prefer_on_target(figures):
    lanes = [{"label": "P1", "amplicons": [
        {"length": 200, "on_target": True, "single_primer": False, "fwd_mm": 0, "rev_mm": 0, "start": 10, "end": 210},
        {"length": 200, "on_target": False, "single_primer": False, "fwd_mm": 1, "rev_mm": 0, "start": 900, "end": 1100}]}]
    regs = figures.gel_regions({"lanes": lanes})
    assert len(regs) == 1                                   # one hit-box per migrated band
    assert regs[0]["amplicon"]["on_target"] is True          # right-click targets the intended product
