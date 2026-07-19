"""Golden tests: real, committed TE fixtures pinned to their published signatures.
A regression in the detectors will fail these loudly."""
from teagle_core import structural, sequtil
from helpers import fixture_seq


def _feature(seq, kind):
    for e in structural.detect_all(seq):
        if e["type"].startswith(kind):
            return e
    return None


def test_golden_copia_ltr():
    seq = fixture_seq("M11240")                         # Drosophila copia, 5146 bp
    assert len(seq) == 5146
    ltr = _feature(seq, "LTR")
    assert ltr is not None
    assert 266 <= ltr["ltr_len"] <= 286                 # published copia LTR = 276 bp
    assert ltr["identity"] >= 99
    assert max(o["length_aa"] for o in sequtil.find_orfs(seq)) >= 1000   # gag-pol polyprotein


def test_golden_tc1_tir():
    seq = fixture_seq("X01005")                         # C. elegans Tc1, 1610 bp
    tir = _feature(seq, "TIR")
    assert tir is not None
    assert 48 <= tir["tir_len"] <= 58                   # published Tc1 TIR = 54 bp
    assert _feature(seq, "LTR") is None                 # DNA transposon, not an LTR element


def test_golden_l1_line():
    seq = fixture_seq("M80343")                         # human LINE-1 (L1.2), 6050 bp
    assert _feature(seq, "LTR") is None                 # non-LTR
    assert _feature(seq, "TIR") is None
    assert max(o["length_aa"] for o in sequtil.find_orfs(seq)) >= 1000   # ORF2 (~1275 aa)
