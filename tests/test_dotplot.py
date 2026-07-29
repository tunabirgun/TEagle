"""Self-similarity dot plot / heat map (roadmap 4.2).

Validated against five real elements rather than synthetic constructs, because the question this panel
answers — "is there a repeat the detectors missed" — can only be judged against sequences whose repeat
content is independently known.

The measurement that shaped the design: the maize Ac element's TIR is 11 bp, and a repeat SHORTER than
the word size cannot produce a single exact k-mer match. At the default k=13 the Ac TIR is invisible
(reverse signal 2); at k=11 it appears (32); at k=8 it is strong (600). Tc1's 54 bp TIR is visible at
every k. Hence the adaptive word size — without it the panel would have silently shown nothing for an
entire superfamily.
"""
import os
import sys

import pytest

_BE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "app", "backend")
if _BE not in sys.path:
    sys.path.insert(0, _BE)
from teagle_core import dotplot, examples, structural                # noqa: E402


def _seq(acc):
    return "".join(l for l in examples.load(acc).splitlines() if not l.startswith(">"))


def _signal(m, which):
    return sum(sum(row) for row in m[which])


def _off_diagonal(m, which="forward", skip=2):
    b = m["bins"]
    return sum(m[which][i][j] for i in range(b) for j in range(b) if abs(i - j) > skip)


# ---------------- the panel shows what the element actually contains ----------------
@pytest.mark.parametrize("acc,label", [("M11240", "copia"), ("M12927", "gypsy")])
def test_ltr_elements_show_off_diagonal_direct_repeats(acc, label):
    """An LTR pair is two copies of the same sequence at opposite ends — an off-diagonal block."""
    m = dotplot.self_matrix(_seq(acc), k=13, bins=120)
    assert _off_diagonal(m, "forward") > 100, f"{label}: no off-diagonal direct-repeat signal"


def test_line_without_terminal_repeats_shows_essentially_only_the_diagonal():
    """L1.2 has no terminal repeat; a strong off-diagonal block here would mean the panel invents one."""
    m = dotplot.self_matrix(_seq("M80343"), k=13, bins=120)
    assert _off_diagonal(m, "forward") < 60


def test_tc1_shows_inverted_repeat_signal():
    m = dotplot.self_matrix(_seq("X01005"), k=13, bins=120)
    assert _signal(m, "reverse") > 20, "Tc1's 54 bp TIR must appear on the reverse layer"


# ---------------- the word-size floor, measured on Ac ----------------
def test_a_repeat_shorter_than_k_is_invisible_and_the_word_size_adapts():
    seq = _seq("X05424")
    tir = structural.find_tir(seq)
    assert tir["tir_len"] < 13, "this test assumes Ac's TIR is shorter than the default word size"
    hidden = dotplot.self_matrix(seq, k=13, bins=120)
    assert _signal(hidden, "reverse") < 10, "expected the short TIR to be invisible at k=13"
    k = dotplot.suggest_k(structural.detect_all(seq))
    assert k < 13, "suggest_k must drop below the default when a short repeat was measured"
    shown = dotplot.self_matrix(seq, k=k, bins=120)
    assert _signal(shown, "reverse") > 100, "the adaptive word size must reveal the short TIR"


def test_suggest_k_keeps_the_default_when_no_repeat_was_measured():
    assert dotplot.suggest_k([]) == dotplot.DEFAULT_K
    assert dotplot.suggest_k(None) == dotplot.DEFAULT_K


# ---------------- guards ----------------
def test_microsatellite_is_masked_rather_than_hanging_the_panel():
    """Every occurrence pairs with every other, so an unmasked (AT)n input is quadratic work."""
    m = dotplot.self_matrix("AT" * 4000, k=13, bins=120)
    assert m["masked_kmers"] > 0
    assert m["masked_positions"] > 1000


def test_short_input_returns_an_empty_matrix_not_an_error():
    m = dotplot.self_matrix("ACGT", k=13, bins=64)
    assert m["forward_max"] == 0 and m["reverse_max"] == 0


def test_matrix_is_square_and_matches_the_declared_bin_count():
    m = dotplot.self_matrix(_seq("X01005"), k=13, bins=64)
    assert len(m["forward"]) == 64 and all(len(r) == 64 for r in m["forward"])
    assert len(m["reverse"]) == 64


# ---------------- chance floor ----------------
def test_chance_floor_is_derived_from_the_sequence_own_composition():
    """An AT-rich element collides far more often than a uniform null predicts, and TEs are often
    AT-rich — so the floor must come from the observed composition, not from 4^-k."""
    at_rich = dotplot.self_matrix("ATATATATGC" * 300, k=8, bins=64)
    balanced = dotplot.self_matrix(("ACGT" * 750), k=8, bins=64)
    assert at_rich["base_collision_p"] > balanced["base_collision_p"]


def test_read_threshold_rises_when_the_word_size_makes_chance_matches_common():
    """At k=13 chance is negligible, so every match counts. At k=8 on a 4.5 kb element roughly 370 of
    19,600 cells hold a chance singleton, so a count of 1 carries no information."""
    seq = _seq("X05424")
    assert dotplot.above_chance(dotplot.self_matrix(seq, k=13, bins=140)) == 1
    assert dotplot.above_chance(dotplot.self_matrix(seq, k=8, bins=140)) >= 2


def test_scope_note_states_the_limits_that_matter():
    m = dotplot.self_matrix(_seq("M11240"), k=13, bins=64)
    note = dotplot.scope_note(m).lower()
    assert "not a local-alignment" in note or "self-blast" in note
    assert "shorter than" in note                       # the k floor
    assert "chance" in note                             # the noise floor
    assert "is not evidence" in note                    # a missing diagonal cannot refute a call


def test_masking_is_reported_in_the_scope_note():
    note = dotplot.scope_note(dotplot.self_matrix("AT" * 4000, k=13, bins=64))
    assert "masked" in note.lower()
