"""Scientific contracts for the whole-genome repeat-composition SVG."""
from __future__ import annotations

import os
import sys
import xml.etree.ElementTree as ET

import pytest


NATIVE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "app", "native")
if NATIVE not in sys.path:
    sys.path.insert(0, NATIVE)

import figures


def _overlap_result():
    return {
        "genome_bp": 100,
        # Backend interval unions: categories overlap, while total_bp is their cross-category union.
        "total_bp": 80,
        "te_bp": 70,
        "tandem_bp": 40,
        "other_bp": 20,
        "library_kind": "user-supplied FASTA library",
        "custom_library": {"name": "lab&reference.fa", "sha256": "a5" * 32},
        # These per-name/model unions sum to 160 bp and must have no effect on the figure geometry or remainder.
        "repeat_names": [
            {"repeat_name": "GypsyModel", "class_family": "LTR/Gypsy", "kind": "TE", "bp": 70},
            {"repeat_name": "(TA)n", "class_family": "Simple_repeat", "kind": "tandem", "bp": 50},
            {"repeat_name": "UnknownModel", "class_family": "Unknown", "kind": "other", "bp": 40},
        ],
    }


def _assert_union_figure_contract(svg: str) -> None:
    text = " ".join(ET.fromstring(svg).itertext())
    assert "All repeats (union) 80.00%" in text
    assert "Transposable elements (union) 70.00%" in text
    assert "Tandem repeats (union) 40.00%" in text
    assert "Other/unclassified annotations (union) 20.00%" in text
    assert "Not annotated as repeat 20.00%" in text
    assert "can overlap and must not be summed" in text
    assert "All repeats and not annotated are complementary" in text


def test_composition_uses_backend_unions_not_summed_repeat_name_rows():
    result = _overlap_result()
    svg = figures.svg_te_composition(result)
    _assert_union_figure_contract(svg)

    # A radical change to repeat-name/model rows is an injected version of the old defect. The union figure must not move.
    mutated = dict(result)
    mutated["repeat_names"] = [{"repeat_name": "hATModel", "class_family": "DNA/hAT", "kind": "TE",
                                 "bp": 9_999_999}]
    assert figures.svg_te_composition(mutated) == svg
    assert "Gypsy" not in svg and "9,999,999" not in svg


def test_custom_library_footer_names_fasta_and_checksum_without_dfam():
    result = _overlap_result()
    svg = figures.svg_te_composition(result)
    text = " ".join(ET.fromstring(svg).itertext())
    assert "lab&reference.fa" in text
    assert result["custom_library"]["sha256"] in text
    assert "Dfam" not in text
    assert "TE coverage requires recognized #class/family headers" in text
    assert "simple/low-complexity is separate" in text


@pytest.mark.parametrize("theme", ["light", "dark"])
def test_union_figure_renders_in_both_display_themes(theme):
    pytest.importorskip("PySide6")
    from PySide6.QtCore import QByteArray
    from PySide6.QtSvg import QSvgRenderer

    svg = figures.svg_te_composition(_overlap_result(), theme=theme)
    renderer = QSvgRenderer(QByteArray(svg.encode("utf-8")))
    assert renderer.isValid()
    assert renderer.defaultSize().width() == 900
    assert renderer.defaultSize().height() == 244


def test_installed_library_footer_names_selected_dfam_scope():
    result = _overlap_result()
    result["library_kind"] = "installed Dfam partitions, curated + uncurated"
    result["dfam_version"] = "4.0"
    result["custom_library"] = None
    text = " ".join(ET.fromstring(figures.svg_te_composition(result)).itertext())
    assert "installed Dfam partitions, curated + uncurated" in text
    assert "user-supplied FASTA" not in text


def test_figure_gates_reject_old_stacking_and_hardcoded_dfam_claims():
    result = _overlap_result()
    svg = figures.svg_te_composition(result)
    _assert_union_figure_contract(svg)
    with pytest.raises(AssertionError):
        _assert_union_figure_contract(svg.replace("can overlap and must not be summed", "are stacked and summed"))

    text = " ".join(ET.fromstring(svg).itertext())
    assert "Dfam" not in text
    injected = text + " Homology to the installed Dfam library."
    with pytest.raises(AssertionError):
        assert "Dfam" not in injected
