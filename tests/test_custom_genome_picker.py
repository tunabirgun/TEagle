"""Custom pairs reach the existing genome scan without changing its request."""
import os
import sys
from pathlib import Path
from types import SimpleNamespace

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "app/native"), str(ROOT / "app/backend")]

from PySide6.QtWidgets import QApplication, QComboBox, QPushButton
from panel_genome import GenomeScanMixin

APP = QApplication.instance() or QApplication([])


class Window(GenomeScanMixin):
    def __init__(self):
        self.state = {"candidates": [], "pcrPairs": []}
        self.scanPicker = QComboBox()
        self.scanBtn = QPushButton()
        self.pcrStageAll = QPushButton()
        self.genomeOrg = QComboBox()
        self.genomeOrg.addItem("Choose", None)
        self.genomeOrg.addItem("Synthetic genome", "synthetic")
        self.genomeOrg.setCurrentIndex(1)
        self._genome_inflight = self._genome_prep_inflight = False
        self._specimen_gen = 1
        self.card_genome = SimpleNamespace(expand=lambda: None)
        self.requests = []
        self.messages = []
        self.engine = SimpleNamespace(submit=lambda *a, **kw: self.requests.append((a, kw)))

    _pair_label = staticmethod(lambda c: c["id"])
    _clear_banner = _refresh_genome_manager = lambda self: None
    _render_genome_busy = lambda self, text: None
    _design_locus_for = lambda self, cand, org: None
    _banner = lambda self, text, level: self.messages.append((text, level))


def pair(name="custom1", **extra):
    return {"id": name, "left_seq": "ACGTACGTACGTACGTAC", "right_seq": "TGCATGCATGCATGCATG", **extra}


def test_custom_pair_is_selectable_without_a_design():
    w = Window()
    c = pair(custom=True)
    w.state["pcrPairs"] = [c]
    w._refresh_scan_picker()
    assert w.scanPicker.currentData() == c
    assert w.scanBtn.isEnabled()
    assert not w.pcrStageAll.isEnabled()
    w._scan_from_picker()
    args, kwargs = w.requests[-1]
    assert args == ("genome_pcr", {"fwd": c["left_seq"], "rev": c["right_seq"], "organism": "synthetic", "design_locus": None,
                                  "params": {"min_perfect": 15, "min_good": 15, "prod_max": 4000, "prod_min": 0}})
    assert kwargs == {"key": "genome_pcr"}


def test_only_current_designs_and_queued_custom_pairs_are_offered():
    w = Window()
    designed, custom, stale = pair("P1"), pair("lab_pair", custom=True), pair("old_design")
    w.state.update(candidates=[designed], pcrPairs=[designed, custom, stale])
    w._refresh_scan_picker()
    assert [w.scanPicker.itemData(i) for i in range(w.scanPicker.count())] == [designed, custom]
    assert w.pcrStageAll.isEnabled()
    w.state.update(candidates=[], pcrPairs=[])
    w._refresh_scan_picker()
    assert w.scanPicker.currentData() is None
    assert not w.scanBtn.isEnabled()


def test_custom_and_designed_pairs_use_identical_search_controls():
    requests = []
    for c in (pair("P1", product_size=224), pair("lab_pair", custom=True)):
        w = Window()
        w._scan_genome(c, "synthetic")
        requests.append(w.requests[-1])
    assert requests[0] == requests[1]


def test_ambiguous_custom_pair_is_not_silently_sent_to_genome_scan():
    w = Window()
    c = pair(custom=True)
    c["left_seq"] = "NCGTACGTACGTACGTAC"
    w._scan_genome(c, "synthetic")
    assert not w.requests
    assert "A/C/G/T" in w.messages[0][0]
    assert not w._genome_inflight
