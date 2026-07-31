"""Which Dfam families a single-sequence family annotation searches, and whether it says so.

RepeatMasker reads the CURATED families only unless it is asked for both. Panel 03 never asked, so the
optional uncurated partitions — a 3.9 GB download the installer offers, and 22.6 GiB on disk — were
unreachable there: installing them changed nothing, and a blank result meant "not searched" while
reading as "not present". Live check on this backend, Ty1 (M18706, S. cerevisiae): curated-only returns
0 hits, including uncurated returns 1 hit classed LTR/Copia, which is what Ty1 is.

The family set searched is part of a result's identity, so it is sealed, reported in the result header,
and named in the message shown when nothing matched.
"""
import os
import sys

import pytest
from teagle_core import wsl


def _capture_script(monkeypatch):
    """Run wsl.annotate far enough to see the RepeatMasker command line it builds, without WSL."""
    seen = {}

    def fake_script(script, timeout=90):
        seen.setdefault("scripts", []).append(script)
        return 0, "", ""

    def fake_wsl(script, stdin=None, timeout=600):
        seen.setdefault("scripts", []).append(script)
        return 0, "", ""

    monkeypatch.setattr(wsl, "_wsl_script", fake_script)
    monkeypatch.setattr(wsl, "_wsl", fake_wsl)
    monkeypatch.setattr(wsl, "resolve_species", lambda s: {"ok": True})
    monkeypatch.setattr(wsl, "env_status", lambda: {"ready": True, "repeatmasker": "4.2.4", "dfam": True})
    return seen


@pytest.mark.parametrize("unc,expect", [(False, False), (True, True)])
def test_uncurated_flag_reaches_the_repeatmasker_command(monkeypatch, unc, expect):
    seen = _capture_script(monkeypatch)
    wsl.annotate(">x\nACGT", species="Saccharomyces cerevisiae", include_uncurated=unc)
    joined = "\n".join(seen.get("scripts", []))
    assert ('-species "Saccharomyces cerevisiae"' in joined)
    assert ("-uncurated" in joined) is expect


def test_no_species_never_smuggles_the_flag_in(monkeypatch):
    """The flag rides on the -species argument; without a lineage there is no argument to carry it."""
    seen = _capture_script(monkeypatch)
    wsl.annotate(">x\nACGT", species=None, include_uncurated=True)
    assert "-uncurated" not in "\n".join(seen.get("scripts", []))


def test_engine_threads_the_choice_through_and_seals_it(monkeypatch):
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                    "app", "backend"))
    import engine
    calls = {}

    def fake_annotate(fasta, species=None, timeout=600, include_uncurated=False, **kw):
        calls["include_uncurated"] = include_uncurated
        return {"ok": True, "hits": [], "n_hits": 0, "repeatmasker_version": "4.2.4",
                "dfam_version": "4.0", "dfam_library": {"version": "4.0", "partitions": []},
                "include_uncurated": include_uncurated,
                "library_kind": ("installed Dfam partitions, curated + uncurated" if include_uncurated
                                 else "installed Dfam partitions, curated families only"),
                "species": species}

    monkeypatch.setattr(engine.wsl, "annotate", fake_annotate)
    for unc in (False, True):
        r = engine.run_annotate({"sequence": ">x\nACGTACGTAC", "species": "Homo sapiens",
                                 "include_uncurated": unc})
        assert calls["include_uncurated"] is unc
        params = (r.get("provenance") or {}).get("parameters") or {}
        # the searched family universe is sealed: the same library answers a different question either way
        assert params.get("library") == r["library_kind"]
        assert ("curated + uncurated" in params["library"]) is unc


def _win():
    pytest.importorskip("PySide6")
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    native = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "app", "native")
    if native not in sys.path:
        sys.path.insert(0, native)
    from PySide6.QtWidgets import QApplication
    QApplication.instance() or QApplication([])
    import main as m
    return m.MainWindow()


def _result(hits, unc, partitions):
    return {"ok": True, "hits": hits, "n_hits": len(hits), "repeatmasker_version": "4.2.4",
            "species": "Saccharomyces cerevisiae", "include_uncurated": unc,
            "dfam_library": {"version": "4.0", "partitions": partitions}}


_HIT = [{"class_family": "LTR/Copia", "family": "DR003680865", "q_start": 0, "q_end": 100,
         "strand": "+", "divergence": 1.0, "score": 500, "pct_del": None, "pct_ins": None,
         "cons_start": None, "cons_end": None, "cons_length": None, "cons_coverage_pct": None}]
_PARTS = ["0 [root]: curated", "1 [Eukaryota]: uncurated"]


def _texts(widget):
    from PySide6.QtWidgets import QLabel
    return " ".join(l.text() for l in widget.findChildren(QLabel))


def test_result_header_names_what_was_searched_not_what_is_installed():
    """A machine with the uncurated partitions installed still searches curated-only by default;
    labelling that result "curated + uncurated" would credit a search it did not perform."""
    w = _win()
    w._render_family(_result(_HIT, False, _PARTS))
    assert "curated only" in _texts(w.wslBody.parentWidget())
    w._render_family(_result(_HIT, True, _PARTS))
    assert "curated + uncurated" in _texts(w.wslBody.parentWidget())


def test_blank_result_distinguishes_not_installed_from_not_searched():
    w = _win()
    w._render_family(_result([], False, ["0 [root]: curated"]))       # nothing extra on the machine
    t = _texts(w.wslBody.parentWidget())
    assert "only the CURATED Dfam" in t and "Install the optional uncurated" in t
    w._render_family(_result([], False, _PARTS))                      # installed, but not read
    t = _texts(w.wslBody.parentWidget())
    assert "searched the CURATED families only" in t and "were not read" in t
    w._render_family(_result([], True, _PARTS))                       # both searched: a real no-match
    t = _texts(w.wslBody.parentWidget())
    assert "genuine no-match" in t


def test_library_choice_is_offered_only_once_the_partitions_are_installed():
    w = _win()
    base = {"wsl2": True, "ready": True, "distro": "Ubuntu-24.04", "repeatmasker": "4.2.4"}
    w._on_wsl_status({**base, "dfam_library": {"partitions": ["0 [root]: curated"]}})
    # tracked explicitly: isVisible() is False for a widget in an unshown or collapsed card either way
    assert w._uncurated_available is False        # offering it would promise a search the machine cannot run
    w._on_wsl_status({**base, "dfam_library": {"partitions": _PARTS}})
    assert w._uncurated_available is True
    assert [w.wslLibrary.itemData(i) for i in range(w.wslLibrary.count())] == [False, True]
    assert w.wslLibrary.currentData() is False    # curated-only stays the default, as RepeatMasker's own is


def test_the_choice_is_read_from_state_not_from_widget_visibility():
    """A widget inside a collapsed card reports not-visible; reading the choice from that would send
    curated-only while the user had selected otherwise."""
    w = _win()
    base = {"wsl2": True, "ready": True, "distro": "Ubuntu-24.04", "repeatmasker": "4.2.4"}
    w._on_wsl_status({**base, "dfam_library": {"partitions": _PARTS}})
    w.wslLibrary.setCurrentIndex(1)                       # "Include uncurated families"
    w.card_wsl.collapse() if hasattr(w.card_wsl, "collapse") else w.wslLibraryRow.setVisible(False)
    sent = {}
    w.engine.submit = lambda op, body=None, key=None: sent.update(body or {})
    w.state["seq"] = "ACGT" * 40
    w._annotate()
    assert sent.get("include_uncurated") is True
