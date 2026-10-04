"""4.1 UI onboarding and comprehension lane: header and panel naming, record selection, one state vocabulary,
caution strip, row actions, quick start and glossary. Drives the real MainWindow offscreen; nothing leaves the process
and the registry is never written (QSettings is replaced by an in-memory store)."""
import copy, functools, os, sys

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
import pytest

pytest.importorskip("PySide6")
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_APP = os.environ.get("M2_APP_ROOT") or os.path.join(_ROOT, "app")      # a pre-change copy can be pointed at for failing-before runs
sys.path.insert(0, os.path.join(_APP, "native"))
sys.path.insert(0, os.path.join(_APP, "backend"))
sys.path.insert(0, os.path.join(_ROOT, "tests"))

from PySide6.QtCore import Qt, QEvent, QPoint, QTimer, QEventLoop
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication, QLabel, QPushButton, QFrame

_app = QApplication.instance() or QApplication([])


class FakeSettings:
    store = {}

    def __init__(self, *a):
        pass

    def value(self, key, default=None):
        return self.store.get(key, default)

    def setValue(self, key, value):
        self.store[key] = value


def _spin(pred, timeout=5.0):
    import time
    loop, t, t0 = QEventLoop(), QTimer(), time.time()
    t.setInterval(20)
    t.timeout.connect(lambda: (pred() or time.time() - t0 > timeout) and loop.quit())
    t.start(); loop.exec(); t.stop()
    return pred()


def _make_win(monkeypatch):
    import fonts, main
    fonts.load_fonts()
    monkeypatch.setattr(main, "QSettings", FakeSettings)
    w = main.MainWindow()
    _app.processEvents()
    return w


def _dispose(w):
    import gc
    w.engine.wait(5000)
    w.close(); w.deleteLater()
    QApplication.sendPostedEvents(None, QEvent.DeferredDelete)
    gc.collect()


@pytest.fixture
def win(monkeypatch):
    FakeSettings.store = {"quickstart_hidden": True}
    w = _make_win(monkeypatch)
    yield w
    _dispose(w)


@pytest.fixture
def fresh_win(monkeypatch):
    """A window as a first-time user sees it (no stored preferences)."""
    FakeSettings.store = {}
    w = _make_win(monkeypatch)
    yield w
    _dispose(w)


@functools.lru_cache(maxsize=None)
def _analysis(name):
    import engine, helpers
    from teagle_core import examples, sequtil
    if name in ("M11240", "M12927", "M80343", "X01005", "X05424"):
        return engine.analyze(examples.load(name))
    if name == "sample":
        from sample import make_sample
        return engine.analyze(make_sample())
    if name == "iupac":
        return engine.analyze(">iupac\n" + "ACGTRYKMSWNN" * 8)
    if name == "array230":
        s = sequtil.parse_fasta(examples.load("M11240"))[0][1]
        return engine.analyze(">arr\n" + (s * 50)[:230000])
    if name == "copia2":
        s = sequtil.parse_fasta(examples.load("M11240"))[0][1]
        return engine.analyze(">dimer\n" + s + s)
    if name == "multi":
        recs = [sequtil.parse_fasta(examples.load(n))[0][1] for n in ("M11240", "X01005", "M80343")]
        return engine.analyze("\n".join(f">rec{i}\n{s}" for i, s in enumerate(recs, 1)))
    raise KeyError(name)


def _accept(win, name, result=None, copy_result=True):
    result = copy.deepcopy(result or _analysis(name)) if copy_result else result
    sequence = "\n".join(f">{r['id']}\n{r['seq']}" for r in result.get("records", []) if r.get("seq"))
    win.seq.setPlainText(sequence)
    win._analyze_gen += 1
    generation = win._analyze_gen
    win._active_analyze_gen = generation
    win._analysis_requests[generation] = {"sequence": sequence, "source": {}}
    win._analyze_inflight = True
    win._on_analyze(result, generation)
    return win.state["records"][win.state.get("record_index", 0)]


def _labels(widget):
    return [l.text() for l in widget.findChildren(QLabel)]


def _texts(widget):
    return " ".join(_labels(widget))


def _by_name(widget, name):
    return [w for w in widget.findChildren(QLabel) if w.objectName() == name]


# ================================================================ header and panel naming (B23, B18 names, B30 casing)
_WSL_STATES = [      # payload, summary, family (RepeatMasker + Dfam), splice (minimap2), WSL2
    ({"error": "boom"}, "status unknown", "status unknown", "status unknown", "status unknown"),
    ({"wsl2": False}, "not installed", "not installed", "not installed", "not installed"),
    ({"wsl2": True, "ready": False}, "not set up", "not set up", "not set up", "WSL2 found"),
    ({"wsl2": True, "ready": True, "minimap2": None}, "partly installed", "ready", "not set up", "WSL2 found"),
    ({"wsl2": True, "ready": True, "minimap2": "2.28"}, "ready", "ready", "ready", "WSL2 found"),
]


@pytest.mark.parametrize("health_first", [True, False])
@pytest.mark.parametrize("payload,expected,_rm,_mm,_wsl", _WSL_STATES)
def test_status_chip_combines_core_and_backend_states(win, payload, expected, _rm, _mm, _wsl, health_first):
    assert win.statusTxt.text() == "connecting…" and "checking" in win.backendTxt.text()
    steps = [lambda: win._on_done("health", {"primer3": "2.3.0"}), lambda: win._on_wsl_status(dict(payload))]
    for s in (steps if health_first else steps[::-1]):
        s()
    assert win.statusTxt.text() == "engine ready · primer3 2.3.0"
    assert win.backendTxt.text() == f"Linux backend: {expected}"
    assert "always available" in win._chip.toolTip() and "Optional Linux (WSL) backend" in win._chip.toolTip()


def test_chip_says_when_the_core_engine_failed(win):
    win._on_failed("health", "x", "")
    assert win.statusTxt.text() == "engine unavailable" and win.led.property("live") is False
    assert "not available" in win._chip.toolTip()


@pytest.mark.parametrize("payload,_summary,rm,mm,wsl", _WSL_STATES)
def test_backend_tags_on_optional_panels(win, payload, _summary, rm, mm, wsl):
    win._on_wsl_status(dict(payload))
    assert f"optional · Linux backend: {rm}" in win.card_wsl.hdr.text()
    assert f"optional · Linux backend: {mm}" in win.card_splice.hdr.text()
    assert f"optional · Linux backend: {wsl}" in win.card_genome.hdr.text()
    assert "ready" not in win.card_genome.hdr.text().split("optional")[1].split("·")[1]       # the scan tool is never certified here
    for card in (win.card_struct, win.card_primer, win.card_pcr, win.card_prov):
        assert "Linux backend" not in card.hdr.text()


def test_set_meta_keeps_the_collapse_state_and_restores_the_default(win):
    card = win.card_primer
    card.expand()
    default = card.hdr.text()
    card.set_meta("x")
    assert card.hdr.text().startswith("▾") and card.hdr.text().endswith("x") and not card.body.isHidden()
    card.set_meta(None)
    assert card.hdr.text() == default
    card.set_collapsed(True)
    card.set_meta("y")
    assert card.hdr.text().startswith("▸") and card.body.isHidden()


def test_theme_and_text_size_buttons_show_state(win):
    import main, theme as theme_mod
    orig = theme_mod.UI_SCALE
    try:
        assert win.theme == "light" and win.themeBtn.text() == "DARK MODE"
        win._toggle_theme()
        assert win.themeBtn.text() == "LIGHT MODE" and "light" in win.themeBtn.accessibleName()
        win._set_ui_scale(1.25)
        assert win.scaleBtn.text() == "TEXT SIZE 125%" and "125%" in win.scaleBtn.accessibleName()
    finally:
        theme_mod.UI_SCALE = orig


def test_rail_toggle_names_what_it_hides_and_keeps_its_casing(win):
    assert win.railToggle.text() == "◂ HIDE INPUT" and win.railToggle.accessibleName() == "Hide input panel"
    assert "panel 01" in win.railToggle.toolTip()
    win._toggle_rail()
    assert win.railToggle.text() == "▸ SHOW INPUT" and win.railToggle.accessibleName() == "Show input panel"


@pytest.mark.parametrize("scale", [1.0, 1.5])
def test_header_never_raises_window_floor(monkeypatch, scale):
    import theme as theme_mod
    orig = theme_mod.UI_SCALE
    theme_mod.UI_SCALE = scale
    try:
        FakeSettings.store = {"quickstart_hidden": True}
        w = _make_win(monkeypatch)
        try:
            w._on_done("health", {"primer3": "2.3.0"})
            w._on_wsl_status({"wsl2": True, "ready": True, "minimap2": None})
            w._toggle_rail(); w._toggle_rail()
            header = w.centralWidget().layout().itemAt(0).widget()
            assert header.minimumSizeHint().width() <= round(820 * scale)
        finally:
            _dispose(w)
    finally:
        theme_mod.UI_SCALE = orig


def test_initial_focus_is_the_accession_field(win):
    win.show()
    _app.processEvents()
    QTest.qWaitForWindowActive(win, 2000)
    assert QApplication.focusWidget() is win.acc


def test_panel_numbers_are_one_table_and_the_splice_panel_is_03b(win):
    import main, re
    assert main.PANEL["splice"] == "03b" and len(set(main.PANEL.values())) == len(main.PANEL)
    assert win.card_splice.hdr.text().startswith("▸ 03b  Splice structure (transcript alignment)")
    _accept(win, "M11240")
    win._render_primers({"candidates": [{"id": "P1", "left_seq": "ACGT" * 5, "right_seq": "TGCA" * 5, "product_size": 200,
                                         "left_tm": 60, "right_tm": 60, "left_gc": 50, "right_gc": 50, "penalty": 0.1,
                                         "qc": {}}], "oligoqc_engines": {}})
    win._render_splice({"counts": {"exons": 2, "introns": 1}, "canonical_introns": 1, "strand": "+",
                        "exons": [{"start": 0, "end": 3}, {"start": 22, "end": 25}],
                        "introns": [{"start": 3, "end": 22, "donor": "GT", "acceptor": "AG", "canonical": True}]})
    seen = []
    for w in win.findChildren(QLabel) + win.findChildren(QPushButton):
        seen += [w.text(), w.toolTip()]
    known = set(main.PANEL.values())
    for text in seen:
        for m in re.finditer(r"[Pp]anel (\d\d[a-z]?)", text):
            assert m.group(1) in known, (m.group(0), text[:80])


def test_splice_panel_text_has_no_de_novo(win):
    texts = [_texts(win.card_splice), win.card_splice.hdr.text()]
    win._on_wsl_status({"wsl2": False})
    texts.append(_texts(win.card_splice))
    win.state["splice_seq"] = "AAAGT" + "C" * 15 + "AGTTT"
    win.state["last_rec"] = {"id": "splice_record"}
    win._render_splice({"counts": {"exons": 2, "introns": 1}, "canonical_introns": 1, "strand": "+",
                        "exons": [{"start": 0, "end": 3}, {"start": 22, "end": 25}],
                        "introns": [{"start": 3, "end": 22, "donor": "GT", "acceptor": "AG", "canonical": True}]})
    root = win.card_splice
    from widgets import DataTable
    for w in root.findChildren(QPushButton):
        texts.append(w.toolTip())
    for t in root.findChildren(DataTable):
        texts += list(getattr(t, "export_notes", None) or [])
    texts.append(_texts(win.card_splice))
    joined = " ".join(texts).lower()
    assert "de novo" not in joined and "de-novo" not in joined
    assert "detected de novo" in _texts(win.rail)                 # the rail caption describes terminal-repeat heuristics and stays


# ================================================================ visible record selection (B6, U9)
def _cand(i=1):
    return {"id": f"P{i}", "left_seq": "ACGTACGTAGCTAGCTAGC" + "ACGT"[i % 4], "right_seq": "TGCATGCATCGATCGATCG" + "ACGT"[i % 4], "product_size": 200,
            "left_tm": 60.0, "right_tm": 60.0, "left_gc": 50.0, "right_gc": 50.0, "penalty": 0.1, "qc": {}}


def _title(win):
    return [l.text() for l in _by_name(win.card_struct, "sech") if l.text().startswith(("Record ", "record_"))]


def test_record_row_selection_binds_without_rebuilding_table(win):
    import shiboken6, main
    _accept(win, "multi")
    recs = win.state["records"]
    t = win._recTable
    col = win._RECORD_HEADERS.index("Length")
    t.sortByColumn(col, Qt.DescendingOrder)
    visual = next(r for r in range(t.rowCount()) if t._orig(r) == 2)
    t.setCurrentCell(visual, 0)
    assert _spin(lambda: win.state["record_index"] == 2)
    assert win._recTable is t and shiboken6.isValid(t)
    assert t.horizontalHeader().sortIndicatorSection() == col
    assert win.state["seq"] == recs[2]["seq"]
    assert _title(win)[0].startswith("Record 3 of 3 · ")
    rid = main._record_public_id(recs[2])
    for text in (win.card_primer.hdr.text(), win.card_pcr.hdr.text(), win.designHint.text()):
        assert rid in text
    assert not win.mScope.isHidden() and "3 of 3" in win.mScope.text() and "this record only" in win.mScope.text()
    assert t.current_orig() == 2


def test_record_row_mouse_and_keyboard_paths(win):
    win.show(); _app.processEvents()
    _accept(win, "multi")
    t = win._recTable
    t.setFocus()
    QTest.mouseClick(t.viewport(), Qt.LeftButton, pos=t.visualRect(t.model().index(1, 1)).center())
    assert _spin(lambda: win.state["record_index"] == t._orig(1))
    before = win.state["record_index"]
    QTest.keyClick(t, Qt.Key_Down)
    assert _spin(lambda: win.state["record_index"] != before)
    assert win.state["record_index"] == t.current_orig()


def test_a_burst_of_row_changes_renders_once_and_the_bound_row_wins(win, monkeypatch):
    _accept(win, "multi")
    t = win._recTable
    shown = []
    real = win._show_record
    monkeypatch.setattr(win, "_show_record", lambda i, keep_table=False: (shown.append(i), real(i, keep_table=keep_table))[1])
    for r in (1, 2, 0):
        t.setCurrentCell(r, 0)
    _spin(lambda: False, 0.4)
    assert shown == [] and win.state["record_index"] == 0           # moving back to the bound row cancels the pending switch
    t.setCurrentCell(2, 0)
    assert _spin(lambda: shown == [2])


def test_invalidation_clears_the_record_binding_and_a_pending_selection(win):
    _accept(win, "multi")
    win._recTable.setCurrentCell(2, 0)
    win._invalidate_analysis("x")
    _spin(lambda: False, 0.4)
    assert win.state.get("record_index") is None and win._recTable is None and win._recDetail is None
    assert "x" in _texts(win.card_struct)
    assert win.card_primer.hdr.text().endswith("designs PCR primers (Primer3)")
    assert win.designHint.text() == "run analysis first" and win.mScope.isHidden()


def test_a_single_record_shows_its_id_without_a_scope_label(win):
    import main
    rec = _accept(win, "M11240")
    rid = main._record_public_id(rec)
    assert win.mScope.isHidden()
    assert _title(win)[0].startswith(rid + " · ")
    assert rid in win.designHint.text() and "record 1 of" not in win.designHint.text()


def test_primers_label_their_design_record_after_a_switch(win):
    _accept(win, "multi")
    recs = win.state["records"]
    win._design_tmpl = recs[0]["seq"]
    win._render_primers({"candidates": [_cand()], "oligoqc_engines": {}})
    assert any("Designed on record 1 of 3" in t for t in _labels(win.primBody.parentWidget()))
    win._show_record(1, keep_table=True)
    warn = [l for l in win.primBody.parentWidget().findChildren(QLabel) if l.property("level") == "warn"]
    assert warn and "designed on record 1 of 3" in warn[0].text() and "record 2 of 3" in warn[0].text()


def test_the_context_menu_route_moves_the_row_highlight_too(win):
    _accept(win, "multi")
    t = win._recTable
    dict(t._row_menu(2))                                             # builder still offers the action
    label, fn = t._row_menu(2)[0]
    assert label == "Show this record below"
    fn()
    assert win.state["record_index"] == 2 and t.current_orig() == 2 and win._recTable is t


# ================================================================ one state vocabulary, call confidence, rank line (B7, B10, B26, U10, U12)
def _tables(win):
    from widgets import DataTable
    return win.card_struct.findChildren(DataTable)


def _summary_row(win):
    t = _tables(win)[0]
    r = next(i for i in range(t.rowCount()) if t._orig(i) == win.state["record_index"])
    return {t.horizontalHeaderItem(j).text(): t.item(r, j).text() for j in range(t.columnCount())}


def _chips(win):
    return [l.text() for l in _by_name(win.card_struct, "cf")]


def _arch_line(win):
    text = _texts(win.card_struct)
    if "Domain architecture:" not in text:
        return None
    return text.split("Domain architecture:")[1].split("Expected modules")[0].split("Coding structure")[0].strip()


@pytest.mark.parametrize("name", ["M11240", "M80343", "X01005", "sample", "iupac"])
def test_summary_row_and_card_share_one_vocabulary(win, name):
    _accept(win, name)
    row, chips = _summary_row(win), _chips(win)
    conf = [c for c in chips if c.startswith("CALL CONFIDENCE:")]
    held = [c for c in chips if "WITHHELD" in c or c == "NOT ASSESSED"]
    if conf:
        assert row["Call confidence"].upper() == conf[0].split(": ")[1].split(" · ")[0]
    else:
        assert row["Call confidence"].startswith("n/a")
    assert (row["Call status"] == "withheld") == bool(held)
    if _arch_line(win) is not None:
        assert row["Domain architecture"] in _arch_line(win)


def test_confidence_chip_has_a_noun_a_tooltip_and_an_accessible_name(win):
    _accept(win, "M80343")
    chip = next(l for l in _by_name(win.card_struct, "cf"))
    assert chip.text() == "CALL CONFIDENCE: MODERATE" and chip.accessibleName() == "Call confidence: Moderate"
    assert "not a probability" in chip.toolTip() and "that class's ceiling" in chip.toolTip()
    assert "separate" in chip.toolTip() and "Hit confidence" in chip.toolTip()


def test_no_signature_input_reads_not_assessable_and_has_no_chip(win):
    _accept(win, "iupac")
    assert "No TE signature detected" in _texts(win.card_struct) and not _chips(win)
    row = _summary_row(win)
    assert row["Domain evidence"] == "not assessable (no ORF ≥ 40 aa)" and row["Call confidence"] == "n/a — no TE signature"
    assert "no open reading frame" in _texts(win.card_struct) and "not a finding" in _texts(win.card_struct)


def test_a_failed_orf_stage_does_not_read_as_not_assessable(win, monkeypatch):
    import engine
    from teagle_core import sequtil
    real = sequtil.find_orfs

    def boom(seq, min_aa=40):                                  # same signature: the seal reads the detector defaults
        raise RuntimeError("forced")
    monkeypatch.setattr(sequtil, "find_orfs", boom)
    res = engine.analyze(">x\n" + "ACGTRYKMSWNN" * 8)
    monkeypatch.setattr(sequtil, "find_orfs", real)
    _accept(win, "x", result=res)
    row = _summary_row(win)
    assert row["Domain evidence"] == "not assessed (ORF detection failed)"      # not "no ORF to search": the ORF stage never completed


def test_every_summary_header_has_a_gloss(win):
    _accept(win, "M11240")
    t = _tables(win)[0]
    assert [t.horizontalHeaderItem(i).toolTip() for i in range(t.columnCount())].count("") == 0
    assert win._RECORD_HEADERS[:7] == ["#", "Record", "Length", "Classification", "Call confidence", "Completeness", "Call status"]


def test_domain_table_says_hit_confidence(win):
    _accept(win, "M11240")
    t = next(t for t in _tables(win) if str(t.export_base).endswith("_domains"))
    heads = [t.horizontalHeaderItem(i).text() for i in range(t.columnCount())]
    assert "Hit confidence" in heads and "Conf" not in heads
    assert t.horizontalHeaderItem(heads.index("Hit confidence")).toolTip()
    assert "Hit confidence</b>" in _texts(win.card_struct)


def test_the_summary_table_export_uses_plain_words_and_names_the_machine_tokens(win):
    _accept(win, "M80343")
    t = _tables(win)[0]
    out = dict(zip(t.export_headers(), t.rows_data()[0]))
    assert out["Domain evidence"].startswith("hits found (partial: ") and out["Call status"] == "called"
    assert any("tested_with_hits" in n for n in t.export_notes)
    raw = dict(zip(win._RECORD_HEADERS, win._record_summary_rows(win.state["records"])[0]))
    assert raw["Domain evidence"].startswith("tested_with_hits")          # the GFF3/JSON vocabulary is unchanged


def test_call_confidence_ceiling_matches_the_tooltip():
    import engine
    from teagle_core import examples
    got = {}
    for acc in ("M11240", "M12927", "M80343", "X01005", "X05424"):
        cl = engine.analyze(examples.load(acc))["records"][0]["classification"]
        got[acc] = (cl["te_class"], cl["confidence"])
    assert got["M11240"][1] == "High" and got["M12927"][1] == "High"
    assert all(conf != "High" for acc, (te, conf) in got.items() if not te.startswith("LTR/"))
    assert got["M80343"][1] == "Moderate" and got["X01005"][1] == "Moderate"


def test_render_does_not_mutate_the_sealed_result(win):
    import engine
    for name in ("M11240", "array230", "M80343"):
        res = copy.deepcopy(_analysis(name))
        before = copy.deepcopy(res)
        _accept(win, name, result=res, copy_result=False)
        win._render_struct_card(win.state["records"][0], res)
        assert res == before
        engine.analysis_result_document(res)


@pytest.mark.parametrize("name,code,rm", [("M11240", "RLC", "LTR/Copia"), ("M12927", "RLG", "LTR/Gypsy"),
                                          ("M80343", "RIX", "LINE"), ("X01005", "DTT", "DNA/TcMar"),
                                          ("X05424", "DTA", "DNA/hAT")])
def test_rank_line_on_card_and_table_agree(win, name, code, rm):
    _accept(win, name)
    row = _summary_row(win)
    line = next(l for l in _labels(win.card_struct) if l.startswith("Wicker 2007:"))
    assert f"({code})" in line and f"RepeatMasker class: {rm}" in line
    assert row["Wicker code"] == code and row["RM class"] == rm


_RM_CLASSES = {"LTR/Copia", "LTR/Gypsy", "LTR/DIRS", "LTR", "LINE", "SINE", "DNA", "DNA/hAT", "DNA/TcMar", "DNA/TcMar-Tc1",
               "DNA/CMC-EnSpm", "DNA/MULE-MuDR", "DNA/PiggyBac", "RC/Helitron"}


def test_the_mapped_repeatmasker_classes_are_strings_that_library_uses():
    import main
    shown = {v[3] for v in main._WICKER.values() if v[3]}
    assert shown <= _RM_CLASSES, shown - _RM_CLASSES
    assert main._rank_view(_matrix()["mule"])["rm"] == "DNA/MULE-MuDR"
    for te, v in main._WICKER.items():                                 # a rank not called stops at the order, never a made-up subclass
        if v[1] is None and v[0] == "TIR":
            assert v[3] == "DNA", te


def test_the_no_orf_input_shows_no_rank_line(win):
    _accept(win, "iupac")
    assert not [l for l in _labels(win.card_struct) if l.startswith("Wicker 2007:")]
    assert _summary_row(win)["Wicker code"] == "—"


def _score_module():
    import importlib.util
    spec = importlib.util.spec_from_file_location("bench_score_m2", os.path.join(_ROOT, "benchmarks", "score.py"))
    mod = importlib.util.module_from_spec(spec)
    sys.modules["bench_score_m2"] = mod
    spec.loader.exec_module(mod)
    return mod


def _ltr(identity=95.0):
    return [{"type": "LTR (terminal direct repeat)", "five_prime": [0, 200], "three_prime": [2800, 3000],
             "element_span": [0, 3000], "ltr_len": 200, "identity": identity}]


def _tir():
    return [{"type": "TIR (terminal inverted repeat)", "five_prime": [0, 40], "three_prime": [960, 1000],
             "element_span": [0, 1000], "tir_len": 40, "identity": 95.0}]


def _dom(code, a, b, hmm="", cls=""):
    return {"domain": code, "hmm": hmm, "class": cls, "nt": [a, b], "strand": "+", "score": 100.0, "orf": 0}


def _matrix():
    from teagle_core import classify
    c = classify.classify
    return {
        "copia": c(_ltr(), [_dom("INT", 500, 900), _dom("RT", 1000, 1500)]),
        "gypsy": c(_ltr(), [_dom("RT", 500, 900), _dom("INT", 1000, 1500)]),
        "dirs": c(_ltr(), [_dom("RT", 300, 600), _dom("YR", 650, 800)]),
        "cacta": c([], [_dom("TPase", 100, 900, cls="dna:CACTA")]),
        "mule": c([], [_dom("TPase", 100, 900, cls="dna:MULE")]),
        "is4": c([], [_dom("TPase", 100, 900, cls="dna:IS4")]),
        "dde": c([], [_dom("TPase", 100, 900, cls="dna:DDE")]),
        "generic_dna": c([], [_dom("TPase", 100, 900)]),
        "non_autonomous": c(_tir(), []),
        "ltr_partial": c(_ltr(), [_dom("GAG", 1000, 1400, hmm="Gag_p24"), _dom("PR", 1500, 1900, hmm="RVP")]),
        "ltr_unclassified": c(_ltr(), [_dom("RT", 500, 900), _dom("RNaseH", 950, 1200)]),
        "erv": c(_ltr(), [_dom("GAG", 500, 900), _dom("RT", 1000, 1500), _dom("INT", 1600, 2000), _dom("ENV", 2100, 2500)]),
        "tdr": c(_ltr(), []),
        "lone_rt": c([], [_dom("RT", 100, 400)]),
        "retro_unassigned": c([], [_dom("RT", 100, 300), _dom("PR", 320, 450, hmm="RVP")], structural_ok=False),
        "helitron": c([], [_dom("HEL", 100, 500)]),
        "none": c([], []),
    }


def test_rank_class_and_order_agree_with_the_benchmark_scorer():
    import main
    score = _score_module()
    for name, cl in _matrix().items():
        view = main._rank_view(cl)
        c, order, _abstain = score.observed(cl)
        assert (view is None) == (c == "NO_CALL"), name
        if view:
            assert view["class"] == {"CLASS_I": "Class I", "CLASS_II": "Class II"}[c], name
            if order in ("LTR", "LINE", "SINE", "TIR"):
                assert view["order"] == order, name
            else:                                                      # ABSTAIN / OTHER: no order is named, except DIRS, which is itself an order
                assert view["order"] == ("DIRS" if cl["te_class"] == "DIRS" else None), name


def test_a_superfamily_is_named_only_where_the_scorer_commits_to_one():
    import main
    score = _score_module()
    order_only, owner_decision = {"DDE_UNPLACED", "DNA_UNPLACED"}, {"PIGGYBAC", "DIRS"}
    exercised = set()
    for name, cl in _matrix().items():
        view = main._rank_view(cl)
        if not view:
            continue
        try:
            token = score.observed_superfamily(cl)
        except ValueError:                                         # a label the frozen scorer cannot read (unmapped refusals): nothing to compare
            assert view["superfamily"] is None, name
            continue
        exercised.add(token)
        if token in score.REPRESENTABLE - order_only - owner_decision and not view["erv"]:
            assert view["superfamily"], (name, token)
        else:
            assert view["superfamily"] is None, (name, token)          # order-only and owner-decision tokens: the display is more conservative
    assert exercised & order_only == order_only and exercised & owner_decision == owner_decision


def test_wicker_table_covers_the_te_classes_the_classifier_emits():
    import main
    for name, cl in _matrix().items():
        te = cl["te_class"]
        if main._rank_view(cl) and not cl.get("is_erv"):
            assert te in main._WICKER, (name, te)


def test_lone_rt_tdr_and_helitron_get_no_rank_line_and_the_erv_line_says_why():
    import main
    m = _matrix()
    for name in ("tdr", "lone_rt", "helitron", "none"):
        assert main._rank_view(m[name]) is None, name
    erv = main._rank_line(main._rank_view(m["erv"]))
    assert erv.startswith("Wicker 2007: Class I › LTR › superfamily not resolved — env-bearing (RLX)")
    assert "(RXX)" in main._rank_line(main._rank_view(m["retro_unassigned"]))


# ================================================================ published-as, scope caution, explanation (B19, B5, B13, B10, U2, U36)
def _notes_with(win, needle):
    return [l for l in win.card_struct.findChildren(QLabel) if needle in l.text()]


def _flanked(name, left="AC" * 150, right="AG" * 150):
    import engine
    seq = _analysis(name)["records"][0]["seq"]
    return engine.analyze(f">flanked\n{left}{seq}{right}")


def test_bundled_example_carries_its_published_class_and_the_verdict_in_words(win):
    import main
    from teagle_core import examples
    _accept(win, "M11240")
    note = _notes_with(win, "Published as:")[0]
    assert "Published as: LTR retrotransposon · Copia (Ty1)" in note.text() and "is consistent with the published class" in note.text()
    assert note.property("level") == "success" and "not whether the published element is complete" in note.text()
    rec = copy.deepcopy(win.state["records"][0])
    rec["classification"]["te_class"] = "LTR/Gypsy"
    win._render_struct_card(rec, {})
    note = _notes_with(win, "Published as:")[0]
    assert "differs from the published class (LTR/Copia)" in note.text() and note.property("level") == "warn"
    rec["classification"]["outcome"] = "withheld"
    win._render_struct_card(rec, {})
    note = _notes_with(win, "Published as:")[0]
    assert "not compared" in note.text() and note.property("level") == "info"
    assert set(main._PUBLISHED_CLASS) == {e[0] for e in examples.CATALOG}


def test_published_class_constants_equal_the_golden_calls():
    import main
    for acc, code in main._PUBLISHED_CLASS.items():
        assert _analysis(acc)["records"][0]["classification"]["te_class"] == code


def test_an_edited_or_reverse_complemented_example_is_not_the_example(win):
    import main
    from teagle_core import sequtil
    seq = _analysis("M11240")["records"][0]["seq"]
    assert main._example_match(seq) and not main._example_match(seq[:-1] + ("A" if seq[-1] != "A" else "C"))
    assert not main._example_match(sequtil.reverse_complement(seq)) and not main._example_match("")
    _accept(win, "M11240_edit", result=_edited_copia())
    assert not _notes_with(win, "Published as:")


def _edited_copia():
    import engine
    seq = _analysis("M11240")["records"][0]["seq"]
    return engine.analyze(">e\n" + seq[:2000] + ("A" if seq[2000] != "A" else "C") + seq[2001:])


def test_large_input_confirmation_gates_the_run(win, monkeypatch):
    import main
    calls, submits = [], []
    answer = {"v": "cancel"}
    monkeypatch.setattr(win, "_confirm_large_input", lambda *a: (calls.append(a), answer["v"])[1])
    win.engine.submit = lambda op, body=None, key=None: submits.append(op)
    seq = ("ACGT" * 60000)[:main.LARGE_INPUT_BP + 1]
    win.seq.setPlainText(">big\n" + seq)
    gen = win._spec_gen
    win._run_analysis()
    assert calls == [(main.LARGE_INPUT_BP + 1, 1, 1)] and submits == [] and win.runBtn.isEnabled()
    assert win._spec_gen == gen                                       # Cancel invalidated nothing
    answer["v"] = "panel03"
    win._run_analysis()
    assert submits == [] and not win.card_wsl.body.isHidden() and win.wslSource.currentIndex() == 0
    answer["v"] = "run"
    win._run_analysis()
    assert submits == ["analyze"]


def test_exactly_the_threshold_asks_nothing_and_the_second_record_is_named(win, monkeypatch):
    import main
    calls, submits = [], []
    monkeypatch.setattr(win, "_confirm_large_input", lambda *a: (calls.append(a), "cancel")[1])
    win.engine.submit = lambda op, body=None, key=None: submits.append(op)
    win.seq.setPlainText(">r\n" + ("ACGT" * 60000)[:main.LARGE_INPUT_BP])
    win._run_analysis()
    assert calls == [] and submits == ["analyze"]
    win.seq.setPlainText(">a\nACGT\n>b\n" + ("ACGT" * 60000)[:main.LARGE_INPUT_BP + 5])
    win._analyze_inflight = False; win._run_analysis()
    assert calls[-1][1:] == (2, 2)


def test_the_engine_large_input_note_follows_the_named_constant(monkeypatch):
    import engine, main
    assert main.LARGE_INPUT_BP == engine.LARGE_INPUT_BP == 200_000
    monkeypatch.setattr(engine, "LARGE_INPUT_BP", 1000)        # same rule, fast input: the note keys off the constant
    for n, expect in ((1001, True), (1000, False)):
        rec = engine.analyze(">r\n" + ("ACGT" * 300)[:n])["records"][0]
        assert any("large input" in x for x in rec["notes"]) == expect, n


def test_caution_strip_on_a_tandem_array(win):
    rec = _accept(win, "array230")
    scan = rec["processing"]["domains"]
    tied = sum(1 for o in rec["orfs"][scan["orfsScanned"]:] if o["length_aa"] == rec["orfs"][scan["orfsScanned"] - 1]["length_aa"])
    banner = next(f for f in win.card_struct.findChildren(QFrame) if f.objectName() == "classbn")
    strip = next(l for l in banner.findChildren(QLabel) if l.accessibleName() == "Scope caution")
    assert strip.property("level") == "warn" and strip.text().startswith("<b>Scope caution</b>")
    text = strip.text()
    assert "Large input (230,000 bp)" in text and "Multiple copies likely: 12" in text and "INT + RT" in text
    assert f"{tied} unsearched ORF(s) are exactly as long as the shortest one searched" in text and tied > 0
    assert "TEagle display heuristic" in text and "the 479 unsearched ORFs may hold more" in text
    chip = next(c for c in _by_name(win.card_struct, "cf") if c.text().startswith("CALL CONFIDENCE"))
    assert chip.text().endswith("SEE CAUTION") and chip.toolTip().startswith("Scope caution below")
    lay = banner.layout()
    idx = {lay.itemAt(i).widget(): i for i in range(lay.count()) if lay.itemAt(i).widget()}
    expl = next(l for l in banner.findChildren(QLabel) if l.objectName() == "classexp")
    assert idx[strip] < idx[expl]


@pytest.mark.parametrize("name", ["M11240", "M12927", "M80343", "X01005", "X05424"])
def test_no_caution_on_a_single_reference(win, name):
    _accept(win, name)
    assert not [l for l in win.card_struct.findChildren(QLabel) if l.accessibleName() == "Scope caution"]
    assert not any("SEE CAUTION" in c for c in _chips(win))


def test_a_copia_dimer_is_flagged_as_more_than_one_copy(win):
    _accept(win, "copia2")
    assert any("Multiple copies likely: 2" in l.text() for l in win.card_struct.findChildren(QLabel))


def test_the_orf_budget_note_is_reworded_only_when_the_budget_cut_a_tie(win):
    rec = _accept(win, "array230")
    tie = __import__("main")._orf_budget_tie(rec)
    notes = [l.text() for l in win.card_struct.findChildren(QLabel) if l.text().startswith("• protein-domain search")]
    assert notes and f"{tie[2]} further ORF(s) were not searched — {tie[3]} of them exactly as long as the shortest searched ORF" in notes[0]
    _accept(win, "M12927")
    notes = [l.text() for l in win.card_struct.findChildren(QLabel) if l.text().startswith("• protein-domain search")]
    assert notes and "12 shorter ORF(s) were not searched" in notes[0]


def test_explanation_rewrites_leave_the_stored_text_alone(win):
    import main
    for name, needles in (("M11240", ["Reverse-transcriptase domain present", "Copia/Ty1-type domain order"]),
                          ("M80343", ["ORF1 ×3 (Tnp_22_trimer, Transposase_22, Tnp_22_dsRBD) · EN · RT."]),
                          ("X01005", ["Domains: TPase ×2 (HTH_Tnp_Tc3_2, DDE_3)."]),
                          ("array230", ["Domains: INT ×12 · RT ×12.", "further ORF(s) were outside"])):
        rec = copy.deepcopy(_analysis(name)["records"][0])
        raw = rec["classification"]["explanation"]
        view = main._explanation_view(rec)
        for n in needles:
            assert n in view, (name, n)
        assert rec["classification"]["explanation"] == raw
        assert "ORF1 · ORF1" not in view and "TPase · TPase" not in view and "INT · RT · INT" not in view
    _accept(win, "M80343")
    assert "ORF1 ×3" in _texts(win.card_struct)


def test_every_explanation_rewrite_rule_fires_and_a_drifted_phrase_is_left_alone():
    import main
    fired = set()
    for name in ("M11240", "M12927", "M80343", "X01005", "array230", "iupac", "sample"):
        rec = _analysis(name)["records"][0]
        raw, view = rec["classification"]["explanation"], main._explanation_view(rec)
        if "×" in view:
            fired.add("collapse")
        if "not assessable (no ORF" in view and "none detected in tested panel" in raw:
            fired.add("no-orf")
        if "Reverse-transcriptase domain present" in view and "reverse-transcriptase domain present" in raw:
            fired.add("capital")
        if "Structural: TDR" in view and "Structural: LTR" in raw:
            fired.add("tdr")
        if "further ORF(s) were outside" in view and "shorter ORF(s) were outside" in raw:
            fired.add("further")
    assert fired == {"collapse", "no-orf", "capital", "tdr", "further"}
    rec = copy.deepcopy(_analysis("array230")["records"][0])
    rec["classification"]["explanation"] = rec["classification"]["explanation"].replace("shorter ORF(s) were outside", "fewer ORFs were outside")
    assert "further" not in main._explanation_view(rec)                # the rule needs its source phrase; without it the text is verbatim


def test_coding_sentence_counts_orfs_and_names_the_multi_copy_reading(win):
    import main
    rec = _analysis("array230")["records"][0]
    note = main._coding_view(rec)
    assert "12 separate ORFs" in note and "12 of these ORFs carry the same INT + RT domain set" in note
    assert "degenerate" not in note and note[0].isupper()
    line = main._coding_view(_analysis("M80343")["records"][0])
    assert "2 separate ORFs" in line and "canonical LINE ORF1/ORF2 organisation" in line and "degenerate" not in line
    assert "reading frame(s)" not in line


def test_reference_sentence_same_record_and_a_flanked_copy(win):
    _accept(win, "M11240")
    text = _texts(win.card_struct)
    assert "This record is the bundled copia reference (M11240)" in text and "even for this published element" in text
    _accept(win, "flanked", result=_flanked("M11240"))
    text = _texts(win.card_struct)
    assert "The bundled copia reference (M11240, Drosophila melanogaster) gives the same result" in text


def test_no_reference_sentence_when_the_reference_reads_intact(win):
    _accept(win, "M80343")
    assert "bundled LINE-1" not in _texts(win.card_struct) and "gives the same result" not in _texts(win.card_struct)


def test_a_failed_domain_scan_gives_no_reference_sentence(win, monkeypatch):
    import engine, inspect
    from teagle_core import domains

    def boom(*a, **k):
        raise RuntimeError("forced")
    boom.__signature__ = inspect.signature(domains.scan_domains)         # the seal reads the detector defaults
    with monkeypatch.context() as m:
        m.setattr(domains, "scan_domains", boom)
        res = engine.analyze(">c\n" + _analysis("M11240")["records"][0]["seq"])
    _accept(win, "failed", result=res)
    text = _texts(win.card_struct)
    assert "did not run" in text and "Profiles tested" not in text and "gives the same result" not in text


# ================================================================ visible row actions (B21, U20)
def _bar_buttons(root, text):
    return [b for b in root.findChildren(QPushButton) if b.text().upper() == text.upper()]


def _capture_menus(monkeypatch):
    """Replace the table's menu class by one that records its items instead of running a modal loop."""
    import widgets
    seen = []

    class _Menu(widgets.QMenu):
        def exec(self, *a):
            seen.append([x.text() for x in self.actions() if x.text()])
    monkeypatch.setattr(widgets, "QMenu", _Menu)
    return seen


def _primer_win(win):
    _accept(win, "M11240")
    win._design_tmpl = "ACGT" * 60
    win._render_primers({"candidates": [_cand(1), _cand(2)], "oligoqc_engines": {}})
    from widgets import DataTable
    return next(t for t in win.primBody.parentWidget().findChildren(DataTable) if t.export_base == "TEagle_primers")


def test_primer_row_actions_match_the_context_menu(win, monkeypatch):
    t = _primer_win(win)
    sent = []
    monkeypatch.setattr(win, "_add_pcr_pair", lambda c: sent.append(c))
    assert t.current_orig() == 0                                   # row 0 is current as soon as the table is built
    t.select_orig(1)
    send = _bar_buttons(win.primBody.parentWidget(), "Send to in-silico PCR")[0]
    t.current_row_changed.emit(1)                                  # what a click on the row does
    assert send.isEnabled()
    send.click()
    assert sent == [win.state["candidates"][1]]
    label = next(f for l, f in t.row_items(1) if l == "→ send to in-silico PCR")
    label()
    assert sent[1] is sent[0]
    lab = [l.text() for l in win.primBody.parentWidget().findChildren(QLabel) if l.text().startswith("Pair ")]
    assert lab and lab[0].endswith("P2:")


def test_shift_f10_and_the_menu_key_open_the_current_row_menu(win, monkeypatch):
    t = _primer_win(win)
    seen = _capture_menus(monkeypatch)
    win.show(); t.setFocus(); _app.processEvents()
    QTest.keyClick(t, Qt.Key_F10, Qt.ShiftModifier)
    QTest.keyClick(t, Qt.Key_Menu)
    assert len(seen) == 2
    for items in seen:
        assert "→ send to in-silico PCR" in items and "Copy row" in items
    QTest.keyClick(t, Qt.Key_F10)                                  # without Shift it is just a key
    assert len(seen) == 2


def test_row_actions_button_pops_exactly_the_right_click_menu(win, monkeypatch):
    t = _primer_win(win)
    seen = _capture_menus(monkeypatch)
    btn = next(b for b in win.primBody.parentWidget().findChildren(QPushButton) if b.accessibleName() == "Actions for the selected row")
    btn.click()
    t._menu(QPoint(5, t.rowViewportPosition(0) + 5))
    assert len(seen) == 2 and seen[0] == seen[1]


def test_domain_and_orf_rows_get_a_design_button_that_routes_like_the_menu(win, monkeypatch):
    rec = _accept(win, "M11240")
    calls = []
    monkeypatch.setattr(win, "_design_for_domain", lambda *a, **k: calls.append((a, tuple(sorted(k.items())))))
    from widgets import DataTable
    dom = next(t for t in win.card_struct.findChildren(DataTable) if str(t.export_base).endswith("_domains"))
    buttons = _bar_buttons(win.card_struct, "Design primers here")
    assert len(buttons) == 2                                       # the ORF table and the domain table each get one
    dom_btn = buttons[-1]                                          # the ORF bar is built first, the domain bar last
    assert dom_btn.isEnabled()
    dom_btn.click()
    first = calls[-1]
    next(f for l, f in dom.row_items(dom.current_orig()) if l == "Design primer here")()
    assert calls[-1] == first


def test_row_actions_are_disabled_after_invalidation(win):
    t = _primer_win(win)
    controls = [b for b in win.primBody.parentWidget().findChildren(QPushButton) if "ACTIONS" in b.text().upper() or "SEND TO" in b.text().upper()]
    assert len(controls) == 2 and all(b.isEnabled() for b in controls)
    win._invalidate_analysis("x")
    assert not any(b.isEnabled() for b in controls)


def test_structural_rows_get_actions_but_no_primary(win):
    _accept(win, "M11240")
    from widgets import DataTable
    t = next(t for t in win.card_struct.findChildren(DataTable) if str(t.export_base).endswith("_structural"))
    assert t.current_orig() == 0
    names = [b.accessibleName() for b in t.parentWidget().findChildren(QPushButton)]
    assert "Actions for the selected row" in names
    lab = [l.text() for l in win.card_struct.findChildren(QLabel) if l.text().startswith(("LTR ", "TDR "))]
    assert lab and lab[0].endswith(":")


def test_the_hints_no_longer_claim_right_click_is_the_only_way(win):
    _primer_win(win)
    text = _texts(win.primBody.parentWidget())
    assert "Row actions ▾" in text and "Right-click a pair" not in text and "collapse the specimen panel" not in text


# ================================================================ quick start, next steps, scroll, help, glossary (B20, U23)
def test_quickstart_examples_load_and_run(fresh_win, monkeypatch):
    from teagle_core import examples, sequtil
    win = fresh_win
    assert not win._quickstart.isHidden()
    links = [b.text() for b in win._quickstart_links]
    assert links == [e[1] for e in examples.available()] + ["negative control"]
    submits = []
    win.engine.submit = lambda op, body=None, key=None: submits.append(op)
    win._quickstart_links[0].click()
    first = examples.available()[0][0]
    assert sequtil.parse_fasta(win.seq.toPlainText())[0][1] == sequtil.parse_fasta(examples.load(first))[0][1]
    assert submits[-1] == "analyze"
    win._analyze_inflight = False
    win._quickstart_links[-1].click()
    assert submits[-1] == "analyze" and win.seq.toPlainText().startswith(">")


def test_quickstart_hide_persists_and_help_reopens(monkeypatch):
    FakeSettings.store = {}
    w = _make_win(monkeypatch)
    try:
        hide = next(b for b in w._quickstart.findChildren(QPushButton) if b.accessibleName() == "Hide the quick-start guide")
        hide.click()
        assert w._quickstart.isHidden() and FakeSettings.store["quickstart_hidden"] is True
    finally:
        _dispose(w)
    w2 = _make_win(monkeypatch)
    try:
        assert w2._quickstart.isHidden()
        w2._show_quickstart()
        assert not w2._quickstart.isHidden() and FakeSettings.store["quickstart_hidden"] is False
    finally:
        _dispose(w2)


def test_a_new_analysis_scrolls_to_its_card_and_a_reselection_does_not(win, monkeypatch):
    from PySide6.QtCore import QPoint
    win.resize(1366, 768); win.show(); _app.processEvents()
    _accept(win, "M11240")
    bar = win.resultsScroll.verticalScrollBar()
    bar.setValue(bar.maximum())
    _accept(win, "X01005")
    inner = win.resultsScroll.widget()
    assert _spin(lambda: bar.value() == win.card_struct.mapTo(inner, QPoint(0, 0)).y())
    calls = []
    monkeypatch.setattr(win, "_scroll_to_top", lambda w: calls.append(w))
    _accept(win, "multi")
    assert len(calls) == 1
    win._show_record(1, keep_table=True)
    _spin(lambda: False, 0.3)
    assert len(calls) == 1                                          # re-selecting a record never moves the page


def test_next_steps_gate_on_flanks(win, monkeypatch):
    import engine, helpers
    from teagle_core import sequtil
    asked = []
    monkeypatch.setattr(win, "_flank_picker", lambda *a, **k: asked.append(a))
    flanked = engine.analyze(">f\n" + _analysis("M11240")["records"][0]["seq"].join(["AC" * 200, "AG" * 200]))
    rec = _accept(win, "flanked", result=flanked)
    span = next(e for e in rec["structural"] if e["type"].startswith("LTR"))["element_span"]
    link = next(b for b in win.card_struct.findChildren(QPushButton) if b.text().startswith("Flanking sequence around the LTR element"))
    link.click()
    assert asked == [(span[0], span[1], "+", "LTR_element")]
    _accept(win, "M11240")
    assert not [b for b in win.card_struct.findChildren(QPushButton) if b.text().startswith("Flanking sequence around")]
    assert "0 bp before and 0 bp after" in _texts(win.card_struct) and "Fetch by coordinate" in _texts(win.card_struct)
    assert [b for b in win.card_struct.findChildren(QPushButton) if b.text().startswith("Design primers for this record")]


def test_genotype_button_needs_a_flank_long_enough_for_a_primer(win):
    from teagle_core import primers
    geno = lambda: [b for b in win.card_struct.findChildren(QPushButton) if b.text().upper().startswith("GENOTYPE THIS INSERTION")]
    _accept(win, "M11240")                                          # the record ends at the element
    assert not geno() and "Fetch by coordinate" in _texts(win.card_struct)
    short = "AC" * (primers.MIN_PRIMER_SIZE // 2 - 1)               # one base under the minimum primer length, both sides
    _accept(win, "flanked", result=_flanked("M11240", left=short, right=short))
    assert not geno() and "Fetch by coordinate" in _texts(win.card_struct)
    _accept(win, "flanked", result=_flanked("M11240", left="AC" * primers.MIN_PRIMER_SIZE, right=""))
    assert geno() and "Fetch by coordinate" in _texts(win.card_struct)   # a primer fits, but not a full product range
    _accept(win, "flanked", result=_flanked("M11240"))
    assert geno() and not [l for l in _labels(win.card_struct) if "Fetch by coordinate" in l]


def test_design_next_step_focuses_the_design_button(win):
    _accept(win, "M11240")
    win.show(); _app.processEvents()
    seen = []
    real = win._scroll_to
    win._scroll_to = lambda w: (seen.append(w), real(w))[1]
    next(b for b in win.card_struct.findChildren(QPushButton) if b.text().startswith("Design primers for this record")).click()
    assert seen == [win.card_primer] and not win.card_primer.body.isHidden()
    assert win.designBtn.hasFocus() or QApplication.focusWidget() is win.designBtn


def test_primer_next_step_stages_all_pairs(win):
    _primer_win(win)
    next(b for b in win.primBody.parentWidget().findChildren(QPushButton) if b.text().startswith("Stage all pairs")).click()
    assert len(win.state["pcrPairs"]) == 2


def test_next_step_buttons_are_disabled_after_invalidation(win):
    _accept(win, "M11240")
    steps = [b for b in win.card_struct.findChildren(QPushButton) if "→" in b.text()]
    assert steps and all(b.isEnabled() for b in steps)
    win._invalidate_analysis("x")
    assert not any(b.isEnabled() for b in steps)


def test_glossary_content_and_f1(win):
    import main
    from teagle_core import domains
    win._show_glossary()
    text = win._glossary_dlg._body.text()
    assert win._glossary_dlg.isVisible() or not win._glossary_dlg.isHidden()
    assert main.GLOSS["Hit confidence"].split("(")[0].strip() in text.replace("&amp;", "&") or "Hit confidence" in text
    assert str(domains.MAX_ORFS_SCANNED) in text and "Call status" in text and "0-based" in text
    win._glossary_dlg.hide()
    win.show(); _app.processEvents()
    QTest.keyClick(win, Qt.Key_F1)
    assert not win._glossary_dlg.isHidden()


def test_manual_opens_the_versioned_release_page(win, monkeypatch):
    import main
    got = []
    monkeypatch.setattr(main, "open_external_url", lambda u: got.append(u) or True)
    win._open_manual()
    assert got == [f"https://github.com/tunabirgun/TEagle/releases/tag/v{main.APP_VERSION}"]


# ================================================================ casing, notices, numbers, theme re-tint, labels (B30, B18, U8, U32, U35, B11, B22)
def test_caps_keeps_units_and_is_idempotent():
    import widgets
    assert widgets.caps("Prod min (bp)") == "PROD MIN (bp)" and widgets.caps("Opt Tm (°C)") == "OPT Tm (°C)"
    assert widgets.caps("ORFs (≥40 aa) — 11") == "ORFS (≥40 aa) — 11"
    assert widgets.caps("Organism — local isPcr") == "ORGANISM — LOCAL isPcr"
    for s in ("Prod min (bp)", "Run analysis", "◴ designing…"):
        assert widgets.caps(widgets.caps(s)) == widgets.caps(s)


def _action_buttons(win):
    import widgets
    return [b for b in win.findChildren(QPushButton) if widgets._is_action_button(b) and b.text()]


def test_action_buttons_keep_one_casing_through_state_changes(win, monkeypatch):
    import widgets
    from sample import make_sample

    def check(stage):
        bad = [b.text() for b in _action_buttons(win) if b.text() != widgets.caps(b.text())]
        assert not bad, (stage, bad)
    check("construction")
    win._set_seq(make_sample()); check("sample loaded")
    submits = []
    win.engine.submit = lambda op, body=None, key=None: submits.append(op)
    win._run_analysis(); check("analysis in flight")
    win._invalidate_analysis("x", cancel_request=True); check("cancelled")
    _accept(win, "M11240"); check("accepted")
    win._toggle_rail(); win._toggle_rail(); check("rail toggled twice")
    win._design_inflight = False
    win._design(); check("design in flight")
    win._on_primers({"candidates": [_cand(1)]}); check("design finished")
    win._add_pcr_pair(win.state["candidates"][0])
    win._run_pcr(); check("pcr in flight")


@pytest.mark.parametrize("handler", ["_fetch", "_run_analysis", "_fetch_coord_empty", "_design", "_run_pcr", "_annotate",
                                      "_splice", "_export_manifest", "_export_analysis_result"])
def test_precondition_notices_are_info_and_sentence_case(win, handler):
    if handler == "_fetch_coord_empty":
        win._fetch_coord()
    else:
        getattr(win, handler)()
    notif = win._notif
    assert notif.property("level") == "info"
    text = notif.findChildren(QLabel)[0].text()
    assert text[0].isupper() or text[0].isdigit() or text[0] in "“\"(", text
    assert text.rstrip()[-1] in ".)?", text


def test_counts_use_separators_and_table_cells_do_not(win):
    import widgets
    win.seq.setPlainText("A" * 230000)
    assert win.charCount.text() == "230,000 nt" or win.charCount.text().startswith("230,000")
    rows = win._record_summary_rows([{"id": "r", "valid": True, "composition": {"length": 5146, "gc": 50.0, "n": 0.0,
                                                                                  "ambiguous": 0.0}, "classification": None,
                                      "processing": {}}])
    assert rows[0][win._RECORD_HEADERS.index("Length")] == 5146
    assert widgets._sortkey("1,500") == 1500.0
    t = widgets.DataTable(["n"]); t.set_rows([["1,500"], ["999"]]); t.sortByColumn(0, Qt.AscendingOrder)
    assert [t.item(i, 0).text() for i in range(2)] == ["999", "1,500"]
    rec = _accept(win, "M11240")
    assert win.rOrf.text().count(",") == 0 and "11" in win.rOrf.text()


def test_theme_toggle_retints_rich_text_and_table_ink(win):
    import theme as theme_mod
    _accept(win, "M11240")
    win._render_primers({"candidates": [dict(_cand(1), qc={"ok": True, "worst": "warn", "left": {"hairpin": {"p3": -10.0, "vrna": -10.0, "flag": "warn", "agree": "agree"}, "self_dimer": {}},
                                                       "right": {"hairpin": {}, "self_dimer": {}}, "hetero_dimer": {}, "end_stability": {}})],
                         "oligoqc_engines": {}})
    old, new = win.theme, "dark"
    acc_old, acc_new = theme_mod.ACCENT[old], theme_mod.ACCENT[new]
    win.rRecords.setText(win._accent_link("https://example.org", "X"))
    assert f"color:{acc_old}" in win.rRecords.text()
    from widgets import DataTable
    extra = DataTable(["x"]); extra.set_rows([["a"]], styles=[[theme_mod.BAD[old]]])    # a table the toggle does not re-render
    win.card_prov.bodylay.addWidget(extra)
    win._toggle_theme()
    assert extra.item(0, 0).foreground().color().name().upper() == theme_mod.BAD[new].upper()
    assert f"color:{acc_old}" not in _texts(win) and f"color:{acc_new}" in win.rRecords.text()
    win._toggle_theme()
    assert extra.item(0, 0).foreground().color().name().upper() == theme_mod.BAD[old].upper()
    assert f"color:{acc_new}" not in _texts(win) and f"color:{acc_old}" in win.rRecords.text()


def test_font_floor_and_names_for_the_new_controls(win):
    for b in (win.railToggle, win.backendBtn, win.scaleBtn, win.themeBtn, win.helpBtn):
        assert len(b.accessibleName() or b.text()) >= 2
    assert win.mark.accessibleName() == "TEagle"
    assert win.backendBtn.text().upper() == "BACKEND"             # tests/test_wsl_install_progress.py pins this spelling


def test_the_iupac_tile_keeps_its_verdict_and_puts_the_share_on_a_second_line(win):
    _accept(win, "iupac_valid", result=_iupac_valid())
    assert win.mValid.text().startswith("valid<br>") and "% ambiguous" in win.mValid.text()
    _accept(win, "M11240")
    assert win.mValid.text() == "valid"


def _iupac_valid():
    import engine
    return engine.analyze(">m\n" + ("ACGTTGCARGGTCATCGATC" * 20))


def test_rail_caption_waits_for_a_result_and_the_footnote_has_no_licence_tag(win):
    assert win.rStructCap.isHidden()
    _accept(win, "M11240")
    assert not win.rStructCap.isHidden()
    win._invalidate_analysis("x")
    assert win.rStructCap.isHidden()
    assert "CC0" not in _texts(win.rail)


def test_the_pcr_mechanism_cell_is_readable_and_the_key_is_kept_for_export(win):
    from widgets import DataTable
    _accept(win, "M11240")
    amp = {"pair": "P1", "mechanism_key": "record_abc(input):F+R:4591:4766:175", "source": "record_abc (input)",
           "start": 4591, "end": 4766, "length": 175, "product_length": 175, "primer_orientation": "F+R", "left_primer": "F",
           "right_primer": "R", "fwd_mm": 0, "rev_mm": 0, "fwd_tail5": 0, "rev_tail5": 0, "fwd_anneal_len": 20, "rev_anneal_len": 20,
           "on_target": None, "matches_intended_product": None, "single_primer": False}
    win._render_pcr([{"label": "P1", "amplicons": [amp]}], [amp])
    t = next(t for t in win.pcrBody.parentWidget().findChildren(DataTable) if t.export_base == "TEagle_amplicons")
    col = [t.horizontalHeaderItem(j).text() for j in range(t.columnCount())].index("Mechanism")
    assert t.item(0, col).text() == "F + R, 4591–4766" and "record_abc" in t.item(0, col).toolTip()
    assert t.rows_data()[0][col] == "record_abc(input):F+R:4591:4766:175"


def test_the_pcr_queue_shows_whole_primers(win):
    c = dict(_cand(1), left_seq="ACGTACGTACGTACGTACGTACGTAC", right_seq="TGCATGCATGCATGCATGCATGCATG")
    win._add_pcr_pair(c)
    text = _texts(win.pcrQueueBox.parentWidget())
    assert c["left_seq"] in text and c["right_seq"] in text and "…" not in text
    names = [b.accessibleName() for b in win.pcrQueueBox.parentWidget().findChildren(QPushButton)]
    assert {"Move this pair up", "Move this pair down", "Remove this pair from the queue"} <= set(names)


def test_coordinate_conventions_are_stated_on_the_tables_that_use_them(win):
    _accept(win, "M80343")
    text = _texts(win.card_struct)
    assert "Start and End are 0-based and half-open" in text and "Coords are 0-based and half-open" in text
    assert "Detection methods" in text


def test_ltr_metric_shows_termini_and_the_poly_a_row_states_an_exact_run(win):
    from widgets import DataTable

    def metrics():
        t = next(t for t in win.card_struct.findChildren(DataTable) if str(t.export_base).endswith("_structural"))
        return {t.item(r, 0).text(): t.item(r, 3).text() for r in range(t.rowCount())}
    _accept(win, "M11240")
    ltr = next(v for k, v in metrics().items() if k.startswith("LTR"))
    assert "TG…CA (canonical)" in ltr and "lower bound" in ltr
    _accept(win, "M12927")
    assert "AG…TT (no known motif)" in next(v for k, v in metrics().items() if k.startswith("LTR"))
    _accept(win, "M80343")
    assert next(v for k, v in metrics().items() if k.startswith("poly-A")) == "uninterrupted A run"
