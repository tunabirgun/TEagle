"""4.1 UI lane M3: assay goals, the primer table (location, notes, order sheet), coordinate display and on-screen record
labels. Drives the real MainWindow offscreen; the registry is never written (QSettings is an in-memory store) and no
request leaves the process."""
import copy, csv, functools, io, os, sys

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
import pytest

pytest.importorskip("PySide6")
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_APP = os.environ.get("M3_APP_ROOT") or os.path.join(_ROOT, "app")      # a pre-change copy can be pointed at for failing-before runs
sys.path.insert(0, os.path.join(_APP, "native"))
sys.path.insert(0, os.path.join(_APP, "backend"))
sys.path.insert(0, os.path.join(_ROOT, "tests"))

from PySide6.QtCore import QEvent
from PySide6.QtWidgets import QApplication, QLabel, QPushButton, QMessageBox, QDialog

_app = QApplication.instance() or QApplication([])


class FakeSettings:
    store = {}

    def __init__(self, *a):
        pass

    def value(self, key, default=None):
        return self.store.get(key, default)

    def setValue(self, key, value):
        self.store[key] = value


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


@functools.lru_cache(maxsize=None)
def _analysis(name):
    import engine, helpers
    from teagle_core import examples, sequtil
    if name in ("M11240", "M12927", "M80343", "X01005", "X05424"):
        return engine.analyze(examples.load(name))
    if name == "padded":                                   # an LTR element with 600 bp of genome on each side
        return engine.analyze(">locus\n" + helpers.make_ltr_element(seed=3, flank=600))
    if name == "labelled":
        return engine.analyze(_labelled_text())
    raise KeyError(name)


@functools.lru_cache(maxsize=None)
def _labelled_text():
    from teagle_core import examples, sequtil
    recs = [sequtil.parse_fasta(examples.load(n))[0][1] for n in ("M11240", "X01005")]
    return f">patient_042 fly stock\n{recs[0]}\n>C:\\Users\\person\\x.fa\n{recs[1]}\n"


def _accept(win, name, result=None, sequence=None):
    """Hand the window a finished analysis as _on_analyze receives it, with the request text the engine saw."""
    result = copy.deepcopy(result or _analysis(name))
    if sequence is None:
        sequence = "\n".join(f">{r['id']}\n{r['seq']}" for r in result.get("records", []) if r.get("seq"))
    win.seq.setPlainText(sequence)
    win._analyze_gen += 1
    gen = win._analyze_gen
    win._active_analyze_gen = gen
    win._analysis_requests[gen] = {"sequence": sequence, "source": {}}
    win._analyze_inflight = True
    win._on_analyze(result, gen)
    return win.state["records"][win.state.get("record_index", 0)]


def _fasta_for(name):
    from teagle_core import examples
    return examples.load(name)


def _tables(root):
    from widgets import DataTable
    return {t.export_base: t for t in root.findChildren(DataTable)}


def _struct_table(win):
    return next(t for t in _tables(win.card_struct).values() if str(t.export_base).endswith("_structural"))


def _table(win, suffix):
    return next(t for t in _tables(win.card_struct).values() if str(t.export_base).endswith(suffix))


def _label_text(root):
    return " ".join(l.text() for l in root.findChildren(QLabel))


def _capture(win, ops=("primers", "pcr", "annotate", "splice", "fetch", "genome_pcr")):
    calls = []
    real = win.engine.submit

    def submit(op, body=None, key=None):
        if op in ops:
            calls.append((op, body, key or op))
        else:
            real(op, body, key=key)
    win.engine.submit = submit
    return calls


def _real_primers(body):
    """Run the engine exactly as the worker would for a captured primers request."""
    import engine
    return engine.run_primers(body)


def _deliver(win, calls, key_op="primers"):
    op, body, key = next(c for c in calls if c[0] == key_op)
    res = _real_primers(body)
    win._on_done(key, res)
    return res


# =============================================================== pure helpers (coords.py)
def test_span_converts_only_the_start_and_marks_an_empty_interval():
    import coords
    assert coords.span(0, 276) == "0–276" and coords.span(0, 276, 1) == "1–276"
    assert coords.span(4870, 5146, 1) == "4871–5146"
    assert coords.span(5, 5, 0) == "5–5" and coords.span(5, 5, 1) == "5^6"


def test_pair_token_agrees_with_the_window_helper():
    import main, coords
    for text in ("D2·P3", "P1", "my pair #2", "···"):
        assert main.MainWindow._pair_token(text) == coords.pair_token(text)


def _rec(domains=None):
    return {"composition": {"length": 1000}, "seq": "A" * 1000, "domains": domains or [],
            "structural": [{"type": "LTR (terminal direct repeat)", "identity": 98.0, "five_prime": [0, 200],
                            "three_prime": [800, 1000], "element_span": [0, 1000]}]}


def _pair(i, lp, rp, f="ACGTACGTACGTACGTACGT", r="TGCATGCATGCATGCATGCA", **kw):
    return {"id": f"P{i}", "left_seq": f, "right_seq": r, "left_pos": lp, "right_pos": rp,
            "product_size": rp[1] - lp[0], "intended_product_span": [lp[0], rp[1]], **kw}


def test_location_names_the_segment_of_each_primer_and_the_domain():
    import coords
    rec = _rec([{"domain": "RT", "nt": [300, 700]}])
    loc = coords.pair_location(_pair(1, [310, 330], [670, 690]), rec, True)
    assert loc == "310–690 · internal (RT)"
    assert coords.pair_location(_pair(1, [150, 170], [670, 690]), rec, True) == "150–690 · 5′ LTR → internal (RT)"
    assert coords.pair_location(_pair(1, [190, 210], [670, 690]), rec, True) == "190–690 · 5′ LTR/internal → internal (RT)"
    assert coords.pair_location(_pair(1, [310, 330], [670, 690]), rec, False, template_label="LTR_5prime") == "310–690 on LTR_5prime"
    assert coords.pair_location(_pair(1, [310, 330], [670, 690]), rec, False) == "310–690 on a different record"
    assert coords.pair_location(_pair(1, [310, 330], [670, 690]), rec, True, 1).startswith("311–690")


def test_terminal_repeat_caution_needs_the_whole_3prime_core_inside_one_copy():
    import coords
    rec = _rec()
    inside = coords.pair_notes([_pair(1, [150, 175], [500, 520])], rec, True)[0]
    assert "F in 5′ LTR" in inside[0] and "98.0%" in inside[1] and "18 nt" in inside[1]
    straddle = coords.pair_notes([_pair(1, [190, 215], [500, 520])], rec, True)[0]
    assert straddle[0] == "—"
    assert coords.pair_notes([_pair(1, [150, 175], [500, 520])], rec, False)[0][0] == "—"      # not bound: no segments known
    tir = dict(rec, structural=[{"type": "TIR (terminal inverted repeat)", "five_prime": [0, 40], "three_prime": [960, 1000],
                                 "element_span": [0, 1000]}])
    note = coords.pair_notes([_pair(1, [100, 120], [980, 1000])], tir, True)[0]
    assert "single-primer" in note[0] and "R in 3′ TIR" in note[0]


def test_near_duplicates_and_shared_primers_are_flagged_without_double_tokens():
    import coords
    a = _pair(1, [300, 320], [600, 620])
    b = _pair(2, [301, 321], [601, 621], f="ACGTACGTACGTACGTAAAA")
    c = _pair(3, [400, 420], [700, 720])                                # shares F with pair 1
    d = _pair(4, [450, 470], [800, 820], f="GGGGCCCCGGGGCCCCGGGG", r="CCCCGGGGCCCCGGGGCCCC")
    notes = coords.pair_notes([a, b, c, d], _rec(), False)
    assert "≈P2" in notes[0][0] and "same F" not in notes[0][0].split("≈P2")[0]
    assert notes[1][0] == "≈P1"
    assert "same F as P1" in notes[2][0]
    assert notes[3][0] == "—"


def test_oligo_name_is_vendor_safe_and_prefix_optional():
    import coords
    assert coords.oligo_name("record_d6becc071991", "D1-P2", "F") == "record_d6becc071991_D1-P2_F"
    assert coords.oligo_name(None, "D1-P2", "R") == "D1-P2_R"
    assert coords.oligo_name("a b/c", "P1", "F") == "a_b_c_P1_F"


def test_order_sheet_merges_a_shared_primer_and_keeps_scale_and_purification_empty():
    import coords
    shared = "ACGTACGTACGTACGTACGT"
    p1 = _pair(1, [10, 30], [180, 200], f=shared, left_tm=60.1, right_tm=59.8, left_gc=50.0, right_gc=45.0, _name_prefix="record_x")
    p5 = _pair(5, [10, 30], [170, 190], f=shared.lower(), r="GGGGCCCCGGGGCCCCGGGG", left_tm=60.1, right_tm=61.0,
               left_gc=50.0, right_gc=70.0, _name_prefix="record_x")
    headers, rows, paste = coords.order_sheet([("D1·P1", p1), ("D1·P5", p5)])
    assert headers[:4] == ["Name", "Sequence (5'-3')", "Scale", "Purification"]
    assert len(rows) == 3 and len({r[0] for r in rows}) == 3
    first = rows[0]
    assert first[0] == "record_x_D1-P1_F" and first[1] == shared and first[2] == "" and first[3] == ""
    assert "D1·P1 F (190 bp)" in first[7] and "D1·P5 F (180 bp)" in first[7]
    assert all(line.count("\t") == 1 for line in paste.splitlines()) and len(paste.splitlines()) == 3


def test_order_sheet_tm_is_the_design_value_and_its_header_is_derived_from_the_seal():
    import coords
    conds = {"mv_conc": 50.0, "dv_conc": 1.5, "dntp_conc": 0.6, "dna_conc": 50.0}
    p = _pair(1, [10, 30], [180, 200], left_tm=60.1, right_tm=59.8, left_gc=50.0, right_gc=45.0,
              qc={"ok": True, "left": {"tm": 99.9}}, design_provenance={"parameters": {"tm_conditions": conds}})
    headers, rows, _ = coords.order_sheet([("P1", p)])
    assert rows[0][5] == 60.1 and "99.9" not in repr(rows)
    assert "1.5 mM Mg²⁺" in headers[5] and "50 nM primer" in headers[5]
    q = copy.deepcopy(p); q["design_provenance"]["parameters"]["tm_conditions"]["dv_conc"] = 2.5
    assert "2.5 mM Mg²⁺" in coords.order_sheet([("P1", q)])[0][5]
    assert "differ between pairs" in coords.order_sheet([("P1", p), ("P2", q)])[0][5]
    plain = {k: v for k, v in p.items() if k != "design_provenance"}
    assert coords.order_sheet([("P1", plain)])[0][5] == "Tm (°C; Primer3 design conditions)"


# =============================================================== coordinate display (B8)
def _coord_state(win):
    st = _struct_table(win), _table(win, "_ORFs"), _table(win, "_domains")
    return [(t.export_headers(), t.rows_data()) for t in st]


def test_coord_toggle_changes_the_display_and_never_a_file(win):
    import widgets
    rec = _accept(win, "M11240")
    s0 = _struct_table(win)
    before = _coord_state(win)
    ser0 = [widgets.serialize_table(h, r, ",") for h, r in before]
    ltr = next(e for e in rec["structural"] if e["type"].startswith("LTR"))
    fp, tp = ltr["five_prime"], ltr["three_prime"]
    assert s0.item(0, 1).text() == f"{fp[0]}–{fp[1]}  ·  {tp[0]}–{tp[1]}"
    win._set_coord_base(1)
    s1, o1, d1 = _struct_table(win), _table(win, "_ORFs"), _table(win, "_domains")
    assert win.coord_base == 1 and "1-BASED" in win.coordBtn.text().upper()
    assert s1.item(0, 1).text() == f"{fp[0] + 1}–{fp[1]}  ·  {tp[0] + 1}–{tp[1]}"
    orf = rec["orfs"][0]
    heads = [s1.horizontalHeaderItem(j).text() for j in range(s1.columnCount())]
    assert heads[1] == "Coords (1-based)"
    oheads = [o1.horizontalHeaderItem(j).text() for j in range(o1.columnCount())]
    row = next(r for r in range(o1.rowCount()) if o1.item(r, 0).text() == "1")
    assert o1.item(row, oheads.index("Start (1-based)")).text() == str(orf["start"] + 1)
    assert o1.item(row, oheads.index("End (inclusive)")).text() == str(orf["end"])
    dheads = [d1.horizontalHeaderItem(j).text() for j in range(d1.columnCount())]
    dom = rec["domains"][0]
    shown = [d1.item(r, dheads.index("nt (1-based)")).text() for r in range(d1.rowCount())]
    assert f"{dom['nt'][0] + 1}–{dom['nt'][1]}" in shown
    assert "aa in ORF" in dheads                                      # HMMER's own 1-based residue range is not converted
    after = _coord_state(win)
    assert after == before                                            # export headers and rows: byte for byte what they were
    assert [widgets.serialize_table(h, r, ",") for h, r in after] == ser0


def test_a_capped_orf_table_still_exports_canonical_numbers_in_its_sort_order(win):
    from PySide6.QtCore import Qt
    from teagle_core import domains as domains_mod
    rec = _accept(win, "M11240")
    n = domains_mod.MAX_ORFS_SCANNED + 5
    rec["orfs"] = [dict(rec["orfs"][0], start=1000 + 10 * i, end=1100 + 10 * i) for i in range(n)]
    win._set_coord_base(1)
    t = _table(win, "_ORFs")
    assert t._complete_rows is not None                               # the capped path
    heads = t.export_headers()
    assert heads[3] == "Start (0-based)" and heads[4] == "End (0-based, exclusive)"
    assert sorted(int(r[3]) for r in t.rows_data()) == [1000 + 10 * i for i in range(n)]
    t.sortByColumn(3, Qt.DescendingOrder)
    assert [int(r[3]) for r in t.rows_data()] == [1000 + 10 * i for i in reversed(range(n))]
    assert int(t.item(0, 3).text()) in {1001 + 10 * i for i in range(n)}      # the screen still shows 1-based


def test_copy_follows_the_screen_and_the_menu_names_the_convention(win):
    rec = _accept(win, "M11240")
    win._set_coord_base(1)
    t = _struct_table(win)
    t._copy_row(0)
    assert t.item(0, 1).text() in QApplication.clipboard().text()
    items = dict(win._feat_menu(0, 276, "+", "LTR", src_seq=rec["seq"]))
    label = next(l for l in items if l.startswith("Copy coords"))
    assert label == "Copy coords (1–276 +, 1-based)"
    items[label]()
    assert QApplication.clipboard().text() == "1-276 +"
    win._set_coord_base(0)
    assert any(l == "Copy coords (0–276 +)" for l, _f in win._feat_menu(0, 276, "+", "LTR", src_seq=rec["seq"]))


def test_coord_toggle_never_changes_a_design_coordinate(win):
    rec = _accept(win, "M11240")
    win._set_coord_base(1)
    calls = _capture(win)
    o = rec["orfs"][0]
    end = min(o["end"], o["start"] + 300)
    win._design_for_domain(o["start"], end, "ORF")
    body = calls[0][1]
    assert body["primer3_target"] == [o["start"], end - o["start"]]
    assert body["included_region"] == [0, len(rec["seq"])]


def test_coord_toggle_redraws_pcr_family_and_splice_with_canonical_exports(win):
    from widgets import DataTable
    rec = _accept(win, "M11240")
    win.state["splice_seq"] = rec["seq"]
    d = {"ok": True, "strand": "+", "counts": {"exons": 2, "introns": 1}, "canonical_introns": 1,
         "exons": [{"start": 0, "end": 100}, {"start": 400, "end": 600}],
         "introns": [{"start": 100, "end": 400, "donor": "GT", "acceptor": "AG", "canonical": True}]}
    win._render_splice(d)
    find = lambda body, suffix: next(t for t in body.parentWidget().findChildren(DataTable) if str(t.export_base).endswith(suffix))
    t = find(win.spliceBody, "introns")
    assert t.item(0, 1).text() == "100–400" and t.export_headers()[1] == "Intron span (0-based)"
    win._set_coord_base(1)
    t = find(win.spliceBody, "introns")
    assert t.item(0, 1).text() == "101–400" and t.rows_data()[0][1] == "100–400"
    assert t.export_headers()[1] == "Intron span (0-based)"
    hits = [{"class_family": "LTR/Copia", "family": "Copia_I", "q_start": 10, "q_end": 900, "strand": "+", "divergence": 2.0,
             "score": 500, "cons_start": 5, "cons_end": 880, "cons_length": 1000, "pct_del": 0.0, "pct_ins": 0.0,
             "cons_coverage_pct": 87.0, "n_fragments": 1}]
    win._render_family({"hits": hits, "species": "x", "repeatmasker_version": "4"})
    fam = next(t for t in win.wslBody.parentWidget().findChildren(DataTable) if "Dfam family" in t.export_headers())
    heads = [fam.horizontalHeaderItem(j).text() for j in range(fam.columnCount())]
    assert fam.item(0, heads.index("Coords (1-based)")).text() == "11–900"
    assert fam.rows_data()[0][fam.export_headers().index("Coords (0-based)")] == "10–900"
    assert "Consensus (1-based, inclusive)" in fam.export_headers()
    win._set_coord_base(0)                                            # and back, through the stored result
    fam = next(t for t in win.wslBody.parentWidget().findChildren(DataTable) if "Dfam family" in t.export_headers())
    assert fam.item(0, 3).text() == "10–900"


def test_coord_setting_persists_and_a_bad_value_reads_as_zero_based(monkeypatch):
    import main
    FakeSettings.store = {"quickstart_hidden": True}
    w = _make_win(monkeypatch)
    try:
        assert w.coord_base == 0 and "0-BASED" in w.coordBtn.text().upper()
        w._set_coord_base(1)
        assert FakeSettings.store["coord_base"] == 1
    finally:
        _dispose(w)
    w2 = _make_win(monkeypatch)
    try:
        assert w2.coord_base == 1 and "1-BASED" in w2.coordBtn.text().upper()
    finally:
        _dispose(w2)
    FakeSettings.store["coord_base"] = "x"
    assert main._saved_coord_base() == 0
    FakeSettings.store["coord_base"] = 7
    assert main._saved_coord_base() == 0


def test_the_viewer_export_ignores_the_display_base_but_its_label_does_not(win):
    from coords import CoordGenomePanel
    from widgets import GenomePanel
    _accept(win, "M11240")
    gv = win.card_struct.findChildren(GenomePanel)[0]
    assert isinstance(gv, CoordGenomePanel)
    svg0, pos0 = gv._cur_svg(920), gv.pos.text()
    win._set_coord_base(1)
    gv1 = win.card_struct.findChildren(GenomePanel)[0]
    assert gv1._cur_svg(920) == svg0
    assert gv1.pos.text() != pos0 and "figure text 0-based, half-open" in gv1.pos.text() and gv1.pos.text().startswith("1–")
    assert any("1-based, inclusive" in (r.get("tip") or "") for r in gv1.canvas.regions)


def test_the_1_based_header_tooltip_is_the_record_meaning_not_the_ispcr_one(win):
    _accept(win, "M11240")
    win._set_coord_base(1)
    tip = _struct_table(win).horizontalHeaderItem(1).toolTip()
    assert "1-based, inclusive [start, end]" in tip and "isPcr" not in tip


def test_the_coordinate_caption_follows_the_setting(win):
    _accept(win, "M11240")
    assert "Coords are 0-based and half-open" in _label_text(win.card_struct) and "Coords button" in _label_text(win.card_struct)
    win._set_coord_base(1)
    assert "Coordinates: 1-based, inclusive" in _label_text(win.card_struct)


# =============================================================== primer table (B9) and its conditions and legend (B12)
def _design_real(win, name="M11240"):
    rec = _accept(win, name)
    calls = _capture(win)
    win._design()
    res = _deliver(win, calls)
    return rec, res


def _primer_table(win):
    return next(t for t in _tables(win.primBody.parentWidget()).values() if t.export_base == "TEagle_primers")


def _heads(t):
    return [t.horizontalHeaderItem(j).text() for j in range(t.columnCount())]


def _by_object(root, name):
    return [l for l in root.findChildren(QLabel) if l.objectName() == name]


@functools.lru_cache(maxsize=None)
def _copia_design():
    import engine
    seq = _analysis("M11240")["records"][0]["seq"]
    return engine.run_primers({"sequence": ">r\n" + seq, "params": {}})


def test_primer_table_puts_the_verdict_and_ranking_first_and_keeps_every_cell_in_its_column(win):
    import main
    _design_real(win)
    t = _primer_table(win)
    assert _heads(t)[:5] == ["ID", "Struct", "Penalty", "Location", "Notes"]
    assert [h.replace(" (0-based)", "") for h in _heads(t)] == main.PRIMER_COLS
    cands = win.state["candidates"]
    assert cands
    for r in range(t.rowCount()):
        c = cands[t._orig(r)]
        assert t.item(r, 1).text() == (c["qc"].get("worst") if c["qc"].get("ok") else "n/a")
        assert t.item(r, 5).text() == c["left_seq"] and t.item(r, 6).text() == c["right_seq"]
        assert t.item(r, 2).text() == str(c["penalty"])
    heads = _heads(t)
    j = heads.index("Self-dim")
    assert t.rows_data()[0][j].lstrip("-").replace(".", "", 1).isdigit() or t.rows_data()[0][j] in ("none", "—")   # bare number


def test_location_and_notes_describe_each_pair_on_the_analysed_record(win):
    rec, res = _design_real(win)
    t = _primer_table(win)
    heads = _heads(t)
    locs = [t.item(r, heads.index("Location")).text() for r in range(t.rowCount())]
    assert all("–" in loc for loc in locs)
    assert any("internal" in loc or "LTR" in loc for loc in locs)
    c = win.state["candidates"][t._orig(0)]
    span = c["intended_product_span"]
    assert t.item(0, heads.index("Location")).text().startswith(f"{span[0]}–{span[1]}")
    assert t.rows_data()[0][heads.index("Location")].startswith(f"{span[0]}–{span[1]}")
    assert t.export_headers()[heads.index("Location")] == "Location (0-based, half-open)"
    assert "3′ core" in t.horizontalHeaderItem(heads.index("Notes")).toolTip() or "3′-terminal" in t.horizontalHeaderItem(heads.index("Notes")).toolTip()


def _cands_on(rec):
    F, R = "ACGTACGTACGTACGTACGT", "TGCATGCATGCATGCATGCA"
    mk = lambda i, lp, rp, f=F, r=R: {"id": f"P{i}", "left_seq": f, "right_seq": r, "left_pos": lp, "right_pos": rp,
                                 "product_size": rp[1] - lp[0], "intended_product_span": [lp[0], rp[1]], "left_tm": 60.0,
                                 "right_tm": 60.0, "left_gc": 50.0, "right_gc": 50.0, "penalty": 0.1}
    return [mk(1, [100, 120], [1000, 1020]), mk(2, [101, 121], [1001, 1021], "ACGTACGTACGTACGTAAAA"), mk(3, [2000, 2020], [3000, 3020], "GGGGCCCCGGGGCCCCGGGG", "CCCCGGGGCCCCGGGGCCCC")]


def test_notes_flag_the_terminal_repeat_primer_and_the_near_duplicate_on_screen(win):
    rec = _accept(win, "M11240")
    win._design_tmpl = rec["seq"]
    win._render_primers({"candidates": _cands_on(rec)})
    t = _primer_table(win)
    heads = _heads(t)
    notes = {t.item(r, 0).text(): t.item(r, heads.index("Notes")).text() for r in range(t.rowCount())}
    assert "F in 5′ LTR" in notes["P1"] and "≈P2" in notes["P1"]
    assert notes["P2"].endswith("≈P1") or "≈P1" in notes["P2"]
    assert notes["P3"] == "—"
    tip = t.item(0, heads.index("Notes")).toolTip()
    assert "100.0%" in tip and "solo LTRs" in tip


def test_a_pair_designed_on_another_template_gives_a_location_on_that_template(win):
    rec = _accept(win, "M11240")
    win._design_tmpl = rec["seq"][:3000]                              # candidates index a different sequence
    win._design_label = "LTR_5prime_upstream_0-250"
    win._render_primers({"candidates": _cands_on(rec)[:1]})
    t = _primer_table(win)
    loc = t.item(0, _heads(t).index("Location")).text()
    assert loc == "100–1020 on LTR_5prime_upstream_0-250" and "LTR →" not in loc


def test_the_explanations_fold_into_one_disclosure_and_the_key_lines_stay_visible(win):
    _design_real(win)
    root = win.primBody.parentWidget()
    toggles = [b for b in root.findChildren(QPushButton) if "How to read this table" in b.text()]
    assert len(toggles) == 1
    body = next(l for l in root.findChildren(QLabel) if "Location is the product span" in l.text())
    assert body.isHidden()
    legend = next(l for l in _by_object(root, "orient") if "caution" in l.text() and "warn" in l.text())
    cond = next(l for l in _by_object(root, "orient") if "Reaction conditions" in l.text())
    assert not legend.isHidden() and not cond.isHidden()
    toggles[0].click()
    assert not body.isHidden()


def test_primer_hint_names_row_actions_and_never_calls_right_click_the_only_way(win):
    _design_real(win)
    text = _label_text(win.primBody.parentWidget())
    assert "Row actions ▾" in text and "Right-click a pair" not in text


def test_caution_tiers_are_labelled_teagle_heuristic_and_derived_from_oligoqc(win, monkeypatch):
    from teagle_core import oligoqc
    _design_real(win)
    monkeypatch.setattr(oligoqc, "_HAIRPIN_CAUTION", -3.0)
    win._render_primers(win.state["lastPrimers"])
    root = win.primBody.parentWidget()
    legend = next(l for l in _by_object(root, "orient") if "caution" in l.text() and "warn" in l.text()).text()
    assert "TEagle heuristic" in legend and "hairpin ≤ −3" in legend and "dimers ≤ −5" in legend
    warn_part = legend.split("warn", 1)[1]
    assert "IDT rule of thumb" in warn_part and "IDT" not in legend.split("warn", 1)[0]
    notes = " ".join(_primer_table(win).export_notes)
    assert "hairpin <= -3" in notes and "TEagle heuristic" in notes
    monkeypatch.setattr(oligoqc, "THRESHOLDS", {"hairpin": {"caution": -4.0, "warn": -8.0}, "dimer": {"caution": -5.0, "warn": -9.0},
                                                 "end": {"caution": -6.0, "warn": -9.0}}, raising=False)
    win._render_primers(win.state["lastPrimers"])
    legend = next(l for l in _by_object(win.primBody.parentWidget(), "orient") if "caution" in l.text() and "warn" in l.text()).text()
    assert "hairpin ≤ −4" in legend and "hairpin ≤ −8" in legend       # per-axis warn values when they differ


def test_the_conditions_line_is_derived_from_the_result_and_exported(win):
    rec, res = _design_real(win)
    d = copy.deepcopy(res)
    win._render_primers(d)
    text = _label_text(win.primBody.parentWidget())
    assert "1.5 mM Mg²⁺" in text and "250 nM primer" in text
    d["design_tm_conditions"]["dv_conc"] = 2.5
    d["oligoqc_conditions"]["dna_conc"] = 100.0
    win._render_primers(d)
    text = _label_text(win.primBody.parentWidget())
    assert "2.5 mM Mg²⁺" in text and "100 nM primer" in text
    assert any("2.5 mM Mg²⁺" in n for n in _primer_table(win).export_notes)
    import widgets
    csv_text = widgets.serialize_table(_primer_table(win).export_headers(), _primer_table(win).rows_data(), ",",
                                       _primer_table(win).export_notes)
    assert "Reaction conditions:" in csv_text and "2.5 mM Mg" in csv_text


def test_a_design_without_usable_structure_qc_says_the_checks_are_unavailable(win):
    rec, res = _design_real(win)
    d = copy.deepcopy(res)
    for c in d["candidates"]:
        c["qc"] = {"ok": False, "error": "x"}
    win._render_primers(d)
    assert "ΔG checks unavailable for this design" in _label_text(win.primBody.parentWidget())
    t = _primer_table(win)
    assert t.item(0, 1).text() == "n/a"


def test_the_structure_dialog_states_the_tiers_from_oligoqc(win, monkeypatch):
    rec, res = _design_real(win)
    shown = []
    monkeypatch.setattr(QDialog, "exec", lambda self: shown.append(self) or 0)
    win._structure_detail(win.state["candidates"][0])
    text = _label_text(shown[0])
    assert "TEagle heuristic" in text and "IDT rule of thumb" in text
    assert text.count("IDT rule of thumb") == 1


def test_copy_pair_fasta_names_the_oligos_by_record_and_pair(win):
    import main
    rec, res = _design_real(win)
    t = _primer_table(win)
    items = dict(t._row_menu(0))
    items["Copy pair FASTA"]()
    c = win.state["candidates"][0]
    prefix = main._record_public_id(rec)
    text = QApplication.clipboard().text()
    assert text.startswith(f">{prefix}_{main.MainWindow._pair_token(main.MainWindow._pair_label(c))}_F\n{c['left_seq']}\n>")
    assert prefix.startswith("record_") and "_R\n" in text


def test_names_are_stamped_once_so_a_record_switch_does_not_rename_a_design(win):
    import main
    rec, res = _design_real(win)
    first = win.state["candidates"][0]["_name_prefix"]
    assert first == main._record_public_id(rec)
    win.state["last_rec"] = {"id": "record_ffffffffffff", "identitySource": "generated", "seq": "ACGT"}
    win._render_primers(win.state["lastPrimers"])
    assert win.state["candidates"][0]["_name_prefix"] == first
    assert win._oligo(win.state["candidates"][0], "F").startswith(first + "_")


def test_a_design_on_a_flank_carries_the_record_and_the_flank_in_its_oligo_names(win):
    import main
    rec = _accept(win, "M11240")
    calls = _capture(win)
    win._design_for_domain(0, 300, "LTR_5prime_upstream_0-300", seq=rec["seq"][:300])
    assert win._design_name_prefix == f"{main._record_public_id(rec)}_LTR_5prime_upstream_0-300"


# =============================================================== oligo order sheet (B15)
def _tsv(path):
    with open(path, encoding="utf-8-sig", newline="") as fh:
        return list(csv.reader(fh, delimiter="\t"))


def _patch_save(monkeypatch, path, fmt="tsv"):
    import widgets
    monkeypatch.setattr(widgets, "pick_table_format", lambda *a, **k: fmt)
    monkeypatch.setattr(widgets.QFileDialog, "getSaveFileName", staticmethod(lambda *a, **k: (str(path), "")))


def _button(root, text):
    return next(b for b in root.findChildren(QPushButton) if b.text().lower() == text.lower())


def test_panel_04_order_sheet_has_one_row_per_unique_oligo_and_no_note_lines(win, monkeypatch, tmp_path):
    rec, res = _design_real(win)
    cands = win.state["candidates"]
    uniq = {c[k].upper() for c in cands for k in ("left_seq", "right_seq")}
    out = tmp_path / "sheet.tsv"
    _patch_save(monkeypatch, out)
    _button(win.primBody.parentWidget(), "Export order sheet").click()
    rows = _tsv(out)
    assert rows[0][:4] == ["Name", "Sequence (5'-3')", "Scale", "Purification"]
    assert len(rows) - 1 == len(uniq)
    assert not any(r and r[0].startswith("#") for r in rows)
    assert all(r[2] == "" and r[3] == "" for r in rows[1:])
    assert len({r[0] for r in rows[1:]}) == len(rows) - 1
    assert "Mg²⁺" in rows[0][5]


def test_copy_order_rows_puts_name_tab_sequence_on_the_clipboard(win):
    import re
    _design_real(win)
    _button(win.primBody.parentWidget(), "Copy order rows").click()
    lines = QApplication.clipboard().text().splitlines()
    assert lines and all(re.fullmatch(r"[A-Za-z0-9_.-]+\t[ACGT]+", l) for l in lines)


def test_the_order_sheet_xlsx_carries_no_author_or_software_metadata(win, monkeypatch, tmp_path):
    zipfile = __import__("zipfile")
    rec, res = _design_real(win)
    out = tmp_path / "sheet.xlsx"
    _patch_save(monkeypatch, out, "xlsx")
    _button(win.primBody.parentWidget(), "Export order sheet").click()
    if not out.exists():
        pytest.skip("Excel export is unavailable on this machine")
    with zipfile.ZipFile(out) as z:
        meta = b"".join(z.read(n) for n in z.namelist() if n.startswith("docProps/")).decode("utf-8", "ignore").lower()
    import re
    assert "openpyxl" not in meta
    assert not re.search(r"<(dc:creator|cp:lastmodifiedby|application)>[^<]+<", meta)


def test_order_sheet_follows_the_on_screen_sort_and_is_disabled_after_invalidation(win, monkeypatch):
    from PySide6.QtCore import Qt
    rec, res = _design_real(win)
    t = _primer_table(win)
    t.sortByColumn(2, Qt.DescendingOrder)
    _button(win.primBody.parentWidget(), "Copy order rows").click()
    first_name = QApplication.clipboard().text().splitlines()[0].split("\t")[0]
    top = win.state["candidates"][t._orig(0)]
    assert first_name == win._oligo(top, "F")
    win._invalidate_analysis("x")
    assert not [b for b in win.primBody.parentWidget().findChildren(QPushButton) if "order" in b.text().lower() and b.isEnabled()]


def test_panel_05_order_sheet_follows_the_queue_and_is_disabled_when_it_is_empty(win, monkeypatch, tmp_path):
    rec, res = _design_real(win)
    cands = win.state["candidates"]
    assert not win.pcrOrderExport.isEnabled() and not win.pcrOrderCopy.isEnabled()
    win._add_pcr_pair(cands[2]); win._add_pcr_pair(cands[0])
    assert win.pcrOrderExport.isEnabled()
    win.pcrOrderCopy.click()
    names = [l.split("\t")[0] for l in QApplication.clipboard().text().splitlines()]
    assert names[0] == win._oligo(cands[2], "F")
    win._pcr_clear()
    assert not win.pcrOrderExport.isEnabled() and not win.pcrOrderCopy.isEnabled()
    win.pcrCustomF.setText("ACGTACGTACGTACGTAC"); win.pcrCustomR.setText("TGCATGCATGCATGCATG")
    win._add_custom_pair()
    out = tmp_path / "q.tsv"
    _patch_save(monkeypatch, out)
    win.pcrOrderExport.click()
    rows = _tsv(out)
    assert rows[1][0] == "custom1_F" and rows[1][5] == ""                        # no design Tm for a pair the user typed
    assert "None" not in repr(rows)


# =============================================================== assay goals (B14)
def _goal(win, key):
    win.pGoal.setCurrentIndex(win.pGoal.findData(key))


def _junction(win, key):
    win.pGoalJunction.setCurrentIndex(win.pGoalJunction.findData(key))


def _sent(calls):
    return [c for c in calls if c[0] == "primers"]


def test_goal_defaults_to_detect_with_the_strongest_domain_and_designs_inside_it(win):
    rec = _accept(win, "M11240")
    assert win.pGoal.currentData() == "detect" and win.designBtn.isEnabled()
    assert win.pGoalDomain.currentText().startswith("RT ·")                     # highest HMMER score first
    rt = next(d for d in rec["domains"] if d["domain"] == "RT")
    calls = _capture(win)
    win.designBtn.click()
    body = _sent(calls)[0][1]
    assert body["primer3_target"] is None and body["included_region"] == [rt["nt"][0], rt["nt"][1] - rt["nt"][0]]


@pytest.mark.parametrize("preset", ["standard", "qpcr"])
def test_a_detect_design_returns_pairs_whose_primers_and_product_stay_in_the_domain(win, preset):
    import engine
    rec = _accept(win, "M11240")
    win._apply_preset(preset)
    calls = _capture(win)
    win.designBtn.click()
    body = _sent(calls)[0][1]
    res = engine.run_primers(body)
    assert res["candidates"], "the routed design returned no pair"
    s, n = body["included_region"]
    for c in res["candidates"]:
        assert c["left_pos"][0] >= s and c["right_pos"][1] <= s + n
        assert c["left_pos"][0] >= 276 and c["right_pos"][1] <= 4870               # clear of both LTR copies
        assert c["primer3_target"] is None
    win._on_done(_sent(calls)[0][2], res)
    t = _primer_table(win)
    assert "within the RT domain" in _label_text(win.primBody.parentWidget()) or "inside the RT domain" in _label_text(win.primBody.parentWidget())
    assert any("(RT)" in t.item(r, _heads(t).index("Location")).text() for r in range(t.rowCount()))


def test_the_domain_is_clipped_to_the_interior_between_the_terminal_repeats(win):
    rec = _accept(win, "M11240")
    rec["domains"][0] = dict(rec["domains"][0], nt=[100, 900])            # reaches into the 5′ LTR (0–276)
    win._refresh_assay_goal()
    calls = _capture(win)
    win.pGoalDomain.setCurrentIndex(win.pGoalDomain.findData(0))
    win._design_for_goal()
    assert _sent(calls)[0][1]["included_region"] == [276, 900 - 276]


def test_a_domain_shorter_than_prod_min_is_not_submitted_and_says_why(win):
    _accept(win, "M11240")
    win.pfields["pMin"].setText("400"); win.pfields["pMax"].setText("500")
    calls = _capture(win)
    win.pGoalDomain.setCurrentIndex(win.pGoalDomain.findData(0))           # INT, ~306 bp
    win._design_for_goal()
    assert not calls
    assert "Prod min" in _label_text(win.primBody.parentWidget())
    assert win.designBtn.isEnabled()


def test_genotype_on_a_record_that_ends_at_the_element_says_what_to_fetch(win):
    _accept(win, "M11240")
    _goal(win, "genotype")
    assert not win.designBtn.isEnabled()
    assert "Fetch the insertion by coordinate" in win.goalGenoInfo.text() and "0–5146" in win.goalGenoInfo.text()
    calls = _capture(win)
    win._design_for_goal()                                                 # the button state is not the only guard
    assert not calls


def test_genotype_without_a_terminal_repeat_pair_says_the_edges_are_unknown(win):
    _accept(win, "M80343")
    _goal(win, "genotype")
    assert "No LTR or TIR pair" in win.goalGenoInfo.text() and not win.designBtn.isEnabled()


def test_a_junction_design_targets_the_two_bases_either_side_of_the_edge(win):
    import engine
    rec = _accept(win, "padded")
    assert win.pGoal.currentData() == "custom"                             # no domain: the whole-record custom range is the default
    _goal(win, "genotype")
    assert "Insertion 605–2261" in win.goalGenoInfo.text() and "upstream flank 605 bp" in win.goalGenoInfo.text()
    assert win.designBtn.isEnabled()
    n = len(rec["seq"])
    for which, s in (("5j", 605), ("3j", 2261)):
        _junction(win, which)
        calls = _capture(win)
        win.designBtn.click()
        body = _sent(calls)[0][1]
        assert body["primer3_target"] == [s - 1, 2] and body["included_region"] == [0, n]
        res = engine.run_primers(body)
        assert res["candidates"]
        assert all(c["left_pos"][1] <= s - 1 and c["right_pos"][0] >= s + 1 for c in res["candidates"])
        win._on_done(_sent(calls)[0][2], res)


def test_flank_to_flank_over_a_long_insertion_is_refused_with_the_arithmetic(win):
    _accept(win, "padded")
    _goal(win, "genotype"); _junction(win, "span")
    calls = _capture(win)
    win.designBtn.click()
    assert not calls
    text = _label_text(win.primBody.parentWidget())
    assert "Prod max" in text and "1,656" in text                          # the insertion's own length is stated


def test_custom_range_over_the_whole_record_is_the_untargeted_design_and_a_subrange_is_inside_only(win):
    rec = _accept(win, "padded")
    n = len(rec["seq"])
    _goal(win, "custom")
    assert (win.pGoalFrom.value(), win.pGoalTo.value()) == (1, n)
    calls = _capture(win)
    win.designBtn.click()
    whole = _sent(calls)[0][1]
    calls.clear(); win._design_inflight = False
    win._design()
    assert _sent(calls)[0][1] == whole and whole.get("primer3_target") is None and whole.get("included_region") is None
    calls.clear(); win._design_inflight = False
    win.pGoalFrom.setValue(701); win.pGoalTo.setValue(1500)
    win.designBtn.click()
    sub = _sent(calls)[0][1]
    assert sub["included_region"] == [700, 800] and sub["primer3_target"] is None


def test_the_insertion_button_in_card_02_opens_the_genotype_goal_and_is_absent_without_a_pair(win):
    _accept(win, "padded")
    btn = next(b for b in win.card_struct.findChildren(QPushButton) if b.text().lower().startswith("genotype this insertion"))
    btn.click()
    assert win.pGoal.currentData() == "genotype" and not win.card_primer.body.isHidden()
    _accept(win, "M80343")
    assert not [b for b in win.card_struct.findChildren(QPushButton) if b.text().lower().startswith("genotype this insertion")]


def test_the_goal_panel_follows_the_selected_record_and_is_disabled_without_a_result(win):
    _accept(win, "labelled", sequence=_labelled_text())
    first = [win.pGoalDomain.itemText(i) for i in range(win.pGoalDomain.count())]
    win._show_record(1)
    second = [win.pGoalDomain.itemText(i) for i in range(win.pGoalDomain.count())]
    assert first and second and first != second
    win._invalidate_analysis("x")
    assert not any(w.isEnabled() for w in (win.pGoal, win.pGoalDomain, win.pGoalJunction, win.pGoalFrom, win.pGoalTo,
                                           win.goalFlankBtn))
    assert not win.designBtn.isEnabled()


def test_a_coordinate_toggle_relabels_the_domain_choices_and_keeps_the_selection(win):
    _accept(win, "M11240")
    win.pGoalDomain.setCurrentIndex(1)
    before = win.pGoalDomain.currentData()
    text0 = win.pGoalDomain.currentText()
    win._set_coord_base(1)
    assert win.pGoalDomain.currentData() == before and win.pGoalDomain.currentText() != text0
    assert "nt 1863–2168" in win.pGoalDomain.currentText() or "nt 3189–3923" in win.pGoalDomain.currentText()


# =============================================================== your label, bundled examples, ID map (B16)
def _spin(pred, timeout=30):
    import time
    from PySide6.QtCore import QEventLoop, QTimer
    loop, t, t0 = QEventLoop(), QTimer(), time.time()
    t.setInterval(20)
    t.timeout.connect(lambda: (pred() or time.time() - t0 > timeout) and loop.quit())
    t.start(); loop.exec(); t.stop()
    return pred()


def _summary(win):
    return next(t for t in win.card_struct.findChildren(__import__("widgets").DataTable) if "_summary" in str(t.export_base))


def test_your_label_is_on_screen_only_and_reaches_no_file(win, tmp_path):
    import engine, main, widgets, zipfile
    from teagle_core import provenance
    result = _analysis("labelled")
    _accept(win, "labelled", sequence=_labelled_text())
    t = _summary(win)
    heads = _heads(t)
    assert heads[2] == "Your label" and heads[1] == "Record"
    assert [t.item(r, 2).text() for r in range(t.rowCount())] == ["patient_042", "C:\\Users\\person\\x.fa"]
    assert "Your label" not in t.export_headers() and len(t.export_headers()) == len(main.MainWindow._RECORD_HEADERS)
    secrets = ("patient_042", "person", "x.fa", "fly stock")
    blob = repr(t.rows_data()) + widgets.serialize_table(t.export_headers(), t.rows_data(), ",", t.export_notes)
    assert not any(s in blob for s in secrets)
    for ext in ("csv", "tsv", "xlsx"):
        p = tmp_path / f"s.{ext}"
        written = widgets.write_table(t.export_headers(), t.rows_data(), str(p), t.export_notes)
        raw = open(written, "rb").read()
        if written.endswith(".xlsx"):
            z = zipfile.ZipFile(written)
            raw = b"".join(z.read(n) for n in z.namelist())
        assert not any(s.encode() in raw for s in secrets)
    t._copy_row(0)
    assert not any(s in QApplication.clipboard().text() for s in secrets)
    for payload in (engine.analysis_result_json_bytes(win.state["analysis"]), main._full_records_fasta(win.state["records"]).encode(),
                    provenance.manifest_json_bytes(win.state["analysis"]["provenance"])):
        assert not any(s.encode() in payload for s in secrets)
    assert all("label" not in str(k).lower() for r in win.state["records"] for k in r)
    assert "Your full header (shown on this screen only): patient_042 fly stock" in t.item(0, 2).toolTip()


def test_the_title_of_the_detail_block_names_the_label_for_the_screen_only(win):
    _accept(win, "labelled", sequence=_labelled_text())
    titles = [l.text() for l in win.card_struct.findChildren(QLabel) if l.objectName() == "sech" and l.text().startswith("Record 1 of 2")]
    assert titles and "your label “patient_042”" in titles[0]


def test_a_label_needs_the_header_to_reproduce_the_record_identity(win):
    result = _analysis("labelled")
    wrong = "\n".join(f">other{i}\n{r['seq']}" for i, r in enumerate(result["records"]))
    _accept(win, "labelled", sequence=wrong)
    assert "Your label" not in _heads(_summary(win))
    _accept(win, "M11240", sequence=">\n".join(["", ""]))          # bare/garbled request text: never a wrong label
    assert "Your label" not in _heads(_summary(win))


def test_labels_are_cleared_on_invalidation(win):
    _accept(win, "labelled", sequence=_labelled_text())
    assert win._local_labels(win.state["records"]) == ["patient_042 fly stock", "C:\\Users\\person\\x.fa"]
    win._set_seq("ACGT" * 30)
    assert "record_headers" not in win.state and win._recTable is None


def test_the_capped_summary_path_drops_the_display_only_column_too(win):
    import engine, main
    text = "\n".join(f">sample{i}\n" + "ACGT" * 30 for i in range(main.MainWindow._MAX_RECORDS + 1))
    result = engine.analyze(text)
    _accept(win, None, result=result, sequence=text)
    t = _summary(win)
    assert t._complete_rows is not None and _heads(t)[2] == "Your label"
    rows = t.rows_data()
    assert len(rows) == main.MainWindow._MAX_RECORDS + 1 and len(rows[0]) == len(main.MainWindow._RECORD_HEADERS)
    assert len(t.export_headers()) == len(main.MainWindow._RECORD_HEADERS)
    assert "sample3" not in repr(rows)


def test_a_bundled_example_is_named_on_screen_but_provenance_stays_local_input(win):
    win._load_example("M11240")
    win._run_analysis()
    assert _spin(lambda: win.state.get("analysis_state") == "current")
    assert "bundled example" in win.rRecords.text() and "M11240.1" in win.rRecords.text()
    assert "identical" in win.rRecords.toolTip() and "fetch M11240" in win.rRecords.toolTip()
    assert win.state["analysis"]["provenance"]["input"]["inputType"] == "local_input"       # the deferred schema change
    seq = win.seq.toPlainText()
    win.seq.setPlainText(seq.replace("ACGT", "ACGA", 1) if "ACGT" in seq else seq + "A")      # an edit is not the bundled sequence
    win._run_analysis()
    assert _spin(lambda: win.state.get("analysis_state") == "current")
    assert "bundled example" not in win.rRecords.text() and "user-supplied" in win.rRecords.text()


def _patch_id_map(monkeypatch, accept, path=None):
    import widgets
    calls = []
    monkeypatch.setattr(QMessageBox, "exec", lambda self: 0)
    monkeypatch.setattr(QMessageBox, "clickedButton",
                        lambda self: next(b for b in self.buttons() if b.text().startswith("Export ID map")) if accept else None)
    monkeypatch.setattr(widgets.QFileDialog, "getSaveFileName",
                        staticmethod(lambda *a, **k: calls.append(a) or ((str(path), "") if path else ("", ""))))
    return calls


_ID_MAP = "Export ID map (local)…"


def _id_map_action(win, row=0):
    return dict(_summary(win)._row_menu(row)).get(_ID_MAP)


def _id_map_buttons(win):
    return [b for b in win.card_struct.findChildren(QPushButton) if "id map" in b.text().lower()]


def test_the_id_map_is_offered_only_from_the_record_row_menu_after_confirmation(win, monkeypatch, tmp_path):
    import hashlib
    _accept(win, "labelled", sequence=_labelled_text())
    assert not _id_map_buttons(win)
    assert _summary(win)._row_menu(0)[0][0] == "Show this record below"
    action = _id_map_action(win)
    assert action is not None
    before = copy.deepcopy(win.state.get("prov_manifest"))
    calls = _patch_id_map(monkeypatch, False)
    action()
    assert not calls and not list(tmp_path.iterdir())
    out = tmp_path / "map.csv"
    calls = _patch_id_map(monkeypatch, True, out)
    action()
    with open(out, encoding="utf-8-sig", newline="") as fh:
        rows = [r for r in csv.reader(fh) if r and not r[0].startswith("#")]
    assert rows[0][:3] == ["Record ID", "Record ID SHA-256", "FASTA header (local)"] and len(rows) == 3
    for rid, digest, header, *_ in rows[1:]:
        assert hashlib.sha256(header.encode("utf-8")).hexdigest() == digest and digest.startswith(rid.split("_")[1])
    assert {r[2] for r in rows[1:]} == {"patient_042 fly stock", "C:\\Users\\person\\x.fa"}
    assert win.state.get("prov_manifest") == before


def test_the_id_map_confirmation_names_the_original_sequence_names_and_defaults_to_cancel(win, monkeypatch):
    _accept(win, "labelled", sequence=_labelled_text())
    seen = {}

    def exec_(self):
        seen.update(text=self.text() + " " + self.informativeText(), default=self.defaultButton().text(),
                    buttons=[b.text() for b in self.buttons()])
        return 0
    monkeypatch.setattr(QMessageBox, "exec", exec_)
    _id_map_action(win)()
    assert "original sequence names" in seen["text"] and "sensitive" in seen["text"] and "not share it" in seen["text"]
    assert seen["default"] == "Cancel" and seen["buttons"][0].startswith("Export ID map")


def test_the_id_map_is_absent_for_a_bare_paste_and_not_reachable_once_the_result_is_gone(win):
    from PySide6.QtCore import Qt
    _accept(win, "labelled", sequence=_labelled_text())
    table = _summary(win)
    win._invalidate_analysis("x")
    try:
        assert not table.isEnabled() and table.contextMenuPolicy() == Qt.NoContextMenu
    except RuntimeError:
        pass                                                   # the table was removed with its result
    import engine
    bare = "ACGT" * 60
    _accept(win, None, result=engine.analyze(bare), sequence=bare)
    assert not _id_map_buttons(win) and _id_map_action(win) is None


def test_a_candidate_with_no_qc_field_still_shows_the_conditions_the_result_records(win):
    rec = _accept(win, "M11240")
    win._design_tmpl = rec["seq"]
    d = {"candidates": _cands_on(rec)[:1], "design_tm_conditions": {"mv_conc": 50.0, "dv_conc": 1.5, "dntp_conc": 0.6, "dna_conc": 50.0},
         "oligoqc_conditions": {"mv_conc": 50.0, "dv_conc": 0.0, "dntp_conc": 0.0, "dna_conc": 250.0, "temp_c": 37.0}}
    win._render_primers(d)
    text = _label_text(win.primBody.parentWidget())
    assert "37 °C" in text and "ΔG checks unavailable" not in text


def test_the_detect_note_claims_avoidance_of_terminal_repeats_only_when_there_are_some(win):
    _accept(win, "M11240")
    assert "avoid both copies" in win.goalDetectNote.text()
    _accept(win, "M80343")
    note = win.goalDetectNote.text()
    assert "terminal repeats" not in note and "inside this domain" in note
