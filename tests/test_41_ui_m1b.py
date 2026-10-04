"""4.1 UI results lane: viewer and table menus, record-card honesty, structural-evidence table, exports,
family and splice panels. Drives the real MainWindow offscreen; no engine call leaves the process."""
import os, sys, functools

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
import pytest

pytest.importorskip("PySide6")
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "app", "native"))
sys.path.insert(0, os.path.join(_ROOT, "app", "backend"))
sys.path.insert(0, os.path.join(_ROOT, "tests"))

from PySide6.QtWidgets import QApplication, QLabel, QPushButton

_app = QApplication.instance() or QApplication([])


@pytest.fixture
def win():
    import gc
    import fonts
    import main
    from PySide6.QtCore import QEvent
    fonts.load_fonts()
    w = main.MainWindow()
    _app.processEvents()
    yield w
    w.engine.wait(5000)
    w.close(); w.deleteLater()
    QApplication.sendPostedEvents(None, QEvent.DeferredDelete)
    del w
    gc.collect()


@functools.lru_cache(maxsize=None)
def _analysis(name):
    import engine, helpers
    from teagle_core import examples, sequtil
    if name in ("M11240", "M80343", "X01005", "X05424"):
        return engine.analyze(examples.load(name))
    if name == "M11240_rc":
        s = sequtil.parse_fasta(examples.load("M11240"))[0][1]
        return engine.analyze(">rc\n" + sequtil.reverse_complement(s))
    if name == "M11240_3kb":
        s = sequtil.parse_fasta(examples.load("M11240"))[0][1]
        return engine.analyze(">cut\n" + s[:3000])
    if name == "tdr":
        return engine.analyze(">tdr\n" + helpers.make_ltr_element(ltr_len=300))
    if name == "mixed":
        return engine.analyze(">bad\nACGT-ACGTAC\n>ok\n" + helpers.make_tir_element())
    raise KeyError(name)


def _accept(win, name):
    import copy
    result = copy.deepcopy(_analysis(name))
    sequence = "\n".join(f">{r['id']}\n{r['seq']}" for r in result.get("records", []) if r.get("seq"))
    win.seq.setPlainText(sequence)
    win._analyze_gen += 1
    generation = win._analyze_gen
    win._active_analyze_gen = generation
    win._analysis_requests[generation] = {"sequence": sequence, "source": {}}
    win._analyze_inflight = True
    win._on_analyze(result, generation)
    return win.state["records"][win.state.get("record_index", 0)]


def _texts(widget):
    return " ".join(l.text() for l in widget.findChildren(QLabel))


def _card(win):
    return _texts(win.card_struct)


def _clip():
    return QApplication.clipboard().text()


def _viewer(win, index=0):
    from widgets import GenomePanel
    return win.card_struct.findChildren(GenomePanel)[index]


def _region(gv, kind, **match):
    for r in gv.canvas.regions:
        if r.get("kind") == kind and all(r.get(k) == v for k, v in match.items()):
            return r
    raise AssertionError(f"no {kind} region in {[(r.get('kind'), r.get('start')) for r in gv.canvas.regions]}")


def _tables(root, base_suffix=None):
    from widgets import DataTable
    ts = root.findChildren(DataTable)
    return ts if base_suffix is None else [t for t in ts if str(t.export_base).endswith(base_suffix)]


# ---------------------------------------------------------------- viewer menus (ui_results-1, -4, -13)
def test_minus_strand_domain_glyph_menu_reports_and_copies_its_own_strand(win):
    rec = _accept(win, "M11240_rc")
    dom = next(d for d in rec["domains"] if d["domain"] == "RT")
    assert dom["strand"] == "-"
    gv = _viewer(win)
    items = dict(gv.on_feature_menu(_region(gv, "domain", start=dom["nt"][0])))
    assert f"Copy coords ({dom['nt'][0]}–{dom['nt'][1]} -)" in items
    assert "Copy protein" in items
    items["Copy DNA"]()
    assert _clip() == dom["dna"]


def test_minus_strand_gene_model_region_copies_the_reverse_complement(win):
    from teagle_core import sequtil
    seq = "ACGTTGCAAC" * 12
    win.state["seq"] = win.state["analyzed_seq"] = seq
    model = {"length": len(seq), "tracks": [{"name": "exons", "height": 20, "features": [
        {"start": 10, "end": 70, "color": "#009E73", "label": "exon", "kind": "exon", "strand": "-"}]}]}
    from widgets import GenomePanel
    from figures import svg_genome
    gv = GenomePanel(svg_genome, "t"); gv.set_model(model)
    gv.set_feature_menu(lambda r: win._region_menu(r, seq))
    r = _region(gv, "exon")
    items = dict(gv.on_feature_menu(r))
    assert "Copy coords (10–70 -)" in items
    items["Copy DNA"]()
    assert _clip() == sequtil.reverse_complement(seq[10:70])


def test_repeat_glyph_in_the_viewer_gets_copy_actions_only(win):
    _accept(win, "M11240")
    gv = _viewer(win)
    labels = [l for l, _ in gv.on_feature_menu(_region(gv, "LTR"))]
    assert "Copy DNA" in labels
    assert not any("Design primer" in l or "splice" in l.lower() or "sub-region" in l for l in labels)


def test_advisory_repeat_row_menu_gets_copy_actions_only(win):
    e = {"type": "Sub-threshold terminal direct repeat (advisory)", "five_prime": [0, 100],
         "three_prime": [300, 400], "element_span": [0, 400], "ltr_len": 100, "identity": 70.0, "advisory": True}
    win.state["analyzed_seq"] = "ACGT" * 100
    labels = [l for l, _ in win._struct_menu(e)]
    assert labels and not any("Design primer" in l or "splice" in l.lower() for l in labels)


def test_splice_viewer_menu_uses_the_spliced_sequence_without_a_prior_analysis(win):
    seq = "ACGTACGTAC" + "GGATTCCATG" * 14 + "TTTT"
    win.seq.setPlainText(">insulin like\n" + seq)
    win.state["seq"] = ">insulin like\n" + seq
    win.state["splice_seq"] = seq
    win.state.pop("analyzed_seq", None)
    win._render_splice({"counts": {"exons": 2, "introns": 1}, "canonical_introns": 1, "strand": "+",
                        "exons": [{"start": 0, "end": 60}, {"start": 90, "end": 150}],
                        "introns": [{"start": 60, "end": 90, "donor": "GT", "acceptor": "AG", "canonical": True}]})
    from widgets import GenomePanel
    gv = win.spliceBody.parentWidget().findChildren(GenomePanel)[0]
    items = dict(gv.on_feature_menu(_region(gv, "exon", start=0)))
    items["Copy DNA"]()
    assert _clip() == seq[0:60]
    items["Flanking sequence (upstream / downstream)…"]()
    dlg = win._flank_dlg
    side = dlg.findChildren(__import__("PySide6.QtWidgets", fromlist=["QComboBox"]).QComboBox)[0]
    side.setCurrentIndex(1)                                                   # downstream
    next(b for b in dlg.findChildren(QPushButton) if b.text().upper() == "COPY FASTA").click()
    assert "insulin" not in _clip() and seq[60:160] in _clip().replace("\n", "")


def test_the_three_flank_entries_of_a_repeat_row_name_their_span(win):
    e = {"type": "LTR (terminal direct repeat)", "five_prime": [0, 100], "three_prime": [300, 400],
         "element_span": [0, 400], "ltr_len": 100, "identity": 95.0}
    win.state["analyzed_seq"] = "ACGT" * 100
    flanks = [l for l, _ in win._struct_menu(e) if "lank" in l]
    assert len(flanks) == len(set(flanks)) >= 1
    assert all(("arm" in l or "element" in l) for l in flanks)


# ---------------------------------------------------------------- flank dialog (ui_results-6)
def _flank_dialog(win, seq, start, end):
    win._flank_picker(start, end, "+", "LTR_element", src_seq=seq)
    dlg = win._flank_dlg
    from PySide6.QtWidgets import QComboBox, QSpinBox
    return dlg, dlg.findChildren(QComboBox)[0], dlg.findChildren(QSpinBox)[0]


def _button(dlg, text):
    return next(b for b in dlg.findChildren(QPushButton) if b.text().upper().startswith(text.upper()))


def test_both_flanks_design_spans_the_element_with_primers_in_each_flank(win, monkeypatch):
    seq = "ACGTTGCA" * 200                                     # 1600 bp
    calls = []
    monkeypatch.setattr(win, "_design_for_domain", lambda s, e, lab, seq=None: calls.append((s, e, lab, seq)))
    dlg, side, spin = _flank_dialog(win, seq, 600, 660)
    side.setCurrentIndex(side.findData("both")); spin.setValue(200)
    _button(dlg, "Design").click()
    s, e, lab, sub = calls[-1]
    assert sub == seq[400:860] and (s, e) == (200, 260)         # target = the element, flanks on both sides
    assert "across" in lab or "element" in lab


def test_both_flanks_design_is_refused_when_the_element_exceeds_the_product_range(win, monkeypatch):
    seq = "ACGTTGCA" * 400                                     # 3200 bp
    calls = []
    monkeypatch.setattr(win, "_design_for_domain", lambda *a, **k: calls.append(a))
    dlg, side, spin = _flank_dialog(win, seq, 500, 2500)
    side.setCurrentIndex(side.findData("both"))
    assert not _button(dlg, "Design").isEnabled()
    assert "longer than" in _texts(dlg) and not calls


def test_single_flank_design_says_it_stays_inside_the_flank(win):
    seq = "ACGTTGCA" * 200
    dlg, side, spin = _flank_dialog(win, seq, 600, 660)
    assert "inside" in _button(dlg, "Design").toolTip().lower()


# ---------------------------------------------------------------- invalid and not-run records (ui_results-5)
def test_invalid_record_is_reported_as_not_analysed_never_as_no_domain_found(win):
    _accept(win, "mixed")
    win._show_record(0)
    text = _card(win)
    assert "No panel domain was detected" not in text and "None were detected" not in text
    assert "not analysed" in text.lower()
    assert not win.designBtn.isEnabled()
    rows = {r[1]: r for r in win._record_summary_rows(win.state["records"])}
    row = dict(zip(win._RECORD_HEADERS, win._record_summary_rows(win.state["records"])[0]))
    assert row["Domain evidence"] == "not_run" and row["Domain architecture"] != "unassigned"
    assert not any(b.text().upper().startswith("EXPORT ANNOTATION") for b in win.card_struct.findChildren(QPushButton))


def test_summary_domain_evidence_follows_the_processing_status(win):
    rec = _accept(win, "M80343")
    row = dict(zip(win._RECORD_HEADERS, win._record_summary_rows([rec])[0]))
    assert row["Domain evidence"].startswith("tested_with_hits") and "8 ORFs not searched" in row["Domain evidence"]
    rec = _accept(win, "M11240")
    row = dict(zip(win._RECORD_HEADERS, win._record_summary_rows([rec])[0]))
    assert row["Domain evidence"] == "tested_with_hits"


# ---------------------------------------------------------------- headline and chip (ui_results-19)
def test_lone_reverse_transcriptase_chip_names_class_and_order(win):
    from teagle_core import classify
    dom = {"domain": "RT", "label": "reverse transcriptase", "pfam": "PF00078", "strand": "+", "nt": [100, 400],
           "aa": [1, 100], "orf": 0, "score": 50.0, "evalue": 1e-20, "confidence": "high", "protein": "M",
           "dna": "ATG"}
    cl = classify.classify([], [dom])
    rec = {"id": "lone_rt", "seq": "ACGT" * 150, "valid": True, "classification": cl,
           "completeness": cl["completeness"], "structural": [], "domains": [dom], "orfs": [],
           "composition": {"length": 600, "gc": 50.0, "n": 0.0, "ambiguous": 0.0},
           "processing": {"domains": {"status": "complete"}}}
    win.state["records"] = [rec]
    win._render_struct_card(rec, {})
    chips = [l.text() for l in win.card_struct.findChildren(QLabel) if l.objectName() == "cf"]
    assert "CLASS & ORDER WITHHELD" in chips and "ORDER WITHHELD" not in chips


def test_unavailable_domain_scan_headline_is_not_an_absence_call(win):
    from teagle_core import classify
    cl = classify.classify([], [], domains_ok=False)
    rec = {"id": "x", "seq": "ACGT" * 100, "valid": True, "classification": cl, "completeness": cl.get("completeness"),
           "structural": [], "domains": [], "orfs": [], "composition": {"length": 400, "gc": 50.0, "n": 0.0},
           "processing": {"domains": {"status": "failed"}}}
    win.state["records"] = [rec]
    win._render_struct_card(rec, {})
    text = _card(win)
    assert "Not assessed" in text and "No TE signature detected" not in text and "no clear TE signature" not in text


# ---------------------------------------------------------------- truncation clause (persona-U1)
def test_truncation_clause_follows_the_coding_signal_not_the_note_text(win):
    rec = _accept(win, "M80343")
    rec["notes"] = ["sequence appears 3′-truncated — some ORFs are open-ended"]
    win._render_struct_card(rec, {})
    assert "3′-truncated" not in _texts(win.card_struct).split("Structural completeness:")[1].split("Domain architecture")[0]
    rec["classification"]["completeness"]["coding"]["truncated_frame"] = True
    win._render_struct_card(rec, {})
    line = _texts(win.card_struct).split("Structural completeness:")[1].split("Domain architecture")[0]
    assert "runs off the end of the supplied sequence" in line


# ---------------------------------------------------------------- summary columns (persona-U11)
def test_summary_domain_architecture_matches_the_card_and_wicker_order_is_separate(win):
    rec = _accept(win, "M80343")
    row = dict(zip(win._RECORD_HEADERS, win._record_summary_rows([rec])[0]))
    assert row["Domain architecture"] == "ORF1–EN–RT" and row["Wicker order"] == "LINE"
    rec = _accept(win, "tdr")
    row = dict(zip(win._RECORD_HEADERS, win._record_summary_rows([rec])[0]))
    assert row["Domain architecture"] == "none detected" and row["Wicker order"] == "withheld"
    assert row["Structural"].startswith("TDR")


# ---------------------------------------------------------------- structural table (U17, U28, U29, exports-5)
def test_order_withheld_candidate_is_labelled_tdr_in_table_and_text(win):
    rec = _accept(win, "tdr")
    t = _tables(win.card_struct, "_structural")[0]
    feats = [t.item(i, 0).text() for i in range(t.rowCount())]
    assert "TDR (terminal direct repeat)" in feats and not any(f.startswith("LTR") for f in feats)
    text = _card(win)
    assert "Structural: TDR" in text and "Structural: LTR" not in text


def test_tsd_row_shows_both_copies_and_its_chance_level(win):
    rec = _accept(win, "tdr")
    tsd = next(e for e in rec["structural"] if e["type"].startswith("TSD"))
    t = _tables(win.card_struct, "_structural")[0]
    row = next(i for i in range(t.rowCount()) if t.item(i, 0).text().startswith("TSD"))
    up, dn = tsd["upstream"], tsd["downstream"]
    assert t.item(row, 1).text().replace("  ", " ") == f"{up[0]}–{up[1]} · {dn[0]}–{dn[1]}"
    assert "chance" in t.item(row, 3).text().lower()


def test_lower_bound_repeat_length_is_hedged_in_table_tooltip_and_export_note(win):
    rec = _accept(win, "M11240")
    e = rec["structural"][0]
    assert e.get("length_is_lower_bound")
    row, tips = win._struct_row(e), win._struct_tips(e)
    assert "lower bound" in row[3] and row[2] == e["ltr_len"]
    assert "record ends" in (tips[2] or "")
    t = _tables(win.card_struct, "_structural")[0]
    assert any("lower bound" in n for n in t.export_notes)


def test_missing_tsd_and_pbs_are_reported_as_searched_not_found_with_their_scope(win):
    _accept(win, "X01005")
    assert "TSD — not assessable" in _card(win)
    _accept(win, "M11240")
    text = _card(win)
    assert "PBS" in text and "tRNA-Lys3" in text and "another tRNA" in text


def test_undetermined_cis_orientation_is_explained_on_the_card(win):
    _accept(win, "tdr")
    assert "orientation of the element is undetermined" in _card(win)


def test_minus_strand_cis_row_is_marked_and_copied_on_its_own_strand(win):
    e = {"type": "PPT (polypurine tract)", "pos": [10, 19], "length": 9, "motif": "GAGGGGGCG", "strand": "-"}
    assert "(−)" in win._struct_row(e)[1]
    win.state["analyzed_seq"] = "ACGTTGCAAC" * 10
    items = dict(win._struct_menu(e))
    assert any(l.startswith("Copy coords") and l.endswith("-)") for l in items)


# ---------------------------------------------------------------- self-similarity guides (structural-14)
def test_selfsim_guides_skip_rejected_pairs_and_colour_by_type(win, monkeypatch):
    import widgets
    seen = {}
    monkeypatch.setattr(widgets.DotPanel, "set_matrix", lambda self, m, guides=None, **k: seen.update(guides=guides))
    monkeypatch.setattr(widgets.DotPanel, "_fit", lambda self: None)
    monkeypatch.setattr(type(win), "_uppercase_buttons", lambda self, root=None: None)
    rec = {"id": "r", "seq": "ACGT" * 300, "structural": [
        {"type": "LTR (terminal direct repeat)", "five_prime": [0, 100], "three_prime": [1000, 1100]},
        {"type": "TIR (terminal inverted repeat)", "five_prime": [0, 30], "three_prime": [1170, 1200]},
        {"type": "Sub-threshold terminal direct repeat (advisory)", "five_prime": [200, 300],
         "three_prime": [600, 700], "advisory": True}]}
    win._build_selfsim(rec)
    g = seen["guides"]
    assert [(x["start"], x["end"]) for x in g] == [(0, 100), (1000, 1100), (0, 30), (1170, 1200)]
    assert g[0]["color"] != g[2]["color"]


# ---------------------------------------------------------------- domain and ORF tables (ui_results-10, -16, B22)
def test_domain_table_has_strand_and_orf_and_names_the_residue_convention(win):
    rec = _accept(win, "M11240_rc")
    t = _tables(win.card_struct, "_domains")[0]
    heads = [t.horizontalHeaderItem(c).text() for c in range(t.columnCount())]
    assert "Strand" in heads and "ORF" in heads
    assert "aa" not in heads and any(h.startswith("aa in ORF") for h in heads)
    rows = t.rows_data()
    assert {r[heads.index("Strand")] for r in rows} == {"-"}
    exp = t.export_headers()
    assert any("1-based" in h for h in exp)
    assert "1-based" in t.horizontalHeaderItem([h.startswith("aa in ORF") for h in heads].index(True)).toolTip()


def test_table_exports_carry_the_record_in_name_and_notes(win):
    from teagle_core import sequtil
    rec = _accept(win, "M11240")
    stem = main_stem(rec)
    for suffix in ("_structural", "_ORFs", "_domains"):
        ts = _tables(win.card_struct, suffix)
        assert ts, suffix
        t = ts[0]
        assert stem in t.export_base
        notes = " ".join(t.export_notes)
        assert rec["idSha256"][:12] in notes or rec["id"] in notes
        assert str(rec["composition"]["length"]) in notes.replace(",", "")


def main_stem(rec):
    import main
    return main._record_export_stem(rec)


def test_orf_table_lists_searched_status_and_shows_all_on_request(win):
    rec = _accept(win, "M80343")                               # 20 ORFs, 12 searched
    t = _tables(win.card_struct, "_ORFs")[0]
    heads = [t.horizontalHeaderItem(c).text() for c in range(t.columnCount())]
    assert "Domains" in heads
    assert len(t.rows_data()) == len(rec["orfs"])                  # export keeps every ORF
    assert t.rowCount() == 12
    btn = next(b for b in win.card_struct.findChildren(QPushButton) if "Show all" in b.text())
    btn.click()
    assert t.rowCount() == len(rec["orfs"])


# ---------------------------------------------------------------- methods text (ui_results-17)
def test_methods_text_states_the_real_scan_scope():
    import main
    from teagle_core import domains
    html = main.MainWindow._methods_html(type("S", (), {
        "_domain_panel_html": staticmethod(main.MainWindow._domain_panel_html),
        "_src_html": lambda self, k: ""})())
    assert f"{domains.MAX_ORFS_SCANNED} longest" in html
    assert "i-Evalue &lt; 1e-3" in html or "i-Evalue < 1e-3" in html
    assert "tRNA-Lys3" in html and "only" in html.split("tRNA-Lys3")[1][:120]
    assert "partitions" in html


# ---------------------------------------------------------------- WSL status (ui_family_splice_prov-7, -10, -11)
_REAL_NO_WSL = {"wsl2": False, "distro": None, "error": "WSL not installed / no distro", "ready": False,
                "repeatmasker": None, "dfam": False, "dfam_uncurated": False, "minimap2": None}


def test_missing_wsl_shows_install_guidance_and_updates_the_splice_card(win):
    win._on_wsl_status(dict(_REAL_NO_WSL))
    assert "WSL2 not installed" in win.wslStatus.text() and "error" not in win.wslStatus.text().lower()
    assert "checking" not in win.spliceStatus.text() and "not installed" in win.spliceStatus.text()


def test_unreadable_wsl_status_still_reads_as_an_error_on_both_cards(win):
    win._on_wsl_status({"wsl2": False, "distro": "Ubuntu-24.04", "reason": "timeout",
                        "error": "the distro did not answer in 30 s"})
    assert "status error" in win.wslStatus.text() and "checking" not in win.spliceStatus.text()


def test_reprobe_that_finds_the_backend_gone_disarms_the_run_buttons(win):
    ready = {"wsl2": True, "ready": True, "distro": "Ubuntu-24.04", "repeatmasker": "4.2.4", "minimap2": "2.28",
             "dfam_uncurated": True, "dfam_library": {"partitions": ["a: uncurated", "b: uncurated"]}}
    win._on_wsl_status(ready)
    win.wslLibrary.setCurrentIndex(1)
    assert win.annotateBtn.isEnabled() and win.spliceBtn.isEnabled() and win._uncurated_available
    win._on_wsl_status({"wsl2": True, "ready": False, "distro": "Ubuntu-24.04", "repeatmasker": None,
                        "dfam": False, "minimap2": None})
    assert not win.annotateBtn.isEnabled() and not win.spliceBtn.isEnabled()
    win._on_wsl_status(ready)
    win._on_wsl_status(dict(_REAL_NO_WSL))
    assert not win.annotateBtn.isEnabled() and not win.spliceBtn.isEnabled()
    assert win._uncurated_available is False and win.wslLibraryRow.isHidden()


def test_partial_uncurated_install_does_not_point_at_the_hidden_library_control(win):
    win._on_wsl_status({"wsl2": True, "ready": True, "distro": "Ubuntu-24.04", "repeatmasker": "4.2.4",
                        "dfam_uncurated": False,
                        "dfam_library": {"partitions": ["dfam40.0.h5", "dfam40.uncurated.consensus.0.h5"]}})
    win._render_family({"ok": True, "hits": [], "species": "Arabidopsis thaliana", "species_defaulted": False,
                        "include_uncurated": False,
                        "dfam_library": {"partitions": ["dfam40.0.h5", "dfam40.uncurated.consensus.0.h5"]}})
    text = _texts(win.wslBody.parentWidget())
    assert "Set Library" not in text and "Set <b>Library" not in text and "missing" in text.lower()


# ---------------------------------------------------------------- family panel (ui_family_splice_prov-6, -12, R24)
def _hit(cf, fam, s, e, score, strand="+"):
    return {"class_family": cf, "family": fam, "q_start": s, "q_end": e, "strand": strand, "divergence": 2.0,
            "score": score, "pct_del": 0.0, "pct_ins": 0.0, "cons_start": 1, "cons_end": e - s,
            "cons_length": e - s, "cons_coverage_pct": 99.0, "n_fragments": 1}


def test_family_header_names_the_dominant_te_and_rna_genes_are_not_te_families(win):
    hits = [_hit("SINE/Alu", "AluSx", 0, 300, 2000), _hit("LTR/ERVK", "HERVK-int", 400, 8000, 52000),
            _hit("snRNA", "U2", 8100, 8300, 900), _hit("Satellite/centr", "ALR/Alpha", 8400, 8800, 700)]
    win._render_family({"ok": True, "hits": hits, "species": "Homo sapiens", "repeatmasker_version": "4.2",
                        "include_uncurated": False})
    root = win.wslBody.parentWidget()
    head = next(l.text() for l in root.findChildren(QLabel) if l.text().startswith("<b>Dfam"))
    assert "LTR/ERVK" in head.split("—")[0] and "U2" not in head and "ALR" not in head
    fam = next(t for t in _tables(root) if t.export_base == "TEagle_family")
    assert {r[2] for r in fam.rows_data()} == {"AluSx", "HERVK-int"}
    assert any("snRNA" in str(r) for t in _tables(root) if t is not fam for r in t.rows_data())


def test_family_table_headers_explain_themselves(win):
    win._render_family({"ok": True, "hits": [_hit("LTR/Copia", "Copia_I", 0, 500, 900)], "species": "Drosophila melanogaster",
                        "repeatmasker_version": "4.2", "include_uncurated": True})
    t = next(t for t in _tables(win.wslBody.parentWidget()) if t.export_base == "TEagle_family")
    tips = {t.horizontalHeaderItem(c).text(): t.horizontalHeaderItem(c).toolTip() for c in range(t.columnCount())}
    assert not tips["Score"].startswith("HMMER") and "Smith-Waterman" in tips["Score"]
    assert "Wicker 2007 scheme" not in tips["Class/family"]
    for h in ("Div (raw)", "Del", "Ins", "Consensus", "Cov", "Blocks"):
        assert tips[h], h


def test_defaulted_organism_is_not_presented_as_a_choice(win):
    win._render_family({"ok": True, "hits": [], "species": "Homo sapiens (RepeatMasker default lineage)",
                        "species_requested": None, "species_defaulted": True, "lineage_searched": "Homo sapiens",
                        "include_uncurated": False, "dfam_library": {"partitions": ["0 [root]: curated"]}})
    text = _texts(win.wslBody.parentWidget())
    assert "no organism was set" in text


def test_non_te_only_hits_still_carry_the_curated_only_caveat(win):
    win._render_family({"ok": True, "hits": [_hit("Simple_repeat", "(AT)n", 0, 40, 30)], "species": "Arabidopsis thaliana",
                        "species_defaulted": False, "include_uncurated": False,
                        "dfam_library": {"partitions": ["0 [root]: curated"]}})
    text = _texts(win.wslBody.parentWidget())
    assert "CURATED" in text and "of this locus" in text


# ---------------------------------------------------------------- splice panel (R19, R24, exports-14)
def test_splice_header_separates_alignment_orientation_from_transcript_strand(win):
    seq = "ACGT" * 60
    win.state["splice_seq"] = seq
    base = {"counts": {"exons": 1, "introns": 0}, "canonical_introns": 0, "exons": [{"start": 0, "end": 60}],
            "introns": []}
    win._render_splice({**base, "strand": "-", "alignment_strand": "-", "transcript_strand": None,
                        "strand_determined": False})
    text = _texts(win.spliceBody.parentWidget())
    assert "strand undetermined" in text and "aligned" in text
    win._render_splice({**base, "strand": "+", "alignment_strand": "-", "transcript_strand": "-",
                        "strand_determined": True})
    text = _texts(win.spliceBody.parentWidget())
    assert "strand undetermined" not in text and "reverse-complement" in text


def test_derived_introns_are_not_annotation_evidence_in_the_cross_check(win):
    seq = "ACGT" * 60
    win.state["splice_seq"] = seq
    win.state["source"] = {"accession": "J00265"}
    win.state["features"] = {"exons": [], "cds": [], "introns": [
        {"start": 60, "end": 90, "strand": "+", "derived": True},
        {"start": 120, "end": 150, "strand": "+"}]}
    win._render_splice({"counts": {"exons": 2, "introns": 1}, "canonical_introns": 1, "strand": "+",
                        "exons": [{"start": 0, "end": 120}, {"start": 150, "end": 200}],
                        "introns": [{"start": 120, "end": 150, "donor": "GT", "acceptor": "AG", "canonical": True}]})
    text = _texts(win.spliceBody.parentWidget())
    assert "1/1 intron" in text


def test_minus_strand_exon_fasta_is_numbered_in_transcript_order():
    import main
    out = main._splice_exons_fasta({"strand": "-", "exons": [{"start": 0, "end": 6}, {"start": 12, "end": 18}]},
                                   "AAACCCGTAAAGTTTGGG", "splice_locus")
    recs = [(h, s) for h, s in zip(out.splitlines()[0::2], out.splitlines()[1::2])]
    assert "_exon_1 " in recs[0][0] and "coords=12-18" in recs[0][0] and recs[0][1] == "CCCAAA"
    assert "_exon_2 " in recs[1][0] and "coords=0-6" in recs[1][0] and recs[1][1] == "GGGTTT"
    assert "".join(s for _, s in recs) == "CCCAAAGGGTTT"


# ---------------------------------------------------------------- FASTA headers (exports-17)
def test_fasta_headers_state_the_coordinate_convention_and_stay_ascii():
    import main
    rec = {"id": "r", "seq": "AAACCCGGGTTT" * 5,
           "orfs": [{"start": 0, "end": 6, "strand": "+", "frame": 1}],
           "domains": [{"domain": "RT", "nt": [1, 10], "strand": "+", "pfam": "", "protein": "MPEPTIDE"}]}
    for text in (main._orfs_fasta(rec), main._domains_fasta(rec)):
        assert "0-based-half-open" in text
        text.encode("ascii")


# ---------------------------------------------------------------- annotation export label (exports-16)
def test_anchored_export_names_the_accession_it_writes(win, monkeypatch):
    import main
    from PySide6.QtWidgets import QMenu
    rec = _accept(win, "M11240")
    labels = []
    monkeypatch.setattr(win, "_bound_coord_region", lambda r: {"chrAccession": "NC_000001.11", "chromLabel": "chr1",
                                                               "start": 10, "stop": 20, "strand": 1})
    from PySide6.QtCore import QTimer

    def grab():
        pop = QApplication.activePopupWidget()
        if isinstance(pop, QMenu):
            labels.extend(a.text() for a in pop.actions()); pop.close()
        else:
            QTimer.singleShot(50, grab)
    QTimer.singleShot(50, grab)
    win._export_annotation(rec)
    anchored = [l for l in labels if "anchored" in l]
    assert anchored and all("NC_000001.11" in l for l in anchored)


# ---------------------------------------------------------------- completeness wording (persona-U21, B10)
def test_partial_tier_is_explained_relative_to_the_bundled_panel(win):
    _accept(win, "M11240")
    text = _card(win)
    assert "Expected modules (2 of 5 found)" in text and "GAG not detected" in text
    assert "relative to the bundled" in text and "Gag_p24" in text and "RVP" in text


# ---------------------------------------------------------------- release self-test positive control (ui_family_splice_prov-9)
def test_selftest_positive_control_passes_on_a_healthy_build_and_fails_when_domains_break(monkeypatch):
    import main
    from teagle_core import domains
    assert main._selftest_positive_control() == []

    def broken(seq, max_orfs=domains.MAX_ORFS_SCANNED, evalue=1e-3):      # same signature: the seal reads its defaults
        raise RuntimeError("hmmsearch failed")
    monkeypatch.setattr(domains, "scan_domains", broken)
    problems = main._selftest_positive_control()
    assert any("domains stage" in p for p in problems) and any("RT or INT" in p for p in problems)


def test_text_export_to_a_locked_file_is_reported_not_raised(win, monkeypatch, tmp_path):
    import main, widgets
    seen = []
    monkeypatch.setattr(main.QFileDialog, "getSaveFileName", staticmethod(lambda *a, **k: (str(tmp_path / "x.fasta"), "")))

    def locked(*a, **k):
        raise PermissionError("locked")
    monkeypatch.setattr(widgets, "atomic_write", locked)
    monkeypatch.setattr(widgets, "_export_failed", lambda parent, msg: seen.append(msg))
    win._save_text(">a\nACGT", "x", "fasta")
    assert seen and "open in another program" in seen[0]


def test_record_transcript_fetch_does_not_block_the_window(win, monkeypatch):
    import threading, time
    import main
    gate = threading.Event()

    def slow(acc, **k):
        gate.wait(10)
        return {"fasta": ">NM_000207.3\nACGTACGT"}
    monkeypatch.setattr(main, "retrieve", slow)
    win.state["features"] = {"transcripts": [{"accession": "NM_000207.3", "product": "insulin"}]}
    win._refresh_transcripts()
    t0 = time.time()
    win._use_record_transcript()
    assert time.time() - t0 < 2 and not win.spliceTxUse.isEnabled()      # returned at once, button shows it is busy
    gate.set()
    deadline = time.time() + 10
    while time.time() < deadline and not win.spliceTx.toPlainText():
        _app.processEvents(); time.sleep(0.05)
    assert "ACGTACGT" in win.spliceTx.toPlainText() and win.spliceTxUse.isEnabled()
    assert win._splice_tx_origin == "record"


# ---------------------------------------------------------------- coordinate convention on copied text
def test_copied_feature_fasta_names_the_coordinate_convention_in_both_display_modes(win):
    import main
    rec = _accept(win, "M11240")
    for base in (1, 0):
        win._set_coord_base(base, persist=False)
        items = dict(win._feat_menu(10, 70, "-", "feat", src_seq=rec["seq"]))
        items["Copy FASTA"]()
        header = _clip().splitlines()[0]
        assert header == f">feat_10-70_rev {main._COORD_SYSTEM}"           # the file numbering, whatever the screen shows


def test_gel_menu_names_the_convention_and_its_fasta_carries_the_token(win):
    import main
    amp = {"source": "rec (input)", "start": 10, "end": 210, "length": 200, "seq": "ACGT" * 50,
           "matches_intended_product": True, "on_target": None}
    region = {"amplicon": amp, "pair": "D1·P2", "has_locus": False, "tip": "200 bp · 1 priming site"}
    win._set_coord_base(0, persist=False)
    items = dict(win._gel_menu(region))
    assert any(l.startswith("Copy coords") and "0-based, half-open" in l for l in items)
    next(f for l, f in items.items() if l.startswith("Copy FASTA"))()
    assert _clip().splitlines()[0].endswith(" " + main._COORD_SYSTEM)
    win._set_coord_base(1, persist=False)
    assert any(l.startswith("Copy coords") and "1-based" in l for l, _f in win._gel_menu(region))
