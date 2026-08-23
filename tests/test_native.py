"""Native (PySide6) app regression tests. Run headless via the offscreen Qt platform. Cover the
figure builders, the engine worker's three-tier taxonomy, and the analyze -> design -> in-silico PCR
workflow driven through the same in-process engine the GUI uses. Skipped if PySide6 is unavailable."""
import hashlib, json, os, sys, time

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
import pytest

pytest.importorskip("PySide6")
_NATIVE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "app", "native")
sys.path.insert(0, _NATIVE)

from PySide6.QtWidgets import QApplication
from PySide6.QtCore import QEventLoop, QTimer, QByteArray
from PySide6.QtGui import QImage, QPainter
from PySide6.QtSvg import QSvgRenderer

_app = QApplication.instance() or QApplication([])


def _spin(pred, timeout=30):
    loop = QEventLoop(); t = QTimer(); t.setInterval(25); t0 = time.time()
    t.timeout.connect(lambda: (pred() or time.time() - t0 > timeout) and loop.quit())
    t.start(); loop.exec()
    return pred()


def _renders(svg):
    r = QSvgRenderer(QByteArray(svg.encode("utf-8")))
    img = QImage(300, 200, QImage.Format_ARGB32); img.fill(0)
    p = QPainter(img); r.render(p); p.end()
    return r.isValid() and any(img.pixelColor(x, y).alpha() > 0
                               for x in range(0, 300, 15) for y in range(0, 200, 15))


# ---------- figure builders ----------
def test_gel_svg_renders():
    import figures
    svg = figures.svg_gel({"lanes": [{"label": "P1", "amplicons": [
        {"length": 250, "on_target": True, "source": "x"}, {"length": 80, "on_target": False, "source": "x"}]}]}, "dark")
    assert svg.startswith("<svg") and _renders(svg)


def test_gel_small_amplicon_axis_floor():
    # a 75 bp band must resolve — the log axis floor tracks the smallest band (obs 4680)
    import figures
    svg = figures.svg_gel({"lanes": [{"label": "P1", "amplicons": [{"length": 75, "on_target": True, "source": "x"}]}]}, "dark")
    assert "75 bp" in svg


def test_genome_svg_renders():
    import figures
    rec = {"composition": {"length": 5000},
           "structural": [{"type": "LTR (x)", "five_prime": [100, 300], "three_prime": [4700, 4900]}],
           "domains": [{"domain": "RT", "label": "RT", "nt": [1000, 2000], "score": 120}],
           "orfs": [{"start": 900, "end": 3000, "strand": "+", "frame": 0, "length_aa": 700}]}
    model = figures.gv_tracks_from_rec(rec)
    assert len(model["tracks"]) == 3
    assert _renders(figures.svg_genome(model, {"start": 0, "end": 5000}, 620, "dark"))


def test_png_export_has_alpha(tmp_path):
    import figures
    from widgets import render_png
    svg = figures.svg_gel({"lanes": [{"label": "P1", "amplicons": [{"length": 200, "on_target": True, "source": "x"}]}]}, "transparent")
    out = str(tmp_path / "gel.png")
    render_png(svg, out, scale=2)
    img = QImage(out)
    assert not img.isNull() and img.hasAlphaChannel()


# ---------- engine worker taxonomy ----------
def test_worker_success_and_user_error():
    from engine_worker import Engine
    eng = Engine()
    got = {}
    eng.done.connect(lambda k, r: got.__setitem__(k, ("done", r)))
    eng.user_error.connect(lambda k, m: got.__setitem__(k, ("user_error", m)))
    eng.failed.connect(lambda k, m, t: got.__setitem__(k, ("failed", m)))
    eng.submit("analyze", {"sequence": "ACGTACGTACGTACGTACGT"}, key="ok")
    eng.submit("pcr", {"sequence": "ACGT" * 40}, key="bad")     # missing primers -> BadRequest
    _spin(lambda: "ok" in got and "bad" in got)
    assert got["ok"][0] == "done" and "records" in got["ok"][1]
    assert got["bad"][0] == "user_error"


# ---------- full workflow ----------
@pytest.fixture
def win():
    import main
    w = main.MainWindow()
    yield w
    w.close()


def _accept_analysis(win, result, sequence=None):
    sequence = sequence or "\n".join(f">{r['id']}\n{r['seq']}" for r in result.get("records", []))
    win.seq.setPlainText(sequence)
    win._analyze_gen += 1
    generation = win._analyze_gen
    win._active_analyze_gen = generation
    win._analysis_requests[generation] = {"sequence": sequence, "source": {}}
    win._analyze_inflight = True
    win._on_analyze(result, generation)


def test_rail_toggle_collapses_and_restores(win):
    # the specimen panel collapses to free width and reopens, restoring the split (state tracked explicitly,
    # since isVisible() is unreliable before the window is shown)
    assert win._rail_collapsed is False and win.railToggle.text().upper() == "HIDE"
    win._rail_sizes = [340, 900]
    win._toggle_rail()
    assert win._rail_collapsed is True and win.rail.isHidden() is True and win.railToggle.text().upper() == "SHOW"
    win._toggle_rail()
    assert win._rail_collapsed is False and win.rail.isHidden() is False and win.split.sizes()[0] > 0


def test_ui_scale_applied_from_settings(monkeypatch):
    import main, theme as theme_mod
    monkeypatch.delenv("TEAGLE_UI_SCALE", raising=False)
    class _FakeSettings:
        def __init__(self, *a): pass
        def value(self, k, d=None): return 1.25
    monkeypatch.setattr(main, "QSettings", _FakeSettings)
    orig = theme_mod.UI_SCALE
    try:
        main._apply_saved_ui_scale()
        assert abs(theme_mod.UI_SCALE - 1.25) < 1e-6             # persisted scale drives theme_mod.UI_SCALE live (no QT_SCALE_FACTOR)
        # an explicit TEAGLE_UI_SCALE env override must win (function returns without changing UI_SCALE)
        theme_mod.UI_SCALE = 2.0
        monkeypatch.setenv("TEAGLE_UI_SCALE", "2.0")
        main._apply_saved_ui_scale()
        assert theme_mod.UI_SCALE == 2.0
    finally:
        theme_mod.UI_SCALE = orig


def test_methods_html_states_databases_and_params(win):
    # the methodology surface must explicitly name the database/model/consensus + thresholds for each layer
    html = win._methods_html()
    assert "HMMER" in html and "Pfam" in html and "PF00078" in html and "PF00665" in html   # domains: HMMER vs Pfam profiles
    assert "1e-3" in html                                                                     # domain E-value threshold
    assert "k=13" in html and "80%" in html                                                   # structural heuristic params
    assert "Wicker" in html                                                                   # classification scheme
    assert "Dfam" in html and "RepeatMasker" in html                                          # family naming (WSL)


def test_ui_scale_button_present(win):
    from PySide6.QtWidgets import QPushButton
    assert any("SCALE" in b.text() for b in win.findChildren(QPushButton))


def test_feature_subregion_routes_subset(win, monkeypatch):
    # a feature right-click offers a sub-region picker; a chosen 1-based-inclusive interval routes only that
    # subset (offset into the feature's own sequence) to primer design / splice
    from PySide6.QtWidgets import QSpinBox, QPushButton, QDialog, QApplication
    dna = "ATCGATCGATCGATCGATCGGGGGCCCCCTTTTTAAAAAACGTACGTACGT"     # 51 bp
    labels = [l for l, _ in win._feat_menu(0, len(dna), "+", "feat", dna=dna)]
    assert any("sub-region" in l for l in labels)
    calls = {}
    monkeypatch.setattr(win, "_design_for_domain", lambda s, e, lab, seq=None: calls.setdefault("design", (s, e, lab, seq)))
    monkeypatch.setattr(win, "_send_to_splice", lambda fasta: calls.setdefault("splice", fasta))
    monkeypatch.setattr(QDialog, "exec", lambda self: 0)
    win._subregion(dna, "feat")
    dlg = [w for w in QApplication.topLevelWidgets()
           if isinstance(w, QDialog) and w.windowTitle() == "Select a sub-region"][-1]
    spins = dlg.findChildren(QSpinBox)
    spins[0].setValue(10); spins[1].setValue(30)                     # bases 10..30 -> dna[9:30]
    [b for b in dlg.findChildren(QPushButton) if "DESIGN" in b.text().upper()][0].click()
    assert calls["design"][3] == dna[9:30] and calls["design"][2] == "feat_10-30"


def test_coord_fetch_submits_and_renders(win):
    # curated organism -> fetch_coords op with the right body; render clears accMeta, sets coordMeta + source
    ops = []
    win.engine.submit = lambda op, body=None, key=None: ops.append((op, body, key))
    win.asmSel.setCurrentIndex(win.asmSel.findData("Homo sapiens"))
    win.coord.setPlainText("chr13:33,016,423-33,066,143")
    win.coordStrand.setCurrentIndex(1)                        # minus strand
    win._fetch_coord()
    assert ops[0] == ("fetch_coords", {"regions": "chr13:33,016,423-33,066,143", "strand": "-",
                                       "organism": "Homo sapiens", "customQuery": ""}, "fetch")
    from teagle_core import provenance
    sequence = "A" * 100
    raw_source = {
        "accession": "NC_000013.11", "assemblyAccession": "GCF_000001405.40", "taxid": "9606",
        "regions": [{"chrAccession": "NC_000013.11", "start": 1, "stop": 100, "strand": 2}],
        "coordSystem": "1-based-inclusive", "retrievalType": "coordinate",
        "source": "NCBI E-utilities (efetch seq_start/seq_stop)",
        "endpoint": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
        "retrievedUtc": "2026-08-21T00:00:00+00:00",
    }
    fasta = ">NC_000013.11:1-100(-)\n" + sequence + "\n"
    source_context = provenance.normalize_source_context(
        raw_source, [("NC_000013.11:1-100(-)", sequence)])
    res = {"ok": True, "runType": "coordinate", "fasta": fasta, "organism": "Homo sapiens",
           "assemblyName": "GRCh38.p14", "displayLocus": "chr13:1-100",
           "regions": [{"chrAccession": "NC_000013.11", "start": 1, "stop": 100, "strand": 2, "chromLabel": "chr13"}],
           "source": raw_source, "sourceContext": source_context}
    win._on_fetch(res)
    assert "chr13" in win.coordMeta.text() and win.accMeta.text() == ""
    assert win.state["source"] == source_context and "displayLocus" not in win.state["source"]
    assert win.state["coord_region"]["_record_id"] == "NC_000013.11:1-100"
    assert win._bound_coord_region({"id": "NC_000013.11:1-100", "seq": sequence}) == {}  # minus-strand is fail closed
    win._set_seq(">new\nACGT")
    assert win.state["coord_region"] is None


def test_plus_coordinate_binding_requires_exact_record_hash_and_source(win):
    from teagle_core import provenance
    accession, sequence = "NC_000001.11", "ACGTACGT"
    record_id = f"{accession}:5-12"
    raw_source = {
        "accession": accession, "assemblyAccession": "GCF_000001405.40", "taxid": "9606",
        "regions": [{"chrAccession": accession, "start": 5, "stop": 12, "strand": 1}],
        "coordSystem": "1-based-inclusive", "retrievalType": "coordinate",
        "source": "NCBI E-utilities (efetch seq_start/seq_stop)",
        "endpoint": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
        "retrievedUtc": "2026-08-21T00:00:00+00:00",
    }
    source_context = provenance.normalize_source_context(raw_source, [(record_id, sequence)])
    res = {"ok": True, "runType": "coordinate", "fasta": f">{record_id}\n{sequence}\n", "organism": "test",
           "assemblyName": "asm", "regions": [{"chrAccession": accession, "start": 5, "stop": 12,
                                                  "strand": 1, "chromLabel": "chr1"}],
           "source": raw_source, "sourceContext": source_context}
    win._on_fetch(res)
    rec = {"id": record_id, "identitySource": "validated_fetch", "seq": sequence}
    assert win._bound_coord_region(rec)["start"] == 5
    assert win._bound_coord_region({**rec, "seq": "ACGTACGA"}) == {}
    assert win._bound_coord_region({**rec, "id": "other"}) == {}
    win.state["source"] = {**win.state["source"], "accession": "OTHER"}
    assert win._bound_coord_region(rec) == {}


def test_coord_custom_organism_path(win):
    ops = []
    win.engine.submit = lambda op, body=None, key=None: ops.append((op, body, key))
    win.asmSel.setCurrentIndex(win.asmSel.findData("__custom__"))
    win.coordCustom.setText("GCF_000001405.40"); win.coord.setPlainText("chr7:1-100")
    win._fetch_coord()
    assert ops[0][1]["customQuery"] == "GCF_000001405.40" and ops[0][1]["organism"] == ""
    # empty custom field -> banner, no submit
    win.coordCustom.setText(""); ops.clear(); win._fetch_coord()
    assert not ops


def test_coord_meta_cleared_only_when_in_flight(win):
    # a failed coordinate fetch must clear the 'fetching…' indicator, not leave it stuck
    win.coordMeta.setText("fetching…")
    win._on_failed("fetch", "boom", "trace")
    assert win.coordMeta.text() == ""
    win.coordMeta.setText("fetching…")
    win._on_user_error("fetch", "bad input")
    assert win.coordMeta.text() == ""
    # a prior successful result on the OTHER panel must survive a failed fetch
    win.accMeta.setText("M18706.1 · Tnt1"); win.coordMeta.setText("fetching…")
    win._on_failed("fetch", "boom", "trace")
    assert win.coordMeta.text() == "" and win.accMeta.text() == "M18706.1 · Tnt1"


def test_analyze_to_pcr_workflow(win):
    win._load_sample()
    win._run_analysis()
    assert _spin(lambda: win.state.get("analyzed_clean")), "analyze did not complete"
    assert win.designBtn.isEnabled()
    assert win.mValid.text() == "valid"

    done = {"p": False}
    orig = win._on_primers
    win._on_primers = lambda d: (orig(d), done.__setitem__("p", True))
    win._design()
    assert _spin(lambda: done["p"]), "primer design did not complete"
    cands = win.state.get("candidates", [])
    assert cands, "no primer candidates"

    win._pcr_stage_all()
    assert len(win.state["pcrPairs"]) == len(cands)
    assert win.runPcrBtn.isEnabled()

    orig_pcr = win._pcr_slot
    finished = {"v": False}
    def wrap(k, d):
        orig_pcr(k, d)
        if win._pcr_run and all(r is not None for r in win._pcr_run["results"]):
            finished["v"] = True
    win._pcr_slot = wrap
    win._run_pcr()
    assert _spin(lambda: finished["v"], 40), "in-silico PCR did not complete"
    assert len(win.state["lastPcr"]["lanes"]) == len(cands)
    # A general design has an intended endpoint pair but no Primer3 target, so target status is explicitly unassessed.
    assert any(a.get("matches_intended_product") is True for a in win.state["lastPcr"]["amplicons"])
    assert all(a.get("on_target") is None for a in win.state["lastPcr"]["amplicons"]
               if not a.get("single_primer"))
    assert win.card_prov.bodylay.count() >= 1        # provenance rendered


def test_staleness_guard_blocks_stale_primer(win):
    win._load_sample()
    win._run_analysis()
    assert _spin(lambda: win.state.get("analyzed_clean"))
    win.seq.setPlainText("ACGTACGTACGTACGTACGTACGTACGTACGT")   # change sequence after analysis
    assert win._stale_block() is True                          # design must be blocked
    assert win.state.get("analysis") is None
    assert win.state.get("prov_manifest") is None


def test_analysis_callback_generation_drops_out_of_order_result(win):
    import engine
    submitted = []
    win.engine.submit = lambda op, body=None, key=None: submitted.append((op, body, key))
    win.seq.setPlainText(">old\n" + "ACGT" * 20)
    win._run_analysis()
    old_key = submitted[-1][2]
    win.seq.setPlainText(">new\n" + "TGCA" * 20)
    win._run_analysis()
    new_key = submitted[-1][2]
    win._on_done(old_key, engine.run_analyze({"sequence": ">old\n" + "ACGT" * 20}))
    assert win.state.get("analysis") is None and win.state.get("prov_manifest") is None
    win._on_done(new_key, engine.run_analyze({"sequence": ">new\n" + "TGCA" * 20}))
    from teagle_core.gff3 import public_record_id
    assert win.state["analysis"]["records"][0]["id"] == public_record_id("new")
    assert win.state["prov_manifest"]["runType"] == "analysis"


def test_stale_analysis_failure_cannot_clear_newer_success(win):
    import engine
    submitted = []
    win.engine.submit = lambda op, body=None, key=None: submitted.append((op, body, key))
    win.seq.setPlainText(">old\n" + "ACGT" * 20); win._run_analysis()
    old_key = submitted[-1][2]
    win.seq.setPlainText(">new\n" + "TGCA" * 20); win._run_analysis()
    new_key = submitted[-1][2]
    win._on_done(new_key, engine.run_analyze({"sequence": ">new\n" + "TGCA" * 20}))
    seal = win.state["prov_manifest"]["manifestSha256"]
    win._on_user_error(old_key, "old request failed")
    assert win.state["prov_manifest"]["manifestSha256"] == seal
    assert win.state["analysis_state"] == "current"


def test_new_analysis_and_failure_clear_prior_primer_manifest(win):
    win.state["prov_manifest"] = {"runType": "primer", "manifestSha256": "f" * 64}
    win.state["prov_manifest_all"] = [win.state["prov_manifest"]]
    submitted = []
    win.engine.submit = lambda op, body=None, key=None: submitted.append((op, body, key))
    win.seq.setPlainText(">x\n" + "ACGT" * 20)
    win._run_analysis()
    assert win.state.get("prov_manifest") is None
    win._on_user_error(submitted[-1][2], "bad input")
    assert win.state.get("analysis") is None and win.state.get("prov_manifest_all") is None


def test_analysis_result_json_export_is_current_and_bound(win, tmp_path, monkeypatch):
    import engine, main
    result = engine.run_analyze({"sequence": ">one\n" + "ACGT" * 20 + "\n>two\n" + "TGCA" * 20})
    _accept_analysis(win, result)
    out = str(tmp_path / "analysis.json")
    proposed = []
    monkeypatch.setattr(main.QFileDialog, "getSaveFileName",
                        staticmethod(lambda *a, **k: (proposed.append(a[2]) or out, "")))
    win._export_analysis_result()
    assert proposed == ["TEagle_multiple_records_analysis.json"]
    payload = json.loads(open(out, encoding="utf-8").read())
    assert engine.validate_analysis_result_document(payload)
    from teagle_core.gff3 import public_record_id
    assert [r["id"] for r in payload["records"]] == [public_record_id("one"), public_record_id("two")]


def test_manifest_export_reverifies_every_batch_seal_at_click_time(win, tmp_path, monkeypatch):
    import copy, engine, main

    result = engine.run_analyze({"sequence": ">one\n" + "ACGT" * 20})
    _accept_analysis(win, result)
    valid = copy.deepcopy(result["provenance"])
    win.state["prov_manifest_all"] = [valid, copy.deepcopy(valid)]
    out = tmp_path / "manifests.json"
    monkeypatch.setattr(main.QFileDialog, "getSaveFileName", staticmethod(lambda *a, **k: (str(out), "")))
    win._export_manifest()
    payload = json.loads(out.read_text(encoding="utf-8"))
    assert len(payload["runs"]) == 2

    original = out.read_bytes()
    win.state["prov_manifest_all"][1]["parameters"]["orf_min_aa"] += 1
    monkeypatch.setattr(win, "_banner", lambda *a, **k: None)
    win._export_manifest()
    assert out.read_bytes() == original
    assert not (tmp_path / "manifests.json.part").exists()


@pytest.mark.parametrize("mutation", ["fake_seal", "excluded_nan"])
def test_landscape_manifest_serializer_preserves_destination_on_click_time_failure(
        win, tmp_path, monkeypatch, mutation):
    import copy, engine, main

    manifest = copy.deepcopy(engine.run_analyze({"sequence": ">one\n" + "ACGT" * 20})["provenance"])
    if mutation == "fake_seal":
        manifest["manifestSha256"] = "0" * 64
    else:
        manifest.setdefault("environment", {})["invalid"] = float("nan")
    out = tmp_path / "landscape.json"
    out.write_text("preserve\n", encoding="utf-8")
    monkeypatch.setattr(main.QFileDialog, "getSaveFileName",
                        staticmethod(lambda *a, **k: (str(out), "")))
    monkeypatch.setattr(win, "_banner", lambda *a, **k: None)
    win._export_manifest_payload([manifest], "TEagle_TE_landscape_safe_manifest")
    assert out.read_text(encoding="utf-8") == "preserve\n"
    assert not (tmp_path / "landscape.json.part").exists()


def test_analysis_result_export_reverifies_payload_at_click_time(win, tmp_path, monkeypatch):
    import engine, main

    result = engine.run_analyze({"sequence": ">one\n" + "ACGT" * 20})
    _accept_analysis(win, result)
    out = tmp_path / "analysis.json"
    out.write_text("validated existing result\n", encoding="utf-8")
    monkeypatch.setattr(main.QFileDialog, "getSaveFileName", staticmethod(lambda *a, **k: (str(out), "")))
    monkeypatch.setattr(win, "_banner", lambda *a, **k: None)
    win.state["analysis"]["records"][0]["composition"]["gc"] = 99.0
    win._export_analysis_result()
    assert out.read_text(encoding="utf-8") == "validated existing result\n"
    assert not (tmp_path / "analysis.json.part").exists()


def test_analysis_callback_rejects_tampered_provenance_binding(win):
    import engine
    result = engine.run_analyze({"sequence": ">x\n" + "ACGT" * 20})
    result["provenance"]["parameters"]["orf_min_aa"] += 1
    _accept_analysis(win, result)
    assert win.state.get("analysis") is None
    assert win.state.get("prov_manifest") is None


def test_untracked_legacy_analysis_callbacks_fail_closed(win):
    import engine
    result = engine.run_analyze({"sequence": ">x\n" + "ACGT" * 20})
    win._on_done("analyze", result)
    assert win.state.get("analysis") is None
    win._on_analyze(result)
    assert win.state.get("analysis") is None


@pytest.mark.parametrize("trigger", ["input_edit", "new_request", "failure", "cancel"])
def test_stale_rendered_exports_cannot_open_after_invalidation(win, monkeypatch, trigger):
    import engine, main
    from PySide6.QtCore import QEvent, QPoint, Qt
    from PySide6.QtWidgets import QPushButton
    from widgets import DataTable
    sequence = ">x\n" + "ACGT" * 30
    _accept_analysis(win, engine.run_analyze({"sequence": sequence}), sequence)
    downstream = (win.primBody, win.pcrBody, win.wslBody, win.spliceBody, win.genomeBody)
    helper_controls = []
    for index, layout in enumerate(downstream):
        table = DataTable(["Value"]); table.set_rows([[index]])
        button = main._export_table_btn(table, f"stale_{index}", win)
        layout.addWidget(table); layout.addWidget(button)
        helper_controls.append((table, button))
    roots = (win.card_struct, win.card_prov) + tuple(layout.parentWidget() for layout in downstream)
    stale_buttons = list(dict.fromkeys(b for root in roots
                     for b in root.findChildren(QPushButton)
                     if "EXPORT" in b.text().upper() or b.text().upper() in {"SVG", "PNG", "PDF"}))
    stale_tables = list(dict.fromkeys(t for root in roots for t in root.findChildren(DataTable)))
    assert all(button in stale_buttons and table in stale_tables for table, button in helper_controls)
    assert {"SVG", "PNG"}.issubset({button.text().upper() for button in stale_buttons})
    dialogs = []
    monkeypatch.setattr(main.QFileDialog, "getSaveFileName",
                        staticmethod(lambda *a, **k: (dialogs.append(a) or ("", ""))))
    monkeypatch.setattr(main.widgets, "pick_table_format", lambda *a, **k: (dialogs.append(a) or None))
    monkeypatch.setattr(win, "_export_annotation", lambda *a, **k: dialogs.append(a))
    submitted = []
    win.engine.submit = lambda op, body=None, key=None: submitted.append((op, body, key))
    if trigger == "input_edit":
        win.seq.setPlainText(">edited\n" + "TGCA" * 30)
    else:
        win._run_analysis()
        if trigger == "failure":
            win._on_user_error(submitted[-1][2], "failed")
        elif trigger == "cancel":
            win.seq.setPlainText(">cancelled\n" + "TGCA" * 30)
    for button in stale_buttons:
        assert not button.isEnabled()
        button.click()
    for table in stale_tables:
        assert not table.isEnabled()
        assert table.contextMenuPolicy() == Qt.NoContextMenu
        table.customContextMenuRequested.emit(QPoint(0, 0))
    assert dialogs == []
    QApplication.sendPostedEvents(None, QEvent.DeferredDelete)
    for button in stale_buttons:
        try:
            button.click()
        except RuntimeError:
            pass
    assert dialogs == []
    assert win.state.get("analysis") is None and win.state.get("prov_manifest") is None
    assert not win.card_struct.findChildren(DataTable)
    assert not [b for card in (win.card_struct, win.card_prov)
                for b in card.findChildren(QPushButton)
                if "EXPORT" in b.text().upper() or b.text().upper() in {"SVG", "PNG", "PDF"}]


def test_accession_fetch_preserves_bounded_source_metadata(win):
    accession = "M11240.1"
    win._on_fetch({"ok": True, "fasta": f">{accession}\nACGT", "accession": accession, "organism": "test",
                   "taxid": "9", "source": "NCBI nuccore",
                   "endpoint": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
                   "retrievedUtc": "2026-08-21T00:00:00+00:00",
                   "sourceUrl": f"https://www.ncbi.nlm.nih.gov/nuccore/{accession}", "length": 4})
    assert {k: win.state["source"][k] for k in
            ("sourceContextVersion", "inputType", "accession", "taxid", "source", "endpoint", "retrievedUtc")} == {
        "sourceContextVersion": "1.0", "inputType": "fetched_accession", "accession": accession,
        "taxid": "9", "source": "ncbi_nuccore", "endpoint": "ncbi_efetch",
        "retrievedUtc": "2026-08-21T00:00:00+00:00"}
    assert win.state["source"]["records"] == [{"id": accession, "length": 4,
                                                  "sha256": hashlib.sha256(b"ACGT").hexdigest()}]


def test_native_operation_input_serializes_exact_selected_fetch_binding(win):
    from teagle_core import provenance

    accession = "M11240.1"
    sequence = "ACGT" * 10
    raw = {"accession": accession, "taxid": "7227", "source": "NCBI nuccore",
           "endpoint": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
           "retrievedUtc": "2026-08-21T00:00:00+00:00",
           "sourceUrl": f"https://www.ncbi.nlm.nih.gov/nuccore/{accession}", "length": len(sequence)}
    source = provenance.normalize_source_context(raw, [(accession, sequence)])
    win.state.update(seq=sequence, source=source, last_rec={"id": accession})
    operation_sequence, operation_source = win._operation_input()
    assert operation_sequence == f">{accession}\n{sequence}"
    assert provenance.validate_source_context(
        json.loads(json.dumps(operation_source)), [(accession, sequence)]) == operation_source


def test_native_operation_input_narrows_coordinate_wire_to_selected_region(win):
    from teagle_core import provenance

    sequences = ["ACGT", "TGCA"]
    raw = {
        "accession": "NC_000001.11", "assemblyAccession": "GCF_000001405.40", "taxid": "9606",
        "regions": [
            {"chrAccession": "NC_000001.11", "start": 1, "stop": 4, "strand": 1},
            {"chrAccession": "NC_000002.12", "start": 9, "stop": 12, "strand": 2},
        ],
        "coordSystem": "1-based-inclusive", "retrievalType": "coordinate",
        "source": "NCBI E-utilities (efetch seq_start/seq_stop)",
        "endpoint": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
        "retrievedUtc": "2026-08-21T00:00:00+00:00",
    }
    records = [("NC_000001.11:1-4", sequences[0]), ("NC_000002.12:9-12(-)", sequences[1])]
    source = provenance.normalize_source_context(raw, records)
    selected_id = "NC_000002.12:9-12"
    win.state.update(seq=sequences[1], source=source, last_rec={"id": selected_id})
    operation_sequence, operation_source = win._operation_input()
    assert operation_sequence == f">{selected_id}\n{sequences[1]}"
    assert operation_source["records"] == [source["records"][1]]
    assert operation_source["regions"] == [source["regions"][1]]
    assert operation_source["accession"] == "NC_000002.12"
    assert provenance.validate_source_context(operation_source, [(selected_id, sequences[1])]) == operation_source


def test_multi_record_analysis_does_not_false_trigger_stale_guard(win):
    """A 2+-record FASTA must NOT self-block Design/PCR. Regression: _show_record once overwrote
    state['analyzed_clean'] with the selected record's bases, so the whole-box _stale_block() comparison
    always mismatched for multi-record input and permanently fired the 'sequence changed' banner."""
    from teagle_core import examples
    fa = examples.load("M11240") + "\n" + examples.load("X01005")   # two records, no edits after analysis
    win.seq.setPlainText(fa)
    win._run_analysis()
    assert _spin(lambda: win.state.get("analyzed_clean")), "analyze did not complete"
    assert len(win.state.get("records") or []) == 2
    assert win._stale_block() is False, "record-0 view falsely reads as stale on a multi-record analysis"
    win._show_record(1)                                           # select the second record
    assert win._stale_block() is False, "selecting record 2 falsely trips the stale guard"
    assert win.designBtn.isEnabled()
    # downstream design must now target the SELECTED record's own sequence (the original CRITICAL)
    assert win.state["seq"] == win.state["records"][1]["seq"]


def test_single_record_gets_common_summary_and_full_record_fasta_surface(win):
    import engine, main
    from widgets import DataTable
    from PySide6.QtWidgets import QPushButton

    result = engine.run_analyze({"sequence": ">single specimen\n" + "ACGT" * 30})
    _accept_analysis(win, result)
    summary = next(table for table in win.card_struct.findChildren(DataTable)
                   if table.export_base.endswith("_summary"))
    assert summary.export_base == f"TEagle_{main._record_export_stem(result['records'][0])}_summary"
    assert summary.rowCount() == 1
    assert summary._headers[-1] == "Processing"
    processing = summary.rows_data()[0][-1]
    assert all(f"{stage}=" in processing for stage in engine.ANALYSIS_STAGE_NAMES)
    labels = {button.text().upper() for button in win.card_struct.findChildren(QPushButton)}
    assert "EXPORT FULL RECORDS (FASTA)" in labels and "EXPORT TABLE" in labels


def test_rewritten_record_filename_uses_full_digest_not_collision_prone_prefix():
    import main
    from teagle_core.gff3 import public_record_identity

    first = public_record_identity("patient_042")
    second = {"id": first["id"], "idSha256": first["idSha256"][:12] + "f" * 52,
              "identitySource": "local_input"}
    assert main._record_export_stem({**first, "identitySource": "local_input"}) == f"record_{first['idSha256']}"
    assert main._record_export_stem(second) != main._record_export_stem(
        {**first, "identitySource": "local_input"})


def test_fasta_rejects_resealed_compact_full_digest_mismatch():
    import main
    from teagle_core.gff3 import public_record_identity

    identity = public_record_identity("patient_042")
    with pytest.raises(ValueError, match="does not match"):
        main._fasta_text([(identity["id"], "ACGT", "", "0" * 64)], identifiers_are_public=True)


def test_bulk_fasta_builders_preserve_feature_sequence_and_strand():
    import main
    rec = {
        "id": "locus one",
        "seq": "AAACCCGGGTTT",
        "orfs": [{"start": 0, "end": 6, "strand": "+", "frame": 1},
                 {"start": 6, "end": 12, "strand": "-", "frame": 2}],
        "domains": [{"domain": "RT", "nt": [1, 10], "strand": "+", "pfam": "PF00078",
                     "protein": "MPEPTIDE"}],
    }
    orfs = main._orfs_fasta(rec)
    from teagle_core.gff3 import public_record_id
    rid = public_record_id(rec["id"])
    from teagle_core.gff3 import _identity_sha256
    digest = _identity_sha256(rec["id"])
    assert f">{rid}_ORF_1 coords=0-6 strand=+ frame=1 record_sha256={digest}\nAAACCC" in orfs
    assert f">{rid}_ORF_2 coords=6-12 strand=- frame=2 record_sha256={digest}\nAAACCC" in orfs
    assert "MPEPTIDE" in main._domains_fasta(rec) and digest in main._domains_fasta(rec)
    exons = main._splice_exons_fasta(
        {"strand": "-", "exons": [{"start": 0, "end": 3}, {"start": 9, "end": 12}]},
        rec["seq"], "splice")
    splice_id = public_record_id("splice")
    splice_digest = _identity_sha256("splice")
    assert f">{splice_id}_exon_1 coords=0-3 strand=- record_sha256={splice_digest}\nTTT" in exons
    assert f">{splice_id}_exon_2 coords=9-12 strand=- record_sha256={splice_digest}\nAAA" in exons


@pytest.mark.parametrize("record_id", [
    r"C:\Users\person\private\sample.fa",
    "/home/person/private/sample.fa",
    "../../private/sample.fa",
    "patient_042",
    "donor-A1",
    "sample_tuna_birgun",
    "descriptive specimen",
])
def test_native_summary_and_sequence_exports_never_repeat_local_record_ids(win, record_id):
    import engine, main
    from teagle_core.gff3 import public_record_id

    result = engine.run_analyze({"sequence": f">{record_id}\n" + "ACGT" * 30})
    public_id = public_record_id(record_id)
    assert result["records"][0]["id"] == public_id
    summary = win._record_summary_rows(result["records"])
    assert summary[0][1] == public_id and record_id not in repr(summary)
    rec = {"id": record_id, "seq": "ATG" + "AAA" * 45 + "TAA",
           "orfs": [{"start": 0, "end": 141, "strand": "+", "frame": 1}],
           "domains": [{"domain": "RT", "protein": "MPEPTIDE", "nt": [0, 24], "strand": "+"}]}
    outputs = [main._fasta_text([(record_id, "ACGT", "full record")]),
               main._fasta_text([(public_id, "ACGT", "full record",
                                  result["records"][0]["idSha256"])]),
               main._orfs_fasta(rec), main._domains_fasta(rec),
               main._splice_exons_fasta({"strand": "+", "exons": [{"start": 0, "end": 4}]},
                                         "ACGT", record_id)]
    for output in outputs:
        assert record_id not in output and public_id in output
        assert result["records"][0]["idSha256"] in output


def test_native_summary_and_fasta_preserve_public_accession(win):
    import engine, main
    from teagle_core import provenance

    sequence = "ACGT" * 30
    raw_source = {"accession": "M11240.1", "source": "NCBI nuccore",
        "endpoint": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
        "retrievedUtc": "2026-08-21T00:00:00+00:00",
        "sourceUrl": "https://www.ncbi.nlm.nih.gov/nuccore/M11240.1", "length": len(sequence)}
    source = provenance.normalize_source_context(raw_source, [("M11240.1", sequence)])
    result = engine.run_analyze({"sequence": ">M11240.1\n" + sequence, "source": source})
    assert win._record_summary_rows(result["records"])[0][1] == "M11240.1"
    fasta = main._full_records_fasta(result["records"])
    assert fasta.startswith(">M11240.1 ") and "record_sha256=" not in fasta


def test_native_pcr_and_genome_scan_tables_use_public_sources(win):
    from teagle_core.gff3 import public_record_id
    from widgets import DataTable

    raw = r"C:\Users\person\private\template.fa"
    public_id = public_record_id(raw)
    amp = {"pair": "P1", "mechanism_key": "F+R", "source": raw, "start": 1, "end": 181,
           "length": 180, "product_length": 180, "fwd_mm": 0, "rev_mm": 0, "on_target": True,
           "matches_intended_product": True, "seq": "A" * 180}
    win._render_pcr([{"label": "P1", "amplicons": [amp], "has_locus": True}], [amp],
                    [{"pair": "P1", "status": "complete", "error": None}])
    pcr_table = next(table for table in win.pcrBody.parentWidget().findChildren(DataTable)
                     if table.export_base == "TEagle_amplicons")
    assert pcr_table.rows_data()[0][2] == public_id
    win._render_genome_scan([{"label": "genome", "amplicons": [amp], "has_locus": True}], [amp],
                            {"verdict": "fixture", "has_locus": True, "n_on": 1, "n_off": 0,
                             "n_single": 0, "n_pair": 1, "per_source": []})
    scan_table = next(table for table in win.genomeBody.parentWidget().findChildren(DataTable)
                      if table.export_base == "TEagle_offtarget_scan")
    assert scan_table.rows_data()[0][1] == public_id
    assert raw not in repr(pcr_table.rows_data()) + repr(scan_table.rows_data())


def test_no_te_repeat_exports_use_the_two_exact_displayed_tables(win, monkeypatch):
    import main
    from widgets import DataTable
    from PySide6.QtWidgets import QPushButton

    hits = [
        {"class_family": "Simple_repeat", "family": "(AT)n", "q_start": 2, "q_end": 20,
         "strand": "+", "divergence": 1.0, "score": 40},
        {"class_family": "Unknown", "family": "Mystery", "q_start": 30, "q_end": 60,
         "strand": "-", "divergence": 4.0, "score": 55},
    ]
    win._render_family({"hits": hits, "species": "test"})
    root = win.wslBody.parentWidget()
    tables = root.findChildren(DataTable)
    assert len(tables) == 2
    assert {table.export_base for table in tables} == {"TEagle_non_TE_repeats", "TEagle_unclassified_repeats"}
    assert all(table.rowCount() == 1 and table.export_notes for table in tables)
    captured = []
    monkeypatch.setattr(main.widgets, "pick_table_format", lambda *a, **k: "csv")
    monkeypatch.setattr(main.widgets, "export_table",
                        lambda h, r, b, p, fmt=None, notes=None: captured.append((h, r, b, notes)))
    buttons = [button for button in root.findChildren(QPushButton) if button.text().upper() == "EXPORT TABLE"]
    assert len(buttons) == 2
    for button in buttons:
        button.click()
    expected = {(table.export_base, tuple(tuple(row) for row in table.rows_data()), tuple(table.export_notes))
                for table in tables}
    observed = {(base, tuple(tuple(row) for row in rows), tuple(notes))
                for _headers, rows, base, notes in captured}
    assert observed == expected


def test_mixed_te_repeat_result_renders_and_exports_all_three_exact_partitions(win, monkeypatch):
    import main
    from widgets import DataTable
    from PySide6.QtWidgets import QPushButton

    hits = [
        {"class_family": "LTR/Gypsy", "family": "Gypsy", "q_start": 0, "q_end": 100,
         "strand": "+", "divergence": 2.0, "score": 120},
        {"class_family": "Simple_repeat", "family": "(AT)n", "q_start": 110, "q_end": 130,
         "strand": "+", "divergence": 1.0, "score": 40},
        {"class_family": "Unknown", "family": "Mystery", "q_start": 140, "q_end": 160,
         "strand": "-", "divergence": 4.0, "score": 55},
    ]
    win._render_family({"hits": hits, "species": "test", "repeatmasker_version": "4.2",
                        "include_uncurated": False})
    root = win.wslBody.parentWidget()
    tables = root.findChildren(DataTable)
    assert {table.export_base for table in tables} == {
        "TEagle_family", "TEagle_non_TE_repeats", "TEagle_unclassified_repeats"}
    assert all(table.rowCount() == 1 and table.export_notes for table in tables)
    captured = []
    monkeypatch.setattr(main.widgets, "pick_table_format", lambda *a, **k: "tsv")
    monkeypatch.setattr(main.widgets, "export_table",
                        lambda h, r, b, p, fmt=None, notes=None: captured.append((r, b, notes)))
    for button in root.findChildren(QPushButton):
        if button.text().upper() == "EXPORT TABLE":
            button.click()
    assert {(base, tuple(tuple(row) for row in rows), tuple(notes)) for rows, base, notes in captured} == {
        (table.export_base, tuple(tuple(row) for row in table.rows_data()), tuple(table.export_notes))
        for table in tables}


def test_orf_and_splice_tables_have_visible_exact_export_surfaces(win):
    import engine
    from widgets import DataTable
    from PySide6.QtWidgets import QPushButton

    seq = ">orf\nATG" + "AAA" * 45 + "TAA"
    _accept_analysis(win, engine.run_analyze({"sequence": seq}), seq)
    orf_table = next(table for table in win.card_struct.findChildren(DataTable)
                     if table.export_base == "TEagle_ORFs")
    assert orf_table.rowCount() >= 1
    assert any(button.text().upper() == "EXPORT ORFS (DNA FASTA)"
               for button in win.card_struct.findChildren(QPushButton))

    win.state["splice_seq"] = "AAAGT" + "C" * 15 + "AGTTT"
    win.state["last_rec"] = {"id": "splice_record"}
    win._render_splice({"counts": {"exons": 2, "introns": 1}, "canonical_introns": 1, "strand": "+",
                        "exons": [{"start": 0, "end": 3}, {"start": 22, "end": 25}],
                        "introns": [{"start": 3, "end": 22, "donor": "GT", "acceptor": "AG", "canonical": True}]})
    splice_root = win.spliceBody.parentWidget()
    introns = next(table for table in splice_root.findChildren(DataTable)
                   if table.export_base == "TEagle_splice_introns")
    assert introns.rowCount() == 1 and introns.export_notes
    labels = {button.text().upper() for button in splice_root.findChildren(QPushButton)}
    assert {"EXPORT TABLE", "EXPORT SPLICE EXONS (FASTA)", "EXPORT SPLICE MODEL ▾"}.issubset(labels)


@pytest.mark.parametrize("extension,marker", [
    ("gff3", "##gff-version 3"), ("bed", "\texon_1|"), ("fasta", ".exon1")])
def test_native_gene_model_export_routes_shared_coordinate_writer(win, monkeypatch, tmp_path,
                                                                   extension, marker):
    import main

    out = tmp_path / f"model.{extension}"
    monkeypatch.setattr(main.QFileDialog, "getSaveFileName", staticmethod(lambda *a, **k: (str(out), "")))
    model = {"exons": [{"start": 0, "end": 6, "strand": "+"}],
             "introns": [{"start": 6, "end": 9, "strand": "+", "donor": "GT", "acceptor": "AG"}],
             "cds": [{"start": 0, "end": 6, "strand": "+", "phase": 0}]}
    win._export_gene_model(model, "ACGTACGTA", "safe", "fetched_gene_model", extension)
    assert marker in out.read_text(encoding="utf-8")


def test_native_gene_model_missing_cds_phase_refuses_only_gff3_and_preserves_destination(
        win, monkeypatch, tmp_path):
    import main

    model = {"exons": [{"start": 0, "end": 6, "strand": "+"}],
             "introns": [], "cds": [{"start": 0, "end": 6, "strand": "+"}]}
    gff = tmp_path / "model.gff3"
    gff.write_text("preserve", encoding="utf-8")
    monkeypatch.setattr(main.QFileDialog, "getSaveFileName",
                        staticmethod(lambda *a, **k: (str(gff), "")))
    win._export_gene_model(model, "ACGTAC", "safe", "fetched_gene_model", "gff3")
    assert gff.read_text(encoding="utf-8") == "preserve"

    for extension in ("bed", "fasta"):
        out = tmp_path / f"model.{extension}"
        monkeypatch.setattr(main.QFileDialog, "getSaveFileName",
                            staticmethod(lambda *a, _out=out, **k: (str(_out), "")))
        win._export_gene_model(model, "ACGTAC", "safe", "fetched_gene_model", extension)
        assert out.exists() and out.stat().st_size > 0


def test_native_gene_model_validation_failure_preserves_existing_file(win, monkeypatch, tmp_path):
    import main

    out = tmp_path / "model.gff3"
    out.write_text("validated existing model\n", encoding="utf-8")
    monkeypatch.setattr(main.QFileDialog, "getSaveFileName", staticmethod(lambda *a, **k: (str(out), "")))
    messages = []
    monkeypatch.setattr(win, "_banner", lambda message, level="error": messages.append((message, level)))
    win._export_gene_model({"exons": [{"start": 0, "end": 99}], "introns": [], "cds": []},
                           "ACGT", "safe", "fetched_gene_model", "gff3")
    assert out.read_text(encoding="utf-8") == "validated existing model\n"
    assert not (tmp_path / "model.gff3.part").exists()
    assert messages and messages[-1][1] == "warn"


def test_fetched_gene_model_has_visible_shared_export_hook(win):
    import engine
    from PySide6.QtWidgets import QPushButton

    sequence = ">record\nATG" + "AAA" * 45 + "TAA"
    result = engine.run_analyze({"sequence": sequence})
    win.state["features"] = {"exons": [{"start": 0, "end": 30}], "introns": [],
                             "cds": [{"start": 0, "end": 30}]}
    _accept_analysis(win, result, sequence)
    assert any(button.text().upper() == "EXPORT GENE MODEL ▾"
               for button in win.card_struct.findChildren(QPushButton))


def test_secondary_structure_detail_has_visible_table_export(win, monkeypatch):
    from PySide6.QtWidgets import QApplication, QDialog, QPushButton
    from widgets import DataTable

    metric = {"p3": -2.5, "vrna": -2.0, "flag": "caution", "agree": "agree"}
    candidate = {
        "id": "P1", "left_seq": "ACGTACGTACGT", "right_seq": "TGCATGCATGCA",
        "qc": {"ok": True, "left": {"hairpin": metric, "self_dimer": metric},
               "right": {"hairpin": metric, "self_dimer": metric},
               "hetero_dimer": metric, "end_stability": metric,
               "conditions": {"temp_c": 37, "mv_conc": 50, "dv_conc": 1.5, "dna_conc": 250},
               "engines": {"primer3": "x", "viennarna": "y"}},
    }
    monkeypatch.setattr(QDialog, "exec", lambda self: 0)
    win._structure_detail(candidate)
    dialog = next(widget for widget in QApplication.topLevelWidgets()
                  if isinstance(widget, QDialog) and widget.windowTitle().startswith("Secondary structure"))
    table = dialog.findChild(DataTable)
    assert table.export_base == "TEagle_P1_secondary_structure"
    assert table.rowCount() == 6 and table.export_notes
    assert any(button.text().upper() == "EXPORT TABLE" for button in dialog.findChildren(QPushButton))
    dialog.close()


# ---------- specimen-identity / stale-state provenance integrity (Loop-1 fixes) ----------
def test_genome_scan_submits_local(win):
    win._prepared_genomes = [{"accession": "GCF_000001405.40", "n_seqs": 24, "bytes": 1_000_000}]
    win._refresh_genome_dropdown()                            # dropdown lists ONLY downloaded genomes now
    win.genomeOrg.setCurrentIndex(win.genomeOrg.findData("Homo sapiens"))
    ops = []
    win.engine.submit = lambda op, body=None, key=None: ops.append((op, body, key))
    cand = {"left_seq": "TTCACCACCATGGAGAAGGC", "right_seq": "GGCATGGACTGTGGTCATGAG", "id": "P1"}
    win._scan_genome(cand)
    assert ops[0][0] == "genome_pcr" and ops[0][2] == "genome_pcr"
    assert ops[0][1]["organism"] == "Homo sapiens" and ops[0][1]["params"]["min_perfect"] == 15
    # cand + the org that started the scan are captured, so a later dropdown change can't retarget it
    assert win._pending_scan == {"cand": cand, "org": "Homo sapiens"} and win._genome_inflight is True


def test_genome_scan_requires_organism(win):
    win.genomeOrg.setCurrentIndex(0)                          # the "— select organism —" placeholder (data=None)
    ops = []
    win.engine.submit = lambda op, body=None, key=None: ops.append(op)
    win._scan_genome({"left_seq": "ACGT" * 5, "right_seq": "TTTT" * 5, "id": "P1"})
    assert ops == [] and win._genome_inflight is False        # no wrong-species scan runs without an explicit pick


def test_re_render_struct_card_no_duplicate_headers(win):
    # regression: a widget-only clear left addLayout'd section headers (STRUCTURAL EVIDENCE / PROTEIN DOMAINS)
    # orphaned but visible, so they DUPLICATED on a second Run analysis. _clear_layout must remove them.
    from PySide6.QtWidgets import QLabel
    from helpers import fixture_seq
    import engine
    seq = fixture_seq("M11240"); win.state["seq"] = seq; win.state["analyzed_seq"] = seq
    def n_headers():
        return sum(1 for w in win.card_struct.findChildren(QLabel) if w.text() == "STRUCTURAL EVIDENCE")
    _accept_analysis(win, engine.run_analyze({"sequence": seq}), seq); a = n_headers()
    _accept_analysis(win, engine.run_analyze({"sequence": seq}), seq); b = n_headers()
    assert a == 1 and b == 1                                  # exactly one header, not accumulating across runs


def test_genome_scan_renders_in_its_own_panel(win):
    # the whole-genome scan is now a dedicated card (06), not the in-silico PCR card
    assert hasattr(win, "card_genome") and hasattr(win, "genomeBody")
    win._genome_inflight = True
    win._on_genome_pcr({"ok": True, "organism": "x", "assemblyName": "y", "assemblyAccession": "GCF_9.9", "n_seqs": 1,
        "amplicons": [{"source": "NC_1", "start": 10, "end": 210, "length": 201, "strand": "+",
                       "single_primer": False, "on_target": True}],
        "summary": {"n_on": 1, "n_off": 0, "n_pair": 1, "n_single": 0, "has_locus": True,
                    "verdict": "copy-specific — 1 on-target", "per_source": [("NC_1", 1)],
                    "size_mode": 201, "size_mode_n": 1, "size_min": 201, "size_max": 201},
        "provenance": {}})
    assert win.genomeBody.count() > 1                          # gel + verdict + table + note rendered in the genome panel


def test_on_genome_pcr_renders_sealed_advisory(win):
    from PySide6.QtWidgets import QLabel
    from widgets import DataTable

    win._genome_inflight = True
    win._on_genome_pcr({"ok": True, "organism": "Saccharomyces cerevisiae", "assemblyName": "R64",
        "assemblyAccession": "GCF_000146045.2", "n_seqs": 17,
        "advisory_note": "candidate priming sites — not wet-lab-validated amplicons",
        "amplicons": [{"source": "NC_001133.9", "start": 1001, "end": 1240, "length": 240, "strand": "+",
                       "single_primer": False, "on_target": False, "pair": "pair"}],
        "summary": {"n_pair": 1, "n_on": 0, "n_off": 1, "n_single": 0, "verdict": "1 genomic priming site", "has_locus": False},
        "provenance": {"manifestSha256": "de" * 32,
                       "databases": [{"name": "RefSeq assembly R64", "version": "GCF_000146045.2", "sha256": "abc123"}]}})
    lp = win.state["lastGenomeScan"]                          # genome scan now renders in its OWN panel, not the PCR card
    assert len(lp["amplicons"]) == 1 and lp["lanes"][0]["advisory"] is True
    assert win._genome_inflight is False
    visible = " ".join(label.text() for label in win.genomeBody.parentWidget().findChildren(QLabel))
    whole_card = " ".join(label.text() for label in win.card_genome.findChildren(QLabel))
    scan_table = next(table for table in win.genomeBody.parentWidget().findChildren(DataTable)
                      if table.export_base == "TEagle_offtarget_scan")
    exported = " ".join(scan_table.export_notes)
    assert "genomic priming site" in visible and "genomic priming site" in exported
    assert "on-target first" not in visible and "on-target first" not in exported
    assert "products are listed on-target first" not in whole_card.lower()
    assert "with a design locus" in whole_card.lower()


def test_genome_prepare_resumes_pending_scan(win):
    ops = []
    win.engine.submit = lambda op, body=None, key=None: ops.append((op, body, key))
    cand = {"left_seq": "ACGT" * 5, "right_seq": "TTTT" * 5, "id": "P1"}
    win.genomeOrg.setCurrentIndex(win.genomeOrg.findData("Saccharomyces cerevisiae"))
    win._prepare_genome("Saccharomyces cerevisiae", then_scan={"cand": cand, "org": "Saccharomyces cerevisiae"})
    assert ops[0][0] == "genome_prepare" and win._genome_prep_inflight is True
    ops.clear()
    win._on_genome_prepare({"ok": True, "organism": "Saccharomyces cerevisiae",
                            "assemblyAccession": "GCF_000146045.2", "n_seqs": 17, "bytes": 12_000_000})
    assert win._genome_prep_inflight is False and ops and ops[0][0] == "genome_pcr"   # resumed into the SAME scan (captured org)


def test_genome_prepare_cache_failure_preserves_remove_and_download_action(win):
    messages = []
    win._render_genome_status = lambda msg, level="error": messages.append((msg, level))
    win._banner = lambda msg, level="error": messages.append((msg, level))
    win._on_genome_prepare({"ok": False, "cache_corrupt": True,
                            "error": "Open Manage genomes, remove this cached genome, then download it again."})
    rendered = " ".join(msg for msg, _ in messages)
    assert "Manage genomes" in rendered and "remove" in rendered and "download" in rendered
    assert all(level == "warn" for _, level in messages)


def test_banner_level_distinguishes_success_from_error(win):
    # messages surface as a single, closable notification dialog (not a top-of-panel banner); level drives the
    # colour style (theme.py #notifmsg[level]) and a new message replaces the previous one (never stacks).
    from PySide6.QtWidgets import QLabel
    win._banner("Genome ready", level="success")
    dlg = win._notif
    assert dlg is not None and dlg.property("level") == "success"
    msg = dlg.findChild(QLabel, "notifmsg")
    assert msg is not None and msg.text() == "Genome ready"     # message shown verbatim (no glyph)
    win._banner("bad input")                                   # default level = error — replaces the success dialog
    assert win._notif is not None and win._notif.property("level") == "error"
    assert win._notif is not dlg                               # single active instance: replaced, not stacked
    win._clear_banner()
    assert win._notif is None


def test_busybar_is_indeterminate_and_updatable():
    from widgets import BusyBar
    b = BusyBar("working…")
    assert b.bar.minimum() == 0 and b.bar.maximum() == 0       # 0..0 = indeterminate 'busy' animation
    b.set_text("stage 2")
    assert b.caption.text() == "stage 2"


def test_genome_dropdown_lists_only_prepared(win):
    # the PCR organism dropdown must show ONLY downloaded+verified genomes, never the full catalog
    win._on_genome_list({"ok": True, "genomes": []})
    assert [win.genomeOrg.itemData(i) for i in range(win.genomeOrg.count())] == [None]   # placeholder only
    win._on_genome_list({"ok": True, "genomes": [{"accession": "GCF_000146045.2", "n_seqs": 17, "bytes": 12_000_000}]})
    orgs = [win.genomeOrg.itemData(i) for i in range(win.genomeOrg.count())]
    assert orgs == [None, "Saccharomyces cerevisiae"]        # placeholder + the one downloaded organism


def test_manager_survives_rebuild_while_open(win):
    # a rebuild while the manager is open (rescale, add-organism, download/delete completing, any genome_list refresh)
    # must NOT collapse the fixed-size dialog to a sliver: the in-place layout swap once left re-added widgets measured
    # as hidden so sizeHint()->~20px and setFixedSize pinned it there unrecoverably.
    win._genome_mgr_open = True
    gen = {"ok": True, "genomes": [{"accession": "GCF_000146045.2", "n_seqs": 17, "bytes": 12_000_000}]}
    win._on_genome_list(gen)                                  # build 1
    h1 = win._genome_mgr.height()
    win._on_genome_list(gen)                                  # build 2 — the rebuild-while-open path
    win._on_genome_list(gen)                                  # build 3 — a second refresh
    h3 = win._genome_mgr.height()
    assert h1 > 300 and h3 > 300 and abs(h3 - h1) <= 8        # stays full height, never collapses to a ~20px sliver


def test_manager_download_button_not_clipped(win):
    # Every Action button must sit FULLY INSIDE its cell and CENTRED in it. Asserted in viewport coordinates against
    # visualRect, because the failures were positional: QTableWidget::item padding is applied to cell widgets too, so
    # an unbudgeted row/column offset the button by the padding (+9/+6 at scale 1.0) and sheared its bottom edge.
    # Checking the holder's SIZE alone missed this for two releases — check where the button actually lands.
    from PySide6.QtWidgets import QPushButton
    from PySide6.QtCore import QRect
    win._genome_mgr_open = True
    win._on_genome_list({"ok": True, "genomes": []})
    tbl = win._genome_mgr_table
    act = tbl.columnCount() - 1
    widths = set()
    for r in range(tbl.rowCount()):
        btn = tbl.cellWidget(r, act).findChild(QPushButton)
        assert btn is not None
        cell = tbl.visualRect(tbl.model().index(r, act))
        brect = QRect(btn.mapTo(tbl.viewport(), btn.rect().topLeft()), btn.size())
        assert cell.contains(brect), f"row {r}: button {brect} not contained by cell {cell}"
        assert abs(brect.center().x() - cell.center().x()) <= 2, f"row {r}: off-centre horizontally"
        assert abs(brect.center().y() - cell.center().y()) <= 2, f"row {r}: off-centre vertically"
        assert btn.height() >= 18                            # still a comfortable, clickable target
        widths.add(btn.width())
    assert len(widths) == 1, f"buttons must share ONE width so they align down the column, got {sorted(widths)}"


def test_on_genome_list_refreshes_dropdown_even_when_manager_closed(win):
    # the split handler must update the dropdown on EVERY genome_list, not only when the manager dialog is open
    win._genome_mgr = None; win._genome_mgr_open = False
    win._on_genome_list({"ok": True, "genomes": [{"accession": "GCF_000001405.40", "n_seqs": 24, "bytes": 1_000_000}]})
    assert win._genome_mgr is None                            # no dialog resurrected
    assert win.genomeOrg.findData("Homo sapiens") >= 0        # ...but the dropdown still refreshed


def test_genome_list_failure_keeps_cached_dropdown(win):
    # a transient genome_list failure must NOT blank the dropdown while genomes are still cached on disk
    win._on_genome_list({"ok": True, "genomes": [{"accession": "GCF_000146045.2", "n_seqs": 17, "bytes": 12_000_000}]})
    assert win.genomeOrg.findData("Saccharomyces cerevisiae") >= 0
    win._on_genome_list({"ok": False, "error": "WSL backend hiccup"})    # transient blip
    assert win._prepared_genomes and win.genomeOrg.findData("Saccharomyces cerevisiae") >= 0   # last-good set preserved


def test_scan_completion_refreshes_open_manager(win):
    # a manager opened mid-scan disables its row buttons; scan completion must refresh it so they re-enable
    win._genome_mgr_open = True
    win._on_genome_list({"ok": True, "genomes": []})          # open the (non-modal) manager dialog
    assert win._genome_mgr is not None and win._genome_mgr.isVisible()
    ops = []
    win.engine.submit = lambda op, body=None, key=None: ops.append(op)
    win._genome_inflight = True
    win._on_genome_pcr({"ok": True, "organism": "x", "assemblyName": "y", "assemblyAccession": "GCF_9.9",
                        "n_seqs": 1, "amplicons": [], "summary": {}, "provenance": {}})
    assert "genome_list" in ops                                # scan completion re-lists -> manager rebuilt, buttons re-enabled
    win._genome_mgr.close()


def test_genome_manager_stale_refresh_does_not_reopen(win):
    # a genome_list result landing AFTER the manager was closed must not resurrect the dialog (M1 completeness)
    win._genome_mgr = None
    win._genome_mgr_open = False
    win._on_genome_list({"ok": True, "genomes": []})
    assert win._genome_mgr is None                            # stale refresh -> no-op
    win._genome_mgr_open = True                               # an explicit open does build one
    win._on_genome_list({"ok": True, "genomes": []})
    assert win._genome_mgr is not None
    win._genome_mgr.close()


def test_genome_need_prepare_uses_captured_org(win, monkeypatch):
    # start a scan for organism X, change the dropdown mid-scan; the download must target X, not the new pick
    win._prepared_genomes = [{"accession": "GCF_000146045.2", "n_seqs": 17, "bytes": 12_000_000}]
    win._refresh_genome_dropdown()                            # dropdown lists ONLY downloaded genomes now
    win.genomeOrg.setCurrentIndex(win.genomeOrg.findData("Saccharomyces cerevisiae"))
    ops = []
    win.engine.submit = lambda op, body=None, key=None: ops.append((op, body))
    win._scan_genome({"left_seq": "ACGT" * 5, "right_seq": "TTTT" * 5, "id": "P1"})   # captures org = yeast
    win.genomeOrg.setCurrentIndex(0)                          # user flips the dropdown to the placeholder (data=None)
    from PySide6.QtWidgets import QMessageBox
    monkeypatch.setattr(QMessageBox, "exec", lambda self: QMessageBox.Yes)
    ops.clear()
    win._on_genome_pcr({"ok": False, "need_prepare": True})
    assert ops and ops[0][0] == "genome_prepare"
    assert ops[0][1]["organism"] == "Saccharomyces cerevisiae"   # the captured org, not the live placeholder


def test_corrupt_genome_cache_never_enters_download_prepare_loop(win, monkeypatch):
    from PySide6.QtWidgets import QMessageBox

    win._pending_scan = {"org": "Saccharomyces cerevisiae", "cand": {"id": "P1"}}
    ops = []
    messages = []
    win.engine.submit = lambda op, body=None, key=None: ops.append((op, body))
    win._render_genome_status = lambda msg, level="error": messages.append((msg, level))
    win._banner = lambda msg, level="error": messages.append((msg, level))
    monkeypatch.setattr(QMessageBox, "exec", lambda self: pytest.fail("corrupt cache opened download prompt"))
    win._on_genome_pcr({"ok": False, "cache_corrupt": True,
                        "error": "cached genome failed its integrity check"})
    assert not any(op == "genome_prepare" for op, _ in ops)
    assert win._pending_scan is None
    visible = " ".join(msg for msg, _ in messages)
    assert "Manage genomes" in visible and "remove" in visible and "download" in visible


def test_corrupt_annotation_cache_preserves_remove_and_download_action(win):
    messages = []
    win._refresh_genome_manager = lambda: None
    win._banner = lambda msg, level="error": messages.append((msg, level))
    win._on_genome_annotate({"ok": False, "cache_corrupt": True,
                             "error": "Open Manage genomes, remove this cached genome, then download it again."})
    rendered = " ".join(msg for msg, _ in messages)
    assert "Manage genomes" in rendered and "remove" in rendered and "download" in rendered
    assert messages[-1][1] == "warn"


def test_flank_region_clickable_and_slices_analyzed_seq(win):
    # a flank region from the gene-model viewer offers copy + slices the analyzed snapshot; contextual actions
    # (primer/splice) appear only when the region is long enough to design in (item 14 — no primer on a 10 bp motif)
    win.state["analyzed_seq"] = ">x\n" + "ACGTACGTAC" + "G" * 90          # first 10 bases ACGTACGTAC, 100 total
    short = [lbl for lbl, _ in win._region_menu({"start": 0, "end": 10, "label": "5prime_flank", "strand": "+"})]
    assert any("Copy FASTA" in l for l in short)
    assert not any("primer" in l.lower() for l in short)                 # 10 bp is too short to design a primer in
    long = [lbl for lbl, _ in win._region_menu({"start": 0, "end": 80, "label": "flank", "strand": "+"})]
    assert any("primer" in l.lower() for l in long)                      # an 80 bp region can
    assert win._slice(0, 10) == "ACGTACGTAC"                             # flank coords index analyzed_seq


def test_app_theme_propagates_to_genome_viewers(win):
    import widgets
    win._load_sample(); win._run_analysis()
    assert _spin(lambda: win.state.get("analyzed_clean"))
    gvs = win.findChildren(widgets.GenomePanel)
    assert gvs and win.theme in ("dark", "light")
    start_theme, start_vtheme = win.theme, gvs[0].theme
    gvs[0].view = {"start": 100.0, "end": 500.0}; saved = dict(gvs[0].view)
    win._toggle_theme()                                        # toggle the app theme
    assert win.theme != start_theme and gvs[0].theme != start_vtheme   # viewer follows the app theme
    assert gvs[0].view == saved                                # pan/zoom preserved (in-place re-render)
    after_one = win.theme
    gvs[0]._set_theme("dark", user=True)                       # per-viewer MANUAL override (button click)
    win._toggle_theme(); win._toggle_theme()                   # two app toggles = identity; must NOT stomp a manual pick (H1)
    assert win.theme == after_one and gvs[0].theme == "dark"   # user's viewer choice is kept across app-theme changes


def test_feature_slice_uses_analyzed_snapshot_not_live_box(win):
    # #1 (Loop-3): after an unanalyzed edit, a panel-01 feature copy must slice the ANALYZED sequence,
    # not the live box (which now backs splice/annotate) — else copies return wrong bases under a header
    win.state["analyzed_seq"] = ">A\n" + "AAAACCCCGGGGTTTT"
    win.seq.setPlainText(">B\n" + "TTTTTTTTTTTTTTTT")         # user edits the box after analysis
    assert win._slice(0, 8) == "AAAACCCC"                     # from analyzed A, not the all-T live box


def test_second_primer_design_blocked_while_inflight(win):
    # #4 (Loop-3): concurrent designs would clobber the shared _design_tmpl -> block the second
    win.state["seq"] = ">x\n" + "ACGT" * 60
    win.state["analyzed_clean"] = win._clean_seq(win.state["seq"])
    ops = []
    win.engine.submit = lambda op, body=None, key=None: ops.append(op)
    win._design()
    win._design()                                             # second call while the first is in flight
    assert ops == ["primers"] and win._design_inflight is True


def test_set_seq_clears_prior_gene_model(win):
    # a fetched gene model must not carry over to an unrelated uploaded/sample specimen
    win.state["features"] = {"counts": {"exons": 3, "introns": 2, "cds": 1}}
    win._load_sample()
    assert win.state["features"] is None


def test_pcr_omits_all_design_geometry_for_mismatched_template(win):
    # A pair designed on another template must not assert endpoint identity or target coverage.
    seq = ">x\n" + "ACGT" * 60
    win.state["seq"] = seq
    win.state["analyzed_clean"] = win._clean_seq(seq)          # pass the stale-gate
    win.state["pcrPairs"] = [{"left_seq": "ACGTACGTACGTACGTAC", "right_seq": "TTTTGGGGTTTTGGGGTT",
                              "left_pos": [0, 18], "right_pos": [30, 48], "_tmpl_sig": "SOME-OTHER-TEMPLATE"}]
    bodies = []
    win.engine.submit = lambda op, body=None, key=None: bodies.append(body)
    win._run_pcr()
    assert bodies
    for field in ("primer3_target", "included_region", "intended_product_span", "design_provenance"):
        assert bodies[0][field] is None
    assert bodies[0]["design_assessment"] == "not_run"  # fixture has no sealed candidate ancestry


def test_pcr_marks_sealed_candidate_geometry_withheld_on_template_mismatch(win):
    seq = ">x\n" + "ACGT" * 60
    win.state["seq"] = seq
    win.state["analyzed_clean"] = win._clean_seq(seq)
    win.state["pcrPairs"] = [{"left_seq": "ACGTACGTACGTACGTAC", "right_seq": "TTTTGGGGTTTTGGGGTT",
                              "left_pos": [0, 18], "right_pos": [30, 48],
                              "_tmpl_sig": "OTHER", "design_provenance": {"manifestSha256": "x"},
                              "primer3_target": [10, 5], "included_region": [0, 40],
                              "intended_product_span": [0, 48], "_design_sequence": ">design\n" + "TGCA" * 60}]
    bodies = []
    win.engine.submit = lambda op, body=None, key=None: bodies.append(body)
    win._run_pcr()
    assert bodies[0]["design_assessment"] == "withheld_template_mismatch"
    assert bodies[0]["design_provenance"] == {"manifestSha256": "x"}
    assert bodies[0]["intended_product_span"] == [0, 48]
    assert bodies[0]["design_sequence"].startswith(">design")


def test_domain_design_routes_domain_as_primer3_target(win):
    seq = "ACGT" * 100
    win.state["seq"] = seq
    win.state["analyzed_clean"] = win._clean_seq(seq)
    submitted = []
    win.engine.submit = lambda op, body=None, key=None: submitted.append((op, body, key))
    win._design_for_domain(80, 140, "reverse-transcriptase domain")
    assert submitted and submitted[0][0] == "primers"
    body = submitted[0][1]
    assert body["primer3_target"] == [80, 60]
    assert body["included_region"] == [0, len(seq)]


def test_full_explicit_subsequence_routes_to_general_design_not_whole_template_target(win):
    seq = "ACGT" * 100
    submitted = []
    win.engine.submit = lambda op, body=None, key=None: submitted.append((op, body, key))
    win._design_for_domain(0, len(seq), "selected amplicon", seq=seq)
    body = submitted[0][1]
    assert body["primer3_target"] is None and body["included_region"] is None


def test_invalid_primer_gui_value_fails_closed_and_clears_prior_candidates(win):
    seq = "ACGT" * 100
    win.state["seq"] = seq; win.seq.setPlainText(seq)
    win.state["analyzed_clean"] = win._clean_seq(seq)
    win.state["candidates"] = [{"id": "old"}]; win.state["lastPrimers"] = {"candidates": [{"id": "old"}]}
    win.pcrStageAll.setEnabled(True)
    win.pfields["pMin"].setText("not-a-number")
    submitted = []
    win.engine.submit = lambda *a, **k: submitted.append((a, k))
    win._design()
    assert not submitted and win.state["candidates"] == [] and win.state["lastPrimers"] is None
    assert not win.pcrStageAll.isEnabled()


def test_failed_primer_request_clears_candidates_and_stage_all(win):
    win.state["candidates"] = [{"id": "old"}]; win.state["lastPrimers"] = {"candidates": [{"id": "old"}]}
    win.pcrStageAll.setEnabled(True)
    win._on_user_error("primers", "invalid primer parameters")
    assert win.state["candidates"] == [] and win.state["lastPrimers"] is None
    assert not win.pcrStageAll.isEnabled()


def test_pcr_queue_keeps_same_oligos_with_distinct_sealed_intended_spans(win):
    common = {"left_seq": "A" * 18, "right_seq": "T" * 18, "product_size": 100}
    first = {**common, "intended_product_span": [10, 110],
             "design_provenance": {"manifestSha256": "1" * 64}}
    second = {**common, "intended_product_span": [20, 120],
              "design_provenance": {"manifestSha256": "2" * 64}}
    win._add_pcr_pair(first); win._add_pcr_pair(second)
    assert len(win.state["pcrPairs"]) == 2


@pytest.mark.parametrize("updates", [{"pcrTP": "2.9"}, {"pcrTP": "19"},
                                      {"pcrMM": "-1"}, {"pcrMM": "19"},
                                      {"pcrPmin": "500", "pcrPmax": "100"}])
def test_invalid_pcr_gui_values_do_not_submit(win, updates):
    seq = ">x\n" + "ACGT" * 80
    win.state["seq"] = seq; win.seq.setPlainText(seq)
    win.state["analyzed_clean"] = win._clean_seq(seq)
    win.state["pcrPairs"] = [{"left_seq": "A" * 18, "right_seq": "T" * 18}]
    for field, value in updates.items():
        win.pfields[field].setText(value)
    submitted = []
    win.engine.submit = lambda *a, **k: submitted.append((a, k))
    win._run_pcr()
    assert not submitted


def test_partial_and_failed_pcr_wording_cannot_claim_exhaustive_absence(win):
    from PySide6.QtWidgets import QLabel
    from widgets import DataTable
    amp = {"pair": "P1", "source": "s", "start": 10, "end": 210, "length": 200,
           "product_length": 200, "product_key": "s:10:210:200", "mechanism_key": "s:F+R:10:210:200",
           "fwd_mm": 0, "rev_mm": 0, "single_primer": False,
           "matches_intended_product": None, "covers_primer3_target": None, "on_target": None,
           "seq": "A" * 200}
    lane = {"label": "P1", "amplicons": [amp], "has_locus": False, "status": "partial_capped"}
    win._render_pcr([lane], [amp], [{"pair": "P1", "status": "partial_capped", "error": None}])
    text = " ".join(l.text() for l in win.pcrBody.parentWidget().findChildren(QLabel)).lower()
    table_text = " ".join(t.item(r, c).text() for t in win.pcrBody.parentWidget().findChildren(DataTable)
                          for r in range(t.rowCount()) for c in range(t.columnCount()) if t.item(r, c)).lower()
    assert "lower bound" in text
    for forbidden in ("no additional products", "no off-target", "exhaustive search complete"):
        assert forbidden not in text and forbidden not in table_text
    win._render_pcr([{"label": "P1", "amplicons": [], "has_locus": False,
                      "status": "partial_capped"}], [],
                    [{"pair": "P1", "status": "partial_capped", "error": None}])
    capped_empty = " ".join(l.text() for l in win.pcrBody.parentWidget().findChildren(QLabel)).lower()
    assert "capped search stopped" in capped_empty and "completed search recorded zero" not in capped_empty
    win._render_pcr([{"label": "P1", "amplicons": [], "has_locus": False, "status": "failed"}], [],
                    [{"pair": "P1", "status": "failed", "error": "failure"}])
    failed_text = " ".join(l.text() for l in win.pcrBody.parentWidget().findChildren(QLabel)).lower()
    assert "failed and was not assessed" in failed_text
    assert "no additional products" not in failed_text and "no off-target" not in failed_text


def test_local_amplicon_wording_distinguishes_unassessed_and_additional():
    import main
    assert main._amp_call({"single_primer": False, "matches_intended_product": True,
                           "covers_primer3_target": None, "on_target": None}) == "intended product; target coverage unavailable"
    assert main._amp_call({"single_primer": False, "matches_intended_product": True,
                           "covers_primer3_target": False, "on_target": False}) == "intended product; target not covered"
    assert main._amp_call({"single_primer": False, "matches_intended_product": False,
                           "covers_primer3_target": True, "on_target": False}) == "additional; covers target"
    assert main._amp_call({"single_primer": False, "matches_intended_product": False,
                           "covers_primer3_target": False, "on_target": False}) == "additional product"


def test_label_ink_accepts_3char_hex():
    import figures
    assert figures._label_ink("#888") in ("#fff", "#111")      # the #888 fallback must not raise
    assert figures._label_ink("#888") == figures._label_ink("#888888")


def test_design_blocked_after_box_edit(win):
    # state['seq'] now tracks the live box; the stale-gate must still block design against an unanalyzed edit
    win._load_sample(); win._run_analysis()
    assert _spin(lambda: win.state.get("analyzed_clean"))
    ops = []; win.engine.submit = lambda op, body=None, key=None: ops.append(op)
    win.seq.setPlainText("ACGTACGTACGTACGTACGTACGTACGTACGT")   # edit after analysis
    win._design()
    assert "primers" not in ops                                # blocked, not designed against the unanalyzed sequence


def test_user_edit_drops_fetched_source(win):
    # after a fetch, editing the box means the specimen is no longer that accession -> source must clear
    accession, sequence = "M11240.1", "ACGT" * 30
    win._on_fetch({"ok": True, "accession": accession, "organism": "Drosophila melanogaster",
                   "title": "copia", "length": len(sequence), "fasta": f">{accession}\n{sequence}",
                   "taxid": "7227", "source": "NCBI nuccore",
                   "endpoint": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
                   "retrievedUtc": "2026-08-21T00:00:00+00:00",
                   "sourceUrl": f"https://www.ncbi.nlm.nih.gov/nuccore/{accession}"})
    assert win.state["source"] and win.state["source"].get("accession") == accession
    win.seq.setPlainText("ACGTACGTACGTACGTACGTACGT")           # genuine user edit
    assert win.state["source"] is None and win.accMeta.text() == ""


def test_programmatic_load_keeps_source(win):
    # the _loading guard: a fetch's own setPlainText must NOT null the source it is about to set
    accession, sequence = "X05424.1", "ACGT" * 25
    win._on_fetch({"ok": True, "accession": accession, "organism": "Nicotiana", "title": "t",
                   "length": len(sequence), "fasta": f">{accession}\n{sequence}", "taxid": "4097",
                   "source": "NCBI nuccore",
                   "endpoint": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
                   "retrievedUtc": "2026-08-21T00:00:00+00:00",
                   "sourceUrl": f"https://www.ncbi.nlm.nih.gov/nuccore/{accession}"})
    assert win.state["source"] and win.state["source"].get("accession") == accession


def test_seq_tracks_box_not_stale_after_refetch(win):
    # analyze A, then fetch B without re-analyzing: state['seq'] (what splice/annotate submit) must be B
    win._set_seq(">A\n" + "AAAA" * 20); win.state["seq"] = win.seq.toPlainText().strip()
    win._run_analysis(); _spin(lambda: win.state.get("analyzed_clean"))
    win._on_fetch({"ok": True, "accession": "B", "organism": "o", "title": "t",
                   "length": 80, "fasta": ">B\n" + "CCCC" * 20})
    assert "CCCC" in win.state["seq"] and "AAAA" not in win.state["seq"]   # current box, not the stale analyzed A


def test_pcr_batch_renders_despite_a_failing_pair(win):
    # a batch where one sub-job errors must still render the gel for the good lanes (slot filled, not stalled)
    win._pcr_gen = 5
    win._pcr_run = {"gen": 6, "results": [None, None], "single": False}
    win._pcr_gen = 6
    win.state["lastPcr"] = None
    win._pcr_slot("pcr#6#0", {"amplicons": [{"start": 1, "end": 100, "length": 99, "fwd_mm": 0,
                              "rev_mm": 0, "on_target": True, "source": "s", "seq": "ACGT"}]})
    win._pcr_slot("pcr#6#1", {"error": "boom", "amplicons": []})              # failed pair
    assert win.state["lastPcr"] and len(win.state["lastPcr"]["lanes"]) == 2   # gel rendered, both lanes present
    assert [s["status"] for s in win.state["lastPcr"]["statuses"]] == ["complete", "failed"]
    assert win.runPcrBtn.isEnabled()                                          # button re-enabled after batch


def test_pcr_stale_batch_result_dropped(win):
    # a late result from a superseded batch (older gen) must not touch the current batch
    win._pcr_run = {"gen": 9, "results": [None], "single": True}
    win._pcr_slot("pcr#8#0", {"amplicons": []})                # stale gen 8 -> ignored
    assert win._pcr_run["results"][0] is None


# ---------- WSL / optional-backend degradation (must never disable core science) ----------
def test_wsl_absent_keeps_optional_disabled(win):
    win._on_wsl_status({"wsl2": False})
    assert not win.annotateBtn.isEnabled() and not win.spliceBtn.isEnabled()
    assert "not installed" in win.wslStatus.text()
    assert not win.wslInstallBtn.isHidden()                    # installer stays reachable — it guides WSL setup


def test_genome_scan_table_declares_ispcr_1based(win):
    """The genome-scan table shows isPcr coordinates verbatim (1-based INCLUSIVE), so it must not reuse the
    0-based half-open gloss — under that convention its own Len contradicts its Coords (end-start+1 == Len).
    Its 'Source' cell is an assembly sequence name, not a specimen, and needs its own gloss too."""
    import main
    from teagle_core import genomepcr
    a = genomepcr.parse_ispcr(">chr7:5530601+5530800 pair 200bp ACGTACGTACGTACGTAC TTTTGGGGCCCCAAAATT\nACGT\n")[0]
    assert a["end"] - a["start"] + 1 == a["length"]                  # 1-based inclusive, as genomepcr documents
    from widgets import DataTable
    win._render_genome_scan([{"label": "P1", "amplicons": [a]}], [a], {"verdict": "v", "n_pair": 1})
    tbl = [t for i in range(win.genomeBody.count())
           for t in ([win.genomeBody.itemAt(i).widget()] if isinstance(win.genomeBody.itemAt(i).widget(), DataTable) else [])][0]
    heads = [tbl.horizontalHeaderItem(c).text() for c in range(tbl.columnCount())]
    assert "Coords (1-based)" in heads and "Coords" not in heads     # not the 0-based half-open key
    assert "Source" not in heads and "Assembly seq" in heads
    g = main.GLOSS
    assert "1-based INCLUSIVE" in g["Coords (1-based)"] and "0-based" not in g["Coords (1-based)"]
    assert "contig accession" in g["Assembly seq"] and "specimen" not in g["Assembly seq"]
    assert "0-based half-open" in g["Coords"]                        # panel 05 (simulated PCR) is unchanged


def test_struct_menu_copies_the_span_the_row_shows(win):
    """An LTR/TIR row lists its two arms (Len = one arm), so the row menu must copy an ARM — not silently
    substitute element_span (the whole element; for a TIR, the entire input). Whole element stays available."""
    from teagle_core import structural
    import random
    random.seed(7)
    core = "".join(random.choice("ACGT") for _ in range(2000))
    ltr = "".join(random.choice("ACGT") for _ in range(300))
    seq = "AAAA" + ltr + core + ltr + "AAAA"
    e = [x for x in structural.detect_all(seq) if x["type"].startswith("LTR")][0]
    win.state["analyzed_seq"] = seq
    row = win._struct_row(e)
    labels = [l for l, _ in win._struct_menu(e)]
    assert "Copy whole element FASTA" in labels                      # capability kept, clearly labelled
    m = dict(win._struct_menu(e))
    from PySide6.QtWidgets import QApplication
    m["Copy 5′ arm FASTA"]()
    txt = QApplication.clipboard().text()
    body = "".join(txt.splitlines()[1:])
    assert len(body) == row[2] and str(e["five_prime"][0]) in txt.splitlines()[0]   # copied bp == the row's Len
    m["Copy whole element FASTA"]()
    assert len("".join(QApplication.clipboard().text().splitlines()[1:])) == e["element_span"][1] - e["element_span"][0]


def test_struct_row_fills_method_and_len_and_hedges_pbs(win):
    """Method is glossed 'How TEagle detected this feature' — it must not be empty for PBS/PPT/TSD/poly-A,
    a PBS must show a Len, and a PBS below the confident threshold must read as tentative in the TABLE."""
    from teagle_core import structural
    import random
    random.seed(11)
    div = list(structural._PRIMER_TRNA["tRNA-Lys3"])
    for i in (2, 5, 8, 11, 14, 17):
        div[i] = {"A": "C", "C": "A", "G": "T", "T": "G"}[div[i]]
    leader = "".join(random.choice("ACGT") for _ in range(6)) + "".join(div) + "".join(random.choice("ACGT") for _ in range(20))
    seq = leader + "".join(random.choice("ACGT") for _ in range(600)) + "AGGGAGGGAGGGAAGGGA"
    pbs = structural.find_pbs(seq, 0)
    assert pbs and pbs["confident"] is False
    row, tips = win._struct_row(pbs), win._struct_tips(pbs)
    assert row[2] == pbs["pos"][1] - pbs["pos"][0]                   # Len no longer blank
    assert "tentative" in row[3] and row[4]                          # hedge visible + Method non-empty
    assert "undetermined" in (tips[3] or "")                         # full hedge on hover, not only in the viewer
    for e in [structural.find_ppt(seq, len(seq)),
              structural.find_tsd("AAAA" + "ACGT" * 30 + "AAAA", 4, 124)] + structural.find_polya("ACGT" * 20 + "A" * 12):
        assert e and win._struct_row(e)[4], e["type"]


def test_methods_panel_documents_pbs_and_ppt(win):
    # the panel claims the annotation is never a black box; PBS/PPT are shown as rows + viewer tracks
    html = win._methods_html()
    assert "PBS" in html and "PPT" in html
    assert "72%" in html and "82% purine" in html and "tentative" in html


def test_genome_manager_error_sizes_to_its_own_content(win):
    """The error rebuild returns before the table fit-tail; the dialog is reused, so it must release the
    previous successful build's fixed size instead of stranding one line of red text in a screenful."""
    win._genome_mgr_open = True
    win._on_genome_list({"ok": True, "genomes": []})
    big = win._genome_mgr.size()
    win._on_genome_list({"ok": False, "error": "WSL backend unavailable"})
    err = win._genome_mgr
    assert err.size().height() < big.height() and err.size().width() <= big.width()
    assert err.size().height() <= err.sizeHint().height() + 4        # no dead space below the content
    win._on_genome_list({"ok": True, "genomes": []})                 # and the table build recovers its own size
    assert win._genome_mgr.size().height() == big.height()


def test_wsl_broken_shows_error_not_crash(win):
    win._on_wsl_status({"error": "wsl.exe returned 0xffffffff"})
    assert "error" in win.wslStatus.text().lower()
    assert not win.annotateBtn.isEnabled()


def test_wsl_stack_missing_offers_install(win):
    win._on_wsl_status({"wsl2": True, "ready": False, "distro": "Ubuntu-24.04", "repeatmasker": None, "dfam": False})
    assert not win.wslInstallBtn.isHidden() and not win.annotateBtn.isEnabled()


def test_wsl_ready_enables_optional(win):
    win._on_wsl_status({"wsl2": True, "ready": True, "distro": "Ubuntu-24.04", "repeatmasker": "4.2.4", "minimap2": "2.28"})
    assert win.annotateBtn.isEnabled() and win.spliceBtn.isEnabled()
    assert not win.wslInstallBtn.isHidden()                    # installer stays reachable for repair / integrity checks


def test_partial_uncurated_install_never_exposes_combined_library_choice(win):
    """One optional shard is not a searchable curated+uncurated installation."""
    win._on_wsl_status({
        "wsl2": True,
        "ready": True,
        "distro": "Ubuntu-24.04",
        "repeatmasker": "4.2.4",
        "dfam": True,
        "dfam_uncurated": False,
        "dfam_library": {"partitions": ["dfam40.0.h5", "dfam40.uncurated.consensus.0.h5"]},
    })
    assert win._uncurated_available is False
    assert win.wslLibraryRow.isHidden()
    assert "curated + uncurated available" not in win.wslStatus.text().lower()


def test_core_science_independent_of_wsl():
    """The engine's core path must not touch WSL — analyze/primers/pcr work with WSL absent or broken."""
    import engine
    fx = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "tests", "fixtures", "M11240.fasta")
    r = engine.run_analyze({"sequence": open(fx).read()})
    assert r["records"][0]["classification"]["te_class"] == "LTR/Copia"


def test_annotate_without_stack_returns_soft_error():
    """run_annotate must degrade gracefully — never an unhandled crash. If the WSL annotation stack
    is absent it returns a soft error (ok:False / error); if the stack happens to be installed on
    this machine it succeeds cleanly. Either way it is a well-formed dict, never an exception."""
    import engine
    from engine import BadRequest
    from teagle_core import wsl
    try:
        res = engine.run_annotate({"sequence": ">x\n" + "ACGT" * 50})
    except BadRequest:
        return
    assert isinstance(res, dict)
    if wsl.env_status().get("ready"):
        assert res.get("ok") is True                            # stack present -> clean success
    else:
        assert res.get("ok") in (False, None) or "error" in res  # stack absent -> soft error, no crash


def test_amplicon_call_and_fasta_kind_are_one_disjoint_label():
    """A self-priming product is its OWN call everywhere — table, gel FASTA header, caption counts.
    It used to read 'off-target · single-primer' in the table while the caption counted it only as
    single-primer, so counting off-target rows disagreed with the caption."""
    import main
    off = {"on_target": False, "single_primer": False}
    single = {"on_target": False, "single_primer": True, "start": 100, "end": 400, "length": 300, "seq": "AC"}
    assert main._amp_call(off) == "off-target" and main._amp_call(single) == "single-primer"
    assert main._amp_call(off, has_locus=False) == "priming site"
    assert main._amp_kind(single) == "singleprimer" and main._amp_kind(off) == "offtarget"
    assert main._amp_kind(single, has_locus=False) == "singleprimer"     # never '_offtarget'/'_primingsite'


def test_delta_g_cell_exports_bare_number():
    """The ! / !! severity marks are on-screen cues; an exported ΔG column must stay numeric."""
    import main, widgets
    txt, flag, tip, export = main._metric_cell([("F", {"p3": -5.2, "vrna": -5.4, "flag": "caution"})])
    assert txt == "-5.4!" and float(export) == -5.4
    t = widgets.DataTable(["ΔG"])
    t.set_rows([[txt]], exports=[[export]])
    assert t.item(0, 0).text() == txt and t.rows_data() == [[export]]    # display marked, export numeric


def test_pcr_queue_buttons_are_uppercased_on_every_mutation(win):
    """The queue is rebuilt on add/move/remove/clear; its row buttons must not sit in sentence case
    until the next theme/scale refresh."""
    c = {"id": "P1", "left_seq": "ACGTACGTACGTACGTACGT", "right_seq": "TGCATGCATGCATGCATGCA",
         "product_size": 200, "left_pos": [0, 20], "right_pos": [180, 200]}
    from PySide6.QtWidgets import QPushButton
    win._add_pcr_pair(c)
    labels = [b.text() for b in win.pcrQueueBox.parentWidget().findChildren(QPushButton)]
    assert labels and all(l == l.upper() for l in labels), labels


def test_pcr_queue_missing_display_metadata_stays_recoverable(win):
    """A partial/stale row must render an honest placeholder instead of raising in a startup callback."""
    from PySide6.QtWidgets import QLabel
    win.state["pcrPairs"] = [{"left_seq": "ACGT" * 5, "right_seq": "TGCA" * 5}]
    win._render_pcr_queue()
    labels = [label.text() for label in win.pcrQueueBox.parentWidget().findChildren(QLabel)]
    assert any("product size unavailable" in text for text in labels)


def test_frozen_selftest_exercises_the_certifi_tls_context():
    """The frozen launcher bypasses envcheck, so its self-test must prove the CA bundle is usable."""
    import inspect
    import main
    source = inspect.getsource(main.selftest)
    assert "certifi.where()" in source
    assert "fetch._ssl_ctx()" in source
    assert "get_ca_certs()" in source


def test_frozen_selftest_exercises_the_actual_pdf_figure_export_path():
    """Importability is insufficient: the frozen build must render, clean and read a tiny figure."""
    import inspect
    import main
    source = inspect.getsource(main.selftest)
    for expected in (
            "TemporaryDirectory", "widgets.render_pdf",
            "media_metadata_violations", "PdfReader", "len(reader.pages)",
    ):
        assert expected in source
    assert "from pypdf import PdfReader, PdfWriter" not in source  # a blank PDF bypasses the runtime figure writer


def _budget_result(mode, request_id, families=None):
    result = {"ok": True, "cores": 8, "mem_gb": 16, "avail_gb": 100,
              "recommended_threads": 7, "limited_by": "cores", "disk_needed_gb": 2,
              "disk_ok": True, "library_mode": mode, "budget_request_id": request_id}
    if families is not None:
        result["library_families_for_species"] = families
    return result


def test_annotation_budget_discloses_and_enforces_the_windows_backing_volume_limit(win):
    from PySide6.QtWidgets import QLabel, QPushButton, QComboBox
    cost, warn, start, gsel, ssel = QLabel(), QLabel(), QPushButton(), QComboBox(), QComboBox()
    warn.hide()
    win._annot_cost_widgets = (cost, warn, start, gsel, ssel)
    win._annot_cost_state = {"organism": "Saccharomyces cerevisiae", "custom_library": None,
                             "library_mode": "dfam_curated", "budget_request_id": 10}
    result = _budget_result("dfam_curated", 10, families=9)
    result.update({"avail_gb": 35.0, "wsl_avail_gb": 950.0, "windows_avail_gb": 35.0,
                   "disk_limited_by": "Windows backing volume", "disk_needed_gb": 42,
                   "disk_ok": False,
                   "disk_note": ("Usable backend storage: 35.0 GiB (Windows backing volume). "
                                 "WSL reports 950.0 GiB; its Windows backing volume has 35.0 GiB free.")})
    win._on_annotate_budget(result)
    rendered = (cost.text() + " " + warn.text()).lower()
    assert "wsl reports 950.0 gib" in rendered
    assert "windows backing volume has 35.0 gib free" in rendered
    assert "smaller wsl filesystem and windows backing-volume limit" in rendered
    assert not start.isEnabled() and not warn.isHidden()


def test_custom_annotation_budget_has_no_dfam_copy_or_warning(win):
    """A custom-library estimate must carry only the library that will actually be searched."""
    from PySide6.QtWidgets import QLabel, QPushButton, QComboBox
    cost, warn, start, gsel, ssel = QLabel(), QLabel(), QPushButton(), QComboBox(), QComboBox()
    warn.hide()
    win._annot_cost_widgets = (cost, warn, start, gsel, ssel)
    win._annot_cost_state = {"organism": "Saccharomyces cerevisiae",
                             "custom_library": r"C:\libraries\lab.fa",
                             "library_mode": "custom", "budget_request_id": 11}
    win._on_annotate_budget(_budget_result("custom", 11, families=9))
    rendered = (cost.text() + " " + warn.text()).lower()
    assert "lab.fa" in rendered
    assert "dfam" not in rendered and "family model" not in rendered
    assert warn.isHidden() and start.isEnabled()


def test_stale_dfam_budget_cannot_overwrite_custom_selection(win):
    from PySide6.QtWidgets import QLabel, QPushButton, QComboBox
    cost, warn, start, gsel, ssel = QLabel("estimating custom"), QLabel(), QPushButton(), QComboBox(), QComboBox()
    warn.hide()
    start.setEnabled(False)
    win._annot_cost_widgets = (cost, warn, start, gsel, ssel)
    win._annot_cost_state = {"organism": "Saccharomyces cerevisiae",
                             "custom_library": r"C:\libraries\lab.fa",
                             "library_mode": "custom", "budget_request_id": 12}
    win._on_annotate_budget(_budget_result("dfam_curated", 11, families=9))
    assert cost.text() == "estimating custom"
    assert warn.isHidden() and not start.isEnabled()


def test_dfam_budget_warnings_name_the_selected_search_set(win):
    from PySide6.QtWidgets import QLabel, QPushButton, QComboBox

    def render(mode, request_id):
        cost, warn, start, gsel, ssel = QLabel(), QLabel(), QPushButton(), QComboBox(), QComboBox()
        win._annot_cost_widgets = (cost, warn, start, gsel, ssel)
        win._annot_cost_state = {"organism": "Saccharomyces cerevisiae", "custom_library": None,
                                 "library_mode": mode, "budget_request_id": request_id}
        win._on_annotate_budget(_budget_result(mode, request_id, families=9))
        return (cost.text() + " " + warn.text()).lower()

    curated = render("dfam_curated", 13)
    uncurated = render("dfam_uncurated", 14)
    assert "curated-only" in curated and "uncurated families will not be searched" in curated
    assert "curated + uncurated" in uncurated and "both family sets are already being searched" in uncurated


def test_missing_uncurated_partitions_block_start_and_never_claim_both_sets(win):
    from PySide6.QtWidgets import QLabel, QPushButton, QComboBox
    cost, warn, start, gsel, ssel = QLabel(), QLabel(), QPushButton(), QComboBox(), QComboBox()
    win._annot_cost_widgets = (cost, warn, start, gsel, ssel)
    win._annot_cost_state = {"organism": "Saccharomyces cerevisiae", "custom_library": None,
                             "library_mode": "dfam_uncurated", "budget_request_id": 19}
    result = _budget_result("dfam_uncurated", 19)
    result.update({"uncurated_available": False, "library_ready": False,
                   "library_error": "the optional uncurated Dfam partitions are not installed"})
    win._on_annotate_budget(result)
    rendered = (cost.text() + " " + warn.text()).lower()
    assert "not installed" in rendered
    assert not start.isEnabled() and not warn.isHidden()
    assert "both family sets are already being searched" not in rendered


def test_uncurated_without_lineage_shows_lineage_action_not_install_action(win):
    from PySide6.QtWidgets import QLabel, QPushButton, QComboBox
    cost, warn, start, gsel, ssel = QLabel(), QLabel(), QPushButton(), QComboBox(), QComboBox()
    win._annot_cost_widgets = (cost, warn, start, gsel, ssel)
    win._annot_cost_state = {"organism": None, "custom_library": None,
                             "library_mode": "dfam_uncurated", "budget_request_id": 20}
    result = _budget_result("dfam_uncurated", 20)
    result.update({"library_ready": False,
                   "library_error": "select an organism/lineage before including uncurated Dfam families"})
    win._on_annotate_budget(result)
    rendered = (cost.text() + " " + warn.text()).lower()
    assert "select an organism" in rendered and "not installed" not in rendered
    assert not start.isEnabled()


def test_custom_budget_copy_gate_rejects_injected_dfam_claim():
    """Mutation negative: an old Dfam sentence must fail the custom-library copy invariant."""
    def assert_custom_copy(text):
        lowered = text.lower()
        assert "lab.fa" in lowered and "dfam" not in lowered and "family model" not in lowered

    assert_custom_copy("Selected library: user-supplied FASTA lab.fa. Only represented sequences can match.")
    with pytest.raises(AssertionError):
        assert_custom_copy("Selected library: lab.fa. Dfam holds 9 family models for this lineage.")


def test_annotation_dialog_refreshes_budget_for_each_selected_library(win, monkeypatch, tmp_path):
    """Drive the non-modal dialog offscreen: custom and each Dfam search set issue a fresh typed budget."""
    from PySide6.QtWidgets import QComboBox, QFileDialog, QLabel
    import main
    library = tmp_path / "lab.fa"
    library.write_text(">LabTE\nACGT\n", encoding="utf-8")
    win._prepared_genomes = [{"accession": "GCF_000146045.2", "n_seqs": 17,
                              "bytes": 12_000_000, "sha256": "a" * 64}]
    submitted = []
    win.engine.submit = lambda op, body=None, key=None: submitted.append((op, body, key))
    monkeypatch.setattr(QFileDialog, "getOpenFileName", lambda *a, **k: (str(library), "FASTA"))

    win._annotate_genome_dialog()
    dialog = win._annot_dialog
    try:
        library_selector = next(
            combo for combo in dialog.findChildren(QComboBox)
            if combo.findData("custom") >= 0 and combo.findData("uncurated") >= 0)
        assert submitted[-1][1]["library_mode"] == "dfam_curated"

        library_selector.setCurrentIndex(library_selector.findData("custom"))
        assert submitted[-1][1]["library_mode"] == "custom"
        custom_id = submitted[-1][1]["budget_request_id"]

        library_selector.setCurrentIndex(library_selector.findData("uncurated"))
        assert submitted[-1][1]["library_mode"] == "dfam_uncurated"
        assert submitted[-1][1]["budget_request_id"] > custom_id

        intro = " ".join(label.text() for label in dialog.findChildren(QLabel)
                         if label.objectName() == "orient")
        assert "selected library" in intro.lower()
        assert "selected dfam" not in intro.lower()
    finally:
        dialog.close()
