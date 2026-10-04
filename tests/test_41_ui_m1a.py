"""4.1 UI shell: request generations, fetch races, invalid-input messages, record identity binding, failure
recovery and keyboard/wording fixes. Drives the real MainWindow offscreen; engine calls are captured, never sent."""
import os, sys, time

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
import pytest

pytest.importorskip("PySide6")
_NATIVE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "app", "native")
sys.path.insert(0, _NATIVE)

from PySide6.QtWidgets import QApplication, QLabel, QPushButton
from PySide6.QtCore import QEventLoop, QTimer

_app = QApplication.instance() or QApplication([])


def _spin(pred, timeout=30):
    loop = QEventLoop(); t = QTimer(); t.setInterval(20); t0 = time.time()
    t.timeout.connect(lambda: (pred() or time.time() - t0 > timeout) and loop.quit())
    t.start(); loop.exec()
    return pred()


@pytest.fixture
def win():
    import gc
    import main
    from PySide6.QtCore import QEvent
    w = main.MainWindow()
    _app.processEvents()                                              # run the window's deferred startup now, not inside a later test
    yield w
    w.engine.wait(5000)
    w.close(); w.deleteLater()
    QApplication.sendPostedEvents(None, QEvent.DeferredDelete)
    del w
    gc.collect()


def _capture(win, ops=("annotate", "splice", "primers", "pcr", "fetch", "fetch_coords", "genome_pcr")):
    """Capture the named engine submissions; every other op still runs for real."""
    calls = []
    real = win.engine.submit

    def submit(op, body=None, key=None):
        if op in ops:
            calls.append((op, body, key or op))
        else:
            real(op, body, key=key)
    win.engine.submit = submit
    return calls


def _analyse(win, text):
    win._set_seq(text)
    win._run_analysis()
    assert _spin(lambda: win.state.get("analysis_state") == "current"), "analysis did not complete"


def _label_texts(layout_owner):
    return " ".join(l.text() for l in layout_owner.findChildren(QLabel))


def _two_specimens():
    from teagle_core import examples
    return examples.load("M11240"), examples.load("X01005")


_FAMILY = {"ok": True, "hits": [{"class_family": "LTR/Copia", "family": "Copia_LTR", "q_start": 10, "q_end": 200,
                                  "strand": "+", "divergence": 1.0, "score": 500}],
           "provenance": {"runType": "annotate", "input": {"id": "stale_old_run", "length": 5146}},
           "species": "Drosophila melanogaster", "repeatmasker_version": "4.1"}


def test_late_family_result_after_specimen_change_is_dropped(win):
    a, b = _two_specimens()
    _analyse(win, a)
    calls = _capture(win)
    win._annotate()
    old_key = calls[-1][2]
    _analyse(win, b)
    assert win.annotateBtn.text().upper() != "◴ ANNOTATING…"
    win._on_done(old_key, _FAMILY)
    assert not win.state.get("family")
    assert "Copia_LTR" not in _label_texts(win.wslBody.parentWidget())
    assert win.state["prov_manifest"]["runType"] == "analysis"


def test_late_splice_result_after_specimen_change_is_dropped(win):
    a, b = _two_specimens()
    _analyse(win, a)
    calls = _capture(win)
    win.spliceTx.setPlainText("ACGT" * 30)
    win._splice()
    old_key = calls[-1][2]
    _analyse(win, b)
    assert win.spliceBtn.text().upper() == "DETECT EXONS / INTRONS"
    win._on_done(old_key, {"ok": True, "counts": {"exons": 2, "introns": 1}, "exons": [], "introns": [],
                           "provenance": {"runType": "splice", "input": {"id": "stale_old_run"}}})
    assert not win.state.get("splice")
    assert win.state["prov_manifest"]["runType"] == "analysis"


def test_late_primers_after_specimen_change_never_render_current(win):
    a, b = _two_specimens()
    _analyse(win, a)
    calls = _capture(win)
    win._design()
    old_key = calls[-1][2]
    assert win._design_inflight
    _analyse(win, b)
    assert not win._design_inflight and win.designBtn.text().upper() == "DESIGN PRIMERS"
    win._on_done(old_key, {"candidates": [{"id": "P1", "left_seq": "A" * 20, "right_seq": "C" * 20, "product_size": 200,
                                           "left_tm": 60, "right_tm": 60, "left_gc": 50, "right_gc": 50,
                                           "penalty": 0.1}],
                           "provenance": {"runType": "primers", "input": {"id": "stale_old_run"}}})
    assert not win.state.get("candidates") and not win.state.get("lastPrimers")
    assert win.state["prov_manifest"]["runType"] == "analysis"
    assert not win.pcrStageAll.isEnabled()


def test_stale_primer_failure_does_not_reset_the_newer_request(win):
    a, _ = _two_specimens()
    _analyse(win, a)
    calls = _capture(win)
    win._design()
    old_key = calls[-1][2]
    win._invalidate_analysis("changed", cancel_request=True)
    _analyse(win, a)
    win._design()
    assert win._design_inflight and len(calls) == 2
    win._on_user_error(old_key, "old failure")
    win._on_failed(old_key, "RuntimeError: old", "trace")
    assert win._design_inflight and win.designBtn.text().upper() == "◴ DESIGNING…"


def test_pcr_batch_in_flight_is_dropped_on_specimen_change(win):
    a, b = _two_specimens()
    _analyse(win, a)
    calls = _capture(win)
    win.state["pcrPairs"] = [{"id": "P1", "left_seq": "ACGTACGTACGTACGTACGT", "right_seq": "TTGCAATTGCAATTGCAATT"}]
    win._run_pcr()
    old_key = calls[-1][2]
    assert win.runPcrBtn.text().upper() == "◴ RUNNING…"
    _analyse(win, b)
    assert win.runPcrBtn.text().upper() == "RUN LOADED PAIRS"
    win._on_done(old_key, {"amplicons": [], "status": "complete",
                           "provenance": {"runType": "pcr", "input": {"id": "stale_old_run"}}})
    assert not win.state.get("lastPcr")
    assert win.state["prov_manifest"]["runType"] == "analysis"


_PCR_PAIR = {"id": "P1", "left_seq": "ACGTACGTACGTACGTACGT", "right_seq": "TTGCAATTGCAATTGCAATT"}
_PCR_DONE = {"amplicons": [], "status": "complete", "provenance": {"runType": "pcr", "input": {"id": "x"}}}


def _pcr_on_two_records(win):
    a, b = _two_specimens()
    _analyse(win, a + "\n" + b)
    calls = _capture(win)
    win.state["pcrPairs"] = [dict(_PCR_PAIR)]
    win._run_pcr()
    return calls, calls[-1][2]


def test_pcr_batch_for_another_record_is_dropped_and_the_button_restored(win):
    calls, key = _pcr_on_two_records(win)
    assert win.runPcrBtn.text().upper() == "◴ RUNNING…"
    win._show_record(1)
    assert win.runPcrBtn.text().upper() == "RUN LOADED PAIRS" and win.runPcrBtn.isEnabled()
    win._on_done(key, _PCR_DONE)
    assert not win.state.get("lastPcr")
    assert win.state["prov_manifest"]["runType"] == "analysis"
    n = len(calls)
    win._run_pcr()
    assert len(calls) == n + 1


def test_pcr_batch_survives_reselecting_its_own_record(win):
    calls, key = _pcr_on_two_records(win)
    win._show_record(0, keep_table=True)
    assert win.runPcrBtn.text().upper() == "◴ RUNNING…"
    win._on_done(key, _PCR_DONE)
    assert win.state.get("lastPcr")


def test_a_finished_pcr_for_another_record_does_not_stay_on_screen_as_current(win):
    calls, key = _pcr_on_two_records(win)
    win._on_done(key, _PCR_DONE)
    assert win.state.get("lastPcr")
    win._show_record(1)
    assert not win.state.get("lastPcr")
    assert "gel" not in _label_texts(win.pcrBody.parentWidget()).lower() and "lane(s)" not in _label_texts(win.pcrBody.parentWidget())


def test_genome_scan_started_for_previous_specimen_is_discarded(win):
    a, b = _two_specimens()
    _analyse(win, a)
    calls = _capture(win)
    cand = {"id": "P1", "left_seq": "ACGTACGTACGTACGTACGT", "right_seq": "TTGCAATTGCAATTGCAATT", "_tmpl_sig": "x"}
    win._scan_genome(cand, org="Homo sapiens")
    assert win._genome_inflight
    _analyse(win, b)
    win._on_done(calls[-1][2], {"ok": True, "amplicons": [], "summary": {"has_locus": False, "verdict": "no sites"},
                                "provenance": {"runType": "genome_pcr", "input": {"id": "stale_old_run"}}})
    assert not win._genome_inflight
    assert not win.state.get("lastGenomeScan")
    assert win.state["prov_manifest"]["runType"] == "analysis"


# ---------- fetch races ----------
_FETCH_OK = {"ok": True, "fasta": ">M11240.1\nACGTACGTACGT", "accession": "M11240.1", "organism": "Drosophila melanogaster",
             "length": 12, "source": "NCBI nuccore", "taxid": "7227",
             "endpoint": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
             "retrievedUtc": "2026-08-21T00:00:00+00:00",
             "sourceUrl": "https://www.ncbi.nlm.nih.gov/nuccore/M11240.1"}


def test_late_fetch_does_not_overwrite_a_pasted_sequence(win):
    calls = _capture(win)
    win.acc.setText("M11240")
    win._fetch()
    key = calls[-1][2]
    win.seq.setPlainText(">my_clone\n" + "ACGT" * 20)
    assert win.fetchBtn.isEnabled()
    win._on_done(key, _FETCH_OK)
    assert win.seq.toPlainText().startswith(">my_clone")
    assert win.state.get("source") is None


def test_late_fetch_does_not_overwrite_a_loaded_example(win):
    from teagle_core import examples
    calls = _capture(win)
    win.acc.setText("M11240")
    win._fetch()
    win._apply_example(examples.load("X01005"))
    win._on_done(calls[-1][2], _FETCH_OK)
    assert win.seq.toPlainText().startswith(">X01005")


def test_cancelling_a_fetch_restores_the_controls_and_ignores_its_result(win):
    calls = _capture(win)
    win.seq.setPlainText(">mine\n" + "ACGT" * 10)
    win.accMeta.setText("previous label")
    win.acc.setText("M11240")
    win._fetch()
    assert win.accMeta.text() == "fetching…" and not win.fetchBtn.isEnabled() and not win.fetchCancelBtn.isHidden()
    win.fetchCancelBtn.click()
    assert win.fetchBtn.isEnabled() and win.fetchCancelBtn.isHidden() and win.accMeta.text() == "previous label"
    win._on_done(calls[-1][2], _FETCH_OK)
    assert win.seq.toPlainText().startswith(">mine")


def test_current_fetch_still_loads(win):
    calls = _capture(win)
    win.acc.setText("M11240")
    win._fetch()
    win._on_done(calls[-1][2], _FETCH_OK)
    assert win.seq.toPlainText().startswith(">M11240.1")
    assert win.state["source"]["accession"] == "M11240.1" and win.fetchBtn.isEnabled()


def test_failed_fetch_keeps_the_loaded_specimen_label_and_binding(win):
    win.accMeta.setText("M11240.1 · copia")
    win.state["coord_region"] = {"marker": 1}
    calls = _capture(win)
    win.acc.setText("BAD")
    win._fetch()
    win._on_done(calls[-1][2], {"ok": False, "error": "unknown accession"})
    assert win.accMeta.text() == "M11240.1 · copia"
    assert win.state["coord_region"] == {"marker": 1}


def test_enter_in_the_accession_box_fetches(win):
    calls = _capture(win)
    win.acc.setText("M11240")
    win.acc.returnPressed.emit()
    assert calls and calls[-1][0] == "fetch"


def test_replaced_record_status_is_shown_beside_the_label(win):
    win._on_fetch({**_FETCH_OK, "status": "replaced", "replacedBy": "X04456"})
    assert "replaced by X04456" in win.accMeta.text()


# ---------- stuck label, whitespace edits, tab focus ----------
def test_run_button_label_resets_when_the_specimen_changes_mid_run(win):
    a, b = _two_specimens()
    win.seq.setPlainText(a)
    win._run_analysis()
    assert win.runBtn.text().upper().endswith("ANALYSING")
    win._set_seq(b)
    assert win.runBtn.text().upper() == "RUN ANALYSIS" and win.runBtn.isEnabled()
    win.engine.wait()


def test_whitespace_only_edit_keeps_the_analysis(win):
    a, _ = _two_specimens()
    _analyse(win, a)
    cur = win.seq.textCursor(); cur.movePosition(cur.MoveOperation.End); win.seq.setTextCursor(cur)
    win.seq.insertPlainText("\t\n")
    assert win.state["analysis_state"] == "current" and win.designBtn.isEnabled()
    win.seq.insertPlainText("A")
    assert win.state["analysis_state"] == "stale"


def test_digit_typed_after_analysis_still_invalidates(win):
    a, _ = _two_specimens()
    _analyse(win, a)
    win.seq.insertPlainText("5")
    assert win.state["analysis_state"] == "stale"


def test_multiline_boxes_let_tab_leave(win):
    for box in (win.seq, win.spliceTx, win.wslPaste, win.pcrBg):
        assert box.tabChangesFocus()


# ---------- invalid input ----------
def test_genbank_numbered_paste_names_the_problem_not_a_seal_failure(win):
    win.seq.setPlainText("        1 acgtacgtac gtacgtacgt acgtacgtac\n       31 acgtacgtac gtacgtacgt acgtacgtac")
    win._run_analysis()
    assert _spin(lambda: win.runBtn.text().upper() == "RUN ANALYSIS")
    text = _label_texts(win.card_struct) + " " + _label_texts(win.card_prov)
    assert "provenance binding" not in text and "seal" not in text.lower()
    assert "GenBank" in text and "“1” at" in text
    assert win.state.get("analysis_state") != "current"
    assert not [b for b in win.card_prov.findChildren(QPushButton) if "EXPORT" in b.text().upper()]


def test_protein_paste_says_it_looks_like_protein(win):
    win.seq.setPlainText(">p\nMKVLAAGIVGLLLAQPEFSTEQWYPLIHFEEQKL")
    win._run_analysis()
    assert _spin(lambda: win.runBtn.text().upper() == "RUN ANALYSIS")
    assert "looks like a protein" in _label_texts(win.card_struct)


def test_accession_typed_into_the_sequence_box_points_to_fetch(win):
    win.seq.setPlainText("M11240")
    win._run_analysis()
    assert _spin(lambda: win.runBtn.text().upper() == "RUN ANALYSIS")
    assert "accession" in _label_texts(win.card_struct) and "Fetch" in _label_texts(win.card_struct)


def test_alignment_gaps_and_masking_are_named():
    from main import invalid_input_message
    msg = invalid_input_message([{"valid": False, "seq": "ACGT--ACGT", "invalid": []}])
    assert "alignment gaps" in msg and "“-” at 5" in msg
    msg = invalid_input_message([{"valid": False, "seq": "ACGTXXXXAC", "invalid": []}])
    assert "X usually marks masked" in msg
    msg = invalid_input_message([{"valid": False, "seq": "", "invalid": [],
                                  "invalidReason": "empty record (a header with no sequence)"}])
    assert "no sequence" in msg


def test_preamble_before_the_header_is_reported_with_the_invalid_record(win):
    win.seq.setPlainText("note to self\n>r1\nAC12GT")
    win._run_analysis()
    assert _spin(lambda: win.runBtn.text().upper() == "RUN ANALYSIS")
    assert "before the first FASTA header" in _label_texts(win.card_struct)


def test_valid_records_next_to_an_invalid_one_name_it_in_the_warning(win):
    from teagle_core import examples
    win.seq.setPlainText(examples.load("M11240") + "\n>bad\nAC12GT")
    win._run_analysis()
    assert _spin(lambda: win.state.get("analysis_state") == "current")
    assert "Record 2" in " ".join(l.text() for l in win._notif.findChildren(QLabel))


# ---------- rail ----------
def test_rail_tiles_reset_when_a_run_is_invalidated(win):
    a, _ = _two_specimens()
    _analyse(win, a)
    assert win.mLen.text() not in ("—", "") and win.rOrf.text() != "—"
    win.seq.setPlainText("MKVLAAGIVGLLLAQPEF")
    assert [w.text() for w in (win.mLen, win.mGC, win.mN, win.mValid, win.rRecords, win.rStruct, win.rOrf)] == ["—"] * 7
    assert win.mValid.property("state") in ("", None)


def test_counter_says_characters_until_a_record_validates(win):
    a, _ = _two_specimens()
    win.seq.setPlainText(">x\nMKVL AAGI 12")
    assert win.charCount.text().startswith("8 characters") and "2 non-letters" in win.charCount.text()
    _analyse(win, a)
    assert win.charCount.text().endswith(" nt") and "," in win.charCount.text()


def test_run_button_tooltip_follows_its_state(win):
    assert "first" in win.runBtn.toolTip()
    win.seq.setPlainText("ACGT" * 10)
    assert "first" not in win.runBtn.toolTip()


# ---------- placeholders, scroll ----------
def test_primer_panel_reads_ready_after_a_run_that_followed_an_invalidation(win):
    a, b = _two_specimens()
    _analyse(win, a)
    win._set_seq(b)
    assert "before designing" in _label_texts(win.primBody.parentWidget())
    _analyse(win, b)
    txt = _label_texts(win.primBody.parentWidget())
    assert "before designing" not in txt and "Ready" in txt


def test_new_analysis_scrolls_to_its_headline(win):
    a, b = _two_specimens()
    win.resize(1240, 560); win.show()
    _analyse(win, a)
    for card in (win.card_wsl, win.card_splice, win.card_primer, win.card_pcr, win.card_genome, win.card_prov):
        card.expand()
    _app.processEvents()
    bar = win.resultsScroll.verticalScrollBar()
    bar.setValue(bar.maximum())
    assert bar.value() > 0
    _analyse(win, b)
    _app.processEvents()
    assert bar.value() == win.card_struct.mapTo(win.resultsScroll.widget(), win.card_struct.rect().topLeft()).y()
    win.hide()


# ---------- organism, preset ----------
def test_example_sets_its_own_organism_and_a_paste_clears_it(win):
    win._load_example("X05424")
    assert win._species() == "Zea mays"
    win._load_example("M80343")
    assert win._species() == "Homo sapiens"
    win.seq.setPlainText(">mine\n" + "ACGT" * 20)
    assert win._species() is None


def test_a_chosen_organism_survives_a_specimen_change(win):
    win._load_example("X05424")
    win.wslSpecies.setCurrentIndex(win.wslSpecies.findData("Mus musculus"))
    win.seq.setPlainText(">mine\n" + "ACGT" * 20)
    assert win._species() == "Mus musculus"


def test_fetch_after_a_human_example_does_not_keep_human_for_a_record_without_organism(win):
    win._load_example("M80343")
    assert win._species() == "Homo sapiens"
    win._on_fetch({**_FETCH_OK, "organism": ""})
    assert win._species() is None


def test_reset_restores_standard_and_a_preset_can_be_reapplied(win):
    win.pPreset.setCurrentIndex(1); win.pPreset.activated.emit(1)
    assert win.pfields["pMin"].text() == "70"
    win.pfields["pMin"].setText("200"); win.pfields["pMin"].textEdited.emit("200")
    assert win.pPreset.currentText().startswith("Custom")
    win.pPreset.activated.emit(1)
    assert win.pfields["pMin"].text() == "70" and win.pPreset.currentIndex() == 1
    rb = [b for b in win.card_primer.findChildren(QPushButton) if b.accessibleName().startswith("Reset primer")][0]
    rb.click()
    assert win.pfields["pMin"].text() == "150" and win.pPreset.currentIndex() == 0
    win.pPreset.activated.emit(1)
    assert win.pfields["pMin"].text() == "70"


# ---------- record identity binding (B3) ----------
def _spy(win, ops=("pcr", "primers")):
    """Record the bodies of the named engine submissions and still run them for real."""
    bodies = []
    real = win.engine.submit

    def submit(op, body=None, key=None):
        if op in ops:
            bodies.append((op, body))
        real(op, body, key=key)
    win.engine.submit = submit
    return bodies


def _design_and_run_pcr(win):
    win._design()
    assert _spin(lambda: win.state.get("candidates")), "primer design did not return"
    win._pcr_stage_all()
    win._run_pcr()
    assert _spin(lambda: win.state.get("lastPcr")), "in-silico PCR did not return"


@pytest.mark.parametrize("index", [0, 1])
def test_primer_pcr_and_provenance_carry_the_selected_records_own_id(win, index):
    from main import _amplicon_public_source
    from teagle_core.gff3 import public_record_id
    a, b = _two_specimens()
    _analyse(win, a + "\n" + b)
    rec = win.state["records"][index]
    win._show_record(index)
    bodies = _spy(win)
    _design_and_run_pcr(win)
    d = win.state["lastPrimers"]["provenance"]["input"]
    assert (d["id"], d["idSha256"]) == (rec["id"], rec["idSha256"])
    assert d["id"] != public_record_id("input_sequence")
    pcr_body = [b for op, b in bodies if op == "pcr"][0]
    assert pcr_body["design_assessment"] == "assess"
    sources = {_amplicon_public_source(x) for x in win.state["lastPcr"]["amplicons"] if x.get("role", "input") == "input"
               or str(x.get("source", "")).endswith("(input)")}
    assert sources == {rec["id"] + " (input)"}
    assert win.state["prov_manifest_all"][0]["input"]["id"] == rec["id"]


def test_the_request_sequence_is_never_a_bare_header_that_leaks_into_a_result(win):
    import json
    a, _ = _two_specimens()
    _analyse(win, a)
    bodies = _spy(win)
    _design_and_run_pcr(win)
    assert bodies[0][1]["sequence"].startswith(">M11240.1")           # travels under its own header ...
    public = [{k: v for k, v in c.items() if not k.startswith("_")} for c in win.state["candidates"]]
    shown = (json.dumps(public) + json.dumps(win.state["lastPrimers"]["provenance"]) + json.dumps(win.state["lastPcr"])
             + json.dumps(win.state["prov_manifest_all"]))
    assert "white-apricot" not in shown                                 # ... which the engine hashes and never echoes


def test_bare_paste_keeps_the_generated_input_sequence_id(win):
    from teagle_core import sequtil, examples
    bare = sequtil.parse_fasta(examples.load("X01005"))[0][1]
    _analyse(win, bare)
    assert win.state["records"][0]["id"] == "input_sequence"
    _design_and_run_pcr(win)
    assert win.state["lastPrimers"]["provenance"]["input"]["id"] == "input_sequence"
    assert {x["source"] for x in win.state["lastPcr"]["amplicons"]} == {"input_sequence (input)"}


def test_a_header_that_does_not_reproduce_the_analysed_identity_is_not_trusted(win):
    a, _ = _two_specimens()
    _analyse(win, a)
    win.state["record_headers"] = ["some other header"]
    seq = win.state["records"][0]["seq"]
    assert win._operation_input()[0] == seq
    win.state["record_headers"] = ["M11240.1 D.melanogaster copia element inserted at white-apricot allele"]
    assert win._operation_input()[0].startswith(">M11240.1")
    assert win._operation_input(seq[:-1])[0] == seq[:-1]                # a different sequence never inherits the header


def test_domain_row_design_is_bound_to_the_record_and_assessed(win):
    a, _ = _two_specimens()
    _analyse(win, a)
    rec = win.state["records"][0]
    bodies = _spy(win, ops=("primers",))
    win._design_for_domain(1000, 1100, "INT", seq=rec["seq"])           # what the domain-table row menu does
    assert _spin(lambda: not win._design_inflight) and win.state["candidates"]
    win._design_for_domain(1000, 1100, "INT")                           # the genome-viewer route
    assert _spin(lambda: not win._design_inflight) and len(bodies) == 2
    assert bodies[0][1]["sequence"] == bodies[1][1]["sequence"] and bodies[0][1]["sequence"].startswith(">")
    assert bodies[0][1]["source"] == bodies[1][1]["source"]
    win._pcr_stage_all()
    pcr = _spy(win, ops=("pcr",))
    win._run_pcr()
    assert _spin(lambda: win.state.get("lastPcr"))
    assert pcr[0][1]["design_assessment"] == "assess"


_PRIMERS = {"candidates": [{"id": "P1", "left_seq": "A" * 20, "right_seq": "C" * 20, "product_size": 200,
                            "left_tm": 60, "right_tm": 60, "left_gc": 50, "right_gc": 50, "penalty": 0.1}]}
_SPLICE = {"ok": True, "counts": {"exons": 2, "introns": 1}, "exons": [], "introns": []}
_KINDS = {"primers": ("designBtn", "primBody", "DESIGN PRIMERS", _PRIMERS),
          "domain": ("designBtn", "primBody", "DESIGN PRIMERS", _PRIMERS),
          "annotate": ("annotateBtn", "wslBody", "RUN FAMILY ANNOTATION", _FAMILY),
          "splice": ("spliceBtn", "spliceBody", "DETECT EXONS / INTRONS", _SPLICE)}


def _start(win, kind):
    """Submit one tagged request on the current record; returns the key the engine would echo."""
    win._on_wsl_status = lambda w: None                                 # the real probe must not reset readiness mid-test
    win._wsl_ready.update(annotate=True, splice=True)
    win.annotateBtn.setEnabled(True); win.spliceBtn.setEnabled(True)
    calls = _capture(win)
    if kind == "primers":
        win._design()
    elif kind == "domain":
        win._design_for_domain(1000, 1100, "INT")
    elif kind == "annotate":
        win._annotate()
    else:
        win.spliceTx.setPlainText("ACGT" * 30)
        win._splice()
    assert len(calls) == 1
    return calls, calls[-1][2]


def _assert_idle(win, kind):
    from widgets import BusyBar
    btn, body, label, _ = _KINDS[kind]
    lay = getattr(win, body)
    first = lay.itemAt(0).widget() if lay.count() else None
    state = {"text": getattr(win, btn).text().upper(), "enabled": getattr(win, btn).isEnabled(),
             "busy_body": isinstance(first, BusyBar), "inflight": win._design_inflight}
    assert state == {"text": label, "enabled": True, "busy_body": False, "inflight": False}, state


def _fresh_request_goes_out(win, kind, calls):
    n = len(calls)
    {"primers": win._design, "domain": lambda: win._design_for_domain(1000, 1100, "INT"),
     "annotate": win._annotate, "splice": win._splice}[kind]()
    assert len(calls) == n + 1, "a new request was refused or never submitted"


@pytest.mark.parametrize("kind", list(_KINDS))
def test_a_request_for_another_record_releases_its_control_and_the_next_one_goes_out(win, kind):
    a, b = _two_specimens()
    _analyse(win, a + "\n" + b)
    calls, key = _start(win, kind)
    win._show_record(1)
    _assert_idle(win, kind)                                             # idle at once, not only when the late result lands
    win._on_done(key, _KINDS[kind][3])
    assert not win.state.get("candidates") and not win.state.get("family") and not win.state.get("splice")
    _assert_idle(win, kind)
    _fresh_request_goes_out(win, kind, calls)


@pytest.mark.parametrize("kind", list(_KINDS))
def test_retire_releases_a_request_it_drops_as_stale(win, kind):
    a, b = _two_specimens()
    _analyse(win, a + "\n" + b)
    calls, key = _start(win, kind)
    win.state["record_index"] = 1                                       # the record moved without going through _show_record
    win._on_done(key, _KINDS[kind][3])
    assert not win.state.get("candidates") and not win.state.get("family") and not win.state.get("splice")
    _assert_idle(win, kind)
    _fresh_request_goes_out(win, kind, calls)


@pytest.mark.parametrize("fail", ["user_error", "fault"])
@pytest.mark.parametrize("kind", list(_KINDS))
def test_a_stale_failure_releases_the_control_too(win, kind, fail):
    a, b = _two_specimens()
    _analyse(win, a + "\n" + b)
    calls, key = _start(win, kind)
    win.state["record_index"] = 1
    (win._on_user_error(key, "late failure") if fail == "user_error" else win._on_failed(key, "RuntimeError: late", "trace"))
    _assert_idle(win, kind)


@pytest.mark.parametrize("kind", ["annotate", "splice"])
def test_a_result_arriving_while_the_same_specimen_is_reanalysed_is_kept_and_shown(win, kind):
    a, _ = _two_specimens()
    _analyse(win, a)
    calls, key = _start(win, kind)
    win._run_analysis()                                                 # same box; the analysis is now in flight
    assert win.state.get("analysis_state") == "running" and win.state.get("record_index") is None
    win._on_done(key, _KINDS[kind][3])
    assert not win.state.get("family") and not win.state.get("splice")  # nothing renders against a record that is not shown yet
    assert _spin(lambda: win.state.get("analysis_state") == "current")
    assert bool(win.state.get("family" if kind == "annotate" else "splice")), "the result was lost in the re-analysis window"
    _assert_idle(win, kind)
    _fresh_request_goes_out(win, kind, calls)


def test_a_held_result_for_another_record_is_discarded_when_the_rerun_lands_elsewhere(win):
    a, b = _two_specimens()
    _analyse(win, a + "\n" + b)
    win._show_record(1)
    calls, key = _start(win, "annotate")
    win._run_analysis()                                                 # the re-analysis shows record 1 of 2, the first valid one
    win._on_done(key, _FAMILY)
    assert _spin(lambda: win.state.get("analysis_state") == "current")
    assert win.state.get("record_index") == 0 and not win.state.get("family")
    _assert_idle(win, "annotate")


def test_a_held_result_is_released_when_the_reanalysis_fails(win):
    a, _ = _two_specimens()
    _analyse(win, a)
    calls, key = _start(win, "annotate")
    win._run_analysis()
    win._on_done(key, _FAMILY)
    win._on_user_error(f"analyze#{win._active_analyze_gen}", "analysis failed")
    assert not win.state.get("family")
    _assert_idle(win, "annotate")
    _fresh_request_goes_out(win, "annotate", calls)


def test_a_held_result_dies_with_a_specimen_change(win):
    a, b = _two_specimens()
    _analyse(win, a)
    calls, key = _start(win, "annotate")
    win._run_analysis()
    win._on_done(key, _FAMILY)
    win._set_seq(b)
    _assert_idle(win, "annotate")
    _analyse(win, b)
    assert not win.state.get("family")


def test_primer_result_for_another_record_is_dropped_after_switching_records(win):
    a, b = _two_specimens()
    _analyse(win, a + "\n" + b)
    calls = _capture(win)
    win._design()
    key = calls[-1][2]
    win._show_record(1)
    win._on_done(key, _PRIMERS)
    assert not win.state.get("candidates")


# ---------- flank / self-similarity tools ----------
def test_flank_and_selfsim_windows_close_when_the_specimen_changes(win):
    a, b = _two_specimens()
    _analyse(win, a)
    win._flank_picker(100, 300, "+", "feature")
    win._build_selfsim(win.state["records"][0])
    flank, selfsim = win._flank_dlg, win._selfsim_dlg
    assert flank.isVisible() and selfsim.isVisible()
    win._set_seq(b)
    assert not flank.isVisible() and not selfsim.isVisible()


# ---------- provenance and exports ----------
def test_analysis_result_stays_exportable_after_primers_and_pcr(win, tmp_path, monkeypatch):
    import json, engine
    from PySide6.QtWidgets import QFileDialog
    a, _ = _two_specimens()
    _analyse(win, a)
    _design_and_run_pcr(win)
    assert win.state["prov_manifest"]["runType"] != "analysis"
    labels = [b.text().upper() for b in win.card_prov.findChildren(QPushButton)]
    assert any("ANALYSIS RESULT" in t for t in labels) and any("ANALYSIS MANIFEST" in t for t in labels)
    out = tmp_path / "analysis.json"
    monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *a, **k: (str(out), ""))
    win._export_analysis_result()
    doc = json.loads(out.read_text(encoding="utf-8"))
    assert engine.validate_analysis_result_document(doc)
    man = tmp_path / "manifest.json"
    monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *a, **k: (str(man), ""))
    win._export_analysis_manifest()
    assert json.loads(man.read_text(encoding="utf-8"))["runType"] == "analysis"


def test_a_write_failure_is_not_reported_as_a_seal_failure(win, tmp_path, monkeypatch):
    import widgets
    from PySide6.QtWidgets import QFileDialog
    a, _ = _two_specimens()
    _analyse(win, a)
    monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *a, **k: (str(tmp_path / "m.json"), ""))

    def locked(*a, **k):
        raise PermissionError("locked")
    monkeypatch.setattr(widgets, "atomic_write", locked)
    for fn in (win._export_manifest, win._export_analysis_result):
        fn()
        text = " ".join(l.text() for l in win._notif.findChildren(QLabel))
        assert "could not be written" in text and "seal" not in text and "binding" not in text


# ---------- failure recovery ----------
def test_failures_reset_every_busy_state_without_stderr(win, monkeypatch):
    monkeypatch.setattr(sys, "stderr", None)
    a, _ = _two_specimens()
    win.seq.setPlainText(a)
    calls = _capture(win)
    win._run_analysis()
    run_key = [k for k in [win._active_analyze_gen]][0]
    win._on_failed(f"analyze#{run_key}", "RuntimeError: x", "trace")
    assert win.runBtn.text().upper() == "RUN ANALYSIS" and not win._analyze_inflight
    _analyse(win, a)
    win._design()
    win._on_failed(calls[-1][2], "RuntimeError: x", "trace")
    assert not win._design_inflight and win.designBtn.isEnabled() and win.designBtn.text().upper() == "DESIGN PRIMERS"
    win._genome_prep_inflight = True
    win._on_failed("genome_prepare", "RuntimeError: x", "trace")
    assert not win._genome_prep_inflight
    win._on_failed("env", "RuntimeError: x", "trace")
    assert "unavailable" in win.envBox.text()
    win._on_failed("wsl_status", "RuntimeError: x", "trace")
    assert "could not be read" in win.wslStatus.text()


def test_failure_trace_is_logged_last_when_stderr_exists(win, monkeypatch):
    import io
    buf = io.StringIO()
    monkeypatch.setattr(sys, "stderr", buf)
    win._on_failed("env", "RuntimeError: x", "TRACE-LINE")
    assert "TRACE-LINE" in buf.getvalue()


def test_a_list_failure_does_not_leave_the_genome_manager_armed(win):
    win._genome_mgr_open = True
    win._on_failed("genome_list", "RuntimeError: x", "trace")
    assert win._genome_mgr_open is False


# ---------- external links ----------
def test_external_links_open_outside_the_kill_on_close_job(win, monkeypatch):
    import subprocess, main
    if sys.platform != "win32":
        pytest.skip("Windows job-object breakaway")
    seen = []
    monkeypatch.setattr(subprocess, "Popen", lambda *a, **k: seen.append((a, k)))
    assert main.open_external_url("https://doi.org/10.1038/nrg2165")
    args, kw = seen[0]
    assert args[0][:2] == ["rundll32.exe", "url.dll,FileProtocolHandler"] and kw["creationflags"] & 0x01000000
    assert not main.open_external_url("file:///c:/windows/system32/calc.exe") and len(seen) == 1
    assert not win.accMeta.openExternalLinks()
    win.accMeta.linkActivated.emit("https://www.ncbi.nlm.nih.gov/nuccore/M11240.1")
    assert len(seen) == 2


def test_external_link_falls_back_when_the_job_forbids_breakaway(monkeypatch):
    import subprocess, main
    from PySide6.QtGui import QDesktopServices
    if sys.platform != "win32":
        pytest.skip("Windows job-object breakaway")
    opened = []
    monkeypatch.setattr(subprocess, "Popen", lambda *a, **k: (_ for _ in ()).throw(PermissionError("job")))
    monkeypatch.setattr(QDesktopServices, "openUrl", lambda u: opened.append(u.toString()) or True)
    assert main.open_external_url("https://doi.org/10.1038/nrg2165") and opened == ["https://doi.org/10.1038/nrg2165"]


# ---------- self-test ----------
def test_frozen_selftest_fails_when_bundled_package_metadata_is_incomplete(monkeypatch, capsys):
    import main, envcheck
    monkeypatch.setattr(sys, "frozen", True, raising=False)
    monkeypatch.setattr(sys, "_MEIPASS", os.path.dirname(_NATIVE), raising=False)     # bundled resources resolve beside the sources
    monkeypatch.setattr(envcheck, "check", lambda: {"packages_ok": False, "importFailures": ["pyhmmer"],
                                                    "packages": [{"name": "primer3-py", "ok": False}]})
    assert main.selftest() == 1
    err = capsys.readouterr().err
    assert "bundled package check failed" in err and "primer3-py" in err and "pyhmmer" in err


def test_unfrozen_selftest_does_not_run_the_package_check(monkeypatch, capsys):
    import main, envcheck
    monkeypatch.setattr(envcheck, "check", lambda: (_ for _ in ()).throw(AssertionError("must not run")))
    assert main.selftest() == 0, capsys.readouterr().err
    assert "bundled package check" not in capsys.readouterr().err


# ---------- headline, summary table, wording ----------
def test_a_record_with_no_evidence_does_not_read_as_a_candidate(win):
    import random
    rng = random.Random(7)
    _analyse(win, ">noise\n" + "".join(rng.choice("ACGT") for _ in range(1800)))
    cl = win.state["last_rec"]["classification"]
    if cl.get("te_class") != "none":
        pytest.skip("this random record happened to carry structural evidence")
    names = {l.objectName(): l.text() for l in win.card_struct.findChildren(QLabel)}
    assert names["classbig"] == "No TE signature detected"
    assert not [l for l in win.card_struct.findChildren(QLabel) if l.objectName() == "cf"]


def test_record_summary_length_is_a_number_that_sorts_and_exports_as_one(win):
    from PySide6.QtCore import Qt
    from widgets import DataTable
    recs = [{"id": f"r{i}", "valid": True, "composition": {"length": n, "gc": 50.0, "n": 0.0, "ambiguous": 0.0},
             "classification": None, "processing": {}} for i, n in enumerate((5146, 980, 1200000))]
    rows = win._record_summary_rows(recs)
    assert [r[2] for r in rows] == [5146, 980, 1200000]
    t = DataTable(win._RECORD_HEADERS, {})
    t.set_rows(rows)
    t.sortByColumn(2, Qt.AscendingOrder)
    assert [t.rows_data()[i][2] for i in range(3)] == ["980", "5146", "1200000"]


def test_card_headers_keep_checked_state_and_a_speakable_name(win):
    win.show()                                                        # a card reads its own visibility, so it has to be on a shown window
    _app.processEvents()
    win.card_wsl.expand()
    assert win.card_wsl.hdr.isChecked() and "Collapsed" not in win.card_wsl.hdr.accessibleDescription()
    win.card_wsl.hdr.click()
    assert _spin(lambda: not win.card_wsl.body.isVisible(), 5)       # the human click eases the card closed
    assert not win.card_wsl.hdr.isChecked()
    assert win.card_wsl.hdr.accessibleName().startswith("03")


def test_example_menu_opens_at_its_button_and_parameter_fields_are_labelled(win, monkeypatch):
    import PySide6.QtWidgets as qw
    seen = []

    class _Menu(qw.QMenu):
        def exec(self, pos=None):
            seen.append(pos)
    monkeypatch.setattr(qw, "QMenu", _Menu)
    win._pick_example()
    assert seen == [win.exampleBtn.mapToGlobal(win.exampleBtn.rect().bottomLeft())]
    assert win.pfields["pMin"].accessibleName() == "Prod min (bp)"


def test_delta_g_glosses_state_both_tiers_from_the_qc_constants():
    import main
    from teagle_core import oligoqc
    for col, tier in (("Hairpin", oligoqc._HAIRPIN_CAUTION), ("Self-dim", oligoqc._DIMER_CAUTION),
                      ("Hetero", oligoqc._DIMER_CAUTION)):
        tip = main.GLOSS[col]
        assert f"≤ −{abs(tier):g}" in tip and f"≤ −{abs(oligoqc._HAIRPIN_WARN):g}" in tip and "caution" in tip


def test_wording_matches_what_the_software_does(win):
    import main
    texts = " ".join(l.text() for l in win.findChildren(QLabel))
    assert "every region is fetched, sealed and analysed" in texts and "first region only" in texts
    assert "Order and superfamily" in texts and "LINE / DNA)" not in texts
    note = main._qq_note().text()
    assert "rush mode (-qq)" in note and "quick mode" not in note


def test_status_chip_names_the_engine_not_the_backend(win):
    win._on_done("health", {"primer3": "2.3.0"})
    assert win.statusTxt.text() == "engine ready · primer3 2.3.0"
    win._on_wsl_status({"wsl2": False})
    assert "not installed" in win.backendBtn.toolTip()


# ---------- bring-your-own primer pair ----------
def _pair_from(seq):
    from main import MainWindow
    return seq[100:122], MainWindow._revcomp(seq[500:522])


def test_a_pair_you_already_have_can_be_queued_and_run_with_a_tail(win):
    a, _ = _two_specimens()
    _analyse(win, a)
    seq = win.state["records"][0]["seq"]
    f, r = _pair_from(seq)
    win.pcrCustomF.setText("GAATTC" + f); win.pcrCustomR.setText(r.lower()); win.pcrCustomName.setText("paper pair")
    win.pcrCustomAdd.click()
    pair = win.state["pcrPairs"][0]
    assert pair["custom"] and pair["id"] == "paper pair" and pair["left_seq"] == "GAATTC" + f
    assert not win.pcrCustomF.text()
    bodies = _spy(win, ops=("pcr",))
    win._run_pcr()
    assert _spin(lambda: win.state.get("lastPcr"))
    assert bodies[0][1]["design_assessment"] == "not_run" and bodies[0][1]["design_provenance"] is None
    sizes = {x.get("product_length", x.get("length")) for x in win.state["lastPcr"]["amplicons"]}
    assert 422 + 6 in sizes or 422 in sizes


def test_a_pair_with_a_bad_letter_or_too_short_is_refused_with_a_reason(win):
    win.pcrCustomF.setText("ACGT1ACGTACGTACGTACG"); win.pcrCustomR.setText("ACGTACGTACGTACGTACGT")
    win.pcrCustomAdd.click()
    assert not win.state.get("pcrPairs")
    assert "not IUPAC" in " ".join(l.text() for l in win._notif.findChildren(QLabel))
    win.pcrCustomF.setText("ACGTACGT"); win.pcrCustomAdd.click()
    assert not win.state.get("pcrPairs") and "at least" in " ".join(l.text() for l in win._notif.findChildren(QLabel))


def test_fetched_specimen_shows_its_accession_not_a_hash_in_the_pcr_table(win):
    from main import _amplicon_public_source
    from teagle_core import examples, sequtil
    seq = sequtil.parse_fasta(examples.load("M11240"))[0][1]
    win._on_fetch({**_FETCH_OK, "fasta": ">M11240.1\n" + seq, "length": len(seq)})
    _analyse_current(win)
    _design_and_run_pcr(win)
    labels = {_amplicon_public_source(a) for a in win.state["lastPcr"]["amplicons"]}
    assert labels == {"M11240.1 (input)"}
    assert win.state["prov_manifest_all"][0]["input"]["id"] == "M11240.1"


def _analyse_current(win):
    win._run_analysis()
    assert _spin(lambda: win.state.get("analysis_state") == "current"), "analysis did not complete"


# ---------- real engine round trips for the tagged operations ----------
def _stub_wsl(monkeypatch):
    import engine
    monkeypatch.setattr(engine.wsl, "annotate", lambda *a, **k: {
        "ok": True, "hits": [{"class_family": "LTR/Copia", "family": "Copia_LTR", "q_start": 10, "q_end": 200,
                              "strand": "+", "divergence": 1.0, "score": 500}],
        "dfam_library": {"version": "fixture"}, "repeatmasker_version": "fixture", "library_kind": "curated",
        "species": "Drosophila melanogaster"})
    monkeypatch.setattr(engine.wsl, "splice_align", lambda *a, **k: {
        "ok": True, "counts": {"exons": 1, "introns": 0}, "exons": [], "introns": [], "strand": "+",
        "minimap2_version": "fixture"})


def test_family_annotation_of_the_specimen_renders_and_carries_the_record_id(win, monkeypatch):
    _stub_wsl(monkeypatch)
    a, b = _two_specimens()
    _analyse(win, a + "\n" + b)
    win._show_record(1)
    win._annotate()
    assert _spin(lambda: win.state.get("family")), "a current annotation result was dropped"
    m = win.state["prov_manifest"]
    assert m["runType"] == "annotate" and m["input"]["id"] == win.state["records"][1]["id"]


def test_family_annotation_of_a_pasted_sequence_renders(win, monkeypatch):
    from teagle_core import sequtil
    _stub_wsl(monkeypatch)
    a, b = _two_specimens()
    _analyse(win, a)
    win.wslSource.setCurrentIndex(1)
    win.wslPaste.setPlainText(sequtil.parse_fasta(b)[0][1])
    win._annotate()
    assert _spin(lambda: win.state.get("family")), "a current annotation result was dropped"


def test_splice_alignment_renders_and_carries_the_record_id(win, monkeypatch):
    _stub_wsl(monkeypatch)
    a, _ = _two_specimens()
    _analyse(win, a)
    win.spliceTx.setPlainText(">t\n" + "ACGT" * 30)
    win._splice()
    assert _spin(lambda: win.state.get("splice")), "a current splice result was dropped"
    m = win.state["prov_manifest"]
    assert m["runType"] == "splice" and m["input"]["id"] == win.state["records"][0]["id"]


def test_a_result_whose_provenance_names_a_different_sequence_is_refused(win):
    a, _ = _two_specimens()
    _analyse(win, a)
    calls = _capture(win)
    win._annotate()
    win._on_done(calls[-1][2], {**_FAMILY, "provenance": {"runType": "annotate", "input": {"sha256": "0" * 64}}})
    assert not win.state.get("family")
    assert "different sequence" in " ".join(l.text() for l in win._notif.findChildren(QLabel))
    assert win.annotateBtn.text().upper() != "◴ ANNOTATING…"
    assert win.state["prov_manifest"]["runType"] == "analysis"


# ---------- a re-run on the same specimen is not a specimen change ----------
def test_rerunning_analysis_keeps_a_running_family_annotation(win):
    a, _ = _two_specimens()
    _analyse(win, a)
    calls = _capture(win)
    win._annotate()
    key = calls[-1][2]
    _analyse_current(win)                                             # same box, analysis run again
    assert win.annotateBtn.text().upper() == "◴ ANNOTATING…" and win.state.get("family_seq")
    win._on_done(key, _FAMILY)
    assert win.state.get("family") and win.annotateBtn.text().upper() == "RUN FAMILY ANNOTATION"
    assert "Copia_LTR" in _label_texts(win.wslBody.parentWidget())


def test_family_annotation_started_before_any_analysis_survives_the_first_analysis(win):
    a, _ = _two_specimens()
    win.seq.setPlainText(a)
    calls = _capture(win)
    win._annotate()
    key = calls[-1][2]
    _analyse_current(win)
    win._on_done(key, _FAMILY)
    assert win.state.get("family")


def test_rerunning_analysis_still_retires_a_running_primer_design(win):
    a, _ = _two_specimens()
    _analyse(win, a)
    calls = _capture(win)
    win._design()
    key = calls[-1][2]
    _analyse_current(win)
    assert not win._design_inflight
    win._on_done(key, {"candidates": [{"id": "P1", "left_seq": "A" * 20, "right_seq": "C" * 20, "product_size": 200,
                                       "left_tm": 60, "right_tm": 60, "left_gc": 50, "right_gc": 50, "penalty": 0.1}]})
    assert not win.state.get("candidates")


# ---------- deferred startup is tied to the window's lifetime ----------
def test_a_window_deleted_before_its_deferred_startup_never_runs_it(monkeypatch):
    import gc
    import main
    from PySide6.QtCore import QEvent
    ran = []
    monkeypatch.setattr(main.MainWindow, "_startup", lambda self: ran.append(self))
    w = main.MainWindow()
    w.deleteLater()
    QApplication.sendPostedEvents(None, QEvent.DeferredDelete)
    del w
    gc.collect()
    _app.processEvents()
    assert not ran


# ---------- the fault banner names the step and says what to do ----------
@pytest.mark.parametrize("key,words", [("wsl_status", "backend status"), ("env", "environment check"),
                                       ("pcr#1#0", "in-silico pcr"), ("add_custom_assembly", "assembly"),
                                       ("annotate_budget", "annotation"), ("some_new_step", "some new step")])
def test_a_fault_banner_names_the_step_and_what_to_do(win, key, words):
    win._on_failed(key, "KeyError: 'x'", "trace")
    text = " ".join(l.text() for l in win._notif.findChildren(QLabel))
    assert words in text.lower() and "this step" not in text.lower()
    assert "try it again" in text.lower() and "restart" in text.lower() and "KeyError: 'x'" in text
