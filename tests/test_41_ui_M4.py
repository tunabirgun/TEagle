"""4.1 provenance and reporting lane: structured provenance view, run list, copyable hashes, session export,
classification table, notebook entry and assay report. Drives the real MainWindow offscreen; nothing leaves the process."""
import copy, json, os, re, sys, time

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
import pytest

pytest.importorskip("PySide6")
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_APP = os.environ.get("M4_APP_ROOT") or os.path.join(_ROOT, "app")      # a pre-change copy can be pointed at for failing-before runs
sys.path.insert(0, os.path.join(_APP, "native"))
sys.path.insert(0, os.path.join(_APP, "backend"))
sys.path.insert(0, os.path.join(_ROOT, "tests"))

from PySide6.QtCore import QEvent, QEventLoop, QTimer
from PySide6.QtWidgets import QApplication, QFileDialog, QLabel, QPlainTextEdit, QPushButton

_app = QApplication.instance() or QApplication([])


class FakeSettings:
    store = {"quickstart_hidden": True}

    def __init__(self, *a):
        pass

    def value(self, key, default=None, type=None):
        return self.store.get(key, default)

    def setValue(self, key, value):
        self.store[key] = value


def _spin(pred, timeout=40):
    loop, t, t0 = QEventLoop(), QTimer(), time.time()
    t.setInterval(20)
    t.timeout.connect(lambda: (pred() or time.time() - t0 > timeout) and loop.quit())
    t.start(); loop.exec(); t.stop()
    return pred()


@pytest.fixture
def win(monkeypatch):
    import gc, fonts, main
    fonts.load_fonts()
    monkeypatch.setattr(main, "QSettings", FakeSettings)
    w = main.MainWindow()
    _app.processEvents()
    yield w
    w.engine.wait(5000)
    w.close(); w.deleteLater()
    QApplication.sendPostedEvents(None, QEvent.DeferredDelete)
    del w
    gc.collect()


def _analyse(win, text):
    win._set_seq(text)
    win._run_analysis()
    assert _spin(lambda: win.state.get("analysis_state") == "current"), "analysis did not complete"


def _example(name="M11240"):
    from teagle_core import examples
    return examples.load(name)


def _design_and_pcr(win):
    win._design()
    assert _spin(lambda: win.state.get("candidates")), "primer design did not return"
    win._pcr_stage_all()
    win._run_pcr()
    assert _spin(lambda: win.state.get("lastPcr")), "in-silico PCR did not return"


def _labels(root):
    return [l.text() for l in root.findChildren(QLabel)]


def _buttons(root):
    return [b for b in root.findChildren(QPushButton)]


def _save_to(monkeypatch, path):
    monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *a, **k: (str(path), ""))


@pytest.fixture(scope="module")
def full(request):
    """One analysis, design and PCR on the sample record, built once and reused read-only."""
    import gc, fonts, main
    fonts.load_fonts()
    mp = pytest.MonkeyPatch()
    mp.setattr(main, "QSettings", FakeSettings)
    from sample import make_sample
    w = main.MainWindow()
    w.setProperty("teagle_keep", True)                 # shared across this module; conftest must not sweep it
    _app.processEvents()
    w._set_seq(make_sample()); w._run_analysis()
    assert _spin(lambda: w.state.get("analysis_state") == "current")
    _design_and_pcr(w)
    yield w
    w.engine.wait(5000)
    w.close(); w.deleteLater()
    QApplication.sendPostedEvents(None, QEvent.DeferredDelete)
    mp.undo()
    gc.collect()


# ---------------------------------------------------------------- report.py: sealed parameters
def test_every_user_facing_sealed_key_has_a_label():
    import engine, report
    from teagle_core import primers
    keys = set(engine._detector_parameters()) | set(primers.DESIGN_DEFAULTS)
    user = {k for k in keys if not k.startswith(("structural_", "classify_"))}
    assert user <= set(report.PARAM_DISPLAY), sorted(user - set(report.PARAM_DISPLAY))


def test_param_sections_show_every_key_once_with_units():
    import engine, report
    m = engine.analyze(_example())["provenance"]
    sections = report.param_sections(m)
    shown = [row["key"] for _g, rows in sections for row in rows]
    assert sorted(shown) == sorted(m["parameters"]) and len(shown) == len(set(shown))
    rows = {row["key"]: row for _g, rs in sections for row in rs}
    assert rows["pbs_min_identity"]["text"] == "55 %" and "0.55" in rows["pbs_min_identity"]["tip"]
    assert rows["ltr_min_identity_pct"]["text"] == "80 %"
    assert rows["orf_min_aa"]["text"] == "40 aa" and rows["ltr_seed_window"]["text"] == "6,000 bp"
    groups = dict(sections)
    assert all(r["key"].startswith(("structural_", "classify_")) for r in groups[report.INTERNALS])
    extra = copy.deepcopy(m); extra["parameters"]["zz_new"] = 3
    got = dict(report.param_sections(extra))
    assert any(r["key"] == "zz_new" and r["label"] == "zz_new" for r in got[report.OTHER])


def test_group_hash_keeps_all_64_characters():
    import report
    h = "0123456789abcdef" * 4
    g = report.group_hash(h)
    assert g.replace(" ", "") == h and all(len(b) == 8 for b in g.split())


def test_panel_profiles_match_the_bundled_panel_and_refuse_a_different_file():
    import engine, report
    from teagle_core import domains
    m = engine.analyze(_example())["provenance"]
    profiles, note = report.panel_profiles(m)
    with open(domains.HMM_PATH, "rb") as fh:
        n_models = fh.read().count(b"HMMER3")
    assert note is None and len(profiles) == len(domains.DOMAIN_INFO) == n_models
    assert all(re.fullmatch(r"PF\d{5}\.\d+", p["accession"]) for p in profiles)
    assert all(p["accession"].startswith(domains.DOMAIN_INFO[p["name"]][3] + ".") for p in profiles)
    assert any(p["name"] == "RVT_1" and p["code"] == "RT" for p in profiles)
    bad = copy.deepcopy(m)
    next(d for d in bad["databases"] if "Pfam" in d["name"])["sha256"] = "0" * 64
    none, why = report.panel_profiles(bad)
    assert none == [] and "differs from the one this run sealed" in why


# ---------------------------------------------------------------- provenance card
def test_provenance_card_is_sectioned_with_full_values(full):
    import report
    texts = " ".join(_labels(full.card_prov))
    for heading in ("INPUT", "RUN AND SEAL", "SOFTWARE", "PARAMETERS", "ENVIRONMENT"):
        assert heading in texts, heading
    assert "Mismatches allowed per primer" in texts and "Exact 3′ bases required" in texts
    assert "…" not in texts, "no value on the card is cut with an ellipsis"
    m = full.state["prov_manifest"]
    assert report.group_hash(m["manifestSha256"]) in texts and report.group_hash(m["input"]["sha256"]) in texts
    env = m["environment"]
    assert str(env["os"]) in texts and str(env["python"]) in texts


def test_analysis_run_shows_labelled_parameters_and_pfam_accessions(win):
    _analyse(win, _example())
    texts = " ".join(_labels(win.card_prov))
    for expected in ("Minimum ORF length", "40 aa", "55 %", "Seed search window", "Minimum LTR–LTR identity"):
        assert expected in texts, expected
    assert "ltr_min_identity_pct" not in texts.replace("manifest key", "")        # raw keys live in tooltips
    assert "PF00078.33" in " ".join(_labels(win.card_prov)) + " ".join(
        l.text() for l in win.card_prov.findChildren(QLabel))
    json_boxes = [b for b in win.card_prov.findChildren(QPlainTextEdit)]
    assert all(b.isReadOnly() for b in json_boxes)


def test_run_list_keeps_the_analysis_after_primers_and_pcr(win):
    from PySide6.QtWidgets import QComboBox
    _analyse(win, _example())
    analysis_seal = win.state["analysis"]["provenance"]["manifestSha256"]
    _design_and_pcr(win)
    combos = win.card_prov.findChildren(QComboBox)
    runs = next(c for c in combos if "run" in c.accessibleName().lower() and "pair" not in c.accessibleName().lower())
    kinds = [runs.itemText(i).split(" · ")[1] for i in range(runs.count())]
    assert kinds == ["Analysis", "Primer design", "In-silico PCR"]
    latest = win.state["prov_manifest"]
    runs.setCurrentIndex(0)
    import report
    assert report.group_hash(analysis_seal) in " ".join(_labels(win.card_prov))
    assert "DATABASES" in " ".join(_labels(win.card_prov))
    assert win.state["prov_manifest"] is latest, "choosing a run to read must not change the latest run"
    assert any("ANALYSIS RESULT" in b.text().upper() for b in _buttons(win.card_prov))


def test_pcr_batch_gets_a_pair_chooser_and_exports_every_pair(win, tmp_path, monkeypatch):
    from PySide6.QtWidgets import QComboBox
    _analyse(win, _example())
    _design_and_pcr(win)
    n_pairs = len(win.state["prov_manifest_all"])
    assert n_pairs > 1
    pair = [c for c in win.card_prov.findChildren(QComboBox) if c.accessibleName().lower().startswith("pair")]
    assert pair and pair[0].count() == n_pairs
    out = tmp_path / "run.json"
    _save_to(monkeypatch, out)
    win._export_selected_run()
    assert len(json.loads(out.read_text(encoding="utf-8"))["runs"]) == n_pairs


def test_export_all_runs_reverifies_each_seal(win, tmp_path, monkeypatch):
    _analyse(win, _example())
    _design_and_pcr(win)
    expected = sum(len(op["manifests"]) for op in win._prov_list())
    out = tmp_path / "all.json"
    _save_to(monkeypatch, out)
    win._export_all_runs()
    doc = json.loads(out.read_text(encoding="utf-8"))
    assert len(doc["runs"]) == expected
    assert {r["runType"] for r in doc["runs"]} == {"analysis", "primer", "in-silico-pcr"}
    # a tampered parameter breaks its seal: nothing is written
    out2 = tmp_path / "tampered.json"
    _save_to(monkeypatch, out2)
    win._prov_list()[1]["manifests"][0]["parameters"]["opt_tm"] += 1
    win._export_all_runs()
    assert not out2.exists() and not (tmp_path / "tampered.json.part").exists()


def test_analysis_result_still_exports_with_a_later_run_selected(win, tmp_path, monkeypatch):
    import engine
    _analyse(win, _example())
    _design_and_pcr(win)
    out = tmp_path / "analysis.json"
    _save_to(monkeypatch, out)
    win._export_analysis_result()
    assert engine.validate_analysis_result_document(json.loads(out.read_text(encoding="utf-8")))


def test_every_hash_has_a_copy_button_that_copies_all_64_characters(full):
    m = full.state["prov_manifest"]
    names = {b.accessibleName(): b for b in _buttons(full.card_prov) if b.accessibleName().startswith("Copy ")}
    assert "Copy Sequence SHA-256" in names and "Copy Seal (manifest SHA-256)" in names
    names["Copy Sequence SHA-256"].click()
    got = QApplication.clipboard().text()
    assert got == m["input"]["sha256"] and len(got) == 64
    names["Copy Seal (manifest SHA-256)"].click()
    assert QApplication.clipboard().text() == m["manifestSha256"]


def test_invalidation_clears_the_run_list_and_every_export(win):
    _analyse(win, _example())
    _design_and_pcr(win)
    assert win._prov_list()
    win.seq.setPlainText("ACGT" * 30)                          # an edit: the analysis, and everything built on it, is stale
    assert not win._prov_list()
    labels = [b.text().upper() for b in _buttons(win.card_prov)]
    assert not [t for t in labels if "EXPORT" in t or "COPY NOTEBOOK" in t or "ASSAY REPORT" in t], labels


def test_a_second_analysis_starts_a_fresh_run_list(win):
    _analyse(win, _example("M11240"))
    _design_and_pcr(win)
    assert len(win._prov_list()) == 3
    _analyse(win, _example("X01005"))
    assert [op["manifests"][0]["runType"] for op in win._prov_list()] == ["analysis"]


# ---------------------------------------------------------------- classification table
def test_classification_tsv_matches_engine_results():
    import engine, main, report
    from teagle_core import examples
    from sample import make_sample
    res = [engine.analyze(examples.load(n)) for n in ("M11240", "X01005", "M80343", "X05424")] + [engine.analyze(make_sample())]
    for r in res:
        text = report.classification_tsv(r["records"], main._rank_view, main._call_view)
        lines = text.split("\n")
        assert lines[0].split("\t") == report.CLS_TSV_HEADER and text.endswith("\n") and "\r" not in text
        assert lines[0].startswith("#TE\tOrder\tSuperfamily\tClade\tComplete\tStrand\tDomains\t")
        rows = [l.split("\t") for l in lines[1:] if l]
        assert len(rows) == len(r["records"]) and all(len(x) == len(report.CLS_TSV_HEADER) for x in rows)
        col = {h: i for i, h in enumerate(report.CLS_TSV_HEADER)}
        for x, rec, bound in zip(rows, r["records"], r["provenance"]["input"]["records"]):
            assert x[col["Sequence_SHA256"]] == bound["sha256"] and x[col["#TE"]] == rec["id"]
            present = set((rec["classification"].get("completeness") or {}).get("present") or [])
            if x[col["Superfamily"]] in ("Copia", "Gypsy"):
                assert (x[col["Complete"]] == "yes") == ({"GAG", "PR", "INT", "RT", "RNaseH"} <= present)
    copia = report.classification_tsv(res[0]["records"], main._rank_view, main._call_view).split("\n")[1].split("\t")
    assert copia[1:3] == ["LTR", "Copia"] and copia[8] == "RLC" and copia[7] == "LTR/Copia" and copia[3] == "NA"
    tc1 = report.classification_tsv(res[1]["records"], main._rank_view, main._call_view).split("\n")[1].split("\t")
    assert tc1[2] == "Tc1_Mariner" and tc1[8] == "DTT"                  # TEsorter's REXdb spelling, not the display name
    cacta = main._rank_view({"te_class": "DNA/CACTA", "class": "Class II"})
    assert report._TESORTER_SUPERFAMILY[cacta["superfamily"]] == "EnSpm_CACTA"


def test_classification_tsv_withheld_invalid_and_failed():
    import main, report
    base = {"valid": True, "seq": "ACGT" * 30, "composition": {"length": 120},
            "classification": {"te_class": "LTR/Copia", "outcome": "withheld", "confidence": "Candidate",
                               "completeness": {"tier": "t", "present": []}, "n_domains": 0},
            "processing": {"classification": {"status": "withheld"}}, "domains": []}
    col = {h: i for i, h in enumerate(report.CLS_TSV_HEADER)}
    row = report.cls_row({**base, "id": "r1"}, main._rank_view, main._call_view)
    assert row[col["Wicker_code"]] == "XXX" and row[col["Order"]] == "unknown" and row[col["TEagle_class"]] == "LTR/Copia"
    none = copy.deepcopy(base); none["classification"].update(te_class="none", outcome="complete")
    assert report.cls_row({**none, "id": "r2"}, main._rank_view, main._call_view)[col["Wicker_code"]] == "NA"
    bad = report.cls_row({"id": "r3", "valid": False, "seq": "XX"}, main._rank_view, main._call_view)
    assert bad[col["Outcome"]] == "invalid_input" and bad[col["Wicker_code"]] == "NA" and len(bad) == len(col)
    failed = copy.deepcopy(base); failed["processing"]["classification"]["status"] = "failed"
    assert report.cls_row({**failed, "id": "r4"}, main._rank_view, main._call_view)[col["Outcome"]] == "failed"


def test_classification_tsv_failed_domain_scan(monkeypatch):
    import engine, main, report
    from teagle_core import domains, examples
    import functools
    real = domains.scan_domains

    @functools.wraps(real)
    def boom(*a, **k):
        raise RuntimeError("boom")
    monkeypatch.setattr(domains, "scan_domains", boom)
    res = engine.analyze(examples.load("M11240"))
    col = {h: i for i, h in enumerate(report.CLS_TSV_HEADER)}
    row = report.classification_tsv(res["records"], main._rank_view, main._call_view).split("\n")[1].split("\t")
    assert row[col["Domains"]] == "unavailable" and row[col["Domain_evidence"]] == "unavailable_failed"
    assert row[col["Strand"]] == "NA"


def test_classification_tsv_button_writes_a_bom_free_lf_file(win, tmp_path, monkeypatch):
    _analyse(win, _example("M11240") + "\n" + _example("X01005"))
    out = tmp_path / "c.cls.tsv"
    _save_to(monkeypatch, out)
    btn = next(b for b in _buttons(win.card_prov) if "CLASSIFICATION" in b.text().upper())
    btn.click()
    raw = out.read_bytes()
    assert raw.startswith(b"#TE\t") and b"\r" not in raw
    assert len([l for l in raw.decode("utf-8").split("\n") if l and not l.startswith("#")]) == len(win.state["records"])


# ---------------------------------------------------------------- notebook entry and assay report
def _clip():
    return QApplication.clipboard().mimeData()


def test_notebook_entry_carries_call_pairs_conditions_and_seals(full):
    full._copy_notebook_entry()
    data = _clip()
    text, page = data.text(), data.html()
    rec = full.state["records"][full.state["record_index"]]
    pairs = full.state["pcrPairs"]
    assert rec["id"] in text and rec["classification"]["te_class"] in text and "0-based, half-open" in text
    for c in pairs:
        assert c["left_seq"] in text and c["right_seq"] in text
    prm = full.state["lastPrimers"]["provenance"]["parameters"]
    assert f"{prm['min_tm']:g}–{prm['max_tm']:g}" in text and f"{prm['prod_min']:g}–{prm['prod_max']:g}" in text
    for m in full.state["prov_manifest_all"]:
        assert m["manifestSha256"] in text
    assert full.state["analysis"]["provenance"]["manifestSha256"] in text
    from teagle_core import sequtil
    assert sequtil.sha256(rec["seq"]) in text
    assert "<html" in page.lower() and "<script" not in page.lower()
    for header in ("sample_LTR_element", "demo construct"):
        assert header not in text and header not in page, "the user's own header must not leave the app"


def test_notebook_pair_locations_and_overlaps_come_from_the_design(full):
    full._copy_notebook_entry()
    text = _clip().text()
    c = full.state["pcrPairs"][0]
    assert f"({c['left_pos'][0]}–{c['left_pos'][1]})" in text and f"({c['right_pos'][0]}–{c['right_pos'][1]})" in text
    assert f"amplicon {c['left_pos'][0]}–{c['right_pos'][1]}" in text
    assert f"product {c['product_size']:,} bp" in text
    qc = c["qc"]["worst"]
    assert f"structure QC {qc}" in text


def test_pairs_bind_to_the_record_by_sequence_hash(full):
    from teagle_core import sequtil
    rec = full.state["records"][0]
    sha = sequtil.sha256(rec["seq"])
    assert all(c["design_provenance"]["input"]["sha256"] == sha for c in full.state["candidates"])


def test_assay_report_pdf_is_metadata_clean_and_complete(full, tmp_path):
    from pypdf import PdfReader
    from teagle_core.artifact_metadata import media_metadata_violations
    out = tmp_path / "assay.pdf"
    assert full._write_assay_report(str(out), "pdf") is True
    assert media_metadata_violations(out, expected_format=".pdf", author=None) == []
    assert b"Tuna" not in out.read_bytes()
    assert not (tmp_path / "assay.pdf.part").exists()
    reader = PdfReader(str(out))
    text = "".join("".join(p.extract_text().split()) for p in reader.pages)     # glyph spacing in a PDF is not stable
    rec = full.state["records"][0]
    assert rec["id"] in text and full.state["pcrPairs"][0]["left_seq"] in text and "0-based,half-open" in text
    assert full.state["analysis"]["provenance"]["manifestSha256"] in text
    assert "sample_LTR_element" not in text and len(reader.pages) <= 3


def test_assay_report_html_is_self_contained(full, tmp_path):
    out = tmp_path / "assay.html"
    assert full._write_assay_report(str(out), "html") is True
    page = out.read_text(encoding="utf-8")
    assert page.count("<svg") == 2 and "<script" not in page and "<!--" not in page
    assert "generator" not in page.lower() and str(tmp_path) not in page and "http" not in page.replace("http://www.w3.org/2000/svg", "")
    assert "sample_LTR_element" not in page and "\r" not in page
    import xml.dom.minidom as md
    ids = re.findall(r'\bid="([^"]+)"', page)
    assert len(ids) == len(set(ids)), "inline figures must not share element ids"


def test_report_excludes_pairs_designed_on_another_record(win):
    _analyse(win, _example("M11240") + "\n" + _example("X01005"))
    win._show_record(0)
    _design_and_pcr(win)
    win._show_record(1)
    ctx = win._report_context()
    assert ctx["pairs"] == [] and ctx["excluded"] == len(win.state["candidates"]) > 0
    win._copy_notebook_entry()
    text = _clip().text()
    assert "designed on another record" in text and "No primer pair has been designed on this record." in text
    assert "record 2 of 2" in text


def test_report_fails_closed_after_invalidation(win):
    _analyse(win, _example())
    win._copy_notebook_entry()
    before = _clip().text()
    win.seq.setPlainText("ACGT" * 30)
    QApplication.clipboard().setText("sentinel")
    win._copy_notebook_entry()
    assert QApplication.clipboard().text() == "sentinel" and before != "sentinel"
    assert "No current analysis" in " ".join(l.text() for l in win._notif.findChildren(QLabel))


def test_withheld_call_still_renders_a_report(win, tmp_path):
    _analyse(win, ">x\n" + "ACGTTGCATGCCGTA" * 8)
    ctx = win._report_context()
    assert ctx is not None
    win._copy_notebook_entry()
    text = _clip().text()
    assert "Cautions:" in text and "Structural evidence in this one sequence" in text
    assert "No primer pair has been designed on this record." in text and "In-silico PCR has not been run" in text


def test_invalid_record_has_nothing_to_record(win):
    _analyse(win, ">ok\n" + "ACGT" * 40 + "\n>bad\nMKVLAAGIVGLLLAQPEF")
    win._show_record(1)
    assert win._report_context() is None
    win._copy_notebook_entry()
    assert "not analysed" in " ".join(l.text() for l in win._notif.findChildren(QLabel))


def test_report_buttons_exist_only_while_analysis_is_current(win):
    _analyse(win, _example())
    names = [b.text().upper() for b in _buttons(win.card_prov)]
    assert any("COPY NOTEBOOK ENTRY" in n for n in names) and any("ASSAY REPORT" in n for n in names)
    assert all(b.property("resultExport") for b in _buttons(win.card_prov) if "NOTEBOOK" in b.text().upper())


def test_report_text_builders_state_empty_cases():
    import report
    ctx = {"version": "9.9.9", "generated": "2026-01-01 00:00",
           "specimen": {"id": "record_x", "length": 100, "sha256": "ab" * 32},
           "call": {"headline": "No TE signature detected", "cautions": []}}
    text = report.notebook_text(ctx)
    assert text.startswith("TEagle 9.9.9") and "No primer pair has been designed on this record." in text
    assert "In-silico PCR has not been run for these pairs." in text and "Element map" not in text
    html_text = report.report_html(ctx, {"map": "<svg></svg>"})
    assert "Element map" in html_text and "<svg></svg>" in html_text and "<!--" not in html_text


# ---------------------------------------------------------------- closing under a running job
def test_close_asks_only_while_a_long_job_runs(win, monkeypatch):
    from PySide6.QtGui import QCloseEvent
    asked = []
    monkeypatch.setattr(win, "_confirm_close", lambda jobs: asked.append(jobs) or False)
    ev = QCloseEvent(); win.closeEvent(ev)
    assert ev.isAccepted() and not asked                              # nothing running: no question
    win._genome_inflight = True
    ev = QCloseEvent(); win.closeEvent(ev)
    assert not ev.isAccepted() and asked == [["genome scan"]]         # declined: the window stays
    monkeypatch.setattr(win, "_confirm_close", lambda jobs: True)
    ev = QCloseEvent(); win.closeEvent(ev)
    assert ev.isAccepted()
    win._genome_inflight = False


def test_a_window_that_is_not_shown_never_blocks_on_the_close_question(win):
    win._annot_inflight = True
    assert not win.isVisible() and win._confirm_close(win._running_jobs()) is True
    win._annot_inflight = False


# ---------------------------------------------------------------- optional export footer
def test_with_footer_extends_the_view_box_escapes_text_and_renders():
    import figures, report
    from PySide6.QtCore import QByteArray
    from PySide6.QtSvg import QSvgRenderer
    svg = figures.svg_gel({"lanes": [{"label": "P1", "amplicons": [{"length": 200}]}]}, "white")
    w, h = [float(x) for x in re.search(r'viewBox="0 0 ([\d.]+) ([\d.]+)"', svg).groups()]
    out = report.with_footer(svg, ["TEagle 9 · a<b & c", "second line"])
    w2, h2 = [float(x) for x in re.search(r'viewBox="0 0 ([\d.]+) ([\d.]+)"', out).groups()]
    assert w2 == w and h2 > h and "a&lt;b &amp; c" in out and f'height="{h2:g}"' in out
    assert QSvgRenderer(QByteArray(out.encode("utf-8"))).isValid()
    assert report.with_footer(svg, []) == svg


def test_figure_footer_lines_carry_prefixes_of_the_hashes():
    import report
    seq, panel, seal = "a" * 64, "b" * 64, "c" * 64
    got = report.figure_footer("map", "4.1.0", "record_x", seq, panel, [seal])
    assert seq[:12] in got[0] and panel[:12] in got[0] and seal[:12] in got[0] and seq[:13] not in got[0]
    assert "0-based, half-open" in got[1]
    gel = report.figure_footer("gel", "4.1.0", "record_x", seq, panel, [("P1", seal)])
    assert "P1 " + seal[:12] in gel[0] and "not a measured gel" in gel[1]


def test_a_coordinate_span_is_not_labelled_in_bp():
    import report
    _g, row = report.param_row("intended_product_span", [354, 821], "in-silico-pcr")
    assert row["text"] == "354–821 (0-based, half-open)"
