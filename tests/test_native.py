"""Native (PySide6) app regression tests. Run headless via the offscreen Qt platform. Cover the
figure builders, the engine worker's three-tier taxonomy, and the analyze -> design -> in-silico PCR
workflow driven through the same in-process engine the GUI uses. Skipped if PySide6 is unavailable."""
import os, sys, time

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
    res = {"ok": True, "runType": "coordinate", "fasta": ">x\n" + "A" * 100 + "\n", "organism": "Homo sapiens",
           "assemblyName": "GRCh38.p14", "displayLocus": "chr13:1-100",
           "regions": [{"chrAccession": "NC_000013.11", "start": 1, "stop": 100, "strand": 2, "chromLabel": "chr13"}],
           "source": {"displayLocus": "chr13:1-100", "accession": "NC_000013.11", "runType": "coordinate"}}
    win._on_fetch(res)
    assert "chr13" in win.coordMeta.text() and win.accMeta.text() == ""
    assert win.state["source"].get("displayLocus") == "chr13:1-100"


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
    # #11 happy path: a panel-01 design->stage->PCR must still call the intended amplicon on-target
    assert any(a.get("on_target") for a in win.state["lastPcr"]["amplicons"]), "template-tagged pair lost its on-target call"
    assert win.card_prov.bodylay.count() >= 1        # provenance rendered


def test_staleness_guard_blocks_stale_primer(win):
    win._load_sample()
    win._run_analysis()
    assert _spin(lambda: win.state.get("analyzed_clean"))
    win.seq.setPlainText("ACGTACGTACGTACGTACGTACGTACGTACGT")   # change sequence after analysis
    assert win._stale_block() is True                          # design must be blocked


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


def test_on_genome_pcr_renders_sealed_advisory(win):
    win._genome_inflight = True
    win._on_genome_pcr({"ok": True, "organism": "Saccharomyces cerevisiae", "assemblyName": "R64",
        "assemblyAccession": "GCF_000146045.2", "n_seqs": 17,
        "advisory_note": "candidate priming sites — not wet-lab-validated amplicons",
        "amplicons": [{"source": "NC_001133.9", "start": 1001, "end": 1240, "length": 240, "strand": "+",
                       "single_primer": False, "on_target": False, "pair": "pair"}],
        "provenance": {"manifestSha256": "de" * 32,
                       "databases": [{"name": "RefSeq assembly R64", "version": "GCF_000146045.2", "sha256": "abc123"}]}})
    lp = win.state["lastPcr"]
    assert len(lp["amplicons"]) == 1 and lp["lanes"][0]["advisory"] is True   # advisory lane -> no "no on-target" stamp
    assert win._genome_inflight is False


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


def test_banner_level_distinguishes_success_from_error(win):
    # a success/advisory message must not render in the red error style with a warning triangle
    win._banner("Genome ready", level="success")
    assert win.errbanner.property("level") == "success" and not win.errbanner.isHidden()
    assert win.errbanner.text().startswith("✓")
    win._banner("bad input")                                   # default level = error
    assert win.errbanner.property("level") == "error" and win.errbanner.text().startswith("⚠")


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


def test_flank_region_clickable_and_slices_analyzed_seq(win):
    # a 5' flank region from the gene-model viewer must offer copy/design and slice the analyzed snapshot
    win.state["analyzed_seq"] = ">x\n" + "ACGTACGTAC" + "N" * 190       # first 10 bases are ACGTACGTAC
    items = win._region_menu({"start": 0, "end": 10, "label": "5prime_flank", "strand": "+"})
    labels = [lbl for lbl, _ in items]
    assert any("Copy FASTA" in l for l in labels) and any("primer" in l.lower() for l in labels)
    assert win._slice(0, 10) == "ACGTACGTAC"                             # flank coords index analyzed_seq


def test_app_theme_propagates_to_genome_viewers(win):
    import widgets
    win._load_sample(); win._run_analysis()
    assert _spin(lambda: win.state.get("analyzed_clean"))
    gvs = win.findChildren(widgets.GenomePanel)
    assert gvs and gvs[0].theme == "dark" and win.theme == "dark"
    gvs[0].view = {"start": 100.0, "end": 500.0}; saved = dict(gvs[0].view)
    win._toggle_theme()                                        # app -> light
    assert win.theme == "light" and gvs[0].theme == "white"    # viewer follows
    assert gvs[0].view == saved                                # pan/zoom preserved (in-place re-render)
    gvs[0]._set_theme("dark", user=True)                       # per-viewer MANUAL override (button click)
    win._toggle_theme(); win._toggle_theme()                   # two app toggles must NOT stomp a manual pick (H1)
    assert win.theme == "light" and gvs[0].theme == "dark"     # user's choice is kept across app-theme changes


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


def test_pcr_omits_target_span_for_mismatched_template(win):
    # a pair designed on a different template than the PCR sequence must NOT assert a (false) on-target
    seq = ">x\n" + "ACGT" * 60
    win.state["seq"] = seq
    win.state["analyzed_clean"] = win._clean_seq(seq)          # pass the stale-gate
    win.state["pcrPairs"] = [{"left_seq": "ACGTACGT", "right_seq": "TTTTGGGG",
                              "left_pos": [0, 8], "right_pos": [30, 38], "_tmpl_sig": "SOME-OTHER-TEMPLATE"}]
    bodies = []
    win.engine.submit = lambda op, body=None, key=None: bodies.append(body)
    win._run_pcr()
    assert bodies and bodies[0]["target_span"] is None         # mismatched template -> no target_span (off-target only)


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
    win._on_fetch({"ok": True, "accession": "M11240", "organism": "Drosophila melanogaster",
                   "title": "copia", "length": 5146, "fasta": ">M11240.1\n" + "ACGT" * 30})
    assert win.state["source"] and win.state["source"].get("accession") == "M11240"
    win.seq.setPlainText("ACGTACGTACGTACGTACGTACGT")           # genuine user edit
    assert win.state["source"] is None and win.accMeta.text() == ""


def test_programmatic_load_keeps_source(win):
    # the _loading guard: a fetch's own setPlainText must NOT null the source it is about to set
    win._on_fetch({"ok": True, "accession": "X05424", "organism": "Nicotiana", "title": "t",
                   "length": 100, "fasta": ">X05424.1\n" + "ACGT" * 25})
    assert win.state["source"] and win.state["source"].get("accession") == "X05424"


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
