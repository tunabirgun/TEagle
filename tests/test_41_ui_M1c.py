"""4.1 UI lane M1c: primer design, in-silico PCR, whole-genome scan and genome-management panels.
Drives the real MainWindow offscreen; engine calls are captured, never sent to WSL or the network."""
import os, sys, time

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
import pytest

pytest.importorskip("PySide6")
_NATIVE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "app", "native")
sys.path.insert(0, _NATIVE)

from PySide6.QtWidgets import QApplication, QLabel, QPushButton, QDialog, QMessageBox
from PySide6.QtCore import QEventLoop, QTimer

_app = QApplication.instance() or QApplication([])

F = "TTCACCACCATGGAGAAGGC"
R = "GGCATGGACTGTGGTCATGAG"


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
    _app.processEvents()
    yield w
    w.engine.wait(5000)
    w.close(); w.deleteLater()
    QApplication.sendPostedEvents(None, QEvent.DeferredDelete)
    del w
    gc.collect()


def _capture(win, ops=("annotate", "splice", "primers", "pcr", "fetch", "fetch_coords", "genome_pcr")):
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


def _labels(root):
    return " ".join(l.text() for l in root.findChildren(QLabel))


def _example(name="M11240"):
    from teagle_core import examples
    return examples.load(name)


def _seq_of(fasta):
    return "".join(l for l in fasta.splitlines() if not l.startswith(">"))


def _coord_specimen(win, regions, assembly="GCF_000001405.40"):
    """Load a coordinate-fetched specimen: regions [(start, stop, strand, sequence)] on NC_000001.11."""
    from teagle_core import provenance
    acc = "NC_000001.11"
    raw = {"accession": acc, "assemblyAccession": assembly, "taxid": "9606",
           "regions": [{"chrAccession": acc, "start": s, "stop": e, "strand": st} for s, e, st, _ in regions],
           "coordSystem": "1-based-inclusive", "retrievalType": "coordinate",
           "source": "NCBI E-utilities (efetch seq_start/seq_stop)",
           "endpoint": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
           "retrievedUtc": "2026-08-21T00:00:00+00:00"}
    recs = [(f"{acc}:{s}-{e}" + ("(-)" if st == 2 else ""), q) for s, e, st, q in regions]
    src = provenance.normalize_source_context(raw, recs)
    res = {"ok": True, "runType": "coordinate", "fasta": "".join(f">{h}\n{q}\n" for h, q in recs),
           "organism": "Homo sapiens", "assemblyName": "GRCh38.p14",
           "regions": [{"chrAccession": acc, "start": s, "stop": e, "strand": st, "chromLabel": "chr1"}
                       for s, e, st, _ in regions],
           "source": raw, "sourceContext": src}
    win._on_fetch(res)
    return src


def _designed(win, calls, **_):
    """Run the captured primer-design request through the real engine and hand the result to the window."""
    import engine
    op, body, key = calls[-1]
    assert op == "primers"
    res = engine.run_primers(body)
    win._on_done(key, res)
    return res["candidates"]


def _scan_rows(*rows):
    raw = "\n".join(rows) + "\n"
    return lambda acc, query, **k: {"ok": True, "raw": raw, "isPcr_version": "33x2", "target": "g.2bit",
                                   "sha256": "s", "n_seqs": 2}


def _isp(src, lo, hi, name="pair"):
    return f">{src}:{lo}+{hi} {name} {hi - lo + 1}bp {F} {R}"


# ======================= spec_conformance-1 / genomepcr_gel-1: the locus the GUI sends =======================
def _human(monkeypatch):
    import main
    monkeypatch.setattr(main, "all_assemblies", lambda: {"Homo sapiens": {
        "assemblyAccession": "GCF_000001405.40", "assemblyName": "GRCh38.p14"}})


def _fetched_and_designed(win, monkeypatch, strand=1, start=1001):
    _human(monkeypatch)
    seq = _seq_of(_example("M11240"))[:2000]
    q = seq if strand == 1 else "".join({"A": "T", "C": "G", "G": "C", "T": "A"}.get(c, "N") for c in reversed(seq))
    _coord_specimen(win, [(start, start + 1999, strand, q)])
    _analyse(win, win.seq.toPlainText())
    calls = _capture(win)
    win._design()
    return _designed(win, calls)


def test_design_locus_carries_window_strand_and_designed_span(win, monkeypatch):
    cands = _fetched_and_designed(win, monkeypatch)
    cand = cands[0]
    loci = win._design_locus_for(cand, "Homo sapiens")
    assert isinstance(loci, list) and len(loci) == 1
    L = loci[0]
    assert (L["accession"], L["start"], L["end"], L["strand"]) == ("NC_000001.11", 1001, 3000, 1)
    assert L["intended_product_span"] == cand["intended_product_span"]
    assert L["intended_endpoint_tolerance_bp"] == cand["intended_endpoint_tolerance_bp"]


@pytest.mark.parametrize("strand", [1, 2])
def test_scan_labels_the_designed_product_on_target_and_an_upstream_hit_off_target(win, monkeypatch, strand):
    import engine
    cands = _fetched_and_designed(win, monkeypatch, strand=strand)
    cand = cands[0]
    a, b = cand["intended_product_span"]
    if strand == 1:
        lo, hi = 1001 + a, 1001 + b - 1
    else:
        lo, hi = 3000 - b + 1, 3000 - a
    monkeypatch.setattr(engine.wsl, "genome_scan", _scan_rows(_isp("NC_000001.11", lo, hi),
                                                              _isp("NC_000001.11", 400, 640)))
    calls = _capture(win)
    win._scan_genome(cand, org="Homo sapiens")
    op, body, _key = calls[-1]
    assert op == "genome_pcr"
    res = engine.run_genome_pcr({**body, "assemblyAccession": "GCF_000001405.40", "organism": "Homo sapiens"})
    on = [(x["start"], x["end"]) for x in res["amplicons"] if x.get("on_target")]
    assert on == [(lo, hi)] and res["has_locus"] is True and res["n_on_target"] == 1


def test_scan_sends_no_locus_for_a_pasted_specimen_or_another_assembly(win, monkeypatch):
    cands = _fetched_and_designed(win, monkeypatch)
    cand = cands[0]
    assert win._design_locus_for(cand, "Mus musculus") is None
    import main
    monkeypatch.setattr(main, "all_assemblies", lambda: {"Homo sapiens": {
        "assemblyAccession": "GCF_999999999.1", "assemblyName": "other"}})
    assert win._design_locus_for(cand, "Homo sapiens") is None
    unsealed = {k: v for k, v in cand.items() if k != "design_provenance"}
    _human(monkeypatch)
    assert win._design_locus_for(unsealed, "Homo sapiens") is None


def test_scan_uses_the_region_the_pair_was_designed_on_in_a_multi_region_fetch(win, monkeypatch):
    _human(monkeypatch)
    seq = _seq_of(_example("M11240"))
    _coord_specimen(win, [(5001, 6200, 1, seq[:1200]), (9001, 10200, 1, seq[1200:2400])])
    _analyse(win, win.seq.toPlainText())
    calls = _capture(win)
    win._design()
    cand = _designed(win, calls)[0]
    loci = win._design_locus_for(cand, "Homo sapiens")
    assert [(L["start"], L["end"]) for L in loci] == [(5001, 6200)]


# ======================= whole-genome scan: size cap, export identity, row cap, locus matching =======================
def _scan_reply(amps, **summary):
    """A genome_pcr reply of the shape engine.run_genome_pcr returns, built from isPcr-style rows."""
    from teagle_core import genomepcr
    summary = {"n_pair": len(amps), "n_on": 0, "n_off": 0, "n_single": 0, "has_locus": False,
               "verdict": "v", "per_source": [], **summary}
    return {"ok": True, "organism": "Homo sapiens", "assemblyName": "GRCh38.p14",
            "assemblyAccession": "GCF_000001405.40", "n_seqs": 25, "amplicons": amps, "summary": summary,
            "provenance": {"manifestSha256": "ab" * 32, "runType": "genome-scan",
                           "parameters": {"engine": "isPcr", "isPcr_version": "33x2", "min_perfect": 15,
                                          "min_good": 15, "max_size": 4000, "min_size": 0},
                           "databases": [{"name": "RefSeq genome assembly", "version": "GCF_000001405.40",
                                          "sha256": "cd" * 32}]}}


def _prime_scan(win, cand, org="Homo sapiens"):
    win._prepared_genomes = [{"accession": "GCF_000001405.40", "n_seqs": 25, "bytes": 1_000_000}]
    win._refresh_genome_dropdown()
    calls = _capture(win)
    win._scan_genome(cand, org=org)
    return calls


def test_scan_window_covers_a_designed_product_longer_than_4000_bp_and_says_so(win):
    long_pair = {"id": "P1", "left_seq": F, "right_seq": R, "product_size": 5500}
    calls = _prime_scan(win, long_pair)
    assert calls[-1][1]["params"]["prod_max"] >= 5500
    assert f"{int(5500 * 1.25):,} bp" in _labels(win.genomeBody.parentWidget())     # the busy text names the window in force
    win._genome_inflight = False
    win._pending_scan = None
    calls2 = _prime_scan(win, {"id": "P2", "left_seq": F, "right_seq": R, "product_size": 300})
    assert calls2[-1][1]["params"]["prod_max"] == 4000
    assert "4,000" in _labels(win.genomeBody.parentWidget())      # the busy text names the cap in force


def test_scan_result_discloses_the_size_cap_and_identifies_pair_organism_and_criteria(win):
    from widgets import DataTable
    cand = {"id": "D2·P3", "left_seq": F, "right_seq": R, "product_size": 300}
    _prime_scan(win, cand)
    amp = {"source": "NC_000001.11", "start": 1101, "end": 1400, "length": 300, "strand": "+",
           "single_primer": False, "on_target": False, "sourceIdentitySource": "validated_fetch"}
    win._on_genome_pcr(_scan_reply([amp]))
    root = win.genomeBody.parentWidget()
    assert "4,000" in _labels(root)
    t = next(t for t in root.findChildren(DataTable) if t.export_base.startswith("TEagle_offtarget_scan"))
    notes = " ".join(t.export_notes)
    for needle in ("D2·P3", F, R, "Homo sapiens", "GCF_000001405.40", "33x2", "4,000"):
        assert needle in notes, needle
    assert "ab" * 8 not in notes and "cd" * 8 not in notes        # checksums stay in the manifest, not the table export
    assert "GCF_000001405.40" in t.export_base and "D2" in t.export_base      # the proposed filename tells two scans apart
    heads = t.export_headers()
    assert "Start (1-based)" in heads and "End (1-based)" in heads
    assert t.rows_data()[0][heads.index("Start (1-based)")] == "1101"


def test_unmatched_locus_is_reported_not_called_off_target(win):
    cand = {"id": "P1", "left_seq": F, "right_seq": R, "product_size": 300}
    _prime_scan(win, cand)
    win._scan_ctx = {**getattr(win, "_scan_ctx", {}), "loci": [{"accession": "NC_000001.11"}]}
    amp = {"source": "CM000663.2", "start": 1101, "end": 1400, "length": 300, "strand": "+",
           "single_primer": False, "on_target": False, "sourceIdentitySource": "validated_fetch"}
    win._on_genome_pcr(_scan_reply([amp], has_locus=True, n_off=1,
                                   verdict="1 off-target candidate priming site(s); no product at the intended locus"))
    text = _labels(win.genomeBody.parentWidget())
    assert "NC_000001.11" in text and "could not be matched" in text


def test_a_scan_with_thousands_of_products_shows_a_capped_table_but_exports_every_row(win):
    from widgets import DataTable
    import main
    _prime_scan(win, {"id": "P1", "left_seq": F, "right_seq": R, "product_size": 300})
    amps = [{"source": f"NC_{i % 20:06d}.1", "start": 1000 + i, "end": 1300 + i, "length": 301, "strand": "+",
             "single_primer": False, "on_target": False, "sourceIdentitySource": "validated_fetch"}
            for i in range(1500)]
    t0 = time.time()
    win._on_genome_pcr(_scan_reply(amps, n_off=1500, has_locus=False))
    assert time.time() - t0 < 20
    root = win.genomeBody.parentWidget()
    t = next(t for t in root.findChildren(DataTable) if t.export_base.startswith("TEagle_offtarget_scan"))
    assert t.rowCount() == main.MainWindow._SCAN_ROWS_SHOWN < 1500
    assert len(t.rows_data()) == 1500
    assert "1,500" in _labels(root) and str(main.MainWindow._SCAN_ROWS_SHOWN) in _labels(root)


def test_gel_menu_never_offers_an_empty_fasta_and_names_the_locus(win):
    amp = {"source": "NC_000001.11", "start": 1101, "end": 1400, "length": 300, "strand": "-",
           "sourceIdentitySource": "validated_fetch"}
    region = {"amplicon": amp, "pair": "genome", "has_locus": False, "tip": "300 bp · 1 priming site"}
    labels = [l for l, _f in win._gel_menu(region)]
    assert not any(l.startswith("Copy FASTA") for l in labels)
    locus = next(l for l in labels if l.startswith("Copy locus"))
    assert "1-based" in locus
    fn = dict(win._gel_menu(region))[locus]
    fn()
    assert QApplication.clipboard().text() == "NC_000001.11:1101-1400(-)"
    local = {"amplicon": {"source": "rec (input)", "start": 10, "end": 210, "length": 200, "seq": "ACGT" * 50,
                          "matches_intended_product": True, "on_target": None},
             "pair": "D1·P2", "has_locus": False, "tip": "200 bp · 1 priming site"}
    labels = [l for l, _f in win._gel_menu(local)]
    assert any(l.startswith("Copy FASTA") for l in labels)
    assert any("0-based" in l for l in labels)
    merged = {**local, "tip": "200 bp · 3 priming sites · 3 products drawn as one band (within 3%)"}
    assert any("3" in l and "band" in l for l, _f in win._gel_menu(merged))


# ======================= primer pair identity, local PCR table, gel =======================
def _cand(i, f=F, r=R, size=200, **extra):
    return {"id": f"P{i}", "left_seq": f, "right_seq": r, "product_size": size, "left_tm": 60.0, "right_tm": 60.0,
            "left_gc": 50.0, "right_gc": 50.0, "penalty": 0.1, "left_pos": [10, 30], "right_pos": [size - 20, size],
            **extra}


def _tables(root):
    from widgets import DataTable
    return {t.export_base: t for t in root.findChildren(DataTable)}


def _run_pcr_with(win, results, pairs):
    """Stage `pairs`, press Run (the request is captured), then deliver one engine reply per pair."""
    calls = _capture(win)
    win.state["pcrPairs"] = list(pairs)
    win._run_pcr()
    assert len(calls) == len(pairs), "the run was refused"
    for (_op, _body, key), res in zip(calls, results):
        win._on_done(key, res)
    return calls


def _pcr_hit(start=10, end=210, **kw):
    return {"source": "rec (input)", "start": start, "end": end, "length": end - start, "product_length": end - start,
            "product_key": f"k:{start}:{end}", "mechanism_key": f"k:F+R:{start}:{end}", "fwd_mm": 0, "rev_mm": 0,
            "fwd_tail5": 0, "rev_tail5": 0, "fwd_anneal_len": 20, "rev_anneal_len": 21, "left_primer": "F",
            "right_primer": "R", "primer_orientation": "F+R", "single_primer": False,
            "matches_intended_product": True, "covers_primer3_target": None, "on_target": None,
            "seq": "ACGT" * ((end - start) // 4), **kw}


def test_each_design_gets_a_batch_tag_that_follows_a_pair_into_the_queue_lane_and_exports(win, monkeypatch):
    import main, widgets
    _analyse(win, _example())
    win._on_primers({"candidates": [_cand(1), _cand(2, f="ACGT" * 5), _cand(3, f="TTGA" * 5)]})
    first = win.state["candidates"]
    assert {c["_batch"] for c in first} == {"D1"}
    win._on_primers({"candidates": [_cand(1, f="GGGA" * 5), _cand(2, f="CCCA" * 5), _cand(3, f="AACC" * 5)]})
    second = win.state["candidates"]
    assert {c["_batch"] for c in second} == {"D2"}
    table = _tables(win.primBody.parentWidget())["TEagle_primers"]
    assert [row[0] for row in table.rows_data()] == ["D2·P1", "D2·P2", "D2·P3"]
    win._add_pcr_pair(second[2])
    win._add_pcr_pair(first[0])
    q = _labels(win.pcrQueueBox.parentWidget())
    assert "D2·P3" in q and "D1·P1" in q
    win._move_pair(0, 1)                                          # reordering must not rename a pair
    assert win.state["pcrPairs"][0]["_batch"] == "D1" and "D2·P3" in _labels(win.pcrQueueBox.parentWidget())
    hits = [{"amplicons": [_pcr_hit()], "status": "complete"}, {"amplicons": [_pcr_hit(30, 230)], "status": "complete"}]
    _run_pcr_with(win, hits, win.state["pcrPairs"])
    lanes = [l["label"] for l in win.state["lastPcr"]["lanes"]]
    assert lanes == ["D1·P1", "D2·P3"]
    t = _tables(win.pcrBody.parentWidget())["TEagle_amplicons"]
    assert [r[0] for r in t.rows_data()] == ["D1·P1", "D2·P3"]
    notes = " ".join(t.export_notes)
    assert "D2·P3" in notes and second[2]["left_seq"] in notes and second[2]["right_seq"] in notes
    saved = []
    monkeypatch.setattr(widgets, "save_fasta", lambda text, base, parent=None: saved.append(text))
    btn = next(b for b in win.pcrBody.parentWidget().findChildren(QPushButton) if "FASTA" in b.text().upper())
    btn.click()
    heads = [l for l in saved[0].splitlines() if l.startswith(">")]
    assert len(heads) == 2 and len(set(h.split()[0] for h in heads)) == 2
    assert "D2-P3" in saved[0] and second[2]["left_seq"] in saved[0]


def test_the_picker_in_the_genome_card_and_the_pair_fasta_use_the_tagged_name(win):
    win._on_primers({"candidates": [_cand(1), _cand(2, f="ACGT" * 5)]})
    assert win.scanPicker.itemText(1).startswith("D1·P2")
    pair_fasta = None
    t = _tables(win.primBody.parentWidget())["TEagle_primers"]
    menu = dict(t._row_menu(1))
    menu["Copy pair FASTA"]()
    pair_fasta = QApplication.clipboard().text()
    assert pair_fasta.startswith(">D1-P2_F\n") and "\n>D1-P2_R\n" in pair_fasta


def test_an_untargeted_design_never_shows_a_red_no_on_target_caption(win):
    import engine
    _analyse(win, _example())
    calls = _capture(win)
    win._design()
    cands = _designed(win, calls)
    win.state["pcrPairs"] = []
    win._add_pcr_pair(win.state["candidates"][0])
    calls.clear()
    win._run_pcr()
    op, body, key = calls[-1]
    assert op == "pcr"
    win._on_done(key, engine.run_pcr(body))
    lane = win.state["lastPcr"]["lanes"][0]
    assert any(a.get("matches_intended_product") for a in lane["amplicons"])
    panel = win.pcrBody.parentWidget().findChildren(__import__("widgets").FigurePanel)[0]
    svg = panel.build_fn("white")
    assert "no on-target" not in svg
    tip = " ".join(r["tip"] for r in __import__("figures").gel_regions({"lanes": [lane]}))
    assert "off-target" not in tip


# ======================= stale card-06 state after the specimen changes =======================
def test_changing_the_specimen_clears_the_scan_picker_and_a_rendered_scan_verdict(win):
    _prime_scan(win, {"id": "P1", "left_seq": F, "right_seq": R, "product_size": 300})
    amp = {"source": "NC_000001.11", "start": 1101, "end": 1400, "length": 300, "strand": "+",
           "single_primer": False, "on_target": False, "sourceIdentitySource": "validated_fetch"}
    win._on_genome_pcr(_scan_reply([amp], verdict="UNIQUE-VERDICT-TEXT"))
    win._on_primers({"candidates": [_cand(1), _cand(2, f="ACGT" * 5)]})
    assert win.scanPicker.count() == 2 and win.pcrStageAll.isEnabled()
    assert "UNIQUE-VERDICT-TEXT" in _labels(win.genomeBody.parentWidget())
    win._invalidate_analysis("The specimen changed.", cancel_request=True)
    assert win.scanPicker.count() == 1 and win.scanPicker.currentData() is None
    assert not win.scanBtn.isEnabled() and not win.pcrStageAll.isEnabled()
    assert "UNIQUE-VERDICT-TEXT" not in _labels(win.genomeBody.parentWidget())


# ======================= primer form, routed design, conditions =======================
def _set_fields(win, **kv):
    for k, v in kv.items():
        win.pfields[k].setText(str(v))


@pytest.mark.parametrize("fields,needle", [
    ({"pMaxS": "40", "pOptS": "20"}, "36"),                      # Primer3's own ceiling on primer length
    ({"pMin": "40", "pMaxS": "27"}, "twice"),                    # backend rule: prod_min >= 2 * max_size
    ({"pGcClamp": "30"}, "clamp"),
    ({"pNum": "100000"}, "50"),
    ({"pTm": "250", "pTmMin": "57", "pTmMax": "300"}, "Tm"),
])
def test_primer_form_rejects_values_the_designer_would_refuse_with_a_plain_message(win, fields, needle):
    _analyse(win, _example())
    calls = _capture(win)
    _set_fields(win, **fields)
    win._design()
    assert not calls, "an invalid form must not reach the engine"
    body = _labels(win.primBody.parentWidget())
    assert needle in body and "OSError" not in body and "built-in maximum" not in body


@pytest.mark.parametrize("preset", ["standard", "qpcr", "highspec", "permissive"])
def test_every_preset_passes_the_primer_form_rules(win, preset):
    win._apply_preset(preset)
    win._read_primer_params()


def test_routed_domain_longer_than_the_product_window_is_not_submitted_and_says_what_is_needed(win):
    _analyse(win, _example())
    calls = _capture(win)
    win._design_for_domain(100, 835, "RT domain")                 # 735 bp against the Standard 150-500 window
    assert not calls
    body = _labels(win.primBody.parentWidget())
    need = 735 + 2 * 18
    assert "735" in body and str(need) in body and "Prod max" in body
    assert win.designBtn.isEnabled()


def test_routed_design_states_it_flanks_and_spans_and_shows_primer_positions(win):
    win._on_primers({"candidates": [_cand(1)], "_routed": {"label": "RT domain", "length": 300, "routed": True}})
    text = _labels(win.primBody.parentWidget())
    assert "flank" in text and "RT domain" in text and "300" in text
    t = _tables(win.primBody.parentWidget())["TEagle_primers"]
    heads = t.export_headers()
    assert "Fwd pos (0-based, half-open)" in heads and "Rev pos (0-based, half-open)" in heads
    row = t.rows_data()[0]
    assert row[heads.index("Fwd pos (0-based, half-open)")] == "10–30"


def test_empty_design_names_the_domain_and_never_blames_the_parameters_for_the_template(win):
    win._on_primers({"candidates": [], "explain_pair": "considered 7001556, unacceptable product size 7001556, ok 0",
                     "_routed": {"label": "RT domain", "length": 300, "routed": True, "prod_max": 500, "min_size": 18},
                     "cautions": ["3 ambiguous template base(s) were treated as N for Primer3; no primer covers them."]})
    text = _labels(win.primBody.parentWidget())
    assert "RT domain" in text and "ambiguous template base" in text
    assert "Permissive preset or widen" not in text


def test_primer_conditions_are_read_from_the_result_and_exported(win):
    d = {"candidates": [_cand(1)], "design_tm_conditions": {"mv_conc": 50.0, "dv_conc": 1.5, "dntp_conc": 0.6, "dna_conc": 50.0},
         "oligoqc_conditions": {"mv_conc": 50.0, "dv_conc": 0.0, "dntp_conc": 0.0, "dna_conc": 250.0, "temp_c": 37.0}}
    win._on_primers(d)
    text = _labels(win.primBody.parentWidget())
    assert "1.5 mM Mg" in text and "37 °C" in text and "250 nM" in text
    t = _tables(win.primBody.parentWidget())["TEagle_primers"]
    assert "1.5 mM Mg" in " ".join(t.export_notes) and "250 nM" in " ".join(t.export_notes)
    win._on_primers({"candidates": [_cand(1)]})                   # a result that predates the fields says so, never guesses
    assert "not recorded" in _labels(win.primBody.parentWidget())


def test_positive_delta_g_is_shown_as_no_structure_not_as_a_value():
    import main
    txt, flag, tip, export = main._metric_cell([("", {"p3": 0.4, "vrna": None, "flag": "ok", "agree": "single"})])
    assert txt == "none" and export == "none" and "no stable structure" in tip.lower() and "0.4" in tip
    txt, _f, _t, export = main._metric_cell([("F", {"p3": -3.1, "vrna": -3.0, "flag": "caution", "agree": "agree"})])
    assert txt == "-3.1!" and export == "-3.1"


def test_structure_detail_states_the_qc_conditions_from_the_result(win, monkeypatch):
    metric = {"p3": -2.5, "vrna": None, "flag": "caution", "agree": "none"}
    cand = {"id": "P1", "_batch": "D1", "left_seq": "ACGTACGTACGT", "right_seq": "TGCATGCATGCA",
            "qc": {"ok": True, "left": {"hairpin": metric, "self_dimer": metric},
                   "right": {"hairpin": metric, "self_dimer": metric}, "hetero_dimer": metric, "end_stability": metric,
                   "conditions": {"temp_c": 37, "mv_conc": 50, "dv_conc": 0, "dntp_conc": 0, "dna_conc": 250},
                   "engines": {"primer3": "x", "viennarna": "unavailable"}}}
    monkeypatch.setattr(QDialog, "exec", lambda self: 0)
    win._structure_detail(cand)
    dlg = next(w for w in QApplication.topLevelWidgets() if isinstance(w, QDialog) and w.windowTitle().startswith("Secondary structure"))
    assert "D1·P1" in dlg.windowTitle()
    text = _labels(dlg)
    assert "0 mM Mg" in text and "250 nM" in text
    from widgets import DataTable
    cells = [dlg.findChild(DataTable).item(r, 4).text() for r in range(6)]
    assert set(cells) == {"no structure"}                         # engines 'none' reads as words, not as a bare code
    assert dlg.findChild(DataTable).export_base == "TEagle_D1-P1_secondary_structure"
    assert "Primer3 only" in text
    dlg.close()


# ======================= in-silico PCR window, lanes, tables =======================
def test_a_designed_product_outside_the_pcr_window_widens_it_visibly_instead_of_vanishing(win):
    _analyse(win, _example())
    calls = _capture(win)
    seen = []
    win._banner = lambda m, level="error": seen.append((m, level))
    win.state["pcrPairs"] = [_cand(1, size=1274), _cand(2, size=300)]
    win._run_pcr()
    prm = calls[-1][1]["params"]
    assert prm["prod_min"] <= 300 and prm["prod_max"] >= 1274
    assert win.pfields["pcrPmax"].text() == str(prm["prod_max"])
    said = " ".join(m for m, _l in seen).replace(",", "")
    assert "1274" in said and "70" in said and "1000" in said and any(l == "warn" for _m, l in seen)


def test_a_pair_inside_the_window_leaves_it_alone(win):
    _analyse(win, _example())
    calls = _capture(win)
    win.state["pcrPairs"] = [_cand(1, size=300)]
    win._run_pcr()
    assert calls[-1][1]["params"]["prod_min"] == 70 and calls[-1][1]["params"]["prod_max"] == 1000


def test_mismatch_tail_and_anneal_columns_name_the_primer_not_the_site(win):
    r_then_f = _pcr_hit(10, 210, left_primer="R", right_primer="F", primer_orientation="R+F",
                        fwd_mm=1, rev_mm=0, fwd_tail5=4, rev_tail5=0, fwd_anneal_len=18, rev_anneal_len=22,
                        matches_intended_product=False)
    f_then_f = _pcr_hit(300, 520, left_primer="F", right_primer="F", primer_orientation="F+F", single_primer=True,
                        fwd_mm=2, rev_mm=1, fwd_tail5=0, rev_tail5=0, fwd_anneal_len=20, rev_anneal_len=20,
                        matches_intended_product=False)
    amps = [{**a, "pair": "P1"} for a in (r_then_f, f_then_f)]
    win._render_pcr([{"label": "P1", "amplicons": amps, "has_locus": False}], amps,
                    [{"pair": "P1", "status": "complete", "error": None}])
    t = _tables(win.pcrBody.parentWidget())["TEagle_amplicons"]
    heads = t.export_headers()
    mm, tail, ann = (heads.index(next(h for h in heads if h.startswith(k))) for k in ("Mism", "Tail", "Anneal"))
    rows = t.rows_data()
    assert rows[0][mm] == "0/1"           # the R primer sits at the LEFT site; the F primer's mismatch is the right site's
    assert rows[0][tail] == "0/4" and rows[0][ann] == "22/18"
    assert rows[1][mm] == "2,1/—"         # a self-priming F+F product has two F sites and no R site
    assert "Coords (0-based, half-open)" in heads
    assert "annealed 3" in " ".join(t.export_notes) and "tail" in " ".join(t.export_notes).lower()


def test_a_failed_pair_is_named_with_its_error_in_the_result(win):
    _analyse(win, _example())
    win._on_primers({"candidates": [_cand(1), _cand(2, f="ACGT" * 5)]})
    win._add_pcr_pair(win.state["candidates"][0]); win._add_pcr_pair(win.state["candidates"][1])
    _run_pcr_with(win, [{"amplicons": [_pcr_hit()], "status": "complete"},
                        {"error": "ValueError: boom", "amplicons": []}], win.state["pcrPairs"])
    text = _labels(win.pcrBody.parentWidget())
    assert "D1·P2" in text and "boom" in text


def test_the_pcr_note_reports_the_widened_window(win):
    _analyse(win, _example())
    win.state["pcrPairs"] = [_cand(1, size=1274)]
    win._on_primers({"candidates": [_cand(1, size=1274)]})
    win.state["pcrPairs"] = [win.state["candidates"][0]]
    _run_pcr_with(win, [{"amplicons": [_pcr_hit(10, 1284)], "status": "complete"}], win.state["pcrPairs"])
    assert "widened" in _labels(win.pcrBody.parentWidget())


# ======================= genome manager and downloads =======================
_YEAST = {"accession": "GCF_000146045.2", "n_seqs": 17, "bytes": 12_000_000, "fasta_bytes": 12_500_000, "sha256": "a" * 64}


def _open_manager(win, genomes=(_YEAST,)):
    win._genome_mgr_open = True
    win._on_genome_list({"ok": True, "genomes": list(genomes)})
    assert win._genome_mgr is not None and win._genome_mgr.isVisible()
    return win._genome_mgr


def _row_buttons(win):
    tbl = win._genome_mgr_table
    return [tbl.cellWidget(r, tbl.columnCount() - 1).findChild(QPushButton) for r in range(tbl.rowCount())]


def _no_modal(monkeypatch, answer=QMessageBox.No):
    texts = []
    monkeypatch.setattr(QMessageBox, "exec", lambda self: texts.append(self.text() + " " + self.informativeText()) or answer)
    return texts


def test_an_open_manager_is_locked_when_a_scan_starts_and_unlocked_when_it_ends(win):
    mgr = _open_manager(win)
    ops = []
    win.engine.submit = lambda op, body=None, key=None: ops.append(op)
    win._prepared_genomes = [_YEAST]; win._refresh_genome_dropdown()
    win._scan_genome({"id": "P1", "left_seq": F, "right_seq": R}, org="Saccharomyces cerevisiae")
    assert "genome_list" in ops                                   # the open manager is told to rebuild
    win._on_genome_list({"ok": True, "genomes": [_YEAST]})
    assert not any(b.isEnabled() for b in _row_buttons(win))
    assert "whole-genome scan finishes" in _labels(win._genome_mgr)      # it names what is actually running
    win._on_genome_pcr({"ok": True, "organism": "x", "amplicons": [], "summary": {}, "provenance": {}})
    win._on_genome_list({"ok": True, "genomes": [_YEAST]})
    assert any(b.isEnabled() for b in _row_buttons(win))
    win._genome_mgr.close()


def test_a_failed_download_in_an_open_manager_re_enables_its_buttons(win):
    _open_manager(win)
    ops = []
    win.engine.submit = lambda op, body=None, key=None: ops.append(op)
    win._on_manager_download("Homo sapiens")
    assert ops.count("genome_list") >= 1 and "genome_prepare" in ops
    win._on_genome_list({"ok": True, "genomes": [_YEAST]})
    assert not any(b.isEnabled() for b in _row_buttons(win)) and "genome download finishes" in _labels(win._genome_mgr)
    ops.clear()
    win._on_genome_prepare({"ok": False, "error": "FAIL insufficient usable disk (5G free, need >=8G for a genome)"})
    assert "genome_list" in ops                                   # the failure branch re-lists, so the manager rebuilds
    win._on_genome_list({"ok": True, "genomes": [_YEAST]})
    assert any(b.isEnabled() for b in _row_buttons(win))
    win._genome_mgr.close()


def test_delete_is_refused_while_a_scan_download_or_annotation_runs(win, monkeypatch):
    texts = _no_modal(monkeypatch, QMessageBox.Yes)
    ops = []
    win.engine.submit = lambda op, body=None, key=None: ops.append(op)
    for flag in ("_genome_inflight", "_genome_prep_inflight", "_annot_inflight"):
        setattr(win, flag, True)
        win._remove_genome(_YEAST["accession"], "Saccharomyces cerevisiae")
        setattr(win, flag, False)
    assert "genome_remove" not in ops and not texts


def test_delete_warns_about_annotation_work_and_states_the_real_size(win, monkeypatch):
    texts = _no_modal(monkeypatch)
    win._prepared_genomes = [_YEAST]
    win._remove_genome(_YEAST["accession"], "Saccharomyces cerevisiae")
    t = texts[0]
    assert "12 MB" in t and "~1 GB" not in t
    assert "annotation" in t and "FASTA" in t and "if there is one" in t        # a cached FASTA exists only after an annotation
    assert "12.5" not in t and "13 MB" not in t                                  # the download-time FASTA size is not on disk any more


def test_download_wording_applies_the_disk_floor_to_every_genome_and_success_points_to_the_scan_button(win, monkeypatch):
    texts = _no_modal(monkeypatch)
    win._pending_scan = {"org": "Saccharomyces cerevisiae", "cand": {"id": "P1"}}
    win._genome_inflight = True
    win._on_genome_pcr({"ok": False, "need_prepare": True, "error": "not downloaded"})
    assert "every genome" in texts[0].lower() and "8 GB" in texts[0]
    win._on_genome_prepare({"ok": True, "organism": "Saccharomyces cerevisiae", "assemblyAccession": "GCF_000146045.2",
                            "n_seqs": 17, "bytes": 12_000_000})
    assert "Scan whole genome" in _labels(win.genomeBody.parentWidget())
    mgr = _open_manager(win)
    assert "every genome" in _labels(mgr).lower()
    mgr.close()


# ======================= whole-genome annotation =======================
_HUMAN = {"accession": "GCF_000001405.40", "n_seqs": 25, "bytes": 800_000_000, "fasta_bytes": 3_100_000_000, "sha256": "b" * 64}


def _annot_dialog(win, calls_ops=("annotate_budget", "genome_annotate", "genome_annotate_reset", "genome_annotate_log")):
    win._prepared_genomes = [_HUMAN]
    calls = []
    real = win.engine.submit

    def submit(op, body=None, key=None):
        if op in calls_ops:
            calls.append((op, body, key or op))
        else:
            real(op, body, key=key)
    win.engine.submit = submit
    win._annotate_genome_dialog()
    return calls


def _budget_ok(**kw):
    return {"ok": True, "cores": 8, "mem_gb": 16, "recommended_threads": 4, "limited_by": "memory", "avail_gb": 300,
            "disk_needed_gb": 8, "disk_ok": True, "library_mode": "dfam_curated", "budget_request_id": 1,
            "library_families_for_species": 1439, **kw}


def test_the_cost_probe_sends_the_measured_fasta_size(win):
    calls = _annot_dialog(win)
    op, body, _k = calls[-1]
    assert op == "annotate_budget" and body["fasta_bytes"] == _HUMAN["fasta_bytes"]
    win._annot_dialog.close()


def test_a_failed_budget_probe_shows_the_error_keeps_start_disabled_and_offers_a_retry(win):
    calls = _annot_dialog(win)
    _op, body, _k = calls[-1]
    win._on_annotate_budget({"ok": False, "error": "TimeoutExpired: wsl", "library_mode": body["library_mode"],
                             "budget_request_id": body["budget_request_id"]})
    cost, warn, start, _g, _s = win._annot_cost_widgets
    assert "None" not in cost.text() and "None" not in warn.text()
    assert "TimeoutExpired" in warn.text() and not warn.isHidden()
    assert not start.isEnabled()
    retry = win._annot_retry
    assert not retry.isHidden()
    n = len(calls)
    retry.click()
    assert len(calls) == n + 1 and calls[-1][0] == "annotate_budget"
    win._on_annotate_budget(_budget_ok(budget_request_id=calls[-1][1]["budget_request_id"]))
    assert start.isEnabled() and retry.isHidden()
    win._annot_dialog.close()


def test_a_running_annotation_shows_a_persistent_progress_line_and_hides_it_when_it_ends(win):
    calls = _capture(win, ops=("genome_annotate", "genome_annotate_log"))
    win._annot_cost_state = {"sha256": _HUMAN["sha256"]}
    win._start_genome_annotate(_HUMAN["accession"], "Homo sapiens", "default")
    bar = win._annotStatus
    assert not bar.isHidden() and "Homo sapiens" in bar.caption.text()
    win._on_genome_annotate_log({"log": "chunk 3 of 40 done"})
    assert "chunk 3 of 40" in bar.caption.text()
    assert win.card_genome.isAncestorOf(bar)                      # lives in card 06, outside the panels that get cleared
    win._on_genome_annotate({"ok": False, "error": "the annotation did not complete"})
    assert bar.isHidden()
    win._on_genome_annotate_log({"log": "late line"})
    assert bar.isHidden()                                         # a late log line never revives a finished run's bar


def test_a_crashed_annotation_also_clears_the_progress_line(win):
    _capture(win, ops=("genome_annotate", "genome_annotate_log"))
    win._start_genome_annotate(_HUMAN["accession"], "Homo sapiens", "default")
    bar = win._annotStatus
    win._on_failed("genome_annotate", "RuntimeError: boom", "trace")      # the fault path never reaches _on_genome_annotate
    win._annot_watch.timeout.emit()
    assert bar.isHidden() and not win._annot_inflight


_CLASH = "FAIL this genome has a part-finished annotation run with different settings (a) than the ones requested (b)."


def _answer(monkeypatch, label_part):
    seen = []

    def exec_(self):
        seen.append(self.text() + " " + self.informativeText())
        self._chosen = next((b for b in self.buttons() if label_part in b.text().lower()), None)
        return 0
    monkeypatch.setattr(QMessageBox, "exec", exec_)
    monkeypatch.setattr(QMessageBox, "clickedButton", lambda self: getattr(self, "_chosen", None))
    return seen


def test_discard_and_start_over_starts_the_run_the_user_asked_for(win, monkeypatch):
    calls = _capture(win, ops=("genome_annotate", "genome_annotate_reset", "genome_annotate_log"))
    win._annot_cost_state = {"sha256": _HUMAN["sha256"], "accession": _HUMAN["accession"]}
    win._start_genome_annotate(_HUMAN["accession"], "Homo sapiens", "quick", None, True)
    first = calls[0][1]
    seen = _answer(monkeypatch, "discard")
    win._on_genome_annotate({"ok": False, "error": _CLASH, "annotation_state": "part_finished"})
    assert "part-finished" in seen[0]
    assert calls[-1][0] == "genome_annotate_reset"
    win._on_genome_annotate_reset({"ok": True})
    again = calls[-1]
    assert again[0] == "genome_annotate"
    assert {k: again[1][k] for k in ("assemblyAccession", "species", "sensitivity", "include_uncurated", "sha256")} == \
        {k: first[k] for k in ("assemblyAccession", "species", "sensitivity", "include_uncurated", "sha256")}
    assert win._annot_inflight


def test_a_finished_annotation_is_not_described_as_part_finished(win, monkeypatch):
    _capture(win, ops=("genome_annotate", "genome_annotate_reset", "genome_annotate_log"))
    win._annot_cost_state = {"accession": _HUMAN["accession"]}
    win._start_genome_annotate(_HUMAN["accession"], "Homo sapiens", "default")
    seen = _answer(monkeypatch, "keep")
    win._on_genome_annotate({"ok": False, "annotation_state": "complete",
                             "error": "FAIL a completed annotation of this genome exists with different settings (a) than (b)."})
    assert "part-finished" not in seen[0] and "finished annotation" in seen[0]


def test_a_declined_discard_starts_nothing(win, monkeypatch):
    calls = _capture(win, ops=("genome_annotate", "genome_annotate_reset", "genome_annotate_log"))
    win._start_genome_annotate(_HUMAN["accession"], "Homo sapiens", "default")
    n = len(calls)
    _answer(monkeypatch, "keep")
    win._on_genome_annotate({"ok": False, "error": _CLASH, "annotation_state": "part_finished"})
    assert len(calls) == n


# ======================= cautions, failure wording, family-level note =======================
def test_a_failed_design_does_not_blame_the_parameters(win):
    _analyse(win, _example())
    win._design_inflight = True
    win._reset_buttons("primers")
    text = _labels(win.primBody.parentWidget())
    assert "Primer design" in text and "parameters" not in text.lower()
    assert win.designBtn.isEnabled()


def test_design_cautions_from_the_engine_are_shown_with_the_pairs(win):
    win._on_primers({"candidates": [_cand(1)],
                     "cautions": ["2 ambiguous template base(s) were treated as N for Primer3; no primer covers them."]})
    assert "ambiguous template base" in _labels(win.primBody.parentWidget())
    t = _tables(win.primBody.parentWidget())["TEagle_primers"]
    assert "ambiguous template base" in " ".join(t.export_notes)


def test_pcr_cautions_from_the_engine_reach_the_result_and_its_export(win):
    _analyse(win, _example())
    win._on_primers({"candidates": [_cand(1)]})
    win._add_pcr_pair(win.state["candidates"][0])
    msg = "3′ exact bases (25) exceed the shortest annealed core (20): primer bases 5′ of it are reported as a tail."
    _run_pcr_with(win, [{"amplicons": [_pcr_hit()], "status": "complete", "cautions": [msg]}], win.state["pcrPairs"])
    assert "exceed the shortest annealed core" in _labels(win.pcrBody.parentWidget())
    assert "exceed the shortest annealed core" in " ".join(_tables(win.pcrBody.parentWidget())["TEagle_amplicons"].export_notes)


def test_the_primer_table_says_a_pair_inside_an_element_can_amplify_every_copy(win):
    win._on_primers({"candidates": [_cand(1)]})
    text = _labels(win.primBody.parentWidget())
    assert "every copy" in text and "whole-genome scan" in text
    assert "locus-specific" not in text and "family-generic" not in text      # the release guardrail's own vocabulary


def test_primer_and_pcr_fields_carry_units_and_plain_language_tooltips(win):
    assert win.pfields["pMaxS"].toolTip() and "36" in win.pfields["pMaxS"].toolTip()
    labels = {l.text() for l in win.findChildren(QLabel)}
    for expected in ("MAX SIZE (NT)", "TM MIN (°C)", "PAIRS TO RETURN", "GC CLAMP (BASES)", "MAX POLY-X (BASES)"):
        assert expected in labels, expected
    assert "annealed 3′ core" in win.pfields["pcrMM"].toolTip()
    assert "polymerase" in win.pfields["pcrTP"].toolTip()


def test_the_scan_manifest_is_what_export_manifest_writes_after_a_pcr_run(win, monkeypatch):
    _analyse(win, _example())
    win._on_primers({"candidates": [_cand(1)]})
    win._add_pcr_pair(win.state["candidates"][0])
    _run_pcr_with(win, [{"amplicons": [_pcr_hit()], "status": "complete",
                         "provenance": {"runType": "in-silico-pcr", "input": {"id": "x"}}}], win.state["pcrPairs"])
    _prime_scan(win, {"id": "P1", "left_seq": F, "right_seq": R, "product_size": 300})
    win._on_genome_pcr(_scan_reply([]))
    sent = []
    monkeypatch.setattr(win, "_export_manifest_payload", lambda manifests, base: sent.append([m["runType"] for m in manifests]))
    win._export_manifest()
    assert sent == [["genome-scan"]]


def test_long_pair_names_are_shortened_on_the_gel_only_and_stay_unique(win):
    _analyse(win, _example())
    for name in ("my_long_qPCR_pair_name", "my_long_qPCR_pair_name_2"):
        win.pcrCustomF.setText(F); win.pcrCustomR.setText(R); win.pcrCustomName.setText(name)
        win._add_custom_pair()
    pairs = win.state["pcrPairs"]
    assert len(pairs) == 1                                        # the same two primers are one queued pair, however named
    win.state["pcrPairs"] = [{**pairs[0], "id": "my_long_qPCR_pair_name"}, {**pairs[0], "id": "my_long_qPCR_pair_name_2",
                                                                         "left_seq": "ACGT" * 5}]
    _run_pcr_with(win, [{"amplicons": [_pcr_hit()], "status": "complete"}] * 2, win.state["pcrPairs"])
    lanes = [l["label"] for l in win.state["lastPcr"]["lanes"]]
    assert all(len(x) <= 10 for x in lanes) and len(set(lanes)) == 2
    t = _tables(win.pcrBody.parentWidget())["TEagle_amplicons"]
    assert [r[0] for r in t.rows_data()] == ["my_long_qPCR_pair_name", "my_long_qPCR_pair_name_2"]


def _unready_goal(win):
    _analyse(win, _example())
    win.pGoal.setCurrentIndex(win.pGoal.findData("custom"))
    win.pGoalTo.setValue(100); win.pGoalFrom.setValue(300)
    assert not win.designBtn.isEnabled() and not win._goal_ready()[0]


def test_a_design_that_returns_leaves_design_primers_disabled_when_the_goal_is_not_ready(win):
    _unready_goal(win)
    win._design_inflight = True
    win._on_primers({"candidates": [_cand(1)]})
    assert not win.designBtn.isEnabled() and win.designBtn.text().upper() == "DESIGN PRIMERS"


def test_a_failed_design_leaves_design_primers_disabled_when_the_goal_is_not_ready(win):
    _unready_goal(win)
    win._design_inflight = True
    win._reset_buttons("primers")
    assert not win.designBtn.isEnabled()
