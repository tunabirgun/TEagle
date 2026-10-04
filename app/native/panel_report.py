"""Assay record report and classification export."""
from __future__ import annotations

import os
from PySide6.QtWidgets import QApplication, QFileDialog, QHBoxLayout, QLabel, QPushButton
import figures
from figures import svg_genome
import widgets
from teagle_core import sequtil, oligoqc, __version__ as APP_VERSION
from ui_common import (
    PANEL, _amp_call, _call_view, _caution_items, _coding_view, _explanation_view, _note_view, _rank_line,
    _rank_view, _record_export_stem, _record_public_id, _sl, _te_headline)


class ReportMixin:
    # =================== assay record ===================
    def _prov_report_block(self, layout):
        """Record-keeping actions for the record on show: a notebook entry, a one-page assay report and the
        classification table. Present only while the analysis is current."""
        if self.state.get("analysis_state") != "current" or not isinstance(self.state.get("analysis"), dict):
            return
        layout.addWidget(_sl("Record this assay"))
        row = QHBoxLayout()
        copy = QPushButton("Copy notebook entry"); copy.setProperty("sm", True); copy.setProperty("resultExport", True)
        copy.setToolTip("Copies a plain-text and formatted summary (call and cautions, primer pairs with locations and "
                        "QC, in-silico PCR result, design and QC conditions, run seals) for a lab notebook. Figures are "
                        "in the assay report.")
        copy.clicked.connect(self._copy_notebook_entry)
        rep = QPushButton("Assay report ▾"); rep.setProperty("sm", True); rep.setProperty("resultExport", True)
        rep.setToolTip("One-page record: call and cautions, element map, primer pairs, predicted gel, conditions and "
                       "provenance.")
        rep.clicked.connect(lambda _=False, b=rep: self._assay_report_menu(b))
        tsv = QPushButton("Export classification (TSV)"); tsv.setProperty("sm", True); tsv.setProperty("resultExport", True)
        tsv.setToolTip("One row per record in TEsorter's .cls.tsv column order (#TE, Order, Superfamily, Clade, Complete, "
                       "Strand, Domains), followed by the RepeatMasker-style class, the Wicker code and TEagle's own call, "
                       "outcome, confidence and completeness tier, for joining against TEsorter, EDTA or RepeatMasker "
                       "results. Clade is NA because TEagle does not resolve clades. When record IDs differ between tools, "
                       "join on Sequence_SHA256.")
        tsv.clicked.connect(self._export_classification)
        for b in (copy, rep, tsv):
            row.addWidget(b)
        row.addStretch(1)
        layout.addLayout(row)
        cap = QLabel(f"Covers the record selected in panel {PANEL['struct']} and the primer pairs designed on it. "
                     "Coordinates are 0-based, half-open [start, end).")
        cap.setObjectName("orient"); cap.setWordWrap(True)
        layout.addWidget(cap)

    def _export_classification(self):
        import report
        recs = self.state.get("records") or []
        if self.state.get("analysis_state") != "current" or not recs:
            return self._banner("No current analysis — run analysis again before exporting the classification table.", "info")
        stem = "multiple_records" if len(recs) > 1 else _record_export_stem(recs[0], "record")
        self._save_text(report.classification_tsv(recs, _rank_view, _call_view), f"TEagle_{stem}_classification.cls", "tsv")

    def _call_display(self, rec):
        """The call as the results card words it, from the same module-level views, for the notebook and the report."""
        cl = rec.get("classification") or {}
        te = cl.get("te_class", "—")
        held = cl.get("outcome") == "withheld"
        unassessed = bool(cl.get("domains_unavailable")) and te in ("none", "unclassified")
        headline = ("Not assessed — domain scan unavailable" if unassessed else
                    "Terminal direct-repeat candidate" if held and te == "repeat/structural-only" else _te_headline(te))
        chip = (("NOT ASSESSED" if unassessed else "CLASS & ORDER WITHHELD" if te in ("unclassified", "none")
                 else "ORDER WITHHELD") if held or unassessed else None)
        family = ("no call made — the protein-domain scan did not run" if unassessed and te == "none"
                  else str(cl.get("superfamily") or ""))
        if family.lower().startswith(headline.lower()) or family == "—":
            family = ""
        comp = cl.get("completeness") or {}
        status, _ = self._scan_state(rec, cl)
        _, arch = self._orders(cl, status)
        tier = comp.get("tier") or ""
        if tier and (comp.get("coding") or {}).get("truncated_frame") is True and "intact" in tier.lower():
            tier += (" at the domain level — but a domain-bearing ORF runs off the end of the supplied sequence, so the "
                     "element may be truncated")
        rank = _rank_view(cl)
        return {"headline": headline, "chip": chip, "superfamily": family,
                "confidence": _call_view(rec)["confidence"], "rank_line": _rank_line(rank) if rank else "",
                "why": _explanation_view(rec) if cl.get("explanation") else "",
                "completeness": (f"{tier} · domain order {arch} · not detected: "
                                 f"{', '.join(comp.get('missing') or []) or 'none'}") if tier else "",
                "coding": _coding_view(rec) if comp.get("coding", {}).get("note") else "",
                "cautions": _caution_items(rec) + [_note_view(n, rec) for n in rec.get("notes") or []]}

    def _report_refusal(self):
        recs = self.state.get("records") or []
        i = self.state.get("record_index")
        if self.state.get("analysis_state") == "current" and i is not None and 0 <= i < len(recs):
            return ("The selected record was not analysed, so there is nothing to record. Select another record in "
                    f"panel {PANEL['struct']}.")
        return "No current analysis — run analysis again before recording this assay."

    def _report_context(self):
        """Everything the notebook entry and the report say, as plain data, or None when there is no current analysis of
        a selected, classified record. Only named fields are read: a candidate carries the text it was designed on, and
        that text (which may hold the user's own header) is never copied out."""
        import engine, report
        from datetime import datetime, timezone
        res = self.state.get("analysis")
        recs = self.state.get("records") or []
        i = self.state.get("record_index")
        if self.state.get("analysis_state") != "current" or not isinstance(res, dict) or i is None or not 0 <= i < len(recs):
            return None
        try:
            engine.analysis_result_document(res)
        except ValueError:
            return None
        rec = recs[i]
        if not rec.get("valid", True) or not rec.get("classification") or _call_view(rec)["status"] == "failed":
            return None
        manifest = res.get("provenance") or {}
        sha = sequtil.sha256(rec.get("seq") or "")
        line = report.source_line(manifest.get("input"))
        source = ("Source: " + line) if line else (
            "Source: pasted or uploaded sequence (header replaced by a content-derived ID)"
            if rec.get("identitySource") == "local_input" else "Source: pasted or uploaded sequence")
        comp = rec.get("composition") or {}
        ctx = {"version": APP_VERSION, "generated": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M"),
               "specimen": {"id": _record_public_id(rec), "position": (i + 1, len(recs)),
                            "length": comp.get("length", len(rec.get("seq") or "")),
                            "gc": comp.get("gc"), "sha256": sha, "source": source},
               "call": self._call_display(rec)}
        ctx.update(self._report_assay(rec, sha, manifest))
        ctx["_rec"] = rec
        return ctx

    def _report_assay(self, rec, sha, analysis_manifest):
        """Pairs, in-silico PCR, genome scan, conditions and seals for the record `sha` identifies. A pair or PCR lane
        is included only when its sealed input hash equals the record's, so a result computed on another record,
        flank or sub-region never reaches this record's report."""
        import report
        feats = report.feature_spans(rec)

        def bound(c):
            return ((c.get("design_provenance") or {}).get("input") or {}).get("sha256") == sha

        cands = self.state.get("candidates") or []
        designed = [c for c in cands if bound(c)]
        last = self.state.get("lastPcr") or {}
        tested = {(p["fwd"], p["rev"]): k for k, p in enumerate(last.get("pairs") or [])}
        staged = [c for c in self.state.get("pcrPairs") or [] if c.get("custom") or bound(c)]
        chosen, seen, note = [], set(), None
        for c in staged + [c for c in designed if (c["left_seq"], c["right_seq"]) in tested]:
            if (c["left_seq"], c["right_seq"]) not in seen:
                seen.add((c["left_seq"], c["right_seq"])); chosen.append(c)
        if not chosen and designed:
            chosen, note = designed, "No pair has been sent to in-silico PCR yet; every designed pair is listed."
        pcr_man = {}
        for run in self._prov_list():
            for m in run["manifests"]:
                if m.get("runType") == "in-silico-pcr":
                    p = m.get("parameters") or {}
                    pcr_man[(p.get("fwd"), p.get("rev"))] = m

        def man_bound(m):
            return ((m.get("parameters") or {}).get("assay_input") or m.get("input") or {}).get("sha256") == sha

        statuses = last.get("statuses") or []
        words = {"complete": "complete", "complete_zero_products": "complete, no product found",
                 "partial_capped": "stopped early (search cap)", "failed": "failed — not assessed",
                 "not_run": "not run"}
        pairs, lines, lanes, runs, stopped = [], [], [], [], False
        for c in chosen:
            key = (c["left_seq"], c["right_seq"])
            m, k = pcr_man.get(key), tested.get(key)
            if k is not None and m is not None and not man_bound(m):
                k = None                                      # that PCR searched another sequence
            if k is not None and k >= len(last.get("lanes") or []):
                k = None
            label = last["pairs"][k]["label"] if k is not None else self._pair_label(c)
            d = {"label": label, "fwd": c["left_seq"], "rev": c["right_seq"]}
            if bound(c):
                lp, rp = c.get("left_pos"), c.get("right_pos")
                span = (lp[0], rp[1]) if lp and rp else None
                d.update(fwd_span=tuple(lp) if lp else None, rev_span=tuple(rp) if rp else None,
                         product=c.get("product_size"), amplicon_span=span, spans=report.overlaps(span, feats),
                         tm=(c.get("left_tm"), c.get("right_tm")), gc=(c.get("left_gc"), c.get("right_gc")),
                         qc=(c.get("qc") or {}).get("worst") or "n/a")
            else:
                d["note"] = "your own pair (no design record)"
            if k is None:
                d["note"] = (d["note"] + " · " if d.get("note") else "") + "not yet tested in silico"
            pairs.append(d)
            if k is None:
                continue
            lane = last["lanes"][k]
            status = (statuses[k]["status"] if k < len(statuses) else lane.get("status")) or "complete"
            prods = [f"{a.get('product_length', a.get('length')):,} bp ({_amp_call(a, lane.get('has_locus', True))})"
                     for a in lane.get("amplicons") or []]
            stopped = stopped or status == "partial_capped"
            lines.append(f"{label}: {words.get(status, status)}"
                         + (f" · {len(prods)} product(s): " + "; ".join(prods) if prods else ""))
            lanes.append(lane)
            if m is not None:
                runs.append((f"In-silico PCR {label}", m.get("manifestSha256")))
        if stopped:
            lines.append("In-silico PCR stopped early (search cap): the product list is a lower bound.")
        excluded = len(cands) - len(designed)
        scan = None
        gs = self.state.get("lastGenomeScan") or {}
        gp = gs.get("pair") or {}
        hit = next((d for d in pairs if (d["fwd"], d["rev"]) == (gp.get("fwd"), gp.get("rev"))), None)
        if gs and hit:
            crit, asm = gs.get("criteria") or {}, gs.get("assembly") or {}
            scan = [f"{hit['label']} on {gs.get('organism') or 'the scanned genome'} ({asm.get('accession') or '—'}): "
                    f"{(gs.get('summary') or {}).get('verdict') or 'scan complete'}",
                    f"isPcr search: at least {report.plain(crit.get('min_perfect'))} bp perfect 3′ match and "
                    f"{report.plain(crit.get('min_good'))} bp good match, products up to {report.plain(crit.get('max_size'))} bp. "
                    "Sites outside that rule are not reported; this is a prediction from the assembly, not a measured result."]
        conditions = []
        base = next((c for c in chosen if bound(c)), None)
        if base:
            p = (base.get("design_provenance") or {}).get("parameters") or {}
            tgt = p.get("primer3_target")
            conditions.append(
                f"Design: Tm {report.plain(p.get('min_tm'))}–{report.plain(p.get('max_tm'))} °C (opt {report.plain(p.get('opt_tm'))}) · "
                f"length {report.plain(p.get('min_size'))}–{report.plain(p.get('max_size'))} nt (opt {report.plain(p.get('opt_size'))}) · "
                f"GC {report.plain(p.get('min_gc'))}–{report.plain(p.get('max_gc'))} % · "
                f"product {report.plain(p.get('prod_min'))}–{report.plain(p.get('prod_max'))} bp · target "
                + (f"start {tgt[0]}, length {tgt[1]} nt (0-based)" if isinstance(tgt, (list, tuple)) and len(tgt) == 2
                   else "none (untargeted design)"))
            tm = report.tm_text(p.get("tm_conditions") or {})
            conditions.append("Primer3 Tm salt and oligo conditions: " +
                              (tm if tm else "Primer3 defaults (not recorded in this run's seal)"))
            qc = (base.get("qc") or {}).get("conditions") or {}
            if qc:
                warn = sorted({abs(w) for w in (oligoqc._HAIRPIN_WARN, oligoqc._DIMER_WARN, oligoqc._END_WARN)})
                conditions.append(
                    f"Structure QC: ΔG at {report.plain(qc.get('temp_c'))} °C, {report.plain(qc.get('mv_conc'))} mM Na⁺, "
                    f"{report.plain(qc.get('dv_conc'))} mM Mg²⁺, {report.plain(qc.get('dna_conc'))} nM oligo; warn at ΔG ≤ "
                    f"−{' / −'.join(f'{w:g}' for w in warn)} kcal/mol (advisory)")
        mm = next((pcr_man[(d["fwd"], d["rev"])] for d in pairs if (d["fwd"], d["rev"]) in pcr_man and
                   man_bound(pcr_man[(d["fwd"], d["rev"])])), None)
        if mm:
            p = mm.get("parameters") or {}
            conditions.append(f"In-silico PCR search: up to {report.plain(p.get('max_mm'))} mismatches per primer, "
                              f"{report.plain(p.get('tp'))} exact 3′ bases, products {report.plain(p.get('prod_min'))}–"
                              f"{report.plain(p.get('prod_max'))} bp")
        prim = (self.state.get("lastPrimers") or {}).get("provenance") or {}
        seals = [("Analysis", analysis_manifest.get("manifestSha256"))]
        if pairs and (prim.get("input") or {}).get("sha256") == sha:
            seals.append(("Primer design", prim.get("manifestSha256")))
        seals += runs
        if scan:
            gm = next((m for run in reversed(self._prov_list()) for m in run["manifests"]
                       if m.get("runType") == "genome-scan"), None)
            if gm:
                seals.append(("Genome scan", gm.get("manifestSha256")))
        panel = next((d for d in analysis_manifest.get("databases") or [] if "pfam" in str(d.get("name", "")).lower()), None)
        profiles, _why = report.panel_profiles(analysis_manifest)
        return {"pairs": pairs, "excluded": excluded, "pairs_note": note, "pcr": {"lines": lines}, "scan": scan,
                "conditions": conditions, "_lanes": lanes,
                "provenance": {"runs": [(a, b) for a, b in seals if b],
                               "software": "; ".join(f"{s.get('name')} {s.get('version')}"
                                                     for s in analysis_manifest.get("software") or []),
                               "panel": ("Pfam TE-domain panel: " + (f"{len(profiles)} models, " if profiles else "")
                                         + f"SHA-256 {panel['sha256']}") if panel and panel.get("sha256") else None}}

    def _copy_notebook_entry(self):
        import report
        ctx = self._report_context()
        if ctx is None:
            return self._banner(self._report_refusal(), "info")
        from PySide6.QtCore import QMimeData
        data = QMimeData()
        data.setText(report.notebook_text(ctx)); data.setHtml(report.report_html(ctx))
        QApplication.clipboard().setMimeData(data)
        self._banner("Notebook entry copied. Paste it into your notebook or a document.", "success")

    def _assay_report_menu(self, anchor):
        from PySide6.QtWidgets import QMenu
        menu = QMenu(self)
        menu.addAction("PDF — print or attach (figures as images)", lambda: self._export_assay_report("pdf"))
        menu.addAction("HTML — self-contained page (vector figures)", lambda: self._export_assay_report("html"))
        menu.exec(anchor.mapToGlobal(anchor.rect().bottomLeft()))

    def _export_assay_report(self, fmt):
        ctx = self._report_context()
        if ctx is None:
            return self._banner(self._report_refusal(), "info")
        path, _ = QFileDialog.getSaveFileName(
            self, "Save assay report", f"TEagle_{_record_export_stem(ctx['_rec'], 'record')}_assay.{fmt}",
            f"{fmt.upper()} (*.{fmt})")
        if not path:
            return
        if not path.lower().endswith("." + fmt):
            path += "." + fmt
        self._write_assay_report(path, fmt)

    def _write_assay_report(self, path, fmt):
        """Write the assay record as 'pdf' (figures as images) or 'html' (vector figures). Returns True when written; a
        missing or stale analysis writes nothing."""
        import report
        ctx = self._report_context()
        if ctx is None:
            self._banner(self._report_refusal(), "info")
            return False
        rec, lanes = ctx["_rec"], ctx["_lanes"]
        svgs = {}
        model = figures.gv_tracks_from_rec(rec)
        if model["tracks"]:
            svgs["map"] = svg_genome(model, {"start": 0, "end": model["length"]}, 920, "white", for_export=True, legend=True)
        if lanes:
            svgs["gel"] = figures.svg_gel({"lanes": lanes}, "white")

        def write(p):
            if fmt == "pdf":
                images, tags = {}, {}
                for name, svg in svgs.items():
                    images[f"teagle://{name}"] = report.svg_image(svg)
                    width = min(680, int(widgets._svg_size(svg)[0]))
                    tags[name] = f'<img src="teagle://{name}" width="{width}">'
                report.write_document_pdf(report.report_html(ctx, tags, qt=True), images, p)
            else:
                with widgets.atomic_write(p, "w", encoding="utf-8", newline="\n") as fh:
                    fh.write(report.report_html(ctx, svgs))
        return widgets.safe_export(self, path, write, ok=f"Assay report written to {os.path.basename(path)}.") is not None
