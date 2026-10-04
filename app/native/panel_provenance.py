"""Provenance panel and run-manifest exports."""
from __future__ import annotations

import html
import json
import os
from PySide6.QtCore import QTimer, Qt
from PySide6.QtWidgets import (
    QApplication, QFileDialog, QGridLayout, QHBoxLayout, QLabel, QPlainTextEdit, QPushButton, QVBoxLayout, QWidget)
import widgets
import theme as theme_mod
from ui_common import (
    _Combo, _clear_layout, _empty, _manifest_export_stem, _note, _sl)


class ProvenanceMixin:
    # =================== provenance ===================
    def _prov_list(self):
        """Runs since the current analysis was accepted (or, with no analysis, since the specimen last changed), oldest
        first; each entry holds the manifests of one run (a PCR batch holds one per pair). The list restarts whenever
        the analysis it belongs to is replaced or invalidated, so it never mixes runs that no longer describe the screen."""
        analysis = self.state.get("analysis") if self.state.get("analysis_state") == "current" else None
        anchor = getattr(self, "_prov_analysis", None)
        if not hasattr(self, "_prov_runs") or (anchor is not None and anchor is not analysis):
            self._prov_runs = []
            self._prov_sel = None
        self._prov_analysis = analysis
        return self._prov_runs

    def _pcr_group(self, m):
        """The manifests of the whole PCR batch that `m` heads: the card is drawn for the first pair as soon as the
        batch is in, and the batch has one manifest per pair."""
        run = getattr(self, "_pcr_run", None)
        if m.get("runType") != "in-silico-pcr" or not run:
            return None
        provs = [r["provenance"] for r in run.get("results") or [] if isinstance(r, dict) and r.get("provenance")]
        return provs if provs and provs[0] is m else None

    def _render_provenance(self, m, group=None):
        self.card_prov.expand()
        runs = self._prov_list()
        members = list(group) if group else (self._pcr_group(m) or [m])
        runs.append({"manifests": members})
        # the latest run is what 'Export this run' and the other panels' callers see; choosing another run in the list
        # only changes what is shown
        self.state["prov_manifest"] = m
        self.state["prov_manifest_all"] = list(members)
        self._prov_sel = (len(runs) - 1, 0)
        self._build_prov_card()
        self._uppercase_buttons()

    def _build_prov_card(self):
        import report
        card = self.card_prov
        card.clear_body()
        body = card.bodylay
        self._prov_report_block(body)
        runs = self._prov_list()
        op, member = self._prov_sel if self._prov_sel else (len(runs) - 1, 0)
        body.addWidget(_sl("Runs"))
        combo = _Combo(); combo.setAccessibleName("Run to show")
        for i, run in enumerate(runs):
            k = len(run["manifests"])
            combo.addItem(report.run_title(i + 1, run["manifests"][0], f" · {k} pairs" if k > 1 else ""), i)
        combo.setCurrentIndex(op)
        combo.setToolTip("Every run since the specimen or its analysis last changed. Earlier runs are cleared because "
                         "they no longer describe the current result.")
        pair = _Combo(); pair.setAccessibleName("Pair to show")
        pair.setToolTip("A batch of in-silico PCR runs has one sealed manifest per primer pair.")
        holder = QWidget(); QVBoxLayout(holder).setContentsMargins(0, 0, 0, 0); holder.layout().setSpacing(theme_mod.sp(9))
        body.addWidget(combo); body.addWidget(pair); body.addWidget(holder)
        self._prov_view = holder

        def fill_pairs(op_index, keep=0):
            pair.blockSignals(True); pair.clear()
            mans = runs[op_index]["manifests"]
            for j, name in enumerate(self._pair_names(mans)):
                pair.addItem(name, j)
            pair.setCurrentIndex(min(keep, pair.count() - 1))
            pair.setVisible(len(mans) > 1)
            pair.blockSignals(False)

        def show(op_index):
            fill_pairs(op_index)
            self._render_run_view(op_index, 0)
            self._uppercase_buttons(holder)

        combo.currentIndexChanged.connect(show)
        pair.currentIndexChanged.connect(lambda j: (self._render_run_view(combo.currentIndex(), j),
                                                    self._uppercase_buttons(holder)))
        fill_pairs(op, member)
        self._render_run_view(op, member)

    def _pair_names(self, manifests):
        """Name of each manifest's primer pair as the gel and the tables call it (P1, P2 … when it is not known)."""
        known = {(p.get("fwd"), p.get("rev")): p.get("label") for p in (self.state.get("lastPcr") or {}).get("pairs") or []}
        out = []
        for k, m in enumerate(manifests, 1):
            prm = m.get("parameters") or {}
            out.append(known.get((prm.get("fwd"), prm.get("rev"))) or f"P{k}")
        return out

    def _render_run_view(self, op_index, member_index=0):
        import report
        runs = self._prov_list()
        if not (0 <= op_index < len(runs)):
            return
        mans = runs[op_index]["manifests"]
        m = mans[max(0, min(member_index, len(mans) - 1))]
        self._prov_sel = (op_index, member_index)
        lay = self._prov_view.layout()
        _clear_layout(lay)
        lay.addWidget(self._kv_section("Input", report.input_rows(m)))
        lay.addWidget(self._kv_section("Run and seal", report.seal_rows(m)))
        sw = m.get("software") or []
        if sw:
            lay.addWidget(self._kv_section("Software", [(s.get("name", "—"), s.get("version") or "—", None) for s in sw]))
        dbs = m.get("databases") or []
        if dbs:
            self._db_sections(lay, m, dbs)
        lay.addWidget(_sl("Parameters"))
        sections = report.param_sections(m)
        if not sections:
            lay.addWidget(_empty("This run sealed no parameters (defaults only)."))
        raw = {}
        for group, rows in sections:
            for r in rows:
                if r["raw"]:
                    raw[r["key"]] = json.loads(r["raw"])
            table = self._param_grid(rows)
            if group == report.INTERNALS:
                self._disclosure(lay, f"{group} · {len(rows)} values", table, expanded=False)
            else:
                box = QWidget(); bl = QVBoxLayout(box); bl.setContentsMargins(0, 0, 0, 0); bl.setSpacing(theme_mod.sp(3))
                head = QLabel(f"<b>{html.escape(group)}</b>"); head.setObjectName("cardnote"); head.setTextFormat(Qt.RichText)
                bl.addWidget(head); bl.addWidget(table)
                lay.addWidget(box)
        if raw:
            box = QPlainTextEdit(json.dumps(raw, indent=2, ensure_ascii=False, sort_keys=True)); box.setReadOnly(True)
            box.setMaximumHeight(theme_mod.sp(220)); box.setAccessibleName("Raw sealed parameters (JSON)")
            self._disclosure(lay, "Raw sealed parameters (JSON)", box, expanded=False)
        env = m.get("environment") or {}
        if env:
            lay.addWidget(self._kv_section(f"Environment ({report.RESULT_NOTE})",
                                           [(k, str(env[k]), None) for k in ("os", "machine", "python") if k in env]))
        notrun = m.get("notRun") or []
        if notrun:
            n = QLabel("<br>".join("• " + html.escape(str(x)) for x in notrun)); n.setObjectName("cardnote")
            n.setTextFormat(Qt.RichText); n.setWordWrap(True); n.setTextInteractionFlags(Qt.TextSelectableByMouse)
            lay.addWidget(_sl("Not run in this step")); lay.addWidget(n)
        refs = m.get("references") or []
        if refs:
            r = QLabel("<br>".join(f"<b>{html.escape(str(x.get('name', '')))}</b> — {html.escape(str(x.get('citation', '')))}"
                                   + (f" doi:{html.escape(str(x['doi']))}" if x.get("doi") else "") for x in refs))
            r.setObjectName("cardnote"); r.setTextFormat(Qt.RichText); r.setWordWrap(True)
            r.setTextInteractionFlags(Qt.TextSelectableByMouse)
            lay.addWidget(_sl("References (source-verified)")); lay.addWidget(r)
        self._prov_export_row(lay, op_index, len(mans))

    def _db_sections(self, lay, m, dbs):
        import report
        rows = []
        for d in dbs:
            name = str(d.get("name") or "—")
            if d.get("version"):
                rows.append((name, str(d["version"]), None))
            elif d.get("file") and not d.get("sha256"):
                rows.append((name, str(d["file"]), None))
            if d.get("sha256"):
                rows.append((f"{name} SHA-256", d["sha256"], "hash"))
        profiles, note = report.panel_profiles(m)
        if profiles:
            rows.append(("Pfam models", f"{len(profiles)} (accessions below)", None))
        lay.addWidget(self._kv_section("Databases", rows))
        if note:
            lay.addWidget(_note(note, "info"))
        if profiles:
            grid = QGridLayout(); grid.setHorizontalSpacing(theme_mod.sp(14)); grid.setVerticalSpacing(theme_mod.sp(2))
            for c, head in enumerate(("Code", "Model", "Accession")):
                h = QLabel(head); h.setObjectName("kdim"); grid.addWidget(h, 0, c)
            for r, p in enumerate(profiles, 1):
                for c, text in enumerate((p["code"], p["name"], p["accession"])):
                    v = QLabel(text); v.setObjectName("cardmeta"); v.setTextInteractionFlags(Qt.TextSelectableByMouse)
                    v.setToolTip(p["label"]); grid.addWidget(v, r, c)
            grid.setColumnStretch(3, 1)
            box = QWidget(); box.setLayout(grid)
            self._disclosure(lay, f"Profile accessions ({len(profiles)})", box, expanded=False)

    def _param_grid(self, rows):
        box = QWidget(); grid = QGridLayout(box); grid.setContentsMargins(0, 0, 0, 0)
        grid.setHorizontalSpacing(theme_mod.sp(10)); grid.setVerticalSpacing(theme_mod.sp(3)); grid.setColumnStretch(1, 1)
        for r, row in enumerate(rows):
            k = QLabel(row["label"]); k.setObjectName("kdim"); k.setWordWrap(True); k.setToolTip(row["tip"])
            v = QLabel(row["text"]); v.setObjectName("cardmeta"); v.setWordWrap(True); v.setToolTip(row["tip"])
            v.setTextInteractionFlags(Qt.TextSelectableByMouse)
            grid.addWidget(k, r, 0, Qt.AlignTop); grid.addWidget(v, r, 1, Qt.AlignTop)
        grid.setColumnMinimumWidth(0, theme_mod.sp(230))      # one label column width for every group
        return box

    def _kv_section(self, title, rows):
        """A labelled block of key/value rows. A row is (label, value, kind); a 'hash' value is shown in full, in blocks
        of eight characters so it wraps, beside a Copy button that puts the exact 64 characters on the clipboard."""
        import report
        box = QWidget(); lay = QVBoxLayout(box); lay.setContentsMargins(0, 0, 0, 0); lay.setSpacing(theme_mod.sp(4))
        lay.addWidget(_sl(title))
        grid = QGridLayout(); grid.setHorizontalSpacing(theme_mod.sp(10)); grid.setVerticalSpacing(theme_mod.sp(3))
        grid.setColumnStretch(1, 1)
        for r, (label, value, kind) in enumerate(rows):
            k = QLabel(label); k.setObjectName("kdim"); k.setWordWrap(True)
            v = QLabel(report.group_hash(value) if kind == "hash" else str(value)); v.setObjectName("cardmeta")
            v.setWordWrap(True); v.setTextInteractionFlags(Qt.TextSelectableByMouse)
            grid.addWidget(k, r, 0, Qt.AlignTop); grid.addWidget(v, r, 1, Qt.AlignTop)
            if kind == "hash":
                b = QPushButton("Copy"); b.setProperty("sm", True); b.setAccessibleName(f"Copy {label}")
                b.setToolTip(f"Copy the full 64-character value of {label}")
                b.clicked.connect(lambda _=False, text=value, btn=b: self._copy_value(text, btn))
                grid.addWidget(b, r, 2, Qt.AlignTop)
        grid.setColumnMinimumWidth(0, theme_mod.sp(150))
        lay.addLayout(grid)
        return box

    def _copy_value(self, text, button=None):
        QApplication.clipboard().setText(str(text))
        if button is not None:
            old = button.text()
            button.setText(widgets.caps("Copied"))

            def restore(b=button, t=old):
                try:
                    b.setText(t)
                except RuntimeError:
                    pass
            QTimer.singleShot(1500, restore)

    def _prov_export_row(self, lay, op_index, n_members):
        """Exports for the run on show. The analysis' own manifest and complete result stay reachable whichever run is
        selected, while the analysis is current."""
        runs = self._prov_list()
        shown = runs[op_index]["manifests"][0]
        analysis = self.state.get("analysis") if self.state.get("analysis_state") == "current" else None
        row = QHBoxLayout(); row.addStretch(1)
        b = QPushButton("Export this run (.json)"); b.setProperty("sm", True); b.setProperty("resultExport", True)
        b.setToolTip(f"Write all {n_members} pair manifests of this batch as one JSON file." if n_members > 1 else
                     "Write the full run-provenance manifest (input hash, tool and database versions, parameters, "
                     "checksums) as a JSON file so the result stays reproducible outside TEagle.")
        b.clicked.connect(self._export_selected_run); row.addWidget(b)
        if len(runs) > 1:
            a = QPushButton("Export all runs (.json)"); a.setProperty("sm", True); a.setProperty("resultExport", True)
            a.setToolTip("Write every sealed manifest listed above, in order, as one JSON file.")
            a.clicked.connect(self._export_all_runs); row.addWidget(a)
        lay.addLayout(row)
        if isinstance(analysis, dict):
            row2 = QHBoxLayout(); row2.addStretch(1)
            if shown.get("runType") != "analysis":
                ab = QPushButton("Export analysis manifest (.json)"); ab.setProperty("sm", True)
                ab.setProperty("resultExport", True)
                ab.setToolTip("Write the sealed provenance manifest of the analysis itself.")
                ab.clicked.connect(self._export_analysis_manifest); row2.addWidget(ab)
            rb = QPushButton("Export analysis result (.json)"); rb.setProperty("sm", True)
            rb.setProperty("resultExport", True)
            rb.setToolTip("Write the complete versioned analysis result: every record, stage status, "
                          "scientific output, reference and verified provenance binding.")
            rb.clicked.connect(self._export_analysis_result); row2.addWidget(rb)
            lay.addLayout(row2)

    def _export_manifest_payload(self, manifests, base):
        """Strictly verify, serialize, and atomically write one or more sealed manifests."""
        import json
        from teagle_core import provenance
        allm = [m for m in manifests if isinstance(m, dict)]
        if not allm:
            return self._banner("No run provenance to export — run the analysis first.", "info")
        path, _ = QFileDialog.getSaveFileName(self, "Export provenance manifest",
                                              f"{base}.json", "JSON (*.json)")
        if not path:
            return
        if not path.lower().endswith(".json"):
            path += ".json"
        try:
            encoded_runs = [provenance.manifest_json_bytes(manifest) for manifest in allm]
            payload = (json.dumps({"runs": [json.loads(encoded) for encoded in encoded_runs]},
                                  indent=2, ensure_ascii=False, allow_nan=False, sort_keys=True).encode("utf-8")
                       if len(encoded_runs) > 1 else encoded_runs[0])
        except (TypeError, ValueError):
            return self._banner("The provenance manifest failed its seal or strict JSON check and was not exported.")
        if not self._write_bytes(path, payload):
            return
        self._banner(f"Provenance manifest written to {os.path.basename(path)}.", level="success")

    def _write_bytes(self, path, payload):
        """Write an export atomically; a disk problem is reported as one, never as a seal failure."""
        try:
            with widgets.atomic_write(path, "wb") as fh:
                fh.write(payload)
        except OSError:
            self._banner(f"“{os.path.basename(path)}” could not be written. The file may be open in another "
                         "program, or the folder may be read-only or no longer available. The provenance itself "
                         "is intact; choose another location and export again.")
            return False
        return True

    def _export_analysis_manifest(self):
        res = self.state.get("analysis")
        manifest = (res or {}).get("provenance") if isinstance(res, dict) else None
        if self.state.get("analysis_state") != "current" or not isinstance(manifest, dict):
            return self._banner("No current analysis to export — run analysis again first.", "info")
        stem = ("multiple_records" if len(res.get("records") or []) > 1
                else _manifest_export_stem(manifest.get("input", {}), "analysis"))
        return self._export_manifest_payload([manifest], f"TEagle_{stem}_manifest")

    def _export_manifest(self):
        allm = self.state.get("prov_manifest_all") or ([self.state["prov_manifest"]]
                                                        if self.state.get("prov_manifest") else [])
        return self._export_run_set(allm)

    def _export_run_set(self, allm):
        """Write the manifests of one run (or batch); each is re-verified against its seal first."""
        if not allm:
            return self._banner("No run provenance to export — run an analysis first.", "info")
        if allm[0].get("runType") == "analysis" and self.state.get("analysis_state") != "current":
            return self._banner("That analysis provenance is stale — run analysis again before exporting.", "info")
        if len(allm) > 1:
            stem = "batch"
        elif allm[0].get("runType") == "analysis" and len(self.state.get("records") or []) > 1:
            stem = "multiple_records"
        else:
            stem = _manifest_export_stem(allm[0].get("input", {}), "run")
        return self._export_manifest_payload(allm, f"TEagle_{stem}_manifest")

    def _export_selected_run(self):
        """'Export this run': the run chosen in the provenance list, which is the latest one until another is picked."""
        runs = self._prov_list()
        sel = self._prov_sel[0] if self._prov_sel else len(runs) - 1
        return self._export_run_set(runs[sel]["manifests"] if 0 <= sel < len(runs) else [])

    def _export_all_runs(self):
        """The session's provenance bundle: every sealed manifest since the analysis, in the order the runs happened."""
        allm = [m for run in self._prov_list() for m in run["manifests"]]
        if not allm:
            return self._banner("No run provenance to export — run an analysis first.", "info")
        recs = self.state.get("records") or []
        stem = "multiple_records" if len(recs) > 1 else _manifest_export_stem(allm[0].get("input", {}), "run")
        return self._export_manifest_payload(allm, f"TEagle_{stem}_all_runs_manifest")

    def _export_analysis_result(self):
        import engine
        result = self.state.get("analysis")
        manifest = result.get("provenance") if isinstance(result, dict) else None      # the analysis' own seal, whatever step ran since
        if (self.state.get("analysis_state") != "current" or not isinstance(manifest, dict) or
                manifest.get("runType") != "analysis"):
            return self._banner("No current analysis result to export — run analysis again first.", "info")
        rid = ("multiple_records" if len(result.get("records") or []) > 1 else
               _manifest_export_stem(manifest.get("input", {}), "analysis"))
        path, _ = QFileDialog.getSaveFileName(self, "Export complete analysis result",
                                              f"TEagle_{rid}_analysis.json", "JSON (*.json)")
        if not path:
            return
        if not path.lower().endswith(".json"):
            path += ".json"
        try:
            payload = engine.analysis_result_json_bytes(result)
        except ValueError:
            return self._banner("The analysis result failed its provenance binding or strict JSON check and was not exported.")
        if not self._write_bytes(path, payload):
            return
        self._banner(f"Complete analysis result written to {os.path.basename(path)}.", level="success")
