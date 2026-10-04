"""Whole-genome TE annotation."""
from __future__ import annotations

import html
import os
from PySide6.QtCore import QTimer, Qt
from PySide6.QtWidgets import QDialog, QHBoxLayout, QLabel, QMessageBox, QPushButton, QVBoxLayout
import figures
from widgets import BusyBar, DataTable, FigurePanel
import theme as theme_mod
from ui_common import (
    _export_table_btn, _safe_export_stem)


class GenomeAnnotateMixin:
    # =================== whole-genome TE annotation ===================
    # Wording taken verbatim from RepeatMasker's own -h text, not paraphrased: an earlier draft turned
    # "5-10% less sensitive" into "5-10x faster", which would have misled the user about both axes.
    _ANNOT_LOW_LIBRARY = 25          # see _on_annotate_budget for how this was chosen (measured, not assumed)
    _ANNOT_SENS = [("default", "Default — the balance RepeatMasker ships with"),
                   ("quick", "Quick (-q) — 2–5× faster, 5–10% less sensitive"),
                   ("slow", "Slow (-s) — 2–3× slower, 0–5% more sensitive")]

    def _on_annotate_budget(self, b):
        """Fill the cost panel from the measured WSL budget + the library's coverage for this lineage."""
        w = getattr(self, "_annot_cost_widgets", None)
        if not w:
            return
        cost, warn, start, gsel, ssel = w
        try:
            st = self._annot_cost_state
            # Budget jobs run concurrently. A slow Dfam lineage query can finish after the user has
            # switched to a custom FASTA (or another Dfam search set); never render that stale response
            # under the new selection.
            if (b.get("budget_request_id") != st.get("budget_request_id")
                    or b.get("library_mode") != st.get("library_mode")):
                return
            retry = getattr(self, "_annot_retry", None)
            if b.get("ok") is False:
                # the probe itself failed (WSL slow or down): none of the promised cost or disk checks ran, so the
                # run must not be armed, and the dialog says why instead of printing a table of 'None'
                cost.setText("TEagle could not read this machine's compute and disk budget, so the cost of this "
                             "annotation is unknown.")
                warn.setText("Budget check failed — " + str(b.get("error") or "the WSL backend did not answer") +
                             ". Check that WSL is running, then press Check again.")
                warn.show(); start.setEnabled(False)
                if retry is not None:
                    retry.show()
                return
            if retry is not None:
                retry.hide()
            # keep the MEASURED thread count so the run uses the number the dialog promised
            self._annot_threads = int(b.get("recommended_threads") or 4)
            fams = b.get("library_families_for_species")
            gb = b.get("disk_needed_gb")
            custom_path = st.get("custom_library")
            library_mode = st.get("library_mode")
            lines = [
                f"<b>This machine (WSL compute):</b> {b.get('cores')} cores · {b.get('mem_gb')} GB RAM. "
                f"TEagle will use {b.get('recommended_threads')} parallel "
                f"RepeatMasker jobs (limited by {b.get('limited_by')}).",
                f"<b>Storage:</b> {b.get('disk_note') or (str(b.get('avail_gb')) + ' GiB usable')}",
                f"<b>Disk needed:</b> about {gb} GB — the genome sequence is fetched once and kept, and the "
                "work chunks are removed as they finish.",
                "<b>Time:</b> a fly-sized genome (~140 Mb) takes minutes to about an hour here; a mammalian "
                "genome (~3 Gb) takes many hours to a day. The run continues in the background and can be "
                "re-started later — finished chunks are not repeated.",
            ]
            if custom_path:
                lines.append(f"<b>Selected library:</b> user-supplied FASTA {html.escape(os.path.basename(custom_path))}. "
                             "This file supplies the interspersed-repeat homology library, so only references "
                             "represented in it can receive homology-based TE calls. RepeatMasker’s separate "
                             "simple- and low-complexity screen can still produce non-TE repeat rows. To enter "
                             "the TE percentage, a custom reference also needs a recognized class/family in "
                             "its header, for example >name#LTR/Gypsy; unclassified matches remain visible but "
                             "are not silently counted as TEs.")
            else:
                selected_set = ("curated + uncurated" if library_mode == "dfam_uncurated"
                                else "curated-only")
                lines.append(f"<b>Library coverage:</b> the selected {selected_set} Dfam search set contains "
                             f"<b>{fams if fams is not None else 'an unknown number of'}</b> family models for "
                             f"{st.get('organism') or 'this lineage'}.")
            cost.setText("<br>".join(lines))
            if library_mode == "dfam_uncurated" and b.get("library_ready") is False:
                warn.setText(b.get("library_error") or
                             "The selected curated + uncurated Dfam search cannot be run as configured.")
                warn.show(); start.setEnabled(False); return
            if b.get("disk_ok") is False:
                warn.setText(f"Not enough usable disk — about {gb} GB is needed after applying the smaller "
                             "WSL filesystem and Windows backing-volume limit. Free space, or remove a cached "
                             "genome from this manager, then try again.")
                warn.show(); start.setEnabled(False); return
            # Threshold rationale: measured on this backend, a lineage Dfam covers usefully has hundreds
            # of models (D. melanogaster 399, H. sapiens 1439) while an uncovered one has single digits
            # (S. cerevisiae 9, and its run found zero TEs). 25 sits in the empty gap between those two
            # regimes; it is a warning trigger only — it never blocks a run or changes a result.
            if not custom_path and fams is not None and fams < self._ANNOT_LOW_LIBRARY:
                lineage = st.get("organism") or "this lineage"
                if library_mode == "dfam_uncurated":
                    warn.setText(f"The selected curated + uncurated Dfam search set contains only {fams} "
                                 f"family model(s) for {lineage}. Both family sets are already being searched, "
                                 "so a blank result would still be limited by this library's lineage coverage. "
                                 "Use a suitable organism-specific FASTA library before treating absence as "
                                 "evidence that the genome lacks transposable elements.")
                else:
                    warn.setText(f"The selected curated-only Dfam search contains only {fams} family model(s) "
                                 f"for {lineage}; uncurated families will not be searched. A blank result would "
                                 "therefore be inconclusive. Select “include uncurated families” (and install "
                                 "those partitions if needed), or use a suitable organism-specific FASTA library, "
                                 "before interpreting absence.")
                warn.show()
            start.setEnabled(True)
        except RuntimeError:
            pass

    def _start_genome_annotate(self, acc, organism, sensitivity, custom_lib=None, include_uncurated=False, sha256=None):
        if self._annot_inflight or not acc:
            return
        self._annot_inflight = True
        sha256 = sha256 or (getattr(self, "_annot_cost_state", None) or {}).get("sha256")
        # kept so 'Discard and start over' can start the run the user chose, not just discard the old one
        self._annot_request = {"acc": acc, "organism": organism, "sensitivity": sensitivity, "custom_lib": custom_lib,
                               "include_uncurated": bool(include_uncurated), "sha256": sha256}
        self._banner(f"Annotating {organism or acc} — this runs in the background. Its progress is shown at the top of "
                     "panel 06; finished chunks are kept if the run is interrupted.", "info")
        self._refresh_genome_manager()
        self._set_status_line(f"Annotating {organism or acc} — started; reading the first chunk. A mammalian genome "
                              "takes many hours. There is no in-app cancel yet.")
        self.card_genome.expand()                              # once, so the line is seen; later updates never re-open a card the user closed
        self.engine.submit("genome_annotate",
                           {"assemblyAccession": acc, "species": organism, "sensitivity": sensitivity,
                            # the thread count the budget probe MEASURED for this machine — the cost dialog
                            # showed it, so the run must actually use it rather than a fixed default
                            "threads": int(getattr(self, "_annot_threads", 0) or 4),
                            "chunk_mb": int(getattr(self, "_annot_chunk_mb", 0) or 40),
                            "custom_library": custom_lib, "include_uncurated": bool(include_uncurated),
                            "sha256": sha256},
                           key="genome_annotate")
        self._annot_timer = QTimer(self); self._annot_timer.setInterval(6000)
        self._annot_timer.timeout.connect(lambda: self.engine.submit("genome_annotate_log", {"tail": 1},
                                                                     key="genome_annotate_log"))
        self._annot_timer.start()
        # A fault ends the run through _reset_buttons, which does not know this bar: its own watcher hides the bar
        # once the run is no longer in flight, so a dead run never keeps reporting progress.
        self._annot_watch = QTimer(self); self._annot_watch.setInterval(2000)
        self._annot_watch.timeout.connect(self._annot_watch_tick)
        self._annot_watch.start()

    def _annot_watch_tick(self):
        if not self._annot_inflight:
            self._hide_annot_status()

    def _hide_annot_status(self):
        w = getattr(self, "_annot_watch", None)
        if w is not None:
            w.stop()
        bar = getattr(self, "_annotStatus", None)
        if bar is not None:
            try:
                bar.hide()
            except RuntimeError:
                pass

    def _on_genome_annotate_log(self, d):
        line = (d or {}).get("log") or ""
        if line and self._annot_inflight:
            org = (getattr(self, "_annot_request", None) or {}).get("organism") or "genome"
            self._set_status_line(f"Annotating {org} — {line}")

    def _set_status_line(self, text):
        """The persistent progress line of a running whole-genome annotation. It lives at the top of card 06, outside
        the panels that are cleared and rebuilt, so it is on screen for the hours the run takes."""
        bar = getattr(self, "_annotStatus", None)
        try:
            if bar is None:
                bar = self._annotStatus = BusyBar(text)
                self.card_genome.bodylay.insertWidget(0, bar)
            bar.set_text(text); bar.show()
        except RuntimeError:
            self._annotStatus = None

    def _on_genome_annotate(self, r):
        self._annot_inflight = False
        self._hide_annot_status()
        t = getattr(self, "_annot_timer", None)
        if t is not None:
            t.stop()
        # The backend reports a refused or failed run as {"ok": False, ...} — a normal RETURN, not an
        # exception — so it arrives here on the success signal. Without this guard the failure was
        # announced as a finished annotation and then crashed the results window on the missing fields.
        if not r or not r.get("ok"):
            msg = (r or {}).get("error") or "the annotation did not complete"
            self._refresh_genome_manager()
            if (r or {}).get("cache_corrupt"):
                self._banner("Genome annotation did not finish — " + msg, "warn")
                return
            # A settings clash is the one failure the user can clear, so offer the action instead of
            # telling them to do something the window gives them no way to do.
            if "different settings" in msg:
                req = getattr(self, "_annot_request", None) or {}
                acc = req.get("acc") or (getattr(self, "_annot_cost_state", None) or {}).get("accession")
                box = QMessageBox(self)
                box.setWindowTitle("Previous annotation found")
                box.setIcon(QMessageBox.Question)
                if (r or {}).get("annotation_state") == "complete":
                    box.setText("This genome already has a finished annotation that used different settings.")
                    box.setInformativeText(
                        "That finished result was made with the earlier settings and cannot stand in for the new ones. "
                        "Discarding it deletes the earlier result from this computer, and its per-hit alignment rows "
                        "cannot be exported first. Discard it and run again with the settings you just chose?")
                else:
                    box.setText("This genome has a part-finished annotation that used different settings.")
                    box.setInformativeText(
                        "Chunks already finished were searched with the earlier settings, so they cannot be "
                        "combined with new ones. Discard that unfinished work and start again with the "
                        "settings you just chose?")
                discard = box.addButton("Discard and start over", QMessageBox.AcceptRole)
                box.addButton("Keep it", QMessageBox.RejectRole)
                box.exec()
                if box.clickedButton() is discard and acc:
                    self._annot_restart = dict(req) if req else None   # the reset reply starts THIS run
                    self.engine.submit("genome_annotate_reset", {"assemblyAccession": acc},
                                       key="genome_annotate_reset")
                return
            self._banner("Genome annotation did not finish — " + msg, "warn")
            return
        # SUCCESS. This block used to live in _on_genome_annotate_reset, so the app's most expensive
        # operation — a multi-hour, uncancellable whole-genome run — finished with no banner and no results
        # window, while DISCARDING an annotation announced "transposable elements cover None% of the
        # assembly across None repeat models" and opened a results window on the reset reply's empty fields.
        # Each handler now reports its own outcome.
        self._annot_result = r
        self._refresh_genome_manager()
        if r.get("coverage_warning"):
            self._banner(r["coverage_warning"], "warn")
        else:
            self._banner(f"Genome annotation finished — transposable elements cover "
                         f"{r.get('te_percent')}% of the assembly across "
                         f"{r.get('te_repeat_name_count')} distinct repeat names/models.",
                         "success")
        self._show_annot_window(r)

    def _on_genome_annotate_reset(self, r):
        redo, self._annot_restart = getattr(self, "_annot_restart", None), None
        if r.get("ok"):
            self._refresh_genome_manager()
            if redo:                                      # 'Discard and start over' means exactly that
                self._banner("Previous annotation discarded. Starting the annotation you chose.", "info")
                return self._start_genome_annotate(redo["acc"], redo["organism"], redo["sensitivity"],
                                                   redo.get("custom_lib"), redo.get("include_uncurated", False),
                                                   redo.get("sha256"))
            self._banner("Previous annotation discarded. Open “Annotate TE landscape” to start again.",
                         "success")
        else:
            self._banner(r.get("error") or "Could not clear the previous annotation.", "warn")
        self._refresh_genome_manager()

    _ANNOT_COLS = [("Repeat name/model", "The matched RepeatMasker repeat name or consensus model (output column 10); "
                                        "this is not necessarily a biological family"),
                   ("Class/family", "RepeatMasker's broader class/family label. TEagle uses this label to "
                                    "separate TE, tandem and other coverage."),
                   ("Kind", "TE = transposable element; tandem = simple/low-complexity/satellite; other = "
                            "rRNA/tRNA; unclassified/other = no recognized class was available"),
                   ("Copies", "RepeatMasker alignment rows; one fragmented insertion can produce several rows"),
                   ("Bases", "Bases in this repeat name/model's within-name interval union. Different rows may "
                             "overlap, so rows must not be summed to recover a headline total."),
                   ("% genome", "Share of all submitted FASTA bases, including ambiguous bases"),
                   ("Mean divergence %", "Alignment-span-weighted mean RepeatMasker divergence from consensus. "
                                          "Raw alignment rows, including overlaps and fragments, supply the "
                                          "weights. It is not Kimura- or CpG-corrected, a mutation rate, or an age.")]

    @staticmethod
    def _annot_library_context(r):
        """Library-specific copy shared by the landscape window and its Markdown export."""
        kind = str(r.get("library_kind") or "")
        raw_custom = r.get("custom_library")
        custom = raw_custom if isinstance(raw_custom, dict) else {}
        is_custom = bool(custom) or kind.strip().lower() == "user" or kind.lower().startswith("user-supplied")
        rm_version = html.escape(str(r.get("repeatmasker_version") or "unresolved"))

        if is_custom:
            name = str(custom.get("name") or "unnamed user FASTA")
            sha = str(custom.get("sha256") or "unavailable")
            sequences = custom.get("sequences")
            sequence_md = f"; {sequences} reference sequences" if sequences is not None else ""
            sequence_ui = f"; {html.escape(str(sequences))} reference sequences" if sequences is not None else ""
            md_name = (name.replace("\\", "\\\\").replace("|", "\\|").replace("*", "\\*")
                       .replace("_", "\\_").replace("`", "\\`"))
            safe_name = html.escape(name)
            safe_sha = html.escape(sha)
            return {
                "report": (f"- Searched with RepeatMasker {rm_version} against the user-supplied FASTA "
                           f"{md_name}; SHA-256 `{sha}`{sequence_md}"),
                "display": (f"Searched with RepeatMasker {rm_version} against user-supplied FASTA "
                            f"<b>{safe_name}</b> (SHA-256 <code>{safe_sha}</code>{sequence_ui})"),
                "tooltip": (f"The run searched the user-supplied FASTA {safe_name}, identified by SHA-256 "
                            f"{safe_sha}. The reference-sequence count describes that file; it is not a count "
                             "of detected repeat names/models."),
                "scope": ("What this is: interspersed-repeat matches to references represented in the "
                          "supplied FASTA, placed by homology. It is not de novo repeat discovery. A reference "
                          "absent from that FASTA cannot be recovered through library homology, so absence "
                          "is not evidence of absence. RepeatMasker's separate simple/low-complexity screen can "
                          "still contribute non-TE rows. A custom reference contributes to the TE percentage only "
                          "when its FASTA header supplies a recognized class/family after '#', for example "
                          ">name#LTR/Gypsy; otherwise its matches remain visible as unclassified/other and are "
                          "excluded because classification is unavailable, not because they were shown to be non-TE. "
                          "Percentages use all submitted FASTA bases, including ambiguous bases and unplaced scaffolds."),
                "report_scope": ("- Interspersed-repeat matches to references represented in the supplied "
                                 "FASTA, placed by homology. This is not de novo repeat discovery. A reference "
                                 "absent from that FASTA cannot be recovered through library homology; "
                                 "RepeatMasker's separate simple/low-complexity screen can still contribute non-TE "
                                 "rows. Custom references without a recognized `#class/family` header remain "
                                 "visible as unclassified/other and do not enter the TE percentage because their "
                                 "classification is unavailable, not because they were demonstrated to be non-TE."),
            }

        library_kind = kind or "installed Dfam partitions"
        dfam_version = r.get("dfam_version")
        models = r.get("library_families_for_species")
        version_md = f" (Dfam {dfam_version})" if dfam_version else ""
        models_md = f"; {models} family models available for this lineage" if models is not None else ""
        safe_kind = html.escape(library_kind)
        version_ui = f" (Dfam {html.escape(str(dfam_version))})" if dfam_version else ""
        models_ui = f"; {html.escape(str(models))} family models available for this lineage" if models is not None else ""
        return {
            "report": (f"- Searched with RepeatMasker {rm_version} against {library_kind}{version_md}{models_md}"),
            "display": (f"Searched with RepeatMasker {rm_version} against {safe_kind}{version_ui}{models_ui}"),
            "tooltip": ("The exact tool and installed Dfam family set this result came from. The lineage model "
                        "count describes profiles available to the search, not families detected in this genome. "
                        "Chunks are assembly pieces used for progress and resumability."),
            "scope": ("What this is: matches to repeat names/models represented in the installed Dfam library, "
                      "placed by homology. It is not de novo repeat discovery. A model absent from the library "
                      "cannot appear here, so absence is not evidence of absence. Percentages use all submitted "
                      "FASTA bases, including ambiguous bases and unplaced scaffolds."),
            "report_scope": ("- Matches to repeat names/models represented in the searched Dfam library, placed "
                             "by homology. This is not de novo repeat discovery. A model absent from the library "
                             "cannot appear here; absence is not evidence of absence."),
        }

    @staticmethod
    def _annot_kind_label(row):
        """An unknown class is unavailable evidence, not a demonstrated non-TE classification."""
        kind = str(row.get("kind") or "")
        class_family = str(row.get("class_family") or "").strip().lower()
        if kind == "other" and class_family in ("", "unknown", "unclassified", "na", "n/a"):
            return "unclassified/other"
        return kind

    def _annot_row_menu(self, row):
        """Context menu for one genome-landscape row — aggregate-appropriate actions only."""
        kind = row.get("kind")
        items = [(f"Copy repeat name/model ({row.get('repeat_name')})", lambda: self._copy(str(row.get("repeat_name")))),
                 ("Copy this row (TSV)", lambda: self._copy("\t".join(
                     str(x) for x in (row.get("repeat_name"), row.get("class_family"), self._annot_kind_label(row), row.get("n"), row.get("bp"),
                                      row.get("percent"), row.get("divergence")))))]
        if kind != "TE":
            if self._annot_kind_label(row) == "unclassified/other":
                reason = ("This row has no recognized RepeatMasker class/family label. It is excluded from "
                          "TE coverage because its classification is unavailable, not because it was demonstrated "
                          "to be non-TE. For a custom library, use headers such as >name#LTR/Gypsy when supported.")
            else:
                reason = ("Simple repeats, low-complexity sequence and satellites are tandem repeats; rRNA and "
                          "tRNA entries are other repeat annotations. TEagle reports these categories separately "
                          "from transposable-element coverage.")
            items.append(("Why is this not counted as a TE?",
                          lambda _reason=reason: self._banner(_reason, "info")))
        else:
            # Same field as the "Mean divergence %" column, so it must carry the same hedges. This text
            # previously said "higher means older" and dropped the raw/uncorrected and saturation caveats
            # the column tooltip states — two different claims about one number, and the looser one was the
            # one a user reached by asking for an explanation.
            items.append(("What does mean divergence mean?",
                          lambda: self._banner(
                              "Alignment-span-weighted mean RepeatMasker divergence from the consensus. Raw "
                              "alignment rows, including overlapping and fragmented rows, supply the weights. "
                              "It is not Kimura- or CpG-corrected, a mutation rate, or an insertion-age estimate.",
                              "info")))
        return items

    def _show_annot_window(self, r):
        """The genome-wide result, in its own window: what was found, what was searched, and what that does
        and does not license the reader to conclude."""
        w = QDialog(self); w.setWindowTitle("Genome TE landscape"); w.setModal(False)
        w.setAttribute(Qt.WA_DeleteOnClose, True)
        lay = QVBoxLayout(w)
        lay.setContentsMargins(theme_mod.sp(14), theme_mod.sp(12), theme_mod.sp(14), theme_mod.sp(12))
        lay.setSpacing(theme_mod.sp(8))
        accession_stem = _safe_export_stem(r.get("accession"), "assembly", trust="validated_fetch")
        head = QLabel(f"{r.get('species')} · {r.get('accession')}")
        hf = head.font(); hf.setBold(True); hf.setPointSizeF(hf.pointSizeF() + 2); head.setFont(hf); lay.addWidget(head)
        # headline numbers: the TE share is stated separately from every repeat, because they are different
        # claims and merging them would overstate transposable-element content.
        tot = QLabel(f"<b>Transposable elements: {r.get('te_percent')}%</b> of "
                     f"{(r.get('genome_bp') or 0)/1e6:.1f} Mb, across "
                     f"{r.get('te_repeat_name_count')} distinct repeat names/models &nbsp;·&nbsp; "
                     f"all repeat classes together: {r.get('masked_percent')}% &nbsp;·&nbsp; "
                     f"{r.get('total_hits'):,} alignment rows")
        tot.setTextFormat(Qt.RichText); tot.setWordWrap(True)
        tot.setToolTip("Transposable elements = the share of the assembly covered by TE-classified matches (LTR, LINE, "
                       "SINE, DNA transposons, rolling-circle, retroposons, Penelope-like elements).\n"
                       "All repeat classes = the separate contig-wise union of every repeat row, including "
                       "tandem repeats (simple repeats, low-complexity sequence, satellites) and non-TE "
                       "entries such as rRNA and tRNA genes.\n"
                       "The categories are independently unioned and non-additive. Only the TE union is a statement "
                       "about transposable-element content. Each headline uses a separate cross-model interval union, "
                       "so a genomic base contributes at most once to that percentage; repeat-name/model rows use "
                       "within-name/model unions and therefore must not be summed. Percentages use all submitted "
                       "FASTA bases, including ambiguous bases.")
        lay.addWidget(tot)
        library = self._annot_library_context(r)
        prov = QLabel(f"{library['display']} · sensitivity {r.get('sensitivity')} · "
                      f"{r.get('chunks')} chunks · {r.get('elapsed_s')} s")
        prov.setObjectName("orient"); prov.setWordWrap(True)
        prov.setToolTip(library["tooltip"])
        lay.addWidget(prov)
        scope = QLabel(library["scope"])
        scope.setObjectName("orient"); scope.setWordWrap(True); lay.addWidget(scope)
        if r.get("coverage_warning"):
            cw = QLabel(r["coverage_warning"]); cw.setObjectName("errbanner"); cw.setWordWrap(True); lay.addWidget(cw)
        complete_note = ("The run completed all submitted assembly chunks. Library scope and sensitivity still "
                         "limit absence claims." if r.get("complete") else
                         "The run is incomplete: all counts, percentages, and rows are partial lower bounds. "
                         "No assembly-wide absence, coverage, or exhaustive-family claim is supported.")
        completion = QLabel(complete_note)
        completion.setObjectName("orient" if r.get("complete") else "errbanner")
        completion.setWordWrap(True); lay.addWidget(completion)
        repeat_names = r.get("repeat_names") or []
        t = DataTable([c[0] for c in self._ANNOT_COLS], dict(self._ANNOT_COLS))
        t.set_rows([[row["repeat_name"], row.get("class_family"), self._annot_kind_label(row), row["n"], row["bp"],
                     row.get("percent"), row["divergence"]] for row in repeat_names])
        t.export_notes = [library["scope"], complete_note,
                          "TE, tandem, other, and unclassified categories are independently unioned and non-additive. "
                          "Repeat-name/model rows use within-name unions and must not be summed to recover a headline. "
                          "Divergence is an alignment-span-weighted raw percentage, not an age estimate."]
        if r.get("coverage_warning"):
            t.export_notes.append(str(r["coverage_warning"]))
        # Row menu matched to what the row IS. A landscape row is a per-name/model AGGREGATE, not a locus, so
        # it offers no "design primer here" or "copy sequence" — there is no single sequence behind it.
        # A tandem/other row additionally says why it is not counted as a transposable element.
        t.set_row_menu(lambda r, _rows=repeat_names: self._annot_row_menu(_rows[r]) if r < len(_rows) else [])
        t.setMinimumHeight(round(260 * theme_mod.UI_SCALE)); lay.addWidget(t)
        # FigurePanel re-renders on theme/zoom, so it takes a builder (bg -> svg), not a rendered string.
        fig = FigurePanel(lambda bg, _r=r: figures.svg_te_composition(_r, theme=bg),
                          f"TEagle_TE_landscape_{accession_stem}")
        fig.setMinimumHeight(round(190 * theme_mod.UI_SCALE)); lay.addWidget(fig)
        row = QHBoxLayout(); row.addStretch(1)
        row.addWidget(_export_table_btn(t, f"TEagle_TE_landscape_{accession_stem}", self))
        # A genome run is sealed like every other result, so the seal has to be reachable from the result
        # — otherwise the run is reproducible in principle and not in practice. The report carries the
        # same hedges the screen shows, so an exported file cannot read as an unqualified claim.
        rep = QPushButton("Export report (.md)"); rep.setProperty("sm", True)
        rep.setToolTip("The summary table plus what was searched, what the numbers mean, and the limits "
                       "stated on this screen — the same hedges, in the file.")
        rep.clicked.connect(lambda: self._save_text(self._annot_report_md(r),
                                                    f"TEagle_TE_landscape_{accession_stem}", "md"))
        row.addWidget(rep)
        manb = QPushButton("Export manifest (.json)"); manb.setProperty("sm", True)
        manb.setToolTip("The sealed provenance record: tool and library versions, every parameter, the "
                        "genome checksum, and whether the run covered the whole assembly.")
        manb.setEnabled(bool(r.get("provenance")))
        manb.clicked.connect(lambda: self._export_manifest_payload(
            [r.get("provenance")], f"TEagle_TE_landscape_{accession_stem}_manifest"))
        row.addWidget(manb)
        close = QPushButton("Close"); close.setProperty("sm", True); close.clicked.connect(w.accept)
        row.addWidget(close); lay.addLayout(row)
        self._uppercase_buttons(w)
        w.resize(round(1000 * theme_mod.UI_SCALE), round(760 * theme_mod.UI_SCALE))
        w.show(); w.raise_()
        self._annot_window = w
