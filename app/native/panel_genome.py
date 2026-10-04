"""Cached-genome manager and whole-genome off-target scan."""
from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QHBoxLayout, QLabel, QMessageBox, QPushButton, QVBoxLayout, QWidget
import figures
from widgets import BusyBar, DataTable, FigurePanel
import theme as theme_mod
from ui_common import (
    CollapsibleCard, GLOSS, PANEL, _Combo, _amp_call, _amplicon_public_source, _clear_layout, _empty,
    _export_table_btn, _note, _sl)


class GenomeScanMixin:
    def _build_genome_card(self):
        # dedicated whole-genome off-target scan — LOCAL isPcr against a downloaded RefSeq assembly. The organism
        # dropdown lists ONLY downloaded+verified genomes (from genome_list). On-target = the product at the
        # specimen's own genome locus when it sits in the scanned assembly; the rest are off-target candidate sites.
        card = CollapsibleCard(PANEL["genome"], "Whole-genome off-target scan",
                               "optional · Linux backend: checking… · a designed pair vs a downloaded RefSeq genome (isPcr)")
        self.card_genome = card
        card.bodylay.addWidget(_sl("Organism — local isPcr against a downloaded genome"))
        gorow = QHBoxLayout()
        self.genomeOrg = _Combo()
        gorow.addWidget(self.genomeOrg, 1)
        self.genomeManageBtn = QPushButton("Manage genomes"); self.genomeManageBtn.setProperty("sm", True)
        self.genomeManageBtn.clicked.connect(self._open_genome_manager)
        gorow.addWidget(self.genomeManageBtn)
        card.bodylay.addLayout(gorow)
        self.genomeOrgHint = QLabel()
        self.genomeOrgHint.setObjectName("orient"); self.genomeOrgHint.setWordWrap(True)
        card.bodylay.addWidget(self.genomeOrgHint)
        # designed-pair picker + primary scan button — the discoverable in-card entry point (no right-click needed).
        # Routes to the SAME _scan_genome handler as the right-click, so the sealed isPcr job stays identical.
        srow = QHBoxLayout()
        self.scanPicker = _Combo(); srow.addWidget(self.scanPicker, 1)
        self.scanBtn = QPushButton("Scan whole genome"); self.scanBtn.setProperty("primary", True)
        self.scanBtn.clicked.connect(self._scan_from_picker); srow.addWidget(self.scanBtn)
        card.bodylay.addLayout(srow)
        self._refresh_genome_dropdown()                       # render from cached prepared set (empty at first build)
        gnote = QLabel("Off-target check: runs a local isPcr scan of a designed pair against a downloaded RefSeq "
                       "genome. Results are sealed with the assembly version + checksum. With a design locus, "
                       "intended-locus sites are listed before additional sites.")
        gnote.setObjectName("orient"); gnote.setWordWrap(True); card.bodylay.addWidget(gnote)
        self.genomeBody = QVBoxLayout(); gb = QWidget(); gb.setLayout(self.genomeBody)
        self.genomeBody.addWidget(_empty("Design a primer pair (panel 04) and download a genome, then scan here."))
        card.bodylay.addWidget(gb)
        return card

    _SCAN_PROD_MAX = 4000            # isPcr's product-size window for a genome scan; shown with every scan result
    _SCAN_ROWS_SHOWN = 500           # rows drawn in the scan table; an export always writes every product

    def _scan_genome(self, cand, org=None):
        """Whole-genome scan of one primer pair via LOCAL isPcr against a downloaded RefSeq assembly. Runtime
        depends on genome size, the installed isPcr build and cached RefSeq assembly, local hardware, and the
        fixed candidate-product-size cap. If the genome is not cached, the backend replies need_prepare and we
        offer a one-time download. Products are candidate priming sites, not validated bands; on/off-target
        labels are applied only when the specimen supplies a matching intended locus.
        `org` is captured at scan start so a later dropdown change can't retarget an in-flight scan/download."""
        self._clear_banner()
        if self._genome_inflight or self._genome_prep_inflight:
            return self._banner("A whole-genome scan / download is already running — wait for it to finish.", "info")
        run = getattr(self, "_pcr_run", None)
        if run and any(r is None for r in run.get("results", [])):
            return self._banner("An in-silico PCR run is still in progress — wait for it before a genome scan.", "info")
        if org is None:
            org = self.genomeOrg.currentData()
        if not org:
            return self._banner("Select an organism in the whole-genome off-target scan card, or download one first.", "info")
        self._pending_scan = {"cand": cand, "org": org}        # cand + org captured so a download resumes the SAME target
        self._scan_gen = self._specimen_gen                    # the specimen this scan was started for (see _scan_is_stale)
        self._genome_inflight = True
        self._refresh_genome_manager()                         # an open manager must lock its row actions while the scan runs
        self.card_genome.expand()
        # fixed wide window (isPcr defaults) — decoupled from the local-PCR size fields the user never sets for a
        # right-click genome scan. Widened only when the designed product itself is longer than the usual cap, so a
        # pair designed for a long amplicon is not reported absent; the window in force is sealed and shown.
        size = cand.get("product_size")
        prod_max = int(size * 1.25) if isinstance(size, int) and size > self._SCAN_PROD_MAX else self._SCAN_PROD_MAX
        loci = self._design_locus_for(cand, org)
        self._scan_ctx = {"prod_max": prod_max, "loci": loci or [], "designed_size": size}
        self._render_genome_busy(f"Scanning the {org} genome locally (isPcr) for priming sites. Runtime depends "
                                 "on genome size, the installed isPcr build and cached RefSeq assembly, and local "
                                 f"hardware. Only products up to {prod_max:,} bp are considered…")
        p = {"min_perfect": 15, "min_good": 15, "prod_max": prod_max, "prod_min": 0}
        self.engine.submit("genome_pcr", {"fwd": cand["left_seq"], "rev": cand["right_seq"], "organism": org,
                                          "design_locus": loci, "params": p}, key="genome_pcr")

    def _poll_prepare_log(self):
        if not self._genome_prep_inflight:
            self._prep_timer.stop(); return
        self.engine.submit("genome_prepare_log", {}, key="genome_prepare_log")

    def _on_genome_prepare_log(self, d):
        if not (self._genome_prep_inflight and d.get("log")):
            return
        org = getattr(self, "_prep_org", "")
        txt = f"Downloading the {org} genome — {d['log']}…"    # stream the milestone into the busy bar's caption
        b = getattr(self, "_genome_busy", None)
        try:
            if b is not None:
                b.set_text(txt)
            else:
                self._render_genome_busy(txt)
        except RuntimeError:                                  # busy bar was replaced (C++ object gone) — re-render
            self._render_genome_busy(txt)

    def _on_genome_prepare(self, d):
        self._genome_prep_inflight = False
        if getattr(self, "_prep_timer", None):
            self._prep_timer.stop()
        if not d.get("ok"):
            msg = "Genome download failed — " + d.get("error", "unknown error")
            level = "warn" if d.get("cache_corrupt") else "error"
            self._render_genome_status(msg, level)
            self._refresh_genome_manager()                    # download settled (failed): an open manager unlocks, so another genome can be deleted to free space
            return self._banner(msg, level)
        mb = (d.get("bytes", 0) or 0) / 1e6
        self._banner(f"Genome ready · {d.get('organism','')} ({d.get('assemblyAccession','')}) · "
                     f"{d.get('n_seqs','?')} sequences · {mb:.0f} MB cached.", level="success")
        ps = self._pending_scan
        if ps and self._scan_is_stale():
            self._pending_scan = None                         # the pair belonged to a specimen that has since changed
            self._render_genome_status(f"{d.get('organism','')} genome downloaded. The scan that was waiting for it "
                                       "belonged to a previous specimen and was not run.", "info")
        elif ps:                                              # resume straight into the SAME scan the download was for
            self._pending_scan = None
            self._scan_genome(ps["cand"], org=ps["org"])
        else:
            self._render_genome_status(f"{d.get('organism','')} genome downloaded — choose a pair above and press Scan whole genome.", "success")
        # refresh unconditionally: the just-downloaded organism must appear in the PCR dropdown even if the
        # manager is closed. The split _on_genome_list updates the dropdown always, the manager only if open.
        self.engine.submit("genome_list", {}, key="genome_list")

    # =================== cached-genome manager ===================
    def _scan_from_picker(self):
        """Run the whole-genome scan for the pair chosen in card 06 — the discoverable entry point. Routes to the
        SAME _scan_genome handler as the right-click, so the sealed isPcr job is byte-identical across entry points."""
        cand = self.scanPicker.currentData() if hasattr(self, "scanPicker") else None
        if not cand:
            return self._banner("Design a primer pair in panel 04 first, then pick it here to scan.", "warn")
        org = self.genomeOrg.currentData()
        if not org:
            return self._banner("Select a downloaded genome above (download one via Manage genomes).", "warn")
        self._scan_genome(cand, org=org)

    def _refresh_scan_picker(self):
        """Populate card 06's designed-pair picker from the current primer candidates, and enable/disable the scan
        button with an honest hint when there are no pairs or no downloaded genome."""
        picker = getattr(self, "scanPicker", None)
        if picker is None:
            return
        picker.blockSignals(True); picker.clear()
        cands = self.state.get("candidates") or []
        if not cands:
            picker.addItem("Design a primer pair in panel 04 first", None)   # self-describing empty state, not a blank combo
        for c in cands:
            picker.addItem(f"{self._pair_label(c)}  ·  {c['left_seq'][:12]}… / {c['right_seq'][:12]}…", c)
        picker.blockSignals(False)
        has_pairs = bool(cands)
        if hasattr(self, "pcrStageAll"):
            self.pcrStageAll.setEnabled(has_pairs)            # 'stage all' only ever stages the current designs
        has_genome = getattr(self, "genomeOrg", None) is not None and self.genomeOrg.count() > 1   # past the placeholder
        picker.setEnabled(has_pairs)
        self.scanBtn.setEnabled(has_pairs and has_genome)
        self.scanBtn.setToolTip("Design a primer pair in panel 04 first" if not has_pairs
                                else "Download a genome via Manage genomes first" if not has_genome else "")

    def _refresh_genome_manager(self):
        """Rebuild a LIVE (open) genome manager. Its row buttons are disabled while a scan/download runs;
        this re-enables them when the work settles. No-op when the manager is closed (never resurrects it)."""
        mgr = getattr(self, "_genome_mgr", None)
        if mgr is not None and mgr.isVisible():
            self.engine.submit("genome_list", {}, key="genome_list")

    def _open_genome_manager(self):
        """Show cached genomes (size, contigs) with delete + pre-download, so the user can manage disk."""
        self._genome_mgr_open = True                          # explicit open — distinguishes from a background refresh
        self.engine.submit("genome_list", {}, key="genome_list")

    def _add_custom_assembly(self):
        """Resolve a typed organism name / assembly accession ONCE, pin its versioned accession to the user store,
        then rebuild the manager + dropdowns so it is immediately scannable. The versioned accession is the seal
        anchor (a bare name can be promoted to a new RefSeq build over time); the name is a display label only."""
        if self._add_asm_inflight:                            # a rapid second click / Enter must not queue a duplicate resolve
            return
        edit = getattr(self, "_addAsmEdit", None)
        q = edit.text().strip() if edit is not None else ""
        if not q:
            return
        self._add_asm_inflight = True
        edit.setEnabled(False)                                # off-thread resolve: never block the UI (the app's
        btn = getattr(self, "_addAsmBtn", None)               # only network call that used to freeze the window)
        if btn is not None:
            try:
                btn.setEnabled(False)                         # the button, not just the field — a disabled field still returns .text()
            except RuntimeError:
                pass
        self._banner(f"Resolving “{q}” against NCBI…", "info")   # an in-progress status, not an error dialog
        self.engine.submit("add_custom_assembly", {"query": q}, key="add_custom_assembly")

    def _on_add_custom_assembly(self, entry):
        self._add_asm_inflight = False                        # cleared before the rebuild so the fresh Add control is enabled
        e = getattr(self, "_addAsmEdit", None)
        try:                                                  # the manager dialog may have been closed mid-resolve
            if e is not None:
                e.setEnabled(True); e.clear()
        except RuntimeError:
            pass
        self._banner(f"Added {entry['organism']} · {entry['assemblyName']} · {entry['assemblyAccession']} — "
                     "download it in the table below to scan.", "success")
        self.engine.submit("genome_list", {}, key="genome_list")   # rebuild manager + dropdowns with the new organism

    def _on_manager_download(self, org):
        """Download a genome from the manager table without closing the dialog; disable every row button
        immediately for feedback (the authoritative rebuild lands when _on_genome_prepare fires genome_list)."""
        self._prepare_genome(org)
        tbl = getattr(self, "_genome_mgr_table", None)
        if tbl is not None:
            for r in range(tbl.rowCount()):
                w = tbl.cellWidget(r, 6)
                if w is not None:
                    w.setEnabled(False)

    def _remove_genome(self, accession, org=""):
        # deleting the cache is destructive and costly to undo (rm -rf, then a fresh multi-minute NCBI
        # download) — confirm first, naming the organism and the re-download cost. Default = No.
        if self._genome_inflight or self._genome_prep_inflight or self._annot_inflight:
            return self._banner("A genome scan, download or annotation is running — wait for it to finish before "
                                "deleting a genome.", "warn")
        name = org or accession
        g = next((x for x in (getattr(self, "_prepared_genomes", None) or []) if x.get("accession") == accession), {})
        size = f"The cached search copy is {(g['bytes'] or 0) / 1e6:.0f} MB. " if g.get("bytes") else ""
        box = QMessageBox(self)
        box.setWindowTitle("Delete genome?")
        box.setText(f"Delete the cached {name} genome ({accession}) from this machine?\n\n{size}"
                    "Deleting also removes any whole-genome annotation saved with it, including every RepeatMasker chunk "
                    "already finished, and the genome FASTA an annotation cached, if there is one. An annotation would "
                    "then start again from the beginning.\n\n"
                    "Scanning it again needs the genome downloaded again from NCBI: several minutes, and at least 8 GB "
                    "free disk in WSL for any genome.")
        box.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
        box.setDefaultButton(QMessageBox.No)
        if box.exec() != QMessageBox.Yes:
            return
        self.engine.submit("genome_remove", {"assemblyAccession": accession}, key="genome_remove")

    def _on_genome_remove(self, d):
        # A refusal has to be visible. genome_remove can now decline (an annotation still holds the genome),
        # and this handler previously discarded its reply entirely — so a declined delete looked exactly like
        # a successful one: the list refreshed and the genome was simply still there.
        if isinstance(d, dict) and not d.get("ok"):
            self._banner(d.get("error") or "The genome could not be deleted.", level="warn")
        # refresh BOTH the dropdown and the manager (if open); the split _on_genome_list self-selects what to update
        self.engine.submit("genome_list", {}, key="genome_list")

    def _reset_genome_scan(self, why):
        """Drop a rendered whole-genome scan once the primers it describes are gone.

        The panel is only ever cleared by the scan paths themselves, so a verdict like "no off-target site
        in this assembly" survived a primer redesign and a new specimen — sitting under panel 06 as though
        it described whatever is on screen now. A stale specificity verdict is the worst kind of wrong
        answer here, because it reads as a clean result. No-op when nothing has been scanned.

        `_invalidate_analysis` drops the stored scan BEFORE it calls this, so the stored scan alone cannot say
        whether a verdict is still on screen; the panel's own flag does."""
        self._refresh_scan_picker()                       # the pairs on offer follow state['candidates'], which a new specimen drops
        if not (self.state.get("lastGenomeScan") or getattr(self, "_genome_shows_result", False)):
            return
        self.state["lastGenomeScan"] = None
        if getattr(self, "_genome_inflight", False) or getattr(self, "_genome_prep_inflight", False):
            return                                        # a running scan owns the panel; it will render its own result
        _clear_layout(self.genomeBody)
        self._genome_shows_result = False
        self.genomeBody.addWidget(_empty(why))

    def _render_genome_status(self, msg, level="error"):
        # every caller reports an OUTCOME (failed / cancelled / downloaded / not-yet-downloaded), never an
        # empty panel — so it gets the errbanner vocabulary at the caller's level, not the grey empty style
        _clear_layout(self.genomeBody)
        self._genome_shows_result = False
        self.genomeBody.addWidget(_note(msg, level))

    def _render_genome_busy(self, text):
        """Long-op liveness in the whole-genome scan panel: an animated indeterminate bar (download / scan),
        so a multi-minute WSL call visibly reads as working, not hung. set_text() updates it from the log poll."""
        _clear_layout(self.genomeBody)
        self._genome_shows_result = False
        self._genome_busy = BusyBar(text)
        self.genomeBody.addWidget(self._genome_busy)

    def _render_genome_scan(self, lanes, amps, summary):
        """Render a whole-genome off-target scan into its own panel: gel + interpretation (verdict,
        per-assembly-sequence spread, size cluster) + a FULL match table, ON-TARGET FIRST. On-target = the product
        at the specimen's own genome locus (when it sits in the scanned assembly); the rest are off-target candidate priming sites,
        or — with no design locus — neutral genomic priming sites."""
        _clear_layout(self.genomeBody)
        self._genome_shows_result = True                      # a verdict is now on screen (see _reset_genome_scan)
        _last = self.state.get("lastGenomeScan") or {}
        gel_base = "_".join(["TEagle_genome_gel"] + [self._pair_token(x) for x in (
            (_last.get("pair") or {}).get("id"), (_last.get("assembly") or {}).get("accession")) if x])
        gel = FigurePanel(lambda bg, L=lanes: figures.svg_gel({"lanes": L}, bg),
                          gel_base, modes=("dark", "white", "uv", "mono"),
                          hit_regions=figures.gel_regions({"lanes": lanes}), on_menu=self._gel_menu)
        # Name the pair and the assembly the verdict belongs to, ABOVE the verdict. Without it the panel
        # states a specificity conclusion with no subject, and a user who redesigns primers can read the
        # old sentence as describing the new ones.
        _last = self.state.get("lastGenomeScan") or {}
        _pair, _org = (_last.get("pair") or {}), _last.get("organism", "")
        if _pair.get("fwd"):
            _pid = _pair.get("id") or "pair"
            _sub = QLabel(f"Scanned pair <b>{_pid}</b>"
                          + (f" against <b>{_org}</b>" if _org else "")
                          + f"<br>F 5′-{_pair['fwd']}-3′<br>R 5′-{_pair.get('rev','')}-3′")
            _sub.setObjectName("orient"); _sub.setWordWrap(True); _sub.setTextFormat(Qt.RichText)
            self.genomeBody.addWidget(_sub)
        gel.apply_app_theme(self.theme); gel.setMinimumHeight(round(420 * theme_mod.UI_SCALE))
        self.genomeBody.addWidget(gel)
        verdict = QLabel(summary.get("verdict", ""))
        verdict.setObjectName("orient"); verdict.setWordWrap(True)
        vf = verdict.font(); vf.setBold(True); verdict.setFont(vf); self.genomeBody.addWidget(verdict)
        per = summary.get("per_source", [])
        if per:
            spread = " · ".join(f"{src}: {n}" for src, n in per[:10]) + (f" · … (+{len(per)-10} more)" if len(per) > 10 else "")
            lab = QLabel("per sequence — " + spread); lab.setObjectName("orient"); lab.setWordWrap(True); self.genomeBody.addWidget(lab)
        mode = summary.get("size_mode")
        if mode is not None:
            lo, hi = summary.get("size_min"), summary.get("size_max")
            rng = f"; range {lo}–{hi} bp" if lo != hi else ""
            sz = QLabel(f"product size clusters at ~{mode} bp ({summary.get('size_mode_n', 0)}/{summary.get('n_pair', 0)} pair hits){rng}")
            sz.setObjectName("orient"); sz.setWordWrap(True); self.genomeBody.addWidget(sz)
        has_locus = summary.get("has_locus")
        crit, asm, loci = _last.get("criteria") or {}, _last.get("assembly") or {}, _last.get("loci") or []
        cap = crit.get("max_size") or self._SCAN_PROD_MAX
        seqs = {a.get("source") for a in amps if not a.get("single_primer")}
        if has_locus and loci and seqs and not any(a.get("on_target") for a in amps) and \
                not seqs & {L.get("accession") for L in loci}:
            # the fetched locus is named on one accession (RefSeq) while the scanned FASTA may name the same
            # chromosome differently (a GenBank GCA assembly), so no product can ever match it by name
            self.genomeBody.addWidget(_note(
                f"The fetched locus lies on {loci[0].get('accession')}, but no product was found on a sequence of "
                "that name. This assembly may name its chromosomes differently (for example a GenBank GCA assembly "
                "against a RefSeq locus), so the intended product could not be matched and the on-target / "
                "off-target split below is not assessed. Treat every product as an unassessed priming site.", "warn"))
        def _call(a):
            return _amp_call(a, has_locus)              # same label the local-PCR table and the gel use
        order = sorted(range(len(amps)), key=lambda i: (amps[i].get("single_primer", False),
                                                        not amps[i].get("on_target", False),
                                                        amps[i]["source"], amps[i]["start"]))
        samps = [amps[i] for i in order]
        scan_table = None
        shown = min(len(samps), self._SCAN_ROWS_SHOWN)
        if samps:
            # isPcr coordinates are 1-based inclusive and are kept verbatim — the headers say so, so Len (end-start+1)
            # does not read as contradicting the span. 'Source' here is the assembly's own sequence name, not a specimen.
            # Start and End are separate numeric columns so an exported file loads as numbers, not 'a–b' text.
            gl = {**GLOSS, "Start (1-based)": "First base of this product in the scanned assembly, 1-based inclusive (isPcr's own).",
                  "End (1-based)": "Last base of this product, 1-based inclusive: Len = End − Start + 1."}
            t = DataTable(["Call", "Assembly seq", "Start (1-based)", "End (1-based)", "Len", "Strand"], gl)
            scan_table = t
            rows = [[_call(a), _amplicon_public_source(a), a["start"], a["end"], a["length"], a.get("strand", "?")]
                    for a in samps]
            if shown < len(rows):                       # a high-copy pair can hit the genome hundreds of thousands of times
                t.set_complete_rows(rows, shown)        # draw the first rows; an export still writes every product
            else:
                t.set_rows(rows)
            t.set_row_menu(lambda r: [("Copy locus (1-based)", lambda a=samps[r]: self._copy(
                f"{_amplicon_public_source(a)}:{a['start']}-{a['end']}({a.get('strand', '?')})"))])
            t.setMaximumHeight(round(260 * theme_mod.UI_SCALE))
            self.genomeBody.addWidget(t)
            pid = _pair.get("id")
            base = "_".join(["TEagle_offtarget_scan"] + [self._pair_token(x) for x in (pid, asm.get("accession")) if x])
            xrow = QHBoxLayout(); xrow.addStretch(1); xrow.addWidget(_export_table_btn(t, base, self))
            self.genomeBody.addLayout(xrow)
            if shown < len(samps):
                self.genomeBody.addWidget(_note(
                    f"Showing the first {shown:,} of {len(samps):,} products (on-target first, then by sequence and "
                    f"position). Sort the table or export it to reach the rest: an export lists all {len(samps):,}.", "info"))
        else:
            self.genomeBody.addWidget(_note("No genome-wide product for this pair under the criteria.", "warn"))
        n_on, n_off, n_single = summary.get("n_on", 0), summary.get("n_off", 0), summary.get("n_single", 0)
        head = f"{n_on} on-target + {n_off} off-target site(s)" if has_locus else f"{summary.get('n_pair', 0)} genomic priming site(s)"
        sp = f" · {n_single} single-primer artefact(s)" if n_single else ""
        comig = (" An on-target that shares a band size with an off-target is drawn in the off-target colour — a "
                 "co-migration cannot be resolved on a gel — with the full split kept in the table above." if has_locus else "")
        listing = "listed above (on-target first)" if has_locus else "listed above"
        export_listing = "listed on-target first" if has_locus else "listed"
        note = QLabel(f"{head}{sp}, {listing}.{comig} Every product is a candidate under isPcr's default ≥15 bp "
                      "minimum perfect 3′-end match requirement — a specificity screen, not wet-lab-validated bands. "
                      f"Only products up to {cap:,} bp are considered, so a longer product would not appear. The verdict is a "
                      "heuristic read of the count/spread; the numbers carry the claim.")
        note.setObjectName("orient"); note.setWordWrap(True); self.genomeBody.addWidget(note)
        if scan_table is not None:
            notes = [f"{head}{sp}, {export_listing}.{comig} Every product is a candidate under isPcr's "
                     "default >=15 bp minimum perfect 3-prime-end match requirement: this is a computational specificity screen, not wet-lab-validated "
                     "bands. The verdict is a heuristic interpretation of the displayed counts and spread."]
            if _pair.get("fwd"):
                where = (f"{_org}" if _org else "the selected genome") + (
                    f", assembly {asm['accession']}" + (f" ({asm['name']})" if asm.get("name") else "") if asm.get("accession") else "")
                notes.append(f"Primer pair {_pair.get('id') or 'pair'}: forward 5'-{_pair['fwd']}-3', reverse "
                             f"5'-{_pair.get('rev', '')}-3'. Scanned against {where}.")
            bits = [f"UCSC isPcr{' ' + str(crit['isPcr']) if crit.get('isPcr') else ''}"]
            if crit.get("min_perfect") is not None:
                bits.append(f"minimum perfect 3-prime match {crit['min_perfect']} bp")
            if crit.get("min_good") is not None:
                bits.append(f"minimum good match {crit['min_good']} bp")
            notes.append("Search: " + ", ".join(bits) + f", product sizes up to {cap:,} bp (longer products are not "
                         "reported). Coordinates are 1-based and inclusive, as isPcr reports them. The sealed run manifest "
                         "(assembly accession, file checksum, parameters) is exported from the Provenance card.")
            if shown < len(samps):
                notes.append(f"The on-screen table shows the first {shown:,} of {len(samps):,} products; this file lists all {len(samps):,}.")
            scan_table.export_notes = notes
        self._uppercase_buttons(self.genomeBody.parentWidget())   # this path has no post-render theme refresh
