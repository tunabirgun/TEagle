"""Dfam and RepeatMasker family annotation panel."""
from __future__ import annotations

import html
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QHBoxLayout, QLabel, QLineEdit, QPushButton, QTextEdit, QVBoxLayout, QWidget
import widgets
from widgets import BusyBar
import theme as theme_mod
from teagle_core.wsl import curated_coverage_sentence as _curated_coverage
from ui_common import (
    CollapsibleCard, PANEL, WSL_ORGANISMS, _Combo, _empty, _export_table_btn, _external, _kb_links, _note,
    _qq_note, _sl)


class FamilyMixin:
    # ---------- 03 Dfam / RepeatMasker family (WSL) ----------
    def _build_wsl_card(self):
        card = CollapsibleCard(PANEL["family"], "Dfam / RepeatMasker family",
                               "optional · Linux backend: checking… · names the TE family (Dfam)")
        self.card_wsl = card
        self.wslStatus = QLabel("checking WSL backend…"); self.wslStatus.setObjectName("cardmeta")
        self.wslStatus.setWordWrap(True); self.wslStatus.setTextFormat(Qt.RichText)
        card.bodylay.addWidget(self.wslStatus)
        srcrow = QHBoxLayout()
        srcrow.addWidget(QLabel("Sequence source"))
        self.wslSource = _Combo()
        self.wslSource.addItems(["Loaded specimen (panel 01)", "Paste sequence…"])
        srcrow.addWidget(self.wslSource); srcrow.addStretch(1)
        card.bodylay.addLayout(srcrow)
        self.wslPaste = QTextEdit(); self.wslPaste.setPlaceholderText("Paste a FASTA or raw DNA sequence to annotate against Dfam")
        self.wslPaste.setTabChangesFocus(True)
        self.wslPaste.setMaximumHeight(round(70 * theme_mod.UI_SCALE)); self.wslPaste.setVisible(False)
        self.wslSource.currentIndexChanged.connect(lambda i: self.wslPaste.setVisible(i == 1))
        card.bodylay.addWidget(self.wslPaste)
        row = QHBoxLayout()
        row.addWidget(QLabel("Organism"))
        self.wslSpecies = _Combo()                     # dropdown of common organisms + 'Other…'
        self.wslSpecies.addItem("— select organism —", None)
        for common, sci in WSL_ORGANISMS:
            self.wslSpecies.addItem(f"{common} · {sci}", sci)
        self.wslSpecies.addItem("Other…", "__other__")
        self.wslSpecies.currentIndexChanged.connect(lambda _i: self._on_species_changed())
        row.addWidget(self.wslSpecies, 1)
        self.wslSpeciesOther = QLineEdit(); self.wslSpeciesOther.setPlaceholderText("type organism / species")
        self.wslSpeciesOther.setVisible(False)
        row.addWidget(self.wslSpeciesOther, 1)
        # the panel's run action: primary, like RUN ANALYSIS / DESIGN PRIMERS / SCAN WHOLE GENOME
        self.annotateBtn = QPushButton("Run family annotation"); self.annotateBtn.setProperty("primary", True)
        self.annotateBtn.setEnabled(False); self.annotateBtn.clicked.connect(self._annotate)
        row.addWidget(self.annotateBtn)
        card.bodylay.addLayout(row)
        # Which families are searched. RepeatMasker reads the CURATED families only unless it is asked
        # for both, so without this the optional uncurated partitions sit on disk unreachable and a blank
        # result means "not searched" rather than "not present". Shown only once they are installed —
        # offering a library the machine does not have would promise a search it cannot run.
        librow = QHBoxLayout()
        librow.addWidget(QLabel("Library"))
        self.wslLibrary = _Combo()
        self.wslLibrary.addItem("Curated families only", False)
        self.wslLibrary.addItem("Include uncurated families", True)
        self.wslLibrary.setToolTip(
            "Which Dfam families RepeatMasker searches. Curated families are reviewed but cover few "
            "lineages deeply; outside a handful of intensively studied species almost every family is "
            "uncurated, so a curated-only search can return nothing for a real element. Including the "
            "uncurated families searches both and is recorded with the result.")
        librow.addWidget(self.wslLibrary, 1); librow.addStretch(1)
        self.wslLibraryRow = QWidget(); self.wslLibraryRow.setLayout(librow)
        self.wslLibraryRow.setVisible(False)      # revealed by _on_wsl_status once the partitions are there
        card.bodylay.addWidget(self.wslLibraryRow)
        self.wslInstallBtn = QPushButton("Backend installer — install · repair · check integrity")
        self.wslInstallBtn.setProperty("sm", True)
        self.wslInstallBtn.clicked.connect(self._open_installer)
        irow = QHBoxLayout(); irow.addWidget(self.wslInstallBtn); irow.addStretch(1)   # secondary: intrinsic width, not card-wide
        card.bodylay.addLayout(irow)
        self.wslBody = QVBoxLayout(); wb = QWidget(); wb.setLayout(self.wslBody)
        self.wslBody.addWidget(_empty("Run RepeatMasker against Dfam to name the TE family. Family naming is the Linux (WSL) backend."))
        card.bodylay.addWidget(wb)
        return card

    # =================== WSL family annotation ===================
    def _init_wsl(self):
        self.engine.submit("wsl_status", key="wsl_status")

    @staticmethod
    def _backend_states(w):
        """Plain state of the OPTIONAL Linux backend from a wsl_status payload (None = no answer yet). 'ready' needs
        both RepeatMasker + Dfam and minimap2; the whole-genome scan is never reported as ready here, because the
        status probe does not test isPcr."""
        if w is None:
            c = "checking…"
            return {"summary": c, "rm": c, "mm": c, "wsl": c}
        no_distro = w.get("wsl2") is False and not w.get("distro")
        if w.get("error") and not no_distro:
            u = "status unknown"
            return {"summary": u, "rm": u, "mm": u, "wsl": u}
        if not w.get("wsl2"):
            n = "not installed"
            return {"summary": n, "rm": n, "mm": n, "wsl": n}
        rm = "ready" if w.get("ready") else "not set up"
        mm = "ready" if w.get("minimap2") else "not set up"
        return {"summary": rm if rm == mm else "partly installed", "rm": rm, "mm": mm, "wsl": "WSL2 found"}

    def _render_status_chip(self):
        """Compose the header chip from the two independent states, whichever arrived first."""
        core, st = self._core, self._backend_states(self._wsl)
        ok = bool(core and core.get("ok"))
        self.statusTxt.setText("connecting…" if core is None else
                               f"engine ready · primer3 {core.get('primer3')}" if ok else "engine unavailable")
        self.backendTxt.setText(f"Linux backend: {st['summary']}")
        tip = ((f"Core engine (in-process; primer3 {core.get('primer3')}): classification, ORFs, protein domains, primer "
                "design and in-silico PCR — always available." if ok else
                "Core engine: not available. Classification, ORFs, protein domains, primer design and in-silico PCR "
                "need it." if core is not None else "Core engine: starting.")
               + f"\nOptional Linux (WSL) backend: family naming (panel {PANEL['family']}) {'ready' if st['rm'] == 'ready' else 'not ready'}; "
               f"splice structure (panel {PANEL['splice']}) {'ready' if st['mm'] == 'ready' else 'not ready'}; whole-genome scan "
               f"(panel {PANEL['genome']}) needs isPcr, which is checked when a scan runs.")
        for w in (self._chip, self.statusTxt, self.backendTxt):
            w.setToolTip(tip)
        self.led.setProperty("live", ok); self._repolish(self.led)

    def _tag_backend_cards(self):
        """Panels that need the optional Linux backend say so in their header, with its live state."""
        st = self._backend_states(self._wsl)
        tag = "optional · Linux backend: "
        self.card_wsl.set_meta(f"{tag}{st['rm']} · names the TE family (Dfam)")
        self.card_splice.set_meta(f"{tag}{st['mm']} · exon–intron structure from an aligned transcript")
        self.card_genome.set_meta(f"{tag}{st['wsl']} · a designed pair vs a downloaded RefSeq genome (isPcr)")

    def _on_wsl_status(self, w):
        self._wsl = w
        self._render_status_chip()
        self._tag_backend_cards()
        self._wsl_ready = {"annotate": False, "splice": False}
        # Every probe starts from "nothing is available" and enables only what THIS status supports: a re-probe
        # after a clean reinstall (or a lost WSL) must not leave the run buttons armed or a stale library choice
        # selectable. A run already in flight keeps its own busy state.
        for kind, btn in (("annotate", self.annotateBtn), ("splice", self.spliceBtn)):
            if kind not in getattr(self, "_reqs", {}):
                btn.setEnabled(False)
        self._uncurated_available = False
        self.wslLibraryRow.setVisible(False)
        self._wsl_probed = True                           # the flags below now come from a real probe, not a default
        self._wsl_partial_uncurated = (not w.get("dfam_uncurated") and any(
            "uncurated" in str(p).lower() for p in ((w.get("dfam_library") or {}).get("partitions") or [])))
        # wsl.available() reports a MISSING distro as wsl2=False with an 'error' text; that is the install
        # guidance case, not a failed probe. A distro that exists but did not answer is the error case.
        no_distro = w.get("wsl2") is False and not w.get("distro")
        unreadable = bool(w.get("error")) and not no_distro
        stack = ("status could not be read" if unreadable else "not installed" if not w.get("wsl2") else
                 "ready" if w.get("ready") else "installed, annotation stack incomplete")
        self.backendBtn.setToolTip("Linux backend — install, repair and check the optional Linux (WSL) stack: Dfam "
                                   f"family naming, splice structure, whole-genome scan. Optional stack: {stack}.")
        if unreadable:
            self.wslStatus.setText(f"<span style='color:{theme_mod.BAD[self.theme]}'>WSL status error: {html.escape(str(w['error']))}</span>")
            self.spliceStatus.setText("WSL status could not be read, so splice structure is unavailable. Open "
                                      "<b>Backend</b> in the header to check the Linux stack.")
            return
        if not w.get("wsl2"):
            found = w.get("found") or []
            self.wslStatus.setText(("<b>WSL2 has no usable Linux distribution</b>" + f" (found: {html.escape(', '.join(found))})"
                                    if found else "<b>WSL2 not installed</b>")
                                   + " — this optional step names the Dfam family. "
                                   "The domain-based superfamily above works without it. Install it in one click with "
                                   "<b>Backend installer</b> below (it runs the elevated <code>wsl --install</code> for you; "
                                   "a Windows restart may be required), or run <code>wsl --install</code> in an Administrator PowerShell.")
            self.spliceStatus.setText("<b>WSL2 not installed</b> — splice-structure alignment needs it (optional).")
            return
        self.engine.submit("genome_list", {}, key="genome_list")   # WSL is up — populate the downloaded-genome dropdown
        # offer the library choice only once the uncurated partitions are actually on the machine —
        # a curated-only backend has nothing extra to search, so the option would promise a search it
        # cannot run. Read from the library famdb reports, not from an assumption about the install.
        # tracked explicitly, not read back off isVisible(): a widget inside a collapsed card reports
        # not-visible, and reading the choice from that would silently discard what the user selected
        # The backend sets this only when BOTH optional Dfam shards passed its file probe. Looking for an
        # "uncurated" partition name here would expose a search choice after a partial install that the
        # backend must reject.
        self._uncurated_available = bool(w.get("dfam_uncurated"))
        self.wslLibraryRow.setVisible(self._uncurated_available)
        if w.get("ready"):
            # name the partitions actually on the machine — "Dfam curated" was fixed text and read as a
            # statement about the library even where the uncurated partitions were installed alongside it
            _lib = "Dfam curated + uncurated available" if self._uncurated_available else "Dfam curated"
            self.wslStatus.setText(f"<span style='color:{theme_mod.GOOD[self.theme]}'>● ready</span> · RepeatMasker {w.get('repeatmasker')} "
                                   f"· {_lib} · distro {w.get('distro')}")
            self.annotateBtn.setEnabled(True); self._wsl_ready["annotate"] = True
        else:
            self.wslStatus.setText(f"WSL2 ok ({w.get('distro')}); annotation stack not installed "
                                   f"(RepeatMasker {w.get('repeatmasker') or 'missing'}, Dfam {'ok' if w.get('dfam') else 'missing'}).")
            self.wslInstallBtn.setVisible(True)
        if w.get("minimap2"):
            self.spliceStatus.setText(f"<span style='color:{theme_mod.GOOD[self.theme]}'>● ready</span> · minimap2 {w.get('minimap2')} "
                                      "· align a transcript to resolve exon–intron structure")
            self.spliceBtn.setEnabled(True); self._wsl_ready["splice"] = True
        elif w.get("wsl2"):
            self.spliceStatus.setText("minimap2 not installed in the WSL backend — it ships with the managed install (Backend in the header).")

    def _on_species_changed(self):
        """Show the free-text field only for 'Other…'; clear it when switching to a listed organism
        so a stale lineage can't silently drive the next annotation."""
        other = self.wslSpecies.currentData() == "__other__"
        self.wslSpeciesOther.setVisible(other)
        if not other:
            self.wslSpeciesOther.clear()
        if not getattr(self, "_species_busy", False):
            self._species_auto = False                        # the user chose it: it outlives a specimen change

    def _species(self):
        """Selected organism for RepeatMasker (-species), or None. 'Other…' uses the free-text field."""
        if self.wslSpecies.currentData() == "__other__":
            return self.wslSpeciesOther.text().strip() or None
        return self.wslSpecies.currentData() or None

    def _set_species(self, name):
        """Auto-select the dropdown entry matching a fetched organism; fall back to 'Other…' free-text."""
        name = (name or "").strip()
        if not name:
            return
        self._species_busy = True
        try:
            for i in range(self.wslSpecies.count()):
                d = self.wslSpecies.itemData(i)
                if isinstance(d, str) and d != "__other__" and d.lower() == name.lower():
                    self.wslSpecies.setCurrentIndex(i); break
            else:
                oi = self.wslSpecies.findData("__other__")    # not a listed organism -> Other + free text
                if oi >= 0:
                    self.wslSpecies.setCurrentIndex(oi)
                    self.wslSpeciesOther.setText(name)
        finally:
            self._species_busy = False
        self._species_auto = True                             # taken from the specimen, so the next specimen replaces it

    def _drop_specimen_species(self):
        """The organism came from the previous specimen (fetch or bundled example): clear it rather than run
        RepeatMasker with a lineage that does not describe the new sequence. A choice the user made stays."""
        if getattr(self, "_species_auto", False) and hasattr(self, "wslSpecies"):
            self._species_busy = True
            try:
                self.wslSpecies.setCurrentIndex(0)
            finally:
                self._species_busy = False
            self._species_auto = False

    def _annotate(self):
        self._clear_banner()                                  # drop any stale error before a retry (matches _design/_run_pcr)
        if self.wslSource.currentIndex() == 1:
            seq = self.wslPaste.toPlainText().strip()
            src = None
        else:
            seq = self.state.get("seq") or self.seq.toPlainText().strip()
            seq, src = self._operation_input(seq)
        if not seq:
            return self._banner("Load or paste a specimen first — there is no sequence to annotate.", "info")
        self.state["family_seq"] = self._norm_seq(seq)    # hit coords index THIS sequence (backend-normalized), not always panel-01
        self.annotateBtn.setEnabled(False); widgets.set_caps(self.annotateBtn, "◴ annotating…")
        self._set_body(self.wslBody, BusyBar("Running RepeatMasker against Dfam — this can take a minute or two…"))
        self.engine.submit("annotate", {"sequence": seq, "species": self._species(), "source": src,
                                        "include_uncurated": bool(self.wslLibrary.currentData())
                                        if getattr(self, "_uncurated_available", False) else False},
                           key=self._issue("annotate", seq, bound=self.wslSource.currentIndex() != 1))

    def _on_annotate(self, d):
        self.annotateBtn.setEnabled(True); widgets.set_caps(self.annotateBtn, "Run family annotation")
        self.card_wsl.expand()
        if not d.get("ok"):
            self._set_body(self.wslBody, _note(d.get("error", "annotation failed")))
            return
        self._render_family(d)
        if d.get("provenance"):
            self._render_provenance(d["provenance"])

    # A satellite or simple repeat is a POSITIVE annotation — RepeatMasker identified what the locus is,
    # and the answer is "not a transposable element". "Unknown" is a different statement: a repeat whose
    # class could not be assigned. Neither is an absence of result, so neither is silently dropped.
    _NONTE = {"Low_complexity", "Simple_repeat", "Satellite"}
    _UNCLASSED = {"Unknown", "Unspecified"}

    # The family tables are RepeatMasker output, not the HMMER domain panel the shared glossary describes.
    _FAMILY_GLOSS = {
        "Score": "RepeatMasker's Smith-Waterman alignment score (RMBlast) for this match. It is not a HMMER bit "
                 "score and cannot be compared with the domain table's Score.",
        "Class/family": "RepeatMasker class/family label of the matched Dfam entry, e.g. LTR/Copia. It also covers "
                        "non-TE classes (satellites, simple repeats, RNA genes), which are listed apart.",
        "Div (raw)": "Raw % substitutions between the query and the Dfam consensus over the aligned span; not "
                     "Kimura-corrected. Lower is closer.",
        "Del": "% of consensus bases missing from the query inside the alignment (RepeatMasker %del).",
        "Ins": "% of query bases with no counterpart in the consensus inside the alignment (RepeatMasker %ins).",
        "Consensus": "Aligned span on the Dfam consensus (begin–end) over the consensus length, in bp.",
        "Cov": "Percent of the Dfam consensus covered by the alignment (summed over merged blocks).",
        "Blocks": "Number of RepeatMasker alignment rows merged into this row because they share a fragment ID "
                  "(one interrupted element).",
    }

    def _repeat_table(self, hits):
        row = lambda i, h, b: [i + 1, h["class_family"], h["family"], self._span(h["q_start"], h["q_end"], b),
                               h["strand"], f"{h['divergence']}%", h["score"]]
        t = self._coord_table(["#", "Class/family", "Dfam family", self._ch("Coords"), "Str", "Div", "Score"],
                              self._gloss(self._FAMILY_GLOSS), [row(i, h, 0) for i, h in enumerate(hits)])
        t.export_header_map = {"Coords (1-based)": "Coords (0-based)"}
        t.set_rows([row(i, h, self.coord_base) for i, h in enumerate(hits)])
        return t

    @classmethod
    def _partition_hits(cls, hits):
        """(te, tandem, unclassified, other) by the backend's own classifier (wsl.repeat_kind), so this panel and
        the whole-genome landscape agree. 'other' holds RNA genes, artefacts and other non-TE annotations."""
        from teagle_core.wsl import repeat_kind
        te, tandem, unclassed, other = [], [], [], []
        for h in hits:
            kind = repeat_kind(h["class_family"])
            if kind == "TE":
                te.append(h)
            elif kind == "tandem":
                tandem.append(h)
            elif h["class_family"].split("/")[0] in cls._UNCLASSED:
                unclassed.append(h)
            else:
                other.append(h)
        return te, tandem, unclassed, other

    def _render_family(self, d):
        self.state["familyResult"] = d                        # kept so the coordinate display can redraw it
        hits = d.get("hits", [])
        te, nonte, unclassed, other = self._partition_hits(hits)
        if not te:
            cont = QWidget(); cl = QVBoxLayout(cont); cl.setContentsMargins(0, 0, 0, 0); cl.setSpacing(6)
            if nonte:
                kinds = sorted({h["class_family"].replace("_", " ").lower() for h in nonte})
                nonte_note = (f"RepeatMasker annotated {len(nonte)} region(s) of this locus as "
                              f"{' and '.join(kinds)}. That is a result, not a failure: those regions match a "
                              f"tandem/low-complexity repeat, and no transposable-element family was found in "
                              f"them.")
                cl.addWidget(_note(nonte_note, "info"))
                nonte_table = self._repeat_table(nonte)
                nonte_table.export_notes = [nonte_note]
                cl.addWidget(nonte_table)
                nonte_row = QHBoxLayout(); nonte_row.addStretch(1)
                nonte_row.addWidget(_export_table_btn(nonte_table, "TEagle_non_TE_repeats", self))
                cl.addLayout(nonte_row)
            if unclassed:
                unclassed_note = (f"{len(unclassed)} repeat region(s) matched a Dfam entry whose class is "
                                  f"unassigned ('Unknown'). A repeat is present; its type is not established.")
                cl.addWidget(_note(unclassed_note, "info"))
                unclassed_table = self._repeat_table(unclassed)
                unclassed_table.export_notes = [unclassed_note]
                cl.addWidget(unclassed_table)
                unclassed_row = QHBoxLayout(); unclassed_row.addStretch(1)
                unclassed_row.addWidget(_export_table_btn(unclassed_table, "TEagle_unclassified_repeats", self))
                cl.addLayout(unclassed_row)
            if other:
                other_note = (f"{len(other)} region(s) of this locus were annotated as RNA genes or other non-TE "
                              "features (for example snRNA, tRNA or an artefact). They are listed apart and are "
                              "not transposable-element family calls.")
                cl.addWidget(_note(other_note, "info"))
                other_table = self._repeat_table(other)
                other_table.export_notes = [other_note]
                cl.addWidget(other_table)
                other_row = QHBoxLayout(); other_row.addStretch(1)
                other_row.addWidget(_export_table_btn(other_table, "TEagle_other_annotations", self))
                cl.addLayout(other_row)
            if True:                                      # the blank-result caveats apply whenever no TE family matched
                # Two independent conditions decide whether a family CAN be named, and the message names
                # whichever is actually missing rather than guessing. Measured on Drosophila copia: with
                # the uncurated partitions installed AND a species it resolves to Copia_LTR and Copia_I at
                # 100% consensus coverage; with the partitions but NO species it returns only
                # low-complexity, because RepeatMasker searches a limited default set without a lineage.
                # How much the curated-only case costs depends entirely on the lineage, so the message
                # below no longer names a figure for one organism — see _CURATED_COVERAGE.
                parts = ((d.get("dfam_library") or {}).get("partitions")) or []
                searched_unc = bool(d.get("include_uncurated"))
                sp = d.get("species") or ""
                defaulted = d.get("species_defaulted")
                has_sp = (not defaulted) if defaulted is not None else (bool(sp) and not str(sp).startswith("(all"))
                if getattr(self, "_wsl_probed", False):   # what the backend probe established, not partition names
                    has_unc, partial_unc = bool(self._uncurated_available), bool(self._wsl_partial_uncurated)
                else:
                    has_unc, partial_unc = any("uncurated" in str(p).lower() for p in parts), False
                blank = "No Dfam family matched" if not hits else "No transposable-element family matched"
                if not has_sp:
                    lineage = d.get("lineage_searched") or "Homo sapiens"
                    cl.addWidget(_note(f"{blank} — and no organism was set, so RepeatMasker searched its default "
                                       f"lineage ({lineage}). It searches only a limited set without the right "
                                       "lineage, so an organism is usually what turns a blank result into a named "
                                       "family. Set the organism above and run again.", "warn"))
                elif not searched_unc and partial_unc:
                    cl.addWidget(_note(f"{blank} for “{sp}” — this run searched the CURATED families only. Only one of "
                                       "the two optional uncurated Dfam partitions is present, so the Library choice "
                                       "is not offered: one partition is missing. Repair it in BACKEND in the "
                                       "header, then run again with <b>Library</b> set to include uncurated "
                                       "families.", "warn"))
                elif not searched_unc and not has_unc:
                    cl.addWidget(_note(f"{blank} for “{sp}” — and only the CURATED Dfam "
                                       "partitions are installed. How much that costs depends entirely on "
                                       "the lineage: " + _curated_coverage() + " A blank result here may "
                                       "therefore be a limit of the installed library rather than a "
                                       "property of the sequence. Install the optional uncurated "
                                       "partitions from BACKEND in the header, then run again with "
                                       "<b>Library</b> set to include them.", "warn"))
                elif not searched_unc:
                    cl.addWidget(_note(f"{blank} for “{sp}” — but this run searched the "
                                       "CURATED families only. The uncurated partitions are installed on "
                                       "this machine and were not read: RepeatMasker searches curated "
                                       "families unless it is asked for both. Set <b>Library</b> to "
                                       "“Include uncurated families” and run again before concluding "
                                       "anything from this blank result.", "warn"))
                else:
                    cl.addWidget(_note(f"{blank} in this locus for “{sp}”. The curated and "
                                       "uncurated partitions were both searched with a lineage set, so "
                                       "this is a genuine no-match rather than a gap in the installed "
                                       "library. The structural and protein-domain evidence above does "
                                       "not depend on this search.", "warn"))
            self._set_body(self.wslBody, cont)
            return
        self.state["family"] = te
        # Name the family set actually SEARCHED, which is decided by the run, not by what is on disk.
        # RepeatMasker reads the curated families unless it is asked for both, so a machine with the
        # uncurated partitions installed still searches curated-only by default — labelling such a result
        # "curated + uncurated" would credit it with a search it did not perform. Conversely a hit named
        # from the uncurated set must not be reported as curated: curated means manually reviewed,
        # uncurated means automatically built and less vetted.
        _dfam_lbl = ("Dfam 4.0 curated + uncurated" if d.get("include_uncurated")
                     else "Dfam 4.0 curated only")
        # the header names the class of the strongest TE match, not of the leftmost one
        top = max(te, key=lambda h: h.get("score") or 0)
        extra = len({h["class_family"] for h in te} - {top["class_family"]})
        species = (f"not set — RepeatMasker searched its default lineage ({d.get('lineage_searched') or 'Homo sapiens'})"
                   if d.get("species_defaulted") else html.escape(str(d.get("species", ""))))
        head = QLabel(f"<b>Dfam · {top['class_family']}</b>"
                      + (f" (strongest match; {extra} other TE class{'es' if extra > 1 else ''} below)" if extra else "")
                      + f" — {' · '.join(sorted({h['family'] for h in te}))} "
                      f"· {_dfam_lbl}{self._src_html('Dfam')} · RepeatMasker {d.get('repeatmasker_version','')}"
                      f"{self._src_html('RepeatMasker')} · species: {species}")
        head.setTextFormat(Qt.RichText); head.setWordWrap(True); _external(head); _kb_links(head)
        # %div/%del/%ins and the consensus-side coverage are RepeatMasker's own alignment evidence; they
        # were parsed and discarded before. Divergence is RAW (not Kimura-corrected) and the header says so.
        headers = ["#", "Class/family", "Dfam family", self._ch("Coords"), "Str",
                   "Div (raw)", "Del", "Ins", "Consensus", "Cov", "Blocks", "Score"]

        def _cons(h):
            if h.get("cons_start") is None or h.get("cons_length") is None:
                return "—"
            return f"{h['cons_start']}–{h['cons_end']} / {h['cons_length']}"

        fam_row = lambda i, h, b: [i + 1, h["class_family"], h["family"], self._span(h["q_start"], h["q_end"], b), h["strand"],
                                   f"{h['divergence']}%",
                                   ("—" if h.get("pct_del") is None else f"{h['pct_del']}%"),
                                   ("—" if h.get("pct_ins") is None else f"{h['pct_ins']}%"),
                                   _cons(h),
                                   ("—" if h.get("cons_coverage_pct") is None else f"{h['cons_coverage_pct']}%"),
                                   h.get("n_fragments", 1), h["score"]]
        t = self._coord_table(headers, self._gloss(self._FAMILY_GLOSS), [fam_row(i, h, 0) for i, h in enumerate(te)])
        t.export_header_map = {"Coords (1-based)": "Coords (0-based)", "Consensus": "Consensus (1-based, inclusive)"}
        t.set_rows([fam_row(i, h, self.coord_base) for i, h in enumerate(te)])
        family_notes = [
            f"Searched {_dfam_lbl} with RepeatMasker {d.get('repeatmasker_version','')} for species "
            f"{d.get('species','') or 'not set'}.",
            "Search ran in RepeatMasker's rush mode (-qq): about 10% less sensitive than a default search and 4 to 10 "
            "times faster (RepeatMasker's own description). Divergence is the raw substitution percentage, not Kimura-corrected; consensus coverage and "
            "divergence are alignment evidence, not a family-membership test.",
        ]
        t.export_notes = family_notes
        t.set_row_menu(lambda r: self._feat_menu(te[r]["q_start"], te[r]["q_end"], te[r]["strand"],
                                                 te[r]["family"], src_seq=self.state.get("family_seq")))
        cont = QWidget(); cl = QVBoxLayout(cont); cl.setContentsMargins(0, 0, 0, 0); cl.setSpacing(6)
        cl.addWidget(head); cl.addWidget(t)
        merged = [h for h in te if h.get("n_fragments", 1) > 1]
        if merged:
            merged_note = ("; ".join(h["fragment_note"] for h in merged) +
                           ". Merged so an interrupted element is one row, with the outermost query span "
                           "and the summed consensus coverage.")
            cl.addWidget(_note(merged_note, "info"))
            t.export_notes.append(merged_note)
        cl.addWidget(_qq_note())
        if nonte or unclassed or other:            # reported, never silently dropped
            also = [f"{len(nonte)} tandem / low-complexity region(s)"] if nonte else []
            also += [f"{len(unclassed)} repeat(s) of unassigned class"] if unclassed else []
            also += [f"{len(other)} RNA-gene or other non-TE region(s)"] if other else []
            cl.addWidget(_note("Also annotated in this locus: " + " and ".join(also) +
                               ", listed apart from the TE families above.", "info"))
        def add_partition(partition, title, base, note_text):
            if not partition:
                return
            cl.addWidget(_sl(title))
            pt = self._repeat_table(partition)
            pt.export_notes = [note_text] + family_notes
            cl.addWidget(pt)
            prow = QHBoxLayout(); prow.addStretch(1)
            prow.addWidget(_export_table_btn(pt, base, self)); cl.addLayout(prow)

        add_partition(nonte, "Tandem / low-complexity repeats", "TEagle_non_TE_repeats",
                      "These rows are positive tandem/low-complexity repeat annotations reported separately; "
                      "they are not transposable-element family calls.")
        add_partition(unclassed, "Unclassified repeats", "TEagle_unclassified_repeats",
                      "These rows matched repeat entries whose class is unassigned. A repeat is present, but its "
                      "type is not established; exclusion from TE coverage is not evidence that it is non-TE.")
        add_partition(other, "RNA genes and other non-TE annotations", "TEagle_other_annotations",
                      "These rows are annotations of RNA genes, artefacts or other non-TE features, reported "
                      "separately; they are not transposable-element family calls.")
        frow = QHBoxLayout(); frow.addStretch(1); frow.addWidget(_export_table_btn(t, "TEagle_family", self))
        cl.addLayout(frow)
        self._set_body(self.wslBody, cont)

    def _open_installer(self):
        """Open the dedicated component-wise installer dialog (per-package status, repair, integrity).
        Re-probe WSL status when it closes so the panel reflects any newly-installed backend."""
        from install_dialog import InstallDialog
        dlg = InstallDialog(self)
        # Destroy it on close. A QDialog with a parent is owned by Qt, and accept() only hides it, so
        # every open used to leave behind a live dialog with its own Engine and QThreadPool. The header
        # BACKEND button makes reopening cheap enough that a user checking on a long download would
        # accumulate them for the life of the session.
        dlg.setAttribute(Qt.WA_DeleteOnClose)
        dlg.finished.connect(lambda _=0: self._init_wsl())
        dlg.exec()
