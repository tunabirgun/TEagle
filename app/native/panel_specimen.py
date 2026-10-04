"""Panel 01 specimen input, request routing and the analysis run flow."""
from __future__ import annotations

import gzip
import html
import json
import os
import sys
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QFileDialog, QFrame, QGridLayout, QHBoxLayout, QLabel, QLineEdit, QMessageBox, QPlainTextEdit, QPushButton,
    QScrollArea, QSizePolicy, QTextEdit, QVBoxLayout, QWidget)
import widgets
from widgets import BusyBar
from sample import make_sample
import theme as theme_mod
from teagle_core import sequtil
from engine import LARGE_INPUT_BP
from ui_common import (
    PANEL, _Combo, _MetaLabel, _clear_layout, _disable_result_exports, _empty, _external, _kb_links, _note,
    invalid_input_message)


class SpecimenMixin:
    # ---------- rail ----------
    def _build_rail(self):
        rail = QFrame(); rail.setObjectName("rail")
        rail.setMinimumWidth(round(300 * theme_mod.UI_SCALE)); rail.setMaximumWidth(round(430 * theme_mod.UI_SCALE))
        lay = QVBoxLayout(rail); lay.setContentsMargins(*(theme_mod.sp(12),) * 4); lay.setSpacing(theme_mod.sp(8))
        lay.addWidget(self._sec(PANEL["specimen"], "Specimen"))
        accrow = QHBoxLayout()
        self.acc = QLineEdit(); self.acc.setPlaceholderText("accession — e.g. M11240, X04456")
        self.acc.setAccessibleName("Accession number")
        self.acc.returnPressed.connect(self._fetch)
        accrow.addWidget(self.acc)
        self.fetchBtn = QPushButton("Fetch"); self.fetchBtn.setProperty("sm", True); self.fetchBtn.clicked.connect(self._fetch)
        accrow.addWidget(self.fetchBtn)
        self.fetchCancelBtn = QPushButton("Cancel"); self.fetchCancelBtn.setProperty("sm", True)
        self.fetchCancelBtn.setToolTip("Stop waiting for the fetch. The request cannot be recalled, but whatever it "
                                       "returns is ignored and your current sequence stays as it is.")
        self.fetchCancelBtn.setVisible(False); self.fetchCancelBtn.clicked.connect(self._cancel_fetch)
        accrow.addWidget(self.fetchCancelBtn)
        lay.addLayout(accrow)
        self.accMeta = _MetaLabel(""); self.accMeta.setObjectName("cardmeta"); self.accMeta.setWordWrap(True)
        self.accMeta.setTextFormat(Qt.RichText); _external(self.accMeta); _kb_links(self.accMeta)
        lay.addWidget(self.accMeta)

        # coordinate fetch (collapsed) — organism + chr:start-end, like the UCSC browser position box
        self.coordToggle = QPushButton("▸ Fetch by coordinate"); self.coordToggle.setProperty("link", True)
        self.coordToggle.setAccessibleName("Fetch by coordinate")   # the ▸/▾ glyph is not speakable
        self.coordToggle.setAccessibleDescription("Disclosure toggle — shows or hides the coordinate-fetch fields.")
        self.coordToggle.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Fixed)
        self.coordToggle.clicked.connect(self._toggle_coord); lay.addWidget(self.coordToggle)
        # 2px insets stay literal: sp(2) snaps to 0 on the (0,6,10,16) ladder, which would delete the gap
        # rather than scale it. Only the spacing, which is on-ladder, goes through sp().
        self.coordBox = QWidget(); cb = QVBoxLayout(self.coordBox); cb.setContentsMargins(0, 2, 0, 2); cb.setSpacing(theme_mod.sp(5))
        # the rail grants this box less width than its combos ask for; unclamped, the rail's layout would
        # evaluate the box's heightForWidth at that unreachable minimum width and hand back one line too
        # few, clipping the wrapped hint below. Let it be as narrow as the rail actually makes it.
        self.coordBox.setMinimumWidth(1)
        orow = QHBoxLayout()
        self.asmSel = _Combo()
        self._rebuild_coord_asm_dropdown()                    # curated + any user-added assemblies
        orow.addWidget(self.asmSel, 1)
        self.coordStrand = _Combo(); self.coordStrand.addItems(["+ strand", "− strand"])
        self.coordStrand.setMaximumWidth(round(104 * theme_mod.UI_SCALE)); orow.addWidget(self.coordStrand)
        cb.addLayout(orow)
        self.coordCustom = QLineEdit(); self.coordCustom.setPlaceholderText("organism name or assembly accession (e.g. GCF_000001405.40)")
        self.coordCustom.setVisible(False); cb.addWidget(self.coordCustom)
        self.asmSel.currentIndexChanged.connect(lambda _=0: self.coordCustom.setVisible(self.asmSel.currentData() == "__custom__"))
        self.coord = QPlainTextEdit(); self.coord.setMaximumHeight(round(66 * theme_mod.UI_SCALE))
        self.coord.setPlaceholderText("chr13:33,016,423-33,066,143   (one region per line for multi-region)")
        cb.addWidget(self.coord)
        crow = QHBoxLayout()
        self.coordFetchBtn = QPushButton("Fetch region(s)"); self.coordFetchBtn.setProperty("sm", True)
        self.coordFetchBtn.clicked.connect(self._fetch_coord)
        crow.addWidget(self.coordFetchBtn); crow.addStretch(1); cb.addLayout(crow)
        self.coordMeta = _MetaLabel(""); self.coordMeta.setObjectName("cardmeta"); self.coordMeta.setWordWrap(True)
        self.coordMeta.setTextFormat(Qt.RichText); _external(self.coordMeta)
        _kb_links(self.coordMeta); cb.addWidget(self.coordMeta)
        cnote = QLabel("UCSC-style, 1-based (same numbers as the browser). Multi-region: every region is fetched, "
                       "sealed and analysed; choose the record to view in the summary table. Genome-anchored "
                       "export is available for the first region only.")
        cnote.setObjectName("orient"); cnote.setWordWrap(True); cb.addWidget(cnote)
        self.coordBox.setVisible(False); lay.addWidget(self.coordBox)

        ub = QPushButton("Upload FASTA (.fa / .fasta / .gz)"); ub.setProperty("sm", True)
        ub.clicked.connect(self._upload); lay.addWidget(ub)
        self.seq = QTextEdit()
        # A QTextEdit placeholder is painted on ONE line and clipped at the viewport — it never wraps. The
        # old string ran to "…Real IUPAC validation runs on analyze." and lost its last four words at 1.0x
        # and its last seven at 1.5x, so the only part a user ever read was the half that says nothing.
        # The clause is NOT relegated to a tooltip: a tooltip has no keyboard or touch path, so a scientific
        # statement parked there is unreachable for some users. It keeps a visible home in the IUPAC readout
        # below, which reads "—" before a run and reports the verdict (and any ambiguous-base content) after.
        self.seq.setPlaceholderText("…or paste DNA (FASTA or raw)")
        self.seq.setMinimumHeight(round(120 * theme_mod.UI_SCALE)); self.seq.textChanged.connect(self._seq_changed)
        self.seq.setTabChangesFocus(True)                     # Tab leaves the box; it must not type a character and drop the analysis
        self.seq.setAccessibleName("Specimen sequence")
        lay.addWidget(self.seq)
        row = QHBoxLayout()
        ls = QPushButton("Load example element ▾"); ls.setProperty("link", True); ls.clicked.connect(self._pick_example)
        self.exampleBtn = ls
        ls.setToolTip("Five published elements — copia, gypsy, LINE-1, Tc1, Ac — plus a synthetic negative control")
        row.addWidget(ls); row.addStretch(1)
        self.charCount = QLabel("0 nt"); self.charCount.setObjectName("kdim"); row.addWidget(self.charCount)
        lay.addLayout(row)
        self.runBtn = QPushButton("Run analysis"); self.runBtn.setProperty("primary", True)
        self.runBtn.clicked.connect(self._run_analysis); lay.addWidget(self.runBtn)
        # The loudest control on the first screen used to be armed with nothing to act on: clicking it before
        # loading anything only ever produced "paste, upload, or fetch a sequence first". A primary button
        # whose single guaranteed outcome is an error is pointing at the wrong step. It now waits for a
        # sequence, the same disabled-until-ready contract designBtn/annotateBtn/spliceBtn already use, which
        # leaves Fetch / Upload / Load example as the live affordances on an empty panel.
        self.runBtn.setEnabled(False)
        self.runBtn.setToolTip("Paste, upload, or fetch a sequence first — or use “Load example element”.")

        # readout gauges (2×2)
        lay.addSpacing(6)
        self.mScope = QLabel(); self.mScope.setObjectName("cardmeta"); self.mScope.setWordWrap(True)
        self.mScope.setVisible(False)                          # multi-record runs only: which record the tiles describe
        lay.addWidget(self.mScope)
        grid = QGridLayout(); grid.setHorizontalSpacing(6); grid.setVerticalSpacing(6)
        grid.setColumnStretch(0, 1); grid.setColumnStretch(1, 1)
        self.mLen = self._readout(grid, 0, "Length"); self.mGC = self._readout(grid, 1, "GC")
        self.mN = self._readout(grid, 2, "N content"); self.mValid = self._readout(grid, 3, "IUPAC")
        lay.addLayout(grid)
        lay.addSpacing(4)
        self.rRecords = _MetaLabel("—"); self.rStruct = _MetaLabel("—"); self.rOrf = _MetaLabel("—")
        for val in (self.rRecords, self.rStruct, self.rOrf):
            val.setObjectName("cardmeta"); val.setTextFormat(Qt.RichText)
        _external(self.rRecords)              # RECORDS -> the real fetched source accession (external DB link)
        self.rStruct.setOpenExternalLinks(False); self.rOrf.setOpenExternalLinks(False)   # in-app scroll, never an external link
        for _l in (self.rRecords, self.rStruct, self.rOrf):   # Tab-reachable (the in-app ones activate _scroll_to below)
            _kb_links(_l)
        self.rStruct.linkActivated.connect(lambda _: self._scroll_to(self.card_struct))
        self.rOrf.linkActivated.connect(lambda _: self._scroll_to(self.card_struct))
        for lbl, val, cap in (("RECORDS", self.rRecords, None),
                              ("STRUCTURAL EVIDENCE", self.rStruct,
                               "detected de novo (terminal-repeat heuristics) — not database-retrieved"),
                              ("ORFS (≥40 aa)", self.rOrf, None)):
            r = QHBoxLayout(); k = QLabel(lbl); k.setObjectName("kdim"); r.addWidget(k); r.addStretch(1); r.addWidget(val)
            lay.addLayout(r)
            if cap:
                c = QLabel(cap); c.setObjectName("gvpos"); c.setWordWrap(True); lay.addWidget(c)
                c.setVisible(False); self.rStructCap = c          # a method caption describes a result, so it appears with the first one

        lay.addSpacing(6)
        self.envHdr = QPushButton("▸ Environment"); self.envHdr.setObjectName("cardhdr")
        self.envHdr.setAccessibleName("Environment")
        self.envHdr.setAccessibleDescription("Disclosure toggle — shows or hides the backend environment status.")
        self.envHdr.clicked.connect(self._toggle_env); lay.addWidget(self.envHdr)
        self.envBox = QLabel("checking…"); self.envBox.setObjectName("cardmeta"); self.envBox.setWordWrap(True)
        self.envBox.setVisible(False); self.envBox.setTextFormat(Qt.RichText)
        lay.addWidget(self.envBox)
        lay.addStretch(1)
        note = QLabel("Order and superfamily (for example LTR/Copia, LTR/Gypsy, LINE, DNA/hAT) are called from "
                      "protein-domain architecture (HMMER with Pfam profiles) plus structural hallmarks, and are "
                      "withheld when the evidence is insufficient. Dfam / RepeatMasker family naming "
                      "runs in the managed WSL backend (panel 03).")
        note.setObjectName("orient"); note.setWordWrap(True); lay.addWidget(note)

        wrap = QScrollArea(); wrap.setWidgetResizable(True); wrap.setWidget(rail)
        wrap.setMinimumWidth(round(320 * theme_mod.UI_SCALE)); wrap.setMaximumWidth(round(440 * theme_mod.UI_SCALE))
        return wrap

    def _sec(self, n, title):
        w = QWidget(); r = QHBoxLayout(w); r.setContentsMargins(0, theme_mod.sp(4), 0, 2); r.setSpacing(theme_mod.sp(9))
        num = QLabel(n); num.setObjectName("secn"); num.setAlignment(Qt.AlignCenter); r.addWidget(num)
        h = QLabel(title); h.setObjectName("sech"); r.addWidget(h); r.addStretch(1)
        return w

    def _readout(self, grid, idx, label):
        cell = QFrame(); cell.setObjectName("cell")
        cl = QVBoxLayout(cell); cl.setContentsMargins(theme_mod.sp(10), theme_mod.sp(8), theme_mod.sp(10), theme_mod.sp(9)); cl.setSpacing(2)
        k = QLabel(widgets.caps(label)); k.setObjectName("kdim")
        v = QLabel("—"); v.setObjectName("value")
        cl.addWidget(k); cl.addWidget(v)
        r, c = divmod(idx, 2)
        grid.addWidget(cell, r, c)
        return v

    # ---------- input handling ----------
    def _set_seq(self, text):
        """Programmatic specimen load (fetch / upload / sample). Wrapped in the _loading guard so the
        textChanged handler keeps state['seq'] in sync but does NOT treat the load as a user edit that
        would wipe the source the loader is about to set."""
        self._invalidate_analysis("The specimen changed — run analysis to create current results and provenance.",
                                  cancel_request=True)
        self._drop_specimen_species()                         # a fetch or example sets its own organism right after
        self.state["features"] = None                         # drop a prior accession's gene model (fetch re-sets it after)
        self.state["coord_region"] = None                     # only a successful coordinate fetch may re-bind this
        self._splice_tx_origin = "external"                   # a new specimen invalidates any prior record-transcript flag
        self._loading = True
        try:
            self.seq.setPlainText(text)
        finally:
            self._loading = False

    def _seq_changed(self):
        txt = self.seq.toPlainText()
        body = "".join(l for l in txt.splitlines() if not l.startswith(">"))
        sig = "".join(txt.split())
        edited = sig != self._box_sig                          # a Tab, blank line or re-wrap changes no base and no header
        self._box_sig = sig
        if not self._loading and edited and (self.state.get("analysis") is not None or
                                             self.state.get("prov_manifest") is not None or
                                             self._active_analyze_gen is not None or self._reqs):
            self._invalidate_analysis("The sequence changed — the previous result is stale. Run analysis again.",
                                      cancel_request=True)
        self._set_char_count(txt)
        # Arm Run analysis exactly when there is something to analyse — never while a run is already in
        # flight, or a keystroke mid-analysis would re-enable the button _run_analysis just disabled.
        self.runBtn.setEnabled(bool(body.strip()) and not self._analyze_inflight)
        self.runBtn.setToolTip("Analyse the sequence in the box: structure, protein domains, classification."
                               if self.runBtn.isEnabled() else
                               "Paste, upload, or fetch a sequence first — or use “Load example element”.")
        self.state["seq"] = txt.strip()                       # specimen tracks the box: splice/annotate read this
        if not self._loading and edited:
            self._drop_specimen_species()
            self.state["coord_region"] = None                 # typed/pasted content has no verified genome-space binding
            if self.state.get("source") is not None:
                self.state["source"] = None                   # a genuine edit no longer matches the fetched identity
                self.accMeta.setText(""); self.coordMeta.setText("")
                self.state["features"] = None
            self._update_splice_ref()                         # refresh the 'Genomic reference' label to the edited specimen

    def _invalidate_analysis(self, message, cancel_request=False):
        """Remove result/export state that no longer describes the specimen or active request."""
        if cancel_request and self._active_analyze_gen is not None:
            self._analysis_requests.pop(self._active_analyze_gen, None)
            self._active_analyze_gen = None
            self._analyze_inflight = False
            widgets.set_caps(self.runBtn, "Run analysis")               # the cancelled run's busy label must not outlive it
        # Re-running analysis on an unchanged specimen leaves a running family or splice job alone, and a result
        # that returns meanwhile is held until the record is shown again; only a change of specimen retires it.
        if self._active_analyze_gen is None:
            self._release_held()                              # no re-analysis is left to show the record the held result needs
        live = set() if cancel_request else {k for k in ("annotate", "splice") if k in self._reqs}
        gone = self._drop_requests(specimen=cancel_request)
        for key in ("analysis", "records", "last_rec", "record_index", "selected_record", "record_headers",
                    "analyzed_seq", "analyzed_clean", "prov_manifest", "prov_manifest_all",
                    "lastPrimers", "candidates", "pcrPairs", "lastPcr", "lastGenomeScan",
                    "family", "familyResult", "family_seq", "splice", "splice_seq"):
            if (key in ("family_seq", "familyResult") and "annotate" in live) or (key == "splice_seq" and "splice" in live):
                continue                                      # the running job's hit coordinates still index this sequence
            self.state.pop(key, None)
        self.state["analysis_state"] = "stale"
        if hasattr(self, "mLen"):
            self._reset_rail_readouts()
        for dlg_name in ("_flank_dlg", "_selfsim_dlg"):       # tools built on the old record must not outlive it
            dlg = getattr(self, dlg_name, None)
            try:
                if dlg is not None:
                    dlg.close()                                # the self-similarity window is reused, so it is hidden, not freed
            except RuntimeError:
                pass
            if dlg_name == "_flank_dlg":
                self._flank_dlg = None
        # A retained Python wrapper remains callable after setParent(None). Every result container is
        # sanitized now and again by _clear_layout before detachment, so neither clicks nor table context
        # menus from a previous analysis generation can export after any invalidation path.
        roots = [getattr(self, "card_struct", None), getattr(self, "card_prov", None)]
        for name in ("primBody", "pcrBody", "wslBody", "spliceBody", "genomeBody"):
            layout = getattr(self, name, None)
            roots.append(layout.parentWidget() if layout is not None else None)
        for root in roots:
            _disable_result_exports(root)
        self._recTimer.stop(); self._rec_pending = None
        if self._recDetail is not None:
            _clear_layout(self._recDetail.layout())           # each result widget is detached (and its exports disabled) on its own
        self._recTable = self._recHead = self._recDetail = None
        if hasattr(self, "designBtn"):
            self.designBtn.setEnabled(False)
            self._sync_record_binding()
        if hasattr(self, "card_struct"):
            self.card_struct.clear_body()
            self.card_struct.bodylay.addWidget(_empty(message))
        if hasattr(self, "primBody"):
            self._set_body(self.primBody, _empty("Run analysis before designing primers."))
        if hasattr(self, "pcrBody"):
            self._set_body(self.pcrBody, _empty("Load a current primer pair, then run pair-aware amplicon search."))
        if hasattr(self, "pcrQueueBox"):
            self._render_pcr_queue()
        if hasattr(self, "wslBody") and "annotate" not in live:
            self._set_body(self.wslBody, _empty("Run family annotation on the current specimen."))
        if hasattr(self, "spliceBody") and "splice" not in live:
            self._set_body(self.spliceBody, _empty("Align a transcript to the current specimen to resolve splicing."))
        if hasattr(self, "genomeBody"):
            self._reset_genome_scan("The specimen changed — design a current primer pair before scanning a genome.")
        if hasattr(self, "card_prov"):
            self.card_prov.clear_body()
            self.card_prov.bodylay.addWidget(_empty(message))
        if gone:
            self._banner("A request that was still running when the results it belonged to were cleared was "
                         "discarded, and its result will not be shown: " + ", ".join(gone) + ".", "info")
        self._uppercase_buttons()

    # ---------- request generations ----------
    _TAGGED = ("annotate", "splice", "primers", "fetch")

    def _issue(self, kind, sequence=None, bound=True):
        """Key for a new engine request. The key carries a counter, and the request records the specimen
        generation and selected record it was sent under plus a digest of the sequence it carries, so a result
        that returns after the specimen or record changed is recognised and dropped instead of rendering as
        current. `bound=False` is for text that does not come from the specimen (a pasted sequence)."""
        self._req_n += 1
        key = f"{kind}#{self._req_n}"
        digest = sequtil.sha256(self._norm_seq(sequence)) if sequence else None
        self._reqs[kind] = {"key": key, "gen": self._spec_gen, "spec": self._specimen_gen, "digest": digest,
                            "bound": bound, "rec": self.state.get("record_index")}
        return key

    def _retire(self, key):
        """(operation, request) for a tagged key, or (None, None) when the request is no longer current: the
        specimen or record changed, it was cancelled, or a newer one replaced it. Untagged keys (direct calls,
        other operations) pass through unchanged."""
        kind, sep, _ = key.partition("#")
        if not sep or kind not in self._TAGGED:
            return key, None
        req = self._reqs.get(kind)
        if req is None or req["key"] != key:
            return None, None
        del self._reqs[kind]
        if kind == "fetch":
            return kind, req
        current = req["gen"] == self._spec_gen if kind == "primers" else req["spec"] == self._specimen_gen
        if not current or (req["bound"] and req["rec"] is not None and req["rec"] != self.state.get("record_index")):
            self._release(kind)                               # the request is gone from _reqs, so nothing else will free its control
            return None, None
        return kind, req

    _IDLE = {"annotate": ("annotateBtn", "Run family annotation", "wslBody"),
             "splice": ("spliceBtn", "Detect exons / introns", "spliceBody")}
    _DISCARDED = "This result was discarded because the analysis or record it belonged to was replaced while it ran. Run it again."

    def _release(self, kind):
        """Free the busy control of a primers, annotate or splice request that was retired: label and enabled state
        back to idle, and a still-spinning body replaced by plain text."""
        if kind == "primers":
            self._design_inflight = False; self._pending_domain = None
            if hasattr(self, "designBtn"):
                widgets.set_caps(self.designBtn, "Design primers")
                self.designBtn.setEnabled(self._design_ready())
                self._unbusy(self.primBody)
            return
        btn, label, body = self._IDLE[kind]
        if hasattr(self, btn):
            widgets.set_caps(getattr(self, btn), label)
            getattr(self, btn).setEnabled(self._wsl_ready[kind])
            self._unbusy(getattr(self, body))

    def _unbusy(self, layout):
        if layout.count() and isinstance(layout.itemAt(0).widget(), BusyBar):
            self._set_body(layout, _empty(self._DISCARDED))

    def _scan_is_stale(self):
        return getattr(self, "_scan_gen", self._specimen_gen) != self._specimen_gen

    def _drop_requests(self, specimen=False):
        """Retire every in-flight request whose result would describe results that are being cleared, and
        release its busy control. Primer design and in-silico PCR belong to the analysis, so any invalidation
        retires them; family annotation, splice, a fetch and a whole-genome scan belong to the specimen and
        are retired only when it changes. Returns plain-language names of what was discarded. A running scan
        cannot be recalled: it keeps its busy flag until the job returns, and then its result is dropped
        (_scan_is_stale)."""
        self._spec_gen += 1
        if specimen:
            self._specimen_gen += 1
        gone = []
        if self._reqs.pop("primers", None):
            gone.append("primer design")
            self._release("primers")
        for kind, name in (("annotate", "family annotation"), ("splice", "splice detection")):
            if specimen and self._reqs.pop(kind, None):
                gone.append(name)
                self._release(kind)
        run = getattr(self, "_pcr_run", None)
        if run and any(r is None for r in run["results"]):
            gone.append("in-silico PCR")
            if hasattr(self, "runPcrBtn"):
                widgets.set_caps(self.runPcrBtn, "Run loaded pairs")
        self._pcr_gen += 1
        self._pcr_run = None
        if specimen and self._genome_inflight:
            gone.append("whole-genome scan")
        if specimen and self._reqs.pop("fetch", None):
            gone.append("accession fetch")
            self._restore_fetch_ui()
        return gone

    def _reset_rail_readouts(self):
        """The rail tiles describe the analysed record; with no current analysis they must not keep showing
        the previous one's length, GC or counts."""
        for w in (self.mLen, self.mGC, self.mN, self.mValid):
            w.setText("—"); w.setToolTip("")
        self.mValid.setProperty("state", ""); self._repolish(self.mValid)
        for w in (self.rRecords, self.rStruct, self.rOrf):
            w.setText("—"); w.setToolTip("")
        self.rStructCap.setVisible(False)

    def _set_char_count(self, txt):
        """Letters in the box, counted as 'characters' until a record has validated as nucleotide sequence."""
        body = "".join(l for l in txt.splitlines() if not l.startswith(">"))
        n = sum(c.isalpha() for c in body)
        other = sum(not (c.isalpha() or c.isspace()) for c in body)
        valid = self.state.get("analysis_state") == "current" and (self.state.get("last_rec") or {}).get("valid")
        self.charCount.setText(f"{n:,} {'nt' if valid else 'characters'}" + (f" · {other:,} non-letters" if other else ""))
        self.charCount.setToolTip("Letters only: digits, spaces and line breaks are not counted. The count reads "
                                  "“nt” once a record has been validated as nucleotide sequence.")

    def _apply_example(self, fasta):
        self._set_seq(fasta)
        self.state["source"] = None
        self.accMeta.setText(""); self.coordMeta.setText("")
        self._update_splice_ref()                             # _set_seq runs under the _loading guard, so refresh the label here

    def _load_sample(self):
        """Load the generated demo construct directly.

        Kept as a NON-BLOCKING, no-argument action because it is the programmatic entry point: the test
        suite calls it, and turning it into a modal menu made every one of those calls block forever —
        tests/test_native.py went from seconds to a hard hang. The user-facing picker is _pick_example."""
        self._apply_example(make_sample())

    def _pick_example(self):
        """Offer the real published specimens first; the generated construct stays, labelled for what it
        is. Its ORF is random sequence, so it cannot match a Pfam profile — landing a first-time user on
        a structural-only null result and making the domain panel look broken when it is working."""
        from PySide6.QtWidgets import QMenu
        from teagle_core import examples
        menu = QMenu(self)
        menu.setToolTipsVisible(True)                      # QMenu hides action tooltips by default (Qt 5.1+)
        for acc, label, organism, expected in examples.available():
            act = menu.addAction(f"{label} — {organism}  ({acc})")
            act.setToolTip(f"Published as: {expected}")
            act.triggered.connect(lambda _=False, a=acc: self._load_example(a))
        if not menu.isEmpty():
            menu.addSeparator()
        neg = menu.addAction("Synthetic construct (negative control)")
        neg.setToolTip("Seeded random LTR-like sequence. Structure is detectable; the ORF is random, so "
                       "no protein domain should be found — that null result is the expected outcome.")
        neg.triggered.connect(lambda _=False: self._apply_example(make_sample()))
        menu.exec(self.exampleBtn.mapToGlobal(self.exampleBtn.rect().bottomLeft()))   # at the button, not wherever the mouse idles

    def _load_example(self, accession):
        from teagle_core import examples
        fasta = examples.load(accession)
        if not fasta:
            return self._banner("This build does not include that example sequence.", "warn")
        self._apply_example(fasta)
        self._set_species(next((o for a, _l, o, _e in examples.CATALOG if a == accession), ""))   # the example's own organism
        return True

    def _upload(self):
        path, _ = QFileDialog.getOpenFileName(self, "Open FASTA", "",
                                              "FASTA (*.fa *.fasta *.fna *.txt *.gz);;All files (*)")
        if not path:
            return
        if os.path.isdir(path):
            # Windows raises PermissionError (errno 13) for open() on a directory, not IsADirectoryError, so
            # the exception type cannot tell these apart here — ask the filesystem instead of the errno.
            return self._banner(f"“{os.path.basename(path.rstrip(os.sep))}” is a folder, not a sequence "
                                "file. Choose a .fa, .fasta, .fna, .txt or .gz file.")
        try:
            if path.lower().endswith(".gz"):
                with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
                    data = f.read()
            else:
                with open(path, "r", encoding="utf-8", errors="replace") as f:
                    data = f.read()
        except (gzip.BadGzipFile, EOFError):
            # BEFORE the OSError arm: BadGzipFile subclasses OSError, so it would otherwise be swallowed by
            # the generic message. EOFError (truncated stream) is not an OSError at all.
            return self._banner(f"“{os.path.basename(path)}” looks like a .gz file but could not be "
                                "decompressed — it may be incomplete or not really gzip-compressed.")
        except OSError as e:
            # An error says what went wrong and what to do next, and never puts a raw
            # traceback or errno in front of the user. str(OSError) is "[Errno 13] Permission denied:
            # 'C:\\Users\\...\\x.fa'" — an errno and a full filesystem path, which is neither.
            name = os.path.basename(path)
            if isinstance(e, PermissionError):
                msg = (f"“{name}” could not be opened — it may be open in another program, or you may not "
                       "have permission to read it. Close it elsewhere and try again.")
            elif isinstance(e, FileNotFoundError):
                msg = f"“{name}” is no longer there — it may have been moved or renamed since you picked it."
            else:
                msg = (f"“{name}” could not be read. If it is on a network or removable drive, copy it to "
                       "this computer and try again.")
            return self._banner(msg)
        self._set_seq(data)
        self.state["source"] = None
        self.accMeta.setText(f"loaded {os.path.basename(path)}"); self.coordMeta.setText("")
        self._update_splice_ref()                             # _set_seq runs under the _loading guard, so refresh the label here

    def _set_fetch_enabled(self, on):
        self.fetchBtn.setEnabled(on); self.coordFetchBtn.setEnabled(on)   # both disabled in-flight: no overlapping fetch race
        self.fetchCancelBtn.setVisible(not on)

    def _start_fetch(self, op, body, label):
        """Submit a fetch under its own key. Until it returns, loading, pasting or editing a specimen cancels
        it (see _drop_requests), so a slow download can never replace what the user put there afterwards."""
        self._fetch_prior = (self.accMeta.text(), self.coordMeta.text())
        label.setText("fetching…")
        self._set_fetch_enabled(False)
        self.engine.submit(op, body, key=self._issue("fetch"))

    def _restore_fetch_ui(self):
        """Fetch ended without loading anything (failed, cancelled, or superseded): re-enable the controls and
        put back the label of the specimen that is still loaded."""
        self._set_fetch_enabled(True)
        prior = getattr(self, "_fetch_prior", ("", ""))
        for lbl, old in zip((self.accMeta, self.coordMeta), prior):
            if lbl.text() == "fetching…":
                lbl.setText(old)

    def _cancel_fetch(self):
        if self._reqs.pop("fetch", None):
            self._restore_fetch_ui()
            self._banner("Fetch cancelled — your current sequence is unchanged. If the download still finishes, "
                         "its result is ignored.", "info")

    def _fetch(self):
        if "fetch" in self._reqs:
            return self._banner("A fetch is already running — wait for it, or press Cancel.", "info")
        acc = self.acc.text().strip()
        if not acc:
            return self._banner("Enter an accession first (e.g. M11240).", "info")
        self._start_fetch("fetch", {"accession": acc}, self.accMeta)

    def _confirm_large_input(self, length, index, n_records):
        """Ask before analysing a record longer than the engine's large-input threshold: the element detectors assume one
        locus, so a contig still gets ONE classification. Returns 'run', 'panel03' or 'cancel'. The single seam tests and
        screenshot harnesses replace; nothing here changes the analysis or the thresholds."""
        box = QMessageBox(self); box.setIcon(QMessageBox.Information); box.setWindowTitle("Large sequence")
        box.setText(f"Record {index} is {length:,} bp." if n_records > 1 else f"{length:,} bp.")
        box.setInformativeText(
            "TEagle's element analysis is tuned for one element at a time: the terminal-repeat and ORF detectors assume a "
            "single locus. On a contig, BAC or chromosome segment it still returns ONE classification for the whole "
            "sequence, which may span many elements. To annotate each element in this sequence, use panel "
            f"{PANEL['family']} (RepeatMasker against Dfam; needs the Linux backend).")
        run = box.addButton("Analyse anyway", QMessageBox.AcceptRole)
        p3 = box.addButton(f"Go to panel {PANEL['family']}", QMessageBox.ActionRole)
        cancel = box.addButton("Cancel", QMessageBox.RejectRole)
        box.setDefaultButton(run); box.setEscapeButton(cancel)
        widgets.uppercase_buttons(box)
        box.exec()
        return "run" if box.clickedButton() is run else "panel03" if box.clickedButton() is p3 else "cancel"

    def _run_analysis(self):
        seq = self.seq.toPlainText().strip()
        if not seq:
            return self._banner("Paste, upload or fetch a sequence first.", "info")
        parsed = sequtil.parse_fasta(seq)
        big = next(((i, len(s)) for i, (_h, s) in enumerate(parsed, 1) if len(s) > LARGE_INPUT_BP), None)
        if big:                                               # before anything is invalidated, so Cancel keeps the current result
            choice = self._confirm_large_input(big[1], big[0], len(parsed))
            if choice == "cancel":
                return
            if choice == "panel03":
                self.wslSource.setCurrentIndex(0)             # 'Loaded specimen (panel 01)'
                return self._scroll_to(self.card_wsl)
        self._invalidate_analysis("Analysis is running — provenance will appear when this request completes.")
        self._analyze_gen += 1
        generation = self._analyze_gen
        self._active_analyze_gen = generation
        self._analysis_requests[generation] = {"sequence": seq, "source": dict(self.state["source"] or {})}
        self.state["seq"] = seq
        self.state["analyzed_seq"] = seq                  # snapshot the sequence the reported feature coords index
        self.state["analysis_state"] = "running"
        self._analyze_inflight = True
        self.runBtn.setEnabled(False); widgets.set_caps(self.runBtn, "… analysing")
        self.engine.submit("analyze", {"sequence": seq, "source": self.state["source"]},
                           key=f"analyze#{generation}")

    # ---------- results routing ----------
    def _hold_result(self, key, res):
        """A family or splice result for the unchanged specimen that returns while it is being re-analysed has no
        selected record to bind to yet. Keep it on its live request; _show_record replays it onto the record it was
        run on, or frees the control if that record is not the one shown."""
        kind = key.partition("#")[0]
        req = self._reqs.get(kind)
        if (kind in self._IDLE and req and req["key"] == key and req["bound"] and req["rec"] is not None
                and req["spec"] == self._specimen_gen and self.state.get("record_index") is None
                and self._active_analyze_gen is not None):
            req["held"] = res
            return True
        return False

    def _release_held(self):
        for kind in [k for k, r in self._reqs.items() if "held" in r]:
            del self._reqs[kind]
            self._release(kind)

    def _on_done(self, key, res):
        if self._hold_result(key, res):
            return
        key, req = self._retire(key)
        if key is None:
            return
        sent = (req or {}).get("digest")
        echoed = ((res.get("provenance") or {}).get("input") or {}).get("sha256") if isinstance(res, dict) else None
        if sent and isinstance(echoed, str) and echoed != sent:
            self._reset_buttons(key)                          # the result describes some other sequence than the one sent
            return self._banner("A result came back for a different sequence than the one that was sent, so it "
                                "was not shown. Run the step again.", "warn")
        if key == "health":
            self._core = {"ok": True, "primer3": res.get("primer3")}
            self._render_status_chip()
            if res.get("core"):
                self.ver.setText("v" + res["core"])
        elif key == "env":
            self._render_env(res)
        elif key == "fetch":
            self._on_fetch(res)
        elif key.startswith("analyze#"):
            try:
                generation = int(key.split("#", 1)[1])
            except ValueError:
                return
            self._on_analyze(res, generation)
        elif key == "wsl_status":
            self._on_wsl_status(res)
        elif key == "annotate":
            self._on_annotate(res)
        elif key == "splice":
            self._on_splice(res)
        elif key == "primers":
            self._on_primers(res)
        elif key.startswith("pcr#"):
            self._pcr_slot(key, res)
        elif key == "genome_pcr":
            self._on_genome_pcr(res)
        elif key == "genome_prepare":
            self._on_genome_prepare(res)
        elif key == "genome_prepare_log":
            self._on_genome_prepare_log(res)
        elif key == "genome_list":
            self._on_genome_list(res)
        elif key == "genome_remove":
            self._on_genome_remove(res)
        elif key == "add_custom_assembly":
            self._on_add_custom_assembly(res)
        elif key == "annotate_budget":
            self._on_annotate_budget(res)
        elif key == "genome_annotate":
            self._on_genome_annotate(res)
        elif key == "genome_annotate_log":
            self._on_genome_annotate_log(res)
        elif key == "genome_annotate_reset":
            self._on_genome_annotate_reset(res)

    def _reset_buttons(self, key):
        if key == "analyze":
            self._analyze_inflight = False
            # re-arm only if the box still holds something — an analysis that failed on an empty/cleared box
            # must not leave the primary button armed with nothing to act on
            self.runBtn.setEnabled(bool(self.seq.toPlainText().strip())); widgets.set_caps(self.runBtn, "Run analysis")
        elif key == "primers":
            self._design_inflight = False
            self._pending_domain = None                       # a failed routed design must not leave a stale re-anchor or busy body
            self._invalidate_primer_results()
            self.designBtn.setEnabled(self._design_ready()); widgets.set_caps(self.designBtn, "Design primers")
            self._set_body(self.primBody, _note("Primer design did not complete — the message explains why. The cause can "
                                                "be the sequence or region as well as the settings; fix it, then try again."))
        elif key == "annotate":
            self.annotateBtn.setEnabled(self._wsl_ready["annotate"]); widgets.set_caps(self.annotateBtn, "Run family annotation")
            # _reset_buttons is error-only (success re-enables in _on_annotate), so clear the stuck busy body here —
            # a BadRequest/fault never reaches _on_annotate, which would otherwise leave "Running RepeatMasker…" spinning
            self._set_body(self.wslBody, _empty("Run RepeatMasker against Dfam to name the TE family. Family naming is the Linux (WSL) backend."))
        elif key == "splice":
            self.spliceBtn.setEnabled(self._wsl_ready["splice"]); widgets.set_caps(self.spliceBtn, "Detect exons / introns")
            self._set_body(self.spliceBody, _empty("Align a transcript (mRNA / cDNA) to the loaded sequence to infer its exon–intron structure."))
        elif key == "fetch":
            self._restore_fetch_ui()                          # fetch failed — re-enable so the user can retry, keep the prior specimen label
        elif key in ("genome_list", "genome_remove"):
            if getattr(self, "_genome_mgr", None) is None:
                self._genome_mgr_open = False                 # a failed list must not leave the manager armed to pop up on the next refresh
        elif key == "env":
            self.envBox.setText("environment check unavailable")
        elif key == "wsl_status":
            self._wsl = {"error": "status could not be read"}
            self._render_status_chip(); self._tag_backend_cards()
            for lab in (self.wslStatus, self.spliceStatus):
                lab.setText("Linux (WSL) backend status could not be read — open Backend in the header to check it.")
        elif key in ("genome_annotate", "annotate_budget"):
            # a failed annotation (or budget probe) must clear the guard, stop the progress poll and
            # re-enable the manager's action — otherwise the feature is dead until the app restarts.
            self._annot_inflight = False
            t = getattr(self, "_annot_timer", None)
            if t is not None:
                try:
                    t.stop()
                except RuntimeError:
                    pass
            self._refresh_genome_manager()
        elif key == "add_custom_assembly":
            self._add_asm_inflight = False                    # bad query / fault — clear the guard and re-enable so the user can retry
            e = getattr(self, "_addAsmEdit", None)
            b = getattr(self, "_addAsmBtn", None)
            for w in (e, b):
                try:
                    if w is not None:
                        w.setEnabled(True)
                except RuntimeError:
                    pass

    def _on_user_error(self, key, msg):
        key, _ = self._retire(key)
        if key is None:
            return                                            # a superseded request: nothing on screen belongs to it
        if key.startswith("analyze#"):
            try:
                generation = int(key.split("#", 1)[1])
            except ValueError:
                return
            self._analysis_requests.pop(generation, None)
            if generation != self._active_analyze_gen:
                return
            self._active_analyze_gen = None
            self._invalidate_analysis("Analysis failed — no result or provenance from this request is current.")
            key = "analyze"
        if key.startswith("pcr#"):                            # a failed pair fills its slot so the batch still renders
            self._pcr_slot(key, {"error": msg, "amplicons": []})
        elif key == "genome_pcr":
            if self._scan_is_stale():
                return self._discard_stale_scan()
            self._genome_inflight = False; self._render_genome_status("Genome scan failed — " + msg)
            self._refresh_genome_manager()                    # scan settled (failed) — re-enable a manager opened mid-scan
            self._banner("Genome scan failed — " + msg, "warn"); return
        elif key == "genome_prepare":
            self._genome_prep_inflight = False; self._pending_scan = None
            self._render_genome_status("Genome download failed — " + msg)
            self._refresh_genome_manager()                    # download settled (failed) — re-enable a manager opened during it
            self._banner("Genome download failed — " + msg, "warn"); return
        elif key == "genome_prepare_log":
            return                                            # a background poll blip: never banner
        else:
            self._reset_buttons(key)
        self._banner(msg)

    @staticmethod
    def _log_trace(trace):
        """Best-effort diagnostic log. A windowed (frozen) build has no stderr, so recovery never depends on it."""
        try:
            if sys.stderr is not None:
                sys.stderr.write(str(trace) + "\n")
        except Exception:
            pass

    def _on_failed(self, key, msg, trace):
        try:
            self._recover_from_fault(key, msg)
        finally:
            self._log_trace(trace)                            # last, so the interface is already usable again

    def _discard_stale_scan(self):
        """A genome scan that was running when the specimen changed has returned: say so, show nothing from it."""
        self._genome_inflight = False; self._pending_scan = None
        self._set_body(self.genomeBody, _empty("The specimen changed while this scan was running, so its result "
                                               "was discarded. Design a pair for the current specimen (panel 04), "
                                               "then scan again."))
        self._refresh_genome_manager()

    def _recover_from_fault(self, key, msg):
        if key == "genome_prepare_log":
            return                                            # a background poll blip: never banner, never spam a traceback
        key, _ = self._retire(key)
        if key is None:
            return
        if key.startswith("analyze#"):
            try:
                generation = int(key.split("#", 1)[1])
            except ValueError:
                return
            self._analysis_requests.pop(generation, None)
            if generation != self._active_analyze_gen:
                return
            self._active_analyze_gen = None
            self._invalidate_analysis("Analysis failed — no result or provenance from this request is current.")
            key = "analyze"
        if key == "health":                                   # state the engine is down IN WORDS — the chip must not
            self._core = {"ok": False}                        # rely on the LED hue alone, nor keep saying "connecting…"
            self._render_status_chip()
        if key.startswith("pcr#"):
            self._pcr_slot(key, {"error": msg, "amplicons": []})
        elif key == "genome_pcr":                             # a genome scan fault is a handled outcome, not an "unexpected error"
            if self._scan_is_stale():
                return self._discard_stale_scan()
            self._genome_inflight = False; self._render_genome_status("Genome scan failed — " + msg)
            self._refresh_genome_manager()                    # scan settled (failed) — re-enable a manager opened mid-scan
            self._banner("Genome scan failed — " + msg, "warn"); return
        elif key == "genome_prepare":                         # likewise a download fault: one clear warning, not a red crash dialog
            self._genome_prep_inflight = False; self._pending_scan = None
            self._render_genome_status("Genome download failed — " + msg)
            self._refresh_genome_manager()                    # download settled (failed) — re-enable a manager opened during it
            self._banner("Genome download failed — " + msg, "warn"); return
        else:
            self._reset_buttons(key)
        # An error never puts a raw traceback in front of the user as the primary message.
        # `msg` is "{ExceptionType}: {str(e)}" straight from engine_worker, so "unexpected error: KeyError:
        # 'accession'" was the whole sentence a bench scientist got. Lead with what failed and what to do;
        # keep the technical string, clearly subordinate, because it is what they would quote in a bug report.
        # The full traceback goes to the log (_on_failed) and is not shown here.
        # No "nothing was changed" reassurance here: this branch is reached by keys that DO change state
        # (genome_annotate, genome_remove), so that would be a claim the handler cannot make.
        step = key.partition("#")[0]
        name = self._STEP_NAMES.get(step) or f"The “{step.replace('_', ' ')}” step"
        self._banner(f"{name} could not be completed. Try it again. If it keeps happening, restart TEagle and "
                     f"report the detail below.\n\nDetail: {msg}")

    # Plain names for the engine keys, so a failure says which step failed in the words the panel uses.
    _STEP_NAMES = {
        "analyze": "The analysis", "fetch": "The accession fetch", "coord_fetch": "The coordinate fetch",
        "primers": "Primer design", "pcr": "The in-silico PCR", "genome_pcr": "The whole-genome scan",
        "genome_prepare": "The genome download", "genome_list": "The genome list",
        "genome_remove": "Deleting the genome", "annotate": "Dfam/RepeatMasker family naming",
        "genome_annotate": "Whole-genome annotation", "splice": "Splice-structure alignment",
        "oligoqc": "Primer secondary-structure QC", "health": "The backend check",
        "wsl_status": "The Linux (WSL) backend status check", "env": "The environment check",
        "annotate_budget": "The whole-genome annotation size check", "genome_annotate_log": "The annotation progress check",
        "genome_annotate_reset": "Resetting the genome annotation", "add_custom_assembly": "Adding the assembly",
    }

    def _toggle_coord(self):
        vis = not self.coordBox.isVisible()
        self.coordBox.setVisible(vis)
        self.coordToggle.setText("▾ Fetch by coordinate" if vis else "▸ Fetch by coordinate")

    def _fetch_coord(self):
        self._clear_banner()
        org = self.asmSel.currentData()
        custom = self.coordCustom.text().strip() if org == "__custom__" else ""
        if org == "__custom__" and not custom:
            return self._banner("Enter an organism name or an assembly accession first (e.g. GCF_000001405.40).", "info")
        regions = self.coord.toPlainText().strip()
        if not regions:
            return self._banner("Enter at least one region first, e.g. chr13:33,016,423-33,066,143.", "info")
        strand = "-" if self.coordStrand.currentIndex() == 1 else "+"
        self._start_fetch("fetch_coords", {"regions": regions, "strand": strand,
                                           "organism": "" if org == "__custom__" else org,
                                           "customQuery": custom}, self.coordMeta)

    def _render_coord_fetch(self, res):
        self.accMeta.setText("")                          # only one specimen identity shows at a time
        regions = res.get("regions", [])
        cached = " · cached (local)" if res.get("fromCache") else ""
        ncbi = (self._src_html("NCBI", "https://www.ncbi.nlm.nih.gov/nuccore/" + regions[0]["chrAccession"])
                if regions else "")
        lines = [f"{r.get('chromLabel','')}:{r.get('start'):,}-{r.get('stop'):,} · {r.get('chrAccession','')} · "
                 f"{r.get('stop',0)-r.get('start',0)+1:,} bp" + ("  (−)" if r.get('strand') == 2 else "")
                 for r in regions]
        self.coordMeta.setText(f"{res.get('assemblyName','')} · {res.get('organism','')}{cached}{ncbi}<br>" + "<br>".join(lines))
        _kb_links(self.coordMeta)                             # text now carries a link -> re-arm the tab stop
        from teagle_core import provenance
        parsed = sequtil.parse_fasta(self.state.get("seq", ""))
        source_context = (provenance.validate_source_context(res.get("sourceContext"), parsed)
                          if res.get("sourceContext") is not None else
                          provenance.normalize_source_context(res.get("source"), parsed))
        self.state["source"] = dict(source_context)
        # Keep the fetched region so an annotation export can be anchored back into genome space. Only the
        # first region is anchored (stated in the panel); every region is still analysed.
        records = sequtil.parse_fasta(self.state.get("seq", ""))
        bindings = self.state["source"].get("records") or []
        if (regions and records and bindings and
                self.state["source"].get("inputType") == "fetched_coordinate"):
            _record_id, sequence = records[0]
            bound = dict(regions[0])
            bound["_record_id"] = bindings[0]["id"]
            bound["_sequence_sha256"] = sequtil.sha256(sequence)
            bound["_source_binding"] = json.dumps(self.state["source"], sort_keys=True, ensure_ascii=False,
                                                    allow_nan=False, default=str)
            self.state["coord_region"] = bound
        else:
            self.state["coord_region"] = None
        self.state["features"] = None

    def _on_fetch(self, res):
        self._set_fetch_enabled(True)                         # fetch settled — allow the next one
        if not res.get("ok"):
            self._restore_fetch_ui()                          # the specimen that is still loaded keeps its label and genome binding
            return self._banner(res.get("error", "The fetch failed."))
        self._clear_banner()
        seqtext = res.get("fasta") or res.get("sequence") or ""
        if seqtext:
            self._set_seq(seqtext)
        org = res.get("organism", "")
        if res.get("runType") == "coordinate":
            self._render_coord_fetch(res)
        else:
            length = res.get("length") or res.get("seq_length") or ""
            cached = " · cached (local)" if res.get("fromCache") else ""
            acc = res.get("accession", "")
            src_label = "ENA" if str(res.get("source", "")).startswith("ENA") else "NCBI"   # name the DB that served it
            src_url = res.get("sourceUrl") or ("https://www.ncbi.nlm.nih.gov/nuccore/" + acc)
            link = self._src_html(src_label, src_url) if acc else ""
            self.accMeta.setText(f"{acc} · {org} · {length} bp{cached}{link}<br>{res.get('title','')}")
            status = str(res.get("status") or "")
            if status and status != "live":                   # a withdrawn or replaced record is still loaded; say so beside it
                note = f"NCBI status: {status}" + (f" — replaced by {res['replacedBy']}" if res.get("replacedBy") else "")
                self.accMeta.setText(self.accMeta.text() + f"<br><b style='color:{theme_mod.FLAG[self.theme]['caution']}'>"
                                     f"{html.escape(note)}</b>")
            _kb_links(self.accMeta)                           # text now carries a link -> re-arm the tab stop
            self.coordMeta.setText("")                    # clear the other specimen identity
            from teagle_core import provenance
            parsed = sequtil.parse_fasta(seqtext)
            source_context = (provenance.validate_source_context(res.get("sourceContext"), parsed)
                              if res.get("sourceContext") is not None else
                              provenance.normalize_source_context(res, parsed))
            self.state["source"] = dict(source_context)
            self.state["features"] = res.get("features")
        # auto-fill species for WSL family annotation if present
        if hasattr(self, "wslSpecies") and org:
            self._set_species(org)
        self._update_splice_ref()

    def _on_analyze(self, res, generation=None):
        if generation is None:
            return
        request = None
        request = self._analysis_requests.pop(generation, None)
        if generation != self._active_analyze_gen or request is None:
            return
        self._active_analyze_gen = None
        self._analyze_inflight = False
        self.runBtn.setEnabled(bool(self.seq.toPlainText().strip())); widgets.set_caps(self.runBtn, "Run analysis")
        self._clear_banner()
        recs = res.get("records", [])
        if not recs:
            if res.get("warning"):
                self._banner(res["warning"], level="warn")
            self._invalidate_analysis("Analysis returned no records — no result or provenance is current.")
            return
        if res.get("status") == "no_valid_records":           # nothing to analyse: say what is wrong with the input, and where
            return self._reject_input(res, request["sequence"])
        invalid = [r for r in recs if not r.get("valid")]
        warning = res.get("warning") or ""
        if invalid and warning:                               # the engine says how many were skipped; this says which, and why
            warning += "\n\n" + invalid_input_message(recs, request["sequence"])
        if warning:
            self._banner(warning, level="warn")
        manifest = res.get("provenance")
        if not isinstance(manifest, dict) or manifest.get("runType") != "analysis":
            self._invalidate_analysis("Analysis returned no valid provenance — its result was not accepted.")
            return
        import engine
        try:
            engine.analysis_result_document(res)
        except ValueError:
            self._invalidate_analysis("Analysis did not match its provenance binding — its result was not accepted.")
            return
        self.state["record_headers"] = [h for h, _ in sequtil.parse_fasta(request["sequence"])]
        self.state["analysis"] = res
        self.state["records"] = recs
        snapshot = request["sequence"] if request else (self.state.get("analyzed_seq") or "".join(r.get("seq", "") for r in recs))
        self.state["analyzed_seq"] = snapshot
        self.state["analyzed_clean"] = self._clean_seq(snapshot)   # snapshot, not the live box — a keystroke mid-analysis must not defeat the stale-block guard
        self.state["analysis_state"] = "current"
        self._reset_genome_scan("New specimen analysed — the previous whole-genome scan described different "
                                "primers. Design a pair (panel 04), then scan here.")
        self._show_record(next((i for i, rec in enumerate(recs) if rec.get("valid")), 0))
        self._render_provenance(manifest)
        self._replay_held()                                   # after the analysis seal, so a held family or splice result keeps the current provenance
        self._scroll_to_top(self.card_struct)                 # land on the new call, not wherever the last result left the page

    def _reject_input(self, res, raw_text):
        """Every record was unusable: show what is wrong with it in the results card and leave no result,
        export or provenance standing. This is a problem with the input, never a seal failure."""
        text = invalid_input_message(res.get("records") or [], raw_text, res.get("warning"))
        self._invalidate_analysis("Nothing was analysed.")
        note = _note("Nothing was analysed — no sequence in the box could be read as DNA or RNA.\n\n" + text +
                     "\n\nFix the sequence in the box and run the analysis again.", "warn")
        note.setTextFormat(Qt.PlainText)
        note.setTextInteractionFlags(Qt.TextSelectableByMouse)
        self.card_struct.clear_body()
        self.card_struct.bodylay.addWidget(note)
        self.card_struct.expand()
        self.card_prov.clear_body()
        self.card_prov.bodylay.addWidget(_empty("No result, so no provenance was produced."))
        self._scroll_to_top(self.card_struct)

    def _replay_held(self):
        """Deliver each result held through a re-analysis to the record now shown, or free its control when the
        record is not the one it was run on."""
        rec = self.state.get("last_rec") or {}
        for kind, req in [(k, r) for k, r in self._reqs.items() if "held" in r]:
            if req["digest"] in (None, sequtil.sha256(self._norm_seq(rec.get("seq") or ""))):
                self._on_done(req["key"], req.pop("held"))
            else:
                del self._reqs[kind]
                self._release(kind)
