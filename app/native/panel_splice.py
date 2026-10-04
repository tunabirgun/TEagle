"""Splice structure panel."""
from __future__ import annotations

import html
import re
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QHBoxLayout, QLabel, QPushButton, QSizePolicy, QTextEdit, QVBoxLayout, QWidget
import figures
from figures import svg_genome
import widgets
from widgets import BusyBar
from coords import CoordGenomePanel
import theme as theme_mod
from teagle_core.fetch import complete_gene_model, cross_check_models
from ui_common import (
    CollapsibleCard, PANEL, _Combo, _empty, _export_fasta_btn, _export_table_btn, _external, _kb_links, _note,
    _record_public_identity, _splice_exons_fasta)


class SpliceMixin:
    # ---------- 03b Splice structure (transcript alignment) ----------
    def _build_splice_card(self):
        card = CollapsibleCard(PANEL["splice"], "Splice structure (transcript alignment)",
                               "optional · Linux backend: checking… · exon–intron structure from an aligned transcript")
        self.card_splice = card
        self.spliceStatus = QLabel("checking splice-alignment backend…"); self.spliceStatus.setObjectName("cardmeta")
        self.spliceStatus.setWordWrap(True); self.spliceStatus.setTextFormat(Qt.RichText)
        card.bodylay.addWidget(self.spliceStatus)
        # the genomic reference is always the specimen loaded in panel 01 (fetched/uploaded/pasted)
        self.spliceRef = QLabel("Genomic reference: none loaded yet — load a specimen in panel 01.")
        self.spliceRef.setObjectName("cardmeta"); self.spliceRef.setWordWrap(True); self.spliceRef.setTextFormat(Qt.RichText)
        card.bodylay.addWidget(self.spliceRef)
        # back-pointer to the record's annotation view (the provenance ladder: ANNOTATED gene model vs MEASURED splice)
        back = QPushButton("← The record's own annotation is the Gene model in panel 02")
        back.setProperty("link", True); back.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Fixed)   # shrink, never force the card wide
        back.clicked.connect(lambda: self._scroll_to(self.card_struct))
        card.bodylay.addWidget(back)
        # transcript picker — the record's OWN annotated mRNAs, so the user need not hunt for the matching transcript
        self.spliceTxRow = QWidget(); _tr = QHBoxLayout(self.spliceTxRow); _tr.setContentsMargins(0, 0, 0, 0)
        self.spliceTxSel = _Combo(); _tr.addWidget(self.spliceTxSel, 1)
        self.spliceTxUse = QPushButton("Use for splice"); self.spliceTxUse.setProperty("sm", True)
        self.spliceTxUse.clicked.connect(self._use_record_transcript); _tr.addWidget(self.spliceTxUse)
        card.bodylay.addWidget(self.spliceTxRow)
        self.spliceTxNote = QLabel(); self.spliceTxNote.setObjectName("cardmeta")
        self.spliceTxNote.setWordWrap(True); self.spliceTxNote.setTextFormat(Qt.RichText)
        card.bodylay.addWidget(self.spliceTxNote)
        self.spliceTx = QTextEdit()
        self.spliceTx.setPlaceholderText("Paste a transcript / cDNA / mRNA. minimap2 -x splice maps it to the loaded sequence; "
                                         "introns are the alignment gaps, checked against canonical GT–AG splice sites.")
        self.spliceTx.setMaximumHeight(round(80 * theme_mod.UI_SCALE))
        self.spliceTx.setTabChangesFocus(True)
        self._splice_tx_origin = "external"                   # a pasted/right-clicked transcript is treated as external
        self.spliceTx.textChanged.connect(self._on_splice_tx_changed)
        card.bodylay.addWidget(self.spliceTx)
        self.spliceBtn = QPushButton("Detect exons / introns"); self.spliceBtn.setProperty("primary", True)
        self.spliceBtn.setEnabled(False); self.spliceBtn.clicked.connect(self._splice)
        srow = QHBoxLayout(); srow.addWidget(self.spliceBtn); srow.addStretch(1)   # the panel's run action, at its own width
        card.bodylay.addLayout(srow)
        self.spliceBody = QVBoxLayout(); sb = QWidget(); sb.setLayout(self.spliceBody)
        self.spliceBody.addWidget(_empty("Align a transcript (mRNA / cDNA) to the loaded sequence to infer its exon–intron structure."))
        card.bodylay.addWidget(sb)
        return card

    def _send_to_splice(self, fasta):
        """Load a right-clicked subsequence as the transcript in the splice card and reveal it —
        mirrors 'send to in-silico PCR'. It is aligned to the loaded genomic reference."""
        self.spliceTx.setPlainText(fasta)
        self._pending_domain = None                       # a newer explicit navigation cancels any pending routed re-anchor
        self.card_splice.expand()
        try:
            self.resultsScroll.ensureWidgetVisible(self.card_splice)
        except Exception:
            pass
        self.spliceTx.setFocus()

    def _on_splice_tx_changed(self):
        # a manual paste / right-click send is an EXTERNAL transcript; only the record-transcript picker marks 'record'
        if not getattr(self, "_tx_programmatic", False):
            self._splice_tx_origin = "external"

    def _refresh_transcripts(self):
        """Populate the splice card's transcript picker with the record's OWN annotated mRNAs, or show an honest
        empty-state. Aligning a record's own transcript back to its locus is a consistency check, not independent
        confirmation — the note says so, and _use_record_transcript flags the origin for the cross-check relabel."""
        if not hasattr(self, "spliceTxSel"):
            return
        feats = self.state.get("features")
        txs = feats.get("transcripts") if isinstance(feats, dict) else None
        self.spliceTxSel.clear()
        if txs:
            for t in txs:
                self.spliceTxSel.addItem(t["accession"] + (f" · {t['product']}" if t.get("product") else ""), t["accession"])
            self.spliceTxRow.setVisible(True); self.spliceTxUse.setEnabled(True)
            self.spliceTxNote.setText("Transcripts annotated on this record. Using one aligns it back to its own "
                                      "locus — a <b>consistency check (same annotation source)</b>, not independent confirmation.")
        else:
            self.spliceTxRow.setVisible(False)
            self.spliceTxNote.setText("No annotated transcript on this specimen — paste a transcript / cDNA / mRNA "
                                      "below, or fetch its accession in panel 01.")

    def _update_splice_ref(self):
        """Show which specimen splice will align a transcript against (the genomic reference)."""
        if not hasattr(self, "spliceRef"):
            return
        self._refresh_transcripts()
        seq = self.state.get("seq") or ""
        if not seq:
            self.spliceRef.setText("Genomic reference: none loaded yet — load a specimen in panel 01.")
            return
        src = self.state.get("source") or {}
        n = len(self._clean_seq(seq))
        name = src.get("displayLocus") or src.get("accession")
        bound = None if name else self._record_binding_text()
        who = name or html.escape(bound or "pasted / uploaded specimen")
        org = f" · {src.get('organism')}" if src.get("organism") else ""
        self.spliceRef.setText(f"Genomic reference: <b>{who}</b>{org} · {n:,} bp "
                               f"({'selected in panel 02' if bound else 'from panel 01'})")

    # =================== splice detection ===================
    def _splice(self):
        self._clear_banner()                                  # drop any stale error before a retry (matches _design/_run_pcr)
        genomic = self.state.get("seq") or self.seq.toPlainText().strip()
        tx = self.spliceTx.toPlainText().strip()
        if not genomic.strip():
            return self._banner("Load a genomic sequence first (fetch, upload, or paste, then Run analysis).", "info")
        if not tx:
            return self._banner("Paste a transcript / cDNA / mRNA to align.", "info")
        genomic, source = self._operation_input(genomic)
        self.state["splice_seq"] = self._norm_seq(genomic)    # intron menus slice THIS (backend-normalized) sequence, not a later box edit
        self.spliceBtn.setEnabled(False); widgets.set_caps(self.spliceBtn, "◴ aligning…")
        self._set_body(self.spliceBody, BusyBar("Aligning the transcript to the genomic sequence (minimap2 -x splice)…"))
        self.engine.submit("splice", {"sequence": genomic, "transcript": tx, "source": source,
                                      "timeout": 300}, key=self._issue("splice", genomic))

    def _on_splice(self, d):
        self.spliceBtn.setEnabled(True); widgets.set_caps(self.spliceBtn, "Detect exons / introns")
        self.card_splice.expand()
        if not d.get("ok"):
            self._set_body(self.spliceBody, _note(d.get("error", "splice alignment failed")))
            return
        self._render_splice(d)
        if d.get("provenance"):
            self._render_provenance(d["provenance"])

    def _render_splice(self, d):
        self.state["splice"] = d
        cont = QWidget(); cl = QVBoxLayout(cont); cl.setContentsMargins(0, 0, 0, 0); cl.setSpacing(6)
        # two different facts: how the supplied transcript lies on this sequence, and which strand the gene is on
        aln, ts, gene = d.get("alignment_strand"), d.get("transcript_strand"), d.get("strand", "")
        orient = "as reverse-complement" if aln == "-" else "forward"
        if d.get("strand_determined") is False:
            strand_txt = (f"strand undetermined (the transcript aligned {orient}; no splice-site motif fixed the gene "
                          "strand, so the alignment orientation is shown)")
        elif aln and ts:
            strand_txt = (f"strand {gene} (the transcript aligned {orient}; its splice motifs read {ts} on the "
                          "transcript)")
        else:
            strand_txt = f"strand {gene}"
        head = QLabel(f"<b>{d['counts']['exons']} exon(s) · {d['counts']['introns']} intron(s)</b> — transcript alignment · "
                      f"{d.get('canonical_introns',0)}/{d['counts']['introns']} canonical splice site(s) · {strand_txt}"
                      f"{self._src_html('minimap2')}")
        head.setTextFormat(Qt.RichText); head.setWordWrap(True); _external(head)
        _kb_links(head); cl.addWidget(head)
        # independent cross-check: does this transcript-aligned (external-transcript) model agree with the record's
        # own annotation? Only meaningful when the genomic IS the fetched record (source accession present).
        src, feats = self.state.get("source"), self.state.get("features")
        if src and src.get("accession") and isinstance(feats, dict) and feats.get("introns"):
            # an intron inferred from the CDS gap is not an annotation, so it cannot be evidence here
            cc = cross_check_models([i for i in complete_gene_model(feats).get("introns", []) if not i.get("derived")],
                                    d.get("introns", []))
            if cc["annotation_total"]:
                extra = ((f" · {len(cc['aligned_only'])} alignment-only" if cc["aligned_only"] else "")
                         + (f" · {len(cc['annotation_only'])} annotation-only" if cc["annotation_only"] else ""))
                # CIRCULARITY GATE (content-based): a transcript that IS one of THIS record's own annotated transcripts,
                # aligned back to its locus, is a consistency check by construction — regardless of HOW it was loaded
                # (the picker, a separate accession fetch, or a paste). Never call that 'independent confirmation'.
                same_source = getattr(self, "_splice_tx_origin", "external") == "record"
                undetermined = False
                if not same_source:
                    _txt = self.spliceTx.toPlainText().lstrip()
                    _hdr = _txt.splitlines()[0] if _txt.startswith(">") else ""
                    _m = re.search(r"[A-Z]{2}_\d+", _hdr)
                    if _m:
                        _own = {re.sub(r"\..*", "", t.get("accession", "")) for t in (feats.get("transcripts") or [])}
                        same_source = _m.group(0) in _own
                    else:
                        undetermined = True                   # a headerless paste — provenance cannot be established
                if same_source:
                    note = QLabel(f"<b>Consistency check vs {src['accession']} annotation (same annotation source):</b> "
                                  f"{cc['matched']}/{cc['annotation_total']} intron(s) reproduced{extra}. "
                                  "This transcript is the record's OWN annotated mRNA aligned back to its locus, so "
                                  "agreement confirms internal consistency — it is NOT independent confirmation — and "
                                  "is not part of the sealed result.")
                elif undetermined:
                    note = QLabel(f"<b>Cross-check vs {src['accession']} annotation (advisory):</b> "
                                  f"{cc['matched']}/{cc['annotation_total']} intron(s) reproduced{extra}. "
                                  "Transcript provenance undetermined (no accession header) — if this is the record's own "
                                  "annotated mRNA, agreement is a same-source consistency check, not independent confirmation. "
                                  "Not part of the sealed result.")
                else:
                    note = QLabel(f"<b>Cross-check vs {src['accession']} annotation (advisory):</b> "
                                  f"{cc['matched']}/{cc['annotation_total']} intron(s) confirmed{extra}. "
                                  "Independent comparison of this alignment against the record's current annotation — "
                                  "not part of the sealed result (the annotation may be revised).")
                note.setTextFormat(Qt.RichText); note.setWordWrap(True); note.setObjectName("orient"); cl.addWidget(note)
                if cc["matched"] == 0:                        # a genomic slice pasted as the 'transcript' aligns gaplessly
                    hint = QLabel("No introns matched — if you aligned a genomic subsequence rather than a spliced "
                                  "transcript / mRNA / cDNA, that is expected (a genomic slice has no spliced-out introns).")
                    hint.setWordWrap(True); hint.setObjectName("orient"); cl.addWidget(hint)
        # measure the sequence splice ACTUALLY ran on — intron coords are in splice_seq's frame. last_rec is a
        # prior analysis that _seq_changed never clears, so preferring it drew an old element's ruler under new results.
        length = len(self.state.get("splice_seq") or "") or len(self._clean_seq(self.state.get("seq", ""))) or 1
        gene_strand = d.get("strand") if d.get("strand") in ("+", "-") else "+"
        model = figures.gv_tracks_from_gene(
            {"exons": [{"strand": gene_strand, **e} for e in d.get("exons", [])],
             "introns": [{"strand": gene_strand, **i} for i in d.get("introns", [])], "cds": []}, length)
        if model["tracks"]:
            gv = CoordGenomePanel(svg_genome, "TEagle_splice"); gv.apply_app_theme(self.theme)
            gv.coord_base = self.coord_base; gv.set_model(model)
            # coordinates index the sequence the alignment ran on, which exists even before any analysis
            gv.set_feature_menu(lambda region, s=self.state.get("splice_seq") or None: self._region_menu(region, s))
            gv.setMinimumHeight(round(220 * theme_mod.UI_SCALE))
            cl.addWidget(gv)
        introns = d.get("introns", [])
        if introns:
            headers = ["#", self._ch("Intron span"), "Len", "Splice site", "Canonical"]
            irow = lambda k, i, b: [k + 1, self._span(i["start"], i["end"], b), i["end"] - i["start"], f"{i['donor']}…{i['acceptor']}",
                                    "canonical" if i.get("canonical") else "non-canonical"]
            t = self._coord_table(headers, self._gloss(), [irow(k, i, 0) for k, i in enumerate(introns)])
            t.export_header_map = {"Intron span (1-based)": "Intron span (0-based)"}
            t.export_notes = ["Introns inferred by aligning the supplied transcript to the genomic sequence (minimap2 -x splice). Coordinates are 0-based "
                              "half-open on the genomic sequence supplied to this run; canonical labels report "
                              "the observed donor and acceptor motif, not independent transcript confirmation."]
            t.set_rows([irow(k, i, self.coord_base) for k, i in enumerate(introns)])
            t.set_row_menu(lambda r: self._feat_menu(introns[r]["start"], introns[r]["end"], d.get("strand", "+"),
                                                     f"intron_{r+1}", src_seq=self.state.get("splice_seq")))
            cl.addWidget(t)
            xrow = QHBoxLayout(); xrow.addStretch(1)
            xrow.addWidget(_export_table_btn(t, "TEagle_splice_introns", self))
            cl.addLayout(xrow)
        else:
            cl.addWidget(_empty("single exon — no introns detected"))
            # foot-gun guard for the common novice case (works for any input, not only fetched+annotated records):
            # a gapless alignment reads as a real single-exon finding but usually means a genomic slice was pasted.
            slice_note = QLabel("A gapless alignment (0 introns) is consistent with either a genuine single-exon "
                                "transcript OR a genomic slice pasted into the transcript box (a genomic slice has no "
                                "spliced-out introns). To resolve splicing, align an mRNA / cDNA / EST, not genomic DNA.")
            slice_note.setWordWrap(True); slice_note.setObjectName("orient"); cl.addWidget(slice_note)
        if d.get("exons") and self.state.get("splice_seq"):
            splice_record = self.state.get("last_rec") or {}
            splice_identity = (_record_public_identity(splice_record) if splice_record else
                               {"id": "splice_locus", "identitySource": "generated"})
            rid = splice_identity["id"].split()[0]
            rid_sha256 = splice_identity.get("idSha256")
            rid_source = splice_identity["identitySource"]
            erow = QHBoxLayout(); erow.addStretch(1)
            erow.addWidget(_export_fasta_btn(
                "Export splice exons (FASTA)",
                lambda result=d, sequence=self.state.get("splice_seq"), record_id=rid, digest=rid_sha256,
                       source=rid_source: _splice_exons_fasta(result, sequence, record_id, digest, source),
                "TEagle_splice_exons", self))
            model_export = QPushButton("Export splice model ▾"); model_export.setProperty("sm", True)
            model_export.setProperty("resultExport", True)
            model_export.setToolTip("Export the displayed transcript-aligned exon/intron model as GFF3, BED, or exon FASTA")
            model_export.clicked.connect(
                lambda _=False, result=d, sequence=self.state.get("splice_seq"), record_id=rid,
                       digest=rid_sha256, source=rid_source:
                self._export_gene_model({"exons": result.get("exons") or [],
                                         "introns": result.get("introns") or [], "cds": [],
                                         "strand": result.get("strand", ".")},
                                        sequence, record_id, "splice_alignment",
                                        identity_source=source, record_id_sha256=digest))
            erow.addWidget(model_export)
            cl.addLayout(erow)
        self._set_body(self.spliceBody, cont)
