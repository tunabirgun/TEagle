"""Record table, record summary and annotation exports."""
from __future__ import annotations

import json
import os
from PySide6.QtCore import Qt
from PySide6.QtGui import QCursor
from PySide6.QtWidgets import QDialog, QFileDialog, QHBoxLayout, QLabel, QMessageBox, QVBoxLayout
import widgets
from coords import CoordTable
import theme as theme_mod
from teagle_core import sequtil, domains as domains_mod, __version__ as APP_VERSION
from ui_common import (
    GLOSS, _MIN_ORF_AA, _call_view, _empty, _export_fasta_btn, _export_table_btn, _full_records_fasta,
    _no_orf_to_search, _panel_references, _rank_view, _record_export_stem, _record_public_id,
    _record_public_identity, _safe_export_stem, _te_headline)


class RecordsMixin:
    # engine.analyze already classifies EVERY record; the UI used to take records[0] and discard the rest,
    # so a user arriving with a shortlist saw one answer and no indication the others had been computed.
    _MAX_RECORDS = 200

    def _on_record_row(self, orig):
        """The user moved to another row of the record table: bind that record once the movement settles."""
        if orig < 0 or orig == self.state.get("record_index"):
            self._recTimer.stop()                             # back on the bound row: a pending switch is withdrawn
            return
        self._rec_pending = orig
        self._recTimer.start()

    def _apply_record_selection(self):
        if self.state.get("analysis_state") == "current" and self._rec_pending is not None:
            self._show_record(self._rec_pending, keep_table=True)

    def _record_binding_text(self):
        """'record 2 of 3 · <id>' (or just the id for one record) for the record every downstream step uses; None
        when no analysis is current."""
        recs = self.state.get("records") or []
        if self.state.get("analysis_state") != "current" or not recs or self.state.get("record_index") is None:
            return None
        rec = recs[min(self.state["record_index"], len(recs) - 1)]
        rid = _record_public_id(rec)
        return f"record {self.state['record_index'] + 1} of {len(recs)} · {rid}" if len(recs) > 1 else rid

    def _sync_record_binding(self):
        """Name the bound record on every surface that acts on it: the rail tiles, the panel 04/05 headers and the
        Design hint. With no current analysis they return to their defaults."""
        b = self._record_binding_text()
        n = len(self.state.get("records") or [])
        self.mScope.setVisible(bool(b) and n > 1)
        if b and n > 1:
            self.mScope.setText(f"Selected {b} — the readouts below describe this record only.")
        self.card_primer.set_meta(f"designs PCR primers (Primer3) · template: {b}" if b else None)
        self.card_pcr.set_meta(f"predicts which fragments amplify · template: {b}" if b else None)
        valid = (self.state.get("last_rec") or {}).get("valid", True)
        if not b:
            self.designHint.setText("run analysis first")
        elif valid:
            self.designHint.setText(f"template: {b}")
        self._refresh_assay_goal()
        self._update_splice_ref()

    def _show_record(self, index: int, keep_table=False):
        """Render one record's detail. Every downstream step (primers, PCR, splice) stays bound to the
        selected record, so the choice is explicit rather than implied by file order."""
        res = self.state.get("analysis") or {}
        recs = self.state.get("records") or []
        if not recs:
            return
        index = max(0, min(index, len(recs) - 1))
        rec = recs[index]
        self.state["last_rec"] = rec
        self.state["record_index"] = index
        bound = ((res.get("provenance") or {}).get("input") or {}).get("records") or []
        self.state["selected_record"] = (dict(bound[index], index=index) if index < len(bound)
                                          else {"id": rec.get("id"), "index": index,
                                                "length": len(rec.get("seq", ""))})
        # Bind every downstream step (primer design, in-silico PCR, splice, annotate, feature slicing) to
        # THIS record's own sequence. state["seq"]/["analyzed_seq"] otherwise hold the raw textbox, which
        # for a multi-FASTA paste is all records concatenated — so selecting a non-first record would
        # silently run design/PCR/annotate/splice on record 1 (engine stores each record's seq since v3.2.0).
        if rec.get("seq"):
            self.state["seq"] = rec["seq"]
            self.state["analyzed_seq"] = rec["seq"]        # _slice() keys off this -> per-record feature coords
            # NB: do NOT touch state["analyzed_clean"] here. It is the whole-textbox snapshot _stale_block()
            # compares the live box against; overwriting it with one record's bases made the box (all records)
            # always differ from it, so the "sequence changed" guard permanently blocked Design/PCR for every
            # multi-FASTA analysis. The stale guard tracks box edits; per-record targeting rides on seq/analyzed_seq.
        self._recTimer.stop()
        usable = bool(rec.get("valid", True))
        self.designBtn.setEnabled(usable and not self._design_inflight)
        self._sync_record_binding()
        for kind in ("primers", "annotate", "splice"):        # a request bound to the record just left can only return as stale: free its control now
            req = self._reqs.get(kind)
            if req and req["bound"] and req["rec"] not in (None, index):
                del self._reqs[kind]
                self._release(kind)
        run = self._pcr_run
        if run and run.get("rec") not in (None, index) and any(r is None for r in run["results"]):
            self._pcr_run = None; self._pcr_gen += 1          # a batch run on the record just left would otherwise land as this record's
            widgets.set_caps(self.runPcrBtn, "Run loaded pairs"); self.runPcrBtn.setEnabled(bool(self.state.get("pcrPairs")))
        lp = self.state.get("lastPcr")
        if lp and lp.get("digest") not in (None, self._selected_digest()):
            self.state.pop("lastPcr", None)                   # a finished run on another sequence is not this record's result
            self._set_body(self.pcrBody, _empty("The previous in-silico PCR result was for another record. Run the loaded "
                                                "pairs on this record to see its products."))
            self._render_pcr_queue()
        if not usable:
            self.designHint.setText("This record was not analysed (invalid input); primers cannot be "
                                    "designed on it. Select a valid record.")
        if not self.state.get("lastPrimers"):                  # the "run analysis first" placeholder is now wrong
            self._set_body(self.primBody, _empty("Ready — choose a preset and press Design primers, or select a domain or ORF row and use "
                                                 "Design primers here to target it."))
        if not (self.state.get("lastPcr") or self.state.get("pcrPairs")):
            self._set_body(self.pcrBody, _empty("Load a pair, then run pair-aware amplicon search."))
        comp = rec.get("composition", {})
        self.mLen.setText(f"{comp.get('length', 0):,}")
        self.mGC.setText(f"{comp.get('gc', 0)}%")
        self.mN.setText(f"{comp.get('n', 0)}%")
        # "valid" alone was read as "clean, unambiguous DNA": a sequence a third R/Y/W is perfectly valid
        # IUPAC and showed "valid" beside "N content 0.0%", because degenerate bases are counted by neither.
        # Say how much of it is degenerate whenever any is, so the reading is not more confident than the
        # sequence. The gauge stays "valid"/"invalid" — the verdict is unchanged, only its scope is stated.
        amb = comp.get("ambiguous", 0.0)
        ok = rec.get("valid")
        if ok and amb:                                        # the share on its own small line: a long value widened its column and clipped mid-word
            self.mValid.setTextFormat(Qt.RichText)
            self.mValid.setText(f"valid<br><span style='font-size:11px'>{amb}% ambiguous</span>")
        else:
            self.mValid.setTextFormat(Qt.PlainText)
            self.mValid.setText("valid" if ok else "invalid")
        if amb:
            codes = ", ".join(f"{b}×{c}" for b, c in sorted((comp.get("ambiguous_counts") or {}).items()))
            self.mValid.setToolTip(f"IUPAC degenerate bases present ({codes}). They are excluded from the GC "
                                   f"denominator and are not counted as N.")
        else:
            self.mValid.setToolTip("")
        self.mValid.setProperty("state", "good" if ok else "bad"); self._repolish(self.mValid)
        self._set_trace_counts(recs, rec)
        self._set_char_count(self.seq.toPlainText())          # reads "nt" now that a record has validated
        self._render_struct_card(rec, res, keep_table=keep_table)
        if keep_table and self._recTable is not None:
            self._recTable.select_orig(index)                 # a menu or code route moves the highlight too (silent)
        if self.state.get("lastPrimers") and not self._design_inflight:
            self._render_primers(self.state["lastPrimers"])   # its "designed on" line is re-checked against the record now shown
        self._uppercase_buttons()

    def _selected_digest(self):
        seq = (self.state.get("last_rec") or {}).get("seq")
        return sequtil.sha256(seq) if seq else None

    def _record_summary_rows(self, recs, limit=None):
        """One row per record: what was actually computed for each, not a promise that it was.

        `limit=None` means every record. The on-screen table passes _MAX_RECORDS because a 200-row widget is
        already unwieldy, but the EXPORT passes None — engine.analyze classified all of them, so capping the
        file would drop computed results without saying so."""
        from teagle_core import gff3
        rows = []
        for i, r in enumerate(recs if limit is None else recs[:limit]):
            cl = r.get("classification") or {}
            comp = r.get("completeness") or cl.get("completeness") or {}
            doms = sorted({d["domain"] for d in (r.get("domains") or [])})
            tdr = gff3.is_neutral_tdr(cl)
            struct = sorted({("TDR" if tdr and e["type"].startswith("LTR") else e["type"].split(" (")[0])
                             for e in (r.get("structural") or [])})
            domain_status, unsearched = self._scan_state(r, cl)
            if unsearched:
                domain_status += f" (partial: {unsearched} ORFs not searched)"
            wicker, arch = self._orders(cl, domain_status)
            processing = "; ".join(f"{stage}={state.get('status', '—')}"
                                   for stage, state in (r.get("processing") or {}).items()) or "—"
            composition = r.get("composition") or {}
            call = _call_view(r)
            rank = _rank_view(cl)
            row = {"#": i + 1, "Record": _record_public_id(r), "Classification": _te_headline(cl.get("te_class") or "—"),
                   "Call confidence": call["cell"], "Completeness": comp.get("tier") or "—", "Call status": call["status"],
                   "Wicker order": wicker, "Domain architecture": arch, "Domain evidence": domain_status,
                   "Wicker code": rank["code"] if rank else "—", "RM class": (rank or {}).get("rm") or "—",
                   "Length": composition.get("length", 0),                       # the number itself: it sorts and exports as one
                   "GC": f"{composition.get('gc', 0)}%", "N": f"{composition.get('n', 0)}%",
                   "Ambiguous": f"{composition.get('ambiguous', 0)}%", "Input": "valid" if r.get("valid") else "INVALID",
                   "Domains": "·".join(doms) or "—", "Structural": ", ".join(struct) or "—",
                   "Rationale": cl.get("explanation") or (f"not analysed: {r['invalidReason']}" if r.get("invalidReason")
                                                          else "not analysed: the sequence is not valid nucleotide input"
                                                          if not r.get("valid", True) else "—"),
                   "Processing": processing}
            rows.append([row[h] for h in self._RECORD_HEADERS])
        return rows

    @staticmethod
    def _scan_state(r, cl):
        """(domain-search status, ORFs left unsearched). 'not_run' when no search happened at all (invalid input,
        or no ORF of 40 aa or more): that is a different statement from 'tested_empty'."""
        from teagle_core import gff3
        if not cl or ((r.get("processing") or {}).get("domains") or {}).get("status") == "not_run":
            return "not_run", 0
        return gff3.domain_scan_state(cl)

    @staticmethod
    def _orders(cl, domain_status):
        """(Wicker order, domain architecture) for the summary, kept apart because the two are different things:
        the first is a rank in the classification scheme, the second the detected domain order."""
        if not cl:
            return "not analysed", "not analysed"
        wicker = cl.get("wicker_order") or ("withheld" if cl.get("outcome") == "withheld" else "unassigned")
        arch = cl.get("domain_architecture") or cl.get("order")
        if domain_status == "not_run":
            arch = "not assessable (no ORF ≥ 40 aa to search)"
        elif cl.get("domains_unavailable"):
            arch = "not assessed (domain scan failed)"
        elif not arch or arch == "unavailable":
            arch = "none detected"
        return wicker, arch

    @staticmethod
    def _profiles_by_code():
        """{reported domain code: [profile names]} grouped from the panel the scan actually loads."""
        out = {}
        for hmm, (code, _label, _cls, _pfam) in domains_mod.DOMAIN_INFO.items():
            out.setdefault(code, []).append(hmm)
        return out

    def _module_lines(self, rec, cl, comp):
        """Which expected modules were found, and what 'not detected' is relative to. A 'partial' tier is a
        verdict on the bundled panel: for copia and gypsy it cannot match a Ty1/Ty3-type gag or protease at all,
        so the profiles that WERE tested are listed (derived from the panel, never retyped)."""
        expected = comp.get("expected") or []
        if not expected:
            return ""
        profiles = self._profiles_by_code()
        present, unassessed = set(comp.get("present") or []), set(comp.get("unassessed") or [])
        has_tir = any(e["type"].startswith("TIR") for e in rec.get("structural") or [])
        states = []
        for m in expected:
            if m in present:
                s = "found"
            elif cl.get("domains_unavailable") and m in profiles:
                s = "not searched (domain scan failed)"
            elif m in unassessed:
                s = "not assessed (structural scan did not run)"
            elif m.startswith("TIR") and has_tir:
                s = "found but not credited (the pair does not enclose the transposase)"
            else:
                s = "not detected"
            states.append((m, s))
        scan = (rec.get("processing") or {}).get("domains") or {}
        scope = (f" — domains searched in the {scan.get('orfsScanned')} longest of {len(rec.get('orfs') or [])} ORFs"
                 if scan.get("orfsNotScanned") else "")
        out = (f"<b>Expected modules ({sum(s == 'found' for _, s in states)} of {len(states)} found)</b>{scope}: "
               + " · ".join(f"{m} {s}" for m, s in states))
        absent = [m for m, s in states if s == "not detected" and m in profiles]
        if absent:
            out += (f"<br>“Not detected” is relative to the bundled {len(domains_mod.DOMAIN_INFO)}-model Pfam panel: "
                    "a divergent or unmodelled domain reads as not detected, not as element decay. Profiles tested "
                    "for the modules not detected — "
                    + "; ".join(f"{m}: {', '.join(profiles[m])}" for m in absent) + ".")
        tier = comp.get("tier") or ""
        ref = _panel_references().get(cl.get("te_class"))
        if absent and ref and ref[3] == tier and not tier.startswith("intact"):
            acc, label, organism, _t, missing, digest = ref   # measured: the bundled reference of this class reads the same
            out += ("<br>" + (f"This record is the bundled {label} reference ({acc}): the panel reports “{tier}” even for "
                              "this published element, so the tier is not by itself evidence of decay."
                              if digest == sequtil.sha256(rec.get("seq") or "") else
                              f"The bundled {label} reference ({acc}, {organism}) gives the same result with this panel "
                              f"({tier}; {', '.join(missing)} not detected), so this tier is not by itself evidence that the "
                              "copy is decayed."))
        return out

    # The columns that qualify the call come right after the record and its length: on a 1366 px window the old order
    # left Confidence, Completeness and Domain architecture off-screen, so a reader saw 'complete' beside a class and
    # took it for element completeness. Length stays third (the column tests and exports address it there).
    _RECORD_HEADERS = ["#", "Record", "Length", "Classification", "Call confidence", "Completeness", "Call status",
                       "Wicker order", "Domain architecture", "Domain evidence", "Wicker code", "RM class", "GC", "N",
                       "Ambiguous", "Input", "Domains", "Structural", "Rationale", "Processing"]
    _SUMMARY_NOTES = [
        "Call status: called = GFF3 classification_outcome=complete; withheld = withheld; failed and not run have no classification.",
        "Domain evidence uses plain words here: hits found = tested_with_hits; none found in panel = tested_empty; scan failed = "
        "unavailable_failed; not assessable = no ORF of at least 40 aa was available to search. GFF3 and JSON keep the "
        "machine tokens.",
        "Call confidence is n/a where the call was withheld or no TE signature was found; GFF3 and JSON keep the classifier's "
        "raw confidence field."]
    _DOMAIN_SEARCH_WORDS = {"tested_with_hits": "hits found", "tested_empty": "none found in panel",
                            "unavailable_failed": "scan failed (not assessed)"}

    def _summary_table_rows(self, recs):
        """The summary rows as the table shows and exports them: the stage tokens of the Domain evidence column spelled
        out in plain words. _record_summary_rows keeps the tokens (they are the GFF3/JSON vocabulary)."""
        col, status_col = self._RECORD_HEADERS.index("Domain evidence"), self._RECORD_HEADERS.index("Call status")
        out = []
        for rec, row in zip(recs, self._record_summary_rows(recs, limit=None)):
            row = list(row)
            token, _, rest = str(row[col]).partition(" (")
            word = (self._DOMAIN_SEARCH_WORDS.get(token) or ("not run" if row[status_col] == "not run" else
                    f"not assessable (no ORF ≥ {_MIN_ORF_AA} aa)" if _no_orf_to_search(rec) else
                    "not assessed (ORF detection failed)"))
            row[col] = word + (f" ({rest}" if rest else "")
            out.append(row)
        return out

    _RESULT_GLOSS = {
        "Wicker order": "Rank in the Wicker 2007 scheme (LTR, LINE, TIR …). 'withheld' means TEagle declined to name "
                        "it from this evidence; 'unassigned' means no TE order applies.",
        "aa in ORF": "Residue range of the match inside the ORF's translated protein, 1-based with both ends "
                     "included (HMMER's envelope). It is a position, not the protein's length; the nt column is the "
                     "same hit on the record, 0-based half-open.",
        "ORF": "Number of the ORF this domain lies in, as listed in the ORF table; '—' when the hit cannot be placed "
               "on one ORF.",
        "Domains": f"Protein domains detected inside this ORF. Only the {domains_mod.MAX_ORFS_SCANNED} longest ORFs are "
                   "searched; 'not searched' marks the rest and 'not assessed' means the search did not run.",
        "Domain evidence": "What the protein-domain search did on this record: tested_with_hits (hits found), tested_empty "
                           "(searched, nothing in the panel), unavailable_failed (the search broke — not an absence), or "
                           "not_run (nothing was searched, e.g. invalid input or no ORF of 40 aa or more). A partial search "
                           "names the ORFs it skipped. The same words are used in the GFF3 and JSON exports.",
    }

    def _build_selfsim(self, rec):
        """Open the self-similarity plot for the analysed record in its own resizable window.

        A dedicated window (rather than an inline card panel) keeps a near-square matrix from dominating
        the results column, and gives the plot room to be zoomed and panned. The word size follows the
        shortest terminal repeat TEagle actually measured: a repeat SHORTER than the word size cannot
        produce a single exact match, so a fixed k=13 hides an 11 bp hAT TIR entirely (measured on maize
        Ac: reverse signal 2 at k=13, 600 at k=8)."""
        from teagle_core import dotplot
        # this record's OWN sequence only — never fall back to analyzed_clean, which for a multi-FASTA paste
        # concatenates every record and would plot a chimeric matrix. A record with an empty sequence
        # (malformed FASTA: back-to-back headers) instead hits the empty-seq banner below.
        seq = rec.get("seq") or ""
        if not seq:
            return self._banner("Run an analysis first — the self-similarity plot needs the analysed sequence.")
        struct = rec.get("structural") or []
        k = dotplot.suggest_k(struct)
        m = dotplot.self_matrix(seq, k=k)
        # only pairs the detector accepted are "detected terminal repeats"; a rejected advisory pair stays
        # off the plot, and LTR and TIR bands take their own colours
        guides = [{"start": e[q][0], "end": e[q][1],
                   "color": theme_mod.OK.get("TIR" if e["type"].startswith("TIR") else "LTR", "#0072B2")}
                  for e in struct if not e.get("advisory") and not e["type"].startswith("Sub-threshold")
                  for q in ("five_prime", "three_prime") if e.get(q)]
        rid = _record_public_id(rec).split()[0]
        # ONE reused window (hidden, not destroyed, on close) whose body is re-populated per record — never a
        # rebuilt widget, which is what collapsed a reused dialog in an earlier release.
        dlg = getattr(self, "_selfsim_dlg", None)
        if dlg is None:
            dlg = QDialog(self)
            dlg.setObjectName("selfsimdlg")
            dlg.resize(940, 820)
            _l = QVBoxLayout(dlg)
            _l.setContentsMargins(*(theme_mod.sp(12),) * 4)
            dlg._panel = widgets.DotPanel(base_name="TEagle_selfsim", parent=dlg)
            _l.addWidget(dlg._panel)
            self._selfsim_dlg = dlg
            self._uppercase_buttons()
        panel = dlg._panel
        dlg.setWindowTitle(f"Self-similarity — {rid}")
        panel.base_name = f"TEagle_{_record_export_stem(rec, 'locus')}_selfsim"
        panel.apply_app_theme(self.theme)
        panel.zoom = None                                 # fit the new record to the window
        panel.set_matrix(m, guides=guides, threshold=dotplot.above_chance(m), scope=dotplot.scope_note(m))
        dlg.show(); dlg.raise_(); dlg.activateWindow()
        panel._fit()
        self._uppercase_buttons()

    def _bound_coord_region(self, rec):
        """Return the coordinate region only when it is still bound to this exact record and source."""
        reg = self.state.get("coord_region") or {}
        rec_id = _record_public_id(rec)
        rec_seq = rec.get("seq") or ""
        source = self.state.get("source") or {}
        if source.get("inputType") != "fetched_coordinate":
            return {}
        from teagle_core import provenance
        if provenance.validate_source_context(source, [(rec_id, rec_seq)]) != source:
            return {}
        try:
            source_binding = json.dumps(source, sort_keys=True, ensure_ascii=False,
                                        allow_nan=False, default=str)
        except (TypeError, ValueError):
            return {}
        if (reg.get("chrAccession") and reg.get("strand") != 2 and rec_seq and
                reg.get("_record_id") == rec_id and
                reg.get("_sequence_sha256") == sequtil.sha256(rec_seq) and
                reg.get("_source_binding") == source_binding):
            return reg
        return {}

    def _export_annotation(self, rec):
        """Write the record's annotation as GFF3 or BED.

        Coordinates are LOCUS-RELATIVE and the seqid is the record's own id: exporting locus offsets
        under a chromosome name would silently misplace every feature in a browser. GFF3 embeds the
        sequence after ##FASTA so the file stands alone."""
        from PySide6.QtWidgets import QMenu
        from teagle_core import gff3
        from teagle_core import __version__ as teagle_ver
        teagle_version = lambda: teagle_ver
        menu = QMenu(self)
        a_gff = menu.addAction("GFF3 (.gff3) — ontology terms + sub-features + sequence")
        a_bed = menu.addAction("BED (.bed) — intervals only")
        # Genome-anchored variants, offered ONLY when this record demonstrably came from a coordinate fetch
        # of the plus strand. A locus-relative file intersected against a chromosome-space track returns
        # zero overlaps and says nothing about it, which is the failure a bioinformatician hits first.
        stale_reg = self.state.get("coord_region") or {}
        reg = self._bound_coord_region(rec)
        rec_id = _record_public_id(rec)
        anchorable = bool(reg)
        a_gff_g = a_bed_g = None
        if anchorable:
            menu.addSeparator()
            # the seqid written is the RefSeq accession; a UCSC-named track (chr1) will not match it by name
            named = reg["chrAccession"] + (f" (UCSC name {reg['chromLabel']}; not written)"
                                           if reg.get("chromLabel") and reg["chromLabel"] != reg["chrAccession"] else "")
            a_gff_g = menu.addAction(f"GFF3, genome-anchored on {named}")
            a_bed_g = menu.addAction(f"BED, genome-anchored on {named}")
        elif stale_reg.get("strand") == 2:
            # a minus-strand fetch is reverse-complemented, so locus coordinates run against the chromosome.
            # Mapping them back is a flip, not an offset — refuse rather than ship a subtly wrong file.
            a_note = menu.addAction("Genome-anchored export unavailable (minus-strand fetch)")
            a_note.setEnabled(False)
        chosen = menu.exec(QCursor.pos())
        if chosen is None:
            return
        anchored = chosen in (a_gff_g, a_bed_g) and anchorable
        fmt = "gff3" if chosen in (a_gff, a_gff_g) else "bed"
        seqid = rec_id.split()[0]
        # 0-based genome coordinate of the locus's first base; the fetch reports 1-based inclusive starts
        offset = int(reg.get("start", 1)) - 1 if anchored else 0
        if anchored:
            seqid = reg["chrAccession"]
        base_stem = (_safe_export_stem(seqid, "locus", trust="validated_fetch") if anchored else
                     _record_export_stem(rec, "locus"))
        base = f"TEagle_{base_stem}"
        path, _ = QFileDialog.getSaveFileName(self, f"Export {fmt.upper()}", base + "." + fmt,
                                              f"{fmt.upper()} (*.{fmt})")
        if not path:
            return
        if not path.lower().endswith("." + fmt):
            path += "." + fmt
        # the record's OWN sequence (not analyzed_clean, which concatenates every record of a multi-FASTA
        # paste — embedding that after ##FASTA would contradict this record's own coordinates)
        seq = rec.get("seq") or None
        span = ((int(reg["start"]), int(reg["stop"])) if anchored else None)
        identity = ({"id": seqid, "identitySource": "validated_fetch"} if anchored else
                    _record_public_identity(rec))
        writer_trust = identity["identitySource"]
        text = (gff3.to_gff3(rec, seqid=seqid, sequence=seq,
                             source_note=f"TEagle {teagle_version()}", offset=offset, region_span=span,
                             identity_trust=writer_trust, record_id_sha256=identity.get("idSha256"))
                if fmt == "gff3" else gff3.to_bed(
                    rec, seqid=seqid, offset=offset, identity_trust=writer_trust,
                    record_id_sha256=identity.get("idSha256")))
        def write(p):
            with widgets.atomic_write(p, "w", encoding="utf-8", newline="\n") as fh:
                fh.write(text)
        widgets.safe_export(self, path, write, ok=f"Annotation written to {os.path.basename(path)}.")

    def _export_gene_model(self, model, sequence, record_id, model_kind, export_format=None,
                           identity_source=None, record_id_sha256=None):
        """Export the exact fetched/splice model shown by its GenomePanel through the shared writers."""
        from PySide6.QtWidgets import QMenu
        from teagle_core import gff3

        extension = export_format
        if extension is None:
            menu = QMenu(self)
            actions = {
                menu.addAction("GFF3 (.gff3) — exon/intron/CDS coordinates + source sequence"): "gff3",
                menu.addAction("BED (.bed) — intervals only"): "bed",
                menu.addAction("Exon FASTA (.fasta)" if model_kind == "splice_alignment"
                               else "Exon/CDS FASTA (.fasta)"): "fasta",
            }
            extension = actions.get(menu.exec(QCursor.pos()))
            if extension is None:
                return
        if extension not in {"gff3", "bed", "fasta"}:
            raise ValueError(f"unsupported model export format: {extension}")
        identity = _record_public_identity({"id": record_id, "identitySource": identity_source,
                                            "idSha256": record_id_sha256})
        writer_trust = identity["identitySource"]
        record_id = identity["id"]
        stem_trust = "generated" if identity.get("idSha256") else writer_trust
        stem = gff3.safe_filename_stem(record_id, trust=stem_trust)
        base = f"TEagle_{stem}_{'splice_model' if model_kind == 'splice_alignment' else 'gene_model'}"
        path, _ = QFileDialog.getSaveFileName(
            self, f"Export {extension.upper()}", base + "." + extension,
            f"{extension.upper()} (*.{extension})")
        if not path:
            return
        if not path.lower().endswith("." + extension):
            path += "." + extension
        try:
            if extension == "gff3":
                text = gff3.gene_model_to_gff3(model, record_id, sequence=sequence,
                                                source=f"TEagle {APP_VERSION}", model_kind=model_kind,
                                                identity_trust=writer_trust,
                                                record_id_sha256=identity.get("idSha256"))
            elif extension == "bed":
                text = gff3.gene_model_to_bed(model, record_id, sequence_length=len(self._norm_seq(sequence)),
                                               source=f"TEagle {APP_VERSION}", model_kind=model_kind,
                                               identity_trust=writer_trust,
                                               record_id_sha256=identity.get("idSha256"))
            else:
                feature_types = ("exon",) if model_kind == "splice_alignment" else ("exon", "CDS")
                text = gff3.gene_model_to_fasta(
                    model, sequence, record_id, feature_types=feature_types,
                    identity_trust=writer_trust, record_id_sha256=identity.get("idSha256"))
            with widgets.atomic_write(path, "w", encoding="utf-8", newline="\n") as fh:
                fh.write(text)
        except (OSError, ValueError) as exc:
            self._banner(f"The {extension.upper()} model export was not written: {exc}", level="warn")
            return
        self._banner(f"Model export written to {os.path.basename(path)}.", level="success")

    def _render_record_table(self, card):
        """Common summary for one or many records, above the selected record's detail."""
        recs = self.state.get("records") or []
        if not recs:
            return
        rows = self._summary_table_rows(recs)
        labels = self._local_labels(recs)
        heads, lab_col, tips = list(self._RECORD_HEADERS), None, None
        if any(labels):                                       # on screen only: the header text is never a file's content
            lab_col = heads.index("Record") + 1
            heads.insert(lab_col, "Your label")
            rows = [r[:lab_col] + [(lab.split()[0] if lab else "—")] + r[lab_col:] for r, lab in zip(rows, labels)]
            tips = [[None] * lab_col + [f"Your full header (shown on this screen only): {lab}" if lab else None]
                    + [None] * (len(heads) - lab_col - 1) for lab in labels]
        t = CoordTable(heads, {**GLOSS, **self._RESULT_GLOSS, "Domains": "Distinct domain codes detected, alphabetical.",
                               "Record": ("Neutral record ID used on screen and in every export. A pasted or uploaded header is replaced "
                                          "by record_ plus the first 12 hex digits of the SHA-256 of the whole header line, so sample "
                                          "names, file paths or patient codes never reach an exported file. It is a pseudonym, not "
                                          "encryption: anyone who already knows the exact header can recompute it. Accessions fetched "
                                          "through TEagle stay readable."),
                               "Your label": ("The first word of your FASTA header, shown on this screen only. It is never written to an "
                                              "exported table, FASTA, GFF3/BED, figure, filename or provenance manifest — those use the "
                                              "Record ID. Hover a cell for the full header line.")})
        if lab_col is not None:
            t.display_only_columns = frozenset({lab_col})
        if len(recs) > self._MAX_RECORDS:
            t.set_complete_rows(rows, self._MAX_RECORDS)
            t.export_notes = [f"Screen shows the first {self._MAX_RECORDS} of {len(recs)} records under the current "
                              "column sort and direction. The export contains all records in that same order; "
                              "the displayed rows are its exact prefix."]
        else:
            t.set_rows(rows, tips=tips)
            t.export_notes = ["The export contains the complete record summary in the table's current sort order."]
        t.export_notes += self._SUMMARY_NOTES
        id_map = [("Export ID map (local)…", lambda: self._export_id_map(recs, labels))] if lab_col is not None else []
        t.set_row_menu(lambda r: [("Show this record below", lambda i=r: self._show_record(i, keep_table=True))] + id_map)
        t.select_orig(self.state.get("record_index", 0))
        t.current_row_changed.connect(self._on_record_row)
        tail = ("The summary keeps composition, classification, completeness and every processing-stage outcome; the "
                "complete versioned result is the JSON export in Run provenance.")
        head = QLabel((f"<b>{len(recs)} records analysed.</b> Click a row, or move to it with the arrow keys, to show that "
                       "record below — primer design, in-silico PCR, splice alignment and family naming all use the record "
                       f"shown. {tail}" if len(recs) > 1 else f"<b>1 record analysed.</b> {tail}")
                      + (f" Showing the first {self._MAX_RECORDS} of {len(recs)} under the current sort; the export "
                         "contains all records in that order, with these displayed rows as its exact prefix."
                         if len(recs) > self._MAX_RECORDS else "")
                      + (" The Your label column shows the first word of the FASTA header you supplied, on this screen only: every "
                         "table, FASTA, GFF3/BED, figure and file uses the Record ID instead (hover the Record column for why)."
                         if lab_col is not None else ""))
        head.setTextFormat(Qt.RichText); head.setWordWrap(True); head.setObjectName("cardnote")
        self._recTable, self._recHead = t, head
        card.bodylay.addWidget(head); card.bodylay.addWidget(t)
        frow = QHBoxLayout(); frow.addStretch(1)
        record_base = ("TEagle_multiple_records" if len(recs) > 1 else
                       f"TEagle_{_record_export_stem(recs[0], 'record')}")
        frow.addWidget(_export_fasta_btn(
            "Export full records (FASTA)",
            lambda rs=recs: _full_records_fasta(rs),
            record_base, self))
        frow.addWidget(_export_table_btn(t, record_base + "_summary", self))
        card.bodylay.addLayout(frow)

    def _local_labels(self, recs):
        """Per record, the FASTA header the user supplied, or None. A header is returned only when it reproduces that record's
        public ID and digest, so a label can never sit on the wrong row; fetched records, bare pastes and any mismatch carry
        none. These strings are for the screen only: nothing here is written to a file or a manifest."""
        from teagle_core.gff3 import public_record_identity
        heads = self.state.get("record_headers") or []
        if len(heads) != len(recs):
            return [None] * len(recs)
        out = []
        for h, r in zip(heads, recs):
            ident = public_record_identity(h)
            out.append(h if (r.get("identitySource") == "local_input" and ident.get("idSha256")
                             and ident == {"id": r.get("id"), "idSha256": r.get("idSha256")}) else None)
        return out

    def _export_id_map(self, recs, labels):
        """Opt-in local key file: each neutral record ID beside the header the user supplied. It is the one export that
        carries header text, so it asks first (Cancel is the default), is CSV only, and sits outside every manifest."""
        box = QMessageBox(self)
        box.setIcon(QMessageBox.Warning); box.setWindowTitle("Export a local ID map?")
        box.setText("This file pairs each record ID with the original sequence names from your FASTA headers (sample names, file "
                    "paths or patient codes included). Anyone holding it can link TEagle's exports back to them.")
        box.setInformativeText("Do not share it if those names are sensitive. Keep it with your private data and never "
                               "send it with the exported results.")
        go = box.addButton("Export ID map…", QMessageBox.AcceptRole)
        cancel = box.addButton("Cancel", QMessageBox.RejectRole)
        box.setDefaultButton(cancel); box.setEscapeButton(cancel)
        box.exec()
        if box.clickedButton() is not go:
            return
        rows = [[r.get("id"), r.get("idSha256"), h, (r.get("composition") or {}).get("length", len(r.get("seq", ""))),
                 sequtil.sha256(r.get("seq", ""))] for r, h in zip(recs, labels) if h]
        widgets.export_table(
            ["Record ID", "Record ID SHA-256", "FASTA header (local)", "Length (bp)", "Sequence SHA-256"], rows,
            "TEagle_ID_map_LOCAL", self, fmt="csv",
            notes=["LOCAL KEY FILE — links the neutral record IDs in TEagle exports to the FASTA headers you supplied. Keep it "
                   "with your private data; do not deposit or share it with exported results.",
                   "Record ID SHA-256 is the SHA-256 of the exact header line in column 3; recomputing it verifies each pairing. "
                   "A header beginning with = + - or @ is written with a leading apostrophe by the spreadsheet-injection guard; "
                   "remove it before recomputing."])
