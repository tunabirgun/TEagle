"""Structure and classification card, feature menus and sequence helpers."""
from __future__ import annotations

import html
import os
import re
from PySide6.QtCore import Qt
from PySide6.QtGui import QCursor
from PySide6.QtWidgets import (
    QApplication, QComboBox, QDialog, QFileDialog, QFrame, QHBoxLayout, QLabel, QPushButton, QSizePolicy, QSpinBox,
    QToolTip, QVBoxLayout, QWidget)
import figures
from figures import svg_genome
import widgets
import coords
from coords import CoordGenomePanel
import theme as theme_mod
from teagle_core import (
    classify, sequtil, domains as domains_mod, structural as structural_mod, primers as primers_mod)
from teagle_core.fetch import complete_gene_model
from ui_common import (
    DOMAIN_COLS, EXPORT_HEADERS, PANEL, _CONFIDENCE_TIP, _COORD_SYSTEM, _ERV_TIP, _RANK_TIP, _SIG, _call_view,
    _caution_items, _clear_layout, _coding_view, _domains_fasta, _explanation_view, _export_fasta_btn,
    _export_table_btn, _external, _invalid_record_text, _kb_links, _note, _note_view, _orfs_fasta, _published_note,
    _rank_line, _rank_view, _record_export_stem, _record_public_id, _record_public_identity, _repeat_word, _sl,
    _te_headline)


class StructureMixin:
    def _next_step_row(self, items):
        """'Next step' followed by link buttons (text, callback). Result controls: disabled with the result they sit under."""
        row = QHBoxLayout(); row.setSpacing(theme_mod.sp(12))
        row.addWidget(_sl("Next step"))
        for text, fn in items:
            b = QPushButton(text); b.setProperty("link", True); b.setProperty("resultExport", True)
            b.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Fixed)
            b.clicked.connect(lambda _=False, f=fn: f()); row.addWidget(b, 1)
        return row

    def _go_design(self):
        self._scroll_to(self.card_primer)
        self.designBtn.setFocus()

    def _flank_next_step(self, rec):
        """(link or None, note or None): outside primers need flanking DNA; offer the flank picker only when the record has
        at least the current minimum product length of it on a side, otherwise say what to fetch."""
        spans = [(e["element_span"], e["type"]) for e in rec.get("structural") or []
                 if e.get("element_span") and not e.get("advisory") and str(e["type"]).startswith(("LTR", "TIR"))]
        if not spans:
            return None, None
        s, t = max(spans, key=lambda x: x[0][1] - x[0][0])
        word = _repeat_word(rec) or t.split(" ")[0]
        n = (rec.get("composition") or {}).get("length") or len(rec.get("seq") or "")
        try:
            need = int(float(self.pfields["pMin"].text()))
        except (ValueError, KeyError):
            need = self.PRESETS["standard"]["pMin"]
        left, right = s[0], n - s[1]
        if max(left, right) >= need:
            return (f"Flanking sequence around the {word} element (for primers outside it) →",
                    lambda: self._flank_picker(s[0], s[1], "+", f"{word}_element")), None
        return None, (f"Primers outside the element need flanking DNA; this record has {left:,} bp before and {right:,} bp "
                      f"after the {word} element ({self._span(*s)} of {n:,} bp). To genotype the insertion, fetch the locus with "
                      f"flanking sequence (panel {PANEL['specimen']} · Fetch by coordinate) and analyse that.")

    def _row_action_bar(self, table, describe, primary=None):
        """Visible controls for the current row of `table`, to sit at the left of its button row: which row is current,
        its main action, and 'Row actions' — exactly the right-click menu (same builder, so no drift). `describe(orig)`
        names the row; `primary` = (button text, menu-label prefix) binds the button to the first menu item with that prefix
        and disables it where the menu offers none. Returns the widgets to add to the row."""
        label = QLabel("Select a row:"); label.setObjectName("cardmeta")
        label.setMinimumWidth(1); label.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)
        acts = QPushButton("Row actions ▾"); acts.setProperty("sm", True); acts.setProperty("resultExport", True)
        acts.setAccessibleName("Actions for the selected row")
        acts.setToolTip("The same menu as right-clicking the row, or pressing Shift+F10 / the Menu key on it.")
        main = None
        if primary:
            main = QPushButton(primary[0]); main.setProperty("sm", True); main.setProperty("resultExport", True)
        hint = QLabel("Also: right-click a row, or Shift+F10."); hint.setObjectName("cardmeta")
        hint.setMinimumWidth(1); hint.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)
        fn = {}

        def sync(_orig=None):
            orig = table.current_orig()
            label.setText(f"{describe(orig)}:" if orig >= 0 else "Select a row:")
            live = orig >= 0 and table.isEnabled()
            acts.setEnabled(live)
            if main is not None:
                fn["go"] = next((f for l, f in table.row_items(orig) if l.startswith(primary[1])), None) if live else None
                main.setEnabled(fn["go"] is not None)
        acts.clicked.connect(lambda _=False: table.show_row_menu(acts.mapToGlobal(acts.rect().bottomLeft())))
        if main is not None:
            main.clicked.connect(lambda _=False: fn.get("go") and fn["go"]())
        table.current_row_changed.connect(sync)
        if table.rowCount():
            table.select_orig(table._orig(0))
        sync()
        return [w for w in (label, main, acts, hint) if w is not None]

    def _render_struct_card(self, rec, res, keep_table=False):
        """Record table (persistent across selections) above a detail block that is rebuilt per record. Selecting a
        row rebuilds only the detail: the table that emitted the signal keeps its sort order and scroll position."""
        card = self.card_struct
        if keep_table and getattr(self, "_recDetail", None) is not None:
            _clear_layout(self._recDetail.layout())
        else:
            card.clear_body()
            card.expand()
            self._render_record_table(card)
            self._recDetail = QWidget()
            dl = QVBoxLayout(self._recDetail); dl.setContentsMargins(0, 0, 0, 0); dl.setSpacing(card.bodylay.spacing())
            card.bodylay.addWidget(self._recDetail)
        body = self._recDetail.layout()
        recs = self.state.get("records") or []
        k = (self.state.get("record_index") or 0) + 1
        who = f"{_record_public_id(rec)} · {(rec.get('composition') or {}).get('length', len(rec.get('seq') or '')):,} bp"
        mine = self._local_labels(recs)
        mine = mine[k - 1] if 0 < k <= len(mine) else None
        if mine:                                              # the screen shows the user's own first word; no file ever does
            who += f" · your label “{mine.split()[0]}”"
        title = QLabel(f"Record {k} of {len(recs)} · {who}" if len(recs) > 1 else who)
        title.setObjectName("sech"); title.setTextInteractionFlags(Qt.TextSelectableByMouse)
        body.addWidget(title)
        from teagle_core import gff3
        cl = rec.get("classification") or {}
        banner = QFrame(); banner.setObjectName("classbn")
        bl = QVBoxLayout(banner)
        bl.setContentsMargins(theme_mod.sp(16), theme_mod.sp(14), theme_mod.sp(16), theme_mod.sp(14)); bl.setSpacing(theme_mod.sp(9))
        if not rec.get("valid", True) or not cl:
            # nothing was searched: say so, and state no domain, structural or class result at all
            top = QHBoxLayout(); top.setSpacing(theme_mod.sp(10))
            failed = rec.get("valid", True) and _call_view(rec)["status"] == "failed"
            big = QLabel("Classification failed" if failed else "Not analysed"); big.setObjectName("classbig"); top.addWidget(big)
            held = QLabel("FAILED" if failed else "NOT ANALYSED"); held.setObjectName("cf"); held.setProperty("level", "Candidate")
            top.addWidget(held); top.addStretch(1); bl.addLayout(top)
            why = (_invalid_record_text(rec, self.state.get("record_index", 0) + 1) if not rec.get("valid", True)
                   else "The classification step raised an error for this record, so no call was made. The analysis "
                        "notes below name it." if failed else "The classification step did not run for this record.")
            wl = QLabel(why); wl.setObjectName("classexp"); wl.setWordWrap(True); bl.addWidget(wl)
            nl = QLabel("No structural, reading-frame or protein-domain search was run on this record, so nothing "
                        "here says whether it is a transposable element. Fix the sequence, or select another record.")
            nl.setObjectName("cardnote"); nl.setWordWrap(True); bl.addWidget(nl)
            body.addWidget(banner)
            return
        tdr = gff3.is_neutral_tdr(cl)
        base = f"TEagle_{_record_export_stem(rec, 'record')}"      # every table export names its record
        domain_status, unsearched = self._scan_state(rec, cl)
        te = cl.get("te_class", "—")
        # TIER 1 — CALL: te_class + confidence chip + superfamily
        top = QHBoxLayout(); top.setSpacing(theme_mod.sp(10))
        order_withheld = cl.get("outcome") == "withheld"
        unassessed = bool(cl.get("domains_unavailable")) and te in ("none", "unclassified")
        headline = ("Not assessed — domain scan unavailable" if unassessed else
                    "Terminal direct-repeat candidate" if order_withheld and te == "repeat/structural-only"
                    else _te_headline(te))
        big = QLabel(headline); big.setObjectName("classbig")
        top.addWidget(big)
        if order_withheld or unassessed:
            # name what was withheld: a lone RT has neither class nor order, an unavailable scan has no call at all
            chip = ("NOT ASSESSED" if unassessed else
                    "CLASS & ORDER WITHHELD" if te in ("unclassified", "none") else "ORDER WITHHELD")
            held = QLabel(chip); held.setObjectName("cf"); held.setProperty("level", "Candidate")
            top.addWidget(held)
        call = _call_view(rec)
        cautions = _caution_items(rec)
        conf = call["confidence"]
        if conf:                                              # a chip beside "no signature" reads as a candidate call, so none is drawn
            cfl = QLabel(f"CALL CONFIDENCE: {conf.upper()}" + (" · SEE CAUTION" if cautions else ""))
            cfl.setObjectName("cf"); cfl.setProperty("level", conf)
            cfl.setAccessibleName(f"Call confidence: {conf}")
            cfl.setToolTip(("Scope caution below — the confidence value is unchanged and does not account for the "
                            "conditions listed.\n" if cautions else "") + _CONFIDENCE_TIP)
            top.addWidget(cfl)
        top.addStretch(1)
        bl.addLayout(top)
        kls_text = ("no call made — the protein-domain scan did not run"
                    if unassessed and te == "none" else cl.get("superfamily", "—"))
        kls = QLabel(kls_text + self._src_html("Wicker2007"))
        kls.setObjectName("classkls"); kls.setWordWrap(True); kls.setTextFormat(Qt.RichText)
        _external(kls); _kb_links(kls); bl.addWidget(kls)
        rank = _rank_view(cl)
        if rank:
            rk = QLabel(_rank_line(rank)); rk.setObjectName("classkls"); rk.setWordWrap(True)
            rk.setToolTip(_RANK_TIP + (_ERV_TIP if rank["erv"] else ""))
            bl.addWidget(rk)
        published = _published_note(rec)
        if published:                                         # a bundled example is a positive control: show the published class beside the call
            bl.addWidget(_note(published[0], published[1]))
        if cautions:
            strip = _note("<b>Scope caution</b>" + "".join(f"<br>• {html.escape(c)}" for c in cautions), "warn")
            strip.setTextFormat(Qt.RichText); strip.setAccessibleName("Scope caution")
            bl.addWidget(strip)
        # TIER 2 — WHY: the reasoning, with reading line-height. The display wording (leading "Classified as …"
        # trimmed, repeats collapsed, a capital after the domain list) is built from the record by _explanation_view;
        # cl["explanation"] itself stays byte-stable for the exports.
        if cl.get("explanation"):
            ex = QLabel(f"<div style='line-height:150%'>{html.escape(_explanation_view(rec), quote=False)}</div>")
            ex.setObjectName("classexp"); ex.setWordWrap(True); ex.setTextFormat(Qt.RichText)
            bl.addWidget(ex)
        # TIER 3 — SCOPE: reliability & completeness. An honest limitation → stays default-VISIBLE, grouped + labelled.
        comp = cl.get("completeness")                         # scoped domain completeness (Axis 2 of reliability)
        if comp:
            rl = QLabel("RELIABILITY & COMPLETENESS"); rl.setObjectName("sectionlabel"); bl.addWidget(rl)
            _cod = comp.get("coding") or {}
            # A record flagged truncated must not also read "intact". The signal is the classifier's own
            # (a domain-bearing ORF reaches the record end), not the wording of an engine note: the note also
            # fires for domain-free ORFs and says "3′" for a minus-strand ORF open at the 5′ end.
            tier = comp["tier"]
            if _cod.get("truncated_frame") is True and "intact" in tier.lower():
                tier += (" at the domain level — but a domain-bearing ORF runs off the end of the supplied sequence, "
                         "so the element may be truncated (see note below)")
            _dead = cl.get("domains_unavailable")
            _, arch = self._orders(cl, domain_status)
            line = (f"<b>Structural completeness:</b> {tier}  ·  {comp.get('kind','')}"
                    f"<br><b>Domain architecture:</b> {arch}")
            modules = self._module_lines(rec, cl, comp)
            if modules:
                line += "<br>" + modules
            # The coding axis, stated beside the tier. "partial" is a domain-ledger verdict, and on its own a
            # reader takes it as "this copy is decayed" — which is why LTR_retriever and TEsorter can say
            # "intact"/"Complete" about the same element. Saying whether the detected domains actually sit in
            # one uninterrupted reading frame separates "modules not detected" from "copy broken".
            if _cod.get("note"):
                line += f"<br><b>Coding structure:</b> {html.escape(_coding_view(rec), quote=False)}"
            cw = QLabel(f"<div style='line-height:150%'>{line}</div>"); cw.setObjectName("classexp"); cw.setTextFormat(Qt.RichText); cw.setWordWrap(True)
            bl.addWidget(cw)
            # the essential caveat stays visible; the full methodology collapses behind a disclosure (less crowding)
            caveat = QLabel("Structural evidence in this one sequence, scored against the tested Pfam panel — "
                            "not a claim about expression, transposition competence, or any individual genome.")
            caveat.setObjectName("cardnote"); caveat.setWordWrap(True); bl.addWidget(caveat)
            scope = QLabel(f"Domains tested: {comp.get('scope','')}. "
                           "(Completeness after Wicker 2007 / TEsorter / LTR_retriever.) The tier reports how much of "
                           "the expected architecture is present at the domain level in THIS SEQUENCE. It says "
                           "nothing about three separate questions it is easily confused with: whether the element "
                           "is transcribed, whether it retains transposition or infection competence, and whether "
                           "any given individual or population carries this insertion. Answering those needs "
                           "expression data, a functional assay, and population genotyping respectively "
                           "(Lanciano &amp; Cristofari 2020, Nat Rev Genet 21:721-736).")
            scope.setObjectName("cardnote"); scope.setWordWrap(True)
            self._disclosure(bl, "Scope and methods", scope, expanded=False)
        else:
            # a no-evidence / structural-only call still needs its SCOPE stated — the essential caveat stays visible
            rl = QLabel("RELIABILITY & COMPLETENESS"); rl.setObjectName("sectionlabel"); bl.addWidget(rl)
            # Same not-detected-vs-not-tested distinction as the completeness branch above: when the scan
            # raised, "no panel domain was detected" is a false statement about the sequence.
            if cl.get("domains_unavailable"):
                caveat = QLabel("The protein-domain scan did not run for this record (see note below), so no "
                                "domain was tested — this is not a finding that none is present.")
                scope_txt = (f"Domains that WOULD have been tested: {classify.DOMAINS_TESTED}. None of them were "
                             "searched on this record, so nothing here says whether the sequence is a "
                             "transposable element.")
            elif domain_status == "not_run":
                caveat = QLabel("No protein-domain search was possible: the sequence has no open reading frame of 40 "
                                "amino acids or more to search. This is not a finding that no domain is present.")
                scope_txt = (f"Domains that WOULD have been tested: {classify.DOMAINS_TESTED}. None of them were "
                             "searched on this record.")
            else:
                caveat = QLabel("No panel domain was detected — scoped to the tested Pfam panel, not proof the sequence is not a TE.")
                scope_txt = (f"Domains tested: {classify.DOMAINS_TESTED}. None were detected — this is relative to this "
                             "bundled Pfam profile panel; a divergent or unmodelled domain reads as not-detected, not as "
                             "proof the sequence is not a transposable element.")
            caveat.setObjectName("cardnote"); caveat.setWordWrap(True); bl.addWidget(caveat)
            scope = QLabel(scope_txt)
            scope.setObjectName("cardnote"); scope.setWordWrap(True)
            self._disclosure(bl, "Scope and methods", scope, expanded=False)
        body.addWidget(banner)

        # genome viewer
        model = figures.gv_tracks_from_rec(rec)
        if model["tracks"]:
            gv = CoordGenomePanel(svg_genome, f"TEagle_{_record_export_stem(rec, 'locus')}_genome")   # the file names its record
            gv.apply_app_theme(self.theme)                # open in the current app theme
            gv.coord_base = self.coord_base
            gv.set_model(model)
            gv.set_feature_menu(lambda region, s=rec.get("seq") or None: self._region_menu(region, s))
            gv.setMinimumHeight(round(260 * theme_mod.UI_SCALE))
            body.addWidget(gv)

        # retroviral transcript architecture (ERV) — the correct coding-organisation model + cis-element legend
        arch = rec.get("retroviral")
        if arch:
            # swatches derive from the figure palette (theme.ARCHCOL / theme.CISCOL) instead of repeating the
            # hex here — a legend that duplicated the hues could drift out of step with the bands it labels.
            _cis = [(k, lbl) for k, pre, lbl in (("PBS", "PBS", "PBS"), ("PPT", "PPT", "PPT"),
                                                 ("PAS", "polyA-signal", "PAS · polyA signal (motif)"))
                    if any(e["type"].startswith(pre) for e in rec.get("structural", []))]
            leg = (f"<span style='color:{theme_mod.ARCHCOL['exon']}'>■</span> env exon &nbsp; "
                   f"<span style='color:{theme_mod.ARCHCOL['intron']}'>■</span> gag–pro–pol interval (env intron)"
                   + "".join(f" &nbsp; <span style='color:{theme_mod.CISCOL[k]}'>■</span> {lbl}" for k, lbl in _cis))
            legw = QLabel(leg); legw.setObjectName("orient"); legw.setTextFormat(Qt.RichText); legw.setWordWrap(True)
            body.addWidget(legw)
            note = QLabel("<b>Endogenous retrovirus — transcript architecture.</b> " + arch["note"] +
                          " For the exact splice bases, send a real env transcript to the Splice structure panel (03b). " +
                          arch.get("subsplice_note", ""))
            note.setObjectName("orient"); note.setWordWrap(True); note.setTextFormat(Qt.RichText)
            body.addWidget(note)

        # structural table (right-click a row → copy FASTA/DNA/coords, design primer here)
        struct = rec.get("structural", [])
        if struct:
            hdr = QHBoxLayout(); hdr.addWidget(_sl("Structural evidence")); hdr.addStretch(1)
            hdr.addWidget(self._src_link("Wicker2007")); body.addLayout(hdr)
            cap = QLabel(self._coord_caption() + " A paired feature lists both copies, upstream · downstream.")
            cap.setObjectName("cardmeta"); cap.setWordWrap(True); body.addWidget(cap)
            t = self._coord_table(["Feature", self._ch("Coords"), "Len", "Metric", "Method"], self._gloss(),
                                  [self._struct_row(e, tdr=tdr, cl=cl, base=0) for e in struct])
            t.export_header_map = {**EXPORT_HEADERS, "Coords (1-based)": "Coords (0-based)"}
            missing_notes = self._searched_not_found(rec, cl, tdr)
            # the caveats shown beneath this table are part of what the numbers mean, so they go in the file
            t.export_notes = [
                self._record_note(rec, res),
                "Structural evidence detected de novo by heuristic terminal-repeat detectors — not retrieved "
                "from a database. Coordinates are 0-based half-open [start, end). A minus-strand row (−) is "
                "still given in the record's own coordinates.",
                "The polyA-signal row is a sequence motif with its downstream element — advisory context only. "
                "It does not locate the U3-R-U5 boundaries, the cleavage site, or the transcript end, which "
                "need RNA evidence TEagle does not use.",
                f"An LTR pair is accepted at >= {structural_mod.MIN_LTR_IDENTITY:g}% identity. This is a "
                "detection threshold, not a calibrated age boundary: divergence can make older copies more "
                "likely to fall below it, but identity also depends on lineage, gene conversion, selection, "
                "sequence quality and detector sensitivity. Sub-threshold pairs are reported as advisory candidates.",
            ]
            if any(e.get("length_is_lower_bound") for e in struct):
                t.export_notes.append(
                    "A repeat marked 'length is a lower bound' still matched where the record ends or where the scan "
                    "stops, so its Len is a minimum and the real repeat may be longer. Len itself is the measured number.")
            if tdr:
                t.export_notes.append("TDR = terminal direct repeat whose LTR origin is unconfirmed (order withheld).")
            t.export_notes += missing_notes
            t.set_rows([self._struct_row(e, tdr=tdr, cl=cl) for e in struct],
                       tips=[self._struct_tips(e, cl=cl) for e in struct])
            t.set_row_menu(lambda r: self._struct_menu(struct[r]))
            t.fit_rows()
            body.addWidget(t)
            methods, seen = [], set()                         # the Method column is elided on a narrow card: the full text, wrapped, once per kind
            for s_ in struct:
                shown, _c, _l, _m, how = self._struct_row(s_, tdr=tdr, cl=cl)
                if (shown, how) in seen or not how:
                    continue
                seen.add((shown, how))
                hedge = self._struct_tips(s_, cl=cl)[0]
                methods.append(f"{html.escape(shown)} — {html.escape(how)}" + (f" <i>{html.escape(hedge.splitlines()[0])}</i>" if hedge else ""))
            if methods:
                mn = QLabel("<b>Detection methods</b><br>" + "<br>".join(methods)); mn.setObjectName("cardnote")
                mn.setTextFormat(Qt.RichText); mn.setWordWrap(True); body.addWidget(mn)
            if missing_notes:                                 # what was searched and not found, stated beside the table
                nf = QLabel("<b>Searched, not found</b><br>" + "<br>".join(html.escape(n) for n in missing_notes))
                nf.setObjectName("orient"); nf.setTextFormat(Qt.RichText); nf.setWordWrap(True)
                body.addWidget(nf)
            # The polyA-signal row is a MOTIF, not a located cleavage site. Its caveat rides in the row's
            # tooltip, but a tooltip is not a disclosure a user is guaranteed to see — so when the row is
            # present the limit is stated in the card itself, next to the table that shows it.
            if any(e.get("type", "").startswith("polyA-signal") for e in struct):
                pas_note = QLabel("The polyA-signal row is a sequence motif with its downstream element — advisory "
                                  "context only. It does not locate the U3–R–U5 boundaries, the cleavage site, or "
                                  "the transcript end, which need RNA evidence TEagle does not use.")
                pas_note.setObjectName("orient"); pas_note.setWordWrap(True)
                body.addWidget(pas_note)
            srow = QHBoxLayout()                              # exportable like every sibling results table
            def _struct_name(o, s=struct):
                r = self._struct_row(s[o], tdr=tdr, cl=cl)
                return f"{r[0].split(' ')[0]} {r[1]}"
            for w in self._row_action_bar(t, _struct_name):
                srow.addWidget(w)
            srow.addStretch(1)
            ib = coords.insertion_bounds(rec)
            if ib and max(ib["up"], ib["down"]) >= primers_mod.MIN_PRIMER_SIZE:
                geno = QPushButton("Genotype this insertion…"); geno.setProperty("sm", True); geno.setProperty("resultExport", True)
                geno.setToolTip("Open primer design (panel 04) with the genotyping options for this insertion.")
                geno.clicked.connect(lambda _=False: self._goto_assay("genotype"))
                srow.addWidget(geno)
            srow.addWidget(_export_table_btn(t, f"{base}_structural", self))
            body.addLayout(srow)

        # ORFs
        orfs = rec.get("orfs", [])
        if orfs:
            body.addWidget(_sl(f"ORFs (≥40 aa) — {len(orfs):,}"))
            cap = QLabel("Start and End are 1-based and inclusive: the first base of the record is 1 and End is the last base of "
                         "the ORF, so length = End − Start + 1. Exports keep the 0-based half-open numbers." if self.coord_base else
                         "Start and End are 0-based and half-open: the first base of the record is 0 and End is not included, "
                         "so length = End − Start. (GenBank, SnapGene and IGV count from 1 and include the end; "
                         "switch with the Coords button in the header.)")
            cap.setObjectName("cardmeta"); cap.setWordWrap(True); body.addWidget(cap)
            ends = "End (inclusive)" if self.coord_base else "End (exclusive)"
            orf_row = lambda i, o, b: [i + 1, o["strand"], o["frame"], o["start"] + b, o["end"], o["length_aa"], _orf_domains(i)]
            t = self._coord_table(["#", "Strand", "Frame", self._ch("Start"), ends, "aa", "Domains"],
                                  self._gloss({**self._RESULT_GLOSS,
                                               "#": "ORF number: its position in the engine's list (longest first). The domain table's "
                                                    "ORF column uses the same number."}))
            t.export_header_map = {**EXPORT_HEADERS, "Start (0-based)": "Start (0-based)", "Start (1-based)": "Start (0-based)",
                                   "End (exclusive)": "End (0-based, exclusive)", "End (inclusive)": "End (0-based, exclusive)"}
            limit = domains_mod.MAX_ORFS_SCANNED
            t.export_notes = [
                self._record_note(rec, res),
                "Open reading frames >= 40 aa, 0-based half-open coordinates on the record's "
                "own sequence. An ORF is a reading frame, not evidence that it is translated.",
                f"Domains: codes detected in that ORF; only the {limit} longest ORFs are searched ('not searched' "
                "marks the rest, 'not assessed' means the domain search did not run). '—' means searched, none in the panel."]
            scan = (rec.get("processing") or {}).get("domains") or {}
            in_orf = figures._domain_orf_map(rec)
            searched = scan.get("orfsScanned") if scan.get("status") in ("complete", "partial") else 0

            def _orf_domains(i):
                if cl.get("domains_unavailable") or searched is None or not scan.get("status") or scan["status"] in ("failed", "not_run"):
                    return "not assessed"
                return " · ".join(in_orf[i]) if i in in_orf else ("—" if i < searched else "not searched")
            rows = [orf_row(i, o, self.coord_base) for i, o in enumerate(orfs)]
            if self.coord_base:                               # a file keeps the canonical numbers, in the table's sort order
                t.canon_rows = [orf_row(i, o, 0) for i, o in enumerate(orfs)]
            if len(rows) > limit:
                t.set_complete_rows(rows, limit)
            else:
                t.set_rows(rows)
            t.set_row_menu(lambda r: self._feat_menu(orfs[r]["start"], orfs[r]["end"], orfs[r]["strand"],
                                                     f"ORF_{orfs[r]['strand']}{orfs[r]['frame']}",
                                                     src_seq=rec.get("seq")))
            t.fit_rows()
            body.addWidget(t)
            orow = QHBoxLayout()
            for w in self._row_action_bar(
                    t, lambda o, os_=orfs: f"ORF {os_[o]['strand']}{os_[o]['frame']} {self._span(os_[o]['start'], os_[o]['end'])}",
                    ("Design primers here", "Design primer here")):
                orow.addWidget(w)
            if len(rows) > limit:
                more = QPushButton(f"Show all {len(rows)} ORFs"); more.setProperty("link", True)
                more.setToolTip(f"Rows follow the current sort. Only the {limit} longest ORFs are searched for protein domains.")
                more.setAccessibleName(f"Show all {len(rows)} ORFs")

                def _toggle_orfs(_=False, b=more, tt=t, n=len(rows)):
                    full = tt._complete_limit < n
                    tt.set_complete_rows(rows, n if full else limit); tt.fit_rows()
                    b.setText(f"Show only the first {limit} ORFs" if full else f"Show all {n} ORFs")
                    b.setAccessibleName(b.text())
                more.clicked.connect(_toggle_orfs)
                orow.addWidget(more)
            orow.addStretch(1)
            orow.addWidget(_export_fasta_btn("Export ORFs (DNA FASTA)", lambda r=rec: _orfs_fasta(r),
                                             f"{base}_ORFs_DNA", self))
            orow.addWidget(_export_table_btn(t, f"{base}_ORFs", self))
            body.addLayout(orow)

        # domains (right-click → copy protein/DNA/FASTA/coords, design primer here)
        doms = rec.get("domains", [])
        if doms:
            hdr = QHBoxLayout(); hdr.addWidget(_sl("Protein domains (HMMER)")); hdr.addStretch(1)
            hdr.addWidget(self._src_link("Pfam")); body.addLayout(hdr)
            dom_row = lambda d, b: [d["domain"], d.get("label", ""), d.get("pfam", ""), d.get("strand") or "+",
                                    self._orf_number(d, orfs), f"{d['aa'][0]}–{d['aa'][1]}", self._span(*d["nt"], b),
                                    d.get("score"), f"{d.get('evalue'):.1e}" if d.get("evalue") is not None else "",
                                    d.get("confidence", "")]
            aa_gloss = self._RESULT_GLOSS["aa in ORF"].replace("0-based half-open", "1-based inclusive" if self.coord_base else "0-based half-open")
            t = self._coord_table([self._ch("nt") if h == "nt" else h for h in DOMAIN_COLS], self._gloss({**self._RESULT_GLOSS, "aa in ORF": aa_gloss}))
            t.export_header_map = {**EXPORT_HEADERS, "aa in ORF": "aa in ORF (1-based, inclusive)",
                                   "nt (0-based)": "nt (0-based, half-open)", "nt (1-based)": "nt (0-based, half-open)",
                                   "ORF": "ORF # (as numbered in the ORF table)"}
            t.export_notes = [
                self._record_note(rec, res),
                f"Protein domains from the bundled {len(domains_mod.DOMAIN_INFO)}-model Pfam panel via "
                "HMMER/pyhmmer. Hit confidence is the per-domain confidence from the HMMER i-Evalue — the "
                "reliability of THIS domain hit, separate from the element's Call confidence and from domain completeness.",
                f"Only the {domains_mod.MAX_ORFS_SCANNED} longest ORFs are searched. A domain absent here was "
                "not necessarily absent from the whole sequence.",
                f"Domains tested: {classify.DOMAINS_TESTED}. 'Not detected' is relative to this panel — a "
                "divergent or unmodelled domain reads as not-detected, not as element decay.",
            ]
            t.set_rows([dom_row(d, self.coord_base) for d in doms])
            if self.coord_base:
                t.canon_rows = [dom_row(d, 0) for d in doms]
            t.set_row_menu(lambda r: self._feat_menu(doms[r]["nt"][0], doms[r]["nt"][1], doms[r].get("strand", "+"),
                                                     doms[r]["domain"], protein=doms[r].get("protein"),
                                                     src_seq=rec.get("seq")))
            t.fit_rows()
            body.addWidget(t)
            dhint = QLabel("The last column is <b>Hit confidence</b> (per-domain HMMER confidence — not the element's Call "
                           "confidence). <b>aa in ORF</b> counts residues of the "
                           f"ORF's protein from 1, both ends included; <b>nt</b> is {'1-based, inclusive' if self.coord_base else '0-based half-open'} on the record. On a "
                           "narrow window, scroll the table sideways — or hide the input panel (Ctrl+B) — to reach "
                           "every column.")
            dhint.setObjectName("orient"); dhint.setTextFormat(Qt.RichText); dhint.setWordWrap(True)
            body.addWidget(dhint)
            drow = QHBoxLayout()
            for w in self._row_action_bar(t, lambda o, ds=doms: f"{ds[o]['domain']} {self._span(*ds[o]['nt'])}",
                                          ("Design primers here", "Design primer here")):
                drow.addWidget(w)
            drow.addStretch(1)
            if _domains_fasta(rec):
                drow.addWidget(_export_fasta_btn("Export domain proteins (FASTA)", lambda r=rec: _domains_fasta(r),
                                                 f"{base}_domain_proteins", self))
            drow.addWidget(_export_table_btn(t, f"{base}_domains", self))
            body.addLayout(drow)

        # Self-similarity: the one view that looks for repeats WITHOUT being told what to expect, so it
        # can catch what a targeted detector missed. Generated on request rather than with every
        # analysis — it is the only panel whose cost grows with the square of the repeat content.
        srow = QHBoxLayout(); srow.addStretch(1)
        selfsim = QPushButton("Self-similarity plot…")
        selfsim.setProperty("sm", True)
        selfsim.setToolTip(
            "Open a resizable window comparing this sequence with itself, plotting every exact word match.\n\n"
            "Finds repeats the targeted detectors were not looking for: a block where nothing was "
            "called is worth investigating. Direct repeats (LTRs) appear off the main diagonal; "
            "inverted repeats (TIRs) appear on the anti-diagonal.\n\n"
            "It cannot do the reverse — exact word matching is not an alignment, so a diverged repeat "
            "fades and an absent block is not evidence that no repeat exists.")
        selfsim.clicked.connect(lambda _=False, r=rec: self._build_selfsim(r))
        srow.addWidget(selfsim)
        body.addLayout(srow)

        # Annotation export: the whole record as GFF3/BED so the call survives into a genome browser.
        # Placed on the structure card because that is where the annotation being exported is shown.
        arow = QHBoxLayout(); arow.addStretch(1)
        ann = QPushButton("Export annotation ▾"); ann.setProperty("sm", True)
        ann.setToolTip("GFF3 (Sequence Ontology terms, sub-features, embedded FASTA) or BED6 intervals")
        ann.clicked.connect(lambda: self._export_annotation(rec))
        arow.addWidget(ann)
        body.addLayout(arow)

        # gene model (exon/intron/CDS) — only when a fetched accession carries feature annotation.
        # For an ERV the host-style CDS/exon view is the misleading one (it shows the env CDS as a single
        # "exon"); the retroviral transcript architecture above is the correct model, so the raw gene model is
        # DE-EMPHASISED behind a collapsed toggle. For non-ERV TEs (a real TE-in-host-gene) it stays visible.
        gm = self.state.get("features")
        if isinstance(gm, dict) and (gm.get("exons") or gm.get("cds")):
            gm = complete_gene_model(gm)                       # fill CDS-implied exons (idempotent; covers old caches)
            derived = any(e.get("derived") for e in gm.get("exons", []))
            _tc = (cl.get("te_class") or "")
            demote = bool(cl.get("is_erv") and rec.get("retroviral"))   # ERV with a transcript-architecture model
            gmbox = QWidget(); gmlay = QVBoxLayout(gmbox); gmlay.setContentsMargins(0, 0, 0, 0); gmlay.setSpacing(6)
            title = "Gene model — NCBI annotation" + (" (feature table + CDS-inferred exons)" if derived else " (feature table)")
            gmlay.addWidget(_sl(title))
            if not demote and (cl.get("is_erv") or _tc.startswith(("LTR", "LINE", "retro", "DNA", "DIRS"))):
                cav = QLabel("This is a transposable element, not a host gene: its coding organisation is the domain "
                             "architecture above, not a host exon–intron structure. The blocks below are the record's own CDS annotation.")
                cav.setObjectName("orient"); cav.setWordWrap(True); gmlay.addWidget(cav)
            legend = ("<span style='color:#009E73'>■</span> exon · "
                      "<span style='color:#8792a0'>■</span> intron · "
                      "<span style='color:#D55E00'>■</span> CDS · "
                      "<span style='color:#5b6b7a'>■</span> flank · "
                      "<span style='color:#c3ccd6'>■</span> gap")
            if derived:
                legend += (" · <span style='color:#7fd3b8'>■</span> <b>exon*</b> = derived from the record's "
                           "CDS/mRNA, not a separate exon annotation")
            leg = QLabel(legend); leg.setTextFormat(Qt.RichText); leg.setWordWrap(True); leg.setObjectName("orient")
            gmlay.addWidget(leg)
            length = rec.get("composition", {}).get("length") or 1
            gmodel = figures.gv_tracks_from_gene(gm, length, include_flanks=True)   # flanks + gaps clickable too
            if gmodel["tracks"]:
                gvg = CoordGenomePanel(svg_genome, "TEagle_genemodel"); gvg.apply_app_theme(self.theme)
                gvg.coord_base = self.coord_base; gvg.set_model(gmodel)
                gvg.set_feature_menu(lambda region, s=rec.get("seq") or None: self._region_menu(region, s))
                gvg.setMinimumHeight(round(200 * theme_mod.UI_SCALE))
                gmlay.addWidget(gvg)
            gmrow = QHBoxLayout(); gmrow.addStretch(1)
            gmexport = QPushButton("Export gene model ▾"); gmexport.setProperty("sm", True)
            gmexport.setProperty("resultExport", True)
            gmexport.setToolTip("Export the displayed exon/intron/CDS model as GFF3, BED, or feature FASTA")
            gmexport.clicked.connect(lambda _=False, model=gm, sequence=rec.get("seq") or "",
                                     record_id=rec.get("id") or "record", source=rec.get("identitySource"),
                                     digest=rec.get("idSha256"):
                                     self._export_gene_model(model, sequence, record_id, "fetched_gene_model",
                                                             identity_source=source, record_id_sha256=digest))
            gmrow.addWidget(gmexport); gmlay.addLayout(gmrow)
            ptr = QPushButton("For exon–intron structure from a supplied transcript, use Splice structure (03b) →")   # ANNOTATED vs MEASURED
            ptr.setProperty("link", True); ptr.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Fixed)
            ptr.clicked.connect(lambda: self._scroll_to(self.card_splice))
            gmlay.addWidget(ptr)
            if demote:                                         # collapse the host-style view for an ERV
                gmbox.setVisible(False)
                gmtog = QPushButton("▸ Record's raw CDS annotation (host-style)")
                gmtog.setAccessibleName("Record's raw CDS annotation (host-style)")
                gmtog.setAccessibleDescription("Disclosure toggle — shows or hides the record's host-style CDS annotation.")
                gmtog.setToolTip("The record's host-style CDS/exon annotation — de-emphasised because the retroviral "
                                 "transcript architecture above is the correct coding model for an ERV.")
                gmtog.setProperty("link", True); gmtog.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Fixed)
                def _tgm(_=False, b=gmbox, t=gmtog):
                    v = not b.isVisible(); b.setVisible(v)
                    t.setText(("▾" if v else "▸") + t.text()[1:])
                gmtog.clicked.connect(_tgm)
                body.addWidget(gmtog)
            body.addWidget(gmbox)

        for note in rec.get("notes", []):
            n = QLabel("• " + _note_view(note, rec)); n.setObjectName("orient"); n.setWordWrap(True)
            body.addWidget(n)

        # explicit methodology — which database / consensus / parameters define each evidence layer
        mtoggle = QPushButton("▸ Methods && databases"); mtoggle.setProperty("link", True)   # '&&' -> a literal '&' (QPushButton eats a lone '&' as a mnemonic)
        mtoggle.setAccessibleName("Methods & databases")
        mtoggle.setAccessibleDescription("Disclosure toggle — shows or hides the methodology and database provenance.")
        mtoggle.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Fixed)
        mbox = QLabel(self._methods_html()); mbox.setObjectName("orient"); mbox.setWordWrap(True)
        mbox.setTextFormat(Qt.RichText); _external(mbox); _kb_links(mbox); mbox.setVisible(False)
        def _tgl():
            v = not mbox.isVisible(); mbox.setVisible(v)
            mtoggle.setText("▾ Methods && databases" if v else "▸ Methods && databases")
        mtoggle.clicked.connect(_tgl)
        body.addWidget(mtoggle); body.addWidget(mbox)
        flank_link, flank_note = self._flank_next_step(rec)
        steps = [(f"Design primers for this record → panel {PANEL['primer']}", self._go_design)]
        if flank_link:
            steps.append(flank_link)
        body.addLayout(self._next_step_row(steps))
        if flank_note:
            fn = QLabel(flank_note); fn.setObjectName("cardnote"); fn.setWordWrap(True); body.addWidget(fn)

    def _methods_html(self):
        """A plain statement of exactly what defines each evidence layer — the database/model, the consensus
        source, and the thresholds — so the annotation is never a black box."""
        p, h, w, dfam = (self._src_html(k) for k in ("Pfam", "HMMER", "Wicker2007", "Dfam"))
        return (
            "<b>Protein domains</b> — profile-HMM search (HMMER" + h + ", run in-process via pyhmmer) of the "
            f"{domains_mod.MAX_ORFS_SCANNED} longest 6-frame ORFs (≥ 40 aa; shorter or later ORFs are not searched) "
            "against a bundled Pfam-A" + p + " TE-domain profile set "
            f"({self._domain_panel_html()}). A hit is kept at i-Evalue &lt; 1e-3 (strictly below, after the "
            "whole-sequence E-value filter at the same cutoff); the gag + env models let TEagle recover the full "
            "GAG–POL–ENV architecture of ERVs (HERV-K, "
            "-W, -L …), not just the pol enzymes.<br>"
            "<b>Structural evidence</b> — heuristic terminal-repeat detectors (no external database): LTR by k-mer "
            f"seed + diagonal cluster (k={_SIG(structural_mod.find_ltr, 'k')}, "
            f"≥ {_SIG(structural_mod.find_ltr, 'min_ltr')} bp, "
            f"≥ {structural_mod.MIN_LTR_IDENTITY:g}% identity, "
            f"≥ {_SIG(structural_mod.find_ltr, 'min_anchors')} anchors — <b>this is a detection threshold, not "
            "a calibrated age boundary</b>. Paired LTRs generally diverge after insertion, so older or more "
            "diverged copies may be missed preferentially; identity also depends on lineage-specific rates, "
            "gene conversion, selection, sequence quality and detector sensitivity. A pair that fails only that floor is "
            "listed separately as an advisory candidate with its measured identity, so 'below threshold' is "
            "distinguishable from 'nothing found'. Further down, the k-mer seeding itself stops finding the "
            "pair at all, and no threshold change recovers that: where this happens depends on how long the "
            "repeat is, because a longer repeat carries more exact seed k-mers at the same identity — "
            "measured on a 574 bp LTR, candidates were still reported at 65% identity and none at 60%); "
            "TIR by a terminal inverted-repeat "
            f"scan plus a k-mer-vs-reverse-complement search, accepted at ≥ {structural_mod.MIN_TIR_IDENTITY:g}% "
            "identity — the same kind of floor as the LTR one, gating every DNA-transposon call the same way; "
            f"poly-A/poly-T tail ≥ 8 bp; TSD as a "
            f"{structural_mod.MIN_TSD}–{structural_mod.MAX_TSD} bp exact "
            "flanking direct repeat. For an LTR element the two retroviral cis-elements are also scanned: <b>PBS</b>, "
            f"the best reverse-complement match to the bundled primer-tRNA template — {', '.join(structural_mod._PRIMER_TRNA)} "
            "only, so an element primed by another tRNA (for example the initiator tRNA-Met of Ty1/copia) is not "
            f"assessed and no PBS row there is not 'no PBS' — within the {_SIG(structural_mod.find_pbs, 'search')} nt "
            "leader after the 5′ LTR, reported at ≥ 55% identity but only NAMED at ≥ 72% — a weaker match is reported as priming tRNA "
            "undetermined and marked <i>tentative</i>, because an endogenised PBS is diverged and a 55% match on random "
            "sequence is not rare; and <b>PPT</b>, a purine-dense run abutting the 3′ LTR (≥ 9 bp, ≥ 82% purine, ≤ 2 "
            "pyrimidine defects, 30 nt window). Coordinates are 0-based half-open, and each row lists its own detection method.<br>"
            "<b>Class / order / superfamily</b> — the Wicker&nbsp;2007" + w + " scheme, derived from the domain architecture and "
            "structural context; Copia vs Gypsy is called from the strand-aware integrase-vs-RT translation order, not "
            "ORF length. An ERV call requires paired terminal direct repeats together with RT and ENV. Without RT, only "
            "core capsid/matrix GAG plus PR, INT, or ENV preserves a hedged LTR fragment; paired direct repeats with an "
            "isolated domain remain a terminal direct-repeat candidate with <b>ORDER WITHHELD</b>.<br>"
            "<b>Reliability</b> — reported on two independent, citable axes rather than one fabricated number: (1) a "
            "per-domain call confidence from the HMMER i-Evalue (Eddy 2011); (2) a categorical domain-completeness "
            "tier — intact / near-complete / partial / structural candidate — mapped to the autonomous/intact criteria "
            "of Wicker 2007, TEsorter and LTR_retriever, and always scoped to the domain models actually tested. A "
            "completed search with no accepted hit says <i>none detected in the tested panel</i>; a failed search is "
            "<i>unavailable</i> and unassessed, never an absence, and a dependent repeat-only call is withheld.<br>"
            "<b>Family naming</b> (optional, WSL backend) — RepeatMasker (RMBLAST) against the Dfam&nbsp;4.0" + dfam +
            " library (curated families by default; the uncurated partitions only when you set Library to include "
            "them — the result header names which set was searched, and Dfam coverage depends on the lineage); "
            "this is the only step that makes a database family call, and it is absent from the offline path.")

    @staticmethod
    def _domain_panel_html():
        """The panel sentence, DERIVED from the profile table the scan actually loads.

        This text used to be hand-written and said "21 models" long after the panel had grown to 30 — the
        app under-reported its own method, which no test could catch because the claim lived in a UI string.
        Built from domains.DOMAIN_INFO, it cannot disagree with the models that were searched."""
        from collections import OrderedDict
        groups = OrderedDict()
        for _hmm, (code, _label, _cls, pfam) in domains_mod.DOMAIN_INFO.items():
            groups.setdefault(code, []).append(pfam)
        parts = [f"{c}&nbsp;{'/'.join(p)}" for c, p in groups.items()]
        return (f"{len(domains_mod.DOMAIN_INFO)} models, all CC0, in {len(groups)} reported groups: "
                + ", ".join(parts))

    # Detector parameters for the evidence types whose backend record carries no "method" key (the values are
    # structural.py's own defaults, which detect_all uses). Without this the Method column — glossed "How TEagle
    # detected this feature" — rendered EMPTY for PBS, PPT, TSD and the poly-A/T tails.
    _STRUCT_METHOD = {
        "PBS": (f"reverse-complement match to the bundled primer-tRNA template ({', '.join(structural_mod._PRIMER_TRNA)} "
                f"only; {len(next(iter(structural_mod._PRIMER_TRNA.values())))} nt) in the "
                f"{_SIG(structural_mod.find_pbs, 'search')} nt leader (≥ {_SIG(structural_mod.find_pbs, 'min_ident') * 100:g}% identity)"),
        "PPT": (f"purine-run extension from the 3′-LTR boundary (≥ {_SIG(structural_mod.find_ppt, 'min_len')} bp, "
                f"≥ {_SIG(structural_mod.find_ppt, 'min_purine') * 100:g}% purine, "
                f"≤ {_SIG(structural_mod.find_ppt, 'max_defects')} defects, "
                f"{_SIG(structural_mod.find_ppt, 'window')} nt window)"),
        "TSD": f"exact {structural_mod.MIN_TSD}–{structural_mod.MAX_TSD} bp direct repeat flanking the element",
        "poly-A": f"terminal homopolymer run (≥ {_SIG(structural_mod.find_polya, 'min_run')} bp, no mismatches)",
        "poly-T": f"terminal homopolymer run (≥ {_SIG(structural_mod.find_polya, 'min_run')} bp, no mismatches)",
        "polyA-signal": ("poly(A)-signal hexamer panel in the 3′ LTR, gated on a GU/U-rich downstream "
                         "element 20–60 nt past the hexamer (≥ 65% G+T, ≥ 4 T)"),
    }

    def _struct_row(self, e, tdr=False, cl=None, base=None):
        """One table row. `tdr` relabels an LTR row of an order-withheld call (the call does not claim LTR origin);
        `cl` (the record's classification) lets a TSD row say whether its length fits the called superfamily; `base`
        overrides the on-screen coordinate numbering (0 gives the canonical row an export writes)."""
        fp, tp = e.get("five_prime") or e.get("upstream"), e.get("three_prime") or e.get("downstream")
        sp = None
        if fp and tp:                                     # a PAIR (LTR/TIR arms, or the two TSD copies): show BOTH,
            coords = f"{self._span(*fp, base)}  ·  {self._span(*tp, base)}"   # matching the blocks drawn in the genome viewer
        else:
            sp = e.get("pos") or e.get("upstream") or e.get("element_span") or [None, None]
            coords = self._span(sp[0], sp[1], base) if sp[0] is not None else ""
        if e.get("strand") == "-" and coords:             # record coordinates, element on the other strand
            coords += "  (−)"
        arm = e.get("ltr_len") or e.get("tir_len") or e.get("length") or ""
        if arm == "" and sp and sp[0] is not None:        # PBS carries no length key — take Len from its own span
            arm = sp[1] - sp[0]
        metric = (f"{e['identity']}%" if e.get("identity") is not None else e.get("motif", ""))
        tm = e.get("termini") or {}
        if tm and e["type"].startswith("LTR"):               # the terminal motif was computed but only ever shown as a tooltip
            metric += (f" · {tm.get('five_start', '')}…{tm.get('five_end', '')} ("
                       + {"canonical": "canonical", "non-canonical": "non-canonical motif"}.get(tm.get("motif_tier"), "no known motif") + ")")
        if e["type"].startswith(("poly-A", "poly-T")):
            metric = f"uninterrupted {e['type'][5]} run"       # find_polya counts an exact run, so there is no purity to report
        if e.get("note") and not e.get("confident", True):   # a hedged call (non-confident PBS) must not read as a
            metric += " · tentative"                          # hard number here; the full hedge is the cell tooltip
        if e.get("length_is_lower_bound"):                # Len stays the measured number; the hedge sits beside it
            metric += " · length is a lower bound"
        if e["type"].startswith("TSD"):
            if e.get("chance_level"):
                metric += " · at chance level"
            elif cl and e.get("length"):                  # a chance-level match is not diagnostic, so no badge for it
                cg = structural_mod.tsd_congruence(e["length"], cl.get("superfamily"), e.get("motif"))
                sf = (cl.get("superfamily") or "").split(" (")[0]
                if cg["verdict"] == "congruent":
                    metric += f" · congruent length ({cg['expected']} bp for {sf})"
                elif cg["verdict"] == "incongruent":
                    metric += f" · incongruent: {sf} duplicates {cg['expected']} bp"
        shown = "TDR (terminal direct repeat)" if tdr and e["type"].startswith("LTR") else e["type"]
        return [shown, coords, arm, metric,
                e.get("method") or self._STRUCT_METHOD.get(e["type"].split(" ")[0], "")]

    def _record_note(self, rec, res=None):
        """Which record an exported table describes: its public id, length and sequence digest."""
        ident = _record_public_identity(rec)
        seq = rec.get("seq") or ""
        return (f"Record {ident['id']} · {(rec.get('composition') or {}).get('length', len(seq)):,} bp"
                + (f" · sequence SHA-256 {sequtil.sha256(seq)}" if seq else "")
                + " · coordinates index this record's own sequence.")

    @staticmethod
    def _orf_number(d, orfs):
        """1-based ORF number (as in the ORF table) a domain hit lies in, or '—' when the hit cannot be placed on
        one (fail closed: strand and interval must both agree)."""
        i, nt = d.get("orf"), d.get("nt")
        if not isinstance(i, int) or isinstance(i, bool) or not 0 <= i < len(orfs) or not nt:
            return "—"
        o = orfs[i]
        return i + 1 if d.get("strand") == o.get("strand") and o["start"] <= nt[0] and nt[1] <= o["end"] else "—"

    def _searched_not_found(self, rec, cl, tdr=False):
        """Plain statements of what was searched for and not found, for hallmarks the called class would carry:
        found, absent and not-searched are three different results and a table row alone cannot tell them apart."""
        s = rec.get("structural") or []
        out = []
        ltr = next((e for e in s if e["type"].startswith("LTR")), None)
        core = ltr or next((e for e in s if e["type"].startswith("TIR")), None)
        if core and not any(e["type"].startswith("TSD") for e in s):
            state = str(core.get("tsd_state") or "")
            if state.startswith("not assessable"):
                out.append("TSD — not assessable: the record has no sequence beyond the element on at least one side, "
                           "so no target-site duplication could be looked for.")
            elif state.startswith("absent"):
                out.append(f"TSD — absent: flanking sequence is present on both sides, but no exact "
                           f"{structural_mod.MIN_TSD}–{structural_mod.MAX_TSD} bp direct repeat flanks the element.")
            elif state:
                out.append("TSD — " + state + ".")
        if ltr:
            word = "TDR" if tdr else "LTR"
            if ltr.get("cis_orientation") == "undetermined":
                out.append(ltr.get("cis_note") or "The orientation of the element is undetermined, so PBS, PPT and "
                                                  "polyA-signal are not reported.")
            else:
                if not any(e["type"].startswith("PBS") for e in s):
                    scope = ltr.get("pbs_scope") or (
                        f"PBS search covers {', '.join(structural_mod._PRIMER_TRNA)} only; an element primed by another "
                        "tRNA is not assessed, so no PBS row is not 'no PBS'")
                    out.append(f"PBS — no match in the {_SIG(structural_mod.find_pbs, 'search')} nt leader after the 5′ "
                               f"{word}. {scope[0].upper() + scope[1:]}.")
                if not any(e["type"].startswith("PPT") for e in s):
                    out.append(f"PPT — no purine tract met the criteria (≥ {_SIG(structural_mod.find_ppt, 'min_len')} bp, "
                               f"≥ {_SIG(structural_mod.find_ppt, 'min_purine'):.0%} purine, ≤ "
                               f"{_SIG(structural_mod.find_ppt, 'max_defects')} pyrimidine defects within "
                               f"{_SIG(structural_mod.find_ppt, 'window')} nt of the 3′ {word}).")
        return out

    def _struct_tips(self, e, cl=None):
        """Per-cell tooltip overrides for one structural row: carry the detector's own hedge (e.g. a PBS below
        the confident threshold) into the table, where it previously showed only in the genome-viewer hover tip."""
        # the polyA-signal note is a standing scope limit, not a per-call hedge, so it rides along even when
        # the motif passed its gate (confident=True) — a confident MOTIF is still not a cleavage site.
        always = str(e.get("type", "")).startswith("polyA-signal")
        note = e.get("note") if (e.get("note") and (always or not e.get("confident", True))) else None
        # The terminal-motif badge was computed but never surfaced anywhere the user could see it. It
        # belongs on the LTR row it describes: canonical TG…CA, one of the documented non-canonical
        # termini, or neither — with the standing caveat that absence is not evidence against the call.
        tm = e.get("termini") or {}
        if tm:
            tier = tm.get("motif_tier")
            obs = f"{tm.get('five_start','')}…{tm.get('five_end','')}"
            lead = ("termini " + obs + " — " +
                    {"canonical": "the canonical TG…CA integrase att motif",
                     "non-canonical": f"the non-canonical motif {tm.get('noncanonical_motif')}",
                     }.get(tier, "no known terminal motif"))
            note = (note + "\n\n" if note else "") + lead + ". " + (tm.get("note") or "")
        metric = note
        if e.get("length_is_lower_bound"):                # the Len cell is the measured value; this says what it is
            bound = "Length is a lower bound: " + (e.get("bound_reason") or "the repeat reaches the end of what was scanned")
        else:
            bound = None
        if str(e.get("type", "")).startswith("TSD") and e.get("chance_probability") is not None:
            chance = (f"Chance probability {e['chance_probability']:g} after {e.get('n_tests', 1)} test(s) of length and "
                      "boundary placement" + ("; at or above the chance threshold, so this match is not evidence of an "
                                              "insertion site" if e.get("chance_level") else "") + ".")
            if e.get("selected_for_motif"):
                chance += " The arm boundary was trimmed to reach this motif, so it is not an independent observation."
            metric = (note + "\n\n" if note else "") + chance
        return [note, None, bound, metric, None]


    # =================== sequence helpers / staleness ===================
    _RC = {"A": "T", "T": "A", "G": "C", "C": "G", "U": "A", "R": "Y", "Y": "R", "S": "S", "W": "W",
           "K": "M", "M": "K", "B": "V", "D": "H", "H": "D", "V": "B", "N": "N"}

    def _clean_seq(self, text=None):
        t = self.seq.toPlainText() if text is None else text
        body = "".join(l for l in t.splitlines() if not l.startswith(">"))
        return "".join(c for c in body if c.isalpha()).upper()

    def _norm_seq(self, text=None):
        """Normalize exactly as the backend (sequtil.parse_fasta record 0 + _norm) does, so a re-sliced
        feature indexes the same positions the engine reported. With a '>' header present, take ONLY the
        first record's sequence — drop any bases before the first header and any later records — then strip
        all whitespace, uppercase, U->T (keep gaps/'*'/digits). Without a header, the whole text is the seq."""
        t = self.seq.toPlainText() if text is None else text
        if ">" in t:
            body, started = [], False
            for line in t.replace("\r\n", "\n").replace("\r", "\n").split("\n"):
                if line.startswith(">"):
                    if started:                               # a second header ends record 0
                        break
                    started = True
                elif started:
                    body.append(line)
            body = "".join(body)
        else:
            body = t
        return "".join(body.split()).upper().replace("U", "T")

    def _operation_input(self, text=None):
        """Bind a downstream operation to the exact serialized fetched-source record, when available."""
        raw = self.state.get("seq", "") if text is None else text
        sequence = self._norm_seq(raw)
        source = dict(self.state.get("source") or {})
        if source.get("inputType") not in {"fetched_accession", "fetched_coordinate"}:
            return self._local_operation_text(raw, sequence), source
        digest = sequtil.sha256(sequence)
        record = self.state.get("last_rec") or {}
        matches = [(index, binding) for index, binding in enumerate(source.get("records") or [])
                   if binding.get("length") == len(sequence) and binding.get("sha256") == digest and
                   (not record.get("id") or binding.get("id") == record.get("id"))]
        if len(matches) != 1:
            return raw, source
        index, binding = matches[0]
        source["records"] = [dict(binding)]
        if source["inputType"] == "fetched_coordinate":
            regions = source.get("regions") or []
            if index >= len(regions):
                return raw, source
            source["regions"] = [dict(regions[index])]
            source["accession"] = source["regions"][0]["chrAccession"]
            source["sourceUrl"] = f"https://www.ncbi.nlm.nih.gov/nuccore/{source['accession']}"
        return f">{binding['id']}\n{sequence}", source

    def _local_operation_text(self, raw, sequence):
        """Text to submit for a pasted or uploaded record. When `sequence` is exactly the analysed record, it
        travels under that record's own header: the engine turns the header into the same public ID the
        analysis card shows (the header itself is never copied into a result), so primers, PCR, family and
        splice results are bound to the record they were run on. The header is accepted only if it reproduces
        the analysed identity and the bases match the sealed digest; anything else travels as given."""
        rec = self.state.get("selected_record") or {}
        heads = self.state.get("record_headers") or []
        i = self.state.get("record_index")
        if (rec.get("identitySource") == "local_input" and isinstance(i, int) and 0 <= i < len(heads) and
                rec.get("length") == len(sequence) and rec.get("sha256") == sequtil.sha256(sequence)):
            from teagle_core.gff3 import public_record_identity
            ident = public_record_identity(heads[i])
            if ident["id"] == rec.get("id") and ident.get("idSha256") == rec.get("idSha256"):
                return f">{heads[i]}\n{sequence}"
        return raw

    def _slice(self, s, e, seq=None):
        # panel-01 feature coords index the ANALYZED snapshot, not the live box (which splice/annotate submit);
        # an explicit seq (family/splice/amplicon) is used verbatim.
        base = self._norm_seq(seq) if seq is not None else self._norm_seq(self.state.get("analyzed_seq", ""))
        return base[s:e]

    def _stale_block(self):
        """Block primer / PCR when the sequence in the box differs from the analysed one."""
        cur = self._clean_seq()
        if not self.state.get("analyzed_clean"):
            self._banner("Run analysis first.", "info"); return True
        if cur and cur != self.state["analyzed_clean"]:
            self._banner("Sequence changed since analysis — Run analysis again before designing primers.", "info")
            return True
        return False

    def _copy(self, text):
        QApplication.clipboard().setText(text)
        QToolTip.showText(QCursor.pos(), "copied", self)     # brief feedback at the cursor, not the status chip

    def _save_text(self, text, base, ext):
        """Write a text payload to a user-chosen file, with the same confirm-and-extension behaviour the
        annotation export uses, so every 'export' in the app behaves the same way."""
        if not text:
            self._banner("There is nothing to export.", "info"); return
        path, _ = QFileDialog.getSaveFileName(self, f"Export {ext.upper()}", f"{base}.{ext}",
                                              f"{ext.upper()} (*.{ext})")
        if not path:
            return
        if not path.lower().endswith("." + ext):
            path += "." + ext
        def write(p):
            with widgets.atomic_write(p, "w", encoding="utf-8", newline="\n") as fh:
                fh.write(text if text.endswith("\n") else text + "\n")
        widgets.safe_export(self, path, write, ok=f"Written to {os.path.basename(path)}.")

    def _feat_menu(self, start, end, strand, label, protein=None, dna=None, src_seq=None, kind=None):
        """Right-click menu for a feature — CONTEXTUAL: only actions valid for the clicked item. Coordinates
        address `src_seq` (default: the currently selected record's analysed sequence); pass `dna` when the exact sequence is already known
        (an amplicon carries its own seq), so copies never re-slice the wrong template. `kind` (the feature type)
        gates the action items so a short structural motif is not offered primer design or splice routing."""
        explicit = dna is not None
        if not explicit:
            raw = self._slice(start, end, src_seq)
            dna = self._revcomp(raw) if strand == "-" else raw
        rev = "_rev" if strand == "-" else ""
        def _design():
            if explicit:                                  # amplicon: design within its own sequence
                self._design_for_domain(0, len(dna), label, seq=dna)
            elif src_seq is not None:                     # feature on a non-panel-01 sequence
                self._design_for_domain(start, end, label, seq=src_seq)
            else:
                self._design_for_domain(start, end, label)
        fid = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(label)).strip("_") or "feature"   # FASTA-safe id (label may hold '5′ flank')
        items = [(f"Copy FASTA", lambda: self._copy(f">{fid}_{start}-{end}{rev} {_COORD_SYSTEM}\n{dna}")),
                 (f"Copy DNA", lambda: self._copy(dna)),
                 (f"Copy coords ({self._span(start, end)} {strand}{', 1-based' if self.coord_base else ''})",
                  lambda: self._copy(f"{self._span(start, end).replace('–', '-')} {strand}"))]
        if protein:
            items.append(("Copy protein", lambda: self._copy(protein)))
        # action items are contextual: a structural motif (LTR / TIR / TSD / PBS / PPT / poly-A) gets Copy actions only —
        # a primer designed INSIDE a repeat would be non-specific, and "send to splice" on a motif is meaningless. A
        # coding/transcript-like feature that is long enough (room for two primers + a product) gets primer + splice + sub-region.
        long_enough = len(dna) >= 50
        structural_motif = str(kind or "").split(" ")[0] in ("LTR", "TIR", "TDR", "TSD", "PBS", "PPT", "PAS", "poly-A",
                                                             "poly-T", "polyA-signal", "Sub-threshold")
        if long_enough and not structural_motif:
            items.append(("Design primer here", _design))
            items.append(("Send to splice structure (03b)",
                          lambda: self._send_to_splice(f">{fid}_{start}-{end}{rev}\n{dna}")))
            items.append(("Select a sub-region → primer / splice…", lambda: self._subregion(dna, fid)))
        # Flanking sequence is offered for EVERY feature, including the structural motifs above: designing
        # primers in the flanks to amplify ACROSS an insertion is the standard way a bench scientist
        # genotypes one, and it is the case where designing inside the element would be wrong.
        if not explicit:
            items.append(("Flanking sequence (upstream / downstream)…",
                          lambda: self._flank_picker(start, end, strand, label, src_seq)))
        return items

    def _flank_picker(self, start, end, strand, label, src_seq=None):
        """Take the sequence upstream and/or downstream of a feature and copy, export or design on it.

        Sides are named in RECORD orientation (upstream = the 5' side of the sequence as loaded), which is
        unambiguous regardless of the feature's strand; a minus-strand feature says so on screen rather
        than silently swapping the meaning of the two words."""
        # the box text carries FASTA headers and line breaks; normalise so coordinates index the bases
        seq = self._norm_seq(src_seq if src_seq is not None else self.state.get("seq", ""))
        if not seq:
            self._banner("No sequence is loaded to take flanks from.", "info"); return
        n = len(seq)
        fid = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(label)).strip("_") or "feature"
        dlg = QDialog(self); dlg.setWindowTitle("Flanking sequence"); dlg.setModal(False)
        dlg.setAttribute(Qt.WA_DeleteOnClose, True)
        lay = QVBoxLayout(dlg)
        lay.setContentsMargins(theme_mod.sp(14), theme_mod.sp(12), theme_mod.sp(14), theme_mod.sp(12))
        lay.setSpacing(theme_mod.sp(8))
        head = QLabel(f"Flanks of {label} ({self._span(start, end)}, {strand} strand)")
        hf = head.font(); hf.setBold(True); head.setFont(hf); lay.addWidget(head)
        note = QLabel("“Upstream” is the 5′ side of the sequence as loaded and “downstream” the 3′ side, so the "
                      "two names do not change meaning with the feature's strand."
                      + ("  This feature is on the minus strand, so its own biological upstream is the "
                         "downstream side here." if strand == "-" else ""))
        note.setObjectName("orient"); note.setWordWrap(True); lay.addWidget(note)
        row = QHBoxLayout()
        side = QComboBox()
        side.addItem("Upstream only", "up"); side.addItem("Downstream only", "down")
        side.addItem("Both flanks (separate records)", "both")
        side.setToolTip("Which side of the feature to take. Both gives two FASTA records, never a "
                        "concatenation — joining them would create a junction that does not exist. Design primers "
                        "with Both places one primer in each flank so the product spans the feature.")
        spin = QSpinBox(); spin.setRange(1, 100000); spin.setValue(500); spin.setSuffix(" bp")
        spin.setToolTip("How many bases to take on each chosen side. Clipped at the ends of the record; "
                        "the actual length taken is shown below.")
        row.addWidget(QLabel("Side")); row.addWidget(side, 1)
        row.addWidget(QLabel("Length")); row.addWidget(spin, 1); lay.addLayout(row)
        info = QLabel(""); info.setObjectName("orient"); info.setWordWrap(True); lay.addWidget(info)

        def parts():
            L = spin.value(); which = side.currentData(); out = []
            if which in ("up", "both"):
                s = max(0, start - L)
                if s < start:
                    out.append((f"{fid}_upstream_{s}-{start}", seq[s:start]))
            if which in ("down", "both"):
                e = min(n, end + L)
                if end < e:
                    out.append((f"{fid}_downstream_{end}-{e}", seq[end:e]))
            return out

        def span_window():
            L = spin.value()
            return max(0, start - L), min(n, end + L)

        def both_problem():
            """Why a pair spanning the feature cannot be designed, or None. One primer must fit in each flank
            and the product (flank + feature + flank) must fit the product range set in panel 04."""
            s0, e0 = span_window()
            if s0 >= start or e0 <= end:
                return ("The feature reaches the end of the record on one side, so there is no flank for one "
                        "of the two primers.")
            try:
                prm = self._read_primer_params()
            except ValueError:
                return None                                   # the design path reports a bad parameter itself
            need = (end - start) + 2 * prm["min_size"]
            if need > prm["prod_max"]:
                return (f"The feature is {end - start:,} bp, longer than the maximum product size "
                        f"({prm['prod_max']:,} bp) minus room for two primers, so no pair can span it. Raise "
                        "the maximum product size in panel 04, or take one flank.")
            return None

        def refresh():
            p = parts()
            both = side.currentData() == "both"
            why = both_problem() if both else None
            dbtn.setEnabled(bool(p) and why is None)
            dbtn.setToolTip(
                "Place one primer in the upstream flank and one in the downstream flank, so the product spans "
                "the feature." if both else
                "Primers are placed inside this flank, so the product lies wholly in flanking sequence and "
                "amplifies the same whether or not the feature is inserted. Choose Both to span the feature.")
            if not p:
                info.setText("No flanking sequence is available on that side — the feature reaches the end "
                             "of the record."); return
            info.setText(" · ".join(f"{nm}: {len(s)} bp" for nm, s in p)
                         + ("" if all(len(s) == spin.value() for _, s in p)
                            else "  (clipped by the end of the record)")
                         + (f"  {why}" if why else ""))

        def fasta():
            return "\n".join(f">{nm}\n{s}" for nm, s in parts())

        def design():
            if side.currentData() == "both":
                if both_problem():
                    return
                s0, e0 = span_window()
                self._design_for_domain(start - s0, end - s0, f"{fid} (across the element)", seq=seq[s0:e0])
            else:
                p = parts()
                if not p:
                    return
                nm, s = p[0]
                self._design_for_domain(0, len(s), nm, seq=s)
            dlg.accept()

        brow = QHBoxLayout(); brow.addStretch(1)
        dbtn = None
        for txt, fn in (("Copy FASTA", lambda: self._copy(fasta())),
                        ("Export FASTA…", lambda: self._save_text(fasta(), fid + "_flank", "fasta")),
                        ("Design primers", design)):
            b = QPushButton(txt); b.setProperty("sm", True); b.clicked.connect(fn); brow.addWidget(b)
            dbtn = b
        side.currentIndexChanged.connect(lambda *_: refresh()); spin.valueChanged.connect(lambda *_: refresh())
        refresh()
        close = QPushButton("Close"); close.setProperty("sm", True); close.clicked.connect(dlg.reject)
        brow.addWidget(close); lay.addLayout(brow)
        self._uppercase_buttons(dlg)
        dlg.adjustSize(); dlg.show(); dlg.raise_()
        self._flank_dlg = dlg

    def _subregion(self, dna, fid):
        """Pick a sub-interval WITHIN a feature (1-based inclusive, offsets into the feature's own sequence) and
        send only that subset to primer design or splice detection. Offset-based, so it is strand- and
        coordinate-space-safe: `dna` is already the feature's sequence in biological orientation."""
        n = len(dna)
        dlg = QDialog(self); dlg.setWindowTitle("Select a sub-region")
        dlg.resize(round(460 * theme_mod.UI_SCALE), round(200 * theme_mod.UI_SCALE))
        lay = QVBoxLayout(dlg)
        lay.addWidget(QLabel(f"<b>{fid}</b> is {n} bp. Choose the sub-interval (1-based, inclusive) to use:"))
        row = QHBoxLayout()
        s = QSpinBox(); s.setRange(1, n); s.setValue(1)
        e = QSpinBox(); e.setRange(1, n); e.setValue(n)
        lenlab = QLabel(); prev = QLabel(); prev.setObjectName("cardmeta"); prev.setWordWrap(True)
        def _upd():
            if e.value() < s.value():
                e.setValue(s.value())
            i, j = s.value() - 1, e.value()
            lenlab.setText(f"{j - i} bp  ·  0-based [{i}, {j})")     # echo the 0-based half-open span used everywhere else
            sub = dna[i:j]
            prev.setText(f"5′ {sub[:48]}{'…' if len(sub) > 48 else ''}")
        s.valueChanged.connect(_upd); e.valueChanged.connect(_upd); _upd()
        row.addWidget(QLabel("from")); row.addWidget(s); row.addWidget(QLabel("to")); row.addWidget(e)
        row.addWidget(lenlab); row.addStretch(1); lay.addLayout(row); lay.addWidget(prev)
        def _go(route):
            i, j = s.value() - 1, e.value()
            sub = dna[i:j]
            if len(sub) < 1:
                return
            sid = f"{fid}_{s.value()}-{e.value()}"
            dlg.accept()
            if route == "primer":
                self._design_for_domain(0, len(sub), sid, seq=sub)
            else:
                self._send_to_splice(f">{sid}\n{sub}")
        brow = QHBoxLayout(); brow.addStretch(1)
        bp = QPushButton("Design primers"); bp.setProperty("sm", True); bp.clicked.connect(lambda: _go("primer"))
        bs = QPushButton("Send to splice"); bs.setProperty("sm", True); bs.clicked.connect(lambda: _go("splice"))
        bc = QPushButton("Cancel"); bc.setProperty("sm", True); bc.clicked.connect(dlg.reject)
        brow.addWidget(bp); brow.addWidget(bs); brow.addWidget(bc); lay.addLayout(brow)
        self._uppercase_buttons()
        dlg.exec()

    def _struct_menu(self, e):
        """Right-click menu for a Structural-evidence row. What is copied must be what the ROW shows: an LTR/TIR
        row lists its two arms (Len = one arm), so each arm is copied on its own instead of silently substituting
        element_span — which is the whole element, and for a TIR the whole input. The whole element stays
        available as its own labelled group of items; the other rows address the single span they display."""
        fp, tp, span = e.get("five_prime"), e.get("three_prime"), e.get("element_span")
        up, dn = e.get("upstream"), e.get("downstream")
        lab = e["type"].split(" ")[0]
        strand = e.get("strand") or "+"                   # a minus-strand cis row addresses record coordinates

        def _grp(name, sp, fid, flanks=False):            # relabel every item so each names its span
            out = []
            for l, f in self._feat_menu(sp[0], sp[1], strand, fid, kind=e["type"]):
                if l.startswith("Flanking"):              # the arms' and copies' own flanks are the element's or internal
                    if flanks:
                        out.append((f"Flanks of the {name}…", f))
                else:
                    out.append((l.replace("Copy ", f"Copy {name} ", 1), f))
            return out
        if fp and tp:
            items = _grp("5′ arm", fp, f"{lab}_5prime") + _grp("3′ arm", tp, f"{lab}_3prime")
            if span and span[0] is not None:
                items += _grp("whole element", span, f"{lab}_element", flanks=True)
            return items
        if up and dn:                                     # a TSD is a pair: each copy is addressable
            return _grp("upstream copy", up, f"{lab}_upstream") + _grp("downstream copy", dn, f"{lab}_downstream")
        sp = e.get("pos") or up or span or [None, None]   # mirrors _struct_row's coord precedence
        if sp[0] is None:
            return [("Copy type", lambda: self._copy(e["type"]))]
        return self._feat_menu(sp[0], sp[1], strand, lab, kind=e["type"])

    def _region_menu(self, region, src_seq=None):
        """Right-click menu for a feature glyph in the genome viewer. The glyph's own strand, kind and protein
        are used (a minus-strand domain copies the coding strand) and coordinates address `src_seq`, the
        sequence the viewer was drawn on; region['dna'] is deliberately not passed, since an explicit sequence
        would switch off the flank picker and restrict primer design to the feature itself."""
        return self._feat_menu(region["start"], region["end"], region.get("strand") or "+",
                               region.get("label") or "feature", protein=region.get("protein"),
                               src_seq=src_seq, kind=region.get("kind"))
