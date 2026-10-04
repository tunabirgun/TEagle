"""In-silico PCR panel."""
from __future__ import annotations

import json
import re
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QGridLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton, QTextEdit, QVBoxLayout, QWidget)
import figures
import widgets
from widgets import FigurePanel
import coords
import theme as theme_mod
from teagle_core import sequtil, primers as primers_mod
from ui_common import (
    CollapsibleCard, EXPORT_HEADERS, GLOSS, PANEL, _COORD_SYSTEM, _amp_call, _amp_kind, _amplicon_public_source,
    _clear_layout, _empty, _export_table_btn, _note, _sl)


class PcrMixin:
    # ---------- 05 In-silico PCR ----------
    def _build_pcr_card(self):
        card = CollapsibleCard(PANEL["pcr"], "In-silico PCR", "predicts which fragments amplify")
        self.card_pcr = card
        card.bodylay.addWidget(_sl("Loaded primer pairs — one gel lane each (in order)"))
        self.pcrQueueBox = QVBoxLayout(); qb = QWidget(); qb.setLayout(self.pcrQueueBox)
        card.bodylay.addWidget(qb)
        srow = QHBoxLayout()
        self.pcrStageAll = QPushButton("+ stage all designed"); self.pcrStageAll.setProperty("sm", True)
        self.pcrStageAll.setEnabled(False); self.pcrStageAll.clicked.connect(self._pcr_stage_all)
        self.pcrClear = QPushButton("Clear"); self.pcrClear.setProperty("sm", True)
        self.pcrClear.setAccessibleName("Clear the in-silico PCR primer-pair queue")
        self.pcrClear.setEnabled(False); self.pcrClear.clicked.connect(self._pcr_clear)
        srow.addWidget(self.pcrStageAll); srow.addWidget(self.pcrClear); srow.addStretch(1)
        self.pcrOrderCopy, self.pcrOrderExport = self._order_buttons(
            lambda: [(self._pair_label(c), c) for c in self.state.get("pcrPairs", [])])
        for b in (self.pcrOrderCopy, self.pcrOrderExport):
            b.setEnabled(False); srow.addWidget(b)
        card.bodylay.addLayout(srow)
        card.bodylay.addWidget(_sl("Or check a primer pair you already have"))
        cgrid = QGridLayout(); cgrid.setHorizontalSpacing(8); cgrid.setVerticalSpacing(4)
        self.pcrCustomF = QLineEdit(); self.pcrCustomR = QLineEdit(); self.pcrCustomName = QLineEdit()
        for ed, ph, nm in ((self.pcrCustomF, "Forward primer, 5′→3′", "Forward primer sequence"),
                           (self.pcrCustomR, "Reverse primer, 5′→3′", "Reverse primer sequence"),
                           (self.pcrCustomName, "Name (optional)", "Primer pair name")):
            ed.setPlaceholderText(ph); ed.setAccessibleName(nm)
        cgrid.addWidget(self.pcrCustomF, 0, 0, 1, 2); cgrid.addWidget(self.pcrCustomR, 1, 0, 1, 2)
        cgrid.addWidget(self.pcrCustomName, 2, 0)
        self.pcrCustomAdd = QPushButton("+ add pair"); self.pcrCustomAdd.setProperty("sm", True)
        self.pcrCustomAdd.clicked.connect(self._add_custom_pair); self.pcrCustomR.returnPressed.connect(self._add_custom_pair)
        cgrid.addWidget(self.pcrCustomAdd, 2, 1)
        card.bodylay.addLayout(cgrid)
        card.bodylay.addWidget(_note("Primers from a paper or your own lab, in any IUPAC letters. A 5′ tail (a restriction "
                                     "site, for example) is fine: it is reported separately from the part that anneals. "
                                     "These pairs were not designed here, so the intended-product check does not apply.",
                                     "info"))
        grid = QGridLayout(); grid.setHorizontalSpacing(8)
        self._grid_fields(grid, [("pcrMM", "Max mismatches", "2"), ("pcrTP", "3′ exact bases", "5"),
                                 ("pcrPmin", "Prod min (bp)", "70"), ("pcrPmax", "Prod max (bp)", "1000")], cols=4)
        card.bodylay.addLayout(grid)
        card.bodylay.addWidget(_sl("Optional custom background (FASTA) — additional-product search"))
        self.pcrBg = QTextEdit(); self.pcrBg.setMaximumHeight(round(60 * theme_mod.UI_SCALE))
        self.pcrBg.setPlaceholderText("Optional: paste extra background sequence(s) to reveal additional predicted products.")
        self.pcrBg.setTabChangesFocus(True)
        card.bodylay.addWidget(self.pcrBg)
        row = QHBoxLayout()
        self.runPcrBtn = QPushButton("Run loaded pairs"); self.runPcrBtn.setProperty("primary", True)
        self.runPcrBtn.setEnabled(False); self.runPcrBtn.clicked.connect(self._run_pcr)
        row.addWidget(self.runPcrBtn)
        self.pcrHint = QLabel("load one or more pairs, then run"); self.pcrHint.setObjectName("kdim")
        row.addWidget(self.pcrHint); row.addStretch(1)
        card.bodylay.addLayout(row)
        self.pcrBody = QVBoxLayout(); pb = QWidget(); pb.setLayout(self.pcrBody)
        self.pcrBody.addWidget(_empty("Load a pair, then run pair-aware amplicon search."))
        card.bodylay.addWidget(pb)
        return card

    def _gel_menu(self, region):
        """Right-click menu for a gel band → copy the amplicon FASTA / coordinates. A genome-scan product has no
        sequence (isPcr's headers are all that is kept), so it offers only its locus, with the sequence name and
        strand. A band that merges several near-equal products copies its representative one, and says so."""
        a = region.get("amplicon") or {}
        pair = region.get("pair", "")
        start, end, seq = a.get("start"), a.get("end"), a.get("seq", "")
        # same call as the band tooltip, the amplicon table and its FASTA export (single-primer is its own
        # bucket; on/off only mean something with a design locus)
        kind = _amp_kind(a, region.get("has_locus", True))
        n = re.search(r"(\d+) products drawn as one band", region.get("tip") or "")
        band = f" — 1 of {n.group(1)} in this band" if n else ""
        items = []
        if seq:
            label = (f"amplicon_{self._pair_token(pair)}_{self._pair_token(_amplicon_public_source(a))}_{start}-{end}_"
                     f"{a.get('product_length', a.get('length',''))}bp_{kind}")
            items.append((f"Copy FASTA{band}", lambda: self._copy(f">{label} {_COORD_SYSTEM}\n{seq}")))
        if start is not None and end is not None:
            if pair == "genome":
                locus = f"{_amplicon_public_source(a)}:{start}-{end}({a.get('strand', '?')})"
                items.append((f"Copy locus, 1-based ({locus}){band}", lambda: self._copy(locus)))
            else:
                items.append((f"Copy coords ({self._span(start, end)}, {'1-based' if self.coord_base else '0-based, half-open'}){band}",
                              lambda: self._copy(self._span(start, end).replace("–", "-"))))
        return items

    # =================== PCR queue ===================
    @staticmethod
    def _pcr_key(c):
        seal = (c.get("design_provenance") or {}).get("manifestSha256", "")
        span = c.get("intended_product_span")
        return c["left_seq"] + "|" + c["right_seq"] + "|" + seal + "|" + json.dumps(span)

    def _add_pcr_pair(self, c):
        pairs = self.state.setdefault("pcrPairs", [])
        if not any(self._pcr_key(p) == self._pcr_key(c) for p in pairs):
            pairs.append(c)
        self._render_pcr_queue()
        self._pending_domain = None                       # a newer explicit navigation cancels any pending routed re-anchor
        self._scroll_to(self.card_pcr)                    # send-to-PCR lands the user on the in-silico PCR panel

    def _add_custom_pair(self):
        """Queue a primer pair typed or pasted by the user. It carries no design seal, so the run assesses the
        products it finds on the specimen but makes no claim about an intended target."""
        try:
            f = primers_mod.normalize_primer(self.pcrCustomF.text(), "Forward primer")
            r = primers_mod.normalize_primer(self.pcrCustomR.text(), "Reverse primer")
        except ValueError as e:
            return self._banner(str(e))
        if min(len(f), len(r)) < primers_mod.MIN_PRIMER_SIZE:
            return self._banner(f"Each primer needs at least {primers_mod.MIN_PRIMER_SIZE} bases, tail included.")
        name = self.pcrCustomName.text().strip()
        n = sum(bool(p.get("custom")) for p in self.state.get("pcrPairs", [])) + 1
        self._add_pcr_pair({"id": name or f"custom{n}", "left_seq": f, "right_seq": r, "custom": True})
        for ed in (self.pcrCustomF, self.pcrCustomR, self.pcrCustomName):
            ed.clear()

    def _pcr_stage_all(self):
        for c in self.state.get("candidates", []):
            self._add_pcr_pair(c)

    def _pcr_clear(self):
        # Clear empties the QUEUE only. A completed run's gel, amplicon table and both export buttons stay —
        # discarding them here destroyed finished work; the caption below is what makes the state honest instead.
        self.state["pcrPairs"] = []
        self._render_pcr_queue()

    def _render_pcr_queue(self):
        while self.pcrQueueBox.count():
            w = self.pcrQueueBox.takeAt(0).widget()
            if w:
                w.setParent(None)
        pairs = self.state.get("pcrPairs", [])
        self.pcrClear.setEnabled(bool(pairs))
        self.pcrOrderCopy.setEnabled(bool(pairs)); self.pcrOrderExport.setEnabled(bool(pairs))
        self.runPcrBtn.setEnabled(bool(pairs))
        if not pairs:
            # never claim "nothing here" while a finished run is still rendered below — say which it is
            stale = bool(self.state.get("lastPcr"))
            self.pcrQueueBox.addWidget(_empty(
                "No pairs staged — the results below are from the previous run. Stage pairs and run again to replace them."
                if stale else
                "No pairs loaded. Design primers, then “send to in-silico PCR”, or stage all."))
            self.pcrHint.setText("load one or more pairs, then run")
            return
        for i, c in enumerate(pairs):
            roww = QWidget(); rl = QHBoxLayout(roww); rl.setContentsMargins(0, 0, 0, 0)
            # A queued row normally comes from Primer3 and has every field. Keep rendering recoverable if
            # an older session/test fixture supplies only the primer sequences: the queue is user-facing
            # state, so one absent display field must not turn the startup callback into a traceback.
            product = c.get("product_size")
            product_text = (f"{product} bp" if product is not None else
                            f"your own pair ({c.get('id')})" if c.get("custom") else "product size unavailable")
            lab = QLabel(f"{self._pair_label(c)}  F {c.get('left_seq', '')}  R {c.get('right_seq', '')}  · {product_text}")
            lab.setObjectName("cardmeta"); lab.setWordWrap(True); lab.setTextInteractionFlags(Qt.TextSelectableByMouse)
            lab.setMinimumWidth(1); rl.addWidget(lab, 1)
            for glyph, name, step, fn in (("↑", "Move this pair up", -1, self._move_pair), ("↓", "Move this pair down", 1, self._move_pair),
                                           ("×", "Remove this pair from the queue", None, None)):
                b = QPushButton(glyph); b.setProperty("sm", True); b.setAccessibleName(name); b.setToolTip(name)
                b.clicked.connect((lambda _=False, k=i, d=step, f=fn: f(k, d)) if fn else (lambda _=False, k=i: self._remove_pair(k)))
                rl.addWidget(b)
            self.pcrQueueBox.addWidget(roww)
        self._uppercase_buttons(self.pcrQueueBox.parentWidget())   # queue rows are rebuilt on every mutation
        self.pcrHint.setText(f"{len(pairs)} pair(s) loaded · run to search")

    def _move_pair(self, i, d):
        a = self.state.get("pcrPairs", [])
        j = i + d
        if 0 <= j < len(a):
            a[i], a[j] = a[j], a[i]
            self._render_pcr_queue()

    def _remove_pair(self, i):
        a = self.state.get("pcrPairs", [])
        if 0 <= i < len(a):
            a.pop(i)
            self._render_pcr_queue()                      # queue only; a finished run's gel + exports survive (see _pcr_clear)

    # =================== run in-silico PCR ===================
    def _run_pcr(self):
        self._clear_banner()
        if self._genome_inflight or self._genome_prep_inflight:
            return self._banner("A whole-genome scan / download is running — wait for it before an in-silico PCR run.", "info")
        if self._stale_block():
            return
        pairs = self.state.get("pcrPairs", [])
        if not pairs:
            return self._banner("Load at least one primer pair first.", "info")
        try:
            p = {"max_mm": self._primer_number("pcrMM", "Maximum mismatches", integer=True, minimum=0),
                 "tp": self._primer_number("pcrTP", "3′ exact bases", integer=True, minimum=0),
                 "prod_min": self._primer_number("pcrPmin", "Minimum PCR product size", integer=True, minimum=1),
                 "prod_max": self._primer_number("pcrPmax", "Maximum PCR product size", integer=True, minimum=1)}
            if p["prod_min"] > p["prod_max"]:
                raise ValueError("Minimum PCR product size cannot exceed maximum PCR product size")
            if any(len(str(c.get(k, ""))) < primers_mod.MIN_PRIMER_SIZE
                   for c in pairs for k in ("left_seq", "right_seq")):
                raise ValueError(f"Each primer must contain at least {primers_mod.MIN_PRIMER_SIZE} bases")
            shortest_primer = min(len(str(c[k])) for c in pairs for k in ("left_seq", "right_seq"))
            if p["tp"] > shortest_primer or p["max_mm"] > shortest_primer:
                raise ValueError("3′ exact bases and maximum mismatches cannot exceed the shortest primer length")
        except ValueError as e:
            return self._banner(str(e))
        # The size window filters the products found, so a designed product outside it would silently vanish
        # ('zero products') even though the pair amplifies its own template. Widen the window to hold every
        # designed product, in the open: the fields show the new values and the result says what was done.
        sizes = [(self._pair_label(c), c["product_size"]) for c in pairs if isinstance(c.get("product_size"), int)]
        window_note = None
        if sizes:
            lo, hi = min(p["prod_min"], min(s for _, s in sizes)), max(p["prod_max"], max(s for _, s in sizes))
            if (lo, hi) != (p["prod_min"], p["prod_max"]):
                out = [f"{n} ({s:,} bp)" for n, s in sizes if not p["prod_min"] <= s <= p["prod_max"]]
                window_note = (f"The designed product of {', '.join(out)} lies outside the product-size window "
                               f"{p['prod_min']:,}–{p['prod_max']:,} bp, which would have hidden it, so the window was "
                               f"widened to {lo:,}–{hi:,} bp for this run.")
                p["prod_min"], p["prod_max"] = lo, hi
                self.pfields["pcrPmin"].setText(str(lo)); self.pfields["pcrPmax"].setText(str(hi))
        labels = []
        for c in pairs:                                        # one unique, design-derived label per lane
            base = lab = self._pair_label(c); n = 2
            while lab in labels:
                lab = f"{base}#{n}"; n += 1
            labels.append(lab)
        lanes_short = []
        for lab in labels:                                     # a gel lane has room for about ten characters: the table and
            s = lab if len(lab) <= 10 else lab[:9] + "…"       # exports keep the full name, only the lane caption is shortened
            base, n = s, 2
            while s in lanes_short:
                s = f"{base[:8]}#{n}"; n += 1
            lanes_short.append(s)
        bg = self.pcrBg.toPlainText()
        # _tmpl_sig uses _norm_seq while the stale-gate uses _clean_seq; a post-analysis edit adding a gap/digit
        # could pass the gate yet flip this sig -> spurious OFF-target (fails safe, never a false on-target).
        operation_sequence, operation_source = self._operation_input(self.state["seq"])
        cur_sig = self._norm_seq(operation_sequence)           # on-target only when the pair was designed on THIS template
        assay_records = sequtil.parse_fasta(operation_sequence)
        assay_record = assay_records[0] if assay_records else None
        self._pcr_gen += 1                                     # new batch id: results from a superseded batch are dropped
        self._pcr_run = {"gen": self._pcr_gen, "results": [None] * len(pairs), "window_note": window_note,
                         "rec": self.state.get("record_index"), "digest": self._selected_digest(),
                         "pairs": [{"label": lab, "lane": sh, "fwd": c["left_seq"], "rev": c["right_seq"]}
                                   for lab, sh, c in zip(labels, lanes_short, pairs)]}
        self.runPcrBtn.setEnabled(False); widgets.set_caps(self.runPcrBtn, "◴ running…")
        if window_note:
            self._banner(window_note, "warn")
        for i, c in enumerate(pairs):
            # Geometry is assessed only on the exact record sealed by this candidate's design manifest.
            # A family/splice/amplicon candidate does not index state["seq"], so all design geometry is withheld.
            sealed = isinstance(c.get("design_provenance"), dict)
            design_records = sequtil.parse_fasta(c.get("_design_sequence", "")) if sealed else []
            bound = bool(sealed and design_records and assay_record == design_records[0] and
                         c.get("_tmpl_sig") == cur_sig)
            body = {"sequence": operation_sequence, "background": bg, "fwd": c["left_seq"], "rev": c["right_seq"],
                    "primer3_target": (c.get("primer3_target") if sealed else None),
                    "included_region": (c.get("included_region") if sealed else None),
                    "intended_product_span": (c.get("intended_product_span") if sealed else None),
                    "intended_endpoint_tolerance_bp": c.get(
                        "intended_endpoint_tolerance_bp", 0),
                    "design_provenance": (c.get("design_provenance") if sealed else None),
                    "design_sequence": (c.get("_design_sequence") if sealed else None),
                    "design_assessment": ("assess" if bound else
                                          ("withheld_template_mismatch"
                                           if sealed else "not_run")),
                    "source": operation_source, "params": p}
            self.engine.submit("pcr", body, key=f"pcr#{self._pcr_gen}#{i}")

    def _tag_fetched_sources(self, amplicons):
        """A fetched specimen is identified by its public accession, which the engine rewrote as a hash when it
        did not mark the amplicon's source as validated. The accession is the one in the sealed source context,
        so mark exactly those sources and the table can show it."""
        src = self.state.get("source") or {}
        if src.get("inputType") not in {"fetched_accession", "fetched_coordinate"}:
            return
        ids = {b.get("id") for b in src.get("records") or []}
        for a in amplicons:
            base = str(a.get("source") or "")
            if base.endswith(" (input)") and base[:-len(" (input)")] in ids and not a.get("sourceIdentitySource"):
                a["sourceIdentitySource"] = "validated_fetch"

    def _pcr_slot(self, key, result):
        """Fill one batch slot (a success dict OR an error placeholder) and render once the whole batch
        is in. A late result from a superseded batch (different gen) is dropped, so the button re-enabling
        after one lane can never let a stale sibling corrupt the next run."""
        parts = key.split("#")                                # pcr#<gen>#<i>
        gen, i = int(parts[1]), int(parts[2])
        run = getattr(self, "_pcr_run", None)
        if not run or gen != run["gen"] or not (0 <= i < len(run["results"])):
            return                                            # stale / superseded batch
        run["results"][i] = result
        if any(r is None for r in run["results"]):
            return
        self.runPcrBtn.setEnabled(True); widgets.set_caps(self.runPcrBtn, "Run loaded pairs")
        self.card_pcr.expand()
        results = run["results"]
        # A capped in-silico search is not an exhaustive one, and the difference is invisible in the
        # product list. Say it before the user reads the gel, naming the pair, since "no other product"
        # is precisely the conclusion a truncated run cannot support.
        names = [x["label"] for x in run.get("pairs") or []] or [f"P{k+1}" for k in range(len(results))]
        cut = [names[i] for i, dd in enumerate(results) if isinstance(dd, dict) and dd.get("truncated")]
        if cut:
            self._banner(f"In-silico PCR stopped early on {', '.join(cut)} — the product list is a lower "
                         "bound on priming, not an exhaustive product search. Repetitive templates hit "
                         "the search cap.", level="warn")
        lanes, amps, provs = [], [], []
        statuses, cautions = [], []
        for idx, dd in enumerate(results):
            lane = names[idx]
            mechanisms = dd.get("amplicons", [])
            self._tag_fetched_sources(mechanisms)
            status = "failed" if dd.get("error") else dd.get("status", "complete")
            statuses.append({"pair": lane, "status": status, "error": dd.get("error")})
            cautions.extend(f"{lane}: {c}" for c in dd.get("cautions") or [])
            physical = {}
            for a in mechanisms:
                key = (a.get("product_key") or
                       f"{a.get('source')}:{a.get('start')}:{a.get('end')}:"
                       f"{a.get('product_length', a.get('length'))}")
                prior_mechanisms = physical.get(key, {}).get("mechanism_keys", [])
                if (key not in physical or
                        (a.get("matches_intended_product") is True and
                         physical[key].get("matches_intended_product") is not True) or
                        (physical[key].get("single_primer") and not a.get("single_primer"))):
                    physical[key] = dict(a)
                    physical[key]["mechanism_keys"] = prior_mechanisms
                physical[key]["mechanism_keys"].append(a.get("mechanism_key", "unknown"))
            products = list(physical.values())
            for a in products:
                a["mechanism_count"] = len(a["mechanism_keys"])
            lanes.append({"label": ((run.get("pairs") or [{}] * len(results))[idx].get("lane") or lane), "amplicons": products,
                          "has_locus": any(a.get("on_target") is not None for a in products),
                          "status": status})
            for a in mechanisms:
                amps.append({**a, "pair": lane})
            if dd.get("provenance"):
                provs.append(dd["provenance"])
        self.state["lastPcr"] = {"lanes": lanes, "amplicons": amps, "statuses": statuses,
                                 "pairs": run.get("pairs") or [], "window_note": run.get("window_note"),
                                 "cautions": cautions, "digest": run.get("digest")}
        self._render_pcr(lanes, amps, statuses)
        self._uppercase_buttons()
        if provs:
            self._render_provenance(provs[0])
            # the card shows the first pair's seal, but a batch has one manifest per pair — the EXPORT must
            # cover every lane on the gel, not just pair 1. Store the full list for _export_manifest.
            self.state["prov_manifest_all"] = provs

    def _render_pcr(self, lanes, amps, statuses=None):
        _clear_layout(self.pcrBody)
        amplicon_table = None
        gel = FigurePanel(lambda bg, L=lanes: figures.svg_gel({"lanes": L}, bg),
                          "TEagle_gel", modes=("dark", "white", "uv", "mono"),
                          hit_regions=figures.gel_regions({"lanes": lanes}), on_menu=self._gel_menu)
        gel.apply_app_theme(self.theme)                   # gel opens in the current app theme (uv/mono stay manual)
        gel.setMinimumHeight(round(420 * theme_mod.UI_SCALE))
        self.pcrBody.addWidget(gel)
        if amps:
            headers = ["Pair", "Mechanism", "Source", self._ch("Coords"), "Template span", "Product len", "Mism F/R",
                       "Tail F/R", "Anneal F/R", "Call"]
            fc = theme_mod.FLAG[self.theme]
            rows, canon, styles, tips, exports = [], [], [], [], []
            for a in amps:
                tails = (a.get("fwd_tail5") or 0, a.get("rev_tail5") or 0)
                mm = self._per_primer(a, "fwd_mm", "rev_mm")
                orient = str(a.get("primer_orientation") or "F+R").replace("+", " + ")
                row = lambda b: [a["pair"], f"{orient}, {self._span(a['start'], a['end'], b)}", _amplicon_public_source(a),
                                 self._span(a["start"], a["end"], b), a["length"], a.get("product_length", a["length"]),
                                 mm, self._per_primer(a, "fwd_tail5", "rev_tail5"),
                                 self._per_primer(a, "fwd_anneal_len", "rev_anneal_len"), _amp_call(a)]
                rows.append(row(self.coord_base))
                canon.append([a.get("mechanism_key") or r if j == 1 else r for j, r in enumerate(row(0))])
                styles.append([None] * 7 + [fc["caution"] if any(tails) else None, None, None])
                exports.append([None, a.get("mechanism_key") or None] + [None] * 8)       # the raw key is the audit trail: export and tooltip
                tips.append([None, f"Mechanism key: {a.get('mechanism_key', '—')}"] + [None] * 4 + [
                    f"Mismatches inside the annealed 3′ core, forward / reverse primer: {mm}. Bases 5′ of the core are the "
                    "tail and are not counted here." + (" This product has unmatched 5′ bases: unless the primer was "
                                                       "designed with a tail, they are mismatches." if any(tails) else ""),
                    None, None, None])
            gl = {**GLOSS,
                  "Pair": "The designed primer pair that produced this amplicon: design batch and Primer3 pair number "
                          "(D2·P3), or the name you gave a pair of your own.",
                  "Mism F/R": "Mismatches between each primer and its binding site, forward / reverse primer, counted "
                              "over the annealed 3′ core only (the 3′ end is kept exact). Bases 5′ of the core are "
                              "reported as the tail, not as mismatches. A self-priming product lists both sites of "
                              "the one primer, e.g. 2,1/—.",
                  "Tail F/R": "Unmatched bases at the 5′ end of the forward / reverse primer, outside the annealed core. "
                              "Expected for a primer with a designed tail (a restriction site, an adapter); for a "
                              "primer without one they are mismatches.",
                  "Anneal F/R": "Length (nt) of the 3′ core of the forward / reverse primer that matches the template."}
            t = self._coord_table(headers, {**gl, **coords.gloss(self.coord_base)}, canon)
            t.export_header_map = {**EXPORT_HEADERS, "Mism F/R": "Mism F/R (annealed 3' core)",
                                   "Tail F/R": "Tail F/R (unmatched 5' nt)", "Anneal F/R": "Anneal F/R (nt)",
                                   "Coords (0-based)": "Coords (0-based, half-open)", "Coords (1-based)": "Coords (0-based, half-open)"}
            amplicon_table = t
            t.set_rows(rows, styles=styles, tips=tips, exports=exports)
            t.set_row_menu(lambda r: self._feat_menu(amps[r]["start"], amps[r]["end"], "+",
                                                     f"amplicon_{amps[r]['pair']}", dna=amps[r].get("seq", "")))
            t.setMaximumHeight(round(200 * theme_mod.UI_SCALE))
            self.pcrBody.addWidget(t)
            arow = QHBoxLayout(); arow.addStretch(1)
            arow.addWidget(_export_table_btn(t, "TEagle_amplicons", self))   # export the notebook table
            if any(a.get("seq") for a in amps):
                prim = {x["label"]: x for x in (self.state.get("lastPcr") or {}).get("pairs") or []}
                def _amps_fasta(_=False, aa=amps, tbl=t):
                    order = [tbl._orig(r) for r in range(tbl.rowCount())]      # follow the table's current (sorted) order
                    seq_amps = [aa[i] for i in order if 0 <= i < len(aa)] or aa
                    seen, out = {}, []
                    for a in seq_amps:
                        # source and orientation are in the ID: a product found in a custom background must not
                        # read as a specimen product, and two products at one position must not share an ID
                        ident = (f"amplicon_{self._pair_token(a['pair'])}_{self._pair_token(_amplicon_public_source(a))}_"
                                 f"{a['start']}-{a['end']}_{a.get('product_length', a['length'])}bp_"
                                 f"{self._pair_token(a.get('primer_orientation', 'F+R'))}_{_amp_kind(a)}")
                        seen[ident] = seen.get(ident, 0) + 1
                        ident += f"_{seen[ident]}" if seen[ident] > 1 else ""
                        pr = prim.get(a["pair"]) or {}
                        desc = (f" F={pr['fwd']} R={pr['rev']}" if pr else "")
                        out.append(f">{ident}{desc}\n{a.get('seq','')}")
                    widgets.save_fasta("\n".join(out), "TEagle_amplicons", self)
                cp = QPushButton("Export amplicons (FASTA)"); cp.setProperty("sm", True)   # export the sequences
                cp.clicked.connect(_amps_fasta)
                arow.addWidget(cp)
            self.pcrBody.addLayout(arow)
        else:
            partial_empty = any(s.get("status") == "partial_capped" for s in (statuses or []))
            failed = sum(1 for s in (statuses or []) if s.get("status") == "failed")
            if partial_empty:
                msg = "No product was recorded before the capped search stopped; unsearched candidates remain."
            elif failed:
                msg = f"No products are available; {failed} pair(s) failed and were not assessed."
            else:
                msg = "The completed search recorded zero products under the stated computational criteria."
            self.pcrBody.addWidget(_note(msg, "warn"))
        products = [a for lane in lanes for a in lane.get("amplicons", [])]
        intended_n = sum(1 for a in products if a.get("matches_intended_product") is True)
        on_n = sum(1 for a in products if a.get("on_target") is True)
        single_n = sum(1 for a in products if a.get("single_primer"))
        unassessed_n = sum(1 for a in products if a.get("on_target") is None and not a.get("single_primer"))
        partial = any(s.get("status") == "partial_capped" for s in (statuses or []))
        failed = sum(1 for s in (statuses or []) if s.get("status") == "failed")
        parts = (f"{len(lanes)} lane(s) · {len(amps)} amplification mechanism(s) · "
                 f"{len(products)} unique physical product(s) · {intended_n} intended product call(s) · "
                 f"{on_n} on-target")
        if single_n:
            parts += f" · {single_n} single-primer"
        if unassessed_n:
            parts += f" · {unassessed_n} target-unassessed"
        if failed:
            parts += f" · {failed} failed/not assessed"
        if partial:
            scope = (" Search caps were reached, so this is a lower bound on amplification mechanisms and no "
                     "specificity conclusion is supported.")
        elif failed:
            scope = (" At least one pair failed and was not assessed, so no batch-wide absence or specificity "
                     "conclusion is supported.")
        else:
            scope = " A complete local search does not establish experimental or genome-wide specificity."
        note = QLabel(f"{parts} · ladder lane “L”. Equal-sized physical products co-migrate into one gel band; "
                      f"the table retains each auditable primer-orientation mechanism. Table coordinates are 0-based, "
                      f"half-open.{scope}" if not self.coord_base else
                      f"{parts} · ladder lane “L”. Equal-sized physical products co-migrate into one gel band; "
                      f"the table retains each auditable primer-orientation mechanism. Coordinates are shown 1-based; "
                      f"exports stay 0-based, half-open.{scope}")
        note.setObjectName("orient"); note.setWordWrap(True); self.pcrBody.addWidget(note)
        lp = self.state.get("lastPcr") or {}
        broken = [f"{s['pair']} ({s.get('error') or 'no detail'})" for s in (statuses or []) if s.get("status") == "failed"]
        if broken:                                        # name each failed pair and why, not just how many
            self.pcrBody.addWidget(_note("Not assessed: " + "; ".join(broken) + ".", "warn"))
        if lp.get("window_note"):
            self.pcrBody.addWidget(_note(lp["window_note"], "info"))
        if lp.get("cautions"):
            self.pcrBody.addWidget(_note(" ".join(lp["cautions"]), "warn"))
        if amplicon_table is not None:
            amplicon_table.export_notes = [
                f"{parts}. Equal-sized physical products co-migrate into one gel band; the table retains each "
                f"auditable primer-orientation mechanism.{scope}",
                "Coordinates are 0-based, half-open positions on the searched sequence (the first base is 0; the end "
                "is not included). Mism, Tail and Anneal are listed forward / reverse primer; mismatches are counted "
                "over the annealed 3-prime core, and unmatched 5-prime bases are reported as the tail."]
            amplicon_table.export_notes += [f"Pair {x['label']}: forward 5'-{x['fwd']}-3', reverse 5'-{x['rev']}-3'."
                                            for x in lp.get("pairs") or []]
            if lp.get("window_note"):
                amplicon_table.export_notes.append(lp["window_note"])
            amplicon_table.export_notes += list(lp.get("cautions") or [])
        self.pcrBody.addLayout(self._next_step_row([
            (f"Scan a pair genome-wide → panel {PANEL['genome']}", lambda: self._scroll_to(self.card_genome))]))
        self._uppercase_buttons(self.pcrBody.parentWidget())   # export buttons are rebuilt with the results
