"""Primer design panel."""
from __future__ import annotations

import html
import math
import re
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QDialog, QGridLayout, QHBoxLayout, QLabel, QPushButton, QSizePolicy, QSpinBox, QVBoxLayout, QWidget)
import widgets
from widgets import BusyBar, DataTable
import coords
import theme as theme_mod
from teagle_core import sequtil, primers as primers_mod
from ui_common import (
    CollapsibleCard, GLOSS, PANEL, PRIMER_COLS, REFLINKS, _Combo, _clear_layout, _dg_clauses, _export_table_btn,
    _external, _kb_links, _metric_cell, _note, _record_public_id, _repeat_word, _safe_export_stem, _sl)


class PrimersMixin:
    # ---------- 04 Primer design ----------
    def _build_primer_card(self):
        card = CollapsibleCard(PANEL["primer"], "Primer design", "designs PCR primers (Primer3)")
        self.card_primer = card
        self._build_assay_goal(card)
        prow = QHBoxLayout()
        prow.addWidget(QLabel("Preset"))
        self.pPreset = _Combo()
        self._preset_keys = ["standard", "qpcr", "highspec", "permissive", "custom"]
        self.pPreset.addItems(["Standard PCR", "qPCR (short amplicon)", "High-specificity", "Permissive (hard targets)",
                               "Custom (edited values)"])
        self.pPreset.setAccessibleName("Primer design preset")
        # `activated` fires on every user pick, so choosing the preset that is already shown re-applies it
        self.pPreset.activated.connect(lambda i: self._apply_preset(self._preset_keys[i]))
        prow.addWidget(self.pPreset); prow.addStretch(1)
        card.bodylay.addLayout(prow)
        # basic params
        self.pfields = {}
        grid = QGridLayout(); grid.setHorizontalSpacing(8); grid.setVerticalSpacing(4)
        basic = [("pMin", "Prod min (bp)", "150"), ("pMax", "Prod max (bp)", "500"), ("pTm", "Opt Tm (°C)", "60"),
                 ("pMinS", "Min size", "18"), ("pMaxS", "Max size", "27"), ("pNum", "Return", "5")]
        self._grid_fields(grid, basic, cols=3)
        card.bodylay.addLayout(grid)
        # advanced params
        adv = [("pOptS", "Opt size", "20"), ("pTmMin", "Tm min", "57"), ("pTmMax", "Tm max", "63"),
               ("pGcMin", "GC min %", "40"), ("pGcMax", "GC max %", "60"), ("pPolyX", "Max poly-X", "4"),
               ("pGcClamp", "GC clamp", "0")]
        self.advBox = QWidget(); ag = QGridLayout(self.advBox); ag.setContentsMargins(0, 0, 0, 0)
        ag.setHorizontalSpacing(8); ag.setVerticalSpacing(4)
        self._grid_fields(ag, adv, cols=4)
        self.advToggle = QPushButton("▸ Advanced parameters"); self.advToggle.setProperty("link", True)
        self.advToggle.setAccessibleName("Advanced parameters")
        self.advToggle.setAccessibleDescription("Disclosure toggle — shows or hides the advanced primer parameters.")
        self.advToggle.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Fixed)
        self.advBox.setVisible(False)
        self.advToggle.clicked.connect(lambda: (self.advBox.setVisible(not self.advBox.isVisible()),
                                                self.advToggle.setText(("▾" if self.advBox.isVisible() else "▸") + " Advanced parameters")))
        card.bodylay.addWidget(self.advToggle)
        card.bodylay.addWidget(self.advBox)
        row = QHBoxLayout()
        self.designBtn = QPushButton("Design primers"); self.designBtn.setProperty("primary", True)
        self.designBtn.setEnabled(False); self.designBtn.clicked.connect(self._design_for_goal)
        row.addWidget(self.designBtn)
        rb = QPushButton("Reset"); rb.setProperty("sm", True); rb.clicked.connect(self._reset_primer_params)
        rb.setAccessibleName("Reset primer parameters to the standard preset")
        row.addWidget(rb)
        self.designHint = QLabel("run analysis first"); self.designHint.setObjectName("kdim")
        row.addWidget(self.designHint); row.addStretch(1)
        card.bodylay.addLayout(row)
        self.primBody = QVBoxLayout(); pb = QWidget(); pb.setLayout(self.primBody)
        card.bodylay.addWidget(pb)
        return card

    @staticmethod
    def _pair_label(c):
        """The name a primer pair keeps everywhere it appears: its design batch tag plus Primer3's own ID
        ('D2·P3'), so a lane, a table row and an export can be traced to one designed pair. A pair the user typed
        keeps its own name."""
        cid = str(c.get("id") or "pair")
        return f"{c['_batch']}·{cid}" if c.get("_batch") else cid

    @staticmethod
    def _pair_token(text):
        """A label as a filename / FASTA-header token ('D2·P3' -> 'D2-P3')."""
        return re.sub(r"[^A-Za-z0-9_.]+", "-", str(text)).strip("-") or "pair"

    @staticmethod
    def _per_primer(a, left, right):
        """'F/R' text for a measure the engine records per SITE (left, right). The primer at each site is named by
        left_primer / right_primer, so an R+F product's right site is reported under F, and a self-priming F+F
        product lists both of F's sites ('2,1/—') instead of passing one primer off as two."""
        got = {"F": [], "R": []}
        for site, dflt, key in (("left_primer", "F", left), ("right_primer", "R", right)):
            v = a.get(key)
            got.setdefault(a.get(site) or dflt, []).append("—" if v is None else str(v))
        return "/".join(",".join(got[k]) or "—" for k in "FR")

    def _src_link(self, key, url=None, label="source"):
        """A small clickable citation label that opens the verified DOI in a browser."""
        r = REFLINKS.get(key)
        lab = QLabel()
        if not r:
            return lab
        accent = theme_mod.ACCENT[self.theme]
        lab.setTextFormat(Qt.RichText)
        lab.setText(f'<a href="{url or r["url"]}" style="color:{accent};text-decoration:underline">'
                    f'{html.escape(label)}</a>')
        _external(lab)
        _kb_links(lab)
        lab.setObjectName("srclink")
        lab.setToolTip("Source — " + r["cite"])
        lab.setAccessibleName("Source citation — " + r["cite"])
        return lab

    def _src_html(self, key, url=None):
        """Inline 'source' anchor for RichText labels (set openExternalLinks on the label)."""
        r = REFLINKS.get(key)
        if not r:
            return ""
        accent = theme_mod.ACCENT[self.theme]
        return (f' <a href="{url or r["url"]}" style="color:{accent};text-decoration:underline" '
                f'title="Source — {r["cite"]}">source</a>')

    # =================== primer params ===================
    PRESETS = {
        "standard":   dict(pMin=150, pMax=500, pTm=60, pMinS=18, pMaxS=27, pNum=5, pOptS=20, pTmMin=57, pTmMax=63, pGcMin=40, pGcMax=60, pPolyX=4, pGcClamp=0),
        "qpcr":       dict(pMin=70, pMax=150, pTm=60, pMinS=18, pMaxS=24, pNum=5, pOptS=20, pTmMin=58, pTmMax=62, pGcMin=40, pGcMax=60, pPolyX=4, pGcClamp=1),
        "highspec":   dict(pMin=150, pMax=500, pTm=62, pMinS=20, pMaxS=26, pNum=8, pOptS=22, pTmMin=60, pTmMax=64, pGcMin=45, pGcMax=60, pPolyX=3, pGcClamp=2),
        "permissive": dict(pMin=100, pMax=1000, pTm=58, pMinS=18, pMaxS=30, pNum=10, pOptS=20, pTmMin=52, pTmMax=65, pGcMin=30, pGcMax=70, pPolyX=5, pGcClamp=0),
    }

    def _apply_preset(self, name):
        p = self.PRESETS.get(name)
        if not p:
            return
        for k, v in p.items():
            if k in self.pfields:
                self.pfields[k].setText(str(v))
        self._show_preset(name)

    # Plain-language help and units for the primer and in-silico PCR fields (label, tooltip). Applied by
    # _decorate_primer_fields; the limits quoted here are the ones _read_primer_params enforces.
    _FIELD_HELP = {
        "pMin": ("Prod min (bp)", "Shortest PCR product to design, in base pairs. It must be at least twice the maximum primer size."),
        "pMax": ("Prod max (bp)", "Longest PCR product to design, in base pairs. Designing around a domain needs it to be at least the domain length plus two primers."),
        "pTm": ("Opt Tm (°C)", "Melting temperature Primer3 aims for, in °C, calculated for the reaction conditions listed under the table."),
        "pMinS": ("Min size (nt)", "Shortest primer, in nucleotides (18 or more)."),
        "pMaxS": ("Max size (nt)", "Longest primer, in nucleotides. Primer3 allows at most 36."),
        "pNum": ("Pairs to return", "How many primer pairs to return (1–50)."),
        "pOptS": ("Opt size (nt)", "Primer length Primer3 aims for, in nucleotides."),
        "pTmMin": ("Tm min (°C)", "Lowest acceptable primer melting temperature, in °C."),
        "pTmMax": ("Tm max (°C)", "Highest acceptable primer melting temperature, in °C."),
        "pGcMin": ("GC min %", "Lowest acceptable G+C content of a primer, in percent."),
        "pGcMax": ("GC max %", "Highest acceptable G+C content of a primer, in percent."),
        "pPolyX": ("Max poly-X (bases)", "Longest run of one repeated base (for example AAAAA) a primer may contain."),
        "pGcClamp": ("GC clamp (bases)", "Number of G or C bases required at the 3′ end of each primer, where polymerase starts extending. 0 means no requirement; it cannot exceed the minimum primer size."),
        "pcrMM": ("Max mismatches", "Mismatches allowed between a primer and its binding site, counted over the annealed 3′ core. Unmatched bases 5′ of that core are reported as a tail, not counted here."),
        "pcrTP": ("3′ exact bases", "The last bases at a primer's 3′ end must match the template exactly, because polymerase cannot extend a mismatched 3′ end."),
        "pcrPmin": ("Prod min (bp)", "Smallest product the search reports, in base pairs. A designed pair outside the window widens it for the run."),
        "pcrPmax": ("Prod max (bp)", "Largest product the search reports, in base pairs. A designed pair outside the window widens it for the run."),
    }

    def _decorate_primer_fields(self):
        """Give every primer and in-silico PCR field its unit and a plain-language tooltip."""
        for fid, (label, tip) in self._FIELD_HELP.items():
            ed = self.pfields.get(fid)
            if ed is None:
                continue
            ed.setToolTip(tip); ed.setAccessibleName(label); ed.setAccessibleDescription(tip)
            cell = ed.parentWidget()
            for lab in (cell.findChildren(QLabel) if cell is not None else []):
                if lab.buddy() is ed:
                    lab.setText(label.upper()); lab.setToolTip(tip)   # tests/test_41_ui_M1c.py pins the all-caps unit spelling here

    def _show_preset(self, name):
        """Name the preset the fields now hold ('custom' once any of them is edited), without re-applying it."""
        self.pPreset.blockSignals(True)
        self.pPreset.setCurrentIndex(self._preset_keys.index(name))
        self.pPreset.blockSignals(False)

    def _reset_primer_params(self):
        self._apply_preset("standard")

    def _primer_number(self, fid, label, *, integer=False, minimum=None, maximum=None, unit="", why=""):
        raw = self.pfields[fid].text().strip()
        try:
            value = float(raw)
        except ValueError:
            raise ValueError(f"{label} must be a number")
        if not math.isfinite(value):
            raise ValueError(f"{label} must be a finite number")
        if integer and not value.is_integer():
            raise ValueError(f"{label} must be a whole number")
        value = int(value) if integer else value
        if minimum is not None and value < minimum:
            raise ValueError(f"{label} must be at least {minimum:g}{unit}{why}")
        if maximum is not None and value > maximum:
            raise ValueError(f"{label} must be at most {maximum:g}{unit}{why}")
        return value

    # Guard rails on the primer form. The first three are Primer3's / the designer's own rules, stated here so a
    # wrong value is caught in plain words before a request is sent; the last two only keep the form sane.
    _PRIMER_LEN_MAX = 36                 # Primer3's built-in maximum primer length
    _PAIRS_MAX = 50                      # pairs per design (each pair is also folded and tabulated)
    _TM_RANGE = (30, 90)                 # °C; outside this a 18-36 nt primer is not a PCR primer

    def _read_primer_params(self):
        nt, c, lo, hi = " nt", " °C", *self._TM_RANGE
        p = {
            "prod_min": self._primer_number("pMin", "Minimum product size", integer=True, minimum=1, unit=" bp"),
            "prod_max": self._primer_number("pMax", "Maximum product size", integer=True, minimum=1, unit=" bp"),
            "opt_tm": self._primer_number("pTm", "Optimal Tm", minimum=lo, maximum=hi, unit=c),
            "min_size": self._primer_number("pMinS", "Minimum primer size", integer=True,
                                             minimum=primers_mod.MIN_PRIMER_SIZE, maximum=self._PRIMER_LEN_MAX, unit=nt),
            "max_size": self._primer_number("pMaxS", "Maximum primer size", integer=True,
                                             minimum=primers_mod.MIN_PRIMER_SIZE, maximum=self._PRIMER_LEN_MAX, unit=nt,
                                             why=" (Primer3's limit)" ),
            "num_return": self._primer_number("pNum", "Number of pairs", integer=True, minimum=1,
                                              maximum=self._PAIRS_MAX, why=" per design (more are slow and rarely useful)"),
            "opt_size": self._primer_number("pOptS", "Optimal primer size", integer=True,
                                             minimum=primers_mod.MIN_PRIMER_SIZE, maximum=self._PRIMER_LEN_MAX, unit=nt),
            "min_tm": self._primer_number("pTmMin", "Minimum Tm", minimum=lo, maximum=hi, unit=c),
            "max_tm": self._primer_number("pTmMax", "Maximum Tm", minimum=lo, maximum=hi, unit=c),
            "min_gc": self._primer_number("pGcMin", "Minimum GC", minimum=0, maximum=100, unit=" %"),
            "max_gc": self._primer_number("pGcMax", "Maximum GC", minimum=0, maximum=100, unit=" %"),
            "max_poly_x": self._primer_number("pPolyX", "Maximum homopolymer", integer=True, minimum=1, unit=" bases"),
            "gc_clamp": self._primer_number("pGcClamp", "GC clamp", integer=True, minimum=0, unit=" bases"),
        }
        if p["prod_min"] > p["prod_max"]:
            raise ValueError("Minimum product size cannot exceed maximum product size")
        if p["prod_min"] < 2 * p["max_size"]:
            raise ValueError(f"Minimum product size ({p['prod_min']} bp) must be at least twice the maximum primer size "
                             f"({2 * p['max_size']} bp), otherwise the two primers could overlap")
        if not p["min_size"] <= p["opt_size"] <= p["max_size"]:
            raise ValueError("Primer sizes must satisfy minimum ≤ optimal ≤ maximum")
        if not p["min_tm"] <= p["opt_tm"] <= p["max_tm"]:
            raise ValueError("Tm values must satisfy minimum ≤ optimal ≤ maximum")
        if p["min_gc"] > p["max_gc"]:
            raise ValueError("Minimum GC cannot exceed maximum GC")
        if p["gc_clamp"] > p["min_size"]:
            raise ValueError(f"GC clamp ({p['gc_clamp']} bases) cannot exceed the minimum primer size ({p['min_size']} nt)")
        return p

    # =================== primer design ===================
    def _invalidate_primer_results(self):
        self.state["lastPrimers"] = None
        self.state["candidates"] = []
        self.pcrStageAll.setEnabled(False)
        self._refresh_scan_picker()

    # ---------- assay goal (panel 04): which region the one Design button designs on ----------
    def _build_assay_goal(self, card):
        """Assay-goal chooser above the preset. Every goal is the same Primer3 call on the analysed record with a different
        region argument, so the preset, the seal and the in-silico PCR binding are unchanged."""
        grow = QHBoxLayout(); grow.addWidget(QLabel("Assay goal"))
        self.pGoal = _Combo()
        for text, key in (("Detect or quantify this element", "detect"), ("Genotype this insertion", "genotype"),
                          ("Custom range", "custom")):
            self.pGoal.addItem(text, key)
        self.pGoal.setAccessibleName("Assay goal")
        self.pGoal.setToolTip("What the primers are for. Detect or quantify: primers and product inside one protein domain of "
                              "the element. Genotype: primers across an insertion junction, from the flanking DNA. Custom: "
                              "primers and product inside a range you choose.")
        grow.addWidget(self.pGoal, 1); card.bodylay.addLayout(grow)
        self.goalPanels = {}
        for key in ("detect", "genotype", "custom"):
            w = QWidget(); lay = QVBoxLayout(w); lay.setContentsMargins(0, 0, 0, 0); lay.setSpacing(theme_mod.sp(4))
            self.goalPanels[key] = (w, lay); card.bodylay.addWidget(w)
        _, dl = self.goalPanels["detect"]
        row = QHBoxLayout(); row.addWidget(QLabel("Domain")); self.pGoalDomain = _Combo(); row.addWidget(self.pGoalDomain, 1)
        self.pGoalDomain.setAccessibleName("Protein domain to design primers inside"); dl.addLayout(row)
        self.goalDetectNote = QLabel(); dl.addWidget(self.goalDetectNote)
        _, gl_ = self.goalPanels["genotype"]
        row = QHBoxLayout(); row.addWidget(QLabel("Design across")); self.pGoalJunction = _Combo(); row.addWidget(self.pGoalJunction, 1)
        for text, key in (("5′ junction — upstream flank into the element", "5j"),
                          ("3′ junction — element into the downstream flank", "3j"),
                          ("Whole insertion — flank to flank (short insertions)", "span")):
            self.pGoalJunction.addItem(text, key)
        self.pGoalJunction.setAccessibleName("Insertion junction to design across")
        self.goalFlankBtn = QPushButton("Flanking sequence…"); self.goalFlankBtn.setProperty("sm", True)
        self.goalFlankBtn.setToolTip("Copy or export the sequence on either side of the insertion (opens the flank picker).")
        self.goalFlankBtn.clicked.connect(self._goal_flanks)
        row.addWidget(self.goalFlankBtn); gl_.addLayout(row)
        self.goalGenoInfo = QLabel(); gl_.addWidget(self.goalGenoInfo)
        _, cl_ = self.goalPanels["custom"]
        row = QHBoxLayout()
        self.pGoalFrom, self.pGoalTo = QSpinBox(), QSpinBox()
        for sp, nm, lab in ((self.pGoalFrom, "Range start, 1-based", "From"), (self.pGoalTo, "Range end, 1-based inclusive", "to")):
            sp.setRange(1, 1); sp.setAccessibleName(nm); row.addWidget(QLabel(lab)); row.addWidget(sp, 1)
        cl_.addLayout(row)
        self.goalCustomInfo = QLabel(); cl_.addWidget(self.goalCustomInfo)
        for lab in (self.goalDetectNote, self.goalGenoInfo, self.goalCustomInfo):
            lab.setObjectName("orient"); lab.setWordWrap(True)
        self.pGoal.currentIndexChanged.connect(lambda *_: self._goal_changed())
        for c in (self.pGoalDomain, self.pGoalJunction):
            c.currentIndexChanged.connect(lambda *_: self._goal_changed())
        for sp in (self.pGoalFrom, self.pGoalTo):
            sp.valueChanged.connect(lambda *_: self._goal_changed())
        self._goal_rec = None
        for w in (self.pGoal, self.pGoalDomain, self.pGoalJunction, self.pGoalFrom, self.pGoalTo, self.goalFlankBtn):
            w.setEnabled(False)
        self._goal_changed()

    def _goal_ready(self):
        """(ok, why): whether the selected goal can be designed on the selected record."""
        rec = self._goal_rec
        if rec is None:
            return False, "run analysis first"
        goal = self.pGoal.currentData()
        if goal == "detect" and self.pGoalDomain.currentData() is None:
            return False, "no domain to target"
        if goal == "genotype":
            model = self.pGoalJunction.model()
            item = model.item(self.pGoalJunction.currentIndex()) if self.pGoalJunction.currentIndex() >= 0 else None
            if item is None or not item.isEnabled():
                return False, "no flank to place a primer in"
        if goal == "custom" and self.pGoalFrom.value() > self.pGoalTo.value():
            return False, "the range start is after its end"
        return True, ""

    def _design_ready(self):
        return self._goal_rec is not None and self._goal_ready()[0]

    def _goal_changed(self):
        goal = self.pGoal.currentData()
        for key, (w, _l) in self.goalPanels.items():
            w.setVisible(key == goal)
        rec = self._goal_rec
        if rec is not None:
            i, j = self.pGoalFrom.value() - 1, self.pGoalTo.value()
            self.goalCustomInfo.setText(
                f"0-based [{i}, {j}) · {max(0, j - i):,} bp — primers and product stay inside this range; the whole record gives "
                "the untargeted design.")
        ok, why = self._goal_ready()
        if hasattr(self, "designBtn") and rec is not None and not self._design_inflight:
            self.designBtn.setEnabled(ok)
        if hasattr(self, "designHint") and rec is not None:
            b = self._record_binding_text()
            self.designHint.setText(f"template: {b}" if ok and b else why if not ok else self.designHint.text())

    def _refresh_assay_goal(self):
        """Fill the goal panels from the selected record (domains, insertion bounds, range limits) and pick the default goal."""
        rec = self.state.get("last_rec")
        usable = bool(rec) and rec.get("valid", True) and bool(rec.get("seq"))
        widgets_ = (self.pGoal, self.pGoalDomain, self.pGoalJunction, self.pGoalFrom, self.pGoalTo, self.goalFlankBtn)
        for w in widgets_:
            w.setEnabled(usable)
        new = (rec if usable else None) is not self._goal_rec
        prev = self.pGoalDomain.currentData()
        self._goal_rec = rec if usable else None
        for w in widgets_:
            w.blockSignals(True)
        try:
            self.pGoalDomain.clear()
            if not usable:
                for lab in (self.goalDetectNote, self.goalGenoInfo, self.goalCustomInfo):
                    lab.setText("")
                self._goal_changed()
                return
            n = (rec.get("composition") or {}).get("length") or len(rec["seq"])
            choices = coords.domain_choices(rec)
            for i, d in choices:
                s, e = d["nt"]
                self.pGoalDomain.addItem(f"{d['domain']} · {d.get('label', '')} · nt {self._span(s, e)} ({e - s:,} bp)", i)
            if not new and prev is not None and self.pGoalDomain.findData(prev) >= 0:
                self.pGoalDomain.setCurrentIndex(self.pGoalDomain.findData(prev))
            q = self.PRESETS["qpcr"]
            self.goalDetectNote.setText(
                "Primers and product stay inside this domain" + (", clipped to the interior between the terminal repeats, so they "
                "avoid both copies" if coords.interior(rec) else "") + ". No Primer3 target is set: "
                "in-silico PCR checks the intended product but leaves target coverage unassessed. For qPCR choose the qPCR "
                f"preset ({q['pMin']}–{q['pMax']} bp)." if choices else
                "No protein domain was detected in this record — use Custom range.")
            b = coords.insertion_bounds(rec)
            k = primers_mod.MIN_PRIMER_SIZE
            item = self.pGoalJunction.model().item
            for idx, ok_ in enumerate((bool(b) and b["up"] >= k, bool(b) and b["down"] >= k, bool(b) and b["up"] >= k and b["down"] >= k)):
                item(idx).setEnabled(ok_)
            if new or not item(self.pGoalJunction.currentIndex()).isEnabled():
                self.pGoalJunction.setCurrentIndex(next((x for x in range(3) if item(x).isEnabled()), 0))
            if not b:
                self.goalGenoInfo.setText("No LTR or TIR pair bounds this element, so its insertion edges are unknown and "
                                          "genotyping primers cannot be placed here.")
            elif b["up"] < k and b["down"] < k:
                self.goalGenoInfo.setText(
                    f"The record ends at the element ({b['kind']} pair spans {self._span(b['start'], b['end'])} of {n:,} bp), so "
                    "there is no genomic flank to place a primer in. Fetch the insertion by coordinate with extra sequence on "
                    f"each side — panel {PANEL['specimen']} › Fetch by coordinate, for example 1,000 bp more on each side — then "
                    "run analysis again.")
            else:
                self.goalGenoInfo.setText(
                    f"Insertion {self._span(b['start'], b['end'])} (from the {b['kind']} pair) · upstream flank {b['up']:,} bp · "
                    f"downstream flank {b['down']:,} bp. A junction pair amplifies only from the insertion allele; a flank-to-flank "
                    "pair spans the whole insertion, and without the insertion the same pair gives a shorter product.")
            self.goalFlankBtn.setEnabled(bool(b))
            for sp in (self.pGoalFrom, self.pGoalTo):
                sp.setRange(1, max(1, n))
            if new:
                self.pGoalFrom.setValue(1); self.pGoalTo.setValue(max(1, n))
                self.pGoal.setCurrentIndex(0 if choices else 2)
        finally:
            for w in widgets_:
                w.blockSignals(False)
        self._goal_changed()

    def _goto_assay(self, goal):
        """Open panel 04 with an assay goal selected ('detect', 'genotype' or 'custom')."""
        idx = self.pGoal.findData(goal)
        if idx >= 0 and self.pGoal.isEnabled():
            self.pGoal.setCurrentIndex(idx)
        self._scroll_to(self.card_primer)

    def _goal_flanks(self):
        b = coords.insertion_bounds(self.state.get("last_rec") or {})
        if b:
            self._flank_picker(b["start"], b["end"], "+", f"{b['kind']}_element")

    def _design_for_goal(self, _=False):
        """Run the selected assay goal. Re-checks feasibility, because a finished design or a failed run re-enables the button."""
        ok, why = self._goal_ready()
        if not ok:
            return self._banner(f"This assay goal cannot be designed here: {why}.", "info")
        rec, goal = self._goal_rec, self.pGoal.currentData()
        n = (rec.get("composition") or {}).get("length") or len(rec["seq"])
        if goal == "custom":
            i, j = self.pGoalFrom.value() - 1, self.pGoalTo.value()
            if i == 0 and j == n:
                return self._design()
            return self._design_for_domain(i, j, f"range {i + 1}–{j}", within=True)
        if goal == "genotype":
            b = coords.insertion_bounds(rec)
            which = self.pGoalJunction.currentData()
            if which == "5j":
                return self._design_for_domain(b["start"] - 1, b["start"] + 1, "5′ insertion junction")
            if which == "3j":
                return self._design_for_domain(b["end"] - 1, b["end"] + 1, "3′ insertion junction")
            return self._design_for_domain(b["start"], b["end"], "whole insertion")
        d = rec["domains"][self.pGoalDomain.currentData()]
        s, e = d["nt"]
        inner = coords.interior(rec)
        if inner:
            s, e = max(s, inner[0]), min(e, inner[1])
        if e <= s:
            return self._banner(f"The {d['domain']} domain lies inside a terminal repeat, so primers inside it would also bind "
                                "the other copy. Pick another domain, or use Custom range.", "warn")
        self._design_for_domain(s, e, f"{d['domain']} domain", within=True)

    def _stamp_design(self, label=None):
        """At submit time, remember what to call this design's oligos and which template its positions index. `label` names
        an explicit sub-sequence (flank, sub-region, amplicon); None means the selected record itself. A record-less call
        leaves the prefix empty rather than naming the placeholder record."""
        rec = self.state.get("last_rec")
        fid = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(label)).strip("_") if label else ""
        self._design_label = label
        self._design_name_prefix = "_".join(p for p in ((_record_public_id(rec) if rec else ""), fid) if p) or None

    def _oligo(self, c, side):
        """Name of one primer of pair `c`, '<record id>_<pair>_<F|R>': the same wherever it is copied or exported."""
        return coords.oligo_name(c.get("_name_prefix"), self._pair_token(self._pair_label(c)), side)

    def _order_buttons(self, entries_fn, result=False):
        """('Copy order rows', 'Export order sheet') for the pairs `entries_fn()` lists as [(label, candidate)]: one row per
        unique oligo. `result` ties them to the analysis result, so an invalidation disables them with the table."""
        cp, ex = QPushButton("Copy order rows"), QPushButton("Export order sheet")
        for b in (cp, ex):
            b.setProperty("sm", True)
            if result:
                b.setProperty("resultExport", True)
        cp.setToolTip("Copy the name and 5′→3′ sequence of each unique oligo, one per line, tab-separated, for a vendor's "
                      "bulk-entry box. Set scale and purification in the vendor form.")
        ex.setToolTip("Save one row per unique oligo: Name, Sequence, Scale, Purification, Length, Tm, GC, Used in. Scale and "
                      "Purification are left blank for you to fill before ordering.")

        def stem():
            ents = entries_fn()
            return _safe_export_stem((ents[0][1].get("_name_prefix") if ents else None) or "primers", "primers", trust="generated")

        def copy(_=False):
            ents = entries_fn()
            if ents:
                self._copy(coords.order_sheet(ents)[2])

        def export(_=False):
            ents = entries_fn()
            if not ents:
                return
            headers, rows, _t = coords.order_sheet(ents)
            fmt = widgets.pick_table_format(ex, ex.mapToGlobal(ex.rect().bottomLeft()))
            if fmt:
                widgets.export_table(headers, rows, f"TEagle_{stem()}_order_sheet", self, fmt=fmt, notes=None)
        cp.clicked.connect(copy); ex.clicked.connect(export)
        return cp, ex

    def _design_block(self, seq=None):
        """ONE guard for BOTH design entry points: one design at a time, a non-empty template, and — when the
        template IS the panel-01 specimen (seq is None) — a non-stale one. The stale test is deliberately not
        applied to an explicit template (amplicon / family hit / splice product): those do not come from the
        specimen box, so the box's edit state says nothing about them and gating on it would block a working path."""
        if self._design_inflight:
            self._banner("A primer design is already running — wait for it to finish.", "info"); return True
        if seq is None and self._stale_block():
            return True
        if not self._clean_seq(seq if seq is not None else self.state.get("seq", "")):
            self._banner("No sequence to design a primer on.", "info"); return True
        return False

    def _design(self):
        self._clear_banner()
        if self._design_block():
            return
        self._invalidate_primer_results()
        try:
            params = self._read_primer_params()
        except ValueError as e:
            self._set_body(self.primBody, _note(f"Primer design not run — {e}.", "warn"))
            return self._banner(str(e))
        self._design_inflight = True
        self._design_ctx = None                            # an untargeted design: no routed domain to name in the result
        self.designBtn.setEnabled(False); widgets.set_caps(self.designBtn, "◴ designing…")
        self._reset_genome_scan("Primers redesigned — the previous whole-genome scan described the earlier "
                                "pair. Re-scan the pair you want to check.")
        operation_sequence, operation_source = self._operation_input(self.state["seq"])
        self._design_tmpl = operation_sequence             # remember the exact bound template submitted for PCR on-target
        self._stamp_design()
        self.engine.submit("primers", {"sequence": operation_sequence, "params": params,
                                        "source": operation_source},
                           key=self._issue("primers", operation_sequence))

    def _is_selected_record(self, seq):
        """True when `seq` is exactly the selected record of the current analysis."""
        rec = self.state.get("selected_record") or {}
        clean = self._norm_seq(seq)
        return bool(rec.get("sha256")) and rec.get("length") == len(clean) and rec.get("sha256") == sequtil.sha256(clean)

    def _design_for_domain(self, start, end, label, seq=None, *, within=False):
        """Design on a region of the selected record (or an explicit `seq`). Routed (default): the region is Primer3's
        required target, so the product spans it. `within`: the region is Primer3's included region instead, so primers and
        product both lie inside it and no target is set."""
        self._clear_banner()
        whole_record = seq is not None and self._is_selected_record(seq)
        full_explicit_template = seq is not None and start == 0 and end == len(self._norm_seq(seq))
        if whole_record:
            seq = None                                        # the analysed record itself: bind it like panel 04 does (same ID, assessed)
        if self._design_block(seq):
            return
        tmpl = seq if seq is not None else self.state.get("seq", "")
        clean_tmpl = self._norm_seq(tmpl)
        self._invalidate_primer_results()
        if (isinstance(start, bool) or isinstance(end, bool) or not isinstance(start, int) or
                not isinstance(end, int) or start < 0 or end <= start or end > len(clean_tmpl)):
            return self._banner("The selected domain has invalid coordinates; primers were not designed.")
        target = None if (full_explicit_template or within) else [start, end - start]
        included = [start, end - start] if within else (None if full_explicit_template else [0, len(clean_tmpl)])
        try:
            params = self._read_primer_params()
        except ValueError as e:
            self._set_body(self.primBody, _note(f"Primer design not run — {e}.", "warn"))
            return self._banner(str(e))
        length = end - start
        if within and length < params["prod_min"]:
            msg = (f"The {label} is {length:,} bp, shorter than Minimum product size (Prod min, {params['prod_min']:,} bp), so no "
                   "product can lie inside it. Lower Prod min (the qPCR preset allows "
                   f"{self.PRESETS['qpcr']['pMin']}), choose a longer range, or design across a junction instead.")
            self._set_body(self.primBody, _note("Primer design not run — " + msg, "warn"))
            return self._banner(msg, "warn")
        if target is not None and length + 2 * params["min_size"] > params["prod_max"]:
            # Primer3 places one primer on each side of a required target, so the product is at least the target
            # plus two primers. Below that no pair can exist, whatever the other settings: say so and what to change.
            need = length + 2 * params["min_size"]
            msg = (f"The {label} is {length:,} bp long. Primers are placed on each side of it, so the product must be at "
                   f"least {length:,} + 2 × {params['min_size']} = {need:,} bp, but Maximum product size (Prod max) is "
                   f"{params['prod_max']:,} bp. Raise Prod max to at least {need:,} (the Permissive preset allows 1,000), "
                   "or pick a shorter stretch of the domain: right-click it → Select a sub-region.")
            self._set_body(self.primBody, _note("Primer design not run — " + msg, "warn"))
            return self._banner(msg, "warn")
        self._design_ctx = {"label": label, "length": length, "routed": target is not None, "within": within,
                            "prod_max": params["prod_max"], "min_size": params["min_size"]}
        self._scroll_to(self.card_primer)                 # bring the primer designer up so the user lands on the next phase
        self.designBtn.setEnabled(False); widgets.set_caps(self.designBtn, "◴ designing…")   # same busy cue as the toolbar design path
        action = "inside" if within else "on" if target is None else "that flank and span"
        self._set_body(self.primBody, BusyBar(f"Designing primers {action} {label}…"))
        self._pending_domain = label
        self._design_inflight = True
        # Same invalidation as the toolbar design path. Deliberately at REQUEST time, not in _render_primers:
        # that is also called by the theme toggle, so resetting there would wipe a valid scan on every
        # light/dark switch. Over-eager here costs a re-scan prompt; under-eager leaves a stale verdict.
        self._reset_genome_scan("Primers redesigned — the previous whole-genome scan described the earlier "
                                "pair. Re-scan the pair you want to check.")
        if seq is None:
            operation_sequence, operation_source = self._operation_input(tmpl)
        else:
            operation_sequence, operation_source = tmpl, None
        self._design_tmpl = operation_sequence           # candidates' left/right_pos index THIS template, not always panel-01
        self._stamp_design(label if seq is not None else None)
        self.engine.submit("primers", {"sequence": operation_sequence, "params": params,
                                        "primer3_target": target,
                                        "included_region": included,
                                        "source": operation_source},
                           key=self._issue("primers", operation_sequence, bound=seq is None))

    def _on_primers(self, d):
        self._design_inflight = False
        self.designBtn.setEnabled(self._design_ready()); widgets.set_caps(self.designBtn, "Design primers")
        self.card_primer.expand()
        cands = d.get("candidates") or []
        if cands:                                         # one tag per design, so P1 of two designs stay two pairs
            self._batch_n = getattr(self, "_batch_n", 0) + 1
            for c in cands:
                c["_batch"] = f"D{self._batch_n}"
        d["_routed"] = d.get("_routed") or getattr(self, "_design_ctx", None)
        self._design_ctx = None
        self._render_primers(d)
        self._uppercase_buttons()
        if d.get("provenance"):
            self._render_provenance(d["provenance"])
        if getattr(self, "_pending_domain", None):        # a routed 'design here' — re-anchor now the card has its results
            self._scroll_to(self.card_primer)
            self._pending_domain = None

    def _render_primers(self, d):
        _clear_layout(self.primBody)                          # recursive: also removes the addLayout'd "Primer pairs" header row
        self.state["lastPrimers"] = d                         # keep for a theme toggle -> re-render with the new palette's flag colours
        cands = d.get("candidates", [])
        design_tmpl = getattr(self, "_design_tmpl", self.state.get("seq", ""))
        tmpl_sig = self._norm_seq(design_tmpl)
        for c in cands:                                       # tag each pair with the template its coords index
            c["_tmpl_sig"] = tmpl_sig
            c["_design_sequence"] = design_tmpl
            c.setdefault("_design_label", getattr(self, "_design_label", None))      # stamped once: a re-render after a record
            c.setdefault("_name_prefix", getattr(self, "_design_name_prefix", None))  # switch must not rename a design's oligos
        self.state["candidates"] = cands
        self._refresh_scan_picker()                           # keep card-06's in-card scan picker in sync with designed pairs
        self.pcrStageAll.setEnabled(bool(cands))
        routed = d.get("_routed") or {}
        if not cands:
            self.primBody.addWidget(_note(self._no_pair_text(d), "warn"))
            return
        srow = QHBoxLayout(); srow.addWidget(_sl(f"Primer pairs — {len(cands):,}")); srow.addStretch(1)
        srow.addWidget(self._src_link("Primer3")); self.primBody.addLayout(srow)
        recs = self.state.get("records") or []
        design_rec = next((i for i, r in enumerate(recs) if tmpl_sig and self._norm_seq(r.get("seq") or "") == tmpl_sig), None)
        if len(recs) > 1 and design_rec is not None:          # a derived template (flank, sub-region, amplicon) matches no record and gets no line
            sel = self.state.get("record_index")
            b_design = f"record {design_rec + 1} of {len(recs)} · {_record_public_id(recs[design_rec])}"
            if design_rec == sel:
                lab = QLabel(f"Designed on {b_design}."); lab.setObjectName("cardmeta"); self.primBody.addWidget(lab)
            elif self._record_binding_text():
                self.primBody.addWidget(_note(
                    f"These pairs were designed on {b_design}. The record selected in panel 02 is now "
                    f"{self._record_binding_text()}: Design primers will use the selected record, and in-silico PCR withholds "
                    "on-target geometry for pairs designed on another record.", "warn"))
        if routed.get("routed") and routed.get("length"):
            self.primBody.addWidget(_note(
                f"Designed to flank and span the {routed['label']} ({routed['length']:,} bp): each product contains the "
                "domain, and the primers sit in the sequence on either side of it, not inside it. Fwd pos / Rev pos "
                "give where each primer lies.", "info"))
        if routed.get("within") and routed.get("length"):
            self.primBody.addWidget(_note(
                f"Designed inside the {routed['label']} ({routed['length']:,} bp): primers and product both lie within it. No "
                "Primer3 target is set, so in-silico PCR checks the intended product but leaves target coverage unassessed.", "info"))
        cautions = [str(x) for x in d.get("cautions") or []]
        if cautions:
            self.primBody.addWidget(_note(" ".join(cautions), "warn"))
        tm_sentence, dg_sentence = self._conditions(d)
        rec = self.state.get("last_rec") or {}
        bound = bool(rec) and bool(tmpl_sig) and tmpl_sig == self._norm_seq(rec.get("seq") or "")   # coordinates index the record shown
        notes = coords.pair_notes(cands, rec, bound, _repeat_word(rec) if rec else None)
        b0 = self.coord_base
        headers = [self._ch(h) if h in ("Fwd pos", "Rev pos") else h for h in PRIMER_COLS]
        gl = {**self._gloss(),
              "ID": "The design batch (D1, D2, … in the order designs were run) and Primer3's pair number. This is "
                    "the name the pair keeps in panel 05, on the gel, in the genome scan and in every export.",
              "Location": coords.LOCATION_GLOSS + (" Numbers are 1-based, inclusive on screen; exports are 0-based, half-open."
                                                  if b0 else " Numbers are 0-based, half-open [start, end)."),
              "Notes": coords.NOTES_GLOSS,
              "Fwd pos (0-based)": "Where the forward primer lies on the designed sequence, 0-based half-open [start, end).",
              "Rev pos (0-based)": "Where the reverse primer lies on the designed sequence (on the forward strand), 0-based half-open [start, end).",
              "Fwd pos (1-based)": "Where the forward primer lies on the designed sequence, 1-based, inclusive [start, end]. Display only.",
              "Rev pos (1-based)": "Where the reverse primer lies on the designed sequence (on the forward strand), 1-based, inclusive. Display only.",
              "Tm F/R": GLOSS["Tm F/R"] + " " + tm_sentence}
        fc = theme_mod.FLAG[self.theme]                       # per-theme flag colours (WCAG-tuned dark vs light)

        def cells(c, n, base):
            """{column: text, colour, tooltip, export} for one pair, projected over PRIMER_COLS below so the four lists
            cannot fall out of step with the header order."""
            pos = lambda p: self._span(*p, base) if p else "—"
            r = {"ID": self._pair_label(c), "Penalty": c["penalty"], "Notes": notes[n][0],
                 "Location": coords.pair_location(c, rec, bound, base, c.get("_design_label")),
                 "Forward (5'→3')": c["left_seq"], "Reverse (5'→3')": c["right_seq"], "Product": c["product_size"],
                 "Fwd pos": pos(c.get("left_pos")), "Rev pos": pos(c.get("right_pos")),
                 "Tm F/R": f"{c['left_tm']}/{c['right_tm']}", "GC% F/R": f"{c['left_gc']}/{c['right_gc']}"}
            st, tp, ex = {}, {"Notes": notes[n][1] or None}, {}
            qc = c.get("qc") or {}
            if qc.get("ok"):
                L, R = qc["left"], qc["right"]
                hp = _metric_cell([("F", L["hairpin"]), ("R", R["hairpin"])])
                sd = _metric_cell([("F", L["self_dimer"]), ("R", R["self_dimer"])])
                het = _metric_cell([("", qc["hetero_dimer"])])
                end = _metric_cell([("", qc["end_stability"])])
                worst = qc.get("worst", "ok")
                r.update({"Struct": worst, "Hairpin": hp[0], "Self-dim": sd[0], "Hetero": het[0], "3′-end": end[0]})   # Struct = the flag as TEXT
                st = {"Hairpin": fc.get(hp[1]), "Self-dim": fc.get(sd[1]), "Hetero": fc.get(het[1]), "3′-end": fc.get(end[1]),
                      "Struct": fc.get(worst)}
                tp.update({"Hairpin": hp[2], "Self-dim": sd[2], "Hetero": het[2], "3′-end": end[2],
                           "Struct": f"worst secondary-structure flag: {worst}"})
                ex = {"Hairpin": hp[3], "Self-dim": sd[3], "Hetero": het[3], "3′-end": end[3]}   # ΔG exports bare, marks are display-only
            else:                                             # QC unavailable for this pair (never drops the pair)
                r.update({"Struct": "n/a", "Hairpin": "—", "Self-dim": "—", "Hetero": "—", "3′-end": "—"})
            return r, st, tp, ex

        rows, styles, tips, exports, canon = [], [], [], [], []
        for n, c in enumerate(cands):
            r, st, tp, ex = cells(c, n, b0)
            rows.append([r[h] for h in PRIMER_COLS]); styles.append([st.get(h) for h in PRIMER_COLS])
            tips.append([tp.get(h) for h in PRIMER_COLS]); exports.append([ex.get(h) for h in PRIMER_COLS])
            if b0:                                            # the file keeps 0-based half-open positions and the bare ΔG numbers
                r0, _s, _t, ex0 = cells(c, n, 0)
                canon.append([ex0.get(h) if ex0.get(h) is not None else r0[h] for h in PRIMER_COLS])
        t = self._coord_table(headers, gl, canon)
        t.export_header_map = {"Location": "Location (0-based, half-open)",
                               **{f"{k} ({b}-based)": f"{k} (0-based, half-open)" for k in ("Fwd pos", "Rev pos") for b in (0, 1)}}
        t.set_rows(rows, styles=styles, tips=tips, exports=exports)
        caution_cl, warn_cl = (s.replace("3′-end", "3-prime end").replace("ΔG", "delta-G") for s in _dg_clauses("-", "<="))
        if any("qc" in c for c in cands) and not any((c.get("qc") or {}).get("ok") for c in cands):
            dg_sentence = "ΔG checks unavailable for this design."
        t.export_notes = [
            f"Secondary-structure values are delta-G (kcal/mol); more negative is worse. Warn: {warn_cl}. Caution tiers "
            f"({caution_cl.replace('TEagle heuristic: ', '')}) are a TEagle heuristic, not an IDT value. "
            "Primer3 is independently cross-checked with ViennaRNA when available. Exported delta-G cells "
            "contain the bare numeric value; display-only warning marks are omitted. A delta-G of zero or above "
            "means no stable structure and is written as the word none. These predictions are advisory and are "
            "not a wet-lab guarantee.",
            "A routed domain design requires the pair to flank and span that target; an untargeted general "
            "design leaves target status unassessed. Fwd pos and Rev pos are 0-based, half-open positions on the "
            "designed sequence. Location is the product span (0-based, half-open) followed by where each primer lies on the "
            "analysed record; Notes flags primers inside a terminal repeat, pairs with overlapping sites and shared primers.",
            "Reaction conditions: " + tm_sentence + " " + dg_sentence,
        ] + cautions
        t.set_row_menu(lambda r: [("→ send to in-silico PCR", lambda rr=r: self._add_pcr_pair(cands[rr])),
                                  ("Secondary-structure detail", lambda rr=r: self._structure_detail(cands[rr])),
                                  ("Scan whole genome for off-targets", lambda rr=r: self._scan_genome(cands[rr])),
                                  ("Copy pair FASTA", lambda rr=r: self._copy(
                                      f">{self._oligo(cands[rr], 'F')}\n{cands[rr]['left_seq']}\n"
                                      f">{self._oligo(cands[rr], 'R')}\n{cands[rr]['right_seq']}"))])
        t.setMaximumHeight(round(220 * theme_mod.UI_SCALE))
        pxr = QHBoxLayout()
        for w in self._row_action_bar(t, lambda o: f"Pair {self._pair_label(cands[o])}",
                                      ("Send to in-silico PCR", "→ send to in-silico PCR")):
            pxr.addWidget(w)
        pxr.addStretch(1)
        for b in self._order_buttons(lambda: [(self._pair_label(cands[o]), cands[o])
                                              for o in (t._orig(i) for i in range(t.rowCount()))], True):
            pxr.addWidget(b)                                      # the order sheet follows the on-screen sort
        pxr.addWidget(_export_table_btn(t, "TEagle_primers", self))
        self.primBody.addLayout(pxr)                              # above the wide table (a beside button scrolls out)
        self.primBody.addWidget(t)
        # visible key for the ΔG colour flags + the ‡ marker (colour is never the only signal)
        caution_cl, warn_cl = _dg_clauses()
        legend = QLabel(f'<span style="color:{fc["caution"]}">■</span>&nbsp;<b>!</b>&nbsp;caution — {html.escape(caution_cl)} &nbsp;&nbsp;'
                        f'<span style="color:{fc["warn"]}">■</span>&nbsp;<b>!!</b>&nbsp;warn — {html.escape(warn_cl)} &nbsp;&nbsp;'
                        f'‡&nbsp;the two engines disagree')
        legend.setObjectName("orient"); legend.setTextFormat(Qt.RichText); legend.setWordWrap(True)
        self.primBody.addWidget(legend)
        eng = d.get("oligoqc_engines", {})
        cross = (f"cross-checked with ViennaRNA {eng.get('viennarna_version')}"
                 if eng.get("viennarna") else "ViennaRNA cross-check unavailable — showing Primer3 only")
        cond = QLabel(f"<b>Reaction conditions.</b> {html.escape(tm_sentence)} {html.escape(dg_sentence)} Both are calculated "
                      "values for these conditions; recalculate Tm if your PCR buffer differs.")
        cond.setObjectName("orient"); cond.setTextFormat(Qt.RichText); cond.setWordWrap(True); self.primBody.addWidget(cond)
        hint = QLabel("Select a pair, then Row actions ▾ (also right-click or Shift+F10): send to in-silico PCR, secondary-structure "
                      "detail, a genome scan or Copy pair FASTA. Stage all pairs in panel 05.")
        hint.setObjectName("orient"); hint.setWordWrap(True); self.primBody.addWidget(hint)
        read = QLabel(
            "Structure columns are ΔG (kcal/mol) of the most stable hairpin / self-dimer / cross-dimer and 3′-end anneal — more "
            "negative = worse; the legend above gives the flag thresholds. The Struct column states the worst flag as text. "
            "3′-end is its own axis because it isolates the interaction that blocks polymerase extension. Primer3 "
            f"(SantaLucia 1998) is cross-checked against an independent engine ({cross}); Row actions ▾ → “Secondary-structure "
            "detail” shows both engines side by side. Advisory — not a wet-lab guarantee. An exported table carries the bare ΔG "
            "number; the marks are on-screen cues.\n\n"
            "Location is the product span, then where each primer sits on the record: flank, terminal repeat or internal "
            "(from panel 02), with the protein domain when a primer lies inside one. Notes flags a primer whose 3′ end lies inside "
            "a terminal repeat (it can also bind the other copy and solo copies), pairs whose sites overlap (nearly the same "
            "product) and primers shared by two pairs (order them once). It flags; it never removes a pair.\n\n"
            "A routed domain design uses the domain as Primer3's required target, so pairs flank and span it; a design inside a "
            "domain or range sets no target, so target status is unassessed; the general design has neither. A pair designed "
            "inside a transposable element can anneal to every copy of that element, and this table does not say how many places "
            "a pair amplifies: the whole-genome scan (panel 06) does. Scroll the table sideways — or hide the input panel "
            "(Ctrl+B) — to see every column.")
        read.setObjectName("orient"); read.setWordWrap(True)
        self._disclosure(self.primBody, "How to read this table", read, expanded=False)
        self.primBody.addLayout(self._next_step_row([
            (f"Stage all pairs for in-silico PCR → panel {PANEL['pcr']}", self._pcr_stage_all),
            (f"Scan a pair against a downloaded genome → panel {PANEL['genome']}",
             lambda: self._scroll_to(self.card_genome))]))

    @staticmethod
    def _conditions(d):
        """(Tm sentence, ΔG sentence) read from the design result's own condition fields; a field the result lacks
        is reported as not recorded, never filled in from memory."""
        mix = coords.conditions_mix
        tm, qc = mix(d.get("design_tm_conditions") or {}), mix(d.get("oligoqc_conditions") or {}, temp=True)
        return (f"Tm (Primer3 design, SantaLucia 1998 nearest-neighbour): {tm}." if tm else
                "Tm conditions: not recorded in this result.",
                f"ΔG checks (Primer3 thal and the ViennaRNA cross-check): {qc}." if qc else
                "ΔG-check conditions: not recorded in this result.")

    @staticmethod
    def _no_pair_text(d):
        """Why a design returned no pair, from what the result and the request say: a routed domain too long for the
        product window, bases the template could not offer, and Primer3's own count. The settings are named only
        when they are the plausible cause."""
        r = d.get("_routed") or {}
        masked = (d.get("template_ambiguity") or {}).get("masked_to_n") or 0
        if r.get("routed") and r.get("length"):
            need = r["length"] + 2 * (r.get("min_size") or primers_mod.MIN_PRIMER_SIZE)
            out = [f"No primer pair flanks and spans the {r['label']} ({r['length']:,} bp)."]
            if r.get("prod_max") and need > r["prod_max"]:
                out.append(f"A product containing it must be at least {need:,} bp, but Prod max is {r['prod_max']:,}.")
        elif r.get("within") and r.get("length"):
            out = [f"No primer pair fits inside the {r['label']} ({r['length']:,} bp) under the current settings."]
        else:
            out = ["No primer pair met the criteria."]
        out += [str(x) for x in d.get("cautions") or []]
        if masked and not (d.get("cautions") or []):
            out.append(f"{masked} ambiguous template base(s) were treated as N, and no primer covers them.")
        if d.get("explain_pair"):
            out.append(f"Primer3 reports: {d['explain_pair']}.")
        out.append("Choose a region without ambiguous or repetitive bases, widen the product-size range, loosen the Tm or GC "
                   "limits, or try the Permissive preset." if not masked else
                   "Choose a region free of ambiguous bases, or widen the product-size range.")
        return " ".join(out)

    def _structure_detail(self, c):
        """IDT OligoAnalyzer-style detail: both engines' ΔG for every structure of one primer pair, side by side."""
        qc = c.get("qc") or {}
        pname = self._pair_label(c)
        dlg = QDialog(self); dlg.setWindowTitle(f"Secondary structure — {pname}")
        dlg.resize(round(640 * theme_mod.UI_SCALE), round(460 * theme_mod.UI_SCALE))
        lay = QVBoxLayout(dlg)
        head = QLabel(f"<b>{html.escape(pname)}</b>  ·  F 5′-{c['left_seq']}-3′  ·  R 5′-{c['right_seq']}-3′")
        head.setWordWrap(True); lay.addWidget(head)
        if not qc.get("ok"):
            lay.addWidget(QLabel("Secondary-structure QC is unavailable for this pair (" + str(qc.get("error", "")) + ")."))
        else:
            L, R = qc["left"], qc["right"]
            def dg(v):                                        # a ΔG of zero or above is no favourable structure
                return "—" if v is None else (f"{v:.2f} (none)" if v >= 0 else f"{v:.2f}")
            def row(name, m):
                return [name, dg(m.get("p3")), dg(m.get("vrna")), m.get("flag", "ok"),
                        {"none": "no structure", "single": "one engine"}.get(m.get("agree", "single"), m.get("agree", "single"))]
            data = [row("Hairpin (F)", L["hairpin"]), row("Hairpin (R)", R["hairpin"]),
                    row("Self-dimer (F)", L["self_dimer"]), row("Self-dimer (R)", R["self_dimer"]),
                    row("Hetero-dimer (F×R)", qc["hetero_dimer"]), row("3′-end stability", qc["end_stability"])]
            t = DataTable(["Structure", "Primer3 ΔG", "ViennaRNA ΔG", "Flag", "Engines"], GLOSS)
            fc = theme_mod.FLAG[self.theme]
            styles = [[None, None, None, fc.get(r[3]), (fc.get("warn") if r[4] == "disagree" else None)] for r in data]
            t.set_rows(data, styles=styles)
            t.setMinimumHeight(round(230 * theme_mod.UI_SCALE)); lay.addWidget(t)
            cd = qc.get("conditions", {}); eng = qc.get("engines", {})
            vr = eng.get("viennarna")
            cross = (f"cross-checked against ViennaRNA {vr} (DNA Mathews-2004, independent)" if vr and vr != "unavailable"
                     else "with no second engine: the ViennaRNA cross-check was unavailable, so these are Primer3 only")
            caution_cl, warn_cl = _dg_clauses()
            meta_text = (f"ΔG in kcal/mol at {cd.get('temp_c','?')} °C, {cd.get('mv_conc','?')} mM Na⁺, "
                         f"{cd.get('dv_conc','?')} mM Mg²⁺, {cd.get('dna_conc','?')} nM oligo (IDT-comparable). "
                         f"Primer3 {eng.get('primer3','?')} (SantaLucia 1998, thal) {cross}. A ΔG of zero or above means "
                         "no stable structure. Dimer ΔG is the intermolecular "
                         f"binding energy. More negative = more stable = worse; warn is {warn_cl}; caution is a "
                         f"{caution_cl}. The 3′-end is reported as "
                         "its own axis because it isolates the interaction that blocks polymerase extension. Neither "
                         "reproduces IDT’s mfold/UNAFold numbers exactly — this is an "
                         "independent second opinion; agreement within ~1–2 kcal/mol is the useful signal.")
            t.export_notes = [meta_text]
            meta = QLabel(meta_text)
            meta.setObjectName("orient"); meta.setWordWrap(True); lay.addWidget(meta)
            srow = QHBoxLayout(); srow.addWidget(QLabel("Methods:"))
            for key, label in (("SantaLucia1998", "SantaLucia 1998"),
                               ("Owczarzy2008", "OligoAnalyzer convention"),
                               ("ViennaRNA", "ViennaRNA")):
                srow.addWidget(self._src_link(key, label=label))
            srow.addStretch(1)
            srow.addWidget(_export_table_btn(t, f"TEagle_{self._pair_token(pname)}_secondary_structure", self))
            lay.addLayout(srow)
        row = QHBoxLayout(); row.addStretch(1)
        copy = QPushButton("Copy"); copy.clicked.connect(lambda: self._copy(self._structure_text(c)))
        copy.setAccessibleName("Copy the secondary-structure ΔG values as text")
        close = QPushButton("Close"); close.clicked.connect(dlg.accept)
        close.setAccessibleName("Close the secondary-structure detail window")
        row.addWidget(copy); row.addWidget(close); lay.addLayout(row)
        self._uppercase_buttons()
        dlg.exec()
