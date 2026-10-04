"""Main-window chrome: header, rail toggle, theme, UI scale and coordinate display, help, banners and shutdown."""
from __future__ import annotations

import html
from PySide6.QtCore import QPoint, QTimer, Qt
from PySide6.QtGui import QColor, QFont, QKeySequence, QShortcut
from PySide6.QtWidgets import (
    QApplication, QDialog, QFrame, QHBoxLayout, QLabel, QLineEdit, QMessageBox, QPushButton, QScrollArea,
    QSizePolicy, QVBoxLayout, QWidget)
import widgets
from widgets import DataTable
import coords
from coords import CoordTable
from sample import make_sample
import theme as theme_mod
from teagle_core import domains as domains_mod, structural as structural_mod, __version__ as APP_VERSION
from ui_common import (
    CollapsibleCard, GLOSS, MARK_H, PANEL, UI_SCALES, WORD_H, _MIN_ORF_AA, _clear_layout, _empty,
    _example_accession, _example_match, _kb_links, _mark_pixmap, _word_pixmap)


class ChromeMixin:
    def _build_header(self):
        wrap = QWidget()
        col = QVBoxLayout(wrap); col.setContentsMargins(6, 0, 2, 0); col.setSpacing(0)
        h = QHBoxLayout(); h.setContentsMargins(0, 2, 0, theme_mod.sp(8)); h.setSpacing(theme_mod.sp(10))
        self.railToggle = QPushButton(); self.railToggle.setProperty("sm", True)
        self._rail_collapsed = False
        self._paint_rail_toggle()
        self.railToggle.clicked.connect(self._toggle_rail)
        h.addWidget(self.railToggle)
        self.mark = QLabel()                                  # eagle brand mark; pixmap set per-theme in _apply_theme
        self.mark.setObjectName("mark")
        self.mark.setToolTip("TEagle")
        self.mark.setAccessibleName("TEagle")
        h.addWidget(self.mark)
        self.word = QLabel()                                  # Cascadia Code wordmark; pixmap set per-theme in _apply_theme
        self.word.setObjectName("word")
        h.addWidget(self.word)
        self.ver = QLabel("v" + APP_VERSION); self.ver.setObjectName("ver")
        h.addWidget(self.ver)
        tag = QLabel("TRANSPOSABLE ELEMENTS ASSAY TERMINAL"); tag.setObjectName("tagline")
        tf = tag.font(); tf.setLetterSpacing(QFont.AbsoluteSpacing, 1.5); tag.setFont(tf)
        # Preferred (NOT Ignored) so the label's own sizeHint is honoured — under Ignored the layout gave it 0px
        # at every window width and the tagline never rendered. minimumWidth(1) keeps the intent Ignored had: the
        # decorative label never imposes its full text width as the window's hard minimum, so a small screen at
        # ≥125% UI scale can still shrink the window inside the display (it is hidden below 1080*scale anyway).
        tag.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)
        tag.setMinimumWidth(1)
        self._tagline = tag                                   # hidden on a narrow window (resizeEvent) so it never overhangs
        h.addWidget(tag)
        h.addStretch(1)
        chip = QFrame(); chip.setObjectName("statuschip")
        cl = QHBoxLayout(chip); cl.setContentsMargins(theme_mod.sp(10), theme_mod.sp(5), theme_mod.sp(11), theme_mod.sp(5)); cl.setSpacing(theme_mod.sp(8))
        self.led = QLabel(); self.led.setObjectName("led")   # size via QSS #led min/max-width (scales with UI_SCALE)
        # Two labels, because they describe two different things: the in-process core engine (classification, ORFs,
        # domains, primers, in-silico PCR — always there) and the OPTIONAL Linux backend (family naming, splice
        # structure, whole-genome scan). One "backend live" chip read as the second being available.
        self.statusTxt = QLabel("connecting…"); self.statusTxt.setObjectName("statusTxt")
        self.backendTxt = QLabel("Linux backend: checking…"); self.backendTxt.setObjectName("statusTxt")
        for lab in (self.statusTxt, self.backendTxt):         # Preferred + a 1 px floor: clips (full text in the tooltip) on a narrow window, never widens it
            lab.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred); lab.setMinimumWidth(1)
        dot = QLabel("·"); dot.setObjectName("statusTxt")
        cl.addWidget(self.led); cl.addWidget(self.statusTxt); cl.addWidget(dot); cl.addWidget(self.backendTxt)
        self._chip = chip
        h.addWidget(chip)
        self._core, self._wsl = None, None                    # health / wsl_status payloads, in whichever order they arrive
        # Top-level entry to the backend installer. It used to be reachable only from a secondary button
        # inside panel 03, where a user who had not yet opened that card never saw it — yet it is what
        # turns on Dfam family naming, splice structure and whole-genome scans.
        self.backendBtn = QPushButton("Backend"); self.backendBtn.setProperty("sm", True)
        self.backendBtn.setToolTip("Linux backend — install, repair and check the optional Linux (WSL) stack: Dfam "
                                   "family naming, splice structure, whole-genome scan")
        self.backendBtn.setAccessibleName("Open the Linux backend installer")
        self.backendBtn.clicked.connect(self._open_installer)
        h.addWidget(self.backendBtn)
        self.scaleBtn = QPushButton(); self.scaleBtn.setProperty("sm", True)
        self.scaleBtn.setToolTip("Scale the whole interface — text and spacing. Applied live; window geometry is exact "
                                 "from the next launch.")
        self.scaleBtn.clicked.connect(lambda: self._ui_scale_menu(self.scaleBtn))
        h.addWidget(self.scaleBtn)
        self.coordBtn = QPushButton(); self.coordBtn.setProperty("sm", True)
        self.coordBtn.setToolTip("How coordinates are numbered on screen. Display only: exported tables, FASTA names, GFF3/BED "
                                 "files and primer-design coordinates keep their documented conventions.")
        self.coordBtn.setAccessibleName("Coordinate numbering on screen")
        self.coordBtn.setMinimumWidth(1)                      # clips on a very narrow window instead of raising its floor
        self.coordBtn.clicked.connect(lambda: self._coord_menu(self.coordBtn))
        self._paint_coord_btn()
        h.addWidget(self.coordBtn)
        self.themeBtn = QPushButton(); self.themeBtn.setProperty("sm", True); self.themeBtn.clicked.connect(self._toggle_theme)
        h.addWidget(self.themeBtn)
        self.helpBtn = QPushButton("Help ▾"); self.helpBtn.setProperty("sm", True)
        self.helpBtn.setAccessibleName("Help: quick start, glossary and user manual")
        self.helpBtn.clicked.connect(lambda: self._help_menu(self.helpBtn))
        h.addWidget(self.helpBtn)
        col.addLayout(h)
        self.headrule = QFrame(); self.headrule.setObjectName("headrule"); self.headrule.setFixedHeight(2)
        col.addWidget(self.headrule)
        QShortcut(QKeySequence("Ctrl+B"), self, activated=self._toggle_rail)   # toggle the specimen panel
        QShortcut(QKeySequence(Qt.Key_F1), self, activated=self._show_glossary)
        return wrap

    def _toggle_rail(self):
        """Collapse/reopen the left specimen panel — hidden, the results panel takes the full width, so wide
        tables (primer QC, off-target scan) need far less horizontal scroll. The header button always shows.
        State is tracked explicitly (isVisible() is unreliable before the window is shown)."""
        collapsed = not getattr(self, "_rail_collapsed", False)
        self._rail_collapsed = collapsed
        if collapsed:
            self._rail_sizes = self.split.sizes()             # remember the split so reopening restores it
            self.rail.setVisible(False)
        else:
            self.rail.setVisible(True)
            total = sum(self.split.sizes()) or 1240
            sizes = self._rail_sizes if sum(self._rail_sizes) <= total else [340, max(300, total - 340)]
            self.split.setSizes(sizes)
        self._paint_rail_toggle()

    def _paint_rail_toggle(self):
        """Name what the toggle hides: the whole of panel 01 (specimen input and readouts), not 'the panel'."""
        hide = not self._rail_collapsed
        widgets.set_caps(self.railToggle, "◂ Hide input" if hide else "▸ Show input")
        self.railToggle.setToolTip(("Hide panel 01 (specimen input and readouts) for more result width — Ctrl+B" if hide else
                                    "Show panel 01 (specimen input and readouts) — Ctrl+B"))
        self.railToggle.setAccessibleName("Hide input panel" if hide else "Show input panel")
        self.railToggle.setAccessibleDescription("Shortcut Ctrl+B")

    def _ui_scale_menu(self, anchor):
        """Pick a global UI scale — applied LIVE (font + spacing) with no restart; the active value is tagged
        'current'. Code-set window/rail geometry is pixel-exact from the next launch."""
        from PySide6.QtWidgets import QMenu
        cur = theme_mod.UI_SCALE
        m = QMenu(self)
        for f in UI_SCALES:
            act = m.addAction(f"{int(f * 100)}%" + ("   (current)" if abs(f - cur) < 1e-3 else ""))
            act.triggered.connect(lambda _=False, x=f: self._set_ui_scale(x))
        m.exec(anchor.mapToGlobal(anchor.rect().bottomLeft()))

    # ---------- coordinate display (screen only; stored values, exports and design coordinates stay 0-based half-open) ----------
    def _paint_coord_btn(self):
        """'Coords: 0-based', or just '0-based' on a window too narrow for the prefix (the tooltip names the setting)."""
        full = self.width() >= round(980 * theme_mod.UI_SCALE)
        widgets.set_caps(self.coordBtn, f"Coords: {self.coord_base}-based" if full else f"{self.coord_base}-based")

    def _coord_menu(self, anchor):
        from PySide6.QtWidgets import QMenu
        m = QMenu(self)
        for b, text in coords.BASES.items():
            m.addAction(text + ("   (current)" if b == self.coord_base else "")).triggered.connect(
                lambda _=False, x=b: self._set_coord_base(x))
        m.exec(anchor.mapToGlobal(anchor.rect().bottomLeft()))

    def _span(self, s, e, base=None):
        return coords.span(s, e, self.coord_base if base is None else base)

    def _ch(self, name):
        return coords.col(name, self.coord_base)

    def _gloss(self, extra=None):
        return {**GLOSS, **(extra or {}), **coords.gloss(self.coord_base)}

    def _coord_caption(self):
        if self.coord_base:
            return ("Coordinates: 1-based, inclusive [start, end] (GenBank, GFF3, SnapGene). Display only — exported tables, "
                    "FASTA names and primer-design coordinates stay 0-based, half-open.")
        return ("Coords are 0-based and half-open [start, end) — the convention of TEagle's exported tables. Protein-domain aa "
                "and Dfam consensus ranges are 1-based, inclusive, as HMMER and RepeatMasker report them. Change with the "
                "Coords button in the header.")

    def _coord_table(self, headers, tips, canon=None):
        """A DataTable showing the chosen numbering whose exports stay canonical: `canon` is the rows as they are at base 0
        (by original index), written to a file whatever the screen shows."""
        t = CoordTable(headers, tips)
        if canon is not None and self.coord_base:
            t.canon_rows = canon
        return t

    # ---------- results column ----------
    def _build_results(self):
        wrap = QScrollArea(); wrap.setWidgetResizable(True)
        # Sideways scrolling ON DEMAND, not never. The invariant was "the page body never scrolls sideways;
        # wide tables scroll inside themselves" — which holds for tables, and fails for everything that is
        # not one. At 1.5x UI scale in a window under ~1230 px the figure toolbar (bg modes + zoom + FIT +
        # SVG + PNG) exceeds the viewport, and with the bar ALWAYS off Qt squeezes the inner widget below its
        # own minimum and CLIPS the children that cannot shrink: the SVG and PNG export buttons became
        # unreachable, with no scrollbar to reveal them and no way to get them back short of resizing the
        # window. AsNeeded keeps the intended behaviour whenever the content fits — which is every normal
        # window — and degrades to a scrollbar instead of to lost controls when it does not.
        wrap.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        self.resultsScroll = wrap                          # inside their own viewport (per the crowding invariant)
        inner = QWidget(); inner.setObjectName("central")
        self.results = QVBoxLayout(inner); self.results.setContentsMargins(theme_mod.sp(4), 0, theme_mod.sp(8), 0); self.results.setSpacing(theme_mod.sp(9))

        self.results.addWidget(self._build_quickstart())
        self.card_struct = CollapsibleCard(PANEL["struct"], "Classification & structure",
                                           "LTR/TIR repeats, ORFs, protein domains")
        self.card_struct.bodylay.addWidget(_empty("Run analysis to detect terminal repeats, ORFs and tails."))
        self.results.addWidget(self.card_struct)

        self.results.addWidget(self._build_wsl_card())
        self.results.addWidget(self._build_splice_card())
        self.results.addWidget(self._build_primer_card())
        self.results.addWidget(self._build_pcr_card())
        self.results.addWidget(self._build_genome_card())
        self.card_prov = CollapsibleCard(PANEL["prov"], "Run provenance", "versions + checksums for reproducibility")
        self.card_prov.bodylay.addWidget(_empty("Populated from the last computation — travels with every result."))
        self.results.addWidget(self.card_prov)

        self.results.addStretch(1)
        wrap.setWidget(inner)
        return wrap

    def _build_quickstart(self):
        """A dismissible 'Start here' strip above the cards: the three steps, the bundled examples as one-click runs,
        and the glossary. Shown until the user hides it; Help > Quick start brings it back."""
        from teagle_core import examples
        box = QFrame(); box.setObjectName("quickstart"); self._quickstart = box
        lay = QVBoxLayout(box)
        lay.setContentsMargins(theme_mod.sp(14), theme_mod.sp(10), theme_mod.sp(14), theme_mod.sp(12)); lay.setSpacing(theme_mod.sp(7))
        top = QHBoxLayout()
        h = QLabel("Start here"); h.setObjectName("sech"); top.addWidget(h); top.addStretch(1)
        hide = QPushButton("Hide"); hide.setProperty("link", True); hide.setAccessibleName("Hide the quick-start guide")
        hide.clicked.connect(lambda _=False: self._set_quickstart_hidden(True))
        top.addWidget(hide); lay.addLayout(top)
        steps = QLabel(f"<b>1</b> Load a sequence in panel {PANEL['specimen']} — an accession, a FASTA file, pasted DNA, or an "
                       f"example below. <b>2</b> Run analysis — card {PANEL['struct']} shows the call and the evidence behind it. "
                       f"<b>3</b> Design primers (panel {PANEL['primer']}), then check them by in-silico PCR "
                       f"(panel {PANEL['pcr']}).")
        steps.setObjectName("cardnote"); steps.setWordWrap(True); steps.setTextFormat(Qt.RichText); lay.addWidget(steps)
        row = QHBoxLayout(); row.setSpacing(theme_mod.sp(8))
        k = QLabel("Analyse a published element:"); k.setObjectName("cardmeta"); row.addWidget(k)
        self._quickstart_links = []
        for acc, label, organism, expected in examples.available():
            b = QPushButton(label); b.setProperty("link", True)
            b.setToolTip(f"{organism} · {acc}. Published as: {expected}. Loads it into panel {PANEL['specimen']} and runs the analysis.")
            b.clicked.connect(lambda _=False, a=acc: self._quickstart_example(a))
            row.addWidget(b); self._quickstart_links.append(b)
        neg = QPushButton("negative control"); neg.setProperty("link", True)
        neg.setToolTip("Seeded random LTR-like sequence. Structure is detectable; the ORF is random, so no protein domain "
                       "should be found — that null result is the expected outcome.")
        neg.clicked.connect(lambda _=False: self._quickstart_example(None))
        row.addWidget(neg); self._quickstart_links.append(neg)
        row.addStretch(1); lay.addLayout(row)
        row2 = QHBoxLayout(); row2.setSpacing(theme_mod.sp(14))
        for text, fn in (("Glossary (F1)", self._show_glossary), ("User manual (online)", self._open_manual)):
            b = QPushButton(text); b.setProperty("link", True); b.clicked.connect(lambda _=False, f=fn: f()); row2.addWidget(b)
        row2.addStretch(1); lay.addLayout(row2)
        box.setVisible(not self._quickstart_hidden())
        self._style_quickstart()
        return box

    def _style_quickstart(self):
        p = theme_mod.PAL[self.theme]
        self._quickstart.setStyleSheet(f"QFrame#quickstart {{ background: {p['panel2']}; border: 1px solid {p['line2']}; "
                                       f"border-left: 3px solid {p['accent']}; border-radius: 2px; }}")

    def _quickstart_example(self, accession):
        """One click: load the example and analyse it. Never goes through the modal example menu."""
        if accession is None:
            self._apply_example(make_sample())
        elif not self._load_example(accession):
            return
        self._run_analysis()

    def _help_menu(self, anchor):
        from PySide6.QtWidgets import QMenu
        m = QMenu(self)
        m.addAction("Quick start", self._show_quickstart)
        m.addAction("Glossary…   F1", self._show_glossary)
        a = m.addAction("User manual (online) ↗", self._open_manual)
        a.setToolTip(f"Opens the TEagle {APP_VERSION} release page, where the illustrated user manual (PDF) is attached.")
        m.setToolTipsVisible(True)
        m.exec(anchor.mapToGlobal(anchor.rect().bottomLeft()))

    def _show_quickstart(self):
        self._set_quickstart_hidden(False)
        self.resultsScroll.verticalScrollBar().setValue(0)

    def _glossary_html(self):
        e = html.escape
        entries = [
            ("Call status", "the “Call status” column of the record table. " + GLOSS["Call status"] +
             " Each analysis stage's own status (complete, partial, failed, not_run) is in the Processing column."),
            ("Call confidence (High, Moderate, Candidate)", "how strongly the combination of protein domains and structural "
             "hallmarks supports the order or superfamily call under TEagle's classification rules (Wicker et al. 2007 "
             "scheme). It is a category, not a probability, and it is separate from the completeness tier."),
            ("Hit confidence", GLOSS["Hit confidence"]),
            ("Completeness tier", "how much of the protein-domain architecture expected for the called order was detected in "
             "this sequence (for example intact / autonomous-consistent, near-complete, partial, fragment, structural "
             "candidate). It is scoped to the tested Pfam panel and reported separately from call confidence: a High-confidence "
             "call can carry a partial tier. It does not establish transcription or transposition competence."),
            ("Not detected is not absent", "‘not detected’ means a search ran and found nothing above its threshold in this "
             f"sequence. A diverged or unmodelled domain reads as not detected, and terminal repeats diverged past the detectors' "
             f"identity floors (LTR {structural_mod.MIN_LTR_IDENTITY:g}%, TIR {structural_mod.MIN_TIR_IDENTITY:g}%) are not "
             f"called; sub-threshold LTR pairs are listed only as advisory candidates. The domain search covers only the "
             f"{domains_mod.MAX_ORFS_SCANNED} longest ORFs of at least {_MIN_ORF_AA} aa. ‘Not assessed’ or ‘unavailable’ means "
             "the search did not run, which says nothing about the sequence."),
            ("Coordinates", "positions in TEagle's tables and FASTA headers are 0-based and half-open, [start, end): the first "
             "base is 0, the end is excluded, and length = end − start. Exceptions are labelled where they appear: Fetch by "
             "coordinate takes UCSC-style 1-based inclusive positions; the sub-region picker takes 1-based inclusive offsets "
             "within the feature; the whole-genome scan table shows isPcr's 1-based inclusive coordinates; the aa column of the "
             "domain table is HMMER's 1-based residue range on the translated ORF. GFF3 exports are 1-based inclusive, as the "
             "format requires; BED exports are 0-based half-open."),
            ("Optional Linux backend", f"family naming (panel {PANEL['family']}), splice structure (panel {PANEL['splice']}) and "
             f"the whole-genome scan (panel {PANEL['genome']}) run in an optional Linux (WSL) backend, installed from the "
             "Backend button in the header. Classification, ORFs, protein domains, primer design and in-silico PCR run in "
             "the core engine and need nothing else."),
            ("Bundled examples", "five published elements. After analysing one, card " + PANEL["struct"] +
             " shows its published class beside TEagle's call."),
        ]
        return ("".join(f"<p><b>{e(t)}</b> — {e(d)}</p>" for t, d in entries)
                + "<p>Each table column has its own definition: hover over the column header.</p>")

    def _show_glossary(self):
        dlg = getattr(self, "_glossary_dlg", None)
        if dlg is None:
            dlg = QDialog(self); dlg.setObjectName("glossarydlg"); dlg.setWindowTitle("TEagle glossary")
            dlg.resize(round(640 * theme_mod.UI_SCALE), round(560 * theme_mod.UI_SCALE))
            lay = QVBoxLayout(dlg); lay.setContentsMargins(*(theme_mod.sp(12),) * 4)
            area = QScrollArea(); area.setWidgetResizable(True); lay.addWidget(area)
            body = QLabel(); body.setObjectName("cardnote"); body.setWordWrap(True); body.setTextFormat(Qt.RichText)
            body.setTextInteractionFlags(Qt.TextSelectableByMouse); body.setAlignment(Qt.AlignTop | Qt.AlignLeft)
            body.setMargin(theme_mod.sp(8)); area.setWidget(body)
            dlg._body = body
            self._glossary_dlg = dlg
        dlg._body.setText(self._glossary_html())
        dlg.show(); dlg.raise_(); dlg.activateWindow()

    def _grid_fields(self, grid, specs, cols):
        for idx, (fid, label, default) in enumerate(specs):
            r, c = divmod(idx, cols)
            cell = QWidget(); cl = QVBoxLayout(cell); cl.setContentsMargins(0, 0, 0, 0); cl.setSpacing(1)
            lab = QLabel(widgets.caps(label)); lab.setObjectName("kdim")
            ed = QLineEdit(default); ed.setMaximumWidth(round(90 * theme_mod.UI_SCALE))
            lab.setBuddy(ed); ed.setAccessibleName(label)
            if fid in self.PRESETS["standard"]:                # a hand edit means the fields no longer hold a named preset
                ed.textEdited.connect(lambda _t="": self._show_preset("custom"))
            cl.addWidget(lab); cl.addWidget(ed)
            grid.addWidget(cell, r, c)
            self.pfields[fid] = ed

    # ---------- theme ----------
    def _apply_theme(self):
        QApplication.instance().setStyleSheet(theme_mod.qss(self.theme))
        accent = theme_mod.ACCENT[self.theme]
        dpr = self.devicePixelRatioF(); z = theme_mod.UI_SCALE
        self.mark.setPixmap(_mark_pixmap(accent, round(MARK_H * z), dpr))
        word_te = "#FFFFFF" if self.theme == "dark" else theme_mod.TEXT[self.theme]  # TE white on dark, dark ink on light
        self.word.setPixmap(_word_pixmap(word_te, accent, round(WORD_H * z), dpr))   # AGLE = accent = eagle-mark colour
        self.headrule.setStyleSheet(f"QFrame#headrule {{ background: {theme_mod.HEADRULE[self.theme]}; }}")
        pct = round(theme_mod.UI_SCALE * 100)
        widgets.set_caps(self.scaleBtn, f"Text size {pct}%")
        self.scaleBtn.setAccessibleName(f"Interface scale, currently {pct}%")
        dark = self.theme == "dark"
        widgets.set_caps(self.themeBtn, "Light mode" if dark else "Dark mode")
        self.themeBtn.setToolTip("Switch to the light theme" if dark else "Switch to the dark theme")
        self.themeBtn.setAccessibleName("Switch to light theme" if dark else "Switch to dark theme")
        if getattr(self, "_quickstart", None) is not None:
            self._style_quickstart()
        prev = getattr(self, "_tinted_theme", self.theme)
        if prev != self.theme:
            self._retint(prev, self.theme)
        self._tinted_theme = self.theme
        self._uppercase_buttons()
        self._sync_viewer_themes()

    def _retint(self, prev, new):
        """Rich-text links and table ink carry the colour of the theme they were written in (a link's colour is fixed when Qt
        parses the HTML, and the app stylesheet cannot reach it), so after a toggle the old theme's accent stayed behind: 1.8:1
        on the other panel. Swap each of our own palette colours for its counterpart in every label and table cell."""
        pal = theme_mod.PAL
        swap = {pal[prev][k]: pal[new][k] for k in ("accent", "good", "bad", "warn")}
        for lab in self.findChildren(QLabel):
            text = lab.text()
            fixed = text
            for old, now in swap.items():
                fixed = fixed.replace(f"color:{old}", f"color:{now}")
            if fixed != text:                                 # changed text re-parses; identical text would not
                lab.setText(fixed)
        up = {k.upper(): QColor(v) for k, v in swap.items()}
        for table in self.findChildren(DataTable):
            for r in range(table.rowCount()):
                for c in range(table.columnCount()):
                    it = table.item(r, c)
                    now = it is not None and up.get(it.foreground().color().name().upper())
                    if now:
                        it.setForeground(now)

    def _sync_viewer_themes(self):
        """Push the app theme into every live genome viewer + figure panel (gel), by default, so they follow
        dark/light. A live findChildren walk (no persistent registry -> no dangling C++ refs across the
        constant card rebuilds); each viewer re-renders in place, preserving its pan/zoom. Per-viewer buttons
        still override until the next app-theme toggle."""
        from widgets import GenomePanel, FigurePanel, DotPanel
        for v in self.findChildren(GenomePanel) + self.findChildren(FigurePanel) + self.findChildren(DotPanel):
            v.apply_app_theme(self.theme)

    def _uppercase_buttons(self, root=None):
        """Uppercase action buttons under `root` (default: the whole window, dialogs parented to it included).
        A freshly rendered panel/dialog passes itself so it never sits in sentence case until the next refresh."""
        widgets.uppercase_buttons(root if root is not None else self)

    def _repolish(self, w):
        w.style().unpolish(w); w.style().polish(w)

    def _toggle_env(self):
        vis = not self.envBox.isVisible()
        self.envBox.setVisible(vis)
        self.envHdr.setText(("▾" if vis else "▸") + " Environment")

    # ---------- startup ----------
    def _startup(self):
        self.engine.submit("health")
        self.engine.submit("env")
        self._init_wsl()
        self._render_pcr_queue()

    def _banner(self, msg, level="error"):
        # transient messages surface as a small, closable notification dialog centred over the window — never a
        # banner wedged above the panels. Non-modal (trailing callback logic must keep running); single active
        # instance (replace, never stack); warn/error persist until dismissed, info/success auto-close.
        old = getattr(self, "_notif", None)
        if old is not None:
            try:
                old.close()
            except Exception:
                pass
        dlg = QDialog(self); dlg.setObjectName("notif"); dlg.setProperty("level", level); dlg.setModal(False)
        dlg.setAttribute(Qt.WA_DeleteOnClose, True)        # a dismissed dialog is freed, not left parented+alive for the session
        dlg.setWindowTitle({"error": "Error", "warn": "Warning", "success": "Done", "info": "Notice"}.get(level, "Notice"))
        self._notif = dlg
        dlg.finished.connect(lambda *_, d=dlg: (setattr(self, "_notif", None) if getattr(self, "_notif", None) is d else None))
        lay = QVBoxLayout(dlg)
        lay.setContentsMargins(theme_mod.sp(16), theme_mod.sp(16), theme_mod.sp(16), theme_mod.sp(14)); lay.setSpacing(theme_mod.sp(12))
        lab = QLabel(msg); lab.setObjectName("notifmsg"); lab.setProperty("level", level); lab.setWordWrap(True)
        lab.setTextInteractionFlags(Qt.TextSelectableByMouse)   # accessions/seals/hashes in a message can be selected + copied
        lab.setMinimumWidth(round(300 * theme_mod.UI_SCALE)); lab.setMaximumWidth(round(460 * theme_mod.UI_SCALE))
        lay.addWidget(lab)
        brow = QHBoxLayout(); brow.addStretch(1)
        ok = QPushButton("Close"); ok.setProperty("sm", True); ok.clicked.connect(dlg.close)
        ok.setAccessibleName("Close this notification")
        brow.addWidget(ok); lay.addLayout(brow)
        self._uppercase_buttons(dlg)                      # CLOSE reads the same here as in every other dialog
        self._repolish(dlg); self._repolish(lab)
        dlg.adjustSize()
        try:                                              # centre over the window, then clamp ALL four edges onto the screen
            g = self.frameGeometry(); dr = dlg.frameGeometry()
            av = self.screen().availableGeometry()
            x = min(max(g.center().x() - dr.width() // 2, av.left()), max(av.left(), av.right() - dr.width()))
            y = min(max(g.center().y() - dr.height() // 2, av.top()), max(av.top(), av.bottom() - dr.height()))
            dlg.move(x, y)
        except Exception:
            pass
        if level in ("info", "success"):                  # confirmations float in without stealing focus, then auto-dismiss
            dlg.setAttribute(Qt.WA_ShowWithoutActivating, True)
            t = QTimer(dlg); t.setSingleShot(True); t.timeout.connect(dlg.close); t.start(4500)   # parented to dlg -> dies with it, never fires on a freed object
            dlg.show()
        else:
            dlg.show(); dlg.raise_(); ok.setFocus()

    def _clear_banner(self):
        old = getattr(self, "_notif", None)
        if old is not None:
            try:
                old.close()
            except Exception:
                pass
            self._notif = None

    def _render_env(self, e):
        if e.get("error"):
            self.envBox.setText(f"<span style='color:{theme_mod.BAD[self.theme]}'>{e['error'][:60]}</span>"); return
        pkgs = "<br>".join(f"{p['name']} {'ok' if p.get('ok') else str(p.get('installed','missing'))}"
                           for p in e.get("packages", []))
        st = ("install needed" if e.get("needs_install") else "up to date")
        bw = e.get("backends", {})
        self.envBox.setText(
            f"<b>state</b> {st}{' · first run' if e.get('first_run') else ''}<br>"
            f"<b>python</b> {e.get('python','?')} {'ok' if e.get('python_ok') else 'old'}<br>"
            f"<b>packages</b><br>{pkgs}<br>"
            f"<b>wsl2</b> {bw.get('wsl2','—')}<br>"
            f"<b>signature</b> {e.get('signature','—')}")

    def _accent_link(self, href, text):
        return f"<a href='{href}' style='color:{theme_mod.ACCENT[self.theme]};text-decoration:none'>{text}</a>"

    def _scroll_to(self, widget):
        """Expand a result card and bring it into view (in-app navigation for the traceable rail links)."""
        try:
            widget.expand(); self.resultsScroll.ensureWidgetVisible(widget)
        except Exception:
            pass

    def _scroll_to_top(self, widget):
        """Expand a result card and put its top edge at the top of the results viewport, so the headline of a
        new analysis is what the user sees rather than wherever the previous result left the page."""
        def align():
            try:
                inner = self.resultsScroll.widget()
                self.resultsScroll.verticalScrollBar().setValue(widget.mapTo(inner, QPoint(0, 0)).y())
            except RuntimeError:
                pass
        try:
            widget.expand()
        except RuntimeError:
            return
        align()
        QTimer.singleShot(0, align)                           # again once the expanded card has been laid out

    def _disclosure(self, layout, summary, body, expanded=False):
        """Progressive disclosure: a one-line clickable summary that expands/collapses a deep-detail widget, so the
        card shows only the essential readout by default and the full methodology stays one click away (less crowding)."""
        body.setVisible(expanded)
        tog = QPushButton(("▾ " if expanded else "▸ ") + summary)
        tog.setProperty("link", True); tog.setObjectName("disclose"); tog.setCursor(Qt.PointingHandCursor)
        tog.setAccessibleName(summary)                    # the ▸/▾ glyph is not speakable
        tog.setAccessibleDescription("Disclosure toggle — shows or hides this section's detail.")
        def _toggle():
            vis = not body.isVisible()
            body.setVisible(vis)
            tog.setText(("▾ " if vis else "▸ ") + summary)
        tog.clicked.connect(_toggle)
        layout.addWidget(tog); layout.addWidget(body)
        return tog

    def _set_trace_counts(self, recs, rec):
        """Make the rail readouts traceable, split by provenance class. RECORDS -> the fetched source accession as
        an external DB link (or 'user-supplied' for pasted input, never a fabricated link). STRUCTURAL EVIDENCE and
        ORFS -> in-app links that scroll to their card-02 tables; these are de-novo heuristic detections with no
        external record, so they are NEVER rendered as a database accession/link."""
        src = self.state.get("source") or {}
        acc = src.get("accession") or src.get("displayLocus")
        if acc and src.get("sourceUrl"):
            extra = f"  +{len(recs) - 1:,} more" if len(recs) > 1 else ""
            self.rRecords.setText(self._accent_link(src["sourceUrl"], str(acc)) + extra)
            _kb_links(self.rRecords)                          # text now carries a link -> re-arm the tab stop
        elif acc:
            self.rRecords.setText(str(acc))
        elif len(recs) == 1 and _example_match(rec.get("seq")):
            acc = _example_accession(_example_match(rec.get("seq"))[0])
            self.rRecords.setText("bundled example · GenBank " + self._accent_link(f"https://www.ncbi.nlm.nih.gov/nuccore/{acc}", acc))
            _kb_links(self.rRecords)
            self.rRecords.setToolTip(
                f"Published sequence shipped with TEagle (GenBank {acc}); the sequence is identical to it. Exports and provenance "
                f"use the neutral record ID because this copy was not fetched; fetch {acc.split('.')[0]} with the accession box to "
                "carry the accession into exports and provenance.")
        else:
            self.rRecords.setText(f"user-supplied ({len(recs):,})" if recs else "—")
        self.rStruct.setText(self._accent_link("#struct", f"{len(rec.get('structural', [])):,}"))
        self.rStructCap.setVisible(True)
        from teagle_core import gff3
        tdr = gff3.is_neutral_tdr(rec.get("classification"))
        self.rStruct.setToolTip(" · ".join(("TDR" if tdr and e["type"].startswith("LTR") else e["type"].split(" (")[0])
                                           for e in rec.get("structural", [])) or "none detected")
        self.rOrf.setText(self._accent_link("#orfs", f"{len(rec.get('orfs', [])):,}"))

    def _set_body(self, layout, widget):
        _clear_layout(layout)                                 # recursive: also removes any addLayout'd sub-layouts
        layout.addWidget(widget)

    def _running_jobs(self):
        """Names of the long jobs still running: they outlive the window (the engine's thread pool is torn down only
        when they finish), so a user who closes TEagle under them should know."""
        return [name for flag, name in (("_annot_inflight", "whole-genome annotation"),
                                        ("_genome_prep_inflight", "genome download"),
                                        ("_genome_inflight", "genome scan")) if getattr(self, flag, False)]

    def _confirm_close(self, jobs):
        """True to close. The one seam harnesses patch; a window that is not on screen never asks."""
        if not self.isVisible():
            return True
        box = QMessageBox(self)
        box.setIcon(QMessageBox.Warning); box.setWindowTitle("Close TEagle?")
        box.setText("Still running: " + ", ".join(jobs) + ".")
        box.setInformativeText("Closing the window does not stop the job. TEagle may keep running in the background "
                               "until it finishes. Close anyway?")
        keep = box.addButton("Keep open", QMessageBox.RejectRole)
        box.addButton("Close anyway", QMessageBox.AcceptRole)
        box.setDefaultButton(keep)
        box.exec()
        return box.clickedButton() is not keep
