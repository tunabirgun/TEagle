"""TEagle native desktop app (PySide6). QMainWindow shell with a specimen rail and a scrollable
column of collapsible result cards. All science runs in-process through the shared engine, off the
GUI thread via engine_worker.Engine. This module wires the analyze workflow, primer/PCR, WSL family
annotation, splice detection, provenance and exports."""
from __future__ import annotations
import gzip, os, re, sys

from PySide6.QtCore import Qt, QTimer, QByteArray, QSettings
from PySide6.QtGui import QGuiApplication, QFont, QPixmap, QPainter, QIcon, QCursor, QShortcut, QKeySequence
from PySide6.QtSvg import QSvgRenderer
from PySide6.QtWidgets import (QApplication, QMainWindow, QWidget, QFrame, QVBoxLayout, QHBoxLayout,
                               QGridLayout, QLabel, QLineEdit, QTextEdit, QPlainTextEdit, QPushButton, QComboBox,
                               QScrollArea, QSplitter, QFileDialog, QSizePolicy, QToolTip, QMessageBox, QDialog,
                               QTableWidget, QTableWidgetItem, QHeaderView, QAbstractItemView, QSpinBox)

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.join(os.path.dirname(_HERE), "backend"))

from engine_worker import Engine
import fonts
import figures
from figures import svg_genome
import widgets
from widgets import FigurePanel, GenomePanel, DataTable, BusyBar
from sample import make_sample
import theme as theme_mod
from teagle_core import appdirs
from teagle_core.fetch import COORD_ASSEMBLIES, complete_gene_model, cross_check_models   # pinned assemblies + gene model
from teagle_core import __version__ as APP_VERSION    # single source of truth (never hardcode a duplicate version)
# common model organisms for RepeatMasker/Dfam lineage (display, value passed to -species).
# 'Other…' at the end reveals a free-text field for anything not listed.
WSL_ORGANISMS = [
    ("Human", "Homo sapiens"), ("Mouse", "Mus musculus"), ("Rat", "Rattus norvegicus"),
    ("Zebrafish", "Danio rerio"), ("Chicken", "Gallus gallus"), ("Frog", "Xenopus tropicalis"),
    ("Fruit fly", "Drosophila melanogaster"), ("Mosquito", "Anopheles gambiae"),
    ("Nematode", "Caenorhabditis elegans"), ("Honey bee", "Apis mellifera"),
    ("Thale cress", "Arabidopsis thaliana"), ("Rice", "Oryza sativa"), ("Maize", "Zea mays"),
    ("Wheat", "Triticum aestivum"), ("Cow", "Bos taurus"), ("Dog", "Canis lupus familiaris"),
    ("Budding yeast", "Saccharomyces cerevisiae"),
]
MARK_H = 36                                               # header brand-mark height (px)
WORD_H = 24                                               # header wordmark height (px); a touch smaller than the eagle mark
ICON_TEAL = "#12B39A"                                     # mid-teal for OS chrome (reads on light + dark taskbars)

def _load_asset(name: str) -> str:
    path = appdirs.resource("native", "assets", name) or os.path.join(_HERE, "assets", name)
    with open(path, encoding="utf-8") as f:
        return f.read()

def _svg_pixmap(svg: str, height: int, dpr: float = 1.0) -> QPixmap:
    r = QSvgRenderer(QByteArray(svg.encode("utf-8")))
    s = r.defaultSize(); w = max(1, round(height * s.width() / s.height()))
    pm = QPixmap(round(w * dpr), round(height * dpr)); pm.fill(Qt.transparent)
    p = QPainter(pm); r.render(p); p.end(); pm.setDevicePixelRatio(dpr)
    return pm

# brand mark: single-color eagle logo (fill=currentColor), recolored per theme
_MARK_SVG = None
def _mark_pixmap(color: str, height: int, dpr: float = 1.0) -> QPixmap:
    global _MARK_SVG
    if _MARK_SVG is None:
        _MARK_SVG = _load_asset("teagle-mark.svg")
    return _svg_pixmap(_MARK_SVG.replace("currentColor", color), height, dpr)

# wordmark: clean Cascadia Code Bold 'TEagle' frozen to static paths, TE/agle recolored per theme
_WORD_SVG = None
def _word_pixmap(te: str, agle: str, height: int, dpr: float = 1.0) -> QPixmap:
    global _WORD_SVG
    if _WORD_SVG is None:
        _WORD_SVG = _load_asset("teagle-wordmark.svg")
    return _svg_pixmap(_WORD_SVG.replace("{TE}", te).replace("{AGLE}", agle), height, dpr)

def _app_icon() -> QIcon:                                 # OS window/taskbar icon: the bundled multi-frame ICO
    # use the same crisp 16..256 teagle.ico the exe + installer shortcut embed, so the running window icon is
    # pixel-identical to the shortcut icon (which is what Windows needs to bind the taskbar button via the AUMID)
    path = appdirs.resource("teagle.ico") or os.path.join(_HERE, "..", "..", "installer", "teagle.ico")
    return QIcon(path)

STRUCT_COLS = ["Feature", "Coords (0-based)", "Len", "Metric", "Method"]
ORF_COLS = ["Strand", "Frame", "Start", "End", "aa"]
DOMAIN_COLS = ["Domain", "Label", "Pfam", "aa", "nt", "Score", "E-value", "Conf"]
# plain-language glossary for every table header (hover to learn the abbreviation) — mirrors web GLOSSARY
GLOSS = {
    "Feature": "Structural hallmark found in the sequence — e.g. LTR, TIR, target-site duplication, poly-A tail.",
    "Coords (0-based)": "Location in the sequence, 0-based half-open [start, end).",
    "Coords": "Location of this amplicon in the searched sequence, 0-based half-open [start, end).",
    "Len": "Length in base pairs.",
    "Metric": "Feature-specific measure — terminal-repeat identity %, a motif, or a length.",
    "Method": "How TEagle detected this feature (the algorithm/heuristic used).",
    "Strand": "Strand of the feature — + forward, − reverse complement.",
    "Frame": "Reading frame (1–3) the ORF is translated in.",
    "Start": "Start position in the sequence, 0-based.",
    "End": "End position in the sequence, 0-based half-open.",
    "aa": "Length of the predicted protein (open reading frame) in amino acids.",
    "Domain": "Detected protein domain code (RT, INT, RNaseH, TPase, …).",
    "Label": "Human-readable name of the protein domain.",
    "Pfam": "Pfam accession for the domain profile that matched.",
    "nt": "Nucleotide span of the domain in the sequence, 0-based.",
    "Score": "HMMER bit score — how strongly this region matches the domain profile; higher is stronger.",
    "E-value": "Expected number of matches this good by chance — lower is more significant (e.g. 1e-30 is highly significant).",
    "Conf": "Per-domain call confidence from the HMMER i-Evalue (Eddy 2011): high (≤ 1e-10) or moderate. This is the per-tool reliability of THIS domain call — separate from the element-level structural-completeness tier.",
    "ID": "Identifier of this designed primer pair.",
    "Forward (5'→3')": "Forward primer sequence, written 5′→3′.",
    "Reverse (5'→3')": "Reverse primer sequence, written 5′→3′.",
    "Product": "Predicted amplicon (product) size in base pairs.",
    "Tm F/R": "Melting temperature (°C) of the forward / reverse primer — matched so both anneal together.",
    "GC% F/R": "Percent G+C content of the forward / reverse primer; ~40–60% is typical.",
    "Penalty": "Primer3's overall penalty for the pair — lower is better; it rises as primers depart from target Tm, size and GC.",
    "Hairpin": "Most stable self-folding (hairpin) ΔG of the forward/reverse primer, kcal/mol (worst of the two). "
               "More negative = more stable structure = worse; ΔG ≤ −9 is flagged. ‡ = the two engines disagree — hover for both.",
    "Self-dim": "Most stable self-dimer (primer with a copy of itself) ΔG, kcal/mol (worst of F/R). "
                "More negative = worse; ≤ −9 flagged. ‡ marks an engine disagreement.",
    "Hetero": "Cross-dimer ΔG between the forward and reverse primer, kcal/mol. More negative = worse; ≤ −9 flagged. "
              "A 3′-end cross-dimer is the classic primer-dimer. ‡ marks an engine disagreement.",
    "3′-end": "3′-end anneal stability (last bases) of the pair, kcal/mol — the strongest single predictor of primer-dimer, "
              "because only a base-paired 3′ end is extended by polymerase. More negative = worse.",
    "Struct": "Worst secondary-structure flag for the pair — ok / caution / warn — across hairpin, self-dimer, hetero-dimer and 3′-end.",
    "Pair": "Which designed primer pair produced this amplicon.",
    "Source": "The sequence that was searched — your specimen or a custom background.",
    "Mism F/R": "Mismatches in the forward / reverse primer binding site (the 3′ end is kept exact).",
    "Call": "On-target = amplicon at the intended locus; off-target = amplified elsewhere; priming site = a "
            "genome-wide product with no single intended locus (a bare-consensus whole-genome scan).",
    "Class/family": "TE class and superfamily (Wicker 2007 scheme), e.g. LTR/Copia.",
    "Dfam family": "The specific named family in the Dfam library, e.g. Copia_I or L1HS.",
    "Str": "Strand of the match — + forward, − reverse complement.",
    "Div": "Divergence — % difference between your sequence and the Dfam family consensus (lower = closer).",
    "Intron span (0-based)": "Intron location in the loaded sequence, 0-based half-open [start, end).",
    "Splice site": "The two bases at each intron boundary (donor…acceptor); canonical introns are GT…AG (or GC–AG / AT–AC).",
    "Canonical": "Whether the intron's donor…acceptor matches a canonical eukaryotic splice motif.",
    "#": "Row number.",
}

_FLAG_ORDER = {"ok": 0, "caution": 1, "warn": 2}       # colours come from the per-theme palette (theme_mod.FLAG), applied at render


def _metric_cell(parts):
    """Fold one or more ΔG metric dicts ({p3, vrna, flag, agree}) into a table cell:
    (text, worst_flag, tooltip). Shows the worst (most negative) ΔG across engines/primers, marks an engine
    disagreement with ‡; the caller maps worst_flag -> a per-theme colour. The Struct column carries the flag
    as TEXT, so the call is legible without relying on colour."""
    dgs, flags, disagree, tips = [], [], False, []
    for lab, m in parts:
        if not m:
            continue
        p3, vr = m.get("p3"), m.get("vrna")
        vals = [v for v in (p3, vr) if v is not None]
        if vals:
            dgs.append(min(vals))
        flags.append(m.get("flag", "ok"))
        if m.get("agree") == "disagree":
            disagree = True
        pre = (lab + ": ") if lab else ""
        tips.append(f"{pre}primer3 {p3 if p3 is not None else '—'} / ViennaRNA {vr if vr is not None else '—'} kcal/mol"
                    + (f" ({m.get('agree')})" if m.get("agree") and m.get("agree") not in ("none",) else ""))
    if not dgs:
        return ("—", "ok", "no structure predicted")
    worst_flag = max(flags, key=lambda f: _FLAG_ORDER.get(f, 0)) if flags else "ok"
    txt = f"{min(dgs):.1f}" + ("‡" if disagree else "")
    return (txt, worst_flag, " · ".join(tips))


# clickable source citations (verified DOIs — mirror backend refs.py and the web REFLINKS)
REFLINKS = {
    "Wicker2007":   {"url": "https://doi.org/10.1038/nrg2165", "cite": "Wicker T, et al. (2007) A unified classification system for eukaryotic transposable elements. Nat Rev Genet 8:973-982."},
    "Pfam":         {"url": "https://www.ebi.ac.uk/interpro/", "cite": "Mistry J, et al. (2021) Pfam: the protein families database in 2021. Nucleic Acids Res 49:D412-D419."},
    "HMMER":        {"url": "https://doi.org/10.1371/journal.pcbi.1002195", "cite": "Eddy SR (2011) Accelerated Profile HMM Searches. PLoS Comput Biol 7:e1002195."},
    "Dfam":         {"url": "https://doi.org/10.1186/s13100-020-00230-y", "cite": "Storer J, et al. (2021) The Dfam community resource of transposable element families. Mob DNA 12:2."},
    "RepeatMasker": {"url": "https://www.repeatmasker.org/", "cite": "Smit AFA, Hubley R, Green P. RepeatMasker Open-4.0."},
    "NCBI":         {"url": "https://www.ncbi.nlm.nih.gov/nuccore/", "cite": "NCBI Entrez / E-utilities (Sayers E, NCBI)."},
    "ENA":          {"url": "https://www.ebi.ac.uk/ena/browser/view/", "cite": "European Nucleotide Archive (EMBL-EBI) — sequence fallback source."},
    "Primer3":      {"url": "https://doi.org/10.1093/nar/gks596", "cite": "Untergasser A, et al. (2012) Primer3 — new capabilities and interfaces. Nucleic Acids Res 40:e115."},
    "minimap2":     {"url": "https://doi.org/10.1093/bioinformatics/bty191", "cite": "Li H (2018) Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics 34:3094-3100."},
    "SantaLucia1998": {"url": "https://doi.org/10.1073/pnas.95.4.1460", "cite": "SantaLucia J Jr (1998) A unified view of … DNA nearest-neighbor thermodynamics. PNAS 95(4):1460-1465."},
    "Owczarzy2008": {"url": "https://doi.org/10.1093/nar/gkn198", "cite": "Owczarzy R, et al. (2008) IDT SciTools. Nucleic Acids Res 36(Web Server):W163-W169. — comparability only: TEagle matches OligoAnalyzer's ΔG convention + −9 kcal/mol threshold, it does not run IDT SciTools (ΔG computed with Primer3 + ViennaRNA)."},
    "ViennaRNA":    {"url": "https://doi.org/10.1186/1748-7188-6-26", "cite": "Lorenz R, et al. (2011) ViennaRNA Package 2.0. Algorithms Mol Biol 6:26."},
}


class CollapsibleCard(QFrame):
    """A titled result card that expands/collapses on header click. Starts collapsed; reveal_on_data
    auto-expands the first time content is set, mirroring the web UI's progressive reveal."""
    def __init__(self, number: str, title: str, meta: str = "", collapsed=True):
        super().__init__()
        self.setObjectName("card")
        self._lay = QVBoxLayout(self)
        self._lay.setContentsMargins(0, 0, 0, 0)
        self._lay.setSpacing(0)
        self.hdr = QPushButton()
        self.hdr.setObjectName("cardhdr")
        self.hdr.setText(f"  {number}   {title}    ·  {meta}" if meta else f"  {number}   {title}")
        self.hdr.setCheckable(True)
        self.hdr.clicked.connect(self._toggle)
        self._lay.addWidget(self.hdr)
        self.body = QWidget()
        self.bodylay = QVBoxLayout(self.body)
        self.bodylay.setContentsMargins(11, 4, 11, 12)
        self.bodylay.setSpacing(9)
        self._lay.addWidget(self.body)
        self._number, self._title, self._meta = number, title, meta
        self.set_collapsed(collapsed)

    def _toggle(self):
        self.set_collapsed(self.body.isVisible())

    def set_collapsed(self, collapsed: bool):
        self.body.setVisible(not collapsed)
        arrow = "▸" if collapsed else "▾"
        title = self._title.replace("&", "&&")           # QPushButton eats a lone '&' as a mnemonic
        meta = self._meta.replace("&", "&&") if self._meta else ""
        txt = f"{arrow} {self._number}  {title}" + (f"    ·  {meta}" if meta else "")
        self.hdr.setText(txt.upper())                     # stay ALL CAPS whether collapsed or expanded; keep the meta gloss

    def expand(self):
        self.set_collapsed(False)

    def clear_body(self):
        _clear_layout(self.bodylay)


def _clear_layout(layout):
    """Recursively remove every widget AND nested sub-layout from `layout`. A widget-only clear
    (setParent(None) on it.widget()) leaves addLayout'd sub-layouts — the 'Structural evidence' /
    'Protein domains' section headers and their source links — orphaned but still parented to the body
    widget, so they accumulate and DUPLICATE on every re-render (e.g. a second Run analysis)."""
    while layout.count():
        it = layout.takeAt(0)
        w = it.widget()
        if w is not None:
            w.setParent(None)
        else:
            sub = it.layout()
            if sub is not None:
                _clear_layout(sub)


def _hline():
    f = QFrame(); f.setObjectName("hline"); f.setFixedHeight(1); return f


def _empty(text):
    l = QLabel(text); l.setObjectName("empty"); l.setWordWrap(True); return l


def _sl(text):
    """A section label — mono, uppercase, tracked (the web UI's `.lbl`)."""
    l = QLabel(text.upper()); l.setObjectName("sectionlabel"); l.setWordWrap(True); return l


def _export_table_btn(table, base, parent):
    """Visible Excel/CSV/TSV export for a DataTable: the button pops a format menu, then a save dialog
    pre-set to the chosen type. Exports in the table's current on-screen (sorted) order."""
    b = QPushButton("⭳ Export table ▾"); b.setProperty("sm", True)
    def pop():
        fmt = widgets.pick_table_format(b, b.mapToGlobal(b.rect().bottomLeft()))
        if fmt:
            widgets.export_table(table._headers, table.rows_data(), base, parent, fmt=fmt)
    b.clicked.connect(lambda _=False: pop())
    return b


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("TEagle")
        self.setWindowIcon(_app_icon())                   # also set on the top-level window (some Qt/Windows paths key the taskbar HICON off WM_SETICON)
        # open at the design size but never larger than the screen — a 1366x768 laptop must not get an
        # 860px-tall window taller than its display; keep a usable floor so the split layout stays coherent.
        # Clamp the MINIMUM to the screen too: setMinimumSize is logical px, so an 820x560 floor becomes
        # 840px physical at 1.5x UI scale and would overhang a 768px screen if the min itself is not clamped.
        w, h = round(1240 * theme_mod.UI_SCALE), round(860 * theme_mod.UI_SCALE)
        minw, minh = 820, 560
        scr = self.screen().availableGeometry() if self.screen() else None
        if scr is not None:
            aw, ah = scr.width() - 40, scr.height() - 60
            w, h = min(w, aw), min(h, ah)
            minw, minh = min(minw, aw), min(minh, ah)
        self.setMinimumSize(minw, minh)
        self.resize(max(w, minw), max(h, minh))
        self.theme = "dark"
        self.state = {"seq": "", "source": None, "last_rec": None}
        self._loading = False                                 # True while a programmatic load writes the specimen box
        self._pcr_gen = 0                                     # monotonic in-silico-PCR batch id (drops stale sibling results)
        self._design_inflight = False                         # one primer design at a time (self._design_tmpl is shared state)
        self._genome_inflight = False                         # one whole-genome isPcr scan at a time
        self._genome_prep_inflight = False                    # one genome download/prepare at a time (large, one-time)
        self._pending_scan = None                             # a scan queued behind a just-started genome download
        self._prepared_genomes = []                           # downloaded+verified (.done) genomes; the ONLY source for the PCR organism dropdown

        self.engine = Engine(self)
        self.engine.done.connect(self._on_done)
        self.engine.user_error.connect(self._on_user_error)
        self.engine.failed.connect(self._on_failed)

        central = QWidget(); central.setObjectName("central")
        self.setCentralWidget(central)
        outer = QVBoxLayout(central)
        outer.setContentsMargins(12, 10, 12, 10)
        outer.setSpacing(8)
        outer.addWidget(self._build_header())

        split = QSplitter(Qt.Horizontal)
        self.rail = self._build_rail()
        split.addWidget(self.rail)
        split.addWidget(self._build_results())
        split.setStretchFactor(0, 0)
        split.setStretchFactor(1, 1)
        split.setSizes([340, 900])
        self.split = split
        self._rail_sizes = [340, 900]                          # remembered split for reopening the collapsed rail
        self._rail_collapsed = False                           # explicit state (isVisible() is unreliable pre-show)
        outer.addWidget(split, 1)

        self._apply_theme()
        QTimer.singleShot(0, self._startup)

    # ---------- header ----------
    def resizeEvent(self, e):
        # hide the decorative tagline on a narrow window so it never contributes to an overhang; the Ignored
        # size policy already lets the window shrink, this just keeps the header clean at small widths.
        super().resizeEvent(e)
        tag = getattr(self, "_tagline", None)
        if tag is not None:
            tag.setVisible(self.width() >= round(1080 * theme_mod.UI_SCALE))

    def _build_header(self):
        wrap = QWidget()
        col = QVBoxLayout(wrap); col.setContentsMargins(6, 0, 2, 0); col.setSpacing(0)
        h = QHBoxLayout(); h.setContentsMargins(0, 2, 0, 8); h.setSpacing(10)
        self.railToggle = QPushButton("◧"); self.railToggle.setProperty("sm", True)
        self.railToggle.setToolTip("Hide the specimen panel for more analysis width (Ctrl+B)")
        self.railToggle.clicked.connect(self._toggle_rail)
        h.addWidget(self.railToggle)
        self.mark = QLabel()                                  # eagle brand mark; pixmap set per-theme in _apply_theme
        self.mark.setObjectName("mark")
        self.mark.setToolTip("TEagle")
        h.addWidget(self.mark)
        self.word = QLabel()                                  # Cascadia Code wordmark; pixmap set per-theme in _apply_theme
        self.word.setObjectName("word")
        h.addWidget(self.word)
        self.ver = QLabel("v" + APP_VERSION); self.ver.setObjectName("ver")
        h.addWidget(self.ver)
        tag = QLabel("TRANSPOSABLE ELEMENTS ASSAY TERMINAL"); tag.setObjectName("tagline")
        tf = tag.font(); tf.setLetterSpacing(QFont.AbsoluteSpacing, 1.5); tag.setFont(tf)
        # Ignored horizontal policy so the decorative tagline never imposes its full text width as the window's
        # hard minimum — otherwise a small screen at ≥125% UI scale cannot shrink the window inside the display.
        tag.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Preferred)
        self._tagline = tag                                   # hidden on a narrow window (resizeEvent) so it never overhangs
        h.addWidget(tag)
        h.addStretch(1)
        chip = QFrame(); chip.setObjectName("statuschip")
        cl = QHBoxLayout(chip); cl.setContentsMargins(10, 5, 11, 5); cl.setSpacing(8)
        self.led = QLabel(); self.led.setObjectName("led"); self.led.setFixedSize(8, 8)
        self.statusTxt = QLabel("connecting…"); self.statusTxt.setObjectName("statusTxt")
        cl.addWidget(self.led); cl.addWidget(self.statusTxt)
        h.addWidget(chip)
        sc = QPushButton("⤢ SCALE"); sc.setProperty("sm", True)
        sc.setToolTip("Global UI scale — shrink or enlarge the whole interface (applied on restart)")
        sc.clicked.connect(lambda: self._ui_scale_menu(sc))
        h.addWidget(sc)
        tb = QPushButton("◐ THEME"); tb.setProperty("sm", True); tb.clicked.connect(self._toggle_theme)
        h.addWidget(tb)
        col.addLayout(h)
        self.headrule = QFrame(); self.headrule.setObjectName("headrule"); self.headrule.setFixedHeight(2)
        col.addWidget(self.headrule)
        QShortcut(QKeySequence("Ctrl+B"), self, activated=self._toggle_rail)   # toggle the specimen panel
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
            self.railToggle.setText("▣")
            self.railToggle.setToolTip("Show the specimen panel (Ctrl+B)")
        else:
            self.rail.setVisible(True)
            total = sum(self.split.sizes()) or 1240
            sizes = self._rail_sizes if sum(self._rail_sizes) <= total else [340, max(300, total - 340)]
            self.split.setSizes(sizes)
            self.railToggle.setText("◧")
            self.railToggle.setToolTip("Hide the specimen panel for more analysis width (Ctrl+B)")

    def _ui_scale_menu(self, anchor):
        """Pick a global UI scale. Persisted to QSettings and applied on the next launch (QT_SCALE_FACTOR must
        be set before the QApplication exists), so a small-screen user can make the whole UI fit."""
        from PySide6.QtWidgets import QMenu
        cur = 1.0
        try:
            cur = float(QSettings("TEagle", "TEagle").value("ui_scale", 1.0))
        except Exception:
            pass
        m = QMenu(self)
        for f in UI_SCALES:
            act = m.addAction(f"{int(f * 100)}%" + ("   ✓" if abs(f - cur) < 1e-3 else ""))
            act.triggered.connect(lambda _=False, x=f: self._set_ui_scale(x))
        m.exec(anchor.mapToGlobal(anchor.rect().bottomLeft()))

    def _set_ui_scale(self, factor):
        QSettings("TEagle", "TEagle").setValue("ui_scale", float(factor))
        box = QMessageBox(self)
        box.setWindowTitle("UI scale")
        box.setText(f"UI scale set to {int(factor * 100)}%.\n\nThe new scale applies on the next launch — restart now, or later.")
        restart = box.addButton("Restart now", QMessageBox.AcceptRole)
        box.addButton("Later", QMessageBox.RejectRole)
        box.exec()
        if box.clickedButton() is restart:
            self._restart_app()

    def _restart_app(self):
        """Relaunch TEagle so a new UI scale takes effect (QT_SCALE_FACTOR is read before the QApplication)."""
        try:
            import subprocess
            if getattr(sys, "frozen", False):                 # packaged build: the exe is sys.executable
                subprocess.Popen([sys.executable])
            else:
                subprocess.Popen([sys.executable, os.path.abspath(sys.argv[0])] + sys.argv[1:])
            QApplication.instance().quit()
        except Exception as e:
            QMessageBox.information(self, "Restart", f"Please restart TEagle manually to apply the scale.\n({type(e).__name__}: {e})")

    # ---------- rail ----------
    def _build_rail(self):
        rail = QFrame(); rail.setObjectName("rail")
        rail.setMinimumWidth(300); rail.setMaximumWidth(430)
        lay = QVBoxLayout(rail); lay.setContentsMargins(12, 12, 12, 12); lay.setSpacing(8)
        lay.addWidget(self._sec("01", "Specimen"))
        accrow = QHBoxLayout()
        self.acc = QLineEdit(); self.acc.setPlaceholderText("accession — e.g. M11240, NC_003075.7")
        accrow.addWidget(self.acc)
        self.fetchBtn = QPushButton("↓ Fetch"); self.fetchBtn.setProperty("sm", True); self.fetchBtn.clicked.connect(self._fetch)
        accrow.addWidget(self.fetchBtn)
        lay.addLayout(accrow)
        self.accMeta = QLabel(""); self.accMeta.setObjectName("cardmeta"); self.accMeta.setWordWrap(True)
        self.accMeta.setTextFormat(Qt.RichText); self.accMeta.setOpenExternalLinks(True)
        lay.addWidget(self.accMeta)

        # coordinate fetch (collapsed) — organism + chr:start-end, like the UCSC browser position box
        self.coordToggle = QPushButton("⌖ Fetch by coordinate ▾"); self.coordToggle.setProperty("link", True)
        self.coordToggle.clicked.connect(self._toggle_coord); lay.addWidget(self.coordToggle)
        self.coordBox = QWidget(); cb = QVBoxLayout(self.coordBox); cb.setContentsMargins(0, 2, 0, 2); cb.setSpacing(5)
        orow = QHBoxLayout()
        self.asmSel = QComboBox()
        for org in sorted(COORD_ASSEMBLIES):
            self.asmSel.addItem(f"{org} · {COORD_ASSEMBLIES[org]['assemblyName']}", org)
        self.asmSel.addItem("Other organism / assembly…", "__custom__")
        orow.addWidget(self.asmSel, 1)
        self.coordStrand = QComboBox(); self.coordStrand.addItems(["+ strand", "− strand"])
        self.coordStrand.setMaximumWidth(104); orow.addWidget(self.coordStrand)
        cb.addLayout(orow)
        self.coordCustom = QLineEdit(); self.coordCustom.setPlaceholderText("organism name or assembly accession (e.g. GCF_000001405.40)")
        self.coordCustom.setVisible(False); cb.addWidget(self.coordCustom)
        self.asmSel.currentIndexChanged.connect(lambda _=0: self.coordCustom.setVisible(self.asmSel.currentData() == "__custom__"))
        self.coord = QPlainTextEdit(); self.coord.setMaximumHeight(66)
        self.coord.setPlaceholderText("chr13:33,016,423-33,066,143   (one region per line for multi-region)")
        cb.addWidget(self.coord)
        crow = QHBoxLayout()
        self.coordFetchBtn = QPushButton("↓ Fetch region(s)"); self.coordFetchBtn.setProperty("sm", True)
        self.coordFetchBtn.clicked.connect(self._fetch_coord)
        crow.addWidget(self.coordFetchBtn); crow.addStretch(1); cb.addLayout(crow)
        self.coordMeta = QLabel(""); self.coordMeta.setObjectName("cardmeta"); self.coordMeta.setWordWrap(True)
        self.coordMeta.setTextFormat(Qt.RichText); self.coordMeta.setOpenExternalLinks(True); cb.addWidget(self.coordMeta)
        cnote = QLabel("UCSC-style, 1-based (same numbers as the browser). Multi-region: all fetched + recorded; "
                       "analysis runs on the first region.")
        cnote.setObjectName("orient"); cnote.setWordWrap(True); cb.addWidget(cnote)
        self.coordBox.setVisible(False); lay.addWidget(self.coordBox)

        ub = QPushButton("⭱ Upload FASTA (.fa / .fasta / .gz)"); ub.setProperty("sm", True)
        ub.clicked.connect(self._upload); lay.addWidget(ub)
        self.seq = QTextEdit(); self.seq.setPlaceholderText("…or paste DNA (FASTA or raw). Real IUPAC validation runs on analyze.")
        self.seq.setMinimumHeight(120); self.seq.textChanged.connect(self._seq_changed)
        lay.addWidget(self.seq)
        row = QHBoxLayout()
        ls = QPushButton("↳ load sample element"); ls.setProperty("link", True); ls.clicked.connect(self._load_sample)
        row.addWidget(ls); row.addStretch(1)
        self.charCount = QLabel("0 nt"); self.charCount.setObjectName("kdim"); row.addWidget(self.charCount)
        lay.addLayout(row)
        self.runBtn = QPushButton("▶ Run analysis"); self.runBtn.setProperty("primary", True)
        self.runBtn.clicked.connect(self._run_analysis); lay.addWidget(self.runBtn)

        # readout gauges (2×2)
        lay.addSpacing(6)
        grid = QGridLayout(); grid.setHorizontalSpacing(6); grid.setVerticalSpacing(6)
        grid.setColumnStretch(0, 1); grid.setColumnStretch(1, 1)
        self.mLen = self._readout(grid, 0, "Length"); self.mGC = self._readout(grid, 1, "GC")
        self.mN = self._readout(grid, 2, "N content"); self.mValid = self._readout(grid, 3, "IUPAC")
        lay.addLayout(grid)
        lay.addSpacing(4)
        self.rRecords = QLabel("—"); self.rStruct = QLabel("—"); self.rOrf = QLabel("—")
        for val in (self.rRecords, self.rStruct, self.rOrf):
            val.setObjectName("cardmeta")
        for lbl, val in (("RECORDS", self.rRecords), ("STRUCTURAL EVIDENCE", self.rStruct), ("ORFS (≥40 aa)", self.rOrf)):
            r = QHBoxLayout(); k = QLabel(lbl); k.setObjectName("kdim"); r.addWidget(k); r.addStretch(1); r.addWidget(val)
            lay.addLayout(r)

        lay.addSpacing(6)
        self.envHdr = QPushButton("▸ ◆  Environment"); self.envHdr.setObjectName("cardhdr")
        self.envHdr.clicked.connect(self._toggle_env); lay.addWidget(self.envHdr)
        self.envBox = QLabel("checking…"); self.envBox.setObjectName("cardmeta"); self.envBox.setWordWrap(True)
        self.envBox.setVisible(False); self.envBox.setTextFormat(Qt.RichText)
        lay.addWidget(self.envBox)
        lay.addStretch(1)
        note = QLabel("Superfamily (Copia / Gypsy / LINE / DNA) is called from protein-domain "
                      "architecture (HMMER + CC0 Pfam profiles). Dfam / RepeatMasker family naming "
                      "runs in the managed WSL backend (panel 03).")
        note.setObjectName("orient"); note.setWordWrap(True); lay.addWidget(note)

        wrap = QScrollArea(); wrap.setWidgetResizable(True); wrap.setWidget(rail)
        wrap.setMinimumWidth(320); wrap.setMaximumWidth(440)
        return wrap

    def _sec(self, n, title):
        w = QWidget(); r = QHBoxLayout(w); r.setContentsMargins(0, 4, 0, 2); r.setSpacing(9)
        num = QLabel(n); num.setObjectName("secn"); num.setAlignment(Qt.AlignCenter); r.addWidget(num)
        h = QLabel(title); h.setObjectName("sech"); r.addWidget(h); r.addStretch(1)
        return w

    def _readout(self, grid, idx, label):
        cell = QFrame(); cell.setObjectName("cell")
        cl = QVBoxLayout(cell); cl.setContentsMargins(10, 8, 10, 9); cl.setSpacing(2)
        k = QLabel(label.upper()); k.setObjectName("kdim")
        v = QLabel("—"); v.setObjectName("value")
        cl.addWidget(k); cl.addWidget(v)
        r, c = divmod(idx, 2)
        grid.addWidget(cell, r, c)
        return v

    # ---------- results column ----------
    def _build_results(self):
        wrap = QScrollArea(); wrap.setWidgetResizable(True)
        wrap.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)   # narrow window: scroll, never hard-clip content
        self.resultsScroll = wrap                          # so send-to-splice etc. can scroll a card into view
        inner = QWidget(); inner.setObjectName("central")
        self.results = QVBoxLayout(inner); self.results.setContentsMargins(4, 0, 8, 0); self.results.setSpacing(9)
        self.errbanner = QLabel(""); self.errbanner.setObjectName("errbanner"); self.errbanner.setWordWrap(True)
        self.errbanner.setVisible(False); self.results.addWidget(self.errbanner)

        self.card_struct = CollapsibleCard("02", "Classification & structure",
                                           "LTR/TIR repeats, ORFs, protein domains")
        self.card_struct.bodylay.addWidget(_empty("Run analysis to detect terminal repeats, ORFs and tails."))
        self.results.addWidget(self.card_struct)

        self.results.addWidget(self._build_wsl_card())
        self.results.addWidget(self._build_splice_card())
        self.results.addWidget(self._build_primer_card())
        self.results.addWidget(self._build_pcr_card())
        self.results.addWidget(self._build_genome_card())
        self.card_prov = CollapsibleCard("07", "Run provenance", "versions + checksums for reproducibility")
        self.card_prov.bodylay.addWidget(_empty("Populated from the last computation — travels with every result."))
        self.results.addWidget(self.card_prov)

        self.results.addStretch(1)
        wrap.setWidget(inner)
        return wrap

    # ---------- 03 Dfam / RepeatMasker family (WSL) ----------
    def _build_wsl_card(self):
        card = CollapsibleCard("03", "Dfam / RepeatMasker family", "names the TE family (Dfam)")
        self.card_wsl = card
        self.wslStatus = QLabel("checking WSL backend…"); self.wslStatus.setObjectName("cardmeta")
        self.wslStatus.setWordWrap(True); self.wslStatus.setTextFormat(Qt.RichText)
        card.bodylay.addWidget(self.wslStatus)
        srcrow = QHBoxLayout()
        srcrow.addWidget(QLabel("Sequence source"))
        self.wslSource = QComboBox()
        self.wslSource.addItems(["Loaded specimen (panel 01)", "Paste sequence…"])
        srcrow.addWidget(self.wslSource); srcrow.addStretch(1)
        card.bodylay.addLayout(srcrow)
        self.wslPaste = QTextEdit(); self.wslPaste.setPlaceholderText("Paste a FASTA or raw DNA sequence to annotate against Dfam")
        self.wslPaste.setMaximumHeight(70); self.wslPaste.setVisible(False)
        self.wslSource.currentIndexChanged.connect(lambda i: self.wslPaste.setVisible(i == 1))
        card.bodylay.addWidget(self.wslPaste)
        row = QHBoxLayout()
        row.addWidget(QLabel("Organism"))
        self.wslSpecies = QComboBox()                     # dropdown of common organisms + 'Other…'
        self.wslSpecies.addItem("— select organism —", None)
        for common, sci in WSL_ORGANISMS:
            self.wslSpecies.addItem(f"{common} · {sci}", sci)
        self.wslSpecies.addItem("Other…", "__other__")
        self.wslSpecies.currentIndexChanged.connect(lambda _i: self._on_species_changed())
        row.addWidget(self.wslSpecies, 1)
        self.wslSpeciesOther = QLineEdit(); self.wslSpeciesOther.setPlaceholderText("type organism / species")
        self.wslSpeciesOther.setVisible(False)
        row.addWidget(self.wslSpeciesOther, 1)
        self.annotateBtn = QPushButton("▶ Run family annotation"); self.annotateBtn.setProperty("sm", True)
        self.annotateBtn.setEnabled(False); self.annotateBtn.clicked.connect(self._annotate)
        row.addWidget(self.annotateBtn)
        card.bodylay.addLayout(row)
        self.wslInstallBtn = QPushButton("⚙ Backend installer — install · repair · check integrity")
        self.wslInstallBtn.setProperty("sm", True)
        self.wslInstallBtn.clicked.connect(self._open_installer)
        card.bodylay.addWidget(self.wslInstallBtn)
        self.wslBody = QVBoxLayout(); wb = QWidget(); wb.setLayout(self.wslBody)
        self.wslBody.addWidget(_empty("Run RepeatMasker against Dfam to name the TE family. Family naming is the Linux (WSL) backend."))
        card.bodylay.addWidget(wb)
        return card

    # ---------- Splice detection ----------
    def _build_splice_card(self):
        card = CollapsibleCard("◧", "Splice detection (de novo)", "finds exons & introns from a transcript")
        self.card_splice = card
        self.spliceStatus = QLabel("checking splice-alignment backend…"); self.spliceStatus.setObjectName("cardmeta")
        self.spliceStatus.setWordWrap(True); self.spliceStatus.setTextFormat(Qt.RichText)
        card.bodylay.addWidget(self.spliceStatus)
        # the genomic reference is always the specimen loaded in panel 01 (fetched/uploaded/pasted)
        self.spliceRef = QLabel("Genomic reference: none loaded yet — load a specimen in panel 01.")
        self.spliceRef.setObjectName("cardmeta"); self.spliceRef.setWordWrap(True); self.spliceRef.setTextFormat(Qt.RichText)
        card.bodylay.addWidget(self.spliceRef)
        tip = QLabel("Aligns a transcript / cDNA / mRNA to that reference; introns are the alignment gaps. "
                     "Tip: right-click any feature → “Send to splice detection” to use it as the transcript.")
        tip.setObjectName("cardmeta"); tip.setWordWrap(True)
        card.bodylay.addWidget(tip)
        self.spliceTx = QTextEdit()
        self.spliceTx.setPlaceholderText("Paste a transcript / cDNA / mRNA. minimap2 -x splice maps it to the loaded sequence; "
                                         "introns are the alignment gaps, checked against canonical GT–AG splice sites.")
        self.spliceTx.setMaximumHeight(80)
        card.bodylay.addWidget(self.spliceTx)
        self.spliceBtn = QPushButton("▶ Detect exons / introns"); self.spliceBtn.setProperty("sm", True)
        self.spliceBtn.setEnabled(False); self.spliceBtn.clicked.connect(self._splice)
        card.bodylay.addWidget(self.spliceBtn)
        self.spliceBody = QVBoxLayout(); sb = QWidget(); sb.setLayout(self.spliceBody)
        self.spliceBody.addWidget(_empty("Align a transcript to the loaded sequence to resolve exon–intron structure de novo."))
        card.bodylay.addWidget(sb)
        return card

    # ---------- 04 Primer design ----------
    def _build_primer_card(self):
        card = CollapsibleCard("04", "Primer design", "designs PCR primers (Primer3)")
        self.card_primer = card
        prow = QHBoxLayout()
        prow.addWidget(QLabel("Preset"))
        self.pPreset = QComboBox()
        self._preset_keys = ["standard", "qpcr", "highspec", "permissive"]
        self.pPreset.addItems(["Standard PCR", "qPCR (short amplicon)", "High-specificity", "Permissive (hard targets)"])
        self.pPreset.currentIndexChanged.connect(lambda i: self._apply_preset(self._preset_keys[i]))
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
        self.advBox.setVisible(False)
        self.advToggle.clicked.connect(lambda: (self.advBox.setVisible(not self.advBox.isVisible()),
                                                self.advToggle.setText(("▾" if self.advBox.isVisible() else "▸") + " Advanced parameters")))
        card.bodylay.addWidget(self.advToggle)
        card.bodylay.addWidget(self.advBox)
        row = QHBoxLayout()
        self.designBtn = QPushButton("⋯ Design primers"); self.designBtn.setProperty("primary", True)
        self.designBtn.setEnabled(False); self.designBtn.clicked.connect(self._design)
        row.addWidget(self.designBtn)
        rb = QPushButton("↺ reset"); rb.setProperty("sm", True); rb.clicked.connect(lambda: self._apply_preset("standard"))
        row.addWidget(rb)
        self.designHint = QLabel("run analysis first"); self.designHint.setObjectName("kdim")
        row.addWidget(self.designHint); row.addStretch(1)
        card.bodylay.addLayout(row)
        self.primBody = QVBoxLayout(); pb = QWidget(); pb.setLayout(self.primBody)
        card.bodylay.addWidget(pb)
        return card

    # ---------- 05 In-silico PCR ----------
    def _build_pcr_card(self):
        card = CollapsibleCard("05", "In-silico PCR", "predicts which fragments amplify")
        self.card_pcr = card
        card.bodylay.addWidget(_sl("Loaded primer pairs — one gel lane each (in order)"))
        self.pcrQueueBox = QVBoxLayout(); qb = QWidget(); qb.setLayout(self.pcrQueueBox)
        card.bodylay.addWidget(qb)
        srow = QHBoxLayout()
        self.pcrStageAll = QPushButton("+ stage all designed"); self.pcrStageAll.setProperty("sm", True)
        self.pcrStageAll.setEnabled(False); self.pcrStageAll.clicked.connect(self._pcr_stage_all)
        self.pcrClear = QPushButton("✕ clear"); self.pcrClear.setProperty("sm", True)
        self.pcrClear.setEnabled(False); self.pcrClear.clicked.connect(self._pcr_clear)
        srow.addWidget(self.pcrStageAll); srow.addWidget(self.pcrClear); srow.addStretch(1)
        card.bodylay.addLayout(srow)
        grid = QGridLayout(); grid.setHorizontalSpacing(8)
        self._grid_fields(grid, [("pcrMM", "Max mismatches", "2"), ("pcrTP", "3′ exact bases", "5"),
                                 ("pcrPmin", "Prod min (bp)", "70"), ("pcrPmax", "Prod max (bp)", "1000")], cols=4)
        card.bodylay.addLayout(grid)
        card.bodylay.addWidget(_sl("Optional custom background (FASTA) — off-target search"))
        self.pcrBg = QTextEdit(); self.pcrBg.setMaximumHeight(60)
        self.pcrBg.setPlaceholderText("Optional: paste extra background sequence(s) to reveal off-target amplicons.")
        card.bodylay.addWidget(self.pcrBg)
        row = QHBoxLayout()
        self.runPcrBtn = QPushButton("▶ Run loaded pairs"); self.runPcrBtn.setProperty("primary", True)
        self.runPcrBtn.setEnabled(False); self.runPcrBtn.clicked.connect(self._run_pcr)
        row.addWidget(self.runPcrBtn)
        self.pcrHint = QLabel("load one or more pairs, then run"); self.pcrHint.setObjectName("kdim")
        row.addWidget(self.pcrHint); row.addStretch(1)
        card.bodylay.addLayout(row)
        self.pcrBody = QVBoxLayout(); pb = QWidget(); pb.setLayout(self.pcrBody)
        self.pcrBody.addWidget(_empty("Load a pair, then run pair-aware amplicon search."))
        card.bodylay.addWidget(pb)
        return card

    def _build_genome_card(self):
        # dedicated whole-genome off-target scan — LOCAL isPcr against a downloaded RefSeq assembly. The organism
        # dropdown lists ONLY downloaded+verified genomes (from genome_list). On-target = the product at the
        # specimen's own genome locus when it sits in the scanned assembly; the rest are off-target paralogs.
        card = CollapsibleCard("06", "Whole-genome off-target scan", "a designed pair vs a downloaded RefSeq genome (isPcr)")
        self.card_genome = card
        card.bodylay.addWidget(_sl("Organism — local isPcr against a downloaded genome"))
        gorow = QHBoxLayout()
        self.genomeOrg = QComboBox()
        gorow.addWidget(self.genomeOrg, 1)
        self.genomeManageBtn = QPushButton("⚙ Manage genomes"); self.genomeManageBtn.setProperty("sm", True)
        self.genomeManageBtn.clicked.connect(self._open_genome_manager)
        gorow.addWidget(self.genomeManageBtn)
        card.bodylay.addLayout(gorow)
        self.genomeOrgHint = QLabel()
        self.genomeOrgHint.setObjectName("orient"); self.genomeOrgHint.setWordWrap(True)
        card.bodylay.addWidget(self.genomeOrgHint)
        self._refresh_genome_dropdown()                       # render from cached prepared set (empty at first build)
        gnote = QLabel("Right-click a designed pair → “⊕ Scan whole genome” to run a local isPcr scan against a "
                       "downloaded genome; download / manage genomes via ⚙ Manage genomes. Products are listed below, "
                       "on-target first — sealed with the assembly version + checksum.")
        gnote.setObjectName("orient"); gnote.setWordWrap(True); card.bodylay.addWidget(gnote)
        self.genomeBody = QVBoxLayout(); gb = QWidget(); gb.setLayout(self.genomeBody)
        self.genomeBody.addWidget(_empty("Design a primer pair, then right-click it → Scan whole genome."))
        card.bodylay.addWidget(gb)
        return card

    def _grid_fields(self, grid, specs, cols):
        for idx, (fid, label, default) in enumerate(specs):
            r, c = divmod(idx, cols)
            cell = QWidget(); cl = QVBoxLayout(cell); cl.setContentsMargins(0, 0, 0, 0); cl.setSpacing(1)
            lab = QLabel(label.upper()); lab.setObjectName("kdim")
            ed = QLineEdit(default); ed.setMaximumWidth(90)
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
        self._uppercase_buttons()
        self._sync_viewer_themes()

    def _sync_viewer_themes(self):
        """Push the app theme into every live genome viewer + figure panel (gel), by default, so they follow
        dark/light. A live findChildren walk (no persistent registry -> no dangling C++ refs across the
        constant card rebuilds); each viewer re-renders in place, preserving its pan/zoom. Per-viewer buttons
        still override until the next app-theme toggle."""
        from widgets import GenomePanel, FigurePanel
        for v in self.findChildren(GenomePanel) + self.findChildren(FigurePanel):
            v.apply_app_theme(self.theme)

    def _uppercase_buttons(self):
        """Mono uppercase action buttons (the web '.btn' look). Skips link buttons and leaves glyphs/
        symbols untouched — str.upper() only affects letters. Idempotent."""
        for b in self.findChildren(QPushButton):
            if b.property("link"):
                continue
            t = b.text()
            if t and t != t.upper():
                b.setText(t.upper())

    def _repolish(self, w):
        w.style().unpolish(w); w.style().polish(w)

    def _toggle_theme(self):
        self.theme = "light" if self.theme == "dark" else "dark"
        self._apply_theme()
        d = self.state.get("lastPrimers")                     # QC ΔG cell colours are per-theme -> re-render on toggle
        if isinstance(d, dict) and d.get("candidates"):
            self._render_primers(d)

    def _toggle_env(self):
        vis = not self.envBox.isVisible()
        self.envBox.setVisible(vis)
        self.envHdr.setText(("▾" if vis else "▸") + " ◆  Environment")

    # ---------- startup ----------
    def _startup(self):
        self.engine.submit("health")
        self.engine.submit("env")
        self._init_wsl()
        self._render_pcr_queue()

    # ---------- input handling ----------
    def _set_seq(self, text):
        """Programmatic specimen load (fetch / upload / sample). Wrapped in the _loading guard so the
        textChanged handler keeps state['seq'] in sync but does NOT treat the load as a user edit that
        would wipe the source the loader is about to set."""
        self.state["features"] = None                         # drop a prior accession's gene model (fetch re-sets it after)
        self._loading = True
        try:
            self.seq.setPlainText(text)
        finally:
            self._loading = False

    def _seq_changed(self):
        txt = self.seq.toPlainText()
        body = "".join(l for l in txt.splitlines() if not l.startswith(">"))
        self.charCount.setText(f"{len(body)} nt")
        self.state["seq"] = txt.strip()                       # specimen tracks the box: splice/annotate read this
        if not self._loading:
            if self.state.get("source") is not None:
                self.state["source"] = None                   # a genuine edit no longer matches the fetched identity
                self.accMeta.setText(""); self.coordMeta.setText("")
                self.state["features"] = None
            self._update_splice_ref()                         # refresh the 'Genomic reference' label to the edited specimen

    def _load_sample(self):
        self._set_seq(make_sample())
        self.state["source"] = None
        self.accMeta.setText(""); self.coordMeta.setText("")

    def _upload(self):
        path, _ = QFileDialog.getOpenFileName(self, "Open FASTA", "",
                                              "FASTA (*.fa *.fasta *.fna *.txt *.gz);;All files (*)")
        if not path:
            return
        try:
            if path.lower().endswith(".gz"):
                with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
                    data = f.read()
            else:
                with open(path, "r", encoding="utf-8", errors="replace") as f:
                    data = f.read()
        except OSError as e:
            return self._banner(f"could not read file: {e}")
        self._set_seq(data)
        self.state["source"] = None
        self.accMeta.setText(f"loaded {os.path.basename(path)}"); self.coordMeta.setText("")

    def _set_fetch_enabled(self, on):
        self.fetchBtn.setEnabled(on); self.coordFetchBtn.setEnabled(on)   # both disabled in-flight: no overlapping fetch race

    def _fetch(self):
        acc = self.acc.text().strip()
        if not acc:
            return self._banner("enter an accession first")
        self.accMeta.setText("fetching…")
        self._set_fetch_enabled(False)
        self.engine.submit("fetch", {"accession": acc}, key="fetch")

    def _run_analysis(self):
        seq = self.seq.toPlainText().strip()
        if not seq:
            return self._banner("paste, upload, or fetch a sequence first")
        self.state["seq"] = seq
        self.state["analyzed_seq"] = seq                  # snapshot the sequence the reported feature coords index
        self.runBtn.setEnabled(False); self.runBtn.setText("… analysing")
        self.engine.submit("analyze", {"sequence": seq, "source": self.state["source"]}, key="analyze")

    # ---------- results routing ----------
    def _on_done(self, key, res):
        if key == "health":
            self.statusTxt.setText(f"backend live · primer3 {res.get('primer3')}")
            self.led.setProperty("live", True); self._repolish(self.led)
            if res.get("core"):
                self.ver.setText("v" + res["core"])
        elif key == "env":
            self._render_env(res)
        elif key == "fetch":
            self._on_fetch(res)
        elif key == "analyze":
            self._on_analyze(res)
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

    def _reset_buttons(self, key):
        if key == "analyze":
            self.runBtn.setEnabled(True); self.runBtn.setText("▶ Run analysis")
        elif key == "primers":
            self._design_inflight = False
            self.designBtn.setEnabled(True); self.designBtn.setText("⋯ Design primers")
        elif key == "annotate":
            self.annotateBtn.setEnabled(True); self.annotateBtn.setText("▶ Run family annotation")
            # _reset_buttons is error-only (success re-enables in _on_annotate), so clear the stuck busy body here —
            # a BadRequest/fault never reaches _on_annotate, which would otherwise leave "Running RepeatMasker…" spinning
            self._set_body(self.wslBody, _empty("Run RepeatMasker against Dfam to name the TE family. Family naming is the Linux (WSL) backend."))
        elif key == "splice":
            self.spliceBtn.setEnabled(True); self.spliceBtn.setText("▶ Detect exons / introns")
            self._set_body(self.spliceBody, _empty("Align a transcript to the loaded sequence to resolve exon–intron structure de novo."))
        elif key == "fetch":
            self._set_fetch_enabled(True)                     # fetch failed — re-enable so the user can retry
            for lbl in (self.accMeta, self.coordMeta):        # clear only the in-flight indicator, keep a prior result
                if lbl.text() == "fetching…":
                    lbl.setText("")

    def _on_user_error(self, key, msg):
        if key.startswith("pcr#"):                            # a failed pair fills its slot so the batch still renders
            self._pcr_slot(key, {"error": msg, "amplicons": []})
        elif key == "genome_pcr":
            self._genome_inflight = False; self._render_genome_status("Genome scan failed — " + msg)
            self._refresh_genome_manager()                    # scan settled (failed) — re-enable a manager opened mid-scan
        elif key == "genome_prepare":
            self._genome_prep_inflight = False; self._pending_scan = None
            self._render_genome_status("Genome download failed — " + msg)
            self._refresh_genome_manager()                    # download settled (failed) — re-enable a manager opened during it
        elif key == "genome_prepare_log":
            return                                            # a background poll blip: never banner
        else:
            self._reset_buttons(key)
        self._banner(msg)

    def _on_failed(self, key, msg, trace):
        if key.startswith("pcr#"):
            self._pcr_slot(key, {"error": msg, "amplicons": []})
        elif key == "genome_pcr":                             # clear the in-flight guard + stuck status on a hard fault
            self._genome_inflight = False; self._render_genome_status("Genome scan failed — " + msg)
            self._refresh_genome_manager()                    # scan settled (failed) — re-enable a manager opened mid-scan
        elif key == "genome_prepare":
            self._genome_prep_inflight = False; self._pending_scan = None
            self._render_genome_status("Genome download failed — " + msg)
            self._refresh_genome_manager()                    # download settled (failed) — re-enable a manager opened during it
        elif key == "genome_prepare_log":
            return                                            # a background poll blip: never banner
        else:
            self._reset_buttons(key)
        sys.stderr.write(trace + "\n")
        self._banner(f"unexpected error: {msg}")

    def _banner(self, msg, level="error"):
        # level drives colour + glyph so a success/advisory message is not painted as a red error (theme.py
        # styles #errbanner[level=...]). Default stays 'error' so every existing call keeps its red styling.
        glyph = {"error": "⚠", "warn": "⚠", "success": "✓", "info": "ℹ"}.get(level, "⚠")
        self.errbanner.setText(f"{glyph} {msg}")
        self.errbanner.setProperty("level", level)
        self._repolish(self.errbanner)                    # re-evaluate the level-aware QSS
        self.errbanner.setVisible(True)

    def _clear_banner(self):
        self.errbanner.setVisible(False)

    def _render_env(self, e):
        if e.get("error"):
            self.envBox.setText(f"<span style='color:#E06A5A'>{e['error'][:60]}</span>"); return
        pkgs = "<br>".join(f"{p['name']} {'✓' if p.get('ok') else '✗ ' + str(p.get('installed','missing'))}"
                           for p in e.get("packages", []))
        st = ("install needed" if e.get("needs_install") else "up to date")
        bw = e.get("backends", {})
        self.envBox.setText(
            f"<b>state</b> {st}{' · first run' if e.get('first_run') else ''}<br>"
            f"<b>python</b> {e.get('python','?')} {'ok' if e.get('python_ok') else 'old'}<br>"
            f"<b>packages</b><br>{pkgs}<br>"
            f"<b>wsl2</b> {bw.get('wsl2','—')}<br>"
            f"<b>signature</b> {e.get('signature','—')}")

    def _toggle_coord(self):
        vis = not self.coordBox.isVisible()
        self.coordBox.setVisible(vis)
        self.coordToggle.setText("⌖ Fetch by coordinate ▴" if vis else "⌖ Fetch by coordinate ▾")

    def _fetch_coord(self):
        self._clear_banner()
        org = self.asmSel.currentData()
        custom = self.coordCustom.text().strip() if org == "__custom__" else ""
        if org == "__custom__" and not custom:
            return self._banner("enter a custom organism name or an assembly accession (e.g. GCF_000001405.40)")
        regions = self.coord.toPlainText().strip()
        if not regions:
            return self._banner("enter at least one region, e.g. chr13:33,016,423-33,066,143")
        strand = "-" if self.coordStrand.currentIndex() == 1 else "+"
        self.coordMeta.setText("fetching…")
        self._set_fetch_enabled(False)
        self.engine.submit("fetch_coords", {"regions": regions, "strand": strand,
                           "organism": "" if org == "__custom__" else org, "customQuery": custom}, key="fetch")

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
        self.state["source"] = res.get("source", {})
        self.state["features"] = None

    def _on_fetch(self, res):
        self._set_fetch_enabled(True)                         # fetch settled — allow the next one
        if not res.get("ok"):
            self.accMeta.setText(""); self.coordMeta.setText("")
            return self._banner(res.get("error", "fetch failed"))
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
            self.coordMeta.setText("")                    # clear the other specimen identity
            self.state["source"] = {k: res.get(k) for k in ("accession", "organism", "title", "length", "moltype") if res.get(k) is not None}
            self.state["features"] = res.get("features")
        # auto-fill species for WSL family annotation if present
        if hasattr(self, "wslSpecies") and org:
            self._set_species(org)
        self._update_splice_ref()

    def _on_analyze(self, res):
        self.runBtn.setEnabled(True); self.runBtn.setText("▶ Run analysis")
        self._clear_banner()
        if res.get("warning"):
            self._banner(res["warning"], level="warn")
        recs = res.get("records", [])
        if not recs:
            return
        rec = recs[0]
        self.state["last_rec"] = rec
        self.state["analyzed_clean"] = self._clean_seq(self.state.get("analyzed_seq", ""))   # snapshot, not the live box — a keystroke mid-analysis must not defeat the stale-block guard
        self._update_splice_ref()
        self.designBtn.setEnabled(True); self.designHint.setText("")
        comp = rec.get("composition", {})
        self.mLen.setText(f"{comp.get('length', 0):,}")
        self.mGC.setText(f"{comp.get('gc', 0)}%")
        self.mN.setText(f"{comp.get('n', 0)}%")
        self.mValid.setText("valid" if rec.get("valid") else "invalid")
        self.mValid.setProperty("state", "good" if rec.get("valid") else "bad"); self._repolish(self.mValid)
        self.rRecords.setText(str(len(recs)))
        self.rStruct.setText(str(len(rec.get("structural", []))))
        self.rOrf.setText(str(len(rec.get("orfs", []))))
        self._render_struct_card(rec, res)
        self._uppercase_buttons()

    def _render_struct_card(self, rec, res):
        card = self.card_struct
        card.clear_body()
        card.expand()
        cl = rec.get("classification") or {}
        banner = QFrame(); banner.setObjectName("classbn")
        bl = QVBoxLayout(banner); bl.setContentsMargins(14, 12, 14, 12); bl.setSpacing(6)
        top = QHBoxLayout(); top.setSpacing(12)
        big = QLabel(cl.get("te_class", "—")); big.setObjectName("classbig")
        top.addWidget(big)
        conf = cl.get("confidence", "")
        if conf:
            cfl = QLabel(conf.upper()); cfl.setObjectName("cf"); cfl.setProperty("level", conf)
            top.addWidget(cfl)
        top.addStretch(1)
        bl.addLayout(top)
        kls = QLabel(cl.get("superfamily", "—") + self._src_html("Wicker2007"))
        kls.setObjectName("classkls"); kls.setWordWrap(True); kls.setTextFormat(Qt.RichText)
        kls.setOpenExternalLinks(True); bl.addWidget(kls)
        if cl.get("explanation"):
            ex = QLabel(cl["explanation"]); ex.setObjectName("classexp"); ex.setWordWrap(True)
            bl.addWidget(ex)
        comp = cl.get("completeness")                         # scoped structural-completeness (Axis 2 of reliability)
        if comp:
            arch = cl.get("order") or " – ".join(comp.get("present", []))
            miss = comp.get("missing") or []
            line = (f"<b>Structural completeness:</b> {comp['tier']}  ·  {comp.get('kind','')}"
                    + (f"<br><b>Domain architecture:</b> {arch}" if arch else "")
                    + (f"  ·  not detected: {', '.join(miss)}" if miss else ""))
            cw = QLabel(line); cw.setObjectName("classexp"); cw.setTextFormat(Qt.RichText); cw.setWordWrap(True)
            bl.addWidget(cw)
            scope = QLabel(f"Domains tested: {comp.get('scope','')}. “Not detected” is relative to this profile "
                           "panel — a divergent or unmodelled domain reads as not-detected, not as element decay "
                           "(completeness after Wicker 2007 / TEsorter / LTR_retriever). The tier reports how much of "
                           "the expected architecture is present at the domain level; it is not a claim that the ORFs "
                           "are intact or that the element is transposition- or infection-competent.")
            scope.setObjectName("cardmeta"); scope.setWordWrap(True); bl.addWidget(scope)
        card.bodylay.addWidget(banner)

        # genome viewer
        model = figures.gv_tracks_from_rec(rec)
        if model["tracks"]:
            gv = GenomePanel(svg_genome, "TEagle_genome")
            gv.apply_app_theme(self.theme)                # open in the current app theme
            gv.set_model(model)
            gv.set_feature_menu(self._region_menu)
            gv.setMinimumHeight(260)
            card.bodylay.addWidget(gv)

        # retroviral transcript architecture (ERV) — the correct coding-organisation model + cis-element legend
        arch = rec.get("retroviral")
        if arch:
            has_cis = any(e["type"].startswith(("PBS", "PPT")) for e in rec.get("structural", []))
            leg = ("<span style='color:#009E73'>■</span> env exon &nbsp; "
                   "<span style='color:#B0752E'>■</span> gag–pro–pol intron (fused polyprotein)"
                   + (" &nbsp; <span style='color:#8459C4'>■</span> PBS &nbsp; "
                      "<span style='color:#2C7FB8'>■</span> PPT" if has_cis else ""))
            legw = QLabel(leg); legw.setObjectName("orient"); legw.setTextFormat(Qt.RichText); legw.setWordWrap(True)
            card.bodylay.addWidget(legw)
            note = QLabel("<b>Endogenous retrovirus — transcript architecture.</b> " + arch["note"] +
                          " For the exact splice bases, send a real env transcript to the splice-detection card. " +
                          arch.get("subsplice_note", ""))
            note.setObjectName("orient"); note.setWordWrap(True); note.setTextFormat(Qt.RichText)
            card.bodylay.addWidget(note)

        # structural table (right-click a row → copy FASTA/DNA/coords, design primer here)
        struct = rec.get("structural", [])
        if struct:
            hdr = QHBoxLayout(); hdr.addWidget(_sl("Structural evidence")); hdr.addStretch(1)
            hdr.addWidget(self._src_link("Wicker2007")); card.bodylay.addLayout(hdr)
            t = DataTable(STRUCT_COLS, GLOSS)
            t.set_rows([self._struct_row(e) for e in struct])
            t.set_row_menu(lambda r: self._struct_menu(struct[r]))
            t.setMaximumHeight(180)
            card.bodylay.addWidget(t)

        # ORFs
        orfs = rec.get("orfs", [])
        if orfs:
            card.bodylay.addWidget(_sl(f"ORFs (≥40 aa) — {len(orfs)}"))
            t = DataTable(ORF_COLS, GLOSS)
            t.set_rows([[o["strand"], o["frame"], o["start"], o["end"], o["length_aa"]] for o in orfs])
            t.set_row_menu(lambda r: self._feat_menu(orfs[r]["start"], orfs[r]["end"], orfs[r]["strand"],
                                                     f"ORF_{orfs[r]['strand']}{orfs[r]['frame']}"))
            t.setMaximumHeight(160)
            card.bodylay.addWidget(t)

        # domains (right-click → copy protein/DNA/FASTA/coords, design primer here)
        doms = rec.get("domains", [])
        if doms:
            hdr = QHBoxLayout(); hdr.addWidget(_sl("Protein domains (HMMER)")); hdr.addStretch(1)
            hdr.addWidget(self._src_link("Pfam")); card.bodylay.addLayout(hdr)
            t = DataTable(DOMAIN_COLS, GLOSS)
            t.set_rows([[d["domain"], d.get("label", ""), d.get("pfam", ""),
                         f"{d['aa'][0]}–{d['aa'][1]}", f"{d['nt'][0]}–{d['nt'][1]}",
                         d.get("score"), f"{d.get('evalue'):.1e}" if d.get("evalue") is not None else "",
                         d.get("confidence", "")]
                        for d in doms])
            t.set_row_menu(lambda r: self._feat_menu(doms[r]["nt"][0], doms[r]["nt"][1], doms[r].get("strand", "+"),
                                                     doms[r]["domain"], protein=doms[r].get("protein")))
            t.setMaximumHeight(180)
            card.bodylay.addWidget(t)
            dhint = QLabel("The last column is <b>Conf</b> (per-domain confidence). On a narrow window, scroll the table "
                           "sideways — or collapse the specimen panel (Ctrl+B / ◧) — to reach every column.")
            dhint.setObjectName("orient"); dhint.setTextFormat(Qt.RichText); dhint.setWordWrap(True)
            card.bodylay.addWidget(dhint)
            drow = QHBoxLayout(); drow.addStretch(1); drow.addWidget(_export_table_btn(t, "TEagle_domains", self))
            card.bodylay.addLayout(drow)

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
            title = "Gene model (NCBI feature table" + (" + CDS-inferred exons)" if derived else ")")
            gmlay.addWidget(_sl(title))
            if not demote and (cl.get("is_erv") or _tc.startswith(("LTR", "LINE", "retro", "DNA"))):
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
                gvg = GenomePanel(svg_genome, "TEagle_genemodel"); gvg.apply_app_theme(self.theme); gvg.set_model(gmodel)
                gvg.set_feature_menu(self._region_menu); gvg.setMinimumHeight(200)
                gmlay.addWidget(gvg)
            if demote:                                         # collapse the host-style view for an ERV
                gmbox.setVisible(False)
                gmtog = QPushButton("▸ Record's raw CDS annotation (host-style — the transcript architecture above is the correct model)")
                gmtog.setProperty("link", True)
                def _tgm(_=False, b=gmbox, t=gmtog):
                    v = not b.isVisible(); b.setVisible(v)
                    t.setText(("▾" if v else "▸") + t.text()[1:])
                gmtog.clicked.connect(_tgm)
                card.bodylay.addWidget(gmtog)
            card.bodylay.addWidget(gmbox)

        for note in rec.get("notes", []):
            n = QLabel("• " + note); n.setObjectName("orient"); n.setWordWrap(True)
            card.bodylay.addWidget(n)

        # explicit methodology — which database / consensus / parameters define each evidence layer
        mtoggle = QPushButton("ⓘ Methods & databases ▾"); mtoggle.setProperty("link", True)
        mbox = QLabel(self._methods_html()); mbox.setObjectName("orient"); mbox.setWordWrap(True)
        mbox.setTextFormat(Qt.RichText); mbox.setOpenExternalLinks(True); mbox.setVisible(False)
        def _tgl():
            v = not mbox.isVisible(); mbox.setVisible(v)
            mtoggle.setText("ⓘ Methods & databases ▴" if v else "ⓘ Methods & databases ▾")
        mtoggle.clicked.connect(_tgl)
        card.bodylay.addWidget(mtoggle); card.bodylay.addWidget(mbox)

    def _methods_html(self):
        """A plain statement of exactly what defines each evidence layer — the database/model, the consensus
        source, and the thresholds — so the annotation is never a black box."""
        p, h, w, dfam = (self._src_html(k) for k in ("Pfam", "HMMER", "Wicker2007", "Dfam"))
        return (
            "<b>Protein domains</b> — profile-HMM search (HMMER" + h + ", run in-process via pyhmmer) of the "
            "6-frame ORFs (≥ 40 aa) against a bundled Pfam-A" + p + " TE-domain profile set (21 models, all CC0): "
            "<b>POL</b> — RT PF00078/PF07727/PF13456, integrase PF00665, RNase&nbsp;H PF00075, protease PF00077; "
            "<b>GAG</b> — matrix PF02337, capsid PF00607/PF19317, nucleocapsid PF14787, retrotransposon-gag PF03732; "
            "<b>ENV</b> — envelope glycoprotein PF13804, transmembrane PF00517, surface PF00429; chromodomain PF00385; "
            "and transposases PF01498/PF03184/PF13358/PF01359/PF05699/PF14372. A hit is kept at per-domain "
            "E-value ≤ 1e-3; the gag + env models let TEagle recover the full GAG–POL–ENV architecture of ERVs (HERV-K, "
            "-W, -L …), not just the pol enzymes.<br>"
            "<b>Structural evidence</b> — heuristic terminal-repeat detectors (no external database): LTR by k-mer "
            "seed + diagonal cluster (k=13, ≥ 80 bp, ≥ 80% identity, ≥ 4 anchors); TIR by a terminal inverted-repeat "
            "scan plus a k-mer-vs-reverse-complement search; poly-A/poly-T tail ≥ 8 bp; TSD as a 4–12 bp exact "
            "flanking direct repeat. Coordinates are 0-based half-open, and each row lists its own detection method.<br>"
            "<b>Superfamily / class</b> — the Wicker&nbsp;2007" + w + " scheme, derived from the domain architecture and "
            "structural context; Copia vs Gypsy is called from the strand-aware integrase-vs-RT translation order, not "
            "ORF length; an env domain with paired LTRs flags an endogenous retrovirus (ERV).<br>"
            "<b>Reliability</b> — reported on two independent, citable axes rather than one fabricated number: (1) a "
            "per-domain call confidence from the HMMER i-Evalue (Eddy 2011); (2) a categorical structural-completeness "
            "tier — intact / near-complete / partial / structural-only — mapped to the autonomous/intact criteria of "
            "Wicker 2007, TEsorter and LTR_retriever, and always scoped to the domain models actually tested.<br>"
            "<b>Family naming</b> (optional, WSL backend) — RepeatMasker (RMBLAST) against the curated Dfam&nbsp;4.0" + dfam +
            " library; this is the only step that makes a database family call, and it is absent from the offline path.")

    def _struct_row(self, e):
        fp, tp = e.get("five_prime"), e.get("three_prime")
        if fp and tp:                                     # a terminal-repeat PAIR (LTR/TIR): show BOTH copies,
            coords = f"{fp[0]}–{fp[1]}  ·  {tp[0]}–{tp[1]}"   # matching the two blocks drawn in the genome viewer
        else:
            sp = e.get("pos") or e.get("upstream") or e.get("element_span") or [None, None]
            coords = f"{sp[0]}–{sp[1]}" if sp[0] is not None else ""
        arm = e.get("ltr_len") or e.get("tir_len") or e.get("length") or ""
        metric = (f"{e['identity']}%" if e.get("identity") is not None else e.get("motif", ""))
        return [e["type"], coords, arm, metric, e.get("method", "")]


    # =================== sequence helpers / staleness ===================
    _RC = {"A": "T", "T": "A", "G": "C", "C": "G", "U": "A", "R": "Y", "Y": "R", "S": "S", "W": "W",
           "K": "M", "M": "K", "B": "V", "D": "H", "H": "D", "V": "B", "N": "N"}

    @staticmethod
    def _revcomp(s):
        return "".join(MainWindow._RC.get(c, "N") for c in reversed(s.upper()))

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

    def _slice(self, s, e, seq=None):
        # panel-01 feature coords index the ANALYZED snapshot, not the live box (which splice/annotate submit);
        # an explicit seq (family/splice/amplicon) is used verbatim.
        base = self._norm_seq(seq) if seq is not None else self._norm_seq(self.state.get("analyzed_seq", ""))
        return base[s:e]

    def _stale_block(self):
        """Block primer / PCR when the sequence in the box differs from the analysed one."""
        cur = self._clean_seq()
        if not self.state.get("analyzed_clean"):
            self._banner("Run analysis first."); return True
        if cur and cur != self.state["analyzed_clean"]:
            self._banner("Sequence changed since analysis — Run analysis again before designing primers.")
            return True
        return False

    def _copy(self, text):
        QApplication.clipboard().setText(text)
        QToolTip.showText(QCursor.pos(), "copied", self)     # brief feedback at the cursor, not the status chip

    def _feat_menu(self, start, end, strand, label, protein=None, dna=None, src_seq=None):
        """Right-click menu for a feature. Coordinates address `src_seq` (default: the panel-01
        specimen); pass `dna` when the exact sequence is already known (an amplicon carries its own
        seq relative to its own background), so copies never re-slice the wrong template."""
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
        items = [(f"⧉ Copy FASTA", lambda: self._copy(f">{fid}_{start}-{end}{rev}\n{dna}")),
                 (f"⧉ Copy DNA", lambda: self._copy(dna)),
                 (f"⧉ Copy coords ({start}–{end} {strand})", lambda: self._copy(f"{start}-{end} {strand}"))]
        if protein:
            items.append(("⧉ Copy protein", lambda: self._copy(protein)))
        items.append(("⌖ Design primer here", _design))
        items.append(("◧ Send to splice detection",
                      lambda: self._send_to_splice(f">{fid}_{start}-{end}{rev}\n{dna}")))
        if len(dna) >= 2:                                 # let the user narrow to a sub-interval before routing
            items.append(("◫ Select a sub-region → primer / splice…", lambda: self._subregion(dna, fid)))
        return items

    def _subregion(self, dna, fid):
        """Pick a sub-interval WITHIN a feature (1-based inclusive, offsets into the feature's own sequence) and
        send only that subset to primer design or splice detection. Offset-based, so it is strand- and
        coordinate-space-safe: `dna` is already the feature's sequence in biological orientation."""
        n = len(dna)
        dlg = QDialog(self); dlg.setWindowTitle("Select a sub-region"); dlg.resize(460, 200)
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
        bp = QPushButton("⌖ Design primers"); bp.setProperty("sm", True); bp.clicked.connect(lambda: _go("primer"))
        bs = QPushButton("◧ Send to splice"); bs.setProperty("sm", True); bs.clicked.connect(lambda: _go("splice"))
        bc = QPushButton("Cancel"); bc.setProperty("sm", True); bc.clicked.connect(dlg.reject)
        brow.addWidget(bp); brow.addWidget(bs); brow.addWidget(bc); lay.addLayout(brow)
        self._uppercase_buttons()
        dlg.exec()

    def _send_to_splice(self, fasta):
        """Load a right-clicked subsequence as the transcript in the splice card and reveal it —
        mirrors 'send to in-silico PCR'. It is aligned to the loaded genomic reference."""
        self.spliceTx.setPlainText(fasta)
        self.card_splice.expand()
        try:
            self.resultsScroll.ensureWidgetVisible(self.card_splice)
        except Exception:
            pass
        self.spliceTx.setFocus()

    def _update_splice_ref(self):
        """Show which specimen splice will align a transcript against (the genomic reference)."""
        if not hasattr(self, "spliceRef"):
            return
        seq = self.state.get("seq") or ""
        if not seq:
            self.spliceRef.setText("Genomic reference: none loaded yet — load a specimen in panel 01.")
            return
        src = self.state.get("source") or {}
        n = len(self._clean_seq(seq))
        who = src.get("displayLocus") or src.get("accession") or "pasted / uploaded specimen"
        org = f" · {src.get('organism')}" if src.get("organism") else ""
        self.spliceRef.setText(f"Genomic reference: <b>{who}</b>{org} · {n:,} bp (from panel 01)")

    def _struct_menu(self, e):
        sp = e.get("element_span") or e.get("five_prime") or e.get("pos") or e.get("upstream") or [None, None]
        if sp[0] is None:
            return [("⧉ Copy type", lambda: self._copy(e["type"]))]
        return self._feat_menu(sp[0], sp[1], "+", e["type"].split(" ")[0])

    def _region_menu(self, region):
        """Right-click menu for a feature glyph in the genome viewer (copy FASTA / design primer)."""
        return self._feat_menu(region["start"], region["end"], region.get("strand", "+"),
                               region.get("label") or "feature")

    def _gel_menu(self, region):
        """Right-click menu for a gel band → copy the amplicon FASTA / coordinates."""
        a = region.get("amplicon") or {}
        pair = region.get("pair", "")
        start, end, seq = a.get("start"), a.get("end"), a.get("seq", "")
        if region.get("has_locus", True):                     # on/off only means something with a design locus
            kind = "ontarget" if a.get("on_target") else "offtarget"
        else:
            kind = "primingsite"                              # no locus -> neutral, matching the band/tooltip
        label = f"amplicon_{pair}_{start}-{end}_{a.get('length','')}bp_{kind}"
        items = [("⧉ Copy FASTA", lambda: self._copy(f">{label}\n{seq}"))]
        if start is not None and end is not None:
            items.append((f"⧉ Copy coords ({start}–{end})", lambda: self._copy(f"{start}-{end}")))
        return items

    def _src_link(self, key, url=None):
        """A small clickable 'source ↗' citation label that opens the verified DOI in a browser."""
        r = REFLINKS.get(key)
        lab = QLabel()
        if not r:
            return lab
        accent = theme_mod.ACCENT[self.theme]
        lab.setTextFormat(Qt.RichText)
        lab.setText(f'<a href="{url or r["url"]}" style="color:{accent};text-decoration:none">source ↗</a>')
        lab.setOpenExternalLinks(True)
        lab.setObjectName("srclink")
        lab.setToolTip("Source — " + r["cite"])
        return lab

    def _src_html(self, key, url=None):
        """Inline 'source ↗' anchor for RichText labels (set openExternalLinks on the label)."""
        r = REFLINKS.get(key)
        if not r:
            return ""
        accent = theme_mod.ACCENT[self.theme]
        return (f' <a href="{url or r["url"]}" style="color:{accent};text-decoration:none" '
                f'title="Source — {r["cite"]}">source ↗</a>')

    # =================== primer params ===================
    PRESETS = {
        "standard":   dict(pMin=150, pMax=500, pTm=60, pMinS=18, pMaxS=27, pNum=5, pOptS=20, pTmMin=57, pTmMax=63, pGcMin=40, pGcMax=60, pPolyX=4, pGcClamp=0),
        "qpcr":       dict(pMin=70, pMax=150, pTm=60, pMinS=18, pMaxS=24, pNum=5, pOptS=20, pTmMin=58, pTmMax=62, pGcMin=40, pGcMax=60, pPolyX=4, pGcClamp=1),
        "highspec":   dict(pMin=150, pMax=500, pTm=62, pMinS=20, pMaxS=26, pNum=8, pOptS=22, pTmMin=60, pTmMax=64, pGcMin=45, pGcMax=60, pPolyX=3, pGcClamp=2),
        "permissive": dict(pMin=100, pMax=1000, pTm=58, pMinS=17, pMaxS=30, pNum=10, pOptS=20, pTmMin=52, pTmMax=65, pGcMin=30, pGcMax=70, pPolyX=5, pGcClamp=0),
    }

    def _apply_preset(self, name):
        p = self.PRESETS.get(name)
        if not p:
            return
        for k, v in p.items():
            if k in self.pfields:
                self.pfields[k].setText(str(v))

    def _numfield(self, fid, default):
        t = self.pfields[fid].text().strip()
        if t == "":
            return default
        try:
            f = float(t)
            return int(f) if f.is_integer() else f
        except ValueError:
            return default

    def _read_primer_params(self):
        p = {"prod_min": self._numfield("pMin", 150), "prod_max": self._numfield("pMax", 500),
             "opt_tm": self._numfield("pTm", 60), "min_size": self._numfield("pMinS", 18),
             "max_size": self._numfield("pMaxS", 27), "num_return": self._numfield("pNum", 5),
             "opt_size": self._numfield("pOptS", 20), "min_tm": self._numfield("pTmMin", 57),
             "max_tm": self._numfield("pTmMax", 63), "min_gc": self._numfield("pGcMin", 40),
             "max_gc": self._numfield("pGcMax", 60), "max_poly_x": self._numfield("pPolyX", 4)}
        clamp = self._numfield("pGcClamp", 0)
        if clamp and clamp > 0:
            p["gc_clamp"] = clamp
        return p

    # =================== primer design ===================
    def _design(self):
        self._clear_banner()
        if self._design_inflight:
            return self._banner("A primer design is already running — wait for it to finish.")
        if self._stale_block():
            return
        self._design_inflight = True
        self.designBtn.setEnabled(False); self.designBtn.setText("◴ designing…")
        self._design_tmpl = self.state["seq"]             # remember the template these candidates index (for PCR on-target)
        self.engine.submit("primers", {"sequence": self.state["seq"], "params": self._read_primer_params()}, key="primers")

    def _design_for_domain(self, start, end, label, seq=None):
        self._clear_banner()
        if self._design_inflight:
            return self._banner("A primer design is already running — wait for it to finish.")
        if seq is None and self._stale_block():           # stale guard only applies to the panel-01 specimen
            return
        tmpl = seq if seq is not None else self.state.get("seq", "")
        if not self._clean_seq(tmpl):
            return self._banner("No sequence to design a primer on.")
        inc = [start, max(60, end - start)]
        self.card_primer.expand()
        self._pending_domain = label
        self._design_inflight = True
        self._design_tmpl = tmpl                          # candidates' left/right_pos index THIS template, not always panel-01
        self.engine.submit("primers", {"sequence": tmpl, "params": self._read_primer_params(),
                                        "included": inc}, key="primers")

    def _on_primers(self, d):
        self._design_inflight = False
        self.designBtn.setEnabled(True); self.designBtn.setText("⋯ Design primers")
        self.card_primer.expand()
        self._render_primers(d)
        self._uppercase_buttons()
        if d.get("provenance"):
            self._render_provenance(d["provenance"])

    def _render_primers(self, d):
        _clear_layout(self.primBody)                          # recursive: also removes the addLayout'd "Primer pairs" header row
        self.state["lastPrimers"] = d                         # keep for a theme toggle -> re-render with the new palette's flag colours
        cands = d.get("candidates", [])
        tmpl_sig = self._norm_seq(getattr(self, "_design_tmpl", self.state.get("seq", "")))
        for c in cands:                                       # tag each pair with the template its coords index
            c["_tmpl_sig"] = tmpl_sig
        self.state["candidates"] = cands
        self.pcrStageAll.setEnabled(bool(cands))
        if not cands:
            self.primBody.addWidget(_empty("No primer pair met the criteria — try the Permissive preset or widen the product range."))
            return
        srow = QHBoxLayout(); srow.addWidget(_sl(f"Primer pairs — {len(cands)}")); srow.addStretch(1)
        srow.addWidget(self._src_link("Primer3")); self.primBody.addLayout(srow)
        headers = ["ID", "Forward (5'→3')", "Reverse (5'→3')", "Product", "Tm F/R", "GC% F/R",
                   "Hairpin", "Self-dim", "Hetero", "3′-end", "Struct", "Penalty"]
        t = DataTable(headers, GLOSS)
        fc = theme_mod.FLAG[self.theme]                       # per-theme flag colours (WCAG-tuned dark vs light)
        rows, styles, tips = [], [], []
        for c in cands:
            base = [c["id"], c["left_seq"], c["right_seq"], c["product_size"],
                    f"{c['left_tm']}/{c['right_tm']}", f"{c['left_gc']}/{c['right_gc']}"]
            qc = c.get("qc") or {}
            if qc.get("ok"):
                L, R = qc["left"], qc["right"]
                hp = _metric_cell([("F", L["hairpin"]), ("R", R["hairpin"])])
                sd = _metric_cell([("F", L["self_dimer"]), ("R", R["self_dimer"])])
                het = _metric_cell([("", qc["hetero_dimer"])])
                end = _metric_cell([("", qc["end_stability"])])
                worst = qc.get("worst", "ok")
                rows.append(base + [hp[0], sd[0], het[0], end[0], worst, c["penalty"]])   # Struct col = the flag as TEXT
                styles.append([None] * 6 + [fc.get(hp[1]), fc.get(sd[1]), fc.get(het[1]), fc.get(end[1]), fc.get(worst), None])
                tips.append([None] * 6 + [hp[2], sd[2], het[2], end[2], f"worst secondary-structure flag: {worst}", None])
            else:                                             # QC unavailable for this pair (never drops the pair)
                rows.append(base + ["—", "—", "—", "—", "n/a", c["penalty"]])
                styles.append([None] * 12); tips.append([None] * 12)
        t.set_rows(rows, styles=styles, tips=tips)
        t.set_row_menu(lambda r: [("→ send to in-silico PCR", lambda rr=r: self._add_pcr_pair(cands[rr])),
                                  ("🧬 Secondary-structure detail", lambda rr=r: self._structure_detail(cands[rr])),
                                  ("⊕ Scan whole genome for off-targets", lambda rr=r: self._scan_genome(cands[rr])),
                                  ("⧉ Copy pair FASTA", lambda rr=r: self._copy(
                                      f">{cands[rr]['id']}_F\n{cands[rr]['left_seq']}\n>{cands[rr]['id']}_R\n{cands[rr]['right_seq']}"))])
        t.setMaximumHeight(220)
        self.primBody.addWidget(t)
        # visible key for the ΔG colour flags + the ‡ marker + the F/R fold (colour is never the only signal)
        legend = QLabel(f'<span style="color:{fc["caution"]}">■</span>&nbsp;caution (moderately stable; threshold varies by structure) &nbsp;&nbsp;'
                        f'<span style="color:{fc["warn"]}">■</span>&nbsp;warn (ΔG ≤ −9) &nbsp;&nbsp;'
                        f'‡&nbsp;the two engines disagree &nbsp;·&nbsp; Hairpin / Self-dim show the worst of the forward and reverse primer')
        legend.setObjectName("orient"); legend.setTextFormat(Qt.RichText); legend.setWordWrap(True)
        self.primBody.addWidget(legend)
        eng = d.get("oligoqc_engines", {})
        cross = (f"cross-checked with ViennaRNA {eng.get('viennarna_version')}"
                 if eng.get("viennarna") else "ViennaRNA cross-check unavailable — showing Primer3 only")
        note = QLabel("Structure columns are ΔG (kcal/mol) of the most stable hairpin / self-dimer / cross-dimer and 3′-end "
                      "anneal — more negative = worse. Red (warn) is ΔG ≤ −9 on any column (the IDT rule of thumb); amber "
                      "(caution) thresholds vary by structure type (hairpin ≤ −2, dimers ≤ −5, 3′-end ≤ −6, the last flagging "
                      "an abnormally stable 3′-end cross-dimer). The Struct column states the worst flag as text. 3′-end is its own axis because it "
                      f"isolates the interaction that blocks polymerase extension. Primer3 (SantaLucia 1998) is cross-checked "
                      f"against an independent engine ({cross}). Right-click a pair → “Secondary-structure detail” for both "
                      "engines side by side. Advisory — not a wet-lab guarantee.")
        note.setObjectName("orient"); note.setWordWrap(True); self.primBody.addWidget(note)
        hint = QLabel("Right-click a pair → “send to in-silico PCR”, or stage all in panel 05. "
                      "Scroll the table sideways — or collapse the specimen panel (Ctrl+B / ◧) — to see every column.")
        hint.setObjectName("orient"); hint.setWordWrap(True); self.primBody.addWidget(hint)

    def _structure_detail(self, c):
        """IDT OligoAnalyzer-style detail: both engines' ΔG for every structure of one primer pair, side by side."""
        qc = c.get("qc") or {}
        dlg = QDialog(self); dlg.setWindowTitle(f"Secondary structure — {c.get('id','pair')}")
        dlg.resize(640, 460)
        lay = QVBoxLayout(dlg)
        head = QLabel(f"<b>{c.get('id','pair')}</b>  ·  F 5′-{c['left_seq']}-3′  ·  R 5′-{c['right_seq']}-3′")
        head.setWordWrap(True); lay.addWidget(head)
        if not qc.get("ok"):
            lay.addWidget(QLabel("Secondary-structure QC is unavailable for this pair (" + str(qc.get("error", "")) + ")."))
        else:
            L, R = qc["left"], qc["right"]
            def row(name, m):
                p3 = m.get("p3"); vr = m.get("vrna")
                return [name, f"{p3:.2f}" if p3 is not None else "—",
                        f"{vr:.2f}" if vr is not None else "—", m.get("flag", "ok"), m.get("agree", "single")]
            data = [row("Hairpin (F)", L["hairpin"]), row("Hairpin (R)", R["hairpin"]),
                    row("Self-dimer (F)", L["self_dimer"]), row("Self-dimer (R)", R["self_dimer"]),
                    row("Hetero-dimer (F×R)", qc["hetero_dimer"]), row("3′-end stability", qc["end_stability"])]
            t = DataTable(["Structure", "Primer3 ΔG", "ViennaRNA ΔG", "Flag", "Engines"], GLOSS)
            fc = theme_mod.FLAG[self.theme]
            styles = [[None, None, None, fc.get(r[3]), (fc.get("warn") if r[4] == "disagree" else None)] for r in data]
            t.set_rows(data, styles=styles)
            t.setMinimumHeight(230); lay.addWidget(t)
            cd = qc.get("conditions", {}); eng = qc.get("engines", {})
            meta = QLabel(f"ΔG in kcal/mol at {cd.get('temp_c','?')} °C, {cd.get('mv_conc','?')} mM Na⁺, "
                          f"{cd.get('dv_conc','?')} mM Mg²⁺, {cd.get('dna_conc','?')} nM oligo (IDT-comparable). "
                          f"Primer3 {eng.get('primer3','?')} (SantaLucia 1998, thal) cross-checked against ViennaRNA "
                          f"{eng.get('viennarna','?')} (DNA Mathews-2004, independent). Dimer ΔG is the intermolecular "
                          "binding energy. More negative = more stable = worse; warn is ΔG ≤ −9 on any structure (IDT rule of "
                          "thumb), caution varies by type (hairpin ≤ −2, dimers ≤ −5, 3′-end ≤ −6). The 3′-end is reported as "
                          "its own axis because it isolates the interaction that blocks polymerase extension. Neither "
                          "reproduces IDT’s mfold/UNAFold numbers exactly — this is an "
                          "independent second opinion; agreement within ~1–2 kcal/mol is the useful signal.")
            meta.setObjectName("orient"); meta.setWordWrap(True); lay.addWidget(meta)
            srow = QHBoxLayout(); srow.addWidget(QLabel("Methods:"))
            for k in ("SantaLucia1998", "Owczarzy2008", "ViennaRNA"):
                srow.addWidget(self._src_link(k))
            srow.addStretch(1); lay.addLayout(srow)
        row = QHBoxLayout(); row.addStretch(1)
        copy = QPushButton("⧉ Copy"); copy.clicked.connect(lambda: self._copy(self._structure_text(c)))
        close = QPushButton("Close"); close.clicked.connect(dlg.accept)
        row.addWidget(copy); row.addWidget(close); lay.addLayout(row)
        self._uppercase_buttons()
        dlg.exec()

    @staticmethod
    def _structure_text(c):
        qc = c.get("qc") or {}
        if not qc.get("ok"):
            return f"{c.get('id','pair')}: secondary-structure QC unavailable"
        L, R = qc["left"], qc["right"]
        def line(name, m):
            return f"{name}\tprimer3 {m.get('p3')}\tViennaRNA {m.get('vrna')}\t{m.get('flag')}\t{m.get('agree')}"
        return "\n".join([f"# {c.get('id','pair')} secondary structure (ΔG kcal/mol)",
                          f"# F {c['left_seq']} / R {c['right_seq']}",
                          line("hairpin_F", L["hairpin"]), line("hairpin_R", R["hairpin"]),
                          line("selfdimer_F", L["self_dimer"]), line("selfdimer_R", R["self_dimer"]),
                          line("heterodimer", qc["hetero_dimer"]), line("end_stability", qc["end_stability"])])

    # =================== PCR queue ===================
    @staticmethod
    def _pcr_key(c):
        return c["left_seq"] + "|" + c["right_seq"]

    def _add_pcr_pair(self, c):
        pairs = self.state.setdefault("pcrPairs", [])
        if not any(self._pcr_key(p) == self._pcr_key(c) for p in pairs):
            pairs.append(c)
        self._render_pcr_queue()
        self.card_pcr.expand()

    def _pcr_stage_all(self):
        for c in self.state.get("candidates", []):
            self._add_pcr_pair(c)

    def _pcr_clear(self):
        self.state["pcrPairs"] = []
        self._render_pcr_queue()

    def _render_pcr_queue(self):
        while self.pcrQueueBox.count():
            w = self.pcrQueueBox.takeAt(0).widget()
            if w:
                w.setParent(None)
        pairs = self.state.get("pcrPairs", [])
        self.pcrClear.setEnabled(bool(pairs))
        self.runPcrBtn.setEnabled(bool(pairs))
        if not pairs:
            self.pcrQueueBox.addWidget(_empty("No pairs loaded. Design primers, then “send to in-silico PCR”, or stage all."))
            self.pcrHint.setText("load one or more pairs, then run")
            return
        for i, c in enumerate(pairs):
            roww = QWidget(); rl = QHBoxLayout(roww); rl.setContentsMargins(0, 0, 0, 0)
            lab = QLabel(f"P{i+1}  {c['left_seq'][:16]}… / {c['right_seq'][:16]}…  · {c['product_size']} bp")
            lab.setObjectName("cardmeta"); rl.addWidget(lab); rl.addStretch(1)
            up = QPushButton("↑"); up.setProperty("sm", True); up.clicked.connect(lambda _=False, k=i: self._move_pair(k, -1))
            dn = QPushButton("↓"); dn.setProperty("sm", True); dn.clicked.connect(lambda _=False, k=i: self._move_pair(k, 1))
            rm = QPushButton("✕"); rm.setProperty("sm", True); rm.clicked.connect(lambda _=False, k=i: self._remove_pair(k))
            rl.addWidget(up); rl.addWidget(dn); rl.addWidget(rm)
            self.pcrQueueBox.addWidget(roww)
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
            self._render_pcr_queue()

    # =================== run in-silico PCR ===================
    def _run_pcr(self):
        self._clear_banner()
        if self._genome_inflight or self._genome_prep_inflight:
            return self._banner("A whole-genome scan / download is running — wait for it before an in-silico PCR run.")
        if self._stale_block():
            return
        pairs = self.state.get("pcrPairs", [])
        if not pairs:
            return self._banner("Load at least one primer pair first.")
        p = {"max_mm": self._numfield("pcrMM", 2), "tp": self._numfield("pcrTP", 5),
             "prod_min": self._numfield("pcrPmin", 70), "prod_max": self._numfield("pcrPmax", 1000)}
        bg = self.pcrBg.toPlainText()
        # _tmpl_sig uses _norm_seq while the stale-gate uses _clean_seq; a post-analysis edit adding a gap/digit
        # could pass the gate yet flip this sig -> spurious OFF-target (fails safe, never a false on-target).
        cur_sig = self._norm_seq(self.state["seq"])           # on-target only when the pair was designed on THIS template
        self._pcr_gen += 1                                     # new batch id: results from a superseded batch are dropped
        self._pcr_run = {"gen": self._pcr_gen, "results": [None] * len(pairs)}
        self.runPcrBtn.setEnabled(False); self.runPcrBtn.setText("◴ running…")
        for i, c in enumerate(pairs):
            # a pair designed on a different template (family/splice/amplicon) has coords that do not index
            # state["seq"]; passing its target_span here would fabricate a false on-target call -> omit it
            ts = [c["left_pos"][0], c["right_pos"][1]] if c.get("_tmpl_sig") == cur_sig else None
            body = {"sequence": self.state["seq"], "background": bg, "fwd": c["left_seq"], "rev": c["right_seq"],
                    "target_span": ts, "params": p}
            self.engine.submit("pcr", body, key=f"pcr#{self._pcr_gen}#{i}")

    def _scan_genome(self, cand, org=None):
        """Whole-genome scan of one primer pair via LOCAL isPcr against a downloaded RefSeq assembly. Fast
        once the genome is cached; if it is not yet downloaded the backend replies need_prepare and we offer
        a one-time download. Every genome-wide amplicon is an off-target candidate (not a validated band).
        `org` is captured at scan start so a later dropdown change can't retarget an in-flight scan/download."""
        self._clear_banner()
        if self._genome_inflight or self._genome_prep_inflight:
            return self._banner("A whole-genome scan / download is already running — wait for it to finish.")
        run = getattr(self, "_pcr_run", None)
        if run and any(r is None for r in run.get("results", [])):
            return self._banner("An in-silico PCR run is still in progress — wait for it before a genome scan.")
        if org is None:
            org = self.genomeOrg.currentData()
        if not org:
            return self._banner("Select an organism in the whole-genome off-target scan card, or download one first.")
        self._pending_scan = {"cand": cand, "org": org}        # cand + org captured so a download resumes the SAME target
        self._genome_inflight = True
        self.card_genome.expand()
        self._render_genome_busy(f"Scanning the {org} genome locally (isPcr) for priming sites "
                                 "(up to ~1–2 min on a large genome)…")
        # fixed wide window (isPcr defaults) — decoupled from the local-PCR size fields the user never sets for a
        # right-click genome scan (coupling them silently capped genome products at the pcrPmax default of 1 kb)
        p = {"min_perfect": 15, "min_good": 15, "prod_max": 4000, "prod_min": 0}
        self.engine.submit("genome_pcr", {"fwd": cand["left_seq"], "rev": cand["right_seq"], "organism": org,
                                          "design_locus": self._design_locus_for(cand, org), "params": p}, key="genome_pcr")

    def _design_locus_for(self, cand, org):
        """The specimen's OWN genome locus/loci — but ONLY when (a) this candidate was designed on the CURRENT
        specimen (its _tmpl_sig matches, so state['source'] truly describes it; a stale row from an earlier
        specimen must never fabricate a false on-target) and (b) that specimen was fetched by coordinate in the
        SAME assembly being scanned. Returns a LIST of {accession,start,stop}, one per fetched region (a
        multi-region specimen's pair may sit in any region), else None -> the scan reports neutral 'genomic
        priming sites' (a pasted/consensus specimen has no genome position)."""
        if self._norm_seq(self.state.get("seq", "")) != cand.get("_tmpl_sig"):
            return None                                       # candidate not designed on the currently loaded specimen
        src = self.state.get("source") or {}
        scanned_asm = (COORD_ASSEMBLIES.get(org, {}) or {}).get("assemblyAccession", "")
        if not scanned_asm or src.get("assemblyAccession") != scanned_asm:
            return None
        loci = [{"accession": r["chrAccession"], "start": r.get("start"), "stop": r.get("stop")}
                for r in (src.get("regions") or []) if r.get("chrAccession")]
        return loci or None

    def _prepare_genome(self, org, then_scan=None):
        """Download + cache an organism's RefSeq genome for local scanning (one-time, non-blocking)."""
        if self._genome_prep_inflight or self._genome_inflight:   # also block while a scan runs (reachable via the manager)
            return self._banner("A whole-genome scan / download is already running — wait for it to finish.")
        if not org:
            return self._banner("Select an organism to download.")
        self._pending_scan = then_scan
        self._genome_prep_inflight = True
        self.card_genome.expand()
        acc = (COORD_ASSEMBLIES.get(org, {}) or {}).get("assemblyAccession", "")
        self._prep_org = org
        self._render_genome_busy(f"Downloading the {org} genome ({acc}) — one-time, kept for future scans. "
                                 "This can take several minutes (larger for mammalian genomes)…")
        self.engine.submit("genome_prepare", {"organism": org}, key="genome_prepare")
        if getattr(self, "_prep_timer", None) is None:        # poll the prepare log so a long download shows liveness
            self._prep_timer = QTimer(self); self._prep_timer.setInterval(2500)
            self._prep_timer.timeout.connect(self._poll_prepare_log)
        self._prep_timer.start()

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
            self._render_genome_status("Genome download failed — " + d.get("error", "unknown error"))
            return self._banner("Genome download failed — " + d.get("error", ""))
        mb = (d.get("bytes", 0) or 0) / 1e6
        self._banner(f"Genome ready · {d.get('organism','')} ({d.get('assemblyAccession','')}) · "
                     f"{d.get('n_seqs','?')} sequences · {mb:.0f} MB cached.", level="success")
        ps = self._pending_scan
        if ps:                                                # resume straight into the SAME scan the download was for
            self._pending_scan = None
            self._scan_genome(ps["cand"], org=ps["org"])
        else:
            self._render_genome_status(f"{d.get('organism','')} genome downloaded — right-click a pair to scan it.")
        # refresh unconditionally: the just-downloaded organism must appear in the PCR dropdown even if the
        # manager is closed. The split _on_genome_list updates the dropdown always, the manager only if open.
        self.engine.submit("genome_list", {}, key="genome_list")

    # =================== cached-genome manager ===================
    def _refresh_genome_dropdown(self):
        """Rebuild the PCR organism dropdown from the prepared (downloaded+verified) genome set only —
        an organism cannot be scanned until it is downloaded. userData stays the ORGANISM string (what
        _scan_genome and the need_prepare fallback consume via currentData()); a downloaded accession
        outside the curated map is skipped. Selection is preserved across the rebuild."""
        box = getattr(self, "genomeOrg", None)
        if box is None:                                       # genome_list can resolve before the PCR card is built
            return
        prev = box.currentData()
        acc2org = {v["assemblyAccession"]: k for k, v in COORD_ASSEMBLIES.items()}
        box.blockSignals(True)
        box.clear()
        box.addItem("— select organism —", None)              # placeholder so no wrong-species scan runs by default
        n = 0
        for g in self._prepared_genomes:
            org = acc2org.get(g.get("accession"))
            if not org:
                continue
            box.addItem(f"{org} · {COORD_ASSEMBLIES[org]['assemblyName']}", org)
            n += 1
        if prev is not None:
            i = box.findData(prev)
            if i >= 0:
                box.setCurrentIndex(i)
        box.blockSignals(False)
        hint = getattr(self, "genomeOrgHint", None)
        if hint is not None:
            hint.setText("" if n else "No genomes downloaded yet — open ⚙ Manage genomes to download one for scanning.")

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

    def _on_genome_list(self, d):
        # This handler is the SINGLE fan-out for every genome_list submit (WSL-ready, download, delete,
        # manager-open). ALWAYS update the dropdown first; the dialog rebuild below is the only part guarded
        # by manager-open — so a download finishing with the manager closed still refreshes the dropdown.
        if d.get("ok"):                                       # only overwrite on success — a transient WSL blip must not
            self._prepared_genomes = d.get("genomes", [])     # blank the dropdown while genomes are still cached on disk
        self._refresh_genome_dropdown()

        dlg = getattr(self, "_genome_mgr", None)
        if dlg is None:
            if not getattr(self, "_genome_mgr_open", False):  # a refresh that landed AFTER the user closed the manager
                return                                        # -> only skip the DIALOG; the dropdown is already updated
            dlg = QDialog(self); dlg.setAttribute(Qt.WA_DeleteOnClose)
            dlg.finished.connect(lambda *_: (setattr(self, "_genome_mgr", None),
                                             setattr(self, "_genome_mgr_open", False)))
            self._genome_mgr = dlg
        dlg.setWindowTitle("Manage genomes")
        dlg.setMinimumWidth(760)
        old = dlg.layout()
        if old is not None:                                   # rebuild the body in place on refresh
            QWidget().setLayout(old)
        lay = QVBoxLayout(dlg); lay.setSpacing(8)
        if not d.get("ok"):
            lay.addWidget(QLabel("Could not list genomes — " + d.get("error", "WSL backend unavailable")))
            close = QPushButton("Close"); close.clicked.connect(dlg.accept)
            lay.addWidget(close); dlg.show(); dlg.raise_(); return
        lay.addWidget(QLabel("Download a genome once to enable its whole-genome off-target scan. It is kept locally "
                             "and then appears in the in-silico PCR organism dropdown. Mammalian genomes are large "
                             "(~1 GB, a few minutes)."))
        prepared = {g["accession"]: g for g in d.get("genomes", [])}
        orgs = sorted(COORD_ASSEMBLIES)
        busy = self._genome_prep_inflight or self._genome_inflight   # lock every row while any download/scan runs
        headers = ["Organism", "Assembly", "Accession", "Status", "Size (MB)", "Contigs", "Action"]
        tbl = QTableWidget(len(orgs), len(headers))
        tbl.setHorizontalHeaderLabels(headers)
        tbl.verticalHeader().setVisible(False)
        tbl.setEditTriggers(QAbstractItemView.NoEditTriggers)
        tbl.setSelectionMode(QAbstractItemView.NoSelection)
        tbl.setSortingEnabled(False)                          # static catalog; sorting conflicts with setCellWidget buttons
        tbl.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        tbl.horizontalHeader().setStretchLastSection(True)
        for r, org in enumerate(orgs):
            meta = COORD_ASSEMBLIES[org]
            acc = meta["assemblyAccession"]
            g = prepared.get(acc)
            if g:
                cells = [org, meta["assemblyName"], acc, "Downloaded ✓",
                         f"{(g.get('bytes', 0) or 0) / 1e6:.0f}", str(g.get("n_seqs", "?"))]
            else:
                cells = [org, meta["assemblyName"], acc, "not downloaded", "—", "—"]
            for c, txt in enumerate(cells):
                it = QTableWidgetItem(txt); it.setTextAlignment(Qt.AlignCenter)
                tbl.setItem(r, c, it)
            btn = QPushButton("🗑 Delete" if g else "⭳ Download"); btn.setProperty("sm", True)
            if g:
                btn.clicked.connect(lambda _=False, a=acc: self._remove_genome(a))
            else:
                btn.clicked.connect(lambda _=False, o=org: self._on_manager_download(o))
            btn.setEnabled(not busy)
            tbl.setCellWidget(r, 6, btn)
        tbl.resizeColumnsToContents()
        lay.addWidget(tbl)
        self._genome_mgr_table = tbl
        if busy:
            lay.addWidget(QLabel("A download or scan is running — actions are disabled until it finishes."))
        close = QPushButton("Close"); close.clicked.connect(dlg.accept)
        lay.addWidget(close)
        dlg.show(); dlg.raise_()

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

    def _remove_genome(self, accession):
        self.engine.submit("genome_remove", {"assemblyAccession": accession}, key="genome_remove")

    def _on_genome_remove(self, d):
        # refresh BOTH the dropdown and the manager (if open); the split _on_genome_list self-selects what to update
        self.engine.submit("genome_list", {}, key="genome_list")

    def _render_genome_status(self, msg):
        _clear_layout(self.genomeBody)
        self.genomeBody.addWidget(_empty(msg))

    def _render_genome_busy(self, text):
        """Long-op liveness in the whole-genome scan panel: an animated indeterminate bar (download / scan),
        so a multi-minute WSL call visibly reads as working, not hung. set_text() updates it from the log poll."""
        _clear_layout(self.genomeBody)
        self._genome_busy = BusyBar(text)
        self.genomeBody.addWidget(self._genome_busy)

    def _on_genome_pcr(self, d):
        self._genome_inflight = False
        self._refresh_genome_manager()                        # scan settled — re-enable a manager opened mid-scan
        if not d.get("ok"):
            if d.get("need_prepare"):                         # genome not downloaded yet — offer the one-time download
                ps = self._pending_scan or {}
                org = ps.get("org") or self.genomeOrg.currentData()   # the org that STARTED this scan, not the live dropdown
                acc = (COORD_ASSEMBLIES.get(org, {}) or {}).get("assemblyAccession", "")
                self._render_genome_status(f"The {org} genome ({acc}) is not downloaded yet.")
                box = QMessageBox(self)
                box.setWindowTitle("Download genome?")
                box.setText(f"The {org} genome ({acc}) is not on this machine yet.\n\n"
                            "Download it once now? It is kept locally so future scans are fast. "
                            "Mammalian genomes are large (~1 GB download, a few minutes).")
                box.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
                box.setDefaultButton(QMessageBox.Yes)
                if org and box.exec() == QMessageBox.Yes:
                    self._prepare_genome(org, then_scan=ps)
                else:
                    self._pending_scan = None
                    self._render_genome_status("Genome scan cancelled — no genome downloaded.")
                return
            self._render_genome_status("No genome scan result — " + d.get("error", "scan failed"))
            return self._banner(d.get("error", "genome scan failed"))
        self._pending_scan = None
        amps = [{**a, "pair": "genome", "remote": True} for a in d.get("amplicons", [])]
        summary = d.get("summary") or {}
        has_locus = bool(summary.get("has_locus"))        # gel: on/off colours with a locus; neutral 'priming site' colour without
        lanes = [{"label": "genome", "amplicons": d.get("amplicons", []), "has_locus": has_locus, "advisory": not has_locus}]
        self.state["lastGenomeScan"] = {"lanes": lanes, "amplicons": amps, "summary": summary}
        self.card_genome.expand()
        self._render_genome_scan(lanes, amps, summary)
        m = d.get("provenance", {})
        db = (m.get("databases") or [{}])[0]
        seal = m.get("manifestSha256", "")
        prov = (f" · assembly {d.get('assemblyAccession','')}"
                + (f" · sha256 {db.get('sha256','')[:12]}…" if db.get("sha256") else "")
                + (f" · seal {seal[:12]}…" if seal else ""))
        self._banner(f"Whole-genome scan · {d.get('organism','')} ({d.get('assemblyName','')}): "
                     f"{summary.get('verdict', 'scan complete')} — a floor at ≥15 bp 3′-perfect match "
                     f"(more-diverged copies are not counted){prov}.", level="info")
        if m:
            self._render_provenance(m)

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
        self.runPcrBtn.setEnabled(True); self.runPcrBtn.setText("▶ Run loaded pairs")
        self.card_pcr.expand()
        results = run["results"]
        lanes, amps, provs = [], [], []
        for idx, dd in enumerate(results):
            lane = f"P{idx+1}"
            lanes.append({"label": lane, "amplicons": dd.get("amplicons", [])})
            for a in dd.get("amplicons", []):
                amps.append({**a, "pair": lane})
            if dd.get("provenance"):
                provs.append(dd["provenance"])
        self.state["lastPcr"] = {"lanes": lanes, "amplicons": amps}
        self._render_pcr(lanes, amps)
        self._uppercase_buttons()
        if provs:
            self._render_provenance(provs[0])

    def _render_pcr(self, lanes, amps):
        _clear_layout(self.pcrBody)
        gel = FigurePanel(lambda bg, L=lanes: figures.svg_gel({"lanes": L}, bg),
                          "TEagle_gel", modes=("dark", "white", "uv", "mono"),
                          hit_regions=figures.gel_regions({"lanes": lanes}), on_menu=self._gel_menu)
        gel.apply_app_theme(self.theme)                   # gel opens in the current app theme (uv/mono stay manual)
        gel.setMinimumHeight(420)
        self.pcrBody.addWidget(gel)
        if amps:
            headers = ["Pair", "Source", "Coords", "Len", "Mism F/R", "Call"]
            rows = [[a["pair"], a["source"], f"{a['start']}–{a['end']}", a["length"],
                     f"{a.get('fwd_mm', '—')}/{a.get('rev_mm', '—')}",
                     ("on-target" if a.get("on_target") else "off-target")
                     + (" · single-primer" if a.get("single_primer") else "")] for a in amps]
            t = DataTable(headers, GLOSS)
            t.set_rows(rows)
            t.set_row_menu(lambda r: self._feat_menu(amps[r]["start"], amps[r]["end"], "+",
                                                     f"amplicon_{amps[r]['pair']}", dna=amps[r].get("seq", "")))
            t.setMaximumHeight(200)
            self.pcrBody.addWidget(t)
            if any(a.get("seq") for a in amps):
                def _amps_fasta(_=False, aa=amps, tbl=t):
                    order = [tbl._orig(r) for r in range(tbl.rowCount())]      # follow the table's current (sorted) order
                    seq_amps = [aa[i] for i in order if 0 <= i < len(aa)] or aa
                    fasta = "\n".join(
                        f">amplicon_{a['pair']}_{a['start']}-{a['end']}_{a['length']}bp_"
                        f"{'on' if a.get('on_target') else 'off'}target"
                        f"{'_singleprimer' if a.get('single_primer') else ''}\n{a.get('seq','')}" for a in seq_amps)
                    widgets.save_fasta(fasta, "TEagle_amplicons", self)
                cp = QPushButton("⭳ Export amplicons → FASTA"); cp.setProperty("sm", True)
                cp.clicked.connect(_amps_fasta)
                self.pcrBody.addWidget(cp)
        else:
            self.pcrBody.addWidget(_empty("No amplicon predicted for any pair under the criteria."))
        onN = sum(1 for a in amps if a.get("on_target"))
        spN = sum(1 for a in amps if a.get("single_primer"))
        offN = len(amps) - onN - spN
        parts = f"{len(lanes)} lane(s) · {onN} on-target · {offN} off-target" + (f" · {spN} single-primer" if spN else "")
        note = QLabel(f"{parts} product(s) · ladder lane “L”. Products of equal size co-migrate into one band (a real "
                      "gel cannot separate them), so a band can carry more than one product — a band carrying an "
                      "off-target alongside the on-target is drawn in the off-target colour (not a clean on-target). "
                      "Every product is listed in the table above. Intensity tracks priming efficiency. Not a claim of "
                      "experimental specificity.")
        note.setObjectName("orient"); note.setWordWrap(True); self.pcrBody.addWidget(note)

    def _render_genome_scan(self, lanes, amps, summary):
        """Render a whole-genome off-target scan into its own panel: gel + interpretation (verdict,
        per-chromosome spread, size cluster) + a FULL match table, ON-TARGET FIRST. On-target = the product
        at the specimen's own genome locus (when it sits in the scanned assembly); the rest are off-targets,
        or — with no design locus — neutral genomic priming sites."""
        _clear_layout(self.genomeBody)
        gel = FigurePanel(lambda bg, L=lanes: figures.svg_gel({"lanes": L}, bg),
                          "TEagle_genome_gel", modes=("dark", "white", "uv", "mono"),
                          hit_regions=figures.gel_regions({"lanes": lanes}), on_menu=self._gel_menu)
        gel.apply_app_theme(self.theme); gel.setMinimumHeight(420)
        self.genomeBody.addWidget(gel)
        verdict = QLabel("▸ " + summary.get("verdict", ""))
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
        def _call(a):
            if a.get("single_primer"):
                return "single-primer"
            if a.get("on_target"):
                return "on-target"
            return "off-target" if has_locus else "priming site"
        order = sorted(range(len(amps)), key=lambda i: (amps[i].get("single_primer", False),
                                                        not amps[i].get("on_target", False),
                                                        amps[i]["source"], amps[i]["start"]))
        samps = [amps[i] for i in order]
        if samps:
            t = DataTable(["Call", "Source", "Coords", "Len", "Strand"], GLOSS)
            t.set_rows([[_call(a), a["source"], f"{a['start']}–{a['end']}", a["length"], a.get("strand", "?")] for a in samps])
            t.set_row_menu(lambda r: [("⧉ Copy locus", lambda a=samps[r]: self._copy(f"{a['source']}:{a['start']}-{a['end']}"))])
            t.setMaximumHeight(260)
            self.genomeBody.addWidget(t)
        else:
            self.genomeBody.addWidget(_empty("No genome-wide product for this pair under the criteria."))
        n_on, n_off, n_single = summary.get("n_on", 0), summary.get("n_off", 0), summary.get("n_single", 0)
        head = f"{n_on} on-target + {n_off} off-target site(s)" if has_locus else f"{summary.get('n_pair', 0)} genomic priming site(s)"
        sp = f" · {n_single} single-primer artefact(s)" if n_single else ""
        comig = (" An on-target that shares a band size with an off-target is drawn in the off-target colour — a "
                 "co-migration cannot be resolved on a gel — with the full split kept in the table above." if has_locus else "")
        note = QLabel(f"{head}{sp}, listed above (on-target first).{comig} Every product is a candidate under isPcr's ≥15 bp "
                      "3′-perfect rule — a specificity screen, not wet-lab-validated bands. The verdict is a heuristic "
                      "read of the count/spread; the numbers carry the claim.")
        note.setObjectName("orient"); note.setWordWrap(True); self.genomeBody.addWidget(note)

    # =================== WSL family annotation ===================
    def _init_wsl(self):
        self.engine.submit("wsl_status", key="wsl_status")

    def _on_wsl_status(self, w):
        if w.get("error"):
            self.wslStatus.setText(f"<span style='color:#E06A5A'>WSL status error: {w['error']}</span>"); return
        if not w.get("wsl2"):
            self.wslStatus.setText("<b>WSL2 not installed</b> — this optional step names the Dfam family. "
                                   "The domain-based superfamily above works without it. Install it in one click with "
                                   "<b>⚙ Backend installer</b> below (it runs the elevated <code>wsl --install</code> for you; "
                                   "a Windows restart may be required), or run <code>wsl --install</code> in an Administrator PowerShell.")
            self.spliceStatus.setText("<b>WSL2 not installed</b> — de-novo splice detection needs it (optional).")
            return
        self.engine.submit("genome_list", {}, key="genome_list")   # WSL is up — populate the downloaded-genome dropdown
        if w.get("ready"):
            self.wslStatus.setText(f"<span style='color:#33D6B8'>● ready</span> · RepeatMasker {w.get('repeatmasker')} "
                                   f"· Dfam curated · distro {w.get('distro')}")
            self.annotateBtn.setEnabled(True)
        else:
            self.wslStatus.setText(f"WSL2 ok ({w.get('distro')}); annotation stack not installed "
                                   f"(RepeatMasker {w.get('repeatmasker') or 'missing'}, Dfam {'ok' if w.get('dfam') else 'missing'}).")
            self.wslInstallBtn.setVisible(True)
        if w.get("minimap2"):
            self.spliceStatus.setText(f"<span style='color:#33D6B8'>● ready</span> · minimap2 {w.get('minimap2')} "
                                      "· align a transcript to resolve exon–intron structure")
            self.spliceBtn.setEnabled(True)
        elif w.get("wsl2"):
            self.spliceStatus.setText("minimap2 not installed in the WSL backend — it ships with the managed install (panel 03).")

    def _on_species_changed(self):
        """Show the free-text field only for 'Other…'; clear it when switching to a listed organism
        so a stale lineage can't silently drive the next annotation."""
        other = self.wslSpecies.currentData() == "__other__"
        self.wslSpeciesOther.setVisible(other)
        if not other:
            self.wslSpeciesOther.clear()

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
        for i in range(self.wslSpecies.count()):
            d = self.wslSpecies.itemData(i)
            if isinstance(d, str) and d != "__other__" and d.lower() == name.lower():
                self.wslSpecies.setCurrentIndex(i); return
        oi = self.wslSpecies.findData("__other__")        # not a listed organism -> Other + free text
        if oi >= 0:
            self.wslSpecies.setCurrentIndex(oi)
            self.wslSpeciesOther.setText(name)

    def _annotate(self):
        self._clear_banner()                                  # drop any stale error before a retry (matches _design/_run_pcr)
        if self.wslSource.currentIndex() == 1:
            seq = self.wslPaste.toPlainText().strip()
            src = None
        else:
            seq = self.state.get("seq") or self.seq.toPlainText().strip()
            src = self.state.get("source")
        if not seq:
            return self._banner("no sequence to annotate — load a specimen or paste one")
        self.state["family_seq"] = self._norm_seq(seq)    # hit coords index THIS sequence (backend-normalized), not always panel-01
        self.annotateBtn.setEnabled(False); self.annotateBtn.setText("◴ annotating…")
        self._set_body(self.wslBody, BusyBar("Running RepeatMasker against Dfam — this can take a minute or two…"))
        self.engine.submit("annotate", {"sequence": seq, "species": self._species(),
                                        "source": src}, key="annotate")

    def _on_annotate(self, d):
        self.annotateBtn.setEnabled(True); self.annotateBtn.setText("▶ Run family annotation")
        self.card_wsl.expand()
        if not d.get("ok"):
            self._set_body(self.wslBody, _empty(d.get("error", "annotation failed")))
            return
        self._render_family(d)
        if d.get("provenance"):
            self._render_provenance(d["provenance"])

    _LC = {"Low_complexity", "Simple_repeat", "Satellite", "Unknown", "Unspecified"}

    def _render_family(self, d):
        te = [h for h in d.get("hits", []) if h["class_family"] not in self._LC]
        lc = len(d.get("hits", [])) - len(te)
        if not te:
            no_species = not d.get("species") or str(d.get("species")).startswith("(all")
            hint = (" Set the organism/species above — RepeatMasker needs a lineage." if no_species
                    else " The family may need an additional Dfam taxon partition.")
            msg = f"No TE family named under the current criteria" + (f" ({lc} low-complexity region(s) found)" if lc else "") + "." + hint
            self._set_body(self.wslBody, _empty(msg))
            return
        self.state["family"] = te
        head = QLabel(f"<b>Dfam · {te[0]['class_family']}</b> — {' · '.join(sorted({h['family'] for h in te}))} "
                      f"· Dfam 4.0 curated{self._src_html('Dfam')} · RepeatMasker {d.get('repeatmasker_version','')}"
                      f"{self._src_html('RepeatMasker')} · species: {d.get('species','')}")
        head.setTextFormat(Qt.RichText); head.setWordWrap(True); head.setOpenExternalLinks(True)
        headers = ["#", "Class/family", "Dfam family", "Coords (0-based)", "Str", "Div", "Score"]
        t = DataTable(headers, GLOSS)
        t.set_rows([[i + 1, h["class_family"], h["family"], f"{h['q_start']}–{h['q_end']}",
                     h["strand"], f"{h['divergence']}%", h["score"]] for i, h in enumerate(te)])
        t.set_row_menu(lambda r: self._feat_menu(te[r]["q_start"], te[r]["q_end"], te[r]["strand"],
                                                 te[r]["family"], src_seq=self.state.get("family_seq")))
        cont = QWidget(); cl = QVBoxLayout(cont); cl.setContentsMargins(0, 0, 0, 0); cl.setSpacing(6)
        cl.addWidget(head); cl.addWidget(t)
        frow = QHBoxLayout(); frow.addStretch(1); frow.addWidget(_export_table_btn(t, "TEagle_family", self))
        cl.addLayout(frow)
        self._set_body(self.wslBody, cont)

    def _open_installer(self):
        """Open the dedicated component-wise installer dialog (per-package status, repair, integrity).
        Re-probe WSL status when it closes so the panel reflects any newly-installed backend."""
        from install_dialog import InstallDialog
        dlg = InstallDialog(self)
        dlg.finished.connect(lambda _=0: self._init_wsl())
        dlg.exec()

    # =================== splice detection ===================
    def _splice(self):
        self._clear_banner()                                  # drop any stale error before a retry (matches _design/_run_pcr)
        genomic = self.state.get("seq") or self.seq.toPlainText().strip()
        tx = self.spliceTx.toPlainText().strip()
        if not genomic.strip():
            return self._banner("Load a genomic sequence first (fetch, upload, or paste, then Run analysis).")
        if not tx:
            return self._banner("Paste a transcript / cDNA / mRNA to align.")
        self.state["splice_seq"] = self._norm_seq(genomic)    # intron menus slice THIS (backend-normalized) sequence, not a later box edit
        self.spliceBtn.setEnabled(False); self.spliceBtn.setText("◴ aligning…")
        self._set_body(self.spliceBody, BusyBar("Aligning the transcript to the genomic sequence (minimap2 -x splice)…"))
        self.engine.submit("splice", {"sequence": genomic, "transcript": tx, "source": self.state.get("source"),
                                      "timeout": 300}, key="splice")

    def _on_splice(self, d):
        self.spliceBtn.setEnabled(True); self.spliceBtn.setText("▶ Detect exons / introns")
        self.card_splice.expand()
        if not d.get("ok"):
            self._set_body(self.spliceBody, _empty(d.get("error", "splice alignment failed")))
            return
        self._render_splice(d)
        if d.get("provenance"):
            self._render_provenance(d["provenance"])

    def _render_splice(self, d):
        self.state["splice"] = d
        cont = QWidget(); cl = QVBoxLayout(cont); cl.setContentsMargins(0, 0, 0, 0); cl.setSpacing(6)
        head = QLabel(f"<b>{d['counts']['exons']} exon(s) · {d['counts']['introns']} intron(s)</b> — de novo · "
                      f"{d.get('canonical_introns',0)}/{d['counts']['introns']} canonical splice site(s) · strand {d.get('strand','')}"
                      f"{self._src_html('minimap2')}")
        head.setTextFormat(Qt.RichText); head.setWordWrap(True); head.setOpenExternalLinks(True); cl.addWidget(head)
        # independent cross-check: does this de-novo (external-transcript) alignment agree with the record's
        # own annotation? Only meaningful when the genomic IS the fetched record (source accession present).
        src, feats = self.state.get("source"), self.state.get("features")
        if src and src.get("accession") and isinstance(feats, dict) and feats.get("introns"):
            cc = cross_check_models(complete_gene_model(feats).get("introns", []), d.get("introns", []))
            if cc["annotation_total"]:
                extra = ((f" · {len(cc['aligned_only'])} alignment-only" if cc["aligned_only"] else "")
                         + (f" · {len(cc['annotation_only'])} annotation-only" if cc["annotation_only"] else ""))
                note = QLabel(f"<b>Cross-check vs {src['accession']} annotation (advisory):</b> "
                              f"{cc['matched']}/{cc['annotation_total']} intron(s) confirmed{extra}. "
                              "Independent comparison of this alignment against the record's current annotation — "
                              "not part of the sealed result (the annotation may be revised).")
                note.setTextFormat(Qt.RichText); note.setWordWrap(True); note.setObjectName("orient"); cl.addWidget(note)
                if cc["matched"] == 0:                        # a genomic slice pasted as the 'transcript' aligns gaplessly
                    hint = QLabel("No introns matched — if you aligned a genomic subsequence rather than a spliced "
                                  "transcript / mRNA / cDNA, that is expected (a genomic slice has no spliced-out introns).")
                    hint.setWordWrap(True); hint.setObjectName("orient"); cl.addWidget(hint)
        length = (self.state.get("last_rec") or {}).get("composition", {}).get("length") or len(self._clean_seq(self.state.get("seq", ""))) or 1
        model = figures.gv_tracks_from_gene({"exons": d.get("exons", []), "introns": d.get("introns", []), "cds": []}, length)
        if model["tracks"]:
            gv = GenomePanel(svg_genome, "TEagle_splice"); gv.apply_app_theme(self.theme); gv.set_model(model)
            gv.set_feature_menu(self._region_menu); gv.setMinimumHeight(220)
            cl.addWidget(gv)
        introns = d.get("introns", [])
        if introns:
            headers = ["#", "Intron span (0-based)", "Len", "Splice site", "Canonical"]
            t = DataTable(headers, GLOSS)
            t.set_rows([[k + 1, f"{i['start']}–{i['end']}", i["end"] - i["start"], f"{i['donor']}…{i['acceptor']}",
                         "canonical" if i.get("canonical") else "non-canonical"] for k, i in enumerate(introns)])
            t.set_row_menu(lambda r: self._feat_menu(introns[r]["start"], introns[r]["end"], d.get("strand", "+"),
                                                     f"intron_{r+1}", src_seq=self.state.get("splice_seq")))
            cl.addWidget(t)
        else:
            cl.addWidget(_empty("single exon — no introns detected"))
            # foot-gun guard for the common novice case (works for any input, not only fetched+annotated records):
            # a gapless alignment reads as a real single-exon finding but usually means a genomic slice was pasted.
            slice_note = QLabel("A gapless alignment (0 introns) is consistent with either a genuine single-exon "
                                "transcript OR a genomic slice pasted into the transcript box (a genomic slice has no "
                                "spliced-out introns). To resolve splicing, align an mRNA / cDNA / EST, not genomic DNA.")
            slice_note.setWordWrap(True); slice_note.setObjectName("orient"); cl.addWidget(slice_note)
        self._set_body(self.spliceBody, cont)

    # =================== provenance ===================
    def _render_provenance(self, m):
        self.card_prov.expand()
        self.card_prov.clear_body()
        inp = m.get("input", {})
        sw = "<br>".join(f"{s['name']} · {s['version']}" for s in m.get("software", []))
        pr = "<br>".join(f"{k} · {'—' if v is None else v}" for k, v in (m.get("parameters") or {}).items()) or "defaults"
        db = "<br>".join(f"{d.get('name','—')} · {d.get('version') or (d.get('sha256','')[:12]+'…' if d.get('sha256') else d.get('file','—'))}"
                         for d in m.get("databases", []))
        env = m.get("environment", {})
        nr = "<br>".join("✕ " + n for n in m.get("notRun", []))
        refs = "<br>".join(f"<b>{r['name']}</b> — {r['citation']}" + (f" doi:{r['doi']}" if r.get('doi') else "")
                           for r in m.get("references", []))
        html = (f"<b>Input</b><br>id · {inp.get('id','')}<br>length · {inp.get('length','')} bp<br>"
                f"sha256 · {str(inp.get('sha256',''))[:16]}…<br>run type · {m.get('runType','')}<br><br>"
                f"<b>Software</b><br>{sw}<br><br><b>Parameters</b><br>{pr}<br><br>")
        if db:
            html += f"<b>Databases</b><br>{db}<br><br>"
        html += (f"<b>Environment</b><br>os · {str(env.get('os',''))[:40]}<br>python · {env.get('python','')}<br>"
                 f"manifest · {str(m.get('manifestSha256',''))[:14]}…")
        if nr:
            html += f"<br><br><b>Not run</b><br>{nr}"
        if refs:
            html += f"<br><br><b>References (source-verified)</b><br>{refs}"
        lab = QLabel(html); lab.setTextFormat(Qt.RichText); lab.setWordWrap(True); lab.setObjectName("cardmeta")
        self.card_prov.bodylay.addWidget(lab)

    def _set_body(self, layout, widget):
        _clear_layout(layout)                                 # recursive: also removes any addLayout'd sub-layouts
        layout.addWidget(widget)

    def closeEvent(self, e):
        it = getattr(self, "_prep_timer", None)               # the genome-download poll timer (the only QTimer here)
        if it is not None:
            it.stop()
        super().closeEvent(e)


def selftest():
    """Bundle self-test (TEAGLE_SELFTEST=1). Proves the packaged build imports the compiled
    scientific stack (pyhmmer, primer3), ships the HMM profiles, renders figures through QtSvg,
    and runs an end-to-end analysis — the checks a double-click launch cannot report. Exit 0/1."""
    app = QApplication.instance() or QApplication([])
    import engine
    from teagle_core import domains, primers
    import figures
    from PySide6.QtSvg import QSvgRenderer
    from PySide6.QtGui import QImage, QPainter
    from PySide6.QtCore import QByteArray, Qt
    problems = []
    missing_fonts = fonts.load_fonts()                    # proves the bundled TTFs actually shipped (offscreen still loads explicit fonts)
    if missing_fonts:
        problems.append("bundled fonts missing from build: " + ", ".join(missing_fonts))
    if primers.PRIMER3_VERSION == "unavailable":
        problems.append(f"primer3 failed to load ({primers.PRIMER3_ERROR})")
    from teagle_core import oligoqc                          # the advertised secondary-structure cross-check engine
    if not oligoqc.available().get("viennarna"):
        problems.append(f"ViennaRNA cross-check engine missing from bundle ({oligoqc.available().get('viennarna_error')})")
    if domains.PYHMMER_VERSION == "unavailable":
        problems.append("pyhmmer failed to load")
    if not domains.HMM_SHA256:
        problems.append("bundled Pfam HMM profiles not found")
    else:                                                 # the profile set must LOAD and cover the full GAG-POL-ENV panel
        try:
            n_hmm = len(domains._hmms())
            codes = {v[0] for v in domains.DOMAIN_INFO.values()}
            miss = {"GAG", "PR", "RT", "RNaseH", "INT", "ENV", "CHR", "TPase"} - codes
            if miss:
                problems.append(f"domain profile set missing expected codes: {sorted(miss)}")
            if n_hmm < 21:
                problems.append(f"domain profile set loaded {n_hmm} profiles, expected >= 21 (gag/env models may be missing)")
        except Exception as e:
            problems.append(f"domain profile set failed to load: {type(e).__name__}: {e}")
    # end-to-end science through the shared engine (also exercises the fixture-free sample)
    try:
        r = engine.run_analyze({"sequence": make_sample()})
        if not r.get("records"):
            problems.append("analyze produced no records")
    except Exception as e:
        problems.append(f"analyze crashed: {type(e).__name__}: {e}")
    # QtSvg must render (the figure layer's single point of failure in a frozen build)
    svg = figures.svg_gel({"lanes": [{"label": "P1", "amplicons": [{"length": 200, "on_target": True, "source": "x"}]}]}, "dark")
    rd = QSvgRenderer(QByteArray(svg.encode("utf-8")))
    img = QImage(200, 150, QImage.Format_ARGB32); img.fill(Qt.transparent)
    p = QPainter(img); rd.render(p); p.end()
    if not (rd.isValid() and any(img.pixelColor(x, y).alpha() > 0 for x in range(0, 200, 10) for y in range(0, 150, 10))):
        problems.append("QtSvg did not render (figure plugin missing from bundle)")
    # XLSX table export must be importable AND functional (save() pulls openpyxl's lazy writer submodules)
    try:
        import io
        if not widgets._HAS_XLSX:
            problems.append("openpyxl (XLSX table export) missing from bundle")
        else:
            wb = widgets._XlWorkbook(); wb.active.append(["h", 1]); wb.save(io.BytesIO())
    except Exception as e:
        problems.append(f"XLSX export self-check failed: {type(e).__name__}: {e}")
    # the installer dialog must construct offscreen (it ships in the frozen build)
    try:
        from install_dialog import InstallDialog
        dlg = InstallDialog()
        dlg._render_components({"wsl2": True, "installing": False, "ready": False, "disk_free_gb": "50",
                                "components": [{"key": "micromamba", "name": "micromamba", "desc": "x",
                                                "ok": False, "detail": "missing", "repairable": True}]})
        if "micromamba" not in dlg._rows:
            problems.append("install dialog did not build component rows")
    except Exception as e:
        problems.append(f"install dialog crashed: {type(e).__name__}: {e}")
    if problems:
        sys.stderr.write("TEAGLE SELFTEST FAILED:\n  - " + "\n  - ".join(problems) + "\n")
        return 1
    print(f"TEAGLE SELFTEST OK · primer3 {primers.PRIMER3_VERSION} · ViennaRNA {oligoqc.VIENNARNA_VERSION} "
          f"· pyhmmer {domains.PYHMMER_VERSION} · HMM {domains.HMM_SHA256[:12]} ({len(domains._hmms())} profiles) "
          f"· QtSvg ok · install dialog ok")
    return 0


UI_SCALES = [0.75, 0.85, 1.0, 1.1, 1.25, 1.5]     # user-selectable global UI scale (persisted, applied at startup)


def _apply_saved_ui_scale():
    """Apply the persisted global UI scale via QT_SCALE_FACTOR BEFORE the QApplication is created, so a user
    on a small screen can shrink the whole UI to fit. An explicit env override always wins."""
    if os.environ.get("QT_SCALE_FACTOR") or os.environ.get("TEAGLE_UI_SCALE"):
        return
    try:
        f = float(QSettings("TEagle", "TEagle").value("ui_scale", 1.0))
        if 0.5 <= f <= 3.0 and abs(f - 1.0) > 1e-3:
            os.environ["QT_SCALE_FACTOR"] = f"{f:.2f}"
    except Exception:
        pass


def main():
    if os.environ.get("TEAGLE_SELFTEST"):
        return selftest()
    _apply_saved_ui_scale()                               # must precede QApplication creation
    if sys.platform == "win32":                           # taskbar groups under our icon, not pythonw's, in dev runs
        try:
            import ctypes
            ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID("TEagle.desktop.2")
        except Exception:
            pass
    app = QApplication.instance() or QApplication(sys.argv)
    app.setApplicationName("TEagle")                      # title stays exactly "TEagle" (no auto-appended display name)
    fonts.load_fonts()                                    # bundled Cascadia Mono (UI) — no dependence on installed fonts
    app.setWindowIcon(_app_icon())
    w = MainWindow()
    w.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
