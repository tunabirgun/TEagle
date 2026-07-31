"""TEagle native desktop app (PySide6). QMainWindow shell with a specimen rail and a scrollable
column of collapsible result cards. All science runs in-process through the shared engine, off the
GUI thread via engine_worker.Engine. This module wires the analyze workflow, primer/PCR, WSL family
annotation, splice detection, provenance and exports."""
from __future__ import annotations
import gzip, json, os, re, sys

from PySide6.QtCore import Qt, QTimer, QByteArray, QSettings
from PySide6.QtGui import (QGuiApplication, QFont, QPixmap, QPainter, QIcon, QCursor, QShortcut,
                           QKeySequence, QColor)
from PySide6.QtSvg import QSvgRenderer
from PySide6.QtWidgets import (QApplication, QMainWindow, QWidget, QFrame, QVBoxLayout, QHBoxLayout,
                               QGridLayout, QLabel, QLineEdit, QTextEdit, QPlainTextEdit, QPushButton, QComboBox,
                               QScrollArea, QSplitter, QFileDialog, QSizePolicy, QToolTip, QMessageBox, QDialog,
                               QTableWidget, QTableWidgetItem, QHeaderView, QAbstractItemView, QSpinBox,
                               QStyledItemDelegate)

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
from teagle_core import appdirs, classify           # classify.DOMAINS_TESTED = the tested-profile panel (scope caveat)
from teagle_core.wsl import curated_coverage_sentence as _curated_coverage   # measured, never retyped on a panel
from teagle_core import domains as domains_mod      # DOMAIN_INFO: the methods panel is derived from it, never retyped
from teagle_core.fetch import (COORD_ASSEMBLIES, all_assemblies,                        # pinned + user-added assemblies
                               complete_gene_model, cross_check_models, retrieve)        # + gene model / transcript fetch
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
    # the genome scan reports isPcr's own coordinates verbatim — a DIFFERENT convention from every other table here
    "Coords (1-based)": "Location of this product in the scanned assembly, exactly as isPcr reports it: "
                        "1-based INCLUSIVE [start, end] — so Len = end − start + 1, not end − start.",
    "Assembly seq": "The assembly sequence this product was found on — the chromosome / contig accession "
                    "as it is named in the downloaded RefSeq assembly.",
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
_FLAG_MARK = {"ok": "", "caution": "!", "warn": "!!"}  # non-colour severity mark, so each ΔG axis is legible without hue


def _amp_call(a, has_locus=True):
    """The ONE call label for an amplicon (table, gel tooltip, FASTA header, caption counts all read it).
    on/off/single are disjoint buckets in the engine (primers.py forces on=False for a self-priming
    product), so single-primer is its own call and never reads as 'off-target'."""
    if a.get("single_primer"):
        return "single-primer"
    if a.get("on_target"):
        return "on-target"
    return "off-target" if has_locus else "priming site"


def _amp_kind(a, has_locus=True):
    """The same call as a punctuation-free FASTA-header token (ontarget / offtarget / singleprimer / primingsite)."""
    return _amp_call(a, has_locus).replace("-", "").replace(" ", "")


def _metric_cell(parts):
    """Fold one or more ΔG metric dicts ({p3, vrna, flag, agree}) into a table cell:
    (text, worst_flag, tooltip, export_value). Shows the worst (most negative) ΔG across engines/primers, appends
    a severity mark (! caution, !! warn) so the per-axis call never rests on colour alone, and marks an engine
    disagreement with ‡; the caller maps worst_flag -> a per-theme colour. The Struct column carries the flag as
    TEXT. export_value is the BARE number (the marks are display cues) so an exported ΔG column stays numeric."""
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
        return ("—", "ok", "no structure predicted", "—")
    worst_flag = max(flags, key=lambda f: _FLAG_ORDER.get(f, 0)) if flags else "ok"
    num = f"{min(dgs):.1f}"
    txt = num + _FLAG_MARK.get(worst_flag, "") + ("‡" if disagree else "")
    return (txt, worst_flag, " · ".join(tips), num)


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
        self.hdr.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Fixed)   # shrink on a narrow window (the meta gloss clips), never force the card wide
        # no setText here — set_collapsed() below owns the header text and escapes '&' for us. Setting the raw
        # title on a QPushButton registered a stray Alt-mnemonic from any '&' in it.
        self.hdr.setCheckable(True)
        self.hdr.clicked.connect(self._toggle)
        self._lay.addWidget(self.hdr)
        self.body = QWidget()
        self.bodylay = QVBoxLayout(self.body)
        self.bodylay.setContentsMargins(theme_mod.sp(12), theme_mod.sp(6), theme_mod.sp(12), theme_mod.sp(12))
        self.bodylay.setSpacing(theme_mod.sp(9))
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
        self.hdr.setText(txt)                             # sentence-case titles (calmer); ALL-CAPS reserved for micro-labels
        self.hdr.setToolTip(f"{self._number}  {self._title}" + (f" · {self._meta}" if self._meta else ""))  # full title when the gloss clips on a narrow window

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


def _kb_links(lab):
    """Make a rich-text label's links Tab-reachable and Enter-activatable. A QLabel's default
    interaction flags are mouse-only, which leaves focusPolicy at NoFocus, so keyboard users could
    never reach the citation/provenance links. Mouse behaviour is unchanged."""
    lab.setTextInteractionFlags(lab.textInteractionFlags()
                                | Qt.LinksAccessibleByMouse | Qt.LinksAccessibleByKeyboard)
    lab.setFocusPolicy(_link_focus(lab.text()))
    return lab


def _link_focus(text):
    # LinksAccessibleByKeyboard raises focusPolicy to StrongFocus unconditionally, so an empty or bare "—" label
    # becomes an invisible dead tab stop. Take focus only while there is actually a link to reach.
    return Qt.StrongFocus if "<a " in (text or "").lower() else Qt.NoFocus


class _MetaLabel(QLabel):
    """Rich-text meta label whose links are Tab-reachable. Re-arms focusPolicy on EVERY setText — evaluating it
    once at construction let later setText calls strand StrongFocus on link-less text (dead tab stop)."""

    def setText(self, t):
        super().setText(t)
        self.setFocusPolicy(_link_focus(t))


def _empty(text):
    l = QLabel(text); l.setObjectName("empty"); l.setWordWrap(True); return l


def _note(text, level="error"):
    """A panel message that is NOT an empty state — a failure (default), a no-result ("warn"), or an
    outcome ("info"/"success"). Same #errbanner / notify vocabulary as the notification dialog, so
    "it failed" never reads as the grey "nothing here yet". In-flight uses BusyBar, not this."""
    l = QLabel(text); l.setObjectName("errbanner"); l.setWordWrap(True)
    if level != "error":                                  # QSS defaults #errbanner to the error look
        l.setProperty("level", level)
    l.style().unpolish(l); l.style().polish(l)
    return l


class _Combo(QComboBox):
    """A combo whose POPUP carries a real keyboard cue. Under the windows11 style the popup's own current-row
    fill measures 1.20:1 (dark) / 1.01:1 (light) against the list — effectively invisible — and the
    stylesheet's selection-background-color never reaches the row. An ::item:selected rule on the popup view
    IS honoured (measured), so build one from the selection colour the theme already resolved into that
    view's palette: no colour literals here, and it re-derives on every theme change."""
    def __init__(self, parent=None):
        super().__init__(parent)
        v = self.view()
        v.setItemDelegate(QStyledItemDelegate(v))          # styles that would use the menu delegate ignore ::item rules

    def showPopup(self):
        v = self.view()
        c = v.palette().highlight().color()                # = the QSS selection-background-color, alpha included
        ring = QColor(c); ring.setAlpha(255)               # same hue at full strength = the theme accent
        # 2px accent ring = the idiom already proven on windows11 for QTableWidget::item:focus. padding-left
        # replaces the native item indent the QSS box model drops, so the row's label does not jump 4px left
        # as the cue lands on it (measured: ink starts at x=11 on every row, current or not).
        v.setStyleSheet(f"QAbstractItemView::item:selected {{ "
                        f"background: rgba({c.red()},{c.green()},{c.blue()},{c.alphaF():.3f}); "
                        f"border: 2px solid {ring.name()}; padding-left: 4px; }}")
        super().showPopup()


def _sl(text):
    """A section label — mono, uppercase, tracked (the web UI's `.lbl`)."""
    l = QLabel(text.upper()); l.setObjectName("sectionlabel"); l.setWordWrap(True); return l


def _export_table_btn(table, base, parent):
    """Visible Excel/CSV/TSV export for a DataTable: the button pops a format menu, then a save dialog
    pre-set to the chosen type. Exports in the table's current on-screen (sorted) order."""
    b = QPushButton("Export table"); b.setProperty("sm", True)
    table.export_base = base                          # the right-click Export… proposes the SAME filename as this button
    def pop():
        fmt = widgets.pick_table_format(b, b.mapToGlobal(b.rect().bottomLeft()))
        if fmt:
            widgets.export_table(table._headers, table.rows_data(), base, parent, fmt=fmt)
    b.clicked.connect(lambda _=False: pop())
    return b


def _qq_note():
    """RepeatMasker runs in -qq (quick) mode. That was an invisible speed choice while only a family NAME
    was shown; once its divergence and coverage are displayed as curation evidence the search sensitivity
    becomes part of the claim, so it is stated beside the numbers rather than buried in the manifest."""
    return _note("Search run in RepeatMasker's quick mode (-qq): faster and less sensitive than a default "
                 "search. Divergence is the raw substitution percentage, not Kimura-corrected, and coverage "
                 "is of the Dfam consensus — both are alignment evidence, not a family-membership test.",
                 "info")


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("TEagle")
        self.setWindowIcon(_app_icon())                   # also set on the top-level window (some Qt/Windows paths key the taskbar HICON off WM_SETICON)
        # open at the design size but never larger than the screen — a 1366x768 laptop must not get an
        # 860px-tall window taller than its display; keep a usable floor so the split layout stays coherent.
        # The floor is scaled by UI_SCALE to match the scaled rail/toolbar minimums (else at 1.5x the results
        # toolbar clips), then clamped to the screen below so a large scale never overhangs a small display.
        w, h = round(1240 * theme_mod.UI_SCALE), round(860 * theme_mod.UI_SCALE)
        minw, minh = round(820 * theme_mod.UI_SCALE), round(560 * theme_mod.UI_SCALE)
        scr = self.screen().availableGeometry() if self.screen() else None
        if scr is not None:
            aw, ah = scr.width() - 40, scr.height() - 60
            w, h = min(w, aw), min(h, ah)
            minw, minh = min(minw, aw), min(minh, ah)
        self.setMinimumSize(minw, minh)
        self.resize(max(w, minw), max(h, minh))
        _saved_theme = str(QSettings("TEagle", "TEagle").value("theme", "light"))   # first run -> light; then remembered
        self.theme = _saved_theme if _saved_theme in ("dark", "light") else "light"
        self.state = {"seq": "", "source": None, "last_rec": None}
        self._loading = False                                 # True while a programmatic load writes the specimen box
        self._pcr_gen = 0                                     # monotonic in-silico-PCR batch id (drops stale sibling results)
        self._design_inflight = False                         # one primer design at a time (self._design_tmpl is shared state)
        self._genome_inflight = False                         # one whole-genome isPcr scan at a time
        self._genome_prep_inflight = False                    # one genome download/prepare at a time (large, one-time)
        self._add_asm_inflight = False                        # one custom-assembly resolve at a time (survives a dialog close/reopen)
        self._annot_inflight = False                          # one whole-genome TE annotation at a time (hours; huge disk/CPU)
        self._annot_result = None                             # last completed landscape (kept for the results window + exports)
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
        self.railToggle = QPushButton("Hide"); self.railToggle.setProperty("sm", True)
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
        cl = QHBoxLayout(chip); cl.setContentsMargins(10, 5, 11, 5); cl.setSpacing(8)
        self.led = QLabel(); self.led.setObjectName("led")   # size via QSS #led min/max-width (scales with UI_SCALE)
        self.statusTxt = QLabel("connecting…"); self.statusTxt.setObjectName("statusTxt")
        cl.addWidget(self.led); cl.addWidget(self.statusTxt)
        h.addWidget(chip)
        # Top-level entry to the backend installer. It used to be reachable only from a secondary button
        # inside panel 03, where a user who had not yet opened that card never saw it — yet it is what
        # turns on Dfam family naming, splice detection and whole-genome scans.
        self.backendBtn = QPushButton("BACKEND"); self.backendBtn.setProperty("sm", True)
        self.backendBtn.setToolTip("Backend installer — install, repair and check the optional Linux (WSL) "
                                   "stack: Dfam family naming, splice detection, whole-genome scans")
        self.backendBtn.clicked.connect(self._open_installer)
        h.addWidget(self.backendBtn)
        sc = QPushButton("SCALE"); sc.setProperty("sm", True)
        sc.setToolTip("Global UI scale — shrink or enlarge the whole interface (applied live; pixel-exact on next launch)")
        sc.clicked.connect(lambda: self._ui_scale_menu(sc))
        h.addWidget(sc)
        tb = QPushButton("THEME"); tb.setProperty("sm", True); tb.clicked.connect(self._toggle_theme)
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
            self.railToggle.setText("Show")
            self.railToggle.setToolTip("Show the specimen panel (Ctrl+B)")
        else:
            self.rail.setVisible(True)
            total = sum(self.split.sizes()) or 1240
            sizes = self._rail_sizes if sum(self._rail_sizes) <= total else [340, max(300, total - 340)]
            self.split.setSizes(sizes)
            self.railToggle.setText("Hide")
            self.railToggle.setToolTip("Hide the specimen panel for more analysis width (Ctrl+B)")

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

    def _set_ui_scale(self, factor):
        """Apply a new UI scale live: re-scale fonts + spacing immediately by re-applying the _scale_px'd QSS.
        No restart, no QT_SCALE_FACTOR (retired). theme_mod.UI_SCALE is seeded before MainWindow so a persisted
        scale is pixel-exact for code-set geometry from the next launch."""
        if abs(factor - theme_mod.UI_SCALE) < 1e-3:
            return
        QSettings("TEagle", "TEagle").setValue("ui_scale", float(factor))
        theme_mod.UI_SCALE = float(factor)
        self._apply_theme()                                   # re-applies scaled QSS + redraws wordmark/mark at the new scale
        z = theme_mod.UI_SCALE                                # re-scale the code-set rail width live so it does not crowd
        self.rail.setMinimumWidth(round(320 * z)); self.rail.setMaximumWidth(round(440 * z))
        inner = self.rail.widget() if hasattr(self.rail, "widget") else None
        if inner is not None:
            inner.setMinimumWidth(round(300 * z)); inner.setMaximumWidth(round(430 * z))
        if getattr(self, "_genome_mgr", None) is not None:    # an open manager was built at the old scale — rebuild it live
            self.engine.submit("genome_list", {}, key="genome_list")
        self._banner(f"UI scale set to {int(factor * 100)}% — applied live (font + spacing); pixel-exact on the next launch.", "info")

    # ---------- rail ----------
    def _build_rail(self):
        rail = QFrame(); rail.setObjectName("rail")
        rail.setMinimumWidth(round(300 * theme_mod.UI_SCALE)); rail.setMaximumWidth(round(430 * theme_mod.UI_SCALE))
        lay = QVBoxLayout(rail); lay.setContentsMargins(12, 12, 12, 12); lay.setSpacing(8)
        lay.addWidget(self._sec("01", "Specimen"))
        accrow = QHBoxLayout()
        self.acc = QLineEdit(); self.acc.setPlaceholderText("accession — e.g. M11240, NC_003075.7")
        accrow.addWidget(self.acc)
        self.fetchBtn = QPushButton("Fetch"); self.fetchBtn.setProperty("sm", True); self.fetchBtn.clicked.connect(self._fetch)
        accrow.addWidget(self.fetchBtn)
        lay.addLayout(accrow)
        self.accMeta = _MetaLabel(""); self.accMeta.setObjectName("cardmeta"); self.accMeta.setWordWrap(True)
        self.accMeta.setTextFormat(Qt.RichText); self.accMeta.setOpenExternalLinks(True); _kb_links(self.accMeta)
        lay.addWidget(self.accMeta)

        # coordinate fetch (collapsed) — organism + chr:start-end, like the UCSC browser position box
        self.coordToggle = QPushButton("▸ Fetch by coordinate"); self.coordToggle.setProperty("link", True)
        self.coordToggle.setAccessibleName("Fetch by coordinate")   # the ▸/▾ glyph is not speakable
        self.coordToggle.setAccessibleDescription("Disclosure toggle — shows or hides the coordinate-fetch fields.")
        self.coordToggle.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Fixed)
        self.coordToggle.clicked.connect(self._toggle_coord); lay.addWidget(self.coordToggle)
        self.coordBox = QWidget(); cb = QVBoxLayout(self.coordBox); cb.setContentsMargins(0, 2, 0, 2); cb.setSpacing(5)
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
        self.coordMeta.setTextFormat(Qt.RichText); self.coordMeta.setOpenExternalLinks(True)
        _kb_links(self.coordMeta); cb.addWidget(self.coordMeta)
        cnote = QLabel("UCSC-style, 1-based (same numbers as the browser). Multi-region: all fetched + recorded; "
                       "analysis runs on the first region.")
        cnote.setObjectName("orient"); cnote.setWordWrap(True); cb.addWidget(cnote)
        self.coordBox.setVisible(False); lay.addWidget(self.coordBox)

        ub = QPushButton("Upload FASTA (.fa / .fasta / .gz)"); ub.setProperty("sm", True)
        ub.clicked.connect(self._upload); lay.addWidget(ub)
        self.seq = QTextEdit(); self.seq.setPlaceholderText("…or paste DNA (FASTA or raw). Real IUPAC validation runs on analyze.")
        self.seq.setMinimumHeight(round(120 * theme_mod.UI_SCALE)); self.seq.textChanged.connect(self._seq_changed)
        lay.addWidget(self.seq)
        row = QHBoxLayout()
        ls = QPushButton("Load example element ▾"); ls.setProperty("link", True); ls.clicked.connect(self._pick_example)
        ls.setToolTip("Five published elements — copia, gypsy, LINE-1, Tc1, Ac — plus a synthetic negative control")
        row.addWidget(ls); row.addStretch(1)
        self.charCount = QLabel("0 nt"); self.charCount.setObjectName("kdim"); row.addWidget(self.charCount)
        lay.addLayout(row)
        self.runBtn = QPushButton("Run analysis"); self.runBtn.setProperty("primary", True)
        self.runBtn.clicked.connect(self._run_analysis); lay.addWidget(self.runBtn)

        # readout gauges (2×2)
        lay.addSpacing(6)
        grid = QGridLayout(); grid.setHorizontalSpacing(6); grid.setVerticalSpacing(6)
        grid.setColumnStretch(0, 1); grid.setColumnStretch(1, 1)
        self.mLen = self._readout(grid, 0, "Length"); self.mGC = self._readout(grid, 1, "GC")
        self.mN = self._readout(grid, 2, "N content"); self.mValid = self._readout(grid, 3, "IUPAC")
        lay.addLayout(grid)
        lay.addSpacing(4)
        self.rRecords = _MetaLabel("—"); self.rStruct = _MetaLabel("—"); self.rOrf = _MetaLabel("—")
        for val in (self.rRecords, self.rStruct, self.rOrf):
            val.setObjectName("cardmeta"); val.setTextFormat(Qt.RichText)
        self.rRecords.setOpenExternalLinks(True)              # RECORDS -> the real fetched source accession (external DB link)
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

        lay.addSpacing(6)
        self.envHdr = QPushButton("▸ Environment"); self.envHdr.setObjectName("cardhdr")
        self.envHdr.setAccessibleName("Environment")
        self.envHdr.setAccessibleDescription("Disclosure toggle — shows or hides the backend environment status.")
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
        wrap.setMinimumWidth(round(320 * theme_mod.UI_SCALE)); wrap.setMaximumWidth(round(440 * theme_mod.UI_SCALE))
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
        wrap.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)  # page body never scrolls sideways; wide tables scroll
        self.resultsScroll = wrap                          # inside their own viewport (per the crowding invariant)
        inner = QWidget(); inner.setObjectName("central")
        self.results = QVBoxLayout(inner); self.results.setContentsMargins(4, 0, 8, 0); self.results.setSpacing(9)

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
        self.wslSource = _Combo()
        self.wslSource.addItems(["Loaded specimen (panel 01)", "Paste sequence…"])
        srcrow.addWidget(self.wslSource); srcrow.addStretch(1)
        card.bodylay.addLayout(srcrow)
        self.wslPaste = QTextEdit(); self.wslPaste.setPlaceholderText("Paste a FASTA or raw DNA sequence to annotate against Dfam")
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

    # ---------- Splice detection ----------
    def _build_splice_card(self):
        card = CollapsibleCard("SP", "Splice detection (de novo)", "exon–intron inferred by aligning a supplied transcript")
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
        self._splice_tx_origin = "external"                   # a pasted/right-clicked transcript is treated as external
        self.spliceTx.textChanged.connect(self._on_splice_tx_changed)
        card.bodylay.addWidget(self.spliceTx)
        self.spliceBtn = QPushButton("Detect exons / introns"); self.spliceBtn.setProperty("primary", True)
        self.spliceBtn.setEnabled(False); self.spliceBtn.clicked.connect(self._splice)
        srow = QHBoxLayout(); srow.addWidget(self.spliceBtn); srow.addStretch(1)   # the panel's run action, at its own width
        card.bodylay.addLayout(srow)
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
        self.pPreset = _Combo()
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
        self.designBtn.setEnabled(False); self.designBtn.clicked.connect(self._design)
        row.addWidget(self.designBtn)
        rb = QPushButton("Reset"); rb.setProperty("sm", True); rb.clicked.connect(lambda: self._apply_preset("standard"))
        rb.setAccessibleName("Reset primer parameters to the standard preset")
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
        self.pcrClear = QPushButton("Clear"); self.pcrClear.setProperty("sm", True)
        self.pcrClear.setAccessibleName("Clear the in-silico PCR primer-pair queue")
        self.pcrClear.setEnabled(False); self.pcrClear.clicked.connect(self._pcr_clear)
        srow.addWidget(self.pcrStageAll); srow.addWidget(self.pcrClear); srow.addStretch(1)
        card.bodylay.addLayout(srow)
        grid = QGridLayout(); grid.setHorizontalSpacing(8)
        self._grid_fields(grid, [("pcrMM", "Max mismatches", "2"), ("pcrTP", "3′ exact bases", "5"),
                                 ("pcrPmin", "Prod min (bp)", "70"), ("pcrPmax", "Prod max (bp)", "1000")], cols=4)
        card.bodylay.addLayout(grid)
        card.bodylay.addWidget(_sl("Optional custom background (FASTA) — off-target search"))
        self.pcrBg = QTextEdit(); self.pcrBg.setMaximumHeight(round(60 * theme_mod.UI_SCALE))
        self.pcrBg.setPlaceholderText("Optional: paste extra background sequence(s) to reveal off-target amplicons.")
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

    def _build_genome_card(self):
        # dedicated whole-genome off-target scan — LOCAL isPcr against a downloaded RefSeq assembly. The organism
        # dropdown lists ONLY downloaded+verified genomes (from genome_list). On-target = the product at the
        # specimen's own genome locus when it sits in the scanned assembly; the rest are off-target paralogs.
        card = CollapsibleCard("06", "Whole-genome off-target scan", "a designed pair vs a downloaded RefSeq genome (isPcr)")
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
                       "genome. Products are listed on-target first and sealed with the assembly version + checksum.")
        gnote.setObjectName("orient"); gnote.setWordWrap(True); card.bodylay.addWidget(gnote)
        self.genomeBody = QVBoxLayout(); gb = QWidget(); gb.setLayout(self.genomeBody)
        self.genomeBody.addWidget(_empty("Design a primer pair (panel 04) and download a genome, then scan here."))
        card.bodylay.addWidget(gb)
        return card

    def _grid_fields(self, grid, specs, cols):
        for idx, (fid, label, default) in enumerate(specs):
            r, c = divmod(idx, cols)
            cell = QWidget(); cl = QVBoxLayout(cell); cl.setContentsMargins(0, 0, 0, 0); cl.setSpacing(1)
            lab = QLabel(label.upper()); lab.setObjectName("kdim")
            ed = QLineEdit(default); ed.setMaximumWidth(round(90 * theme_mod.UI_SCALE))
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

    def _uppercase_buttons(self, root=None):
        """Uppercase action buttons under `root` (default: the whole window, dialogs parented to it included).
        A freshly rendered panel/dialog passes itself so it never sits in sentence case until the next refresh."""
        widgets.uppercase_buttons(root if root is not None else self)

    def _repolish(self, w):
        w.style().unpolish(w); w.style().polish(w)

    def _toggle_theme(self):
        self.theme = "light" if self.theme == "dark" else "dark"
        QSettings("TEagle", "TEagle").setValue("theme", self.theme)   # remember the choice for next launch
        self._apply_theme()
        d = self.state.get("lastPrimers")                     # QC ΔG cell colours are per-theme -> re-render on toggle
        if isinstance(d, dict) and d.get("candidates"):
            self._render_primers(d)

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

    # ---------- input handling ----------
    def _set_seq(self, text):
        """Programmatic specimen load (fetch / upload / sample). Wrapped in the _loading guard so the
        textChanged handler keeps state['seq'] in sync but does NOT treat the load as a user edit that
        would wipe the source the loader is about to set."""
        self.state["features"] = None                         # drop a prior accession's gene model (fetch re-sets it after)
        self._splice_tx_origin = "external"                   # a new specimen invalidates any prior record-transcript flag
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
        menu.exec(QCursor.pos())

    def _load_example(self, accession):
        from teagle_core import examples
        fasta = examples.load(accession)
        if not fasta:
            return self._banner("This build does not include that example sequence.", "warn")
        self._apply_example(fasta)

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
        self._update_splice_ref()                             # _set_seq runs under the _loading guard, so refresh the label here

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
            self.runBtn.setEnabled(True); self.runBtn.setText("Run analysis")
        elif key == "primers":
            self._design_inflight = False
            self._pending_domain = None                       # a failed routed design must not leave a stale re-anchor or busy body
            self.designBtn.setEnabled(True); self.designBtn.setText("Design primers")
            self._set_body(self.primBody, _note("Primer design failed — adjust the parameters or region, then try again."))
        elif key == "annotate":
            self.annotateBtn.setEnabled(True); self.annotateBtn.setText("Run family annotation")
            # _reset_buttons is error-only (success re-enables in _on_annotate), so clear the stuck busy body here —
            # a BadRequest/fault never reaches _on_annotate, which would otherwise leave "Running RepeatMasker…" spinning
            self._set_body(self.wslBody, _empty("Run RepeatMasker against Dfam to name the TE family. Family naming is the Linux (WSL) backend."))
        elif key == "splice":
            self.spliceBtn.setEnabled(True); self.spliceBtn.setText("Detect exons / introns")
            self._set_body(self.spliceBody, _empty("Align a transcript to the loaded sequence to resolve exon–intron structure de novo."))
        elif key == "fetch":
            self._set_fetch_enabled(True)                     # fetch failed — re-enable so the user can retry
            for lbl in (self.accMeta, self.coordMeta):        # clear only the in-flight indicator, keep a prior result
                if lbl.text() == "fetching…":
                    lbl.setText("")
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
        if key.startswith("pcr#"):                            # a failed pair fills its slot so the batch still renders
            self._pcr_slot(key, {"error": msg, "amplicons": []})
        elif key == "genome_pcr":
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

    def _on_failed(self, key, msg, trace):
        if key == "genome_prepare_log":
            return                                            # a background poll blip: never banner, never spam a traceback
        sys.stderr.write(trace + "\n")                        # log the trace for diagnosis (all real faults)
        if key == "health":                                   # state the backend is down IN WORDS — the chip must not
            self.statusTxt.setText("backend unavailable")     # rely on the LED hue alone, nor keep saying "connecting…"
            self.led.setProperty("live", False); self._repolish(self.led)
        if key.startswith("pcr#"):
            self._pcr_slot(key, {"error": msg, "amplicons": []})
        elif key == "genome_pcr":                             # a genome scan fault is a handled outcome, not an "unexpected error"
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
        self._banner(f"unexpected error: {msg}")

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

    def _toggle_coord(self):
        vis = not self.coordBox.isVisible()
        self.coordBox.setVisible(vis)
        self.coordToggle.setText("▾ Fetch by coordinate" if vis else "▸ Fetch by coordinate")

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
        _kb_links(self.coordMeta)                             # text now carries a link -> re-arm the tab stop
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
            _kb_links(self.accMeta)                           # text now carries a link -> re-arm the tab stop
            self.coordMeta.setText("")                    # clear the other specimen identity
            self.state["source"] = {k: res.get(k) for k in ("accession", "organism", "title", "length", "moltype") if res.get(k) is not None}
            self.state["source"]["sourceUrl"] = src_url        # traceable RECORDS link (display/provenance label, not sealed)
            self.state["features"] = res.get("features")
        # auto-fill species for WSL family annotation if present
        if hasattr(self, "wslSpecies") and org:
            self._set_species(org)
        self._update_splice_ref()

    def _on_analyze(self, res):
        self.runBtn.setEnabled(True); self.runBtn.setText("Run analysis")
        self._clear_banner()
        if res.get("warning"):
            self._banner(res["warning"], level="warn")
        recs = res.get("records", [])
        if not recs:
            return
        self.state["analysis"] = res
        self.state["records"] = recs
        self.state["analyzed_clean"] = self._clean_seq(self.state.get("analyzed_seq", ""))   # snapshot, not the live box — a keystroke mid-analysis must not defeat the stale-block guard
        self._show_record(0)

    # engine.analyze already classifies EVERY record; the UI used to take records[0] and discard the rest,
    # so a user arriving with a shortlist saw one answer and no indication the others had been computed.
    _MAX_RECORDS = 200

    def _show_record(self, index: int):
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
        self._update_splice_ref()
        self.designBtn.setEnabled(True); self.designHint.setText("")
        comp = rec.get("composition", {})
        self.mLen.setText(f"{comp.get('length', 0):,}")
        self.mGC.setText(f"{comp.get('gc', 0)}%")
        self.mN.setText(f"{comp.get('n', 0)}%")
        self.mValid.setText("valid" if rec.get("valid") else "invalid")
        self.mValid.setProperty("state", "good" if rec.get("valid") else "bad"); self._repolish(self.mValid)
        self._set_trace_counts(recs, rec)
        self._render_struct_card(rec, res)
        self._uppercase_buttons()

    def _record_summary_rows(self, recs):
        """One row per record: what was actually computed for each, not a promise that it was."""
        rows = []
        for i, r in enumerate(recs[:self._MAX_RECORDS]):
            cl = r.get("classification") or {}
            comp = r.get("completeness") or cl.get("completeness") or {}
            doms = sorted({d["domain"] for d in (r.get("domains") or [])})
            struct = sorted({e["type"].split(" (")[0] for e in (r.get("structural") or [])})
            rows.append([i + 1, r.get("id", ""), f"{(r.get('composition') or {}).get('length', 0):,}",
                         f"{(r.get('composition') or {}).get('gc', 0)}%",
                         "valid" if r.get("valid") else "INVALID",
                         cl.get("te_class") or "—", cl.get("confidence") or "—",
                         "·".join(doms) or "—", ", ".join(struct) or "—",
                         comp.get("tier") or "—"])
        return rows

    _RECORD_HEADERS = ["#", "Record", "Length", "GC", "Input", "Class", "Confidence",
                       "Domains", "Structural", "Completeness"]

    def _accent_link(self, href, text):
        return f"<a href='{href}' style='color:{theme_mod.ACCENT[self.theme]};text-decoration:none'>{text}</a>"

    def _scroll_to(self, widget):
        """Expand a result card and bring it into view (in-app navigation for the traceable rail links)."""
        try:
            widget.expand(); self.resultsScroll.ensureWidgetVisible(widget)
        except Exception:
            pass

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
            extra = f"  +{len(recs) - 1} more" if len(recs) > 1 else ""
            self.rRecords.setText(self._accent_link(src["sourceUrl"], str(acc)) + extra)
            _kb_links(self.rRecords)                          # text now carries a link -> re-arm the tab stop
        elif acc:
            self.rRecords.setText(str(acc))
        else:
            self.rRecords.setText(f"user-supplied ({len(recs)})" if recs else "—")
        self.rStruct.setText(self._accent_link("#struct", str(len(rec.get("structural", [])))))
        self.rStruct.setToolTip(" · ".join(e["type"].split(" (")[0] for e in rec.get("structural", [])) or "none detected")
        self.rOrf.setText(self._accent_link("#orfs", str(len(rec.get("orfs", [])))))

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
        guides = [{"start": e[q][0], "end": e[q][1], "color": theme_mod.OK.get("LTR", "#0072B2")}
                  for e in struct for q in ("five_prime", "three_prime") if e.get(q)]
        rid = (rec.get("id") or "locus").split()[0]
        # ONE reused window (hidden, not destroyed, on close) whose body is re-populated per record — never a
        # rebuilt widget, which is what collapsed a reused dialog in an earlier release.
        dlg = getattr(self, "_selfsim_dlg", None)
        if dlg is None:
            dlg = QDialog(self)
            dlg.setObjectName("selfsimdlg")
            dlg.resize(940, 820)
            _l = QVBoxLayout(dlg)
            _l.setContentsMargins(12, 12, 12, 12)
            dlg._panel = widgets.DotPanel(base_name="TEagle_selfsim", parent=dlg)
            _l.addWidget(dlg._panel)
            self._selfsim_dlg = dlg
            self._uppercase_buttons()
        panel = dlg._panel
        dlg.setWindowTitle(f"Self-similarity — {rid}")
        panel.base_name = f"TEagle_{rid}_selfsim"
        panel.apply_app_theme(self.theme)
        panel.zoom = None                                 # fit the new record to the window
        panel.set_matrix(m, guides=guides, threshold=dotplot.above_chance(m), scope=dotplot.scope_note(m))
        dlg.show(); dlg.raise_(); dlg.activateWindow()
        panel._fit()
        self._uppercase_buttons()

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
        chosen = menu.exec(QCursor.pos())
        if chosen is None:
            return
        fmt = "gff3" if chosen is a_gff else "bed"
        seqid = (rec.get("id") or "locus").split()[0]
        base = f"TEagle_{seqid}"
        path, _ = QFileDialog.getSaveFileName(self, f"Export {fmt.upper()}", base + "." + fmt,
                                              f"{fmt.upper()} (*.{fmt})")
        if not path:
            return
        if not path.lower().endswith("." + fmt):
            path += "." + fmt
        # the record's OWN sequence (not analyzed_clean, which concatenates every record of a multi-FASTA
        # paste — embedding that after ##FASTA would contradict this record's own coordinates)
        seq = rec.get("seq") or None
        text = (gff3.to_gff3(rec, seqid=seqid, sequence=seq,
                             source_note=f"TEagle {teagle_version()}")
                if fmt == "gff3" else gff3.to_bed(rec, seqid=seqid))
        with open(path, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(text)
        self._banner(f"Annotation written to {os.path.basename(path)}.", level="success")

    def _render_record_table(self, card):
        """Per-record summary for a multi-record input, above the selected record's detail."""
        recs = self.state.get("records") or []
        if len(recs) < 2:
            return
        rows = self._record_summary_rows(recs)
        t = DataTable(self._RECORD_HEADERS, GLOSS)
        t.set_rows(rows)
        t.set_row_menu(lambda r: [("Show this record below", lambda i=r: self._show_record(i))])
        sel = self.state.get("record_index", 0)
        head = QLabel(f"<b>{len(recs)} records analysed</b> — every one was classified; the detail below is "
                      f"record {sel + 1}. Right-click a row to show it, or export the table."
                      + (f" Showing the first {self._MAX_RECORDS} of {len(recs)}."
                         if len(recs) > self._MAX_RECORDS else ""))
        head.setTextFormat(Qt.RichText); head.setWordWrap(True); head.setObjectName("cardmeta")
        wrap = QWidget(); wl = QVBoxLayout(wrap)
        wl.setContentsMargins(0, 0, 0, 0); wl.setSpacing(theme_mod.sp(6))
        wl.addWidget(head); wl.addWidget(t)
        frow = QHBoxLayout(); frow.addStretch(1)
        frow.addWidget(_export_table_btn(t, "TEagle_records", self))
        wl.addLayout(frow)
        card.bodylay.addWidget(wrap)

    def _render_struct_card(self, rec, res):
        card = self.card_struct
        card.clear_body()
        card.expand()
        self._render_record_table(card)
        cl = rec.get("classification") or {}
        banner = QFrame(); banner.setObjectName("classbn")
        bl = QVBoxLayout(banner)
        bl.setContentsMargins(theme_mod.sp(16), theme_mod.sp(14), theme_mod.sp(16), theme_mod.sp(14)); bl.setSpacing(theme_mod.sp(9))
        # TIER 1 — CALL: te_class + confidence chip + superfamily
        top = QHBoxLayout(); top.setSpacing(theme_mod.sp(10))
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
        kls.setOpenExternalLinks(True); _kb_links(kls); bl.addWidget(kls)
        # TIER 2 — WHY: the reasoning, with reading line-height. Drop the leading "Classified as … (confidence)."
        # sentence — it just repeats the CALL above (display-only trim; classify.py.explanation stays byte-stable).
        if cl.get("explanation"):
            why = re.sub(r"^\s*Classified as [^.]+\.\s*", "", cl["explanation"])
            ex = QLabel(f"<div style='line-height:150%'>{why}</div>"); ex.setObjectName("classexp")
            ex.setWordWrap(True); ex.setTextFormat(Qt.RichText)
            bl.addWidget(ex)
        # TIER 3 — SCOPE: reliability & completeness. An honest limitation → stays default-VISIBLE, grouped + labelled.
        comp = cl.get("completeness")                         # scoped structural-completeness (Axis 2 of reliability)
        if comp:
            rl = QLabel("RELIABILITY & COMPLETENESS"); rl.setObjectName("sectionlabel"); bl.addWidget(rl)
            arch = cl.get("order") or " – ".join(comp.get("present", []))
            miss = comp.get("missing") or []
            # A record flagged 3′-truncated must not also read "intact": the card was asserting completeness
            # directly above its own truncation warning. Presentation only — classify.py's tier value is
            # untouched (fixing the tier logic would move benchmark rows and is the user's call).
            tier = comp["tier"]
            if any("truncat" in str(n) for n in (rec.get("notes") or [])) and "intact" in tier.lower():
                tier += " at the domain level — but the sequence is 3′-truncated (see note below)"
            line = (f"<b>Structural completeness:</b> {tier}  ·  {comp.get('kind','')}"
                    + (f"<br><b>Domain architecture:</b> {arch}" if arch else "")
                    + (f"  ·  not detected: {', '.join(miss)}" if miss else ""))
            cw = QLabel(f"<div style='line-height:150%'>{line}</div>"); cw.setObjectName("classexp"); cw.setTextFormat(Qt.RichText); cw.setWordWrap(True)
            bl.addWidget(cw)
            # the essential caveat stays visible; the full methodology collapses behind a disclosure (less crowding)
            caveat = QLabel("Structural evidence in this one sequence, scored against the tested Pfam panel — "
                            "not a claim about expression, transposition competence, or any individual genome.")
            caveat.setObjectName("cardmeta"); caveat.setWordWrap(True); bl.addWidget(caveat)
            scope = QLabel(f"Domains tested: {comp.get('scope','')}. “Not detected” is relative to this profile "
                           "panel — a divergent or unmodelled domain reads as not-detected, not as element decay "
                           "(completeness after Wicker 2007 / TEsorter / LTR_retriever). The tier reports how much of "
                           "the expected architecture is present at the domain level in THIS SEQUENCE. It says "
                           "nothing about three separate questions it is easily confused with: whether the element "
                           "is transcribed, whether it retains transposition or infection competence, and whether "
                           "any given individual or population carries this insertion. Answering those needs "
                           "expression data, a functional assay, and population genotyping respectively "
                           "(Lanciano &amp; Cristofari 2020, Nat Rev Genet 21:721-736).")
            scope.setObjectName("cardmeta"); scope.setWordWrap(True)
            self._disclosure(bl, "Scope and methods", scope, expanded=False)
        else:
            # a no-evidence / structural-only call still needs its SCOPE stated — the essential caveat stays visible
            rl = QLabel("RELIABILITY & COMPLETENESS"); rl.setObjectName("sectionlabel"); bl.addWidget(rl)
            caveat = QLabel("No panel domain was detected — scoped to the tested Pfam panel, not proof the sequence is not a TE.")
            caveat.setObjectName("cardmeta"); caveat.setWordWrap(True); bl.addWidget(caveat)
            scope = QLabel(f"Domains tested: {classify.DOMAINS_TESTED}. None were detected — this is relative to this "
                           "bundled Pfam profile panel; a divergent or unmodelled domain reads as not-detected, not as "
                           "proof the sequence is not a transposable element.")
            scope.setObjectName("cardmeta"); scope.setWordWrap(True)
            self._disclosure(bl, "Scope and methods", scope, expanded=False)
        card.bodylay.addWidget(banner)

        # genome viewer
        model = figures.gv_tracks_from_rec(rec)
        if model["tracks"]:
            gv = GenomePanel(svg_genome, "TEagle_genome")
            gv.apply_app_theme(self.theme)                # open in the current app theme
            gv.set_model(model)
            gv.set_feature_menu(self._region_menu)
            gv.setMinimumHeight(round(260 * theme_mod.UI_SCALE))
            card.bodylay.addWidget(gv)

        # retroviral transcript architecture (ERV) — the correct coding-organisation model + cis-element legend
        arch = rec.get("retroviral")
        if arch:
            # swatches derive from the figure palette (theme.ARCHCOL / theme.CISCOL) instead of repeating the
            # hex here — a legend that duplicated the hues could drift out of step with the bands it labels.
            _cis = [(k, lbl) for k, pre, lbl in (("PBS", "PBS", "PBS"), ("PPT", "PPT", "PPT"),
                                                 ("PAS", "polyA-signal", "PAS · polyA signal (motif)"))
                    if any(e["type"].startswith(pre) for e in rec.get("structural", []))]
            leg = (f"<span style='color:{theme_mod.ARCHCOL['exon']}'>■</span> env exon &nbsp; "
                   f"<span style='color:{theme_mod.ARCHCOL['intron']}'>■</span> gag–pro–pol intron (fused polyprotein)"
                   + "".join(f" &nbsp; <span style='color:{theme_mod.CISCOL[k]}'>■</span> {lbl}" for k, lbl in _cis))
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
            t.set_rows([self._struct_row(e) for e in struct], tips=[self._struct_tips(e) for e in struct])
            t.set_row_menu(lambda r: self._struct_menu(struct[r]))
            t.setMaximumHeight(round(180 * theme_mod.UI_SCALE))
            card.bodylay.addWidget(t)
            # The polyA-signal row is a MOTIF, not a located cleavage site. Its caveat rides in the row's
            # tooltip, but a tooltip is not a disclosure a user is guaranteed to see — so when the row is
            # present the limit is stated in the card itself, next to the table that shows it.
            if any(e.get("type", "").startswith("polyA-signal") for e in struct):
                pas_note = QLabel("The polyA-signal row is a sequence motif with its downstream element — advisory "
                                  "context only. It does not locate the U3–R–U5 boundaries, the cleavage site, or "
                                  "the transcript end, which need RNA evidence TEagle does not use.")
                pas_note.setObjectName("orient"); pas_note.setWordWrap(True)
                card.bodylay.addWidget(pas_note)
            srow = QHBoxLayout(); srow.addStretch(1)          # exportable like every sibling results table
            srow.addWidget(_export_table_btn(t, "TEagle_structural", self))
            card.bodylay.addLayout(srow)

        # ORFs
        orfs = rec.get("orfs", [])
        if orfs:
            card.bodylay.addWidget(_sl(f"ORFs (≥40 aa) — {len(orfs)}"))
            t = DataTable(ORF_COLS, GLOSS)
            t.set_rows([[o["strand"], o["frame"], o["start"], o["end"], o["length_aa"]] for o in orfs])
            t.set_row_menu(lambda r: self._feat_menu(orfs[r]["start"], orfs[r]["end"], orfs[r]["strand"],
                                                     f"ORF_{orfs[r]['strand']}{orfs[r]['frame']}",
                                                     src_seq=rec.get("seq")))
            t.setMaximumHeight(round(160 * theme_mod.UI_SCALE))
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
                                                     doms[r]["domain"], protein=doms[r].get("protein"),
                                                     src_seq=rec.get("seq")))
            t.setMaximumHeight(round(180 * theme_mod.UI_SCALE))
            card.bodylay.addWidget(t)
            dhint = QLabel("The last column is <b>Conf</b> (per-domain confidence). On a narrow window, scroll the table "
                           "sideways — or collapse the specimen panel (Ctrl+B) — to reach every column.")
            dhint.setObjectName("orient"); dhint.setTextFormat(Qt.RichText); dhint.setWordWrap(True)
            card.bodylay.addWidget(dhint)
            drow = QHBoxLayout(); drow.addStretch(1); drow.addWidget(_export_table_btn(t, "TEagle_domains", self))
            card.bodylay.addLayout(drow)

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
        card.bodylay.addLayout(srow)

        # Annotation export: the whole record as GFF3/BED so the call survives into a genome browser.
        # Placed on the structure card because that is where the annotation being exported is shown.
        arow = QHBoxLayout(); arow.addStretch(1)
        ann = QPushButton("Export annotation ▾"); ann.setProperty("sm", True)
        ann.setToolTip("GFF3 (Sequence Ontology terms, sub-features, embedded FASTA) or BED6 intervals")
        ann.clicked.connect(lambda: self._export_annotation(rec))
        arow.addWidget(ann)
        card.bodylay.addLayout(arow)

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
                gvg = GenomePanel(svg_genome, "TEagle_genemodel"); gvg.apply_app_theme(self.theme); gvg.set_model(gmodel)
                gvg.set_feature_menu(self._region_menu); gvg.setMinimumHeight(round(200 * theme_mod.UI_SCALE))
                gmlay.addWidget(gvg)
            ptr = QPushButton("For exon–intron from a supplied transcript, use Splice detection →")   # ANNOTATED vs MEASURED
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
                card.bodylay.addWidget(gmtog)
            card.bodylay.addWidget(gmbox)

        for note in rec.get("notes", []):
            n = QLabel("• " + note); n.setObjectName("orient"); n.setWordWrap(True)
            card.bodylay.addWidget(n)

        # explicit methodology — which database / consensus / parameters define each evidence layer
        mtoggle = QPushButton("▸ Methods && databases"); mtoggle.setProperty("link", True)   # '&&' -> a literal '&' (QPushButton eats a lone '&' as a mnemonic)
        mtoggle.setAccessibleName("Methods & databases")
        mtoggle.setAccessibleDescription("Disclosure toggle — shows or hides the methodology and database provenance.")
        mtoggle.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Fixed)
        mbox = QLabel(self._methods_html()); mbox.setObjectName("orient"); mbox.setWordWrap(True)
        mbox.setTextFormat(Qt.RichText); mbox.setOpenExternalLinks(True); _kb_links(mbox); mbox.setVisible(False)
        def _tgl():
            v = not mbox.isVisible(); mbox.setVisible(v)
            mtoggle.setText("▾ Methods && databases" if v else "▸ Methods && databases")
        mtoggle.clicked.connect(_tgl)
        card.bodylay.addWidget(mtoggle); card.bodylay.addWidget(mbox)

    def _methods_html(self):
        """A plain statement of exactly what defines each evidence layer — the database/model, the consensus
        source, and the thresholds — so the annotation is never a black box."""
        p, h, w, dfam = (self._src_html(k) for k in ("Pfam", "HMMER", "Wicker2007", "Dfam"))
        return (
            "<b>Protein domains</b> — profile-HMM search (HMMER" + h + ", run in-process via pyhmmer) of the "
            "6-frame ORFs (≥ 40 aa) against a bundled Pfam-A" + p + " TE-domain profile set "
            f"({self._domain_panel_html()}). A hit is kept at per-domain "
            "E-value ≤ 1e-3; the gag + env models let TEagle recover the full GAG–POL–ENV architecture of ERVs (HERV-K, "
            "-W, -L …), not just the pol enzymes.<br>"
            "<b>Structural evidence</b> — heuristic terminal-repeat detectors (no external database): LTR by k-mer "
            "seed + diagonal cluster (k=13, ≥ 80 bp, ≥ 80% identity, ≥ 4 anchors); TIR by a terminal inverted-repeat "
            "scan plus a k-mer-vs-reverse-complement search; poly-A/poly-T tail ≥ 8 bp; TSD as a 4–12 bp exact "
            "flanking direct repeat. For an LTR element the two retroviral cis-elements are also scanned: <b>PBS</b>, "
            "the best reverse-complement match to a bundled 18 nt primer-tRNA panel within the 44 nt leader after the "
            "5′ LTR, reported at ≥ 55% identity but only NAMED at ≥ 72% — a weaker match is reported as priming tRNA "
            "undetermined and marked <i>tentative</i>, because an endogenised PBS is diverged and a 55% match on random "
            "sequence is not rare; and <b>PPT</b>, a purine-dense run abutting the 3′ LTR (≥ 9 bp, ≥ 82% purine, ≤ 2 "
            "pyrimidine defects, 30 nt window). Coordinates are 0-based half-open, and each row lists its own detection method.<br>"
            "<b>Superfamily / class</b> — the Wicker&nbsp;2007" + w + " scheme, derived from the domain architecture and "
            "structural context; Copia vs Gypsy is called from the strand-aware integrase-vs-RT translation order, not "
            "ORF length; an env domain with paired LTRs flags an endogenous retrovirus (ERV).<br>"
            "<b>Reliability</b> — reported on two independent, citable axes rather than one fabricated number: (1) a "
            "per-domain call confidence from the HMMER i-Evalue (Eddy 2011); (2) a categorical structural-completeness "
            "tier — intact / near-complete / partial / structural-only — mapped to the autonomous/intact criteria of "
            "Wicker 2007, TEsorter and LTR_retriever, and always scoped to the domain models actually tested.<br>"
            "<b>Family naming</b> (optional, WSL backend) — RepeatMasker (RMBLAST) against the curated Dfam&nbsp;4.0" + dfam +
            " library; this is the only step that makes a database family call, and it is absent from the offline path.")

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
        "PBS": "reverse-complement match to a bundled 18 nt primer-tRNA panel in the 44 nt leader (≥ 55% identity)",
        "PPT": "purine-run extension from the 3′-LTR boundary (≥ 9 bp, ≥ 82% purine, ≤ 2 defects, 30 nt window)",
        "TSD": "exact 4–12 bp direct repeat flanking the element",
        "poly-A": "terminal homopolymer run (≥ 8 bp)",
        "poly-T": "terminal homopolymer run (≥ 8 bp)",
        "polyA-signal": ("poly(A)-signal hexamer panel in the 3′ LTR, gated on a GU/U-rich downstream "
                         "element 20–60 nt past the hexamer (≥ 65% G+T, ≥ 4 T)"),
    }

    def _struct_row(self, e):
        fp, tp = e.get("five_prime"), e.get("three_prime")
        sp = None
        if fp and tp:                                     # a terminal-repeat PAIR (LTR/TIR): show BOTH copies,
            coords = f"{fp[0]}–{fp[1]}  ·  {tp[0]}–{tp[1]}"   # matching the two blocks drawn in the genome viewer
        else:
            sp = e.get("pos") or e.get("upstream") or e.get("element_span") or [None, None]
            coords = f"{sp[0]}–{sp[1]}" if sp[0] is not None else ""
        arm = e.get("ltr_len") or e.get("tir_len") or e.get("length") or ""
        if arm == "" and sp and sp[0] is not None:        # PBS carries no length key — take Len from its own span
            arm = sp[1] - sp[0]
        metric = (f"{e['identity']}%" if e.get("identity") is not None else e.get("motif", ""))
        if e.get("note") and not e.get("confident", True):   # a hedged call (non-confident PBS) must not read as a
            metric += " · tentative"                          # hard number here; the full hedge is the cell tooltip
        return [e["type"], coords, arm, metric,
                e.get("method") or self._STRUCT_METHOD.get(e["type"].split(" ")[0], "")]

    def _struct_tips(self, e):
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
        return [note, None, None, note, None]


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
        with open(path, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(text if text.endswith("\n") else text + "\n")
        self._banner(f"Written to {os.path.basename(path)}.", "success")

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
        items = [(f"Copy FASTA", lambda: self._copy(f">{fid}_{start}-{end}{rev}\n{dna}")),
                 (f"Copy DNA", lambda: self._copy(dna)),
                 (f"Copy coords ({start}–{end} {strand})", lambda: self._copy(f"{start}-{end} {strand}"))]
        if protein:
            items.append(("Copy protein", lambda: self._copy(protein)))
        # action items are contextual: a structural motif (LTR / TIR / TSD / PBS / PPT / poly-A) gets Copy actions only —
        # a primer designed INSIDE a repeat would be non-specific, and "send to splice" on a motif is meaningless. A
        # coding/transcript-like feature that is long enough (room for two primers + a product) gets primer + splice + sub-region.
        long_enough = len(dna) >= 50
        structural_motif = str(kind or "").split(" ")[0] in ("LTR", "TIR", "TSD", "PBS", "PPT", "poly-A", "poly-T",
                                                             "polyA-signal")
        if long_enough and not structural_motif:
            items.append(("Design primer here", _design))
            items.append(("Send to splice detection",
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
        seq = src_seq if src_seq is not None else self.state.get("seq", "")
        if not seq:
            self._banner("No sequence is loaded to take flanks from.", "info"); return
        n = len(seq)
        fid = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(label)).strip("_") or "feature"
        dlg = QDialog(self); dlg.setWindowTitle("Flanking sequence"); dlg.setModal(False)
        dlg.setAttribute(Qt.WA_DeleteOnClose, True)
        lay = QVBoxLayout(dlg)
        lay.setContentsMargins(theme_mod.sp(14), theme_mod.sp(12), theme_mod.sp(14), theme_mod.sp(12))
        lay.setSpacing(theme_mod.sp(8))
        head = QLabel(f"Flanks of {label} ({start}–{end}, {strand} strand)")
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
                        "concatenation — joining them would create a junction that does not exist.")
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

        def refresh():
            p = parts()
            if not p:
                info.setText("No flanking sequence is available on that side — the feature reaches the end "
                             "of the record."); return
            info.setText(" · ".join(f"{nm}: {len(s)} bp" for nm, s in p)
                         + ("" if all(len(s) == spin.value() for _, s in p)
                            else "  (clipped by the end of the record)"))
        side.currentIndexChanged.connect(lambda *_: refresh()); spin.valueChanged.connect(lambda *_: refresh())
        refresh()

        def fasta():
            return "\n".join(f">{nm}\n{s}" for nm, s in parts())

        def design():
            p = parts()
            if not p:
                return
            nm, s = p[0]
            self._design_for_domain(0, len(s), nm, seq=s)
            dlg.accept()

        brow = QHBoxLayout(); brow.addStretch(1)
        for txt, fn in (("Copy FASTA", lambda: self._copy(fasta())),
                        ("Export FASTA…", lambda: self._save_text(fasta(), fid + "_flank", "fasta")),
                        ("Design primers", design)):
            b = QPushButton(txt); b.setProperty("sm", True); b.clicked.connect(fn); brow.addWidget(b)
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

    def _use_record_transcript(self):
        """Fetch the record's own selected mRNA into the transcript box and flag origin='record' so _render_splice
        reports a consistency check, never 'independent confirmation' (the check is circular by construction)."""
        acc = self.spliceTxSel.currentData() if hasattr(self, "spliceTxSel") else None
        if not acc:
            return
        self.spliceTxUse.setEnabled(False); QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            fasta = retrieve(acc).get("fasta", "")
        except Exception as e:
            QApplication.restoreOverrideCursor(); self.spliceTxUse.setEnabled(True)
            self._banner(f"Could not fetch {acc}: {e.args[0] if getattr(e,'args',None) else e}", "warn"); return
        QApplication.restoreOverrideCursor(); self.spliceTxUse.setEnabled(True)
        self._tx_programmatic = True
        self.spliceTx.setPlainText(fasta)
        self._tx_programmatic = False
        self._splice_tx_origin = "record"                     # same annotation source -> consistency check, not independent
        self.card_splice.expand(); self.spliceTx.setFocus()

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
        who = src.get("displayLocus") or src.get("accession") or "pasted / uploaded specimen"
        org = f" · {src.get('organism')}" if src.get("organism") else ""
        self.spliceRef.setText(f"Genomic reference: <b>{who}</b>{org} · {n:,} bp (from panel 01)")

    def _struct_menu(self, e):
        """Right-click menu for a Structural-evidence row. What is copied must be what the ROW shows: an LTR/TIR
        row lists its two arms (Len = one arm), so each arm is copied on its own instead of silently substituting
        element_span — which is the whole element, and for a TIR the whole input. The whole element stays
        available as its own labelled group of items; the other rows address the single span they display."""
        fp, tp, span = e.get("five_prime"), e.get("three_prime"), e.get("element_span")
        lab = e["type"].split(" ")[0]
        if fp and tp:
            def _grp(name, sp, fid):                      # relabel the shared Copy items so each names its span
                return [(l.replace("Copy ", f"Copy {name} ", 1), f)
                        for l, f in self._feat_menu(sp[0], sp[1], "+", fid, kind=e["type"])]
            items = _grp("5′ arm", fp, f"{lab}_5prime") + _grp("3′ arm", tp, f"{lab}_3prime")
            if span and span[0] is not None:
                items += _grp("whole element", span, f"{lab}_element")
            return items
        sp = e.get("pos") or e.get("upstream") or span or [None, None]   # mirrors _struct_row's coord precedence
        if sp[0] is None:
            return [("Copy type", lambda: self._copy(e["type"]))]
        return self._feat_menu(sp[0], sp[1], "+", lab, kind=e["type"])

    def _region_menu(self, region):
        """Right-click menu for a feature glyph in the genome viewer (copy FASTA / design primer)."""
        return self._feat_menu(region["start"], region["end"], region.get("strand", "+"),
                               region.get("label") or "feature", kind=region.get("kind"))

    def _gel_menu(self, region):
        """Right-click menu for a gel band → copy the amplicon FASTA / coordinates."""
        a = region.get("amplicon") or {}
        pair = region.get("pair", "")
        start, end, seq = a.get("start"), a.get("end"), a.get("seq", "")
        # same call as the band tooltip, the amplicon table and its FASTA export (single-primer is its own
        # bucket; on/off only mean something with a design locus)
        kind = _amp_kind(a, region.get("has_locus", True))
        label = f"amplicon_{pair}_{start}-{end}_{a.get('length','')}bp_{kind}"
        items = [("Copy FASTA", lambda: self._copy(f">{label}\n{seq}"))]
        if start is not None and end is not None:
            items.append((f"Copy coords ({start}–{end})", lambda: self._copy(f"{start}-{end}")))
        return items

    def _src_link(self, key, url=None):
        """A small clickable 'source' citation label that opens the verified DOI in a browser."""
        r = REFLINKS.get(key)
        lab = QLabel()
        if not r:
            return lab
        accent = theme_mod.ACCENT[self.theme]
        lab.setTextFormat(Qt.RichText)
        lab.setText(f'<a href="{url or r["url"]}" style="color:{accent};text-decoration:underline">source</a>')
        lab.setOpenExternalLinks(True)
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
    def _design_block(self, seq=None):
        """ONE guard for BOTH design entry points: one design at a time, a non-empty template, and — when the
        template IS the panel-01 specimen (seq is None) — a non-stale one. The stale test is deliberately not
        applied to an explicit template (amplicon / family hit / splice product): those do not come from the
        specimen box, so the box's edit state says nothing about them and gating on it would block a working path."""
        if self._design_inflight:
            self._banner("A primer design is already running — wait for it to finish."); return True
        if seq is None and self._stale_block():
            return True
        if not self._clean_seq(seq if seq is not None else self.state.get("seq", "")):
            self._banner("No sequence to design a primer on."); return True
        return False

    def _design(self):
        self._clear_banner()
        if self._design_block():
            return
        self._design_inflight = True
        self.designBtn.setEnabled(False); self.designBtn.setText("◴ designing…")
        self._design_tmpl = self.state["seq"]             # remember the template these candidates index (for PCR on-target)
        self.engine.submit("primers", {"sequence": self.state["seq"], "params": self._read_primer_params()}, key="primers")

    def _design_for_domain(self, start, end, label, seq=None):
        self._clear_banner()
        if self._design_block(seq):
            return
        tmpl = seq if seq is not None else self.state.get("seq", "")
        inc = [start, max(60, end - start)]
        self._scroll_to(self.card_primer)                 # bring the primer designer up so the user lands on the next phase
        self.designBtn.setEnabled(False); self.designBtn.setText("◴ designing…")   # same busy cue as the toolbar design path
        self._set_body(self.primBody, BusyBar(f"Designing primers for {label}…"))  # in-flight, not empty (same cue as annotate/splice)
        self._pending_domain = label
        self._design_inflight = True
        self._design_tmpl = tmpl                          # candidates' left/right_pos index THIS template, not always panel-01
        self.engine.submit("primers", {"sequence": tmpl, "params": self._read_primer_params(),
                                        "included": inc}, key="primers")

    def _on_primers(self, d):
        self._design_inflight = False
        self.designBtn.setEnabled(True); self.designBtn.setText("Design primers")
        self.card_primer.expand()
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
        tmpl_sig = self._norm_seq(getattr(self, "_design_tmpl", self.state.get("seq", "")))
        for c in cands:                                       # tag each pair with the template its coords index
            c["_tmpl_sig"] = tmpl_sig
        self.state["candidates"] = cands
        self._refresh_scan_picker()                           # keep card-06's in-card scan picker in sync with designed pairs
        self.pcrStageAll.setEnabled(bool(cands))
        if not cands:
            self.primBody.addWidget(_note("No primer pair met the criteria — try the Permissive preset or widen the product range.", "warn"))
            return
        srow = QHBoxLayout(); srow.addWidget(_sl(f"Primer pairs — {len(cands)}")); srow.addStretch(1)
        srow.addWidget(self._src_link("Primer3")); self.primBody.addLayout(srow)
        headers = ["ID", "Forward (5'→3')", "Reverse (5'→3')", "Product", "Tm F/R", "GC% F/R",
                   "Hairpin", "Self-dim", "Hetero", "3′-end", "Struct", "Penalty"]
        t = DataTable(headers, GLOSS)
        fc = theme_mod.FLAG[self.theme]                       # per-theme flag colours (WCAG-tuned dark vs light)
        rows, styles, tips, exports = [], [], [], []
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
                exports.append([None] * 6 + [hp[3], sd[3], het[3], end[3], None, None])   # ΔG exports bare, marks are display-only
            else:                                             # QC unavailable for this pair (never drops the pair)
                rows.append(base + ["—", "—", "—", "—", "n/a", c["penalty"]])
                styles.append([None] * 12); tips.append([None] * 12); exports.append([None] * 12)
        t.set_rows(rows, styles=styles, tips=tips, exports=exports)
        t.set_row_menu(lambda r: [("→ send to in-silico PCR", lambda rr=r: self._add_pcr_pair(cands[rr])),
                                  ("Secondary-structure detail", lambda rr=r: self._structure_detail(cands[rr])),
                                  ("Scan whole genome for off-targets", lambda rr=r: self._scan_genome(cands[rr])),
                                  ("Copy pair FASTA", lambda rr=r: self._copy(
                                      f">{cands[rr]['id']}_F\n{cands[rr]['left_seq']}\n>{cands[rr]['id']}_R\n{cands[rr]['right_seq']}"))])
        t.setMaximumHeight(round(220 * theme_mod.UI_SCALE))
        pxr = QHBoxLayout(); pxr.addStretch(1); pxr.addWidget(_export_table_btn(t, "TEagle_primers", self))
        self.primBody.addLayout(pxr)                              # above the 12-col table (a beside button scrolls out)
        self.primBody.addWidget(t)
        # visible key for the ΔG colour flags + the ‡ marker + the F/R fold (colour is never the only signal)
        legend = QLabel(f'<span style="color:{fc["caution"]}">■</span>&nbsp;<b>!</b>&nbsp;caution (moderately stable; threshold varies by structure) &nbsp;&nbsp;'
                        f'<span style="color:{fc["warn"]}">■</span>&nbsp;<b>!!</b>&nbsp;warn (ΔG ≤ −9) &nbsp;&nbsp;'
                        f'‡&nbsp;the two engines disagree &nbsp;·&nbsp; Hairpin / Self-dim show the worst of the forward and reverse primer '
                        f'&nbsp;·&nbsp; an exported table carries the bare ΔG number (the marks are on-screen cues)')
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
                      "Scroll the table sideways — or collapse the specimen panel (Ctrl+B) — to see every column.")
        hint.setObjectName("orient"); hint.setWordWrap(True); self.primBody.addWidget(hint)

    def _structure_detail(self, c):
        """IDT OligoAnalyzer-style detail: both engines' ΔG for every structure of one primer pair, side by side."""
        qc = c.get("qc") or {}
        dlg = QDialog(self); dlg.setWindowTitle(f"Secondary structure — {c.get('id','pair')}")
        dlg.resize(round(640 * theme_mod.UI_SCALE), round(460 * theme_mod.UI_SCALE))
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
            t.setMinimumHeight(round(230 * theme_mod.UI_SCALE)); lay.addWidget(t)
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
        copy = QPushButton("Copy"); copy.clicked.connect(lambda: self._copy(self._structure_text(c)))
        copy.setAccessibleName("Copy the secondary-structure ΔG values as text")
        close = QPushButton("Close"); close.clicked.connect(dlg.accept)
        close.setAccessibleName("Close the secondary-structure detail window")
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
        self._pending_domain = None                       # a newer explicit navigation cancels any pending routed re-anchor
        self._scroll_to(self.card_pcr)                    # send-to-PCR lands the user on the in-silico PCR panel

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
            lab = QLabel(f"P{i+1}  {c['left_seq'][:16]}… / {c['right_seq'][:16]}…  · {c['product_size']} bp")
            lab.setObjectName("cardmeta"); rl.addWidget(lab); rl.addStretch(1)
            up = QPushButton("Up"); up.setProperty("sm", True); up.clicked.connect(lambda _=False, k=i: self._move_pair(k, -1))
            dn = QPushButton("Down"); dn.setProperty("sm", True); dn.clicked.connect(lambda _=False, k=i: self._move_pair(k, 1))
            rm = QPushButton("Remove"); rm.setProperty("sm", True); rm.clicked.connect(lambda _=False, k=i: self._remove_pair(k))
            rl.addWidget(up); rl.addWidget(dn); rl.addWidget(rm)
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
        scanned_asm = (all_assemblies().get(org, {}) or {}).get("assemblyAccession", "")
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
        acc = (all_assemblies().get(org, {}) or {}).get("assemblyAccession", "")
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
            self._render_genome_status(f"{d.get('organism','')} genome downloaded — right-click a pair to scan it.", "success")
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
            picker.addItem(f"{c['id']}  ·  {c['left_seq'][:12]}… / {c['right_seq'][:12]}…", c)
        picker.blockSignals(False)
        has_pairs = bool(cands)
        has_genome = getattr(self, "genomeOrg", None) is not None and self.genomeOrg.count() > 1   # past the placeholder
        picker.setEnabled(has_pairs)
        self.scanBtn.setEnabled(has_pairs and has_genome)
        self.scanBtn.setToolTip("Design a primer pair in panel 04 first" if not has_pairs
                                else "Download a genome via Manage genomes first" if not has_genome else "")

    def _rebuild_coord_asm_dropdown(self):
        """(Re)populate the coordinate-fetch organism dropdown from all_assemblies(), so a newly added custom
        organism appears in-session (preserving the current selection). The '__custom__' free-text item stays last."""
        box = getattr(self, "asmSel", None)
        if box is None:
            return
        prev = box.currentData()
        box.blockSignals(True); box.clear()
        _asm = all_assemblies()
        for org in sorted(_asm):
            box.addItem(f"{org} · {_asm[org]['assemblyName']}", org)
        box.addItem("Other organism / assembly…", "__custom__")
        if prev is not None:
            i = box.findData(prev)
            if i >= 0:
                box.setCurrentIndex(i)
        box.blockSignals(False)

    def _refresh_genome_dropdown(self):
        """Rebuild the PCR organism dropdown from the prepared (downloaded+verified) genome set only —
        an organism cannot be scanned until it is downloaded. userData stays the ORGANISM string (what
        _scan_genome and the need_prepare fallback consume via currentData()); a downloaded accession
        outside the curated map is skipped. Selection is preserved across the rebuild."""
        box = getattr(self, "genomeOrg", None)
        if box is None:                                       # genome_list can resolve before the PCR card is built
            return
        prev = box.currentData()
        asm = all_assemblies()                                # curated + user-added, so a downloaded custom genome shows here too
        acc2org = {v["assemblyAccession"]: k for k, v in asm.items()}
        box.blockSignals(True)
        box.clear()
        box.addItem("— select organism —", None)              # placeholder so no wrong-species scan runs by default
        n = 0
        for g in self._prepared_genomes:
            org = acc2org.get(g.get("accession"))
            if not org:
                continue
            box.addItem(f"{org} · {asm[org]['assemblyName']}", org)
            n += 1
        if prev is not None:
            i = box.findData(prev)
            if i >= 0:
                box.setCurrentIndex(i)
        box.blockSignals(False)
        hint = getattr(self, "genomeOrgHint", None)
        if hint is not None:
            hint.setText("" if n else "No genomes downloaded yet — open Manage genomes to download one for scanning.")
        self._refresh_scan_picker()                           # genome availability gates the scan button

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
        self._rebuild_coord_asm_dropdown()                    # a newly added custom organism must also appear in coord-fetch

        dlg = getattr(self, "_genome_mgr", None)
        if dlg is None:
            if not getattr(self, "_genome_mgr_open", False):  # a refresh that landed AFTER the user closed the manager
                return                                        # -> only skip the DIALOG; the dropdown is already updated
            dlg = QDialog(self); dlg.setAttribute(Qt.WA_DeleteOnClose)
            dlg.finished.connect(lambda *_: (setattr(self, "_genome_mgr", None),
                                             setattr(self, "_genome_mgr_open", False)))
            self._genome_mgr = dlg
        from PySide6.QtGui import QColor
        dlg.setWindowTitle("Manage genomes")
        z = theme_mod.UI_SCALE
        try:                                                  # the dialog is sized to fit the screen + fixed (non-resizable) at the end
            scrg = self.screen().availableGeometry(); scrw, scrh = scrg.width(), scrg.height()
        except Exception:
            scrw, scrh = 1440, 900
        old = dlg.layout()
        if old is not None:                                   # rebuild the body in place on refresh
            dlg.hide()                                        # hide while swapping the layout, else the re-added widgets measure as hidden and sizeHint()
            QWidget().setLayout(old)                          # collapses to ~20px -> setFixedSize would pin the whole dialog to a sliver (the final dlg.show() re-shows it)
        lay = QVBoxLayout(dlg)
        lay.setContentsMargins(theme_mod.sp(12), theme_mod.sp(10), theme_mod.sp(12), theme_mod.sp(10)); lay.setSpacing(theme_mod.sp(7))
        title = QLabel("Manage genomes"); tf = title.font(); tf.setBold(True); tf.setPointSizeF(tf.pointSizeF() + 2)
        title.setFont(tf); lay.addWidget(title)
        if not d.get("ok"):                                   # WSL-down error card (styled, not a bare label)
            err = QLabel("Could not list genomes — " + d.get("error", "WSL backend unavailable"))
            err.setObjectName("errbanner"); err.setWordWrap(True); lay.addWidget(err)
            close = QPushButton("Close"); close.clicked.connect(dlg.accept)
            close.setAccessibleName("Close the genome manager"); lay.addWidget(close)
            # this branch returns BEFORE the table fit-tail that owns setFixedWidth/adjustSize/setFixedSize, and the
            # dialog is REUSED — so release the size a previous successful build pinned on it, else the error card is
            # stranded inside the old table's geometry with a screenful of dead space.
            self._uppercase_buttons(dlg)                  # uppercase BEFORE the fit tail measures the dialog
            ew = min(round(560 * z), scrw - 56)
            dlg.setMinimumSize(0, 0); dlg.setMaximumSize(16777215, 16777215)   # QWIDGETSIZE_MAX
            dlg.setFixedWidth(ew); dlg.adjustSize()
            dlg.setFixedSize(ew, min(dlg.sizeHint().height(), scrh - 56))
            dlg.show(); dlg.raise_(); return
        intro = QLabel("Download a genome once to enable its whole-genome off-target scan; it is kept locally and then "
                       "appears in the organism menus. Mammalian genomes are large (~1 GB) and the download needs "
                       "≥8 GB free disk in WSL: extraction peaks near 4 GB (FASTA ~3 GB + 2bit ~0.8 GB) before the "
                       "temporary files are removed.")
        intro.setObjectName("orient"); intro.setWordWrap(True); lay.addWidget(intro)
        # lock every action while any download/scan/annotation runs. The annotation belongs in this guard:
        # it reads the cached genome for hours, and Delete does `rm -rf` on that same directory.
        busy = self._genome_prep_inflight or self._genome_inflight or self._annot_inflight
        # ADD-ORGANISM row — one input, auto-detects an organism name vs a GCF/GCA accession; resolved + pinned once.
        arow = QHBoxLayout()
        self._addAsmEdit = QLineEdit()
        self._addAsmEdit.setPlaceholderText("Add an organism by name or assembly accession — e.g. Danio rerio or GCF_000002035.6")
        addBtn = QPushButton("Add"); addBtn.setProperty("sm", True); addBtn.setEnabled(not busy and not self._add_asm_inflight)
        addBtn.setAccessibleName("Add the typed organism or assembly accession to the genome list")
        addBtn.clicked.connect(self._add_custom_assembly); self._addAsmEdit.returnPressed.connect(self._add_custom_assembly)
        self._addAsmBtn = addBtn                               # kept so the add handlers can disable/re-enable it across a dialog rebuild
        if self._add_asm_inflight:                            # a resolve started in a since-closed dialog is still running
            self._addAsmEdit.setEnabled(False)
        arow.addWidget(self._addAsmEdit, 1); arow.addWidget(addBtn); lay.addLayout(arow)
        prepared = {g["accession"]: g for g in d.get("genomes", [])}
        asm = all_assemblies(); orgs = sorted(asm)            # curated + user-added
        pal = theme_mod._DARK if self.theme == "dark" else theme_mod._LIGHT
        rows = len(orgs)
        # size every row so the WHOLE table (all rows + columns + buttons) shows with NO scrolling in a fixed dialog.
        # rows stay comfortable on a normal screen and shrink (with the font) only if a small screen forces it. The final
        # row height is computed in the tail after the real chrome (title+intro+add-row+close) is measured at the true width.
        avail_h = scrh - round(56); avail_w = scrw - round(56)
        # ONE fixed, comfortable row height, chosen so the natural-height 'sm' button + a small margin always fit (the
        # button is NEVER forced short or shrunk -> never squeezed/clipped). Set BEFORE the cells are built so each cell
        # widget is laid out at the right height. A screen too small to show every row scrolls vertically (button stays
        # comfortable) rather than squeezing the button — a readable button matters more than avoiding a scroll.
        rowH = round(34 * z)
        target_w = min(round(1180 * z), avail_w)              # bigger, so full organism + assembly names show without eliding
        headers = ["Organism", "Assembly", "Accession", "Status", "Size (MB)", "Contigs", "Action"]
        tbl = QTableWidget(len(orgs), len(headers))
        tbl.setHorizontalHeaderLabels(headers)
        tbl.verticalHeader().setVisible(False)
        tbl.setEditTriggers(QAbstractItemView.NoEditTriggers)
        tbl.setSelectionMode(QAbstractItemView.NoSelection)
        tbl.setSortingEnabled(False)                          # static catalog; sorting conflicts with setCellWidget buttons
        tbl.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        # Organism AND Assembly (the two widest text columns) stretch, so on a narrow dialog THEY shrink+elide while the
        # short fixed columns — including Action — keep their content width. The Download button can never be the column
        # that overflows off the right edge into a horizontal scrollbar (the small-screen clip a single stretch column left).
        tbl.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        tbl.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
        tbl.setTextElideMode(Qt.ElideRight)                   # a shrunk Organism/Assembly elides ("Drosophila melan…"), never clips Action
        tbl.setWordWrap(False)                                # (h-scrollbar stays as-needed — a pathological tiny screen keeps Action reachable by scroll rather than silently clipped)
        tbl.horizontalHeader().setStretchLastSection(False)
        tbl.horizontalHeader().setMinimumSectionSize(round(56 * z))
        tbl.horizontalHeader().setDefaultAlignment(Qt.AlignCenter)   # header + cells both centred (they disagreed before)
        tbl.verticalHeader().setDefaultSectionSize(rowH)      # provisional row height; the real value is set in the tail
        aligns = [Qt.AlignHCenter] * 6                        # every column centred, matching the centred header
        row_btns = []                                         # collected so every row's button gets ONE uniform size
        for r, org in enumerate(orgs):
            meta = asm[org]; acc = meta["assemblyAccession"]; g = prepared.get(acc); downloaded = bool(g)
            cells = [org, meta.get("assemblyName", ""), acc, "Downloaded" if downloaded else "Not downloaded",
                     (f"{(g.get('bytes', 0) or 0) / 1e6:.0f}" if downloaded else "—"),
                     (str(g.get("n_seqs", "?")) if downloaded else "—")]
            for c, txt in enumerate(cells):
                it = QTableWidgetItem(txt); it.setTextAlignment(aligns[c] | Qt.AlignVCenter)
                if c in (0, 1):                               # Organism/Assembly may elide on a narrow dialog — tooltip recovers the full text
                    it.setToolTip(txt)
                if c == 3:                                    # status carried by colour + word (green = downloaded)
                    it.setForeground(QColor(pal["good"] if downloaded else pal["faint"]))
                tbl.setItem(r, c, it)
            btn = QPushButton("Delete" if downloaded else "Download")
            btn.setProperty("cellbtn", True)                  # the ONLY style property on it (see theme.py) -> no specificity clash
            btn.ensurePolished()                              # apply the stylesheet BEFORE sizeHint is read, else it measures the base font
            btn.clicked.connect((lambda _=False, a=acc, o=org: self._remove_genome(a, o)) if downloaded
                                else (lambda _=False, o=org: self._on_manager_download(o)))
            btn.setEnabled(not busy)
            holder = QWidget()                                # full-cell holder; AlignCenter centres the button both ways in whatever
            hl = QHBoxLayout(holder)                          # geometry the view hands the cell -> no top/bottom bias, no clipping
            hl.setContentsMargins(0, 0, 0, 0); hl.setSpacing(0)
            hl.addWidget(btn, 0, Qt.AlignCenter)
            tbl.setCellWidget(r, 6, holder)
            row_btns.append(btn)
        # ONE uniform button size for every row, from the widest label ("Download" > "Delete"), so the buttons align
        # with each other instead of each sizing to its own text. Fixed size => the holder can never squeeze or stretch it.
        # Qt applies QTableWidget::item padding to CELL WIDGETS as well as text: it offsets the widget by the
        # padding's top-left AND shrinks its usable box by twice the padding. So the row/column must be sized to
        # (button + 2x padding), or the button is both pushed off-centre and sheared. Read the padding from the
        # theme so this can never drift from the stylesheet.
        pv, ph = round(theme_mod.TABLE_ITEM_PAD_V * z), round(theme_mod.TABLE_ITEM_PAD_H * z)
        if row_btns:
            for b in row_btns:                                # _uppercase_buttons() runs LATER (any theme/scale/result
                b.setText(b.text().upper())                   # refresh reaches this dialog) — measure the FINAL text or
            bw = max(b.sizeHint().width() for b in row_btns)  # the pinned width is the sentence-case one and "DOWNLOAD"
            bh = max(b.sizeHint().height() for b in row_btns) # gets clipped to "OWNLOA"
            for b in row_btns:
                b.setFixedSize(bw, bh)                        # uniform size => the buttons align with each other
            rowH = max(rowH, bh + 2 * pv + round(4 * z))      # inset height (rowH - 2*pv) always exceeds the button
            tbl.verticalHeader().setDefaultSectionSize(rowH)
        tbl.resizeColumnsToContents()
        if row_btns:
            # resizeColumnsToContents() measures ITEMS only — it is blind to setCellWidget widgets, which is why the
            # Action column used to size to the word "Action" and let the button spill past the table's right edge.
            tbl.horizontalHeader().setSectionResizeMode(6, QHeaderView.Fixed)
            tbl.setColumnWidth(6, bw + 2 * ph + round(6 * z))
        for r in range(rows):                                 # the view owns the holder's geometry (cell minus the
            tbl.setRowHeight(r, rowH)                          # symmetric item padding) -> AlignCenter is true centring
        lay.addWidget(tbl)
        self._genome_mgr_table = tbl
        if busy:
            b = QLabel("A download or scan is running — actions are disabled until it finishes.")
            b.setObjectName("orient"); b.setWordWrap(True); lay.addWidget(b)
        close = QPushButton("Close"); close.setProperty("sm", True); close.clicked.connect(dlg.accept)
        close.setAccessibleName("Close the genome manager")
        # Whole-genome TE annotation lives in the footer rather than a per-row button: the Action column's
        # width and row height are pinned from ONE measured button, and a second per-row control would
        # reopen the squeeze/clipping problems that sizing was written to prevent.
        annot = QPushButton("Annotate TE landscape…"); annot.setProperty("sm", True)
        annot.setAccessibleName("Annotate every transposable element in a downloaded genome")
        annot.setToolTip("Run RepeatMasker over a downloaded genome and summarise its transposable-element content")
        annot.setEnabled(bool(prepared) and not busy and not self._annot_inflight)
        annot.clicked.connect(self._annotate_genome_dialog)
        crow = QHBoxLayout(); crow.addStretch(1); crow.addWidget(annot); crow.addWidget(close)
        crow.addStretch(1); lay.addLayout(crow)
        self._uppercase_buttons(dlg)                          # ADD/CLOSE match the row DOWNLOAD/DELETE; measured below
        # ---- fit: full table if it fits the screen, else cap to whole rows + a vertical scroll (button never squeezed) ----
        dlg.setFixedWidth(target_w)                           # set the real width first so the intro wraps as it truly will
        hdr_h = tbl.horizontalHeader().sizeHint().height() or round(34 * z)
        tbl.setFixedHeight(hdr_h + rows * rowH + round(4 * z))
        tbl.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        dlg.adjustSize()
        if dlg.sizeHint().height() > avail_h:                 # too tall for this screen -> show as many WHOLE rows as fit + scroll
            over = dlg.sizeHint().height() - avail_h
            vis = max(1, (rows * rowH - over) // rowH)
            tbl.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)
            tbl.setFixedHeight(hdr_h + vis * rowH + round(4 * z))
            dlg.adjustSize()
        dlg.setFixedSize(target_w, min(dlg.sizeHint().height(), avail_h))
        dlg.show(); dlg.raise_()

    # =================== whole-genome TE annotation ===================
    # Wording taken verbatim from RepeatMasker's own -h text, not paraphrased: an earlier draft turned
    # "5-10% less sensitive" into "5-10x faster", which would have misled the user about both axes.
    _ANNOT_LOW_LIBRARY = 25          # see _on_annotate_budget for how this was chosen (measured, not assumed)
    _ANNOT_SENS = [("default", "Default — the balance RepeatMasker ships with"),
                   ("quick", "Quick (-q) — 2–5× faster, 5–10% less sensitive"),
                   ("slow", "Slow (-s) — 2–3× slower, 0–5% more sensitive")]

    def _annotate_genome_dialog(self):
        """Choose a downloaded genome and the run's parameters, and state the cost BEFORE anything starts.

        The cost panel is the scientific gate as much as the practical one: it reports how many family
        models the installed Dfam partitions actually hold for the chosen lineage. A lineage with none
        cannot yield a TE result at any runtime, and the user is told that here rather than after hours."""
        if self._annot_inflight:
            self._banner("A genome annotation is already running.", "info"); return
        # _prepared_genomes is the app's single store of what is cached on disk (kept current by the one
        # genome_list fan-out), so the picker can never disagree with the manager table beside it.
        prepared = [g for g in (getattr(self, "_prepared_genomes", None) or []) if g.get("accession")]
        if not prepared:
            self._banner("Download a genome first — annotation runs on a genome already on this machine.", "info")
            return
        asm = all_assemblies()
        acc_to_org = {m["assemblyAccession"]: o for o, m in asm.items()}
        dlg = QDialog(self); dlg.setWindowTitle("Annotate TE landscape"); dlg.setModal(False)
        dlg.setAttribute(Qt.WA_DeleteOnClose, True)
        lay = QVBoxLayout(dlg)
        lay.setContentsMargins(theme_mod.sp(14), theme_mod.sp(12), theme_mod.sp(14), theme_mod.sp(12))
        lay.setSpacing(theme_mod.sp(8))
        t = QLabel("Annotate the TE landscape of a downloaded genome")
        tf = t.font(); tf.setBold(True); tf.setPointSizeF(tf.pointSizeF() + 1); t.setFont(tf); lay.addWidget(t)
        intro = QLabel("TEagle runs RepeatMasker over the whole assembly and summarises which transposable-element "
                       "families it places, how much of the genome they cover, and how diverged each is. This finds "
                       "copies of families that are already in the installed Dfam library; it does not discover new "
                       "families, so coverage depends on that library.")
        intro.setObjectName("orient"); intro.setWordWrap(True); lay.addWidget(intro)

        grid = QHBoxLayout()
        gsel = QComboBox()
        for g in prepared:
            org = acc_to_org.get(g["accession"], g["accession"])
            gsel.addItem(f"{org} · {g['accession']} · {(g.get('bytes') or 0)/1e6:.0f} MB", g["accession"])
        gsel.setToolTip("A genome already downloaded to this machine. Only downloaded genomes can be "
                        "annotated — the whole assembly is read locally, nothing is sent to a server.")
        _gl = QLabel("Genome"); _gl.setToolTip(gsel.toolTip())
        grid.addWidget(_gl); grid.addWidget(gsel, 1)
        ssel = QComboBox()
        for key, label in self._ANNOT_SENS:
            ssel.addItem(label, key)
        ssel.setToolTip("How hard RepeatMasker looks. The figures are RepeatMasker's own: Quick (-q) is "
                        "2–5× faster and 5–10% less sensitive; Slow (-s) is 2–3× slower and 0–5% more "
                        "sensitive. Sensitivity here means recovering older, more diverged copies — the "
                        "recent ones are found either way.")
        _sl2 = QLabel("Sensitivity"); _sl2.setToolTip(ssel.toolTip())
        grid.addWidget(_sl2); grid.addWidget(ssel, 1)
        lay.addLayout(grid)
        # LIBRARY CHOICE — the single setting that decides which families can be found at all, so it is
        # offered explicitly rather than assumed. The two options are alternatives, not additions:
        # RepeatMasker's -lib replaces the database search, so a run uses one or the other.
        lrow = QHBoxLayout()
        lsel = QComboBox()
        lsel.addItem("Installed Dfam library — curated families only", "")
        lsel.addItem("Installed Dfam library — include uncurated families", "uncurated")
        lsel.addItem("My own repeat library (FASTA)…", "custom")
        lsel.setToolTip(
            "Dfam families are either curated (reviewed) or uncurated (automatically built). RepeatMasker "
            "searches curated families only unless you ask for both.\n\n"
            "Outside a few heavily studied species most families are uncurated, so curated-only can find "
            "nothing at all: on baker's yeast it sees 9 families and reports no transposable element, while "
            "including uncurated sees 421 more and finds the Ty elements the genome really has.\n\n"
            "Including them needs the optional uncurated partitions from the backend installer, and the "
            "extra families are lower confidence — the result says which setting produced it.\n\n"
            "A custom FASTA library replaces Dfam entirely for that run.")
        lrow.addWidget(QLabel("Library")); lrow.addWidget(lsel, 1)
        libpath = QLabel(""); libpath.setObjectName("orient"); libpath.setWordWrap(True)
        lrow.addWidget(libpath, 2); lay.addLayout(lrow)
        self._annot_custom_lib = None
        self._annot_uncurated = False

        def pick_library(idx):
            self._annot_uncurated = (lsel.currentData() == "uncurated")
            if lsel.currentData() != "custom":
                self._annot_custom_lib = None
                libpath.setText("Curated + uncurated families. Needs the uncurated Dfam partitions installed; "
                                "uncurated families are automatically built and lower confidence."
                                if self._annot_uncurated else "")
                return
            fn, _ = QFileDialog.getOpenFileName(dlg, "Choose a repeat library (FASTA)", "",
                                                "FASTA library (*.fa *.fasta *.lib *.txt);;All files (*)")
            if not fn:
                lsel.setCurrentIndex(0); self._annot_custom_lib = None; libpath.setText(""); return
            self._annot_custom_lib = fn
            libpath.setText(f"{os.path.basename(fn)} — searched instead of Dfam; TEagle records its checksum "
                            "but cannot vouch for its contents.")
        lsel.currentIndexChanged.connect(pick_library)
        cost = QLabel("Estimating…"); cost.setObjectName("orient"); cost.setWordWrap(True)
        cost.setTextInteractionFlags(Qt.TextSelectableByMouse); lay.addWidget(cost)
        warn = QLabel(""); warn.setObjectName("errbanner"); warn.setWordWrap(True); warn.hide(); lay.addWidget(warn)
        start = QPushButton("Start annotation"); start.setEnabled(False)
        cancel = QPushButton("Close"); cancel.setProperty("sm", True); cancel.clicked.connect(dlg.reject)
        brow = QHBoxLayout(); brow.addStretch(1); brow.addWidget(start); brow.addWidget(cancel); lay.addLayout(brow)
        self._annot_cost_state = {}

        def refresh_cost():
            acc = gsel.currentData()
            g = next((x for x in prepared if x["accession"] == acc), {})
            org = acc_to_org.get(acc, "")
            cost.setText("Estimating cost — reading the WSL budget and the installed Dfam library…")
            start.setEnabled(False); warn.hide()
            self._annot_cost_state = {"accession": acc, "organism": org, "sha256": g.get("sha256")}
            self.engine.submit("annotate_budget",
                               {"species": org, "genome_bytes": int(g.get("bytes") or 0)}, key="annotate_budget")

        self._annot_cost_widgets = (cost, warn, start, gsel, ssel)
        gsel.currentIndexChanged.connect(lambda *_: refresh_cost())
        start.clicked.connect(lambda: (self._start_genome_annotate(gsel.currentData(), acc_to_org.get(gsel.currentData(), ""),
                                                                   ssel.currentData(), self._annot_custom_lib,
                                                                   getattr(self, "_annot_uncurated", False)), dlg.accept()))
        self._uppercase_buttons(dlg)
        refresh_cost()
        dlg.adjustSize(); dlg.show(); dlg.raise_()
        self._annot_dialog = dlg

    def _on_annotate_budget(self, b):
        """Fill the cost panel from the measured WSL budget + the library's coverage for this lineage."""
        w = getattr(self, "_annot_cost_widgets", None)
        if not w:
            return
        cost, warn, start, gsel, ssel = w
        try:
            st = self._annot_cost_state
            # keep the MEASURED thread count so the run uses the number the dialog promised
            self._annot_threads = int(b.get("recommended_threads") or 4)
            fams = b.get("library_families_for_species")
            gb = b.get("disk_needed_gb")
            lines = [
                f"<b>This machine (WSL):</b> {b.get('cores')} cores · {b.get('mem_gb')} GB RAM · "
                f"{b.get('avail_gb')} GB free. TEagle will use {b.get('recommended_threads')} parallel "
                f"RepeatMasker jobs (limited by {b.get('limited_by')}).",
                f"<b>Disk needed:</b> about {gb} GB — the genome sequence is fetched once and kept, and the "
                "work chunks are removed as they finish.",
                "<b>Time:</b> a fly-sized genome (~140 Mb) takes minutes to about an hour here; a mammalian "
                "genome (~3 Gb) takes many hours to a day. The run continues in the background and can be "
                "re-started later — finished chunks are not repeated.",
                f"<b>Library coverage:</b> the installed Dfam partitions hold "
                f"<b>{fams if fams is not None else 'an unknown number of'}</b> family models for "
                f"{st.get('organism') or 'this lineage'}.",
            ]
            cost.setText("<br>".join(lines))
            if b.get("disk_ok") is False:
                warn.setText(f"Not enough free disk in WSL — about {gb} GB is needed. Free space, or remove a "
                             "cached genome from this manager, then try again.")
                warn.show(); start.setEnabled(False); return
            # Threshold rationale: measured on this backend, a lineage Dfam covers usefully has hundreds
            # of models (D. melanogaster 399, H. sapiens 1439) while an uncovered one has single digits
            # (S. cerevisiae 9, and its run found zero TEs). 25 sits in the empty gap between those two
            # regimes; it is a warning trigger only — it never blocks a run or changes a result.
            if fams is not None and fams < self._ANNOT_LOW_LIBRARY:
                warn.setText(f"The installed Dfam library holds only {fams} family model(s) for "
                             f"{st.get('organism') or 'this lineage'}. A run would report tandem and "
                             "low-complexity repeats but could not find transposable elements, however long it "
                             "ran — that would be a limit of the library, not a property of the genome. Install "
                             "the additional Dfam partitions from the backend installer before relying on this.")
                warn.show()
            start.setEnabled(True)
        except RuntimeError:
            pass

    def _start_genome_annotate(self, acc, organism, sensitivity, custom_lib=None, include_uncurated=False):
        if self._annot_inflight or not acc:
            return
        self._annot_inflight = True
        self._banner(f"Annotating {organism or acc} — this runs in the background; progress appears here.", "info")
        self._refresh_genome_manager()
        self.engine.submit("genome_annotate",
                           {"assemblyAccession": acc, "species": organism, "sensitivity": sensitivity,
                            # the thread count the budget probe MEASURED for this machine — the cost dialog
                            # showed it, so the run must actually use it rather than a fixed default
                            "threads": int(getattr(self, "_annot_threads", 0) or 4),
                            "chunk_mb": int(getattr(self, "_annot_chunk_mb", 0) or 40),
                            "custom_library": custom_lib, "include_uncurated": bool(include_uncurated),
                            "sha256": (self._annot_cost_state or {}).get("sha256")},
                           key="genome_annotate")
        self._annot_timer = QTimer(self); self._annot_timer.setInterval(6000)
        self._annot_timer.timeout.connect(lambda: self.engine.submit("genome_annotate_log", {"tail": 1},
                                                                     key="genome_annotate_log"))
        self._annot_timer.start()

    def _on_genome_annotate_log(self, d):
        line = (d or {}).get("log") or ""
        if line and self._annot_inflight:
            self._set_status_line(f"Genome annotation · {line}")

    def _set_status_line(self, text):
        lbl = getattr(self, "_annotStatus", None)
        if lbl is not None:
            try:
                lbl.setText(text)
            except RuntimeError:
                pass

    def _on_genome_annotate(self, r):
        self._annot_inflight = False
        t = getattr(self, "_annot_timer", None)
        if t is not None:
            t.stop()
        # The backend reports a refused or failed run as {"ok": False, ...} — a normal RETURN, not an
        # exception — so it arrives here on the success signal. Without this guard the failure was
        # announced as a finished annotation and then crashed the results window on the missing fields.
        if not r or not r.get("ok"):
            msg = (r or {}).get("error") or "the annotation did not complete"
            self._refresh_genome_manager()
            # A settings clash is the one failure the user can clear, so offer the action instead of
            # telling them to do something the window gives them no way to do.
            if "different settings" in msg:
                acc = (self._annot_cost_state or {}).get("accession")
                box = QMessageBox(self)
                box.setWindowTitle("Previous annotation found")
                box.setIcon(QMessageBox.Question)
                box.setText("This genome has a part-finished annotation that used different settings.")
                box.setInformativeText(
                    "Chunks already finished were searched with the earlier settings, so they cannot be "
                    "combined with new ones. Discard that unfinished work and start again with the "
                    "settings you just chose?")
                discard = box.addButton("Discard and start over", QMessageBox.AcceptRole)
                box.addButton("Keep it", QMessageBox.RejectRole)
                box.exec()
                if box.clickedButton() is discard and acc:
                    self.engine.submit("genome_annotate_reset", {"assemblyAccession": acc},
                                       key="genome_annotate_reset")
                return
            self._banner("Genome annotation did not finish — " + msg, "warn")
            return

    def _on_genome_annotate_reset(self, r):
        if r.get("ok"):
            self._banner("Previous annotation discarded. Open “Annotate TE landscape” to start again.",
                         "success")
        else:
            self._banner(r.get("error") or "Could not clear the previous annotation.", "warn")
        self._annot_result = r
        self._refresh_genome_manager()
        if r.get("coverage_warning"):
            self._banner(r["coverage_warning"], "warn")
        else:
            self._banner(f"Genome annotation finished — transposable elements cover "
                         f"{r.get('te_percent')}% of the assembly across {r.get('te_family_count')} families.",
                         "success")
        self._show_annot_window(r)

    _ANNOT_COLS = [("Family", "The repeat family RepeatMasker assigned (class/subfamily)"),
                   ("Kind", "TE = transposable element; tandem = simple/low-complexity/satellite; other = "
                            "rRNA, tRNA and unclassified repeats"),
                   ("Copies", "How many separate hits were placed"),
                   ("Bases", "Total bases covered by those hits"),
                   ("% genome", "Share of the whole assembly covered"),
                   ("Mean divergence %", "RepeatMasker's RAW percent mismatch to the family consensus, "
                                         "averaged weighted by hit length. It is NOT Kimura- or "
                                         "CpG-corrected, so it saturates for old copies and is inflated "
                                         "where CpG sites decay fast. Higher means more decayed, but it "
                                         "is a similarity measure, not an age.")]

    @staticmethod
    def _annot_report_md(r):
        """The landscape as a self-contained Markdown report.

        Every hedge the window shows travels with it: what was searched, what the two percentages mean,
        and why an empty TE result may be a library limit. An export that carried only the table would
        read as an unqualified statement about the genome."""
        fams = r.get("families") or []
        te = [f for f in fams if f.get("kind") == "TE"]
        L = [f"# Transposable-element landscape — {r.get('species')} ({r.get('accession')})", "",
             f"- **Transposable elements: {r.get('te_percent')}%** of "
             f"{(r.get('genome_bp') or 0)/1e6:.1f} Mb, across {r.get('te_family_count')} families",
             f"- All repeat classes together: {r.get('masked_percent')}% "
             f"({r.get('total_hits'):,} alignment rows)",
             f"- Searched with RepeatMasker {r.get('repeatmasker_version')} against "
             f"{r.get('library_kind')} (Dfam {r.get('dfam_version')}); "
             f"{r.get('library_families_for_species')} family models available for this lineage",
             f"- Sensitivity {r.get('sensitivity')} · {r.get('chunks')} chunks · "
             f"complete: {bool(r.get('complete'))}", ""]
        if r.get("coverage_warning"):
            L += ["> **Coverage warning.** " + r["coverage_warning"], ""]
        L += ["## Scope and limits", "",
              "- Copies of families already present in the searched library, placed by homology. New",
              "  families are not discovered, and a family absent from the library cannot appear here —",
              "  absence is not evidence of absence.",
              "- Transposable elements are counted separately from tandem repeats (simple repeats,",
              "  low-complexity sequence, satellites) and non-TE entries; only the first is a statement",
              "  about transposable-element content.",
              "- Coverage is computed on merged intervals, so overlapping alignments are not double-counted.",
              "- Divergence is RepeatMasker's raw percent mismatch to the consensus, length-weighted. It is",
              "  not Kimura- or CpG-corrected and is not an age.",
              "- Percentages are of the whole assembly, including unplaced scaffolds.", "",
              "## Families", "",
              "| Family | Kind | Copies (alignment rows) | Bases (merged) | % genome | Mean divergence % |",
              "|---|---|---|---|---|---|"]
        for f in fams:
            L.append(f"| {f.get('family')} | {f.get('kind')} | {f.get('n'):,} | {f.get('bp'):,} | "
                     f"{f.get('percent')} | {f.get('divergence')} |")
        L += ["", f"_{len(te)} of {len(fams)} rows are transposable-element families._"]
        return "\n".join(L)

    def _annot_row_menu(self, f):
        """Context menu for one genome-landscape row — aggregate-appropriate actions only."""
        kind = f.get("kind")
        items = [(f"Copy family name ({f.get('family')})", lambda: self._copy(str(f.get("family")))),
                 ("Copy this row (TSV)", lambda: self._copy("\t".join(
                     str(x) for x in (f.get("family"), f.get("kind"), f.get("n"), f.get("bp"),
                                      f.get("percent"), f.get("divergence")))))]
        if kind != "TE":
            items.append(("Why is this not counted as a TE?",
                          lambda: self._banner(
                              "Simple repeats, low-complexity sequence and satellites are tandem repeats, and "
                              "rRNA/tRNA/unclassified entries are not transposable elements. RepeatMasker places "
                              "them alongside TEs, so TEagle counts them separately — they are included in the "
                              "all-repeat percentage but never in the transposable-element percentage.", "info")))
        else:
            items.append(("What does mean divergence mean?",
                          lambda: self._banner(
                              "Average percent difference between the placed copies and their family consensus, "
                              "weighted by how long each copy is. Higher means older, more decayed copies; lower "
                              "means recent or still-active ones. It is a homology measure, not an age in years.",
                              "info")))
        return items

    def _show_annot_window(self, r):
        """The genome-wide result, in its own window: what was found, what was searched, and what that does
        and does not license the reader to conclude."""
        w = QDialog(self); w.setWindowTitle("Genome TE landscape"); w.setModal(False)
        w.setAttribute(Qt.WA_DeleteOnClose, True)
        lay = QVBoxLayout(w)
        lay.setContentsMargins(theme_mod.sp(14), theme_mod.sp(12), theme_mod.sp(14), theme_mod.sp(12))
        lay.setSpacing(theme_mod.sp(8))
        head = QLabel(f"{r.get('species')} · {r.get('accession')}")
        hf = head.font(); hf.setBold(True); hf.setPointSizeF(hf.pointSizeF() + 2); head.setFont(hf); lay.addWidget(head)
        # headline numbers: the TE share is stated separately from every repeat, because they are different
        # claims and merging them would overstate transposable-element content.
        tot = QLabel(f"<b>Transposable elements: {r.get('te_percent')}%</b> of "
                     f"{(r.get('genome_bp') or 0)/1e6:.1f} Mb, in {r.get('te_family_count')} families &nbsp;·&nbsp; "
                     f"all repeat classes together: {r.get('masked_percent')}% &nbsp;·&nbsp; "
                     f"{r.get('total_hits'):,} hits")
        tot.setTextFormat(Qt.RichText); tot.setWordWrap(True)
        tot.setToolTip("Transposable elements = the share of the assembly covered by TE families (LTR, LINE, "
                       "SINE, DNA transposons, rolling-circle, retroposons).\n"
                       "All repeat classes = that plus tandem repeats (simple repeats, low-complexity "
                       "sequence, satellites) and non-TE entries such as rRNA and tRNA genes.\n"
                       "The two are shown separately because only the first is a statement about "
                       "transposable-element content.")
        lay.addWidget(tot)
        prov = QLabel(f"Searched with RepeatMasker {r.get('repeatmasker_version')} against Dfam "
                      f"{r.get('dfam_version')} ({r.get('library_families_for_species')} family models available "
                      f"for this lineage) · sensitivity {r.get('sensitivity')} · "
                      f"{r.get('chunks')} chunks · {r.get('elapsed_s')} s")
        prov.setObjectName("orient"); prov.setWordWrap(True)
        prov.setToolTip("The exact tool and library this result came from. 'Family models available for this "
                        "lineage' is how many family profiles the installed Dfam partitions hold for this "
                        "organism — a small number means most of its TEs cannot be found, whatever the "
                        "settings. Chunks are the pieces the assembly was split into so progress could be "
                        "shown and the run resumed.")
        lay.addWidget(prov)
        scope = QLabel("What this is: copies of families already present in the installed Dfam library, placed by "
                       "homology. It is not a discovery of new families, and a family absent from the library "
                       "cannot appear here — so absence is not evidence of absence. Percentages are of the whole "
                       "assembly, including unplaced scaffolds.")
        scope.setObjectName("orient"); scope.setWordWrap(True); lay.addWidget(scope)
        if r.get("coverage_warning"):
            cw = QLabel(r["coverage_warning"]); cw.setObjectName("errbanner"); cw.setWordWrap(True); lay.addWidget(cw)
        fams = r.get("families") or []
        t = DataTable([c[0] for c in self._ANNOT_COLS], dict(self._ANNOT_COLS))
        t.set_rows([[f["family"], f["kind"], f["n"], f["bp"], f.get("percent"), f["divergence"]] for f in fams])
        # Row menu matched to what the row IS. A landscape row is a per-family AGGREGATE, not a locus, so
        # it offers no "design primer here" or "copy sequence" — there is no single sequence behind it.
        # A tandem/other row additionally says why it is not counted as a transposable element.
        t.set_row_menu(lambda r, _f=fams: self._annot_row_menu(_f[r]) if r < len(_f) else [])
        t.setMinimumHeight(round(260 * theme_mod.UI_SCALE)); lay.addWidget(t)
        # FigurePanel re-renders on theme/zoom, so it takes a builder (bg -> svg), not a rendered string.
        fig = FigurePanel(lambda bg, _r=r: figures.svg_te_composition(_r, theme=bg), "TEagle_TE_landscape")
        fig.setMinimumHeight(round(190 * theme_mod.UI_SCALE)); lay.addWidget(fig)
        row = QHBoxLayout(); row.addStretch(1)
        row.addWidget(_export_table_btn(t, f"TEagle_TE_landscape_{r.get('accession')}", self))
        # A genome run is sealed like every other result, so the seal has to be reachable from the result
        # — otherwise the run is reproducible in principle and not in practice. The report carries the
        # same hedges the screen shows, so an exported file cannot read as an unqualified claim.
        rep = QPushButton("Export report (.md)"); rep.setProperty("sm", True)
        rep.setToolTip("The summary table plus what was searched, what the numbers mean, and the limits "
                       "stated on this screen — the same hedges, in the file.")
        rep.clicked.connect(lambda: self._save_text(self._annot_report_md(r),
                                                    f"TEagle_TE_landscape_{r.get('accession')}", "md"))
        row.addWidget(rep)
        manb = QPushButton("Export manifest (.json)"); manb.setProperty("sm", True)
        manb.setToolTip("The sealed provenance record: tool and library versions, every parameter, the "
                        "genome checksum, and whether the run covered the whole assembly.")
        manb.setEnabled(bool(r.get("provenance")))
        manb.clicked.connect(lambda: self._save_text(
            json.dumps(r.get("provenance") or {}, indent=2, sort_keys=True),
            f"TEagle_TE_landscape_{r.get('accession')}_manifest", "json"))
        row.addWidget(manb)
        close = QPushButton("Close"); close.setProperty("sm", True); close.clicked.connect(w.accept)
        row.addWidget(close); lay.addLayout(row)
        self._uppercase_buttons(w)
        w.resize(round(1000 * theme_mod.UI_SCALE), round(760 * theme_mod.UI_SCALE))
        w.show(); w.raise_()
        self._annot_window = w

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
        name = org or accession
        box = QMessageBox(self)
        box.setWindowTitle("Delete genome?")
        box.setText(f"Delete the cached {name} genome ({accession}) from this machine?\n\n"
                    f"Its whole-genome off-target scan then needs the genome downloaded again — "
                    "~1 GB from NCBI, several minutes, and ≥8 GB free disk.")
        box.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
        box.setDefaultButton(QMessageBox.No)
        if box.exec() != QMessageBox.Yes:
            return
        self.engine.submit("genome_remove", {"assemblyAccession": accession}, key="genome_remove")

    def _on_genome_remove(self, d):
        # refresh BOTH the dropdown and the manager (if open); the split _on_genome_list self-selects what to update
        self.engine.submit("genome_list", {}, key="genome_list")

    def _render_genome_status(self, msg, level="error"):
        # every caller reports an OUTCOME (failed / cancelled / downloaded / not-yet-downloaded), never an
        # empty panel — so it gets the errbanner vocabulary at the caller's level, not the grey empty style
        _clear_layout(self.genomeBody)
        self.genomeBody.addWidget(_note(msg, level))

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
                acc = (all_assemblies().get(org, {}) or {}).get("assemblyAccession", "")
                self._render_genome_status(f"The {org} genome ({acc}) is not downloaded yet.", "info")
                box = QMessageBox(self)
                box.setWindowTitle("Download genome?")
                box.setText(f"The {org} genome ({acc}) is not on this machine yet.\n\n"
                            "Download it once now? It is kept locally so future scans are fast. "
                            "Mammalian genomes are large (~1 GB download, a few minutes) and need ≥8 GB free disk "
                            "in WSL — extraction peaks near 4 GB (FASTA ~3 GB + 2bit ~0.8 GB) before the temporary "
                            "files are removed.")
                box.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
                box.setDefaultButton(QMessageBox.Yes)
                if org and box.exec() == QMessageBox.Yes:
                    self._prepare_genome(org, then_scan=ps)
                else:
                    self._pending_scan = None
                    self._render_genome_status("Genome scan cancelled — no genome downloaded.", "warn")
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
        self.runPcrBtn.setEnabled(True); self.runPcrBtn.setText("Run loaded pairs")
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
            # the card shows the first pair's seal, but a batch has one manifest per pair — the EXPORT must
            # cover every lane on the gel, not just pair 1. Store the full list for _export_manifest.
            self.state["prov_manifest_all"] = provs

    def _render_pcr(self, lanes, amps):
        _clear_layout(self.pcrBody)
        gel = FigurePanel(lambda bg, L=lanes: figures.svg_gel({"lanes": L}, bg),
                          "TEagle_gel", modes=("dark", "white", "uv", "mono"),
                          hit_regions=figures.gel_regions({"lanes": lanes}), on_menu=self._gel_menu)
        gel.apply_app_theme(self.theme)                   # gel opens in the current app theme (uv/mono stay manual)
        gel.setMinimumHeight(round(420 * theme_mod.UI_SCALE))
        self.pcrBody.addWidget(gel)
        if amps:
            headers = ["Pair", "Source", "Coords", "Len", "Mism F/R", "Call"]
            rows = [[a["pair"], a["source"], f"{a['start']}–{a['end']}", a["length"],
                     f"{a.get('fwd_mm', '—')}/{a.get('rev_mm', '—')}",
                     _amp_call(a)] for a in amps]
            t = DataTable(headers, GLOSS)
            t.set_rows(rows)
            t.set_row_menu(lambda r: self._feat_menu(amps[r]["start"], amps[r]["end"], "+",
                                                     f"amplicon_{amps[r]['pair']}", dna=amps[r].get("seq", "")))
            t.setMaximumHeight(round(200 * theme_mod.UI_SCALE))
            self.pcrBody.addWidget(t)
            arow = QHBoxLayout(); arow.addStretch(1)
            arow.addWidget(_export_table_btn(t, "TEagle_amplicons", self))   # export the notebook table
            if any(a.get("seq") for a in amps):
                def _amps_fasta(_=False, aa=amps, tbl=t):
                    order = [tbl._orig(r) for r in range(tbl.rowCount())]      # follow the table's current (sorted) order
                    seq_amps = [aa[i] for i in order if 0 <= i < len(aa)] or aa
                    fasta = "\n".join(
                        f">amplicon_{a['pair']}_{a['start']}-{a['end']}_{a['length']}bp_"
                        f"{_amp_kind(a)}\n{a.get('seq','')}" for a in seq_amps)
                    widgets.save_fasta(fasta, "TEagle_amplicons", self)
                cp = QPushButton("Export amplicons (FASTA)"); cp.setProperty("sm", True)   # export the sequences
                cp.clicked.connect(_amps_fasta)
                arow.addWidget(cp)
            self.pcrBody.addLayout(arow)
        else:
            self.pcrBody.addWidget(_note("No amplicon predicted for any pair under the criteria.", "warn"))
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
        self._uppercase_buttons(self.pcrBody.parentWidget())   # export buttons are rebuilt with the results

    def _render_genome_scan(self, lanes, amps, summary):
        """Render a whole-genome off-target scan into its own panel: gel + interpretation (verdict,
        per-chromosome spread, size cluster) + a FULL match table, ON-TARGET FIRST. On-target = the product
        at the specimen's own genome locus (when it sits in the scanned assembly); the rest are off-targets,
        or — with no design locus — neutral genomic priming sites."""
        _clear_layout(self.genomeBody)
        gel = FigurePanel(lambda bg, L=lanes: figures.svg_gel({"lanes": L}, bg),
                          "TEagle_genome_gel", modes=("dark", "white", "uv", "mono"),
                          hit_regions=figures.gel_regions({"lanes": lanes}), on_menu=self._gel_menu)
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
        def _call(a):
            return _amp_call(a, has_locus)              # same label the local-PCR table and the gel use
        order = sorted(range(len(amps)), key=lambda i: (amps[i].get("single_primer", False),
                                                        not amps[i].get("on_target", False),
                                                        amps[i]["source"], amps[i]["start"]))
        samps = [amps[i] for i in order]
        if samps:
            # isPcr coordinates are 1-based inclusive and are kept verbatim — the header says so, so Len (end-start+1)
            # does not read as contradicting the span. 'Source' here is the assembly's own sequence name, not a specimen.
            t = DataTable(["Call", "Assembly seq", "Coords (1-based)", "Len", "Strand"], GLOSS)
            t.set_rows([[_call(a), a["source"], f"{a['start']}–{a['end']}", a["length"], a.get("strand", "?")] for a in samps])
            t.set_row_menu(lambda r: [("Copy locus", lambda a=samps[r]: self._copy(f"{a['source']}:{a['start']}-{a['end']}"))])
            t.setMaximumHeight(round(260 * theme_mod.UI_SCALE))
            self.genomeBody.addWidget(t)
            xrow = QHBoxLayout(); xrow.addStretch(1); xrow.addWidget(_export_table_btn(t, "TEagle_offtarget_scan", self))
            self.genomeBody.addLayout(xrow)
        else:
            self.genomeBody.addWidget(_note("No genome-wide product for this pair under the criteria.", "warn"))
        n_on, n_off, n_single = summary.get("n_on", 0), summary.get("n_off", 0), summary.get("n_single", 0)
        head = f"{n_on} on-target + {n_off} off-target site(s)" if has_locus else f"{summary.get('n_pair', 0)} genomic priming site(s)"
        sp = f" · {n_single} single-primer artefact(s)" if n_single else ""
        comig = (" An on-target that shares a band size with an off-target is drawn in the off-target colour — a "
                 "co-migration cannot be resolved on a gel — with the full split kept in the table above." if has_locus else "")
        note = QLabel(f"{head}{sp}, listed above (on-target first).{comig} Every product is a candidate under isPcr's ≥15 bp "
                      "3′-perfect rule — a specificity screen, not wet-lab-validated bands. The verdict is a heuristic "
                      "read of the count/spread; the numbers carry the claim.")
        note.setObjectName("orient"); note.setWordWrap(True); self.genomeBody.addWidget(note)
        self._uppercase_buttons(self.genomeBody.parentWidget())   # this path has no post-render theme refresh

    # =================== WSL family annotation ===================
    def _init_wsl(self):
        self.engine.submit("wsl_status", key="wsl_status")

    def _on_wsl_status(self, w):
        if w.get("error"):
            self.wslStatus.setText(f"<span style='color:{theme_mod.BAD[self.theme]}'>WSL status error: {w['error']}</span>"); return
        if not w.get("wsl2"):
            self.wslStatus.setText("<b>WSL2 not installed</b> — this optional step names the Dfam family. "
                                   "The domain-based superfamily above works without it. Install it in one click with "
                                   "<b>Backend installer</b> below (it runs the elevated <code>wsl --install</code> for you; "
                                   "a Windows restart may be required), or run <code>wsl --install</code> in an Administrator PowerShell.")
            self.spliceStatus.setText("<b>WSL2 not installed</b> — de-novo splice detection needs it (optional).")
            return
        self.engine.submit("genome_list", {}, key="genome_list")   # WSL is up — populate the downloaded-genome dropdown
        # offer the library choice only once the uncurated partitions are actually on the machine —
        # a curated-only backend has nothing extra to search, so the option would promise a search it
        # cannot run. Read from the library famdb reports, not from an assumption about the install.
        _parts = ((w.get("dfam_library") or {}).get("partitions")) or []
        # tracked explicitly, not read back off isVisible(): a widget inside a collapsed card reports
        # not-visible, and reading the choice from that would silently discard what the user selected
        self._uncurated_available = any("uncurated" in str(p).lower() for p in _parts)
        self.wslLibraryRow.setVisible(self._uncurated_available)
        if w.get("ready"):
            # name the partitions actually on the machine — "Dfam curated" was fixed text and read as a
            # statement about the library even where the uncurated partitions were installed alongside it
            _lib = "Dfam curated + uncurated available" if self._uncurated_available else "Dfam curated"
            self.wslStatus.setText(f"<span style='color:{theme_mod.GOOD[self.theme]}'>● ready</span> · RepeatMasker {w.get('repeatmasker')} "
                                   f"· {_lib} · distro {w.get('distro')}")
            self.annotateBtn.setEnabled(True)
        else:
            self.wslStatus.setText(f"WSL2 ok ({w.get('distro')}); annotation stack not installed "
                                   f"(RepeatMasker {w.get('repeatmasker') or 'missing'}, Dfam {'ok' if w.get('dfam') else 'missing'}).")
            self.wslInstallBtn.setVisible(True)
        if w.get("minimap2"):
            self.spliceStatus.setText(f"<span style='color:{theme_mod.GOOD[self.theme]}'>● ready</span> · minimap2 {w.get('minimap2')} "
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
        self.engine.submit("annotate", {"sequence": seq, "species": self._species(), "source": src,
                                        "include_uncurated": bool(self.wslLibrary.currentData())
                                        if getattr(self, "_uncurated_available", False) else False},
                           key="annotate")

    def _on_annotate(self, d):
        self.annotateBtn.setEnabled(True); self.annotateBtn.setText("Run family annotation")
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

    def _repeat_table(self, hits):
        t = DataTable(["#", "Class/family", "Dfam family", "Coords (0-based)", "Str", "Div", "Score"], GLOSS)
        t.set_rows([[i + 1, h["class_family"], h["family"], f"{h['q_start']}–{h['q_end']}",
                     h["strand"], f"{h['divergence']}%", h["score"]] for i, h in enumerate(hits)])
        return t

    def _render_family(self, d):
        hits = d.get("hits", [])
        te = [h for h in hits if h["class_family"] not in self._NONTE | self._UNCLASSED]
        nonte = [h for h in hits if h["class_family"] in self._NONTE]
        unclassed = [h for h in hits if h["class_family"] in self._UNCLASSED]
        if not te:
            cont = QWidget(); cl = QVBoxLayout(cont); cl.setContentsMargins(0, 0, 0, 0); cl.setSpacing(6)
            if nonte:
                kinds = sorted({h["class_family"].replace("_", " ").lower() for h in nonte})
                cl.addWidget(_note(f"RepeatMasker annotated this locus as {' and '.join(kinds)} — "
                                   f"{len(nonte)} region(s). That is a result, not a failure: the sequence "
                                   f"matches a tandem/low-complexity repeat, and no transposable-element "
                                   f"family was found in it.", "info"))
                cl.addWidget(self._repeat_table(nonte))
            if unclassed:
                cl.addWidget(_note(f"{len(unclassed)} repeat region(s) matched a Dfam entry whose class is "
                                   f"unassigned ('Unknown'). A repeat is present; its type is not established.",
                                   "info"))
                cl.addWidget(self._repeat_table(unclassed))
            if not hits:
                # Two independent conditions decide whether a family CAN be named, and the message names
                # whichever is actually missing rather than guessing. Measured on Drosophila copia: with
                # the uncurated partitions installed AND a species it resolves to Copia_LTR and Copia_I at
                # 100% consensus coverage; with the partitions but NO species it returns only
                # low-complexity, because RepeatMasker searches a limited default set without a lineage.
                # How much the curated-only case costs depends entirely on the lineage, so the message
                # below no longer names a figure for one organism — see _CURATED_COVERAGE.
                parts = ((d.get("dfam_library") or {}).get("partitions")) or []
                has_unc = any("uncurated" in str(p).lower() for p in parts)
                searched_unc = bool(d.get("include_uncurated"))
                sp = d.get("species") or ""
                has_sp = bool(sp) and not str(sp).startswith("(all")
                if not has_sp:
                    cl.addWidget(_note("No Dfam family matched — and no organism was set. RepeatMasker "
                                       "searches only a limited default set without a lineage, so an "
                                       "organism is usually what turns a blank result into a named family. "
                                       "Set the organism above and run again.", "warn"))
                elif not searched_unc and not has_unc:
                    cl.addWidget(_note(f"No Dfam family matched for “{sp}” — and only the CURATED Dfam "
                                       "partitions are installed. How much that costs depends entirely on "
                                       "the lineage: " + _curated_coverage() + " A blank result here may "
                                       "therefore be a limit of the installed library rather than a "
                                       "property of the sequence. Install the optional uncurated "
                                       "partitions from BACKEND in the header, then run again with "
                                       "<b>Library</b> set to include them.", "warn"))
                elif not searched_unc:
                    cl.addWidget(_note(f"No Dfam family matched for “{sp}” — but this run searched the "
                                       "CURATED families only. The uncurated partitions are installed on "
                                       "this machine and were not read: RepeatMasker searches curated "
                                       "families unless it is asked for both. Set <b>Library</b> to "
                                       "“Include uncurated families” and run again before concluding "
                                       "anything from this blank result.", "warn"))
                else:
                    cl.addWidget(_note(f"No Dfam family matched this locus for “{sp}”. The curated and "
                                       "uncurated partitions were both searched with a lineage set, so "
                                       "this is a genuine no-match rather than a gap in the installed "
                                       "library. The structural and protein-domain evidence above does "
                                       "not depend on this search.", "warn"))
            frow = QHBoxLayout(); frow.addStretch(1)
            if hits:
                frow.addWidget(_export_table_btn(self._repeat_table(nonte + unclassed), "TEagle_repeats", self))
            cl.addLayout(frow)
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
        head = QLabel(f"<b>Dfam · {te[0]['class_family']}</b> — {' · '.join(sorted({h['family'] for h in te}))} "
                      f"· {_dfam_lbl}{self._src_html('Dfam')} · RepeatMasker {d.get('repeatmasker_version','')}"
                      f"{self._src_html('RepeatMasker')} · species: {d.get('species','')}")
        head.setTextFormat(Qt.RichText); head.setWordWrap(True); head.setOpenExternalLinks(True); _kb_links(head)
        # %div/%del/%ins and the consensus-side coverage are RepeatMasker's own alignment evidence; they
        # were parsed and discarded before. Divergence is RAW (not Kimura-corrected) and the header says so.
        headers = ["#", "Class/family", "Dfam family", "Coords (0-based)", "Str",
                   "Div (raw)", "Del", "Ins", "Consensus", "Cov", "Blocks", "Score"]
        t = DataTable(headers, GLOSS)

        def _cons(h):
            if h.get("cons_start") is None or h.get("cons_length") is None:
                return "—"
            return f"{h['cons_start']}–{h['cons_end']} / {h['cons_length']}"

        t.set_rows([[i + 1, h["class_family"], h["family"], f"{h['q_start']}–{h['q_end']}", h["strand"],
                     f"{h['divergence']}%",
                     ("—" if h.get("pct_del") is None else f"{h['pct_del']}%"),
                     ("—" if h.get("pct_ins") is None else f"{h['pct_ins']}%"),
                     _cons(h),
                     ("—" if h.get("cons_coverage_pct") is None else f"{h['cons_coverage_pct']}%"),
                     h.get("n_fragments", 1), h["score"]] for i, h in enumerate(te)])
        t.set_row_menu(lambda r: self._feat_menu(te[r]["q_start"], te[r]["q_end"], te[r]["strand"],
                                                 te[r]["family"], src_seq=self.state.get("family_seq")))
        cont = QWidget(); cl = QVBoxLayout(cont); cl.setContentsMargins(0, 0, 0, 0); cl.setSpacing(6)
        cl.addWidget(head); cl.addWidget(t)
        merged = [h for h in te if h.get("n_fragments", 1) > 1]
        if merged:
            cl.addWidget(_note("; ".join(h["fragment_note"] for h in merged) +
                               ". Merged so an interrupted element is one row, with the outermost query span "
                               "and the summed consensus coverage.", "info"))
        cl.addWidget(_qq_note())
        if nonte or unclassed:                     # reported, never silently dropped
            extra = [f"{len(nonte)} tandem / low-complexity region(s)"] if nonte else []
            extra += [f"{len(unclassed)} repeat(s) of unassigned class"] if unclassed else []
            cl.addWidget(_note("Also annotated in this locus: " + " and ".join(extra) +
                               ", listed apart from the TE families above.", "info"))
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
        self.spliceBtn.setEnabled(True); self.spliceBtn.setText("Detect exons / introns")
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
        head = QLabel(f"<b>{d['counts']['exons']} exon(s) · {d['counts']['introns']} intron(s)</b> — de novo · "
                      f"{d.get('canonical_introns',0)}/{d['counts']['introns']} canonical splice site(s) · strand {d.get('strand','')}"
                      f"{self._src_html('minimap2')}")
        head.setTextFormat(Qt.RichText); head.setWordWrap(True); head.setOpenExternalLinks(True)
        _kb_links(head); cl.addWidget(head)
        # independent cross-check: does this de-novo (external-transcript) alignment agree with the record's
        # own annotation? Only meaningful when the genomic IS the fetched record (source accession present).
        src, feats = self.state.get("source"), self.state.get("features")
        if src and src.get("accession") and isinstance(feats, dict) and feats.get("introns"):
            cc = cross_check_models(complete_gene_model(feats).get("introns", []), d.get("introns", []))
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
        model = figures.gv_tracks_from_gene({"exons": d.get("exons", []), "introns": d.get("introns", []), "cds": []}, length)
        if model["tracks"]:
            gv = GenomePanel(svg_genome, "TEagle_splice"); gv.apply_app_theme(self.theme); gv.set_model(model)
            gv.set_feature_menu(self._region_menu); gv.setMinimumHeight(round(220 * theme_mod.UI_SCALE))
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
        nr = "<br>".join("· " + n for n in m.get("notRun", []))
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
        lab.setTextInteractionFlags(Qt.TextSelectableByMouse)   # the hashes/versions must be selectable to copy
        self.card_prov.bodylay.addWidget(lab)
        # the seal is the reproducibility record — it must be able to leave the app, not only be read on screen.
        self.state["prov_manifest"] = m
        self.state["prov_manifest_all"] = [m]              # single by default; a PCR batch overrides with all pairs
        erow = QHBoxLayout(); erow.addStretch(1)
        eb = QPushButton("Export manifest (.json)"); eb.setProperty("sm", True)
        eb.setToolTip("Write the full run-provenance manifest (input hash, tool + database versions, "
                      "parameters, checksums) as a JSON file so the result stays reproducible outside TEagle.")
        eb.clicked.connect(self._export_manifest)
        erow.addWidget(eb); self.card_prov.bodylay.addLayout(erow)
        self._uppercase_buttons()

    def _export_manifest(self):
        import json
        allm = self.state.get("prov_manifest_all") or ([self.state["prov_manifest"]]
                                                        if self.state.get("prov_manifest") else [])
        if not allm:
            return self._banner("No run provenance to export — run an analysis first.")
        rid = (str(allm[0].get("input", {}).get("id") or "run")).split()[0]
        path, _ = QFileDialog.getSaveFileName(self, "Export provenance manifest",
                                              f"TEagle_{rid}_manifest.json", "JSON (*.json)")
        if not path:
            return
        if not path.lower().endswith(".json"):
            path += ".json"
        # a batch (multi-pair in-silico PCR) seals one manifest PER pair — export them all so the file
        # covers every lane shown, not just the first; a single run writes its one manifest as before.
        payload = {"runs": allm} if len(allm) > 1 else allm[0]
        with open(path, "w", encoding="utf-8", newline="\n") as fh:
            json.dump(payload, fh, indent=2, ensure_ascii=False, sort_keys=True)
        self._banner(f"Provenance manifest written to {os.path.basename(path)}.", level="success")

    def _set_body(self, layout, widget):
        _clear_layout(layout)                                 # recursive: also removes any addLayout'd sub-layouts
        layout.addWidget(widget)

    def closeEvent(self, e):
        # Stop EVERY polling timer, not a named one. This used to stop only the genome-download poll, so
        # each timer added later (the annotation progress poll, for one) could still fire while the window
        # was being torn down. Walking the children keeps that correct without anyone remembering to.
        for t in self.findChildren(QTimer):
            try:
                t.stop()
            except RuntimeError:                              # already destroyed by Qt — nothing to stop
                pass
        # Off-thread engine jobs are deliberately NOT waited on: a whole-genome annotation runs for hours
        # inside WSL and resumes on the next launch, so blocking the close would trade a slow exit for no
        # benefit. _Job._emit already drops results whose receiver has gone.
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
    from teagle_core import oligoqc                          # secondary-structure QC; primer3 required, ViennaRNA optional
    if not oligoqc.available().get("primer3"):               # the PRIMARY engine must be in the bundle
        problems.append(f"primer3 secondary-structure engine missing from bundle "
                        f"({oligoqc.available().get('primer3_error')})")
    # ViennaRNA is deliberately NOT bundled — its licence forbids redistribution inside an AGPL work, so it
    # is a user-installed optional second engine. Its absence is the expected shipped state, never a build
    # failure; the UI already reports which engines ran and the manifest records it.
    if domains.PYHMMER_VERSION == "unavailable":
        problems.append("pyhmmer failed to load")
    if not domains.HMM_SHA256:
        problems.append("bundled Pfam HMM profiles not found")
    else:                                                 # the profile set must LOAD and cover the full GAG-POL-ENV panel
        try:
            n_hmm = len(domains._hmms())
            codes = {v[0] for v in domains.DOMAIN_INFO.values()}
            # every class the 30-model panel claims must be represented, incl. the LINE ORF1p/EN, the DIRS
            # tyrosine recombinase (YR) and the Helitron (HEL) modules added this release
            miss = {"GAG", "PR", "RT", "RNaseH", "INT", "ENV", "CHR", "TPase", "ORF1", "EN", "YR", "HEL"} - codes
            if miss:
                problems.append(f"domain profile set missing expected codes: {sorted(miss)}")
            if n_hmm < len(domains.DOMAIN_INFO):          # the bundled .hmm must load the full declared panel, not a truncated subset
                problems.append(f"domain profile set loaded {n_hmm} profiles, expected {len(domains.DOMAIN_INFO)} "
                                "(bundled .hmm may be truncated)")
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
            # go through the lazy accessor, which is the path a real export takes — this check exists to
            # prove openpyxl actually LOADS and writes inside the frozen bundle, not merely that it is present
            wb = widgets._xl()["Workbook"](); wb.active.append(["h", 1]); wb.save(io.BytesIO())
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
    _vr = oligoqc.VIENNARNA_VERSION if oligoqc.available().get("viennarna") else "not installed (optional)"
    print(f"TEAGLE SELFTEST OK · primer3 {primers.PRIMER3_VERSION} · ViennaRNA {_vr} "
          f"· pyhmmer {domains.PYHMMER_VERSION} · HMM {domains.HMM_SHA256[:12]} ({len(domains._hmms())} profiles) "
          f"· QtSvg ok · install dialog ok")
    return 0


UI_SCALES = [0.75, 0.85, 1.0, 1.1, 1.25, 1.5]     # user-selectable global UI scale (persisted, applied at startup)


def _apply_saved_ui_scale():
    """Seed the live UI scale from QSettings into theme_mod.UI_SCALE BEFORE MainWindow is built, so first
    paint is correct. Scale is driven at runtime (theme_mod.UI_SCALE + _apply_theme) — no QT_SCALE_FACTOR,
    no restart. An explicit TEAGLE_UI_SCALE env override (screenshot harness / dev) always wins."""
    if os.environ.get("TEAGLE_UI_SCALE"):
        return
    try:
        f = float(QSettings("TEagle", "TEagle").value("ui_scale", 1.0))
        if 0.5 <= f <= 3.0:
            theme_mod.UI_SCALE = f
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
