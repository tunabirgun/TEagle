"""TEagle native desktop app (PySide6). QMainWindow shell with a specimen rail and a scrollable
column of collapsible result cards. All science runs in-process through the shared engine, off the
GUI thread via engine_worker.Engine. This module wires the analyze workflow, primer/PCR, WSL family
annotation, splice detection, provenance and exports."""
from __future__ import annotations
import gzip, os, sys

from PySide6.QtCore import Qt, QTimer, QByteArray
from PySide6.QtGui import QGuiApplication, QFont, QPixmap, QPainter, QIcon
from PySide6.QtSvg import QSvgRenderer
from PySide6.QtWidgets import (QApplication, QMainWindow, QWidget, QFrame, QVBoxLayout, QHBoxLayout,
                               QGridLayout, QLabel, QLineEdit, QTextEdit, QPushButton, QComboBox,
                               QScrollArea, QSplitter, QFileDialog, QSizePolicy)

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.join(os.path.dirname(_HERE), "backend"))

from engine_worker import Engine
import fonts
import figures
from figures import svg_genome
import widgets
from widgets import FigurePanel, GenomePanel, DataTable
from sample import make_sample
import theme as theme_mod
from teagle_core import appdirs

APP_VERSION = "2.0.0"
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

# wordmark: clean Cascadia Code Bold 'TEAGLE' frozen to static paths, TE/AGLE recolored per theme
_WORD_SVG = None
def _word_pixmap(te: str, agle: str, height: int, dpr: float = 1.0) -> QPixmap:
    global _WORD_SVG
    if _WORD_SVG is None:
        _WORD_SVG = _load_asset("teagle-wordmark.svg")
    return _svg_pixmap(_WORD_SVG.replace("{TE}", te).replace("{AGLE}", agle), height, dpr)

def _app_icon() -> QIcon:                                 # window/taskbar icon from the mark, mid-teal on transparent
    icon = QIcon()
    for s in (16, 20, 24, 32, 40, 48, 64, 128, 256):
        big = _mark_pixmap(ICON_TEAL, s * 4)              # supersample x4 then smooth-scale -> clean small frames
        icon.addPixmap(big.scaled(s, s, Qt.KeepAspectRatio, Qt.SmoothTransformation))
    return icon

STRUCT_COLS = ["Feature", "Coords (0-based)", "Len", "Metric", "Method"]
ORF_COLS = ["Strand", "Frame", "Start", "End", "aa"]
DOMAIN_COLS = ["Domain", "Label", "Pfam", "aa", "nt", "Score", "E-value"]
# plain-language glossary for every table header (hover to learn the abbreviation) — mirrors web GLOSSARY
GLOSS = {
    "Feature": "Structural hallmark found in the sequence — e.g. LTR, TIR, target-site duplication, poly-A tail.",
    "Coords (0-based)": "Location in the sequence, 0-based half-open [start, end).",
    "Coords": "Location of this amplicon in the searched sequence.",
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
    "ID": "Identifier of this designed primer pair.",
    "Forward (5'→3')": "Forward primer sequence, written 5′→3′.",
    "Reverse (5'→3')": "Reverse primer sequence, written 5′→3′.",
    "Product": "Predicted amplicon (product) size in base pairs.",
    "Tm F/R": "Melting temperature (°C) of the forward / reverse primer — matched so both anneal together.",
    "GC% F/R": "Percent G+C content of the forward / reverse primer; ~40–60% is typical.",
    "Penalty": "Primer3's overall penalty for the pair — lower is better; it rises as primers depart from target Tm, size and GC.",
    "Pair": "Which designed primer pair produced this amplicon.",
    "Source": "The sequence that was searched — your specimen or a custom background.",
    "Mism F/R": "Mismatches in the forward / reverse primer binding site (the 3′ end is kept exact).",
    "Call": "On-target = amplicon at the intended locus; off-target = amplified elsewhere.",
    "Class/family": "TE class and superfamily (Wicker 2007 scheme), e.g. LTR/Copia.",
    "Dfam family": "The specific named family in the Dfam library, e.g. Copia_I or L1HS.",
    "Str": "Strand of the match — + forward, − reverse complement.",
    "Div": "Divergence — % difference between your sequence and the Dfam family consensus (lower = closer).",
    "Intron span (0-based)": "Intron location in the loaded sequence, 0-based half-open [start, end).",
    "Splice site": "The two bases at each intron boundary (donor…acceptor); canonical introns are GT…AG (or GC–AG / AT–AC).",
    "Canonical": "Whether the intron's donor…acceptor matches a canonical eukaryotic splice motif.",
    "#": "Row number.",
}

# clickable source citations (verified DOIs — mirror backend refs.py and the web REFLINKS)
REFLINKS = {
    "Wicker2007":   {"url": "https://doi.org/10.1038/nrg2165", "cite": "Wicker T, et al. (2007) A unified classification system for eukaryotic transposable elements. Nat Rev Genet 8:973-982."},
    "Pfam":         {"url": "https://www.ebi.ac.uk/interpro/", "cite": "Mistry J, et al. (2021) Pfam: the protein families database in 2021. Nucleic Acids Res 49:D412-D419."},
    "HMMER":        {"url": "https://doi.org/10.1371/journal.pcbi.1002195", "cite": "Eddy SR (2011) Accelerated Profile HMM Searches. PLoS Comput Biol 7:e1002195."},
    "Dfam":         {"url": "https://doi.org/10.1186/s13100-020-00230-y", "cite": "Storer J, et al. (2021) The Dfam community resource of transposable element families. Mob DNA 12:2."},
    "RepeatMasker": {"url": "https://www.repeatmasker.org/", "cite": "Smit AFA, Hubley R, Green P. RepeatMasker Open-4.0."},
    "NCBI":         {"url": "https://www.ncbi.nlm.nih.gov/nuccore/", "cite": "NCBI Entrez / E-utilities (Sayers E, NCBI)."},
    "Primer3":      {"url": "https://doi.org/10.1093/nar/gks596", "cite": "Untergasser A, et al. (2012) Primer3 — new capabilities and interfaces. Nucleic Acids Res 40:e115."},
    "minimap2":     {"url": "https://doi.org/10.1093/bioinformatics/bty191", "cite": "Li H (2018) Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics 34:3094-3100."},
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
        self._number, self._title = number, title
        self.set_collapsed(collapsed)

    def _toggle(self):
        self.set_collapsed(self.body.isVisible())

    def set_collapsed(self, collapsed: bool):
        self.body.setVisible(not collapsed)
        arrow = "▸" if collapsed else "▾"
        title = self._title.replace("&", "&&")           # QPushButton eats a lone '&' as a mnemonic
        self.hdr.setText(f"{arrow} {self._number}  {title}")

    def expand(self):
        self.set_collapsed(False)

    def clear_body(self):
        while self.bodylay.count():
            it = self.bodylay.takeAt(0)
            w = it.widget()
            if w:
                w.setParent(None)


def _hline():
    f = QFrame(); f.setObjectName("hline"); f.setFixedHeight(1); return f


def _empty(text):
    l = QLabel(text); l.setObjectName("empty"); l.setWordWrap(True); return l


def _sl(text):
    """A section label — mono, uppercase, tracked (the web UI's `.lbl`)."""
    l = QLabel(text.upper()); l.setObjectName("sectionlabel"); l.setWordWrap(True); return l


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("TEagle")
        self.resize(round(1240 * theme_mod.UI_SCALE), round(860 * theme_mod.UI_SCALE))
        self.theme = "dark"
        self.state = {"seq": "", "source": None, "last_rec": None}

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
        split.addWidget(self._build_rail())
        split.addWidget(self._build_results())
        split.setStretchFactor(0, 0)
        split.setStretchFactor(1, 1)
        split.setSizes([340, 900])
        outer.addWidget(split, 1)

        self._apply_theme()
        QTimer.singleShot(0, self._startup)

    # ---------- header ----------
    def _build_header(self):
        wrap = QWidget()
        col = QVBoxLayout(wrap); col.setContentsMargins(6, 0, 2, 0); col.setSpacing(0)
        h = QHBoxLayout(); h.setContentsMargins(0, 2, 0, 8); h.setSpacing(10)
        self.mark = QLabel()                                  # eagle brand mark; pixmap set per-theme in _apply_theme
        self.mark.setObjectName("mark")
        self.mark.setToolTip("TEagle")
        h.addWidget(self.mark)
        self.word = QLabel()                                  # Cascadia Code wordmark; pixmap set per-theme in _apply_theme
        self.word.setObjectName("word")
        h.addWidget(self.word)
        self.ver = QLabel("v" + APP_VERSION); self.ver.setObjectName("ver")
        h.addWidget(self.ver)
        tag = QLabel("TRANSPOSABLE ELEMENTS"); tag.setObjectName("tagline")
        tf = tag.font(); tf.setLetterSpacing(QFont.AbsoluteSpacing, 1.5); tag.setFont(tf)
        h.addWidget(tag)
        h.addStretch(1)
        chip = QFrame(); chip.setObjectName("statuschip")
        cl = QHBoxLayout(chip); cl.setContentsMargins(10, 5, 11, 5); cl.setSpacing(8)
        self.led = QLabel(); self.led.setObjectName("led"); self.led.setFixedSize(8, 8)
        self.statusTxt = QLabel("connecting…"); self.statusTxt.setObjectName("statusTxt")
        cl.addWidget(self.led); cl.addWidget(self.statusTxt)
        h.addWidget(chip)
        tb = QPushButton("◐ THEME"); tb.setProperty("sm", True); tb.clicked.connect(self._toggle_theme)
        h.addWidget(tb)
        col.addLayout(h)
        self.headrule = QFrame(); self.headrule.setObjectName("headrule"); self.headrule.setFixedHeight(2)
        col.addWidget(self.headrule)
        return wrap

    # ---------- rail ----------
    def _build_rail(self):
        rail = QFrame(); rail.setObjectName("rail")
        rail.setMinimumWidth(300); rail.setMaximumWidth(430)
        lay = QVBoxLayout(rail); lay.setContentsMargins(12, 12, 12, 12); lay.setSpacing(8)
        lay.addWidget(self._sec("01", "Specimen"))
        accrow = QHBoxLayout()
        self.acc = QLineEdit(); self.acc.setPlaceholderText("accession — e.g. M11240, NC_003075.7")
        accrow.addWidget(self.acc)
        fb = QPushButton("↓ Fetch"); fb.setProperty("sm", True); fb.clicked.connect(self._fetch)
        accrow.addWidget(fb)
        lay.addLayout(accrow)
        self.accMeta = QLabel(""); self.accMeta.setObjectName("cardmeta"); self.accMeta.setWordWrap(True)
        self.accMeta.setTextFormat(Qt.RichText); self.accMeta.setOpenExternalLinks(True)
        lay.addWidget(self.accMeta)
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
        self.card_prov = CollapsibleCard("06", "Run provenance", "versions + checksums for reproducibility")
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
        self.wslSpecies = QLineEdit(); self.wslSpecies.setPlaceholderText("species (e.g. drosophila melanogaster) — from accession if fetched")
        row.addWidget(self.wslSpecies)
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
    def _seq_changed(self):
        txt = self.seq.toPlainText()
        n = len(txt.replace("\n", "").replace("\r", ""))
        # rough nt count excluding header lines
        body = "".join(l for l in txt.splitlines() if not l.startswith(">"))
        self.charCount.setText(f"{len(body)} nt")

    def _load_sample(self):
        self.seq.setPlainText(make_sample())
        self.state["source"] = None
        self.accMeta.setText("")

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
        self.seq.setPlainText(data)
        self.state["source"] = None
        self.accMeta.setText(f"loaded {os.path.basename(path)}")

    def _fetch(self):
        acc = self.acc.text().strip()
        if not acc:
            return self._banner("enter an accession first")
        self.accMeta.setText("fetching…")
        self.engine.submit("fetch", {"accession": acc}, key="fetch")

    def _run_analysis(self):
        seq = self.seq.toPlainText().strip()
        if not seq:
            return self._banner("paste, upload, or fetch a sequence first")
        self.state["seq"] = seq
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
            self._on_pcr(key, res)

    def _reset_buttons(self, key):
        if key == "analyze":
            self.runBtn.setEnabled(True); self.runBtn.setText("▶ Run analysis")
        elif key == "primers":
            self.designBtn.setEnabled(True); self.designBtn.setText("⋯ Design primers")
        elif key.startswith("pcr#"):
            self.runPcrBtn.setEnabled(True); self.runPcrBtn.setText("▶ Run loaded pairs")
        elif key == "annotate":
            self.annotateBtn.setEnabled(True); self.annotateBtn.setText("▶ Run family annotation")
        elif key == "splice":
            self.spliceBtn.setEnabled(True); self.spliceBtn.setText("▶ Detect exons / introns")
        elif key == "fetch":
            self.accMeta.setText("")

    def _on_user_error(self, key, msg):
        self._reset_buttons(key)
        self._banner(msg)

    def _on_failed(self, key, msg, trace):
        self._reset_buttons(key)
        sys.stderr.write(trace + "\n")
        self._banner(f"unexpected error: {msg}")

    def _banner(self, msg):
        self.errbanner.setText("⚠ " + msg)
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

    def _on_fetch(self, res):
        if not res.get("ok"):
            self.accMeta.setText("")
            return self._banner(res.get("error", "fetch failed"))
        self._clear_banner()
        seqtext = res.get("fasta") or res.get("sequence") or ""
        if seqtext:
            self.seq.setPlainText(seqtext)
        org = res.get("organism", "")
        length = res.get("length") or res.get("seq_length") or ""
        cached = " · cached (local)" if res.get("fromCache") else ""
        acc = res.get("accession", "")
        ncbi = self._src_html("NCBI", "https://www.ncbi.nlm.nih.gov/nuccore/" + acc) if acc else ""
        self.accMeta.setText(f"{acc} · {org} · {length} bp{cached}{ncbi}<br>{res.get('title','')}")
        self.state["source"] = {k: res.get(k) for k in ("accession", "organism", "title", "length", "moltype") if res.get(k) is not None}
        self.state["features"] = res.get("features")
        # auto-fill species for WSL family annotation if present
        if hasattr(self, "wslSpecies") and org:
            self.wslSpecies.setText(org.lower())

    def _on_analyze(self, res):
        self.runBtn.setEnabled(True); self.runBtn.setText("▶ Run analysis")
        self._clear_banner()
        if res.get("warning"):
            self._banner(res["warning"])
        recs = res.get("records", [])
        if not recs:
            return
        rec = recs[0]
        self.state["last_rec"] = rec
        self.state["analyzed_clean"] = self._clean_seq(self.state.get("seq", ""))
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
        card.bodylay.addWidget(banner)

        # genome viewer
        model = figures.gv_tracks_from_rec(rec)
        if model["tracks"]:
            gv = GenomePanel(svg_genome, "TEagle_genome")
            gv.set_model(model)
            gv.set_feature_menu(self._region_menu)
            gv.setMinimumHeight(260)
            card.bodylay.addWidget(gv)

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
                         d.get("score"), f"{d.get('evalue'):.1e}" if d.get("evalue") is not None else ""]
                        for d in doms])
            t.set_row_menu(lambda r: self._feat_menu(doms[r]["nt"][0], doms[r]["nt"][1], doms[r].get("strand", "+"),
                                                     doms[r]["domain"], protein=doms[r].get("protein")))
            t.setMaximumHeight(180)
            card.bodylay.addWidget(t)

        # gene model (exon/intron/CDS) — only when a fetched accession carries feature annotation
        gm = self.state.get("features")
        if isinstance(gm, dict) and (gm.get("exons") or gm.get("cds")):
            card.bodylay.addWidget(_sl("Gene model (NCBI feature table)"))
            length = rec.get("composition", {}).get("length") or 1
            gmodel = figures.gv_tracks_from_gene(gm, length)
            if gmodel["tracks"]:
                gvg = GenomePanel(svg_genome, "TEagle_genemodel"); gvg.set_model(gmodel)
                gvg.set_feature_menu(self._region_menu); gvg.setMinimumHeight(200)
                card.bodylay.addWidget(gvg)

        for note in rec.get("notes", []):
            n = QLabel("• " + note); n.setObjectName("orient"); n.setWordWrap(True)
            card.bodylay.addWidget(n)

    def _struct_row(self, e):
        sp = e.get("element_span") or e.get("five_prime") or e.get("pos") or e.get("upstream") or [None, None]
        arm = e.get("ltr_len") or e.get("tir_len") or e.get("length") or ""
        metric = (f"{e['identity']}%" if e.get("identity") is not None else e.get("motif", ""))
        coords = f"{sp[0]}–{sp[1]}" if sp[0] is not None else ""
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

    def _slice(self, s, e):
        return self._clean_seq(self.state.get("seq", ""))[s:e]

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
        self.statusTxt.setText("copied to clipboard")
        QTimer.singleShot(1500, lambda: self.statusTxt.setText(self.statusTxt.text()))

    def _feat_menu(self, start, end, strand, label, protein=None):
        dna = self._revcomp(self._slice(start, end)) if strand == "-" else self._slice(start, end)
        items = [(f"⧉ Copy FASTA", lambda: self._copy(f">{label}_{start}-{end}{'_rev' if strand=='-' else ''}\n{dna}")),
                 (f"⧉ Copy DNA", lambda: self._copy(dna)),
                 (f"⧉ Copy coords ({start}–{end} {strand})", lambda: self._copy(f"{start}-{end} {strand}"))]
        if protein:
            items.append(("⧉ Copy protein", lambda: self._copy(protein)))
        items.append(("⌖ Design primer here", lambda: self._design_for_domain(start, end, label)))
        return items

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
        onoff = "on" if a.get("on_target") else "off"
        label = f"amplicon_{pair}_{start}-{end}_{a.get('length','')}bp_{onoff}target"
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
        if self._stale_block():
            return
        self.designBtn.setEnabled(False); self.designBtn.setText("◴ designing…")
        self.engine.submit("primers", {"sequence": self.state["seq"], "params": self._read_primer_params()}, key="primers")

    def _design_for_domain(self, start, end, label):
        self._clear_banner()
        if self._stale_block():
            return
        inc = [start, max(60, end - start)]
        self.card_primer.expand()
        self._pending_domain = label
        self.engine.submit("primers", {"sequence": self.state["seq"], "params": self._read_primer_params(),
                                        "included": inc}, key="primers")

    def _on_primers(self, d):
        self.designBtn.setEnabled(True); self.designBtn.setText("⋯ Design primers")
        self.card_primer.expand()
        self._render_primers(d)
        self._uppercase_buttons()
        if d.get("provenance"):
            self._render_provenance(d["provenance"])

    def _render_primers(self, d):
        while self.primBody.count():
            w = self.primBody.takeAt(0).widget()
            if w:
                w.setParent(None)
        cands = d.get("candidates", [])
        self.state["candidates"] = cands
        self.pcrStageAll.setEnabled(bool(cands))
        if not cands:
            self.primBody.addWidget(_empty("No primer pair met the criteria — try the Permissive preset or widen the product range."))
            return
        srow = QHBoxLayout(); srow.addWidget(_sl(f"Primer pairs — {len(cands)}")); srow.addStretch(1)
        srow.addWidget(self._src_link("Primer3")); self.primBody.addLayout(srow)
        headers = ["ID", "Forward (5'→3')", "Reverse (5'→3')", "Product", "Tm F/R", "GC% F/R", "Penalty"]
        t = DataTable(headers, GLOSS)
        t.set_rows([[c["id"], c["left_seq"], c["right_seq"], c["product_size"],
                     f"{c['left_tm']}/{c['right_tm']}", f"{c['left_gc']}/{c['right_gc']}", c["penalty"]] for c in cands])
        t.set_row_menu(lambda r: [("→ send to in-silico PCR", lambda rr=r: self._add_pcr_pair(cands[rr])),
                                  ("⧉ Copy pair FASTA", lambda rr=r: self._copy(
                                      f">{cands[rr]['id']}_F\n{cands[rr]['left_seq']}\n>{cands[rr]['id']}_R\n{cands[rr]['right_seq']}"))])
        t.setMaximumHeight(200)
        self.primBody.addWidget(t)
        hint = QLabel("Right-click a pair → “send to in-silico PCR”, or stage all in panel 05.")
        hint.setObjectName("orient"); self.primBody.addWidget(hint)

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
        if self._stale_block():
            return
        pairs = self.state.get("pcrPairs", [])
        if not pairs:
            return self._banner("Load at least one primer pair first.")
        p = {"max_mm": self._numfield("pcrMM", 2), "tp": self._numfield("pcrTP", 5),
             "prod_min": self._numfield("pcrPmin", 70), "prod_max": self._numfield("pcrPmax", 1000)}
        bg = self.pcrBg.toPlainText()
        self._pcr_run = {"expected": len(pairs), "results": [None] * len(pairs), "single": len(pairs) == 1}
        self.runPcrBtn.setEnabled(False); self.runPcrBtn.setText("◴ running…")
        for i, c in enumerate(pairs):
            body = {"sequence": self.state["seq"], "background": bg, "fwd": c["left_seq"], "rev": c["right_seq"],
                    "target_span": [c["left_pos"][0], c["right_pos"][1]], "params": p}
            self.engine.submit("pcr", body, key=f"pcr#{i}")

    def _on_pcr(self, key, d):
        i = int(key.split("#")[1])
        run = getattr(self, "_pcr_run", None)
        if not run or i >= len(run["results"]):
            return
        run["results"][i] = d
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
        self._render_pcr(lanes, amps, results if run["single"] else None)
        self._uppercase_buttons()
        if provs:
            self._render_provenance(provs[0])

    def _render_pcr(self, lanes, amps, single_results):
        while self.pcrBody.count():
            w = self.pcrBody.takeAt(0).widget()
            if w:
                w.setParent(None)
        gel = FigurePanel(lambda bg, L=lanes: figures.svg_gel({"lanes": L}, bg),
                          "TEagle_gel", modes=("dark", "white", "uv", "mono"),
                          hit_regions=figures.gel_regions({"lanes": lanes}), on_menu=self._gel_menu)
        gel.setMinimumHeight(420)
        self.pcrBody.addWidget(gel)
        if amps:
            headers = ["Pair", "Source", "Coords", "Len", "Mism F/R", "Call"]
            t = DataTable(headers, GLOSS)
            t.set_rows([[a["pair"], a["source"], f"{a['start']}–{a['end']}", a["length"],
                         f"{a['fwd_mm']}/{a['rev_mm']}", "on-target" if a.get("on_target") else "off-target"] for a in amps])
            t.set_row_menu(lambda r: self._feat_menu(amps[r]["start"], amps[r]["end"], "+", f"amplicon_{amps[r]['pair']}"))
            t.setMaximumHeight(200)
            self.pcrBody.addWidget(t)
            cp = QPushButton("⧉ amplicons → FASTA"); cp.setProperty("sm", True)
            cp.clicked.connect(lambda: self._copy("\n".join(
                f">amplicon_{a['pair']}_{a['start']}-{a['end']}_{a['length']}bp_{'on' if a.get('on_target') else 'off'}target\n{a.get('seq','')}" for a in amps)))
            self.pcrBody.addWidget(cp)
        else:
            self.pcrBody.addWidget(_empty("No amplicon predicted for any pair under the criteria."))
        onN = sum(1 for a in amps if a.get("on_target"))
        note = QLabel(f"{len(lanes)} lane(s) · {onN} on-target band(s) · ladder lane “L”. Not a claim of experimental specificity.")
        note.setObjectName("orient"); self.pcrBody.addWidget(note)

    # =================== WSL family annotation ===================
    def _init_wsl(self):
        self.engine.submit("wsl_status", key="wsl_status")

    def _on_wsl_status(self, w):
        if w.get("error"):
            self.wslStatus.setText(f"<span style='color:#E06A5A'>WSL status error: {w['error']}</span>"); return
        if not w.get("wsl2"):
            self.wslStatus.setText("<b>WSL2 not installed</b> — this optional step names the Dfam family. "
                                   "The domain-based superfamily above works without it. Install: open PowerShell as "
                                   "Administrator and run <code>wsl --install</code>, then restart Windows.")
            self.spliceStatus.setText("<b>WSL2 not installed</b> — de-novo splice detection needs it (optional).")
            return
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

    def _annotate(self):
        if self.wslSource.currentIndex() == 1:
            seq = self.wslPaste.toPlainText().strip()
            src = None
        else:
            seq = self.state.get("seq") or self.seq.toPlainText().strip()
            src = self.state.get("source")
        if not seq:
            return self._banner("no sequence to annotate — load a specimen or paste one")
        self.annotateBtn.setEnabled(False); self.annotateBtn.setText("◴ annotating…")
        self.engine.submit("annotate", {"sequence": seq, "species": self.wslSpecies.text().strip() or None,
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
        t.set_row_menu(lambda r: self._feat_menu(te[r]["q_start"], te[r]["q_end"], te[r]["strand"], te[r]["family"]))
        cont = QWidget(); cl = QVBoxLayout(cont); cl.setContentsMargins(0, 0, 0, 0); cl.setSpacing(6)
        cl.addWidget(head); cl.addWidget(t)
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
        genomic = self.state.get("seq") or self.seq.toPlainText().strip()
        tx = self.spliceTx.toPlainText().strip()
        if not genomic.strip():
            return self._banner("Load a genomic sequence first (fetch, upload, or paste, then Run analysis).")
        if not tx:
            return self._banner("Paste a transcript / cDNA / mRNA to align.")
        self.spliceBtn.setEnabled(False); self.spliceBtn.setText("◴ aligning…")
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
        length = (self.state.get("last_rec") or {}).get("composition", {}).get("length") or len(self._clean_seq(self.state.get("seq", ""))) or 1
        model = figures.gv_tracks_from_gene({"exons": d.get("exons", []), "introns": d.get("introns", []), "cds": []}, length)
        if model["tracks"]:
            gv = GenomePanel(svg_genome, "TEagle_splice"); gv.set_model(model)
            gv.set_feature_menu(self._region_menu); gv.setMinimumHeight(220)
            cl.addWidget(gv)
        introns = d.get("introns", [])
        if introns:
            headers = ["#", "Intron span (0-based)", "Len", "Splice site", "Canonical"]
            t = DataTable(headers, GLOSS)
            t.set_rows([[k + 1, f"{i['start']}–{i['end']}", i["end"] - i["start"], f"{i['donor']}…{i['acceptor']}",
                         "canonical" if i.get("canonical") else "non-canonical"] for k, i in enumerate(introns)])
            t.set_row_menu(lambda r: self._feat_menu(introns[r]["start"], introns[r]["end"], d.get("strand", "+"), f"intron_{r+1}"))
            cl.addWidget(t)
        else:
            cl.addWidget(_empty("single exon — no introns detected"))
        self._set_body(self.spliceBody, cont)

    # =================== provenance ===================
    def _render_provenance(self, m):
        self.card_prov.expand()
        self.card_prov.clear_body()
        inp = m.get("input", {})
        sw = "<br>".join(f"{s['name']} · {s['version']}" for s in m.get("software", []))
        pr = "<br>".join(f"{k} · {'—' if v is None else v}" for k, v in (m.get("parameters") or {}).items()) or "defaults"
        db = "<br>".join(f"{d['name']} · {d.get('version') or (d.get('sha256','')[:12]+'…' if d.get('sha256') else d.get('file','—'))}"
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
        while layout.count():
            w = layout.takeAt(0).widget()
            if w:
                w.setParent(None)
        layout.addWidget(widget)

    def closeEvent(self, e):
        it = getattr(self, "_install_timer", None)
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
    if domains.PYHMMER_VERSION == "unavailable":
        problems.append("pyhmmer failed to load")
    if not domains.HMM_SHA256:
        problems.append("bundled Pfam HMM profiles not found")
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
    print(f"TEAGLE SELFTEST OK · primer3 {primers.PRIMER3_VERSION} · pyhmmer {domains.PYHMMER_VERSION} "
          f"· HMM {domains.HMM_SHA256[:12]} · QtSvg ok · install dialog ok")
    return 0


def main():
    if os.environ.get("TEAGLE_SELFTEST"):
        return selftest()
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
