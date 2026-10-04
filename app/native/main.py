"""TEagle native desktop app (PySide6). QMainWindow shell with a specimen rail and a scrollable
column of collapsible result cards. All science runs in-process through the shared engine, off the
GUI thread via engine_worker.Engine. This module wires the analyze workflow, primer/PCR, WSL family
annotation, splice detection, provenance and exports."""
from __future__ import annotations
import functools, gzip, html, json, math, os, re, sys

from PySide6.QtCore import Qt, QTimer, QByteArray, QSettings, QPropertyAnimation, QEasingCurve, QPoint
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
import coords
from coords import CoordTable, CoordGenomePanel
from sample import make_sample
import theme as theme_mod
from teagle_core import appdirs, classify, sequtil  # classify.DOMAINS_TESTED = the tested-profile panel (scope caveat)
from teagle_core.wsl import (curated_coverage_sentence as _curated_coverage,
                             dfam_coverage_counts as _dfam_coverage_counts)  # measured, never retyped on a panel
from teagle_core import domains as domains_mod      # DOMAIN_INFO: the methods panel is derived from it, never retyped
from teagle_core import structural as structural_mod   # MIN_TSD/MAX_TSD: the TSD method text states the real window
from teagle_core import primers as primers_mod
from teagle_core import oligoqc                      # ΔG tier thresholds: tooltips read them, never retype them


def _safe_export_stem(value, default="record", *, trust="local_input"):
    """Filesystem-safe deterministic stem for any identifier supplied by a record or backend."""
    from teagle_core.gff3 import safe_filename_stem
    return safe_filename_stem(value, default=default, trust=trust)


def _public_record_id(value):
    """Neutral outward record ID, retaining only application-owned PCR role suffixes."""
    from teagle_core.gff3 import public_record_label
    return public_record_label(value)


def _record_public_identity(record):
    """Return an already engine-validated record identity, or fail closed to a local-input identity."""
    from teagle_core.gff3 import public_record_identity
    raw = str((record or {}).get("id") or "record")
    source = (record or {}).get("identitySource")
    digest = (record or {}).get("idSha256")
    if source == "local_input" and isinstance(digest, str):
        digest = digest.lower()
        if re.fullmatch(r"[0-9a-f]{64}", digest) and raw == f"record_{digest[:12]}":
            return {"id": raw, "idSha256": digest, "identitySource": source}
    if source in {"validated_fetch", "generated"} and digest is None:
        expected = public_record_identity(raw, trust=source)
        if expected == {"id": raw}:
            return {"id": raw, "identitySource": source}
    return {**public_record_identity(raw), "identitySource": "local_input"}


def _record_public_id(record):
    return _record_public_identity(record)["id"]


def _record_export_stem(record, default="record"):
    from teagle_core.gff3 import safe_filename_stem
    identity = _record_public_identity(record)
    if identity.get("idSha256"):
        return safe_filename_stem(f"record_{identity['idSha256']}", default=default, trust="generated")
    trust = "generated" if identity.get("idSha256") else identity["identitySource"]
    return safe_filename_stem(identity["id"], default=default, trust=trust)


def _amplicon_public_source(amplicon):
    from teagle_core.gff3 import bound_public_record_label, public_record_label
    trust = amplicon.get("sourceIdentitySource")
    if trust in {"validated_fetch", "generated"}:
        return public_record_label(amplicon.get("source"), trust=trust)
    return bound_public_record_label(amplicon.get("source"), amplicon.get("sourceIdSha256"))


_ACCESSION_LIKE = re.compile(r"[A-Za-z]{1,2}_?\d{5,9}(\.\d+)?")
_NON_DNA_LETTERS = set("EFIJLOPQZ")                       # letters no nucleotide code uses (X is ambiguous: masking or protein)


def _invalid_record_text(rec, number, raw_text=""):
    """Plain-language reason one record was not analysed, naming what is wrong and where (1-based, counted
    after spaces and line breaks are removed). Built from the record's own sequence, so it is exact even when
    the engine's list of offending characters is capped."""
    label = f"Record {number}"
    reason = rec.get("invalidReason")
    if reason:
        return (f"{label} has a header but no sequence." if str(reason).startswith("empty") else
                f"{label} holds only N and IUPAC ambiguity codes, so there is no A, C, G or T to analyse.")
    seq = rec.get("seq") or ""
    bad = [(i + 1, c) for i, c in enumerate(seq) if c not in sequtil.IUPAC]
    chars = {c for _, c in bad}
    letters = [c for c in seq if c.isalpha()]
    out = [f"{label}: {len(bad):,} of {len(seq):,} characters are not nucleotide letters (A C G T U or an IUPAC "
           "code). First: " + ", ".join(f"“{c}” at {p}" for p, c in bad[:6]) + (" …" if len(bad) > 6 else "") + "."]
    protein = sum(c in _NON_DNA_LETTERS for c in letters)
    token = "".join(str(raw_text).split())
    if protein >= 5 and protein / max(1, len(letters)) >= 0.1:
        out.append("This looks like a protein sequence. TEagle analyses nucleotide sequence; paste the DNA instead.")
    elif ">" not in str(raw_text) and _ACCESSION_LIKE.fullmatch(token):
        out.append("This looks like an accession number, not a sequence. Type it in the accession box and press Fetch.")
    elif any(c.isdigit() for c in chars):
        out.append("Digits usually mean position numbers copied from a GenBank flat file. Delete the numbers, or "
                   "fetch the record by accession.")
    elif chars & set("-.~"):
        out.append("Hyphens and dots usually mark alignment gaps. Delete them to analyse the sequence.")
    elif "X" in chars:
        out.append("X usually marks masked or unknown bases. Replace it with N if the base is unknown.")
    elif "*" in chars:
        out.append("An asterisk usually marks a stop codon in a protein sequence.")
    return " ".join(out)


def invalid_input_message(records, raw_text="", warning=""):
    """Everything wrong with the records that were not analysed, in plain words. `warning` is the engine's own
    text; only its statement about characters before the first FASTA header is repeated here."""
    lines = [_invalid_record_text(r, i, raw_text) for i, r in enumerate(records, 1) if not r.get("valid")]
    shown = lines[:5] + ([f"… and {len(lines) - 5} more records."] if len(lines) > 5 else [])
    pre = [s.strip() for s in re.split(r"(?<=\.) ", str(warning or "")) if "before the first FASTA header" in s]
    return "\n".join(shown + pre)


_TE_CLASS_TEXT = {"none": "No TE signature detected"}


def _te_headline(te_class):
    """Plain headline for a classification code; a record with no evidence must not read as a bare 'none'."""
    return _TE_CLASS_TEXT.get(te_class, te_class)


def _manifest_export_stem(manifest_input, default="run"):
    inp = manifest_input or {}
    return _record_export_stem({"id": inp.get("id"), "idSha256": inp.get("idSha256"),
                                "identitySource": inp.get("identitySource")}, default)


def _SIG(fn, name):
    """A detector's actual default, for panels that DESCRIBE that detector. The methods disclosure restated
    the LTR thresholds as literal prose, so it could drift from the code silently — the TSD line already had,
    claiming a 4 bp floor against a detector searching from 2."""
    import inspect
    return inspect.signature(fn).parameters[name].default
from teagle_core.fetch import (COORD_ASSEMBLIES, all_assemblies,                        # pinned + user-added assemblies
                               complete_gene_model, cross_check_models, retrieve)        # + gene model / transcript fetch
from teagle_core import __version__ as APP_VERSION    # single source of truth (never hardcode a duplicate version)
_DFAM_YEAST_CURATED, _DFAM_YEAST_TOTAL, _DFAM_YEAST_ADDITIONAL = \
    _dfam_coverage_counts("Saccharomyces cerevisiae")
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
PANEL = {"specimen": "01", "struct": "02", "family": "03", "splice": "03b", "primer": "04", "pcr": "05",
         "genome": "06", "prov": "07"}                    # the one numbering every card header and "panel NN" pointer uses
REPO_URL = "https://github.com/tunabirgun/TEagle"        # the constant is the specification (same URL as CITATION.cff)
from engine import LARGE_INPUT_BP                         # the engine owns the documented single-locus threshold
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
# Export-only header renames. On screen the convention is a GLOSS tooltip and the columns are sized to
# contents, so spelling it out would push the table sideways at every window size; in a FILE the tooltip is
# gone and an unlabelled "Start" is ambiguous between the two conventions this app legitimately uses
# (0-based half-open for locus features, 1-based inclusive for isPcr genome hits).
EXPORT_HEADERS = {"Start": "Start (0-based)", "End": "End (0-based, exclusive)",
                  "nt": "nt (0-based)", "Coords": "Coords (0-based, half-open)"}
DOMAIN_COLS = ["Domain", "Label", "Pfam", "Strand", "ORF", "aa in ORF", "nt", "Score", "E-value", "Hit confidence"]


def _dg_tiers(caution, warn):
    """Both ΔG tiers for a structure type, read from oligoqc so a tooltip never disagrees with the cell marks."""
    return (f"! caution at ΔG ≤ −{abs(caution):g} (TEagle heuristic), !! warn at ΔG ≤ −{abs(warn):g} "
            "(IDT rule of thumb).")


def _dg_limits():
    """{axis: (caution, warn)} in kcal/mol, read from oligoqc on every call (never cached or retyped here)."""
    t = getattr(oligoqc, "THRESHOLDS", None)
    if t:
        return {k: (v["caution"], v["warn"]) for k, v in t.items()}
    return {"hairpin": (oligoqc._HAIRPIN_CAUTION, oligoqc._HAIRPIN_WARN),
            "dimer": (oligoqc._DIMER_CAUTION, oligoqc._DIMER_WARN), "end": (oligoqc._END_CAUTION, oligoqc._END_WARN)}


def _dg_clauses(minus="−", le="≤"):
    """(caution clause, warn clause) for the legend and notes, derived from _dg_limits: the caution tiers are a TEagle
    heuristic, and only the warn threshold is attributed to IDT's rule of thumb."""
    lim = _dg_limits()
    n = lambda v: f"{minus}{abs(v):g}"
    names = {"hairpin": "hairpin", "dimer": "dimers", "end": "3′-end"}
    caution = ", ".join(f"{names[k]} {le} {n(lim[k][0])}" for k in names)
    warns = {lim[k][1] for k in names}
    warn = (f"ΔG {le} {n(warns.pop())} kcal/mol on any column" if len(warns) == 1 else
            ", ".join(f"{names[k]} {le} {n(lim[k][1])}" for k in names) + " kcal/mol")
    return f"TEagle heuristic: {caution} kcal/mol", f"{warn} (IDT rule of thumb)"


PRIMER_COLS = ["ID", "Struct", "Penalty", "Location", "Notes", "Forward (5'→3')", "Reverse (5'→3')", "Product", "Fwd pos",
               "Rev pos", "Tm F/R", "GC% F/R", "Hairpin", "Self-dim", "Hetero", "3′-end"]


# plain-language glossary for every table header (hover to learn the abbreviation)
GLOSS = {
    "Record": "Record identifier used in every panel and export; a pasted header is shown as record_<SHA-256 prefix>.",
    "Classification": "TEagle class code, e.g. LTR/Copia, DNA/hAT, LINE. repeat/structural-only = terminal-repeat candidate with the "
                      "order withheld; unclassified = class not assigned; none = no TE signature detected.",
    "Call confidence": "Strength of this classification under TEagle's fixed rule — High, Moderate or Candidate; not a probability. "
                       "High is reachable only for LTR/Copia and LTR/Gypsy calls; Moderate is the top rung for every other class. "
                       "n/a when the call was withheld or no TE signature was found.",
    "Completeness": "Domain-completeness tier scored against the bundled Pfam panel. partial = some expected modules were not "
                    "detected by these profiles, not that the copy is decayed.",
    "Call status": "Whether the classification step made a call: called; withheld (the evidence does not license naming the order "
                   "or class); failed (the step raised an error); not run (invalid input). It describes the call, not the "
                   "element: a called record can still be a partial copy.",
    "Rationale": "The classifier's explanation, verbatim as written to the GFF3 and JSON exports.",
    "Structural": "Structural hallmarks detected (LTR or TDR, TIR, TSD, PBS, PPT, poly-A/T tail, polyA-signal motif).",
    "Processing": "Status of each analysis stage. domains=partial means only the longest ORFs were searched.",
    "Length": "Record length in base pairs.",
    "GC": "G+C percentage of the unambiguous A/C/G/T bases.",
    "N": "Percentage of N (unknown) bases.",
    "Ambiguous": "Percentage of IUPAC degenerate bases (R, Y, K, M, S, W, B, D, H, V); not counted as N.",
    "Input": "Whether the record passed IUPAC nucleotide validation; invalid records are not analysed.",
    "Wicker code": "Three-letter code of Wicker et al. 2007 mapped from this call (class, order, superfamily; X = rank not called), "
                   "e.g. RLC = Copia, DTT = Tc1–Mariner.",
    "RM class": "RepeatMasker/Dfam-style class string mapped from this call, only to the level TEagle called (e.g. LTR/Copia, "
                "DNA/TcMar, LINE). A mapping, not a Dfam search; — = no class asserted.",
    "Domain architecture": "Detected protein-domain codes in their inferred biological order; this is not a taxonomic rank. "
                           "none detected = the search ran and matched nothing; not assessable = no ORF of 40 aa or more to "
                           "search; not assessed = the scan failed.",
    "Feature": "Structural hallmark found in the sequence — e.g. LTR, TIR, target-site duplication, poly-A tail.",
    "Coords (0-based)": "Location in the sequence, 0-based half-open [start, end): the first base is 0, the end is not "
                        "included, length = end − start. Add 1 to the start (and keep the end) for a 1-based, inclusive "
                        "position as GenBank, SnapGene and IGV show it.",
    "Coords": "Location of this amplicon in the searched sequence, 0-based half-open [start, end).",
    "Template span": "Length of the half-open interval copied from the searched template.",
    "Product len": "Full predicted PCR-product length, including incorporated non-templated 5′ primer tails.",
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
    "Hit confidence": "Per-domain confidence from the HMMER i-Evalue (Eddy 2011): high (≤ 1e-10) or moderate. It applies to this one "
                      "hit, separately from the element's Call confidence and from domain completeness.",
    "ID": "Identifier of this designed primer pair.",
    "Forward (5'→3')": "Forward primer sequence, written 5′→3′.",
    "Reverse (5'→3')": "Reverse primer sequence, written 5′→3′.",
    "Product": "Predicted amplicon (product) size in base pairs.",
    "Tm F/R": "Melting temperature (°C) of the forward / reverse primer — matched so both anneal together.",
    "GC% F/R": "Percent G+C content of the forward / reverse primer; ~40–60% is typical.",
    "Penalty": "Primer3's overall penalty for the pair — lower is better; it rises as primers depart from target Tm, size and GC.",
    "Hairpin": "Most stable self-folding (hairpin) ΔG of the forward/reverse primer, kcal/mol (worst of the two). "
               f"More negative = more stable structure = worse. {_dg_tiers(oligoqc._HAIRPIN_CAUTION, oligoqc._HAIRPIN_WARN)} "
               "‡ = the two engines disagree — hover for both.",
    "Self-dim": "Most stable self-dimer (primer with a copy of itself) ΔG, kcal/mol (worst of F/R). "
                f"More negative = worse. {_dg_tiers(oligoqc._DIMER_CAUTION, oligoqc._DIMER_WARN)} ‡ marks an engine disagreement.",
    "Hetero": "Cross-dimer ΔG between the forward and reverse primer, kcal/mol. More negative = worse. "
              f"{_dg_tiers(oligoqc._DIMER_CAUTION, oligoqc._DIMER_WARN)} "
              "A 3′-end cross-dimer is the classic primer-dimer. ‡ marks an engine disagreement.",
    "3′-end": "3′-end anneal stability (last bases) of the pair, kcal/mol — the strongest single predictor of primer-dimer, "
              "because only a base-paired 3′ end is extended by polymerase. More negative = worse.",
    "Struct": "Worst secondary-structure flag for the pair — ok / caution / warn — across hairpin, self-dimer, hetero-dimer and 3′-end.",
    "Pair": "Which designed primer pair produced this amplicon.",
    "Mechanism": "Contributing primer identities and orientation plus the 0-based half-open product endpoints. "
                 "Several mechanisms may describe one physical product.",
    "Source": "The sequence that was searched — your specimen or a custom background.",
    "Mism F/R": "Mismatches in the forward / reverse primer binding site (the 3′ end is kept exact).",
    "Call": "For a design-bound local PCR, on-target requires the intended F+R endpoints and coverage of the "
            "Primer3 target. Without a Primer3 target, target status is unassessed; other products are additional "
            "products, not experimentally established off-targets. A genome-wide priming site has no bound locus.",
    "Class/family": "RepeatMasker/Dfam class/family of the matched Dfam family, e.g. LTR/Copia or DNA/TcMar-Tc1 — RepeatMasker's "
                    "vocabulary, which differs from Wicker 2007 at several ranks.",
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
    """One evidence-bounded call label for local PCR and genome-scan products."""
    if a.get("single_primer"):
        return "single-primer"
    if "matches_intended_product" in a:
        if a.get("matches_intended_product") is True:
            if a.get("covers_primer3_target") is True:
                return "intended product; target covered"
            if a.get("covers_primer3_target") is False:
                return "intended product; target not covered"
            return "intended product; target coverage unavailable"
        if a.get("covers_primer3_target") is True:
            return "additional; covers target"
        return ("additional product" if a.get("matches_intended_product") is False
                else "product; intended/target assessment unavailable")
    if a.get("on_target"):
        return "on-target"
    return "off-target" if has_locus else "priming site"


def _amp_kind(a, has_locus=True):
    """The same call as a punctuation-free FASTA-header token (ontarget / offtarget / singleprimer / primingsite)."""
    return re.sub(r"[^A-Za-z0-9]", "", _amp_call(a, has_locus))


def _metric_cell(parts):
    """Fold one or more ΔG metric dicts ({p3, vrna, flag, agree}) into a table cell:
    (text, worst_flag, tooltip, export_value). Shows the worst (most negative) ΔG across engines/primers, appends
    a severity mark (! caution, !! warn) so the per-axis call never rests on colour alone, and marks an engine
    disagreement with ‡; the caller maps worst_flag -> a per-theme colour. The Struct column carries the flag as
    TEXT. export_value is the BARE number (the marks are display cues) so an exported ΔG column stays numeric."""
    dgs, flags, disagree, tips, seen = [], [], False, [], False
    for lab, m in parts:
        if not m:
            continue
        p3, vr = m.get("p3"), m.get("vrna")
        vals = [v for v in (p3, vr) if v is not None and v < 0]   # a ΔG of zero or above is no favourable structure
        seen = seen or any(v is not None for v in (p3, vr))
        if vals:
            dgs.append(min(vals))
        flags.append(m.get("flag", "ok"))
        if m.get("agree") == "disagree":
            disagree = True
        pre = (lab + ": ") if lab else ""
        tips.append(f"{pre}primer3 {p3 if p3 is not None else '—'} / ViennaRNA {vr if vr is not None else '—'} kcal/mol"
                    + (f" ({m.get('agree')})" if m.get("agree") and m.get("agree") not in ("none",) else ""))
    if not dgs:
        if seen:                                          # the engine returned a value, but it is not a stable structure
            return ("none", "ok", "no stable structure (ΔG ≥ 0) — " + " · ".join(tips), "none")
        return ("—", "ok", "no structure predicted", "—")
    worst_flag = max(flags, key=lambda f: _FLAG_ORDER.get(f, 0)) if flags else "ok"
    num = f"{min(dgs):.1f}"
    txt = num + _FLAG_MARK.get(worst_flag, "") + ("‡" if disagree else "")
    return (txt, worst_flag, " · ".join(tips), num)


# clickable source citations (verified DOIs — mirror backend refs.py)
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


_QWIDGETSIZE_MAX = 16777215          # Qt's "no maximum"; a height animation must hand this back or the body freezes
CARD_ANIM_MS = 160                   # long enough to read as motion, short enough that it never gates the next click


def motion_enabled() -> bool:
    """Whether UI motion should play, honouring the OS accessibility setting.

    Qt exposes no `prefers-reduced-motion` equivalent, so read Windows' "Show animations in Windows"
    flag (SPI_GETCLIENTAREAANIMATION = 0x1042) directly. Cached: the setting changes about once in a
    user's life, and a per-click syscall is not worth it. A platform or call that cannot answer leaves
    motion ON, matching the Windows default — this is a 160 ms height ease, not something whose absence
    would be safer than its presence."""
    cached = getattr(motion_enabled, "_v", None)
    if cached is not None:
        return cached
    v = True
    try:
        import ctypes
        got = ctypes.c_int(0)
        if ctypes.windll.user32.SystemParametersInfoW(0x1042, 0, ctypes.byref(got), 0):
            v = bool(got.value)
    except Exception:                 # non-Windows, or a locked-down user32 — fall back to the default
        pass
    motion_enabled._v = v
    return v


class CollapsibleCard(QFrame):
    """A titled result card that expands/collapses on header click. Starts collapsed; reveal_on_data
    auto-expands the first time content is set (progressive reveal).

    Only a HUMAN header click animates. set_collapsed()/expand() — the programmatic API every result
    handler, test and screenshot harness drives — stays instantaneous, so a captured figure can never
    catch a card mid-ease and no test has to wait on a timer."""
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
        self._meta_default = meta
        self._collapsed = collapsed
        self.set_collapsed(collapsed)

    def _paint_hdr(self):
        arrow = "▸" if self._collapsed else "▾"
        title = self._title.replace("&", "&&")           # QPushButton eats a lone '&' as a mnemonic
        meta = self._meta.replace("&", "&&") if self._meta else ""
        self.hdr.setText(f"{arrow} {self._number}  {title}" + (f"    ·  {meta}" if meta else ""))   # sentence-case titles (calmer); ALL-CAPS reserved for micro-labels
        self.hdr.setToolTip(f"{self._number}  {self._title}" + (f" · {self._meta}" if self._meta else ""))  # full title when the gloss clips on a narrow window

    def set_meta(self, text=None):
        """Replace the header's gloss (None restores the default). Repaints the header only: never the body, the
        collapse state or a running ease."""
        self._meta = self._meta_default if text is None else text
        self._paint_hdr()

    def _toggle(self):
        self.set_collapsed(self.body.isVisible(), animate=True)

    def _ease(self, collapsed: bool):
        """Ease the body's height. Animates maximumHeight rather than the widget's own height because the
        card sits in a QVBoxLayout, which owns height directly and would fight a height animation."""
        prev = getattr(self, "_anim", None)
        if prev is not None:
            prev.stop()                                   # a fast double-click must not leave two animations racing
        if collapsed:
            start, end = self.body.height(), 0
        else:
            self.body.setVisible(True)
            self.body.setMaximumHeight(0)
            start, end = 0, max(self.body.sizeHint().height(), 1)
        anim = QPropertyAnimation(self.body, b"maximumHeight", self)
        anim.setDuration(CARD_ANIM_MS)
        anim.setStartValue(start)
        anim.setEndValue(end)
        anim.setEasingCurve(QEasingCurve.OutCubic)        # decelerate into rest; entering motion, not bouncing

        def settled():
            if collapsed:
                self.body.setVisible(False)
            # Hand the cap back unconditionally: a body pinned at its opening sizeHint could not grow when
            # a later result made it taller, which would clip the card's own content.
            self.body.setMaximumHeight(_QWIDGETSIZE_MAX)

        anim.finished.connect(settled)
        self._anim = anim
        anim.start()

    def set_collapsed(self, collapsed: bool, animate: bool = False):
        if animate and motion_enabled():
            self._ease(collapsed)
        else:
            # Cancel any in-flight human-triggered ease FIRST. Without this, a programmatic reveal that
            # lands mid-animation (a result handler calling expand(), or _scroll_to navigating to the card
            # the user just clicked) set the final state and then let the running animation keep driving
            # maximumHeight — and QPropertyAnimation.stop() does not emit `finished`, so the body could be
            # left clamped at an interpolated height with no callback to hand the cap back.
            anim = getattr(self, "_anim", None)
            if anim is not None:
                anim.stop()
                self._anim = None
            self.body.setVisible(not collapsed)
            self.body.setMaximumHeight(_QWIDGETSIZE_MAX)
        self._collapsed = collapsed                       # explicit: body.isVisible() is unreliable before show and mid-ease
        self._paint_hdr()
        self.hdr.blockSignals(True); self.hdr.setChecked(not collapsed); self.hdr.blockSignals(False)   # checked = open, however it got there
        self.hdr.setAccessibleName(f"{self._number} {self._title}")     # the ▸/▾ glyph is not speakable
        self.hdr.setAccessibleDescription("Expanded section — activate to collapse." if not collapsed
                                          else "Collapsed section — activate to expand.")

    def expand(self):
        self.set_collapsed(False)

    def clear_body(self):
        _clear_layout(self.bodylay)


def _disable_result_exports(root):
    """Fail closed every export surface retained from a result widget being removed."""
    if root is None:
        return
    buttons = ([root] if isinstance(root, QPushButton) else []) + root.findChildren(QPushButton)
    for button in buttons:
        label = button.text().strip().upper()
        if not (button.property("resultExport") or "EXPORT" in label or label in {"SVG", "PNG", "PDF"}):
            continue
        button.setEnabled(False)
        if button.property("resultExportDisabled"):
            continue
        try:
            button.clicked.disconnect()
        except (RuntimeError, TypeError):
            pass
        button.setProperty("resultExportDisabled", True)
    tables = ([root] if isinstance(root, DataTable) else []) + root.findChildren(DataTable)
    for table in tables:
        table.setEnabled(False)
        table.setContextMenuPolicy(Qt.NoContextMenu)
        if table.property("resultExportDisabled"):
            continue
        try:
            table.customContextMenuRequested.disconnect()
        except (RuntimeError, TypeError):
            pass
        table.setProperty("resultExportDisabled", True)


def _clear_layout(layout):
    """Recursively remove every widget AND nested sub-layout from `layout`. A widget-only clear
    (setParent(None) on it.widget()) leaves addLayout'd sub-layouts — the 'Structural evidence' /
    'Protein domains' section headers and their source links — orphaned but still parented to the body
    widget, so they accumulate and DUPLICATE on every re-render (e.g. a second Run analysis)."""
    while layout.count():
        it = layout.takeAt(0)
        w = it.widget()
        if w is not None:
            _disable_result_exports(w)
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


def open_external_url(url):
    """Open an http(s) link in the default browser. On Windows the launcher starts outside TEagle's
    kill-on-close job (CREATE_BREAKAWAY_FROM_JOB), so closing TEagle does not close the browser it opened;
    where the job forbids breakaway, or off Windows, the Qt desktop-services call is used instead."""
    url = str(url)
    if not re.match(r"https?://", url, re.IGNORECASE):
        return False                                          # never hand another scheme to the shell
    if sys.platform == "win32":
        import subprocess
        try:
            subprocess.Popen(["rundll32.exe", "url.dll,FileProtocolHandler", url], close_fds=True,
                             creationflags=0x01000000 | 0x08000000)   # CREATE_BREAKAWAY_FROM_JOB | CREATE_NO_WINDOW
            return True
        except OSError:
            pass
    from PySide6.QtCore import QUrl
    from PySide6.QtGui import QDesktopServices
    return QDesktopServices.openUrl(QUrl(url))


def _external(lab):
    """Route a label's links through open_external_url instead of Qt's own opener."""
    lab.setOpenExternalLinks(False)
    lab.linkActivated.connect(open_external_url)
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
    """A section label — mono, uppercase, tracked."""
    l = QLabel(widgets.caps(text)); l.setObjectName("sectionlabel"); l.setWordWrap(True); return l


def _export_table_btn(table, base, parent):
    """Visible Excel/CSV/TSV export for a DataTable: the button pops a format menu, then a save dialog
    pre-set to the chosen type. Exports in the table's current on-screen (sorted) order."""
    b = QPushButton("Export table"); b.setProperty("sm", True); b.setProperty("resultExport", True)
    table.export_base = base                          # the right-click Export… proposes the SAME filename as this button
    def pop():
        fmt = widgets.pick_table_format(b, b.mapToGlobal(b.rect().bottomLeft()))
        if fmt:
            widgets.export_table(table.export_headers(), table.rows_data(), base, parent, fmt=fmt,
                                 notes=getattr(table, "export_notes", None))
    b.clicked.connect(lambda _=False: pop())
    return b


def _fasta_text(records, line_width=60, identifiers_are_public=False):
    """Deterministic FASTA for already-derived records: (identifier, sequence, optional description)."""
    chunks = []
    for item in records:
        ident, sequence = item[0], item[1]
        description = item[2] if len(item) > 2 else ""
        supplied_digest = item[3] if len(item) > 3 else None
        ident = str(ident)
        ident = ident[1:] if ident.startswith(">") else ident
        if identifiers_are_public or supplied_digest:
            if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._:-]*", ident):
                raise ValueError("derived FASTA identifier is not safe")
        else:
            from teagle_core.gff3 import public_record_identity
            identity = public_record_identity(ident)
            ident = identity["id"]
            supplied_digest = supplied_digest or identity.get("idSha256")
        digest = None
        if supplied_digest:
            digest = str(supplied_digest).lower()
            if (not re.fullmatch(r"[0-9a-f]{64}", digest) or
                    ident != f"record_{digest[:12]}" and not ident.startswith(f"record_{digest[:12]}_")):
                raise ValueError("record identity SHA-256 does not match the public FASTA identifier")
        ident = re.sub(r"\s+", "_", ident) or "record"
        sequence = re.sub(r"\s+", "", str(sequence)).upper()
        if not sequence:
            continue
        if supplied_digest:
            description = (str(description).strip() + f" record_sha256={digest}").strip()
        header = ident + ((" " + str(description).strip()) if str(description).strip() else "")
        chunks.append(">" + header)
        chunks.extend(sequence[i:i + line_width] for i in range(0, len(sequence), line_width))
    return "\n".join(chunks) + ("\n" if chunks else "")


_COORD_SYSTEM = "coordinate_system=0-based-half-open"       # the convention every coords=/nt= header token uses


def _orfs_fasta(rec):
    seq = rec.get("seq") or ""
    identity = _record_public_identity(rec)
    rid, digest = identity["id"], identity.get("idSha256")
    records = []
    for i, feature in enumerate(rec.get("orfs") or [], 1):
        start, end = feature["start"], feature["end"]
        dna = seq[start:end]
        if feature.get("strand") == "-":
            dna = sequtil.reverse_complement(dna)
        records.append((f"{rid}_ORF_{i}", dna,
                        f"coords={start}-{end} strand={feature.get('strand', '+')} frame={feature.get('frame', '')} "
                        f"{_COORD_SYSTEM}", digest))
    return _fasta_text(records, identifiers_are_public=True)


def _domains_fasta(rec):
    identity = _record_public_identity(rec)
    rid, digest = identity["id"], identity.get("idSha256")
    return _fasta_text([
        (f"{rid}_domain_{i}_{domain.get('domain') or 'domain'}", domain.get("protein") or "",
         f"nt={domain.get('nt', ['?', '?'])[0]}-{domain.get('nt', ['?', '?'])[1]} "
         f"strand={domain.get('strand', '+')} pfam={domain.get('pfam') or 'NA'} {_COORD_SYSTEM}", digest)
        for i, domain in enumerate(rec.get("domains") or [], 1)
        if domain.get("protein")
    ], identifiers_are_public=True)


def _splice_exons_fasta(result, genomic_sequence, record_id="splice_locus", record_id_sha256=None,
                         identity_source="generated"):
    seq = re.sub(r"\s+", "", genomic_sequence or "").upper()
    identity = _record_public_identity({"id": record_id, "idSha256": record_id_sha256,
                                        "identitySource": identity_source})
    record_id, digest = identity["id"], identity.get("idSha256")
    strand = result.get("strand", "+")
    # transcript order: on the minus strand exon_1 is the 3'-most block on the record, the 5' exon of the mRNA
    exons = sorted(result.get("exons") or [], key=lambda x: x["start"], reverse=(strand == "-"))
    records = []
    for i, exon in enumerate(exons, 1):
        start, end = exon["start"], exon["end"]
        dna = seq[start:end]
        if strand == "-":
            dna = sequtil.reverse_complement(dna)
        records.append((f"{record_id}_exon_{i}", dna,
                        f"coords={start}-{end} strand={strand} transcript_rank={i} {_COORD_SYSTEM}", digest))
    return _fasta_text(records, identifiers_are_public=True)


def _full_records_fasta(records):
    bound = []
    for record in records or []:
        identity = _record_public_identity(record)
        sequence = record.get("seq") or ""
        bound.append((identity["id"], sequence, f"length={len(sequence)}", identity.get("idSha256")))
    return _fasta_text(bound, identifiers_are_public=True)


def _selftest_positive_control():
    """Problems found when the bundled copia reference is analysed end to end; an empty list is healthy. The
    release self-test used a random-ORF negative control, which still passes when the domain search is broken:
    a positive control must find its known RT/INT domains with every stage complete."""
    import engine
    from teagle_core import examples, provenance
    if not examples.available():
        return ["no bundled example sequences were found in the build"]
    problems = []
    r = engine.analyze(examples.load("M11240"))
    rec = (r.get("records") or [{}])[0]
    for stage in ("structural", "orfs", "domains", "classification"):
        status = ((rec.get("processing") or {}).get(stage) or {}).get("status")
        if status != "complete":
            problems.append(f"positive control: the {stage} stage reported {status!r}, expected 'complete'")
    te_class = (rec.get("classification") or {}).get("te_class")
    if te_class != "LTR/Copia":
        problems.append(f"positive control: the bundled copia reference was called {te_class!r}, expected 'LTR/Copia'")
    if not {d.get("domain") for d in rec.get("domains") or []} & {"RT", "INT"}:
        problems.append("positive control: no RT or INT domain was detected in the bundled copia reference")
    try:
        engine.analysis_result_json_bytes(r)
        provenance.manifest_json_bytes(r.get("provenance") or {})
    except Exception as e:
        problems.append(f"positive control: result or manifest did not serialise: {type(e).__name__}: {e}")
    return problems


def _export_fasta_btn(label, fasta_fn, base, parent):
    """Visible, invalidation-aware FASTA export backed by a callable so it cannot drift from current data."""
    b = QPushButton(label); b.setProperty("sm", True); b.setProperty("resultExport", True)
    b.clicked.connect(lambda _=False: widgets.save_fasta(fasta_fn(), base, parent))
    return b


def _qq_note():
    """RepeatMasker runs in -qq (rush) mode. That was an invisible speed choice while only a family NAME
    was shown; once its divergence and coverage are displayed as curation evidence the search sensitivity
    becomes part of the claim, so it is stated beside the numbers rather than buried in the manifest."""
    return _note("Search run in RepeatMasker's rush mode (-qq): about 10% less sensitive than a default search "
                 "and 4 to 10 times faster (RepeatMasker's own description). Divergence is the raw substitution percentage, not Kimura-corrected, and coverage "
                 "is of the Dfam consensus — both are alignment evidence, not a family-membership test.",
                 "info")


# ---------------------------------------------------------------------------------------------------------------
# Display views of one analysed record. Pure functions of the record: they return new strings and never write
# into the sealed result, so analysis_result_document(res) still passes after any render.
# ---------------------------------------------------------------------------------------------------------------
_STATUS_WORD = {"complete": "called", "withheld": "withheld", "failed": "failed", "not_run": "not run"}
_MIN_ORF_AA = _SIG(sequtil.find_orfs, "min_aa")


def _call_view(rec):
    """One state vocabulary for the record summary and the banner. status: called | withheld | failed | not run
    (the classification STEP, not a statement about the element); confidence: the chip level or None; cell: what the
    summary's 'Call confidence' column shows. The chip and the cell cannot disagree because both come from here."""
    cl = rec.get("classification") or {}
    stage = ((rec.get("processing") or {}).get("classification") or {}).get("status")
    if not rec.get("valid", True):
        status = "not run"
    else:
        status = _STATUS_WORD.get(stage) or (_STATUS_WORD.get(cl.get("outcome")) if cl else None) or \
                 ("called" if cl else "not run")
    te, conf = cl.get("te_class"), cl.get("confidence")
    chip = conf if (status == "called" and te != "none" and conf) else None
    cell = (conf if chip else "n/a — call withheld" if status == "withheld" else
            "n/a — no TE signature" if status == "called" and te == "none" else "n/a")
    return {"status": status, "confidence": chip, "cell": cell}


_CONFIDENCE_TIP = (
    "Call confidence is TEagle's fixed evidence rule for this classification, not a probability.\n"
    "High — only an LTR/Copia or LTR/Gypsy call with paired LTRs, at least two domain types and a readable "
    "integrase-vs-RT order (and no co-occurring transposase).\n"
    "Moderate — the highest rung for every other class: a LINE with a 3′ tail; a DNA transposon whose transposase lies "
    "between its TIRs, or that carries two domain types; DIRS; an RT with terminal repeats; a Copia/Gypsy call whose "
    "INT/RT order cannot be read.\n"
    "Candidate — a single line of evidence, conflicting transposase families or INT/RT order, or a withheld or "
    "unassessed call.\n"
    "MODERATE on a LINE or DNA transposon is that class's ceiling, not doubt about the call. The domain table's Hit "
    "confidence is a separate, per-hit measure.")


# Wicker et al. 2007 (doi:10.1038/nrg2165) three-letter codes (X = rank not called; TREP lists RLC, RLG, RLX, RYD, RIX,
# RXX, DTT, DTA, DTM, DTC, DTX, DXX) and the RepeatMasker/Dfam Type/SubType strings of the Dfam 4.0 curated library.
# The table IS the specification of a vocabulary mapping, so it is written out; a test cross-checks it against
# benchmarks/score.py. Keyed on te_class: (order, superfamily, code, RepeatMasker class).
_WICKER = {
    "LTR/Copia": ("LTR", "Copia", "RLC", "LTR/Copia"),
    "LTR/Gypsy": ("LTR", "Gypsy", "RLG", "LTR/Gypsy"),
    "LTR/unclassified": ("LTR", None, "RLX", "LTR"),
    "LTR/partial": ("LTR", None, "RLX", "LTR"),
    "DIRS": ("DIRS", None, "RYX", "LTR/DIRS"),
    "LINE": ("LINE", None, "RIX", "LINE"),
    "retro/unassigned": (None, None, "RXX", None),
    "DNA/hAT": ("TIR", "hAT", "DTA", "DNA/hAT"),
    "DNA/Tc1-Mariner": ("TIR", "Tc1–Mariner", "DTT", "DNA/TcMar"),
    "DNA/CACTA": ("TIR", "CACTA", "DTC", "DNA/CMC-EnSpm"),
    "DNA/MULE": ("TIR", "Mutator", "DTM", "DNA/MULE-MuDR"),
    "DNA/IS4-like": ("TIR", None, "DTX", "DNA"),
    "DNA/DDE": ("TIR", None, "DTX", "DNA"),
    "DNA/unclassified": ("TIR", None, "DTX", "DNA"),
    "DNA/non-autonomous": ("TIR", None, "DTX", "DNA"),
}


def _rank_view(cl):
    """TEagle's call on Wicker 2007's ranks and in RepeatMasker's vocabulary, or None when no class is derived. A display
    mapping only: it goes no deeper than the rank TEagle called, and an unknown te_class gets an X code and no string."""
    cl = cl or {}
    te = str(cl.get("te_class") or "").strip()
    if not te or te.lower() in ("unclassified", "none", "-"):
        return None
    head = te.split("/")[0].strip().lower()
    klass = str(cl.get("class") or "").lower()
    order = ("LTR" if head == "ltr" else "LINE" if head == "line" else "SINE" if head == "sine" else
             "TIR" if head in ("tir", "dna", "mite", "helitron") else "ABSTAIN" if head == "retro" else
             "NO_CALL" if head == "repeat" else "OTHER")
    if "class ii" in klass or order == "TIR":
        rank = "Class II"
    elif "class i" in klass or order in ("LTR", "LINE", "SINE"):
        rank = "Class I"
    else:
        return None
    if cl.get("is_erv"):
        return {"class": rank, "order": "LTR", "superfamily": None, "code": "RLX", "rm": "LTR", "erv": True}
    got = _WICKER.get(te)
    if got is None:
        return {"class": rank, "order": None, "superfamily": None, "code": "DXX" if rank == "Class II" else "RXX",
                "rm": None, "erv": False}
    return {"class": rank, "order": got[0], "superfamily": got[1], "code": got[2], "rm": got[3], "erv": False}


def _rank_line(view):
    """The banner text for a _rank_view."""
    if view["erv"]:
        sf = "superfamily not resolved — env-bearing"
    elif view["order"] is None:
        sf = None
    else:
        sf = f"{view['superfamily']}" if view["superfamily"] else "superfamily not called"
    path = f"{view['class']} › " + (f"{view['order']} › {sf}" if view["order"] else "order not called")
    rm = view["rm"] or "— (none asserted)"
    return f"Wicker 2007: {path} ({view['code']})   ·   RepeatMasker class: {rm}"


_RANK_TIP = ("A display mapping of TEagle's call onto Wicker et al. 2007 (three-letter code; X = rank not called) and "
             "RepeatMasker's class/family vocabulary. No new call is made. For a Dfam family name use panel 03.")
_ERV_TIP = (" Env-bearing LTR elements fall in Wicker's Gypsy superfamily (errantiviruses) or its Retrovirus/ERV "
            "superfamilies (RLR/RLE); TEagle cannot tell these apart, so the INT/RT-order call is not mapped to a "
            "superfamily here.")


# The five bundled examples carry their published class in TEagle's te_class vocabulary. The constants are the
# specification (the published identity); a test pins their keys to examples.CATALOG and the golden classifications.
_PUBLISHED_CLASS = {"M11240": "LTR/Copia", "M12927": "LTR/Gypsy", "M80343": "LINE", "X01005": "DNA/Tc1-Mariner",
                    "X05424": "DNA/hAT"}


@functools.lru_cache(maxsize=None)
def _example_digests():
    """{sha256 of a bundled example's record sequence: catalog entry + published class}. The binding is exact sequence
    identity: an edited, trimmed or reverse-complemented copy of an example is not that example."""
    from teagle_core import examples
    out = {}
    for entry in examples.available():
        parsed = sequtil.parse_fasta(examples.load(entry[0]) or "")
        if parsed and entry[0] in _PUBLISHED_CLASS:
            out[sequtil.sha256(parsed[0][1])] = (*entry[:4], _PUBLISHED_CLASS[entry[0]])
    return out


def _example_match(seq):
    return _example_digests().get(sequtil.sha256(seq)) if seq else None


@functools.lru_cache(maxsize=None)
def _example_accession(acc):
    """Versioned GenBank accession of a bundled example (the first word of its FASTA header), or `acc` itself."""
    from teagle_core import examples
    parsed = sequtil.parse_fasta(examples.load(acc) or "")
    tok = parsed[0][0].split()[0] if parsed else ""
    return tok if tok == acc or tok.startswith(acc + ".") else acc


def _published_note(rec):
    """(text, level) for a record that IS a bundled example: the published class beside TEagle's call, or None."""
    entry = _example_match(rec.get("seq"))
    cl = rec.get("classification") or {}
    if not entry:
        return None
    acc, label, organism, expected, code = entry
    head = f"Bundled example {acc} — {label}, {organism}. Published as: {expected}. "
    te = cl.get("te_class")
    if cl.get("outcome") == "withheld" or not cl or te in (None, "unclassified", "none"):
        return head + "TEagle withheld the order call for this record, so it is not compared.", "info"
    if te == code:
        return (head + f"TEagle's call {te} is consistent with the published class. The completeness tier below reports "
                "what the tested Pfam panel detected in this sequence, not whether the published element is complete."), "success"
    return head + (f"TEagle's call {te} differs from the published class ({code}). Check the evidence below before "
                   "relying on this call."), "warn"


@functools.lru_cache(maxsize=None)
def _panel_references():
    """{te_class: (accession, label, organism, tier, missing modules, sequence digest)} for the bundled examples, run once
    per session, so the statement 'the bundled reference reads the same' is measured, not asserted."""
    import engine
    from teagle_core import examples
    out = {}
    try:
        for acc, label, organism, _expected in examples.available():
            rec = engine.analyze(examples.load(acc))["records"][0]
            comp = (rec.get("classification") or {}).get("completeness") or {}
            out[(rec.get("classification") or {}).get("te_class")] = (
                acc, label, organism, comp.get("tier"), tuple(comp.get("missing") or ()), sequtil.sha256(rec["seq"]))
    except Exception:
        return {}
    return out


def _no_orf_to_search(rec):
    """True only when the ORF stage completed and found none, so the domain search had nothing to read."""
    p = rec.get("processing") or {}
    return ((p.get("orfs") or {}).get("status") == "complete" and (p.get("orfs") or {}).get("count") == 0 and
            (p.get("domains") or {}).get("status") in ("complete", "partial", "not_run"))


def _orf_budget_tie(rec):
    """(scanned, n_orfs, unsearched, tied, aa) when ORFs beyond the domain search's budget are exactly as long as the
    shortest one it read, or None. Equal lengths are cut by search order (strand, frame, position), not by content."""
    scan = (rec.get("processing") or {}).get("domains") or {}
    orfs, scanned = rec.get("orfs") or [], scan.get("orfsScanned")
    if not scanned or not scan.get("orfsNotScanned") or len(orfs) <= scanned:
        return None
    aa = orfs[scanned - 1].get("length_aa")
    tied = sum(1 for o in orfs[scanned:] if o.get("length_aa") == aa)
    return (scanned, len(orfs), len(orfs) - scanned, tied, aa) if tied else None


def _repeated_domain_set(rec):
    """(k, codes) when k >= 2 separate, non-overlapping ORFs each carry the same set of two or more domain codes that
    includes RT, INT or TPase; None otherwise. A display heuristic for 'probably more than one element copy': a single-code
    set is excluded on purpose (a frameshifted transposase also splits into one-code ORFs), accepting false negatives."""
    orfs, hosts = rec.get("orfs") or [], {}
    for d in rec.get("domains") or []:
        i = d.get("orf")
        if isinstance(i, int) and not isinstance(i, bool) and 0 <= i < len(orfs):
            hosts.setdefault(i, []).append(d)
    groups = {}
    for i, ds in hosts.items():
        codes = frozenset(d["domain"] for d in ds)
        if len(codes) >= 2 and codes & {"RT", "INT", "TPase"}:
            groups.setdefault(codes, []).append(i)
    best = None
    for codes, idx in groups.items():
        k, end = 0, -1
        for i in sorted(idx, key=lambda j: orfs[j]["start"]):
            if orfs[i]["start"] >= end:
                k, end = k + 1, orfs[i]["end"]
        if k >= 2 and (best is None or k > best[0]):
            first = hosts[min(idx, key=lambda j: orfs[j]["start"])]
            order = list(dict.fromkeys(d["domain"] for d in sorted(first, key=lambda d: d["nt"][0])))
            best = (k, order)
    return best


def _repeat_word(rec):
    """'TDR' for an order-withheld direct repeat, 'LTR' or 'TIR' for an accepted pair, None when the record has neither."""
    from teagle_core import gff3
    types = [str(e.get("type", "")) for e in rec.get("structural") or [] if not e.get("advisory")]
    if any(t.startswith("LTR") for t in types):
        return "TDR" if gff3.is_neutral_tdr(rec.get("classification")) else "LTR"
    return "TIR" if any(t.startswith("TIR") for t in types) else None


def _caution_items(rec):
    """Conditions that limit what the single classification above can mean, from structured fields only. The engine's
    note text, the confidence value and the thresholds are unchanged: this lifts what is already known into view."""
    out = []
    n = (rec.get("composition") or {}).get("length") or len(rec.get("seq") or "")
    if n > LARGE_INPUT_BP:
        out.append(f"Large input ({n:,} bp). The detectors are tuned for a single element, so this one call may span "
                   f"several elements. To annotate each element, use panel {PANEL['family']} (RepeatMasker against Dfam).")
    scan = (rec.get("processing") or {}).get("domains") or {}
    rep = _repeated_domain_set(rec)
    if rep:
        k, codes = rep
        word = _repeat_word(rec)
        s = (f"Multiple copies likely: {k} separate, non-overlapping ORFs each carry the same {' + '.join(codes)} domain "
             "set, so this sequence probably holds more than one element"
             + (f" and the {word} coordinates describe only the outermost repeat pair" if word else "")
             + ". (TEagle display heuristic — the classification itself is unchanged.)")
        if scan.get("orfsNotScanned"):
            s += (f" {k} counted among the {scan.get('orfsScanned')} ORFs searched; the {scan['orfsNotScanned']} "
                  "unsearched ORFs may hold more.")
        out.append(s)
    tie = _orf_budget_tie(rec)
    if tie:
        scanned, total, _u, t, aa = tie
        out.append(f"The domain search read the {scanned} longest of {total} ORFs, but {t} unsearched ORF(s) are exactly as "
                   f"long as the shortest one searched ({aa} aa). Equal lengths are cut by search order (strand, frame, "
                   "position), not by content, so domains in those ORFs were not looked for.")
    coding = ((rec.get("classification") or {}).get("completeness") or {}).get("coding") or {}
    if coding.get("truncated_frame") is True:
        out.append("A domain-bearing ORF runs to the end of the supplied sequence — the element may continue beyond it.")
    return out


def _note_view(note, rec):
    """The engine's ORF-budget note, reworded at display time only when the budget cut a tie."""
    tie = _orf_budget_tie(rec)
    old = f"{tie[2]} shorter ORF(s) were not searched" if tie else None
    if not tie or old not in note:
        return note
    return note.replace(old, f"{tie[2]} further ORF(s) were not searched — {tie[3]} of them exactly as long as the shortest "
                             f"searched ORF ({tie[4]} aa); equal lengths are cut by search order, not by content")


def _cap_first(s):
    return s[:1].upper() + s[1:]


def _domain_segment(rec):
    """'INT ×12 · RT ×12' for the domains found: repeats collapsed in first-appearance order, with the profile names
    added when one code is hit by different profiles (the three ORF1p models, Tc1's HTH and DDE)."""
    seen = {}
    for d in rec.get("domains") or []:
        seen.setdefault(d["domain"], []).append(d.get("hmm") or "")
    parts = []
    for code, hmms in seen.items():
        part = code + (f" ×{len(hmms)}" if len(hmms) > 1 else "")
        distinct = list(dict.fromkeys(h for h in hmms if h))
        if len(hmms) > 1 and len(distinct) > 1:
            part += f" ({', '.join(distinct)})"
        parts.append(part)
    return " · ".join(parts)


def _explanation_view(rec):
    """The classifier's explanation, rewritten for reading: collapsed repeats, a sentence that starts with a capital, the
    no-ORF case stated as not assessable, TDR for an order-withheld repeat, and 'further' on a tied ORF budget. Each
    rule is an exact-substring replacement built from the record's own data; when the expected text is not there (the
    classifier's wording drifted) the sentence is left as written. cl['explanation'] itself is never touched."""
    cl = rec.get("classification") or {}
    why = re.sub(r"^\s*Classified as [^.]+\.\s*", "", str(cl.get("explanation") or ""))
    from teagle_core import gff3
    if gff3.is_neutral_tdr(cl):
        why = re.sub(r"Structural: [^.]*\.", lambda m: re.sub(r"\bLTR\b", "TDR", m.group(0)), why, count=1)
    raw = " · ".join(d["domain"] for d in rec.get("domains") or [])
    if raw and f"Domains: {raw}." in why:
        why = why.replace(f"Domains: {raw}.", f"Domains: {_domain_segment(rec)}.", 1)
    elif not raw and "Domains: none detected in tested panel." in why and _no_orf_to_search(rec):
        why = why.replace("Domains: none detected in tested panel.",
                          f"Domains: not assessable (no ORF ≥ {_MIN_ORF_AA} aa to search).", 1)
    m = re.search(r"Domains: [^.]*\. ([a-z])", why)
    if m:
        why = why[:m.start(1)] + m.group(1).upper() + why[m.end(1):]
    if _orf_budget_tie(rec):
        why = why.replace("shorter ORF(s) were outside", "further ORF(s) were outside", 1)
    return why


_COORD_TAIL = ("consistent with an interrupted, rearranged or degenerate copy, though a genuine multi-ORF architecture "
               "(LINE ORF1/ORF2, some ERVs) also looks like this")


def _coding_view(rec):
    """The classifier's coding-structure sentence with a capital, and — when the ORFs carrying the detected domains hold the
    same domain set more than once — the multi-copy reading instead of the 'degenerate copy' one. Verbatim otherwise."""
    note = str(((rec.get("classification") or {}).get("completeness") or {}).get("coding", {}).get("note") or "")
    rep = _repeated_domain_set(rec)
    if rep and _COORD_TAIL in note:
        note = note.replace(_COORD_TAIL, f"{rep[0]} of these ORFs carry the same {' + '.join(rep[1])} domain set, which "
                                         "suggests more than one element copy (see the scope caution above)")
    return _cap_first(note)


def _saved_coord_base():
    """The remembered display convention for coordinates: 1 (1-based, inclusive) or 0; anything else reads as 0."""
    try:
        return 1 if int(QSettings("TEagle", "TEagle").value("coord_base", 0)) == 1 else 0
    except (TypeError, ValueError):
        return 0


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
        self.coord_base = _saved_coord_base()                 # on-screen coordinate convention only; files keep 0-based half-open
        self.state = {"seq": "", "source": None, "last_rec": None}
        self._loading = False                                 # True while a programmatic load writes the specimen box
        self._pcr_run = None
        self._pcr_gen = 0                                    # monotonic in-silico-PCR batch id (drops stale sibling results)
        self._spec_gen = 0                                    # bumped by every invalidation: a primer/PCR result sent under an older one is dropped
        self._specimen_gen = 0                                # bumped only when the specimen itself changes (family, splice, scan, fetch)
        self._recTimer = QTimer(self); self._recTimer.setSingleShot(True); self._recTimer.setInterval(150)
        self._recTimer.timeout.connect(self._apply_record_selection)   # debounces a run of arrow-key row changes
        self._recTable = self._recHead = self._recDetail = None        # the persistent record table, its heading, the rebuilt detail block
        self._rec_pending = None
        self._req_n = 0
        self._reqs = {}                                       # kind -> newest live tagged request {key, gen, digest}
        self._box_sig = ""                                    # whitespace-free text of the specimen box (a pure whitespace edit is no change)
        self._wsl_ready = {"annotate": False, "splice": False}
        self._analyze_gen = 0                                 # monotonic analysis id (drops stale/out-of-order callbacks)
        self._active_analyze_gen = None
        self._analysis_requests = {}
        self._analyze_inflight = False                        # gates Run analysis alongside "is there a sequence yet?"
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
        outer.setContentsMargins(theme_mod.sp(12), theme_mod.sp(10), theme_mod.sp(12), theme_mod.sp(10))
        outer.setSpacing(theme_mod.sp(8))
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
        self.acc.setFocus(Qt.OtherFocusReason)                 # the first thing to do is name a sequence, not hide a panel
        QTimer.singleShot(0, self, self._startup)             # context object: a window freed first never runs it

    # ---------- header ----------
    def resizeEvent(self, e):
        # hide the decorative tagline on a narrow window so it never contributes to an overhang; the Ignored
        # size policy already lets the window shrink, this just keeps the header clean at small widths.
        super().resizeEvent(e)
        tag = getattr(self, "_tagline", None)
        if tag is not None:
            tag.setVisible(self.width() >= round(1080 * theme_mod.UI_SCALE))
        if hasattr(self, "coordBtn"):
            self._paint_coord_btn()

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

    def _set_coord_base(self, base, persist=True):
        """Switch the on-screen numbering and redraw what is showing, without rebinding the record or touching a running
        design (so it re-renders the cards, never _show_record)."""
        base = 1 if base == 1 else 0
        if base == self.coord_base:
            return
        self.coord_base = base
        if persist:
            QSettings("TEagle", "TEagle").setValue("coord_base", base)
        self._paint_coord_btn()
        st = self.state
        if st.get("last_rec") and st.get("analysis_state") == "current":
            self._render_struct_card(st["last_rec"], st.get("analysis") or {}, keep_table=True)
        if st.get("lastPrimers") and not self._design_inflight:
            self._render_primers(st["lastPrimers"])
        if st.get("lastPcr"):
            lp = st["lastPcr"]
            self._render_pcr(lp["lanes"], lp["amplicons"], lp["statuses"])
        if st.get("splice"):
            self._render_splice(st["splice"])
        if st.get("familyResult"):
            self._render_family(st["familyResult"])
        self._refresh_assay_goal()
        self._uppercase_buttons()
        self._banner("Coordinates are now shown 1-based, inclusive. Exported tables, FASTA names and primer-design coordinates "
                     "are unchanged (0-based, half-open)." if base else
                     "Coordinates are now shown 0-based, half-open — the convention TEagle's exported tables use.", "info")

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

    # ---------- quick start, help and glossary ----------
    def _quickstart_hidden(self):
        return str(QSettings("TEagle", "TEagle").value("quickstart_hidden", False)).lower() in ("true", "1")

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

    def _set_quickstart_hidden(self, flag):
        QSettings("TEagle", "TEagle").setValue("quickstart_hidden", bool(flag))
        self._quickstart.setVisible(not flag)

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

    def _open_manual(self):
        return open_external_url(f"{REPO_URL}/releases/tag/v{APP_VERSION}")

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

    # ---------- 03 Dfam / RepeatMasker family (WSL) ----------
    def _build_wsl_card(self):
        card = CollapsibleCard(PANEL["family"], "Dfam / RepeatMasker family",
                               "optional · Linux backend: checking… · names the TE family (Dfam)")
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
        self.wslPaste.setTabChangesFocus(True)
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

    def _build_genome_card(self):
        # dedicated whole-genome off-target scan — LOCAL isPcr against a downloaded RefSeq assembly. The organism
        # dropdown lists ONLY downloaded+verified genomes (from genome_list). On-target = the product at the
        # specimen's own genome locus when it sits in the scanned assembly; the rest are off-target candidate sites.
        card = CollapsibleCard(PANEL["genome"], "Whole-genome off-target scan",
                               "optional · Linux backend: checking… · a designed pair vs a downloaded RefSeq genome (isPcr)")
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
                       "genome. Results are sealed with the assembly version + checksum. With a design locus, "
                       "intended-locus sites are listed before additional sites.")
        gnote.setObjectName("orient"); gnote.setWordWrap(True); card.bodylay.addWidget(gnote)
        self.genomeBody = QVBoxLayout(); gb = QWidget(); gb.setLayout(self.genomeBody)
        self.genomeBody.addWidget(_empty("Design a primer pair (panel 04) and download a genome, then scan here."))
        card.bodylay.addWidget(gb)
        return card

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
        # A cache hit can still trigger a feature re-fetch with retry back-off, so the fetch runs off the GUI thread
        # and the window stays responsive; the result is applied here on the GUI thread when it lands.
        import threading
        self.spliceTxUse.setEnabled(False)
        out, specimen = {}, self._specimen_gen

        def work():
            try:
                out["fasta"] = retrieve(acc).get("fasta", "")
            except Exception as e:
                out["error"] = e.args[0] if getattr(e, "args", None) else e
        threading.Thread(target=work, daemon=True).start()
        poll = QTimer(self); poll.setInterval(80)

        def done():
            if not out:
                return
            poll.stop(); poll.deleteLater()
            self.spliceTxUse.setEnabled(True)
            if specimen != self._specimen_gen:                # the specimen changed while this was in flight
                return
            if "error" in out:
                return self._banner(f"Could not fetch {acc}: {out['error']}", "warn")
            self._tx_programmatic = True
            self.spliceTx.setPlainText(out["fasta"])
            self._tx_programmatic = False
            self._splice_tx_origin = "record"                 # same annotation source -> consistency check, not independent
            self.card_splice.expand(); self.spliceTx.setFocus()
        poll.timeout.connect(done); poll.start()

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

    @staticmethod
    def _structure_text(c):
        qc = c.get("qc") or {}
        if not qc.get("ok"):
            return f"{MainWindow._pair_label(c)}: secondary-structure QC unavailable"
        L, R = qc["left"], qc["right"]
        def line(name, m):
            return f"{name}\tprimer3 {m.get('p3')}\tViennaRNA {m.get('vrna')}\t{m.get('flag')}\t{m.get('agree')}"
        return "\n".join([f"# {MainWindow._pair_label(c)} secondary structure (ΔG kcal/mol)",
                          f"# F {c['left_seq']} / R {c['right_seq']}",
                          line("hairpin_F", L["hairpin"]), line("hairpin_R", R["hairpin"]),
                          line("selfdimer_F", L["self_dimer"]), line("selfdimer_R", R["self_dimer"]),
                          line("heterodimer", qc["hetero_dimer"]), line("end_stability", qc["end_stability"])])

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

    _SCAN_PROD_MAX = 4000            # isPcr's product-size window for a genome scan; shown with every scan result
    _SCAN_ROWS_SHOWN = 500           # rows drawn in the scan table; an export always writes every product

    def _scan_genome(self, cand, org=None):
        """Whole-genome scan of one primer pair via LOCAL isPcr against a downloaded RefSeq assembly. Runtime
        depends on genome size, the installed isPcr build and cached RefSeq assembly, local hardware, and the
        fixed candidate-product-size cap. If the genome is not cached, the backend replies need_prepare and we
        offer a one-time download. Products are candidate priming sites, not validated bands; on/off-target
        labels are applied only when the specimen supplies a matching intended locus.
        `org` is captured at scan start so a later dropdown change can't retarget an in-flight scan/download."""
        self._clear_banner()
        if self._genome_inflight or self._genome_prep_inflight:
            return self._banner("A whole-genome scan / download is already running — wait for it to finish.", "info")
        run = getattr(self, "_pcr_run", None)
        if run and any(r is None for r in run.get("results", [])):
            return self._banner("An in-silico PCR run is still in progress — wait for it before a genome scan.", "info")
        if org is None:
            org = self.genomeOrg.currentData()
        if not org:
            return self._banner("Select an organism in the whole-genome off-target scan card, or download one first.", "info")
        self._pending_scan = {"cand": cand, "org": org}        # cand + org captured so a download resumes the SAME target
        self._scan_gen = self._specimen_gen                    # the specimen this scan was started for (see _scan_is_stale)
        self._genome_inflight = True
        self._refresh_genome_manager()                         # an open manager must lock its row actions while the scan runs
        self.card_genome.expand()
        # fixed wide window (isPcr defaults) — decoupled from the local-PCR size fields the user never sets for a
        # right-click genome scan. Widened only when the designed product itself is longer than the usual cap, so a
        # pair designed for a long amplicon is not reported absent; the window in force is sealed and shown.
        size = cand.get("product_size")
        prod_max = int(size * 1.25) if isinstance(size, int) and size > self._SCAN_PROD_MAX else self._SCAN_PROD_MAX
        loci = self._design_locus_for(cand, org)
        self._scan_ctx = {"prod_max": prod_max, "loci": loci or [], "designed_size": size}
        self._render_genome_busy(f"Scanning the {org} genome locally (isPcr) for priming sites. Runtime depends "
                                 "on genome size, the installed isPcr build and cached RefSeq assembly, and local "
                                 f"hardware. Only products up to {prod_max:,} bp are considered…")
        p = {"min_perfect": 15, "min_good": 15, "prod_max": prod_max, "prod_min": 0}
        self.engine.submit("genome_pcr", {"fwd": cand["left_seq"], "rev": cand["right_seq"], "organism": org,
                                          "design_locus": loci, "params": p}, key="genome_pcr")

    def _design_locus_for(self, cand, org):
        """Where the DESIGNED product sits in the genome, as the one-element list the engine expects:
        {accession, start, end, strand, intended_product_span, intended_endpoint_tolerance_bp}, the fetched
        window (1-based inclusive) plus the designed product's span inside that record. Sent only for a pair
        whose own sealed design record is a coordinate fetch of the SAME assembly being scanned and is still the
        record the specimen on screen holds; anything else returns None and the scan stays neutral ('genomic
        priming sites'). The region is taken from the pair's sealed design input, so in a multi-region fetch it is
        the region the pair was designed on, never another one."""
        inp = (cand.get("design_provenance") or {}).get("input") or {}
        src = self.state.get("source") or {}
        scanned_asm = (all_assemblies().get(org, {}) or {}).get("assemblyAccession", "")
        regions = inp.get("regions") or []
        if (not scanned_asm or inp.get("inputType") != "fetched_coordinate" or len(regions) != 1 or
                inp.get("assemblyAccession") != scanned_asm or src.get("assemblyAccession") != scanned_asm or
                cand.get("intended_product_span") is None):
            return None
        r = regions[0]
        try:
            i = (src.get("regions") or []).index(r)
            on_screen = (src.get("records") or [])[i].get("sha256") == inp.get("sha256")
        except (ValueError, IndexError, AttributeError):
            return None                                       # the specimen on screen is not the one it was designed on
        if not on_screen or not r.get("chrAccession"):
            return None
        return [{"accession": r["chrAccession"], "start": r.get("start"), "end": r.get("end", r.get("stop")),
                 "strand": r.get("strand"), "intended_product_span": cand["intended_product_span"],
                 "intended_endpoint_tolerance_bp": cand.get("intended_endpoint_tolerance_bp", 0)}]

    def _prepare_genome(self, org, then_scan=None):
        """Download + cache an organism's RefSeq genome for local scanning (one-time, non-blocking)."""
        if self._genome_prep_inflight or self._genome_inflight:   # also block while a scan runs (reachable via the manager)
            return self._banner("A whole-genome scan / download is already running — wait for it to finish.", "info")
        if not org:
            return self._banner("Select an organism to download.", "info")
        self._pending_scan = then_scan
        self._genome_prep_inflight = True
        self._refresh_genome_manager()                        # an open manager must lock its row actions while this runs
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
            msg = "Genome download failed — " + d.get("error", "unknown error")
            level = "warn" if d.get("cache_corrupt") else "error"
            self._render_genome_status(msg, level)
            self._refresh_genome_manager()                    # download settled (failed): an open manager unlocks, so another genome can be deleted to free space
            return self._banner(msg, level)
        mb = (d.get("bytes", 0) or 0) / 1e6
        self._banner(f"Genome ready · {d.get('organism','')} ({d.get('assemblyAccession','')}) · "
                     f"{d.get('n_seqs','?')} sequences · {mb:.0f} MB cached.", level="success")
        ps = self._pending_scan
        if ps and self._scan_is_stale():
            self._pending_scan = None                         # the pair belonged to a specimen that has since changed
            self._render_genome_status(f"{d.get('organism','')} genome downloaded. The scan that was waiting for it "
                                       "belonged to a previous specimen and was not run.", "info")
        elif ps:                                              # resume straight into the SAME scan the download was for
            self._pending_scan = None
            self._scan_genome(ps["cand"], org=ps["org"])
        else:
            self._render_genome_status(f"{d.get('organism','')} genome downloaded — choose a pair above and press Scan whole genome.", "success")
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
            picker.addItem(f"{self._pair_label(c)}  ·  {c['left_seq'][:12]}… / {c['right_seq'][:12]}…", c)
        picker.blockSignals(False)
        has_pairs = bool(cands)
        if hasattr(self, "pcrStageAll"):
            self.pcrStageAll.setEnabled(has_pairs)            # 'stage all' only ever stages the current designs
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
        if not getattr(self, "_fields_decorated", False) and "pcrPmax" in getattr(self, "pfields", {}):
            # cards 04 and 05 exist by the time card 06 builds, and this runs once while it does: the primer and
            # in-silico PCR fields get their units and tooltips without editing the card builders
            self._fields_decorated = True
            self._decorate_primer_fields()
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
                       "appears in the organism menus. Every genome download needs at least 8 GB free disk in WSL, "
                       "however small the genome. Mammalian genomes are large (~1 GB): extraction peaks near 4 GB "
                       "(FASTA ~3 GB + 2bit ~0.8 GB) before the temporary files are removed.")
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
                b.setText(widgets.caps(b.text()))                   # refresh reaches this dialog) — measure the FINAL text or
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
            what = ("genome download" if self._genome_prep_inflight else "whole-genome scan" if self._genome_inflight
                    else "genome annotation")
            b = QLabel(f"Actions are disabled until the {what} finishes.")
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
        intro = QLabel("TEagle runs RepeatMasker over the whole assembly and summarises which repeat names/models "
                       "match, how much of the genome their intervals cover, and their alignment-span-weighted "
                       "divergence. It matches references represented in the selected library; it does not discover "
                       "new repeats, so coverage depends on that library.")
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
            f"nothing at all: on baker's yeast it sees {_DFAM_YEAST_CURATED} of {_DFAM_YEAST_TOTAL} family "
            f"models and reports no transposable element, while including uncurated adds "
            f"{_DFAM_YEAST_ADDITIONAL} models and finds the Ty elements the genome really has.\n\n"
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
                refresh_cost()
                return
            fn, _ = QFileDialog.getOpenFileName(dlg, "Choose a repeat library (FASTA)", "",
                                                "FASTA library (*.fa *.fasta *.lib *.txt);;All files (*)")
            if not fn:
                lsel.setCurrentIndex(0); self._annot_custom_lib = None; libpath.setText(""); return
            self._annot_custom_lib = fn
            libpath.setText(f"{os.path.basename(fn)} — used as the sole repeat library; TEagle records its "
                            "checksum but cannot vouch for its contents.")
            refresh_cost()
        lsel.currentIndexChanged.connect(pick_library)
        cost = QLabel("Estimating…"); cost.setObjectName("orient"); cost.setWordWrap(True)
        cost.setTextInteractionFlags(Qt.TextSelectableByMouse); lay.addWidget(cost)
        warn = QLabel(""); warn.setObjectName("errbanner"); warn.setWordWrap(True); warn.hide(); lay.addWidget(warn)
        start = QPushButton("Start annotation"); start.setEnabled(False)
        retry = QPushButton("Check again"); retry.setProperty("sm", True); retry.hide()
        retry.setToolTip("Read this machine's compute and disk budget again")
        cancel = QPushButton("Close"); cancel.setProperty("sm", True); cancel.clicked.connect(dlg.reject)
        brow = QHBoxLayout(); brow.addStretch(1); brow.addWidget(retry); brow.addWidget(start); brow.addWidget(cancel)
        lay.addLayout(brow)
        self._annot_retry = retry
        self._annot_cost_state = {}

        def refresh_cost():
            acc = gsel.currentData()
            g = next((x for x in prepared if x["accession"] == acc), {})
            org = acc_to_org.get(acc, "")
            custom_path = self._annot_custom_lib
            library_mode = ("custom" if custom_path else
                            ("dfam_uncurated" if self._annot_uncurated else "dfam_curated"))
            cost.setText("Estimating cost — reading the WSL compute and disk budget…" if custom_path else
                         "Estimating cost — reading the WSL budget and the installed Dfam library…")
            start.setEnabled(False); warn.hide(); retry.hide()
            self._annot_budget_generation = int(getattr(self, "_annot_budget_generation", 0)) + 1
            self._annot_cost_state = {"accession": acc, "organism": org, "sha256": g.get("sha256"),
                                      "custom_library": custom_path, "library_mode": library_mode,
                                      "budget_request_id": self._annot_budget_generation}
            body = {"species": org, "genome_bytes": int(g.get("bytes") or 0), "library_mode": library_mode,
                    "budget_request_id": self._annot_budget_generation}
            if g.get("fasta_bytes"):                     # the source FASTA is what the disk estimate is sized from
                body["fasta_bytes"] = int(g["fasta_bytes"])
            self.engine.submit("annotate_budget", body, key="annotate_budget")

        self._annot_cost_widgets = (cost, warn, start, gsel, ssel)
        retry.clicked.connect(lambda *_: refresh_cost())
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
            # Budget jobs run concurrently. A slow Dfam lineage query can finish after the user has
            # switched to a custom FASTA (or another Dfam search set); never render that stale response
            # under the new selection.
            if (b.get("budget_request_id") != st.get("budget_request_id")
                    or b.get("library_mode") != st.get("library_mode")):
                return
            retry = getattr(self, "_annot_retry", None)
            if b.get("ok") is False:
                # the probe itself failed (WSL slow or down): none of the promised cost or disk checks ran, so the
                # run must not be armed, and the dialog says why instead of printing a table of 'None'
                cost.setText("TEagle could not read this machine's compute and disk budget, so the cost of this "
                             "annotation is unknown.")
                warn.setText("Budget check failed — " + str(b.get("error") or "the WSL backend did not answer") +
                             ". Check that WSL is running, then press Check again.")
                warn.show(); start.setEnabled(False)
                if retry is not None:
                    retry.show()
                return
            if retry is not None:
                retry.hide()
            # keep the MEASURED thread count so the run uses the number the dialog promised
            self._annot_threads = int(b.get("recommended_threads") or 4)
            fams = b.get("library_families_for_species")
            gb = b.get("disk_needed_gb")
            custom_path = st.get("custom_library")
            library_mode = st.get("library_mode")
            lines = [
                f"<b>This machine (WSL compute):</b> {b.get('cores')} cores · {b.get('mem_gb')} GB RAM. "
                f"TEagle will use {b.get('recommended_threads')} parallel "
                f"RepeatMasker jobs (limited by {b.get('limited_by')}).",
                f"<b>Storage:</b> {b.get('disk_note') or (str(b.get('avail_gb')) + ' GiB usable')}",
                f"<b>Disk needed:</b> about {gb} GB — the genome sequence is fetched once and kept, and the "
                "work chunks are removed as they finish.",
                "<b>Time:</b> a fly-sized genome (~140 Mb) takes minutes to about an hour here; a mammalian "
                "genome (~3 Gb) takes many hours to a day. The run continues in the background and can be "
                "re-started later — finished chunks are not repeated.",
            ]
            if custom_path:
                lines.append(f"<b>Selected library:</b> user-supplied FASTA {html.escape(os.path.basename(custom_path))}. "
                             "This file supplies the interspersed-repeat homology library, so only references "
                             "represented in it can receive homology-based TE calls. RepeatMasker’s separate "
                             "simple- and low-complexity screen can still produce non-TE repeat rows. To enter "
                             "the TE percentage, a custom reference also needs a recognized class/family in "
                             "its header, for example >name#LTR/Gypsy; unclassified matches remain visible but "
                             "are not silently counted as TEs.")
            else:
                selected_set = ("curated + uncurated" if library_mode == "dfam_uncurated"
                                else "curated-only")
                lines.append(f"<b>Library coverage:</b> the selected {selected_set} Dfam search set contains "
                             f"<b>{fams if fams is not None else 'an unknown number of'}</b> family models for "
                             f"{st.get('organism') or 'this lineage'}.")
            cost.setText("<br>".join(lines))
            if library_mode == "dfam_uncurated" and b.get("library_ready") is False:
                warn.setText(b.get("library_error") or
                             "The selected curated + uncurated Dfam search cannot be run as configured.")
                warn.show(); start.setEnabled(False); return
            if b.get("disk_ok") is False:
                warn.setText(f"Not enough usable disk — about {gb} GB is needed after applying the smaller "
                             "WSL filesystem and Windows backing-volume limit. Free space, or remove a cached "
                             "genome from this manager, then try again.")
                warn.show(); start.setEnabled(False); return
            # Threshold rationale: measured on this backend, a lineage Dfam covers usefully has hundreds
            # of models (D. melanogaster 399, H. sapiens 1439) while an uncovered one has single digits
            # (S. cerevisiae 9, and its run found zero TEs). 25 sits in the empty gap between those two
            # regimes; it is a warning trigger only — it never blocks a run or changes a result.
            if not custom_path and fams is not None and fams < self._ANNOT_LOW_LIBRARY:
                lineage = st.get("organism") or "this lineage"
                if library_mode == "dfam_uncurated":
                    warn.setText(f"The selected curated + uncurated Dfam search set contains only {fams} "
                                 f"family model(s) for {lineage}. Both family sets are already being searched, "
                                 "so a blank result would still be limited by this library's lineage coverage. "
                                 "Use a suitable organism-specific FASTA library before treating absence as "
                                 "evidence that the genome lacks transposable elements.")
                else:
                    warn.setText(f"The selected curated-only Dfam search contains only {fams} family model(s) "
                                 f"for {lineage}; uncurated families will not be searched. A blank result would "
                                 "therefore be inconclusive. Select “include uncurated families” (and install "
                                 "those partitions if needed), or use a suitable organism-specific FASTA library, "
                                 "before interpreting absence.")
                warn.show()
            start.setEnabled(True)
        except RuntimeError:
            pass

    def _start_genome_annotate(self, acc, organism, sensitivity, custom_lib=None, include_uncurated=False, sha256=None):
        if self._annot_inflight or not acc:
            return
        self._annot_inflight = True
        sha256 = sha256 or (getattr(self, "_annot_cost_state", None) or {}).get("sha256")
        # kept so 'Discard and start over' can start the run the user chose, not just discard the old one
        self._annot_request = {"acc": acc, "organism": organism, "sensitivity": sensitivity, "custom_lib": custom_lib,
                               "include_uncurated": bool(include_uncurated), "sha256": sha256}
        self._banner(f"Annotating {organism or acc} — this runs in the background. Its progress is shown at the top of "
                     "panel 06; finished chunks are kept if the run is interrupted.", "info")
        self._refresh_genome_manager()
        self._set_status_line(f"Annotating {organism or acc} — started; reading the first chunk. A mammalian genome "
                              "takes many hours. There is no in-app cancel yet.")
        self.card_genome.expand()                              # once, so the line is seen; later updates never re-open a card the user closed
        self.engine.submit("genome_annotate",
                           {"assemblyAccession": acc, "species": organism, "sensitivity": sensitivity,
                            # the thread count the budget probe MEASURED for this machine — the cost dialog
                            # showed it, so the run must actually use it rather than a fixed default
                            "threads": int(getattr(self, "_annot_threads", 0) or 4),
                            "chunk_mb": int(getattr(self, "_annot_chunk_mb", 0) or 40),
                            "custom_library": custom_lib, "include_uncurated": bool(include_uncurated),
                            "sha256": sha256},
                           key="genome_annotate")
        self._annot_timer = QTimer(self); self._annot_timer.setInterval(6000)
        self._annot_timer.timeout.connect(lambda: self.engine.submit("genome_annotate_log", {"tail": 1},
                                                                     key="genome_annotate_log"))
        self._annot_timer.start()
        # A fault ends the run through _reset_buttons, which does not know this bar: its own watcher hides the bar
        # once the run is no longer in flight, so a dead run never keeps reporting progress.
        self._annot_watch = QTimer(self); self._annot_watch.setInterval(2000)
        self._annot_watch.timeout.connect(self._annot_watch_tick)
        self._annot_watch.start()

    def _annot_watch_tick(self):
        if not self._annot_inflight:
            self._hide_annot_status()

    def _hide_annot_status(self):
        w = getattr(self, "_annot_watch", None)
        if w is not None:
            w.stop()
        bar = getattr(self, "_annotStatus", None)
        if bar is not None:
            try:
                bar.hide()
            except RuntimeError:
                pass

    def _on_genome_annotate_log(self, d):
        line = (d or {}).get("log") or ""
        if line and self._annot_inflight:
            org = (getattr(self, "_annot_request", None) or {}).get("organism") or "genome"
            self._set_status_line(f"Annotating {org} — {line}")

    def _set_status_line(self, text):
        """The persistent progress line of a running whole-genome annotation. It lives at the top of card 06, outside
        the panels that are cleared and rebuilt, so it is on screen for the hours the run takes."""
        bar = getattr(self, "_annotStatus", None)
        try:
            if bar is None:
                bar = self._annotStatus = BusyBar(text)
                self.card_genome.bodylay.insertWidget(0, bar)
            bar.set_text(text); bar.show()
        except RuntimeError:
            self._annotStatus = None

    def _on_genome_annotate(self, r):
        self._annot_inflight = False
        self._hide_annot_status()
        t = getattr(self, "_annot_timer", None)
        if t is not None:
            t.stop()
        # The backend reports a refused or failed run as {"ok": False, ...} — a normal RETURN, not an
        # exception — so it arrives here on the success signal. Without this guard the failure was
        # announced as a finished annotation and then crashed the results window on the missing fields.
        if not r or not r.get("ok"):
            msg = (r or {}).get("error") or "the annotation did not complete"
            self._refresh_genome_manager()
            if (r or {}).get("cache_corrupt"):
                self._banner("Genome annotation did not finish — " + msg, "warn")
                return
            # A settings clash is the one failure the user can clear, so offer the action instead of
            # telling them to do something the window gives them no way to do.
            if "different settings" in msg:
                req = getattr(self, "_annot_request", None) or {}
                acc = req.get("acc") or (getattr(self, "_annot_cost_state", None) or {}).get("accession")
                box = QMessageBox(self)
                box.setWindowTitle("Previous annotation found")
                box.setIcon(QMessageBox.Question)
                if (r or {}).get("annotation_state") == "complete":
                    box.setText("This genome already has a finished annotation that used different settings.")
                    box.setInformativeText(
                        "That finished result was made with the earlier settings and cannot stand in for the new ones. "
                        "Discarding it deletes the earlier result from this computer, and its per-hit alignment rows "
                        "cannot be exported first. Discard it and run again with the settings you just chose?")
                else:
                    box.setText("This genome has a part-finished annotation that used different settings.")
                    box.setInformativeText(
                        "Chunks already finished were searched with the earlier settings, so they cannot be "
                        "combined with new ones. Discard that unfinished work and start again with the "
                        "settings you just chose?")
                discard = box.addButton("Discard and start over", QMessageBox.AcceptRole)
                box.addButton("Keep it", QMessageBox.RejectRole)
                box.exec()
                if box.clickedButton() is discard and acc:
                    self._annot_restart = dict(req) if req else None   # the reset reply starts THIS run
                    self.engine.submit("genome_annotate_reset", {"assemblyAccession": acc},
                                       key="genome_annotate_reset")
                return
            self._banner("Genome annotation did not finish — " + msg, "warn")
            return
        # SUCCESS. This block used to live in _on_genome_annotate_reset, so the app's most expensive
        # operation — a multi-hour, uncancellable whole-genome run — finished with no banner and no results
        # window, while DISCARDING an annotation announced "transposable elements cover None% of the
        # assembly across None repeat models" and opened a results window on the reset reply's empty fields.
        # Each handler now reports its own outcome.
        self._annot_result = r
        self._refresh_genome_manager()
        if r.get("coverage_warning"):
            self._banner(r["coverage_warning"], "warn")
        else:
            self._banner(f"Genome annotation finished — transposable elements cover "
                         f"{r.get('te_percent')}% of the assembly across "
                         f"{r.get('te_repeat_name_count')} distinct repeat names/models.",
                         "success")
        self._show_annot_window(r)

    def _on_genome_annotate_reset(self, r):
        redo, self._annot_restart = getattr(self, "_annot_restart", None), None
        if r.get("ok"):
            self._refresh_genome_manager()
            if redo:                                      # 'Discard and start over' means exactly that
                self._banner("Previous annotation discarded. Starting the annotation you chose.", "info")
                return self._start_genome_annotate(redo["acc"], redo["organism"], redo["sensitivity"],
                                                   redo.get("custom_lib"), redo.get("include_uncurated", False),
                                                   redo.get("sha256"))
            self._banner("Previous annotation discarded. Open “Annotate TE landscape” to start again.",
                         "success")
        else:
            self._banner(r.get("error") or "Could not clear the previous annotation.", "warn")
        self._refresh_genome_manager()

    _ANNOT_COLS = [("Repeat name/model", "The matched RepeatMasker repeat name or consensus model (output column 10); "
                                        "this is not necessarily a biological family"),
                   ("Class/family", "RepeatMasker's broader class/family label. TEagle uses this label to "
                                    "separate TE, tandem and other coverage."),
                   ("Kind", "TE = transposable element; tandem = simple/low-complexity/satellite; other = "
                            "rRNA/tRNA; unclassified/other = no recognized class was available"),
                   ("Copies", "RepeatMasker alignment rows; one fragmented insertion can produce several rows"),
                   ("Bases", "Bases in this repeat name/model's within-name interval union. Different rows may "
                             "overlap, so rows must not be summed to recover a headline total."),
                   ("% genome", "Share of all submitted FASTA bases, including ambiguous bases"),
                   ("Mean divergence %", "Alignment-span-weighted mean RepeatMasker divergence from consensus. "
                                          "Raw alignment rows, including overlaps and fragments, supply the "
                                          "weights. It is not Kimura- or CpG-corrected, a mutation rate, or an age.")]

    @staticmethod
    def _annot_library_context(r):
        """Library-specific copy shared by the landscape window and its Markdown export."""
        kind = str(r.get("library_kind") or "")
        raw_custom = r.get("custom_library")
        custom = raw_custom if isinstance(raw_custom, dict) else {}
        is_custom = bool(custom) or kind.strip().lower() == "user" or kind.lower().startswith("user-supplied")
        rm_version = html.escape(str(r.get("repeatmasker_version") or "unresolved"))

        if is_custom:
            name = str(custom.get("name") or "unnamed user FASTA")
            sha = str(custom.get("sha256") or "unavailable")
            sequences = custom.get("sequences")
            sequence_md = f"; {sequences} reference sequences" if sequences is not None else ""
            sequence_ui = f"; {html.escape(str(sequences))} reference sequences" if sequences is not None else ""
            md_name = (name.replace("\\", "\\\\").replace("|", "\\|").replace("*", "\\*")
                       .replace("_", "\\_").replace("`", "\\`"))
            safe_name = html.escape(name)
            safe_sha = html.escape(sha)
            return {
                "report": (f"- Searched with RepeatMasker {rm_version} against the user-supplied FASTA "
                           f"{md_name}; SHA-256 `{sha}`{sequence_md}"),
                "display": (f"Searched with RepeatMasker {rm_version} against user-supplied FASTA "
                            f"<b>{safe_name}</b> (SHA-256 <code>{safe_sha}</code>{sequence_ui})"),
                "tooltip": (f"The run searched the user-supplied FASTA {safe_name}, identified by SHA-256 "
                            f"{safe_sha}. The reference-sequence count describes that file; it is not a count "
                             "of detected repeat names/models."),
                "scope": ("What this is: interspersed-repeat matches to references represented in the "
                          "supplied FASTA, placed by homology. It is not de novo repeat discovery. A reference "
                          "absent from that FASTA cannot be recovered through library homology, so absence "
                          "is not evidence of absence. RepeatMasker's separate simple/low-complexity screen can "
                          "still contribute non-TE rows. A custom reference contributes to the TE percentage only "
                          "when its FASTA header supplies a recognized class/family after '#', for example "
                          ">name#LTR/Gypsy; otherwise its matches remain visible as unclassified/other and are "
                          "excluded because classification is unavailable, not because they were shown to be non-TE. "
                          "Percentages use all submitted FASTA bases, including ambiguous bases and unplaced scaffolds."),
                "report_scope": ("- Interspersed-repeat matches to references represented in the supplied "
                                 "FASTA, placed by homology. This is not de novo repeat discovery. A reference "
                                 "absent from that FASTA cannot be recovered through library homology; "
                                 "RepeatMasker's separate simple/low-complexity screen can still contribute non-TE "
                                 "rows. Custom references without a recognized `#class/family` header remain "
                                 "visible as unclassified/other and do not enter the TE percentage because their "
                                 "classification is unavailable, not because they were demonstrated to be non-TE."),
            }

        library_kind = kind or "installed Dfam partitions"
        dfam_version = r.get("dfam_version")
        models = r.get("library_families_for_species")
        version_md = f" (Dfam {dfam_version})" if dfam_version else ""
        models_md = f"; {models} family models available for this lineage" if models is not None else ""
        safe_kind = html.escape(library_kind)
        version_ui = f" (Dfam {html.escape(str(dfam_version))})" if dfam_version else ""
        models_ui = f"; {html.escape(str(models))} family models available for this lineage" if models is not None else ""
        return {
            "report": (f"- Searched with RepeatMasker {rm_version} against {library_kind}{version_md}{models_md}"),
            "display": (f"Searched with RepeatMasker {rm_version} against {safe_kind}{version_ui}{models_ui}"),
            "tooltip": ("The exact tool and installed Dfam family set this result came from. The lineage model "
                        "count describes profiles available to the search, not families detected in this genome. "
                        "Chunks are assembly pieces used for progress and resumability."),
            "scope": ("What this is: matches to repeat names/models represented in the installed Dfam library, "
                      "placed by homology. It is not de novo repeat discovery. A model absent from the library "
                      "cannot appear here, so absence is not evidence of absence. Percentages use all submitted "
                      "FASTA bases, including ambiguous bases and unplaced scaffolds."),
            "report_scope": ("- Matches to repeat names/models represented in the searched Dfam library, placed "
                             "by homology. This is not de novo repeat discovery. A model absent from the library "
                             "cannot appear here; absence is not evidence of absence."),
        }

    @staticmethod
    def _annot_kind_label(row):
        """An unknown class is unavailable evidence, not a demonstrated non-TE classification."""
        kind = str(row.get("kind") or "")
        class_family = str(row.get("class_family") or "").strip().lower()
        if kind == "other" and class_family in ("", "unknown", "unclassified", "na", "n/a"):
            return "unclassified/other"
        return kind

    @staticmethod
    def _annot_report_md(r):
        """The landscape as a self-contained Markdown report.

        Every hedge the window shows travels with it: what was searched, what the two percentages mean,
        and why an empty TE result may be a library limit. An export that carried only the table would
        read as an unqualified statement about the genome."""
        repeat_names = r.get("repeat_names") or []
        te = [row for row in repeat_names if row.get("kind") == "TE"]
        library = MainWindow._annot_library_context(r)
        L = [f"# Transposable-element landscape — {r.get('species')} ({r.get('accession')})", "",
             f"- **Transposable elements: {r.get('te_percent')}%** of "
             f"{(r.get('genome_bp') or 0)/1e6:.1f} Mb, across "
             f"{r.get('te_repeat_name_count')} distinct repeat names/models",
             f"- All repeat classes together: {r.get('masked_percent')}% "
             f"({r.get('total_hits'):,} alignment rows)",
             library["report"],
             f"- Sensitivity {r.get('sensitivity')} · {r.get('chunks')} chunks · "
             f"complete: {bool(r.get('complete'))}", ""]
        if not r.get("complete"):
            L += ["> **Incomplete run.** Counts, percentages, and rows are partial lower bounds. No assembly-wide "
                  "absence, coverage, or exhaustive-family claim is supported.", ""]
        if r.get("coverage_warning"):
            L += ["> **Coverage warning.** " + r["coverage_warning"], ""]
        L += ["## Scope and limits", "",
              library["report_scope"],
              "- Transposable elements are counted separately from tandem repeats (simple repeats,",
              "  low-complexity sequence, satellites) and non-TE entries; only the first is a statement",
              "  about transposable-element content.",
              "- Each repeat name/model's Bases value is its within-name interval union. The all-repeat and TE headline "
              "percentages use separate cross-model unions, so a genomic base contributes at most once to each "
              "headline. Summing rows can exceed a headline because repeat models can overlap.",
              "- TE, tandem and other category percentages are independently unioned and non-additive.",
              "- Divergence is the alignment-span-weighted mean RepeatMasker divergence from consensus. Raw alignment",
              "  rows, including overlaps and fragments, supply the weights. It is not Kimura- or CpG-corrected,",
              "  a mutation rate, or an age.",
              "- Percentages use all submitted FASTA bases, including ambiguous bases and unplaced scaffolds.", "",
              "## Repeat names/models", "",
              "| Repeat name/model | Class/family | Kind | Copies (alignment rows) | Bases (within-name union) | % genome | Mean divergence % |",
              "|---|---|---|---|---|---|---|"]
        for row in repeat_names:
            L.append(f"| {row.get('repeat_name')} | {row.get('class_family')} | {MainWindow._annot_kind_label(row)} | {row.get('n'):,} | {row.get('bp'):,} | "
                     f"{row.get('percent')} | {row.get('divergence')} |")
        L += ["", f"_{len(te)} of {len(repeat_names)} rows have a recognized transposable-element class._"]
        return "\n".join(L)

    def _annot_row_menu(self, row):
        """Context menu for one genome-landscape row — aggregate-appropriate actions only."""
        kind = row.get("kind")
        items = [(f"Copy repeat name/model ({row.get('repeat_name')})", lambda: self._copy(str(row.get("repeat_name")))),
                 ("Copy this row (TSV)", lambda: self._copy("\t".join(
                     str(x) for x in (row.get("repeat_name"), row.get("class_family"), self._annot_kind_label(row), row.get("n"), row.get("bp"),
                                      row.get("percent"), row.get("divergence")))))]
        if kind != "TE":
            if self._annot_kind_label(row) == "unclassified/other":
                reason = ("This row has no recognized RepeatMasker class/family label. It is excluded from "
                          "TE coverage because its classification is unavailable, not because it was demonstrated "
                          "to be non-TE. For a custom library, use headers such as >name#LTR/Gypsy when supported.")
            else:
                reason = ("Simple repeats, low-complexity sequence and satellites are tandem repeats; rRNA and "
                          "tRNA entries are other repeat annotations. TEagle reports these categories separately "
                          "from transposable-element coverage.")
            items.append(("Why is this not counted as a TE?",
                          lambda _reason=reason: self._banner(_reason, "info")))
        else:
            # Same field as the "Mean divergence %" column, so it must carry the same hedges. This text
            # previously said "higher means older" and dropped the raw/uncorrected and saturation caveats
            # the column tooltip states — two different claims about one number, and the looser one was the
            # one a user reached by asking for an explanation.
            items.append(("What does mean divergence mean?",
                          lambda: self._banner(
                              "Alignment-span-weighted mean RepeatMasker divergence from the consensus. Raw "
                              "alignment rows, including overlapping and fragmented rows, supply the weights. "
                              "It is not Kimura- or CpG-corrected, a mutation rate, or an insertion-age estimate.",
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
        accession_stem = _safe_export_stem(r.get("accession"), "assembly", trust="validated_fetch")
        head = QLabel(f"{r.get('species')} · {r.get('accession')}")
        hf = head.font(); hf.setBold(True); hf.setPointSizeF(hf.pointSizeF() + 2); head.setFont(hf); lay.addWidget(head)
        # headline numbers: the TE share is stated separately from every repeat, because they are different
        # claims and merging them would overstate transposable-element content.
        tot = QLabel(f"<b>Transposable elements: {r.get('te_percent')}%</b> of "
                     f"{(r.get('genome_bp') or 0)/1e6:.1f} Mb, across "
                     f"{r.get('te_repeat_name_count')} distinct repeat names/models &nbsp;·&nbsp; "
                     f"all repeat classes together: {r.get('masked_percent')}% &nbsp;·&nbsp; "
                     f"{r.get('total_hits'):,} alignment rows")
        tot.setTextFormat(Qt.RichText); tot.setWordWrap(True)
        tot.setToolTip("Transposable elements = the share of the assembly covered by TE-classified matches (LTR, LINE, "
                       "SINE, DNA transposons, rolling-circle, retroposons, Penelope-like elements).\n"
                       "All repeat classes = the separate contig-wise union of every repeat row, including "
                       "tandem repeats (simple repeats, low-complexity sequence, satellites) and non-TE "
                       "entries such as rRNA and tRNA genes.\n"
                       "The categories are independently unioned and non-additive. Only the TE union is a statement "
                       "about transposable-element content. Each headline uses a separate cross-model interval union, "
                       "so a genomic base contributes at most once to that percentage; repeat-name/model rows use "
                       "within-name/model unions and therefore must not be summed. Percentages use all submitted "
                       "FASTA bases, including ambiguous bases.")
        lay.addWidget(tot)
        library = self._annot_library_context(r)
        prov = QLabel(f"{library['display']} · sensitivity {r.get('sensitivity')} · "
                      f"{r.get('chunks')} chunks · {r.get('elapsed_s')} s")
        prov.setObjectName("orient"); prov.setWordWrap(True)
        prov.setToolTip(library["tooltip"])
        lay.addWidget(prov)
        scope = QLabel(library["scope"])
        scope.setObjectName("orient"); scope.setWordWrap(True); lay.addWidget(scope)
        if r.get("coverage_warning"):
            cw = QLabel(r["coverage_warning"]); cw.setObjectName("errbanner"); cw.setWordWrap(True); lay.addWidget(cw)
        complete_note = ("The run completed all submitted assembly chunks. Library scope and sensitivity still "
                         "limit absence claims." if r.get("complete") else
                         "The run is incomplete: all counts, percentages, and rows are partial lower bounds. "
                         "No assembly-wide absence, coverage, or exhaustive-family claim is supported.")
        completion = QLabel(complete_note)
        completion.setObjectName("orient" if r.get("complete") else "errbanner")
        completion.setWordWrap(True); lay.addWidget(completion)
        repeat_names = r.get("repeat_names") or []
        t = DataTable([c[0] for c in self._ANNOT_COLS], dict(self._ANNOT_COLS))
        t.set_rows([[row["repeat_name"], row.get("class_family"), self._annot_kind_label(row), row["n"], row["bp"],
                     row.get("percent"), row["divergence"]] for row in repeat_names])
        t.export_notes = [library["scope"], complete_note,
                          "TE, tandem, other, and unclassified categories are independently unioned and non-additive. "
                          "Repeat-name/model rows use within-name unions and must not be summed to recover a headline. "
                          "Divergence is an alignment-span-weighted raw percentage, not an age estimate."]
        if r.get("coverage_warning"):
            t.export_notes.append(str(r["coverage_warning"]))
        # Row menu matched to what the row IS. A landscape row is a per-name/model AGGREGATE, not a locus, so
        # it offers no "design primer here" or "copy sequence" — there is no single sequence behind it.
        # A tandem/other row additionally says why it is not counted as a transposable element.
        t.set_row_menu(lambda r, _rows=repeat_names: self._annot_row_menu(_rows[r]) if r < len(_rows) else [])
        t.setMinimumHeight(round(260 * theme_mod.UI_SCALE)); lay.addWidget(t)
        # FigurePanel re-renders on theme/zoom, so it takes a builder (bg -> svg), not a rendered string.
        fig = FigurePanel(lambda bg, _r=r: figures.svg_te_composition(_r, theme=bg),
                          f"TEagle_TE_landscape_{accession_stem}")
        fig.setMinimumHeight(round(190 * theme_mod.UI_SCALE)); lay.addWidget(fig)
        row = QHBoxLayout(); row.addStretch(1)
        row.addWidget(_export_table_btn(t, f"TEagle_TE_landscape_{accession_stem}", self))
        # A genome run is sealed like every other result, so the seal has to be reachable from the result
        # — otherwise the run is reproducible in principle and not in practice. The report carries the
        # same hedges the screen shows, so an exported file cannot read as an unqualified claim.
        rep = QPushButton("Export report (.md)"); rep.setProperty("sm", True)
        rep.setToolTip("The summary table plus what was searched, what the numbers mean, and the limits "
                       "stated on this screen — the same hedges, in the file.")
        rep.clicked.connect(lambda: self._save_text(self._annot_report_md(r),
                                                    f"TEagle_TE_landscape_{accession_stem}", "md"))
        row.addWidget(rep)
        manb = QPushButton("Export manifest (.json)"); manb.setProperty("sm", True)
        manb.setToolTip("The sealed provenance record: tool and library versions, every parameter, the "
                        "genome checksum, and whether the run covered the whole assembly.")
        manb.setEnabled(bool(r.get("provenance")))
        manb.clicked.connect(lambda: self._export_manifest_payload(
            [r.get("provenance")], f"TEagle_TE_landscape_{accession_stem}_manifest"))
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
        if self._genome_inflight or self._genome_prep_inflight or self._annot_inflight:
            return self._banner("A genome scan, download or annotation is running — wait for it to finish before "
                                "deleting a genome.", "warn")
        name = org or accession
        g = next((x for x in (getattr(self, "_prepared_genomes", None) or []) if x.get("accession") == accession), {})
        size = f"The cached search copy is {(g['bytes'] or 0) / 1e6:.0f} MB. " if g.get("bytes") else ""
        box = QMessageBox(self)
        box.setWindowTitle("Delete genome?")
        box.setText(f"Delete the cached {name} genome ({accession}) from this machine?\n\n{size}"
                    "Deleting also removes any whole-genome annotation saved with it, including every RepeatMasker chunk "
                    "already finished, and the genome FASTA an annotation cached, if there is one. An annotation would "
                    "then start again from the beginning.\n\n"
                    "Scanning it again needs the genome downloaded again from NCBI: several minutes, and at least 8 GB "
                    "free disk in WSL for any genome.")
        box.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
        box.setDefaultButton(QMessageBox.No)
        if box.exec() != QMessageBox.Yes:
            return
        self.engine.submit("genome_remove", {"assemblyAccession": accession}, key="genome_remove")

    def _on_genome_remove(self, d):
        # A refusal has to be visible. genome_remove can now decline (an annotation still holds the genome),
        # and this handler previously discarded its reply entirely — so a declined delete looked exactly like
        # a successful one: the list refreshed and the genome was simply still there.
        if isinstance(d, dict) and not d.get("ok"):
            self._banner(d.get("error") or "The genome could not be deleted.", level="warn")
        # refresh BOTH the dropdown and the manager (if open); the split _on_genome_list self-selects what to update
        self.engine.submit("genome_list", {}, key="genome_list")

    def _reset_genome_scan(self, why):
        """Drop a rendered whole-genome scan once the primers it describes are gone.

        The panel is only ever cleared by the scan paths themselves, so a verdict like "no off-target site
        in this assembly" survived a primer redesign and a new specimen — sitting under panel 06 as though
        it described whatever is on screen now. A stale specificity verdict is the worst kind of wrong
        answer here, because it reads as a clean result. No-op when nothing has been scanned.

        `_invalidate_analysis` drops the stored scan BEFORE it calls this, so the stored scan alone cannot say
        whether a verdict is still on screen; the panel's own flag does."""
        self._refresh_scan_picker()                       # the pairs on offer follow state['candidates'], which a new specimen drops
        if not (self.state.get("lastGenomeScan") or getattr(self, "_genome_shows_result", False)):
            return
        self.state["lastGenomeScan"] = None
        if getattr(self, "_genome_inflight", False) or getattr(self, "_genome_prep_inflight", False):
            return                                        # a running scan owns the panel; it will render its own result
        _clear_layout(self.genomeBody)
        self._genome_shows_result = False
        self.genomeBody.addWidget(_empty(why))

    def _render_genome_status(self, msg, level="error"):
        # every caller reports an OUTCOME (failed / cancelled / downloaded / not-yet-downloaded), never an
        # empty panel — so it gets the errbanner vocabulary at the caller's level, not the grey empty style
        _clear_layout(self.genomeBody)
        self._genome_shows_result = False
        self.genomeBody.addWidget(_note(msg, level))

    def _render_genome_busy(self, text):
        """Long-op liveness in the whole-genome scan panel: an animated indeterminate bar (download / scan),
        so a multi-minute WSL call visibly reads as working, not hung. set_text() updates it from the log poll."""
        _clear_layout(self.genomeBody)
        self._genome_shows_result = False
        self._genome_busy = BusyBar(text)
        self.genomeBody.addWidget(self._genome_busy)

    def _on_genome_pcr(self, d):
        if self._scan_is_stale():
            return self._discard_stale_scan()
        self._genome_inflight = False
        self._refresh_genome_manager()                        # scan settled — re-enable a manager opened mid-scan
        if not d.get("ok"):
            if d.get("cache_corrupt"):
                ps = self._pending_scan or {}
                org = ps.get("org") or self.genomeOrg.currentData()
                acc = (all_assemblies().get(org, {}) or {}).get("assemblyAccession", "")
                self._pending_scan = None
                msg = (f"The cached {org or 'selected'} genome ({acc or 'unknown accession'}) failed its "
                       "integrity check. Open Manage genomes, remove that cached genome, then download it "
                       "again before scanning. TEagle did not run isPcr on the altered cache.")
                self._render_genome_status(msg, "warn")
                return self._banner(msg, "warn")
            if d.get("need_prepare"):                         # genome not downloaded yet — offer the one-time download
                ps = self._pending_scan or {}
                org = ps.get("org") or self.genomeOrg.currentData()   # the org that STARTED this scan, not the live dropdown
                acc = (all_assemblies().get(org, {}) or {}).get("assemblyAccession", "")
                self._render_genome_status(f"The {org} genome ({acc}) is not downloaded yet.", "info")
                box = QMessageBox(self)
                box.setWindowTitle("Download genome?")
                box.setText(f"The {org} genome ({acc}) is not on this machine yet.\n\n"
                            "Download it once now? It is kept locally for future scans. Every genome download needs at "
                            "least 8 GB free disk in WSL, however small the genome. Mammalian genomes are large (~1 GB "
                            "download, a few minutes) and extraction peaks near 4 GB (FASTA ~3 GB + 2bit ~0.8 GB) "
                            "before the temporary files are removed.")
                box.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
                box.setDefaultButton(QMessageBox.Yes)
                if org and box.exec() == QMessageBox.Yes:
                    self._prepare_genome(org, then_scan=ps)
                else:
                    self._pending_scan = None
                    self._render_genome_status("Genome scan cancelled — no genome downloaded.", "warn")
                return
            self._render_genome_status("No genome scan result — " + d.get("error", "scan failed"))
            return self._banner(d.get("error", "The genome scan failed."))
        # WHICH pair this verdict describes, captured before _pending_scan is dropped. A scan result says
        # things like "no off-target site in this assembly"; without the pair on the panel that sentence
        # reads as a property of the specimen rather than of two specific oligos.
        _ps = self._pending_scan or {}
        _cand = _ps.get("cand") or {}
        self._pending_scan = None
        amps = [{**a, "pair": "genome", "remote": True} for a in d.get("amplicons", [])]
        summary = d.get("summary") or {}
        has_locus = bool(summary.get("has_locus"))        # gel: on/off colours with a locus; neutral 'priming site' colour without
        lanes = [{"label": "genome", "amplicons": d.get("amplicons", []), "has_locus": has_locus, "advisory": not has_locus}]
        ctx = getattr(self, "_scan_ctx", None) or {}
        prm = (d.get("provenance") or {}).get("parameters") or {}
        self.state["lastGenomeScan"] = {"lanes": lanes, "amplicons": amps, "summary": summary,
                                        "pair": {"id": self._pair_label(_cand) if _cand else "", "fwd": _cand.get("left_seq", ""),
                                                 "rev": _cand.get("right_seq", "")},
                                        "organism": _ps.get("org", ""),
                                        "assembly": {"accession": d.get("assemblyAccession", ""), "name": d.get("assemblyName", "")},
                                        "criteria": {"isPcr": prm.get("isPcr_version"), "min_perfect": prm.get("min_perfect"),
                                                     "min_good": prm.get("min_good"),
                                                     "max_size": prm.get("max_size") or ctx.get("prod_max") or self._SCAN_PROD_MAX},
                                        "loci": ctx.get("loci") or []}
        self.card_genome.expand()
        self._render_genome_scan(lanes, amps, summary)
        m = d.get("provenance", {})
        db = (m.get("databases") or [{}])[0]
        seal = m.get("manifestSha256", "")
        prov = (f" · assembly {d.get('assemblyAccession','')}"
                + (f" · sha256 {db.get('sha256','')[:12]}…" if db.get("sha256") else "")
                + (f" · seal {seal[:12]}…" if seal else ""))
        self._banner(f"Whole-genome scan · {d.get('organism','')} ({d.get('assemblyName','')}): "
                     f"{summary.get('verdict', 'scan complete')} — under isPcr's ≥15 bp minimum perfect "
                     f"3′-end match requirement (sites outside that rule are not reported){prov}.", level="info")
        if m:
            self._render_provenance(m)

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

    def _render_genome_scan(self, lanes, amps, summary):
        """Render a whole-genome off-target scan into its own panel: gel + interpretation (verdict,
        per-assembly-sequence spread, size cluster) + a FULL match table, ON-TARGET FIRST. On-target = the product
        at the specimen's own genome locus (when it sits in the scanned assembly); the rest are off-target candidate priming sites,
        or — with no design locus — neutral genomic priming sites."""
        _clear_layout(self.genomeBody)
        self._genome_shows_result = True                      # a verdict is now on screen (see _reset_genome_scan)
        _last = self.state.get("lastGenomeScan") or {}
        gel_base = "_".join(["TEagle_genome_gel"] + [self._pair_token(x) for x in (
            (_last.get("pair") or {}).get("id"), (_last.get("assembly") or {}).get("accession")) if x])
        gel = FigurePanel(lambda bg, L=lanes: figures.svg_gel({"lanes": L}, bg),
                          gel_base, modes=("dark", "white", "uv", "mono"),
                          hit_regions=figures.gel_regions({"lanes": lanes}), on_menu=self._gel_menu)
        # Name the pair and the assembly the verdict belongs to, ABOVE the verdict. Without it the panel
        # states a specificity conclusion with no subject, and a user who redesigns primers can read the
        # old sentence as describing the new ones.
        _last = self.state.get("lastGenomeScan") or {}
        _pair, _org = (_last.get("pair") or {}), _last.get("organism", "")
        if _pair.get("fwd"):
            _pid = _pair.get("id") or "pair"
            _sub = QLabel(f"Scanned pair <b>{_pid}</b>"
                          + (f" against <b>{_org}</b>" if _org else "")
                          + f"<br>F 5′-{_pair['fwd']}-3′<br>R 5′-{_pair.get('rev','')}-3′")
            _sub.setObjectName("orient"); _sub.setWordWrap(True); _sub.setTextFormat(Qt.RichText)
            self.genomeBody.addWidget(_sub)
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
        crit, asm, loci = _last.get("criteria") or {}, _last.get("assembly") or {}, _last.get("loci") or []
        cap = crit.get("max_size") or self._SCAN_PROD_MAX
        seqs = {a.get("source") for a in amps if not a.get("single_primer")}
        if has_locus and loci and seqs and not any(a.get("on_target") for a in amps) and \
                not seqs & {L.get("accession") for L in loci}:
            # the fetched locus is named on one accession (RefSeq) while the scanned FASTA may name the same
            # chromosome differently (a GenBank GCA assembly), so no product can ever match it by name
            self.genomeBody.addWidget(_note(
                f"The fetched locus lies on {loci[0].get('accession')}, but no product was found on a sequence of "
                "that name. This assembly may name its chromosomes differently (for example a GenBank GCA assembly "
                "against a RefSeq locus), so the intended product could not be matched and the on-target / "
                "off-target split below is not assessed. Treat every product as an unassessed priming site.", "warn"))
        def _call(a):
            return _amp_call(a, has_locus)              # same label the local-PCR table and the gel use
        order = sorted(range(len(amps)), key=lambda i: (amps[i].get("single_primer", False),
                                                        not amps[i].get("on_target", False),
                                                        amps[i]["source"], amps[i]["start"]))
        samps = [amps[i] for i in order]
        scan_table = None
        shown = min(len(samps), self._SCAN_ROWS_SHOWN)
        if samps:
            # isPcr coordinates are 1-based inclusive and are kept verbatim — the headers say so, so Len (end-start+1)
            # does not read as contradicting the span. 'Source' here is the assembly's own sequence name, not a specimen.
            # Start and End are separate numeric columns so an exported file loads as numbers, not 'a–b' text.
            gl = {**GLOSS, "Start (1-based)": "First base of this product in the scanned assembly, 1-based inclusive (isPcr's own).",
                  "End (1-based)": "Last base of this product, 1-based inclusive: Len = End − Start + 1."}
            t = DataTable(["Call", "Assembly seq", "Start (1-based)", "End (1-based)", "Len", "Strand"], gl)
            scan_table = t
            rows = [[_call(a), _amplicon_public_source(a), a["start"], a["end"], a["length"], a.get("strand", "?")]
                    for a in samps]
            if shown < len(rows):                       # a high-copy pair can hit the genome hundreds of thousands of times
                t.set_complete_rows(rows, shown)        # draw the first rows; an export still writes every product
            else:
                t.set_rows(rows)
            t.set_row_menu(lambda r: [("Copy locus (1-based)", lambda a=samps[r]: self._copy(
                f"{_amplicon_public_source(a)}:{a['start']}-{a['end']}({a.get('strand', '?')})"))])
            t.setMaximumHeight(round(260 * theme_mod.UI_SCALE))
            self.genomeBody.addWidget(t)
            pid = _pair.get("id")
            base = "_".join(["TEagle_offtarget_scan"] + [self._pair_token(x) for x in (pid, asm.get("accession")) if x])
            xrow = QHBoxLayout(); xrow.addStretch(1); xrow.addWidget(_export_table_btn(t, base, self))
            self.genomeBody.addLayout(xrow)
            if shown < len(samps):
                self.genomeBody.addWidget(_note(
                    f"Showing the first {shown:,} of {len(samps):,} products (on-target first, then by sequence and "
                    f"position). Sort the table or export it to reach the rest: an export lists all {len(samps):,}.", "info"))
        else:
            self.genomeBody.addWidget(_note("No genome-wide product for this pair under the criteria.", "warn"))
        n_on, n_off, n_single = summary.get("n_on", 0), summary.get("n_off", 0), summary.get("n_single", 0)
        head = f"{n_on} on-target + {n_off} off-target site(s)" if has_locus else f"{summary.get('n_pair', 0)} genomic priming site(s)"
        sp = f" · {n_single} single-primer artefact(s)" if n_single else ""
        comig = (" An on-target that shares a band size with an off-target is drawn in the off-target colour — a "
                 "co-migration cannot be resolved on a gel — with the full split kept in the table above." if has_locus else "")
        listing = "listed above (on-target first)" if has_locus else "listed above"
        export_listing = "listed on-target first" if has_locus else "listed"
        note = QLabel(f"{head}{sp}, {listing}.{comig} Every product is a candidate under isPcr's default ≥15 bp "
                      "minimum perfect 3′-end match requirement — a specificity screen, not wet-lab-validated bands. "
                      f"Only products up to {cap:,} bp are considered, so a longer product would not appear. The verdict is a "
                      "heuristic read of the count/spread; the numbers carry the claim.")
        note.setObjectName("orient"); note.setWordWrap(True); self.genomeBody.addWidget(note)
        if scan_table is not None:
            notes = [f"{head}{sp}, {export_listing}.{comig} Every product is a candidate under isPcr's "
                     "default >=15 bp minimum perfect 3-prime-end match requirement: this is a computational specificity screen, not wet-lab-validated "
                     "bands. The verdict is a heuristic interpretation of the displayed counts and spread."]
            if _pair.get("fwd"):
                where = (f"{_org}" if _org else "the selected genome") + (
                    f", assembly {asm['accession']}" + (f" ({asm['name']})" if asm.get("name") else "") if asm.get("accession") else "")
                notes.append(f"Primer pair {_pair.get('id') or 'pair'}: forward 5'-{_pair['fwd']}-3', reverse "
                             f"5'-{_pair.get('rev', '')}-3'. Scanned against {where}.")
            bits = [f"UCSC isPcr{' ' + str(crit['isPcr']) if crit.get('isPcr') else ''}"]
            if crit.get("min_perfect") is not None:
                bits.append(f"minimum perfect 3-prime match {crit['min_perfect']} bp")
            if crit.get("min_good") is not None:
                bits.append(f"minimum good match {crit['min_good']} bp")
            notes.append("Search: " + ", ".join(bits) + f", product sizes up to {cap:,} bp (longer products are not "
                         "reported). Coordinates are 1-based and inclusive, as isPcr reports them. The sealed run manifest "
                         "(assembly accession, file checksum, parameters) is exported from the Provenance card.")
            if shown < len(samps):
                notes.append(f"The on-screen table shows the first {shown:,} of {len(samps):,} products; this file lists all {len(samps):,}.")
            scan_table.export_notes = notes
        self._uppercase_buttons(self.genomeBody.parentWidget())   # this path has no post-render theme refresh

    # =================== WSL family annotation ===================
    def _init_wsl(self):
        self.engine.submit("wsl_status", key="wsl_status")

    @staticmethod
    def _backend_states(w):
        """Plain state of the OPTIONAL Linux backend from a wsl_status payload (None = no answer yet). 'ready' needs
        both RepeatMasker + Dfam and minimap2; the whole-genome scan is never reported as ready here, because the
        status probe does not test isPcr."""
        if w is None:
            c = "checking…"
            return {"summary": c, "rm": c, "mm": c, "wsl": c}
        no_distro = w.get("wsl2") is False and not w.get("distro")
        if w.get("error") and not no_distro:
            u = "status unknown"
            return {"summary": u, "rm": u, "mm": u, "wsl": u}
        if not w.get("wsl2"):
            n = "not installed"
            return {"summary": n, "rm": n, "mm": n, "wsl": n}
        rm = "ready" if w.get("ready") else "not set up"
        mm = "ready" if w.get("minimap2") else "not set up"
        return {"summary": rm if rm == mm else "partly installed", "rm": rm, "mm": mm, "wsl": "WSL2 found"}

    def _render_status_chip(self):
        """Compose the header chip from the two independent states, whichever arrived first."""
        core, st = self._core, self._backend_states(self._wsl)
        ok = bool(core and core.get("ok"))
        self.statusTxt.setText("connecting…" if core is None else
                               f"engine ready · primer3 {core.get('primer3')}" if ok else "engine unavailable")
        self.backendTxt.setText(f"Linux backend: {st['summary']}")
        tip = ((f"Core engine (in-process; primer3 {core.get('primer3')}): classification, ORFs, protein domains, primer "
                "design and in-silico PCR — always available." if ok else
                "Core engine: not available. Classification, ORFs, protein domains, primer design and in-silico PCR "
                "need it." if core is not None else "Core engine: starting.")
               + f"\nOptional Linux (WSL) backend: family naming (panel {PANEL['family']}) {'ready' if st['rm'] == 'ready' else 'not ready'}; "
               f"splice structure (panel {PANEL['splice']}) {'ready' if st['mm'] == 'ready' else 'not ready'}; whole-genome scan "
               f"(panel {PANEL['genome']}) needs isPcr, which is checked when a scan runs.")
        for w in (self._chip, self.statusTxt, self.backendTxt):
            w.setToolTip(tip)
        self.led.setProperty("live", ok); self._repolish(self.led)

    def _tag_backend_cards(self):
        """Panels that need the optional Linux backend say so in their header, with its live state."""
        st = self._backend_states(self._wsl)
        tag = "optional · Linux backend: "
        self.card_wsl.set_meta(f"{tag}{st['rm']} · names the TE family (Dfam)")
        self.card_splice.set_meta(f"{tag}{st['mm']} · exon–intron structure from an aligned transcript")
        self.card_genome.set_meta(f"{tag}{st['wsl']} · a designed pair vs a downloaded RefSeq genome (isPcr)")

    def _on_wsl_status(self, w):
        self._wsl = w
        self._render_status_chip()
        self._tag_backend_cards()
        self._wsl_ready = {"annotate": False, "splice": False}
        # Every probe starts from "nothing is available" and enables only what THIS status supports: a re-probe
        # after a clean reinstall (or a lost WSL) must not leave the run buttons armed or a stale library choice
        # selectable. A run already in flight keeps its own busy state.
        for kind, btn in (("annotate", self.annotateBtn), ("splice", self.spliceBtn)):
            if kind not in getattr(self, "_reqs", {}):
                btn.setEnabled(False)
        self._uncurated_available = False
        self.wslLibraryRow.setVisible(False)
        self._wsl_probed = True                           # the flags below now come from a real probe, not a default
        self._wsl_partial_uncurated = (not w.get("dfam_uncurated") and any(
            "uncurated" in str(p).lower() for p in ((w.get("dfam_library") or {}).get("partitions") or [])))
        # wsl.available() reports a MISSING distro as wsl2=False with an 'error' text; that is the install
        # guidance case, not a failed probe. A distro that exists but did not answer is the error case.
        no_distro = w.get("wsl2") is False and not w.get("distro")
        unreadable = bool(w.get("error")) and not no_distro
        stack = ("status could not be read" if unreadable else "not installed" if not w.get("wsl2") else
                 "ready" if w.get("ready") else "installed, annotation stack incomplete")
        self.backendBtn.setToolTip("Linux backend — install, repair and check the optional Linux (WSL) stack: Dfam "
                                   f"family naming, splice structure, whole-genome scan. Optional stack: {stack}.")
        if unreadable:
            self.wslStatus.setText(f"<span style='color:{theme_mod.BAD[self.theme]}'>WSL status error: {html.escape(str(w['error']))}</span>")
            self.spliceStatus.setText("WSL status could not be read, so splice structure is unavailable. Open "
                                      "<b>Backend</b> in the header to check the Linux stack.")
            return
        if not w.get("wsl2"):
            found = w.get("found") or []
            self.wslStatus.setText(("<b>WSL2 has no usable Linux distribution</b>" + f" (found: {html.escape(', '.join(found))})"
                                    if found else "<b>WSL2 not installed</b>")
                                   + " — this optional step names the Dfam family. "
                                   "The domain-based superfamily above works without it. Install it in one click with "
                                   "<b>Backend installer</b> below (it runs the elevated <code>wsl --install</code> for you; "
                                   "a Windows restart may be required), or run <code>wsl --install</code> in an Administrator PowerShell.")
            self.spliceStatus.setText("<b>WSL2 not installed</b> — splice-structure alignment needs it (optional).")
            return
        self.engine.submit("genome_list", {}, key="genome_list")   # WSL is up — populate the downloaded-genome dropdown
        # offer the library choice only once the uncurated partitions are actually on the machine —
        # a curated-only backend has nothing extra to search, so the option would promise a search it
        # cannot run. Read from the library famdb reports, not from an assumption about the install.
        # tracked explicitly, not read back off isVisible(): a widget inside a collapsed card reports
        # not-visible, and reading the choice from that would silently discard what the user selected
        # The backend sets this only when BOTH optional Dfam shards passed its file probe. Looking for an
        # "uncurated" partition name here would expose a search choice after a partial install that the
        # backend must reject.
        self._uncurated_available = bool(w.get("dfam_uncurated"))
        self.wslLibraryRow.setVisible(self._uncurated_available)
        if w.get("ready"):
            # name the partitions actually on the machine — "Dfam curated" was fixed text and read as a
            # statement about the library even where the uncurated partitions were installed alongside it
            _lib = "Dfam curated + uncurated available" if self._uncurated_available else "Dfam curated"
            self.wslStatus.setText(f"<span style='color:{theme_mod.GOOD[self.theme]}'>● ready</span> · RepeatMasker {w.get('repeatmasker')} "
                                   f"· {_lib} · distro {w.get('distro')}")
            self.annotateBtn.setEnabled(True); self._wsl_ready["annotate"] = True
        else:
            self.wslStatus.setText(f"WSL2 ok ({w.get('distro')}); annotation stack not installed "
                                   f"(RepeatMasker {w.get('repeatmasker') or 'missing'}, Dfam {'ok' if w.get('dfam') else 'missing'}).")
            self.wslInstallBtn.setVisible(True)
        if w.get("minimap2"):
            self.spliceStatus.setText(f"<span style='color:{theme_mod.GOOD[self.theme]}'>● ready</span> · minimap2 {w.get('minimap2')} "
                                      "· align a transcript to resolve exon–intron structure")
            self.spliceBtn.setEnabled(True); self._wsl_ready["splice"] = True
        elif w.get("wsl2"):
            self.spliceStatus.setText("minimap2 not installed in the WSL backend — it ships with the managed install (Backend in the header).")

    def _on_species_changed(self):
        """Show the free-text field only for 'Other…'; clear it when switching to a listed organism
        so a stale lineage can't silently drive the next annotation."""
        other = self.wslSpecies.currentData() == "__other__"
        self.wslSpeciesOther.setVisible(other)
        if not other:
            self.wslSpeciesOther.clear()
        if not getattr(self, "_species_busy", False):
            self._species_auto = False                        # the user chose it: it outlives a specimen change

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
        self._species_busy = True
        try:
            for i in range(self.wslSpecies.count()):
                d = self.wslSpecies.itemData(i)
                if isinstance(d, str) and d != "__other__" and d.lower() == name.lower():
                    self.wslSpecies.setCurrentIndex(i); break
            else:
                oi = self.wslSpecies.findData("__other__")    # not a listed organism -> Other + free text
                if oi >= 0:
                    self.wslSpecies.setCurrentIndex(oi)
                    self.wslSpeciesOther.setText(name)
        finally:
            self._species_busy = False
        self._species_auto = True                             # taken from the specimen, so the next specimen replaces it

    def _drop_specimen_species(self):
        """The organism came from the previous specimen (fetch or bundled example): clear it rather than run
        RepeatMasker with a lineage that does not describe the new sequence. A choice the user made stays."""
        if getattr(self, "_species_auto", False) and hasattr(self, "wslSpecies"):
            self._species_busy = True
            try:
                self.wslSpecies.setCurrentIndex(0)
            finally:
                self._species_busy = False
            self._species_auto = False

    def _annotate(self):
        self._clear_banner()                                  # drop any stale error before a retry (matches _design/_run_pcr)
        if self.wslSource.currentIndex() == 1:
            seq = self.wslPaste.toPlainText().strip()
            src = None
        else:
            seq = self.state.get("seq") or self.seq.toPlainText().strip()
            seq, src = self._operation_input(seq)
        if not seq:
            return self._banner("Load or paste a specimen first — there is no sequence to annotate.", "info")
        self.state["family_seq"] = self._norm_seq(seq)    # hit coords index THIS sequence (backend-normalized), not always panel-01
        self.annotateBtn.setEnabled(False); widgets.set_caps(self.annotateBtn, "◴ annotating…")
        self._set_body(self.wslBody, BusyBar("Running RepeatMasker against Dfam — this can take a minute or two…"))
        self.engine.submit("annotate", {"sequence": seq, "species": self._species(), "source": src,
                                        "include_uncurated": bool(self.wslLibrary.currentData())
                                        if getattr(self, "_uncurated_available", False) else False},
                           key=self._issue("annotate", seq, bound=self.wslSource.currentIndex() != 1))

    def _on_annotate(self, d):
        self.annotateBtn.setEnabled(True); widgets.set_caps(self.annotateBtn, "Run family annotation")
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

    # The family tables are RepeatMasker output, not the HMMER domain panel the shared glossary describes.
    _FAMILY_GLOSS = {
        "Score": "RepeatMasker's Smith-Waterman alignment score (RMBlast) for this match. It is not a HMMER bit "
                 "score and cannot be compared with the domain table's Score.",
        "Class/family": "RepeatMasker class/family label of the matched Dfam entry, e.g. LTR/Copia. It also covers "
                        "non-TE classes (satellites, simple repeats, RNA genes), which are listed apart.",
        "Div (raw)": "Raw % substitutions between the query and the Dfam consensus over the aligned span; not "
                     "Kimura-corrected. Lower is closer.",
        "Del": "% of consensus bases missing from the query inside the alignment (RepeatMasker %del).",
        "Ins": "% of query bases with no counterpart in the consensus inside the alignment (RepeatMasker %ins).",
        "Consensus": "Aligned span on the Dfam consensus (begin–end) over the consensus length, in bp.",
        "Cov": "Percent of the Dfam consensus covered by the alignment (summed over merged blocks).",
        "Blocks": "Number of RepeatMasker alignment rows merged into this row because they share a fragment ID "
                  "(one interrupted element).",
    }

    def _repeat_table(self, hits):
        row = lambda i, h, b: [i + 1, h["class_family"], h["family"], self._span(h["q_start"], h["q_end"], b),
                               h["strand"], f"{h['divergence']}%", h["score"]]
        t = self._coord_table(["#", "Class/family", "Dfam family", self._ch("Coords"), "Str", "Div", "Score"],
                              self._gloss(self._FAMILY_GLOSS), [row(i, h, 0) for i, h in enumerate(hits)])
        t.export_header_map = {"Coords (1-based)": "Coords (0-based)"}
        t.set_rows([row(i, h, self.coord_base) for i, h in enumerate(hits)])
        return t

    @classmethod
    def _partition_hits(cls, hits):
        """(te, tandem, unclassified, other) by the backend's own classifier (wsl.repeat_kind), so this panel and
        the whole-genome landscape agree. 'other' holds RNA genes, artefacts and other non-TE annotations."""
        from teagle_core.wsl import repeat_kind
        te, tandem, unclassed, other = [], [], [], []
        for h in hits:
            kind = repeat_kind(h["class_family"])
            if kind == "TE":
                te.append(h)
            elif kind == "tandem":
                tandem.append(h)
            elif h["class_family"].split("/")[0] in cls._UNCLASSED:
                unclassed.append(h)
            else:
                other.append(h)
        return te, tandem, unclassed, other

    def _render_family(self, d):
        self.state["familyResult"] = d                        # kept so the coordinate display can redraw it
        hits = d.get("hits", [])
        te, nonte, unclassed, other = self._partition_hits(hits)
        if not te:
            cont = QWidget(); cl = QVBoxLayout(cont); cl.setContentsMargins(0, 0, 0, 0); cl.setSpacing(6)
            if nonte:
                kinds = sorted({h["class_family"].replace("_", " ").lower() for h in nonte})
                nonte_note = (f"RepeatMasker annotated {len(nonte)} region(s) of this locus as "
                              f"{' and '.join(kinds)}. That is a result, not a failure: those regions match a "
                              f"tandem/low-complexity repeat, and no transposable-element family was found in "
                              f"them.")
                cl.addWidget(_note(nonte_note, "info"))
                nonte_table = self._repeat_table(nonte)
                nonte_table.export_notes = [nonte_note]
                cl.addWidget(nonte_table)
                nonte_row = QHBoxLayout(); nonte_row.addStretch(1)
                nonte_row.addWidget(_export_table_btn(nonte_table, "TEagle_non_TE_repeats", self))
                cl.addLayout(nonte_row)
            if unclassed:
                unclassed_note = (f"{len(unclassed)} repeat region(s) matched a Dfam entry whose class is "
                                  f"unassigned ('Unknown'). A repeat is present; its type is not established.")
                cl.addWidget(_note(unclassed_note, "info"))
                unclassed_table = self._repeat_table(unclassed)
                unclassed_table.export_notes = [unclassed_note]
                cl.addWidget(unclassed_table)
                unclassed_row = QHBoxLayout(); unclassed_row.addStretch(1)
                unclassed_row.addWidget(_export_table_btn(unclassed_table, "TEagle_unclassified_repeats", self))
                cl.addLayout(unclassed_row)
            if other:
                other_note = (f"{len(other)} region(s) of this locus were annotated as RNA genes or other non-TE "
                              "features (for example snRNA, tRNA or an artefact). They are listed apart and are "
                              "not transposable-element family calls.")
                cl.addWidget(_note(other_note, "info"))
                other_table = self._repeat_table(other)
                other_table.export_notes = [other_note]
                cl.addWidget(other_table)
                other_row = QHBoxLayout(); other_row.addStretch(1)
                other_row.addWidget(_export_table_btn(other_table, "TEagle_other_annotations", self))
                cl.addLayout(other_row)
            if True:                                      # the blank-result caveats apply whenever no TE family matched
                # Two independent conditions decide whether a family CAN be named, and the message names
                # whichever is actually missing rather than guessing. Measured on Drosophila copia: with
                # the uncurated partitions installed AND a species it resolves to Copia_LTR and Copia_I at
                # 100% consensus coverage; with the partitions but NO species it returns only
                # low-complexity, because RepeatMasker searches a limited default set without a lineage.
                # How much the curated-only case costs depends entirely on the lineage, so the message
                # below no longer names a figure for one organism — see _CURATED_COVERAGE.
                parts = ((d.get("dfam_library") or {}).get("partitions")) or []
                searched_unc = bool(d.get("include_uncurated"))
                sp = d.get("species") or ""
                defaulted = d.get("species_defaulted")
                has_sp = (not defaulted) if defaulted is not None else (bool(sp) and not str(sp).startswith("(all"))
                if getattr(self, "_wsl_probed", False):   # what the backend probe established, not partition names
                    has_unc, partial_unc = bool(self._uncurated_available), bool(self._wsl_partial_uncurated)
                else:
                    has_unc, partial_unc = any("uncurated" in str(p).lower() for p in parts), False
                blank = "No Dfam family matched" if not hits else "No transposable-element family matched"
                if not has_sp:
                    lineage = d.get("lineage_searched") or "Homo sapiens"
                    cl.addWidget(_note(f"{blank} — and no organism was set, so RepeatMasker searched its default "
                                       f"lineage ({lineage}). It searches only a limited set without the right "
                                       "lineage, so an organism is usually what turns a blank result into a named "
                                       "family. Set the organism above and run again.", "warn"))
                elif not searched_unc and partial_unc:
                    cl.addWidget(_note(f"{blank} for “{sp}” — this run searched the CURATED families only. Only one of "
                                       "the two optional uncurated Dfam partitions is present, so the Library choice "
                                       "is not offered: one partition is missing. Repair it in BACKEND in the "
                                       "header, then run again with <b>Library</b> set to include uncurated "
                                       "families.", "warn"))
                elif not searched_unc and not has_unc:
                    cl.addWidget(_note(f"{blank} for “{sp}” — and only the CURATED Dfam "
                                       "partitions are installed. How much that costs depends entirely on "
                                       "the lineage: " + _curated_coverage() + " A blank result here may "
                                       "therefore be a limit of the installed library rather than a "
                                       "property of the sequence. Install the optional uncurated "
                                       "partitions from BACKEND in the header, then run again with "
                                       "<b>Library</b> set to include them.", "warn"))
                elif not searched_unc:
                    cl.addWidget(_note(f"{blank} for “{sp}” — but this run searched the "
                                       "CURATED families only. The uncurated partitions are installed on "
                                       "this machine and were not read: RepeatMasker searches curated "
                                       "families unless it is asked for both. Set <b>Library</b> to "
                                       "“Include uncurated families” and run again before concluding "
                                       "anything from this blank result.", "warn"))
                else:
                    cl.addWidget(_note(f"{blank} in this locus for “{sp}”. The curated and "
                                       "uncurated partitions were both searched with a lineage set, so "
                                       "this is a genuine no-match rather than a gap in the installed "
                                       "library. The structural and protein-domain evidence above does "
                                       "not depend on this search.", "warn"))
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
        # the header names the class of the strongest TE match, not of the leftmost one
        top = max(te, key=lambda h: h.get("score") or 0)
        extra = len({h["class_family"] for h in te} - {top["class_family"]})
        species = (f"not set — RepeatMasker searched its default lineage ({d.get('lineage_searched') or 'Homo sapiens'})"
                   if d.get("species_defaulted") else html.escape(str(d.get("species", ""))))
        head = QLabel(f"<b>Dfam · {top['class_family']}</b>"
                      + (f" (strongest match; {extra} other TE class{'es' if extra > 1 else ''} below)" if extra else "")
                      + f" — {' · '.join(sorted({h['family'] for h in te}))} "
                      f"· {_dfam_lbl}{self._src_html('Dfam')} · RepeatMasker {d.get('repeatmasker_version','')}"
                      f"{self._src_html('RepeatMasker')} · species: {species}")
        head.setTextFormat(Qt.RichText); head.setWordWrap(True); _external(head); _kb_links(head)
        # %div/%del/%ins and the consensus-side coverage are RepeatMasker's own alignment evidence; they
        # were parsed and discarded before. Divergence is RAW (not Kimura-corrected) and the header says so.
        headers = ["#", "Class/family", "Dfam family", self._ch("Coords"), "Str",
                   "Div (raw)", "Del", "Ins", "Consensus", "Cov", "Blocks", "Score"]

        def _cons(h):
            if h.get("cons_start") is None or h.get("cons_length") is None:
                return "—"
            return f"{h['cons_start']}–{h['cons_end']} / {h['cons_length']}"

        fam_row = lambda i, h, b: [i + 1, h["class_family"], h["family"], self._span(h["q_start"], h["q_end"], b), h["strand"],
                                   f"{h['divergence']}%",
                                   ("—" if h.get("pct_del") is None else f"{h['pct_del']}%"),
                                   ("—" if h.get("pct_ins") is None else f"{h['pct_ins']}%"),
                                   _cons(h),
                                   ("—" if h.get("cons_coverage_pct") is None else f"{h['cons_coverage_pct']}%"),
                                   h.get("n_fragments", 1), h["score"]]
        t = self._coord_table(headers, self._gloss(self._FAMILY_GLOSS), [fam_row(i, h, 0) for i, h in enumerate(te)])
        t.export_header_map = {"Coords (1-based)": "Coords (0-based)", "Consensus": "Consensus (1-based, inclusive)"}
        t.set_rows([fam_row(i, h, self.coord_base) for i, h in enumerate(te)])
        family_notes = [
            f"Searched {_dfam_lbl} with RepeatMasker {d.get('repeatmasker_version','')} for species "
            f"{d.get('species','') or 'not set'}.",
            "Search ran in RepeatMasker's rush mode (-qq): about 10% less sensitive than a default search and 4 to 10 "
            "times faster (RepeatMasker's own description). Divergence is the raw substitution percentage, not Kimura-corrected; consensus coverage and "
            "divergence are alignment evidence, not a family-membership test.",
        ]
        t.export_notes = family_notes
        t.set_row_menu(lambda r: self._feat_menu(te[r]["q_start"], te[r]["q_end"], te[r]["strand"],
                                                 te[r]["family"], src_seq=self.state.get("family_seq")))
        cont = QWidget(); cl = QVBoxLayout(cont); cl.setContentsMargins(0, 0, 0, 0); cl.setSpacing(6)
        cl.addWidget(head); cl.addWidget(t)
        merged = [h for h in te if h.get("n_fragments", 1) > 1]
        if merged:
            merged_note = ("; ".join(h["fragment_note"] for h in merged) +
                           ". Merged so an interrupted element is one row, with the outermost query span "
                           "and the summed consensus coverage.")
            cl.addWidget(_note(merged_note, "info"))
            t.export_notes.append(merged_note)
        cl.addWidget(_qq_note())
        if nonte or unclassed or other:            # reported, never silently dropped
            also = [f"{len(nonte)} tandem / low-complexity region(s)"] if nonte else []
            also += [f"{len(unclassed)} repeat(s) of unassigned class"] if unclassed else []
            also += [f"{len(other)} RNA-gene or other non-TE region(s)"] if other else []
            cl.addWidget(_note("Also annotated in this locus: " + " and ".join(also) +
                               ", listed apart from the TE families above.", "info"))
        def add_partition(partition, title, base, note_text):
            if not partition:
                return
            cl.addWidget(_sl(title))
            pt = self._repeat_table(partition)
            pt.export_notes = [note_text] + family_notes
            cl.addWidget(pt)
            prow = QHBoxLayout(); prow.addStretch(1)
            prow.addWidget(_export_table_btn(pt, base, self)); cl.addLayout(prow)

        add_partition(nonte, "Tandem / low-complexity repeats", "TEagle_non_TE_repeats",
                      "These rows are positive tandem/low-complexity repeat annotations reported separately; "
                      "they are not transposable-element family calls.")
        add_partition(unclassed, "Unclassified repeats", "TEagle_unclassified_repeats",
                      "These rows matched repeat entries whose class is unassigned. A repeat is present, but its "
                      "type is not established; exclusion from TE coverage is not evidence that it is non-TE.")
        add_partition(other, "RNA genes and other non-TE annotations", "TEagle_other_annotations",
                      "These rows are annotations of RNA genes, artefacts or other non-TE features, reported "
                      "separately; they are not transposable-element family calls.")
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

    # =================== provenance ===================
    def _prov_list(self):
        """Runs since the current analysis was accepted (or, with no analysis, since the specimen last changed), oldest
        first; each entry holds the manifests of one run (a PCR batch holds one per pair). The list restarts whenever
        the analysis it belongs to is replaced or invalidated, so it never mixes runs that no longer describe the screen."""
        analysis = self.state.get("analysis") if self.state.get("analysis_state") == "current" else None
        anchor = getattr(self, "_prov_analysis", None)
        if not hasattr(self, "_prov_runs") or (anchor is not None and anchor is not analysis):
            self._prov_runs = []
            self._prov_sel = None
        self._prov_analysis = analysis
        return self._prov_runs

    def _pcr_group(self, m):
        """The manifests of the whole PCR batch that `m` heads: the card is drawn for the first pair as soon as the
        batch is in, and the batch has one manifest per pair."""
        run = getattr(self, "_pcr_run", None)
        if m.get("runType") != "in-silico-pcr" or not run:
            return None
        provs = [r["provenance"] for r in run.get("results") or [] if isinstance(r, dict) and r.get("provenance")]
        return provs if provs and provs[0] is m else None

    def _render_provenance(self, m, group=None):
        self.card_prov.expand()
        runs = self._prov_list()
        members = list(group) if group else (self._pcr_group(m) or [m])
        runs.append({"manifests": members})
        # the latest run is what 'Export this run' and the other panels' callers see; choosing another run in the list
        # only changes what is shown
        self.state["prov_manifest"] = m
        self.state["prov_manifest_all"] = list(members)
        self._prov_sel = (len(runs) - 1, 0)
        self._build_prov_card()
        self._uppercase_buttons()

    def _build_prov_card(self):
        import report
        card = self.card_prov
        card.clear_body()
        body = card.bodylay
        self._prov_report_block(body)
        runs = self._prov_list()
        op, member = self._prov_sel if self._prov_sel else (len(runs) - 1, 0)
        body.addWidget(_sl("Runs"))
        combo = _Combo(); combo.setAccessibleName("Run to show")
        for i, run in enumerate(runs):
            k = len(run["manifests"])
            combo.addItem(report.run_title(i + 1, run["manifests"][0], f" · {k} pairs" if k > 1 else ""), i)
        combo.setCurrentIndex(op)
        combo.setToolTip("Every run since the specimen or its analysis last changed. Earlier runs are cleared because "
                         "they no longer describe the current result.")
        pair = _Combo(); pair.setAccessibleName("Pair to show")
        pair.setToolTip("A batch of in-silico PCR runs has one sealed manifest per primer pair.")
        holder = QWidget(); QVBoxLayout(holder).setContentsMargins(0, 0, 0, 0); holder.layout().setSpacing(theme_mod.sp(9))
        body.addWidget(combo); body.addWidget(pair); body.addWidget(holder)
        self._prov_view = holder

        def fill_pairs(op_index, keep=0):
            pair.blockSignals(True); pair.clear()
            mans = runs[op_index]["manifests"]
            for j, name in enumerate(self._pair_names(mans)):
                pair.addItem(name, j)
            pair.setCurrentIndex(min(keep, pair.count() - 1))
            pair.setVisible(len(mans) > 1)
            pair.blockSignals(False)

        def show(op_index):
            fill_pairs(op_index)
            self._render_run_view(op_index, 0)
            self._uppercase_buttons(holder)

        combo.currentIndexChanged.connect(show)
        pair.currentIndexChanged.connect(lambda j: (self._render_run_view(combo.currentIndex(), j),
                                                    self._uppercase_buttons(holder)))
        fill_pairs(op, member)
        self._render_run_view(op, member)

    def _pair_names(self, manifests):
        """Name of each manifest's primer pair as the gel and the tables call it (P1, P2 … when it is not known)."""
        known = {(p.get("fwd"), p.get("rev")): p.get("label") for p in (self.state.get("lastPcr") or {}).get("pairs") or []}
        out = []
        for k, m in enumerate(manifests, 1):
            prm = m.get("parameters") or {}
            out.append(known.get((prm.get("fwd"), prm.get("rev"))) or f"P{k}")
        return out

    def _render_run_view(self, op_index, member_index=0):
        import report
        runs = self._prov_list()
        if not (0 <= op_index < len(runs)):
            return
        mans = runs[op_index]["manifests"]
        m = mans[max(0, min(member_index, len(mans) - 1))]
        self._prov_sel = (op_index, member_index)
        lay = self._prov_view.layout()
        _clear_layout(lay)
        lay.addWidget(self._kv_section("Input", report.input_rows(m)))
        lay.addWidget(self._kv_section("Run and seal", report.seal_rows(m)))
        sw = m.get("software") or []
        if sw:
            lay.addWidget(self._kv_section("Software", [(s.get("name", "—"), s.get("version") or "—", None) for s in sw]))
        dbs = m.get("databases") or []
        if dbs:
            self._db_sections(lay, m, dbs)
        lay.addWidget(_sl("Parameters"))
        sections = report.param_sections(m)
        if not sections:
            lay.addWidget(_empty("This run sealed no parameters (defaults only)."))
        raw = {}
        for group, rows in sections:
            for r in rows:
                if r["raw"]:
                    raw[r["key"]] = json.loads(r["raw"])
            table = self._param_grid(rows)
            if group == report.INTERNALS:
                self._disclosure(lay, f"{group} · {len(rows)} values", table, expanded=False)
            else:
                box = QWidget(); bl = QVBoxLayout(box); bl.setContentsMargins(0, 0, 0, 0); bl.setSpacing(theme_mod.sp(3))
                head = QLabel(f"<b>{html.escape(group)}</b>"); head.setObjectName("cardnote"); head.setTextFormat(Qt.RichText)
                bl.addWidget(head); bl.addWidget(table)
                lay.addWidget(box)
        if raw:
            box = QPlainTextEdit(json.dumps(raw, indent=2, ensure_ascii=False, sort_keys=True)); box.setReadOnly(True)
            box.setMaximumHeight(theme_mod.sp(220)); box.setAccessibleName("Raw sealed parameters (JSON)")
            self._disclosure(lay, "Raw sealed parameters (JSON)", box, expanded=False)
        env = m.get("environment") or {}
        if env:
            lay.addWidget(self._kv_section(f"Environment ({report.RESULT_NOTE})",
                                           [(k, str(env[k]), None) for k in ("os", "machine", "python") if k in env]))
        notrun = m.get("notRun") or []
        if notrun:
            n = QLabel("<br>".join("• " + html.escape(str(x)) for x in notrun)); n.setObjectName("cardnote")
            n.setTextFormat(Qt.RichText); n.setWordWrap(True); n.setTextInteractionFlags(Qt.TextSelectableByMouse)
            lay.addWidget(_sl("Not run in this step")); lay.addWidget(n)
        refs = m.get("references") or []
        if refs:
            r = QLabel("<br>".join(f"<b>{html.escape(str(x.get('name', '')))}</b> — {html.escape(str(x.get('citation', '')))}"
                                   + (f" doi:{html.escape(str(x['doi']))}" if x.get("doi") else "") for x in refs))
            r.setObjectName("cardnote"); r.setTextFormat(Qt.RichText); r.setWordWrap(True)
            r.setTextInteractionFlags(Qt.TextSelectableByMouse)
            lay.addWidget(_sl("References (source-verified)")); lay.addWidget(r)
        self._prov_export_row(lay, op_index, len(mans))

    def _db_sections(self, lay, m, dbs):
        import report
        rows = []
        for d in dbs:
            name = str(d.get("name") or "—")
            if d.get("version"):
                rows.append((name, str(d["version"]), None))
            elif d.get("file") and not d.get("sha256"):
                rows.append((name, str(d["file"]), None))
            if d.get("sha256"):
                rows.append((f"{name} SHA-256", d["sha256"], "hash"))
        profiles, note = report.panel_profiles(m)
        if profiles:
            rows.append(("Pfam models", f"{len(profiles)} (accessions below)", None))
        lay.addWidget(self._kv_section("Databases", rows))
        if note:
            lay.addWidget(_note(note, "info"))
        if profiles:
            grid = QGridLayout(); grid.setHorizontalSpacing(theme_mod.sp(14)); grid.setVerticalSpacing(theme_mod.sp(2))
            for c, head in enumerate(("Code", "Model", "Accession")):
                h = QLabel(head); h.setObjectName("kdim"); grid.addWidget(h, 0, c)
            for r, p in enumerate(profiles, 1):
                for c, text in enumerate((p["code"], p["name"], p["accession"])):
                    v = QLabel(text); v.setObjectName("cardmeta"); v.setTextInteractionFlags(Qt.TextSelectableByMouse)
                    v.setToolTip(p["label"]); grid.addWidget(v, r, c)
            grid.setColumnStretch(3, 1)
            box = QWidget(); box.setLayout(grid)
            self._disclosure(lay, f"Profile accessions ({len(profiles)})", box, expanded=False)

    def _param_grid(self, rows):
        box = QWidget(); grid = QGridLayout(box); grid.setContentsMargins(0, 0, 0, 0)
        grid.setHorizontalSpacing(theme_mod.sp(10)); grid.setVerticalSpacing(theme_mod.sp(3)); grid.setColumnStretch(1, 1)
        for r, row in enumerate(rows):
            k = QLabel(row["label"]); k.setObjectName("kdim"); k.setWordWrap(True); k.setToolTip(row["tip"])
            v = QLabel(row["text"]); v.setObjectName("cardmeta"); v.setWordWrap(True); v.setToolTip(row["tip"])
            v.setTextInteractionFlags(Qt.TextSelectableByMouse)
            grid.addWidget(k, r, 0, Qt.AlignTop); grid.addWidget(v, r, 1, Qt.AlignTop)
        grid.setColumnMinimumWidth(0, theme_mod.sp(230))      # one label column width for every group
        return box

    def _kv_section(self, title, rows):
        """A labelled block of key/value rows. A row is (label, value, kind); a 'hash' value is shown in full, in blocks
        of eight characters so it wraps, beside a Copy button that puts the exact 64 characters on the clipboard."""
        import report
        box = QWidget(); lay = QVBoxLayout(box); lay.setContentsMargins(0, 0, 0, 0); lay.setSpacing(theme_mod.sp(4))
        lay.addWidget(_sl(title))
        grid = QGridLayout(); grid.setHorizontalSpacing(theme_mod.sp(10)); grid.setVerticalSpacing(theme_mod.sp(3))
        grid.setColumnStretch(1, 1)
        for r, (label, value, kind) in enumerate(rows):
            k = QLabel(label); k.setObjectName("kdim"); k.setWordWrap(True)
            v = QLabel(report.group_hash(value) if kind == "hash" else str(value)); v.setObjectName("cardmeta")
            v.setWordWrap(True); v.setTextInteractionFlags(Qt.TextSelectableByMouse)
            grid.addWidget(k, r, 0, Qt.AlignTop); grid.addWidget(v, r, 1, Qt.AlignTop)
            if kind == "hash":
                b = QPushButton("Copy"); b.setProperty("sm", True); b.setAccessibleName(f"Copy {label}")
                b.setToolTip(f"Copy the full 64-character value of {label}")
                b.clicked.connect(lambda _=False, text=value, btn=b: self._copy_value(text, btn))
                grid.addWidget(b, r, 2, Qt.AlignTop)
        grid.setColumnMinimumWidth(0, theme_mod.sp(150))
        lay.addLayout(grid)
        return box

    def _copy_value(self, text, button=None):
        QApplication.clipboard().setText(str(text))
        if button is not None:
            old = button.text()
            button.setText(widgets.caps("Copied"))

            def restore(b=button, t=old):
                try:
                    b.setText(t)
                except RuntimeError:
                    pass
            QTimer.singleShot(1500, restore)

    def _prov_export_row(self, lay, op_index, n_members):
        """Exports for the run on show. The analysis' own manifest and complete result stay reachable whichever run is
        selected, while the analysis is current."""
        runs = self._prov_list()
        shown = runs[op_index]["manifests"][0]
        analysis = self.state.get("analysis") if self.state.get("analysis_state") == "current" else None
        row = QHBoxLayout(); row.addStretch(1)
        b = QPushButton("Export this run (.json)"); b.setProperty("sm", True); b.setProperty("resultExport", True)
        b.setToolTip(f"Write all {n_members} pair manifests of this batch as one JSON file." if n_members > 1 else
                     "Write the full run-provenance manifest (input hash, tool and database versions, parameters, "
                     "checksums) as a JSON file so the result stays reproducible outside TEagle.")
        b.clicked.connect(self._export_selected_run); row.addWidget(b)
        if len(runs) > 1:
            a = QPushButton("Export all runs (.json)"); a.setProperty("sm", True); a.setProperty("resultExport", True)
            a.setToolTip("Write every sealed manifest listed above, in order, as one JSON file.")
            a.clicked.connect(self._export_all_runs); row.addWidget(a)
        lay.addLayout(row)
        if isinstance(analysis, dict):
            row2 = QHBoxLayout(); row2.addStretch(1)
            if shown.get("runType") != "analysis":
                ab = QPushButton("Export analysis manifest (.json)"); ab.setProperty("sm", True)
                ab.setProperty("resultExport", True)
                ab.setToolTip("Write the sealed provenance manifest of the analysis itself.")
                ab.clicked.connect(self._export_analysis_manifest); row2.addWidget(ab)
            rb = QPushButton("Export analysis result (.json)"); rb.setProperty("sm", True)
            rb.setProperty("resultExport", True)
            rb.setToolTip("Write the complete versioned analysis result: every record, stage status, "
                          "scientific output, reference and verified provenance binding.")
            rb.clicked.connect(self._export_analysis_result); row2.addWidget(rb)
            lay.addLayout(row2)

    def _export_manifest_payload(self, manifests, base):
        """Strictly verify, serialize, and atomically write one or more sealed manifests."""
        import json
        from teagle_core import provenance
        allm = [m for m in manifests if isinstance(m, dict)]
        if not allm:
            return self._banner("No run provenance to export — run the analysis first.", "info")
        path, _ = QFileDialog.getSaveFileName(self, "Export provenance manifest",
                                              f"{base}.json", "JSON (*.json)")
        if not path:
            return
        if not path.lower().endswith(".json"):
            path += ".json"
        try:
            encoded_runs = [provenance.manifest_json_bytes(manifest) for manifest in allm]
            payload = (json.dumps({"runs": [json.loads(encoded) for encoded in encoded_runs]},
                                  indent=2, ensure_ascii=False, allow_nan=False, sort_keys=True).encode("utf-8")
                       if len(encoded_runs) > 1 else encoded_runs[0])
        except (TypeError, ValueError):
            return self._banner("The provenance manifest failed its seal or strict JSON check and was not exported.")
        if not self._write_bytes(path, payload):
            return
        self._banner(f"Provenance manifest written to {os.path.basename(path)}.", level="success")

    def _write_bytes(self, path, payload):
        """Write an export atomically; a disk problem is reported as one, never as a seal failure."""
        try:
            with widgets.atomic_write(path, "wb") as fh:
                fh.write(payload)
        except OSError:
            self._banner(f"“{os.path.basename(path)}” could not be written. The file may be open in another "
                         "program, or the folder may be read-only or no longer available. The provenance itself "
                         "is intact; choose another location and export again.")
            return False
        return True

    def _export_analysis_manifest(self):
        res = self.state.get("analysis")
        manifest = (res or {}).get("provenance") if isinstance(res, dict) else None
        if self.state.get("analysis_state") != "current" or not isinstance(manifest, dict):
            return self._banner("No current analysis to export — run analysis again first.", "info")
        stem = ("multiple_records" if len(res.get("records") or []) > 1
                else _manifest_export_stem(manifest.get("input", {}), "analysis"))
        return self._export_manifest_payload([manifest], f"TEagle_{stem}_manifest")

    def _export_manifest(self):
        allm = self.state.get("prov_manifest_all") or ([self.state["prov_manifest"]]
                                                        if self.state.get("prov_manifest") else [])
        return self._export_run_set(allm)

    def _export_run_set(self, allm):
        """Write the manifests of one run (or batch); each is re-verified against its seal first."""
        if not allm:
            return self._banner("No run provenance to export — run an analysis first.", "info")
        if allm[0].get("runType") == "analysis" and self.state.get("analysis_state") != "current":
            return self._banner("That analysis provenance is stale — run analysis again before exporting.", "info")
        if len(allm) > 1:
            stem = "batch"
        elif allm[0].get("runType") == "analysis" and len(self.state.get("records") or []) > 1:
            stem = "multiple_records"
        else:
            stem = _manifest_export_stem(allm[0].get("input", {}), "run")
        return self._export_manifest_payload(allm, f"TEagle_{stem}_manifest")

    def _export_selected_run(self):
        """'Export this run': the run chosen in the provenance list, which is the latest one until another is picked."""
        runs = self._prov_list()
        sel = self._prov_sel[0] if self._prov_sel else len(runs) - 1
        return self._export_run_set(runs[sel]["manifests"] if 0 <= sel < len(runs) else [])

    def _export_all_runs(self):
        """The session's provenance bundle: every sealed manifest since the analysis, in the order the runs happened."""
        allm = [m for run in self._prov_list() for m in run["manifests"]]
        if not allm:
            return self._banner("No run provenance to export — run an analysis first.", "info")
        recs = self.state.get("records") or []
        stem = "multiple_records" if len(recs) > 1 else _manifest_export_stem(allm[0].get("input", {}), "run")
        return self._export_manifest_payload(allm, f"TEagle_{stem}_all_runs_manifest")

    def _export_analysis_result(self):
        import engine
        result = self.state.get("analysis")
        manifest = result.get("provenance") if isinstance(result, dict) else None      # the analysis' own seal, whatever step ran since
        if (self.state.get("analysis_state") != "current" or not isinstance(manifest, dict) or
                manifest.get("runType") != "analysis"):
            return self._banner("No current analysis result to export — run analysis again first.", "info")
        rid = ("multiple_records" if len(result.get("records") or []) > 1 else
               _manifest_export_stem(manifest.get("input", {}), "analysis"))
        path, _ = QFileDialog.getSaveFileName(self, "Export complete analysis result",
                                              f"TEagle_{rid}_analysis.json", "JSON (*.json)")
        if not path:
            return
        if not path.lower().endswith(".json"):
            path += ".json"
        try:
            payload = engine.analysis_result_json_bytes(result)
        except ValueError:
            return self._banner("The analysis result failed its provenance binding or strict JSON check and was not exported.")
        if not self._write_bytes(path, payload):
            return
        self._banner(f"Complete analysis result written to {os.path.basename(path)}.", level="success")

    # =================== assay record ===================
    def _prov_report_block(self, layout):
        """Record-keeping actions for the record on show: a notebook entry, a one-page assay report and the
        classification table. Present only while the analysis is current."""
        if self.state.get("analysis_state") != "current" or not isinstance(self.state.get("analysis"), dict):
            return
        layout.addWidget(_sl("Record this assay"))
        row = QHBoxLayout()
        copy = QPushButton("Copy notebook entry"); copy.setProperty("sm", True); copy.setProperty("resultExport", True)
        copy.setToolTip("Copies a plain-text and formatted summary (call and cautions, primer pairs with locations and "
                        "QC, in-silico PCR result, design and QC conditions, run seals) for a lab notebook. Figures are "
                        "in the assay report.")
        copy.clicked.connect(self._copy_notebook_entry)
        rep = QPushButton("Assay report ▾"); rep.setProperty("sm", True); rep.setProperty("resultExport", True)
        rep.setToolTip("One-page record: call and cautions, element map, primer pairs, predicted gel, conditions and "
                       "provenance.")
        rep.clicked.connect(lambda _=False, b=rep: self._assay_report_menu(b))
        tsv = QPushButton("Export classification (TSV)"); tsv.setProperty("sm", True); tsv.setProperty("resultExport", True)
        tsv.setToolTip("One row per record in TEsorter's .cls.tsv column order (#TE, Order, Superfamily, Clade, Complete, "
                       "Strand, Domains), followed by the RepeatMasker-style class, the Wicker code and TEagle's own call, "
                       "outcome, confidence and completeness tier, for joining against TEsorter, EDTA or RepeatMasker "
                       "results. Clade is NA because TEagle does not resolve clades. When record IDs differ between tools, "
                       "join on Sequence_SHA256.")
        tsv.clicked.connect(self._export_classification)
        for b in (copy, rep, tsv):
            row.addWidget(b)
        row.addStretch(1)
        layout.addLayout(row)
        cap = QLabel(f"Covers the record selected in panel {PANEL['struct']} and the primer pairs designed on it. "
                     "Coordinates are 0-based, half-open [start, end).")
        cap.setObjectName("orient"); cap.setWordWrap(True)
        layout.addWidget(cap)

    def _export_classification(self):
        import report
        recs = self.state.get("records") or []
        if self.state.get("analysis_state") != "current" or not recs:
            return self._banner("No current analysis — run analysis again before exporting the classification table.", "info")
        stem = "multiple_records" if len(recs) > 1 else _record_export_stem(recs[0], "record")
        self._save_text(report.classification_tsv(recs, _rank_view, _call_view), f"TEagle_{stem}_classification.cls", "tsv")

    def _call_display(self, rec):
        """The call as the results card words it, from the same module-level views, for the notebook and the report."""
        cl = rec.get("classification") or {}
        te = cl.get("te_class", "—")
        held = cl.get("outcome") == "withheld"
        unassessed = bool(cl.get("domains_unavailable")) and te in ("none", "unclassified")
        headline = ("Not assessed — domain scan unavailable" if unassessed else
                    "Terminal direct-repeat candidate" if held and te == "repeat/structural-only" else _te_headline(te))
        chip = (("NOT ASSESSED" if unassessed else "CLASS & ORDER WITHHELD" if te in ("unclassified", "none")
                 else "ORDER WITHHELD") if held or unassessed else None)
        family = ("no call made — the protein-domain scan did not run" if unassessed and te == "none"
                  else str(cl.get("superfamily") or ""))
        if family.lower().startswith(headline.lower()) or family == "—":
            family = ""
        comp = cl.get("completeness") or {}
        status, _ = self._scan_state(rec, cl)
        _, arch = self._orders(cl, status)
        tier = comp.get("tier") or ""
        if tier and (comp.get("coding") or {}).get("truncated_frame") is True and "intact" in tier.lower():
            tier += (" at the domain level — but a domain-bearing ORF runs off the end of the supplied sequence, so the "
                     "element may be truncated")
        rank = _rank_view(cl)
        return {"headline": headline, "chip": chip, "superfamily": family,
                "confidence": _call_view(rec)["confidence"], "rank_line": _rank_line(rank) if rank else "",
                "why": _explanation_view(rec) if cl.get("explanation") else "",
                "completeness": (f"{tier} · domain order {arch} · not detected: "
                                 f"{', '.join(comp.get('missing') or []) or 'none'}") if tier else "",
                "coding": _coding_view(rec) if comp.get("coding", {}).get("note") else "",
                "cautions": _caution_items(rec) + [_note_view(n, rec) for n in rec.get("notes") or []]}

    def _report_refusal(self):
        recs = self.state.get("records") or []
        i = self.state.get("record_index")
        if self.state.get("analysis_state") == "current" and i is not None and 0 <= i < len(recs):
            return ("The selected record was not analysed, so there is nothing to record. Select another record in "
                    f"panel {PANEL['struct']}.")
        return "No current analysis — run analysis again before recording this assay."

    def _report_context(self):
        """Everything the notebook entry and the report say, as plain data, or None when there is no current analysis of
        a selected, classified record. Only named fields are read: a candidate carries the text it was designed on, and
        that text (which may hold the user's own header) is never copied out."""
        import engine, report
        from datetime import datetime, timezone
        res = self.state.get("analysis")
        recs = self.state.get("records") or []
        i = self.state.get("record_index")
        if self.state.get("analysis_state") != "current" or not isinstance(res, dict) or i is None or not 0 <= i < len(recs):
            return None
        try:
            engine.analysis_result_document(res)
        except ValueError:
            return None
        rec = recs[i]
        if not rec.get("valid", True) or not rec.get("classification") or _call_view(rec)["status"] == "failed":
            return None
        manifest = res.get("provenance") or {}
        sha = sequtil.sha256(rec.get("seq") or "")
        line = report.source_line(manifest.get("input"))
        source = ("Source: " + line) if line else (
            "Source: pasted or uploaded sequence (header replaced by a content-derived ID)"
            if rec.get("identitySource") == "local_input" else "Source: pasted or uploaded sequence")
        comp = rec.get("composition") or {}
        ctx = {"version": APP_VERSION, "generated": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M"),
               "specimen": {"id": _record_public_id(rec), "position": (i + 1, len(recs)),
                            "length": comp.get("length", len(rec.get("seq") or "")),
                            "gc": comp.get("gc"), "sha256": sha, "source": source},
               "call": self._call_display(rec)}
        ctx.update(self._report_assay(rec, sha, manifest))
        ctx["_rec"] = rec
        return ctx

    def _report_assay(self, rec, sha, analysis_manifest):
        """Pairs, in-silico PCR, genome scan, conditions and seals for the record `sha` identifies. A pair or PCR lane
        is included only when its sealed input hash equals the record's, so a result computed on another record,
        flank or sub-region never reaches this record's report."""
        import report
        feats = report.feature_spans(rec)

        def bound(c):
            return ((c.get("design_provenance") or {}).get("input") or {}).get("sha256") == sha

        cands = self.state.get("candidates") or []
        designed = [c for c in cands if bound(c)]
        last = self.state.get("lastPcr") or {}
        tested = {(p["fwd"], p["rev"]): k for k, p in enumerate(last.get("pairs") or [])}
        staged = [c for c in self.state.get("pcrPairs") or [] if c.get("custom") or bound(c)]
        chosen, seen, note = [], set(), None
        for c in staged + [c for c in designed if (c["left_seq"], c["right_seq"]) in tested]:
            if (c["left_seq"], c["right_seq"]) not in seen:
                seen.add((c["left_seq"], c["right_seq"])); chosen.append(c)
        if not chosen and designed:
            chosen, note = designed, "No pair has been sent to in-silico PCR yet; every designed pair is listed."
        pcr_man = {}
        for run in self._prov_list():
            for m in run["manifests"]:
                if m.get("runType") == "in-silico-pcr":
                    p = m.get("parameters") or {}
                    pcr_man[(p.get("fwd"), p.get("rev"))] = m

        def man_bound(m):
            return ((m.get("parameters") or {}).get("assay_input") or m.get("input") or {}).get("sha256") == sha

        statuses = last.get("statuses") or []
        words = {"complete": "complete", "complete_zero_products": "complete, no product found",
                 "partial_capped": "stopped early (search cap)", "failed": "failed — not assessed",
                 "not_run": "not run"}
        pairs, lines, lanes, runs, stopped = [], [], [], [], False
        for c in chosen:
            key = (c["left_seq"], c["right_seq"])
            m, k = pcr_man.get(key), tested.get(key)
            if k is not None and m is not None and not man_bound(m):
                k = None                                      # that PCR searched another sequence
            if k is not None and k >= len(last.get("lanes") or []):
                k = None
            label = last["pairs"][k]["label"] if k is not None else self._pair_label(c)
            d = {"label": label, "fwd": c["left_seq"], "rev": c["right_seq"]}
            if bound(c):
                lp, rp = c.get("left_pos"), c.get("right_pos")
                span = (lp[0], rp[1]) if lp and rp else None
                d.update(fwd_span=tuple(lp) if lp else None, rev_span=tuple(rp) if rp else None,
                         product=c.get("product_size"), amplicon_span=span, spans=report.overlaps(span, feats),
                         tm=(c.get("left_tm"), c.get("right_tm")), gc=(c.get("left_gc"), c.get("right_gc")),
                         qc=(c.get("qc") or {}).get("worst") or "n/a")
            else:
                d["note"] = "your own pair (no design record)"
            if k is None:
                d["note"] = (d["note"] + " · " if d.get("note") else "") + "not yet tested in silico"
            pairs.append(d)
            if k is None:
                continue
            lane = last["lanes"][k]
            status = (statuses[k]["status"] if k < len(statuses) else lane.get("status")) or "complete"
            prods = [f"{a.get('product_length', a.get('length')):,} bp ({_amp_call(a, lane.get('has_locus', True))})"
                     for a in lane.get("amplicons") or []]
            stopped = stopped or status == "partial_capped"
            lines.append(f"{label}: {words.get(status, status)}"
                         + (f" · {len(prods)} product(s): " + "; ".join(prods) if prods else ""))
            lanes.append(lane)
            if m is not None:
                runs.append((f"In-silico PCR {label}", m.get("manifestSha256")))
        if stopped:
            lines.append("In-silico PCR stopped early (search cap): the product list is a lower bound.")
        excluded = len(cands) - len(designed)
        scan = None
        gs = self.state.get("lastGenomeScan") or {}
        gp = gs.get("pair") or {}
        hit = next((d for d in pairs if (d["fwd"], d["rev"]) == (gp.get("fwd"), gp.get("rev"))), None)
        if gs and hit:
            crit, asm = gs.get("criteria") or {}, gs.get("assembly") or {}
            scan = [f"{hit['label']} on {gs.get('organism') or 'the scanned genome'} ({asm.get('accession') or '—'}): "
                    f"{(gs.get('summary') or {}).get('verdict') or 'scan complete'}",
                    f"isPcr search: at least {report.plain(crit.get('min_perfect'))} bp perfect 3′ match and "
                    f"{report.plain(crit.get('min_good'))} bp good match, products up to {report.plain(crit.get('max_size'))} bp. "
                    "Sites outside that rule are not reported; this is a prediction from the assembly, not a measured result."]
        conditions = []
        base = next((c for c in chosen if bound(c)), None)
        if base:
            p = (base.get("design_provenance") or {}).get("parameters") or {}
            tgt = p.get("primer3_target")
            conditions.append(
                f"Design: Tm {report.plain(p.get('min_tm'))}–{report.plain(p.get('max_tm'))} °C (opt {report.plain(p.get('opt_tm'))}) · "
                f"length {report.plain(p.get('min_size'))}–{report.plain(p.get('max_size'))} nt (opt {report.plain(p.get('opt_size'))}) · "
                f"GC {report.plain(p.get('min_gc'))}–{report.plain(p.get('max_gc'))} % · "
                f"product {report.plain(p.get('prod_min'))}–{report.plain(p.get('prod_max'))} bp · target "
                + (f"start {tgt[0]}, length {tgt[1]} nt (0-based)" if isinstance(tgt, (list, tuple)) and len(tgt) == 2
                   else "none (untargeted design)"))
            tm = report.tm_text(p.get("tm_conditions") or {})
            conditions.append("Primer3 Tm salt and oligo conditions: " +
                              (tm if tm else "Primer3 defaults (not recorded in this run's seal)"))
            qc = (base.get("qc") or {}).get("conditions") or {}
            if qc:
                warn = sorted({abs(w) for w in (oligoqc._HAIRPIN_WARN, oligoqc._DIMER_WARN, oligoqc._END_WARN)})
                conditions.append(
                    f"Structure QC: ΔG at {report.plain(qc.get('temp_c'))} °C, {report.plain(qc.get('mv_conc'))} mM Na⁺, "
                    f"{report.plain(qc.get('dv_conc'))} mM Mg²⁺, {report.plain(qc.get('dna_conc'))} nM oligo; warn at ΔG ≤ "
                    f"−{' / −'.join(f'{w:g}' for w in warn)} kcal/mol (advisory)")
        mm = next((pcr_man[(d["fwd"], d["rev"])] for d in pairs if (d["fwd"], d["rev"]) in pcr_man and
                   man_bound(pcr_man[(d["fwd"], d["rev"])])), None)
        if mm:
            p = mm.get("parameters") or {}
            conditions.append(f"In-silico PCR search: up to {report.plain(p.get('max_mm'))} mismatches per primer, "
                              f"{report.plain(p.get('tp'))} exact 3′ bases, products {report.plain(p.get('prod_min'))}–"
                              f"{report.plain(p.get('prod_max'))} bp")
        prim = (self.state.get("lastPrimers") or {}).get("provenance") or {}
        seals = [("Analysis", analysis_manifest.get("manifestSha256"))]
        if pairs and (prim.get("input") or {}).get("sha256") == sha:
            seals.append(("Primer design", prim.get("manifestSha256")))
        seals += runs
        if scan:
            gm = next((m for run in reversed(self._prov_list()) for m in run["manifests"]
                       if m.get("runType") == "genome-scan"), None)
            if gm:
                seals.append(("Genome scan", gm.get("manifestSha256")))
        panel = next((d for d in analysis_manifest.get("databases") or [] if "pfam" in str(d.get("name", "")).lower()), None)
        profiles, _why = report.panel_profiles(analysis_manifest)
        return {"pairs": pairs, "excluded": excluded, "pairs_note": note, "pcr": {"lines": lines}, "scan": scan,
                "conditions": conditions, "_lanes": lanes,
                "provenance": {"runs": [(a, b) for a, b in seals if b],
                               "software": "; ".join(f"{s.get('name')} {s.get('version')}"
                                                     for s in analysis_manifest.get("software") or []),
                               "panel": ("Pfam TE-domain panel: " + (f"{len(profiles)} models, " if profiles else "")
                                         + f"SHA-256 {panel['sha256']}") if panel and panel.get("sha256") else None}}

    def _copy_notebook_entry(self):
        import report
        ctx = self._report_context()
        if ctx is None:
            return self._banner(self._report_refusal(), "info")
        from PySide6.QtCore import QMimeData
        data = QMimeData()
        data.setText(report.notebook_text(ctx)); data.setHtml(report.report_html(ctx))
        QApplication.clipboard().setMimeData(data)
        self._banner("Notebook entry copied. Paste it into your notebook or a document.", "success")

    def _assay_report_menu(self, anchor):
        from PySide6.QtWidgets import QMenu
        menu = QMenu(self)
        menu.addAction("PDF — print or attach (figures as images)", lambda: self._export_assay_report("pdf"))
        menu.addAction("HTML — self-contained page (vector figures)", lambda: self._export_assay_report("html"))
        menu.exec(anchor.mapToGlobal(anchor.rect().bottomLeft()))

    def _export_assay_report(self, fmt):
        ctx = self._report_context()
        if ctx is None:
            return self._banner(self._report_refusal(), "info")
        path, _ = QFileDialog.getSaveFileName(
            self, "Save assay report", f"TEagle_{_record_export_stem(ctx['_rec'], 'record')}_assay.{fmt}",
            f"{fmt.upper()} (*.{fmt})")
        if not path:
            return
        if not path.lower().endswith("." + fmt):
            path += "." + fmt
        self._write_assay_report(path, fmt)

    def _write_assay_report(self, path, fmt):
        """Write the assay record as 'pdf' (figures as images) or 'html' (vector figures). Returns True when written; a
        missing or stale analysis writes nothing."""
        import report
        ctx = self._report_context()
        if ctx is None:
            self._banner(self._report_refusal(), "info")
            return False
        rec, lanes = ctx["_rec"], ctx["_lanes"]
        svgs = {}
        model = figures.gv_tracks_from_rec(rec)
        if model["tracks"]:
            svgs["map"] = svg_genome(model, {"start": 0, "end": model["length"]}, 920, "white", for_export=True, legend=True)
        if lanes:
            svgs["gel"] = figures.svg_gel({"lanes": lanes}, "white")

        def write(p):
            if fmt == "pdf":
                images, tags = {}, {}
                for name, svg in svgs.items():
                    images[f"teagle://{name}"] = report.svg_image(svg)
                    width = min(680, int(widgets._svg_size(svg)[0]))
                    tags[name] = f'<img src="teagle://{name}" width="{width}">'
                report.write_document_pdf(report.report_html(ctx, tags, qt=True), images, p)
            else:
                with widgets.atomic_write(p, "w", encoding="utf-8", newline="\n") as fh:
                    fh.write(report.report_html(ctx, svgs))
        return widgets.safe_export(self, path, write, ok=f"Assay report written to {os.path.basename(path)}.") is not None

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

    def closeEvent(self, e):
        jobs = self._running_jobs()
        if jobs and not self._confirm_close(jobs):
            e.ignore()
            return
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
    scientific stack (pyhmmer, primer3), ships the HMM profiles and TLS roots, renders figures through
    QtSvg, and runs an end-to-end analysis — the checks a double-click launch cannot report. Exit 0/1."""
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
    # the frozen launcher bypasses the startup environment check, so a dist-info folder or a module missing
    # from the bundle would only surface as a crash mid-analysis; fail the release gate here instead
    if getattr(sys, "frozen", False):
        try:
            import envcheck
            env = envcheck.check()
            if not env.get("packages_ok"):
                bad = [p["name"] for p in env.get("packages", []) if not p.get("ok")] + list(env.get("importFailures") or [])
                problems.append("bundled package check failed (missing metadata or import): " + ", ".join(sorted(set(bad))))
        except Exception as e:
            problems.append(f"bundled package check crashed: {type(e).__name__}: {e}")
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
    # Network retrieval is intentionally not exercised here, but its trust store must be functional. The
    # frozen launcher bypasses envcheck, and package metadata alone cannot prove that certifi's bundled CA
    # file is present or that OpenSSL can parse it into the context every HTTPS fetch uses.
    try:
        import certifi
        from teagle_core import fetch
        fetch._SSL_CTX_CACHE.clear()
        ca_path = certifi.where()
        tls_context = fetch._ssl_ctx()
        if not os.path.isfile(ca_path):
            problems.append("certifi CA bundle missing from build")
        elif tls_context is None or not tls_context.get_ca_certs():
            problems.append("certifi CA bundle did not load into the TLS context")
    except Exception as e:
        problems.append(f"TLS trust-store self-check failed: {type(e).__name__}: {e}")
    # Vector PDF export relies on the shared SVG -> QPdfWriter path and metadata minimisation. Exercise
    # the actual runtime writer so a frozen build cannot pass merely because the dependencies import.
    try:
        import tempfile
        from pypdf import PdfReader
        from teagle_core.artifact_metadata import media_metadata_violations
        with tempfile.TemporaryDirectory(prefix="teagle-selftest-") as temporary:
            path = os.path.join(temporary, "figure-export.pdf")
            widgets.render_pdf(svg, path)
            violations = media_metadata_violations(path, author=None)   # end-user exports carry no author
            reader = PdfReader(path)
            if violations:
                problems.append("PDF figure export metadata check failed: " + "; ".join(violations))
            if len(reader.pages) != 1:
                problems.append("PDF figure export did not preserve its page")
    except Exception as e:
        problems.append(f"PDF figure export self-check crashed: {type(e).__name__}: {e}")
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
          f"· TLS roots ok · QtSvg ok · install dialog ok")
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
