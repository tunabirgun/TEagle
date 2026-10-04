"""Helpers, constants and shared widgets used by main.py and the panel mixins."""
from __future__ import annotations

import functools
import os
import re
import sys
from PySide6.QtCore import QByteArray, QEasingCurve, QPropertyAnimation, Qt
from PySide6.QtGui import QColor, QIcon, QPainter, QPixmap
from PySide6.QtSvg import QSvgRenderer
from PySide6.QtWidgets import (
    QComboBox, QFrame, QLabel, QPushButton, QSizePolicy, QStyledItemDelegate, QVBoxLayout, QWidget)
import widgets
from widgets import DataTable
import theme as theme_mod
from teagle_core import appdirs, sequtil, oligoqc
from teagle_core.wsl import dfam_coverage_counts as _dfam_coverage_counts
from engine import LARGE_INPUT_BP

_HERE = os.path.dirname(os.path.abspath(__file__))

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

UI_SCALES = [0.75, 0.85, 1.0, 1.1, 1.25, 1.5]     # user-selectable global UI scale (persisted, applied at startup)
