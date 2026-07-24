"""Primer / oligo secondary-structure QC — hairpin, self-dimer (homodimer), hetero-dimer, and 3'-end
stability — the way IDT OligoAnalyzer reports them, computed from real nearest-neighbor thermodynamics.

Two independent engines cross-check each other (the user's requirement: Primer3 alone can disagree with IDT):
  - PRIMARY: primer3-py calc_hairpin / calc_homodimer / calc_heterodimer / calc_end_stability
    (libprimer3 `thal`, SantaLucia 1998 NN params). dg is returned in cal/mol -> divided by 1000 here.
  - CROSS-CHECK: ViennaRNA (RNAfold / RNAcofold) with the DNA Mathews-2004 parameters — an INDEPENDENT
    minimum-free-energy implementation. For dimers the comparable quantity is the BINDING free energy
    (ΔG of the AB duplex minus each strand's own MFE), not the raw cofold MFE.

Neither is an IDT clone: IDT OligoAnalyzer runs the mfold/UNAFold engine, which is academic-non-commercial
and not redistributable. When the optional UNAFold backend component is installed (WSL, bioconda), the
engine module `oligoqc_unafold` adds that IDT-matching column; this module stays pure/in-process (no WSL,
no network) so it is unit-testable and always available.

Thresholds follow IDT's published rule of thumb: a structure with ΔG ≤ -9 kcal/mol on ANY axis is flagged
warn (SantaLucia 1998; Owczarzy 2008). The caution boundary varies by structure type — hairpin ≤ -2,
dimers ≤ -5, and 3'-end ≤ -6. The 3'-end axis measures the ΔG of one primer's 3' end annealing to the
OTHER primer (a 3'-end cross-dimer), NOT the primer's stability against its intended template; its ~-6
cutoff flags only an abnormally stable 3'-end cross-dimer. It is shown separately because a base-paired
3' end is the one geometry a polymerase can extend, so it is the single most decisive dimer-failure
signal. The tiering is a TEagle design choice anchored on those sources, not a verbatim IDT table.
"""
from __future__ import annotations

try:                                                 # primer3 is already a hard dep of the primer module;
    import primer3                                   # a broken build must degrade, not crash the engine
    _P3_ERR = None
except Exception as _e:                              # pragma: no cover
    primer3 = None
    _P3_ERR = f"{type(_e).__name__}: {_e}"

try:                                                 # ViennaRNA is optional: absent -> primer3-only, cross-check "n/a"
    import RNA
    _VRNA_VER = RNA.__version__
    _VRNA_ERR = None
except Exception as _e:                              # pragma: no cover
    RNA = None
    _VRNA_VER = None
    _VRNA_ERR = f"{type(_e).__name__}: {_e}"

PRIMER3_VERSION = getattr(primer3, "__version__", "unavailable") if primer3 is not None else "unavailable"
VIENNARNA_VERSION = _VRNA_VER or "unavailable"

# IDT OligoAnalyzer-comparable reaction conditions (documented + sealed). Mg2+=0 and 0.25 uM oligo match
# IDT's OligoAnalyzer defaults; ViennaRNA has no divalent model, so both engines run at matched monovalent
# salt + temperature so their numbers are compared like-for-like.
CONDITIONS = {"mv_conc": 50.0, "dv_conc": 0.0, "dntp_conc": 0.0, "dna_conc": 250.0, "temp_c": 37.0}

# ΔG (kcal/mol) tiers — TEagle heuristic anchored on IDT's -9 kcal/mol figure + the 3'-end/internal asymmetry.
_HAIRPIN_WARN, _HAIRPIN_CAUTION = -9.0, -2.0
_DIMER_WARN, _DIMER_CAUTION = -9.0, -5.0
_END_WARN, _END_CAUTION = -9.0, -6.0            # 3'-end cross-dimer ΔG (last bases annealing to the other primer) — its own axis
_AGREE_BAND = 2.0                                # two independent NN engines within ~1-2 kcal/mol => concordant


def available() -> dict:
    return {"primer3": primer3 is not None, "primer3_version": PRIMER3_VERSION,
            "viennarna": RNA is not None, "viennarna_version": VIENNARNA_VERSION,
            "primer3_error": _P3_ERR, "viennarna_error": _VRNA_ERR}


def _round(x, n=2):
    return round(x, n) if isinstance(x, (int, float)) else None


# ---------------- primer3 (primary engine) ----------------
def _p3_kcal(tr):
    """A primer3 ThermoResult -> (dg_kcal_per_mol or None if no structure, tm)."""
    if tr is None or not getattr(tr, "structure_found", False):
        return None, None
    return tr.dg / 1000.0, tr.tm                     # dg is cal/mol -> kcal/mol


def _p3_hairpin(seq):
    return _p3_kcal(primer3.calc_hairpin(seq, **_p3_args()))


def _p3_homodimer(seq):
    return _p3_kcal(primer3.calc_homodimer(seq, **_p3_args()))


def _p3_heterodimer(a, b):
    return _p3_kcal(primer3.calc_heterodimer(a, b, **_p3_args()))


def _p3_end_stability(a, b):
    # gate on structure_found (like _p3_kcal): no favorable 3'-end cross-dimer -> None -> renders "—", not 0.0
    tr = primer3.calc_end_stability(a, b, **_p3_args())
    return (tr.dg / 1000.0 if (tr is not None and getattr(tr, "structure_found", False)) else None)


def _p3_args():
    c = CONDITIONS
    return {"mv_conc": c["mv_conc"], "dv_conc": c["dv_conc"], "dntp_conc": c["dntp_conc"],
            "dna_conc": c["dna_conc"], "temp_c": c["temp_c"]}


# ---------------- ViennaRNA (independent cross-check) ----------------
_VRNA_READY = False


def _vrna_md():
    """A ViennaRNA model_details set to DNA Mathews-2004 params at matched T + monovalent salt.
    NB: params_load_DNA_Mathews2004() mutates ViennaRNA's PROCESS-GLOBAL energy table irreversibly (there is no
    load-default counterpart called here). oligoqc is currently the ONLY module that folds RNA, so there is no
    collision; if a future feature folds RNA under default (RNA) parameters in the same process, load the params
    it needs explicitly before folding rather than assuming defaults. Loaded once (idempotent) per process."""
    global _VRNA_READY
    if not _VRNA_READY:
        RNA.params_load_DNA_Mathews2004()
        _VRNA_READY = True
    md = RNA.md()
    md.temperature = CONDITIONS["temp_c"]
    try:
        md.salt = CONDITIONS["mv_conc"] / 1000.0     # mol/L monovalent (ViennaRNA >=2.6); ignored if unsupported
    except Exception:                                # pragma: no cover
        pass
    return md


def _vrna_mfe(seq):
    fc = RNA.fold_compound(seq, _vrna_md())
    return fc.mfe()[1]


def _vrna_hairpin(seq):
    if RNA is None:
        return None
    dg = _vrna_mfe(seq)
    return dg if dg < 0 else None                    # >=0 -> no favourable fold


def _vrna_binding(a, b):
    """Intermolecular binding ΔG = MFE(A&B duplex) - MFE(A) - MFE(B). This is the quantity comparable to
    primer3's calc_homodimer/heterodimer (which report the duplex, not each strand's own hairpin)."""
    if RNA is None:
        return None
    md = _vrna_md()
    ab = RNA.fold_compound(a + "&" + b, md).mfe()[1]
    ga = RNA.fold_compound(a, md).mfe()[1]
    gb = RNA.fold_compound(b, md).mfe()[1]
    dg = ab - ga - gb
    return dg if dg < 0 else None


# ---------------- flag logic ----------------
def _tier(dg, warn, caution):
    """ok / caution / warn from a ΔG (more negative = more stable = worse). None -> ok (no structure)."""
    if dg is None:
        return "ok"
    if dg <= warn:
        return "warn"
    if dg <= caution:
        return "caution"
    return "ok"


def _agree(p3, vrna):
    """Concordance label between the two engines. Both None -> 'none'. One value missing: 'n/a' when the
    ViennaRNA cross-check engine is not installed (couldn't run), else 'single' (a genuine one-engine result,
    e.g. ViennaRNA predicted no structure). Both present -> agree/disagree on the ±band."""
    if vrna is None and p3 is None:
        return "none"
    if p3 is None or vrna is None:
        return "n/a" if RNA is None else "single"
    return "agree" if abs(p3 - vrna) <= _AGREE_BAND else "disagree"


def _metric(p3_dg, vrna_dg, warn, caution):
    tier = _tier(p3_dg, warn, caution)
    vt = _tier(vrna_dg, warn, caution)
    # worst-of the two engines drives the shown flag, so a concerning call from EITHER is not hidden
    order = {"ok": 0, "caution": 1, "warn": 2}
    flag = tier if order[tier] >= order[vt] else vt
    return {"p3": _round(p3_dg), "vrna": _round(vrna_dg), "flag": flag, "agree": _agree(p3_dg, vrna_dg)}


# ---------------- sequence-only descriptors ----------------
def gc_clamp(seq: str, n: int = 5) -> int:
    """Count G/C in the last n (3') bases — a 1-3 clamp stabilises the 3' end, >3 risks mispriming."""
    tail = seq.upper()[-n:]
    return sum(1 for b in tail if b in "GC")


def longest_poly_x(seq: str) -> int:
    seq = seq.upper()
    best = run = 0
    prev = ""
    for b in seq:
        run = run + 1 if b == prev else 1
        prev = b
        best = max(best, run)
    return best


def _oligo_tm(seq):
    if primer3 is None:
        return None
    try:
        return _round(primer3.calc_tm(seq, mv_conc=CONDITIONS["mv_conc"], dv_conc=CONDITIONS["dv_conc"],
                                      dntp_conc=CONDITIONS["dntp_conc"], dna_conc=CONDITIONS["dna_conc"]), 1)
    except Exception:                                # pragma: no cover
        return None


# ---------------- public API ----------------
def qc_oligo(seq: str) -> dict:
    """Per-oligo QC: hairpin + self-dimer ΔG (both engines) + sequence descriptors."""
    seq = (seq or "").upper()
    if primer3 is None or not seq:
        return {"ok": False, "error": _P3_ERR or "empty oligo"}
    hp_p3, _ = _p3_hairpin(seq)
    hd_p3, _ = _p3_homodimer(seq)
    return {
        "ok": True, "seq": seq, "len": len(seq),
        "tm": _oligo_tm(seq), "gc": _round(100.0 * sum(1 for b in seq if b in "GC") / len(seq), 1),
        "gc_clamp": gc_clamp(seq), "poly_x": longest_poly_x(seq),
        "hairpin": _metric(hp_p3, _vrna_hairpin(seq), _HAIRPIN_WARN, _HAIRPIN_CAUTION),
        "self_dimer": _metric(hd_p3, _vrna_binding(seq, seq), _DIMER_WARN, _DIMER_CAUTION),
    }


def qc_pair(left: str, right: str) -> dict:
    """Full primer-pair QC: each primer's hairpin/self-dimer, the cross (hetero) dimer, and 3'-end stability."""
    left, right = (left or "").upper(), (right or "").upper()
    if primer3 is None:
        return {"ok": False, "error": _P3_ERR or "primer3 unavailable"}
    if not left or not right:                            # guard empty primers (mirrors qc_oligo), never call calc_* on ""
        return {"ok": False, "error": "empty primer sequence"}
    het_p3, _ = _p3_heterodimer(left, right)
    het = _metric(het_p3, _vrna_binding(left, right), _DIMER_WARN, _DIMER_CAUTION)
    # 3'-end stability is a primer3-only metric (last-bases anneal ΔG); worst of F-on-R / R-on-F
    ends = [d for d in (_p3_end_stability(left, right), _p3_end_stability(right, left)) if d is not None]
    end_dg = min(ends) if ends else None
    end = {"p3": _round(end_dg), "vrna": None, "flag": _tier(end_dg, _END_WARN, _END_CAUTION), "agree": "single"}
    left_qc, right_qc = qc_oligo(left), qc_oligo(right)
    worst = _worst_flag([left_qc["hairpin"]["flag"], left_qc["self_dimer"]["flag"],
                         right_qc["hairpin"]["flag"], right_qc["self_dimer"]["flag"], het["flag"], end["flag"]])
    return {"ok": True, "left": left_qc, "right": right_qc, "hetero_dimer": het, "end_stability": end,
            "worst": worst, "conditions": dict(CONDITIONS),
            "engines": {"primer3": PRIMER3_VERSION, "viennarna": VIENNARNA_VERSION}}


def _worst_flag(flags):
    order = {"ok": 0, "caution": 1, "warn": 2}
    return max(flags, key=lambda f: order.get(f, 0)) if flags else "ok"
