"""Score TEagle's raw corpus output against the literature-derived labels. Reads only files that were
written by executing the tool; scores nothing that was not run.

Two levels are scored separately, because collapsing them would hide the behaviour this paper is about:

  CLASS   - Class I retrotransposon vs Class II DNA transposon vs not-a-TE. A call at this level is made
            whenever any diagnostic evidence is found.
  ORDER   - LTR / LINE / SINE / TIR. TEagle frequently declines here while still making a class call
            ('retro/partial'), and an abstention is NOT scored as an error. It is counted in its own
            column, so accuracy and abstention can be read together. A tool that abstains often will show
            high accuracy on a small denominator, and hiding that would be a misrepresentation.
  SUPERFAMILY - Copia / Gypsy / hAT / Tc1-Mariner / CACTA / MULE / piggyBac, the rank a bench user
            actually asks about. Scored on the same abstention-aware footing as the two ranks above; see
            the block of tables and commentary below, which carries the judgement calls in full.

Cases whose analysed input was a whole containing record are stratified out of the primary accuracy
figures: analysing a 160 kb contig when the corpus names a 4 kb element inside it tests something else.
They are reported separately, as their own stratum, rather than dropped.

    python benchmarks/score.py                                   # broader corpus
    python benchmarks/score.py --corpus benchmarks/corpus_holdout.tsv \
        --rawdir benchmarks/raw/teagle_holdout --out benchmarks/raw/scores_holdout.json

Output: benchmarks/raw/scores.json
"""
from __future__ import annotations
import argparse, csv, glob, json, math, os, re, sys
from collections import Counter, defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CORPUS = os.path.join(ROOT, "benchmarks", "corpus.tsv")
RAWDIR = os.path.join(ROOT, "benchmarks", "raw", "teagle")
OUT = os.path.join(ROOT, "benchmarks", "raw", "scores.json")

# The two permitted values of the corpus `record_scope` column. `element` means the deposit IS the
# transposable element (or the negative-control gene) essentially in full; `containing_record` means the
# deposit is a clone, contig, chromosome, assembly, vector or genomic region that carries the element
# among other sequence.
ELEMENT, CONTAINING = "element", "containing_record"
SCOPES = (ELEMENT, CONTAINING)

# A corpus row is identified by these three columns together. Not by accession alone: several rows share
# AF391808 and several more share AF123535, one deposited record per group. Not by row order either, because a
# scorer that joined on position would silently re-label every later case the first time a row was added
# or removed. This triple is checked for uniqueness at load, so an ambiguous corpus fails before anything
# is scored.
CASE_KEY = ("accession", "coords", "expected_superfamily")


def case_key(row: dict) -> tuple:
    return tuple((row.get(k) or "").strip() for k in CASE_KEY)


def load_record_scope(corpus_path: str):
    """Map each corpus case to its curated record scope, and give the panel list the corpus defines.

    The scope is read from the corpus rather than from the raw result files, because it is a curation
    decision about the deposit and not an observation made at run time; the raw files predate the column.

    Nothing here defaults. The predecessor of this function inferred scope from the deposit title with two
    keyword lists combined by OR, which could only ever re-admit a row and never exclude one, and which
    fell back on a coordinate test that was vacuously true for any row carrying no coordinate -- so the
    stratification it documented never once fired. A missing or unrecognised value is now a hard failure.
    """
    with open(corpus_path, encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        raise ValueError(f"{corpus_path}: no rows")
    if "record_scope" not in rows[0]:
        raise ValueError(f"{corpus_path}: no record_scope column; every case must declare whether the "
                         f"deposit is the element ({ELEMENT}) or carries it ({CONTAINING})")
    scope, panels = {}, []
    for n, r in enumerate(rows, 1):
        v = (r.get("record_scope") or "").strip()
        if v not in SCOPES:
            raise ValueError(f"{corpus_path} row {n} ({r['accession']}): record_scope is {v!r}, "
                             f"which is not one of {SCOPES}")
        k = case_key(r)
        if k in scope:
            raise ValueError(f"{corpus_path} row {n}: {k} is not unique, so its result cannot be "
                             f"attributed to one ground truth")
        scope[k] = v
        if r["panel"] not in panels:
            panels.append(r["panel"])
    return scope, panels


def wilson(k, n, z=1.959963985):
    if n == 0:
        return (None, None, None)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return round(p, 4), round(max(0.0, c - h), 4), round(min(1.0, c + h), 4)


def expected_class(s: str) -> str:
    t = (s or "").lower()
    if t.startswith("not a te") or "falls outside the wicker" in t:
        return "NOT_TE"
    if "class ii" in t:
        return "CLASS_II"
    if "class i" in t:
        return "CLASS_I"
    return "OTHER"


def expected_order(s: str) -> str:
    t = (s or "").lower()
    if t.startswith("not a te") or "falls outside the wicker" in t:
        return "NOT_TE"
    if "sine" in t:
        return "SINE"
    if "line" in t or "non-ltr" in t:
        return "LINE"
    if "ltr" in t:
        return "LTR"
    if "tir" in t or "class ii" in t:
        return "TIR"
    return "OTHER"


def observed(cl: dict):
    """Map TEagle's output onto the same vocabulary. Returns (class, order, abstained_at_order)."""
    te = (cl.get("te_class") or "").strip()
    klass = (cl.get("class") or "").strip().lower()
    if not te or te.lower() in ("unclassified", "none", "-"):
        return ("NO_CALL", "NO_CALL", True)
    head = te.split("/")[0].strip().lower()
    tail = te.split("/")[1].strip().lower() if "/" in te else ""

    if head in ("ltr", "line", "sine"):
        order = head.upper()
    elif head in ("tir", "dna", "mite", "helitron"):
        order = "TIR"
    elif head == "retro":
        order = "ABSTAIN"                     # class called, order withheld ('retro/partial')
    elif head == "repeat":
        # 'repeat/structural-only' is the tool saying, in its own words, "terminal inverted repeat, class
        # unassigned". That is an abstention, and scoring it as a wrong answer charged the tool for a call
        # it explicitly declined to make -- the opposite of what an abstention-aware benchmark should do.
        order = "NO_CALL"
    else:
        order = "OTHER"

    if "class ii" in klass or order == "TIR":
        c = "CLASS_II"
    elif "class i" in klass or order in ("LTR", "LINE", "SINE", "ABSTAIN"):
        c = "CLASS_I"
    else:
        c = "NO_CALL"
    # Abstention at ORDER level means the order itself was withheld. An earlier version also treated the
    # tail token as an abstention, so "LTR/unclassified" and "LTR/partial" were scored as declining to
    # answer -- but both DO answer the order question (LTR) and withhold only the superfamily beneath it.
    # That conflated two different levels and understated how often the tool commits to an order.
    abstain = order in ("ABSTAIN", "NO_CALL")
    return (c, order, abstain)


# =====================================================================================================
# SUPERFAMILY
# =====================================================================================================
# Class and order are the ranks at which this tool is strongest, and reporting only those two reports the
# levels at which a known error cannot register. AF104899 (Gmr1) is labelled Gypsy by its describing
# authors and carries a documented Copia-like pol order, so the integrase-vs-RT rule MUST call it Copia;
# at order level both labels are LTR and the error is invisible. This block makes it visible.
#
# Everything below is an EXPLICIT table keyed on the exact literal label. It is deliberately not a
# substring or fuzzy match: a reviewer has to be able to disagree with one line of it. Any label that is
# not in the table raises rather than silently becoming "other", so adding a new superfamily to
# classify.py or a new row to a corpus breaks scoring instead of quietly leaving a case ungraded.

ABSTAIN = "ABSTAIN"          # the tool declined to name a superfamily (see SF_OBSERVED)
NG_NO_LABEL = "no superfamily label in the corpus"
NG_RETROVIRUS = "Wicker Retrovirus/ERV superfamily: not in the classifier's answer space"
NOT_GRADABLE = (NG_NO_LABEL, NG_RETROVIRUS)

# ---- What TEagle can say -----------------------------------------------------------------------------
# The classifier's closed answer space at superfamily rank, read off app/backend/teagle_core/classify.py.
# Committed calls map to a token; every other label is the tool explicitly declining at this rank.
#
# WHY "LINE (non-LTR)" IS AN ABSTENTION AND NOT A WRONG ANSWER. Against an expected L1 / CR1 / Jockey the
# tool has emitted an ORDER name into the superfamily field (classify.py:208). It has named the rank above
# and declined the one below - the same act as "LTR retrotransposon (superfamily undetermined)", which is
# already an abstention here. Scoring it wrong would charge the tool for a call it explicitly refused to
# make, which is the opposite of what an abstention-aware benchmark does, and it is inconsistent with how
# the manuscript treats declining everywhere else: an abstention is counted in its own column and never as
# an error. The cost of this choice is paid in the abstention column, where it belongs - these cases stay
# in the denominator, so "TEagle never names a LINE superfamily" shows up as a rate rather than vanishing.
SF_OBSERVED = {
    # committed superfamily calls
    "Copia (Ty1)": "COPIA",                                      # classify.py:148
    "Gypsy (Ty3)": "GYPSY",                                      # classify.py:151
    "hAT": "HAT",                                                # classify.py:254
    "Tc1/Mariner": "TCMAR",                                      # classify.py:256
    "CACTA (En/Spm)": "CACTA",                                   # classify.py:258
    "MULE (Mutator)": "MULE",                                    # classify.py:260
    "IS4-like DDE (piggyBac group)": "PIGGYBAC",                 # classify.py:264
    "DDE transposon": "DDE_UNPLACED",                            # classify.py:266
    "DNA transposon": "DNA_UNPLACED",                            # classify.py:268
    "DIRS-group (tyrosine-recombinase retroelement)": "DIRS",    # classify.py:116
    # explicit refusals at superfamily rank
    "LTR retrotransposon (superfamily undetermined)": ABSTAIN,   # classify.py:143/183/228/242
    "LTR retrotransposon (no coding domains detected)": ABSTAIN,  # classify.py:375
    "retrotransposon (partial)": ABSTAIN,                        # classify.py:202
    "LINE (non-LTR)": ABSTAIN,                                   # classify.py:208 - order named, see above
    "non-autonomous TIR element (superfamily undetermined)": ABSTAIN,   # classify.py:394
    "terminal inverted repeat, class unassigned": ABSTAIN,       # classify.py:405
    "no clear TE signature": ABSTAIN,                            # classify.py:411
}

# These labels are composed at run time rather than written out, so they cannot appear in a literal table.
# Each rule is anchored to BOTH ends of the generated string and to the line that generates it; none of
# them is a bare substring test. Every one is a refusal whose variable part names the domains that WERE
# found, which is why that part has to be allowed to vary.
SF_OBSERVED_COMPOSED = (
    ("LTR retroelement fragment (", "; RT/pol not detected)", ABSTAIN),        # classify.py:358
    ("coding fragment (", "; RT/transposase not detected)", ABSTAIN),          # classify.py:362
    ("", " domain(s) detected; no transposable-element assignment", ABSTAIN),  # classify.py:367
    ("Helitron-family helicase present", "superfamily not called", ABSTAIN),   # classify.py:333
)

# classify.py:153 appends this to "Gypsy (Ty3)" when a chromodomain is present. The chromoviruses are a
# Ty3/Gypsy lineage, so the suffix refines the call without changing the superfamily; it is stripped
# before lookup rather than given its own row, which would invite the two to drift apart.
CHROMOVIRUS_SUFFIX = " · chromovirus"

# Tokens the classifier can actually produce. DERIVED from the table above, not restated, so that
# deleting a branch from classify.py and its row here automatically narrows what counts as representable.
REPRESENTABLE = frozenset(v for v in SF_OBSERVED.values() if v != ABSTAIN)


def observed_superfamily(cl: dict) -> str:
    """Map TEagle's superfamily string onto the token vocabulary, or ABSTAIN. Raises on anything else."""
    s = (cl.get("superfamily") or "").strip()
    if not s:
        return ABSTAIN                     # nothing was classified at all; a refusal, not a wrong answer
    if s.endswith(CHROMOVIRUS_SUFFIX):
        s = s[:-len(CHROMOVIRUS_SUFFIX)]
    if s in SF_OBSERVED:
        return SF_OBSERVED[s]
    for pre, suf, tok in SF_OBSERVED_COMPOSED:
        if s.startswith(pre) and s.endswith(suf):
            return tok
    raise ValueError(f"superfamily {s!r} is not in the scoring vocabulary. classify.py emits a label this "
                     f"scorer has never been told how to read, so leaving it ungraded would understate or "
                     f"overstate accuracy silently. Add it to SF_OBSERVED with a comment.")


# ---- What the corpora ask -----------------------------------------------------------------------------
# Curated labels normalised to the same tokens. Where a label carries a family or lineage qualifier
# ("Ty3/gypsy (family Huck)"), the superfamily is the leading token and the qualifier is below the rank
# being scored, so it is dropped - this is a rank projection, not a loosening of the ground truth.
#
# Tokens that are NOT in REPRESENTABLE (L1, CR1, JOCKEY, ...) are still graded. They are the LINE, SINE,
# Bel-Pao, P and PIF/Harbinger superfamilies, for every one of which the tool declines; keeping them in
# the denominator is what turns "cannot name a LINE superfamily" into a reported abstention rate instead
# of a silent omission. check_declines_where_unrepresentable() asserts that the tool really does decline
# on every one of them, so this cannot quietly become a source of unfair wrong answers.
SF_EXPECTED = {
    # --- Ty1/Copia --------------------------------------------------------------------------------
    "Copia (Ty1)": "COPIA",
    "Copia (Ty1-copia)": "COPIA",
    "Copia (Ty1) LTR retrotransposon": "COPIA",
    "Ty1/copia (Pseudoviridae, family copia)": "COPIA",
    "Ty1/copia (family Fourf)": "COPIA",
    "Ty1/copia (family Opie)": "COPIA",
    "Ty1/copia (RIRE1-related lineage)": "COPIA",
    "Ty1/copia (PREM family, element Prem1)": "COPIA",
    "Ty1/copia (PREM-2 family, element Prem2-a)": "COPIA",
    "Ty1/copia (PREM-2 family, element Prem2-b)": "COPIA",
    # --- Ty3/Gypsy --------------------------------------------------------------------------------
    # Errantivirus / Metavirus are ICTV genera of Metaviridae, but these curated labels lead with the
    # Wicker superfamily and give the genus as a qualifier, so the superfamily is stated and gradable.
    # Contrast the ERV rows below, where ONLY viral taxonomy is given and no Wicker superfamily is named.
    "Gypsy (Ty3)": "GYPSY",
    "Gypsy (Ty3-gypsy)": "GYPSY",
    "Tf2 (Ty3/Gypsy) LTR retrotransposon": "GYPSY",
    "Ty3/gypsy (Errantivirus, family 17.6)": "GYPSY",
    "Ty3/gypsy (Errantivirus, family 297)": "GYPSY",
    "Ty3/gypsy (Metavirus, family 412)": "GYPSY",
    "Ty3/gypsy (Cinful-Zeon family, element Zeon1)": "GYPSY",
    "Ty3/gypsy (family Huck)": "GYPSY",
    "Ty3/gypsy-clade (family blood)": "GYPSY",
    # --- Bel-Pao ----------------------------------------------------------------------------------
    # A Wicker LTR superfamily alongside Copia and Gypsy. The classifier has no Bel-Pao branch, so it can
    # never be answered correctly; it is graded anyway because the tool declines on both rows.
    "Pao/BEL (Semotivirus, family roo)": "BELPAO",
    # --- LINE superfamilies (Wicker order LINE) ---------------------------------------------------
    "L1": "L1",
    "CR1": "CR1",
    "Jockey": "JOCKEY",
    "Jockey (non-LTR LINE retrotransposon)": "JOCKEY",
    "Tad1": "TAD1",
    "R1": "R1",
    "R2": "R2",
    "I": "I",
    "RTE": "RTE",
    "RTE/BovB": "RTE",              # BovB is the mammalian RTE lineage, not a superfamily of its own
    # --- SINE superfamilies -----------------------------------------------------------------------
    # Wicker names SINE superfamilies by the RNA the element derives from; Alu and B1 are both 7SL.
    "7SL-RNA-derived (Alu family)": "SINE_7SL",
    "7SL-RNA-derived (B1 family)": "SINE_7SL",
    # --- TIR superfamilies ------------------------------------------------------------------------
    "hAT": "HAT",
    "Tc1/Mariner": "TCMAR",
    "Tc1/Mariner (Tc1-like)": "TCMAR",
    "Tc1 (Tc1/mariner) DNA transposon": "TCMAR",
    # AM493772.1 (Bari1-family MITE, transposase annotated /pseudo). Autonomy and superfamily are separate
    # questions: a MITE is a deletion derivative OF a superfamily, and the corpus names which one. Grading
    # it at superfamily rank is therefore legitimate; whether the tool should have called it autonomous is
    # a different axis and is not what this block measures.
    "Tc1/Mariner-related MITE": "TCMAR",
    "CACTA": "CACTA",
    "CACTA (En/Spm)": "CACTA",
    "MULE (Mutator)": "MULE",
    "MULE (Mutator, atypical/FHY3-FAR1-related)": "MULE",
    "MULE-related MITE": "MULE",                 # same reasoning as the Bari1 MITE above
    "IS4-like/piggyBac": "PIGGYBAC",
    "PIF/Harbinger": "PIFHARB",                  # Wicker superfamily; no branch in classify.py
    # X06779.1 / X06590.1. "generic DDE (P-element family)" is read at SUPERFAMILY rank as P: the
    # parenthetical names the Wicker superfamily and "generic DDE" describes the transposase model that
    # would have to match it, not a coarser ground truth. Mapping it to the classifier's own
    # "DDE transposon" bin instead would grade the corpus against the tool's vocabulary rather than
    # against the literature, which is backwards.
    "generic DDE (P-element family)": "P_ELEMENT",
    "P-element family": "P_ELEMENT",
    # --- Not gradable: no superfamily in the ground truth ------------------------------------------
    # Blank covers the negative controls (no superfamily exists for a host gene or a satellite) and four
    # elements the curators deliberately left unresolved at this rank. With no ground truth there is
    # nothing to be right or wrong about, so they leave the denominator entirely.
    "": NG_NO_LABEL,
    # AB305073 (MusD). The corpus note says in as many words that the superfamily was left unresolved.
    # That is ground truth withheld by the curator - a different thing from the ERV rows below, where the
    # ground truth is stated but the classifier has no bin for it. Both leave the denominator; they are
    # counted under separate reasons so the two are not conflated.
    "MusD element, ERVL-related (superfamily left unresolved, see note)": NG_NO_LABEL,
    # --- Not gradable: Retrovirus / ERV ------------------------------------------------------------
    # THE ONE CARVE-OUT, AND THE ARGUMENT FOR IT.
    #
    # Wicker et al. 2007 place Retrovirus and ERV as superfamilies of the LTR order ALONGSIDE Copia,
    # Gypsy and Bel-Pao - not beneath Gypsy. classify.py's LTR branch can emit exactly three things:
    # Copia, Gypsy, or "superfamily undetermined" (classify.py:143-155). There is NO Retrovirus or ERV
    # token anywhere in its output vocabulary, and REPRESENTABLE above is derived from that vocabulary,
    # so it contains neither. The Wicker Retrovirus superfamily is therefore NOT representable by this
    # classifier at all.
    #
    # Given that, "Gypsy (Ty3)" on an ERV is:
    #   - not CORRECT. Retroviral RTs do group with the Ty3/Gypsy clade and ERVs share the Gypsy pol
    #     order, which is exactly why the integrase-vs-RT rule fires - but a phylogenetic affinity is not
    #     a superfamily assignment. The manuscript states class and superfamily follow Wicker, and under
    #     Wicker these are sibling superfamilies. Scoring it correct would silently redefine the rank.
    #   - not an ABSTENTION. Unlike "LINE (non-LTR)", "Gypsy (Ty3)" is not the rank above with the rank
    #     below withheld; it is a positive, committed, and false superfamily call.
    #   - so either WRONG or NOT GRADABLE. Counting it wrong would report the absence of a bin as a
    #     discrimination failure, and would pool it with the two genuine Copia/Gypsy inversions, which
    #     are a different defect. These rows are therefore excluded from the denominator, counted under
    #     their own reason, and the accuracy that WOULD have followed from grading them incorrect is
    #     reported alongside as a sensitivity figure, so the exclusion cannot flatter the result unseen.
    #
    # The exclusion is decided on the EXPECTED label - the whole lineage leaves together. Roughly half
    # these rows draw an abstention rather than a Gypsy call, and splitting the panel on what the tool
    # happened to emit would put one ground-truth category on both sides of the denominator boundary.
    "HML-2": NG_RETROVIRUS,
    "HML-2 (LTR5Hs)": NG_RETROVIRUS,
    "HERV-K (HML-2), Betaretrovirus-like / ERVK": NG_RETROVIRUS,
    "HERV-K (HML-2) LTR retrotransposon / endogenous retrovirus": NG_RETROVIRUS,
    "HERV-H, Gammaretrovirus-like / ERV1": NG_RETROVIRUS,
    "HERV-W, Gammaretrovirus-like / ERV1": NG_RETROVIRUS,
    "IAP element, Betaretrovirus-like / ERVK": NG_RETROVIRUS,
    "Gammaretrovirus-like / ERV1 (KoRV, GALV-related)": NG_RETROVIRUS,
    "Gammaretrovirus-like / ERV1 (PERV lineage)": NG_RETROVIRUS,
    "Gammaretrovirus-like / ERV1 (RD-114/BaEV lineage)": NG_RETROVIRUS,
    "Gammaretrovirus-like / ERV1 (murine leukemia virus lineage)": NG_RETROVIRUS,
    "Betaretrovirus-like / ERVK (JSRV/enJSRV lineage)": NG_RETROVIRUS,
    "Betaretrovirus-like / ERVK (MMTV lineage)": NG_RETROVIRUS,
    "Betaretrovirus-like / ERVK (type-D lineage, Mason-Pfizer/BaEV-related)": NG_RETROVIRUS,
    "Retroviridae, Alpharetrovirus (endogenous ALV subgroup E; chicken ev-3 locus)": NG_RETROVIRUS,
    "Retroviridae, Betaretrovirus (type D; enJSRV lineage, clone enJS56A1)": NG_RETROVIRUS,
    "Retroviridae, Epsilonretrovirus (RMERV)": NG_RETROVIRUS,
    "Retroviridae, Gammaretrovirus (type C; PERV-A)": NG_RETROVIRUS,
    "Retroviridae, gammaretrovirus (MLV)-related clade; genus deliberately left unassigned by the "
    "describing authors (ZFERV)": NG_RETROVIRUS,
    "Retroviridae, genus unassignable per the describing authors (PyERV) - NOT betaretrovirus, "
    "NOT gammaretrovirus": NG_RETROVIRUS,
}


def expected_superfamily(s: str) -> str:
    """Normalise a curated superfamily label. Raises on anything the table does not name."""
    t = (s or "").strip()
    if t in SF_EXPECTED:
        return SF_EXPECTED[t]
    raise ValueError(f"expected_superfamily {t!r} is not in SF_EXPECTED. A corpus row carries a label this "
                     f"scorer cannot place, and guessing its rank would put an unreviewed mapping into a "
                     f"published accuracy figure. Add it with a comment saying which Wicker superfamily "
                     f"it names, or mark it not gradable.")


def check_declines_where_unrepresentable(cases):
    """Assert the premise the LINE/Bel-Pao/P/PIF grading decision rests on.

    Those cases are kept in the denominator on the grounds that the tool DECLINES on them, so keeping them
    costs it abstention rather than accuracy. That is an empirical claim about this corpus and this build,
    not a law, so it is checked rather than assumed: if a future build ever answers "Copia" for a Bel-Pao
    element or "Tc1/Mariner" for a PIF/Harbinger, the case would be scored wrong for a distinction the
    classifier was never built to make - exactly the unfairness the Retrovirus carve-out exists to avoid -
    and the scorer must stop and force that decision to be made explicitly rather than absorb it.
    """
    bad = [c for c in cases if c["superfamily_gradable"]
           and c["expected_superfamily_token"] not in REPRESENTABLE
           and not c["abstained_at_superfamily"]]
    if bad:
        d = "; ".join(f"{c['accession']} expected {c['expected_superfamily_token']} "
                      f"-> called {c['observed_superfamily']!r}" for c in bad)
        raise ValueError(
            f"{len(bad)} case(s) whose expected superfamily is outside the classifier's answer space drew "
            f"a committed superfamily call: {d}. These would be scored incorrect for a superfamily the "
            f"tool cannot express. Decide explicitly whether they join the Retrovirus carve-out.")


def acc_superfamily(sub):
    """Superfamily accuracy, in the same shape as acc() at class and order.

    n_gradable counts cases with a curated superfamily the classifier could in principle be asked for;
    n_answered counts those where it committed to one. Not-gradable cases are outside n_gradable entirely
    and are reported with their reason, because an excluded row that is not counted is an excluded row
    nobody can audit.
    """
    gradable = [c for c in sub if c["superfamily_gradable"]]
    answered = [c for c in gradable if not c["abstained_at_superfamily"]]
    correct = [c for c in answered
               if c["expected_superfamily_token"] == c["observed_superfamily_token"]]
    p, lo, hi = wilson(len(correct), len(answered))
    return {"n_gradable": len(gradable), "n_answered": len(answered), "n_correct": len(correct),
            "accuracy": p, "ci95_low": lo, "ci95_high": hi,
            "abstention_rate": round(1 - len(answered) / len(gradable), 4) if gradable else None,
            "n_not_gradable": len(sub) - len(gradable),
            "not_gradable_by_reason": dict(Counter(c["superfamily_not_gradable_reason"] for c in sub
                                                   if not c["superfamily_gradable"]))}


def acc_superfamily_erv_as_incorrect(sub):
    """Sensitivity: the same figures if the Retrovirus/ERV rows were graded rather than excluded.

    Exists so the carve-out cannot flatter the headline unseen. An ERV row that drew a committed call
    joins the answered denominator as an error; one that drew an abstention joins as an abstention, which
    is how every other declined case is already treated.
    """
    sub2 = [dict(c) for c in sub]
    for c in sub2:
        if c["superfamily_not_gradable_reason"] == NG_RETROVIRUS:
            c["superfamily_gradable"] = True
            c["expected_superfamily_token"] = "RETROVIRUS"      # not in REPRESENTABLE -> never matches
    return acc_superfamily(sub2)


def acc_superfamily_representable_only(sub):
    """Sensitivity in the other direction: restricted to cases the classifier could answer correctly.

    Answers the opposite objection to the one above - that the abstention rate is inflated by rows where
    declining was the tool's only honest option. Accuracy is unchanged by construction; the abstention
    rate is the number that moves.
    """
    return acc_superfamily([c for c in sub if not c["superfamily_gradable"]
                            or c["expected_superfamily_token"] in REPRESENTABLE])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--corpus", default=CORPUS, help="corpus TSV carrying the record_scope column")
    ap.add_argument("--rawdir", default=RAWDIR, help="directory of raw per-case results")
    ap.add_argument("--out", default=OUT, help="scores JSON to write")
    args = ap.parse_args()

    files = sorted(f for f in glob.glob(os.path.join(args.rawdir, "*.json"))
                   if not os.path.basename(f).startswith("_"))
    if not files:
        print(f"no raw output in {args.rawdir} - run benchmarks/run_teagle.py first")
        return 1
    record_scope, panels = load_record_scope(args.corpus)

    cases, stratified = [], []
    for f in files:
        d = json.load(open(f, encoding="utf-8"))
        row = d["corpus_row"]
        recs = d["result"].get("records") or [{}]
        cl = recs[0].get("classification") or {}
        applied = (d.get("fetch_source") or {}).get("coords_applied")
        oc, oo, ab = observed(cl)
        esf = expected_superfamily(row.get("expected_superfamily"))
        osf = observed_superfamily(cl)
        entry = {
            "accession": d["accession"], "panel": row.get("panel"), "organism": row.get("organism"),
            "expected_class": expected_class(row.get("expected_class")),
            "expected_order": expected_order(row.get("expected_class")),
            "expected_superfamily": row.get("expected_superfamily"),
            "observed_class": oc, "observed_order": oo, "observed_te_class": cl.get("te_class"),
            "observed_superfamily": cl.get("superfamily"), "confidence": cl.get("confidence"),
            # Raw labels are kept above; these are the normalised tokens the superfamily block scores on,
            # carried per case so any figure in the tables can be traced back to the row that produced it.
            "expected_superfamily_token": None if esf in NOT_GRADABLE else esf,
            "observed_superfamily_token": osf,
            "superfamily_gradable": esf not in NOT_GRADABLE,
            "superfamily_not_gradable_reason": esf if esf in NOT_GRADABLE else None,
            "abstained_at_superfamily": osf == ABSTAIN,
            "completeness": (cl.get("completeness") or {}).get("tier")
            if isinstance(cl.get("completeness"), dict) else cl.get("completeness"),
            "abstained_at_order": ab, "orfs_unscanned": cl.get("orfs_unscanned"),
            "n_domains": cl.get("n_domains"), "input_length": d.get("input_length"),
            # Carried so the count of distinct inputs can be derived rather than assumed equal to the
            # count of cases: some rows of this corpus are several ground truths over one analysed record.
            "input_sha256": d.get("input_sha256"),
            "coords_applied": applied, "corpus_confidence": row.get("confidence"),
            "citation": row.get("citation"),
        }
        # Is the ANALYSED input one element? Two things decide it, and they are different things.
        #
        # The first is the curated scope of the deposit, read from the corpus. It is curated rather than
        # inferred because no mechanical test of the DEFINITION line survives contact with the records:
        # Drosophila P1 clones are titled "... DNA sequence (P1s ...), complete sequence" and one
        # Ty3 deposit leads with the tRNA-Cys gene it sits beside, so a keyword rule either misses the
        # clones or throws out the element. Deposit length fails for the same reason in the other
        # direction -- element size and record size overlap across 0.4 to 20 kb.
        key = case_key(row)
        if key not in record_scope:
            raise ValueError(f"{os.path.basename(f)}: {key} is not a case in {args.corpus}; the corpus and "
                             f"the raw results are out of step, so no result can be attributed safely")
        entry["record_scope"] = record_scope[key]

        # The second is whether a coordinate was actually applied. A containing record narrowed to a plain
        # span before analysis IS a single-element input: the engine never saw the rest of the record.
        # Negative controls are commonly deposited as a chromosome or a whole bacterial genome that
        # supplies the control feature as a sub-range, and stratifying those out on deposit scope alone
        # would delete the cases most likely to expose a false positive -- the wrong direction for a
        # scorer to err in. run_teagle.py records the applied span as "<start>-<end>"; anything else
        # ("NOT APPLICABLE ...", "OUT OF RANGE ...") means the whole record was analysed.
        span = re.fullmatch(r"(\d+)-(\d+)", applied or "")
        if span:
            expect = int(span.group(2)) - int(span.group(1)) + 1
            if entry["input_length"] != expect:
                raise ValueError(f"{os.path.basename(f)}: span {applied} is {expect} bp but "
                                 f"{entry['input_length']} bp were analysed")
        entry["span_applied"] = bool(span)
        entry["single_element_input"] = entry["record_scope"] == ELEMENT or entry["span_applied"]
        if not entry["single_element_input"]:
            entry["excluded_because"] = (
                f"deposit is a containing record and no coordinate was applied, so the whole "
                f"{entry['input_length']} bp record was analysed"
                + (f" ({applied})" if applied else ""))
        entry["deposit_title"] = str((d.get("fetch_source") or {}).get("title") or "")[:120]
        (cases if entry["single_element_input"] else stratified).append(entry)

    # The premise behind grading the non-representable superfamilies is checked over EVERY case, retained
    # and stratified alike, before any figure is computed from them.
    check_declines_where_unrepresentable(cases + stratified)

    def acc(sub, level):
        # Superfamily is scored by its own function rather than through the branch below: its gradability
        # is a per-case decision from an explicit table, not the single "OTHER" sentinel that class and
        # order use, and folding it in would have meant editing the path those two ranks run on.
        if level == "superfamily":
            return acc_superfamily(sub)
        exp, obs = f"expected_{level}", f"observed_{level}"
        gradable = [c for c in sub if c[exp] not in ("OTHER",)]
        if level == "order":
            answered = [c for c in gradable if not c["abstained_at_order"]]
        else:
            answered = [c for c in gradable if c[obs] != "NO_CALL"]
        correct = [c for c in answered if c[exp] == c[obs]]
        p, lo, hi = wilson(len(correct), len(answered))
        return {"n_gradable": len(gradable), "n_answered": len(answered), "n_correct": len(correct),
                "accuracy": p, "ci95_low": lo, "ci95_high": hi,
                "abstention_rate": round(1 - len(answered) / len(gradable), 4) if gradable else None}

    result = {
        "cases_scored": len(cases),
        "cases_excluded_not_single_element": len(stratified),
        "overall": {"class": acc(cases, "class"), "order": acc(cases, "order"),
                    "superfamily": acc(cases, "superfamily")},
        # Both sensitivity views of the superfamily denominator, so neither exclusion argument has to be
        # taken on trust. Neither is a table column; they exist to be quoted against the headline.
        "superfamily_sensitivity": {
            "erv_graded_incorrect": acc_superfamily_erv_as_incorrect(cases),
            "representable_expected_only": acc_superfamily_representable_only(cases),
        },
        # Every case the superfamily block scores as an error, named. Derived from the same tokens the
        # accuracy figure uses, so the list and the count cannot disagree.
        "superfamily_incorrect": [
            {"accession": c["accession"], "panel": c["panel"], "organism": c["organism"],
             "expected": c["expected_superfamily"], "observed": c["observed_superfamily"],
             "expected_token": c["expected_superfamily_token"],
             "observed_token": c["observed_superfamily_token"], "confidence": c["confidence"]}
            for c in cases if c["superfamily_gradable"] and not c["abstained_at_superfamily"]
            and c["expected_superfamily_token"] != c["observed_superfamily_token"]],
        # Every case excluded from the superfamily denominator, named, with the reason. An exclusion that
        # is only a count is an exclusion nobody can check.
        "superfamily_not_gradable": [
            {"accession": c["accession"], "panel": c["panel"], "expected": c["expected_superfamily"],
             "observed": c["observed_superfamily"], "reason": c["superfamily_not_gradable_reason"]}
            for c in cases if not c["superfamily_gradable"]],
        # Panels come from the CORPUS, not from the retained cases. A panel every one of whose rows is a
        # containing record -- refusal-supply is exactly that, three views of one maize adh1 contig --
        # would otherwise disappear from this table and from Table 2 without leaving a trace, which is the
        # same silent omission this stratifier exists to end. It now reports n = 0 and its stratified count.
        "by_panel": {p: {"n": len([c for c in cases if c["panel"] == p]),
                         "n_stratified": len([c for c in stratified if c["panel"] == p]),
                         "class": acc([c for c in cases if c["panel"] == p], "class"),
                         "order": acc([c for c in cases if c["panel"] == p], "order"),
                         "superfamily": acc([c for c in cases if c["panel"] == p], "superfamily")}
                     for p in sorted(panels)},
        "by_confidence_tier": {t: {"n": len([c for c in cases if c["confidence"] == t]),
                                   "class": acc([c for c in cases if c["confidence"] == t], "class"),
                                   "order": acc([c for c in cases if c["confidence"] == t], "order"),
                                   "superfamily": acc([c for c in cases if c["confidence"] == t],
                                                      "superfamily")}
                               for t in sorted({str(c["confidence"]) for c in cases})},
        "negative_controls": {
            "n": len([c for c in cases if c["expected_class"] == "NOT_TE"]),
            "correctly_not_called": len([c for c in cases if c["expected_class"] == "NOT_TE"
                                         and c["observed_class"] == "NO_CALL"]),
            "false_positives": [{"accession": c["accession"], "called": c["observed_te_class"],
                                 "confidence": c["confidence"]}
                                for c in cases if c["expected_class"] == "NOT_TE"
                                and c["observed_class"] != "NO_CALL"],
        },
        "confusion_order": dict(Counter(f"{c['expected_order']}->{c['observed_order']}" for c in cases)),
        # The stratum, reported rather than dropped. `incorrect_calls` is the number that says whether
        # removing these rows flattered the tool: a stratified row that carried a WRONG call was
        # suppressing accuracy, so excluding it raises the headline and that has to be visible.
        "stratum_containing_record": {
            "n": len(stratified),
            "class": acc(stratified, "class"), "order": acc(stratified, "order"),
            "superfamily": acc(stratified, "superfamily"),
            "incorrect_calls": [
                {"accession": c["accession"], "panel": c["panel"], "level": lvl,
                 "expected": c[f"expected_{lvl}"], "observed": c[f"observed_{lvl}"]}
                for c in stratified for lvl in ("class", "order")
                if c[f"expected_{lvl}"] != "OTHER" and c[f"observed_{lvl}"] not in ("NO_CALL", "ABSTAIN")
                and c[f"expected_{lvl}"] != c[f"observed_{lvl}"]],
            # Superfamily is listed separately rather than added to the loop above, which compares the
            # raw expected/observed columns; at this rank the comparison is between normalised tokens.
            "superfamily_incorrect": [
                {"accession": c["accession"], "panel": c["panel"],
                 "expected": c["expected_superfamily"], "observed": c["observed_superfamily"]}
                for c in stratified if c["superfamily_gradable"] and not c["abstained_at_superfamily"]
                and c["expected_superfamily_token"] != c["observed_superfamily_token"]],
            "distinct_inputs": len({c["input_sha256"] for c in stratified}),
        },
        "distinct_inputs_scored": len({c["input_sha256"] for c in cases}),
        "excluded_cases": [{"accession": c["accession"], "panel": c["panel"],
                            "reason": c.get("excluded_because")} for c in stratified],
        "cases": cases,
    }
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    json.dump(result, open(args.out, "w", encoding="utf-8"), indent=1, sort_keys=True)

    o = result["overall"]
    print(f"scored {len(cases)} cases from {result['distinct_inputs_scored']} distinct inputs; "
          f"{len(stratified)} stratified out as containing records\n")
    for lvl in ("class", "order", "superfamily"):
        a = o[lvl]
        print(f"{lvl.upper():11s} answered {a['n_answered']}/{a['n_gradable']}  "
              f"correct {a['n_correct']}  accuracy {a['accuracy']} "
              f"[{a['ci95_low']}, {a['ci95_high']}]  abstention {a['abstention_rate']}"
              + (f"  not gradable {a['n_not_gradable']}" if "n_not_gradable" in a else ""))
    for r, n in sorted(o["superfamily"]["not_gradable_by_reason"].items()):
        print(f"            not gradable: {n:3d}  {r}")
    for k, a in sorted(result["superfamily_sensitivity"].items()):
        print(f"            sensitivity {k}: {a['n_correct']}/{a['n_answered']} of {a['n_gradable']} "
              f"= {a['accuracy']}  abstention {a['abstention_rate']}")
    print("\nby panel:")
    for p, v in result["by_panel"].items():
        print(f"  {p:20s} n={v['n']:3d} (+{v['n_stratified']:2d} stratified)  "
              f"class {v['class']['accuracy']}  order {v['order']['accuracy']} "
              f"(abstained {v['order']['abstention_rate']})  "
              f"superfamily {v['superfamily']['n_correct']}/{v['superfamily']['n_answered']}"
              f" of {v['superfamily']['n_gradable']} = {v['superfamily']['accuracy']} "
              f"(abstained {v['superfamily']['abstention_rate']})")
    print("\nsuperfamily errors:")
    for w in result["superfamily_incorrect"]:
        print(f"   INCORRECT {w['accession']:12s} {w['panel']:20s} expected {w['expected']!r} "
              f"-> called {w['observed']!r} ({w['confidence']})")
    st = result["stratum_containing_record"]
    print(f"\ncontaining-record stratum: {st['n']} cases over {st['distinct_inputs']} distinct inputs; "
          f"class {st['class']['n_correct']}/{st['class']['n_answered']} answered, "
          f"order {st['order']['n_correct']}/{st['order']['n_answered']} answered")
    for w in st["incorrect_calls"]:
        print(f"   INCORRECT (stratified) {w['accession']} {w['level']}: "
              f"expected {w['expected']}, called {w['observed']}")
    nc = result["negative_controls"]
    print(f"\nnegative controls: {nc['correctly_not_called']}/{nc['n']} correctly not called")
    for fp in nc["false_positives"]:
        print(f"   FALSE POSITIVE {fp['accession']}: {fp['called']} ({fp['confidence']})")
    print(f"\n-> {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
