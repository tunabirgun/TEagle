"""Score TEagle's raw corpus output against the literature-derived labels. Reads only files that were
written by executing the tool; scores nothing that was not run.

Two levels are scored separately, because collapsing them would hide the behaviour this paper is about:

  CLASS   - Class I retrotransposon vs Class II DNA transposon vs not-a-TE. A call at this level is made
            whenever any diagnostic evidence is found.
  ORDER   - LTR / LINE / SINE / TIR, and below it the superfamily. TEagle frequently declines here while
            still making a class call ('retro/partial'), and an abstention is NOT scored as an error. It
            is counted in its own column, so accuracy and abstention can be read together. A tool that
            abstains often will show high accuracy on a small denominator, and hiding that would be a
            misrepresentation.

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
        entry = {
            "accession": d["accession"], "panel": row.get("panel"), "organism": row.get("organism"),
            "expected_class": expected_class(row.get("expected_class")),
            "expected_order": expected_order(row.get("expected_class")),
            "expected_superfamily": row.get("expected_superfamily"),
            "observed_class": oc, "observed_order": oo, "observed_te_class": cl.get("te_class"),
            "observed_superfamily": cl.get("superfamily"), "confidence": cl.get("confidence"),
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

    def acc(sub, level):
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
        "overall": {"class": acc(cases, "class"), "order": acc(cases, "order")},
        # Panels come from the CORPUS, not from the retained cases. A panel every one of whose rows is a
        # containing record -- refusal-supply is exactly that, three views of one maize adh1 contig --
        # would otherwise disappear from this table and from Table 2 without leaving a trace, which is the
        # same silent omission this stratifier exists to end. It now reports n = 0 and its stratified count.
        "by_panel": {p: {"n": len([c for c in cases if c["panel"] == p]),
                         "n_stratified": len([c for c in stratified if c["panel"] == p]),
                         "class": acc([c for c in cases if c["panel"] == p], "class"),
                         "order": acc([c for c in cases if c["panel"] == p], "order")}
                     for p in sorted(panels)},
        "by_confidence_tier": {t: {"n": len([c for c in cases if c["confidence"] == t]),
                                   "class": acc([c for c in cases if c["confidence"] == t], "class"),
                                   "order": acc([c for c in cases if c["confidence"] == t], "order")}
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
            "incorrect_calls": [
                {"accession": c["accession"], "panel": c["panel"], "level": lvl,
                 "expected": c[f"expected_{lvl}"], "observed": c[f"observed_{lvl}"]}
                for c in stratified for lvl in ("class", "order")
                if c[f"expected_{lvl}"] != "OTHER" and c[f"observed_{lvl}"] not in ("NO_CALL", "ABSTAIN")
                and c[f"expected_{lvl}"] != c[f"observed_{lvl}"]],
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
    for lvl in ("class", "order"):
        a = o[lvl]
        print(f"{lvl.upper():6s} answered {a['n_answered']}/{a['n_gradable']}  "
              f"correct {a['n_correct']}  accuracy {a['accuracy']} "
              f"[{a['ci95_low']}, {a['ci95_high']}]  abstention {a['abstention_rate']}")
    print("\nby panel:")
    for p, v in result["by_panel"].items():
        print(f"  {p:20s} n={v['n']:3d} (+{v['n_stratified']:2d} stratified)  "
              f"class {v['class']['accuracy']}  order {v['order']['accuracy']} "
              f"(abstained {v['order']['abstention_rate']})")
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
