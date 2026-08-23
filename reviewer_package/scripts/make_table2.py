"""Table 2 - classification accuracy by panel, built from benchmarks/raw/scores.json.

Accuracy and abstention are reported side by side for every panel, because either alone is misleading: a
tool that abstains freely posts a high accuracy on a shrinking denominator, and a tool that always answers
posts a low one while being more useful. The denominator for each figure is printed rather than implied.

    python benchmarks/make_table2.py
"""
from __future__ import annotations
import json, math, os, sys
from collections import Counter

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [os.path.join(ROOT, "benchmarks")]
from make_tables import write_csv, write_xlsx, write_docx           # noqa: E402
# Imported rather than restated, so the caption cannot describe an exclusion category under a name the
# scorer has since renamed.
from score import NG_RETROVIRUS, NG_NO_LABEL                        # noqa: E402

SRC = os.path.join(ROOT, "benchmarks", "raw", "scores.json")
OUT = os.path.join(ROOT, "manuscript", "tables")

PANEL_LABEL = {
    "ltr-copia-gypsy": "LTR retrotransposons (Copia / Gypsy)",
    "erv-architecture": "Endogenous retroviruses",
    "line-nonltr": "LINEs (non-LTR retrotransposons)",
    "tir-dna": "DNA transposons (TIR)",
    "divergence-gradient": "LTR divergence gradient",
    "negative-controls": "Negative controls (non-TE)",
    "refusal-supply": "Unreadable domain order",
    "primers-copynumber": "Assay design",
}


def fmt(a):
    # Abstention is reported even when accuracy is not. Blanking it alongside accuracy emptied the
    # column at exactly the rows where the tool abstained on everything — the column the paper spends a
    # paragraph justifying. "n/a" means no gradable case; a number means the rate was computed.
    acc = ("n/a", "n/a") if a.get("accuracy") is None else (
        f"{a['accuracy']:.3f}", f"{a['ci95_low']:.3f}–{a['ci95_high']:.3f}")
    ab = "n/a" if a.get("abstention_rate") is None else f"{a['abstention_rate']:.3f}"
    return acc[0], acc[1], ab


def _wilson(k, n, z=1.959963985):
    if not n:
        return None, None, None
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return round(p, 4), round(max(0.0, c - h), 4), round(min(1.0, c + h), 4)


def _rank_accuracy(cases, level):
    """Independently recompute one stored score block from case records."""
    if level == "superfamily":
        gradable = [c for c in cases if c["superfamily_gradable"]]
        answered = [c for c in gradable if not c["abstained_at_superfamily"]]
        correct = [c for c in answered
                   if c["expected_superfamily_token"] == c["observed_superfamily_token"]]
        p, lo, hi = _wilson(len(correct), len(answered))
        return {"n_gradable": len(gradable), "n_answered": len(answered), "n_correct": len(correct),
                "accuracy": p, "ci95_low": lo, "ci95_high": hi,
                "abstention_rate": round(1 - len(answered) / len(gradable), 4) if gradable else None,
                "n_not_gradable": len(cases) - len(gradable),
                "not_gradable_by_reason": dict(Counter(c["superfamily_not_gradable_reason"] for c in cases
                                                        if not c["superfamily_gradable"]))}
    expected, observed = f"expected_{level}", f"observed_{level}"
    gradable = [c for c in cases if c[expected] != "OTHER"]
    answered = ([c for c in gradable if not c["abstained_at_order"]]
                if level == "order" else [c for c in gradable if c[observed] != "NO_CALL"])
    correct = [c for c in answered if c[expected] == c[observed]]
    p, lo, hi = _wilson(len(correct), len(answered))
    return {"n_gradable": len(gradable), "n_answered": len(answered), "n_correct": len(correct),
            "accuracy": p, "ci95_low": lo, "ci95_high": hi,
            "abstention_rate": round(1 - len(answered) / len(gradable), 4) if gradable else None}


def _collapse_unique_inputs(cases, level):
    fields = {
        "class": ("expected_class", "observed_class"),
        "order": ("expected_order", "observed_order", "abstained_at_order"),
        "superfamily": ("expected_superfamily_token", "observed_superfamily_token",
                        "superfamily_gradable", "superfamily_not_gradable_reason",
                        "abstained_at_superfamily"),
    }[level]
    grouped = {}
    for case in cases:
        sha = (case.get("input_sha256") or "").strip()
        if not sha:
            raise ValueError(f"scores.json case {case.get('accession')} lacks input_sha256")
        grouped.setdefault(sha, []).append(case)
    representatives = []
    for sha, rows in sorted(grouped.items()):
        signatures = {tuple(row.get(field) for field in fields) for row in rows}
        if len(signatures) != 1:
            raise ValueError(f"scores.json input_sha256 {sha} has conflicting {level} rows")
        representatives.append(rows[0])
    return representatives


def assert_score_integrity(d):
    """Reject stale or hand-edited aggregates before Table 2 can reproduce them."""
    cases = d["cases"]
    stratified = d.get("stratified_cases")
    if stratified is None:
        raise ValueError("scores.json lacks stratified_cases; regenerate it with score.py")
    if d.get("cases_scored") != len(cases):
        raise ValueError("scores.json cases_scored disagrees with cases")
    for level in ("class", "order", "superfamily"):
        actual = _rank_accuracy(cases, level)
        if d["overall"].get(level) != actual:
            raise ValueError(f"scores.json overall {level} disagrees with its case records")
        unique_actual = _rank_accuracy(_collapse_unique_inputs(cases, level), level)
        if (d.get("overall_unique_inputs") or {}).get(level) != unique_actual:
            raise ValueError(f"scores.json overall_unique_inputs {level} disagrees with its case records")
    distinct_inputs = len(_collapse_unique_inputs(cases, "class"))
    if d.get("distinct_inputs_scored") != distinct_inputs:
        raise ValueError("scores.json distinct_inputs_scored disagrees with case hashes")
    element_cases = [c for c in cases if c["expected_class"] != "NOT_TE"
                     and c["expected_order"] not in (None, "", "OTHER")]
    element_unique = _collapse_unique_inputs(element_cases, "order")
    withheld = [c for c in element_unique
                if c["abstained_at_order"] or c["observed_order"] == "NO_CALL"]
    element_actual = {"n_gradable": len(element_unique), "n_withheld": len(withheld),
                      "n_answered": len(element_unique) - len(withheld),
                      "withholding_rate": round(len(withheld) / len(element_unique), 4)
                                          if element_unique else None}
    stored_element = ((d.get("independence") or {}).get("element_only_order_withholding"))
    if stored_element != element_actual:
        raise ValueError("scores.json element-only unique-input withholding disagrees with case records")
    for panel, stored in d["by_panel"].items():
        panel_cases = [c for c in cases if c["panel"] == panel]
        if stored.get("n") != len(panel_cases):
            raise ValueError(f"scores.json panel {panel} n disagrees with its case records")
        if stored.get("n_stratified") != sum(c["panel"] == panel for c in stratified):
            raise ValueError(f"scores.json panel {panel} n_stratified disagrees with its case records")
        for level in ("class", "order", "superfamily"):
            if stored.get(level) != _rank_accuracy(panel_cases, level):
                raise ValueError(f"scores.json panel {panel} {level} disagrees with its case records")


def main():
    if not os.path.exists(SRC):
        print("missing scores.json - run benchmarks/score.py first")
        return 1
    d = json.load(open(SRC, encoding="utf-8"))
    assert_score_integrity(d)

    header = ["Panel", "Unit", "n", "Class: answered/gradable", "Class accuracy", "Class 95% CI",
              "Class abstention", "Order: answered/gradable", "Order accuracy", "Order 95% CI",
              "Order abstention", "Superfamily: answered/gradable", "Superfamily accuracy",
              "Superfamily 95% CI", "Superfamily abstention"]
    rows = []
    for panel, v in sorted(d["by_panel"].items(), key=lambda kv: -kv[1]["n"]):
        c, o, s = v["class"], v["order"], v["superfamily"]
        ca, cci, cab = fmt(c)
        oa, oci, oab = fmt(o)
        sa, sci, sab = fmt(s)
        rows.append([PANEL_LABEL.get(panel, panel), "curated case rows", v["n"],
                     f"{c['n_answered']}/{c['n_gradable']}", ca, cci, cab,
                     f"{o['n_answered']}/{o['n_gradable']}", oa, oci, oab,
                     f"{s['n_answered']}/{s['n_gradable']}", sa, sci, sab])
    # The combined headline must not pretend that duplicated annotations are independent observations.
    # Panel rows stay row-level because panels are curation strata; the visibly labelled final row retains
    # the historical row-level summary as a sensitivity analysis.
    ov = d["overall_unique_inputs"]
    ca, cci, cab = fmt(ov["class"])
    oa, oci, oab = fmt(ov["order"])
    sa, sci, sab = fmt(ov["superfamily"])
    rows.append(["All panels combined (primary)", "distinct input hashes", d["distinct_inputs_scored"],
                 f"{ov['class']['n_answered']}/{ov['class']['n_gradable']}", ca, cci, cab,
                 f"{ov['order']['n_answered']}/{ov['order']['n_gradable']}", oa, oci, oab,
                 f"{ov['superfamily']['n_answered']}/{ov['superfamily']['n_gradable']}", sa, sci, sab])
    ov = d["overall"]
    ca, cci, cab = fmt(ov["class"])
    oa, oci, oab = fmt(ov["order"])
    sa, sci, sab = fmt(ov["superfamily"])
    rows.append(["All panels combined (row-level sensitivity)", "curated case rows", d["cases_scored"],
                 f"{ov['class']['n_answered']}/{ov['class']['n_gradable']}", ca, cci, cab,
                 f"{ov['order']['n_answered']}/{ov['order']['n_gradable']}", oa, oci, oab,
                 f"{ov['superfamily']['n_answered']}/{ov['superfamily']['n_gradable']}", sa, sci, sab])

    nc = d["negative_controls"]
    # Which panels lost rows, and how many, is derived from the scores file rather than restated here. A
    # panel can now legitimately show n = 0 -- refusal-supply is three views of one maize contig and all
    # three are containing records -- and a zero with no explanation beside it would be worse than the
    # silent disappearance it replaced.
    lost = ", ".join(f"{PANEL_LABEL.get(p, p)} {v['n_stratified']}"
                     for p, v in sorted(d["by_panel"].items(), key=lambda kv: -kv[1]["n_stratified"])
                     if v["n_stratified"])
    sf = ov["superfamily"]
    ng = sf["not_gradable_by_reason"]
    sens = d["superfamily_sensitivity"]["erv_graded_incorrect"]
    # The superfamily errors are named in the caption rather than counted, because two inversions in
    # opposite directions is the finding this row exists to expose and a bare "2" would hide it.
    wrong = "; ".join(f"{w['accession']} ({w['organism']}), labelled {w['expected']}, called "
                      f"{w['observed']}" for w in d["superfamily_incorrect"])
    caption = (f"Table 2. Classification accuracy by panel. A class call distinguishes Class I from "
               f"Class II from not-a-transposable-element; an order call names LTR, LINE or TIR; "
               f"a superfamily call names Copia, Gypsy, hAT, Tc1/Mariner, CACTA, MULE or piggyBac. "
               f"Abstentions are counted in their own column, not scored as errors, so accuracy and "
               f"abstention must be read together; naming the order while withholding the superfamily "
               f"beneath it, as TEagle does on every LINE, is an abstention at superfamily rank. "
               f"Combined primary figures use distinct analysed input hashes; panel rows and the explicitly "
               f"labelled sensitivity row use curated case rows. Intervals are Wilson 95%. "
               f"{sf['n_not_gradable']} cases lie outside the superfamily "
               f"denominator: {ng.get(NG_NO_LABEL, 0)} carry no superfamily label in the corpus, and "
               f"{ng.get(NG_RETROVIRUS, 0)} are endogenous retroviruses, whose Wicker superfamily "
               f"(Retrovirus / ERV) has no token in this classifier's output vocabulary — it calls them "
               f"Gypsy (Ty3), the superfamily whose pol order they share, so neither verdict would "
               f"measure its discrimination. Had those cases been graded incorrect instead, superfamily "
               f"accuracy would be {sens['accuracy']:.3f} on {sens['n_answered']} answered. The "
               f"{len(d['superfamily_incorrect'])} superfamily errors are inversions in both directions: "
               f"{wrong}. Negative controls: "
               f"{nc['correctly_not_called']} of {nc['n']} correctly returned no call. "
               f"{d['cases_excluded_not_single_element']} cases were stratified out of these figures "
               f"because the deposit carries the element among other sequence and no coordinate narrowed "
               f"the analysed input ({lost}); they are reported as their own stratum. The supplied "
               f"corpus-scope description gives the rule and counts.")

    os.makedirs(OUT, exist_ok=True)
    write_csv(os.path.join(OUT, "table2_classification_accuracy.csv"), header, rows)
    write_xlsx(os.path.join(OUT, "table2_classification_accuracy.xlsx"), header, rows, "Table 2")
    write_docx(os.path.join(OUT, "table2_classification_accuracy.docx"), header, rows, caption)
    print(f"Table 2 : {len(rows)} rows x {len(header)} cols -> .csv .xlsx .docx")
    for r in rows:
        print(f"  {str(r[0])[:34]:36s} n={r[2]:<4} class {r[4]:>6}  order {r[8]:>6}  "
              f"superfamily {r[11]:>7} {r[12]:>6}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
