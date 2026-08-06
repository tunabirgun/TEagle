"""Table 2 - classification accuracy by panel, built from benchmarks/raw/scores.json.

Accuracy and abstention are reported side by side for every panel, because either alone is misleading: a
tool that abstains freely posts a high accuracy on a shrinking denominator, and a tool that always answers
posts a low one while being more useful. The denominator for each figure is printed rather than implied.

    python benchmarks/make_table2.py
"""
from __future__ import annotations
import json, os, sys

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


def main():
    if not os.path.exists(SRC):
        print("missing scores.json - run benchmarks/score.py first")
        return 1
    d = json.load(open(SRC, encoding="utf-8"))

    header = ["Panel", "n", "Class: answered/gradable", "Class accuracy", "Class 95% CI",
              "Class abstention", "Order: answered/gradable", "Order accuracy", "Order 95% CI",
              "Order abstention", "Superfamily: answered/gradable", "Superfamily accuracy",
              "Superfamily 95% CI", "Superfamily abstention"]
    rows = []
    for panel, v in sorted(d["by_panel"].items(), key=lambda kv: -kv[1]["n"]):
        c, o, s = v["class"], v["order"], v["superfamily"]
        ca, cci, cab = fmt(c)
        oa, oci, oab = fmt(o)
        sa, sci, sab = fmt(s)
        rows.append([PANEL_LABEL.get(panel, panel), v["n"],
                     f"{c['n_answered']}/{c['n_gradable']}", ca, cci, cab,
                     f"{o['n_answered']}/{o['n_gradable']}", oa, oci, oab,
                     f"{s['n_answered']}/{s['n_gradable']}", sa, sci, sab])
    ov = d["overall"]
    ca, cci, cab = fmt(ov["class"])
    oa, oci, oab = fmt(ov["order"])
    sa, sci, sab = fmt(ov["superfamily"])
    rows.append(["All panels combined", d["cases_scored"],
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
               f"Intervals are Wilson 95%. {sf['n_not_gradable']} cases lie outside the superfamily "
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
               f"the analysed input ({lost}); they are reported as their own stratum. See "
               f"benchmarks/CORPUS_SCOPE.md.")

    os.makedirs(OUT, exist_ok=True)
    write_csv(os.path.join(OUT, "table2_classification_accuracy.csv"), header, rows)
    write_xlsx(os.path.join(OUT, "table2_classification_accuracy.xlsx"), header, rows, "Table 2")
    write_docx(os.path.join(OUT, "table2_classification_accuracy.docx"), header, rows, caption)
    print(f"Table 2 : {len(rows)} rows x {len(header)} cols -> .csv .xlsx .docx")
    for r in rows:
        print(f"  {str(r[0])[:34]:36s} n={r[1]:<4} class {r[3]:>6}  order {r[7]:>6}  "
              f"superfamily {r[10]:>7} {r[11]:>6}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
