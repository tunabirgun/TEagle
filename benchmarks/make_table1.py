"""Table 1 - classification accuracy on the independent validation corpus, from
benchmarks/raw/scores_holdout.json.

The manuscript cited a table of these numbers that was never produced, while a second table carried the
same number. This builds the missing one. Its structure matches Table 2 so the two corpora can be read
against each other without re-learning the columns.

    python benchmarks/make_table1.py
"""
from __future__ import annotations
import json, os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [os.path.join(ROOT, "benchmarks")]
from make_tables import write_csv, write_xlsx, write_docx           # noqa: E402

SRC = os.path.join(ROOT, "benchmarks", "raw", "scores_holdout.json")
OUT = os.path.join(ROOT, "manuscript", "tables")

PANEL_LABEL = {
    "ltr-copia": "LTR retrotransposons (Ty1/Copia)",
    "ltr-gypsy": "LTR retrotransposons (Ty3/Gypsy)",
    "erv": "Endogenous retroviruses",
    "line": "LINEs (non-LTR retrotransposons)",
    "tir": "DNA transposons (TIR)",
    "negative": "Negative controls (non-TE)",
}


def fmt(a):
    if a.get("accuracy") is None:
        return "—", "—", "—"
    return (f"{a['accuracy']:.3f}", f"{a['ci95_low']:.3f}–{a['ci95_high']:.3f}",
            f"{a['abstention_rate']:.3f}")


def main():
    if not os.path.exists(SRC):
        print("missing scores_holdout.json - run the validation scoring first")
        return 1
    d = json.load(open(SRC, encoding="utf-8"))

    header = ["Panel", "n", "Class: answered/gradable", "Class accuracy", "Class 95% CI",
              "Class abstention", "Order: answered/gradable", "Order accuracy", "Order 95% CI",
              "Order abstention"]
    rows = []
    for panel, v in sorted(d["by_panel"].items(), key=lambda kv: -kv[1]["n"]):
        c, o = v["class"], v["order"]
        ca, cci, cab = fmt(c)
        oa, oci, oab = fmt(o)
        rows.append([PANEL_LABEL.get(panel, panel), v["n"],
                     f"{c['n_answered']}/{c['n_gradable']}", ca, cci, cab,
                     f"{o['n_answered']}/{o['n_gradable']}", oa, oci, oab])
    ov = d["overall"]
    ca, cci, cab = fmt(ov["class"])
    oa, oci, oab = fmt(ov["order"])
    # Element-only totals, because pooling the negative controls with the elements mixes two questions:
    # on an element an abstention is a withheld answer, on a negative control it is the right answer.
    el = [c for c in d["cases"] if c["expected_class"] != "NOT_TE"]
    el_ans = [c for c in el if not c["abstained_at_order"] and c["observed_order"] != "NO_CALL"]
    el_ok = [c for c in el_ans if c["expected_order"] == c["observed_order"]]
    rows.append(["Transposable elements only", len(el), "—", "—", "—", "—",
                 f"{len(el_ans)}/{len(el)}",
                 f"{len(el_ok) / len(el_ans):.3f}" if el_ans else "—", "—",
                 f"{1 - len(el_ans) / len(el):.3f}" if el else "—"])
    rows.append(["All panels combined", len(d["cases"]),
                 f"{ov['class']['n_answered']}/{ov['class']['n_gradable']}", ca, cci, cab,
                 f"{ov['order']['n_answered']}/{ov['order']['n_gradable']}", oa, oci, oab])

    nc = d["negative_controls"]
    fp = nc["false_positives"]
    caption = (
        "Table 1. Classification accuracy on the independent validation corpus. Each case carries a "
        "class and superfamily label taken from the publication that described the element, and each "
        "record is the element itself rather than a clone or contig containing it, so no coordinate "
        "arithmetic intervenes between the deposit and the analysis. A class call distinguishes Class I "
        "from Class II from not-a-transposable-element; an order call names LTR, LINE or TIR and, below "
        "it, the superfamily. Abstentions are counted in their own column rather than scored as errors, "
        "because a tool that abstains freely would otherwise post a high accuracy on a shrinking "
        "denominator; accuracy and abstention must therefore be read together. Intervals are Wilson 95%. "
        "The two summary rows are given separately because declining is the correct answer on a negative "
        "control and a withheld answer on a real element. Negative controls: "
        f"{nc['correctly_not_called']} of {nc['n']} drew no call of any kind"
        + (f"; {fp[0]['accession']} was assigned a class at {fp[0]['confidence']} confidence with the "
           "order withheld." if fp else "."))

    os.makedirs(OUT, exist_ok=True)
    write_csv(os.path.join(OUT, "table1_validation_accuracy.csv"), header, rows)
    write_xlsx(os.path.join(OUT, "table1_validation_accuracy.xlsx"), header, rows, "Table 1")
    write_docx(os.path.join(OUT, "table1_validation_accuracy.docx"), header, rows, caption)
    print(f"Table 1 : {len(rows)} rows x {len(header)} cols -> .csv .xlsx .docx")
    for r in rows:
        print(f"  {str(r[0])[:34]:36s} n={r[1]:<4} class {r[3]:>6}  order {r[7]:>6}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
