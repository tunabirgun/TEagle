"""Table 3 - predicted against published amplicon size, from benchmarks/raw/assay_scores.json.

Reference numbers are not written here. Each corpus entry carries the DOI of its source, and the
manuscript's reference list carries the same DOIs, so the number is looked up by DOI at build time. A
renumbered bibliography therefore renumbers this table, and a source cited in the table but missing from
the bibliography is an error rather than a wrong-looking bracket.

    python benchmarks/make_table3.py
"""
from __future__ import annotations
import json, os, re, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [os.path.join(ROOT, "benchmarks")]
from make_tables import write_csv, write_xlsx, write_docx           # noqa: E402

SRC = os.path.join(ROOT, "benchmarks", "raw", "assay_scores.json")
MS = os.path.join(ROOT, "manuscript", "manuscript.md")
OUT = os.path.join(ROOT, "manuscript", "tables")

DOI_RE = re.compile(r"10\.\d{4,9}/[^\s;,)\]]+")


def reference_numbers():
    """DOI -> bibliography number, read from the manuscript's own reference list."""
    if not os.path.exists(MS):
        return {}
    text = open(MS, encoding="utf-8").read()
    tail = text.split("## References", 1)[-1]
    out = {}
    for m in re.finditer(r"^(\d+)\.\s+(.*)$", tail, re.M):
        d = DOI_RE.search(m.group(2))
        if d:
            out[d.group(0).rstrip(".").lower()] = int(m.group(1))
    return out


def cite(citation, numbers):
    # A table whose purpose is checkability must not carry a row a reader cannot trace, so a citation
    # with no DOI is a failure rather than an em dash.
    d = DOI_RE.search(citation or "")
    if not d:
        return "NO DOI IN CORPUS"
    n = numbers.get(d.group(0).rstrip(".").lower())
    return f"[{n}]" if n else "NOT IN BIBLIOGRAPHY"


def main():
    if not os.path.exists(SRC):
        print("missing assay_scores.json - run benchmarks/run_assay.py first")
        return 1
    d = json.load(open(SRC, encoding="utf-8"))
    numbers = reference_numbers()

    header = ["Assay", "Target", "Organism", "Template", "Published (bp)", "Predicted (bp)",
              "Difference (bp)", "Products", "Loci", "Size stated in source", "Source"]
    rows = []
    for r in sorted(d["rows"], key=lambda x: (x["panel"] or "", x["name"])):
        diff = r["difference_bp"]
        rows.append([
            # The target is not truncated: a cut cell is indistinguishable from a whole one, and eight
            # of these seventeen ended mid-word at the previous 60-character limit.
            r["name"], (r["element"] or "").split(".")[0], r["organism"],
            r["accession"], r["expected_bp"],
            r["closest_bp"] if r["closest_bp"] is not None else "no product",
            f"{diff:+d}" if diff is not None else "—",
            r["n_products"],
            "single" if r["loci_on_template"] == "single" else "multiple",
            "yes" if r["ground_truth"] == "stated_in_source" else "gel estimate",
            cite(r["citation"], numbers),
        ])

    st = d["by_ground_truth"].get("stated_in_source", {})
    caption = (
        "Table 3. Predicted against published amplicon size for assays whose product size was determined "
        "at the bench. Each pair was run on the template its source names, through the same in-silico PCR "
        "the application exposes, with one size window applied to every case "
        f"({d['search_window_bp'][0]}–{d['search_window_bp'][1]} bp: half the smallest published "
        "amplicon to three times the largest) rather than a window centred on the expected value, which "
        "would let the search decide the answer. Predicted size is the full product, templated span plus "
        "any incorporated non-templated 5′ tail, because that is the band a gel reports. The "
        "'Size stated in source' column separates sizes the source prints as an integer from those it "
        "gives only as a gel estimate; the concordance statistic in the text is taken over the former. "
        "The 'Loci' column marks the one assay that is multi-locus by design, for which a product count "
        "is not a specificity result. Concordance is scored within a tolerance of 10% or 20 bp, whichever "
        "is larger, since a size read off a gel against a ladder carries its own error. "
        f"Of the {d['single_locus_cases']} single-locus assays, {d['single_locus_concordant']} fall within "
        f"tolerance and {d['exact_matches']} reproduce the published size exactly; median absolute "
        f"difference {d['median_abs_difference_bp']} bp, maximum "
        f"{st.get('max_abs_difference_bp', '—')} bp.")

    os.makedirs(OUT, exist_ok=True)
    write_csv(os.path.join(OUT, "table3_assay_validation.csv"), header, rows)
    write_xlsx(os.path.join(OUT, "table3_assay_validation.xlsx"), header, rows, "Table 3")
    write_docx(os.path.join(OUT, "table3_assay_validation.docx"), header, rows, caption)
    print(f"Table 3 : {len(rows)} rows x {len(header)} cols -> .csv .xlsx .docx")
    bad = [r for r in rows if not r[-1].startswith("[")]
    for r in rows:
        print(f"  {str(r[0])[:38]:40s} {r[4]:>6} -> {str(r[5]):>6}  {r[6]:>5}  {r[-1]}")
    if bad:
        print(f"\n{len(bad)} row(s) could not be traced to a numbered reference")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
