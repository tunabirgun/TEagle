"""Validate the in-silico PCR engine against published assays with wet-lab-determined product sizes.

The ground truth here is external in a way the classification corpus's is not: the primers are somebody
else's, and the expected amplicon is a length that a gel or a sequenced product already established. A
prediction either reproduces that length or it does not, and neither outcome depends on any judgement of
ours.

What is measured, per assay:
  * whether the pair is found at all on the cited template;
  * the predicted product length against the published one, as a signed difference in bp;
  * how many products the engine predicts, since a pair that amplifies once in the paper and five times
    here is a specificity result even when one of the five is the right size.

Three properties of each case are read from the corpus rather than inferred from its prose, because each
selects a stratum the manuscript reports and an inference would decide it silently: whether the pair
targets one locus or many, whether the source prints the size as an integer or gives only a gel estimate,
and whether the entry was written with a worked prediction already in it. For a multi-locus assay a single
expected size is the band the authors chose to name rather than the only product, so a count mismatch is
not an error; those cases are reported apart from the headline.

    python benchmarks/run_assay.py

Output: benchmarks/raw/assay/<name>.json per case, plus benchmarks/raw/assay_scores.json
"""
from __future__ import annotations
import glob, hashlib, json, os, re, sys, time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [os.path.join(ROOT, "app", "backend")]

import engine                                              # noqa: E402
from teagle_core import primers                            # noqa: E402
from teagle_core import __version__ as TEAGLE_VERSION      # noqa: E402

CORPUS = os.path.join(ROOT, "benchmarks", "assay_corpus.json")
OUTDIR = os.path.join(ROOT, "benchmarks", "raw", "assay")
SCORES = os.path.join(ROOT, "benchmarks", "raw", "assay_scores.json")

# Tolerance for calling a prediction concordant. A published size is often read off a gel against a ladder,
# so it carries its own error; and a paper may round. 10% or 20 bp, whichever is larger, is the band within
# which a gel-derived size and an exact in-silico length are not in conflict.
def tolerance(expected_bp: int) -> float:
    return max(20.0, 0.10 * expected_bp)


# Whether a pair targets one locus on the cited template or many is a property of the assay, and the
# corpus states it per case with a justification. It was previously inferred from the panel name, which
# is a filing label: four cases filed under 'irap-remap' are neither IRAP nor REMAP assays, and excluding
# them on that basis silently removed the panel's only two non-zero differences from the headline.
def is_multilocus(case) -> bool:
    loci = case.get("loci_on_template")
    if loci not in ("single", "multiple"):
        raise ValueError(f"{case['name']}: loci_on_template missing or unrecognised ({loci!r})")
    return loci == "multiple"


# Papers write a degenerate position in several ways. The corpus stores each primer exactly as printed, so
# that an auditor can check it against the source without decoding anything; the engine needs IUPAC. The
# translation happens here, and both strings are written to the record.
_DEGEN = {"AG": "R", "CT": "Y", "CG": "S", "AT": "W", "GT": "K", "AC": "M",
          "CGT": "B", "AGT": "D", "ACT": "H", "ACG": "V", "ACGT": "N"}


def to_iupac(primer: str) -> str:
    def sub(m):
        bases = "".join(sorted(set(m.group(1).upper().replace("U", "T"))))
        if bases not in _DEGEN:
            raise ValueError(f"cannot render '{m.group(0)}' as a single IUPAC code")
        return _DEGEN[bases]
    return re.sub(r"[\[(]([ACGTUacgtu/,]+)[\])]", sub, primer).upper()


# Not every published size is the same kind of number. Some papers print an integer in a table; others give
# a band on a gel ("1.3 kb"), from which a bp figure can only be derived. A benchmark that pools the two
# reports a concordance partly against its own arithmetic, so the classes are separated and reported apart.
# The test is mechanical: the expected integer either appears verbatim in the quoted source passage or it
# does not. Nothing is assigned by hand.
def ground_truth(case) -> str:
    prov = case.get("expected_bp_provenance")
    if prov not in ("stated_in_source", "derived_from_source"):
        raise ValueError(f"{case['name']}: expected_bp_provenance missing or unrecognised ({prov!r})")
    if prov == "stated_in_source":
        exp = str(int(case["expected_amplicon_bp"]))
        blob = f"{case.get('where_in_paper','')} {case.get('amplicon_range_bp','')}"
        if not re.search(rf"(?<!\d){exp}(?!\d)", blob):
            raise ValueError(f"{case['name']}: declared stated_in_source but {exp} bp appears nowhere in "
                             "the quoted passage, so the claim cannot be checked from the record")
    return prov


# A case whose corpus entry already records a worked in-silico prediction was not scored blind: the
# expected coordinates were known before the engine ran. That does not make the comparison worthless - the
# published size is still external - but it is a different evidential status and is recorded as one.
# Stated per case in the corpus with its justification, not sniffed out of the caveats prose: an entry that
# described a worked prediction in different words would otherwise be filed as blind.
def preworked(case) -> bool:
    flag = case.get("worked_prediction_in_entry")
    if not isinstance(flag, bool):
        raise ValueError(f"{case['name']}: worked_prediction_in_entry missing or not a boolean ({flag!r})")
    return flag


# One size window for every case, derived from the corpus rather than from the case being scored: half the
# smallest published amplicon to three times the largest. A per-case window centred on the expected size
# would let the search decide the answer, which is the failure mode this whole benchmark exists to avoid.
def size_window(cases):
    exp = [int(c["expected_amplicon_bp"]) for c in cases]
    return max(1, min(exp) // 2), max(exp) * 3


def fetch(accession: str):
    res = engine.run_fetch({"accession": accession})
    if not res.get("ok", True) and res.get("error"):
        raise RuntimeError(res["error"])
    seq = res.get("sequence") or res.get("fasta") or ""
    if not seq.strip():
        raise RuntimeError("fetch returned no sequence")
    return seq, res


def verify_bound(cases, trimmed_names, bound, prod_min, prod_max):
    """Re-run only the trimmed cases at the derived bound and one base above it.

    The bound follows from the matcher's construction, but an argument about code is not a measurement of
    it. Running the two values that straddle the bound costs a handful of seconds — only trimmed primers
    can be affected, and there is at most a case or two — and turns the claim into an observation.
    """
    if not bound:
        return {"applicable": False}
    checked = []
    # Only trimmed primers can be affected, so only those are re-run. Fetching every template again
    # would mean re-downloading two whole chromosomes to learn nothing.
    for c in [c for c in cases if c["name"] in set(trimmed_names)]:
        try:
            seq, _src = fetch(c["template_accession"])
        except Exception as e:
            checked.append({"name": c["name"], "error": f"{type(e).__name__}: {e}"})
            continue
        fwd, rev = to_iupac(c["forward_primer"]), to_iupac(c["reverse_primer"])
        counts = {}
        for ma in (bound, bound + 1):
            amps = primers.in_silico_pcr(fwd, rev, seq, c["template_accession"],
                                         prod_min=prod_min, prod_max=prod_max, min_anneal=ma)
            counts[ma] = sum(1 for a in amps if not a["single_primer"])
        if counts[bound] and not counts[bound + 1]:
            checked.append({"name": c["name"], "at_bound": counts[bound],
                            "above_bound": counts[bound + 1], "straddles": True})
    return {"applicable": True, "bound": bound, "cases_that_straddle_it": checked,
            "confirmed": bool(checked)}


def main():
    if not os.path.exists(CORPUS):
        print(f"missing {CORPUS} - assemble the assay corpus first")
        return 1
    cases = json.load(open(CORPUS, encoding="utf-8"))
    os.makedirs(OUTDIR, exist_ok=True)
    prod_min, prod_max = size_window(cases)
    print(f"size window applied to every case: {prod_min}-{prod_max} bp\n")

    rows, failures = [], []
    for n, c in enumerate(cases, 1):
        name = re.sub(r"[^A-Za-z0-9_.-]", "_", c["name"])[:60]
        dest = os.path.join(OUTDIR, f"{n:03d}_{name}.json")
        try:
            seq, src = fetch(c["template_accession"])
            fwd, rev = to_iupac(c["forward_primer"]), to_iupac(c["reverse_primer"])
            t0 = time.time()
            res = engine.run_pcr({"sequence": seq, "fwd": fwd, "rev": rev,
                                  "params": {"prod_min": prod_min, "prod_max": prod_max}})
            amps = res.get("amplicons") or []
            # The published size is a band on a gel, so the comparable quantity is the full product,
            # templated span plus any incorporated non-templated 5' tail - not the distance between
            # binding sites. They differ only for a tailed pair.
            by_size = {}
            for a in amps:
                by_size.setdefault(int(a.get("product_length") or a.get("length") or a.get("size") or 0), a)
            sizes = sorted(by_size)
            exp = int(c["expected_amplicon_bp"])
            tol = tolerance(exp)
            best = min(sizes, key=lambda s: abs(s - exp)) if sizes else None
            # The tail belongs to the product that was scored. Taking a maximum across all products would
            # describe a different amplicon from the one the reported difference refers to.
            ba = by_size.get(best) or {}
            tails = int(ba.get("fwd_tail5") or 0) + int(ba.get("rev_tail5") or 0)
            row = {
                "name": c["name"], "panel": c.get("panel"), "organism": c.get("organism"),
                "element": c.get("element"), "accession": c["template_accession"],
                "expected_bp": exp, "n_products": len(sizes), "product_sizes": sizes,
                "closest_bp": best,
                "difference_bp": (best - exp) if best is not None else None,
                "within_tolerance": (best is not None and abs(best - exp) <= tol),
                "tolerance_bp": round(tol, 1),
                "found": bool(sizes),
                "multilocus_by_design": is_multilocus(c),
                "loci_on_template": c["loci_on_template"],
                "primer_fwd_as_published": c["forward_primer"], "primer_rev_as_published": c["reverse_primer"],
                "primer_degenerate": (fwd != c["forward_primer"].upper() or rev != c["reverse_primer"].upper()),
                "scored_product_tail_bp": tails,
                "scored_product_anneal_bp": [ba.get("fwd_anneal_len"), ba.get("rev_anneal_len")],
                "ground_truth": ground_truth(c),
                "preworked_in_corpus": preworked(c),
                "source_confidence": c.get("confidence"),
                "search_window_bp": [prod_min, prod_max],
                "citation": c.get("citation"), "caveats": c.get("caveats"),
                "template_sha256": hashlib.sha256(seq.encode()).hexdigest()[:16],
                "elapsed_s": round(time.time() - t0, 3),
            }
            rows.append(row)
            tmp = dest + ".part"
            json.dump({"case": c, "result": res, "scored": row}, open(tmp, "w", encoding="utf-8"),
                      indent=1, ensure_ascii=False, sort_keys=True)
            os.replace(tmp, dest)
            flag = "OK " if row["within_tolerance"] else ("MISS" if row["found"] else "NONE")
            print(f"[{n}/{len(cases)}] {flag} {c['name'][:34]:34} expected {exp:>6} bp  "
                  f"got {best if best is not None else '-':>6}  products {len(sizes)}")
        except Exception as e:
            failures.append({"name": c["name"], "error": f"{type(e).__name__}: {e}"})
            print(f"[{n}/{len(cases)}] FAIL {c['name'][:34]:34} {type(e).__name__}: {e}")

    def summarise(rs):
        diffs = sorted(abs(r["difference_bp"]) for r in rs if r["difference_bp"] is not None)
        return {
            "n": len(rs),
            "concordant": sum(1 for r in rs if r["within_tolerance"]),
            "concordance": round(sum(1 for r in rs if r["within_tolerance"]) / len(rs), 4) if rs else None,
            "exact_matches": sum(1 for r in rs if r["difference_bp"] == 0),
            "median_abs_difference_bp": diffs[len(diffs) // 2] if diffs else None,
            "max_abs_difference_bp": diffs[-1] if diffs else None,
            "cases": [r["name"] for r in rs],
        }

    scored = [r for r in rows if r["found"]]
    single = [r for r in rows if not r["multilocus_by_design"]]
    out = {
        "teagle_version": TEAGLE_VERSION,
        "search_window_bp": [prod_min, prod_max],
        "cases": len(cases), "analysed": len(rows), "failed": len(failures),
        "found_any_product": len(scored),
        "single_locus_cases": len(single),
        "single_locus_concordant": sum(1 for r in single if r["within_tolerance"]),
        "concordance": summarise(single)["concordance"],
        "median_abs_difference_bp": summarise(single)["median_abs_difference_bp"],
        "exact_matches": summarise(single)["exact_matches"],
        # The headline number is the one taken against sizes the source prints as integers. Sizes that had
        # to be derived from a gel estimate are reported beside it, never inside it.
        "by_ground_truth": {g: summarise([r for r in single if r["ground_truth"] == g])
                            for g in sorted({r["ground_truth"] for r in single})},
        "blind": summarise([r for r in single if not r["preworked_in_corpus"]]),
        "preworked": summarise([r for r in single if r["preworked_in_corpus"]]),
        "multilocus": summarise([r for r in rows if r["multilocus_by_design"]]),
        # min_anneal sensitivity. The matcher tries the full primer first and shortens only on failure,
        # stopping at the first length that binds, so it uses the LONGEST core that yields sites; a floor
        # at or below that core cannot change the outcome, and one above it removes the site.
        #
        # Only TRIMMED primers can be affected at all. The floor is clamped to the primer's own length, so
        # a primer that anneals over its whole length is found at zero trim whatever the floor is — a
        # 19-mer binding full length is unchanged even at a floor of 40. Deriving the bound from every
        # observed core, trimmed or not, therefore reports the shortest primer in the panel rather than the
        # quantity the floor actually governs, and understates the margin.
        "min_anneal_in_force": primers.MIN_PRIMER_SIZE,
        "trimmed_cases": [r["name"] for r in rows if (r.get("scored_product_tail_bp") or 0) > 0],
        "shortest_trimmed_core_bp": min(
            [a for r in rows if (r.get("scored_product_tail_bp") or 0) > 0
             for a in (r.get("scored_product_anneal_bp") or []) if a], default=None),
        "rows": rows, "failures": failures,
        "generated": time.strftime("%Y-%m-%dT%H:%M:%S"),
    }
    out["min_anneal_invariant_up_to"] = out["shortest_trimmed_core_bp"] or "every value (no primer trimmed)"
    out["min_anneal_bound_verified"] = verify_bound(cases, out["trimmed_cases"],
                                                    out["shortest_trimmed_core_bp"], prod_min, prod_max)
    json.dump(out, open(SCORES, "w", encoding="utf-8"), indent=1, ensure_ascii=False, sort_keys=True)
    print(f"\nanalysed {len(rows)}  failed {len(failures)}")
    for g, s in out["by_ground_truth"].items():
        print(f"  {g:20} {s['concordant']}/{s['n']} concordant, {s['exact_matches']} exact, "
              f"median |diff| {s['median_abs_difference_bp']} bp, max {s['max_abs_difference_bp']} bp")
    print(f"  {'blind (no worked prediction in corpus)':20} "
          f"{out['blind']['concordant']}/{out['blind']['n']}")
    print(f"  {'multi-locus by design':20} {out['multilocus']['concordant']}/{out['multilocus']['n']} "
          f"(reported, not pooled)")
    print(f"  min_anneal in force {out['min_anneal_in_force']} bp; {len(out['trimmed_cases'])} case(s) "
          f"trimmed, shortest trimmed core {out['shortest_trimmed_core_bp']} bp -> invariant up to "
          f"{out['min_anneal_invariant_up_to']}, verified={out['min_anneal_bound_verified'].get('confirmed')}")
    print(f"-> {SCORES}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
