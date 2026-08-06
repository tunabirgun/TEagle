"""Paired comparison of TEagle against TEsorter on byte-identical input.

The pre-registered analysis. TEsorter is reported on BOTH passes: its own default rexdb, which is the
like-for-like comparison since TEagle is also run untuned, and the lineage-matched database its authors
recommend, which is the comparator at its best. Quoting only the pass that flatters TEagle would not be an
evaluation.

WHAT IS PAIRED, AND WHAT IS NOT. Only cases whose corpus label maps onto the order vocabulary
(LTR / LINE / SINE / TIR) enter the paired statistics: negative controls (`NOT_TE`), labels outside the
vocabulary (`OTHER`) and records that `score.py` stratified out as not-single-element cannot be graded at
order level on the same scale as the rest. That filter used to be applied silently, which mattered, because
the negative controls are exactly where the order-level errors live. Every excluded case is now enumerated
in the output under `excluded`, with the reason and with BOTH tools' calls on it, and the order-error count
is emitted twice — over the paired set and over the paired set plus the negative controls — so a statement
about "incorrect order calls" has to name its denominator.

BYTE IDENTITY IS VERIFIED, NOT ASSERTED. Identity holds by construction (`run_tesorter.py` takes the
sequence from TEagle's stored record rather than re-fetching), so a check that merely restates the
construction proves nothing. What is compared here is a hash recomputed from the bytes each tool actually
received: the analysed sequence inside TEagle's raw record, against the residues parsed back out of the
FASTA file `run_tesorter.py` wrote and handed to TEsorter (`raw/tesorter/_work/<accession>.fa`). A pair
whose two hashes differ is excluded and counted. Separately, the recorded fetch-level SHA-256 carried by
both raw records is compared, which catches TEagle being re-run after TEsorter. The two hashes cover
different byte strings — fetched FASTA text with header and wrapping, versus residues only — and are named
apart in the output for that reason. Those input FASTA files are large intermediates and are not committed;
if they have been cleaned up the check reports the affected pairs as unverifiable rather than as verified,
so the count, not a boolean, is what the manuscript should bind.

Discordant pairs are compared with McNemar's test (exact binomial, two-sided), which uses only the cases
where the two tools disagree - the correct denominator for "do these tools differ". Cases where TEagle
abstains and TEsorter answers correctly are counted and reported separately, because that is the direction
in which abstention costs a user information and it is the number the abstention argument lives or dies on.

    python benchmarks/head_to_head.py

Output: benchmarks/raw/head_to_head.json
"""
from __future__ import annotations
import glob, hashlib, json, math, os, sys
from collections import Counter

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [os.path.join(ROOT, "benchmarks")]
from score import expected_order, observed                       # noqa: E402  same mapping both tools

SCORES = os.path.join(ROOT, "benchmarks", "raw", "scores.json")
TES = os.path.join(ROOT, "benchmarks", "raw", "tesorter")
TES_WORK = os.path.join(TES, "_work")
TEAGLE_RAW = os.path.join(ROOT, "benchmarks", "raw", "teagle")
OUT = os.path.join(ROOT, "benchmarks", "raw", "head_to_head.json")

# The four orders both tools are scored on. A call outside this set is not an order call and is never
# counted as an order-level error; it is reported in its own field instead.
ORDER_VOCAB = ("LTR", "LINE", "SINE", "TIR")


def tesorter_order(rec):
    """Map TEsorter's Order column onto the same vocabulary score.py uses for TEagle."""
    cls = rec.get("classifications") or []
    if not cls:
        return "NO_CALL"
    order = str(cls[0].get("order") or "").strip().lower()
    if not order or order in ("unknown", "unclassified", "mixture", "na"):
        return "NO_CALL"
    if order.startswith("ltr"):
        return "LTR"
    if order.startswith("line"):
        return "LINE"
    if order.startswith("sine"):
        return "SINE"
    if order.startswith(("tir", "dna", "mite", "helitron", "maverick")):
        return "TIR"
    return "OTHER"


def mcnemar_exact(b, c):
    """Two-sided exact binomial on the discordant pairs. Exact rather than the chi-square approximation
    because the discordant count here is small, where the approximation is known to be anticonservative."""
    n = b + c
    if n == 0:
        return 1.0
    k = min(b, c)
    tail = sum(math.comb(n, i) for i in range(0, k + 1)) / (2 ** n)
    return min(1.0, 2 * tail)


def fasta_residue_sha256(path):
    """SHA-256 of the residues in a FASTA file, header and line wrapping removed, hashed a line at a time
    so a 150 Mb chromosome does not have to be held as one string. This is the hash of what TEsorter was
    actually given, read back off disk rather than taken from any record that claims what it was given."""
    h = hashlib.sha256()
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                continue
            h.update(line.strip().encode("utf-8"))
    return h.hexdigest()


_TEAGLE_CACHE = {}


def teagle_input_hashes(accession):
    """What TEagle analysed for this accession, as hashes.

    Returns (residue_shas, fetch_shas, n_records). Both are SETS because fourteen corpus rows share an
    accession with another row, so an accession can carry more than one raw record; pairing by accession is
    only legitimate while those records hold the SAME input, and the size of the set is what says so. The
    residue hash is computed over the sequence the engine analysed, which is the same string
    `run_tesorter.py` wrote into the FASTA it passed to TEsorter. The fetch hash is the one recorded at
    fetch time and covers the FASTA text including its header, so it is not comparable to the first."""
    if accession in _TEAGLE_CACHE:
        return _TEAGLE_CACHE[accession]
    safe = accession.replace("/", "_")
    paths = sorted(glob.glob(os.path.join(TEAGLE_RAW, f"*_{safe}.json")))
    legacy = os.path.join(TEAGLE_RAW, safe + ".json")
    if not paths and os.path.exists(legacy):
        paths = [legacy]
    residues, fetch = set(), set()
    for p in paths:
        d = json.load(open(p, encoding="utf-8"))
        recs = (d.get("result") or {}).get("records") or []
        seq = (recs[0] or {}).get("seq") if recs else None
        if seq:
            residues.add(hashlib.sha256(seq.encode("utf-8")).hexdigest())
        if d.get("input_sha256"):
            fetch.add(d["input_sha256"])
    out = (residues, fetch, len(paths))
    _TEAGLE_CACHE[accession] = out
    return out


def check_input_identity(accession, tes_rec):
    """Compare the bytes the two tools received. Returns a dict with a verdict of 'verified',
    'mismatch' or 'unverifiable', the last carrying the reason it could not be checked."""
    residues, fetch, n_recs = teagle_input_hashes(accession)
    fa = os.path.join(TES_WORK, f"{accession.replace('/', '_')}.fa")
    rec = {"accession": accession, "teagle_records": n_recs,
           "teagle_distinct_input_residue_sha256": len(residues)}

    if not residues:
        rec["verdict"] = "unverifiable"
        rec["reason"] = "no TEagle raw record with an analysed sequence for this accession"
    elif not os.path.exists(fa):
        rec["verdict"] = "unverifiable"
        rec["reason"] = "the FASTA handed to TEsorter was not retained (raw/tesorter/_work is an " \
                        "uncommitted intermediate)"
    else:
        got = fasta_residue_sha256(fa)
        rec["tesorter_input_residue_sha256"] = got
        rec["verdict"] = "verified" if got in residues else "mismatch"

    # Independent of the residue check: the fetch-level SHA-256 both raw records carry. Agreement here is
    # weak evidence on its own -- run_tesorter.py copied the value across -- but disagreement is strong,
    # because it means one tool's raw record was regenerated after the other's.
    t_fetch = tes_rec.get("input_sha256")
    if not fetch or not t_fetch:
        rec["recorded_fetch_sha256"] = "unavailable"
    elif t_fetch in fetch:
        rec["recorded_fetch_sha256"] = "agrees"
    else:
        rec["recorded_fetch_sha256"] = "differs"
    return rec


def order_call_error(call, expected):
    """True when a tool committed to an order and the order is wrong; False when it committed and was
    right; None when it made no order call, or made one outside the four-order vocabulary, or the case
    carries no gradable expected order. None means 'not an order call to grade', never 'correct'."""
    if expected not in ORDER_VOCAB and expected != "NOT_TE":
        return None
    if call not in ORDER_VOCAB:
        return None
    if expected == "NOT_TE":
        return True                                   # any order call on a negative control is wrong
    return call != expected


def main():
    if not os.path.exists(SCORES):
        print("run benchmarks/score.py first")
        return 1
    sc = json.load(open(SCORES, encoding="utf-8"))
    by_acc = {}
    for c in sc["cases"]:
        by_acc.setdefault(c["accession"], c)
    # score.py already records why each case it stratified out was stratified out; reuse its wording rather
    # than inventing a second, weaker explanation here.
    stratified = {e["accession"]: e.get("reason") for e in sc.get("excluded_cases", [])}

    tes = {}
    for f in glob.glob(os.path.join(TES, "*.json")):
        if os.path.basename(f).startswith("_"):
            continue
        r = json.load(open(f, encoding="utf-8"))
        tes.setdefault(r["mode"], {})[r["accession"]] = r

    result = {"modes": {}}
    for mode in sorted(tes):
        pairs, excluded, identity = [], [], []
        for acc in sorted(tes[mode]):
            t = tes[mode][acc]
            a = by_acc.get(acc)
            ts = tesorter_order(t)
            idc = check_input_identity(acc, t)
            idc["mode"] = mode
            identity.append(idc)

            # ---- decide whether this case is pairable, and if not, say why -------------------------
            if a is None:
                reason = "not_scored_as_single_element"
            elif a["expected_order"] == "NOT_TE":
                reason = "negative_control"
            elif a["expected_order"] == "OTHER":
                reason = "expected_order_outside_vocabulary"
            elif idc["verdict"] == "mismatch" or idc["recorded_fetch_sha256"] == "differs":
                reason = "input_bytes_differ_between_tools"
            else:
                reason = None

            if reason is not None:
                exp = a["expected_order"] if a else None
                te_call = a["observed_order"] if a else None
                excluded.append({
                    "accession": acc, "reason": reason,
                    "stratification_note": stratified.get(acc) if a is None else None,
                    "expected_order": exp,
                    "teagle": te_call, "teagle_te_class": a["observed_te_class"] if a else None,
                    "teagle_confidence": a["confidence"] if a else None,
                    "teagle_abstained": a["abstained_at_order"] if a else None,
                    "tesorter": ts, "database": t.get("database"),
                    "teagle_order_call_incorrect": order_call_error(te_call, exp) if a else None,
                    "tesorter_order_call_incorrect": order_call_error(ts, exp) if a else None,
                })
                continue

            exp = a["expected_order"]
            te_ok = (not a["abstained_at_order"]) and a["observed_order"] == exp
            te_abst = a["abstained_at_order"]
            ts_ok = ts not in ("NO_CALL", "OTHER") and ts == exp
            pairs.append({"accession": acc, "expected": exp, "teagle": a["observed_order"],
                          "teagle_correct": te_ok, "teagle_abstained": te_abst,
                          "tesorter": ts, "tesorter_correct": ts_ok,
                          "database": t.get("database"),
                          "teagle_order_call_incorrect": order_call_error(a["observed_order"], exp),
                          "tesorter_order_call_incorrect": order_call_error(ts, exp)})

        b = sum(1 for p in pairs if p["teagle_correct"] and not p["tesorter_correct"])
        c = sum(1 for p in pairs if p["tesorter_correct"] and not p["teagle_correct"])
        both = sum(1 for p in pairs if p["teagle_correct"] and p["tesorter_correct"])
        neither = sum(1 for p in pairs if not p["teagle_correct"] and not p["tesorter_correct"])
        cost = [p for p in pairs if p["teagle_abstained"] and p["tesorter_correct"]]
        gain = [p for p in pairs if p["teagle_abstained"] and not p["tesorter_correct"]
                and p["tesorter"] not in ("NO_CALL",)]

        # ---- the exclusion ledger --------------------------------------------------------------
        by_reason = Counter(e["reason"] for e in excluded)
        negatives = [e for e in excluded if e["reason"] == "negative_control"]
        outside = [e for e in excluded if e["reason"] == "expected_order_outside_vocabulary"]
        unscored = [e for e in excluded if e["reason"] == "not_scored_as_single_element"]
        te_neg_err = [e for e in negatives if e["teagle_order_call_incorrect"]]
        ts_neg_err = [e for e in negatives if e["tesorter_order_call_incorrect"]]
        te_pair_err = [p for p in pairs if p["teagle_order_call_incorrect"]]
        ts_pair_err = [p for p in pairs if p["tesorter_order_call_incorrect"]]
        # A call of 'OTHER' is a call, but not one of the four orders, so `order_call_error` grades it as
        # neither right nor wrong. Counting those is what licenses reading an error count of 0 as "no
        # incorrect order call" rather than "no gradable order call": if a tool started emitting OTHER the
        # error count would otherwise stay at 0 while the tool was in fact answering.
        oov = {"teagle": sum(1 for p in pairs if p["teagle"] == "OTHER")
                         + sum(1 for e in negatives if e["teagle"] == "OTHER"),
               "tesorter": sum(1 for p in pairs if p["tesorter"] == "OTHER")
                           + sum(1 for e in negatives if e["tesorter"] == "OTHER")}

        ver = Counter(i["verdict"] for i in identity)
        fetch_agree = Counter(i["recorded_fetch_sha256"] for i in identity)

        result["modes"][mode] = {
            "n_considered": len(tes[mode]),
            "n_paired": len(pairs),
            "n_excluded": len(excluded),
            # ---- byte identity, measured -----------------------------------------------------
            "input_identity": {
                "method": "SHA-256 over the residues of the FASTA run_tesorter.py handed to TEsorter "
                          "(raw/tesorter/_work/<accession>.fa, re-read from disk) compared with SHA-256 "
                          "over the sequence TEagle analysed (raw/teagle/*.json -> result.records[0].seq)",
                "n_checked": len(identity),
                "n_input_residues_verified_identical": ver.get("verified", 0),
                "n_input_residues_mismatched": ver.get("mismatch", 0),
                "n_unverifiable": ver.get("unverifiable", 0),
                "unverifiable_reasons": dict(Counter(
                    i.get("reason") for i in identity if i["verdict"] == "unverifiable")),
                "recorded_fetch_sha256_agrees": fetch_agree.get("agrees", 0),
                "recorded_fetch_sha256_differs": fetch_agree.get("differs", 0),
                "recorded_fetch_sha256_unavailable": fetch_agree.get("unavailable", 0),
                # >0 would mean an accession's raw records hold different inputs, which would make pairing
                # by accession ambiguous. It is reported because the 0 is what licenses the pairing.
                "n_accessions_with_multiple_distinct_inputs": sum(
                    1 for i in identity if i["teagle_distinct_input_residue_sha256"] > 1),
                "per_case": identity,
            },
            "input_sha_mismatches_excluded": by_reason.get("input_bytes_differ_between_tools", 0),
            # ---- what the pairing left out ---------------------------------------------------
            "excluded": {
                "n_excluded": len(excluded),
                "by_reason": dict(by_reason),
                "negative_controls": {
                    "n": len(negatives),
                    "teagle_incorrect_order_calls": [
                        {"accession": e["accession"], "called": e["teagle"],
                         "te_class": e["teagle_te_class"], "confidence": e["teagle_confidence"]}
                        for e in te_neg_err],
                    "tesorter_incorrect_order_calls": [
                        {"accession": e["accession"], "called": e["tesorter"],
                         "database": e["database"]} for e in ts_neg_err],
                    "n_teagle_incorrect": len(te_neg_err),
                    "n_tesorter_incorrect": len(ts_neg_err),
                },
                # Labels outside LTR/LINE/SINE/TIR cannot be scored right or wrong at order level, so a
                # positive call on one is neither. The calls are listed here so that "no incorrect calls
                # among the excluded cases" cannot be read as "no calls among the excluded cases".
                "expected_order_outside_vocabulary": {
                    "n": len(outside),
                    "ungradable_order_calls": [
                        {"accession": e["accession"], "teagle": e["teagle"], "tesorter": e["tesorter"]}
                        for e in outside if e["teagle"] in ORDER_VOCAB or e["tesorter"] in ORDER_VOCAB],
                },
                "not_scored_as_single_element": {
                    "n": len(unscored),
                    "cases": [{"accession": e["accession"], "note": e["stratification_note"],
                               "tesorter": e["tesorter"]} for e in unscored],
                },
                "cases": excluded,
            },
            # ---- order-level errors, with both denominators ----------------------------------
            "order_errors": {
                "paired_denominator": len(pairs),
                "teagle_paired": len(te_pair_err),
                "tesorter_paired": len(ts_pair_err),
                "with_negative_controls_denominator": len(pairs) + len(negatives),
                "teagle_with_negative_controls": len(te_pair_err) + len(te_neg_err),
                "tesorter_with_negative_controls": len(ts_pair_err) + len(ts_neg_err),
                "n_out_of_vocabulary_order_calls": oov,
                "teagle_cases": [p["accession"] for p in te_pair_err] + [e["accession"] for e in te_neg_err],
                "tesorter_cases": ([p["accession"] for p in ts_pair_err]
                                   + [e["accession"] for e in ts_neg_err]),
            },
            "both_correct": both, "teagle_only": b, "tesorter_only": c, "neither": neither,
            "teagle_accuracy": round((both + b) / len(pairs), 4) if pairs else None,
            "tesorter_accuracy": round((both + c) / len(pairs), 4) if pairs else None,
            "mcnemar_discordant": b + c,
            "mcnemar_p_two_sided_exact": float(f"{mcnemar_exact(b, c):.3g}"),
            "abstention_cost_teagle_silent_tesorter_right": len(cost),
            "abstention_vindicated_teagle_silent_tesorter_wrong": len(gain),
            "cost_cases": [p["accession"] for p in cost],
            "tesorter_no_call": sum(1 for p in pairs if p["tesorter"] == "NO_CALL"),
            "pairs": pairs,
        }
        # The ledger has to account for every case that was considered, or it is decoration.
        m = result["modes"][mode]
        assert m["n_considered"] == m["n_paired"] + sum(m["excluded"]["by_reason"].values()), \
            f"{mode}: exclusion ledger does not account for every considered case"

    json.dump(result, open(OUT, "w", encoding="utf-8"), indent=1)

    for mode, m in result["modes"].items():
        e, oe, ii = m["excluded"], m["order_errors"], m["input_identity"]
        print(f"\n=== TEsorter pass: {mode} ===  considered {m['n_considered']} = "
              f"paired {m['n_paired']} + excluded {m['n_excluded']}")
        print(f"  both correct {m['both_correct']:3d}   TEagle only {m['teagle_only']:3d}   "
              f"TEsorter only {m['tesorter_only']:3d}   neither {m['neither']:3d}")
        print(f"  accuracy over paired cases: TEagle {m['teagle_accuracy']}   "
              f"TEsorter {m['tesorter_accuracy']}")
        print(f"  McNemar exact two-sided p = {m['mcnemar_p_two_sided_exact']} "
              f"on {m['mcnemar_discordant']} discordant pairs")
        print(f"  TEagle abstained AND TEsorter was right : {m['abstention_cost_teagle_silent_tesorter_right']}"
              f"  <- what abstention COSTS")
        print(f"  TEagle abstained AND TEsorter was wrong : "
              f"{m['abstention_vindicated_teagle_silent_tesorter_wrong']}  <- what it BUYS")
        print(f"  TEsorter returned no call              : {m['tesorter_no_call']}")
        print(f"  excluded from pairing: " + ", ".join(f"{k}={v}" for k, v in e["by_reason"].items()))
        print(f"  incorrect ORDER calls  over {oe['paired_denominator']} paired: "
              f"TEagle {oe['teagle_paired']}, TEsorter {oe['tesorter_paired']}")
        print(f"  incorrect ORDER calls  over {oe['with_negative_controls_denominator']} paired + negative "
              f"controls: TEagle {oe['teagle_with_negative_controls']} "
              f"{oe['teagle_cases'] or ''}, TEsorter {oe['tesorter_with_negative_controls']} "
              f"{oe['tesorter_cases'] or ''}")
        if e["expected_order_outside_vocabulary"]["ungradable_order_calls"]:
            print(f"  ungradable order calls (label outside LTR/LINE/SINE/TIR): "
                  f"{e['expected_order_outside_vocabulary']['ungradable_order_calls']}")
        print(f"  input bytes: verified identical {ii['n_input_residues_verified_identical']}/"
              f"{ii['n_checked']}, mismatched {ii['n_input_residues_mismatched']}, "
              f"unverifiable {ii['n_unverifiable']}; recorded fetch SHA-256 agrees "
              f"{ii['recorded_fetch_sha256_agrees']}, differs {ii['recorded_fetch_sha256_differs']}")
    print(f"\n-> {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
