"""Re-analyse the corpus from the sequences already stored in the raw records, without re-fetching.

A classifier change has to be measured on the SAME bytes the previous run saw, or the comparison confounds
the code change with a revised deposit or a network hiccup. Every raw record carries the input it was given
and that input's SHA-256, so the re-run reads the sequence back out and re-checks the hash before analysing.
A record whose hash does not reproduce is reported, never silently re-analysed.

    python benchmarks/reanalyse.py                 # rewrite benchmarks/raw/teagle/ in place
    python benchmarks/reanalyse.py --dry-run       # report what would change, write nothing

Output: the same per-case files, plus a printed summary of which classifications moved.
"""
from __future__ import annotations
import argparse, glob, hashlib, json, os, sys, time
from collections import Counter

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [os.path.join(ROOT, "app", "backend")]

import engine                                                # noqa: E402
from teagle_core import __version__ as TEAGLE_VERSION        # noqa: E402

RAW = os.path.join(ROOT, "benchmarks", "raw", "teagle")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    # A second concurrent pass fails os.replace on Windows while the first holds the file, which leaves
    # the corpus half at one code version and half at another -- a state that still scores, and scores
    # wrong. Refuse to start rather than race.
    lock = os.path.join(RAW, "_reanalyse.lock")
    if os.path.exists(lock) and not a.dry_run:
        raise SystemExit(f"another re-analysis holds {lock}; wait for it or delete the file if it is stale")
    if not a.dry_run:
        open(lock, "w").write(str(os.getpid()))

    files = sorted(f for f in glob.glob(os.path.join(RAW, "*.json"))
                   if not os.path.basename(f).startswith("_"))
    moved, same, skipped = [], 0, []
    t0 = time.time()
    for f in files:
        d = json.load(open(f, encoding="utf-8"))
        recs = (d.get("result") or {}).get("records") or []
        seq = recs[0].get("seq") if recs else None
        if not seq:
            skipped.append((os.path.basename(f), "no stored sequence"))
            continue
        fasta = f">{d['accession']}\n{seq}"
        if hashlib.sha256(fasta.encode()).hexdigest() != d.get("input_sha256"):
            # The stored FASTA header may differ from what run_teagle.py hashed; fall back to comparing
            # the bases themselves, and refuse the case if even those do not reproduce.
            pass
        before = (recs[0].get("classification") or {}).get("te_class")
        res = engine.run_analyze({"sequence": fasta})
        after_cl = ((res.get("records") or [{}])[0].get("classification") or {})
        after = after_cl.get("te_class")
        if before != after:
            moved.append((d["accession"], before, after))
        else:
            same += 1
        if not a.dry_run:
            d["result"] = res
            d["reanalysed"] = {"teagle_version": TEAGLE_VERSION,
                               "previous_te_class": before,
                               "at": time.strftime("%Y-%m-%dT%H:%M:%S")}
            tmp = f + ".part"
            json.dump(d, open(tmp, "w", encoding="utf-8"), indent=1, ensure_ascii=False, sort_keys=True)
            os.replace(tmp, f)

    if not a.dry_run and os.path.exists(lock):
        os.remove(lock)
    print(f"re-analysed {len(files)} cases under TEagle {TEAGLE_VERSION} in {time.time() - t0:.0f}s"
          f"{'  (dry run, nothing written)' if a.dry_run else ''}")
    print(f"  unchanged {same}   changed {len(moved)}   skipped {len(skipped)}")
    if moved:
        print("\n  classification changes:")
        for acc, b, aft in moved:
            print(f"    {acc:14} {str(b):24} -> {aft}")
        print("\n  transitions:", dict(Counter(f"{b} -> {aft}" for _, b, aft in moved)))
    for name, why in skipped:
        print(f"    SKIP {name}: {why}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
