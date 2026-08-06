"""Run TEagle over the benchmark corpus and write RAW output, one JSON per case.

Nothing here scores anything. This step exists so that every number in the paper is traceable to a file
produced by executing the shipped tool on a fetched accession — never to a figure someone remembered,
transcribed, or asked a model for.

    python benchmarks/run_teagle.py                 # fetch + analyse every case, resumable
    python benchmarks/run_teagle.py --limit 5       # smoke test
    python benchmarks/run_teagle.py --only AY129008 # one case

Resumable: a case whose raw file already exists is skipped, so a network failure costs only the
unfinished cases. Delete the raw file to force a re-run.

Outputs
  benchmarks/raw/teagle/<accession>.json   the full engine result plus the fetched input's sha256
  benchmarks/raw/teagle/_run.json          tool versions, panel hash, corpus hash, start/end time
"""
from __future__ import annotations
import argparse, csv, hashlib, json, os, re, sys, time, traceback

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [os.path.join(ROOT, "app", "backend")]

import engine                                        # noqa: E402
from teagle_core import __version__ as TEAGLE_VERSION  # noqa: E402
from teagle_core import domains, provenance          # noqa: E402

CORPUS = os.path.join(ROOT, "benchmarks", "corpus.tsv")
OUTDIR = os.path.join(ROOT, "benchmarks", "raw", "teagle")


def _sha(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def load_corpus(path):
    """One dict per row. The corpus is the single source of truth for what is benchmarked; this reader
    deliberately does no filtering, so a case cannot be silently dropped between corpus and results."""
    with open(path, encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for i, r in enumerate(rows, 1):
        if not r.get("accession"):
            raise SystemExit(f"corpus row {i} has no accession")
    return rows


def fetch(accession, coords=""):
    """Fetch the case's sequence through the SAME engine path the application uses, so the benchmark
    measures the shipped fetch-and-analyse chain rather than a bespoke harness."""
    # The corpus records coordinates as they appear in the source publication, which is not one format:
    # GenBank join() expressions, UCSC browser ranges, cytogenetic bands and plain spans all occur. Only a
    # plain span on the fetched accession can be applied without an assembly context, so that is the only
    # case narrowed here; everything else is analysed as the whole deposited record and SAYS SO, because a
    # silently unapplied coordinate would make the case look narrower than it was.
    span = re.match(r"^\s*(\d[\d,]*)\s*(?:\.\.|-|–)\s*(\d[\d,]*)", coords or "")
    applied = None
    res = engine.run_fetch({"accession": accession})
    if not res.get("ok", True) and res.get("error"):
        raise RuntimeError(res["error"])
    seq = res.get("sequence") or res.get("fasta") or ""
    if not seq.strip():
        raise RuntimeError("fetch returned no sequence")

    if span:
        head, _, body = seq.partition("\n") if seq.startswith(">") else ("", "", seq)
        bases = "".join(l.strip() for l in body.splitlines() if not l.startswith(">"))
        a, b = (int(span.group(1).replace(",", "")), int(span.group(2).replace(",", "")))
        if 0 < a <= b <= len(bases):
            seq = (head + "\n" if head else "") + bases[a - 1:b]
            applied = f"{a}-{b}"
        else:
            applied = f"OUT OF RANGE {a}-{b} on {len(bases)} bp - whole record used"
    elif coords:
        applied = f"NOT APPLICABLE without an assembly ({coords[:60]}) - whole record used"
    res = dict(res)
    res["coords_applied"] = applied
    return seq, res


def run_case(row):
    acc = row["accession"].strip()
    seq, fetch_res = fetch(acc, (row.get("coords") or "").strip())
    t0 = time.time()
    result = engine.run_analyze({"sequence": seq})
    elapsed = round(time.time() - t0, 3)
    return {
        "accession": acc,
        "corpus_row": row,
        "input_sha256": _sha(seq),
        "input_length": len("".join(l for l in seq.splitlines() if not l.startswith(">"))),
        "fetch_source": {k: fetch_res.get(k) for k in
                         ("accession", "organism", "title", "length", "coords_applied") if k in fetch_res},
        "elapsed_s": elapsed,
        "result": result,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--corpus", default=CORPUS)
    ap.add_argument("--outdir", default=OUTDIR)
    ap.add_argument("--limit", type=int, default=0, help="stop after N cases (smoke test)")
    ap.add_argument("--only", default="", help="run a single accession")
    ap.add_argument("--force", action="store_true", help="re-run cases that already have raw output")
    a = ap.parse_args()

    rows = load_corpus(a.corpus)
    if a.only:
        rows = [r for r in rows if r["accession"].strip() == a.only]
        if not rows:
            raise SystemExit(f"{a.only} is not in the corpus")
    os.makedirs(a.outdir, exist_ok=True)

    corpus_hash = _sha(open(a.corpus, encoding="utf-8").read())
    started = time.strftime("%Y-%m-%dT%H:%M:%S")
    ok = skipped = failed = 0
    failures = []

    for n, row in enumerate(rows, 1):
        if a.limit and ok + failed >= a.limit:
            break
        acc = row["accession"].strip()
        # Keyed by row, not by accession. Several corpus cases share an accession with another case -
        # different elements within one deposited record - and keying on the accession alone silently
        # overwrote them, so the run appeared complete while losing a tenth of the corpus. The count is
        # deliberately not written here; it is a property of the corpus and changes when the corpus does.
        dest = os.path.join(a.outdir, f"{n:03d}_{acc.replace('/', '_')}.json")
        if os.path.exists(dest) and not a.force:
            skipped += 1
            continue
        try:
            rec = run_case(row)
            tmp = dest + ".part"                      # same atomic discipline the app uses for exports
            with open(tmp, "w", encoding="utf-8") as fh:
                json.dump(rec, fh, indent=1, ensure_ascii=False, sort_keys=True)
            os.replace(tmp, dest)
            cl = (rec["result"]["records"][0].get("classification") or {}) if rec["result"].get("records") else {}
            print(f"[{n}/{len(rows)}] {acc:14} {cl.get('te_class') or '-':22} "
                  f"{cl.get('confidence') or '-':10} {rec['elapsed_s']}s")
            ok += 1
        except Exception as e:                        # a failed case is DATA, not a reason to stop
            failed += 1
            failures.append({"accession": acc, "error": f"{type(e).__name__}: {e}",
                             "traceback": traceback.format_exc()[-800:]})
            print(f"[{n}/{len(rows)}] {acc:14} FAILED  {type(e).__name__}: {e}")

    meta = {
        "teagle_version": TEAGLE_VERSION,
        "panel_profiles": len(domains.DOMAIN_INFO),
        "panel_sha256": getattr(domains, "HMM_SHA256", None),
        "detector_parameters": engine._detector_parameters(),
        "provenance_schema": getattr(provenance, "SCHEMA_VERSION", None),
        "corpus_path": os.path.relpath(a.corpus, ROOT).replace("\\", "/"),
        "corpus_sha256": corpus_hash,
        "corpus_rows": len(rows),
        "started": started,
        "finished": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "counts": {"analysed": ok, "skipped_existing": skipped, "failed": failed},
        "failures": failures,
        "python": sys.version.split()[0],
    }
    with open(os.path.join(a.outdir, "_run.json"), "w", encoding="utf-8") as fh:
        json.dump(meta, fh, indent=1, ensure_ascii=False, sort_keys=True)

    print(f"\nanalysed {ok} · skipped {skipped} · failed {failed}")
    print(f"raw output -> {a.outdir}")
    if failures:
        print("FAILED CASES (these are reported in the paper, not dropped):")
        for f in failures:
            print("  ", f["accession"], f["error"][:90])


if __name__ == "__main__":
    main()
