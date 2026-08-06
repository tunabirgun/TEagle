"""Run TEsorter over the same benchmark corpus TEagle saw, and write RAW output per case.

Fairness is the point of this file, so the choices are stated rather than buried:

* Both tools receive the SAME fetched FASTA — the file `run_teagle.py` already wrote and hashed. TEsorter
  is never given a cleaner or differently-trimmed input.
* Two passes are run and reported separately:
    - `default`  : `-db rexdb`, TEsorter's own default, no per-case tuning. This is the head-to-head,
                   because TEagle is also run untuned.
    - `matched`  : the lineage-appropriate database (`rexdb-plant` / `rexdb-metazoa`) chosen from the
                   corpus row's organism. TEsorter's authors recommend this, so reporting only the
                   default pass would understate it. The matched pass is the comparator AT ITS BEST.
  A paper that quotes only the pass flattering to TEagle is not an evaluation.
* TEsorter runs in its own micromamba environment (`tesorter`), never the shipped `te` environment that
  serves the RepeatMasker/Dfam backend.

    python benchmarks/run_tesorter.py            # both passes, resumable
    python benchmarks/run_tesorter.py --pass default
"""
from __future__ import annotations
import argparse, csv, json, os, re, subprocess, sys, time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CORPUS = os.path.join(ROOT, "benchmarks", "corpus.tsv")
TEAGLE_RAW = os.path.join(ROOT, "benchmarks", "raw", "teagle")
OUTDIR = os.path.join(ROOT, "benchmarks", "raw", "tesorter")
MAMBA = "~/.local/bin/micromamba"

# Lineage databases TEsorter ships. The mapping is by organism kingdom only — never by expected answer,
# which would leak the ground truth into the comparator's configuration.
PLANT = re.compile(r"arabidopsis|oryza|zea |maize|rice|nicotiana|tobacco|solanum|triticum|wheat|"
                   r"hordeum|barley|glycine|soybean|vitis|populus|brassica|medicago|sorghum", re.I)
METAZOA = re.compile(r"drosophila|homo |human|mus |mouse|danio|zebrafish|caenorhabditis|anopheles|"
                     r"aedes|bombyx|gallus|xenopus|rattus|macaca|pan troglodytes", re.I)


def wsl(cmd: str, timeout=1800):
    """Run one command inside WSL. Returns (rc, stdout, stderr) — never raises on a non-zero rc, because a
    tool failing on a case is a RESULT to record, not an exception to swallow."""
    p = subprocess.run(["wsl", "-e", "bash", "-lc", cmd], capture_output=True, text=True,
                       timeout=timeout, encoding="utf-8", errors="replace")
    return p.returncode, p.stdout, p.stderr


def db_for(organism: str, mode: str) -> str:
    if mode == "default":
        return "rexdb"
    if PLANT.search(organism or ""):
        return "rexdb-plant"
    if METAZOA.search(organism or ""):
        return "rexdb-metazoa"
    return "rexdb"                                    # unknown lineage -> the general database


def to_wsl_path(win_path: str) -> str:
    p = os.path.abspath(win_path).replace("\\", "/")
    return "/mnt/" + p[0].lower() + p[2:]


def sequence_for(accession: str, row_index: int | None = None):
    """The EXACT sequence TEagle analysed, taken from its raw output — not re-fetched. Re-fetching could
    return a revised record and would silently make the two tools' inputs differ.

    Raw files are keyed by CORPUS ROW, not by accession, because several corpus cases share an accession
    with another case (different elements inside one deposited record). Looking up by accession alone
    would hand TEsorter whichever of those happened to be written last, so the row index is used when it
    is known and an accession match is only a fallback."""
    p = None
    if row_index is not None:
        cand = os.path.join(TEAGLE_RAW, f"{row_index:03d}_{accession.replace('/', '_')}.json")
        if os.path.exists(cand):
            p = cand
    if p is None:
        import glob as _glob
        hits = sorted(_glob.glob(os.path.join(TEAGLE_RAW, f"*_{accession.replace('/', '_')}.json")))
        legacy = os.path.join(TEAGLE_RAW, accession.replace("/", "_") + ".json")
        if hits:
            p = hits[0]
        elif os.path.exists(legacy):
            p = legacy
    if p is None or not os.path.exists(p):
        return None, None
    rec = json.load(open(p, encoding="utf-8"))
    recs = (rec.get("result") or {}).get("records") or []
    if not recs:
        return None, rec.get("input_sha256")
    return recs[0].get("seq"), rec.get("input_sha256")


def parse_cls(text: str):
    """TEsorter's .cls.tsv: one row per sequence, columns #TE Order Superfamily Clade Complete Strand Domains."""
    rows = []
    for line in (text or "").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        rows.append({"id": f[0] if len(f) > 0 else "", "order": f[1] if len(f) > 1 else "",
                     "superfamily": f[2] if len(f) > 2 else "", "clade": f[3] if len(f) > 3 else "",
                     "complete": f[4] if len(f) > 4 else "", "strand": f[5] if len(f) > 5 else "",
                     "domains": f[6] if len(f) > 6 else ""})
    return rows


def run_one(acc, organism, seq, mode, workdir):
    db = db_for(organism, mode)
    os.makedirs(workdir, exist_ok=True)
    fa = os.path.join(workdir, f"{acc}.fa")
    with open(fa, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(f">{acc}\n")
        for i in range(0, len(seq), 60):
            fh.write(seq[i:i + 60] + "\n")
    wfa = to_wsl_path(fa)
    wdir = to_wsl_path(workdir)
    t0 = time.time()
    rc, out, err = wsl(f"cd {wdir} && {MAMBA} run -n tesorter TEsorter {wfa} -db {db} -p 4 -pre {acc}.{mode}")
    elapsed = round(time.time() - t0, 3)
    cls_path = os.path.join(workdir, f"{acc}.{mode}.cls.tsv")
    cls_text = open(cls_path, encoding="utf-8").read() if os.path.exists(cls_path) else ""
    return {"accession": acc, "organism": organism, "mode": mode, "database": db,
            "returncode": rc, "elapsed_s": elapsed,
            "classifications": parse_cls(cls_text),
            "cls_tsv_raw": cls_text,
            "stdout_tail": (out or "")[-1500:], "stderr_tail": (err or "")[-1500:]}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--corpus", default=CORPUS)
    ap.add_argument("--outdir", default=OUTDIR)
    ap.add_argument("--pass", dest="passes", default="both", choices=["default", "matched", "both"])
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()

    rows = list(csv.DictReader(open(a.corpus, encoding="utf-8", newline=""), delimiter="\t"))
    modes = ["default", "matched"] if a.passes == "both" else [a.passes]
    workdir = os.path.join(a.outdir, "_work")
    os.makedirs(a.outdir, exist_ok=True)

    rc, ver, _ = wsl(f"{MAMBA} run -n tesorter TEsorter --version")
    version = (ver or "").strip().splitlines()[-1] if ver else "unknown"
    print(f"TEsorter: {version}")

    done = failed = skipped = noseq = 0
    for n, row in enumerate(rows, 1):
        if a.limit and done + failed >= a.limit:
            break
        acc = row["accession"].strip()
        organism = (row.get("organism") or "").strip()
        seq, sha = sequence_for(acc, n)
        if not seq:
            noseq += 1
            print(f"[{n}/{len(rows)}] {acc:14} SKIP - no TEagle raw output (run run_teagle.py first)")
            continue
        for mode in modes:
            dest = os.path.join(a.outdir, f"{acc}.{mode}.json")
            if os.path.exists(dest) and not a.force:
                skipped += 1
                continue
            try:
                rec = run_one(acc, organism, seq, mode, workdir)
                rec["input_sha256"] = sha              # proves both tools saw the same bytes
                tmp = dest + ".part"
                json.dump(rec, open(tmp, "w", encoding="utf-8"), indent=1, ensure_ascii=False, sort_keys=True)
                os.replace(tmp, dest)
                c = rec["classifications"]
                call = f"{c[0]['order']}/{c[0]['superfamily']}" if c else "(no call)"
                print(f"[{n}/{len(rows)}] {acc:14} {mode:8} {rec['database']:16} {call:28} {rec['elapsed_s']}s")
                done += 1
            except Exception as e:
                failed += 1
                print(f"[{n}/{len(rows)}] {acc:14} {mode:8} FAILED {type(e).__name__}: {e}")

    json.dump({"tesorter_version": version, "modes": modes,
               "corpus_rows": len(rows), "counts": {"run": done, "skipped_existing": skipped,
                                                    "failed": failed, "no_teagle_input": noseq},
               "finished": time.strftime("%Y-%m-%dT%H:%M:%S")},
              open(os.path.join(a.outdir, "_run.json"), "w", encoding="utf-8"),
              indent=1, ensure_ascii=False, sort_keys=True)
    print(f"\nrun {done} · skipped {skipped} · failed {failed} · no-input {noseq}")


if __name__ == "__main__":
    main()
