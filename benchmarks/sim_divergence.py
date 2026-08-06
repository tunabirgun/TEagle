"""Measure TEagle's terminal-repeat detection boundary as a function of LTR-LTR divergence.

Why simulate. The two LTRs of a retrotransposon are identical at insertion and diverge with time, so
LTR-LTR identity is an age proxy. To measure where a detector stops seeing the pair you need cases whose
true identity is known exactly, across a controlled gradient. Real deposits do not come with an exact,
independently measured identity at every step. Earl Grey (Baril et al. 2024, Mol Biol Evol 41:msae068)
established the same approach for the same reason: simulated sequence with known ground truth for the
quantitative panels, real genomes for validation.

What is simulated and what is not. Only the DIVERGENCE is simulated. The element scaffold is a real
deposited retrotransposon; one LTR copy is mutated at a controlled per-base rate while the internal
coding region and the other copy are left untouched. Ground truth is therefore exact: the pair IS a
terminal direct repeat at a known identity, and any non-detection is a false negative attributable to the
detector, not to an ambiguous input.

Substitutions only, no indels: an indel would change the alignment problem as well as the identity, and
this panel measures one variable.

    python benchmarks/sim_divergence.py                     # full series
    python benchmarks/sim_divergence.py --replicates 3      # fewer replicates, quicker

Output: benchmarks/raw/sim_divergence.json - one record per (identity target, replicate).
"""
from __future__ import annotations
import argparse, hashlib, json, os, random, sys, time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [os.path.join(ROOT, "app", "backend")]

import engine                                              # noqa: E402
from teagle_core import __version__ as TEAGLE_VERSION      # noqa: E402
from teagle_core import structural                         # noqa: E402

OUT = os.path.join(ROOT, "benchmarks", "raw", "sim_divergence.json")
BASES = "ACGT"

# Identity targets, dense across the two boundaries this panel exists to locate:
#   the acceptance floor (structural.MIN_LTR_IDENTITY, 80%) and the lower limit of k-mer seeding.
# The seeding limit is deliberately not given a number here. An earlier comment put it at ~72% from a
# single rough test; this panel is what retracted that, finding candidates still reported at 65% and
# the limit length-dependent rather than constant. Sampling stays dense through that region.
TARGETS = [100, 99, 98, 96, 94, 92, 90, 88, 86, 85, 84, 83, 82, 81, 80.5, 80, 79.5, 79, 78,
           77, 76, 75, 74, 73, 72, 71, 70, 68, 65, 60]


def mutate(seq: str, pct_identity: float, rng: random.Random) -> str:
    """Substitute a fraction of positions so the copy retains `pct_identity` percent of the original.
    Every substitution is to a DIFFERENT base, so the realised identity is exactly the target (a random
    replacement would restore the original base a quarter of the time and inflate identity)."""
    s = list(seq)
    n_sub = int(round(len(s) * (100.0 - pct_identity) / 100.0))
    for i in rng.sample(range(len(s)), min(n_sub, len(s))):
        s[i] = rng.choice([b for b in BASES if b != s[i].upper()])
    return "".join(s)


def build_element(ltr: str, internal: str, ltr3: str) -> str:
    return ltr + internal + ltr3


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--accession", default="D83003",
                    help="scaffold element; default Tto1 (tobacco copia, Hirochika 1993) - NOT a bundled example")
    ap.add_argument("--replicates", type=int, default=5)
    ap.add_argument("--seed", type=int, default=20260804)
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args()

    print(f"fetching scaffold {a.accession} ...")
    fetched = engine.run_fetch({"accession": a.accession})
    seq = (fetched.get("sequence") or fetched.get("fasta") or "")
    seq = "".join(l.strip() for l in seq.splitlines() if not l.startswith(">")).upper()
    if len(seq) < 1000:
        raise SystemExit(f"fetch returned {len(seq)} bp - cannot build a series from that")

    # Locate the element's own LTRs with the detector, then rebuild the element from those coordinates so
    # the scaffold is a REAL element with REAL termini, not a construct.
    base = structural.find_ltr(seq)
    if not base:
        raise SystemExit(f"{a.accession}: no LTR pair found in the scaffold; choose another accession")
    l5s, l5e = base["five_prime"]
    l3s, l3e = base["three_prime"]
    ltr5, internal, ltr3 = seq[l5s:l5e], seq[l5e:l3s], seq[l3s:l3e]
    print(f"scaffold: LTR {len(ltr5)} bp, internal {len(internal)} bp, "
          f"native identity {base['identity']}%")

    records = []
    t0 = time.time()
    for target in TARGETS:
        for rep in range(a.replicates):
            rng = random.Random(a.seed + int(target * 100) * 1000 + rep)
            mutated3 = mutate(ltr3, target, rng)
            elem = build_element(ltr5, internal, mutated3)
            realised = structural._identity(ltr5, mutated3)
            res = engine.run_analyze({"sequence": f">sim_{target}_{rep}\n{elem}"})
            rec = (res.get("records") or [{}])[0]
            ev = rec.get("structural") or []
            accepted = next((e for e in ev if e["type"].startswith("LTR")), None)
            advisory = next((e for e in ev if e.get("advisory")), None)
            cl = rec.get("classification") or {}
            records.append({
                "target_identity": target,
                "realised_identity": realised,
                "replicate": rep,
                "detected_as_ltr": accepted is not None,
                "reported_identity": (accepted or {}).get("identity"),
                "reported_ltr_len": (accepted or {}).get("ltr_len"),
                "advisory_reported": advisory is not None,
                "advisory_identity": (advisory or {}).get("identity"),
                "te_class": cl.get("te_class"),
                "superfamily": cl.get("superfamily"),
                "confidence": cl.get("confidence"),
                "element_sha256": hashlib.sha256(elem.encode()).hexdigest()[:16],
            })
        det = sum(1 for r in records[-a.replicates:] if r["detected_as_ltr"])
        adv = sum(1 for r in records[-a.replicates:] if r["advisory_reported"])
        print(f"  identity {target:5}%  detected {det}/{a.replicates}  advisory {adv}/{a.replicates}")

    meta = {
        "teagle_version": TEAGLE_VERSION,
        "scaffold_accession": a.accession,
        "scaffold_ltr_len": len(ltr5),
        "scaffold_internal_len": len(internal),
        "scaffold_native_identity": base["identity"],
        "min_ltr_identity_floor": structural.MIN_LTR_IDENTITY,
        "seed": a.seed,
        "replicates": a.replicates,
        "targets": TARGETS,
        "elapsed_s": round(time.time() - t0, 1),
        "generated": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "records": records,
    }
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    tmp = a.out + ".part"
    json.dump(meta, open(tmp, "w", encoding="utf-8"), indent=1, sort_keys=True)
    os.replace(tmp, a.out)
    print(f"\n{len(records)} runs -> {a.out}")


if __name__ == "__main__":
    main()
