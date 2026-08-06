"""How long a terminal inverted repeat arises by chance, and whether MIN_TIR_STANDALONE clears it.

`structural.MIN_TIR_STANDALONE` decides whether a terminal inverted repeat may name a Class II element on
its own, with no transposase to corroborate it. A floor like that is only defensible if it sits above what
coincidence produces, and the only way to know what coincidence produces is to measure it: run the shipped
detector over random sequence and record what it finds. The constant's own comment reports such a
measurement, but the measurement lived in prose, and the shipped test (tests/test_structural_bounds.py)
runs a deliberately cheap 400 trials and says so. Neither reproduces the number that set the constant, so
this script does — with its seed recorded, so the run is repeatable rather than remembered.

The load-bearing result is NOT the single longest hit. That is a maximum over 4,500 draws from a steep
tail, so it moves between seeds; the ceiling is reported, but the distribution and the per-length hit rates
are what a reader should judge the floor against, and the verdict is stated as a comparison against the
constant imported from the module rather than against a number restated here.

Two conditions the output records, because they shape the distribution and would otherwise be read as
properties of chance itself: the lengths reported are conditioned on find_tir's own detection floor
(`min_tir`), its identity gate and its low-complexity guard — a shorter chance repeat exists in the
sequence and is simply not reported — and the sequences are uniform in base composition, whereas real
genomes are usually not, and an AT-rich background yields longer chance repeats than this measures.

    python benchmarks/chance_tir.py                 # the seeded 4,500-sequence measurement
    python benchmarks/chance_tir.py --seeds 5       # and how far the ceiling moves across seeds

The committed output was produced by the second form, so it carries the `seed_sensitivity` block; the
first form writes the same file without it.

Output: benchmarks/raw/chance_tir.json
"""
from __future__ import annotations
import inspect, json, os, random, sys, time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [os.path.join(ROOT, "app", "backend")]

from teagle_core import structural                         # noqa: E402
from teagle_core import __version__ as TEAGLE_VERSION      # noqa: E402

OUT = os.path.join(ROOT, "benchmarks", "raw", "chance_tir.json")

# The three lengths span the range over which the detector is actually used: a short non-autonomous
# element, a typical complete DNA transposon, and a long element or a small contig. Chance repeat length
# grows with template length, so a single length would not bound the quantity the floor has to clear.
LENGTHS = (1500, 4000, 9000)
TRIALS_PER_LENGTH = 1500
SEED = 20260806                       # arbitrary but fixed and recorded


def _defaults(fn, *names):
    """Read a function's own defaults, so the conditions reported are the ones actually applied."""
    sig = inspect.signature(fn)
    return {n: sig.parameters[n].default for n in names}


def run(lengths=LENGTHS, trials_per_length: int = TRIALS_PER_LENGTH, seed: int = SEED) -> dict:
    """Run the measurement and return the result record. Writes nothing."""
    rng = random.Random(seed)
    t0 = time.time()
    per_length = {}
    all_lens = []
    for L in lengths:
        hits = []
        for _ in range(trials_per_length):
            seq = "".join(rng.choices("ACGT", k=L))
            ev = structural.find_tir(seq)
            if ev:
                hits.append((ev["tir_len"], ev["identity"]))
        lens = sorted(n for n, _i in hits)
        all_lens.extend(lens)
        hist = {}
        for n in lens:
            hist[n] = hist.get(n, 0) + 1
        per_length[str(L)] = {
            "trials": trials_per_length,
            "hits": len(hits),
            "hit_rate": round(len(hits) / trials_per_length, 5),
            "longest": max(lens) if lens else None,
            "median_when_found": lens[len(lens) // 2] if lens else None,
            "length_histogram": {str(k): hist[k] for k in sorted(hist)},
            # the quantity the floor is set against: how often chance reaches AT LEAST a given length
            "at_least": {str(k): round(sum(1 for n in lens if n >= k) / trials_per_length, 5)
                         for k in sorted(hist)},
        }

    overall_max = max(all_lens) if all_lens else None
    floor = structural.MIN_TIR_STANDALONE
    total = len(lengths) * trials_per_length
    supports = overall_max is not None and overall_max < floor
    if overall_max is None:                               # no chance hit at all: the floor is unopposed,
        supports = True                                   # but the margin is unmeasured, not infinite
    hist_all = {}
    for n in all_lens:
        hist_all[n] = hist_all.get(n, 0) + 1

    return {
        "trials": total,
        "trials_per_length": trials_per_length,
        "lengths_bp": list(lengths),
        "seed": seed,
        "hits": len(all_lens),
        "hit_rate": round(len(all_lens) / total, 5),
        "longest_chance_tir_bp": overall_max,
        "length_histogram": {str(k): hist_all[k] for k in sorted(hist_all)},
        "at_least": {str(k): round(sum(1 for n in all_lens if n >= k) / total, 5)
                     for k in sorted(hist_all)},
        "per_length": per_length,
        "min_tir_standalone": floor,
        "margin_bp": (floor - overall_max) if overall_max is not None else None,
        "supports_min_tir_standalone": supports,
        "conditions": {
            "find_tir_defaults": _defaults(structural.find_tir, "k", "min_tir", "max_tir", "min_anchors"),
            "identity_floor": structural.MIN_TIR_IDENTITY,
            "note": "reported lengths are conditioned on find_tir's own detection floor, its identity "
                    "gate and its low-complexity guard; sequences are uniform in base composition, and "
                    "an AT-rich background would yield longer chance repeats than this measures",
        },
        "teagle_version": TEAGLE_VERSION,
        "python": sys.version.split()[0],
        "elapsed_s": round(time.time() - t0, 1),
        "generated": time.strftime("%Y-%m-%dT%H:%M:%S"),
    }


def seed_sensitivity(n_seeds: int, lengths=LENGTHS, trials_per_length: int = TRIALS_PER_LENGTH,
                     seed: int = SEED) -> dict:
    """Repeat the whole measurement on consecutive seeds and report how far the ceiling moves.

    The ceiling is a maximum over a steep tail, so one seed fixes it only to within its own sampling
    noise. Anyone quoting the ceiling as a bound should know that width, and a committed procedure for
    measuring it is better than an assertion that it is small. Off by default: the headline run is the
    single seeded measurement above.
    """
    t0 = time.time()
    tops = []
    for s in range(seed, seed + n_seeds):
        r = run(lengths, trials_per_length, s)
        tops.append(r["longest_chance_tir_bp"])
    return {"n_seeds": n_seeds, "first_seed": seed, "ceilings_bp": tops,
            "sequences_total": n_seeds * len(lengths) * trials_per_length,
            "min_bp": min(tops), "max_bp": max(tops), "elapsed_s": round(time.time() - t0, 1)}


def main():
    n_seeds = 1
    if "--seeds" in sys.argv:
        n_seeds = int(sys.argv[sys.argv.index("--seeds") + 1])
    r = run()
    if n_seeds > 1:
        r["seed_sensitivity"] = seed_sensitivity(n_seeds)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(r, open(OUT, "w", encoding="utf-8"), indent=1, sort_keys=True)
    print(f"{r['trials']} random sequences ({r['trials_per_length']} at each of "
          f"{', '.join(str(x) for x in r['lengths_bp'])} bp), seed {r['seed']}")
    print(f"  chance inverted repeat found in {r['hits']} ({r['hit_rate']:.2%})")
    for L in r["lengths_bp"]:
        d = r["per_length"][str(L)]
        print(f"  {L:>5} bp: {d['hits']:>4} hits ({d['hit_rate']:.2%}), longest {d['longest']} bp")
    print(f"  longest chance repeat overall: {r['longest_chance_tir_bp']} bp")
    print(f"  MIN_TIR_STANDALONE = {r['min_tir_standalone']} bp, margin {r['margin_bp']} bp")
    print("  distribution: " + ", ".join(f"{k} bp {v/r['trials']:.2%}"
                                         for k, v in r["length_histogram"].items()))
    if "seed_sensitivity" in r:
        s = r["seed_sensitivity"]
        print(f"  ceiling over {s['n_seeds']} seeds: {s['min_bp']}-{s['max_bp']} bp {s['ceilings_bp']}")
    print(f"-> {OUT}")
    if not r["supports_min_tir_standalone"]:
        print(f"chance reached {r['longest_chance_tir_bp']} bp, at or above the standalone floor of "
              f"{r['min_tir_standalone']} bp: the shipped constant does not clear what this measures")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
