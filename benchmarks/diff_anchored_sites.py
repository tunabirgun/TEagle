"""Differential test: the single-pass 3'-anchored matcher against the descending-scan semantics it replaced.

`primers.anchored_sites` finds the longest 3'-proximal core of a primer that anneals. The literal way to
do that is to try the whole primer, then progressively longer 5' tails, and stop at the first core length
that yields sites; it costs one full template scan per tail length, which is unusable on a chromosome. The
shipped implementation scans once, at the shortest permitted core, and extends each hit 5' while the
mismatch budget and the strict 3' rule hold. That rewrite is admissible only if it returns the same sites,
and "only if" is a claim about behaviour, not about intent -- so it is measured here rather than argued.

The check exists in the repository because the manuscript reports its result. An equivalence measured in a
scratch file is a recollection; a reader cannot re-run it, and cannot tell a run that agreed from a run
that was never made. The first attempt at the replacement did disagree, in the regime where the strict 3'
window is wider than the shortest permitted core (`tp > min_anneal`), which is exactly the regime a
uniform sample of the stated parameter ranges almost never reaches: with the ranges below it accounts for
about 2% of draws. A quarter of the trials are therefore drawn from that regime by construction, and the
two stratum counts are written to the output, so the headline number cannot be read as covering an axis it
did not exercise.

What the reference reimplements is the SCAN STRUCTURE, which is what changed. It reuses the module's
base-level IUPAC predicate and reverse-complement helper: retyping the ambiguity table would add a
typo-shaped failure mode that manufactures disagreements the code does not have, and the matcher itself is
covered by tests/test_primers.py.

Compared per site: strand-appropriate position (`left` or `right`), mismatch count, mismatch positions,
5' tail length and core length. Any disagreement is a defect in one of the two implementations, not a
tolerance, so the script exits non-zero.

    python benchmarks/diff_anchored_sites.py

Output: benchmarks/raw/diff_anchored_sites.json
"""
from __future__ import annotations
import json, os, random, sys, time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [os.path.join(ROOT, "app", "backend")]

from teagle_core import primers                            # noqa: E402
from teagle_core.sequtil import reverse_complement         # noqa: E402
from teagle_core import __version__ as TEAGLE_VERSION      # noqa: E402

OUT = os.path.join(ROOT, "benchmarks", "raw", "diff_anchored_sites.json")

TRIALS = 6000
SEED = 20260806                       # arbitrary but fixed and recorded: the run is reproducible verbatim

# Parameter ranges, inclusive. `min_anneal` brackets the shipped default (primers.MIN_PRIMER_SIZE) on both
# sides so the default is not a boundary of the tested space; the primer-length range spans a short
# diagnostic oligo to a long tailed consensus primer.
MAX_MM_RANGE = (0, 3)
TP_RANGE = (0, 8)
MIN_ANNEAL_RANGE = (6, 22)
PRIMER_LEN_RANGE = (10, 42)
SEQ_LEN_RANGE = (90, 500)

# Share of trials drawn so that the strict 3' window is wider than the shortest permitted core. This is the
# only regime in which the two implementations' clamps could differ -- the shipped code clamps once, at the
# shortest core, and the reference clamps per core length -- and it is where the first replacement failed.
STRICT_SHARE = 0.25

# A site scan that hits primers.MAX_SITES stops early; the reference has no such cap, so a capped scan is a
# legitimate divergence rather than a defect. Templates are kept short enough that the cap cannot bite:
# the candidate count is bounded by the number of template positions.
assert SEQ_LEN_RANGE[1] < primers.MAX_SITES, "template range can reach the site cap; trials would not be comparable"
assert MIN_ANNEAL_RANGE[0] < primers.MIN_PRIMER_SIZE < MIN_ANNEAL_RANGE[1], \
    "min_anneal range no longer brackets the shipped minimum primer size"

# Degenerate codes are read off the module's own ambiguity table rather than restated, so a code added
# there is exercised here without an edit. A code is degenerate when it covers more than one base, which
# excludes U (an RNA spelling of T, not an ambiguity).
DEGENERATE_CODES = tuple(sorted(c for c, allowed in primers._IUPAC.items() if len(allowed) > 1))

# (min_anneal, tp) pairs from the declared ranges for which the strict window exceeds the shortest core.
# Derived from the ranges: change a range and the stratum follows, or fails loudly if it becomes empty.
_STRICT_PAIRS = [(ma, tp)
                 for ma in range(MIN_ANNEAL_RANGE[0], MIN_ANNEAL_RANGE[1] + 1)
                 for tp in range(TP_RANGE[0], TP_RANGE[1] + 1)
                 if tp > ma and ma <= PRIMER_LEN_RANGE[0]]     # ma <= shortest primer, so floor == ma
assert _STRICT_PAIRS, "no (min_anneal, tp) pair in the declared ranges reaches tp > min_anneal"


# --------------------------------------------------------------------------------------------------
# Reference implementation: the descending scan the single-pass matcher replaced.
# --------------------------------------------------------------------------------------------------
def descending_sites(primer: str, seq: str, max_mm: int, tp: int, min_anneal: int, on_reverse: bool,
                     clamp_at_floor: bool = False):
    """Binding sites for the longest 3'-proximal core that anneals, found by rescanning per core length.

    Tries the full primer first and shortens one base at a time; the first core length that yields any
    surviving site wins, and the sites at that length are the answer. The strict 3' rule is applied with
    t = min(tp, n) at EACH core length n, because the window cannot extend past the core it constrains.
    Getting that clamp wrong (fixing t once, at the shortest core) is the historical defect this harness
    was built to catch, so it is stated explicitly here rather than inherited.

    `clamp_at_floor=True` reintroduces that defect deliberately. It is not used for the equivalence
    measurement; it is the sensitivity control, and a trial set that cannot tell the two apart has not
    established anything about the clamp.
    """
    L = len(primer)
    if L == 0:
        return []
    floor = max(1, min(min_anneal, L))
    for n in range(L, floor - 1, -1):
        core = primer[L - n:]
        pat = reverse_complement(core) if on_reverse else core
        t = min(tp, floor) if clamp_at_floor else min(tp, n)
        # the primer's 3' end maps to the left of the window on the reverse strand and to the right of it
        # on the forward strand, because reverse-complementing swaps the ends
        three = frozenset(range(0, t)) if on_reverse else frozenset(range(n - t, n))
        hits = []
        for i in range(len(seq) - n + 1):
            mm, pos, ok = 0, [], True
            for k in range(n):
                if not primers._base_ok(pat[k], seq[i + k]):
                    mm += 1
                    if mm > max_mm or k in three:
                        ok = False
                        break
                    pos.append(k)
            if not ok:
                continue
            site = {"mm": mm, "mm_pos": pos, "tail5": L - n, "anneal_len": n}
            site["right" if on_reverse else "left"] = (i + n) if on_reverse else i
            hits.append(site)
        if hits:
            return hits
    return []


def canonical(sites):
    """Every field the two implementations must agree on, in a form that ignores enumeration order."""
    out = []
    for s in sites:
        side = "right" if "right" in s else "left"
        out.append((side, s[side], s["mm"], tuple(s["mm_pos"]), s["tail5"], s["anneal_len"]))
    return sorted(out)


# --------------------------------------------------------------------------------------------------
# Trial generation
# --------------------------------------------------------------------------------------------------
def _draw_params(rng, strict: bool):
    """One trial's parameters. `strict` forces tp > floor, the regime the clamps could differ in."""
    L = rng.randint(*PRIMER_LEN_RANGE)
    max_mm = rng.randint(*MAX_MM_RANGE)
    on_reverse = rng.random() < 0.5
    degenerate = rng.random() < 0.5
    if strict:
        min_anneal, tp = rng.choice(_STRICT_PAIRS)
    else:
        min_anneal = rng.randint(*MIN_ANNEAL_RANGE)
        tp = rng.randint(*TP_RANGE)
    return {"primer_len": L, "max_mm": max_mm, "tp": tp, "min_anneal": min_anneal,
            "on_reverse": on_reverse, "degenerate": degenerate,
            "seq_len": rng.randint(*SEQ_LEN_RANGE)}


def _random_primer(rng, L: int, degenerate: bool) -> str:
    """A concrete primer, or one carrying IUPAC ambiguity codes at a minority of positions."""
    out = []
    for _ in range(L):
        if degenerate and rng.random() < 0.2:
            out.append(rng.choice(DEGENERATE_CODES))
        else:
            out.append(rng.choice("ACGT"))
    return "".join(out)


def _instantiate(rng, core: str) -> str:
    """Resolve each (possibly degenerate) primer base to one template base it covers."""
    return "".join(rng.choice(primers._IUPAC[b]) for b in core)


def _plant(rng, chars: list, primer: str, p: dict) -> int:
    """Write annealing copies of 3'-proximal cores into the template.

    Purely random templates almost never carry a site long enough to test anything, so the comparison
    would be between two empty lists. Several copies are planted, at different core lengths and different
    mutation loads, so the candidates a trial produces reach different core lengths: that exercises the
    single-pass code's selection of the longest reachable core and its discard of the shorter candidates,
    which a single planted copy would not. Mutations are placed by distance from the 3' end, half the time
    biased into the strict window, so the 3' rule is exercised on both sides of its own boundary.
    """
    L = p["primer_len"]
    floor = max(1, min(p["min_anneal"], L))
    planted = 0
    for _ in range(rng.randint(1, 3)):
        n = rng.randint(max(1, floor - 2), L)          # occasionally below the floor: that copy yields no site
        frag = list(_instantiate(rng, primer[L - n:]))
        for _m in range(rng.randint(0, p["max_mm"] + 1)):
            d = rng.randint(0, min(p["tp"], n - 1)) if (rng.random() < 0.5 and p["tp"] > 0) \
                else rng.randint(0, n - 1)             # d = distance from the primer's 3' end
            k = d if p["on_reverse"] else n - 1 - d
            frag[k] = rng.choice([b for b in "ACGT" if b != frag[k]])
        frag = "".join(frag)
        if p["on_reverse"]:
            frag = reverse_complement(frag)
        if len(frag) > len(chars):
            continue
        at = rng.randint(0, len(chars) - len(frag))
        chars[at:at + len(frag)] = list(frag)
        planted += 1
    return planted


def run(trials: int = TRIALS, seed: int = SEED, strict_share: float = STRICT_SHARE) -> dict:
    """Run the comparison and return the result record. Writes nothing."""
    rng = random.Random(seed)
    t0 = time.time()
    n_strict = int(round(trials * strict_share))
    schedule = [True] * n_strict + [False] * (trials - n_strict)
    rng.shuffle(schedule)

    disagreements, examples = 0, []
    with_sites = both_empty = sites_compared = 0
    strict_trials = strict_with_sites = 0
    capped_trials = trials_with_tail = 0
    sensitivity_caught = strict_sensitivity_caught = tp_gt_floor = 0
    anneal_lens, tail_lens = {}, {}
    # drawn-parameter coverage, so "spanning both strands, degenerate and concrete primers, and the full
    # parameter range" is a countable property of the run rather than a property of the declared ranges
    cover = {"reverse_strand": 0, "forward_strand": 0, "degenerate_primer": 0, "concrete_primer": 0,
             "by_max_mm": {}, "by_tp": {}}

    for idx, strict in enumerate(schedule):
        p = _draw_params(rng, strict)
        cover["reverse_strand" if p["on_reverse"] else "forward_strand"] += 1
        cover["degenerate_primer" if p["degenerate"] else "concrete_primer"] += 1
        cover["by_max_mm"][str(p["max_mm"])] = cover["by_max_mm"].get(str(p["max_mm"]), 0) + 1
        cover["by_tp"][str(p["tp"])] = cover["by_tp"].get(str(p["tp"]), 0) + 1
        primer = _random_primer(rng, p["primer_len"], p["degenerate"])
        chars = [rng.choice("ACGT") for _ in range(p["seq_len"])]
        _plant(rng, chars, primer, p)
        seq = "".join(chars)

        got, capped = primers.anchored_sites(primer, seq, p["max_mm"], p["tp"],
                                             p["min_anneal"], p["on_reverse"])
        want = descending_sites(primer, seq, p["max_mm"], p["tp"], p["min_anneal"], p["on_reverse"])
        capped_trials += int(capped)

        floor = max(1, min(p["min_anneal"], p["primer_len"]))
        if strict:
            strict_trials += 1
            if want:
                strict_with_sites += 1
        if want or got:
            with_sites += 1
        else:
            both_empty += 1
        sites_compared += len(want)
        for s in want:
            anneal_lens[s["anneal_len"]] = anneal_lens.get(s["anneal_len"], 0) + 1
            tail_lens[s["tail5"]] = tail_lens.get(s["tail5"], 0) + 1
        if any(s["tail5"] > 0 for s in want):          # the winning core was shorter than the whole primer,
            trials_with_tail += 1                      # so the descent/extension path did work in this trial

        # Sensitivity control: the same trial run against a reference whose strict window is clamped once,
        # at the shortest core. Trials that catch that error are trials capable of catching a clamp defect.
        if p["tp"] > floor:
            tp_gt_floor += 1
            broken = descending_sites(primer, seq, p["max_mm"], p["tp"], p["min_anneal"],
                                      p["on_reverse"], clamp_at_floor=True)
            if canonical(broken) != canonical(want):
                sensitivity_caught += 1
                if strict:
                    strict_sensitivity_caught += 1

        a, b = canonical(got), canonical(want)
        if a != b:
            disagreements += 1
            if len(examples) < 20:                     # enough to see whether they share a regime
                examples.append({"trial": idx, "strict_stratum": strict, "floor": floor,
                                 "tp_exceeds_floor": p["tp"] > floor, "params": dict(p),
                                 "primer": primer, "seq": seq,
                                 "single_pass": [list(x) for x in a],
                                 "descending": [list(x) for x in b]})

    # A generator that stops planting would leave both implementations returning empty lists and the
    # harness green for the wrong reason, so the site yield is a reported result and a failure condition.
    yield_frac = with_sites / trials if trials else 0.0
    return {
        "trials": trials,
        "disagreements": disagreements,
        "agreed": trials - disagreements,
        "seed": seed,
        "strata": {"strict_window_wider_than_shortest_core": strict_trials,
                   "uniform": trials - strict_trials,
                   "strict_share_requested": strict_share},
        "informative": {"trials_with_sites": with_sites,
                        "both_empty": both_empty,
                        "site_yield": round(yield_frac, 4),
                        "sites_compared": sites_compared,
                        "trials_where_the_winning_core_was_shorter_than_the_primer": trials_with_tail,
                        "strict_stratum_trials_with_sites": strict_with_sites},
        "sensitivity_control": {"description": "trials on which a reference clamping the strict 3' window "
                                               "once, at the shortest permitted core, gives a different "
                                               "answer from the correct per-core-length clamp",
                                "trials_tested": tp_gt_floor,
                                "trials_caught": sensitivity_caught,
                                "strict_stratum_caught": strict_sensitivity_caught},
        "core_length_histogram": {str(k): anneal_lens[k] for k in sorted(anneal_lens)},
        "tail_length_histogram": {str(k): tail_lens[k] for k in sorted(tail_lens)},
        "parameter_ranges": {"max_mm": list(MAX_MM_RANGE), "tp": list(TP_RANGE),
                             "min_anneal": list(MIN_ANNEAL_RANGE),
                             "primer_len": list(PRIMER_LEN_RANGE),
                             "seq_len": list(SEQ_LEN_RANGE),
                             "strands": ["forward", "reverse"],
                             "primers": ["concrete", "IUPAC-degenerate"]},
        "coverage": {**cover, "by_max_mm": {k: cover["by_max_mm"][k] for k in sorted(cover["by_max_mm"])},
                     "by_tp": {k: cover["by_tp"][k] for k in sorted(cover["by_tp"], key=int)}},
        "site_cap_hits": capped_trials,
        "module_constants": {"MIN_PRIMER_SIZE": primers.MIN_PRIMER_SIZE,
                             "MAX_SITES": primers.MAX_SITES},
        "degenerate_codes": list(DEGENERATE_CODES),
        "examples": examples,
        "site_yield_ok": yield_frac >= 0.5,
        "sensitivity_ok": sensitivity_caught > 0,
        "teagle_version": TEAGLE_VERSION,
        "python": sys.version.split()[0],
        "elapsed_s": round(time.time() - t0, 1),
        "generated": time.strftime("%Y-%m-%dT%H:%M:%S"),
    }


def main():
    r = run()
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(r, open(OUT, "w", encoding="utf-8"), indent=1, sort_keys=True)
    inf = r["informative"]
    print(f"{r['trials']} trials, seed {r['seed']}: {r['disagreements']} disagreements")
    print(f"  strata: {r['strata']['strict_window_wider_than_shortest_core']} with the strict window wider "
          f"than the shortest core ({inf['strict_stratum_trials_with_sites']} of them yielded sites), "
          f"{r['strata']['uniform']} uniform")
    print(f"  {inf['trials_with_sites']} trials yielded sites ({inf['site_yield']:.1%}), "
          f"{inf['sites_compared']} sites compared; {inf['both_empty']} trials empty on both sides")
    sc = r["sensitivity_control"]
    print(f"  sensitivity control: a strict window clamped once at the shortest core is caught on "
          f"{sc['trials_caught']} of the {sc['trials_tested']} trials with tp > floor")
    print(f"  site cap reached in {r['site_cap_hits']} trials")
    print(f"  {r['elapsed_s']} s")
    for e in r["examples"]:
        print(f"  DISAGREE trial {e['trial']}: floor={e['floor']} tp={e['params']['tp']} "
              f"max_mm={e['params']['max_mm']} L={e['params']['primer_len']} "
              f"on_reverse={e['params']['on_reverse']} tp>floor={e['tp_exceeds_floor']}")
    print(f"-> {OUT}")
    bad = 0
    if not r["site_yield_ok"]:
        print("site yield collapsed: the trials are not testing anything; treat agreement as unproven")
        bad = 1
    if not r["sensitivity_ok"]:
        print("no trial distinguishes a mis-clamped strict window: the equivalence is unproven on the "
              "one axis the historical defect lived on")
        bad = 1
    return 1 if (bad or r["disagreements"]) else 0


if __name__ == "__main__":
    sys.exit(main())
