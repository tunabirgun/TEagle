"""Real primer design (Primer3 via primer3-py) + pair-aware in-silico PCR.
No hand-rolled thermodynamics: candidate generation and Tm come from Primer3."""
from __future__ import annotations
from .sequtil import reverse_complement

try:                                                # a broken/missing Primer3 must not crash the engine —
    import primer3                                  # primer design reports unavailable; in-silico PCR (pure Python) still runs
    PRIMER3_VERSION = primer3.__version__
    PRIMER3_ERROR = None
except Exception as _e:
    primer3 = None
    PRIMER3_VERSION = "unavailable"
    PRIMER3_ERROR = f"{type(_e).__name__}: {_e}"


def design_primers(template: str, params: dict | None = None, target: list | None = None,
                   included: list | None = None):
    """Design primer pairs with Primer3. target=[start,len] forces the product to span a locus;
    included=[start,len] restricts primers to a region (e.g. a protein domain).
    Returns list of candidate dicts with real Tm/GC/product size/penalty."""
    if primer3 is None:
        raise RuntimeError("Primer3 is unavailable in this environment (" + (PRIMER3_ERROR or "import failed") + ")")
    p = params or {}
    global_args = {
        "PRIMER_OPT_SIZE": p.get("opt_size", 20),
        "PRIMER_MIN_SIZE": p.get("min_size", 18),
        "PRIMER_MAX_SIZE": p.get("max_size", 27),
        "PRIMER_OPT_TM": p.get("opt_tm", 60.0),
        "PRIMER_MIN_TM": p.get("min_tm", 57.0),
        "PRIMER_MAX_TM": p.get("max_tm", 63.0),
        "PRIMER_MIN_GC": p.get("min_gc", 40.0),
        "PRIMER_MAX_GC": p.get("max_gc", 60.0),
        "PRIMER_PRODUCT_SIZE_RANGE": [p.get("prod_min", 120), p.get("prod_max", 700)],
        "PRIMER_NUM_RETURN": p.get("num_return", 5),
        "PRIMER_MAX_POLY_X": p.get("max_poly_x", 4),
        "PRIMER_GC_CLAMP": p.get("gc_clamp", 0),
    }
    seq_args = {"SEQUENCE_ID": "teagle", "SEQUENCE_TEMPLATE": template}
    if target:
        seq_args["SEQUENCE_TARGET"] = target
    if included:
        seq_args["SEQUENCE_INCLUDED_REGION"] = included
    r = primer3.bindings.design_primers(seq_args, global_args)
    n = r.get("PRIMER_PAIR_NUM_RETURNED", 0)
    out = []
    for i in range(n):
        lpos = r[f"PRIMER_LEFT_{i}"]           # (start, len)
        rpos = r[f"PRIMER_RIGHT_{i}"]
        out.append({
            "id": f"P{i+1}",
            "left_seq": r[f"PRIMER_LEFT_{i}_SEQUENCE"],
            "right_seq": r[f"PRIMER_RIGHT_{i}_SEQUENCE"],
            "left_pos": [lpos[0], lpos[0] + lpos[1]],
            "right_pos": [rpos[0] - rpos[1] + 1, rpos[0] + 1],
            "left_tm": round(r[f"PRIMER_LEFT_{i}_TM"], 1),
            "right_tm": round(r[f"PRIMER_RIGHT_{i}_TM"], 1),
            "left_gc": round(r[f"PRIMER_LEFT_{i}_GC_PERCENT"], 1),
            "right_gc": round(r[f"PRIMER_RIGHT_{i}_GC_PERCENT"], 1),
            "product_size": r[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"],
            "penalty": round(r[f"PRIMER_PAIR_{i}_PENALTY"], 2),
        })
    return {"candidates": out, "explain_left": r.get("PRIMER_LEFT_EXPLAIN", ""),
            "explain_right": r.get("PRIMER_RIGHT_EXPLAIN", ""),
            "explain_pair": r.get("PRIMER_PAIR_EXPLAIN", "")}


def _scan(pattern: str, seq: str, max_mm: int, tp: int):
    """Return match positions of `pattern` on `seq` with <=max_mm mismatches.
    Each hit: (start, mismatches, mm_positions, three_prime_end_exact_flag_left, ...)."""
    L = len(pattern)
    hits = []
    for i in range(len(seq) - L + 1):
        mm, pos = 0, []
        for k in range(L):
            if pattern[k] != seq[i + k]:
                mm += 1
                pos.append(k)
                if mm > max_mm:
                    break
        if mm <= max_mm:
            hits.append((i, mm, pos))
    return hits


def in_silico_pcr(fwd: str, rev: str, seq: str, seq_id: str = "template",
                  max_mm: int = 2, tp: int = 5, prod_min: int = 70, prod_max: int = 1000,
                  target_span: list | None = None):
    """Pair-aware in-silico PCR. Searches both strands, applies a strict 3' rule
    (zero mismatches in the terminal `tp` bases), builds amplicons from inward-facing
    sites within [prod_min, prod_max] — both two-primer (F+R) products and single-primer
    (F+F / R+R) self-priming products across inverted repeats (marked single_primer).
    Returns real amplicon list."""
    fwd, rev, seq = fwd.upper(), rev.upper(), seq.upper()     # case-insensitive: a lowercase primer must still bind
    max_mm, tp = max(0, int(max_mm)), max(0, int(tp))         # non-negative ints; tp<0 must not silently disable the 3' rule
    MAX_SITES, MAX_AMPS = 4000, 4000                          # bound work + memory on repetitive templates
    amps = []
    # forward-capable binding: primer matches top strand 5'->3'; 3' end = right side
    def forward_sites(primer):
        t = min(tp, len(primer))                              # clamp so a huge tp cannot build an enormous index set
        sites = []
        for i, mm, pos in _scan(primer, seq, max_mm, t):
            three = set(range(len(primer) - t, len(primer)))
            if any(pp in three for pp in pos):        # strict 3': no mismatch in last t
                continue
            sites.append({"left": i, "mm": mm, "mm_pos": pos})
            if len(sites) >= MAX_SITES:
                break
        return sites
    # reverse-capable binding: revcomp(primer) matches top strand; primer 3' = left side
    def reverse_sites(primer):
        t = min(tp, len(primer))
        rc = reverse_complement(primer)
        sites = []
        for i, mm, pos in _scan(rc, seq, max_mm, t):
            three = set(range(0, t))                  # primer 3' maps to left of rc window
            if any(pp in three for pp in pos):
                continue
            sites.append({"right": i + len(rc), "mm": mm, "mm_pos": pos})
            if len(sites) >= MAX_SITES:
                break
        return sites

    fset = [("F", s) for s in forward_sites(fwd)] + [("R", s) for s in forward_sites(rev)]
    rset = [("F", s) for s in reverse_sites(fwd)] + [("R", s) for s in reverse_sites(rev)]
    capped = False
    for fo, fs in fset:
        for ro, rs in rset:
            left, right = fs["left"], rs["right"]
            plen = right - left
            if left < right and prod_min <= plen <= prod_max:
                single = fo == ro                     # same primer at both ends: self-priming across a TIR/LTR
                on = bool(target_span and left >= target_span[0] - 5 and right <= target_span[1] + 5)
                amps.append({
                    "source": seq_id, "start": left, "end": right, "length": plen,
                    "fwd_primer": fo, "rev_primer": ro, "single_primer": single,
                    "fwd_mm": fs["mm"], "rev_mm": rs["mm"],
                    "on_target": on,
                    "amplicon_5p": seq[left:left + 30] + ("…" if plen > 60 else ""),
                    "seq": seq[left:right],
                })
                if len(amps) >= MAX_AMPS:
                    capped = True
                    break
        if capped:
            break
    # on-target first, then strongest priming (fewest mismatches), then two-primer before self-priming, then position
    amps.sort(key=lambda a: (not a["on_target"], a["fwd_mm"] + a["rev_mm"], a["single_primer"], a["start"]))
    return amps
