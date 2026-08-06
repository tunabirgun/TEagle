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


# The shortest oligo the tool will accept as a primer. It is the Primer3 minimum used for design and, in
# in-silico PCR, the shortest 3'-proximal core that may be called a binding site. One constant for both, so
# the tool cannot report a product from an annealing footprint shorter than a primer it would design.
MIN_PRIMER_SIZE = 18


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
        "PRIMER_MIN_SIZE": p.get("min_size", MIN_PRIMER_SIZE),
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


# IUPAC ambiguity: a primer base (possibly degenerate) matches a template base when the template base is one
# of the bases that code allows. Degenerate consensus primers are standard in TE work, and the genome-scan
# path (isPcr) is ambiguity-aware — local in-silico PCR must be too, or the two disagree for the same pair.
_IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "U": "T",
          "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
          "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}


def _base_ok(p: str, s: str) -> bool:
    if p == s:
        return True
    allowed = _IUPAC.get(p)                                # a degenerate primer base matches any template base it covers;
    return allowed is not None and s in allowed           # an ambiguous/unknown TEMPLATE base stays a mismatch (conservative)


def _scan(pattern: str, seq: str, max_mm: int, tp: int):
    """Return match positions of `pattern` on `seq` with <=max_mm mismatches. `pattern` may carry IUPAC
    ambiguity codes (degenerate primers); the template is plain ACGT. Each hit: (start, mismatches, mm_positions)."""
    L = len(pattern)
    hits = []
    for i in range(len(seq) - L + 1):
        mm, pos = 0, []
        for k in range(L):
            if not _base_ok(pattern[k], seq[i + k]):
                mm += 1
                pos.append(k)
                if mm > max_mm:
                    break
        if mm <= max_mm:
            hits.append((i, mm, pos))
    return hits


MAX_SITES, MAX_AMPS = 4000, 4000              # bound work + memory on repetitive templates


def _match_at(pattern, seq, at):
    """Mismatch count and positions for `pattern` placed at `at`, or None if it runs off the template."""
    if at < 0 or at + len(pattern) > len(seq):
        return None
    mm, pos = 0, []
    for k, pb in enumerate(pattern):
        if not _base_ok(pb, seq[at + k]):
            mm += 1
            pos.append(k)
    return mm, pos


def anchored_sites(primer, seq, max_mm, tp, min_anneal, on_reverse):
    """Binding sites for the longest 3'-proximal core of `primer` that anneals, found in ONE pass.

    Trying the full primer, then one base shorter, and so on until something binds is the obvious
    implementation and costs a whole template scan per tail length — twenty-five scans of a 58 Mb
    chromosome for a 42-mer. It is also unnecessary: shortening a match cannot introduce a mismatch, so
    every site of a long core is already a site of its own 3' suffix. The shortest permitted core is
    therefore scanned once and each hit extended 5' while the mismatch budget holds. The longest
    extension reached anywhere is the core the assay uses, which is what the descending scan returns.

    Returns (sites, capped). `on_reverse` selects the strand the primer primes from; on that strand a 5'
    extension of the primer appends complemented bases to the RIGHT of the matched window, because
    reverse-complementing swaps the ends.
    """
    L = len(primer)
    if L == 0:
        return [], False
    floor = max(1, min(min_anneal, L))
    t = min(tp, floor)                                # the strict 3' window lies inside every core
    core_min = primer[L - floor:]
    pat = reverse_complement(core_min) if on_reverse else core_min
    three = set(range(0, t)) if on_reverse else set(range(floor - t, floor))
    capped = False
    cands = []
    for i, _mm, pos in _scan(pat, seq, max_mm, t):
        if any(pp in three for pp in pos):            # strict 3': no mismatch in the terminal t bases
            continue
        cands.append(i)
        if len(cands) >= MAX_SITES:
            capped = True
            break
    if not cands:
        return [], capped

    # For each candidate, the longest core it supports. Working in DISTANCE FROM THE 3' END makes both
    # constraints closed-form: a core of length n is acceptable when it carries at most `max_mm`
    # mismatches (all mismatches at distance < n) and none inside its strict 3' window, which spans
    # min(tp, n) bases. The window widens with n only while n < tp, so the second constraint caps n at
    # the nearest mismatch whenever that mismatch falls inside the window at all.
    reach = []
    for i in cands:
        dists = sorted((p if on_reverse else floor - 1 - p) for p in _match_at(pat, seq, i)[1])
        limit = L                                     # how far the template allows the core to extend
        for n in range(floor, L):
            b = primer[L - n - 1]
            q = (i + n) if on_reverse else (i - (n - floor) - 1)
            if q < 0 or q >= len(seq):
                limit = n
                break
            if not _base_ok(reverse_complement(b) if on_reverse else b, seq[q]):
                dists.append(n)                       # the base just added lies n from the 3' end
        n_budget = L if len(dists) <= max_mm else sorted(dists)[max_mm]
        nearest = min(dists) if dists else None
        n_strict = L if (nearest is None or nearest >= tp) else nearest
        reach.append(min(limit, n_budget, n_strict))

    viable = [(i, n) for i, n in zip(cands, reach) if n >= floor]
    if not viable:                                    # every candidate failed once the exact 3' rule was
        return [], capped                             # applied at its own core length
    best = max(n for _i, n in viable)
    core = primer[L - best:]
    pat_best = reverse_complement(core) if on_reverse else core
    out = []
    for i, n in viable:
        if n != best:
            continue
        at = i if on_reverse else i - (best - floor)
        m = _match_at(pat_best, seq, at)
        if m is None:                                 # the extension check bounds this; kept as a guard
            continue
        mm, pos = m
        site = {"mm": mm, "mm_pos": pos, "tail5": L - best, "anneal_len": best}
        site["right" if on_reverse else "left"] = (at + best) if on_reverse else at
        out.append(site)
    return out, capped


def in_silico_pcr(fwd: str, rev: str, seq: str, seq_id: str = "template",
                  max_mm: int = 2, tp: int = 5, prod_min: int = 70, prod_max: int = 1000,
                  target_span: list | None = None, stats: dict | None = None,
                  min_anneal: int = MIN_PRIMER_SIZE):
    """Pair-aware in-silico PCR. Searches both strands, applies a strict 3' rule
    (zero mismatches in the terminal `tp` bases), builds amplicons from inward-facing
    sites within [prod_min, prod_max] — both two-primer (F+R) products and single-primer
    (F+F / R+R) self-priming products across inverted repeats (marked single_primer).
    on_target and single_primer are mutually exclusive, so on/off/single counts are disjoint
    and exhaustive downstream (off = total - on - single can never go negative).
    A primer may carry a non-templated 5' tail - a restriction site, an adapter, a barcode, a promoter -
    which does not anneal on the first cycle but is incorporated into the product. Requiring the whole
    primer to match therefore rejects a large and ordinary class of real assays. Binding is instead
    anchored at the 3' end: the longest 3'-proximal suffix that matches is used, down to `min_anneal`
    bases, and any unmatched 5' remainder is reported as a tail rather than silently ignored. Product
    length is reported both as the templated span and as the full product including both tails, because
    those differ for a tailed pair and a user comparing against a gel needs the latter.

    Returns real amplicon list."""
    fwd, rev, seq = fwd.upper(), rev.upper(), seq.upper()     # case-insensitive: a lowercase primer must still bind
    max_mm, tp = max(0, int(max_mm)), max(0, int(tp))         # non-negative ints; tp<0 must not silently disable the 3' rule
    amps = []
    sites_capped = False                                      # a binding-site scan that hit MAX_SITES stopped early

    def sites(primer, on_reverse):
        nonlocal sites_capped
        found, capped = anchored_sites(primer, seq, max_mm, tp, min_anneal, on_reverse)
        sites_capped = sites_capped or capped
        return found

    fset = [("F", s) for s in sites(fwd, False)] + [("R", s) for s in sites(rev, False)]
    rset = [("F", s) for s in sites(fwd, True)] + [("R", s) for s in sites(rev, True)]
    capped = False
    for fo, fs in fset:
        for ro, rs in rset:
            left, right = fs["left"], rs["right"]
            plen = right - left
            t5f, t5r = fs.get("tail5", 0), rs.get("tail5", 0)
            # The size window is about the band, so it is applied to the product that is actually made -
            # templated span plus both incorporated 5' tails - not to the distance between binding sites.
            # For an untailed pair the two are identical and the behaviour is unchanged.
            if left < right and prod_min <= plen + t5f + t5r <= prod_max:
                single = fo == ro                     # same primer at both ends: self-priming across a TIR/LTR
                # a self-priming product is an artefact, never the intended amplicon (that is always F+R), so it
                # stays in the single-primer bucket even inside the target window -> on/off/single stay disjoint
                on = (not single) and bool(target_span and left >= target_span[0] - 5 and right <= target_span[1] + 5)
                # A non-templated 5' tail is incorporated into the product, so the band a user sees on a
                # gel is the templated span PLUS both tails. Both are reported: `length` stays the
                # templated span so coordinates and slices remain consistent, and `product_length` is
                # what the assay actually yields.
                amps.append({
                    "source": seq_id, "start": left, "end": right, "length": plen,
                    "product_length": plen + t5f + t5r,
                    "fwd_tail5": t5f, "rev_tail5": t5r,
                    "fwd_anneal_len": fs.get("anneal_len"), "rev_anneal_len": rs.get("anneal_len"),
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
    if stats is not None:
        # A truncated search must be distinguishable from an exhaustive one. Both caps bite exactly on the
        # templates where the answer matters most — a repetitive element, where "no further off-target
        # product" would otherwise be read as a specificity result rather than as the search giving up.
        stats["amplicons_capped"] = capped
        stats["sites_capped"] = sites_capped
        stats["max_amplicons"] = MAX_AMPS
        stats["max_sites"] = MAX_SITES
    return amps
