"""Real primer design (Primer3 via primer3-py) + pair-aware in-silico PCR.
No hand-rolled thermodynamics: candidate generation and Tm come from Primer3."""
from __future__ import annotations
import math
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
INTENDED_ENDPOINT_TOLERANCE_BP = 0
DESIGN_DEFAULTS = {
    "opt_size": 20, "min_size": MIN_PRIMER_SIZE, "max_size": 27,
    "opt_tm": 60.0, "min_tm": 57.0, "max_tm": 63.0,
    "min_gc": 40.0, "max_gc": 60.0,
    "prod_min": 120, "prod_max": 700, "num_return": 5,
    "max_poly_x": 4, "gc_clamp": 0,
}


def resolved_design_parameters(params: dict | None = None) -> dict:
    submitted = params or {}
    resolved = {name: submitted.get(name, default) for name, default in DESIGN_DEFAULTS.items()}
    integer_fields = ("opt_size", "min_size", "max_size", "prod_min", "prod_max",
                      "num_return", "max_poly_x", "gc_clamp")
    for name, value in resolved.items():
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
            raise ValueError(f"{name} must be a finite number")
        if name in integer_fields and not isinstance(value, int):
            raise ValueError(f"{name} must be an integer")
    if (resolved["min_size"] < MIN_PRIMER_SIZE or
            not resolved["min_size"] <= resolved["opt_size"] <= resolved["max_size"]):
        raise ValueError(f"primer sizes must satisfy {MIN_PRIMER_SIZE} <= minimum <= optimal <= maximum")
    if not resolved["min_tm"] <= resolved["opt_tm"] <= resolved["max_tm"]:
        raise ValueError("Tm values must satisfy minimum <= optimal <= maximum")
    if not 0 <= resolved["min_gc"] <= resolved["max_gc"] <= 100:
        raise ValueError("GC bounds must satisfy 0 <= minimum <= maximum <= 100")
    if resolved["prod_min"] < 1 or resolved["prod_max"] < resolved["prod_min"]:
        raise ValueError("product-size bounds must be positive and prod_min <= prod_max")
    if resolved["num_return"] < 1 or resolved["max_poly_x"] < 1 or resolved["gc_clamp"] < 0:
        raise ValueError("num_return/max_poly_x must be positive and gc_clamp non-negative")
    return resolved


def _coordinate(value, name):
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{name} coordinates must be integers")
    return value


def start_length_region(value, template_length: int, name: str):
    """Validate and normalize a 0-based ``[start, length]`` region."""
    if value is None:
        return None
    if not isinstance(value, (list, tuple)) or len(value) != 2:
        raise ValueError(f"{name} must be a [start, length] pair")
    start = _coordinate(value[0], name)
    length = _coordinate(value[1], name)
    if start < 0 or length <= 0 or start + length > template_length:
        raise ValueError(f"{name} must be a positive region within the template")
    return [start, length]


def endpoint_span(value, template_length: int, name: str):
    """Validate and normalize a 0-based half-open ``[start, end]`` span."""
    if value is None:
        return None
    if not isinstance(value, (list, tuple)) or len(value) != 2:
        raise ValueError(f"{name} must be a [start, end] pair")
    start = _coordinate(value[0], name)
    end = _coordinate(value[1], name)
    if start < 0 or end <= start or end > template_length:
        raise ValueError(f"{name} must be a positive half-open span within the template")
    return [start, end]


def endpoint_tolerance(value):
    value = _coordinate(value, "intended endpoint tolerance")
    if value < 0:
        raise ValueError("intended endpoint tolerance must be non-negative")
    return value


def design_primers(template: str, params: dict | None = None, target: list | None = None,
                   included: list | None = None):
    """Design primer pairs with Primer3. target=[start,len] forces the product to span a locus;
    included=[start,len] restricts primers to a region (e.g. a protein domain).
    Returns list of candidate dicts with real Tm/GC/product size/penalty."""
    if primer3 is None:
        raise RuntimeError("Primer3 is unavailable in this environment (" + (PRIMER3_ERROR or "import failed") + ")")
    target = start_length_region(target, len(template), "primer3_target")
    included = start_length_region(included, len(template), "included_region")
    if target is not None and included is not None:
        target_end = target[0] + target[1]
        included_end = included[0] + included[1]
        if target[0] < included[0] or target_end > included_end:
            raise ValueError("primer3_target must be contained in included_region")
    p = resolved_design_parameters(params)
    global_args = {
        "PRIMER_OPT_SIZE": p["opt_size"],
        "PRIMER_MIN_SIZE": p["min_size"],
        "PRIMER_MAX_SIZE": p["max_size"],
        "PRIMER_OPT_TM": p["opt_tm"],
        "PRIMER_MIN_TM": p["min_tm"],
        "PRIMER_MAX_TM": p["max_tm"],
        "PRIMER_MIN_GC": p["min_gc"],
        "PRIMER_MAX_GC": p["max_gc"],
        "PRIMER_PRODUCT_SIZE_RANGE": [p["prod_min"], p["prod_max"]],
        "PRIMER_NUM_RETURN": p["num_return"],
        "PRIMER_MAX_POLY_X": p["max_poly_x"],
        "PRIMER_GC_CLAMP": p["gc_clamp"],
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
        left_pos = [lpos[0], lpos[0] + lpos[1]]
        right_pos = [rpos[0] - rpos[1] + 1, rpos[0] + 1]
        intended_product_span = endpoint_span(
            [left_pos[0], right_pos[1]], len(template), "intended_product_span")
        product_size = r[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"]
        if product_size != intended_product_span[1] - intended_product_span[0]:
            raise RuntimeError("Primer3 product size disagrees with its primer coordinates")
        if template[left_pos[0]:left_pos[1]].upper() != r[f"PRIMER_LEFT_{i}_SEQUENCE"].upper():
            raise RuntimeError("Primer3 forward-primer sequence disagrees with its template coordinates")
        right_template = template[right_pos[0]:right_pos[1]]
        if reverse_complement(right_template).upper() != r[f"PRIMER_RIGHT_{i}_SEQUENCE"].upper():
            raise RuntimeError("Primer3 reverse-primer sequence disagrees with its template coordinates")
        if target is not None:
            target_end = target[0] + target[1]
            if intended_product_span[0] > target[0] or intended_product_span[1] < target_end:
                raise RuntimeError("Primer3 candidate does not span primer3_target")
        if included is not None:
            included_end = included[0] + included[1]
            if intended_product_span[0] < included[0] or intended_product_span[1] > included_end:
                raise RuntimeError("Primer3 candidate leaves included_region")
        out.append({
            "id": f"P{i+1}",
            "left_seq": r[f"PRIMER_LEFT_{i}_SEQUENCE"],
            "right_seq": r[f"PRIMER_RIGHT_{i}_SEQUENCE"],
            "left_pos": left_pos,
            "right_pos": right_pos,
            "left_tm": round(r[f"PRIMER_LEFT_{i}_TM"], 1),
            "right_tm": round(r[f"PRIMER_RIGHT_{i}_TM"], 1),
            "left_gc": round(r[f"PRIMER_LEFT_{i}_GC_PERCENT"], 1),
            "right_gc": round(r[f"PRIMER_RIGHT_{i}_GC_PERCENT"], 1),
            "product_size": product_size,
            "penalty": round(r[f"PRIMER_PAIR_{i}_PENALTY"], 2),
            "primer3_target": target,
            "included_region": included,
            "intended_product_span": intended_product_span,
            "intended_endpoint_tolerance_bp": INTENDED_ENDPOINT_TOLERANCE_BP,
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
                  primer3_target: list | None = None, intended_product_span: list | None = None,
                  intended_endpoint_tolerance_bp: int = INTENDED_ENDPOINT_TOLERANCE_BP,
                  assess_design_geometry: bool = True,
                  stats: dict | None = None,
                  min_anneal: int = MIN_PRIMER_SIZE):
    """Pair-aware in-silico PCR. Searches both strands, applies a strict 3' rule
    (zero mismatches in the terminal `tp` bases), builds amplicons from inward-facing
    sites within [prod_min, prod_max] — both two-primer (F+R) products and single-primer
    (F+F / R+R) self-priming products across inverted repeats (marked single_primer).
    ``matches_intended_product`` identifies the designed endpoint pair within the explicit endpoint
    tolerance. ``covers_primer3_target`` separately reports whether a two-primer product spans the
    Primer3 target. ``on_target`` is their conjunction when a Primer3 target exists; without a target
    it is unassessed (None), not an off-target call.
    A primer may carry a non-templated 5' tail - a restriction site, an adapter, a barcode, a promoter -
    which does not anneal on the first cycle but is incorporated into the product. Requiring the whole
    primer to match therefore rejects a large and ordinary class of real assays. Binding is instead
    anchored at the 3' end: the longest 3'-proximal suffix that matches is used, down to `min_anneal`
    bases, and any unmatched 5' remainder is reported as a tail rather than silently ignored. Product
    length is reported both as the templated span and as the full product including both tails, because
    those differ for a tailed pair and a user comparing against a gel needs the latter.

    Returns real amplicon list."""
    fwd, rev, seq = fwd.upper(), rev.upper(), seq.upper()     # case-insensitive: a lowercase primer must still bind
    for name, value in (("max_mm", max_mm), ("tp", tp), ("min_anneal", min_anneal),
                        ("prod_min", prod_min), ("prod_max", prod_max)):
        if isinstance(value, bool) or not isinstance(value, int):
            raise ValueError(f"{name} must be an integer")
    if max_mm < 0 or tp < 0:
        raise ValueError("max_mm and tp must be non-negative")
    if min_anneal < MIN_PRIMER_SIZE:
        raise ValueError(f"min_anneal must be at least {MIN_PRIMER_SIZE}")
    if prod_min < 1 or prod_max < prod_min:
        raise ValueError("product-size bounds must be positive and prod_min <= prod_max")
    if len(fwd) < min_anneal or len(rev) < min_anneal:
        raise ValueError("each primer must be at least min_anneal bases long")
    shortest_primer = min(len(fwd), len(rev))
    if tp > shortest_primer or max_mm > shortest_primer:
        raise ValueError("tp and max_mm cannot exceed the shorter primer length")
    primer3_target = start_length_region(primer3_target, len(seq), "primer3_target")
    intended_product_span = endpoint_span(intended_product_span, len(seq), "intended_product_span")
    tolerance = endpoint_tolerance(intended_endpoint_tolerance_bp)
    if not isinstance(assess_design_geometry, bool):
        raise ValueError("assess_design_geometry must be boolean")
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
                full_product_length = plen + t5f + t5r
                left_oligo = fwd if fo == "F" else rev
                right_oligo = rev if ro == "R" else fwd
                product_seq = (left_oligo[:t5f] + seq[left:right] +
                               reverse_complement(right_oligo[:t5r]))
                matches_intended = None
                if assess_design_geometry and intended_product_span is not None:
                    matches_intended = (fo == "F" and ro == "R" and
                                        abs(left - intended_product_span[0]) <= tolerance and
                                        abs(right - intended_product_span[1]) <= tolerance)
                covers_target = None
                if assess_design_geometry and primer3_target is not None:
                    target_end = primer3_target[0] + primer3_target[1]
                    covers_target = not single and left <= primer3_target[0] and right >= target_end
                on = ((matches_intended is True and covers_target is True)
                      if assess_design_geometry and primer3_target is not None else None)
                # A non-templated 5' tail is incorporated into the product, so the band a user sees on a
                # gel is the templated span PLUS both tails. Both are reported: `length` stays the
                # templated span so coordinates and slices remain consistent, and `product_length` is
                # what the assay actually yields.
                amps.append({
                    "source": seq_id, "start": left, "end": right, "length": plen,
                    "product_key": f"{seq_id}:{left}:{right}:{full_product_length}",
                    "mechanism_key": f"{seq_id}:{fo}+{ro}:{left}:{right}:{full_product_length}",
                    "product_length": full_product_length,
                    "fwd_tail5": t5f, "rev_tail5": t5r,
                    "left_tail5": t5f, "right_tail5": t5r,
                    "fwd_anneal_len": fs.get("anneal_len"), "rev_anneal_len": rs.get("anneal_len"),
                    "fwd_primer": fo, "rev_primer": ro, "left_primer": fo, "right_primer": ro,
                    "primer_orientation": f"{fo}+{ro}", "single_primer": single,
                    "fwd_mm": fs["mm"], "rev_mm": rs["mm"],
                    "matches_intended_product": matches_intended,
                    "covers_primer3_target": covers_target,
                    "on_target": on,
                    "amplicon_5p": seq[left:left + 30] + ("…" if plen > 60 else ""),
                    "templated_seq": seq[left:right], "seq": product_seq,
                })
                if len(amps) >= MAX_AMPS:
                    capped = True
                    break
        if capped:
            break
    # intended endpoint pair first, then target-covering products, then strongest priming and position
    amps.sort(key=lambda a: (a["matches_intended_product"] is not True,
                             a["covers_primer3_target"] is not True,
                             a["fwd_mm"] + a["rev_mm"], a["single_primer"], a["start"]))
    if stats is not None:
        # A truncated search must be distinguishable from an exhaustive one. Both caps bite exactly on the
        # templates where the answer matters most — a repetitive element, where "no further off-target
        # product" would otherwise be read as a specificity result rather than as the search giving up.
        stats["amplicons_capped"] = capped
        stats["sites_capped"] = sites_capped
        stats["max_amplicons"] = MAX_AMPS
        stats["max_sites"] = MAX_SITES
    return amps
