"""Real structural TE evidence (Layer D): terminal direct repeats (LTR), terminal
inverted repeats (TIR), poly-A tails, flanking target-site duplications (TSD).
Heuristic detectors — reported as candidate structural evidence, never as family/DB calls."""
from __future__ import annotations
from collections import Counter
from .sequtil import reverse_complement


def _is_simple(s: str) -> bool:
    """Low-complexity guard: a homopolymer or short-period tandem repeat satisfies any identity test
    trivially and is not terminal-repeat evidence. Flags a copy dominated by one base or built from
    very few distinct 4-mers (a satellite/microsatellite signature, not an LTR)."""
    if len(s) < 8:
        return True
    u = s.upper()
    if max(u.count(c) for c in "ACGT") / len(u) > 0.7:
        return True
    quads = {u[i:i + 4] for i in range(len(u) - 3)}
    # absolute floor as well as proportional: a short (TA)n tract yields 2 distinct 4-mers, which clears
    # any percentage of a small window but is still a microsatellite, not a terminal repeat.
    return len(quads) < max(4, 0.05 * (len(u) - 3))


def _kmer_index(a: str, k: int):
    idx = {}
    for i in range(len(a) - k + 1):
        kmer = a[i:i + k]
        if "N" in kmer or "n" in kmer:              # an N-run is not repeat evidence; never seed anchors on it
            continue
        idx.setdefault(kmer, []).append(i)
    return idx


def _identity(r5: str, r3: str) -> float:
    m = min(len(r5), len(r3))
    if m == 0:
        return 0.0
    # a match needs two KNOWN equal bases; N is unknown, so N==N is not base-pairing evidence
    return round(100 * sum(1 for x, y in zip(r5, r3) if x == y and x not in "Nn") / m, 1)


def _termini_motif(seq: str, l5s: int, l5e: int, l3s: int, l3e: int) -> dict:
    """Report whether the LTR termini carry the canonical 5'-TG…CA-3' dinucleotides.

    A BADGE, never a boundary mover. The motif is characteristic of many LTR retrotransposons but is not
    universal — measured here, copia (M11240) satisfies it exactly while gypsy/mdg4 (M12927) reads AG…TT
    and no offset within +/-10 bp satisfies it on both copies. Snapping boundaries onto the motif would
    therefore have MOVED a correct gypsy boundary to satisfy a rule that element does not follow.
    Tolerance of one mismatch across the four positions follows LTRharvest's -motifmis 1
    (Ellinghaus et al. 2008 BMC Bioinformatics 9:18)."""
    obs = {"five_start": seq[l5s:l5s + 2].upper(), "five_end": seq[l5e - 2:l5e].upper(),
           "three_start": seq[l3s:l3s + 2].upper(), "three_end": seq[l3e - 2:l3e].upper()}
    want = {"five_start": "TG", "five_end": "CA", "three_start": "TG", "three_end": "CA"}
    mism = sum(1 for k, v in want.items() if obs[k] != v)
    return {**obs, "canonical": mism == 0, "mismatches": mism,
            "within_ltrharvest_tolerance": mism <= 1,
            "note": ("termini match the canonical TG…CA" if mism == 0 else
                     "termini do not match TG…CA; the motif is common but not universal among LTR "
                     "retrotransposons, so its absence is not evidence against the call")}


def _extend(seq: str, i: int, j: int, step: int, limit: int, x_drop: int = 30) -> int:
    """Ungapped X-drop extension of a repeat pair along a fixed diagonal (BLAST-style: Altschul 1990).
    Walks (i, j) outward by `step` (+1 right, -1 left) and returns the offset of the best-scoring
    prefix — how far the pair genuinely stays a repeat, not how far it can be walked. Match +1 /
    mismatch -1 separates an 80%-identity LTR pair (+0.6 per base) from random (-0.5 per base)."""
    score = best = best_off = 0
    for off in range(max(0, limit)):
        a, b = seq[i + step * off], seq[j + step * off]
        if a == b and a not in "Nn":                 # N is unknown, never base-pairing evidence
            score += 1
        else:
            score -= 2                               # break-even at 67% identity: a diverged pair still
        if score > best:                             # extends, a chance match past the end cannot set a best
            best, best_off = score, off + 1
        elif score < best - x_drop:
            break
    return best_off


# Seeding window for the terminal-repeat scan. It must EXCEED the repeat being sought: the two windows
# only share k-mers where they overlap in repeat-internal coordinates, so a window at or below the true
# LTR length seeds nothing. 6 kb covers the longest LTRs described (multi-kb plant Ty3/Ogre LTRs).
_LTR_SEED_WINDOW = 6000


def find_ltr(seq: str, k: int = 13, min_ltr: int = 80, min_anchors: int = 4):
    """Detect a terminal DIRECT repeat pair (LTR candidate). Coords 0-based half-open.
    Anchors only SEED the repeat — its extent is measured by X-drop extension along the seeded
    diagonal over the whole sequence, so the reported length is not an artefact of the window."""
    n = len(seq)
    W = min(n // 2, _LTR_SEED_WINDOW)
    if W < min_ltr:
        return None
    a, b = seq[:W], seq[n - W:]                     # b offset by (n-W) in full seq
    idx = _kmer_index(a, k)
    diag, pts = Counter(), {}
    for j in range(len(b) - k + 1):
        for i in idx.get(b[j:j + k], ()):
            p5, p3 = i, (n - W) + j
            if p3 <= p5 + min_ltr:                   # 3' copy must be downstream
                continue
            d = p3 - p5
            diag[d] += 1
            pts.setdefault(d, []).append((p5, p3))
    if not diag:
        return None
    best_d, _ = max(diag.items(), key=lambda kv: kv[1])
    grp = [pt for dd, pl in pts.items() if abs(dd - best_d) <= 3 for pt in pl]
    if len(grp) < min_anchors:
        return None
    p5s = [p for p, _ in grp]
    l5s, l5e = min(p5s), max(p5s) + k
    L = l5e - l5s
    if L < min_ltr:
        return None
    # Anchor density floor. This only has to reject a handful of chance k-mers sharing a diagonal; the
    # real guard is the identity of the EXTENDED copies below. The old 40% floor demanded that 40% of
    # positions start an exact k-mer, which needs ~93% identity (0.93^13 = 0.4) — so it silently dropped
    # every LTR pair older than a few million years and made the 80% identity gate unreachable. At 75%
    # identity exact 13-mers still start at 2.4% of positions, so 1% keeps a wide margin.
    if len(set(p5s)) < max(min_anchors, 0.01 * (L - k + 1)):
        return None
    l3s, l3e = best_d + l5s, best_d + l5e
    if l3e > n:
        l3e = n
    # Measure the repeat by extending the seeded diagonal across the FULL sequence. Without this the
    # length is whatever fraction of the LTR happened to fall inside both seeding windows, reported as
    # if it were measured. Neither copy may run into the other, nor past the sequence ends.
    gap = l3s - l5e                                  # internal region: neither copy may extend into the other
    lim_l, lim_r = min(l5s, gap), min(n - l3e, gap)
    ext_l = _extend(seq, l5s - 1, l3s - 1, -1, lim_l)
    ext_r = _extend(seq, l5e, l3e, +1, lim_r)
    l5s, l3s = l5s - ext_l, l3s - ext_l
    l5e, l3e = l5e + ext_r, l3e + ext_r
    L = l5e - l5s
    # A repeat that reaches a record edge, or that was stopped only by running into its own other copy,
    # may continue past what was supplied: the length is then a LOWER BOUND and must never be presented
    # as a measurement. Touching the edge counts even when no extension was needed to get there.
    bounded = (l5s == 0 or l3e == n
               or (ext_l > 0 and ext_l == lim_l) or (ext_r > 0 and ext_r == lim_r))
    ident = _identity(seq[l5s:l5e], seq[l3s:l3e])
    if ident < 80 or _is_simple(seq[l5s:l5e]):       # identity now carries the weight the density floor used to
        return None
    ev = {"type": "LTR (terminal direct repeat)", "ltr_len": L, "identity": ident,
          "five_prime": [l5s, l5e], "three_prime": [l3s, l3e],
          "element_span": [l5s, l3e],
          "method": "k-mer seed (k=%d, %d bp window) + X-drop diagonal extension" % (k, W)}
    ev["termini"] = _termini_motif(seq, l5s, l5e, l3s, l3e)
    if bounded:
        ev["length_is_lower_bound"] = True
        ev["bound_reason"] = ("the repeat still matched where the record ends, so the copies may extend "
                              "beyond it — the length is a lower bound, not a measurement")
    return ev


def _terminal_tir(seq: str, min_tir: int, max_tir: int, min_ident: int = 80):
    """Direct scan for a terminal inverted repeat anchored at the true element ends [0,L)/[n-L,n).
    A real TIR is often imperfect, so the length is the LONGEST window whose identity stays above
    threshold before the identity 'cliff' where the match becomes random — not the shortest, most
    perfect core. Catches short canonical TIRs (11 bp Ac) and longer imperfect ones (28 bp mariner)."""
    n = len(seq)
    hi = min(max_tir, n // 2)
    if hi < min_tir:
        return None
    b_all = reverse_complement(seq[n - hi:])             # revcomp of the 3' window; its first L bases = revcomp(seq[n-L:])
    # `matches` accumulates across L (b_all[:L] is a true prefix for every L), so the whole scan is one
    # pass instead of a re-count per window — which is what makes a realistic max_tir affordable.
    best_L, best_score, best_ident = None, 0, 0.0
    matches = 0
    for L in range(1, hi + 1):
        x, y = seq[L - 1], b_all[L - 1]
        if x == y and x not in "Nn":                     # N is not base-pairing evidence
            matches += 1
        if L < min_tir:
            continue
        # match +1 / mismatch -2, matching the LTR extension. The old 2*matches - L (match +1 /
        # mismatch -1) breaks even at 50% identity, so chance matches past the true end could still
        # out-score it — a perfect 54 bp TIR scored lower than a 61 bp window at 95%. Only the 60 bp
        # scan limit hid that. At -2 the true boundary wins, and a real 80%-identity TIR still climbs.
        score = 3 * matches - 2 * L
        ident = 100.0 * matches / L
        if ident >= min_ident and score > best_score:    # extend while each added base pair keeps the score climbing
            best_L, best_score, best_ident = L, score, round(ident, 1)
    if best_L and not _is_simple(seq[:best_L]):      # (AT)n is its own reverse complement — not a TIR
        ev = {"type": "TIR (terminal inverted repeat)", "tir_len": best_L, "identity": best_ident,
              "five_prime": [0, best_L], "three_prime": [n - best_L, n], "element_span": [0, n],
              "method": "terminal inverted-repeat scan (element termini, %d bp limit)" % hi}
        if best_L == hi:                                 # the scan limit stopped it, not the repeat ending
            ev["length_is_lower_bound"] = True
            ev["bound_reason"] = ("the inverted repeat still matched at the scan limit — the length is a "
                                  "lower bound, not a measurement")
        return ev
    return None


def find_tir(seq: str, k: int = 11, min_tir: int = 10, max_tir: int = 1000, min_anchors: int = 3):
    """Detect terminal INVERTED repeats (TIR candidate, DNA transposons).
    The 1 kb limit spans the described range: ~11-30 bp for hAT and Tc1/Mariner, ~13 bp for CACTA,
    but hundreds of bp for MULE/Mutator and Foldback-type elements. A repeat still matching at the
    limit is reported as a lower bound rather than silently cut to the limit."""
    n = len(seq)
    term = _terminal_tir(seq, min_tir, max_tir)     # canonical TIR at the true termini wins
    if term:
        return term
    W = min(n // 2, 1200)
    if W < min_tir:
        return None
    a = seq[:W]
    b = reverse_complement(seq[n - W:])             # revcomp of 3' window
    idx = _kmer_index(a, k)
    inv, pts = Counter(), {}
    for j in range(len(b) - k + 1):
        for i in idx.get(b[j:j + k], ()):
            # forward 3' position of this k-mer: revcomp index j -> forward [n-W + (W-k-j)]
            p3s = (n - W) + (W - k - j)
            s = i + (p3s + k)                        # invariant ~ n for a true TIR pair
            if p3s <= i:
                continue
            inv[s] += 1
            pts.setdefault(s, []).append((i, p3s))
    if not inv:
        return None
    best_s, cnt = max(inv.items(), key=lambda kv: kv[1])
    if cnt < min_anchors:
        return None
    grp = [pt for ss, pl in pts.items() if abs(ss - best_s) <= 4 for pt in pl]
    i5 = [p for p, _ in grp]
    t5s, t5e = min(i5), max(i5) + k
    if (t5e - t5s) > max_tir or (t5e - t5s) < min_tir:
        return None
    r5 = seq[t5s:t5e]
    r3 = seq[best_s - t5e:best_s - t5s]              # 3' TIR forward region
    ident = _identity(r5, reverse_complement(r3))
    if ident < 80 or _is_simple(r5):                 # a simple repeat is not terminal-repeat evidence
        return None
    return {"type": "TIR (terminal inverted repeat)", "tir_len": t5e - t5s, "identity": ident,
            "five_prime": [t5s, t5e], "three_prime": [best_s - t5e, best_s - t5s],
            "element_span": [t5s, best_s - t5s],
            "method": "k-mer seed vs reverse-complement (k=%d)" % k}


def find_polya(seq: str, min_run: int = 8):
    """Poly-A (3') or poly-T (5') homopolymer tail — LINE/retro signature."""
    out = []
    m = len(seq)
    r = 0
    while r < m and seq[m - 1 - r] == "A":
        r += 1
    if r >= min_run:
        out.append({"type": "poly-A tail (3')", "length": r, "pos": [m - r, m]})
    r = 0
    while r < m and seq[r] == "T":
        r += 1
    if r >= min_run:
        out.append({"type": "poly-T (5')", "length": r, "pos": [0, r]})
    return out


# Expected TSD length by superfamily, restricted to values that can be attributed and to superfamilies
# classify.py can actually emit. Anything absent here is "not diagnostic" — deliberately, rather than
# filled with plausible numbers. LTR (5 bp) and hAT (8 bp): Ou et al. 2019 Genome Biol 20:275. The
# Tc1/Mariner TA dinucleotide is the defining insertion site of the superfamily (Wicker et al. 2007).
TSD_EXPECT = {
    "Copia": (5, "Ou et al. 2019 Genome Biol 20:275"),
    "Gypsy": (5, "Ou et al. 2019 Genome Biol 20:275"),
    "hAT": (8, "Ou et al. 2019 Genome Biol 20:275"),
    "Tc1/Mariner": (2, "Wicker et al. 2007 Nat Rev Genet 8:973-982 (TA target site)"),
}


def tsd_congruence(length, superfamily):
    """Compare an observed TSD length with the superfamily's expected one.

    Only INCONGRUENT carries information. An exact L-mer direct repeat at a fixed offset arises by chance
    with probability ~4^-L under a uniform null — common for short L, and commoner still in the AT-rich
    flanks TEs favour — so a matching length corroborates weakly while a clashing one is a real signal
    that the boundary, the superfamily call, or both are wrong."""
    key = (superfamily or "").split(" ")[0]
    for name, (exp, cite) in TSD_EXPECT.items():
        if key and (key == name.split("/")[0] or superfamily.startswith(name)):
            if length is None:
                return {"verdict": "not assessed", "expected": exp, "basis": cite}
            return {"verdict": "congruent" if length == exp else "incongruent",
                    "expected": exp, "observed": length, "basis": cite}
    return {"verdict": "not diagnostic", "expected": None,
            "basis": "no attributable expected TSD length for this superfamily"}


def find_tsd(seq: str, elem_start: int, elem_end: int, min_tsd: int = 2, max_tsd: int = 12, expect=None):
    """Target-site duplication: a short direct repeat immediately flanking the element.
    Requires flanking sequence beyond [elem_start, elem_end]; else returns None (honest).

    `expect` (a length, when the superfamily is known) is tested FIRST. Plain longest-first search lets a
    chance 7-mer in an AT-rich flank outrank the real 2 bp TA of a Tc1/Mariner element, which is also why
    min_tsd is 2 rather than 4 — at 4 the diagnostic TA was unreachable."""
    order = [L for L in ([expect] if expect else []) if min_tsd <= L <= max_tsd]
    order += [L for L in range(max_tsd, min_tsd - 1, -1) if L not in order]
    for L in order:
        up_s, up_e = elem_start - L, elem_start
        dn_s, dn_e = elem_end, elem_end + L
        if up_s < 0 or dn_e > len(seq):
            continue
        left, right = seq[up_s:up_e], seq[dn_s:dn_e]
        if left == right and "N" not in left:
            return {"type": "TSD (target-site duplication)", "length": L,
                    "motif": left, "upstream": [up_s, up_e], "downstream": [dn_s, dn_e],
                    "matched_expected": bool(expect and L == expect)}
    return None


# Primer-binding-site reference: PBS (viral +strand, 5'->3') = reverse complement of the 3'-terminal 18 nt
# of the priming tRNA. Endogenised proviruses carry a DIVERGED PBS, so a match is often partial — the tRNA
# identity is reported hedged, never as a hard call. Panel is Lys-anchored (HML-2 primes tRNA-Lys); the Lys3
# entry equals the canonical HIV-1 tRNA-Lys3 PBS (independently verifiable) and anchors the panel's accuracy.
_PRIMER_TRNA = {
    "tRNA-Lys3": "TGGCGCCCGAACAGGGAC",     # Lys-primed lentiviruses / HML-2-class ERVs
}


def find_pbs(seq: str, ltr_five_prime_end: int, search: int = 44, min_ident: float = 0.55):
    """Primer-binding site: the ~18 nt just 3' of the 5' LTR that the priming tRNA anneals to. Detected
    by best reverse-complement match to a bundled primer-tRNA panel, at the canonical leader position.
    tRNA identity is hedged when the match is weak (expected for endogenised, diverged proviruses)."""
    leader = seq[ltr_five_prime_end:ltr_five_prime_end + search]
    best = None
    for name, pbs in _PRIMER_TRNA.items():
        L = len(pbs)
        for i in range(len(leader) - L + 1):
            w = leader[i:i + L]
            ident = sum(1 for a, b in zip(w, pbs) if a == b and a not in "Nn") / L
            if best is None or ident > best["_id"]:
                best = {"_id": ident, "trna": name, "pos": [ltr_five_prime_end + i, ltr_five_prime_end + i + L], "motif": w}
    if not best or best["_id"] < min_ident:            # no credible PBS in the leader window
        return None
    # the panel is Lys-anchored, so only NAME the priming tRNA on a strong match; a weak (diverged) match
    # is reported as undetermined with the closest panel match, never a hard call for a non-Lys genus.
    strong = best["_id"] >= 0.72
    return {"type": "PBS (primer-binding site)", "pos": best["pos"],
            "priming_trna": best["trna"] if strong else "undetermined",
            "best_match": best["trna"], "identity": round(100 * best["_id"], 1),
            "confident": strong, "motif": best["motif"],
            "note": "" if strong else (f"priming tRNA undetermined — closest panel match {best['trna']} "
                                       f"({round(100 * best['_id'], 1)}%), below the confident threshold; "
                                       "endogenised PBS is often diverged")}


def find_ppt(seq: str, ltr_three_prime_start: int, window: int = 30, min_len: int = 9,
             min_purine: float = 0.82, max_defects: int = 2):
    """Polypurine tract: the run of purines (A/G) abutting the 3' LTR that primes plus-strand synthesis.
    Extended backward from the LTR boundary while it stays a dense purine run (bounded pyrimidine defects,
    no N) and TRIMMED so the reported tract starts on a purine — never a window with leading pyrimidines."""
    reg_s = max(0, ltr_three_prime_start - window)
    region = seq[reg_s:ltr_three_prime_start]
    best_start, defects = None, 0
    for i in range(len(region) - 1, -1, -1):           # extend 5' from the 3'-LTR boundary
        c = region[i]
        if c in "Nn":
            break
        if c not in "AG":
            defects += 1
            if defects > max_defects:
                break
        sub = region[i:]                               # candidate tract must START on two purines and stay dense
        starts_clean = c in "AG" and (i + 1 >= len(region) or region[i + 1] in "AG")
        if starts_clean and len(sub) >= min_len and sum(1 for x in sub if x in "AG") / len(sub) >= min_purine:
            best_start = i
    if best_start is None:
        return None
    sub = region[best_start:]
    return {"type": "PPT (polypurine tract)", "pos": [reg_s + best_start, ltr_three_prime_start],
            "length": len(sub), "purine_frac": round(sum(1 for x in sub if x in "AG") / len(sub), 2), "motif": sub}


def _tsd_state(seq: str, start: int, end: int, found, min_tsd: int = 2) -> str:
    """Three outcomes, not two. A TSD can only be looked for when the record extends beyond the element
    on BOTH sides; a bare pasted element — the normal case — has no flanks, so 'no TSD' there means the
    question was never asked. Reporting that as 'absent' would be evidence the run does not have."""
    if found:
        return "found"
    if start < min_tsd or end + min_tsd > len(seq):
        return "not assessable (no flanking sequence in the record)"
    return "absent (flanks present, no direct repeat found)"


def detect_all(seq: str):
    """Run all structural detectors. Returns list of evidence dicts (empty if none)."""
    ev = []
    ltr = find_ltr(seq)
    if ltr:
        ev.append(ltr)
        tsd = find_tsd(seq, ltr["element_span"][0], ltr["element_span"][1])
        ltr["tsd_state"] = _tsd_state(seq, ltr["element_span"][0], ltr["element_span"][1], tsd)
        if tsd:
            ev.append(tsd)
        pbs = find_pbs(seq, ltr["five_prime"][1])       # LTR-class cis-elements: PBS (leader) + PPT (before 3' LTR)
        if pbs:
            ev.append(pbs)
        ppt = find_ppt(seq, ltr["three_prime"][0])
        if ppt:
            ev.append(ppt)
    # LTR (direct) and TIR (inverted) terminal architectures are mutually exclusive:
    # only look for a TIR when no LTR was found, so an LTR element never reports a spurious TIR.
    tir = find_tir(seq) if not ltr else None
    if tir:
        ev.append(tir)
        tsd = find_tsd(seq, tir["five_prime"][0], tir["three_prime"][1])
        tir["tsd_state"] = _tsd_state(seq, tir["five_prime"][0], tir["three_prime"][1], tsd)
        if tsd:
            ev.append(tsd)
    ev.extend(find_polya(seq))
    return ev
