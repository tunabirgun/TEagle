"""Real structural TE evidence (Layer D): terminal direct repeats (LTR), terminal
inverted repeats (TIR), poly-A tails, flanking target-site duplications (TSD).
Heuristic detectors — reported as candidate structural evidence, never as family/DB calls."""
from __future__ import annotations
import math
from bisect import bisect_left, bisect_right
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
    return len(quads) < max(4, min(192, 0.05 * (len(u) - 3)))     # 192 = 0.75 x the 256 possible 4-mers: longer repeats stay judgeable


# A seed word seen more than this often in one search window is a microsatellite, a homopolymer or a repeat
# family, not terminal-repeat evidence, and is masked (the word masking of BLAST/DUST). Unmasked, its occurrences
# pair with each other quadratically (a 3-kb (AT)n window gives millions of seed pairs, 45-90 s in find_tir).
# A random 11-mer is seen about once in 4 Mb, a random 6-mer about once in a 400-bp window, and a word that a
# real terminal repeat carries is seen in each window about once, so 20 bounds the pairs at 20 x window without
# touching ordinary sequence.
_SEED_MAX_OCC = 20
# Diverged microsatellites escape a count cap (every mutated copy of the word is rare) and then chain into a
# gapped alignment of up to _MAX_BAND_CELLS, so the LTR seed words are also masked where the sequence is
# low-complexity by the DUST score (Morgulis et al. 2006 J Comput Biol 13:1028): in a 64-bp window, sum over
# triplets of c(c-1)/2, divided by (triplets - 1), above 2 (dustmasker's default level 20, window 64; random
# sequence scores about 0.5). The TIR seed search uses the count cap only: DUST also masks 9-20% of the Ac and
# Tc1 specimens (their terminal and subterminal repeats), which would change real calls.
_DUST_WINDOW, _DUST_SCORE = 64, 2.0


def _dust(s: str) -> bytearray:
    """1 at every base of a window whose DUST score exceeds _DUST_SCORE, else 0."""
    n = len(s)
    bad = bytearray(n)
    if n < 4:
        return bad
    w = min(_DUST_WINDOW, n)
    limit = _DUST_SCORE * (w - 3)
    counts, score = {}, 0
    for i in range(n - 2):
        t = s[i:i + 3]
        c = counts.get(t, 0)
        score += c
        counts[t] = c + 1
        if i >= w - 2:
            r = s[i + 2 - w:i + 5 - w]
            c = counts[r] - 1
            counts[r] = c
            score -= c
        if i >= w - 3 and score > limit:
            bad[i + 3 - w:i + 3] = b"\x01" * w
    return bad


def _kmer_index(a: str, k: int, max_occ: int | None = None, mask: bool = False):
    idx = {}
    clean = "N" not in a and "n" not in a
    bad = _dust(a) if mask else b""
    if 1 not in bad:
        bad = b""
    for i in range(len(a) - k + 1):
        kmer = a[i:i + k]
        if not clean and ("N" in kmer or "n" in kmer):   # an N-run is not repeat evidence; never seed anchors on it
            continue
        if bad and bad.find(1, i, i + k) >= 0:           # a word on low-complexity sequence seeds nothing
            continue
        hit = idx.get(kmer)
        if hit is None:
            idx[kmer] = [i]
        else:
            hit.append(i)
    if max_occ:
        idx = {w: p for w, p in idx.items() if len(p) <= max_occ}
    return idx


_ACGT = frozenset("ACGT")


def _identity_stats(r5: str, r3: str):
    """(matches, known columns, ambiguous columns) of an ungapped pair. A column counts only when BOTH
    bases are unambiguous: N and the IUPAC codes are unknown, so they are neither a match nor a mismatch."""
    m = known = 0
    for x, y in zip(r5, r3):
        if x in _ACGT and y in _ACGT:
            known += 1
            m += x == y
    return m, known, min(len(r5), len(r3)) - known


def _identity(r5: str, r3: str) -> float:
    m, known, _ = _identity_stats(r5, r3)
    return round(100 * m / known, 1) if known else 0.0


# Terminal-motif panel. A motif is a 4-mer: the LTR's first two bases + its last two ("TGCA" = 5'-TG…CA-3').
# TGCA is the canonical integrase att pair; the other seven are the non-canonical motifs LTR_retriever
# searches when TGCA is absent — verified against its manual (Ou & Jiang 2018 Plant Physiol 176:1410;
# LTR_retriever 2.9.0 manual Table 1: -motif [TCCA TGCT TACA TACT TGGA TATA TGTA TGCA]). The list is
# copied from that default, not derived: it is prior knowledge about which termini occur in nature.
_CANONICAL_MOTIF = "TGCA"
_NONCANONICAL_MOTIFS = ("TCCA", "TGCT", "TACA", "TACT", "TGGA", "TATA", "TGTA")


def _motif_mismatches(obs: dict, motif: str) -> int:
    want = {"five_start": motif[:2], "five_end": motif[2:],
            "three_start": motif[:2], "three_end": motif[2:]}
    return sum(1 for k, v in want.items() if obs[k] != v)


def _termini_motif(seq: str, l5s: int, l5e: int, l3s: int, l3e: int) -> dict:
    """Report whether the LTR termini carry the canonical 5'-TG…CA-3' dinucleotides, and — only when they
    do not — whether they match one of the known non-canonical terminal motifs instead.

    A BADGE, never a boundary mover. The motif is characteristic of many LTR retrotransposons but is not
    universal — measured here, copia (M11240) satisfies it exactly while gypsy/mdg4 (M12927) reads AG…TT
    and no offset within +/-10 bp satisfies it on both copies. Snapping boundaries onto the motif would
    therefore have MOVED a correct gypsy boundary to satisfy a rule that element does not follow.
    Tolerance of one mismatch across the four positions follows LTRharvest's -motifmis 1
    (Ellinghaus et al. 2008 BMC Bioinformatics 9:18). The non-canonical panel is only consulted after the
    canonical test fails, and only an EXACT match counts — a mismatch-tolerant search across eight motifs
    would match almost any terminus by chance and turn the badge into noise."""
    obs = {"five_start": seq[l5s:l5s + 2].upper(), "five_end": seq[l5e - 2:l5e].upper(),
           "three_start": seq[l3s:l3s + 2].upper(), "three_end": seq[l3e - 2:l3e].upper()}
    mism = _motif_mismatches(obs, _CANONICAL_MOTIF)
    alt = next((m for m in _NONCANONICAL_MOTIFS if _motif_mismatches(obs, m) == 0), None) if mism else None
    if mism == 0:
        note = "termini match the canonical TG…CA"
    elif alt:
        note = (f"termini do not match TG…CA but match the non-canonical motif {alt[:2]}…{alt[2:]}, one of "
                "the seven non-TGCA termini LTR_retriever searches (Ou & Jiang 2018 Plant Physiol 176:1410); "
                "the motif is corroborating context only and does not set the boundaries")
    else:
        note = ("termini match neither TG…CA nor any of the seven known non-canonical motifs; the motif is "
                "common but not universal among LTR retrotransposons, so its absence is not evidence "
                "against the call")
    return {**obs, "canonical": mism == 0, "mismatches": mism,
            "within_ltrharvest_tolerance": mism <= 1,
            "noncanonical_motif": alt,                    # None unless an exact non-TGCA panel match
            "motif_tier": "canonical" if mism == 0 else ("non-canonical" if alt else "none"),
            "note": note}


def _extend(seq: str, i: int, j: int, step: int, limit: int, x_drop: int = 30) -> int:
    """Ungapped X-drop extension of a repeat pair along a fixed diagonal (BLAST-style: Altschul 1990).
    Walks (i, j) outward by `step` (+1 right, -1 left) and returns the offset of the best-scoring
    prefix — how far the pair genuinely stays a repeat, not how far it can be walked. Match +1 /
    mismatch -1 separates an 80%-identity LTR pair (+0.6 per base) from random (-0.5 per base)."""
    score = best = best_off = 0
    for off in range(max(0, limit)):
        a, b = seq[i + step * off], seq[j + step * off]
        if a not in _ACGT or b not in _ACGT:         # unknown base: neither evidence for nor against the pair
            continue
        if a == b:
            score += 1
        else:
            score -= 2                               # break-even at 67% identity: a diverged pair still
        if score > best:                             # extends, a chance match past the end cannot set a best
            best, best_off = score, off + 1
        elif score < best - x_drop:
            break
    return best_off


# Gapped measurement of a terminal direct repeat. LTR copies diverge by indels as well as substitutions, and
# an ungapped comparison puts everything downstream of a 1-bp indel out of phase (a 99.6%-identical copia LTR
# pair reads 67%). Scoring matches the extension above (+1 / -2); a gap of g columns costs _GAP_OPEN + g.
_GAP_OPEN, _GAP_EXT = 4, 1
_LTR_MAX_INDEL = 32            # anchors up to this many diagonals from the best one may be chained across an indel
_GAP_BAND_MARGIN = 8           # slack either side of the diagonals the chain visits
_GAPPED_BELOW_IDENTITY = 97.0  # the gapped alignment is only attempted when the ungapped pair is below this
_NEG = -10 ** 9
_MAX_BAND_CELLS = 250_000      # a gapped alignment bigger than this is not attempted: the ungapped measurement stands


def _gapped_extend(seq: str, i0: int, j0: int, step: int, lim_a: int, lim_b: int, w: int, x_drop: int = 30):
    """Banded affine-gap X-drop extension of a repeat pair. Consumes a bases of copy A (from i0, direction
    `step`, at most lim_a) and b bases of copy B (from j0, at most lim_b) with |a - b| <= w; returns
    (best score, a, b) of the best-scoring prefix pair."""
    oe = _GAP_OPEN + _GAP_EXT
    best = bi = bj = 0
    Hp = {0: 0}
    for j in range(1, min(w, lim_b) + 1):
        Hp[j] = -(_GAP_OPEN + j * _GAP_EXT)
    Fp = {}
    for r in range(1, lim_a + 1):
        a = seq[i0 + step * (r - 1)]
        H, F = {}, {}
        e = _NEG
        rowmax = _NEG
        for j in range(max(0, r - w), min(lim_b, r + w) + 1):
            h = _NEG
            if j >= 1 and (j - 1) in Hp:
                b = seq[j0 + step * (j - 1)]
                h = Hp[j - 1] + (0 if a not in _ACGT or b not in _ACGT else (1 if a == b else -2))
            f = max(Hp[j] - oe, Fp.get(j, _NEG) - _GAP_EXT) if j in Hp else _NEG
            e = max(H[j - 1] - oe, e - _GAP_EXT) if (j - 1) in H else _NEG
            v = max(h, f, e)
            H[j], F[j] = v, f
            if v > rowmax:
                rowmax = v
            if v > best:
                best, bi, bj = v, r, j
        if rowmax < best - x_drop:
            break
        Hp, Fp = H, F
    return best, bi, bj


def _band_align(a: str, b: str, lo: int, hi: int):
    """Global affine-gap alignment of a and b restricted to the band lo <= j - i <= hi.
    Returns (matches, mismatches, gap columns, gap events, ambiguous columns) or None when the band
    cannot reach the far corner or the matrix would be unreasonably large."""
    la, lb = len(a), len(b)
    lo, hi = min(lo, 0, lb - la), max(hi, 0, lb - la)
    if la * (hi - lo + 1) > _MAX_BAND_CELLS:
        return None
    oe = _GAP_OPEN + _GAP_EXT
    tabs = ([], [], [])                      # M (diagonal), E (consumes b), F (consumes a): per row {j: (score, back-state)}
    M, E, F = tabs

    def sc(t, i, j):
        v = t[i].get(j)
        return v[0] if v else _NEG

    def pick(*cands):
        best = max(cands)
        return best if best[0] > _NEG // 2 else None

    for i in range(la + 1):
        for t in tabs:
            t.append({})
        for j in range(max(0, i + lo), min(lb, i + hi) + 1):
            if i == 0 and j == 0:
                M[0][0] = (0, 3)
                continue
            if i and j:
                x, y = a[i - 1], b[j - 1]
                s = 0 if x not in _ACGT or y not in _ACGT else (1 if x == y else -2)
                v = pick(*((sc(t, i - 1, j - 1) + s, st) for st, t in enumerate(tabs)))
                if v:
                    M[i][j] = v
            if j:
                v = pick((sc(M, i, j - 1) - oe, 0), (sc(E, i, j - 1) - _GAP_EXT, 1), (sc(F, i, j - 1) - oe, 2))
                if v:
                    E[i][j] = v
            if i:
                v = pick((sc(M, i - 1, j) - oe, 0), (sc(E, i - 1, j) - oe, 1), (sc(F, i - 1, j) - _GAP_EXT, 2))
                if v:
                    F[i][j] = v
    top, state = max((sc(t, la, lb), st) for st, t in enumerate(tabs))
    if top <= _NEG // 2:
        return None
    i, j = la, lb
    mt = mm = gc = ge = amb = 0
    while i or j:
        back = tabs[state][i][j][1]
        if state == 0:
            x, y = a[i - 1], b[j - 1]
            if x not in _ACGT or y not in _ACGT:
                amb += 1
            elif x == y:
                mt += 1
            else:
                mm += 1
            i, j = i - 1, j - 1
        else:
            gc += 1
            ge += back != state
            if state == 1:
                j -= 1
            else:
                i -= 1
        state = back
    return mt, mm, gc, ge, amb


# Seeding window for the terminal-repeat scan. It must EXCEED the repeat being sought: the two windows
# only share k-mers where they overlap in repeat-internal coordinates, so a window at or below the true
# LTR length seeds nothing. 6 kb covers the longest LTRs described (multi-kb plant Ty3/Ogre LTRs).
_LTR_SEED_WINDOW = 6000


MIN_LTR_IDENTITY = 80.0      # named so the panels that describe this detector can state the real floor
# The TIR floor is the LTR floor's twin and gates every DNA-transposon call the same way, so it is named
# for the same reason: a 79%-identity terminal inverted repeat is silently rejected by a number the user
# cannot see on screen and a reviewer cannot recover from the sealed manifest.
MIN_TIR_IDENTITY = 80

# The length a terminal inverted repeat must reach before it can carry a Class II call ON ITS OWN, with no
# transposase to corroborate it. MIN_TIR_IDENTITY above governs whether a repeat is detected at all; this
# governs whether a detected repeat is strong enough to name an element by itself.
#
# Derived, not chosen. The check that runs with the suite is live: tests/test_41_structural.py::live_chance_tir
# runs find_tir over seeded random sequences (seed 41) on every run, so a find_tir change that lifts chance
# repeats toward this floor fails the suite. The larger measurement is a committed script rather than this comment: run
# `python benchmarks/chance_tir.py`, which runs find_tir over 4,500 random sequences (1.5, 4 and 9 kb,
# 1,500 each) at a recorded seed and writes benchmarks/raw/chance_tir.json. The chance hit rate is about
# 1%, and the length distribution falls away steeply — the great majority of chance hits are 13 bp.
#
# The ceiling is itself a sampled extreme, so it moves with the seed: the recorded seed gives 15 bp, and
# the script's own seed-sensitivity pass over five seeds gives a maximum of 17 bp. Twenty is set above the
# across-seed maximum rather than above any single run, which is the conservative reading, and still sits
# well below the 31 bp repeats of the real non-autonomous elements in the benchmark corpus. Chance repeat
# length grows with template length, so this ceiling holds for the lengths sampled and not in general. 4.1
# re-ran the same measurement after find_tir gained its flanked-record search: ceilings 15, 17, 16, 17, 15 bp
# over the same five seeds, the same 15-17 bp range as before (an uncorroborated candidate is reported as its
# exact core, a corroborated one is measured at _TIR_WEAK_IDENTITY, precisely so that this stays true).
#
# The floor also has to survive the case that motivated the conservative branch it guards: a solo copia 5'
# LTR (276 bp) fed alone trips the terminal-inverted-repeat scan and used to be filed as a DNA transposon.
# Measured, that solo LTR yields a 12 bp repeat, so this floor rejects it. tests/test_structural_bounds.py
# re-derives both numbers, so a change to find_tir that moves either one fails the suite rather than
# silently loosening the gate.
MIN_TIR_STANDALONE = 20


def _ltr_stats(l5s, l5e, l3s, l3e, mode, mt, mism, amb, events=0, gap_cols=0):
    cols = mt + mism + gap_cols
    return {"l5s": l5s, "l5e": l5e, "l3s": l3s, "l3e": l3e, "mode": mode, "amb": amb, "events": events,
            "ident": round(100 * mt / cols, 1) if cols else 0.0,
            "score": mt - 2 * mism - _GAP_OPEN * events - _GAP_EXT * gap_cols}


def _ltr_ungapped(seq, grp, best_d, k, min_ltr, min_anchors):
    """Seeded diagonal pair measured by ungapped X-drop extension. Returns a candidate dict or None."""
    n = len(seq)
    p5s = [p for p, _ in grp]
    l5s, l5e = min(p5s), max(p5s) + k
    L = l5e - l5s
    if L < min_ltr:
        return None
    # Anchor density floor. This only has to reject a handful of chance k-mers sharing a diagonal; the
    # real guard is the identity of the EXTENDED copies. The old 40% floor demanded that 40% of
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
    mt, known, amb = _identity_stats(seq[l5s:l5e], seq[l3s:l3e])
    c = _ltr_stats(l5s, l5e, l3s, l3e, "ungapped", mt, known - mt, amb)
    # A repeat that reaches a record edge, or that was stopped only by running into its own other copy,
    # may continue past what was supplied: the length is then a LOWER BOUND and must never be presented
    # as a measurement. Touching the edge counts even when no extension was needed to get there.
    c["edge"] = l5s == 0 or l3e == n
    c["meet"] = (ext_l > 0 and ext_l == lim_l and lim_l == gap) or (ext_r > 0 and ext_r == lim_r and lim_r == gap)
    return c


def _chain(anchors):
    """Longest collinear subset of seed anchors (both coordinates strictly increasing)."""
    anchors = sorted(anchors, key=lambda t: (t[0], -t[1]))
    tails, idxs, prev = [], [], [-1] * len(anchors)
    for x, (_, p3) in enumerate(anchors):
        pos = bisect_left(tails, p3)
        if pos == len(tails):
            tails.append(p3)
            idxs.append(x)
        else:
            tails[pos], idxs[pos] = p3, x
        prev[x] = idxs[pos - 1] if pos else -1
    out, x = [], (idxs[-1] if idxs else -1)
    while x >= 0:
        out.append(anchors[x])
        x = prev[x]
    return out[::-1]


def _ltr_gapped(seq, pts, best_d, k, min_ltr, min_anchors):
    """Seeded pair measured by a gapped alignment: anchors on neighbouring diagonals are chained across an
    indel, the ends are extended with a banded affine-gap X-drop, and identity is matches over aligned
    columns of one banded global alignment. Returns a candidate dict or None."""
    n = len(seq)
    chain = _chain([pt for d, pl in pts.items() if abs(d - best_d) <= _LTR_MAX_INDEL for pt in pl])
    if len(chain) < min_anchors:
        return None
    a5, a3 = chain[0]
    b5, b3 = chain[-1][0] + k, chain[-1][1] + k
    L = b5 - a5
    if L < min_ltr or b5 > a3 or len(chain) < max(min_anchors, 0.01 * (L - k + 1)):
        return None
    half = (a3 - b5) // 2                            # internal region: neither copy may run into the other
    w = _GAP_BAND_MARGIN + 4
    _, ial, jbl = _gapped_extend(seq, a5 - 1, a3 - 1, -1, a5, half, w)
    _, iar, jbr = _gapped_extend(seq, b5, b3, +1, half, n - b3, w)
    l5s, l3s, l5e, l3e = a5 - ial, a3 - jbl, b5 + iar, b3 + jbr
    rel = [q - p - (l3s - l5s) for p, q in chain]
    res = _band_align(seq[l5s:l5e], seq[l3s:l3e], min(rel) - _GAP_BAND_MARGIN, max(rel) + _GAP_BAND_MARGIN)
    if res is None:
        return None
    mt, mm, gc, ge, amb = res
    c = _ltr_stats(l5s, l5e, l3s, l3e, "gapped", mt, mm, amb, ge, gc)
    c["edge"] = l5s == 0 or l3e == n
    c["meet"] = (jbl > 0 and jbl == half) or (iar > 0 and iar == half)
    return c


def find_ltr(seq: str, k: int = 13, min_ltr: int = 80, min_anchors: int = 4, near_miss: list | None = None):
    """Detect a terminal DIRECT repeat pair (LTR candidate). Coords 0-based half-open.
    Anchors only SEED the repeat — its extent is measured by X-drop extension along the seeded
    diagonal over the whole sequence, so the reported length is not an artefact of the window. A pair
    that is not near-identical under the ungapped comparison is re-measured with a gapped alignment, so an
    indel does not read as divergence; the better-scoring measurement is reported."""
    seq = seq.upper()
    n = len(seq)
    W = min(n // 2, _LTR_SEED_WINDOW)
    if W < min_ltr:
        return None
    a, b = seq[:W], seq[n - W:]                     # b offset by (n-W) in full seq
    idx = _kmer_index(a, k, _SEED_MAX_OCC, True)
    idx_b = _kmer_index(b, k, _SEED_MAX_OCC, True)         # a word masked in either window seeds nothing
    diag, pts = Counter(), {}
    for j in range(len(b) - k + 1):
        w = b[j:j + k]
        if w not in idx_b:
            continue
        for i in idx.get(w, ()):
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
    cand = _ltr_ungapped(seq, grp, best_d, k, min_ltr, min_anchors) if len(grp) >= min_anchors else None
    if (cand is None or cand["ident"] < _GAPPED_BELOW_IDENTITY
            or any(3 < abs(d - best_d) <= _LTR_MAX_INDEL for d in pts)):
        g = _ltr_gapped(seq, pts, best_d, k, min_ltr, min_anchors)
        if g and (cand is None or g["score"] > cand["score"]):
            cand = g
    if cand is None:
        return None
    l5s, l5e, l3s, l3e = cand["l5s"], cand["l5e"], cand["l3s"], cand["l3e"]
    L, ident = l5e - l5s, cand["ident"]
    if ident < MIN_LTR_IDENTITY or _is_simple(seq[l5s:l5e]):   # identity carries the weight the density floor used to
        # A candidate pair that fails ONLY the identity floor is not the same as no candidate at all, and the
        # difference matters scientifically: LTR-LTR identity decays with insertion age, so a hard floor makes
        # the detector silently blind to older elements and reports that blindness as absence. The measured
        # identity is handed back as an advisory near-miss (same idiom as the PBS 55-72% reported-but-not-named
        # band) so "below threshold" is visibly different from "nothing there". Simple/low-complexity repeats
        # are NOT reported: those are a true negative, not a diverged element.
        if near_miss is not None and not _is_simple(seq[l5s:l5e]):
            near_miss.append({
                # The type MUST NOT begin with "LTR". Four consumers dispatch on that prefix
                # (classify.has_ltr, retroviral.py's ltr lookup, figures.py's band colour, main.py's
                # evidence roll-up), so an advisory row named "LTR ..." was read as a CONFIRMED terminal
                # repeat: a pair rejected for being below the identity floor came back as LTR/Copia at HIGH
                # confidence — strictly worse than the silent discard this row was added to replace.
                "type": "Sub-threshold terminal direct repeat (advisory)",
                "ltr_len": L, "identity": ident, "threshold": MIN_LTR_IDENTITY,
                "five_prime": [l5s, l5e], "three_prime": [l3s, l3e], "element_span": [l5s, l3e],
                "advisory": True,
                "method": ("k-mer seed (k=%d, %d bp window) + X-drop diagonal extension; reported because the "
                           "pair fell below the %.0f%% identity floor, NOT accepted as an LTR" % (k, W, MIN_LTR_IDENTITY)),
            })
        return None
    ev = {"type": "LTR (terminal direct repeat)", "ltr_len": L, "identity": ident,
          "five_prime": [l5s, l5e], "three_prime": [l3s, l3e],
          "element_span": [l5s, l3e],
          "method": "k-mer seed (k=%d, %d bp window) + %s" % (
              k, W, "banded gapped alignment" if cand["mode"] == "gapped" else "X-drop diagonal extension")}
    if cand["events"]:
        ev["indels"] = cand["events"]
    if cand["amb"]:
        ev["ambiguous_columns"] = cand["amb"]
    ev["termini"] = _termini_motif(seq, l5s, l5e, l3s, l3e)
    causes = [c for c, f in (("the repeat still matched where the record ends", cand["edge"]),
                             ("the two copies met each other", cand["meet"])) if f]
    if causes:
        ev["length_is_lower_bound"] = True
        ev["bound_reason"] = (" and ".join(causes) + ", so the copies may extend beyond it — the length is a "
                              "lower bound, not a measurement")
    return _snap_ltr(seq, ev)


# Where a terminal inverted repeat is looked for. The anchored scan covers a record that starts and ends on
# the element. A flanked record (any genomic window) needs the arms located, and the search space grows with
# the square of the flank, so two seed searches run, each bounded and named so the bound is on record:
#   _TIR_SEED_WINDOW   exact 11-mer seeds (>= 3 on one diagonal, i.e. a ~13 bp perfect core) within this many
#                      bp of each record end. Strong enough to stand on its own: chance gives ~1% (chance_tir).
#   _TIR_JOINT_WINDOW  exact 6-mer seeds (two on one diagonal) within this many bp. Short or imperfect TIRs (hAT's 11 bp at 91%) have
#                      no 11-mer, and at k=6 the arms are indistinguishable from chance on their own, so such a
#                      candidate is accepted only when a flanking TSD of >= _TIR_WEAK_TSD bp corroborates it.
# A record flanked by more than the window on either side is not searched beyond it; scan_scope() says so.
_TIR_SEED_WINDOW = 1200
_TIR_JOINT_WINDOW = 400
_TIR_JOINT_K = 6
_TIR_JOINT_ANCHORS = 2          # 6-mer seeds needed on one diagonal: an exact run of 7+, or two shorter ones
MIN_TIR_LEN = 10               # find_tir's default min_tir; arms are never trimmed below it
_TIR_XDROP = 60                # score below the best at which a seeded candidate's scan stops
_TIR_MAX_SEEDS = 64             # candidate diagonals kept per seed search; each costs a scan, and ordinary sequence has a handful
_TIR_JOINT_MAX = 60             # longest arm a 6-mer-seeded candidate is measured to; longer arms carry 11-mer seeds
_TIR_WEAK_TSD = 6
_TIR_WEAK_IDENTITY = 90.0        # identity the arms of a TSD-corroborated candidate are measured at
_TIR_OUT_MAX = 2                # a terminal pair that does not pair is cut off by the seeded extension; the TSD can restore up to this many
_TIR_TRIM_MAX = 6               # arms found by extension may have run through the target site; trim at most this many bp


def scan_scope(seq: str) -> dict:
    """What the terminal-repeat scans can see of a record, for a caller that must say so when nothing is found."""
    n = len(seq)
    ltr, tir, joint = min(n // 2, _LTR_SEED_WINDOW), min(n // 2, _TIR_SEED_WINDOW), min(n // 2, _TIR_JOINT_WINDOW)
    tail = min(n, _TAIL_WINDOW)
    return {"record_bp": n, "ltr_seed_window_bp": ltr, "tir_seed_window_bp": tir, "tir_joint_window_bp": joint,
            "tail_window_bp": tail, "covers_whole_record": 2 * max(ltr, tir) >= n and tail >= n}


# GC range over which the random-sequence TIR ceiling quoted with a call was checked: the benchmark runs at
# 50% GC and the seeded tests at 35% and 65%. Chance inverted repeats lengthen as composition departs from
# uniform, so the ceiling is not a statement about a record outside this range.
CHANCE_TIR_GC_RANGE = (0.35, 0.65)


def chance_tir(seq: str, length: int):
    """Chance that random sequence of this record's composition holds a perfect inverted repeat of `length` bp
    between the two outer seed windows, or None without ACGT bases. Per placement a base pairs with probability
    q = 2(fA fT + fC fG); the two windows give W x W placements, so p = 1 - exp(-W^2 q^length). Exact arms only:
    an imperfect arm of the same length is likelier by chance, so this understates it. It bounds the live
    detector from above for arms of 14 bp and longer (tests/test_41_structural.py)."""
    s = seq.upper()
    f = {b: s.count(b) for b in "ACGT"}
    tot = sum(f.values())
    if not tot:
        return None
    q = 2 * (f["A"] * f["T"] + f["C"] * f["G"]) / tot ** 2
    w = scan_scope(seq)["tir_seed_window_bp"]
    gc = (f["C"] + f["G"]) / tot
    return {"gc": round(gc, 3), "placements": w * w, "per_placement": float(f"{q ** length:.3g}"),
            "chance_probability": float(f"{-math.expm1(-w * w * q ** length):.3g}"),
            "in_measured_gc_range": CHANCE_TIR_GC_RANGE[0] <= gc <= CHANCE_TIR_GC_RANGE[1]}


def _scan_tir(seq: str, s: int, e: int, min_tir: int, max_tir: int, min_ident: float = MIN_TIR_IDENTITY,
              x_drop: int | None = None, reach: int = 16):
    """Inverted-repeat arms of an element whose outer ends are [s, e): returns (length, identity, matches,
    scan limit) of the best window, or None. A real TIR is often imperfect, so the length is the LONGEST
    window whose identity stays above threshold before the identity 'cliff' where the match becomes random —
    not the shortest, most perfect core. Catches short canonical TIRs (11 bp Ac) and longer imperfect ones
    (28 bp mariner)."""
    hi = min(max_tir, (e - s) // 2)
    if hi < min_tir:
        return None
    b_all = reverse_complement(seq[e - hi:e])            # revcomp of the 3' window; its first L bases = revcomp(seq[e-L:e])
    # `matches` accumulates across L (b_all[:L] is a true prefix for every L), so the whole scan is one
    # pass instead of a re-count per window — which is what makes a realistic max_tir affordable.
    best = None
    best_score = 0
    matches = known = 0
    for L in range(1, hi + 1):
        x, y = seq[s + L - 1], b_all[L - 1]
        if x in _ACGT and y in _ACGT:                    # an ambiguity code is not base-pairing evidence either way
            known += 1
            matches += x == y
        if L < min_tir or not known:
            continue
        # match +1 / mismatch -2, matching the LTR extension. The old 2*matches - L (match +1 /
        # mismatch -1) breaks even at 50% identity, so chance matches past the true end could still
        # out-score it — a perfect 54 bp TIR scored lower than a 61 bp window at 95%. Only the 60 bp
        # scan limit hid that. At -2 the true boundary wins, and a real 80%-identity TIR still climbs.
        score = 3 * matches - 2 * known
        ident = 100.0 * matches / known
        if ident >= min_ident and score > best_score:    # extend while each added base pair keeps the score climbing
            best, best_score = (L, round(ident, 1), matches, hi), score
        elif x_drop and best and score < best_score - x_drop:
            break                                        # past the cliff: a seeded candidate need not be walked to the limit
        elif x_drop and not best and L >= 2 * min_tir and ident < min_ident - 10:
            break                                        # never reached the identity floor and now far below it
    if best and not _simple_arms(seq, s, e, best[0], reach):    # (AT)n is its own reverse complement — not a TIR
        return best
    return None


def _period_run(w: str) -> bool:
    """w is a tandem run of period 1-2, (AT)n or A(n), or from 20 bp of period 3, (CAG)n, or of period 4-6 with a
    unit that is its own reverse complement, (GATC)n, (AAATTT)n, or repeats a shorter unit, with at most 15% of its
    comparisons broken (30% from 20 bp). Such a repeat satisfies an inverted-repeat identity test against itself (a
    palindromic unit pairs with itself outright), so it is no terminal-repeat evidence. A longer unit that is neither,
    (TCTA)n, cannot pair with itself, and an arm made of it is left to the identity test."""
    tol = 0.15 if len(w) < 20 else 0.3
    return any(sum(w[i] != w[i + p] for i in range(len(w) - p)) <= tol * (len(w) - p) and (p < 4 or _short_or_paired(w, p))
               for p in range(1, 7 if len(w) >= 20 else 3))


def _simple_arms(seq: str, s: int, e: int, L: int, reach: int = 16) -> bool:
    """The arms seq[s:s+L] and seq[e-L:e] are a microsatellite, alone or with the `reach` bp outside them or the same
    length inside them. A tandem repeat pairs with its own reverse complement, so it is no inverted-repeat evidence,
    and a few substitutions do not make it one (a 35-bp (TA)n arm with one C, or a 10-bp CTGCAGCAGC, is not
    _is_simple). The inside context is what exposes a short arm at the edge of a record cut from an array, which has
    no outside; `reach` is 0 for arms from which a run outside has been cut."""
    ws = (seq[s:s + L], seq[e - L:e], seq[max(0, s - reach):s + L], seq[e - L:e + reach],
          seq[s:s + L + 16], seq[max(0, e - L - 16):e])
    return _is_simple(ws[0]) or any(_period_run(w) for w in ws)


def _tsd_in_microsatellite(seq: str, s: int, e: int, t: int, reach: int = 8) -> bool:
    """The t-bp TSD flanking [s, e) lies in a tandem run that goes on for `reach` more bases outside it."""
    return any(len(w) >= t + reach and _period_run(w) for w in (seq[max(0, s - t - reach):s], seq[e:e + t + reach]))


_RUN_TRIM_SCORE = 10             # net matching comparisons (+1 match / -2 mismatch) that make a stretch of sequence a tandem run
_RUN_TRIM_SLACK = 10             # a seeded extension can stop short of, or overshoot, the run's edge: the run may start this far in
_RUN_TRIM_OUT = 40               # bases outside the arm that are read with it, so that a few run bases left inside the arm are seen


def _run_extent(seq: str, i: int, step: int, limit: int, x_drop: int = 14) -> int:
    """Bases of the arm that starts at seq[i] and runs inward by `step` (+1 from a 5' end, -1 from a 3' end) that belong
    to a tandem run of period 1-6 continuing from the arm's outer end. Each base is compared with the one a period
    further in; the run is a stretch of those comparisons that scores _RUN_TRIM_SCORE (+1 / -2, as _outer_start walks
    a diagonal) and starts within _RUN_TRIM_SLACK bases of the arm's end, and it ends where the score falls x_drop below
    its peak, so a few substitutions do not end it. The walk starts _RUN_TRIM_OUT bases outside the arm, so that the
    last few bases of a run that the arm only just entered are seen. Returns (bases in the arm, start of the run, which
    is negative when it begins outside, period); (0, 0, 0) when no such run enters the arm, as random sequence and a
    TA site do not reach the score."""
    n = len(seq)
    out = min(_RUN_TRIM_OUT, i if step > 0 else n - 1 - i)
    limit = min(limit, n - i if step > 0 else i + 1)
    line = seq[i - out:i + limit + 6] if step > 0 else seq[max(0, i - limit - 5):i + out + 1][::-1]    # line[k] is x = k - out
    clean = set(line) <= _ACGT
    runs = []
    for p in range(1, 7):
        cur = best = 0
        start, best_end, held = -out, -p, False
        for k in range(min(len(line) - p, out + limit)):
            u, v = line[k], line[k + p]
            if clean or (u in _ACGT and v in _ACGT):
                cur += 1 if u == v else -2
            x = k - out
            if held:
                if cur > best:
                    best, best_end = cur, x + 1
                elif cur < best - x_drop:
                    break
            elif cur <= 0:
                cur, start = 0, x + 1
                if start > _RUN_TRIM_SLACK:
                    break                                  # a run that starts here no longer reaches the arm's end
            elif cur >= _RUN_TRIM_SCORE:
                best, best_end, held = cur, x + 1, True
        if held and best_end + p > 0:
            runs.append((best, p, best_end + p, start))
            if best_end >= limit:
                break                                      # the run fills the walk: a longer period would only agree
    top = max((r[0] for r in runs), default=0)
    return next(((end, st, p) for score, p, end, st in runs if score >= 0.9 * top), (0, 0, 0))   # the smallest period that fits as well as the best: a multiple of it, or chance matches with a substituted base, would overshoot


def _unit(w: str, p: int) -> str:
    """The repeat unit of the period-p run w: the majority base at each phase."""
    return "".join(Counter(w[j::p]).most_common(1)[0][0] for j in range(p))


def _self_paired(w: str, p: int) -> bool:
    """The unit of the period-p run w is its own reverse complement up to rotation: (TA)n, (GATC)n and (AATT)n are,
    A(n), (CA)n and (CAG)n are not, and no odd period can be."""
    u = _unit(w, p)
    return reverse_complement(u) in u + u


def _short_or_paired(w: str, p: int) -> bool:
    """A period-4-6 reading of w counts as a run when its unit is palindromic or repeats a unit of period 1-3."""
    u = _unit(w, p)
    return reverse_complement(u) in u + u or any(u == u[:q] * (p // q) for q in (1, 2, 3) if p % q == 0)


def _run_trim(seq: str, s: int, e: int, cap: int) -> int:
    """Bases to cut from both ends of the arms [s, e) because their outer ends are a tandem run that pairs with its own
    reverse complement, or -1 when such a run fills every base an arm can have (`cap`, or half the span). A transposon
    inserted into a (TA)n run splits it, and its arms extend through the target site into the two halves, which pair
    with each other; that run is no part of the element, so it is cut rather than the whole candidate discarded. A
    run of any other unit, A(n) or (CA)n, cannot pair across a split, so at both arm ends it is the element's own.
    The longer of the two runs is cut from both ends, since the arms pair base for base."""
    limit = min(cap, (e - s) // 2 - MIN_TIR_LEN)
    if limit <= 0:
        return 0
    cut = 0
    for i, step in ((s, 1), (e - 1, -1)):
        x, st, p = _run_extent(seq, i, step, limit)
        if not x or not _self_paired(seq[max(0, i + st):i + x] if step > 0 else seq[max(0, i - x + 1):i - st + 1], p):
            continue
        if x >= limit:
            return -1
        cut = max(cut, x)
    return cut


def _tir_row(s: int, e: int, L: int, ident: float, method: str, lower_bound: bool = False) -> dict:
    ev = {"type": "TIR (terminal inverted repeat)", "tir_len": L, "identity": ident,
          "five_prime": [s, s + L], "three_prime": [e - L, e], "element_span": [s, e], "method": method}
    if lower_bound:                                      # the scan limit stopped it, not the repeat ending
        ev["length_is_lower_bound"] = True
        ev["bound_reason"] = ("the inverted repeat still matched at the scan limit — the length is a "
                              "lower bound, not a measurement")
    return ev


def _terminal_tir(seq: str, min_tir: int, max_tir: int, min_ident: float = MIN_TIR_IDENTITY):
    """Direct scan for a terminal inverted repeat anchored at the true element ends [0,L)/[n-L,n)."""
    n = len(seq)
    r = _scan_tir(seq, 0, n, min_tir, max_tir, min_ident)
    if not r:
        return None
    L, ident, _, hi = r
    return _tir_row(0, n, L, ident, "terminal inverted-repeat scan (element termini, %d bp limit)" % hi, L == hi)


def _outer_start(seq: str, rc: str, i: int, j: int, x_drop: int = 6):
    """Walk a seeded diagonal outward (toward both record ends) with an ungapped X-drop; (i0, j0) is the
    outermost pair still in the repeat. rc[j] pairs with seq[i] when they are equal."""
    score = best = best_off = off = 0
    while i - 1 - off >= 0 and j - 1 - off >= 0:
        a, b = seq[i - 1 - off], rc[j - 1 - off]
        if a in _ACGT and b in _ACGT:
            score += 1 if a == b else -2
        off += 1
        if score > best:
            best, best_off = score, off
        elif score < best - x_drop:
            break
    return i - best_off, j - best_off


def _core_run(seq: str, rc: str, s: int, j0: int, L: int):
    """(length, offset from the outer end) of the longest run of exactly paired bases in the arms."""
    run = top = off = 0
    for t in range(L):
        a = seq[s + t]
        run = run + 1 if a in _ACGT and a == rc[j0 + t] else 0
        if run > top:
            top, off = run, t - run + 1
    return top, off


def palindromic_target_site(seq: str, tir: dict, motif: str = "TA", max_trim: int = _TIR_TRIM_MAX):
    """Recover a palindromic target-site duplication that the inverted-repeat arms ran through.

    TA is its own reverse complement, so on a flanked Tc1/Mariner insertion the arms continue 2 bp through the
    TSD and tir_len comes out as TIR + 2. Trimming t outer bases from a TIR leaves a direct repeat whenever
    those bases are palindromic, so this cannot be done for every TIR (about one in four would gain a
    "TSD"): the caller must know the superfamily (classify, for Tc1/Mariner) and pass its motif.

    Returns (tir, tsd) with the arms trimmed to the smallest trim that exposes `motif` at both ends, or None.
    The trimmed row records the adjustment; nothing is changed in place."""
    if not tir or not tir.get("tir_len"):
        return None
    s, e = tir["element_span"]
    m, L = len(motif), tir["tir_len"]
    for a in range(0, min(max_trim, L - MIN_TIR_LEN) + 1):
        s2, e2 = s + a, e - a
        if s2 >= m and e2 + m <= len(seq) and seq[s2 - m:s2] == motif and seq[e2:e2 + m] == motif:
            row = dict(tir)
            if a:
                row.update({"tir_len": L - a, "five_prime": [s2, s2 + L - a], "three_prime": [e2 - L + a, e2],
                            "element_span": [s2, e2], "boundary_adjustment": [a, -a],
                            "identity": _identity(seq[s2:s2 + L - a], reverse_complement(seq[e2 - L + a:e2]))})
                row.pop("boundary_uncertain", None)
            return row, find_tsd(seq, s2, e2, expect=m, placements=a + 1)
    return None


def _tir_candidates(seq, rc, k, min_tir, max_tir, min_anchors):
    """Unanchored arm candidates from the two seed searches, one per diagonal, each as a dict with its
    outer ends, length and identity before any TSD-based trimming."""
    n = len(seq)
    seeds = {}
    for kk, window, need, limit in ((k, _TIR_SEED_WINDOW, min_anchors, max_tir), (_TIR_JOINT_K, _TIR_JOINT_WINDOW, _TIR_JOINT_ANCHORS, _TIR_JOINT_MAX)):
        W = min(n // 2, window)
        if W < kk or W < min_tir:
            continue
        idx = _kmer_index(seq[:W], kk, _SEED_MAX_OCC)
        idx_q = _kmer_index(rc[:W], kk, _SEED_MAX_OCC)
        hits = {}
        for j in range(W - kk + 1):
            w = rc[j:j + kk]
            if w not in idx_q:
                continue
            for i in idx.get(w, ()):
                hits.setdefault(j - i, []).append(i)
        found = {}
        for d, il in hits.items():
            if len(il) >= need and d not in seeds:   # seed the densest cluster, not a chance hit elsewhere on the diagonal
                srt = sorted(il)
                dens = {x: bisect_right(srt, x + kk + 4) - bisect_left(srt, x - kk - 4) for x in il}
                i = max(il, key=dens.get)
                if dens[i] >= need:
                    found[d] = (dens[i], ((i, i + d), limit))
        if len(found) > _TIR_MAX_SEEDS:              # repeat-rich window: keep the densest diagonals, in discovery order
            keep = set(sorted(found, key=lambda d: -found[d][0])[:_TIR_MAX_SEEDS])
            found = {d: v for d, v in found.items() if d in keep}
        seeds.update({d: v for d, (_, v) in found.items()})
    out = []
    for d, ((i, j), limit) in seeds.items():
        i0, j0 = _outer_start(seq, rc, i, j)
        s, e = i0, n - j0
        cut = _run_trim(seq, s, e, max_tir)
        if cut < 0:
            continue                                       # nothing but a tandem run to pair: not a terminal repeat
        s, e, j0 = s + cut, e - cut, j0 + cut
        r = _scan_tir(seq, s, e, min_tir, limit, x_drop=_TIR_XDROP, reach=0 if cut else 16)
        if r:
            L, ident, m, hi = r
            core, off = _core_run(seq, rc, s, j0, L)
            out.append({"s": s, "e": e, "L": L, "ident": ident, "hi": hi, "j0": j0, "core": core, "off": off,
                        "reach": 0 if cut else 16})
    return out


def find_tir(seq: str, k: int = 11, min_tir: int = MIN_TIR_LEN, max_tir: int = 1000, min_anchors: int = 3):
    """Detect terminal INVERTED repeats (TIR candidate, DNA transposons).
    The 1 kb limit spans the described range: ~11-30 bp for hAT and Tc1/Mariner, ~13 bp for CACTA,
    but hundreds of bp for MULE/Mutator and Foldback-type elements. A repeat still matching at the
    limit is reported as a lower bound rather than silently cut to the limit.

    A record whose ends already carry an inverted repeat is read at its ends, as in 4.0.0 (a chance repeat at the
    ends of a flanked record then stands in for the real one, about 1% of the time). Otherwise the arms are located
    from seeds (see _TIR_SEED_WINDOW), the arms are trimmed back when they ran through a palindromic target
    site, and the TSD that fixes the boundary is what lets a short or imperfect TIR be accepted at all. When
    several diagonals qualify the arms are chosen by corroboration, not by the number of anchors: a
    subterminal inverted repeat can out-count the true termini."""
    seq = seq.upper()
    n = len(seq)
    term = _terminal_tir(seq, min_tir, max_tir)     # a record that starts and ends on the element: read at its ends
    if term:
        return term
    rc = reverse_complement(seq)
    best, best_key, placements = None, None, 0
    for c in _tir_candidates(seq, rc, k, min_tir, max_tir, min_anchors):
        s, e, L = c["s"], c["e"], c["L"]
        shifts = [(a_, -a_) for a_ in range(-_TIR_OUT_MAX, min(_TIR_TRIM_MAX, L - min_tir) + 1) if a_]
        placements += len(shifts) + 1                      # every candidate's boundary placements are searched for a TSD
        pick = _snap_choice(seq, s, e, shifts, cap=MAX_TSD, min_tsd=5, paired=True)
        a = pick[0] if pick else 0
        s2, e2, L2 = s + a, e - a, L - a
        t = _tsd_at(seq, s2, e2, MIN_TSD, MAX_TSD)[0] if (s2 > 0 and e2 < n) else 0
        if t >= _TIR_WEAK_TSD and _tsd_in_microsatellite(seq, s2, e2, t):
            t = 0                                         # a TSD inside an (AT)n run is no evidence of an insertion site
        supported = t >= _TIR_WEAK_TSD
        strong = c["core"] >= min_anchors + k - 1
        if not (supported or strong):
            continue
        if a:
            ident = _identity(seq[s2:s2 + L2], reverse_complement(seq[e2 - L2:e2]))
            if ident < MIN_TIR_IDENTITY or L2 < min_tir or _simple_arms(seq, s2, e2, L2, c["reach"]):
                s2, e2, L2, a = s, e, L, 0
        if supported:
            # A candidate accepted on the strength of its TSD has its arms measured at _TIR_WEAK_IDENTITY:
            # at the 80% floor, chance TIR + chance TSD together extended arms to 19 bp across the chance_tir
            # seeds, one bp under MIN_TIR_STANDALONE, which would erode that floor's margin.
            r2 = _scan_tir(seq, s2, e2, min_tir, max_tir if strong else _TIR_JOINT_MAX, _TIR_WEAK_IDENTITY, _TIR_XDROP)
            if r2:
                L2 = r2[0]
            else:                                         # nothing reaches that identity: report the exact core only
                s2, e2, L2, a = s + c["off"], e - c["off"], c["core"], 0
        if not supported:
            # Arms placed by repeat strength alone are reported as the exact core that qualified them, not as
            # an extension: a seeded search over every start offset of a long flank reaches 20+ bp of chance
            # repeat once mismatches are allowed (chance_tir.py), so an extended length here would no longer
            # sit below MIN_TIR_STANDALONE. A corroborated candidate carries its TSD and is extended.
            s2, e2, L2, a = s + c["off"], e - c["off"], c["core"], 0
        if L2 < min_tir:
            continue                                      # an exact core shorter than the detection floor is not a TIR call
        if _simple_arms(seq, s2, e2, L2, c["reach"]):
            continue                                      # whichever path set the arms: a microsatellite pairs with its own reverse complement
        long_ = strong and L2 >= MIN_TIR_STANDALONE       # longer than chance reaches: outranks a corroborated short candidate
        key = (supported and long_, long_, supported, t if supported else 0, strong, L2)
        if best_key is None or key > best_key:
            best_key, best = key, (c, s2, e2, L2, a, supported)
    if not best:
        return None
    c, s, e, L, a, supported = best
    ident = _identity(seq[s:s + L], reverse_complement(seq[e - L:e])) if (a or s != c["s"] or L != c["L"]) else c["ident"]
    ev = _tir_row(s, e, L, ident, "k-mer seed vs reverse-complement (k=%d, %d bp window)" % (k, _TIR_SEED_WINDOW),
                  L == c["hi"] and L == c["L"])
    if a:
        ev["boundary_adjustment"] = [a, -a]
    if not supported:
        ev["boundary_uncertain"] = True              # arms placed by repeat strength alone; no flanking TSD fixes them
    else:
        ev["boundary_tests"] = placements            # accepted because of its TSD: that TSD's chance is adjusted for these
    return ev


# Poly-A tail search. A TPRT-made tail is not always the last thing in a record: a genomic window carries
# flank (and the target-site duplication) after it, and a real tail has the odd non-A. So the tail is the
# best-scoring A-rich tract in the terminal window (A +1, anything else -_TAIL_PENALTY: about one
# interruption per ten bases at the purity floor), and it counts when either
#   - it ends within _TAIL_SLACK bases of the record end (a bare element), or
#   - it is followed by a direct repeat of the sequence before the element, >= _TAIL_TSD_MIN bp: the
#     target-site duplication that frames a TPRT insertion. That corroboration is what licenses a tail that
#     is not at the end, since an A-rich tract of 8+ bp turns up by chance in a few hundred bp.
# Only the last _TAIL_WINDOW bases of the record are searched (the first, for the poly-T of a reverse-oriented
# element); a tail further from the end is not seen, and scan_scope() reports the window.
_TAIL_PENALTY = 9
_TAIL_SLACK = 3
_TAIL_WINDOW = 600
_TAIL_TSD_MIN = 8


def _a_tracts(seq: str, lo: int, hi: int, min_run: int):
    """A-rich tracts of seq[lo:hi] as (start, end), each starting and ending on an A, maximal under the
    +1 / -_TAIL_PENALTY score, with at least min_run A's."""
    out, start, score, best, best_end = [], None, 0, 0, 0
    for i in range(lo, hi):
        if seq[i] == "A":
            if start is None:
                start, score, best = i, 0, 0
            score += 1
            if score > best:
                best, best_end = score, i + 1
        elif start is not None:
            score -= _TAIL_PENALTY
            if score <= 0:
                out.append((start, best_end))
                start = None
    if start is not None:
        out.append((start, best_end))
    return [(a, b) for a, b in out if seq.count("A", a, b) >= min_run]


def _tail_tsd(seq: str, heads: dict, tail_start: int, end: int, n_tests: int = 1):
    """The direct repeat that frames a tail ending at `end`: a copy of seq[end:end+t] (t >= _TAIL_TSD_MIN)
    lying before the tail, as a TSD row, or None. It is accepted only when its family-wise chance over the
    `n_tests` placements tried (left copy x tail end x tract) is at or below CHANCE_LEVEL_P: an 8-mer match
    somewhere in ~3,000 placements is routine in random sequence, a 12-mer is not. The search is capped at
    MAX_TSD (12 bp) and the test is per strand, so in practice only a 12 bp copy corroborates a tail, and not in
    flanks below about 30% GC, where a perfect 12-mer still exceeds the chance level."""
    probe = seq[end:end + _TAIL_TSD_MIN]
    if len(probe) < _TAIL_TSD_MIN or not _ACGT.issuperset(probe):
        return None
    for p in heads.get(probe, ()):
        for t in range(MAX_TSD, _TAIL_TSD_MIN - 1, -1):       # the left copy begins at p; its length fixes the element start
            if p + t < tail_start and seq[p:p + t] == seq[end:end + t] and end + t <= len(seq):
                fam = 1 - (1 - _chance_p(seq, t)) ** n_tests
                if fam > CHANCE_LEVEL_P:
                    break                                      # a shorter copy at this position is weaker still
                return {"type": "TSD (target-site duplication)", "length": t, "motif": seq[p:p + t],
                        "upstream": [p, p + t], "downstream": [end, end + t], "matched_expected": False,
                        "chance_probability": float(f"{fam:.3g}"), "n_tests": n_tests, "chance_level": fam > CHANCE_LEVEL_P,
                        "context": "flanks a poly-A tail (target-primed reverse transcription)"}
    return None


def _find_tails(seq: str, min_run: int):
    """3' poly-A tails as (row, tsd row or None) pairs."""
    n = len(seq)
    out = []
    lo = max(0, n - _TAIL_WINDOW)
    heads = None
    tracts = _a_tracts(seq, lo, n, min_run)
    n_tracts = sum(1 for _, b in tracts if b < n - _TAIL_SLACK) or 1
    n_left = max(1, min(n, _TAIL_WINDOW) - _TAIL_TSD_MIN + 1)      # left-copy start positions searched
    for a, b in tracts:
        row = {"type": "poly-A tail (3')", "length": b - a, "pos": [a, b]}
        tsd = None
        if b >= n - _TAIL_SLACK:
            pass
        else:
            if heads is None:
                w = min(n, _TAIL_WINDOW)
                heads = {}
                for p in range(max(0, w - _TAIL_TSD_MIN) + 1):
                    heads.setdefault(seq[p:p + _TAIL_TSD_MIN], []).append(p)
            # the tract may have absorbed leading A's of the target site, so the true tail end is tried at each
            # of the last few A positions and the longest duplication wins (ties go to the longer tail)
            ends = range(b, max(a + min_run, b - _TAIL_SLACK - 2) - 1, -1)
            for end in ends:
                cand = _tail_tsd(seq, heads, a, end, len(ends) * n_left * n_tracts)
                if cand and (tsd is None or cand["length"] > tsd["length"]):
                    tsd, b = cand, end
            if not tsd:
                continue
            row.update({"length": b - a, "pos": [a, b]})
            if not tsd:
                continue
            row["tsd_corroborated"] = True
        out.append((row, tsd))
    return out[-1:] if out else []                           # one tail: the 3'-most that qualifies


def _tails_both(seq: str, min_run: int):
    """(rows, tsd rows) for the 3' poly-A tail and the 5' poly-T tract of a reverse-oriented element."""
    n = len(seq)
    rows, tsds = [], []
    for row, tsd in _find_tails(seq, min_run):
        rows.append(row)
        tsds += [tsd] if tsd else []
    rc = reverse_complement(seq)
    for row, tsd in _find_tails(rc, min_run):
        a, b = row["pos"]
        t = {"type": "poly-T (5')", "length": row["length"], "pos": [n - b, n - a]}
        if row.get("tsd_corroborated"):
            t["tsd_corroborated"] = True
        rows.append(t)
        if tsd:
            tsds.append({**tsd, "motif": reverse_complement(tsd["motif"]), "upstream": [n - tsd["downstream"][1], n - tsd["downstream"][0]],
                         "downstream": [n - tsd["upstream"][1], n - tsd["upstream"][0]]})
    return rows, tsds


POLYA_MIN_RUN = 8


def find_polya(seq: str, min_run: int = POLYA_MIN_RUN):
    """Poly-A (3') or poly-T (5') tail — LINE/retro signature. See the note above for what counts."""
    return _tails_both(seq, min_run)[0]


# Expected TSD length by superfamily, restricted to values that can be attributed and to superfamilies
# classify.py can actually emit. Anything absent here is "not diagnostic" — deliberately, rather than
# filled with plausible numbers. Each entry is (preferred length, accepted range, basis, required motif).
# LTR retrotransposons duplicate 4-6 bp (Wicker et al. 2007 Nat Rev Genet 8:973-982; the same range is
# reported for HERV insertions, e.g. the HML-2 survey in BMC Genomics 2017, doi:10.1186/s12864-017-3872-6),
# with 5 bp the preferred value in Ou et al. 2019 Genome Biol 20:275. hAT
# 8 bp: Ou et al. 2019. The Tc1/Mariner TA dinucleotide is the defining insertion site of the superfamily
# (Wicker et al. 2007), so for that superfamily the motif is part of the expectation, not only the length.
TSD_EXPECT = {
    "Copia": (5, (4, 6), "Ou et al. 2019 Genome Biol 20:275; 4-6 bp range of LTR retrotransposons, Wicker et al. 2007 Nat Rev Genet 8:973-982", None),
    "Gypsy": (5, (4, 6), "Ou et al. 2019 Genome Biol 20:275; 4-6 bp range of LTR retrotransposons, Wicker et al. 2007 Nat Rev Genet 8:973-982", None),
    "hAT": (8, (8, 8), "Ou et al. 2019 Genome Biol 20:275", None),
    "Tc1/Mariner": (2, (2, 2), "Wicker et al. 2007 Nat Rev Genet 8:973-982 (TA target site)", "TA"),
}


def tsd_congruence(length, superfamily, motif=None):
    """Compare an observed TSD with the superfamily's expected length range and, where one is defined,
    motif. `expected` stays the single preferred length (the value find_tsd is asked to prefer);
    `expected_range` is what counts as congruent. `motif` is optional: without it only length is judged.

    Only INCONGRUENT carries information. An exact L-mer direct repeat at a fixed offset arises by chance
    with probability ~4^-L under a uniform null — common for short L, and commoner still in the AT-rich
    flanks TEs favour — so a matching length corroborates weakly while a clashing one is uncertain evidence
    that the boundary, the superfamily call, or both are wrong."""
    key = (superfamily or "").split(" ")[0]
    for name, (exp, rng, cite, want_motif) in TSD_EXPECT.items():
        if key and (key == name.split("/")[0] or superfamily.startswith(name)):
            base = {"expected": exp, "expected_range": list(rng), "basis": cite}
            if length is None:
                return {"verdict": "not assessed", **base}
            ok = rng[0] <= length <= rng[1]
            out = {"verdict": "congruent" if ok else "incongruent", **base, "observed": length}
            if want_motif and motif is not None:
                out["expected_motif"], out["observed_motif"] = want_motif, motif
                if motif != want_motif:
                    out["verdict"] = "incongruent"
            return out
    return {"verdict": "not diagnostic", "expected": None,
            "basis": "no attributable expected TSD length for this superfamily"}


# The searched TSD length window, named so the panels that DESCRIBE this detector can interpolate the
# real bounds instead of restating them. They disagreed: the Methods disclosure and the Method column both
# claimed "4–12 bp" while the search has run from 2 bp since 3.2.1, so the screen told the user the exact
# Tc1/Mariner TA case this detector was tuned to catch could not be found. tests/test_docs_track_build.py
# now fails if a panel string and these constants drift apart again.
MIN_TSD, MAX_TSD = 2, 12

# A chance TSD is not evidence of an insertion site: P(an exact L-mer direct repeat) = (sum of squared base
# frequencies)^L at one placement, so 2-3 bp matches occur in ~6% and ~2% of random flank pairs. A row whose
# chance probability, adjusted for every length and boundary placement searched, exceeds this is flagged
# `chance_level`; a single pre-specified 5 bp match is at 1e-3 under uniform composition, and an unspecified
# one (any length searched) or a snapped boundary needs more.
CHANCE_LEVEL_P = 1e-3


def _tsd_at(seq: str, elem_start: int, elem_end: int, min_tsd: int, max_tsd: int, expect=None):
    """Longest exact ACGT-only direct repeat flanking [elem_start, elem_end) as (length, motif), the
    expected length tried first; (0, '') when none. Ambiguity codes are unknown bases, never a match."""
    order = [L for L in ([expect] if expect else []) if min_tsd <= L <= max_tsd]
    order += [L for L in range(max_tsd, min_tsd - 1, -1) if L not in order]
    for L in order:
        up_s, dn_e = elem_start - L, elem_end + L
        if up_s < 0 or dn_e > len(seq):
            continue
        left = seq[up_s:elem_start]
        if left == seq[elem_end:dn_e] and _ACGT.issuperset(left):
            return L, left
    return 0, ""


def _chance_p(seq: str, length: int) -> float:
    acgt = [seq.count(b) for b in "ACGT"]
    tot = sum(acgt)
    return float(sum((c / tot) ** 2 for c in acgt) ** length) if tot else 1.0


def find_tsd(seq: str, elem_start: int, elem_end: int, min_tsd: int = MIN_TSD, max_tsd: int = MAX_TSD, expect=None,
             placements: int = 1):
    """Target-site duplication: a short direct repeat immediately flanking the element.
    Requires flanking sequence beyond [elem_start, elem_end]; else returns None (honest).

    `expect` (a length, when the superfamily is known) is tested FIRST. Plain longest-first search lets a
    chance 7-mer in an AT-rich flank outrank the real 2 bp TA of a Tc1/Mariner element, which is also why
    min_tsd is 2 rather than 4 — at 4 the diagnostic TA was unreachable. The detect_all call sites run
    BEFORE the superfamily is known, so they pass no `expect`; classify.classify re-runs find_tsd with the
    resolved superfamily's expected length once it has one, correcting a coincidental longer pick in place.

    The row carries `chance_probability`, `n_tests` and `chance_level`. The probability is the chance, from the
    record's own composition, that a search this wide finds a repeat at least this long: the union over the
    lengths from L to max_tsd (one length when `expect` matched first) and over `placements`, the boundary
    placements searched to arrive at this element span (1 for a fixed span; an LTR or TIR boundary that was
    snapped jointly with its TSD passes the number of shifts it chose among). `n_tests` is lengths x placements."""
    L, left = _tsd_at(seq, elem_start, elem_end, min_tsd, max_tsd, expect)
    if not L:
        return None
    lengths = [L] if expect and L == expect else range(L, max_tsd + 1)
    miss = 1.0
    for l in lengths:
        miss *= 1 - _chance_p(seq, l)
    p = 1 - miss ** placements
    return {"type": "TSD (target-site duplication)", "length": L,
            "motif": left, "upstream": [elem_start - L, elem_start], "downstream": [elem_end, elem_end + L],
            "matched_expected": bool(expect and L == expect),
            "chance_probability": float(f"{p:.3g}"), "n_tests": len(lengths) * placements, "chance_level": p > CHANCE_LEVEL_P}


# Boundary snapping. An X-drop extension past the true end of a repeat keeps one chance-matching base about
# a quarter of the time, so a flanked element's boundary lands 1-4 bp too wide and its real TSD is then
# searched at the wrong place. The boundary is therefore re-chosen jointly with the TSD, within a few bp,
# the way LTRharvest (-vic) and LTR_retriever treat the terminal vicinity. Chance protection: a shifted
# boundary is taken only when its evidence reaches _SNAP_MIN_BP and beats the original boundary's, where
# evidence is the TSD length (counted up to _SNAP_CAP bp, the LTR range, so a lucky long match cannot buy a
# large move) plus _CANONICAL_BP when both LTR ends read the canonical TG..CA (a 1-in-256 chance event,
# worth about 4 bp of TSD), less half a bp per base moved.
_SNAP_START = range(-1, 5)       # shift of the 5' outer end: + moves it inward (the usual overshoot)
_SNAP_END = range(-4, 2)         # shift of the 3' outer end: - moves it inward
_SNAP_MIN_BP = 5
_SNAP_MIN_TSD = 4                  # the bottom of the LTR range: a shorter match is not worth moving a boundary for
_SNAP_CAP = 6
_CANONICAL_BP = 3
_SHIFT_COST = 0.5
# A boundary that already carries a TSD is not moved to buy one more base: a chance match at the neighbouring
# position extends a TSD by one about a quarter of the time. The move must gain this much beyond its cost.
_SNAP_MARGIN = 1.5


def _snap_choice(seq, s, e, shifts, canon=None, prefer=None, cap=_SNAP_CAP, min_tsd=_SNAP_MIN_TSD,
                 min_bp=_SNAP_MIN_BP, paired=False, max_tsd=MAX_TSD):
    """Best (a, b) among `shifts` ((a, b) pairs: a moves the start, b the end) or None to keep (s, e).
    `prefer` is the family's modal TSD length: one overshooting end leaves a prefix of the TSD matching, so a
    4 bp match beside a 5 bp candidate is explained by the shift, and the modal length breaks that tie."""
    def evidence(a, b):
        t = min(_tsd_at(seq, s + a, e + b, MIN_TSD, max_tsd)[0], cap)
        return t + (t == prefer) + (_CANONICAL_BP if canon and canon(a, b) else 0), t

    base = best_score = evidence(0, 0)[0]
    best = None
    for a, b in shifts:
        c = abs(a) if paired else abs(a) + abs(b)          # paired: one move of the arms, not two
        if not c:
            continue
        ev, t = evidence(a, b)
        score = ev - _SHIFT_COST * c
        if t >= min_tsd and ev >= min_bp and score >= base + _SNAP_MARGIN and score > best_score:
            best, best_score = (a, b), score
    return best


def _snap_ltr(seq: str, ev: dict) -> dict:
    """Re-choose an accepted LTR pair's outer boundaries jointly with the flanking TSD (see above)."""
    (l5s, l5e), (l3s, l3e) = ev["five_prime"], ev["three_prime"]
    canon = lambda a, b: _termini_motif(seq, l5s + a, l5e + b, l3s + a, l3e + b)["canonical"]
    shifts = [(a, b) for a in _SNAP_START for b in _SNAP_END]
    prefer = TSD_EXPECT["Copia"][0]
    pick = _snap_choice(seq, l5s, l3e, shifts, canon, prefer=prefer)
    if not pick:
        return ev
    bonus = pick != _snap_choice(seq, l5s, l3e, shifts, None, prefer=prefer)     # the canonical bonus decided the move
    a, b = pick
    l5s, l3s, l5e, l3e = l5s + a, l3s + a, l5e + b, l3e + b
    ev.update({"ltr_len": l5e - l5s, "five_prime": [l5s, l5e], "three_prime": [l3s, l3e],
               "element_span": [l5s, l3e], "termini": _termini_motif(seq, l5s, l5e, l3s, l3e),
               "boundary_adjustment": [a, b], "boundary_tests": len(shifts)})
    if bonus:
        ev["termini"]["selected_for_motif"] = True
    if "gapped" not in ev["method"]:
        ev["identity"] = _identity(seq[l5s:l5e], seq[l3s:l3e])
    return ev


# Primer-binding-site reference: PBS (viral +strand, 5'->3') = reverse complement of the 3'-terminal 18 nt
# of the priming tRNA. Endogenised proviruses carry a DIVERGED PBS, so a match is often partial — the tRNA
# identity is reported hedged, never as a hard call. Panel is Lys-anchored (HML-2 primes tRNA-Lys); the Lys3
# entry equals the canonical HIV-1 tRNA-Lys3 PBS (independently verifiable) and anchors the panel's accuracy.
_PRIMER_TRNA = {
    "tRNA-Lys3": "TGGCGCCCGAACAGGGAC",     # Lys-primed lentiviruses / HML-2-class ERVs
}


def _pbs_chance(seq: str, pbs: str, matches: int, windows: int) -> float:
    """P(at least one of `windows` leader windows matches the panel PBS in >= `matches` of its positions by
    chance), from the record's own base composition (Poisson-binomial over the 18 positions). The reporting
    floor (55%) is crossed by chance in about one random 200 nt leader in seven (more in GC-rich sequence);
    this puts that rate beside the identity instead of leaving it to the Methods text."""
    tot = sum(seq.count(b) for b in "ACGT") or 1
    freq = {b: seq.count(b) / tot for b in "ACGT"}
    dist = [1.0]
    for ch in pbs:
        q = freq.get(ch, 0.25)
        dist = [(dist[i] if i < len(dist) else 0.0) * (1 - q) + (dist[i - 1] * q if i else 0.0) for i in range(len(dist) + 1)]
    p = sum(dist[matches:])
    return 1 - (1 - p) ** max(1, windows)


def find_pbs(seq: str, ltr_five_prime_end: int, search: int = 44, min_ident: float = 0.55):
    """Primer-binding site: the ~18 nt just 3' of the 5' LTR that the priming tRNA anneals to. Detected
    by best reverse-complement match to a bundled primer-tRNA panel, at the canonical leader position.
    tRNA identity is hedged when the match is weak (expected for endogenised, diverged proviruses).
    The panel holds tRNA-Lys3 only, so an element primed by another tRNA (initiator Met in Ty1/copia) is not
    assessed; every row carries `chance_probability`, the composition-adjusted chance of its match."""
    leader = seq[ltr_five_prime_end:ltr_five_prime_end + search]
    best = None
    for name, pbs in _PRIMER_TRNA.items():
        L = len(pbs)
        for i in range(len(leader) - L + 1):
            w = leader[i:i + L]
            m = sum(1 for a, b in zip(w, pbs) if a == b and a in _ACGT)
            if best is None or m / L > best["_id"]:
                best = {"_id": m / L, "m": m, "L": L, "trna": name, "pos": [ltr_five_prime_end + i, ltr_five_prime_end + i + L],
                        "motif": w, "pbs": pbs, "windows": len(leader) - L + 1}
    if not best or best["_id"] < min_ident:            # no credible PBS in the leader window
        return None
    # the panel is Lys-anchored, so only NAME the priming tRNA on a strong match; a weak (diverged) match
    # is reported as undetermined with the closest panel match, never a hard call for a non-Lys genus.
    strong = best["_id"] >= 0.72
    chance = _pbs_chance(seq, best["pbs"], best["m"], best["windows"])
    return {"type": "PBS (primer-binding site)", "pos": best["pos"],
            "priming_trna": best["trna"] if strong else "undetermined",
            "best_match": best["trna"], "identity": round(100 * best["_id"], 1),
            "confident": strong, "motif": best["motif"], "chance_probability": float(f"{chance:.3g}"),
            "panel": list(_PRIMER_TRNA),
            "note": "" if strong else (f"priming tRNA undetermined — closest panel match {best['trna']} "
                                       f"({round(100 * best['_id'], 1)}%), below the confident threshold; "
                                       f"a match this good arises by chance in about {100 * chance:.0f}% of leaders "
                                       "of this composition; endogenised PBS is often diverged")}


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


# Polyadenylation-signal hexamer panel. AATAAA and ATTAAA together account for roughly three quarters of
# MAMMALIAN poly(A) sites; the remainder are single-base variants of AATAAA (Beaudoing et al. 2000 Genome
# Res 10:1001-1010). Tiers are reported, never summed into a score: a hexamer's identity says how COMMON
# the variant is, not whether THIS copy is used.
# SCOPE, measured: this panel is derived from mammalian 3'-end data, but the detector runs on whatever
# element the user loads. Benchmarking across organisms returned hexamer hits inside yeast (Ty1) and plant
# (Ta1, Tnt1) LTRs, where 3'-end formation does NOT follow the mammalian AATAAA consensus — budding-yeast
# and plant signals are degenerate and A-rich rather than a fixed hexamer. The panel's provenance therefore
# travels with every call (see the note built in find_polya_signal); it is never presented as a validated
# signal outside mammals.
_POLYA_HEXAMERS = {
    "AATAAA": "canonical",
    "ATTAAA": "principal variant",
    # single-base variants of AATAAA catalogued by Beaudoing et al. 2000; each is individually uncommon
    "AGTAAA": "minor variant", "TATAAA": "minor variant", "CATAAA": "minor variant",
    "GATAAA": "minor variant", "AATATA": "minor variant", "AATACA": "minor variant",
    "AATAGA": "minor variant", "ACTAAA": "minor variant", "AAGAAA": "minor variant",
    "AATGAA": "minor variant",
}


def _dse_state(seq: str, hex_end: int, lo: int, hi: int, min_gt_frac: float, min_t: int) -> tuple:
    """Downstream sequence element (DSE): the GU/U-rich stretch that follows a real poly(A) signal.

    Three outcomes, like _tsd_state. A hexamer alone is not evidence — AATAAA-type words occur about once
    per 4 kb by chance and roughly half of genomic occurrences are not used as poly(A) signals, so the
    DSE requirement is what makes a hit worth reporting. When the record ends before the window, the
    question was never asked: that is 'not assessable', never 'absent'.

    WINDOW, derived rather than guessed: the hexamer lies 10-30 nt UPSTREAM of the cleavage site, and the
    DSE begins roughly 15 nt DOWNSTREAM of that site and runs on for tens of bases, U/GU-rich (Zhao,
    Hyman & Moore 1999 Microbiol Mol Biol Rev 63:405; Legendre & Gautheret 2003 NAR 31:1375; Proudfoot
    2011 Genes Dev 25:1770). Measured from the hexamer's 3' end that places the DSE at about +25 to +60,
    which is what `lo`/`hi` default to. An earlier draft used +10 to +30 — the gap between the hexamer and
    the cleavage site — and therefore scored the wrong stretch of sequence entirely."""
    w_s, w_e = hex_end + lo, hex_end + hi
    if w_e > len(seq):
        return "not assessable (record ends before the downstream window)", None
    win = seq[w_s:w_e].upper()
    if "N" in win:
        return "not assessable (ambiguous bases in the downstream window)", None
    t = win.count("T")
    gt = (t + win.count("G")) / len(win)
    det = {"window": [w_s, w_e], "gt_frac": round(gt, 2), "t_count": t}
    if gt >= min_gt_frac and t >= min_t:
        return "found", det
    return "absent (window present, not GU/U-rich)", det


def find_polya_signal(seq: str, ltr: dict, lo: int = 20, hi: int = 60,
                      min_gt_frac: float = 0.65, min_t: int = 4):
    """Advisory poly(A)-signal motif inside the 3' LTR copy — the copy whose signal is used in vivo.

    ADVISORY ONLY. This reports a MOTIF, not a function: it never locates the R/U5 boundary, the cleavage
    site, or a transcription end, all of which are defined by processing events that sequence alone cannot
    place. Reported only when a hexamer from the panel is followed by a GU/U-rich downstream element
    (`lo`-`hi` nt past the hexamer, >= `min_gt_frac` G+T and >= `min_t` T) — the hexamer alone is not
    evidence. See _dse_state for where the window comes from.

    Because the two LTR copies are duplicates, the same relative offset is checked in the 5' copy: agreement
    corroborates, disagreement is reported as drift and lowers confidence rather than being averaged away.
    Returns None when no hexamer is present at all."""
    l3s, l3e = ltr["three_prime"]
    l5s, l5e = ltr["five_prime"]
    best = None
    for i in range(l3s, max(l3s, l3e - 6) + 1):
        hexa = seq[i:i + 6].upper()
        tier = _POLYA_HEXAMERS.get(hexa)
        if not tier:
            continue
        state, det = _dse_state(seq, i + 6, lo, hi, min_gt_frac, min_t)
        rank = (0 if state == "found" else 1, 0 if tier == "canonical" else (1 if tier == "principal variant" else 2))
        if best is None or rank < best["_rank"]:
            best = {"_rank": rank, "pos": [i, i + 6], "motif": hexa, "variant": tier,
                    "dse_state": state, "dse": det, "offset_in_ltr": i - l3s}
    if best is None or best["dse_state"].startswith("absent"):
        # The gate is the whole point: a hexamer with a assessable-but-not-GU/U-rich downstream window is
        # the chance match this detector exists to reject, so nothing is reported. A window that could not
        # be assessed at all is still reported (hedged) — that is "the question could not be asked", the
        # same three-outcome discipline _tsd_state uses, not a negative result.
        return None
    best.pop("_rank", None)                               # sort key, not evidence
    mirror = seq[l5s + best["offset_in_ltr"]: l5s + best["offset_in_ltr"] + 6].upper()
    both = mirror == best["motif"]
    confident = best["dse_state"] == "found" and both
    # first token is the export Name and must not read like the LINE "poly-A" TAIL — a different feature
    # with a different SO term. "polyA-signal" keeps the two distinguishable in a genome browser.
    return {"type": "polyA-signal (motif, advisory)", **{k: v for k, v in best.items() if not k.startswith("_")},
            "in_both_ltr_copies": both, "five_prime_copy_motif": mirror, "confident": confident,
            "method": ("poly(A)-signal hexamer panel (Beaudoing et al. 2000 Genome Res 10:1001) in the 3' LTR, "
                       f"gated on a GU/U-rich downstream element {lo}-{hi} nt downstream "
                       f"(>= {int(100 * min_gt_frac)}% G+T, >= {min_t} T)"),
            "note": ("motif only — advisory. The hexamer panel is derived from MAMMALIAN poly(A) sites; yeast, "
                     "plant and insect 3'-end formation uses degenerate signals this panel does not model, so "
                     "a hit outside mammals is a sequence coincidence unless independently supported. A poly(A) "
                     "hexamer occurs by chance about once per 4 kb and about "
                     "half of genomic occurrences are never used, so this does NOT locate the R/U5 boundary, "
                     "the cleavage site, or the transcript end; those need RNA evidence this tool does not use."
                     + ("" if both else " The 5' LTR copy does not carry the same motif at the same offset, "
                        "which is consistent with post-integration drift — treat the call as weaker."))}


def _tsd_state(seq: str, start: int, end: int, found, min_tsd: int = MIN_TSD, uncertain: bool = False) -> str:
    """Four outcomes, not two. A TSD can only be looked for when the record extends beyond the element
    on BOTH sides; a bare pasted element — the normal case — has no flanks, so 'no TSD' there means the
    question was never asked. Reporting that as 'absent' would be evidence the run does not have. Flanks
    shorter than the longest TSD searched can hold only part of the question, and a boundary that was placed
    by repeat strength alone is not a boundary a TSD can be tested against."""
    if found:
        return "found"
    room = min(start, len(seq) - end)
    if room < min_tsd:
        return "not assessable (no flanking sequence in the record)"
    if uncertain:
        return "not assessed (element boundary uncertain: the repeat arms were placed by repeat strength alone)"
    if room < MAX_TSD:
        return f"partially assessable (flanks of {room} bp allow only TSDs of up to {room} bp)"
    return "absent (flanks present, no direct repeat found)"


def _rc_ltr(ltr: dict, n: int) -> dict:
    """The LTR pair as seen on the reverse-complement strand: copies swap roles."""
    (l5s, l5e), (l3s, l3e) = ltr["five_prime"], ltr["three_prime"]
    return {"five_prime": [n - l3e, n - l3s], "three_prime": [n - l5e, n - l5s], "element_span": [n - l3e, n - l5s]}


def _cis_on(seq: str, ltr: dict, strand: str) -> list:
    """PBS, PPT and polyA-signal rows for one element orientation, in record coordinates. On '-' the
    detectors run on the reverse complement (where the element reads 5'->3') and are mapped back; the
    motifs stay as read on the element's own strand."""
    n = len(seq)
    s, l = (seq, ltr) if strand == "+" else (reverse_complement(seq), _rc_ltr(ltr, n))
    rows = [r for r in (find_pbs(s, l["five_prime"][1]), find_ppt(s, l["three_prime"][0]), find_polya_signal(s, l)) if r]
    if strand == "-":
        flip = lambda p: [n - p[1], n - p[0]]
        for r in rows:
            r["pos"] = flip(r["pos"])
            if r.get("dse"):
                r["dse"] = {**r["dse"], "window": flip(r["dse"]["window"])}
    for r in rows:
        r["strand"] = strand
    return rows


def _cis_evidence(rows: list) -> int:
    """How much one orientation's cis-elements say about the element's orientation: a polypurine tract, a
    confident PBS, a poly(A) signal with its downstream element."""
    return sum(1 for r in rows if r["type"].startswith("PPT") or (r["type"].startswith("PBS") and r.get("confident"))
               or (r["type"].startswith("polyA") and r.get("dse_state") == "found"))


def _cis_elements(seq: str, ltr: dict, orientation=None) -> list:
    """LTR cis-elements for the element's orientation, which a record does not state. A caller that knows it
    (from the coding strand) passes it; otherwise both strands are read and the one whose cis-elements say
    more wins. If neither does, or they tie, the elements are withheld rather than reported on a strand the
    evidence does not support, and the LTR row says so."""
    ltr["pbs_scope"] = ("PBS search covers " + ", ".join(_PRIMER_TRNA) + " only; an element primed by another "
                        "tRNA (e.g. initiator tRNA-Met in Ty1/copia) is not assessed, so no PBS row is not 'no PBS'")
    if orientation in ("+", "-"):
        ltr["cis_orientation"], ltr["cis_orientation_source"] = orientation, "caller"
        return _cis_on(seq, ltr, orientation)
    fwd, rev = _cis_on(seq, ltr, "+"), _cis_on(seq, ltr, "-")
    ef, er = _cis_evidence(fwd), _cis_evidence(rev)
    if ef == er:
        ltr["cis_orientation"] = "undetermined"
        ltr["cis_note"] = ("orientation of the element is undetermined (no polypurine tract, confident PBS or "
                           "poly(A) signal favours either strand), so PBS, PPT and polyA-signal are not reported")
        return []
    ltr["cis_orientation"], ltr["cis_orientation_source"] = ("+" if ef > er else "-"), "cis-elements"
    return fwd if ef > er else rev


def detect_all(seq: str, orientation=None):
    """Run all structural detectors. Returns list of evidence dicts (empty if none).
    `orientation` ('+' or '-', the strand the element reads on) is optional; see _cis_elements."""
    seq = seq.upper()
    ev = []
    near = []
    ltr = find_ltr(seq, near_miss=near)
    if not ltr and near:
        # no accepted LTR, but a real diagonal pair was found and rejected on identity alone — say so rather
        # than let the panel read "no terminal repeats detected"
        ev.extend(near)
    if ltr:
        ev.append(ltr)
        span = ltr["element_span"]
        tsd = find_tsd(seq, span[0], span[1], placements=ltr.get("boundary_tests", 1))
        ltr["tsd_state"] = _tsd_state(seq, span[0], span[1], tsd)
        if tsd:
            ev.append(tsd)
        ev.extend(_cis_elements(seq, ltr, orientation))  # advisory; never a boundary or a cleavage site
    # LTR (direct) and TIR (inverted) terminal architectures are mutually exclusive:
    # only look for a TIR when no LTR was found, so an LTR element never reports a spurious TIR.
    tir = find_tir(seq) if not ltr else None
    if tir:
        ev.append(tir)
        span = tir["element_span"]
        tsd = find_tsd(seq, span[0], span[1], placements=tir.get("boundary_tests", 1))
        tir["tsd_state"] = _tsd_state(seq, span[0], span[1], tsd, uncertain=bool(tir.get("boundary_uncertain")))
        if tsd:
            ev.append(tsd)
    tails, tail_tsds = _tails_both(seq, POLYA_MIN_RUN)
    ev.extend(tails)
    if not (ltr or tir):                              # a TSD that frames a tail is for a non-LTR element; never let it stand in for an LTR/TIR one
        ev.extend(tail_tsds)
    return ev
