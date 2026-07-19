"""Real structural TE evidence (Layer D): terminal direct repeats (LTR), terminal
inverted repeats (TIR), poly-A tails, flanking target-site duplications (TSD).
Heuristic detectors — reported as candidate structural evidence, never as family/DB calls."""
from __future__ import annotations
from collections import Counter
from .sequtil import reverse_complement


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


def find_ltr(seq: str, k: int = 13, min_ltr: int = 80, min_anchors: int = 4):
    """Detect a terminal DIRECT repeat pair (LTR candidate). Coords 0-based half-open."""
    n = len(seq)
    W = min(n // 2, 1800)
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
    if len(set(p5s)) < 0.4 * (L - k + 1):            # require a contiguous repeat, not scattered matches
        return None
    l3s, l3e = best_d + l5s, best_d + l5e
    if l3e > n:
        l3e = n
    ident = _identity(seq[l5s:l5e], seq[l3s:l3e])
    if ident < 80:
        return None
    return {"type": "LTR (terminal direct repeat)", "ltr_len": L, "identity": ident,
            "five_prime": [l5s, l5e], "three_prime": [l3s, l3e],
            "element_span": [l5s, l3e], "method": "k-mer seed + diagonal cluster (k=%d)" % k}


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
    best_L, best_score, best_ident = None, 0, 0.0
    for L in range(min_tir, hi + 1):
        a, b = seq[:L], b_all[:L]
        matches = sum(1 for x, y in zip(a, b) if x == y and x not in "Nn")   # N is not base-pairing evidence
        score = 2 * matches - L                          # cumulative match minus mismatch; peaks at the TIR boundary
        ident = 100.0 * matches / L
        if ident >= min_ident and score > best_score:    # extend while each added base pair keeps the score climbing
            best_L, best_score, best_ident = L, score, round(ident, 1)
    if best_L:
        return {"type": "TIR (terminal inverted repeat)", "tir_len": best_L, "identity": best_ident,
                "five_prime": [0, best_L], "three_prime": [n - best_L, n], "element_span": [0, n],
                "method": "terminal inverted-repeat scan (element termini)"}
    return None


def find_tir(seq: str, k: int = 11, min_tir: int = 10, max_tir: int = 60, min_anchors: int = 3):
    """Detect terminal INVERTED repeats (TIR candidate, DNA transposons)."""
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
    if ident < 80:
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


def find_tsd(seq: str, elem_start: int, elem_end: int, min_tsd: int = 4, max_tsd: int = 12):
    """Target-site duplication: a short direct repeat immediately flanking the element.
    Requires flanking sequence beyond [elem_start, elem_end]; else returns None (honest)."""
    for L in range(max_tsd, min_tsd - 1, -1):
        up_s, up_e = elem_start - L, elem_start
        dn_s, dn_e = elem_end, elem_end + L
        if up_s < 0 or dn_e > len(seq):
            continue
        left, right = seq[up_s:up_e], seq[dn_s:dn_e]
        if left == right and "N" not in left:
            return {"type": "TSD (target-site duplication)", "length": L,
                    "motif": left, "upstream": [up_s, up_e], "downstream": [dn_s, dn_e]}
    return None


def detect_all(seq: str):
    """Run all structural detectors. Returns list of evidence dicts (empty if none)."""
    ev = []
    ltr = find_ltr(seq)
    if ltr:
        ev.append(ltr)
        tsd = find_tsd(seq, ltr["element_span"][0], ltr["element_span"][1])
        if tsd:
            ev.append(tsd)
    # LTR (direct) and TIR (inverted) terminal architectures are mutually exclusive:
    # only look for a TIR when no LTR was found, so an LTR element never reports a spurious TIR.
    tir = find_tir(seq) if not ltr else None
    if tir:
        ev.append(tir)
        tsd = find_tsd(seq, tir["five_prime"][0], tir["three_prime"][1])
        if tsd:
            ev.append(tsd)
    ev.extend(find_polya(seq))
    return ev
