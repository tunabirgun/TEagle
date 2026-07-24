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


def detect_all(seq: str):
    """Run all structural detectors. Returns list of evidence dicts (empty if none)."""
    ev = []
    ltr = find_ltr(seq)
    if ltr:
        ev.append(ltr)
        tsd = find_tsd(seq, ltr["element_span"][0], ltr["element_span"][1])
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
        if tsd:
            ev.append(tsd)
    ev.extend(find_polya(seq))
    return ev
