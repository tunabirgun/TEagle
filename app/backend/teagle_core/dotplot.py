"""Self-similarity of one locus: exact k-mer matches binned into a square matrix.

What this is for. Every structural detector in TEagle looks for a specific thing — a terminal direct
repeat, an inverted repeat, a tail. Each can therefore miss. A self-similarity matrix looks for repeats
of ANY arrangement without being told what to expect, so a strong block where nothing was called is a
positive observation: it is the one view that can catch a false negative in a capped detector.

What it is NOT. Exact k-mer matching, not a local-alignment self-BLAST. At 80% identity an exact 13-mer
survives with probability 0.8^13 ~ 5.5%, so a diverged repeat produces a faint diagonal or none at all.
A missing diagonal therefore cannot refute a call — only a present one is evidence. The asymmetry is
deliberate and must be stated wherever the panel is shown.

Two signals, computed together because they answer different questions:
  - FORWARD matches (i, j) with i != j lie on diagonals parallel to the main one — direct repeats, the
    LTR signature.
  - REVERSE-COMPLEMENT matches lie on ANTI-diagonals — inverted repeats, the TIR signature.

Pure stdlib: numpy, PIL and matplotlib are excluded from the frozen bundle by the build guard, and
pulling any of them back in would undo the v3.0.0 installer reduction.

Reimplemented from the published algorithm (a dot matrix is Gibbs & McIntyre 1970); TE-Aid renders an
equivalent panel and is MIT-licensed, and no GPL EMBOSS code was consulted or translated.
"""
from __future__ import annotations

from .sequtil import reverse_complement

# A microsatellite makes a k-mer occur thousands of times, and every occurrence pairs with every other:
# the work is quadratic in the occurrence count, so an unmasked (AT)n input hangs the interface. Masking
# above a cap bounds the work AND is the scientifically right call — a k-mer present that many times
# carries no positional information.
DEFAULT_MAX_OCC = 120
DEFAULT_K = 13
DEFAULT_BINS = 300
# hard ceiling on pair work, so a pathological input degrades to a truncated picture rather than a hang
MAX_PAIRS = 4_000_000


def _index(seq: str, k: int):
    idx = {}
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        if "N" in kmer:
            continue
        idx.setdefault(kmer, []).append(i)
    return idx


def self_matrix(seq: str, k: int = DEFAULT_K, bins: int = DEFAULT_BINS,
                max_occ: int = DEFAULT_MAX_OCC, include_diagonal: bool = True) -> dict:
    """Binned self-similarity counts for `seq`.

    Returns a dict with `forward` and `reverse` matrices (lists of rows, counts per bin), the binning
    metadata, and an honest account of what was masked or truncated — a picture that silently dropped
    half its input would be read as an absence of repeats.
    """
    seq = (seq or "").upper()
    n = len(seq)
    bins = max(16, min(int(bins), 1200))
    out = {"k": k, "bins": bins, "length": n, "max_occ": max_occ,
           "forward": [[0] * bins for _ in range(bins)],
           "reverse": [[0] * bins for _ in range(bins)],
           "masked_kmers": 0, "masked_positions": 0, "truncated": False,
           "forward_max": 0, "reverse_max": 0}
    if n < k + 1:
        return out
    idx = _index(seq, k)
    scale = bins / float(n)

    def _bin(p):
        b = int(p * scale)
        return b if b < bins else bins - 1

    fwd, rev = out["forward"], out["reverse"]
    pairs = 0
    # ---- forward: direct repeats ----
    for kmer, pos in idx.items():
        if len(pos) > max_occ:
            out["masked_kmers"] += 1
            out["masked_positions"] += len(pos)
            continue
        for a in range(len(pos)):
            ba = _bin(pos[a])
            start = a if include_diagonal else a + 1
            for b in range(start, len(pos)):
                bb = _bin(pos[b])
                fwd[ba][bb] += 1
                if ba != bb:
                    fwd[bb][ba] += 1
                pairs += 1
        if pairs > MAX_PAIRS:
            out["truncated"] = True
            break
    # ---- reverse complement: inverted repeats ----
    if not out["truncated"]:
        rc = reverse_complement(seq)
        for i in range(n - k + 1):
            kmer = rc[i:i + k]
            if "N" in kmer:
                continue
            hits = idx.get(kmer)
            if not hits or len(hits) > max_occ:
                continue
            # position i in the reverse complement corresponds to forward position n-k-i
            fpos = n - k - i
            bf = _bin(fpos)
            for p in hits:
                rev[bf][_bin(p)] += 1
                pairs += 1
            if pairs > MAX_PAIRS:
                out["truncated"] = True
                break
    out["forward_max"] = max((max(r) for r in fwd), default=0)
    out["reverse_max"] = max((max(r) for r in rev), default=0)
    out.update(chance_floor(seq, k, bins, n))
    return out


def chance_floor(seq: str, k: int, bins: int, n: int) -> dict:
    """How much signal this sequence would show if it contained no repeat at all.

    Two random positions share a k-mer with probability (Σ p_b²)^k under the sequence's OWN base
    composition — an AT-rich element collides far more often than a uniform null would predict, and TEs
    are frequently AT-rich. With (n/bins)² position pairs per bin, that gives an expected count per bin.
    Reporting it is what separates 'a repeat is present' from 'this is what noise looks like at this word
    size': at k=8 on a 4.5 kb element the floor is high enough to fill the panel with speckle, which a
    reader would otherwise mistake for structure."""
    counts = {b: seq.count(b) for b in "ACGT"}
    tot = sum(counts.values()) or 1
    p_same = sum((c / tot) ** 2 for c in counts.values())
    per_pair = p_same ** k
    pairs_per_bin = (n / float(bins)) ** 2
    exp = per_pair * pairs_per_bin
    return {"chance_per_bin": round(exp, 3),
            "chance_total": round(per_pair * n * n / 2.0, 1),
            "base_collision_p": round(p_same, 4)}


def suggest_k(structural_ev, default: int = DEFAULT_K) -> int:
    """Pick a word size the detected repeats can actually produce.

    A repeat SHORTER than k cannot generate a single exact k-mer match, so it is invisible however
    strong it is. Measured on the maize Ac element, whose TIR is 11 bp at 90.9% identity: the reverse
    signal is 600 at k=8, 32 at k=11 and 2 at k=13 — the default word size hides it completely, while
    Tc1's 54 bp TIR is visible at every k. When TEagle has already measured a short terminal repeat, the
    word size follows it instead of the other way round."""
    lens = [ev.get("tir_len") or ev.get("ltr_len") for ev in (structural_ev or [])]
    lens = [L for L in lens if L]
    if not lens:
        return default
    shortest = min(lens)
    # keep k well below the repeat so a mismatch inside it cannot erase every window
    return max(6, min(default, shortest - 3))


def scope_note(m: dict) -> str:
    """The sentence that must travel with the picture."""
    bits = [f"Exact {m['k']}-mer self-similarity, not a local-alignment self-BLAST — at 80% identity an "
            f"exact {m['k']}-mer survives with probability ~{0.8 ** m['k'] * 100:.1f}%, so a faint or "
            f"absent diagonal is not evidence that no repeat exists.",
            f"A repeat shorter than {m['k']} bp cannot appear at all, however strong: lower the word "
            f"size to see short terminal inverted repeats."]
    if m.get("chance_per_bin") is not None:
        cpb = m["chance_per_bin"]
        bits.append(f"Expected by chance at this word size and base composition: {cpb:.2f} matches per "
                    f"bin. " + ("Most of the scattered signal here is that noise floor, not structure — "
                                "raise the word size to suppress it."
                                if cpb >= 0.5 else
                                "Any visible block is therefore above the noise floor."))
    if m.get("masked_kmers"):
        bits.append(f"{m['masked_kmers']} k-mer(s) occurring more than {m['max_occ']} times were masked "
                    f"({m['masked_positions']} positions) — simple/satellite repeat, which carries no "
                    f"positional information and would otherwise dominate the picture.")
    if m.get("truncated"):
        bits.append("The comparison was truncated at the work limit; the picture is incomplete.")
    return " ".join(bits)


def above_chance(m: dict, expected_false: float = 1.0) -> int:
    """Smallest per-bin count that should not appear anywhere by chance.

    A per-bin significance test is the wrong instrument here: the matrix has bins² cells, so an event
    that is rare in one cell is commonplace across the panel. On the maize Ac element at k=8 the mean is
    0.019 matches per bin, which makes a count of 1 look 'significant' cell-by-cell — yet 19,600 cells
    produce roughly 370 of them by chance, which is exactly the speckle that fills the picture. This
    returns the count c for which the WHOLE matrix is expected to show fewer than `expected_false` cells,
    so a reader knows which cells carry information.

    Returned, never applied: filtering the picture silently would hide the noise a reader needs in order
    to judge the signal."""
    mu = m.get("chance_per_bin") or 0.0
    cells = (m.get("bins") or 1) ** 2
    if mu <= 0:
        return 1
    # survival function of a Poisson with small mean, accumulated term by term
    term = pow(2.718281828459045, -mu)          # P(X = 0)
    tail = 1.0 - term
    c = 1
    while c < 64:
        if cells * tail < expected_false:
            return c
        term *= mu / c                           # P(X = c)
        tail -= term
        c += 1
    return c


def density_rows(m: dict, which: str = "forward"):
    """Matrix rows normalised to 0..1 for the heat map, with the peak reported separately so a caller
    can label the scale with real counts rather than an unlabelled gradient."""
    mat = m[which]
    peak = m[f"{which}_max"] or 1
    return [[c / peak for c in row] for row in mat], m[f"{which}_max"]
