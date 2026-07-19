"""Real sequence utilities: FASTA parse, IUPAC validation, composition, revcomp, ORFs."""
from __future__ import annotations
import re, hashlib

IUPAC = set("ACGTURYSWKMBDHVN")
_COMP = str.maketrans("ACGTURYSWKMBDHVNacgturyswkmbdhvn",
                      "TGCAAYRSWMKVHDBNtgcaayrswmkvhdbn")

def parse_fasta(text):
    """Return list of (header, sequence). Bare sequence (no '>') becomes one record.
    Tolerant of non-string input (None -> [], other types coerced) so a malformed
    request never crashes the parser."""
    if text is None:
        return []
    if not isinstance(text, str):
        text = str(text)
    text = text.replace("\r\n", "\n").replace("\r", "\n")
    recs = []
    if ">" not in text:
        seq = _norm(text)
        if seq:
            recs.append(("input_sequence", seq))
        return recs
    header, buf = None, []
    for line in text.split("\n"):
        if line.startswith(">"):
            if header is not None:
                recs.append((header, _norm("".join(buf))))
            header, buf = line[1:].strip() or "record", []
        elif header is not None:
            buf.append(re.sub(r"\s+", "", line))
    if header is not None:
        recs.append((header, _norm("".join(buf))))
    return recs


def _norm(s: str) -> str:
    # strip whitespace, uppercase, and normalize RNA to its DNA equivalent (U -> T) so composition,
    # ORF finding, translation, and primer design all operate on DNA
    return re.sub(r"\s+", "", s).upper().replace("U", "T")

def validate_iupac(seq: str):
    """Return (ok, [(0-based position, char), ...]) for characters outside the IUPAC set."""
    bad = [(i, c) for i, c in enumerate(seq) if c not in IUPAC]
    return (len(bad) == 0, bad[:50])

def composition(seq: str):
    n = len(seq)
    if n == 0:
        return {"length": 0, "gc": 0.0, "n": 0.0, "counts": {}}
    counts = {b: seq.count(b) for b in "ACGTN"}
    gc = counts["G"] + counts["C"]
    acgt = counts["A"] + counts["C"] + counts["G"] + counts["T"]
    return {
        "length": n,
        "gc": round(100 * gc / acgt, 1) if acgt else 0.0,
        "n": round(100 * counts["N"] / n, 1),
        "counts": counts,
    }

def reverse_complement(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


_CODON = {}
for _i, _c1 in enumerate("TCAG"):
    for _j, _c2 in enumerate("TCAG"):
        for _k, _c3 in enumerate("TCAG"):
            _CODON[_c1 + _c2 + _c3] = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"[_i*16 + _j*4 + _k]


def translate(nt: str) -> str:
    """Standard genetic code; unknown/ambiguous codons -> 'X'."""
    nt = nt.upper()
    return "".join(_CODON.get(nt[i:i+3], "X") for i in range(0, len(nt) - 2, 3))

def sha256(seq: str) -> str:
    return hashlib.sha256(seq.encode()).hexdigest()

_CODON_STOP = {"TAA", "TAG", "TGA"}

def find_orfs(seq: str, min_aa: int = 40):
    """Real 6-frame ORF finder (ATG..stop). Returns list of dicts, 0-based half-open coords."""
    orfs = []
    rc = reverse_complement(seq)
    n = len(seq)
    for strand, s in (("+", seq), ("-", rc)):
        m = len(s)
        for frame in range(3):
            atg = None                                    # start of the currently-open ORF, if any
            j = frame
            while j < m - 2:
                cod = s[j:j+3]
                if atg is None:
                    if cod == "ATG":
                        atg = j
                elif cod in _CODON_STOP:                  # single left-to-right pass per frame: no per-ATG rescan (O(n), not O(n^2))
                    aa = (j - atg) // 3
                    if aa >= min_aa:
                        if strand == "+":
                            gs, ge = atg, j + 3
                        else:  # map back to forward coords
                            gs, ge = n - (j + 3), n - atg
                        orfs.append({"strand": strand, "frame": frame + 1,
                                     "start": gs, "end": ge,
                                     "length_nt": j + 3 - atg, "length_aa": aa, "open_end": False})
                    atg = None
                j += 3
            if atg is not None:                           # 3'-truncated ORF: reached the end with an open ATG and no stop.
                end = atg + ((m - atg) // 3) * 3           # keep it (flagged) so a degraded/truncated element still yields a scannable ORF
                aa = (end - atg) // 3
                if aa >= min_aa:
                    if strand == "+":
                        gs, ge = atg, end
                    else:
                        gs, ge = n - end, n - atg
                    orfs.append({"strand": strand, "frame": frame + 1,
                                 "start": gs, "end": ge,
                                 "length_nt": end - atg, "length_aa": aa, "open_end": True})
    orfs.sort(key=lambda o: (-o["length_aa"]))
    return orfs
