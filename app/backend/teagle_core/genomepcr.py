"""Local whole-genome in-silico PCR — parse UCSC isPcr output into amplicons.

The scan itself runs in the WSL backend (wsl.genome_scan) with isPcr against a downloaded,
checksummed RefSeq assembly — no remote query (a local safety timeout still applies). isPcr applies the 3' perfect-match
priming rule and amplicon assembly natively (it is the engine behind UCSC In-Silico PCR), so this
module only has to parse its FASTA output and shape amplicon records. Kept pure (no WSL, no network)
so the parser is unit-testable in isolation. Because the assembly is a fixed local file (accession +
sha256), a genome scan is reproducible and IS sealed — unlike the retired remote path.
"""
from __future__ import annotations
import re
from collections import Counter

from .fetch import FetchError


class GenomePcrError(FetchError):
    """A whole-genome-scan failure surfaced to the user (never a 500)."""


# isPcr FASTA header:  >CHROM:START(+|-)END  name  SIZEbp  FWD  REV
# START/END are 1-based inclusive; the sign between them is the subject strand of the product.
_HDR = re.compile(r"^>(\S+):(\d+)([+-])(\d+)\s+(\S+)\s+(\d+)bp\s+(\S+)\s+(\S+)\s*$")


def parse_ispcr(fasta_text: str) -> list:
    """Parse isPcr FASTA output into amplicon dicts. Query rows are named 'pair' (fwd+rev),
    'fwdonly' (fwd+fwd) and 'revonly' (rev+rev); the last two are single-primer products. Every
    genome-wide product is off-target by definition (a TE primer pair has no one designed locus)."""
    amps = []
    for line in fasta_text.splitlines():
        line = line.rstrip()
        if not line.startswith(">"):
            continue
        m = _HDR.match(line)
        if not m:
            continue
        chrom, s, sign, e, name, size, fwd, rev = m.groups()
        start, end = int(s), int(e)
        lo, hi = (start, end) if end >= start else (end, start)
        amps.append({
            "source": chrom, "start": lo, "end": hi, "length": int(size),
            "strand": "+" if sign == "+" else "-",
            "fwd_primer": fwd, "rev_primer": rev,
            "single_primer": name in ("fwdonly", "revonly") or fwd.upper() == rev.upper(),
            "pair": name, "on_target": False,
        })
    # stable order: fewest surprises first — by chromosome then position
    amps.sort(key=lambda a: (a["single_primer"], a["source"], a["start"]))
    return amps


def summarize(amps: list, has_locus: bool = False) -> dict:
    """Interpret a whole-genome isPcr scan of ONE primer pair. 'pair' products (fwd+rev) are the real
    amplicons; fwdonly/revonly rows are single-primer artefacts, counted separately.

    On/off-target is defined relative to an intended locus. When the primers were designed on a locus that
    sits in the scanned assembly (has_locus=True), the product overlapping that locus is the ON-TARGET
    (marked on_target=True upstream) and the rest are OFF-target paralogs; the verdict is then a specificity
    call — how many off-target sites accompany the intended one. When the primers are a bare consensus with
    no genome position (has_locus=False), there is no single intended locus, so a genome-wide product is
    neither on- nor off-target — it is a neutral GENOMIC PRIMING SITE (candidate), and the verdict reads
    family-generic vs locus-specific. isPcr's perfect-3' rule makes the count a conservative floor. Pure."""
    pair = [a for a in amps if not a.get("single_primer")]
    single = [a for a in amps if a.get("single_primer")]
    on = [a for a in pair if a.get("on_target")]
    off = [a for a in pair if not a.get("on_target")]
    src_counts = Counter(a["source"] for a in pair)
    per_source = sorted(src_counts.items(), key=lambda kv: (-kv[1], kv[0]))   # busiest sequence first
    lengths = [a["length"] for a in pair]
    size_mode = size_mode_n = size_min = size_max = None
    if lengths:
        size_mode, size_mode_n = Counter(lengths).most_common(1)[0]
        size_min, size_max = min(lengths), max(lengths)
    n_pair, n_on, n_off, n_sources = len(pair), len(on), len(off), len(src_counts)
    if has_locus:
        # n_on is 0 or 1 (exactly one best-overlapping product is the on-target, chosen upstream)
        if n_on == 0 and n_off == 0:
            tier, verdict = "none", "no genome-wide product — the primers do not amplify anywhere in this assembly"
        elif n_on == 0:
            tier, verdict = "off-target-only", (f"{n_off} off-target site(s) and NO product at the intended locus — "
                                                "the pair amplifies elsewhere, not the designed target")
        elif n_off == 0:
            tier, verdict = "specific", "copy-specific — 1 on-target product, no off-target site in this assembly"
        elif n_off <= 5:
            tier, verdict = "low-off-target", f"1 on-target + {n_off} off-target site(s) — low-copy / paralogous"
        else:
            tier, verdict = "family-generic", (f"1 on-target + {n_off} off-target site(s) across {n_sources} "
                                               "sequence(s) — family-generic (expected for a TE-consensus pair)")
    else:
        if n_pair == 0:
            tier, verdict = "none", "no genome-wide priming site (in this assembly, at the current 3' match length)"
        elif n_pair == 1:
            tier, verdict = "locus-specific", "locus-specific in this assembly (1 genomic priming site)"
        elif n_pair <= 5:
            tier, verdict = "low-copy", f"low-copy / paralogous ({n_pair} genomic priming sites)"
        else:
            tier, verdict = "family-generic", (f"family-generic — {n_pair} genomic priming sites across {n_sources} "
                                               "sequence(s) (expected for a TE-consensus pair)")
    return {"n_total": len(amps), "n_pair": n_pair, "n_on": n_on, "n_off": n_off, "n_single": len(single),
            "n_sources": n_sources, "per_source": per_source, "has_locus": has_locus,
            "size_mode": size_mode, "size_mode_n": size_mode_n, "size_min": size_min, "size_max": size_max,
            "tier": tier, "verdict": verdict}


def query_rows(fwd: str, rev: str) -> str:
    """isPcr query file: the pair plus both single-primer combinations (F+F, R+R), one per line.
    Names are self-generated (never user text), so the file carries only validated primer bases."""
    fwd, rev = fwd.upper().strip(), rev.upper().strip()
    return f"pair\t{fwd}\t{rev}\nfwdonly\t{fwd}\t{fwd}\nrevonly\t{rev}\t{rev}\n"
