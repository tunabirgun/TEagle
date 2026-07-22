"""Local whole-genome in-silico PCR — parse UCSC isPcr output into amplicons.

The scan itself runs in the WSL backend (wsl.genome_scan) with isPcr against a downloaded,
checksummed RefSeq assembly — no remote query, no timeouts. isPcr applies the 3' perfect-match
priming rule and amplicon assembly natively (it is the engine behind UCSC In-Silico PCR), so this
module only has to parse its FASTA output and shape amplicon records. Kept pure (no WSL, no network)
so the parser is unit-testable in isolation. Because the assembly is a fixed local file (accession +
sha256), a genome scan is reproducible and IS sealed — unlike the retired remote path.
"""
from __future__ import annotations
import re

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


def query_rows(fwd: str, rev: str) -> str:
    """isPcr query file: the pair plus both single-primer combinations (F+F, R+R), one per line.
    Names are self-generated (never user text), so the file carries only validated primer bases."""
    fwd, rev = fwd.upper().strip(), rev.upper().strip()
    return f"pair\t{fwd}\t{rev}\nfwdonly\t{fwd}\t{fwd}\nrevonly\t{rev}\t{rev}\n"
