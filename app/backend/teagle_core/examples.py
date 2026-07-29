"""Bundled example elements — real, published transposable elements, not synthetic sequence.

The generated demo construct (`sample.py`) is a seeded random LTR element: its 210-codon ORF is random,
so it cannot hit a Pfam profile at E=1e-3 and a first-time user's first result is a structural-only null.
These five are the canonical specimens already pinned in the benchmark suite, chosen to exercise every
branch of `classify.py`: Copia and Gypsy (resolved by integrase-vs-RT order), a LINE, and both a
Tc1/Mariner and a hAT DNA transposon. Each carries its published identity so the classification can be
checked against the literature rather than trusted.
"""
from __future__ import annotations
import os

from . import appdirs

# (accession, label, organism, expected class — as published, for the user to check the call against)
CATALOG = [
    ("M11240", "copia", "Drosophila melanogaster",
     "LTR retrotransposon · Copia (Ty1) — integrase N-terminal to RT"),
    ("M12927", "gypsy (mdg4)", "Drosophila melanogaster",
     "LTR retrotransposon · Gypsy (Ty3) — integrase C-terminal to RT; env-bearing errantivirus"),
    ("M80343", "LINE-1 (L1.2)", "Homo sapiens",
     "Non-LTR retrotransposon · LINE — RT, no integrase, 3′ poly-A tail"),
    ("X01005", "Tc1", "Caenorhabditis elegans",
     "DNA transposon · Tc1/Mariner — transposase between terminal inverted repeats"),
    ("X05424", "Activator (Ac)", "Zea mays",
     "DNA transposon · hAT — the element McClintock described"),
]


def _dir() -> str:
    return appdirs.resource("teagle_core", "data", "examples") or \
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "examples")


def load(accession: str) -> str | None:
    """FASTA text for a bundled example, or None when it is not present in this build."""
    if not any(accession == a for a, *_ in CATALOG):          # never read an arbitrary path from the UI
        return None
    path = os.path.join(_dir(), accession + ".fasta")
    try:
        with open(path, encoding="utf-8") as fh:
            return fh.read()
    except OSError:
        return None


def available() -> list:
    """The catalog entries this build can actually serve, so the UI never offers a missing example."""
    return [e for e in CATALOG if load(e[0]) is not None]
