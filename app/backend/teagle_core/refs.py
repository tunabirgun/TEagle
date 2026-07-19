"""Verified database/tool references. Every DOI here was confirmed against a primary
source (CrossRef / the audited resource matrix, Deliverables/01_resource_matrix.md).
Citations are attached to results so a database call always carries its reference.
Never fabricate a citation — if a source is unverified it does not go here."""
from __future__ import annotations

try:
    import pyhmmer as _ph
    _PH = _ph.__version__
except Exception:
    _PH = "unknown"
try:
    import primer3 as _p3
    _P3 = _p3.__version__
except Exception:
    _P3 = "unknown"

# Note: version strings state only what is verifiable. The bundled Pfam-A profiles were
# fetched from InterPro on 2026-07-17 (exact Pfam release not reconciled — see matrix §3.2);
# HMMER runs via the engine bundled in pyhmmer; Primer3 via the primer3-py wheel.
REFS = {
    "HMMER": {
        "name": f"HMMER (via pyhmmer {_PH})",
        "citation": "Eddy SR (2011) Accelerated Profile HMM Searches. PLoS Comput Biol 7(10):e1002195.",
        "doi": "10.1371/journal.pcbi.1002195", "license": "BSD-3-Clause", "url": "http://hmmer.org"},
    "Pfam": {
        "name": "Pfam-A profiles (via InterPro, retrieved 2026-07-17)",
        "citation": "Mistry J, Chuguransky S, Williams L, et al. (2021) Pfam: The protein families database in 2021. Nucleic Acids Res 49(D1):D412-D419.",
        "doi": "10.1093/nar/gkaa913", "license": "CC0", "url": "https://www.ebi.ac.uk/interpro/"},
    "Primer3": {
        "name": f"Primer3 (via primer3-py {_P3})",
        "citation": "Untergasser A, Cutcutache I, Koressaar T, et al. (2012) Primer3-new capabilities and interfaces. Nucleic Acids Res 40(15):e115.",
        "doi": "10.1093/nar/gks596", "license": "GPL-2.0", "url": "https://primer3.org"},
    "NCBI-Eutilities": {
        "name": "NCBI E-utilities (Entrez)",
        "citation": "Sayers E. A General Introduction to the E-utilities. Entrez Programming Utilities Help. Bethesda (MD): NCBI.",
        "doi": "", "license": "Public Domain", "url": "https://www.ncbi.nlm.nih.gov/books/NBK25497/"},
    "Dfam": {
        "name": "Dfam 4.0 (curated)",
        "citation": "Storer J, Hubley R, Rosen J, Wheeler TJ, Smit AF (2021) The Dfam community resource of transposable element families, sequence models, and genome annotations. Mob DNA 12:2.",
        "doi": "10.1186/s13100-020-00230-y", "license": "CC0", "url": "https://www.dfam.org/"},
    "RepeatMasker": {
        "name": "RepeatMasker 4.2.4 (RMBLAST engine)",
        "citation": "Smit AFA, Hubley R, Green P. RepeatMasker Open-4.0. http://www.repeatmasker.org",
        "doi": "", "license": "OSL-2.1", "url": "https://www.repeatmasker.org/"},
    "Wicker2007": {
        "name": "TE classification framework",
        "citation": "Wicker T, Sabot F, Hua-Van A, et al. (2007) A unified classification system for eukaryotic transposable elements. Nat Rev Genet 8(12):973-982.",
        "doi": "10.1038/nrg2165", "license": "", "url": "https://doi.org/10.1038/nrg2165"},
    "minimap2": {
        "name": "minimap2 (splice-aware alignment)",
        "citation": "Li H (2018) Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics 34(18):3094-3100.",
        "doi": "10.1093/bioinformatics/bty191", "license": "MIT", "url": "https://github.com/lh3/minimap2"},
}


def for_run(run_type: str, domains=None, fetched: bool = False):
    """Return the references actually used to produce this result."""
    keys = []
    if run_type == "analysis":
        keys += ["Wicker2007"]                 # classification + structural-hallmark framework
        if domains:
            keys += ["HMMER", "Pfam"]          # domain evidence only when domains were found
    if run_type == "primer":                   # Primer3 designs primers; in-silico PCR is TEagle's own pure-Python scanner,
        keys += ["Primer3"]                     # so its seal must not carry (and depend on) the primer3-py version
    if run_type == "annotate":
        keys += ["RepeatMasker", "Dfam", "Wicker2007"]
    if run_type == "splice":
        keys += ["minimap2"]
    if fetched:
        keys += ["NCBI-Eutilities"]
    out, seen = [], set()
    for k in keys:
        if k not in seen:
            seen.add(k)
            out.append({"key": k, **REFS[k]})
    return out
