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
    "TEsorter": {
        "name": "Domain-architecture completeness (structural-completeness tier basis)",
        "citation": "Zhang RG, Li GY, Wang XL, et al. (2022) TEsorter: an accurate and fast method to classify LTR-retrotransposons in plant genomes. Horticulture Research 9:uhac017.",
        "doi": "10.1093/hr/uhac017", "license": "", "url": "https://doi.org/10.1093/hr/uhac017"},
    "LTRretriever": {
        "name": "Intact LTR-retrotransposon criteria (completeness tier basis)",
        "citation": "Ou S, Jiang N (2018) LTR_retriever: a highly accurate and sensitive program for identification of long terminal repeat retrotransposons. Plant Physiology 176(2):1410-1422.",
        "doi": "10.1104/pp.17.01310", "license": "", "url": "https://doi.org/10.1104/pp.17.01310"},
    "minimap2": {
        "name": "minimap2 (splice-aware alignment)",
        "citation": "Li H (2018) Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics 34(18):3094-3100.",
        "doi": "10.1093/bioinformatics/bty191", "license": "MIT", "url": "https://github.com/lh3/minimap2"},
    "miniprot": {
        "name": "miniprot (protein-to-genome spliced alignment)",
        "citation": "Li H (2023) Protein-to-genome alignment with miniprot. Bioinformatics 39(1):btad014.",
        "doi": "10.1093/bioinformatics/btad014", "license": "MIT", "url": "https://github.com/lh3/miniprot"},
    "isPcr": {                                     # UCSC In-Silico PCR; shares the BLAT tiling method (no separate isPcr paper)
        "name": "isPcr (UCSC In-Silico PCR)",
        "citation": "Kent WJ (2002) BLAT—the BLAST-like alignment tool. Genome Res 12(4):656-664.",
        "doi": "10.1101/gr.229202", "license": "Non-commercial (UCSC)", "url": "https://genome.ucsc.edu/cgi-bin/hgPcr"},
    "NCBI-RefSeq": {
        "name": "NCBI RefSeq genome assemblies (retrieved via NCBI Datasets)",
        "citation": "O'Leary NA, Wright MW, Brister JR, et al. (2016) Reference sequence (RefSeq) database at NCBI: current status, taxonomic expansion, and functional annotation. Nucleic Acids Res 44(D1):D733-D745.",
        "doi": "10.1093/nar/gkv1189", "license": "Public Domain", "url": "https://www.ncbi.nlm.nih.gov/refseq/"},
    # Primer secondary-structure QC (oligoqc). Every DOI below was independently confirmed against CrossRef +
    # PubMed by a 3-verifier pass. These are ADVISORY-method citations (attached to the result, never sealed —
    # the QC does not change the Primer3-designed primers, so it must not change the primer design seal).
    "SantaLucia1998": {
        "name": "Nearest-neighbor DNA thermodynamics (Tm / ΔG basis)",
        "citation": "SantaLucia J Jr (1998) A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics. Proc Natl Acad Sci USA 95(4):1460-1465.",
        "doi": "10.1073/pnas.95.4.1460", "license": "", "url": "https://doi.org/10.1073/pnas.95.4.1460"},
    "SantaLucia2004": {
        "name": "DNA structural-motif thermodynamic parameters",
        "citation": "SantaLucia J Jr, Hicks D (2004) The thermodynamics of DNA structural motifs. Annu Rev Biophys Biomol Struct 33:415-440.",
        "doi": "10.1146/annurev.biophys.32.110601.141800", "license": "", "url": "https://doi.org/10.1146/annurev.biophys.32.110601.141800"},
    "Owczarzy2008": {
        "name": "IDT OligoAnalyzer / SciTools methodology",
        "citation": "Owczarzy R, Tataurov AV, Wu Y, et al. (2008) IDT SciTools: a suite for analysis and design of nucleic acid oligomers. Nucleic Acids Res 36(Web Server issue):W163-W169.",
        "doi": "10.1093/nar/gkn198", "license": "", "url": "https://doi.org/10.1093/nar/gkn198"},
    "ViennaRNA": {
        "name": "ViennaRNA (independent secondary-structure cross-check; DNA params)",
        "citation": "Lorenz R, Bernhart SH, Höner zu Siederdissen C, et al. (2011) ViennaRNA Package 2.0. Algorithms Mol Biol 6:26.",
        "doi": "10.1186/1748-7188-6-26", "license": "custom (free for academic use)", "url": "https://www.tbi.univie.ac.at/RNA/"},
}


def oligoqc_refs():
    """References for the primer secondary-structure QC (hairpin/dimer/3'-end ΔG). Attached to a primer
    result as advisory-method provenance; NEVER folded into the sealed primer-design references, since the
    QC does not alter the Primer3-designed primers."""
    keys = ["SantaLucia1998", "SantaLucia2004", "Owczarzy2008", "ViennaRNA"]
    return [{"key": k, **REFS[k]} for k in keys]


def for_run(run_type: str, domains=None, fetched: bool = False):
    """Return the references actually used to produce this result."""
    keys = []
    if run_type == "analysis":
        keys += ["Wicker2007", "TEsorter", "LTRretriever"]    # classification framework + structural-completeness-tier basis
        if domains:
            keys += ["HMMER", "Pfam"]          # domain evidence only when domains were found
    if run_type == "primer":                   # Primer3 designs primers; in-silico PCR is TEagle's own pure-Python scanner,
        keys += ["Primer3"]                     # so its seal must not carry (and depend on) the primer3-py version
    if run_type == "annotate":
        keys += ["RepeatMasker", "Dfam", "Wicker2007"]
    if run_type == "splice":
        keys += ["minimap2"]
    if run_type == "homology":                 # protein-to-genome evidence; Dfam added when the bundled library is used
        keys += ["miniprot"]
    if run_type == "genome-scan":              # local whole-genome in-silico PCR against a downloaded RefSeq assembly
        keys += ["isPcr", "NCBI-RefSeq"]
    if fetched:
        keys += ["NCBI-Eutilities"]
    out, seen = [], set()
    for k in keys:
        if k not in seen:
            seen.add(k)
            out.append({"key": k, **REFS[k]})
    return out
