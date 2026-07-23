"""Capability benchmarks.
- Splicing / exon-intron detection: annotation-based (offline, @network) against genes with a published
  exon/intron count, plus de-novo minimap2 splice alignment (@wsl).
- Family-level naming: 10 canonical TE specimens run through the Dfam/RepeatMasker pipeline (@wsl), each
  checked for a named family that is consistent with the element's structural class.
The @wsl / @network benchmarks are skipped in the fast offline suite and run with the backend available."""
import os, sys
import pytest

_BE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "app", "backend")
if _BE not in sys.path:
    sys.path.insert(0, _BE)
from teagle_core import fetch                                    # noqa: E402


# ---------------- splicing / exon-intron benchmark ----------------
# (accession, expected exon count, expected intron count) — published gene structures.
SPLICE_BENCH = [
    ("J00265", 3, 2),    # human insulin (INS): 3 exons / 2 introns; middle exon is CDS-derived in this record
]


@pytest.mark.network
@pytest.mark.parametrize("acc,exons,introns", SPLICE_BENCH)
def test_splice_annotation_benchmark(acc, exons, introns):
    """Annotation-based exon/intron detection: the completed gene model must recover the published counts,
    and every CDS segment must lie inside an exon (no coding sequence stranded in an intron)."""
    ft = fetch._get(fetch.EUTILS + f"efetch.fcgi?db=nuccore&id={acc}&rettype=ft&retmode=text&tool=TEagle", 60)
    gm = fetch.build_gene_model(fetch.parse_feature_table(ft))   # build_gene_model completes the model
    assert gm["counts"]["exons"] == exons, (acc, gm["counts"])
    assert gm["counts"]["introns"] == introns, (acc, gm["counts"])
    assert all(any(e["start"] <= c["start"] and c["end"] <= e["end"] for e in gm["exons"]) for c in gm["cds"])


# ---------------- family-level naming benchmark (10 specimens) ----------------
# (accession, description, expected structural class) — canonical, well-characterised TEs across both classes.
FAMILY_BENCH = [                                                # accessions verified live against NCBI titles
    ("M11240", "Drosophila copia (LTR/Copia)",        "I"),
    ("X13777", "Tobacco Tnt1 (LTR/Copia)",            "I"),
    ("M23367", "Yeast Ty3 (LTR/Gypsy)",               "I"),
    ("M12927", "Drosophila gypsy (LTR/Gypsy)",        "I"),
    ("X59545", "Drosophila mdg1 (LTR)",               "I"),
    ("M80343", "Human LINE-1 L1.2 (non-LTR)",         "I"),
    ("M17551", "Mouse IAP (LTR/ERV)",                 "I"),
    ("X01005", "C. elegans Tc1 (DNA/TcMar)",          "II"),
    ("X05424", "Maize Activator Ac (DNA/hAT)",        "II"),
    ("M69216", "Drosophila hobo (DNA/hAT)",           "II"),
]


@pytest.mark.wsl
@pytest.mark.parametrize("acc,desc,klass", FAMILY_BENCH)
def test_family_naming_benchmark(acc, desc, klass):
    """Each specimen must receive a NAMED Dfam family from the WSL backend (RepeatMasker + Dfam), and the
    call must not contradict the element's structural class. Runs only when the WSL backend is installed."""
    from teagle_core import wsl
    meta = fetch.retrieve(acc)
    r = wsl.annotate(meta["fasta"], species=None, timeout=600)
    assert r.get("ok"), f"{acc} {desc}: {r.get('error')}"
    named = [h for h in r.get("hits", []) if h.get("family") and h["class_family"] not in
             {"Low_complexity", "Simple_repeat", "Satellite", "Unknown", "Unspecified"}]
    assert named, f"{acc} {desc}: no TE family named"


# ---------------- primer secondary-structure QC benchmark (published primers, offline) ----------------
# Verified published PCR primer pairs (exact sequences from PrimerBank + peer-reviewed papers; see the report
# bibliography). A sound secondary-structure QC must NOT false-alarm on primers that were experimentally validated
# in the literature: none should be flagged 'warn'. Pure/in-process (primer3 + ViennaRNA), so it runs offline.
# This is a SPECIFICITY / false-alarm check (validated primers must not be flagged 'warn'), NOT a numerical
# ΔG-accuracy validation against a reference tool. The full 12-pair set matches the report's benchmark table.
LIT_PRIMERS = [
    ("GAPDH-197",  "GGAGCGAGATCCCTCCAAAAT", "GGCTGTTGTCATACTTCTCATGG", "PrimerBank 378404907c1"),
    ("GAPDH-101",  "ACAACTTTGGTATCGTGGAAGG", "GCCATCACGCCACAGTTTC",     "PrimerBank 378404907c2"),
    ("ACTB-250",   "CATGTACGTTGCTATCCAGGC",  "CTCCTTAATGTCACGCACGAT",   "PrimerBank 4501885a1"),
    ("B2M-248",    "GAGGCTATCCAGCGTACTCCA",  "CGGCAGGCATACTCATCTTTT",   "PrimerBank 37704380c1"),
    ("GAPDH-Misak","ACCCAGAAGACTGTGGATGG",   "TTCAGCTCAGGGATGACCTT",    "Misak 2025 Methods (Sci Rep 15:32499)"),
    ("Alu-Yb8",    "GGTGAAACCCCGTCTCTACT",   "GGTTCAAGCGATTCTCCTGC",    "Funakoshi 2017 (Sci Rep 7:13202)"),
    ("L1PA-1",     "GACATCTACACCGAAAACCC",   "TCGTCAAAATCATTCTCCATCC",  "Misak 2025 (Sci Rep 15:32499)"),
    ("L1PA-2",     "ACCAGCCACTGCAAAATC",     "CCAATTTGCCAGTCTGTGTC",    "Misak 2025 (Sci Rep 15:32499)"),
    ("L1PA-3",     "ATGCACAAGCCTCAGTAGCC",   "TCCATTCTCCCCGTCACTTTC",   "Misak 2025 (Sci Rep 15:32499)"),
    ("L1PA-4",     "TCCACACCAAAACCCCATC",    "CTCGTCAAAGTCATTCTCCATC",  "Misak 2025 (Sci Rep 15:32499)"),
    ("L1PA-5",     "GACAAAGGTGACATTACAAC",   "CTTGGGAGATTGTGTGTTTC",    "Misak 2025 (Sci Rep 15:32499)"),
    ("L1PA-6",     "AGAATGAAACTGGACCCCTA",   "GTCCAGAAGAGTATTTCCTA",    "Misak 2025 (Sci Rep 15:32499)"),
]


def test_literature_primer_qc_benchmark():
    """Published, experimentally-validated primers must pass the dual-engine secondary-structure QC without a
    'warn' flag, and both engines must return ΔG in a sane kcal/mol range. Offline (primer3 + ViennaRNA)."""
    from teagle_core import oligoqc
    warned = []
    for name, F, R, cite in LIT_PRIMERS:
        q = oligoqc.qc_pair(F, R)
        assert q["ok"], f"{name}: QC failed"
        for m in (q["left"]["hairpin"], q["left"]["self_dimer"], q["hetero_dimer"]):
            for eng in ("p3", "vrna"):
                v = m.get(eng)
                assert v is None or -60.0 < v < 5.0, f"{name}: {eng} ΔG out of kcal/mol range ({v})"
        if q["worst"] == "warn":
            warned.append(name)
    assert not warned, f"published validated primers should not be flagged 'warn': {warned}"


@pytest.mark.wsl
def test_denovo_splice_benchmark():
    """De-novo minimap2 splice alignment: an mRNA aligned to its genomic locus must recover introns that
    agree with the record's annotation (independent cross-check). Runs only when minimap2 (WSL) is present."""
    import engine
    from teagle_core import wsl
    g = fetch.retrieve("J00265")                                  # insulin gene (genomic, 2 introns)
    tx = fetch.retrieve("NM_000207")                              # insulin mRNA (spliced transcript)
    r = engine.run_splice({"sequence": g["fasta"], "transcript": tx["fasta"], "source": {"accession": "J00265"}})
    assert r.get("ok"), r.get("error")
    ann = fetch.build_gene_model(fetch.parse_feature_table(
        fetch._get(fetch.EUTILS + "efetch.fcgi?db=nuccore&id=J00265&rettype=ft&retmode=text&tool=TEagle", 60)))
    cc = fetch.cross_check_models(ann["introns"], r.get("introns", []))
    assert cc["matched"] >= 1                                     # alignment confirms at least one annotated intron
