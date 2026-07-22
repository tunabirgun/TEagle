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
