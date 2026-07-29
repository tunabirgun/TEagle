"""Capability benchmarks.
- Splicing / exon-intron detection: annotation-based (offline, @network) against genes with a published
  exon/intron count, plus de-novo minimap2 splice alignment (@wsl).
- Family-level naming: canonical TE specimens run through the Dfam/RepeatMasker pipeline (@wsl), each
  checked for a named family consistent with the element's structural class. Naming needs BOTH a lineage
  and the optional uncurated Dfam partitions. Three canonical elements that Dfam 4.0 simply does not
  contain are held separately, as an executable record of that coverage boundary.
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
    ("M12927", "Drosophila gypsy (LTR/Gypsy)",        "I"),
    ("X59545", "Drosophila mdg1 (LTR)",               "I"),
    ("M80343", "Human LINE-1 L1.2 (non-LTR)",         "I"),
    ("M17551", "Mouse IAP (LTR/ERV)",                 "I"),
    ("X01005", "C. elegans Tc1 (DNA/TcMar)",          "II"),
    ("M69216", "Drosophila hobo (DNA/hAT)",           "II"),
]

# Canonical elements that Dfam 4.0 CANNOT name, measured 2026-07-28 with the curated AND both uncurated
# consensus partitions installed and the organism supplied. Not a TEagle defect and not an installation
# gap — every consensus partition these taxa need was present; the families are simply absent from the
# database. Kept as an executable record of the coverage boundary rather than deleted, because a user
# who runs one of these needs to know the blank result is Dfam's limit and not their mistake.
#   maize Ac      -> 4 hits, all Low_complexity/Simple_repeat, no TE family
#   tobacco Tnt1  -> 0 hits
#   yeast Ty3     -> 0 hits  (S. cerevisiae has 32 families in the whole of Dfam)
# Dfam's coverage is strongest for metazoa and vertebrates; plant and fungal repeats are largely held in
# other resources. The structural and protein-domain layers, which do not consult Dfam, classify all
# three correctly.
DFAM_UNNAMEABLE = [
    ("X13777", "Tobacco Tnt1 (LTR/Copia)",            "I"),
    ("M23367", "Yeast Ty3 (LTR/Gypsy)",               "I"),
    ("X05424", "Maize Activator Ac (DNA/hAT)",        "II"),
]


@pytest.mark.wsl
@pytest.mark.parametrize("acc,desc,klass", FAMILY_BENCH)
def test_family_naming_benchmark(acc, desc, klass):
    """Each specimen must receive a NAMED Dfam family from the WSL backend (RepeatMasker + Dfam), and the
    call must not contradict the element's structural class.

    The organism comes from the fetched RECORD, not from a table written here: NCBI already states it
    authoritatively, and a hand-copied species is one more thing that can silently drift. It has to be
    supplied at all because RepeatMasker searches only a limited default set without a lineage —
    measured on Drosophila copia, which returns nothing but low-complexity with no species and resolves
    to Copia_LTR + Copia_I at 100% consensus coverage with one. Naming a family therefore needs BOTH a
    lineage and the uncurated Dfam partitions, which are an optional backend component.

    Requires the uncurated partitions; skipped rather than failed when only the curated set is present,
    since that is a legitimate installation and not a regression."""
    from teagle_core import wsl
    lib = (wsl.env_status().get("dfam_library") or {})
    if not any("uncurated" in str(p).lower() for p in lib.get("partitions") or []):
        pytest.skip("uncurated Dfam partitions not installed — most clade families cannot be named")
    meta = fetch.retrieve(acc)
    organism = meta.get("organism") or None
    r = wsl.annotate(meta["fasta"], species=organism, timeout=900)
    assert r.get("ok"), f"{acc} {desc} [{organism}]: {r.get('error')}"
    named = [h for h in r.get("hits", []) if h.get("family") and h["class_family"] not in
             {"Low_complexity", "Simple_repeat", "Satellite", "Unknown", "Unspecified"}]
    assert named, f"{acc} {desc} [{organism}]: no TE family named"
    # the call must not contradict the element's published class
    fams = " ".join(h["class_family"] for h in named)
    if klass == "I":
        assert "LTR" in fams or "LINE" in fams or "SINE" in fams or "Retro" in fams, (acc, desc, fams)
    else:
        assert "DNA" in fams or "RC" in fams, (acc, desc, fams)


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


# ---------------- HERV GAG-POL-ENV domain benchmark (v2.9.0; @network — fetch + offline pyhmmer) ----------------
# Verified full/partial HERV proviruses across families (efetch-fetched). A sound tool must (a) recover the full
# GAG-POL-ENV of an intact HML-2 provirus, and (b) NOT invent domains that are genuinely absent (HERV-L is env-less).
HERV_BENCH = [
    ("AY037928", "HERV-K113", {"GAG", "PR", "RT", "RNaseH", "INT", "ENV"}, True, None),   # full HML-2 -> intact ERV
    ("AF164615", "HERV-K109", {"GAG", "PR", "RT", "RNaseH", "INT", "ENV"}, True, None),   # full HML-2
    ("AJ289709", "HERV-H",    {"RT", "INT", "ENV"},                        True, None),   # env present, gag degraded
    ("X89211",   "HERV-L",    {"RT", "INT"},                               False, "ENV"), # env-less lineage: env must NOT appear
]


@pytest.mark.network
@pytest.mark.parametrize("acc,name,need,erv,forbid", HERV_BENCH)
def test_herv_domain_architecture_benchmark(acc, name, need, erv, forbid):
    """Each HERV specimen must show at least its expected domain modules; an intact HML-2 provirus is called an ERV;
    and a genuinely env-less lineage (HERV-L) must not be given a spurious env domain."""
    import sys as _s
    _n = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "app", "backend")
    if _n not in _s.path:
        _s.path.insert(0, _n)
    from teagle_core import structural, domains, classify
    seq = fetch.retrieve(acc)["fasta"]
    from teagle_core.sequtil import parse_fasta
    s = parse_fasta(seq)[0][1]
    dm = domains.scan_domains(s)
    codes = {d["domain"] for d in dm}
    assert need <= codes, f"{acc} {name}: expected {need}, got {codes}"
    if forbid:
        assert forbid not in codes, f"{acc} {name}: {forbid} should NOT be detected (genuinely absent)"
    cl = classify.classify(structural.detect_all(s), dm)
    assert cl.get("is_erv") is erv, f"{acc} {name}: is_erv expected {erv}"
    assert all("confidence" in d for d in dm)                # per-domain reliability (Axis 1) attached


@pytest.mark.wsl
@pytest.mark.parametrize("acc,desc,klass", DFAM_UNNAMEABLE)
def test_dfam_coverage_boundary_is_still_where_it_was(acc, desc, klass):
    """These canonical elements are NOT in Dfam 4.0. The test asserts the boundary, so that if a future
    Dfam release adds them the failure tells us to promote them into FAMILY_BENCH — and so that nobody
    re-investigates a blank result that has already been traced to the database rather than the tool.

    The structural and protein-domain layers do not consult Dfam and must classify them regardless."""
    from teagle_core import wsl, structural, domains, classify
    lib = (wsl.env_status().get("dfam_library") or {})
    if not any("uncurated" in str(p).lower() for p in lib.get("partitions") or []):
        pytest.skip("uncurated Dfam partitions not installed")
    meta = fetch.retrieve(acc)
    r = wsl.annotate(meta["fasta"], species=meta.get("organism"), timeout=900)
    assert r.get("ok"), f"{acc}: {r.get('error')}"
    named = [h for h in r.get("hits", []) if h.get("family") and h["class_family"] not in
             {"Low_complexity", "Simple_repeat", "Satellite", "Unknown", "Unspecified"}]
    assert not named, (f"{acc} {desc} is now nameable in Dfam ({[h['family'] for h in named]}) — move it "
                       f"into FAMILY_BENCH and update the coverage note")
    # the layers that do not depend on Dfam must still work
    seq = "".join(l for l in meta["fasta"].splitlines() if not l.startswith(">"))
    cl = classify.classify(structural.detect_all(seq), domains.scan_domains(seq))
    assert cl["te_class"] not in (None, "none"), f"{acc}: structural/domain classification also failed"
