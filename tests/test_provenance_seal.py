"""Regression: the content-addressed provenance seal must be invariant to volatile inputs —
the fetch wall-clock timestamp and the version/availability of tools a run does not use —
and DB provenance must be recorded, not hardcoded empty."""
from teagle_core import provenance


def test_seal_invariant_to_retrieved_utc():
    a_src = {"accession": "M11240", "organism": "Drosophila", "taxid": "7227",
             "source": "NCBI nuccore", "endpoint": "efetch", "retrievedUtc": "2026-07-18T10:00:00+00:00"}
    b_src = dict(a_src, retrievedUtc="2026-07-18T23:59:59+00:00")
    a = provenance.build_manifest("analysis", "ACGTACGTACGT", "M11240", {"orf_min_aa": 40}, source=a_src)
    b = provenance.build_manifest("analysis", "ACGTACGTACGT", "M11240", {"orf_min_aa": 40}, source=b_src)
    assert a["input"]["retrievedUtc"] != b["input"]["retrievedUtc"]     # recorded verbatim...
    assert a["manifestSha256"] == b["manifestSha256"]                   # ...but never sealed


def test_analysis_seal_excludes_primer3():
    m = provenance.build_manifest("analysis", "ACGT", "x", {})
    names = " ".join(s["name"].lower() for s in m["software"])
    assert "primer3" not in names                                       # unused tool must not enter the seal
    assert "pyhmmer" in names or "hmmer" in names


def test_primer_seal_includes_primer3():
    names = " ".join(s["name"].lower() for s in provenance.build_manifest("primer", "ACGT", "x", {})["software"])
    assert "primer3" in names


def test_seal_deterministic_same_input():
    a = provenance.build_manifest("primer", "ACGTACGT", "x", {"opt_size": 20})
    b = provenance.build_manifest("primer", "ACGTACGT", "x", {"opt_size": 20})
    assert a["manifestSha256"] == b["manifestSha256"]


def test_databases_recorded_when_passed():
    dbs = [{"name": "Pfam TE-domain profiles", "sha256": "deadbeef"}]
    assert provenance.build_manifest("analysis", "ACGT", "x", {}, databases=dbs)["databases"] == dbs


def test_annotate_manifest_not_self_contradictory():
    m = provenance.build_manifest("annotate", "ACGT", "x", {"engine": "RMBLAST"},
                                  not_run=["External NCBI Primer-BLAST"],
                                  databases=[{"name": "Dfam (curated)", "version": "4.0"}])
    assert m["databases"], "an annotate run must record the DB it used"
    joined = " ".join(m["notRun"]).lower()
    assert "dfam" not in joined and "repeatmasker" not in joined        # cannot claim the step it ran did not run


# ---------- coordinate-fetch provenance: identity sealed, display labels excluded ----------
def _coord_src(**over):
    src = {"accession": "NC_000013.11", "organism": "Homo sapiens", "taxid": "9606",
           "assemblyAccession": "GCF_000001405.40",
           "regions": [{"chrAccession": "NC_000013.11", "start": 33016423, "stop": 33066143, "strand": 1}],
           "coordSystem": "1-based-inclusive", "retrievalType": "coordinate",
           "assemblyName": "GRCh38.p14", "displayLocus": "chr13:33,016,423-33,066,143", "chromName": "chr13",
           "source": "NCBI E-utilities", "endpoint": "efetch", "retrievedUtc": "2026-07-21T10:00:00+00:00"}
    src.update(over)
    return src


def _seal(src):
    return provenance.build_manifest("analysis", "ACGTACGT", "coord", {"orf_min_aa": 40}, source=src)["manifestSha256"]


def test_coord_seal_records_identity_and_labels():
    m = provenance.build_manifest("analysis", "ACGTACGT", "coord", {}, source=_coord_src())
    inp = m["input"]
    for k in ("assemblyAccession", "regions", "coordSystem", "retrievalType",   # sealed identity...
              "assemblyName", "displayLocus", "chromName"):                      # ...and recorded labels
        assert k in inp, k                                                       # everything is recorded verbatim


def test_coord_seal_invariant_to_display_labels():
    # cosmetic labels (assembly display name, human locus string, chrom label) must NOT change the seal
    base = _seal(_coord_src())
    assert base == _seal(_coord_src(assemblyName="Genome Reference Consortium h38"))
    assert base == _seal(_coord_src(displayLocus="chr13 : 33016423 - 33066143"))
    assert base == _seal(_coord_src(chromName="13"))
    assert base == _seal(_coord_src(retrievedUtc="2026-01-01T00:00:00+00:00"))   # wall-clock never sealed


def test_coord_seal_changes_with_identity():
    # the pinned assembly, taxid, coordinates, and strand are the reproducible identity — each flips the seal
    base = _seal(_coord_src())
    assert base != _seal(_coord_src(assemblyAccession="GCF_009914755.1"))
    assert base != _seal(_coord_src(taxid="10090"))          # organism pin is sealed (curated must carry the real taxid)
    moved = [{"chrAccession": "NC_000013.11", "start": 33016424, "stop": 33066143, "strand": 1}]
    assert base != _seal(_coord_src(regions=moved))
    minus = [{"chrAccession": "NC_000013.11", "start": 33016423, "stop": 33066143, "strand": 2}]
    assert base != _seal(_coord_src(regions=minus))


def test_seal_invariant_to_serving_database():
    # the DB that served the FASTA (NCBI vs the ENA fallback) is recorded but must NOT change the seal
    a = {"accession": "M11240.1", "organism": "Drosophila", "taxid": "7227", "source": "NCBI nuccore",
         "endpoint": "https://eutils.ncbi.nlm.nih.gov/.../efetch.fcgi", "sourceUrl": "https://www.ncbi.nlm.nih.gov/nuccore/M11240.1",
         "retrievedUtc": "2026-07-21T10:00:00+00:00"}
    b = dict(a, source="ENA (EMBL-EBI)", endpoint="https://www.ebi.ac.uk/ena/browser/api/fasta/M11240.1",
             sourceUrl="https://www.ebi.ac.uk/ena/browser/view/M11240.1")
    ma = provenance.build_manifest("analysis", "ACGTACGT", "M11240.1", {}, source=a)
    mb = provenance.build_manifest("analysis", "ACGTACGT", "M11240.1", {}, source=b)
    assert ma["input"]["source"] != mb["input"]["source"]              # recorded verbatim...
    assert ma["manifestSha256"] == mb["manifestSha256"]                # ...never sealed (refetch-invariant)


def test_splice_seal_includes_transcript():
    import hashlib
    g = "ACGTACGTACGTACGT"
    m1 = provenance.build_manifest("splice", g, "gid", {"tool": "minimap2",
         "transcript_sha256": hashlib.sha256(b"transcript-one").hexdigest()})
    m2 = provenance.build_manifest("splice", g, "gid", {"tool": "minimap2",
         "transcript_sha256": hashlib.sha256(b"transcript-two").hexdigest()})
    assert m1["manifestSha256"] != m2["manifestSha256"]                # same genomic, different transcript -> different seal


def test_coord_and_accession_seals_are_distinct():
    # a coordinate run and a bare-accession run of the same bytes must not collide
    acc_src = {"accession": "NC_000013.11", "organism": "Homo sapiens", "taxid": "9606",
               "source": "NCBI nuccore", "endpoint": "efetch", "retrievedUtc": "2026-07-21T10:00:00+00:00"}
    assert _seal(_coord_src()) != _seal(acc_src)
