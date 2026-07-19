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
