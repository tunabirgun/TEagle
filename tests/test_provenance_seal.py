"""Regression: the content-addressed provenance seal must be invariant to volatile inputs —
the fetch wall-clock timestamp and the version/availability of tools a run does not use —
and DB provenance must be recorded, not hardcoded empty."""
import copy
import hashlib
import json
import pytest

from teagle_core import provenance
from teagle_core.gff3 import public_record_identity


_NCBI_ENDPOINT = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def _validated_accession_source(accession="M11240.1", sequence="ACGT", *, ena=False, retrieved=None):
    label = "ENA (EMBL-EBI)" if ena else "NCBI nuccore"
    endpoint = (f"https://www.ebi.ac.uk/ena/browser/api/fasta/{accession}" if ena else _NCBI_ENDPOINT)
    source_url = (f"https://www.ebi.ac.uk/ena/browser/view/{accession}" if ena else
                  f"https://www.ncbi.nlm.nih.gov/nuccore/{accession}")
    raw = {"accession": accession, "taxid": "7227", "source": label, "endpoint": endpoint,
           "retrievedUtc": retrieved or "2026-07-18T10:00:00+00:00", "sourceUrl": source_url,
           "length": len(sequence)}
    return provenance.normalize_source_context(raw, [(accession, sequence)])


def test_accession_source_wire_json_round_trip_revalidates_against_actual_record():
    sequence = "ACGTACGT"
    source = _validated_accession_source("M11240.1", sequence)
    wire = json.loads(json.dumps(source))
    assert provenance.validate_source_context(wire, [("M11240.1", sequence)]) == source
    assert set(wire) == {"sourceContextVersion", "inputType", "source", "endpoint", "accession",
                         "retrievedUtc", "sourceUrl", "records", "taxid"}
    assert wire["records"] == [{"id": "M11240.1", "length": len(sequence),
                                 "sha256": hashlib.sha256(sequence.encode()).hexdigest()}]


def test_raw_fetch_response_is_accepted_only_by_normalizer_boundary():
    sequence = "ACGTACGT"
    raw = {
        "accession": "M11240.1", "taxid": "7227", "source": "NCBI nuccore",
        "endpoint": _NCBI_ENDPOINT, "retrievedUtc": "2026-07-18T10:00:00+00:00",
        "sourceUrl": "https://www.ncbi.nlm.nih.gov/nuccore/M11240.1", "length": len(sequence),
    }
    records = [("M11240.1", sequence)]
    assert provenance.normalize_source_context(raw, records)["inputType"] == "fetched_accession"
    assert provenance.validate_source_context(raw, records) == {
        "sourceContextVersion": "1.0", "inputType": "local_input"}


@pytest.mark.parametrize("mutate", [
    lambda d: d.update(provider="ncbi"),
    lambda d: d.update(endpoint="efetch"),
    lambda d: d.update(accession="M11240.2"),
    lambda d: d.update(retrievedUtc="2026-07-18T10:00:00.123+00:00"),
    lambda d: d.update(sourceUrl="https://www.ncbi.nlm.nih.gov.attacker.example/nuccore/M11240.1"),
    lambda d: d.update(taxid=7227),
    lambda d: d.update(identitySource="validated_fetch"),
    lambda d: d["records"][0].update(id="M11240.2"),
    lambda d: d["records"][0].update(length=9),
    lambda d: d["records"][0].update(length=True),
    lambda d: d["records"][0].update(sha256="0" * 64),
])
def test_accession_source_wire_forgery_downgrades_atomically_to_local(mutate):
    sequence = "ACGTACGT"
    source = _validated_accession_source("M11240.1", sequence)
    mutate(source)
    assert provenance.validate_source_context(source, [("M11240.1", sequence)]) == {
        "sourceContextVersion": "1.0", "inputType": "local_input"}


def test_seal_invariant_to_retrieved_utc():
    sequence = "ACGTACGTACGT"
    a_src = _validated_accession_source("M11240", sequence, retrieved="2026-07-18T10:00:00+00:00")
    b_src = _validated_accession_source("M11240", sequence, retrieved="2026-07-18T23:59:59+00:00")
    identity = {**public_record_identity("M11240", trust="validated_fetch"),
                "identitySource": "validated_fetch"}
    a = provenance.build_manifest("analysis", sequence, "M11240", {"orf_min_aa": 40}, source=a_src,
                                  input_identity=identity)
    b = provenance.build_manifest("analysis", sequence, "M11240", {"orf_min_aa": 40}, source=b_src,
                                  input_identity=identity)
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


def test_seal_contract_is_versioned_utf8_and_verifiable():
    m = provenance.build_manifest("analysis", "ACGT", "é–record", {"note": "not run — unavailable"})
    assert m["schemaVersion"] == "0.2"
    assert m["canonicalization"] == "teagle-json-utf8-v1"
    assert provenance.verify_manifest(m)
    encoded = provenance.canonical_manifest_bytes(m)
    assert "not run — unavailable".encode("utf-8") in encoded
    assert b"\\u00e9" not in encoded and b"\\u2013" not in encoded
    assert provenance.verify_manifest(json.loads(json.dumps(m, ensure_ascii=False)))
    encoded_document = provenance.manifest_json_bytes(m)
    assert provenance.verify_manifest(json.loads(encoded_document))


def test_manifest_strict_json_rejects_nonfinite_or_tampered_input():
    manifest = provenance.build_manifest("analysis", "ACGT", "x", {"threshold": 1})
    tampered = copy.deepcopy(manifest)
    tampered["parameters"]["threshold"] = 2
    with pytest.raises(ValueError, match="seal verification"):
        provenance.manifest_json_bytes(tampered)
    nonfinite = provenance.build_manifest("analysis", "ACGT", "x", {"threshold": 1})
    nonfinite["environment"]["probe"] = float("nan")
    with pytest.raises(ValueError):
        provenance.manifest_json_bytes(nonfinite)


def test_canonical_utf8_golden_vector():
    probe = {
        "schemaVersion": "0.2", "canonicalization": "teagle-json-utf8-v1",
        "createdUtc": "volatile", "environment": {"os": "volatile"},
        "input": {"id": "record_aaaaaaaaaaaa", "idSha256": "a" * 64,
                  "identitySource": "local_input", "sourceContextVersion": "1.0",
                  "inputType": "local_input", "length": 4,
                  "sha256": "1" * 64,
                  "records": [{"id": "record_bbbbbbbbbbbb", "idSha256": "b" * 64,
                               "identitySource": "local_input", "length": 4, "sha256": "0" * 64}]},
        "parameters": {"note": "not run — unavailable"}, "manifestSha256": "0" * 64,
    }
    expected = (b'{"canonicalization":"teagle-json-utf8-v1","input":{"idSha256":"' + b"a" * 64 +
                b'","identitySource":"local_input","inputType":"local_input","length":4,"records":'
                b'[{"id":"record_bbbbbbbbbbbb","idSha256":"' + b"b" * 64 +
                b'","identitySource":"local_input","length":4,"sha256":"' + b"0" * 64 +
                b'"}],"sha256":"' + b"1" * 64 + b'","sourceContextVersion":"1.0"},"parameters":{"note":"not run '
                b'\xe2\x80\x94 unavailable"},"schemaVersion":"0.2"}')
    assert provenance.canonical_manifest_bytes(probe) == expected
    probe["manifestSha256"] = hashlib.sha256(expected).hexdigest()
    assert provenance.verify_manifest(probe)


def test_seal_verification_rejects_every_sealed_mutation():
    m = provenance.build_manifest("analysis", "ACGT", "x", {"threshold": 1})
    for path, value in [
        (("input", "sha256"), "0" * 64),
        (("parameters", "threshold"), 2),
        (("notRun",), ["different"]),
        (("canonicalization",), "unknown"),
    ]:
        changed = copy.deepcopy(m)
        target = changed
        for key in path[:-1]:
            target = target[key]
        target[path[-1]] = value
        assert not provenance.verify_manifest(changed)


def test_ordered_record_boundaries_and_ids_are_sealed():
    one = [{"id": "whole", "sequence": "ACGT" * 6}]
    two = [{"id": "left", "sequence": "ACGT" * 3},
           {"id": "right", "sequence": "ACGT" * 3}]
    a = provenance.build_manifest("analysis", "ACGT" * 6, "whole", {}, input_records=one)
    b = provenance.build_manifest("analysis", "ACGT" * 6, "left", {}, input_records=two)
    assert a["input"]["sha256"] == b["input"]["sha256"]
    assert a["manifestSha256"] != b["manifestSha256"]
    assert [r["id"] for r in b["input"]["records"]] == [
        provenance.public_record_identity("left")["id"],
        provenance.public_record_identity("right")["id"],
    ]
    assert all("idSha256" in r for r in b["input"]["records"])
    reversed_records = copy.deepcopy(b)
    reversed_records["input"]["records"].reverse()
    assert not provenance.verify_manifest(reversed_records)


def test_analysis_result_descriptor_is_inside_manifest_seal():
    descriptor = {"resultSchemaVersion": "1.0", "resultPayloadSha256": "a" * 64}
    manifest = provenance.build_manifest("analysis", "ACGT", "x", {}, result_descriptor=descriptor)
    assert manifest["resultDescriptor"] == descriptor and provenance.verify_manifest(manifest)
    manifest["resultDescriptor"]["resultPayloadSha256"] = "b" * 64
    assert not provenance.verify_manifest(manifest)


@pytest.mark.parametrize("descriptor", [
    {},
    {"resultSchemaVersion": "1.0", "resultPayloadSha256": "a" * 63},
    {"resultSchemaVersion": "", "resultPayloadSha256": "a" * 64},
    {"resultSchemaVersion": "1.0", "resultPayloadSha256": "g" * 64},
    {"resultSchemaVersion": "1.0", "resultPayloadSha256": "a" * 64, "extra": True},
])
def test_invalid_analysis_result_descriptor_is_rejected(descriptor):
    with pytest.raises(ValueError):
        provenance.build_manifest("analysis", "ACGT", "x", {}, result_descriptor=descriptor)


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
           "regions": [{"chrAccession": "NC_000013.11", "start": 1, "stop": 8, "strand": 1}],
           "coordSystem": "1-based-inclusive", "retrievalType": "coordinate",
           "assemblyName": "GRCh38.p14", "displayLocus": "chr13:33,016,423-33,066,143", "chromName": "chr13",
           "source": "NCBI E-utilities (efetch seq_start/seq_stop)", "endpoint": _NCBI_ENDPOINT,
           "retrievedUtc": "2026-07-21T10:00:00+00:00"}
    src.update(over)
    return src


def _validated_coord(src):
    region = src["regions"][0]
    header = f"{region['chrAccession']}:{region['start']}-{region['stop']}" + ("(-)" if region["strand"] == 2 else "")
    return provenance.normalize_source_context(src, [(header, "ACGTACGT")])


def _seal(src):
    validated = _validated_coord(src)
    region = validated["regions"][0]
    rid = f"{region['chrAccession']}:{region['start']}-{region['end']}"
    identity = {**public_record_identity(rid, trust="validated_fetch"), "identitySource": "validated_fetch"}
    return provenance.build_manifest("analysis", "ACGTACGT", rid, {"orf_min_aa": 40}, source=validated,
                                     input_identity=identity)["manifestSha256"]


def test_coord_seal_records_identity_and_labels():
    source = _validated_coord(_coord_src())
    identity = {**public_record_identity("NC_000013.11:1-8", trust="validated_fetch"),
                "identitySource": "validated_fetch"}
    m = provenance.build_manifest("analysis", "ACGTACGT", identity["id"], {}, source=source,
                                  input_identity=identity)
    inp = m["input"]
    for k in ("assemblyAccession", "regions", "coordSystem", "source", "endpoint", "sourceUrl"):
        assert k in inp, k
    for omitted in ("organism", "assemblyName", "displayLocus", "chromName", "retrievalType"):
        assert omitted not in inp


def test_coordinate_source_exports_only_exact_bound_region_schema():
    raw = _coord_src(organism="patient_042", assemblyName="donor-A1", endpointPath=r"C:\private\endpoint")
    raw["regions"][0]["unknownPath"] = "/private/region"
    validated = _validated_coord(raw)
    assert validated["inputType"] == "fetched_coordinate"
    assert set(validated) == {"sourceContextVersion", "inputType", "source", "endpoint", "accession", "assemblyAccession",
                              "retrievedUtc", "coordSystem", "regions", "sourceUrl", "records", "taxid"}
    assert set(validated["regions"][0]) == {"chrAccession", "start", "end", "strand"}
    assert all(value not in json.dumps(validated) for value in
               ("patient_042", "donor-A1", r"C:\private\endpoint", "/private/region"))


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
    moved = [{"chrAccession": "NC_000013.11", "start": 2, "stop": 9, "strand": 1}]
    assert base != _seal(_coord_src(regions=moved))
    minus = [{"chrAccession": "NC_000013.11", "start": 1, "stop": 8, "strand": 2}]
    assert base != _seal(_coord_src(regions=minus))


def test_seal_invariant_to_serving_database():
    # the DB that served the FASTA (NCBI vs the ENA fallback) is recorded but must NOT change the seal
    sequence = "ACGTACGT"
    a = _validated_accession_source("M11240.1", sequence)
    b = _validated_accession_source("M11240.1", sequence, ena=True)
    identity = {**public_record_identity("M11240.1", trust="validated_fetch"),
                "identitySource": "validated_fetch"}
    ma = provenance.build_manifest("analysis", sequence, "M11240.1", {}, source=a, input_identity=identity)
    mb = provenance.build_manifest("analysis", sequence, "M11240.1", {}, source=b, input_identity=identity)
    assert ma["input"]["source"] != mb["input"]["source"]              # recorded verbatim...
    assert ma["input"]["sourceUrl"] == a["sourceUrl"]
    assert mb["input"]["sourceUrl"] == b["sourceUrl"]
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
    sequence = "ACGTACGT"
    acc_src = _validated_accession_source("NC_000013.11", sequence)
    identity = {**public_record_identity("NC_000013.11", trust="validated_fetch"),
                "identitySource": "validated_fetch"}
    acc = provenance.build_manifest("analysis", sequence, "NC_000013.11", {}, source=acc_src,
                                    input_identity=identity)["manifestSha256"]
    assert _seal(_coord_src()) != acc


@pytest.mark.parametrize("source_url", [
    "file:///C:/Users/person/private/sample.fa",
    r"C:\Users\person\private\sample.fa",
    r"\\server\share\sample.fa",
    "https://example.com/path\nprivate",
    "https://user:secret@example.com/record",
    "http://localhost/record",
    "http://127.0.0.1/record",
    "http://10.0.0.8/record",
    "http://192.168.1.8/record",
    "http://[::1]/record",
    "https://www.ncbi.nlm.nih.gov.attacker.example/nuccore/M11240.1",
    "https://www.ncbi.nlm.nih.gov/nuccore/M11240.2",
    "https://www.ncbi.nlm.nih.gov/nuccore/M11240.1?download=1",
    "https://www.ncbi.nlm.nih.gov/nuccore/M11240.1#fragment",
    "https://www.ncbi.nlm.nih.gov/nuccore/%2e%2e/M11240.1",
    "https://www.ncbi.nlm.nih.gov/nuccore/M11240.1%0aextra",
    "https://www.ebi.ac.uk/ena/browser/view/M11240.2",
])
def test_manifest_omits_nonpublic_or_credentialed_source_urls(source_url):
    manifest = provenance.build_manifest(
        "analysis", "ACGT", "M11240.1", {}, source={
            "accession": "M11240.1", "source": "NCBI nuccore", "sourceUrl": source_url,
        })
    assert "sourceUrl" not in manifest["input"]
    assert source_url not in json.dumps(manifest, ensure_ascii=False)
    assert provenance.verify_manifest(manifest)


@pytest.mark.parametrize("source_url", [
    "https://www.ncbi.nlm.nih.gov/nuccore/M11240.1",
    "https://www.ebi.ac.uk/ena/browser/view/M11240.1",
])
def test_manifest_retains_public_ncbi_and_ena_source_urls(source_url):
    ena = "ebi.ac.uk" in source_url
    source = _validated_accession_source("M11240.1", "ACGT", ena=ena)
    identity = {**public_record_identity("M11240.1", trust="validated_fetch"),
                "identitySource": "validated_fetch"}
    manifest = provenance.build_manifest("analysis", "ACGT", "M11240.1", {}, source=source,
                                         input_identity=identity)
    assert manifest["input"]["sourceUrl"] == source_url
    assert provenance.verify_manifest(manifest)


def test_manifest_source_url_is_canonical_accession_bound_and_semantically_verified():
    source_url = "https://www.ncbi.nlm.nih.gov/nuccore/M11240.1"
    identity = {**public_record_identity("M11240.1", trust="validated_fetch"),
                "identitySource": "validated_fetch"}
    manifest = provenance.build_manifest(
        "analysis", "ACGT", "M11240.1", {}, source=_validated_accession_source(),
        input_identity=identity)
    for mutation in (
        "https://www.ncbi.nlm.nih.gov/nuccore/M11240.2",
        "https://www.ebi.ac.uk/ena/browser/view/M11240.1",
        "https://attacker.example/nuccore/M11240.1",
    ):
        changed = copy.deepcopy(manifest)
        changed["input"]["sourceUrl"] = mutation
        assert not provenance.verify_manifest(changed)
        with pytest.raises(ValueError, match="seal verification"):
            provenance.manifest_json_bytes(changed)


def test_resealed_fetch_source_and_public_identity_mismatch_is_rejected():
    identity = {**public_record_identity("M11240.1", trust="validated_fetch"),
                "identitySource": "validated_fetch"}
    manifest = provenance.build_manifest(
        "analysis", "ACGT", "M11240.1", {}, source=_validated_accession_source(),
        input_identity=identity)
    changed = copy.deepcopy(manifest)
    changed["input"]["id"] = "M11241.1"
    changed["manifestSha256"] = hashlib.sha256(provenance.canonical_manifest_bytes(changed)).hexdigest()
    assert not provenance.verify_manifest(changed)


def test_resealed_coordinate_region_and_record_identity_mismatch_is_rejected():
    source = _validated_coord(_coord_src())
    identity = {**public_record_identity("NC_000013.11:1-8", trust="validated_fetch"),
                "identitySource": "validated_fetch"}
    manifest = provenance.build_manifest(
        "analysis", "ACGTACGT", identity["id"], {}, source=source, input_identity=identity,
        input_records=[{**identity, "sequence": "ACGTACGT"}])
    changed = copy.deepcopy(manifest)
    changed["input"]["records"][0]["id"] = "NC_000013.11:2-9"
    changed["manifestSha256"] = hashlib.sha256(provenance.canonical_manifest_bytes(changed)).hexdigest()
    assert not provenance.verify_manifest(changed)


@pytest.mark.parametrize("record_id", [
    "patient_042", "donor-A1", "sample_tuna_birgun", "descriptive specimen",
    " specimen ", "specimen\nprivate", r"C:\private\sample.fa", "/private/sample.fa",
])
def test_manifest_rewrites_nonpublic_record_identity_without_raw_label(record_id):
    manifest = provenance.build_manifest(
        "analysis", "ACGT", record_id, {}, input_records=[{"id": record_id, "sequence": "ACGT"}])
    identity = provenance.public_record_identity(record_id)
    assert manifest["input"]["id"] == identity["id"]
    assert manifest["input"]["idSha256"] == identity["idSha256"]
    assert manifest["input"]["records"][0]["idSha256"] == identity["idSha256"]
    assert record_id not in provenance.manifest_json_bytes(manifest).decode("utf-8")


def test_compact_display_collision_keeps_full_digest_and_order_bindings(monkeypatch):
    from teagle_core import gff3

    monkeypatch.setattr(gff3, "_identity_digest", lambda _value: "deadbeefdead")
    real_sha = gff3._identity_sha256
    monkeypatch.setattr(gff3, "_identity_sha256",
                        lambda value: "deadbeefdead" + real_sha(value)[12:])
    raw = ["patient_042", "donor-A1"]
    top = {**gff3.public_record_identity("multiple_records"), "identitySource": "local_input"}
    manifest = provenance.build_manifest("analysis", "ACGTACGT", "multiple_records", {}, input_identity=top,
                                         input_records=[{"id": raw[0], "sequence": "ACGT"},
                                                        {"id": raw[1], "sequence": "ACGT"}])
    records = manifest["input"]["records"]
    assert [record["id"] for record in records] == ["record_deadbeefdead"] * 2
    assert len({record["idSha256"] for record in records}) == 2
    assert provenance.verify_manifest(manifest)
    changed = copy.deepcopy(manifest)
    changed["input"]["records"].reverse()
    assert not provenance.verify_manifest(changed)


def test_resealed_compact_identity_prefix_mismatch_is_rejected():
    manifest = provenance.build_manifest("analysis", "ACGT", "patient_042", {})
    changed = copy.deepcopy(manifest)
    changed["input"]["id"] = "record_000000000000"
    changed["manifestSha256"] = hashlib.sha256(provenance.canonical_manifest_bytes(changed)).hexdigest()
    assert not provenance.verify_manifest(changed)


@pytest.mark.parametrize("field", ["selected_record", "assay_input"])
def test_resealed_nested_compact_identity_prefix_mismatch_is_rejected(field):
    identity = public_record_identity("patient_042")
    manifest = provenance.build_manifest("primer", "ACGT", "patient_042", {},)
    changed = copy.deepcopy(manifest)
    changed["parameters"][field] = {**identity, "length": 4, "sha256": hashlib.sha256(b"ACGT").hexdigest()}
    changed["parameters"][field]["id"] = "record_000000000000"
    changed["manifestSha256"] = hashlib.sha256(provenance.canonical_manifest_bytes(changed)).hexdigest()
    assert not provenance.verify_manifest(changed)


def test_resealed_pcr_source_role_digest_mismatch_is_rejected():
    identity = public_record_identity("patient_042")
    manifest = provenance.build_manifest("pcr", "ACGT", "patient_042", {})
    changed = copy.deepcopy(manifest)
    changed["parameters"]["amplicon"] = {
        "source": identity["id"] + " (input)", "sourceIdSha256": "0" * 64}
    changed["manifestSha256"] = hashlib.sha256(provenance.canonical_manifest_bytes(changed)).hexdigest()
    assert not provenance.verify_manifest(changed)


def test_manifest_public_record_binding_never_serializes_local_identifier():
    raw = r"C:\Users\person\private\sample.fa"
    manifest = provenance.build_manifest(
        "analysis", "ACGT", raw, {}, input_records=[{"id": raw, "sequence": "ACGT"}])
    expected = provenance.public_record_identity(raw)
    assert manifest["input"]["id"] == expected["id"]
    assert manifest["input"]["idSha256"] == expected["idSha256"]
    assert manifest["input"]["records"][0]["idSha256"] == expected["idSha256"]
    assert raw not in provenance.manifest_json_bytes(manifest).decode("utf-8")
