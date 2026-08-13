import envcheck
from teagle_core import provenance


def test_parse_requirements():
    reqs = dict(envcheck.parse_requirements())
    assert "primer3-py" in reqs and "PySide6" in reqs         # native app deps (Qt replaced pywebview)
    assert reqs["primer3-py"]                    # pinned version present


def test_every_bundled_direct_dependency_has_a_functional_import_probe():
    requirements = {name for name, _ in envcheck.parse_requirements()}
    assert requirements <= set(envcheck._IMPORT_NAME)
    assert envcheck._IMPORT_NAME["openpyxl"] == "openpyxl"
    assert envcheck._IMPORT_NAME["certifi"] == "certifi"
    assert envcheck._IMPORT_NAME["pypdf"] == "pypdf"


def test_metadata_without_a_working_import_requires_repair(monkeypatch):
    reqs = [("openpyxl", "3.1.5"), ("certifi", "2026.7.22")]
    monkeypatch.setattr(envcheck, "parse_requirements", lambda: reqs)
    monkeypatch.setattr(envcheck, "_installed", lambda name: dict(reqs)[name])
    monkeypatch.setattr(envcheck, "_importable", lambda name: name != "openpyxl")
    monkeypatch.setattr(envcheck, "probe_backends", lambda: {})
    rep = envcheck.check()
    packages = {p["name"]: p for p in rep["packages"]}
    assert packages["openpyxl"]["ok"] is False
    assert packages["certifi"]["ok"] is True
    assert rep["packages_ok"] is False
    assert rep["needs_install"] is True
    assert rep["importFailures"] == ["openpyxl"]


def test_check_shape_and_python_ok():
    rep = envcheck.check()
    for k in ("python", "python_ok", "packages", "packages_ok", "signature",
              "needs_install", "first_run", "backends"):
        assert k in rep
    assert rep["python_ok"] is True
    assert isinstance(rep["packages"], list) and rep["packages"]


def test_signature_stable():
    reqs = envcheck.parse_requirements()
    assert envcheck._signature(reqs) == envcheck._signature(reqs)


def test_manifest_real_versions_and_stable_checksum():
    m1 = provenance.build_manifest("analysis", "ACGTACGT", "seqA", {"x": 1})
    m2 = provenance.build_manifest("analysis", "ACGTACGT", "seqA", {"x": 1})
    assert m1["input"]["sha256"] == m2["input"]["sha256"]     # deterministic input hash
    names = {s["name"].lower() for s in m1["software"]}
    assert any("hmmer" in n for n in names)                   # an analysis run uses HMMER (pyhmmer)...
    assert not any("primer3" in n for n in names)             # ...and must not seal an unused tool
    assert m1["software"][0]["version"]                       # real version string
    assert m1["notRun"]                                       # deferred items disclosed
    assert len(m1["manifestSha256"]) == 64


def test_manifest_source_provenance():
    src = {"accession": "M11240.1", "organism": "Drosophila melanogaster", "taxid": 7227,
           "source": "NCBI nuccore", "endpoint": "https://x", "retrievedUtc": "2026-01-01T00:00:00Z"}
    m = provenance.build_manifest("analysis", "ACGT", "M11240.1", {}, source=src)
    assert m["input"]["accession"] == "M11240.1"
    assert m["input"]["organism"] == "Drosophila melanogaster"


def test_manifest_seal_is_reproducible():
    # identical scientific input must yield an identical seal, regardless of wall-clock time
    a = provenance.build_manifest("primer", "ACGTACGT", "seq1", {"prod_min": 100})
    b = provenance.build_manifest("primer", "ACGTACGT", "seq1", {"prod_min": 100})
    assert a["manifestSha256"] == b["manifestSha256"]
    # a different input must change the seal
    c = provenance.build_manifest("primer", "ACGTACGA", "seq1", {"prod_min": 100})
    assert c["manifestSha256"] != a["manifestSha256"]
