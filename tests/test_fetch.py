import pytest
from teagle_core import fetch


def test_validate_accession_valid():
    for acc in ["M11240", "NC_003075.7", "X05424", "AB087615"]:
        assert fetch.validate_accession(acc) == acc


def test_validate_accession_invalid():
    for bad in ["", "   ", "not!an!acc", "12345", "!!!", "M112 40"]:
        with pytest.raises(fetch.FetchError):
            fetch.validate_accession(bad)


def test_retrieve_caches_and_reuses(tmp_path, monkeypatch):
    # the same accession must be served from disk on the second call, not re-downloaded
    monkeypatch.setattr(fetch, "CACHE_DIR", str(tmp_path))
    calls = {"n": 0}
    def fake_resolve(acc):
        calls["n"] += 1
        return {"accession": "M11240.1", "organism": "Drosophila melanogaster", "taxid": 7227,
                "title": "copia", "length": 5146, "moltype": "DNA", "source": "NCBI nuccore"}
    monkeypatch.setattr(fetch, "resolve", fake_resolve)
    monkeypatch.setattr(fetch, "fetch_fasta", lambda acc: ">M11240.1\n" + "ACGT" * 20)
    first = fetch.retrieve("M11240")
    assert first["fromCache"] is False and calls["n"] == 1
    second = fetch.retrieve("M11240")
    assert second["fromCache"] is True and calls["n"] == 1          # no second network hit
    assert second["fasta"] == first["fasta"]
    third = fetch.retrieve("M11240", refresh=True)                  # force re-download
    assert third["fromCache"] is False and calls["n"] == 2


def test_feature_table_parse_and_gene_model():
    # NCBI feature-table (rettype=ft) with explicit exon/intron + a spliced CDS join
    ft = ("\n".join([
        ">Feature gb|X.1|Y",
        "2186\t2227\texon", "\t\t\tnumber\t1",
        "2228\t2406\tintron", "\t\t\tnumber\t1",
        "2424\t2610\tCDS", "3397\t3542", "\t\t\tproduct\tinsulin",
        "2611\t3396\tintron", "\t\t\tnumber\t2",
        "3397\t>3615\texon", "\t\t\tnumber\t2",
    ]) + "\n")
    gm = fetch.build_gene_model(fetch.parse_feature_table(ft))
    assert gm["counts"] == {"exons": 2, "introns": 2, "cds": 2}
    assert gm["derived_introns"] is False and gm["derived_exons"] is False
    assert gm["exons"][0] == {"start": 2185, "end": 2227, "strand": "+", "note": "1"}   # 1-based -> 0-based
    assert (gm["cds"][0]["start"], gm["cds"][0]["end"]) == (2423, 2610)                  # spliced CDS interval


def test_gene_model_derives_introns_from_cds_join():
    # only a spliced CDS (no explicit exon/intron) -> exons + introns derived from the join
    ft = ">Feature gb|X.1|Y\n100\t200\tCDS\n300\t400\n\t\t\tproduct\tp\n"
    gm = fetch.build_gene_model(fetch.parse_feature_table(ft))
    assert gm["derived_exons"] and gm["derived_introns"]
    assert gm["counts"]["exons"] == 2 and gm["counts"]["introns"] == 1
    assert (gm["introns"][0]["start"], gm["introns"][0]["end"]) == (200, 299)            # gap between exon ends


@pytest.mark.network
def test_retrieve_real_accession():
    meta = fetch.retrieve("M11240")
    assert "Drosophila" in meta["organism"]
    assert meta["seq_length"] == 5146
    assert meta["fasta"].startswith(">")
    assert meta["endpoint"].startswith("https://")
