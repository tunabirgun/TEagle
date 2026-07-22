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
    monkeypatch.setattr(fetch, "fetch_fasta", lambda acc, served=None: ">M11240.1\n" + "ACGT" * 20)
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
    # build_gene_model now completes the model: the middle exon (only 2 exon features annotated, but the
    # first CDS segment implies a third) is derived -> 3 exons, and derived_exons flips True.
    assert gm["counts"] == {"exons": 3, "introns": 2, "cds": 2}
    assert gm["derived_introns"] is False and gm["derived_exons"] is True
    assert gm["exons"][0] == {"start": 2185, "end": 2227, "strand": "+", "note": "1"}   # 1-based -> 0-based
    assert any(e.get("derived") and (e["start"], e["end"]) == (2406, 2610) for e in gm["exons"])   # middle exon
    assert (gm["cds"][0]["start"], gm["cds"][0]["end"]) == (2423, 2610)                  # spliced CDS interval


def test_complete_gene_model_fills_cds_implied_exon():
    # J00265-shaped: exon features for exon 1 + the last exon only, a 2-segment CDS, 2 explicit introns.
    # The middle CDS segment is in no exon -> a derived exon must be synthesized at its true (inter-intron) extent.
    gm = {"exons": [{"start": 2185, "end": 2227, "strand": "+"}, {"start": 3396, "end": 3615, "strand": "+"}],
          "introns": [{"start": 2227, "end": 2406}, {"start": 2610, "end": 3396}],
          "cds": [{"start": 2423, "end": 2610, "strand": "+"}, {"start": 3396, "end": 3542, "strand": "+"}],
          "counts": {"exons": 2, "introns": 2, "cds": 2}, "derived_exons": False, "derived_introns": False}
    out = fetch.complete_gene_model(gm)
    assert out["counts"]["exons"] == 3 and out["derived_exons"] is True
    mid = [e for e in out["exons"] if e.get("derived")]
    assert len(mid) == 1 and (mid[0]["start"], mid[0]["end"]) == (2406, 2610)   # spans between the two introns
    # every CDS segment now lies inside an exon
    assert all(any(e["start"] <= c["start"] and c["end"] <= e["end"] for e in out["exons"]) for c in out["cds"])


def test_cross_check_models_matches_within_tolerance():
    ann = [{"start": 100, "end": 200}, {"start": 400, "end": 500}]
    aln = [{"start": 101, "end": 199}, {"start": 402, "end": 500}]     # both within ±2 bp
    cc = fetch.cross_check_models(ann, aln, tol=2)
    assert cc["matched"] == 2 and not cc["annotation_only"] and not cc["aligned_only"]


def test_cross_check_models_flags_discrepancies():
    ann = [{"start": 100, "end": 200}, {"start": 400, "end": 500}]     # annotation has 2 introns
    aln = [{"start": 100, "end": 200}, {"start": 700, "end": 800}]     # alignment: 1 shared, 1 novel
    cc = fetch.cross_check_models(ann, aln, tol=2)
    assert cc["matched"] == 1
    assert len(cc["annotation_only"]) == 1 and len(cc["aligned_only"]) == 1


def test_complete_gene_model_recomputes_derived_introns():
    # partial exons + a CDS-implied middle exon + DERIVED introns (no explicit intron features): inserting the
    # middle exon must NOT leave a stale intron overlapping it — the derived introns are recomputed.
    gm = {"exons": [{"start": 0, "end": 100}, {"start": 500, "end": 600}],
          "introns": [{"start": 100, "end": 500}],                # derived from the 2-exon gap
          "cds": [{"start": 200, "end": 300}],
          "counts": {"exons": 2, "introns": 1, "cds": 1}, "derived_exons": False, "derived_introns": True}
    out = fetch.complete_gene_model(gm)
    assert out["counts"]["exons"] == 3 and out["counts"]["introns"] == 2
    for i in out["introns"]:                                       # no intron may overlap any exon
        for e in out["exons"]:
            assert not (i["start"] < e["end"] and e["start"] < i["end"]), ("overlap", i, e)


def test_explicit_exons_derived_introns_flag_no_derived_exon():
    # a normal multi-exon record (explicit exons, introns derived from the gaps) must NOT flag any exon as
    # derived — otherwise the gene-model title/legend would claim CDS-inferred exons the figure never draws
    # (Loop-6 MED: gate the exon* legend on the per-exon flag, not the model-level derived_introns).
    ft = ">Feature gb|X.1|Y\n100\t200\texon\n300\t400\texon\n"
    gm = fetch.build_gene_model(fetch.parse_feature_table(ft))
    assert gm["derived_introns"] is True and gm["derived_exons"] is False
    assert not any(e.get("derived") for e in gm["exons"])


def test_complete_gene_model_idempotent():
    gm = {"exons": [{"start": 0, "end": 50}], "introns": [], "cds": [{"start": 100, "end": 200}],
          "counts": {"exons": 1, "introns": 0, "cds": 1}, "derived_exons": False, "derived_introns": False}
    once = fetch.complete_gene_model(gm)
    twice = fetch.complete_gene_model(once)
    assert once["counts"] == twice["counts"] and len(once["exons"]) == len(twice["exons"])


def test_gene_model_derives_introns_from_cds_join():
    # only a spliced CDS (no explicit exon/intron) -> exons + introns derived from the join
    ft = ">Feature gb|X.1|Y\n100\t200\tCDS\n300\t400\n\t\t\tproduct\tp\n"
    gm = fetch.build_gene_model(fetch.parse_feature_table(ft))
    assert gm["derived_exons"] and gm["derived_introns"]
    assert gm["counts"]["exons"] == 2 and gm["counts"]["introns"] == 1
    assert (gm["introns"][0]["start"], gm["introns"][0]["end"]) == (200, 299)            # gap between exon ends
    # honesty invariant: model-level derived_exons implies EVERY wholesale-derived exon is flagged, so the
    # viewer marks them exon* and the legend's derived/annotated distinction matches the figure (Loop-5 MED)
    assert all(e.get("derived") for e in gm["exons"])


@pytest.mark.network
def test_retrieve_real_accession():
    meta = fetch.retrieve("M11240")
    assert "Drosophila" in meta["organism"]
    assert meta["seq_length"] == 5146
    assert meta["fasta"].startswith(">")
    assert meta["endpoint"].startswith("https://")


# ---------- coordinate fetch (UCSC-style) ----------
def test_parse_regions_single():
    r = fetch.parse_regions("chr13:33,016,423-33,066,143")
    assert len(r) == 1
    assert r[0] == {"chromKey": "13", "chromLabel": "chr13", "start": 33016423, "end": 33066143}


def test_parse_regions_comma_is_thousands_not_delimiter():
    # a comma inside the coordinate is a thousands separator, never a region delimiter
    r = fetch.parse_regions("chr1:1,000-2,000")
    assert len(r) == 1 and r[0]["start"] == 1000 and r[0]["end"] == 2000


def test_parse_regions_multi_newline_and_semicolon():
    r = fetch.parse_regions("chr1:1-10\nchr2:5-9 ; chrX:100-200")
    assert [x["chromLabel"] for x in r] == ["chr1", "chr2", "chrX"]
    assert r[1] == {"chromKey": "2", "chromLabel": "chr2", "start": 5, "end": 9}


def test_parse_regions_organism_specific_names_kept_verbatim():
    # roman numerals (yeast), arm names (fly), no 'chr' prefix — passed through for the resolver
    assert fetch.parse_regions("II:1-100")[0]["chromKey"] == "II"
    assert fetch.parse_regions("2L:1-100")[0]["chromKey"] == "2L"
    assert fetch.parse_regions("chrMT:1-5")[0]["chromKey"] == "MT"


@pytest.mark.parametrize("bad", ["", "   ", "chr1", "chr1:100", "chr1:abc-def", "chr1:200-100",
                                 "chr1:0-100", "chr1:,-100", "chr1:100-,", "chr1: , - , "])
def test_parse_regions_rejects_malformed(bad):
    # a comma/space-only coordinate group must raise CoordError, never a bare int('') ValueError
    with pytest.raises(fetch.CoordError):
        fetch.parse_regions(bad)


def test_fetch_fasta_records_ena_fallback_source(monkeypatch):
    # efetch returns a non-FASTA body, ENA has the sequence: the served-by must be reported as ENA, not NCBI
    def fake_get(url, *a, **k):
        if "efetch" in url:
            return "<html>error</html>"                       # NCBI can't serve it as FASTA
        return ">X\n" + "ACGT" * 10                           # ENA fallback serves it
    monkeypatch.setattr(fetch, "_get", fake_get)
    served = []
    fetch.fetch_fasta("M11240", served)
    assert served and served[0][0].startswith("ENA")


def test_fetch_fasta_ena_fallback_on_ncbi_request_error(monkeypatch):
    # an NCBI HTTP/URL error (not just a non-FASTA 200) must still fall through to ENA
    def fake_get(url, *a, **k):
        if "efetch" in url:
            raise fetch.FetchError("HTTP 503 from source")
        return ">X\n" + "ACGT" * 10
    monkeypatch.setattr(fetch, "_get", fake_get)
    served = []
    fetch.fetch_fasta("M11240", served)
    assert served and served[0][0].startswith("ENA")


def test_fetch_fasta_records_ncbi_when_efetch_serves(monkeypatch):
    monkeypatch.setattr(fetch, "_get", lambda url, *a, **k: ">X\n" + "ACGT" * 10)
    served = []
    fetch.fetch_fasta("M11240", served)
    assert served and served[0][0] == "NCBI nuccore"


def test_resolve_non_json_body_is_fetcherror(monkeypatch):
    # accession metadata path: a 200 with a non-JSON body must be a clean FetchError, not JSONDecodeError/500
    monkeypatch.setattr(fetch, "_get", lambda *a, **k: "<html>maintenance</html>")
    with pytest.raises(fetch.FetchError):
        fetch.resolve("M11240")


def test_datasets_json_non_json_body_is_coorderror(monkeypatch):
    # NCBI Datasets 200 with an HTML/maintenance body must surface as a clean CoordError, not JSONDecodeError
    monkeypatch.setattr(fetch, "_get", lambda *a, **k: "<html>maintenance</html>")
    with pytest.raises(fetch.CoordError):
        fetch._datasets_json("genome/taxon/foo/dataset_report")
    with pytest.raises(fetch.FetchError):                     # and it propagates as FetchError through resolve_assembly
        fetch.resolve_assembly("Homo sapiens")


_FAKE_MAP = {"assemblyAccession": "GCF_TEST.1", "molecules": [
    {"chrName": "13", "ucscStyleName": "chr13", "refseqAccession": "NC_000013.11", "length": 114364328},
    {"chrName": "I", "ucscStyleName": "chrI", "refseqAccession": "NC_001133.9", "length": 230218},
    {"chrName": "MT", "ucscStyleName": "chrM", "refseqAccession": "NC_012920.1", "length": 16569},
]}


def test_resolve_chrom_matches_and_mt_alias(monkeypatch):
    monkeypatch.setattr(fetch, "_load_assembly_map", lambda acc, refresh=False: _FAKE_MAP)
    assert fetch.resolve_chrom("GCF_TEST.1", "13")["refseqAccession"] == "NC_000013.11"
    assert fetch.resolve_chrom("GCF_TEST.1", "I")["refseqAccession"] == "NC_001133.9"    # roman numeral
    assert fetch.resolve_chrom("GCF_TEST.1", "M")["refseqAccession"] == "NC_012920.1"    # M -> MT alias
    assert fetch.resolve_chrom("GCF_TEST.1", "chrMT")["refseqAccession"] == "NC_012920.1"
    with pytest.raises(fetch.CoordError):
        fetch.resolve_chrom("GCF_TEST.1", "99")


def test_fetch_fasta_range_length_guard(monkeypatch):
    # efetch silently clamps a stop past the sequence end — the length assert must catch it
    monkeypatch.setattr(fetch, "_get", lambda *a, **k: ">x:1-100\n" + "A" * 50 + "\n")
    with pytest.raises(fetch.FetchError):
        fetch.fetch_fasta_range("NC_000013.11", 1, 100, "+")


def _stub_coord_backend(monkeypatch, tmp_path):
    monkeypatch.setattr(fetch, "CACHE_DIR", str(tmp_path))
    monkeypatch.setattr(fetch, "resolve_chrom",
                        lambda acc, key, refresh=False: {"refseqAccession": "NC_000013.11",
                        "length": 114364328, "ucscStyleName": "chr13", "chrName": "13"})
    calls = {"n": 0}
    def fake_range(acc, start, end, strand="+"):
        calls["n"] += 1
        return f">{acc}:{start}-{end}\n" + "ACGT" * ((end - start + 1 + 3) // 4)
    monkeypatch.setattr(fetch, "fetch_fasta_range", fake_range)
    return calls


def test_retrieve_coords_caches_and_strand_keys_differ(tmp_path, monkeypatch):
    calls = _stub_coord_backend(monkeypatch, tmp_path)
    a = fetch.retrieve_coords("chr13:1-100", "GCF_TEST.1", "GRCh38", "Homo sapiens")
    assert a["fromCache"] is False and calls["n"] == 1
    assert a["runType"] == "coordinate" and a["source"]["retrievalType"] == "coordinate"
    b = fetch.retrieve_coords("chr13:1-100", "GCF_TEST.1", "GRCh38", "Homo sapiens")
    assert b["fromCache"] is True and calls["n"] == 1                    # served from disk
    c = fetch.retrieve_coords("chr13:1-100", "GCF_TEST.1", "GRCh38", "Homo sapiens", strand="-")
    assert c["fromCache"] is False and calls["n"] == 2                   # strand is part of the cache identity


def test_retrieve_coords_assembly_and_taxid_in_cache_key(tmp_path, monkeypatch):
    # the stub resolves every assembly to the SAME chr accession/coords, so only the assembly+taxid in the
    # cache key can keep two genuinely different runs from colliding (patch assemblies share primary-chr RefSeq IDs)
    calls = _stub_coord_backend(monkeypatch, tmp_path)
    a = fetch.retrieve_coords("chr13:1-100", "GCF_000001405.40", "GRCh38.p14", "Homo sapiens")
    assert a["fromCache"] is False and calls["n"] == 1
    b = fetch.retrieve_coords("chr13:1-100", "GCF_000001405.39", "GRCh38.p13", "Homo sapiens")   # different assembly
    assert b["fromCache"] is False and calls["n"] == 2                                            # must NOT collide
    c = fetch.retrieve_coords("chr13:1-100", "GCF_000001405.40", "GRCh38.p14", "Homo sapiens", taxid="9606")
    assert c["fromCache"] is False and calls["n"] == 3                                            # different taxid path
    d = fetch.retrieve_coords("chr13:1-100", "GCF_000001405.40", "GRCh38.p14", "Homo sapiens")   # identical identity
    assert d["fromCache"] is True and calls["n"] == 3                                             # DOES hit cache


def test_retrieve_coords_multi_region_seal_fields(tmp_path, monkeypatch):
    _stub_coord_backend(monkeypatch, tmp_path)
    m = fetch.retrieve_coords("chr13:1-100\nchr13:200-300", "GCF_TEST.1", "GRCh38", "Homo sapiens")
    assert len(m["source"]["regions"]) == 2
    assert m["source"]["regions"][0] == {"chrAccession": "NC_000013.11", "start": 1, "stop": 100, "strand": 1}
    assert m["source"]["coordSystem"] == "1-based-inclusive"
    assert "+1 more" in m["displayLocus"]


def test_retrieve_coords_rejects_out_of_bounds(tmp_path, monkeypatch):
    monkeypatch.setattr(fetch, "CACHE_DIR", str(tmp_path))
    monkeypatch.setattr(fetch, "resolve_chrom",
                        lambda acc, key, refresh=False: {"refseqAccession": "NC_x.1", "length": 500,
                        "ucscStyleName": "chr13", "chrName": "13"})
    with pytest.raises(fetch.CoordError):
        fetch.retrieve_coords("chr13:1-999", "GCF_TEST.1", "GRCh38", "Homo sapiens")


def test_coord_assemblies_are_pinned_accessions():
    # every curated organism pins a versioned GCF/GCA assembly accession + numeric taxid (complete organism pin)
    for org, a in fetch.COORD_ASSEMBLIES.items():
        assert fetch._ASM_ACC_RE.match(a["assemblyAccession"]), org
        assert a["assemblyName"]
        assert a["taxid"].isdigit(), org                      # taxid recorded so curated runs seal the same identity


def test_retrieve_coords_records_taxid(tmp_path, monkeypatch):
    _stub_coord_backend(monkeypatch, tmp_path)
    m = fetch.retrieve_coords("chr13:1-100", "GCF_000001405.40", "GRCh38.p14", "Homo sapiens", taxid="9606")
    assert m["source"]["taxid"] == "9606"                     # sealed organism identity, not left blank


@pytest.mark.network
def test_retrieve_coords_real_no_off_by_one(tmp_path, monkeypatch):
    monkeypatch.setattr(fetch, "CACHE_DIR", str(tmp_path))
    a = fetch.COORD_ASSEMBLIES["Homo sapiens"]
    m = fetch.retrieve_coords("chr13:33,016,423-33,066,143", a["assemblyAccession"], a["assemblyName"], "Homo sapiens")
    assert m["seq_length"] == 33066143 - 33016423 + 1                    # UCSC display == efetch, no conversion
