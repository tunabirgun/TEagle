"""Local whole-genome in-silico PCR — offline tests of the pure isPcr parser + query builder.
The isPcr run itself (WSL) is exercised by a @wsl test; the FASTA parsing is tested here."""
from teagle_core import genomepcr


def test_query_rows_emits_pair_and_single_primer_rows():
    q = genomepcr.query_rows("acgtACGT", "ttttGGGG")
    rows = [r for r in q.splitlines() if r.strip()]
    assert rows == ["pair\tACGTACGT\tTTTTGGGG",
                    "fwdonly\tACGTACGT\tACGTACGT",
                    "revonly\tTTTTGGGG\tTTTTGGGG"]


def test_parse_ispcr_plus_strand_amplicon():
    out = ">NC_001133.9:1001+1240 pair 240bp TACAATTATATCTTATTTCC ATGAAGTGAGACAATATCGT\nACGT\n"
    amps = genomepcr.parse_ispcr(out)
    assert len(amps) == 1
    a = amps[0]
    assert a["source"] == "NC_001133.9" and (a["start"], a["end"], a["length"]) == (1001, 1240, 240)
    assert a["strand"] == "+" and a["on_target"] is False and a["single_primer"] is False


def test_parse_ispcr_minus_strand_orders_coords_low_high():
    out = ">NC_1:5000-5200 pair 201bp AAA CCC\nACGT\n"
    a = genomepcr.parse_ispcr(out)[0]
    assert a["strand"] == "-" and (a["start"], a["end"]) == (5000, 5200)   # lo..hi regardless of report order


def test_parse_ispcr_flags_single_primer_products():
    out = "\n".join([
        ">NC_1:100+300 pair 201bp AAA CCC",
        ">NC_1:900+1000 fwdonly 101bp AAA AAA",
        ">NC_2:10+80 revonly 71bp CCC CCC",
    ]) + "\n"
    amps = genomepcr.parse_ispcr(out)
    by = {a["pair"]: a for a in amps}
    assert by["pair"]["single_primer"] is False
    assert by["fwdonly"]["single_primer"] is True and by["revonly"]["single_primer"] is True


def test_parse_ispcr_ignores_sequence_and_malformed_lines():
    out = "not a header\n>bad header without coords\nACGTACGT\n>NC_1:1+50 pair 50bp AAA CCC\nACGT\n"
    amps = genomepcr.parse_ispcr(out)
    assert len(amps) == 1 and amps[0]["length"] == 50


# --- summarize(): off-target interpretation (pair-vs-single split, per-chrom spread, verdict, size cluster) ---
def _amps(*headers):
    return genomepcr.parse_ispcr("".join(h + "\n" for h in headers))


def test_summarize_splits_pair_and_single_primer_counts():
    s = genomepcr.summarize(_amps(
        ">NC_1:100+300 pair 201bp AAA CCC",
        ">NC_2:10+210 pair 201bp AAA CCC",
        ">NC_1:5+55 fwdonly 51bp AAA AAA",
        ">NC_2:5+55 revonly 51bp CCC CCC"))
    assert s["n_total"] == 4 and s["n_pair"] == 2 and s["n_single"] == 2


def test_summarize_groups_per_source_pair_only_busiest_first():
    # single-primer products must NOT inflate the per-chromosome pair-product counts
    s = genomepcr.summarize(_amps(
        ">NC_1:100+300 pair 201bp AAA CCC",
        ">NC_1:900+1100 pair 201bp AAA CCC",
        ">NC_2:10+210 pair 201bp AAA CCC",
        ">NC_1:5+55 fwdonly 51bp AAA AAA"))
    assert s["per_source"] == [("NC_1", 2), ("NC_2", 1)] and s["n_sources"] == 2


def test_summarize_verdict_locus_specific_for_single_pair_product():
    s = genomepcr.summarize(_amps(">NC_1:100+300 pair 201bp AAA CCC"))
    assert s["tier"] == "locus-specific" and s["n_pair"] == 1


def test_summarize_verdict_family_generic_for_many_copies():
    s = genomepcr.summarize(_amps(*[f">NC_{i}:100+300 pair 201bp AAA CCC" for i in range(1, 9)]))
    assert s["tier"] == "family-generic" and s["n_pair"] == 8 and s["n_sources"] == 8


def test_summarize_size_cluster_over_pair_products():
    s = genomepcr.summarize(_amps(
        ">NC_1:1+200 pair 200bp AAA CCC",
        ">NC_1:9+208 pair 200bp AAA CCC",
        ">NC_2:1+250 pair 250bp AAA CCC"))
    assert s["size_mode"] == 200 and s["size_mode_n"] == 2
    assert s["size_min"] == 200 and s["size_max"] == 250


def test_summarize_empty_scan_is_clean_not_error():
    # an empty isPcr result is the legitimate 'specific primer' case, not a crash
    s = genomepcr.summarize([])
    assert s["n_total"] == 0 and s["n_pair"] == 0 and s["per_source"] == []
    assert s["tier"] == "none" and s["size_mode"] is None


def test_summarize_has_locus_splits_on_and_off_target():
    # with a design locus, the product marked on_target is the on-target; the rest are off-target paralogs
    amps = _amps(">NC_1:100+300 pair 201bp AAA CCC", ">NC_2:10+210 pair 201bp AAA CCC")
    amps[0]["on_target"] = True
    s = genomepcr.summarize(amps, has_locus=True)
    assert s["has_locus"] and s["n_on"] == 1 and s["n_off"] == 1
    assert "on-target" in s["verdict"] and "off-target" in s["verdict"]


def test_summarize_no_locus_is_neutral_priming_sites():
    # no design locus -> no single intended target, so products are neutral genomic priming sites (no on/off)
    s = genomepcr.summarize(_amps(">NC_1:100+300 pair 201bp AAA CCC"), has_locus=False)
    assert s["has_locus"] is False and s["n_on"] == 0 and "priming site" in s["verdict"]
