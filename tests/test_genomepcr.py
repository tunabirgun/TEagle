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
