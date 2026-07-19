from teagle_core import sequtil


def test_parse_fasta_bare_and_records():
    assert sequtil.parse_fasta("acgt\nAC GT")[0] == ("input_sequence", "ACGTACGT")
    recs = sequtil.parse_fasta(">a desc\nACGT\nGGGG\n>b\nTTTT")
    assert [r[0] for r in recs] == ["a desc", "b"]
    assert recs[0][1] == "ACGTGGGG"
    assert recs[1][1] == "TTTT"


def test_parse_fasta_crlf_and_empty():
    assert sequtil.parse_fasta(">x\r\nAC\r\nGT\r\n")[0][1] == "ACGT"
    assert sequtil.parse_fasta("") == []


def test_validate_iupac():
    assert sequtil.validate_iupac("ACGTRYSWKMN") == (True, [])
    ok, bad = sequtil.validate_iupac("ACZG*T")
    assert ok is False
    assert (2, "Z") in bad and (4, "*") in bad


def test_composition():
    c = sequtil.composition("GGCCAATTNN")
    assert c["length"] == 10
    assert c["gc"] == 50.0            # 4 GC of 8 non-N ACGT
    assert c["n"] == 20.0


def test_reverse_complement_involution_and_iupac():
    for s in ["ACGT", "AACCGGTTN", "ATGCRYSWKMBDHVN"]:
        assert sequtil.reverse_complement(sequtil.reverse_complement(s)) == s
    assert sequtil.reverse_complement("AAAAC") == "GTTTT"
    assert sequtil.reverse_complement("N") == "N"


def test_find_orfs_known():
    # ATG + 50 codons (no stop) + TAA on the + strand
    body = ("AAC" * 50)
    seq = "GG" + "ATG" + body + "TAA" + "GG"
    orfs = sequtil.find_orfs(seq, min_aa=40)
    plus = [o for o in orfs if o["strand"] == "+"]
    assert plus, "expected a + strand ORF"
    top = plus[0]
    assert top["start"] == 2 and top["end"] == 2 + 3 + len(body) + 3
    assert top["length_aa"] == 51           # ATG + 50 body codons (initiator M counted)


def test_sha256_stable():
    assert sequtil.sha256("ACGT") == sequtil.sha256("ACGT")
    assert len(sequtil.sha256("ACGT")) == 64


def test_parse_fasta_tolerates_bad_input():
    # malformed request payloads must never crash the parser
    assert sequtil.parse_fasta(None) == []
    assert sequtil.parse_fasta("") == []
    assert sequtil.parse_fasta(12345) == [("input_sequence", "12345")]   # coerced, not crashed
    assert sequtil.parse_fasta("   \n\t ") == []
