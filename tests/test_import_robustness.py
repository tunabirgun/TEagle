"""Import robustness: malformed input must be rejected clearly, never mis-parsed silently.

The failure that matters here is not a crash — it is a file that parses into something plausible but
wrong, because every downstream number is then computed from the wrong sequence and nothing announces
it. A pasted protein read as nucleotide, digits from a copied alignment absorbed into the sequence, or a
BOM swallowed into the first base all produce a result that looks ordinary and is not.
"""
import gzip
import io
import os
import sys

import pytest

_BE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "app", "backend")
if _BE not in sys.path:
    sys.path.insert(0, _BE)
from teagle_core import sequtil                                  # noqa: E402
import engine                                                    # noqa: E402


# ---------------- record splitting ----------------
def test_bare_sequence_becomes_one_record():
    recs = sequtil.parse_fasta("ACGTACGT")
    assert recs == [("input_sequence", "ACGTACGT")]


def test_multi_record_split_is_exact():
    recs = sequtil.parse_fasta(">a desc\nACGT\nACGT\n>b\nTTTT\n")
    assert [r[0] for r in recs] == ["a desc", "b"]
    assert [r[1] for r in recs] == ["ACGTACGT", "TTTT"]


def test_crlf_and_lone_cr_are_normalised():
    for text in (">a\r\nACGT\r\nACGT\r\n", ">a\rACGT\rACGT\r"):
        assert sequtil.parse_fasta(text) == [("a", "ACGTACGT")]


def test_duplicate_ids_are_kept_as_separate_records():
    """Silently merging two records that share a name would fabricate a chimeric sequence."""
    recs = sequtil.parse_fasta(">dup\nAAAA\n>dup\nCCCC\n")
    assert len(recs) == 2
    assert [r[1] for r in recs] == ["AAAA", "CCCC"]


def test_empty_record_is_preserved_not_dropped():
    recs = sequtil.parse_fasta(">empty\n>next\nACGT\n")
    assert [r[0] for r in recs] == ["empty", "next"]
    assert recs[0][1] == ""


def test_headerless_blank_input_yields_nothing():
    assert sequtil.parse_fasta("") == []
    assert sequtil.parse_fasta("   \n \n") == []
    assert sequtil.parse_fasta(None) == []


def test_internal_whitespace_and_case_are_normalised():
    recs = sequtil.parse_fasta(">a\nac gt\n  AC\tGT \n")
    assert recs[0][1] == "ACGTACGT"


def test_rna_is_normalised_to_dna_but_protein_keeps_selenocysteine():
    assert sequtil.parse_fasta(">a\nACGU\n")[0][1] == "ACGT"
    assert sequtil.parse_protein(">p\nMSTU\n")[0][1] == "MSTU"      # U = selenocysteine, never U->T


# ---------------- characters that must be caught, not absorbed ----------------
def test_digits_from_a_copied_alignment_are_reported_as_invalid():
    """Sequence copied out of an alignment viewer carries column numbers. Absorbing them shifts every
    downstream coordinate silently."""
    ok, bad = sequtil.validate_iupac(sequtil.parse_fasta(">a\n1 ACGT 60\n")[0][1])
    assert not ok
    assert [c for _, c in bad] == ["1", "6", "0"]


def test_invalid_characters_report_their_positions():
    ok, bad = sequtil.validate_iupac("ACGT@ACGT")
    assert not ok and bad and bad[0][0] == 4


def test_iupac_ambiguity_codes_are_accepted():
    ok, bad = sequtil.validate_iupac("ACGTRYSWKMBDHVN")
    assert ok and bad == []


# ---------------- alphabet guards at the engine boundary ----------------
def test_accession_pasted_into_the_sequence_box_is_rejected():
    with pytest.raises(engine.BadRequest) as e:
        engine.run_primers({"sequence": "NM_001301717.2"})
    assert "accession" in str(e.value).lower()


def test_protein_where_nucleotide_expected_is_rejected():
    with pytest.raises(engine.BadRequest):
        engine.run_primers({"sequence": ">p\nMSTAVLENPGLGRKLSDFGQETSYIEDNCNQNGAISLIFSLKEEVGALAKVLRLFEE\n"})


def test_nucleotide_where_protein_expected_is_rejected():
    with pytest.raises(engine.BadRequest) as e:
        engine.run_miniprot({"sequence": ">g\n" + "ACGT" * 50, "protein": ">p\n" + "ACGT" * 30})
    assert "amino acid" in str(e.value).lower() or "nucleotide" in str(e.value).lower()


def test_non_nucleotide_primer_is_rejected_rather_than_scanned_as_empty():
    with pytest.raises(engine.BadRequest):
        engine.run_pcr({"sequence": ">a\n" + "ACGT" * 50, "fwd": "not a primer", "rev": "ACGTACGTACGTACGT"})


def test_empty_sequence_is_a_clear_user_error_not_a_crash():
    with pytest.raises(engine.BadRequest):
        engine.run_primers({"sequence": ""})


# ---------------- gzip upload ----------------
def test_gzipped_fasta_round_trips():
    raw = b">a\nACGTACGTAA\n"
    blob = gzip.compress(raw)
    text = gzip.decompress(blob).decode()
    assert sequtil.parse_fasta(text) == [("a", "ACGTACGTAA")]


def test_utf8_bom_is_not_absorbed_into_the_first_base():
    """A file saved by Excel or Notepad carries a BOM; read naively it becomes a leading invalid char."""
    text = b"\xef\xbb\xbf>a\nACGT\n".decode("utf-8-sig")
    recs = sequtil.parse_fasta(text)
    assert recs == [("a", "ACGT")]
    ok, _ = sequtil.validate_iupac(recs[0][1])
    assert ok


# ---------------- analysis-level reporting ----------------
def test_multi_record_input_is_announced():
    """The wording changed with the record summary table: all records ARE analysed, and the warning now
    says the downstream steps follow the selected one rather than implying the rest were ignored."""
    res = engine.analyze(">a\nACGTACGTAC\n>b\nTTTTTTTTTT\n")
    assert len(res["records"]) == 2
    assert res.get("warning") and "2 records" in res["warning"]


def test_rna_input_is_announced_not_silently_converted():
    res = engine.analyze(">a\n" + "ACGU" * 40)
    notes = " ".join(res["records"][0]["notes"]).lower()
    assert "rna" in notes and "u was read as t" in notes


def test_invalid_record_is_flagged_and_does_not_abort_the_run():
    res = engine.analyze(">bad\nACGT@@@ACGT\n")
    rec = res["records"][0]
    assert rec["valid"] is False and rec["invalid"]


# ---------------- multi-record analysis (roadmap 4.3) ----------------
def test_every_record_is_classified_not_only_the_first():
    """engine.analyze always computed all records; the UI discarded all but the first. The summary table
    exists because that work was already done and thrown away."""
    from teagle_core import examples
    fa = "".join(examples.load(a) for a in ("M11240", "M80343", "X01005"))
    res = engine.analyze(fa)
    assert len(res["records"]) == 3
    classes = [r["classification"]["te_class"] for r in res["records"]]
    assert classes == ["LTR/Copia", "LINE", "DNA/Tc1"], classes
    for r in res["records"]:                       # each record carries its OWN evidence, not the first's
        assert r["classification"]["completeness"] is not None
        assert r["composition"]["length"] > 0


def test_multi_record_warning_states_that_all_were_analysed():
    from teagle_core import examples
    res = engine.analyze("".join(examples.load(a) for a in ("M11240", "X01005")))
    w = res["warning"].lower()
    assert "were analysed" in w or "all" in w
    assert "selected" in w, "the warning must say downstream steps follow the SELECTED record"
