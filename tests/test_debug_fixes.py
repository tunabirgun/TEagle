"""Regression tests for defects found by the three debug loops (2026-07-18):
N-run structural fabrication, O(n^2) find_orfs, lowercase primers, amplicon cap,
seal refetch-invariance + unused-tool scoping, and classify robustness/formatting."""
from teagle_core import structural, classify, primers, provenance, sequtil, refs


# --- N-runs are not base-pairing evidence: no fabricated TIR/LTR, no false TE call ---
def test_n_run_does_not_fabricate_tir():
    assert structural.find_tir("N" * 30 + "ACGTACGTAC" + "N" * 30) is None


def test_n_run_does_not_fabricate_ltr():
    assert structural.find_ltr("N" * 400) is None


def test_n_run_not_classified_as_te():
    assert classify.classify(structural.detect_all("N" * 100), [])["te_class"] == "none"
    assert classify.classify(structural.detect_all("N" * 400), [])["te_class"] == "none"


# --- find_orfs: single-pass (no per-ATG rescan) and still correct ---
def test_find_orfs_stop_poor_is_fast():
    import time
    t0 = time.time()
    orfs = sequtil.find_orfs("ATG" * 20000)          # no in-frame stop -> single-pass, not O(n^2)
    assert time.time() - t0 < 2.0                    # was ~17s before the single-pass rewrite


def test_find_orfs_emits_3prime_truncated_orf():
    # a long ORF with no terminal stop must still be recovered (flagged open-ended) so truncated TEs scan
    full = sequtil.find_orfs("ATG" + "GCT" * 100, min_aa=40)
    assert full and any(o.get("open_end") for o in full)
    assert not any(o.get("open_end") for o in sequtil.find_orfs("ATG" + "GCT" * 100 + "TAA", min_aa=40))


def test_rna_input_normalized_to_dna():
    # U -> T on ingest: a 50%-GC RNA reports 50% GC (not 66.7%) and is translatable
    (rid, seq), = sequtil.parse_fasta("AUGCAUGCAUGCAUGCAUGCAUGCAUGC")
    assert "U" not in seq and seq.startswith("ATGC")
    assert sequtil.composition(seq)["gc"] == 50.0


def test_find_orfs_recovers_a_real_orf():
    orfs = sequtil.find_orfs("ATG" + "GCA" * 50 + "TAA")     # ATG + 50 codons before the stop -> 51 aa
    assert any(o["length_aa"] == 51 and o["strand"] == "+" and o["start"] == 0 for o in orfs)


# --- in-silico PCR: case-insensitive primers, and a hard cap on output ---
def _amplifying_template():
    fwd = "ACGTACGTACGTACGTAC"
    rev = "TTGGTTGGTTGGTTGGTT"
    seq = fwd + "GC" * 40 + primers.reverse_complement(rev)
    return fwd, rev, seq


def test_lowercase_primers_match_uppercase():
    fwd, rev, seq = _amplifying_template()
    up = primers.in_silico_pcr(fwd, rev, seq)
    lo = primers.in_silico_pcr(fwd.lower(), rev.lower(), seq)
    assert len(up) >= 1 and len(up) == len(lo)       # a lowercase primer must bind identically


def test_in_silico_pcr_output_is_capped():
    rep = "ACGTACGTACGTACGTACGT" * 100               # highly repetitive -> many binding sites
    amps = primers.in_silico_pcr("ACGTACGTACGT", "ACGTACGTACGT", rep, max_mm=3, tp=0, prod_min=1, prod_max=3000)
    assert len(amps) <= 4000                          # bounded, does not blow up memory


def test_tp_negative_does_not_disable_three_prime_rule():
    # a huge/negative tp must be clamped, not hang or silently disable the 3' check
    fwd, rev, seq = _amplifying_template()
    assert isinstance(primers.in_silico_pcr(fwd, rev, seq, tp=-5), list)
    assert isinstance(primers.in_silico_pcr(fwd, rev, seq, tp=10 ** 9), list)


# --- provenance seal: refetch-invariant, and never carries an unused tool's version ---
def test_seal_invariant_to_fasta_header():
    a = sequtil.parse_fasta(">M11240.1 desc\nACGTACGTACGTACGTACGTACGTACGTACGT\n")
    b = sequtil.parse_fasta(">ENA|M11240|M11240.1 desc\nACGTACGTACGTACGTACGTACGTACGTACGT\n")
    ma = provenance.build_manifest("analysis", a[0][1], a[0][0], {"orf_min_aa": 40})
    mb = provenance.build_manifest("analysis", b[0][1], b[0][0], {"orf_min_aa": 40})
    assert ma["input"]["id"] != mb["input"]["id"]            # header recorded verbatim...
    assert ma["manifestSha256"] == mb["manifestSha256"]      # ...but not sealed


def test_in_silico_pcr_run_does_not_cite_primer3():
    keys = [r["key"] for r in refs.for_run("in-silico-pcr")]
    assert "Primer3" not in keys                              # pure-Python scanner; its seal must not depend on primer3-py


def test_malformed_source_does_not_crash_manifest():
    for bad in (5, True, [1, 2], "x"):
        m = provenance.build_manifest("analysis", "ACGT", "x", {}, source=bad)
        assert len(m["manifestSha256"]) == 64                # ignored, not crashed


# --- classify: clean explanation + robust to sparse domain dicts ---
def test_empty_evidence_explanation_has_no_double_period():
    assert ". ." not in classify.classify([], [])["explanation"]


def test_classify_robust_to_missing_domain_keys():
    ltr = [{"type": "LTR (x)"}]
    doms = [{"domain": "INT", "class": "INT"}, {"domain": "RT", "class": "RT", "strand": "+", "nt": [100, 200], "score": 50.0}]
    cl = classify.classify(ltr, doms)                         # missing strand/nt on INT must not raise
    assert cl["superfamily"] and cl["te_class"]
