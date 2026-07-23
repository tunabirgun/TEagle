"""WSL backend: hermetic .out-parser tests + a WSL-gated live annotation test."""
import pytest
from teagle_core import wsl

SAMPLE_OUT = """   SW   perc perc perc  query     position in query    matching  repeat       position in repeat
score   div. del. ins.  sequence  begin end   (left)   repeat    class/family begin  end    (left)  ID

 2405    0.0  0.0  0.0  M11240.1      1   276 (4870) + Copia_LTR LTR/Copia       1    276    (0)   1
40720    0.5  0.0  0.0  M11240.1    277  4870  (276) + Copia_I   LTR/Copia       1   4593    (0)   1
  283   14.9  0.0  0.0  Xyz             5   120 (900)  C Gypsy-7  LTR/Gypsy    (0)  200    50   2
"""


def test_parse_out_basic():
    hits = wsl.parse_out(SAMPLE_OUT)
    assert len(hits) == 3
    h = hits[0]
    assert h["family"] == "Copia_LTR" and h["class_family"] == "LTR/Copia"
    assert h["q_start"] == 0 and h["q_end"] == 276     # 1-based .out -> 0-based half-open
    assert h["strand"] == "+" and h["score"] == 2405


def test_parse_out_complement_strand():
    hits = wsl.parse_out(SAMPLE_OUT)
    assert hits[2]["strand"] == "-"                     # 'C' in the .out maps to minus strand
    assert hits[2]["class_family"] == "LTR/Gypsy"


def test_parse_out_skips_headers_and_blanks():
    assert wsl.parse_out("garbage\n\n   header line\n") == []


def test_prep_script_extracts_with_python3_zipfile():
    # a fresh minimal Ubuntu WSL ships python3 (zipfile built into CPython) but NOT unzip; extraction must
    # use python3 stdlib first so a first-time genome download does not fail with 'genome preparation failed: unzip'
    script = wsl._PREP_SCRIPT
    assert "python3 -m zipfile -e" in script            # primary extractor: guaranteed-present stdlib
    assert "FAIL unzip" not in script                   # the old bare-unzip failure path is gone
    if "unzip -oq" in script:                           # a fallback unzip is allowed ONLY if guarded on the system PATH
        assert "command -v unzip" in script


def test_prep_script_keeps_zip_free_and_fna_discovery():
    # the extraction rewrite must not lose the peak-disk (free the zip) or FNA-discovery invariants
    script = wsl._PREP_SCRIPT
    assert "rm -f dl.zip" in script
    assert "_genomic.fna" in script


def test_prep_script_guards_empty_contig_count():
    # a broken/empty FASTA must never be sealed as a valid 0-contig genome
    assert "no sequences" in wsl._PREP_SCRIPT


def test_integrity_probe_verifies_genome_scan_tools():
    # integrity_check must verify isPcr + NCBI datasets, not certify the backend healthy while the scan is broken
    # (same dependency-binding-truthfulness class as the fixed unzip gap)
    p = wsl._INTEGRITY_PROBE
    assert "=SCAN=" in p and "isPcr" in p and "datasets" in p


def test_genome_list_delivered_via_stdin_not_inline():
    # inline `bash -lc` mangles the loop variable $d (wsl.exe command-line rebuild), so an inline for-loop
    # reports EVERY cached genome as missing; genome_list MUST deliver its loop via STDIN (_wsl_script)
    import inspect
    src = inspect.getsource(wsl.genome_list)
    assert "_wsl_script(" in src and "for d in" in src and "_wsl(" not in src.replace("_wsl_script(", "")


def test_genome_scan_distinguishes_missing_ispcr_from_missing_genome():
    # a broken isPcr binary must not be misattributed to a missing genome (which would send the user to
    # re-download a multi-GB assembly) — the run script guards the binary with a distinct exit code
    import inspect
    src = inspect.getsource(wsl.genome_scan)
    assert "exit 8" in src and "isPcr tool is missing" in src


def test_species_validation_rejects_injection():
    # invalid species tokens must be rejected before reaching the shell
    r = wsl.annotate(">x\nACGT", species="dros; rm -rf /")
    assert r["ok"] is False and "invalid species" in r["error"]


def test_start_install_executes_script_not_reaped_nohup(monkeypatch):
    # Regression: WSL2 reaps a detached (nohup) orphan when the launching session exits,
    # so the installer must EXECUTE the script inside a session-holding call, never nohup+return.
    import time
    calls = []
    def fake_wsl(script, stdin=None, timeout=600):
        calls.append(script)
        return 0, "FREE", ""
    def fake_wsl_script(script, timeout=90):        # lock-liveness pre-check goes through _wsl_script
        calls.append(script)
        return 0, "FREE", ""                        # no live install -> FREE, so start proceeds
    monkeypatch.setattr(wsl, "_wsl", fake_wsl)
    monkeypatch.setattr(wsl, "_wsl_script", fake_wsl_script)
    monkeypatch.setattr(wsl, "available", lambda: {"wsl2": True})
    wsl._install_thread = None
    r = wsl.start_install()
    assert r["started"] is True
    for _ in range(100):                            # wait for the session-holding thread to run the script
        if any('bash "$HOME/teagle_wsl_install.sh"' in c for c in calls):
            break
        time.sleep(0.02)
    if wsl._install_thread:
        wsl._install_thread.join(timeout=2)
    assert any('bash "$HOME/teagle_wsl_install.sh"' in c for c in calls), calls
    assert not any("nohup" in c for c in calls), "installer must not use a reaped nohup launch"


def test_parse_sam_splice_cigar():
    # a spliced alignment: 50M 200N 50M at 1-based POS 101 -> two exons, one intron (0-based)
    sam = "@HD\tVN:1.6\nq1\t0\tref\t101\t60\t50M200N50M\t*\t0\t0\t*\t*\n"
    r = wsl._parse_sam_splice(sam)
    assert r["strand"] == "+"
    assert r["exons"] == [{"start": 100, "end": 150}, {"start": 350, "end": 400}]
    assert r["introns"] == [{"start": 150, "end": 350}]


def test_parse_sam_splice_skips_unmapped():
    sam = "q1\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n"    # FLAG 4 = unmapped
    assert wsl._parse_sam_splice(sam) is None


def test_start_install_blocks_concurrent(monkeypatch):
    # a second install must not start while one is running
    monkeypatch.setattr(wsl, "available", lambda: {"wsl2": True})
    class _Alive:
        def is_alive(self): return True
    monkeypatch.setattr(wsl, "_install_thread", _Alive())
    r = wsl.start_install()
    assert r["started"] is False and "already running" in r["error"]


@pytest.mark.wsl
def test_annotate_copia_family_live():
    seq = open("tests/fixtures/M11240.fasta").read()
    r = wsl.annotate(seq, species="drosophila melanogaster", timeout=400)
    assert r["ok"], r.get("error")
    assert any(h["class_family"] == "LTR/Copia" for h in r["hits"])
    assert any("Copia" in h["family"] for h in r["hits"])


@pytest.mark.wsl
def test_annotate_human_l1_family_live():
    seq = open("tests/fixtures/M80343.fasta").read()
    r = wsl.annotate(seq, species="Homo sapiens", timeout=400)
    assert r["ok"], r.get("error")
    assert any(h["class_family"] == "LINE/L1" for h in r["hits"])   # real Dfam family, 2nd organism
