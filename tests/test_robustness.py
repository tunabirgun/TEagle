"""Robustness / broken-environment tests (hermetic, mocked): no WSL, broken WSL,
malformed inputs. Verifies failures are explicit and never a silent success."""
import pytest
from teagle_core import wsl, sequtil, structural, domains, classify


# ---- WSL environment scenarios (mocked, no real WSL needed) ----

def test_no_wsl_at_all(monkeypatch):
    monkeypatch.setattr(wsl, "resolve_distro", lambda: None)
    av = wsl.available()
    assert av["wsl2"] is False and av["error"]
    r = wsl.annotate(">x\nACGTACGT", species=None)
    assert r["ok"] is False                              # explicit failure, not a fake result


def test_broken_wsl_stack(monkeypatch):
    # distro present, but RepeatMasker / Dfam not installed -> not ready, must not annotate
    monkeypatch.setattr(wsl, "resolve_distro", lambda: "Ubuntu-24.04")
    def fake_wsl(script, stdin=None, timeout=600):
        return (0, "ok\n", "") if "echo ok" in script else (0, "", "")
    monkeypatch.setattr(wsl, "_wsl", fake_wsl)
    st = wsl.env_status()
    assert st["wsl2"] is True and st["ready"] is False
    r = wsl.annotate(">x\nACGT")
    assert r["ok"] is False and "not ready" in r["error"]


def test_wsl_species_injection_blocked_without_env():
    # untrusted species is rejected before any WSL call (hermetic)
    r = wsl.annotate(">x\nACGT", species="drosophila; curl evil")
    assert r["ok"] is False and "invalid species" in r["error"]


# ---- malformed / degenerate inputs must not crash the core ----

@pytest.mark.parametrize("seq", ["", "N" * 50, "ACGT", "acgtACGTnnnn", ">only header\n", "X" * 5])
def test_core_handles_degenerate_inputs(seq):
    recs = sequtil.parse_fasta(seq)
    for _, s in recs:
        sequtil.composition(s)
        structural.detect_all(s)            # must not raise
        sequtil.find_orfs(s)


def test_domains_on_short_sequence_no_crash():
    assert domains.scan_domains("ACGT" * 5) == []      # too short for any ORF/domain


def test_classify_empty_evidence():
    c = classify.classify([], [])
    assert c["te_class"] == "none" and c["confidence"] == "Candidate"


def test_invalid_iupac_reported_not_silent():
    ok, bad = sequtil.validate_iupac("ACGT!@#Z")
    assert ok is False and len(bad) >= 3                # bad chars surfaced, not dropped
