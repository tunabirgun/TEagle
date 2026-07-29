"""Regression tests for species resolution before RepeatMasker (2026-07-28).

RepeatMasker delegates its lineage lookup to famdb, and famdb rejects an ambiguous name by printing its
own USAGE text and exiting non-zero — RepeatMasker then exits 255. The user's error message for
"drosophila" was therefore a wall of famdb help beginning
`% ./famdb.py lineage "Heterodontus japonicus" -ad -f semicolon`, which says nothing about the real
problem: the name matches 141 taxa and is not specific enough.

Live behaviour is covered by the @wsl-marked test; the parsing of famdb's replies is hermetic.
"""
import pytest

from teagle_core import wsl


class _FakeReply:
    """Stand in for _wsl_script so the message logic is testable without a backend."""

    def __init__(self, text):
        self.text = text

    def __call__(self, script, timeout=90):
        return 0, self.text, ""


def test_ambiguous_name_is_reported_as_ambiguity(monkeypatch):
    monkeypatch.setattr(wsl, "_wsl_script", _FakeReply(
        "Ambiguous search term 'drosophila' (found 141 results, 2 exact).\n"
        "Please use a more specific name or taxa ID"))
    r = wsl.resolve_species("drosophila")
    assert r["ok"] is False and r["ambiguous"] is True
    assert "141" in r["error"]
    assert "Drosophila melanogaster" in r["error"]      # tells the user what to type instead
    assert "famdb" not in r["error"].lower()            # never surface the tool's own usage text


def test_unknown_name_is_distinguished_from_ambiguous(monkeypatch):
    monkeypatch.setattr(wsl, "_wsl_script", _FakeReply("No species found matching that name"))
    r = wsl.resolve_species("Nonexistent taxon")
    assert r["ok"] is False and not r.get("ambiguous")
    assert "not found" in r["error"].lower()


def test_resolvable_name_passes(monkeypatch):
    monkeypatch.setattr(wsl, "_wsl_script", _FakeReply(
        "# Format: <NCBI tax ID> <scientific name>\n1 root [9]\n└─9606 Homo sapiens [1200]"))
    assert wsl.resolve_species("Homo sapiens")["ok"] is True


def test_invalid_token_is_rejected_without_touching_wsl(monkeypatch):
    def _boom(*a, **k):
        raise AssertionError("must not reach WSL for a malformed token")
    monkeypatch.setattr(wsl, "_wsl_script", _boom)
    assert wsl.resolve_species("../../etc/passwd; rm -rf /")["ok"] is False


def test_probe_failure_does_not_block_the_run(monkeypatch):
    """A broken probe must not turn into a refusal — let RepeatMasker try and report its own outcome."""
    def _raise(*a, **k):
        raise RuntimeError("wsl unavailable")
    monkeypatch.setattr(wsl, "_wsl_script", _raise)
    r = wsl.resolve_species("Homo sapiens")
    assert r["ok"] is True and "unchecked" in r


@pytest.mark.wsl
def test_species_scoped_annotate_succeeds_live():
    """The path that returned exit 255 for every species before this fix."""
    from helpers import fixture_seq
    r = wsl.annotate(">q\n" + fixture_seq("M80343"), species="Homo sapiens", timeout=600)
    assert r.get("ok") is True, r.get("error")


@pytest.mark.wsl
def test_ambiguous_species_is_rejected_live_with_guidance():
    r = wsl.annotate(">q\nACGT" * 100, species="drosophila", timeout=300)
    assert r.get("ok") is False
    assert r.get("ambiguous_species") is True
    assert "specific" in r["error"] or "full scientific name" in r["error"]
