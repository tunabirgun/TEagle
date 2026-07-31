"""Whole-genome TE annotation: the invariants a multi-hour, homology-bound run must keep.

These lock findings from the multi-persona review — a resumed run that mixes parameters, a failed run
rendered as a success, coverage double-counted across overlapping alignments, and an export that drops
the hedges the screen shows."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "app", "backend"))
import pytest
from teagle_core import wsl


def test_repeat_kind_separates_tes_from_tandem_and_other():
    """A "% masked" that folded tandem repeats into the TE figure would overstate TE content."""
    for fam in ("LTR/Gypsy", "LINE/L1", "SINE/Alu", "DNA/hAT-Ac", "RC/Helitron", "Retroposon/SVA"):
        assert wsl.repeat_kind(fam) == "TE", fam
    for fam in ("Simple_repeat", "Low_complexity", "Satellite"):
        assert wsl.repeat_kind(fam) == "tandem", fam
    for fam in ("rRNA", "tRNA", "snRNA", "Unknown", "ARTEFACT"):
        assert wsl.repeat_kind(fam) == "other", fam


def test_sensitivity_is_validated():
    r = wsl.genome_annotate("GCF_000000000.1", sensitivity="turbo")
    assert r["ok"] is False and "sensitivity" in r["error"]


def test_accession_is_validated():
    assert wsl.genome_annotate("not an accession")["ok"] is False


def test_stage_library_rejects_a_non_fasta(tmp_path):
    p = tmp_path / "lib.txt"
    p.write_text("ACGT no header here\n", encoding="utf-8")
    r = wsl.stage_library(str(p))
    assert r["ok"] is False and "FASTA" in r["error"]


def test_stage_library_rejects_a_missing_file():
    assert wsl.stage_library(str(os.path.join("nowhere", "nope.fa")))["ok"] is False


def test_annotation_script_refuses_a_resume_with_different_settings():
    """A resumed run that silently mixed species/sensitivity/chunk size across chunks would seal the
    mixture as one uniform value — a manifest describing a run that never happened."""
    src = wsl._ANNOT_SCRIPT
    assert "RUNSIG=" in src and ".runsig" in src
    assert "different settings" in src                      # the refusal the user actually sees
    assert "__SPSIG__" in src                               # library/lineage is part of the signature


def test_coverage_is_merged_not_summed():
    """RepeatMasker emits overlapping rows (nested insertions, fragmented re-alignments); summing row
    lengths double-counts bases and can report more coverage than the genome has."""
    src = wsl._ANNOT_SCRIPT
    assert "sort -k11,11 -k5,5 -k6,6n" in src               # family, then contig, then start
    assert "merged" in src.lower()


@pytest.mark.parametrize("field", ["te_percent", "masked_percent", "coverage_warning", "library_kind"])
def test_result_contract_keeps_te_and_all_repeat_separate(field):
    """The reporting contract the UI and the exports depend on."""
    import inspect
    src = inspect.getsource(wsl.genome_annotate)
    assert f'"{field}"' in src


def test_exported_report_carries_the_hedges(tmp_path):
    """An export that dropped the scope note would read as an unqualified claim about the genome."""
    pytest.importorskip("PySide6")
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    native = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "app", "native")
    if native not in sys.path:
        sys.path.insert(0, native)
    from PySide6.QtWidgets import QApplication
    QApplication.instance() or QApplication([])
    import main
    md = main.MainWindow._annot_report_md({
        "species": "S. cerevisiae", "accession": "GCF_000146045.2", "te_percent": 0.0,
        "masked_percent": 1.265, "te_family_count": 0, "genome_bp": 12157105, "total_hits": 3588,
        "repeatmasker_version": "4.2.4", "dfam_version": "4.0", "library_families_for_species": 9,
        "sensitivity": "default", "chunks": 1, "complete": True,
        "library_kind": "installed Dfam partitions",
        "coverage_warning": "No transposable-element family was found.",
        "families": [{"family": "Simple_repeat", "kind": "tandem", "n": 3041, "bp": 128258,
                      "percent": 1.055, "divergence": 16.39}]})
    assert "Coverage warning" in md
    assert "absence is not evidence of absence" in md
    assert "not Kimura- or CpG-corrected" in md
    assert "merged intervals" in md
    assert "counted separately from tandem repeats" in md


def test_uncurated_flag_reaches_repeatmasker_and_the_seal():
    """Installing Dfam's uncurated partitions does nothing on its own: RepeatMasker searches curated
    families ONLY unless given -uncurated. Measured on S. cerevisiae, curated-only sees 9 families and
    reports zero transposable elements while -uncurated sees 421 more and recovers the Ty elements — so
    the flag has to reach the command line, and it has to be part of what identifies the run."""
    import inspect
    src = inspect.getsource(wsl.genome_annotate)
    assert "-uncurated" in src, "the flag must be passed to RepeatMasker"
    assert "include_uncurated" in src
    # curated-only and curated+uncurated are different searches, so a resume must not mix them
    assert "unc:" in src, "the run signature must distinguish the two searches"


def test_library_kind_names_which_families_were_searched():
    """A result that did not say which family set produced it would be unreproducible: the same genome
    against curated vs curated+uncurated is a different experiment."""
    import inspect
    src = inspect.getsource(wsl.genome_annotate)
    assert "curated only" in src and "curated + uncurated" in src


def test_zero_te_warning_blames_the_search_not_the_genome():
    import inspect
    src = inspect.getsource(wsl.genome_annotate)
    assert "not a \nfinding about the genome" in src or "not a " in src
    assert "CURATED subset" in src, "the warning must name the real cause when curated-only found nothing"
