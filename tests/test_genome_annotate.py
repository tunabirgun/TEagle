"""Whole-genome TE annotation: the invariants a multi-hour, homology-bound run must keep.

These lock findings from the multi-persona review — a resumed run that mixes parameters, a failed run
rendered as a success, coverage double-counted across overlapping alignments, and an export that drops
the hedges the screen shows."""
import hashlib
import os
import re
import shutil
import subprocess
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "app", "backend"))
import pytest
from teagle_core import provenance, wsl


def _bash_path():
    """Prefer Git Bash on Windows; System32/bash.exe is only a WSL launcher."""
    candidates = []
    if os.name == "nt":
        candidates.extend((r"C:\Program Files\Git\bin\bash.exe",
                           r"C:\Program Files (x86)\Git\bin\bash.exe"))
    candidates.append(shutil.which("bash"))
    path = next((p for p in candidates if p and os.path.isfile(p)), None)
    if not path:
        pytest.skip("bash is required to execute the production WSL scripts")
    return path


def _run_bash(script, cwd, *, home=None):
    env = os.environ.copy()
    if home is not None:
        env["HOME"] = os.fspath(home).replace("\\", "/")
    return subprocess.run([_bash_path(), "-s"], input=script, text=True, cwd=cwd,
                          env=env, capture_output=True, timeout=60)


def _coverage_block():
    match = re.search(r"# BEGIN_COVERAGE_SUMMARY.*?# END_COVERAGE_SUMMARY",
                      wsl._ANNOT_SCRIPT, flags=re.S)
    assert match, "production coverage block markers are missing"
    return match.group(0)


_OVERLAPPING_OUT = """\
100 1.0 0.0 0.0 chr1 1 50 (50) + teA LTR/Gypsy
100 1.5 0.0 0.0 chr1 20 30 (70) + teC LTR/Gypsy
100 2.0 0.0 0.0 chr1 41 80 (20) + teB LINE/L1
100 2.5 0.0 0.0 chr1 75 90 (10) + pen PLE/Penelope
100 3.0 0.0 0.0 chr1 70 100 (0) + low Simple_repeat
100 4.0 0.0 0.0 chr2 1 10 (0) + rib rRNA
"""


def _parse_coverage(text):
    parsed = {"repeat_names": {}, "kinds": {}}
    for line in text.splitlines():
        if line.startswith("REP\t"):
            _, repeat_name, class_family, _n, bp, _div = line.split("\t")
            parsed["repeat_names"][repeat_name] = {"class_family": class_family, "bp": int(bp)}
        elif line.startswith("KIND\t"):
            _, kind, bp = line.split("\t")
            parsed["kinds"][kind] = int(bp)
        elif line.startswith(("TOTAL_BP ", "TOTAL_HITS ", "GENOME_BP ")):
            key, value = line.split()
            parsed[key.lower()] = int(value)
    return parsed


def _assert_overlap_summary(summary):
    assert summary["repeat_names"] == {
        "teA": {"class_family": "LTR/Gypsy", "bp": 50},
        "teC": {"class_family": "LTR/Gypsy", "bp": 11},
        "teB": {"class_family": "LINE/L1", "bp": 40},
        "pen": {"class_family": "PLE/Penelope", "bp": 16},
        "low": {"class_family": "Simple_repeat", "bp": 31},
        "rib": {"class_family": "rRNA", "bp": 10},
    }
    assert summary["total_hits"] == 6
    assert summary["total_bp"] == 110      # chr1 1..100 plus chr2 1..10, regardless of repeat name/model
    assert summary["kinds"] == {"TE": 90, "tandem": 31, "other": 10}


def test_repeat_kind_separates_tes_from_tandem_and_other():
    """A "% masked" that folded tandem repeats into the TE figure would overstate TE content."""
    for fam in ("LTR/Gypsy", "LINE/L1", "SINE/Alu", "DNA/hAT-Ac", "RC/Helitron", "Retroposon/SVA",
                "PLE/Penelope"):
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
    assert "genome_sha256=$GENOME_SHA" in src               # input identity is part of the signature too


def _resume_fixture(tmp_path, sequence="A" * 200):
    accession = "GCF_000000000.1"
    genome_dir = tmp_path / "teagle_genomes" / accession
    annot_dir = genome_dir / "annot"
    annot_dir.mkdir(parents=True)
    fasta = ">chr1\n" + sequence + "\n"
    sha = hashlib.sha256(fasta.encode()).hexdigest()
    (genome_dir / "genome.fna").write_bytes(fasta.encode())
    (genome_dir / "meta.txt").write_text(
        f"accession={accession}\ntarget=genome.fna\nsha256={sha}\n"
        f"target_sha256={sha}\nn_seqs=1\nbytes={len(fasta.encode())}\n", encoding="utf-8")
    (genome_dir / ".done").touch()
    (annot_dir / ".split_done").touch()
    output = annot_dir / "chunk_0001.out"
    output.write_text(_OVERLAPPING_OUT, encoding="utf-8")
    (annot_dir / "chunk_0001.out.sha256").write_text(
        hashlib.sha256(output.read_bytes()).hexdigest() + "\n", encoding="utf-8")
    (annot_dir / ".runsig").write_text(
        f"genome_sha256={sha}|species=test-library|sens=|chunkbp=1000000", encoding="utf-8")
    script = (wsl._ANNOT_SCRIPT.replace("__ACC__", accession).replace("__PA__", "1")
              .replace("__CHUNKBP__", "1000000").replace("__SENS__", "")
              .replace("__SPSIG__", "test-library").replace("__SP__", ""))
    return genome_dir, annot_dir, sha, script


def test_annotation_rehashes_cached_fasta_before_reusing_chunks(tmp_path):
    genome_dir, _annot_dir, _sha, script = _resume_fixture(tmp_path)
    intact = _run_bash(script, tmp_path, home=tmp_path)
    assert intact.returncode == 0, intact.stdout + intact.stderr
    assert "ANNOT_OK" in intact.stdout

    (genome_dir / "genome.fna").write_bytes((">chr1\n" + "C" * 200 + "\n").encode())
    corrupt = _run_bash(script, tmp_path, home=tmp_path)
    assert corrupt.returncode != 0
    assert "cached genome FASTA failed its sealed checksum" in corrupt.stdout
    assert "summarising" not in corrupt.stdout


def test_annotation_resume_signature_rejects_a_changed_verified_genome(tmp_path):
    genome_dir, annot_dir, old_sha, script = _resume_fixture(tmp_path)
    changed = ">chr1\n" + "C" * 200 + "\n"
    new_sha = hashlib.sha256(changed.encode()).hexdigest()
    (genome_dir / "genome.fna").write_bytes(changed.encode())
    meta = (genome_dir / "meta.txt").read_text(encoding="utf-8")
    (genome_dir / "meta.txt").write_text(meta.replace(old_sha, new_sha), encoding="utf-8")

    resumed = _run_bash(script, tmp_path, home=tmp_path)
    assert resumed.returncode != 0
    assert "part-finished annotation run with different settings" in resumed.stdout
    assert old_sha in (annot_dir / ".runsig").read_text(encoding="utf-8")


def test_annotation_resume_rejects_a_mutated_completed_output(tmp_path):
    _genome_dir, annot_dir, _sha, script = _resume_fixture(tmp_path)
    with (annot_dir / "chunk_0001.out").open("a", encoding="utf-8") as handle:
        handle.write("100 0 0 0 chr1 1 2 (0) + injected LTR/Gypsy\n")
    resumed = _run_bash(script, tmp_path, home=tmp_path)
    assert resumed.returncode != 0
    assert "completed annotation output failed its integrity seal" in resumed.stdout
    assert "summarising" not in resumed.stdout


def test_annotation_resume_reconciles_output_committed_before_input_cleanup(tmp_path):
    """A crash between output commit and input deletion must not count one chunk twice."""
    genome_dir, annot_dir, _sha, script = _resume_fixture(tmp_path)
    chunk = annot_dir / "chunk_0001.fa"
    chunk.write_bytes((genome_dir / "genome.fna").read_bytes())
    input_seal = annot_dir / "chunk_0001.input.sha256"
    input_seal.write_text(hashlib.sha256(chunk.read_bytes()).hexdigest() + "\n", encoding="utf-8")

    resumed = _run_bash(script, tmp_path, home=tmp_path)
    assert resumed.returncode == 0, resumed.stdout + resumed.stderr
    assert "chunks: 1 total, 1 already finished" in resumed.stdout
    assert "TOTAL_CHUNKS 1" in resumed.stdout
    assert not chunk.exists()
    assert not input_seal.exists()


def _pending_annotation_fixture(tmp_path):
    genome_dir, annot_dir, sha, script = _resume_fixture(tmp_path)
    (annot_dir / "chunk_0001.out").unlink()
    (annot_dir / "chunk_0001.out.sha256").unlink()
    chunk = annot_dir / "chunk_0001.fa"
    chunk.write_bytes((genome_dir / "genome.fna").read_bytes())
    (annot_dir / "chunk_0001.input.sha256").write_text(
        hashlib.sha256(chunk.read_bytes()).hexdigest() + "\n", encoding="utf-8")
    return genome_dir, annot_dir, chunk, sha, script


def test_annotation_resume_rejects_a_mutated_pending_chunk_before_repeatmasker(tmp_path):
    _genome_dir, _annot_dir, chunk, _sha, script = _pending_annotation_fixture(tmp_path)
    chunk.write_bytes(b">chr1\nALTERED\n")
    called = tmp_path / "repeatmasker-called"
    binary = tmp_path / "bin" / "micromamba"
    binary.parent.mkdir()
    binary.write_text(f"#!/usr/bin/env bash\ntouch '{called.as_posix()}'\nexit 0\n", encoding="utf-8")
    os.chmod(binary, 0o755)

    resumed = _run_bash(script, tmp_path, home=tmp_path)
    assert resumed.returncode != 0
    assert "pending work chunk failed its integrity seal" in resumed.stdout
    assert not called.exists()


def test_annotation_final_genome_rehash_rejects_a_midrun_mutation(tmp_path):
    genome_dir, _annot_dir, _chunk, _sha, script = _pending_annotation_fixture(tmp_path)
    binary = tmp_path / "bin" / "micromamba"
    binary.parent.mkdir()
    binary.write_text(
        "#!/usr/bin/env bash\n"
        "chunk=\"${@: -1}\"\n"
        "printf '' > \"$chunk.out\"\n"
        f"printf '>chr1\\n{'C' * 200}\\n' > '{(genome_dir / 'genome.fna').as_posix()}'\n"
        "exit 0\n", encoding="utf-8")
    os.chmod(binary, 0o755)

    run = _run_bash(script, tmp_path, home=tmp_path)
    assert run.returncode != 0
    assert "genome FASTA changed during annotation" in run.stdout
    assert "ANNOT_OK" not in run.stdout


def test_annotation_integrity_gate_rejects_missing_output_seal(tmp_path):
    _genome_dir, annot_dir, _sha, script = _resume_fixture(tmp_path)
    (annot_dir / "chunk_0001.out.sha256").unlink()
    run = _run_bash(script, tmp_path, home=tmp_path)
    assert run.returncode != 0
    assert "completed annotation output lacks its integrity seal" in run.stdout


def test_prepare_rehashes_an_existing_cached_target(tmp_path):
    accession = "GCF_000000000.1"
    genome_dir = tmp_path / "teagle_genomes" / accession
    genome_dir.mkdir(parents=True)
    fasta = ">chr1\nACGT\n"
    sha = hashlib.sha256(fasta.encode()).hexdigest()
    (genome_dir / "genome.fna").write_bytes(fasta.encode())
    # Legacy 3.8.0 FASTA cache: no target_sha256, so the source-FASTA seal is the valid fallback.
    (genome_dir / "meta.txt").write_text(
        f"accession={accession}\ntarget=genome.fna\nsha256={sha}\nn_seqs=1\nbytes={len(fasta)}\n",
        encoding="utf-8")
    (genome_dir / ".done").touch()
    script = wsl._PREP_SCRIPT.replace("__ACC__", accession)
    intact = _run_bash(script, tmp_path, home=tmp_path)
    assert intact.returncode == 0 and "PREPARED" in intact.stdout

    (genome_dir / "genome.fna").write_bytes(b">chr1\nTGCA\n")
    corrupt = _run_bash(script, tmp_path, home=tmp_path)
    assert corrupt.returncode != 0
    assert "failed its sealed checksum" in corrupt.stdout


@pytest.mark.parametrize("target,checksum,body,expected", [
    ("unexpected.bin", "a" * 64, b"x", "invalid target"),
    ("genome.2bit", "", b"x", "no sealed target checksum"),
    ("genome.2bit", "a" * 64, b"changed", "failed its sealed checksum"),
])
def test_prepare_cache_failures_are_classified_before_download_tools(tmp_path, target, checksum, body, expected):
    accession = "GCF_000000000.1"
    genome_dir = tmp_path / "teagle_genomes" / accession
    genome_dir.mkdir(parents=True)
    (genome_dir / target).write_bytes(body)
    (genome_dir / "meta.txt").write_text(
        f"accession={accession}\ntarget={target}\nsha256={'b' * 64}\ntarget_sha256={checksum}\n",
        encoding="utf-8")
    (genome_dir / ".done").touch()
    called = tmp_path / "datasets-called"
    datasets = tmp_path / "micromamba" / "envs" / "te" / "bin" / "datasets"
    datasets.parent.mkdir(parents=True)
    datasets.write_text(f"#!/usr/bin/env bash\ntouch '{called.as_posix()}'\nexit 1\n", encoding="utf-8")
    os.chmod(datasets, 0o755)
    run = _run_bash(wsl._PREP_SCRIPT.replace("__ACC__", accession), tmp_path, home=tmp_path)
    assert run.returncode != 0 and "CACHE_CORRUPT" in run.stdout and expected in run.stdout
    assert not called.exists()


@pytest.mark.parametrize("detail", [
    "cached genome metadata has an invalid target",
    "cached genome has no sealed target checksum",
    "cached genome target failed its sealed checksum",
])
def test_prepare_cache_failure_returns_recoverable_contract(monkeypatch, detail):
    monkeypatch.setattr(wsl, "available", lambda: {"wsl2": True})
    monkeypatch.setattr(wsl, "_wsl", lambda *a, **k: (0, "yes\n", ""))
    monkeypatch.setattr(wsl, "_wsl_script",
                        lambda *a, **k: (1, f"FAIL CACHE_CORRUPT {detail}\n", ""))
    result = wsl.genome_prepare("GCF_000000000.1")
    assert result["ok"] is False and result["cache_corrupt"] is True
    assert "Manage genomes" in result["error"] and "remove" in result["error"] and "download" in result["error"]


def test_genome_scan_stops_on_cache_integrity_failure(monkeypatch):
    calls = []
    monkeypatch.setattr(wsl, "available", lambda: {"wsl2": True})
    monkeypatch.setattr(wsl, "_wsl_script",
                        lambda script, **kw: (calls.append(script) or (0, "CACHE_BAD checksum mismatch\n", "")))
    result = wsl.genome_scan("GCF_000000000.1", "q\tACGT\tTGCA\n")
    assert result["ok"] is False and result["cache_corrupt"] is True
    assert "integrity" in result["error"] and "remove it" in result["error"] and len(calls) == 1


def test_genome_scan_executes_the_production_cache_probe_before_ispcr(tmp_path, monkeypatch):
    """Run the exact probe generated by genome_scan against an intact then mutated cache.

    This is a second gate over the mocked CACHE_BAD contract above: it proves the production shell reads
    and hashes the cached target, and proves a mutation returns a distinct corrupt-cache response before
    isPcr is staged or run.
    """
    accession = "GCF_000000000.1"
    genome_dir = tmp_path / "teagle_genomes" / accession
    genome_dir.mkdir(parents=True)
    fasta = b">chr1\nACGTACGT\n"
    sha = hashlib.sha256(fasta).hexdigest()
    (genome_dir / "genome.fna").write_bytes(fasta)
    (genome_dir / "meta.txt").write_text(
        f"accession={accession}\ntarget=genome.fna\nsha256={sha}\n"
        f"target_sha256={sha}\nn_seqs=1\nbytes={len(fasta)}\n", encoding="utf-8")
    (genome_dir / ".done").touch()

    probes, downstream = [], []
    monkeypatch.setattr(wsl, "available", lambda: {"wsl2": True})

    def execute_probe(script, **_kwargs):
        if "echo CACHE_OK" in script:
            probes.append(script)
            returncode = _run_bash(script, tmp_path, home=tmp_path)
            return returncode.returncode, returncode.stdout, returncode.stderr
        downstream.append(script)
        return 8, "", ""

    monkeypatch.setattr(wsl, "_wsl_script", execute_probe)
    monkeypatch.setattr(wsl, "_wsl", lambda *a, **k: (downstream.append(a[0]) or (0, "", "")))

    intact = wsl.genome_scan(accession, "q\tACGT\tACGT\n")
    assert intact["ok"] is False                         # the fake downstream cannot execute real isPcr
    assert downstream, "intact cache never reached query staging or isPcr"

    downstream.clear()
    (genome_dir / "genome.fna").write_bytes(b">chr1\nTGCATGCA\n")
    mutated = wsl.genome_scan(accession, "q\tACGT\tACGT\n")
    assert mutated["ok"] is False and mutated["cache_corrupt"] is True
    assert "integrity" in mutated["error"]
    assert downstream == [], "mutated cache reached query staging or isPcr"
    assert len(probes) == 2


def test_generated_genome_scripts_are_valid_bash(tmp_path):
    scripts = [wsl._PREP_SCRIPT.replace("__ACC__", "GCF_000000000.1"),
               (wsl._ANNOT_SCRIPT.replace("__ACC__", "GCF_000000000.1").replace("__PA__", "1")
                .replace("__CHUNKBP__", "1000000").replace("__SENS__", "")
                .replace("__SPSIG__", "test-library").replace("__SP__", ""))]
    for script in scripts:
        check = subprocess.run([_bash_path(), "-n", "-s"], input=script, text=True,
                               capture_output=True, cwd=tmp_path, timeout=60)
        assert check.returncode == 0, check.stderr


def test_coverage_uses_cross_family_and_kind_unions(tmp_path):
    """Execute the production summary on cross-family and TE/non-TE overlaps.

    Repeat-name/model rows retain their own unions, but the all-repeat and TE headlines must each be a contig-wise
    union. Therefore the rows intentionally sum to 158 bp while only 110 distinct genomic bases are
    masked, and the tandem overlap cannot inflate the 90-bp TE-only union.
    """
    (tmp_path / "all.out").write_text(_OVERLAPPING_OUT, encoding="utf-8")
    run = _run_bash("set -euo pipefail\nGENOME_BP=200\n" + _coverage_block(), tmp_path)
    assert run.returncode == 0, run.stderr
    _assert_overlap_summary(_parse_coverage(run.stdout))


def test_coverage_gate_rejects_the_old_per_row_family_sum(tmp_path):
    """Inject the 3.8.0 defect and watch the numerical gate fail.

    In this corpus no family overlaps itself, so summing raw rows exactly reproduces the former sum of
    per-name/model unions: 158 bp rather than the 110-bp cross-model union.
    """
    (tmp_path / "all.out").write_text(_OVERLAPPING_OUT, encoding="utf-8")
    block = _coverage_block()
    needle = 'printf "TOTAL_BP %d\\n" "$TOTAL_BP"'
    mutant = ('TOTAL_BP=$(awk \'$1 ~ /^[0-9]+$/ { total += (($7+0)-($6+0)+1) } '
              'END { print total+0 }\' all.out)\n' + needle)
    assert needle in block
    run = _run_bash("set -euo pipefail\nGENOME_BP=200\n" + block.replace(needle, mutant), tmp_path)
    assert run.returncode == 0, run.stderr
    assert _parse_coverage(run.stdout)["total_bp"] == 158
    with pytest.raises(AssertionError):
        _assert_overlap_summary(_parse_coverage(run.stdout))


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
        "masked_percent": 1.265, "te_repeat_name_count": 0, "genome_bp": 12157105, "total_hits": 3588,
        "repeatmasker_version": "4.2.4", "dfam_version": "4.0", "library_families_for_species": 9,
        "sensitivity": "default", "chunks": 1, "complete": True,
        "library_kind": "installed Dfam partitions",
        "coverage_warning": "No TE-classified repeat name/model was found.",
        "repeat_names": [{"repeat_name": "(TA)n", "class_family": "Simple_repeat", "kind": "tandem",
                          "n": 3041, "bp": 128258,
                          "percent": 1.055, "divergence": 16.39}]})
    assert "Coverage warning" in md
    assert "absence is not evidence of absence" in md
    assert "not Kimura- or CpG-corrected" in md
    assert "within-name interval union" in md
    assert "separate cross-model unions" in md
    assert "Summing rows can exceed a headline" in md
    assert "independently unioned and non-additive" in md
    assert "including ambiguous bases" in md
    assert "alignment-span-weighted mean RepeatMasker divergence" in md
    assert "counted separately from tandem repeats" in md


def _assert_custom_library_copy(texts, name, sha):
    joined = "\n".join(texts)
    assert name in joined
    assert sha in joined
    assert "dfam" not in joined.lower()
    assert "family models available" not in joined.lower()
    assert "installed partitions" not in joined.lower()


def test_custom_library_native_copy_and_report_name_only_the_supplied_fasta():
    pytest.importorskip("PySide6")
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    native = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "app", "native")
    if native not in sys.path:
        sys.path.insert(0, native)
    from PySide6.QtWidgets import QApplication, QLabel
    app = QApplication.instance() or QApplication([])
    import main

    name = "curated-lab-library.fa"
    sha = "9f" * 32
    result = {
        "species": "Example species", "accession": "GCF_000000001.1", "te_percent": 1.25,
        "masked_percent": 1.5, "te_repeat_name_count": 1, "genome_bp": 1_000_000, "total_hits": 2,
        "repeatmasker_version": "4.2.4", "dfam_version": "4.0", "library_families_for_species": 1439,
        "sensitivity": "default", "chunks": 1, "complete": True, "library_kind": "user",
        "custom_library": {"name": name, "sha256": sha, "sequences": 7},
        "repeat_names": [{"repeat_name": "LabTE", "class_family": "LTR/Gypsy", "kind": "TE", "n": 2, "bp": 100,
                          "percent": 0.01, "divergence": 1.2}],
    }
    context = main.MainWindow._annot_library_context(result)
    _assert_custom_library_copy(
        [context["display"], context["tooltip"], context["scope"], context["report_scope"]], name, sha)

    report = main.MainWindow._annot_report_md(result)
    _assert_custom_library_copy([report], name, sha)
    assert "7 reference sequences" in report
    assert "Bases (within-name union)" in report
    assert "Repeat name/model" in report and "distinct repeat names/models" in report
    assert "Class/family" in report and "LTR/Gypsy" in report

    window = main.MainWindow()
    try:
        result["elapsed_s"] = 1.2
        window._show_annot_window(result)
        rendered = "\n".join(
            label.text() + "\n" + label.toolTip()
            for label in window._annot_window.findChildren(QLabel)
        )
        _assert_custom_library_copy([rendered], name, sha)
        assert "separate cross-model interval union" in rendered
        assert "separate contig-wise union of every repeat row" in rendered
        assert "that plus" not in rendered.lower()
        copies = dict(main.MainWindow._ANNOT_COLS)["Copies"]
        assert "alignment rows" in copies and "fragmented insertion" in copies
    finally:
        if getattr(window, "_annot_window", None) is not None:
            window._annot_window.close()
        window.close()


def test_custom_library_copy_gate_rejects_injected_dfam_claim():
    """The reporting gate must fail if an old Dfam/model-count sentence is reintroduced."""
    name = "lab.fa"
    sha = "ab" * 32
    valid = [f"User-supplied FASTA {name}; SHA-256 {sha}"]
    _assert_custom_library_copy(valid, name, sha)
    with pytest.raises(AssertionError):
        _assert_custom_library_copy(
            valid + ["Dfam 4.0; 1439 family models available for this lineage"], name, sha)


def test_uncurated_flag_reaches_repeatmasker_and_the_seal():
    """Installing Dfam's uncurated partitions does nothing on its own: RepeatMasker searches curated
    families ONLY unless given -uncurated. Measured on S. cerevisiae, curated-only sees 9 families and
    reports zero transposable elements while -uncurated adds the remaining models and recovers Ty elements — so
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


def _annotation_output(genome_sha="a" * 64):
    return ("TOTAL_CHUNKS 1\n" + f"GENOME_SHA256 {genome_sha}\n"
            "GENOME_BP 200\nTOTAL_HITS 4\nTOTAL_BP 110\n"
            "KIND\tTE\t80\nKIND\ttandem\t31\nKIND\tother\t10\n"
            "REP\tteA\tLTR/Gypsy\t1\t50\t1.00\n"
            "REP\tteB\tLINE/L1\t1\t40\t2.00\n"
            "REP\tlow\tSimple_repeat\t1\t31\t3.00\n"
            "REP\trib\trRNA\t1\t10\t4.00\nANNOT_OK\n")


def test_custom_library_does_not_require_or_bind_unused_dfam(monkeypatch):
    scripts = []
    staged = {"ok": True, "path": "$HOME/teagle_libs/abc.fa", "sha256": "c" * 64,
              "sequences": 2, "bytes": 20, "name": "custom.fa"}
    monkeypatch.setattr(wsl, "stage_library", lambda _path: dict(staged))
    states = iter((
        {"ready": False, "repeatmasker": "4.2.4", "dfam": False,
         "dfam_library": {"version": "4.0", "partitions": ["old"]}},
        {"ready": True, "repeatmasker": "4.2.4", "dfam": True,
         "dfam_library": {"version": "4.1", "partitions": ["new", "extra"]}},
    ))
    monkeypatch.setattr(wsl, "env_status", lambda: next(states))
    monkeypatch.setattr(wsl, "_wsl_script",
                        lambda script, **kw: (scripts.append(script) or (0, _annotation_output(), "")))

    first = wsl.genome_annotate("GCF_000000000.1", custom_lib="custom.fa", include_uncurated=True)
    second = wsl.genome_annotate("GCF_000000000.1", custom_lib="custom.fa", include_uncurated=True)
    assert first["ok"] and second["ok"]
    assert first["dfam_library"] is None and first["dfam_version"] is None
    assert first["include_uncurated"] is False and first["species"] is None
    assert scripts[0] == scripts[1], "an unused Dfam mutation changed the custom-library run signature"
    assert '-lib "$HOME/teagle_libs/abc.fa"' in scripts[0]
    runsig = next(line for line in scripts[0].splitlines() if line.startswith("RUNSIG="))
    assert "lib:" + "c" * 64 in runsig and "|parts:" not in runsig and "|unc:" not in runsig


def test_dfam_backed_genome_annotation_still_requires_dfam(monkeypatch):
    monkeypatch.setattr(wsl, "env_status", lambda: {
        "ready": False, "repeatmasker": "4.2.4", "dfam": False, "dfam_library": None})
    called = []
    monkeypatch.setattr(wsl, "_wsl_script", lambda *a, **k: (called.append(True) or (0, "", "")))
    result = wsl.genome_annotate("GCF_000000000.1")
    assert result["ok"] is False and "not ready" in result["error"]
    assert called == []


def test_budget_custom_skips_dfam_and_selected_search_set_controls_famdb(monkeypatch):
    """The pre-run estimate must count the exact library selected, and no unrelated library for custom."""
    import engine
    monkeypatch.setattr(engine.wsl, "annotate_budget", lambda _n: {"ok": True, "recommended_threads": 2})
    calls = []
    monkeypatch.setattr(engine.wsl, "dfam_lineage_families",
                        lambda species, include_uncurated=False:
                        (calls.append((species, include_uncurated)) or {"ok": True, "families": 9}))

    custom = engine.run_annotate_budget({"species": "Saccharomyces cerevisiae",
                                         "library_mode": "custom", "budget_request_id": 1})
    curated = engine.run_annotate_budget({"species": "Saccharomyces cerevisiae",
                                          "library_mode": "dfam_curated", "budget_request_id": 2})
    all_families = engine.run_annotate_budget({"species": "Saccharomyces cerevisiae",
                                               "library_mode": "dfam_uncurated", "budget_request_id": 3})
    assert "library_families_for_species" not in custom
    assert calls == [("Saccharomyces cerevisiae", False),
                     ("Saccharomyces cerevisiae", True)]
    assert curated["library_mode"] == "dfam_curated" and curated["budget_request_id"] == 2
    assert all_families["library_mode"] == "dfam_uncurated" and all_families["budget_request_id"] == 3


def test_uncurated_budget_reports_missing_partitions_instead_of_claiming_the_search(monkeypatch):
    import engine
    monkeypatch.setattr(engine.wsl, "annotate_budget", lambda _n: {"ok": True, "recommended_threads": 2})
    monkeypatch.setattr(engine.wsl, "dfam_lineage_families", lambda *a, **k: {
        "ok": False, "uncurated_available": False,
        "error": "the optional uncurated Dfam partitions are not installed",
    })
    result = engine.run_annotate_budget({"species": "Saccharomyces cerevisiae",
                                         "library_mode": "dfam_uncurated", "budget_request_id": 9})
    assert result["uncurated_available"] is False
    assert result["library_ready"] is False
    assert result["library_families_for_species"] is None
    assert "not installed" in result["library_error"]


def test_uncurated_budget_without_lineage_blocks_start_before_famdb(monkeypatch):
    import engine
    monkeypatch.setattr(engine.wsl, "annotate_budget", lambda _n: {"ok": True, "recommended_threads": 2})
    calls = []
    monkeypatch.setattr(engine.wsl, "dfam_lineage_families",
                        lambda *a, **k: calls.append((a, k)))
    result = engine.run_annotate_budget({"species": None, "library_mode": "dfam_uncurated",
                                         "budget_request_id": 10})
    assert result["library_ready"] is False
    assert "select an organism" in result["library_error"]
    assert calls == []


def test_genome_uncurated_without_lineage_fails_before_annotation(monkeypatch):
    monkeypatch.setattr(wsl, "env_status", lambda: {
        "ready": True, "repeatmasker": "4.2.4", "dfam": True, "dfam_uncurated": True,
        "dfam_library": {"version": "4.0", "partitions": ["curated", "uncurated"]},
    })
    calls = []
    monkeypatch.setattr(wsl, "_wsl_script", lambda *a, **k: (calls.append(True) or (0, "", "")))
    result = wsl.genome_annotate("GCF_000000000.1", species=None, include_uncurated=True)
    assert result["ok"] is False and "select an organism" in result["error"]
    assert calls == []


@pytest.mark.parametrize("detail", [
    "prepared genome has no sealed source FASTA checksum",
    "cached genome FASTA failed its sealed checksum",
    "genome FASTA changed during annotation; no result was sealed",
])
def test_annotation_cache_failure_returns_recoverable_contract(monkeypatch, detail):
    monkeypatch.setattr(wsl, "env_status", lambda: {
        "ready": True, "repeatmasker": "4.2.4", "dfam": True, "dfam_uncurated": False,
        "dfam_library": {"version": "4.0", "partitions": ["curated"]},
    })
    monkeypatch.setattr(wsl, "_wsl_script",
                        lambda *a, **k: (1, f"FAIL CACHE_CORRUPT {detail}\n", ""))
    result = wsl.genome_annotate("GCF_000000000.1")
    assert result["ok"] is False and result["cache_corrupt"] is True and result["partial"] is True
    assert "Manage genomes" in result["error"] and "remove" in result["error"] and "download" in result["error"]


def test_famdb_count_command_matches_the_selected_search_set(monkeypatch):
    scripts = []
    monkeypatch.setattr(wsl, "_wsl_script",
                        lambda script, **_kwargs: (scripts.append(script) or (0, "9\n", "")))
    assert wsl.dfam_lineage_families("Saccharomyces cerevisiae", False)["ok"]
    assert wsl.dfam_lineage_families("Saccharomyces cerevisiae", True)["ok"]
    assert " families -a -c -f summary " in scripts[0]
    assert " families -a -f summary " in scripts[1]
    assert " -c " not in scripts[1].split("|", 1)[0]
    assert "dfam40.uncurated.consensus.0.h5" in scripts[1]
    assert "dfam40.uncurated.consensus.1.h5" in scripts[1]


def test_uncurated_lineage_count_fails_closed_when_partitions_are_absent(monkeypatch):
    monkeypatch.setattr(wsl, "_wsl_script",
                        lambda script, **_kwargs: (3, "TEAGLE_UNCURATED_MISSING\n", ""))
    result = wsl.dfam_lineage_families("Saccharomyces cerevisiae", True)
    assert result["ok"] is False and result["uncurated_available"] is False
    assert "not installed" in result["error"]


def test_dfam_annotation_rejects_requested_uncurated_set_when_not_installed(monkeypatch):
    monkeypatch.setattr(wsl, "resolve_species", lambda _species: {"ok": True})
    monkeypatch.setattr(wsl, "env_status", lambda: {
        "ready": True, "repeatmasker": "4.2.4", "dfam": True, "dfam_uncurated": False,
        "dfam_library": {"version": "4.0", "partitions": ["root", "curated"]},
    })
    called = []
    monkeypatch.setattr(wsl, "_wsl_script", lambda *a, **k: (called.append(True) or (0, "", "")))
    result = wsl.genome_annotate("GCF_000000000.1", species="Saccharomyces cerevisiae",
                                 include_uncurated=True)
    assert result["ok"] is False and result["uncurated_available"] is False
    assert "not installed" in result["error"]
    assert called == []


def test_selected_search_set_count_gate_rejects_a_mutated_curated_command(monkeypatch):
    """Inject omission of -c and watch the exact-command gate reject the old all-family count."""
    scripts = []
    monkeypatch.setattr(wsl, "_wsl_script",
                        lambda script, **_kwargs: (scripts.append(script) or (0, "9\n", "")))
    wsl.dfam_lineage_families("Saccharomyces cerevisiae", False)
    assert " families -a -c -f summary " in scripts[0]
    mutant = scripts[0].replace(" -a -c -f ", " -a -f ")
    with pytest.raises(AssertionError):
        assert " families -a -c -f summary " in mutant


def _engine_annotation(custom=True, genome_sha="a" * 64):
    base = {"ok": True, "accession": "GCF_000000000.1", "repeat_names": [],
            "genome_sha256": genome_sha, "genome_bp": 200, "chunks": 1, "complete": True,
            "total_hits": 0, "total_bp": 0, "te_bp": 0, "repeatmasker_version": "4.2.4",
            "include_uncurated": False, "sensitivity": "default", "threads": 1, "chunk_mb": 40}
    if custom:
        base.update({"custom_library": {"name": "custom.fa", "sha256": "c" * 64,
                                        "sequences": 2, "bytes": 20},
                     "library_kind": "user-supplied FASTA library", "dfam_library": None,
                     "dfam_version": None, "species": None})
    else:
        base.update({"custom_library": None, "library_kind": "installed Dfam partitions, curated only",
                     "dfam_library": {"version": "4.0", "partitions": ["0 [root]: curated"]},
                     "dfam_version": "4.0", "species": "Homo sapiens"})
    return base


def test_genome_annotation_preserves_repeat_names_separately_from_class_labels(monkeypatch):
    staged = {"ok": True, "path": "$HOME/teagle_libs/abc.fa", "sha256": "c" * 64,
              "sequences": 2, "bytes": 20, "name": "custom.fa"}
    monkeypatch.setattr(wsl, "stage_library", lambda _path: dict(staged))
    monkeypatch.setattr(wsl, "env_status", lambda: {
        "ready": False, "repeatmasker": "4.2.4", "dfam": False, "dfam_uncurated": False,
        "dfam_library": None,
    })
    monkeypatch.setattr(wsl, "_wsl_script", lambda *a, **k: (0, _annotation_output(), ""))
    result = wsl.genome_annotate("GCF_000000000.1", custom_lib="custom.fa")
    assert result["ok"] is True
    assert {(row["repeat_name"], row["class_family"], row["kind"]) for row in result["repeat_names"]} == {
        ("teA", "LTR/Gypsy", "TE"), ("teB", "LINE/L1", "TE"),
        ("low", "Simple_repeat", "tandem"), ("rib", "rRNA", "other"),
    }
    assert result["te_repeat_name_count"] == 2
    assert "families" not in result and "te_family_count" not in result


def test_unclassified_custom_matches_remain_visible_but_do_not_inflate_te_coverage(monkeypatch):
    staged = {"ok": True, "path": "$HOME/teagle_libs/abc.fa", "sha256": "c" * 64,
              "sequences": 1, "bytes": 10, "name": "custom.fa"}
    output = ("TOTAL_CHUNKS 1\nGENOME_SHA256 " + "a" * 64
              + "\nGENOME_BP 200\nTOTAL_HITS 1\nTOTAL_BP 50\nKIND\tother\t50\n"
              + "REP\tLabTE\tUnknown\t1\t50\t1.00\nANNOT_OK\n")
    monkeypatch.setattr(wsl, "stage_library", lambda _path: dict(staged))
    monkeypatch.setattr(wsl, "env_status", lambda: {
        "ready": False, "repeatmasker": "4.2.4", "dfam": False, "dfam_uncurated": False,
        "dfam_library": None,
    })
    monkeypatch.setattr(wsl, "_wsl_script", lambda *a, **k: (0, output, ""))
    result = wsl.genome_annotate("GCF_000000000.1", custom_lib="custom.fa")
    assert result["ok"] is True and result["te_bp"] == 0 and result["te_repeat_name_count"] == 0
    assert result["repeat_names"][0]["repeat_name"] == "LabTE" and result["repeat_names"][0]["kind"] == "other"
    assert "#LTR/Gypsy" in result["coverage_warning"] and "unclassified/other" in result["coverage_warning"]
    assert "classification is unavailable" in result["coverage_warning"]


def test_custom_library_manifest_records_only_the_library_actually_searched(monkeypatch):
    import copy
    import engine
    result = _engine_annotation(custom=True)
    monkeypatch.setattr(engine.wsl, "genome_annotate", lambda *a, **k: copy.deepcopy(result))
    manifest_result = engine.run_genome_annotate({
        "assemblyAccession": "GCF_000000000.1", "custom_library": "custom.fa",
        "sha256": "f" * 64})  # hostile caller value: the backend-verified hash must win
    manifest = manifest_result["provenance"]
    assert manifest["input"]["sha256"] == "a" * 64
    assert manifest["input"]["length"] == 200
    assert manifest["input"]["sha256"] != hashlib.sha256(b"").hexdigest()
    joined_refs = " ".join(ref["key"] for ref in manifest["references"])
    joined_dbs = " ".join(db["name"] for db in manifest["databases"])
    assert "Dfam" not in joined_refs and "Dfam" not in joined_dbs
    assert manifest["databases"][0] == {
        "name": "User-supplied repeat library", "sha256": "c" * 64,
        "sequences": 2, "bytes": 20}
    assert any("Dfam database search" in item for item in manifest["notRun"])


def test_custom_library_manifest_is_invariant_to_renaming_same_bytes(monkeypatch):
    import copy
    import engine
    names = iter(("first-name.fa", "renamed-library.fa"))
    def fake(*_args, **_kwargs):
        result = _engine_annotation(custom=True)
        result["custom_library"]["name"] = next(names)
        return copy.deepcopy(result)
    monkeypatch.setattr(engine.wsl, "genome_annotate", fake)
    body = {"assemblyAccession": "GCF_000000000.1", "custom_library": "chosen.fa"}
    first = engine.run_genome_annotate(body)
    second = engine.run_genome_annotate(body)
    assert first["provenance"]["manifestSha256"] == second["provenance"]["manifestSha256"]
    assert all("fileName" not in db for db in first["provenance"]["databases"])


@pytest.mark.parametrize("custom", [
    None,
    {"name": "custom.fa", "sequences": 1, "bytes": 10},
    {"name": "custom.fa", "sha256": "invalid", "sequences": 1, "bytes": 10},
    {"name": "custom.fa", "sha256": "c" * 64, "sequences": 0, "bytes": 10},
    {"name": "custom.fa", "sha256": "c" * 64, "sequences": True, "bytes": 10},
])
def test_requested_custom_library_requires_verified_identity(monkeypatch, custom):
    import copy
    import engine
    result = _engine_annotation(custom=True)
    result["custom_library"] = custom
    monkeypatch.setattr(engine.wsl, "genome_annotate", lambda *a, **k: copy.deepcopy(result))
    rejected = engine.run_genome_annotate({"assemblyAccession": "GCF_000000000.1",
                                           "custom_library": "chosen.fa"})
    assert rejected["ok"] is False and rejected["partial"] is True
    assert "custom library" in rejected["error"] and "not sealed" in rejected["error"]


def test_default_dfam_manifest_is_unchanged_in_scope(monkeypatch):
    import copy
    import engine
    result = _engine_annotation(custom=False)
    monkeypatch.setattr(engine.wsl, "genome_annotate", lambda *a, **k: copy.deepcopy(result))
    manifest = engine.run_genome_annotate({"assemblyAccession": "GCF_000000000.1"})["provenance"]
    assert any(ref["key"] == "Dfam" for ref in manifest["references"])
    assert manifest["databases"][0]["name"] == "Dfam (via famdb)"
    assert manifest["databases"][0]["version"] == "4.0"


def test_genome_manifest_requires_paired_verified_identity_and_changes_with_hash():
    empty_hash = hashlib.sha256(b"").hexdigest()
    with pytest.raises(ValueError, match="supplied together"):
        provenance.build_manifest("analysis", "ACGT", "x", {}, input_sha256="a" * 64)
    with pytest.raises(ValueError, match="supplied together"):
        provenance.build_manifest("analysis", "ACGT", "x", {}, input_length=4)
    with pytest.raises(ValueError, match="requires"):
        provenance.build_manifest("genome_annotate", "", "GCF_000000000.1", {})
    with pytest.raises(ValueError, match="positive"):
        provenance.build_manifest("genome_annotate", "", "GCF_000000000.1", {},
                                  input_sha256="a" * 64, input_length=0)
    with pytest.raises(ValueError, match="non-negative"):
        provenance.build_manifest("genome_annotate", "", "GCF_000000000.1", {},
                                  input_sha256="a" * 64, input_length=True)
    with pytest.raises(ValueError, match="64-character"):
        provenance.build_manifest("genome_annotate", "", "GCF_000000000.1", {},
                                  input_sha256="not-a-hash", input_length=200)

    a = provenance.build_manifest("genome_annotate", "", "GCF_000000000.1", {"x": 1},
                                  input_sha256="a" * 64, input_length=200)
    b = provenance.build_manifest("genome_annotate", "", "GCF_000000000.1", {"x": 1},
                                  input_sha256="b" * 64, input_length=200)
    public_digest = hashlib.sha256(b"GCF_000000000.1").hexdigest()
    assert a["input"] == {
        "id": f"record_{public_digest[:12]}",
        "idSha256": public_digest,
        "identitySource": "local_input",
        "inputType": "local_input",
        "length": 200,
        "sha256": "a" * 64,
        "sourceContextVersion": "1.0",
    }
    assert a["input"]["sha256"] != empty_hash
    assert a["manifestSha256"] != b["manifestSha256"]


def test_engine_refuses_to_seal_missing_or_invalid_backend_genome_identity(monkeypatch):
    import engine
    for sha, bp in ((None, 200), ("bad", 200), ("a" * 64, 0), ("a" * 64, True)):
        result = _engine_annotation(custom=False, genome_sha=sha)
        result["genome_bp"] = bp
        monkeypatch.setattr(engine.wsl, "genome_annotate", lambda *a, _r=result, **k: dict(_r))
        rejected = engine.run_genome_annotate({"assemblyAccession": "GCF_000000000.1",
                                               "sha256": "f" * 64})
        assert rejected["ok"] is False and rejected["partial"] is True
        assert "not sealed" in rejected["error"]
