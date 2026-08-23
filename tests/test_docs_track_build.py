"""Documentation claims that must be derived from the build, not restated by hand.

A number written into prose drifts the moment the thing it describes changes: the README claimed a
"21-model Pfam panel" after the panel had grown to 30. These tests fail when a stated figure and the
artefact it describes disagree, so the rot is caught by the suite rather than by a reader.
"""
import csv
import io
import json
import os
import re
import zipfile

import pytest

from teagle_core import domains

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _readme():
    with io.open(os.path.join(_ROOT, "README.md"), encoding="utf-8") as fh:
        return fh.read()


def _read(relpath):
    path = os.path.join(_ROOT, relpath)
    if not os.path.exists(path) and relpath.replace("\\", "/").startswith(("manuscript/", "report/")):
        pytest.skip(f"{relpath} is an ignored local publication artifact")
    with io.open(path, encoding="utf-8") as fh:
        return fh.read()


def _require_manuscript():
    if not os.path.isdir(os.path.join(_ROOT, "manuscript")):
        pytest.skip("manuscript working sources are not distributed in public source checkouts")


def _holdout_rows():
    with io.open(os.path.join(_ROOT, "benchmarks", "corpus_holdout.tsv"), encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _holdout_scores():
    with io.open(os.path.join(_ROOT, "benchmarks", "raw", "scores_holdout.json"), encoding="utf-8") as fh:
        return json.load(fh)


def _broader_scores():
    with io.open(os.path.join(_ROOT, "benchmarks", "raw", "scores.json"), encoding="utf-8") as fh:
        return json.load(fh)


def test_readme_panel_size_matches_the_bundled_profiles():
    stated = re.search(r"(\d+)-model Pfam panel", _readme())
    assert stated, "README no longer states the panel size — update this test if that is intended"
    assert int(stated.group(1)) == len(domains._hmms())


def test_readme_broader_benchmark_summary_tracks_current_scores():
    scores = _broader_scores()
    readme = _readme()
    primary = scores["overall_unique_inputs"]
    row = scores["overall"]
    element = scores["independence"]["element_only_order_withholding"]
    order = primary["order"]
    sf = primary["superfamily"]
    assert f"{order['n_gradable']} gradable content hashes at class and order rank" in readme
    assert f"answered {order['n_answered']} and all {order['n_correct']} were correct" in readme
    assert (f"accuracy {order['accuracy']:.3f}; Wilson 95% CI {order['ci95_low']:.3f}–"
            f"{order['ci95_high']:.3f}") in readme
    assert (f"Among {element['n_gradable']} gradable element hashes, it answered "
            f"{element['n_answered']} and withheld the order on {element['n_withheld']}") in readme
    assert (f"row-level sensitivity summary is {row['order']['n_correct']}/"
            f"{row['order']['n_answered']} answered among {row['order']['n_gradable']} gradable") in readme
    assert (f"{sf['n_correct']} of {sf['n_answered']} answered calls were correct among "
            f"{sf['n_gradable']} gradable hashes") in readme
    assert (f"accuracy {sf['accuracy']:.3f}; Wilson 95% CI {sf['ci95_low']:.3f}–"
            f"{sf['ci95_high']:.3f}") in readme


def test_active_comparison_scripts_make_no_preregistration_timing_claim():
    for path in ("benchmarks/head_to_head.py", "reviewer_package/scripts/head_to_head.py"):
        text = _read(path).lower()
        for forbidden in ("pre-registered", "preregistered", "pre-register", "preregister"):
            assert forbidden not in text, f"{path} contains timing claim {forbidden!r}"


# ---------------- report/manual capability enumeration (gitignored working files) ----------------
def _doc(name):
    p = os.path.join(_ROOT, "report", name)
    if not os.path.exists(p):
        return None
    with io.open(p, encoding="utf-8") as fh:
        return fh.read()


@pytest.mark.parametrize("doc", ["teagle_manual.tex", "teagle_report.tex"])
def test_doc_states_the_correct_panel_size(doc):
    """The report and manual both enumerate the panel; the stated model count must equal the bundled set,
    so growing the panel without updating the docs fails the build (README once claimed 21 after it hit 30)."""
    text = _doc(doc)
    if text is None:
        pytest.skip(f"{doc} not present (report/ is a gitignored working dir)")
    n = len(domains.DOMAIN_INFO)
    assert re.search(rf"\b{n}[ ~]?(?:public-domain )?Pfam", text) or re.search(rf"\b{n} Pfam profile", text), \
        f"{doc} does not state the current panel size of {n} Pfam models"


@pytest.mark.parametrize("doc", ["teagle_manual.tex", "teagle_report.tex"])
def test_doc_lists_every_superfamily_the_classifier_can_name(doc):
    """Every superfamily/lineage name classify() can emit must appear in the enumeration, or the docs
    under-state what the tool does. Helitron is deliberately NOT here: its helicase is a detected DOMAIN
    (checked by the accession test) but classify() emits no Helitron superfamily call, so listing it as a
    superfamily would over-state the tool."""
    text = _doc(doc)
    if text is None:
        pytest.skip(f"{doc} not present")
    for name in ("Copia", "Gypsy", "LINE", "DIRS", "Tc1/Mariner", "hAT", "CACTA", "MULE", "IS4"):
        assert name in text, f"{doc} does not name the {name!r} superfamily/group in its enumeration"


@pytest.mark.parametrize("doc", ["teagle_manual.tex", "teagle_report.tex"])
def test_doc_does_not_claim_a_helitron_superfamily_call(doc):
    """classify() detects the Helitron helicase domain but emits no Helitron superfamily call. The docs must
    say so, not imply the tool assigns a Helitron superfamily (guards against re-introducing that overclaim)."""
    import sys as _sys
    _be = os.path.join(_ROOT, "app", "backend")
    if _be not in _sys.path:
        _sys.path.insert(0, _be)
    from teagle_core import classify
    cl = classify.classify([], [{"domain": "HEL", "class": "dna:Helitron",
                                 "nt": [100, 900], "strand": "+", "score": 150.0}])
    assert not cl["te_class"].startswith("HEL"), "classify now emits a Helitron call — update the docs and this test"
    text = _doc(doc)
    if text is None:
        pytest.skip(f"{doc} not present")
    assert "does not currently emit" in text or "no dedicated" in text, \
        f"{doc} must state that no dedicated Helitron superfamily call is emitted"


@pytest.mark.parametrize("doc", ["teagle_manual.tex", "teagle_report.tex"])
def test_doc_has_a_decision_methodology_table(doc):
    """Each capability doc must carry the 'how each call is made' table (label tab:calls), so a reader sees
    the decision rule behind each class/superfamily/family/domain call, not only the vocabulary."""
    text = _doc(doc)
    if text is None:
        pytest.skip(f"{doc} not present")
    assert r"\label{tab:calls}" in text, f"{doc} is missing the decision-methodology table (tab:calls)"
    for rule in ("N-terminal", "C-terminal", "tyrosine recombinase", "80--80--80", "strongest-scoring"):
        assert rule in text, f"{doc} decision table omits the {rule!r} rule"


def test_report_oligoqc_table_row_count_matches_the_benchmark():
    """The report's primer-QC table (tab:oligoqc) must have one data row per benchmarked pair. It had
    drifted to a stale 6-row set (different primers, one non-reproducing value) while the benchmark held 12;
    this ties the table size to LIT_PRIMERS so the primer set cannot silently diverge again. (Row COUNT
    only — the per-cell ΔG values are ViennaRNA-dependent and so not hermetically assertable.)"""
    text = _doc("teagle_report.tex")
    if text is None:
        pytest.skip("report not present")
    import sys as _sys, re as _re
    _t = os.path.join(_ROOT, "tests")
    if _t not in _sys.path:
        _sys.path.insert(0, _t)
    from test_benchmarks import LIT_PRIMERS
    body = text.split(r"\label{tab:oligoqc}", 1)[1].split(r"\end{tabularx}", 1)[0]
    rows = [l for l in body.splitlines() if l.strip().endswith(r"\\") and ("ok" in l or "caution" in l or "warn" in l)]
    assert len(rows) == len(LIT_PRIMERS), \
        f"report oligoqc table has {len(rows)} data rows but the benchmark has {len(LIT_PRIMERS)} pairs"


@pytest.mark.parametrize("doc", ["teagle_manual.tex", "teagle_report.tex"])
def test_doc_lists_every_reported_domain_code_pfam_accession(doc):
    """Each Pfam accession in the bundled panel must appear in the doc's domain table — the table is the
    ground-truth enumeration the user asked for, so a model absent from it is a documentation defect."""
    text = _doc(doc)
    if text is None:
        pytest.skip(f"{doc} not present")
    for hmm, (code, label, cls, pfam) in domains.DOMAIN_INFO.items():
        assert pfam in text, f"{doc} omits {pfam} ({hmm}) from its domain-panel table"


def test_readme_version_badge_matches_the_source_version():
    from teagle_core import __version__
    stated = re.search(r"badge/version-([0-9.]+)-", _readme())
    assert stated, "README version badge missing"
    assert stated.group(1) == __version__


@pytest.mark.parametrize("relpath", ("report/teagle_manual.tex", "report/teagle_report.tex"))
def test_release_document_version_matches_the_source_version(relpath):
    """Ignored release documents must not retain an older version than the shipped build."""
    from teagle_core import __version__
    path = os.path.join(_ROOT, relpath)
    if not os.path.exists(path):
        pytest.skip(f"{relpath} not present in this checkout")
    stated = re.search(r"\\newcommand\{\\teagleversion\}\{([0-9.]+)\}", _read(relpath))
    assert stated, f"{relpath} does not define \\teagleversion"
    assert stated.group(1) == __version__


@pytest.mark.parametrize("doc", ["teagle_manual.tex", "teagle_report.tex"])
def test_genome_scan_docs_do_not_claim_a_fixed_runtime(doc):
    text = _doc(doc)
    if text is None:
        pytest.skip(f"{doc} not present")
    flat = " ".join(text.split())
    assert not re.search(r"(?:one|1)\s*(?:to|[-–])\s*2\s*(?:minutes?|mins?)", text, re.IGNORECASE)
    assert "typically scanned" not in text.lower()
    assert ("Runtime depends on genome size, the installed isPcr build and cached RefSeq assembly, local "
            "hardware, and the fixed candidate-product-size cap.") in flat
    assert "in-progress state while the scan runs and reports a timeout as a failure" in flat


def test_native_genome_scan_wording_matches_runtime_behavior():
    native = _read("app/native/main.py")
    wsl = _read("app/backend/teagle_core/wsl.py")
    flat = " ".join(native.split())
    assert not re.search(r"(?:one|1)\s*(?:to|[-–~])\s*2\s*(?:minutes?|mins?)", native, re.IGNORECASE)
    assert "fast once the genome is cached" not in native.lower()
    assert "future scans are fast" not in native.lower()
    assert ("Runtime depends on genome size, the installed isPcr build and cached RefSeq assembly, local "
            "hardware, and the fixed candidate-product-size cap") in flat
    assert "_render_genome_busy" in native
    assert "max_size: int = 4000" in wsl
    assert "timeout: int = 600" in wsl
    assert '"error": f"isPcr timed out after {timeout}s"' in wsl


def test_current_release_text_tracks_version_and_schema_contracts():
    """Current shipped text must distinguish the 3.9 contracts from historical 3.8 evidence."""
    from teagle_core import __version__, provenance
    import engine

    citation = _read("CITATION.cff")
    stated = re.search(r"(?m)^version:\s*([0-9.]+)\s*$", citation)
    assert stated and stated.group(1) == __version__
    released = re.search(r"(?m)^date-released:\s*(\d{4}-\d{2}-\d{2})\s*$", citation)
    assert released, "a published release must carry its release date"

    readme = _read("README.md")
    app_readme = _read("app/README.md")
    changelog = _read("CHANGELOG.md")
    try_it = _read("docs/TRY-IT-YOURSELF.md")
    combined = "\n".join((readme, app_readme, changelog, try_it))
    assert f"## [{__version__}] — {released.group(1)}" in changelog
    assert f"Manifest schema {provenance.SCHEMA_VERSION}" in combined
    assert (f"result schema {engine.ANALYSIS_RESULT_SCHEMA}" in combined or
            f"analysis-result schema {engine.ANALYSIS_RESULT_SCHEMA}" in combined)
    assert f"source-context schema {provenance.SOURCE_CONTEXT_VERSION}" in combined
    assert "unkeyed integrity" in combined
    assert "clean-machine" in combined
    assert "No new 3.9.0 classification, assay, or prospective benchmark result" in combined
    for stale in (
        "under release validation", "has not yet been published",
        "No 3.9.0 installer or manual asset is public yet",
        "still awaiting rendering and page-by-page release QA",
    ):
        assert stale not in combined


def test_current_docs_state_native_and_legacy_export_boundaries():
    readme = _read("README.md")
    app_readme = _read("app/README.md")
    for required in (
        "XLSX/CSV/TSV", "SVG/PNG/vector PDF", "GFF3/BED/FASTA",
        "provenance JSON", "complete-result JSON",
    ):
        assert required in readme
    assert "legacy browser" in readme.lower()
    assert "disables XLSX, PDF, provenance JSON, and complete-result JSON" in app_readme


def test_screenshot_index_records_completed_39_capture_contract():
    text = _read("docs/img/README.md")
    for required in (
        "native TEagle 3.9.0 app", "light theme at 1920x1080",
        "1500x1150 window and UI scale 1.0", "Every accepted figure is a real run",
        "No figure contains synthesised results", "recapture sidecar",
        "dimensions and SHA-256 hashes",
    ):
        assert required in text
    assert "must be recaptured" not in text


def test_manual_cites_subregion_picker_with_coordinate_contract():
    manual_path = os.path.join(_ROOT, "report", "teagle_manual.tex")
    if not os.path.exists(manual_path):
        pytest.skip("report/teagle_manual.tex not present in this checkout")
    manual = _read("report/teagle_manual.tex")
    assert "{subregion_picker.png}" in manual
    assert "1-based, inclusive interval" in manual
    assert "0-based, half-open span" in manual
    assert "primer design or splice detection" in manual


def test_manual_keeps_ltr_preposition_with_its_subject():
    manual = _read("report/teagle_manual.tex")

    def assert_phrase_kept(text):
        assert "For~LTR~retrotransposons" in text
        assert "For~LTR retrotransposons" not in text
        assert "For LTR retrotransposons" not in text

    assert_phrase_kept(manual)
    with pytest.raises(AssertionError):
        assert_phrase_kept(manual.replace("For~LTR~retrotransposons", "For~LTR retrotransposons", 1))


def test_manual_keeps_qpcr_gc_clamp_together():
    manual = _read("report/teagle_manual.tex")

    def assert_phrase_kept(text):
        assert "$3'$ GC~clamp." in text
        assert "$3'$ GC clamp." not in text

    assert_phrase_kept(manual)
    with pytest.raises(AssertionError):
        assert_phrase_kept(manual.replace("$3'$ GC~clamp.", "$3'$ GC clamp.", 1))


def test_manual_keeps_sensitivity_item_clear_of_the_footer():
    manual = _read("report/teagle_manual.tex")
    guard = "\\needspace{6\\baselineskip}\n\\item[Sensitivity]"

    def assert_guard(text):
        assert guard in text

    assert_guard(manual)
    with pytest.raises(AssertionError):
        assert_guard(manual.replace(guard, r"\item[Sensitivity]", 1))


def test_report_keeps_terminal_run_length_with_its_criterion():
    report = _read("report/teagle_report.tex")
    kept = r"length~$r \ge r_{\min}$"
    breakable = r"length $r \ge r_{\min}$"

    def assert_phrase_kept(text):
        assert kept in text
        assert breakable not in text

    assert_phrase_kept(report)
    with pytest.raises(AssertionError):
        assert_phrase_kept(report.replace(kept, breakable, 1))


def test_report_keeps_no_longer_together():
    report = _read("report/teagle_report.tex")
    kept = "no~longer"
    breakable = "no longer"

    def assert_phrase_kept(text):
        assert kept in text
        assert breakable not in text

    assert_phrase_kept(report)
    with pytest.raises(AssertionError):
        assert_phrase_kept(report.replace(kept, breakable, 1))


def test_report_keeps_stacked_rows_together():
    report = _read("report/teagle_report.tex")
    kept = "stacked~rows"
    breakable = "stacked rows"

    def assert_phrase_kept(text):
        assert kept in text
        assert breakable not in text

    assert_phrase_kept(report)
    with pytest.raises(AssertionError):
        assert_phrase_kept(report.replace(kept, breakable, 1))


def test_report_keeps_hue_is_together():
    report = _read("report/teagle_report.tex")
    kept = "Hue~is"
    breakable = "Hue is"

    def assert_phrase_kept(text):
        assert kept in text
        assert breakable not in text

    assert_phrase_kept(report)
    with pytest.raises(AssertionError):
        assert_phrase_kept(report.replace(kept, breakable, 1))


def test_manual_contents_card_navigation_has_no_gap_or_duplicate():
    manual = _read("report/teagle_manual.tex")
    navigation = re.compile(
        r"\\section\{Loading a sequence\}"
        r"|\\section\{Card (?P<section>\d{2}) --- (?P<section_title>[^{}]+)\}"
    )
    expected = [
        (1, "Loading a sequence"),
        (2, "Classification and structure"),
        (3, "Searching for a family match (optional)"),
        (4, "Designing PCR primers"),
        (5, "Predicting amplification (in-silico PCR)"),
        (6, "Whole-genome off-target scan"),
        (7, "The run provenance"),
    ]
    target = "\n".join((
        r"\section{Card 06 --- Whole-genome off-target scan}",
        r"\label{sec:genome-scan}",
    ))

    def assert_navigation(text):
        entries = []
        for match in navigation.finditer(text):
            if match.group(0) == r"\section{Loading a sequence}":
                entries.append((1, "Loading a sequence"))
            else:
                number = int(match.group("section"))
                title = match.group("section_title")
                entries.append((number, title))
        assert entries == expected
        assert target in text
        assert r"\addcontentsline{toc}" not in text

    card_06 = r"\section{Card 06 --- Whole-genome off-target scan}"
    assert_navigation(manual)
    with pytest.raises(AssertionError):
        assert_navigation(manual.replace(card_06 + "\n", "", 1))
    with pytest.raises(AssertionError):
        assert_navigation(manual.replace(card_06, card_06 + "\n" + card_06, 1))
    with pytest.raises(AssertionError):
        assert_navigation(manual + "\n" + r"\addcontentsline{toc}{section}{Card 06}")


def test_version_locked_benchmark_docs_match_their_derived_evidence_version():
    """Benchmark evidence keeps the version that generated it; an application patch must not relabel it."""
    with io.open(os.path.join(_ROOT, "reviewer_package", "results", "provenance_summary.json"),
                 encoding="utf-8") as handle:
        summary = json.load(handle)
    versions = summary["version_coverage"]["teagle_versions"]
    assert len(versions) == 1
    evidence_version, evidence_cases = next(iter(versions.items()))
    assert evidence_cases == summary["case_count"]
    assert f"TEagle {evidence_version}." in _read("reviewer_package/README.md")
    spec = _read("benchmarks/protocol/BENCHMARK_SPEC_v1.md")
    assert f"current execution pin is **{evidence_version}**" in spec
    # manuscript/ is ignored in public checkouts; if locally maintained, its evidence caption is still pinned.
    figure_index = os.path.join(_ROOT, "manuscript", "figures", "figure-index.md")
    if os.path.exists(figure_index):
        assert f"TEagle {evidence_version}." in _read("manuscript/figures/figure-index.md")


def test_documented_holdout_summary_is_derived_from_corpus_and_scores():
    """Validation-corpus prose must follow the current rows and their scored denominator."""
    rows = _holdout_rows()
    scores = _holdout_scores()
    n_records = len(rows)
    n_organism_values = len({row["organism"] for row in rows})
    n_answered = scores["overall"]["order"]["n_answered"]
    n_correct = scores["overall"]["order"]["n_correct"]
    n_negative = sum(
        row["expected_class"].strip().lower().startswith("not a te") for row in rows
    )
    n_te = sum(
        not row["expected_class"].strip().lower().startswith("not a te") for row in rows
    )
    assert scores["negative_controls"]["n"] == n_negative
    assert scores["overall"]["order"]["n_gradable"] == n_records

    summary = re.search(
        r"Measured on a development-example-excluded retrospective literature-derived corpus: "
        r"(?P<records>\d+) records\s+carrying (?P<organism_values>\d+) distinct literal organism annotations.*?"
        r"named the Wicker order for (?P<answered>\d+) of the (?P<te_records>\d+) transposable elements\s+"
        r"and was\s+correct on all (?P<correct>\d+), and declined to name an order on all "
        r"(?P<negative_controls>\d+) negative controls\.",
        _readme(),
        re.S,
    )
    assert summary, "README validation summary is missing or no longer machine-checkable"
    assert int(summary["records"]) == n_records
    assert int(summary["organism_values"]) == n_organism_values
    assert int(summary["answered"]) == n_answered
    assert int(summary["te_records"]) == n_te
    assert int(summary["correct"]) == n_correct
    assert int(summary["negative_controls"]) == n_negative

    reviewer = _read("reviewer_package/README.md")
    assert f"{n_records} records carrying {n_organism_values} distinct literal organism annotations" in reviewer


def test_manuscript_uses_literal_holdout_organism_annotations_and_cautious_specification_language():
    """The manuscript must not silently turn literal corpus annotations into a species count or timing claim."""
    _require_manuscript()
    import sys as _sys
    manuscript_dir = os.path.join(_ROOT, "manuscript")
    if manuscript_dir not in _sys.path:
        _sys.path.insert(0, manuscript_dir)
    import values

    rows = _holdout_rows()
    rendered_values = values.build()
    assert "holdout_species" not in rendered_values
    assert rendered_values["holdout_organism_annotations"] == len(
        {row["organism"].strip() for row in rows}
    )

    text = _read("manuscript/manuscript.md")
    assert re.search(
        r"{{holdout_organism_annotations}} distinct\s+literal organism annotations", text
    )
    for stale in (
        "{{holdout_species}}",
        "independent corpus",
        "assembled independently of the software",
        "pre-registers",
        "pre-registered",
        "pre-registering",
        "before the search",
        "before these panels",
        "available from the corresponding author",
        "available on request",
    ):
        assert stale not in text
    flat = " ".join(text.split())
    assert "no prospective preregistration is claimed" in flat
    assert flat.lower().count("preregistration") == 1
    assert "dated, versioned specification" not in flat
    assert "amendments made during execution" not in flat


def test_manuscript_primary_units_duplicate_and_comparator_accounting_are_data_bound():
    """Primary units and comparison accounting must be derived from current case-level evidence."""
    _require_manuscript()
    import sys as _sys
    manuscript_dir = os.path.join(_ROOT, "manuscript")
    if manuscript_dir not in _sys.path:
        _sys.path.insert(0, manuscript_dir)
    import values

    rendered_values = values.build()
    scores = _broader_scores()
    source = _read("manuscript/manuscript.md")
    flat = " ".join(source.split())
    abstract = flat.split("## Abstract", 1)[1].split("## Abbreviations", 1)[0]

    primary = scores["overall_unique_inputs"]["order"]
    assert rendered_values["broad_unique_order_gradable"] == primary["n_gradable"]
    assert rendered_values["broad_unique_order_answered"] == primary["n_answered"]
    assert rendered_values["broad_unique_order_correct"] == primary["n_correct"]
    assert rendered_values["broad_unique_order_withheld"] == (
        primary["n_gradable"] - primary["n_answered"]
    )
    assert "{{broad_unique_order_gradable}} distinct gradable input hashes" in abstract
    assert ("among {{broad_element_unique_gradable}} element hashes it answered "
            "{{broad_element_unique_answered}}") in abstract
    assert "withheld the order on {{broad_element_unique_withheld}}" in abstract
    assert "any of the {{broad_neg_n}} broader control hashes" in abstract
    assert "{{broad_order_answered}}" not in abstract
    assert "{{broad_gradable}}" not in abstract

    gradable = [case for case in scores["cases"]
                if not (case["expected_class"] != "NOT_TE" and case["expected_order"] == "OTHER")]
    rows_by_hash = {}
    for case in gradable:
        rows_by_hash.setdefault(case["input_sha256"], []).append(case)
    duplicates = [rows for rows in rows_by_hash.values() if len(rows) > 1]
    assert len(duplicates) == 1
    assert rendered_values["broad_duplicate_rows"] == len(duplicates[0])
    assert rendered_values["broad_duplicate_accession"] == duplicates[0][0]["accession"]
    assert "{{broad_duplicate_accession}}" in source
    assert "sole duplicated gradable input" in flat

    sf_cases = [case for case in scores["cases"] if case.get("superfamily_gradable")]
    sf_order_answered = [case for case in sf_cases
                         if not case.get("abstained_at_order")
                         and case.get("observed_order") != "NO_CALL"]
    assert rendered_values["broad_sf_same_order_gradable"] == len(sf_cases)
    assert rendered_values["broad_sf_same_order_withheld"] == len(sf_cases) - len(sf_order_answered)
    assert rendered_values["broad_sf_withheld"] == len(sf_cases) - scores["overall"]["superfamily"]["n_answered"]
    assert ("Both rates use the same case set" in flat
            and "{{broad_sf_same_order_gradable}} cases" in flat)
    assert ("superfamily on {{broad_sf_withheld}} of {{broad_sf_same_order_gradable}}"
            in flat)

    with io.open(os.path.join(_ROOT, "benchmarks", "raw", "head_to_head.json"),
                 encoding="utf-8") as fh:
        comparison = json.load(fh)
    historical = comparison["historical_vs_current"]
    assert rendered_values["cmp_historical_tesorter_rows"] == historical["historical_tesorter_corpus_rows"]
    assert rendered_values["cmp_current_teagle_rows"] == historical["current_teagle_corpus_rows"]
    for mode in ("default", "matched"):
        evidence = comparison["modes"][mode]
        accounting = evidence["accounting"]
        assert rendered_values[f"cmp_{mode}_n"] == evidence["n_paired"]
        assert rendered_values[f"cmp_{mode}_both_correct"] == evidence["both_correct"]
        assert rendered_values[f"cmp_{mode}_teagle_only"] == evidence["teagle_only"]
        assert rendered_values[f"cmp_{mode}_tesorter_only"] == evidence["tesorter_only"]
        assert rendered_values[f"cmp_{mode}_neither"] == evidence["neither"]
        assert rendered_values[f"cmp_{mode}_successful_accessions"] == accounting["current_successful_unique_accessions"]
        assert rendered_values[f"cmp_{mode}_comparator_outputs"] == accounting["comparator_outputs"]
        assert rendered_values[f"cmp_{mode}_excluded_with_output"] == accounting["excluded_with_output"]
        assert rendered_values[f"cmp_{mode}_missing_comparator_output"] == accounting["missing_comparator_output"]["n"]
        assert rendered_values[f"cmp_{mode}_extra_comparator_outputs"] == accounting["extra_comparator_outputs"]["n"]
    assert "All comparative inferences below are conditional on the {{cmp_default_n}} paired cases" in flat


def test_manuscript_divergence_boundary_observations_are_case_bound():
    """Observed boundary values must come from the replicate ledger, not typed prose."""
    _require_manuscript()
    import sys as _sys
    manuscript_dir = os.path.join(_ROOT, "manuscript")
    if manuscript_dir not in _sys.path:
        _sys.path.insert(0, manuscript_dir)
    import values

    with io.open(os.path.join(_ROOT, "benchmarks", "raw", "sim_divergence.json"),
                 encoding="utf-8") as fh:
        raw = json.load(fh)
    by_target = {}
    for record in raw["records"]:
        by_target.setdefault(record["target_identity"], []).append(record)
    fully = [target for target, rows in by_target.items()
             if all(row["detected_as_ltr"] for row in rows)]
    partial = [target for target, rows in by_target.items()
               if 0 < sum(row["detected_as_ltr"] for row in rows) < len(rows)]
    advisory = [target for target, rows in by_target.items()
                if any(row["advisory_reported"] for row in rows)]
    never = [target for target, rows in by_target.items()
             if not any(row["detected_as_ltr"] for row in rows)]
    rendered_values = values.build()
    assert rendered_values["div_highest_fully_accepted"] == f"{max(fully):g}%"
    assert rendered_values["div_transition_low"] == f"{min(partial):g}"
    assert rendered_values["div_transition_high"] == f"{max(partial):g}"
    assert rendered_values["div_highest_advisory"] == f"{max(advisory):g}"
    highest_never_rows = by_target[max(never)]
    assert rendered_values["div_highest_never_accepted_n"] == sum(
        row["detected_as_ltr"] for row in highest_never_rows
    )

    source = _read("manuscript/manuscript.md")
    flat = " ".join(source.split())
    for key in (
        "div_highest_fully_accepted", "div_highest_never_accepted_n", "div_transition_low",
        "div_transition_high", "div_highest_advisory",
    ):
        assert f"{{{{{key}}}}}" in source
    for stale in (
        "from 100% down to {{div_lowest_accepted}}", "absent (0/{{div_replicates}})",
        "transition at 79–79.5%", "candidate pairs across 79.5% to",
    ):
        assert stale not in flat

    # These are protocol/scaffold constants, not observed performance summaries.
    assert raw["scaffold_native_identity"] == 100
    assert max(raw["targets"]) == 100 and min(raw["targets"]) == 60
    assert "terminal repeat pair at 100% identity" in flat
    assert "identity targets from 100% to 60%" in flat


def test_preprint_layout_source_requires_landscape_native_tables_and_caption_attachment():
    """Source-level guard for the DOCX layout policy; the artifact is checked below when present locally."""
    source = _read("manuscript/build_docx.py")
    for required in (
        "WD_ORIENT.LANDSCAPE",
        "WD_ORIENT.PORTRAIT",
        "w:tblHeader",
        "w:cantSplit",
        "keep_with_next",
        "keep_together",
        "WD_ALIGN_PARAGRAPH.LEFT",
        "WD_ALIGN_PARAGRAPH.JUSTIFY",
        "restart\"), \"continuous\"",
        "w:instrText",
        'instruction.text = " PAGE "',
        "_page_number_footer(section)",
    ):
        assert required in source


def test_manuscript_scientific_typography_source_is_explicit():
    """Scientific notation and biological names must carry unambiguous source markup."""
    text = _read("manuscript/manuscript.md")
    flat = " ".join(text.split())
    assert "10<sup>−6</sup>" in text
    assert "10<sup>4</sup>" in text
    assert not re.search(r"10(?:\^|\*\*)[-−]?\d", text)
    assert "3′" in text and "5′" in text
    assert not re.search(r"[35]['’]", text)
    assert text.count("*P* =") == 2
    assert "*p* =" not in text
    assert text.count("*n* =") == 5
    assert "Pol polyprotein" in text and "pol polyprotein" not in text
    assert "Gag, protease or envelope" in text and "gag, protease or envelope" not in text
    assert text.count("*IS711*") == 2
    assert "IS*711*" not in text
    for styled in (
        "*Euplotes aediculatus*", "*Gadus morhua*", "*Aspergillus nidulans*",
        "*Fusarium oxysporum*", "*Brucella melitensis*", "*Brucella ovis*",
        "*Homo sapiens*", "*Cichorium intybus*", "*niaD*", "*ACE*", "*pol*",
    ):
        assert styled in flat


def test_manuscript_claim_bearing_prose_uses_direct_scoped_wording():
    """Guard the reviewed claim passages without banning necessary scientific contrasts globally."""
    text = _read("manuscript/manuscript.md")
    flat = " ".join(text.split())
    for stale in (
        "Beneath this practical gap lies an epistemic one",
        "What is scarce is not the idea but its reach",
        "Two consequences follow",
        "More informative is that",
        "One property of this comparison has to be stated plainly",
        "The result that governs how the rest of this comparison must be read",
        "One qualification belongs with that statement",
        "we say so here rather than leave it to be inferred",
        "The decomposition isolates the observed cost",
        "The second observation concerns scope",
        "What the survey does support is narrower",
        "The general point survives whatever the count",
        "Both push in the same direction, making the recovery reported here optimistic",
        "The ground truth for the assay layer is external in a way a classification label is not",
        "[8, software documentation]",
        "DANTE_LTR reports each element",
        "The gap between the two is the honest measure",
        "informative rather than random",
        "A provenance claim is worth only as much as its weakest property",
        "Stability: does the same input yield the same seal",
        "Sensitivity: do different inputs yield different seals",
        "Coverage: does the seal record every threshold",
        "Classification is compared with literature-derived labels",
        "The idea of withholding rather than guessing is not new",
        "Withholding a call where the evidence cannot support one is not a new idea",
    ):
        assert stale not in flat
    for required in (
        "DANTE reports the lowest common ancestor of ambiguous top hits rather than forcing a lineage call [8]",
        "Evidence-dependent withholding predates TEagle",
        "We assessed classification against labels curated from the publications describing each element",
        "we infer that this rule can censor older elements; this panel does not estimate how often that occurs",
        "The assay targets and reported product sizes were externally sourced; scoring nevertheless required documented choices about template, locus and tolerance",
        "The panel was not scored blind",
        "The paired set is restricted to cases carrying a gradable order label",
        "The zero-error result does not extend to the whole corpus",
        "The survey supports a narrower conclusion",
        "The comparison shows why abstention requires external evaluation",
        "The primary superfamily denominator excludes {{broad_sf_erv_excluded}} endogenous-retrovirus rows",
        "The difference between the two estimates quantifies the result's sensitivity to that decision",
        "TEagle made two superfamily errors on real elements",
        "We assessed seal stability, sensitivity and parameter coverage separately",
        "Parameter coverage was complete in this audit",
        "TEagle analysed {{prov_inputs}} synthetic inputs {{prov_replicates}} times each",
        "We assessed classification against literature-derived labels",
        "We measured terminal-repeat detection on a simulated divergence series",
        "Earlier tools established evidence-dependent withholding",
    ):
        assert required in flat


def test_optional_local_manuscript_reports_all_three_reliability_axes():
    """The manuscript must not omit the HMMER call-confidence axis used by the shipped build."""
    _require_manuscript()
    flat = " ".join(_read("manuscript/manuscript.md").split())
    assert "Classification and its three reliability axes" in flat
    assert "per-domain call confidence" in flat
    assert "domain completeness" in flat
    assert "coding-structure" in flat
    assert "Classification and its two reliability axes" not in flat
    assert "Reliability is reported on two axes" not in flat


def test_docx_inline_renderer_uses_true_superscript_and_subscript_runs():
    """HTML-like scientific markup must become Word run properties, not smaller baseline text."""
    _require_manuscript()
    import sys as _sys
    from docx import Document

    manuscript_dir = os.path.join(_ROOT, "manuscript")
    if manuscript_dir not in _sys.path:
        _sys.path.insert(0, manuscript_dir)
    import build_docx

    paragraph = Document().add_paragraph()
    build_docx.add_inline(paragraph, "10<sup>−6</sup> and H<sub>2</sub>O")
    exponent = next(run for run in paragraph.runs if run.text == "−6")
    subscript = next(run for run in paragraph.runs if run.text == "2")
    assert exponent.font.superscript is True
    assert exponent.font.subscript is False
    assert subscript.font.subscript is True
    assert subscript.font.superscript is False


def test_optional_local_preprint_artifact_has_readable_table_and_figure_structure():
    """Optional local-artifact check: generated manuscript files are gitignored in some CI checkouts."""
    path = os.path.join(_ROOT, "manuscript", "manuscript_preprint.docx")
    if not os.path.exists(path):
        pytest.skip("preprint artifact not present (manuscript/ is a gitignored working directory)")
    from docx import Document
    from docx.oxml.ns import qn

    document = Document(path)
    landscape_sections = [section for section in document.sections
                          if section.page_width > section.page_height]
    portrait_sections = [section for section in document.sections
                         if section.page_width < section.page_height]
    assert len(document.tables) == 4
    assert len(landscape_sections) == len(document.tables)
    base_portrait = portrait_sections[0]
    base_margins = (
        base_portrait.left_margin, base_portrait.right_margin,
        base_portrait.top_margin, base_portrait.bottom_margin,
    )
    assert all((section.left_margin, section.right_margin,
                section.top_margin, section.bottom_margin) == base_margins
               for section in portrait_sections)
    for section in document.sections:
        numbering = section._sectPr.find(qn("w:lnNumType"))
        assert numbering is not None
        assert numbering.get(qn("w:restart")) == "continuous"
        footer_xml = section.footer._element.xml
        assert re.search(r"<w:instrText[^>]*>\s*PAGE\s*</w:instrText>", footer_xml)
        field_types = set(re.findall(r'w:fldCharType="([^"]+)"', footer_xml))
        assert {"begin", "separate", "end"}.issubset(field_types)
    for table in document.tables:
        assert table.rows[0]._tr.get_or_add_trPr().find(qn("w:tblHeader")) is not None
        for row in table.rows:
            assert row._tr.get_or_add_trPr().find(qn("w:cantSplit")) is not None
    for table_number in (1, 2, 4):
        table = document.tables[table_number - 1]
        frame_width = (landscape_sections[table_number - 1].page_width
                       - landscape_sections[table_number - 1].left_margin
                       - landscape_sections[table_number - 1].right_margin)
        assert sum(cell.width for cell in table.rows[0].cells) < frame_width
        assert table.rows[0].cells[0].width.mm >= 28
    figure_paragraphs = [paragraph for paragraph in document.paragraphs
                         if "<w:drawing" in paragraph._p.xml]
    assert len(figure_paragraphs) == 3
    assert all(paragraph.paragraph_format.keep_with_next for paragraph in figure_paragraphs)
    captions = [paragraph for paragraph in document.paragraphs
                if re.match(r"^Figure \d+\.", paragraph.text)]
    assert len(captions) == 3
    assert all(paragraph.paragraph_format.keep_together for paragraph in captions)
    table_captions = [paragraph for paragraph in document.paragraphs
                      if re.match(r"^Table \d+\.", paragraph.text)]
    assert len(table_captions) == 4
    assert all(paragraph.paragraph_format.keep_together for paragraph in table_captions)
    assert all(paragraph.paragraph_format.keep_with_next for paragraph in table_captions)


def test_optional_local_manuscript_table_captions_do_not_split_across_pages():
    """Every manuscript rendition keeps each complete table caption in one paragraph box."""
    from docx import Document

    for filename in ("manuscript.docx", "manuscript_preprint.docx"):
        path = os.path.join(_ROOT, "manuscript", filename)
        if not os.path.exists(path):
            pytest.skip(f"{filename} not present")
        captions = [paragraph for paragraph in Document(path).paragraphs
                    if re.match(r"^Table \d+\.", paragraph.text)]
        assert len(captions) == 4
        assert all(paragraph.paragraph_format.keep_together for paragraph in captions)


def test_optional_local_manuscript_artifacts_preserve_scientific_run_properties():
    """Inspect the generated OOXML-level run styling for notation, taxa, genes and proteins."""
    from docx import Document

    for filename in ("manuscript.docx", "manuscript_preprint.docx"):
        path = os.path.join(_ROOT, "manuscript", filename)
        if not os.path.exists(path):
            pytest.skip(f"{filename} not present")
        document = Document(path)
        paragraphs = list(document.paragraphs)
        cells = [cell for table in document.tables for row in table.rows for cell in row.cells]
        paragraph_runs = [run for paragraph in paragraphs for run in paragraph.runs]
        runs = list(paragraph_runs)
        runs.extend(run for cell in cells for paragraph in cell.paragraphs for run in paragraph.runs)
        full_text = " ".join([paragraph.text for paragraph in paragraphs] + [cell.text for cell in cells])

        assert "<sup>" not in full_text and "<sub>" not in full_text
        assert "3′" in full_text and "5′" in full_text
        assert not re.search(r"[35]['’]", full_text)
        for exponent in ("−6", "4"):
            matching = [run for run in paragraph_runs if run.text == exponent]
            assert matching and all(run.font.superscript is True for run in matching)
        assert len([run for run in runs if run.text == "P" and run.font.italic is True]) == 2
        assert any(run.text == "E" and run.font.italic is True for run in runs)
        assert len([run for run in runs if run.text == "n" and run.font.italic is True]) == 5

        for gene_or_taxon in (
            "IS711", "niaD", "ACE", "pol", "Euplotes aediculatus", "Gadus morhua",
            "Aspergillus nidulans", "Fusarium oxysporum", "Brucella melitensis",
            "Brucella ovis", "Homo sapiens", "Cichorium intybus",
        ):
            matching = [run for run in runs if run.text == gene_or_taxon]
            assert matching and all(run.font.italic is True for run in matching), (filename, gene_or_taxon)
        assert any("Pol polyprotein" in run.text and run.font.italic is False for run in runs)
        assert any("Gag, protease or envelope" in run.text and run.font.italic is False for run in runs)

        if document.tables:
            table_runs = [run for cell in cells for paragraph in cell.paragraphs for run in paragraph.runs]
            for taxon in (
                "Salmo salar", "Cichorium intybus", "Brucella melitensis", "Brucella ovis",
                "Brucella suis", "Homo sapiens", "Hordeum vulgare", "Aspergillus nidulans",
                "B. melitensis", "B. ovis", "B. suis",
            ):
                matching = [run for run in table_runs if run.text == taxon]
                assert matching and all(run.font.italic is True for run in matching), taxon
            for gene in (
                "GAPDH", "EF1AA", "TUB", "ndhF", "ACE", "DCP1", "DEFB107B", "niaD",
                "gag", "mimp1", "impala", "Alu", "IS711", "IS6501",
            ):
                matching = [run for run in table_runs if run.text == gene]
                assert matching and all(run.font.italic is True for run in matching), gene


def test_optional_local_manuscript_and_supplement_word_packages_use_only_times_new_roman():
    """Every Word text role, theme fallback, bullet and equation must resolve to Times New Roman."""
    paths = [
        os.path.join(_ROOT, "manuscript", "manuscript.docx"),
        os.path.join(_ROOT, "manuscript", "manuscript_preprint.docx"),
        os.path.join(_ROOT, "manuscript", "tables", "table4_tool_comparison.docx"),
        os.path.join(_ROOT, "manuscript", "tables", "tableS1_tool_comparison_sourced.docx"),
    ]
    present = [path for path in paths if os.path.exists(path)]
    if not present:
        pytest.skip("local manuscript artifacts are not present")
    for path in present:
        with zipfile.ZipFile(path) as package:
            word_xml = b"\n".join(
                package.read(name) for name in package.namelist()
                if name.startswith("word/") and name.endswith(".xml")
            )
        rfonts = re.findall(rb'<w:rFonts\b([^>]*)/?>', word_xml)
        fonts = {
            value
            for attributes in rfonts
            for value in re.findall(
                rb'w:(?:ascii|hAnsi|cs|eastAsia)="([^"]+)"', attributes
            )
        }
        theme_fonts = set(re.findall(rb'typeface="([^"]*)"', word_xml))
        math_fonts = set(re.findall(rb'<m:mathFont\b[^>]*m:val="([^"]+)"', word_xml))
        assert fonts == {b"Times New Roman"}, (os.path.basename(path), fonts)
        assert theme_fonts <= {b"Times New Roman"}, (os.path.basename(path), theme_fonts)
        assert math_fonts <= {b"Times New Roman"}, (os.path.basename(path), math_fonts)
        assert "\uf0b7".encode("utf-8") not in word_xml, os.path.basename(path)


def test_optional_local_manuscript_artifacts_fully_justify_body_prose():
    """Both Word variants justify body prose while leaving captions and references readable."""
    from docx import Document
    from docx.enum.text import WD_ALIGN_PARAGRAPH

    for filename in ("manuscript.docx", "manuscript_preprint.docx"):
        path = os.path.join(_ROOT, "manuscript", filename)
        if not os.path.exists(path):
            pytest.skip(f"{filename} not present")
        document = Document(path)
        in_body = False
        in_references = False
        justified = 0
        references = 0
        for paragraph in document.paragraphs:
            if paragraph.text == "Abstract" and paragraph.style.name.startswith("Heading"):
                in_body = True
                continue
            if paragraph.text == "References" and paragraph.style.name.startswith("Heading"):
                in_references = True
                continue
            if not in_body or not paragraph.text.strip():
                continue
            is_heading = paragraph.style.name.startswith("Heading")
            is_list = paragraph.style.name == "List Bullet"
            is_caption = bool(re.match(r"^(Figure|Table) \d+\.", paragraph.text))
            is_figure = "<w:drawing" in paragraph._p.xml
            if is_heading or is_list or is_caption or is_figure:
                continue
            if in_references:
                assert paragraph.alignment == WD_ALIGN_PARAGRAPH.LEFT, (
                    filename, paragraph.text[:80]
                )
                references += 1
                continue
            assert paragraph.alignment == WD_ALIGN_PARAGRAPH.JUSTIFY, (
                filename, paragraph.text[:80]
            )
            justified += 1
        assert justified > 100, (filename, justified)
        assert references == 34, (filename, references)
        title = next(p for p in document.paragraphs if p.text.startswith("TEagle: an evidence-traceable"))
        assert title.alignment != WD_ALIGN_PARAGRAPH.JUSTIFY
        captions = [p for p in document.paragraphs
                    if re.match(r"^(Figure|Table) \d+\.", p.text)]
        assert captions
        assert all(p.alignment == WD_ALIGN_PARAGRAPH.LEFT for p in captions)


def test_optional_local_manuscript_artifacts_are_double_spaced_throughout():
    """Mobile DNA asks for double spacing throughout the submitted manuscript."""
    from docx import Document

    for filename in ("manuscript.docx", "manuscript_preprint.docx"):
        path = os.path.join(_ROOT, "manuscript", filename)
        if not os.path.exists(path):
            pytest.skip(f"{filename} not present")
        document = Document(path)
        for paragraph in document.paragraphs:
            if not paragraph.text.strip():
                continue
            spacing = paragraph.paragraph_format.line_spacing
            if spacing is None:
                spacing = paragraph.style.paragraph_format.line_spacing
            assert spacing == 2.0, (filename, paragraph.text[:80], spacing)


def test_reviewer_corpus_scope_is_an_exact_mirror():
    """The reviewer copy must not silently diverge from the benchmark corpus-scope statement."""
    assert _read("reviewer_package/corpora/CORPUS_SCOPE.md") == _read("benchmarks/CORPUS_SCOPE.md")


def test_report_benchmark_spec_is_an_exact_mirror_when_maintained():
    """The optional, gitignored report mirror is exact when present in a local working tree."""
    report_spec = os.path.join(_ROOT, "report", "BENCHMARK_SPEC_v1.md")
    if not os.path.exists(report_spec):
        pytest.skip("report benchmark-spec mirror not present")
    assert _read("report/BENCHMARK_SPEC_v1.md") == _read("benchmarks/protocol/BENCHMARK_SPEC_v1.md")


def test_domains_tested_string_covers_every_emitted_domain_code():
    """DOMAINS_TESTED is the scope statement shown beside every 'not detected'. If a code can be emitted
    but is absent from that sentence, the honesty invariant is false."""
    from teagle_core import classify
    described = classify.DOMAINS_TESTED.lower()
    emitted = {info[0] for info in domains.DOMAIN_INFO.values()}
    expected_words = {"ORF1": "orf1", "EN": "endonuclease", "YR": "tyrosine recombinase",
                      "HEL": "helitron", "RT": "rt", "INT": "integrase", "ENV": "envelope",
                      "GAG": "gag", "PR": "protease", "RNaseH": "rnase h", "CHR": "chromodomain",
                      "TPase": "transposase"}
    for code in emitted:
        word = expected_words.get(code)
        assert word, f"domain code {code} has no expected wording — extend this test"
        assert word in described, f"{code} is emitted but not named in DOMAINS_TESTED"


# ---------------- curated assembly table ----------------
def test_every_curated_assembly_entry_is_well_formed():
    """Shape check, offline: a malformed accession or taxid would fail only when a user picked that
    organism and waited for a multi-GB download to fail."""
    import re
    from teagle_core import fetch
    for organism, meta in fetch.COORD_ASSEMBLIES.items():
        assert re.match(r"^GC[AF]_\d+\.\d+$", meta["assemblyAccession"]), (organism, meta)
        assert meta["taxid"].isdigit(), (organism, meta)
        assert meta["assemblyName"].strip(), organism


@pytest.mark.network
def test_every_curated_accession_resolves_to_the_organism_it_claims():
    """An accession that looks right and resolves to a different species is worse than none: it would
    silently scan the wrong genome. Verified against NCBI Datasets."""
    import json
    import urllib.request
    from teagle_core import fetch
    for organism, meta in fetch.COORD_ASSEMBLIES.items():
        acc = meta["assemblyAccession"]
        url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{acc}/dataset_report"
        report = (json.load(urllib.request.urlopen(url, timeout=60)).get("reports") or [{}])[0]
        got = report.get("organism", {})
        assert str(got.get("tax_id")) == meta["taxid"], (organism, acc, got.get("tax_id"))


def test_methods_panel_text_is_derived_from_the_profile_table():
    """The in-app methods disclosure once claimed a 21-model Pfam panel after it had grown to 30 — the app
    under-reporting its own method, in a UI string no docs test covered. It is now derived from
    domains.DOMAIN_INFO, so the count and the Pfam list cannot disagree with what the scan loads."""
    import os as _os, sys as _sys
    _os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    _native = _os.path.join(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))), "app", "native")
    if _native not in _sys.path:
        _sys.path.insert(0, _native)
    import pytest as _pt
    _pt.importorskip("PySide6")
    from PySide6.QtWidgets import QApplication
    QApplication.instance() or QApplication([])
    import main
    from teagle_core import domains
    html = main.MainWindow._domain_panel_html()
    assert f"{len(domains.DOMAIN_INFO)} models" in html
    for _hmm, (_code, _label, _cls, pfam) in domains.DOMAIN_INFO.items():
        assert pfam in html, f"{pfam} missing from the methods panel text"


def test_tsd_method_text_states_the_window_the_detector_actually_searches():
    """Both TSD method strings claimed "4-12 bp" while find_tsd has searched from 2 bp since 3.2.1 — so the
    screen (and the exported Method column) told the user that the diagnostic 2 bp Tc1/Mariner TA, the exact
    case the detector was tuned to catch, could not be found. Both strings are now interpolated from
    structural.MIN_TSD/MAX_TSD; this fails if either drifts from the search bounds again."""
    import os as _os, sys as _sys
    _os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    _native = _os.path.join(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))), "app", "native")
    if _native not in _sys.path:
        _sys.path.insert(0, _native)
    import pytest as _pt
    _pt.importorskip("PySide6")
    from PySide6.QtWidgets import QApplication
    QApplication.instance() or QApplication([])
    import inspect
    import main
    from teagle_core import structural

    # the constants must be what find_tsd actually defaults to, or deriving from them proves nothing
    sig = inspect.signature(structural.find_tsd)
    assert sig.parameters["min_tsd"].default == structural.MIN_TSD
    assert sig.parameters["max_tsd"].default == structural.MAX_TSD

    window = f"{structural.MIN_TSD}–{structural.MAX_TSD} bp"
    assert window in main.MainWindow._STRUCT_METHOD["TSD"], main.MainWindow._STRUCT_METHOD["TSD"]

    # _methods_html reads only self._src_html (citation links), so a stub `self` exercises the real string
    # without standing up a MainWindow.
    class _Stub:
        _domain_panel_html = staticmethod(main.MainWindow._domain_panel_html)

        def _src_html(self, _key):
            return ""

    html = main.MainWindow._methods_html(_Stub())
    assert window in html, "the Methods & databases disclosure states a stale TSD window"


def test_domains_tested_names_every_transposase_group_the_panel_can_call():
    """The scope sentence named five transposase groups while the panel also bundles Pfam's generic DDE
    family (PF03184, class dna:DDE) — which classify names outright as "DDE transposon". A user reading
    that line was told a family the classifier CAN call is not tested, which is the not-detected-vs-
    not-tested confusion the line exists to prevent. The phrase is now derived from DOMAIN_INFO; this
    fails if a transposase class is added to the panel and the sentence does not follow."""
    from teagle_core import classify, domains

    panel = {cls for code, _l, cls, _p in domains.DOMAIN_INFO.values()
             if code == "TPase" and cls.startswith("dna:")}
    assert panel, "no dna: transposase classes in the panel — update this test if that is intended"
    assert len(classify._tpase_groups()) == len(panel), (classify._tpase_groups(), panel)
    for group in classify._tpase_groups():
        assert group in classify.DOMAINS_TESTED, f"{group} missing from the tested-panel sentence"


def test_public_release_text_has_no_automated_author_attribution():
    forbidden = tuple("".join(parts).lower() for parts in (
        ("clau", "de"), ("chat", "gpt"), ("open", "ai"), ("co", "dex"),
    ))
    for relpath in ("README.md", "CHANGELOG.md", "app/README.md", "THIRD-PARTY-NOTICES.md"):
        path = os.path.join(_ROOT, relpath)
        if not os.path.exists(path):
            continue
        lower = _read(relpath).lower()
        assert not any(term in lower for term in forbidden), relpath


def test_mobile_dna_manuscript_submission_guardrails():
    """Keep the submitted prose aligned with the explicit Mobile DNA requirements and scope limits."""
    text = _read("manuscript/manuscript.md")
    flat = " ".join(text.split())
    lower = flat.lower()
    for heading in (
        "## Background", "## Implementation", "## Results", "## Discussion", "## Conclusions",
        "## Availability and requirements", "## Declarations", "## Abbreviations",
    ):
        assert heading in text
    for declaration in (
        "Ethics approval and consent to participate", "Consent for publication",
        "Availability of data and materials", "Competing interests", "Funding",
        "Authors' contributions", "Acknowledgements",
    ):
        assert declaration in text
    acknowledgement = (
        "The author Tuna Birgün acknowledges financial support from the Scientific and Technological "
        "Research Council of Türkiye (TÜBİTAK) through the 2210-A National MSc Scholarship Programme "
        "(Application No. 1649B022502808)."
    )
    assert text.count(acknowledgement) == 1
    assert "**Acknowledgements.** Not applicable." not in text
    assert "The support disclosed in Acknowledgements was a personal scholarship to Tuna Birgün." in text
    assert "This research received no separate project-specific grant." in text
    assert "16. UCSC Genome Browser. In-Silico PCR. https://genome.ucsc.edu/cgi-bin/hgPcr " \
           "(accessed 18 August 2026)." in text
    assert "Abd el Wahab H" in text
    assert "Abd El Wahab H" not in text
    assert "GNU Affero General Public License v3.0 or later" in text
    assert text.count("Sibel Yılmaz") >= 2
    assert "Sibel Yilmaz" not in text
    notes_path = os.path.join(_ROOT, "manuscript", "SUBMISSION_NOTES.md")
    if os.path.exists(notes_path):
        notes = _read("manuscript/SUBMISSION_NOTES.md")
        assert "Sibel Yılmaz" in notes
        assert "Sibel Yilmaz" not in notes
    for item in (
        "Project name", "Project home page", "Archived version / persistent identifier",
        "Operating system", "Programming language", "Other requirements", "License",
        "Any restrictions to use by non-academics",
    ):
        assert f"*{item}:*" in text
    source_lines = text.splitlines()
    assert "Binding is anchored at the 3′ end rather than requiring the whole oligo to match: the longest" in source_lines
    assert (
        "- *Archived version / persistent identifier:* Every quantitative result was generated with TEagle "
        "`v{{corpus_run_version}}`, commit `642b271e2cd4d496a5277f92dd4abe9d6e8578b3`: "
        "https://github.com/tunabirgun/TEagle/releases/tag/v3.8.0. The version-locked benchmark evidence is archived "
        "at https://doi.org/10.5281/zenodo.21848814. A repository-independent software-source DOI has not been assigned."
    ) in source_lines
    assert any(
        "https://doi.org/10.5281/zenodo.21848814 contains curated corpus metadata" in line
        for line in source_lines
    )
    stale_archive_placeholders = (
        "Submission-time archive identifier required",
        "No DOI is claimed in this manuscript",
        "persistent identifier must be inserted",
        "submission archive will contain",
    )
    assert not any(phrase in text for phrase in stale_archive_placeholders)
    assert ("- *Operating system:* Windows 10/11, 64-bit. The application installs optional components "
            "in a WSL2 environment.") in source_lines
    assert "binding is anchored at the 3′ end rather than requiring the whole oligo to match: the longest" not in source_lines
    for forbidden in (
        "{{corpus_species}}", "{{holdout_species}}", "{{assay_species}}",
        "upper bound on performance", "proprietary", "benchmarks/", "benchmarks\\",
        "reviewer_package/", "reviewer_package\\", "STATISTICAL_ANALYSIS_PLAN.md",
        "repository test suite", "fails the build", "raw-output directory",
        "raw output by a script", "summary header", "version bump", "widely used",
        "A distinct problem remains unaddressed",
    ):
        assert forbidden not in text
    assert not re.search(r"\b(?:35|36|50)\s+species\b", text, re.I)
    assert not re.search(r"\b[\w.-]+\.(?:py|json|md)\b", text, re.I)
    assert not re.search(r"`[^`\n]*(?:[\\/]|\.(?:py|json|md))[^`\n]*`", text, re.I)
    for unsupported_claim in (
        "general utility", "broad utility", "establishes population-level accuracy",
        "population-level accuracy is", "population-level accuracy was",
    ):
        assert unsupported_claim not in lower
    assert "selection may bias the reported results in either direction" in lower
    assert "not population-level accuracy" in lower
    assert "integrated single-locus workflow" in lower
    assert "In the selected workflows, these steps are not assembled into a single workflow" in flat
    assert "A distinct problem is not covered by the selected workflows above" not in flat
    assert "Neither annotation count is a species count" in flat
    assert "primary independent unit is the distinct input SHA-256" in flat
    assert "assembled retrospectively from external records after excluding five bundled examples used during software development" in flat
    assert "no prospective preregistration is claimed" in flat
    assert "identification and parameter coverage only" in lower
    assert "cannot reconstruct a run by itself" in lower
    assert "not external biological ground truth" in lower
    assert ("does not replace prospective wet-lab validation of primers designed by teagle and "
            "does not establish de novo primer success") in lower
    assert "Mob DNA. 2026 Jul 16. Online ahead of print. doi:10.1186/s13100-026-00405-z" in text
    assert "Biology (Basel). 2026;15(1):59. doi:10.3390/biology15010059" in text
    assert "Row-level sensitivity of order accuracy by the confidence tier" in flat
    assert "Accession {{broad_duplicate_accession}} contributes a duplicated row" in flat
    # Corpus-level totals must not be presented as one performance denominator because controls and
    # element withholding have different scientific meanings.
    assert "{{broad_unique_order_answered}}" not in flat
    assert "{{broad_unique_order_withheld}}" not in flat
    assert re.search(r"Among the \{\{broad_element_unique_gradable\}\} element hashes, TEagle "
                     r"answered \{\{broad_element_unique_answered\}\} and withheld "
                     r"\{\{broad_element_unique_withheld\}\}", flat)

    spec = " ".join(_read("benchmarks/protocol/BENCHMARK_SPEC_v1.md").split()).lower()
    assert "selection can bias measured performance in either direction" in spec
    assert "upper bound" not in spec


def test_submission_assets_use_reader_facing_evidence_and_current_status():
    """Submission assets must not expose repository paths or repeat superseded review claims."""
    figure_sources = (
        "benchmarks/fig_workflow.py",
        "reviewer_package/scripts/fig_workflow.py",
        "manuscript/figures/figure-index.md",
        "manuscript/figures/fig_workflow.svg",
        "reviewer_package/figures/fig_workflow.svg",
    )
    for relpath in figure_sources:
        text = _read(relpath)
        assert "permissively licensed" not in text
    assert "Bundled — runs offline" in _read("manuscript/figures/fig_workflow.svg")

    table_paths = (
        "manuscript/tables/tableS1_tool_comparison_sourced.csv",
        "reviewer_package/results/tableS1_tool_comparison_sourced.csv",
    )
    forbidden_paths = ("app/backend/", "app/native/", "::_detector_parameters")
    for relpath in table_paths:
        text = _read(relpath)
        for forbidden in forbidden_paths:
            assert forbidden not in text
        assert "source and verification evidence supplied with the submission" in text

    for relpath in (
        "manuscript/tables/tableS1_tool_comparison_sourced.docx",
        "manuscript/tables/tableS1_tool_comparison_sourced.xlsx",
    ):
        with zipfile.ZipFile(os.path.join(_ROOT, relpath)) as package:
            payload = b"\n".join(package.read(name) for name in package.namelist())
        for forbidden in (b"app/backend/", b"app/native/", b"::_detector_parameters"):
            assert forbidden not in payload

    for relpath in (
        "manuscript/tables/table2_classification_accuracy.docx",
        "reviewer_package/results/table2_classification_accuracy.docx",
    ):
        with zipfile.ZipFile(os.path.join(_ROOT, relpath)) as package:
            payload = b"\n".join(package.read(name) for name in package.namelist())
        assert b"benchmarks/CORPUS_SCOPE.md" not in payload

    for relpath in (
        "manuscript/tables/table3_assay_validation.csv",
        "reviewer_package/results/table3_assay_validation.csv",
    ):
        header = _read(relpath).splitlines()[0]
        assert "Reference size (bp)" in header
        assert "Ground truth (bp)" not in header

    manuscript = " ".join(_read("manuscript/manuscript.md").split())
    assert "Clean-machine installation and uninstallation, and SmartScreen behaviour, were not tested" in manuscript
    assert "all required core dependencies are bundled" in manuscript
    assert "retrospectively replayed {{assay_cases}} published assays" in manuscript
    assert "without blinding" in manuscript
    assert "survive on a locked-down institutional laptop" not in manuscript
    assert "which capabilities are bundled and which require the optional WSL2 environment" in manuscript
    assert "Activation of tobacco retrotransposons during tissue culture" in manuscript
    assert "https://www.repeatmasker.org/ (accessed 7 August 2026" in manuscript

    reviewer = " ".join(_read("reviewer_package/README.md").split())
    assert "Phase 3 post-run" not in reviewer
    assert "Independent verification and submission readiness remain unclaimed" not in reviewer
    assert "benchmark evidence has persistent public access" in reviewer.lower()
