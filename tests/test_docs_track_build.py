"""Documentation claims that must be derived from the build, not restated by hand.

A number written into prose drifts the moment the thing it describes changes: the README claimed a
"21-model Pfam panel" after the panel had grown to 30. These tests fail when a stated figure and the
artefact it describes disagree, so the rot is caught by the suite rather than by a reader.
"""
import io
import os
import re

import pytest

from teagle_core import domains

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _readme():
    with io.open(os.path.join(_ROOT, "README.md"), encoding="utf-8") as fh:
        return fh.read()


def test_readme_panel_size_matches_the_bundled_profiles():
    stated = re.search(r"(\d+)-model Pfam panel", _readme())
    assert stated, "README no longer states the panel size — update this test if that is intended"
    assert int(stated.group(1)) == len(domains._hmms())


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
