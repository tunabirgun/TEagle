"""The raw benchmark directory must contain exactly what the last run wrote, and nothing else.

`benchmarks/run_teagle.py` names each case `<row-index>_<accession>.json`. The row index comes from the
case's position in `corpus.tsv`, so inserting or removing a corpus row shifts the indices of every case
below it. `--force` then rewrites the cases under their NEW names and leaves the files written under the
OLD names in place. `benchmarks/score.py` globs the directory rather than looking cases up by accession,
so both copies are read.

That is not hypothetical. The directory reached 232 files for 135 corpus rows before anyone noticed, and
the orphans were stale enough to predate a classifier fix — so records produced by two different versions
of the decision rules were scored together, the same case was counted more than once, and published panel
sizes were inflated. The symptom was subtle: duplicated lines in the scorer's error list, which reads as
two similar cases rather than one case counted twice.

An orphan is invisible by construction: it is a well-formed file, in the right directory, with the right
schema, and every value inside it was true when it was written. Nothing about the file is wrong. Only its
continued presence is. That is exactly the sort of defect that has to be caught mechanically.
"""
from __future__ import annotations
import csv
import glob
import json
import os

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "benchmarks", "raw", "teagle")
RUN = os.path.join(RAW, "_run.json")
CORPUS = os.path.join(ROOT, "benchmarks", "corpus.tsv")


def _case_files():
    return sorted(glob.glob(os.path.join(RAW, "[0-9]*.json")))


def _run_record():
    if not os.path.exists(RUN):
        pytest.skip("no benchmarks/raw/teagle/_run.json — the corpus has not been run in this checkout")
    with open(RUN, encoding="utf-8") as fh:
        return json.load(fh)


def test_case_file_count_matches_what_the_run_reported():
    """The run record states how many cases it analysed. Any excess is an orphan from an earlier run."""
    run = _run_record()
    analysed = run["counts"]["analysed"]
    files = _case_files()
    assert len(files) == analysed, (
        f"benchmarks/raw/teagle holds {len(files)} case files but the last run analysed {analysed}. "
        f"The surplus was written by an earlier run under different corpus row indices and is still being "
        f"scored. Move the files older than the run's start time out of the directory — "
        f"benchmarks/outdated/ is where superseded output goes — and rescore.")


def test_each_accession_has_exactly_as_many_files_as_corpus_rows():
    """The orphans take the shape 035_AF018167.json beside 036_AF018167.json — but an accession under
    several indices is NOT itself the defect, because the corpus repeats accessions on purpose: 135 rows
    cover 122 distinct accessions, and AF391808 alone carries eight, one per element in the deposit.

    So the invariant is not "one file per accession"; it is that the number of files for an accession
    equals the number of corpus rows citing it. Anything above that count is a file the current corpus
    cannot account for."""
    if not os.path.exists(CORPUS):
        pytest.skip("no corpus.tsv")
    with open(CORPUS, encoding="utf-8") as fh:
        expected = {}
        for r in csv.DictReader(fh, delimiter="\t"):
            acc = (r.get("accession") or "").strip().split(".")[0]
            if acc:
                expected[acc] = expected.get(acc, 0) + 1

    found = {}
    for path in _case_files():
        acc = os.path.basename(path).partition("_")[2][:-len(".json")].split(".")[0]
        found[acc] = found.get(acc, 0) + 1

    # Only a SURPLUS is an orphan. A shortfall means the case failed to fetch or analyse, which the run
    # record already reports as a failure and the manuscript reports rather than dropping.
    surplus = {a: (n, expected.get(a, 0)) for a, n in found.items() if n > expected.get(a, 0)}
    assert not surplus, (
        "accessions with more raw files than the corpus has rows for them, so the scorer counts those "
        f"cases more than once — {{accession: (files, corpus rows)}}: {surplus}")


def test_every_case_file_belongs_to_a_corpus_row():
    """An accession no longer in corpus.tsv cannot be scored against a label, so it should not be here."""
    if not os.path.exists(CORPUS):
        pytest.skip("no corpus.tsv")
    with open(CORPUS, encoding="utf-8") as fh:
        corpus = {r["accession"].strip() for r in csv.DictReader(fh, delimiter="\t") if r.get("accession")}
    # The runner strips a version suffix from some accessions when naming files, so compare on the stem.
    stems = {a.split(".")[0] for a in corpus}
    orphans = sorted(
        os.path.basename(p) for p in _case_files()
        if os.path.basename(p).partition("_")[2][:-len(".json")].split(".")[0] not in stems)
    assert not orphans, (
        f"case files whose accession is not in corpus.tsv, so they were left by a corpus that has since "
        f"changed: {orphans}")
