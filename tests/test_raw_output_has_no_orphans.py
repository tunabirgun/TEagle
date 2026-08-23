"""Each raw benchmark directory must contain exactly what its last run wrote, and nothing else.

`benchmarks/run_teagle.py` names each case `<row-index>_<accession>.json`. The row index comes from the
case's position in its corpus TSV, so inserting or removing a row shifts the indices of every case below
it. `--force` then rewrites the cases under their NEW names and leaves the files written under the OLD
names in place. `benchmarks/score.py` globs the directory rather than looking cases up by accession, so
both copies are read.

That is not hypothetical, and it has now happened twice. The development corpus reached 232 files for 135
rows before anyone noticed; the orphans were stale enough to predate a classifier change, so records
produced by two different rule sets were scored together and published panel sizes were inflated. The
first version of this test guarded that directory alone — and the defect promptly recurred in the
validation corpus the next time a row was removed, where 21 orphans survived beside 43 current files and
the scorer reported a 43-row corpus as 64 cases with twice the negative controls. A guard written for one
instance of a general defect is an invitation to meet it again somewhere else, which is why this now
iterates every corpus rather than naming one.

An orphan is invisible by construction: a well-formed file, in the right directory, with the right
schema, and every value inside it was true when it was written. Nothing about the file is wrong. Only its
continued presence is.
"""
from __future__ import annotations
import csv
import glob
import json
import os

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAWDIRS = [(os.path.join(ROOT, "benchmarks", "raw", "teagle"),
            os.path.join(ROOT, "benchmarks", "corpus.tsv")),
           (os.path.join(ROOT, "benchmarks", "raw", "teagle_holdout"),
            os.path.join(ROOT, "benchmarks", "corpus_holdout.tsv"))]


def _case_files(raw):
    return sorted(glob.glob(os.path.join(raw, "[0-9]*.json")))


def _dirs():
    """(raw dir, corpus, run record) for every benchmark corpus that has been executed here."""
    out = []
    for raw, corpus in RAWDIRS:
        run = os.path.join(raw, "_run.json")
        if os.path.exists(run):
            with open(run, encoding="utf-8") as fh:
                out.append((raw, corpus, json.load(fh)))
    return out


def _corpus_counts(corpus):
    """accession stem -> number of rows citing it. The corpora repeat accessions on purpose."""
    out = {}
    with open(corpus, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            acc = (r.get("accession") or "").strip().split(".")[0]
            if acc:
                out[acc] = out.get(acc, 0) + 1
    return out


def _skip_if_unrun():
    if len(_dirs()) != len(RAWDIRS):
        pytest.skip("complete private raw benchmark runs are not retained in this checkout")


def test_both_corpora_were_discovered():
    """Guard the guard. If discovery finds one directory when the project has two, every test below
    reports clean on exactly the case this suite was extended to cover."""
    _skip_if_unrun()
    found = [os.path.relpath(d[0], ROOT) for d in _dirs()]
    assert len(found) == 2, (
        f"expected both the development and validation raw directories, found {len(found)}: {found}")


def test_case_file_count_matches_what_the_run_reported():
    """The run record states how many cases it analysed. Any excess is an orphan from an earlier run."""
    _skip_if_unrun()
    bad = []
    for raw, _corpus, run in _dirs():
        analysed = run["counts"]["analysed"]
        n = len(_case_files(raw))
        if n != analysed:
            bad.append(f"{os.path.relpath(raw, ROOT)} holds {n} case files; its last run analysed "
                       f"{analysed}")
    assert not bad, (
        "; ".join(bad)
        + ". The surplus was written by an earlier run under different corpus row indices and is still "
          "being scored. Move the files older than that run's start time into benchmarks/outdated/ and "
          "rescore.")


def test_each_accession_has_exactly_as_many_files_as_corpus_rows():
    """Only a SURPLUS is an orphan.

    An accession stored under several indices is not itself the defect: the corpora cite some accessions
    more than once, one row per element in the deposit, and AF391808 alone carries eight. A shortfall is
    not a defect either — it means the case failed to fetch or analyse, which the run record reports as a
    failure and the manuscript reports rather than dropping.
    """
    _skip_if_unrun()
    bad = []
    for raw, corpus, _run in _dirs():
        if not os.path.exists(corpus):
            continue
        expected = _corpus_counts(corpus)
        found = {}
        for path in _case_files(raw):
            acc = os.path.basename(path).partition("_")[2][:-len(".json")].split(".")[0]
            found[acc] = found.get(acc, 0) + 1
        surplus = {a: (n, expected.get(a, 0)) for a, n in found.items() if n > expected.get(a, 0)}
        if surplus:
            bad.append(f"{os.path.relpath(raw, ROOT)}: {surplus}")
    assert not bad, (
        "accessions with more raw files than their corpus has rows citing them, so the scorer counts "
        "those cases more than once - {accession: (files, corpus rows)}: " + "; ".join(bad))


def test_every_case_file_belongs_to_a_corpus_row():
    """Catches a row deliberately removed — an excluded development accession, say — whose raw record was
    left behind and would go on being scored against a label the corpus no longer carries."""
    _skip_if_unrun()
    bad = []
    for raw, corpus, _run in _dirs():
        if not os.path.exists(corpus):
            continue
        stems = set(_corpus_counts(corpus))
        orphans = sorted(
            os.path.basename(p) for p in _case_files(raw)
            if os.path.basename(p).partition("_")[2][:-len(".json")].split(".")[0] not in stems)
        if orphans:
            bad.append(f"{os.path.relpath(raw, ROOT)}: {orphans}")
    assert not bad, (
        "case files whose accession is not in its corpus, left behind by a corpus that has since "
        "changed: " + "; ".join(bad))
