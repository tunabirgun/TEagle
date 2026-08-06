"""The reviewer package must be a copy, not a fork.

`reviewer_package/` exists so a reviewer can check a number without navigating the working tree, and its
README states the arrangement plainly: nothing in it is unique to the folder, and its contents correspond
to the manuscript as submitted. Both sentences are promises about byte equality with files that live
elsewhere, and nothing enforced them.

They drifted the first time a generator changed. `benchmarks/make_tables.py` was edited to fix a truncated
table cell; the manuscript's Table 4 was regenerated; the package's copy of the script and of the table
were not. The shipped package then contained a script that, run as its own README instructs, reproduces a
table different from the one in the paper — on the very table a reviewer would open it to check. Four cells
disagreed, one of them visibly mangled.

A copy that can silently stop being a copy is worse than no copy, because the README asserts it is one.
These tests fail the moment the equality the package claims stops holding.
"""
from __future__ import annotations
import filecmp
import os

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PKG = os.path.join(ROOT, "reviewer_package")
SCRIPTS = os.path.join(PKG, "scripts")
RESULTS = os.path.join(PKG, "results")

# Where each copied file's original lives. Directories rather than a file list, so a file added to the
# package later is covered without anyone remembering to extend this.
SOURCE_DIRS = {
    SCRIPTS: [os.path.join(ROOT, "benchmarks")],
    RESULTS: [os.path.join(ROOT, "benchmarks", "raw"), os.path.join(ROOT, "manuscript", "tables")],
}


def _copied_files(pkg_dir):
    """(package copy, original) for every file in pkg_dir that has an original elsewhere.

    A file with no counterpart is not an error here — `corpora/` and the README are assembled for the
    package and have no single upstream — but a file that DOES have one must match it.
    """
    if not os.path.isdir(pkg_dir):
        return []
    out = []
    for name in sorted(os.listdir(pkg_dir)):
        copy = os.path.join(pkg_dir, name)
        if not os.path.isfile(copy):
            continue
        for src_dir in SOURCE_DIRS[pkg_dir]:
            original = os.path.join(src_dir, name)
            if os.path.isfile(original):
                out.append((copy, original))
                break
    return out


ALL_PAIRS = _copied_files(SCRIPTS) + _copied_files(RESULTS)


def test_the_package_actually_copies_something():
    """Guard the guard. If the discovery above silently found nothing — a directory renamed, the package
    restructured — every parametrised test below would vacuously pass and the check would be gone."""
    assert len(ALL_PAIRS) >= 10, (
        f"only {len(ALL_PAIRS)} reviewer-package files were matched to an original; the discovery rule "
        f"in {os.path.basename(__file__)} has probably gone stale against the current layout")


@pytest.mark.parametrize("copy,original", ALL_PAIRS,
                         ids=[os.path.basename(c) for c, _o in ALL_PAIRS])
def test_reviewer_package_copy_matches_its_source(copy, original):
    rel_c = os.path.relpath(copy, ROOT)
    rel_o = os.path.relpath(original, ROOT)
    assert filecmp.cmp(original, copy, shallow=False), (
        f"{rel_c} has drifted from {rel_o}. The reviewer package promises to be a copy, so regenerate "
        f"or re-copy it in the same commit as the change that made them differ.")


def test_the_scripts_the_readme_names_are_present():
    """The README lists the commands that reproduce `results/` from `corpora/`. A command naming a script
    the package does not ship sends a reviewer to a file that is not there."""
    readme = os.path.join(PKG, "README.md")
    if not os.path.exists(readme):
        pytest.skip("no reviewer_package/README.md")
    with open(readme, encoding="utf-8") as fh:
        text = fh.read()
    named = {line.split("benchmarks/", 1)[1].split()[0]
             for line in text.splitlines() if "python benchmarks/" in line}
    have = set(os.listdir(SCRIPTS)) if os.path.isdir(SCRIPTS) else set()
    missing = sorted(named - have)
    assert not missing, f"reviewer_package/README.md names scripts the package does not ship: {missing}"
