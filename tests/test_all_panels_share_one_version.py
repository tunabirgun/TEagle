"""Keep locally present benchmark panels internally version-consistent.

Application releases and sanitized benchmark-evidence deposits are controlled separately. The default
tool suite therefore checks agreement among locally present panels without treating them as application
release assets. A benchmark release must also set ``TEAGLE_VERIFY_BENCHMARK_RELEASE=1``; that explicit
gate requires every panel to match the application source in the checkout being deposited.
"""
from __future__ import annotations
import glob
import json
import os

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "benchmarks", "raw")
SOURCE = os.path.join(ROOT, "app", "backend", "teagle_core", "__init__.py")


def _shipped_version():
    import re
    with open(SOURCE, encoding="utf-8") as fh:
        m = re.search(r'__version__\s*=\s*["\']([^"\']+)["\']', fh.read())
    assert m, "could not parse __version__ from teagle_core/__init__.py"
    return m.group(1)


def _panel_versions():
    """(file, version) for every raw artefact that records the version that produced it."""
    out = []
    if not os.path.isdir(RAW):
        return out
    paths = sorted(glob.glob(os.path.join(RAW, "*.json")))
    paths += sorted(glob.glob(os.path.join(RAW, "*", "_run.json")))
    for p in paths:
        try:
            with open(p, encoding="utf-8") as fh:
                d = json.load(fh)
        except (json.JSONDecodeError, OSError):
            continue
        if isinstance(d, dict) and d.get("teagle_version"):
            out.append((os.path.relpath(p, ROOT), d["teagle_version"]))
    return out


PANELS = _panel_versions()


def test_panels_were_found():
    """Guard the guard. An empty list would make every assertion below pass without checking anything,
    and this test exists precisely because a check that inspects a subset is worse than none."""
    if not os.path.isdir(RAW):
        pytest.skip("no benchmarks/raw — benchmarks have not been run in this checkout")
    assert len(PANELS) >= 4, (
        f"only {len(PANELS)} raw artefacts record a teagle_version; the discovery rule has gone stale "
        f"against the current layout, so this test is not checking what it claims to check")


def test_every_panel_shares_one_version():
    versions = {}
    for path, ver in PANELS:
        versions.setdefault(ver, []).append(path)
    assert len(versions) == 1, (
        "benchmark panels were produced by more than one TEagle version, so the manuscript cannot "
        "report a single version as describing all of its results. Re-run the panels that lag:\n"
        + "\n".join(f"  {v}: {', '.join(f)}" for v, f in sorted(versions.items())))


def test_panels_match_the_shipped_source():
    """The explicit benchmark-release gate binds all local panels to the checked-out source."""
    if os.environ.get("TEAGLE_VERIFY_BENCHMARK_RELEASE") != "1":
        pytest.skip("tool-only release; set TEAGLE_VERIFY_BENCHMARK_RELEASE=1 for the benchmark gate")
    shipped = _shipped_version()
    stale = sorted({v for _p, v in PANELS if v != shipped})
    assert not stale, (
        f"the source ships {shipped} but benchmark panels record {stale}. The release would carry "
        f"results describing a different build. Re-run the benchmarks, or bump the version after them.")
