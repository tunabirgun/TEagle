"""A numeric constant that decides a call must appear in the provenance manifest, or be provably unused.

`engine._detector_parameters()` reads each threshold from `inspect.signature` of the function that applies
it, so a default changed in a detector changes the seal automatically and the manifest cannot drift from
the code. That idiom is good and it has one blind spot: it can only see values that are *signature
defaults*. A threshold defined at module scope and referenced directly from another module is invisible
to it, and invisible in a way no audit keyed on detector functions will report, because such an audit
enumerates the functions and never learns the constant exists.

Three constants sat outside it. `MIN_TIR_STANDALONE` is the length at or above which a terminal inverted
repeat licenses a Class II call with no transposase — the one bound this project measures a chance
ceiling against, and therefore the single threshold a reader is most likely to want to reconstruct.
`_LTR_SEED_WINDOW` caps the region searched for terminal-repeat seeds, so on a long record it decides
what can be found at all. Both were applied and neither was sealed, while the manuscript claimed the
manifest reconstructs every bound applied and the pre-registered specification set "manifest sufficient
to reconstruct every threshold" as a gate.

The third, `_MIN_LINE_RECORD_BP`, is genuinely unused: its own comment records that record length is not
used to license a LINE call. Sealing it would be its own dishonesty — a manifest entry implying a bound
that decides nothing. So the rule this test enforces is not "seal every constant" but "seal every
constant that is read outside the module that defines it, and keep the unused ones on an explicit list
that says why they are exempt."
"""
from __future__ import annotations
import os
import re
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CORE = os.path.join(ROOT, "app", "backend", "teagle_core")
sys.path.insert(0, os.path.join(ROOT, "app", "backend"))

# Constants that are defined but do not decide anything, each with the reason it is exempt. A constant may
# only be added here with evidence that nothing reads it — not because sealing it is inconvenient.
KNOWN_UNUSED = {
    "_MIN_LINE_RECORD_BP": "retained for reference; record length does not license a LINE call, and the "
                           "module comment says so. Read by nothing.",
}


def _numeric_module_constants(path):
    """Module-scope numeric constants, by name, from source rather than by import, so a constant is found
    whether or not importing the module happens to expose it."""
    with open(path, encoding="utf-8") as fh:
        src = fh.read()
    out = {}
    for m in re.finditer(r"^([A-Z_][A-Z0-9_]*|_[A-Za-z][A-Za-z0-9_]*)\s*=\s*(-?\d+(?:\.\d+)?)\s*(?:#.*)?$",
                         src, re.M):
        out[m.group(1)] = m.group(2)
    return out


def _read_outside_defining_module(name, defining):
    """Is this constant referenced from any module other than the one that defines it?"""
    for fn in os.listdir(CORE):
        if not fn.endswith(".py") or fn == defining:
            continue
        with open(os.path.join(CORE, fn), encoding="utf-8") as fh:
            body = fh.read()
        # Strip comments so a constant merely discussed in prose does not count as applied.
        body = re.sub(r"#.*$", "", body, flags=re.M)
        if re.search(rf"\b{re.escape(name)}\b", body):
            return fn
    return None


def test_every_cross_module_threshold_is_sealed_or_declared_unused():
    import engine
    sealed_values = {str(v) for v in engine._detector_parameters().values()}
    sealed_names = {k.lower().replace("_", "") for k in engine._detector_parameters()}

    unsealed = []
    for fn in ("structural.py", "classify.py", "domains.py", "sequtil.py"):
        path = os.path.join(CORE, fn)
        if not os.path.exists(path):
            continue
        for name, value in _numeric_module_constants(path).items():
            if name in KNOWN_UNUSED:
                continue
            user = _read_outside_defining_module(name, fn)
            if user is None:
                continue                       # defined and used only where it is defined
            key = name.lower().lstrip("_").replace("_", "")
            if key in sealed_names or any(key in s for s in sealed_names) or value in sealed_values:
                continue
            unsealed.append(f"{fn}::{name} = {value}, applied in {user}")

    assert not unsealed, (
        "numeric thresholds are applied across module boundaries but do not appear in the provenance "
        "manifest, so a run cannot be reconstructed from its own seal:\n  "
        + "\n  ".join(unsealed)
        + "\nAdd each to engine._detector_parameters(), or add it to KNOWN_UNUSED with evidence that "
          "nothing reads it.")


def test_known_unused_constants_really_are_unused():
    """Guard the exemption list. An entry that starts being read would otherwise stay exempt forever,
    which turns the escape hatch into exactly the blind spot this test exists to close."""
    still_unused = []
    for name, reason in KNOWN_UNUSED.items():
        for fn in os.listdir(CORE):
            if not fn.endswith(".py"):
                continue
            with open(os.path.join(CORE, fn), encoding="utf-8") as fh:
                body = re.sub(r"#.*$", "", fh.read(), flags=re.M)
            # A single hit is the definition itself; two or more means something reads it.
            if len(re.findall(rf"\b{re.escape(name)}\b", body)) > 1:
                still_unused.append(f"{name} is exempt as unused ({reason}) but is read in {fn}")
    assert not still_unused, "\n  ".join(still_unused)
