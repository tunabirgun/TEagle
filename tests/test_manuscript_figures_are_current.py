"""Every figure the manuscript embeds must be byte-identical to what its generator last produced.

The generators write to `benchmarks/figures/`. The manuscript build embeds from `manuscript/figures/`.
Nothing connected the two, so regenerating a figure updated the analysis copy and left the submitted copy
untouched — and the failure is invisible from the Markdown, which carries only a filename.

It reached a submission. The confidence-calibration figure embedded in `manuscript_preprint.docx` plotted
the Moderate tier at n = 22 and roughly 0.955 accuracy, from a run that predated a corpus de-duplication,
while the caption printed directly beneath it read n = 20 and 1.000 from the current scored output. A
referee extracted the embedded raster, hashed it, and found it identical to the superseded file. That
figure was the paper's flagship negative result, so the one asset a reader would scrutinise most closely
was the one that disagreed with its own caption.

Modification time cannot be trusted for this: copying, checkout and archive extraction all rewrite it
without changing content, and a regenerated figure that happens to be identical is not a defect. The
comparison is therefore on content.
"""
from __future__ import annotations
import glob
import hashlib
import os

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GENERATED = os.path.join(ROOT, "benchmarks", "figures")
EMBEDDED = os.path.join(ROOT, "manuscript", "figures")


def _sha(path):
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def _pairs():
    """(name, generated, embedded) for every figure asset present in both directories."""
    if not (os.path.isdir(GENERATED) and os.path.isdir(EMBEDDED)):
        return []
    out = []
    for src in sorted(glob.glob(os.path.join(GENERATED, "fig_*.*"))):
        if src.endswith(".json"):          # summary values, not an asset the manuscript embeds
            continue
        dst = os.path.join(EMBEDDED, os.path.basename(src))
        if os.path.exists(dst):
            out.append((os.path.basename(src), src, dst))
    return out


PAIRS = _pairs()


def test_there_are_figures_to_compare():
    """Guard the guard: if the discovery returns nothing, the parametrised test below would pass
    vacuously and the check would be silently gone."""
    if not os.path.isdir(EMBEDDED):
        pytest.skip("no manuscript/figures — manuscript is gitignored and absent from this checkout")
    assert len(PAIRS) >= 3, (
        f"only {len(PAIRS)} figure assets matched between benchmarks/figures and manuscript/figures; "
        f"the manuscript has three figures in three formats each, so the discovery rule has gone stale")


@pytest.mark.parametrize("name,generated,embedded", PAIRS, ids=[p[0] for p in PAIRS])
def test_embedded_figure_matches_its_generator_output(name, generated, embedded):
    g, e = _sha(generated), _sha(embedded)
    assert g == e, (
        f"manuscript/figures/{name} is not what the generator last produced. The manuscript would embed "
        f"a superseded figure while its caption prints current numbers. Copy benchmarks/figures/{name} "
        f"over it and rebuild the renditions.\n"
        f"  generated {g[:16]}  {os.path.getsize(generated):,} B\n"
        f"  embedded  {e[:16]}  {os.path.getsize(embedded):,} B")
