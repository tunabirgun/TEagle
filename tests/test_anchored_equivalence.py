"""The single-pass 3'-anchored matcher must keep returning what the descending scan returned.

`primers.anchored_sites` replaced a literal implementation that rescanned the template once per candidate
core length. The replacement is a performance rewrite, so its whole justification is that it changes
nothing observable; the full measurement is benchmarks/diff_anchored_sites.py, which runs 6,000 seeded
trials and writes its result to benchmarks/raw/diff_anchored_sites.json.

This runs a small sample of the same comparison inside the ordinary suite, so a future edit to the
matcher fails here rather than only in a benchmark someone has to remember to run. It is deliberately
cheap and deliberately not a replay: a different seed, so a change that happens to survive the recorded
run still has to survive fresh trials.

Half the sample is drawn from the regime where the strict 3' window is wider than the shortest permitted
core, which is where the first attempt at the replacement failed and which a uniform draw over the stated
parameter ranges reaches about 2% of the time. The harness's own sensitivity control — a reference that
clamps the strict window once instead of per core length — is asserted to be caught, so a sample that
could not detect that defect fails instead of passing quietly.
"""
import os
import sys

# appended, never inserted: the benchmarks directory carries ordinary module names (score, figures,
# testsuite) and must not shadow anything the rest of the suite imports
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "benchmarks")))

import diff_anchored_sites as diff                         # noqa: E402

TRIALS = 300
SEED = 8112026                     # not the benchmark's seed: this is a fresh sample, not a replay


def test_single_pass_matcher_matches_the_descending_scan():
    r = diff.run(trials=TRIALS, seed=SEED, strict_share=0.5)

    assert r["disagreements"] == 0, r["examples"][:3]

    # A trial set that produced no sites would agree trivially, and so would one that never reached the
    # regime the clamps can differ in. Both are asserted, so the zero above cannot be vacuous.
    assert r["site_yield_ok"], r["informative"]
    assert r["informative"]["trials_where_the_winning_core_was_shorter_than_the_primer"] > 0
    sc = r["sensitivity_control"]
    assert sc["trials_caught"] > 0, sc
    assert r["site_cap_hits"] == 0, "a capped site scan is not comparable with the uncapped reference"
