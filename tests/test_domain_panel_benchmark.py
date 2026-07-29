"""Domain-panel benchmark: every bundled Pfam model must fire on a real specimen of the element it models.

A profile that never fires is not coverage — it is an unfalsifiable claim in the scope statement. Each
model added on 2026-07-28 was validated against a named specimen before being enabled, and this pins that
validation so a panel edit, a threshold change, or a pyhmmer upgrade cannot silently break it.

The E-value floors are deliberately loose (an order of magnitude or more below what was measured): they
catch a model going dark, not normal variation. Measured values on 2026-07-28 are recorded beside each
row so a drift is visible rather than merely a pass.
"""
import os
import sys

import pytest

_BE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "app", "backend")
if _BE not in sys.path:
    sys.path.insert(0, _BE)
from teagle_core import domains, fetch                           # noqa: E402

# (accession, specimen, hmm profile that must fire, loose E-value ceiling, E-value measured 2026-07-28)
PANEL_BENCH = [
    ("M80343", "Human LINE-1 L1.2",              "Transposase_22",  1e-30, 4.3e-68),
    ("M80343", "Human LINE-1 L1.2",              "Tnp_22_trimer",   1e-15, 2.1e-35),
    ("M80343", "Human LINE-1 L1.2",              "Tnp_22_dsRBD",    1e-20, 3.6e-51),
    ("M80343", "Human LINE-1 L1.2",              "Exo_endo_phos",   1e-10, 1.6e-23),
    ("M11340", "Dictyostelium DIRS-1",           "Phage_integrase", 1e-04, 1.6e-08),
    ("AB267078", "Ipomoea Hel-It1 (Helitron)",   "Helitron_like_N", 1e-20, 6.1e-49),
    ("M25427", "Maize En-1 (CACTA)",             "Transposase_24",  1e-05, 3.1e-10),
    ("M76978", "Maize MuDR (MULE)",              "MULE",            1e-08, 1.4e-17),
    ("OQ718454", "Bactrocera piggyBac-like",     "DDE_Tnp_1_7",     1e-20, 1.8e-44),
    ("M11240", "Drosophila copia",               "RVT_2",           1e-30, 1.7e-78),
    ("M11240", "Drosophila copia",               "rve",             1e-05, 2.0e-11),
    ("M12927", "Drosophila gypsy",               "RVT_1",           1e-20, 2.2e-42),
    ("X01005", "C. elegans Tc1",                 "HTH_Tnp_Tc3_2",   1e-15, 1.1e-26),
    ("X05424", "Maize Ac",                       "Dimer_Tnp_hAT",   1e-20, 4.9e-35),
]

# Specificity: the LINE modules must NOT fire on elements that are not LINEs. Before ORF1p/EN were added
# a full-length L1 and a dead fragment scored identically; adding them must not buy that at the cost of
# calling every retroelement a LINE.
LINE_ONLY = {"Transposase_22", "Tnp_22_trimer", "Tnp_22_dsRBD"}
NON_LINE_SPECIMENS = ["M11240", "M12927", "X01005", "X05424"]


def _hits(accession):
    """Best (lowest-E) hit per profile. A plain dict comprehension keeps the LAST occurrence instead,
    which reported a weaker duplicate hit and made a healthy model look near its threshold."""
    seq = "".join(l for l in fetch.retrieve(accession)["fasta"].splitlines() if not l.startswith(">"))
    best = {}
    for h in domains.scan_domains(seq):
        cur = best.get(h["hmm"])
        if cur is None or (h.get("evalue") is not None and h["evalue"] < cur.get("evalue", float("inf"))):
            best[h["hmm"]] = h
    return best


@pytest.mark.network
@pytest.mark.parametrize("acc,specimen,profile,ceiling,measured", PANEL_BENCH)
def test_every_panel_model_fires_on_its_specimen(acc, specimen, profile, ceiling, measured):
    hits = _hits(acc)
    assert profile in hits, f"{profile} did not fire on {specimen} ({acc}) — the model is dark"
    ev = hits[profile].get("evalue")
    if ev is not None:
        assert ev <= ceiling, f"{profile} on {specimen}: E={ev:.1e} exceeds the {ceiling:.0e} floor"


@pytest.mark.network
@pytest.mark.parametrize("acc", NON_LINE_SPECIMENS)
def test_line_specific_models_do_not_fire_on_non_lines(acc):
    fired = LINE_ONLY & set(_hits(acc))
    assert not fired, f"LINE-specific model(s) {fired} fired on non-LINE specimen {acc}"


@pytest.mark.network
def test_full_length_line_reaches_intact_and_a_fragment_does_not():
    """The defect this panel expansion existed to fix: before ORF1p and the ORF2p endonuclease were
    modelled, a full-length L1 and a 5'-truncated fragment returned byte-identical completeness."""
    from teagle_core import classify, structural
    seq = "".join(l for l in fetch.retrieve("M80343")["fasta"].splitlines() if not l.startswith(">"))

    def verdict(s):
        return classify.classify(structural.detect_all(s), domains.scan_domains(s))["completeness"]

    full, frag = verdict(seq), verdict(seq[-3000:])
    assert full != frag
    assert full["missing"] == [] and "intact" in full["tier"]
    assert frag["missing"], "a 5'-truncated LINE must record something missing"


def test_panel_size_matches_the_documented_scope():
    """DOMAINS_TESTED is shown to the user beside every 'not detected'; it must not describe a panel
    larger or smaller than the one that ran."""
    from teagle_core import classify
    assert len(domains._hmms()) == len(domains.DOMAIN_INFO), \
        "every bundled profile needs a DOMAIN_INFO entry or it reports under its raw HMM name"
    assert classify.DOMAINS_TESTED.strip()
