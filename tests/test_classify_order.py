"""Regression tests for the strand-aware INT-vs-RT ordering that decides Copia vs Gypsy.
Synthetic domain hits place INT and RT in different ORFs / strands without a network fetch —
the split-ORF case the golden fixtures cannot exercise (their INT+RT share one ORF)."""
from teagle_core import classify

LTR = [{"type": "LTR (276 bp, 98% id)"}]


def _dom(code, nt, strand="+", score=50.0):
    return {"domain": code, "nt": list(nt), "strand": strand, "score": score,
            "aa": [1, max(1, (nt[1] - nt[0]) // 3)], "class": "retro"}


def test_copia_when_integrase_upstream_plus_strand():
    # INT genomically upstream of RT on + strand -> Copia, regardless of ORF length rank
    cl = classify.classify(LTR, [_dom("INT", (100, 300), "+", 40), _dom("RT", (500, 1400), "+", 90)])
    assert cl["superfamily"].startswith("Copia"), cl["superfamily"]
    assert cl["confidence"] == "High"


def test_protease_evidence_recorded_from_pr_code():
    # domains.py emits code "PR" for the aspartic protease (RVP profile); the evidence line must fire on "PR"
    cl = classify.classify(LTR, [_dom("INT", (100, 300), "+", 40), _dom("RT", (500, 1400), "+", 90),
                                 _dom("PR", (350, 480), "+", 30)])
    assert any("aspartic-protease" in e for e in cl["evidence"]), cl["evidence"]


def test_line_polyt_not_mislabeled_as_polya():
    # a LINE with only a 5' poly-T tract must not be reported as a "3' poly-A tail"
    struct = [{"type": "poly-T (5')", "length": 20, "pos": [0, 20]}]
    cl = classify.classify(struct, [_dom("RT", (100, 1000), "+", 90)])
    ev = " ".join(cl["evidence"])
    assert "poly-T" in ev and "poly-A tail consistent" not in ev


def test_domain_dedup_keeps_opposite_strand():
    from teagle_core import domains
    # same code + same nt region but OPPOSITE strands -> both kept (different features)
    hits = [{"domain": "INT", "nt": [100, 300], "strand": "+", "score": 90},
            {"domain": "INT", "nt": [100, 300], "strand": "-", "score": 40}]
    kept = domains._dedup_domains(hits)
    assert len(kept) == 2 and {h["strand"] for h in kept} == {"+", "-"}


def test_domain_dedup_drops_same_strand_lower_score():
    from teagle_core import domains
    hits = [{"domain": "RT", "nt": [100, 300], "strand": "+", "score": 90},
            {"domain": "RT", "nt": [150, 350], "strand": "+", "score": 40}]   # same code+strand, overlap
    kept = domains._dedup_domains(hits)
    assert len(kept) == 1 and kept[0]["score"] == 90


def test_gypsy_when_integrase_downstream_plus_strand():
    cl = classify.classify(LTR, [_dom("RT", (100, 1000), "+", 90), _dom("INT", (1200, 1400), "+", 40)])
    assert cl["superfamily"].startswith("Gypsy"), cl["superfamily"]
    assert cl["confidence"] == "High"


def test_minus_strand_translation_order_respected():
    # minus-strand pol reads right->left: INT at LARGER nt is N-terminal -> Copia
    cl = classify.classify(LTR, [_dom("INT", (1200, 1400), "-", 40), _dom("RT", (100, 1000), "-", 90)])
    assert cl["superfamily"].startswith("Copia"), cl["superfamily"]
    # mirror: INT at SMALLER nt on minus strand is C-terminal -> Gypsy
    cl2 = classify.classify(LTR, [_dom("INT", (100, 300), "-", 40), _dom("RT", (500, 1400), "-", 90)])
    assert cl2["superfamily"].startswith("Gypsy"), cl2["superfamily"]


def test_order_indeterminate_downgrades_confidence():
    # INT and RT on different strands -> order not resolvable -> must not be High
    cl = classify.classify(LTR, [_dom("INT", (100, 300), "+", 40), _dom("RT", (500, 1400), "-", 90)])
    assert cl["confidence"] != "High"
    assert cl["superfamily"].split(" ")[0] in ("Copia", "Gypsy")


def test_split_orf_does_not_flip_on_orf_length():
    # historical bug: INT in a SHORTER ORF (higher length-rank) but genomically upstream of a long
    # RT ORF must stay Copia by translation position, not flip to Gypsy by ORF-length rank
    cl = classify.classify(LTR, [_dom("RT", (2000, 3500), "+", 95), _dom("INT", (200, 600), "+", 45)])
    assert cl["superfamily"].startswith("Copia"), cl["superfamily"]
    assert cl["confidence"] == "High"
