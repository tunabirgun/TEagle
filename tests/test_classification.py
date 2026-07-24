"""Golden classification tests: real fixtures pinned to their published TE class/superfamily.
Exercises the domain (pyhmmer) + classify pipeline end-to-end."""
from teagle_core import structural, domains, classify
from helpers import fixture_seq


def _classify(acc):
    seq = fixture_seq(acc)
    st = structural.detect_all(seq)
    dm = domains.scan_domains(seq)
    return classify.classify(st, dm), {d["domain"] for d in dm}


def test_copia_is_copia_superfamily():
    cl, dcodes = _classify("M11240")
    assert cl["te_class"] == "LTR/Copia"
    assert cl["superfamily"].startswith("Copia")
    assert {"RT", "INT"} <= dcodes                # RT + integrase detected
    assert cl["confidence"] == "High"


def test_gypsy_is_gypsy_superfamily():
    cl, dcodes = _classify("M12927")
    assert cl["te_class"] == "LTR/Gypsy"
    assert cl["superfamily"].startswith("Gypsy")
    assert {"RT", "INT"} <= dcodes
    assert cl["confidence"] == "High"


def test_l1_is_line():
    cl, dcodes = _classify("M80343")
    assert cl["te_class"] == "LINE"
    assert "RT" in dcodes
    assert "INT" not in dcodes                    # LINEs lack integrase


def test_tc1_is_dna_transposon():
    cl, dcodes = _classify("X01005")
    assert cl["te_class"].startswith("DNA/")
    assert "TPase" in dcodes
    assert "Class II" in cl["class"]


def test_classify_flags_composite_rt_plus_tpase():
    # a transposase co-occurring with RT (nested/composite locus) must be surfaced + confidence capped —
    # not silently dropped by the `elif tpase` branch that is unreachable once rt is present
    st = [{"type": "LTR (5')"}, {"type": "LTR (3')"}]
    dm = [{"domain": "RT", "nt": [1000, 1500], "strand": "+", "score": 100.0},
          {"domain": "INT", "nt": [500, 900], "strand": "+", "score": 90.0},     # INT before RT -> Copia, order resolvable
          {"domain": "TPase", "nt": [2000, 2400], "strand": "+", "score": 80.0, "class": "hAT"}]
    cl = classify.classify(st, dm)
    assert cl["superfamily"].startswith("Copia")                                  # primary call unchanged
    assert any("transposase domain also present" in e for e in cl["evidence"])    # composite surfaced
    assert cl["confidence"] != "High"                                             # capped from High by the composite signal
    assert "Class I" in cl["class"]                                               # still a Class I retrotransposon primary call


def test_ac_is_hat():
    cl, dcodes = _classify("X05424")
    assert cl["te_class"] == "DNA/hAT"
    assert "TPase" in dcodes


def test_domain_fields_present():
    dm = domains.scan_domains(fixture_seq("M11240"))
    assert dm, "copia should have detectable domains"
    d = dm[0]
    for k in ("domain", "label", "pfam", "score", "evalue", "aa", "nt", "dna", "protein"):
        assert k in d, f"domain missing field {k}"
    assert d["pfam"].startswith("PF")
    assert set(d["dna"]) <= set("ACGTN")
    assert len(d["protein"]) >= 20
    assert (d["nt"][1] - d["nt"][0]) == len(d["dna"])


def test_copia_gypsy_distinguished_by_integrase_order():
    copia, _ = _classify("M11240")
    gypsy, _ = _classify("M12927")
    assert copia["superfamily"].split(" ")[0] != gypsy["superfamily"].split(" ")[0]


# ---------------- ERV detection + structural completeness (v2.9.0; pure unit, no fetch) ----------------
def _mock_ltr(span=(0, 9000), ident=98.0):
    return [{"type": "LTR (terminal direct repeat)", "element_span": list(span), "identity": ident}]


def _dom(code, s, e, score=100.0, hmm=None):
    # hmm defaults to a CORE gag profile for GAG (capsid) so the mock reads as a real capsid hit;
    # pass hmm="zf-CCHC_5" to model a nucleocapsid-only gag hit.
    return {"domain": code, "label": code, "class": "retro", "score": score, "strand": "+", "nt": [s, e],
            "hmm": hmm if hmm is not None else ("Gag_p24" if code == "GAG" else code)}


def test_erv_full_architecture_is_intact():
    # env + paired LTRs + full gag-pol -> ERV, 'intact / autonomous-consistent', GAG..ENV architecture
    doms = [_dom("GAG", 2000, 2400), _dom("PR", 3600, 3760), _dom("RT", 4100, 4600),
            _dom("RNaseH", 5200, 5600), _dom("INT", 5800, 6100), _dom("ENV", 7900, 8500)]
    cl = classify.classify(_mock_ltr(), doms)
    assert cl["is_erv"] is True
    assert cl["completeness"]["tier"].startswith("intact")
    assert not cl["completeness"]["missing"]
    assert "GAG" in (cl["order"] or "") and "ENV" in (cl["order"] or "")   # env closes the architecture
    assert cl["completeness"]["scope"]                                     # honesty-scope note present


def test_nucleocapsid_only_gag_is_not_intact():
    # a promiscuous zf-CCHC nucleocapsid hit alone is not core gag evidence: still shown as a detected
    # module, but must not earn the intact / near-complete tier without a capsid/matrix hit
    doms = [_dom("GAG", 2000, 2400, hmm="zf-CCHC_5"), _dom("PR", 3600, 3760), _dom("RT", 4100, 4600),
            _dom("RNaseH", 5200, 5600), _dom("INT", 5800, 6100), _dom("ENV", 7900, 8500)]
    comp = classify.classify(_mock_ltr(), doms)["completeness"]
    assert not comp["tier"].startswith("intact")
    assert not comp["tier"].startswith("near-complete")
    assert "GAG" in comp["present"]                                        # detected, honestly shown


def test_env_detected_without_ltr_is_kept_and_kind_is_honest():
    # RT+INT+ENV but no detected LTR: ENV must appear in present (not silently dropped), and the
    # completeness kind must not assert an LTR that was never detected
    doms = [_dom("RT", 4100, 4600), _dom("INT", 5800, 6100), _dom("ENV", 7900, 8500)]
    comp = classify.classify([], doms)["completeness"]
    assert "ENV" in comp["present"]
    assert not comp["kind"].startswith("LTR")                              # no false LTR claim


def test_ltr_without_env_is_not_erv():
    doms = [_dom("GAG", 2000, 2400), _dom("RT", 4100, 4600), _dom("INT", 5800, 6100)]
    cl = classify.classify(_mock_ltr(), doms)
    assert cl["is_erv"] is False
    assert "ENV" not in (cl["order"] or "")


def test_partial_tier_never_claims_present_env_or_gag_missing():
    # HERV-H-like: full gag-pol-env recovered but the LTR heuristic missed the terminal repeats.
    # The tier must NOT say gag/env are incomplete while they sit in present (Loop-2 fix).
    doms = [_dom("GAG", 2000, 2400), _dom("PR", 3600, 3760), _dom("RT", 4100, 4600),
            _dom("RNaseH", 5200, 5600), _dom("INT", 5800, 6100), _dom("ENV", 7900, 8500)]
    comp = classify.classify([], doms)["completeness"]              # no structural LTR
    assert "ENV" in comp["present"] and "GAG" in comp["present"]
    assert "incomplete" not in comp["tier"]                         # nothing present is called incomplete


def test_plain_ltr_retrotransposon_tier_does_not_blame_env():
    # a canonical env-less Gypsy/Copia (RT+INT+LTR, no gag core, no env) must not be labelled env-decayed
    doms = [_dom("RT", 4100, 4600), _dom("INT", 5800, 6100)]
    comp = classify.classify(_mock_ltr(), doms)["completeness"]
    assert "env" not in comp["tier"].lower()
    assert "ENV" not in (comp["missing"] or [])


def test_partial_tier_does_not_embed_missing_list_the_banner_renders():
    # the banner appends 'not detected: <missing>' itself; the tier must not embed the same list (double-render)
    doms = [_dom("RT", 4100, 4600), _dom("INT", 5800, 6100)]        # pol-only + LTR -> partial, missing gag/PR/RNaseH
    comp = classify.classify(_mock_ltr(), doms)["completeness"]
    assert comp["tier"].startswith("partial")
    assert "not detected" not in comp["tier"]                       # list lives only in comp['missing']
    assert "GAG" in comp["missing"]


def test_ltr_with_env_only_domain_is_not_called_no_coding_detected():
    # an ERV relic that kept env but lost RT/pol: env is detected, so neither the superfamily, the tier,
    # nor the present list may claim 'no coding domain detected' (would contradict the domain table)
    doms = [_dom("ENV", 7900, 8500)]
    cl = classify.classify(_mock_ltr(), doms)
    assert "no coding domains detected" not in cl["superfamily"]
    comp = cl["completeness"]
    assert "ENV" in comp["present"]
    assert "no coding domain" not in comp["tier"]


def test_missing_gag_is_partial_not_intact():
    # pol core but no gag detected -> 'partial', gag listed as not-detected (scoped, honest)
    doms = [_dom("RT", 4100, 4600), _dom("INT", 5800, 6100)]
    cl = classify.classify(_mock_ltr(), doms)
    assert cl["completeness"]["tier"].startswith("partial")
    assert "GAG" in cl["completeness"]["missing"]
