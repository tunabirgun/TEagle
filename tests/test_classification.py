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
