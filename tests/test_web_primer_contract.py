import json
from pathlib import Path
import shutil
import subprocess

import pytest


APP = Path(__file__).parents[1] / "app" / "web" / "app.js"


def _source():
    return APP.read_text(encoding="utf-8")


def _node(script):
    node = shutil.which("node")
    if not node:
        pytest.skip("Node.js is unavailable")
    result = subprocess.run([node, "-e", script], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    return json.loads(result.stdout)


def test_web_primer_numeric_parser_rejects_invalid_values():
    src = _source()
    helper = src[src.index("function requiredNumber"):src.index("function readPrimerParams")]
    script = f'''const values = {{}};
const $ = id => ({{value: values[id.slice(1)]}});
{helper}
const cases = [
  ["", {{integer:true,min:0}}], ["x", {{integer:true,min:0}}],
  ["2.9", {{integer:true,min:0}}], ["-1", {{integer:true,min:0}}]
];
const rejected = cases.map(([value, options]) => {{
  values.x = value; try {{ requiredNumber("x", "x", options); return false; }} catch(e) {{ return true; }}
}});
values.x = "18"; const accepted = requiredNumber("x", "x", {{integer:true,min:18}});
console.log(JSON.stringify({{rejected, accepted}}));'''
    observed = _node(script)
    assert observed == {"rejected": [True, True, True, True], "accepted": 18}


def test_web_primer_ranges_and_minimum_are_fail_closed():
    src = _source()
    block = src[src.index("function requiredNumber"):src.index('if($("#pPreset"))')]
    defaults = {"pMin": "150", "pMax": "500", "pTm": "60", "pMinS": "18",
                "pMaxS": "27", "pNum": "5", "pOptS": "20", "pTmMin": "57",
                "pTmMax": "63", "pGcMin": "40", "pGcMax": "60", "pPolyX": "4",
                "pGcClamp": "0"}
    script = f'''const values = {json.dumps(defaults)};
const $ = id => ({{value: values[id.slice(1)]}});
{block}
function rejected(changes) {{ Object.assign(values, changes); try {{ readPrimerParams(); return false; }} catch(e) {{ return true; }} }}
const observations = [rejected({{pMinS:"17"}}), rejected({{pMinS:"18",pMin:"501",pMax:"500"}}),
  rejected({{pMin:"150",pMax:"500",pTmMin:"64",pTm:"60",pTmMax:"63"}}),
  rejected({{pTmMin:"57",pTm:"60",pTmMax:"63",pGcMin:"61",pGcMax:"60"}})];
console.log(JSON.stringify(observations));'''
    assert _node(script) == [True, True, True, True]
    assert "pMinS:18" in src and "pMinS:17" not in src


def test_web_primer_geometry_provenance_and_wording_contract():
    src = _source()
    pcr = src[src.index("/* ---------- in-silico PCR engine"):src.index("/* ---------- provenance ----------")]
    assert "primer3_target: [start, end - start]" in src
    assert "included_region: [0, templateLength]" in src
    assert "target_span:" not in pcr
    for field in ("design_provenance:c.design_provenance", "design_sequence:c._design_sequence",
                  "primer3_target:c.primer3_target", "included_region:c.included_region",
                  "intended_product_span:c.intended_product_span",
                  "intended_endpoint_tolerance_bp:c.intended_endpoint_tolerance_bp"):
        assert field in pcr
    assert "manifestSha256" in pcr and "intended_product_span" in pcr
    assert "partial result; product counts are lower bounds" in pcr
    assert "absence of further products was not assessed" in pcr
    assert "target not covered" in pcr and "target coverage unavailable" in pcr
    assert "No additional amplicon" not in pcr
    assert "off-target</span>" not in pcr
    assert "product_length ?? a.length" in src
    assert "product_key" in src
    assert "const shortestPrimer = Math.min" in pcr
    assert "p.tp > shortestPrimer || p.max_mm > shortestPrimer" in pcr


def test_web_javascript_parses():
    node = shutil.which("node")
    if not node:
        pytest.skip("Node.js is unavailable")
    result = subprocess.run([node, "--check", str(APP)], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
