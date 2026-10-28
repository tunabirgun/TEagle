import copy

import pytest
import report


@pytest.mark.parametrize("run_type", ("primer", "in-silico-pcr"))
@pytest.mark.parametrize("region, expected", (
    ([125, 75], "125–200 (0-based, half-open)"),
    ([0, 80], "0–80 (0-based, half-open)"),
))
def test_included_region_displays_start_plus_length(run_type, region, expected):
    manifest = {"runType": run_type, "parameters": {"included_region": region}}
    original = copy.deepcopy(manifest)
    rows = [row for _, section in report.param_sections(manifest) for row in section]
    assert rows[0]["text"] == expected
    assert manifest == original


def test_endpoint_span_and_target_keep_their_own_contracts():
    assert report.param_row("intended_product_span", [125, 200])[1]["text"] == "125–200 (0-based, half-open)"
    assert report.param_row("primer3_target", [125, 75])[1]["text"] == "start 125, length 75 nt (0-based)"
    assert report.param_row("included_region", None)[1]["text"] == "the whole record"
