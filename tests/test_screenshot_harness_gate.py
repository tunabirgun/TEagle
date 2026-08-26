import copy
import json
import os
import sys
import time
from pathlib import Path

import pytest
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
VERIFY = ROOT / "verification"
if not (VERIFY / "screenshot_gate.py").is_file():
    pytest.skip(
        "local screenshot-harness helper is gitignored and unavailable in this checkout",
        allow_module_level=True,
    )
if str(VERIFY) not in sys.path:
    sys.path.insert(0, str(VERIFY))

import screenshot_gate as gate
from teagle_core.artifact_metadata import clean_png_metadata


def _png(path: Path, color="white"):
    path.parent.mkdir(parents=True, exist_ok=True)
    Image.new("RGB", (2, 2), color).save(path)
    clean_png_metadata(path)


def _recapture_tree(tmp_path):
    started = time.time_ns()
    outputs = {}
    for rel in sorted(gate.RECAP_ALLOWED):
        path = tmp_path / rel
        _png(path)
        record = gate.output_record(tmp_path, path)
        outputs[rel] = {key: record[key] for key in ("sha256", "width", "height")}
    return {
        "_capture": {"schemaVersion": "1.0", "appVersion": gate.RELEASE_VERSION, "status": "ok",
                     "startedNs": started, "completedNs": time.time_ns(), "failures": [], "outputs": outputs},
        "shots": {}, "img": {}, "man": {},
    }


def test_release_version_font_and_exit_gates_fail_closed():
    assert gate.version_issues("3.9.0") and not gate.version_issues("4.0.0")
    assert gate.font_issues(["CascadiaMono.ttf"]) and not gate.font_issues([])
    assert gate.exit_status(["capture failed"]) == 1 and gate.exit_status([]) == 0


def test_capture_environment_pin_mutation_fails_closed():
    good = type("Capture", (), {"SCREEN_W": 1920, "SCREEN_H": 1080, "WIN_W": 1500,
                                 "WIN_H": 1150, "THEME": "light", "UI_SCALE": 1.0})
    assert gate.capture_environment_issues(good) == []
    bad = type("Capture", (), {"SCREEN_W": 1920, "SCREEN_H": 1080, "WIN_W": 1499,
                                "WIN_H": 1150, "THEME": "light", "UI_SCALE": 1.0})
    assert gate.capture_environment_issues(bad)


def test_recapture_sidecar_accepts_only_fresh_exact_content_bindings(tmp_path):
    sidecar = _recapture_tree(tmp_path)
    exact = {"docs/img/overview.png": (2, 2), "report/manual_figures/overview.png": (2, 2)}
    assert gate.validate_sidecar(tmp_path, sidecar, exact_dimensions=exact) == []

    cases = []
    stale_version = copy.deepcopy(sidecar); stale_version["_capture"]["appVersion"] = "3.8.1"; cases.append(stale_version)
    failed = copy.deepcopy(sidecar); failed["_capture"]["status"] = "failed"; cases.append(failed)
    missing = copy.deepcopy(sidecar); missing["_capture"]["outputs"].pop(next(iter(gate.RECAP_REQUIRED))); cases.append(missing)
    unknown = copy.deepcopy(sidecar); unknown["_capture"]["outputs"]["docs/img/private.png"] = {}; cases.append(unknown)
    mismatch = copy.deepcopy(sidecar); mismatch["_capture"]["outputs"][next(iter(gate.RECAP_REQUIRED))]["sha256"] = "0" * 64; cases.append(mismatch)
    for mutation in cases:
        assert gate.validate_sidecar(tmp_path, mutation, exact_dimensions=exact), mutation

    rel = next(iter(gate.RECAP_REQUIRED))
    old = sidecar["_capture"]["startedNs"] - 1_000_000
    os.utime(tmp_path / rel, ns=(old, old))
    assert any("STALE" in issue for issue in gate.validate_sidecar(tmp_path, sidecar, exact_dimensions=exact))


def test_recapture_sidecar_rejects_rebound_wrong_pinned_dimension(tmp_path):
    sidecar = _recapture_tree(tmp_path)
    rel = "docs/img/overview.png"
    Image.new("RGB", (3, 2), "white").save(tmp_path / rel)
    clean_png_metadata(tmp_path / rel)
    record = gate.output_record(tmp_path, tmp_path / rel)
    sidecar["_capture"]["outputs"][rel] = {
        key: record[key] for key in ("sha256", "width", "height")
    }
    issues = gate.validate_sidecar(
        tmp_path, sidecar,
        exact_dimensions={rel: (2, 2), "report/manual_figures/overview.png": (2, 2)})
    assert any("dimension" in issue for issue in issues)


def test_png_gate_mutations_cover_dimensions_freshness_and_decodability(tmp_path):
    rel = "docs/img/overview.png"
    path = tmp_path / rel
    _png(path)
    started = path.stat().st_mtime_ns
    assert gate.validate_png_set(tmp_path, [rel], started_ns=started, exact_dimensions={rel: (2, 2)}) == []
    assert any("dimension" in issue for issue in
               gate.validate_png_set(tmp_path, [rel], exact_dimensions={rel: (3, 2)}))
    assert any("STALE" in issue for issue in gate.validate_png_set(tmp_path, [rel], started_ns=started + 1))
    path.write_bytes(b"not a png")
    assert any("undecodable" in issue for issue in gate.validate_png_set(tmp_path, [rel]))


def test_duplicate_hash_gate_detects_mutation(tmp_path):
    left, right = gate.DUPLICATE_PAIRS[0]
    _png(tmp_path / left)
    _png(tmp_path / right)
    assert gate.validate_duplicate_pairs(tmp_path, [(left, right)]) == []
    _png(tmp_path / right, "black")
    assert any("hash mismatch" in issue for issue in gate.validate_duplicate_pairs(tmp_path, [(left, right)]))


def test_inventory_and_uncited_gates_detect_mutations(tmp_path):
    shots = tmp_path / "docs/img/screenshots"; manual = tmp_path / "report/manual_figures"; img = tmp_path / "docs/img"
    shots.mkdir(parents=True); manual.mkdir(parents=True); img.mkdir(exist_ok=True)
    for name in gate.RUN_SCREENSHOTS | gate.RECAP_SCREENSHOTS:
        _png(shots / name)
    for name in gate.MANUAL_FIGURES:
        _png(manual / name)
    for name in gate.README_FIGURES:
        _png(img / name)
    for name in gate.STATIC_ASSETS:
        path = img / name; path.write_bytes(b"asset")
    (tmp_path / "report/teagle_manual.tex").write_text(
        "\n".join(f"\\includegraphics{{{name}}}" for name in gate.MANUAL_FIGURES - gate.CROP_SOURCES),
        encoding="utf-8")
    assert gate.validate_inventory(tmp_path) == []
    assert gate.validate_tex_citations(tmp_path) == []
    _png(shots / "old_panel.png")
    assert any("UNDOC" in issue for issue in gate.validate_inventory(tmp_path))
    (tmp_path / "report/teagle_manual.tex").write_text("", encoding="utf-8")
    assert any("UNCITED" in issue for issue in gate.validate_tex_citations(tmp_path))


def test_current_tex_cites_every_expected_manual_figure():
    if not (ROOT / "report/teagle_manual.tex").exists():
        pytest.skip("report/teagle_manual.tex not present in this checkout")
    assert gate.validate_tex_citations(ROOT) == []


@pytest.mark.parametrize("marker", ["FAIL", "MISS", "UNDOC", "STALE", "UNCITED"])
def test_manifest_gate_rejects_every_failure_marker(marker):
    required = ("CSV", "PDF")
    good = "TEagle documentation screenshots — TEagle v4.0.0\nOK    export_CSV\nOK    export_PDF\n"
    assert gate.validate_manifest_text(good, required_exports=required) == []
    assert gate.validate_manifest_text(good + marker + " mutation\n", required_exports=required)


def test_manifest_gate_rejects_old_version_and_missing_export():
    text = "TEagle documentation screenshots — TEagle v3.8.1\nOK    export_CSV\n"
    issues = gate.validate_manifest_text(text, required_exports=("CSV", "PDF"))
    assert any("version" in issue for issue in issues)
    assert any("export_PDF" in issue for issue in issues)


def test_manifest_gate_rejects_false_off_target_wording_for_unsealed_ltr_products():
    prefix = "TEagle documentation screenshots — TEagle v4.0.0\n"
    good = (prefix + "OK    pcr_off_target [runs] In-silico PCR card: repeated-LTR priming sites; "
            "three target-unassessed products because no intended locus was sealed\n")
    bad = prefix + "OK    pcr_off_target [runs] In-silico PCR card: off-target bands\n"
    assert gate.validate_manifest_text(good) == []
    assert any("overstates unsealed products" in issue for issue in gate.validate_manifest_text(bad))


def test_capture_scripts_use_fail_closed_gate_and_exclude_legacy_harness():
    recap = (VERIFY / "recapture_shots.py").read_text(encoding="utf-8")
    runs = (VERIFY / "screenshot_runs.py").read_text(encoding="utf-8")
    joined = recap + runs
    assert "font_issues" in joined
    assert "validate_sidecar" in runs and "validate_manifest_text" in runs
    assert "sha256=" in runs and "_image_evidence" in runs
    assert "NamedTemporaryFile" in joined and "os.replace" in joined
    assert "raise SystemExit" in recap and "raise SystemExit" in runs
    assert "capture_screens.py" not in joined
    assert "capture_screens" not in json.dumps({"recap": recap, "runs": runs})


def test_current_export_probe_and_family_partition_are_complete():
    runs = (VERIFY / "screenshot_runs.py").read_text(encoding="utf-8")
    required = {
        "CSV", "TSV", "XLSX", "SVG", "PNG", "PDF", "FASTA_SEQUENCE",
        "GFF3_ANNOTATION", "BED_ANNOTATION", "GFF3_GENE_MODEL", "BED_GENE_MODEL",
        "FASTA_GENE_MODEL", "GFF3_SPLICE", "BED_SPLICE", "FASTA_SPLICE",
        "RESULT_JSON", "MANIFEST_JSON",
    }
    for name in required:
        assert f'"{name}"' in runs
    assert "TemporaryDirectory" in runs
    assert "strict_json" in runs
    assert '"target_span"' not in runs
    assert '"intended_product_span"' in runs and '"design_provenance"' in runs
    assert "_analysis_requests[generation]" in runs
    assert "all categories displayed and exported" in runs
    assert "TE-classified" in runs and "non-TE repeat" in runs and "unclassified" in runs
