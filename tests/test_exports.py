"""Image-export integrity: every figure must export to valid SVG and a correctly-sized PNG in every
background mode (dark / light / uv / mono / transparent), single-row and wrapped multi-row gels alike.
Runs headless via the same save_svg / render_png the export buttons call."""
import os, sys, re, tempfile
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
import pytest

pytest.importorskip("PySide6")
_NATIVE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "app", "native")
if _NATIVE not in sys.path:
    sys.path.insert(0, _NATIVE)

from PySide6.QtWidgets import QApplication
from PySide6.QtGui import QImage
from PySide6.QtSvg import QSvgRenderer
from PySide6.QtCore import QByteArray
import figures, widgets

_app = QApplication.instance() or QApplication([])


def _lanes(n):
    return [{"label": f"P{i+1}", "amplicons": [{"length": 150 + i * 11, "on_target": i % 3 == 0,
             "single_primer": i % 4 == 0, "fwd_mm": i % 3, "rev_mm": 0, "source": "src"}]} for i in range(n)]


def _check_export(svg, tmp, name, transparent):
    assert QSvgRenderer(QByteArray(svg.encode())).isValid(), f"{name}: invalid SVG"
    p_svg = os.path.join(tmp, name + ".svg")
    widgets.save_svg(svg, p_svg)
    assert os.path.getsize(p_svg) > 200
    p_png = os.path.join(tmp, name + ".png")
    widgets.render_png(svg, p_png, scale=3)
    img = QImage(p_png)
    assert not img.isNull(), f"{name}: null PNG"
    vb = re.search(r'viewBox="0 0 ([\d.]+) ([\d.]+)"', svg)
    assert img.width() == int(float(vb.group(1)) * 3) and img.height() == int(float(vb.group(2)) * 3)
    assert any(img.pixelColor(x, y).alpha() > 0 for x in range(0, img.width(), 40)
               for y in range(0, img.height(), 40)), f"{name}: blank PNG"
    corner = img.pixelColor(2, 2).alpha()
    assert (corner < 255) if transparent else (corner == 255), f"{name}: wrong background alpha ({corner})"


@pytest.mark.parametrize("mode", ["dark", "white", "uv", "mono", "transparent"])
@pytest.mark.parametrize("n", [2, 23])          # single-row and wrapped multi-row
def test_gel_export_every_mode(mode, n, tmp_path):
    svg = figures.svg_gel({"lanes": _lanes(n)}, mode)
    _check_export(svg, str(tmp_path), f"gel_{mode}_{n}", transparent=(mode == "transparent"))


def test_xlsx_val_typing_and_injection_guard():
    assert widgets._xlsx_val("5") == 5 and isinstance(widgets._xlsx_val("5"), int)
    assert widgets._xlsx_val("1.5") == 1.5 and isinstance(widgets._xlsx_val("1.5"), float)
    assert widgets._xlsx_val("-5") == -5                      # negative number stays numeric
    assert widgets._xlsx_val("-") == "-" and widgets._xlsx_val("+") == "+"   # bare strand cell literal
    assert widgets._xlsx_val("=SUM(A1)") == "'=SUM(A1)"        # formula neutralised
    assert widgets._xlsx_val("120–450") == "120–450"          # en-dash coords stay text


def test_csv_escape_keeps_bare_sign():
    assert widgets._csv_escape("-", ",") == "-" and widgets._csv_escape("+", ",") == "+"
    assert widgets._csv_escape("-5", ",") == "'-5"            # multi-char leading sign guarded
    assert widgets._csv_escape("=cmd", ",") == "'=cmd"
    assert widgets._csv_escape("a,b", ",") == '"a,b"'


def test_export_xlsx_strand_column_not_corrupted(tmp_path):
    load_workbook = pytest.importorskip("openpyxl").load_workbook
    p = str(tmp_path / "t.xlsx")
    widgets._export_xlsx(["Dfam family", "Str", "Score"], [["Copia_I", "-", 40720], ["Copia_LTR", "+", 2405]], p)
    ws = load_workbook(p).active
    assert ws["A1"].font.bold and ws.freeze_panes == "A2"
    rows = list(ws.iter_rows(min_row=2, values_only=True))
    assert rows[0][1] == "-" and rows[1][1] == "+"           # strand not turned into '-/'+
    assert rows[0][2] == 40720 and isinstance(rows[0][2], int)


def test_gene_flanks_and_gaps_regions():
    # 5' + 3' flanks surface; an intron-covered middle is NOT a gap
    fg = figures._flanks_and_gaps([{"start": 100, "end": 200}, {"start": 400, "end": 500}],
                                  [{"start": 200, "end": 400}], 700)
    got = {(r["label"], r["start"], r["end"]) for r in fg}
    assert got == {("5prime_flank", 0, 100), ("3prime_flank", 500, 700)}
    # an uncovered interior hole IS a gap
    fg2 = figures._flanks_and_gaps([{"start": 100, "end": 200}, {"start": 400, "end": 500}], [], 500)
    assert ("gap", 200, 400) in {(r["label"], r["start"], r["end"]) for r in fg2}


def test_gene_model_track_includes_flanks_when_requested():
    m0 = figures.gv_tracks_from_gene({"exons": [{"start": 100, "end": 200}], "introns": [], "cds": []}, 500)
    m1 = figures.gv_tracks_from_gene({"exons": [{"start": 100, "end": 200}], "introns": [], "cds": []}, 500,
                                     include_flanks=True)
    labels0 = [f.get("label") for t in m0["tracks"] for f in t["features"]]
    labels1 = [f.get("label") for t in m1["tracks"] for f in t["features"]]
    assert "5′ flank" not in labels0                          # off by default (splice viewer must stay clean)
    assert "5′ flank" in labels1 and "3′ flank" in labels1    # readable on-glyph names (FASTA id sanitised in _feat_menu)


def test_derived_exon_marked_distinctly():
    # a CDS-inferred exon must be visually distinct from an annotated one (honesty: never cite an
    # inferred coord as a GenBank annotation)
    gm = {"exons": [{"start": 0, "end": 100}, {"start": 300, "end": 400, "derived": True}],
          "introns": [], "cds": []}
    feats = figures.gv_tracks_from_gene(gm, 500)["tracks"][0]["features"]
    ann = next(f for f in feats if f["label"] == "exon")
    der = next(f for f in feats if f["label"] == "exon*")
    assert ann["color"] == figures.GENECOL["exon"]
    assert der["color"] == figures.GENECOL["exon_derived"] != figures.GENECOL["exon"]
    assert "derived from the record's CDS/mRNA" in der["tip"]


def test_gap_and_flank_use_distinct_colors():
    # an interior gap and a terminal flank must not render as one ambiguous colour
    gm = {"exons": [{"start": 100, "end": 150}, {"start": 300, "end": 350}], "introns": [], "cds": []}
    feats = figures.gv_tracks_from_gene(gm, 500, include_flanks=True)["tracks"][0]["features"]
    gap = next(f for f in feats if f["label"] == "gap")       # hole between the two exons
    flank = next(f for f in feats if f["label"] == "5′ flank")
    assert gap["color"] == figures.GENECOL["gap"]
    assert flank["color"] == figures.GENECOL["flank"]
    assert gap["color"] != flank["color"]


def test_available_formats_lists_all_when_xlsx_present():
    fmts = [f for _, f in widgets._available_formats()]
    if widgets._HAS_XLSX:
        assert fmts == ["xlsx", "csv", "tsv"]
    else:
        assert fmts == ["csv", "tsv"]


@pytest.mark.parametrize("fmt,ext,sep", [("csv", ".csv", ","), ("tsv", ".tsv", "\t")])
def test_export_table_fmt_writes_chosen_delimiter(fmt, ext, sep, tmp_path, monkeypatch):
    out = str(tmp_path / ("t" + ext))
    monkeypatch.setattr(widgets.QFileDialog, "getSaveFileName", staticmethod(lambda *a, **k: (out, "")))
    widgets.export_table(["type", "start"], [["LTR", "1,000"]], "base", None, fmt=fmt)
    text = open(out, encoding="utf-8-sig").read()
    assert text.splitlines()[0] == sep.join(["type", "start"])
    if fmt == "csv":
        assert '"1,000"' in text                              # comma cell quoted only in CSV
    else:
        assert "1,000" in text and '"1,000"' not in text      # TSV needs no quoting


def test_export_table_fmt_appends_missing_extension(tmp_path, monkeypatch):
    noext = str(tmp_path / "table")
    monkeypatch.setattr(widgets.QFileDialog, "getSaveFileName", staticmethod(lambda *a, **k: (noext, "")))
    widgets.export_table(["c"], [["v"]], "base", None, fmt="tsv")
    assert os.path.isfile(noext + ".tsv")


def test_export_flat_action_routes_chosen_format(monkeypatch):
    # the DataTable Export is now a single FLAT action (no submenu -> no native ▶ arrow to overlap the label);
    # it pops the flat format picker then routes the chosen format to export_table
    from PySide6.QtCore import QPoint
    calls = []
    monkeypatch.setattr(widgets, "export_table", lambda h, r, b, p, fmt=None: calls.append(fmt))
    monkeypatch.setattr(widgets, "pick_table_format", lambda parent, pos: "csv")
    t = widgets.DataTable(["c"], {}); t.set_rows([["v"]])
    t._do_export(QPoint(0, 0))
    assert calls == ["csv"]
    assert not hasattr(widgets, "add_export_submenu")         # the arrow-prone submenu helper is gone


@pytest.mark.parametrize("theme", ["dark", "white"])
def test_genome_export(theme, tmp_path):
    sys.path.insert(0, os.path.join(os.path.dirname(_NATIVE), "backend"))
    import engine
    seq = open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                            "tests", "fixtures", "M11240.fasta")).read()
    rec = engine.run_analyze({"sequence": seq})["records"][0]
    model = figures.gv_tracks_from_rec(rec)
    svg = figures.svg_genome(model, {"start": 0.0, "end": float(model["length"])}, 920, theme, False)
    _check_export(svg, str(tmp_path), f"genome_{theme}", transparent=False)


def _erv_model():
    """Model carrying the longest track name in the app: the ERV 'env mRNA (predicted)' architecture track."""
    rec = {"composition": {"length": 9472},
           "structural": [{"type": "LTR (terminal direct repeat)", "five_prime": [0, 968],
                           "three_prime": [8504, 9472]},
                          {"type": "PPT (polypurine tract)", "pos": [8486, 8504], "purine_frac": 0.9},
                          {"type": "poly-A tail", "pos": [9400, 9430]}],
           "domains": [{"domain": "ENV", "label": "envelope", "nt": [6500, 7800], "score": 90.0}],
           "retroviral": {"leader_exon": [968, 1141], "intron": [1141, 6400], "env_exon": [6400, 8000]},
           "orfs": [{"start": 1111, "end": 3112, "strand": "+", "frame": 2, "length_aa": 666}]}
    return figures.gv_tracks_from_rec(rec)


def test_track_label_gutter_fits_longest_name():
    """The label gutter must hold the widest track name — 'env mRNA (predicted)' was clipped to
    'nv mRNA (predicted)' at the fixed 96 px gutter. Mono metric = upper bound on the drawn width."""
    model = _erv_model()
    names = [t["name"] for t in model["tracks"]]
    assert "env mRNA (predicted)" in names
    ML = figures.gv_gutter(model)
    for n in names:
        start_x = ML - figures.GV_NAME_GAP - len(n) * figures.GV_NAME_PX * figures.GV_CHAR_EM
        assert start_x >= 0, f"{n!r} starts at x={start_x:.1f}, clipped by the figure's left edge"
    assert figures.gv_gutter({"tracks": [{"name": "CDS"}]}) == figures.GV_ML_MIN   # short names keep the old layout


@pytest.mark.parametrize("W", [320, 620, 920])
def test_genome_labels_clear_left_edge_when_rasterised(W):
    """Rasterise and require ink-free columns at the figure's left edge, so no name is cut off."""
    model = _erv_model()
    svg = figures.svg_genome(model, {"start": 0.0, "end": float(model["length"])}, W, "white")
    w, h = widgets._svg_size(svg)
    sc = 3
    img = QImage(int(w * sc), int(h * sc), QImage.Format_ARGB32)
    img.fill(0xFFFFFFFF)
    from PySide6.QtGui import QPainter
    p = QPainter(img); QSvgRenderer(QByteArray(svg.encode())).render(p); p.end()
    for x in range(0, sc):                                    # first logical pixel column must be blank
        assert all(img.pixelColor(x, y).lightness() > 230 for y in range(img.height())), \
            f"ink at the left edge (W={W}) — a track name is clipped"


def test_canvas_gutter_matches_the_svg():
    """The canvas' bp<->pixel map reads the same gutter the SVG drew, or zoom/pan anchors drift."""
    model = _erv_model()
    panel = widgets.GenomePanel(figures.svg_genome, "t")
    panel.set_model(model)
    c = panel.canvas
    svg = figures.svg_genome(model, {"start": 0.0, "end": float(model["length"])}, c._w, "white")
    ml = float(re.search(r'<rect x="([\d.]+)" y="[\d.]+" width="[\d.]+" height="13"', svg).group(1))
    assert abs(c.ML - ml) < 0.01                              # overview rect starts at the gutter
    panel.deleteLater()
