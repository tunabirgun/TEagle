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


def test_export_submenu_actions_route_to_each_format(monkeypatch):
    from PySide6.QtWidgets import QMenu, QWidget
    calls = []
    monkeypatch.setattr(widgets, "export_table", lambda h, r, b, p, fmt=None: calls.append(fmt))
    host = QWidget(); m = QMenu(host)
    sub = widgets.add_export_submenu(m, ["c"], lambda: [["v"]], "base", host)
    import gc; gc.collect()                                    # the parented submenu must outlive a GC pass
    for a in list(sub.actions()):
        a.trigger()
    assert calls == [f for _, f in widgets._available_formats()]


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
