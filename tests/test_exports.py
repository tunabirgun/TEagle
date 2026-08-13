"""Image-export integrity: every figure must export to valid SVG and a correctly-sized PNG in every
background mode (dark / light / uv / mono / transparent), single-row and wrapped multi-row gels alike.
Runs headless via the same save_svg / render_png the export buttons call."""
import builtins, os, re, struct, sys, tempfile, zlib
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
_PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


def _png_chunk(chunk_type, data):
    crc = zlib.crc32(chunk_type + data) & 0xFFFFFFFF
    return struct.pack(">I", len(data)) + chunk_type + data + struct.pack(">I", crc)


def _png_records(payload):
    assert payload.startswith(_PNG_SIGNATURE)
    offset = len(_PNG_SIGNATURE)
    records = []
    while offset < len(payload):
        length = struct.unpack(">I", payload[offset:offset + 4])[0]
        end = offset + 12 + length
        records.append((payload[offset + 4:offset + 8], payload[offset + 8:offset + 8 + length], offset, end))
        offset = end
    assert offset == len(payload)
    return records


def _png_bytes(records):
    return _PNG_SIGNATURE + b"".join(_png_chunk(kind, data) for kind, data in records)


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
    assert widgets._csv_escape("=cmd", ",") == "'=cmd"
    assert widgets._csv_escape("a,b", ",") == '"a,b"'


def test_csv_and_xlsx_agree_on_negative_numbers():
    """The same value must not export as a number in one format and text in another.

    _xlsx_val has always kept "-5" numeric (see the test above), while _csv_escape quoted it to
    "'-5" — so every negative delta-G in the primer-QC tables reached a CSV as text, and the same
    table exported to XLSX and CSV disagreed about its own type. A CSV read into R or pandas came
    back with a character column for every thermodynamic field. Numbers are now exempt from the
    formula guard in both paths; a non-numeric signed cell is still neutralised."""
    for v in ("-5", "-9.4", "-26.61", "+1.5e-3", "0.0"):
        assert widgets._csv_escape(v, ",") == v, v
        assert widgets._xlsx_val(v) == float(v), v
    # the injection guard still bites where it should
    assert widgets._csv_escape("-1+cmd", ",") == "'-1+cmd"
    assert widgets._csv_escape("+cmd", ",") == "'+cmd"
    assert widgets._csv_escape("@ref", ",") == "'@ref"


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
    # capture headers and notes too: the context-menu path silently dropped the on-screen caveats that the
    # visible Export button carried, so which entry point a user happened to take decided whether their
    # exported table kept its scientific qualifications. Both paths must now pass the same things.
    monkeypatch.setattr(widgets, "export_table",
                        lambda h, r, b, p, fmt=None, notes=None: calls.append((fmt, h, notes)))
    monkeypatch.setattr(widgets, "pick_table_format", lambda parent, pos: "csv")
    t = widgets.DataTable(["c"], {}); t.set_rows([["v"]])
    t.export_notes = ["a caveat that must reach the file"]
    t._do_export(QPoint(0, 0))
    assert [c[0] for c in calls] == ["csv"]
    assert calls[0][1] == ["c"]                               # export headers, not the raw widget headers
    assert calls[0][2] == ["a caveat that must reach the file"]
    assert not hasattr(widgets, "add_export_submenu")         # the arrow-prone submenu helper is gone


def test_data_table_copy_cell_uses_export_value():
    t = widgets.DataTable(["display"])
    t.set_rows([["12.0 !"]], exports=[["12.0"]])
    t._copy_cell(0, 0)
    assert _app.clipboard().text() == "12.0"


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


# ---------------- self-similarity plot: custom colours + PDF export ----------------
def _dot_matrix():
    sys.path.insert(0, os.path.join(os.path.dirname(_NATIVE), "backend"))
    from teagle_core import dotplot, examples
    return dotplot.self_matrix(examples.load("M11240"), k=dotplot.suggest_k([]))


def test_dot_plot_export_honours_custom_mark_colours():
    """The colour the user picks must reach the exported figure (screen and export share one path); the
    Okabe–Ito defaults must not leak in once overridden. The dot plot draws both layers as flat fills."""
    m = _dot_matrix()
    svg = figures.svg_dotplot(m, W=900, theme="white", fwd="#7B1FA2", rev="#00897B")
    assert "#7B1FA2" in svg and "#00897B" in svg
    assert figures._DOT_FWD not in svg and figures._DOT_REV not in svg


def test_heat_map_export_honours_the_custom_layer_colour():
    """The heat map ramps the active layer's base colour toward white, so the base appears at peak
    intensity and the default base must be gone. Only the shown layer's colour is threaded."""
    m = _dot_matrix()
    svg = figures.svg_dotheat(m, W=900, theme="white", which="forward", fwd="#7B1FA2", rev="#00897B")
    assert figures._heat_ramp(1.0, "#7B1FA2") == "#7B1FA2"    # the ramp's ceiling is the base itself
    assert "#7B1FA2" in svg                                    # a peak-intensity cell renders at the base
    assert figures._DOT_FWD not in svg                        # the default blue must not survive the override


def test_dot_export_carries_its_own_caveat_only_on_export():
    """The honesty caveat is a separate on-screen widget, so the exported figure must embed its own one-line
    limit — but only on export, never in the interactive SVG (which shows the widget beneath it)."""
    m = _dot_matrix()
    assert "not local alignment" in figures.svg_dotplot(m, W=900, for_export=True)
    assert "not local alignment" not in figures.svg_dotplot(m, W=900, for_export=False)


def test_dot_pdf_export_is_a_valid_vector_pdf(tmp_path):
    m = _dot_matrix()
    svg = figures.svg_dotplot(m, W=920, for_export=True, fwd="#7B1FA2", rev="#00897B")
    p = str(tmp_path / "selfsim.pdf")
    widgets.render_pdf(svg, p)
    assert os.path.getsize(p) > 1000
    with open(p, "rb") as fh:
        assert fh.read(5) == b"%PDF-"
    from teagle_core.artifact_metadata import media_metadata_violations
    assert media_metadata_violations(p) == []


def test_runtime_figure_exports_have_author_only_metadata(tmp_path):
    svg = '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 100 50"><text x="5" y="25">x</text></svg>'
    outputs = {
        "svg": tmp_path / "figure.svg",
        "png": tmp_path / "figure.png",
        "pdf": tmp_path / "figure.pdf",
    }
    widgets.save_svg(svg, str(outputs["svg"]))
    widgets.render_png(svg, str(outputs["png"]), scale=1)
    widgets.render_pdf(svg, str(outputs["pdf"]))
    from teagle_core.artifact_metadata import media_metadata_violations
    for path in outputs.values():
        assert media_metadata_violations(path) == [], path


def test_png_metadata_cleanup_needs_no_pillow_and_rejects_extra_text(tmp_path, monkeypatch):
    from teagle_core import artifact_metadata

    path = tmp_path / "figure.png"
    svg = '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 100 50"><rect width="100" height="50"/></svg>'
    widgets.render_png(svg, str(path), scale=1)
    pristine = path.read_bytes()
    pristine_records = [(kind, data) for kind, data, _start, _end in _png_records(pristine)]
    before = QImage(str(path))
    before_pixels = [before.pixel(x, y) for y in range(before.height()) for x in range(before.width())]
    dirty_records = pristine_records[:-1] + [(b"tEXt", b"Software\0synthetic generator")] + pristine_records[-1:]
    path.write_bytes(_png_bytes(dirty_records))
    assert artifact_metadata.media_metadata_violations(path)

    real_import = builtins.__import__
    def reject_pillow(name, *args, **kwargs):
        if name == "PIL" or name.startswith("PIL."):
            raise AssertionError("PNG metadata handling must not import Pillow")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", reject_pillow)
    artifact_metadata.clean_png_metadata(path)
    assert artifact_metadata.media_metadata_violations(path) == []
    cleaned_records = [(kind, data) for kind, data, _start, _end in _png_records(path.read_bytes())]
    assert [item for item in cleaned_records if item[0] == b"tEXt"] == [(b"tEXt", b"Author\0Tuna Birgun")]
    assert [item for item in cleaned_records if item[0] != b"tEXt"] == [
        item for item in pristine_records if item[0] != b"tEXt"
    ]
    after = QImage(str(path))
    assert (after.width(), after.height()) == (before.width(), before.height())
    assert [after.pixel(x, y) for y in range(after.height()) for x in range(after.width())] == before_pixels


@pytest.mark.parametrize("chunk", [
    (b"vpAg", b"generator=synthetic;tracking=123"),
    (b"iTXt", b"Author\0\0\0en\0Translated Author\0Tuna Birgun"),
])
def test_png_metadata_cleanup_removes_private_and_extended_text_chunks(tmp_path, chunk):
    from teagle_core import artifact_metadata

    path = tmp_path / "figure.png"
    widgets.render_png('<svg xmlns="http://www.w3.org/2000/svg" width="8" height="8"/>', str(path), scale=1)
    records = [(kind, data) for kind, data, _start, _end in _png_records(path.read_bytes())]
    path.write_bytes(_png_bytes(records[:-1] + [chunk] + records[-1:]))
    assert artifact_metadata.media_metadata_violations(path)
    artifact_metadata.clean_png_metadata(path)
    assert artifact_metadata.media_metadata_violations(path) == []
    assert chunk not in [(kind, data) for kind, data, _start, _end in _png_records(path.read_bytes())]


def test_png_metadata_gate_rejects_bad_crc_without_replacing_source(tmp_path):
    from teagle_core import artifact_metadata

    path = tmp_path / "figure.png"
    widgets.render_png('<svg xmlns="http://www.w3.org/2000/svg" width="8" height="8"/>', str(path), scale=1)
    corrupt = bytearray(path.read_bytes())
    idat = next(record for record in _png_records(corrupt) if record[0] == b"IDAT")
    corrupt[idat[3] - 1] ^= 1
    path.write_bytes(corrupt)
    original = path.read_bytes()
    assert any("invalid CRC" in item for item in artifact_metadata.media_metadata_violations(path))
    with pytest.raises(ValueError, match="invalid CRC"):
        artifact_metadata.clean_png_metadata(path)
    assert path.read_bytes() == original


@pytest.mark.parametrize("mutation", ["duplicate-ihdr", "nonconsecutive-idat"])
def test_png_metadata_gate_rejects_invalid_chunk_structure(tmp_path, mutation):
    from teagle_core import artifact_metadata

    path = tmp_path / "figure.png"
    widgets.render_png('<svg xmlns="http://www.w3.org/2000/svg" width="8" height="8"/>', str(path), scale=1)
    records = [(kind, data) for kind, data, _start, _end in _png_records(path.read_bytes())]
    if mutation == "duplicate-ihdr":
        records.insert(1, records[0])
    else:
        index = next(index for index, item in enumerate(records) if item[0] == b"IDAT")
        kind, data = records[index]
        split = max(1, len(data) // 2)
        records[index:index + 1] = [
            (kind, data[:split]),
            (b"tEXt", b"Author\0Tuna Birgun"),
            (kind, data[split:]),
        ]
    path.write_bytes(_png_bytes(records))
    original = path.read_bytes()
    assert artifact_metadata.media_metadata_violations(path)
    with pytest.raises(ValueError):
        artifact_metadata.clean_png_metadata(path)
    assert path.read_bytes() == original


def test_pdf_export_failure_preserves_existing_file(tmp_path, monkeypatch):
    """Postprocessing failure must not replace an existing journal figure with a partial PDF."""
    from teagle_core import artifact_metadata

    path = tmp_path / "figure.pdf"
    original = b"existing validated PDF"
    path.write_bytes(original)
    monkeypatch.setattr(
        artifact_metadata,
        "clean_pdf_metadata",
        lambda _path: (_ for _ in ()).throw(OSError("synthetic metadata-cleaning failure")),
    )
    svg = '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 100 50"/>'
    with pytest.raises(OSError, match="synthetic metadata-cleaning failure"):
        widgets.render_pdf(svg, str(path))
    assert path.read_bytes() == original
    assert not (tmp_path / "figure.pdf.part").exists()
