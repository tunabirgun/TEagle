"""4.1 UI lane V1: theme tokens, figure builders and reusable widgets (behavioural, offscreen).

Run in its own pytest session: benchmarks/figures.py shadows native/figures.py when a benchmark test has
already put benchmarks/ on sys.path."""
import os, re, sys, math, types

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
import pytest

pytest.importorskip("PySide6")
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for _p in (os.path.join(_ROOT, "app", "native"), os.path.join(_ROOT, "app", "backend")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from PySide6.QtCore import Qt, QPoint, QPointF, QByteArray, QEvent
from PySide6.QtGui import QImage, QPainter, QColor, QMouseEvent
from PySide6.QtSvg import QSvgRenderer
from PySide6.QtTest import QTest
from PySide6.QtWidgets import (QApplication, QComboBox, QFrame, QLabel, QLineEdit, QPlainTextEdit, QPushButton,
                               QTextEdit, QVBoxLayout, QWidget)

import fonts
import theme
import figures
import widgets

_app = QApplication.instance() or QApplication([])
fonts.load_fonts()
FIXTURES = os.path.join(_ROOT, "tests", "fixtures")


# ---------------------------------------------------------------- colour helpers
def _rgba(s):
    m = re.match(r"rgba\((\d+),\s*(\d+),\s*(\d+),\s*([\d.]+)\)", s)
    if m:
        return tuple(int(v) for v in m.groups()[:3]), float(m.group(4))
    h = s.lstrip("#")
    return tuple(int(h[i:i + 2], 16) for i in (0, 2, 4)), 1.0


def _over(fg, bg):
    """Composite a (possibly translucent) token over an opaque token -> (r, g, b)."""
    (c, a), (b, _) = _rgba(fg), _rgba(bg)
    return tuple(round(c[i] * a + b[i] * (1 - a)) for i in range(3))


def _lum(c):
    f = lambda v: (v / 255) / 12.92 if v / 255 <= 0.03928 else ((v / 255 + 0.055) / 1.055) ** 2.4
    return 0.2126 * f(c[0]) + 0.7152 * f(c[1]) + 0.0722 * f(c[2])


def _cr(a, b):
    la, lb = sorted((_lum(a), _lum(b)), reverse=True)
    return (la + 0.05) / (lb + 0.05)


def _tok(P, name, over=None):
    c, a = _rgba(P[name])
    return _over(P[name], P[over]) if (a < 1 and over) else c


@pytest.fixture
def scale():
    saved = theme.UI_SCALE
    yield
    theme.UI_SCALE = saved
    _app.setStyleSheet("")


# ================================================================= theme: contrast tokens (B18, U26, -18)
@pytest.mark.parametrize("name", ["light", "dark"])
def test_text_tokens_meet_aa_on_every_surface_they_sit_on(name):
    P = theme._PAL[name]
    for fg in ("text", "text2", "dim", "faint", "accent", "good", "bad", "warn", "info"):
        for surf in ("bg", "panel", "panel2"):
            assert _cr(_tok(P, fg), _tok(P, surf)) >= 4.5, (name, fg, surf)
    # the status chips and notices are translucent tints over the surface they are drawn on
    for tone in ("good", "warn"):
        assert _cr(_tok(P, tone), _over(P[tone + "soft"], P["panel2"])) >= 4.5, (name, "chip", tone)   # classbn
    for tone in ("good", "warn", "bad", "info"):
        for surf in ("panel", "panel2"):                                                                 # card, dialog, banner
            assert _cr(_tok(P, tone), _over(P[tone + "soft"], P[surf])) >= 4.5, (name, "notice", tone, surf)
    assert _cr(_tok(P, "dim"), _tok(P, "panel3")) >= 4.5                                                  # Moderate chip
    assert _cr(_tok(P, "accentink"), _tok(P, "accent")) >= 4.5                                            # primary button


def test_pal_is_the_public_name_of_the_palette_table():
    assert theme.PAL is theme._PAL and set(theme.PAL) == {"light", "dark"}


def test_light_chip_inks_were_below_aa_before_and_now_clear_it():
    P = theme._PAL["light"]
    bg = _over(P["goodsoft"], P["panel2"])
    assert _cr(_rgba("#14774F")[0], bg) < 4.5                    # the originals (measured 4.28)
    assert _cr(_rgba(P["good"])[0], bg) >= 5.0
    assert _cr(_rgba("#7C621F")[0], _over(P["warnsoft"], P["panel2"])) < 4.5
    assert _cr(_rgba(P["warn"])[0], _over(P["warnsoft"], P["panel2"])) >= 5.0


@pytest.mark.parametrize("name", ["light", "dark"])
def test_placeholder_text_is_opaque_and_meets_aa(name, scale):
    _app.setStyleSheet(theme.qss(name))
    P = theme._PAL[name]
    for cls in (QLineEdit, QTextEdit, QPlainTextEdit):
        w = cls()
        w.ensurePolished()
        c = w.palette().placeholderText().color()
        assert c.alpha() == 255, cls.__name__                    # Qt's default is the text colour at 50% alpha
        assert _cr((c.red(), c.green(), c.blue()), _tok(P, "bg")) >= 4.5, (cls.__name__, name)


# ================================================================= theme: focus rings (B18, U30)
def _buttons(card):
    out = []
    for props in ({}, {"primary": True}, {"sm": True}, {"sm": True, "primary": True}, {"cellbtn": True}):
        b = QPushButton("RUN ANALYSIS")
        for k, v in props.items():
            b.setProperty(k, v)
        card.layout().addWidget(b)
        out.append(b)
    return out


@pytest.mark.parametrize("ui", [0.75, 0.85, 1.0, 1.1, 1.25, 1.5])
@pytest.mark.parametrize("name", ["light", "dark"])
def test_button_box_does_not_move_when_it_takes_focus(name, ui, scale):
    theme.UI_SCALE = ui
    _app.setStyleSheet(theme.qss(name))
    card = QFrame(); card.setObjectName("card"); QVBoxLayout(card)
    btns = _buttons(card)
    card.show(); _app.processEvents()
    before = [b.sizeHint() for b in btns]
    for b, sz in zip(btns, before):
        b.setFocus(Qt.TabFocusReason); _app.processEvents()
        assert b.sizeHint() == sz, (b.property("primary"), b.property("sm"), b.property("cellbtn"), ui)
        b.clearFocus(); _app.processEvents()
    card.close()


@pytest.mark.parametrize("name", ["light", "dark"])
def test_focused_button_ring_contrasts_3_to_1_with_the_panel(name, scale):
    theme.UI_SCALE = 1.0
    _app.setStyleSheet(theme.qss(name))
    P = theme._PAL[name]
    card = QFrame(); card.setObjectName("card"); QVBoxLayout(card)
    btns = _buttons(card)
    card.show(); QTest.qWaitForWindowActive(card)
    panel = _tok(P, "panel")
    for b in btns:
        b.setFocus(Qt.TabFocusReason); _app.processEvents()
        assert b.hasFocus()
        img = b.grab().toImage()
        for x in (0, 1):                                         # the 2px ring: outermost two pixel columns
            c = img.pixelColor(x, img.height() // 2)
            assert _cr((c.red(), c.green(), c.blue()), panel) >= 3.0, (b.property("primary"), b.property("sm"), x)
        b.clearFocus()
    card.close()


# ================================================================= theme: size floors, arrow, active theme
def test_label_size_floors_and_prose_face(scale):
    theme.UI_SCALE = 1.0
    _app.setStyleSheet(theme.qss("light"))
    for obj in ("kdim", "ver", "tagline", "sectionlabel", "gvpos", "cf", "cardmeta", "queuerow", "statusTxt", "secn"):
        lab = QLabel("x"); lab.setObjectName(obj); lab.ensurePolished()
        assert lab.font().pixelSize() >= 11, obj
    for obj in ("orient", "cardnote", "classexp", "empty", "errbanner", "notifmsg"):
        lab = QLabel("x"); lab.setObjectName(obj); lab.ensurePolished()
        assert lab.font().pixelSize() >= 12, obj
        assert lab.font().family() == fonts.BODY_FAMILY, obj
    from PySide6.QtWidgets import QHeaderView
    t = widgets.DataTable(["a", "b"]); t.ensurePolished(); t.horizontalHeader().ensurePolished()
    assert t.horizontalHeader().sizeHint().height() >= 11 + 12


def test_combo_box_draws_a_drop_down_arrow(scale):
    theme.UI_SCALE = 1.0
    _app.setStyleSheet(theme.qss("light"))
    cb = QComboBox(); cb.addItems(["qPCR", "Standard"]); cb.resize(160, 34); cb.show(); _app.processEvents()
    img = cb.grab().toImage()
    bg = QColor(theme._PAL["light"]["bg"])
    dark_in_arrow_zone = sum(1 for x in range(img.width() - 22, img.width() - 4) for y in range(8, img.height() - 8)
                             if img.pixelColor(x, y).lightness() < bg.lightness() - 60)
    assert dark_in_arrow_zone >= 8                                # a visible glyph where the drop-down sits
    cb.close()


def test_qss_records_the_active_theme_and_new_figure_panels_open_in_it(scale):
    theme.qss("light")
    assert theme.ACTIVE == "light"
    panel = widgets.FigurePanel(lambda bg: '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 10 10"/>', "f")
    assert panel.bg == "white"                                    # the app's light look, not modes[0]
    panel2 = widgets.FigurePanel(lambda bg: '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 10 10"/>', "f",
                                 app_theme="dark")
    assert panel2.bg == "dark"
    theme.qss("dark")
    assert widgets.FigurePanel(lambda bg: '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 10 10"/>',
                               "f").bg == "dark"
    for p in (panel, panel2):
        p.deleteLater()


# ================================================================= theme: domain hues (-13)
def test_every_domain_code_has_its_own_figure_hue_with_contrast_on_both_tracks():
    from teagle_core import domains
    codes = {v[0] for v in domains.DOMAIN_INFO.values()}
    assert codes <= set(theme.OK), codes - set(theme.OK)
    fallback = "#888"
    hues = {c: theme.OK[c].lower() for c in codes}
    assert fallback not in hues.values() and "#888888" not in hues.values()
    for c in ("EN", "ORF1", "HEL", "YR"):                        # the four that used to fall back to #888
        for track in (theme.GV_THEME["dark"]["track"], theme.GV_THEME["white"]["track"]):
            assert figures._contrast(theme.OK[c], track) >= 3.0, (c, track)
    # codes that can share an element must not share a hue (TPase and CHR share vermillion by design: DNA vs LTR)
    for a, b in (("EN", "ORF1"), ("EN", "RT"), ("ORF1", "RT"), ("YR", "RT"), ("YR", "RNaseH"), ("HEL", "TPase")):
        assert hues[a] != hues[b], (a, b)


def test_new_domain_hues_stay_apart_under_deuteranomaly_and_protanomaly():
    cs = pytest.importorskip("colorspacious")
    import numpy as np
    new = {c: theme.OK[c] for c in ("EN", "ORF1", "HEL", "YR")}
    others = {k: v for k, v in {**theme.OK, **theme.GENECOL, **theme.CISCOL, **theme.ARCHCOL}.items()
              if isinstance(v, str) and v.startswith("#") and k not in new and k not in ("ladder", "on", "off")}
    others["fallback"] = "#888888"
    rgb = lambda h: np.array([int(h[i:i + 2], 16) for i in (1, 3, 5)]) / 255.0
    for vis in (None, {"name": "sRGB1+CVD", "cvd_type": "deuteranomaly", "severity": 100},
                {"name": "sRGB1+CVD", "cvd_type": "protanomaly", "severity": 100}):
        def u(h):
            c = rgb(h)
            if vis:
                c = np.clip(cs.cspace_convert(c, vis, "sRGB1"), 0, 1)
            return cs.cspace_convert(c, "sRGB1", "CAM02-UCS")
        for code, h in new.items():
            for other, oh in {**others, **{k: v for k, v in new.items() if k != code}}.items():
                if {code, other} == {"EN", "HEL"}:                # never share an element; weak under deuteranomaly
                    continue
                assert np.linalg.norm(u(h) - u(oh)) >= 8.0, (vis and vis["cvd_type"], code, other)


# ================================================================= figures: shared helpers
def _rec(acc="M11240"):
    import engine
    seq = open(os.path.join(FIXTURES, acc + ".fasta")).read()
    return engine.run_analyze({"sequence": seq})["records"][0]


def _texts(svg):
    """[(x, y, size, anchor, text)] of every <text> element."""
    import xml.etree.ElementTree as ET
    root = ET.fromstring(svg)
    ns = "{http://www.w3.org/2000/svg}"
    return [(float(t.get("x")), float(t.get("y")), float(t.get("font-size", 12)), t.get("text-anchor", "start"),
             (t.text or "") + "".join(c.tail or "" for c in t)) for t in root.iter(ns + "text")]   # <title> excluded


def _text_box(t):
    x, y, size, anchor, text = t
    w = len(text) * size * 0.6
    x0 = x - w if anchor == "end" else (x - w / 2 if anchor == "middle" else x)
    return (x0, y - size * 0.8, x0 + w, y + size * 0.2)


def _overlaps(boxes):
    return [(a, b) for i, a in enumerate(boxes) for b in boxes[i + 1:]
            if a[0] < b[2] - 0.5 and b[0] < a[2] - 0.5 and a[1] < b[3] - 0.5 and b[1] < a[3] - 0.5]


def _raster(svg, scale=2, ground=0x00FFFFFF):
    w, h = widgets._svg_size(svg)
    img = QImage(int(w * scale), int(h * scale), QImage.Format_ARGB32)
    img.fill(ground)
    p = QPainter(img)
    QSvgRenderer(QByteArray(svg.encode())).render(p)
    p.end()
    return img


# ================================================================= genome viewer: export fills (-2)
def test_export_track_and_lane_fills_are_not_opaque_black():
    model = figures.gv_tracks_from_rec(_rec())
    svg, regions = figures.svg_genome(model, {"start": 0.0, "end": float(model["length"])}, 760, "dark",
                                      for_export=True, return_regions=True)
    img = _raster(svg, 2, 0x00FFFFFF)                            # transparent paper, as render_png does
    ML = figures.gv_gutter(model, 760)
    track_y = 34 + 13 + 10 + 24 + 8                              # inside the first track's band
    row = [img.pixelColor(x * 2, track_y * 2) for x in range(int(ML) + 3, 740)
           if not any(r["x0"] - 3 <= x <= r["x1"] + 3 and r["y0"] - 3 <= track_y <= r["y1"] + 3 for r in regions)]
    assert len(row) > 300
    opaque = sum(1 for c in row if c.alpha() > 100)              # grid lines only; a black track would be all of them
    assert opaque < 0.10 * len(row), (opaque, len(row))
    dark = sum(1 for x in range(0, img.width(), 5) for y in range(0, img.height(), 5)
               if img.pixelColor(x, y).alpha() > 200 and img.pixelColor(x, y).lightness() < 40)
    assert dark < 600, dark                                      # ink text only; four opaque black tracks are ~2000


def test_screen_themes_keep_opaque_tracks():
    model = figures.gv_tracks_from_rec(_rec())
    svg = figures.svg_genome(model, {"start": 0.0, "end": float(model["length"])}, 760, "white")
    assert 'fill="#f6f8fa" fill-opacity="1.0"' in svg


# ================================================================= genome viewer: contract (strand, kind, protein, dna, index)
def test_domain_features_and_regions_carry_strand_kind_protein_dna_and_index():
    rec = _rec()
    model = figures.gv_tracks_from_rec(rec)
    dom = next(t for t in model["tracks"] if t["name"] == "protein domains")["features"]
    assert len(dom) == len(rec["domains"])
    for i, (f, d) in enumerate(zip(dom, rec["domains"])):
        assert (f["kind"], f["strand"], f["protein"], f["dna"], f["index"]) == ("domain", d["strand"],
                                                                                d["protein"], d["dna"], i)
    svg, regions = figures.svg_genome(model, {"start": 0.0, "end": float(model["length"])}, 760, "white",
                                      return_regions=True)
    drs = [r for r in regions if r["kind"] == "domain"]
    assert len(drs) == len(dom)
    for r, d in zip(drs, rec["domains"]):
        assert r["strand"] == d["strand"] and r["protein"] == d["protein"] and r["dna"] == d["dna"]
        assert r["index"] in range(len(dom)) and r["label"] == d["domain"]
    assert all("kind" in r and "strand" in r for r in regions)


def test_minus_strand_domains_keep_their_strand_and_get_an_arrowhead():
    from teagle_core.primers import reverse_complement
    import engine
    seq = "".join(l for l in open(os.path.join(FIXTURES, "M11240.fasta")).read().splitlines() if l[:1] != ">")
    rec = engine.run_analyze({"sequence": reverse_complement(seq)})["records"][0]
    model = figures.gv_tracks_from_rec(rec)
    dom = next(t for t in model["tracks"] if t["name"] == "protein domains")["features"]
    assert dom and all(f["strand"] == "-" for f in dom)
    svg = figures.svg_genome(model, {"start": 0.0, "end": float(model["length"])}, 760, "white")
    assert svg.count("<path") >= len(dom)                        # one arrowhead per wide-enough domain glyph
    assert all("strand -" in f["tip"] for f in dom)


def test_gene_model_features_carry_strand_and_derived_introns_render_as_intron_star():
    gm = {"exons": [{"start": 0, "end": 100, "strand": "-"}, {"start": 300, "end": 400, "strand": "-"}],
          "introns": [{"start": 100, "end": 300, "strand": "-", "derived": True}],
          "cds": [{"start": 20, "end": 380, "strand": "-"}]}
    model = figures.gv_tracks_from_gene(gm, 500)
    feats = model["tracks"][0]["features"]
    intron = next(f for f in feats if f.get("intron"))
    assert intron["derived"] is True and intron["label"] == "intron*" and intron["strand"] == "-"
    assert all(f["strand"] == "-" for f in feats)
    assert model["tracks"][1]["features"][0]["strand"] == "-"
    svg, regions = figures.svg_genome(model, {"start": 0.0, "end": 500.0}, 620, "white", return_regions=True)
    r = next(r for r in regions if r["kind"] == "intron")
    assert r["label"] == "intron*" and r["derived"] is True and r["strand"] == "-"
    annotated = figures.gv_tracks_from_gene({"exons": [{"start": 0, "end": 10}],
                                             "introns": [{"start": 10, "end": 50}], "cds": []}, 60)
    assert annotated["tracks"][0]["features"][1]["label"] == "intron"


def test_overview_label_reads_whole_record():
    model = figures.gv_tracks_from_rec(_rec())
    svg = figures.svg_genome(model, {"start": 0.0, "end": float(model["length"])}, 760, "white")
    assert "whole record · 5,146 bp" in svg and "whole element" not in svg


# ================================================================= genome viewer: TSD and the advisory pair (-11)
def test_tsd_copies_are_drawn_with_their_own_kind():
    rec = {"composition": {"length": 3000},
           "structural": [{"type": "LTR (terminal direct repeat)", "five_prime": [105, 400], "three_prime": [2700, 2995]},
                          {"type": "TSD (target-site duplication)", "length": 5, "motif": "ACGTA",
                           "upstream": [100, 105], "downstream": [2995, 3000], "chance_level": False}],
           "domains": [], "orfs": []}
    model = figures.gv_tracks_from_rec(rec)
    tsd = next(t for t in model["tracks"] if t["name"] == "TSD")["features"]
    assert [(f["start"], f["end"], f["kind"]) for f in tsd] == [(100, 105, "TSD"), (2995, 3000, "TSD")]
    svg, regions = figures.svg_genome(model, {"start": 0.0, "end": 3000.0}, 760, "white", return_regions=True)
    assert [r["kind"] for r in regions if r["kind"] == "TSD"] == ["TSD", "TSD"]


def test_sub_threshold_pair_is_drawn_open_and_cannot_pass_for_an_ltr():
    rec = {"composition": {"length": 3000},
           "structural": [{"type": "Sub-threshold terminal direct repeat (advisory)", "ltr_len": 300, "identity": 61.0,
                           "threshold": 80, "five_prime": [0, 300], "three_prime": [2700, 3000], "advisory": True}],
           "domains": [], "orfs": []}
    reps = figures.gv_tracks_from_rec(rec)["tracks"][0]["features"]
    assert len(reps) == 2 and all(f["outline"] == "sub" and f["label"] == "sub-threshold" for f in reps)
    assert all(f["kind"].split(" ")[0] == "LTR" for f in reps)   # the menu offers copy actions only, never primer design
    svg = figures.svg_genome({"length": 3000, "tracks": [{"name": "terminal repeats", "height": 20, "features": reps}]},
                             {"start": 0.0, "end": 3000.0}, 760, "white")
    assert 'stroke-dasharray="3 2"' in svg and 'fill-opacity="0.15"' in svg


def test_withheld_terminal_repeat_candidate_is_labelled_tdr_not_ltr():
    rec = {"composition": {"length": 3000},
           "classification": {"te_class": "repeat/structural-only",
                              "completeness": {"kind": "terminal direct-repeat candidate (structure-only)"}},
           "structural": [{"type": "LTR (terminal direct repeat)", "five_prime": [0, 300], "three_prime": [2700, 3000]}],
           "domains": [], "orfs": []}
    reps = figures.gv_tracks_from_rec(rec)["tracks"][0]["features"]
    assert {f["label"] for f in reps} == {"TDR"} and {f["kind"] for f in reps} == {"LTR"}


# ================================================================= genome viewer: legibility (B25, U25, -12)
def test_orf_rows_are_marked_plus_and_minus_and_domain_bearing_orfs_are_outlined():
    model = figures.gv_tracks_from_rec(_rec())
    orfs = next(t for t in model["tracks"] if t.get("stranded"))
    assert orfs["name"] == "ORFs"
    svg = figures.svg_genome(model, {"start": 0.0, "end": float(model["length"])}, 760, "white")
    texts = _texts(svg)
    plus, minus = [t for t in texts if t[4] == "+"], [t for t in texts if t[4] == "−"]
    assert len(plus) == 1 and len(minus) == 1 and plus[0][1] < minus[0][1]
    name = next(t for t in texts if t[4] == "ORFs")
    assert _text_box(name)[2] <= _text_box(plus[0])[0] - 2        # the name clears the mark
    emph = [f for f in orfs["features"] if f["emph"]]
    assert emph and all("carries" in f["tip"] for f in emph)
    assert svg.count('stroke-width="1.2"') >= len(emph)


def test_orf_attribution_fails_closed_on_a_strand_or_interval_mismatch():
    rec = {"orfs": [{"start": 0, "end": 600, "strand": "+", "frame": 0, "length_aa": 200},
                    {"start": 700, "end": 1300, "strand": "+", "frame": 1, "length_aa": 200}],
           "domains": [{"domain": "RT", "orf": 0, "strand": "+", "nt": [100, 400]},
                       {"domain": "INT", "orf": 1, "strand": "+", "nt": [100, 400]},        # outside ORF 1
                       {"domain": "RT", "orf": 1, "strand": "-", "nt": [800, 900]},         # wrong strand
                       {"domain": "RT", "orf": 9, "strand": "+", "nt": [800, 900]}]}        # no such ORF
    assert figures._domain_orf_map(rec) == {0: ["RT"]}


def test_short_terminal_repeats_are_never_narrower_than_three_pixels_and_stay_inside_the_plot():
    rec = {"composition": {"length": 4565},
           "structural": [{"type": "TIR (terminal inverted repeat)", "five_prime": [0, 11], "three_prime": [4554, 4565]}],
           "domains": [], "orfs": []}
    model = figures.gv_tracks_from_rec(rec)
    W = 620
    svg, regions = figures.svg_genome(model, {"start": 0.0, "end": 4565.0}, W, "white", return_regions=True)
    ML = figures.gv_gutter(model, W)
    plot = max(120, W - ML - figures.GV_MR)
    assert len(regions) == 2
    for r in regions:
        assert r["x1"] - r["x0"] >= figures.MIN_GLYPH_PX - 1e-6
        assert r["x0"] >= ML - 1e-6 and r["x1"] <= ML + plot + 1e-6
    glyphs = re.findall(r'class="gvglyph" x="([\d.]+)" y="[\d.]+" width="([\d.]+)"', svg)
    assert [(round(float(x), 1), round(float(w), 1)) for x, w in glyphs] == [
        (round(r["x0"], 1), round(r["x1"] - r["x0"], 1)) for r in regions]      # hover boxes match what is drawn


@pytest.mark.parametrize("export", [False, True])
def test_figure_text_floor_and_no_overlapping_text(export):
    floor = figures.FIG_PX_EXPORT if export else figures.FIG_PX_SCREEN
    for acc in ("M11240", "X05424", "M80343"):
        model = figures.gv_tracks_from_rec(_rec(acc))
        for W in (320, 620, 920):
            svg = figures.svg_genome(model, {"start": 0.0, "end": float(model["length"])}, W, "white", export,
                                     legend=export)
            assert min(t[2] for t in _texts(svg)) >= floor, (acc, W)
            assert not _overlaps([_text_box(t) for t in _texts(svg)]), (acc, W, export)


def test_ruler_tick_count_follows_the_width_so_labels_never_collide():
    model = {"length": 230000, "tracks": []}
    for W in (320, 620, 920):
        svg = figures.svg_genome(model, {"start": 0.0, "end": 230000.0}, W, "white")
        ticks = [t for t in _texts(svg) if re.fullmatch(r"[\d.]+[kM]?", t[4])]
        assert ticks and not _overlaps([_text_box(t) for t in ticks]), W


def test_narrow_feature_captions_wrap_and_are_never_truncated():
    feats = [{"start": 100 * i, "end": 100 * i + 2, "color": "#0072B2", "label": f"F{i}"} for i in range(1, 30)]
    model = {"length": 3000, "tracks": [{"name": "t", "height": 18, "features": feats}]}
    svg = figures.svg_genome(model, {"start": 0.0, "end": 3000.0}, 620, "white")
    assert "more" not in svg
    joined = " ".join(t[4] for t in _texts(svg))
    for i in range(1, 30):
        assert f"F{i} {100 * i}–{100 * i + 2}" in joined, i       # every label, with its plain-digit interval
    assert all(_text_box(t)[2] <= 620 - figures.GV_MR + 1 for t in _texts(svg) if t[4].startswith("F"))   # wrapped


def test_caption_coordinates_are_plain_digits_matching_the_table_form():
    model = figures.gv_tracks_from_rec(_rec())
    svg = figures.svg_genome(model, {"start": 0.0, "end": float(model["length"])}, 760, "white")
    caps = [t[4] for t in _texts(svg) if "PPT" in t[4]]
    assert caps and "4861–4870" in caps[0] and not re.search(r"\d,\d{3}–", caps[0])
    assert figures.fmt_interval(4861, 4870) == "4861–4870"


def _glyph_attrs(svg):
    return [dict(re.findall(r'([\w-]+)="([^"]*)"', tag)) for tag in re.findall(r'<rect class="gvglyph"[^>]*>', svg)]


def test_light_and_export_glyphs_have_an_outline_of_their_own_hue_that_holds_contrast():
    model = figures.gv_tracks_from_rec(_rec())
    for theme_name, export, ground in (("white", False, theme.GV_THEME["white"]["track"]), ("dark", True, "#FFFFFF")):
        svg = figures.svg_genome(model, {"start": 0.0, "end": float(model["length"])}, 760, theme_name, export)
        glyphs = _glyph_attrs(svg)
        assert glyphs
        for g in glyphs:
            stroke = g.get("stroke")
            assert stroke, (theme_name, g["fill"])
            if g.get("stroke-width") == "1.2":                    # an ORF that carries a domain: figure-ink outline
                continue
            assert figures._contrast(stroke, ground) >= 4.5, (theme_name, g["fill"], stroke)
            fh = figures.colorsys.rgb_to_hls(*[v / 255 for v in figures._rgb(g["fill"])])
            sh = figures.colorsys.rgb_to_hls(*[v / 255 for v in figures._rgb(stroke)])
            assert abs(fh[0] - sh[0]) < 0.02 or fh[2] < 0.05, (g["fill"], stroke)    # same hue
    dark = _glyph_attrs(figures.svg_genome(model, {"start": 0.0, "end": float(model["length"])}, 760, "dark"))
    assert dark and all(g.get("stroke-width") in (None, "1.2") for g in dark)      # the dark viewer is unchanged


def test_every_exported_genome_figure_carries_a_legend_of_what_it_draws_and_the_screen_does_not():
    model = figures.gv_tracks_from_rec(_rec())
    view = {"start": 0.0, "end": float(model["length"])}
    shown = figures.svg_genome(model, view, 760, "white")
    exported = figures.svg_genome(model, view, 760, "dark", True, legend=True)
    assert "positions in bp from the sequence start (0-based)" not in shown
    assert "positions in bp from the sequence start (0-based)" in exported
    text = " ".join(t[4] for t in _texts(exported))
    for name in ("LTR", "PPT", "INT integrase", "RT reverse transcriptase", "ORF carrying a detected domain", "ORF"):
        assert name in text, name
    swatches = re.findall(r'<rect x="[\d.]+" y="[\d.]+" width="9" height="9"', exported)
    drawn = {(f["color"], f["legend"]) for t in model["tracks"] for f in t["features"]}
    assert len(swatches) == len(drawn)
    w, _h = widgets._svg_size(exported)
    assert all(_text_box(t)[2] <= w for t in _texts(exported))


def test_panel_exports_use_the_legend_and_the_publication_palette_by_default():
    panel = widgets.GenomePanel(figures.svg_genome, "g")
    panel.set_model(figures.gv_tracks_from_rec(_rec()))
    tgt = panel._export_target()
    assert tgt["for_export"] is True and tgt["legend"] is True
    svg = panel._cur_svg(920, **tgt)
    assert "positions in bp from the sequence start (0-based)" in svg
    panel.deleteLater()


# ================================================================= gel: lane status, wording, bands, ladder
def _amp(length, on=None, intended=None, single=False, start=None, **kw):
    a = {"length": length, "product_length": length, "on_target": on, "single_primer": single,
         "fwd_mm": 0, "rev_mm": 0, "source": "chrT"}
    if intended is not None:
        a["matches_intended_product"] = intended
    if start is not None:
        a["start"], a["end"] = start, start + length
    a.update(kw)
    return a


def _gel_text(svg):
    return [t[4] for t in _texts(svg)]


@pytest.mark.parametrize("bg", ["dark", "white", "uv", "mono", "transparent"])
def test_failed_lane_is_hatched_and_captioned_failed_never_a_bare_dash(bg):
    lanes = [{"label": "P1", "status": "failed", "amplicons": []},
             {"label": "P2", "status": "complete_zero_products", "amplicons": []}]
    svg = figures.svg_gel({"lanes": lanes}, bg)
    texts = _gel_text(svg)
    assert "failed" in texts and "—" in texts                    # P1 says failed; P2 is the genuine empty lane
    assert texts.count("—") == 1                                 # failed is NOT drawn as an empty lane
    assert "<path" in svg and "stroke-dasharray=\"3 2\"" in svg   # hatched
    assert "failed (not assessed)" in texts                      # the legend carries the state
    assert QSvgRenderer(QByteArray(svg.encode())).isValid()


def test_capped_lane_gets_a_badge_and_never_claims_no_on_target():
    lanes = [{"label": "P1", "status": "partial_capped", "has_locus": True,
              "amplicons": [_amp(300, False, False, start=10)]}]
    svg = figures.svg_gel({"lanes": lanes}, "white")
    texts = _gel_text(svg)
    assert "capped" in texts and "no on-target" not in texts
    assert "capped (lower bound)" in texts
    lane_tip = re.search(r"<title>(The search hit its cap[^<]*)</title>", svg)
    assert lane_tip and "no specificity conclusion is supported" in lane_tip.group(1)


def test_lane_without_a_designed_target_never_prints_a_red_no_on_target():
    lanes = [{"label": "P1", "has_locus": False, "amplicons": [_amp(300, None, start=10)]}]
    svg = figures.svg_gel({"lanes": lanes}, "white")
    assert "no on-target" not in _gel_text(svg)
    assert "No designed target was set for this lane" in svg      # the reason is in the lane tooltip
    assert "priming site" in _gel_text(svg)


def test_design_bound_lane_without_a_target_says_designed_product_not_a_genomic_site():
    lanes = [{"label": "P1", "has_locus": False,
              "amplicons": [_amp(300, None, True, start=10), _amp(500, None, False, start=900)]}]
    svg = figures.svg_gel({"lanes": lanes}, "white")
    tips = " ".join(re.findall(r"<title>([^<]*)</title>", svg)) + " " + " ".join(r["tip"] for r in figures.gel_regions({"lanes": lanes}))
    legend = _gel_text(svg)
    assert "designed product (target not assessed)" in legend and "priming site" in legend
    assert "set no target" in tips
    only = [{"label": "P1", "has_locus": False, "amplicons": [_amp(300, None, True, start=10)]}]
    assert "priming site" not in _gel_text(figures.svg_gel({"lanes": only}, "white"))
    assert "genomic" not in " ".join(re.findall(r"<title>([^<]*)</title>", figures.svg_gel({"lanes": only}, "white")))
    assert any("designed product (target not assessed)" in r["tip"] for r in figures.gel_regions({"lanes": lanes}))
    scan = [{"label": "S", "has_locus": False, "amplicons": [_amp(300, None, start=10)]}]
    assert "priming site" in _gel_text(figures.svg_gel({"lanes": scan}, "white"))
    assert "genomic priming sites" in figures.svg_gel({"lanes": scan}, "white")


def test_a_complete_lane_whose_bands_are_all_unintended_still_says_no_on_target():
    lanes = [{"label": "P1", "status": "complete", "has_locus": True, "amplicons": [_amp(300, False, False, start=1)]}]
    assert "no on-target" in _gel_text(figures.svg_gel({"lanes": lanes}, "white"))


def test_legend_lists_only_the_states_the_gel_uses():
    only_on = [{"label": "P1", "has_locus": True, "amplicons": [_amp(300, True, True, start=1)]}]
    t = _gel_text(figures.svg_gel({"lanes": only_on}, "white"))
    assert "on-target" in t and "failed (not assessed)" not in t and "capped (lower bound)" not in t
    assert "single-primer" not in t
    mixed = only_on + [{"label": "P2", "status": "failed", "amplicons": []},
                       {"label": "P3", "status": "partial_capped", "has_locus": True,
                        "amplicons": [_amp(400, False, False, True, start=1)]}]
    t = _gel_text(figures.svg_gel({"lanes": mixed}, "white"))
    assert {"failed (not assessed)", "capped (lower bound)", "single-primer"} <= set(t)


def test_design_bound_extras_are_additional_products_not_off_targets():
    lanes = [{"label": "P1", "has_locus": True,
              "amplicons": [_amp(300, True, True, start=1), _amp(500, False, False, start=900)]}]
    svg = figures.svg_gel({"lanes": lanes}, "white")
    regions = figures.gel_regions({"lanes": lanes})
    tips = " ".join(r["tip"] for r in regions) + " " + " ".join(re.findall(r"<title>([^<]*)</title>", svg))
    assert "1 additional product" in tips and "off-target" not in tips
    legend = _gel_text(svg)
    assert "additional product" in legend and "off-target" not in legend
    # a genome-scan product (no engine verdict on designed endpoints) keeps the wording the app itself defines
    scan = [{"label": "S", "has_locus": True, "amplicons": [_amp(300, True, start=1), _amp(900, False, start=800)]}]
    assert "1 off-target" in figures.gel_regions({"lanes": scan})[1]["tip"]
    assert "off-target" in _gel_text(figures.svg_gel({"lanes": scan}, "white"))


def test_unassessed_products_are_neutral_not_off_target():
    g = [_amp(300, None, start=1), _amp(300, None, start=500)]
    n = figures._call_split(g, True)
    assert n == {"on": 0, "additional": 0, "off": 0, "single": 0, "site": 2}
    bands = figures._lane_bands(g, figures.GELPAL["dark"], True)
    assert bands[0]["call"] == "site" and bands[0]["on"] is False


def test_near_equal_products_merge_into_one_band_drawn_in_the_worst_case_colour():
    P = figures.GELPAL["white"]
    lane = [_amp(502, True, True, start=0), _amp(498, False, False, start=900), _amp(496, False, False, start=2000)]
    bands = figures._lane_bands(lane, P, True)
    assert len(bands) == 1 and bands[0]["call"] == "off" and bands[0]["color"] == P["off"]
    assert bands[0]["size"] == 502 and bands[0]["size_text"] == "496–502 bp" and bands[0]["count"] == 3
    regs = figures.gel_regions({"lanes": [{"label": "P1", "has_locus": True, "amplicons": lane}]})
    assert len(regs) == 1 and "496–502 bp" in regs[0]["tip"] and "3 products drawn as one band" in regs[0]["tip"]
    assert regs[0]["amplicon"]["on_target"] is True               # the menu still addresses the designed product
    far = figures._lane_bands([_amp(500, True, True, start=0), _amp(530, False, False, start=900)], P, True)
    assert len(far) == 2                                          # 6% apart: two bands


def test_gel_ladder_extends_to_cover_large_genome_scan_products():
    lanes = [{"label": "S", "has_locus": False, "amplicons": [_amp(3800, None, start=0)]}]
    G = figures._gel_geometry({"lanes": lanes})
    assert 3000 in G["LADDER"] and 2000 in G["LADDER"] and G["maxbp"] == 4000
    svg = figures.svg_gel({"lanes": lanes}, "white")
    assert "4000" in _gel_text(svg) and "L = size ladder, 100–4000 bp" in _gel_text(svg)
    small = figures._gel_geometry({"lanes": [{"label": "P", "amplicons": [_amp(300, True)]}]})
    assert small["maxbp"] == 1600 and max(small["LADDER"]) == 1500   # unchanged for ordinary products


def test_gel_tooltips_keep_a_validated_chromosome_accession_readable():
    a = _amp(300, True, True, start=0, source="NC_000001.11", sourceIdentitySource="validated_fetch")
    assert figures._public_amplicon_source(a) == "NC_000001.11"
    assert "NC_000001.11" in figures.gel_regions({"lanes": [{"label": "P", "has_locus": True, "amplicons": [a]}]})[0]["tip"]
    local = _amp(300, True, True, start=0, source="/home/person/private/sample.fa")
    assert "person" not in figures._public_amplicon_source(local)   # a local path is still rewritten


@pytest.mark.parametrize("bg", ["dark", "white", "uv", "mono", "transparent"])
def test_gel_text_floor_pitch_and_no_overlap_in_every_palette(bg):
    lanes = [{"label": f"P{i+1}", "status": "partial_capped" if i == 3 else "complete", "has_locus": True,
              "amplicons": [_amp(300 + 60 * i, False, False, start=0)]} for i in range(10)]
    for case in (lanes, lanes[:2], [lanes[0]]):
        svg = figures.svg_gel({"lanes": case}, bg)
        texts = _texts(svg)
        assert min(t[2] for t in texts) >= figures.GEL_PX
        assert not _overlaps([_text_box(t) for t in texts]), (bg, len(case))
        w, _h = widgets._svg_size(svg)
        assert all(_text_box(t)[2] <= w + 0.5 for t in texts)
    # a one-lane gel is as wide as its legend needs
    one = widgets._svg_size(figures.svg_gel({"lanes": [lanes[0]]}, bg))[0]
    assert one >= 300


def test_gel_footer_height_comes_from_the_tallest_caption_so_rows_do_not_collide():
    lanes = [{"label": "P1", "amplicons": [_amp(300, True, start=0)], "has_locus": True}]
    g1 = figures._gel_geometry({"lanes": lanes})
    assert g1["foot"] >= 1 and g1["ROW_BOT_PAD"] >= 20
    regs = figures.gel_regions({"lanes": lanes * 11})              # wrapped rows keep matching hit-boxes
    assert len(regs) == 11 and max(r["y0"] for r in regs) - min(r["y0"] for r in regs) > 100


# ================================================================= TE landscape (-3, -15)
def _te_result(**kw):
    r = {"genome_bp": 1_000_000, "total_bp": 460_000, "te_bp": 451_600, "tandem_bp": 1000, "other_bp": 12_900,
         "library_kind": "installed Dfam partitions, curated + uncurated", "dfam_version": "3.8"}
    r.update(kw)
    return r


def test_te_percentage_is_never_painted_over_by_its_own_bar():
    r = _te_result()
    svg = figures.svg_te_composition(r)
    import xml.etree.ElementTree as ET
    root = ET.fromstring(svg)
    ns = "{http://www.w3.org/2000/svg}"
    bar_x = min(float(e.get("x")) for e in root.iter(ns + "rect") if e.get("stroke") and e.get("height") == "18")
    labels = [t for t in root.iter(ns + "text") if re.search(r"\d\.\d\d%$", "".join(t.itertext()))]
    assert len(labels) == 5
    for t in labels:
        label = "".join(t.itertext())
        assert float(t.get("x")) + len(label) * float(t.get("font-size")) * 0.6 <= bar_x, label
    assert figures.svg_te_composition(_te_result(other_bp=1_000_000, total_bp=1_000_000, te_bp=1_000_000,
                                                 tandem_bp=1_000_000))   # builds at 100% for every row


def test_te_footers_wrap_inside_the_figure_even_for_a_long_custom_library():
    r = _te_result(library_kind="user-supplied FASTA library",
                   custom_library={"name": "a_really_long_lab_reference_library_name_" * 3 + ".fa", "sha256": "ab" * 32})
    svg = figures.svg_te_composition(r)
    w, h = widgets._svg_size(svg)
    texts = _texts(svg)
    assert all(_text_box(t)[2] <= w for t in texts), [t[4][-30:] for t in texts if _text_box(t)[2] > w]
    assert all(_text_box(t)[3] <= h for t in texts)
    joined = " ".join(t[4] for t in texts)
    assert "ab" * 32 in joined and "simple/low-complexity is separate." in joined      # nothing cut off
    assert not _overlaps([_text_box(t) for t in texts])
    default = figures.svg_te_composition(_te_result(library_kind="user-supplied FASTA library",
                                                    custom_library={"name": "lab.fa", "sha256": "a5" * 32}))
    assert widgets._svg_size(default) == (900.0, 244.0)           # the ordinary custom-library figure is unchanged


# ================================================================= dot plot and heat map (-8, -9, B25)
def _dot_matrix():
    from teagle_core import dotplot, examples
    return dotplot.self_matrix(examples.load("M11240"), k=dotplot.suggest_k([]))


def test_dot_plot_ticks_use_the_ruler_formatter_and_text_sizes_hold_the_floor():
    m = _dot_matrix()
    for export in (False, True):
        svg = figures.svg_dotplot(m, W=920 if export else 620, for_export=export)
        texts = _texts(svg)
        assert min(t[2] for t in texts) >= (figures.FIG_PX_EXPORT if export else figures.FIG_PX_SCREEN)
        assert any(re.fullmatch(r"\d(\.\d)?k", t[4]) for t in texts) and not any("," in t[4] and t[4][0].isdigit()
                                                                                for t in texts if "bp" not in t[4])
        assert not _overlaps([_text_box(t) for t in texts]), export


def test_dot_export_footer_wraps_to_the_plot_width_and_the_page_grows_to_hold_it():
    m = _dot_matrix()
    m = {**m, "masked_kmers": 12, "masked_positions": 400, "max_occ": 30, "truncated": True}
    on = figures.svg_dotplot(m, W=920, for_export=True)
    off = figures.svg_dotplot(m, W=920, for_export=False)
    w, h = widgets._svg_size(on)
    foot = [t for t in _texts(on) if t[1] > 46 + (920 - 58 - 16) + 40]
    assert foot and min(t[2] for t in foot) >= 10
    assert all(_text_box(t)[2] <= 58 + (920 - 58 - 16) + 1 for t in foot)      # wrapped to the plot width
    assert all(_text_box(t)[3] <= h for t in _texts(on))
    assert widgets._svg_size(on)[1] > widgets._svg_size(off)[1]
    assert "Comparison truncated at the work limit" in " ".join(t[4] for t in foot)


def test_dot_legend_second_entry_starts_after_the_first_ends():
    m = _dot_matrix()
    for px in (9, 10):
        svg = figures.svg_dotplot(m, W=920 if px == 10 else 620, for_export=px == 10)
        texts = {t[4]: t for t in _texts(svg)}
        first, second = texts["direct repeat (forward)"], texts["inverted repeat (reverse complement)"]
        assert _text_box(first)[2] < second[0] - 10


def test_heat_map_on_the_dark_theme_ramps_from_the_background_so_dense_bins_are_brightest():
    m = _dot_matrix()
    svg = figures.svg_dotheat(m, W=620, theme="dark", which="forward")
    bg = figures._DOT_THEME["dark"]["bg"]
    fills = re.findall(r'<rect x="[\d.]+" y="[\d.]+" width="[\d.]+" height="[\d.]+" fill="(#[0-9A-F]{6})"/>', svg)
    cells = [f for f in fills if f.upper() != bg.upper() and f.upper() != "#0F1519"]
    assert len(cells) > 20
    # every drawn bin lies between the background and the signal hue; none (a lone match included) is brighter than
    # the peak colour, which a near-white floor on the dark theme made the most salient mark of all
    assert max(figures._lum(f) for f in cells) <= figures._lum(figures._DOT_FWD) + 0.01
    low = figures._heat_ramp(1 / max(m["forward_max"], 2), figures._DOT_FWD, bg)
    assert figures._contrast(low, bg) < figures._contrast(figures._DOT_FWD, bg)
    light = figures.svg_dotheat(m, W=620, theme="white", which="forward")
    assert figures._HEAT_FLOOR.upper() in light.upper()            # light and export keep the near-white floor
    assert figures._heat_ramp(1.0, "#7B1FA2") == "#7B1FA2"


def test_dot_canvas_hover_readout_is_not_replaced_by_the_help_tooltip(monkeypatch):
    panel = widgets.DotPanel()
    panel.set_matrix(_dot_matrix(), threshold=1)
    assert panel.canvas.toolTip() == ""                            # nothing to overwrite the per-cell readout
    assert panel.helpBtn.toolTip() == widgets.DotPanel.HELP and panel.helpBtn.accessibleName()
    panel.resize(900, 900); panel.show(); _app.processEvents()
    c = panel.canvas
    sc = c._sc()
    shown = []
    monkeypatch.setattr(widgets, "QToolTip", types.SimpleNamespace(
        showText=lambda pos, text, w=None: shown.append(text), hideText=lambda: None))
    pos = QPointF((figures.DOT_ML + 100) * sc, (figures.DOT_MT + 100) * sc)
    c.mouseMoveEvent(QMouseEvent(QEvent.MouseMove, pos, c.mapToGlobal(pos), Qt.NoButton, Qt.NoButton, Qt.NoModifier))
    assert shown and "bp" in shown[-1] and "," not in shown[-1].splitlines()[0]
    panel.close()


# ================================================================= widgets: casing (B30, U34)
def test_caps_keeps_units_and_is_idempotent():
    assert widgets.caps("Prod min (bp)") == "PROD MIN (bp)"
    assert widgets.caps("Opt Tm (°C)") == "OPT Tm (°C)"
    assert widgets.caps("ORFs (≥40 aa) — 1,234") == "ORFS (≥40 aa) — 1,234"
    assert widgets.caps("Organism — local isPcr") == "ORGANISM — LOCAL isPcr"
    assert widgets.caps("ΔG kcal/mol, 5 nt and 150bp") == "ΔG kcal/mol, 5 nt AND 150bp"
    assert widgets.caps("run analysis") == "RUN ANALYSIS" and widgets.caps("−") == "−"
    for t in ("Prod min (bp)", "Opt Tm (°C)", "run analysis", "minimap2 and mRNA", "x"):
        assert widgets.caps(widgets.caps(t)) == widgets.caps(t)
    assert {"bp", "kb", "Mb", "nt", "aa", "Tm", "ΔG", "kcal", "mol", "mRNA", "cDNA", "tRNA", "isPcr",
            "minimap2"} <= widgets.KEEP_CASE


def test_set_caps_applies_to_action_buttons_only_and_uppercase_buttons_includes_the_root():
    act, link, hdr = QPushButton("x"), QPushButton("x"), QPushButton("x")
    link.setProperty("link", True)
    hdr.setObjectName("cardhdr")
    widgets.set_caps(act, "Run analysis (30 bp)")
    widgets.set_caps(link, "Show all 24 ORFs")
    widgets.set_caps(hdr, "Primer design")
    assert (act.text(), link.text(), hdr.text()) == ("RUN ANALYSIS (30 bp)", "Show all 24 ORFs", "Primer design")
    solo = QPushButton("Show")
    widgets.uppercase_buttons(solo)                                # the root itself, with no parent to search under
    assert solo.text() == "SHOW"
    box = QWidget(); lay = QVBoxLayout(box)
    b = QPushButton("Prod min (bp)"); lay.addWidget(b)
    widgets.uppercase_buttons(box)
    assert b.text() == "PROD MIN (bp)"


def test_panel_toolbar_buttons_keep_their_casing_when_their_text_changes():
    panel = widgets.DotPanel()
    panel.set_matrix(_dot_matrix(), threshold=1)
    panel._toggle_guides(); panel._toggle_layer(); panel._set_mode("heat"); panel._toggle_layer()
    for b in (panel.guideBtn, panel.layerBtn):
        assert b.text() == widgets.caps(b.text()), b.text()
    assert panel.guideBtn.text() == "SHOW GUIDES"


def test_sortkey_reads_thousands_separators_as_part_of_the_number():
    assert widgets._sortkey("1,500") == 1500.0 and widgets._sortkey("230,000") == 230000.0
    assert widgets._sortkey("5,146") == 5146.0 and widgets._sortkey("60.1/58.3") == pytest.approx(118.4)
    assert widgets._sortkey("0–5146") == 0.0 and widgets._sortkey("2.0e-11") == pytest.approx(2e-11)
    t = widgets.DataTable(["len"])
    t.set_rows([["1,500"], ["999"], ["12,000"]])
    t.sortItems(0, Qt.AscendingOrder)
    assert [t.item(i, 0).text() for i in range(3)] == ["999", "1,500", "12,000"]


# ================================================================= widgets: DataTable contract
def _table(rows=6, cols=("id", "val")):
    t = widgets.DataTable(list(cols))
    t.set_rows([[f"r{i}", (i * 7) % 5] for i in range(rows)])
    t.resize(300, 300); t.show(); QTest.qWaitForWindowActive(t)
    return t


def test_current_row_changed_reports_the_original_index_for_keyboard_and_mouse_after_a_sort():
    t = _table()
    got = []
    t.current_row_changed.connect(got.append)
    t.sortItems(1, Qt.DescendingOrder)                             # visual order no longer equals original order
    visual = [t.item(i, 0).text() for i in range(t.rowCount())]
    assert visual != [f"r{i}" for i in range(6)]
    t.setFocus(); t.setCurrentCell(0, 0); got.clear()
    QTest.keyClick(t, Qt.Key_Down)
    assert got and got[-1] == int(visual[1][1:])                  # r<orig> of the second visual row
    assert t.current_orig() == got[-1]
    QTest.keyClick(t, Qt.Key_Down)
    assert got[-1] == int(visual[2][1:])
    rect = t.visualItemRect(t.item(4, 0))
    QTest.mouseClick(t.viewport(), Qt.LeftButton, pos=rect.center())
    assert got[-1] == int(visual[4][1:])
    n = len(got)
    QTest.mouseClick(t.viewport(), Qt.LeftButton, pos=rect.center())
    assert len(got) == n                                          # same row again: no duplicate emission
    t.close()


def test_select_orig_selects_without_emitting_and_unknown_clears():
    t = _table()
    got = []
    t.current_row_changed.connect(got.append)
    t.sortItems(1, Qt.AscendingOrder)
    assert t.select_orig(3) is True
    assert got == [] and t.current_orig() == 3
    assert t.item(t.currentRow(), 0).text() == "r3" and t.selectedItems()
    assert t.select_orig(99) is False and t.current_orig() == -1 and not t.selectedItems()
    assert t.select_orig(-1) is False and got == []
    t.close()


def test_selection_is_sticky_across_a_rebuild_and_a_capped_resort():
    t = _table()
    got = []
    t.select_orig(2)
    t.current_row_changed.connect(got.append)
    t.set_rows([[f"s{i}", i] for i in range(6)])                   # same shape rebuilt (a re-render)
    assert t.current_orig() == 2 and t.item(t.currentRow(), 0).text() == "s2" and got == []
    t.set_rows([["only", 1]])                                      # the row is gone: cleared, and said so once
    assert t.current_orig() == -1 and got == [-1]
    capped = widgets.DataTable(["id", "n"])
    rows = [[f"c{i}", i] for i in range(10)]
    capped.set_complete_rows(rows, 4)                              # shows c0..c3
    capped.select_orig(1)
    capped.sortItems(1, Qt.DescendingOrder)                        # now the visible prefix is c9..c6
    assert [capped.item(i, 0).text() for i in range(4)] == ["c9", "c8", "c7", "c6"]
    assert capped.current_orig() == -1                             # c1 left the visible prefix: cleared
    capped.sortItems(1, Qt.AscendingOrder)
    assert [capped._orig(i) for i in range(4)] == [0, 1, 2, 3]
    t.close()


def test_row_activated_stays_the_double_click_signal():
    t = _table()
    got = []
    t.row_activated.connect(got.append)
    t.sortItems(1, Qt.DescendingOrder)
    r = 2
    orig = t._orig(r)
    t.doubleClicked.emit(t.model().index(r, 0))                    # what a double-click raises
    assert got == [orig] and orig != r
    t.close()


def test_ctrl_c_copies_the_selected_rows_as_clean_values_not_one_cell_of_display_text():
    t = widgets.DataTable(["a", "dG"])
    t.set_rows([["x", "-9.4 !"], ["y", "-3.1"]], exports=[["x", "-9.4"], ["y", None]])
    t.resize(300, 200); t.show(); QTest.qWaitForWindowActive(t)
    t.setFocus(); t.selectRow(0)
    QTest.keyClick(t, Qt.Key_C, Qt.ControlModifier)
    assert _app.clipboard().text() == "x\t-9.4"
    t.selectAll()
    QTest.keyClick(t, Qt.Key_C, Qt.ControlModifier)
    assert _app.clipboard().text() == "x\t-9.4\ny\t-3.1"
    t.close()


def test_fit_rows_makes_the_table_exactly_as_tall_as_its_rows_and_keeps_it_that_way(scale):
    t = widgets.DataTable(["a", "b"])
    t.set_rows([[f"r{i}", i] for i in range(24)])
    holder = QWidget(); lay = QVBoxLayout(holder); lay.addWidget(t)
    holder.resize(400, 300); holder.show(); _app.processEvents()
    assert t.verticalScrollBar().maximum() > 0                     # unfitted: an inner scroll box
    t.fit_rows()
    _app.processEvents()
    assert t.height() == t._rows_height(24) and t.verticalScrollBar().maximum() == 0
    t.set_rows([[f"r{i}", i] for i in range(10)])                  # refits on a new row set
    assert t.height() == t._rows_height(10)
    t.fit_rows(4)                                                  # a capped view: four rows, scroll beyond
    t.set_rows([[f"r{i}", i] for i in range(10)])
    assert t.height() == t._rows_height(4) and t.verticalScrollBarPolicy() == Qt.ScrollBarAsNeeded
    t.fit_rows()
    theme.UI_SCALE = 1.5                                           # a live scale change refits on the next event
    _app.setStyleSheet(theme.qss("light"))
    for _ in range(5):
        _app.processEvents()
    assert t.height() == t._rows_height(10) and t.verticalScrollBar().maximum() == 0
    holder.close()


def test_wheel_over_a_fitted_table_is_left_for_the_page_but_an_inner_scroll_box_keeps_it():
    from PySide6.QtGui import QWheelEvent
    from PySide6.QtCore import QPointF, QPoint as QP

    def wheel(table):
        ev = QWheelEvent(QPointF(20, 20), QPointF(20, 20), QP(0, 0), QP(0, -120), Qt.NoButton, Qt.NoModifier,
                         Qt.NoScrollPhase, False)
        QApplication.sendEvent(table.viewport(), ev)
        return ev.isAccepted()
    boxed = widgets.DataTable(["a"])
    boxed.set_rows([[f"r{i}"] for i in range(40)])
    boxed.setMaximumHeight(120); boxed.show(); _app.processEvents()
    assert boxed.verticalScrollBar().maximum() > 0 and wheel(boxed) is True       # the capped box swallows the wheel
    fitted = widgets.DataTable(["a"])
    fitted.set_rows([[f"r{i}"] for i in range(40)])
    fitted.fit_rows(); fitted.show(); _app.processEvents()
    assert fitted.verticalScrollBar().maximum() == 0 and wheel(fitted) is False   # ignored: the page scrolls instead
    boxed.close(); fitted.close()


# ================================================================= widgets: export (CSV notes, '#', fallbacks, failures)
def test_csv_with_notes_round_trips_through_pandas_both_ways_and_the_csv_module(tmp_path):
    pd = pytest.importorskip("pandas")
    import csv, io
    t = widgets.DataTable(["#", "Feature", "Coords"])
    t.set_rows([[1, "LTR", "0–276"], [2, "PPT", "4861–4870"], [3, "TSD", "10–15"]])
    notes = ["Caveat one, with a comma, and a second, and a third, and more.",
             'A "quoted" caveat, with commas, too, and more',
             "A plain caveat with no separator"]
    assert t.export_headers()[0] == "No."                          # '#' would turn the header into a comment
    text = widgets.serialize_table(t.export_headers(), t.rows_data(), ",", notes)
    plain = pd.read_csv(io.StringIO(text))                         # a ragged row used to raise ParserError here
    assert list(plain.columns) == ["No.", "Feature", "Coords"] and plain.iloc[:3, 1].tolist() == ["LTR", "PPT", "TSD"]
    skipped = pd.read_csv(io.StringIO(text), comment="#")          # the header and first row used to be lost here
    assert list(skipped.columns) == ["No.", "Feature", "Coords"] and skipped.iloc[0, 1] == "LTR"
    assert max(len(r) for r in csv.reader(io.StringIO(text))) == 3  # no note is wider than the table
    notes_back = [r[0] for r in csv.reader(io.StringIO(text)) if r and r[0].startswith("# ")]
    assert notes_back == ["# " + n for n in notes]                 # every note comes back verbatim
    assert t.export_headers() == ["No.", "Feature", "Coords"]
    other = widgets.DataTable(["Pair", "#x"])
    assert other.export_headers() == ["Pair", "#x"]


def test_write_table_returns_the_path_it_wrote_and_export_table_says_when_excel_fell_back(tmp_path, monkeypatch):
    p = str(tmp_path / "t.csv")
    assert widgets.write_table(["a"], [[1]], p) == p
    monkeypatch.setattr(widgets, "_HAS_XLSX", True)
    monkeypatch.setattr(widgets, "_export_xlsx", lambda *a, **k: False)   # an install openpyxl cannot import
    got = widgets.write_table(["a"], [[1]], str(tmp_path / "t.xlsx"))
    assert got == str(tmp_path / "t.csv")
    sent = []
    class Owner(QWidget):
        def _banner(self, msg, level="error"):
            sent.append((level, msg))
    owner = Owner()
    out = str(tmp_path / "TEagle_ORFs.xlsx")
    monkeypatch.setattr(widgets.QFileDialog, "getSaveFileName", staticmethod(lambda *a, **k: (out, "")))
    widgets.export_table(["a"], [[1]], "base", owner, fmt="xlsx")
    assert os.path.exists(str(tmp_path / "TEagle_ORFs.csv")) and not os.path.exists(out)
    assert sent and sent[-1][0] == "warn" and "saved as CSV" in sent[-1][1] and "TEagle_ORFs.csv" in sent[-1][1]


def test_a_permission_error_on_any_exporter_is_reported_not_swallowed(tmp_path, monkeypatch):
    sent = []
    class Owner(QWidget):
        def _banner(self, msg, level="error"):
            sent.append((level, msg))
    owner = Owner()
    locked = str(tmp_path / "open_in_acrobat.pdf")
    monkeypatch.setattr(widgets.QFileDialog, "getSaveFileName", staticmethod(lambda *a, **k: (locked, "")))
    def boom(svg, path):
        raise PermissionError(13, "in use")
    monkeypatch.setattr(widgets, "render_pdf", boom)
    panel = widgets.FigurePanel(lambda bg: '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 10 10"/>', "f",
                                parent=owner)
    panel._export_pdf()                                            # used to raise out of the Qt slot
    assert sent and sent[-1][0] == "warn" and "open_in_acrobat.pdf" in sent[-1][1] and "another program" in sent[-1][1]
    sent.clear()
    gp = widgets.GenomePanel(figures.svg_genome, "g", parent=owner)
    gp.set_model(figures.gv_tracks_from_rec(_rec()))
    gp._export_pdf()
    assert sent and sent[-1][0] == "warn"
    sent.clear()
    dp = widgets.DotPanel(parent=owner)
    dp.set_matrix(_dot_matrix())
    dp._export_pdf()
    assert sent and sent[-1][0] == "warn"
    sent.clear()
    fa = str(tmp_path / "locked.fasta")
    monkeypatch.setattr(widgets.QFileDialog, "getSaveFileName", staticmethod(lambda *a, **k: (fa, "")))
    monkeypatch.setattr(widgets, "atomic_write", lambda *a, **k: (_ for _ in ()).throw(PermissionError(13, "ro")))
    widgets.save_fasta(">x\nACGT", "base", owner)
    assert sent and sent[-1][0] == "warn" and "locked.fasta" in sent[-1][1]


def test_successful_exports_confirm_by_name_and_fasta_gets_its_extension(tmp_path, monkeypatch):
    sent = []
    class Owner(QWidget):
        def _banner(self, msg, level="error"):
            sent.append((level, msg))
    owner = Owner()
    noext = str(tmp_path / "amplicons")
    monkeypatch.setattr(widgets.QFileDialog, "getSaveFileName", staticmethod(lambda *a, **k: (noext, "")))
    widgets.save_fasta(">x\nACGT", "base", owner)
    assert os.path.isfile(noext + ".fasta") and sent[-1] == ("success", "Saved amplicons.fasta.")
    keep = str(tmp_path / "keep.fa")
    monkeypatch.setattr(widgets.QFileDialog, "getSaveFileName", staticmethod(lambda *a, **k: (keep, "")))
    widgets.save_fasta(">x\nACGT", "base", owner)
    assert os.path.isfile(keep)                                    # an existing FASTA extension is respected
    svgp = str(tmp_path / "fig")
    monkeypatch.setattr(widgets.QFileDialog, "getSaveFileName", staticmethod(lambda *a, **k: (svgp, "")))
    panel = widgets.FigurePanel(lambda bg: '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 10 10"/>', "f",
                                parent=owner)
    panel._export_svg()
    assert os.path.isfile(svgp + ".svg") and sent[-1] == ("success", "Saved fig.svg.")


def test_xlsx_strings_that_start_with_a_sign_or_at_carry_no_visible_apostrophe(tmp_path):
    openpyxl = pytest.importorskip("openpyxl")
    assert widgets._xlsx_val("+/-") == "+/-" and widgets._xlsx_val("@ref") == "@ref" and widgets._xlsx_val("-1+x") == "-1+x"
    assert widgets._xlsx_val("=SUM(A1)") == "'=SUM(A1)"            # only '=' is a formula in openpyxl: still neutralised
    p = str(tmp_path / "t.xlsx")
    widgets._export_xlsx(["strand", "n"], [["+/-", 1], ["@x", 2], ["-1+x", 3]], p)
    ws = openpyxl.load_workbook(p).active
    assert [ws.cell(row=r, column=1).value for r in (2, 3, 4)] == ["+/-", "@x", "-1+x"]
    assert all(ws.cell(row=r, column=1).data_type == "s" for r in (2, 3, 4))


# ================================================================= widgets: genome viewer interaction
def test_dragging_the_genome_viewer_pans_by_the_distance_dragged_at_every_ui_scale(scale):
    for ui in (1.0, 1.5):
        theme.UI_SCALE = ui
        panel = widgets.GenomePanel(figures.svg_genome, "g")
        panel.resize(int(900 * ui), int(500 * ui))
        panel.set_model({"length": 100000, "tracks": [{"name": "t", "height": 20, "features": [
            {"start": 40000, "end": 41000, "color": "#0072B2", "label": "RT"}]}]})
        panel.show(); _app.processEvents()
        panel.zoom_at(50000, 0.1)                                  # a 10 kb window
        v0 = dict(panel.view)
        c = panel.canvas
        x0 = 400
        QTest.mousePress(c, Qt.LeftButton, pos=QPoint(x0, 50))
        QTest.mouseMove(c, QPoint(x0 + 100, 50))
        QTest.mouseRelease(c, Qt.LeftButton, pos=QPoint(x0 + 100, 50))
        span = v0["end"] - v0["start"]
        expected = 100 / c._sc() / c._plot_w() * span              # 100 widget px, in the figure's logical units
        assert (v0["start"] - panel.view["start"]) == pytest.approx(expected, rel=0.02), ui
        panel.close()


def test_export_palette_is_an_explicit_visible_choice_that_a_background_click_also_sets():
    panel = widgets.GenomePanel(figures.svg_genome, "g")
    panel.set_model(figures.gv_tracks_from_rec(_rec()))
    assert "PUBLICATION" in panel.exportPalBtn.text().upper() and panel._export_target()["for_export"] is True
    panel._th_btns["white"].click()                                # even the already-selected background
    assert "AS SHOWN" in panel.exportPalBtn.text().upper() and panel._export_target()["for_export"] is False
    panel.exportPalBtn.click()
    assert "PUBLICATION" in panel.exportPalBtn.text().upper() and panel._export_target()["for_export"] is True
    assert panel.exportPalBtn.toolTip() and panel.exportPalBtn.accessibleName().startswith("Export palette")
    dot = widgets.DotPanel()
    dot.set_matrix(_dot_matrix())
    assert dot._export_as_shown is False
    dot.exportPalBtn.click()
    assert dot._export_as_shown is True and "AS SHOWN" in dot.exportPalBtn.text().upper()
    panel.deleteLater(); dot.deleteLater()


def test_toolbar_glyph_buttons_have_spoken_names_and_tips():
    panels = [widgets.FigurePanel(lambda bg: '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 10 10"/>', "f",
                                  modes=("dark", "white", "uv", "mono"))]
    gp = widgets.GenomePanel(figures.svg_genome, "g")
    dp = widgets.DotPanel()
    dp.set_matrix(_dot_matrix())
    panels += [gp, dp]
    for panel in panels:
        for b in panel.findChildren(QPushButton):
            name = b.accessibleName() or b.text()
            assert sum(ch.isalpha() for ch in name) >= 2, (type(panel).__name__, b.text())
            assert b.toolTip(), (type(panel).__name__, b.text())
    names = {b.accessibleName() for b in panels[0].findChildren(QPushButton)}
    assert {"Zoom out", "Zoom in", "Fit the figure to the window", "Export the figure as SVG",
            "Dark figure background", "Light figure background", "UV figure background",
            "Mono figure background"} <= names
    for p in panels:
        p.deleteLater()


def test_canvases_take_keyboard_focus_draw_a_ring_and_open_band_actions_from_the_keyboard():
    shown = []
    panel = widgets.FigurePanel(
        lambda bg: '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 100 100"><rect width="100" height="100" fill="#fff"/></svg>',
        "f", hit_regions=[{"x0": 0, "y0": 0, "x1": 10, "y1": 10, "tip": "300 bp"},
                          {"x0": 20, "y0": 0, "x1": 30, "y1": 10, "tip": "500 bp"}],
        on_menu=lambda r: shown.append(r["tip"]) or [])
    c = panel.canvas
    assert c.focusPolicy() == Qt.StrongFocus and c.accessibleName() and "Page Down" in c.accessibleDescription()
    panel.resize(400, 400); panel.show(); QTest.qWaitForWindowActive(panel)
    c.setFocus()
    before = c.grab().toImage()
    tx0 = c.tx
    QTest.keyClick(c, Qt.Key_Right)
    assert c.tx != tx0 and c._user_view                            # arrows pan
    QTest.keyClick(c, Qt.Key_PageDown); QTest.keyClick(c, Qt.Key_PageDown)
    assert c._kb == 1
    QTest.keyClick(c, Qt.Key_Return)
    assert shown == ["500 bp"]                                     # the band stepped to, not the first one
    assert c.grab().toImage() != before                            # a focus ring / band marker is painted
    gp = widgets.GenomePanel(figures.svg_genome, "g")
    assert gp.canvas.focusPolicy() == Qt.StrongFocus and gp.canvas.accessibleName() and gp.canvas.accessibleDescription()
    panel.close(); gp.deleteLater()


def test_intended_flag_marks_the_designed_product_and_a_missed_target_is_neutral_not_off_target():
    g = [{**_amp(300, False, start=1), "intended": True}, {**_amp(500, False, start=900), "intended": False}]
    assert figures._call_split([g[0]], True)["site"] == 1          # the designed product that misses the target
    assert figures._call_split([g[1]], True)["additional"] == 1    # carries `intended`, so it is design-bound
    bands = figures._lane_bands(g, figures.GELPAL["white"], True)
    assert [b["call"] for b in bands] == ["site", "off"]
    assert "designed product, target not covered" in bands[0]["t"] and "1 additional product" in bands[1]["t"]


def test_the_chevron_url_is_quoted_so_a_user_data_path_with_spaces_still_parses(tmp_path, monkeypatch):
    from teagle_core import appdirs
    spaced = tmp_path / "Jane Doe" / "TEagle data"
    spaced.mkdir(parents=True)
    monkeypatch.setattr(appdirs, "user_data_dir", lambda: str(spaced))
    css = theme.qss("light")
    m = re.search(r'down-arrow \{ image: url\("([^"]+)"\)', css)
    assert m and " " in m.group(1) and os.path.isfile(m.group(1))
    theme.UI_SCALE = 1.0
    _app.setStyleSheet(css)
    cb = QComboBox(); cb.addItems(["a", "b"]); cb.resize(160, 34); cb.show(); _app.processEvents()
    img = cb.grab().toImage()
    assert sum(1 for x in range(img.width() - 22, img.width() - 4) for y in range(8, img.height() - 8)
               if img.pixelColor(x, y).lightness() < 150) >= 8
    cb.close(); _app.setStyleSheet("")


def test_genome_viewer_position_readout_uses_plain_digits():
    panel = widgets.GenomePanel(figures.svg_genome, "g")
    panel.set_model({"length": 123456, "tracks": []})
    assert panel.pos.text() == "0–123456 bp · 123.46 kb"
    panel.deleteLater()


# ================================================================= DataTable fit: re-entrancy and lifetime
def test_datatable_apply_fit_is_reentrancy_safe(monkeypatch):
    t = widgets.DataTable(["a"]); t.set_rows([["x"]]); t.fit_rows()
    calls = []
    real = widgets.QTableWidget.setFixedHeight
    def rec(self, h):
        calls.append(h)
        if len(calls) < 50:
            self._apply_fit()                                      # the scrollbar-range change re-enters
        real(self, h)
    monkeypatch.setattr(widgets.QTableWidget, "setFixedHeight", rec)
    t._apply_fit()
    assert len(calls) <= 2                                         # bounded: one pass plus one settle pass
    t.deleteLater()


def test_datatable_apply_fit_ignores_a_deleted_widget():
    import shiboken6
    t = widgets.DataTable(["a"]); t.fit_rows()
    shiboken6.delete(t)
    widgets.DataTable._apply_fit(t)                                # must not raise on a dead C++ object


def test_mixed_gel_legend_covers_every_lane_kind():
    scored = {"label": "P1", "has_locus": True, "amplicons": [_amp(300, True, True, start=1)]}
    neutral = {"label": "P2", "has_locus": False, "amplicons": [_amp(500, None, True, start=900)]}
    texts = [t for _, t in figures._gel_legend_items([scored, neutral])]
    assert "on-target" in texts and "designed product (target not assessed)" in texts
    scan = {"label": "S", "has_locus": False, "amplicons": [_amp(700, None, start=50)]}
    assert "priming site" in [t for _, t in figures._gel_legend_items([scored, scan])]
    legend = _gel_text(figures.svg_gel({"lanes": [scored, neutral]}, "white"))
    assert "designed product (target not assessed)" in legend and "on-target" in legend


def test_gel_legend_names_only_the_call_buckets_present():
    def label(*amps):
        lane = {"label": "P", "has_locus": True, "amplicons": list(amps)}
        return dict(figures._gel_legend_items([lane]))["off"]
    assert label(_amp(300, True, True, start=1)) == "off-target"          # no extra band: the swatch keeps its generic name (test_gel pins it)
    assert label(_amp(300, True, True, start=1), _amp(500, False, False, start=900)) == "additional product"
    both = label(_amp(300, True, True, start=1), _amp(500, False, False, start=900), _amp(800, False, start=2000))
    assert both == "additional product / off-target"
