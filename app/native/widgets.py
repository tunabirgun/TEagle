"""Reusable Qt widgets for the native app: an SVG figure panel (zoom/pan/WYSIWYG export),
an interactive genome viewer (windowed semantic zoom), and a data table with CSV/TSV export
and a copy context menu. All figure rendering goes through QSvgRenderer so on-screen output
matches the exported SVG/PNG (the gel's Gaussian-blur glow is the one SVG-Tiny casualty on screen)."""
from __future__ import annotations
import contextlib, functools, math, os, re

from PySide6.QtCore import Qt, QByteArray, QEvent, QPointF, QRectF, QSize, QTimer, Signal
from PySide6.QtGui import QImage, QPainter, QColor, QPixmap, QIcon, QKeySequence, QPen
from PySide6.QtSvg import QSvgRenderer
from PySide6.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QLabel, QSizePolicy,
                               QTableWidget, QTableWidgetItem, QMenu, QFileDialog, QApplication,
                               QAbstractItemView, QHeaderView, QToolTip, QProgressBar, QScrollArea,
                               QColorDialog)

import theme                       # UI_SCALE: figures are authored in logical units and stretched by it
from figures import gv_gutter, GV_MR, DOT_ML, DOT_MT, DOT_MR   # layout the hit-tests must share with the SVG

MODE_LABEL = {"dark": "dark", "white": "light", "uv": "UV", "mono": "mono", "transparent": "transparent"}

# per-cell export/copy value, when the displayed text carries a UI-only mark (see DataTable.set_rows)
EXPORT_ROLE = Qt.UserRole + 1


# Units and symbols that keep their written case inside an uppercased label: "PROD MIN (bp)", "OPT Tm (°C)",
# "ORFS (≥40 aa)", "LOCAL isPcr". The set IS the specification — uppercasing them is wrong notation.
KEEP_CASE = frozenset({"bp", "kb", "Mb", "nt", "aa", "Tm", "ΔG", "kcal", "mol", "mRNA", "cDNA", "tRNA",
                       "isPcr", "minimap2"})
_WORD = re.compile(r"[^\W_]+")
_NUM_UNIT = re.compile(r"(\d+)(\D.*)")


def caps(text: str) -> str:
    """The app's assay-terminal uppercase for labels, leaving KEEP_CASE tokens (and a number glued to one, "150bp")
    as written. Idempotent."""
    def word(m):
        w = m.group(0)
        if w in KEEP_CASE:
            return w
        d = _NUM_UNIT.fullmatch(w)
        return w if d and d.group(2) in KEEP_CASE else w.upper()
    return _WORD.sub(word, text)


def _is_action_button(b) -> bool:
    return not (b.property("link") or b.objectName() == "cardhdr")


def set_caps(button, text: str):
    """setText for a button whose label is re-set after construction: applies the same casing as
    uppercase_buttons (none for link buttons and card headers), so a state change never leaves it in sentence case."""
    button.setText(caps(text) if _is_action_button(button) else text)


def uppercase_buttons(root):
    """Uppercase every action button under `root` (the web '.btn' look), root itself included. Skips link buttons and
    card headers — their titles stay sentence-case — and leaves glyphs and KEEP_CASE units untouched. Idempotent.
    Lives here so a dialog built outside MainWindow (installer, notification, genome manager) can call it as it is
    built, instead of waiting for the next theme/scale/result refresh to reach it."""
    buttons = root.findChildren(QPushButton)
    if isinstance(root, QPushButton):
        buttons.append(root)
    for b in buttons:
        if not _is_action_button(b):
            continue
        t = b.text()
        if t and t != caps(t):
            b.setText(caps(t))


def _svg_size(svg: str):
    m = None
    import re
    mm = re.search(r'viewBox="0 0 ([\d.]+) ([\d.]+)"', svg)
    return (float(mm.group(1)), float(mm.group(2))) if mm else (800.0, 400.0)


@contextlib.contextmanager
def atomic_write(path: str, mode: str = "w", **kw):
    """Write a user deliverable so an interrupted run cannot leave a half-file wearing the final name.

    A plain open(path,'w') truncates the target the instant it is called, so a crash, a full disk, or a
    closed lid mid-write leaves a file that opens fine and is quietly incomplete — a truncated results
    table is exactly the kind of wrong data that survives into a figure unnoticed. Same temp-then-replace
    the assembly cache already uses (fetch.py:537-540). os.replace is atomic on the same filesystem, and
    the temp file sits beside the target so it never crosses one."""
    tmp = f"{path}.part"
    try:
        with open(tmp, mode, **kw) as fh:
            yield fh
        os.replace(tmp, path)
    except BaseException:
        try:
            os.remove(tmp)                      # never leave a stray .part beside the user's export
        except OSError:
            pass
        raise


def atomic_save(path: str, writer):
    """Same guarantee for writers that take a PATH rather than a handle (QImage.save, openpyxl, QPdfWriter).
    `writer(tmp_path)` must return falsey ONLY to signal failure — anything else counts as written."""
    tmp = f"{path}.part"
    try:
        ok = writer(tmp)
        if ok is False or not os.path.exists(tmp):
            raise OSError(f"writer produced no file for {os.path.basename(path)}")
        os.replace(tmp, path)
    except BaseException:
        try:
            os.remove(tmp)
        except OSError:
            pass
        raise


def save_svg(svg: str, path: str):
    def write_clean(tmp):
        with open(tmp, "w", encoding="utf-8", newline="\n") as f:
            f.write(svg)
        from teagle_core.artifact_metadata import clean_svg_metadata
        clean_svg_metadata(tmp, author=None)

    atomic_save(path, write_clean)


def render_png(svg: str, path: str, scale: int = 3):
    """Rasterise an SVG string to a transparent-background PNG at `scale`x for publication."""
    w, h = _svg_size(svg)
    r = QSvgRenderer(QByteArray(svg.encode("utf-8")))
    img = QImage(int(w * scale), int(h * scale), QImage.Format_ARGB32)
    img.fill(Qt.transparent)
    p = QPainter(img)
    r.render(p)
    p.end()
    def write_clean(tmp):
        if not img.save(tmp, "PNG"):
            raise OSError("Qt could not write the PNG export")
        from teagle_core.artifact_metadata import clean_png_metadata
        clean_png_metadata(tmp, author=None)

    atomic_save(path, write_clean)


def render_pdf(svg: str, path: str):
    """Write an SVG string to a VECTOR PDF (a journal-ready figure): SVG -> QPainter -> QPdfWriter, never
    via a raster QImage, so text stays text and the marks stay resolution-independent. The page is sized to
    the figure's own aspect ratio (px -> mm at 96 dpi) with no margin, so the whole plot fills one page."""
    from PySide6.QtGui import QPdfWriter, QPageSize, QPageLayout
    from PySide6.QtCore import QSizeF, QMarginsF, QRectF
    w, h = _svg_size(svg)
    def write_clean(tmp):
        writer = QPdfWriter(tmp)
        writer.setCreator("")
        writer.setTitle("")
        writer.setResolution(96)                                      # 1 SVG px == 1 device px at 96 dpi
        writer.setPageSize(QPageSize(QSizeF(w * 25.4 / 96.0, h * 25.4 / 96.0), QPageSize.Unit.Millimeter))
        writer.setPageMargins(QMarginsF(0, 0, 0, 0), QPageLayout.Unit.Millimeter)
        r = QSvgRenderer(QByteArray(svg.encode("utf-8")))
        p = QPainter(writer)
        r.render(p, QRectF(0, 0, writer.width(), writer.height()))   # fill the page in device pixels
        p.end()
        from teagle_core.artifact_metadata import clean_pdf_metadata
        clean_pdf_metadata(tmp, author=None)

    atomic_save(path, write_clean)


_BG_NAME = {"dark": "Dark", "white": "Light", "uv": "UV", "mono": "Mono", "transparent": "Transparent"}
_EXPORT_TIP = {"SVG": "Export the figure as SVG (vector; text stays editable text)",
               "PNG": "Export the figure as PNG (raster, 3x resolution)",
               "PDF": "Export the figure as a vector PDF (one page, text stays text)"}


def _name_button(b, name, tip=None):
    """Spoken name and hover tip for a toolbar button whose label is a glyph or an abbreviation."""
    b.setAccessibleName(name)
    b.setToolTip(tip or name)


def _with_ext(path, ext):
    return path if path.lower().endswith(ext) else path + ext


def _ring_color():
    return QColor(theme.PAL["light" if theme.ACTIVE == "light" else "dark"]["accent"])


def _focus_ring(widget, painter):
    """A 2px accent frame around a canvas that has keyboard focus — canvases draw their own content, so Qt paints
    no focus cue on them."""
    pen = QPen(_ring_color()); pen.setWidth(2)
    painter.setPen(pen)
    painter.setBrush(Qt.NoBrush)
    painter.drawRect(widget.rect().adjusted(1, 1, -1, -1))


class SvgCanvas(QWidget):
    """Paints an SVG string with a user scale + pan offset. Wheel = zoom at cursor, drag = pan.
    Keyboard (the canvas takes focus): arrows pan, + / - zoom, 0 or Home fits, Page Up / Page Down step through the
    hover regions (gel bands) and show their detail, Enter or the Menu key opens that region's actions."""
    def __init__(self, parent=None):
        super().__init__(parent)
        self._svg = ""
        self._renderer = None
        self._sw, self._sh = 800.0, 400.0
        self.scale = 1.0
        self.tx = 0.0
        self.ty = 0.0
        self._drag = None
        self._user_view = False       # user zoomed/panned -> stop auto-refitting on resize
        self.regions = []             # [{x0,y0,x1,y1 (svg coords), tip, ...}] for hover / right-click
        self._on_menu = None          # callable(region) -> list[(label, fn)]
        self._kb = -1                 # index of the region the keyboard has stepped to
        self.setMinimumHeight(300)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.setMouseTracking(True)
        self.setFocusPolicy(Qt.StrongFocus)
        self.setAccessibleName("Figure")
        self.setAccessibleDescription("Arrow keys pan, plus and minus zoom, 0 fits the figure. "
                                      "Page Up and Page Down step through the bands; Enter opens the actions for a band.")

    def set_hit_regions(self, regions, on_menu=None):
        self.regions = regions or []
        self._on_menu = on_menu
        self._kb = -1

    def _open_menu(self, region, gpos):
        items = self._on_menu(region) if self._on_menu else None
        if not items:
            return
        m = QMenu(self)
        for label, fn in items:
            m.addAction(label, fn)
        m.exec(gpos)

    def keyPressEvent(self, e):
        k = e.key()
        step = 40.0
        if k in (Qt.Key_Left, Qt.Key_Right, Qt.Key_Up, Qt.Key_Down):
            dx = {Qt.Key_Left: step, Qt.Key_Right: -step}.get(k, 0.0)
            dy = {Qt.Key_Up: step, Qt.Key_Down: -step}.get(k, 0.0)
            self.tx += dx; self.ty += dy; self._user_view = True; self.update()
        elif k in (Qt.Key_Plus, Qt.Key_Equal):
            self._zoom_at(self.width() / 2, self.height() / 2, 1.25)
        elif k in (Qt.Key_Minus, Qt.Key_Underscore):
            self._zoom_at(self.width() / 2, self.height() / 2, 0.8)
        elif k in (Qt.Key_0, Qt.Key_Home):
            self.fit()
        elif k in (Qt.Key_PageDown, Qt.Key_PageUp) and self.regions:
            self._kb = (self._kb + (1 if k == Qt.Key_PageDown else -1)) % len(self.regions)
            r = self.regions[self._kb]
            QToolTip.showText(self.mapToGlobal(self.rect().center()), r.get("tip", ""), self)
            self.update()
        elif k in (Qt.Key_Menu, Qt.Key_Return, Qt.Key_Enter) and self.regions and self._on_menu:
            r = self.regions[self._kb if self._kb >= 0 else 0]
            self._open_menu(r, self.mapToGlobal(self.rect().center()))
        else:
            super().keyPressEvent(e)

    def _svg_at(self, wx, wy):
        if self.scale <= 0:
            return (0.0, 0.0)
        return ((wx - self.tx) / self.scale, (wy - self.ty) / self.scale)

    def _region_at(self, wx, wy):
        sx, sy = self._svg_at(wx, wy)
        for r in self.regions:
            if r["x0"] <= sx <= r["x1"] and r["y0"] <= sy <= r["y1"]:
                return r
        return None

    def set_svg(self, svg: str, refit: bool = True):
        self._svg = svg
        self._renderer = QSvgRenderer(QByteArray(svg.encode("utf-8")))
        self._sw, self._sh = _svg_size(svg)
        if refit:
            self.fit()
        self.update()

    def fit(self):
        vw, vh = max(self.width(), 1), max(self.height(), 1)
        self.scale = min(vw / self._sw, vh / self._sh) * 0.94
        self.tx = (vw - self._sw * self.scale) / 2
        self.ty = (vh - self._sh * self.scale) / 2
        self._user_view = False       # back to auto-fit until the user zooms/pans again
        self.update()

    def resizeEvent(self, e):
        """Refit to the new viewport, like GenomePanel's canvas re-renders on resize; a user's own
        zoom/pan is left alone. fit() only repaints, so this cannot loop."""
        super().resizeEvent(e)
        if not self._user_view:
            self.fit()

    def _zoom_at(self, cx, cy, factor):
        self._user_view = True
        ns = min(16.0, max(0.08, self.scale * factor))
        self.tx = cx - (cx - self.tx) * (ns / self.scale)
        self.ty = cy - (cy - self.ty) * (ns / self.scale)
        self.scale = ns
        self.update()

    def wheelEvent(self, e):
        f = 1.12 if e.angleDelta().y() > 0 else 0.89
        self._zoom_at(e.position().x(), e.position().y(), f)

    def mousePressEvent(self, e):
        if e.button() != Qt.LeftButton:      # right-click is for the context menu, not panning
            return
        self._drag = (e.position().x(), e.position().y(), self.tx, self.ty)
        self.setCursor(Qt.ClosedHandCursor)

    def mouseMoveEvent(self, e):
        if self._drag:
            x0, y0, tx0, ty0 = self._drag
            self.tx = tx0 + (e.position().x() - x0)
            self.ty = ty0 + (e.position().y() - y0)
            self._user_view = True
            self.update()
            return
        r = self._region_at(e.position().x(), e.position().y())    # hover: show the feature detail
        if r:
            QToolTip.showText(e.globalPosition().toPoint(), r.get("tip", ""), self)
        else:
            QToolTip.hideText()

    def mouseReleaseEvent(self, e):
        self._drag = None
        self.unsetCursor()

    def contextMenuEvent(self, e):
        if not self._on_menu:
            return
        r = self._region_at(e.pos().x(), e.pos().y())
        if not r:
            return
        self._open_menu(r, e.globalPos())

    def paintEvent(self, _):
        p = QPainter(self)
        p.setRenderHint(QPainter.Antialiasing)
        p.setRenderHint(QPainter.SmoothPixmapTransform)
        if self._renderer:
            p.save()
            p.translate(self.tx, self.ty)
            p.scale(self.scale, self.scale)
            self._renderer.render(p, QRectF(0, 0, self._sw, self._sh))
            if self.hasFocus() and 0 <= self._kb < len(self.regions):       # the band the keyboard stepped to
                r = self.regions[self._kb]
                pen = QPen(_ring_color()); pen.setWidthF(2.0 / max(self.scale, 0.01))
                p.setPen(pen)
                p.drawRect(QRectF(r["x0"] - 2, r["y0"] - 2, r["x1"] - r["x0"] + 4, r["y1"] - r["y0"] + 4))
            p.restore()
        if self.hasFocus():
            _focus_ring(self, p)
        p.end()


class FigurePanel(QWidget):
    """Toolbar (bg modes + zoom + fit + export) over an SvgCanvas. `build_fn(bg)->svg` supplies
    the figure; export is WYSIWYG — it writes the currently selected background mode."""
    def __init__(self, build_fn, base_name: str, modes=("dark", "white"), parent=None,
                 hit_regions=None, on_menu=None, app_theme=None):
        super().__init__(parent)
        self.build_fn = build_fn
        self.base_name = base_name
        self.modes = list(modes)
        self.bg = self.modes[0]
        # open in the app's look: the given app theme, else the one the stylesheet was last built for
        first = {"light": "white", "dark": "dark"}.get(app_theme or theme.ACTIVE)
        if first in self.modes:
            self.bg = first
        self._theme_locked = False        # a manual bg pick pins this panel; app-theme toggles stop following it
        self._hit_regions = hit_regions
        self._on_menu = on_menu
        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(4)
        bar = QHBoxLayout()
        bar.addWidget(QLabel("Background"))
        self._mode_btns = {}
        for m in self.modes:
            b = QPushButton(caps(MODE_LABEL.get(m, m)))
            b.setProperty("sm", True)
            _name_button(b, f"{_BG_NAME.get(m, m)} figure background",
                         f"{_BG_NAME.get(m, m)} background for the figure and its exports")
            b.clicked.connect(lambda _=False, mm=m: self._set_bg(mm))
            bar.addWidget(b)
            self._mode_btns[m] = b
        bar.addStretch(1)
        for txt, fn, name, tip in (("−", self._zoom_out, "Zoom out", None), ("FIT", self._fit, "Fit the figure to the window", None),
                                   ("+", self._zoom_in, "Zoom in", None),
                                   ("SVG", self._export_svg, "Export the figure as SVG", _EXPORT_TIP["SVG"]),
                                   ("PNG", self._export_png, "Export the figure as PNG", _EXPORT_TIP["PNG"]),
                                   ("PDF", self._export_pdf, "Export the figure as vector PDF", _EXPORT_TIP["PDF"])):
            b = QPushButton(txt)
            b.setProperty("sm", True)
            _name_button(b, name, tip)
            b.clicked.connect(fn)
            bar.addWidget(b)
        lay.addLayout(bar)
        self.canvas = SvgCanvas()
        if self._hit_regions is not None or self._on_menu is not None:
            self.canvas.set_hit_regions(self._hit_regions, self._on_menu)
        lay.addWidget(self.canvas)
        self.render()

    def _set_bg(self, m):
        self._theme_locked = True         # user chose a bg explicitly -> keep it across later app-theme changes
        self.bg = m
        self.render()

    def apply_app_theme(self, app_theme):
        """Follow the app's dark/light theme (app 'light' -> figure bg 'white') UNTIL the user picks a bg
        manually, after which this panel keeps that choice (incl. uv/mono) and ignores app-theme toggles.
        Re-renders WITHOUT refitting, so the user's zoom/pan is preserved. No-op if the mapped mode isn't
        offered by this panel."""
        if self._theme_locked:            # user's manual bg pick wins over app-theme propagation
            return
        m = "white" if app_theme == "light" else "dark"
        if m not in self.modes or m == self.bg:
            return
        self.bg = m
        for mm, b in self._mode_btns.items():
            b.setProperty("primary", mm == self.bg); b.style().unpolish(b); b.style().polish(b)
        self.canvas.set_svg(self.build_fn(self.bg), refit=False)

    def render(self):
        for m, b in self._mode_btns.items():
            b.setProperty("primary", m == self.bg)
            b.style().unpolish(b); b.style().polish(b)
        self.canvas.set_svg(self.build_fn(self.bg))

    def _zoom_in(self):
        self.canvas._zoom_at(self.canvas.width() / 2, self.canvas.height() / 2, 1.25)

    def _zoom_out(self):
        self.canvas._zoom_at(self.canvas.width() / 2, self.canvas.height() / 2, 0.8)

    def _fit(self):
        self.canvas.fit()

    def _export_svg(self):
        path, _ = QFileDialog.getSaveFileName(self, "Export SVG", self.base_name + ".svg", "SVG (*.svg)")
        if path:
            path = _with_ext(path, ".svg")
            safe_export(self, path, lambda p: save_svg(self.build_fn(self.bg), p))   # export what you see (bg mode)

    def _export_png(self):
        path, _ = QFileDialog.getSaveFileName(self, "Export PNG", self.base_name + ".png", "PNG (*.png)")
        if path:
            path = _with_ext(path, ".png")
            safe_export(self, path, lambda p: render_png(self.build_fn(self.bg), p))

    def _export_pdf(self):
        path, _ = QFileDialog.getSaveFileName(self, "Export PDF", self.base_name + ".pdf", "PDF (*.pdf)")
        if path:
            path = _with_ext(path, ".pdf")
            safe_export(self, path, lambda p: render_pdf(self.build_fn(self.bg), p))   # same semantic SVG as SVG/PNG


class GenomePanel(QWidget):
    """Interactive genome viewer: wheel/buttons zoom the *bp window* (semantic zoom), drag pans it,
    export is WYSIWYG of the current window. Re-renders svg_genome each interaction, like the web viewer."""
    def __init__(self, svg_genome_fn, base_name="TEagle_genome", parent=None):
        super().__init__(parent)
        self._svg_genome = svg_genome_fn
        self.base_name = base_name
        self.model = {"length": 1, "tracks": []}
        self.theme = "dark"
        self._theme_locked = False        # a manual bg pick pins this viewer; app-theme toggles stop following it
        self._export_as_shown = False     # False: exports use the publication palette; True: the on-screen palette
        self.view = {"start": 0.0, "end": 1.0}
        self.on_feature_menu = None       # callable(region) -> list[(label, fn)] for right-click
        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(4)
        bar = QHBoxLayout()
        bar.addWidget(QLabel("Background"))
        self._th_btns = {}
        for th, lab in (("dark", "DARK"), ("white", "LIGHT")):
            b = QPushButton(lab); b.setProperty("sm", True)
            _name_button(b, f"{_BG_NAME[th]} figure background", f"{_BG_NAME[th]} background for the viewer")
            b.clicked.connect(lambda _=False, t=th: self._set_theme(t, user=True))
            bar.addWidget(b); self._th_btns[th] = b
        self.pos = QLabel(""); self.pos.setObjectName("gvpos")
        bar.addWidget(self.pos)
        bar.addStretch(1)
        for txt, fn, name in (("−", lambda: self._zoom(1.6), "Zoom out"), ("FIT", self._fit, "Fit the whole record"),
                              ("+", lambda: self._zoom(0.625), "Zoom in")):
            b = QPushButton(txt); b.setProperty("sm", True); _name_button(b, name); b.clicked.connect(fn); bar.addWidget(b)
        self.exportPalBtn = QPushButton(); self.exportPalBtn.setProperty("sm", True)
        self.exportPalBtn.clicked.connect(self._toggle_export_palette)
        bar.addWidget(self.exportPalBtn)
        self._sync_export_palette_btn()
        for txt, fn, name in (("SVG", self._export_svg, "Export the figure as SVG"),
                              ("PNG", self._export_png, "Export the figure as PNG"),
                              ("PDF", self._export_pdf, "Export the figure as vector PDF")):
            b = QPushButton(txt); b.setProperty("sm", True); _name_button(b, name, _EXPORT_TIP[txt])
            b.clicked.connect(fn); bar.addWidget(b)
        lay.addLayout(bar)
        self.canvas = _GenomeCanvas(self)
        self.canvas.setMinimumWidth(round(320 * max(theme.UI_SCALE, 0.1)))   # 320 LOGICAL units — the SVG authoring floor
        lay.addWidget(self.canvas)                         # (below 320 the bp<->pixel map would drift; obs: genome hit-test)

    def set_model(self, model: dict):
        self.model = model
        L = model.get("length", 1) or 1
        self.view = {"start": 0.0, "end": float(L)}
        self._render()

    def set_feature_menu(self, cb):
        """cb(region) -> list[(label, fn)] built on right-click over a feature glyph."""
        self.on_feature_menu = cb

    def _set_theme(self, t, user=False):
        if user:                          # user chose a bg explicitly -> keep it across later app-theme changes
            self._theme_locked = True
            self._export_as_shown = True  # ... and that choice is also what the files carry (shown on the palette button)
            self._sync_export_palette_btn()
        self.theme = t
        self._render()

    def _toggle_export_palette(self):
        self._export_as_shown = not self._export_as_shown
        self._sync_export_palette_btn()

    def _sync_export_palette_btn(self):
        shown = self._export_as_shown
        set_caps(self.exportPalBtn, "Palette: as shown" if shown else "Palette: publication")
        self.exportPalBtn.setAccessibleName("Export palette: " + ("as shown on screen" if shown else "publication"))
        self.exportPalBtn.setToolTip(
            "Palette of the exported SVG, PNG and PDF. 'As shown' writes the viewer exactly as it looks now. "
            "'Publication' writes dark ink on a transparent background, legible on white paper. Click to switch.")

    def apply_app_theme(self, app_theme):
        """Follow the app's dark/light theme (app 'light' -> viewer 'white') UNTIL the user picks a bg
        manually, after which this viewer keeps that choice and ignores app-theme toggles. Re-renders in
        place so the current pan/zoom window (self.view) is preserved."""
        if self._theme_locked:            # user's manual bg pick wins over app-theme propagation
            return
        self._set_theme("white" if app_theme == "light" else "dark")

    def _cur_svg(self, w, for_export=False, theme=None, legend=False):
        kw = {"legend": True} if legend else {}
        return self._svg_genome(self.model, {"start": self.view["start"], "end": self.view["end"]},
                                w, theme or self.theme, for_export, **kw)

    def _render(self):
        for th, b in self._th_btns.items():
            b.setProperty("primary", th == self.theme)
            b.style().unpolish(b); b.style().polish(b)
        # author the SVG in LOGICAL units (device pixels / UI scale); the canvas CSS-stretches it back by the
        # same factor, so the viewer's own text and glyphs grow with the global UI scale like the rest of the chrome
        z = max(theme.UI_SCALE, 0.1)
        w = max((self.canvas.width() or round(620 * z)) / z, 320)
        svg, regions = self._svg_genome(self.model, {"start": self.view["start"], "end": self.view["end"]},
                                        w, self.theme, False, True)
        self.canvas.set_svg(svg)
        self.canvas.regions = regions
        self.canvas._kb = -1
        L = self.model.get("length", 1) or 1
        self.pos.setText(f"{int(self.view['start'])}–{int(self.view['end'])} bp · "
                         f"{(self.view['end']-self.view['start'])/1000:.2f} kb")

    def _clamp(self):
        L = self.model.get("length", 1) or 1
        sp = min(max(self.view["end"] - self.view["start"], 20), L)
        if sp >= L:
            self.view = {"start": 0.0, "end": float(L)}
            return
        st = max(0.0, min(self.view["start"], L - sp))
        self.view = {"start": st, "end": st + sp}

    def zoom_at(self, bp, factor):
        sp = (self.view["end"] - self.view["start"]) * factor
        frac = (bp - self.view["start"]) / max(self.view["end"] - self.view["start"], 1e-9)
        self.view = {"start": bp - frac * sp, "end": bp + (1 - frac) * sp}
        self._clamp(); self._render()

    def _zoom(self, factor):
        mid = (self.view["start"] + self.view["end"]) / 2
        self.zoom_at(mid, factor)

    def pan_bp(self, dbp):
        self.view = {"start": self.view["start"] + dbp, "end": self.view["end"] + dbp}
        self._clamp(); self._render()

    def _fit(self):
        L = self.model.get("length", 1) or 1
        self.view = {"start": 0.0, "end": float(L)}
        self._render()

    def _export_target(self):
        """Which palette a file export should use, and that every export carries its legend.

        The figure leaves the app to become a journal figure, so the DEFAULT export is the publication
        palette — dark ink, legible on white paper — regardless of the on-screen background, which exists
        for reading the viewer rather than for print. Choosing a background, or the 'as shown' palette
        button, is a deliberate instruction and is honoured; the palette button always says which rule applies."""
        return {"for_export": not self._export_as_shown, "theme": self.theme, "legend": True}

    def _export_svg(self):
        path, _ = QFileDialog.getSaveFileName(self, "Export SVG", self.base_name + ".svg", "SVG (*.svg)")
        if path:
            path = _with_ext(path, ".svg")
            safe_export(self, path, lambda p: save_svg(self._cur_svg(920, **self._export_target()), p))

    def _export_png(self):
        path, _ = QFileDialog.getSaveFileName(self, "Export PNG", self.base_name + ".png", "PNG (*.png)")
        if path:
            path = _with_ext(path, ".png")
            safe_export(self, path, lambda p: render_png(self._cur_svg(920, **self._export_target()), p))

    def _export_pdf(self):
        path, _ = QFileDialog.getSaveFileName(self, "Export PDF", self.base_name + ".pdf", "PDF (*.pdf)")
        if path:
            path = _with_ext(path, ".pdf")
            safe_export(self, path, lambda p: render_pdf(self._cur_svg(920, **self._export_target()), p))


class _GenomeCanvas(QWidget):
    MR = GV_MR

    @property
    def ML(self):
        """Left gutter of the SVG being painted — svg_genome sizes it to the widest track name, so the
        bp<->pixel map has to read the same number or wheel-zoom/pan would anchor off by the difference."""
        return gv_gutter(getattr(self.panel, "model", None), self._w or 620.0)

    def __init__(self, panel: GenomePanel):
        super().__init__(panel)
        self.panel = panel
        self._svg = ""
        self._renderer = None
        self._w = 620.0
        self._drag = None
        self.regions = []
        self.setMinimumHeight(220)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)
        self.setFocusPolicy(Qt.StrongFocus)
        self.setMouseTracking(True)
        self._kb = -1                 # feature the keyboard stepped to
        self.setAccessibleName("Genome viewer")
        self.setAccessibleDescription("Left and Right arrows pan, plus and minus zoom, Home or 0 shows the whole record. "
                                      "Page Up and Page Down step through the features; Enter opens the actions for one.")

    def _sc(self):
        """CSS-stretch factor: the SVG is authored in logical units and painted across the widget width,
        so pixel coords divide by this to reach SVG coords (and the figure scales with the UI scale)."""
        return ((self.width() or self._w) / self._w) if self._w else 1.0

    def _region_at(self, wx, wy):
        sc = self._sc()
        if sc <= 0:
            return None
        sx, sy = wx / sc, wy / sc
        for r in self.regions:
            if r["x0"] <= sx <= r["x1"] and r["y0"] <= sy <= r["y1"]:
                return r
        return None

    def _open_menu(self, region, gpos):
        cb = getattr(self.panel, "on_feature_menu", None)
        items = cb(region) if cb else None
        if not items:
            return
        m = QMenu(self)
        for label, fn in items:
            m.addAction(label, fn)
        m.exec(gpos)

    def set_svg(self, svg):
        self._svg = svg
        self._renderer = QSvgRenderer(QByteArray(svg.encode("utf-8")))
        self._w, self._h = _svg_size(svg)
        self.setMinimumHeight(int(self._h * self._sc()))      # painted height = SVG height x the CSS stretch
        self.update()

    def _plot_w(self):
        return max(120.0, (self._w or 620) - self.ML - self.MR)     # logical (SVG) units

    def _x_to_bp(self, px):
        frac = max(0.0, min(1.0, (px / self._sc() - self.ML) / self._plot_w()))
        v = self.panel.view
        return v["start"] + frac * (v["end"] - v["start"])

    def resizeEvent(self, _):
        self.panel._render()

    def wheelEvent(self, e):
        v = self.panel.view
        unit = 0.002
        factor = 2 ** (-e.angleDelta().y() * unit)
        self.panel.zoom_at(self._x_to_bp(e.position().x()), factor)

    def mousePressEvent(self, e):
        if e.button() != Qt.LeftButton:      # right-click opens the feature menu, not a pan
            return
        v = self.panel.view
        self._drag = (e.position().x(), v["start"], v["end"])
        self.setCursor(Qt.ClosedHandCursor)

    def mouseMoveEvent(self, e):
        if self._drag:
            x0, s0, e0 = self._drag
            dbp = (e.position().x() - x0) / self._sc() / self._plot_w() * (e0 - s0)   # px -> logical units, like _x_to_bp
            self.panel.view = {"start": s0 - dbp, "end": e0 - dbp}
            self.panel._clamp(); self.panel._render()
            return
        r = self._region_at(e.position().x(), e.position().y())     # hover: feature detail tooltip
        if r:
            QToolTip.showText(e.globalPosition().toPoint(), r.get("tip", ""), self)
        else:
            QToolTip.hideText()

    def mouseReleaseEvent(self, e):
        self._drag = None
        self.unsetCursor()

    def contextMenuEvent(self, e):
        r = self._region_at(e.pos().x(), e.pos().y())
        if r:
            self._open_menu(r, e.globalPos())

    def keyPressEvent(self, e):
        v = self.panel.view
        sp = v["end"] - v["start"]
        k = e.key()
        if k in (Qt.Key_PageDown, Qt.Key_PageUp) and self.regions:
            self._kb = (self._kb + (1 if k == Qt.Key_PageDown else -1)) % len(self.regions)
            r = self.regions[self._kb]
            QToolTip.showText(self.mapToGlobal(self.rect().center()), r.get("tip", ""), self)
            self.update()
        elif k in (Qt.Key_Menu, Qt.Key_Return, Qt.Key_Enter) and self.regions:
            self._open_menu(self.regions[self._kb if 0 <= self._kb < len(self.regions) else 0],
                            self.mapToGlobal(self.rect().center()))
        elif k == Qt.Key_Left:
            self.panel.pan_bp(-sp * 0.15)
        elif k == Qt.Key_Right:
            self.panel.pan_bp(sp * 0.15)
        elif k in (Qt.Key_Up, Qt.Key_Plus, Qt.Key_Equal):
            self.panel._zoom(0.625)
        elif k in (Qt.Key_Down, Qt.Key_Minus, Qt.Key_Underscore):
            self.panel._zoom(1.6)
        elif k in (Qt.Key_Home, Qt.Key_0):
            self.panel._fit()
        else:
            super().keyPressEvent(e)

    def paintEvent(self, _):
        p = QPainter(self)
        p.setRenderHint(QPainter.Antialiasing)
        if self._renderer:
            # SVG is authored at self._w; CSS-stretch to widget width (height scales to keep aspect)
            sc = self._sc()
            self._renderer.render(p, QRectF(0, 0, self._w * sc, self._h * sc))
            if self.hasFocus() and 0 <= self._kb < len(self.regions):          # the feature the keyboard stepped to
                r = self.regions[self._kb]
                pen = QPen(_ring_color()); pen.setWidth(2)
                p.setPen(pen)
                p.drawRect(QRectF(r["x0"] * sc - 2, r["y0"] * sc - 2, (r["x1"] - r["x0"]) * sc + 4,
                                  (r["y1"] - r["y0"]) * sc + 4))
        if self.hasFocus():
            _focus_ring(self, p)
        p.end()


# ---------- data tables ----------
def _is_number(s):
    try:
        float(s)
        return True
    except (TypeError, ValueError):
        return False


def _csv_escape(v, sep):
    v = "" if v is None else str(v)
    # Neutralise spreadsheet formula injection (CWE-1236). A NUMBER is exempt: "-9.4" is not a formula,
    # and quoting it turned every negative delta-G in the primer-QC tables into a text cell, so a CSV
    # loaded into R or pandas came back as a character column. Only a non-numeric +/- lead is escaped
    # (e.g. "-1+cmd"), which still leaves a bare strand marker intact.
    if v[:1] in ("=", "@", "\t", "\r") or (v[:1] in ("+", "-") and len(v) > 1 and not _is_number(v)):
        v = "'" + v
    if sep in v or '"' in v or "\n" in v:
        v = '"' + v.replace('"', '""') + '"'
    return v


# Excel export is optional: degrade to CSV/TSV if openpyxl is absent OR broken.
#
# The import is DEFERRED. openpyxl pulls numpy and costs ~245 ms — a quarter of a second added to every
# launch, for a feature used only when someone exports a table. `find_spec` locates the package without
# executing it, which is what makes the saving possible.
#
# But locating is a WEAKER test than importing: a present-but-broken install (a shadowing openpyxl.py, a
# truncated package, a submodule that raises) passes find_spec and then fails on the first export. The old
# eager `try: import ... except: _HAS_XLSX = False` caught that at startup and simply hid the Excel option.
# To keep that behaviour, `_xl()` returns None instead of raising when the real import fails, and clears
# `_HAS_XLSX` so every later menu, format list and writer sees the same answer the eager version gave.
import importlib.util as _ilu

_HAS_XLSX = _ilu.find_spec("openpyxl") is not None
_XL = {}


def _xl():
    """Import openpyxl on first use; return None (and disable XLSX) if it cannot be imported."""
    global _HAS_XLSX
    if not _XL:
        try:
            from openpyxl import Workbook
            from openpyxl.styles import Font
        except Exception:
            _HAS_XLSX = False                                # broken install: behave exactly as if absent
            return None
        _XL["Workbook"], _XL["Font"] = Workbook, Font
    return _XL


def _xlsx_val(v):
    """Write a real number when the whole cell is numeric (so Excel sorts/filters it), else guarded
    text. A cell that starts with '=' is neutralised (CWE-1236); openpyxl writes every other string, '+/-' and
    '@ref' included, as a string cell, so those carry no visible apostrophe."""
    s = "" if v is None else str(v)
    t = s.strip()
    if t:
        try:
            if t.lstrip("-").isdigit():
                return int(t)
            f = float(t)
            if math.isfinite(f):                     # never write nan/inf (Excel renders them as an error cell)
                return f
        except ValueError:
            pass
    if s[:1] == "=":                                 # only '=' makes openpyxl store a formula; +, - and @ stay text
        return "'" + s
    return s


def _export_xlsx(headers, rows, path, notes=None):
    """Write the table as a workbook. Returns False when openpyxl cannot be imported, so the caller can
    fall back — a broken install must cost the user their chosen FORMAT, never their data or a traceback."""
    _x = _xl()                                           # first XLSX export pays the import, later ones do not
    if _x is None:
        return False                                     # unusable install — write_table falls back to CSV
    wb = _x["Workbook"]()
    ws = wb.active
    ws.title = "TEagle"
    ws.append([str(h) for h in headers])
    for cell in ws[1]:
        cell.font = _x["Font"](bold=True)
    for r in rows:
        ws.append([_xlsx_val(c) for c in r])
    ws.freeze_panes = "A2"                                # keep the header visible while scrolling
    for n in (notes or []):                               # caveats below the data, after one blank row
        ws.append([])
        for ln in str(n).splitlines():
            if ln.strip():
                ws.append([ln])
    def write_clean(tmp):
        wb.save(tmp)
        from teagle_core.artifact_metadata import clean_ooxml_metadata
        clean_ooxml_metadata(tmp, author=None)

    atomic_save(path, write_clean)


_TABLE_FORMATS = [("Excel workbook (.xlsx)", "xlsx"),
                  ("CSV — comma-separated (.csv)", "csv"),
                  ("TSV — tab-separated (.tsv)", "tsv")]
_FMT_FILTER = {"xlsx": "Excel workbook (*.xlsx)", "csv": "CSV (*.csv)", "tsv": "TSV (*.tsv)"}


def _available_formats():
    return [(lbl, fmt) for lbl, fmt in _TABLE_FORMATS if fmt != "xlsx" or _HAS_XLSX]


def pick_table_format(parent, global_pos):
    """Pop a format menu (Excel/CSV/TSV) at global_pos; return 'xlsx'|'csv'|'tsv' or None if dismissed."""
    m = QMenu(parent)
    amap = {m.addAction(label): fmt for label, fmt in _available_formats()}
    return amap.get(m.exec(global_pos))


def export_table(headers, rows, base, parent=None, fmt=None, notes=None):
    """Write a table to a user-chosen file. `fmt` in {'xlsx','csv','tsv'} pre-selects the format so the
    save dialog offers exactly that type; fmt=None falls back to a multi-filter dialog."""
    if fmt is None:
        filters = (["Excel (*.xlsx)"] if _HAS_XLSX else []) + ["CSV (*.csv)", "TSV (*.tsv)"]
        path, _sel = QFileDialog.getSaveFileName(parent, "Export table",
                                                 base + (".xlsx" if _HAS_XLSX else ".csv"), ";;".join(filters))
    else:
        if fmt == "xlsx" and not _HAS_XLSX:
            fmt = "csv"
        ext = "." + fmt
        path, _sel = QFileDialog.getSaveFileName(parent, f"Export table as {fmt.upper()}",
                                                 base + ext, _FMT_FILTER[fmt])
        if path and not path.lower().endswith(ext):
            path += ext                                       # honor the chosen format if the user omits the extension
    if not path:
        return
    # A failed export must say what to do, not raise. atomic_write correctly cleans up and re-raises, and
    # an unhandled OSError here reached the generic crash banner as a bare exception string — the exact
    # wrong primary message on the most ordinary failure there is (a read-only folder, a
    # disconnected drive, the file open in Excel).
    written = safe_export(parent, path, lambda p: write_table(headers, rows, p, notes), notify=False)
    if written is None:
        return
    name = os.path.basename(written)
    if os.path.normcase(written) != os.path.normcase(path):
        _notify(parent, f"Excel export is unavailable on this computer, so the table was saved as CSV instead: {name}.",
                "warn")
    else:
        _notify(parent, f"Saved {name}.", "success")


def _banner_owner(parent):
    """The nearest object that has a `_banner` (MainWindow), looking at the parent and then its window."""
    for cand in (parent, parent.window() if hasattr(parent, "window") else None):
        if callable(getattr(cand, "_banner", None)):
            return cand
    return None


def _notify(parent, msg, level="success"):
    owner = _banner_owner(parent)
    if owner is not None:
        owner._banner(msg, level)
    return owner is not None


def _export_failed(parent, msg):
    """Report an export failure through the owner's banner when there is one, else a message box. The table
    widgets are reused by dialogs that are not MainWindow, so this cannot assume _banner exists."""
    if _notify(parent, msg, "warn"):
        return
    from PySide6.QtWidgets import QMessageBox
    box = QMessageBox(parent)
    box.setWindowTitle("Export failed")
    box.setIcon(QMessageBox.Warning)
    box.setText(msg)
    box.exec()


def safe_export(parent, path, writer, ok=None, notify=True):
    """Run `writer(path)` for a user-chosen file and say what happened. A PermissionError (the file open in Acrobat
    or Excel, a read-only folder) or other OSError/ValueError is reported through the owner's banner instead of
    escaping a Qt slot — where it prints to a stderr nobody sees and leaves the old file in place with no message.
    Returns the path actually written (`writer` may return a different one, e.g. an Excel fallback to CSV) or None
    on failure. With notify, success is confirmed by name; `ok` overrides the message. Public: every export path
    in the app routes through this."""
    name = os.path.basename(path)
    try:
        out = writer(path)
    except PermissionError:
        _export_failed(parent, f"“{name}” could not be written — it may be open in another program "
                               "(Excel and PDF viewers lock a file while it is open), or the folder may be read-only. "
                               "Close it or choose another folder, then export again.")
        return None
    except (OSError, ValueError) as e:
        _export_failed(parent, f"“{name}” could not be written. If the folder is on a network or removable "
                               f"drive, check it is still connected, or choose a folder on this computer. "
                               f"({type(e).__name__}: {e})")
        return None
    out = out if isinstance(out, str) else path
    if notify:
        _notify(parent, ok or f"Saved {os.path.basename(out)}.", "success")
    return out


def serialize_table(headers, rows, sep=",", notes=None) -> str:
    """The exact text written for a CSV/TSV export. Separated from the file dialog so the contract —
    every displayed row present, in order, escaped so the delimiter survives a round trip — is testable
    without a GUI. An export that silently drops a column or mangles a value is a defect, and a defect
    that only a human reading a spreadsheet can catch is one that ships."""
    lines = [sep.join(_csv_escape(h, sep) for h in headers)]
    lines += [sep.join(_csv_escape(c, sep) for c in r) for r in rows]
    # The scientific qualifications shown beside the table travel with it. An export that
    # silently drops a hedge or a units label is a defect — and these tables become figures and supplementary
    # files in papers, where the screen the numbers came from is long gone. Each note is ONE field, so a plain
    # reader never meets a ragged row: it is written as "# text" and quoted only when the text holds the
    # separator or a quote (a comment-skipping reader, pandas comment='#' or R comment.char, then drops every
    # note that needed no quoting).
    for n in (notes or []):
        for ln in str(n).splitlines():
            if ln.strip():
                t = "# " + ln
                lines.append('"' + t.replace('"', '""') + '"' if (sep in t or '"' in t) else t)
    return "\r\n".join(lines)


def write_table(headers, rows, path, notes=None):
    """Write a table to `path`, choosing the format from its extension. `notes` are the on-screen caveats.
    Returns the path actually written: an .xlsx request whose openpyxl will not import is written as CSV beside
    it, and the caller must be able to say so."""
    low = str(path).lower()
    if low.endswith(".xlsx") and _HAS_XLSX:
        if _export_xlsx(headers, rows, path, notes) is not False:
            return path
        # openpyxl was locatable but would not import: write the same data as CSV under the chosen name's
        # stem rather than failing. _xl() has already cleared _HAS_XLSX, so the option disappears next time.
        path = path[:-5] + ".csv"
        low = path.lower()
    sep = "\t" if low.endswith(".tsv") else ","
    with atomic_write(path, "w", encoding="utf-8-sig", newline="") as f:   # BOM so Excel reads UTF-8
        f.write(serialize_table(headers, rows, sep, notes))
    return path


def save_fasta(fasta: str, base: str, parent=None):
    """Write a ready-made FASTA string to a user-chosen .fasta file. A name without a FASTA extension gets
    .fasta; a write failure is reported, not raised."""
    path, _ = QFileDialog.getSaveFileName(parent, "Export FASTA", base + ".fasta",
                                          "FASTA (*.fasta *.fa *.fna *.txt)")
    if not path:
        return
    if not path.lower().endswith((".fasta", ".fa", ".fna", ".txt")):
        path += ".fasta"
    if not fasta.endswith("\n"):
        fasta += "\n"
    def write(p):
        with atomic_write(p, "w", encoding="utf-8", newline="\n") as f:
            f.write(fasta)
    return safe_export(parent, path, write)


_NUM_RE = re.compile(r'-?\d+\.?\d*(?:[eE][-+]?\d+)?')
_GROUP_COMMA = re.compile(r'(?<=\d),(?=\d{3}(?!\d))')    # a thousands separator: digit , exactly three digits

def _sortkey(text: str):
    """Numeric key for a cell: a slash-composite like '2/1' or '60.1/58.3' (Mism/Tm/GC F/R) sorts by
    the sum of both values; otherwise the leading number (handles 1.2e-30, 45%, '0–5146' → 0). None
    when the cell has no number, so those columns fall back to case-insensitive string order."""
    if not text:
        return None
    text = _GROUP_COMMA.sub("", text)                    # '1,500' is fifteen hundred, not 1 then 500
    nums = _NUM_RE.findall(text)
    if not nums:
        return None
    try:
        if "/" in text and len(nums) >= 2:               # F/R composite → combine both, not just forward
            return float(nums[0]) + float(nums[1])
        return float(nums[0])
    except ValueError:
        return None


class BusyBar(QWidget):
    """Indeterminate progress bar + caption for a long operation with no parseable progress
    (RepeatMasker / minimap2 / isPcr / genome download). The bar animates on its own (Qt-driven), so the
    panel visibly stays alive without a hang-vs-working ambiguity; no manual QTimer to leak. Call set_text()
    to update the caption (e.g. from the genome-download log poll)."""
    def __init__(self, text="working…", parent=None):
        super().__init__(parent)
        lay = QVBoxLayout(self); lay.setContentsMargins(0, 0, 0, 0); lay.setSpacing(4)
        self.caption = QLabel(text); self.caption.setObjectName("orient"); self.caption.setWordWrap(True)
        self.bar = QProgressBar(); self.bar.setRange(0, 0)         # 0..0 = indeterminate 'busy' animation
        self.bar.setTextVisible(False); self.bar.setFixedHeight(6)
        lay.addWidget(self.caption); lay.addWidget(self.bar)

    def set_text(self, text):
        self.caption.setText(text)


class _Cell(QTableWidgetItem):
    """Table cell that sorts numerically when both cells parse as numbers, else alphabetically —
    so Score / E-value / aa / divergence / coords sort by value, not by string."""
    def __lt__(self, other):
        a, b = _sortkey(self.text()), _sortkey(other.text())
        if a is not None and b is not None:
            return a < b
        return self.text().casefold() < other.text().casefold()


# A cell that is ENTIRELY one number: optional sign, thousands separators, decimal, exponent, and an
# optional trailing % or unit-less suffix-free tail. Deliberately stricter than _NUM_RE (which matches a
# number ANYWHERE) because alignment is decided by what the whole cell is, not by what it contains:
# "0–276 · 4870–5146" holds numbers but is a coordinate composite and must not right-align away from its
# neighbours, while "2.0e-11" and "100.0%" must.
_PURE_NUM_RE = re.compile(r'^[+-]?\d{1,3}(?:,\d{3})*(?:\.\d+)?(?:[eE][+-]?\d+)?%?$'
                          r'|^[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?%?$')


def _numeric_columns(rows, ncols):
    """Which columns hold nothing but numbers, so digits can be made to line up.

    DERIVED from the data on every set_rows rather than declared per call site: a column list would go
    stale the moment a table gained a column, and would be wrong for the tables whose columns change
    shape with the result (the Metric column is '100.0%' for an LTR and 'GAGGGGGCG' for a PPT, so it is
    correctly NOT numeric). A column with no non-empty cell is not numeric — an all-blank column has no
    digits to align."""
    out = []
    for j in range(ncols):
        seen = False
        numeric = True
        for r in rows:
            if j >= len(r):
                continue
            t = "" if r[j] is None else str(r[j]).strip()
            if not t:
                continue
            seen = True
            if not _PURE_NUM_RE.match(t):
                numeric = False
                break
        out.append(seen and numeric)
    return out


class DataTable(QTableWidget):
    """A read-only table with CSV/TSV export and a copy-cell/row context menu. Optional per-row
    activation callback (double-click / Enter) and a right-click menu builder for FASTA-style actions.

    Selection contract: `current_row_changed(orig)` fires when the user moves the current row (mouse or keyboard),
    carrying the row's ORIGINAL index (its position in set_rows, whatever the sort) or -1 when there is none.
    `select_orig(orig)` selects that row from code WITHOUT emitting; `current_orig()` reads it back. A rebuild
    (set_rows / set_complete_rows / a re-sort of a capped view) re-selects the same original row silently; when
    that row is no longer shown (a capped view re-sorted past its limit) the selection is cleared and
    current_row_changed(-1) fires once. `row_activated(orig)` stays the double-click signal."""
    row_activated = Signal(int)
    current_row_changed = Signal(int)

    def __init__(self, headers, tooltips=None, parent=None):
        super().__init__(0, len(headers), parent)
        self.setHorizontalHeaderLabels(headers)
        self._headers = headers
        if tooltips:
            for i, h in enumerate(headers):
                if h in tooltips:
                    self.horizontalHeaderItem(i).setToolTip(tooltips[h])
        self.verticalHeader().setVisible(False)
        self.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.setAlternatingRowColors(True)
        self.setContextMenuPolicy(Qt.CustomContextMenu)
        self.customContextMenuRequested.connect(self._menu)
        self.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.horizontalHeader().setStretchLastSection(True)
        self.horizontalHeader().setMaximumSectionSize(360)            # cap wide free-text cols so the last column stays visible
        self.setTextElideMode(Qt.ElideRight)                          # elide overflow; full text is in the cell tooltip
        self.horizontalHeader().setDefaultAlignment(Qt.AlignCenter)   # until set_rows sees the data, headers centre
        self.setSortingEnabled(True)                                  # click a header to sort (numeric-aware via _Cell)
        self.doubleClicked.connect(lambda idx: self.row_activated.emit(self._orig(idx.row())))
        self._row_menu = None            # callable(orig_row_index)->list[(label, fn)]
        self.export_base = "TEagle_table"   # proposed export filename; the visible Export button sets the specific one
        self._complete_rows = None
        self._complete_limit = None
        self._complete_refresh = False
        self.horizontalHeader().sortIndicatorChanged.connect(self._complete_sort_changed)
        self._quiet = 0                  # >0 while code (not the user) moves the current row: no current_row_changed
        self._last_orig = -1             # original index last reported / restored
        self._fit_limit = False          # False = height not fitted; None = fit every row; int = fit that many rows
        self.currentCellChanged.connect(self._on_current_cell)
        self._fitting = self._refit = False
        self.horizontalScrollBar().rangeChanged.connect(self._on_hrange)

    def _rows_height(self, n):
        """Height of the frame, the header and the first n rows, plus the sideways scrollbar when one is needed."""
        hh = max(self.horizontalHeader().height(), self.horizontalHeader().sizeHint().height())
        h = 2 * self.frameWidth() + hh + sum(self.rowHeight(i) for i in range(min(n, self.rowCount())))
        hb = self.horizontalScrollBar()
        if hb.isVisible() or hb.maximum() > 0:            # a sideways-scrolling table must keep room for its bar
            h += hb.sizeHint().height()
        return h

    def sizeHint(self):
        """Height the CURRENT rows need — QAbstractScrollArea's stock 192px hint reserved 80–135px of empty
        grid under a short result. Callers still cap tall tables with setMaximumHeight()."""
        s = super().sizeHint()
        return QSize(s.width(), max(self._rows_height(self.rowCount()), self.minimumSizeHint().height()))

    def fit_rows(self, limit=None):
        """Make the table exactly as tall as its header plus its rows (all of them, or the first `limit`), with no
        inner scrollbar when every row is shown, so the page scrolls instead of a small box capturing the wheel.
        Stays fitted across set_rows, a UI-scale / style change and a sideways-scrollbar change."""
        self._fit_limit = limit
        self._apply_fit()

    def _on_hrange(self, *_):
        self._apply_fit()

    def _apply_fit(self):
        if self._fit_limit is False:
            return
        if self._fitting:                              # setFixedHeight moves the scrollbar range, which re-enters
            self._refit = True
            return
        self._fitting = True
        try:
            for _ in range(2):                         # one settle pass after the scrollbar range changed
                self._refit = False
                n = self.rowCount() if self._fit_limit is None else min(self.rowCount(), self._fit_limit)
                self.setFixedHeight(self._rows_height(n))
                self.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff if n >= self.rowCount() else Qt.ScrollBarAsNeeded)
                if not self._refit:
                    break
        except RuntimeError:                           # a queued fit reached a widget whose C++ object is already gone
            pass
        finally:
            self._fitting = False

    def changeEvent(self, e):
        super().changeEvent(e)
        if getattr(self, "_fit_limit", False) is not False and e.type() in (QEvent.StyleChange, QEvent.FontChange):
            QTimer.singleShot(0, self, self._apply_fit)    # row heights follow the style, which settles after this event

    def _on_current_cell(self, row, _col, _prow, _pcol):
        if self._quiet:
            return
        orig = self._orig(row) if row >= 0 else -1
        if orig != self._last_orig:
            self._last_orig = orig
            self.current_row_changed.emit(orig)

    def current_orig(self):
        """Original index of the current row, or -1."""
        r = self.currentRow()
        return self._orig(r) if r >= 0 and self.item(r, 0) is not None else -1

    def _visual_of(self, orig):
        for r in range(self.rowCount()):
            if self._orig(r) == orig:
                return r
        return -1

    def select_orig(self, orig):
        """Select the row whose original index is `orig` (scrolling it into view) WITHOUT emitting
        current_row_changed; -1, or an index that is not shown, clears the selection. Returns whether a row was selected."""
        self._quiet += 1
        try:
            r = self._visual_of(orig) if orig is not None and orig >= 0 else -1
            if r >= 0:
                self.setCurrentCell(r, 0)
                self.scrollToItem(self.item(r, 0))
            else:
                self.clearSelection()
                self.selectionModel().clearCurrentIndex()
            self._last_orig = orig if r >= 0 else -1
            return r >= 0
        finally:
            self._quiet -= 1

    def keyPressEvent(self, e):
        """Ctrl+C copies the selected rows (or the current row) as tab-separated MACHINE values, the same text as
        'Copy row' — Qt's own copy would put one cell's display text, UI marks (! / ‡) included, on the clipboard."""
        if e.matches(QKeySequence.Copy):
            rows = sorted({i.row() for i in self.selectedIndexes()}) or ([self.currentRow()] if self.currentRow() >= 0 else [])
            if rows:
                QApplication.clipboard().setText(
                    "\n".join("\t".join(self._value(r, j) for j in range(self.columnCount())) for r in rows))
            e.accept()
            return
        if (e.key() == Qt.Key_Menu or (e.key() == Qt.Key_F10 and e.modifiers() & Qt.ShiftModifier)) and self.currentRow() >= 0:
            self.show_row_menu()                     # the keyboard route to the same menu as the right-click
            e.accept()
            return
        super().keyPressEvent(e)

    def _orig(self, visual_row):
        """Map a visual row (post-sort) back to the index it had in set_rows, so row menus /
        activation address the right data record no matter how the user sorted."""
        if visual_row < 0:
            return -1
        it = self.item(visual_row, 0)
        d = it.data(Qt.UserRole) if it else None
        return d if d is not None else visual_row

    def set_rows(self, rows, styles=None, tips=None, exports=None, _from_complete=False,
                 _orig_indices=None):
        """rows: list of row value-lists. Optional styles[i][j] = foreground colour (hex/QColor or None)
        and tips[i][j] = tooltip override — both additive and sort-safe (the item carries its own colour
        and tooltip when a sort moves it). exports[i][j] overrides what export/copy writes for a cell whose
        DISPLAYED text carries a UI-only mark (the ΔG severity ! / ‡), so the exported column stays numeric.
        The previously selected row (by original index) is re-selected without emitting current_row_changed."""
        if not _from_complete:
            self._complete_rows = None
            self._complete_limit = None
        prev = self.current_orig()
        if prev < 0:
            prev = self._last_orig
        self._quiet += 1                             # the rebuild is not the user moving the selection
        try:
            self._fill_rows(rows, styles, tips, exports, _orig_indices)
            r = self._visual_of(prev) if prev >= 0 else -1   # sticky selection: the same original row, if still shown
            if r >= 0:
                self.setCurrentCell(r, 0)
                self._last_orig = prev
        finally:
            self._quiet -= 1
        self._apply_fit()
        if prev >= 0 and r < 0:
            self._last_orig = -1
            self.current_row_changed.emit(-1)

    def _fill_rows(self, rows, styles, tips, exports, _orig_indices):
        self.setSortingEnabled(False)                # never sort mid-insert (it scrambles rows)
        self.setRowCount(0)
        # Numeric columns right-align so their digits form a column; everything else stays centred. The
        # data font is Cascadia Mono precisely so figures are tabular (theme.py: "numeric data stay in
        # Cascadia Mono for column alignment"), and centring threw that away — 1409/110/98/75 sat on four
        # different left edges, which is the one thing a monospaced column is meant to prevent.
        numcols = _numeric_columns(rows, self.columnCount())
        # Stretch the last column only when it holds free text (a Method or Label column genuinely wants the
        # slack). A numeric last column stretched to the viewport strands its digits against the far edge,
        # a column-width away from the value they belong beside — visible on the ORFs table, whose last
        # column is "aa". Text columns keep the old behaviour.
        self.horizontalHeader().setStretchLastSection(not (numcols and numcols[-1]))
        for j, isnum in enumerate(numcols):
            h = self.horizontalHeaderItem(j)
            if h is not None:                        # header follows its column, or the label floats off its digits
                h.setTextAlignment((Qt.AlignRight if isnum else Qt.AlignHCenter) | Qt.AlignVCenter)
        for i, r in enumerate(rows):
            self.insertRow(i)
            for j, c in enumerate(r):
                text = "" if c is None else str(c)
                item = _Cell(text)
                align = Qt.AlignRight if (j < len(numcols) and numcols[j]) else Qt.AlignHCenter
                item.setTextAlignment(align | Qt.AlignVCenter)
                tip = tips[i][j] if (tips and tips[i] and j < len(tips[i]) and tips[i][j]) else text
                item.setToolTip(tip)                     # full value / richer detail on hover
                col = styles[i][j] if (styles and styles[i] and j < len(styles[i])) else None
                if col:
                    item.setForeground(QColor(col) if isinstance(col, str) else col)
                xv = exports[i][j] if (exports and exports[i] and j < len(exports[i]) and exports[i][j] is not None) else None
                if xv is not None:
                    item.setData(EXPORT_ROLE, str(xv))
                if j == 0:
                    item.setData(Qt.UserRole, (_orig_indices[i] if _orig_indices is not None else i))
                self.setItem(i, j, item)
        self.horizontalHeader().setSortIndicator(-1, Qt.AscendingOrder)  # keep engine order until a header is clicked
        self.setSortingEnabled(True)
        self.resizeColumnsToContents()
        if numcols and numcols[-1]:
            # A numeric last column does not stretch (see above), which on a wide card left the table's
            # frame spanning the full width with most of it empty and the row-scroll thumb stranded far
            # from the data. Cap the WIDGET to what its columns actually need, so the table reads as
            # deliberately sized rather than as a broken full-width one. No cap when the last column is
            # text — that column is meant to absorb the slack.
            w = sum(self.columnWidth(j) for j in range(self.columnCount())) + 2 * self.frameWidth()
            if self.verticalScrollBar().isVisible() or self.rowCount() > 8:
                w += self.verticalScrollBar().sizeHint().width()
            self.setMaximumWidth(max(w, 1))
        else:
            self.setMaximumWidth(16777215)                # QWIDGETSIZE_MAX: no cap for a text last column

    def set_complete_rows(self, rows, visible_limit):
        """Bind a capped view to its complete row set; sorting reselects the visible prefix globally."""
        self._complete_rows = [list(row) for row in rows]
        self._complete_limit = max(0, int(visible_limit))
        self.set_rows(self._complete_rows[:self._complete_limit], _from_complete=True,
                      _orig_indices=list(range(min(self._complete_limit, len(self._complete_rows)))))

    @staticmethod
    def _compare_rows(left, right, column):
        a = "" if column >= len(left) or left[column] is None else str(left[column])
        b = "" if column >= len(right) or right[column] is None else str(right[column])
        an, bn = _sortkey(a), _sortkey(b)
        av, bv = ((an, bn) if an is not None and bn is not None else (a.casefold(), b.casefold()))
        return (av > bv) - (av < bv)

    def _sorted_complete_pairs(self):
        pairs = [(i, list(row)) for i, row in enumerate(self._complete_rows or [])]
        section = self.horizontalHeader().sortIndicatorSection()
        if section < 0:
            return pairs
        order = self.horizontalHeader().sortIndicatorOrder()
        return sorted(pairs, key=functools.cmp_to_key(
            lambda a, b: self._compare_rows(a[1], b[1], section)), reverse=order == Qt.DescendingOrder)

    def _sorted_complete_rows(self):
        return [row for _, row in self._sorted_complete_pairs()]

    def _complete_sort_changed(self, section, order):
        if self._complete_rows is None or self._complete_refresh or section < 0:
            return
        pairs = self._sorted_complete_pairs()
        shown = pairs[:self._complete_limit]
        self._complete_refresh = True
        try:
            self.set_rows([row for _, row in shown], _from_complete=True,
                          _orig_indices=[i for i, _ in shown])
            self.sortItems(section, order)
            self.horizontalHeader().setSortIndicator(section, order)
        finally:
            self._complete_refresh = False

    def set_row_menu(self, builder):
        self._row_menu = builder

    def _build_menu(self, row, col, gpos):
        """The row's menu: its own actions, then Copy cell / Copy row / Export table…. One builder serves the
        right-click, the Menu / Shift+F10 key and the visible 'Row actions' button, so they cannot drift apart."""
        m = QMenu(self)
        if row >= 0 and self._row_menu:
            for label, fn in self._row_menu(self._orig(row)):        # pass the original data index, not the sorted row
                m.addAction(label, fn)
            m.addSeparator()
        m.addAction("Copy cell", lambda: self._copy_cell(row, col))
        m.addAction("Copy row", lambda: self._copy_row(row))
        m.addAction("Export table…", lambda: self._do_export(gpos))   # flat action (no submenu -> no native ▶ arrow to collide with)
        return m

    def _menu(self, pos):
        gpos = self.viewport().mapToGlobal(pos)
        self._build_menu(self.rowAt(pos.y()), self.columnAt(pos.x()), gpos).exec(gpos)

    def row_items(self, orig):
        """[(label, callable)] the row menu offers for the row with original index `orig` (empty when none)."""
        return list(self._row_menu(orig)) if self._row_menu and orig is not None and orig >= 0 else []

    def show_row_menu(self, gpos=None):
        """Pop the menu for the CURRENT row (anchored under it unless `gpos` is given). False when no row is current."""
        r = self.currentRow()
        if r < 0:
            return False
        if gpos is None:
            gpos = self.viewport().mapToGlobal(self.visualRect(self.currentIndex()).bottomLeft())
        self._build_menu(r, max(self.currentColumn(), 0), gpos).exec(gpos)
        return True

    def _do_export(self, gpos):
        """Flat table export: pop the arrow-free format picker at gpos, then route the chosen format to export_table."""
        fmt = pick_table_format(self, gpos)
        if fmt:
            # notes too, or the two entry points produce DIFFERENT files: the visible button carried the
            # caveats and this one silently dropped them, so which path a user happened to take decided
            # whether their exported table kept its scientific qualifications.
            export_table(self.export_headers(), self.rows_data(), getattr(self, "export_base", "TEagle_table"),
                         self, fmt=fmt,   # same base as the visible Export button -> identical proposed filename
                         notes=getattr(self, "export_notes", None))

    def _value(self, i, j):
        """The cell's machine value: the export override when one was set, else the displayed text."""
        it = self.item(i, j)
        if it is None:
            return ""
        v = it.data(EXPORT_ROLE)
        return it.text() if v is None else str(v)

    def _copy_row(self, row):
        if row < 0:
            return
        QApplication.clipboard().setText("\t".join(self._value(row, j) for j in range(self.columnCount())))

    def _copy_cell(self, row, col):
        if row < 0 or col < 0:
            return
        QApplication.clipboard().setText(self._value(row, col))

    def export_headers(self):
        """Headers as written to a FILE. `export_header_map` renames a column for export only, so a units
        or coordinate-convention label can reach the file without widening the on-screen column (headers are
        ResizeToContents, so "Start" -> "Start (0-based)" would push the table sideways on every screen).
        On screen the convention is a GLOSS tooltip, which no export can carry."""
        m = getattr(self, "export_header_map", None) or {}
        # a lone '#' header is renamed 'No.': a comment-skipping reader (pandas comment='#', R comment.char) would
        # otherwise read the whole header line as a comment and promote the first data row to column names
        return [m.get(h, "No." if str(h).strip() == "#" else h) for h in self._headers]

    def rows_data(self):
        """Rows for export/copy, in the table's current on-screen (sorted) order.

        `full_rows` overrides this when the table on screen is a capped VIEW of a bigger set. A file that
        silently held the first N of M would be an export defect: what reaches the file
        must not be quietly less than what was computed. Both export entry points — the visible button and
        the right-click menu — read through here, so neither can drift from the other."""
        if self._complete_rows is not None:
            return self._sorted_complete_rows()
        full = getattr(self, "full_rows", None)
        if callable(full):
            return full()

        return [[self._value(i, j) for j in range(self.columnCount())]
                for i in range(self.rowCount())]


# ============================ self-similarity panel (dot plot + heat map) ============================
class _DotCanvas(QWidget):
    """Renders a self-similarity SVG and reports, on hover, exactly which comparison a cell represents.

    A dot matrix is unreadable without a coordinate readout: every cell is a pair of positions, and a
    user cannot be expected to measure one off two rulers. Hovering names both coordinates, the layer
    (direct or inverted) and the match count, so the picture can be interrogated rather than admired."""

    def __init__(self, panel):
        super().__init__()
        self.panel = panel
        self._svg, self._renderer = "", None
        self._w = self._h = 1
        self.setMinimumHeight(300)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)
        self.setMouseTracking(True)

    def _sc(self):
        return ((self.width() or self._w) / self._w) if self._w else 1.0

    def set_svg(self, svg):
        self._svg = svg
        self._renderer = QSvgRenderer(QByteArray(svg.encode("utf-8")))
        self._w, self._h = _svg_size(svg)
        self.update()

    def paintEvent(self, e):
        if self._renderer is not None:
            p = QPainter(self)
            self._renderer.render(p)                       # fills the widget rect -> scaled by the fixed size DotPanel sets
            p.end()

    def wheelEvent(self, e):
        # Ctrl + wheel zooms; a plain wheel scrolls the enclosing scroll area (default handling)
        if e.modifiers() & Qt.ControlModifier:
            self.panel._zoom_by(1.15 if e.angleDelta().y() > 0 else 1.0 / 1.15)
            e.accept()
        else:
            e.ignore()

    def mouseMoveEvent(self, e):
        m = self.panel.matrix
        if not m:
            return
        sc = self._sc() or 1.0
        ML, MT, MR = DOT_ML, DOT_MT, DOT_MR                # one definition, in figures._dot_frame
        plot = max(self._w - ML - MR, 80.0)
        sx, sy = e.position().x() / sc, e.position().y() / sc
        if not (ML <= sx <= ML + plot and MT <= sy <= MT + plot):
            QToolTip.hideText()
            return
        b = m["bins"]
        n = max(m["length"], 1)
        j = min(b - 1, int((sx - ML) / plot * b))
        i = min(b - 1, int((sy - MT) / plot * b))
        bp_x = int((j + 0.5) / b * n)
        bp_y = int((i + 0.5) / b * n)
        fwd, rev = m["forward"][i][j], m["reverse"][i][j]
        thr = self.panel.threshold
        parts = [f"{bp_y} bp  vs  {bp_x} bp"]
        if fwd:
            parts.append(f"direct (forward): {fwd} match{'es' if fwd != 1 else ''}"
                         + ("" if fwd >= thr else "  — at or below the chance level"))
        if rev:
            parts.append(f"inverted (reverse complement): {rev} match{'es' if rev != 1 else ''}"
                         + ("" if rev >= thr else "  — at or below the chance level"))
        if not fwd and not rev:
            parts.append(f"no exact {m['k']}-mer match between these two regions")
        if abs(i - j) <= 1 and fwd:
            parts.append("on the identity diagonal — every sequence matches itself here")
        QToolTip.showText(e.globalPosition().toPoint(), "\n".join(parts), self)


class DotPanel(QWidget):
    """Self-similarity of one locus, as a dot plot or a binned heat map.

    Two views of one computation: the dot plot answers "is there a repeat, and where"; the heat map
    answers "how much of the locus takes part", which a binary mark cannot show once bins saturate."""

    HELP = ("Compares the sequence with itself and marks every position pair that shares an exact word.\n\n"
            "• The solid line corner-to-corner is the identity diagonal — every sequence matches itself.\n"
            "• A block OFF that diagonal is a DIRECT repeat: the two LTRs of a retroelement.\n"
            "• A block on the ANTI-diagonal is an INVERTED repeat: the two TIRs of a DNA transposon.\n"
            "• Shaded bands mark what TEagle itself detected, so the picture can confirm or contradict it.\n\n"
            "This finds repeats the targeted detectors were not looking for — a strong block where nothing "
            "was called is worth investigating. It cannot do the reverse: exact word matching is not an "
            "alignment, so a diverged repeat fades out, and a faint or absent block is NOT evidence that "
            "no repeat exists. A repeat shorter than the word size cannot appear at all.")

    _BASE_W = 900                    # the SVG is authored once at this logical width; zoom scales the display

    def __init__(self, base_name="TEagle_selfsim", parent=None):
        super().__init__(parent)
        import figures
        self.base_name = base_name
        self.matrix = None
        self.guides = []
        self.threshold = 1
        self.mode = "dot"            # 'dot' | 'heat'
        self.layer = "forward"       # heat map layer
        self.show_guides = True
        self.theme = "dark"
        self._theme_locked = False
        self._export_as_shown = False   # False: exports use the publication palette; True: the on-screen palette
        self.zoom = None             # None = fit the whole plot to the window; a float = explicit zoom factor
        self._def_fwd, self._def_rev = figures._DOT_FWD, figures._DOT_REV
        self.fwd_color, self.rev_color = self._def_fwd, self._def_rev
        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(6)

        # row 1 — which view
        bar = QHBoxLayout()
        bar.setSpacing(6)
        self._btns = {}
        for key, label, tip in (
                ("dot", "Dot plot", "One mark per matching position pair — shows WHERE repeats are"),
                ("heat", "Heat map", "Match density per bin — shows HOW MUCH of the locus takes part, "
                                     "which a binary mark cannot once bins saturate")):
            b = QPushButton(label)
            b.setProperty("sm", True)
            b.setToolTip(tip)
            b.clicked.connect(lambda _=False, k=key: self._set_mode(k))
            self._btns[key] = b
            bar.addWidget(b)
        self.layerBtn = QPushButton("Inverted layer")
        self.layerBtn.setProperty("sm", True)
        self.layerBtn.setToolTip("Heat map only: switch between direct (LTR-type) and inverted "
                                 "(TIR-type) match density")
        self.layerBtn.clicked.connect(self._toggle_layer)
        bar.addWidget(self.layerBtn)
        self.guideBtn = QPushButton("Hide guides")
        self.guideBtn.setProperty("sm", True)
        self.guideBtn.setToolTip("Shaded bands mark the terminal repeats TEagle detected. Hide them to "
                                 "read the raw self-similarity without that prompt.")
        self.guideBtn.clicked.connect(self._toggle_guides)
        bar.addWidget(self.guideBtn)
        bar.addStretch(1)
        # the help text is a button's tooltip: on the canvas it was overwritten, within ~20 ms of every mouse move, by
        # the per-cell coordinate readout it shares the tooltip slot with
        self.helpBtn = QPushButton("?")
        self.helpBtn.setProperty("sm", True)
        _name_button(self.helpBtn, "About this plot", self.HELP)
        bar.addWidget(self.helpBtn)
        lay.addLayout(bar)

        # row 2 — zoom · export colours · export formats
        tools = QHBoxLayout()
        tools.setSpacing(6)
        for label, name, tip, slot in (("−", "Zoom out", "Zoom out (or Ctrl + mouse wheel)", lambda: self._zoom_by(1 / 1.25)),
                                       ("Fit", "Fit the figure to the window", "Fit the whole plot to the window", self._fit),
                                       ("+", "Zoom in", "Zoom in (or Ctrl + mouse wheel)", lambda: self._zoom_by(1.25))):
            b = QPushButton(label)
            b.setProperty("sm", True)
            _name_button(b, name, tip)
            b.clicked.connect(slot)
            tools.addWidget(b)
        tools.addSpacing(12)
        self.fwdBtn = QPushButton("Forward")
        self.fwdBtn.setProperty("sm", True)
        self.fwdBtn.setToolTip("Colour of the direct-repeat (forward) marks — applies on screen and in every export")
        self.fwdBtn.clicked.connect(lambda: self._pick_color("fwd"))
        self.revBtn = QPushButton("Reverse")
        self.revBtn.setProperty("sm", True)
        self.revBtn.setToolTip("Colour of the inverted-repeat (reverse-complement) marks")
        self.revBtn.clicked.connect(lambda: self._pick_color("rev"))
        rst = QPushButton("Reset colours")
        rst.setProperty("sm", True)
        rst.setToolTip("Restore the colour-vision-safe Okabe–Ito default pair")
        rst.clicked.connect(self._reset_colors)
        tools.addWidget(self.fwdBtn)
        tools.addWidget(self.revBtn)
        tools.addWidget(rst)
        tools.addStretch(1)
        self.exportPalBtn = QPushButton()
        self.exportPalBtn.setProperty("sm", True)
        self.exportPalBtn.clicked.connect(self._toggle_export_palette)
        tools.addWidget(self.exportPalBtn)
        self._sync_export_palette_btn()
        for label, name, tip, slot in (
                ("SVG", "Export the figure as SVG", "Vector export — text stays text, legible on white", self._export_svg),
                ("PNG", "Export the figure as PNG", "Raster export at publication resolution", self._export_png),
                ("PDF", "Export the figure as vector PDF", "Vector PDF — one journal-ready figure per page", self._export_pdf)):
            b = QPushButton(label)
            b.setProperty("sm", True)
            _name_button(b, name, tip)
            b.clicked.connect(slot)
            tools.addWidget(b)
        lay.addLayout(tools)

        # the canvas lives in a scroll area so a zoomed plot can exceed the window and be panned
        self.canvas = _DotCanvas(self)
        self._scroll = QScrollArea()
        self._scroll.setWidgetResizable(False)
        self._scroll.setAlignment(Qt.AlignCenter)
        self._scroll.setWidget(self.canvas)
        lay.addWidget(self._scroll, 1)

        self.note = QLabel("")
        self.note.setObjectName("cardnote")             # running prose: the sans face, 12px floor
        self.note.setWordWrap(True)
        lay.addWidget(self.note)
        self._sync_swatches()

    def _swatch(self, color):
        pm = QPixmap(12, 12)
        pm.fill(QColor(color))
        return QIcon(pm)

    def _sync_swatches(self):
        self.fwdBtn.setIcon(self._swatch(self.fwd_color))
        self.revBtn.setIcon(self._swatch(self.rev_color))

    def _pick_color(self, which):
        cur = self.fwd_color if which == "fwd" else self.rev_color
        c = QColorDialog.getColor(QColor(cur), self, "Choose mark colour")
        if not c.isValid():
            return
        if which == "fwd":
            self.fwd_color = c.name()
        else:
            self.rev_color = c.name()
        self._sync_swatches()
        self._render()

    def _reset_colors(self):
        self.fwd_color, self.rev_color = self._def_fwd, self._def_rev
        self._sync_swatches()
        self._render()

    def _apply_zoom(self):
        if not self.matrix or self.canvas._w <= 1:
            return
        if self.zoom is None:                             # fit the WHOLE plot (both axes) into the viewport
            vp = self._scroll.viewport()
            z = min((vp.width() - 4) / self.canvas._w, (vp.height() - 4) / self.canvas._h)
            z = max(z, 0.15)
        else:
            z = self.zoom
        self.canvas.setFixedSize(max(int(self.canvas._w * z), 60), max(int(self.canvas._h * z), 60))

    def _zoom_by(self, f):
        cur = (self.canvas.width() / self.canvas._w) if self.canvas._w > 1 else 1.0
        self.zoom = max(0.2, min(cur * f, 8.0))
        self._apply_zoom()

    def _fit(self):
        self.zoom = None
        self._apply_zoom()

    def resizeEvent(self, e):
        super().resizeEvent(e)
        if self.zoom is None:                             # fit follows the window until the user zooms explicitly
            self._apply_zoom()

    def set_matrix(self, matrix, guides=None, threshold=1, scope=""):
        self.matrix = matrix
        self.guides = guides or []
        self.threshold = threshold or 1
        self.note.setText(scope)
        self._render()

    def apply_app_theme(self, app_theme):
        if not self._theme_locked:
            self.theme = "white" if app_theme == "light" else "dark"
            self._render()

    def _set_mode(self, mode):
        self.mode = mode
        self._render()

    def _toggle_layer(self):
        self.layer = "reverse" if self.layer == "forward" else "forward"
        self._render()

    def _toggle_guides(self):
        self.show_guides = not self.show_guides
        set_caps(self.guideBtn, "Show guides" if not self.show_guides else "Hide guides")
        self._render()

    def _toggle_export_palette(self):
        self._export_as_shown = not self._export_as_shown
        self._sync_export_palette_btn()

    def _sync_export_palette_btn(self):
        shown = self._export_as_shown
        set_caps(self.exportPalBtn, "Palette: as shown" if shown else "Palette: publication")
        self.exportPalBtn.setAccessibleName("Export palette: " + ("as shown on screen" if shown else "publication"))
        self.exportPalBtn.setToolTip(
            "Palette of the exported SVG, PNG and PDF. 'As shown' writes the plot exactly as it looks now. "
            "'Publication' writes dark ink on white, legible on paper. Click to switch.")

    def _svg(self, w, for_export=False):
        import figures
        g = self.guides if self.show_guides else None
        if self.mode == "heat":
            return figures.svg_dotheat(self.matrix, W=w, theme=self.theme, for_export=for_export,
                                       guides=g, which=self.layer, fwd=self.fwd_color, rev=self.rev_color)
        return figures.svg_dotplot(self.matrix, W=w, theme=self.theme, for_export=for_export,
                                   guides=g, read_threshold=self.threshold, fwd=self.fwd_color, rev=self.rev_color)

    def _render(self):
        if not self.matrix:
            return
        for key, b in self._btns.items():
            b.setProperty("primary", key == self.mode)
            b.style().unpolish(b)
            b.style().polish(b)
        self.layerBtn.setEnabled(self.mode == "heat")
        set_caps(self.layerBtn, "Direct layer" if self.layer == "reverse" else "Inverted layer")
        # author the SVG once at a fixed logical width; the zoom factor scales the display, not the SVG
        self.canvas.set_svg(self._svg(self._BASE_W))
        self._apply_zoom()

    def _export_svg(self):
        path, _ = QFileDialog.getSaveFileName(self, "Export SVG", self.base_name + ".svg", "SVG (*.svg)")
        if path:
            path = _with_ext(path, ".svg")
            safe_export(self, path, lambda p: save_svg(self._svg(920, for_export=not self._export_as_shown), p))

    def _export_png(self):
        path, _ = QFileDialog.getSaveFileName(self, "Export PNG", self.base_name + ".png", "PNG (*.png)")
        if path:
            path = _with_ext(path, ".png")
            safe_export(self, path, lambda p: render_png(self._svg(920, for_export=not self._export_as_shown), p))

    def _export_pdf(self):
        path, _ = QFileDialog.getSaveFileName(self, "Export PDF", self.base_name + ".pdf", "PDF (*.pdf)")
        if path:
            path = _with_ext(path, ".pdf")
            safe_export(self, path, lambda p: render_pdf(self._svg(920, for_export=not self._export_as_shown), p))
