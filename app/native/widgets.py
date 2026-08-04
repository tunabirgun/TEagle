"""Reusable Qt widgets for the native app: an SVG figure panel (zoom/pan/WYSIWYG export),
an interactive genome viewer (windowed semantic zoom), and a data table with CSV/TSV export
and a copy context menu. All figure rendering goes through QSvgRenderer so on-screen output
matches the exported SVG/PNG (the gel's Gaussian-blur glow is the one SVG-Tiny casualty on screen)."""
from __future__ import annotations
import contextlib, math, os, re

from PySide6.QtCore import Qt, QByteArray, QPointF, QRectF, QSize, Signal
from PySide6.QtGui import QImage, QPainter, QColor, QPixmap, QIcon
from PySide6.QtSvg import QSvgRenderer
from PySide6.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QLabel, QSizePolicy,
                               QTableWidget, QTableWidgetItem, QMenu, QFileDialog, QApplication,
                               QAbstractItemView, QHeaderView, QToolTip, QProgressBar, QScrollArea,
                               QColorDialog)

import theme                       # UI_SCALE: figures are authored in logical units and stretched by it
from figures import gv_gutter, GV_MR   # label gutter is model-dependent; the hit-test must use the same value

MODE_LABEL = {"dark": "dark", "white": "light", "uv": "UV", "mono": "mono", "transparent": "transparent"}

# per-cell export/copy value, when the displayed text carries a UI-only mark (see DataTable.set_rows)
EXPORT_ROLE = Qt.UserRole + 1


def uppercase_buttons(root):
    """Uppercase every action button under `root` (the web '.btn' look). Skips link buttons and card
    headers — their titles stay sentence-case — and leaves glyphs untouched. Idempotent. Lives here so a
    dialog built outside MainWindow (installer, notification, genome manager) can call it as it is built,
    instead of waiting for the next theme/scale/result refresh to reach it."""
    for b in root.findChildren(QPushButton):
        if b.property("link") or b.objectName() == "cardhdr":
            continue
        t = b.text()
        if t and t != t.upper():
            b.setText(t.upper())


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
    with atomic_write(path, "w", encoding="utf-8") as f:
        f.write(svg)


def render_png(svg: str, path: str, scale: int = 3):
    """Rasterise an SVG string to a transparent-background PNG at `scale`x for publication."""
    w, h = _svg_size(svg)
    r = QSvgRenderer(QByteArray(svg.encode("utf-8")))
    img = QImage(int(w * scale), int(h * scale), QImage.Format_ARGB32)
    img.fill(Qt.transparent)
    p = QPainter(img)
    r.render(p)
    p.end()
    atomic_save(path, lambda tmp: img.save(tmp, "PNG"))


def render_pdf(svg: str, path: str):
    """Write an SVG string to a VECTOR PDF (a journal-ready figure): SVG -> QPainter -> QPdfWriter, never
    via a raster QImage, so text stays text and the marks stay resolution-independent. The page is sized to
    the figure's own aspect ratio (px -> mm at 96 dpi) with no margin, so the whole plot fills one page."""
    from PySide6.QtGui import QPdfWriter, QPageSize, QPageLayout
    from PySide6.QtCore import QSizeF, QMarginsF, QRectF
    w, h = _svg_size(svg)
    writer = QPdfWriter(path)
    writer.setResolution(96)                                          # 1 SVG px == 1 device px at 96 dpi
    writer.setPageSize(QPageSize(QSizeF(w * 25.4 / 96.0, h * 25.4 / 96.0), QPageSize.Unit.Millimeter))
    writer.setPageMargins(QMarginsF(0, 0, 0, 0), QPageLayout.Unit.Millimeter)
    r = QSvgRenderer(QByteArray(svg.encode("utf-8")))
    p = QPainter(writer)
    r.render(p, QRectF(0, 0, writer.width(), writer.height()))       # fill the page in device pixels
    p.end()


class SvgCanvas(QWidget):
    """Paints an SVG string with a user scale + pan offset. Wheel = zoom at cursor, drag = pan."""
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
        self.setMinimumHeight(300)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.setMouseTracking(True)

    def set_hit_regions(self, regions, on_menu=None):
        self.regions = regions or []
        self._on_menu = on_menu

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
        items = self._on_menu(r)
        if not items:
            return
        m = QMenu(self)
        for label, fn in items:
            m.addAction(label, fn)
        m.exec(e.globalPos())

    def paintEvent(self, _):
        p = QPainter(self)
        p.setRenderHint(QPainter.Antialiasing)
        p.setRenderHint(QPainter.SmoothPixmapTransform)
        if self._renderer:
            p.translate(self.tx, self.ty)
            p.scale(self.scale, self.scale)
            self._renderer.render(p, QRectF(0, 0, self._sw, self._sh))
        p.end()


class FigurePanel(QWidget):
    """Toolbar (bg modes + zoom + fit + export) over an SvgCanvas. `build_fn(bg)->svg` supplies
    the figure; export is WYSIWYG — it writes the currently selected background mode."""
    def __init__(self, build_fn, base_name: str, modes=("dark", "white"), parent=None,
                 hit_regions=None, on_menu=None):
        super().__init__(parent)
        self.build_fn = build_fn
        self.base_name = base_name
        self.modes = list(modes)
        self.bg = self.modes[0]
        self._theme_locked = False        # a manual bg pick pins this panel; app-theme toggles stop following it
        self._hit_regions = hit_regions
        self._on_menu = on_menu
        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(4)
        bar = QHBoxLayout()
        bar.addWidget(QLabel("bg"))
        self._mode_btns = {}
        for m in self.modes:
            b = QPushButton(MODE_LABEL.get(m, m).upper())
            b.setProperty("sm", True)
            b.clicked.connect(lambda _=False, mm=m: self._set_bg(mm))
            bar.addWidget(b)
            self._mode_btns[m] = b
        bar.addStretch(1)
        for txt, fn in (("−", self._zoom_out), ("FIT", self._fit), ("+", self._zoom_in),
                        ("SVG", self._export_svg), ("PNG", self._export_png)):
            b = QPushButton(txt)
            b.setProperty("sm", True)
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
            save_svg(self.build_fn(self.bg), path)                 # export what you see (selected bg mode)

    def _export_png(self):
        path, _ = QFileDialog.getSaveFileName(self, "Export PNG", self.base_name + ".png", "PNG (*.png)")
        if path:
            render_png(self.build_fn(self.bg), path)               # export what you see (selected bg mode)


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
        self.view = {"start": 0.0, "end": 1.0}
        self.on_feature_menu = None       # callable(region) -> list[(label, fn)] for right-click
        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(4)
        bar = QHBoxLayout()
        bar.addWidget(QLabel("bg"))
        self._th_btns = {}
        for th, lab in (("dark", "DARK"), ("white", "LIGHT")):
            b = QPushButton(lab); b.setProperty("sm", True)
            b.clicked.connect(lambda _=False, t=th: self._set_theme(t, user=True))
            bar.addWidget(b); self._th_btns[th] = b
        self.pos = QLabel(""); self.pos.setObjectName("gvpos")
        bar.addWidget(self.pos)
        bar.addStretch(1)
        for txt, fn in (("−", lambda: self._zoom(1.6)), ("FIT", self._fit), ("+", lambda: self._zoom(0.625)),
                        ("SVG", self._export_svg), ("PNG", self._export_png)):
            b = QPushButton(txt); b.setProperty("sm", True); b.clicked.connect(fn); bar.addWidget(b)
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
        self.theme = t
        self._render()

    def apply_app_theme(self, app_theme):
        """Follow the app's dark/light theme (app 'light' -> viewer 'white') UNTIL the user picks a bg
        manually, after which this viewer keeps that choice and ignores app-theme toggles. Re-renders in
        place so the current pan/zoom window (self.view) is preserved."""
        if self._theme_locked:            # user's manual bg pick wins over app-theme propagation
            return
        self._set_theme("white" if app_theme == "light" else "dark")

    def _cur_svg(self, w, for_export=False, theme=None):
        return self._svg_genome(self.model, {"start": self.view["start"], "end": self.view["end"]},
                                w, theme or self.theme, for_export)

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
        L = self.model.get("length", 1) or 1
        self.pos.setText(f"{int(self.view['start']):,}–{int(self.view['end']):,} bp · "
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
        """Which palette a file export should use.

        The figure leaves the app to become a journal figure, so the DEFAULT export is the publication
        palette — dark ink, legible on white paper — regardless of the on-screen background, which exists
        for reading the viewer rather than for print. An explicitly chosen background still wins: if the
        user picked one from the bg buttons, that is a deliberate instruction and is honoured."""
        return {"for_export": not self._theme_locked, "theme": self.theme}

    def _export_svg(self):
        path, _ = QFileDialog.getSaveFileName(self, "Export SVG", self.base_name + ".svg", "SVG (*.svg)")
        if path:
            save_svg(self._cur_svg(920, **self._export_target()), path)

    def _export_png(self):
        path, _ = QFileDialog.getSaveFileName(self, "Export PNG", self.base_name + ".png", "PNG (*.png)")
        if path:
            render_png(self._cur_svg(920, **self._export_target()), path)


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
            dbp = (e.position().x() - x0) / self._plot_w() * (e0 - s0)
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
        cb = getattr(self.panel, "on_feature_menu", None)
        if not cb:
            return
        r = self._region_at(e.pos().x(), e.pos().y())
        if not r:
            return
        items = cb(r)
        if not items:
            return
        m = QMenu(self)
        for label, fn in items:
            m.addAction(label, fn)
        m.exec(e.globalPos())

    def keyPressEvent(self, e):
        v = self.panel.view
        sp = v["end"] - v["start"]
        k = e.key()
        if k == Qt.Key_Left:
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
    text. A leading formula char is neutralised the same way as the CSV path (CWE-1236)."""
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
    if s[:1] in ("=", "@") or (s[:1] in ("+", "-") and len(s) > 1):    # bare +/- (strand) stays literal
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
    atomic_save(path, lambda tmp: wb.save(tmp))


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
    # thing AGENTS.md rule 2 forbids, on the most ordinary failure there is (a read-only folder, a
    # disconnected drive, the file open in Excel).
    name = os.path.basename(path)
    try:
        write_table(headers, rows, path, notes)
    except PermissionError:
        _export_failed(parent, f"“{name}” could not be written — it may be open in another program "
                               "(Excel locks a file while it is open), or the folder may be read-only. "
                               "Close it or choose another folder, then export again.")
    except OSError as e:
        _export_failed(parent, f"“{name}” could not be written. If the folder is on a network or removable "
                               f"drive, check it is still connected, or choose a folder on this computer. "
                               f"({type(e).__name__})")


def _export_failed(parent, msg):
    """Report an export failure through the parent's banner when it has one, else a message box. The table
    widgets are reused by dialogs that are not MainWindow, so this cannot assume _banner exists."""
    banner = getattr(parent, "_banner", None)
    if callable(banner):
        banner(msg, "warn")
        return
    from PySide6.QtWidgets import QMessageBox
    box = QMessageBox(parent)
    box.setWindowTitle("Export failed")
    box.setIcon(QMessageBox.Warning)
    box.setText(msg)
    box.exec()


def serialize_table(headers, rows, sep=",", notes=None) -> str:
    """The exact text written for a CSV/TSV export. Separated from the file dialog so the contract —
    every displayed row present, in order, escaped so the delimiter survives a round trip — is testable
    without a GUI. An export that silently drops a column or mangles a value is a defect, and a defect
    that only a human reading a spreadsheet can catch is one that ships."""
    lines = [sep.join(_csv_escape(h, sep) for h in headers)]
    lines += [sep.join(_csv_escape(c, sep) for c in r) for r in rows]
    # The scientific qualifications shown beside the table travel WITH it. AGENTS.md: an export that
    # silently drops a hedge or a units label is a defect — and these tables become figures and supplementary
    # files in papers, where the screen the numbers came from is long gone. Comment-prefixed so the data rows
    # stay machine-readable for anything that skips '#'.
    for n in (notes or []):
        lines += ["# " + ln for ln in str(n).splitlines() if ln.strip()]
    return "\r\n".join(lines)


def write_table(headers, rows, path, notes=None):
    """Write a table to `path`, choosing the format from its extension. `notes` are the on-screen caveats."""
    low = str(path).lower()
    if low.endswith(".xlsx") and _HAS_XLSX:
        if _export_xlsx(headers, rows, path, notes) is not False:
            return
        # openpyxl was locatable but would not import: write the same data as CSV under the chosen name's
        # stem rather than failing. _xl() has already cleared _HAS_XLSX, so the option disappears next time.
        path = path[:-5] + ".csv"
        low = path.lower()
    sep = "\t" if low.endswith(".tsv") else ","
    with atomic_write(path, "w", encoding="utf-8-sig", newline="") as f:   # BOM so Excel reads UTF-8
        f.write(serialize_table(headers, rows, sep, notes))


def save_fasta(fasta: str, base: str, parent=None):
    """Write a ready-made FASTA string to a user-chosen .fasta file."""
    path, _ = QFileDialog.getSaveFileName(parent, "Export FASTA", base + ".fasta",
                                          "FASTA (*.fasta *.fa *.txt)")
    if not path:
        return
    if not fasta.endswith("\n"):
        fasta += "\n"
    with atomic_write(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(fasta)


_NUM_RE = re.compile(r'-?\d+\.?\d*(?:[eE][-+]?\d+)?')

def _sortkey(text: str):
    """Numeric key for a cell: a slash-composite like '2/1' or '60.1/58.3' (Mism/Tm/GC F/R) sorts by
    the sum of both values; otherwise the leading number (handles 1.2e-30, 45%, '0–5146' → 0). None
    when the cell has no number, so those columns fall back to case-insensitive string order."""
    if not text:
        return None
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
    activation callback (double-click / Enter) and a right-click menu builder for FASTA-style actions."""
    row_activated = Signal(int)

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

    def sizeHint(self):
        """Height the CURRENT rows need — QAbstractScrollArea's stock 192px hint reserved 80–135px of empty
        grid under a short result. Callers still cap tall tables with setMaximumHeight()."""
        s = super().sizeHint()
        hh = max(self.horizontalHeader().height(), self.horizontalHeader().sizeHint().height())
        h = 2 * self.frameWidth() + hh + sum(self.rowHeight(i) for i in range(self.rowCount()))
        hb = self.horizontalScrollBar()
        if hb.isVisible() or hb.maximum() > 0:            # a sideways-scrolling table must keep room for its bar
            h += hb.sizeHint().height()
        return QSize(s.width(), max(h, self.minimumSizeHint().height()))

    def _orig(self, visual_row):
        """Map a visual row (post-sort) back to the index it had in set_rows, so row menus /
        activation address the right data record no matter how the user sorted."""
        if visual_row < 0:
            return -1
        it = self.item(visual_row, 0)
        d = it.data(Qt.UserRole) if it else None
        return d if d is not None else visual_row

    def set_rows(self, rows, styles=None, tips=None, exports=None):
        """rows: list of row value-lists. Optional styles[i][j] = foreground colour (hex/QColor or None)
        and tips[i][j] = tooltip override — both additive and sort-safe (the item carries its own colour
        and tooltip when a sort moves it). exports[i][j] overrides what export/copy writes for a cell whose
        DISPLAYED text carries a UI-only mark (the ΔG severity ! / ‡), so the exported column stays numeric."""
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
                    item.setData(Qt.UserRole, i)         # remember the original row index for menus
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

    def set_row_menu(self, builder):
        self._row_menu = builder

    def _menu(self, pos):
        row = self.rowAt(pos.y())
        m = QMenu(self)
        if row >= 0 and self._row_menu:
            for label, fn in self._row_menu(self._orig(row)):        # pass the original data index, not the sorted row
                m.addAction(label, fn)
            m.addSeparator()
        m.addAction("Copy row", lambda: self._copy_row(row))
        gpos = self.viewport().mapToGlobal(pos)                   # flat action (no submenu -> no native ▶ arrow to collide with)
        m.addAction("Export table…", lambda: self._do_export(gpos))
        m.exec(gpos)

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

    def export_headers(self):
        """Headers as written to a FILE. `export_header_map` renames a column for export only, so a units
        or coordinate-convention label can reach the file without widening the on-screen column (headers are
        ResizeToContents, so "Start" -> "Start (0-based)" would push the table sideways on every screen).
        On screen the convention is a GLOSS tooltip, which no export can carry."""
        m = getattr(self, "export_header_map", None) or {}
        return [m.get(h, h) for h in self._headers]

    def rows_data(self):
        """Rows for export/copy, in the table's current on-screen (sorted) order.

        `full_rows` overrides this when the table on screen is a capped VIEW of a bigger set. A file that
        silently held the first N of M would be the export defect AGENTS.md names: what reaches the file
        must not be quietly less than what was computed. Both export entry points — the visible button and
        the right-click menu — read through here, so neither can drift from the other."""
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
        # geometry mirrors figures._dot_frame; kept in one place there and read back here
        ML, MT, MR, MB = 58.0, 46.0, 16.0, 66.0
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
        parts = [f"{bp_y:,} bp  vs  {bp_x:,} bp"]
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
        lay.addLayout(bar)

        # row 2 — zoom · export colours · export formats
        tools = QHBoxLayout()
        tools.setSpacing(6)
        for label, tip, slot in (("−", "Zoom out (or Ctrl + mouse wheel)", lambda: self._zoom_by(1 / 1.25)),
                                 ("Fit", "Fit the whole plot to the window", self._fit),
                                 ("+", "Zoom in (or Ctrl + mouse wheel)", lambda: self._zoom_by(1.25))):
            b = QPushButton(label)
            b.setProperty("sm", True)
            b.setToolTip(tip)
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
        for label, tip, slot in (
                ("SVG", "Vector export — text stays text, legible on white", self._export_svg),
                ("PNG", "Raster export at publication resolution", self._export_png),
                ("PDF", "Vector PDF — one journal-ready figure per page", self._export_pdf)):
            b = QPushButton(label)
            b.setProperty("sm", True)
            b.setToolTip(tip)
            b.clicked.connect(slot)
            tools.addWidget(b)
        lay.addLayout(tools)

        # the canvas lives in a scroll area so a zoomed plot can exceed the window and be panned
        self.canvas = _DotCanvas(self)
        self.canvas.setToolTip(self.HELP)
        self._scroll = QScrollArea()
        self._scroll.setWidgetResizable(False)
        self._scroll.setAlignment(Qt.AlignCenter)
        self._scroll.setWidget(self.canvas)
        lay.addWidget(self._scroll, 1)

        self.note = QLabel("")
        self.note.setObjectName("cardmeta")
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
        self.guideBtn.setText("Show guides" if not self.show_guides else "Hide guides")
        self._render()

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
        self.layerBtn.setText("Direct layer" if self.layer == "reverse" else "Inverted layer")
        # author the SVG once at a fixed logical width; the zoom factor scales the display, not the SVG
        self.canvas.set_svg(self._svg(self._BASE_W))
        self._apply_zoom()

    def _export_svg(self):
        path, _ = QFileDialog.getSaveFileName(self, "Export SVG", self.base_name + ".svg", "SVG (*.svg)")
        if path:
            save_svg(self._svg(920, for_export=not self._theme_locked), path)

    def _export_png(self):
        path, _ = QFileDialog.getSaveFileName(self, "Export PNG", self.base_name + ".png", "PNG (*.png)")
        if path:
            render_png(self._svg(920, for_export=not self._theme_locked), path)

    def _export_pdf(self):
        path, _ = QFileDialog.getSaveFileName(self, "Export PDF", self.base_name + ".pdf", "PDF (*.pdf)")
        if path:
            render_pdf(self._svg(920, for_export=not self._theme_locked), path)
