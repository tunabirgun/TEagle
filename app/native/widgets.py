"""Reusable Qt widgets for the native app: an SVG figure panel (zoom/pan/WYSIWYG export),
an interactive genome viewer (windowed semantic zoom), and a data table with CSV/TSV export
and a copy context menu. All figure rendering goes through QSvgRenderer so on-screen output
matches the exported SVG/PNG (the gel's Gaussian-blur glow is the one SVG-Tiny casualty on screen)."""
from __future__ import annotations
import math, re

from PySide6.QtCore import Qt, QByteArray, QPointF, QRectF, Signal
from PySide6.QtGui import QImage, QPainter, QColor
from PySide6.QtSvg import QSvgRenderer
from PySide6.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QLabel, QSizePolicy,
                               QTableWidget, QTableWidgetItem, QMenu, QFileDialog, QApplication,
                               QAbstractItemView, QHeaderView, QToolTip)

MODE_LABEL = {"dark": "dark", "white": "light", "uv": "UV", "mono": "mono", "transparent": "transparent"}


def _svg_size(svg: str):
    m = None
    import re
    mm = re.search(r'viewBox="0 0 ([\d.]+) ([\d.]+)"', svg)
    return (float(mm.group(1)), float(mm.group(2))) if mm else (800.0, 400.0)


def save_svg(svg: str, path: str):
    with open(path, "w", encoding="utf-8") as f:
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
    img.save(path, "PNG")


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
        self.update()

    def _zoom_at(self, cx, cy, factor):
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
                        ("⭳ SVG", self._export_svg), ("⭳ PNG", self._export_png)):
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
        self.bg = m
        self.render()

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
            b.clicked.connect(lambda _=False, t=th: self._set_theme(t))
            bar.addWidget(b); self._th_btns[th] = b
        self.pos = QLabel(""); self.pos.setObjectName("gvpos")
        bar.addWidget(self.pos)
        bar.addStretch(1)
        for txt, fn in (("−", lambda: self._zoom(1.6)), ("FIT", self._fit), ("+", lambda: self._zoom(0.625)),
                        ("⭳ SVG", self._export_svg), ("⭳ PNG", self._export_png)):
            b = QPushButton(txt); b.setProperty("sm", True); b.clicked.connect(fn); bar.addWidget(b)
        lay.addLayout(bar)
        self.canvas = _GenomeCanvas(self)
        self.canvas.setMinimumWidth(320)                  # match the SVG authoring floor so it never CSS-stretches
        lay.addWidget(self.canvas)                         # (below 320 the bp<->pixel map would drift; obs: genome hit-test)

    def set_model(self, model: dict):
        self.model = model
        L = model.get("length", 1) or 1
        self.view = {"start": 0.0, "end": float(L)}
        self._render()

    def set_feature_menu(self, cb):
        """cb(region) -> list[(label, fn)] built on right-click over a feature glyph."""
        self.on_feature_menu = cb

    def _set_theme(self, t):
        self.theme = t
        self._render()

    def _cur_svg(self, w, for_export=False, theme=None):
        return self._svg_genome(self.model, {"start": self.view["start"], "end": self.view["end"]},
                                w, theme or self.theme, for_export)

    def _render(self):
        for th, b in self._th_btns.items():
            b.setProperty("primary", th == self.theme)
            b.style().unpolish(b); b.style().polish(b)
        w = max(self.canvas.width() or 620, 320)
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

    def _export_svg(self):
        path, _ = QFileDialog.getSaveFileName(self, "Export SVG", self.base_name + ".svg", "SVG (*.svg)")
        if path:
            save_svg(self._cur_svg(920, for_export=False, theme=self.theme), path)   # honor the selected bg

    def _export_png(self):
        path, _ = QFileDialog.getSaveFileName(self, "Export PNG", self.base_name + ".png", "PNG (*.png)")
        if path:
            render_png(self._cur_svg(920, for_export=False, theme=self.theme), path)  # honor the selected bg


class _GenomeCanvas(QWidget):
    ML, MR = 96, 16

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

    def _region_at(self, wx, wy):
        sc = ((self.width() or self._w) / self._w) if self._w else 1.0
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
        self.setMinimumHeight(int(self._h))
        self.update()

    def _plot_w(self):
        return max(120.0, (self.width() or 620) - self.ML - self.MR)

    def _x_to_bp(self, px):
        frac = max(0.0, min(1.0, (px - self.ML) / self._plot_w()))
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
            sc = (self.width() or self._w) / self._w
            self._renderer.render(p, QRectF(0, 0, self._w * sc, self._h * sc))
        p.end()


# ---------- data tables ----------
def _csv_escape(v, sep):
    v = "" if v is None else str(v)
    # neutralise spreadsheet formula injection (CWE-1236) but keep a bare +/- (e.g. a strand cell) intact
    if v[:1] in ("=", "@", "\t", "\r") or (v[:1] in ("+", "-") and len(v) > 1):
        v = "'" + v
    if sep in v or '"' in v or "\n" in v:
        v = '"' + v.replace('"', '""') + '"'
    return v


try:                                                     # Excel export is optional: degrade to CSV/TSV if openpyxl is absent
    from openpyxl import Workbook as _XlWorkbook
    from openpyxl.styles import Font as _XlFont
    _HAS_XLSX = True
except Exception:
    _HAS_XLSX = False


def _xlsx_val(v):
    """Write a real number when the whole cell is numeric (so Excel sorts/filters it), else guarded
    text. A leading formula char is neutralised the same way as the CSV path (CWE-1236)."""
    s = "" if v is None else str(v)
    t = s.strip()
    if t:
        try:
            return int(t) if t.lstrip("-").isdigit() else float(t)
        except ValueError:
            pass
    if s[:1] in ("=", "@") or (s[:1] in ("+", "-") and len(s) > 1):    # bare +/- (strand) stays literal
        return "'" + s
    return s


def _export_xlsx(headers, rows, path):
    wb = _XlWorkbook()
    ws = wb.active
    ws.title = "TEagle"
    ws.append([str(h) for h in headers])
    for cell in ws[1]:
        cell.font = _XlFont(bold=True)
    for r in rows:
        ws.append([_xlsx_val(c) for c in r])
    ws.freeze_panes = "A2"                                # keep the header visible while scrolling
    wb.save(path)


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


def add_export_submenu(menu, headers, rows_fn, base, parent):
    """Attach a '⭳ Export table ▸ Excel/CSV/TSV' submenu; rows_fn is called at click time (current order).
    The submenu is parented to `menu` so C++ keeps it alive after this returns (a bare addMenu(str) can
    be GC'd out from under the menu)."""
    sub = QMenu("⭳ Export table", menu)
    for label, fmt in _available_formats():
        sub.addAction(label, lambda _=False, f=fmt: export_table(headers, rows_fn(), base, parent, fmt=f))
    menu.addMenu(sub)
    return sub


def export_table(headers, rows, base, parent=None, fmt=None):
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
    low = path.lower()
    if low.endswith(".xlsx") and _HAS_XLSX:
        _export_xlsx(headers, rows, path)
        return
    sep = "\t" if low.endswith(".tsv") else ","
    lines = [sep.join(_csv_escape(h, sep) for h in headers)]
    lines += [sep.join(_csv_escape(c, sep) for c in r) for r in rows]
    with open(path, "w", encoding="utf-8-sig", newline="") as f:      # BOM so Excel reads UTF-8
        f.write("\r\n".join(lines))


def save_fasta(fasta: str, base: str, parent=None):
    """Write a ready-made FASTA string to a user-chosen .fasta file."""
    path, _ = QFileDialog.getSaveFileName(parent, "Export FASTA", base + ".fasta",
                                          "FASTA (*.fasta *.fa *.txt)")
    if not path:
        return
    if not fasta.endswith("\n"):
        fasta += "\n"
    with open(path, "w", encoding="utf-8", newline="\n") as f:
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


class _Cell(QTableWidgetItem):
    """Table cell that sorts numerically when both cells parse as numbers, else alphabetically —
    so Score / E-value / aa / divergence / coords sort by value, not by string."""
    def __lt__(self, other):
        a, b = _sortkey(self.text()), _sortkey(other.text())
        if a is not None and b is not None:
            return a < b
        return self.text().casefold() < other.text().casefold()


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
        self.horizontalHeader().setDefaultAlignment(Qt.AlignCenter)   # centered headers
        self.setSortingEnabled(True)                                  # click a header to sort (numeric-aware via _Cell)
        self.doubleClicked.connect(lambda idx: self.row_activated.emit(self._orig(idx.row())))
        self._row_menu = None            # callable(orig_row_index)->list[(label, fn)]

    def _orig(self, visual_row):
        """Map a visual row (post-sort) back to the index it had in set_rows, so row menus /
        activation address the right data record no matter how the user sorted."""
        if visual_row < 0:
            return -1
        it = self.item(visual_row, 0)
        d = it.data(Qt.UserRole) if it else None
        return d if d is not None else visual_row

    def set_rows(self, rows):
        self.setSortingEnabled(False)                # never sort mid-insert (it scrambles rows)
        self.setRowCount(0)
        for i, r in enumerate(rows):
            self.insertRow(i)
            for j, c in enumerate(r):
                text = "" if c is None else str(c)
                item = _Cell(text)
                item.setTextAlignment(Qt.AlignCenter)    # centered cells, matching the headers
                item.setToolTip(text)                    # full value on hover, in case the cell elides
                if j == 0:
                    item.setData(Qt.UserRole, i)         # remember the original row index for menus
                self.setItem(i, j, item)
        self.horizontalHeader().setSortIndicator(-1, Qt.AscendingOrder)  # keep engine order until a header is clicked
        self.setSortingEnabled(True)
        self.resizeColumnsToContents()

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
        add_export_submenu(m, self._headers, self.rows_data, "TEagle_table", self)
        m.exec(self.viewport().mapToGlobal(pos))

    def _copy_row(self, row):
        if row < 0:
            return
        vals = [self.item(row, j).text() if self.item(row, j) else "" for j in range(self.columnCount())]
        QApplication.clipboard().setText("\t".join(vals))

    def rows_data(self):
        return [[self.item(i, j).text() if self.item(i, j) else "" for j in range(self.columnCount())]
                for i in range(self.rowCount())]
