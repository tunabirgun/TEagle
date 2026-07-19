"""Brand banner: eagle mark + two-tone TEagle wordmark, transparent, like the in-app header.
Fully self-contained vector SVG (text flattened to paths, no font dependency) + high-res PNG.
Dark variant for dark READMEs/backgrounds, light variant for white. Run: python render_banner.py"""
import io, os
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
from PySide6.QtCore import QByteArray, QRectF, Qt
from PySide6.QtGui import (QImage, QPainter, QColor, QFont, QFontDatabase,
                           QFontMetricsF, QPainterPath)
from PySide6.QtSvg import QSvgRenderer, QSvgGenerator
from PySide6.QtGui import QGuiApplication

_APP = QGuiApplication.instance() or QGuiApplication([])  # required for font/QPainterPath APIs

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
OUT = os.path.join(ROOT, "docs", "img")
MARK = open(os.path.join(HERE, "teagle-mark.svg"), encoding="utf-8").read()

VARIANTS = {
    "dark":  dict(mark="#33D6B8", te="#E6EDF1", agle="#33D6B8"),   # on dark backgrounds
    "light": dict(mark="#0E9E86", te="#141B21", agle="#0E9E86"),   # on white backgrounds
}
H = 260.0                                                # mark height (layout units)
PAD = 24.0                                               # transparent margin around content

_FAMILY = None
def _font():
    global _FAMILY
    if _FAMILY is None:                                  # headless offscreen ships no fonts; load a real TTF
        winf = os.path.join(os.environ.get("WINDIR", "C:/Windows"), "Fonts")
        for fn in ("CascadiaCode.ttf", "CascadiaCodeNF.ttf", "consolab.ttf", "arialbd.ttf"):
            fid = QFontDatabase.addApplicationFont(os.path.join(winf, fn))
            if fid != -1:
                _FAMILY = QFontDatabase.applicationFontFamilies(fid)[0]; break
        else:
            _FAMILY = QFont().defaultFamily()
    f = QFont(_FAMILY); f.setStyleHint(QFont.Monospace)
    f.setBold(True); f.setWeight(QFont.Bold); f.setPixelSize(int(H * 0.52))
    return f

def _path_to_d(path, dx=0.0, dy=0.0):                    # QPainterPath -> SVG path data (freezes the exact glyphs)
    parts = []; started = False; i = 0; n = path.elementCount()
    while i < n:
        e = path.elementAt(i); x = e.x + dx; y = e.y + dy
        if e.type == QPainterPath.MoveToElement:
            if started: parts.append("Z")
            parts.append(f"M{x:.2f} {y:.2f}"); started = True; i += 1
        elif e.type == QPainterPath.LineToElement:
            parts.append(f"L{x:.2f} {y:.2f}"); i += 1
        elif e.type == QPainterPath.CurveToElement:
            e2 = path.elementAt(i + 1); e3 = path.elementAt(i + 2)
            parts.append(f"C{x:.2f} {y:.2f} {e2.x+dx:.2f} {e2.y+dy:.2f} {e3.x+dx:.2f} {e3.y+dy:.2f}"); i += 3
        else:
            i += 1
    if started: parts.append("Z")
    return "".join(parts)

def _layout():
    r = QSvgRenderer(QByteArray(MARK.encode())); s = r.defaultSize()
    mw = H * s.width() / s.height()
    font = _font(); tracking = font.pixelSize() * 0.14
    # two separate paths so TE and AGLE can take different colors
    fm = QFontMetricsF(font); x = 0.0
    te = QPainterPath(); agle = QPainterPath()
    for i, ch in enumerate("TEAGLE"):
        p = QPainterPath(); p.addText(x, 0.0, font, ch)
        (te if i < 2 else agle).addPath(p)
        x += fm.horizontalAdvance(ch) + tracking
    text_w = x - tracking
    gap = H * 0.16
    total_w = mw + gap + text_w
    return r, mw, font, te, agle, text_w, gap, total_w, fm

def _paint(p, colors):
    r, mw, font, te, agle, text_w, gap, total_w, fm = _layout()
    p.setRenderHint(QPainter.Antialiasing, True)
    p.setRenderHint(QPainter.TextAntialiasing, True)
    # eagle mark, recolored, left-aligned, full height
    r2 = QSvgRenderer(QByteArray(MARK.replace("currentColor", colors["mark"]).encode()))
    r2.render(p, QRectF(PAD, PAD, mw, H))
    # wordmark: cap-centered on the mark's vertical midline
    baseline = PAD + H / 2 + fm.capHeight() / 2
    p.translate(PAD + mw + gap, baseline)
    p.setPen(Qt.NoPen)
    p.setBrush(QColor(colors["te"]));   p.drawPath(te)
    p.setBrush(QColor(colors["agle"])); p.drawPath(agle)
    return total_w + 2 * PAD, H + 2 * PAD

def _dims():
    r, mw, font, te, agle, text_w, gap, total_w, fm = _layout()
    return total_w + 2 * PAD, H + 2 * PAD

def render(variant, scale=4):
    colors = VARIANTS[variant]
    W, Hh = _dims()
    os.makedirs(OUT, exist_ok=True)
    # PNG (high-res, transparent)
    img = QImage(int(W * scale), int(Hh * scale), QImage.Format_ARGB32); img.fill(0)
    p = QPainter(img); p.scale(scale, scale); _paint(p, colors); p.end()
    png = os.path.join(OUT, f"teagle-banner-{variant}.png"); img.save(png)
    # SVG (vector, self-contained)
    svg = os.path.join(OUT, f"teagle-banner-{variant}.svg")
    gen = QSvgGenerator(); gen.setFileName(svg)
    gen.setSize(img.size()); gen.setViewBox(QRectF(0, 0, W, Hh))
    gen.setTitle("TEagle"); gen.setDescription("TEagle brand banner")
    p = QPainter(gen); _paint(p, colors); p.end()
    print("wrote", png, "and", svg, f"({int(W)}x{int(Hh)})")

def write_wordmark_asset():
    """Freeze the notched Cascadia 'TEAGLE' into a static two-tone SVG (tokens {TE}/{AGLE})
    for the in-app header — no font dependency, identical to the banner wordmark."""
    r, mw, font, te, agle, text_w, gap, total_w, fm = _layout()
    both = QPainterPath(); both.addPath(te); both.addPath(agle)
    br = both.boundingRect()
    dx, dy = -br.x(), -br.y()
    te_d = _path_to_d(te, dx, dy); ag_d = _path_to_d(agle, dx, dy)
    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {br.width():.2f} {br.height():.2f}" '
           f'role="img" aria-label="TEagle">'
           f'<path d="{te_d}" fill="{{TE}}"/><path d="{ag_d}" fill="{{AGLE}}"/></svg>\n')
    dest = os.path.join(HERE, "teagle-wordmark.svg")
    open(dest, "w", encoding="utf-8").write(svg)
    print("wrote", dest, f"({br.width():.0f}x{br.height():.0f})")

def render_social(w=1280, h=640, bg="#0B0F14"):
    """GitHub social-preview card (2:1). Dark background so it reads on any embed;
    dark-variant banner centered. Transparent cards break on light OG surfaces."""
    from PySide6.QtGui import QColor
    W, Hh = _dims(); scale = min(w / (W * 1.35), h / (Hh * 1.9))
    img = QImage(w, h, QImage.Format_ARGB32); img.fill(QColor(bg))
    p = QPainter(img)
    p.translate((w - W * scale) / 2, (h - Hh * scale) / 2); p.scale(scale, scale)
    _paint(p, VARIANTS["dark"]); p.end()
    os.makedirs(OUT, exist_ok=True)
    dest = os.path.join(OUT, "teagle-social-preview.png"); img.save(dest)
    print("wrote", dest, f"({w}x{h})")

if __name__ == "__main__":
    write_wordmark_asset()
    for v in VARIANTS:
        render(v)
    render_social()
