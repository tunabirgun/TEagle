"""Render brand assets from the single-source mono mark (teagle-mark.svg).
Outputs, all recolored from the same currentColor SVG:
  installer/teagle.ico            multi-size app/installer icon, mid-teal on transparent
  report/figures/teagle-mark-dark.pdf/.png  dark-ink title-page logo (vector + raster)
Run: python app/native/assets/render_brand.py"""
import io, os
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
from PySide6.QtCore import QByteArray, QBuffer, QRectF, QSizeF, QMarginsF
from PySide6.QtGui import QImage, QPainter, QPdfWriter, QPageSize
from PySide6.QtSvg import QSvgRenderer
from PIL import Image, ImageFilter

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
SVG = os.path.join(HERE, "teagle-mark.svg")
ICON_TEAL = "#12B39A"                                    # OS chrome, reads on light + dark taskbars
DARK_INK  = "#141B21"                                    # report title page on white paper

with open(SVG, encoding="utf-8") as f:
    BASE = f.read()

def _renderer(color):
    return QSvgRenderer(QByteArray(BASE.replace("currentColor", color).encode("utf-8")))

def _png(color, h):                                      # aspect-preserved PNG bytes at height h
    r = _renderer(color); s = r.defaultSize()
    w = max(1, round(h * s.width() / s.height()))
    img = QImage(w, h, QImage.Format_ARGB32); img.fill(0)
    p = QPainter(img); r.render(p); p.end()
    buf = QByteArray(); b = QBuffer(buf); b.open(QBuffer.WriteOnly)
    img.save(b, "PNG"); b.close()
    return bytes(buf)

def _square(color, n, ss=12):                            # mark on an n x n transparent canvas, supersampled ss x
    big = Image.open(io.BytesIO(_png(color, round(n * ss * 0.94)))).convert("RGBA")
    canvas = Image.new("RGBA", (n * ss, n * ss), (0, 0, 0, 0))
    canvas.alpha_composite(big, ((n * ss - big.width) // 2, (n * ss - big.height) // 2))
    out = canvas.resize((n, n), Image.LANCZOS)           # downscale = clean anti-aliasing at small sizes
    if n <= 64:                                          # the detailed eagle/DNA trace blurs when reduced;
        out = out.filter(ImageFilter.UnsharpMask(radius=0.7, percent=110, threshold=0))  # restore edge crispness
    return out

def _bmp(color, w, h, fill, ss=4):                       # eagle centered on an opaque white canvas (Inno wizard image, no alpha)
    r = _renderer(color); s = r.defaultSize()
    scale = min(w * fill / s.width(), h * fill / s.height())
    eag = Image.open(io.BytesIO(_png(color, max(1, round(s.height() * scale * ss))))).convert("RGBA")
    canvas = Image.new("RGBA", (w * ss, h * ss), (255, 255, 255, 255))
    canvas.alpha_composite(eag, ((w * ss - eag.width) // 2, (h * ss - eag.height) // 2))
    return canvas.resize((w, h), Image.LANCZOS).convert("RGB")   # 24-bit BMP, supersampled for clean edges

def wizard_images():                                     # Inno Setup wizard eagle logo: top-right small + Welcome/Finish large
    instdir = os.path.join(ROOT, "installer")
    for n in (55, 83, 110, 138):                         # DPI ladder; Inno picks the closest to the current scaling
        p = os.path.join(instdir, "wizard-small" + ("" if n == 55 else f"-{n}") + ".bmp")
        # fill 0.72 → the eagle sits inset from the header corner with an equal ~14% margin on every side
        # (Inno anchors this image top-right, so the centered mark reads equally spaced from top and right)
        _bmp(ICON_TEAL, n, n, fill=0.72).save(p, "BMP"); print("wrote", p)
    for w, hh in ((164, 314), (246, 459), (328, 604), (410, 797)):
        p = os.path.join(instdir, "wizard-large" + ("" if w == 164 else f"-{w}") + ".bmp")
        _bmp(ICON_TEAL, w, hh, fill=0.66).save(p, "BMP"); print("wrote", p)

def main():
    sizes = [16, 20, 24, 32, 40, 48, 64, 96, 128, 256]  # standard Windows icon sizes; exact frame per DPI
    frames = [_square(ICON_TEAL, s) for s in sizes]
    ico = os.path.join(ROOT, "installer", "teagle.ico")
    # save from the LARGEST frame as base + append the rest, so every size embeds its own high-res
    # render. (Saving from the 16px frame made Pillow drop all larger sizes → a blurry upscaled icon.)
    frames[-1].save(ico, format="ICO", sizes=[(s, s) for s in sizes], append_images=frames[:-1])
    print("wrote", ico)

    wizard_images()                                      # Inno Setup wizard eagle logos (BMP)

    figdir = os.path.join(ROOT, "report", "figures"); os.makedirs(figdir, exist_ok=True)
    Image.open(io.BytesIO(_png(DARK_INK, 900))).save(os.path.join(figdir, "teagle-mark-dark.png"))
    print("wrote", os.path.join(figdir, "teagle-mark-dark.png"))

    r = _renderer(DARK_INK); s = r.defaultSize(); ar = s.width() / s.height()
    w_mm = 40.0; h_mm = w_mm / ar
    pdf = os.path.join(figdir, "teagle-mark-dark.pdf")
    pw = QPdfWriter(pdf); pw.setResolution(600)
    pw.setPageSize(QPageSize(QSizeF(w_mm, h_mm), QPageSize.Millimeter))
    pw.setPageMargins(QMarginsF(0, 0, 0, 0))
    p = QPainter(pw)
    r.render(p, QRectF(0, 0, p.viewport().width(), p.viewport().height())); p.end()
    print("wrote", pdf)

if __name__ == "__main__":
    main()
