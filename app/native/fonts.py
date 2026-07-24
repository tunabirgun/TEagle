"""Ship the UI fonts inside the app so text and the wordmark render identically on any machine,
whether or not the fonts are installed. Body/UI text is Roboto (Apache-2.0, see Roboto-LICENSE.txt);
sequences, accessions and numeric data stay in Cascadia Mono for column alignment; the brand
wordmark uses Cascadia Code Bold. Cascadia Code/Mono are SIL OFL 1.1 (see assets/fonts/OFL.txt).
One loader, used by main(), selftest(), the screenshot capture, and the banner generator — a
single source so screenshots and releases reflect exactly what an end user gets."""
import os
from PySide6.QtGui import QFontDatabase
from teagle_core import appdirs

_HERE = os.path.dirname(os.path.abspath(__file__))
# Roboto = body/UI text (Regular/Medium/Bold cover QSS weights 400/500/700).
# Cascadia Mono = sequences/accessions/data tables (Regular/SemiBold/Bold cover 400/600–650/700).
# Cascadia Code Bold = the wordmark generator only; the header ships it as frozen SVG paths.
FONT_FILES = ("Roboto-Regular.ttf", "Roboto-Medium.ttf", "Roboto-Bold.ttf",
              "CascadiaMono-Regular.ttf", "CascadiaMono-SemiBold.ttf",
              "CascadiaMono-Bold.ttf", "CascadiaCode-Bold.ttf")
BODY_FAMILY = "Roboto"
UI_FAMILY = "Cascadia Mono"
WORDMARK_FAMILY = "Cascadia Code"

def font_path(name: str) -> str:
    return appdirs.resource("native", "assets", "fonts", name) or \
        os.path.join(_HERE, "assets", "fonts", name)

def load_fonts() -> list[str]:
    """Register the bundled fonts with Qt. Returns the files that FAILED to load — a file missing
    from the frozen bundle (or unreadable) yields addApplicationFont == -1. An empty list means
    every font shipped. Requires a QGuiApplication to already exist."""
    missing = []
    for fn in FONT_FILES:
        p = font_path(fn)
        if not p or not os.path.exists(p) or QFontDatabase.addApplicationFont(p) == -1:
            missing.append(fn)
    return missing
