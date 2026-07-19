"""Where the installed app keeps its writable state and finds its bundled resources.

From source, writable state lives in <project>/.teagle (unchanged) and resources are
read from the source tree. When packaged as a one-file/one-dir app (PyInstaller sets
sys.frozen), writable state moves to %LOCALAPPDATA%/TEagle so nothing is written into
Program Files, and bundled resources are read from the extraction dir (sys._MEIPASS)."""
from __future__ import annotations
import os, sys


def _frozen() -> bool:
    return bool(getattr(sys, "frozen", False))


def user_data_dir() -> str:
    """User-writable directory for cache, timings, and the environment signature."""
    if _frozen():
        base = os.environ.get("LOCALAPPDATA") or os.path.expanduser("~")
        d = os.path.join(base, "TEagle")
    else:
        d = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", ".teagle"))
    try:
        os.makedirs(d, exist_ok=True)
    except Exception:
        pass
    return d


def resource(*parts: str) -> str | None:
    """Absolute path to a bundled read-only resource (e.g. web/, data/) when frozen; else None."""
    if _frozen():
        return os.path.join(sys._MEIPASS, *parts)
    return None
