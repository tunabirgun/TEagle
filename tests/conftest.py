import sys, os

import pytest

ROOT = os.path.dirname(os.path.abspath(__file__))
# make the app's scientific core and the native UI modules importable, whichever test file runs first
sys.path.insert(0, os.path.abspath(os.path.join(ROOT, "..", "app", "native")))
sys.path.insert(0, os.path.abspath(os.path.join(ROOT, "..", "app", "backend")))


class _MemorySettings:
    """In-memory stand-in for QSettings("TEagle", "TEagle"), fresh for every test."""
    store = {}

    def __init__(self, *a, **k):
        pass

    def value(self, key, default=None, *a, **k):
        return self.store.get(key, default)

    def setValue(self, key, value):
        self.store[key] = value

    def contains(self, key):
        return key in self.store

    def remove(self, key):
        self.store.pop(key, None)

    def sync(self):
        pass


@pytest.fixture(autouse=True)
def _private_settings(monkeypatch):
    """Keep tests off the user's saved TEagle settings (Windows registry): a test that toggles the theme,
    UI scale or coordinate base must neither read the developer's choices nor overwrite them."""
    try:
        import main
    except ImportError:
        return
    _MemorySettings.store = {}
    monkeypatch.setattr(main, "QSettings", _MemorySettings)


@pytest.fixture(autouse=True)
def _dispose_top_level_widgets():
    """Delete top-level widgets a test leaves behind. A closed window is only hidden, and every new
    window applies the application-wide stylesheet, which repolishes every live widget; leftovers made
    each later window slower to build, so a one-session run grew to over an hour. A widget a
    module-scoped fixture still needs carries the dynamic property 'teagle_keep'."""
    yield
    try:
        from PySide6.QtWidgets import QApplication
        from PySide6.QtCore import QEvent
    except ImportError:
        return
    app = QApplication.instance()
    if app is None:
        return
    def kept(w):
        while w is not None:
            if w.property("teagle_keep"):
                return True
            w = w.parentWidget()
        return False

    stale = [w for w in app.topLevelWidgets() if not kept(w)]
    if not stale:
        return
    import gc
    for w in stale:
        w.deleteLater()            # not close(): a window with a long job running would ask a modal question
    QApplication.sendPostedEvents(None, QEvent.DeferredDelete)
    gc.collect()


@pytest.fixture(autouse=True)
def _hermetic_wsl(request, monkeypatch):
    """Present the machine as having no WSL distro, as on CI. Without this, a hermetic test on a machine
    with WSL installed probes the real distro (or runs a real annotation), and every window's thread pool
    waits on that probe at teardown. Tests marked `wsl` and the distro-resolution tests are exempt."""
    if request.node.get_closest_marker("wsl") or request.node.fspath.basename == "test_wsl.py":
        return
    from teagle_core import wsl
    monkeypatch.setattr(wsl, "resolve_distro", lambda *a, **k: None)
    monkeypatch.setattr(wsl, "_distro_cache", None, raising=False)
