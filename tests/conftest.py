import sys, os

import pytest

ROOT = os.path.dirname(os.path.abspath(__file__))
# make the app's scientific core and the native UI modules importable, whichever test file runs first
sys.path.insert(0, os.path.abspath(os.path.join(ROOT, "..", "app", "native")))
sys.path.insert(0, os.path.abspath(os.path.join(ROOT, "..", "app", "backend")))


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
