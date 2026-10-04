"""Deferred callbacks must not outlive the window that scheduled them."""
import gc
import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PySide6")
from PySide6.QtCore import QEvent
from PySide6.QtWidgets import QApplication

_app = QApplication.instance() or QApplication([])


def test_installer_deleted_before_its_deferred_refresh_never_runs_it(monkeypatch):
    import install_dialog
    ran = []
    monkeypatch.setattr(install_dialog.InstallDialog, "_refresh", lambda self: ran.append(self))
    dlg = install_dialog.InstallDialog()
    dlg.deleteLater()
    QApplication.sendPostedEvents(None, QEvent.DeferredDelete)
    del dlg
    gc.collect()
    _app.processEvents()
    assert not ran
