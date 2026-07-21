"""Threaded engine adapter for the Qt UI. Runs blocking engine.run_* calls off the GUI thread
and routes the outcome through the three-tier taxonomy (obs 4702):

  done(key, result_dict)        -> success
  user_error(key, message)      -> BadRequest: a correctable input problem, shown inline
  failed(key, message, trace)   -> unexpected fault: error banner + traceback to the log

`key` lets a panel match a result to the request it issued (and ignore stale ones)."""
from __future__ import annotations
import os, sys, traceback

from PySide6.QtCore import QObject, QRunnable, QThreadPool, Signal

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "backend"))
import engine
from engine import BadRequest


class _Signals(QObject):
    done = Signal(str, object)
    user_error = Signal(str, str)
    failed = Signal(str, str, str)


class _Job(QRunnable):
    def __init__(self, key, fn, signals):
        super().__init__()
        self.key, self.fn, self.signals = key, fn, signals

    def _emit(self, signal, *args):
        try:
            signal.emit(*args)
        except RuntimeError:
            pass                                                 # window closed mid-job: receiver gone, drop the result

    def run(self):
        try:
            res = self.fn()
        except BadRequest as e:
            self._emit(self.signals.user_error, self.key, str(e))
        except Exception as e:                                    # unexpected — surface + log full trace
            self._emit(self.signals.failed, self.key, f"{type(e).__name__}: {e}", traceback.format_exc())
        else:
            self._emit(self.signals.done, self.key, res)


class Engine(QObject):
    """Submit engine operations by key; connect to done / user_error / failed."""
    done = Signal(str, object)
    user_error = Signal(str, str)
    failed = Signal(str, str, str)

    # op name -> engine function taking the request body dict
    OPS = {
        "health": lambda b: engine.run_health(),
        "env": lambda b: engine.run_env(),
        "wsl_status": lambda b: engine.run_wsl_status(),
        "wsl_install": lambda b: engine.run_wsl_install(),
        "wsl_install_log": lambda b: engine.run_wsl_install_log(),
        "wsl_install_wsl2": lambda b: engine.run_wsl_install_wsl2(),
        "wsl2_install_log": lambda b: engine.run_wsl2_install_log(),
        "wsl_components": lambda b: engine.run_wsl_components(),
        "wsl_repair": engine.run_wsl_repair,
        "wsl_integrity": lambda b: engine.run_wsl_integrity(),
        "fetch": engine.run_fetch,
        "fetch_coords": engine.run_fetch_coords,
        "analyze": engine.run_analyze,
        "eta": engine.run_eta,
        "annotate": engine.run_annotate,
        "splice": engine.run_splice,
        "primers": engine.run_primers,
        "pcr": engine.run_pcr,
    }

    def __init__(self, parent=None, max_threads: int | None = None):
        super().__init__(parent)
        self.pool = QThreadPool(self)
        if max_threads:
            self.pool.setMaxThreadCount(max_threads)
        self._sig = _Signals()
        self._sig.done.connect(self.done)
        self._sig.user_error.connect(self.user_error)
        self._sig.failed.connect(self.failed)

    def submit(self, op: str, body: dict | None = None, key: str | None = None):
        """Run engine op `op` with request `body`. `key` (default = op) tags the result."""
        fn = self.OPS[op]
        b = body or {}
        self.pool.start(_Job(key or op, lambda: fn(b), self._sig))

    def wait(self, msecs: int = -1) -> bool:
        """Block until all queued jobs finish (used by headless tests)."""
        return self.pool.waitForDone(msecs)
