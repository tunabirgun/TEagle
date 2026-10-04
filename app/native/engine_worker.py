"""Threaded engine adapter for the Qt UI. Runs blocking engine.run_* calls off the GUI thread
and routes the outcome through the three-tier taxonomy (obs 4702):

  done(key, result_dict)        -> success
  user_error(key, message)      -> BadRequest: a correctable input problem, shown inline
  failed(key, message, trace)   -> unexpected fault: error banner + traceback to the log

`key` lets a panel match a result to the request it issued (and ignore stale ones)."""
from __future__ import annotations
import os, sys, threading, traceback

from PySide6.QtCore import QObject, QThread, Signal

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "backend"))
import engine
from engine import BadRequest


class _Signals(QObject):
    done = Signal(str, object)
    user_error = Signal(str, str)
    failed = Signal(str, str, str)


class _Job:
    def __init__(self, key, fn, signals):
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
        "genome_pcr": engine.run_genome_pcr,
        "genome_prepare": engine.run_genome_prepare,
        "genome_prepare_log": engine.run_genome_prepare_log,
        "genome_list": engine.run_genome_list,
        "genome_remove": engine.run_genome_remove,
        "add_custom_assembly": engine.run_add_custom_assembly,
        "genome_annotate": engine.run_genome_annotate,
        "genome_annotate_log": engine.run_genome_annotate_log,
        "genome_annotate_status": engine.run_genome_annotate_status,
        "genome_annotate_reset": engine.run_genome_annotate_reset,
        "annotate_budget": engine.run_annotate_budget,
    }

    def __init__(self, parent=None, max_threads: int | None = None):
        super().__init__(parent)
        # Daemon threads, not a QThreadPool owned by this object: a pool's destructor waits for every running
        # job, and a whole-genome annotation or an integrity probe can run for minutes to hours inside WSL.
        # Destroying an Engine (closing the backend dialog) or exiting the application must never wait on one.
        self._slots = threading.Semaphore(max_threads or max(2, QThread.idealThreadCount()))
        self._idle = threading.Condition()
        self._active = 0
        self._sig = _Signals()
        self._sig.done.connect(self.done)
        self._sig.user_error.connect(self.user_error)
        self._sig.failed.connect(self.failed)

    def submit(self, op: str, body: dict | None = None, key: str | None = None):
        """Run engine op `op` with request `body`. `key` (default = op) tags the result."""
        fn = self.OPS[op]
        b = body or {}
        job = _Job(key or op, lambda: fn(b), self._sig)
        with self._idle:
            self._active += 1
        threading.Thread(target=self._run, args=(job,), daemon=True, name=f"teagle-{op}").start()

    def _run(self, job):
        try:
            with self._slots:
                job.run()
        finally:
            with self._idle:
                self._active -= 1
                self._idle.notify_all()

    def wait(self, msecs: int = -1) -> bool:
        """Block until all submitted jobs finish (used by headless tests)."""
        with self._idle:
            return self._idle.wait_for(lambda: self._active == 0, None if msecs < 0 else msecs / 1000)
