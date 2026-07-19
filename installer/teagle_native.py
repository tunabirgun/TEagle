"""TEagle native desktop launcher (packaged one-click build).

Opens the PySide6/Qt application window directly — no browser, no WebView2, no local HTTP server.
The engine and every dependency (PySide6, primer3, pyhmmer, the Pfam HMM profiles) are bundled in
the executable; a wetlab user installs one thing and runs it. The Dfam family / splice features are
optional and install with one click from inside the app; nothing here is required for the core science.

Windows hardening: a named mutex enforces single-instance, and a KILL_ON_JOB_CLOSE Job Object binds
every child (any WSL/RepeatMasker/minimap2 subprocess the app spawns) so closing TEagle — or an
in-place upgrade terminating it — never leaves an orphaned process tree.
Set TEAGLE_SELFTEST=1 to run the headless bundle self-test (imports + QtSvg + a real analysis) and exit.
"""
import os, sys

# make the bundled backend + native package importable (frozen: PyInstaller collects them;
# source: add the two source dirs so `python installer/teagle_native.py` also works)
if not getattr(sys, "frozen", False):
    _root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    sys.path.insert(0, os.path.join(_root, "app", "backend"))
    sys.path.insert(0, os.path.join(_root, "app", "native"))

_KEEP = []            # keep OS handles alive for the whole process lifetime


def _msgbox(text, title="TEagle"):
    try:
        import ctypes
        ctypes.windll.user32.MessageBoxW(0, str(text), title, 0x10)   # MB_ICONERROR
    except Exception:
        sys.stderr.write(f"{title}: {text}\n")


def _windows_single_instance() -> bool:
    """Named mutex: a second launch detects the first instance and returns True (should bow out)."""
    try:
        import ctypes
        k = ctypes.windll.kernel32
        m = k.CreateMutexW(None, False, "Global\\TEagle_native_single_instance")
        _KEEP.append(m)
        return k.GetLastError() == 183               # ERROR_ALREADY_EXISTS
    except Exception:
        return False


def _enable_kill_on_close_job():
    """Bind this process (and its children) to a KILL_ON_JOB_CLOSE job. Non-fatal if it fails
    (e.g. already inside a non-nestable job). Must run before any child subprocess is spawned."""
    try:
        import ctypes
        import ctypes.wintypes as wt

        class _BASIC(ctypes.Structure):
            _fields_ = [("PerProcessUserTimeLimit", ctypes.c_longlong),
                        ("PerJobUserTimeLimit", ctypes.c_longlong),
                        ("LimitFlags", ctypes.c_uint32),
                        ("MinimumWorkingSetSize", ctypes.c_size_t),
                        ("MaximumWorkingSetSize", ctypes.c_size_t),
                        ("ActiveProcessLimit", ctypes.c_uint32),
                        ("Affinity", ctypes.c_size_t),
                        ("PriorityClass", ctypes.c_uint32),
                        ("SchedulingClass", ctypes.c_uint32)]

        class _IOC(ctypes.Structure):
            _fields_ = [("r", ctypes.c_ulonglong), ("w", ctypes.c_ulonglong), ("o", ctypes.c_ulonglong),
                        ("rt", ctypes.c_ulonglong), ("wt", ctypes.c_ulonglong), ("ot", ctypes.c_ulonglong)]

        class _EXT(ctypes.Structure):
            _fields_ = [("BasicLimitInformation", _BASIC), ("IoInfo", _IOC),
                        ("ProcessMemoryLimit", ctypes.c_size_t), ("JobMemoryLimit", ctypes.c_size_t),
                        ("PeakProcessMemoryUsed", ctypes.c_size_t), ("PeakJobMemoryUsed", ctypes.c_size_t)]

        k = ctypes.WinDLL("kernel32", use_last_error=True)
        k.CreateJobObjectW.restype = wt.HANDLE
        k.AssignProcessToJobObject.argtypes = [wt.HANDLE, wt.HANDLE]
        job = k.CreateJobObjectW(None, None)
        if not job:
            return
        info = _EXT()
        info.BasicLimitInformation.LimitFlags = 0x2000        # JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE
        k.SetInformationJobObject(job, 9, ctypes.byref(info), ctypes.sizeof(info))   # ExtendedLimitInformation
        k.AssignProcessToJobObject(job, k.GetCurrentProcess())
        _KEEP.append(job)
    except Exception:
        pass                                                  # best-effort; not required for correctness


def main():
    if os.environ.get("TEAGLE_SELFTEST"):
        import main as native_main                            # app/native/main.py
        return native_main.selftest()
    if sys.platform == "win32":
        _enable_kill_on_close_job()                           # before any WSL child is spawned
        if _windows_single_instance():
            return 0                                          # another TEagle is already running
    try:
        import main as native_main
        return native_main.main()
    except Exception as e:
        import traceback
        traceback.print_exc()
        _msgbox(f"TEagle could not start:\n\n{type(e).__name__}: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
