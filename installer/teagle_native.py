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


def _ensure_std_streams():
    """A windowed build started from Explorer has no std handles, so sys.stdout/stderr are None and any
    write raises. Give them a UTF-8 log (devnull if it cannot be opened) before anything imports main."""
    if sys.stdout is not None and sys.stderr is not None:
        return
    try:
        d = os.path.join(os.environ.get("LOCALAPPDATA") or os.path.expanduser("~"), "TEagle", "logs")
        os.makedirs(d, exist_ok=True)
        path = os.path.join(d, "teagle.log")
        f = open(path, "w" if os.path.exists(path) and os.path.getsize(path) > 1_000_000 else "a",
                 encoding="utf-8", buffering=1)
    except OSError:
        f = open(os.devnull, "w", encoding="utf-8")
    _KEEP.append(f)
    if sys.stdout is None:
        sys.stdout = f
    if sys.stderr is None:
        sys.stderr = f


def _msgbox(text, title="TEagle", icon=0x10):
    """icon: MB_ICONERROR (default) or MB_ICONINFORMATION (0x40)."""
    try:
        import ctypes
        ctypes.windll.user32.MessageBoxW(0, str(text), title, icon)
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


def _raise_running_window(title="TEagle") -> bool:
    """Bring the already-running instance's main window forward. True if a Qt window with that title was found."""
    try:
        import ctypes
        u = ctypes.windll.user32
        hwnd = u.FindWindowW(None, title)
        if not hwnd:
            return False
        cls = ctypes.create_unicode_buffer(256)
        u.GetClassNameW(hwnd, cls, 256)
        if not cls.value.startswith("Qt"):                      # a folder or document that happens to be named TEagle
            return False
        if u.IsIconic(hwnd):
            u.ShowWindow(hwnd, 9)                               # SW_RESTORE
        u.SetForegroundWindow(hwnd)
        return True
    except Exception:
        return False


# KILL_ON_JOB_CLOSE reaps wsl.exe children; BREAKAWAY_OK lets a launcher that must outlive TEagle (a browser opened
# from a link) start itself with CREATE_BREAKAWAY_FROM_JOB.
_JOB_LIMIT_FLAGS = 0x2000 | 0x800


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
        info.BasicLimitInformation.LimitFlags = _JOB_LIMIT_FLAGS
        k.SetInformationJobObject(job, 9, ctypes.byref(info), ctypes.sizeof(info))   # ExtendedLimitInformation
        k.AssignProcessToJobObject(job, k.GetCurrentProcess())
        _KEEP.append(job)
    except Exception:
        pass                                                  # best-effort; not required for correctness


def main():
    if os.environ.get("TEAGLE_SELFTEST"):
        _ensure_std_streams()
        import main as native_main                            # app/native/main.py
        return native_main.selftest()
    if sys.platform == "win32":
        _enable_kill_on_close_job()                           # before any WSL child is spawned
        if _windows_single_instance():                        # another TEagle is already running
            if not _raise_running_window():
                _msgbox("TEagle is already running.", icon=0x40)
            return 0
    _ensure_std_streams()                                     # after the mutex: a second launch must not touch the running instance's log
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
