"""TEagle environment check + first-run installer.

On first run — and after any install or upgrade (detected via a signature over
requirements.txt + the app version) — this verifies Python and the pinned Python
packages, installs/repairs them via pip, and records the verified state so it does
not reinstall every launch. It also probes optional backends (WebView2, WSL2) that
gate advanced features. Nothing fails silently: a failed install is reported with
the pip error, never treated as success.
"""
from __future__ import annotations
import sys, os, json, subprocess, hashlib, shutil
import importlib.metadata as md
from teagle_core import appdirs

FROZEN = bool(getattr(sys, "frozen", False))                 # packaged one-click build (deps bundled)
# no console-window flash when the windowed GUI shells out (Windows only; 0 elsewhere)
_NO_WINDOW = getattr(subprocess, "CREATE_NO_WINDOW", 0) if os.name == "nt" else 0
HERE = os.path.dirname(os.path.abspath(__file__))            # app/backend
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))       # project root
REQ = appdirs.resource("requirements.txt") or os.path.join(HERE, "requirements.txt")
STATE_DIR = appdirs.user_data_dir()                          # %LOCALAPPDATA%/TEagle when installed
STATE = os.path.join(STATE_DIR, "env_state.json")
MIN_PY = (3, 10)


def _app_version() -> str:
    try:
        sys.path.insert(0, HERE)
        from teagle_core import __version__
        return __version__
    except Exception:
        return "unknown"


def parse_requirements() -> list[tuple[str, str]]:
    reqs = []
    if not os.path.isfile(REQ):
        return reqs
    for line in open(REQ, encoding="utf-8"):
        line = line.split("#")[0].strip()
        if not line:
            continue
        if "==" in line:
            name, ver = line.split("==", 1)
            reqs.append((name.strip(), ver.strip()))
        else:
            reqs.append((line.strip(), ""))
    return reqs


def _installed(dist: str):
    try:
        return md.version(dist)
    except md.PackageNotFoundError:
        return None
    except Exception:
        return None


def _signature(reqs) -> str:
    raw = _app_version() + "|" + ";".join(f"{n}=={v}" for n, v in reqs)
    return hashlib.sha256(raw.encode()).hexdigest()[:16]


def _read_state() -> dict:
    try:
        return json.load(open(STATE, encoding="utf-8"))
    except Exception:
        return {}


def _write_state(sig: str, pkgs):
    os.makedirs(STATE_DIR, exist_ok=True)
    json.dump({"signature": sig, "app_version": _app_version(),
               "packages": {n: v for n, v, _ in pkgs}},
              open(STATE, "w", encoding="utf-8"), indent=2)


def probe_backends() -> dict:
    """Optional backends that gate advanced features (report, never fail)."""
    wv = os.path.isdir(r"C:\Program Files (x86)\Microsoft\EdgeWebView\Application") or \
         os.path.isdir(r"C:\Program Files\Microsoft\EdgeWebView\Application")
    wsl = None
    exe = shutil.which("wsl")
    if exe:
        try:
            out = subprocess.run([exe, "-l", "-v"], capture_output=True, timeout=8,
                                 text=True, errors="ignore", creationflags=_NO_WINDOW)
            txt = (out.stdout or "").replace("\x00", "")
            if "Running" in txt or "Stopped" in txt:
                wsl = "available"
            else:                                            # wsl.exe present but no distro registered
                wsl = "no distro"
        except Exception:
            wsl = "error"
    return {"webview2": bool(wv), "wsl2": wsl or "not found"}


_IMPORT_NAME = {"primer3-py": "primer3", "pyhmmer": "pyhmmer", "PySide6": "PySide6",
                "ViennaRNA": "RNA"}   # verify they truly import (dist name -> module; ViennaRNA imports as RNA)
_SOURCE_ONLY = set()                                            # every pinned dep ships in the frozen native build


def check() -> dict:
    reqs = parse_requirements()
    pkgs = []            # (name, required, installed)
    all_ok = True
    import_fail = []
    for name, ver in reqs:
        if FROZEN and name in _SOURCE_ONLY:                      # not part of the packaged build -> don't flag "install needed"
            continue
        inst = _installed(name)
        ok = (inst is not None) and (ver == "" or inst == ver)
        if ok and name in _IMPORT_NAME:                      # metadata present != importable
            try:
                __import__(_IMPORT_NAME[name])
            except Exception:
                ok = False
                import_fail.append(name)
        if not ok:
            all_ok = False
        pkgs.append((name, ver, inst))
    sig = _signature(reqs)
    last = _read_state().get("signature")
    py_ok = sys.version_info[:2] >= MIN_PY
    return {
        "python": ".".join(map(str, sys.version_info[:3])),
        "python_ok": py_ok,
        "app_version": _app_version(),
        "packages": [{"name": n, "required": r or "(any)", "installed": i, "ok": (i is not None and (r == "" or i == r))}
                     for n, r, i in pkgs],
        "packages_ok": all_ok,
        "importFailures": import_fail,
        "signature": sig,
        "signature_matches": (sig == last),
        "first_run": last is None,
        "needs_install": (not all_ok) or (sig != last),
        "backends": probe_backends(),
    }


def ensure(auto_install: bool = True, verbose: bool = True) -> dict:
    """Full first-run flow: check, install if needed, re-verify, persist state."""
    rep = check()
    rep["installed_now"] = []
    rep["error"] = None
    if not rep["python_ok"]:
        rep["error"] = f"Python {rep['python']} < required {'.'.join(map(str, MIN_PY))}"
        return rep
    if rep["needs_install"] and auto_install and not FROZEN:   # packaged build ships its deps; nothing to pip-install
        if verbose:
            why = "first run" if rep["first_run"] else ("requirements/version changed"
                   if not rep["signature_matches"] else "missing/outdated packages")
            print(f"[env] {why} — installing pinned dependencies…")
        proc = subprocess.run([sys.executable, "-m", "pip", "install", "-r", REQ],
                              capture_output=True, text=True, errors="ignore", creationflags=_NO_WINDOW)
        if proc.returncode != 0:
            rep["error"] = "pip install failed:\n" + (proc.stderr or proc.stdout or "")[-1500:]
            if verbose:
                print("[env] ERROR — dependency install failed. See message.")
            rep = {**check(), "installed_now": [], "error": rep["error"], "backends": rep["backends"]}
            return rep
        import importlib
        importlib.invalidate_caches()               # so a same-process import sees just-installed packages
        after = check()
        after["installed_now"] = [p["name"] for p in after["packages"] if p["ok"]]
        after["error"] = None
        if after["packages_ok"]:
            _write_state(after["signature"], [(p["name"], p["installed"], True) for p in after["packages"]])
            if verbose:
                print("[env] dependencies verified and recorded.")
        else:
            after["error"] = "some packages still unsatisfied after install"
        return after
    if rep["packages_ok"] and rep["signature"] != _read_state().get("signature"):
        _write_state(rep["signature"], [(p["name"], p["installed"], True) for p in rep["packages"]])
    if verbose and not rep["needs_install"]:
        print(f"[env] up to date (sig {rep['signature']}, {rep['app_version']}).")
    return rep


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="TEagle environment check")
    ap.add_argument("--no-install", action="store_true", help="check only, do not install")
    a = ap.parse_args()
    r = ensure(auto_install=not a.no_install)
    print(json.dumps({k: r[k] for k in ("python", "python_ok", "app_version", "packages_ok",
          "needs_install", "first_run", "installed_now", "error", "backends")}, indent=2))
    sys.exit(0 if (r["python_ok"] and r["packages_ok"] and not r["error"]) else 1)
