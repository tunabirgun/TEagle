"""TEagle launcher — the one entry point (source / dev).

    python app/teagle.py            # first-run env check, then the native PySide6 window
    python app/teagle.py --check    # run the environment check and exit
    python app/teagle.py --selftest # headless bundle self-test (imports + QtSvg + a real analysis)
    python app/teagle.py --server   # legacy web UI over a local browser server (no Qt window)

On first run, and after any install or upgrade, it verifies Python + the pinned
dependencies and installs them before launching. It never starts on a broken
environment and never reports a failed install as success. The packaged one-click
build uses installer/teagle_native.py instead (deps already bundled).
"""
import os, sys, json, argparse, threading, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "backend"))
sys.path.insert(0, os.path.join(HERE, "native"))
import envcheck


def main() -> int:
    ap = argparse.ArgumentParser(description="TEagle launcher")
    ap.add_argument("--server", action="store_true", help="legacy web UI over a local browser server")
    ap.add_argument("--check", action="store_true", help="run environment check and exit")
    ap.add_argument("--selftest", action="store_true", help="headless bundle self-test and exit")
    ap.add_argument("--no-install", action="store_true", help="do not auto-install missing deps")
    ap.add_argument("--port", type=int, default=8765)
    a = ap.parse_args()

    if a.selftest:
        os.environ["TEAGLE_SELFTEST"] = "1"
        import main as native_main
        return native_main.selftest()

    print("TEagle — checking environment…")
    rep = envcheck.ensure(auto_install=not a.no_install)

    if a.check:
        print(json.dumps({k: rep[k] for k in ("python", "python_ok", "app_version",
              "packages_ok", "needs_install", "first_run", "installed_now", "error", "backends")}, indent=2))
        return 0 if (rep["python_ok"] and rep["packages_ok"] and not rep["error"]) else 1

    if rep["error"] or not rep["python_ok"] or not rep["packages_ok"]:
        print("\nx Cannot start — environment not ready:")
        print("  " + (rep["error"] or f"Python {rep['python']} unsupported / packages missing"))
        print("  Fix: python -m pip install -r app/backend/requirements.txt")
        return 1

    if a.server:                                     # legacy browser path (kept for headless / remote use)
        import server
        server.serve(port=a.port)
        return 0

    import main as native_main                        # app/native/main.py — the native PySide6 app
    return native_main.main()


if __name__ == "__main__":
    sys.exit(main())
