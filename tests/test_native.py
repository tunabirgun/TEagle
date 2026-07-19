"""Native (PySide6) app regression tests. Run headless via the offscreen Qt platform. Cover the
figure builders, the engine worker's three-tier taxonomy, and the analyze -> design -> in-silico PCR
workflow driven through the same in-process engine the GUI uses. Skipped if PySide6 is unavailable."""
import os, sys, time

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
import pytest

pytest.importorskip("PySide6")
_NATIVE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "app", "native")
sys.path.insert(0, _NATIVE)

from PySide6.QtWidgets import QApplication
from PySide6.QtCore import QEventLoop, QTimer, QByteArray
from PySide6.QtGui import QImage, QPainter
from PySide6.QtSvg import QSvgRenderer

_app = QApplication.instance() or QApplication([])


def _spin(pred, timeout=30):
    loop = QEventLoop(); t = QTimer(); t.setInterval(25); t0 = time.time()
    t.timeout.connect(lambda: (pred() or time.time() - t0 > timeout) and loop.quit())
    t.start(); loop.exec()
    return pred()


def _renders(svg):
    r = QSvgRenderer(QByteArray(svg.encode("utf-8")))
    img = QImage(300, 200, QImage.Format_ARGB32); img.fill(0)
    p = QPainter(img); r.render(p); p.end()
    return r.isValid() and any(img.pixelColor(x, y).alpha() > 0
                               for x in range(0, 300, 15) for y in range(0, 200, 15))


# ---------- figure builders ----------
def test_gel_svg_renders():
    import figures
    svg = figures.svg_gel({"lanes": [{"label": "P1", "amplicons": [
        {"length": 250, "on_target": True, "source": "x"}, {"length": 80, "on_target": False, "source": "x"}]}]}, "dark")
    assert svg.startswith("<svg") and _renders(svg)


def test_gel_small_amplicon_axis_floor():
    # a 75 bp band must resolve — the log axis floor tracks the smallest band (obs 4680)
    import figures
    svg = figures.svg_gel({"lanes": [{"label": "P1", "amplicons": [{"length": 75, "on_target": True, "source": "x"}]}]}, "dark")
    assert "75 bp" in svg


def test_genome_svg_renders():
    import figures
    rec = {"composition": {"length": 5000},
           "structural": [{"type": "LTR (x)", "five_prime": [100, 300], "three_prime": [4700, 4900]}],
           "domains": [{"domain": "RT", "label": "RT", "nt": [1000, 2000], "score": 120}],
           "orfs": [{"start": 900, "end": 3000, "strand": "+", "frame": 0, "length_aa": 700}]}
    model = figures.gv_tracks_from_rec(rec)
    assert len(model["tracks"]) == 3
    assert _renders(figures.svg_genome(model, {"start": 0, "end": 5000}, 620, "dark"))


def test_png_export_has_alpha(tmp_path):
    import figures
    from widgets import render_png
    svg = figures.svg_gel({"lanes": [{"label": "P1", "amplicons": [{"length": 200, "on_target": True, "source": "x"}]}]}, "transparent")
    out = str(tmp_path / "gel.png")
    render_png(svg, out, scale=2)
    img = QImage(out)
    assert not img.isNull() and img.hasAlphaChannel()


# ---------- engine worker taxonomy ----------
def test_worker_success_and_user_error():
    from engine_worker import Engine
    eng = Engine()
    got = {}
    eng.done.connect(lambda k, r: got.__setitem__(k, ("done", r)))
    eng.user_error.connect(lambda k, m: got.__setitem__(k, ("user_error", m)))
    eng.failed.connect(lambda k, m, t: got.__setitem__(k, ("failed", m)))
    eng.submit("analyze", {"sequence": "ACGTACGTACGTACGTACGT"}, key="ok")
    eng.submit("pcr", {"sequence": "ACGT" * 40}, key="bad")     # missing primers -> BadRequest
    _spin(lambda: "ok" in got and "bad" in got)
    assert got["ok"][0] == "done" and "records" in got["ok"][1]
    assert got["bad"][0] == "user_error"


# ---------- full workflow ----------
@pytest.fixture
def win():
    import main
    w = main.MainWindow()
    yield w
    w.close()


def test_analyze_to_pcr_workflow(win):
    win._load_sample()
    win._run_analysis()
    assert _spin(lambda: win.state.get("analyzed_clean")), "analyze did not complete"
    assert win.designBtn.isEnabled()
    assert win.mValid.text() == "valid"

    done = {"p": False}
    orig = win._on_primers
    win._on_primers = lambda d: (orig(d), done.__setitem__("p", True))
    win._design()
    assert _spin(lambda: done["p"]), "primer design did not complete"
    cands = win.state.get("candidates", [])
    assert cands, "no primer candidates"

    win._pcr_stage_all()
    assert len(win.state["pcrPairs"]) == len(cands)
    assert win.runPcrBtn.isEnabled()

    orig_pcr = win._on_pcr
    finished = {"v": False}
    def wrap(k, d):
        orig_pcr(k, d)
        if win._pcr_run and all(r is not None for r in win._pcr_run["results"]):
            finished["v"] = True
    win._on_pcr = wrap
    win._run_pcr()
    assert _spin(lambda: finished["v"], 40), "in-silico PCR did not complete"
    assert len(win.state["lastPcr"]["lanes"]) == len(cands)
    assert win.card_prov.bodylay.count() >= 1        # provenance rendered


def test_staleness_guard_blocks_stale_primer(win):
    win._load_sample()
    win._run_analysis()
    assert _spin(lambda: win.state.get("analyzed_clean"))
    win.seq.setPlainText("ACGTACGTACGTACGTACGTACGTACGTACGT")   # change sequence after analysis
    assert win._stale_block() is True                          # design must be blocked


# ---------- WSL / optional-backend degradation (must never disable core science) ----------
def test_wsl_absent_keeps_optional_disabled(win):
    win._on_wsl_status({"wsl2": False})
    assert not win.annotateBtn.isEnabled() and not win.spliceBtn.isEnabled()
    assert "not installed" in win.wslStatus.text()
    assert not win.wslInstallBtn.isHidden()                    # installer stays reachable — it guides WSL setup


def test_wsl_broken_shows_error_not_crash(win):
    win._on_wsl_status({"error": "wsl.exe returned 0xffffffff"})
    assert "error" in win.wslStatus.text().lower()
    assert not win.annotateBtn.isEnabled()


def test_wsl_stack_missing_offers_install(win):
    win._on_wsl_status({"wsl2": True, "ready": False, "distro": "Ubuntu-24.04", "repeatmasker": None, "dfam": False})
    assert not win.wslInstallBtn.isHidden() and not win.annotateBtn.isEnabled()


def test_wsl_ready_enables_optional(win):
    win._on_wsl_status({"wsl2": True, "ready": True, "distro": "Ubuntu-24.04", "repeatmasker": "4.2.4", "minimap2": "2.28"})
    assert win.annotateBtn.isEnabled() and win.spliceBtn.isEnabled()
    assert not win.wslInstallBtn.isHidden()                    # installer stays reachable for repair / integrity checks


def test_core_science_independent_of_wsl():
    """The engine's core path must not touch WSL — analyze/primers/pcr work with WSL absent or broken."""
    import engine
    fx = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "tests", "fixtures", "M11240.fasta")
    r = engine.run_analyze({"sequence": open(fx).read()})
    assert r["records"][0]["classification"]["te_class"] == "LTR/Copia"


def test_annotate_without_stack_returns_soft_error():
    """run_annotate must degrade gracefully — never an unhandled crash. If the WSL annotation stack
    is absent it returns a soft error (ok:False / error); if the stack happens to be installed on
    this machine it succeeds cleanly. Either way it is a well-formed dict, never an exception."""
    import engine
    from engine import BadRequest
    from teagle_core import wsl
    try:
        res = engine.run_annotate({"sequence": ">x\n" + "ACGT" * 50})
    except BadRequest:
        return
    assert isinstance(res, dict)
    if wsl.env_status().get("ready"):
        assert res.get("ok") is True                            # stack present -> clean success
    else:
        assert res.get("ok") in (False, None) or "error" in res  # stack absent -> soft error, no crash
