"""Backend installer dialog (PySide6). A dedicated window that shows every component of the WSL
annotation stack — WSL2, micromamba, RepeatMasker, minimap2, the two Dfam libraries and the FamDB
config — each with a live status tick, a per-component Repair button, plus Install-all and a deep
Check-integrity pass and a live log. Runs all WSL work off the GUI thread through its own Engine
worker so the window stays responsive. Designed to make install failures diagnosable and fixable
one component at a time on any PC (no/broken WSL, no/broken conda env, partial Dfam download)."""
from __future__ import annotations
import os, re, sys

from PySide6.QtCore import Qt, QTimer
from PySide6.QtGui import QFont, QDesktopServices
from PySide6.QtCore import QUrl
from PySide6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QGridLayout, QLabel, QPushButton,
                               QPlainTextEdit, QScrollArea, QWidget, QFrame, QApplication)

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.join(os.path.dirname(_HERE), "backend"))
from engine_worker import Engine
from widgets import uppercase_buttons          # one casing rule for every window (see widgets.uppercase_buttons)

_ICON = {"ok": ("✓", "#33D6B8"), "bad": ("✕", "#E06A5A"), "work": ("●", "#E6A23C"), "unknown": ("—", "#8A959D")}


class InstallDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("TEagle — backend installer")
        self.setObjectName("central")
        import theme as _t
        z = _t.UI_SCALE                                       # open larger, capped to the screen, so every component row is comfortably readable
        try:
            _g = self.screen().availableGeometry(); _sw, _sh = _g.width(), _g.height()
        except Exception:
            _sw, _sh = 1440, 900
        self.resize(min(round(900 * z), _sw - 60), min(round(820 * z), _sh - 60))
        self._fitted = False             # the grow-to-content pass runs once, after the rows exist
        self.engine = Engine(self)
        self.engine.done.connect(self._on_done)
        self.engine.user_error.connect(self._on_user_error)
        self.engine.failed.connect(self._on_failed)
        self._rows = {}
        self._log_seen = ""
        self._busy = False
        self._wsl2_ok = True             # last-seen WSL2 availability (drives Install-all routing)
        self._wsl2_installing = False    # a Windows-side elevated WSL install is in progress
        self._wsl2_ticks = 0
        self._wsl2_reboot_pending = False  # a WSL install completed this session but Windows needs a restart

        root = QVBoxLayout(self)
        root.setContentsMargins(16, 14, 16, 14); root.setSpacing(10)
        title = QLabel("Backend installer"); title.setObjectName("sech")
        root.addWidget(title)
        intro = QLabel("This installs the optional Linux (WSL) annotation stack used for Dfam family "
                       "naming and de-novo splice detection. The domain-based superfamily classification "
                       "works without any of this. Each component installs and repairs independently — a "
                       "failure in one never blocks the others.")
        intro.setObjectName("orient"); intro.setWordWrap(True); root.addWidget(intro)

        # component grid inside a scroll area
        holder = QWidget(); self.grid = QGridLayout(holder)
        self.grid.setContentsMargins(2, 2, 2, 2); self.grid.setHorizontalSpacing(10); self.grid.setVerticalSpacing(4)
        self.grid.setColumnStretch(1, 1)
        self.scroll = scroll = QScrollArea(); scroll.setWidgetResizable(True); scroll.setWidget(holder)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)   # rows wrap their text; never scroll sideways
        # the list takes every spare pixel (stretch 1) and _fit_to_content grows the window to the whole
        # stack where the screen allows; the minimum stays low so a 720p screen can still fit the window
        # (a taller floor pushed the action bar and Close off the bottom of the display)
        scroll.setMinimumHeight(round(150 * _t.UI_SCALE)); root.addWidget(scroll, 1)

        # action bar
        bar = QHBoxLayout()
        self.installAllBtn = QPushButton("⭳ Install / update all"); self.installAllBtn.setProperty("primary", True)
        self.installAllBtn.clicked.connect(self._install_all); bar.addWidget(self.installAllBtn)
        self.integBtn = QPushButton("✔ Check integrity"); self.integBtn.setProperty("sm", True)
        self.integBtn.clicked.connect(self._check_integrity); bar.addWidget(self.integBtn)
        self.refreshBtn = QPushButton("↻ Refresh"); self.refreshBtn.setProperty("sm", True)
        self.refreshBtn.clicked.connect(self._refresh); bar.addWidget(self.refreshBtn)
        bar.addStretch(1)
        self.statusLine = QLabel(""); self.statusLine.setObjectName("cardmeta"); bar.addWidget(self.statusLine)
        root.addLayout(bar)

        loglbl = QLabel("INSTALL LOG"); loglbl.setObjectName("kdim"); root.addWidget(loglbl)
        self.log = QPlainTextEdit(); self.log.setReadOnly(True)
        self.log.setMinimumHeight(round(96 * z))     # floor only: it keeps its taller sizeHint wherever there is room
        lf = QFont("Cascadia Code"); lf.setStyleHint(QFont.Monospace); lf.setPointSize(9); self.log.setFont(lf)
        root.addWidget(self.log)

        close = QPushButton("Close"); close.setProperty("sm", True); close.clicked.connect(self.accept)
        crow = QHBoxLayout(); crow.addStretch(1); crow.addWidget(close); root.addLayout(crow)
        uppercase_buttons(self)

        self._poll = QTimer(self); self._poll.setInterval(2500); self._poll.timeout.connect(self._tick)
        # accept()/reject()/X all emit finished (closeEvent fires only for the X) — stop polling on every path
        self.finished.connect(lambda *_: self._poll.stop())
        QTimer.singleShot(0, self._refresh)

    # ---------- sizing ----------
    def _fit_to_content(self):
        """Grow the window once, after the component rows exist, so the whole stack is visible where the
        screen allows — the last row (the ~6.9 GB Dfam root download) used to sit clipped against the
        action bar. Where the stack still does not fit, the list scrolls to it."""
        if self._fitted or not self._rows:
            return
        self._fitted = True
        vp = self.scroll.viewport()
        # ask the grid what the rows need at the width they already have — the scroll area has not resized
        # its holder yet at this point, so the holder's own height still reads as the viewport height
        need = self.grid.heightForWidth(vp.width()) if self.grid.hasHeightForWidth() else self.grid.sizeHint().height()
        extra = need - vp.height()
        if extra <= 0:
            return
        try:
            g = self.screen().availableGeometry()
        except Exception:
            return
        import theme as _t                              # UI_SCALE is live, so read it at call time
        frame = max(self.frameGeometry().height() - self.height(), round(32 * _t.UI_SCALE))   # title bar
        h = min(self.height() + extra, g.height() - frame)
        if h > self.height():
            self.resize(self.width(), h)
        fg = self.frameGeometry()                       # keep the grown window on the screen it grew on
        if fg.bottom() > g.bottom():
            self.move(fg.x(), max(g.top(), fg.y() - (fg.bottom() - g.bottom())))

    # ---------- rows ----------
    def _ensure_row(self, key, name, desc):
        if key in self._rows:
            return self._rows[key]
        r = self.grid.rowCount()
        icon = QLabel("—"); icon.setFont(QFont("", 12)); icon.setFixedWidth(18); icon.setAlignment(Qt.AlignCenter)
        nm = QLabel(name); nf = nm.font(); nf.setBold(True); nm.setFont(nf)
        detail = QLabel("…"); detail.setObjectName("cardmeta")
        detail.setWordWrap(True)     # a long guide line used to widen the grid past the viewport and push every button out of reach
        btn = QPushButton("Repair"); btn.setProperty("sm", True)
        btn.clicked.connect(lambda _=False, k=key: self._repair(k))
        self.grid.addWidget(icon, r, 0)
        cell = QWidget(); cl = QVBoxLayout(cell); cl.setContentsMargins(0, 2, 0, 2); cl.setSpacing(0)
        top = QLabel(name); tf = top.font(); tf.setBold(True); top.setFont(tf)
        d = QLabel(desc); d.setObjectName("orient"); d.setWordWrap(True)
        cl.addWidget(top); cl.addWidget(detail); cl.addWidget(d)
        self.grid.addWidget(cell, r, 1)
        self.grid.addWidget(btn, r, 2, Qt.AlignTop)
        self._rows[key] = {"icon": icon, "detail": detail, "btn": btn, "repairable_now": False}
        return self._rows[key]

    def _set_icon(self, key, state):
        glyph, color = _ICON[state]
        ic = self._rows[key]["icon"]
        ic.setText(glyph); ic.setStyleSheet(f"color:{color}; font-weight:700;")

    # ---------- refresh / render ----------
    def _refresh(self):
        self.statusLine.setText("checking…")
        self.engine.submit("wsl_components", key="components")

    def _render_components(self, res):
        if res.get("error"):
            self.statusLine.setText(res["error"][:70]); return
        installing = res.get("installing")
        self._wsl2_ok = bool(res.get("wsl2"))
        if self._wsl2_ok:
            self._wsl2_reboot_pending = False       # WSL came up -> the reboot is done, clear the sticky state
        for c in res.get("components", []):
            row = self._ensure_row(c["key"], c["name"], c.get("desc", ""))
            ok = c.get("ok")
            # a not-yet-ok repairable component during an active install shows the working glyph
            state = "ok" if ok else ("work" if (installing and c.get("repairable")) else ("unknown" if not c.get("repairable") and not ok else "bad"))
            self._set_icon(c["key"], state)
            det = c.get("detail", "")
            if c.get("guide"):
                det += "  —  " + c["guide"]
            row["detail"].setText(det)
            if c["key"] == "wsl2":
                if self._wsl2_reboot_pending and not ok:
                    # a just-completed install shows as registered-but-won't-start pre-reboot; DON'T advise
                    # the destructive unregister — the distro is fine, Windows just needs a restart
                    row["detail"].setText("installed — restart Windows to finish setup, then reopen this installer")
                    row["repairable_now"] = False
                    row["btn"].setVisible(False)
                else:
                    # WSL itself installs from an elevated helper — the button is live even though WSL is absent
                    actionable = bool(c.get("installable")) and not ok
                    row["repairable_now"] = actionable
                    row["btn"].setEnabled(actionable and not self._busy)
                    row["btn"].setText("Install WSL")
                    row["btn"].setVisible(actionable)
            else:
                repairable = bool(c.get("repairable") and res.get("wsl2"))
                row["repairable_now"] = repairable      # remembered so _set_busy can restore it, not read a clobbered isEnabled()
                row["btn"].setEnabled(repairable and not self._busy)
                row["btn"].setText("Repair" if not ok else "Reinstall")
                row["btn"].setVisible(bool(c.get("repairable")))
        if self._wsl2_installing:
            pass                                        # keep the install-in-progress status line
        elif self._wsl2_reboot_pending:
            self.statusLine.setText("Restart Windows to finish WSL setup, then reopen this installer")
        elif not res.get("wsl2"):
            self.statusLine.setText("WSL2 not installed — click Install WSL on the first row")
        elif res.get("ready"):
            # ready == family naming (RepeatMasker + both Dfam partitions); only claim splice / genome scan
            # when their own components are actually installed, so the line never over-promises a missing tool
            ok = {c.get("key"): c.get("ok") for c in res.get("components", [])}
            feats = ["family naming"]
            if ok.get("minimap2"):
                feats.append("splice detection")
            if ok.get("genomescan"):
                feats.append("whole-genome scan")
            self.statusLine.setText("● ready — " + " & ".join(feats) + " available")
        elif installing:
            self.statusLine.setText("installing…")
        else:
            self.statusLine.setText(f"not ready · {res.get('disk_free_gb','?')} GB free")
        uppercase_buttons(self)                         # rows are created here and their labels re-set per state
        if not self._fitted:
            QTimer.singleShot(0, self._fit_to_content)  # after the new rows have been laid out

    # ---------- operations ----------
    def _set_busy(self, busy, note=""):
        self._busy = busy
        for b in (self.installAllBtn, self.integBtn, self.refreshBtn):
            b.setEnabled(not busy)
        for r in self._rows.values():
            r["btn"].setEnabled((not busy) and r.get("repairable_now", False))
        if note:
            self.statusLine.setText(note)

    def _append_log(self, text):
        if not text:
            return
        # show only the new suffix so the pane doesn't rebuild every poll
        if text.startswith(self._log_seen):
            new = text[len(self._log_seen):]
        else:
            self.log.setPlainText(text); self._log_seen = text
            self.log.verticalScrollBar().setValue(self.log.verticalScrollBar().maximum()); return
        if new:
            self.log.setPlainText(text)
            self.log.verticalScrollBar().setValue(self.log.verticalScrollBar().maximum())
        self._log_seen = text

    def _start(self, op, body, note):
        self._log_seen = ""; self.log.clear()
        self._set_busy(True, note)
        self.engine.submit(op, body, key=op)

    def _install_all(self):
        if not self._wsl2_ok:
            self._repair("wsl2")            # WSL must exist before the Linux stack; route to the elevated installer
            return
        self._start("wsl_install", None, "starting full install…")

    def _repair(self, key):
        if key == "wsl2":
            self._wsl2_installing = True; self._wsl2_ticks = 0
            self._start("wsl_install_wsl2", None, "installing WSL2 — accept the Windows (UAC) prompt…")
        else:
            self._start("wsl_repair", {"component": key}, f"repairing {key}…")

    def _check_integrity(self):
        self._set_busy(True, "running integrity check…")
        self.engine.submit("wsl_integrity", key="wsl_integrity")

    def _tick(self):
        if self._wsl2_installing:
            self._wsl2_ticks += 1
            if self._wsl2_ticks > 480:              # ~20 min guard so a stuck WSL install never polls forever
                self._poll.stop(); self._set_busy(False); self._wsl2_installing = False
                self.statusLine.setText("WSL install timed out — see the log, or run `wsl --install` in an admin PowerShell")
                return
            self.engine.submit("wsl2_install_log", key="wsl2_log")
        else:
            self.engine.submit("wsl_install_log", key="log")
        self.engine.submit("wsl_components", key="components")

    # ---------- results ----------
    def _on_done(self, key, res):
        if key == "components":
            self._render_components(res)
        elif key in ("wsl_install", "wsl_repair", "wsl_install_wsl2"):
            if not res.get("started"):
                self._set_busy(False); self._wsl2_installing = False
                self.statusLine.setText("could not start: " + str(res.get("error", "unknown")))
                return
            self._poll.start(); self._tick()
        elif key == "wsl2_log":
            log = res.get("log", "")
            self._append_log(log)
            m = re.search(r"DONE-WSL\s+(-?\d+)", log)
            if m or "[teagle] FAILED" in log:
                self._poll.stop(); self._set_busy(False); self._wsl2_installing = False
                if "[teagle] FAILED" in log:
                    self.statusLine.setText("WSL install failed — see the log")
                elif m and m.group(1) != "0":          # elevation ran but wsl.exe returned nonzero
                    self._wsl2_reboot_pending = True
                    self.statusLine.setText("WSL install reported errors — see the log; a Windows restart may still be required")
                else:
                    self._wsl2_reboot_pending = True    # success: pre-reboot the distro shows registered-but-unstartable
                    self.statusLine.setText("WSL installed — restart Windows to finish setup, then reopen this installer")
                self._refresh()
        elif key == "log":
            log = res.get("log", "")
            self._append_log(log)
            if "[teagle] DONE" in log or "[teagle] FAILED" in log:
                self._poll.stop(); self._set_busy(False)
                self.statusLine.setText("install finished — verify with Check integrity"
                                        if "[teagle] DONE" in log else "install reported a failure — see log")
                self._refresh()
        elif key == "wsl_integrity":
            self._set_busy(False)
            lines = ["=== integrity check ==="]
            for c in res.get("checks", []):
                lines.append(f"[{'OK' if c['ok'] else 'FAIL'}] {c['name']} — {c['detail']}")
            lines.append("RESULT: " + ("all checks passed" if res.get("ok") else "problems found — repair the failing components"))
            if res.get("error"):
                lines.append("error: " + str(res["error"]))
            self.log.setPlainText("\n".join(lines)); self._log_seen = self.log.toPlainText()
            self.statusLine.setText("integrity: " + ("OK" if res.get("ok") else "problems found"))

    def _on_user_error(self, key, msg):
        self._poll.stop(); self._set_busy(False); self._wsl2_installing = False
        self.statusLine.setText(msg[:70])

    def _on_failed(self, key, msg, trace):
        self._poll.stop(); self._set_busy(False); self._wsl2_installing = False
        self.statusLine.setText("error: " + msg[:60])
        sys.stderr.write(trace + "\n")

    def closeEvent(self, e):
        self._poll.stop()
        super().closeEvent(e)


def _selftest():
    """Offscreen construction check (bundled selftest gate): the dialog builds and lists rows."""
    app = QApplication.instance() or QApplication([])
    d = InstallDialog()
    d._render_components({"wsl2": True, "installing": False, "ready": False, "disk_free_gb": "100",
                          "components": [{"key": "micromamba", "name": "micromamba", "desc": "x", "ok": True, "detail": "installed", "repairable": True},
                                         {"key": "wsl2", "name": "WSL2", "desc": "y", "ok": True, "detail": "Ubuntu", "repairable": False}]})
    assert "micromamba" in d._rows and "wsl2" in d._rows
    return 0


if __name__ == "__main__":
    sys.exit(_selftest())
