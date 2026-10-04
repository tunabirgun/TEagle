"""TEagle native desktop app (PySide6). QMainWindow shell with a specimen rail and a scrollable
column of collapsible result cards. All science runs in-process through the shared engine, off the
GUI thread via engine_worker.Engine. This module wires the analyze workflow, primer/PCR, WSL family
annotation, splice detection, provenance and exports."""
from __future__ import annotations
import functools, gzip, html, json, math, os, re, sys

from PySide6.QtCore import Qt, QTimer, QByteArray, QSettings, QPropertyAnimation, QEasingCurve, QPoint
from PySide6.QtGui import (QGuiApplication, QFont, QPixmap, QPainter, QIcon, QCursor, QShortcut,
                           QKeySequence, QColor)
from PySide6.QtSvg import QSvgRenderer
from PySide6.QtWidgets import (QApplication, QMainWindow, QWidget, QFrame, QVBoxLayout, QHBoxLayout,
                               QGridLayout, QLabel, QLineEdit, QTextEdit, QPlainTextEdit, QPushButton, QComboBox,
                               QScrollArea, QSplitter, QFileDialog, QSizePolicy, QToolTip, QMessageBox, QDialog,
                               QTableWidget, QTableWidgetItem, QHeaderView, QAbstractItemView, QSpinBox,
                               QStyledItemDelegate)

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.join(os.path.dirname(_HERE), "backend"))

from engine_worker import Engine
import fonts
import figures
from figures import svg_genome
import widgets
from widgets import FigurePanel, GenomePanel, DataTable, BusyBar
import coords
from coords import CoordTable, CoordGenomePanel
from sample import make_sample
import theme as theme_mod
from teagle_core import appdirs, classify, sequtil  # classify.DOMAINS_TESTED = the tested-profile panel (scope caveat)
from teagle_core.wsl import (curated_coverage_sentence as _curated_coverage,
                             dfam_coverage_counts as _dfam_coverage_counts)  # measured, never retyped on a panel
from teagle_core import domains as domains_mod      # DOMAIN_INFO: the methods panel is derived from it, never retyped
from teagle_core import structural as structural_mod   # MIN_TSD/MAX_TSD: the TSD method text states the real window
from teagle_core import primers as primers_mod
from teagle_core import oligoqc                      # ΔG tier thresholds: tooltips read them, never retype them
from teagle_core.fetch import (COORD_ASSEMBLIES, all_assemblies,                        # pinned + user-added assemblies
                               complete_gene_model, cross_check_models, retrieve)        # + gene model / transcript fetch
from teagle_core import __version__ as APP_VERSION    # single source of truth (never hardcode a duplicate version)
from engine import LARGE_INPUT_BP                         # the engine owns the documented single-locus threshold
from ui_common import (
    CARD_ANIM_MS, CollapsibleCard, DOMAIN_COLS, EXPORT_HEADERS, GLOSS, ICON_TEAL, MARK_H, ORF_COLS, PANEL,
    PRIMER_COLS, REFLINKS, REPO_URL, STRUCT_COLS, UI_SCALES, WORD_H, WSL_ORGANISMS, _ACCESSION_LIKE,
    _CONFIDENCE_TIP, _COORD_SYSTEM, _COORD_TAIL, _Combo, _DFAM_YEAST_ADDITIONAL, _DFAM_YEAST_CURATED,
    _DFAM_YEAST_TOTAL, _ERV_TIP, _FLAG_MARK, _FLAG_ORDER, _MIN_ORF_AA, _MetaLabel, _NON_DNA_LETTERS,
    _PUBLISHED_CLASS, _QWIDGETSIZE_MAX, _RANK_TIP, _SIG, _STATUS_WORD, _TE_CLASS_TEXT, _WICKER, _amp_call,
    _amp_kind, _amplicon_public_source, _app_icon, _call_view, _cap_first, _caution_items, _clear_layout,
    _coding_view, _dg_clauses, _dg_limits, _dg_tiers, _disable_result_exports, _domain_segment, _domains_fasta,
    _empty, _example_accession, _example_digests, _example_match, _explanation_view, _export_fasta_btn,
    _export_table_btn, _external, _fasta_text, _full_records_fasta, _invalid_record_text, _kb_links, _link_focus,
    _load_asset, _manifest_export_stem, _mark_pixmap, _metric_cell, _no_orf_to_search, _note, _note_view,
    _orf_budget_tie, _orfs_fasta, _panel_references, _public_record_id, _published_note, _qq_note, _rank_line,
    _rank_view, _record_export_stem, _record_public_id, _record_public_identity, _repeat_word,
    _repeated_domain_set, _safe_export_stem, _selftest_positive_control, _sl, _splice_exons_fasta, _svg_pixmap,
    _te_headline, _word_pixmap, invalid_input_message, motion_enabled, open_external_url)
from panel_chrome import ChromeMixin
from panel_specimen import SpecimenMixin
from panel_records import RecordsMixin
from panel_structure import StructureMixin
from panel_splice import SpliceMixin
from panel_primers import PrimersMixin
from panel_pcr import PcrMixin
from panel_genome import GenomeScanMixin
from panel_genome_annotate import GenomeAnnotateMixin
from panel_family import FamilyMixin
from panel_provenance import ProvenanceMixin
from panel_report import ReportMixin


def _saved_coord_base():
    """The remembered display convention for coordinates: 1 (1-based, inclusive) or 0; anything else reads as 0."""
    try:
        return 1 if int(QSettings("TEagle", "TEagle").value("coord_base", 0)) == 1 else 0
    except (TypeError, ValueError):
        return 0


class MainWindow(ChromeMixin, SpecimenMixin, RecordsMixin, StructureMixin, SpliceMixin, PrimersMixin, PcrMixin, GenomeScanMixin, GenomeAnnotateMixin, FamilyMixin, ProvenanceMixin, ReportMixin, QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("TEagle")
        self.setWindowIcon(_app_icon())                   # also set on the top-level window (some Qt/Windows paths key the taskbar HICON off WM_SETICON)
        # open at the design size but never larger than the screen — a 1366x768 laptop must not get an
        # 860px-tall window taller than its display; keep a usable floor so the split layout stays coherent.
        # The floor is scaled by UI_SCALE to match the scaled rail/toolbar minimums (else at 1.5x the results
        # toolbar clips), then clamped to the screen below so a large scale never overhangs a small display.
        w, h = round(1240 * theme_mod.UI_SCALE), round(860 * theme_mod.UI_SCALE)
        minw, minh = round(820 * theme_mod.UI_SCALE), round(560 * theme_mod.UI_SCALE)
        scr = self.screen().availableGeometry() if self.screen() else None
        if scr is not None:
            aw, ah = scr.width() - 40, scr.height() - 60
            w, h = min(w, aw), min(h, ah)
            minw, minh = min(minw, aw), min(minh, ah)
        self.setMinimumSize(minw, minh)
        self.resize(max(w, minw), max(h, minh))
        _saved_theme = str(QSettings("TEagle", "TEagle").value("theme", "light"))   # first run -> light; then remembered
        self.theme = _saved_theme if _saved_theme in ("dark", "light") else "light"
        self.coord_base = _saved_coord_base()                 # on-screen coordinate convention only; files keep 0-based half-open
        self.state = {"seq": "", "source": None, "last_rec": None}
        self._loading = False                                 # True while a programmatic load writes the specimen box
        self._pcr_run = None
        self._pcr_gen = 0                                    # monotonic in-silico-PCR batch id (drops stale sibling results)
        self._spec_gen = 0                                    # bumped by every invalidation: a primer/PCR result sent under an older one is dropped
        self._specimen_gen = 0                                # bumped only when the specimen itself changes (family, splice, scan, fetch)
        self._recTimer = QTimer(self); self._recTimer.setSingleShot(True); self._recTimer.setInterval(150)
        self._recTimer.timeout.connect(self._apply_record_selection)   # debounces a run of arrow-key row changes
        self._recTable = self._recHead = self._recDetail = None        # the persistent record table, its heading, the rebuilt detail block
        self._rec_pending = None
        self._req_n = 0
        self._reqs = {}                                       # kind -> newest live tagged request {key, gen, digest}
        self._box_sig = ""                                    # whitespace-free text of the specimen box (a pure whitespace edit is no change)
        self._wsl_ready = {"annotate": False, "splice": False}
        self._analyze_gen = 0                                 # monotonic analysis id (drops stale/out-of-order callbacks)
        self._active_analyze_gen = None
        self._analysis_requests = {}
        self._analyze_inflight = False                        # gates Run analysis alongside "is there a sequence yet?"
        self._design_inflight = False                         # one primer design at a time (self._design_tmpl is shared state)
        self._genome_inflight = False                         # one whole-genome isPcr scan at a time
        self._genome_prep_inflight = False                    # one genome download/prepare at a time (large, one-time)
        self._add_asm_inflight = False                        # one custom-assembly resolve at a time (survives a dialog close/reopen)
        self._annot_inflight = False                          # one whole-genome TE annotation at a time (hours; huge disk/CPU)
        self._annot_result = None                             # last completed landscape (kept for the results window + exports)
        self._pending_scan = None                             # a scan queued behind a just-started genome download
        self._prepared_genomes = []                           # downloaded+verified (.done) genomes; the ONLY source for the PCR organism dropdown

        self.engine = Engine(self)
        self.engine.done.connect(self._on_done)
        self.engine.user_error.connect(self._on_user_error)
        self.engine.failed.connect(self._on_failed)

        central = QWidget(); central.setObjectName("central")
        self.setCentralWidget(central)
        outer = QVBoxLayout(central)
        outer.setContentsMargins(theme_mod.sp(12), theme_mod.sp(10), theme_mod.sp(12), theme_mod.sp(10))
        outer.setSpacing(theme_mod.sp(8))
        outer.addWidget(self._build_header())

        split = QSplitter(Qt.Horizontal)
        self.rail = self._build_rail()
        split.addWidget(self.rail)
        split.addWidget(self._build_results())
        split.setStretchFactor(0, 0)
        split.setStretchFactor(1, 1)
        split.setSizes([340, 900])
        self.split = split
        self._rail_sizes = [340, 900]                          # remembered split for reopening the collapsed rail
        self._rail_collapsed = False                           # explicit state (isVisible() is unreliable pre-show)
        outer.addWidget(split, 1)

        self._apply_theme()
        self.acc.setFocus(Qt.OtherFocusReason)                 # the first thing to do is name a sequence, not hide a panel
        QTimer.singleShot(0, self, self._startup)             # context object: a window freed first never runs it

    # ---------- header ----------
    def resizeEvent(self, e):
        # hide the decorative tagline on a narrow window so it never contributes to an overhang; the Ignored
        # size policy already lets the window shrink, this just keeps the header clean at small widths.
        super().resizeEvent(e)
        tag = getattr(self, "_tagline", None)
        if tag is not None:
            tag.setVisible(self.width() >= round(1080 * theme_mod.UI_SCALE))
        if hasattr(self, "coordBtn"):
            self._paint_coord_btn()

    def _set_ui_scale(self, factor):
        """Apply a new UI scale live: re-scale fonts + spacing immediately by re-applying the _scale_px'd QSS.
        No restart, no QT_SCALE_FACTOR (retired). theme_mod.UI_SCALE is seeded before MainWindow so a persisted
        scale is pixel-exact for code-set geometry from the next launch."""
        if abs(factor - theme_mod.UI_SCALE) < 1e-3:
            return
        QSettings("TEagle", "TEagle").setValue("ui_scale", float(factor))
        theme_mod.UI_SCALE = float(factor)
        self._apply_theme()                                   # re-applies scaled QSS + redraws wordmark/mark at the new scale
        z = theme_mod.UI_SCALE                                # re-scale the code-set rail width live so it does not crowd
        self.rail.setMinimumWidth(round(320 * z)); self.rail.setMaximumWidth(round(440 * z))
        inner = self.rail.widget() if hasattr(self.rail, "widget") else None
        if inner is not None:
            inner.setMinimumWidth(round(300 * z)); inner.setMaximumWidth(round(430 * z))
        if getattr(self, "_genome_mgr", None) is not None:    # an open manager was built at the old scale — rebuild it live
            self.engine.submit("genome_list", {}, key="genome_list")
        self._banner(f"UI scale set to {int(factor * 100)}% — applied live (font + spacing); pixel-exact on the next launch.", "info")

    def _set_coord_base(self, base, persist=True):
        """Switch the on-screen numbering and redraw what is showing, without rebinding the record or touching a running
        design (so it re-renders the cards, never _show_record)."""
        base = 1 if base == 1 else 0
        if base == self.coord_base:
            return
        self.coord_base = base
        if persist:
            QSettings("TEagle", "TEagle").setValue("coord_base", base)
        self._paint_coord_btn()
        st = self.state
        if st.get("last_rec") and st.get("analysis_state") == "current":
            self._render_struct_card(st["last_rec"], st.get("analysis") or {}, keep_table=True)
        if st.get("lastPrimers") and not self._design_inflight:
            self._render_primers(st["lastPrimers"])
        if st.get("lastPcr"):
            lp = st["lastPcr"]
            self._render_pcr(lp["lanes"], lp["amplicons"], lp["statuses"])
        if st.get("splice"):
            self._render_splice(st["splice"])
        if st.get("familyResult"):
            self._render_family(st["familyResult"])
        self._refresh_assay_goal()
        self._uppercase_buttons()
        self._banner("Coordinates are now shown 1-based, inclusive. Exported tables, FASTA names and primer-design coordinates "
                     "are unchanged (0-based, half-open)." if base else
                     "Coordinates are now shown 0-based, half-open — the convention TEagle's exported tables use.", "info")

    # ---------- quick start, help and glossary ----------
    def _quickstart_hidden(self):
        return str(QSettings("TEagle", "TEagle").value("quickstart_hidden", False)).lower() in ("true", "1")

    def _set_quickstart_hidden(self, flag):
        QSettings("TEagle", "TEagle").setValue("quickstart_hidden", bool(flag))
        self._quickstart.setVisible(not flag)

    def _open_manual(self):
        return open_external_url(f"{REPO_URL}/releases/tag/v{APP_VERSION}")

    def _toggle_theme(self):
        self.theme = "light" if self.theme == "dark" else "dark"
        QSettings("TEagle", "TEagle").setValue("theme", self.theme)   # remember the choice for next launch
        self._apply_theme()
        d = self.state.get("lastPrimers")                     # QC ΔG cell colours are per-theme -> re-render on toggle
        if isinstance(d, dict) and d.get("candidates"):
            self._render_primers(d)

    @staticmethod
    def _revcomp(s):
        return "".join(MainWindow._RC.get(c, "N") for c in reversed(s.upper()))

    def _use_record_transcript(self):
        """Fetch the record's own selected mRNA into the transcript box and flag origin='record' so _render_splice
        reports a consistency check, never 'independent confirmation' (the check is circular by construction)."""
        acc = self.spliceTxSel.currentData() if hasattr(self, "spliceTxSel") else None
        if not acc:
            return
        # A cache hit can still trigger a feature re-fetch with retry back-off, so the fetch runs off the GUI thread
        # and the window stays responsive; the result is applied here on the GUI thread when it lands.
        import threading
        self.spliceTxUse.setEnabled(False)
        out, specimen = {}, self._specimen_gen

        def work():
            try:
                out["fasta"] = retrieve(acc).get("fasta", "")
            except Exception as e:
                out["error"] = e.args[0] if getattr(e, "args", None) else e
        threading.Thread(target=work, daemon=True).start()
        poll = QTimer(self); poll.setInterval(80)

        def done():
            if not out:
                return
            poll.stop(); poll.deleteLater()
            self.spliceTxUse.setEnabled(True)
            if specimen != self._specimen_gen:                # the specimen changed while this was in flight
                return
            if "error" in out:
                return self._banner(f"Could not fetch {acc}: {out['error']}", "warn")
            self._tx_programmatic = True
            self.spliceTx.setPlainText(out["fasta"])
            self._tx_programmatic = False
            self._splice_tx_origin = "record"                 # same annotation source -> consistency check, not independent
            self.card_splice.expand(); self.spliceTx.setFocus()
        poll.timeout.connect(done); poll.start()

    @staticmethod
    def _structure_text(c):
        qc = c.get("qc") or {}
        if not qc.get("ok"):
            return f"{MainWindow._pair_label(c)}: secondary-structure QC unavailable"
        L, R = qc["left"], qc["right"]
        def line(name, m):
            return f"{name}\tprimer3 {m.get('p3')}\tViennaRNA {m.get('vrna')}\t{m.get('flag')}\t{m.get('agree')}"
        return "\n".join([f"# {MainWindow._pair_label(c)} secondary structure (ΔG kcal/mol)",
                          f"# F {c['left_seq']} / R {c['right_seq']}",
                          line("hairpin_F", L["hairpin"]), line("hairpin_R", R["hairpin"]),
                          line("selfdimer_F", L["self_dimer"]), line("selfdimer_R", R["self_dimer"]),
                          line("heterodimer", qc["hetero_dimer"]), line("end_stability", qc["end_stability"])])

    def _design_locus_for(self, cand, org):
        """Where the DESIGNED product sits in the genome, as the one-element list the engine expects:
        {accession, start, end, strand, intended_product_span, intended_endpoint_tolerance_bp}, the fetched
        window (1-based inclusive) plus the designed product's span inside that record. Sent only for a pair
        whose own sealed design record is a coordinate fetch of the SAME assembly being scanned and is still the
        record the specimen on screen holds; anything else returns None and the scan stays neutral ('genomic
        priming sites'). The region is taken from the pair's sealed design input, so in a multi-region fetch it is
        the region the pair was designed on, never another one."""
        inp = (cand.get("design_provenance") or {}).get("input") or {}
        src = self.state.get("source") or {}
        scanned_asm = (all_assemblies().get(org, {}) or {}).get("assemblyAccession", "")
        regions = inp.get("regions") or []
        if (not scanned_asm or inp.get("inputType") != "fetched_coordinate" or len(regions) != 1 or
                inp.get("assemblyAccession") != scanned_asm or src.get("assemblyAccession") != scanned_asm or
                cand.get("intended_product_span") is None):
            return None
        r = regions[0]
        try:
            i = (src.get("regions") or []).index(r)
            on_screen = (src.get("records") or [])[i].get("sha256") == inp.get("sha256")
        except (ValueError, IndexError, AttributeError):
            return None                                       # the specimen on screen is not the one it was designed on
        if not on_screen or not r.get("chrAccession"):
            return None
        return [{"accession": r["chrAccession"], "start": r.get("start"), "end": r.get("end", r.get("stop")),
                 "strand": r.get("strand"), "intended_product_span": cand["intended_product_span"],
                 "intended_endpoint_tolerance_bp": cand.get("intended_endpoint_tolerance_bp", 0)}]

    def _prepare_genome(self, org, then_scan=None):
        """Download + cache an organism's RefSeq genome for local scanning (one-time, non-blocking)."""
        if self._genome_prep_inflight or self._genome_inflight:   # also block while a scan runs (reachable via the manager)
            return self._banner("A whole-genome scan / download is already running — wait for it to finish.", "info")
        if not org:
            return self._banner("Select an organism to download.", "info")
        self._pending_scan = then_scan
        self._genome_prep_inflight = True
        self._refresh_genome_manager()                        # an open manager must lock its row actions while this runs
        self.card_genome.expand()
        acc = (all_assemblies().get(org, {}) or {}).get("assemblyAccession", "")
        self._prep_org = org
        self._render_genome_busy(f"Downloading the {org} genome ({acc}) — one-time, kept for future scans. "
                                 "This can take several minutes (larger for mammalian genomes)…")
        self.engine.submit("genome_prepare", {"organism": org}, key="genome_prepare")
        if getattr(self, "_prep_timer", None) is None:        # poll the prepare log so a long download shows liveness
            self._prep_timer = QTimer(self); self._prep_timer.setInterval(2500)
            self._prep_timer.timeout.connect(self._poll_prepare_log)
        self._prep_timer.start()

    def _rebuild_coord_asm_dropdown(self):
        """(Re)populate the coordinate-fetch organism dropdown from all_assemblies(), so a newly added custom
        organism appears in-session (preserving the current selection). The '__custom__' free-text item stays last."""
        box = getattr(self, "asmSel", None)
        if box is None:
            return
        prev = box.currentData()
        box.blockSignals(True); box.clear()
        _asm = all_assemblies()
        for org in sorted(_asm):
            box.addItem(f"{org} · {_asm[org]['assemblyName']}", org)
        box.addItem("Other organism / assembly…", "__custom__")
        if prev is not None:
            i = box.findData(prev)
            if i >= 0:
                box.setCurrentIndex(i)
        box.blockSignals(False)

    def _refresh_genome_dropdown(self):
        """Rebuild the PCR organism dropdown from the prepared (downloaded+verified) genome set only —
        an organism cannot be scanned until it is downloaded. userData stays the ORGANISM string (what
        _scan_genome and the need_prepare fallback consume via currentData()); a downloaded accession
        outside the curated map is skipped. Selection is preserved across the rebuild."""
        box = getattr(self, "genomeOrg", None)
        if box is None:                                       # genome_list can resolve before the PCR card is built
            return
        if not getattr(self, "_fields_decorated", False) and "pcrPmax" in getattr(self, "pfields", {}):
            # cards 04 and 05 exist by the time card 06 builds, and this runs once while it does: the primer and
            # in-silico PCR fields get their units and tooltips without editing the card builders
            self._fields_decorated = True
            self._decorate_primer_fields()
        prev = box.currentData()
        asm = all_assemblies()                                # curated + user-added, so a downloaded custom genome shows here too
        acc2org = {v["assemblyAccession"]: k for k, v in asm.items()}
        box.blockSignals(True)
        box.clear()
        box.addItem("— select organism —", None)              # placeholder so no wrong-species scan runs by default
        n = 0
        for g in self._prepared_genomes:
            org = acc2org.get(g.get("accession"))
            if not org:
                continue
            box.addItem(f"{org} · {asm[org]['assemblyName']}", org)
            n += 1
        if prev is not None:
            i = box.findData(prev)
            if i >= 0:
                box.setCurrentIndex(i)
        box.blockSignals(False)
        hint = getattr(self, "genomeOrgHint", None)
        if hint is not None:
            hint.setText("" if n else "No genomes downloaded yet — open Manage genomes to download one for scanning.")
        self._refresh_scan_picker()                           # genome availability gates the scan button

    def _on_genome_list(self, d):
        # This handler is the SINGLE fan-out for every genome_list submit (WSL-ready, download, delete,
        # manager-open). ALWAYS update the dropdown first; the dialog rebuild below is the only part guarded
        # by manager-open — so a download finishing with the manager closed still refreshes the dropdown.
        if d.get("ok"):                                       # only overwrite on success — a transient WSL blip must not
            self._prepared_genomes = d.get("genomes", [])     # blank the dropdown while genomes are still cached on disk
        self._refresh_genome_dropdown()
        self._rebuild_coord_asm_dropdown()                    # a newly added custom organism must also appear in coord-fetch

        dlg = getattr(self, "_genome_mgr", None)
        if dlg is None:
            if not getattr(self, "_genome_mgr_open", False):  # a refresh that landed AFTER the user closed the manager
                return                                        # -> only skip the DIALOG; the dropdown is already updated
            dlg = QDialog(self); dlg.setAttribute(Qt.WA_DeleteOnClose)
            dlg.finished.connect(lambda *_: (setattr(self, "_genome_mgr", None),
                                             setattr(self, "_genome_mgr_open", False)))
            self._genome_mgr = dlg
        from PySide6.QtGui import QColor
        dlg.setWindowTitle("Manage genomes")
        z = theme_mod.UI_SCALE
        try:                                                  # the dialog is sized to fit the screen + fixed (non-resizable) at the end
            scrg = self.screen().availableGeometry(); scrw, scrh = scrg.width(), scrg.height()
        except Exception:
            scrw, scrh = 1440, 900
        old = dlg.layout()
        if old is not None:                                   # rebuild the body in place on refresh
            dlg.hide()                                        # hide while swapping the layout, else the re-added widgets measure as hidden and sizeHint()
            QWidget().setLayout(old)                          # collapses to ~20px -> setFixedSize would pin the whole dialog to a sliver (the final dlg.show() re-shows it)
        lay = QVBoxLayout(dlg)
        lay.setContentsMargins(theme_mod.sp(12), theme_mod.sp(10), theme_mod.sp(12), theme_mod.sp(10)); lay.setSpacing(theme_mod.sp(7))
        title = QLabel("Manage genomes"); tf = title.font(); tf.setBold(True); tf.setPointSizeF(tf.pointSizeF() + 2)
        title.setFont(tf); lay.addWidget(title)
        if not d.get("ok"):                                   # WSL-down error card (styled, not a bare label)
            err = QLabel("Could not list genomes — " + d.get("error", "WSL backend unavailable"))
            err.setObjectName("errbanner"); err.setWordWrap(True); lay.addWidget(err)
            close = QPushButton("Close"); close.clicked.connect(dlg.accept)
            close.setAccessibleName("Close the genome manager"); lay.addWidget(close)
            # this branch returns BEFORE the table fit-tail that owns setFixedWidth/adjustSize/setFixedSize, and the
            # dialog is REUSED — so release the size a previous successful build pinned on it, else the error card is
            # stranded inside the old table's geometry with a screenful of dead space.
            self._uppercase_buttons(dlg)                  # uppercase BEFORE the fit tail measures the dialog
            ew = min(round(560 * z), scrw - 56)
            dlg.setMinimumSize(0, 0); dlg.setMaximumSize(16777215, 16777215)   # QWIDGETSIZE_MAX
            dlg.setFixedWidth(ew); dlg.adjustSize()
            dlg.setFixedSize(ew, min(dlg.sizeHint().height(), scrh - 56))
            dlg.show(); dlg.raise_(); return
        intro = QLabel("Download a genome once to enable its whole-genome off-target scan; it is kept locally and then "
                       "appears in the organism menus. Every genome download needs at least 8 GB free disk in WSL, "
                       "however small the genome. Mammalian genomes are large (~1 GB): extraction peaks near 4 GB "
                       "(FASTA ~3 GB + 2bit ~0.8 GB) before the temporary files are removed.")
        intro.setObjectName("orient"); intro.setWordWrap(True); lay.addWidget(intro)
        # lock every action while any download/scan/annotation runs. The annotation belongs in this guard:
        # it reads the cached genome for hours, and Delete does `rm -rf` on that same directory.
        busy = self._genome_prep_inflight or self._genome_inflight or self._annot_inflight
        # ADD-ORGANISM row — one input, auto-detects an organism name vs a GCF/GCA accession; resolved + pinned once.
        arow = QHBoxLayout()
        self._addAsmEdit = QLineEdit()
        self._addAsmEdit.setPlaceholderText("Add an organism by name or assembly accession — e.g. Danio rerio or GCF_000002035.6")
        addBtn = QPushButton("Add"); addBtn.setProperty("sm", True); addBtn.setEnabled(not busy and not self._add_asm_inflight)
        addBtn.setAccessibleName("Add the typed organism or assembly accession to the genome list")
        addBtn.clicked.connect(self._add_custom_assembly); self._addAsmEdit.returnPressed.connect(self._add_custom_assembly)
        self._addAsmBtn = addBtn                               # kept so the add handlers can disable/re-enable it across a dialog rebuild
        if self._add_asm_inflight:                            # a resolve started in a since-closed dialog is still running
            self._addAsmEdit.setEnabled(False)
        arow.addWidget(self._addAsmEdit, 1); arow.addWidget(addBtn); lay.addLayout(arow)
        prepared = {g["accession"]: g for g in d.get("genomes", [])}
        asm = all_assemblies(); orgs = sorted(asm)            # curated + user-added
        pal = theme_mod._DARK if self.theme == "dark" else theme_mod._LIGHT
        rows = len(orgs)
        # size every row so the WHOLE table (all rows + columns + buttons) shows with NO scrolling in a fixed dialog.
        # rows stay comfortable on a normal screen and shrink (with the font) only if a small screen forces it. The final
        # row height is computed in the tail after the real chrome (title+intro+add-row+close) is measured at the true width.
        avail_h = scrh - round(56); avail_w = scrw - round(56)
        # ONE fixed, comfortable row height, chosen so the natural-height 'sm' button + a small margin always fit (the
        # button is NEVER forced short or shrunk -> never squeezed/clipped). Set BEFORE the cells are built so each cell
        # widget is laid out at the right height. A screen too small to show every row scrolls vertically (button stays
        # comfortable) rather than squeezing the button — a readable button matters more than avoiding a scroll.
        rowH = round(34 * z)
        target_w = min(round(1180 * z), avail_w)              # bigger, so full organism + assembly names show without eliding
        headers = ["Organism", "Assembly", "Accession", "Status", "Size (MB)", "Contigs", "Action"]
        tbl = QTableWidget(len(orgs), len(headers))
        tbl.setHorizontalHeaderLabels(headers)
        tbl.verticalHeader().setVisible(False)
        tbl.setEditTriggers(QAbstractItemView.NoEditTriggers)
        tbl.setSelectionMode(QAbstractItemView.NoSelection)
        tbl.setSortingEnabled(False)                          # static catalog; sorting conflicts with setCellWidget buttons
        tbl.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        # Organism AND Assembly (the two widest text columns) stretch, so on a narrow dialog THEY shrink+elide while the
        # short fixed columns — including Action — keep their content width. The Download button can never be the column
        # that overflows off the right edge into a horizontal scrollbar (the small-screen clip a single stretch column left).
        tbl.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        tbl.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
        tbl.setTextElideMode(Qt.ElideRight)                   # a shrunk Organism/Assembly elides ("Drosophila melan…"), never clips Action
        tbl.setWordWrap(False)                                # (h-scrollbar stays as-needed — a pathological tiny screen keeps Action reachable by scroll rather than silently clipped)
        tbl.horizontalHeader().setStretchLastSection(False)
        tbl.horizontalHeader().setMinimumSectionSize(round(56 * z))
        tbl.horizontalHeader().setDefaultAlignment(Qt.AlignCenter)   # header + cells both centred (they disagreed before)
        tbl.verticalHeader().setDefaultSectionSize(rowH)      # provisional row height; the real value is set in the tail
        aligns = [Qt.AlignHCenter] * 6                        # every column centred, matching the centred header
        row_btns = []                                         # collected so every row's button gets ONE uniform size
        for r, org in enumerate(orgs):
            meta = asm[org]; acc = meta["assemblyAccession"]; g = prepared.get(acc); downloaded = bool(g)
            cells = [org, meta.get("assemblyName", ""), acc, "Downloaded" if downloaded else "Not downloaded",
                     (f"{(g.get('bytes', 0) or 0) / 1e6:.0f}" if downloaded else "—"),
                     (str(g.get("n_seqs", "?")) if downloaded else "—")]
            for c, txt in enumerate(cells):
                it = QTableWidgetItem(txt); it.setTextAlignment(aligns[c] | Qt.AlignVCenter)
                if c in (0, 1):                               # Organism/Assembly may elide on a narrow dialog — tooltip recovers the full text
                    it.setToolTip(txt)
                if c == 3:                                    # status carried by colour + word (green = downloaded)
                    it.setForeground(QColor(pal["good"] if downloaded else pal["faint"]))
                tbl.setItem(r, c, it)
            btn = QPushButton("Delete" if downloaded else "Download")
            btn.setProperty("cellbtn", True)                  # the ONLY style property on it (see theme.py) -> no specificity clash
            btn.ensurePolished()                              # apply the stylesheet BEFORE sizeHint is read, else it measures the base font
            btn.clicked.connect((lambda _=False, a=acc, o=org: self._remove_genome(a, o)) if downloaded
                                else (lambda _=False, o=org: self._on_manager_download(o)))
            btn.setEnabled(not busy)
            holder = QWidget()                                # full-cell holder; AlignCenter centres the button both ways in whatever
            hl = QHBoxLayout(holder)                          # geometry the view hands the cell -> no top/bottom bias, no clipping
            hl.setContentsMargins(0, 0, 0, 0); hl.setSpacing(0)
            hl.addWidget(btn, 0, Qt.AlignCenter)
            tbl.setCellWidget(r, 6, holder)
            row_btns.append(btn)
        # ONE uniform button size for every row, from the widest label ("Download" > "Delete"), so the buttons align
        # with each other instead of each sizing to its own text. Fixed size => the holder can never squeeze or stretch it.
        # Qt applies QTableWidget::item padding to CELL WIDGETS as well as text: it offsets the widget by the
        # padding's top-left AND shrinks its usable box by twice the padding. So the row/column must be sized to
        # (button + 2x padding), or the button is both pushed off-centre and sheared. Read the padding from the
        # theme so this can never drift from the stylesheet.
        pv, ph = round(theme_mod.TABLE_ITEM_PAD_V * z), round(theme_mod.TABLE_ITEM_PAD_H * z)
        if row_btns:
            for b in row_btns:                                # _uppercase_buttons() runs LATER (any theme/scale/result
                b.setText(widgets.caps(b.text()))                   # refresh reaches this dialog) — measure the FINAL text or
            bw = max(b.sizeHint().width() for b in row_btns)  # the pinned width is the sentence-case one and "DOWNLOAD"
            bh = max(b.sizeHint().height() for b in row_btns) # gets clipped to "OWNLOA"
            for b in row_btns:
                b.setFixedSize(bw, bh)                        # uniform size => the buttons align with each other
            rowH = max(rowH, bh + 2 * pv + round(4 * z))      # inset height (rowH - 2*pv) always exceeds the button
            tbl.verticalHeader().setDefaultSectionSize(rowH)
        tbl.resizeColumnsToContents()
        if row_btns:
            # resizeColumnsToContents() measures ITEMS only — it is blind to setCellWidget widgets, which is why the
            # Action column used to size to the word "Action" and let the button spill past the table's right edge.
            tbl.horizontalHeader().setSectionResizeMode(6, QHeaderView.Fixed)
            tbl.setColumnWidth(6, bw + 2 * ph + round(6 * z))
        for r in range(rows):                                 # the view owns the holder's geometry (cell minus the
            tbl.setRowHeight(r, rowH)                          # symmetric item padding) -> AlignCenter is true centring
        lay.addWidget(tbl)
        self._genome_mgr_table = tbl
        if busy:
            what = ("genome download" if self._genome_prep_inflight else "whole-genome scan" if self._genome_inflight
                    else "genome annotation")
            b = QLabel(f"Actions are disabled until the {what} finishes.")
            b.setObjectName("orient"); b.setWordWrap(True); lay.addWidget(b)
        close = QPushButton("Close"); close.setProperty("sm", True); close.clicked.connect(dlg.accept)
        close.setAccessibleName("Close the genome manager")
        # Whole-genome TE annotation lives in the footer rather than a per-row button: the Action column's
        # width and row height are pinned from ONE measured button, and a second per-row control would
        # reopen the squeeze/clipping problems that sizing was written to prevent.
        annot = QPushButton("Annotate TE landscape…"); annot.setProperty("sm", True)
        annot.setAccessibleName("Annotate every transposable element in a downloaded genome")
        annot.setToolTip("Run RepeatMasker over a downloaded genome and summarise its transposable-element content")
        annot.setEnabled(bool(prepared) and not busy and not self._annot_inflight)
        annot.clicked.connect(self._annotate_genome_dialog)
        crow = QHBoxLayout(); crow.addStretch(1); crow.addWidget(annot); crow.addWidget(close)
        crow.addStretch(1); lay.addLayout(crow)
        self._uppercase_buttons(dlg)                          # ADD/CLOSE match the row DOWNLOAD/DELETE; measured below
        # ---- fit: full table if it fits the screen, else cap to whole rows + a vertical scroll (button never squeezed) ----
        dlg.setFixedWidth(target_w)                           # set the real width first so the intro wraps as it truly will
        hdr_h = tbl.horizontalHeader().sizeHint().height() or round(34 * z)
        tbl.setFixedHeight(hdr_h + rows * rowH + round(4 * z))
        tbl.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        dlg.adjustSize()
        if dlg.sizeHint().height() > avail_h:                 # too tall for this screen -> show as many WHOLE rows as fit + scroll
            over = dlg.sizeHint().height() - avail_h
            vis = max(1, (rows * rowH - over) // rowH)
            tbl.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)
            tbl.setFixedHeight(hdr_h + vis * rowH + round(4 * z))
            dlg.adjustSize()
        dlg.setFixedSize(target_w, min(dlg.sizeHint().height(), avail_h))
        dlg.show(); dlg.raise_()

    def _annotate_genome_dialog(self):
        """Choose a downloaded genome and the run's parameters, and state the cost BEFORE anything starts.

        The cost panel is the scientific gate as much as the practical one: it reports how many family
        models the installed Dfam partitions actually hold for the chosen lineage. A lineage with none
        cannot yield a TE result at any runtime, and the user is told that here rather than after hours."""
        if self._annot_inflight:
            self._banner("A genome annotation is already running.", "info"); return
        # _prepared_genomes is the app's single store of what is cached on disk (kept current by the one
        # genome_list fan-out), so the picker can never disagree with the manager table beside it.
        prepared = [g for g in (getattr(self, "_prepared_genomes", None) or []) if g.get("accession")]
        if not prepared:
            self._banner("Download a genome first — annotation runs on a genome already on this machine.", "info")
            return
        asm = all_assemblies()
        acc_to_org = {m["assemblyAccession"]: o for o, m in asm.items()}
        dlg = QDialog(self); dlg.setWindowTitle("Annotate TE landscape"); dlg.setModal(False)
        dlg.setAttribute(Qt.WA_DeleteOnClose, True)
        lay = QVBoxLayout(dlg)
        lay.setContentsMargins(theme_mod.sp(14), theme_mod.sp(12), theme_mod.sp(14), theme_mod.sp(12))
        lay.setSpacing(theme_mod.sp(8))
        t = QLabel("Annotate the TE landscape of a downloaded genome")
        tf = t.font(); tf.setBold(True); tf.setPointSizeF(tf.pointSizeF() + 1); t.setFont(tf); lay.addWidget(t)
        intro = QLabel("TEagle runs RepeatMasker over the whole assembly and summarises which repeat names/models "
                       "match, how much of the genome their intervals cover, and their alignment-span-weighted "
                       "divergence. It matches references represented in the selected library; it does not discover "
                       "new repeats, so coverage depends on that library.")
        intro.setObjectName("orient"); intro.setWordWrap(True); lay.addWidget(intro)

        grid = QHBoxLayout()
        gsel = QComboBox()
        for g in prepared:
            org = acc_to_org.get(g["accession"], g["accession"])
            gsel.addItem(f"{org} · {g['accession']} · {(g.get('bytes') or 0)/1e6:.0f} MB", g["accession"])
        gsel.setToolTip("A genome already downloaded to this machine. Only downloaded genomes can be "
                        "annotated — the whole assembly is read locally, nothing is sent to a server.")
        _gl = QLabel("Genome"); _gl.setToolTip(gsel.toolTip())
        grid.addWidget(_gl); grid.addWidget(gsel, 1)
        ssel = QComboBox()
        for key, label in self._ANNOT_SENS:
            ssel.addItem(label, key)
        ssel.setToolTip("How hard RepeatMasker looks. The figures are RepeatMasker's own: Quick (-q) is "
                        "2–5× faster and 5–10% less sensitive; Slow (-s) is 2–3× slower and 0–5% more "
                        "sensitive. Sensitivity here means recovering older, more diverged copies — the "
                        "recent ones are found either way.")
        _sl2 = QLabel("Sensitivity"); _sl2.setToolTip(ssel.toolTip())
        grid.addWidget(_sl2); grid.addWidget(ssel, 1)
        lay.addLayout(grid)
        # LIBRARY CHOICE — the single setting that decides which families can be found at all, so it is
        # offered explicitly rather than assumed. The two options are alternatives, not additions:
        # RepeatMasker's -lib replaces the database search, so a run uses one or the other.
        lrow = QHBoxLayout()
        lsel = QComboBox()
        lsel.addItem("Installed Dfam library — curated families only", "")
        lsel.addItem("Installed Dfam library — include uncurated families", "uncurated")
        lsel.addItem("My own repeat library (FASTA)…", "custom")
        lsel.setToolTip(
            "Dfam families are either curated (reviewed) or uncurated (automatically built). RepeatMasker "
            "searches curated families only unless you ask for both.\n\n"
            "Outside a few heavily studied species most families are uncurated, so curated-only can find "
            f"nothing at all: on baker's yeast it sees {_DFAM_YEAST_CURATED} of {_DFAM_YEAST_TOTAL} family "
            f"models and reports no transposable element, while including uncurated adds "
            f"{_DFAM_YEAST_ADDITIONAL} models and finds the Ty elements the genome really has.\n\n"
            "Including them needs the optional uncurated partitions from the backend installer, and the "
            "extra families are lower confidence — the result says which setting produced it.\n\n"
            "A custom FASTA library replaces Dfam entirely for that run.")
        lrow.addWidget(QLabel("Library")); lrow.addWidget(lsel, 1)
        libpath = QLabel(""); libpath.setObjectName("orient"); libpath.setWordWrap(True)
        lrow.addWidget(libpath, 2); lay.addLayout(lrow)
        self._annot_custom_lib = None
        self._annot_uncurated = False

        def pick_library(idx):
            self._annot_uncurated = (lsel.currentData() == "uncurated")
            if lsel.currentData() != "custom":
                self._annot_custom_lib = None
                libpath.setText("Curated + uncurated families. Needs the uncurated Dfam partitions installed; "
                                "uncurated families are automatically built and lower confidence."
                                if self._annot_uncurated else "")
                refresh_cost()
                return
            fn, _ = QFileDialog.getOpenFileName(dlg, "Choose a repeat library (FASTA)", "",
                                                "FASTA library (*.fa *.fasta *.lib *.txt);;All files (*)")
            if not fn:
                lsel.setCurrentIndex(0); self._annot_custom_lib = None; libpath.setText(""); return
            self._annot_custom_lib = fn
            libpath.setText(f"{os.path.basename(fn)} — used as the sole repeat library; TEagle records its "
                            "checksum but cannot vouch for its contents.")
            refresh_cost()
        lsel.currentIndexChanged.connect(pick_library)
        cost = QLabel("Estimating…"); cost.setObjectName("orient"); cost.setWordWrap(True)
        cost.setTextInteractionFlags(Qt.TextSelectableByMouse); lay.addWidget(cost)
        warn = QLabel(""); warn.setObjectName("errbanner"); warn.setWordWrap(True); warn.hide(); lay.addWidget(warn)
        start = QPushButton("Start annotation"); start.setEnabled(False)
        retry = QPushButton("Check again"); retry.setProperty("sm", True); retry.hide()
        retry.setToolTip("Read this machine's compute and disk budget again")
        cancel = QPushButton("Close"); cancel.setProperty("sm", True); cancel.clicked.connect(dlg.reject)
        brow = QHBoxLayout(); brow.addStretch(1); brow.addWidget(retry); brow.addWidget(start); brow.addWidget(cancel)
        lay.addLayout(brow)
        self._annot_retry = retry
        self._annot_cost_state = {}

        def refresh_cost():
            acc = gsel.currentData()
            g = next((x for x in prepared if x["accession"] == acc), {})
            org = acc_to_org.get(acc, "")
            custom_path = self._annot_custom_lib
            library_mode = ("custom" if custom_path else
                            ("dfam_uncurated" if self._annot_uncurated else "dfam_curated"))
            cost.setText("Estimating cost — reading the WSL compute and disk budget…" if custom_path else
                         "Estimating cost — reading the WSL budget and the installed Dfam library…")
            start.setEnabled(False); warn.hide(); retry.hide()
            self._annot_budget_generation = int(getattr(self, "_annot_budget_generation", 0)) + 1
            self._annot_cost_state = {"accession": acc, "organism": org, "sha256": g.get("sha256"),
                                      "custom_library": custom_path, "library_mode": library_mode,
                                      "budget_request_id": self._annot_budget_generation}
            body = {"species": org, "genome_bytes": int(g.get("bytes") or 0), "library_mode": library_mode,
                    "budget_request_id": self._annot_budget_generation}
            if g.get("fasta_bytes"):                     # the source FASTA is what the disk estimate is sized from
                body["fasta_bytes"] = int(g["fasta_bytes"])
            self.engine.submit("annotate_budget", body, key="annotate_budget")

        self._annot_cost_widgets = (cost, warn, start, gsel, ssel)
        retry.clicked.connect(lambda *_: refresh_cost())
        gsel.currentIndexChanged.connect(lambda *_: refresh_cost())
        start.clicked.connect(lambda: (self._start_genome_annotate(gsel.currentData(), acc_to_org.get(gsel.currentData(), ""),
                                                                   ssel.currentData(), self._annot_custom_lib,
                                                                   getattr(self, "_annot_uncurated", False)), dlg.accept()))
        self._uppercase_buttons(dlg)
        refresh_cost()
        dlg.adjustSize(); dlg.show(); dlg.raise_()
        self._annot_dialog = dlg

    @staticmethod
    def _annot_report_md(r):
        """The landscape as a self-contained Markdown report.

        Every hedge the window shows travels with it: what was searched, what the two percentages mean,
        and why an empty TE result may be a library limit. An export that carried only the table would
        read as an unqualified statement about the genome."""
        repeat_names = r.get("repeat_names") or []
        te = [row for row in repeat_names if row.get("kind") == "TE"]
        library = MainWindow._annot_library_context(r)
        L = [f"# Transposable-element landscape — {r.get('species')} ({r.get('accession')})", "",
             f"- **Transposable elements: {r.get('te_percent')}%** of "
             f"{(r.get('genome_bp') or 0)/1e6:.1f} Mb, across "
             f"{r.get('te_repeat_name_count')} distinct repeat names/models",
             f"- All repeat classes together: {r.get('masked_percent')}% "
             f"({r.get('total_hits'):,} alignment rows)",
             library["report"],
             f"- Sensitivity {r.get('sensitivity')} · {r.get('chunks')} chunks · "
             f"complete: {bool(r.get('complete'))}", ""]
        if not r.get("complete"):
            L += ["> **Incomplete run.** Counts, percentages, and rows are partial lower bounds. No assembly-wide "
                  "absence, coverage, or exhaustive-family claim is supported.", ""]
        if r.get("coverage_warning"):
            L += ["> **Coverage warning.** " + r["coverage_warning"], ""]
        L += ["## Scope and limits", "",
              library["report_scope"],
              "- Transposable elements are counted separately from tandem repeats (simple repeats,",
              "  low-complexity sequence, satellites) and non-TE entries; only the first is a statement",
              "  about transposable-element content.",
              "- Each repeat name/model's Bases value is its within-name interval union. The all-repeat and TE headline "
              "percentages use separate cross-model unions, so a genomic base contributes at most once to each "
              "headline. Summing rows can exceed a headline because repeat models can overlap.",
              "- TE, tandem and other category percentages are independently unioned and non-additive.",
              "- Divergence is the alignment-span-weighted mean RepeatMasker divergence from consensus. Raw alignment",
              "  rows, including overlaps and fragments, supply the weights. It is not Kimura- or CpG-corrected,",
              "  a mutation rate, or an age.",
              "- Percentages use all submitted FASTA bases, including ambiguous bases and unplaced scaffolds.", "",
              "## Repeat names/models", "",
              "| Repeat name/model | Class/family | Kind | Copies (alignment rows) | Bases (within-name union) | % genome | Mean divergence % |",
              "|---|---|---|---|---|---|---|"]
        for row in repeat_names:
            L.append(f"| {row.get('repeat_name')} | {row.get('class_family')} | {MainWindow._annot_kind_label(row)} | {row.get('n'):,} | {row.get('bp'):,} | "
                     f"{row.get('percent')} | {row.get('divergence')} |")
        L += ["", f"_{len(te)} of {len(repeat_names)} rows have a recognized transposable-element class._"]
        return "\n".join(L)

    def _on_genome_pcr(self, d):
        if self._scan_is_stale():
            return self._discard_stale_scan()
        self._genome_inflight = False
        self._refresh_genome_manager()                        # scan settled — re-enable a manager opened mid-scan
        if not d.get("ok"):
            if d.get("cache_corrupt"):
                ps = self._pending_scan or {}
                org = ps.get("org") or self.genomeOrg.currentData()
                acc = (all_assemblies().get(org, {}) or {}).get("assemblyAccession", "")
                self._pending_scan = None
                msg = (f"The cached {org or 'selected'} genome ({acc or 'unknown accession'}) failed its "
                       "integrity check. Open Manage genomes, remove that cached genome, then download it "
                       "again before scanning. TEagle did not run isPcr on the altered cache.")
                self._render_genome_status(msg, "warn")
                return self._banner(msg, "warn")
            if d.get("need_prepare"):                         # genome not downloaded yet — offer the one-time download
                ps = self._pending_scan or {}
                org = ps.get("org") or self.genomeOrg.currentData()   # the org that STARTED this scan, not the live dropdown
                acc = (all_assemblies().get(org, {}) or {}).get("assemblyAccession", "")
                self._render_genome_status(f"The {org} genome ({acc}) is not downloaded yet.", "info")
                box = QMessageBox(self)
                box.setWindowTitle("Download genome?")
                box.setText(f"The {org} genome ({acc}) is not on this machine yet.\n\n"
                            "Download it once now? It is kept locally for future scans. Every genome download needs at "
                            "least 8 GB free disk in WSL, however small the genome. Mammalian genomes are large (~1 GB "
                            "download, a few minutes) and extraction peaks near 4 GB (FASTA ~3 GB + 2bit ~0.8 GB) "
                            "before the temporary files are removed.")
                box.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
                box.setDefaultButton(QMessageBox.Yes)
                if org and box.exec() == QMessageBox.Yes:
                    self._prepare_genome(org, then_scan=ps)
                else:
                    self._pending_scan = None
                    self._render_genome_status("Genome scan cancelled — no genome downloaded.", "warn")
                return
            self._render_genome_status("No genome scan result — " + d.get("error", "scan failed"))
            return self._banner(d.get("error", "The genome scan failed."))
        # WHICH pair this verdict describes, captured before _pending_scan is dropped. A scan result says
        # things like "no off-target site in this assembly"; without the pair on the panel that sentence
        # reads as a property of the specimen rather than of two specific oligos.
        _ps = self._pending_scan or {}
        _cand = _ps.get("cand") or {}
        self._pending_scan = None
        amps = [{**a, "pair": "genome", "remote": True} for a in d.get("amplicons", [])]
        summary = d.get("summary") or {}
        has_locus = bool(summary.get("has_locus"))        # gel: on/off colours with a locus; neutral 'priming site' colour without
        lanes = [{"label": "genome", "amplicons": d.get("amplicons", []), "has_locus": has_locus, "advisory": not has_locus}]
        ctx = getattr(self, "_scan_ctx", None) or {}
        prm = (d.get("provenance") or {}).get("parameters") or {}
        self.state["lastGenomeScan"] = {"lanes": lanes, "amplicons": amps, "summary": summary,
                                        "pair": {"id": self._pair_label(_cand) if _cand else "", "fwd": _cand.get("left_seq", ""),
                                                 "rev": _cand.get("right_seq", "")},
                                        "organism": _ps.get("org", ""),
                                        "assembly": {"accession": d.get("assemblyAccession", ""), "name": d.get("assemblyName", "")},
                                        "criteria": {"isPcr": prm.get("isPcr_version"), "min_perfect": prm.get("min_perfect"),
                                                     "min_good": prm.get("min_good"),
                                                     "max_size": prm.get("max_size") or ctx.get("prod_max") or self._SCAN_PROD_MAX},
                                        "loci": ctx.get("loci") or []}
        self.card_genome.expand()
        self._render_genome_scan(lanes, amps, summary)
        m = d.get("provenance", {})
        db = (m.get("databases") or [{}])[0]
        seal = m.get("manifestSha256", "")
        prov = (f" · assembly {d.get('assemblyAccession','')}"
                + (f" · sha256 {db.get('sha256','')[:12]}…" if db.get("sha256") else "")
                + (f" · seal {seal[:12]}…" if seal else ""))
        self._banner(f"Whole-genome scan · {d.get('organism','')} ({d.get('assemblyName','')}): "
                     f"{summary.get('verdict', 'scan complete')} — under isPcr's ≥15 bp minimum perfect "
                     f"3′-end match requirement (sites outside that rule are not reported){prov}.", level="info")
        if m:
            self._render_provenance(m)

    def closeEvent(self, e):
        jobs = self._running_jobs()
        if jobs and not self._confirm_close(jobs):
            e.ignore()
            return
        # Stop EVERY polling timer, not a named one. This used to stop only the genome-download poll, so
        # each timer added later (the annotation progress poll, for one) could still fire while the window
        # was being torn down. Walking the children keeps that correct without anyone remembering to.
        for t in self.findChildren(QTimer):
            try:
                t.stop()
            except RuntimeError:                              # already destroyed by Qt — nothing to stop
                pass
        # Off-thread engine jobs are deliberately NOT waited on: a whole-genome annotation runs for hours
        # inside WSL and resumes on the next launch, so blocking the close would trade a slow exit for no
        # benefit. _Job._emit already drops results whose receiver has gone.
        super().closeEvent(e)


def selftest():
    """Bundle self-test (TEAGLE_SELFTEST=1). Proves the packaged build imports the compiled
    scientific stack (pyhmmer, primer3), ships the HMM profiles and TLS roots, renders figures through
    QtSvg, and runs an end-to-end analysis — the checks a double-click launch cannot report. Exit 0/1."""
    app = QApplication.instance() or QApplication([])
    import engine
    from teagle_core import domains, primers
    import figures
    from PySide6.QtSvg import QSvgRenderer
    from PySide6.QtGui import QImage, QPainter
    from PySide6.QtCore import QByteArray, Qt
    problems = []
    missing_fonts = fonts.load_fonts()                    # proves the bundled TTFs actually shipped (offscreen still loads explicit fonts)
    if missing_fonts:
        problems.append("bundled fonts missing from build: " + ", ".join(missing_fonts))
    if primers.PRIMER3_VERSION == "unavailable":
        problems.append(f"primer3 failed to load ({primers.PRIMER3_ERROR})")
    from teagle_core import oligoqc                          # secondary-structure QC; primer3 required, ViennaRNA optional
    if not oligoqc.available().get("primer3"):               # the PRIMARY engine must be in the bundle
        problems.append(f"primer3 secondary-structure engine missing from bundle "
                        f"({oligoqc.available().get('primer3_error')})")
    # ViennaRNA is deliberately NOT bundled — its licence forbids redistribution inside an AGPL work, so it
    # is a user-installed optional second engine. Its absence is the expected shipped state, never a build
    # failure; the UI already reports which engines ran and the manifest records it.
    if domains.PYHMMER_VERSION == "unavailable":
        problems.append("pyhmmer failed to load")
    if not domains.HMM_SHA256:
        problems.append("bundled Pfam HMM profiles not found")
    else:                                                 # the profile set must LOAD and cover the full GAG-POL-ENV panel
        try:
            n_hmm = len(domains._hmms())
            codes = {v[0] for v in domains.DOMAIN_INFO.values()}
            # every class the 30-model panel claims must be represented, incl. the LINE ORF1p/EN, the DIRS
            # tyrosine recombinase (YR) and the Helitron (HEL) modules added this release
            miss = {"GAG", "PR", "RT", "RNaseH", "INT", "ENV", "CHR", "TPase", "ORF1", "EN", "YR", "HEL"} - codes
            if miss:
                problems.append(f"domain profile set missing expected codes: {sorted(miss)}")
            if n_hmm < len(domains.DOMAIN_INFO):          # the bundled .hmm must load the full declared panel, not a truncated subset
                problems.append(f"domain profile set loaded {n_hmm} profiles, expected {len(domains.DOMAIN_INFO)} "
                                "(bundled .hmm may be truncated)")
        except Exception as e:
            problems.append(f"domain profile set failed to load: {type(e).__name__}: {e}")
    # the frozen launcher bypasses the startup environment check, so a dist-info folder or a module missing
    # from the bundle would only surface as a crash mid-analysis; fail the release gate here instead
    if getattr(sys, "frozen", False):
        try:
            import envcheck
            env = envcheck.check()
            if not env.get("packages_ok"):
                bad = [p["name"] for p in env.get("packages", []) if not p.get("ok")] + list(env.get("importFailures") or [])
                problems.append("bundled package check failed (missing metadata or import): " + ", ".join(sorted(set(bad))))
        except Exception as e:
            problems.append(f"bundled package check crashed: {type(e).__name__}: {e}")
    # end-to-end science through the shared engine (also exercises the fixture-free sample)
    try:
        r = engine.run_analyze({"sequence": make_sample()})
        if not r.get("records"):
            problems.append("analyze produced no records")
    except Exception as e:
        problems.append(f"analyze crashed: {type(e).__name__}: {e}")
    # QtSvg must render (the figure layer's single point of failure in a frozen build)
    svg = figures.svg_gel({"lanes": [{"label": "P1", "amplicons": [{"length": 200, "on_target": True, "source": "x"}]}]}, "dark")
    rd = QSvgRenderer(QByteArray(svg.encode("utf-8")))
    img = QImage(200, 150, QImage.Format_ARGB32); img.fill(Qt.transparent)
    p = QPainter(img); rd.render(p); p.end()
    if not (rd.isValid() and any(img.pixelColor(x, y).alpha() > 0 for x in range(0, 200, 10) for y in range(0, 150, 10))):
        problems.append("QtSvg did not render (figure plugin missing from bundle)")
    # XLSX table export must be importable AND functional (save() pulls openpyxl's lazy writer submodules)
    try:
        import io
        if not widgets._HAS_XLSX:
            problems.append("openpyxl (XLSX table export) missing from bundle")
        else:
            # go through the lazy accessor, which is the path a real export takes — this check exists to
            # prove openpyxl actually LOADS and writes inside the frozen bundle, not merely that it is present
            wb = widgets._xl()["Workbook"](); wb.active.append(["h", 1]); wb.save(io.BytesIO())
    except Exception as e:
        problems.append(f"XLSX export self-check failed: {type(e).__name__}: {e}")
    # Network retrieval is intentionally not exercised here, but its trust store must be functional. The
    # frozen launcher bypasses envcheck, and package metadata alone cannot prove that certifi's bundled CA
    # file is present or that OpenSSL can parse it into the context every HTTPS fetch uses.
    try:
        import certifi
        from teagle_core import fetch
        fetch._SSL_CTX_CACHE.clear()
        ca_path = certifi.where()
        tls_context = fetch._ssl_ctx()
        if not os.path.isfile(ca_path):
            problems.append("certifi CA bundle missing from build")
        elif tls_context is None or not tls_context.get_ca_certs():
            problems.append("certifi CA bundle did not load into the TLS context")
    except Exception as e:
        problems.append(f"TLS trust-store self-check failed: {type(e).__name__}: {e}")
    # Vector PDF export relies on the shared SVG -> QPdfWriter path and metadata minimisation. Exercise
    # the actual runtime writer so a frozen build cannot pass merely because the dependencies import.
    try:
        import tempfile
        from pypdf import PdfReader
        from teagle_core.artifact_metadata import media_metadata_violations
        with tempfile.TemporaryDirectory(prefix="teagle-selftest-") as temporary:
            path = os.path.join(temporary, "figure-export.pdf")
            widgets.render_pdf(svg, path)
            violations = media_metadata_violations(path, author=None)   # end-user exports carry no author
            reader = PdfReader(path)
            if violations:
                problems.append("PDF figure export metadata check failed: " + "; ".join(violations))
            if len(reader.pages) != 1:
                problems.append("PDF figure export did not preserve its page")
    except Exception as e:
        problems.append(f"PDF figure export self-check crashed: {type(e).__name__}: {e}")
    # the installer dialog must construct offscreen (it ships in the frozen build)
    try:
        from install_dialog import InstallDialog
        dlg = InstallDialog()
        dlg._render_components({"wsl2": True, "installing": False, "ready": False, "disk_free_gb": "50",
                                "components": [{"key": "micromamba", "name": "micromamba", "desc": "x",
                                                "ok": False, "detail": "missing", "repairable": True}]})
        if "micromamba" not in dlg._rows:
            problems.append("install dialog did not build component rows")
    except Exception as e:
        problems.append(f"install dialog crashed: {type(e).__name__}: {e}")
    if problems:
        sys.stderr.write("TEAGLE SELFTEST FAILED:\n  - " + "\n  - ".join(problems) + "\n")
        return 1
    _vr = oligoqc.VIENNARNA_VERSION if oligoqc.available().get("viennarna") else "not installed (optional)"
    print(f"TEAGLE SELFTEST OK · primer3 {primers.PRIMER3_VERSION} · ViennaRNA {_vr} "
          f"· pyhmmer {domains.PYHMMER_VERSION} · HMM {domains.HMM_SHA256[:12]} ({len(domains._hmms())} profiles) "
          f"· TLS roots ok · QtSvg ok · install dialog ok")
    return 0




def _apply_saved_ui_scale():
    """Seed the live UI scale from QSettings into theme_mod.UI_SCALE BEFORE MainWindow is built, so first
    paint is correct. Scale is driven at runtime (theme_mod.UI_SCALE + _apply_theme) — no QT_SCALE_FACTOR,
    no restart. An explicit TEAGLE_UI_SCALE env override (screenshot harness / dev) always wins."""
    if os.environ.get("TEAGLE_UI_SCALE"):
        return
    try:
        f = float(QSettings("TEagle", "TEagle").value("ui_scale", 1.0))
        if 0.5 <= f <= 3.0:
            theme_mod.UI_SCALE = f
    except Exception:
        pass


def main():
    if os.environ.get("TEAGLE_SELFTEST"):
        return selftest()
    _apply_saved_ui_scale()                               # must precede QApplication creation
    if sys.platform == "win32":                           # taskbar groups under our icon, not pythonw's, in dev runs
        try:
            import ctypes
            ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID("TEagle.desktop.2")
        except Exception:
            pass
    app = QApplication.instance() or QApplication(sys.argv)
    app.setApplicationName("TEagle")                      # title stays exactly "TEagle" (no auto-appended display name)
    fonts.load_fonts()                                    # bundled Cascadia Mono (UI) — no dependence on installed fonts
    app.setWindowIcon(_app_icon())
    w = MainWindow()
    w.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
