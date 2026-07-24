"""Qt stylesheet for the TEagle 'assay terminal' look — a monospace-forward precision instrument,
matching the web UI (app/web/app.css): exact near-black + teal-mint palette, sharp 2px edges,
tabular-mono data, uppercase tracked micro-labels, and a teal accent reserved for interactive
affordances. Okabe-Ito data hues live in the figures. Qt QSS has no letter-spacing / text-transform,
so the wordmark spacing and uppercased labels are applied in code (main.py)."""

import os, re

# Whole UI is Cascadia Mono (bundled with the app; see native/fonts.py). SANS is kept as a
# separate name for the few prose labels, but points at the same mono family per the design brief.
MONO = '"Cascadia Mono", "Cascadia Code", "Consolas", "Courier New", monospace'
SANS = MONO

# global UI zoom: scale fonts + padding uniformly (1px/2px borders & radii left untouched).
UI_SCALE = float(os.environ.get("TEAGLE_UI_SCALE", "1.0"))

_COMMON = """
* {{ font-family: {sans}; font-size: 13px; color: {text}; }}
QMainWindow, QWidget#central, QDialog, QMessageBox {{ background: {bg}; }}
QScrollArea {{ border: none; background: {bg}; }}
QScrollArea > QWidget > QWidget {{ background: transparent; }}
QLabel {{ background: transparent; }}

/* header wordmark + chrome (mono, tracked) */
QLabel#word {{ font-family: {mono}; font-weight: 700; font-size: 16px; }}
QLabel#tagline {{ font-family: {mono}; color: {faint}; font-size: 10px; }}
QLabel#ver {{ font-family: {mono}; color: {faint}; font-size: 9px; font-weight: 600;
    border: 1px solid {line}; border-radius: 3px; padding: 1px 5px; }}
QLabel#statusTxt {{ font-family: {mono}; color: {dim}; font-size: 11px; }}
QFrame#statuschip {{ border: 1px solid {line}; background: {panel2}; border-radius: 2px; }}
QLabel#led {{ border-radius: 4px; background: {faint}; min-width: 8px; max-width: 8px; min-height: 8px; max-height: 8px; }}
QLabel#led[live="true"] {{ background: {good}; border: 2px solid {goodsoft}; }}
QFrame#headrule {{ border: none; }}

/* section headers + number badges */
QLabel#secn {{ font-family: {mono}; color: {accent}; font-size: 11px; font-weight: 700;
    border: 1px solid {line2}; border-radius: 2px; padding: 1px 5px; min-width: 15px; }}
QLabel#sech {{ font-family: {sans}; font-weight: 650; font-size: 13px; }}

/* rail + cards (sharp, panel-2 headers) */
QFrame#rail {{ background: {panel}; border: none; border-right: 1px solid {line}; }}
QFrame#card {{ background: {panel}; border: 1px solid {line}; border-radius: 2px; }}
QFrame#hline {{ background: {line}; max-height: 1px; border: none; }}
QPushButton#cardhdr {{ font-family: {mono}; background: {panel2}; border: none;
    border-bottom: 1px solid {line}; text-align: left; padding: 11px 14px; font-size: 12px; font-weight: 700; color: {text}; }}
QPushButton#cardhdr:hover {{ color: {accent}; }}
QPushButton#cardhdr:checked {{ border-bottom: none; }}

/* data-entry fields = mono, near-black well, sharp */
QLineEdit, QTextEdit, QPlainTextEdit {{ background: {bg}; border: 1px solid {line2}; border-radius: 2px;
    padding: 7px 9px; color: {text}; font-family: {mono}; font-size: 12px;
    selection-background-color: {accentsoft}; selection-color: {text}; }}
QComboBox {{ background: {bg}; border: 1px solid {line2}; border-radius: 2px; padding: 6px 9px;
    color: {text}; font-family: {mono}; font-size: 12px; }}
QLineEdit:focus, QTextEdit:focus, QComboBox:focus {{ border: 1px solid {accent}; }}
QComboBox::drop-down {{ border: none; width: 18px; }}
QComboBox QAbstractItemView {{ background: {panel3}; color: {text}; border: 1px solid {line2};
    selection-background-color: {accentsoft}; selection-color: {text}; font-family: {mono}; outline: none; }}

/* buttons = mono, uppercase (applied in code), sharp 2px */
QPushButton {{ background: {panel2}; border: 1px solid {line2}; border-radius: 2px; padding: 8px 13px;
    color: {text}; font-family: {mono}; font-size: 12px; font-weight: 600; }}
QPushButton:hover {{ border: 1px solid {accent}; background: {panel}; }}
QPushButton:pressed {{ background: {panel3}; }}
QPushButton:disabled {{ color: {faint}; border: 1px solid {line}; background: {panel}; }}
QPushButton[primary="true"] {{ background: {accent}; color: {accentink}; border: 1px solid {accent}; font-weight: 700; }}
QPushButton[primary="true"]:hover {{ background: {accent2}; }}
QPushButton[primary="true"]:disabled {{ background: {panel2}; color: {faint}; border: 1px solid {line}; }}
QPushButton[sm="true"] {{ padding: 5px 9px; font-size: 11px; }}
QPushButton[link="true"] {{ background: transparent; border: none; color: {accent};
    font-family: {mono}; font-size: 11px; text-align: left; padding: 2px; font-weight: 600; }}
QPushButton[link="true"]:hover {{ color: {accent2}; }}

/* metric readout gauges + key/value chrome */
QFrame#cell {{ background: {panel}; border: 1px solid {line}; border-radius: 2px; }}
QLabel#kdim {{ font-family: {mono}; color: {faint}; font-size: 9px; font-weight: 600; }}
QLabel#value {{ font-family: {mono}; font-size: 21px; font-weight: 600; color: {text}; }}
QLabel#value[state="good"] {{ color: {good}; }}
QLabel#value[state="bad"] {{ color: {bad}; }}

/* tables = mono data, uppercase tracked headers */
QTableWidget {{ background: {panel}; gridline-color: {line}; border: 1px solid {line}; border-radius: 2px;
    font-family: {mono}; font-size: 12px; alternate-background-color: {panel2}; outline: none; }}
QTableWidget::item {{ padding: 4px 7px; }}
QTableWidget::item:selected {{ background: {accentsoft}; color: {text}; }}
QHeaderView::section {{ background: {panel2}; color: {faint}; padding: 6px 8px; border: none;
    border-bottom: 1px solid {line}; font-family: {mono}; font-size: 10px; font-weight: 700; }}
QTableCornerButton::section {{ background: {panel2}; border: none; }}

/* classification banner (left accent border, big title) */
QFrame#classbn {{ background: {panel2}; border: 1px solid {line2}; border-left: 3px solid {accent}; border-radius: 2px; }}
QLabel#classbig {{ font-family: {sans}; font-size: 18px; font-weight: 700; color: {text}; }}
QLabel#classkls {{ font-family: {mono}; font-size: 11px; color: {dim}; }}
QLabel#classexp {{ font-family: {sans}; font-size: 12px; color: {dim}; }}
QLabel#cf {{ font-family: {mono}; font-size: 10px; font-weight: 700; padding: 3px 9px; border-radius: 2px; }}
QLabel#cf[level="High"] {{ background: {goodsoft}; color: {good}; }}
QLabel#cf[level="Moderate"] {{ background: {panel3}; color: {dim}; }}
QLabel#cf[level="Candidate"] {{ background: {warnsoft}; color: {warn}; }}

/* prose + misc chrome */
QLabel#orient {{ font-family: {sans}; color: {dim}; font-size: 11px; background: {panel2};
    border-left: 2px solid {line2}; border-radius: 2px; padding: 8px 11px; }}
QLabel#sectionlabel {{ font-family: {mono}; color: {faint}; font-size: 10px; font-weight: 700; }}
QLabel#cardmeta {{ font-family: {mono}; color: {dim}; font-size: 10.5px; }}
QLabel#empty {{ font-family: {mono}; color: {faint}; font-size: 12px; }}
QLabel#gvpos {{ font-family: {mono}; color: {faint}; font-size: 10px; }}
QLabel#errbanner {{ font-family: {mono}; background: {badsoft}; color: {bad}; border: 1px solid {bad};
    border-radius: 2px; padding: 9px 12px; font-size: 11.5px; }}
QLabel#errbanner[level="success"] {{ background: {goodsoft}; color: {good}; border: 1px solid {good}; }}
QLabel#errbanner[level="warn"] {{ background: {warnsoft}; color: {warn}; border: 1px solid {warn}; }}
QLabel#errbanner[level="info"] {{ background: {accentsoft}; color: {accent}; border: 1px solid {accent}; }}
QLabel#queuerow {{ font-family: {mono}; color: {dim}; font-size: 11px; }}
QProgressBar {{ background: {panel2}; border: 1px solid {line2}; border-radius: 2px; }}
QProgressBar::chunk {{ background: {accent}; border-radius: 2px; }}

/* scrollbars — thin, terminal */
QScrollBar:vertical {{ background: {bg}; width: 10px; margin: 0; }}
QScrollBar::handle:vertical {{ background: {line2}; border-radius: 2px; min-height: 30px; }}
QScrollBar::handle:vertical:hover {{ background: {accent}; }}
QScrollBar:horizontal {{ background: {bg}; height: 10px; margin: 0; }}
QScrollBar::handle:horizontal {{ background: {line2}; border-radius: 2px; min-width: 30px; }}
QScrollBar::add-line, QScrollBar::sub-line {{ height: 0; width: 0; }}
QScrollBar::add-page, QScrollBar::sub-page {{ background: transparent; }}
QSplitter::handle {{ background: {line}; }}
QMenu {{ background: {panel3}; border: 1px solid {line2}; font-family: {mono}; font-size: 11.5px; }}
QMenu::item {{ padding: 6px 14px; }}
QMenu::item:selected {{ background: {accentsoft}; color: {text}; }}
QToolTip {{ background: {panel3}; color: {text}; border: 1px solid {line2}; font-family: {mono}; font-size: 11px; padding: 4px 6px; }}
"""

_DARK = dict(bg="#0A0D10", panel="#0F1317", panel2="#131A1F", panel3="#19222A", line="#1E272E",
             line2="#2B3740", text="#E6EDF1", dim="#8A959D", faint="#7E8A93",
             accent="#33D6B8", accent2="#1FB89C", accentink="#042420", accentsoft="rgba(51,214,184,0.14)",
             good="#40C088", goodsoft="rgba(64,192,136,0.16)", warn="#D8B368", warnsoft="rgba(216,179,104,0.16)",
             bad="#EC5F49", badsoft="rgba(236,95,73,0.14)", mono=MONO, sans=SANS)
_LIGHT = dict(bg="#EDF1F3", panel="#FFFFFF", panel2="#F2F5F7", panel3="#E7ECEF", line="#DCE3E7",
              line2="#C6D0D6", text="#141B21", dim="#57636B", faint="#5E6A72",
              accent="#0E9E86", accent2="#0B8571", accentink="#FFFFFF", accentsoft="rgba(14,158,134,0.12)",
              good="#178A5C", goodsoft="rgba(23,138,92,0.14)", warn="#8A6D22", warnsoft="rgba(138,109,34,0.16)",
              bad="#C6432E", badsoft="rgba(198,67,46,0.12)", mono=MONO, sans=SANS)

# header accent underline (teal fading to transparent), like the web header::after
HEADRULE = {"dark": "qlineargradient(x1:0, y1:0, x2:1, y2:0, stop:0 #33D6B8, stop:0.42 transparent)",
            "light": "qlineargradient(x1:0, y1:0, x2:1, y2:0, stop:0 #0E9E86, stop:0.42 transparent)"}
ACCENT = {"dark": "#33D6B8", "light": "#0E9E86"}
TEXT = {"dark": "#E6EDF1", "light": "#141B21"}         # wordmark "TE" ink, per theme
# per-theme flag colours for QC ΔG cells — reuse the WCAG-tuned bad/warn palette (dark vs light), so amber/red
# read correctly on both backgrounds instead of a single hardcoded dark-tuned hex
FLAG = {"dark":  {"warn": _DARK["bad"],  "caution": _DARK["warn"]},
        "light": {"warn": _LIGHT["bad"], "caution": _LIGHT["warn"]}}


def _scale_px(css: str, f: float) -> str:
    if abs(f - 1.0) < 1e-3:
        return css
    # scale font/padding px (>=6px) only; keep 1-2px borders and 2px radii crisp
    return re.sub(r"(\d+(?:\.\d+)?)px",
                  lambda m: (f"{float(m.group(1)) * f:.1f}px" if float(m.group(1)) >= 6 else m.group(0)),
                  css)

def qss(theme: str = "dark") -> str:
    return _scale_px(_COMMON.format(**(_LIGHT if theme == "light" else _DARK)), UI_SCALE)
