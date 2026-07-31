"""Qt stylesheet for the TEagle 'assay terminal' look — a monospace-forward precision instrument,
matching the web UI (app/web/app.css): exact near-black + teal-mint palette, sharp 2px edges,
tabular-mono data, uppercase tracked micro-labels, and a teal accent reserved for interactive
affordances. Qt QSS has no letter-spacing / text-transform, so the wordmark spacing and uppercased
labels are applied in code (main.py).

Single source of colour truth for the whole app: the chrome palettes (_DARK/_LIGHT), every token
derived from them (ACCENT/TEXT/HEADRULE/FLAG) and the figure data hues (OKABE_ITO/OK/GENECOL/ARCHCOL/
CISCOL/GV_THEME/GELPAL, re-exported by figures.py) all live here — one definition per colour."""

import math, os, re

import fonts as _fonts

# Body/UI text is Roboto (bundled; see native/fonts.py); sequences, accessions and numeric data
# stay in Cascadia Mono for column alignment. SANS = prose/chrome, MONO = data/brand.
# The bundled families are NAMED ONCE, in fonts.py beside the .ttf list that ships them, and these
# stacks are built from those constants — so a change to what is bundled cannot leave the stylesheet
# asking for a font that is no longer in the build. Only the non-bundled OS fallbacks are literals.
MONO = f'"{_fonts.UI_FAMILY}", "{_fonts.WORDMARK_FAMILY}", "Consolas", "Courier New", monospace'
SANS = f'"{_fonts.BODY_FAMILY}", "Segoe UI", "Helvetica Neue", Arial, sans-serif'

# global UI zoom: scale fonts + padding uniformly (1px/2px borders & radii left untouched).
# UI_SCALE is now driven live (set at runtime + re-apply qss); the env var is only a first-paint seed.
UI_SCALE = float(os.environ.get("TEAGLE_UI_SCALE", "1.0"))

# spacing scale — the one rhythm for code-side margins/gaps. sp() SNAPS its argument onto SCALE, so the
# scale is enforced rather than merely declared: the 7/9/12/14 call sites land on 6/10/10/16 (each moves
# <=2px). 0 stays 0 (an explicit "no margin"); anything above L passes through, the ladder ends there.
S, M, L = 6, 10, 16
SCALE = (0, S, M, L)

def sp(token: float) -> int:
    """A spacing token snapped onto SCALE, then scaled to the current UI_SCALE (contentsMargins/setSpacing)."""
    v = token if token > L else min(SCALE, key=lambda s: (abs(s - token), -s))   # ties round up, never to 0
    return max(0, round(v * UI_SCALE))

# type scale: integer px only. Qt rounds a fractional QSS px size UP to the next whole pixel (10.5px
# renders exactly as 11px — measured), so the old 10.5/11.5/12.5 steps were shadows of 11/12/13 and are
# now written as those integers: same pixels, one step per size.
_COMMON = """
* {{ font-family: {sans}; font-size: 13px; color: {text}; }}
QMainWindow, QWidget#central, QDialog, QMessageBox {{ background: {bg}; }}
QScrollArea {{ border: none; background: {bg}; }}
QScrollArea > QWidget > QWidget {{ background: transparent; }}
QLabel {{ background: transparent; }}

/* header wordmark + chrome (mono, tracked) */
QLabel#word {{ font-family: {mono}; font-weight: 700; font-size: 16px; }}
QLabel#tagline {{ font-family: {sans}; color: {faint}; font-size: 10px; }}
QLabel#ver {{ font-family: {mono}; color: {faint}; font-size: 9px; font-weight: 600;
    border: 1px solid {line}; border-radius: 3px; padding: 1px 5px; }}
QLabel#statusTxt {{ font-family: {sans}; color: {dim}; font-size: 11px; }}
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
QPushButton#cardhdr {{ font-family: {sans}; background: {panel2}; border: none;
    border-bottom: 1px solid {line}; text-align: left; padding: 13px 16px; font-size: 12px; font-weight: 700; color: {text}; }}
QPushButton#cardhdr:hover {{ color: {accent}; }}
QPushButton#cardhdr:checked {{ border-bottom: none; }}

/* data-entry fields = mono, near-black well, sharp */
QLineEdit, QTextEdit, QPlainTextEdit {{ background: {bg}; border: 1px solid {line2}; border-radius: 2px;
    padding: 7px 9px; color: {text}; font-family: {mono}; font-size: 12px;
    selection-background-color: {accentsoft}; selection-color: {text}; }}
QComboBox {{ background: {bg}; border: 1px solid {line2}; border-radius: 2px; padding: 6px 9px;
    color: {text}; font-family: {sans}; font-size: 12px; }}
QLineEdit:focus, QTextEdit:focus, QComboBox:focus {{ border: 1px solid {accent}; }}
QComboBox::drop-down {{ border: none; width: 18px; }}
/* the popup's highlighted row is drawn by QComboBox's own delegate, which honours selection-background-color
   but ignores ::item rules and paints no focus rect — so outline here is inert either way (measured, Qt 6.11) */
QComboBox QAbstractItemView {{ background: {panel3}; color: {text}; border: 1px solid {line2};
    selection-background-color: {accentsoft}; selection-color: {text}; font-family: {sans}; outline: none; }}

/* buttons = mono, uppercase (applied in code), sharp 2px */
QPushButton {{ background: {panel2}; border: 1px solid {line2}; border-radius: 2px; padding: 8px 13px;
    color: {text}; font-family: {sans}; font-size: 12px; font-weight: 600; }}
QPushButton:hover {{ border: 1px solid {accent}; background: {panel}; }}
QPushButton:pressed {{ background: {panel3}; }}
QPushButton:disabled {{ color: {faint}; border: 1px solid {line}; background: {panel}; }}
QPushButton[primary="true"] {{ background: {accent}; color: {accentink}; border: 1px solid {accent}; font-weight: 700; }}
QPushButton[primary="true"]:hover {{ background: {accent2}; }}
QPushButton[primary="true"]:disabled {{ background: {panel2}; color: {faint}; border: 1px solid {line}; }}
QPushButton[sm="true"] {{ padding: 5px 9px; font-size: 11px; }}
/* in-table action button: the ONLY property set on it (never combined with sm -> no selector-specificity
   clash), compact enough to sit inside a table row without crowding it */
QPushButton[cellbtn="true"] {{ padding: 4px 8px; font-size: 11px; font-weight: 600; }}
QPushButton[link="true"] {{ background: transparent; border: none; color: {accent};
    font-family: {sans}; font-size: 11px; text-align: left; padding: 2px; font-weight: 600; }}
QPushButton[link="true"]:hover {{ color: {accent2}; }}
/* keyboard focus (web parity: app.css ::focus-visible). NOT `outline`: Qt paints outline at the CONTENTS rect
   (border box minus padding), so on the 2px-padded link buttons it lands on the glyphs and shears the text.
   Focus therefore recolours the 1px border and tints the fill — same border WIDTH as the base rule, so nothing
   shifts, and the tint keeps focus distinguishable from hover (which recolours the border only). */
QPushButton:focus {{ border: 1px solid {accent}; background: {accentsoft}; }}
/* primary sits on an accent fill: ring in accent ink. The fill MUST be restated — the generic :focus rule above
   has equal specificity and comes later, so without this it overwrites the accent fill with the tint and the
   accentink label drops to ~1.1:1 against it (measured), i.e. the label vanishes on click-focus. */
QPushButton[primary="true"]:focus {{ border: 1px solid {accentink}; background: {accent}; }}
/* link/cardhdr have border:none, so the generic :focus border above would ADD a box and shift them. Give link its
   own ring and pay for it out of the padding (1px border + 1px pad == the base 2px pad), so the box is unchanged
   and the ring sits outside the glyph area. cardhdr keeps its bottom rule and takes a tint only. */
QPushButton[link="true"]:focus {{ border: 1px solid {accent}; padding: 1px; background: {accentsoft}; }}
/* cardhdr is a StrongFocus button and the FIRST tab stop of every panel, so the tint alone (measured 1.09:1
   light / 1.21:1 dark against the unfocused header — invisible) is not a cue. It gets the same accent ring as
   every other button, paid for out of its 13/16px padding so the box and the glyph positions do not move, and
   keeps its bottom rule. */
QPushButton#cardhdr:focus {{ border: 2px solid {accent}; border-bottom: 1px solid {line};
    padding: 11px 14px 13px 14px; background: {accentsoft}; }}

/* metric readout gauges + key/value chrome */
QFrame#cell {{ background: {panel}; border: 1px solid {line}; border-radius: 2px; }}
QLabel#kdim {{ font-family: {mono}; color: {faint}; font-size: 9px; font-weight: 600; }}
QLabel#value {{ font-family: {mono}; font-size: 21px; font-weight: 600; color: {text}; }}
QLabel#value[state="good"] {{ color: {good}; }}
QLabel#value[state="bad"] {{ color: {bad}; }}

/* tables = mono data, uppercase tracked headers */
QTableWidget {{ background: {panel}; gridline-color: {line}; border: 1px solid {line}; border-radius: 2px;
    font-family: {mono}; font-size: 12px; alternate-background-color: {panel2}; outline: none; }}
QTableWidget::item {{ padding: {tpadv}px {tpadh}px; }}
QTableWidget::item:selected {{ background: {accentsoft}; color: {text}; }}
QHeaderView::section {{ background: {panel2}; color: {faint}; padding: 6px 8px; border: none;
    border-bottom: 1px solid {line}; font-family: {mono}; font-size: 10px; font-weight: 700; }}
QTableCornerButton::section {{ background: {panel2}; border: none; }}

/* classification banner (left accent border, big title) */
QFrame#classbn {{ background: {panel2}; border: 1px solid {line2}; border-left: 3px solid {accent}; border-radius: 2px; }}
QLabel#classbig {{ font-family: {sans}; font-size: 18px; font-weight: 700; color: {text}; }}
QLabel#classkls {{ font-family: {mono}; font-size: 11px; color: {dim}; }}
QLabel#classexp {{ font-family: {sans}; font-size: 12px; color: {text2}; }}
QLabel#cf {{ font-family: {sans}; font-size: 10px; font-weight: 700; padding: 3px 9px; border-radius: 2px; }}
QLabel#cf[level="High"] {{ background: {goodsoft}; color: {good}; }}
QLabel#cf[level="Moderate"] {{ background: {panel3}; color: {dim}; }}
QLabel#cf[level="Candidate"] {{ background: {warnsoft}; color: {warn}; }}

/* prose + misc chrome */
QLabel#orient {{ font-family: {sans}; color: {text2}; font-size: 11px; background: {panel2};
    border-left: 2px solid {line2}; border-radius: 2px; padding: 11px 14px; }}
QLabel#sectionlabel {{ font-family: {mono}; color: {faint}; font-size: 10px; font-weight: 700; }}
QLabel#cardmeta {{ font-family: {mono}; color: {dim}; font-size: 11px; }}
QLabel#empty {{ font-family: {sans}; color: {faint}; font-size: 12px; }}
QLabel#gvpos {{ font-family: {mono}; color: {faint}; font-size: 10px; }}
/* Level vocabulary, shared by the panel note (#errbanner) and the notification dialog (#notifmsg). TWO
   channels, because colour alone is a weak one (the gel learned this):
   1. hue — error red / warn amber / success green / info BLUE. info must NOT reuse the teal accent: light
      accentsoft and goodsoft composite to rgb(226,238,235) vs rgb(223,239,232), CAM02-UCS dE 1.5 (1.2 under
      deuteranomaly) — below a JND, so the severity carried no signal at all in the app's default theme.
   2. shape — the three OUTCOME levels (a run failed / finished with a caveat / finished clean) carry a 3px
      left rule, the same emphasis idiom as #classbn/#orient; a plain "info" NOTICE stays a flat 1px box.
      Every level restates border-left AND padding explicitly: the unsuffixed rule below is the error style,
      so whatever it declares leaks into any level that does not override it. The 2px the rule takes is paid
      back out of padding-left, so the ink starts at the same inset in all four. */
QLabel#errbanner {{ font-family: {sans}; background: {badsoft}; color: {bad}; border: 1px solid {bad};
    border-left: 3px solid {bad}; border-radius: 2px; padding: 9px 12px 9px 10px; font-size: 12px; }}
/* level variants below are the styling contract for a level-tagged banner; no caller sets it today */
QLabel#errbanner[level="success"] {{ background: {goodsoft}; color: {good}; border: 1px solid {good};
    border-left: 3px solid {good}; padding: 9px 12px 9px 10px; }}
QLabel#errbanner[level="warn"] {{ background: {warnsoft}; color: {warn}; border: 1px solid {warn};
    border-left: 3px solid {warn}; padding: 9px 12px 9px 10px; }}
QLabel#errbanner[level="info"] {{ background: {infosoft}; color: {info}; border: 1px solid {info};
    border-left: 1px solid {info}; padding: 9px 12px; }}
/* notification dialog — the old top-of-panel banner relocated into a small, closable dialog */
QDialog#notif {{ background: {panel}; }}
QLabel#notifmsg {{ font-family: {sans}; font-size: 13px; border-radius: 2px; padding: 12px 14px 12px 12px;
    background: {badsoft}; color: {bad}; border: 1px solid {bad}; border-left: 3px solid {bad}; }}
QLabel#notifmsg[level="success"] {{ background: {goodsoft}; color: {good}; border: 1px solid {good};
    border-left: 3px solid {good}; padding: 12px 14px 12px 12px; }}
QLabel#notifmsg[level="warn"] {{ background: {warnsoft}; color: {warn}; border: 1px solid {warn};
    border-left: 3px solid {warn}; padding: 12px 14px 12px 12px; }}
QLabel#notifmsg[level="info"] {{ background: {infosoft}; color: {info}; border: 1px solid {info};
    border-left: 1px solid {info}; padding: 12px 14px; }}
QLabel#queuerow {{ font-family: {sans}; color: {dim}; font-size: 11px; }}
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
QMenu {{ background: {panel3}; border: 1px solid {line2}; font-family: {sans}; font-size: 12px; }}
QMenu::item {{ padding: 7px 28px 7px 16px; }}
QMenu::right-arrow {{ width: 10px; margin-right: 8px; }}
QMenu::item:selected {{ background: {accentsoft}; color: {text}; }}
QToolTip {{ background: {panel3}; color: {text}; border: 1px solid {line2}; font-family: {sans}; font-size: 11px; padding: 4px 6px; }}
"""

_DARK = dict(bg="#0A0D10", panel="#0F1317", panel2="#131A1F", panel3="#19222A", line="#1E272E",
             line2="#2B3740", text="#E6EDF1", text2="#B4BEC5", dim="#8A959D", faint="#7E8A93",
             accent="#33D6B8", accent2="#1FB89C", accentink="#042420", accentsoft="rgba(51,214,184,0.14)",
             good="#40C088", goodsoft="rgba(64,192,136,0.16)", warn="#D8B368", warnsoft="rgba(216,179,104,0.16)",
             bad="#EC5F49", badsoft="rgba(236,95,73,0.14)",
             info="#58A6FF", infosoft="rgba(88,166,255,0.16)", mono=MONO, sans=SANS)
_LIGHT = dict(bg="#EDF1F3", panel="#FFFFFF", panel2="#F2F5F7", panel3="#E7ECEF", line="#DCE3E7",
              line2="#C6D0D6", text="#141B21", text2="#3A454D", dim="#57636B", faint="#5E6A72",
              accent="#0A7259", accent2="#086048", accentink="#FFFFFF", accentsoft="rgba(10,114,89,0.12)",
              good="#14774F", goodsoft="rgba(23,138,92,0.14)", warn="#7C621F", warnsoft="rgba(138,109,34,0.16)",
              bad="#BA3F2B", badsoft="rgba(198,67,46,0.12)",
              # info is the ONE level hue that is not the brand teal — see the level-vocabulary note above.
              # Blue also survives the red-green deficiencies: the info/success background pair holds
              # dE 7.6 (light) / 11.3 (dark) at full deuteranomaly AND protanomaly. J*=40.8 C=27.8 puts it
              # in the same lightness/chroma register as good/warn/bad. Re-measure both tracks + CVD before
              # changing it.
              info="#1A5FA8", infosoft="rgba(26,95,168,0.16)", mono=MONO, sans=SANS)

_PAL = {"dark": _DARK, "light": _LIGHT}
# every token below DERIVES from the palettes above — never restate a palette hex here
ACCENT = {k: p["accent"] for k, p in _PAL.items()}
TEXT = {k: p["text"] for k, p in _PAL.items()}         # wordmark "TE" ink, per theme
GOOD = {k: p["good"] for k, p in _PAL.items()}         # per-theme "ready/ok" status ink (WCAG-tuned both tracks)
BAD = {k: p["bad"] for k, p in _PAL.items()}           # per-theme "error/failed" status ink
# header accent underline (teal fading to transparent), like the web header::after
HEADRULE = {k: f"qlineargradient(x1:0, y1:0, x2:1, y2:0, stop:0 {c}, stop:0.42 transparent)"
            for k, c in ACCENT.items()}
# per-theme flag colours for QC ΔG cells — reuse the WCAG-tuned bad/warn palette (dark vs light), so amber/red
# read correctly on both backgrounds instead of a single hardcoded dark-tuned hex
FLAG = {k: {"warn": p["bad"], "caution": p["warn"]} for k, p in _PAL.items()}

# ============ figure data hues (re-exported by figures.py — one definition per colour) ============
# Okabe-Ito colour-blind-safe base (Okabe & Ito 2008). Every figure hue below NAMES one of these rather
# than restating its hex, so a shared hue exists in exactly one place.
OKABE_ITO = {"orange": "#E69F00", "skyblue": "#56B4E9", "green": "#009E73",
             "blue": "#0072B2", "vermillion": "#D55E00", "purple": "#CC79A7"}
_OI = OKABE_ITO
# domain / feature hues — must match the web figure bands (app/web/app.js).
# ENV needs its own hue — without a key it fell back to grey beside GAG's grey. Okabe-Ito yellow #F0E442 was
# tried and rejected: 1.24:1 on the light track (the #888 fallback it replaced managed 3.33:1). #B15928 is the
# measured pick — 3.72:1 dark / 4.56:1 light, and the only candidate that holds up under CVD (worst separation
# 9.5 CAM02-UCS vs 1.9–2.8 for the alternatives). Re-measure both tracks AND deuteranomaly/protanomaly before
# changing it.
# GAG likewise had to leave grey: #7A7A7A sat only 5.63 CAM02-UCS from the #888 no-entry fallback, and both being
# achromatic, CVD moved that by 0.00 — a named domain read as an unnamed one. #6438FC is the measured replacement:
# 3.02:1 on the dark track / 5.61:1 on the light track, worst-case separation 18.4 CAM02-UCS from any required hue
# (#888, RT, INT/TIR, RNaseH, PR, CHR/TPase, LTR, ENV) and 13.8 once PBS/PPT/ORF and the gene-model bands are
# included, under normal, deuteranomalous and protanomalous vision alike. Magenta (#C71585 10.6, #B03A8E 10.1) and
# brown (#8C564B 13.0, but nearest neighbour is ENV) were measurably worse. Re-measure before changing it.
OK = {"RT": _OI["blue"], "INT": _OI["orange"], "RNaseH": _OI["green"], "PR": _OI["purple"], "GAG": "#6438FC",
      "ENV": "#B15928", "CHR": _OI["vermillion"], "TPase": _OI["vermillion"], "LTR": _OI["skyblue"],
      "TIR": _OI["orange"], "tail": _OI["purple"], "ORF": "#4C6C97",
      "on": _OI["green"], "off": _OI["vermillion"], "ladder": "#999999"}
# exon_derived = a lighter tint of the annotated-exon green: same family, reads "inferred, not annotated";
# gap kept distinct (light) from the darker slate flank so filler regions aren't a single ambiguous colour.
GENECOL = {"exon": _OI["green"], "exon_derived": "#7fd3b8", "intron": "#8792a0",
           "cds": _OI["vermillion"], "flank": "#5b6b7a", "gap": "#c3ccd6"}
# retroviral transcript architecture (ERV): env exons green, the removed gag-pro-pol span amber-brown so the
# "single large intron = frameshift-fused polyprotein" reads distinctly from a grey host intron.
ARCHCOL = {"exon": _OI["green"], "intron": "#B0752E"}
# LTR cis-elements: PBS (leader, purple) and PPT (before 3' LTR, blue) — each distinct from the LTR blocks,
# the env-exon green and the intron amber.
# PAS (the advisory polyA-signal motif) needed a third hue that no other band in a shared render already
# owns. #B5316B is the measured pick: worst-case separation 20.6 CAM02-UCS from every required hue (PBS,
# PPT, LTR, tail/PR, ORF, GAG, ENV, RT, INT, RNaseH, TPase, the gene-model bands and the #888 fallback),
# holding 10.1 under deuteranomaly and 12.1 under protanomaly — the best CVD floor of the candidates
# measured — at 3.26:1 on the dark track and 5.81:1 on the light one. Re-measure both tracks AND
# deuteranomaly/protanomaly before changing it.
CISCOL = {"PBS": "#8459C4", "PPT": "#2C7FB8", "PAS": "#B5316B"}

# genome-viewer chrome per render target. 'export' is a print-neutral scheme (transparent paper), NOT the app's
# light theme. Only values that genuinely equal a chrome token derive from one; figure-specific greys stay
# literal here, because snapping them onto the nearest palette token would silently repalette the figure.
_GV_WIN = "#1f6feb"                                    # viewport-window blue, on light/neutral paper
GV_THEME = {
    "export": {"paper": "none", "ink": "#222", "faint": "#555", "grid": "#dcdfe3",
               "track": "#00000000", "lane": "#0000000d", "frame": "#c7ccd2", "win": _GV_WIN},
    "white":  {"paper": _LIGHT["panel"], "ink": _LIGHT["text"], "faint": "#5a6570", "grid": "#eceef1",
               "track": "#f6f8fa", "lane": "#eef1f4", "frame": "#dde1e6", "win": _GV_WIN},
    "dark":   {"paper": "#0b1016", "ink": _DARK["text"], "faint": _DARK["dim"], "grid": "#182029",
               "track": "#10171e", "lane": "#121b23", "frame": "#243039", "win": "#4aa8ff"},
}

# agarose-gel palettes. "site" = a NEUTRAL colour for a whole-genome scan with no design locus: the products are
# neither on- nor off-target, just genomic priming sites, so they must not read as the off-target warning colour.
GELPAL = {
    "transparent": {"paper": "none", "gel": "#0f1316", "well": "#04060a", "stroke": "#2a3138", "ink": "#5a656f",
                    "on": OK["on"], "off": OK["off"], "single": _OI["blue"], "site": _OI["skyblue"],
                    "ladder": OK["ladder"], "glow": 1.4, "band": 2.6},
    "dark":        {"paper": "#0b0e11", "gel": "#0f1316", "well": "#04060a", "stroke": "#232a30", "ink": "#8792a0",
                    "on": OK["on"], "off": OK["off"], "single": _OI["blue"], "site": _OI["skyblue"],
                    "ladder": OK["ladder"], "glow": 1.4, "band": 2.6},
    # light theme keeps near-black vs dark red deliberately: measured against Machado-2009 and Viénot-Brettel-Mollon,
    # this pair separates by 27 L* (greyscale 2.27:1) and beats an Okabe-Ito green/vermillion pair here, which goes
    # near-isoluminant under protanopia on a near-white gel. Do not "fix" it to brand hues without re-measuring.
    "white":       {"paper": "#ffffff", "gel": "#ededed", "well": "#c4c4c4", "stroke": "#cccccc", "ink": "#555555",
                    "on": "#151515", "off": "#992222", "single": "#1f5fa8", "site": "#2a6f97",
                    "ladder": "#9a9a9a", "glow": 0.3, "band": 2.6},
    "uv":          {"paper": "#050310", "gel": "#0a0714", "well": "#000000", "stroke": "#1c1236", "ink": "#9fb4d8",
                    "on": "#5bff6b", "off": "#ffcf47", "single": "#6fb2ff", "site": "#79d0ff",
                    "ladder": "#79d0ff", "glow": 3.2, "band": 3.1},
    "mono":        {"paper": "#0d0d0d", "gel": "#181818", "well": "#000000", "stroke": "#2b2b2b", "ink": "#b2b2b2",
                    "on": "#f2f2f2", "off": "#9a9a9a", "single": "#6f6f6f", "site": "#c0c0c0",
                    "ladder": "#cfcfcf", "glow": 2.0, "band": 2.9},
}


def _scale_px(css: str, f: float) -> str:
    if abs(f - 1.0) < 1e-3:
        return css
    # scale font/padding px (>=6px) only; keep 1-2px borders and 2px radii crisp
    return re.sub(r"(\d+(?:\.\d+)?)px",
                  lambda m: (f"{float(m.group(1)) * f:.1f}px" if float(m.group(1)) >= 6 else m.group(0)),
                  css)

# QTableWidget::item padding, as UNSCALED px. Qt applies this padding to setCellWidget widgets too — it both
# offsets them and shrinks their usable box — so any code putting a widget in a cell must budget for it. Single
# source of truth: the QSS below is built from these, and callers read them instead of re-hardcoding 6/9.
TABLE_ITEM_PAD_V, TABLE_ITEM_PAD_H = 6, 9


def qss(theme: str = "dark") -> str:
    pal = _LIGHT if theme == "light" else _DARK
    css = _scale_px(_COMMON.format(**pal, tpadv=TABLE_ITEM_PAD_V, tpadh=TABLE_ITEM_PAD_H), UI_SCALE)
    # Current-cell keyboard cue. MUST be ::item:focus — a view-level `outline` paints 0px under the windows11
    # style the app actually ships on (its CE_ItemViewItem never issues PE_FrameFocusRect); it only worked on
    # fusion/vista. Independent of selection, so NoSelection tables (the genome manager) still get a cue.
    # The 2px ring is paid for out of the item padding, EXACTLY at every UI scale (columns are sized
    # ResizeToContents = widest text + padding, with zero slack, so a ring that added width would elide the
    # focused cell's own text). Appended post-_scale_px so these already-scaled paddings are not scaled twice.
    fpv = max(0, math.ceil(TABLE_ITEM_PAD_V * UI_SCALE) - 2)   # ceil: Qt rounds a fractional QSS px UP
    fph = max(0, math.ceil(TABLE_ITEM_PAD_H * UI_SCALE) - 2)
    css += (f"\nQTableWidget::item:focus {{ border: 2px solid {pal['accent']}; "
            f"padding: {fpv}px {fph}px; }}\n")
    return css
