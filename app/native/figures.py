"""Publication figure builders: pure functions returning SVG strings. The native app renders them on-screen
via QSvgRenderer and exports the identical string as SVG / rasterises it to PNG / writes it to a vector PDF.

Figures: an interactive genome viewer (ruler + tracks, windowed), a to-scale multi-lane agarose gel (log MW
axis, 5 palettes), the self-similarity dot plot and heat map, and the whole-genome repeat-composition bars."""
from __future__ import annotations
import colorsys
import functools
import math
import re

# every colour used here is defined once, in theme.py (single source of truth); re-exported so
# `figures.OK` / `figures.GELPAL` stay valid for existing callers and tests.
from theme import OKABE_ITO, OK, GENECOL, ARCHCOL, CISCOL, GV_THEME, GELPAL   # noqa: F401
from teagle_core.gff3 import bound_public_record_label, is_neutral_tdr, public_record_label

FIGFONT = "Cascadia Mono, Consolas, monospace"   # bundled UI font; mono digits align on the gel ladder / ruler


def esc(s) -> str:
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def _rgb(hex_color: str):
    """(r, g, b) of a 3- or 6-digit hex (a bare '#888' fallback would otherwise raise and abort the whole render)."""
    h = hex_color.lstrip("#")
    if len(h) == 3:
        h = "".join(c * 2 for c in h)                         # #888 -> #888888
    return int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)


def _lum(hex_color: str) -> float:
    """WCAG relative luminance."""
    def lin(v):
        v /= 255.0
        return v / 12.92 if v <= 0.03928 else ((v + 0.055) / 1.055) ** 2.4
    r, g, b = _rgb(hex_color)
    return 0.2126 * lin(r) + 0.7152 * lin(g) + 0.0722 * lin(b)


def _contrast(a: str, b: str) -> float:
    la, lb = _lum(a), _lum(b)
    return (max(la, lb) + 0.05) / (min(la, lb) + 0.05)


def _label_ink(hex_color: str) -> str:
    """Pick black/white label by best WCAG contrast against a fill colour."""
    L = _lum(hex_color)
    return "#fff" if (1.05 / (L + 0.05)) >= ((L + 0.05) / 0.05) else "#111"


# Figure text and glyph rules (one definition for every figure here). Small text is 9 px on screen and 10 px in a
# file that leaves the app; a glyph is never narrower than MIN_GLYPH_PX (an 11 bp TIR on a 4.5 kb element is 2 px
# wide at fit); a white-paper glyph gets a 1 px outline of its own hue darkened to OUTLINE_MIN_CR against the
# paper. WCAG 1.4.11 asks 3:1 for graphical objects; 4.5 is a deliberate margin (a 1 px anti-aliased stroke
# loses contrast), not a standard value.
FIG_PX_SCREEN = 9
FIG_PX_EXPORT = 10
MIN_GLYPH_PX = 3
OUTLINE_MIN_CR = 4.5


@functools.lru_cache(maxsize=256)
def _outline(fill: str, bg: str) -> str:
    """The fill's own hue, darkened (HLS lightness only) until it holds OUTLINE_MIN_CR against `bg`; the fill itself
    when it already does."""
    if _contrast(fill, bg) >= OUTLINE_MIN_CR:
        return fill if len(fill.lstrip("#")) == 6 else "#%02x%02x%02x" % _rgb(fill)
    r, g, b = (v / 255.0 for v in _rgb(fill))
    h, l, s = colorsys.rgb_to_hls(r, g, b)
    while l > 0.0:
        l = max(0.0, l - 0.005)
        c = "#%02x%02x%02x" % tuple(round(v * 255) for v in colorsys.hls_to_rgb(h, l, s))
        if _contrast(c, bg) >= OUTLINE_MIN_CR:
            return c
    return "#000000"


def _fmt_int(n) -> str:
    return f"{int(round(n)):,}"


def fmt_interval(a, b) -> str:
    """The one coordinate form shared by the tables and every figure caption: plain digits, 'start–end'."""
    return f"{int(a)}–{int(b)}"


def _text_w(s: str, px: float) -> float:
    """Upper-bound width of a mono string (the figure font's advance is 0.6 em; the fallback is narrower)."""
    return len(s) * px * GV_CHAR_EM


def _wrap_mono(items, width_px: float, px: float, sep: str = " · "):
    """Pack text items into lines no wider than `width_px`; an item is never dropped or cut, and one wider than the
    line sits alone on its own line."""
    lines, cur = [], ""
    for it in items:
        cand = it if not cur else cur + sep + it
        if cur and _text_w(cand, px) > width_px:
            lines.append(cur)
            cur = it
        else:
            cur = cand
    if cur:
        lines.append(cur)
    return lines


# ================= genome viewer =================
GV_NAME_PX = 9.5      # track-name font size in the viewer
GV_CHAR_EM = 0.6      # Cascadia Mono advance width (measured 0.600 em; the SVG renderer draws no wider)
GV_NAME_GAP = 10      # gap between the end of a track name and the plot area
GV_NAME_PAD = 4       # breathing room at the figure's left edge
GV_ML_MIN = 96        # gutter floor, so short-named models keep the historical layout
GV_MR = 16            # right margin


def _strand_mark_w(name_px: float) -> float:
    """Width reserved left of a stranded track's plot for its '+' / '−' row marks."""
    return name_px * GV_CHAR_EM + 4


def gv_gutter(model: dict, W: float = None, name_px: float = GV_NAME_PX) -> float:
    """Width of the left track-label gutter: sized to the longest track name in this model, so a
    name like 'env mRNA (predicted)' is never clipped by the figure's left edge. Mono metric, so the
    result is an upper bound on the rendered text width. A stranded track also reserves room for its
    '+' / '−' row marks. With W, the gutter is capped so the plot area keeps its 120 px floor on a very
    narrow figure. One definition — svg_genome draws it and the interactive canvas maps pixels to bp
    with it, so the two can never disagree. `name_px` is the name font size (larger in an export)."""
    w = max((len(t.get("name", "")) * name_px * GV_CHAR_EM + (_strand_mark_w(name_px) if t.get("stranded") else 0)
             for t in (model or {}).get("tracks", []) if t.get("name")), default=0)
    g = float(max(GV_ML_MIN, math.ceil(w) + GV_NAME_GAP + GV_NAME_PAD))
    return g if W is None else min(g, max(float(GV_ML_MIN), W - GV_MR - 120))


def gv_tick_label(bp: float) -> str:
    return f"{bp/1e6:g}M" if bp >= 1e6 else (f"{bp/1000:g}k" if bp >= 1000 else f"{bp:g}")


def gv_nice_step(span: float, ticks: int) -> float:
    raw = max(1.0, span / ticks)
    mag = 10 ** math.floor(math.log10(raw))
    n = raw / mag
    return (1 if n < 1.5 else 2 if n < 3.5 else 5 if n < 7.5 else 10) * mag


_DOMAIN_NAME = {"RT": "reverse transcriptase", "INT": "integrase", "RNaseH": "RNase H", "PR": "protease",
                "GAG": "gag capsid / matrix", "ENV": "envelope", "CHR": "chromodomain", "TPase": "transposase",
                "EN": "endonuclease (APE)", "ORF1": "LINE ORF1p", "HEL": "Helitron helicase",
                "YR": "tyrosine recombinase"}


def _domain_orf_map(rec: dict) -> dict:
    """{ORF index in rec['orfs']: [domain codes detected inside it]}. scan_domains numbers its ORFs from the same
    find_orfs list as rec['orfs']; a hit is attributed only when its strand matches that ORF and its nucleotide
    interval lies inside it. A hit failing either check is attributed to nothing (fail closed), so an outline is
    never drawn on the wrong ORF."""
    orfs = rec.get("orfs") or []
    out = {}
    for d in rec.get("domains") or []:
        i, nt = d.get("orf"), d.get("nt")
        if not isinstance(i, int) or not 0 <= i < len(orfs) or not nt:
            continue
        o = orfs[i]
        if d.get("strand") != o.get("strand") or nt[0] < o["start"] or nt[1] > o["end"]:
            continue
        codes = out.setdefault(i, [])
        if d.get("domain") not in codes:
            codes.append(d.get("domain"))
    return out


def gv_tracks_from_rec(rec: dict) -> dict:
    """Build a genome-viewer model from an analysis record (structural + domains + ORFs + ERV architecture).

    Feature keys (a track is {name, height, features[, stranded]}): start, end (0-based half-open), color, label,
    tip, legend (the export legend text); and, where they apply, kind, strand, intron, derived, outline ('sub'
    for an advisory pair drawn open), emph (an ORF that carries a detected domain), name and, for a domain hit,
    protein, dna and index (its position in rec['domains']). svg_genome copies kind, strand, protein, dna,
    index and derived into the hit regions it returns, which the right-click menu reads."""
    tracks, reps, cis, tails, tsds = [], [], [], [], []
    tdr = is_neutral_tdr(rec.get("classification"))          # a withheld candidate is a TDR, not an LTR claim
    for e in rec.get("structural", []):
        t = e["type"]
        if t.startswith("LTR") or t.startswith("TIR"):
            col = OK["LTR"] if t.startswith("LTR") else OK["TIR"]
            for p in (e.get("five_prime"), e.get("three_prime")):
                if p:
                    kind = t.split(" ")[0]
                    reps.append({"start": p[0], "end": p[1], "color": col,
                                 "label": "TDR" if (tdr and kind == "LTR") else kind, "kind": kind,
                                 "legend": "TDR (terminal direct repeat)" if (tdr and kind == "LTR") else kind,
                                 "tip": f"{t} {p[0]}–{p[1]}"})
        elif t.startswith("Sub-threshold"):
            # an advisory pair below the identity floor: drawn open and labelled so it cannot read as an accepted
            # LTR. Its kind starts with 'LTR' so the right-click menu treats it as a repeat (no primer design).
            for p in (e.get("five_prime"), e.get("three_prime")):
                if p:
                    reps.append({"start": p[0], "end": p[1], "color": OK["LTR"], "label": "sub-threshold",
                                 "kind": "LTR sub-threshold", "outline": "sub", "legend": "sub-threshold repeat pair",
                                 "tip": f"{t} {p[0]}–{p[1]} · below the {e.get('threshold', '')}% identity floor, "
                                        f"not accepted as an LTR"})
        elif t.startswith("TSD"):
            for side, p in (("upstream", e.get("upstream")), ("downstream", e.get("downstream"))):
                if p:
                    tsds.append({"start": p[0], "end": p[1], "color": GENECOL["flank"], "label": "TSD", "kind": "TSD",
                                 "legend": "TSD (target-site duplication)",
                                 "tip": f"TSD {side} copy {p[0]}–{p[1]} · {e.get('length')} bp {e.get('motif', '')}"
                                        + (" · at chance level" if e.get("chance_level") else "")})
        elif t.startswith("polyA-signal"):
            # an LTR cis-element, NOT the LINE poly-A tail — it must not land in the tail track below, where
            # a motif would read as a detected tail. Labelled with the hexamer so the row is self-describing.
            p = e["pos"]
            cis.append({"start": p[0], "end": p[1], "color": CISCOL["PAS"], "kind": "PAS",
                        "label": "PAS·" + (e.get("motif") or "?"), "legend": "PAS (polyA-signal motif)",
                        "tip": (f"polyA-signal motif {e.get('motif')} ({e.get('variant')}) {p[0]}–{p[1]} · "
                                f"downstream element {e.get('dse_state')} · advisory: a motif, not a cleavage site")})
        elif t.startswith("PBS") or t.startswith("PPT"):
            p = e["pos"]; key = t[:3]
            if key == "PBS":                              # name the tRNA only when confident, else "PBS·?"
                lab = "PBS·" + ((e.get("priming_trna") or "").replace("tRNA-", "") if e.get("confident") else "?")
                tip = (f"{t} {p[0]}–{p[1]} · " + (f"{e.get('priming_trna')} {e.get('identity')}%"
                       if e.get("confident") else
                       f"priming tRNA undetermined (closest {e.get('best_match','?')} {e.get('identity')}%)"))
            else:
                lab = "PPT"
                tip = f"{t} {p[0]}–{p[1]} · {int(round(e.get('purine_frac', 0) * 100))}% purine"
            cis.append({"start": p[0], "end": p[1], "color": CISCOL[key], "label": lab, "kind": key, "legend": key,
                        "tip": tip})
        elif e.get("pos"):
            p = e["pos"]
            f = {"start": p[0], "end": p[1], "color": OK["tail"], "label": t.split(" ")[0],
                 "kind": t.split(" ")[0], "tip": f"{t} {p[0]}–{p[1]}"}
            if t.startswith("poly"):
                f["legend"] = "poly-A/T tail"
            (tails if t.startswith("poly") else reps).append(f)   # a poly-A/T tail is not a terminal repeat
    if reps:
        tracks.append({"name": "terminal repeats", "height": 20, "features": reps})
    if tsds:
        tracks.append({"name": "TSD", "height": 18, "features": tsds})
    if tails:
        tracks.append({"name": "poly-A/T tail", "height": 18, "features": tails})
    if cis:
        tracks.append({"name": "cis-elements", "height": 18, "features": cis})
    doms = []
    for di, d in enumerate(rec.get("domains", [])):
        code = d["domain"]
        strand = d.get("strand") or "+"
        doms.append({"start": d["nt"][0], "end": d["nt"][1], "color": OK.get(code, "#888"), "label": code,
                     "kind": "domain", "strand": strand, "protein": d.get("protein"), "dna": d.get("dna"),
                     "index": di, "name": d.get("label", ""),
                     "legend": f"{code} {_DOMAIN_NAME.get(code, d.get('label', ''))}".strip(),
                     "tip": (f"{code} · {d.get('label','')} · nt {d['nt'][0]}–{d['nt'][1]} · strand {strand} · "
                             f"score {d.get('score')}")})
    if doms:
        tracks.append({"name": "protein domains", "height": 22, "features": doms})
    arch = rec.get("retroviral")
    if arch:                                                  # ERV: env expressed from a spliced subgenomic mRNA
        feat = [{"start": arch["leader_exon"][0], "end": arch["leader_exon"][1], "color": ARCHCOL["exon"],
                 "label": "leader", "kind": "exon", "legend": "env mRNA exon",
                 "tip": f"5′ leader exon {arch['leader_exon'][0]}–{arch['leader_exon'][1]}"},
                {"start": arch["intron"][0], "end": arch["intron"][1], "color": ARCHCOL["intron"], "intron": True,
                 "label": "gag–pol", "kind": "intron", "legend": "gag–pol interval (env intron)",
                 "tip": f"gag–pro–pol coding interval removed from env RNA {arch['intron'][0]}–{arch['intron'][1]}"},
                {"start": arch["env_exon"][0], "end": arch["env_exon"][1], "color": ARCHCOL["exon"],
                 "label": "env", "kind": "exon", "legend": "env mRNA exon",
                 "tip": f"env exon {arch['env_exon'][0]}–{arch['env_exon'][1]}"}]
        tracks.append({"name": "env mRNA (predicted)", "height": 20, "features": feat})
    carrying = _domain_orf_map(rec)
    orfs = []
    for oi, o in enumerate(rec.get("orfs", [])):
        codes = carrying.get(oi, [])
        orfs.append({"start": o["start"], "end": o["end"], "color": OK["ORF"], "strand": o["strand"], "kind": "ORF",
                     "index": oi, "emph": bool(codes), "domains": codes,
                     "legend": "ORF carrying a detected domain" if codes else "ORF",
                     "tip": f"ORF {o['strand']}{o['frame']} · {o['length_aa']} aa"
                            + (f" · carries {', '.join(codes)}" if codes else "")})
    if orfs:
        tracks.append({"name": "ORFs", "height": 26, "features": orfs, "stranded": True})
    return {"length": rec.get("composition", {}).get("length", 1) or 1, "tracks": tracks}


def _flanks_and_gaps(exons: list, introns: list, length: int) -> list:
    """Flanking + inter-feature regions that are neither exon nor intron: the 5' upstream flank
    (0 .. first feature), the 3' downstream flank (last feature .. length), and any interior gap not
    covered by an exon or intron. Returned as clickable region dicts so the user can copy/design there
    too. Labels avoid spaces/apostrophes so they read cleanly in a FASTA header."""
    spans = [(f["start"], f["end"]) for f in (exons + introns)]
    if not spans:
        return []
    lo = min(s for s, _ in spans); hi = max(e for _, e in spans)
    out = []
    if lo > 0:
        out.append({"start": 0, "end": lo, "label": "5prime_flank", "kind": "flank", "name": "5′ flank"})
    if length and hi < length:
        out.append({"start": hi, "end": length, "label": "3prime_flank", "kind": "flank", "name": "3′ flank"})
    covered = sorted(spans)                                    # interior gaps = holes in the exon∪intron cover
    cur = lo
    for s, e in covered:
        if s > cur:
            out.append({"start": cur, "end": s, "label": "gap", "kind": "gap", "name": "gap"})
        cur = max(cur, e)
    return [r for r in out if r["end"] > r["start"]]           # skip zero/negative-length regions


def gv_tracks_from_gene(gm: dict, length: int, include_flanks: bool = False) -> dict:
    tracks = []
    # a CDS-implied exon that the record does NOT annotate is marked distinctly (lighter green + 'exon*' +
    # tip) so a tool-inferred coordinate is never mistaken for a GenBank-annotated exon.
    feat = [{"start": e["start"], "end": e["end"],
             "color": GENECOL["exon_derived"] if e.get("derived") else GENECOL["exon"],
             "label": "exon*" if e.get("derived") else "exon", "kind": "exon", "strand": e.get("strand") or "+",
             "derived": bool(e.get("derived")),
             "legend": "exon* (derived from the CDS)" if e.get("derived") else "exon",
             "tip": f"exon {e['start']}–{e['end']} ({e['end']-e['start']} bp)"
                    + (" · derived from the record's CDS/mRNA — not a separate exon annotation" if e.get("derived") else "")}
            for e in gm.get("exons", [])]
    # a derived intron is the gap between two exons of one unit, not an annotated intron: marked 'intron*' so a
    # tool-inferred coordinate is never read as a GenBank annotation (the same rule as a derived exon)
    feat += [{"start": i["start"], "end": i["end"], "color": GENECOL["intron"], "intron": True, "kind": "intron",
              "strand": i.get("strand") or "+", "derived": bool(i.get("derived")),
              "label": "intron*" if i.get("derived") else "intron",
              "legend": "intron* (derived from the exons)" if i.get("derived") else "intron",
              "tip": f"intron {i['start']}–{i['end']}" + (
                  " · derived from the gap between exons — not a separate intron annotation" if i.get("derived") else "") + (
                  f" · {i['donor']}…{i['acceptor']}{' (canonical)' if i.get('canonical') else ''}"
                  if i.get("donor") else "")} for i in gm.get("introns", [])]
    if include_flanks:                                        # 5'/3' flanks + interior gaps, clickable for copy/design
        for r in _flanks_and_gaps(gm.get("exons", []), gm.get("introns", []), length):
            feat.append({"start": r["start"], "end": r["end"],
                         "color": GENECOL["gap"] if r["kind"] == "gap" else GENECOL["flank"],
                         "label": r["name"],                  # readable on-glyph name; FASTA id is sanitised in _feat_menu
                         "kind": r["kind"], "legend": r["kind"],
                         "tip": f"{r['name']} {r['start']}–{r['end']} ({r['end']-r['start']} bp) · not exon/intron"})
    if feat:
        tracks.append({"name": "exons / introns", "height": 22, "features": feat})
    cds = [{"start": c["start"], "end": c["end"], "color": GENECOL["cds"], "label": "CDS", "kind": "CDS",
            "strand": c.get("strand") or "+", "legend": "CDS",
            "tip": f"CDS {c['start']}–{c['end']}"} for c in gm.get("cds", [])]
    if cds:
        tracks.append({"name": "CDS (coding)", "height": 16, "features": cds})
    return {"length": length or 1, "tracks": tracks}


def _gv_theme(theme: str, for_export: bool) -> dict:
    """Genome-viewer chrome for a render target; the palettes live in theme.GV_THEME (copy, so a
    caller mutating the returned dict cannot corrupt the shared one)."""
    return dict(GV_THEME["export" if for_export else ("white" if theme == "white" else "dark")])


def _tick_step(span: float, max_ticks: int) -> float:
    """Smallest 1/2/5 x 10^n step that puts at most `max_ticks` ticks across `span`."""
    raw = max(1.0, span / max(max_ticks, 1))
    mag = 10 ** math.floor(math.log10(raw))
    for m in (1, 2, 5, 10):
        if m * mag >= raw:
            return m * mag
    return 10 * mag


def _gv_legend_entries(plans):
    """[(fill, text, style)] for every colour the visible features draw, in track order, once per (colour, text).
    style: 'plain', 'emph' (an ORF with a detected domain, ink outline) or 'open' (an advisory pair, drawn open)."""
    seen, out = set(), []
    for pl in plans:
        for f, _a, _b, _inside, _arrow in pl["items"]:
            text = f.get("legend") or f.get("label")
            if not text:
                continue
            key = (f["color"], text)
            if key in seen:
                continue
            seen.add(key)
            out.append((f["color"], text, "emph" if f.get("emph") else ("open" if f.get("outline") == "sub" else "plain")))
    return out


def svg_genome(model: dict, view: dict, W: float, theme: str, for_export: bool = False, return_regions: bool = False,
               legend: bool = False):
    """Render the genome viewer SVG. With return_regions=True, also return a list of on-screen
    feature hit-boxes (SVG coords + feature identity) so the interactive canvas can hover/right-click.

    Each region carries x0, y0, x1, y1 (SVG coords), start, end, label, tip, strand and kind, plus protein, dna,
    index, derived and domains when the feature has them. Small text is FIG_PX_SCREEN px (FIG_PX_EXPORT in an
    export); a glyph is at least MIN_GLYPH_PX wide; on white paper (the light viewer or any export) each glyph has a
    1 px outline of its own hue. A feature too narrow for its label is named in a wrapped caption below its track
    (nothing is truncated). legend=True adds a legend of every colour drawn and the coordinate convention."""
    L = model["length"]
    regions = []
    px = FIG_PX_EXPORT if for_export else FIG_PX_SCREEN        # small text: ticks, glyph labels, captions
    name_px = FIG_PX_EXPORT if for_export else GV_NAME_PX      # track names
    MR, ovH, rulerH, MT = GV_MR, 13, 24, 34
    ML = gv_gutter(model, W, name_px)               # gutter fits the widest track name (no clipped label)
    plotW = max(120, W - ML - MR)
    s0, s1 = view["start"], view["end"]
    span = max(1, s1 - s0)
    def bx(bp): return ML + (bp - s0) / span * plotW
    def ox(bp): return ML + bp / L * plotW
    T = _gv_theme(theme, for_export)
    light = for_export or theme == "white"                    # outlined glyphs on a light ground
    ground = "#FFFFFF" if for_export else T["track"]          # what a glyph is read against
    # pass 1: each visible feature's drawn box (never narrower than MIN_GLYPH_PX, kept inside the plot) and, for a
    # label that does not fit inside it, a caption item
    plans = []
    for t in model["tracks"]:
        items, caps = [], []
        for f in t["features"]:
            a = max(bx(f["start"]), ML); b = min(bx(f["end"]), ML + plotW)
            if b < ML - 0.5 or a > ML + plotW + 0.5:
                continue
            if not f.get("intron") and b - a < MIN_GLYPH_PX:
                c = (a + b) / 2
                a, b = c - MIN_GLYPH_PX / 2, c + MIN_GLYPH_PX / 2
                if a < ML:
                    a, b = ML, ML + MIN_GLYPH_PX
                if b > ML + plotW:
                    a, b = ML + plotW - MIN_GLYPH_PX, ML + plotW
            # a domain glyph wide enough also carries an arrowhead at its coding end: Copia (INT->RT) versus Gypsy
            # (RT->INT) is read from the ORDER of domains along the element's own strand, so the direction is drawn
            arrow = 8 if (f.get("kind") == "domain" and f.get("strand") in ("+", "-") and b - a >= 16) else 0
            inside = bool(f.get("label")) and not f.get("intron") and _text_w(f["label"], px) + 6 + arrow <= b - a
            if f.get("label") and not f.get("intron") and not inside:
                caps.append(f'{f["label"]} {fmt_interval(f["start"], f["end"])}')
            items.append((f, a, b, inside, arrow))
        plans.append({"t": t, "items": items, "lines": _wrap_mono(caps, plotW, px)})
    y = MT + ovH + 10 + rulerH
    track_ys = []
    for pl in plans:
        track_ys.append(y)
        cap_h = (px + 3 + (len(pl["lines"]) - 1) * (px + 2)) if pl["lines"] else 0
        y += pl["t"].get("height", 20) + max(20, cap_h + 10)
    tracks_end = y
    entries = _gv_legend_entries(plans) if legend else []
    leg_rows, leg_pos, lx0 = 0, [], GV_NAME_PAD + 4
    if legend:
        x, row = lx0, 0
        for fill, text, style in entries:
            wi = 9 + 4 + _text_w(text, px) + 14
            if x > lx0 and x + wi > W - GV_MR:
                x, row = lx0, row + 1
            leg_pos.append((x, row))
            x += wi
        leg_rows = row + 2                                     # the entries, then the coordinate convention
    leg_row_h = px + 6
    H = tracks_end + 12 + (leg_rows * leg_row_h + 4 if legend else 0)
    s = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" '
         f'font-family="{FIGFONT}">')
    s += (f'<defs><linearGradient id="gvwin" x1="0" x2="0" y1="0" y2="1">'
          f'<stop offset="0" stop-color="{T["win"]}" stop-opacity="0.28"/>'
          f'<stop offset="1" stop-color="{T["win"]}" stop-opacity="0.10"/></linearGradient></defs>')
    if T["paper"] != "none":
        s += f'<rect width="{W}" height="{H}" fill="{T["paper"]}"/>'
    ovY = MT
    s += (f'<text x="{ML}" y="{ovY-5}" fill="{T["faint"]}" font-size="{px}">whole record · {_fmt_int(L)} bp</text>')
    s += (f'<rect x="{ML}" y="{ovY}" width="{plotW}" height="{ovH}" rx="3" fill="{T["lane"]}" '
          f'fill-opacity="{T["lane_op"]}" stroke="{T["frame"]}"/>')
    for t in model["tracks"]:
        for f in t["features"]:
            a = ox(f["start"]); w = max(ox(f["end"]) - ox(f["start"]), 1)
            s += f'<rect x="{a:.1f}" y="{ovY+2}" width="{w:.1f}" height="{ovH-4}" rx="1" fill="{f["color"]}" opacity="0.6"/>'
    s += (f'<rect x="{ox(s0):.1f}" y="{ovY-2:.1f}" width="{max(ox(s1)-ox(s0),2):.1f}" height="{ovH+4}" '
          f'rx="2" fill="url(#gvwin)" stroke="{T["win"]}" stroke-width="1.2"/>')
    ry = MT + ovH + 10 + rulerH - 7
    step = _tick_step(span, int(plotW // (_text_w("000.0k", px) + 10)))   # a tick every label-width plus 10 px, at most
    first = math.ceil(s0 / step) * step
    s += f'<line x1="{ML}" y1="{ry}" x2="{ML+plotW}" y2="{ry}" stroke="{T["frame"]}"/>'
    bp = first
    while bp <= s1 + 1:
        x = bx(bp)
        if not (x < ML - 1 or x > ML + plotW + 1):
            s += (f'<line x1="{x:.1f}" y1="{ry}" x2="{x:.1f}" y2="{ry-5}" stroke="{T["faint"]}"/>'
                  f'<text x="{x:.1f}" y="{ry-8:.1f}" fill="{T["faint"]}" font-size="{px}" text-anchor="middle">'
                  f'{gv_tick_label(bp)}</text>'
                  f'<line x1="{x:.1f}" y1="{ry+2:.1f}" x2="{x:.1f}" y2="{tracks_end+2}" stroke="{T["grid"]}"/>')
        bp += step
    for ti, pl in enumerate(plans):
        t = pl["t"]
        ty = track_ys[ti]; th = t.get("height", 20)
        stranded = bool(t.get("stranded"))
        name_x = ML - GV_NAME_GAP - (_strand_mark_w(name_px) if stranded else 0)
        s += (f'<text x="{name_x}" y="{ty+th/2+3:.1f}" fill="{T["faint"]}" font-size="{name_px}" '
              f'text-anchor="end">{esc(t["name"])}</text>')
        s += (f'<rect x="{ML}" y="{ty}" width="{plotW}" height="{th}" rx="3" fill="{T["track"]}" '
              f'fill-opacity="{T["track_op"]}"/>')
        if stranded:
            s += f'<line x1="{ML}" y1="{ty+th/2:.1f}" x2="{ML+plotW}" y2="{ty+th/2:.1f}" stroke="{T["grid"]}"/>'
            for mark, my in (("+", ty + th / 4), ("−", ty + 3 * th / 4)):    # which half-lane is which strand
                s += (f'<text x="{ML-4}" y="{my+px/3:.1f}" fill="{T["faint"]}" font-size="{px}" '
                      f'text-anchor="end">{mark}</text>')
        for f, a, b, inside, arrow in pl["items"]:
            stroke = _outline(f["color"], ground) if light else None
            tip = f'<title>{esc(f.get("tip",""))}</title>'
            if f.get("intron"):
                mid = (a + b) / 2
                s += (f'<path d="M {a:.1f} {ty+th/2:.1f} L {mid:.1f} {ty+3:.1f} L {b:.1f} {ty+th/2:.1f}" '
                      f'fill="none" stroke="{stroke or f["color"]}" stroke-width="1.3">{tip}</path>')
                y0, y1, hh = ty, ty + th, 0
            else:
                w = b - a
                yy = (ty + 2.5 if f.get("strand") == "+" else ty + th / 2 + 1.5) if stranded else ty + 2.5
                hh = max(th / 2 - 4 if stranded else th - 5, 3)
                if f.get("outline") == "sub":                  # an advisory pair: open, so it cannot pass for a call
                    paint = (f'fill="{f["color"]}" fill-opacity="0.15" stroke="{stroke or f["color"]}" '
                             f'stroke-width="1.2" stroke-dasharray="3 2"')
                elif f.get("emph"):                            # an ORF that carries a detected domain
                    paint = f'fill="{f["color"]}" stroke="{T["ink"]}" stroke-width="1.2"'
                elif stroke:
                    paint = f'fill="{f["color"]}" stroke="{stroke}" stroke-width="1"'
                else:
                    paint = f'fill="{f["color"]}"'
                s += (f'<rect class="gvglyph" x="{a:.1f}" y="{yy:.1f}" width="{w:.1f}" height="{hh:.1f}" '
                      f'rx="2.5" {paint}>{tip}</rect>')
                y0, y1 = yy, yy + hh
                ink = T["ink"] if f.get("outline") == "sub" else _label_ink(f["color"])
                tx = a + 4
                if arrow:                                      # head at the coding end; the label moves off a left head
                    my = yy + hh / 2
                    if f["strand"] == "+":
                        s += (f'<path d="M {b-8:.1f} {my-3.5:.1f} L {b-3:.1f} {my:.1f} L {b-8:.1f} {my+3.5:.1f} Z" '
                              f'fill="{ink}" pointer-events="none"/>')
                    else:
                        s += (f'<path d="M {a+8:.1f} {my-3.5:.1f} L {a+3:.1f} {my:.1f} L {a+8:.1f} {my+3.5:.1f} Z" '
                              f'fill="{ink}" pointer-events="none"/>')
                        tx = a + 12
                if inside:
                    s += (f'<text x="{tx:.1f}" y="{yy+hh-3:.1f}" fill="{ink}" font-size="{px}" font-weight="700" '
                          f'pointer-events="none">{esc(f["label"])}</text>')
                elif f.get("label"):                           # too narrow to letter: leader line + caption below
                    cx = (a + b) / 2
                    s += (f'<line x1="{cx:.1f}" y1="{y1:.1f}" x2="{cx:.1f}" y2="{ty+th+2:.1f}" '
                          f'stroke="{f["color"]}" stroke-width="0.6"/>')
            reg = {"x0": a, "y0": y0, "x1": b, "y1": y1, "start": f["start"], "end": f["end"],
                   "strand": f.get("strand", "+"), "tip": f.get("tip", ""), "kind": f.get("kind", ""),
                   "label": f.get("label", "intron" if f.get("intron") else "")}
            reg.update({k: f[k] for k in ("protein", "dna", "index", "derived", "domains") if k in f})
            regions.append(reg)
        for k, line in enumerate(pl["lines"]):                 # never drop a label (PBS·? hedge is often ~2 px wide)
            s += (f'<text x="{ML}" y="{ty+th+px+3+k*(px+2):.1f}" fill="{T["faint"]}" font-size="{px}">'
                  f'{esc(line)}</text>')
    if legend:
        ly0 = tracks_end + 4
        for (fill, text, style), (lx, row) in zip(entries, leg_pos):
            ly = ly0 + row * leg_row_h
            if style == "emph":
                paint = f'fill="{fill}" stroke="{T["ink"]}" stroke-width="1.2"'
            elif style == "open":
                paint = (f'fill="{fill}" fill-opacity="0.15" stroke="{_outline(fill, ground) if light else fill}" '
                         f'stroke-width="1.2" stroke-dasharray="3 2"')
            elif light:
                paint = f'fill="{fill}" stroke="{_outline(fill, ground)}" stroke-width="1"'
            else:
                paint = f'fill="{fill}"'
            s += (f'<rect x="{lx:.1f}" y="{ly:.1f}" width="9" height="9" rx="2" {paint}/>'
                  f'<text x="{lx+13:.1f}" y="{ly+8:.1f}" fill="{T["ink"]}" font-size="{px}">{esc(text)}</text>')
        s += (f'<text x="{lx0}" y="{ly0+(leg_rows-1)*leg_row_h+8:.1f}" fill="{T["faint"]}" font-size="{px}">'
              f'positions in bp from the sequence start (0-based)</text>')
    svg = s + "</svg>"
    return (svg, regions) if return_regions else svg


# ================= agarose gel =================
# palettes: theme.GELPAL (imported above)


def _call_mark(cx: float, cy: float, call: str, color: str, r: float = 2.4) -> str:
    """Redundant NON-COLOUR cue for a band's call, so on- vs off-target survives greyscale/colour-blind reading:
    on-target = solid triangle, off-target = cross, single-primer = open ring, priming site = filled square."""
    if call == "on":
        return f'<path d="M {cx-r:.1f} {cy-r:.1f} L {cx+r:.1f} {cy:.1f} L {cx-r:.1f} {cy+r:.1f} Z" fill="{color}"/>'
    if call == "off":
        return (f'<path d="M {cx-r:.1f} {cy-r:.1f} L {cx+r:.1f} {cy+r:.1f} M {cx-r:.1f} {cy+r:.1f} '
                f'L {cx+r:.1f} {cy-r:.1f}" stroke="{color}" stroke-width="1.1" stroke-linecap="round" fill="none"/>')
    if call == "single":
        return f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{max(r-0.6,1.0):.1f}" fill="none" stroke="{color}" stroke-width="1"/>'
    return f'<rect x="{cx-r+0.5:.1f}" y="{cy-r+0.5:.1f}" width="{2*r-1:.1f}" height="{2*r-1:.1f}" fill="{color}"/>'


def _band_opacity(total_mm: int) -> float:
    """Priming-efficiency proxy: a perfect match reads bright; each mismatch dims the band.
    All reported mismatches are 5'-proximal (the strict-3' rule forbids 3'-end mismatches)."""
    return round(max(0.4, 1.0 - 0.22 * max(0, total_mm)), 3)


def _call_counts(g):
    """Disjoint, exhaustive split of one band's products -> (n_on, n_off, n_single), n_on+n_off+n_single == len(g).
    A self-priming product is an artefact, never the intended amplicon, so it is counted as single-primer even if
    it falls inside the target window; off-target is the remainder, so no derived count can go negative."""
    n_single = sum(1 for a in g if a.get("single_primer"))
    n_on = sum(1 for a in g if a.get("on_target") and not a.get("single_primer"))
    return n_on, len(g) - n_on - n_single, n_single


def _is_design_bound(a) -> bool:
    """A product of a design-bound local PCR carries the engine's verdict on the DESIGNED endpoints
    (`matches_intended_product`, or the contract's `intended`). A genome-scan product does not."""
    return "matches_intended_product" in a or "intended" in a


def _is_intended(a) -> bool:
    return bool(a["intended"]) if "intended" in a else a.get("matches_intended_product") is True


def _call_split(g, has_locus: bool = True) -> dict:
    """Disjoint, exhaustive split of one band's products: on-target, additional (a design-bound product that is not
    the designed one — the table's 'additional product', not an established off-target), off-target (a genome-scan
    product away from the intended locus, which the app itself defines), single-primer, and site (no on/off call is
    possible: no design locus, a target the run did not assess, or the designed product that misses the Primer3
    target). The counts sum to len(g); an unassessed product is never counted as off-target."""
    n = {"on": 0, "additional": 0, "off": 0, "single": 0, "site": 0}
    for a in g:
        if a.get("single_primer"):
            n["single"] += 1                       # a self-priming artefact is never the intended amplicon
        elif not has_locus or a.get("on_target") is None:
            n["site"] += 1
        elif a.get("on_target"):
            n["on"] += 1
        elif _is_intended(a):
            n["site"] += 1                         # the designed product, but it does not cover the target
        else:
            n["additional" if _is_design_bound(a) else "off"] += 1
    return n


def _call_text(n: dict, g, has_locus: bool = True) -> str:
    """The band tooltip's call list: '1 on-target, 2 additional products, 1 single-primer'."""
    parts = []
    if n["on"]:
        parts.append(f"{n['on']} on-target")
    if n["additional"]:
        parts.append(f"{n['additional']} additional product" + ("s" if n["additional"] != 1 else ""))
    if n["off"]:
        parts.append(f"{n['off']} off-target")
    if n["site"]:
        if not has_locus:
            des = sum(1 for a in g if not a.get("single_primer") and _is_intended(a))
            if des:
                parts.append(f"{des} designed product" + ("s" if des != 1 else "") + " (target not assessed)")
            if n["site"] - des:
                parts.append(f"{n['site'] - des} priming site" + ("s" if n["site"] - des != 1 else ""))
        else:
            intended = sum(1 for a in g if not a.get("single_primer") and _is_intended(a) and not a.get("on_target"))
            parts.append(f"{n['site']} product" + ("s" if n["site"] != 1 else "") + " with no on/off call"
                         + (" (designed product, target not covered)" if intended == n["site"] else ""))
    if n["single"]:
        parts.append(f"{n['single']} single-primer")
    return ", ".join(parts)


def _physical_product_length(amplicon):
    return amplicon.get("product_length", amplicon["length"])


def _public_amplicon_source(amplicon):
    """The record label a tooltip may show. A validated-fetch or generated identity (a verified chromosome
    accession, a generated input id) stays readable, as in the table; anything else is public only through its full
    digest binding."""
    trust = amplicon.get("sourceIdentitySource")
    if trust in {"validated_fetch", "generated"}:
        return public_record_label(amplicon.get("source"), trust=trust)
    return bound_public_record_label(amplicon.get("source"), amplicon.get("sourceIdSha256"))


def _unique_physical_products(amplicons):
    """Collapse repeated mechanisms for one physical product without merging distinct genomic intervals."""
    products = {}
    for index, a in enumerate(amplicons or []):
        size = _physical_product_length(a)
        if a.get("product_key"):
            key = a["product_key"]
        elif a.get("start") is not None and a.get("end") is not None:
            key = (a.get("source"), a["start"], a["end"], size)
        else:
            key = ("unbound-row", index)
        if (key not in products or
                (a.get("matches_intended_product") is True and
                 products[key].get("matches_intended_product") is not True) or
                (products[key].get("single_primer") and not a.get("single_primer"))):
            products[key] = a
    return list(products.values())


BAND_MERGE_TOL = 0.03     # products within 3% of a band's smallest product are drawn as one band


def _band_groups(amplicons):
    """The lane's physical products grouped into the bands a gel would show: sorted by size, a new band starting
    when a product is more than BAND_MERGE_TOL larger than the smallest in the current one. On the 240 px log axis
    3% is about 2.5 px, the height of a band, so products that close are one smear on screen whatever colour is
    drawn last; equal lengths (a real gel cannot resolve them) always merge. A display heuristic, not a statement
    about resolving power: the table lists every product and the tooltip gives the size range and count."""
    groups = []
    for a in sorted(_unique_physical_products(amplicons), key=_physical_product_length):
        size = _physical_product_length(a)
        if groups and size <= _physical_product_length(groups[-1][0]) * (1 + BAND_MERGE_TOL):
            groups[-1].append(a)
        else:
            groups.append([a])
    return groups


def _group_size(g) -> int:
    """Where a band is drawn: the designed product's size when it is in the band, else the middle of the band."""
    on = [a for a in g if a.get("on_target") and not a.get("single_primer")]
    pick = on or g
    return sorted(_physical_product_length(a) for a in pick)[(len(pick) - 1) // 2]


def _size_text(g) -> str:
    lo, hi = min(_physical_product_length(a) for a in g), max(_physical_product_length(a) for a in g)
    return f"{lo} bp" if lo == hi else f"{lo}–{hi} bp"


def _lane_bands(amplicons, P, has_locus=True):
    """Collapse a lane's amplicons into one band per product size (a real gel cannot resolve equal lengths, and
    near-equal ones, within BAND_MERGE_TOL, share one drawn band — see _band_groups). Band intensity follows the
    strongest (fewest-mismatch) product in the band.

    has_locus=True (a designed on-target exists): worst-case colour wins — a band that also carries an additional
    or off-target product (or a single-primer artefact) is NOT a clean on-target, so it reads off-target (a
    specificity warning), never a reassuring on-target band. has_locus=False (whole-genome scan of a bare consensus
    pair): there is no on/off target, so every product is a NEUTRAL 'genomic priming site' and must not read as the
    off-target warning colour. A product with no on/off call (target not assessed) is neutral too. Every product is
    still enumerated in the table below the gel."""
    bands = []
    for g in _band_groups(amplicons):
        n = _call_split(g, has_locus)                     # disjoint buckets: no count can go negative
        if not has_locus:                                 # neutral priming sites (single-primer artefacts still distinct)
            call = "single" if (n["single"] and not n["site"]) else "site"
        elif n["additional"] or n["off"]:                 # any extra product co-migrating here -> the whole band is flagged
            call = "off"
        elif n["single"]:
            call = "single"
        elif n["on"]:                                     # purely on-target in this band
            call = "on"
        else:
            call = "site"
        min_mm = min((a.get("fwd_mm", 0) + a.get("rev_mm", 0)) for a in g)
        raw_source = g[0].get("source", "")
        src = esc(_public_amplicon_source(g[0])) if raw_source else ""
        color = P[call]                                   # colour derives from the call, so the shape cue can never disagree
        bands.append({"size": _group_size(g), "size_text": _size_text(g), "color": color, "call": call,
                      "opacity": _band_opacity(min_mm), "on": bool(n["on"]) and has_locus, "single": bool(n["single"]),
                      "count": len(g), "merged": len({_physical_product_length(a) for a in g}) > 1,
                      "t": _call_text(n, g, has_locus) + (f" · {src}" if src else "")})
    return bands


LANES_PER_ROW = 10        # sample lanes per gel row; more than this wraps onto stacked rows


GEL_PX = 9                # every gel text is at least this size, in every palette
# size marks of the ladder lane, as a 1 kb Plus-style set; marks above 1.5 kb are drawn only when a product needs them
_LADDER_MARKS = (6000, 5000, 4000, 3000, 2000, 1500, 1000, 700, 500, 400, 300, 200, 100, 50)


def _lane_status(lane: dict) -> str:
    return str(lane.get("status") or "complete")


def _lane_footer(lane: dict):
    """The caption lines under a lane, [(text, role)], role in {'warn', 'ink', 'bold'}. Never states a specificity
    claim the lane cannot support: a failed search says 'failed' (it is not an empty lane), a search that hit its
    cap says 'capped' (a lower bound — never 'no on-target'), and a lane with no designed target prints no
    on-target verdict at all (its tooltip says why)."""
    status = _lane_status(lane)
    amps = lane.get("amplicons") or []
    if status == "failed":
        return [("failed", "bold")]
    if status == "not_run":
        return [("not run", "ink")]
    bands = _lane_bands(amps, GELPAL["dark"], lane.get("has_locus", True))
    if status == "partial_capped":
        return [("capped", "bold")]
    if not bands:
        return [("—", "ink")]
    if lane.get("has_locus", True) and not lane.get("advisory", False) and not any(b["on"] for b in bands):
        return [("no on-target", "warn")]                    # bands present but none intended
    return []


def _hatch(x, y, w, h, step=5.0):
    """A path of 45-degree hatch lines clipped analytically to the rectangle (QtSvg has no reliable <pattern>)."""
    segs = []
    c = -h
    while c < w:
        x0, y0 = max(c, 0.0), max(-c, 0.0)                  # line  y = x - c, entering at the left or top edge
        x1 = min(w, h + c); y1 = x1 - c
        if x1 > x0:
            segs.append(f"M {x + x0:.1f} {y + y0:.1f} L {x + x1:.1f} {y + y1:.1f}")
        c += step
    return " ".join(segs)


def _gel_legend_items(lanes):
    """[(kind, text)] for the states this gel actually draws; kind picks the swatch. Lanes with a design locus and
    lanes without one each contribute their own entries, so a mixed gel explains every band colour it uses."""
    neutral = [l for l in lanes if not l.get("has_locus", True)]
    scored_lanes = [l for l in lanes if l.get("has_locus", True)]
    amps = [a for l in lanes for a in (l.get("amplicons") or [])]
    items = []
    plain = [a for l in neutral for a in (l.get("amplicons") or []) if not a.get("single_primer")]
    if neutral and (plain or not scored_lanes):
        if any(_is_intended(a) for a in plain):
            items.append(("site", "designed product (target not assessed)"))
        if not plain or any(not _is_intended(a) for a in plain):
            items.append(("site", "priming site"))
    if scored_lanes or not neutral:
        scored = [a for l in scored_lanes for a in (l.get("amplicons") or [])]
        n = _call_split(scored, True)
        words = [w for w, k in (("additional product", "additional"), ("off-target", "off")) if n[k]]
        items += [("on", "on-target"), ("off", " / ".join(words) or "off-target")]
        if any(a.get("on_target") is None and not a.get("single_primer") for a in scored):
            items.append(("site", "no on/off call"))
    if any(a.get("single_primer") for a in amps):
        items.append(("single", "single-primer"))
    if any(_lane_status(l) == "failed" for l in lanes):
        items.append(("failed", "failed (not assessed)"))
    if any(_lane_status(l) == "partial_capped" for l in lanes):
        items.append(("capped", "capped (lower bound)"))
    return items


def _legend_item_w(text, px=GEL_PX):
    return 6 + 5 + _text_w(text, px) + 14


def _gel_geometry(data: dict):
    """Shared layout math for svg_gel and gel_regions (must stay identical so hit-boxes line up).
    Sample lanes wrap into stacked rows of <=LANES_PER_ROW; each row carries its own ladder + bp scale."""
    lanes = data.get("lanes") or [{"label": data.get("laneLabel", "PCR"), "amplicons": data.get("amplicons", [])}]
    sizes = [_physical_product_length(a) for l in lanes
             for a in _unique_physical_products(l.get("amplicons") or [])]
    smallest = min(sizes) if sizes else 90
    minbp = max(25, min(90, smallest - 10))
    # the top of the scale is the next ladder mark at or above the biggest product (1600 as before when nothing is
    # larger), so a 3.8 kb product sits between marked sizes instead of floating above the last one
    big = max(sizes) if sizes else 0
    maxbp = max(1600, big)
    if big > 1600:
        maxbp = next((m for m in sorted(_LADDER_MARKS) if m >= big), big)
    LADDER = [m for m in _LADDER_MARKS if minbp <= m <= maxbp]
    laneW, gap, x0 = 40, 28, 62                             # pitch 68: a 9 px 'no on-target' caption (65 px) clears its neighbours'
    foot = max([len(_lane_footer(l)) for l in lanes] + [1])
    TOP_MARGIN, ROW_TOP_PAD, BODY, ROW_GAP = 8, 30, 240, 16
    ROW_BOT_PAD = max(20, 12 + foot * (GEL_PX + 2) + 6)     # room for the tallest caption of any lane
    ROW_PITCH = ROW_TOP_PAD + BODY + ROW_BOT_PAD + ROW_GAP
    rows = [list(range(i, min(i + LANES_PER_ROW, len(lanes)))) for i in range(0, len(lanes), LANES_PER_ROW)] or [[]]
    widest = min(LANES_PER_ROW, max(1, len(lanes)))
    ladder_text = f"L = size ladder, {min(LADDER) if LADDER else minbp}–{max(LADDER) if LADDER else maxbp} bp"
    items = _gel_legend_items(lanes) + [("ladder", ladder_text)]
    # the gel is as wide as its lanes need, or as the widest single legend entry; the legend wraps to more rows
    # rather than stretching the figure (a wide empty band beside two lanes)
    W = math.ceil(max(x0 + (1 + widest) * (laneW + gap) + 12, 300,
                      x0 - 3 + max(_legend_item_w(t) for _, t in items) + 12))
    legend_pos, xs, row = [], x0 - 3, 0
    for kind, text in items:
        wi = _legend_item_w(text)
        if xs > x0 - 3 and xs + wi > W - 12:
            xs, row = x0 - 3, row + 1
        legend_pos.append((kind, text, xs, row))
        xs += wi
    leg_rows = row + 1
    last_bot = TOP_MARGIN + (len(rows) - 1) * ROW_PITCH + ROW_TOP_PAD + BODY + ROW_BOT_PAD
    H = last_bot + 22 + (leg_rows - 1) * (GEL_PX + 4)       # legend strip below the last row

    def row_top(r):                                        # band-area top of row r
        return TOP_MARGIN + r * ROW_PITCH + ROW_TOP_PAD

    def y(bp, r):
        t = row_top(r)
        return t + (math.log(maxbp) - math.log(max(bp, minbp))) / (math.log(maxbp) - math.log(minbp)) * BODY

    def laneX(col):
        return x0 + col * (laneW + gap)

    return {"lanes": lanes, "rows": rows, "minbp": minbp, "maxbp": maxbp, "LADDER": LADDER,
            "laneW": laneW, "gap": gap, "x0": x0, "BODY": BODY, "W": W, "H": H, "foot": foot,
            "ROW_BOT_PAD": ROW_BOT_PAD, "legend": items, "legend_pos": legend_pos, "leg_rows": leg_rows,
            "ladder_text": ladder_text,
            "y": y, "laneX": laneX, "row_top": row_top}


def gel_regions(data: dict):
    """Per-band hit-boxes (SVG coords) + the amplicon each band represents, for hover/right-click.
    Follows the wrapped multi-row layout; the ladder column (0) is skipped, sample lanes start at 1."""
    G = _gel_geometry(data)
    lanes, rows, laneW, y, laneX = G["lanes"], G["rows"], G["laneW"], G["y"], G["laneX"]
    regions = []
    for r, idxs in enumerate(rows):
        for j, li in enumerate(idxs):
            l = lanes[li]
            lx = laneX(j + 1)
            has_locus = l.get("has_locus", True)
            for g in _band_groups(l.get("amplicons") or []):
                yy = y(_group_size(g), r)
                # representative for the right-click menu: the intended (on-target) product if any, else the strongest
                rep = min(g, key=lambda a: (not a.get("on_target"), a.get("fwd_mm", 0) + a.get("rev_mm", 0)))
                call = _call_text(_call_split(g, has_locus), g, has_locus)   # same split as the band tooltip
                raw_source = rep.get("source", "")
                public_source = _public_amplicon_source(rep) if raw_source else ""
                tip = f'{_size_text(g)} · {call}' + (f' · {public_source}' if public_source else "")
                if len({_physical_product_length(a) for a in g}) > 1:
                    tip += f" · {len(g)} products drawn as one band (within {BAND_MERGE_TOL:.0%})"
                regions.append({"x0": lx + 3, "y0": yy - 4, "x1": lx + laneW - 3, "y1": yy + 4,
                                "tip": tip, "amplicon": rep, "pair": l.get("label", ""),
                                "has_locus": has_locus, "status": _lane_status(l)})
    return regions


def svg_gel(data: dict, bg: str) -> str:
    """data: {lanes:[{label, amplicons:[{length,on_target,source}]}]} or legacy {amplicons}.
    Lanes beyond LANES_PER_ROW wrap onto stacked gel rows, each with its own ladder + bp scale."""
    G = _gel_geometry(data)
    lanes, rows, LADDER = G["lanes"], G["rows"], G["LADDER"]
    laneW, gap, x0, BODY = G["laneW"], G["gap"], G["x0"], G["BODY"]
    W, H, y, laneX, row_top, foot = G["W"], G["H"], G["y"], G["laneX"], G["row_top"], G["foot"]
    P = GELPAL.get(bg, GELPAL["transparent"])
    px = GEL_PX
    s = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" '
         f'font-family="{FIGFONT}">')
    s += (f'<defs><filter id="glow" x="-40%" y="-140%" width="180%" height="380%">'
          f'<feGaussianBlur stdDeviation="{P["glow"]}" result="b"/><feMerge>'
          f'<feMergeNode in="b"/><feMergeNode in="SourceGraphic"/></feMerge></filter></defs>')
    if P["paper"] != "none":
        s += f'<rect width="{W}" height="{H}" fill="{P["paper"]}"/>'

    _LANE_TIP = {"failed": "This pair's search failed, so nothing was assessed — that is not the same as no product.",
                 "partial_capped": "The search hit its cap: the bands drawn are a lower bound, and no specificity "
                                   "conclusion is supported.",
                 "not_run": "This pair was not run."}

    def draw_lane(col, r, lane, bands, is_ladder):
        nonlocal s
        lx = laneX(col)
        rtop = row_top(r); rbot = rtop + BODY
        label = "L" if is_ladder else lane["label"]
        tip = ""
        if not is_ladder:
            status = _lane_status(lane)
            designed = any(_is_intended(a) for a in lane.get("amplicons") or [])
            msg = _LANE_TIP.get(status) or ("" if lane.get("has_locus", True) else
                                            "The design set no target for this lane, so on-target and off-target are "
                                            "not assessed." if designed else
                                            "No designed target was set for this lane, so its products are genomic "
                                            "priming sites, not on-target or off-target calls.")
            tip = f"<title>{esc(msg)}</title>" if msg else ""
        s += f'<rect x="{lx+4:.1f}" y="{rtop-13}" width="{laneW-8}" height="4" rx="1" fill="{P["well"]}"/>'
        s += (f'<text x="{lx+laneW/2:.1f}" y="{rtop-18}" fill="{P["ink"]}" font-size="{px}" '
              f'text-anchor="middle">{esc(label)}{tip}</text>')
        if not is_ladder and _lane_status(lane) == "failed":   # a failed search is hatched, never an empty lane
            s += (f'<rect x="{lx+3:.1f}" y="{rtop}" width="{laneW-6}" height="{BODY}" fill="none" '
                  f'stroke="{P["ink"]}" stroke-width="0.6" stroke-dasharray="3 2"/>'
                  f'<path d="{_hatch(lx + 3, rtop, laneW - 6, BODY)}" stroke="{P["ink"]}" stroke-width="0.6" '
                  f'opacity="0.45" fill="none"/>')
        for bd in (bands or []):
            yy = y(bd["size"], r); h = 1.6 if is_ladder else P["band"]
            filt = "" if is_ladder else ' filter="url(#glow)"'
            op = "" if is_ladder else f' fill-opacity="{bd.get("opacity", 1.0)}"'
            dash = ' stroke="#ffffff" stroke-width="0.5" stroke-dasharray="2 1.5"' if bd.get("single") else ""
            n = bd.get("count", 1)
            title = (bd.get("size_text") or f'{bd["size"]} bp') + (f' · {bd["t"]}' if bd.get("t") else "") \
                + (f' (×{n})' if n > 1 else "")
            s += (f'<rect x="{lx+3:.1f}" y="{yy-h/2:.1f}" width="{laneW-6}" height="{h}" rx="1" '
                  f'fill="{bd["color"]}"{op}{dash}{filt}><title>{title}</title></rect>')
            if not is_ladder and bd.get("call"):          # shape marker: the call read without colour
                s += _call_mark(lx - 1.0, yy, bd["call"], bd["color"])
        if not is_ladder:
            for k, (text, role) in enumerate(_lane_footer(lane)):
                cx, cy = lx + laneW / 2, rbot + 12 + k * (px + 2)
                fill = P["off"] if role == "warn" else P["ink"]
                weight = ' font-weight="700"' if role == "bold" else ""
                if text == "capped":                       # a boxed badge: this lane is a lower bound
                    bw = _text_w(text, px) + 6
                    s += (f'<rect x="{cx-bw/2:.1f}" y="{cy-px+1:.1f}" width="{bw:.1f}" height="{px+2}" rx="2" '
                          f'fill="none" stroke="{P["ink"]}" stroke-width="0.8"/>')
                s += (f'<text x="{cx:.1f}" y="{cy:.1f}" fill="{fill}" font-size="{px}"{weight} '
                      f'text-anchor="middle">{esc(text)}</text>')

    box_h = 16 + BODY + 12 + (foot - 1) * (px + 2) + 4
    for r, idxs in enumerate(rows):
        rtop = row_top(r)
        ncols = 1 + len(idxs)
        s += (f'<rect x="{x0-7:.1f}" y="{rtop-16}" width="{ncols*(laneW+gap)+2:.1f}" height="{box_h:.1f}" '
              f'rx="2" fill="{P["gel"]}" stroke="{P["stroke"]}"/>')
        s += f'<text x="{x0-14}" y="{rtop-19}" fill="{P["ink"]}" font-size="{px}" text-anchor="end">bp</text>'
        for m in LADDER:
            s += f'<text x="{x0-14}" y="{y(m, r)+3:.1f}" fill="{P["ink"]}" font-size="{px}" text-anchor="end">{m}</text>'
        draw_lane(0, r, {}, [{"size": m, "color": P["ladder"]} for m in LADDER], True)
        for j, li in enumerate(idxs):
            l = lanes[li]
            draw_lane(j + 1, r, l, _lane_bands(l.get("amplicons") or [], P, l.get("has_locus", True)), False)

    # legend: only the states this gel draws, laid out from measured widths (the gel is widened to hold it)
    for kind, text, xs, row in G["legend_pos"]:
        ly = H - 8 - (G["leg_rows"] - 1 - row) * (px + 4)
        cx = xs + 3
        if kind in ("on", "off", "single", "site"):
            s += _call_mark(cx, ly, kind, P[kind], 3)
        elif kind == "failed":
            s += (f'<rect x="{xs:.1f}" y="{ly-4:.1f}" width="7" height="8" fill="none" stroke="{P["ink"]}" '
                  f'stroke-width="0.8"/><path d="{_hatch(xs, ly - 4, 7, 8, 3.0)}" stroke="{P["ink"]}" '
                  f'stroke-width="0.6" fill="none"/>')
        elif kind == "capped":
            s += (f'<rect x="{xs:.1f}" y="{ly-4:.1f}" width="8" height="8" rx="2" fill="none" '
                  f'stroke="{P["ink"]}" stroke-width="0.8"/>')
        tx = xs if kind == "ladder" else xs + 11
        s += f'<text x="{tx:.1f}" y="{ly+3}" fill="{P["ink"]}" font-size="{px}">{esc(text)}</text>'
    return s + "</svg>"


# ============================ self-similarity: dot plot + heat map ============================
# One locus compared with itself. Rendered through the same hand-written SVG path as the genome viewer
# and the gel — no numpy/PIL/matplotlib, which the build guard excludes and which would undo the
# v3.0.0 installer reduction. Forward matches mark direct repeats (LTRs); reverse-complement matches
# mark inverted repeats (TIRs), so the two are drawn as distinguishable layers rather than merged.
_DOT_THEME = {
    "dark":   {"bg": "#0F1519", "grid": "#243038", "ink": "#E6EDF1", "faint": "#7C8B95"},
    "white":  {"bg": "#FFFFFF", "grid": "#DDE4E9", "ink": "#141B21", "faint": "#5B6B76"},
    "export": {"bg": "#FFFFFF", "grid": "#D4DBE0", "ink": "#222222", "faint": "#555555"},
}
# Forward (direct) and reverse (inverted) get different hues AND different marks, so the two signals
# stay separable in greyscale and under colour-vision deficiency. Blue and vermillion are the pair the
# Okabe-Ito set keeps furthest apart under deuteranopia and protanopia.
_DOT_FWD, _DOT_REV = OKABE_ITO["blue"], OKABE_ITO["vermillion"]

# The exported figure must carry its own limit — the on-screen caveat is a separate widget that does not
# travel into an SVG/PNG/PDF a reader drops into a paper. One concise line, added only on export.
_DOT_FOOTER = ("Exact k-mer matching, not local alignment — an absent or faint diagonal is not evidence a "
               "repeat is absent; a visible block is.")


def _dot_run_caveat(m):
    """Run-SPECIFIC caveats (masking / truncation) that the on-screen scope note shows. An exported figure
    must carry them too — omitting them would let a masked or truncated picture read as complete."""
    parts = []
    if m.get("masked_kmers"):
        parts.append("%d k-mer(s) over %d× were masked (%d positions) as simple/satellite repeat and are "
                     "not shown." % (m["masked_kmers"], m.get("max_occ", 0), m.get("masked_positions", 0)))
    if m.get("truncated"):
        parts.append("Comparison truncated at the work limit — the picture is incomplete.")
    return " ".join(parts)


DOT_ML, DOT_MT, DOT_MR, DOT_MB = 58.0, 46.0, 16.0, 66.0   # plot-box margins; the canvas hit-test reads the first three


def _dot_footer_lines(m, W, px):
    """[(text, bold)] — the generic method caveat plus any run-specific masking/truncation caveat, wrapped to the
    plot width. An export carries them; wrapping (not shrinking) keeps the 10 px export floor."""
    width = max(W - DOT_ML - DOT_MR, 80.0)
    out = [(ln, False) for ln in _wrap_mono(_DOT_FOOTER.split(" "), width, px, sep=" ")]
    run = _dot_run_caveat(m)
    if run:
        out += [(ln, True) for ln in _wrap_mono(run.split(" "), width, px, sep=" ")]
    return out


def _dot_footer_svg(ML, y, T, lines, px):
    return "".join('<text x="%.1f" y="%.1f" fill="%s" font-size="%s"%s>%s</text>' % (
        ML, y + i * (px + 2), T["faint"], px, ' font-weight="700"' if bold else "", esc(ln))
        for i, (ln, bold) in enumerate(lines))


def _dot_theme(theme: str, for_export: bool) -> dict:
    return dict(_DOT_THEME["export" if for_export else ("white" if theme == "white" else "dark")])


_HEAT_FLOOR = "#F5F7F9"                                      # the empty end of the ramp on a light ground


def _heat_ramp(t: float, base: str, floor: str = _HEAT_FLOOR) -> str:
    """Sequential ramp from `floor` — the figure's own background colour, so a sparse bin is as quiet as an empty one
    — to the signal hue, monotone in lightness so the scale still reads in greyscale (an unordered rainbow would not)."""
    t = 0.0 if t < 0 else (1.0 if t > 1 else t)
    r0, g0, b0 = _rgb(floor)
    r1, g1, b1 = int(base[1:3], 16), int(base[3:5], 16), int(base[5:7], 16)
    return "#%02X%02X%02X" % (int(r0 + (r1 - r0) * t), int(g0 + (g1 - g0) * t), int(b0 + (b1 - b0) * t))


def _dot_frame(m, W, T, title, sub, px=FIG_PX_SCREEN, footer_lines=0):
    """Shared chrome: a square plot box with a bp ruler on both axes, title and scope line. The bottom margin is the
    on-screen one unless a footer of `footer_lines` lines (an export) needs more."""
    ML, MT, MR = DOT_ML, DOT_MT, DOT_MR
    MB = DOT_MB if not footer_lines else max(DOT_MB, 44 + footer_lines * (px + 2) + 10)
    plot = max(W - ML - MR, 80.0)
    H = MT + plot + MB
    n = max(m["length"], 1)
    s = ('<svg xmlns="http://www.w3.org/2000/svg" width="%.0f" height="%.0f" viewBox="0 0 %.0f %.0f" '
         'font-family="%s">' % (W, H, W, H, FIGFONT))
    s += '<rect width="%.0f" height="%.0f" fill="%s"/>' % (W, H, T["bg"])
    s += ('<text x="%.1f" y="20" fill="%s" font-size="12" font-weight="700">%s</text>'
          '<text x="%.1f" y="34" fill="%s" font-size="%s">%s</text>'
          % (ML, T["ink"], esc(title), ML, T["faint"], px, esc(sub)))
    s += ('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="none" stroke="%s" stroke-width="1"/>'
          % (ML, MT, plot, plot, T["grid"]))
    step = _tick_step(n, int(plot // (_text_w("000.0k", px) + 10)))   # the genome ruler's rule: a label plus 10 px each
    t = 0
    while t <= n:
        x = ML + plot * (t / n)
        y = MT + plot * (t / n)
        s += ('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="0.5" opacity="0.55"/>'
              '<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="0.5" opacity="0.55"/>'
              '<text x="%.1f" y="%.1f" fill="%s" font-size="%s" text-anchor="middle">%s</text>'
              '<text x="%.1f" y="%.1f" fill="%s" font-size="%s" text-anchor="end">%s</text>'
              % (x, MT, x, MT + plot, T["grid"], ML, y, ML + plot, y, T["grid"],
                 x, MT + plot + 13, T["faint"], px, gv_tick_label(t),
                 ML - 6, y + 3, T["faint"], px, gv_tick_label(t)))
        t += step
    return s, ML, MT, plot, H


def _dot_guides(m, ML, MT, plot, T, guides):
    """Translucent bands marking TEagle's own detected features, so the picture is read against the
    calls it exists to corroborate or contradict. Never obscures the data."""
    if not guides:
        return ""
    n = max(m["length"], 1)
    s = ""
    for g in guides:
        a, z = g.get("start", 0), g.get("end", 0)
        if z <= a:
            continue
        x0, x1 = ML + plot * (a / n), ML + plot * (z / n)
        col = g.get("color", T["faint"])
        s += ('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="%s" opacity="0.10"/>'
              '<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="%s" opacity="0.10"/>'
              % (x0, MT, max(x1 - x0, 1), plot, col,
                 ML, MT + plot * (a / n), plot, max(plot * ((z - a) / n), 1), col))
    return s


def _dot_legend(x, y, T, heat=None, fwd=_DOT_FWD, rev=_DOT_REV, px=FIG_PX_SCREEN, floor=_HEAT_FLOOR):
    if heat is None:
        first = "direct repeat (forward)"
        x2 = x + 13 + _text_w(first, px) + 16                  # the second entry starts where the first one ends
        return ('<circle cx="%.1f" cy="%.1f" r="3" fill="%s"/>'
                '<text x="%.1f" y="%.1f" fill="%s" font-size="%s">%s</text>'
                '<rect x="%.1f" y="%.1f" width="6" height="6" fill="%s"/>'
                '<text x="%.1f" y="%.1f" fill="%s" font-size="%s">inverted repeat (reverse complement)</text>'
                % (x + 4, y - 3, fwd, x + 13, y, T["ink"], px, first,
                   x2, y - 6, rev, x2 + 11, y, T["ink"], px))
    base, peak = heat
    s = ""
    for i in range(24):
        s += '<rect x="%.1f" y="%.1f" width="7" height="8" fill="%s"/>' % (
            x + i * 7, y - 8, _heat_ramp(i / 23.0, base, floor))
    s += ('<text x="%.1f" y="%.1f" fill="%s" font-size="%s">0</text>'
          '<text x="%.1f" y="%.1f" fill="%s" font-size="%s">%d matches / bin</text>'
          % (x, y + 11, T["faint"], px, x + 176, y - 1, T["faint"], px, peak))
    return s


def svg_dotplot(m: dict, W: float = 620, theme: str = "dark", for_export: bool = False,
                guides: list = None, title: str = "Self-similarity dot plot",
                read_threshold: int = 1, fwd: str = None, rev: str = None) -> str:
    """Binary dot matrix: a mark wherever the locus matches itself."""
    T = _dot_theme(theme, for_export)
    px = FIG_PX_EXPORT if for_export else FIG_PX_SCREEN
    fwd, rev = fwd or _DOT_FWD, rev or _DOT_REV       # caller may override the mark hues; bg stays white on export
    sub = "exact %d-mer matches · %s bp" % (m["k"], _fmt_int(m["length"]))
    if read_threshold and read_threshold > 1:
        sub += " · cells below %d matches shown faint (chance level)" % read_threshold
    foot = _dot_footer_lines(m, W, px) if for_export else []
    s, ML, MT, plot, H = _dot_frame(m, W, T, title, sub, px, len(foot))
    b = m["bins"]
    cell = plot / b
    r = max(cell * 0.42, 0.5)
    # Cells below the chance-derived read threshold are drawn FAINT, never removed: at a small word size
    # most single-count cells are noise, but hiding them would leave a reader unable to judge how noisy
    # the panel is. Emphasis, not filtering.
    thr = read_threshold or 1
    for name, mat, col in (("fwd", m["forward"], fwd), ("rev", m["reverse"], rev)):
        for i, row in enumerate(mat):
            y = MT + (i + 0.5) * cell
            for j, c in enumerate(row):
                if not c:
                    continue
                x = ML + (j + 0.5) * cell
                op = "" if c >= thr else ' opacity="0.22"'
                if name == "fwd":
                    s += '<circle cx="%.2f" cy="%.2f" r="%.2f" fill="%s"%s/>' % (x, y, r, col, op)
                else:                       # a different MARK, not merely a different hue
                    s += '<rect x="%.2f" y="%.2f" width="%.2f" height="%.2f" fill="%s"%s/>' % (
                        x - r, y - r, 2 * r, 2 * r, col, op)
    s += _dot_guides(m, ML, MT, plot, T, guides)
    s += _dot_legend(ML, MT + plot + 28, T, heat=None, fwd=fwd, rev=rev, px=px)
    if foot:
        s += _dot_footer_svg(ML, MT + plot + 44, T, foot, px)
    return s + "</svg>"


def svg_dotheat(m: dict, W: float = 620, theme: str = "dark", for_export: bool = False,
                guides: list = None, which: str = "forward",
                title: str = "Self-similarity heat map", fwd: str = None, rev: str = None) -> str:
    """Binned match DENSITY. The dot plot answers "is there a repeat"; the heat map answers "how much of
    the locus takes part", which a binary mark cannot show once bins saturate."""
    T = _dot_theme(theme, for_export)
    px = FIG_PX_EXPORT if for_export else FIG_PX_SCREEN
    fwd, rev = fwd or _DOT_FWD, rev or _DOT_REV
    peak = m.get("%s_max" % which) or 1
    kind = "direct (forward)" if which == "forward" else "inverted (reverse-complement)"
    sub = "%s match density · exact %d-mer · peak %d per bin · %s bp" % (
        kind, m["k"], peak, _fmt_int(m["length"]))
    foot = _dot_footer_lines(m, W, px) if for_export else []
    s, ML, MT, plot, H = _dot_frame(m, W, T, title, sub, px, len(foot))
    base = fwd if which == "forward" else rev
    # the ramp starts at the figure's own ground (an empty bin is not drawn), so on the dark theme a sparse bin is
    # dark and the densest is the brightest — a light floor there made the noisiest bins the most salient
    floor = T["bg"] if T["bg"].lstrip("#").upper() != "FFFFFF" else _HEAT_FLOOR
    mat = m[which]
    cell = plot / m["bins"]
    for i, row in enumerate(mat):
        y = MT + i * cell
        for j, c in enumerate(row):
            if not c:
                continue
            s += '<rect x="%.2f" y="%.2f" width="%.2f" height="%.2f" fill="%s"/>' % (
                ML + j * cell, y, cell + 0.4, cell + 0.4, _heat_ramp(c / peak, base, floor))
    s += _dot_guides(m, ML, MT, plot, T, guides)
    s += _dot_legend(ML, MT + plot + 28, T, heat=(base, peak), px=px, floor=floor)
    if foot:
        s += _dot_footer_svg(ML, MT + plot + 44, T, foot, px)
    return s + "</svg>"


# ---------------- genome TE landscape (whole-genome annotation) ----------------
# Class hues reuse the domain palette so a class means the same colour everywhere in the app. Tandem and
# "other" repeats are drawn in neutral greys: they are real repeats but NOT transposable elements, and the
# figure must not let them read as TE content.
_TE_CLASS_COLOUR = {"LTR": OK["LTR"], "LINE": OK["RT"], "SINE": OK["RNaseH"], "DNA": OK["TIR"],
                    "RC": OK["GAG"], "Retroposon": OK["PR"], "PLE": OK["RT"]}
_NONTE_COLOUR = {"tandem": "#9aa3ad", "other": "#c3ccd6"}


def _te_library_footer(r: dict) -> tuple[str, str | None]:
    """Return library-specific figure footnotes; custom-library output must never imply a Dfam search."""
    kind = str(r.get("library_kind") or "").strip()
    raw_custom = r.get("custom_library")
    custom = raw_custom if isinstance(raw_custom, dict) else {}
    is_custom = bool(custom) or kind.lower() == "user" or kind.lower().startswith("user-supplied")
    if is_custom:
        name = str(custom.get("name") or "unnamed user FASTA")
        sha = str(custom.get("sha256") or "unavailable")
        return (f"Interspersed-repeat homology used user-supplied FASTA {name}; absent repeat models cannot be matched.",
                f"Library SHA-256: {sha}. TE coverage requires recognized #class/family headers; simple/low-complexity is separate.")
    label = kind or "installed Dfam partitions"
    version = r.get("dfam_version")
    if version and "dfam" not in label.lower():
        label += f" (Dfam {version})"
    return (f"Homology to {label}; families absent from it cannot appear.", None)


_TE_NOTE = ("Independent union bars: TE, tandem and other/unclassified can overlap and must not be summed. "
            "All repeats and not annotated are complementary.")


def svg_te_composition(r: dict, W: float = 900, theme: str = "light") -> str:
    """Independent union-coverage bars for the whole-genome RepeatMasker result.

    Repeat-name/model aggregates are deliberately not stacked: intervals assigned to different models can overlap.
    The all-repeat union and its unannotated complement form one partition; TE, tandem and other are separate
    kind-specific unions that may overlap one another and therefore must not be added."""
    gbp = float(r.get("genome_bp") or 0)
    # honour the panel's background mode: the figure carries text and a frame, and painting them in a
    # fixed dark ink left the whole figure unreadable on the dark track (and made the DARK/WHITE buttons
    # look broken, because nothing changed when they were pressed).
    dark = str(theme).lower() in ("dark", "uv", "mono-dark")
    paper = "#12161b" if dark else "#ffffff"
    ink = "#e6edf3" if dark else "#333333"
    faint = "#8b929a" if dark else "#666666"
    rule = "#5a636c" if dark else "#8b929a"
    if gbp <= 0:
        return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="60">'
                f'<rect width="{W}" height="60" fill="{paper}"/></svg>')

    def union_bp(field: str) -> float:
        try:
            value = float(r.get(field) or 0)
        except (TypeError, ValueError):
            value = 0.0
        return min(gbp, max(0.0, value))

    total_bp = union_bp("total_bp")
    rest = max(0.0, gbp - total_bp)  # only the all-repeat union has a valid assembly complement
    rows = [
        ("All repeats (union)", total_bp, "#4c78a8", False),
        ("Transposable elements (union)", union_bp("te_bp"), _TE_CLASS_COLOUR["LTR"], False),
        ("Tandem repeats (union)", union_bp("tandem_bp"), _NONTE_COLOUR["tandem"], False),
        ("Other/unclassified annotations (union)", union_bp("other_bp"), _NONTE_COLOUR["other"], False),
        ("Not annotated as repeat", rest, "none", True),
    ]
    footer1, footer2 = _te_library_footer(r)
    barH, gap, top, padx, lab_px, fpx, lead = 18, 8, 42, 12, 10.5, 10, 12
    # the label column is sized to the longest "label pct%" string (mono metric, an upper bound), so the percentage
    # is never painted over by its own bar
    texts = [f"{label} {100.0 * bp / gbp:.2f}%" for label, bp, _c, _d in rows]
    labelW = math.ceil(max(_text_w(t, lab_px) for t in texts)) + 12
    barX, barW = padx + labelW, max(1.0, W - (padx + labelW) - padx)
    # caption lines wrapped to the figure width (a long FASTA name or checksum used to run off the edge); sentences
    # are kept whole where they can be, and the footers sit at the bottom edge as before
    wrapw = W - 2 * padx
    note_lines = _wrap_mono(_TE_NOTE.split(" "), wrapw, fpx, sep=" ")
    foot_lines = _wrap_mono(footer1.split(" "), wrapw, fpx, sep=" ")
    if footer2:
        foot_lines += _wrap_mono(re.split(r"(?<=\.) +", footer2), wrapw, fpx, sep=" ")
    note_y = top + len(rows) * (barH + gap) + 3
    H = max(244 if footer2 else 230, note_y + len(note_lines) * lead + len(foot_lines) * lead + 12)
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
           f'viewBox="0 0 {W} {H}" font-family="{FIGFONT}">',
           f'<rect width="{W}" height="{H}" fill="{paper}"/>',
           f'<text x="{padx}" y="19" font-size="12" fill="{ink}">Share of the {gbp/1e6:.1f} Mb assembly '
           f'covered by independent interval unions</text>']
    for i, (label, bp, colour, dashed) in enumerate(rows):
        y = top + i * (barH + gap)
        pct = 100.0 * bp / gbp
        width = (bp / gbp) * barW
        out.append(f'<text x="{padx}" y="{y+13}" font-size="{lab_px}" fill="{ink}">{esc(label)} {pct:.2f}%</text>')
        out.append(f'<rect x="{barX}" y="{y}" width="{barW:.2f}" height="{barH}" fill="none" stroke="{rule}"/>')
        attrs = (f'fill="none" stroke="{rule}" stroke-dasharray="3 3"' if dashed
                 else f'fill="{colour}"')
        out.append(f'<rect x="{barX}" y="{y}" width="{width:.2f}" height="{barH}" {attrs}>'
                   f'<title>{esc(label)}: {pct:.2f}% ({int(round(bp)):,} bp)</title></rect>')
    for k, ln in enumerate(note_lines):
        out.append(f'<text x="{padx}" y="{note_y + k * lead}" font-size="{fpx}" fill="{faint}">{esc(ln)}</text>')
    for k, ln in enumerate(foot_lines):
        out.append(f'<text x="{padx}" y="{H - 5 - (len(foot_lines) - 1 - k) * lead}" font-size="{fpx}" '
                   f'fill="{faint}">{esc(ln)}</text>')
    out.append("</svg>")
    return "".join(out)
