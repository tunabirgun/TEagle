"""Publication figure builders, ported from the validated web UI (app/web/app.js).
Pure functions returning SVG strings. The native app renders them on-screen via QSvgRenderer
and exports the identical string as SVG / rasterises it to PNG — so figures match the web output.

Two figures: an interactive genome viewer (ruler + tracks, windowed) and a to-scale
multi-lane agarose gel (log MW axis, 5 palettes). Geometry mirrors svgGenome / svgGel exactly."""
from __future__ import annotations
import math

# every colour used here is defined once, in theme.py (single source of truth); re-exported so
# `figures.OK` / `figures.GELPAL` stay valid for existing callers and tests.
from theme import OKABE_ITO, OK, GENECOL, ARCHCOL, CISCOL, GV_THEME, GELPAL   # noqa: F401

FIGFONT = "Cascadia Mono, Consolas, monospace"   # bundled UI font; mono digits align on the gel ladder / ruler


def esc(s) -> str:
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def _label_ink(hex_color: str) -> str:
    """Pick black/white label by best WCAG contrast against a fill colour. Accepts 3- or 6-digit hex
    (a bare '#888' fallback would otherwise raise ValueError and abort the whole figure render)."""
    def lin(v):
        v /= 255.0
        return v / 12.92 if v <= 0.03928 else ((v + 0.055) / 1.055) ** 2.4
    h = hex_color.lstrip("#")
    if len(h) == 3:
        h = "".join(c * 2 for c in h)                         # #888 -> #888888
    r = int(h[0:2], 16); g = int(h[2:4], 16); b = int(h[4:6], 16)
    L = 0.2126 * lin(r) + 0.7152 * lin(g) + 0.0722 * lin(b)
    return "#fff" if (1.05 / (L + 0.05)) >= ((L + 0.05) / 0.05) else "#111"


def _fmt_int(n) -> str:
    return f"{int(round(n)):,}"


def _fit_caption(items: list, plotW: float) -> str:
    """One caption line that fits the plot width (mono ≈4.5 px/char at 7.5 px); overflow is counted, not dropped."""
    maxc = max(24, int(plotW / 4.5))
    out, n = [], 0
    for it in items:
        if out and n + len(it) + 3 > maxc:
            return " · ".join(out) + f" · +{len(items)-len(out)} more"
        out.append(it); n += len(it) + 3
    return " · ".join(out)


# ================= genome viewer =================
GV_NAME_PX = 9.5      # track-name font size in the viewer
GV_CHAR_EM = 0.6      # Cascadia Mono advance width (measured 0.600 em; the SVG renderer draws no wider)
GV_NAME_GAP = 10      # gap between the end of a track name and the plot area
GV_NAME_PAD = 4       # breathing room at the figure's left edge
GV_ML_MIN = 96        # gutter floor, so short-named models keep the historical layout
GV_MR = 16            # right margin


def gv_gutter(model: dict, W: float = None) -> float:
    """Width of the left track-label gutter: sized to the longest track name in this model, so a
    name like 'env mRNA (predicted)' is never clipped by the figure's left edge. Mono metric, so the
    result is an upper bound on the rendered text width. With W, the gutter is capped so the plot
    area keeps its 120 px floor on a very narrow figure. One definition — svg_genome draws it and the
    interactive canvas maps pixels to bp with it, so the two can never disagree."""
    n = max((len(t.get("name", "")) for t in (model or {}).get("tracks", []) if t.get("name")), default=0)
    g = float(max(GV_ML_MIN, math.ceil(n * GV_NAME_PX * GV_CHAR_EM) + GV_NAME_GAP + GV_NAME_PAD))
    return g if W is None else min(g, max(float(GV_ML_MIN), W - GV_MR - 120))


def gv_nice_step(span: float, ticks: int) -> float:
    raw = max(1.0, span / ticks)
    mag = 10 ** math.floor(math.log10(raw))
    n = raw / mag
    return (1 if n < 1.5 else 2 if n < 3.5 else 5 if n < 7.5 else 10) * mag


def gv_tracks_from_rec(rec: dict) -> dict:
    """Build a genome-viewer model from an analysis record (structural + domains + ORFs + ERV architecture)."""
    tracks, reps, cis, tails = [], [], [], []
    for e in rec.get("structural", []):
        t = e["type"]
        if t.startswith("LTR") or t.startswith("TIR"):
            col = OK["LTR"] if t.startswith("LTR") else OK["TIR"]
            for p in (e.get("five_prime"), e.get("three_prime")):
                if p:
                    reps.append({"start": p[0], "end": p[1], "color": col, "label": t.split(" ")[0],
                                 "kind": t.split(" ")[0], "tip": f"{t} {p[0]}–{p[1]}"})
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
            cis.append({"start": p[0], "end": p[1], "color": CISCOL[key], "label": lab, "kind": key, "tip": tip})
        elif e.get("pos"):
            p = e["pos"]
            f = {"start": p[0], "end": p[1], "color": OK["tail"], "label": t.split(" ")[0],
                 "kind": t.split(" ")[0], "tip": f"{t} {p[0]}–{p[1]}"}
            (tails if t.startswith("poly") else reps).append(f)   # a poly-A/T tail is not a terminal repeat
    if reps:
        tracks.append({"name": "terminal repeats", "height": 20, "features": reps})
    if tails:
        tracks.append({"name": "poly-A/T tail", "height": 18, "features": tails})
    if cis:
        tracks.append({"name": "cis-elements", "height": 18, "features": cis})
    doms = [{"start": d["nt"][0], "end": d["nt"][1], "color": OK.get(d["domain"], "#888"),
             "label": d["domain"],
             "tip": f"{d['domain']} · {d.get('label','')} · nt {d['nt'][0]}–{d['nt'][1]} · score {d.get('score')}"}
            for d in rec.get("domains", [])]
    if doms:
        tracks.append({"name": "protein domains", "height": 22, "features": doms})
    arch = rec.get("retroviral")
    if arch:                                                  # ERV: env expressed from a spliced subgenomic mRNA
        feat = [{"start": arch["leader_exon"][0], "end": arch["leader_exon"][1], "color": ARCHCOL["exon"],
                 "label": "leader", "tip": f"5′ leader exon {arch['leader_exon'][0]}–{arch['leader_exon'][1]}"},
                {"start": arch["intron"][0], "end": arch["intron"][1], "color": ARCHCOL["intron"], "intron": True,
                 "label": "gag–pol", "tip": f"gag–pro–pol intron (fused polyprotein, removed) {arch['intron'][0]}–{arch['intron'][1]}"},
                {"start": arch["env_exon"][0], "end": arch["env_exon"][1], "color": ARCHCOL["exon"],
                 "label": "env", "tip": f"env exon {arch['env_exon'][0]}–{arch['env_exon'][1]}"}]
        tracks.append({"name": "env mRNA (predicted)", "height": 20, "features": feat})
    orfs = [{"start": o["start"], "end": o["end"], "color": OK["ORF"], "strand": o["strand"],
             "tip": f"ORF {o['strand']}{o['frame']} · {o['length_aa']} aa"} for o in rec.get("orfs", [])]
    if orfs:
        tracks.append({"name": "ORFs (± strand)", "height": 26, "features": orfs, "stranded": True})
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
             "label": "exon*" if e.get("derived") else "exon",
             "tip": f"exon {e['start']}–{e['end']} ({e['end']-e['start']} bp)"
                    + (" · derived from the record's CDS/mRNA — not a separate exon annotation" if e.get("derived") else "")}
            for e in gm.get("exons", [])]
    feat += [{"start": i["start"], "end": i["end"], "color": GENECOL["intron"], "intron": True,
              "tip": f"intron {i['start']}–{i['end']}" + (
                  f" · {i['donor']}…{i['acceptor']}{' (canonical)' if i.get('canonical') else ''}"
                  if i.get("donor") else "")} for i in gm.get("introns", [])]
    if include_flanks:                                        # 5'/3' flanks + interior gaps, clickable for copy/design
        for r in _flanks_and_gaps(gm.get("exons", []), gm.get("introns", []), length):
            feat.append({"start": r["start"], "end": r["end"],
                         "color": GENECOL["gap"] if r["kind"] == "gap" else GENECOL["flank"],
                         "label": r["name"],                  # readable on-glyph name; FASTA id is sanitised in _feat_menu
                         "tip": f"{r['name']} {r['start']}–{r['end']} ({r['end']-r['start']} bp) · not exon/intron"})
    if feat:
        tracks.append({"name": "exons / introns", "height": 22, "features": feat})
    cds = [{"start": c["start"], "end": c["end"], "color": GENECOL["cds"], "label": "CDS",
            "tip": f"CDS {c['start']}–{c['end']}"} for c in gm.get("cds", [])]
    if cds:
        tracks.append({"name": "CDS (coding)", "height": 16, "features": cds})
    return {"length": length or 1, "tracks": tracks}


def _gv_theme(theme: str, for_export: bool) -> dict:
    """Genome-viewer chrome for a render target; the palettes live in theme.GV_THEME (copy, so a
    caller mutating the returned dict cannot corrupt the shared one)."""
    return dict(GV_THEME["export" if for_export else ("white" if theme == "white" else "dark")])


def svg_genome(model: dict, view: dict, W: float, theme: str, for_export: bool = False, return_regions: bool = False):
    """Render the genome viewer SVG. With return_regions=True, also return a list of on-screen
    feature hit-boxes (SVG coords + feature identity) so the interactive canvas can hover/right-click."""
    L = model["length"]
    regions = []
    MR, ovH, rulerH, MT = GV_MR, 13, 24, 34
    ML = gv_gutter(model, W)                        # gutter fits the widest track name (no clipped label)
    plotW = max(120, W - ML - MR)
    s0, s1 = view["start"], view["end"]
    span = max(1, s1 - s0)
    def bx(bp): return ML + (bp - s0) / span * plotW
    def ox(bp): return ML + bp / L * plotW
    T = _gv_theme(theme, for_export)
    y = MT + ovH + 10 + rulerH
    track_ys = []
    for t in model["tracks"]:
        track_ys.append(y)
        y += (t.get("height", 20)) + 20
    H = y + 12
    s = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" '
         f'font-family="{FIGFONT}">')
    s += (f'<defs><linearGradient id="gvwin" x1="0" x2="0" y1="0" y2="1">'
          f'<stop offset="0" stop-color="{T["win"]}" stop-opacity="0.28"/>'
          f'<stop offset="1" stop-color="{T["win"]}" stop-opacity="0.10"/></linearGradient></defs>')
    if T["paper"] != "none":
        s += f'<rect width="{W}" height="{H}" fill="{T["paper"]}"/>'
    ovY = MT
    s += f'<text x="{ML}" y="{ovY-5}" fill="{T["faint"]}" font-size="8.5">whole element · {_fmt_int(L)} bp</text>'
    s += f'<rect x="{ML}" y="{ovY}" width="{plotW}" height="{ovH}" rx="3" fill="{T["lane"]}" stroke="{T["frame"]}"/>'
    for t in model["tracks"]:
        for f in t["features"]:
            a = ox(f["start"]); w = max(ox(f["end"]) - ox(f["start"]), 1)
            s += f'<rect x="{a:.1f}" y="{ovY+2}" width="{w:.1f}" height="{ovH-4}" rx="1" fill="{f["color"]}" opacity="0.6"/>'
    s += (f'<rect x="{ox(s0):.1f}" y="{ovY-2:.1f}" width="{max(ox(s1)-ox(s0),2):.1f}" height="{ovH+4}" '
          f'rx="2" fill="url(#gvwin)" stroke="{T["win"]}" stroke-width="1.2"/>')
    ry = MT + ovH + 10 + rulerH - 7
    step = gv_nice_step(span, 7)
    first = math.ceil(s0 / step) * step
    s += f'<line x1="{ML}" y1="{ry}" x2="{ML+plotW}" y2="{ry}" stroke="{T["frame"]}"/>'
    bp = first
    while bp <= s1 + 1:
        x = bx(bp)
        if not (x < ML - 1 or x > ML + plotW + 1):
            lab = f"{bp/1e6:g}M" if bp >= 1e6 else (f"{bp/1000:g}k" if bp >= 1000 else f"{bp:g}")
            s += (f'<line x1="{x:.1f}" y1="{ry}" x2="{x:.1f}" y2="{ry-5}" stroke="{T["faint"]}"/>'
                  f'<text x="{x:.1f}" y="{ry-8:.1f}" fill="{T["faint"]}" font-size="9" text-anchor="middle">{lab}</text>'
                  f'<line x1="{x:.1f}" y1="{ry+2:.1f}" x2="{x:.1f}" y2="{H-10}" stroke="{T["grid"]}"/>')
        bp += step
    for ti, t in enumerate(model["tracks"]):
        ty = track_ys[ti]; th = t.get("height", 20)
        s += (f'<text x="{ML-GV_NAME_GAP}" y="{ty+th/2+3:.1f}" fill="{T["faint"]}" font-size="{GV_NAME_PX}" '
              f'text-anchor="end">{esc(t["name"])}</text>')
        s += f'<rect x="{ML}" y="{ty}" width="{plotW}" height="{th}" rx="3" fill="{T["track"]}"/>'
        stranded = t.get("stranded")
        if stranded:
            s += f'<line x1="{ML}" y1="{ty+th/2:.1f}" x2="{ML+plotW}" y2="{ty+th/2:.1f}" stroke="{T["grid"]}"/>'
        tiny = []                                              # labels of features too narrow to letter in place
        for f in t["features"]:
            a = max(bx(f["start"]), ML); b = min(bx(f["end"]), ML + plotW)
            if b < ML - 0.5 or a > ML + plotW + 0.5:
                continue
            if f.get("intron"):
                mid = (a + b) / 2
                s += (f'<path d="M {a:.1f} {ty+th/2:.1f} L {mid:.1f} {ty+3:.1f} L {b:.1f} {ty+th/2:.1f}" '
                      f'fill="none" stroke="{f["color"]}" stroke-width="1.3"><title>{esc(f.get("tip",""))}</title></path>')
                regions.append({"x0": a, "y0": ty, "x1": b, "y1": ty + th, "start": f["start"], "end": f["end"],
                                "strand": f.get("strand", "+"), "label": f.get("label", "intron"), "tip": f.get("tip", "")})
                continue
            w = max(b - a, 1.5)
            yy = (ty + 2.5 if f.get("strand") == "+" else ty + th / 2 + 1.5) if stranded else ty + 2.5
            hh = th / 2 - 4 if stranded else th - 5
            s += (f'<rect class="gvglyph" x="{a:.1f}" y="{yy:.1f}" width="{w:.1f}" height="{max(hh,3):.1f}" '
                  f'rx="2.5" fill="{f["color"]}"><title>{esc(f.get("tip",""))}</title></rect>')
            regions.append({"x0": a, "y0": yy, "x1": a + w, "y1": yy + max(hh, 3), "start": f["start"], "end": f["end"],
                            "strand": f.get("strand", "+"), "label": f.get("label", ""), "tip": f.get("tip", "")})
            if f.get("label"):
                if w > 26:
                    s += (f'<text x="{a+4:.1f}" y="{yy+max(hh,3)-3:.1f}" fill="{_label_ink(f["color"])}" '
                          f'font-size="9" font-weight="700" pointer-events="none">{esc(f["label"])}</text>')
                else:                                          # too narrow to letter: leader line + caption below
                    cx = a + w / 2
                    s += (f'<line x1="{cx:.1f}" y1="{yy+max(hh,3):.1f}" x2="{cx:.1f}" y2="{ty+th+2:.1f}" '
                          f'stroke="{f["color"]}" stroke-width="0.6"/>')
                    tiny.append(f'{f["label"]} {_fmt_int(f["start"])}–{_fmt_int(f["end"])}')
        if tiny:                                               # never drop a label (PBS·? hedge is often ~2 px wide)
            s += (f'<text x="{ML}" y="{ty+th+10:.1f}" fill="{T["faint"]}" font-size="7.5">'
                  f'{esc(_fit_caption(tiny, plotW))}</text>')
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


def _lane_bands(amplicons, P, has_locus=True):
    """Collapse a lane's amplicons into one band per product size (a real gel cannot resolve equal lengths).
    Band intensity follows the strongest (fewest-mismatch) product at that size.

    has_locus=True (a designed on-target exists): worst-case colour wins — a band whose size also carries an
    off-target (or single-primer artefact) is NOT a clean on-target, so it reads off-target (a specificity
    warning), never a reassuring on-target band. has_locus=False (whole-genome scan of a bare consensus pair):
    there is no on/off target, so every product is a NEUTRAL 'genomic priming site' and must not read as the
    off-target warning colour. Every product is still enumerated in the table below the gel."""
    groups = {}
    for a in (amplicons or []):
        groups.setdefault(a["length"], []).append(a)
    bands = []
    for size in sorted(groups):
        g = groups[size]
        n_on, n_off, n_single = _call_counts(g)           # disjoint buckets: no count can go negative
        if not has_locus:                                 # neutral priming sites (single-primer artefacts still distinct)
            n_site = len(g) - n_single
            call = "single" if (n_single and not n_site) else "site"
            parts = (([f"{n_site} priming site" + ("s" if n_site != 1 else "")] if n_site else [])
                     + ([f"{n_single} single-primer"] if n_single else []))
            on = False
        else:
            if n_off:                                   # any off-target co-migrating here -> flag the whole band off-target
                call = "off"
            elif n_single:
                call = "single"
            elif n_on:                                    # purely on-target at this size
                call = "on"
            else:
                call = "off"
            parts = (([f"{n_on} on-target"] if n_on else []) + ([f"{n_off} off-target"] if n_off else [])
                     + ([f"{n_single} single-primer"] if n_single else []))
            on = bool(n_on)
        min_mm = min((a.get("fwd_mm", 0) + a.get("rev_mm", 0)) for a in g)
        src = esc(g[0].get("source", ""))
        color = P[call]                                   # colour derives from the call, so the shape cue can never disagree
        bands.append({"size": size, "color": color, "call": call, "opacity": _band_opacity(min_mm),
                      "on": on, "single": bool(n_single), "count": len(g),
                      "t": ", ".join(parts) + (f" · {src}" if src else "")})
    return bands


LANES_PER_ROW = 10        # sample lanes per gel row; more than this wraps onto stacked rows


def _gel_geometry(data: dict):
    """Shared layout math for svg_gel and gel_regions (must stay identical so hit-boxes line up).
    Sample lanes wrap into stacked rows of <=LANES_PER_ROW; each row carries its own ladder + bp scale."""
    lanes = data.get("lanes") or [{"label": data.get("laneLabel", "PCR"), "amplicons": data.get("amplicons", [])}]
    sizes = [a["length"] for l in lanes for a in (l.get("amplicons") or [])]
    smallest = min(sizes) if sizes else 90
    minbp = max(25, min(90, smallest - 10))
    maxbp = max([1600] + sizes)
    LADDER = [m for m in (1500, 1000, 700, 500, 400, 300, 200, 100, 50) if minbp <= m <= maxbp]
    laneW, gap, x0 = 40, 12, 62
    TOP_MARGIN, ROW_TOP_PAD, BODY, ROW_BOT_PAD, ROW_GAP = 8, 30, 240, 20, 16
    ROW_PITCH = ROW_TOP_PAD + BODY + ROW_BOT_PAD + ROW_GAP
    rows = [list(range(i, min(i + LANES_PER_ROW, len(lanes)))) for i in range(0, len(lanes), LANES_PER_ROW)] or [[]]
    widest = min(LANES_PER_ROW, max(1, len(lanes)))
    any_single = any(a.get("single_primer") for l in lanes for a in (l.get("amplicons") or []))
    legend_w = x0 + (214 if any_single else 132) + 100    # room for the on/off(/single-primer)/ladder legend row
    W = max(x0 + (1 + widest) * (laneW + gap) + 12, legend_w, 300)
    last_bot = TOP_MARGIN + (len(rows) - 1) * ROW_PITCH + ROW_TOP_PAD + BODY + ROW_BOT_PAD
    H = last_bot + 22                                       # legend strip below the last row

    def row_top(r):                                        # band-area top of row r
        return TOP_MARGIN + r * ROW_PITCH + ROW_TOP_PAD

    def y(bp, r):
        t = row_top(r)
        return t + (math.log(maxbp) - math.log(max(bp, minbp))) / (math.log(maxbp) - math.log(minbp)) * BODY

    def laneX(col):
        return x0 + col * (laneW + gap)

    return {"lanes": lanes, "rows": rows, "minbp": minbp, "maxbp": maxbp, "LADDER": LADDER,
            "laneW": laneW, "gap": gap, "x0": x0, "BODY": BODY, "W": W, "H": H,
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
            groups = {}
            for a in (l.get("amplicons") or []):
                groups.setdefault(a["length"], []).append(a)
            for size in sorted(groups):
                g = groups[size]
                yy = y(size, r)
                # representative for the right-click menu: the intended (on-target) product if any, else the strongest
                rep = min(g, key=lambda a: (not a.get("on_target"), a.get("fwd_mm", 0) + a.get("rev_mm", 0)))
                n_on, n_off, n_single = _call_counts(g)   # same disjoint split as the band tooltip
                if l.get("has_locus", True):              # with a design locus -> on/off-target; else neutral priming sites
                    call = ", ".join(([f"{n_on} on-target"] if n_on else []) + ([f"{n_off} off-target"] if n_off else [])
                                     + ([f"{n_single} single-primer"] if n_single else []))
                else:
                    n_site = len(g) - n_single
                    call = ", ".join(([f"{n_site} priming site" + ("s" if n_site != 1 else "")] if n_site else [])
                                     + ([f"{n_single} single-primer"] if n_single else []))
                tip = f'{size} bp · {call}' + (f' · {rep.get("source","")}' if rep.get("source") else "")
                regions.append({"x0": lx + 3, "y0": yy - 4, "x1": lx + laneW - 3, "y1": yy + 4,
                                "tip": tip, "amplicon": rep, "pair": l.get("label", ""),
                                "has_locus": l.get("has_locus", True)})
    return regions


def svg_gel(data: dict, bg: str) -> str:
    """data: {lanes:[{label, amplicons:[{length,on_target,source}]}]} or legacy {amplicons}.
    Lanes beyond LANES_PER_ROW wrap onto stacked gel rows, each with its own ladder + bp scale."""
    G = _gel_geometry(data)
    lanes, rows, LADDER = G["lanes"], G["rows"], G["LADDER"]
    laneW, gap, x0, BODY = G["laneW"], G["gap"], G["x0"], G["BODY"]
    W, H, y, laneX, row_top = G["W"], G["H"], G["y"], G["laneX"], G["row_top"]
    P = GELPAL.get(bg, GELPAL["transparent"])
    s = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" '
         f'font-family="{FIGFONT}">')
    s += (f'<defs><filter id="glow" x="-40%" y="-140%" width="180%" height="380%">'
          f'<feGaussianBlur stdDeviation="{P["glow"]}" result="b"/><feMerge>'
          f'<feMergeNode in="b"/><feMergeNode in="SourceGraphic"/></feMerge></filter></defs>')
    if P["paper"] != "none":
        s += f'<rect width="{W}" height="{H}" fill="{P["paper"]}"/>'

    def draw_lane(col, r, label, bands, is_ladder, advisory=False):
        nonlocal s
        lx = laneX(col)
        rtop = row_top(r); rbot = rtop + BODY
        s += f'<rect x="{lx+4:.1f}" y="{rtop-13}" width="{laneW-8}" height="4" rx="1" fill="{P["well"]}"/>'
        s += f'<text x="{lx+laneW/2:.1f}" y="{rtop-18}" fill="{P["ink"]}" font-size="8.5" text-anchor="middle">{esc(label)}</text>'
        for bd in (bands or []):
            yy = y(bd["size"], r); h = 1.6 if is_ladder else P["band"]
            filt = "" if is_ladder else ' filter="url(#glow)"'
            op = "" if is_ladder else f' fill-opacity="{bd.get("opacity", 1.0)}"'
            dash = ' stroke="#ffffff" stroke-width="0.5" stroke-dasharray="2 1.5"' if bd.get("single") else ""
            n = bd.get("count", 1)
            title = f'{bd["size"]} bp' + (f' · {bd["t"]}' if bd.get("t") else "") + (f' (×{n})' if n > 1 else "")
            s += (f'<rect x="{lx+3:.1f}" y="{yy-h/2:.1f}" width="{laneW-6}" height="{h}" rx="1" '
                  f'fill="{bd["color"]}"{op}{dash}{filt}><title>{title}</title></rect>')
            if not is_ladder and bd.get("call"):          # shape marker: the call read without colour
                s += _call_mark(lx - 1.0, yy, bd["call"], bd["color"])
        if not is_ladder:
            if not (bands or []):
                s += f'<text x="{lx+laneW/2:.1f}" y="{rbot+13}" fill="{P["ink"]}" font-size="7" text-anchor="middle">—</text>'
            elif not advisory and not any(b.get("on") for b in bands):   # bands present but none intended
                s += (f'<text x="{lx+laneW/2:.1f}" y="{rbot+13}" fill="{P["off"]}" font-size="6.5" '   # (a genome scan is
                      f'text-anchor="middle">no on-target</text>')                                     # all-off-target by design)

    for r, idxs in enumerate(rows):
        rtop = row_top(r)
        ncols = 1 + len(idxs)
        s += (f'<rect x="{x0-7:.1f}" y="{rtop-16}" width="{ncols*(laneW+gap)+2:.1f}" height="{BODY+30:.1f}" '
              f'rx="2" fill="{P["gel"]}" stroke="{P["stroke"]}"/>')
        s += f'<text x="{x0-14}" y="{rtop-19}" fill="{P["ink"]}" font-size="8" text-anchor="end">bp</text>'
        for m in LADDER:
            s += f'<text x="{x0-14}" y="{y(m, r)+2.5:.1f}" fill="{P["ink"]}" font-size="8" text-anchor="end">{m}</text>'
        draw_lane(0, r, "L", [{"size": m, "color": P["ladder"]} for m in LADDER], True)
        for j, li in enumerate(idxs):
            l = lanes[li]
            draw_lane(j + 1, r, l["label"], _lane_bands(l.get("amplicons") or [], P, l.get("has_locus", True)),
                      False, l.get("advisory", False))

    any_single = any(a.get("single_primer") for l in lanes for a in (l.get("amplicons") or []))
    all_neutral = bool(lanes) and all(not l.get("has_locus", True) for l in lanes)   # every lane a no-locus scan
    ly = H - 8
    if all_neutral:                          # neutral 'priming site' bands -> a matching swatch, not on/off (which match nothing)
        s += (_call_mark(x0, ly, "site", P["site"], 3)
              + f'<text x="{x0+7}" y="{ly+3}" fill="{P["ink"]}" font-size="8">priming site</text>')
        xnext = x0 + 78
    else:                                    # any locus (incl. local PCR, which shares this gel) -> on/off swatches
        s += (_call_mark(x0, ly, "on", P["on"], 3)            # swatch shapes = the in-lane band markers
              + f'<text x="{x0+7}" y="{ly+3}" fill="{P["ink"]}" font-size="8">on-target</text>'
              + _call_mark(x0 + 64, ly, "off", P["off"], 3)
              + f'<text x="{x0+71}" y="{ly+3}" fill="{P["ink"]}" font-size="8">off-target</text>')
        xnext = x0 + 132
    if any_single:
        s += (_call_mark(xnext, ly, "single", P["single"], 3)
              + f'<text x="{xnext+7}" y="{ly+3}" fill="{P["ink"]}" font-size="8">single-primer</text>')
        xnext += 82
    s += f'<text x="{xnext}" y="{ly+3}" fill="{P["ink"]}" font-size="8">L = MW ladder (bp)</text>'
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


def _dot_theme(theme: str, for_export: bool) -> dict:
    return dict(_DOT_THEME["export" if for_export else ("white" if theme == "white" else "dark")])


def _heat_ramp(t: float, base: str) -> str:
    """Sequential ramp from a near-white floor to the signal hue, monotone in lightness so the scale
    still reads in greyscale — which an unordered rainbow would not."""
    t = 0.0 if t < 0 else (1.0 if t > 1 else t)
    r0, g0, b0 = 245, 247, 249
    r1, g1, b1 = int(base[1:3], 16), int(base[3:5], 16), int(base[5:7], 16)
    return "#%02X%02X%02X" % (int(r0 + (r1 - r0) * t), int(g0 + (g1 - g0) * t), int(b0 + (b1 - b0) * t))


def _dot_frame(m, W, T, title, sub):
    """Shared chrome: a square plot box with a bp ruler on both axes, title and scope line."""
    ML, MT, MR, MB = 58.0, 46.0, 16.0, 66.0
    plot = max(W - ML - MR, 80.0)
    H = MT + plot + MB
    n = max(m["length"], 1)
    s = ('<svg xmlns="http://www.w3.org/2000/svg" width="%.0f" height="%.0f" viewBox="0 0 %.0f %.0f" '
         'font-family="%s">' % (W, H, W, H, FIGFONT))
    s += '<rect width="%.0f" height="%.0f" fill="%s"/>' % (W, H, T["bg"])
    s += ('<text x="%.1f" y="20" fill="%s" font-size="12" font-weight="700">%s</text>'
          '<text x="%.1f" y="34" fill="%s" font-size="9">%s</text>'
          % (ML, T["ink"], esc(title), ML, T["faint"], esc(sub)))
    s += ('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="none" stroke="%s" stroke-width="1"/>'
          % (ML, MT, plot, plot, T["grid"]))
    step = gv_nice_step(n, 6)
    t = 0
    while t <= n:
        x = ML + plot * (t / n)
        y = MT + plot * (t / n)
        s += ('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="0.5" opacity="0.55"/>'
              '<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="0.5" opacity="0.55"/>'
              '<text x="%.1f" y="%.1f" fill="%s" font-size="8" text-anchor="middle">%s</text>'
              '<text x="%.1f" y="%.1f" fill="%s" font-size="8" text-anchor="end">%s</text>'
              % (x, MT, x, MT + plot, T["grid"], ML, y, ML + plot, y, T["grid"],
                 x, MT + plot + 13, T["faint"], _fmt_int(t),
                 ML - 6, y + 3, T["faint"], _fmt_int(t)))
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


def _dot_legend(x, y, T, heat=None, fwd=_DOT_FWD, rev=_DOT_REV):
    if heat is None:
        return ('<circle cx="%.1f" cy="%.1f" r="3" fill="%s"/>'
                '<text x="%.1f" y="%.1f" fill="%s" font-size="9">direct repeat (forward)</text>'
                '<rect x="%.1f" y="%.1f" width="6" height="6" fill="%s"/>'
                '<text x="%.1f" y="%.1f" fill="%s" font-size="9">inverted repeat (reverse complement)</text>'
                % (x + 4, y - 3, fwd, x + 13, y, T["ink"],
                   x + 176, y - 6, rev, x + 187, y, T["ink"]))
    base, peak = heat
    s = ""
    for i in range(24):
        s += '<rect x="%.1f" y="%.1f" width="7" height="8" fill="%s"/>' % (
            x + i * 7, y - 8, _heat_ramp(i / 23.0, base))
    s += ('<text x="%.1f" y="%.1f" fill="%s" font-size="8">0</text>'
          '<text x="%.1f" y="%.1f" fill="%s" font-size="8">%d matches / bin</text>'
          % (x, y + 10, T["faint"], x + 176, y - 1, T["faint"], peak))
    return s


def svg_dotplot(m: dict, W: float = 620, theme: str = "dark", for_export: bool = False,
                guides: list = None, title: str = "Self-similarity dot plot",
                read_threshold: int = 1, fwd: str = None, rev: str = None) -> str:
    """Binary dot matrix: a mark wherever the locus matches itself."""
    T = _dot_theme(theme, for_export)
    fwd, rev = fwd or _DOT_FWD, rev or _DOT_REV       # caller may override the mark hues; bg stays white on export
    sub = "exact %d-mer matches · %s bp" % (m["k"], _fmt_int(m["length"]))
    if read_threshold and read_threshold > 1:
        sub += " · cells below %d matches shown faint (chance level)" % read_threshold
    s, ML, MT, plot, H = _dot_frame(m, W, T, title, sub)
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
    s += _dot_legend(ML, MT + plot + 28, T, heat=None, fwd=fwd, rev=rev)
    if for_export:
        s += '<text x="%.1f" y="%.1f" fill="%s" font-size="7.5">%s</text>' % (
            ML, MT + plot + 44, T["faint"], esc(_DOT_FOOTER))
    return s + "</svg>"


def svg_dotheat(m: dict, W: float = 620, theme: str = "dark", for_export: bool = False,
                guides: list = None, which: str = "forward",
                title: str = "Self-similarity heat map", fwd: str = None, rev: str = None) -> str:
    """Binned match DENSITY. The dot plot answers "is there a repeat"; the heat map answers "how much of
    the locus takes part", which a binary mark cannot show once bins saturate."""
    T = _dot_theme(theme, for_export)
    fwd, rev = fwd or _DOT_FWD, rev or _DOT_REV
    peak = m.get("%s_max" % which) or 1
    kind = "direct (forward)" if which == "forward" else "inverted (reverse-complement)"
    sub = "%s match density · exact %d-mer · peak %d per bin · %s bp" % (
        kind, m["k"], peak, _fmt_int(m["length"]))
    s, ML, MT, plot, H = _dot_frame(m, W, T, title, sub)
    base = fwd if which == "forward" else rev
    mat = m[which]
    cell = plot / m["bins"]
    for i, row in enumerate(mat):
        y = MT + i * cell
        for j, c in enumerate(row):
            if not c:
                continue
            s += '<rect x="%.2f" y="%.2f" width="%.2f" height="%.2f" fill="%s"/>' % (
                ML + j * cell, y, cell + 0.4, cell + 0.4, _heat_ramp(c / peak, base))
    s += _dot_guides(m, ML, MT, plot, T, guides)
    s += _dot_legend(ML, MT + plot + 28, T, heat=(base, peak))
    if for_export:
        s += '<text x="%.1f" y="%.1f" fill="%s" font-size="7.5">%s</text>' % (
            ML, MT + plot + 44, T["faint"], esc(_DOT_FOOTER))
    return s + "</svg>"
