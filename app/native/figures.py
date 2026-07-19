"""Publication figure builders, ported from the validated web UI (app/web/app.js).
Pure functions returning SVG strings. The native app renders them on-screen via QSvgRenderer
and exports the identical string as SVG / rasterises it to PNG — so figures match the web output.

Two figures: an interactive genome viewer (ruler + tracks, windowed) and a to-scale
multi-lane agarose gel (log MW axis, 5 palettes). Geometry mirrors svgGenome / svgGel exactly."""
from __future__ import annotations
import math

FIGFONT = "Cascadia Mono, Consolas, monospace"   # bundled UI font; mono digits align on the gel ladder / ruler

# Okabe-Ito data palette (domain / feature hues) — must match the web figure bands
OK = {"RT": "#0072B2", "INT": "#E69F00", "RNaseH": "#009E73", "PR": "#CC79A7", "GAG": "#7A7A7A",
      "CHR": "#D55E00", "TPase": "#D55E00", "LTR": "#56B4E9", "TIR": "#E69F00", "tail": "#CC79A7",
      "ORF": "#4C6C97", "on": "#009E73", "off": "#D55E00", "ladder": "#999999"}
GENECOL = {"exon": "#009E73", "intron": "#8792a0", "cds": "#D55E00"}


def esc(s) -> str:
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def _label_ink(hex_color: str) -> str:
    """Pick black/white label by best WCAG contrast against a fill colour."""
    def lin(v):
        v /= 255.0
        return v / 12.92 if v <= 0.03928 else ((v + 0.055) / 1.055) ** 2.4
    r = int(hex_color[1:3], 16); g = int(hex_color[3:5], 16); b = int(hex_color[5:7], 16)
    L = 0.2126 * lin(r) + 0.7152 * lin(g) + 0.0722 * lin(b)
    return "#fff" if (1.05 / (L + 0.05)) >= ((L + 0.05) / 0.05) else "#111"


def _fmt_int(n) -> str:
    return f"{int(round(n)):,}"


# ================= genome viewer =================
def gv_nice_step(span: float, ticks: int) -> float:
    raw = max(1.0, span / ticks)
    mag = 10 ** math.floor(math.log10(raw))
    n = raw / mag
    return (1 if n < 1.5 else 2 if n < 3.5 else 5 if n < 7.5 else 10) * mag


def gv_tracks_from_rec(rec: dict) -> dict:
    """Build a genome-viewer model from an analysis record (structural + domains + ORFs)."""
    tracks, reps = [], []
    for e in rec.get("structural", []):
        t = e["type"]
        if t.startswith("LTR") or t.startswith("TIR"):
            col = OK["LTR"] if t.startswith("LTR") else OK["TIR"]
            for p in (e.get("five_prime"), e.get("three_prime")):
                if p:
                    reps.append({"start": p[0], "end": p[1], "color": col, "label": t.split(" ")[0],
                                 "tip": f"{t} {p[0]}–{p[1]}"})
        elif e.get("pos"):
            p = e["pos"]
            reps.append({"start": p[0], "end": p[1], "color": OK["tail"], "label": t.split(" ")[0],
                         "tip": f"{t} {p[0]}–{p[1]}"})
    if reps:
        tracks.append({"name": "terminal repeats", "height": 20, "features": reps})
    doms = [{"start": d["nt"][0], "end": d["nt"][1], "color": OK.get(d["domain"], "#888"),
             "label": d["domain"],
             "tip": f"{d['domain']} · {d.get('label','')} · nt {d['nt'][0]}–{d['nt'][1]} · score {d.get('score')}"}
            for d in rec.get("domains", [])]
    if doms:
        tracks.append({"name": "protein domains", "height": 22, "features": doms})
    orfs = [{"start": o["start"], "end": o["end"], "color": OK["ORF"], "strand": o["strand"],
             "tip": f"ORF {o['strand']}{o['frame']} · {o['length_aa']} aa"} for o in rec.get("orfs", [])]
    if orfs:
        tracks.append({"name": "ORFs (± strand)", "height": 26, "features": orfs, "stranded": True})
    return {"length": rec.get("composition", {}).get("length", 1) or 1, "tracks": tracks}


def gv_tracks_from_gene(gm: dict, length: int) -> dict:
    tracks = []
    feat = [{"start": e["start"], "end": e["end"], "color": GENECOL["exon"], "label": "exon",
             "tip": f"exon {e['start']}–{e['end']} ({e['end']-e['start']} bp)"} for e in gm.get("exons", [])]
    feat += [{"start": i["start"], "end": i["end"], "color": GENECOL["intron"], "intron": True,
              "tip": f"intron {i['start']}–{i['end']}" + (
                  f" · {i['donor']}…{i['acceptor']}{' (canonical)' if i.get('canonical') else ''}"
                  if i.get("donor") else "")} for i in gm.get("introns", [])]
    if feat:
        tracks.append({"name": "exons / introns", "height": 22, "features": feat})
    cds = [{"start": c["start"], "end": c["end"], "color": GENECOL["cds"], "label": "CDS",
            "tip": f"CDS {c['start']}–{c['end']}"} for c in gm.get("cds", [])]
    if cds:
        tracks.append({"name": "CDS (coding)", "height": 16, "features": cds})
    return {"length": length or 1, "tracks": tracks}


def _gv_theme(theme: str, for_export: bool) -> dict:
    if for_export:
        return {"paper": "none", "ink": "#222", "faint": "#555", "grid": "#dcdfe3",
                "track": "#00000000", "lane": "#0000000d", "frame": "#c7ccd2", "win": "#1f6feb"}
    if theme == "white":
        return {"paper": "#ffffff", "ink": "#141b21", "faint": "#5a6570", "grid": "#eceef1",
                "track": "#f6f8fa", "lane": "#eef1f4", "frame": "#dde1e6", "win": "#1f6feb"}
    return {"paper": "#0b1016", "ink": "#e6edf1", "faint": "#8a959d", "grid": "#182029",
            "track": "#10171e", "lane": "#121b23", "frame": "#243039", "win": "#4aa8ff"}


def svg_genome(model: dict, view: dict, W: float, theme: str, for_export: bool = False, return_regions: bool = False):
    """Render the genome viewer SVG. With return_regions=True, also return a list of on-screen
    feature hit-boxes (SVG coords + feature identity) so the interactive canvas can hover/right-click."""
    L = model["length"]
    regions = []
    ML, MR, ovH, rulerH, MT = 96, 16, 13, 24, 34
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
        s += f'<text x="{ML-10}" y="{ty+th/2+3:.1f}" fill="{T["faint"]}" font-size="9.5" text-anchor="end">{esc(t["name"])}</text>'
        s += f'<rect x="{ML}" y="{ty}" width="{plotW}" height="{th}" rx="3" fill="{T["track"]}"/>'
        stranded = t.get("stranded")
        if stranded:
            s += f'<line x1="{ML}" y1="{ty+th/2:.1f}" x2="{ML+plotW}" y2="{ty+th/2:.1f}" stroke="{T["grid"]}"/>'
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
            if f.get("label") and w > 26:
                s += (f'<text x="{a+4:.1f}" y="{yy+max(hh,3)-3:.1f}" fill="{_label_ink(f["color"])}" '
                      f'font-size="9" font-weight="700" pointer-events="none">{esc(f["label"])}</text>')
    svg = s + "</svg>"
    return (svg, regions) if return_regions else svg


# ================= agarose gel =================
GELPAL = {
    "transparent": {"paper": "none", "gel": "#0f1316", "well": "#04060a", "stroke": "#2a3138", "ink": "#5a656f",
                    "on": OK["on"], "off": OK["off"], "ladder": OK["ladder"], "glow": 1.4, "band": 2.6},
    "dark":        {"paper": "#0b0e11", "gel": "#0f1316", "well": "#04060a", "stroke": "#232a30", "ink": "#8792a0",
                    "on": OK["on"], "off": OK["off"], "ladder": OK["ladder"], "glow": 1.4, "band": 2.6},
    "white":       {"paper": "#ffffff", "gel": "#ededed", "well": "#c4c4c4", "stroke": "#cccccc", "ink": "#555555",
                    "on": "#151515", "off": "#992222", "ladder": "#9a9a9a", "glow": 0.3, "band": 2.6},
    "uv":          {"paper": "#050310", "gel": "#0a0714", "well": "#000000", "stroke": "#1c1236", "ink": "#9fb4d8",
                    "on": "#5bff6b", "off": "#ffcf47", "ladder": "#79d0ff", "glow": 3.2, "band": 3.1},
    "mono":        {"paper": "#0d0d0d", "gel": "#181818", "well": "#000000", "stroke": "#2b2b2b", "ink": "#b2b2b2",
                    "on": "#f2f2f2", "off": "#9a9a9a", "ladder": "#cfcfcf", "glow": 2.0, "band": 2.9},
}


def _gel_geometry(data: dict):
    """Shared layout math for svg_gel and gel_regions (must stay identical so hit-boxes line up)."""
    lanes = data.get("lanes") or [{"label": data.get("laneLabel", "PCR"), "amplicons": data.get("amplicons", [])}]
    sizes = [a["length"] for l in lanes for a in (l.get("amplicons") or [])]
    smallest = min(sizes) if sizes else 90
    minbp = max(25, min(90, smallest - 10))
    maxbp = max([1600] + sizes)
    laneW, gap, x0, top, botPad, H = 40, 12, 62, 48, 46, 366
    bot = H - botPad

    def y(bp):
        return top + (math.log(maxbp) - math.log(max(bp, minbp))) / (math.log(maxbp) - math.log(minbp)) * (bot - top)

    def laneX(i):
        return x0 + i * (laneW + gap)
    return lanes, minbp, maxbp, laneW, x0, top, bot, H, y, laneX


def gel_regions(data: dict):
    """Per-band hit-boxes (SVG coords) + the amplicon each band represents, for hover/right-click
    on the gel figure. Lane 0 is the ladder (skipped); sample lanes start at column 1."""
    lanes, minbp, maxbp, laneW, x0, top, bot, H, y, laneX = _gel_geometry(data)
    P = GELPAL["dark"]
    regions = []
    for i, l in enumerate(lanes):
        lx = laneX(i + 1)
        for a in (l.get("amplicons") or []):
            yy = y(a["length"]); h = P["band"]
            call = "on-target" if a.get("on_target") else "off-target"
            tip = f'{a["length"]} bp · {call}' + (f' · {a.get("source","")}' if a.get("source") else "")
            regions.append({"x0": lx + 3, "y0": yy - 4, "x1": lx + laneW - 3, "y1": yy + 4,
                            "tip": tip, "amplicon": a, "pair": l.get("label", "")})
    return regions


def svg_gel(data: dict, bg: str) -> str:
    """data: {lanes:[{label, amplicons:[{length,on_target,source}]}]} or legacy {amplicons}."""
    lanes = data.get("lanes") or [{"label": data.get("laneLabel", "PCR"), "amplicons": data.get("amplicons", [])}]
    sizes = [a["length"] for l in lanes for a in (l.get("amplicons") or [])]
    smallest = min(sizes) if sizes else 90
    minbp = max(25, min(90, smallest - 10))
    maxbp = max([1600] + sizes)
    LADDER = [m for m in (1500, 1000, 700, 500, 400, 300, 200, 100, 50) if minbp <= m <= maxbp]
    laneW, gap, x0, top, botPad, H = 40, 12, 62, 48, 46, 366
    bot = H - botPad
    def y(bp):
        return top + (math.log(maxbp) - math.log(max(bp, minbp))) / (math.log(maxbp) - math.log(minbp)) * (bot - top)
    cols = 1 + len(lanes)
    W = max(x0 + cols * (laneW + gap) + 12, 300)
    P = GELPAL.get(bg, GELPAL["transparent"])
    def laneX(i): return x0 + i * (laneW + gap)
    s = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" '
         f'font-family="{FIGFONT}">')
    s += (f'<defs><filter id="glow" x="-40%" y="-140%" width="180%" height="380%">'
          f'<feGaussianBlur stdDeviation="{P["glow"]}" result="b"/><feMerge>'
          f'<feMergeNode in="b"/><feMergeNode in="SourceGraphic"/></feMerge></filter></defs>')
    if P["paper"] != "none":
        s += f'<rect width="{W}" height="{H}" fill="{P["paper"]}"/>'
    s += (f'<rect x="{x0-7:.1f}" y="{top-16}" width="{cols*(laneW+gap)+2:.1f}" height="{bot-top+30:.1f}" '
          f'rx="2" fill="{P["gel"]}" stroke="{P["stroke"]}"/>')
    s += f'<text x="{x0-14}" y="{top-19}" fill="{P["ink"]}" font-size="8" text-anchor="end">bp</text>'
    for m in LADDER:
        yy = y(m)
        s += f'<text x="{x0-14}" y="{yy+2.5:.1f}" fill="{P["ink"]}" font-size="8" text-anchor="end">{m}</text>'

    def draw_lane(col, label, bands, is_ladder):
        nonlocal s
        lx = laneX(col)
        s += f'<rect x="{lx+4:.1f}" y="{top-13}" width="{laneW-8}" height="4" rx="1" fill="{P["well"]}"/>'
        s += f'<text x="{lx+laneW/2:.1f}" y="{top-18}" fill="{P["ink"]}" font-size="8.5" text-anchor="middle">{esc(label)}</text>'
        for bd in (bands or []):
            yy = y(bd["size"]); h = 1.6 if is_ladder else P["band"]
            filt = "" if is_ladder else ' filter="url(#glow)"'
            title = f'{bd["size"]} bp' + (f' · {bd["t"]}' if bd.get("t") else "")
            s += (f'<rect x="{lx+3:.1f}" y="{yy-h/2:.1f}" width="{laneW-6}" height="{h}" rx="1" '
                  f'fill="{bd["color"]}"{filt}><title>{title}</title></rect>')
        if not is_ladder and not (bands or []):
            s += f'<text x="{lx+laneW/2:.1f}" y="{bot+13}" fill="{P["ink"]}" font-size="7" text-anchor="middle">—</text>'

    draw_lane(0, "L", [{"size": m, "color": P["ladder"]} for m in LADDER], True)
    for i, l in enumerate(lanes):
        bands = [{"size": a["length"], "color": P["on"] if a.get("on_target") else P["off"],
                  "t": ("on-target" if a.get("on_target") else "off-target") + " · " + esc(a.get("source", ""))}
                 for a in (l.get("amplicons") or [])]
        draw_lane(i + 1, l["label"], bands, False)
    ly = H - 12
    s += (f'<circle cx="{x0}" cy="{ly}" r="3" fill="{P["on"]}"/>'
          f'<text x="{x0+7}" y="{ly+3}" fill="{P["ink"]}" font-size="8">on-target</text>'
          f'<circle cx="{x0+64}" cy="{ly}" r="3" fill="{P["off"]}"/>'
          f'<text x="{x0+71}" y="{ly+3}" fill="{P["ink"]}" font-size="8">off-target</text>'
          f'<text x="{x0+140}" y="{ly+3}" fill="{P["ink"]}" font-size="8">L = MW ladder (bp)</text>')
    return s + "</svg>"
