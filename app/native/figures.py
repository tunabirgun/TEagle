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
# exon_derived = a lighter tint of the annotated-exon green: same family, reads "inferred, not annotated";
# gap kept distinct (light) from the darker slate flank so filler regions aren't a single ambiguous colour.
GENECOL = {"exon": "#009E73", "exon_derived": "#7fd3b8", "intron": "#8792a0",
           "cds": "#D55E00", "flank": "#5b6b7a", "gap": "#c3ccd6"}
# retroviral transcript architecture (ERV): env exons green, the removed gag-pro-pol span amber-brown so the
# "single large intron = frameshift-fused polyprotein" reads distinctly from a grey host intron.
ARCHCOL = {"exon": "#009E73", "intron": "#B0752E"}
# LTR cis-elements: PBS (leader, purple) and PPT (before 3' LTR, blue) — each distinct from the LTR blocks,
# the env-exon green and the intron amber.
CISCOL = {"PBS": "#8459C4", "PPT": "#2C7FB8"}


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


# ================= genome viewer =================
def gv_nice_step(span: float, ticks: int) -> float:
    raw = max(1.0, span / ticks)
    mag = 10 ** math.floor(math.log10(raw))
    n = raw / mag
    return (1 if n < 1.5 else 2 if n < 3.5 else 5 if n < 7.5 else 10) * mag


def gv_tracks_from_rec(rec: dict) -> dict:
    """Build a genome-viewer model from an analysis record (structural + domains + ORFs + ERV architecture)."""
    tracks, reps, cis = [], [], []
    for e in rec.get("structural", []):
        t = e["type"]
        if t.startswith("LTR") or t.startswith("TIR"):
            col = OK["LTR"] if t.startswith("LTR") else OK["TIR"]
            for p in (e.get("five_prime"), e.get("three_prime")):
                if p:
                    reps.append({"start": p[0], "end": p[1], "color": col, "label": t.split(" ")[0],
                                 "tip": f"{t} {p[0]}–{p[1]}"})
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
            cis.append({"start": p[0], "end": p[1], "color": CISCOL[key], "label": lab, "tip": tip})
        elif e.get("pos"):
            p = e["pos"]
            reps.append({"start": p[0], "end": p[1], "color": OK["tail"], "label": t.split(" ")[0],
                         "tip": f"{t} {p[0]}–{p[1]}"})
    if reps:
        tracks.append({"name": "terminal repeats", "height": 20, "features": reps})
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
# "site" = a NEUTRAL colour for a whole-genome scan with no design locus: the products are neither on- nor
# off-target, just genomic priming sites, so they must not read as the off-target warning colour.
GELPAL = {
    "transparent": {"paper": "none", "gel": "#0f1316", "well": "#04060a", "stroke": "#2a3138", "ink": "#5a656f",
                    "on": OK["on"], "off": OK["off"], "single": "#0072B2", "site": "#56B4E9", "ladder": OK["ladder"], "glow": 1.4, "band": 2.6},
    "dark":        {"paper": "#0b0e11", "gel": "#0f1316", "well": "#04060a", "stroke": "#232a30", "ink": "#8792a0",
                    "on": OK["on"], "off": OK["off"], "single": "#0072B2", "site": "#56B4E9", "ladder": OK["ladder"], "glow": 1.4, "band": 2.6},
    "white":       {"paper": "#ffffff", "gel": "#ededed", "well": "#c4c4c4", "stroke": "#cccccc", "ink": "#555555",
                    "on": "#151515", "off": "#992222", "single": "#1f5fa8", "site": "#2a6f97", "ladder": "#9a9a9a", "glow": 0.3, "band": 2.6},
    "uv":          {"paper": "#050310", "gel": "#0a0714", "well": "#000000", "stroke": "#1c1236", "ink": "#9fb4d8",
                    "on": "#5bff6b", "off": "#ffcf47", "single": "#6fb2ff", "site": "#79d0ff", "ladder": "#79d0ff", "glow": 3.2, "band": 3.1},
    "mono":        {"paper": "#0d0d0d", "gel": "#181818", "well": "#000000", "stroke": "#2b2b2b", "ink": "#b2b2b2",
                    "on": "#f2f2f2", "off": "#9a9a9a", "single": "#6f6f6f", "site": "#c0c0c0", "ladder": "#cfcfcf", "glow": 2.0, "band": 2.9},
}


def _band_opacity(total_mm: int) -> float:
    """Priming-efficiency proxy: a perfect match reads bright; each mismatch dims the band.
    All reported mismatches are 5'-proximal (the strict-3' rule forbids 3'-end mismatches)."""
    return round(max(0.4, 1.0 - 0.22 * max(0, total_mm)), 3)


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
        n_single = sum(1 for a in g if a.get("single_primer"))
        if not has_locus:                                 # neutral priming sites (single-primer artefacts still distinct)
            n_site = len(g) - n_single
            color = P["single"] if (n_single and not n_site) else P["site"]
            parts = (([f"{n_site} priming site" + ("s" if n_site != 1 else "")] if n_site else [])
                     + ([f"{n_single} single-primer"] if n_single else []))
            on = False
        else:
            n_on = sum(1 for a in g if a.get("on_target"))
            n_off = len(g) - n_on - n_single
            if n_off:                                     # any off-target co-migrating here -> flag the whole band off-target
                color = P["off"]
            elif n_single:
                color = P["single"]
            elif n_on:                                    # purely on-target at this size
                color = P["on"]
            else:
                color = P["off"]
            parts = (([f"{n_on} on-target"] if n_on else []) + ([f"{n_off} off-target"] if n_off else [])
                     + ([f"{n_single} single-primer"] if n_single else []))
            on = bool(n_on)
        min_mm = min((a.get("fwd_mm", 0) + a.get("rev_mm", 0)) for a in g)
        src = esc(g[0].get("source", ""))
        bands.append({"size": size, "color": color, "opacity": _band_opacity(min_mm),
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
                n_single = sum(1 for a in g if a.get("single_primer"))
                if l.get("has_locus", True):              # with a design locus -> on/off-target; else neutral priming sites
                    n_on = sum(1 for a in g if a.get("on_target"))
                    n_off = len(g) - n_on - n_single
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
        s += (f'<circle cx="{x0}" cy="{ly}" r="3" fill="{P["site"]}"/>'
              f'<text x="{x0+7}" y="{ly+3}" fill="{P["ink"]}" font-size="8">priming site</text>')
        xnext = x0 + 78
    else:                                    # any locus (incl. local PCR, which shares this gel) -> on/off swatches
        s += (f'<circle cx="{x0}" cy="{ly}" r="3" fill="{P["on"]}"/>'
              f'<text x="{x0+7}" y="{ly+3}" fill="{P["ink"]}" font-size="8">on-target</text>'
              f'<circle cx="{x0+64}" cy="{ly}" r="3" fill="{P["off"]}"/>'
              f'<text x="{x0+71}" y="{ly+3}" fill="{P["ink"]}" font-size="8">off-target</text>')
        xnext = x0 + 132
    if any_single:
        s += (f'<circle cx="{xnext}" cy="{ly}" r="3" fill="{P["single"]}"/>'
              f'<text x="{xnext+7}" y="{ly+3}" fill="{P["ink"]}" font-size="8">single-primer</text>')
        xnext += 82
    s += f'<text x="{xnext}" y="{ly+3}" fill="{P["ink"]}" font-size="8">L = MW ladder (bp)</text>'
    return s + "</svg>"
