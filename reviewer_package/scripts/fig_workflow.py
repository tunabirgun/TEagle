"""Figure 1 - the workflow schematic.

Drawn rather than screenshotted, so it stays legible at column width and regenerates with the rest of the
figure set. Follows the convention TE-Seq (Mobile DNA 2025, 16:44) and TEtrimmer (Nat Commun 2025, 16:8429)
both use independently: name the third-party tool behind every step on the diagram, because a reader
deciding whether to trust a domain call needs to know it is HMMER's.

The one architectural claim the figure has to make legible is that the classifier sits downstream of the
structural and domain layers ONLY. It cannot see the family-naming layer, which is why a database name can
never inflate a call's confidence and why a failed lookup degrades the record rather than the call. That
point is carried by a note rather than by a crossed-out arrow: the arrow has to cross the whole diagram to
reach anything, and every routing of it collided with a box.

Layout is a named grid. Every box is placed from a column and row constant, and none is positioned by eye,
so resizing one cannot silently push it into a neighbour.

    python benchmarks/fig_workflow.py
"""
from __future__ import annotations
import os, sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                            # noqa: E402
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [os.path.join(ROOT, "app", "backend")]
from teagle_core import __version__ as TEAGLE_VERSION                      # noqa: E402
from teagle_core import domains                                            # noqa: E402

OUT = os.path.join(ROOT, "benchmarks", "figures")
MM = 1 / 25.4
DOUBLE = 183 * MM

OI = {"orange": "#E69F00", "skyblue": "#56B4E9", "green": "#009E73",
      "blue": "#0072B2", "vermillion": "#D55E00", "purple": "#CC79A7", "black": "#000000"}

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "pdf.fonttype": 42, "svg.fonttype": "none",
    "savefig.bbox": "tight", "savefig.pad_inches": 0.03,
})

N_PROFILES = len(domains.DOMAIN_INFO)      # derived, so the panel cannot drift from the bundled panel

# ---- grid -------------------------------------------------------------------------------------------
CL, CLW = 0.045, 0.150          # column A: input, then the assay boxes
CM, CMW = 0.245, 0.290          # column B: evidence layers and the classifier
CR, CRW = 0.575, 0.160          # column C: provenance seal
CO, COW = 0.775, 0.205          # column D: optional WSL2 components
BH = 0.115                      # standard box height

R_TOP = 0.795                   # input and structural detectors share a row
R_DOM = 0.635
R_CLASS = 0.450
R_ASSAY = 0.262
R_EXPORT = 0.100
AW = 0.150                      # assay box width


def box(ax, x, y, w, h, title, sub=None, face="#FFFFFF", edge=OI["black"], lw=0.9, fs=6.6, tool=None):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.006,rounding_size=0.010",
                                linewidth=lw, edgecolor=edge, facecolor=face, zorder=3))
    lines = [(title, fs, "bold", OI["black"], "normal")]
    lines += [(b, fs - 1.0, "normal", "#333333", "normal") for b in (sub or "").split("\n") if b]
    if tool:
        lines.append((tool, fs - 1.4, "normal", OI["blue"], "italic"))
    step = 0.026
    top = y + h / 2 + (len(lines) - 1) * step / 2
    for i, (txt, size, weight, col, sty) in enumerate(lines):
        ax.text(x + w / 2, top - i * step, txt, ha="center", va="center", fontsize=size,
                fontweight=weight, color=col, style=sty, zorder=4)


def arrow(ax, p, q, colour=OI["black"], lw=0.9, ls="-"):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=7, linewidth=lw,
                                 color=colour, linestyle=ls, shrinkA=2, shrinkB=2, zorder=2))


def main():
    fig, ax = plt.subplots(figsize=(DOUBLE, DOUBLE * 0.62))
    ax.set_xlim(0, 1); ax.set_ylim(0.02, 1.0); ax.axis("off")

    # bundled / optional boundary -------------------------------------------------------------------
    ax.add_patch(Rectangle((0.020, 0.070), 0.720, 0.855, fill=False, linestyle=(0, (5, 3)),
                           linewidth=0.9, edgecolor=OI["green"], zorder=1))
    ax.text(0.024, 0.935, "Bundled — runs offline, permissively licensed", fontsize=6.0,
            color=OI["green"], fontweight="bold", ha="left", va="bottom")
    ax.text(CO, 0.935, "Optional, installed WSL2 environment", fontsize=6.0,
            color=OI["vermillion"], fontweight="bold", ha="left", va="bottom")

    # input and evidence layers ---------------------------------------------------------------------
    box(ax, CL, R_TOP, CLW, BH, "One sequence", sub="accession · file · pasted\nup to tens of kb",
        face="#F2F7FA")
    box(ax, CM, R_TOP, CMW, BH, "Structural detectors",
        sub="LTR · TIR · TSD\nPBS · PPT · poly-A", edge=OI["blue"])
    box(ax, CM, R_DOM, CMW, BH, "Protein-domain panel",
        sub=f"{N_PROFILES} Pfam-A profiles\n12 longest ORFs", tool="HMMER via pyhmmer", edge=OI["blue"])
    arrow(ax, (CL + CLW, R_TOP + BH * 0.68), (CM, R_TOP + BH * 0.5), colour="#888888", lw=0.8)
    arrow(ax, (CL + CLW, R_TOP + BH * 0.30), (CM, R_DOM + BH * 0.72), colour="#888888", lw=0.8)

    # optional layers, outside the boundary ----------------------------------------------------------
    box(ax, CO, R_TOP, COW, BH, "Family naming", sub="Dfam library lookup",
        tool="RepeatMasker", edge=OI["vermillion"], face="#FFF8F4")
    box(ax, CO, R_DOM, COW, BH, "Splice alignment", sub="exons and introns\nfrom a transcript",
        tool="minimap2", edge=OI["vermillion"], face="#FFF8F4")
    arrow(ax, (CM + CMW, R_TOP + BH * 0.75), (CO, R_TOP + BH * 0.75),
          colour="#BBBBBB", lw=0.8, ls=(0, (3, 2)))

    # classifier and seal ----------------------------------------------------------------------------
    # No italic line on this box or on the seal: the italic slot names the third-party tool behind a
    # step, and both steps are TEagle's own. What the classifier does when the evidence will not support
    # a call is explanation, and explanation belongs in the caption rather than inside the image.
    box(ax, CM, R_CLASS, CMW, BH + 0.012, "Classifier", sub="Wicker class · order · superfamily",
        edge=OI["blue"], lw=1.6, face="#EAF3F8")
    # BOTH evidence layers feed the classifier. Only the domain arrow was drawn, so the figure showed one
    # input where the caption and the prose claim two. The structural arrow is routed down the left of the
    # column so it does not overlap the domain box.
    arrow(ax, (CM + CMW / 2, R_DOM), (CM + CMW / 2, R_CLASS + BH + 0.012), colour=OI["blue"], lw=1.3)
    ax.plot([CM - 0.022, CM - 0.022], [R_TOP + BH * 0.35, R_CLASS + BH * 0.55],
            color=OI["blue"], lw=1.1, zorder=2)
    ax.plot([CM, CM - 0.022], [R_TOP + BH * 0.35, R_TOP + BH * 0.35], color=OI["blue"], lw=1.1, zorder=2)
    arrow(ax, (CM - 0.022, R_CLASS + BH * 0.55), (CM, R_CLASS + BH * 0.55), colour=OI["blue"], lw=1.1)

    box(ax, CR, R_CLASS, CRW, BH + 0.012, "Provenance seal",
        sub="versions · checksums\nthresholds · input hash",
        edge=OI["green"], lw=1.5, face="#EDF7F3")
    mid = R_CLASS + (BH + 0.012) / 2
    arrow(ax, (CM + CMW, mid), (CR, mid), colour=OI["green"], lw=1.0)

    # assay layer ------------------------------------------------------------------------------------
    box(ax, CL, R_ASSAY, AW, BH, "Primer design", sub="QC on every pair", tool="Primer3",
        edge=OI["purple"], face="#FBF4F8")
    box(ax, CL + AW + 0.030, R_ASSAY, AW, BH, "In-silico PCR",
        sub="predicted amplicons\nto-scale gel", edge=OI["purple"], face="#FBF4F8")
    box(ax, CO, R_ASSAY, COW, BH, "Off-target scan", sub="against a downloaded\nassembly",
        tool="isPcr", edge=OI["vermillion"], face="#FFF8F4")
    arrow(ax, (CM + 0.015, R_CLASS), (CL + AW * 0.80, R_ASSAY + BH), colour=OI["purple"], lw=1.0)
    arrow(ax, (CL + AW, R_ASSAY + BH / 2), (CL + AW + 0.030, R_ASSAY + BH / 2),
          colour=OI["purple"], lw=1.0)
    arrow(ax, (CL + 2 * AW + 0.030, R_ASSAY + BH / 2), (CO, R_ASSAY + BH / 2),
          colour="#BBBBBB", lw=0.9, ls=(0, (3, 2)))
    # The assay run is sealed too, so this arrow points INTO the seal. Drawn from the in-silico PCR box up
    # the clear channel to the seal's lower edge.
    arrow(ax, (CL + 2 * AW + 0.045, R_ASSAY + BH * 0.75), (CR + 0.020, R_CLASS),
          colour=OI["green"], lw=0.8, ls=(0, (2, 2)))

    # export -----------------------------------------------------------------------------------------
    box(ax, CL, R_EXPORT, 2 * AW + 0.030, 0.082, "Export",
        sub="GFF3 · BED · FASTA · XLSX/CSV/TSV\nJSON manifest · SVG/PNG/PDF",
        face="#F5F5F5", edge="#999999", fs=6.2)
    arrow(ax, (CL + AW * 0.5, R_ASSAY), (CL + AW * 0.5, R_EXPORT + 0.082), colour="#999999", lw=0.8)


    os.makedirs(OUT, exist_ok=True)
    for ext in ("pdf", "svg", "png"):
        fig.savefig(os.path.join(OUT, f"fig_workflow.{ext}"), dpi=600 if ext == "png" else None)
    plt.close(fig)
    print(f"wrote fig_workflow.pdf / .svg / .png  ({N_PROFILES} profiles, TEagle {TEAGLE_VERSION})")


if __name__ == "__main__":
    main()
