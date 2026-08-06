"""Generate manuscript figures from the RAW benchmark output. No number is typed in here.

Every value plotted is read from a file under benchmarks/raw/ that was written by executing a tool. If a
panel's raw file is absent the figure is skipped with a message rather than drawn from a placeholder.

Conventions, so the figures read as one set and survive a journal's production pipeline:
  * Okabe-Ito colour-blind-safe palette, the same hues the application uses (app/native/theme.py)
  * colour is never the only carrier - every series also differs by marker and line style
  * vector PDF and SVG for the journal, 600 dpi PNG for drafts; text stays text, never outlined
  * sized for a single (89 mm) or double (183 mm) column, so nothing is rescaled by the typesetter

    python benchmarks/figures.py
"""
from __future__ import annotations
import json, math, os, sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                    # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "benchmarks", "raw")
OUT = os.path.join(ROOT, "benchmarks", "figures")

# Okabe & Ito (2008), as used throughout the application
OI = {"orange": "#E69F00", "skyblue": "#56B4E9", "green": "#009E73",
      "blue": "#0072B2", "vermillion": "#D55E00", "purple": "#CC79A7", "black": "#000000"}

MM = 1 / 25.4
SINGLE, DOUBLE = 89 * MM, 183 * MM

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 9,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 150, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
    "pdf.fonttype": 42, "svg.fonttype": "none",      # keep text as text, not outlines
})


def wilson(k: int, n: int, z: float = 1.959963985):
    """Wilson score interval. Used rather than the normal approximation because n is small and the
    proportions sit at 0 and 1, where the normal interval is degenerate."""
    if n == 0:
        return (float("nan"),) * 3
    p = k / n
    d = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return p, max(0.0, centre - half), min(1.0, centre + half)


def save(fig, stem):
    os.makedirs(OUT, exist_ok=True)
    for ext in ("pdf", "svg", "png"):
        fig.savefig(os.path.join(OUT, f"{stem}.{ext}"), dpi=600 if ext == "png" else None)
    plt.close(fig)
    print(f"  wrote {stem}.pdf / .svg / .png")


def _pearson(x, y):
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    sx = math.sqrt(sum((a - mx) ** 2 for a in x))
    sy = math.sqrt(sum((b - my) ** 2 for b in y))
    return sum((a - mx) * (b - my) for a, b in zip(x, y)) / (sx * sy) if sx and sy else float("nan")


def _ranks(v):
    """Average ranks, ties shared - so Spearman is Pearson on these."""
    order = sorted(range(len(v)), key=lambda i: v[i])
    out = [0.0] * len(v)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
            j += 1
        shared = (i + j) / 2 + 1
        for k in range(i, j + 1):
            out[order[k]] = shared
        i = j + 1
    return out


def _spearman(x, y):
    return _pearson(_ranks(x), _ranks(y))


def fig_divergence():
    """Figure: (a) reporting route against LTR-LTR identity; (b) whether the identity the detector reports
    is exact, resolved against the length of the repeat window it actually recovered.

    Panel b exists because panel a alone invites the wrong reading. A detector that reports an identity
    slightly above the truth would look like a miscalibrated instrument. The raw data show instead that the
    reported value is exact whenever the whole repeat is recovered, and departs from the truth only where
    seeding recovers a shortened, better-conserved sub-window. That is the difference between a bias and a
    bounded-evidence statement, and it is not visible in panel a.
    """
    src = os.path.join(RAW, "sim_divergence.json")
    if not os.path.exists(src):
        print("  SKIP divergence - benchmarks/raw/sim_divergence.json not present")
        return
    d = json.load(open(src, encoding="utf-8"))
    by = {}
    for r in d["records"]:
        by.setdefault(r["target_identity"], []).append(r)
    xs = sorted(by)
    floor = d["min_ltr_identity_floor"]
    n_rep = d["replicates"]
    true_len = d["scaffold_ltr_len"]

    fig, (ax, bx) = plt.subplots(1, 2, figsize=(DOUBLE, SINGLE * 0.86),
                                 gridspec_kw={"width_ratios": [1.25, 1.0], "wspace": 0.42})

    # ---- panel a: reporting route against simulated identity -------------------------------------
    det = [wilson(sum(1 for r in by[x] if r["detected_as_ltr"]), len(by[x])) for x in xs]
    adv = [wilson(sum(1 for r in by[x] if r["advisory_reported"]), len(by[x])) for x in xs]
    either = [wilson(sum(1 for r in by[x] if r["detected_as_ltr"] or r["advisory_reported"]), len(by[x]))
              for x in xs]

    # Shaded operating region rather than a bare line - the presentation habit of TEtrimmer Fig 3C. The
    # band is a reading aid here, NOT evidence for the floor: x is a property of the input, not a tunable
    # parameter, so this is not a parameter sweep and must not be read as one.
    ax.axvspan(floor, max(xs) + 2, color=OI["blue"], alpha=0.07, lw=0, zorder=0)
    ax.axvline(floor, color=OI["black"], lw=0.7, ls=(0, (4, 2)), zorder=1)

    # The union carries no independent information outside the 79-79.5 transition: above the floor it
    # equals the accepted route, below it equals the advisory route. Drawn as a wide translucent envelope
    # underneath rather than as a third equal line, so it cannot be mistaken for a separate measurement
    # and cannot hide the two series it is composed of.
    ax.plot(xs, [s[0] for s in either], color=OI["green"], lw=3.4, alpha=0.32, solid_capstyle="round",
            zorder=2, label="Either route (envelope)")
    for series, colour, marker, style, label, z in (
            (adv, OI["orange"], "s", "--", "Sub-threshold candidate", 3),
            (det, OI["blue"], "o", "-", "Accepted as a terminal repeat", 4)):
        p = [s[0] for s in series]
        lo = [s[0] - s[1] for s in series]
        hi = [s[2] - s[0] for s in series]
        ax.errorbar(xs, p, yerr=[lo, hi], color=colour, marker=marker, ls=style, ms=3.2, lw=1.1,
                    elinewidth=0.6, capsize=1.6, label=label, zorder=z)
    h, lab = ax.get_legend_handles_labels()
    order = [2, 1, 0]                                     # accepted / advisory / envelope

    ax.set_xlabel("Simulated LTR–LTR identity (%)")
    ax.set_ylabel("Proportion of replicates reported")
    ax.set_xlim(max(xs) + 2, min(xs) - 2)                 # descending: left = young, right = old
    ax.set_ylim(-0.04, 1.10)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    # Placed OUTSIDE the axes. Every in-axes position collides with something: the accepted route sits at
    # 1.0 across the whole upper left, the advisory route sits at 0.0 across the whole lower left, and the
    # transition fills the right. An external legend cannot occlude data at any figure scale.
    ax.legend([h[i] for i in order], [lab[i] for i in order], loc="upper center",
              bbox_to_anchor=(0.5, -0.215), ncol=3, frameon=False, fontsize=6.2,
              columnspacing=1.15, handlelength=1.7, handletextpad=0.45)

    # ---- panel b: is the reported identity exact? ------------------------------------------------
    acc = [r for r in d["records"] if r["detected_as_ltr"] and r["reported_ltr_len"]]
    lens = [r["reported_ltr_len"] for r in acc]
    delta = [r["reported_identity"] - r["realised_identity"] for r in acc]
    whole = [(L, D) for L, D in zip(lens, delta) if L >= true_len]
    part = [(L, D) for L, D in zip(lens, delta) if L < true_len]

    bx.axhline(0.0, color=OI["black"], lw=0.7, ls=(0, (4, 2)), zorder=1)
    bx.scatter([L for L, _ in part], [D for _, D in part], s=9, facecolors="none",
               edgecolors=OI["vermillion"], linewidths=0.6, marker="o", zorder=3,
               label="Repeat partly recovered")
    # Every whole-recovery replicate lands on the same coordinate, so one visible marker stands for all
    # of them. The count and the coordinate go in the legend label; an arrow annotation into this corner
    # collides with the dense partly-recovered cloud.
    bx.scatter([L for L, _ in whole], [D for _, D in whole], s=15, color=OI["blue"], marker="D",
               zorder=4, label="Whole repeat recovered")

    rho = _spearman(lens, delta)
    bx.set_xlabel(f"Repeat length recovered (bp of {true_len})")
    bx.set_ylabel("Reported − true identity (percentage points)")
    bx.set_ylim(-0.26, max(delta) * 1.34)
    bx.legend(loc="upper right", frameon=False, handletextpad=0.4, fontsize=6.4, labelspacing=0.45)

    for axis, letter in ((ax, "A"), (bx, "B")):            # bold CAPS panel letters
        axis.annotate(letter, xy=(-0.185, 1.14), xycoords="axes fraction", ha="left", va="top",
                      fontsize=10, fontweight="bold")

    save(fig, "fig_divergence_sensitivity")

    # the numbers the manuscript text will quote, written out so the prose cannot drift from the figure
    summary = {
        "floor": floor,
        "lowest_identity_fully_accepted": min((x for x in xs if all(r["detected_as_ltr"] for r in by[x])),
                                              default=None),
        "highest_identity_never_accepted": max((x for x in xs if not any(r["detected_as_ltr"] for r in by[x])),
                                               default=None),
        "lowest_identity_any_advisory": min((x for x in xs if any(r["advisory_reported"] for r in by[x])),
                                            default=None),
        "highest_identity_no_report_at_all": max((x for x in xs if not any(
            r["detected_as_ltr"] or r["advisory_reported"] for r in by[x])), default=None),
        "replicates": n_rep,
        "total_runs": len(d["records"]),
        "scaffold": d["scaffold_accession"],
        "scaffold_ltr_len": true_len,
        # panel b - reported identity against the length of the window actually recovered
        "accepted_n": len(acc),
        "whole_repeat_recovered_n": len(whole),
        "whole_repeat_max_abs_error": max((abs(D) for _, D in whole), default=None),
        "partly_recovered_n": len(part),
        "partly_recovered_mean_error": round(sum(D for _, D in part) / len(part), 4) if part else None,
        "partly_recovered_max_error": round(max(D for _, D in part), 4) if part else None,
        "shortest_recovered_bp": min(lens),
        "shortest_recovered_pct_of_true": round(100 * min(lens) / true_len, 1),
        "spearman_rho_len_vs_error": round(rho, 4),
    }
    json.dump(summary, open(os.path.join(OUT, "fig_divergence_summary.json"), "w"), indent=1)
    print("   summary:", json.dumps(summary))


def main():
    print("generating figures from benchmarks/raw/ ...")
    fig_divergence()
    print("done ->", OUT)


if __name__ == "__main__":
    main()
