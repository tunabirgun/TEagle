"""Figure 3 - does TEagle's reported confidence track correctness?

Reads benchmarks/raw/scores.json only. Every value plotted is computed here from the scored cases; none is
typed in. Tier order on the axis is the tool's own ordering, not the ordering of the result, so a
non-monotone outcome is visible as non-monotone rather than hidden by sorting.

    python benchmarks/fig_calibration.py
"""
from __future__ import annotations
import json, math, os, sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "benchmarks", "raw", "scores.json")
OUT = os.path.join(ROOT, "benchmarks", "figures")

OI = {"orange": "#E69F00", "green": "#009E73", "blue": "#0072B2",
      "vermillion": "#D55E00", "black": "#000000"}
MM = 1 / 25.4
SINGLE = 89 * MM

# The tool's own ordering, strongest first. Fixed here so the panel cannot be silently re-sorted into
# looking monotone.
TIER_ORDER = ["High", "Moderate", "Candidate"]

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 7.5, "ytick.labelsize": 7,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "axes.spines.top": False, "axes.spines.right": False,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
    "pdf.fonttype": 42, "svg.fonttype": "none",
})


def wilson(k, n, z=1.959963985):
    if n == 0:
        return (float("nan"),) * 3
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return p, max(0.0, c - h), min(1.0, c + h)


def main():
    if not os.path.exists(SRC):
        print("  SKIP calibration - benchmarks/raw/scores.json not present")
        return 1
    data = json.load(open(SRC, encoding="utf-8"))

    # EVERY answered call counts, including the negative controls. An earlier version excluded them, and
    # the case it excluded was the single Moderate-tier false positive - a telomerase reverse transcriptase
    # called a LINE. Dropping it made Moderate read 46/46 rather than 46/47, which is precisely the
    # direction of error a calibration panel exists to expose. A call on something that is not a
    # transposable element is wrong at any tier, so it is scored wrong here.
    tally = defaultdict(lambda: [0, 0])
    for c in data["cases"]:
        if c["abstained_at_order"] or c["observed_order"] == "NO_CALL":
            continue                               # no call made, so nothing to score at any tier
        if c["expected_class"] != "NOT_TE" and c["expected_order"] == "OTHER":
            continue                               # no gradable order label; out of numerator AND denominator
        t = str(c["confidence"])
        tally[t][1] += 1
        # Negative controls ARE counted, and a call on one is an error at any tier. Excluding them would
        # drop precisely the cases a calibration panel exists to expose.
        if c["expected_class"] != "NOT_TE" and c["expected_order"] == c["observed_order"]:
            tally[t][0] += 1

    tiers = [t for t in TIER_ORDER if t in tally] + [t for t in tally if t not in TIER_ORDER]
    xs = list(range(len(tiers)))
    stats = [wilson(*reversed(tally[t])) if False else wilson(tally[t][0], tally[t][1]) for t in tiers]

    fig, ax = plt.subplots(figsize=(SINGLE * 1.30, SINGLE * 0.92))
    ax.axhline(1.0, color="#CCCCCC", lw=0.6, zorder=1)
    for x, t, (p, lo, hi) in zip(xs, tiers, stats):
        n = tally[t][1]
        colour = OI["blue"] if p >= 0.9 else OI["vermillion"]
        # Clamped at zero: at k = n the Wilson upper bound can land a floating-point ulp below 1.0 while
        # p is exactly 1.0, which makes the upper error bar a negative 1e-16 and matplotlib refuses it.
        ax.errorbar([x], [p], yerr=[[max(0.0, p - lo)], [max(0.0, hi - p)]], color=colour, marker="o", ms=5.0,
                    lw=0, elinewidth=1.0, capsize=3.0, zorder=3)
        # n only. The accuracy value is already on the y-axis, and every word of explanation belongs in
        # the manuscript caption rather than baked into the image.
        ax.annotate(f"n = {n}", xy=(x, lo), xytext=(0, -10), textcoords="offset points",
                    ha="center", va="top", fontsize=6.6, color="#333333")

    ax.set_xticks(xs)
    ax.set_xticklabels(tiers)
    ax.set_xlim(-0.55, len(tiers) - 0.45)
    ax.set_ylim(0.20, 1.10)
    ax.set_yticks([0.25, 0.5, 0.75, 1.0])
    ax.set_xlabel("Confidence tier reported by TEagle")
    ax.set_ylabel("Order-level accuracy")



    os.makedirs(OUT, exist_ok=True)
    for ext in ("pdf", "svg", "png"):
        fig.savefig(os.path.join(OUT, f"fig_confidence_calibration.{ext}"),
                    dpi=600 if ext == "png" else None)
    plt.close(fig)

    summary = {t: {"n": tally[t][1], "correct": tally[t][0],
                   "accuracy": round(stats[i][0], 4),
                   "ci95": [round(stats[i][1], 4), round(stats[i][2], 4)]}
               for i, t in enumerate(tiers)}
    json.dump(summary, open(os.path.join(OUT, "fig_calibration_summary.json"), "w"), indent=1)
    print("wrote fig_confidence_calibration.pdf / .svg / .png")
    print("  ", json.dumps(summary))
    return 0


if __name__ == "__main__":
    sys.exit(main())
