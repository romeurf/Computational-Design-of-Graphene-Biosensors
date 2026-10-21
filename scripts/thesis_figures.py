"""
Schematic figures for the dissertation that are not produced by scripts/analysis.py.

Three figures, each drawn from the committed result tables so that the numbers printed on
them cannot drift from the text:

  pipeline_funnel.png   the stages of the method and the candidate windows surviving each
                        filter, from the screening funnel
  reference_panel.png   the reference panel under a single global criterion and under the
                        criteria of each probe's own species and type
  debye_geometries.png  the three geometries in which the Debye coverage is measured, with
                        the measured phosphate counts

Usage:
    python scripts/thesis_figures.py
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle

BASE_DIR = Path(__file__).resolve().parent.parent
OUT = BASE_DIR / "thesis" / "images"
FIGDIR = BASE_DIR / "output" / "figures"

# validated single-hue ordinal ramp (light -> dark) and neutral ink
RAMP = ["#cde2fb", "#86b6ef", "#5598e7", "#256abf", "#104281"]
PALE, MID, DARK = "#cde2fb", "#86b6ef", "#256abf"
INK, INK2, MUTED, LINE = "#0b0b0b", "#52514e", "#898781", "#c3c2b7"
FONT = "DejaVu Sans"
plt.rcParams.update({"font.family": FONT, "font.size": 9,
                     "pdf.fonttype": 42, "svg.fonttype": "none"})


def save(fig, name):
    """Vector PDF for the dissertation, raster copy for quick viewing."""
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(FIGDIR / f"{name}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  {name}.pdf")


# ----------------------------------------------------------------- figure 1
def pipeline_funnel():
    """Method stages, with a bar for the windows surviving each filter."""
    stages = [
        ("Sequence retrieval (NCBI Entrez)", None),
        ("Length clustering", None),
        ("Multiple sequence alignment (MAFFT)", None),
        ("Conservation scoring (PPI)", None),
        ("Sliding-window candidates", 7661),
        ("Melting temperature", 6052),
        ("GC content", 4365),
        ("Hairpin $\\Delta G$", 3978),
        ("Homodimer $\\Delta G$", 2575),
        ("Secondary structure (seqfold)", 2439),
        ("Quality ranking: top 5 per gene", 30),
        ("3D prediction (Boltz-2)", 30),
        ("Assembly on graphene + Debye coverage", 30),
        ("Pareto front: final selection", 11),
    ]
    n = len(stages)
    fig, ax = plt.subplots(figsize=(7.2, 7.4), dpi=200)
    box_w, x0 = 0.40, 0.03
    bar_x0, bar_w = 0.47, 0.40
    top = 1.0
    h = top / n
    vmax = 7661

    for i, (label, count) in enumerate(stages):
        y = top - (i + 0.5) * h
        filt = count is not None and i >= 5
        face = PALE if count is None else ("#eef4fd" if not filt else "#ffffff")
        ax.add_patch(FancyBboxPatch((x0, y - h * 0.33), box_w, h * 0.66,
                                    boxstyle="round,pad=0.004,rounding_size=0.012",
                                    linewidth=0.8, edgecolor=LINE, facecolor=face))
        ax.text(x0 + 0.012, y, label, va="center", ha="left", fontsize=8.6, color=INK)
        if i < n - 1:
            ax.add_patch(FancyArrowPatch((x0 + box_w / 2, y - h * 0.34),
                                         (x0 + box_w / 2, y - h * 0.66),
                                         arrowstyle="-|>", mutation_scale=7,
                                         linewidth=0.8, color=MUTED))
        if count is None:
            continue
        frac = count / vmax
        shade = DARK if i >= 10 else MID if i >= 5 else PALE
        ax.add_patch(Rectangle((bar_x0, y - h * 0.20), max(bar_w * frac, 0.004), h * 0.40,
                               linewidth=0, facecolor=shade))
        ax.text(bar_x0 + max(bar_w * frac, 0.004) + 0.012, y, f"{count:,}".replace(",", " "),
                va="center", ha="left", fontsize=8.6, color=INK)

    ax.text(x0, top + 0.022, "Stage", fontsize=8.6, color=INK2, fontweight="bold")
    ax.text(bar_x0, top + 0.022, "Candidates remaining", fontsize=8.6, color=INK2, fontweight="bold")
    ax.set_xlim(0, 1.0)
    ax.set_ylim(-0.02, top + 0.055)
    ax.axis("off")
    fig.tight_layout()
    save(fig, "pipeline_funnel")
    plt.close(fig)


# ----------------------------------------------------------------- figure 2
def reference_panel():
    """Reference probes passing under a global criterion vs under their own profile."""
    rows = [("Bacteria", 67, 31, 40), ("Virus", 4, 1, 1), ("Fungus", 1, 1, 0),
            ("Protozoan", 1, 0, 1), ("Host (human)", 1, 0, 0)]
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(7.8, 3.1), dpi=200,
                                  gridspec_kw={"width_ratios": [1.5, 1]})
    y = np.arange(len(rows))[::-1]
    bh = 0.34
    for i, (_, _, g, o) in zip(y, rows):
        ax.barh(i + bh / 2 + 0.02, g, height=bh, color=MID, linewidth=0)
        ax.barh(i - bh / 2 - 0.02, o, height=bh, color=DARK, linewidth=0)
        for v, yy in ((g, i + bh / 2 + 0.02), (o, i - bh / 2 - 0.02)):
            ax.text(v + 0.7, yy, str(v), va="center", ha="left", fontsize=8.5, color=INK)
    ax.set_yticks(y, [f"{lab}\n($n$ = {n})" for lab, n, _, _ in rows], fontsize=8.5)
    ax.set_xlabel("reference probes passing the basic screen", color=INK2, fontsize=8.8)
    ax.set_xlim(0, 48)
    ax.tick_params(axis="y", length=0, colors=INK)
    ax.tick_params(axis="x", colors=INK2, labelsize=8.2)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.spines["bottom"].set_color(LINE)
    ax.legend([Rectangle((0, 0), 1, 1, color=MID), Rectangle((0, 0), 1, 1, color=DARK)],
              ["single global criterion", "own species and type"],
              loc="lower right", frameon=False, fontsize=8.2)

    moves = [("Pseudomonas aeruginosa", "+6"), ("Klebsiella pneumoniae", "+2"),
             ("Staphylococcus aureus", "+1"), ("Plasmodium falciparum", "+1"),
             ("Aspergillus fumigatus", "-1")]
    ax2.text(0, 1.02, "Probes that change outcome", fontsize=8.8, color=INK2,
             fontweight="bold", va="top")
    for k, (sp, lab) in enumerate(moves):
        yy = 0.84 - k * 0.135
        col = DARK if lab.startswith("+") else "#9a3a35"
        ax2.text(0.0, yy, lab, fontsize=9.3, color=col, fontweight="bold", va="center")
        ax2.text(0.15, yy, sp, fontsize=8.3, color=INK, va="center", style="italic")
    ax2.plot([0, 0.95], [0.24, 0.24], color=LINE, lw=0.8)
    ax2.text(0, 0.17, "Ten probes are recovered whose composition departs\n"
                      "from the balanced bacterial case; one is rejected,\n"
                      "the A. fumigatus profile being narrower than the\n"
                      "global one.",
             fontsize=7.9, color=INK2, va="top", linespacing=1.5)
    ax2.set_xlim(0, 1); ax2.set_ylim(-0.12, 1.06); ax2.axis("off")
    fig.tight_layout()
    save(fig, "reference_panel")
    plt.close(fig)


# ----------------------------------------------------------------- figure 3
def debye_geometries():
    """The three measured geometries, with the measured phosphate counts."""
    cov = pd.read_csv(BASE_DIR / "results" / "debye_coverage.csv", sep=";")
    med = {"rigid": cov["P_1x PBS"].median(),
           "anchored": cov["P_anchored_1x PBS"].median(),
           "duplex": cov["P_target_1x PBS"].median()}
    rng = np.random.default_rng(0)
    LAM = 7.6                                    # lambda_D at 1x PBS, in Angstrom

    fig, axes = plt.subplots(1, 3, figsize=(7.8, 3.5), dpi=200, sharey=True)
    titles = ["Rigid placement", "Anchored by the $3\'$ end", "Anchored duplex"]
    notes = [f"median {int(med['rigid'])} of the 18\u201326 probe\nphosphates within $\\lambda_D$",
             f"median {int(med['anchored'])} probe phosphate\nwithin $\\lambda_D$",
             f"median {int(med['duplex'])} target phosphates\nwithin $\\lambda_D$"]

    for ax, title, note, kind in zip(axes, titles, notes, ("rigid", "anchored", "duplex")):
        ax.add_patch(Rectangle((-1.0, 0), 2.0, LAM, linewidth=0, facecolor="#eaf2fd"))
        ax.plot([-1.0, 1.0], [LAM, LAM], color=INK2, lw=0.9, ls=(0, (4, 3)))
        ax.plot([-1.0, 1.0], [0, 0], color=INK, lw=2.4, solid_capstyle="butt")

        if kind == "duplex":                     # two strands twisting about one axis
            t = np.linspace(0, 1, 26)
            ph = 2 * np.pi * t * 2.0
            hp = 3.2 + t * 72                    # probe, tethered at its 3' end
            # the target runs antiparallel, so its first phosphate sits one rise higher
            ht = hp + 3.4
            xa, xb = 0.30 * np.sin(ph), 0.30 * np.sin(ph + np.pi)
            ax.plot(xa, hp, color=PALE, lw=1.1, zorder=2)
            ax.plot(xb, ht, color=MID, lw=1.1, zorder=2)
            ax.scatter(xa, hp, s=11, color=PALE, zorder=3, linewidths=0)
            ax.scatter(xb, ht, s=14, color=MID, zorder=3, linewidths=0)
            ax.text(-0.95, 68, "probe", fontsize=7.6, color=MID)
            ax.text(-0.95, 61, "target", fontsize=7.6, color=DARK)
        else:
            if kind == "rigid":                  # compact, lowest atom at the contact distance
                x = rng.normal(0, 0.26, 22)
                h = 3.2 + np.abs(rng.normal(0, 13, 22))
            else:                                # tethered at the 3' end, extending away
                t = np.linspace(0, 1, 22)
                x = 0.34 * np.sin(2 * np.pi * t * 1.4) * t
                h = 3.2 + t * 52
            inside = h <= LAM
            ax.plot(x, h, color=PALE, lw=1.0, zorder=2)
            ax.scatter(x[~inside], h[~inside], s=16, color=MID, zorder=3, linewidths=0)
            ax.scatter(x[inside], h[inside], s=26, color=DARK, zorder=4, linewidths=0)

        ax.set_title(title, fontsize=9, color=INK, pad=8)
        ax.text(0, -9, note, ha="center", va="top", fontsize=8.1, color=INK2, linespacing=1.5)
        ax.set_xlim(-1.0, 1.0)
        ax.set_ylim(-24, 84)
        ax.set_xticks([])
        for sp in ("top", "right", "bottom"):
            ax.spines[sp].set_visible(False)
        ax.spines["left"].set_color(LINE)

    axes[0].set_ylabel("height above the graphene ($\\AA$)", color=INK2, fontsize=8.6)
    axes[0].tick_params(axis="y", colors=INK2, labelsize=8)
    axes[0].annotate("$\\lambda_D$ = 7.6 $\\AA$", xy=(-0.55, LAM), xytext=(-0.9, 26),
                     fontsize=8.0, color=INK2,
                     arrowprops=dict(arrowstyle="-", color=INK2, lw=0.7))
    axes[0].text(-0.9, 20, "(1$\\times$ PBS)", fontsize=8.0, color=INK2)
    fig.legend([plt.Line2D([], [], marker="o", ls="", ms=5, color=DARK),
                plt.Line2D([], [], marker="o", ls="", ms=4.5, color=MID)],
               ["phosphate within $\\lambda_D$", "phosphate beyond it"],
               loc="lower center", ncol=2, frameon=False, fontsize=8.2,
               bbox_to_anchor=(0.5, -0.03))
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    save(fig, "debye_geometries")
    plt.close(fig)


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    FIGDIR.mkdir(parents=True, exist_ok=True)
    print("writing figures:")
    pipeline_funnel()
    reference_panel()
    debye_geometries()
