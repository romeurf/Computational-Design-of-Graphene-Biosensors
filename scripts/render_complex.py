"""
Render a probe-graphene complex PDB as a two-panel figure (side and top view).

Draws the actual atomic coordinates of the assembled model: the graphene sheet as its
carbon lattice and the probe coloured by element, with the probe-sheet separation
annotated. Used to produce the figures in the dissertation.

Usage:
    python scripts/render_complex.py <complex.pdb> <out.png> [title]
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser

pdb, out = sys.argv[1], sys.argv[2]
title = sys.argv[3] if len(sys.argv) > 3 else ""

model = PDBParser(QUIET=True).get_structure("c", pdb)[0]
G = np.array([a.coord for a in model["G"].get_atoms()])
patoms = list(model["D"].get_atoms())
P = np.array([a.coord for a in patoms])
elem = [a.element.strip().upper() for a in patoms]

# CPK-like colours for the probe; graphene stays neutral grey
COL = {"C": "#4d4d4d", "N": "#2b6cb0", "O": "#c0392b", "P": "#d68910"}
pc = [COL.get(e, "#7f8c8d") for e in elem]
gap = P[:, 2].min() - G[:, 2].mean()

fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))

# --- side view (x-z): shows the probe standing above the plane -------------
ax[0].scatter(G[:, 0], G[:, 2], s=4, c="#9aa0a6", linewidths=0)
ax[0].scatter(P[:, 0], P[:, 2], s=14, c=pc, linewidths=0)
ax[0].axhline(G[:, 2].mean(), color="#9aa0a6", lw=0.8, zorder=0)
ax[0].annotate("", xy=(G[:, 0].min() + 6, G[:, 2].mean()), xytext=(G[:, 0].min() + 6, P[:, 2].min()),
               arrowprops=dict(arrowstyle="<->", color="#c0392b", lw=1.2))
ax[0].text(G[:, 0].min() + 8, (P[:, 2].min()) / 2, f"{gap:.1f} Å", color="#c0392b", fontsize=9, va="center")
ax[0].set(xlabel="x (Å)", ylabel="z (Å)", title="Side view")

# legend: graphene + the elements actually present in the probe
from matplotlib.lines import Line2D
present = [e for e in ("C", "N", "O", "P") if e in set(elem)]
handles = [Line2D([], [], marker="o", ls="", ms=5, color="#9aa0a6", label="graphene (C)")] + \
          [Line2D([], [], marker="o", ls="", ms=5, color=COL[e], label=f"probe {e}") for e in present]
ax[0].legend(handles=handles, loc="upper right", fontsize=7.5, frameon=False, ncol=2)

# --- top view (x-y): shows the probe footprint inside the sheet ------------
ax[1].scatter(G[:, 0], G[:, 1], s=3, c="#d5d8dc", linewidths=0)
ax[1].scatter(P[:, 0], P[:, 1], s=14, c=pc, linewidths=0)
ax[1].set(xlabel="x (Å)", ylabel="y (Å)", title="Top view")
ax[1].set_aspect("equal")

for a in ax:
    a.spines[["top", "right"]].set_visible(False)
    a.tick_params(labelsize=8)
if title:
    fig.suptitle(title, fontsize=11)
fig.tight_layout()
fig.savefig(out, dpi=200)
print(f"wrote {out}  (probe {len(P)} atoms, sheet {len(G)} atoms, gap {gap:.2f} A)")
