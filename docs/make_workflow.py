#!/usr/bin/env python3
"""Draw the toolkit workflow diagram -> docs/workflow.svg and docs/workflow.png.
Reproducible source for the figure shown in the README."""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

IN, MOD, OUT = "#eef2f7", "#2c7fb8", "#fdf0e6"
INE, MODE, OUTE = "#9aa7b4", "#1c5a84", "#d9995f"
EDGE = "#6b7a88"
fig, ax = plt.subplots(figsize=(9.6, 5.0))
ax.set_xlim(0, 10); ax.set_ylim(0, 5); ax.axis("off")

def box(x, y, w, h, text, fc, ec, tc, bold=False, fs=10.5):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.08",
                                fc=fc, ec=ec, lw=1.4, mutation_aspect=1))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", color=tc,
            fontsize=fs, fontweight="bold" if bold else "normal", zorder=5)
    return (x, y, w, h)

def arrow(a, b, rad=0.0, label="", lx=0, ly=0):
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=12, lw=1.4,
                                 color=EDGE, connectionstyle="arc3,rad=%s" % rad, zorder=1))
    if label:
        ax.text(lx, ly, label, fontsize=8.5, style="italic", color=EDGE, ha="center", va="center")

# column headers
for x, t in ((1.45, "INPUTS"), (4.95, "MODULES  (config-selectable)"), (8.45, "OUTPUTS  (CSV)")):
    ax.text(x, 4.86, t, ha="center", color="#8a97a4", fontsize=10, fontweight="bold")

W, H = 2.3, 0.6
# inputs
cC = box(0.3, 4.05, W, H, "Count matrix", IN, INE, "#1b2733")
cS = box(0.3, 3.30, W, H, "Sample sheet", IN, INE, "#1b2733")
cD = box(0.3, 2.20, W, H, "Dilution series", IN, INE, "#1b2733")
cR = box(0.3, 0.55, W, H + 0.1, "Resources\nPPI / GMT / CollecTRI", IN, INE, "#1b2733", fs=9)
# modules
m2 = box(3.8, 3.65, W, H, "2 · QC + diff. expr.", MOD, MODE, "white", bold=True, fs=10)
m1 = box(3.8, 2.20, W, H, "1 · Cross-reactivity", MOD, MODE, "white", bold=True, fs=10)
m3 = box(3.8, 1.25, W, H, "3 · Decomposition", MOD, MODE, "white", bold=True, fs=10)
m4 = box(3.8, 0.45, W, H, "4 · Entropy + pathways", MOD, MODE, "white", bold=True, fs=9.5)
# outputs
o2 = box(7.3, 3.65, W, H, "TMM CPM · DE · QC", OUT, OUTE, "#5a3410", fs=9.5)
o1 = box(7.3, 2.20, W, H, "host fraction + mask", OUT, OUTE, "#5a3410", fs=9.5)
o3 = box(7.3, 1.25, W, H, "heritable / environmental", OUT, OUTE, "#5a3410", fs=9)
o4 = box(7.3, 0.45, W, H, "diversity · entropy · TF", OUT, OUTE, "#5a3410", fs=9)

def R(b):  # right-mid
    return (b[0] + b[2], b[1] + b[3] / 2)
def L(b):
    return (b[0], b[1] + b[3] / 2)
def T(b):
    return (b[0] + b[2] / 2, b[1] + b[3])
def B(b):
    return (b[0] + b[2] / 2, b[1])

# inputs -> modules
arrow(R(cC), (L(m2)[0], L(m2)[1] + 0.08), rad=-0.15)
arrow(R(cS), (L(m2)[0], L(m2)[1] - 0.08), rad=0.12)
arrow(R(cD), L(m1))
arrow(R(cR), L(m4))
# internal
arrow(T(m1), B(m2), label="exclusion\nmask", lx=5.3, ly=2.95)
arrow((L(m2)[0], L(m2)[1] - 0.05), (L(m3)[0], L(m3)[1] + 0.05), rad=0.32,
      label="normalized\nmatrix + fits", lx=3.15, ly=2.55)
arrow((L(m2)[0], L(m2)[1] - 0.12), (L(m4)[0], L(m4)[1] + 0.05), rad=0.42)
# modules -> outputs
for m, o in ((m2, o2), (m1, o1), (m3, o3), (m4, o4)):
    arrow(R(m), L(o))

here = os.path.dirname(os.path.abspath(__file__))
fig.savefig(os.path.join(here, "workflow.svg"), bbox_inches="tight")
fig.savefig(os.path.join(here, "workflow.png"), dpi=150, bbox_inches="tight")
print("wrote workflow.svg and workflow.png")
