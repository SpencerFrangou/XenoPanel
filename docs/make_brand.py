#!/usr/bin/env python3
"""Build docs/banner.png (README header) and docs/social_preview.png
(the card GitHub shows when the repo is shared). Clean, typographic, no clip art."""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

here = os.path.dirname(os.path.abspath(__file__))
NAVY, BLUE, LIGHT, PAPER, INK = "#16222e", "#2c7fb8", "#eef2f7", "#fbfcfd", "#1b2733"
MODS = ["Cross-reactivity", "QC + diff. expr.", "Decomposition", "Entropy + pathways"]


def chips(ax, y, x0=0.055, fs=12, gap=0.015):
    x = x0
    for m in MODS:
        w = 0.013 * len(m) + 0.03
        ax.add_patch(FancyBboxPatch((x, y), w, 0.1, boxstyle="round,pad=0.004,rounding_size=0.03",
                                    transform=ax.transAxes, fc="#ffffff22", ec=BLUE, lw=1.2))
        ax.text(x + w / 2, y + 0.05, m, transform=ax.transAxes, ha="center", va="center",
                color=LIGHT, fontsize=fs)
        x += w + gap


def banner():
    fig = plt.figure(figsize=(12.8, 3.0), dpi=100)
    ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off")
    ax.add_patch(plt.Rectangle((0, 0), 1, 1, transform=ax.transAxes, fc=NAVY, ec="none"))
    ax.add_patch(plt.Rectangle((0, 0.0), 0.012, 1, transform=ax.transAxes, fc=BLUE, ec="none"))
    ax.text(0.055, 0.66, "XenoPanel", transform=ax.transAxes, color="#ffffff",
            fontsize=40, fontweight="bold", va="center")
    ax.text(0.056, 0.40, "Pathway, network and species-resolved analysis for targeted transcriptomic panels",
            transform=ax.transAxes, color="#a9c3d6", fontsize=15, va="center")
    chips(ax, 0.12)
    fig.savefig(os.path.join(here, "banner.png"), dpi=100)
    plt.close(fig)


def social():
    fig = plt.figure(figsize=(12.8, 6.4), dpi=100)
    ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off")
    ax.add_patch(plt.Rectangle((0, 0), 1, 1, transform=ax.transAxes, fc=NAVY, ec="none"))
    ax.add_patch(plt.Rectangle((0, 0.86), 1, 0.14, transform=ax.transAxes, fc=BLUE, ec="none"))
    ax.text(0.5, 0.93, "TARGETED TRANSCRIPTOMIC PANEL ANALYSIS", transform=ax.transAxes,
            ha="center", va="center", color="#eaf3fa", fontsize=15, fontweight="bold",
            fontfamily="monospace")
    ax.text(0.5, 0.60, "XenoPanel", transform=ax.transAxes, ha="center", va="center",
            color="#ffffff", fontsize=58, fontweight="bold")
    ax.text(0.5, 0.45, "Cross-reactivity  ·  QC + differential expression  ·  decomposition  ·  entropy & pathways",
            transform=ax.transAxes, ha="center", va="center", color="#a9c3d6", fontsize=18)
    ax.text(0.5, 0.17, "Configuration-driven  ·  R core with an optional Python layer  ·  MIT",
            transform=ax.transAxes, ha="center", va="center", color="#7f9bb0", fontsize=15)
    ax.plot([0.33, 0.67], [0.29, 0.29], transform=ax.transAxes, color=BLUE, lw=2)
    fig.savefig(os.path.join(here, "social_preview.png"), dpi=100)
    plt.close(fig)


banner()
social()
print("wrote docs/banner.png and docs/social_preview.png")
