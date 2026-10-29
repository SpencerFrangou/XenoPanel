"""Standard plots for the toolkit's outputs.

Plotting was deliberately left out of the R pipeline, which writes tables
only. These functions draw the figures a user most often wants, from those
tables, using matplotlib. Each takes a loaded frame (or a Results object) and
returns a matplotlib Figure, so they compose into scripts or notebooks.
"""
from __future__ import annotations
import numpy as np
import matplotlib
matplotlib.use("Agg")  # safe default for headless runs; override before import for interactive use
import matplotlib.pyplot as plt

_CLASS_COLORS = {
    "heritable": "#2c7fb8", "environmental": "#d95f0e",
    "codominant": "#756bb1", "undetermined": "#969696", "not_changed": "#d9d9d9",
}


def volcano(de, contrast="", fdr=0.05, lfc=np.log2(1.5), ax=None):
    """logFC against -log10(TREAT FDR) for one contrast."""
    fig = None
    if ax is None:
        fig, ax = plt.subplots(figsize=(5, 4.5))
    q = de["TREAT_FDR"].clip(lower=1e-300)
    sig = (de["TREAT_FDR"] < fdr)
    ax.scatter(de.loc[~sig, "logFC"], -np.log10(q[~sig]), s=6, c="#cccccc", linewidths=0, label="n.s.")
    ax.scatter(de.loc[sig, "logFC"], -np.log10(q[sig]), s=8, c="#d95f0e", linewidths=0, label="TREAT FDR<%.2g" % fdr)
    ax.axhline(-np.log10(fdr), color="#888888", lw=0.7, ls="--")
    for x in (-lfc, lfc):
        ax.axvline(x, color="#888888", lw=0.7, ls="--")
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("-log10 TREAT FDR")
    ax.set_title("Volcano: %s" % contrast if contrast else "Volcano")
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    _tidy(ax)
    return fig if fig is not None else ax.figure


def decomposition_scatter(dec, ax=None):
    """Heritable vs environmental component, coloured by class."""
    fig = None
    if ax is None:
        fig, ax = plt.subplots(figsize=(5, 4.6))
    d = dec[dec["class"] != "not_changed"]
    for cls, sub in d.groupby("class"):
        ax.scatter(sub["heritable"], sub["environmental"], s=10, linewidths=0,
                   c=_CLASS_COLORS.get(cls, "#444444"), label="%s (%d)" % (cls, len(sub)))
    lim = np.nanmax(np.abs(np.r_[d["heritable"], d["environmental"]])) * 1.05
    ax.plot([-lim, lim], [lim, -lim], color="#bbbbbb", lw=0.7)  # total-constant line
    ax.axhline(0, color="#dddddd", lw=0.6); ax.axvline(0, color="#dddddd", lw=0.6)
    ax.set_xlabel("heritable component (log2)")
    ax.set_ylabel("environmental component (log2)")
    ax.set_title("Decomposition of the program")
    ax.legend(frameon=False, fontsize=8, loc="best")
    _tidy(ax)
    return fig if fig is not None else ax.figure


def entropy_by_group(ne, measure="CCAT", ax=None):
    """Strip of a network-entropy measure (CCAT or SR) by group."""
    fig = None
    if ax is None:
        fig, ax = plt.subplots(figsize=(5, 4))
    if measure not in ne.columns:
        raise ValueError("measure '%s' not in network_entropy table" % measure)
    groups = list(dict.fromkeys(ne["group"]))
    for i, g in enumerate(groups):
        y = ne.loc[ne["group"] == g, measure].values
        x = np.random.default_rng(0).normal(i, 0.05, size=len(y))
        ax.scatter(x, y, s=24, c="#2c7fb8", linewidths=0)
        ax.hlines(np.median(y), i - 0.2, i + 0.2, color="#d95f0e", lw=2)
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels(groups, rotation=30, ha="right", fontsize=8)
    ax.set_ylabel(measure)
    ax.set_title("%s by group" % measure)
    _tidy(ax)
    return fig if fig is not None else ax.figure


def crossreactivity_hist(xr, tolerance=0.10, ax=None):
    """Distribution of host-derived signal fraction, with the mask threshold."""
    fig = None
    if ax is None:
        fig, ax = plt.subplots(figsize=(5, 4))
    v = xr.loc[xr["expressed"], "host_fraction_op"].clip(0, 1)
    ax.hist(v, bins=40, color="#9ecae1", edgecolor="white")
    ax.axvline(tolerance, color="#d95f0e", lw=1.2, ls="--", label="mask > %.0f%%" % (100 * tolerance))
    ax.set_yscale("log")
    ax.set_xlabel("host-derived fraction of signal (expressed features)")
    ax.set_ylabel("features (log scale)")
    ax.set_title("Cross-reactivity")
    ax.legend(frameon=False, fontsize=8)
    _tidy(ax)
    return fig if fig is not None else ax.figure


def _tidy(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
