#!/usr/bin/env python3
"""Draw every applicable plot from a toolkit results directory.

Usage:
    python plot_results.py <results_dir> [--out <figures_dir>] [--fdr 0.05]

Writes one PNG per plot. Only the plots whose inputs exist are produced, so it
works whatever subset of modules was run.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from xenopanel import Results, plots


def main():
    ap = argparse.ArgumentParser(description="Plot DriverMap toolkit results.")
    ap.add_argument("results_dir", help="the toolkit output directory (holds the CSV results)")
    ap.add_argument("--out", default=None, help="where to write PNGs (default: <results_dir>/figures)")
    ap.add_argument("--fdr", type=float, default=0.05)
    args = ap.parse_args()

    res = Results(args.results_dir)
    out = args.out or os.path.join(args.results_dir, "figures")
    os.makedirs(out, exist_ok=True)
    made = []

    xr = res.crossreactivity()
    if xr is not None:
        fig = plots.crossreactivity_hist(xr)
        f = os.path.join(out, "crossreactivity.png"); fig.savefig(f, dpi=150, bbox_inches="tight"); made.append(f)

    for cn in res.contrasts():
        fig = plots.volcano(res.de(cn), contrast=cn, fdr=args.fdr)
        f = os.path.join(out, "volcano_%s.png" % cn); fig.savefig(f, dpi=150, bbox_inches="tight"); made.append(f)

    dec = res.decomposition()
    if dec is not None:
        fig = plots.decomposition_scatter(dec)
        f = os.path.join(out, "decomposition.png"); fig.savefig(f, dpi=150, bbox_inches="tight"); made.append(f)

    ne = res.network_entropy()
    if ne is not None:
        for meas in ("CCAT", "SR"):
            if meas in ne.columns and ne[meas].notna().any():
                fig = plots.entropy_by_group(ne, measure=meas)
                f = os.path.join(out, "entropy_%s.png" % meas); fig.savefig(f, dpi=150, bbox_inches="tight"); made.append(f)

    print("wrote %d figure(s) to %s:" % (len(made), out))
    for f in made:
        print("  " + os.path.basename(f))
    if not made:
        print("  (no plottable outputs found -- run the pipeline first)")


if __name__ == "__main__":
    main()
