#!/usr/bin/env python3
"""Build docs/example_outputs.png: a 2x2 montage of representative outputs,
drawn from the synthetic example run with the toolkit's own plot functions."""
import os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, "..", "python"))
from xenopanel import Results, plots

res = Results(os.path.join(here, "..", "examples", "results"))
fig, axs = plt.subplots(2, 2, figsize=(11, 9))
plots.crossreactivity_hist(res.crossreactivity(), ax=axs[0, 0])
plots.volcano(res.de("MetLung_vs_Parental"), contrast="MetLung vs Parental", ax=axs[0, 1])
plots.decomposition_scatter(res.decomposition(), ax=axs[1, 0])
plots.entropy_by_group(res.network_entropy(), measure="CCAT", ax=axs[1, 1])
fig.suptitle("Representative outputs (synthetic example)", fontsize=13, y=0.995)
fig.tight_layout(rect=[0, 0, 1, 0.98])
fig.savefig(os.path.join(here, "example_outputs.png"), dpi=140, bbox_inches="tight")
print("wrote docs/example_outputs.png")
