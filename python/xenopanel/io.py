"""Load the toolkit's CSV outputs into pandas DataFrames.

The R pipeline writes every result as plain CSV; these readers turn a results
directory into tidy frames so downstream work (plotting, further analysis)
can happen in Python. Nothing here recomputes anything -- it only reads.
"""
from __future__ import annotations
import os
import glob
import pandas as pd


class Results:
    """A thin view over a toolkit output directory."""

    def __init__(self, results_dir: str):
        if not os.path.isdir(results_dir):
            raise FileNotFoundError(f"results directory not found: {results_dir}")
        self.dir = results_dir

    def _path(self, *parts: str) -> str:
        return os.path.join(self.dir, *parts)

    def has(self, *parts: str) -> bool:
        return os.path.exists(self._path(*parts))

    # ---- per-module loaders -------------------------------------------------
    def crossreactivity(self):
        p = self._path("crossreactivity_per_feature.csv")
        return pd.read_csv(p) if os.path.exists(p) else None

    def de_summary(self):
        p = self._path("DE_summary.csv")
        return pd.read_csv(p) if os.path.exists(p) else None

    def contrasts(self):
        files = glob.glob(self._path("DE_tables", "DE_*.csv"))
        return sorted(os.path.basename(f)[3:-4] for f in files)

    def de(self, contrast):
        p = self._path("DE_tables", "DE_%s.csv" % contrast)
        if not os.path.exists(p):
            raise FileNotFoundError("no DE table for contrast '%s' in %s" % (contrast, self.dir))
        return pd.read_csv(p)

    def decomposition(self):
        p = self._path("decomposition.csv")
        return pd.read_csv(p) if os.path.exists(p) else None

    def diversity(self):
        p = self._path("diversity_per_sample.csv")
        return pd.read_csv(p) if os.path.exists(p) else None

    def network_entropy(self):
        p = self._path("network_entropy.csv")
        return pd.read_csv(p) if os.path.exists(p) else None

    def signatures(self):
        p = self._path("signature_singscore.csv")
        return pd.read_csv(p) if os.path.exists(p) else None

    def progeny(self):
        p = self._path("progeny_activity.csv")
        return pd.read_csv(p) if os.path.exists(p) else None

    def tf_activity(self):
        p = self._path("tf_activity.csv")
        return pd.read_csv(p) if os.path.exists(p) else None

    def __repr__(self):
        return "Results(%r)" % self.dir
