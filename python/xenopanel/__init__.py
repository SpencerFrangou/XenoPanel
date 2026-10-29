"""xenopanel -- Python companion to the XenoPanel toolkit.

The analysis runs in R (see the repository root). This package loads the R
pipeline's CSV outputs into pandas and draws the standard plots, so results
can be taken into a Python workflow for figures and further analysis. It does
not reimplement any of the statistics.
"""
from .io import Results
from . import plots

__all__ = ["Results", "plots"]
__version__ = "0.1.0"
