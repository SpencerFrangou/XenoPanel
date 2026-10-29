# Requirements

## Software

**R 4.3 or newer.** The toolkit was developed and validated on R 4.6.0.

Packages, by the module that needs them. A module stops with a clear message
if one of its packages is missing, so you only need to install what you use.

| Package      | Source          | Needed by                          |
|--------------|-----------------|------------------------------------|
| yaml         | CRAN            | always (reads the config)          |
| edgeR        | Bioconductor    | qc_de, decomposition               |
| limma        | Bioconductor    | pulled in by edgeR                 |
| igraph       | CRAN            | entropy_pathways (network entropy) |
| SCENT        | GitHub (aet21)  | entropy_pathways (signalling entropy rate; CCAT works without it) |
| singscore    | Bioconductor    | entropy_pathways (signature scoring) |
| decoupleR    | Bioconductor    | entropy_pathways (PROGENy, TF activity) |
| progeny      | Bioconductor    | entropy_pathways (PROGENy)         |
| vegan        | CRAN            | optional diversity extras          |

Versions used in validation: edgeR 4.10.0, limma 3.68.2, SCENT 1.0.3,
igraph 2.3.3, singscore 1.32.0, GSVA 2.6.6, decoupleR 2.17.0, progeny 1.34.0,
vegan 2.7.6, yaml 2.3.12. The exact versions of any run are written to
`sessionInfo.txt` in the output folder.

Install the CRAN and Bioconductor packages with:

```r
install.packages(c("yaml", "igraph", "vegan"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("edgeR", "limma", "singscore", "decoupleR", "progeny"))
# SCENT is on GitHub:
# remotes::install_github("aet21/SCENT")
```

## Hardware

Targeted panels are small (8,000–20,000 features), so a run fits in the memory
of an ordinary laptop; no GPU or cluster is needed. The synthetic example
(2,000 features, 15 samples) runs in under a minute. The signalling entropy
rate is the slowest step and completes in minutes for a panel-sized network.

## Python (optional plotting layer only)

Python 3.8+ with `pandas`, `matplotlib` and `numpy`
(`pip install -r python/requirements.txt`). Not needed for the analysis.

## Operating system

Any platform with R: macOS, Linux or Windows. Paths in the config may use
either separator.
