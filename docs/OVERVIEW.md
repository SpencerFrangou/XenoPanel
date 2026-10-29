# Overview: premise, design, and the methods the toolkit builds on

## The premise

Targeted transcriptomic panels (DriverMap and similar assays) measure a defined
set of genes using multiplexed, primer-based amplification rather than
whole-transcriptome sequencing. In practice, that makes them relatively cheap,
sensitive at low input, and useful for small or otherwise difficult samples.
The trade-off is that they introduce a couple of analysis problems that a
general RNA-seq pipeline is unlikely to solve on its own.

The first is cross-species signal. A panel designed for one species does not
simply stop amplifying because RNA from a second species is present. In a
xenograft — for example, a human tumour growing in a mouse — some primers may
also amplify mouse template. Part of a measured value can therefore come from
the host rather than the tumour. Before making biological comparisons, that
host-derived fraction needs to be estimated and the features with the worst
cross-reactivity removed.

The second issue comes after the signal has been cleaned up. At that point, the
more interesting question is often not just *which* genes change, but *why*
they change. If the same cells are profiled in culture and in tissue, a gene's
change between those states can be separated into a component the cells carry
with them (heritable) and a component imposed by the tissue environment
(environmental). Provided the relevant states were profiled, that split is an
exact piece of arithmetic. It can also be more informative than a conventional
differential-expression list on its own.

The toolkit brings those two ideas together with standard differential
expression and several transcriptome-state read-outs. They are organized as
four modules that can be switched on or off in a config file. The pipeline is
intentionally input-agnostic: it accepts any count matrix and sample sheet that
follow the documented format, along with whatever contrasts the analysis
requires.

## The four modules

**Cross-reactivity.** This module starts from a dilution series in which
target-species RNA is titrated into a constant mass of host RNA. For each
feature, it estimates the fraction of signal that appears to be host-derived
and masks features that exceed a user-defined tolerance. The estimator uses two
non-overlapping calibration sets from the titration: features with no host
signal and features dominated by host signal. That setup means the
decomposition does not need an assumption about PCR yield.

**QC and differential expression.** The workflow applies expression filtering
and TMM normalization, followed by per-sample quality checks such as library
size, detection, and replicate correlation. Differential expression is then
tested using quasi-likelihood F-tests and TREAT against a fold-change threshold,
for whichever contrasts are defined in the configuration.

**Decomposition.** For designs with in-vitro and in-vivo arms derived from the
same cells, each gene's total change is separated as
`total = heritable + environmental`. This is an exact additive identity rather
than a fitted approximation. Genes are then assigned to classes using bands,
with an optional site-specific split.

**Entropy and pathway activity.** This module collects several measures that
describe transcriptome state from different angles: transcriptome diversity
(Shannon entropy, Hill numbers, N50, and Gini), network entropy (signalling
entropy rate and CCAT), gene-signature scores, pathway activity, and
transcription-factor activity. These measures are related, but they are not
interchangeable; each captures a somewhat different property of the data.

## The established methods it uses, and their own descriptions

The toolkit does not reinvent the statistics. The differential-expression and
activity-inference steps call published, peer-reviewed packages; what the
toolkit adds is the panel-specific logic (cross-reactivity, decomposition) and
the configuration layer that ties everything together. The packages, with
their own descriptions:

- **edgeR** (Bioconductor) — "Differential expression analysis of sequence
  count data. Implements a range of statistical methodology based on the
  negative binomial distributions, including empirical Bayes estimation, exact
  tests, generalized linear models, quasi-likelihood, and gene set
  enrichment." Used for filtering, normalization and all differential
  expression. Robinson, McCarthy & Smyth 2010; Chen, Lun & Smyth 2016.

- **SCENT** (Single Cell ENTropy; github.com/aet21/SCENT) — estimates
  signalling entropy, a systems measure of differentiation potency, by
  integrating an expression profile with a protein–protein interaction
  network. Used for the signalling entropy rate (SR) and CCAT. Banerji et al.
  2013; Teschendorff & Enver 2017; Chen et al. 2019 (CCAT).

- **decoupleR** (Bioconductor) — "Ensemble of computational methods to infer
  biological activities from omics data." Used to run the PROGENy and
  transcription-factor activity models. Badia-i-Mompel et al. 2022.

- **PROGENy** (Bioconductor) — "PROGENy is a resource that leverages a large
  compendium of publicly available signaling perturbation experiments to yield
  a common core of pathway responsive genes for human and mouse. These,
  coupled with any statistical method, can be used to infer pathway activities
  from bulk or single-cell transcriptomics." Schubert et al. 2018.

- **singscore** (Bioconductor) — "A simple single-sample gene signature
  scoring method that uses rank-based statistics to analyze the sample's gene
  expression profile. It scores the expression activities of gene sets at a
  single-sample level." Foroutan et al. 2018; Bhuva et al. 2020.

- **igraph** (CRAN) — network analysis and manipulation routines, used to
  build the protein-interaction network for the entropy measures.

Transcription-factor regulons are read from a CollecTRI-format network, and
gene sets from a standard GMT file (for example MSigDB), both supplied by the
user so the toolkit carries no third-party gene-set data of its own.

## How to run it

```
Rscript run.R your_config.yaml
```

That command runs the modules enabled in the config, in order, and writes the
results as CSV files under the output directory. A sensible first check is to
run the bundled synthetic example before using a real dataset:

```
Rscript examples/make_example_data.R     # make the example inputs
Rscript run.R config.example.yaml        # run the whole pipeline
Rscript tests/test_recovery.R            # check it recovers the planted truth
```

`USAGE.md` provides the step-by-step workflow, and `docs/INPUT_FORMAT.md`
defines the expected input files.
