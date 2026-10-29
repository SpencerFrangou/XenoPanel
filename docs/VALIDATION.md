# Validation

Two checks. Neither ships research data.

## 1. Recovery of planted truth (synthetic data, included)

`examples/make_example_data.R` simulates a dataset with a known structure: 150
genes heritable-up, 150 environmental-up, and a dilution series in which 150
genes carry a constant host-species template. `tests/test_recovery.R` runs the
pipeline and asserts recovery. Result:

- cross-reactivity: 100% recall of expressed planted host features; ~10% of
  masked features non-planted (the tolerance margin).
- decomposition: 92% recall of heritable genes, 99% of environmental genes;
  additive identity to 4e-16 log2; no meaningful cross-class confusion.
- entropy: diversity and network-entropy (CCAT) tables for every sample.

Run with `Rscript tests/test_recovery.R`.

## 2. Reproduction of a reference pipeline (data not included)

The engine was run on the inputs of the study it was abstracted from, in a
scratch directory outside this repository, and compared against that study's
independent, study-specific scripts. On the 5,677 shared genes:

- differential expression: total-contrast log fold-change, Pearson r =
  1.000000, maximum absolute difference 0.0000.
- decomposition: 100% class agreement, reproducing 131 heritable, 320
  codominant, 493 environmental, 11 undetermined genes, and 610
  in-vivo-general, 198 site-specific, 5 site-undetermined.

The scratch data were removed after the comparison.
