# Usage guide

## 1. Install the dependencies

See `docs/REQUIREMENTS.md`. At minimum you need R (4.3 or newer) with the
`yaml` and `edgeR` packages. The other packages are needed only by the
modules that use them, and a module whose package is missing stops with a
message naming it.

## 2. Prepare your inputs

Put your count matrix and sample sheet where the config can find them, in the
format in `docs/INPUT_FORMAT.md`. The `examples/` folder holds a complete
working set you can copy and adapt.

## 3. Write a config

Copy `config.example.yaml` and edit it. The parts you will usually change:

- `project_dir` / `output_dir` — where inputs live and outputs go.
- `inputs:` — the count matrix and sample sheet paths.
- `modules:` — set each of the four to `true` or `false`.
- `qc_de: contrasts:` — the comparisons you want, as `"GroupA - GroupB"`.
- `decomposition:` — which three contrasts form the additive identity
  `total = heritable + environmental`. Leave the module off if your design
  has no in-vitro / in-vivo fork from the same cells.
- thresholds — `treat_fold`, `fdr`, `expressed_floor_cpm`, the
  cross-reactivity `host_tolerance`, and the decomposition bands.

## 4. Run

```
Rscript run.R your_config.yaml
```

The modules run in order. Cross-reactivity (if on) builds a mask that the
qc_de module applies before normalization. Decomposition and entropy both read
the normalized matrix and fits from qc_de, so those two modules require qc_de
to be enabled in the same run.

## 5. Outputs

Everything is written as CSV under `output_dir`:

- `crossreactivity_per_feature.csv`, `exclusion_mask.txt`
- `TMM_logCPM.csv`, `TMM_CPM.csv`, `qc_sample_metrics.csv`,
  `DE_summary.csv`, and `DE_tables/DE_<contrast>.csv`
- `decomposition.csv` (per-gene components, residual, class, site class)
- `diversity_per_sample.csv`, `network_entropy.csv`,
  `signature_singscore.csv`, `progeny_activity.csv`, `tf_activity.csv`
- `sessionInfo.txt` recording the exact package versions of the run

## 6. Check it works

`Rscript tests/test_recovery.R` runs the whole pipeline on the synthetic
example and asserts that it recovers the planted cross-reactive features and
the planted heritable / environmental gene classes. It exits non-zero if any
check fails, so it also serves as a smoke test after installing.
