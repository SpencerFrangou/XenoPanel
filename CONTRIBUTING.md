# Contributing to XenoPanel

Bug reports, questions about an unexpected result, requests for an input format or a method that is not
handled, and pull requests are welcome.

## Reporting a problem

Open an issue at <https://github.com/SpencerFrangou/XenoPanel/issues> using the bug report template. Include
the XenoPanel version (`Rscript run.R --version`, or the first line of `sessionInfo.txt` in the output
folder), the operating system and R version, the configuration file with paths removed, what you ran, what
happened and what you expected. Do not attach count matrices that you are not free to share; the synthetic
generator in `examples/make_example_data.R` can usually be adapted to reproduce a problem.

## Asking for help

Questions about configuring a run, choosing thresholds, or interpreting the decomposition or entropy outputs
can be opened as issues (label `question`). Expect a reply within a week.

## Proposing a change

1. Open an issue first for anything larger than a typo, so the design can be discussed before work starts.
2. Fork the repository and create a branch from `main`.
3. Install the R packages listed in `docs/REQUIREMENTS.md`, then run the recovery test, which generates the
   synthetic data, runs the whole pipeline and checks that the planted truth is recovered:
   ```
   Rscript tests/test_recovery.R
   ```
   For the Python plotting layer: `pip install -r python/requirements.txt && python python/plot_results.py examples/results`.
4. Keep the analysis deterministic (fixed seeds, no network access during a run) and keep research data out
   of the repository.
5. Add or extend a check in `tests/test_recovery.R` for the change, add a line under `Unreleased` in
   `CHANGELOG.md`, and update `USAGE.md`, `docs/INPUT_FORMAT.md` or `config.example.yaml` if the behaviour,
   an input or a configuration key changed.
6. Open a pull request. The GitHub Actions workflow runs the recovery test on Linux, macOS and Windows.

A change that alters results for existing data (a filtering threshold, the model formula, the decomposition
rule) is a minor version bump and must be described in the changelog with the reason and the expected size of
the effect.

## Code style

The R code is organised as one file per module, each exposing a single `run_<module>()` entry point called
from `run.R`, with shared helpers in `R/io.R`. Follow the existing style, keep configuration keys documented
in `config.example.yaml`, and fail early with a clear message (`die()`) when an input or a package is missing.

## Licence

By contributing you agree that your contribution is licensed under the MIT licence of this repository.
