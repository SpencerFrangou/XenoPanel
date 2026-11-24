# Changelog

All notable changes to XenoPanel are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow [Semantic Versioning](https://semver.org/).
A patch release fixes a bug without changing results, a minor release adds a feature or changes a default in a
way that can change results, a major release changes the configuration format or the input layout.

## [Unreleased]

### Added
- `Rscript run.R --version`; the version is printed in the run banner and written, with the run time, the
  configuration path and its checksum and the enabled modules, at the top of `sessionInfo.txt` in every output
  folder.
- GitHub Actions workflow that installs the R dependencies, generates the synthetic example, runs the
  pipeline and the recovery test on Linux, macOS and Windows, exercises the Python plotting layer, and
  repeats monthly against current CRAN and Bioconductor releases.
- `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`, issue templates, `docs/statement_of_need.md`.

### Fixed
- The run banner and the Python plotting script still said "DriverMap toolkit" after the rename to XenoPanel.

## [1.0.0] - 2026-10-29

First public release: cross-reactivity estimation and exclusion mask from a dilution series; QC, TMM
normalisation and edgeR quasi-likelihood differential expression with TREAT; heritable/environmental
decomposition over three states; transcriptome diversity, network entropy (SCENT, CCAT), singscore, PROGENy
and CollecTRI activity; YAML configuration; synthetic example with recovery test; optional Python plotting
layer; GitHub Pages site. Archived at Zenodo, DOI 10.5281/zenodo.23768221.

[Unreleased]: https://github.com/SpencerFrangou/XenoPanel/compare/v1.0.0...HEAD
[1.0.0]: https://github.com/SpencerFrangou/XenoPanel/releases/tag/v1.0.0
