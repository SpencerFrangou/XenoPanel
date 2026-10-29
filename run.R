#!/usr/bin/env Rscript
# ============================================================================
# XenoPanel toolkit -- main driver
# Usage:  Rscript run.R [config.yaml]        (default: config.example.yaml)
# Reads the config, runs the enabled modules in order, writes all outputs
# under <project_dir>/<output_dir>. See README.md and docs/INPUT_FORMAT.md.
# ============================================================================
suppressWarnings(suppressMessages({
  here <- tryCatch(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                   error = function(e) ".")
  if (!length(here) || !nzchar(here)) here <- "."
  for (f in c("io.R", "mod_crossreactivity.R", "mod_qc_de.R",
              "mod_decomposition.R", "mod_entropy_pathways.R"))
    source(file.path(here, "R", f))
}))
`%||%` <- function(a, b) if (is.null(a)) b else a

args <- commandArgs(trailingOnly = TRUE)
cfg_path <- if (length(args) >= 1) args[1] else file.path(here, "config.example.yaml")
if (!file.exists(cfg_path)) die("config file not found: %s", cfg_path)
if (!have("yaml")) die("the 'yaml' package is required to read the config")
cfg <- yaml::read_yaml(cfg_path)

project_dir <- as.character(cfg$project_dir %||% ".")
if (is.na(project_dir) || !nzchar(project_dir)) project_dir <- "."
if (!grepl("^(/|~|[A-Za-z]:)", project_dir))
  project_dir <- file.path(dirname(normalizePath(cfg_path)), project_dir)
out_dir <- ensure_dir(resolve_path(cfg$output_dir %||% "results", project_dir))
mods <- cfg$modules %||% list()

say("============================================================")
say("DriverMap toolkit")
say("  config : %s", normalizePath(cfg_path))
say("  project: %s", project_dir)
say("  output : %s", out_dir)
say("  modules: %s", paste(names(Filter(isTRUE, mods)), collapse = ", "))
say("============================================================")

# ---- shared inputs (needed by qc_de / decomposition / entropy) -------------
m <- NULL; samples <- NULL
if (any(unlist(mods[c("qc_de", "decomposition", "entropy_pathways")]))) {
  samples <- read_samples(resolve_path(cfg$inputs$samples, project_dir))
  m <- read_counts(resolve_path(cfg$inputs$counts, project_dir),
                   collapse = isTRUE(cfg$qc_de$collapse_duplicate_symbols))
  m <- align_counts(m, samples)
  say("loaded counts: %d features x %d samples", nrow(m), ncol(m))
}

# ---- 1. cross-reactivity ---------------------------------------------------
mask <- NULL
if (isTRUE(mods$crossreactivity)) {
  mask <- run_crossreactivity(cfg, project_dir, out_dir)
  if (!isTRUE(cfg$crossreactivity$apply_mask_in_de)) mask <- NULL
}

# ---- 2. qc + differential expression --------------------------------------
de <- NULL
if (isTRUE(mods$qc_de)) {
  de <- run_qc_de(cfg, m, samples, out_dir, mask = mask)
} else if (any(unlist(mods[c("decomposition", "entropy_pathways")]))) {
  die("modules 'decomposition' and 'entropy_pathways' need 'qc_de' enabled (they use its normalized matrix and fits)")
}

# ---- 3. decomposition ------------------------------------------------------
if (isTRUE(mods$decomposition)) run_decomposition(cfg, de, out_dir)

# ---- 4. entropy and pathway activity --------------------------------------
if (isTRUE(mods$entropy_pathways)) run_entropy_pathways(cfg, de, samples, out_dir, project_dir)

say("============================================================")
say("done. outputs in %s", out_dir)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
