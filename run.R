#!/usr/bin/env Rscript
# ============================================================================
# XenoPanel toolkit -- main driver
# Usage:  Rscript run.R [config.yaml]        (default: config.example.yaml)
#         Rscript run.R --version
# Reads the config, runs the enabled modules in order, writes all outputs
# under <project_dir>/<output_dir>. See README.md and docs/INPUT_FORMAT.md.
# The same pipeline is available as a function after installing the package:
#   XenoPanel::xenopanel_run("my_config.yaml")
# ============================================================================
suppressWarnings(suppressMessages({
  here <- tryCatch(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                   error = function(e) ".")
  if (!length(here) || !nzchar(here)) here <- "."
  for (f in c("io.R", "mod_crossreactivity.R", "mod_qc_de.R",
              "mod_decomposition.R", "mod_entropy_pathways.R", "pipeline.R"))
    source(file.path(here, "R", f))
}))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1 && args[1] %in% c("--version", "-V")) { cat(sprintf("XenoPanel %s\n", XENOPANEL_VERSION)); quit(status = 0) }
cfg_path <- if (length(args) >= 1) args[1] else file.path(here, "config.example.yaml")
xenopanel_run(cfg_path)
