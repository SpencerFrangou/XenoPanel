# pipeline.R -- the whole XenoPanel run as one function, so the toolkit can be used as an installed
# package (XenoPanel::xenopanel_run("config.yaml")) as well as through run.R.

`%||%` <- function(a, b) if (is.null(a)) b else a

XENOPANEL_VERSION <- "1.1.0.dev"

#' Version string of the toolkit
xenopanel_version <- function() XENOPANEL_VERSION

#' Run the XenoPanel pipeline from a YAML configuration file
#'
#' Reads the configuration, runs the enabled modules in order (cross-reactivity, QC and
#' differential expression, decomposition, entropy and pathway activity) and writes every
#' output under <project_dir>/<output_dir>. Returns the output directory invisibly.
xenopanel_run <- function(cfg_path) {
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
  say("XenoPanel %s", XENOPANEL_VERSION)
  say("  config : %s", normalizePath(cfg_path))
  say("  project: %s", project_dir)
  say("  output : %s", out_dir)
  say("  modules: %s", paste(names(Filter(isTRUE, mods)), collapse = ", "))
  say("============================================================")

  # ---- shared inputs (needed by qc_de / decomposition / entropy) -----------
  m <- NULL; samples <- NULL
  if (any(unlist(mods[c("qc_de", "decomposition", "entropy_pathways")]))) {
    samples <- read_samples(resolve_path(cfg$inputs$samples, project_dir))
    m <- read_counts(resolve_path(cfg$inputs$counts, project_dir),
                     collapse = isTRUE(cfg$qc_de$collapse_duplicate_symbols))
    m <- align_counts(m, samples)
    say("loaded counts: %d features x %d samples", nrow(m), ncol(m))
  }

  # ---- 1. cross-reactivity -------------------------------------------------
  mask <- NULL
  if (isTRUE(mods$crossreactivity)) {
    mask <- run_crossreactivity(cfg, project_dir, out_dir)
    if (!isTRUE(cfg$crossreactivity$apply_mask_in_de)) mask <- NULL
  }

  # ---- 2. qc + differential expression ------------------------------------
  de <- NULL
  if (isTRUE(mods$qc_de)) {
    de <- run_qc_de(cfg, m, samples, out_dir, mask = mask)
  } else if (any(unlist(mods[c("decomposition", "entropy_pathways")]))) {
    die("modules 'decomposition' and 'entropy_pathways' need 'qc_de' enabled (they use its normalized matrix and fits)")
  }

  # ---- 3. decomposition ----------------------------------------------------
  if (isTRUE(mods$decomposition)) run_decomposition(cfg, de, out_dir)

  # ---- 4. entropy and pathway activity ------------------------------------
  if (isTRUE(mods$entropy_pathways)) run_entropy_pathways(cfg, de, samples, out_dir, project_dir)

  say("============================================================")
  say("done. outputs in %s", out_dir)
  writeLines(c(sprintf("XenoPanel %s", XENOPANEL_VERSION), sprintf("run: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
               sprintf("config: %s", normalizePath(cfg_path)),
               sprintf("config md5: %s", tryCatch(as.character(tools::md5sum(cfg_path)), error = function(e) "n/a")),
               sprintf("modules: %s", paste(names(Filter(isTRUE, mods)), collapse = ", ")), "",
               utils::capture.output(utils::sessionInfo())),
             file.path(out_dir, "sessionInfo.txt"))
  invisible(out_dir)
}
