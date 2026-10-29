# io.R -- input reading, validation, and small shared helpers.
# No analysis lives here; every reader fails loudly on malformed input.

say <- function(...) cat(sprintf(...), "\n", sep = "")
die <- function(...) stop(sprintf(...), call. = FALSE)

# Resolve a path against project_dir unless it is already absolute.
resolve_path <- function(path, project_dir) {
  if (is.null(path) || !nzchar(path)) return(NULL)
  if (grepl("^(/|~|[A-Za-z]:)", path)) path else file.path(project_dir, path)
}

# Read a count matrix. Column 1 is the feature ID; the rest are samples.
# A "LibID|SYMBOL" feature ID is reduced to SYMBOL. Duplicate symbols are an
# error unless `collapse = TRUE`, in which case their counts are summed.
read_counts <- function(path, collapse = FALSE) {
  if (!file.exists(path)) die("count file not found: %s", path)
  d <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(d) < 2) die("count file needs a feature column and >=1 sample column")
  ids <- as.character(d[[1]])
  sym <- sub("^[0-9]+\\|", "", ids)          # strip DriverMap LibID prefix if present
  m <- as.matrix(d[, -1, drop = FALSE])
  mode(m) <- "numeric"
  if (any(is.na(m))) die("count matrix has %d non-numeric/NA cells", sum(is.na(m)))
  if (any(m < 0)) die("count matrix has negative values")
  rownames(m) <- sym
  dup <- unique(sym[duplicated(sym)])
  if (length(dup)) {
    if (!collapse)
      die("%d duplicate gene symbols (e.g. %s). Set collapse_duplicate_symbols: true to sum them.",
          length(dup), paste(head(dup, 3), collapse = ", "))
    m <- rowsum(m, group = sym)
    say("  collapsed %d duplicate symbols by summation", length(dup))
  }
  storage.mode(m) <- "double"
  m
}

# Read and validate a sample sheet.
read_samples <- function(path) {
  if (!file.exists(path)) die("sample sheet not found: %s", path)
  s <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  need <- c("sample", "group")
  miss <- setdiff(need, names(s))
  if (length(miss)) die("sample sheet missing column(s): %s", paste(miss, collapse = ", "))
  if (anyDuplicated(s$sample)) die("duplicate sample names in sample sheet")
  if (!"context" %in% names(s)) s$context <- NA_character_
  if (!"replicate" %in% names(s)) s$replicate <- ave(seq_len(nrow(s)), s$group, FUN = seq_along)
  s
}

# Align a count matrix to a sample sheet, honouring source_column if present.
align_counts <- function(m, samples) {
  if ("source_column" %in% names(samples) && !all(is.na(samples$source_column))) {
    miss <- setdiff(samples$source_column, colnames(m))
    if (length(miss)) die("source_column value(s) not found in count file: %s",
                          paste(head(miss, 5), collapse = ", "))
    m <- m[, samples$source_column, drop = FALSE]
  } else {
    miss <- setdiff(samples$sample, colnames(m))
    if (length(miss)) die("sample(s) not found as count columns: %s (give source_column to map)",
                          paste(head(miss, 5), collapse = ", "))
    m <- m[, samples$sample, drop = FALSE]
  }
  colnames(m) <- samples$sample
  m
}

# Load an optional package, returning TRUE if available.
have <- function(pkg) requireNamespace(pkg, quietly = TRUE)

# Require a package or stop with an actionable message.
need_pkg <- function(pkg, what) {
  if (!have(pkg)) die("module '%s' needs the R package '%s' (install it or disable the module)", what, pkg)
}

ensure_dir <- function(p) { dir.create(p, showWarnings = FALSE, recursive = TRUE); p }
