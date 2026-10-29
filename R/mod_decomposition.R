# mod_decomposition.R
# Split each gene's total change into a heritable part and an environmental
# part using an exact additive identity, then attribute each gene to a class
# in bands. This requires the three contrasts to satisfy
#   total = heritable + environmental
# by construction (e.g. (MetLung-Parental) = (MetLungCL-Parental) + (MetLung-MetLungCL)).
# The residual is reported; it should be ~0 to numerical precision.
#
# Banded attribution (component evidence from the unthresholded QL test, so
# the fold-change threshold is applied once, to the total):
#   heritable      heritable fraction > band_high, heritable component significant
#   environmental  heritable fraction < band_low,  environmental component significant
#   codominant     fraction in between, both components significant
#   undetermined   otherwise

run_decomposition <- function(cfg, de, out_dir) {
  d <- cfg$decomposition
  for (k in c("total", "heritable", "environmental"))
    if (is.null(d[[k]]) || !d[[k]] %in% names(de$LFC))
      die("decomposition: contrast '%s' (%s) not found among qc_de contrasts", d[[k]], k)

  g <- de$genes
  tot <- de$LFC[[d$total]][g]; her <- de$LFC[[d$heritable]][g]; env <- de$LFC[[d$environmental]][g]
  out <- data.frame(gene = g, total = tot, heritable = her, environmental = env,
                    residual = tot - (her + env),
                    heritable_fraction = her / tot, row.names = NULL)
  say("[decomposition] max |residual| of total = heritable + environmental: %.2e log2",
      max(abs(out$residual), na.rm = TRUE))

  # expressed floor from the two endpoint states of the total contrast
  parts <- strsplit(cfg$qc_de$contrasts[[d$total]], "\\s*-\\s*")[[1]]
  floor <- de$floor
  hi_ok <- rep(TRUE, length(g))
  if (floor > 0 && all(parts %in% colnames(de$gmean))) {
    hi <- pmax(de$gmean[g, parts[1]], de$gmean[g, parts[2]])
    hi_ok <- hi >= floor
    say("  expressed floor %g CPM: %d of %d genes pass", floor, sum(hi_ok), length(g))
  }

  sig <- de$TFDR[[d$total]][g] < de$fdr & hi_ok
  her_sig <- de$QLFDR[[d$heritable]][g] < de$fdr & sign(her) == sign(tot)
  env_sig <- de$QLFDR[[d$environmental]][g] < de$fdr & sign(env) == sign(tot)
  lo <- d$band_low %||% 0.25; hi_b <- d$band_high %||% 0.75
  fr <- out$heritable_fraction
  out$class <- ifelse(!sig | is.na(sig), "not_changed",
    ifelse(fr > hi_b & her_sig, "heritable",
    ifelse(fr < lo & env_sig, "environmental",
    ifelse(fr >= lo & fr <= hi_b & her_sig & env_sig, "codominant", "undetermined"))))

  # optional site attribution
  if (!is.null(d$reimplant) && nzchar(d$reimplant) && d$reimplant %in% names(de$LFC)) {
    reimp <- de$LFC[[d$reimplant]][g]; reimp_fdr <- de$QLFDR[[d$reimplant]][g]
    gen <- reimp_fdr < de$fdr & sign(reimp) == sign(tot)
    lung <- !gen
    if (!is.null(d$site) && nzchar(d$site) && d$site %in% names(de$LFC)) {
      site <- de$LFC[[d$site]][g]; site_fdr <- de$QLFDR[[d$site]][g]
      lung <- !gen & site_fdr < de$fdr & sign(site) == sign(tot)
    }
    out$site_class <- ifelse(out$class %in% c("environmental", "codominant"),
      ifelse(gen, "in_vivo_general", ifelse(lung, "site_specific", "site_undetermined")), NA)
  }

  write.csv(out, file.path(out_dir, "decomposition.csv"), row.names = FALSE)
  say("[decomposition] class counts:")
  print(table(out$class[out$class != "not_changed"]))
  if ("site_class" %in% names(out)) print(table(out$site_class, useNA = "no"))
  say("  wrote decomposition.csv")
  invisible(out)
}
