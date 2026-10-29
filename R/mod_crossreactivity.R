# mod_crossreactivity.R
# Estimate the host-species contribution to each feature from a dilution
# series, and build an exclusion mask. The estimator is the one validated in
# the source study: counts from a PCR-targeted library do not measure input
# on an absolute scale, so the pure arms cannot be compared directly, but
# within the titration both templates compete in one tube and two disjoint
# calibration sets make the decomposition exact.
#
#   observed(d) = s(d)*pureTarget + m(d)*pureHost
#   s(d)  from host-negative features (zero host signal) -> target scaling
#   m(d)  from host-dominated features (host >> target)  -> host scaling
#   host fraction = m(d)*pureHost / observed(d)
#
# A feature is masked if its host fraction at the operating-point mixture
# exceeds `host_tolerance`, in union with any curated list supplied.

run_crossreactivity <- function(cfg, project_dir, out_dir) {
  x <- cfg$crossreactivity
  path <- resolve_path(x$dilution_counts, project_dir)
  say("[crossreactivity] reading %s", path)
  m <- read_counts(path, collapse = TRUE)

  # Mean CPM per mixture, from the replicate columns that share each prefix.
  mix_cpm <- list(); masses <- data.frame()
  for (mx in x$mixtures) {
    cols <- grep(mx$column_prefix, colnames(m), fixed = TRUE)
    if (!length(cols)) die("no columns match mixture prefix '%s'", mx$column_prefix)
    sub <- m[, cols, drop = FALSE]
    cpm <- sweep(sub, 2, colSums(sub), "/") * 1e6
    mix_cpm[[mx$column_prefix]] <- rowMeans(cpm)
    masses <- rbind(masses, data.frame(name = mx$column_prefix, target = mx$target_ng,
                                       host = mx$host_ng, nrep = length(cols)))
  }
  M <- do.call(cbind, mix_cpm)
  pure_t <- masses$name[masses$host == 0][1]
  pure_h <- masses$name[masses$target == 0][1]
  if (is.na(pure_t) || is.na(pure_h)) die("need one pure-target (host_ng 0) and one pure-host (target_ng 0) mixture")
  Tg <- M[, pure_t]; Ho <- M[, pure_h]
  say("  features: %d | pure-target arm: %s | pure-host arm: %s", nrow(M), pure_t, pure_h)

  dil <- masses[masses$target > 0 & masses$host > 0, ]
  host_frac <- setNames(dil$host / (dil$host + dil$target), dil$name)

  # calibration sets
  tneg <- Ho == 0 & Tg >= 20                 # host-negative, target-expressed -> target scaling
  hdom <- Ho >= 10 & Ho / (Tg + 0.01) >= 20  # host-dominated -> host scaling
  say("  calibration features: %d host-negative, %d host-dominated", sum(tneg), sum(hdom))
  if (sum(tneg) < 50 || sum(hdom) < 10)
    say("  WARNING: few calibration features; host-fraction estimates may be noisy")

  s_d <- sapply(dil$name, function(d) median(M[tneg, d] / Tg[tneg]))
  m_d <- sapply(dil$name, function(d) median(M[hdom, d] / Ho[hdom]))

  contam <- sapply(dil$name, function(d) {
    hp <- s_d[d] * Tg; mp <- m_d[d] * Ho
    ifelse(hp + mp > 0, mp / (hp + mp), 0)
  })
  rownames(contam) <- rownames(M)

  op <- x$operating_point
  if (!op %in% colnames(contam)) die("operating_point '%s' is not a target+host mixture", op)
  expr <- Tg >= (x$min_target_cpm %||% 5)
  tol <- x$host_tolerance %||% 0.10
  mask <- contam[, op] > tol
  curated <- x$curated_exclusions
  if (length(curated)) mask <- mask | (rownames(contam) %in% curated)

  res <- data.frame(feature = rownames(M), cpm_pure_target = round(Tg, 3),
                    cpm_pure_host = round(Ho, 3),
                    host_fraction_op = round(contam[, op], 4),
                    expressed = expr, masked = mask, row.names = NULL)
  write.csv(res, file.path(out_dir, "crossreactivity_per_feature.csv"), row.names = FALSE)
  masked_ids <- rownames(contam)[mask]
  writeLines(masked_ids, file.path(out_dir, "exclusion_mask.txt"))

  say("  operating point %s = %.1f%% host RNA by mass", op, 100 * host_frac[op])
  say("  masked at >%.0f%% host signal: %d features (%d also expressed in target)",
      100 * tol, sum(mask), sum(mask & expr))
  say("  wrote crossreactivity_per_feature.csv and exclusion_mask.txt")
  invisible(masked_ids)
}

`%||%` <- function(a, b) if (is.null(a)) b else a
