# mod_qc_de.R
# Expression filtering (filterByExpr), TMM normalization, per-sample QC, and
# differential expression by edgeR quasi-likelihood F-test plus TREAT against
# a fold-change threshold. Contrasts are whatever the config defines.
# Returns the objects the decomposition module needs.

run_qc_de <- function(cfg, m, samples, out_dir, mask = NULL) {
  need_pkg("edgeR", "qc_de")
  suppressPackageStartupMessages(library(edgeR))
  q <- cfg$qc_de
  grp <- factor(samples$group, levels = unique(samples$group))

  if (length(mask)) {
    drop <- intersect(rownames(m), mask)
    m <- m[setdiff(rownames(m), drop), , drop = FALSE]
    say("[qc_de] applied cross-reactivity mask: removed %d features, %d remain", length(drop), nrow(m))
  }

  y0 <- DGEList(m, group = grp)
  keep <- filterByExpr(y0, group = grp, min.count = q$min_count %||% 10)
  say("[qc_de] filterByExpr kept %d of %d features", sum(keep), nrow(m))
  y <- DGEList(m[keep, , drop = FALSE], group = grp)
  y <- normLibSizes(y, method = "TMM")
  logcpm <- cpm(y, log = TRUE, prior.count = 2)
  lincpm <- cpm(y, normalized.lib.sizes = TRUE)

  write.csv(data.frame(gene = rownames(logcpm), logcpm, check.names = FALSE),
            file.path(out_dir, "TMM_logCPM.csv"), row.names = FALSE)
  write.csv(data.frame(gene = rownames(lincpm), lincpm, check.names = FALSE),
            file.path(out_dir, "TMM_CPM.csv"), row.names = FALSE)

  # per-sample QC
  qc <- data.frame(sample = colnames(y), group = as.character(grp),
                   lib_size = colSums(m[keep, ]), detected_ge10 = colSums(m[keep, ] >= 10),
                   frac_zero = round(colMeans(m[keep, ] == 0), 4),
                   norm_factor = round(y$samples$norm.factors, 4))
  write.csv(qc, file.path(out_dir, "qc_sample_metrics.csv"), row.names = FALSE)
  # replicate correlation
  for (g in levels(grp)) {
    idx <- which(grp == g); if (length(idx) < 2) next
    cc <- cor(logcpm[, idx], method = "spearman")
    say("  %-14s replicate Spearman mean %.3f (min %.3f)", g, mean(cc[lower.tri(cc)]), min(cc[lower.tri(cc)]))
  }

  design <- model.matrix(~ 0 + grp); colnames(design) <- levels(grp)
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  say("  genes: %d | BCV: %.4f", nrow(y), sqrt(y$common.dispersion))

  # parse contrasts of the form "A - B"
  con_defs <- q$contrasts
  if (!length(con_defs)) die("qc_de needs at least one contrast")
  con <- makeContrasts(contrasts = unlist(con_defs), levels = design)
  colnames(con) <- names(con_defs)

  lfc_thr <- log2(q$treat_fold %||% 1.5); fdr <- q$fdr %||% 0.05
  LFC <- list(); QLFDR <- list(); TFDR <- list(); summ <- data.frame()
  de_dir <- ensure_dir(file.path(out_dir, "DE_tables"))
  for (cn in colnames(con)) {
    qlf <- glmQLFTest(fit, contrast = con[, cn])
    tt <- topTags(qlf, n = Inf, sort.by = "none")$table
    tr <- glmTreat(fit, contrast = con[, cn], lfc = lfc_thr)
    tt2 <- topTags(tr, n = Inf, sort.by = "none")$table
    out <- data.frame(gene = rownames(tt), logFC = tt$logFC, logCPM = tt$logCPM,
                      F = tt$F, PValue = tt$PValue, FDR = tt$FDR,
                      TREAT_PValue = tt2[rownames(tt), "PValue"], TREAT_FDR = tt2[rownames(tt), "FDR"])
    write.csv(out[order(out$FDR), ], file.path(de_dir, paste0("DE_", cn, ".csv")), row.names = FALSE)
    LFC[[cn]] <- setNames(out$logFC, out$gene)
    QLFDR[[cn]] <- setNames(out$FDR, out$gene)
    TFDR[[cn]] <- setNames(out$TREAT_FDR, out$gene)
    summ <- rbind(summ, data.frame(contrast = cn, tested = nrow(out),
      QL_sig = sum(out$FDR < fdr), TREAT_sig = sum(out$TREAT_FDR < fdr),
      up = sum(out$TREAT_FDR < fdr & out$logFC > 0), down = sum(out$TREAT_FDR < fdr & out$logFC < 0)))
  }
  write.csv(summ, file.path(out_dir, "DE_summary.csv"), row.names = FALSE)
  say("[qc_de] wrote %d contrast tables and DE_summary.csv", length(con_defs))
  print(summ, row.names = FALSE)

  # group-mean linear CPM for the expressed-floor rule used downstream
  gmean <- sapply(levels(grp), function(g) rowMeans(lincpm[, grp == g, drop = FALSE]))
  list(fit = fit, grp = grp, logcpm = logcpm, lincpm = lincpm, gmean = gmean,
       LFC = LFC, QLFDR = QLFDR, TFDR = TFDR, fdr = fdr, treat_lfc = lfc_thr,
       floor = q$expressed_floor_cpm %||% 0, genes = rownames(logcpm))
}
