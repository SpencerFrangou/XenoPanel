#!/usr/bin/env Rscript
# test_recovery.R -- run the toolkit on the synthetic example and check it
# recovers the planted truth. Exits non-zero on any failure.
here <- tryCatch(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                 error = function(e) ".")
root <- normalizePath(file.path(here, ".."))
ex   <- file.path(root, "examples")

system2("Rscript", c(file.path(ex, "make_example_data.R")), stdout = FALSE)
status <- system2("Rscript", c(file.path(root, "run.R"), file.path(root, "config.example.yaml")),
                  stdout = FALSE, stderr = FALSE)
res <- file.path(ex, "results")
fail <- 0; ok <- function(lab, cond) { cat(sprintf("  [%s] %s\n", ifelse(cond, "PASS", "FAIL"), lab)); if (!cond) fail <<- fail + 1 }

# cross-reactivity: recover the planted host-cross-reactive genes
truth <- readLines(file.path(ex, "TRUTH_host_crossreactive_genes.txt"))
mask  <- readLines(file.path(res, "exclusion_mask.txt"))
sens  <- mean(truth %in% mask)               # of planted host genes, fraction masked
xr <- read.csv(file.path(res, "crossreactivity_per_feature.csv"))
expressed_truth <- intersect(truth, xr$feature[xr$expressed])
sens_expr <- mean(expressed_truth %in% mask)
fdr_mask <- mean(!(mask %in% truth))
cat(sprintf("crossreactivity: %d masked; recall(all)=%.2f recall(expressed)=%.2f falsemask=%.2f\n",
            length(mask), sens, sens_expr, fdr_mask))
ok("masks most expressed host-cross-reactive genes (recall>=0.9)", sens_expr >= 0.9)
ok("few non-host genes masked (<0.15)", fdr_mask < 0.15)

# decomposition: planted genes 1-150 heritable-up, 151-300 environmental-up
dec <- read.csv(file.path(res, "decomposition.csv"))
her_truth <- sprintf("GENE%04d", 1:150); env_truth <- sprintf("GENE%04d", 151:300)
cls <- setNames(dec$class, dec$gene)
her_recall <- mean(cls[intersect(her_truth, dec$gene)] == "heritable", na.rm = TRUE)
env_recall <- mean(cls[intersect(env_truth, dec$gene)] == "environmental", na.rm = TRUE)
cat(sprintf("decomposition: heritable recall=%.2f environmental recall=%.2f max|resid|=%.1e\n",
            her_recall, env_recall, max(abs(dec$residual), na.rm = TRUE)))
ok("heritable genes classed heritable (recall>=0.8)", her_recall >= 0.8)
ok("environmental genes classed environmental (recall>=0.8)", env_recall >= 0.8)
ok("additive identity holds (max|residual|<1e-6)", max(abs(dec$residual), na.rm = TRUE) < 1e-6)

# cross-contamination between the two classes should be low
her_as_env <- mean(cls[intersect(her_truth, dec$gene)] == "environmental", na.rm = TRUE)
ok("heritable rarely misclassed environmental (<0.1)", her_as_env < 0.1)

# entropy outputs exist
ok("diversity table written", file.exists(file.path(res, "diversity_per_sample.csv")))
ok("network entropy table written", file.exists(file.path(res, "network_entropy.csv")))
ne <- read.csv(file.path(res, "network_entropy.csv"))
ok("CCAT computed for every sample", all(is.finite(ne$CCAT)))

cat(sprintf("\n%s (%d failure%s)\n", ifelse(fail == 0, "ALL TESTS PASSED", "TESTS FAILED"),
            fail, ifelse(fail == 1, "", "s")))
quit(status = if (fail == 0) 0 else 1)
