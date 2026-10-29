#!/usr/bin/env Rscript
# make_example_data.R -- generate synthetic, fully simulated inputs so the
# toolkit can be run end to end with no real data. Nothing here comes from any
# experiment; the numbers are drawn from a negative binomial with a planted
# five-state structure and a planted host cross-reactivity pattern, so the
# modules have something true to recover.
set.seed(1)
here <- tryCatch(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                 error = function(e) ".")
if (!length(here) || !nzchar(here)) here <- "."

G <- 2000L                                   # features
genes <- sprintf("GENE%04d", seq_len(G))
states <- c("Parental", "PrimaryTumor", "MetLung", "MetSubQ", "MetLungCL")
reps <- 3L
samples <- expand.grid(rep = seq_len(reps), group = states)
samples$sample <- sprintf("%s_%d", samples$group, samples$rep)
context <- c(Parental = "in_vitro", PrimaryTumor = "in_vivo", MetLung = "in_vivo",
             MetSubQ = "in_vivo", MetLungCL = "in_vitro")

base <- 2^runif(G, 2, 12)                     # baseline expression per gene
# planted effects (log2) relative to Parental, per state, for two gene blocks:
#  genes 1-150  heritable up    (present already in the cultured line)
#  genes 151-300 environmental up (only in tissue, lost on explant)
eff <- matrix(0, G, length(states), dimnames = list(genes, states))
her <- 1:150; env <- 151:300
eff[her, c("PrimaryTumor","MetLung","MetSubQ","MetLungCL")] <- 2.0
eff[env, c("PrimaryTumor","MetLung","MetSubQ")] <- 2.5          # requires tissue; zero in MetLungCL
mu <- base * 2^eff

counts <- matrix(0L, G, nrow(samples), dimnames = list(genes, samples$sample))
size <- 20                                    # NB dispersion ~ 1/size
for (j in seq_len(nrow(samples))) {
  lambda <- mu[, as.character(samples$group[j])]
  counts[, j] <- rnbinom(G, mu = lambda, size = size)
}
write.csv(data.frame(feature = genes, counts, check.names = FALSE),
          file.path(here, "counts.csv"), row.names = FALSE)
write.csv(data.frame(sample = samples$sample, group = as.character(samples$group),
                     context = context[as.character(samples$group)], replicate = samples$rep),
          file.path(here, "samples.csv"), row.names = FALSE)

# ---- dilution series with a planted host pattern ---------------------------
# 150 features are "host cross-reactive": they carry a constant host template
# that shows up when target RNA is diluted. The rest respond only to target.
host_genes <- sample(genes, 150)
host_level <- setNames(rep(0, G), genes); host_level[host_genes] <- 2^runif(150, 6, 12)
mix <- data.frame(name = c("T50_H0","T10_H50","T2.5_H50","T0.64_H50","T0_H50"),
                  target = c(50,10,2.5,0.64,0), host = c(0,50,50,50,50))
dil <- data.frame(feature = genes)
for (i in seq_len(nrow(mix))) {
  tgt <- base * (mix$target[i] / 50)          # target signal scales with target mass
  hst <- host_level * (mix$host[i] / 50)      # host signal scales with host mass
  for (r in 1:3) {
    lambda <- pmax(tgt + hst, 0.01)
    dil[[sprintf("%s_R%d", mix$name[i], r)]] <- rnbinom(G, mu = lambda, size = size)
  }
}
write.csv(dil, file.path(here, "dilution_series.csv"), row.names = FALSE)

writeLines(host_genes, file.path(here, "TRUTH_host_crossreactive_genes.txt"))
cat("wrote counts.csv, samples.csv, dilution_series.csv and TRUTH_host_crossreactive_genes.txt\n")
cat(sprintf("planted: %d heritable-up, %d environmental-up, %d host cross-reactive\n",
            length(her), length(env), length(host_genes)))

# ---- synthetic PPI edge list (so network_entropy/CCAT has a network) -------
# A random geometric-ish graph over the simulated genes. Not biology; just a
# connected network so SR/CCAT have something to run on.
ne <- 6000L
e1 <- sample(genes, ne, replace = TRUE); e2 <- sample(genes, ne, replace = TRUE)
keep <- e1 != e2
edges <- unique(data.frame(a = e1[keep], b = e2[keep]))
write.table(edges, file.path(here, "ppi_edges.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# ---- synthetic gene sets (GMT) so signature scoring runs -------------------
gmt <- c(
  paste(c("HERITABLE_UP_SET", "synthetic", genes[her]), collapse = "\t"),
  paste(c("ENVIRONMENTAL_UP_SET", "synthetic", genes[env]), collapse = "\t"),
  paste(c("RANDOM_SET", "synthetic", sample(genes, 100)), collapse = "\t"))
writeLines(gmt, file.path(here, "gene_sets.gmt"))
cat("wrote ppi_edges.tsv and gene_sets.gmt\n")
