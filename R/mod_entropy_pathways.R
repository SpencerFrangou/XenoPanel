# mod_entropy_pathways.R
# Four optional read-outs on the normalized matrix. Each sub-part checks for
# the resources and packages it needs and is skipped with a note if they are
# absent, so the module runs with whatever is available.

shannon <- function(x) { p <- x[x > 0] / sum(x); -sum(p * log2(p)) }
gini    <- function(x) { x <- sort(as.numeric(x)); n <- length(x); sum((2 * seq_len(n) - n - 1) * x) / (n * sum(x)) }
n50     <- function(x) { s <- sort(x, decreasing = TRUE); which(cumsum(s) >= 0.5 * sum(s))[1] }

# CCAT: Pearson correlation between a sample's expression and node degree,
# over genes shared by the matrix and the network. Implemented directly (the
# single published definition) so it runs on panels below SCENT's 5,000-gene
# guard. Expression is log-transformed if it is not already.
ccat_vec <- function(E, deg) {
  if (max(E, na.rm = TRUE) > 100) E <- log2(E + 1)
  common <- intersect(names(deg), rownames(E))
  k <- deg[common]
  as.vector(cor(E[common, , drop = FALSE], k))
}

read_edges <- function(path) {
  e <- read.table(path, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                  col.names = c("a", "b"), fill = TRUE)
  e[nzchar(e$a) & nzchar(e$b) & e$a != e$b, ]
}

run_entropy_pathways <- function(cfg, de, samples, out_dir, project_dir) {
  e <- cfg$entropy_pathways
  grp <- de$grp
  lincpm <- de$lincpm; logcpm <- de$logcpm

  # ---- diversity ----------------------------------------------------------
  if (isTRUE(e$diversity)) {
    div <- data.frame(sample = colnames(lincpm), group = as.character(grp),
                      H_bits = apply(lincpm, 2, shannon),
                      evenness = apply(lincpm, 2, function(x) shannon(x) / log2(sum(x > 0))),
                      gini = apply(lincpm, 2, gini),
                      N50_genes = apply(lincpm, 2, n50))
    div$hill_q1 <- 2^div$H_bits
    write.csv(div, file.path(out_dir, "diversity_per_sample.csv"), row.names = FALSE)
    say("[entropy] diversity: wrote diversity_per_sample.csv")
  }

  # ---- network entropy: SR (SCENT if present) and CCAT --------------------
  if (isTRUE(e$network_entropy)) {
    ppi <- resolve_path(e$ppi_edges, project_dir)
    if (is.null(ppi) || !file.exists(ppi)) {
      say("[entropy] network_entropy skipped: no ppi_edges file")
    } else if (!have("igraph")) {
      say("[entropy] network_entropy skipped: package 'igraph' not installed")
    } else {
      suppressPackageStartupMessages(library(igraph))
      ed <- read_edges(ppi)
      G <- simplify(graph_from_data_frame(ed, directed = FALSE))
      Gp <- induced_subgraph(G, intersect(rownames(logcpm), V(G)$name))
      comp <- components(Gp); Gm <- induced_subgraph(Gp, V(Gp)$name[comp$membership == which.max(comp$csize)])
      deg <- setNames(degree(Gm), V(Gm)$name)
      say("[entropy] network: %d genes in the largest connected component", length(deg))
      ccat <- ccat_vec(logcpm, deg)
      res <- data.frame(sample = colnames(logcpm), group = as.character(grp), CCAT = ccat)
      if (have("SCENT")) {
        suppressPackageStartupMessages(library(SCENT))
        A <- as.matrix(as_adjacency_matrix(Gm)); X <- logcpm[rownames(A), ]
        sr <- tryCatch(as.numeric(CompSRana(list(expMC = X, adjMC = A), local = FALSE)$SR),
                       error = function(err) { say("  SR failed: %s", conditionMessage(err)); NA })
        res$SR <- sr
      } else say("  SR skipped: package 'SCENT' not installed (CCAT still computed)")
      write.csv(res, file.path(out_dir, "network_entropy.csv"), row.names = FALSE)
      say("  wrote network_entropy.csv")
    }
  }

  # ---- signature scoring --------------------------------------------------
  if (isTRUE(e$signatures)) {
    gmt <- resolve_path(e$gene_sets, project_dir)
    if (is.null(gmt) || !file.exists(gmt)) say("[entropy] signatures skipped: no gene_sets GMT file")
    else run_signatures(gmt, logcpm, grp, out_dir)
  }

  # ---- PROGENy ------------------------------------------------------------
  if (isTRUE(e$progeny))
    tryCatch(run_progeny(logcpm, grp, out_dir),
             error = function(err) say("[entropy] progeny skipped: %s", conditionMessage(err)))

  # ---- CollecTRI TF activity ---------------------------------------------
  if (isTRUE(e$tf_activity)) {
    net <- resolve_path(e$collectri, project_dir)
    if (is.null(net) || !file.exists(net)) say("[entropy] tf_activity skipped: no collectri network file")
    else tryCatch(run_tf_activity(net, logcpm, grp, e$min_regulon_size %||% 5, out_dir),
                  error = function(err) say("[entropy] tf_activity skipped: %s", conditionMessage(err)))
  }
}

read_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE)
  sets <- lapply(strsplit(lines, "\t"), function(f) f[-c(1, 2)])
  names(sets) <- sapply(strsplit(lines, "\t"), `[`, 1)
  lapply(sets, function(g) g[nzchar(g)])
}

run_signatures <- function(gmt, logcpm, grp, out_dir) {
  sets <- read_gmt(gmt)
  sets <- lapply(sets, function(g) intersect(g, rownames(logcpm)))
  sets <- sets[sapply(sets, length) >= 5]
  if (!length(sets)) { say("[entropy] signatures skipped: no gene set has >=5 genes on the panel"); return() }
  if (have("singscore") && have("SummarizedExperiment")) {
    suppressPackageStartupMessages({library(singscore); library(SummarizedExperiment)})
    rk <- rankGenes(logcpm)
    sc <- t(sapply(sets, function(s) simpleScore(rk, upSet = s)$TotalScore))
    colnames(sc) <- colnames(logcpm)
    write.csv(data.frame(signature = rownames(sc), sc, check.names = FALSE),
              file.path(out_dir, "signature_singscore.csv"), row.names = FALSE)
    say("[entropy] signatures: scored %d sets with singscore", length(sets))
  } else say("[entropy] signatures skipped: package 'singscore' not installed")
}

run_progeny <- function(logcpm, grp, out_dir) {
  if (!have("progeny") || !have("decoupleR")) { say("[entropy] progeny skipped: needs 'progeny' and 'decoupleR'"); return() }
  suppressPackageStartupMessages({library(progeny); library(decoupleR)})
  mod <- progeny::getModel("Human", top = 500)
  net <- data.frame(source = rep(colnames(mod), each = nrow(mod)),
                    target = rep(rownames(mod), times = ncol(mod)),
                    weight = as.vector(as.matrix(mod)))
  net <- net[net$weight != 0 & net$target %in% rownames(logcpm), ]
  keep <- names(which(table(net$source) >= 30))
  net <- net[net$source %in% keep, ]
  acts <- decoupleR::run_mlm(mat = logcpm, network = net, .source = "source", .target = "target", .mor = "weight", minsize = 30)
  acts <- acts[acts$statistic == "mlm", ]
  w <- reshape(as.data.frame(acts[, c("source", "condition", "score")]), idvar = "source", timevar = "condition", direction = "wide")
  colnames(w) <- sub("^score\\.", "", colnames(w))
  write.csv(w, file.path(out_dir, "progeny_activity.csv"), row.names = FALSE)
  say("[entropy] progeny: scored %d pathways", nrow(w))
}

run_tf_activity <- function(path, logcpm, grp, minsize, out_dir) {
  if (!have("decoupleR")) { say("[entropy] tf_activity skipped: needs 'decoupleR'"); return() }
  suppressPackageStartupMessages(library(decoupleR))
  x <- read.csv(path, stringsAsFactors = FALSE)
  if (!all(c("source", "target", "weight") %in% names(x)))
    die("collectri file needs columns source,target,weight")
  net <- unique(data.frame(source = x$source, target = x$target, mor = as.numeric(x$weight)))
  net <- net[net$target %in% rownames(logcpm) & is.finite(net$mor), ]
  keep <- names(which(table(net$source) >= minsize)); net <- net[net$source %in% keep, ]
  if (!nrow(net)) { say("[entropy] tf_activity skipped: no regulon has >=%d targets on the panel", minsize); return() }
  acts <- run_ulm(mat = logcpm, net = net, .source = "source", .target = "target", .mor = "mor", minsize = minsize)
  w <- reshape(as.data.frame(acts[, c("source", "condition", "score")]), idvar = "source", timevar = "condition", direction = "wide")
  colnames(w) <- sub("^score\\.", "", colnames(w))
  write.csv(w, file.path(out_dir, "tf_activity.csv"), row.names = FALSE)
  say("[entropy] tf_activity: scored %d regulons", nrow(w))
}
