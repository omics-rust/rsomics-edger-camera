#!/usr/bin/env Rscript
# Generate a log-expression matrix, a two-group design, and a GMT of gene sets
# for compat/perf fixtures. Some sets carry a real DE signal (up or down), some
# carry a shared latent factor (positive inter-gene correlation), the rest are
# random background draws.
#
# Usage: mkfixture.R <ngenes> <nsamples_per_group> <nsets> <setsize> \
#                    <expr_out> <design_out> <gmt_out> [seed]
args <- commandArgs(trailingOnly = TRUE)
ng <- as.integer(args[1])
nper <- as.integer(args[2])
nsets <- as.integer(args[3])
setsize <- as.integer(args[4])
expr_out <- args[5]
design_out <- args[6]
gmt_out <- args[7]
seed <- if (length(args) >= 8) as.integer(args[8]) else 1L
set.seed(seed)

n <- 2L * nper
group <- c(rep(0L, nper), rep(1L, nper))
mu <- rnorm(ng, mean = 6, sd = 2)
sd_gene <- sqrt(0.2 + rgamma(ng, shape = 2, rate = 4))
effect <- rnorm(ng, 0, 0.15)

E <- matrix(0, nrow = ng, ncol = n)
for (i in seq_len(ng)) {
  E[i, ] <- mu[i] + effect[i] * group + rnorm(n, 0, sd_gene[i])
}
rownames(E) <- sprintf("g%06d", seq_len(ng))
colnames(E) <- sprintf("s%03d", seq_len(n))

gene_ids <- rownames(E)
sets <- vector("list", nsets)
for (s in seq_len(nsets)) {
  members <- sample(gene_ids, setsize)
  rows <- match(members, gene_ids)
  kind <- s %% 4
  if (kind == 1) {
    E[rows, ] <- E[rows, ] + outer(rep(1, setsize), group) * 1.0
  } else if (kind == 2) {
    E[rows, ] <- E[rows, ] - outer(rep(1, setsize), group) * 1.0
  } else if (kind == 3) {
    latent <- rnorm(n) * 0.7
    E[rows, ] <- E[rows, ] + matrix(latent, setsize, n, byrow = TRUE)
  }
  sets[[s]] <- members
}
names(sets) <- sprintf("set%05d", seq_len(nsets))

design <- cbind(Intercept = 1, group = group)
rownames(design) <- colnames(E)

write_tsv <- function(mat, path, corner) {
  con <- file(path, "w")
  writeLines(paste(c(corner, colnames(mat)), collapse = "\t"), con)
  g <- function(x) formatC(x, digits = 10, format = "g", flag = "")
  for (i in seq_len(nrow(mat))) {
    writeLines(paste(c(rownames(mat)[i], g(mat[i, ])), collapse = "\t"), con)
  }
  close(con)
}
write_tsv(E, expr_out, "gene")
write_tsv(design, design_out, "sample")

con <- file(gmt_out, "w")
for (s in seq_len(nsets)) {
  writeLines(paste(c(names(sets)[s], "na", sets[[s]]), collapse = "\t"), con)
}
close(con)
