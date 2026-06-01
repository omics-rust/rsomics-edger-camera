#!/usr/bin/env Rscript
# camera oracle: read a log-expression matrix TSV, a design TSV, and a GMT, run
# limma camera, and write GeneSet,(Correlation,)Direction,PValue,FDR matching
# rsomics-edger-camera. NGenes is dropped from the diff (it is an integer count
# both sides agree on by construction).
#
# Usage: camera_oracle.R <expr.tsv> <design.tsv> <sets.gmt> <coef> <out.tsv> [mode]
#   mode = "fixed" (default, inter.gene.cor=0.01) | "estimate" (inter.gene.cor=NA)
suppressMessages(library(limma))

args <- commandArgs(trailingOnly = TRUE)
expr_path <- args[1]
design_path <- args[2]
gmt_path <- args[3]
coef <- as.integer(args[4])
out_path <- args[5]
mode <- if (length(args) >= 6) args[6] else "fixed"

E <- as.matrix(read.delim(expr_path, row.names = 1, check.names = FALSE))
design <- as.matrix(read.delim(design_path, row.names = 1, check.names = FALSE))

read_gmt <- function(path) {
  lines <- readLines(path)
  lines <- lines[nchar(trimws(lines)) > 0]
  sets <- lapply(lines, function(l) {
    f <- strsplit(l, "\t")[[1]]
    unique(f[-(1:2)])
  })
  names(sets) <- vapply(lines, function(l) strsplit(l, "\t")[[1]][1], "")
  sets
}
sets <- read_gmt(gmt_path)
idx <- lapply(sets, function(g) which(rownames(E) %in% g))

if (mode == "estimate") {
  res <- camera(E, idx, design, contrast = coef, inter.gene.cor = NA)
} else {
  res <- camera(E, idx, design, contrast = coef)
}

con <- file(out_path, "w")
g <- function(x) formatC(x, digits = 10, format = "g", flag = "")
has_cor <- "Correlation" %in% colnames(res)
if (has_cor) {
  writeLines("GeneSet\tNGenes\tCorrelation\tDirection\tPValue\tFDR", con)
} else {
  writeLines("GeneSet\tNGenes\tDirection\tPValue\tFDR", con)
}
for (i in seq_len(nrow(res))) {
  fields <- c(rownames(res)[i], as.character(res$NGenes[i]))
  if (has_cor) fields <- c(fields, g(res$Correlation[i]))
  fields <- c(fields, res$Direction[i], g(res$PValue[i]), g(res$FDR[i]))
  writeLines(paste(fields, collapse = "\t"), con)
}
close(con)
