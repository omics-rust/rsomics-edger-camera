# rsomics-edger-camera

CAMERA competitive gene-set test — given a log-expression matrix, a design, and
a collection of gene sets, ask whether each set is ranked highly in a
differential-expression analysis relative to the other genes, correcting the
set-mean variance for inter-gene correlation. A single-binary Rust
reimplementation of limma's `camera`.

Each gene is scored by its empirical-Bayes moderated *t* (lmFit + eBayes),
mapped to a standard-normal *z*. For a set of *m* genes out of *G*, the variance
inflation factor `VIF = 1 + (m-1)ρ̄` widens the set-mean standard error; a
modified two-sample test then compares the set's mean *z* against the rest.

## Usage

```
rsomics-edger-camera expr.tsv --design design.tsv --gene-sets sets.gmt [--coef N | --contrast c.tsv] [-o camera.tsv]
```

- `expr.tsv` — header row of sample ids, first column gene ids, log-expression values.
- `--design` — header row of coefficient names, first column sample ids (the model matrix).
- `--gene-sets` — GMT: set name, a description (ignored), then gene ids, tab-separated.
- `--coef N` — 1-based coefficient (or contrast) to rank on; defaults to the last column.
- `--contrast` — optional contrast matrix; applies `contrasts.fit` before ranking.
- `--inter-gene-cor ρ` — preset inter-gene correlation shared by all sets (default 0.01, limma's default).
- `--estimate-cor` — estimate the correlation per set instead of using the preset.
- `--allow-neg-cor` — with `--estimate-cor`, allow a negative correlation (VIF below 1).
- `--no-sort` — keep input set order instead of sorting by PValue.

Output columns: `GeneSet NGenes [Correlation] Direction PValue FDR`, where the
`Correlation` column appears only when the correlation is estimated. With a
preset correlation the test uses the full-data degrees of freedom (`G-2`); with
an estimated correlation it uses the residual degrees of freedom.

```
rsomics-edger-camera E.tsv --design design.tsv --gene-sets sets.gmt --coef 2 -o camera.tsv
rsomics-edger-camera E.tsv --design design.tsv --gene-sets sets.gmt --estimate-cor > camera.tsv
```

## Origin

This crate is an independent Rust reimplementation of limma's `camera` based on:

- The published method: Wu, D. and Smyth, G.K. (2012), "Camera: a competitive
  gene set test accounting for inter-gene correlation", Nucleic Acids Research
  40(17):e133, doi:10.1093/nar/gks461 — the variance inflation factor, the
  QR-rotated residual estimator of the inter-gene correlation, the moderated-*t*
  to *z* conversion, and the modified two-sample competitive test.
- The empirical-Bayes moderation it builds on: Smyth, G.K. (2004), Statistical
  Applications in Genetics and Molecular Biology 3(1):3,
  doi:10.2202/1544-6115.1027.
- Black-box behaviour testing against the limma binary via an R oracle
  (`Rscript` + limma), diffed field-by-field in `tests/compat.rs`.

No source code from limma (GPL) was used as reference during implementation.
Output is value-exact against limma `camera` (relative deviation < 1e-6 for
Correlation, PValue, and FDR; Direction matches exactly) for both the preset and
the per-set estimated correlation modes.

License: MIT OR Apache-2.0.
Upstream credit: limma (https://bioconductor.org/packages/limma/), GPL (>=2).
