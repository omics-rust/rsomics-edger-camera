//! limma camera: a competitive gene-set test that asks whether a set ranks high
//! in a differential-expression analysis, correcting the set-mean variance for
//! inter-gene correlation (Wu & Smyth 2012, NAR, doi:10.1093/nar/gks461).
//!
//! Clean-room: the method follows the published paper and is validated
//! black-box against the limma binary. No limma (GPL) source was consulted.

mod camera;
mod ebayes;
mod fit;
mod matrix;
mod special;

use std::collections::HashMap;
use std::io::{BufWriter, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

use camera::{Cor, Engine, bh_adjust};
pub use matrix::{read_contrast, read_design, read_expr, read_gmt};

pub struct Options<'a> {
    pub expr: &'a Path,
    pub design: &'a Path,
    pub gene_sets: &'a Path,
    pub contrast: Option<&'a Path>,
    /// 1-based coefficient (or contrast) to rank on; defaults to the last.
    pub coef: Option<usize>,
    /// fixed preset inter-gene correlation; None estimates it per set.
    pub inter_gene_cor: Option<f64>,
    pub allow_neg_cor: bool,
    pub sort: bool,
}

pub struct Row {
    pub name: String,
    pub n_genes: usize,
    pub correlation: Option<f64>,
    pub direction: &'static str,
    pub p_value: f64,
    pub fdr: f64,
}

pub struct Results {
    pub coef_name: String,
    pub rows: Vec<Row>,
    /// true when an estimated (per-set) correlation column is present.
    pub has_correlation: bool,
}

pub fn run(opts: &Options) -> Result<Results> {
    let expr = read_expr(opts.expr)?;
    let design = read_design(opts.design)?;
    if design.x.len() != expr.samples.len() {
        return Err(RsomicsError::InvalidInput(format!(
            "design has {} rows, expression has {} samples",
            design.x.len(),
            expr.samples.len()
        )));
    }

    let (base, qr) = fit::lm_fit(&expr.y, &expr.genes, &design.x, &design.coef_names)?;
    let fit = if let Some(cpath) = opts.contrast {
        let contrast = read_contrast(cpath, &design.coef_names)?;
        let xtx_inv = qr.xtx_inv();
        fit::contrasts_fit(&base, &contrast, &xtx_inv)
    } else {
        base
    };

    let m = ebayes::ebayes(&fit);

    let nc = fit.coef_names.len();
    let coef = opts.coef.unwrap_or(nc);
    if coef < 1 || coef > nc {
        return Err(RsomicsError::InvalidInput(format!(
            "--coef {coef} out of range 1..={nc}"
        )));
    }
    let k = coef - 1;
    let modt: Vec<f64> = (0..fit.genes.len()).map(|gi| m.t[gi][k]).collect();

    let engine = Engine::new(&modt, m.df_total, &qr, &expr.y);

    let gene_index: HashMap<&str, usize> = expr
        .genes
        .iter()
        .enumerate()
        .map(|(i, g)| (g.as_str(), i))
        .collect();

    let cor = match opts.inter_gene_cor {
        Some(rho) => Cor::Preset(rho),
        None => Cor::Estimate {
            allow_neg: opts.allow_neg_cor,
        },
    };
    let has_correlation = opts.inter_gene_cor.is_none();

    let sets = read_gmt(opts.gene_sets)?;
    let mut results = Vec::with_capacity(sets.len());
    for set in &sets {
        let idx: Vec<usize> = set
            .members
            .iter()
            .filter_map(|g| gene_index.get(g.as_str()).copied())
            .collect();
        if idx.len() < 2 || idx.len() >= expr.genes.len() - 1 {
            continue;
        }
        results.push(engine.test(&set.name, &idx, &cor));
    }
    if results.is_empty() {
        return Err(RsomicsError::InvalidInput(
            "no gene set had between 2 and G-2 genes present in the matrix".into(),
        ));
    }

    let pvals: Vec<f64> = results.iter().map(|r| r.p_value).collect();
    let fdr = bh_adjust(&pvals);

    let mut order: Vec<usize> = (0..results.len()).collect();
    if opts.sort {
        order.sort_by(|&a, &b| results[a].p_value.partial_cmp(&results[b].p_value).unwrap());
    }

    let rows = order
        .into_iter()
        .map(|i| Row {
            name: results[i].name.clone(),
            n_genes: results[i].n_genes,
            correlation: results[i].correlation,
            direction: if results[i].up { "Up" } else { "Down" },
            p_value: results[i].p_value,
            fdr: fdr[i],
        })
        .collect();

    Ok(Results {
        coef_name: fit.coef_names[k].clone(),
        rows,
        has_correlation,
    })
}

pub fn write_results(res: &Results, out: &mut dyn Write) -> Result<()> {
    let mut w = BufWriter::with_capacity(1 << 20, out);
    if res.has_correlation {
        writeln!(w, "GeneSet\tNGenes\tCorrelation\tDirection\tPValue\tFDR")
            .map_err(RsomicsError::Io)?;
    } else {
        writeln!(w, "GeneSet\tNGenes\tDirection\tPValue\tFDR").map_err(RsomicsError::Io)?;
    }
    let mut fmt = ryu::Buffer::new();
    let mut line = String::with_capacity(128);
    let mut ni = itoa::Buffer::new();
    for row in &res.rows {
        line.clear();
        line.push_str(&row.name);
        line.push('\t');
        line.push_str(ni.format(row.n_genes));
        if let Some(c) = row.correlation {
            line.push('\t');
            line.push_str(fmt.format(c));
        }
        line.push('\t');
        line.push_str(row.direction);
        line.push('\t');
        line.push_str(fmt.format(row.p_value));
        line.push('\t');
        line.push_str(fmt.format(row.fdr));
        line.push('\n');
        w.write_all(line.as_bytes()).map_err(RsomicsError::Io)?;
    }
    w.flush().map_err(RsomicsError::Io)?;
    Ok(())
}
