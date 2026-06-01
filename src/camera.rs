//! CAMERA competitive gene-set test (Wu & Smyth 2012, NAR, doi:10.1093/nar/gks461).
//!
//! Genes are scored by their moderated-t mapped to standard-normal z. For each
//! set, the inter-gene correlation inflates the variance of the set mean: the
//! modified two-sample test compares the set's mean z to the rest, with the
//! set's standard error scaled by sqrt(VIF). A fixed preset correlation uses
//! the full-data df (G-2); an estimated per-set correlation uses the residual
//! df and floors the VIF at 1 unless negative correlation is allowed.

use crate::fit::Qr;
use crate::special::{t_pvalue_two_sided, t_to_z};

pub struct SetResult {
    pub name: String,
    pub n_genes: usize,
    /// estimated inter-gene correlation; None when a preset value was used.
    pub correlation: Option<f64>,
    pub up: bool,
    pub p_value: f64,
}

/// limma interGeneCorrelation: VIF = m · Σ_k (ū·k)² / d over the row-standardized
/// trailing residual rotation Q2'y; correlation solves VIF = 1 + (m-1)ρ.
pub fn inter_gene_correlation(set_rows: &[&[f64]], qr: &Qr) -> (f64, f64) {
    let (n, p) = qr.dim();
    let d = n - p;
    let m = set_rows.len();
    let mut col_sums = vec![0.0f64; d];
    let mut u = vec![0.0f64; d];
    for &row in set_rows {
        qr.rotated_residual(row, &mut u);
        let ss: f64 = u.iter().map(|&v| v * v).sum();
        let scale = (ss / d as f64).sqrt();
        for (cs, &v) in col_sums.iter_mut().zip(&u) {
            *cs += v / scale;
        }
    }
    let mean_sq: f64 = col_sums.iter().map(|&s| (s / m as f64).powi(2)).sum();
    let vif = m as f64 * mean_sq / d as f64;
    let correlation = (vif - 1.0) / (m as f64 - 1.0);
    (vif, correlation)
}

pub enum Cor {
    /// fixed preset shared by every set
    Preset(f64),
    /// estimate per set; bool = allow negative (VIF below 1)
    Estimate { allow_neg: bool },
}

/// Run the competitive test for one set given the genewise moderated t (column
/// already selected), the QR for residual rotation, the expression rows in gene
/// order, and a mapping from set member id to row index.
pub struct Engine<'a> {
    pub z: Vec<f64>,
    qr: &'a Qr,
    expr: &'a [Vec<f64>],
}

impl<'a> Engine<'a> {
    pub fn new(modt: &[f64], df_total: f64, qr: &'a Qr, expr: &'a [Vec<f64>]) -> Engine<'a> {
        let z = modt.iter().map(|&t| t_to_z(t, df_total)).collect();
        Engine { z, qr, expr }
    }

    pub fn test(&self, name: &str, idx: &[usize], cor: &Cor) -> SetResult {
        let g = self.z.len();
        let m1 = idx.len();
        let m2 = g - m1;

        let in_set: f64 = idx.iter().map(|&i| self.z[i]).sum();
        let zbar1 = in_set / m1 as f64;
        let total: f64 = self.z.iter().sum();
        let zbar2 = (total - in_set) / m2 as f64;

        let mut in_mask = vec![false; g];
        for &i in idx {
            in_mask[i] = true;
        }
        let mut ss = 0.0;
        for (i, &zi) in self.z.iter().enumerate() {
            let mean = if in_mask[i] { zbar1 } else { zbar2 };
            ss += (zi - mean).powi(2);
        }
        let sigma2 = ss / (g as f64 - 2.0);

        let (correlation, vif_used, df) = match cor {
            Cor::Preset(rho) => {
                let vif = 1.0 + (m1 as f64 - 1.0) * rho;
                (None, vif, g as f64 - 2.0)
            }
            Cor::Estimate { allow_neg } => {
                let rows: Vec<&[f64]> = idx.iter().map(|&i| self.expr[i].as_slice()).collect();
                let (vif, rho) = inter_gene_correlation(&rows, self.qr);
                let used = if *allow_neg { vif } else { vif.max(1.0) };
                let (n, p) = self.qr.dim();
                (Some(rho), used, (n - p) as f64)
            }
        };

        let se = (sigma2 * (vif_used / m1 as f64 + 1.0 / m2 as f64)).sqrt();
        let t = (zbar1 - zbar2) / se;
        let p_value = t_pvalue_two_sided(t, df);

        SetResult {
            name: name.to_string(),
            n_genes: m1,
            correlation,
            up: zbar1 > zbar2,
            p_value,
        }
    }
}

/// Benjamini-Hochberg adjusted p-values, returned in input order.
pub fn bh_adjust(p: &[f64]) -> Vec<f64> {
    let n = p.len();
    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_by(|&a, &b| p[b].partial_cmp(&p[a]).unwrap());
    let mut adj = vec![0.0; n];
    let mut cummin = f64::INFINITY;
    for (rank, &i) in idx.iter().enumerate() {
        let m = (n - rank) as f64;
        let v = (n as f64 / m * p[i]).min(1.0);
        cummin = cummin.min(v);
        adj[i] = cummin;
    }
    adj
}
