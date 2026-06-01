//! Per-gene least-squares fit (limma lmFit, method="ls") and contrasts.fit.
//!
//! Householder QR of the design X once, then solve every gene. The QR also
//! yields R for (X'X)^-1 and the trailing rotated residuals Q2'y that
//! interGeneCorrelation needs.

use rsomics_common::{Result, RsomicsError};

pub struct Qr {
    n: usize,
    p: usize,
    qr: Vec<Vec<f64>>,
    rdiag: Vec<f64>,
}

impl Qr {
    pub fn new(x: &[Vec<f64>]) -> Result<Qr> {
        let n = x.len();
        let p = x[0].len();
        if n < p {
            return Err(RsomicsError::InvalidInput(format!(
                "design has {n} samples < {p} coefficients (rank-deficient)"
            )));
        }
        let mut qr: Vec<Vec<f64>> = x.to_vec();
        let mut rdiag = vec![0.0; p];
        for k in 0..p {
            let mut nrm = 0.0f64;
            for row in qr.iter().take(n).skip(k) {
                nrm = nrm.hypot(row[k]);
            }
            if nrm == 0.0 {
                return Err(RsomicsError::InvalidInput(
                    "design matrix is rank-deficient".into(),
                ));
            }
            if qr[k][k] < 0.0 {
                nrm = -nrm;
            }
            for row in qr.iter_mut().take(n).skip(k) {
                row[k] /= nrm;
            }
            qr[k][k] += 1.0;
            for j in (k + 1)..p {
                let mut s = 0.0;
                for row in qr.iter().take(n).skip(k) {
                    s += row[k] * row[j];
                }
                s = -s / qr[k][k];
                for row in qr.iter_mut().take(n).skip(k) {
                    let add = s * row[k];
                    row[j] += add;
                }
            }
            rdiag[k] = -nrm;
        }
        Ok(Qr { n, p, qr, rdiag })
    }

    pub fn dim(&self) -> (usize, usize) {
        (self.n, self.p)
    }

    /// Apply Q' to y in place (length n).
    #[allow(clippy::needless_range_loop)]
    fn qty(&self, y: &mut [f64]) {
        for k in 0..self.p {
            let mut s = 0.0;
            for i in k..self.n {
                s += self.qr[i][k] * y[i];
            }
            s = -s / self.qr[k][k];
            for i in k..self.n {
                y[i] += s * self.qr[i][k];
            }
        }
    }

    /// (beta[p], rss) for one gene.
    pub fn solve(&self, y: &[f64]) -> (Vec<f64>, f64) {
        let mut qty = y.to_vec();
        self.qty(&mut qty);
        let rss: f64 = qty[self.p..].iter().map(|&e| e * e).sum();
        let mut beta = vec![0.0; self.p];
        for j in (0..self.p).rev() {
            beta[j] = qty[j];
            for k in (j + 1)..self.p {
                beta[j] -= self.qr[j][k] * beta[k];
            }
            beta[j] /= self.rdiag[j];
        }
        (beta, rss)
    }

    /// The trailing d = n-p rotated residuals (Q2'y) for one gene.
    pub fn rotated_residual(&self, y: &[f64], out: &mut [f64]) {
        let mut qty = y.to_vec();
        self.qty(&mut qty);
        out.copy_from_slice(&qty[self.p..]);
    }

    /// (X'X)^-1 = R^-1 R^-T, the p×p unscaled coefficient covariance.
    #[allow(clippy::needless_range_loop)]
    pub fn xtx_inv(&self) -> Vec<Vec<f64>> {
        let p = self.p;
        let r_at =
            |i: usize, j: usize| -> f64 { if i == j { self.rdiag[i] } else { self.qr[i][j] } };
        let mut rinv = vec![vec![0.0; p]; p];
        for i in 0..p {
            rinv[i][i] = 1.0 / r_at(i, i);
        }
        for j in 0..p {
            for i in (0..j).rev() {
                let mut s = 0.0;
                for k in (i + 1)..=j {
                    s += r_at(i, k) * rinv[k][j];
                }
                rinv[i][j] = -s / r_at(i, i);
            }
        }
        let mut cov = vec![vec![0.0; p]; p];
        for i in 0..p {
            for j in 0..p {
                let mut s = 0.0;
                for (ri, rj) in rinv[i].iter().zip(&rinv[j]) {
                    s += ri * rj;
                }
                cov[i][j] = s;
            }
        }
        cov
    }
}

pub struct Fit {
    pub coef_names: Vec<String>,
    /// [gene][coef]
    pub coefficients: Vec<Vec<f64>>,
    /// per-coef unscaled sd = sqrt(diag (X'X)^-1), shared across genes
    pub stdev_unscaled: Vec<f64>,
    pub sigma: Vec<f64>,
    pub df_residual: f64,
    pub amean: Vec<f64>,
    pub genes: Vec<String>,
}

pub fn lm_fit(
    y: &[Vec<f64>],
    genes: &[String],
    x: &[Vec<f64>],
    coef_names: &[String],
) -> Result<(Fit, Qr)> {
    let n = x.len();
    let p = x[0].len();
    if y.iter().any(|row| row.len() != n) {
        return Err(RsomicsError::InvalidInput(
            "expression samples do not match design rows".into(),
        ));
    }
    let df_residual = (n - p) as f64;
    if df_residual < 1.0 {
        return Err(RsomicsError::InvalidInput(
            "residual degrees of freedom < 1 (need more samples than coefficients)".into(),
        ));
    }
    let qr = Qr::new(x)?;
    let cov = qr.xtx_inv();
    let stdev_unscaled: Vec<f64> = (0..p).map(|j| cov[j][j].sqrt()).collect();

    let mut coefficients = Vec::with_capacity(y.len());
    let mut sigma = Vec::with_capacity(y.len());
    let mut amean = Vec::with_capacity(y.len());
    for row in y {
        let (beta, rss) = qr.solve(row);
        coefficients.push(beta);
        sigma.push((rss / df_residual).sqrt());
        amean.push(row.iter().sum::<f64>() / n as f64);
    }

    Ok((
        Fit {
            coef_names: coef_names.to_vec(),
            coefficients,
            stdev_unscaled,
            sigma,
            df_residual,
            amean,
            genes: genes.to_vec(),
        },
        qr,
    ))
}

/// contrasts.fit: transform a coefficient-space fit into contrast space.
#[allow(clippy::needless_range_loop)]
pub fn contrasts_fit(fit: &Fit, contrast: &crate::matrix::Contrast, xtx_inv: &[Vec<f64>]) -> Fit {
    let p = fit.coef_names.len();
    let q = contrast.names.len();
    let cmat = &contrast.c;

    let mut new_coef = Vec::with_capacity(fit.coefficients.len());
    for beta in &fit.coefficients {
        let mut nb = vec![0.0; q];
        for (col, item) in nb.iter_mut().enumerate() {
            let mut s = 0.0;
            for (i, &b) in beta.iter().enumerate() {
                s += cmat[i][col] * b;
            }
            *item = s;
        }
        new_coef.push(nb);
    }

    let mut stdev = vec![0.0; q];
    for (col, sd) in stdev.iter_mut().enumerate() {
        let mut acc = 0.0;
        for i in 0..p {
            for j in 0..p {
                acc += cmat[i][col] * xtx_inv[i][j] * cmat[j][col];
            }
        }
        *sd = acc.sqrt();
    }

    Fit {
        coef_names: contrast.names.clone(),
        coefficients: new_coef,
        stdev_unscaled: stdev,
        sigma: fit.sigma.clone(),
        df_residual: fit.df_residual,
        amean: fit.amean.clone(),
        genes: fit.genes.clone(),
    }
}
