//! Empirical-Bayes variance moderation (Smyth 2004, Stat Appl Genet Mol Biol
//! 3:1, doi:10.2202/1544-6115.1027): moment-fit the scaled-inverse-chisquare
//! prior, shrink each gene's residual variance toward it, and form the
//! moderated t on df.total = df.residual + df.prior.

use crate::fit::Fit;
use crate::special::{digamma, trigamma, trigamma_inverse};

pub struct Moderated {
    /// moderated t per [gene][coef]
    pub t: Vec<Vec<f64>>,
    pub df_total: f64,
}

fn median(x: &[f64]) -> f64 {
    let mut s = x.to_vec();
    s.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = s.len();
    if n % 2 == 1 {
        s[n / 2]
    } else {
        0.5 * (s[n / 2 - 1] + s[n / 2])
    }
}

/// limma fitFDist (scalar df1): moment-fit the scaled-inverse-chisquare prior.
/// Zero and near-zero variances are offset away from zero (floor at 1e-5·median)
/// before the moment estimation, matching limma's "Zero sample variances
/// detected, have been offset away from zero".
fn fit_f_dist(x: &[f64], df1: f64) -> (f64, f64) {
    let n = x.len();
    if n == 0 {
        return (f64::NAN, f64::NAN);
    }
    if n == 1 {
        return (x[0], 0.0);
    }
    let half = df1 / 2.0;
    let df1_ok = df1.is_finite() && df1 > 1e-15;
    let ok: Vec<f64> = if df1_ok {
        x.iter()
            .copied()
            .filter(|v| v.is_finite() && *v > -1e-15)
            .collect()
    } else {
        vec![]
    };
    let nok = ok.len();
    if nok == 0 {
        return (f64::NAN, f64::NAN);
    }
    if nok == 1 {
        return (ok[0], 0.0);
    }
    let mut xs: Vec<f64> = ok.iter().map(|v| v.max(0.0)).collect();
    let m = median(&xs);
    let m = if m == 0.0 { 1.0 } else { m };
    let floor = 1e-5 * m;
    for v in &mut xs {
        *v = v.max(floor);
    }
    let nf = nok as f64;
    let e: Vec<f64> = xs
        .iter()
        .map(|v| v.ln() - digamma(half) + half.ln())
        .collect();
    let emean: f64 = e.iter().sum::<f64>() / nf;
    let evar: f64 =
        e.iter().map(|&v| (v - emean).powi(2)).sum::<f64>() / (nf - 1.0) - trigamma(half);
    if evar > 0.0 {
        let df2 = 2.0 * trigamma_inverse(evar);
        let s20 = (emean + digamma(df2 / 2.0) - (df2 / 2.0).ln()).exp();
        (s20, df2)
    } else {
        (xs.iter().sum::<f64>() / nf, f64::INFINITY)
    }
}

fn squeeze_var(sigma2: &[f64], df: f64) -> (Vec<f64>, f64) {
    let (s20, df0) = fit_f_dist(sigma2, df);
    let s2_post: Vec<f64> = if df0.is_infinite() {
        sigma2.iter().map(|_| s20).collect()
    } else {
        sigma2
            .iter()
            .map(|&s2| (df0 * s20 + df * s2) / (df0 + df))
            .collect()
    };
    (s2_post, df0)
}

pub fn ebayes(fit: &Fit) -> Moderated {
    let sigma2: Vec<f64> = fit.sigma.iter().map(|s| s * s).collect();
    let (s2_post, df_prior) = squeeze_var(&sigma2, fit.df_residual);

    let ng = fit.coefficients.len();
    let q = fit.coef_names.len();
    let df_pooled = ng as f64 * fit.df_residual;
    let df_total = (fit.df_residual + df_prior).min(df_pooled);

    let mut t = vec![vec![0.0; q]; ng];
    for (gi, row) in t.iter_mut().enumerate() {
        let post_sd = s2_post[gi].sqrt();
        for (cj, tv) in row.iter_mut().enumerate() {
            *tv = fit.coefficients[gi][cj] / (fit.stdev_unscaled[cj] * post_sd);
        }
    }

    Moderated { t, df_total }
}
