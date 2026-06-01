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

fn fit_f_dist(x: &[f64], df1: f64) -> (f64, f64) {
    let half = df1 / 2.0;
    let mut e = Vec::with_capacity(x.len());
    for &xi in x {
        if xi > 0.0 {
            e.push(xi.ln() - digamma(half) + half.ln());
        }
    }
    let m = e.len() as f64;
    let emean: f64 = e.iter().sum::<f64>() / m;
    let evar: f64 = e.iter().map(|&v| (v - emean).powi(2)).sum::<f64>() / (m - 1.0);
    let evar = evar - trigamma(half);
    if evar > 0.0 {
        let df2 = 2.0 * trigamma_inverse(evar);
        let s20 = (emean + digamma(df2 / 2.0) - (df2 / 2.0).ln()).exp();
        (s20, df2)
    } else {
        (emean.exp(), f64::INFINITY)
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
