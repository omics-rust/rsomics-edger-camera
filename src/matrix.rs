use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

fn open(path: &Path) -> Result<BufReader<File>> {
    let f = File::open(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
    Ok(BufReader::new(f))
}

fn parse_f64(s: &str) -> Result<f64> {
    let t = s.trim();
    t.parse::<f64>()
        .map_err(|_| RsomicsError::InvalidInput(format!("non-numeric value '{t}'")))
}

pub struct Expr {
    pub samples: Vec<String>,
    pub genes: Vec<String>,
    /// row-major [gene][sample]
    pub y: Vec<Vec<f64>>,
}

pub fn read_expr(path: &Path) -> Result<Expr> {
    let mut lines = open(path)?.lines();
    let header = lines
        .next()
        .ok_or_else(|| RsomicsError::InvalidInput("empty expression matrix".into()))?
        .map_err(RsomicsError::Io)?;
    let samples: Vec<String> = header.split('\t').skip(1).map(str::to_string).collect();
    if samples.is_empty() {
        return Err(RsomicsError::InvalidInput(
            "expression matrix needs at least one sample column".into(),
        ));
    }
    let mut genes = Vec::new();
    let mut y = Vec::new();
    for line in lines {
        let line = line.map_err(RsomicsError::Io)?;
        if line.is_empty() {
            continue;
        }
        let mut f = line.split('\t');
        let gene = f
            .next()
            .ok_or_else(|| RsomicsError::InvalidInput("missing gene id".into()))?;
        let row: Vec<f64> = f.map(parse_f64).collect::<Result<_>>()?;
        if row.len() != samples.len() {
            return Err(RsomicsError::InvalidInput(format!(
                "gene '{gene}' has {} values, header declares {} samples",
                row.len(),
                samples.len()
            )));
        }
        genes.push(gene.to_string());
        y.push(row);
    }
    if genes.is_empty() {
        return Err(RsomicsError::InvalidInput("no genes in matrix".into()));
    }
    Ok(Expr { samples, genes, y })
}

pub struct Design {
    pub coef_names: Vec<String>,
    /// row-major [sample][coef]
    pub x: Vec<Vec<f64>>,
    pub row_ids: Vec<String>,
}

/// Design TSV: first column = sample id, remaining columns are coefficient
/// names with the numeric model-matrix entries, one row per sample.
pub fn read_design(path: &Path) -> Result<Design> {
    let mut lines = open(path)?.lines();
    let header = lines
        .next()
        .ok_or_else(|| RsomicsError::InvalidInput("empty design matrix".into()))?
        .map_err(RsomicsError::Io)?;
    let coef_names: Vec<String> = header.split('\t').skip(1).map(str::to_string).collect();
    if coef_names.is_empty() {
        return Err(RsomicsError::InvalidInput(
            "design matrix needs at least one coefficient column".into(),
        ));
    }
    let mut x = Vec::new();
    let mut row_ids = Vec::new();
    for line in lines {
        let line = line.map_err(RsomicsError::Io)?;
        if line.is_empty() {
            continue;
        }
        let mut f = line.split('\t');
        let id = f
            .next()
            .ok_or_else(|| RsomicsError::InvalidInput("missing design row id".into()))?;
        let row: Vec<f64> = f.map(parse_f64).collect::<Result<_>>()?;
        if row.len() != coef_names.len() {
            return Err(RsomicsError::InvalidInput(format!(
                "design row '{id}' has {} values, header declares {} coefficients",
                row.len(),
                coef_names.len()
            )));
        }
        row_ids.push(id.to_string());
        x.push(row);
    }
    Ok(Design {
        coef_names,
        x,
        row_ids,
    })
}

pub struct Contrast {
    pub names: Vec<String>,
    /// row-major [coef][contrast]
    pub c: Vec<Vec<f64>>,
}

/// Contrast TSV: first column = original coefficient name, header = contrast
/// names; entry [coef][contrast] is the contrast weight.
pub fn read_contrast(path: &Path, coef_names: &[String]) -> Result<Contrast> {
    let mut lines = open(path)?.lines();
    let header = lines
        .next()
        .ok_or_else(|| RsomicsError::InvalidInput("empty contrast matrix".into()))?
        .map_err(RsomicsError::Io)?;
    let names: Vec<String> = header.split('\t').skip(1).map(str::to_string).collect();
    if names.is_empty() {
        return Err(RsomicsError::InvalidInput(
            "contrast matrix needs at least one contrast column".into(),
        ));
    }
    let mut by_coef = std::collections::HashMap::new();
    for line in lines {
        let line = line.map_err(RsomicsError::Io)?;
        if line.is_empty() {
            continue;
        }
        let mut f = line.split('\t');
        let coef = f
            .next()
            .ok_or_else(|| RsomicsError::InvalidInput("missing contrast row id".into()))?;
        let row: Vec<f64> = f.map(parse_f64).collect::<Result<_>>()?;
        if row.len() != names.len() {
            return Err(RsomicsError::InvalidInput(format!(
                "contrast row '{coef}' width mismatch"
            )));
        }
        by_coef.insert(coef.to_string(), row);
    }
    let mut c = Vec::with_capacity(coef_names.len());
    for name in coef_names {
        let row = by_coef
            .get(name)
            .cloned()
            .unwrap_or_else(|| vec![0.0; names.len()]);
        c.push(row);
    }
    Ok(Contrast { names, c })
}

pub struct GeneSet {
    pub name: String,
    /// gene ids declared in the set, in file order, deduplicated.
    pub members: Vec<String>,
}

/// GMT: one set per line, tab-separated — set name, a description (ignored),
/// then the gene ids. Set names must be unique.
pub fn read_gmt(path: &Path) -> Result<Vec<GeneSet>> {
    let mut sets = Vec::new();
    let mut seen = std::collections::HashSet::new();
    for line in open(path)?.lines() {
        let line = line.map_err(RsomicsError::Io)?;
        if line.trim().is_empty() {
            continue;
        }
        let mut f = line.split('\t');
        let name = f
            .next()
            .ok_or_else(|| RsomicsError::InvalidInput("gene set with no name".into()))?
            .to_string();
        let _description = f.next();
        let mut members = Vec::new();
        let mut local = std::collections::HashSet::new();
        for g in f {
            let g = g.trim();
            if !g.is_empty() && local.insert(g.to_string()) {
                members.push(g.to_string());
            }
        }
        if members.is_empty() {
            return Err(RsomicsError::InvalidInput(format!(
                "gene set '{name}' has no gene ids"
            )));
        }
        if !seen.insert(name.clone()) {
            return Err(RsomicsError::InvalidInput(format!(
                "duplicate gene set name '{name}'"
            )));
        }
        sets.push(GeneSet { name, members });
    }
    if sets.is_empty() {
        return Err(RsomicsError::InvalidInput("no gene sets in GMT".into()));
    }
    Ok(sets)
}
