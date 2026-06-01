//! Differential compat against limma camera.
//!
//! - `golden_*` always runs: ours vs a committed R-captured camera table.
//! - `live_r_*` runs only when an Rscript with limma is found; it regenerates
//!   the oracle and diffs against ours (loud-skip otherwise).

use std::collections::BTreeMap;
use std::path::PathBuf;
use std::process::Command;

const EPS: f64 = 1e-6;

fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-edger-camera"))
}

fn golden(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

fn manifest(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(rel)
}

struct Table {
    header: Vec<String>,
    /// set name -> (numeric fields, direction string)
    rows: BTreeMap<String, (Vec<f64>, String)>,
}

fn parse(text: &str) -> Table {
    let mut lines = text.lines();
    let header: Vec<String> = lines
        .next()
        .unwrap()
        .split('\t')
        .map(str::to_string)
        .collect();
    let dir_col = header.iter().position(|h| h == "Direction").unwrap();
    let mut rows = BTreeMap::new();
    for line in lines {
        if line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        let name = f[0].to_string();
        let direction = f[dir_col].to_string();
        let nums: Vec<f64> = f
            .iter()
            .enumerate()
            .skip(1)
            .filter(|(i, _)| *i != dir_col)
            .map(|(_, s)| s.trim().parse().unwrap())
            .collect();
        rows.insert(name, (nums, direction));
    }
    Table { header, rows }
}

fn assert_close(a: &Table, b: &Table, label: &str) {
    assert_eq!(a.header, b.header, "{label}: header mismatch");
    assert_eq!(a.rows.len(), b.rows.len(), "{label}: row count mismatch");
    let mut max_rel = 0.0f64;
    for (name, (x, dx)) in &a.rows {
        let (y, dy) = b
            .rows
            .get(name)
            .unwrap_or_else(|| panic!("{label}: missing set {name}"));
        assert_eq!(dx, dy, "{label}: {name} direction {dx} vs {dy}");
        assert_eq!(x.len(), y.len(), "{label}: {name} width mismatch");
        for (vx, vy) in x.iter().zip(y) {
            let rel = (vx - vy).abs() / vy.abs().max(1e-9);
            max_rel = max_rel.max(rel);
            assert!(rel < EPS, "{label}: {name} ours={vx} ref={vy} rel={rel:e}");
        }
    }
    eprintln!("{label}: max relative deviation = {max_rel:e}");
}

fn run_ours(coef: usize, estimate: bool) -> String {
    let mut cmd = Command::new(ours());
    cmd.arg(golden("expr.tsv"))
        .args(["--design", golden("design.tsv").to_str().unwrap()])
        .args(["--gene-sets", golden("sets.gmt").to_str().unwrap()])
        .args(["--coef", &coef.to_string()]);
    if estimate {
        cmd.arg("--estimate-cor");
    }
    let out = cmd.output().unwrap();
    assert!(
        out.status.success(),
        "ours failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).unwrap()
}

#[test]
fn golden_fixed() {
    let ours_out = run_ours(2, false);
    let expected = std::fs::read_to_string(golden("camera.fixed.expected.tsv")).unwrap();
    assert_close(
        &parse(&ours_out),
        &parse(&expected),
        "camera fixed (golden)",
    );
}

#[test]
fn golden_estimate() {
    let ours_out = run_ours(2, true);
    let expected = std::fs::read_to_string(golden("camera.estimate.expected.tsv")).unwrap();
    assert_close(
        &parse(&ours_out),
        &parse(&expected),
        "camera estimate (golden)",
    );
}

fn rscript() -> Option<String> {
    let conda = format!(
        "{}/miniconda3/envs/r-bioc/bin/Rscript",
        std::env::var("HOME").unwrap_or_default()
    );
    for cand in [conda.as_str(), "Rscript"] {
        let ok = Command::new(cand)
            .args([
                "-e",
                "if(!requireNamespace('limma',quietly=TRUE)) quit(status=1)",
            ])
            .output()
            .map(|o| o.status.success())
            .unwrap_or(false);
        if ok {
            return Some(cand.to_string());
        }
    }
    None
}

fn live(mode: &str, estimate: bool, label: &str) {
    let Some(rs) = rscript() else {
        eprintln!("SKIP {label}: no Rscript with limma found");
        return;
    };
    let scratch = std::env::temp_dir();
    let r_out = scratch.join(format!("camera_r_{}_{}.tsv", mode, std::process::id()));
    let oracle = Command::new(&rs)
        .arg(manifest("tests/camera_oracle.R"))
        .arg(golden("expr.tsv"))
        .arg(golden("design.tsv"))
        .arg(golden("sets.gmt"))
        .arg("2")
        .arg(&r_out)
        .arg(mode)
        .output()
        .unwrap();
    assert!(
        oracle.status.success(),
        "oracle failed: {}",
        String::from_utf8_lossy(&oracle.stderr)
    );
    let ours_out = run_ours(2, estimate);
    let r_ref = std::fs::read_to_string(&r_out).unwrap();
    let _ = std::fs::remove_file(&r_out);
    assert_close(&parse(&ours_out), &parse(&r_ref), label);
}

#[test]
fn live_r_fixed() {
    live("fixed", false, "camera fixed (live R)");
}

#[test]
fn live_r_estimate() {
    live("estimate", true, "camera estimate (live R)");
}
