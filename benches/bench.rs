use std::hint::black_box;
use std::path::PathBuf;

use criterion::{Criterion, criterion_group, criterion_main};
use rsomics_edger_camera::{Options, run};

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

fn bench_camera(c: &mut Criterion) {
    let expr = fixture("expr.tsv");
    let design = fixture("design.tsv");
    let sets = fixture("sets.gmt");
    if !expr.exists() {
        return;
    }
    c.bench_function("camera", |b| {
        b.iter(|| {
            let opts = Options {
                expr: &expr,
                design: &design,
                gene_sets: &sets,
                contrast: None,
                coef: Some(2),
                inter_gene_cor: Some(0.01),
                allow_neg_cor: false,
                sort: true,
            };
            black_box(run(&opts).unwrap());
        })
    });
}

criterion_group!(benches, bench_camera);
criterion_main!(benches);
