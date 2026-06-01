use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, Tool, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

use rsomics_edger_camera::{Options, run, write_results};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-edger-camera", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    /// log-expression matrix TSV: header = sample ids, col 1 = gene ids.
    pub expr: PathBuf,
    /// Design matrix TSV: header = coefficient names, col 1 = sample ids.
    #[arg(long)]
    design: PathBuf,
    /// Gene sets in GMT: set name, description, then gene ids (tab-separated).
    #[arg(long)]
    gene_sets: PathBuf,
    /// Contrast matrix TSV (contrasts.fit): col 1 = coefficient names.
    #[arg(long)]
    contrast: Option<PathBuf>,
    /// 1-based coefficient (or contrast) to rank on; default = last.
    #[arg(long)]
    coef: Option<usize>,
    /// Preset inter-gene correlation shared by all sets.
    #[arg(long, default_value_t = 0.01)]
    inter_gene_cor: f64,
    /// Estimate the inter-gene correlation per set instead of using a preset.
    #[arg(long)]
    estimate_cor: bool,
    /// With --estimate-cor, allow a negative correlation (VIF below 1).
    #[arg(long)]
    allow_neg_cor: bool,
    /// Keep input set order instead of sorting by PValue.
    #[arg(long)]
    no_sort: bool,
    /// Results TSV destination; "-" is stdout.
    #[arg(short = 'o', long, default_value = "-")]
    output: String,
    #[command(flatten)]
    pub common: CommonFlags,
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }
    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        let opts = Options {
            expr: &self.expr,
            design: &self.design,
            gene_sets: &self.gene_sets,
            contrast: self.contrast.as_deref(),
            coef: self.coef,
            inter_gene_cor: if self.estimate_cor {
                None
            } else {
                Some(self.inter_gene_cor)
            },
            allow_neg_cor: self.allow_neg_cor,
            sort: !self.no_sort,
        };
        let res = run(&opts)?;

        let mut out: Box<dyn std::io::Write> = if self.output == "-" {
            Box::new(std::io::stdout().lock())
        } else {
            Box::new(std::fs::File::create(&self.output).map_err(RsomicsError::Io)?)
        };
        write_results(&res, &mut out)?;

        if !self.common.quiet {
            eprintln!(
                "{} sets tested, ranked on coef '{}'",
                res.rows.len(),
                res.coef_name
            );
        }
        Ok(())
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
    tagline: "Competitive gene-set test accounting for inter-gene correlation (limma camera).",
    origin: Some(Origin {
        upstream: "limma camera",
        upstream_license: "GPL (>=2)",
        our_license: "MIT OR Apache-2.0",
        paper_doi: Some("10.1093/nar/gks461"),
    }),
    usage_lines: &[
        "<expr.tsv> --design <design.tsv> --gene-sets <sets.gmt> [--coef N] [-o out.tsv]",
    ],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: None,
                long: "design",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("PathBuf"),
                required: true,
                default: None,
                description: "Design matrix TSV (header = coefficient names, col 1 = sample ids).",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "gene-sets",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("PathBuf"),
                required: true,
                default: None,
                description: "Gene sets in GMT format.",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "contrast",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("PathBuf"),
                required: false,
                default: None,
                description: "Contrast matrix TSV; applies contrasts.fit before ranking.",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "coef",
                aliases: &[],
                value: Some("<N>"),
                type_hint: Some("usize"),
                required: false,
                default: None,
                description: "1-based coefficient/contrast to rank on.",
                why_default: Some("Last coefficient — the typical treatment effect."),
            },
            FlagSpec {
                short: None,
                long: "inter-gene-cor",
                aliases: &[],
                value: Some("<rho>"),
                type_hint: Some("f64"),
                required: false,
                default: Some("0.01"),
                description: "Preset inter-gene correlation shared by all sets.",
                why_default: Some("limma's default since 3.29.6."),
            },
            FlagSpec {
                short: None,
                long: "estimate-cor",
                aliases: &[],
                value: None,
                type_hint: None,
                required: false,
                default: None,
                description: "Estimate the correlation per set instead of using the preset.",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "allow-neg-cor",
                aliases: &[],
                value: None,
                type_hint: None,
                required: false,
                default: None,
                description: "With --estimate-cor, allow a negative correlation (VIF below 1).",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "no-sort",
                aliases: &[],
                value: None,
                type_hint: None,
                required: false,
                default: None,
                description: "Keep input set order instead of sorting by PValue.",
                why_default: None,
            },
            FlagSpec {
                short: Some('o'),
                long: "output",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("String"),
                required: false,
                default: Some("-"),
                description: "Results TSV destination; \"-\" is stdout.",
                why_default: None,
            },
        ],
    }],
    examples: &[
        Example {
            description: "Default competitive test (preset correlation 0.01)",
            command: "rsomics-edger-camera E.tsv --design design.tsv --gene-sets sets.gmt --coef 2 -o camera.tsv",
        },
        Example {
            description: "Estimate the inter-gene correlation per set",
            command: "rsomics-edger-camera E.tsv --design design.tsv --gene-sets sets.gmt --estimate-cor > camera.tsv",
        },
    ],
    json_result_schema_doc: None,
};

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn cli_debug_assert() {
        Cli::command().debug_assert();
    }
}
