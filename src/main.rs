//! `fracsync` — sketch references and estimate containment of a query sample.
//!
//! ```sh
//! fracsync sketch refs/*.fa -o refs.sig          # k=21, closed syncmer s=11, scale=1
//! fracsync contain refs.sig reads.fastq.gz       # which references are present, by containment
//! fracsync info refs.sig                          # parameters and per-reference sizes
//! ```

use std::io::{self, Write};
use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::{Parser, ValueEnum};

use fracsync::contain;
use fracsync::select::{SelectConfig, SelectorKind};
use fracsync::sketch::{self, SignatureFile};

#[derive(Parser)]
#[command(
    name = "fracsync",
    version,
    about = "K-mer selection (strand-symmetric syncmers + FracMinHash) and containment estimation"
)]
enum Cli
{
    /// Sketch each reference FASTA/FASTQ into a signature file (one reference per file)
    Sketch
    {
        /// Reference files (FASTA/FASTQ, optionally gzipped)
        #[arg(required = true)]
        paths: Vec<PathBuf>,
        /// Output signature file
        #[arg(long, short)]
        out: PathBuf,
        #[command(flatten)]
        select: SelectArgs,
    },
    /// Estimate containment of each reference in a query sample of reads
    Contain
    {
        /// Signature file from `sketch`
        db: PathBuf,
        /// Query reads (FASTA/FASTQ, optionally gzipped); all records form one sample
        #[arg(required = true)]
        queries: Vec<PathBuf>,
        /// Only report references with at least this containment
        #[arg(long, default_value_t = 0.0)]
        min_containment: f64,
    },
    /// Print the parameters and per-reference sizes of a signature file
    Info
    {
        /// Signature file from `sketch`
        db: PathBuf,
    },
}

#[derive(clap::Args)]
struct SelectArgs
{
    /// k-mer length (1..=32)
    #[arg(long, short, default_value_t = 21)]
    k: usize,
    /// FracMinHash scale D: keep ~1/D of hashes (1 keeps all)
    #[arg(long, default_value_t = 1)]
    scale: u64,
    /// Selection scheme
    #[arg(long, value_enum, default_value_t = Selector::Syncmer)]
    selector: Selector,
    /// Syncmer s-mer length (must be < k)
    #[arg(long, default_value_t = 11)]
    s: u8,
    /// Syncmer offset (and its mirror): 0 is closed; a single centre offset is open
    #[arg(long, default_value_t = 0)]
    offset: u8,
    /// Minimizer window length (number of k-mers)
    #[arg(long, default_value_t = 11)]
    w: u8,
}

#[derive(Copy, Clone, ValueEnum)]
enum Selector
{
    /// Syncmer: error-robust, strand-symmetric (default)
    Syncmer,
    /// Minimizer: for comparison
    Minimizer,
    /// Every k-mer under the scale threshold (plain FracMinHash): for comparison
    All,
}

impl SelectArgs
{
    fn config(&self) -> SelectConfig
    {
        let kind = match self.selector
        {
            Selector::Syncmer => SelectorKind::OpenSyncmer {
                s: self.s,
                offset: self.offset,
            },
            Selector::Minimizer => SelectorKind::Minimizer { w: self.w },
            Selector::All => SelectorKind::All,
        };
        SelectConfig {
            k: self.k,
            scale: self.scale,
            kind,
        }
    }
}

fn selector_desc(cfg: &SelectConfig) -> String
{
    match cfg.kind
    {
        SelectorKind::OpenSyncmer { s, offset } =>
        {
            format!("syncmer s={s} offset={offset} mirror={}", cfg.k - s as usize - offset as usize)
        }
        SelectorKind::Minimizer { w } => format!("minimizer w={w}"),
        SelectorKind::All => "all k-mers (FracMinHash)".to_string(),
    }
}

fn run() -> Result<()>
{
    match Cli::parse()
    {
        Cli::Sketch { paths, out, select } =>
        {
            let cfg = select.config();
            cfg.validate().map_err(anyhow::Error::msg)?;
            let signatures = sketch::sketch_files(&paths, &cfg)?;
            let total: usize = signatures.iter().map(|s| s.hashes.len()).sum();
            let reference_count = signatures.len();
            SignatureFile { signatures }.write(&out)?;
            eprintln!(
                "sketched {} references ({} selected hashes) at k={} scale={} {} into {}",
                reference_count,
                total,
                cfg.k,
                cfg.scale,
                selector_desc(&cfg),
                out.display()
            );
        }
        Cli::Contain {
            db,
            queries,
            min_containment,
        } =>
        {
            let sigs = SignatureFile::read(&db)?;
            let cfg = sigs
                .config()
                .context("signature file is empty; nothing to compare against")?;

            let sample = sketch::select_files(&queries, &cfg)?;

            let mut scores: Vec<_> = sigs
                .signatures
                .iter()
                .map(|s| contain::score(&s.name, &s.hashes, &sample, cfg.k))
                .filter(|s| s.containment >= min_containment)
                .collect();
            // Most-contained first; ties broken by the more specific match
            // (larger fraction of the sample explained).
            scores.sort_by(|a, b| {
                b.containment
                    .total_cmp(&a.containment)
                    .then(b.query_containment.total_cmp(&a.query_containment))
            });

            eprintln!(
                "sample: {} distinct selected hashes (k={} scale={} {})",
                sample.len(),
                cfg.k,
                cfg.scale,
                selector_desc(&cfg)
            );
            let mut out = io::BufWriter::new(io::stdout().lock());
            writeln!(out, "#reference\tref_kmers\tshared\tcontainment\tani\tsample_frac")?;
            for s in &scores
            {
                writeln!(
                    out,
                    "{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}",
                    s.name, s.ref_size, s.shared, s.containment, s.ani, s.query_containment
                )?;
            }
            out.flush()?;
        }
        Cli::Info { db } =>
        {
            let sigs = SignatureFile::read(&db)?;
            match sigs.config()
            {
                Some(cfg) => eprintln!(
                    "{}: {} references, k={} scale={} {}",
                    db.display(),
                    sigs.signatures.len(),
                    cfg.k,
                    cfg.scale,
                    selector_desc(&cfg)
                ),
                None => eprintln!("{}: empty signature file", db.display()),
            }
            let mut out = io::BufWriter::new(io::stdout().lock());
            writeln!(out, "#reference\tselected_hashes")?;
            for s in &sigs.signatures
            {
                writeln!(out, "{}\t{}", s.name, s.hashes.len())?;
            }
            out.flush()?;
        }
    }
    Ok(())
}

fn main()
{
    if let Err(e) = run()
    {
        eprintln!("error: {e:#}");
        std::process::exit(1);
    }
}
