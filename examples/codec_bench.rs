//! Compare FRACSYN2 with framed bitcode, keeping production I/O unchanged.
//! cargo run --release --example codec_bench -- --out /tmp/fracsync-codecs
//! Optional: --input existing.sig; --repeats 5. Worker processes isolate RSS.
#[path = "support/blocked.rs"]
mod blocked;

use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::Instant;

use anyhow::{Context, Result, ensure};
use clap::Parser;
use fracsync::select::{SelectConfig, SelectorKind};
use fracsync::sketch::{Signature, SignatureFile};

#[derive(Parser)]
struct Args
{
    #[arg(long, default_value = "/tmp/fracsync-codecs")]
    out: PathBuf,
    /// Benchmark an existing FRACSYN2 database instead of generated fixtures.
    #[arg(long)]
    input: Option<PathBuf>,
    #[arg(long, default_value_t = 5)]
    repeats: usize,
    #[arg(long, hide = true)]
    worker: Option<String>,
    #[arg(long, hide = true)]
    codec: Option<String>,
}

fn next(state: &mut u64) -> u64
{
    *state ^= *state << 13;
    *state ^= *state >> 7;
    *state ^= *state << 17;
    *state
}

fn signature(name: &str, hashes: Vec<u64>, scale: u64) -> Signature
{
    Signature {
        name: name.into(),
        k: 21,
        scale,
        kind: SelectorKind::OpenSyncmer { s: 11, offset: 0 },
        hashes,
    }
}

fn make_fixture(label: &str) -> SignatureFile
{
    let mut rng = 42;
    let signatures = match label
    {
        "dna-scale1" | "dna-scale10" =>
        {
            let seq: Vec<_> = (0..5_000_000)
                .map(|_| b"ACGT"[(next(&mut rng) % 4) as usize])
                .collect();
            let scale = if label == "dna-scale1" { 1 } else { 10 };
            let cfg = SelectConfig {
                k: 21,
                scale,
                kind: SelectorKind::OpenSyncmer { s: 11, offset: 0 },
            };
            let mut hashes = Vec::new();
            cfg.select(&seq, &mut |h| hashes.push(h)).unwrap();
            hashes.sort_unstable();
            hashes.dedup();
            vec![signature(label, hashes, scale)]
        }
        "many-small" => (0..1024)
            .map(|i| {
                let mut hashes: Vec<_> = (0..(1 + i % 127)).map(|_| next(&mut rng)).collect();
                hashes.sort_unstable();
                hashes.dedup();
                signature(&format!("ref-{i}"), hashes, 1)
            })
            .collect(),
        "large-uniform" =>
        {
            let mut hashes: Vec<_> = (0..8_000_000).map(|_| next(&mut rng)).collect();
            hashes.sort_unstable();
            hashes.dedup();
            vec![signature(label, hashes, 1)]
        }
        _ => unreachable!(),
    };
    SignatureFile { signatures }
}

fn read(codec: &str, path: &Path) -> Result<SignatureFile>
{
    if codec == "bincode"
    {
        SignatureFile::read(path)
    }
    else
    {
        blocked::read(path)
    }
}
fn write(codec: &str, sf: &SignatureFile, path: &Path) -> Result<()>
{
    match codec
    {
        "bincode" => sf.write(path),
        "bitcode-raw" => blocked::write(sf, path, false),
        "bitcode-delta" => blocked::write(sf, path, true),
        _ => anyhow::bail!("unknown codec"),
    }
}
fn equal(a: &SignatureFile, b: &SignatureFile) -> bool
{
    a.signatures.len() == b.signatures.len()
        && a.signatures.iter().zip(&b.signatures).all(|(a, b)| {
            a.name == b.name
                && a.k == b.k
                && a.scale == b.scale
                && a.kind == b.kind
                && a.hashes == b.hashes
        })
}
fn peak_kib() -> u64
{
    std::fs::read_to_string("/proc/self/status")
        .ok()
        .and_then(|s| {
            s.lines()
                .find(|s| s.starts_with("VmHWM:"))
                .and_then(|line| line.split_whitespace().nth(1)?.parse().ok())
        })
        .unwrap_or(0)
}

fn worker(args: &Args, operation: &str) -> Result<()>
{
    let input = args.input.as_ref().context("worker needs input")?;
    let codec = args.codec.as_deref().context("worker needs codec")?;
    let (elapsed, peak, initial);
    if operation == "write"
    {
        let sf = SignatureFile::read(input)?;
        initial = peak_kib();
        let start = Instant::now();
        write(codec, &sf, &args.out)?;
        elapsed = start.elapsed();
        peak = peak_kib();
        // Full equality outside timing and after the RSS sample.
        ensure!(equal(&sf, &read(codec, &args.out)?), "write roundtrip mismatch");
    }
    else
    {
        initial = peak_kib();
        let start = Instant::now();
        let sf = read(codec, &args.out)?;
        elapsed = start.elapsed();
        peak = peak_kib();
        ensure!(equal(&sf, &SignatureFile::read(input)?), "read mismatch");
    }
    println!("{} {} {} {}", elapsed.as_nanos(), std::fs::metadata(&args.out)?.len(), peak, initial);
    Ok(())
}

#[derive(Default)]
struct Measurements
{
    nanos: Vec<u128>,
    size: u64,
    peak: u64,
    extra: u64,
}
impl Measurements
{
    fn add(&mut self, text: &str) -> Result<()>
    {
        let fields: Vec<_> = text.split_whitespace().collect();
        ensure!(fields.len() == 4, "invalid worker result: {text}");
        self.nanos.push(fields[0].parse()?);
        self.size = fields[1].parse()?;
        let peak: u64 = fields[2].parse()?;
        let initial: u64 = fields[3].parse()?;
        self.peak = self.peak.max(peak);
        self.extra = self.extra.max(peak.saturating_sub(initial));
        Ok(())
    }
    fn median_ms(&mut self) -> f64
    {
        self.nanos.sort_unstable();
        self.nanos[self.nanos.len() / 2] as f64 / 1e6
    }
}

fn main() -> Result<()>
{
    let args = Args::parse();
    if let Some(operation) = args.worker.as_deref()
    {
        return worker(&args, operation);
    }
    ensure!(args.repeats > 0, "repeats must be positive");
    std::fs::create_dir_all(&args.out)?;
    let mut inputs = Vec::new();
    if let Some(path) = args.input
    {
        inputs.push(("user-input".to_string(), path));
    }
    else
    {
        for label in ["dna-scale1", "dna-scale10", "many-small", "large-uniform"]
        {
            let sf = make_fixture(label);
            let path = args.out.join(format!("{label}.sig"));
            eprintln!(
                "{label}: {} references, {} hashes",
                sf.signatures.len(),
                sf.signatures.iter().map(|s| s.hashes.len()).sum::<usize>()
            );
            sf.write(&path)?;
            inputs.push((label.into(), path));
        }
    }
    println!(
        "dataset\tcodec\tbytes\twrite_ms\tread_ms\twrite_peak_kib\twrite_extra_peak_kib\tread_peak_kib"
    );
    let codecs = ["bincode", "bitcode-raw", "bitcode-delta"];
    let exe = std::env::current_exe()?;
    for (label, input) in inputs
    {
        let mut results: Vec<_> = (0..3)
            .map(|_| (Measurements::default(), Measurements::default()))
            .collect();
        for repeat in 0..args.repeats
        {
            for step in 0..3
            {
                let idx = (repeat + step) % 3;
                let codec = codecs[idx];
                let output = args.out.join(format!("{label}.{codec}.bench"));
                for operation in ["write", "read"]
                {
                    let run = Command::new(&exe)
                        .arg("--worker")
                        .arg(operation)
                        .arg("--codec")
                        .arg(codec)
                        .arg("--input")
                        .arg(&input)
                        .arg("--out")
                        .arg(&output)
                        .output()?;
                    ensure!(
                        run.status.success(),
                        "worker failed: {}",
                        String::from_utf8_lossy(&run.stderr)
                    );
                    let m = if operation == "write"
                    {
                        &mut results[idx].0
                    }
                    else
                    {
                        &mut results[idx].1
                    };
                    m.add(std::str::from_utf8(&run.stdout)?)?;
                }
            }
        }
        for (codec, (mut write, mut read)) in codecs.iter().zip(results)
        {
            let write_ms = write.median_ms();
            let read_ms = read.median_ms();
            println!(
                "{label}\t{codec}\t{}\t{write_ms:.3}\t{read_ms:.3}\t{}\t{}\t{}",
                write.size, write.peak, write.extra, read.peak
            );
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests
{
    use super::*;
    #[test]
    fn block_boundaries_and_both_modes_roundtrip()
    {
        let path =
            std::env::temp_dir().join(format!("fracsync_blocked_{}.bench", std::process::id()));
        for count in [0, 1, 3, 63, 64, 65, 4095, 4096, 4097, 8193]
        {
            let sf = SignatureFile {
                signatures: vec![signature(
                    "référence",
                    (0..count).map(|i| i * 100_000_000).collect(),
                    1,
                )],
            };
            for delta in [false, true]
            {
                blocked::write(&sf, &path, delta).unwrap();
                assert!(equal(&sf, &blocked::read(&path).unwrap()));
            }
        }
        for hashes in [vec![0, u64::MAX], vec![u64::MAX]]
        {
            let sf = SignatureFile {
                signatures: vec![signature("edge", hashes, 1)],
            };
            blocked::write(&sf, &path, true).unwrap();
            assert!(equal(&sf, &blocked::read(&path).unwrap()));
        }
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn fixed_codec_fixtures_decode_and_reencode_identically()
    {
        let path = std::env::temp_dir()
            .join(format!("fracsync_codec_golden_{}.bench", std::process::id()));
        let sf = SignatureFile {
            signatures: vec![signature("golden", vec![0, 1, 42, u64::MAX], 1)],
        };
        for (delta, bytes) in [
            (false, include_bytes!("../tests/fixtures/bitcode-raw-v1.bin").as_slice()),
            (true, include_bytes!("../tests/fixtures/bitcode-delta-v1.bin").as_slice()),
        ]
        {
            std::fs::write(&path, bytes).unwrap();
            assert!(equal(&sf, &blocked::read(&path).unwrap()));
            blocked::write(&sf, &path, delta).unwrap();
            assert_eq!(std::fs::read(&path).unwrap(), bytes);
        }
        std::fs::remove_file(path).unwrap();
    }
}
