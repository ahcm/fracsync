//! Reference signatures and record-at-a-time FASTA/FASTQ sketching.
//!
//! Selectors use bounded rolling state; the parser still buffers a record and
//! the distinct selected hash set grows with the input's selected content.

use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Take, Write};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use bincode::{Decode, Encode};
use fastx::FastX::{self, FastARecord, FastQRead, FastQRecord, FastXFormat, FastXRead};
use flate2::read::MultiGzDecoder;

use crate::select::{SelectConfig, SelectorKind};

/// Magic + format version, so a `.sig` file is self-describing and mismatched
/// formats fail loudly instead of decoding into garbage.
const SIG_MAGIC: &[u8; 8] = b"FRACSYN2";

/// One reference reduced to its selected hashes.
#[derive(Encode, Decode, Clone)]
pub struct Signature
{
    pub name: String,
    pub k: u8,
    pub scale: u64,
    pub kind: SelectorKind,
    /// Sorted, deduplicated selected canonical hashes.
    pub hashes: Vec<u64>,
}

/// A `.sig` file holds one or more signatures sharing the same parameters.
#[derive(Encode, Decode, Clone)]
pub struct SignatureFile
{
    pub signatures: Vec<Signature>,
}

impl SignatureFile
{
    /// The selection parameters shared by every signature (from the first one).
    /// A file with mixed parameters is rejected on read, so this is well-defined.
    pub fn config(&self) -> Option<SelectConfig>
    {
        self.signatures.first().map(|s| SelectConfig {
            k: s.k as usize,
            scale: s.scale,
            kind: s.kind,
        })
    }

    /// Checks the invariants needed by selection and merge-based scoring.
    pub fn validate(&self) -> Result<()>
    {
        if let Some(cfg) = self.config()
        {
            cfg.validate().map_err(anyhow::Error::msg)?;
            let threshold = u64::MAX / cfg.scale;
            for s in &self.signatures
            {
                anyhow::ensure!(
                    s.k as usize == cfg.k && s.scale == cfg.scale && s.kind == cfg.kind,
                    "signatures have mixed selection parameters"
                );
                anyhow::ensure!(
                    s.hashes.windows(2).all(|w| w[0] < w[1]),
                    "signature {:?}: hashes must be sorted and unique",
                    s.name
                );
                anyhow::ensure!(
                    s.hashes.last().is_none_or(|&h| h <= threshold),
                    "signature {:?}: hash exceeds scale threshold",
                    s.name
                );
            }
        }
        Ok(())
    }

    pub fn write(&self, path: &Path) -> Result<()>
    {
        // Validate before opening so invalid data cannot truncate an existing file.
        self.validate()?;
        let mut out = BufWriter::new(
            File::create(path).with_context(|| format!("creating {}", path.display()))?,
        );
        out.write_all(SIG_MAGIC)?;
        bincode::encode_into_std_write(self, &mut out, bincode::config::standard())
            .context("encoding signatures")?;
        out.flush()
            .with_context(|| format!("writing {}", path.display()))?;
        Ok(())
    }

    pub fn read(path: &Path) -> Result<Self>
    {
        let mut input = BufReader::new(
            File::open(path).with_context(|| format!("reading {}", path.display()))?,
        );
        let mut magic = [0; 8];
        input
            .read_exact(&mut magic)
            .context("reading signature header")?;
        if &magic == b"FRACSYN1"
        {
            anyhow::bail!(
                "{}: legacy FRACSYN1 hashes; re-sketch the source sequences with this version",
                path.display()
            );
        }
        anyhow::ensure!(
            &magic == SIG_MAGIC,
            "{}: unsupported fracsync signature format",
            path.display()
        );
        // Decode collections incrementally: derived Vec decoding trusts length
        // prefixes enough to allocate before discovering a truncated payload.
        let remaining = input.get_ref().metadata()?.len().saturating_sub(8);
        let mut body = input.take(remaining);
        let sf = read_body(&mut body).with_context(|| format!("decoding {}", path.display()))?;
        anyhow::ensure!(
            body.read(&mut [0u8; 1])? == 0,
            "{}: trailing data after signatures",
            path.display()
        );
        sf.validate()
            .with_context(|| format!("validating {}", path.display()))?;
        Ok(sf)
    }
}

// These helpers read the same standard-bincode field layout that Encode emits,
// but do not preallocate whole collections from untrusted length prefixes.
fn read_value<T: Decode<()>>(input: &mut impl Read) -> Result<T>
{
    Ok(bincode::decode_from_std_read(input, bincode::config::standard())?)
}

fn read_count(input: &mut Take<impl Read>) -> Result<usize>
{
    let len: usize = read_value(input)?;
    // Every element used by this format occupies at least one encoded byte.
    anyhow::ensure!(
        len as u64 <= input.limit(),
        "collection length exceeds remaining signature bytes"
    );
    Ok(len)
}

fn read_values<T: Decode<()>>(input: &mut Take<impl Read>) -> Result<Vec<T>>
{
    let len = read_count(input)?;
    let mut values = Vec::with_capacity(len.min(4096));
    for _ in 0..len
    {
        values.push(read_value(input)?);
    }
    Ok(values)
}

fn read_body(input: &mut Take<impl Read>) -> Result<SignatureFile>
{
    let count = read_count(input)?;
    let mut signatures = Vec::new();
    for _ in 0..count
    {
        signatures.push(Signature {
            name: String::from_utf8(read_values(input)?)?,
            k: read_value(input)?,
            scale: read_value(input)?,
            kind: read_value(input)?,
            hashes: read_values(input)?,
        });
    }
    Ok(SignatureFile { signatures })
}

/// Derives a reference name from a file path (stem, sans common extensions).
fn name_from_path(path: &Path) -> String
{
    let mut name = path
        .file_name()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| path.display().to_string());
    for ext in [".gz", ".fasta", ".fa", ".fna", ".fastq", ".fq"]
    {
        if let Some(stripped) = name.strip_suffix(ext)
        {
            name = stripped.to_string();
        }
    }
    name
}

/// Selects sorted, distinct hashes across all records in one file.
/// Syncmers are content-defined; minimizers depend on the record's context.
pub fn select_file(path: &Path, cfg: &SelectConfig) -> Result<Vec<u64>>
{
    select_paths(std::iter::once(path), cfg)
}

/// Selects a union across files, deduplicating as records arrive so repeated
/// query files do not accumulate duplicate vectors in memory.
pub fn select_files(paths: &[PathBuf], cfg: &SelectConfig) -> Result<Vec<u64>>
{
    select_paths(paths.iter().map(PathBuf::as_path), cfg)
}

fn select_paths<'a>(
    paths: impl IntoIterator<Item = &'a Path>,
    cfg: &SelectConfig,
) -> Result<Vec<u64>>
{
    cfg.validate().map_err(anyhow::Error::msg)?;
    let mut set = HashSet::new();
    for path in paths
    {
        let mut reader =
            sequence_reader(path).with_context(|| format!("opening {}", path.display()))?;
        // fastx 0.6.1's peek indexes the first byte; guard empty input before
        // invoking it, including an empty decompressed gzip stream.
        anyhow::ensure!(
            !reader
                .fill_buf()
                .with_context(|| format!("reading {}", path.display()))?
                .is_empty(),
            "{}: expected FASTA or FASTQ input",
            path.display()
        );
        let (format, _) = FastX::peek(&mut reader)
            .with_context(|| format!("detecting sequence format in {}", path.display()))?;
        match format
        {
            FastXFormat::FASTA =>
            {
                let mut record = FastARecord::default();
                while record
                    .read(&mut reader)
                    .with_context(|| format!("parsing FASTA record in {}", path.display()))?
                    > 0
                {
                    // Join wrapped lines within this record without another
                    // sequence allocation; read() reuses the record buffer.
                    record.raw_seq.retain(|&b| b != b'\n' && b != b'\r');
                    cfg.select(record.seq_raw(), &mut |h| {
                        set.insert(h);
                    })
                    .map_err(anyhow::Error::msg)?;
                }
            }
            FastXFormat::FASTQ =>
            {
                let mut record = FastQRecord::default();
                while record
                    .read(&mut reader)
                    .with_context(|| format!("parsing FASTQ record in {}", path.display()))?
                    > 0
                {
                    anyhow::ensure!(
                        record.seq_raw().len() == record.qual().len(),
                        "{}: FASTQ sequence and quality lengths differ for {:?}",
                        path.display(),
                        record.id()
                    );
                    cfg.select(record.seq_raw(), &mut |h| {
                        set.insert(h);
                    })
                    .map_err(anyhow::Error::msg)?;
                }
            }
            FastXFormat::EOF | FastXFormat::UNKNOWN =>
            {
                anyhow::bail!("{}: expected FASTA or FASTQ input", path.display());
            }
        }
    }
    let mut hashes: Vec<u64> = set.into_iter().collect();
    hashes.sort_unstable();
    Ok(hashes)
}

// FastX's path helper uses a 600 MiB buffer. Supply a smaller reader to its
// record API and detect gzip by magic, preserving support regardless of suffix.
fn sequence_reader(path: &Path) -> std::io::Result<Box<dyn BufRead>>
{
    const BUFFER_SIZE: usize = 64 * 1024;
    let mut input = BufReader::with_capacity(BUFFER_SIZE, File::open(path)?);
    if input.fill_buf()?.starts_with(&[0x1f, 0x8b])
    {
        Ok(Box::new(BufReader::with_capacity(BUFFER_SIZE, MultiGzDecoder::new(input))))
    }
    else
    {
        Ok(Box::new(input))
    }
}

/// Sketches each input file into one [`Signature`] (one reference per file).
pub fn sketch_files(paths: &[PathBuf], cfg: &SelectConfig) -> Result<Vec<Signature>>
{
    cfg.validate().map_err(anyhow::Error::msg)?;
    let mut sigs = Vec::with_capacity(paths.len());
    for path in paths
    {
        let hashes = select_file(path, cfg)?;
        sigs.push(Signature {
            name: name_from_path(path),
            k: cfg.k as u8,
            scale: cfg.scale,
            kind: cfg.kind,
            hashes,
        });
    }
    Ok(sigs)
}

#[cfg(test)]
mod tests
{
    use super::*;

    #[test]
    fn signature_roundtrips()
    {
        let dir = std::env::temp_dir();
        let path = dir.join(format!("fracsync_sig_{}.sig", std::process::id()));
        let sf = SignatureFile {
            signatures: vec![Signature {
                name: "ref".into(),
                k: 21,
                scale: 1,
                kind: SelectorKind::OpenSyncmer { s: 11, offset: 0 },
                hashes: vec![1, 2, 3, 9, 42],
            }],
        };
        sf.write(&path).unwrap();
        let back = SignatureFile::read(&path).unwrap();
        assert_eq!(back.signatures.len(), 1);
        assert_eq!(back.signatures[0].hashes, vec![1, 2, 3, 9, 42]);
        assert_eq!(back.signatures[0].k, 21);
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn rejects_foreign_file()
    {
        let dir = std::env::temp_dir();
        let path = dir.join(format!("fracsync_bad_{}.sig", std::process::id()));
        std::fs::write(&path, b"not a signature").unwrap();
        assert!(SignatureFile::read(&path).is_err());
        std::fs::remove_file(&path).ok();
    }

    fn fixture() -> SignatureFile
    {
        SignatureFile {
            signatures: vec![Signature {
                name: "ref".into(),
                k: 21,
                scale: 1,
                kind: SelectorKind::OpenSyncmer { s: 11, offset: 0 },
                hashes: vec![1, 2, 3],
            }],
        }
    }

    #[test]
    fn rejects_invalid_signature_invariants_on_read_and_write()
    {
        let path =
            std::env::temp_dir().join(format!("fracsync_invariants_{}.sig", std::process::id()));
        let mut cases = Vec::new();
        for hashes in [vec![3, 1, 2], vec![1, 1, 2]]
        {
            let mut sf = fixture();
            sf.signatures[0].hashes = hashes;
            cases.push(sf);
        }
        for (k, scale, kind) in [
            (0, 1, SelectorKind::Minimizer { w: 1 }),
            (33, 1, SelectorKind::Minimizer { w: 1 }),
            (21, 0, SelectorKind::Minimizer { w: 1 }),
            (21, 1, SelectorKind::Minimizer { w: 0 }),
            (21, 1, SelectorKind::OpenSyncmer { s: 22, offset: 0 }),
            (21, 1, SelectorKind::OpenSyncmer { s: 11, offset: 11 }),
        ]
        {
            let mut sf = fixture();
            let sig = &mut sf.signatures[0];
            sig.k = k;
            sig.scale = scale;
            sig.kind = kind;
            cases.push(sf);
        }
        let mut mixed = fixture();
        let mut second = mixed.signatures[0].clone();
        second.scale = 2;
        mixed.signatures.push(second);
        cases.push(mixed);
        let mut threshold = fixture();
        threshold.signatures[0].scale = 2;
        threshold.signatures[0].hashes = vec![u64::MAX];
        cases.push(threshold);
        for sf in cases
        {
            std::fs::write(&path, b"preserve me").unwrap();
            assert!(sf.write(&path).is_err());
            assert_eq!(std::fs::read(&path).unwrap(), b"preserve me");
            let mut bytes = SIG_MAGIC.to_vec();
            bytes.extend(bincode::encode_to_vec(&sf, bincode::config::standard()).unwrap());
            std::fs::write(&path, bytes).unwrap();
            assert!(SignatureFile::read(&path).is_err());
        }
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn rejects_legacy_truncated_and_trailing_data()
    {
        let path = std::env::temp_dir().join(format!("fracsync_format_{}.sig", std::process::id()));
        fixture().write(&path).unwrap();
        let bytes = std::fs::read(&path).unwrap();
        assert_eq!(&bytes[..8], b"FRACSYN2");
        for len in 0..bytes.len()
        {
            std::fs::write(&path, &bytes[..len]).unwrap();
            assert!(SignatureFile::read(&path).is_err());
        }
        let mut trailing = bytes.clone();
        trailing.push(0);
        std::fs::write(&path, trailing).unwrap();
        assert!(SignatureFile::read(&path).is_err());
        let mut legacy = bytes;
        legacy[..8].copy_from_slice(b"FRACSYN1");
        std::fs::write(&path, legacy).unwrap();
        let error = SignatureFile::read(&path).err().unwrap().to_string();
        assert!(error.contains("re-sketch"), "{error}");
        let empty = SignatureFile { signatures: vec![] };
        empty.write(&path).unwrap();
        assert!(SignatureFile::read(&path).unwrap().config().is_none());
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn selection_deduplicates_files_and_never_joins_records()
    {
        let path = std::env::temp_dir().join(format!("fracsync_records_{}.fa", std::process::id()));
        let cfg = SelectConfig {
            k: 21,
            scale: 1,
            kind: SelectorKind::Minimizer { w: 1 },
        };
        std::fs::write(
            &path,
            b">one\nACGTTGCAACGTACGTAACCGGTTACGT\n>two\nACGTTGCAACGTACGTAACCGGTTACGT\n",
        )
        .unwrap();
        let single = select_file(&path, &cfg).unwrap();
        assert!(!single.is_empty());
        assert_eq!(single, select_files(&[path.clone(), path.clone()], &cfg).unwrap());
        std::fs::write(&path, b">one\nACGTTGCAACGT\n>two\nACGTAACCGGTT\n").unwrap();
        assert!(select_file(&path, &cfg).unwrap().is_empty());
        let invalid = SelectConfig { scale: 0, ..cfg };
        assert!(select_file(&path, &invalid).is_err());
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn rejects_impossible_collection_lengths_without_allocating()
    {
        let path =
            std::env::temp_dir().join(format!("fracsync_lengths_{}.sig", std::process::id()));
        let huge = bincode::encode_to_vec(usize::MAX, bincode::config::standard()).unwrap();
        let mut count = SIG_MAGIC.to_vec();
        count.extend_from_slice(&huge);
        let mut name = SIG_MAGIC.to_vec();
        name.push(1);
        name.extend_from_slice(&huge);
        let mut hashes = SIG_MAGIC.to_vec();
        // One signature; empty name, k=21, scale=1, Minimizer (variant 1), w=1.
        hashes.extend_from_slice(&[1, 0, 21, 1, 1, 1]);
        hashes.extend_from_slice(&huge);
        for bytes in [count, name, hashes]
        {
            std::fs::write(&path, bytes).unwrap();
            let error = SignatureFile::read(&path).err().unwrap();
            assert!(format!("{error:#}").contains("collection length"));
        }
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn streaming_codec_roundtrips_varints_and_unicode()
    {
        let path =
            std::env::temp_dir().join(format!("fracsync_varints_{}.sig", std::process::id()));
        let mut sf = fixture();
        sf.signatures[0].name = "référence".repeat(40);
        sf.signatures[0].hashes = vec![0, 250, 251, 65535, 65536, u32::MAX as u64, u64::MAX];
        sf.signatures = vec![sf.signatures[0].clone(); 300];
        sf.write(&path).unwrap();
        let decoded = SignatureFile::read(&path).unwrap();
        assert_eq!(decoded.signatures.len(), 300);
        for sig in decoded.signatures
        {
            assert_eq!(sig.name, sf.signatures[0].name);
            assert_eq!(sig.hashes, sf.signatures[0].hashes);
        }
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn fastx_parses_wrapped_fasta_fastq_and_concatenated_gzip()
    {
        let path = std::env::temp_dir().join(format!("fracsync_fastx_{}.data", std::process::id()));
        let cfg = SelectConfig {
            k: 21,
            scale: 1,
            kind: SelectorKind::Minimizer { w: 1 },
        };
        let seq = b"ACGTTGCAACGTACGTAACCGGTTACGTTTGGCCAATGCATGCAACGT";
        let mut expected = Vec::new();
        cfg.select(seq, &mut |h| expected.push(h)).unwrap();
        expected.sort_unstable();
        expected.dedup();
        assert!(!expected.is_empty());
        let fasta = format!(
            ">wrapped\r\n{}\r\n{}",
            std::str::from_utf8(&seq[..12]).unwrap(),
            std::str::from_utf8(&seq[12..]).unwrap()
        );
        let fastq = format!(
            "@read\r\n{}\r\n+\r\n{}",
            std::str::from_utf8(seq).unwrap(),
            "I".repeat(seq.len())
        );
        for content in [fasta, fastq]
        {
            std::fs::write(&path, &content).unwrap();
            assert_eq!(select_file(&path, &cfg).unwrap(), expected);
            let mut compressed = Vec::new();
            for part in [content.as_bytes(), b"\n"]
            {
                let mut encoder =
                    flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
                encoder.write_all(part).unwrap();
                compressed.extend(encoder.finish().unwrap());
            }
            // The .data suffix deliberately tests magic-based gzip detection.
            std::fs::write(&path, &compressed).unwrap();
            assert_eq!(select_file(&path, &cfg).unwrap(), expected);
            compressed.truncate(compressed.len() - 4);
            std::fs::write(&path, compressed).unwrap();
            assert!(select_file(&path, &cfg).is_err());
        }
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn fastx_errors_include_input_path()
    {
        let path =
            std::env::temp_dir().join(format!("fracsync_fastx_bad_{}.fq", std::process::id()));
        let cfg = SelectConfig {
            k: 21,
            scale: 1,
            kind: SelectorKind::Minimizer { w: 1 },
        };
        for input in [
            b"".as_slice(),
            b"invalid",
            b"@r\nACGT\n+\n",
            b"@r\nACGT\n-\nIIII\n",
            b"@r\nACGT\n+\nIII\n",
        ]
        {
            std::fs::write(&path, input).unwrap();
            let error = select_file(&path, &cfg).unwrap_err();
            assert!(format!("{error:#}").contains(path.to_str().unwrap()), "{error:#}");
        }
        std::fs::remove_file(path).unwrap();
    }
}
