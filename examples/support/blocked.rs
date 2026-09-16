//! Experimental format for codec measurements, not a supported signature format.
//! Bitcode 0.6.9; length-prefixed metadata and <=4096-hash blocks. Arrays bound
//! decoded hash allocation even if a block contains hostile compressed values.
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;

use anyhow::{Context, Result, ensure};
use fracsync::select::SelectorKind;
use fracsync::sketch::{Signature, SignatureFile};

const MAGIC: &[u8; 8] = b"FSBENCH1";
const BLOCK: usize = 4096;
const MAX_FRAME: usize = 64 * 1024;

#[derive(bitcode::Encode, bitcode::Decode)]
struct Metadata
{
    name: String,
    k: u8,
    scale: u64,
    kind: [u8; 3],
    hashes: u64,
}

// Use power-of-two arrays to bound decoding. Split the final remainder into
// smaller blocks so short references never pay for zero padding.
macro_rules! dispatch {
    ($n:expr, $f:ident, $($arg:expr),*) => {
        match $n {
            1 => $f::<1>($($arg),*), 2 => $f::<2>($($arg),*),
            4 => $f::<4>($($arg),*), 8 => $f::<8>($($arg),*),
            16 => $f::<16>($($arg),*), 32 => $f::<32>($($arg),*),
            64 => $f::<64>($($arg),*), 128 => $f::<128>($($arg),*),
            256 => $f::<256>($($arg),*), 512 => $f::<512>($($arg),*),
            1024 => $f::<1024>($($arg),*), 2048 => $f::<2048>($($arg),*),
            4096 => $f::<4096>($($arg),*), _ => unreachable!(),
        }
    };
}

fn put_frame(out: &mut impl Write, bytes: &[u8]) -> Result<()>
{
    ensure!(bytes.len() <= MAX_FRAME, "frame too large");
    out.write_all(&(bytes.len() as u32).to_le_bytes())?;
    out.write_all(bytes)?;
    Ok(())
}

fn get_frame(input: &mut impl Read, bytes: &mut Vec<u8>) -> Result<()>
{
    let mut size = [0; 4];
    input.read_exact(&mut size)?;
    let size = u32::from_le_bytes(size) as usize;
    ensure!(size > 0 && size <= MAX_FRAME, "invalid frame size");
    bytes.resize(size, 0);
    input.read_exact(bytes)?;
    Ok(())
}

fn encode_block<'a, const N: usize>(
    hashes: &[u64],
    delta: bool,
    codec: &'a mut bitcode::Buffer,
) -> &'a [u8]
{
    let mut values = [0u64; N];
    if delta
    {
        // First hash is stored outside the codec; each block is independent.
        for (i, pair) in hashes.windows(2).enumerate()
        {
            values[i + 1] = pair[1] - pair[0];
        }
    }
    else
    {
        values[..hashes.len()].copy_from_slice(hashes);
    }
    codec.encode(&values)
}

fn decode_block<const N: usize>(
    bytes: &[u8],
    codec: &mut bitcode::Buffer,
    count: usize,
    base: u64,
    delta: bool,
    hashes: &mut Vec<u64>,
) -> Result<()>
{
    let mut values: [u64; N] = codec.decode(bytes).context("decoding hash block")?;
    ensure!(values[count..].iter().all(|&v| v == 0), "nonzero padding");
    if delta
    {
        ensure!(values[0] == 0, "invalid first delta");
        values[0] = base;
        for i in 1..count
        {
            ensure!(values[i] > 0, "zero delta");
            values[i] = values[i - 1]
                .checked_add(values[i])
                .context("delta overflow")?;
        }
    }
    else
    {
        ensure!(base == 0, "unexpected raw block base");
    }
    hashes.extend_from_slice(&values[..count]);
    Ok(())
}

pub fn write(sf: &SignatureFile, path: &Path, delta: bool) -> Result<()>
{
    sf.validate()?;
    let mut out = BufWriter::new(File::create(path)?);
    out.write_all(MAGIC)?;
    out.write_all(&[u8::from(delta)])?;
    out.write_all(&(sf.signatures.len() as u64).to_le_bytes())?;
    let mut codec = bitcode::Buffer::new();
    for sig in &sf.signatures
    {
        let kind = match sig.kind
        {
            SelectorKind::OpenSyncmer { s, offset } => [0, s, offset],
            SelectorKind::Minimizer { w } => [1, w, 0],
        };
        ensure!(sig.name.len() <= 4096, "benchmark reference name too long");
        let metadata = Metadata {
            name: sig.name.clone(),
            k: sig.k,
            scale: sig.scale,
            kind,
            hashes: sig.hashes.len() as u64,
        };
        put_frame(&mut out, codec.encode(&metadata))?;
        let mut rest = sig.hashes.as_slice();
        while !rest.is_empty()
        {
            let available = rest.len().min(BLOCK);
            let count = 1usize << available.ilog2();
            let (hashes, tail) = rest.split_at(count);
            rest = tail;
            out.write_all(&(hashes.len() as u32).to_le_bytes())?;
            out.write_all(&if delta { hashes[0] } else { 0 }.to_le_bytes())?;
            let encoded = dispatch!(hashes.len(), encode_block, hashes, delta, &mut codec);
            put_frame(&mut out, encoded)?;
        }
    }
    out.flush()?;
    Ok(())
}

pub fn read(path: &Path) -> Result<SignatureFile>
{
    let mut input = BufReader::new(File::open(path)?);
    let mut header = [0; 17];
    input.read_exact(&mut header)?;
    ensure!(&header[..8] == MAGIC && header[8] <= 1, "invalid prototype header");
    let delta = header[8] == 1;
    let count = u64::from_le_bytes(header[9..].try_into()?);
    let mut signatures = Vec::new();
    let mut codec = bitcode::Buffer::new();
    let mut bytes = Vec::new();
    for _ in 0..count
    {
        get_frame(&mut input, &mut bytes)?;
        let m: Metadata = codec.decode(&bytes).context("decoding metadata")?;
        ensure!(m.name.len() <= 4096, "reference name too long");
        let kind = match m.kind
        {
            [0, s, offset] => SelectorKind::OpenSyncmer { s, offset },
            [1, w, 0] => SelectorKind::Minimizer { w },
            _ => anyhow::bail!("invalid selector"),
        };
        let mut hashes = Vec::new();
        while (hashes.len() as u64) < m.hashes
        {
            let mut block = [0; 12];
            input.read_exact(&mut block)?;
            let n = u32::from_le_bytes(block[..4].try_into()?) as usize;
            ensure!(
                n.is_power_of_two() && n <= BLOCK && n as u64 <= m.hashes - hashes.len() as u64,
                "invalid block count"
            );
            let base = u64::from_le_bytes(block[4..].try_into()?);
            get_frame(&mut input, &mut bytes)?;
            dispatch!(n, decode_block, &bytes, &mut codec, n, base, delta, &mut hashes)?;
        }
        signatures.push(Signature {
            name: m.name,
            k: m.k,
            scale: m.scale,
            kind,
            hashes,
        });
    }
    ensure!(input.read(&mut [0])? == 0, "trailing data");
    let sf = SignatureFile { signatures };
    sf.validate()?;
    Ok(sf)
}

#[cfg(test)]
mod tests
{
    use super::*;

    #[test]
    fn rejects_truncated_oversized_and_invalid_frames()
    {
        let path =
            std::env::temp_dir().join(format!("fracsync_blocked_bad_{}.bench", std::process::id()));
        let sf = SignatureFile {
            signatures: vec![Signature {
                name: "test".into(),
                k: 21,
                scale: 1,
                kind: SelectorKind::OpenSyncmer { s: 11, offset: 0 },
                hashes: vec![0, 1, 42],
            }],
        };
        write(&sf, &path, true).unwrap();
        let bytes = std::fs::read(&path).unwrap();
        for len in 0..bytes.len()
        {
            std::fs::write(&path, &bytes[..len]).unwrap();
            assert!(read(&path).is_err(), "accepted truncation at {len}");
        }
        let mut trailing = bytes.clone();
        trailing.push(0);
        std::fs::write(&path, trailing).unwrap();
        assert!(read(&path).is_err());
        let block = 21 + u32::from_le_bytes(bytes[17..21].try_into().unwrap()) as usize;
        for (offset, value) in [
            (17, u32::MAX),
            (block, 0),
            (block, 3),
            (block, 8192),
            (block + 12, u32::MAX),
        ]
        {
            let mut bad = bytes.clone();
            bad[offset..offset + 4].copy_from_slice(&value.to_le_bytes());
            std::fs::write(&path, bad).unwrap();
            assert!(read(&path).is_err());
        }
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn rejects_invalid_delta_arithmetic()
    {
        let mut codec = bitcode::Buffer::new();
        let mut hashes = Vec::new();
        for (values, base) in [([1, 1], 0), ([0, 0], 0), ([0, u64::MAX], 1)]
        {
            let bytes = bitcode::encode(&values);
            assert!(decode_block::<2>(&bytes, &mut codec, 2, base, true, &mut hashes).is_err());
            assert!(hashes.is_empty());
        }
    }
}
