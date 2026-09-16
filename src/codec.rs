//! Wire primitives for bitcode 0.6.9: little-endian u64 counts, u32-length
//! frames, and greedy power-of-two arrays of at most 4096 raw u64 hashes.
//! Callers supply format magic, metadata schema, and semantic validation.
use anyhow::{Context, Result, ensure};
use std::io::{Read, Write};

/// Maximum UTF-8 byte length of a reference or group name.
pub const MAX_NAME: usize = 4096;
const MAX_FRAME: usize = 64 * 1024;
const BLOCK: usize = 4096;

pub fn write_count(out: &mut impl Write, count: u64) -> Result<()>
{
    out.write_all(&count.to_le_bytes())?;
    Ok(())
}
pub fn read_count(input: &mut impl Read) -> Result<u64>
{
    let mut bytes = [0; 8];
    input.read_exact(&mut bytes)?;
    Ok(u64::from_le_bytes(bytes))
}
pub fn write_frame(out: &mut impl Write, bytes: &[u8]) -> Result<()>
{
    ensure!(!bytes.is_empty() && bytes.len() <= MAX_FRAME, "invalid frame size");
    out.write_all(&(bytes.len() as u32).to_le_bytes())?;
    out.write_all(bytes)?;
    Ok(())
}
pub fn read_frame(input: &mut impl Read) -> Result<Vec<u8>>
{
    let mut size = [0; 4];
    input.read_exact(&mut size)?;
    let size = u32::from_le_bytes(size) as usize;
    ensure!(size > 0 && size <= MAX_FRAME, "invalid frame size");
    let mut bytes = vec![0; size];
    input.read_exact(&mut bytes)?;
    Ok(bytes)
}
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
fn encode_block<const N: usize>(
    out: &mut impl Write,
    hashes: &[u64],
    codec: &mut bitcode::Buffer,
) -> Result<()>
{
    let values: &[u64; N] = hashes.try_into()?;
    write_frame(out, codec.encode(values))
}
fn decode_block<const N: usize>(
    bytes: &[u8],
    hashes: &mut Vec<u64>,
    codec: &mut bitcode::Buffer,
) -> Result<()>
{
    let values: [u64; N] = codec.decode(bytes).context("decoding hash block")?;
    hashes.extend_from_slice(&values);
    Ok(())
}
/// Encodes without copying the full hash vector into a second buffer.
pub fn write_hashes(out: &mut impl Write, mut hashes: &[u64]) -> Result<()>
{
    let mut codec = bitcode::Buffer::new();
    while !hashes.is_empty()
    {
        let count = 1usize << hashes.len().min(BLOCK).ilog2();
        dispatch!(count, encode_block, out, &hashes[..count], &mut codec)?;
        hashes = &hashes[count..];
    }
    Ok(())
}
/// Grows the output only after successfully decoding bounded, fixed-size arrays.
pub fn read_hashes(input: &mut impl Read, mut count: u64) -> Result<Vec<u64>>
{
    let mut hashes = Vec::new();
    let mut codec = bitcode::Buffer::new();
    while count > 0
    {
        let n = 1usize << count.min(BLOCK as u64).ilog2();
        let bytes = read_frame(input)?;
        dispatch!(n, decode_block, &bytes, &mut hashes, &mut codec)?;
        count -= n as u64;
    }
    Ok(hashes)
}
pub fn finish(input: &mut impl Read) -> Result<()>
{
    ensure!(input.read(&mut [0])? == 0, "trailing data after signatures");
    Ok(())
}

#[cfg(test)]
mod tests
{
    use super::*;
    #[test]
    fn roundtrips_block_boundaries()
    {
        for n in [0, 1, 2, 3, 4095, 4096, 4097, 8193]
        {
            let hashes: Vec<u64> = (0..n).map(|i| u64::MAX / n.max(1) * i).collect();
            let mut bytes = Vec::new();
            write_hashes(&mut bytes, &hashes).unwrap();
            let mut input = bytes.as_slice();
            assert_eq!(read_hashes(&mut input, n).unwrap(), hashes);
            finish(&mut input).unwrap();
        }
    }
    #[test]
    fn rejects_bad_frames_and_counts()
    {
        for size in [0u32, 65537, u32::MAX]
        {
            assert!(read_frame(&mut size.to_le_bytes().as_slice()).is_err());
        }
        assert!(read_hashes(&mut [].as_slice(), u64::MAX).is_err());
        let mut bytes = Vec::new();
        write_hashes(&mut bytes, &[42]).unwrap();
        assert!(read_hashes(&mut bytes.as_slice(), 4096).is_err());
    }
}
