# Signature format 3

`FRACSYN3` is an eight-byte ASCII magic followed by a little-endian `u64`
signature count. For each signature, write one metadata frame followed by its
hash blocks. No padding or trailing bytes are permitted.

Each frame consists of a little-endian `u32` byte length and exactly that many
bytes encoded with bitcode **0.6.9**. Frame lengths must be 1–65536 bytes.
Metadata encodes this struct, in declaration order:

```rust
struct Metadata {
    name: String,
    k: u8,
    scale: u64,
    kind: SelectorKind,
    hashes: u64,
}
```

`SelectorKind` variants, in declaration order, are
`OpenSyncmer { s: u8, offset: u8 }`, `Minimizer { w: u8 }`, and `All`.
Names are limited to 4096 UTF-8 bytes. Selection parameters must be valid and
identical across signatures in a file.

For each hash block, choose the largest power of two no greater than both the
remaining hash count and 4096. Encode exactly `[u64; N]` as one bitcode frame.
Repeat until all hashes are written. Empty signatures have no hash blocks.
Hash counts and block sizes are inferred from metadata; blocks have no extra
count, base, padding, or delta field. Hashes must be strictly increasing and
no greater than `u64::MAX / scale`.

Writers buffer only one block of encoded hashes at a time. Readers decode fixed
arrays and grow the result incrementally, without allocating a hash vector from
an untrusted count. Metadata is frame bounded; total decoded signature storage
still grows with the file's contents. This format provides validation, not a
checksum or authentication.

Version 3 preserves version 2's hash values but changes the encoding. Versions
1 and 2 require re-sketching. Codec/schema changes require a format-version
review and compatible fixture checks. NCSI uses the same block primitives with
its own magic and metadata (including group and marker fields); its signatures
are not interchangeable with fracsync signatures.
