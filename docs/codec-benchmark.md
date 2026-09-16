# Bounded bitcode prototype: 2026-09-16

**Historical measurements:** this report compares the former `FRACSYN2`
bincode implementation with experimental bitcode codecs. Production now uses
`FRACSYN3` bounded raw bitcode blocks; bincode has been removed. The recorded
numbers below have not been remeasured for the production format.

The current benchmark's `production` baseline uses `FRACSYN3`, so running it now
compares production and prototype bitcode framing; it does not reproduce the old
bincode timings. Delta encoding remains experimental.

## Reproduce

```sh
cargo run --release --example codec_bench -- --out /tmp/fracsync-codecs --repeats 5
cargo run --release --example codec_bench -- --input refs.sig --out /tmp/fracsync-real-codecs --repeats 5
cargo test --example codec_bench
```

The first command creates deterministic input sketches and `.bench` outputs in
the output directory. The second uses an existing `FRACSYN3` database. The program
prints TSV to stdout; [recorded raw measurements](codec-benchmark.tsv) include
file bytes, timing, and Linux memory measurements.

Measured on x86-64 Linux, AMD Ryzen 9 7945HX3D, rustc 1.98.0
(`88d9e12ae`, 2026-08-18), release profile with thin LTO. No CPU pinning or
frequency control was applied.

## Method

At measurement time, the baseline called the **then-production** `SignatureFile::write/read`, including
its validation and incremental bincode reader. Each encode and decode runs in a
new child process. Timing includes file open/close, buffered I/O and validation;
input preparation, fixture generation and round-trip comparison are outside the
timed interval. Every child checks **full equality** of every reference name,
parameter and hash after sampling its peak memory. A mismatch aborts the run.

Values below are median times from five repetitions, with codec order rotated
between repetitions. These are page-cache-backed measurements: writes flush
userspace buffers but do not call `fsync`; reads immediately follow writes.
They measure codec and buffered-I/O performance, not durable disk throughput or
cold-storage latency. Very short timings are especially sensitive to scheduling.

Fixtures use deterministic xorshift-generated data, seed 42:

- **DNA scale 1:** the actual default selector applied to 5 million generated
  A/C/G/T bases, then sorted and deduplicated: 866,794 hashes.
- **DNA scale 10:** the same sequence and selector at scale 10: 86,438 hashes.
- **Many small references:** 1,024 references with 1–127 sorted random hashes
  each (65,060 total). Varying lengths exercise final-block overhead.
- **Large uniform set:** 8 million sorted random 64-bit hashes in one reference.
  This stresses I/O and memory, not the sequence-selection algorithm.

These are generated inputs, not real genomes or a biological accuracy benchmark.
They are useful for the approximately uniform hashes produced by our sampler,
but they do not establish performance on every reference collection.

## Results

Sizes are decimal MB; times are milliseconds.

| Input | Codec | File MB | Write ms | Read ms |
|-------|-------|--------:|---------:|--------:|
| DNA scale 1 | Historical bincode | 7.801 | 5.400 | 13.312 |
| DNA scale 1 | Bitcode raw blocks | 6.938 | 2.643 | 2.676 |
| DNA scale 1 | Bitcode delta blocks | 6.938 | 2.901 | 3.184 |
| DNA scale 10 | Historical bincode | 0.778 | 0.507 | 1.473 |
| DNA scale 10 | Bitcode raw blocks | 0.692 | 0.238 | 0.321 |
| DNA scale 10 | Bitcode delta blocks | 0.692 | 0.259 | 0.358 |
| Many small references | Historical bincode | 0.600 | 0.391 | 1.162 |
| Many small references | Bitcode raw blocks | 0.602 | 0.395 | 0.510 |
| Many small references | Bitcode delta blocks | 0.598 | 0.440 | 0.585 |
| Large uniform set | Historical bincode | 72.000 | 47.368 | 104.136 |
| Large uniform set | Bitcode raw blocks | 64.033 | 26.052 | 22.667 |
| Large uniform set | Bitcode delta blocks | 64.033 | 25.983 | 27.125 |

For the large inputs, raw bitcode was approximately **11% smaller**, **1.8–2.1×
faster to write**, and **4.6–5.0× faster to read**. On many small references it
was **0.37% larger**, essentially tied on write time, and **2.3× faster to read**.
Framing overhead cancels the size saving there.

Much of the large-file size improvement is removing bincode's integer tag:
large random `u64` values occupy nine bytes with the current variable-integer
encoding and about eight bytes in these bitcode blocks. This is not substantial
compression of random hash information, nor evidence that bitcode beats a custom
fixed-width little-endian format; that alternative was not measured.

The delta variant produced byte-for-byte equal **file sizes** to the raw variant
on all three large fixtures (not equal contents). Inspecting bitcode 0.6.9's
integer packing explains why: it uses integer widths 8/16/32/64/128 bits, and
the gaps in these sets still exceed 32 bits. They therefore remain 64-bit values.
Delta reconstruction adds work without reducing those payloads. Tiny blocks can
occasionally benefit, as in the small-reference case, but that saved only 0.23%
relative to bincode and decoded more slowly than raw bitcode.

## Memory

Linux `/proc/self/status` `VmHWM` is sampled before the operation and immediately
after it, before the full-equality check. The TSV reports the maximum sampled
peak across the five runs; write time starts with the source database already
loaded. `write_extra_peak_kib` is the increase over that process's **previous
high-water mark**, not a precise accounting of codec allocations. Allocator
reuse can hide allocations; small differences are noise. On systems without
`/proc/self/status`, the program reports zero for memory, not an actual zero use.

| Large uniform set | Write peak KiB | Additional write high-water KiB | Read peak KiB |
|-------------------|---------------:|-------------------------------:|--------------:|
| Historical bincode | 67,552 | 68 | 67,548 |
| Bitcode raw | 66,284 | 96 | 67,608 |
| Bitcode delta | 67,588 | 96 | 67,540 |

Peak memory is comparable and dominated by the decoded hash vectors. There is
no whole-file encoded copy. The design, rather than these RSS numbers alone,
provides the bound on hash-codec working space: at most 4,096 hashes per block,
a 64 KiB encoded-frame limit, and reusable codec buffers for a finite set of
fixed array types. The complete returned signature collection still grows with
the database.

## Experimental framing and checks

The prototype uses a separate `FSBENCH1` magic, never `FRACSYN3`:

1. Eight-byte magic, one-byte raw/delta mode, little-endian `u64` reference count.
2. Per reference, a little-endian `u32` frame length and bitcode metadata
   `(name, k, scale, selector-tag/parameters, hash-count)`.
3. Per hash block, little-endian `u32` hash count, `u64` base, `u32` payload
   length, then bitcode-encoded fixed-array values. Counts are powers of two,
   at most 4,096. The last remainder is split into smaller exact arrays.
4. Raw mode uses base zero. Delta mode stores the first hash as base, zero as
   the first encoded element, then positive differences within that block.

Decoding bounds frame lengths and hash counts before allocation, checks delta
overflow, rejects truncated/trailing data, and applies the production signature
invariants. Fixed arrays prevent a malicious vector length from expanding a hash
block arbitrarily. These checks are tested, but the experimental reader has not
undergone fuzzing or a complete hostile-input review.

The tests cover block boundaries, empty sketches, Unicode names, extreme hash
values, invalid block lengths, truncated frames and invalid delta arithmetic.
Two checked-in 85-byte fixtures encode one reference named `golden`, with
`k=21`, `scale=1`, endpoint syncmers (`s=11`, `offset=0`) and hashes
`[0, 1, 42, u64::MAX]`. Tests both decode them and reproduce their exact bytes.
They detect accidental changes to this prototype; they are not a commitment to
support its format in the production CLI.

Bitcode's [documented API](https://docs.rs/bitcode/0.6.9/bitcode/) is buffer-based
and does not promise a stable format across major versions. Its
[reusable Buffer](https://docs.rs/bitcode/0.6.9/bitcode/struct.Buffer.html) avoids
reallocating codec scratch space for every frame. The production migration uses
explicit versioning, fixed fixtures, bounded hash blocks, a documented schema,
and re-sketch instructions for older formats. Hash selection is unchanged.
