# fracsync

Strand-symmetric **syncmer selection + FracMinHash sampling**, with a Rust library
and CLI for sketching FASTA/FASTQ references and measuring their containment in a
sample. The default selects **closed syncmers** (`k=21, s=11, offset=0`).

FracSync originated as the selection and containment substrate of my **NCSI**
project for classifying pathogens. The NCSI manuscript is **unpublished** and its
code is **unreleased**. This crate implements the reusable core, with the engineering
changes described below; it is not the full NCSI classifier. Syncmers, minimizers,
ntHash and FracMinHash are established techniques, credited below.

## Build and use

```sh
cargo build --release
cargo test

# One reference per input file; all records in that file form one reference.
fracsync sketch refs/*.fa -o refs.sig

# All query files/records form one sample.
fracsync contain refs.sig reads.fastq.gz
fracsync contain refs.sig reads1.fastq.gz reads2.fastq.gz --min-containment 0.05

# Retain approximately 1/10 of the syncmer-selected distinct hashes.
fracsync sketch refs/*.fa -o refs.sig --scale 10

# Central open syncmers: k-s=10, so the mirror of offset 5 is itself.
fracsync sketch refs/*.fa -o open.sig -k 21 --s 11 --offset 5

# Window-dependent minimizers, for comparison.
fracsync sketch refs/*.fa -o min.sig --selector minimizer --w 11

# Every k-mer under the scale threshold (plain FracMinHash), for comparison:
# the same content-defined robustness as a syncmer, no spacing guarantee.
fracsync sketch refs/*.fa -o all.sig --selector all --scale 11

fracsync info refs.sig
```

Input is parsed with [fastx](https://docs.rs/fastx/0.6.1/), using reusable record
buffers and a 64 KiB input buffer. Wrapped FASTA and four-line FASTQ are supported;
gzip (including concatenated members) is detected by content, not filename. Bases are
case-insensitive and `U` is treated as `T` in both hashing and complexity filtering.
K-mers containing other bases are skipped. Records are never concatenated across
boundaries. Parameters must satisfy `1 <= k <= 32`, `scale >= 1`, and either
`1 <= s < k`, `0 <= offset <= k-s`, or `1 <= w <= 255` for minimizers.

`contain` writes TSV, sorted by decreasing containment, then decreasing
`sample_frac`:

| Column | Meaning |
|--------|---------|
| `reference` | Input filename with common sequence/compression extensions removed |
| `ref_kmers` | Number of distinct hashes in the reference **sketch**, not all genomic k-mers |
| `shared` | Size of the reference/sample sketch intersection |
| `containment` | `shared / ref_kmers` |
| `ani` | Uncorrected point estimate `containment^(1/k)` |
| `sample_frac` | `shared / sample_kmers`, using the distinct hashes in the sample sketch |

Containment measures how much of the reference's **selected, complexity-filtered
content** appears in the sample. It is not automatically the fraction of all bases
or all k-mers covered. `sample_frac` is not abundance or unique attribution:
related references may explain overlapping hashes.

The ANI transform assumes that intact k-mer survival reflects nucleotide identity.
Incomplete coverage and sequencing errors also lower containment, so this value
can substantially underestimate identity. There is no depth correction, mutation
confidence interval, or correction for empty-sketch bias. Empty reference/query
sketches produce zero for the corresponding ratios; zero from an empty sketch is
not evidence of biological absence. Small sketches have high sampling uncertainty.
These limitations matter particularly at large scales. See Blanca et al., Hera
et al. and Sylph in the references below.

## How our syncmers relate to the literature

[Edgar (2021)](https://doi.org/10.7717/peerj.10805) introduced syncmers: select a
k-mer according to where its smallest constituent s-mer occurs. This decision
uses only the k-mer itself. An intact selected k-mer is therefore selected again
in another read, even if its flanking sequence has changed.

FracSync orders s-mers by their canonical ntHash-style value and accepts the
zero-based offset pair:

```text
{offset, k - s - offset}
```

Reverse complementation reverses the s-mer positions, while their canonical
hashes remain the same. Accepting mirrored offsets makes the decision
strand-symmetric. Equal minimum hashes prefer the position nearest the centre;
if two such positions remain tied, they are mirrors and have the same acceptance
decision. This is the rule inherited from the NCSI draft.

| Offsets (`k=21, s=11`) | Relation to published schemes |
|------------------------|--------------------------------|
| `{0, 10}` (default) | Edgar's **closed** endpoint scheme, with our tie policy |
| `{5}` | Central **open** syncmer |
| `{o, 10-o}` | A symmetric two-offset instance of parameterized syncmers |

Per-k-mer conservation under error is a property of any rule that decides from
the k-mer alone, so a syncmer and plain FracMinHash (`--selector all`) at the
same density reselect the same intact k-mers; what the closed syncmer adds is
bounded spacing, which matters for short, noisy reads. The library exposes
`SelectConfig::expected_density()` — selector density times `1/D` — so a
chance model can be sized by the universe of *selectable* k-mers.

[Dutta, Pellow and Shamir (2022)](https://doi.org/10.1371/journal.pcbi.1010638)
study parameterized syncmer schemes and how offset choices affect spacing and
conservation. FracSync exposes a mirrored pair, not an arbitrary offset set or
their complete mapping implementation. The historical Rust variant name
`SelectorKind::OpenSyncmer` remains for compatibility; it covers both open and
closed settings.

With random, distinct s-mer orderings, the endpoint pair has approximate density
`2/(k-s+1)` and a single central position has approximate density `1/(k-s+1)`.
Classical closed syncmers have a window guarantee. **The final FracSync sketch
has no unconditional maximum-gap guarantee**: our centre-preferring tie policy
can suppress endpoint selections, low-complexity masking removes candidates,
and fractional sampling can remove any remaining selection. Ambiguous bases
also break valid sequence runs.

[Shaw and Yu (2022)](https://doi.org/10.1093/bioinformatics/btab790) provide the
local-selection framework underlying these conservation arguments. Reselecting
an intact k-mer is not unique to syncmers: plain hash-threshold sampling has that
property too. Syncmers provide a way to shape spacing before thinning; we do not
claim superior per-k-mer conservation to plain FracMinHash or transfer NCSI's
draft benchmark numbers to this crate.

Minimizers follow [Roberts et al. (2004)](https://doi.org/10.1093/bioinformatics/bth408):
choose the smallest raw canonical k-mer hash in each window of `w` valid k-mers,
with leftmost tie-breaking. They depend on neighbouring k-mers and read
boundaries. A record shorter than `k+w-1` emits no minimizers, and cutting a
reference into reads can change which k-mers are selected. The syncmer
boundary-independence claim does not apply to this comparison mode.

## Hashing, scaling and compatibility

The rolling recurrence follows [ntHash (Mohamadi et al., 2016)](https://doi.org/10.1093/bioinformatics/btw397).
The raw canonical value is `min(forward, reverse)`. This minimum is not uniform:
thresholding it directly retains roughly `1-(1-1/D)^2`, about `2/D` for large `D`,
under the independent-uniform approximation. The NCSI draft already identifies
this issue; the original FracSync format used that rule.

This crate applies the **SplitMix64 bijective finalizer** to the canonical
k-mer value before emitting or thresholding it. The mixed value is retained when
`hash <= u64::MAX / D`, giving approximately `1/D` retention in the tested random
sequence cases. Mixing preserves reverse-complement equality and introduces no
additional hash collisions. It is not a proof of independent sampling for every
sequence distribution. Selector ordering still uses the raw hashes, so changing
scale only removes selected hashes; it does not change the underlying selector.

The fractional-threshold approach follows the sourmash/FracMinHash work of
[Brown and Irber (2016)](https://doi.org/10.21105/joss.00027) and
[Irber et al. (2022)](https://doi.org/10.1101/2022.01.11.475838).
Under ideal uniform sampling, conditioning on a nonempty reference sketch gives
an unbiased estimate of containment of the underlying selected sets. Reporting
zero for empty sketches introduces bias, as discussed by
[Hera et al. (2023)](https://doi.org/10.1101/gr.277651.123).

**Signature format is now `FRACSYN2`. Re-sketch existing references.**
`FRACSYN1` files are rejected with an explicit migration message, because stored
hashes and scale semantics changed; merely changing the header is not a valid
conversion. This crate does not promise hash or file compatibility with NCSI,
sourmash, or other ntHash implementations. Reads and writes validate selection
parameters, parameter agreement across references, strictly increasing hashes,
and the scale threshold. Readers also reject truncated files and trailing data.

## Memory and runtime

Selectors use rolling hashes and bounded monotone queues: `O(k-s+1)` working
space for syncmers and `O(w)` for minimizers, independent of record length.
Minima are maintained in amortized linear time; syncmer ties scan at most
`k-s+1` entries, and complexity checks inspect at most `k <= 32` bases.
There are no full-record k-mer/s-mer hash arrays in the selection path.

The parser still buffers sequence records, and sketching keeps distinct selected
hashes in a hash set before sorting. The CLI retains reference signatures in RAM;
query memory grows with distinct selected sample hashes. Thus this is
**record-at-a-time processing, not a fixed-memory or out-of-core pipeline**.
Increasing scale reduces the expected hash-set size but does not cap memory or
shrink the parser's largest-record buffer.

Signature I/O is buffered without whole-file encoding copies. Query files share
one deduplication set. Intersections use a linear merge for similarly sized sets
and binary searches when one set is much smaller, avoiding a complete scan of a
large sample for every tiny reference. There is no inverted reference index or
parallel query engine here.

## Experimental signature-codec benchmark

A separate prototype compares the current `FRACSYN2` bincode I/O with bounded
bitcode blocks, with and without hash deltas. Bitcode is pinned to 0.6.9 as a
**development-only dependency**; the library and CLI still use `FRACSYN2`.

```sh
cargo run --release --example codec_bench -- --out /tmp/fracsync-codecs --repeats 5
# Or benchmark a real signature database:
cargo run --release --example codec_bench -- --input refs.sig --out /tmp/fracsync-codecs
cargo test --example codec_bench
```

The [benchmark report](docs/codec-benchmark.md) records the method, raw results,
limitations and recommendation. On the generated large sketches, raw bitcode
blocks were about 11% smaller and substantially faster to read. Delta encoding
gave no additional size reduction on those large cases. Small-reference framing
overhead matters, so these are not universal compression claims. Prototype
`.bench` files are experimental and are not accepted by the production CLI.

## As a library

```rust
use fracsync::select::{SelectConfig, SelectorKind};

let cfg = SelectConfig {
    k: 21,
    scale: 10,
    kind: SelectorKind::OpenSyncmer { s: 11, offset: 0 }, // closed endpoints
};
let seq = b"ACGTTGCAACGTACGTAACCGGTTACGT";
let mut hashes = Vec::new();
cfg.select(seq, &mut |h| hashes.push(h)).expect("valid selection parameters");
hashes.sort_unstable();
hashes.dedup();
```

`SelectConfig::select` now returns `Result<(), String>` and validates before
emitting. Callers must handle that result. The callback can receive duplicates.
`hash::for_each_canonical` still exposes the raw ntHash-style values; those are
not the mixed values stored in version-2 signatures. Use `select` to reproduce
a signature's hashes.

| Module | Role |
|--------|------|
| `hash` | Canonical rolling ntHash-style hashing |
| `select` | Mirrored syncmers / minimizers, complexity masking, mixed-hash sampling |
| `sketch` | Validated signature I/O and FASTA/FASTQ sketching |
| `contain` | Sorted-set intersection, containment and uncorrected ANI |

Library scoring expects sorted, deduplicated sets produced under the same
configuration. `SignatureFile::read` enforces file invariants; callers assembling
raw slices themselves must uphold those preconditions.

## Relation to NCSI and references

The source work is *NCSI: containment-based open-set sequence
identification for noisy long reads*, unpublished manuscript draft dated
2026-09-11, with unreleased accompanying code. FracSync extracts and develops its
selection/containment core. The draft's count-spectrum depth model,
coverage-corrected ANI, abundance estimates, open-set reporting, mmap index,
contamination masking and strain-marker tier are outside this crate. In
particular, we do not implement the coverage correction discussed in
[Sylph (Shaw and Yu)](https://doi.org/10.1038/s41587-024-02412-y).

Relevant references from the manuscript, grouped by their role here:

1. **Syncmers:** Edgar R. (2021). *Syncmers are more sensitive than minimizers for
   selecting conserved k-mers in biological sequences.* PeerJ 9:e10805.
   [DOI](https://doi.org/10.7717/peerj.10805).
2. **Offset generalization:** Dutta A., Pellow D., Shamir R. (2022).
   *Parameterized syncmer schemes improve long-read mapping.* PLOS Computational
   Biology 18(10):e1010638. [DOI](https://doi.org/10.1371/journal.pcbi.1010638).
3. **Selection theory:** Shaw J., Yu Y.W. (2022). *Theory of local k-mer selection
   with applications to long-read alignment.* Bioinformatics 38(20):4659–4669.
   [DOI](https://doi.org/10.1093/bioinformatics/btab790).
4. **Minimizers:** Roberts M. et al. (2004). *Reducing storage requirements for
   biological sequence comparison.* Bioinformatics 20(18):3363–3369.
   [DOI](https://doi.org/10.1093/bioinformatics/bth408).
5. **Rolling hash:** Mohamadi H. et al. (2016). *ntHash: recursive nucleotide
   hashing.* Bioinformatics 32(22):3492–3494.
   [DOI](https://doi.org/10.1093/bioinformatics/btw397).
6. **Sketching software:** Brown C.T., Irber L. (2016). *sourmash: a library for
   MinHash sketching of DNA.* JOSS 1(5):27.
   [DOI](https://doi.org/10.21105/joss.00027).
7. **FracMinHash:** Irber L. et al. (2022). *Lightweight compositional analysis of
   metagenomes with FracMinHash and minimum metagenome covers.* bioRxiv preprint.
   [DOI](https://doi.org/10.1101/2022.01.11.475838).
8. **Mutation model:** Blanca A. et al. (2022). *The statistics of k-mers from a
   sequence undergoing a simple mutation process without spurious matches.*
   Journal of Computational Biology 29(2):155–168.
   [DOI](https://doi.org/10.1089/cmb.2021.0431).
9. **FracMinHash statistics and ANI:** Hera M.R., Pierce-Ward N.T., Koslicki D.
   (2023). *Deriving confidence intervals for mutation rates across a wide range
   of evolutionary distances using FracMinHash.* Genome Research.
   [DOI](https://doi.org/10.1101/gr.277651.123).
10. **Coverage correction (not implemented):** Shaw J., Yu Y.W. (online 2024;
    volume publication 2025). *Rapid species-level metagenome profiling and
    containment estimation with sylph.* Nature Biotechnology 43:1348–1359.
    [DOI](https://doi.org/10.1038/s41587-024-02412-y).

## License

EUPL-1.2
