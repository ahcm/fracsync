//! **FracSync** — strand-symmetric syncmers, FracMinHash sampling and containment.
//!
//! Developed from the selection core of NCSI, an unpublished manuscript and
//! unreleased classifier. See the crate README for method references and scope.
//!
//! - [`hash`] — canonical ntHash-style rolling hashes (raw, unmixed values).
//! - [`select`] — mirrored syncmers (closed endpoints by default) or
//!   window-dependent minimizers, low-complexity masking and mixed-hash sampling.
//! - [`sketch`] — validated version-2 signatures and record-at-a-time sketching.
//! - [`contain`] — selected-set containment and an ANI point estimate without
//!   coverage correction or confidence intervals.
//!
//! Syncmers reselect intact k-mers independently of flanking sequence;
//! minimizers do not. Scale D retains approximately 1/D of selected hashes.
//! Selectors have bounded working state, but records and distinct hash sets
//! still consume memory proportional to their size.
//!
//! ```
//! use fracsync::select::{SelectConfig, SelectorKind};
//! let cfg = SelectConfig {
//!     k: 21, scale: 1,
//!     kind: SelectorKind::OpenSyncmer { s: 11, offset: 0 }, // closed syncmers
//! };
//! let mut hashes = Vec::new();
//! cfg.select(b"ACGTTGCAACGTACGTAACCGGTTACGT", &mut |h| hashes.push(h))?;
//! hashes.sort_unstable();
//! hashes.dedup();
//! # Ok::<(), String>(())
//! ```

pub mod contain;
pub mod hash;
pub mod select;
pub mod sketch;
