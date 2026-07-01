# Changelog

All notable changes to DIMPLE are documented here. The format loosely follows
[Keep a Changelog](https://keepachangelog.com/); versions above the published
`v1.0` manuscript snapshot begin at `1.1.0`.

## [1.1.0] - 2026-06-30

Reliability and cross-entrypoint consistency release. No change to the core
oligo-design algorithm; the shipped `maxfrag` formula and default codon choices
are unchanged. Regression goldens were regenerated (see *maxfrag* below).

### Changed
- **Unified defaults across all four entrypoints (CLI, GUI, web, notebook).**
  Previously the "same" defaults produced different physical oligos depending on
  how you launched DIMPLE.
  - **Oligo length** now defaults to `250` everywhere (was `230` on CLI/web/
    notebook, `250` in the GUI). This is the *total synthesized* oligo length.
  - **Random seed** is now deterministic-by-default: `1848` everywhere, so the
    same gene and settings yield identical oligos on every entrypoint. Pass any
    integer (including `0`) to vary the library. The effective seed is logged on
    every run.
- **DMS is now on by default** in the GUI (it already was in the web app and
  notebook). The CLI keeps `-DMS` as an explicit opt-in for scripted use.
- **`-include_substitutions` is now an alias of `-DMS`** on the CLI (it was
  declared but never read). Both set the DMS toggle.
- **`requires-python`** capped to `>=3.12,<3.13` and a `.python-version` (3.12)
  committed, so installs no longer resolve onto the biopython-incompatible 3.13.

### Fixed
- **Non-DMS runs no longer emit an empty `mutations.csv`.** The file is written
  only for DMS runs (restoring pre-1.0.x behavior); DIS/insertion/deletion-only
  runs produce no spurious zero-byte file.
- **GUI runs now complete.** The GUI passed a `config=` keyword to
  `generate_DMS_fragments`, `post_qc`, and `print_all` — none of which accept it
  — so every real run raised `TypeError` after fragment generation. Removed;
  these read `pool.config` internally. The DMS checkbox also now actually
  launches selected (a `deselect()` at widget creation had been overriding the
  default-on behavior).
- **GUI pipeline errors now surface.** The whole run pipeline is wrapped so
  errors (e.g. the common "unwanted restriction site" `ValueError`) appear in a
  message box and the output log instead of vanishing to stderr behind a blank
  window. A full-pipeline GUI smoke test now guards against this class of
  call-signature drift.
- **Custom mutations no longer silently drop residues.** Overlapping range and
  single-position lines now accumulate and de-duplicate consistently regardless
  of line order (previously a range could overwrite an earlier single).
- **Deletion "extends beyond ORF" warning now fires correctly.** The guard
  compared mismatched units (amino-acid index vs nucleotide count vs a
  primer-buffer-inflated length) and never triggered; it is now unit-correct and
  advisory (the variant is still emitted).
- **`align_genevariation` (homolog linking, `-matchSequences match`) no longer
  crashes at the default overlap.** Its shared-breaklist snap grid was aligned
  to multiples of 3 while the breaksites setter requires alignment relative to
  `primer_buffer`; the two disagreed at overlap=4, raising a `ValueError`. Now
  frame-aligned at any overlap.
- **`doublefrag=1`** now raises a clear `NotImplementedError` at config time
  instead of a cryptic mid-run `AttributeError`. The two-fragments-per-oligo
  layout remains unsupported pending a real use case.

### Tests
- Regression tests aligned to the shipped `maxfrag = synth_len - 64 -
  overlap_l - overlap_r` formula at the real default overlap (4, was 3), and the
  byte goldens regenerated. Oligo and variant counts are unchanged; fragment
  boundaries (and therefore oligo sequences, barcodes, and primer choices)
  shifted to match what the entrypoints actually emit.

### Security
- **biopython `Bio.Entrez` XXE advisory (CVE-2025-68463 / GHSA-x3vf-39hj-gxr4)
  does not affect DIMPLE.** The codebase never imports `Bio.Entrez`, so the
  vector is unreachable. The `biopython==1.84` pin is retained for this release:
  upgrading past the advisory requires a coupled `pydna` major-minor upgrade
  (pydna 5.4.0 caps `biopython<1.85`; only pydna ≥5.5.x permits ≥1.87), which
  touches the QC amplification path and is deferred to a dedicated dependency
  pass.
