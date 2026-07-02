# verification/

One-off scripts to independently validate DIMPLE library output (not part of the
test suite). Built 2026-07-01 to check that the designed oligos actually encode
the intended mutations and assemble correctly. Untracked (git-excluded).

## Scripts

### `validate_library.py` — independent oligo-level check
Reads `All_Oligos.fasta` + the WT gene FASTA (independent of DIMPLE's own CSVs),
digests each oligo with the Type IIS enzyme, reconstructs the mutant CDS by
flank-offset splicing, translates, and compares the observed edit to the intent
parsed from the oligo *name*.

```
python verification/validate_library.py \
    --oligos tests/expected/All_Oligos.fasta \
    --gene   tests/expected/Kir.fa \
    --enzyme BsmBI
```

Reconstructs substitutions, insertions and deletions. Confirms edits away from
subpool edges cleanly; edits *at* a subpool 5'/3' edge have one WT flank too
short to reconstruct from the oligo alone (their context lives across the
assembly junction) — those are covered by the pydna cross-check instead.

### `build_pool.py` — rebuild the Kir DMS pool
Rebuilds the pool exactly like `tests/regression/test_dms_pipeline_kir.py` but
with `enzyme="BsmBI"` so the assembly machinery is available. Import `build()`.

### `pydna_crosscheck.py` — assembly-level cross-check (method 2)
For each variant: PCR the WT gene with its subpool gene primers -> BsmBI cut ->
backbone; PCR the oligo with its oligo primers -> cut -> insert; ligate ->
circular product; translate the ORF (coordinate-free, best-ORF by block-diff);
compare to intent. Uses pydna's real sticky-end ligation, so it is independent
of `validate_library.py`. Runs the method-1 residual + a sample of passing
variants and reports agreement.

```
python verification/pydna_crosscheck.py
```

## Findings (Kir DMS golden run, 11,729 variants)

- **Method 1: 11,525 / 11,729 confirmed**, incl. **9119/9119 substitutions**.
- **Completeness is full**: every codon position (including boundary positions)
  carries all substitutions + all 3 deletion sizes + all 3 insertion sizes; only
  the final codon lacks del/ins (correct).
- **Method 2 (pydna)** confirmed 165/204 of method-1's boundary blind-spots as
  correctly assembled; **83/83 agreement** with method 1 on assemblable passing
  variants. An apparent "HD instead of DH" boundary insertion was a
  translation-frame artifact of the cross-check harness (disproved by reading the
  assembled DNA in-frame).
- **MAJOR confirmed defect**: **fragment/subpool 7 (codons 293-340) — all 1293
  variants (100%) fail to assemble**; every other fragment is 0%. Flagged by
  DIMPLE's OWN `qc.check_final_assembly` ("DNA cannot be circularized / sticky
  ends not compatible"), not just this harness. Fragment 7's 5' junction overhangs
  are incompatible (oligo 5' `aaat` vs backbone 3' `gaga`, needs `attt`); the 3'
  junction is fine. Scaffold-level (identical across all 1293 variants), likely
  tied to the `Non specific Fragment:1` fragment-size reshuffle during layout.
  Reproduce: run `check_final_assembly` on a pool from `build_pool.py`, or
  `pydna_crosscheck.py`.
- Minor: deletion variants at `wt_start=1` delete the initiator ATG
  (`Kir_delete-1_{3,6,9}-1`) — per the maintainer, acceptable.
- Completeness is full (no missing classes at any boundary); the issue is
  assembly of subpool 7, and the encoding of every variant is correct.
