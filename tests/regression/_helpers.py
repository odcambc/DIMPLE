"""Shared post-``print_all`` invariant checks for regression tests.

Asserts: no duplicate oligo IDs, no duplicate variant rows, row-count parity
between the oligo FASTA and the variants CSV, every oligo within
``config.synth_len``, every gene primer of plausible length.
"""

from __future__ import annotations

import csv
from pathlib import Path

from Bio import SeqIO

from DIMPLE.pool import DimpleRuntimeConfig

# find_geneprimer's design window is ~20-25 nt with a ±12 nt adjust loop;
# anything past 60 nt means the loop ran off the rails.
_GENE_PRIMER_MAX_LEN = 60


def assert_outputs_consistent(
    work_dir,
    gene_id: str,
    config: DimpleRuntimeConfig,
) -> None:
    """Assert invariants on per-gene DIMPLE output files in *work_dir*."""
    work = Path(work_dir)
    oligos_path = work / f"{gene_id}_DMS_Oligos.fasta"
    variants_path = work / f"{gene_id}_designed_variants.csv"
    primers_path = work / f"{gene_id}_DMS_Gene_Primers.fasta"

    for p in (oligos_path, variants_path, primers_path):
        assert p.is_file(), f"missing expected output: {p}"

    oligo_records = list(SeqIO.parse(oligos_path, "fasta"))
    oligo_ids = [r.id for r in oligo_records]

    dup_oligo_ids = _duplicates(oligo_ids)
    assert not dup_oligo_ids, (
        f"{len(dup_oligo_ids)} duplicate oligo IDs in {oligos_path.name}: "
        f"{sorted(dup_oligo_ids)[:3]}"
    )

    overlong = [(r.id, len(r.seq)) for r in oligo_records if len(r.seq) > config.synth_len]
    assert not overlong, (
        f"{len(overlong)} oligos exceed synth_len={config.synth_len}; "
        f"first few: {overlong[:3]}"
    )

    with variants_path.open(newline="") as f:
        rows = list(csv.DictReader(f))
    variant_names = [row["name"] for row in rows]

    dup_variants = _duplicates(variant_names)
    assert not dup_variants, (
        f"{len(dup_variants)} duplicate variant names in {variants_path.name}: "
        f"{sorted(dup_variants)[:3]}"
    )

    assert len(oligo_records) == len(rows), (
        f"oligo/variant count mismatch: {len(oligo_records)} oligos in "
        f"{oligos_path.name} vs {len(rows)} rows in {variants_path.name}"
    )

    primer_records = list(SeqIO.parse(primers_path, "fasta"))
    bad_primers = [
        (r.id, len(r.seq))
        for r in primer_records
        if len(r.seq) == 0 or len(r.seq) > _GENE_PRIMER_MAX_LEN
    ]
    assert not bad_primers, (
        f"{len(bad_primers)} gene primers have implausible length "
        f"(empty or > {_GENE_PRIMER_MAX_LEN} nt): {bad_primers[:3]}"
    )


def _duplicates(items):
    seen: set = set()
    dups: set = set()
    for item in items:
        if item in seen:
            dups.add(item)
        seen.add(item)
    return dups
