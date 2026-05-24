"""Regression test: full DMS pipeline on the Kir gene.

Runs generate_DMS_fragments → post_qc → print_all with the canonical Kir FASTA
and compares five output files byte-for-byte against golden files in
tests/expected/. The five golden comparisons are:

    Kir_DMS_Gene_Primers.fasta
    Kir_DMS_Oligo_Primers.fasta
    Kir_DMS_Oligos.fasta
    Kir_mutations.csv
    Kir_designed_variants.csv

Pass --update-golden to regenerate the golden files instead of asserting.
"""

import shutil
from pathlib import Path

import pytest
from Bio.Seq import Seq

from DIMPLE.DIMPLE import addgene, generate_DMS_fragments, post_qc, print_all
from DIMPLE.pool import DimpleRuntimeConfig

EXPECTED = Path(__file__).parent.parent / "expected"

_GOLDEN_FILES = [
    "Kir_DMS_Gene_Primers.fasta",
    "Kir_DMS_Oligo_Primers.fasta",
    "Kir_DMS_Oligos.fasta",
    "Kir_mutations.csv",
    "Kir_designed_variants.csv",
]

_OVERLAP = 3


@pytest.mark.slow
def test_dms_pipeline_kir(tmp_path, dimple_human_usage, kir_fa, update_golden):
    """Full DMS pipeline on the Kir gene; byte-compares outputs to golden files."""
    gene_file = tmp_path / kir_fa.name
    shutil.copy(kir_fa, gene_file)
    wDir = str(tmp_path) + "/"

    # Run configuration -- mirrors run_dimple.py for the Kir DMS scan.
    # Restriction "CGTCTC(G)1/5": cutsite CGTCTC, buffer G, overhang 5 - 1 = 4.
    config = DimpleRuntimeConfig(
        handle="",
        synth_len=230,
        maxfrag=230 - 62 - _OVERLAP,  # 165
        primer_buffer=30 + _OVERLAP,  # PRIMER_BUFFER_BASE + overlap = 33
        dms=True,
        stop_codon=True,
        make_double=False,
        maximize_nucleotide_change=False,
        cutsite=Seq("CGTCTC"),
        cutsite_buffer=Seq("G"),
        cutsite_overhang=4,
        enzyme=None,  # skip golden-gate assembly QC
        avoid_sequence=[Seq("CGTCTC"), Seq("GGTCTC")],
        random_seed=1848,
        usage=dimple_human_usage,
    )

    pool = addgene(str(gene_file), config)

    generate_DMS_fragments(
        pool,
        _OVERLAP,
        _OVERLAP,
        True,  # synonymous
        None,  # custom_mutations
        True,  # dms
        ["GAC", "GACCAT", "GACCATGTA"],  # insert
        [3, 6, 9],  # delete
        False,  # dis
        wDir,
    )
    post_qc(pool)
    print_all(pool, wDir)

    if update_golden:
        for name in _GOLDEN_FILES:
            shutil.copy(tmp_path / name, EXPECTED / name)
        return

    for name in _GOLDEN_FILES:
        actual = (tmp_path / name).read_text()
        expected = (EXPECTED / name).read_text()
        assert actual == expected, f"Output mismatch in {name}"
