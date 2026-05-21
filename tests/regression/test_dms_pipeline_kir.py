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

import pytest
from Bio.Seq import Seq

from DIMPLE.DIMPLE import (
    DIMPLE,
    addgene,
    generate_DMS_fragments,
    post_qc,
    print_all,
)

from pathlib import Path

EXPECTED = Path(__file__).parent.parent / "expected"

_GOLDEN_FILES = [
    "Kir_DMS_Gene_Primers.fasta",
    "Kir_DMS_Oligo_Primers.fasta",
    "Kir_DMS_Oligos.fasta",
    "Kir_mutations.csv",
    "Kir_designed_variants.csv",
]


@pytest.mark.slow
def test_dms_pipeline_kir(tmp_path, dimple_human_usage, kir_fa, update_golden):
    """Full DMS pipeline on the Kir gene; byte-compares outputs to golden files."""
    gene_file = tmp_path / kir_fa.name
    shutil.copy(kir_fa, gene_file)

    # Mirror the parameter setup from the original test / run_dimple.py.
    # dimple_state (via dimple_human_usage) snapshots and restores all of these.
    DIMPLE.handle = ""
    DIMPLE.overlap = 3
    DIMPLE.synth_len = 230          # shadows the property descriptor with a plain int
    DIMPLE.maxfrag = DIMPLE.synth_len - 62 - DIMPLE.overlap  # 165
    DIMPLE.primerBuffer += DIMPLE.overlap                     # 30 + 3 = 33

    DIMPLE.dms = True
    DIMPLE.stop_codon = True
    DIMPLE.make_double = False
    DIMPLE.maximize_nucleotide_change = False

    restriction_sequence = "CGTCTC(G)1/5"
    tmp_cutsite = restriction_sequence.split("(")
    DIMPLE.cutsite = Seq(tmp_cutsite[0])
    DIMPLE.cutsite_buffer = Seq(tmp_cutsite[1].split(")")[0])
    tmp_overhang = tmp_cutsite[1].split(")")[1].split("/")
    DIMPLE.cutsite_overhang = int(tmp_overhang[1]) - int(tmp_overhang[0])
    DIMPLE.enzyme = None  # skip golden-gate assembly QC; doesn't affect output files

    DIMPLE.avoid_sequence = [Seq(x) for x in ["CGTCTC", "GGTCTC"]]
    DIMPLE.random_seed = 1848

    wDir = str(tmp_path) + "/"
    pool = addgene(str(gene_file))

    generate_DMS_fragments(
        pool,
        DIMPLE.overlap,
        DIMPLE.overlap,
        True,                                   # synonymous
        None,                                   # custom_mutations
        DIMPLE.dms,
        ["GAC", "GACCAT", "GACCATGTA"],         # insert
        [3, 6, 9],                              # delete
        False,                                  # dis
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
