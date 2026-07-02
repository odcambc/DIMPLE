"""Regression test: Golden-Gate assembly correctness on the Kir gene.

This guards the class of bug where a fragment's Type IIS overhang (or an
edit-induced coding site) reconstitutes the enzyme recognition site, so oligos
gain a spurious internal cut and cannot assemble. It caught fragment 7 of the
Kir DMS run (1293/1293 variants unassemblable) -- invisible to the byte-compare
regressions because those run with ``enzyme=None``, which disables the assembly
QC entirely.

Unlike the byte-compare tests, this one runs with ``enzyme="BsmBI"`` and asserts
that every SHIPPED oligo assembles cleanly and that only a small handful of
variants are dropped for unavoidable edit-induced internal sites.
"""

import logging
import shutil

import pytest
from Bio.Seq import Seq

from DIMPLE.DIMPLE import addgene, check_final_assembly, generate_DMS_fragments, post_qc
from DIMPLE.pool import DimpleRuntimeConfig

_OVERLAP = 4
# A few variants inherently carry the enzyme site (fixed insert/delete handle, or
# no synonymous codon avoids it) and are legitimately dropped. A whole subpool
# failing (the fragment-7 bug) would blow far past this.
_MAX_EXPECTED_DROPS = 25


@pytest.mark.slow
def test_kir_oligos_assemble(tmp_path, dimple_human_usage, kir_fa, caplog):
    gene_file = tmp_path / kir_fa.name
    shutil.copy(kir_fa, gene_file)
    wDir = str(tmp_path) + "/"

    config = DimpleRuntimeConfig(
        handle="",
        synth_len=230,
        maxfrag=230 - 64 - _OVERLAP - _OVERLAP,
        primer_buffer=30 + _OVERLAP,
        dms=True,
        stop_codon=True,
        make_double=False,
        maximize_nucleotide_change=False,
        cutsite=Seq("CGTCTC"),
        cutsite_buffer=Seq("G"),
        cutsite_overhang=4,
        enzyme="BsmBI",  # ENABLE the assembly QC that the byte-compare tests skip
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

    gene = pool[0]

    # Only a small handful of variants should be dropped for internal sites; a
    # large or subpool-concentrated drop means a boundary overhang reconstituted
    # the cut site (the fragment-7 class) and must fail here.
    assert len(gene.dropped_oligos) <= _MAX_EXPECTED_DROPS, (
        f"{len(gene.dropped_oligos)} oligos dropped for internal restriction sites "
        f"(expected <= {_MAX_EXPECTED_DROPS}): {gene.dropped_oligos[:20]}"
    )

    # Every oligo that ships must assemble. check_final_assembly logs an error per
    # unassemblable variant; assert none fire.
    with caplog.at_level(logging.ERROR, logger="DIMPLE.qc"):
        check_final_assembly(gene)
    unassemblable = [
        r.getMessage()
        for r in caplog.records
        if "does not assemble with template" in r.getMessage()
    ]
    assert (
        not unassemblable
    ), f"{len(unassemblable)} shipped oligos fail to assemble: {unassemblable[:10]}"

    # post_qc must not hard-fail on this handful of drops (stays under the fraction
    # limit) and its report exercises the drop-summary path.
    post_qc(pool)
