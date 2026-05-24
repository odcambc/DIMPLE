"""Regression test: combined DMS + DIS scan on the Kir gene.

Guards the bug where ``generate_DMS_fragments`` overflowed ``synth_len`` when
both DMS and DIS were enabled. ``DMS`` and ``DIS`` push fragments into the
same ``dms_sequences`` list per gene-fragment, but DIS fragments are
``len(handle)`` longer than DMS substitutions. The barcode pair sized once
from the smallest fragment got concatenated whole onto the longer DIS
sequences, producing oligos `handle_len - 2*overhang` bytes over ``synth_len``::

    Exception: Oligo too long: 246 is longer than 230

Fixed by extending the per-sequence trim branch in the oligo-assembly loop
to also fire when ``dis=True``. This test runs DMS+DIS and asserts the run
completes, every oligo is at most ``synth_len`` long, and the
``designed_variants`` table includes both mutation types.
"""

import shutil

import pytest
from Bio.Seq import Seq

from DIMPLE.DIMPLE import addgene, generate_DMS_fragments
from DIMPLE.pool import DimpleRuntimeConfig

_HANDLE = "AGCGGGAGACCGGGGTCTCTGAGC"
_OVERLAP = 3


@pytest.mark.slow
def test_dms_dis_combined_kir(tmp_path, dimple_human_usage, kir_fa):
    """DMS+DIS together stays within synth_len for every emitted oligo."""
    gene_file = tmp_path / kir_fa.name
    shutil.copy(kir_fa, gene_file)
    wDir = str(tmp_path) + "/"

    config = DimpleRuntimeConfig(
        handle=_HANDLE,
        synth_len=230,
        maxfrag=230 - 62 - _OVERLAP,
        primer_buffer=30 + _OVERLAP,
        dms=True,
        stop_codon=False,
        make_double=False,
        maximize_nucleotide_change=False,
        cutsite=Seq("CGTCTC"),
        cutsite_buffer=Seq("G"),
        cutsite_overhang=4,
        enzyme=None,
        # The handle carries a BsaI site; only avoid the BsmBI cassette site.
        avoid_sequence=[Seq("CGTCTC")],
        random_seed=1848,
        usage=dimple_human_usage,
    )

    pool = addgene(str(gene_file), config)

    generate_DMS_fragments(
        pool,
        _OVERLAP,
        _OVERLAP,
        False,  # synonymous
        None,  # custom_mutations
        True,  # dms
        False,  # insert
        False,  # delete
        True,  # dis
        wDir,
    )

    gene = pool[0]

    # The crash being guarded: every emitted oligo must fit in synth_len.
    overlong = [o for o in gene.oligos if len(o.seq) > config.synth_len]
    assert not overlong, (
        f"{len(overlong)} oligos exceed synth_len={config.synth_len}: "
        f"max={max(len(o.seq) for o in overlong)}"
    )

    # Both mutation types should be represented in designed_variants.
    types = {v["mutation_type"] for v in gene.designed_variants.values()}
    assert "M" in types or "S" in types, f"no DMS substitution variants registered (types: {types})"
    assert "DI" in types, f"no DI variants registered (types: {types})"
