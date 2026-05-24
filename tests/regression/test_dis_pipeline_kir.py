"""Regression test: domain-insertion scan (DIS) on the Kir gene.

Guards the bug where the ``dis`` branch of ``generate_DMS_fragments`` built
oligo records but never registered them in ``gene.designed_variants``. The
later barcode-assembly loop then did::

    gene.designed_variants[sequence.id]['oligo_sequence'] = combined_sequence

which raised ``KeyError: 'Kir_DIS-1_1'`` for any run including DIS.

This test runs the pipeline with ``dis=True`` and asserts it completes and
that every emitted oligo has a matching ``designed_variants`` entry carrying
an ``oligo_sequence``.
"""

import shutil

import pytest
from Bio.Seq import Seq

from DIMPLE.DIMPLE import addgene, generate_DMS_fragments
from DIMPLE.pool import DimpleRuntimeConfig

# Default domain-insertion handle from run_dimple.py (BsaI-based, 24 nt).
_HANDLE = "AGCGGGAGACCGGGGTCTCTGAGC"
_OVERLAP = 3


@pytest.mark.slow
def test_dis_pipeline_kir(tmp_path, dimple_human_usage, kir_fa):
    """DIS scan on Kir registers a designed_variants entry per oligo."""
    gene_file = tmp_path / kir_fa.name
    shutil.copy(kir_fa, gene_file)
    wDir = str(tmp_path) + "/"

    # Run configuration mirrors test_dms_pipeline_kir, with dis-specific values.
    config = DimpleRuntimeConfig(
        handle=_HANDLE,
        synth_len=230,
        maxfrag=230 - 62 - _OVERLAP,
        primer_buffer=30 + _OVERLAP,
        dms=False,
        stop_codon=True,
        make_double=False,
        maximize_nucleotide_change=False,
        cutsite=Seq("CGTCTC"),
        cutsite_buffer=Seq("G"),
        cutsite_overhang=4,
        enzyme=None,
        # The handle deliberately carries a BsaI site; only avoid the BsmBI
        # cassette site so the DIS restriction-site guard does not fire.
        avoid_sequence=[Seq("CGTCTC")],
        random_seed=1848,
        usage=dimple_human_usage,
    )

    pool = addgene(str(gene_file), config)

    generate_DMS_fragments(
        pool,
        _OVERLAP,
        _OVERLAP,
        True,    # synonymous
        None,    # custom_mutations
        False,   # dms
        False,   # insert
        False,   # delete
        True,    # dis
        wDir,
    )

    gene = pool[0]
    # The crash being guarded: designed_variants is populated for DIS.
    dis_variants = {
        k: v for k, v in gene.designed_variants.items() if "_DIS-" in k
    }
    assert dis_variants, "DIS run produced no designed_variants entries"

    # Every emitted oligo must have a registered variant carrying its sequence.
    for oligo in gene.oligos:
        assert oligo.id in gene.designed_variants, oligo.id
        variant = gene.designed_variants[oligo.id]
        assert variant["mutation_type"] == "DI"
        assert variant["oligo_sequence"] == oligo.seq
