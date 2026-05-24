"""Regression test: long deletions (5 aa and 10 aa) on the Kir gene.

Guards upstream issue #23: the pipeline crashed in ``find_geneprimer`` with
``IndexError: index out of range`` whenever the inflated ``overlap_r``
(``max(deletions) + overlap - 3``) pushed the gene-primer end position below
the fixed start (``15``). The fix grows ``primer_buffer`` in
``compute_overlaps_and_maxfrag`` so the primer-design window stays at least one
base wide; this test exercises a 5 aa (15 nt) and 10 aa (30 nt) deletion to
cover both the originally-reported sizes.
"""

import shutil

import pytest
from Bio.Seq import Seq

from DIMPLE.DIMPLE import addgene, generate_DMS_fragments
from DIMPLE.pool import DimpleRuntimeConfig
from DIMPLE.run_settings import compute_overlaps_and_maxfrag

_DELETIONS = [15, 30]
_OVERLAP = 4
_OLIGO_LEN = 230


@pytest.mark.slow
def test_long_deletions_kir(tmp_path, dimple_human_usage, kir_fa):
    """5 aa and 10 aa deletion scan on Kir completes without crashing."""
    gene_file = tmp_path / kir_fa.name
    shutil.copy(kir_fa, gene_file)
    wDir = str(tmp_path) + "/"

    config = DimpleRuntimeConfig(
        handle="",
        dms=False,
        stop_codon=True,
        make_double=False,
        maximize_nucleotide_change=False,
        cutsite=Seq("CGTCTC"),
        cutsite_buffer=Seq("G"),
        cutsite_overhang=4,
        enzyme=None,
        avoid_sequence=[Seq("CGTCTC"), Seq("GGTCTC")],
        random_seed=1848,
        usage=dimple_human_usage,
    )
    overlap_l, overlap_r = compute_overlaps_and_maxfrag(
        _OLIGO_LEN, 0, _OVERLAP, _DELETIONS, config=config
    )

    pool = addgene(str(gene_file), config)

    # Should not raise: pre-fix this crashed with IndexError in find_geneprimer.
    generate_DMS_fragments(
        pool,
        overlap_l,
        overlap_r,
        True,  # synonymous
        None,  # custom_mutations
        False,  # dms
        False,  # insert
        _DELETIONS,  # delete
        False,  # dis
        wDir,
    )

    gene = pool[0]

    # At least one D variant exists for each requested deletion size.
    for delete_n in _DELETIONS:
        length = delete_n // 3
        d_variants = [
            v
            for v in gene.designed_variants.values()
            if v["mutation_type"] == "D" and v["length"] == length
        ]
        assert d_variants, f"no D variants emitted for deletion length {delete_n} nt"

    # Every emitted oligo fits within synth_len.
    assert gene.oligos, "no oligos emitted"
    over_length = [(o.id, len(o.seq)) for o in gene.oligos if len(o.seq) > config.synth_len]
    assert not over_length, f"oligos exceed synth_len={config.synth_len}: {over_length[:5]}"
