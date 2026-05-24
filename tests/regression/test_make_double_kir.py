"""Regression test: DMS + ``make_double=True`` on the Kir gene.

Guards two bugs that together made ``make_double=True`` raise on any single-gene
run for ~12 months:

1. **Off-by-one between ``positions`` and mutation-name AA numbers.** The DMS
   path built ``positions`` with ``+ 3`` but mutation names with ``+ 6``, so
   ``positions.index(int(re.findall(r'\\d+', combi[N])[0]))`` raised
   ``ValueError: <N> is not in list`` for the last AA of every fragment.

2. **Double-mutant SeqRecords never registered in ``designed_variants``.** The
   shared barcode-assembly loop writes ``designed_variants[sequence.id][
   "oligo_sequence"] = ...`` for every emitted sequence, so the missing
   registration crashed with ``KeyError: 'Kir_DMS-N_AAA1+BBB2'``.

Uses a 2-AA amino-acid set (``Ala``, ``Cys``) to keep the pair-count tractable
(full 20-AA scan generates ~480k double-mutant pairs per fragment, which OOMs
on small CI workers). The bugs reproduce identically at this scale.
"""

import shutil

import pytest
from Bio.Seq import Seq

from DIMPLE.DIMPLE import addgene, generate_DMS_fragments
from DIMPLE.pool import DimpleRuntimeConfig
from DIMPLE.run_settings import apply_instance_settings

_OVERLAP = 3


@pytest.mark.slow
def test_make_double_kir(tmp_path, dimple_human_usage, kir_fa):
    gene_file = tmp_path / kir_fa.name
    shutil.copy(kir_fa, gene_file)
    wDir = str(tmp_path) + "/"

    config = DimpleRuntimeConfig(
        handle="",
        synth_len=230,
        maxfrag=230 - 62 - _OVERLAP,
        primer_buffer=30 + _OVERLAP,
        dms=True,
        stop_codon=False,
        make_double=True,
        maximize_nucleotide_change=False,
        cutsite=Seq("CGTCTC"),
        cutsite_buffer=Seq("G"),
        cutsite_overhang=4,
        avoid_sequence=[Seq("CGTCTC"), Seq("GGTCTC")],
        random_seed=1848,
        usage=dimple_human_usage,
    )

    pool = addgene(str(gene_file), config)
    # Tractable scan: 2 amino acids → ~98 singles/fragment → ~5k double-mutant
    # pairs/fragment instead of ~480k.
    apply_instance_settings(pool, config=config, aminoacids=["Ala", "Cys"])

    generate_DMS_fragments(
        pool,
        _OVERLAP,
        _OVERLAP,
        False,  # synonymous
        None,  # custom_mutations
        True,  # dms
        False,  # insert
        False,  # delete
        False,  # dis
        wDir,
    )

    gene = pool[0]

    # Bug 2 guard: doubles registered in designed_variants.
    doubles = [v for v in gene.designed_variants.values() if v["mutation_type"] == "MM"]
    assert doubles, "no double-mutant variants registered"

    # Bug 1 guard: every emitted oligo has a registered variant carrying its sequence.
    # (Without the off-by-one fix the run would crash before reaching this assertion.)
    for oligo in gene.oligos:
        assert oligo.id in gene.designed_variants, oligo.id
        assert gene.designed_variants[oligo.id]["oligo_sequence"] == oligo.seq, oligo.id

    # A double-mutant variant looks structurally right.
    sample = doubles[0]
    assert sample["length"] == 2
    assert "+" in sample["name"]
    assert sample["hgvs"].startswith("p.([") and sample["hgvs"].endswith("])")
