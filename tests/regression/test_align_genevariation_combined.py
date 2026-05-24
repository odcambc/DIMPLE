"""Regression test for ``align_genevariation``.

``align_genevariation`` is opt-in (via ``-matchSequences match`` in the CLI; the
GUI checkbox is commented out at ``run_dimple_gui.py:502-503``) and has no
default test coverage. A 2024 commit (``3a7503bf``) silently broke its
breaksite remap by wrapping a list-comprehension in parens — fixed
2026-05-23 on the architecture branch, but until coverage exists future
regressions are easy.

Uses ``tests/data/Kir_pair.fa`` — a synthetic homologous pair built by
copying ``Kir.fa`` with 5 single-nt substitutions in non-ORF regions.
Kir vs. Shaker (the natural multi-gene fixture in ``combined_fasta.fa``)
scores below the ``> 1.5`` alignment threshold so doesn't exercise the
linking branch.

Asserts the key invariants:

  * Both genes end up in each other's ``gene.linked`` set.
  * Linked genes share a common ``breaklist`` (so they get the same fragment
    boundaries — the whole point of the function).
  * ``unique_Frag`` is populated per-gene with one boolean per fragment,
    flagging True where the fragment differs between linked genes.
"""

import shutil
from pathlib import Path

import pytest
from Bio.Seq import Seq

from DIMPLE.DIMPLE import addgene, align_genevariation
from DIMPLE.pool import DimpleRuntimeConfig

_OVERLAP = 3
_FIXTURE = Path(__file__).parent.parent / "data" / "Kir_pair.fa"


@pytest.mark.slow
def test_align_genevariation_links_homologs(tmp_path, dimple_human_usage):
    """Two near-identical genes get linked, share breaklist, populate unique_Frag."""
    assert _FIXTURE.exists(), f"missing fixture: {_FIXTURE}"
    gene_file = tmp_path / _FIXTURE.name
    shutil.copy(_FIXTURE, gene_file)

    config = DimpleRuntimeConfig(
        handle="",
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
        avoid_sequence=[Seq("CGTCTC"), Seq("GGTCTC")],
        random_seed=1848,
        usage=dimple_human_usage,
        # Force linking without prompting; matches what -matchSequences match passes.
        link_policy="always",
        non_interactive=True,
    )

    pool = addgene(str(gene_file), config)
    assert len(pool) == 2, "fixture should yield two genes"

    align_genevariation(pool)

    kir, homolog = pool[0], pool[1]

    # Both must be cross-linked.
    assert hasattr(kir, "linked"), "Kir.linked never assigned by align_genevariation"
    assert hasattr(homolog, "linked"), "homolog.linked never assigned"
    assert 1 in kir.linked, f"Kir not linked to homolog; kir.linked={kir.linked}"
    assert 0 in homolog.linked, f"homolog not linked to Kir; homolog.linked={homolog.linked}"

    # Shared fragment boundaries.
    assert (
        kir.breaklist == homolog.breaklist
    ), "linked genes should share breaklist after align_genevariation"

    # unique_Frag populated and well-shaped per gene.
    assert hasattr(kir, "unique_Frag") and hasattr(homolog, "unique_Frag")
    assert len(kir.unique_Frag) == len(
        kir.breaklist
    ), f"unique_Frag length {len(kir.unique_Frag)} != breaklist length {len(kir.breaklist)}"
    assert len(homolog.unique_Frag) == len(homolog.breaklist)
    # Each entry is bool-like (True or False).
    assert all(isinstance(x, bool) for x in kir.unique_Frag)
    # The two synthetic mutations live within the ORF region the pipeline
    # mutates, so we expect at least one fragment-pair to differ.
    differing = [a != b for a, b in zip(kir.unique_Frag, homolog.unique_Frag)]
    assert any(kir.unique_Frag) or any(differing), (
        f"expected at least one unique fragment between Kir and its homolog; "
        f"kir.unique_Frag={kir.unique_Frag}, homolog.unique_Frag={homolog.unique_Frag}"
    )
