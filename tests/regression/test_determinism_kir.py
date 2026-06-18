"""Determinism contract: same config + seed → byte-identical outputs.

Runs the full DMS pipeline on Kir twice with an identical configuration and a
pinned ``random_seed``, into two separate work directories, then byte-compares
every output file. Asserts the two runs are identical.

This is independent of the golden fixtures: the goldens prove "this config
produces *these* bytes"; this test proves "this config produces *the same* bytes
every time", which would still catch a newly-introduced nondeterministic source
(e.g. an unseeded stdlib ``random.choice`` slipping into codon selection) even
if someone regenerated the goldens to match the drift.

Each run builds a fresh config + ``addgene`` -- the config carries ``Seq``
cutsites that are mutated in place during fragment generation and ``addgene``
builds stateful per-gene objects, so reusing either across runs would not be a
clean-start comparison.
"""

import shutil
from pathlib import Path

import pytest
from Bio.Seq import Seq

from DIMPLE.DIMPLE import addgene, generate_DMS_fragments, post_qc, print_all
from DIMPLE.pool import DimpleRuntimeConfig

_OVERLAP = 3
_SEED = 1848


def _run(work_dir: Path, gene_file: Path, usage: dict) -> None:
    """Run the full DMS pipeline into *work_dir* with a pinned seed."""
    local_gene = work_dir / gene_file.name
    shutil.copy(gene_file, local_gene)
    wDir = str(work_dir) + "/"

    config = DimpleRuntimeConfig(
        handle="",
        synth_len=230,
        maxfrag=230 - 62 - _OVERLAP,
        primer_buffer=30 + _OVERLAP,
        dms=True,
        stop_codon=True,
        make_double=False,
        maximize_nucleotide_change=False,
        cutsite=Seq("CGTCTC"),
        cutsite_buffer=Seq("G"),
        cutsite_overhang=4,
        enzyme=None,
        avoid_sequence=[Seq("CGTCTC"), Seq("GGTCTC")],
        random_seed=_SEED,
        usage=dict(usage),
    )

    pool = addgene(str(local_gene), config)
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


def _output_files(work_dir: Path) -> dict:
    """Map output filename -> bytes for every .fasta/.csv the run produced."""
    files = {}
    for path in sorted(work_dir.iterdir()):
        if path.suffix in (".fasta", ".csv"):
            files[path.name] = path.read_bytes()
    return files


@pytest.mark.slow
def test_dms_run_is_deterministic(tmp_path, dimple_human_usage, kir_fa):
    """Two identical seeded runs produce byte-identical output files."""
    work_a = tmp_path / "run_a"
    work_b = tmp_path / "run_b"
    work_a.mkdir()
    work_b.mkdir()

    _run(work_a, kir_fa, dimple_human_usage)
    _run(work_b, kir_fa, dimple_human_usage)

    out_a = _output_files(work_a)
    out_b = _output_files(work_b)

    # Guard against a vacuous pass if outputs silently stopped being written.
    assert out_a, "first run produced no .fasta/.csv outputs"
    assert out_a.keys() == out_b.keys(), (
        f"runs produced different file sets: "
        f"only in A={sorted(set(out_a) - set(out_b))}, "
        f"only in B={sorted(set(out_b) - set(out_a))}"
    )

    mismatched = [name for name in out_a if out_a[name] != out_b[name]]
    assert not mismatched, f"nondeterministic output in: {mismatched}"
