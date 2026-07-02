"""Build the Kir DMS pool exactly like the regression test, but enzyme=BsmBI so
check_final_assembly's PCR/cut/ligate machinery is available for cross-checking."""

import sys, shutil, tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from Bio.Seq import Seq
from DIMPLE.DIMPLE import addgene, generate_DMS_fragments
from DIMPLE.pool import DimpleRuntimeConfig
from tests.conftest import _HUMAN_USAGE

_OVERLAP = 4


def build():
    tmp = Path(tempfile.mkdtemp())
    gene_file = tmp / "Kir.fa"
    shutil.copy(REPO / "tests/expected/Kir.fa", gene_file)
    wDir = str(tmp) + "/"
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
        enzyme="BsmBI",
        avoid_sequence=[Seq("CGTCTC"), Seq("GGTCTC")],
        random_seed=1848,
        usage=dict(_HUMAN_USAGE),
    )
    pool = addgene(str(gene_file), config)
    generate_DMS_fragments(
        pool,
        _OVERLAP,
        _OVERLAP,
        True,
        None,
        True,
        ["GAC", "GACCAT", "GACCATGTA"],
        [3, 6, 9],
        False,
        wDir,
    )
    return pool


if __name__ == "__main__":
    pool = build()
    g = pool[0]
    print("gene", g.geneid, "n genePrimer", len(g.genePrimer), "n barPrimer", len(g.barPrimer))
    dv = g.designed_variants
    print("n designed_variants", len(dv))
    keys = list(dv)[:6]
    print("sample keys:", keys)
    k0 = keys[0]
    print("entry fields:", list(dv[k0].keys()))
    print(
        "entry:",
        {
            f: (
                str(v.seq)[:36]
                if hasattr(v, "seq")
                else (str(v)[:36] if not isinstance(v, (int, float)) else v)
            )
            for f, v in dv[k0].items()
        },
    )
