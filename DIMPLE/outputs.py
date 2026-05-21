"""Output aggregation helpers for the DIMPLE pipeline.

Writes the combined oligo and primer FASTA files, migrated out of
``DIMPLE.DIMPLE``.
"""

from __future__ import annotations

import os

from Bio import SeqIO

from DIMPLE.core import DIMPLE


def print_all(pool, folder="", config=None):
    """Writes oligos and primers to files."""
    if not isinstance(pool[0], DIMPLE):
        raise TypeError("Not an instance of the DIMPLE class")
    alloligos = []
    allprimers = []
    for obj in pool:
        try:
            alloligos.extend(obj.oligos)
            allprimers.extend(obj.barPrimer)
            allprimers.extend(obj.genePrimer)
        except AttributeError:
            print(obj.geneid + " has not been processed")
    # Remove redundant sequences?
    SeqIO.write(
        alloligos, os.path.join(folder.replace("\\", ""), "All_Oligos.fasta"), "fasta"
    )
    SeqIO.write(
        allprimers, os.path.join(folder.replace("\\", ""), "All_Primers.fasta"), "fasta"
    )
