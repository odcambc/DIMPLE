"""The oligo ``Pool`` -- one DIMPLE run's genes plus their shared run config.

A DIMPLE run takes a batch of genes and produces a single **oligo pool** (the
synthesis-vendor deliverable; the wet-lab variant library is generated
downstream from the pool). The pipeline previously passed this around as a bare
list of genes; ``Pool`` gives it a class and a home for the run configuration.

This module owns three cohesive things: ``DimpleRuntimeConfig`` (the typed run
config), ``Pool``, and ``addgene`` (the Pool factory). It is a near-leaf in the
import graph -- it imports only ``DIMPLE.core`` -- so ``run_settings`` and the
pipeline modules can depend on it without a cycle.
"""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass, field
from typing import List, Optional

from Bio import SeqIO
from Bio.Seq import Seq

from DIMPLE.core import DIMPLE

logger = logging.getLogger(__name__)

# Base primer buffer before overlap extension (matches DIMPLE.primerBuffer default).
PRIMER_BUFFER_BASE: int = 30

# Bundled barcode primer sets, loaded once. Each DimpleRuntimeConfig gets its
# own copy because barcodes are consumed (popped) as a run assigns subpools.
_DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
try:
    _BARCODES_F = list(SeqIO.parse(os.path.join(_DATA_DIR, "forward_finalprimers.fasta"), "fasta"))
    _BARCODES_R = list(SeqIO.parse(os.path.join(_DATA_DIR, "reverse_finalprimers.fasta"), "fasta"))
except FileNotFoundError as exc:
    raise ValueError(
        "Could not find barcode files. Please upload your own or place standard "
        "barcodes in the data file."
    ) from exc


@dataclass
class DimpleRuntimeConfig:
    """Explicit runtime settings used by CLI/GUI and core pipeline."""

    handle: str = ""
    usage: Optional[dict] = None
    cutsite: Optional[Seq] = None
    cutsite_buffer: Optional[Seq] = None
    cutsite_overhang: int = 0
    enzyme: Optional[str] = None
    synth_len: Optional[int] = None
    maxfrag: Optional[int] = None
    primer_buffer: int = PRIMER_BUFFER_BASE
    random_seed: Optional[int] = 0
    avoid_sequence: Optional[List[Seq]] = None
    barcode_f: list = field(default_factory=lambda: list(_BARCODES_F))
    barcode_r: list = field(default_factory=lambda: list(_BARCODES_R))
    dms: bool = False
    stop_codon: bool = False
    make_double: bool = False
    maximize_nucleotide_change: bool = False
    gene_primer_tm: tuple[int, int] = (58, 62)
    non_interactive: bool = False
    preferred_orf_index: Optional[int] = None
    link_policy: str = "prompt"
    breaksite_change_policy: str = "prompt"


class Pool:
    """An oligo pool: the ``DIMPLE`` genes of one run plus their shared config.

    Iterating / indexing a ``Pool`` yields its gene instances, so it substitutes
    transparently for the bare gene list the pipeline previously passed around
    while config reads are migrated onto ``pool.config``.
    """

    def __init__(self, config: DimpleRuntimeConfig):
        self.config = config
        self.genes: list = []

    def append(self, gene) -> None:
        self.genes.append(gene)

    def extend(self, genes) -> None:
        self.genes.extend(genes)

    def __iter__(self):
        return iter(self.genes)

    def __getitem__(self, index):
        return self.genes[index]

    def __len__(self) -> int:
        return len(self.genes)


def addgene(genefile, config: Optional[DimpleRuntimeConfig] = None, start=None, end=None):
    """Build an oligo :class:`Pool` from a FASTA file containing one or more genes.

    Each gene becomes a :class:`DIMPLE` instance carrying a ``pool``
    back-reference. ``config`` is transitional-optional: when omitted it falls
    back to the ``run_settings`` runtime-config singleton (Phase 3 of the Pool
    migration makes it required).
    """
    if config is None:
        # Transitional fallback until all callers pass an explicit config.
        from DIMPLE.run_settings import get_runtime_config

        config = get_runtime_config()
    pool = Pool(config)
    print("Barcode: " + str(config.barcode_f[0].seq))
    print("Number of barcodes: " + str(len(config.barcode_f)))
    if start is None:
        start = []
    if end is None:
        end = []
    tmpgene = list(SeqIO.parse(genefile.replace("\\", ""), "fasta"))
    tmpgene[0].seq = tmpgene[0].seq.upper()
    for gene in tmpgene:
        if "start:" in gene.description and "end:" in gene.description:
            start = int(gene.description.split("start:")[1].split(" ")[0]) - 1
            end = int(gene.description.split("end:")[1].split(" ")[0])
            gene.filename = genefile.replace("\\", "")
            logger.info("Found start: " + str(start) + " and end: " + str(end))
            logger.info("Inferred ORF sequence: " + str(gene.seq[start:end]))
            logger.info("ORF translation: " + str(gene.seq[start:end].translate()))
            instance = DIMPLE(gene, start, end, pool)
        else:
            gene.filename = genefile.replace("\\", "")
            instance = DIMPLE(gene, start, end, pool)
        pool.append(instance)
    return pool
