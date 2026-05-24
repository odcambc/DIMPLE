"""Core domain model for DIMPLE: the per-gene ``DIMPLE`` class.

A ``DIMPLE`` instance models one gene of an oligo pool -- its sequence,
fragments, oligos, and primers. Run-wide configuration lives on the pool's
:class:`~DIMPLE.pool.DimpleRuntimeConfig`, reachable from any gene via
``gene.pool.config``. This module is a leaf in the package import graph -- it
depends only on ``DIMPLE.utilities`` and biopython -- so ``pool``, the
functional modules, and the ``DIMPLE.DIMPLE`` orchestrator can import the class
from here without an import cycle.
"""

from __future__ import annotations

import logging

import numpy as np

from DIMPLE.utilities import findORF

logger = logging.getLogger(__name__)


class DIMPLE:
    """A single gene of a DIMPLE oligo pool."""

    # Fixed constants (not run configuration -- never change per run).
    maxfrag_offset = 64  # nt reserved for barcodes, cut sites, and handle
    minfrag = 24  # smallest size for golden gate fragment efficiency
    primerTm = (56.5, 60)  # melting temperature limits for fragment primers

    def __init__(self, gene, start=None, end=None, pool=None):

        self.pool = pool
        cfg = pool.config

        # Set up random number generator
        self.rng = np.random.default_rng(cfg.random_seed)

        #  Search for ORF
        maxfrag = cfg.maxfrag
        if maxfrag is None:
            maxfrag = (
                cfg.synth_len - self.maxfrag_offset
            )  # based on space for barcodes, cut sites, handle. Doesn't need to be exact
            cfg.maxfrag = maxfrag
        # Per-gene maxfrag: seeded from config, but switch_fragmentsize may
        # decrement it per gene while resolving fragment-boundary conflicts.
        self.maxfrag = maxfrag

        self.geneid = gene.name
        self.linked = set()
        self.genePrimer = []
        self.oligos = []
        self.barPrimer = []
        self.fullGene = gene.seq.upper()
        self.split = 0
        self.num_frag_per_oligo = 1
        self.doublefrag = 0
        self.filename = gene.filename
        # Set up variables. Could have this as user input in the class
        self.SynonymousCodons = {
            "Cys": ["TGT", "TGC"],
            "Asp": ["GAT", "GAC"],
            "Ser": ["TCT", "TCG", "TCA", "TCC", "AGC", "AGT"],
            "Gln": ["CAA", "CAG"],
            "Met": ["ATG"],
            "Asn": ["AAC", "AAT"],
            "Pro": ["CCT", "CCG", "CCA", "CCC"],
            "Lys": ["AAG", "AAA"],
            "STOP": ["TAG", "TGA", "TAA"],
            "Thr": ["ACC", "ACA", "ACG", "ACT"],
            "Phe": ["TTT", "TTC"],
            "Ala": ["GCA", "GCC", "GCG", "GCT"],
            "Gly": ["GGT", "GGG", "GGA", "GGC"],
            "Ile": ["ATC", "ATA", "ATT"],
            "Leu": ["TTA", "TTG", "CTC", "CTT", "CTG", "CTA"],
            "His": ["CAT", "CAC"],
            "Arg": ["CGA", "CGC", "CGG", "CGT", "AGG", "AGA"],
            "Trp": ["TGG"],
            "Val": ["GTA", "GTC", "GTG", "GTT"],
            "Glu": ["GAG", "GAA"],
            "Tyr": ["TAT", "TAC"],
        }
        self.aminoacids = [
            "Cys",
            "Asp",
            "Ser",
            "Gln",
            "Met",
            "Asn",
            "Pro",
            "Lys",
            "Thr",
            "Phe",
            "Ala",
            "Gly",
            "Ile",
            "Leu",
            "His",
            "Arg",
            "Trp",
            "Val",
            "Glu",
            "Tyr"
        ]
        if cfg.stop_codon:
            self.aminoacids.append("STOP")
        self.complement = {"A": "T", "C": "G", "G": "C", "T": "A"}

        self.designed_variants = {}


        # First check for unwanted cutsites (BsaI sites and BsmBI sites)
        match_sites = [
                gene.seq.upper().count(cut)
                + gene.seq.upper().count(cut.reverse_complement())
                for cut in cfg.avoid_sequence
            ]
        if any(match_sites):
            raise ValueError(
                "Unwanted Restriction cut sites found. Please input plasmids with these removed."
                + str([cfg.avoid_sequence[i] for i, x in enumerate(match_sites) if bool(x)])
            )  # change codon

        # Check for ORF specification and record start and end
        logger.info("Checking for ORF specification")
        logger.info("Start: " + str(start) + " End: " + str(end))
        if start is not None and end is not None:
            logger.info("Using user-specified ORF")
            logger.info("Start: " + str(start) + " End: " + str(end))
            logger.info("ORF length: " + str(end - start))
            if (end - start) % 3 != 0:
                print("Gene length is not divisible by 3. Resetting and attempting to identify ORF.")
                logger.warning("Gene length is not divisible by 3. Resetting and attempting to identify ORF.")
                start = None
                end = None
        if start is None or end is None:
            logger.info("Start and end of ORF were not provided. Manually identifying ORF.")
            start, end = findORF(
                gene,
                non_interactive=cfg.non_interactive,
                preferred_orf_index=cfg.preferred_orf_index,
            )
            logger.info("Found the following positions: Start: " + str(start) + " End: " + str(end))

        self.aacount = int((end - start) / 3)
        self.start = start
        self.end = end
        logger.info("Using the following ORF positions:")
        logger.info("Start: " + str(start) + " End: " + str(end))


        # record sequence with extra bp to account for primer. for plasmids (circular) we can rearrange linear sequence)
        if start - cfg.primer_buffer < 0:
            self.seq = (
                gene.seq[start + 3 - cfg.primer_buffer:]
                + gene.seq[: end + cfg.primer_buffer]
            )
        elif end + cfg.primer_buffer > len(gene.seq):
            self.seq = (
                gene.seq[start + 3 - cfg.primer_buffer:]
                + gene.seq[: end + cfg.primer_buffer - len(gene.seq)]
            )
        else:
            self.seq = gene.seq[start + 3 - cfg.primer_buffer: end + cfg.primer_buffer]
        self.seq = self.seq.upper()

        # Determine Fragment Size and store beginning and end of each fragment
        num = int(
            round(((end - start - 3) / float(maxfrag)) + 0.499999999)
        )  # total bins needed (rounded up)

        insertionsites = range(start + 3, end, 3)
        fragsize = [len(insertionsites[i::num]) * 3 for i in list(range(num))]

        # if any(x<144 for x in fragsize):
        #     raise ValueError('Fragment size too low')
        print("Initial Fragment Sizes for:" + self.geneid)
        print(fragsize)

        total = cfg.primer_buffer
        breaksites = [cfg.primer_buffer]

        for x in fragsize:
            total += x
            breaksites.extend([total])
        self.breaklist = [
            [x, x + fragsize[idx]] for idx, x in enumerate(breaksites[:-1])
        ]  # insertion site to insertion site
        self.problemsites = set()
        self.unique_Frag = [True] * len(fragsize)
        self.fragsize = fragsize
        self.__breaksites = breaksites


    # Update Breaksites
    @property
    def breaksites(self):
        return self.__breaksites

    @breaksites.setter
    def breaksites(self, value):
        cfg = self.pool.config
        if isinstance(value, list):
            if any([(x - cfg.primer_buffer) % 3 != 0 for x in value]):
                raise ValueError("New Breaksites are not divisible by 3")
            if (
                value[0] != self.breaksites[0] or value[-1] != self.breaksites[-1]
            ) and not cfg.dms:
                if cfg.non_interactive:
                    if cfg.breaksite_change_policy == "error":
                        raise ValueError(
                            "Beginning and end of gene changed during non-interactive run. "
                            "Set breaksite_change_policy to 'warn' to continue."
                        )
                    logger.warning(
                        "Beginning and end of gene changed in non-interactive mode; continuing."
                    )
                else:
                    if (
                        input(
                            "Beginning and End of gene have changed. Are you sure you want to continue? (y/n)"
                        )
                        != "y"
                    ):
                        raise Exception("Canceled user set break sites")
            self.__breaksites = value
            fragsize = [j - i for i, j in zip(value[:-1], value[1:])]
            self.fragsize = fragsize
            self.breaklist = [
                [x, x + fragsize[idx]] for idx, x in enumerate(value[:-1])
            ]  # insertion site to insertion site
            print("New Fragment Sizes for: " + self.geneid)
            print(fragsize)
            # fragment_genes(self)

        else:
            raise ValueError("Breaklist input is not a list")
