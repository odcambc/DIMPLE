"""Core domain model for DIMPLE: the ``DIMPLE`` class and the ``addgene`` loader.

The ``DIMPLE`` class holds both per-gene state and class-level library-design
configuration (mutable class attributes set by callers / ``run_settings``).
That conflation is a known smell; separating per-gene state from global
configuration is tracked as a follow-up refactor. This module is a leaf in the
package import graph -- it depends only on ``DIMPLE.utilities`` and biopython,
so functional modules (``primers``, ``qc``, ...) and the ``DIMPLE.DIMPLE``
orchestrator can import the class from here without an import cycle.
"""

from __future__ import annotations

import logging
import os

import numpy as np
from Bio import SeqIO

from DIMPLE.utilities import findORF

logger = logging.getLogger(__name__)


def addgene(genefile, start=None, end=None):
    """Generate a list of DIMPLE classes from a fasta file containing genes."""
    print("Barcode: " + str(DIMPLE.barcodeF[0].seq))
    print("Number of barcodes: " + str(len(DIMPLE.barcodeF)))
    if start is None:
        start = []
    if end is None:
        end = []
    tmpgene = list(SeqIO.parse(genefile.replace("\\", ""), "fasta"))
    tmpgene[0].seq = tmpgene[0].seq.upper()
    tmpOLS = []
    for gene in tmpgene:
        if "start:" in gene.description and "end:" in gene.description:
            start = int(gene.description.split("start:")[1].split(" ")[0]) - 1
            end = int(gene.description.split("end:")[1].split(" ")[0])
            gene.filename = genefile.replace("\\", "")
            logger.info("Found start: " + str(start) + " and end: " + str(end))
            logger.info("Inferred ORF sequence: " + str(gene.seq[start:end]))
            logger.info("ORF translation: " + str(gene.seq[start:end].translate()))
            tmpOLS.append(DIMPLE(gene, start, end))
        else:
            gene.filename = genefile.replace("\\", "")
            tmpOLS.append(DIMPLE(gene, start, end))
    return tmpOLS  # only return the class object itself if one gene is given


class DIMPLE:
    """Class for generating indel mutagenic scanning libraries."""

    # Calculate and update maxfrag - Max number of nucleotides that a fragment can carry
    @property
    def synth_len(self):
        return self._synth_len

    @synth_len.setter
    def synth_len(self, value):
        self._synth_len = value
        self.maxfrag = value - self.maxfrag_offset

    random_seed = 0
    non_interactive = False
    preferred_orf_index = None
    link_policy = "prompt"  # prompt|always|never
    breaksite_change_policy = "prompt"  # prompt|warn|error

    # Shared variables for all genes
    # Number of nucleotides in synthesis length to preserve for cutsites and primers. Cutsites are
    # composed of the cutsite, the cutsite buffer, and the cutsite overhang
    # Length of cutsite is
    # len_cutsite = len(DIMPLE.cutsite) + len(DIMPLE.cutsite_buffer) + DIMPLE.cutsite_overhang
    # Max oligo primer pair length is 2*21 = 42
    # len_cutsite = len(self.cutsite) + len(self.cutsite_buffer) + self.cutsite_overhang = 22
    maxfrag_offset = 64
    minfrag = 24  # Picked based on smallest size for golden gate fragment efficiency
    primerBuffer = 30  # This extends the sequence beyond the ORF for a primer. Must be greater than 30
    allhangF = []
    allhangR = []
    primerTm = (56.5, 60)  # Melting temperature limits for primers
    gene_primerTm = (58, 62)  # Help gene primer amplification
    # BsaI / BsmBI / None; set by run_settings or callers before pipeline
    enzyme = None
    # Load Barcodes
    dataDirectory = os.path.abspath(os.path.dirname(__file__))
    try:
        barcodeF = list(
            SeqIO.parse(dataDirectory + "/data/forward_finalprimers.fasta", "fasta")
        )
        barcodeR = list(
            SeqIO.parse(dataDirectory + "/data/reverse_finalprimers.fasta", "fasta")
        )
    except FileNotFoundError as exc:
        raise ValueError(
            "Could not find barcode files. Please upload your own or place standard barcodes in the data file."
        ) from exc

    def __init__(self, gene, start=None, end=None):

        # Set up random number generator
        self.rng = np.random.default_rng(DIMPLE.random_seed)

        #  Search for ORF
        try:
            DIMPLE.maxfrag  # if DIMPLE.maxfrag doesn't exist, create it
        except AttributeError:
            DIMPLE.maxfrag = (
                self.synth_len - self.maxfrag_offset
            )  # based on space for barcodes, cut sites, handle. Doesn't need to be exact

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
        if self.stop_codon:
            self.aminoacids.append("STOP")
        self.complement = {"A": "T", "C": "G", "G": "C", "T": "A"}

        self.designed_variants = {}


        # First check for unwanted cutsites (BsaI sites and BsmBI sites)
        match_sites = [
                gene.seq.upper().count(cut)
                + gene.seq.upper().count(cut.reverse_complement())
                for cut in DIMPLE.avoid_sequence
            ]
        if any(match_sites):
            raise ValueError(
                "Unwanted Restriction cut sites found. Please input plasmids with these removed."
                + str([DIMPLE.avoid_sequence[i] for i, x in enumerate(match_sites) if bool(x)])
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
                non_interactive=DIMPLE.non_interactive,
                preferred_orf_index=DIMPLE.preferred_orf_index,
            )
            logger.info("Found the following positions: Start: " + str(start) + " End: " + str(end))

        self.aacount = int((end - start) / 3)
        self.start = start
        self.end = end
        logger.info("Using the following ORF positions:")
        logger.info("Start: " + str(start) + " End: " + str(end))


        # record sequence with extra bp to account for primer. for plasmids (circular) we can rearrange linear sequence)
        if start - self.primerBuffer < 0:
            self.seq = (
                gene.seq[start + 3 - self.primerBuffer:]
                + gene.seq[: end + self.primerBuffer]
            )
        elif end + self.primerBuffer > len(gene.seq):
            self.seq = (
                gene.seq[start + 3 - self.primerBuffer:]
                + gene.seq[: end + self.primerBuffer - len(gene.seq)]
            )
        else:
            self.seq = gene.seq[start + 3 - self.primerBuffer: end + self.primerBuffer]
        self.seq = self.seq.upper()

        # Determine Fragment Size and store beginning and end of each fragment
        num = int(
            round(((end - start - 3) / float(DIMPLE.maxfrag)) + 0.499999999)
        )  # total bins needed (rounded up)

        insertionsites = range(start + 3, end, 3)
        fragsize = [len(insertionsites[i::num]) * 3 for i in list(range(num))]

        # if any(x<144 for x in fragsize):
        #     raise ValueError('Fragment size too low')
        print("Initial Fragment Sizes for:" + self.geneid)
        print(fragsize)

        total = DIMPLE.primerBuffer
        breaksites = [DIMPLE.primerBuffer]

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


    def ochre(self):
        if len(self.SynonymousCodons["STOP"]) < 2:
            raise Exception("You have removed all stop codons")
        self.usage_ecoli["TAG"] = 1
        self.usage_human["TAG"] = 1
        # SynonymousCodons['STOP'] = ['TGA','TAA']
        del self.SynonymousCodons["STOP"][0]
        self.SynonymousCodons["OCHRE"] = ["TAG"]
        self.aminoacids.extend("OCHRE")

    def amber(self):
        if len(self.SynonymousCodons["STOP"]) < 2:
            raise Exception("You have removed all stop codons")
        self.usage_ecoli["TAA"] = 1
        self.usage_human["TAA"] = 1
        del self.SynonymousCodons["STOP"][2]
        self.SynonymousCodons["AMBER"] = ["TAA"]
        self.aminoacids.extend("AMBER")

    def opal(self):
        if len(self.SynonymousCodons["STOP"]) < 2:
            raise Exception("You have removed all stop codons")
        self.usage_ecoli["TGA"] = 1
        self.usage_human["TGA"] = 1
        del self.SynonymousCodons["STOP"][1]
        self.SynonymousCodons["OPAL"] = ["TGA"]
        self.aminoacids.extend("OPAL")

    def __getitem__(self):
        return

    # Update Breaksites
    @property
    def breaksites(self):
        return self.__breaksites

    @breaksites.setter
    def breaksites(self, value):
        if isinstance(value, list):
            if any([(x - DIMPLE.primerBuffer) % 3 != 0 for x in value]):
                raise ValueError("New Breaksites are not divisible by 3")
            if (
                value[0] != self.breaksites[0] or value[-1] != self.breaksites[-1]
            ) and not DIMPLE.dms:
                if DIMPLE.non_interactive:
                    if DIMPLE.breaksite_change_policy == "error":
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
