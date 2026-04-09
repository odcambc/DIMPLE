# RUN DIMPLE
# script for usage with command line

import argparse
import os
import warnings
from datetime import datetime

import logging

from DIMPLE.DIMPLE import (
    DIMPLE,
    addgene,
    align_genevariation,
    generate_DMS_fragments,
    post_qc,
    print_all,
)
from DIMPLE.run_settings import (
    apply_barcode_start,
    apply_handle,
    apply_instance_settings,
    apply_random_seed,
    apply_restriction_settings,
    compute_overlaps_and_maxfrag,
    configure_dimple_logging,
    normalize_avoid_list,
    resolve_codon_usage,
    validate_insertions,
)
from DIMPLE.utilities import parse_custom_mutations

log_file = "logs/DIMPLE-{:%Y-%m-%d-%s}.log".format(datetime.now())
configure_dimple_logging(log_file)
logger = logging.getLogger(__name__)
logger.info("Started logging")

parser = argparse.ArgumentParser(
    description="DIMPLE: Deep Indel Missense Programmable Library Engineering"
)
parser.add_argument(
    "-wDir",
    help="Working directory for fasta files and output folder",
)
parser.add_argument(
    "-geneFile",
    required=True,
    help=(
        "Input all gene sequences including backbone in a fasta format. "
        "Place all in one fasta file. Name description can include start and end "
        "points (>gene1 start:1 end:2)"
    ),
)
parser.add_argument(
    "-handle",
    default="AGCGGGAGACCGGGGTCTCTGAGC",
    help=(
        "Genetic handle for domain insertion. This is important for defining the "
        "linker. Currently uses BsaI (4 base overhang), but this can be swapped "
        "for SapI (3 base overhang)."
    ),
)
parser.add_argument(
    "-dis",
    default=False,
    help="use the handle to insert domains at every position in POI",
)
parser.add_argument(
    "-matchSequences",
    action="store_const",
    const="match",
    default="nomatch",
    help=(
        "Find similar sequences between genes to avoid printing the same oligos "
        "multiple times. Default: No matching"
    ),
)
parser.add_argument(
    "-oligoLen",
    type=int,
    default=230,
    help="Synthesized oligo length",
)
parser.add_argument(
    "-fragmentLen",
    default=0,
    type=int,
    help="Maximum length of gene fragment",
)
parser.add_argument(
    "-overlap",
    default=4,
    type=int,
    help=(
        "Enter number of bases to extend each fragment for overlap. This will help "
        "with insertions close to fragment boundary"
    ),
)
parser.add_argument(
    "-DMS",
    action="store_const",
    const=True,
    default=False,
    help="Choose if you will run deep deep mutation scan",
)
parser.add_argument(
    "-custom_mutations",
    default=None,
    help="Path to file that includes custom mutations with the format position:AA",
)
parser.add_argument(
    "-usage",
    default="human",
    help='Default is "human" or "ecoli", or pass a path to a codon usage file',
)
parser.add_argument(
    "-insertions",
    default=False,
    nargs="+",
    help=(
        "Enter a list of insertions (nucleotides) to make at every position. Note, "
        "you should enter multiples of 3 nucleotides to maintain reading frame"
    ),
)
parser.add_argument(
    "-deletions",
    default=False,
    nargs="+",
    help=(
        "Enter a list of deletions (number of nucleotides) to symmetrically delete "
        "(it will make deletions in multiples of 2x). Note you should enter "
        "multiples of 3 to maintain reading frame"
    ),
)
parser.add_argument(
    "-include_substitutions",
    default=False,
    help="If you are running DMS but only want to insert or delete AA",
)
parser.add_argument(
    "-barcode_start",
    default=0,
    help=(
        "To run DIMPLE multiple times, you will need to avoid using the same "
        "barcodes. This allows you to start at a different barcode."
    ),
)
parser.add_argument(
    "-restriction_sequence",
    default="CGTCTC(G)1/5",
    help=(
        "Recommended using BsmBI - CGTCTC(G)1/5 or BsaI - GGTCTC(G)1/5. Do not use N"
    ),
)
parser.add_argument(
    "-avoid_sequence",
    nargs="+",
    default=["CGTCTC", "GGTCTC"],
    help=(
        "Avoid these sequences in the backbone - BsaI and BsmBI. For multiple "
        "sequences use a space between inputs. Example -avoid_sequence CGTCTC GGTCTC"
    ),
)
parser.add_argument(
    "-include_stop_codons",
    help="Include stop codons in the list of scanning mutations.",
    default=False,
    const=True,
    action="store_const",
)
parser.add_argument(
    "-include_synonymous",
    help="Include synonymous codons in the list of scanning mutations.",
    default=False,
    const=True,
    action="store_const",
)
parser.add_argument(
    "-make_double",
    help="Make each combination of mutations within a fragment",
    default=False,
    const=True,
    action="store_const",
)
parser.add_argument(
    "-maximize_nucleotide_change",
    help=(
        "Maximize the number of nucleotide changes in each codon for easier "
        "detection in NGS"
    ),
    default=False,
    const=True,
    action="store_const",
)
parser.add_argument(
    "-seed",
    help="Seed for random number generation",
    default=None,
)
args = parser.parse_args()

if args.wDir is None:
    if "/" in args.geneFile:
        args.wDir = args.geneFile.rsplit("/", 1)[0] + "/"
        args.geneFile = args.geneFile.rsplit("/", 1)[1]
    else:
        args.wDir = ""

apply_handle(args.handle)

deletions_for_overlap = args.deletions if args.deletions else False
overlap_l, overlap_r = compute_overlaps_and_maxfrag(
    args.oligoLen,
    args.fragmentLen,
    args.overlap,
    deletions_for_overlap,
    logger=logger,
)

apply_barcode_start(int(args.barcode_start))

apply_restriction_settings(args.restriction_sequence)

normalize_avoid_list(args.avoid_sequence, logger=logger)

DIMPLE.dms = args.DMS
DIMPLE.stop_codon = args.include_stop_codons
DIMPLE.make_double = args.make_double
DIMPLE.maximize_nucleotide_change = args.maximize_nucleotide_change

if args.custom_mutations:
    with open(args.custom_mutations, encoding="utf-8") as f:
        custom_mutations = f.readlines()
    custom_mutations = parse_custom_mutations(custom_mutations)
else:
    custom_mutations = None

if args.seed:
    apply_random_seed(int(args.seed))
else:
    apply_random_seed(None)

resolve_codon_usage(args.usage)

OLS = addgene(os.path.join(args.wDir, args.geneFile).strip())

apply_instance_settings(OLS)

if args.matchSequences == "match":
    align_genevariation(OLS)
if args.deletions:
    args.deletions = [int(x) for x in args.deletions]
if not any([DIMPLE.dms, args.insertions, args.deletions]):
    raise ValueError("Didn't select any mutations to generate")

if args.insertions:
    validate_insertions(list(args.insertions))

logger.info("Generating DMS fragments")

generate_DMS_fragments(
    OLS,
    overlap_l,
    overlap_r,
    args.include_synonymous,
    custom_mutations,
    DIMPLE.dms,
    args.insertions,
    args.deletions,
    args.dis,
    args.wDir,
)

post_qc(OLS)
print_all(OLS, args.wDir)

logger.info("Finished")
