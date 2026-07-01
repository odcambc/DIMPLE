"""QC helpers for the DIMPLE pipeline.

Contains post-run QC checks for barcode primer specificity and simulated
final-assembly verification.
"""

from __future__ import annotations

import logging
import re
from collections import Counter

from Bio.Restriction import BsaI, BsmBI
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import MeltingTemp as mt
from pydna.amplify import pcr
from pydna.dseqrecord import Dseqrecord

from DIMPLE.core import DIMPLE

logger = logging.getLogger(__name__)


_DROP_FRACTION_LIMIT = 0.02  # >2% dropped signals a systematic layout problem


def report_dropped_oligos(pool):
    """Summarise oligos dropped for edit-induced internal restriction sites.

    A handful of drops is expected (some variants inherently carry the enzyme
    site). Make it loud regardless, and hard-fail if more than a small fraction
    is lost -- a large or fragment-concentrated drop means a systematic problem
    (e.g. a whole subpool whose overhang reconstitutes the site) rather than a few
    unavoidable variants, and must not slip through as a routine warning.
    """
    dropped = [(g.geneid, vid) for g in pool for vid in getattr(g, "dropped_oligos", [])]
    if not dropped:
        return
    total = sum(len(g.designed_variants) + len(getattr(g, "dropped_oligos", [])) for g in pool)
    by_subpool = Counter(re.search(r"-(\d+)_", vid).group(1) for _, vid in dropped if re.search(r"-(\d+)_", vid))
    msg = (
        f"Dropped {len(dropped)} of {total} oligos with an internal restriction site "
        f"(cannot assemble). By subpool: {dict(sorted(by_subpool.items()))}. "
        f"Examples: {[vid for _, vid in dropped[:8]]}"
    )
    logger.warning("=" * 78)
    logger.warning(msg)
    logger.warning("=" * 78)
    if total and len(dropped) / total > _DROP_FRACTION_LIMIT:
        raise Exception(
            f"Too many oligos dropped for internal restriction sites: {len(dropped)}/{total} "
            f"(> {_DROP_FRACTION_LIMIT:.0%}). This indicates a systematic layout problem "
            f"(likely a subpool overhang reconstituting {pool.config.cutsite}), not a few "
            f"unavoidable variants. By subpool: {dict(sorted(by_subpool.items()))}."
        )


def post_qc(pool):
    logger.info("Running post QC")
    if not isinstance(pool[0], DIMPLE):
        raise TypeError("Not an instance of the DIMPLE class")
    report_dropped_oligos(pool)
    # Post QC
    all_oligos = []
    all_barPrimers = []
    for obj in pool:
        logger.info(f"Running QC for {obj.geneid}")
        try:
            all_oligos.extend(obj.oligos)
            all_barPrimers.extend(obj.barPrimer)
        except AttributeError:
            print(obj.geneid + " has not been processed")

        if pool.config.enzyme is not None:
            logger.info(f"Checking oligo assembly for {obj.geneid}")
            check_final_assembly(obj)
        else:
            logger.info(
                f"Skipping oligo assembly check for {obj.geneid} as no enzyme is specified."
            )

    print("Running QC for barcode primer specificity")
    cassetteSet = set(all_oligos[0].id[:-6])
    uCassette = [SeqRecord(all_oligos[0].seq, id=all_oligos[0].id[:-6])]
    for idx in range(len(all_oligos)):
        if all_oligos[idx].id[:-6] not in cassetteSet:
            uCassette.append(SeqRecord(all_oligos[idx].seq, id=all_oligos[idx].id[:-6]))
        cassetteSet.add(all_oligos[idx].id[:-6])
    grouped = iter(all_barPrimers)
    grouped = zip(grouped, grouped)  # create the combinatorial comparisons
    nonspecific = []
    # iterate over every barcode primer pair and match to each oligo to check for
    # nonspecific amplification
    for idxPrime, primers in enumerate(grouped):  # iterate over every barcode primer pair
        print("Checking primer set:" + primers[0].id[:-2])
        for idxCassette, fragment in enumerate(uCassette):  # iterate over every pool oligo
            if (
                primers[0].id.split("_")[2] != all_oligos[idxCassette].id.split("_")[2]
            ):  # ignore designed annealing (same name)
                fragname = fragment.id
                fragment = fragment.seq
                non = [[False], [False]]
                for idxDirection, primer in enumerate(primers):
                    primername = primer.id
                    primer = primer.seq
                    for i in range(
                        len(fragment) - len(primer)
                    ):  # iterate over the length of the oligo for a binding site
                        match = [
                            primer[j].lower() == fragment[i + j].lower() for j in range(len(primer))
                        ]
                        first = 10
                        for k in range(len(match) - 3):
                            if (match[k] and match[k + 1] and match[k + 3]) or (
                                match[k] and match[k + 1] and match[k + 2]
                            ):
                                first = k
                                break
                        if (
                            sum(match[first:]) > len(primer[first:]) * 0.8
                            and sum(match[first:]) > 6
                            and match[-1]
                        ):  # string compare - sum of matched nt is greater than 80%
                            try:
                                melt = mt.Tm_NN(
                                    primer[first:],
                                    c_seq=fragment[i + first : i + len(primer)].complement(),
                                    nn_table=mt.DNA_NN2,
                                    de_table=mt.DNA_DE1,
                                    imm_table=mt.DNA_IMM1,
                                )
                                # if melt>20:
                                #     print('Found non-specific match:'
                                #           + fragment[i:i+len(primer)])
                                #     print('                 primer:'
                                #           + primer + ' Tm:' + str(round(melt, 1)))
                                if melt > 35:
                                    non[0].append(True)
                            except ValueError:
                                # print(str(valerr) + ". Please check position manually:"
                                #       + str(i+1) + " forward")
                                # print('Primer:'+primer)
                                # print('Match: '+fragment[i:i+len(primer)])
                                pass
                    fragment = fragment.reverse_complement()
                    for i in range(len(fragment) - len(primer)):
                        match = [
                            primer[j].lower() == fragment[i + j].lower() for j in range(len(primer))
                        ]
                        first = 10
                        for k in range(0, len(match) - 3, 1):
                            if match[k] and match[k + 1] and match[k + 3]:
                                first = k
                                break
                        if (
                            sum(match[first:]) > len(primer[first:]) * 0.8
                            and sum(match[first:]) > 6
                            and match[-1]
                        ):  # string compare - sum of matched nt is greater than 80%
                            try:
                                melt = mt.Tm_NN(
                                    primer[first:],
                                    c_seq=fragment[i + first : i + len(primer)].complement(),
                                    nn_table=mt.DNA_NN2,
                                    de_table=mt.DNA_DE1,
                                    imm_table=mt.DNA_IMM1,
                                )
                                # if melt > 20:
                                #     print('Found non-specific match:'+fragment[i:i+len(primer)])
                                #     print('                 primer:'+primer+' Tm:'+str(melt))
                                if melt > 35:
                                    non[1].append(True)
                            except ValueError:
                                # print(str(valerr) + ". Please check position manually:"
                                #       + str(i+1) + " reverse")
                                # print('Primer:'+primer)
                                # print('Match: '+fragment[i:i+len(primer)])
                                pass
                    if sum(non[0]) == 0 and sum(non[1]) == 0:
                        break
                    if sum(non[0]) > 0 and sum(non[1]) > 0:
                        nonspecific.append([primername, fragname])
                        print("Found Non-specific Amplification")
    if nonspecific:
        print("Nonspecific Primers: (Manually changing primer sequence recommended)")
        print(nonspecific)
    else:
        print("No non-specific primers detected")


def _pcr_or_friendly_error(fwd, rev, template, gene_id, context):
    """Simulate PCR, re-raising pydna's raw failure as actionable DIMPLE guidance.

    pydna raises a bare ``PCR not specific! ... anneals ... at N`` (or a no-product)
    error that gives a DIMPLE user no lever. This names the gene/fragment and points
    at the fix -- a designed primer that anneals at more than one site usually means
    the oligo is too short to place a unique primer; a longer oligo (or adjusted
    gene-primer Tm) resolves it.
    """
    try:
        product = pcr(fwd, rev, template)
    except ValueError as exc:  # pydna raises ValueError for non-specific / no product
        raise ValueError(
            f"Primer specificity check failed for {gene_id} ({context}): a designed "
            f"primer anneals at more than one site, so the fragment cannot be "
            f"amplified cleanly. This usually means the oligo length is too short to "
            f"place a unique primer -- try a longer oligo length, or adjust the "
            f"gene-primer Tm bounds. (pydna: {exc})"
        ) from exc
    return Dseqrecord(product)


def check_final_assembly(gene):
    """Test that each oligo assembles properly and contains the designed mutation."""

    cfg = gene.pool.config
    # Check whether the enzyme is set.
    if cfg.enzyme:
        if cfg.enzyme == "BsmBI":
            enzyme = BsmBI
        elif cfg.enzyme == "BsaI":
            enzyme = BsaI
        else:
            logger.warning("Enzyme not recognized. Not performing assembly check.")
            return None
    else:
        logger.warning("No enzyme set. Not performing assembly check.")
        return None

    n_fragments = len(gene.genePrimer) // 2
    full_template = Dseqrecord(gene.seq, circular=True)

    backbones = []
    oligo_primer_dseqs = []
    logger.info(f"Testing assembly for {gene.geneid}")
    logger.info(f"Using enzyme: {cfg.enzyme}")
    logger.info(f"Number of fragments: {n_fragments}")

    for frag in range(0, n_fragments):
        fwd_primer = Dseqrecord(gene.genePrimer[frag * 2])
        rev_primer = Dseqrecord(gene.genePrimer[frag * 2 + 1])
        template_pcr_product = _pcr_or_friendly_error(
            fwd_primer, rev_primer, full_template, gene.geneid, f"gene primer, fragment {frag + 1}"
        )
        cut_template_product = max(template_pcr_product.cut(enzyme), key=len)
        backbones.append(cut_template_product)

        fwd_oligo_primer = Dseqrecord(gene.barPrimer[frag * 2])
        rev_oligo_primer = Dseqrecord(gene.barPrimer[frag * 2 + 1])
        oligo_primer_dseqs.append((fwd_oligo_primer, rev_oligo_primer))

    for variant in gene.designed_variants:
        variant_dict = gene.designed_variants[variant]
        # Get the fragment that the variant is in.
        fragment = variant_dict["fragment"]
        sequence = variant_dict["xfrag"]
        oligo_sequence = variant_dict["oligo_sequence"]
        # Simulate PCR of oligo with oligo primers.
        fwd_oligo_primer, rev_oligo_primer = oligo_primer_dseqs[fragment - 1]
        oligo_pcr_product = _pcr_or_friendly_error(
            fwd_oligo_primer,
            rev_oligo_primer,
            oligo_sequence,
            gene.geneid,
            f"oligo primer, variant {variant}",
        )

        try:
            cut_oligo_product = max(oligo_pcr_product.cut(enzyme), key=len)
        except ValueError as error:
            logger.error(f"Oligo cut site issue: {variant}")
            logger.error(str(error))

        try:
            assembled = (cut_oligo_product + backbones[fragment - 1]).looped()
            if str(sequence[4:-4]) not in str(assembled.seq):
                logger.error(f"Assembly product incorrect: {variant}.")
                logger.error(f"Expected variant to contain: {str(sequence[4:-4])}")
                logger.error(f"Predicted assembly product: {str(assembled.seq)}")

        except TypeError as error:
            logger.error(f"Oligo does not assemble with template: {variant}")
            logger.error(str(error))
