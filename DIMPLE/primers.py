"""Primer design helpers for the DIMPLE pipeline.

Gene-primer and fragment-primer design plus a nonspecific-annealing check,
migrated out of ``DIMPLE.DIMPLE`` so primer logic owns its own module.
"""

from __future__ import annotations

import logging

from Bio.SeqUtils import MeltingTemp as mt

from DIMPLE.core import DIMPLE

logger = logging.getLogger(__name__)


def find_geneprimer(genefrag, start, end):
    # 3' end of primer is variable to adjust melting temperature
    # 5' end of primer is fixed, with restriction site added
    # Also add space for maximum deletion on 5' end
    primer = (
        genefrag[start:end].complement() + DIMPLE.cutsite[::-1] + "ATA"
    )  # added ATA for cleavage close to end of DNA fragment
    # Check melting temperature
    # find complementary sequences
    comp = 0  # compensate for bases that align with bsmbi
    while primer.complement()[end - start + comp] == genefrag[end + comp]:
        comp += 1
    # comp += 1 # This is important for single basepair overhang
    tm2 = mt.Tm_NN(primer[0: end - start + comp], nn_table=mt.DNA_NN2)
    tm4 = mt.Tm_NN(primer[0: end - start + comp], nn_table=mt.DNA_NN4)
    count = 0
    while (
        tm2 < DIMPLE.gene_primerTm[0]
        or tm2 > DIMPLE.gene_primerTm[1]
        or tm4 < DIMPLE.gene_primerTm[0]
        or tm4 > DIMPLE.gene_primerTm[1]
    ):
        if tm2 < DIMPLE.gene_primerTm[0] or tm4 < DIMPLE.gene_primerTm[0]:
            start += -1
            primer = (
                genefrag[start:end].complement() + DIMPLE.cutsite[::-1] + "ATA"
            )  # cut site addition
            tm2 = mt.Tm_NN(primer[0: end - start + comp], nn_table=mt.DNA_NN2)
            tm4 = mt.Tm_NN(primer[0: end - start + comp], nn_table=mt.DNA_NN4)
        if (
            count > 12 or start == 0
        ):  # stop if caught in inf loop or if linker is at max (31 + 7 = 38 bases)
            break
        if tm2 > DIMPLE.gene_primerTm[1] and tm4 > DIMPLE.gene_primerTm[1]:
            start += 1
            primer = genefrag[start:end].complement() + DIMPLE.cutsite[::-1] + "ATA"
            # tm = mt.Tm_NN(primer[0:e-s+comp],c_seq=genefrag[s:e+comp],nn_table=mt.DNA_NN2)
            tm2 = mt.Tm_NN(primer[0: end - start + comp], nn_table=mt.DNA_NN2)
            tm4 = mt.Tm_NN(primer[0: end - start + comp], nn_table=mt.DNA_NN4)
        count += 1
    # optional - force first nucleotide to a C or G
    # while primer[0]=="T" or primer[0]=="A" or primer[0]=="t" or primer[0]=="a":
    #     s += -1
    #     primer = genefrag[s:e].complement()+"CTCTGCA"
    #     tm = mt.Tm_NN(primer[0:e-s+comp],nn_table=mt.DNA_NN2)
    # return final primer with tm
    print(
        "Generated primers: ",
        primer.complement().reverse_complement(),
        round(tm2, 1),
        start,
    )
    return primer.complement().reverse_complement(), round(tm2, 1), start


def find_fragment_primer(fragment, stop):
    # This function finds optimal primer for OLS subpool by changing 3' end of primer
    start = 0  # starts at maximum length (5' is fixed)
    if stop > 25:  # limit primer to 25 bases to begin with
        end = 25
    else:
        end = stop
    count = 0
    primer = fragment[start:end]
    tm2 = mt.Tm_NN(
        primer, nn_table=mt.DNA_NN2
    )  # Two methods of finding melting temperature seems more consistent
    tm4 = mt.Tm_NN(primer, nn_table=mt.DNA_NN4)
    while (
        tm2 < DIMPLE.primerTm[0]
        or tm2 > DIMPLE.primerTm[1]
        or tm4 < DIMPLE.primerTm[0]
        or tm4 > DIMPLE.primerTm[1]
        or len(primer) < 16
    ):
        count += 1
        if (
            count > 12 or end > stop
        ):  # stop if caught in inf loop or if primer is larger than the barcode
            end = stop
            primer = fragment[start:end]
            break

        if tm2 < DIMPLE.primerTm[0] or tm4 < DIMPLE.primerTm[0]:
            if start == 0:
                break
            end += 1
            primer = fragment[start:end]
            tm2 = mt.Tm_NN(primer, nn_table=mt.DNA_NN2)
            tm4 = mt.Tm_NN(primer, nn_table=mt.DNA_NN4)

        if tm2 > DIMPLE.primerTm[1] or tm4 > DIMPLE.primerTm[1]:
            end += -1
            primer = fragment[start:end]
            tm2 = mt.Tm_NN(primer, nn_table=mt.DNA_NN2)
            tm4 = mt.Tm_NN(primer, nn_table=mt.DNA_NN4)
    return primer, round(tm2, 1)


def check_nonspecific(primer, fragment, point):
    non = []
    # fragment is the entire gene sequence plus the buffer sequence on each side
    # point is the position of the primer in the fragment
    # Forward
    for i in range(len(fragment) - len(primer)):  # Scan each position
        # first check if the primer binds at each position in the fragment
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
        # if the primer binds to 80% of the first ... bases
        # and more than 6 bases
        # and the 3' matches
        # then check melting temperature
        if (
            sum(match[first:]) > len(primer[first:]) * 0.8
            and sum(match[first:]) > 6
            and match[-1]
            and point != i
        ):  # string compare - sum of matched nt is greater than 80%
            try:
                # check the melting temperature of the primer
                melt = mt.Tm_NN(
                    primer[first:],
                    c_seq=fragment[i + first : i + len(primer)].complement(),
                    nn_table=mt.DNA_NN2,
                    de_table=mt.DNA_DE1,
                    imm_table=mt.DNA_IMM1,
                )
                if melt > 25:
                    print("Found non-specific match at " + str(i + 1) + "bp:")
                    print("match: " + fragment[i: i + len(primer)])
                    print("primer:" + primer + " Tm:" + str(round(melt, 1)))
                    logger.warning("Found non-specific match at " + str(i + 1) + "bp:")
                    logger.warning("match: " + fragment[i: i + len(primer)])
                    logger.warning("primer:" + primer + " Tm:" + str(round(melt, 1)))
                if melt > 35:
                    non.append(True)
            except ValueError as valerr:
                print(
                    str(valerr)
                    + ". Please check position manually:"
                    + str(i + 1)
                    + " forward"
                )
                print("Primer:" + primer)
                print("Match: " + fragment[i : i + len(primer)])
                non.append(False)
    # Reverse
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
            and point != -i
        ):  # string compare - sum of matched nt is greater than 80%
            try:
                melt = mt.Tm_NN(
                    primer[first:],
                    c_seq=fragment[i + first : i + len(primer)].complement(),
                    nn_table=mt.DNA_NN2,
                    de_table=mt.DNA_DE1,
                    imm_table=mt.DNA_IMM1,
                )
                if melt > 20:
                    print("Found non-specific match at " + str(i + 1) + "bp:")
                    print(" match:" + fragment[i : i + len(primer)])
                    print("primer:" + primer + " Tm:" + str(melt))
                if melt > 35:
                    non.append(True)
            except ValueError as valerr:
                print(
                    str(valerr)
                    + ". Please check position manually:"
                    + str(i + 1)
                    + " reverse"
                )
                print("Primer:" + primer)
                print("Match: " + fragment[i : i + len(primer)])
                non.append(False)
    return sum(non)
