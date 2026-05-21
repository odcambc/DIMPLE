"""Fragment layout helpers for the DIMPLE pipeline.

This module owns the internal cluster of functions responsible for deciding
how a gene is split into fragments and for resolving fragment-boundary
conflicts (non-specific primers, matching overhangs).
"""

from __future__ import annotations

from DIMPLE.core import DIMPLE


def recalculate_num_fragments(gene):
    num = int(
        round(((gene.end - gene.start) / float(gene.maxfrag)) + 0.499999999)
    )  # total bins needed (rounded up)
    insertionsites = range(gene.start + 3, gene.end, 3)
    gene.fragsize = [len(insertionsites[i::num]) * 3 for i in list(range(num))]
    total = DIMPLE.primerBuffer
    breaksites = [DIMPLE.primerBuffer]
    for x in gene.fragsize:
        total += x
        breaksites.extend([total])
    # if DIMPLE.dms:
    #     tmpbreaklist = []
    #     for idx, x in enumerate(breaksites[:-1]):
    #         if idx:
    #             tmpbreaklist.append([x, x + gene.fragsize[idx]])
    #         else:
    #             tmpbreaklist.append([x + 3, x + gene.fragsize[idx] + 3])
    #     gene.breaklist = tmpbreaklist
    # else:
    gene.breaklist = [
        [x, x + gene.fragsize[idx]] for idx, x in enumerate(breaksites[:-1])
    ]  # insertion site to insertion site
    #gene.problemsites = set()
    gene.breaksites = breaksites
    gene.unique_Frag = [True] * len(gene.fragsize)
    return gene


def switch_fragmentsize(gene, detectedsite, pool):
    """TODO:
    Docstring
    """

    if not isinstance(gene, DIMPLE):
        raise TypeError("Not an instance of the DIMPLE class")
    skip = False
    count = 0
    count2 = 0
    print("Non specific Fragment:" + str(detectedsite))
    if (
        len(gene.fragsize) * gene.maxfrag < len(gene.seq) - gene.primerBuffer * 2
    ):  # if the maxfrag has changed and it is impossible to split the gene into x number of fragments it should recalculate the number of fragments
        gene = recalculate_num_fragments(gene)
    else:
        gene.problemsites.add(gene.breaksites[detectedsite])
    if all(item == gene.maxfrag for item in gene.fragsize) or any(
        item > gene.maxfrag for item in gene.fragsize
    ):
        gene.maxfrag += -1
    while True:
        if count > len(gene.breaksites):
            # Randomly shift a fragment
            count = 0
            count2 += 1
            detectedsite = gene.rng.integers(
                1, len(gene.breaksites) - 1, dtype=int
            )  # dont change beginning or end
            if gene.fragsize[detectedsite - 1] == gene.maxfrag:
                shift = -3
            else:
                if gene.rng.integers(0, 2, dtype=int):
                    shift = 3
                else:
                    shift = -3
            gene.breaksites[detectedsite] = gene.breaksites[detectedsite] + shift
            gene.fragsize = [
                j - i for i, j in zip(gene.breaksites[:-1], gene.breaksites[1:])
            ]
            # if DIMPLE.dms:
            #     tmpbreaklist = []
            #     for idx, x in enumerate(gene.breaksites[:-1]):
            #         if idx:
            #             tmpbreaklist.append([x, x + gene.fragsize[idx]])
            #         else:
            #             tmpbreaklist.append([x + 3, x + gene.fragsize[idx] + 3])
            #     gene.breaklist = tmpbreaklist
            # else:
            gene.breaklist = [
                [x, x + gene.fragsize[idx]]
                for idx, x in enumerate(gene.breaksites[:-1])
            ]
            if count2 > len(gene.breaklist) * 3:
                gene.maxfrag += -1  # try to change for only this gene...
                if len(gene.fragsize) * gene.maxfrag < len(gene.seq):
                    gene = recalculate_num_fragments(gene)
                    count = 0
                    count2 = 0
        count += 1
        # Find connecting Fragments
        if detectedsite == 0 or detectedsite == len(gene.fragsize):
            print("Issue with primer on end of gene")
            skip = True
            break
        if (
            gene.fragsize[detectedsite] == gene.fragsize[detectedsite - 1]
            and gene.fragsize[detectedsite] >= gene.maxfrag
        ):
            if all(
                item >= gene.maxfrag for item in gene.fragsize[detectedsite + 1 :]
            ) and not all(
                item >= gene.maxfrag for item in gene.fragsize[: detectedsite - 1]
            ):
                shift = 3
                while gene.breaksites[detectedsite] + shift in gene.problemsites:
                    shift += 3
            if all(
                item >= gene.maxfrag for item in gene.fragsize[: detectedsite - 1]
            ) and not all(
                item >= gene.maxfrag for item in gene.fragsize[detectedsite + 1 :]
            ):
                shift = -3
                while gene.breaksites[detectedsite] + shift in gene.problemsites:
                    shift += -3
            else:
                if (
                    detectedsite < len(gene.fragsize) / 2
                ):  # should be based on problemsites not where it is located in the gene
                    shift = 3
                    while gene.breaksites[detectedsite] + shift in gene.problemsites:
                        shift += 3
                else:
                    shift = -3
                    while gene.breaksites[detectedsite] + shift in gene.problemsites:
                        shift += -3
        elif gene.fragsize[detectedsite] > gene.fragsize[detectedsite - 1]:
            shift = 3
            while gene.breaksites[detectedsite] + shift in gene.problemsites:
                shift += 3
        elif gene.fragsize[detectedsite] < gene.fragsize[detectedsite - 1]:
            shift = -3
            while gene.breaksites[detectedsite] + shift in gene.problemsites:
                shift += -3
        elif (
            gene.fragsize[detectedsite] == gene.fragsize[detectedsite - 1]
            and gene.fragsize[detectedsite] < gene.maxfrag
        ):
            shift = -3
            while gene.breaksites[detectedsite] + shift in gene.problemsites:
                shift = -shift
                if shift < 0:
                    shift += -3
        # Process shift and reprocess fragments
        gene.breaksites[detectedsite] = gene.breaksites[detectedsite] + shift
        gene.fragsize = [
            j - i for i, j in zip(gene.breaksites[:-1], gene.breaksites[1:])
        ]
        # if DIMPLE.dms:
        #     tmpbreaklist = []
        #     for idx, x in enumerate(gene.breaksites[:-1]):
        #         if idx:
        #             tmpbreaklist.append([x, x + gene.fragsize[idx] + 3])
        #         else:
        #             tmpbreaklist.append([x + 3, x + gene.fragsize[idx] + 3])
        #     gene.breaklist = tmpbreaklist
        # else:
        gene.breaklist = [
            [x, x + gene.fragsize[idx]]
            for idx, x in enumerate(gene.breaksites[:-1])
        ]
        # recheck for size limit issues
        tmpsite = [
            topidx for topidx, item in enumerate(gene.fragsize) if item > gene.maxfrag
        ]
        if tmpsite:
            # pick which side to adjust
            if tmpsite[0] == len(gene.fragsize):
                detectedsite = tmpsite[0]
            elif tmpsite[0] == 0:
                detectedsite = tmpsite[0] + 1
            elif tmpsite[0] == detectedsite and tmpsite[0] + 1 < len(gene.fragsize):
                detectedsite = tmpsite[0] + 1
            else:
                detectedsite = tmpsite[0]
        else:
            break
    print(gene.fragsize)
    # align all linked genes to the same breaksites
    for tmp in gene.linked:
        pool[tmp].breaksites = gene.breaksites
        pool[tmp].fragsize = gene.fragsize
        pool[tmp].breaklist = gene.breaklist
    return skip


def check_overhangs(gene, pool, overlap_l, overlap_r):
    """TODO:
    Docstring
    """
    # Force all overhangs to be different within a gene (no more than 2 matching in a row)
    switched = False
    if not isinstance(gene, DIMPLE):
        raise TypeError("Not an instance of the DIMPLE class")
    while True:
        detectedsites = set()  # stores matching overhangs
        for idx, y in enumerate(gene.breaklist):
            overhang_F = gene.seq[
                y[0] - DIMPLE.cutsite_overhang - overlap_l: y[0] - overlap_r
            ]  # Forward overhang
            overhang_R = gene.seq[
                y[1] + overlap_l: y[1] + DIMPLE.cutsite_overhang + overlap_r
            ]  # Reverse overhang
            if (
                overhang_F == overhang_R
                or overhang_F == overhang_R.reverse_complement()
            ):
                detectedsites.update([idx])
        # overhang = []
        # for idx, y in enumerate(gene.breaklist):
        #     overhang.append([gene.seq[y[0] - DIMPLE.cutsite_overhang - overlapL: y[0] - overlapR], idx])  # Forward overhang
        #     overhang.append([gene.seq[y[1] + overlapL: y[1] + DIMPLE.cutsite_overhang + overlapR], idx + 1])  # Reverse overhang
        # detectedsites = set()  # stores matching overhangs
        # for i in range(len(overhang)):  # check each overhang for matches
        #     for j in [x for x in range(len(overhang)) if x != i]:  # permutate over every overhang combination to find matches
        #         #if overhang[i][0] == overhang[j][0] or overhang[i][0][:3] == overhang[j][0][:3] or overhang[i][0][1:] == overhang[j][0][1:] or overhang[i][0] == overhang[i][0].reverse_complement():  # no 3 matching bases
        #         if overhang[i][0] == overhang[j][0] or overhang[i][0] == overhang[i][0].reverse_complement():  # no 3 matching bases
        #             detectedsites.update([overhang[i][1]])
        for detectedsite in detectedsites:
            switched = True
            if detectedsite == 0:
                detectedsite = 1  # don't mess with the first cut site
            print(
                "------------------ Fragment size swapped due to matching overhangs ------------------"
            )
            skip = switch_fragmentsize(gene, detectedsite, pool)
        else:  # if no detected sites
            break
    return switched
