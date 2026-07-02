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
    total = gene.pool.config.primer_buffer
    breaksites = [gene.pool.config.primer_buffer]
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
    # gene.problemsites = set()
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
        len(gene.fragsize) * gene.maxfrag < len(gene.seq) - gene.pool.config.primer_buffer * 2
    ):  # if the maxfrag has changed and it is impossible to split the gene into x number of
        # fragments it should recalculate the number of fragments
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
            gene.fragsize = [j - i for i, j in zip(gene.breaksites[:-1], gene.breaksites[1:])]
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
                [x, x + gene.fragsize[idx]] for idx, x in enumerate(gene.breaksites[:-1])
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
            if all(item >= gene.maxfrag for item in gene.fragsize[detectedsite + 1 :]) and not all(
                item >= gene.maxfrag for item in gene.fragsize[: detectedsite - 1]
            ):
                shift = 3
                while gene.breaksites[detectedsite] + shift in gene.problemsites:
                    shift += 3
            if all(item >= gene.maxfrag for item in gene.fragsize[: detectedsite - 1]) and not all(
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
        gene.fragsize = [j - i for i, j in zip(gene.breaksites[:-1], gene.breaksites[1:])]
        # if DIMPLE.dms:
        #     tmpbreaklist = []
        #     for idx, x in enumerate(gene.breaksites[:-1]):
        #         if idx:
        #             tmpbreaklist.append([x, x + gene.fragsize[idx] + 3])
        #         else:
        #             tmpbreaklist.append([x + 3, x + gene.fragsize[idx] + 3])
        #     gene.breaklist = tmpbreaklist
        # else:
        gene.breaklist = [[x, x + gene.fragsize[idx]] for idx, x in enumerate(gene.breaksites[:-1])]
        # recheck for size limit issues
        tmpsite = [topidx for topidx, item in enumerate(gene.fragsize) if item > gene.maxfrag]
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


def _overlapping_count(haystack, needle) -> int:
    """Count occurrences of ``needle`` in ``haystack`` allowing overlaps.

    ``str.count`` is non-overlapping and therefore undercounts sites that overlap
    (e.g. ``"CGTCTCGTCTC".count("CGTCTC") == 1`` even though a Type IIS enzyme cuts
    at both position 0 and position 5).
    """
    hay = str(haystack).upper()
    sub = str(needle).upper()
    if not sub:
        return 0
    return sum(1 for i in range(len(hay) - len(sub) + 1) if hay[i : i + len(sub)] == sub)


def _overhang_recreates_cutsite(cutsite, cutsite_buffer, overhang) -> bool:
    """True if an overhang reconstitutes the enzyme recognition site.

    The oligo lays down the synthetic 5'->3' junction ``cutsite + buffer + overhang``.
    If the overhang (or its reverse complement) makes that junction contain the
    recognition site anywhere beyond the single intended copy at its 5' end -- on
    either strand -- the enzyme gains a spurious internal cut site and the fragment
    assembles with the wrong overhang (or not at all). Example: overhang ``TCTC``
    with BsmBI ``CGTCTC`` + ``G`` buffer -> ``CGTCTCGTCTC``, which contains ``CGTCTC``
    twice. The coding-side ``avoid_sequence`` screen never sees this because the
    spurious site spans the scaffold/overhang boundary, not the coding alone.
    """
    cs = str(cutsite).upper()
    csr = str(cutsite.reverse_complement()).upper()
    for oh in (overhang, overhang.reverse_complement()):
        junction = cs + str(cutsite_buffer).upper() + str(oh).upper()
        if _overlapping_count(junction, cs) > 1 or _overlapping_count(junction, csr) > 0:
            return True
    return False


def _internal_cutsite(oligo_seq, coding_core, cutsite) -> bool:
    """True if a Type IIS recognition site falls INSIDE the oligo's coding region.

    A well-formed oligo has its two designed sites flanking the coding; a site
    *inside* the coding span means the enzyme cuts the insert internally, so the
    fragment cannot assemble. A spurious site elsewhere (in a barcode arm) is
    benign -- it only trims a discarded piece. This distinguishes genuinely broken
    oligos from benign extra sites without a full assembly simulation.

    ``coding_core`` is the coding embedded in the oligo (mutation carrier, without
    the flanking overhangs), used only to locate the coding span in ``oligo_seq``.
    """
    oseq = str(oligo_seq).upper()
    core = str(coding_core).upper()
    start = oseq.find(core)
    if start < 0 or not core:
        return False
    end = start + len(core)
    for site in (str(cutsite).upper(), str(cutsite.reverse_complement()).upper()):
        j = oseq.find(site, start)
        if 0 <= j < end:
            return True
    return False


def check_overhangs(gene, pool, overlap_l, overlap_r):
    """TODO:
    Docstring
    """
    # Force all overhangs to be different within a gene (no more than 2 matching in a row)
    switched = False
    if not isinstance(gene, DIMPLE):
        raise TypeError("Not an instance of the DIMPLE class")
    cfg = gene.pool.config
    while True:
        detectedsites = set()  # stores matching overhangs
        for idx, y in enumerate(gene.breaklist):
            overhang_F = gene.seq[
                y[0] - gene.pool.config.cutsite_overhang - overlap_l : y[0] - overlap_r
            ]  # Forward overhang
            overhang_R = gene.seq[
                y[1] + overlap_l : y[1] + gene.pool.config.cutsite_overhang + overlap_r
            ]  # Reverse overhang
            if overhang_F == overhang_R or overhang_F == overhang_R.reverse_complement():
                detectedsites.update([idx])
            # Reject overhangs that reconstitute the Type IIS site through the
            # cutsite buffer (creates a spurious internal cut -> dead fragment).
            if (
                cfg.cutsite is not None
                and cfg.cutsite_buffer is not None
                and (
                    _overhang_recreates_cutsite(cfg.cutsite, cfg.cutsite_buffer, overhang_F)
                    or _overhang_recreates_cutsite(cfg.cutsite, cfg.cutsite_buffer, overhang_R)
                )
            ):
                detectedsites.update([idx])
        for detectedsite in detectedsites:
            switched = True
            if detectedsite == 0:
                detectedsite = 1  # don't mess with the first cut site
            print(
                "------------------ Fragment size swapped due to matching overhangs "
                "------------------"
            )
            switch_fragmentsize(gene, detectedsite, pool)
        else:  # if no detected sites
            break
    return switched
