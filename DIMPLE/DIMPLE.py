"""
DIMPLE: Deep Indel Missense Programmable Library Engineering

Python 3.7 package for generating oligo fragments and respective primers for
scanning Indel/Missense mutations

Written By: David Nedrud

Requires installation of Biopython
Simple installation command: pip install biopython

File input must be .fasta/.fa format and must include the whole plasmid for
primer specificity and binding
File output will also be .fasta format

Genes with variable sections can be aligned to save library space (avoid
synthesizing the same sequence multiple times)
Use align_genevariation()

"""

import csv
import itertools
import logging
import os
import re
import warnings
from difflib import SequenceMatcher
from math import ceil

from Bio import Align, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import seq1, seq3

from DIMPLE.core import DIMPLE
from DIMPLE.pool import DimpleRuntimeConfig, Pool, addgene

logger = logging.getLogger(__name__)

# Public API. `DIMPLE/DIMPLE.py` re-exports module-split helpers so existing
# `from DIMPLE.DIMPLE import ...` callers keep working. Listed here so ruff F401
# doesn't strip them as unused.
__all__ = [
    "DIMPLE",
    "addgene",
    "Pool",
    "DimpleRuntimeConfig",
    "align_genevariation",
    "generate_DMS_fragments",
    "check_overhangs",
    "recalculate_num_fragments",
    "switch_fragmentsize",
    "combine_fragments",
    "print_all",
    "check_nonspecific",
    "find_fragment_primer",
    "find_geneprimer",
    "post_qc",
    "check_final_assembly",
]


# This function is not used in the current version of the code
# This will find genes that share the same sequence and avoid synthesizing the same oligos
# multiple times
def align_genevariation(pool):
    if not isinstance(pool[0], DIMPLE):
        raise TypeError("Not an instance of the DIMPLE class")

    def should_link_genes() -> bool:
        if pool.config.link_policy == "always":
            return True
        if pool.config.link_policy == "never" or pool.config.non_interactive:
            return False
        return input("Are these genes linked? (y/n):") == "y"

    match = []
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    print("------------Finding homologous regions------------")
    # First find genes with matching sequences
    for m in range(len(pool)):
        remlist = range(len(pool))[m + 1 :]
        for p in remlist:
            alignments = aligner.align(pool[m].seq, pool[p].seq)
            for alignment in alignments:
                score = alignment.score
                score_len = alignment.indices.shape[1]
                break
            if score / score_len > 1.5:  # Threshold for a matched gene set
                index = [
                    x for x, geneset in enumerate(match) if m in geneset or p in geneset
                ]  # Determine if aligned genes are in any of the previously matched sets
                if not index:  # Create a new set if not
                    print(pool[m].geneid)
                    print(pool[p].geneid)
                    if should_link_genes():
                        match.append(set([m, p]))
                else:
                    if p not in match[index[0]] or m not in match[index[0]]:
                        for items in match[index[0]].union(set([p, m])):
                            print(pool[items].geneid)
                        if should_link_genes():
                            match[index[0]].add(p)
                            match[index[0]].add(m)
    # Create fragments for each match
    if match:
        for tmpset in match:
            matchset = list(tmpset)
            print(
                "Determining Gene Variation for genes:"
                + ",".join([pool[i].geneid for i in matchset])
            )
            max_gene_len = 0
            variablesites = set()
            for i, j in itertools.combinations(matchset, 2):
                max_gene_len = max(
                    max_gene_len,
                    len(pool[i].seq) - 2 * pool.config.primer_buffer,
                    len(pool[j].seq) - 2 * pool.config.primer_buffer,
                )
                seq_match = SequenceMatcher(None, pool[i].seq, pool[j].seq)
                # Determine variable regions
                variablesites.update(
                    [
                        x.size
                        for x in seq_match.get_matching_blocks()
                        if x.size != len(pool[i].seq) and x.size != len(pool[j].seq) and x.size != 0
                    ]
                )  # not sure how to account for zero
            problemsites = set()
            for kk in variablesites:
                problemsites.update(
                    range(kk - pool.config.primer_buffer, kk + pool.config.primer_buffer)
                )  # Add space for primers to bind
            # Determine Fragment Size while avoiding variable regions - must be same for all genes
            num = int(
                round(((max_gene_len) / float(pool.config.maxfrag)) + 0.499999999)
            )  # total bins needed (rounded up)
            insertionsites = range(
                pool.config.primer_buffer,
                max_gene_len + pool.config.primer_buffer - 6,
                3,
            )  # all genes start with a buffer
            fragsize = [len(insertionsites[i::num]) * 3 for i in list(range(num))]
            total = pool.config.primer_buffer
            breaksites = [
                pool.config.primer_buffer
            ]  # first site is always the max primer length (adjusted at beginning)
            for x in fragsize:
                total += x
                breaksites.extend([total])
            available_sites = [
                xsite
                for xsite in range(0, max_gene_len + pool.config.primer_buffer + 1, 3)
                if xsite not in problemsites
            ]
            breaksites = [
                (
                    site
                    if site in available_sites
                    else min(available_sites, key=lambda x: abs(x - site))
                )
                for site in breaksites
            ]  # remove problemsites?
            if any(x < DIMPLE.minfrag or x > pool.config.maxfrag for x in fragsize):
                print(fragsize)
                raise ValueError(
                    "Fragment size too low"
                )  # this was decided by author. could be changed
            fragsize = [j - i for i, j in zip(breaksites[:-1], breaksites[1:])]
            breaklist = [
                [x, x + fragsize[idx]] for idx, x in enumerate(breaksites[:-1])
            ]  # insertion site to insertion site
            unique_Frag = [
                [] for x in range(max(matchset) + 1)
            ]  # a list of fragments that do not match
            for x in breaklist:
                sequences = [str(pool[i].seq[x[0] : x[1]]) for i in matchset]
                index = [matchset[i] for i, x in enumerate(sequences) if i == sequences.index(x)]
                for idx in matchset:
                    if idx in index:
                        unique_Frag[idx].extend([True])
                    else:
                        unique_Frag[idx].extend([False])

            print("Finished Alignment. Fragment Sizes for combined genes:")
            print(fragsize)
            for (
                idx
            ) in (
                matchset
            ):  # setting these to the same variable should link them for processing later
                pool[idx].problemsites = (
                    # add gap range to problemsites variable to avoid breaking in a gap
                    problemsites
                )
                pool[idx].breaklist = breaklist
                pool[idx].fragsize = fragsize
                pool[idx].breaksites = breaksites
                pool[idx].linked.update(matchset)
                pool[idx].unique_Frag = unique_Frag[idx]
    else:
        print(
            "No redundant sequences found. Matching sequences may be too short or not aligned "
            "to reduce number of oligos synthesized"
        )


def generate_DMS_fragments(
    pool,
    overlapL,
    overlapR,
    synonymous,
    custom_mutations,
    dms=True,
    insert=False,
    delete=False,
    dis=False,
    folder="",
    config=None,
):
    """Generates the mutagenic oligos and writes the output to files."""

    # For each variant, also add an entry for the designed variants
    # sheet used in Dumpling. This is a list of dicts, where each dict
    # corresponds to a variant and contains the following keys:
    # - 'count': the number of reads observed for the variant (0 here)
    # - 'pos': the (codon) position of the variant
    # - 'mutation_type': the class of mutation (e.g. 'M', 'S', 'I', 'D', 'X')
    # - 'name': the simple name of the variant (e.g. 'M1A', 'E15del', 'N7_Y9del', 'N28_G29insGSG')
    # - 'codon': the variant codon sequence
    # - 'wt_codon': the wt codon sequence
    # - 'mutation': the specific type of mutation (e.g. 'D_1', 'I_3', 'R', 'A')
    # - 'length': the number of codons changed by the variant
    # - 'hgvs': the hgvs notation for the variant
    # - 'sequence': the full sequence of the variant

    # dms set to true for subsitition mutations
    # insert set to a list of insertions
    # delete set to a list of numbers of symmetrical deletions
    if not isinstance(pool[0], DIMPLE):
        raise TypeError("Not an instance of the DIMPLE class")
    # Loop through each gene or gene variation
    finishedGenes = []
    # Adjust fragments to account for variable sized fragments with the same subpool
    # barcodes/primers
    if insert or delete or dis:
        insert_list = []
        if insert:
            insert_list.extend(insert)
        if dis:
            insert_list.append(pool.config.handle)
        if insert or dis:
            pool.config.maxfrag = (
                pool.config.synth_len
                - DIMPLE.maxfrag_offset
                - max([len(x) for x in insert_list])
                - overlapL
                - overlapR
            )  # increase barcode space to allow for variable sized fragments within an oligo
        if delete and not insert and not dis:
            pool.config.maxfrag = (
                pool.config.synth_len - DIMPLE.maxfrag_offset - overlapL - overlapR
            )
        print("New max fragment:" + str(pool.config.maxfrag))
        for gene in pool:
            gene.maxfrag = pool.config.maxfrag
            switch_fragmentsize(gene, 1, pool)

    # Generate oligos for each gene
    for ii, gene in enumerate(pool):
        print(gene.breaklist)
        print(
            "--------------------------------- Analyzing Gene:"
            + gene.geneid
            + " ---------------------------------"
        )
        # gene.breaklist[0][0] += 0  # Do not mutate first codon
        # gene.fragsize[0] += -3  # Adjust size to match breaklist
        gene.maxfrag = pool.config.maxfrag
        if not any(
            [tmp in finishedGenes for tmp in gene.linked]
        ):  # only run analysis for one of the linked genes
            # Quality Control for overhangs from the same gene
            check_overhangs(gene, pool, overlapL, overlapR)
        # Generate oligos and Primers
        idx = 0  # index for fragment
        # storage for unused barcodes
        compileF = []
        compileR = []
        missingSites = []
        # missingTable = [[1]*gene.aacount]*gene.aacount
        missingFragments = []
        all_grouped_oligos = []
        # Loop through each fragment
        while idx < len(gene.breaklist):
            if idx == 0:
                gene.oligos = []
                gene.barPrimer = []
                gene.genePrimer = []
            frag = gene.breaklist[idx]
            grouped_oligos = []
            # AA range for fragment (need to subtract beginning primer buffer)
            fragstart = str(int((frag[0] - pool.config.primer_buffer) / 3) + 2)
            fragend = str(int((frag[1] - pool.config.primer_buffer) / 3) + 1)
            print(
                "Creating Fragment:"
                + gene.geneid
                + " --- Fragment #"
                + str(idx + 1)
                + " AA:"
                + fragstart
                + "-"
                + fragend
            )
            # only run analysis for one of the linked genes
            if not any([tmp in finishedGenes for tmp in gene.linked]):
                # Primers for gene amplification with addition of restriction enzyme site
                genefrag_R = gene.seq[
                    frag[0] - pool.config.primer_buffer : frag[0] + pool.config.primer_buffer
                ]
                reverse, tmR, sR = find_geneprimer(
                    genefrag_R, 15, pool.config.primer_buffer + 1 - overlapL, pool
                )  # 15 is just a starting point
                genefrag_F = gene.seq[
                    frag[1] - pool.config.primer_buffer : frag[1] + pool.config.primer_buffer
                ]
                forward, tmF, sF = find_geneprimer(
                    genefrag_F.reverse_complement(),
                    15,
                    pool.config.primer_buffer + 1 - overlapR,
                    pool,
                )
                # negative numbers look for reverse primers
                # 10 bases is the buffer overhang on the primer (ATA + (N))
                tmpr = check_nonspecific(
                    reverse,
                    gene.seq,
                    frag[0]
                    - len(gene.seq)
                    + 3
                    + len(pool.config.cutsite_buffer)
                    + len(pool.config.cutsite)
                    - overlapL,
                )
                tmpf = check_nonspecific(
                    forward,
                    gene.seq,
                    frag[1]
                    - 3
                    - len(pool.config.cutsite_buffer)
                    - len(pool.config.cutsite)
                    + overlapR,
                )
                if tmpf or tmpr:
                    # swap size with another fragment
                    print(
                        "------------------ Fragment size swapped due to non-specific primers "
                        "------------------"
                    )
                    if tmpf:
                        idx = idx + 1
                        print("Non specific primer F: " + forward)
                    else:
                        print("Non specific primer R: " + reverse)
                    # swap size with another fragment
                    skip = switch_fragmentsize(gene, idx, pool)
                    if skip:
                        # if end of gene, try to extend primer to make it more specific?
                        if tmpr:
                            reverse += gene.complement[genefrag_R[sR - 1]]
                            warnings.warn(
                                "Gene primer at the end of gene has non specific annealing. "
                                "Please Check this primer manually: " + str(reverse)
                            )
                            logger.warning(
                                "Gene primer at the end of gene has non specific annealing. "
                                "Please Check this primer manually: " + str(reverse)
                            )
                        if tmpf:
                            idx -= 1
                            forward += Seq(
                                genefrag_F.reverse_complement()[sF - 1]
                            ).reverse_complement()
                            warnings.warn(
                                "Gene primer at the end of gene has non specific annealing. "
                                "Please Check this primer manually: " + str(forward)
                            )
                            logger.warning(
                                "Gene primer at the end of gene has non specific annealing. "
                                "Please Check this primer manually: " + str(forward)
                            )
                    else:
                        # Quality Control for overhangs from the same gene
                        # check_overhangs(gene, pool)
                        pool.config.barcode_f.extend(compileF)  # return unused barcodes
                        pool.config.barcode_r.extend(compileR)
                        compileF = []  # reset unused primers
                        compileR = []
                        gene.genePrimer = []  # reset gene all primers due to nonspecific primer
                        gene.barPrimer = []
                        idx = 0
                        continue  # return to the beginning
                elif check_overhangs(gene, pool, overlapL, overlapR):
                    pool.config.barcode_f.extend(compileF)  # return unused barcodes
                    pool.config.barcode_r.extend(compileR)
                    compileF = []  # reset unused primers
                    compileR = []
                    gene.genePrimer = []  # reset gene all primers due to nonspecific primer
                    gene.barPrimer = []
                    idx = 0
                    continue  # return to the beginning
                # Store
                gene.genePrimer.append(
                    SeqRecord(
                        reverse,
                        id=gene.geneid + "_geneP_Mut-" + str(idx + 1) + "_R",
                        description="Frag" + fragstart + "-" + fragend + " " + str(tmR) + "C",
                    )
                )
                gene.genePrimer.append(
                    SeqRecord(
                        forward,
                        id=gene.geneid + "_geneP_Mut-" + str(idx + 1) + "_F",
                        description="Frag" + fragstart + "-" + fragend + " " + str(tmF) + "C",
                    )
                )
                # Determine missing double mutations
                beginning = int(
                    (frag[0] - pool.config.primer_buffer - sR) / 3
                )  # Region missing double mutations
                if beginning < 1:
                    beginning = 1
                end = ceil((frag[1] - pool.config.primer_buffer + sF) / 3)
                if end > ceil((gene.breaksites[-1] - pool.config.primer_buffer) / 3):
                    end = ceil((gene.breaksites[-1] - pool.config.primer_buffer) / 3)
                missingTmp = set()
                for site in range(beginning, end):  # Record these missing double mutations
                    for site2 in range(site + 1, end):
                        missingSites.append([site, site2])
                        missingTmp.add(site)
                        missingTmp.add(site2)
                        # missingTable[site][site2] = 0
                missingFragments.append([(frag[0] - 30) / 3, (frag[1] - 30) / 3, list(missingTmp)])
            if gene.unique_Frag[idx]:  # only for unique sequences
                # Create gene fragments with insertions
                count = 0
                tmpseq = gene.seq[
                    frag[0]
                    - pool.config.cutsite_overhang
                    - overlapL : frag[1]
                    + pool.config.cutsite_overhang
                    + overlapR
                ].replace(
                    "-", ""
                )  # extract sequence for oligo fragment include an extra 4 bases for BsmBI
                #    cut site and overlap
                offset = pool.config.cutsite_overhang + overlapL
                ## Create the mutations
                dms_sequences = []
                dms_sequences_double = []
                # list positions to mutate
                if custom_mutations:
                    # find custom mutations in the fragment range
                    tmp_positions = list(custom_mutations.keys())
                    tmp_tmp_positions = [
                        x * 3 - 3 + pool.config.primer_buffer for x in list(custom_mutations.keys())
                    ]
                    tmp_mut_positions = [
                        [i, x + offset - frag[0]]
                        for i, x in enumerate(tmp_tmp_positions)
                        if frag[0] <= x + 3 <= frag[1] - 3
                    ]
                    mut_positions = [x for i, x in tmp_mut_positions]
                    positions = [tmp_positions[i] for i, x in tmp_mut_positions]
                else:
                    mut_positions = range(offset, offset + frag[1] - frag[0], 3)
                    positions = [
                        int((frag[0] + x + 3 - offset - pool.config.primer_buffer) / 3)
                        for x in mut_positions
                    ]
                ### Deep Mutational Scanning
                if dms:
                    mutations = {}
                    for i in mut_positions:
                        wt_codon = tmpseq[i : i + 3].upper()
                        wt = [
                            name
                            for name, codon in gene.SynonymousCodons.items()
                            if wt_codon in codon
                        ]
                        if custom_mutations:
                            mutations_to_make = [
                                seq3(x)
                                for x in custom_mutations[positions[mut_positions.index(i)]].split(
                                    ","
                                )
                            ]
                        else:
                            mutations_to_make = gene.aminoacids
                        for jk in (x for x in mutations_to_make):
                            # check if synonymous and if user wants these mutations
                            if jk not in wt[0] or synonymous:
                                if jk == wt[0]:
                                    is_synonymous = True
                                elif jk == "STOP":
                                    is_stop = True
                                else:
                                    is_stop = False
                                    is_synonymous = False
                                codons = [
                                    aa for aa in gene.SynonymousCodons[jk] if aa not in wt_codon
                                ]
                                p = [
                                    pool.config.usage[aa] for aa in codons
                                ]  # Find probabilities but not wild type codon
                                p = [
                                    xp if xp > 0.1 else 0 for xp in p
                                ]  # Remove probabilities below 0.1
                                p = [xp / sum(p) for xp in p]  # Normalize to 1
                                if not p:
                                    continue
                                # if the user wants to maximize the number of nucleotide changes
                                synonymous_mutation = []
                                synonymous_position = 0
                                if pool.config.maximize_nucleotide_change:
                                    # remove codons with only one change compared to wt_codon
                                    max_codons = [
                                        x
                                        for x in codons
                                        if sum([x[i] != wt_codon[i] for i in range(3)]) > 1
                                    ]
                                    if max_codons:
                                        # if there are codons with more than one base change
                                        mutation = gene.rng.choice(
                                            max_codons, 1, p
                                        )  # Pick one codon
                                        xfrag = (
                                            tmpseq[0:i] + mutation[0] + tmpseq[i + 3 :]
                                        )  # Add mutation to fragment
                                    else:
                                        # no codons with more than one base change.
                                        # Creating synonymous mutation in neighboring codon.
                                        mutation = gene.rng.choice(codons, 1, p)  # Pick one codon
                                        # find neighboring codon
                                        tmp_synonymous = [
                                            name
                                            for name, codon in gene.SynonymousCodons.items()
                                            if tmpseq[i - 3 : i] in codon
                                        ]
                                        synonymous_codons = gene.SynonymousCodons[tmp_synonymous[0]]
                                        max_synonymous = [
                                            x
                                            for x in synonymous_codons
                                            if sum([x[c] != tmpseq[i - 3 : i][c] for c in range(3)])
                                            > 0
                                        ]
                                        if max_synonymous and not (
                                            idx == 0 and mut_positions.index(i) == 0
                                        ):
                                            synonymous_mutation = gene.rng.choice(max_synonymous, 1)
                                            xfrag = (
                                                tmpseq[0 : i - 3]
                                                + synonymous_mutation[0]
                                                + mutation[0]
                                                + tmpseq[i + 3 :]
                                            )  # Add mutation to fragment
                                            synonymous_position = -1
                                        else:
                                            tmp_synonymous = [
                                                name
                                                for name, codon in gene.SynonymousCodons.items()
                                                if tmpseq[i + 3 : i + 6] in codon
                                            ]
                                            synonymous_codons = gene.SynonymousCodons[
                                                tmp_synonymous[0]
                                            ]
                                            max_synonymous = [
                                                x
                                                for x in synonymous_codons
                                                if sum(
                                                    [
                                                        x[c] != tmpseq[i + 3 : i + 6][c]
                                                        for c in range(3)
                                                    ]
                                                )
                                                > 0
                                            ]
                                            if max_synonymous:
                                                synonymous_mutation = gene.rng.choice(
                                                    max_synonymous, 1
                                                )
                                                xfrag = (
                                                    tmpseq[0:i]
                                                    + mutation[0]
                                                    + synonymous_mutation[0]
                                                    + tmpseq[i + 6 :]
                                                )  # Add mutation to fragment
                                                synonymous_position = +1
                                            else:
                                                print(
                                                    "Unable to create synonymous mutation in "
                                                    "neighboring codon. Continuing with single "
                                                    "nucleotide change"
                                                )
                                                xfrag = tmpseq[0:i] + mutation[0] + tmpseq[i + 3 :]
                                                print(xfrag)
                                else:
                                    mutation = gene.rng.choice(codons, 1, p)  # Pick one codon
                                    xfrag = (
                                        tmpseq[0:i] + mutation[0] + tmpseq[i + 3 :]
                                    )  # Add mutation to fragment
                                # Check each cassette for more than 2 BsmBI and 2 BsaI sites
                                avoid_count = 0
                                while any(
                                    [
                                        (
                                            xfrag.upper().count(x)
                                            + xfrag.upper().count(x.reverse_complement())
                                        )
                                        > 0
                                        for x in pool.config.avoid_sequence
                                    ]
                                ):
                                    mutation = gene.rng.choice(
                                        gene.SynonymousCodons[jk], 1, p
                                    )  # Pick one codon
                                    avoid_count += 1
                                    xfrag = tmpseq[0:i] + mutation[0] + tmpseq[i + 3 :]
                                    if avoid_count > 10:
                                        warnings.warn(
                                            "Unwanted restriction site found within "
                                            f"substitution fragment: {str(xfrag)}"
                                        )
                                        logger.error(
                                            "Unwanted restriction site found within "
                                            f"substitution fragment: {str(xfrag)}"
                                        )
                                        break
                                mutations[
                                    ">"
                                    + wt[0]
                                    + str(
                                        int(
                                            (frag[0] + i + 6 - offset - pool.config.primer_buffer)
                                            / 3
                                        )
                                    )
                                    + jk
                                ] = mutation[0]
                                # if there was a synonymous mutation added then add the
                                # synonymous mutation to the mutation list
                                if synonymous_mutation:
                                    mutations[
                                        ">"
                                        + wt[0]
                                        + str(
                                            int(
                                                (
                                                    frag[0]
                                                    + i
                                                    + 6
                                                    - offset
                                                    - pool.config.primer_buffer
                                                )
                                                / 3
                                            )
                                        )
                                        + jk
                                    ] += (str(synonymous_position) + "_" + synonymous_mutation[0])
                                oligo_id = (
                                    gene.geneid
                                    + "_DMS-"
                                    + str(idx + 1)
                                    + "_"
                                    + wt[0]
                                    + str(
                                        int(
                                            (frag[0] + i + 6 - offset - pool.config.primer_buffer)
                                            / 3
                                        )
                                    )
                                    + jk
                                )
                                dms_sequences.append(
                                    SeqRecord(
                                        xfrag,
                                        id=oligo_id,
                                        description="Frag " + fragstart + "-" + fragend,
                                    )
                                )
                                if is_synonymous:
                                    mutation_type = "S"
                                elif is_stop:
                                    mutation_type = "X"
                                else:
                                    mutation_type = "M"
                                aa_pos = int(
                                    (frag[0] + i + 6 - offset - pool.config.primer_buffer) / 3
                                )
                                name = f"{seq1(wt[0])}{aa_pos}{seq1(jk)}"
                                gene.designed_variants[oligo_id] = {
                                    "count": 0,
                                    "pos": int(
                                        (frag[0] + i + 6 - offset - pool.config.primer_buffer) / 3
                                    ),
                                    "mutation_type": mutation_type,
                                    "name": name,
                                    "codon": mutation[0],
                                    "wt_codon": wt_codon,
                                    "mutation": seq1(jk),
                                    "length": 1,
                                    "hgvs": f"p.({name})",
                                    "fragment": idx + 1,
                                    "xfrag": xfrag,
                                }
                        # if double mutations are selected then make every possible double mutation
                        if pool.config.make_double:
                            # select every permutation of mut_positions order doesn't matter
                            for combi in itertools.combinations(mutations.keys(), 2):
                                # extract number from mutation name
                                if "STOP" not in combi[0] and "STOP" not in combi[1]:
                                    pos1 = mut_positions[
                                        positions.index(int(re.findall(r"\d+", combi[0])[0]))
                                    ]
                                    pos2 = mut_positions[
                                        positions.index(int(re.findall(r"\d+", combi[1])[0]))
                                    ]
                                    if pos1 != pos2:
                                        xfrag = (
                                            tmpseq[0:pos1]
                                            + mutations[combi[0]]
                                            + tmpseq[pos1 + 3 : pos2]
                                            + mutations[combi[1]]
                                            + tmpseq[pos2 + 3 :]
                                        )
                                        dms_sequences_double.append(
                                            SeqRecord(
                                                xfrag,
                                                id=gene.geneid
                                                + "_DMS-"
                                                + str(idx + 1)
                                                + "_"
                                                + combi[0].strip(">")
                                                + "+"
                                                + combi[1].strip(">"),
                                                description="Frag " + fragstart + "-" + fragend,
                                            )
                                        )
                    # record mutation for analysis with NGS
                    # TODO: Don't append.
                    with open(
                        os.path.join(folder.replace("\\", ""), gene.geneid + "_mutations.csv"),
                        "a",
                    ) as file:
                        for mut in mutations.keys():
                            file.write(mut + "\n")
                            file.write(mutations[mut] + "\n")
                ### Scanning Insertions
                if insert:
                    insert_translations = {}
                    for insertion_sequence in insert:
                        if len(insertion_sequence) % 3 == 0:
                            insert_translations[insertion_sequence] = Seq(
                                insertion_sequence
                            ).translate()
                        else:
                            logger.warning(
                                f"Insertion sequence {insertion_sequence} is not a multiple "
                                "of 3. Will not translate in output."
                            )
                            insert_translations[insertion_sequence] = f"({insertion_sequence})"

                    # insertion

                    for i in range(offset, offset + frag[1] - frag[0], 3):
                        pos = int((frag[0] + i + 3 - offset - pool.config.primer_buffer) / 3)
                        wt_pre_codon = tmpseq[i : i + 3].upper()
                        wt_post_codon = tmpseq[i + 3 : i + 6].upper()
                        wt_pre_aa = [
                            name
                            for name, codon in gene.SynonymousCodons.items()
                            if wt_pre_codon in codon
                        ]
                        wt_post_aa = [
                            name
                            for name, codon in gene.SynonymousCodons.items()
                            if wt_post_codon in codon
                        ]

                        for insert_n in insert:
                            xfrag = tmpseq[0:i] + insert_n + tmpseq[i:]  # Add mutation to fragment
                            # Check each cassette for more than 2 BsmBI and 2 BsaI sites
                            while any(
                                [
                                    (
                                        xfrag.upper().count(x)
                                        + xfrag.upper().count(x.reverse_complement())
                                    )
                                    > 0
                                    for x in pool.config.avoid_sequence
                                ]
                            ):
                                warnings.warn(
                                    "Unwanted restriction site found within insertion fragment: "
                                    + str(xfrag)
                                )
                                logger.warning(
                                    "Unwanted restriction site found within insertion fragment: "
                                    + str(xfrag)
                                )
                                break
                                # not sure how to solve this issue
                                # mutation?
                                # xfrag = tmpseq[0:i] + mutation + tmpseq[i + 3:]
                            oligo_id = (
                                gene.geneid
                                + "_insert-"
                                + str(idx + 1)
                                + "_"
                                + insert_n
                                + "-"
                                + str(pos)
                            )
                            dms_sequences.append(
                                SeqRecord(
                                    xfrag,
                                    id=oligo_id,
                                    description="Frag " + fragstart + "-" + fragend,
                                )
                            )
                            # Translate insert_n
                            insert_name = insert_translations[insert_n]
                            name = (
                                f"{seq1(wt_pre_aa)}{pos}_{seq1(wt_post_aa)}{pos+1}_ins{insert_name}"
                            )
                            # TODO: Insert length assumes that the insert is a multiple of 3
                            # (i.e. codon insertions). Make more flexible.
                            gene.designed_variants[oligo_id] = {
                                "count": 0,
                                "pos": pos,
                                "mutation_type": "I",
                                "name": name,
                                "codon": insert_n,
                                "wt_codon": "",
                                "mutation": f"I_{len(insert_n)//3}",
                                "length": f"{len(insert_n)//3}",
                                "hgvs": f"p.({name})",
                                "fragment": idx + 1,
                                "xfrag": xfrag,
                            }
                ### Scanning Deletions
                if delete:
                    # deletion
                    # TODO: failing here, for some reason. i becomes too large.
                    # fragment lengths are too high? no.
                    # overlaps are too small for larger deletion sizes. why?

                    # Iterate over codon boundaries in the fragment
                    # Shifted down by 3 to avoid long deletions running into the primer
                    # binding region
                    for i in range(offset - 3, offset + frag[1] - frag[0] - 3, 3):
                        # Calculate the amino acid position too
                        pos = int((frag[0] + i + 6 - offset - pool.config.primer_buffer) / 3)
                        # List of wt codons for each position in the range of deletion lengths
                        wt_codons = [
                            tmpseq[i + j : i + j + 3].upper() for j in range(0, max(delete), 3)
                        ]
                        wt_aas = [
                            [
                                name
                                for name, codon in gene.SynonymousCodons.items()
                                if wt_codon in codon
                            ]
                            for wt_codon in wt_codons
                        ]

                        for delete_n in delete:
                            # Check if deletion extends beyond ORF.
                            if pos + delete_n > len(gene.seq) / 3:
                                logger.warning(
                                    "Deletion extends beyond ORF: " + f"D{pos}_{delete_n}"
                                )
                                pass
                            # Check if deletion extends beyond the fragment.
                            if delete_n + i > len(tmpseq):
                                print("overlap: ", overlapL)
                                print("offset: ", offset)
                                print("frag: ", frag)
                                print("tmpseq: ", tmpseq)
                                print("delete_n: ", delete_n)
                                print("length: ", len(tmpseq))
                                print("max i: ", offset + frag[1] - frag[0] + 3)
                                print("i: ", i)
                                raise ValueError(
                                    "deletions cannot be larger than fragment itself: "
                                    "adjust settings and retry."
                                )
                            else:
                                xfrag = (
                                    tmpseq[0:i] + tmpseq[i + delete_n :]
                                )  # delete forward from position only

                            # Make sure that the 3' end has sufficient sequence to trim for cutsites
                            # Number of bases trimmed is pool.config.cutsite_overhang (usually 4)
                            # Add dummy bases to 3' end if not enough sequence.
                            if len(tmpseq[i + delete_n :]) < pool.config.cutsite_overhang:
                                pass
                                # buffer_length = (
                                #     pool.config.cutsite_overhang - len(tmpseq[i + delete_n :])
                                # )
                                # xfrag = xfrag + "N" * buffer_length

                            # Check each cassette for more than 2 BsmBI and 2 BsaI sites

                            while any(
                                [
                                    (
                                        xfrag.upper().count(x)
                                        + xfrag.upper().count(x.reverse_complement())
                                    )
                                    > 0
                                    for x in pool.config.avoid_sequence
                                ]
                            ):
                                warnings.warn(
                                    "Unwanted restriction site found within deletion fragment: "
                                    + str(xfrag)
                                )
                                logger.warning(
                                    "Unwanted restriction site found within deletion fragment: "
                                    + str(xfrag)
                                )
                                break
                                # xfrag = tmpseq[0:i-delete_n-3] + tmpseq[i+delete_n:]
                                # iteratively shift deletion to avoid cut sites? or mutate
                                # codons of near by aa?
                            oligo_id = (
                                gene.geneid
                                + "_delete-"
                                + str(idx + 1)
                                + "_"
                                + str(delete_n)
                                + "-"
                                + str(pos)
                            )
                            dms_sequences.append(
                                SeqRecord(
                                    xfrag,
                                    id=oligo_id,
                                    description="Frag " + fragstart + "-" + fragend,
                                )
                            )
                            length = int(delete_n / 3)
                            if length == 1:
                                name = f"{seq1(wt_aas[0][0])}{pos}del"
                            else:
                                first_aa = seq1(wt_aas[0][0])
                                last_aa = seq1(wt_aas[length - 1][0])
                                name = f"{first_aa}{pos}_{last_aa}{pos + length - 1}del"

                            gene.designed_variants[oligo_id] = {
                                "count": 0,
                                "pos": pos,
                                "mutation_type": "D",
                                "name": name,
                                "codon": "",
                                "wt_codon": tmpseq[i : i + delete_n].upper(),
                                "mutation": f"D_{length}",
                                "length": length,
                                "hgvs": f"p.({name})",
                                "fragment": idx + 1,
                                "xfrag": xfrag,
                            }
                ### Scanning Domain Insertions
                if dis:
                    # Translate the domain insertion handle for naming in the output
                    if len(pool.config.handle) % 3 == 0:
                        handle_name = Seq(str(pool.config.handle)).translate()
                    else:
                        logger.warning(
                            f"Domain insertion handle {pool.config.handle} is not a multiple "
                            "of 3. Will not translate in output."
                        )
                        handle_name = "(handle)"
                    # insertion
                    for i in range(offset, offset + frag[1] - frag[0], 3):
                        # if idx == 0:
                        #    continue
                        pos = int((frag[0] + i + 3 - offset - pool.config.primer_buffer) / 3)
                        wt_pre_codon = tmpseq[i : i + 3].upper()
                        wt_post_codon = tmpseq[i + 3 : i + 6].upper()
                        wt_pre_aa = [
                            name
                            for name, codon in gene.SynonymousCodons.items()
                            if wt_pre_codon in codon
                        ]
                        wt_post_aa = [
                            name
                            for name, codon in gene.SynonymousCodons.items()
                            if wt_post_codon in codon
                        ]
                        xfrag = (
                            tmpseq[0:i] + pool.config.handle + tmpseq[i:]
                        )  # Add mutation to fragment
                        # Check each cassette for more than 2 BsmBI and 2 BsaI sites
                        while any(
                            [
                                (
                                    xfrag.upper().count(x)
                                    + xfrag.upper().count(x.reverse_complement())
                                )
                                > 2
                                for x in pool.config.avoid_sequence
                            ]
                        ):
                            warnings.warn(
                                "Unwanted restriction site found within domain insertion fragment: "
                                + str(xfrag)
                            )
                            logger.warning(
                                "Unwanted restriction site found within domain insertion fragment: "
                                + str(xfrag)
                            )
                            break
                            # not sure how to solve this issue
                            # mutation?
                            # xfrag = tmpseq[0:i] + mutation + tmpseq[i + 3:]
                        oligo_id = gene.geneid + "_DIS-" + str(idx + 1) + "_" + str(pos)
                        dms_sequences.append(
                            SeqRecord(
                                xfrag,
                                id=oligo_id,
                                description="Frag " + fragstart + "-" + fragend,
                            )
                        )
                        name = f"{seq1(wt_pre_aa)}{pos}_{seq1(wt_post_aa)}{pos+1}_ins{handle_name}"
                        gene.designed_variants[oligo_id] = {
                            "count": 0,
                            "pos": pos,
                            "mutation_type": "DI",
                            "name": name,
                            "codon": str(pool.config.handle),
                            "wt_codon": "",
                            "mutation": f"DI_{len(pool.config.handle) // 3}",
                            "length": len(pool.config.handle) // 3,
                            "hgvs": f"p.({name})",
                            "fragment": idx + 1,
                            "xfrag": xfrag,
                        }
                for idx_type, dms_sequence_list in enumerate([dms_sequences, dms_sequences_double]):
                    if dms_sequence_list:  # are there any sequences to write?
                        tmF = 0
                        tmR = 0
                        if gene.num_frag_per_oligo > 1:
                            dms_sequence_list = combine_fragments(
                                dms_sequence_list,
                                gene.num_frag_per_oligo,
                                gene.split,
                                pool,
                            )
                        len_cutsite = (
                            len(pool.config.cutsite)
                            + len(pool.config.cutsite_buffer)
                            + pool.config.cutsite_overhang
                        )
                        # determine barcodes for subpool amplification based on smallest size
                        frag_sizes = [len(xf.seq) for xf in dms_sequence_list]
                        smallest_frag = dms_sequence_list[frag_sizes.index(min(frag_sizes))].seq
                        while (
                            tmF < DIMPLE.primerTm[0] or tmR < DIMPLE.primerTm[0]
                        ):  # swap out barcode if tm is low
                            difference = pool.config.synth_len - (
                                len(smallest_frag) + len_cutsite * 2
                            )  # 14 bases is the length of the restriction sites with overhangs
                            #    (7 bases each)
                            try:
                                barF = pool.config.barcode_f.pop(0)
                                barR = pool.config.barcode_r.pop(0)
                            except IndexError:
                                raise Exception("Ran out of barcodes.")
                            count += 1  # How many barcodes used
                            compileF.append(barF)
                            compileR.append(barR)
                            while (difference / 2) > len(barF):
                                tmpF = pool.config.barcode_f.pop(0)
                                tmpR = pool.config.barcode_r.pop(0)
                                compileF.append(tmpF)
                                compileR.append(tmpR)
                                barF += tmpF
                                barR += tmpR
                                count += 1  # How many barcodes used
                            tmpfrag_1 = (
                                barF.seq[0 : int(difference / 2)]
                                + pool.config.cutsite
                                + pool.config.cutsite_buffer
                                + tmpseq[0 : pool.config.cutsite_overhang]
                            )  # include recognition site and the 4 base overhang
                            tmpfrag_2 = (
                                tmpseq[-pool.config.cutsite_overhang :]
                                + pool.config.cutsite_buffer.reverse_complement()
                                + pool.config.cutsite.reverse_complement()
                                + barR.seq.reverse_complement()[
                                    0 : difference - int(difference / 2)
                                ]
                            )
                            # primers for amplifying subpools
                            offset = (
                                int(difference / 2) + len_cutsite
                            )  # add 11 bases for type 2 restriction
                            primerF, tmF = find_fragment_primer(tmpfrag_1, 25)
                            if len(primerF) > 21:
                                tmF = 0
                            primerR, tmR = find_fragment_primer(tmpfrag_2.reverse_complement(), 25)
                            if len(primerR) > 21:
                                tmR = 0
                        for (
                            sequence
                        ) in dms_sequence_list:  # add barcodes to the fragments to make the oligos
                            # Per-sequence trim path is required whenever the list can hold
                            # variable-length entries: insert/delete vary by definition; dis
                            # inserts a handle so its sequences are longer than concurrent
                            # DMS substitutions in the same list. The else-branch's blind
                            # concat only stays at synth_len when every sequence matches the
                            # smallest_frag the barcodes were sized for.
                            if insert or delete or dis:
                                cutsite_overhang = pool.config.cutsite_overhang
                                difference = (
                                    pool.config.synth_len
                                    - len(sequence.seq[cutsite_overhang:-cutsite_overhang])
                                    - len_cutsite * 2
                                )  # how many bases need to be added to make oligo correct length
                                offset = int(difference / 2)  # force it to be a integer

                                combined_sequence = (
                                    tmpfrag_1[:offset]
                                    + tmpfrag_1[-len_cutsite:]
                                    + sequence.seq[
                                        pool.config.cutsite_overhang : -pool.config.cutsite_overhang
                                    ]
                                    + tmpfrag_2[:len_cutsite]
                                    + tmpfrag_2[-(difference - offset) :]
                                )
                            else:
                                combined_sequence = (
                                    tmpfrag_1
                                    + sequence.seq[
                                        pool.config.cutsite_overhang : -pool.config.cutsite_overhang
                                    ]
                                    + tmpfrag_2
                                )
                            if (
                                primerF not in combined_sequence
                                or primerR.reverse_complement() not in combined_sequence
                            ):
                                print(primerF)
                                print(combined_sequence)
                                print("---")
                                print(combined_sequence.reverse_complement())
                                print(primerR)
                                logger.error(
                                    "Primers no longer bind to oligo. Was not able to add "
                                    "barcode to oligo. Try adjusting fragment length or "
                                    "synthesis length and try again."
                                )
                                raise Exception(
                                    "Primers no longer bind to oligo. Was not able to add "
                                    "barcode to oligo. Try adjusting fragment length or "
                                    "synthesis length and try again."
                                )
                            if (
                                combined_sequence.upper().count(pool.config.cutsite)
                                + combined_sequence.upper().count(
                                    pool.config.cutsite.reverse_complement()
                                )
                                < 2
                            ):
                                raise Exception("Oligo does not have 2 cutsites")
                            if len(combined_sequence) > pool.config.synth_len:
                                raise Exception(
                                    f"Oligo too long: {str(len(combined_sequence))} is longer "
                                    f"than {str(pool.config.synth_len)}"
                                )
                            if gene.doublefrag == 0:
                                gene.oligos.append(
                                    SeqRecord(
                                        combined_sequence,
                                        id=sequence.id,
                                        description="",
                                    )
                                )
                            else:
                                grouped_oligos.append(
                                    SeqRecord(
                                        combined_sequence,
                                        id=sequence.id,
                                        description="",
                                    )
                                )
                            gene.designed_variants[sequence.id][
                                "oligo_sequence"
                            ] = combined_sequence

                        # Store primers for gene fragment
                        if idx_type == 0:
                            gene.barPrimer.append(
                                SeqRecord(
                                    primerF,
                                    id=gene.geneid + "_oligoP_DMS-" + str(idx + 1) + "_F",
                                    description="Frag"
                                    + fragstart
                                    + "-"
                                    + fragend
                                    + "_"
                                    + str(tmF)
                                    + "C",
                                )
                            )
                            gene.barPrimer.append(
                                SeqRecord(
                                    primerR,
                                    id=gene.geneid + "_oligoP_DMS-" + str(idx + 1) + "_R",
                                    description="Frag"
                                    + fragstart
                                    + "-"
                                    + fragend
                                    + "_"
                                    + str(tmR)
                                    + "C",
                                )
                            )
                        else:
                            gene.barPrimer.append(
                                SeqRecord(
                                    primerF,
                                    id=gene.geneid + "_oligoP_DMS-double-" + str(idx + 1) + "_F",
                                    description="Frag"
                                    + fragstart
                                    + "-"
                                    + fragend
                                    + "_"
                                    + str(tmF)
                                    + "C",
                                )
                            )
                            gene.barPrimer.append(
                                SeqRecord(
                                    primerR,
                                    id=gene.geneid + "_oligoP_DMS-double-" + str(idx + 1) + "_R",
                                    description="Frag"
                                    + fragstart
                                    + "-"
                                    + fragend
                                    + "_"
                                    + str(tmR)
                                    + "C",
                                )
                            )
                        print("Barcodes tested:" + str(count))
                        # return unused barcodes
                        pool.config.barcode_f.extend(compileF[:-2])
                        pool.config.barcode_r.extend(compileR[:-2])
                        print("Barcodes Remaining:" + str(len(pool.config.barcode_f)))
                        compileF = []  # reset unused primers
                        compileR = []
            if gene.doublefrag == 1:
                all_grouped_oligos.append(grouped_oligos)
            idx += 1
        # Resolve Double Fragment
        if gene.doublefrag == 1:
            while len(all_grouped_oligos) > 1:
                listOne = all_grouped_oligos.pop(0)
                listTwo = all_grouped_oligos.pop(0)
                while listOne and listTwo:
                    one = listOne.pop(0)
                    two = listTwo.pop(0)
                    combined_sequence = one.seq + two.seq.reverse_complement()
                    combined_id = one.id + two.id
                    gene.oligos.append(SeqRecord(combined_sequence, id=combined_id, description=""))
                if listOne or listTwo:
                    if listOne:
                        sequence = listOne.pop(0)
                    if listTwo:
                        sequence = listTwo.pop(0)
                    combined_id = sequence.id
                    combined_sequence = sequence.seq
                    difference = 230 - len(combined_sequence)
                    # print(len(tmpseq))
                    barF2 = pool.config.barcode_f.pop(0)
                    barR2 = pool.config.barcode_r.pop(0)
                    while difference / 2 > len(barF2):
                        barF2 += pool.config.barcode_f.pop(0)
                        barR2 += pool.config.barcode_r.pop(0)
                    combined_sequence2 = (
                        barF2.seq[0 : int(difference / 2)]
                        + combined_sequence
                        + barR2.seq.reverse_complement()[0 : difference - int(difference / 2)]
                    )
                    gene.oligos.append(
                        SeqRecord(combined_sequence2, id=combined_id, description="")
                    )
            if all_grouped_oligos:
                one = all_grouped_oligos
                while one:
                    sequence_one = one.pop(0)
                    combined_id = sequence_one.id
                    combined_sequence = sequence_one.seq
                    difference = 230 - len(combined_sequence)
                    # print(len(tmpseq))
                    barF2 = pool.config.barcode_f.pop(0)
                    barR2 = pool.config.barcode_r.pop(0)
                    while difference / 2 > len(barF2):
                        barF2 += pool.config.barcode_f.pop(0)
                        barR2 += pool.config.barcode_r.pop(0)
                    combined_sequence2 = (
                        barF2.seq[0 : int(difference / 2)]
                        + combined_sequence
                        + barR2.seq.reverse_complement()[0 : difference - int(difference / 2)]
                    )
                    gene.oligos.append(
                        SeqRecord(combined_sequence2, id=combined_id, description="")
                    )
        # Export files (fasta)
        # Missing Mutation Pairs
        # import csv
        # missing2_path = os.path.join(folder.replace('\\', ''),
        #                              gene.geneid + "_missing2Mutations.csv")
        # with open(missing2_path, 'w') as csvfile:
        #    mutationwriter = csv.writer(csvfile, delimiter=',')
        #    mutationwriter.writerows(missingSites)
        #    mutationwriter.writerow('Fragment Info')
        #    mutationwriter.writerows(missingFragments)
        # Print table?
        # from tabulate import tabulate
        # print('Missing Double Mutation Table:')
        # print(tabulate(missingTable))
        # Fragments
        SeqIO.write(
            gene.oligos,
            os.path.join(folder.replace("\\", ""), gene.geneid + "_DMS_Oligos.fasta"),
            "fasta",
        )
        # Barcode Primers
        SeqIO.write(
            gene.barPrimer,
            os.path.join(folder.replace("\\", ""), gene.geneid + "_DMS_Oligo_Primers.fasta"),
            "fasta",
        )
        # Amplification Primers
        SeqIO.write(
            gene.genePrimer,
            os.path.join(folder.replace("\\", ""), gene.geneid + "_DMS_Gene_Primers.fasta"),
            "fasta",
        )

        # Designed Variants
        with open(
            os.path.join(folder.replace("\\", ""), gene.geneid + "_designed_variants.csv"),
            "w",
        ) as csvfile:
            fieldnames = [
                "count",
                "pos",
                "mutation_type",
                "name",
                "codon",
                "wt_codon",
                "mutation",
                "length",
                "hgvs",
            ]
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames, extrasaction="ignore")
            writer.writeheader()
            for variant in gene.designed_variants:
                writer.writerow(gene.designed_variants[variant])

        # Record finished gene for aligned genes
        finishedGenes.extend([ii])


from DIMPLE.fragment_layout import (  # noqa: E402
    check_overhangs,
    recalculate_num_fragments,
    switch_fragmentsize,
)
from DIMPLE.oligo_assembly import combine_fragments  # noqa: E402
from DIMPLE.outputs import print_all  # noqa: E402
from DIMPLE.primers import (  # noqa: E402
    check_nonspecific,
    find_fragment_primer,
    find_geneprimer,
)
from DIMPLE.qc import check_final_assembly, post_qc  # noqa: E402
