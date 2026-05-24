def parse_custom_mutations(mutation_text):
    custom_mutations = {}
    for set in mutation_text:
        set = set.split(":")
        if set[1] == "All":
            if "-" in set[0]:
                for i in range(int(set[0].split("-")[0]), int(set[0].split("-")[1]) + 1):
                    custom_mutations[i] = "A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y"
            else:
                custom_mutations[int(set[0])] = "A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y"
        else:
            if "-" in set[0]:
                for i in range(int(set[0].split("-")[0]), int(set[0].split("-")[1]) + 1):
                    custom_mutations[i] = set[1]
            else:
                # if mutation exists, add to it
                if int(set[0]) in custom_mutations.keys():
                    custom_mutations[int(set[0])] = custom_mutations[int(set[0])] + "," + set[1]
                else:
                    custom_mutations[int(set[0])] = set[1]
    return custom_mutations


def codon_usage(usage):
    if usage == "ecoli":
        usage_table = {
            "TTT": 0.58,
            "TTC": 0.42,
            "TTA": 0.14,
            "TTG": 0.13,
            "TAT": 0.59,
            "TAC": 0.41,
            "TAA": 0.61,
            "TAG": 0.09,
            "CTT": 0.12,
            "CTC": 0.1,
            "CTA": 0.04,
            "CTG": 0.47,
            "CAT": 0.57,
            "CAC": 0.43,
            "CAA": 0.34,
            "CAG": 0.66,
            "ATT": 0.49,
            "ATC": 0.39,
            "ATA": 0.11,
            "ATG": 1,
            "AAT": 0.49,
            "AAC": 0.51,
            "AAA": 0.74,
            "AAG": 0.26,
            "GTT": 0.28,
            "GTC": 0.2,
            "GTA": 0.17,
            "GTG": 0.35,
            "GAT": 0.63,
            "GAC": 0.37,
            "GAA": 0.68,
            "GAG": 0.32,
            "TCT": 0.17,
            "TCC": 0.15,
            "TCA": 0.14,
            "TCG": 0.14,
            "TGT": 0.46,
            "TGC": 0.54,
            "TGA": 0.3,
            "TGG": 1,
            "CCT": 0.18,
            "CCC": 0.13,
            "CCA": 0.2,
            "CCG": 0.49,
            "CGT": 0.36,
            "CGC": 0.36,
            "CGA": 0.07,
            "CGG": 0.11,
            "ACT": 0.19,
            "ACC": 0.4,
            "ACA": 0.17,
            "ACG": 0.25,
            "AGT": 0.16,
            "AGC": 0.25,
            "AGA": 0.07,
            "AGG": 0.04,
            "GCT": 0.18,
            "GCC": 0.26,
            "GCA": 0.23,
            "GCG": 0.33,
            "GGT": 0.35,
            "GGC": 0.37,
            "GGA": 0.13,
            "GGG": 0.15,
        }  # E.coli codon usage table
    elif usage == "human":
        usage_table = {
            "TTT": 0.45,
            "TTC": 0.55,
            "TTA": 0.07,
            "TTG": 0.13,
            "TAT": 0.43,
            "TAC": 0.57,
            "TAA": 0.28,
            "TAG": 0.2,
            "CTT": 0.13,
            "CTC": 0.2,
            "CTA": 0.07,
            "CTG": 0.41,
            "CAT": 0.41,
            "CAC": 0.59,
            "CAA": 0.25,
            "CAG": 0.75,
            "ATT": 0.36,
            "ATC": 0.48,
            "ATA": 0.16,
            "ATG": 1,
            "AAT": 0.46,
            "AAC": 0.54,
            "AAA": 0.42,
            "AAG": 0.58,
            "GTT": 0.18,
            "GTC": 0.24,
            "GTA": 0.11,
            "GTG": 0.47,
            "GAT": 0.46,
            "GAC": 0.54,
            "GAA": 0.42,
            "GAG": 0.58,
            "TCT": 0.18,
            "TCC": 0.22,
            "TCA": 0.15,
            "TCG": 0.06,
            "TGT": 0.45,
            "TGC": 0.55,
            "TGA": 0.52,
            "TGG": 1,
            "CCT": 0.28,
            "CCC": 0.33,
            "CCA": 0.27,
            "CCG": 0.11,
            "CGT": 0.08,
            "CGC": 0.19,
            "CGA": 0.11,
            "CGG": 0.21,
            "ACT": 0.24,
            "ACC": 0.36,
            "ACA": 0.28,
            "ACG": 0.12,
            "AGT": 0.15,
            "AGC": 0.24,
            "AGA": 0.2,
            "AGG": 0.2,
            "GCT": 0.26,
            "GCC": 0.4,
            "GCA": 0.23,
            "GCG": 0.11,
            "GGT": 0.16,
            "GGC": 0.34,
            "GGA": 0.25,
            "GGG": 0.25,
        }
    else:
        usage_table = usage
    return usage_table


def findORF(gene, non_interactive=False, preferred_orf_index=None):
    # Scan through all strands and frames for open reading frames
    min_protein_len = 100  # Support for finding gene position in vector
    genestart = []
    geneend = []
    genestrand = []
    print("Analyzing Gene:" + gene.name)
    for strand, nuc in [(+1, gene.seq), (-1, gene.seq.reverse_complement())]:
        for frame in range(3):
            length = 3 * ((len(gene) - frame) // 3)  # Multiple of three
            translated = nuc[frame : frame + length].translate()
            for protein in translated.split("*"):
                if len(protein) >= min_protein_len:
                    if len(protein.split("M", 1)) > 1:
                        ORF = "M" + protein.split("M", 1)[1]
                        genestart.append(translated.find(ORF) * 3 + frame + 1)
                        geneend.append(genestart[-1] + len(ORF) * 3 - 1)
                        genestrand.append(strand)
                        print(
                            "ORF#%i %s...%s - length %i, strand %i, frame %i"
                            % (
                                len(genestart),
                                ORF[:25],
                                ORF[-10:],
                                len(ORF),
                                strand,
                                frame + 1,
                            )
                        )
    if not genestart:
        raise ValueError(
            f"No valid ORF candidates found for {gene.name}. "
            "Provide explicit start/end coordinates in the FASTA header."
        )

    # Select Gene Frame
    if non_interactive:
        if preferred_orf_index is not None:
            genenum = int(preferred_orf_index)
            if genenum < 1 or genenum > len(genestart):
                raise ValueError(
                    f"Preferred ORF index {genenum} is out of range for {gene.name}. "
                    f"Found {len(genestart)} ORF candidates."
                )
        elif len(genestart) == 1:
            genenum = 1
        else:
            raise ValueError(
                f"Multiple ORF candidates found for {gene.name}. "
                "Use FASTA header start/end or provide an explicit ORF index."
            )
        if genestrand[genenum - 1] == -1:
            gene.seq = gene.seq.reverse_complement()
    else:
        while True:
            try:
                genenum = int(input("Which ORF are you targeting? (number):"))  # userinput
                if genestrand[genenum - 1] == -1:
                    gene.seq = gene.seq.reverse_complement()
            except (ValueError, IndexError):
                print("Please enter number")
                continue
            else:
                break
    start = genestart[genenum - 1] - 1  # subtract 1 to account for 0 indexing
    end = geneend[genenum - 1]
    if not non_interactive:
        print(gene.seq[start:end].translate()[:10])
        quest = "g"  # holding place

        while quest != "n" and quest != "y":
            quest = input("Is this the beginning of your gene?(position %i) (y/n):" % (start))
        while quest == "n":
            start = int(input("Enter the starting position your gene:"))
            print(gene.seq[start:end].translate()[:10])
            quest = input("Is this the beginning of your gene?(position %i) (y/n):" % (start))
        print(gene.seq[start:end].translate()[-10:])
        quest = "g"
        while quest != "n" and quest != "y":
            quest = input("Is the size of your gene %ibp? (y/n):" % (end - start))
        while quest == "n":
            try:
                end = int(input("Enter nucleotide length of your gene:")) + start
                if (end - start) % 3 != 0:
                    print("Length is not divisible by 3")
                    continue
                print(gene.seq[start:end].translate()[-10:])
                quest = input("Is this end correct? (y/n):")
            except ValueError:
                print("Please enter a number")
                quest = "n"
    return start, end
