"""Generate SnapGene-openable GenBank (.gb) files for visual inspection.

SnapGene opens GenBank natively with full feature + primer rendering (there is no
robust open-source writer for the native .dna binary). Produces:

  Kir_gene_template.gb        -- the gene template with all 9 fragments + gene primers
  <variant>_oligo.gb          -- one designed oligo, annotated (barcodes/BsmBI/overhang/coding/mutation)
  <variant>_assembled.gb      -- that variant assembled back into the template, mutation marked

Run: python verification/make_snapgene.py   (writes to verification/snapgene_examples/)
"""

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from build_pool import build

OUT = HERE / "snapgene_examples"
CUTSITE = "CGTCTC"
CUTSITE_RC = "GAGACG"
OVERHANG = 4
OVERLAP = 4

# a substitution in the (previously broken) fragment 7, an insertion, a deletion
EXAMPLES = ["Kir_DMS-7_Gln293Cys", "Kir_insert-2_GAC-60", "Kir_delete-2_9-60"]


def feat(start, end, kind, label, strand=1):
    return SeqFeature(
        FeatureLocation(int(start), int(end), strand=strand),
        type=kind,
        qualifiers={"label": label},
    )


def anneal_region(primer_seq):
    """Return the 3' annealing tail of a gene primer (drop the ATA+CGTCTC+buffer)."""
    s = str(primer_seq).upper()
    i = s.find(CUTSITE)
    return s[i + len(CUTSITE) + 1 + OVERHANG :] if i >= 0 else s


def add_primer_binds(record, template, primers):
    tmpl = str(template).upper()
    for p in primers:
        tail = anneal_region(p.seq)[-18:]
        if not tail:
            continue
        j = tmpl.find(tail)
        if j >= 0:
            record.features.append(feat(j, j + len(tail), "primer_bind", p.id, strand=1))
            continue
        rc = str(Seq(tail).reverse_complement())
        j = tmpl.find(rc)
        if j >= 0:
            record.features.append(feat(j, j + len(rc), "primer_bind", p.id, strand=-1))


def gb_record(seq, name, desc, circular=False):
    rec = SeqRecord(Seq(str(seq)), id=name[:16], name=name[:16], description=desc)
    rec.annotations["molecule_type"] = "ds-DNA"
    rec.annotations["topology"] = "circular" if circular else "linear"
    return rec


def make_template(gene):
    rec = gb_record(gene.seq, "Kir_template", "Kir gene template: fragments + gene primers")
    rec.features.append(
        feat(gene.breaklist[0][0], gene.breaklist[-1][1], "CDS", "ORF (mutable region)")
    )
    for idx, (a, b) in enumerate(gene.breaklist):
        rec.features.append(feat(a, b, "misc_feature", f"Fragment {idx + 1}"))
    add_primer_binds(rec, gene.seq, gene.genePrimer)
    SeqIO.write(rec, OUT / "Kir_gene_template.gb", "genbank")


def make_oligo(gene, vid):
    e = gene.designed_variants[vid]
    oseq = str(e["oligo_sequence"]).upper()
    rec = gb_record(oseq, vid.replace("Kir_", ""), f"{vid} designed oligo")
    p1 = oseq.find(CUTSITE)
    p2 = oseq.rfind(CUTSITE_RC)
    if p1 >= 0 and p2 > p1:
        rec.features.append(feat(0, p1, "misc_feature", "barcode primer F"))
        rec.features.append(feat(p1, p1 + 6, "protein_bind", "BsmBI"))
        rec.features.append(feat(p2, p2 + 6, "protein_bind", "BsmBI"))
        rec.features.append(feat(p2 + 6, len(oseq), "misc_feature", "barcode primer R"))
        # coding insert = between the two cut+buffer offsets; overhangs are the 4 nt each side
        ins_start, ins_end = p1 + 6 + 1, p2 - 1
        rec.features.append(feat(ins_start, ins_start + OVERHANG, "misc_feature", "5' overhang"))
        rec.features.append(feat(ins_end - OVERHANG, ins_end, "misc_feature", "3' overhang"))
        rec.features.append(feat(ins_start, ins_end, "CDS", f"insert ({vid.split('_')[-1]})"))
    # mark the edited bases: diff the mutant coding core vs the WT coding core
    a, b = gene.breaklist[e["fragment"] - 1]
    wt_frag = str(gene.seq)[a - OVERHANG - OVERLAP : b + OVERHANG + OVERLAP].upper()
    xfrag = str(e["xfrag"]).upper()
    core, wt_core = xfrag[OVERHANG:-OVERHANG], wt_frag[OVERHANG:-OVERHANG]
    cpos = oseq.find(core[:24])
    if cpos >= 0 and core != wt_core:
        p = 0
        while p < min(len(core), len(wt_core)) and core[p] == wt_core[p]:
            p += 1
        s = 0
        while s < (min(len(core), len(wt_core)) - p) and core[-1 - s] == wt_core[-1 - s]:
            s += 1
        rec.features.append(
            feat(
                cpos + p, cpos + p + max(1, len(core) - p - s), "modified_base", vid.split("_")[-1]
            )
        )
    SeqIO.write(rec, OUT / f"{vid}_oligo.gb", "genbank")


def make_assembled(gene, vid):
    e = gene.designed_variants[vid]
    frag = e["fragment"]
    a, b = gene.breaklist[frag - 1]
    A = a - OVERHANG - OVERLAP
    B = b + OVERHANG + OVERLAP
    xfrag = str(e["xfrag"]).upper()
    wt = str(gene.seq).upper()
    mutant = wt[:A] + xfrag + wt[B:]
    rec = gb_record(
        mutant, vid.replace("Kir_", "") + "_asm", f"{vid} assembled into template", circular=True
    )
    # fragment span (shifts downstream bases by the indel, but the insert sits at [A:A+len(xfrag)])
    rec.features.append(feat(A, A + len(xfrag), "misc_feature", f"Fragment {frag} (assembled)"))
    # locate the edit by trimming the common prefix/suffix of mutant vs wt
    p = 0
    while p < min(len(mutant), len(wt)) and mutant[p] == wt[p]:
        p += 1
    s = 0
    while s < (min(len(mutant), len(wt)) - p) and mutant[-1 - s] == wt[-1 - s]:
        s += 1
    edit_len = max(1, len(mutant) - p - s)
    rec.features.append(feat(p, p + edit_len, "modified_base", vid.split("_")[-1]))
    add_primer_binds(rec, mutant, gene.genePrimer)
    SeqIO.write(rec, OUT / f"{vid}_assembled.gb", "genbank")


def main():
    OUT.mkdir(exist_ok=True)
    gene = build()[0]
    make_template(gene)
    for vid in EXAMPLES:
        if vid in gene.designed_variants:
            make_oligo(gene, vid)
            make_assembled(gene, vid)
        else:
            print(f"(skip {vid}: not in library)")
    print("wrote:")
    for f in sorted(OUT.glob("*.gb")):
        print("  ", f)


if __name__ == "__main__":
    main()
