"""pydna assembly cross-check (method 2).

For each variant: PCR the WT gene with its subpool gene primers -> BsmBI cut ->
backbone; PCR the oligo with its oligo primers -> BsmBI cut -> insert; ligate ->
circular product; rotate to WT frame; translate the ORF; compare the observed
protein edit to the intent. This uses pydna's real sticky-end ligation, so it is
independent of the flank-offset reconstruction in validate_library.

Run: builds backbones once, then assembles the 204 residual + a 100-variant
sample of already-passing variants, and reports agreement with method 1.
"""

import sys, re, random
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO))

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Restriction import BsmBI
from pydna.dseqrecord import Dseqrecord
from pydna.amplify import pcr

import validate_library as V
from build_pool import build


def max_piece(pieces):
    return max(pieces, key=len)


def assemble_and_translate(gene, backbones, dv_entry, wt_prot, edit_pos):
    """Ligate the oligo insert into its subpool backbone and return the mutant
    protein, located coordinate-free: translate the (circular) product in all 3
    frames on both strands, and among every Met-initiated ORF pick the one whose
    block-diff against WT is the smallest single local edit (self-validating -- the
    correct ORF start yields a tiny diff, wrong starts yield garbage)."""
    frag = dv_entry["fragment"]
    oligo = Dseqrecord(dv_entry["oligo_sequence"])
    fwd, rev = gene.barPrimer[(frag - 1) * 2], gene.barPrimer[(frag - 1) * 2 + 1]
    insert = max_piece(Dseqrecord(pcr(fwd, rev, oligo)).cut(BsmBI))
    assembled = (insert + backbones[frag - 1]).looped()

    L = len(wt_prot)
    best = None  # (diff_score, orf)
    seq = str(assembled.seq).upper()
    for cand in (seq, str(Seq(seq).reverse_complement())):
        dbl = cand + cand
        for frame in (0, 1, 2):
            prot = str(Seq(dbl[frame : frame + (len(dbl) - frame) // 3 * 3]).translate())
            for mpos in (i for i, c in enumerate(prot) if c == "M"):
                orf = prot[mpos:].split("*")[0]
                if abs(len(orf) - L) > 30 or len(orf) < L - 30:
                    continue
                p, wtb, mutb = V.block_diff(wt_prot, orf)
                score = len(wtb) + len(mutb)
                if best is None or score < best[0]:
                    best = (score, orf)
    return best[1] if best else None


def main():
    gene = build()[0]
    wt_gene = str(gene.seq).upper() if hasattr(gene, "seq") else None
    rec = next(SeqIO.parse(str(REPO / "tests/expected/Kir.fa"), "fasta"))
    m = re.search(r"start:(\d+)\s+end:(\d+)", rec.description)
    orf_start, orf_end = int(m.group(1)) - 1, int(m.group(2))
    wt_gene = str(rec.seq).upper()
    wt_cds = wt_gene[orf_start:orf_end]
    wt_prot = str(Seq(wt_cds).translate(to_stop=True))

    # build the 9 backbones once
    full_template = Dseqrecord(gene.seq, circular=True)
    n_frag = len(gene.genePrimer) // 2
    backbones = []
    for f in range(n_frag):
        prod = Dseqrecord(pcr(gene.genePrimer[f * 2], gene.genePrimer[f * 2 + 1], full_template))
        backbones.append(max_piece(prod.cut(BsmBI)))
    print(f"built {len(backbones)} backbones", flush=True)

    # method-1 verdicts on the golden oligos -> residual + passing sets
    residual, passing = [], []
    for r in SeqIO.parse(str(REPO / "tests/expected/All_Oligos.fasta"), "fasta"):
        it = V.parse_intent(r.id)
        if not it:
            continue
        ob = V.observe(V.digest_insert(str(r.seq).upper(), BsmBI), wt_gene, orf_start, orf_end)
        ok, _ = V.verdict(it, ob, wt_cds)
        (residual if ok is not True else passing).append(r.id)
    rng = random.Random(1848)
    sample = rng.sample(passing, 100)
    print(f"method-1: {len(residual)} residual, {len(passing)} passing; sampling 100", flush=True)

    dv = gene.designed_variants

    def pydna_verdict(vid):
        it = V.parse_intent(vid)
        try:
            prot = assemble_and_translate(gene, backbones, dv[vid], wt_prot, it.pos or 0)
        except Exception as exc:  # noqa
            return None, f"assembly error: {type(exc).__name__}: {str(exc)[:80]}"
        if prot is None:
            return None, "could not locate ORF in assembled product"
        p, wtb, mutb = V.block_diff(wt_prot, prot)
        clean = wt_prot[:p] == prot[:p] and wt_prot[p + len(wtb) :] == prot[p + len(mutb) :]
        # reuse the same policy as method 1 by synthesizing an Observed
        obs = V.Observed(
            placed=True,
            frame_ok=(len(prot) != len(wt_prot)) or True,
            delta_nt=(len(prot) - len(wt_prot)) * 3,
            wt_protein=wt_prot,
            mut_protein=prot,
            edit_res=p + 1,
            wt_block=wtb,
            mut_block=mutb,
            clean=clean,
            mut_cds=None,
        )
        # protein-level judgement (indels compared at protein level here)
        if it.kind == "sub":
            return V.verdict(it, obs, wt_cds)
        # indel: check clean in-frame block of right size/residues near the position
        if not clean:
            return False, "downstream altered"
        if it.kind == "insert":
            ok = wtb == "" and len(mutb) == it.ncodons and str(Seq(it.ins_seq).translate()) == mutb
            return ok, f"+{mutb} (want {str(Seq(it.ins_seq).translate())})"
        ok = mutb == "" and len(wtb) == it.ncodons
        return ok, f"-{wtb}"

    def run(ids, label):
        agree = disagree = err = 0
        details = []
        for i, vid in enumerate(ids):
            pk, preason = pydna_verdict(vid)
            it = V.parse_intent(vid)
            ob = V.observe(
                (
                    V.digest_insert(dv[vid]["oligo_sequence"], BsmBI)
                    if False
                    else V.digest_insert(dv[vid]["oligo_sequence"], BsmBI)
                ),
                wt_gene,
                orf_start,
                orf_end,
            )
            mk, _ = V.verdict(it, ob, wt_cds)
            if pk is None:
                err += 1
                details.append((vid, "PYDNA-ERR", preason))
            elif pk is True:
                agree += 1  # pydna says correct
                if mk is not True:
                    details.append((vid, "PYDNA-OK / method1-FAIL", preason))
            else:
                disagree += 1
                details.append((vid, "PYDNA-FAIL", preason))
            if (i + 1) % 25 == 0:
                print(f"  [{label}] {i+1}/{len(ids)}", flush=True)
        print(
            f"\n=== {label}: pydna OK={agree}  pydna FAIL={disagree}  errors={err} (n={len(ids)}) ==="
        )
        for d in details[:60]:
            print("   ", d)
        return agree, disagree, err

    print("\n########## RESIDUAL (method-1 could not confirm) ##########", flush=True)
    run(residual, "residual")
    print("\n########## SAMPLE of method-1 PASSING (agreement check) ##########", flush=True)
    run(sample, "sample")


if __name__ == "__main__":
    main()
