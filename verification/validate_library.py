#!/usr/bin/env python3
"""One-off validator for a DIMPLE oligo library.

Independent of DIMPLE's own CSV bookkeeping: it reads the oligos you would
actually order (All_Oligos.fasta) + the WT gene FASTA, digests each oligo with
the Type IIS enzyme, reconstructs the mutant CDS by aligning the insert against
the WT gene, translates, and reports what actually changed.

The `verdict()` function -- the policy for what counts as "correct" -- is left
for you to write (see the TODO).
"""
from __future__ import annotations

import argparse
import re
from dataclasses import dataclass

from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Restriction import BsaI, BsmBI
from Bio.Seq import Seq

ENZYMES = {"BsmBI": BsmBI, "BsaI": BsaI}

_ALIGNER = PairwiseAligner()
_ALIGNER.mode = "local"
_ALIGNER.match_score = 2
_ALIGNER.mismatch_score = -1
_ALIGNER.open_gap_score = -6
_ALIGNER.extend_gap_score = -0.5

AA3TO1 = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Cys": "C",
    "Gln": "Q",
    "Glu": "E",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Phe": "F",
    "Pro": "P",
    "Ser": "S",
    "Thr": "T",
    "Trp": "W",
    "Tyr": "Y",
    "Val": "V",
    "Ter": "*",
    "End": "*",
    "Stop": "*",
    "STOP": "*",
}


# --------------------------------------------------------------------------- #
# Intent: what the oligo NAME claims the edit should be.
# --------------------------------------------------------------------------- #
@dataclass
class Intent:
    kind: str  # "sub" | "insert" | "delete"
    raw: str  # original name
    wt_aa: str | None = None  # sub only, 1-letter
    pos: int | None = None  # sub: residue; indel: wt_start from the name
    mut_aa: str | None = None  # sub only, 1-letter (or "*")
    ncodons: int | None = None  # indel only: number of codons added/removed
    ins_seq: str | None = None  # insert only: inserted nucleotides
    ins_aa: str | None = None  # insert only: translated inserted residues


# name shapes:
#   DMS:    Kir_DMS-<subpool>_<Wt3><pos><Mut3|STOP>
#   insert: Kir_insert-<subpool>_<inserted_nt>-<wt_start>
#   delete: Kir_delete-<subpool>_<n_nucleotides>-<wt_start>
_DMS = re.compile(r"_DMS-\d+_([A-Za-z]{3})(\d+)([A-Za-z]{3,4})$")
_INS = re.compile(r"_insert-\d+_([ACGT]+)-(\d+)$")
_DEL = re.compile(r"_delete-\d+_(\d+)-(\d+)$")


def parse_intent(name: str) -> Intent | None:
    m = _DMS.search(name)
    if m:
        wt, pos, mut = m.group(1), int(m.group(2)), m.group(3)
        return Intent("sub", name, AA3TO1.get(wt), pos, AA3TO1.get(mut))
    m = _INS.search(name)
    if m:
        seq = m.group(1)
        aa = str(Seq(seq).translate()) if len(seq) % 3 == 0 else None
        return Intent(
            "insert", name, pos=int(m.group(2)), ncodons=len(seq) // 3, ins_seq=seq, ins_aa=aa
        )
    m = _DEL.search(name)
    if m:
        return Intent("delete", name, pos=int(m.group(2)), ncodons=int(m.group(1)) // 3)
    return None


# --------------------------------------------------------------------------- #
# Observed: what the DNA actually encodes, derived independently.
# --------------------------------------------------------------------------- #
@dataclass
class Observed:
    placed: bool  # could we locate the insert in the WT CDS?
    note: str = ""
    frame_ok: bool | None = None  # is the net length change a multiple of 3?
    delta_nt: int | None = None  # mutant_cds_len - wt_cds_len
    wt_protein: str | None = None
    mut_protein: str | None = None
    # protein-level block diff (longest common prefix + suffix):
    edit_res: int | None = None  # 1-based residue where the edit block starts
    wt_block: str | None = None  # WT residues removed at the edit
    mut_block: str | None = None  # mutant residues present at the edit
    clean: bool | None = None  # do the flanks outside the block match WT exactly?
    # nucleotide-level, left-normalized (canonical) indel description:
    ins_nt: str | None = None  # nucleotides inserted (5'->3'), "" if none
    del_nt: str | None = None  # nucleotides deleted, "" if none
    norm_res: int | None = None  # 1-based residue the normalized edit starts at
    codon_aligned: bool | None = None  # does the normalized edit land on a codon boundary?
    mut_cds: str | None = None  # reconstructed mutant CDS (for whole-CDS indel checks)


def block_diff(wt: str, mut: str):
    """Return (start_0based, wt_block, mut_block): the minimal differing middle
    after stripping the longest common prefix and suffix."""
    n = min(len(wt), len(mut))
    p = 0
    while p < n and wt[p] == mut[p]:
        p += 1
    s = 0
    while s < (n - p) and wt[-1 - s] == mut[-1 - s]:
        s += 1
    return p, wt[p : len(wt) - s], mut[p : len(mut) - s]


def left_normalize(ref: str, pos: int, block: str) -> tuple[int, str]:
    """Left-align a pure indel block against ref (VCF-style).

    `block` is the inserted-or-deleted run occurring at `pos` in `ref` context;
    slide it as far 5' as the sequence allows so the representation is canonical.
    Returns (canonical_pos, canonical_block).
    """
    block = list(block)
    while pos > 0 and block[-1] == ref[pos - 1]:
        block = [ref[pos - 1]] + block[:-1]
        pos -= 1
    return pos, "".join(block)


def largest_fragment(dseq_pieces):
    return max(dseq_pieces, key=len)


def digest_insert(oligo_seq: str, enzyme) -> str:
    """Cut the oligo with the enzyme; return the largest internal piece (the insert)."""
    pieces = enzyme.catalyze(Seq(oligo_seq))
    return str(largest_fragment(pieces))


def _unique_pos(wt_gene: str, kmer: str) -> int | None:
    p = wt_gene.find(kmer)
    if p == -1 or wt_gene.find(kmer, p + 1) != -1:
        return None
    return p


def _exact_run(insert: str, wt_gene: str, off: int, start: int, step: int) -> int:
    """Walk from `start` in `step` direction while insert[i] == wt[i+off]; return
    the index just past the last match (step=+1) or just before it (step=-1)."""
    n = len(insert)
    i = start
    while 0 <= i < n and 0 <= i + off < len(wt_gene) and insert[i] == wt_gene[i + off]:
        i += step
    return i


def reconstruct(insert: str, wt_gene: str, k: int = 15, max_indel: int = 60):
    """Splice the insert into the WT top-strand gene via its two WT flank offsets.

    The coding insert is [synthetic?][5' WT flank][edit][3' WT flank]. Each flank
    matches WT at an offset (wt_pos - insert_index); the 5' flank fixes the left
    splice point, the 3' flank the right one, and the difference between the two
    offsets *is* the edit -- equal => substitution, 3' offset smaller => insertion,
    larger => deletion. Both flanks are recovered by exact extension, so this works
    down to the 4-nt Type IIS overhang (edit at a fragment edge). If the 3' flank
    never resumes in WT, that is a real mis-generation. Returns (mutant, note) or
    (None, reason).
    """
    from collections import Counter

    n, W = len(insert), len(wt_gene)
    anchors = [
        (i, p - i)
        for i in range(n - k + 1)
        if (p := _unique_pos(wt_gene, insert[i : i + k])) is not None
    ]
    if not anchors:
        return None, "no unique anchor"

    def first_match(off):
        i = 0
        while i < n and not (0 <= i + off < W and insert[i] == wt_gene[i + off]):
            i += 1
        return i

    # 5' flank offset: the offset of the earliest unique anchor (skips any
    # synthetic Kozak/handle, which yields no WT match). NOTE: this assumes the 5'
    # flank is long enough to hold a unique k-mer; edits within a few codons of a
    # fragment's 5' edge have too short a 5' flank and fall to the pydna assembly
    # cross-check instead (their upstream context lives in the backbone, not the
    # oligo).
    off5 = min(anchors, key=lambda a: a[0])[1]
    a5 = first_match(off5)

    # Substitution vs indel: a substitution changes at most one codon, so the whole
    # insert (from a5) matches WT at off5 with <=3 mismatches; an indel frameshifts
    # off5 and shows many. This cleanly separates them even when the edit is at the
    # 3' fragment edge (where only the 4-nt overhang remains as a 3' flank).
    mismatches = sum(
        1 for i in range(a5, n) if not (0 <= i + off5 < W) or insert[i] != wt_gene[i + off5]
    )
    if mismatches <= 3:
        off3 = off5  # substitution / synonymous / no indel
        best_run = n - a5
    else:
        # indel: the 3' flank sits at a different offset. Pick the offset with the
        # longest exact suffix run (works down to the 4-nt overhang).
        off3, best_run = off5, -1
        for off in range(off5 - max_indel, off5 + max_indel + 1):
            if off == off5:
                continue
            run = (n - 1) - _exact_run(insert, wt_gene, off, n - 1, -1)
            if run > best_run or (run == best_run and abs(off - off5) < abs(off3 - off5)):
                best_run, off3 = run, off

    A, B = a5 + off5, n + off3
    if best_run < 4 or not (0 <= A <= B <= W):
        return None, "3' flank does not resume in WT (possible mis-generation)"
    return wt_gene[:A] + insert[a5:] + wt_gene[B:], ""


def translate_cds(gene: str, orf_start: int, orf_end: int) -> str:
    """Translate from orf_start (0-based) to the first stop or orf_end."""
    cds = gene[orf_start:orf_end]
    prot = str(Seq(cds).translate(to_stop=False))
    return prot


def observe(insert: str, wt_gene: str, orf_start: int, orf_end: int) -> Observed:
    mutant, note = reconstruct(insert, wt_gene)
    if mutant is None:
        return Observed(placed=False, note=note)

    delta = len(mutant) - len(wt_gene)
    frame_ok = delta % 3 == 0

    wt_prot = translate_cds(wt_gene, orf_start, orf_end).split("*")[0]
    # mutant ORF end shifts by delta (indels move the stop codon downstream)
    mut_prot = str(Seq(mutant[orf_start : orf_end + delta]).translate(to_stop=True))

    p, wt_block, mut_block = block_diff(wt_prot, mut_prot)
    edit_res = p + 1
    # "clean" = the block diff fully explains the difference (flanks identical),
    # i.e. no frameshift scrambling the downstream sequence.
    clean = (
        wt_prot[:p] == mut_prot[:p]
        and wt_prot[p + len(wt_block) :] == mut_prot[p + len(mut_block) :]
    )

    # Nucleotide-level, left-normalized indel description (canonical position).
    wt_cds = wt_gene[orf_start:orf_end]
    mut_cds = mutant[orf_start : orf_end + delta]
    np_, del_nt, ins_nt = block_diff(wt_cds, mut_cds)
    if ins_nt and not del_nt:  # pure insertion
        np_, ins_nt = left_normalize(wt_cds, np_, ins_nt)
    elif del_nt and not ins_nt:  # pure deletion
        np_, del_nt = left_normalize(wt_cds, np_, del_nt)
    norm_res = np_ // 3 + 1
    codon_aligned = np_ % 3 == 0

    return Observed(
        placed=True,
        frame_ok=frame_ok,
        delta_nt=delta,
        wt_protein=wt_prot,
        mut_protein=mut_prot,
        edit_res=edit_res,
        wt_block=wt_block,
        mut_block=mut_block,
        clean=clean,
        ins_nt=ins_nt,
        del_nt=del_nt,
        norm_res=norm_res,
        codon_aligned=codon_aligned,
        mut_cds=mut_cds,
    )


# --------------------------------------------------------------------------- #
# VERDICT
# --------------------------------------------------------------------------- #
def verdict(intent: Intent, obs: Observed, wt_cds: str) -> tuple[bool | None, str]:
    """Decide whether the observed edit matches the intent. Returns (ok, reason).

    Policy (agreed):
      * substitutions -- protein-level: right residue change at the named position,
        synonymous => no protein change, nonsense => truncation at the position.
      * indels -- whole-CDS: reconstruct the expected mutant CDS from the name
        (WT with the named handle spliced in after codon wt_start, or the named
        codon range deleted) and require an exact nucleotide match. This subsumes
        position + frame + inserted-nucleotide identity in one ambiguity-free test.
    """
    if not obs.placed:
        return None, f"unplaced ({obs.note})"
    if not obs.frame_ok:
        return False, f"frameshift: delta_nt={obs.delta_nt}"

    if intent.kind == "sub":
        if not obs.clean:
            return False, "downstream sequence altered (not a clean point change)"
        if intent.wt_aa == intent.mut_aa:  # synonymous
            ok = obs.wt_block == "" and obs.mut_block == ""
            return ok, (
                "synonymous: no protein change"
                if ok
                else f"synonymous but protein changed: {obs.wt_block}->{obs.mut_block}"
            )
        if intent.mut_aa == "*":  # nonsense
            ok = obs.edit_res == intent.pos and obs.mut_block == ""
            return ok, (
                f"nonsense: truncated at {intent.pos}"
                if ok
                else f"nonsense mismatch: edit_res={obs.edit_res} mut_block={obs.mut_block!r}"
            )
        ok = (
            obs.edit_res == intent.pos
            and obs.wt_block == intent.wt_aa
            and obs.mut_block == intent.mut_aa
        )
        return ok, (
            f"{intent.wt_aa}{intent.pos}{intent.mut_aa} confirmed"
            if ok
            else f"want {intent.wt_aa}{intent.pos}{intent.mut_aa}, "
            f"got {obs.wt_block}{obs.edit_res}{obs.mut_block}"
        )

    if intent.kind == "insert":
        exp = wt_cds[: 3 * intent.pos] + intent.ins_seq + wt_cds[3 * intent.pos :]
        ok = exp == obs.mut_cds
        return ok, (
            f"+{intent.ins_seq} after codon {intent.pos} (exact CDS)"
            if ok
            else f"+{intent.ins_seq}@{intent.pos}: CDS mismatch"
        )

    if intent.kind == "delete":
        ndel = intent.ncodons * 3
        cut = 3 * (intent.pos - 1)
        exp = wt_cds[:cut] + wt_cds[cut + ndel :]
        ok = exp == obs.mut_cds
        return ok, (
            f"-{ndel}nt at codon {intent.pos} (exact CDS)"
            if ok
            else f"-{ndel}nt@{intent.pos}: CDS mismatch"
        )

    return None, "unknown intent kind"


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--oligos", required=True)
    ap.add_argument("--gene", required=True)
    ap.add_argument("--enzyme", default="BsmBI", choices=list(ENZYMES))
    ap.add_argument("--limit", type=int, default=0, help="cap records for a quick look")
    ap.add_argument("--show", type=int, default=8, help="print N example observations")
    args = ap.parse_args()

    enzyme = ENZYMES[args.enzyme]
    rec = next(SeqIO.parse(args.gene, "fasta"))
    header = rec.description
    m = re.search(r"start:(\d+)\s+end:(\d+)", header)
    orf_start = int(m.group(1)) - 1 if m else 0
    orf_end = int(m.group(2)) if m else len(rec.seq)
    wt_gene = str(rec.seq).upper()

    wt_cds = wt_gene[orf_start:orf_end]
    n_res = len(wt_cds) // 3

    counts = {"sub": 0, "insert": 0, "delete": 0, "unparsed": 0}
    passed = {"sub": 0, "insert": 0, "delete": 0}
    unplaced = {"sub": 0, "insert": 0, "delete": 0}
    fails = []  # (id, kind, reason, pos, near_edge)
    unparsed_examples = []

    for o in SeqIO.parse(args.oligos, "fasta"):
        if args.limit and sum(counts.values()) >= args.limit:
            break
        intent = parse_intent(o.id)
        if intent is None:
            counts["unparsed"] += 1
            if len(unparsed_examples) < 8:
                unparsed_examples.append(o.id)
            continue
        counts[intent.kind] += 1
        insert = digest_insert(str(o.seq).upper(), enzyme)
        obs = observe(insert, wt_gene, orf_start, orf_end)
        ok, reason = verdict(intent, obs, wt_cds)
        if ok is True:
            passed[intent.kind] += 1
        elif ok is None:
            unplaced[intent.kind] += 1
        else:
            # is the edit within 1 residue of a subpool/CDS boundary?
            pos = intent.pos or 0
            near_edge = pos <= 1 or pos >= n_res - 1
            fails.append((o.id, intent.kind, reason, pos, near_edge))

    print("\n=== validation summary ===")
    for k in ("sub", "insert", "delete"):
        print(f"  {k:7s}: {passed[k]:5d}/{counts[k]:<5d} confirmed" f"  ({unplaced[k]} unplaced)")
    print(f"  unparsed names: {counts['unparsed']}")
    if unparsed_examples:
        print("    e.g.", unparsed_examples[:5])

    if fails:
        edge = [f for f in fails if f[4]]
        interior = [f for f in fails if not f[4]]
        print(
            f"\n=== {len(fails)} verdict FAILURES "
            f"({len(edge)} at a CDS/subpool edge, {len(interior)} interior) ==="
        )
        print("  -- interior (most concerning) --")
        for fid, kind, reason, pos, _ in interior[:15]:
            print(f"    {fid:30s} {reason}")
        print("  -- edge (likely reconstruction boundary artifacts) --")
        for fid, kind, reason, pos, _ in edge[:10]:
            print(f"    {fid:30s} {reason}")
    else:
        print("\nAll parsed variants confirmed. No mismatches.")


if __name__ == "__main__":
    main()
