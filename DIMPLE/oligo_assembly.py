"""Oligo assembly helpers for DIMPLE pipeline."""

from __future__ import annotations

from Bio.SeqRecord import SeqRecord

from DIMPLE.core import DIMPLE


def combine_fragments(tandem, num_frag_per_oligo, split):
    """TODO:
    Docstring
    """
    tandem_seq = []
    barcodes = []
    if split:
        tmpF = DIMPLE.barcodeF.pop(0)
        tmpR = DIMPLE.barcodeR.pop(0)
    direction = -1
    while len(tandem) > num_frag_per_oligo:
        tmp = tandem.pop(0)
        tmp_tandem = tmp.seq
        tandem_id = tmp.id
        for x in range(num_frag_per_oligo - 1):
            if split:
                name = tmp.id
                tmp = tandem.pop(0)
                if direction == 1:
                    tmp_tandem += (
                        "G"
                        + DIMPLE.cutsite.reverse_complement()
                        + tmpR
                        + tmpF
                        + DIMPLE.cutsite
                        + "C"
                        + tmp.seq
                    )  # concatenate and add cut sites with buffer
                    direction = -1
                else:
                    tmp_tandem += (
                        "G"
                        + DIMPLE.cutsite.reverse_complement()
                        + tmpR
                        + tmpF
                        + DIMPLE.cutsite
                        + "C"
                        + tmp.seq.reverse_complement()
                    )  # concatenate and add cut sites with buffer
                    direction = 1
                tandem_id += "+" + tmp.id
                barcodes.append(SeqRecord(tmpR, id=name))
                barcodes.append(SeqRecord(tmpF, id=tmp.id))
            else:
                tmp = tandem.pop(0)
                if direction == 1:
                    tmp_tandem += (
                        "G"
                        + DIMPLE.cutsite.reverse_complement()
                        + "ACGT"
                        + DIMPLE.cutsite
                        + "C"
                        + tmp.seq
                    )  # concatenate and add cut sites with buffer
                    direction = -1
                else:
                    tmp_tandem += (
                        "G"
                        + DIMPLE.cutsite.reverse_complement()
                        + "ACGT"
                        + DIMPLE.cutsite
                        + "C"
                        + tmp.seq.reverse_complement()
                    )
                    direction = 1
                tandem_id += "+" + tmp.id
        tandem_seq.append(SeqRecord(tmp_tandem, id=tandem_id, description=""))

    if tandem:
        direction = -1
        print(len(tandem))
        tmp = tandem.pop(0)
        tmp_tandem = tmp.seq
        tandem_id = tmp.id
        while tandem:
            if split:
                name = tmp.id
                tmp = tandem.pop(0)
                if direction == 1:
                    tmp_tandem += (
                        "G"
                        + DIMPLE.cutsite.reverse_complement()
                        + "ACGT"
                        + DIMPLE.cutsite
                        + "C"
                        + tmp.seq
                    )  # concatenate and add cut sites with buffer
                    direction = 1
                else:
                    tmp_tandem += (
                        "G"
                        + DIMPLE.cutsite.reverse_complement()
                        + "ACGT"
                        + DIMPLE.cutsite
                        + "C"
                        + tmp.seq.reverse_complement()
                    )  # concatenate and add cut sites with buffer
                    direction = -1
                tandem_id += "+" + tmp.id
                barcodes.append(SeqRecord(tmpR, id=name))
                barcodes.append(SeqRecord(tmpF, id=tmp.id))
            else:
                tmp = tandem.pop(0)
                if direction == -1:
                    tmp_tandem += (
                        "G"
                        + DIMPLE.cutsite.reverse_complement()
                        + "ACGT"
                        + DIMPLE.cutsite
                        + "C"
                        + tmp.seq
                    )  # concatenate and add cut sites with buffer
                    direction = 1
                else:
                    tmp_tandem += (
                        "G"
                        + DIMPLE.cutsite.reverse_complement()
                        + "ACGT"
                        + DIMPLE.cutsite
                        + "C"
                        + tmp.seq
                    )
                    direction = -1
                tandem_id += "+" + tmp.id
        difference = len(tandem_seq[-1].seq) - len(tmp_tandem)
        barF = DIMPLE.barcodeF.pop(0)
        barR = DIMPLE.barcodeR.pop(0)
        while difference / 2 > len(barF):
            barF += DIMPLE.barcodeF.pop(0)
            barR += DIMPLE.barcodeR.pop(0)
        tmpfrag = (
            barF.seq[0 : int(difference / 2)]
            + tmp_tandem
            + barR.seq.reverse_complement()[0 : difference - int(difference / 2)]
        )
        tandem_seq.append(SeqRecord(tmpfrag, id=tandem_id, description=""))
        print("Partial sequence" + str(len(tmpfrag)))
    return tandem_seq
