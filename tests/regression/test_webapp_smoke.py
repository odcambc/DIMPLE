"""Smoke test for the Gradio web layer (``webapp.app.run_dimple_job``).

Exercises the same wrapper the browser form calls -- no Gradio server, no
browser automation. Confirms a Kir DMS run returns a results zip with the
expected members, that the raw uploaded gene is cleaned up, and that the cheap
input-validation paths short-circuit before any pipeline work.

Skipped entirely when the ``web`` optional extra (gradio) isn't installed.
"""

import zipfile
from pathlib import Path

import pytest

gradio = pytest.importorskip("gradio")  # noqa: F841 - skip module if web extra absent

from webapp import app as webapp_app  # noqa: E402 - must follow importorskip

# Keyword args shared by the happy-path calls; mirrors the form defaults.
_BASE_KWARGS = dict(
    oligo_len=230,
    fragment_len="auto",
    overlap=4,
    dis=False,
    deletions_raw="",
    insertions_raw="",
    make_double=False,
    stop_codon=False,
    include_synonymous=False,
    maximize_nucleotide_change=False,
    restriction_sequence="CGTCTC(G)1/5",
    avoid_sequence="CGTCTC, GGTCTC",
    handle="AGCGGGAGACCGGGGTCTCTGAGC",
    codon_usage="human",
    barcode_start=0,
    tm_low=58,
    tm_high=62,
    orf_index="",
    random_seed="1",
)

_EXPECTED_ZIP_MEMBERS = {
    "All_Oligos.fasta",
    "All_Primers.fasta",
    "Kir_DMS_Gene_Primers.fasta",
    "Kir_DMS_Oligo_Primers.fasta",
    "Kir_DMS_Oligos.fasta",
    "Kir_designed_variants.csv",
    "Kir_mutations.csv",
}


def test_run_job_rejects_missing_file():
    """No upload -> error message, no crash, no zip."""
    zip_path, log = webapp_app.run_dimple_job(None, dms=True, **_BASE_KWARGS)
    assert zip_path is None
    assert "FASTA" in log


def test_run_job_rejects_no_mutation_type():
    """A file but no mutation type -> error before any pipeline work."""
    zip_path, log = webapp_app.run_dimple_job("dummy.fa", **{**_BASE_KWARGS, "dms": False})
    assert zip_path is None
    assert "mutation type" in log


@pytest.mark.slow
def test_run_job_kir_dms_produces_zip(tmp_path, kir_fa, monkeypatch):
    """End-to-end: Kir DMS run returns a zip with the expected members."""
    # Confine job dirs to the test's tmp_path so nothing leaks into the shared
    # app root; run_dimple_job reads JOBS_DIR at call time.
    jobs_dir = tmp_path / "jobs"
    monkeypatch.setattr(webapp_app, "JOBS_DIR", jobs_dir)

    zip_path, log = webapp_app.run_dimple_job(str(kir_fa), dms=True, **_BASE_KWARGS)

    assert zip_path is not None, log
    assert "Total oligos" in log
    assert "passed" in log  # advisory consistency gate ran and passed

    with zipfile.ZipFile(zip_path) as zf:
        members = set(zf.namelist())
    assert members == _EXPECTED_ZIP_MEMBERS

    # The raw uploaded gene must not survive in the job dir or the zip.
    assert kir_fa.name not in members
    assert not (Path(zip_path).parent / kir_fa.name).exists()
