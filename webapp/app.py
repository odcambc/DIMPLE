"""Gradio front-end for DIMPLE.

Wraps the same pipeline the CLI and Tk GUI drive
(:func:`DIMPLE.runner.build_runtime_config` + :func:`DIMPLE.runner.run_pipeline`)
behind a browser form. A user uploads a FASTA, sets the design knobs, and gets
back the oligo / primer FASTAs + variant / mutation CSVs as a zip.

Design notes (see ``tasks/webapp.md``):

* **Always non-interactive.** ``non_interactive=True`` plus ``link_policy="never"``
  and ``breaksite_change_policy="error"`` so the ``input()`` paths in
  ``findORF`` / the breaksite setter never block the server. ORF disambiguation
  comes from the FASTA ``start:/end:`` header or the explicit ORF-index field.
* **Per-request temp dir.** Each run gets a fresh ``tempfile.mkdtemp`` work dir;
  ``run_pipeline`` reads/writes files there and we zip the results out.
* This is the fast-demo layer. Out-of-process execution, upload-size caps, and a
  hard wall-clock timeout are the step-4 robustness items, not done here.
"""

from __future__ import annotations

import logging
import os
import shutil
import tempfile
import threading
import time
import zipfile
from pathlib import Path

# Disable Gradio's network telemetry before importing gradio -- the analytics
# env var is read at import time, so setting it later (or only on Blocks) leaves
# the import-time version check firing.
os.environ.setdefault("GRADIO_ANALYTICS_ENABLED", "False")

import gradio as gr  # noqa: E402  (must follow the env var above)

from DIMPLE.runner import build_runtime_config, run_pipeline  # noqa: E402

logger = logging.getLogger(__name__)

# --- Multi-user data isolation -------------------------------------------------
# DIMPLE itself holds no process-global mutable state (per-instance RNG,
# config-scoped everything), so in-memory leakage between concurrent requests is
# not a concern. The leak surface is the filesystem + Gradio's file-serving
# layer, which this module owns. Layout:
#
#   APP_ROOT/jobs/<random>/   one isolated work dir per request; reaped on a TTL
#   APP_ROOT/gradio/          Gradio's own cache (uploads + served copies)
#
# Each run's uploaded gene FASTA is deleted as soon as the pipeline finishes, and
# whole job dirs are removed once older than JOB_TTL_SECONDS. Gradio's cache (the
# raw upload it stores, plus the served copy of the result) is reaped by
# delete_cache on the Blocks. No work dir is ever added to Gradio's allowed
# paths, so there is no /gradio_api/file= route into it.
APP_ROOT = Path(os.environ.get("DIMPLE_WEB_TMP", Path(tempfile.gettempdir()) / "dimple_web"))
JOBS_DIR = APP_ROOT / "jobs"
GRADIO_CACHE_DIR = APP_ROOT / "gradio"

# How long a finished job's files remain fetchable before the reaper deletes
# them. Bounds the window in which a result could be served at all.
JOB_TTL_SECONDS = int(os.environ.get("DIMPLE_WEB_JOB_TTL", "1800"))  # 30 min

# Cap upload size; oligo generation is CPU-bound and a huge FASTA is a DoS.
MAX_UPLOAD_MB = int(os.environ.get("DIMPLE_WEB_MAX_UPLOAD_MB", "10"))

# Output artifacts live in the work dir alongside the uploaded FASTA; collect
# everything DIMPLE emits (FASTA + CSV) and exclude the input we wrote in.
_OUTPUT_GLOBS = ("*.fasta", "*.csv")


def _reap_old_jobs(ttl_seconds: int = JOB_TTL_SECONDS) -> None:
    """Delete job dirs whose mtime is older than *ttl_seconds*. Best-effort."""
    if not JOBS_DIR.is_dir():
        return
    cutoff = time.time() - ttl_seconds
    for job_dir in JOBS_DIR.iterdir():
        try:
            if job_dir.is_dir() and job_dir.stat().st_mtime < cutoff:
                shutil.rmtree(job_dir, ignore_errors=True)
        except OSError:
            pass


def _start_reaper(interval_seconds: int = 300) -> None:
    """Start a daemon thread that periodically reaps expired job dirs."""

    def _loop() -> None:
        while True:
            time.sleep(interval_seconds)
            _reap_old_jobs()

    threading.Thread(target=_loop, name="dimple-web-reaper", daemon=True).start()


def _parse_int_list(raw: str) -> list[int]:
    return [int(x.strip()) for x in raw.split(",") if x.strip()]


def _parse_str_list(raw: str) -> list[str]:
    return [x.strip() for x in raw.split(",") if x.strip()]


def run_dimple_job(
    fasta_file,
    *,
    oligo_len: int,
    fragment_len: str,
    overlap: int,
    dms: bool,
    dis: bool,
    deletions_raw: str,
    insertions_raw: str,
    make_double: bool,
    stop_codon: bool,
    include_synonymous: bool,
    maximize_nucleotide_change: bool,
    restriction_sequence: str,
    avoid_sequence: str,
    handle: str,
    codon_usage: str,
    barcode_start: int,
    tm_low: int,
    tm_high: int,
    orf_index,
    random_seed,
) -> tuple[str | None, str]:
    """Run one DIMPLE design and return ``(zip_path_or_None, status_log)``.

    Returns the path to a results zip on success, or ``None`` plus an error
    message in the status log on failure. Never raises out to Gradio -- all
    errors are surfaced as text so the UI shows a message, not a traceback.
    """
    log: list[str] = []

    if fasta_file is None:
        return None, "Error: please upload a FASTA gene file."

    work_dir: Path | None = None
    try:
        deletions = _parse_int_list(deletions_raw) if deletions_raw.strip() else False
        insertions = _parse_str_list(insertions_raw) if insertions_raw.strip() else False

        if not any([dms, dis, deletions, insertions]):
            return None, (
                "Error: select at least one mutation type " "(DMS, DIS, deletions, or insertions)."
            )

        # Isolated per-request work dir under the app-owned root; the uploaded FASTA
        # is copied in under its own name so run_pipeline's
        # os.path.join(work_dir, target_file) resolves.
        JOBS_DIR.mkdir(parents=True, exist_ok=True)
        work_dir = Path(tempfile.mkdtemp(prefix="job_", dir=JOBS_DIR))
        src = Path(fasta_file)
        gene_file = src.name
        shutil.copy(src, work_dir / gene_file)

        frag = 0 if str(fragment_len).strip().lower() in ("", "auto") else int(fragment_len)
        seed = int(random_seed) if str(random_seed).strip() not in ("", "none") else None
        orf = int(orf_index) if str(orf_index).strip() not in ("", "none") else None

        config, overlap_l, overlap_r = build_runtime_config(
            oligo_len=int(oligo_len),
            fragment_len=frag,
            overlap=int(overlap),
            handle=handle,
            restriction_sequence=restriction_sequence,
            avoid_sequence=avoid_sequence,
            codon_usage=codon_usage,
            barcode_start=int(barcode_start),
            deletions=deletions,
            dms=bool(dms),
            dis=bool(dis),
            make_double=bool(make_double),
            stop_codon=bool(stop_codon),
            maximize_nucleotide_change=bool(maximize_nucleotide_change),
            # --- web layer always forces non-interactive ---
            non_interactive=True,
            preferred_orf_index=orf,
            link_policy="never",
            breaksite_change_policy="error",
            random_seed=seed,
            logger=logger,
        )

        if config.enzyme is not None:
            log.append(f"Restriction enzyme: {config.enzyme}")
        log.append(f"Restriction sequence: {config.cutsite}")

        pool = run_pipeline(
            gene_file,
            str(work_dir),
            config,
            overlap_l,
            overlap_r,
            include_synonymous=bool(include_synonymous),
            insertions=insertions,
            deletions=deletions,
            dis=bool(dis),
            gene_primer_tm=(int(tm_low), int(tm_high)),
        )
    except Exception as exc:  # noqa: BLE001 - surface every failure as UI text
        logger.exception("DIMPLE web run failed")
        if work_dir is not None:
            shutil.rmtree(work_dir, ignore_errors=True)
        return None, "\n".join(log + [f"Error: {type(exc).__name__}: {exc}"])

    # The pipeline has read the upload; remove the raw gene FASTA now so the
    # plaintext input does not linger in the job dir for the TTL window.
    (work_dir / gene_file).unlink(missing_ok=True)

    genes = [obj.geneid for obj in pool]
    n_oligos = sum(len(getattr(obj, "oligos", [])) for obj in pool)
    log.append(f"Genes processed: {', '.join(genes)}")
    log.append(f"Total oligos: {n_oligos}")

    # Advisory consistency gate (DMS-shape filenames only; don't fail the run).
    if dms:
        try:
            from tests.regression._helpers import assert_outputs_consistent

            for gene_id in genes:
                assert_outputs_consistent(work_dir, gene_id, config)
            log.append("Output consistency check: passed.")
        except AssertionError as exc:
            log.append(f"Warning - output consistency check failed: {exc}")
        except Exception:  # noqa: BLE001 - helper is best-effort here
            pass

    out_files = sorted(
        p for pat in _OUTPUT_GLOBS for p in work_dir.glob(pat) if p.name != gene_file
    )
    if not out_files:
        return None, "\n".join(log + ["Error: run produced no output files."])

    zip_path = work_dir / "DIMPLE_results.zip"
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for p in out_files:
            zf.write(p, p.name)

    log.append("")
    log.append("Output files:")
    log.extend(f"  - {p.name}" for p in out_files)
    return str(zip_path), "\n".join(log)


def build_interface() -> gr.Blocks:
    # delete_cache=(frequency, age): Gradio sweeps its own cache (the raw upload
    # it stores + the served copy of each result) every `frequency` seconds,
    # removing entries older than `age`. Bounds how long any served file lives.
    with gr.Blocks(
        title="DIMPLE", delete_cache=(300, JOB_TTL_SECONDS), analytics_enabled=False
    ) as demo:
        gr.Markdown(
            "# DIMPLE\n"
            "Deep Indel Missense Programmable Library Engineering. Upload a gene "
            "FASTA, choose your mutation types, and download the oligo / primer "
            "library.\n\n"
            "*Headless mode: ambiguous ORFs must be resolved via `start:`/`end:` "
            "in the FASTA header or the ORF index field below.*"
        )

        with gr.Row():
            with gr.Column(scale=1):
                fasta = gr.File(label="Gene FASTA", file_types=[".fa", ".fasta", ".txt"])

                gr.Markdown("**Mutation types**")
                dms = gr.Checkbox(label="DMS (deep mutational scan)", value=True)
                dis = gr.Checkbox(label="DIS (domain insertion scan)", value=False)
                deletions_raw = gr.Textbox(
                    label="Deletions (nt, comma-separated; blank = none)", value=""
                )
                insertions_raw = gr.Textbox(
                    label="Insertions (nt sequences, comma-separated; blank = none)", value=""
                )

                with gr.Accordion("Core parameters", open=True):
                    oligo_len = gr.Number(label="Oligo length", value=230, precision=0)
                    overlap = gr.Number(label="Fragment overlap (nt)", value=4, precision=0)
                    fragment_len = gr.Textbox(label="Fragment length (or 'auto')", value="auto")
                    restriction_sequence = gr.Textbox(
                        label="Type IIS restriction sequence", value="CGTCTC(G)1/5"
                    )

                with gr.Accordion("Advanced", open=False):
                    avoid_sequence = gr.Textbox(
                        label="Sequences to avoid (comma-separated)", value="CGTCTC, GGTCTC"
                    )
                    handle = gr.Textbox(
                        label="Domain-insertion handle",
                        value="AGCGGGAGACCGGGGTCTCTGAGC",
                    )
                    codon_usage = gr.Dropdown(
                        label="Codon usage", choices=["human", "ecoli"], value="human"
                    )
                    make_double = gr.Checkbox(label="Make double mutants", value=False)
                    stop_codon = gr.Checkbox(label="Include stop codons", value=False)
                    include_synonymous = gr.Checkbox(label="Include synonymous", value=False)
                    maximize_nucleotide_change = gr.Checkbox(
                        label="Maximize nucleotide change", value=False
                    )
                    barcode_start = gr.Number(label="Barcode start", value=0, precision=0)
                    tm_low = gr.Number(label="Gene primer Tm (low)", value=58, precision=0)
                    tm_high = gr.Number(label="Gene primer Tm (high)", value=62, precision=0)
                    orf_index = gr.Textbox(label="Preferred ORF index (blank = auto)", value="")
                    random_seed = gr.Textbox(label="Random seed (blank = none)", value="1")

                run_btn = gr.Button("Design library", variant="primary")

            with gr.Column(scale=1):
                out_zip = gr.File(label="Results (zip)")
                status = gr.Textbox(label="Status", lines=20, max_lines=40)

        def _on_run(*vals):
            (
                fasta_v,
                dms_v,
                dis_v,
                dels_v,
                ins_v,
                oligo_v,
                overlap_v,
                frag_v,
                restr_v,
                avoid_v,
                handle_v,
                codon_v,
                double_v,
                stop_v,
                syn_v,
                maxnt_v,
                barcode_v,
                tmlo_v,
                tmhi_v,
                orf_v,
                seed_v,
            ) = vals
            return run_dimple_job(
                fasta_v,
                oligo_len=oligo_v,
                fragment_len=frag_v,
                overlap=overlap_v,
                dms=dms_v,
                dis=dis_v,
                deletions_raw=dels_v,
                insertions_raw=ins_v,
                make_double=double_v,
                stop_codon=stop_v,
                include_synonymous=syn_v,
                maximize_nucleotide_change=maxnt_v,
                restriction_sequence=restr_v,
                avoid_sequence=avoid_v,
                handle=handle_v,
                codon_usage=codon_v,
                barcode_start=barcode_v,
                tm_low=tmlo_v,
                tm_high=tmhi_v,
                orf_index=orf_v,
                random_seed=seed_v,
            )

        run_btn.click(
            _on_run,
            inputs=[
                fasta,
                dms,
                dis,
                deletions_raw,
                insertions_raw,
                oligo_len,
                overlap,
                fragment_len,
                restriction_sequence,
                avoid_sequence,
                handle,
                codon_usage,
                make_double,
                stop_codon,
                include_synonymous,
                maximize_nucleotide_change,
                barcode_start,
                tm_low,
                tm_high,
                orf_index,
                random_seed,
            ],
            outputs=[out_zip, status],
        )

    return demo


def _parse_auth(raw: str | None):
    """Parse ``DIMPLE_WEB_AUTH`` ("user:pass[,user2:pass2]") into Gradio auth."""
    if not raw:
        return None
    pairs = [tuple(p.split(":", 1)) for p in raw.split(",") if ":" in p]
    return pairs or None


def main() -> None:
    logging.basicConfig(level=logging.INFO)

    # Point Gradio's cache at an app-owned dir (not the shared system temp root)
    # so its served files are confined and reaped by delete_cache. Must be set
    # before launch reads it.
    GRADIO_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("GRADIO_TEMP_DIR", str(GRADIO_CACHE_DIR))

    _reap_old_jobs()  # clear anything stale from a previous process
    _start_reaper()

    # No broad allowed_paths: nothing outside Gradio's own cache is servable, so
    # there is no /gradio_api/file= route into a job dir. auth (if configured)
    # gives a true per-user boundary; without it, results are protected only by
    # unguessable paths + the TTL reaper.
    build_interface().launch(
        server_name=os.environ.get("DIMPLE_WEB_HOST", "127.0.0.1"),
        max_file_size=f"{MAX_UPLOAD_MB}mb",
        auth=_parse_auth(os.environ.get("DIMPLE_WEB_AUTH")),
    )


if __name__ == "__main__":
    main()
