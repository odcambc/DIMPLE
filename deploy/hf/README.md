---
title: DIMPLE
emoji: 🧬
colorFrom: blue
colorTo: green
sdk: gradio
app_file: app.py
python_version: "3.12"
pinned: false
short_description: Design oligo libraries for deep mutational scanning
---

# DIMPLE

Deep Indel Missense Programmable Library Engineering — a browser front-end for
designing oligonucleotide libraries for deep mutational scanning. Upload a gene
FASTA, choose your mutation types (DMS / DIS / insertions / deletions), and
download the oligo + primer FASTAs and variant / mutation CSVs.

This Space wraps the same pipeline the CLI and Tk GUI drive
(`DIMPLE.runner.build_runtime_config` + `run_pipeline`); see the
[DIMPLE repository](https://github.com/coywil26/DIMPLE) for the science and the
full tool.

> **Public demo.** Runs are isolated per request (isolated temp dirs, the raw
> upload is deleted as soon as the pipeline consumes it, and results are reaped
> on a TTL), but there is no per-user login — anyone with the URL can use it.
> Don't upload anything you wouldn't want a stranger to run.

This README's YAML header is the HuggingFace Spaces config. Notes:

- `app_file: app.py` is the repo-root entrypoint.
- `sdk_version` is intentionally unpinned so the Space uses HF's current Gradio
  (the app needs `delete_cache` ≥ 4.19 and `max_file_size` ≥ 4.x). Pin it here
  if you want a reproducible build.
