"""HuggingFace Spaces / standalone entrypoint for the DIMPLE web app.

On HF Spaces this is the configured ``app_file``: Spaces runs ``python app.py``,
which sets up the isolation scaffolding (app-owned Gradio cache dir + job
reaper) and launches the interface. Gradio auto-binds to ``0.0.0.0:7860`` under
Spaces, so no ``server_name`` is forced here (unlike ``webapp.app.main`` which
defaults to localhost for local dev).

Locally you can run either this file (``python app.py``) or the module
(``python -m webapp.app``).
"""

from __future__ import annotations

import os

from webapp.app import (
    GRADIO_CACHE_DIR,
    MAX_UPLOAD_MB,
    _reap_old_jobs,
    _start_reaper,
    build_interface,
)

GRADIO_CACHE_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("GRADIO_TEMP_DIR", str(GRADIO_CACHE_DIR))

_reap_old_jobs()
_start_reaper()

demo = build_interface()

if __name__ == "__main__":
    demo.launch(max_file_size=f"{MAX_UPLOAD_MB}mb")
