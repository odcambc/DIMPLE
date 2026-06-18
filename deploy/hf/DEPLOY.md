# Deploy the DIMPLE web app to HuggingFace Spaces

Target: a **public, no-login** Gradio Space. The Space is its own git repo; we
push a curated subset of the project to it (app + code + a Space-specific
`requirements.txt`/`README.md`). The main repo is uv/`pyproject.toml`-managed and
ships no flat `requirements.txt`, so the Space carries its own (core deps +
`gradio`) for HF's pip-based builder.

## What ships to the Space

At the Space repo root:

- `app.py` — repo-root entrypoint (from the DIMPLE root).
- `requirements.txt` — from `deploy/hf/requirements.txt` (adds `gradio`).
- `README.md` — from `deploy/hf/README.md` (carries the HF Spaces YAML config).
- `DIMPLE/` and `webapp/` — the packages the app imports.

Not shipped: `tests/`, goldens, tracking files. (The app's *advisory* output
consistency check imports `tests.regression._helpers`; on the Space that import
is simply absent and the check is skipped — the run still completes and zips.)

No Git LFS needed — nothing large ships.

## One-time setup

```bash
pip install -U "huggingface_hub[cli]"
hf auth login        # paste a WRITE token from https://huggingface.co/settings/tokens
```

## Create + push the Space

```bash
# 1. Create the Space (Gradio SDK). Pick your own <user>/<name>.
hf repo create dimple-demo --repo-type space --space_sdk gradio
#    (or: huggingface.co -> New -> Space -> Gradio SDK)

# 2. Clone the new (near-empty) Space repo OUTSIDE the DIMPLE checkout.
git clone https://huggingface.co/spaces/<user>/dimple-demo /tmp/dimple-space

# 3. Assemble the Space contents from the DIMPLE checkout.
DIMPLE_DIR="$HOME/Projects/DIMPLE"          # adjust to your path
cd /tmp/dimple-space
cp -R "$DIMPLE_DIR/DIMPLE" "$DIMPLE_DIR/webapp" "$DIMPLE_DIR/app.py" .
cp "$DIMPLE_DIR/deploy/hf/requirements.txt" requirements.txt
cp "$DIMPLE_DIR/deploy/hf/README.md" README.md

# 4. Push — HF builds and serves automatically.
git add -A
git commit -m "Deploy DIMPLE Gradio demo"
git push
```

The first build installs deps (a couple of minutes); then the app is live at
`https://huggingface.co/spaces/<user>/dimple-demo`. Watch the build/run logs in
the Space's **Logs** tab.

## Smoke-test the live Space

Upload `tests/data/Kir.fa`, leave **DMS** checked, click **Design library**.
Expect a `DIMPLE_results.zip` (7 files) and a status log ending in the output
file list. (Kir DMS is ~8k oligos and runs in well under a minute on the free
CPU tier.)

## Updating later

Re-run step 3's copies in the Space clone, then `git add -A && git commit && git
push`. To pin Gradio for reproducible builds, add `sdk_version: "<x.y.z>"` to
the README front-matter and match it in `requirements.txt`.

## If you'd rather not maintain a separate Space repo

Add the Space as a remote on a throwaway deploy branch of the main repo and push
that — but the main repo ships no root `requirements.txt` (deps live in
`pyproject.toml` / `uv.lock`), and its `README.md` has no HF front-matter. HF's
Gradio builder wants a flat `requirements.txt`, so you'd have to add one (core
deps + `gradio`) and prepend the YAML front-matter to the root README — i.e.
re-create exactly what `deploy/hf/` already holds. The curated-copy flow above
avoids polluting the repo root with deploy-only files.
```
