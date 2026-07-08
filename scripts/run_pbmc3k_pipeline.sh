#!/usr/bin/env bash
# Integration test: run the cellconsensus Nextflow pipeline on scanpy's pbmc3k dataset.
#
# Uses the `local` profile so processes run on the host against the current
# cellconsensus source tree — no container involved.
#
# Assumes cellconsensus (with its CLI) and scanpy are already installed in the
# active Python environment.
#
# Usage:
#   scripts/run_pbmc3k_pipeline.sh [OUTDIR]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
OUTDIR="${1:-${REPO_ROOT}/tests/integration/pbmc3k_results}"
WORKDIR="${REPO_ROOT}/tests/integration/pbmc3k_work"
DATA_DIR="${REPO_ROOT}/tests/integration/pbmc3k_data"
H5AD_PATH="${DATA_DIR}/pbmc3k.h5ad"
SAMPLESHEET="${DATA_DIR}/samplesheet.csv"

mkdir -p "${DATA_DIR}" "${OUTDIR}" "${WORKDIR}"

# 1. Materialize the pbmc3k dataset via scanpy.
if [[ ! -f "${H5AD_PATH}" ]]; then
    echo "[pbmc3k] Downloading pbmc3k via scanpy → ${H5AD_PATH}"
    python - <<PY
import scanpy as sc

adata = sc.datasets.pbmc3k()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.write_h5ad("${H5AD_PATH}")
print(f"[pbmc3k] wrote {adata.n_obs} cells x {adata.n_vars} genes")
PY
else
    echo "[pbmc3k] Reusing cached ${H5AD_PATH}"
fi

cat > "${SAMPLESHEET}" <<CSV
sample_id,adata_path
pbmc3k,${H5AD_PATH}
CSV

echo "[pbmc3k] Samplesheet:"
cat "${SAMPLESHEET}"

echo "[pbmc3k] Launching pipeline (local profile) → ${OUTDIR}"
nextflow run "${REPO_ROOT}/nextflow/main.nf" \
    -profile local \
    -work-dir "${WORKDIR}" \
    --input "${SAMPLESHEET}" \
    --outdir "${OUTDIR}"

echo "[pbmc3k] Verifying outputs..."
for f in \
    "${OUTDIR}/annotate/pbmc3k_annotate/predictions.csv" \
    "${OUTDIR}/annotate/pbmc3k_annotate/annotated.h5ad"
do
    if [[ ! -s "${f}" ]]; then
        echo "[pbmc3k] MISSING expected output: ${f}" >&2
        exit 1
    fi
    echo "[pbmc3k]   ok: ${f}"
done

echo "[pbmc3k] Integration test passed."
