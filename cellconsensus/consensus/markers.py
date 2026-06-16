"""Load consensus marker caches and build reference matrices."""
import json
import pickle

import numpy as np
from pathlib import Path
from scipy import sparse

from .gene_alias import resolve_gene_names

DATA_DIR = Path(__file__).parent.parent / "data"
DEFAULT_CACHES = {
    1: DATA_DIR / "consensus_cache_meta.pkl",
    2: DATA_DIR / "consensus_cache_meta2.pkl",
    3: DATA_DIR / "consensus_cache_meta3.pkl",
}

_cl_to_meta = None


def load_cl_to_meta():
    """Load the cell type to meta group mapping."""
    global _cl_to_meta
    if _cl_to_meta is not None:
        return _cl_to_meta
    with open(DATA_DIR / "cl_to_meta.json") as f:
        _cl_to_meta = json.load(f)
    return _cl_to_meta


def load_meta_groups(level=1):
    """Load meta group definitions at a given level."""
    data = load_cl_to_meta()
    if level == 1:
        return data["meta_groups"]
    elif level == 2:
        return data["meta2_groups"]
    return data["meta3_groups"]


def load_meta_mapping(level=2):
    """Load child -> parent mapping for a given level."""
    data = load_cl_to_meta()
    if level == 2:
        return {v["meta2"]: v["meta"] for v in data["cl_to_meta"].values()}
    elif level == 3:
        return {v["meta3"]: v["meta2"] for v in data["cl_to_meta"].values()}
    return {}


def load_consensus(level=1):
    """Load consensus cache at a given level."""
    cache_path = DEFAULT_CACHES[level]
    with open(cache_path, "rb") as f:
        return pickle.load(f)


def get_meta_keys(level=1, include_cancer=False):
    """Get the list of meta keys to use for scoring.

    A key is returned only if it appears both in the consensus cache and in
    ``meta_groups`` for that level (so removing a key from ``meta_groups``
    is enough to drop it from level scoring, no cache rebuild needed). The
    catch-all ``"other"`` and any ``cancer*`` keys are always excluded.
    """
    cache = load_consensus(level)
    groups = load_meta_groups(level=level)
    all_cancer_keys = {k for k in cache["consensus"] if k.startswith("cancer")}
    exclude = {"other"} | all_cancer_keys
    return [k for k in sorted(cache["consensus"].keys())
            if k not in exclude and k in groups]


def load_cell_type(level=1):
    """Return the cell-type keys available for scoring at a given level.

    Use this to discover valid arguments for ``CellConsensus.predict_score``.
    Cancer keys are tracked separately — see ``list_cancer_types()``.

    Parameters
    ----------
    level : int
        1 (broad), 2 (subtypes), or 3 (fine-grained).

    Returns
    -------
    dict[str, str]
        Mapping ``key -> human-readable name`` for every cell type scorable
        at this level.

    Examples
    --------
    >>> from cellconsensus import load_cell_type
    >>> load_cell_type(level=1)
    {'adipocyte': 'adipocyte', 'b_plasma': 'B cell / plasma cell', ...}
    >>> list(load_cell_type(level=2).keys())[:3]
    ['absorptive_epithelial', 'adipocyte', 'astrocyte']
    """
    if level not in (1, 2, 3):
        raise ValueError(f"level must be 1, 2, or 3 (got {level}).")
    groups = load_meta_groups(level=level)
    keys = get_meta_keys(level=level)
    return {k: groups[k] for k in keys}


def build_reference_matrix(var_names, meta_keys, ref_top_k=200, level=1,
                           consensus=None):
    """Build reference matrix R (n_genes x n_meta_types), L1-normalized columns.

    Parameters
    ----------
    var_names : array-like
        Gene names from adata.var_names.
    meta_keys : list of str
        Meta type keys to include as columns.
    ref_top_k : int
        Number of top marker genes per meta type.
    level : int
        Which consensus cache to use (ignored if `consensus` is provided).
    consensus : dict or None
        Optional explicit consensus dict (e.g. from the cancer cache).

    Returns
    -------
    R : sparse CSC matrix (n_genes x len(meta_keys)) with each column summing
    to 1 — paired with L1-normalized rows of Q, the resulting Q @ R score is
    panel-size-invariant.
    """
    if consensus is None:
        consensus = load_consensus(level)["consensus"]

    gene_names = resolve_gene_names(var_names)
    gene_to_idx = {g: i for i, g in enumerate(gene_names)}
    n_genes = len(gene_names)

    row_idx, col_idx, vals = [], [], []
    for j, mk in enumerate(meta_keys):
        scores = consensus.get(mk, {})
        # Filter to genes in data, then take top-k
        ranked = sorted(scores.items(), key=lambda x: x[1], reverse=True)
        top_genes = [(g, s) for g, s in ranked if g in gene_to_idx][:ref_top_k]
        for g, s in top_genes:
            row_idx.append(gene_to_idx[g])
            col_idx.append(j)
            vals.append(np.sqrt(s))   # EXPERIMENT: sqrt-compress weights before L1-norm

    R = sparse.csc_matrix((vals, (row_idx, col_idx)),
                          shape=(n_genes, len(meta_keys)))
    # L1-normalize each column (sum to 1)
    for j in range(R.shape[1]):
        col_data = R.data[R.indptr[j]:R.indptr[j + 1]]
        norm = np.sum(col_data)
        if norm > 0:
            col_data /= norm

    return R
