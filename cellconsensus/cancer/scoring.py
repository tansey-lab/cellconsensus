"""Standalone scoring helpers (cancer + generic gene-set)."""
import pandas as pd

from ..consensus import build_reference_matrix
from ..clustering.ccc import compute_scores
from .cache import load_cancer_cache, validate_cancer_types


def _score_consensus(Q, var_names, consensus, keys, reduce_fn, ref_top_k=200):
    """Cosine similarity of Q against a custom consensus dict, then reduced.

    Parameters
    ----------
    Q : sparse matrix (n_cells x n_genes)
        Quantile-normalized cell matrix.
    var_names : array-like
        Gene names.
    consensus : dict[str, dict[str, float]]
        One inner gene->weight dict per key.
    keys : list of str
        Columns in the returned matrix; must be a subset of ``consensus``.
    reduce_fn : callable (S) -> S
        Reduces the per-cell score matrix (NN smoothing or cluster averaging).
    ref_top_k : int

    Returns
    -------
    ndarray (n_cells x len(keys))
    """
    R = build_reference_matrix(var_names, keys, ref_top_k, consensus=consensus)
    S = compute_scores(Q, R, ref_top_k=ref_top_k)
    return reduce_fn(S)


def score_cancer(Q, var_names, cancer_types, reduce_fn, ref_top_k=200):
    """Compute cancer scores for all cells.

    Parameters
    ----------
    Q : sparse matrix
        Quantile-normalized cell matrix.
    var_names : array-like
        Gene names.
    cancer_types : list of str or None
        Specific cancer types (e.g. ["melanoma", "breast_carcinoma"]).
        Pass None for pan-cancer (key "cancer").
        Import ``list_cancer_types()`` to see all valid keys.
    reduce_fn : callable (S) -> S
        Reduces the per-cell score matrix (NN smoothing or cluster averaging).
    ref_top_k : int

    Returns
    -------
    DataFrame with cancer scores per cell, one column per cancer type.

    Raises
    ------
    ValueError
        If ``cancer_types`` contains unknown keys.
    """
    validate_cancer_types(cancer_types)

    cache = load_cancer_cache()
    consensus = cache["consensus"]

    if cancer_types is None:
        cancer_keys = ["cancer"] if "cancer" in consensus else []
    else:
        cancer_keys = list(cancer_types)

    if not cancer_keys:
        return pd.DataFrame()

    S = _score_consensus(Q, var_names, consensus, cancer_keys, reduce_fn,
                         ref_top_k=ref_top_k)
    return pd.DataFrame(S, columns=cancer_keys)
