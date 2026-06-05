"""Standalone scoring helper: cosine similarity against a consensus dict."""
from ..consensus import build_reference_matrix
from ..clustering.ccc import compute_scores


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
    S = compute_scores(Q, R)
    return reduce_fn(S)

