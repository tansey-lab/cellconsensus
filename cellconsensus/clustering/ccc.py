"""CellConsensus Clustering (CCC): NN-based unsupervised cell typing."""
import numpy as np
from scipy import sparse
from scipy.stats import rankdata
from pynndescent import NNDescent


def build_nn_graph(adata, n_top_genes=2000, n_pcs=30, n_neighbors=20):
    """Build k-NN graph from PCA of normalized data.

    Modifies adata in-place (normalize_total, log1p, HVG, PCA).
    Returns (nn_indices, nn_distances) of shape (n_cells, n_neighbors).
    """
    import scanpy as sc

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)
    sc.tl.pca(adata, n_comps=n_pcs, use_highly_variable=True)

    X_pca = adata.obsm["X_pca"]
    index = NNDescent(X_pca, metric="cosine")
    nn_indices, nn_distances = index.query(X_pca, k=n_neighbors + 1)
    return nn_indices[:, 1:], nn_distances[:, 1:]


def quantile_normalize(adata):
    """Per-gene quantile normalization + per-cell L1 normalization.

    Returns sparse CSR matrix Q (n_cells x n_genes).
    Zeros stay zero. Nonzero values are ranked per gene and scaled to (0, 1].
    Each row is then L1-normalized (sums to 1) so cell-type scores are
    panel-size-invariant.
    """
    X = adata.X
    if not sparse.issparse(X):
        X = sparse.csc_matrix(X)
    else:
        X = X.tocsc()

    Q = X.copy().astype(np.float64)
    for j in range(X.shape[1]):
        nz_vals = Q.data[Q.indptr[j]:Q.indptr[j + 1]]
        if len(nz_vals) == 0:
            continue
        ranks = rankdata(nz_vals, method="average")
        Q.data[Q.indptr[j]:Q.indptr[j + 1]] = ranks / ranks.max()

    Q = Q.tocsr()
    for i in range(Q.shape[0]):
        row_data = Q.data[Q.indptr[i]:Q.indptr[i + 1]]
        norm = np.sum(row_data)
        if norm > 0:
            row_data /= norm

    return Q


def compute_scores(Q, R, ref_top_k=None):
    """L1-normalized cosine similarity S = Q @ R, optionally rescaled.

    With L1-normalized rows of Q and L1-normalized columns of R, the score is
    panel-size-invariant. The raw max of ``Q @ R`` is bounded by ``1 / k``
    where k is the per-column support of R, so multiplying by ``ref_top_k``
    lifts scores onto roughly ``[0, 1]`` and makes them comparable across
    panels of different sizes. Pass ``ref_top_k=None`` to skip the rescale
    (e.g. when assignment is by argmax and absolute scale does not matter).
    """
    S = Q @ R
    if sparse.issparse(S):
        S = S.toarray()
    if ref_top_k is not None:
        S = S * float(ref_top_k)
    return S
