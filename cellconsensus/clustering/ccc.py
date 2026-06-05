"""CellConsensus Clustering (CCC): NN-based unsupervised cell typing."""
import numpy as np
from scipy import sparse
from scipy.stats import rankdata
from pynndescent import NNDescent


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


def _quantile_axis(Q):
    """Zero-aware quantile-normalize the nonzero values of each CSR row.

    Zeros stay 0; nonzeros are ranked within the full row and scaled to
    (frac_zero, 1] — the ranks start right after the zero mass. For column
    normalization, pass the transpose (CSR of the matrix's T).
    """
    n = Q.shape[1]
    for i in range(Q.shape[0]):
        s, e = Q.indptr[i], Q.indptr[i + 1]
        seg = Q.data[s:e]
        if len(seg) == 0:
            continue
        r = rankdata(seg, method="average")
        Q.data[s:e] = (n - len(seg) + r) / n
    return Q


def double_quantile_normalize(X):
    """Quantile-normalize rows -> quantile-normalize cols -> L1-normalize rows.

    Count-free, scanpy-free normalization for the score-space graph. Both
    quantile passes are zero-aware (sparsity preserved). The final per-cell L1
    (each row sums to 1) pairs with the L1-normalized reference columns so
    ``Q @ R`` scores are panel-size-invariant. Returns CSR.
    """
    X = X.tocsr().astype(np.float64).copy() if sparse.issparse(X) \
        else sparse.csr_matrix(X, dtype=np.float64)
    Qr = _quantile_axis(X)                      # quantile-norm rows
    Qc = _quantile_axis(Qr.T.tocsr())           # quantile-norm cols (rows of T)
    Q = Qc.T.tocsr()
    rs = np.asarray(Q.sum(1)).ravel()           # L1 rows (per-cell, sums to 1)
    rs[rs == 0] = 1.0
    return (sparse.diags(1.0 / rs) @ Q).tocsr()


def build_score_graph(adata, R, n_neighbors=20):
    """Cell-type score-space NN graph: double-QN, project onto R, cosine kNN.

    ``Q = double_quantile_normalize(counts)``; ``A = Q @ R`` embeds each cell in
    cell-type score space (one coordinate per reference program); the cosine
    kNN graph is built on ``A``. No normalize_total / log1p / HVG / PCA.

    Returns (nn_indices, Q, A): nn_indices (n_cells x n_neighbors), the
    double-QN matrix Q (reused for cell-type scoring), and the embedding A
    (n_cells x n_programs).
    """
    Q = double_quantile_normalize(adata.X)
    A = Q @ R
    A = np.asarray(A.todense()) if sparse.issparse(A) \
        else np.ascontiguousarray(A)
    index = NNDescent(A, metric="cosine")
    nn_indices, _ = index.query(A, k=n_neighbors + 1)
    return nn_indices[:, 1:], Q, A


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
