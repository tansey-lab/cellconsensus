"""CellConsensus Clustering (CCC): NN-based unsupervised cell typing."""
import numpy as np
from scipy import sparse
from scipy.stats import rankdata
from pynndescent import NNDescent


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


def compute_scores(Q, R):
    """L1-normalized cosine similarity S = Q @ R.

    With L1-normalized rows of Q and L1-normalized columns of R, each entry is
    the fraction of the cell's expression mass captured by that reference
    column's markers — already panel-size-invariant, so the same scale holds
    whether R has one column (a gene set) or many (the full cell-type panel).
    """
    S = Q @ R
    if sparse.issparse(S):
        S = S.toarray()
    return S
