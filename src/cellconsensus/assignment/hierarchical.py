"""Hierarchical assignment: level 1 -> 2 -> 3 cascading via argmax.

Both clustering modes share the same machinery and differ only in how the
per-cell score matrix is reduced before taking the argmax:

* ``ccc``         — smooth scores over the NN graph (``smooth``)
* ``precomputed`` — average scores within each given cluster (``cluster_average``)

A cell is labelled by the cell type with the highest (reduced) score. Within
each level-1 label, the same argmax is repeated over that type's level-2
children, then over the level-3 children of the chosen level-2 label, so the
label only ever gains precision down the hierarchy.
"""

import numpy as np
from scipy import sparse

from ..clustering.ccc import compute_scores
from ..consensus import (
    build_reference_matrix,
    load_consensus,
    load_meta_groups,
    load_meta_mapping,
)


# --------------------------------------------------------------- reducers
def _build_adjacency(nn_indices, n_cells):
    """Row-normalized adjacency matrix with self-loops."""
    row_idx, col_idx = [], []
    for i in range(n_cells):
        row_idx.append(i)
        col_idx.append(i)
        for j in nn_indices[i]:
            row_idx.append(i)
            col_idx.append(j)
    A = sparse.csr_matrix(
        (np.ones(len(row_idx)), (row_idx, col_idx)),
        shape=(n_cells, n_cells),
    )
    row_sums = np.array(A.sum(axis=1)).ravel()
    A = sparse.diags(1.0 / row_sums) @ A
    return A


def smooth(S, nn_indices, n_smooth):
    """Smooth a score matrix over the NN graph (iterative neighbor averaging)."""
    if n_smooth <= 0:
        return S
    A = _build_adjacency(nn_indices, S.shape[0])
    for _ in range(n_smooth):
        S = A @ S
    return S


def cluster_average(S, clusters):
    """Replace each cell's score row with the mean over its cluster."""
    out = np.zeros_like(S)
    for c in np.unique(clusters):
        mask = clusters == c
        out[mask] = S[mask].mean(axis=0)
    return out


def build_sub_nn(nn_indices, cell_idx, sub_n_neighbors):
    """Restrict the global NN graph to neighbors within ``cell_idx``.

    Returns an (len(cell_idx) x sub_n_neighbors) array of *local* indices,
    padded with self-loops when a cell has too few in-subset neighbors.
    """
    idx_set = set(int(i) for i in cell_idx)
    global_to_local = {int(g): loc for loc, g in enumerate(cell_idx)}
    sub_nn = np.zeros((len(cell_idx), sub_n_neighbors), dtype=int)
    for li, gi in enumerate(cell_idx):
        local = [
            global_to_local[int(nb)]
            for nb in nn_indices[gi]
            if int(nb) in idx_set and int(nb) != int(gi)
        ]
        while len(local) < sub_n_neighbors:
            local.append(li)
        sub_nn[li] = local[:sub_n_neighbors]
    return sub_nn


# --------------------------------------------------------------- assignment
def assign_argmax(S, keys):
    """Assign each cell to the highest-scoring type.

    Parameters
    ----------
    S : ndarray (n_cells, n_keys)
        (Reduced) score matrix.
    keys : list of str
        Type keys, one per column of S.

    Returns
    -------
    labels : ndarray of object
        Best key per cell; "unassigned" where every score is <= 0.
    scores : ndarray of float
        Winning score per cell.
    """
    keys = np.asarray(keys, dtype=object)
    best = np.argmax(S, axis=1)
    scores = S[np.arange(S.shape[0]), best]
    labels = keys[best].astype(object)
    labels[scores <= 0] = "unassigned"
    return labels, scores


def run_sublevel(
    Q,
    cell_labels,
    child_level,
    var_names,
    reduce_fn,
    ref_top_k=200,
    skip_labels=("unassigned",),
    verbose=True,
):
    """Refine parent-level labels into their child-level subtypes.

    For each parent label, score its cells against that parent's sibling
    subtypes, reduce the scores with ``reduce_fn``, and re-assign by argmax.
    Cells whose parent has fewer than two subtypes (or is in ``skip_labels``)
    keep their parent label.

    Parameters
    ----------
    Q : sparse matrix (n_cells x n_genes)
        Quantile-normalized cell matrix.
    cell_labels : ndarray of str
        Parent-level assignment per cell.
    child_level : int
        2 or 3.
    var_names : array-like
        Gene names.
    reduce_fn : callable (S, cell_idx) -> S
        Reduces a subset's score matrix (smoothing or cluster averaging).
    ref_top_k : int
        Top consensus markers per subtype.
    skip_labels : tuple of str
        Parent labels to leave untouched.

    Returns
    -------
    (labels, scores) : ndarray of object, ndarray of float
    """
    n_cells = len(cell_labels)
    child_mapping = load_meta_mapping(level=child_level)
    child_groups = load_meta_groups(level=child_level)
    child_consensus = load_consensus(level=child_level)["consensus"]

    new_labels = np.array(cell_labels, dtype=object)
    new_scores = np.zeros(n_cells)

    skip_set = set(skip_labels)
    parents = [p for p in sorted(set(cell_labels)) if p not in skip_set]

    for parent in parents:
        siblings = sorted(
            [
                k
                for k in child_consensus
                if child_mapping.get(k) == parent
                and k != parent
                and not k.endswith("_other")
            ]
        )
        if len(siblings) < 2:
            continue  # keep parent label, no subtypes to choose from

        cell_idx = np.where(cell_labels == parent)[0]
        if verbose:
            print(f"  {parent}: {len(cell_idx)} cells, {len(siblings)} subtypes")

        R = build_reference_matrix(var_names, siblings, ref_top_k, level=child_level)
        S = compute_scores(Q[cell_idx], R)
        S = reduce_fn(S, cell_idx)
        labels, scores = assign_argmax(S, siblings)
        new_labels[cell_idx] = labels
        new_scores[cell_idx] = scores

    if verbose:
        print(f"\n  Level {child_level} summary:")
        for mk in sorted(set(new_labels)):
            count = (new_labels == mk).sum()
            if count > 0:
                avg = new_scores[new_labels == mk].mean()
                print(
                    f"    {child_groups.get(mk, mk):<35} {count:>6} cells  "
                    f"(avg score: {avg:.4f})"
                )

    return new_labels, new_scores
