"""CellConsensus: hierarchical cell type annotation."""
import pickle
import warnings

import numpy as np
import pandas as pd
from scipy import sparse

from .consensus import (
    get_meta_keys,
    build_reference_matrix,
    load_meta_groups,
    load_consensus,
    load_cell_type,
)
from .clustering.ccc import (compute_scores, double_quantile_normalize,
                             build_score_graph)
from .assignment import (
    smooth,
    cluster_average,
    build_sub_nn,
    assign_argmax,
    run_sublevel,
)
from .cancer import (
    _score_consensus,
    load_cancer_cache,
    validate_cancer_types,
    cancer_key_to_ncit,
    cancer_key_to_name,
    is_cancer_key,
)


class CellConsensus:
    """Hierarchical unsupervised cell type annotation.

    CellConsensus assigns cell types at three levels of granularity using
    consensus marker genes aggregated from multiple databases and publications.
    A cell (or cluster) is labelled by the cell type with the highest score;
    the label is then refined into that type's level-2 children and finally its
    level-3 children, so precision only ever increases down the hierarchy.

    Two clustering strategies are available, differing only in how per-cell
    scores are pooled before the argmax:

    **"ccc" (CellConsensus Clustering)** — fully unsupervised, per-cell:
        1. Build a k-NN graph in PCA space (pynndescent, cosine distance)
        2. Quantile-normalize each gene across cells, then L1-normalize cells
        3. Score each cell against consensus marker references (cosine)
        4. Smooth scores over the NN graph (iterative neighbor averaging)
        5. Assign each cell to its argmax type, then refine within levels 2/3
           using progressively smaller NN graphs (20 -> 10 -> 5 neighbors)

    **"precomputed"** — user provides cluster labels:
        Same scoring, but instead of NN smoothing the scores are averaged
        across all cells in each given cluster. Every cluster therefore gets a
        single label per level; clusters are never split.

    Parameters
    ----------
    clustering : str
        "ccc" or "precomputed".
    n_neighbors : int
        Number of nearest neighbors for level 1 (ccc mode).
    n_neighbors_lvl2 : int
        Number of nearest neighbors for level 2 refinement (ccc mode).
    n_neighbors_lvl3 : int
        Number of nearest neighbors for level 3 refinement (ccc mode).
    n_smooth : int
        Number of smoothing iterations over the NN graph (ccc mode).
    ref_top_k : int
        Number of top consensus markers per cell type.
    graph_level : int
        Consensus level whose programs define the score-space embedding the
        ccc kNN graph is built on (ccc mode).
    cluster_key : str
        Column in adata.obs holding cluster labels (precomputed mode only).

    Examples
    --------
    >>> cc = CellConsensus()
    >>> cc.fit(adata)
    >>> labels = cc.predict(level=2)
    >>> cancer = cc.predict_score("cancer")
    >>> # Bring your own clusters
    >>> cc = CellConsensus(clustering="precomputed", cluster_key="leiden")
    >>> cc.fit(adata)
    >>> labels = cc.predict(level=3)
    >>> # Include cancer as a class during fit
    >>> cc.fit(adata, include_cancer=True)
    >>> # User-supplied gene signature
    >>> sig = cc.predict_gene_set(["CD3D", "CD3E"], name="t_sig")
    >>> # Persist a fit
    >>> cc.save("model.pkl")
    >>> cc2 = CellConsensus.load("model.pkl")
    """

    def __init__(self, clustering="ccc", n_neighbors=20,
                 n_neighbors_lvl2=10, n_neighbors_lvl3=5,
                 n_smooth=1, ref_top_k=200, graph_level=3,
                 cluster_key=None):
        if clustering not in ("ccc", "precomputed"):
            raise ValueError(f"Unknown clustering: {clustering}. "
                             f"Use 'ccc' or 'precomputed'.")
        self.clustering = clustering
        self.n_neighbors = n_neighbors
        self.n_neighbors_lvl2 = n_neighbors_lvl2
        self.n_neighbors_lvl3 = n_neighbors_lvl3
        self.n_smooth = n_smooth
        self.ref_top_k = ref_top_k
        self.graph_level = graph_level
        self.cluster_key = cluster_key

        # Populated by fit()
        self.obs_names_ = None
        self.var_names_ = None
        self.nn_indices_ = None      # ccc only
        self.clusters_ = None        # precomputed only
        self.Q_ = None
        self.S_ = None               # level-1 (reduced) score matrix
        self.meta_keys_ = None       # level-1 keys (+ cancer keys)

        # Cancer-in-fit metadata (set by fit when include_cancer=True)
        self.include_cancer_ = False
        self.cancer_keys_ = None

        # level -> (labels, scores)
        self._predictions = {}

    # ------------------------------------------------------------------ fit
    def fit(self, adata, include_cancer=False, cancer_types=None, verbose=True):
        """Fit the model on an AnnData object.

        Assigns cell types at all three levels. For ccc mode this builds a NN
        graph and smooths per-cell scores; for precomputed mode it averages
        scores within each cluster from ``adata.obs[cluster_key]``.

        Parameters
        ----------
        adata : AnnData
            Raw count matrix. Modified in-place for normalization/PCA.
        include_cancer : bool
            If True, add cancer key(s) as extra columns of the level-1 score
            matrix. Cells/clusters whose cancer score wins are labeled with the
            cancer key (e.g. "cancer"); refinement skips them.
        cancer_types : list of str or None
            Cancer keys to include when ``include_cancer=True``. None defaults
            to the pan-cancer key "cancer". See ``list_cancer_types()``.
        verbose : bool
            Print progress.

        Returns
        -------
        self
        """
        self.obs_names_ = np.asarray(adata.obs_names)
        self.var_names_ = np.asarray(adata.var_names)
        self._predictions = {}

        self.include_cancer_ = bool(include_cancer)
        if self.include_cancer_:
            validate_cancer_types(cancer_types)
            cache = load_cancer_cache()
            if cancer_types is None:
                self.cancer_keys_ = (["cancer"]
                                     if "cancer" in cache["consensus"] else [])
            else:
                self.cancer_keys_ = list(cancer_types)
            if not self.cancer_keys_:
                # No cancer keys available — silently disable
                self.include_cancer_ = False
                self.cancer_keys_ = None
        else:
            self.cancer_keys_ = None

        # Level-1 cell-type keys (cancer keys are appended in _assign_all_levels)
        self.meta_keys_ = list(get_meta_keys(level=1))

        if self.clustering == "ccc":
            self._prepare_ccc(adata, verbose)
        else:
            self._prepare_precomputed(adata, verbose)

        self._assign_all_levels(verbose)
        return self

    def _prepare_ccc(self, adata, verbose):
        """Build the score-space NN graph and double-quantile-normalize.

        Cells are double-quantile-normalized (zero-aware QN rows -> QN cols ->
        L1 rows) and projected onto the level-``graph_level`` reference, giving
        each cell a coordinate per cell-type program. The cosine kNN graph is
        built in that score space — no normalize_total / log1p / HVG / PCA.
        """
        keys = list(load_consensus(self.graph_level)["consensus"].keys())
        R = build_reference_matrix(np.asarray(adata.var_names), keys,
                                   self.ref_top_k, level=self.graph_level)
        if verbose:
            print(f"Score-space {self.n_neighbors}-NN graph "
                  f"({adata.n_obs} cells, {len(keys)} level-"
                  f"{self.graph_level} programs, double quantile-norm)...")
        nn_indices, Q, A = build_score_graph(adata, R, self.n_neighbors)
        self.nn_indices_ = nn_indices
        self.clusters_ = None
        self.Q_ = Q
        adata.obsm["X_cc_embed"] = A

    def _prepare_precomputed(self, adata, verbose):
        """Validate clusters and double-quantile-normalize (precomputed mode).

        Uses the same data processing as ccc mode (double quantile-norm with
        L1-normalized rows); the only difference is that scores are reduced by
        averaging within the supplied clusters rather than NN smoothing.
        """
        if self.cluster_key is None:
            raise ValueError("cluster_key must be set for precomputed mode.")
        if self.cluster_key not in adata.obs.columns:
            raise ValueError(
                f"Cluster key '{self.cluster_key}' not found in adata.obs. "
                f"Available columns: {list(adata.obs.columns)}"
            )
        self.clusters_ = adata.obs[self.cluster_key].astype(str).values
        self.nn_indices_ = None

        if verbose:
            n_clusters = len(np.unique(self.clusters_))
            print(f"Using {n_clusters} precomputed clusters "
                  f"({self.cluster_key}); double quantile-norm...")
        self.Q_ = double_quantile_normalize(adata.X)

    # ----------------------------------------------------------- reducers
    def _level1_reduce(self, S):
        """Reduce a full-matrix score: smooth (ccc) or cluster-average."""
        if self.clustering == "ccc":
            return smooth(S, self.nn_indices_, self.n_smooth)
        return cluster_average(S, self.clusters_)

    def _sublevel_reduce_fn(self, sub_n_neighbors):
        """Build a reduce_fn(S, cell_idx) for level-2/3 refinement."""
        if self.clustering == "ccc":
            nn_indices = self.nn_indices_
            n_smooth = self.n_smooth

            def fn(S, cell_idx):
                sub_nn = build_sub_nn(nn_indices, cell_idx, sub_n_neighbors)
                return smooth(S, sub_nn, n_smooth)
            return fn

        clusters = self.clusters_

        def fn(S, cell_idx):
            return cluster_average(S, clusters[cell_idx])
        return fn

    def _assign_all_levels(self, verbose):
        """Score and assign levels 1, 2, 3 with the mode's reducer."""
        meta_keys = list(self.meta_keys_)  # cell-type keys only at this point
        cancer_keys = list(self.cancer_keys_) if self.include_cancer_ else []
        all_keys = meta_keys + cancer_keys

        if verbose:
            extra = f" + {len(cancer_keys)} cancer" if cancer_keys else ""
            print(f"Scoring level 1 ({len(meta_keys)} types{extra})...")

        R = build_reference_matrix(
            self.var_names_, meta_keys, self.ref_top_k, level=1)
        if cancer_keys:
            cache = load_cancer_cache()
            R_cancer = build_reference_matrix(
                self.var_names_, cancer_keys, self.ref_top_k,
                consensus=cache["consensus"])
            R = sparse.hstack([R, R_cancer]).tocsc()

        S = compute_scores(self.Q_, R)
        S = self._level1_reduce(S)
        self.S_ = S
        self.meta_keys_ = all_keys

        labels1, scores1 = assign_argmax(S, all_keys)
        self._predictions[1] = (labels1, scores1)
        if verbose:
            groups = load_meta_groups(level=1)
            for mk in all_keys:
                count = int((labels1 == mk).sum())
                if count > 0:
                    display = (groups.get(mk, mk) if mk not in cancer_keys
                               else cancer_key_to_name(mk))
                    print(f"  {display:<35} {count:>6} cells")

        # Refinement skips "unassigned" and all cancer keys (no taxonomy below).
        skip_labels = tuple(["unassigned"] + cancer_keys)

        if verbose:
            print("\n--- Level 2 ---")
        labels2, scores2 = run_sublevel(
            self.Q_, labels1, 2, self.var_names_,
            self._sublevel_reduce_fn(self.n_neighbors_lvl2),
            ref_top_k=self.ref_top_k, skip_labels=skip_labels, verbose=verbose)
        self._predictions[2] = (labels2, scores2)

        if verbose:
            print("\n--- Level 3 ---")
        labels3, scores3 = run_sublevel(
            self.Q_, labels2, 3, self.var_names_,
            self._sublevel_reduce_fn(self.n_neighbors_lvl3),
            ref_top_k=self.ref_top_k, skip_labels=skip_labels, verbose=verbose)
        self._predictions[3] = (labels3, scores3)

    # -------------------------------------------------------------- predict
    def predict(self, level=3, output="name", include_cancer=True):
        """Return cell type predictions.

        Parameters
        ----------
        level : int
            Granularity level.
            1 = broad (T cell, myeloid, ...)
            2 = subtypes (CD4 T, NK, monocyte, ...)
            3 = fine-grained (naive CD4, effector CD8, cDC2, ...)
        output : str
            "name"  - human-readable names (default)
            "cl_id" - Cell Ontology IDs (e.g. CL:0000624) or NCIT codes for
                      cancer-labeled cells (e.g. NCIT:C3262 for "cancer")
            "key"   - internal meta keys (e.g. cd4_t, naive_cd4, cancer)
        include_cancer : bool
            If True (default), return cancer labels as-is for cells assigned to
            a cancer class during fit. If False and the model was fit with
            include_cancer=True, those cells are re-assigned post-hoc to their
            best non-cancer column. No-op when no cancer keys were fit.

        Returns
        -------
        pandas Series indexed by obs_names.
        """
        from .consensus.cl_mapping import meta_to_cl, meta_to_name

        self._check_fitted()
        if level not in self._predictions:
            raise RuntimeError(f"Level {level} predictions not available.")

        cancer_keys = set(self.cancer_keys_ or [])
        raw_labels = np.asarray(self._predictions[level][0])

        if cancer_keys and not include_cancer:
            raw_labels = self._strip_cancer_labels(raw_labels)

        if output == "key":
            values = raw_labels
        elif output == "cl_id":
            values = np.array([
                cancer_key_to_ncit(mk) if mk in cancer_keys else meta_to_cl(mk)
                for mk in raw_labels
            ])
        else:  # name
            values = np.array([
                cancer_key_to_name(mk) if mk in cancer_keys else meta_to_name(mk)
                for mk in raw_labels
            ])

        return pd.Series(values, index=self.obs_names_, name="cellconsensus")

    def _strip_cancer_labels(self, raw_labels):
        """Re-assign cancer-labeled cells to their best non-cancer column in S_.

        Returns a new array; does not mutate ``self._predictions``.
        """
        cancer_keys = set(self.cancer_keys_ or [])
        if not cancer_keys:
            return raw_labels

        meta_keys = self.meta_keys_
        non_cancer_idx = [j for j, mk in enumerate(meta_keys)
                          if mk not in cancer_keys]
        if not non_cancer_idx:
            return raw_labels
        sub_keys = [meta_keys[j] for j in non_cancer_idx]
        sub_S = self.S_[:, non_cancer_idx]

        out = np.array(raw_labels, dtype=object)
        cancer_mask = np.array([lbl in cancer_keys for lbl in raw_labels])
        if cancer_mask.any():
            best = np.argmax(sub_S[cancer_mask], axis=1)
            out[cancer_mask] = np.array(sub_keys, dtype=object)[best]
        return out

    def predict_score(self, cell_types, level=1, smooth=True, verbose=False):
        """Score cells against any cell type(s) or cancer type(s).

        Pass any consensus key (cell type at
        the chosen ``level``, or any cancer key) and get the per-cell score,
        reduced the same way as fit (NN smoothing for ccc, cluster averaging
        for precomputed).

        Parameters
        ----------
        cell_types : str or list of str
            One or more keys. Cell-type keys are looked up in the level-N
            consensus cache; cancer keys are auto-detected and pulled from
            the cancer cache (the ``level`` argument is irrelevant for
            those). See ``load_cell_type(level=N)`` and ``list_cancer_types()``
            for valid keys.
        level : int
            1, 2, or 3 — which cell-type cache to draw non-cancer references
            from. Ignored for cancer keys.
        smooth : bool
            If True (default), reduce scores like fit does. If False, return
            raw per-cell cosine scores without smoothing.
        verbose : bool

        Returns
        -------
        pandas DataFrame (n_cells x len(cell_types)), indexed by obs_names.

        Examples
        --------
        >>> cc.predict_score("t_cell")                       # L1 T-cell score
        >>> cc.predict_score(["nk", "monocyte"], level=2)    # two L2 scores
        >>> cc.predict_score(["t_cell", "melanoma"])         # cell type + cancer
        """
        self._check_fitted()

        if isinstance(cell_types, str):
            keys = [cell_types]
        else:
            keys = list(cell_types)
        if not keys:
            raise ValueError("`cell_types` must be non-empty.")

        ct_consensus = load_consensus(level=level)["consensus"]
        cancer_consensus = load_cancer_cache()["consensus"]

        merged = {}
        unknown = []
        for k in keys:
            if is_cancer_key(k):
                merged[k] = cancer_consensus[k]
            elif k in ct_consensus:
                merged[k] = ct_consensus[k]
            else:
                unknown.append(k)

        if unknown:
            raise ValueError(
                f"Unknown key(s) at level={level}: {unknown}. "
                f"See `load_cell_type(level={level})` for cell types and "
                f"`list_cancer_types()` for cancer keys."
            )

        reduce_fn = self._level1_reduce if smooth else (lambda S: S)
        if verbose:
            print(f"Scoring {len(keys)} type(s) at level={level} "
                  f"({self.clustering}, reduce={'yes' if smooth else 'no'})...")
        S = _score_consensus(
            self.Q_, self.var_names_, merged, keys,
            reduce_fn, ref_top_k=self.ref_top_k,
        )
        return pd.DataFrame(S, columns=keys, index=self.obs_names_)

    def predict_gene_set(self, genes, weights=None, name="gene_set",
                         smooth=True, verbose=False):
        """Score cells against a user-supplied gene signature.

        Parameters
        ----------
        genes : list of str
            Gene symbols. Resolved through the gene alias map. Genes absent
            from ``var_names_`` after resolution are silently dropped.
        weights : list of float or dict[str, float] or None
            Per-gene weights. None => uniform 1.0. If a list, must match
            ``genes`` 1-to-1. If a dict, keys must be a subset of ``genes``;
            missing entries default to 1.0. Negative weights are allowed.
        name : str
            Column name in the returned DataFrame.
        smooth : bool
            Reduce scores like fit does (NN smoothing for ccc, cluster
            averaging for precomputed). If False, return raw per-cell scores.
        verbose : bool

        Returns
        -------
        pandas DataFrame with a single column ``name``, indexed by obs_names.
        """
        self._check_fitted()

        gene_list = list(genes)
        if len(gene_list) == 0:
            raise ValueError("`genes` must be non-empty.")

        if weights is None:
            weight_map = {g: 1.0 for g in gene_list}
        elif isinstance(weights, dict):
            weight_map = {g: float(weights.get(g, 1.0)) for g in gene_list}
        else:
            weights_list = list(weights)
            if len(weights_list) != len(gene_list):
                raise ValueError(
                    f"`weights` length ({len(weights_list)}) does not match "
                    f"`genes` length ({len(gene_list)})."
                )
            weight_map = {g: float(w) for g, w in zip(gene_list, weights_list)}

        consensus = {name: weight_map}
        reduce_fn = self._level1_reduce if smooth else (lambda S: S)

        if verbose:
            print(f"Scoring gene set '{name}' ({len(gene_list)} genes, "
                  f"reduce={'yes' if smooth else 'no'})...")
        S = _score_consensus(
            self.Q_, self.var_names_, consensus, [name],
            reduce_fn, ref_top_k=len(gene_list),
        )
        if S.shape[1] == 0 or np.allclose(S, 0):
            warnings.warn(
                f"Gene set '{name}': no input genes matched var_names "
                "after alias resolution; returning zeros."
            )
        return pd.DataFrame(S, columns=[name], index=self.obs_names_)

    def score_matrix(self, level=1):
        """Return the level-1 score matrix.

        Parameters
        ----------
        level : int

        Returns
        -------
        DataFrame (n_cells x n_types).
        """
        self._check_fitted()
        if level == 1:
            return pd.DataFrame(
                self.S_, columns=self.meta_keys_, index=self.obs_names_)
        raise NotImplementedError("Score matrix for level > 1 not yet supported.")

    # ---------------------------------------------------------- persistence
    def save(self, path):
        """Persist fitted state to ``path`` (pickle).

        Saves what is required to call ``predict``, ``predict_score``, and
        ``predict_gene_set`` after ``load()``. The original AnnData is not
        stored (users keep that separately).
        """
        self._check_fitted()
        from . import __version__ as _cc_version
        state = {
            "version": _cc_version,
            "hyperparameters": {
                "clustering": self.clustering,
                "n_neighbors": self.n_neighbors,
                "n_neighbors_lvl2": self.n_neighbors_lvl2,
                "n_neighbors_lvl3": self.n_neighbors_lvl3,
                "n_smooth": self.n_smooth,
                "ref_top_k": self.ref_top_k,
                "graph_level": self.graph_level,
                "cluster_key": self.cluster_key,
            },
            "fitted": {
                "obs_names": self.obs_names_,
                "var_names": self.var_names_,
                "nn_indices": self.nn_indices_,
                "clusters": self.clusters_,
                "Q": self.Q_,
                "S": self.S_,
                "meta_keys": self.meta_keys_,
                "predictions": self._predictions,
                "include_cancer": self.include_cancer_,
                "cancer_keys": self.cancer_keys_,
            },
        }
        with open(path, "wb") as f:
            pickle.dump(state, f, protocol=pickle.HIGHEST_PROTOCOL)

    @classmethod
    def load(cls, path):
        """Reload a saved fit. Returns a CellConsensus ready for predict_*."""
        with open(path, "rb") as f:
            state = pickle.load(f)

        from . import __version__ as _cc_version
        saved_version = state.get("version", "unknown")
        if saved_version != _cc_version:
            warnings.warn(
                f"Loading CellConsensus model saved with version "
                f"{saved_version}; installed version is {_cc_version}."
            )

        obj = cls(**state["hyperparameters"])
        f = state["fitted"]
        obj.obs_names_      = f["obs_names"]
        obj.var_names_      = f["var_names"]
        obj.nn_indices_     = f["nn_indices"]
        obj.clusters_       = f.get("clusters")
        obj.Q_              = f["Q"]
        obj.S_              = f["S"]
        obj.meta_keys_      = f["meta_keys"]
        obj._predictions    = f["predictions"]
        obj.include_cancer_ = f.get("include_cancer", False)
        obj.cancer_keys_    = f.get("cancer_keys")
        return obj

    def _check_fitted(self):
        if self.S_ is None:
            raise RuntimeError("Call fit() before predict().")
