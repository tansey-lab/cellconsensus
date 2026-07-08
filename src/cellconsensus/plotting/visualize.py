"""Plotting helpers for :mod:`cellconsensus`.

Builds a multi-page PDF summarising a CellConsensus annotation:

* bar plots of the cell-type breakdown at every level,
* spatial scatter ("choropleth") of ``obsm['spatial']`` coloured by cell type,
* histograms of confidence scores and any other numeric score columns.

Matplotlib is imported lazily so importing :mod:`cellconsensus` stays cheap and
the base package keeps no hard plotting dependency.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def _discrete_colors(categories):
    """Return an ordered ``{category: rgba}`` map that scales past 20 classes."""
    import matplotlib.pyplot as plt

    cats = list(categories)
    n = len(cats)
    if n <= 10:
        cmap = plt.get_cmap("tab10")
        colors = [cmap(i) for i in range(n)]
    elif n <= 20:
        cmap = plt.get_cmap("tab20")
        colors = [cmap(i) for i in range(n)]
    else:
        cmap = plt.get_cmap("hsv")
        colors = [cmap(i / n) for i in range(n)]
    return dict(zip(cats, colors))


def _category_counts(series):
    """Ordered value counts (descending) for a label column, as a Series."""
    s = pd.Series(series).astype("object").fillna("NA")
    return s.value_counts(dropna=False)


def plot_celltype_barplot(ax, labels, title=None):
    """Horizontal bar plot of a cell-type breakdown onto ``ax``."""
    counts = _category_counts(labels)
    colors = _discrete_colors(counts.index)
    y = np.arange(len(counts))[::-1]  # largest at top
    ax.barh(y, counts.values, color=[colors[c] for c in counts.index], edgecolor="none")
    ax.set_yticks(y)
    ax.set_yticklabels(counts.index, fontsize=8)
    ax.set_xlabel("cells")
    total = int(counts.sum())
    for yi, v in zip(y, counts.values):
        ax.text(
            v,
            yi,
            f" {int(v)} ({100 * v / total:.1f}%)",
            va="center",
            ha="left",
            fontsize=7,
        )
    ax.margins(x=0.15)
    if title:
        ax.set_title(title, fontsize=11)
    return ax


def plot_spatial(ax, coords, labels, title=None, point_size=6, legend=True):
    """Spatial scatter of ``coords`` (n x 2) coloured by ``labels``."""
    coords = np.asarray(coords)
    labels = pd.Series(labels).astype("object").fillna("NA").values
    # Colour by frequency so the legend order matches the bar plots.
    order = list(_category_counts(labels).index)
    colors = _discrete_colors(order)
    for cat in order:
        mask = labels == cat
        ax.scatter(
            coords[mask, 0],
            coords[mask, 1],
            s=point_size,
            c=[colors[cat]],
            label=str(cat),
            linewidths=0,
            rasterized=True,
        )
    ax.set_aspect("equal")
    ax.invert_yaxis()  # image convention: y grows downward
    ax.set_xlabel("spatial_1")
    ax.set_ylabel("spatial_2")
    ax.set_xticks([])
    ax.set_yticks([])
    if title:
        ax.set_title(title, fontsize=11)
    if legend:
        ax.legend(
            loc="center left",
            bbox_to_anchor=(1.01, 0.5),
            fontsize=6,
            markerscale=2,
            frameon=False,
            ncol=1 if len(order) <= 30 else 2,
        )
    return ax


def plot_score_histogram(ax, values, title=None, bins=50):
    """Histogram of a numeric score column onto ``ax``."""
    v = np.asarray(pd.to_numeric(pd.Series(values), errors="coerce"), dtype=float)
    v = v[np.isfinite(v)]
    if v.size:
        ax.hist(v, bins=bins, color="#4C72B0", edgecolor="white", linewidth=0.3)
        med = float(np.median(v))
        ax.axvline(med, color="#C44E52", linestyle="--", linewidth=1)
        ax.text(
            0.98,
            0.95,
            f"n={v.size}\nmedian={med:.3g}",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=7,
        )
    else:
        ax.text(
            0.5,
            0.5,
            "no finite values",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=8,
        )
    ax.set_ylabel("cells")
    if title:
        ax.set_title(title, fontsize=10)
    return ax


def visualize(
    adata,
    output_path,
    level_cols=None,
    score_cols=None,
    spatial_key="spatial",
    point_size=6,
    verbose=True,
):
    """Write a multi-page PDF summary of a CellConsensus annotation.

    Parameters
    ----------
    adata : AnnData
        Annotated object. Cell-type labels are read from ``adata.obs``.
    output_path : str
        Destination ``.pdf`` path.
    level_cols : list of str or None
        Categorical ``obs`` columns to treat as cell-type levels. If None,
        auto-detect columns starting with ``cellconsensus_level_`` that do not
        end in ``_score``.
    score_cols : list of str or None
        Numeric ``obs`` columns to histogram. If None, auto-detect columns
        ending in ``_score`` plus any other numeric ``cellconsensus_`` columns.
    spatial_key : str
        Key in ``adata.obsm`` holding 2-D coordinates. Skipped if absent.
    point_size : float
        Marker size for the spatial scatter.
    verbose : bool

    Returns
    -------
    dict with the pages actually rendered (for logging/testing).
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    obs = adata.obs
    if level_cols is None:
        level_cols = [
            c
            for c in obs.columns
            if c.startswith("cellconsensus_level_") and not c.endswith("_score")
        ]
    else:
        level_cols = [c for c in level_cols if c in obs.columns]
    level_cols = _natural_sort(level_cols)

    if score_cols is None:
        score_cols = [
            c
            for c in obs.columns
            if c.endswith("_score")
            or (
                c.startswith("cellconsensus_") and pd.api.types.is_numeric_dtype(obs[c])
            )
        ]
        # keep only genuinely numeric columns, preserve order, drop dups
        seen = set()
        score_cols = [
            c
            for c in score_cols
            if pd.api.types.is_numeric_dtype(obs[c]) and not (c in seen or seen.add(c))
        ]
    else:
        score_cols = [c for c in score_cols if c in obs.columns]

    has_spatial = (
        spatial_key in adata.obsm and np.asarray(adata.obsm[spatial_key]).shape[1] >= 2
    )

    if not level_cols and not score_cols and not has_spatial:
        raise ValueError(
            "Nothing to plot: found no cell-type level columns "
            "(cellconsensus_level_*), no score columns (*_score), and no "
            f"obsm['{spatial_key}']. Pass --level-cols / --score-cols / "
            "--spatial-key explicitly, or run cellconsensus-annotate first."
        )

    rendered = {"barplots": [], "spatial": [], "histograms": []}

    with PdfPages(output_path) as pdf:
        # --- Page 1: cell-type breakdown bar plots (one row per level) ------
        if level_cols:
            fig, axes = plt.subplots(
                len(level_cols),
                1,
                figsize=(8.5, max(2.6, 2.6 * len(level_cols))),
                squeeze=False,
            )
            for ax, col in zip(axes[:, 0], level_cols):
                plot_celltype_barplot(ax, obs[col], title=col)
                rendered["barplots"].append(col)
            fig.suptitle("Cell-type breakdown", fontsize=13)
            fig.tight_layout(rect=(0, 0, 1, 0.98))
            pdf.savefig(fig)
            plt.close(fig)

        # --- Spatial choropleth: one page per level -------------------------
        if has_spatial:
            coords = np.asarray(adata.obsm[spatial_key])[:, :2]
            cols_for_spatial = level_cols or [None]
            for col in cols_for_spatial:
                fig, ax = plt.subplots(figsize=(9, 7.5))
                labels = (
                    obs[col] if col is not None else np.full(adata.n_obs, "all cells")
                )
                plot_spatial(
                    ax,
                    coords,
                    labels,
                    point_size=point_size,
                    title=f"Spatial: {col}" if col else "Spatial",
                )
                fig.tight_layout()
                pdf.savefig(fig)
                plt.close(fig)
                if col is not None:
                    rendered["spatial"].append(col)
        elif verbose:
            print(f"  (no obsm['{spatial_key}'] — skipping spatial plots)")

        # --- Score histograms -----------------------------------------------
        if score_cols:
            ncol = 2 if len(score_cols) > 1 else 1
            nrow = int(np.ceil(len(score_cols) / ncol))
            fig, axes = plt.subplots(
                nrow, ncol, figsize=(4.5 * ncol, 3.0 * nrow), squeeze=False
            )
            flat = axes.ravel()
            for ax, col in zip(flat, score_cols):
                plot_score_histogram(ax, obs[col], title=col)
                rendered["histograms"].append(col)
            for ax in flat[len(score_cols) :]:
                ax.axis("off")
            fig.suptitle("Score distributions", fontsize=13)
            fig.tight_layout(rect=(0, 0, 1, 0.97))
            pdf.savefig(fig)
            plt.close(fig)

    if verbose:
        print(
            f"Wrote {output_path}: "
            f"{len(rendered['barplots'])} barplot(s), "
            f"{len(rendered['spatial'])} spatial page(s), "
            f"{len(rendered['histograms'])} histogram(s)."
        )
    return rendered


def _natural_sort(cols):
    """Sort so ``level_2`` precedes ``level_10`` (digits compared as ints)."""
    import re

    def key(c):
        return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", c)]

    return sorted(cols, key=key)
