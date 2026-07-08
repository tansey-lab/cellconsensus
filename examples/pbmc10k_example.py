#!/usr/bin/env python3
"""CellConsensus on 10x's 10k PBMC (v3) — produces the two README figures.

Self-contained: downloads the public 10x filtered matrix (~37 MB) into
``data/`` on first run, then caches it there.

Outputs:
  assets/pbmc_cell_types.png   — UMAP colored by level 1 / 2 / 3 labels
  assets/pbmc_score_lvl1.png   — UMAP per-type level-1 score heatmaps
"""

import os
import urllib.request

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc

from cellconsensus import CellConsensus

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ASSETS = os.path.join(HERE, "assets")
DATA = os.path.join(HERE, "data")
os.makedirs(ASSETS, exist_ok=True)
os.makedirs(DATA, exist_ok=True)

URL = (
    "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/"
    "pbmc_10k_v3_filtered_feature_bc_matrix.h5"
)
H5 = os.path.join(DATA, "pbmc10k_v3.h5")

if not os.path.exists(H5):
    print(f"Downloading 10k PBMC (v3) from 10x Genomics...\n  {URL}")
    urllib.request.urlretrieve(URL, H5)

print("Loading 10k PBMC (v3)...")
adata = sc.read_10x_h5(H5)
adata.var_names_make_unique()
# Light QC: drop near-empty cells and rarely-detected genes.
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
print(f"  {adata.n_obs} cells × {adata.n_vars} genes")

cc = CellConsensus()
cc.fit(adata)

adata.obs["cc_lvl1"] = cc.predict(level=1).values
adata.obs["cc_lvl2"] = cc.predict(level=2).values
adata.obs["cc_lvl3"] = cc.predict(level=3).values
S = cc.score_matrix(level=1)

print("Computing UMAP...")
sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=15)
sc.tl.umap(adata)
emb = adata.obsm["X_umap"]


def scatter_cat(ax, labels, title):
    cats = sorted(np.unique(labels).tolist())
    cmap = plt.get_cmap("tab20")
    for i, c in enumerate(cats):
        m = labels == c
        ax.scatter(
            emb[m, 0],
            emb[m, 1],
            s=5,
            c=[cmap(i % 20)],
            label=f"{c} ({m.sum()})",
            linewidths=0,
            rasterized=True,
        )
    ax.set_title(title, fontsize=13)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.legend(
        loc="center left",
        bbox_to_anchor=(1.0, 0.5),
        fontsize=8,
        markerscale=2.5,
        frameon=False,
    )


# --- Figure 1: CellConsensus levels 1 / 2 / 3 ---
fig, axes = plt.subplots(1, 3, figsize=(22, 6.5))
for ax, lvl in zip(axes, (1, 2, 3)):
    scatter_cat(
        ax, adata.obs[f"cc_lvl{lvl}"].astype(str).values, f"CellConsensus level {lvl}"
    )
fig.suptitle(
    f"10k PBMC (v3) — CellConsensus predictions  (n={adata.n_obs})", fontsize=14, y=1.02
)
fig.tight_layout()
out1 = os.path.join(ASSETS, "pbmc_cell_types.png")
fig.savefig(out1, dpi=130, bbox_inches="tight")
plt.close(fig)
print(f"  wrote {out1}")

# --- Figure 2: level-1 per-type score heatmaps on UMAP ---
cols = list(S.columns)
vals_all = S.values
n = len(cols)
ncol = 4
nrow = int(np.ceil(n / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(4 * ncol, 3.6 * nrow))
axes = np.atleast_2d(axes).ravel()
vmax = float(np.quantile(vals_all, 0.995))
for j, name in enumerate(cols):
    v = vals_all[:, j]
    order = np.argsort(v)
    sc_h = axes[j].scatter(
        emb[order, 0],
        emb[order, 1],
        c=v[order],
        s=5,
        cmap="viridis",
        vmin=0,
        vmax=vmax,
        linewidths=0,
        rasterized=True,
    )
    axes[j].set_title(name, fontsize=11)
    axes[j].set_xticks([])
    axes[j].set_yticks([])
    plt.colorbar(sc_h, ax=axes[j], fraction=0.046, pad=0.02)
for k in range(n, len(axes)):
    axes[k].axis("off")
fig.suptitle(
    "10k PBMC (v3) — CellConsensus level-1 score per cell type", fontsize=14, y=1.0
)
fig.tight_layout()
out2 = os.path.join(ASSETS, "pbmc_score_lvl1.png")
fig.savefig(out2, dpi=130, bbox_inches="tight")
plt.close(fig)
print(f"  wrote {out2}")
