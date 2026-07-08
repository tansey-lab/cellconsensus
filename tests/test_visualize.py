"""Tests for the cellconsensus-visualize CLI and plotting helpers."""

import anndata as ad
import numpy as np
import pandas as pd
import pytest


def _annotated_adata(n=200, spatial=True):
    rng = np.random.RandomState(0)
    obs = pd.DataFrame(index=[f"c{i}" for i in range(n)])
    obs["cellconsensus_level_1"] = rng.choice(["t_cell", "myeloid"], n)
    obs["cellconsensus_level_2"] = rng.choice(["cd4_t", "cd8_t", "nk"], n)
    obs["cellconsensus_level_1_score"] = rng.rand(n)
    obs["cellconsensus_level_2_score"] = rng.rand(n)
    a = ad.AnnData(X=rng.poisson(1.0, size=(n, 15)).astype(float), obs=obs)
    if spatial:
        a.obsm["spatial"] = rng.rand(n, 2) * 100
    return a


def test_visualize_writes_pdf(tmp_path):
    from cellconsensus.plotting import visualize

    adata = _annotated_adata()
    out = tmp_path / "plots.pdf"
    rendered = visualize(adata, str(out), verbose=False)

    assert out.exists() and out.stat().st_size > 0
    assert rendered["barplots"] == ["cellconsensus_level_1",
                                    "cellconsensus_level_2"]
    assert rendered["spatial"] == ["cellconsensus_level_1",
                                   "cellconsensus_level_2"]
    assert set(rendered["histograms"]) == {"cellconsensus_level_1_score",
                                           "cellconsensus_level_2_score"}


def test_visualize_without_spatial(tmp_path):
    from cellconsensus.plotting import visualize

    adata = _annotated_adata(spatial=False)
    out = tmp_path / "plots.pdf"
    rendered = visualize(adata, str(out), verbose=False)

    assert out.exists()
    assert rendered["spatial"] == []
    assert rendered["barplots"]


def test_visualize_nothing_to_plot(tmp_path):
    from cellconsensus.plotting import visualize

    adata = ad.AnnData(X=np.zeros((10, 3)))
    with pytest.raises(ValueError, match="Nothing to plot"):
        visualize(adata, str(tmp_path / "x.pdf"), verbose=False)


def test_visualize_cli(tmp_path):
    from cellconsensus.cli.visualize import main

    h5 = tmp_path / "ann.h5ad"
    _annotated_adata().write_h5ad(str(h5))
    main([str(h5), "-o", str(tmp_path / "out.pdf"), "-q"])
    assert (tmp_path / "out.pdf").exists()

    # Default output path derives from the input name.
    main([str(h5), "-q"])
    assert (tmp_path / "ann.plots.pdf").exists()
