"""``cellconsensus-annotate`` — annotate an .h5ad file with cell types.

Reads an AnnData ``.h5ad`` file, fits :class:`~cellconsensus.core.CellConsensus`,
and writes hierarchical predictions to CSV and/or back into the AnnData.
"""

import argparse
import sys

from .. import __version__


def build_parser():
    p = argparse.ArgumentParser(
        prog="cellconsensus-annotate",
        description="Fit CellConsensus on an AnnData file and write predictions.",
    )
    p.add_argument(
        "--version", action="version", version=f"cellconsensus {__version__}"
    )
    p.add_argument("input", help="Path to input AnnData (.h5ad) file.")
    p.add_argument(
        "-o",
        "--output",
        default=None,
        help="Path to write predictions CSV (default: stdout).",
    )
    p.add_argument(
        "--output-h5ad",
        default=None,
        help="If set, write the AnnData with predictions in .obs to this path.",
    )
    p.add_argument(
        "--clustering",
        choices=("ccc", "precomputed"),
        default="ccc",
        help="Clustering strategy (default: ccc).",
    )
    p.add_argument(
        "--cluster-key",
        default=None,
        help="Column in adata.obs with cluster labels (precomputed mode).",
    )
    p.add_argument(
        "--levels",
        default="1,2,3",
        help="Comma-separated levels to emit, from {1,2,3} (default: 1,2,3).",
    )
    p.add_argument(
        "--output-format",
        choices=("name", "cl_id", "key"),
        default="name",
        help="Label representation (default: name).",
    )
    p.add_argument(
        "--n-neighbors",
        type=int,
        default=20,
        help="Level-1 nearest neighbors (ccc mode, default: 20).",
    )
    p.add_argument(
        "--n-neighbors-lvl2",
        type=int,
        default=10,
        help="Level-2 nearest neighbors (ccc mode, default: 10).",
    )
    p.add_argument(
        "--n-neighbors-lvl3",
        type=int,
        default=5,
        help="Level-3 nearest neighbors (ccc mode, default: 5).",
    )
    p.add_argument(
        "--n-smooth",
        type=int,
        default=1,
        help="Smoothing iterations (ccc mode, default: 1).",
    )
    p.add_argument(
        "--ref-top-k",
        type=int,
        default=200,
        help="Top consensus markers per cell type (default: 200).",
    )
    p.add_argument(
        "--graph-level",
        type=int,
        default=3,
        choices=(1, 2, 3),
        help="Consensus level for the ccc kNN graph (default: 3).",
    )
    p.add_argument(
        "--include-cancer",
        action="store_true",
        help="Add cancer key(s) as level-1 classes during fit.",
    )
    p.add_argument(
        "--cancer-types",
        default=None,
        help="Comma-separated cancer keys for --include-cancer "
        "(default: pan-cancer 'cancer'). See cellconsensus-list-cancer-types.",
    )
    p.add_argument(
        "--save-model", default=None, help="Persist the fitted model to this .pkl path."
    )
    p.add_argument(
        "-q", "--quiet", action="store_true", help="Suppress progress output."
    )
    return p


def _parse_levels(levels_str):
    try:
        levels = [int(x) for x in levels_str.split(",") if x.strip()]
    except ValueError:
        raise SystemExit(f"error: invalid --levels value: {levels_str!r}")
    if not levels or any(lvl not in (1, 2, 3) for lvl in levels):
        raise SystemExit("error: --levels must be a comma-separated subset of {1,2,3}.")
    return levels


def main(argv=None):
    import anndata
    import pandas as pd

    from ..core import CellConsensus

    args = build_parser().parse_args(argv)
    levels = _parse_levels(args.levels)
    verbose = not args.quiet

    if verbose:
        print(f"Reading {args.input}...", file=sys.stderr)
    adata = anndata.read_h5ad(args.input)

    cancer_types = None
    if args.cancer_types:
        cancer_types = [c.strip() for c in args.cancer_types.split(",") if c.strip()]

    cc = CellConsensus(
        clustering=args.clustering,
        n_neighbors=args.n_neighbors,
        n_neighbors_lvl2=args.n_neighbors_lvl2,
        n_neighbors_lvl3=args.n_neighbors_lvl3,
        n_smooth=args.n_smooth,
        ref_top_k=args.ref_top_k,
        graph_level=args.graph_level,
        cluster_key=args.cluster_key,
    )
    cc.fit(
        adata,
        include_cancer=args.include_cancer,
        cancer_types=cancer_types,
        verbose=verbose,
    )

    columns = {}
    for lvl in levels:
        columns[f"level_{lvl}"] = cc.predict(level=lvl, output=args.output_format)
        columns[f"level_{lvl}_score"] = cc.confidence(level=lvl).values
    result = pd.DataFrame(columns, index=adata.obs_names)

    if args.save_model:
        cc.save(args.save_model)
        if verbose:
            print(f"Saved model to {args.save_model}", file=sys.stderr)

    if args.output_h5ad:
        for col in result.columns:
            adata.obs[f"cellconsensus_{col}"] = result[col].values
        adata.write_h5ad(args.output_h5ad)
        if verbose:
            print(f"Wrote annotated AnnData to {args.output_h5ad}", file=sys.stderr)

    if args.output:
        result.to_csv(args.output)
        if verbose:
            print(f"Wrote predictions to {args.output}", file=sys.stderr)
    elif not args.output_h5ad:
        result.to_csv(sys.stdout)


if __name__ == "__main__":
    main()
