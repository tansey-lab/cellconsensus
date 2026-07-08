"""``cellconsensus-visualize`` — render annotation summary plots to PDF.

Reads an annotated AnnData ``.h5ad`` (typically the ``--output-h5ad`` from
``cellconsensus-annotate``) and writes a multi-page PDF with:

* bar plots of the cell-type breakdown at every level,
* a spatial scatter of ``obsm['spatial']`` coloured by cell type (if present),
* histograms of confidence scores and any other numeric score columns.
"""
import argparse
import sys

from .. import __version__


def build_parser():
    p = argparse.ArgumentParser(
        prog="cellconsensus-visualize",
        description="Render cell-type breakdown, spatial, and score plots "
                    "from an annotated AnnData to a PDF.",
    )
    p.add_argument("--version", action="version",
                   version=f"cellconsensus {__version__}")
    p.add_argument("input", help="Path to annotated AnnData (.h5ad) file.")
    p.add_argument(
        "-o", "--output", default=None,
        help="Output PDF path (default: <input>.plots.pdf).",
    )
    p.add_argument(
        "--level-cols", default=None,
        help="Comma-separated obs columns to treat as cell-type levels "
             "(default: auto-detect cellconsensus_level_* columns).",
    )
    p.add_argument(
        "--score-cols", default=None,
        help="Comma-separated numeric obs columns to histogram "
             "(default: auto-detect *_score and numeric cellconsensus_* "
             "columns).",
    )
    p.add_argument(
        "--spatial-key", default="spatial",
        help="Key in adata.obsm with 2-D coordinates (default: spatial).",
    )
    p.add_argument(
        "--point-size", type=float, default=6.0,
        help="Marker size for the spatial scatter (default: 6).",
    )
    p.add_argument("-q", "--quiet", action="store_true",
                   help="Suppress progress output.")
    return p


def _split_cols(value):
    if value is None:
        return None
    return [c.strip() for c in value.split(",") if c.strip()]


def main(argv=None):
    import anndata

    from ..plotting import visualize

    args = build_parser().parse_args(argv)
    verbose = not args.quiet

    output = args.output
    if output is None:
        output = args.input
        for ext in (".h5ad", ".h5", ".ad"):
            if output.endswith(ext):
                output = output[: -len(ext)]
                break
        output += ".plots.pdf"

    if verbose:
        print(f"Reading {args.input}...", file=sys.stderr)
    adata = anndata.read_h5ad(args.input)

    try:
        visualize(
            adata,
            output,
            level_cols=_split_cols(args.level_cols),
            score_cols=_split_cols(args.score_cols),
            spatial_key=args.spatial_key,
            point_size=args.point_size,
            verbose=verbose,
        )
    except ValueError as e:
        raise SystemExit(f"error: {e}")


if __name__ == "__main__":
    main()
