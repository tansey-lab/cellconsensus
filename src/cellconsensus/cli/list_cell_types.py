"""``cellconsensus-list-cell-types`` — list valid cell-type keys by level."""

import argparse

from .. import __version__


def build_parser():
    p = argparse.ArgumentParser(
        prog="cellconsensus-list-cell-types",
        description="List valid cell-type keys at a given consensus level, "
        "one 'key<TAB>name' pair per line.",
    )
    p.add_argument(
        "--version", action="version", version=f"cellconsensus {__version__}"
    )
    p.add_argument(
        "--level",
        type=int,
        default=1,
        choices=(1, 2, 3),
        help="Consensus level (default: 1).",
    )
    return p


def main(argv=None):
    from ..consensus import load_cell_type

    args = build_parser().parse_args(argv)
    for key, name in load_cell_type(level=args.level).items():
        print(f"{key}\t{name}")


if __name__ == "__main__":
    main()
