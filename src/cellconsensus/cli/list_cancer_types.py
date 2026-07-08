"""``cellconsensus-list-cancer-types`` — list valid cancer-type keys."""

import argparse

from .. import __version__


def build_parser():
    p = argparse.ArgumentParser(
        prog="cellconsensus-list-cancer-types",
        description="List valid cancer-type keys, one per line.",
    )
    p.add_argument(
        "--version", action="version", version=f"cellconsensus {__version__}"
    )
    return p


def main(argv=None):
    from ..cancer import list_cancer_types

    build_parser().parse_args(argv)
    for key in list_cancer_types():
        print(key)


if __name__ == "__main__":
    main()
