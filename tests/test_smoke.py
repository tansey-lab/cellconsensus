"""Smoke test: the package imports and exposes its public API."""

import cellconsensus


def test_import():
    assert cellconsensus.__version__
    assert cellconsensus.__all__
