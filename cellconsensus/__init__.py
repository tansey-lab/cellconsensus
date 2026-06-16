"""CellConsensus: hierarchical unsupervised cell type annotation."""

__version__ = "1.2.0"

from .core import CellConsensus
from .cancer import list_cancer_types
from .consensus import load_cell_type

__all__ = ["CellConsensus", "list_cancer_types", "load_cell_type"]
