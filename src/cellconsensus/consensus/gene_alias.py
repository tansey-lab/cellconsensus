"""Gene alias resolution to current HUGO symbols."""

import json
from pathlib import Path

_gene_alias = None
_ALIAS_PATH = Path(__file__).parent.parent / "data" / "gene_alias.json"


def load_gene_alias():
    """Load gene alias mapping (old symbol -> current symbol)."""
    global _gene_alias
    if _gene_alias is not None:
        return _gene_alias

    if _ALIAS_PATH.exists():
        with open(_ALIAS_PATH) as f:
            _gene_alias = json.load(f)
    else:
        _gene_alias = {}

    return _gene_alias


def resolve_gene_names(var_names):
    """Resolve a list of gene names to current HUGO symbols."""
    alias = load_gene_alias()
    return [alias.get(g, g) for g in var_names]
