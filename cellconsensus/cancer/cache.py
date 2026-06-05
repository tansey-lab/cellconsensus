"""Cancer consensus cache loader and helpers.

The cancer cache (consensus_cache_cancer.pkl) contains 119 cancer types
keyed by slugified NCIT term names (e.g. 'melanoma', 'breast_carcinoma').
The root key 'cancer' (NCIT:C3262 Neoplasm) aggregates all descendants.
"""
import json
import pickle
from pathlib import Path

DATA_DIR = Path(__file__).parent.parent / "data"
CANCER_CACHE_PATH = DATA_DIR / "consensus_cache_cancer.pkl"
CANCER_CODES_PATH = DATA_DIR / "malignant_cancer_codes.json"

_cache = None
_codes = None
_key_to_ncit = None
_key_to_display_name = None


def load_cancer_cache():
    """Load the cancer consensus cache (cached after first call)."""
    global _cache
    if _cache is None:
        with open(CANCER_CACHE_PATH, "rb") as f:
            _cache = pickle.load(f)
    return _cache


def load_cancer_codes():
    """Load the NCIT cancer codes JSON (hierarchy + metadata)."""
    global _codes
    if _codes is None:
        with open(CANCER_CODES_PATH) as f:
            _codes = json.load(f)
    return _codes


def list_cancer_types():
    """Return the list of all valid cancer type keys, sorted alphabetically.

    Use this to see which strings can be passed as cancer keys to
    `CellConsensus.predict_score()`.

    Returns
    -------
    list of str
        Alphabetically sorted list of 119 valid cancer type keys.

    Examples
    --------
    >>> from cellconsensus import list_cancer_types
    >>> types = list_cancer_types()
    >>> print(len(types), types[:5])
    119 ['acute_erythroid_leukemia', 'acute_lymphoblastic_leukemia', ...]
    """
    cache = load_cancer_cache()
    return sorted(cache["consensus"].keys())


def _build_key_maps():
    """Build slug -> NCIT code and slug -> display-name maps from the codes JSON."""
    global _key_to_ncit, _key_to_display_name
    if _key_to_ncit is not None:
        return
    codes = load_cancer_codes()
    entries = codes.get("entries", {})
    _key_to_ncit = {info["key"]: ncit for ncit, info in entries.items()}
    _key_to_display_name = {
        info["key"]: info["name"] for ncit, info in entries.items()
    }


def cancer_key_to_ncit(key):
    """Return the NCIT code for a cancer slug, or the slug itself if missing."""
    _build_key_maps()
    return _key_to_ncit.get(key, key)


def cancer_key_to_name(key):
    """Return the display name for a cancer slug (e.g. "Breast Carcinoma").

    Falls back to a title-cased version of the slug if the key is unknown.
    """
    _build_key_maps()
    if key in _key_to_display_name:
        return _key_to_display_name[key]
    return key.replace("_", " ").title()


def is_cancer_key(key):
    """True if `key` is one of the keys in the cancer consensus cache."""
    cache = load_cancer_cache()
    return key in cache["consensus"]


def validate_cancer_types(cancer_types):
    """Check cancer_types against the cache. Raise ValueError if unknown.

    Parameters
    ----------
    cancer_types : list of str or None
        None is always valid (defaults to pan-cancer "cancer").

    Raises
    ------
    ValueError
        With a helpful message pointing to list_cancer_types().
    """
    if cancer_types is None:
        return
    cache = load_cancer_cache()
    valid = set(cache["consensus"].keys())
    unknown = [c for c in cancer_types if c not in valid]
    if unknown:
        raise ValueError(
            f"Unknown cancer type(s): {unknown}. "
            f"Import the full list of valid names with:\n"
            f"  from cellconsensus import list_cancer_types\n"
            f"  list_cancer_types()  # returns {len(valid)} alphabetically-sorted keys"
        )
