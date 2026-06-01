"""Mapping between meta keys, CL IDs, and human-readable names."""

# Meta key -> representative CL ID (one per meta level)
META_TO_CL = {
    # Level 1
    "b_plasma": "CL:0000236", "t_cell": "CL:0000084", "myeloid_cell": "CL:0000763",
    "mast_cell": "CL:0000097", "erythroid_megakaryocyte": "CL:0000764",
    "fibroblast": "CL:0000057", "smc_pericyte": "CL:0000192",
    "endothelial": "CL:0000115", "epithelial": "CL:0000066",
    "neural": "CL:0002319", "stem_progenitor": "CL:0000034",
    "germ_cell": "CL:0000586", "adipocyte": "CL:0000136",
    "parathyroid": "CL:0000446", "muscle": "CL:0000187",
    # Level 2
    "b_cell": "CL:0000236", "plasma_cell": "CL:0000786",
    "cd4_t": "CL:0000624", "cd8_t": "CL:0000625", "gd_t": "CL:0000798",
    "nk": "CL:0000623", "ilc": "CL:0001065",
    "monocyte": "CL:0000576", "macrophage": "CL:0000235",
    "dendritic_cell": "CL:0000451", "neutrophil": "CL:0000775",
    "eosinophil_basophil": "CL:0000771",
    "erythroid": "CL:0000764", "megakaryocyte": "CL:0000556",
    "vascular_endothelial": "CL:0000071",
    "lymphatic_endothelial": "CL:0002138",
    "smooth_muscle": "CL:0000192", "pericyte": "CL:0000669",
    "absorptive_epithelial": "CL:0000584", "secretory_epithelial": "CL:0000160",
    "basal_squamous": "CL:0000646", "endocrine_epithelial": "CL:0000163",
    "ciliated_epithelial": "CL:0000067",
    "neuron": "CL:0000540", "astrocyte": "CL:0000127",
    "oligodendrocyte": "CL:0000128", "schwann_cell": "CL:0002573",
    "hsc": "CL:0000037", "msc": "CL:0000134",
    # Level 3
    "naive_cd4": "CL:0000895", "memory_cd4": "CL:0000897",
    "treg": "CL:0000815", "helper_cd4": "CL:0000492",
    "naive_cd8": "CL:0000900", "effector_cd8": "CL:0000794",
    "memory_cd8": "CL:0000909",
    "classical_mono": "CL:0001054", "nonclassical_mono": "CL:0000875",
    "cdc1": "CL:0002394", "cdc2": "CL:0002399", "pdc": "CL:0000784",
    "naive_b": "CL:0000788", "memory_b": "CL:0000787",
    "germinal_center_b": "CL:0000844",
    "progenitor_erythroid": "CL:0000765", "mature_erythroid": "CL:0000232",
    "excitatory_neuron": "CL:0000679", "inhibitory_neuron": "CL:0000617",
    "hepatocyte": "CL:0000182", "enterocyte": "CL:0000584",
    "renal_tubular": "CL:0002306",
    "at2": "CL:0002063", "goblet": "CL:0000160", "club": "CL:0000158",
    "acinar": "CL:0000622", "ductal": "CL:0002079",
    "at1": "CL:0002062", "multiciliated": "CL:0005012",
    "beta_cell": "CL:0000169", "alpha_cell": "CL:0000171",
    "enteroendocrine": "CL:0000164",
}

# CL ID -> human-readable name (from cell_type table)
CL_TO_NAME = None


def _load_cl_names():
    """Load CL ID -> name mapping from cl_to_meta.json."""
    global CL_TO_NAME
    if CL_TO_NAME is not None:
        return CL_TO_NAME
    from .markers import load_cl_to_meta
    data = load_cl_to_meta()
    CL_TO_NAME = {cl: v["name"] for cl, v in data["cl_to_meta"].items()}
    # Also add from the full consensus cache ct_info
    from .markers import DATA_DIR
    import pickle
    full_path = DATA_DIR / "consensus_cache_full.pkl"
    if full_path.exists():
        with open(full_path, "rb") as f:
            cache = pickle.load(f)
        for cl, info in cache.get("ct_info", {}).items():
            if cl not in CL_TO_NAME:
                CL_TO_NAME[cl] = info.get("name", cl)
    return CL_TO_NAME


def meta_to_cl(meta_key):
    """Convert a meta key to its representative CL ID."""
    return META_TO_CL.get(meta_key, meta_key)


def cl_to_name(cl_id):
    """Convert a CL ID to its human-readable name."""
    names = _load_cl_names()
    return names.get(cl_id, cl_id)


def meta_to_name(meta_key):
    """Convert a meta key to its human-readable name."""
    from .markers import load_meta_groups
    for level in [3, 2, 1]:
        groups = load_meta_groups(level=level)
        if meta_key in groups:
            return groups[meta_key]
    return meta_key
