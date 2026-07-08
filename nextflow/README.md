# nf-core/cellconsensus

## Introduction

**nf-core/cellconsensus** is a bioinformatics pipeline that runs the [CellConsensus](https://github.com/tansey-lab/cellconsensus) hierarchical cell-type annotation on a set of AnnData (`.h5ad`) samples.

The pipeline currently exposes a single step:

1. `CELLCONSENSUS_ANNOTATE`: fits `CellConsensus` on each input `.h5ad`, emits per-cell hierarchical predictions as CSV and writes the annotated AnnData back to disk.

## Usage

Prepare a samplesheet (CSV) with the following columns:

```csv
sample_id,adata_path
SAMPLE_1,/path/to/data/sample_1.h5ad
SAMPLE_2,/path/to/data/sample_2.h5ad
```

Then run the pipeline:

```bash
nextflow run tansey-lab/cellconsensus \
    -profile <docker/singularity/.../institute> \
    --input samplesheet.csv \
    --outdir <OUTDIR>
```

### Common parameters

- `--clustering {ccc,precomputed}` (default: `ccc`)
- `--cluster-key <col>` (required for `precomputed`)
- `--levels 1,2,3`
- `--output-format {name,cl_id,key}`
- `--include-cancer` and `--cancer-types <keys>`

See `nextflow run tansey-lab/cellconsensus --help_full` for the full parameter list.

## Pipeline output

Per-sample results are published under `<outdir>/annotate/<sample_id>_annotate/`:

- `predictions.csv` — hierarchical predictions per cell
- `annotated.h5ad` — original AnnData with `cellconsensus_level_*` columns added
- `model.pkl` — fitted CellConsensus model (only if `--save_model`)

## Credits

nf-core/cellconsensus wraps the [CellConsensus](https://github.com/tansey-lab/cellconsensus) package.

## Citations

See the [`CITATIONS.md`](CITATIONS.md) file for references.
