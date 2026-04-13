# scRNA-seq Batch Correction Benchmarking Pipeline
## Project Context

Exploratory benchmarking framework comparing standard scRNA-seq batch correction/integration methods.
NOT a production clinical pipeline. Designed for modularity, reproducibility, and portfolio demonstration.

## Data

- **Batch 1**: PBMC 3k (`sc.datasets.pbmc3k()`) — 10x Chromium v1, ~2,700 cells
- **Batch 2**: PBMC 68k downsampled to ~3k cells (`sc.datasets.pbmc68k_reduced()`) — different library prep
- **Batch column**: `batch` (values: `"pbmc3k"`, `"pbmc68k"`)
- **Cell type column**: `cell_type` (Louvain-derived, populated during preprocessing)
- **Organism**: Human
- **Input format**: loaded directly via Scanpy API, saved as `.h5ad`

## Project Structure

```
scrna-batch-benchmark/
├── CLAUDE.md                  ← this file
├── environment.yml            ← conda environment spec
├── data/
│   ├── raw/                   ← raw AnnData objects per batch
│   └── processed/             ← merged, preprocessed AnnData
├── results/
│   ├── corrected/             ← one .h5ad per method
│   ├── metrics/               ← scIB scores as CSV
│   └── figures/               ← UMAP PNGs + heatmap
├── scripts/
│   ├── 00_download_data.py
│   ├── 01_preprocess.py
│   ├── 02a_integrate_bbknn.py
│   ├── 02b_integrate_harmony.py
│   ├── 02c_integrate_seurat.R       ← R script (called via subprocess)
│   ├── 02c_integrate_seurat_wrap.py ← Python wrapper for R script
│   ├── 03_benchmark.py
│   └── 04_visualize.py
├── run_pipeline.sh            ← executes all steps in order
└── report/
    └── summary.html           ← auto-generated HTML report
```

## Execution Order

Always run scripts in numeric order:
```bash
bash run_pipeline.sh
```
Or step by step:
```bash
python scripts/00_download_data.py
python scripts/01_preprocess.py
python scripts/02a_integrate_bbknn.py
python scripts/02b_integrate_harmony.py
python scripts/02c_integrate_seurat_wrap.py   # requires R environment
python scripts/03_benchmark.py
python scripts/04_visualize.py
```

## Environment

- **Conda env name**: `scrna-bench`
- **Python**: 3.10
- **Key Python packages**: scanpy, anndata, bbknn, harmonypy, scib-metrics, matplotlib, seaborn, pandas, numpy, leidenalg, python-igraph
- **R packages**: Seurat (v4+), SeuratDisk, Matrix
- **R/Python bridge**: rpy2, anndata2ri (install separately after conda env)

Activate before running anything:
```bash
conda activate scrna-bench
```

## Key Design Decisions

1. **Each integration method produces a standardized AnnData** saved to `results/corrected/<method>.h5ad`
   - Must contain: `.obsm["X_umap"]`, `.obs["batch"]`, `.obs["cell_type"]`
   - Optional: `.obsm["X_emb_corrected"]` for embedding-level methods

2. **Preprocessing is run once** and saved to `data/processed/merged_preprocessed.h5ad`
   - All integration methods load from this single file
   - Never re-run preprocessing per method

3. **scIB metrics to compute** (Phase 1 scope):
   - `iLISI` — batch mixing
   - `NMI` — biological conservation (vs cell_type labels)
   - `ARI` — clustering agreement
   - Full scIB battery deferred to Phase 2

4. **Seurat arm**: R script handles the integration, writes a temporary `.h5ad`, Python wrapper picks it up.
   - If R/rpy2 environment fails, skip Seurat arm and log a warning — do not block pipeline.

5. **UMAP** computed per method independently (not shared layout)

## Integration Methods

| Arm | Method | Script | Output key |
|-----|--------|--------|------------|
| A | BBKNN | 02a | `X_umap` via corrected graph |
| B | Harmony | 02b | `X_umap` via harmony embedding |
| C | Seurat CCA | 02c | `X_umap` via integrated PCA |

## Output Contracts

`03_benchmark.py` expects:
- `results/corrected/*.h5ad` — one file per method
- Each file has `.obs["batch"]` and `.obs["cell_type"]`

`04_visualize.py` expects:
- `results/metrics/scores.csv` — columns: `method, iLISI, NMI, ARI`
- `results/corrected/*.h5ad` — for UMAP plots

## Error Handling Notes

- If `pbmc68k_reduced()` is deprecated in current Scanpy version, use `sc.datasets.pbmc68k_neighbors()` or download from https://figshare.com/articles/dataset/pbmc68k/5765138 and subsample to 3000 cells
- If scib-metrics import fails, try `import scib` (package name varies by version)
- All scripts should write a completion marker: `results/<step>.done`

## What NOT to Do

- Do not hardcode absolute paths — use `pathlib.Path` relative to project root
- Do not re-merge or re-preprocess inside integration scripts
- Do not install packages inside scripts — if something is missing, print the install command and exit cleanly
- Do not run Seurat arm if R is not available — fail gracefully
