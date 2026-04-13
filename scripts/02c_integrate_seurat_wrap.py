"""
02c_integrate_seurat_wrap.py
Python wrapper: calls 02c_integrate_seurat.R via subprocess,
then reads the output PCA embedding and builds a standardized AnnData.
Fails gracefully if R is not available.
Output: results/corrected/seurat.h5ad
"""

import subprocess
import sys
import shutil
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW_DIR  = ROOT / "data" / "raw"
PROC_DIR = ROOT / "data" / "processed"
OUTPUT   = ROOT / "results" / "corrected" / "seurat.h5ad"
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
DONE_MARKER = ROOT / "results" / "02c_seurat.done"
SKIP_MARKER = ROOT / "results" / "02c_seurat.skipped"

RSCRIPT  = ROOT / "scripts" / "02c_integrate_seurat.R"
TMP_DIR  = ROOT / "results" / "corrected" / "_seurat_tmp"
TMP_DIR.mkdir(parents=True, exist_ok=True)

PCA_CSV  = TMP_DIR / "seurat_pca.csv"
META_CSV = TMP_DIR / "seurat_meta.csv"

print("=== Step 02c: Seurat CCA Integration (via R) ===")

# Check R availability
if shutil.which("Rscript") is None:
    print("  WARNING: Rscript not found in PATH. Skipping Seurat arm.")
    SKIP_MARKER.touch()
    sys.exit(0)

# Call R script
cmd = [
    "Rscript", str(RSCRIPT),
    str(RAW_DIR / "pbmc3k_raw.h5ad"),
    str(RAW_DIR / "pbmc68k_raw.h5ad"),
    str(PCA_CSV),
    str(META_CSV),
]

print("  Running R script...")
result = subprocess.run(cmd, capture_output=True, text=True)

if result.returncode != 0:
    print("  ERROR: R script failed.")
    print(result.stderr[-2000:])
    print("  Skipping Seurat arm.")
    SKIP_MARKER.touch()
    sys.exit(0)

print(result.stdout)

# Read outputs
if not PCA_CSV.exists() or not META_CSV.exists():
    print("  ERROR: R script did not produce expected output files. Skipping.")
    SKIP_MARKER.touch()
    sys.exit(0)

pca_df  = pd.read_csv(PCA_CSV, index_col=0)
meta_df = pd.read_csv(META_CSV, index_col=0)

# Load preprocessed AnnData to get expression + obs alignment
adata_pre = sc.read_h5ad(PROC_DIR / "merged_preprocessed.h5ad")

# Align cell barcodes (R may reorder)
common_cells = pca_df.index.intersection(adata_pre.obs_names)
if len(common_cells) < 100:
    print(f"  WARNING: Only {len(common_cells)} matching barcodes. Check barcode format.")

adata_seurat = adata_pre[common_cells].copy()
adata_seurat.obsm["X_pca_seurat"] = pca_df.loc[common_cells].values

# Recompute neighbors + UMAP on Seurat-integrated PCA
sc.pp.neighbors(adata_seurat, use_rep="X_pca_seurat", n_neighbors=15, n_pcs=30)
sc.tl.umap(adata_seurat)
sc.tl.leiden(adata_seurat, resolution=0.5, key_added="leiden_seurat")

# Output contract
assert "X_umap" in adata_seurat.obsm
assert "batch" in adata_seurat.obs
assert "cell_type" in adata_seurat.obs

adata_seurat.write_h5ad(OUTPUT)
print(f"Saved: {OUTPUT}")
DONE_MARKER.touch()
print("=== Step 02c complete ===\n")
