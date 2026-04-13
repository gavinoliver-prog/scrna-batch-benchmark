"""
01_preprocess.py
Merge PBMC 3k + PBMC 68k, run standard Scanpy QC/normalization/HVG/PCA.
Output: data/processed/merged_preprocessed.h5ad
Cell type labels assigned via Leiden clustering + marker-based annotation.
"""

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = ROOT / "data" / "raw"
PROC_DIR = ROOT / "data" / "processed"
PROC_DIR.mkdir(parents=True, exist_ok=True)
DONE_MARKER = ROOT / "results" / "01_preprocess.done"

print("=== Step 01: Preprocessing ===")

# --- Load raw data ---
adata_3k  = sc.read_h5ad(RAW_DIR / "pbmc3k_raw.h5ad")
adata_68k = sc.read_h5ad(RAW_DIR / "pbmc68k_raw.h5ad")

# Align to shared genes
shared_genes = adata_3k.var_names.intersection(adata_68k.var_names)
print(f"Shared genes: {len(shared_genes)}")
adata_3k  = adata_3k[:, shared_genes].copy()
adata_68k = adata_68k[:, shared_genes].copy()

# Merge
adata = ad.concat([adata_3k, adata_68k], label="batch", keys=["pbmc3k", "pbmc68k"],
                  join="inner", merge="same")
adata.obs_names_make_unique()
print(f"Merged: {adata.shape[0]} cells x {adata.shape[1]} genes")

# --- QC ---
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None,
                            log1p=False, inplace=True)

# Diagnostics before filtering
print(f"  Mito % — median: {adata.obs.pct_counts_mt.median():.1f}%, "
      f"p95: {adata.obs.pct_counts_mt.quantile(0.95):.1f}%")

# 20% mito cap keeps both datasets — PBMC 3k is clean (<5%), GSE132044 v3
# cells can run slightly higher due to different library prep; doublets /
# dying cells at the extreme tail are still removed.
mito_thresh = 20
adata = adata[adata.obs.pct_counts_mt < mito_thresh, :]
adata = adata[adata.obs.n_genes_by_counts < 6000, :]
print(f"After QC: {adata.shape[0]} cells")

# Store raw counts BEFORE normalization (needed for seurat_v3 HVG)
adata.layers["counts"] = adata.X.copy()

# --- HVG on raw counts (seurat_v3 requires pre-normalization counts) ---
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="batch",
                             flavor="seurat_v3", layer="counts", subset=False)
print(f"HVGs: {adata.var.highly_variable.sum()}")

# --- Normalize & log-transform ---
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata  # freeze normalized log counts

# --- Scale & PCA on HVGs ---
adata_hvg = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(adata_hvg, max_value=10)
sc.tl.pca(adata_hvg, n_comps=50, svd_solver="arpack")

# Copy PCA back to main object
adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"].copy()
adata.varm["PCs"] = np.zeros((adata.n_vars, 50))  # placeholder for non-HVG vars
adata.uns["pca"] = adata_hvg.uns["pca"].copy()

# --- Uncorrected neighbors + UMAP (baseline) ---
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, use_rep="X_pca")
sc.tl.leiden(adata, resolution=0.5, key_added="leiden_uncorrected")
sc.tl.umap(adata)
adata.obsm["X_umap_uncorrected"] = adata.obsm["X_umap"].copy()

# --- Coarse cell type annotation (marker-based) ---
# PBMC marker genes -> broad cell type map
marker_genes = {
    "CD14+ Mono":  ["CD14", "LYZ"],
    "CD4 T":       ["IL7R", "CD4", "CCR7"],
    "CD8 T":       ["CD8A", "CD8B"],
    "NK":          ["GNLY", "NKG7"],
    "B cell":      ["MS4A1", "CD79A"],
    "DC":          ["FCER1A", "CST3"],
    "Platelet":    ["PPBP"],
}

sc.tl.rank_genes_groups(adata, groupby="leiden_uncorrected", method="wilcoxon",
                         key_added="rank_genes_leiden")

# Assign cell types by scoring
sc.tl.score_genes(adata, gene_list=[g for genes in marker_genes.values() for g in genes],
                  score_name="_marker_score")

# Simple assignment: score each cluster against each cell type
cell_type_map = {}
for cluster in adata.obs["leiden_uncorrected"].unique():
    mask = adata.obs["leiden_uncorrected"] == cluster
    best_type = "Unknown"
    best_score = -np.inf
    for ct, genes in marker_genes.items():
        available = [g for g in genes if g in adata.var_names]
        if not available:
            continue
        score = adata[mask, available].X.mean()
        if hasattr(score, "toarray"):
            score = score.toarray().mean()
        if score > best_score:
            best_score = score
            best_type = ct
    cell_type_map[cluster] = best_type

adata.obs["cell_type"] = adata.obs["leiden_uncorrected"].map(cell_type_map).astype("category")
print("Cell type distribution:")
print(adata.obs["cell_type"].value_counts())

# --- Save ---
adata.write_h5ad(PROC_DIR / "merged_preprocessed.h5ad")
print(f"Saved: {PROC_DIR / 'merged_preprocessed.h5ad'}")

DONE_MARKER.touch()
print("=== Step 01 complete ===\n")
