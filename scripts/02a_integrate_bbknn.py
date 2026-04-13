"""
02a_integrate_bbknn.py
Batch correction using BBKNN (Batch Balanced KNN).
Loads preprocessed AnnData, applies BBKNN on PCA embedding,
computes UMAP on corrected neighbor graph.
Output: results/corrected/bbknn.h5ad
"""

import scanpy as sc
import bbknn
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
INPUT  = ROOT / "data" / "processed" / "merged_preprocessed.h5ad"
OUTPUT = ROOT / "results" / "corrected" / "bbknn.h5ad"
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
DONE_MARKER = ROOT / "results" / "02a_bbknn.done"

print("=== Step 02a: BBKNN Integration ===")

adata = sc.read_h5ad(INPUT)
print(f"Loaded: {adata.shape[0]} cells")

# BBKNN operates on X_pca, corrects the neighbor graph
bbknn.bbknn(adata, batch_key="batch", use_rep="X_pca", n_pcs=30,
            neighbors_within_batch=3)

# UMAP on corrected graph
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5, key_added="leiden_bbknn")

# Output contract: must have X_umap, batch, cell_type
assert "X_umap" in adata.obsm
assert "batch" in adata.obs
assert "cell_type" in adata.obs

adata.write_h5ad(OUTPUT)
print(f"Saved: {OUTPUT}")
DONE_MARKER.touch()
print("=== Step 02a complete ===\n")
