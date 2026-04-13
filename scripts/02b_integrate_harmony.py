"""
02b_integrate_harmony.py
Batch correction using Harmony (via harmonypy).
Corrects PCA embedding, recomputes neighbors + UMAP on corrected embedding.
Output: results/corrected/harmony.h5ad
"""

import scanpy as sc
import harmonypy as hm
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
INPUT  = ROOT / "data" / "processed" / "merged_preprocessed.h5ad"
OUTPUT = ROOT / "results" / "corrected" / "harmony.h5ad"
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
DONE_MARKER = ROOT / "results" / "02b_harmony.done"

print("=== Step 02b: Harmony Integration ===")

adata = sc.read_h5ad(INPUT)
print(f"Loaded: {adata.shape[0]} cells")

# Run Harmony on PCA embedding
ho = hm.run_harmony(
    adata.obsm["X_pca"][:, :30],
    adata.obs,
    vars_use=["batch"],
    max_iter_harmony=20,
    random_state=42,
    verbose=True,
)

adata.obsm["X_pca_harmony"] = ho.Z_corr  # cells x components (harmonypy 0.2+ convention)

# Recompute neighbors on harmony-corrected PCA
sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_neighbors=15, n_pcs=30,
                key_added="neighbors_harmony")

# UMAP
sc.tl.umap(adata, neighbors_key="neighbors_harmony")
sc.tl.leiden(adata, resolution=0.5, key_added="leiden_harmony",
             neighbors_key="neighbors_harmony")

# Output contract
assert "X_umap" in adata.obsm
assert "batch" in adata.obs
assert "cell_type" in adata.obs

adata.write_h5ad(OUTPUT)
print(f"Saved: {OUTPUT}")
DONE_MARKER.touch()
print("=== Step 02b complete ===\n")
