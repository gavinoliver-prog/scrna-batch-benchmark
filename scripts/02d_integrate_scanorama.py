"""
02d_integrate_scanorama.py
Scanorama manifold alignment integration.
Loads merged preprocessed AnnData, runs Scanorama on PCA embeddings per batch,
computes UMAP + Leiden, saves standardised AnnData.
Output: results/corrected/scanorama.h5ad
"""

import sys
import scanpy as sc
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PROC_DIR = ROOT / "data" / "processed"
OUTPUT = ROOT / "results" / "corrected" / "scanorama.h5ad"
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
DONE_MARKER = ROOT / "results" / "02d_scanorama.done"

print("=== Step 02d: Scanorama Integration ===")

try:
    import scanorama  # noqa: F401
except ImportError:
    print("ERROR: scanorama not found.")
    print("Install with: pip install scanorama")
    sys.exit(1)

adata = sc.read_h5ad(PROC_DIR / "merged_preprocessed.h5ad")
print(f"Loaded: {adata.shape[0]} cells × {adata.shape[1]} genes")

# Scanorama integration on per-batch PCA embeddings.
# sc.external.pp.scanorama_integrate takes obsm[basis] sliced by batch key,
# aligns them, and stores corrected low-dim embedding in obsm[adjusted_basis].
print("Running Scanorama (basis=X_pca)...")
sc.external.pp.scanorama_integrate(
    adata,
    key="batch",
    basis="X_pca",
    adjusted_basis="X_scanorama",
    knn=20,
    sigma=15,
    approx=True,
    alpha=0.1,
    batch_size=5000,
)

print(f"  Scanorama embedding shape: {adata.obsm['X_scanorama'].shape}")

# Neighbors + UMAP + Leiden on Scanorama-corrected embedding
sc.pp.neighbors(adata, use_rep="X_scanorama", n_neighbors=15)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5, key_added="leiden_scanorama")

# Output contract
assert "X_umap" in adata.obsm, "X_umap missing"
assert "batch" in adata.obs, "batch missing"
assert "cell_type" in adata.obs, "cell_type missing"

adata.write_h5ad(OUTPUT)
print(f"Saved: {OUTPUT}")
DONE_MARKER.touch()
print("=== Step 02d complete ===\n")
