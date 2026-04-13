"""
03_benchmark.py
Compute scIB batch correction metrics for each integrated method.
Metrics: iLISI (batch mixing), NMI, ARI (biological conservation)
Output: results/metrics/scores.csv
"""

import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path

ROOT       = Path(__file__).resolve().parents[1]
CORR_DIR   = ROOT / "results" / "corrected"
METRICS_DIR = ROOT / "results" / "metrics"
METRICS_DIR.mkdir(parents=True, exist_ok=True)
DONE_MARKER = ROOT / "results" / "03_benchmark.done"

print("=== Step 03: Benchmarking ===")

# Detect which methods completed
method_files = {
    "Uncorrected": ROOT / "data" / "processed" / "merged_preprocessed.h5ad",
    "BBKNN":       CORR_DIR / "bbknn.h5ad",
    "Harmony":     CORR_DIR / "harmony.h5ad",
    "Seurat CCA":  CORR_DIR / "seurat.h5ad",
}

# Try scib-metrics first, fall back to scib
def get_scib():
    try:
        import scib_metrics
        return "scib_metrics"
    except ImportError:
        pass
    try:
        import scib
        return "scib"
    except ImportError:
        return None

scib_backend = get_scib()
if scib_backend is None:
    print("ERROR: Neither scib-metrics nor scib found.")
    print("Install with: pip install scib-metrics")
    raise ImportError("scib not available")

print(f"Using backend: {scib_backend}")

# Map each method to the Leiden key produced by its integration script
LEIDEN_KEY_MAP = {
    "Uncorrected": "leiden_uncorrected",
    "BBKNN":       "leiden_bbknn",
    "Harmony":     "leiden_harmony",
    "Seurat CCA":  "leiden_seurat",
}


def compute_metrics(adata, method_name, use_rep="X_pca"):
    """Compute iLISI, NMI, ARI for a given AnnData."""
    from sklearn.metrics import normalized_mutual_info_score, adjusted_rand_score
    from sklearn.preprocessing import LabelEncoder
    from sklearn.neighbors import NearestNeighbors

    results = {"method": method_name}

    # Use the method-specific Leiden key; recompute on the corrected embedding
    # if the expected key isn't present (e.g. first run or missing).
    leiden_key = LEIDEN_KEY_MAP.get(method_name, "leiden_eval")
    if leiden_key not in adata.obs.columns:
        print(f"  Leiden key '{leiden_key}' not found — recomputing on {use_rep}")
        sc.pp.neighbors(adata, use_rep=use_rep, n_neighbors=15, n_pcs=30,
                        key_added="_eval_neighbors")
        sc.tl.leiden(adata, resolution=0.5, key_added="leiden_eval",
                     neighbors_key="_eval_neighbors")
        leiden_key = "leiden_eval"

    # Encode labels
    le = LabelEncoder()
    true_labels = le.fit_transform(adata.obs["cell_type"].astype(str))
    pred_labels = le.fit_transform(adata.obs[leiden_key].astype(str))

    results["NMI"] = round(normalized_mutual_info_score(true_labels, pred_labels), 4)
    results["ARI"] = round(adjusted_rand_score(true_labels, pred_labels), 4)

    # iLISI — build a NeighborsResults object from the corrected space.
    # BBKNN corrects the graph, not an embedding, so we rebuild a high-k
    # BBKNN graph purely for evaluation rather than using raw PCA.
    try:
        from scib_metrics import ilisi_knn
        from scib_metrics.nearest_neighbors import NeighborsResults

        N_LISI = 90   # neighbors for LISI

        if method_name == "BBKNN":
            # BBKNN corrects graph topology, not the PCA embedding: within-batch
            # cells remain closer in Euclidean distance, so kNN on raw PCA gives
            # the uncorrected result.  The UMAP was projected from the corrected
            # graph, so kNN in UMAP space faithfully captures batch mixing.
            emb = adata.obsm["X_umap"]
            nn = NearestNeighbors(n_neighbors=N_LISI, metric="euclidean",
                                  n_jobs=-1).fit(emb)
            dists, indices = nn.kneighbors(emb)
            nn_results = NeighborsResults(indices=indices, distances=dists)
            ilisi_val = ilisi_knn(nn_results, adata.obs["batch"].values)
            results["iLISI"] = round(float(ilisi_val), 4)
            dist_csr = None
        else:
            # Embedding-based methods: kNN on the corrected embedding
            emb = adata.obsm[use_rep]
            if emb.ndim > 1 and emb.shape[1] > 30:
                emb = emb[:, :30]
            nn = NearestNeighbors(n_neighbors=N_LISI, metric="euclidean",
                                  n_jobs=-1).fit(emb)
            dists, indices = nn.kneighbors(emb)
            nn_results = NeighborsResults(indices=indices, distances=dists)
            ilisi_val = ilisi_knn(nn_results, adata.obs["batch"].values)
            results["iLISI"] = round(float(np.median(ilisi_val)), 4)

    except Exception as e:
        print(f"  iLISI failed ({e}), setting NaN")
        results["iLISI"] = np.nan

    print(f"  {method_name}: NMI={results['NMI']:.4f}, ARI={results['ARI']:.4f}, "
          f"iLISI={results.get('iLISI', 'N/A')}")
    return results

all_scores = []

for method, path in method_files.items():
    if not path.exists():
        print(f"  Skipping {method}: file not found ({path})")
        continue

    print(f"\nEvaluating: {method}")
    adata = sc.read_h5ad(path)

    # Determine representation to use
    if method == "Harmony" and "X_pca_harmony" in adata.obsm:
        use_rep = "X_pca_harmony"
    elif method == "Seurat CCA" and "X_pca_seurat" in adata.obsm:
        use_rep = "X_pca_seurat"
    else:
        use_rep = "X_pca"

    try:
        scores = compute_metrics(adata, method, use_rep=use_rep)
        all_scores.append(scores)
    except Exception as e:
        print(f"  ERROR evaluating {method}: {e}")

scores_df = pd.DataFrame(all_scores)
scores_df = scores_df.sort_values("NMI", ascending=False).reset_index(drop=True)
scores_df.to_csv(METRICS_DIR / "scores.csv", index=False)

print(f"\nScores saved: {METRICS_DIR / 'scores.csv'}")
print(scores_df.to_string(index=False))

DONE_MARKER.touch()
print("=== Step 03 complete ===\n")
