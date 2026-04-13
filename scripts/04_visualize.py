"""
04_visualize.py
Generate UMAP plots per method (colored by batch + cell type)
and a metric summary heatmap.
Output: results/figures/
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from pathlib import Path

ROOT       = Path(__file__).resolve().parents[1]
CORR_DIR   = ROOT / "results" / "corrected"
METRICS_DIR = ROOT / "results" / "metrics"
FIG_DIR    = ROOT / "results" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)
DONE_MARKER = ROOT / "results" / "04_visualize.done"

sc.set_figure_params(dpi=150, figsize=(5, 5), fontsize=12)

print("=== Step 04: Visualization ===")

method_files = {
    "Uncorrected": ROOT / "data" / "processed" / "merged_preprocessed.h5ad",
    "BBKNN":       CORR_DIR / "bbknn.h5ad",
    "Harmony":     CORR_DIR / "harmony.h5ad",
    "Seurat CCA":  CORR_DIR / "seurat.h5ad",
}

available = {k: v for k, v in method_files.items() if v.exists()}
print(f"Methods available: {list(available.keys())}")

# --- Per-method UMAP panels ---
for method, path in available.items():
    adata = sc.read_h5ad(path)

    # Use uncorrected UMAP key for baseline
    umap_key = "X_umap_uncorrected" if method == "Uncorrected" and "X_umap_uncorrected" in adata.obsm else "X_umap"
    if umap_key not in adata.obsm:
        print(f"  No UMAP found for {method}, skipping plots.")
        continue

    adata.obsm["X_umap"] = adata.obsm[umap_key]  # ensure sc.pl uses correct key

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(f"Integration: {method}", fontsize=14, fontweight="bold")

    sc.pl.umap(adata, color="batch",      ax=axes[0], show=False,
               title="Colored by batch",      frameon=False, legend_loc="right margin")
    sc.pl.umap(adata, color="cell_type",  ax=axes[1], show=False,
               title="Colored by cell type",  frameon=False, legend_loc="right margin")

    plt.tight_layout()
    outpath = FIG_DIR / f"umap_{method.replace(' ', '_').lower()}.png"
    plt.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {outpath}")

# --- Side-by-side comparison grid ---
n_methods = len(available)
fig, axes = plt.subplots(n_methods, 2, figsize=(12, 5 * n_methods))
if n_methods == 1:
    axes = axes[np.newaxis, :]  # ensure 2D

for row_idx, (method, path) in enumerate(available.items()):
    adata = sc.read_h5ad(path)
    umap_key = "X_umap_uncorrected" if method == "Uncorrected" and "X_umap_uncorrected" in adata.obsm else "X_umap"
    if umap_key not in adata.obsm:
        continue
    adata.obsm["X_umap"] = adata.obsm[umap_key]

    sc.pl.umap(adata, color="batch",     ax=axes[row_idx, 0], show=False,
               title=f"{method} | batch",     frameon=False, legend_loc="right margin")
    sc.pl.umap(adata, color="cell_type", ax=axes[row_idx, 1], show=False,
               title=f"{method} | cell type", frameon=False, legend_loc="right margin")

plt.tight_layout()
comp_path = FIG_DIR / "umap_comparison_grid.png"
plt.savefig(comp_path, dpi=150, bbox_inches="tight")
plt.close()
print(f"  Saved: {comp_path}")

# --- Metric heatmap ---
scores_path = METRICS_DIR / "scores.csv"
if scores_path.exists():
    scores_df = pd.read_csv(scores_path)
    metric_cols = [c for c in ["iLISI", "NMI", "ARI"] if c in scores_df.columns]

    if metric_cols:
        heatmap_data = scores_df.set_index("method")[metric_cols].astype(float)

        # Normalize each metric to 0–1 for visual comparison
        norm_data = (heatmap_data - heatmap_data.min()) / (heatmap_data.max() - heatmap_data.min() + 1e-8)

        fig, axes = plt.subplots(1, 2, figsize=(14, max(4, len(heatmap_data) * 0.8 + 2)))

        # Raw scores
        sns.heatmap(heatmap_data, annot=True, fmt=".3f", cmap="YlOrRd",
                    ax=axes[0], cbar=True, linewidths=0.5)
        axes[0].set_title("Raw Scores", fontweight="bold")
        axes[0].set_xlabel("")

        # Normalized scores
        sns.heatmap(norm_data, annot=True, fmt=".3f", cmap="YlOrRd",
                    ax=axes[1], cbar=True, linewidths=0.5)
        axes[1].set_title("Normalized (0–1 per metric)", fontweight="bold")
        axes[1].set_xlabel("")

        plt.suptitle("scRNA-seq Batch Correction Benchmarking\nMetric Summary",
                     fontsize=13, fontweight="bold", y=1.02)
        plt.tight_layout()
        heatmap_path = FIG_DIR / "metric_heatmap.png"
        plt.savefig(heatmap_path, dpi=150, bbox_inches="tight")
        plt.close()
        print(f"  Saved: {heatmap_path}")
else:
    print("  No scores.csv found — skipping heatmap.")

# --- Generate simple HTML report ---
html_dir = ROOT / "report"
html_dir.mkdir(exist_ok=True)

scores_table = ""
if scores_path.exists():
    df = pd.read_csv(scores_path)
    scores_table = df.to_html(index=False, classes="scores-table", border=0, float_format="{:.4f}".format)

figure_blocks = ""
for method in available.keys():
    fname = f"umap_{method.replace(' ', '_').lower()}.png"
    fpath = FIG_DIR / fname
    if fpath.exists():
        rel_path = f"../results/figures/{fname}"
        figure_blocks += f"""
        <div class="method-block">
            <h3>{method}</h3>
            <img src="{rel_path}" alt="{method} UMAP" />
        </div>
        """

heatmap_block = ""
if (FIG_DIR / "metric_heatmap.png").exists():
    heatmap_block = '<img src="../results/figures/metric_heatmap.png" alt="Metric Heatmap" style="max-width:800px;" />'

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>scRNA-seq Batch Correction Benchmarking Report</title>
<style>
  body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
          max-width: 1100px; margin: 0 auto; padding: 2rem; color: #222; background: #fafafa; }}
  h1   {{ color: #1a1a2e; border-bottom: 2px solid #e0e0e0; padding-bottom: 0.5rem; }}
  h2   {{ color: #2c3e50; margin-top: 2rem; }}
  h3   {{ color: #34495e; }}
  .method-block {{ margin: 1.5rem 0; padding: 1rem; background: white;
                   border-radius: 8px; box-shadow: 0 1px 4px rgba(0,0,0,0.1); }}
  .method-block img {{ max-width: 100%; border-radius: 4px; }}
  .scores-table {{ border-collapse: collapse; width: 100%; }}
  .scores-table th, .scores-table td {{ padding: 0.6rem 1rem; text-align: left;
                                         border-bottom: 1px solid #e0e0e0; }}
  .scores-table th {{ background: #2c3e50; color: white; }}
  .scores-table tr:hover {{ background: #f0f4f8; }}
  .tag {{ display: inline-block; padding: 2px 8px; border-radius: 12px;
          font-size: 0.8rem; background: #e8f4fd; color: #2980b9; margin: 0 4px; }}
</style>
</head>
<body>
<h1>scRNA-seq Batch Correction Benchmarking Report</h1>
<p>
  <span class="tag">PBMC 3k</span>
  <span class="tag">PBMC 68k</span>
  <span class="tag">Methods: {", ".join(available.keys())}</span>
</p>

<h2>Metric Summary</h2>
{scores_table}

<h2>Metric Heatmap</h2>
{heatmap_block}

<h2>UMAP Projections</h2>
{figure_blocks}

<h2>Comparison Grid</h2>
<img src="../results/figures/umap_comparison_grid.png"
     alt="All methods comparison" style="max-width:100%;" />

<hr>
<p style="color:#888; font-size:0.85rem;">
Generated by scrna-batch-benchmark pipeline &middot;
Metrics: iLISI (batch mixing), NMI &amp; ARI (biological conservation vs cell type labels)
</p>
</body>
</html>
"""

report_path = html_dir / "summary.html"
report_path.write_text(html)
print(f"  Saved report: {report_path}")

DONE_MARKER.touch()
print("=== Step 04 complete ===\n")
