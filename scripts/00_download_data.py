"""
00_download_data.py
Download PBMC 3k (Scanpy built-in, 10x Chromium v1) and a real second-batch PBMC
dataset (GSE132044, 10x Chromium v3 cells; Ding et al. 2020 Nature Biotechnology).
Both datasets have full human-transcriptome gene coverage and a genuine batch effect
arising from different 10x Genomics library-prep chemistry versions.
"""

import gzip
import io
import urllib.request
import scipy.io
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = ROOT / "data" / "raw"
RAW_DIR.mkdir(parents=True, exist_ok=True)
DONE_MARKER = ROOT / "results" / "00_download.done"
DONE_MARKER.parent.mkdir(parents=True, exist_ok=True)

GEO_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE132nnn/GSE132044/suppl"
PLATFORM = "10x-Chromium-v3"   # which cells to keep from GSE132044
N_TARGET = 3000                  # subsample to this many cells
RANDOM_STATE = 42


def download_gz(url: str) -> bytes:
    """Download a .gz URL and return decompressed bytes."""
    print(f"  Downloading {url.split('/')[-1]} ...", flush=True)
    with urllib.request.urlopen(url, timeout=120) as resp:
        compressed = resp.read()
    return gzip.decompress(compressed)


# ──────────────────────────────────────────────
# Batch 1: PBMC 3k (Scanpy built-in)
# ──────────────────────────────────────────────
print("=== Step 00: Downloading data ===")
print("\n[Batch 1] PBMC 3k (10x Chromium v1) ...")
adata_3k = sc.datasets.pbmc3k()
adata_3k.obs["batch"] = "pbmc3k"
adata_3k.write_h5ad(RAW_DIR / "pbmc3k_raw.h5ad")
print(f"  Saved: {adata_3k.shape[0]} cells × {adata_3k.shape[1]} genes")

# ──────────────────────────────────────────────
# Batch 2: GSE132044 – 10x Chromium v3 PBMCs
# Ding et al. 2020, Nature Biotechnology
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132044
# ──────────────────────────────────────────────
print(f"\n[Batch 2] GSE132044 – {PLATFORM} PBMCs (Ding et al. 2020) ...")

# 1. Cell barcodes (small)
cell_bytes = download_gz(f"{GEO_BASE}/GSE132044_pbmc_hg38_cell.tsv.gz")
all_cells = cell_bytes.decode().strip().splitlines()

# 2. Gene list (small)
gene_bytes = download_gz(f"{GEO_BASE}/GSE132044_pbmc_hg38_gene.tsv.gz")
raw_gene_names = gene_bytes.decode().strip().splitlines()
# Format is "ENSGXXX_SYMBOL" — extract the symbol part
gene_symbols = [g.split("_", 1)[1] if "_" in g else g for g in raw_gene_names]
print(f"  Genes: {len(gene_symbols)}")

# 3. Count matrix (genes × cells, MTX format)
mtx_bytes = download_gz(f"{GEO_BASE}/GSE132044_pbmc_hg38_count_matrix.mtx.gz")
matrix = scipy.io.mmread(io.BytesIO(mtx_bytes)).tocsc()   # genes × cells (CSC for column slicing)
print(f"  Full matrix: {matrix.shape[0]} genes × {matrix.shape[1]} cells")

# 4. Filter to chosen platform
platform_mask = np.array([f".{PLATFORM}." in c for c in all_cells])
platform_indices = np.where(platform_mask)[0]
print(f"  {PLATFORM} cells found: {len(platform_indices)}")

# Slice to platform cells (CSC → subset columns efficiently)
mat_platform = matrix[:, platform_indices].T.tocsr()    # cells × genes
cells_platform = [all_cells[i] for i in platform_indices]

# Build AnnData
adata_geo = ad.AnnData(
    X=mat_platform,
    obs=pd.DataFrame(index=cells_platform),
    var=pd.DataFrame(index=gene_symbols),
)
adata_geo.var_names_make_unique()   # handle any duplicate gene symbols
adata_geo.obs_names_make_unique()

# 5. Subsample to N_TARGET cells
np.random.seed(RANDOM_STATE)
n_keep = min(N_TARGET, adata_geo.n_obs)
sc.pp.subsample(adata_geo, n_obs=n_keep, random_state=RANDOM_STATE)
adata_geo.obs["batch"] = "pbmc68k"   # keep original batch-key value per spec

adata_geo.write_h5ad(RAW_DIR / "pbmc68k_raw.h5ad")
print(f"  Saved: {adata_geo.shape[0]} cells × {adata_geo.shape[1]} genes")
print(f"  Source: GSE132044 ({PLATFORM}), doi:10.1038/s41587-020-0465-8")

DONE_MARKER.touch()
print("\n=== Step 00 complete ===\n")
