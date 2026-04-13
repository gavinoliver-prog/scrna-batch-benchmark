#!/usr/bin/env Rscript
# 02c_integrate_seurat.R
# Seurat CCA-based anchor integration.
# Reads two h5ad files (via SeuratDisk), runs FindIntegrationAnchors + IntegrateData,
# writes integrated PCA embedding back to a temp CSV for Python to consume.

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
input_3k  <- args[1]   # path to pbmc3k_raw.h5ad
input_68k <- args[2]   # path to pbmc68k_raw.h5ad
output_embed <- args[3] # path to write integrated PCA CSV
output_meta  <- args[4] # path to write cell metadata CSV

cat("=== Seurat CCA Integration (R) ===\n")

load_h5ad_as_seurat <- function(h5ad_path) {
  # Convert h5ad -> h5seurat -> Seurat object
  h5seurat_path <- sub("\\.h5ad$", ".h5seurat", h5ad_path)
  Convert(h5ad_path, dest = "h5seurat", overwrite = TRUE)
  obj <- LoadH5Seurat(h5seurat_path)
  return(obj)
}

cat("Loading data...\n")
seurat_3k  <- load_h5ad_as_seurat(input_3k)
seurat_68k <- load_h5ad_as_seurat(input_68k)

seurat_3k$batch  <- "pbmc3k"
seurat_68k$batch <- "pbmc68k"

# Preprocess each object independently
preprocess_seurat <- function(obj) {
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst",
                              nfeatures = 2000, verbose = FALSE)
  return(obj)
}

cat("Preprocessing...\n")
seurat_3k  <- preprocess_seurat(seurat_3k)
seurat_68k <- preprocess_seurat(seurat_68k)

obj_list <- list(pbmc3k = seurat_3k, pbmc68k = seurat_68k)

# Find integration anchors (CCA)
cat("Finding integration anchors (CCA)...\n")
anchors <- FindIntegrationAnchors(
  object.list = obj_list,
  anchor.features = 2000,
  reduction = "cca",
  dims = 1:30,
  verbose = FALSE
)

# Integrate
cat("Integrating...\n")
integrated <- IntegrateData(anchorset = anchors, dims = 1:30, verbose = FALSE)
DefaultAssay(integrated) <- "integrated"

# Scale + PCA on integrated assay
integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)

# Write outputs
pca_embed <- Embeddings(integrated, reduction = "pca")
meta      <- integrated@meta.data[, c("batch"), drop = FALSE]

cat(paste("Writing PCA embedding:", output_embed, "\n"))
write.csv(pca_embed, output_embed, quote = FALSE)

cat(paste("Writing metadata:", output_meta, "\n"))
write.csv(meta, output_meta, quote = FALSE)

cat("=== Seurat R script complete ===\n")
