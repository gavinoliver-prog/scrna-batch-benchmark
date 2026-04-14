#!/usr/bin/env Rscript
# 02c_integrate_seurat.R
# Seurat CCA-based anchor integration.
# Reads MTX-format count matrices (no SeuratDisk — incompatible with Seurat v5),
# runs FindIntegrationAnchors + IntegrateData, writes integrated PCA embedding
# and metadata as CSV for the Python wrapper to consume.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
mtx_dir_3k   <- args[1]   # directory: matrix.mtx, barcodes.tsv, features.tsv
mtx_dir_68k  <- args[2]   # directory: matrix.mtx, barcodes.tsv, features.tsv
output_embed <- args[3]   # path to write integrated PCA CSV
output_meta  <- args[4]   # path to write cell metadata CSV

cat("=== Seurat CCA Integration (R) ===\n")

load_mtx_as_seurat <- function(mtx_dir, project_name) {
  counts   <- readMM(file.path(mtx_dir, "matrix.mtx"))   # genes x cells
  barcodes <- read.table(file.path(mtx_dir, "barcodes.tsv"),
                         header = FALSE, stringsAsFactors = FALSE)$V1
  features <- read.table(file.path(mtx_dir, "features.tsv"),
                         header = FALSE, stringsAsFactors = FALSE)$V1

  rownames(counts) <- make.unique(features)
  colnames(counts) <- barcodes

  CreateSeuratObject(counts = counts, project = project_name)
}

cat("Loading data...\n")
seurat_3k  <- load_mtx_as_seurat(mtx_dir_3k,  "pbmc3k")
seurat_68k <- load_mtx_as_seurat(mtx_dir_68k, "pbmc68k")

seurat_3k$batch  <- "pbmc3k"
seurat_68k$batch <- "pbmc68k"

preprocess_seurat <- function(obj) {
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst",
                              nfeatures = 2000, verbose = FALSE)
  obj
}

cat("Preprocessing...\n")
seurat_3k  <- preprocess_seurat(seurat_3k)
seurat_68k <- preprocess_seurat(seurat_68k)

cat("Finding integration anchors (CCA)...\n")
anchors <- FindIntegrationAnchors(
  object.list     = list(pbmc3k = seurat_3k, pbmc68k = seurat_68k),
  anchor.features = 2000,
  reduction       = "cca",
  dims            = 1:30,
  verbose         = FALSE
)

cat("Integrating...\n")
integrated <- IntegrateData(anchorset = anchors, dims = 1:30, verbose = FALSE)
DefaultAssay(integrated) <- "integrated"

integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)

pca_embed <- Embeddings(integrated, reduction = "pca")
meta      <- integrated@meta.data[, "batch", drop = FALSE]

cat(paste("Writing PCA embedding:", output_embed, "\n"))
write.csv(pca_embed, output_embed, quote = FALSE)

cat(paste("Writing metadata:", output_meta, "\n"))
write.csv(meta, output_meta, quote = FALSE)

cat("=== Seurat R script complete ===\n")
