
# Get Integrated Matrix  ---------------------------------------------------------------

library(Seurat)
library(SingleCellExperiment)
library(scDblFinder)
library(preprocessCore)
library(clustree)
library(SingleR)
library(vroom)
library(tibble)
library(tidyverse)
library(patchwork)
library(harmony)

set.seed(1024)

remove_doublets <- function(obj, filename) {
  # Remove doublets using the scDblFinder
  results <- scDblFinder(obj, returnType = 'table') %>%
    as.data.frame() %>%
    dplyr::filter(type == "real")

  # Count how many doublets are
  results %>%
    dplyr::count(class)

  # Keep only the singlets
  keep <- results %>%
    dplyr::filter(class == "singlet") %>%
    rownames()

  med <- Seurat::as.Seurat(obj, counts = "counts", data = "counts")

  med <- med[, keep]

  med[["orig.ident"]] <- med[['sample']]

  saveRDS(med, file = filename)

  return(med)
}

run_qc <- function(obj, mtn_genes, plotname, filename) {

  mtn_genes <- mtn_genes[mtn_genes %in% rownames(GetAssayData(obj))]

  obj[["percent.mt"]] <- PercentageFeatureSet(obj, features = mtn_genes)

  obj[["percent.mt"]][is.na(obj[["percent.mt"]])] <- 0

  med.qc <- FetchData(obj, vars=c("nFeature_originalexp","nCount_originalexp","percent.mt"))

  p <- med.qc %>%
    mutate(keep = if_else(nCount_originalexp > 500 & nFeature_originalexp < 7000 & percent.mt < 10, "keep", "remove")) %>%
    ggplot() +
    geom_point(aes(nCount_originalexp, nFeature_originalexp, colour=keep), alpha=.50) +
    scale_x_log10() +
    scale_y_log10()

  ggsave(p, filename = plotname, width = 15, height = 12)

  obj <- subset(obj, nCount_originalexp > 500 & nFeature_originalexp < 8000 & percent.mt < 10)

  saveRDS(obj, file = "data/mmusculus/filtered_seurat_obj.rds")

  return(obj)
}

integrate_data <- function(obj, plotname, filename) {
  # Data integration using Seurat method

  # Split data according to sample id
  obj.list <- SplitObject(obj, split.by = "sample")

  # Normalize and identify variable features for each sample
  obj.list <- lapply(X = obj.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  })

  # Select features that are variable across different samples
  features <- SelectIntegrationFeatures(object.list = obj.list)

  # Find the anchors across the different samples
  obj.anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features)

  # Create an integrated data assay
  obj.combined <- IntegrateData(anchorset = obj.anchors, k.weight = 50)

  DefaultAssay(obj.combined) <- "integrated"

  # Run the standard workflow for visualization and clustering
  obj.combined <- ScaleData(obj.combined, verbose = FALSE) %>%
    RunPCA(npcs = 30, verbose = FALSE) %>%
    RunUMAP(reduction = "pca", dims = 1:30) %>%
    FindNeighbors(reduction = "pca", dims = 1:30) %>%
    FindClusters(resolution = 0.5)

  # Visualization
  p1 <- DimPlot(obj.combined, reduction = "umap", group.by = "sample")
  p2 <- DimPlot(obj.combined, reduction = "umap", label = TRUE, repel = TRUE)
  p_all <- p1 + p2

  ggsave(p_all, filename = plotname, height = 12, width = 18)

  saveRDS(obj.combined, filename)

  return(obj.combined)
}


# Run Mus musculus ---------------------

mm_mt_gns <- vroom("data/mmusculus/mitchondrial_genes.txt") %>%
  filter(`Gene name` %in% rownames(GetAssayData(mm_seurat))) %>%
  pull(5)

mm <- zellkonverter::readH5AD("data/mmusculus/combined_matrix.h5ad", X_name = "counts")

mm_seurat <- remove_doublets(mm, "results/mmusculus/cleaned_seurat.rds")

mm_filtered <- run_qc(mm_seurat, mm_mt_gns, "results/mmusculus/keep_plot.pdf", "results/mmusculus/filtered_seurat.rds")

mm_integrated <- integrate_data(mm_filtered, "results/mmusculus/integration_plot.pdf", "results/mmusculus/integrated.rds")

# Run Lonchura striata ----------------

ls_mito_genes <-
  c(
    "ATP6",
    "ATP8",
    "COX1",
    "COX2",
    "COX3",
    "CYTB",
    "ND1",
    "ND2",
    "ND3",
    "ND4",
    "ND4L",
    "ND5",
    "ND6"
  )

lstr <- zellkonverter::readH5AD("data/lstriata/combined_matrix.h5ad", X_name = "counts")

ls_seurat <- remove_doublets(lstr, "results/lstriata/cleaned_seurat.rds")

ls_filtered <- run_qc(ls_seurat, ls_mito_genes, "results/lstriata/keep_plot.pdf", "results/lstriata/filtered_seurat.rds")

ls_integrated <- integrate_data(ls_filtered, "results/lstriata/integration_plot.pdf", "results/lstriata/integrated.rds")

# Run Pogona vitticeps ----------------

pv_mito_genes <-
  c(
    "CO1",
    "CO2",
    "CO3",
    "CR1",
    "CR2",
    "AT8",
    "AT6",
    "ND1",
    "ND3",
    "ND4",
    "ND4L",
    "ND5"
  )

pv <- zellkonverter::readH5AD("data/pvitticeps/combined_matrix.h5ad", X_name = "counts")

pv_seurat <- remove_doublets(pv, "results/pvitticeps/cleaned_seurat.rds")

pv_filtered <- run_qc(pv_seurat, pv_mito_genes, "results/pvitticeps/keep_plot.pdf", "results/pvitticeps/filtered_seurat.rds")

pv_integrated <- integrate_data(pv_filtered, "results/pvitticeps/integration_plot.pdf", "results/pvitticeps/integrated.rds")