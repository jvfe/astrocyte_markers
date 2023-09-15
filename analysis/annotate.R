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
library(celldex)
library(org.Mm.eg.db)

sc <- readRDS("data/mmusculus/integrated_data.rds")

sce <- as.SingleCellExperiment(sc)

rownames(sce) <-
  stringr::str_remove(rownames(sce), pattern = "\\.\\d+")

ref <- MouseRNAseqData()

pred <- SingleR(sce, ref = ref, labels = ref$label.main)

table(pred$labels)
plotScoreHeatmap(pred)

my.table <- table(Assigned = pred$pruned.labels,
                  cluster = sce$seurat_clusters)

pheatmap::pheatmap(log2(my.table + 10), color = colorRampPalette(c("white", "blue"))(101))

aggregated <- scater::sumCountsAcrossCells(sce,
                                           pred$pruned.labels,
                                           exprs_values = "logcounts",
                                           average = TRUE)

matrix <- assay(aggregated) %>%
  as.data.frame() %>%
  rownames_to_column()

matrix %>%
  write_csv("results/mmusculus/per_cell_expression.csv")
