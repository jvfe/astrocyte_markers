library(Seurat)
library(SingleR)
library(tidyverse)
library(patchwork)
library(celldex)

annotate_cell_type <- function(object, ref, results_path) {
  sce <- as.SingleCellExperiment(object)

  pred <- SingleR(sce, ref = ref, labels = ref$label.main)

  score_heats <- plotScoreHeatmap(pred)

  ggsave(
    score_heats,
    filename = paste0(results_path, "score_heatmap.pdf"),
    height = 8,
    width = 12
  )

  my.table <- table(Assigned = pred$pruned.labels,
                    cluster = sce$seurat_clusters)

  pheat <-
    pheatmap::pheatmap(log2(my.table + 10), color = colorRampPalette(c("white", "blue"))(101))

  ggsave(
    pheat,
    filename = paste0(results_path, "cell_pheatmap.pdf"),
    height = 8
  )

  aggregated <- scater::sumCountsAcrossCells(sce,
                                             pred$pruned.labels,
                                             exprs_values = "logcounts",
                                             average = TRUE)

  by_cell_matrix <- assay(aggregated) %>%
    as.data.frame() %>%
    rownames_to_column()

  by_cell_matrix %>%
    write_csv(paste0(results_path, "per_cell_expression.csv"))

  return(by_cell_matrix)
}

# Refs
mm_ref <- MouseRNAseqData()
ref <- BlueprintEncodeData()

# Mus musculus
mm <- readRDS("results/mmusculus/integrated.rds")
mm_cells <- annotate_cell_type(mm, mm_ref, "results/mmusculus/")

# Lonchura striata domestica
lstr <- readRDS("results/lstriata/integrated.rds")
lstr_cells <- annotate_cell_type(lstr, ref, "results/lstriata/")

# Pogona vitticeps
pv <- readRDS("results/pvitticeps/integrated.rds")
pv_cells <- annotate_cell_type(pv, ref, "results/pvitticeps/")
