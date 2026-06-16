# Cell annotation, marker analysis, differential expression, and differential accessibility
#
# Purpose:
#   Representative analysis code for annotating WNN clusters, finding marker genes,
#   comparing early hair cells with related populations, and visualizing chromatin
#   accessibility around key hair-cell genes.

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(dplyr)
  library(ggplot2)
  library(openxlsx)
})

project_dir <- Sys.getenv("PROJECT_DIR", unset = ".")
out_dir <- file.path(project_dir, "processed_data")
fig_dir <- file.path(project_dir, "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

obj <- readRDS(file.path(out_dir, "integrated_multiome_seurat.rds"))
Idents(obj) <- obj$wsnn_clusters

# -------------------------------------------------------------------------
# Cluster annotation
# -------------------------------------------------------------------------

cluster_annotations <- c(
  "0" = "Supporting cells",
  "1" = "Otic prosensory progenitors",
  "2" = "Neural progenitors",
  "3" = "Stroma cells",
  "4" = "Mesenchymal cells",
  "5" = "Epithelial cells",
  "6" = "Neural progenitors",
  "7" = "Schwann cells",
  "8" = "Stroma cells",
  "9" = "Myogenic precursors",
  "10" = "Epithelial cells",
  "11" = "Early otic progenitors",
  "12" = "Late hair cells",
  "13" = "Stroma cells",
  "14" = "Neural progenitors",
  "15" = "Dark cells",
  "16" = "Nerve",
  "17" = "Melanocytes",
  "18" = "Epithelial cells",
  "19" = "Early hair cells"
)

obj$cell_type <- cluster_annotations[as.character(obj$wsnn_clusters)]
Idents(obj) <- obj$cell_type

saveRDS(obj, file.path(out_dir, "integrated_multiome_seurat_annotated.rds"))

p_celltype <- DimPlot(
  obj,
  reduction = "wnn.umap",
  group.by = "cell_type",
  label = TRUE,
  repel = TRUE,
  label.size = 3
) + NoLegend()

ggsave(file.path(fig_dir, "wnn_umap_cell_types.pdf"), p_celltype, width = 7, height = 6)

# -------------------------------------------------------------------------
# Marker gene plots and marker detection
# -------------------------------------------------------------------------

marker_genes <- c(
  "PAX2", "SOX2", "SIX1", "JAG1",
  "ATOH1", "POU4F3", "MYO7A",
  "FBXO2", "BRICD5", "OC90"
)

DefaultAssay(obj) <- "RNA"

p_dot <- DotPlot(obj, features = marker_genes, assay = "RNA") +
  RotatedAxis() +
  ggtitle("Canonical marker genes")

ggsave(file.path(fig_dir, "marker_gene_dotplot.pdf"), p_dot, width = 9, height = 5)

p_feature <- FeaturePlot(
  obj,
  features = c("FOXI1", "ATOH1", "POU4F3", "MYO7A"),
  reduction = "wnn.umap",
  min.cutoff = "q5",
  max.cutoff = "q95"
)

ggsave(file.path(fig_dir, "key_gene_featureplots.pdf"), p_feature, width = 10, height = 8)

Idents(obj) <- obj$wsnn_clusters
DefaultAssay(obj) <- "RNA"
markers_by_cluster <- FindAllMarkers(obj, only.pos = TRUE)
write.csv(markers_by_cluster, file.path(out_dir, "gene_markers_by_WNN_cluster.csv"))

# -------------------------------------------------------------------------
# Differential expression analyses
# -------------------------------------------------------------------------

DefaultAssay(obj) <- "SCT"
obj <- PrepSCTFindMarkers(obj, assay = "SCT", verbose = TRUE)
Idents(obj) <- obj$wsnn_clusters

# Cluster definitions used in the manuscript analyses.
otic_cluster <- "1"
ehc_cluster <- "19"
lhc_cluster <- "12"
cluster_4 <- "4"
cluster_8 <- "8"

DEgenes_EHC_vs_Otic <- FindMarkers(
  obj,
  ident.1 = ehc_cluster,
  ident.2 = otic_cluster,
  assay = "SCT",
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.1
)
DEgenes_EHC_vs_Otic$p_val_adj <- p.adjust(DEgenes_EHC_vs_Otic$p_val, method = "BH")

DEgenes_EHC_vs_cluster4 <- FindMarkers(
  obj,
  ident.1 = ehc_cluster,
  ident.2 = cluster_4,
  assay = "SCT",
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.1
)
DEgenes_EHC_vs_cluster4$p_val_adj <- p.adjust(DEgenes_EHC_vs_cluster4$p_val, method = "BH")

DEgenes_EHC_vs_cluster8 <- FindMarkers(
  obj,
  ident.1 = ehc_cluster,
  ident.2 = cluster_8,
  assay = "SCT",
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.1
)
DEgenes_EHC_vs_cluster8$p_val_adj <- p.adjust(DEgenes_EHC_vs_cluster8$p_val, method = "BH")

DEgenes_LHC_vs_Otic <- FindMarkers(
  obj,
  ident.1 = lhc_cluster,
  ident.2 = otic_cluster,
  assay = "SCT",
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.1
)
DEgenes_LHC_vs_Otic$p_val_adj <- p.adjust(DEgenes_LHC_vs_Otic$p_val, method = "BH")

write.xlsx(DEgenes_EHC_vs_Otic, file.path(out_dir, "DEgenes_EHC_vs_Otic.xlsx"), rowNames = TRUE)
write.xlsx(DEgenes_EHC_vs_cluster4, file.path(out_dir, "DEgenes_EHC_vs_cluster4.xlsx"), rowNames = TRUE)
write.xlsx(DEgenes_EHC_vs_cluster8, file.path(out_dir, "DEgenes_EHC_vs_cluster8.xlsx"), rowNames = TRUE)
write.xlsx(DEgenes_LHC_vs_Otic, file.path(out_dir, "DEgenes_LHC_vs_Otic.xlsx"), rowNames = TRUE)

# -------------------------------------------------------------------------
# Differential accessibility analyses
# -------------------------------------------------------------------------

DefaultAssay(obj) <- "ATAC"
Idents(obj) <- obj$wsnn_clusters

DApeaks_EHC_vs_Otic <- FindMarkers(
  object = obj,
  ident.1 = ehc_cluster,
  ident.2 = otic_cluster,
  assay = "ATAC",
  test.use = "LR",
  latent.vars = "nCount_ATAC",
  min.pct = 0.1,
  logfc.threshold = 0.25
)

DApeaks_LHC_vs_Otic <- FindMarkers(
  object = obj,
  ident.1 = lhc_cluster,
  ident.2 = otic_cluster,
  assay = "ATAC",
  test.use = "LR",
  latent.vars = "nCount_ATAC",
  min.pct = 0.1,
  logfc.threshold = 0.25
)

write.csv(DApeaks_EHC_vs_Otic, file.path(out_dir, "DApeaks_EHC_vs_Otic.csv"))
write.csv(DApeaks_LHC_vs_Otic, file.path(out_dir, "DApeaks_LHC_vs_Otic.csv"))

# -------------------------------------------------------------------------
# Chromatin accessibility visualization around selected genes
# -------------------------------------------------------------------------

DefaultAssay(obj) <- "ATAC"
Idents(obj) <- obj$wsnn_clusters

plot_coverage <- function(gene, filename) {
  p <- CoveragePlot(
    object = obj,
    region = gene,
    assay = "ATAC",
    idents = c(otic_cluster, cluster_4, cluster_8, lhc_cluster, ehc_cluster),
    annotation = TRUE,
    peaks = TRUE,
    links = FALSE
  )
  ggsave(file.path(fig_dir, filename), p, width = 10, height = 6)
}

plot_coverage("FOXI1", "coverage_FOXI1.pdf")
plot_coverage("ATOH1", "coverage_ATOH1.pdf")
plot_coverage("POU4F3", "coverage_POU4F3.pdf")
plot_coverage("MYO7A", "coverage_MYO7A.pdf")

# -------------------------------------------------------------------------
# Temporal composition by sampling day
# -------------------------------------------------------------------------

temporal_counts <- as.data.frame(table(obj$experiment_day, obj$cell_type))
colnames(temporal_counts) <- c("experiment_day", "cell_type", "n")

temporal_props <- temporal_counts %>%
  group_by(experiment_day) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  ungroup()

write.csv(temporal_props, file.path(out_dir, "cell_type_composition_by_day.csv"), row.names = FALSE)

p_temporal <- ggplot(temporal_props, aes(x = experiment_day, y = percent, fill = cell_type)) +
  geom_col(position = "fill") +
  ylab("Cell-type proportion") +
  xlab("Sampling day") +
  theme_classic()

ggsave(file.path(fig_dir, "cell_type_composition_by_day.pdf"), p_temporal, width = 8, height = 5)
