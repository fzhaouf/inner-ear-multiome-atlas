# Seurat/Signac multiome preprocessing and integration


suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(EnsDb.Hsapiens.v86)
  library(GenomicRanges)
  library(dplyr)
  library(ggplot2)
})

project_dir <- Sys.getenv("PROJECT_DIR", unset = ".")
data_dir <- file.path(project_dir, "Data")
out_dir <- file.path(project_dir, "processed_data")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

samples <- data.frame(
  sample_id = c("D45", "D50", "D55"),
  day_label = c("day45", "day50", "day55"),
  stringsAsFactors = FALSE
)

# -------------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------------

sample_path <- function(sample_id, ...) {
  file.path(data_dir, sample_id, ...)
}

read_peak_granges <- function(sample_id) {
  peaks <- read.table(
    file = sample_path(sample_id, "atac_peaks.bed"),
    col.names = c("chr", "start", "end")
  )
  makeGRangesFromDataFrame(peaks)
}

make_combined_peaks <- function(sample_ids) {
  peak_granges <- do.call(c, lapply(sample_ids, read_peak_granges))
  combined_peaks <- reduce(peak_granges)
  peakwidths <- width(combined_peaks)
  combined_peaks[peakwidths < 10000 & peakwidths > 20]
}

create_multiome_object <- function(sample_id, day_label, combined_peaks) {
  message("Creating Seurat multiome object for ", day_label)

  fragment_file <- sample_path(sample_id, "atac_fragments.tsv.gz")
  feature_dir <- sample_path(sample_id, "filtered_feature_bc_matrix")

  fragments <- CreateFragmentObject(path = fragment_file)
  atac_counts <- FeatureMatrix(fragments = fragments, features = combined_peaks)

  multiome <- Read10X(data.dir = feature_dir)
  rna_counts <- multiome$`Gene Expression`

  common_cells <- intersect(colnames(rna_counts), colnames(atac_counts))
  rna_counts <- rna_counts[, common_cells]
  atac_counts <- atac_counts[, common_cells]

  obj <- CreateSeuratObject(counts = rna_counts)
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

  # Keep standard chromosomes for ATAC peaks.
  grange_counts <- StringToGRanges(rownames(atac_counts), sep = c("-", "-"))
  grange_use <- seqnames(grange_counts) %in% standardChromosomes(grange_counts)
  atac_counts <- atac_counts[as.vector(grange_use), ]

  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  seqlevelsStyle(annotations) <- "UCSC"
  genome(annotations) <- "hg38"

  chrom_assay <- CreateChromatinAssay(
    counts = atac_counts,
    sep = c("-", "-"),
    genome = "hg38",
    fragments = fragment_file,
    annotation = annotations
  )

  obj[["ATAC"]] <- chrom_assay
  obj$experiment_day <- day_label
  obj <- RenameCells(obj, add.cell.id = day_label)
  obj
}

qc_filter <- function(obj) {
  subset(
    obj,
    subset = nCount_ATAC < 1e5 &
      nCount_ATAC > 1e2 &
      nCount_RNA < 30000 &
      nCount_RNA > 1000 &
      percent.mt < 20
  )
}

preprocess_rna_atac <- function(obj) {
  DefaultAssay(obj) <- "RNA"
  obj <- SCTransform(obj, verbose = TRUE)

  DefaultAssay(obj) <- "ATAC"
  obj <- FindTopFeatures(obj, min.cutoff = 10)
  obj <- RunTFIDF(obj)
  obj <- RunSVD(obj)

  DefaultAssay(obj) <- "SCT"
  obj
}

# -------------------------------------------------------------------------
# Build per-sample objects
# -------------------------------------------------------------------------

combined_peaks <- make_combined_peaks(samples$sample_id)
saveRDS(combined_peaks, file.path(out_dir, "combined_peaks_hg38.rds"))

seurat_list <- Map(
  f = function(sample_id, day_label) {
    obj <- create_multiome_object(sample_id, day_label, combined_peaks)
    obj <- qc_filter(obj)
    obj <- preprocess_rna_atac(obj)
    saveRDS(obj, file.path(out_dir, paste0(day_label, ".multiome.seurat.rds")))
    obj
  },
  sample_id = samples$sample_id,
  day_label = samples$day_label
)
names(seurat_list) <- samples$day_label

# -------------------------------------------------------------------------
# RNA integration using SCT/rPCA
# -------------------------------------------------------------------------

features <- SelectIntegrationFeatures(seurat_list, nfeatures = 3000)
seurat_list <- PrepSCTIntegration(seurat_list, anchor.features = features)
seurat_list <- lapply(seurat_list, RunPCA, features = features, verbose = FALSE)

integration_anchors <- FindIntegrationAnchors(
  object.list = seurat_list,
  normalization.method = "SCT",
  anchor.features = features,
  dims = 1:30,
  reduction = "rpca",
  k.anchor = 20
)

integrated <- IntegrateData(
  anchorset = integration_anchors,
  normalization.method = "SCT",
  new.assay.name = "integratedRNA",
  dims = 1:30
)

# -------------------------------------------------------------------------
# ATAC preprocessing and integrated LSI embedding
# -------------------------------------------------------------------------

DefaultAssay(integrated) <- "ATAC"
integrated <- FindTopFeatures(integrated, min.cutoff = "q25")
integrated <- RunTFIDF(integrated)
integrated <- RunSVD(integrated)

integrated_atac <- IntegrateEmbeddings(
  anchorset = integration_anchors,
  new.reduction.name = "integratedLSI",
  reductions = integrated@reductions$lsi
)
integrated@reductions$integratedLSI <- integrated_atac@reductions$integratedLSI

# -------------------------------------------------------------------------
# WNN UMAP and clustering
# -------------------------------------------------------------------------

DefaultAssay(integrated) <- "integratedRNA"
integrated <- ScaleData(integrated)
integrated <- RunPCA(integrated)

integrated <- FindMultiModalNeighbors(
  integrated,
  reduction.list = list("pca", "integratedLSI"),
  dims.list = list(1:50, 2:50)
)

integrated <- RunUMAP(
  integrated,
  nn.name = "weighted.nn",
  reduction.name = "wnn.umap",
  reduction.key = "wnnUMAP_"
)

integrated <- FindClusters(
  integrated,
  graph.name = "wsnn",
  algorithm = 1,
  resolution = 0.1,
  verbose = FALSE
)

integrated$wsnn_clusters <- integrated$seurat_clusters

# Optional RNA-only UMAP/clustering used for comparison.
DefaultAssay(integrated) <- "integratedRNA"
integrated <- FindNeighbors(integrated, reduction = "pca", dims = 1:30)
integrated <- RunUMAP(
  integrated,
  reduction = "pca",
  dims = 1:50,
  reduction.name = "umap.rna",
  reduction.key = "rnaUMAP_"
)
integrated <- FindClusters(
  integrated,
  graph.name = "integratedRNA_snn",
  resolution = 0.8
)
integrated$rna_clusters <- integrated$seurat_clusters

# Restore WNN cluster identities for downstream analyses.
Idents(integrated) <- integrated$wsnn_clusters

saveRDS(integrated, file.path(out_dir, "integrated_multiome_seurat.rds"))
save(integrated, file = file.path(out_dir, "integrated_multiome_seurat.rda"))

# Basic visualization examples.
p_wnn <- DimPlot(
  integrated,
  reduction = "wnn.umap",
  group.by = "wsnn_clusters",
  label = TRUE,
  repel = TRUE
) + ggtitle("WNN clusters")

ggsave(file.path(out_dir, "wnn_umap_clusters.pdf"), p_wnn, width = 7, height = 6)
