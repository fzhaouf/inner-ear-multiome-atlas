# Prepare CellOracle RNA and ATAC inputs from integrated Seurat object
#
# Purpose:
#   Representative code for preparing the hair-cell lineage subset and peak
#   co-accessibility inputs for CellOracle. The RNA subset is exported to matrix
#   files using utilities/seurat_to_mtx_custom.R and then converted to AnnData
#   by the companion Python code in celloracle_motif_and_grn_analysis.py.

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(monocle3)
  library(cicero)
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(Matrix)
})

project_dir <- Sys.getenv("PROJECT_DIR", unset = ".")
out_dir <- file.path(project_dir, "celloracle", "data")
tmp_export_dir <- file.path(project_dir, "celloracle", "tmp_export")
cicero_dir <- file.path(project_dir, "celloracle", "cicero")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tmp_export_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(cicero_dir, showWarnings = FALSE, recursive = TRUE)

seurat_file <- file.path(project_dir, "processed_data", "integrated_multiome_seurat_annotated.rds")
obj <- readRDS(seurat_file)

# -------------------------------------------------------------------------
# RNA layer: subset otic -> early hair cell -> late hair cell lineage
# -------------------------------------------------------------------------

DefaultAssay(obj) <- "RNA"
Idents(obj) <- obj$wsnn_clusters

# Harmonize lineage labels for CellOracle.
obj$celloracle_cell_type <- "other"
obj$celloracle_cell_type[obj$wsnn_clusters == "1"] <- "otic"
obj$celloracle_cell_type[obj$wsnn_clusters == "19"] <- "early_hair_cell"
obj$celloracle_cell_type[obj$wsnn_clusters == "12"] <- "late_hair_cell"

# Normalize RNA and identify variable genes before export.
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000)

lineage_cells <- WhichCells(
  obj,
  expression = celloracle_cell_type %in% c("otic", "early_hair_cell", "late_hair_cell")
)
lineage_obj <- subset(obj, cells = lineage_cells)

# Keep only the information needed by CellOracle.
lineage_obj[["ATAC"]] <- NULL
lineage_obj[["SCT"]] <- NULL
if ("integratedRNA" %in% names(lineage_obj@assays)) {
  lineage_obj[["integratedRNA"]] <- NULL
}
DefaultAssay(lineage_obj) <- "RNA"

# Keep WNN UMAP as the primary embedding.
keep_reductions <- intersect(names(lineage_obj@reductions), c("wnn.umap", "umap.rna"))
lineage_obj@reductions <- lineage_obj@reductions[keep_reductions]

lineage_obj$cell_type <- lineage_obj$celloracle_cell_type
Idents(lineage_obj) <- "cell_type"

keep_meta <- intersect(
  colnames(lineage_obj@meta.data),
  c("nCount_RNA", "nFeature_RNA", "percent.mt", "experiment_day", "cell_type", "wsnn_clusters")
)
lineage_obj@meta.data <- lineage_obj@meta.data[, keep_meta, drop = FALSE]

lineage_rds <- file.path(out_dir, "lineage_subset_RNAonly.rds")
saveRDS(lineage_obj, lineage_rds)

# Export matrices for conversion to AnnData.
# Equivalent command-line usage:
# Rscript utilities/seurat_to_mtx_custom.R celloracle/data/lineage_subset_RNAonly.rds celloracle/tmp_export
source(file.path(project_dir, "utilities", "seurat_to_mtx_custom.R"))
export_seurat_to_celloracle_mtx(lineage_rds, tmp_export_dir)

# -------------------------------------------------------------------------
# ATAC layer: create Cicero CDS and peak co-accessibility links
# -------------------------------------------------------------------------

obj <- readRDS(seurat_file)
DefaultAssay(obj) <- "ATAC"

# Peak x cell binary counts.
atac_counts <- GetAssayData(obj, slot = "counts", assay = "ATAC")
atac_counts@x[atac_counts@x > 0] <- 1

cellinfo <- data.frame(cells = colnames(obj))
rownames(cellinfo) <- cellinfo$cells

# Convert peak names from chr-start-end to chr_start_end for Cicero metadata.
rownames(atac_counts) <- gsub("-", "_", rownames(atac_counts))
peakinfo <- data.frame(
  site_name = rownames(atac_counts),
  chr = sapply(strsplit(rownames(atac_counts), "_"), `[`, 1),
  bp1 = as.integer(sapply(strsplit(rownames(atac_counts), "_"), `[`, 2)),
  bp2 = as.integer(sapply(strsplit(rownames(atac_counts), "_"), `[`, 3))
)
rownames(peakinfo) <- peakinfo$site_name

colnames(atac_counts) <- rownames(cellinfo)

input_cds <- new_cell_data_set(
  atac_counts,
  cell_metadata = cellinfo,
  gene_metadata = peakinfo
)

input_cds <- monocle3::detect_genes(input_cds)
input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0, ]

set.seed(2017)
input_cds <- estimate_size_factors(input_cds)
input_cds <- preprocess_cds(input_cds, method = "LSI")
input_cds <- reduce_dimension(input_cds, reduction_method = "UMAP", preprocess_method = "LSI")

umap_coords <- reducedDims(input_cds)$UMAP
cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap_coords)
saveRDS(cicero_cds, file.path(out_dir, "cicero_cds.rds"))

# hg38 chromosome sizes for Cicero.
genome <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)
genome <- genome[!grepl("_", names(genome))]
genome <- genome[!grepl("random", names(genome))]
hg38_genome <- data.frame(V1 = names(genome), V2 = as.numeric(genome))

cicero_connections <- run_cicero(cicero_cds, genomic_coords = hg38_genome)
write.csv(cicero_connections, file.path(cicero_dir, "cicero_connections.csv"))

all_peaks <- data.frame(x = rownames(atac_counts))
write.csv(all_peaks, file.path(cicero_dir, "all_peaks.csv"))
