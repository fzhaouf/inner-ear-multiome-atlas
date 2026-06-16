# Export Seurat RNA data and metadata to matrix/CSV files for Python/AnnData.
#
# Can be used either as a sourced function or from command line:
#   Rscript utilities/seurat_to_mtx_custom.R path/to/object.rds output_folder

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(colorspace)
})

export_seurat_to_celloracle_mtx <- function(file_path_seurat_object, outdir) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  so <- readRDS(file_path_seurat_object)
  DefaultAssay(so) <- "RNA"

  # RNA normalized data.
  rna_data <- GetAssayData(so, slot = "data", assay = "RNA")
  Matrix::writeMM(obj = rna_data, file = file.path(outdir, "assay_RNA_data.mtx"))
  write.csv(colnames(rna_data), file = file.path(outdir, "assay_RNA_cells.csv"))
  write.csv(rownames(rna_data), file = file.path(outdir, "assay_RNA_genes.csv"))

  # RNA raw counts.
  rna_counts <- GetAssayData(so, slot = "counts", assay = "RNA")
  Matrix::writeMM(obj = rna_counts, file = file.path(outdir, "assay_RNA_rawdata.mtx"))

  # Metadata.
  meta <- so@meta.data
  meta$active_ident <- Idents(so)
  write.csv(meta, file = file.path(outdir, "meta_data.csv"))

  # Metadata dtypes.
  meta_dtypes <- sapply(meta, class)
  write.csv(data.frame(dtype = meta_dtypes), file = file.path(outdir, "meta_data_dtype.csv"))

  # Variable features.
  var_genes <- VariableFeatures(so)
  write.csv(var_genes, file = file.path(outdir, "var_genes.csv"))

  # Dimensional reductions.
  for (reduction_name in names(so@reductions)) {
    coords <- Embeddings(so[[reduction_name]])
    safe_name <- gsub("\\.", "_", reduction_name)
    write.csv(coords, file = file.path(outdir, paste0("reduction_", safe_name, ".csv")))
  }

  # Cluster colors for plotting, if needed.
  clusters <- Idents(so)
  n_cluster <- length(levels(clusters))
  hues <- seq(15, 375, length = n_cluster + 1)
  colors_hex <- hcl(h = hues, l = 65, c = 100)[1:n_cluster]
  color_df <- data.frame(colors_hex, row.names = levels(clusters))
  write.csv(color_df, file = file.path(outdir, "cluster_color_hex.csv"))

  invisible(outdir)
}

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2) {
    stop("Usage: Rscript seurat_to_mtx_custom.R path/to/object.rds output_folder")
  }
  export_seurat_to_celloracle_mtx(args[1], args[2])
}
