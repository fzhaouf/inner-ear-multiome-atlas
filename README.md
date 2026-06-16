# Single-cell multiome analysis code for developing human vestibular hair cells

This repository contains representative analysis scripts used for the single-cell multiome atlas of developing human vestibular hair cells. The goal of this code release is transparency: the scripts document the major computational analyses used in the manuscript, including Seurat/Signac integration, RNA velocity, CellRank lineage analysis, and CellOracle perturbation analysis.

These scripts are not intended to be a fully automated end-to-end pipeline. File paths, object names, and data-access steps may need to be adapted for local use. Raw sequencing data and processed intermediate objects may be subject to data-use restrictions and are not included here.

## Repository structure

```text
code/
  seurat_multiome_integration.R
  cell_annotation_marker_analysis.R
  rna_velocity_analysis.py
  cellrank_lineage_analysis.py
  celloracle_input_preparation.R
  celloracle_motif_and_grn_analysis.py
  celloracle_perturbation_analysis.py

utilities/
  merge_loom.py
  seurat_to_mtx_custom.R

environment/
  python_environment.yml
  R_sessionInfo_template.R

data/
  README_data_access.md
```

## Main workflow

### Multiome preprocessing and integration

`code/seurat_multiome_integration.R` creates Seurat objects for day 45, day 50, and day 55 10x Multiome samples. It creates a unified ATAC peak set, re-quantifies fragments, performs RNA and ATAC preprocessing, integrates samples using reciprocal PCA/SCT-based integration, integrates ATAC LSI embeddings, and performs WNN UMAP/clustering.

### Cell annotation and marker analyses

`code/cell_annotation_marker_analysis.R` adds cell-type annotations to WNN clusters, generates marker gene plots, calculates cluster markers, runs differential expression and differential accessibility analyses, summarizes temporal composition across sampling days, and produces chromatin accessibility plots for genes such as FOXI1, ATOH1, POU4F3, and MYO7A.

### RNA velocity and CellRank

`code/rna_velocity_analysis.py` merges spliced/unspliced loom files with the integrated AnnData object, computes RNA velocity using scVelo, generates velocity stream plots, and creates a lineage subset for otic progenitor, early hair cell, and late hair cell populations.

`code/cellrank_lineage_analysis.py` uses the velocity-inferred AnnData object to estimate terminal/initial states, compute fate probabilities, reconstruct lineage relationships, and identify lineage driver genes.

### CellOracle regulatory network and perturbation analysis

`code/celloracle_input_preparation.R` prepares RNA and ATAC inputs from the integrated Seurat object, exports RNA matrices/metadata for AnnData conversion, extracts lineage cells, and creates Cicero peak co-accessibility inputs.

`code/celloracle_motif_and_grn_analysis.py` processes Cicero peak links, performs motif scanning with CellOracle, creates the base GRN, imports the lineage AnnData object into CellOracle, and estimates cell-type-specific GRNs.

`code/celloracle_perturbation_analysis.py` performs in silico perturbation simulations, including FOXI1 overexpression, and compares simulated perturbation vectors with the inferred developmental flow.

## Notes

- This code assumes a human hg38 reference genome.
- Input data paths are represented as local placeholders and should be changed before running.
- The code reflects the major analysis tasks used in the study rather than a polished software package.
- For reproducibility, users should record exact package versions after running `sessionInfo()` in R and `conda env export` or equivalent in Python.

