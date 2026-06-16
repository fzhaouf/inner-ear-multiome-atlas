# Single-cell multiome analysis code for developing human vestibular hair cells

This repository contains analysis scripts used for the single-cell multiome atlas of developing human vestibular hair cells. The scripts document the major computational analyses used in the manuscript, including Seurat/Signac integration, RNA velocity, CellRank lineage analysis, and CellOracle perturbation analysis. Raw sequencing data and processed intermediate objects may be subject to data-use restrictions and are not included here.

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

```

## Main workflow

The analysis code is organized around four main tasks:

**Multiome preprocessing and integration**
`code/seurat_multiome_integration.R` contains the Seurat/Signac workflow used to process RNA and ATAC data, integrate day 45, day 50, and day 55 samples, and generate WNN UMAP clusters.

**Cell annotation and marker analyses**
`code/cell_annotation_marker_analysis.R` includes cluster annotation, marker gene analysis, differential expression/accessibility analysis, temporal composition summaries, and figure-related plotting code.

**RNA velocity and CellRank**
`code/rna_velocity_analysis.py` and `code/cellrank_lineage_analysis.py` contain the scVelo and CellRank analyses used to infer hair cell lineage progression and identify lineage-associated driver genes.

**CellOracle regulatory analysis**
`code/celloracle_input_preparation.R`, `code/celloracle_motif_and_grn_analysis.py`, and `code/celloracle_perturbation_analysis.py` contain the CellOracle workflow for GRN inference and in silico perturbation analysis, including FOXI1 perturbation.



