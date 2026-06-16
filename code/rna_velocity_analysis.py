# RNA velocity analysis with scVelo.


from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv

scv.settings.verbosity = 3
scv.settings.set_figure_params("scvelo", facecolor="white", dpi=100, frameon=False)

PROJECT_DIR = Path(".")
ANNDATA_DIR = PROJECT_DIR / "anndata"
FIG_DIR = PROJECT_DIR / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

INTEGRATED_H5AD = ANNDATA_DIR / "integrated_anndata.h5ad"
OUTPUT_H5AD = ANNDATA_DIR / "integrated_anndata_velocity.h5ad"
LINEAGE_OUTPUT_H5AD = ANNDATA_DIR / "integrated_anndata_velocity_lineage_subset.h5ad"

LOOM_FILES = {
    "day45": PROJECT_DIR / "D45.loom",
    "day50": PROJECT_DIR / "D50.loom",
    "day55": PROJECT_DIR / "D55.loom",
}


def clean_loom_barcodes(loom_adata: ad.AnnData, day_label: str) -> ad.AnnData:
    """Convert velocyto loom barcodes to match Seurat/AnnData barcodes."""
    barcodes = []
    for raw_bc in loom_adata.obs_names:
        # Example velocyto barcode pattern: sample:AAAC...x
        bc = raw_bc.split(":")[-1].replace("x", "")
        barcodes.append(f"{day_label}_{bc}-1")
    loom_adata.obs_names = barcodes
    return loom_adata


def load_and_merge_looms() -> ad.AnnData:
    """Load day-specific loom files and concatenate spliced/unspliced layers."""
    loom_objects = []
    for day_label, loom_file in LOOM_FILES.items():
        ldata = scv.read(str(loom_file), cache=True)
        ldata = clean_loom_barcodes(ldata, day_label)
        loom_objects.append(ldata)

    ldata_merged = ad.concat(
        loom_objects,
        join="outer",
        label="loom_sample",
        keys=list(LOOM_FILES.keys()),
        index_unique=None,
    )
    ldata_merged.var_names_make_unique()
    return ldata_merged


# ---------------------------------------------------------------------
# Load integrated AnnData and loom matrices
# ---------------------------------------------------------------------

adata = sc.read_h5ad(INTEGRATED_H5AD)
ldata = load_and_merge_looms()

# Merge spliced/unspliced layers into the integrated AnnData object.
adata = scv.utils.merge(adata, ldata)

if "X_wnn_umap" in adata.obsm and "X_umap" not in adata.obsm:
    adata.obsm["X_umap"] = adata.obsm["X_wnn_umap"]

# ---------------------------------------------------------------------
# scVelo preprocessing and velocity inference
# ---------------------------------------------------------------------

scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=3000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

# Stochastic velocity was used for the full dataset overview.
scv.tl.velocity(adata, mode="stochastic")
scv.tl.velocity_graph(adata, n_jobs=8)
scv.tl.velocity_pseudotime(adata)

adata.write(OUTPUT_H5AD, compression="gzip")

# Overview velocity plots.
scv.pl.velocity_embedding_stream(
    adata,
    basis="umap",
    color="cell_type" if "cell_type" in adata.obs else "wsnn_clusters",
    density=3,
    smooth=0.5,
    title="",
    save="velocity_embedding_stream_full.pdf",
)

scv.pl.velocity_embedding_grid(
    adata,
    basis="umap",
    color="cell_type" if "cell_type" in adata.obs else "wsnn_clusters",
    title="",
    scale=0.25,
    save="velocity_embedding_grid_full.pdf",
)

# Gene-level velocity examples.
for gene in ["JAG1", "POU4F3", "ATOH1", "MYO7A"]:
    if gene in adata.var_names:
        scv.pl.velocity(adata, var_names=[gene], color="wsnn_clusters", save=f"velocity_{gene}.pdf")

# ---------------------------------------------------------------------
# Hair-cell lineage subset for dynamical velocity and CellRank
# ---------------------------------------------------------------------

lineage_labels = [
    "Otic",
    "otic",
    "Otic prosensory progenitors",
    "Early hair cells",
    "early_hair_cell",
    "EHC",
    "Late hair cells",
    "late_hair_cell",
    "LHC",
]

if "cell_type" in adata.obs:
    lineage_mask = adata.obs["cell_type"].isin(lineage_labels)
elif "wsnn_clusters" in adata.obs:
    lineage_mask = adata.obs["wsnn_clusters"].astype(str).isin(["1", "12", "19"])
else:
    raise ValueError("No cell_type or wsnn_clusters column found for lineage subsetting.")

adata_lineage = adata[lineage_mask].copy()

# Dynamical model was used for the focused hair-cell lineage analysis.
scv.tl.recover_dynamics(adata_lineage, n_jobs=8)
scv.tl.velocity(adata_lineage, mode="dynamical")
scv.tl.velocity_graph(adata_lineage, show_progress_bar=False)
scv.tl.velocity_pseudotime(adata_lineage)

scv.pl.velocity_embedding_stream(
    adata_lineage,
    basis="umap",
    color="cell_type" if "cell_type" in adata_lineage.obs else "wsnn_clusters",
    density=2.5,
    smooth=0.8,
    title="",
    save="velocity_embedding_stream_lineage.pdf",
)

adata_lineage.write(LINEAGE_OUTPUT_H5AD, compression="gzip")
