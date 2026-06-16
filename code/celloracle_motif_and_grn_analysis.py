# CellOracle motif scanning and gene regulatory network inference.


from pathlib import Path

import celloracle as co
from celloracle import motif_analysis as ma
from celloracle.applications import Pseudotime_calculator
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc

PROJECT_DIR = Path(".")
CELLORACLE_DIR = PROJECT_DIR / "celloracle"
DATA_DIR = CELLORACLE_DIR / "data"
TMP_EXPORT_DIR = CELLORACLE_DIR / "tmp_export"
CICERO_DIR = CELLORACLE_DIR / "cicero"
FIG_DIR = PROJECT_DIR / "figures"

DATA_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------
# Convert Seurat matrix export to AnnData
# ---------------------------------------------------------------------

def make_anndata_from_seurat_export(folder: Path, output_h5ad: Path) -> sc.AnnData:
    """Create AnnData from matrix/metadata files exported from Seurat."""
    adata = sc.read_mtx(folder / "assay_RNA_data.mtx").T

    genes = pd.read_csv(folder / "assay_RNA_genes.csv", index_col=0).iloc[:, 0].values
    cells = pd.read_csv(folder / "assay_RNA_cells.csv", index_col=0).iloc[:, 0].values
    adata.var_names = genes
    adata.obs_names = cells

    meta = pd.read_csv(folder / "meta_data.csv", index_col=0)
    adata.obs = meta.loc[adata.obs_names]

    raw_counts = sc.read_mtx(folder / "assay_RNA_rawdata.mtx").T.X
    adata.layers["raw_count"] = raw_counts

    var_genes = pd.read_csv(folder / "var_genes.csv", index_col=0).iloc[:, 0].values
    adata.var["variable_gene"] = adata.var_names.isin(var_genes)

    for redfile in folder.glob("reduction_*.csv"):
        redname = redfile.name.replace("reduction_", "").replace(".csv", "")
        coords = pd.read_csv(redfile, index_col=0).loc[adata.obs_names].values
        adata.obsm[f"X_{redname}"] = coords

    # Standardize UMAP key for CellOracle.
    if "X_wnn.umap" in adata.obsm and "X_umap" not in adata.obsm:
        adata.obsm["X_umap"] = adata.obsm["X_wnn.umap"]
    if "X_wnn_umap" in adata.obsm and "X_umap" not in adata.obsm:
        adata.obsm["X_umap"] = adata.obsm["X_wnn_umap"]

    adata.write_h5ad(output_h5ad, compression="gzip")
    return adata


adata_file = DATA_DIR / "lineage_subset_celloracle.h5ad"
adata = make_anndata_from_seurat_export(TMP_EXPORT_DIR, adata_file)

# Optional outlier removal from the otic progenitor region on the UMAP.
if "cell_type" in adata.obs and "X_umap" in adata.obsm:
    emb = adata.obsm["X_umap"]
    otic = (adata.obs["cell_type"] == "otic").values
    if otic.sum() > 0:
        otic_idx = otic.nonzero()[0]
        center = pd.DataFrame(emb[otic]).median(axis=0).values
        distances = ((emb[otic] - center) ** 2).sum(axis=1) ** 0.5
        threshold = pd.Series(distances).quantile(0.95)
        keep_otic_small = distances < threshold
        keep_otic_full = pd.Series(False, index=range(adata.n_obs)).values
        keep_otic_full[otic_idx] = keep_otic_small
        keep = (~otic) | (otic & keep_otic_full)
        adata = adata[keep].copy()

# ---------------------------------------------------------------------
# Add pseudotime for CellOracle development-flow analysis
# ---------------------------------------------------------------------

sc.pp.pca(adata, n_comps=50)
sc.pp.neighbors(adata, n_neighbors=30, use_rep="X_pca")
sc.tl.diffmap(adata)

pt = Pseudotime_calculator(
    adata=adata,
    obsm_key="X_umap",
    cluster_column_name="cell_type",
)

pt.set_lineage(
    lineage_dictionary={
        "Otic_to_Hair": ["otic", "early_hair_cell", "late_hair_cell"]
    }
)

root_candidates = adata.obs_names[adata.obs["cell_type"] == "otic"]
root_cell = root_candidates[min(100, len(root_candidates) - 1)]
pt.set_root_cells(root_cells={"Otic_to_Hair": root_cell})
pt.get_pseudotime_per_each_lineage()
adata.obs["Pseudotime"] = pt.adata.obs["Pseudotime"]

adata.write_h5ad(DATA_DIR / "lineage_subset_celloracle_with_pseudotime.h5ad", compression="gzip")

# ---------------------------------------------------------------------
# Process Cicero peak co-accessibility data
# ---------------------------------------------------------------------

ref_genome = "hg38"
peaks = pd.read_csv(CICERO_DIR / "all_peaks.csv", index_col=0).x.values
cicero_connections = pd.read_csv(CICERO_DIR / "cicero_connections.csv", index_col=0)

tss_annotated = ma.get_tss_info(peak_str_list=peaks, ref_genome=ref_genome)
integrated = ma.integrate_tss_peak_with_cicero(
    tss_peak=tss_annotated,
    cicero_connections=cicero_connections,
)

peak_gene_links = integrated[integrated.coaccess >= 0.8]
peak_gene_links = peak_gene_links[["peak_id", "gene_short_name"]].reset_index(drop=True)
peak_gene_links.to_csv(DATA_DIR / "processed_peak_file.csv")

# ---------------------------------------------------------------------
# Motif scan and base GRN creation
# ---------------------------------------------------------------------

genome_installation = ma.is_genome_installed(ref_genome=ref_genome, genomes_dir=None)
if not genome_installation:
    import genomepy

    genomepy.install_genome(name=ref_genome, provider="UCSC", genomes_dir=None)

peak_data = pd.read_csv(DATA_DIR / "processed_peak_file.csv", index_col=0)
peak_data = ma.check_peak_format(peak_data, ref_genome, genomes_dir=None)

tfi = ma.TFinfo(peak_data_frame=peak_data, ref_genome=ref_genome, genomes_dir=None)
tfi.scan(fpr=0.02, motifs=None, verbose=True)
tfi.to_hdf5(file_path=str(DATA_DIR / "motif_scan.celloracle.tfinfo"))

tfi.reset_filtering()
tfi.filter_motifs_by_score(threshold=10)
tfi.make_TFinfo_dataframe_and_dictionary(verbose=True)
base_grn = tfi.to_dataframe()
base_grn.to_parquet(DATA_DIR / "base_GRN_dataframe.parquet")

# ---------------------------------------------------------------------
# CellOracle object and cell-type-specific GRNs
# ---------------------------------------------------------------------

adata = sc.read_h5ad(DATA_DIR / "lineage_subset_celloracle_with_pseudotime.h5ad")
adata = adata[:, adata.var["variable_gene"]].copy()
adata.X = adata.layers["raw_count"].copy()

oracle = co.Oracle()
oracle.import_anndata_as_raw_count(
    adata=adata,
    cluster_column_name="cell_type",
    embedding_name="X_umap",
)
oracle.import_TF_data(TF_info_matrix=base_grn)

oracle.perform_PCA()
explained = oracle.pca.explained_variance_ratio_
n_comps = min(50, max(10, int((explained.cumsum() < 0.9).sum())))

n_cell = oracle.adata.shape[0]
k = int(0.025 * n_cell)
k = max(k, 10)
oracle.knn_imputation(
    n_pca_dims=n_comps,
    k=k,
    balanced=True,
    b_sight=k * 8,
    b_maxl=k * 4,
    n_jobs=4,
)

oracle.to_hdf5(str(DATA_DIR / "TANG.celloracle.oracle"))

links = oracle.get_links(
    cluster_name_for_GRN_unit="cell_type",
    alpha=10,
    verbose_level=10,
)
links.to_hdf5(file_path=str(DATA_DIR / "links.celloracle.links"))

# Quick QC plots.
sc.settings.figdir = str(FIG_DIR)
sc.pl.umap(oracle.adata, color=["cell_type", "FOXI1", "ATOH1"], save="celloracle_qc_umap.pdf")
