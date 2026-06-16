"""CellOracle perturbation simulation.

Purpose:
    Representative analysis code for CellOracle in silico perturbation analysis.
    The main manuscript result focused on FOXI1 overexpression and comparison of
    perturbation vectors with the inferred developmental flow along the otic to
    hair-cell trajectory.
"""

from pathlib import Path

import celloracle as co
from celloracle.applications import Gradient_calculator, Oracle_development_module
import matplotlib.pyplot as plt
import scanpy as sc

PROJECT_DIR = Path(".")
DATA_DIR = PROJECT_DIR / "celloracle" / "data"
FIG_DIR = PROJECT_DIR / "figures" / "celloracle_perturbation"
FIG_DIR.mkdir(parents=True, exist_ok=True)

oracle = co.load_hdf5(str(DATA_DIR / "TANG.celloracle.oracle"))
links = co.load_hdf5(file_path=str(DATA_DIR / "links.celloracle.links"))

# Prepare GRN simulation model using cluster-specific TF dictionaries.
links.filter_links()
oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
oracle.fit_GRN_for_simulation(alpha=10, use_cluster_specific_TFdict=True)

# ---------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------

def plot_gene_expression(gene: str) -> None:
    """Plot imputed expression of selected genes."""
    if gene not in oracle.adata.var_names:
        print(f"{gene} not found in oracle object")
        return
    sc.settings.figdir = str(FIG_DIR)
    sc.pl.umap(
        oracle.adata,
        color=[gene, oracle.cluster_column_name],
        layer="imputed_count",
        use_raw=False,
        cmap="viridis",
        save=f"_{gene}_expression.pdf",
    )


def run_perturbation(
    gene: str,
    value: float,
    label: str,
    n_propagation: int = 3,
    n_grid: int = 40,
    min_mass: float = 5.8,
    scale_simulation: float = 10,
) -> None:
    """Run CellOracle perturbation and save vector-field plots."""
    print(f"Running perturbation: {gene} -> {value} ({label})")

    oracle.simulate_shift(perturb_condition={gene: value}, n_propagation=n_propagation)
    oracle.estimate_transition_prob(n_neighbors=200, knn_random=True, sampled_fraction=1)
    oracle.calculate_embedding_shift(sigma_corr=0.05)

    # Cell-level perturbation vector.
    fig, ax = plt.subplots(1, 2, figsize=(13, 6))
    oracle.plot_quiver(scale=20, ax=ax[0])
    ax[0].set_title(f"Simulated cell identity shift vector: {gene} {label}")
    oracle.plot_quiver_random(scale=20, ax=ax[1])
    ax[1].set_title("Randomized simulation vector")
    fig.savefig(FIG_DIR / f"{gene}_{label}_cell_vectors.pdf", bbox_inches="tight")
    plt.close(fig)

    # Grid-level vector field.
    oracle.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
    oracle.calculate_mass_filter(min_mass=min_mass, plot=True)

    fig, ax = plt.subplots(1, 2, figsize=(13, 6))
    oracle.plot_simulation_flow_on_grid(scale=scale_simulation, ax=ax[0])
    ax[0].set_title(f"Simulated vector field: {gene} {label}")
    oracle.plot_simulation_flow_random_on_grid(scale=scale_simulation, ax=ax[1])
    ax[1].set_title("Randomized vector field")
    fig.savefig(FIG_DIR / f"{gene}_{label}_grid_vectors.pdf", bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 7))
    oracle.plot_cluster_whole(ax=ax, s=10)
    oracle.plot_simulation_flow_on_grid(scale=scale_simulation, ax=ax, show_background=False)
    ax.set_title(f"{gene} {label} perturbation vectors")
    fig.savefig(FIG_DIR / f"{gene}_{label}_vectors_on_clusters.pdf", bbox_inches="tight")
    plt.close(fig)

    compare_with_development_flow(
        gene=gene,
        label=label,
        n_grid=n_grid,
        min_mass=min_mass,
        scale_simulation=scale_simulation,
    )


def compare_with_development_flow(
    gene: str,
    label: str,
    n_grid: int,
    min_mass: float,
    scale_simulation: float,
) -> None:
    """Compare perturbation vectors to pseudotime-derived developmental flow."""
    gradient = Gradient_calculator(oracle_object=oracle, pseudotime_key="Pseudotime")
    gradient.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
    gradient.calculate_mass_filter(min_mass=min_mass, plot=True)
    gradient.transfer_data_into_grid(args={"method": "polynomial", "n_poly": 3}, plot=True)
    gradient.calculate_gradient()

    dev = Oracle_development_module()
    dev.load_differentiation_reference_data(gradient_object=gradient)
    dev.load_perturb_simulation_data(oracle_object=oracle)
    dev.calculate_inner_product()
    dev.calculate_digitized_ip(n_bins=10)

    vm = 0.02
    fig, ax = plt.subplots(1, 2, figsize=(12, 6))
    dev.plot_inner_product_on_grid(vm=vm, s=50, ax=ax[0])
    ax[0].set_title("Perturbation score")
    dev.plot_inner_product_random_on_grid(vm=vm, s=50, ax=ax[1])
    ax[1].set_title("Randomized perturbation score")
    fig.savefig(FIG_DIR / f"{gene}_{label}_perturbation_score.pdf", bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 6))
    dev.plot_inner_product_on_grid(vm=vm, s=50, ax=ax)
    dev.plot_simulation_flow_on_grid(scale=scale_simulation, show_background=False, ax=ax)
    ax.set_title(f"{gene} {label}: perturbation vs development flow")
    fig.savefig(FIG_DIR / f"{gene}_{label}_development_alignment.pdf", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------
# Main perturbation analyses
# ---------------------------------------------------------------------

for marker in ["FOXI1", "ATOH1", "MYO7A"]:
    plot_gene_expression(marker)

# Manuscript-focused perturbation: FOXI1 overexpression.
run_perturbation(
    gene="FOXI1",
    value=0.6,
    label="OE",
    n_propagation=3,
    n_grid=40,
    min_mass=5.8,
    scale_simulation=10,
)

# Optional exploratory perturbations used during analysis.
# Uncomment if needed.
# run_perturbation(gene="FOXI1", value=0.0, label="KO", n_propagation=3, scale_simulation=25)
# run_perturbation(gene="ATOH1", value=0.0, label="KO", n_propagation=3, scale_simulation=5)
