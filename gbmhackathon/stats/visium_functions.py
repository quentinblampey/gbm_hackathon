"""Script containing functions to perform statistics on Visium data."""

from typing import Optional, Union

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc


def perform_multi_clustering(
    dataset: dict[str, ad.AnnData],  # Change to a dictionary of AnnData objects
    sample_list: Union[list[str], None] = None,
    cl_method: str = "louvain",
    cl_input: str = "X_pca",  ## to clarify
    repres_n_neighbors: int = 15,
    repres_n_pca: int = 50,
    resolution: Optional[list] = None,
    random_state: Optional[int] = 42,
    layer: str = "CPM",
    plot: bool = True,
    save_output: Optional[str] = None,
    **kwargs: dict,
) -> dict[str, ad.AnnData]:  # pylint: disable=PLR0912
    """Perform clustering for each sample in the dataset (dictionary of AnnData)."""
    sample_list = sample_list or list(
        dataset.keys()
    )  # Get the keys of the dictionary as the sample IDs
    resolution = resolution or [0.25, 0.5, 0.75, 1]

    adats = {}

    valid_methods = ["louvain", "leiden"]
    valid_layers = [
        "CPM",
        "log_CPM",
        "raw",
    ]
    valid_inputs = ["X", "X_pca"]

    if cl_method not in valid_methods:
        raise ValueError(f"Method should be one of {valid_methods}")

    if layer not in valid_layers:
        raise ValueError(f"Layer should be one of {valid_layers}")

    if cl_input not in valid_inputs:
        raise ValueError(f"cl_input should be one of {valid_inputs}")

    for samp in sample_list:  # Iterate over sample_list (sample IDs)
        print(f"Processing sample {samp}")
        # Get the AnnData object directly from the dictionary
        adata = dataset[samp]

        # Define layer to use for PCA or Clustering
        adata.X = adata.layers[layer]
        # Computes the nearest neighbors distance matrix and a neighborhood graph of
        # observations
        print("    - Computing neighbors")
        if cl_input == "X":
            sc.pp.neighbors(
                adata,
                n_neighbors=repres_n_neighbors,
                use_rep="X",
                random_state=random_state,
            )
        elif cl_input == "X_pca":
            print(f"Shape of AnnData X {adata.X.shape}")
            n_comps = (
                np.min(adata.X.shape) - 1
                if np.min(adata.X.shape) < repres_n_pca
                else repres_n_pca
            )

            sc.pp.scale(adata)
            sc.pp.pca(adata, n_comps=n_comps, random_state=random_state)
            sc.pp.neighbors(
                adata,
                n_neighbors=repres_n_neighbors,
                random_state=random_state,
                use_rep="X_pca",
            )
        print("    - Performing clustering at each resolution level")
        # Perform clustering for each resolution
        for res in resolution:
            key_added = f"multi_cl_{cl_method}_{res}"
            # Perform clustering
            if cl_method == "louvain":
                sc.tl.louvain(
                    adata,
                    resolution=res,
                    random_state=random_state,
                    key_added=key_added,
                    **kwargs,
                )
            elif cl_method == "leiden":
                sc.tl.leiden(
                    adata,
                    resolution=res,
                    random_state=random_state,
                    key_added=key_added,
                    **kwargs,
                )
        if plot:
            if save_output is not None:
                sc.pl.spatial(
                    adata,
                    color=[i for i in adata.obs.columns if i.startswith("multi_cl_")],
                    save=f"{samp}_{save_output}",
                    show=plot,
                    **kwargs,
                )
            else:
                sc.pl.spatial(
                    adata,
                    color=[i for i in adata.obs.columns if i.startswith("multi_cl_")],
                )
        adats[samp] = adata
    return adats


def quantify_cell_population_activity(
    dataset: dict[str, ad.AnnData],
    biomarker_dict: dict[str, list[str]],
    layer: str = "CPM",
) -> dict[str, ad.AnnData]:
    """Quantify the activity of specific cell populations based on their biomarkers and store the result in the .obsm slot.

    Parameters
    ----------
    dataset : dict[str, ad.AnnData]
        A dictionary where the keys are sample IDs and the values are AnnData objects.
    biomarker_dict : dict[str, list[str]], optional
        A dictionary where the keys are cell types (e.g., 'Fibrosis', 'Lymphocytes', etc.)
        and the values are lists of gene names representing biomarkers for those cell types.
        If None, default biomarker sets for GBM are used.

    Returns
    -------
    dict[str, ad.AnnData]
        Updated dictionary of AnnData objects with a new `.obsm` slot containing the quantified activity for each cell population.
    """
    if biomarker_dict is None:
        print(
            "No biomarker dictionary provided. You need to provide a dictionary with cell types and their biomarkers."
        )

    for samp, adata in dataset.items():
        print(f"Processing sample {samp}")

        # Create a dictionary to store the activity for each cell type
        activity_dict = {}

        # Iterate over each cell type and their corresponding biomarkers
        for cell_type, biomarkers in biomarker_dict.items():
            # Ensure biomarkers exist in the gene names of the dataset
            valid_genes = [gene for gene in biomarkers if gene in adata.var_names]

            if valid_genes:
                # Sum (or average) the expression of the biomarkers for the current cell type
                activity = adata[:, valid_genes].layers[layer].sum(axis=1).A.flatten()
                activity_dict[cell_type] = activity
            else:
                print(f"Warning: No valid genes found for {cell_type} in sample {samp}")

        # Convert activity_dict to a DataFrame and store it in the .obsm slot
        activity_df = pd.DataFrame(activity_dict)

        # Ensure the DataFrame is aligned with the number of cells
        activity_df.index = (
            adata.obs_names
        )  # Make sure the index matches the cells in the dataset

        # Store the activity_df in the .obsm slot of the AnnData object
        adata.obsm["cell_population_activity"] = activity_df
        adata.obsm["cell_population_activity_normalized"] = activity_df.div(
            activity_df.sum(axis=1), axis=0
        )

    return dataset
