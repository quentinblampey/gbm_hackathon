"""Script that contains functions to manipulate Visium data in the starter notebook."""

# Import all necessary packages
from typing import Dict

import anndata as ad
import pandas as pd
import scanpy as sc


def normalize_anndata(adata: ad.AnnData, target_sum: float = 1e6) -> ad.AnnData:
    """
    Normalize an AnnData object using CPM normalization and log1p transformation.

    Parameters
    ----------
        adata (ad.AnnData): The AnnData object to normalize.
        target_sum (float): The target sum for the normalization. Default is 1e6.
            if None, the sum of adata.X will be used
            if 1e6 this is CPM normalization

    Returns
    -------
        ad.AnnData: The normalized AnnData object.
    """
    adata.layers["raw"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=target_sum)
    adata.layers["CPM"] = adata.X.copy()
    sc.pp.log1p(adata)
    adata.layers["log_CPM"] = adata.X.copy()
    return adata


def normalize_anndata_wrapper(
    anndata_dict: Dict[str, ad.AnnData], target_sum: float = 1e6
) -> Dict[str, ad.AnnData]:
    """
    Normalize all AnnData objects in a dictionary.

    Parameters
    ----------
        anndata_dict (Dict[str, ad.AnnData]): A dictionary of AnnData objects.
        target_sum (float): The target sum for the normalization. Default is 1e4.

    Returns
    -------
        Dict[str, ad.AnnData]: The dictionary with normalized AnnData objects.
    """
    for sample_id, adata in anndata_dict.items():
        try:
            anndata_dict[sample_id] = normalize_anndata(adata, target_sum=target_sum)
        except Exception as e:
            print(f"Failed to normalize data for sample {sample_id}: {e}")
    return anndata_dict


def convert_obsm_to_adata(adata: ad.AnnData, obsm_key: str) -> ad.AnnData:
    """Convert obsm matrices to adata.

    This function converts an obsm matrix to an AnnData object in order to use it with liana functions.

    Parameters
    ----------
    adata : ad.AnnData
        AnnData object containing the obsm matrix.
    obsm_key : str
        Key of the obsm matrix to convert to adata.

    Returns
    -------
    ad.AnnData
        AnnData object with obsm matrix as X.

    Raises
    ------
    KeyError
        If obsm_key is not in adata.obsm.
    """
    if obsm_key not in adata.obsm.keys():
        raise KeyError(f"No key with name '{obsm_key}' in adata.obsm")
    obsm = adata.obsm[obsm_key]
    obsm_df = pd.DataFrame(obsm, index=adata.obs.index)
    obsm_adata = ad.AnnData(
        obsm_df,
        # obs=adata.obs,
        var=pd.DataFrame(index=obsm_df.columns),
        uns=adata.uns,
        obsm=adata.obsm,
        obsp=adata.obsp,
    )
    return obsm_adata
