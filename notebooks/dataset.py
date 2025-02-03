import numpy as np
import scanpy as sc
import zarr
from scipy.spatial import cKDTree

from gbmhackathon import MosaicDataset

source_dict_mosaic = MosaicDataset.load_tabular()

df_he = source_dict_mosaic["he"]["H1 features"]
df_he.head()


def load_visium_he(subject_id: str) -> sc.AnnData:
    try:
        subject_id_visium = f"{subject_id}_vis"

        adata = MosaicDataset.load_visium(sample_list=[subject_id_visium])[
            subject_id_visium
        ]

        h1_zarr = zarr.open(df_he.loc[subject_id, "path"], mode="r")

        coords_he = h1_zarr["coords"][:] * 224 + 112
        coords_visium = adata.obsm["spatial"]

        tree = cKDTree(coords_he)
        _, indices = tree.query(coords_visium)

        adata.obsm["emb"] = h1_zarr["emb"][indices]

        print(adata)

        return adata
    except Exception as e:
        print(e)
        return None


adata_dict = {subject_id: load_visium_he(subject_id) for subject_id in df_he.index}
adata_dict = {
    subject_id: adata for subject_id, adata in adata_dict.items() if adata is not None
}
