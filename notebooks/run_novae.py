import novae

from gbmhackathon import MosaicDataset

visium_dict = MosaicDataset.load_visium(
    sample_list=["HK_G_022a_vis", "HK_G_024a_vis", "HK_G_030a_vis"]
)

for sid, adata in visium_dict.items():
    adata.obs["subject_id"] = sid
    adata.obs_names_make_unique()

adatas = list(visium_dict.values())

novae.utils.spatial_neighbors(adatas, technology="visium", slide_key="subject_id")

novae.plot.connectivities(adatas[:12])

model = novae.Novae(adatas, n_hops_local=1, n_hops_view=1)

model.fit(max_epochs=10)

model.compute_representations(adatas)

model.assign_domains(adatas)

novae.plot.domains(adatas, cell_size=200)
