"""Hyperparameters and definitions for the GBM MOSAIC data on abstra."""

import os
from pathlib import Path


def generate_data_paths_mosaic(
    clinical_path: Path,
    bulkrna_path: Path,
    wes_path: Path,
    single_cell_path: Path,
    he_path: Path,
    he_features_path: Path,
    visium_path: Path,
):
    """Generate a dictionary of data paths."""
    return {
        "Processed clinical": clinical_path / "GBM_HK_sample_and_clinical_data.csv",
        "Key events clinical": clinical_path / "GBM_HK_multi_entry_events.csv",
        "Treatments": clinical_path / "GBM_HK_multi_entry_treatments.csv",
        "Data dictionary": clinical_path / "GBM_HK_data_dictionary.csv",
        "Bulk RNA raw counts": bulkrna_path / "counts/raw_counts.tsv",
        "Bulk RNA normalized counts": bulkrna_path / "counts/deseq2norm_counts.tsv",
        "Bulk RNA TPM counts": bulkrna_path / "counts/tpm_counts.tsv",
        "Bulk RNA fpkm counts": bulkrna_path / "counts/fpkm_counts.tsv",
        "WES CNV deletion": wes_path
        / "results/cohort_level/cnv/binary_deletion_status.csv",
        "WES CNV amplification": wes_path
        / "results/cohort_level/cnv/binary_duplication_status.csv",
        "WES CNV oncogenic": wes_path
        / "results/cohort_level/cnv/binary_oncogenic_alteration_status.csv",
        "WES mutations": wes_path
        / "results/cohort_level/snv_indel/binary_oncogenic_alteration_status.csv",
        "Single cell anndata": single_cell_path
        / "preprocessed/sc_merged_annotated.h5ad",
        "HE": he_path,
        "HE_features_H1": he_features_path,
        "Visium anndata": visium_path,
    }

def generate_data_paths_bruce(
    global_path: Path
):
    """Generate a dictionary of data paths."""
    return {
        "Metadata": global_path / "death_confirmed_metadata_owkin.csv",
        "MIBI_immune_images": global_path,
        "MIBI_tumor_images": global_path,
    }

HOME_PATH = Path(os.environ["HOME"])
DATA_PATH = HOME_PATH / "SageMaker" / "data"
MOSAIC_PATH = DATA_PATH / "mosaic_dataset"
BRUCE_PATH = DATA_PATH / "bruce_dataset"

DATA_PATH_FILES_MOSAIC = generate_data_paths_mosaic(
    MOSAIC_PATH / "Clinical",
    MOSAIC_PATH / "RNAseq",
    MOSAIC_PATH / "WES",
    MOSAIC_PATH / "single_cell",
    MOSAIC_PATH / "Visium" / "converted_he",
    DATA_PATH / "h1_bioptimus_features",  # H1 features
    MOSAIC_PATH / "Visium" / "spaceranger_count",
)

DATA_PATH_FILES_BRUCE = generate_data_paths_bruce(
    BRUCE_PATH
)
