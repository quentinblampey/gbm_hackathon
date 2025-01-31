"""Pre-processor of tabular data."""

from abc import ABC, abstractmethod

import numpy as np
import pandas as pd
from pathlib import Path


# Opt-in to the future behavior
pd.set_option("future.no_silent_downcasting", True)


class BaseDataCleaner(ABC):
    """Clean Tabular Data. Class to remove typos, preprocess data columns.

    Attributes
    ----------
    dataset_name : str

    Methods
    -------
    transform: Transform the input data.
    """

    @abstractmethod
    def transform(self, x):
        """Transform the input data.

        Parameters
        ----------
        x : pd.DataFrame
            Input data

        Returns
        -------
        pd.DataFrame
            x with no typos
        """
        pass


class MOSAICLightClinicalDataCleaner(BaseDataCleaner):
    """Add sample_id as index and remove columns with only nan values.

    Useful only for the original clinical data from MOSAIC.

    Methods
    -------
    transform: Transform the input data.
    """

    def transform(self, x):
        """Transform the input data.

        Note that this method acts **inplace** on the input data.

        Parameters
        ----------
        x : pd.DataFrame
            Input data

        Returns
        -------
        pd.DataFrame
            x with no typos
        """
        # Add sample_id as index
        x.set_index("sample_id", inplace=True)

        # Remove columns with only nan values
        x.dropna(axis=1, how="all", inplace=True)

        return x


class MOSAICClinicalDataCleaner(BaseDataCleaner):
    """Clean processed clinical data from MOSAIC.

    Remove typos, preprocess data columns.

    Methods
    -------
    transform: Transform the input data.
    """

    def transform(self, x):
        """Transform the input data.

        Note that this method acts **inplace** on the input data.

        Parameters
        ----------
        x : pd.DataFrame
            Input data

        Returns
        -------
        pd.DataFrame
            x with no typos
        """
        # Add sample_id as index
        x.set_index("sample_id", inplace=True)

        # Replace Unknown values by nan
        x.replace("Unknown", np.nan, inplace=True)

        # Remove columns with only nan values
        x.dropna(axis=1, how="all", inplace=True)

        return x


class MOSAICKeyEventsDataCleaner(BaseDataCleaner):
    """Clean Key Events data file from MOSAIC.

    Remove typos, preprocess data columns.

    Methods
    -------
    transform: Transform the input data.
    """

    def transform(self, x):
        """Transform the input data.

        Note that this method acts **inplace** on the input data.

        Parameters
        ----------
        x : pd.DataFrame
            Input data

        Returns
        -------
        pd.DataFrame
            x with no typos
        """
        # Add Patient ID as index
        x.set_index("patient_id", inplace=True)

        # Replace Unknown values by nan
        x.replace("Unknown", np.nan, inplace=True)

        # Remove columns with only nan values
        x.dropna(axis=1, how="all", inplace=True)

        return x


class MOSAICBulkRNADataCleaner(BaseDataCleaner):
    """Clean Key Events data file from MOSAIC.

    Remove typos, preprocess data columns.

    Methods
    -------
    transform: Transform the input data.
    """

    def transform(self, x):
        """Transform the input data.

        Note that this method acts **inplace** on the input data.

        Parameters
        ----------
        x : pd.DataFrame
            Input data

        Returns
        -------
        pd.DataFrame
            x with no typos
        """
        x.rename(columns={"gene": "EnsemblID"}, inplace=True)
        x.set_index("EnsemblID", inplace=True)

        # If weird ensembl IDs exist, drop them
        # Appears in the deseq normalized counts for instance
        x = x[~x.index.str.contains("g@-ext")]  # corresponds to
        # immunoglobulin super-loci, mainly used to find fusion genes
        # and therefore not necessary for DGEA

        return x


class MOSAICWESDataCleaner(BaseDataCleaner):
    """Clean WES data file from MOSAIC.

    Remove typos, preprocess data columns.

    Methods
    -------
    transform: Transform the input data.
    """

    def transform(self, x):
        """Transform the input data.

        Note that this method acts **inplace** on the input data.

        Parameters
        ----------
        x : pd.DataFrame
            Input data

        Returns
        -------
        pd.DataFrame
            x with no typos
        """
        x = x.dropna(subset=["gene_name"])
        x = x.set_index("gene_name")
        x = x.drop(columns=["ensembl_gene_id"])
        x = x.T  # Transpose the data to have patient id in index
        x.index.name = "sample_id"

        return x


class MOSAICHEDataCleaner(BaseDataCleaner):
    """Clean HE data file from MOSAIC.

    Add patient ID, preprocess data columns.

    Methods
    -------
    transform: Transform the input data.
    """

    def transform(self, x):
        """Transform the input data.

        Note that this method acts **inplace** on the input data.

        Parameters
        ----------
        x : pd.DataFrame
            Input data

        Returns
        -------
        pd.DataFrame
            x with no typos, parsed dates
        """
        x.loc[:, "Subject Id"] = x["path"].apply(lambda v: v.name[:9])
        x.loc[:, "patient id"] = x["Subject Id"].apply(lambda v: v[:-1])
        x = x.set_index("Subject Id")
        return x


class MOSAICHEFeaturesDataCleaner(BaseDataCleaner):
    """Clean HE features data file from MOSAIC.

    Add patient ID, preprocess data columns.

    Methods
    -------
    transform: Transform the input data.
    """

    def transform(self, x):
        """Transform the input data.

        Note that this method acts **inplace** on the input data.

        Parameters
        ----------
        x : pd.DataFrame
            Input data

        Returns
        -------
        pd.DataFrame
            x with no typos, parsed dates
        """
        x.loc[:, "Subject Id"] = x["path"].apply(lambda v: v.name.replace("features_", "")[:9])
        x.loc[:, "patient id"] = x["Subject Id"].apply(lambda v: v[:-1])
        x = x.set_index("Subject Id")
        return x


class BruceLightMetadataDataCleaner(BaseDataCleaner):
    """Clean the medata from Bruce dataset.

    Methods
    -------
    transform: Transform the input data.
    """

    def transform(self, x):
        """Transform the input data.

        Note that this method acts **inplace** on the input data.

        Parameters
        ----------
        x : pd.DataFrame
            Input data

        Returns
        -------
        pd.DataFrame
            x with no typos
        """
        # Drop the "Unnamed: 0" column
        x.drop(columns=["Unnamed: 0"], inplace=True)

        # Add sample ID as index
        x.set_index("sample_id", inplace=True)

        # Replace Unknown values by nan
        x.replace("unknown", np.nan, inplace=True)

        # Remove columns with only nan values
        x.dropna(axis=1, how="all", inplace=True)

        return x


class BruceMIBIImagesCleaner(BaseDataCleaner):
    """Clean MIBI data file from Bruce.

    Add sample ID, preprocess data columns.

    Methods
    -------
    transform: Transform the input data.
    """

    def transform(self, x):
        """Transform the input data.

        Note that this method acts **inplace** on the input data.

        Parameters
        ----------
        x : pd.DataFrame
            Input data

        Returns
        -------
        pd.DataFrame
            x with no typos, parsed dates
        """
        x.loc[:, "sample_id"] = x["path"].apply(lambda v: Path(v).parts[6])
        x.loc[:, "marker"] = x["path"].apply(lambda v: Path(v).stem)
        x = x.set_index("sample_id")
        return x
