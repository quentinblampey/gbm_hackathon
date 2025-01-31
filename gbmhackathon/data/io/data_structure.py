"""Define data structures for tabular data."""

from pathlib import Path
from typing import NamedTuple, Optional, Union

from gbmhackathon.data.io.loaders import BaseDataLoader


class DataFile(NamedTuple):
    """Define the structure of a datafile."""

    name: Optional[Union[str, Path]] = None
    loader: Optional[BaseDataLoader] = None

    def load(self):
        """Load datafile.

        Returns
        -------
        pd.DataFrame
            Loaded and cleaned datafile for a given source.

        Raises
        ------
        ValueError
            In case no loader has been specified.
        """
        if self.loader is not None:
            return self.loader.load_data(self.name)
        raise ValueError("Please specify a loader to load data.")


class DataSheet(NamedTuple):
    """Define the structure of an excel sheet."""

    name: str
    header: Optional[tuple[int]] = None


class DataSource(NamedTuple):
    """Define the structure of a data source (tabular modality)."""

    files: Union[DataFile, dict[str, DataFile]]
