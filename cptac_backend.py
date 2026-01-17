"""
CPTAC Backend - Cancer dataclass and data loading/preprocessing

This module provides the Cancer dataclass with cached data loading for CPTAC datasets.
Data is lazily loaded and cached on first access for efficient reuse.
"""

# ============================================================================
# Patch cptac to handle rate limiting (429) errors gracefully
# This must be done BEFORE importing cptac
# ============================================================================
import os
import sys
import os.path as path

def _patch_cptac_for_offline_fallback():
    """
    Patch requests.get to convert 429/5xx HTTPError into a ConnectionError,
    which cptac handles gracefully by falling back to cached data.
    This must be called BEFORE importing cptac.
    """
    import requests
    from requests.exceptions import ConnectionError as RequestsConnectionError

    _original_get = requests.get

    def _patched_get(*args, **kwargs):
        """Patched requests.get that converts rate limiting errors to ConnectionError."""
        response = _original_get(*args, **kwargs)
        # Check if this is a Zenodo request with rate limiting or server error
        if response.status_code in [429, 500, 502, 503, 504]:
            url = args[0] if args else kwargs.get('url', '')
            if 'zenodo.org' in str(url):
                # Raise ConnectionError which cptac treats as "no internet"
                raise RequestsConnectionError(
                    f"Server unavailable (HTTP {response.status_code}). Using cached data."
                )
        return response

    # Apply the patch
    requests.get = _patched_get

# Apply the patch before importing cptac
_patch_cptac_for_offline_fallback()

import cptac
import pandas as pd
import numpy as np
import logging
from dataclasses import dataclass, field
from typing import Optional, Dict, List, Any


# Registry mapping cancer names to cptac classes (not instances - for lazy loading)
CANCER_REGISTRY: Dict[str, type] = {
    'brca': cptac.Brca,
    'coad': cptac.Coad,
    'hnscc': cptac.Hnscc,
    'luad': cptac.Luad,
    'ovarian': cptac.Ov,
    'ccrcc': cptac.Ccrcc,
    'gbm': cptac.Gbm,
    'lscc': cptac.Lscc,
    'pdac': cptac.Pdac,
}


@dataclass
class Cancer:
    """
    Represents a CPTAC cancer dataset with lazy-loaded and cached data.

    Properties:
        name: Cancer type identifier (e.g., 'brca', 'luad')
        cptac_obj: The underlying cptac cancer object

    Cached Properties (lazily loaded on first access):
        phospho_deduplicated: Deduplicated phosphoproteomics data
        phospho_normalized: Protein-normalized phosphoproteomics data
        proteomics_deduplicated: Deduplicated proteomics data
        clinical_data: Clinical metadata for the cancer cohort (placeholder for future use)
    """
    name: str
    cptac_obj: Any = field(repr=False)

    # Private cache storage
    _phospho_deduplicated: Optional[pd.DataFrame] = field(default=None, repr=False, init=False)
    _phospho_normalized: Optional[pd.DataFrame] = field(default=None, repr=False, init=False)
    _proteomics_deduplicated: Optional[pd.DataFrame] = field(default=None, repr=False, init=False)
    _clinical_data: Optional[pd.DataFrame] = field(default=None, repr=False, init=False)

    @property
    def phospho_deduplicated(self) -> pd.DataFrame:
        """Lazily load and cache deduplicated phosphoproteomics data."""
        if self._phospho_deduplicated is None:
            self._phospho_deduplicated = self._load_deduplicated_phospho()
        return self._phospho_deduplicated

    @property
    def phospho_normalized(self) -> pd.DataFrame:
        """Lazily load and cache protein-normalized phosphoproteomics data."""
        if self._phospho_normalized is None:
            self._phospho_normalized = self._load_normalized_phospho()
        return self._phospho_normalized

    @property
    def proteomics_deduplicated(self) -> pd.DataFrame:
        """Lazily load and cache deduplicated proteomics data."""
        if self._proteomics_deduplicated is None:
            self._proteomics_deduplicated = self._load_deduplicated_proteomics()
        return self._proteomics_deduplicated

    @property
    def clinical_data(self) -> pd.DataFrame:
        """Lazily load and cache clinical data (placeholder for future cohort support)."""
        if self._clinical_data is None:
            self._clinical_data = self._load_clinical_data()
        return self._clinical_data

    def _load_deduplicated_phospho(self) -> pd.DataFrame:
        """
        Load and deduplicate phosphoproteomics data.

        Returns:
            pd.DataFrame: Transposed, deduplicated phosphoproteomics data with samples as columns.
                         Multi-index: (Name/Gene, Site, Peptide, Database_ID)
        """
        logging.info(f"[{self.name}] Loading and deduplicating phosphoproteomics data...")
        phospho = self.cptac_obj.get_phosphoproteomics('bcm').T
        phospho.columns = phospho.columns.values.tolist()
        phospho = phospho.groupby(level=[0, 1, 2, 3]).agg('mean').replace(-np.inf, np.nan).replace(np.inf, np.nan)
        logging.info(f"[{self.name}] Phospho data shape: {phospho.shape} (phosphosites x samples)")
        return phospho

    def _load_normalized_phospho(self) -> pd.DataFrame:
        """
        Load protein-normalized phosphoproteomics data.

        Normalization is performed by subtracting whole-cell proteomics values
        from phosphoproteomics values for each sample.

        Returns:
            pd.DataFrame: Protein-normalized phosphoproteomics data
        """
        logging.info(f"[{self.name}] Calculating protein-normalized phosphoproteomics...")
        phospho = self.phospho_deduplicated  # Uses cached data if available
        whole_cell = self.cptac_obj.get_proteomics('bcm').T
        whole_cell.columns = whole_cell.columns.values.tolist()

        normalized_phospho = phospho.copy()
        samples_normalized = 0
        samples_skipped = 0

        for colname in list(normalized_phospho.columns):
            if colname in whole_cell.columns:
                wc_col = whole_cell[colname]
                phospho_col = phospho[colname]
                merge = pd.merge(phospho_col, wc_col, left_index=True, right_index=True, suffixes=('_phospho', '_wc'))
                merge[colname] = merge[colname + '_phospho'] - merge[colname + '_wc']
                normalized_phospho[colname] = merge[colname]
                samples_normalized += 1
            else:
                normalized_phospho = normalized_phospho.drop(columns=[colname])
                samples_skipped += 1

        logging.info(f"[{self.name}] Normalization complete: {samples_normalized} samples normalized, {samples_skipped} samples removed")
        logging.info(f"[{self.name}] Normalized data shape: {normalized_phospho.shape}")
        return normalized_phospho

    def _load_deduplicated_proteomics(self) -> pd.DataFrame:
        """
        Load and deduplicate proteomics data.

        Returns:
            pd.DataFrame: Transposed, deduplicated proteomics data with samples as columns
        """
        logging.info(f"[{self.name}] Loading and deduplicating proteomics data...")
        proteomics = self.cptac_obj.get_proteomics('bcm').T
        proteomics.columns = proteomics.columns.values.tolist()

        # Check if index is MultiIndex and deduplicate if needed
        if isinstance(proteomics.index, pd.MultiIndex):
            proteomics = proteomics.groupby(level=list(range(proteomics.index.nlevels))).agg('mean').replace(-np.inf, np.nan).replace(np.inf, np.nan)
            logging.info(f"[{self.name}] Deduplicated proteomics data")

        logging.info(f"[{self.name}] Proteomics data shape: {proteomics.shape} (proteins x samples)")
        return proteomics

    def _load_clinical_data(self) -> pd.DataFrame:
        """
        Load clinical metadata from available data sources.

        Tries multiple CPTAC data sources to find clinical data.

        Returns:
            pd.DataFrame: Clinical data for the cancer cohort
        """
        logging.info(f"[{self.name}] Loading clinical data...")

        # Try each potential source for clinical data
        sources = ['mssm', 'bcm', 'broad', 'umich', 'washu', 'harmonized']

        for source in sources:
            try:
                clinical = self.cptac_obj.get_clinical(source=source)
                if clinical is not None and not clinical.empty:
                    logging.info(f"[{self.name}] Clinical data loaded from '{source}': {clinical.shape}")
                    return clinical
            except Exception as e:
                logging.debug(f"[{self.name}] No clinical data from '{source}': {e}")
                continue

        logging.warning(f"[{self.name}] Could not load clinical data from any source")
        return pd.DataFrame()

    def clear_cache(self) -> None:
        """Clear all cached data for this cancer object."""
        self._phospho_deduplicated = None
        self._phospho_normalized = None
        self._proteomics_deduplicated = None
        self._clinical_data = None
        logging.info(f"[{self.name}] Cache cleared")

    # Placeholder methods for future cohort support
    def get_sample_ids(self, tumor_only: bool = False, normal_only: bool = False) -> List[str]:
        """
        Get sample IDs with optional filtering.

        Args:
            tumor_only: If True, return only tumor samples
            normal_only: If True, return only normal samples

        Returns:
            List of sample IDs
        """
        # Get column names from phospho data as sample identifiers
        all_samples = list(self.phospho_deduplicated.columns)

        if tumor_only:
            return [s for s in all_samples if '.N' not in s]
        elif normal_only:
            return [s for s in all_samples if '.N' in s]
        else:
            return all_samples


class CancerManager:
    """
    Manages Cancer instances with lazy loading and caching.

    Cancer objects are only instantiated when first requested and then cached
    for subsequent access. This avoids loading all CPTAC datasets at startup.
    """
    _instances: Dict[str, Cancer] = {}

    @classmethod
    def get(cls, cancer_name: str) -> Cancer:
        """
        Get or create a Cancer instance.

        Args:
            cancer_name: Name of the cancer type (e.g., 'brca', 'luad')

        Returns:
            Cancer instance for the specified cancer type

        Raises:
            ValueError: If the cancer type is not supported
        """
        if cancer_name not in cls._instances:
            if cancer_name not in CANCER_REGISTRY:
                raise ValueError(f"Unknown cancer type: {cancer_name}. Supported types: {list(CANCER_REGISTRY.keys())}")

            logging.info(f"Instantiating CPTAC {cancer_name} dataset...")
            cptac_class = CANCER_REGISTRY[cancer_name]
            cptac_obj = cptac_class()
            cls._instances[cancer_name] = Cancer(name=cancer_name, cptac_obj=cptac_obj)
            logging.info(f"CPTAC {cancer_name} dataset ready")

        return cls._instances[cancer_name]

    @classmethod
    def available_cancers(cls) -> List[str]:
        """Return list of available cancer types."""
        return list(CANCER_REGISTRY.keys())

    @classmethod
    def clear_cache(cls, cancer_name: Optional[str] = None) -> None:
        """
        Clear cached Cancer instances.

        Args:
            cancer_name: If provided, clear only that cancer's instance.
                        If None, clear all cached instances.
        """
        if cancer_name:
            if cancer_name in cls._instances:
                cls._instances[cancer_name].clear_cache()
                del cls._instances[cancer_name]
                logging.info(f"Cleared cache for {cancer_name}")
        else:
            for name in list(cls._instances.keys()):
                cls._instances[name].clear_cache()
            cls._instances.clear()
            logging.info("Cleared all cancer caches")

    @classmethod
    def loaded_cancers(cls) -> List[str]:
        """Return list of currently loaded cancer types."""
        return list(cls._instances.keys())
