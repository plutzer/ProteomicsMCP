"""
CPTAC Mock Data Generators

This module provides functions to generate mock CPTAC-style DataFrames
for testing purposes.
"""
import pandas as pd
import numpy as np
from typing import List, Optional


def create_phospho_dataframe(
    genes: Optional[List[str]] = None,
    sites_per_gene: int = 3,
    tumor_samples: Optional[List[str]] = None,
    normal_samples: Optional[List[str]] = None,
    include_nan: bool = True,
    seed: int = 42,
) -> pd.DataFrame:
    """
    Create a mock phosphoproteomics DataFrame with 4-level MultiIndex.

    The DataFrame mimics CPTAC phospho data structure:
    - MultiIndex: (Name/Gene, Site, Peptide, Database_ID)
    - Columns: Sample IDs (tumor samples + normal samples with .N suffix)
    - Values: Log2 abundance values

    Parameters:
    -----------
    genes : list of str, optional
        Gene names to include. Default: ['AKT1', 'TP53', 'MTOR', 'GSK3B', 'MAPK1', 'EGFR']
    sites_per_gene : int
        Number of phosphosites per gene. Default: 3
    tumor_samples : list of str, optional
        Tumor sample IDs. Default: ['S001', 'S002', 'S003', 'S004', 'S005']
    normal_samples : list of str, optional
        Normal sample IDs (will have .N suffix). Default: ['S001', 'S002', 'S003']
    include_nan : bool
        Whether to include NaN values in the data. Default: True
    seed : int
        Random seed for reproducibility. Default: 42

    Returns:
    --------
    pd.DataFrame
        Mock phosphoproteomics DataFrame with 4-level MultiIndex
    """
    np.random.seed(seed)

    if genes is None:
        genes = ['AKT1', 'TP53', 'MTOR', 'GSK3B', 'MAPK1', 'EGFR']

    if tumor_samples is None:
        tumor_samples = ['S001', 'S002', 'S003', 'S004', 'S005']

    if normal_samples is None:
        normal_samples = ['S001', 'S002', 'S003']

    # Standard phosphorylation sites for common genes
    site_map = {
        'AKT1': ['S473', 'T308', 'S124'],
        'TP53': ['S15', 'S392', 'S46'],
        'MTOR': ['S2448', 'S2481', 'T2446'],
        'GSK3B': ['S9', 'Y216', 'S389'],
        'MAPK1': ['T185', 'Y187', 'T202'],
        'EGFR': ['Y1068', 'Y1173', 'Y992'],
    }

    # Build MultiIndex tuples
    index_tuples = []
    for gene in genes:
        sites = site_map.get(gene, [f'S{i}' for i in range(100, 100 + sites_per_gene)])[:sites_per_gene]
        for i, site in enumerate(sites):
            peptide = f"PEPTIDE_{gene}_{site}"
            db_id = f"NP_{hash(gene) % 10000:05d}_{site}"
            index_tuples.append((gene, site, peptide, db_id))

    # Create MultiIndex
    index = pd.MultiIndex.from_tuples(
        index_tuples,
        names=['Name', 'Site', 'Peptide', 'Database_ID']
    )

    # Create column names
    columns = tumor_samples + [f"{s}.N" for s in normal_samples]

    # Generate random data
    n_rows = len(index_tuples)
    n_cols = len(columns)

    # Generate values with some structure:
    # - Tumor samples have slightly higher values on average
    # - Normal samples cluster together
    data = np.random.randn(n_rows, n_cols) * 2

    # Add tumor vs normal difference
    n_tumor = len(tumor_samples)
    data[:, :n_tumor] += 0.5  # Tumor samples slightly higher

    # Introduce NaN values if requested
    if include_nan:
        nan_mask = np.random.random((n_rows, n_cols)) < 0.1
        data[nan_mask] = np.nan

    df = pd.DataFrame(data, index=index, columns=columns)

    return df


def create_proteomics_dataframe(
    genes: Optional[List[str]] = None,
    tumor_samples: Optional[List[str]] = None,
    normal_samples: Optional[List[str]] = None,
    include_nan: bool = True,
    seed: int = 42,
) -> pd.DataFrame:
    """
    Create a mock proteomics DataFrame with MultiIndex.

    The DataFrame mimics CPTAC proteomics data structure:
    - MultiIndex: (Name/Gene, Database_ID)
    - Columns: Sample IDs (tumor samples + normal samples with .N suffix)
    - Values: Log2 abundance values

    Parameters:
    -----------
    genes : list of str, optional
        Gene names to include. Default: ['AKT1', 'TP53', 'MTOR', 'GSK3B', 'MAPK1', 'EGFR']
    tumor_samples : list of str, optional
        Tumor sample IDs. Default: ['S001', 'S002', 'S003', 'S004', 'S005']
    normal_samples : list of str, optional
        Normal sample IDs (will have .N suffix). Default: ['S001', 'S002', 'S003']
    include_nan : bool
        Whether to include NaN values in the data. Default: True
    seed : int
        Random seed for reproducibility. Default: 42

    Returns:
    --------
    pd.DataFrame
        Mock proteomics DataFrame with MultiIndex
    """
    np.random.seed(seed)

    if genes is None:
        genes = ['AKT1', 'TP53', 'MTOR', 'GSK3B', 'MAPK1', 'EGFR']

    if tumor_samples is None:
        tumor_samples = ['S001', 'S002', 'S003', 'S004', 'S005']

    if normal_samples is None:
        normal_samples = ['S001', 'S002', 'S003']

    # Build MultiIndex tuples
    index_tuples = []
    for gene in genes:
        db_id = f"NP_{hash(gene) % 10000:05d}"
        index_tuples.append((gene, db_id))

    # Create MultiIndex
    index = pd.MultiIndex.from_tuples(
        index_tuples,
        names=['Name', 'Database_ID']
    )

    # Create column names
    columns = tumor_samples + [f"{s}.N" for s in normal_samples]

    # Generate random data
    n_rows = len(index_tuples)
    n_cols = len(columns)

    data = np.random.randn(n_rows, n_cols) * 1.5

    # Add tumor vs normal difference
    n_tumor = len(tumor_samples)
    data[:, :n_tumor] += 0.3  # Tumor samples slightly higher

    # Introduce NaN values if requested
    if include_nan:
        nan_mask = np.random.random((n_rows, n_cols)) < 0.05
        data[nan_mask] = np.nan

    df = pd.DataFrame(data, index=index, columns=columns)

    return df


def create_proteomics_simple_index(
    genes: Optional[List[str]] = None,
    tumor_samples: Optional[List[str]] = None,
    normal_samples: Optional[List[str]] = None,
    include_nan: bool = True,
    seed: int = 42,
) -> pd.DataFrame:
    """
    Create a mock proteomics DataFrame with simple gene name index.

    Some CPTAC datasets have simple index (just gene names) instead of MultiIndex.

    Parameters:
    -----------
    genes : list of str, optional
        Gene names to include. Default: ['AKT1', 'TP53', 'MTOR', 'GSK3B', 'MAPK1', 'EGFR']
    tumor_samples : list of str, optional
        Tumor sample IDs. Default: ['S001', 'S002', 'S003', 'S004', 'S005']
    normal_samples : list of str, optional
        Normal sample IDs (will have .N suffix). Default: ['S001', 'S002', 'S003']
    include_nan : bool
        Whether to include NaN values in the data. Default: True
    seed : int
        Random seed for reproducibility. Default: 42

    Returns:
    --------
    pd.DataFrame
        Mock proteomics DataFrame with simple index
    """
    np.random.seed(seed)

    if genes is None:
        genes = ['AKT1', 'TP53', 'MTOR', 'GSK3B', 'MAPK1', 'EGFR']

    if tumor_samples is None:
        tumor_samples = ['S001', 'S002', 'S003', 'S004', 'S005']

    if normal_samples is None:
        normal_samples = ['S001', 'S002', 'S003']

    # Create column names
    columns = tumor_samples + [f"{s}.N" for s in normal_samples]

    # Generate random data
    n_rows = len(genes)
    n_cols = len(columns)

    data = np.random.randn(n_rows, n_cols) * 1.5

    # Add tumor vs normal difference
    n_tumor = len(tumor_samples)
    data[:, :n_tumor] += 0.3

    # Introduce NaN values if requested
    if include_nan:
        nan_mask = np.random.random((n_rows, n_cols)) < 0.05
        data[nan_mask] = np.nan

    df = pd.DataFrame(data, index=genes, columns=columns)
    df.index.name = 'Name'

    return df


def create_phospho_with_duplicates(
    genes: Optional[List[str]] = None,
    seed: int = 42,
) -> pd.DataFrame:
    """
    Create a phospho DataFrame with duplicate index entries for testing deduplication.

    Parameters:
    -----------
    genes : list of str, optional
        Gene names. Default: ['AKT1', 'TP53']
    seed : int
        Random seed. Default: 42

    Returns:
    --------
    pd.DataFrame
        DataFrame with duplicate MultiIndex entries
    """
    np.random.seed(seed)

    if genes is None:
        genes = ['AKT1', 'TP53']

    # Create duplicated index entries
    index_tuples = [
        ('AKT1', 'S473', 'PEPTIDE_1', 'NP_001'),
        ('AKT1', 'S473', 'PEPTIDE_1', 'NP_001'),  # Duplicate
        ('AKT1', 'S473', 'PEPTIDE_2', 'NP_001'),  # Different peptide
        ('TP53', 'S15', 'PEPTIDE_3', 'NP_002'),
        ('TP53', 'S15', 'PEPTIDE_3', 'NP_002'),  # Duplicate
    ]

    index = pd.MultiIndex.from_tuples(
        index_tuples,
        names=['Name', 'Site', 'Peptide', 'Database_ID']
    )

    columns = ['S001', 'S002', 'S001.N', 'S002.N']
    data = np.random.randn(len(index_tuples), len(columns))

    return pd.DataFrame(data, index=index, columns=columns)
