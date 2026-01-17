"""
Shared pytest fixtures for ProteomicsMCP tests.

This module provides fixtures that can be used across all test modules.
"""
import pytest
import pandas as pd
import numpy as np
from unittest.mock import MagicMock, patch
from typing import Dict, List, Any
import os
import sys

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from tests.fixtures.cptac_mock_data import (
    create_phospho_dataframe,
    create_proteomics_dataframe,
    create_proteomics_simple_index,
    create_phospho_with_duplicates,
)
from tests.fixtures.psp_mock_data import (
    create_kinase_substrate_df,
    create_regulatory_sites_df,
    create_disease_sites_df,
    create_phospho_sites_df,
)


# ============================================================================
# CPTAC Fixtures
# ============================================================================

@pytest.fixture
def sample_phospho_multiindex():
    """
    Create a sample phosphoproteomics DataFrame with 4-level MultiIndex.

    Structure mimics CPTAC phospho data:
    - Index: (Name, Site, Peptide, Database_ID)
    - Columns: Tumor samples + Normal samples (.N suffix)
    """
    return create_phospho_dataframe()


@pytest.fixture
def sample_proteomics_multiindex():
    """
    Create a sample proteomics DataFrame with MultiIndex.

    Structure mimics CPTAC proteomics data with MultiIndex:
    - Index: (Name, Database_ID)
    - Columns: Tumor samples + Normal samples (.N suffix)
    """
    return create_proteomics_dataframe()


@pytest.fixture
def sample_proteomics_simple_index():
    """
    Create a sample proteomics DataFrame with simple gene name index.

    Some CPTAC datasets use simple index instead of MultiIndex.
    """
    return create_proteomics_simple_index()


@pytest.fixture
def sample_phospho_with_duplicates():
    """Create a phospho DataFrame with duplicate index entries for deduplication tests."""
    return create_phospho_with_duplicates()


@pytest.fixture
def mock_cptac_cancer_obj(sample_phospho_multiindex, sample_proteomics_multiindex):
    """
    Create a mock cptac cancer object mimicking cptac.Brca() behavior.

    This mock returns transposed DataFrames from get_phosphoproteomics and get_proteomics
    to match the expected .T operation in the backend.
    """
    mock_obj = MagicMock()

    # Return transposed data (backend does .T on the result)
    mock_obj.get_phosphoproteomics.return_value = sample_phospho_multiindex.T
    mock_obj.get_proteomics.return_value = sample_proteomics_multiindex.T

    # Mock clinical data
    clinical_df = pd.DataFrame({
        'age': [55, 62, 48, 70, 45],
        'sex': ['F', 'M', 'F', 'F', 'M'],
        'stage': ['I', 'II', 'III', 'II', 'I'],
    }, index=['S001', 'S002', 'S003', 'S004', 'S005'])
    mock_obj.get_clinical.return_value = clinical_df

    return mock_obj


@pytest.fixture
def mock_cancer_instance(mock_cptac_cancer_obj, sample_phospho_multiindex, sample_proteomics_multiindex):
    """
    Create a mock Cancer instance with pre-populated cached data.

    This fixture bypasses lazy loading by directly setting cached properties.
    """
    from cptac_backend import Cancer

    cancer = Cancer(name='test_cancer', cptac_obj=mock_cptac_cancer_obj)

    # Pre-populate caches to avoid cptac calls
    cancer._phospho_deduplicated = sample_phospho_multiindex
    cancer._phospho_normalized = sample_phospho_multiindex.copy()  # Use same for simplicity
    cancer._proteomics_deduplicated = sample_proteomics_multiindex

    return cancer


@pytest.fixture
def mock_cancer_manager(mock_cancer_instance):
    """
    Create a patched CancerManager that returns mock Cancer instances.

    Use this fixture to test MCP tools without network calls.
    """
    with patch('cptac_mcp.CancerManager') as mock_manager:
        mock_manager.get.return_value = mock_cancer_instance
        mock_manager.available_cancers.return_value = ['brca', 'coad', 'luad']
        yield mock_manager


# ============================================================================
# PSP Fixtures
# ============================================================================

@pytest.fixture
def sample_kinase_substrate_df():
    """Create a sample kinase-substrate DataFrame matching PSP format."""
    return create_kinase_substrate_df()


@pytest.fixture
def sample_regulatory_sites_df():
    """Create a sample regulatory sites DataFrame matching PSP format."""
    return create_regulatory_sites_df()


@pytest.fixture
def sample_disease_sites_df():
    """Create a sample disease-associated sites DataFrame matching PSP format."""
    return create_disease_sites_df()


@pytest.fixture
def sample_phospho_sites_df():
    """Create a sample phosphorylation sites DataFrame matching PSP format."""
    return create_phospho_sites_df()


@pytest.fixture
def mock_psp_data(
    sample_kinase_substrate_df,
    sample_regulatory_sites_df,
    sample_disease_sites_df,
    sample_phospho_sites_df
):
    """
    Patch PSP global variables with mock data.

    Use this fixture to test PSP MCP tools without loading real data files.
    """
    with patch.multiple(
        'psp_proteomics',
        kinase_substrate_data=sample_kinase_substrate_df,
        regulatory_data=sample_regulatory_sites_df,
        disease_data=sample_disease_sites_df,
        phospho_data=sample_phospho_sites_df,
    ):
        yield {
            'kinase_substrate': sample_kinase_substrate_df,
            'regulatory': sample_regulatory_sites_df,
            'disease': sample_disease_sites_df,
            'phospho': sample_phospho_sites_df,
        }


# ============================================================================
# CSV Parsing Helpers
# ============================================================================

@pytest.fixture
def csv_parser():
    """
    Provide a helper function for parsing CSV strings from MCP tool results.

    Returns:
    --------
    callable
        Function that parses a CSV string into a list of dictionaries
    """
    def parse_csv(csv_string: str) -> List[Dict[str, str]]:
        """Parse CSV string into list of dicts."""
        if not csv_string or not csv_string.strip():
            return []

        lines = csv_string.strip().split('\n')
        if len(lines) < 2:
            return []

        headers = lines[0].split(',')
        rows = []

        for line in lines[1:]:
            values = line.split(',')
            # Pad with empty strings if needed
            while len(values) < len(headers):
                values.append('')
            rows.append(dict(zip(headers, values)))

        return rows

    return parse_csv


@pytest.fixture
def assert_csv_structure():
    """
    Provide a helper function for asserting CSV structure from MCP tool results.

    Returns:
    --------
    callable
        Function that asserts CSV has expected columns
    """
    def _assert_csv_structure(csv_string: str, expected_columns: List[str]):
        """Assert CSV has expected column headers."""
        if not csv_string or not csv_string.strip():
            pytest.fail("CSV string is empty")

        lines = csv_string.strip().split('\n')
        headers = lines[0].split(',')

        for col in expected_columns:
            assert col in headers, f"Expected column '{col}' not found in headers: {headers}"

    return _assert_csv_structure


# ============================================================================
# Test Data Helpers
# ============================================================================

@pytest.fixture
def standard_test_genes():
    """Return standard list of genes used in tests."""
    return ['AKT1', 'TP53', 'MTOR', 'GSK3B', 'MAPK1', 'EGFR']


@pytest.fixture
def standard_test_sites():
    """Return standard list of phosphosites used in tests."""
    return ['AKT1_S473', 'AKT1_T308', 'TP53_S15', 'MTOR_S2448', 'GSK3B_S9']


@pytest.fixture
def standard_tumor_samples():
    """Return standard tumor sample IDs."""
    return ['S001', 'S002', 'S003', 'S004', 'S005']


@pytest.fixture
def standard_normal_samples():
    """Return standard normal sample IDs."""
    return ['S001', 'S002', 'S003']


# ============================================================================
# Real Data Fixtures (for e2e tests)
# ============================================================================

@pytest.fixture
def psp_data_dir():
    """Return path to PSP data directory."""
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(base_dir, 'datasets', 'PSP_2025')


@pytest.fixture
def psp_data_available(psp_data_dir):
    """Check if PSP data files are available."""
    required_files = [
        'Phosphorylation_site_dataset.txt',
        'Kinase_Substrate_Dataset.txt',
        'Regulatory_sites.txt',
        'Disease-associated_sites.txt',
    ]

    for filename in required_files:
        filepath = os.path.join(psp_data_dir, filename)
        if not os.path.exists(filepath):
            return False

    return True


# ============================================================================
# Parametrization Helpers
# ============================================================================

def pytest_generate_tests(metafunc):
    """Generate parametrized test cases."""
    # Query parsing test cases
    if "query_parse_case" in metafunc.fixturenames:
        cases = [
            ("AKT1", ("AKT1", None)),
            ("AKT1_S473", ("AKT1", "S473")),
            ("AKT1_S473_S474", ("AKT1", "S473_S474")),
            ("MTOR_S2448", ("MTOR", "S2448")),
            ("GSK3B", ("GSK3B", None)),
        ]
        metafunc.parametrize("query_parse_case", cases, ids=[c[0] for c in cases])

    # Cancer type test cases
    if "cancer_type" in metafunc.fixturenames:
        cancers = ['brca', 'coad', 'hnscc', 'luad', 'ovarian', 'ccrcc', 'gbm', 'lscc', 'pdac']
        metafunc.parametrize("cancer_type", cancers)


# ============================================================================
# Logging Configuration
# ============================================================================

@pytest.fixture(autouse=True)
def configure_logging():
    """Configure logging for tests."""
    import logging
    logging.basicConfig(level=logging.WARNING)
    # Suppress noisy loggers during tests
    logging.getLogger('cptac').setLevel(logging.ERROR)
    logging.getLogger('cptac_backend').setLevel(logging.WARNING)
    logging.getLogger('cptac_mcp').setLevel(logging.WARNING)
    logging.getLogger('psp_proteomics').setLevel(logging.WARNING)
