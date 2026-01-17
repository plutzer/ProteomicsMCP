"""
Test fixtures for ProteomicsMCP.

This module provides mock data generators for CPTAC and PSP data.
"""
from .cptac_mock_data import (
    create_phospho_dataframe,
    create_proteomics_dataframe,
    create_proteomics_simple_index,
)
from .psp_mock_data import (
    create_kinase_substrate_df,
    create_regulatory_sites_df,
    create_disease_sites_df,
    create_phospho_sites_df,
)
