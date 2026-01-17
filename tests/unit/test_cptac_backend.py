"""
Unit tests for cptac_backend.py

Tests the Cancer dataclass and CancerManager class.
"""
import pytest
import pandas as pd
import numpy as np
from unittest.mock import MagicMock, patch, PropertyMock
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from cptac_backend import Cancer, CancerManager, CANCER_REGISTRY


class TestCancerDataclass:
    """Tests for the Cancer dataclass."""

    def test_cancer_initialization(self, mock_cptac_cancer_obj):
        """Verify Cancer fields initialize correctly."""
        cancer = Cancer(name='test_brca', cptac_obj=mock_cptac_cancer_obj)

        assert cancer.name == 'test_brca'
        assert cancer.cptac_obj is mock_cptac_cancer_obj
        assert cancer._phospho_deduplicated is None
        assert cancer._phospho_normalized is None
        assert cancer._proteomics_deduplicated is None
        assert cancer._clinical_data is None

    def test_phospho_deduplicated_lazy_loading(self, mock_cptac_cancer_obj):
        """Verify lazy load triggers on first access."""
        cancer = Cancer(name='test', cptac_obj=mock_cptac_cancer_obj)

        # Cache should be None initially
        assert cancer._phospho_deduplicated is None

        # Access property to trigger lazy loading
        result = cancer.phospho_deduplicated

        # Should have called get_phosphoproteomics
        mock_cptac_cancer_obj.get_phosphoproteomics.assert_called_once_with('bcm')

        # Result should be a DataFrame
        assert isinstance(result, pd.DataFrame)
        assert result.index.nlevels == 4  # MultiIndex with 4 levels

    def test_phospho_deduplicated_caching(self, mock_cptac_cancer_obj):
        """Verify second access returns cached data."""
        cancer = Cancer(name='test', cptac_obj=mock_cptac_cancer_obj)

        # First access
        result1 = cancer.phospho_deduplicated

        # Second access
        result2 = cancer.phospho_deduplicated

        # Should only call get_phosphoproteomics once
        assert mock_cptac_cancer_obj.get_phosphoproteomics.call_count == 1

        # Both results should be the same object
        assert result1 is result2

    def test_proteomics_deduplicated_lazy_loading(self, mock_cptac_cancer_obj):
        """Verify proteomics lazy loading."""
        cancer = Cancer(name='test', cptac_obj=mock_cptac_cancer_obj)

        assert cancer._proteomics_deduplicated is None

        result = cancer.proteomics_deduplicated

        mock_cptac_cancer_obj.get_proteomics.assert_called_once_with('bcm')
        assert isinstance(result, pd.DataFrame)

    def test_phospho_normalized_uses_cached_phospho(self, mock_cptac_cancer_obj):
        """Verify normalized phospho uses cached deduplicated data."""
        cancer = Cancer(name='test', cptac_obj=mock_cptac_cancer_obj)

        # Access normalized which should use cached deduplicated
        result = cancer.phospho_normalized

        # get_phosphoproteomics called for deduplicated data
        # get_proteomics called for whole cell normalization
        mock_cptac_cancer_obj.get_phosphoproteomics.assert_called()
        mock_cptac_cancer_obj.get_proteomics.assert_called()

        assert isinstance(result, pd.DataFrame)

    def test_clear_cache(self, mock_cptac_cancer_obj, sample_phospho_multiindex):
        """Verify all cached fields reset to None."""
        cancer = Cancer(name='test', cptac_obj=mock_cptac_cancer_obj)

        # Pre-populate caches
        cancer._phospho_deduplicated = sample_phospho_multiindex
        cancer._phospho_normalized = sample_phospho_multiindex.copy()
        cancer._proteomics_deduplicated = sample_phospho_multiindex.copy()
        cancer._clinical_data = pd.DataFrame({'a': [1, 2, 3]})

        # Clear cache
        cancer.clear_cache()

        # All caches should be None
        assert cancer._phospho_deduplicated is None
        assert cancer._phospho_normalized is None
        assert cancer._proteomics_deduplicated is None
        assert cancer._clinical_data is None

    def test_get_sample_ids_all(self, mock_cancer_instance):
        """Verify get_sample_ids returns all samples by default."""
        samples = mock_cancer_instance.get_sample_ids()

        # Should include both tumor and normal samples
        assert 'S001' in samples
        assert 'S001.N' in samples
        assert len(samples) == 8  # 5 tumor + 3 normal

    def test_get_sample_ids_tumor_only(self, mock_cancer_instance):
        """Verify get_sample_ids filters to tumor samples only."""
        samples = mock_cancer_instance.get_sample_ids(tumor_only=True)

        # Should only include tumor samples (no .N suffix)
        assert 'S001' in samples
        assert 'S001.N' not in samples
        assert len(samples) == 5

        # Verify no normal samples
        for s in samples:
            assert '.N' not in s

    def test_get_sample_ids_normal_only(self, mock_cancer_instance):
        """Verify get_sample_ids filters to normal samples only."""
        samples = mock_cancer_instance.get_sample_ids(normal_only=True)

        # Should only include normal samples (with .N suffix)
        assert 'S001' not in samples
        assert 'S001.N' in samples
        assert len(samples) == 3

        # Verify all samples have .N suffix
        for s in samples:
            assert '.N' in s


class TestCancerManager:
    """Tests for the CancerManager class."""

    def test_available_cancers(self):
        """Verify available_cancers returns list from CANCER_REGISTRY."""
        available = CancerManager.available_cancers()

        assert isinstance(available, list)
        assert 'brca' in available
        assert 'coad' in available
        assert 'luad' in available
        assert len(available) == len(CANCER_REGISTRY)

    def test_get_invalid_cancer_raises(self):
        """Verify ValueError is raised for unknown cancer type."""
        with pytest.raises(ValueError) as exc_info:
            CancerManager.get('invalid_cancer_type')

        assert 'Unknown cancer type' in str(exc_info.value)
        assert 'invalid_cancer_type' in str(exc_info.value)

    def test_loaded_cancers_initially_empty(self):
        """Verify loaded_cancers returns empty list initially after clearing."""
        # Clear any existing cache first
        CancerManager.clear_cache()

        loaded = CancerManager.loaded_cancers()

        assert isinstance(loaded, list)
        # After clearing, should be empty
        assert len(loaded) == 0

    @patch('cptac_backend.CANCER_REGISTRY', {'test': MagicMock})
    def test_get_caches_cancer_instance(self):
        """Verify get() caches cancer instances."""
        # Clear cache first
        CancerManager._instances.clear()

        with patch('cptac_backend.CANCER_REGISTRY', {'test': MagicMock}):
            # First access creates instance
            cancer1 = CancerManager.get('test')

            # Second access returns cached instance
            cancer2 = CancerManager.get('test')

            assert cancer1 is cancer2
            assert 'test' in CancerManager.loaded_cancers()

        # Cleanup
        CancerManager._instances.clear()

    def test_clear_cache_single_cancer(self):
        """Verify clear_cache can clear a specific cancer."""
        # Setup: add mock cancer
        mock_cancer = MagicMock()
        CancerManager._instances['test_cancer'] = mock_cancer

        # Clear specific cancer
        CancerManager.clear_cache('test_cancer')

        assert 'test_cancer' not in CancerManager._instances
        mock_cancer.clear_cache.assert_called_once()

    def test_clear_cache_all_cancers(self):
        """Verify clear_cache clears all cancers when no argument given."""
        # Setup: add mock cancers
        mock1 = MagicMock()
        mock2 = MagicMock()
        CancerManager._instances['cancer1'] = mock1
        CancerManager._instances['cancer2'] = mock2

        # Clear all
        CancerManager.clear_cache()

        assert len(CancerManager._instances) == 0
        mock1.clear_cache.assert_called_once()
        mock2.clear_cache.assert_called_once()


class TestDataDeduplication:
    """Tests for data deduplication logic."""

    def test_groupby_mean_deduplication(self, sample_phospho_with_duplicates):
        """Verify duplicate indices are averaged correctly."""
        # The sample_phospho_with_duplicates has duplicate entries
        df = sample_phospho_with_duplicates

        # Apply deduplication (same logic as in Cancer._load_deduplicated_phospho)
        deduplicated = df.groupby(level=[0, 1, 2, 3]).agg('mean')

        # Should have fewer rows after deduplication
        assert len(deduplicated) < len(df)

        # Check specific duplicate case
        # ('AKT1', 'S473', 'PEPTIDE_1', 'NP_001') had 2 rows with same index
        # They should now be averaged
        assert ('AKT1', 'S473', 'PEPTIDE_1', 'NP_001') in deduplicated.index

    def test_inf_values_replaced_with_nan(self, mock_cptac_cancer_obj):
        """Verify infinite values are replaced with NaN during deduplication."""
        # Create DataFrame with inf values
        df = pd.DataFrame(
            [[1.0, np.inf, -np.inf, 2.0]],
            index=pd.MultiIndex.from_tuples([('AKT1', 'S473', 'PEP', 'DB')],
                                             names=['Name', 'Site', 'Peptide', 'Database_ID']),
            columns=['S001', 'S002', 'S003', 'S004']
        )

        mock_cptac_cancer_obj.get_phosphoproteomics.return_value = df.T
        cancer = Cancer(name='test', cptac_obj=mock_cptac_cancer_obj)

        result = cancer.phospho_deduplicated

        # Inf values should be NaN
        assert pd.isna(result.loc[('AKT1', 'S473', 'PEP', 'DB'), 'S002'])
        assert pd.isna(result.loc[('AKT1', 'S473', 'PEP', 'DB'), 'S003'])
        # Valid values should be preserved
        assert result.loc[('AKT1', 'S473', 'PEP', 'DB'), 'S001'] == 1.0
        assert result.loc[('AKT1', 'S473', 'PEP', 'DB'), 'S004'] == 2.0


class TestCancerRegistry:
    """Tests for the CANCER_REGISTRY constant."""

    def test_registry_contains_expected_cancers(self):
        """Verify registry has all expected cancer types."""
        expected_cancers = ['brca', 'coad', 'hnscc', 'luad', 'ovarian', 'ccrcc', 'gbm', 'lscc', 'pdac']

        for cancer in expected_cancers:
            assert cancer in CANCER_REGISTRY

    def test_registry_values_are_classes(self):
        """Verify registry values are cptac classes (not instances)."""
        for name, cls in CANCER_REGISTRY.items():
            assert callable(cls), f"Registry value for {name} should be a class/callable"


class TestNormalSampleIdentification:
    """Tests for identifying normal vs tumor samples."""

    def test_normal_sample_suffix_detection(self, sample_phospho_multiindex):
        """Verify normal samples are correctly identified by .N suffix."""
        columns = sample_phospho_multiindex.columns.tolist()

        normal_cols = [col for col in columns if '.N' in col]
        tumor_cols = [col for col in columns if '.N' not in col]

        assert len(normal_cols) == 3  # S001.N, S002.N, S003.N
        assert len(tumor_cols) == 5   # S001, S002, S003, S004, S005

        # Verify normal samples have corresponding tumor samples
        for normal_col in normal_cols:
            tumor_col = normal_col.replace('.N', '')
            assert tumor_col in tumor_cols

    def test_paired_sample_extraction(self, sample_phospho_multiindex):
        """Verify tumor-normal pairs can be extracted correctly."""
        columns = sample_phospho_multiindex.columns.tolist()
        normal_cols = [col for col in columns if '.N' in col]

        pairs = []
        for normal_col in normal_cols:
            tumor_col = normal_col.replace('.N', '')
            if tumor_col in columns:
                pairs.append((tumor_col, normal_col))

        assert len(pairs) == 3
        assert ('S001', 'S001.N') in pairs
        assert ('S002', 'S002.N') in pairs
        assert ('S003', 'S003.N') in pairs
