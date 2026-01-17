"""
Integration tests for CPTAC MCP tools.

Tests the MCP tool functions with mocked CancerManager to avoid network calls.
"""
import pytest
import pandas as pd
import numpy as np
from unittest.mock import patch, MagicMock
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from cptac_backend import Cancer


class TestGetCancerTypes:
    """Tests for get_cancer_types MCP tool."""

    def test_returns_list(self):
        """Test that get_cancer_types returns a list."""
        from cptac_mcp import get_cancer_types

        result = get_cancer_types()

        assert isinstance(result, list)

    def test_contains_expected_cancers(self):
        """Test that result contains expected cancer types."""
        from cptac_mcp import get_cancer_types

        result = get_cancer_types()

        expected = ['brca', 'coad', 'hnscc', 'luad', 'ovarian', 'ccrcc', 'gbm', 'lscc', 'pdac']
        for cancer in expected:
            assert cancer in result


@pytest.mark.integration
class TestPhosphoTumorVsNormal:
    """Tests for phospho_tumor_vs_normal MCP tool."""

    def test_single_site_query(self, mock_cancer_instance, csv_parser):
        """Test query for a single phosphosite."""
        with patch('cptac_mcp.CancerManager') as mock_manager:
            mock_manager.get.return_value = mock_cancer_instance
            mock_manager.available_cancers.return_value = ['brca']

            from cptac_mcp import phospho_tumor_vs_normal

            result = phospho_tumor_vs_normal(
                cancer='brca',
                query='AKT1_S473',
                normalized='false'
            )

            assert 'data' in result
            assert 'AKT1' in result['data']

            # Parse CSV and check structure
            csv_data = result['data']['AKT1']
            rows = csv_parser(csv_data)

            assert len(rows) > 0
            # Check expected columns exist
            assert 'site' in rows[0]
            assert 'log2_fold_change' in rows[0]
            assert 'p_value' in rows[0]

    def test_gene_query_all_sites(self, mock_cancer_instance, csv_parser):
        """Test query for gene name returns all sites."""
        with patch('cptac_mcp.CancerManager') as mock_manager:
            mock_manager.get.return_value = mock_cancer_instance
            mock_manager.available_cancers.return_value = ['brca']

            from cptac_mcp import phospho_tumor_vs_normal

            result = phospho_tumor_vs_normal(
                cancer='brca',
                query='AKT1',
                normalized='false'
            )

            assert 'data' in result
            assert 'AKT1' in result['data']

            # Should return multiple sites for AKT1
            csv_data = result['data']['AKT1']
            rows = csv_parser(csv_data)

            assert len(rows) >= 3  # AKT1 has S473, T308, S124

    def test_multiple_genes_query(self, mock_cancer_instance):
        """Test query for multiple genes."""
        with patch('cptac_mcp.CancerManager') as mock_manager:
            mock_manager.get.return_value = mock_cancer_instance
            mock_manager.available_cancers.return_value = ['brca']

            from cptac_mcp import phospho_tumor_vs_normal

            result = phospho_tumor_vs_normal(
                cancer='brca',
                query='AKT1,TP53',
                normalized='false'
            )

            assert 'data' in result
            assert 'AKT1' in result['data']
            assert 'TP53' in result['data']

    def test_invalid_normalized_param(self, mock_cancer_instance):
        """Test that invalid normalized parameter raises ValueError."""
        with patch('cptac_mcp.CancerManager') as mock_manager:
            mock_manager.get.return_value = mock_cancer_instance

            from cptac_mcp import phospho_tumor_vs_normal

            with pytest.raises(ValueError) as exc_info:
                phospho_tumor_vs_normal(
                    cancer='brca',
                    query='AKT1_S473',
                    normalized='invalid'
                )

            assert 'normalized' in str(exc_info.value).lower()

    def test_no_matches_returns_error(self, mock_cancer_instance):
        """Test that query with no matches returns error dict."""
        with patch('cptac_mcp.CancerManager') as mock_manager:
            mock_manager.get.return_value = mock_cancer_instance
            mock_manager.available_cancers.return_value = ['brca']

            from cptac_mcp import phospho_tumor_vs_normal

            result = phospho_tumor_vs_normal(
                cancer='brca',
                query='NONEXISTENT_GENE',
                normalized='false'
            )

            assert 'error' in result

    def test_fdr_correction_applied(self, mock_cancer_instance, csv_parser):
        """Test that p_value_adjusted is present (FDR correction)."""
        with patch('cptac_mcp.CancerManager') as mock_manager:
            mock_manager.get.return_value = mock_cancer_instance
            mock_manager.available_cancers.return_value = ['brca']

            from cptac_mcp import phospho_tumor_vs_normal

            result = phospho_tumor_vs_normal(
                cancer='brca',
                query='AKT1',
                normalized='false'
            )

            csv_data = result['data']['AKT1']
            rows = csv_parser(csv_data)

            # Check p_value_adjusted column exists
            assert 'p_value_adjusted' in rows[0]

    def test_normalized_true_uses_normalized_data(self, mock_cancer_instance):
        """Test that normalized='true' uses phospho_normalized property."""
        with patch('cptac_mcp.CancerManager') as mock_manager:
            mock_manager.get.return_value = mock_cancer_instance
            mock_manager.available_cancers.return_value = ['brca']

            from cptac_mcp import phospho_tumor_vs_normal

            result = phospho_tumor_vs_normal(
                cancer='brca',
                query='AKT1_S473',
                normalized='true'
            )

            assert 'data' in result


@pytest.mark.integration
class TestProteinTumorVsNormal:
    """Tests for protein_tumor_vs_normal MCP tool."""

    def test_single_protein_query(self, mock_cancer_instance, csv_parser):
        """Test query for a single protein."""
        with patch('cptac_mcp.CancerManager') as mock_manager:
            mock_manager.get.return_value = mock_cancer_instance
            mock_manager.available_cancers.return_value = ['brca']

            from cptac_mcp import protein_tumor_vs_normal

            result = protein_tumor_vs_normal(
                cancer='brca',
                query='AKT1'
            )

            assert 'data' in result
            rows = csv_parser(result['data'])

            assert len(rows) > 0
            assert rows[0]['gene'] == 'AKT1'

    def test_multiple_proteins_query(self, mock_cancer_instance, csv_parser):
        """Test query for multiple proteins."""
        with patch('cptac_mcp.CancerManager') as mock_manager:
            mock_manager.get.return_value = mock_cancer_instance
            mock_manager.available_cancers.return_value = ['brca']

            from cptac_mcp import protein_tumor_vs_normal

            result = protein_tumor_vs_normal(
                cancer='brca',
                query='AKT1,TP53,MTOR'
            )

            assert 'data' in result
            rows = csv_parser(result['data'])

            genes = [row['gene'] for row in rows]
            assert 'AKT1' in genes
            assert 'TP53' in genes
            assert 'MTOR' in genes

    def test_no_matches_returns_error(self, mock_cancer_instance):
        """Test query with no matches returns error."""
        with patch('cptac_mcp.CancerManager') as mock_manager:
            mock_manager.get.return_value = mock_cancer_instance
            mock_manager.available_cancers.return_value = ['brca']

            from cptac_mcp import protein_tumor_vs_normal

            result = protein_tumor_vs_normal(
                cancer='brca',
                query='NONEXISTENT_GENE'
            )

            assert 'error' in result


@pytest.mark.integration
class TestCorrelationAnalysis:
    """Tests for correlation_analysis MCP tool."""

    def test_phospho_correlation(self, mock_cancer_instance, csv_parser):
        """Test phospho correlation analysis."""
        with patch('cptac_mcp.CancerManager') as mock_manager:
            mock_manager.get.return_value = mock_cancer_instance
            mock_manager.available_cancers.return_value = ['brca']

            from cptac_mcp import correlation_analysis

            result = correlation_analysis(
                cancer='brca',
                query='AKT1_S473,AKT1_T308',
                data_type='phospho',
                normalized='false'
            )

            assert 'correlation_matrix' in result
            assert 'p_value_matrix' in result
            assert 'n_samples' in result

            # Check n_samples is positive
            assert result['n_samples'] > 0

    def test_protein_correlation(self, mock_cancer_instance):
        """Test protein correlation analysis."""
        with patch('cptac_mcp.CancerManager') as mock_manager:
            mock_manager.get.return_value = mock_cancer_instance
            mock_manager.available_cancers.return_value = ['brca']

            from cptac_mcp import correlation_analysis

            result = correlation_analysis(
                cancer='brca',
                query='AKT1,TP53',
                data_type='proteomics'
            )

            assert 'correlation_matrix' in result
            assert 'p_value_matrix' in result

    def test_mixed_correlation_both_datatype(self, mock_cancer_instance):
        """Test mixed phospho and protein correlation with data_type='both'."""
        with patch('cptac_mcp.CancerManager') as mock_manager:
            mock_manager.get.return_value = mock_cancer_instance
            mock_manager.available_cancers.return_value = ['brca']

            from cptac_mcp import correlation_analysis

            result = correlation_analysis(
                cancer='brca',
                query='AKT1_S473,AKT1_protein',
                data_type='both',
                normalized='false'
            )

            assert 'correlation_matrix' in result or 'error' in result

    def test_invalid_data_type_raises(self, mock_cancer_instance):
        """Test that invalid data_type raises ValueError."""
        with patch('cptac_mcp.CancerManager') as mock_manager:
            mock_manager.get.return_value = mock_cancer_instance

            from cptac_mcp import correlation_analysis

            with pytest.raises(ValueError) as exc_info:
                correlation_analysis(
                    cancer='brca',
                    query='AKT1_S473',
                    data_type='invalid_type'
                )

            assert 'data_type' in str(exc_info.value).lower()

    def test_correlation_matrix_structure(self, mock_cancer_instance):
        """Test correlation matrix has proper structure."""
        with patch('cptac_mcp.CancerManager') as mock_manager:
            mock_manager.get.return_value = mock_cancer_instance
            mock_manager.available_cancers.return_value = ['brca']

            from cptac_mcp import correlation_analysis

            result = correlation_analysis(
                cancer='brca',
                query='AKT1_S473,AKT1_T308,TP53_S15',
                data_type='phospho',
                normalized='false'
            )

            # Parse correlation matrix
            lines = result['correlation_matrix'].split('\n')
            header = lines[0].split(',')

            # Should have empty first cell + item labels
            assert header[0] == ''
            assert len(header) >= 2  # At least 2 items

            # Check matrix is square
            n_cols = len(header)
            n_rows = len(lines)
            assert n_cols == n_rows  # Including header

    def test_diagonal_correlation_is_one(self, mock_cancer_instance):
        """Test that diagonal of correlation matrix is 1."""
        with patch('cptac_mcp.CancerManager') as mock_manager:
            mock_manager.get.return_value = mock_cancer_instance
            mock_manager.available_cancers.return_value = ['brca']

            from cptac_mcp import correlation_analysis

            result = correlation_analysis(
                cancer='brca',
                query='AKT1_S473,AKT1_T308',
                data_type='phospho',
                normalized='false'
            )

            lines = result['correlation_matrix'].split('\n')

            # Check diagonal values (should be 1.000)
            for i, line in enumerate(lines[1:], start=1):
                values = line.split(',')
                diagonal_value = values[i]
                assert diagonal_value == '1.000', f"Diagonal at position {i} should be 1.000"


@pytest.mark.integration
class TestCSVOutputFormat:
    """Tests for verifying CSV output format of MCP tools."""

    def test_phospho_csv_columns(self, mock_cancer_instance, assert_csv_structure):
        """Verify phospho tumor vs normal CSV has expected columns."""
        with patch('cptac_mcp.CancerManager') as mock_manager:
            mock_manager.get.return_value = mock_cancer_instance
            mock_manager.available_cancers.return_value = ['brca']

            from cptac_mcp import phospho_tumor_vs_normal

            result = phospho_tumor_vs_normal(
                cancer='brca',
                query='AKT1_S473',
                normalized='false'
            )

            expected_cols = ['site', 'peptide', 'database_id', 'log2_fold_change',
                           'p_value', 'p_value_adjusted', 'n_pairs', 'mean_tumor', 'mean_normal']

            assert_csv_structure(result['data']['AKT1'], expected_cols)

    def test_protein_csv_columns(self, mock_cancer_instance, assert_csv_structure):
        """Verify protein tumor vs normal CSV has expected columns."""
        with patch('cptac_mcp.CancerManager') as mock_manager:
            mock_manager.get.return_value = mock_cancer_instance
            mock_manager.available_cancers.return_value = ['brca']

            from cptac_mcp import protein_tumor_vs_normal

            result = protein_tumor_vs_normal(
                cancer='brca',
                query='AKT1'
            )

            expected_cols = ['gene', 'log2_fold_change', 'p_value',
                           'p_value_adjusted', 'n_pairs', 'mean_tumor', 'mean_normal']

            assert_csv_structure(result['data'], expected_cols)


@pytest.mark.integration
class TestCancerTypeValidation:
    """Tests for cancer type validation in MCP tools."""

    def test_invalid_cancer_type_raises(self):
        """Test that invalid cancer type raises ValueError."""
        from cptac_mcp import phospho_tumor_vs_normal

        with pytest.raises(ValueError) as exc_info:
            phospho_tumor_vs_normal(
                cancer='invalid_cancer',
                query='AKT1_S473'
            )

        assert 'invalid_cancer' in str(exc_info.value).lower() or 'unsupported' in str(exc_info.value).lower()


@pytest.mark.integration
class TestEdgeCases:
    """Tests for edge cases and error handling."""

    def test_empty_query_string(self, mock_cancer_instance):
        """Test handling of empty query string."""
        with patch('cptac_mcp.CancerManager') as mock_manager:
            mock_manager.get.return_value = mock_cancer_instance
            mock_manager.available_cancers.return_value = ['brca']

            from cptac_mcp import phospho_tumor_vs_normal

            result = phospho_tumor_vs_normal(
                cancer='brca',
                query='',
                normalized='false'
            )

            assert 'error' in result

    def test_query_with_extra_whitespace(self, mock_cancer_instance, csv_parser):
        """Test query with extra whitespace is handled."""
        with patch('cptac_mcp.CancerManager') as mock_manager:
            mock_manager.get.return_value = mock_cancer_instance
            mock_manager.available_cancers.return_value = ['brca']

            from cptac_mcp import phospho_tumor_vs_normal

            result = phospho_tumor_vs_normal(
                cancer='brca',
                query='  AKT1_S473  ,  TP53_S15  ',
                normalized='false'
            )

            assert 'data' in result
            assert 'AKT1' in result['data']
            assert 'TP53' in result['data']

    def test_single_paired_sample(self, mock_cptac_cancer_obj):
        """Test handling when only one paired sample exists."""
        # Create DataFrame with only one paired sample
        index = pd.MultiIndex.from_tuples(
            [('AKT1', 'S473', 'PEP', 'DB')],
            names=['Name', 'Site', 'Peptide', 'Database_ID']
        )
        df = pd.DataFrame(
            [[1.0, 2.0]],
            index=index,
            columns=['S001', 'S001.N']
        )

        mock_cptac_cancer_obj.get_phosphoproteomics.return_value = df.T
        cancer = Cancer(name='test', cptac_obj=mock_cptac_cancer_obj)

        # Pre-populate cache
        cancer._phospho_deduplicated = df

        with patch('cptac_mcp.CancerManager') as mock_manager:
            mock_manager.get.return_value = cancer
            mock_manager.available_cancers.return_value = ['test']

            from cptac_mcp import phospho_tumor_vs_normal

            result = phospho_tumor_vs_normal(
                cancer='test',
                query='AKT1_S473',
                normalized='false'
            )

            # Should still work but p_value might be NaN with only 1 pair
            assert 'data' in result
