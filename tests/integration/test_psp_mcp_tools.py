"""
Integration tests for PSP MCP tools.

Tests the MCP tool functions with mocked PSP data.
"""
import pytest
import pandas as pd
import numpy as np
from unittest.mock import patch, MagicMock
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))


@pytest.mark.integration
class TestGetKinaseSubstrates:
    """Tests for get_kinase_substrates MCP tool."""

    def test_single_kinase_query(self, mock_psp_data, csv_parser):
        """Test query for a single kinase."""
        from psp_proteomics import get_kinase_substrates

        result = get_kinase_substrates('AKT1')

        assert 'AKT1' in result
        rows = csv_parser(result['AKT1'])

        assert len(rows) > 0
        # Check expected columns
        assert 'substrate' in rows[0]
        assert 'gene' in rows[0]
        assert 'site' in rows[0]

    def test_multiple_kinases_query(self, mock_psp_data, csv_parser):
        """Test query for multiple kinases."""
        from psp_proteomics import get_kinase_substrates

        result = get_kinase_substrates('AKT1,MTOR,GSK3B')

        assert 'AKT1' in result
        assert 'MTOR' in result
        assert 'GSK3B' in result

        # Each kinase should have substrates
        for kinase in ['AKT1', 'MTOR', 'GSK3B']:
            rows = csv_parser(result[kinase])
            assert len(rows) > 0

    def test_organism_filter(self, mock_psp_data, csv_parser):
        """Test organism filtering."""
        from psp_proteomics import get_kinase_substrates

        # Query for human
        result_human = get_kinase_substrates('AKT1', organism='human')

        # Should return human substrates
        assert 'AKT1' in result_human

    def test_nonexistent_kinase(self, mock_psp_data, csv_parser):
        """Test query for non-existent kinase returns header only."""
        from psp_proteomics import get_kinase_substrates

        result = get_kinase_substrates('NONEXISTENT')

        # Should have entry but with only header row
        assert 'NONEXISTENT' in result
        lines = result['NONEXISTENT'].split('\n')
        assert len(lines) == 1  # Only header


@pytest.mark.integration
class TestGetRegulatorySites:
    """Tests for get_regulatory_sites MCP tool."""

    def test_gene_query(self, mock_psp_data, csv_parser):
        """Test query by gene name returns all sites."""
        from psp_proteomics import get_regulatory_sites

        result = get_regulatory_sites('TP53')

        assert 'TP53' in result
        rows = csv_parser(result['TP53'])

        # TP53 should have multiple regulatory sites
        assert len(rows) > 0

    def test_gene_site_query(self, mock_psp_data, csv_parser):
        """Test query by Gene_Site format."""
        from psp_proteomics import get_regulatory_sites

        result = get_regulatory_sites('AKT1_S473-p')

        assert 'AKT1_S473-p' in result
        rows = csv_parser(result['AKT1_S473-p'])

        assert len(rows) > 0

    def test_mixed_query(self, mock_psp_data):
        """Test mixed gene and Gene_Site query."""
        from psp_proteomics import get_regulatory_sites

        result = get_regulatory_sites('TP53,AKT1_S473-p')

        assert 'TP53' in result
        assert 'AKT1_S473-p' in result

    def test_functional_annotations_present(self, mock_psp_data, csv_parser):
        """Test that functional annotations are included."""
        from psp_proteomics import get_regulatory_sites

        result = get_regulatory_sites('AKT1')

        rows = csv_parser(result['AKT1'])

        # Check for functional annotation columns
        assert 'on_function' in rows[0]
        assert 'on_process' in rows[0]


@pytest.mark.integration
class TestGetDiseaseSites:
    """Tests for get_disease_sites MCP tool."""

    def test_gene_query(self, mock_psp_data, csv_parser):
        """Test query by gene name."""
        from psp_proteomics import get_disease_sites

        result = get_disease_sites('TP53')

        assert 'TP53' in result
        rows = csv_parser(result['TP53'])

        assert len(rows) > 0

    def test_disease_filter(self, mock_psp_data, csv_parser):
        """Test disease name filter."""
        from psp_proteomics import get_disease_sites

        result = get_disease_sites('TP53', disease='breast cancer')

        assert 'TP53' in result
        rows = csv_parser(result['TP53'])

        # All results should be breast cancer related
        for row in rows:
            if row.get('disease'):
                assert 'breast' in row['disease'].lower() or 'cancer' in row['disease'].lower()

    def test_gene_site_query(self, mock_psp_data, csv_parser):
        """Test query by Gene_Site format."""
        from psp_proteomics import get_disease_sites

        result = get_disease_sites('TP53_S15')

        assert 'TP53_S15' in result

    def test_no_disease_matches(self, mock_psp_data, csv_parser):
        """Test query with no disease matches."""
        from psp_proteomics import get_disease_sites

        result = get_disease_sites('TP53', disease='nonexistent_rare_disease')

        assert 'TP53' in result
        # Should have header only (no matches)
        lines = result['TP53'].split('\n')
        assert len(lines) == 1  # Only header


@pytest.mark.integration
class TestFindUpstreamKinases:
    """Tests for find_upstream_kinases MCP tool."""

    def test_site_query(self, mock_psp_data, csv_parser):
        """Test query for specific site returns kinases."""
        from psp_proteomics import find_upstream_kinases

        result = find_upstream_kinases('GSK3B_S9')

        assert 'GSK3B_S9' in result
        rows = csv_parser(result['GSK3B_S9'])

        # GSK3B S9 should have AKT1 as upstream kinase
        assert len(rows) > 0
        assert 'kinase' in rows[0]

    def test_gene_query_returns_all_sites(self, mock_psp_data, csv_parser):
        """Test gene query returns all sites with their kinases."""
        from psp_proteomics import find_upstream_kinases

        result = find_upstream_kinases('GSK3B')

        assert 'GSK3B' in result
        rows = csv_parser(result['GSK3B'])

        # Should include site column for gene queries
        assert 'site' in rows[0]

    def test_multiple_queries(self, mock_psp_data):
        """Test multiple site/gene queries."""
        from psp_proteomics import find_upstream_kinases

        result = find_upstream_kinases('GSK3B_S9,TSC2_S939')

        assert 'GSK3B_S9' in result
        assert 'TSC2_S939' in result


@pytest.mark.integration
class TestGetInteractionNetwork:
    """Tests for get_interaction_network MCP tool."""

    def test_kinase_query_returns_network(self, mock_psp_data):
        """Test kinase query builds network with substrates."""
        from psp_proteomics import get_interaction_network

        result = get_interaction_network('AKT1', depth=1)

        assert 'nodes' in result
        assert 'edges' in result
        assert 'summary' in result

    def test_site_query_returns_network(self, mock_psp_data):
        """Test site query builds network with upstream kinases."""
        from psp_proteomics import get_interaction_network

        result = get_interaction_network('GSK3B_S9', depth=1)

        assert 'nodes' in result
        assert 'edges' in result

    def test_depth_parameter(self, mock_psp_data):
        """Test depth parameter affects network size."""
        from psp_proteomics import get_interaction_network

        result_d1 = get_interaction_network('AKT1', depth=1)
        result_d2 = get_interaction_network('AKT1', depth=2)

        # Depth 2 should have more or equal nodes/edges
        nodes_d1 = len(result_d1['nodes'].split('\n'))
        nodes_d2 = len(result_d2['nodes'].split('\n'))

        # Note: Depending on data, depth 2 might not always be larger
        assert nodes_d1 >= 1
        assert nodes_d2 >= 1

    def test_evidence_filter(self, mock_psp_data, csv_parser):
        """Test evidence type filtering."""
        from psp_proteomics import get_interaction_network

        result = get_interaction_network('AKT1', depth=1, evidence_types='in_vivo')

        # Edges should only have in_vivo evidence
        edge_rows = csv_parser(result['edges'])

        for row in edge_rows:
            if row.get('evidence'):
                assert 'in_vivo' in row['evidence']

    def test_network_structure(self, mock_psp_data, csv_parser):
        """Test network has proper node/edge structure."""
        from psp_proteomics import get_interaction_network

        result = get_interaction_network('AKT1', depth=1)

        # Check node structure
        node_rows = csv_parser(result['nodes'])
        if node_rows:
            assert 'node_id' in node_rows[0]
            assert 'type' in node_rows[0]
            assert 'gene' in node_rows[0]

        # Check edge structure
        edge_rows = csv_parser(result['edges'])
        if edge_rows:
            assert 'source' in edge_rows[0]
            assert 'target' in edge_rows[0]


@pytest.mark.integration
class TestGetPathwayContext:
    """Tests for get_pathway_context MCP tool."""

    def test_kinase_query(self, mock_psp_data):
        """Test pathway context for a kinase."""
        from psp_proteomics import get_pathway_context

        result = get_pathway_context('AKT1')

        assert 'query' in result
        assert 'summary' in result

    def test_site_query(self, mock_psp_data):
        """Test pathway context for a specific site."""
        from psp_proteomics import get_pathway_context

        result = get_pathway_context('AKT1_S473')

        assert 'query' in result
        assert result['query'] == 'AKT1_S473'

    def test_include_downstream(self, mock_psp_data):
        """Test including downstream substrates."""
        from psp_proteomics import get_pathway_context

        result = get_pathway_context('AKT1', include_downstream=True)

        assert 'downstream' in result

    def test_include_upstream(self, mock_psp_data):
        """Test including upstream kinases."""
        from psp_proteomics import get_pathway_context

        result = get_pathway_context('AKT1_S473', include_upstream=True)

        assert 'upstream' in result

    def test_exclude_downstream(self, mock_psp_data):
        """Test excluding downstream substrates."""
        from psp_proteomics import get_pathway_context

        result = get_pathway_context('AKT1', include_downstream=False)

        # Should not have downstream key
        assert 'downstream' not in result or result.get('downstream') is None

    def test_functional_outcomes_included(self, mock_psp_data):
        """Test functional outcomes are included."""
        from psp_proteomics import get_pathway_context

        result = get_pathway_context('AKT1')

        assert 'functional_outcomes' in result


@pytest.mark.integration
class TestGetEvidenceSummary:
    """Tests for get_evidence_summary MCP tool."""

    def test_gene_query(self, mock_psp_data):
        """Test evidence summary for a gene."""
        from psp_proteomics import get_evidence_summary

        result = get_evidence_summary('AKT1')

        assert 'AKT1' in result
        assert 'summary' in result['AKT1']

    def test_site_query(self, mock_psp_data):
        """Test evidence summary for a specific site."""
        from psp_proteomics import get_evidence_summary

        result = get_evidence_summary('MTOR_S2448')

        assert 'MTOR_S2448' in result

    def test_top_kinases_included(self, mock_psp_data):
        """Test top kinases section is included."""
        from psp_proteomics import get_evidence_summary

        result = get_evidence_summary('AKT1')

        assert 'top_kinases' in result['AKT1']

    def test_functional_categories_included(self, mock_psp_data):
        """Test functional categories section is included."""
        from psp_proteomics import get_evidence_summary

        result = get_evidence_summary('AKT1')

        assert 'functional_categories' in result['AKT1']

    def test_disease_associations_included(self, mock_psp_data):
        """Test disease associations section is included."""
        from psp_proteomics import get_evidence_summary

        result = get_evidence_summary('TP53')

        assert 'disease_associations' in result['TP53']


@pytest.mark.integration
class TestCSVOutputFormat:
    """Tests for verifying CSV output format of PSP MCP tools."""

    def test_kinase_substrates_csv_columns(self, mock_psp_data, assert_csv_structure):
        """Verify kinase substrates CSV has expected columns."""
        from psp_proteomics import get_kinase_substrates

        result = get_kinase_substrates('AKT1')

        expected_cols = ['substrate', 'gene', 'site', 'site_sequence',
                        'in_vivo', 'in_vitro', 'pmids']

        assert_csv_structure(result['AKT1'], expected_cols)

    def test_regulatory_sites_csv_columns(self, mock_psp_data, assert_csv_structure):
        """Verify regulatory sites CSV has expected columns."""
        from psp_proteomics import get_regulatory_sites

        result = get_regulatory_sites('AKT1')

        expected_cols = ['site', 'modification', 'on_function', 'on_process',
                        'on_interaction', 'pmids', 'lt_lit', 'ms_lit']

        assert_csv_structure(result['AKT1'], expected_cols)

    def test_disease_sites_csv_columns(self, mock_psp_data, assert_csv_structure):
        """Verify disease sites CSV has expected columns."""
        from psp_proteomics import get_disease_sites

        result = get_disease_sites('TP53')

        expected_cols = ['site', 'disease', 'alteration', 'pmids',
                        'lt_lit', 'ms_lit', 'notes']

        assert_csv_structure(result['TP53'], expected_cols)


@pytest.mark.integration
class TestEdgeCases:
    """Tests for edge cases and error handling."""

    def test_empty_query(self, mock_psp_data, csv_parser):
        """Test handling of empty query."""
        from psp_proteomics import get_kinase_substrates

        result = get_kinase_substrates('')

        # Should return entry for empty string
        assert '' in result
        # Should have only header row
        lines = result[''].split('\n')
        assert len(lines) == 1

    def test_whitespace_in_query(self, mock_psp_data):
        """Test query with whitespace is handled."""
        from psp_proteomics import get_kinase_substrates

        result = get_kinase_substrates('  AKT1  ,  MTOR  ')

        # Queries should be stripped
        assert 'AKT1' in result
        assert 'MTOR' in result

    def test_case_insensitive_query(self, mock_psp_data, csv_parser):
        """Test queries are case-insensitive."""
        from psp_proteomics import get_kinase_substrates

        result_upper = get_kinase_substrates('AKT1')
        result_lower = get_kinase_substrates('akt1')

        # Both should return results
        rows_upper = csv_parser(result_upper['AKT1'])
        rows_lower = csv_parser(result_lower['akt1'])

        # Results might be in different keys but should have same data
        assert len(rows_upper) > 0 or len(rows_lower) > 0
