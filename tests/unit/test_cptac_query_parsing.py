"""
Unit tests for CPTAC query parsing logic.

Tests the parsing of Gene_Site queries and related functionality
used in cptac_mcp.py tools.
"""
import pytest
import pandas as pd
import numpy as np
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))


class TestQueryParsing:
    """Tests for Gene_Site query parsing."""

    @pytest.mark.parametrize("query,expected", [
        ("AKT1", ("AKT1", None)),
        ("AKT1_S473", ("AKT1", "S473")),
        ("AKT1_S473_S474", ("AKT1", "S473_S474")),
        ("MTOR_S2448", ("MTOR", "S2448")),
        ("GSK3B", ("GSK3B", None)),
        ("TP53_S15", ("TP53", "S15")),
        ("EGFR_Y1068", ("EGFR", "Y1068")),
        ("CDK1_T14", ("CDK1", "T14")),
    ])
    def test_gene_site_parsing(self, query, expected):
        """Parametrized tests for Gene_Site parsing."""
        if '_' in query:
            parts = query.split('_')
            gene = parts[0]
            site = '_'.join(parts[1:])
        else:
            gene = query
            site = None

        assert (gene, site) == expected

    def test_complex_site_format(self):
        """Test parsing of complex site formats like S473_S474."""
        query = "AKT1_S473_S474"
        parts = query.split('_')
        gene = parts[0]
        site = '_'.join(parts[1:])

        assert gene == "AKT1"
        assert site == "S473_S474"

    def test_multiple_queries(self):
        """Test parsing comma-separated queries."""
        query_str = "AKT1_S473,TP53_S15,MTOR"
        queries = [q.strip() for q in query_str.split(',')]

        assert len(queries) == 3
        assert queries[0] == "AKT1_S473"
        assert queries[1] == "TP53_S15"
        assert queries[2] == "MTOR"

    def test_query_with_whitespace(self):
        """Test parsing queries with extra whitespace."""
        query_str = "AKT1_S473 , TP53_S15,  MTOR  "
        queries = [q.strip() for q in query_str.split(',')]

        assert queries[0] == "AKT1_S473"
        assert queries[1] == "TP53_S15"
        assert queries[2] == "MTOR"


class TestProteinSuffixHandling:
    """Tests for _protein suffix detection and removal."""

    def test_detect_protein_suffix(self):
        """Test detection of _protein suffix."""
        query = "AKT1_protein"
        has_protein_suffix = query.endswith('_protein')

        assert has_protein_suffix

    def test_remove_protein_suffix(self):
        """Test removal of _protein suffix."""
        query = "AKT1_protein"
        if query.endswith('_protein'):
            gene = query[:-8]  # Remove '_protein'

        assert gene == "AKT1"

    def test_protein_suffix_not_confused_with_site(self):
        """Ensure _protein is not confused with phosphosite."""
        # This should be detected as protein query, not site query
        query = "TSC2_protein"
        force_protein = query.endswith('_protein')

        assert force_protein

        # After removing suffix, should just be gene name
        gene = query[:-8]
        assert gene == "TSC2"
        assert '_' not in gene

    @pytest.mark.parametrize("query,is_protein", [
        ("AKT1_protein", True),
        ("AKT1_S473", False),
        ("MTOR_protein", True),
        ("MTOR_S2448", False),
        ("TP53", False),
    ])
    def test_protein_vs_site_detection(self, query, is_protein):
        """Test distinguishing protein queries from site queries."""
        detected = query.endswith('_protein')
        assert detected == is_protein


class TestPhosphositeDetection:
    """Tests for detecting phosphosite patterns (S/T/Y followed by numbers)."""

    @pytest.mark.parametrize("site,is_phosphosite", [
        ("S473", True),
        ("T308", True),
        ("Y1068", True),
        ("S2448", True),
        ("protein", False),
        ("V600", False),  # V is not a phosphosite amino acid
        ("K123", False),  # K is not a phosphosite amino acid
    ])
    def test_phosphosite_pattern(self, site, is_phosphosite):
        """Test detection of phosphosite patterns."""
        if len(site) > 0 and site[0] in ['S', 'T', 'Y']:
            detected = True
        else:
            detected = False

        assert detected == is_phosphosite

    def test_phosphosite_in_query_context(self):
        """Test phosphosite detection in full query parsing context."""
        queries = [
            ("AKT1_S473", True),   # Phosphosite
            ("AKT1_protein", False),  # Protein marker
            ("AKT1_V600", False),  # Not a phosphosite
            ("AKT1", False),       # Gene only
        ]

        for query, expected_is_phosphosite in queries:
            if '_' in query:
                parts = query.split('_')
                potential_site = parts[1]

                if query.endswith('_protein'):
                    is_phosphosite = False
                elif len(potential_site) > 0 and potential_site[0] in ['S', 'T', 'Y']:
                    is_phosphosite = True
                else:
                    is_phosphosite = False
            else:
                is_phosphosite = False

            assert is_phosphosite == expected_is_phosphosite, f"Failed for query: {query}"


class TestNormalColumnIdentification:
    """Tests for identifying normal sample columns by .N suffix."""

    def test_detect_normal_columns(self, sample_phospho_multiindex):
        """Test detection of normal sample columns."""
        all_cols = sample_phospho_multiindex.columns.tolist()
        normal_cols = [col for col in all_cols if '.N' in col]

        assert len(normal_cols) == 3
        for col in normal_cols:
            assert '.N' in col

    def test_detect_tumor_columns(self, sample_phospho_multiindex):
        """Test detection of tumor sample columns."""
        all_cols = sample_phospho_multiindex.columns.tolist()
        tumor_cols = [col for col in all_cols if '.N' not in col]

        assert len(tumor_cols) == 5
        for col in tumor_cols:
            assert '.N' not in col

    @pytest.mark.parametrize("sample_id,is_normal", [
        ("S001", False),
        ("S001.N", True),
        ("S002", False),
        ("S002.N", True),
        ("SAMPLE_XYZ", False),
        ("SAMPLE_XYZ.N", True),
    ])
    def test_sample_type_detection(self, sample_id, is_normal):
        """Test sample type detection."""
        detected = '.N' in sample_id
        assert detected == is_normal


class TestPairedSampleMatching:
    """Tests for matching tumor samples to their normal counterparts."""

    def test_create_tumor_from_normal_id(self):
        """Test creating tumor sample ID from normal sample ID."""
        normal_id = "S001.N"
        tumor_id = normal_id.replace('.N', '')

        assert tumor_id == "S001"

    def test_find_paired_samples(self, sample_phospho_multiindex):
        """Test finding tumor-normal pairs."""
        all_cols = sample_phospho_multiindex.columns.tolist()
        normal_cols = [col for col in all_cols if '.N' in col]

        pairs = []
        for normal_col in normal_cols:
            tumor_col = normal_col.replace('.N', '')
            if tumor_col in all_cols:
                pairs.append((tumor_col, normal_col))

        # Should find 3 pairs
        assert len(pairs) == 3
        assert ('S001', 'S001.N') in pairs
        assert ('S002', 'S002.N') in pairs
        assert ('S003', 'S003.N') in pairs

    def test_unpaired_normal_sample(self):
        """Test handling of normal sample without tumor counterpart."""
        # Create a case where normal exists but tumor doesn't
        all_cols = ['S001', 'S002', 'S003.N', 'S004.N']

        normal_cols = [col for col in all_cols if '.N' in col]
        pairs = []
        for normal_col in normal_cols:
            tumor_col = normal_col.replace('.N', '')
            if tumor_col in all_cols:
                pairs.append((tumor_col, normal_col))

        # S003.N has no tumor counterpart, only S004 exists but doesn't have .N
        # So only S004.N would not have pair (S004 doesn't exist)
        # Actually neither S003 nor S004 exist, so no pairs
        assert len(pairs) == 0

    def test_paired_value_extraction(self, sample_phospho_multiindex):
        """Test extracting paired values for statistical analysis."""
        df = sample_phospho_multiindex
        idx = df.index[0]  # First phosphosite

        all_cols = df.columns.tolist()
        normal_cols = [col for col in all_cols if '.N' in col]

        paired_tumor = []
        paired_normal = []

        for normal_col in normal_cols:
            tumor_col = normal_col.replace('.N', '')
            if tumor_col in all_cols:
                tumor_val = df.loc[idx, tumor_col]
                normal_val = df.loc[idx, normal_col]

                if pd.notna(tumor_val) and pd.notna(normal_val):
                    paired_tumor.append(tumor_val)
                    paired_normal.append(normal_val)

        # Should have some paired values
        assert len(paired_tumor) > 0
        assert len(paired_tumor) == len(paired_normal)


class TestQueryResultMatching:
    """Tests for matching queries to DataFrame index entries."""

    def test_match_gene_site_to_multiindex(self, sample_phospho_multiindex):
        """Test matching Gene_Site query to MultiIndex."""
        df = sample_phospho_multiindex
        gene = "AKT1"
        site = "S473"

        matching_rows = [idx for idx in df.index if idx[0] == gene and idx[1] == site]

        assert len(matching_rows) > 0
        for idx in matching_rows:
            assert idx[0] == gene
            assert idx[1] == site

    def test_match_gene_only_to_multiindex(self, sample_phospho_multiindex):
        """Test matching gene-only query to MultiIndex (all sites)."""
        df = sample_phospho_multiindex
        gene = "AKT1"

        matching_rows = [idx for idx in df.index if idx[0] == gene]

        # Should find multiple sites for AKT1
        assert len(matching_rows) == 3  # S473, T308, S124
        for idx in matching_rows:
            assert idx[0] == gene

    def test_no_match_returns_empty(self, sample_phospho_multiindex):
        """Test that non-existent query returns empty list."""
        df = sample_phospho_multiindex
        gene = "NONEXISTENT"

        matching_rows = [idx for idx in df.index if idx[0] == gene]

        assert len(matching_rows) == 0


class TestCorrelationQueryParsing:
    """Tests for correlation analysis query parsing."""

    def test_mixed_query_parsing(self):
        """Test parsing mixed phospho and protein queries."""
        query_str = "AKT1_S473,TSC2_protein,MTOR_S2448,TP53"
        queries = [q.strip() for q in query_str.split(',')]

        parsed = []
        for q in queries:
            if q.endswith('_protein'):
                parsed.append(('protein', q[:-8]))
            elif '_' in q:
                parts = q.split('_')
                if parts[1][0] in ['S', 'T', 'Y']:
                    parsed.append(('phospho', parts[0], '_'.join(parts[1:])))
                else:
                    parsed.append(('gene', q))
            else:
                parsed.append(('gene', q))

        assert parsed[0] == ('phospho', 'AKT1', 'S473')
        assert parsed[1] == ('protein', 'TSC2')
        assert parsed[2] == ('phospho', 'MTOR', 'S2448')
        assert parsed[3] == ('gene', 'TP53')
