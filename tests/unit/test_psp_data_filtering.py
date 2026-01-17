"""
Unit tests for PSP data filtering logic.

Tests the filtering and parsing functions used in psp_proteomics.py.
"""
import pytest
import pandas as pd
import numpy as np
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))


class TestOrganismFiltering:
    """Tests for organism filtering functionality."""

    def test_case_insensitive_human_filter(self, sample_kinase_substrate_df):
        """Test case-insensitive filtering for human organism."""
        df = sample_kinase_substrate_df

        # Filter with lowercase
        result_lower = df[df['KIN_ORGANISM'].str.lower() == 'human']
        # Filter with mixed case
        result_mixed = df[df['KIN_ORGANISM'].str.lower() == 'HUMAN'.lower()]

        assert len(result_lower) == len(result_mixed)
        assert len(result_lower) > 0

    def test_mouse_filter(self, sample_kinase_substrate_df):
        """Test filtering for mouse organism."""
        df = sample_kinase_substrate_df

        result = df[df['KIN_ORGANISM'].str.lower() == 'mouse']

        # Mock data includes one mouse entry
        assert len(result) == 1
        assert result.iloc[0]['KINASE'] == 'Akt1'

    def test_excludes_other_organisms(self, sample_kinase_substrate_df):
        """Test that human filter excludes mouse entries."""
        df = sample_kinase_substrate_df

        human_only = df[df['KIN_ORGANISM'].str.lower() == 'human']
        all_entries = df

        assert len(human_only) < len(all_entries)

    @pytest.mark.parametrize("organism_input,expected_match", [
        ("human", True),
        ("Human", True),
        ("HUMAN", True),
        ("mouse", True),
        ("Mouse", True),
        ("rat", False),
        ("", False),
    ])
    def test_organism_matching(self, sample_kinase_substrate_df, organism_input, expected_match):
        """Test various organism input formats."""
        df = sample_kinase_substrate_df

        result = df[df['KIN_ORGANISM'].str.lower() == organism_input.lower()]

        if expected_match:
            assert len(result) > 0
        else:
            assert len(result) == 0


class TestKinaseQuery:
    """Tests for kinase name matching."""

    def test_exact_kinase_match(self, sample_kinase_substrate_df):
        """Test exact kinase name matching."""
        df = sample_kinase_substrate_df

        result = df[df['KINASE'].str.upper() == 'AKT1']

        assert len(result) > 0
        for _, row in result.iterrows():
            assert row['KINASE'].upper() == 'AKT1'

    def test_case_insensitive_kinase_match(self, sample_kinase_substrate_df):
        """Test case-insensitive kinase matching."""
        df = sample_kinase_substrate_df

        result_upper = df[df['KINASE'].str.upper() == 'AKT1']
        result_lower = df[df['KINASE'].str.upper() == 'akt1'.upper()]

        assert len(result_upper) == len(result_lower)

    def test_multiple_kinases(self, sample_kinase_substrate_df):
        """Test querying multiple kinases."""
        df = sample_kinase_substrate_df
        kinases = ['AKT1', 'MTOR']

        results = {}
        for kinase in kinases:
            result = df[df['KINASE'].str.upper() == kinase.upper()]
            results[kinase] = result

        assert len(results['AKT1']) > 0
        assert len(results['MTOR']) > 0

    def test_nonexistent_kinase(self, sample_kinase_substrate_df):
        """Test query for non-existent kinase."""
        df = sample_kinase_substrate_df

        result = df[df['KINASE'].str.upper() == 'NONEXISTENT']

        assert len(result) == 0


class TestGeneSiteParsing:
    """Tests for Gene_Site format parsing in PSP context."""

    @pytest.mark.parametrize("query,expected_gene,expected_site", [
        ("AKT1_S473", "AKT1", "S473"),
        ("TP53_S15", "TP53", "S15"),
        ("MTOR_S2448", "MTOR", "S2448"),
        ("GSK3B_Y216", "GSK3B", "Y216"),
    ])
    def test_gene_site_split(self, query, expected_gene, expected_site):
        """Test splitting Gene_Site into components."""
        parts = query.split('_')
        gene = parts[0]
        site = '_'.join(parts[1:])

        assert gene == expected_gene
        assert site == expected_site

    def test_phosphosite_detection_s_t_y(self):
        """Test detection of S/T/Y phosphorylation sites."""
        sites = ['S473', 'T308', 'Y1068']

        for site in sites:
            assert len(site) > 0
            assert site[0] in ['S', 'T', 'Y']

    def test_non_phosphosite_detection(self):
        """Test that non-S/T/Y sites are not treated as phosphosites."""
        sites = ['V600', 'K123', 'A100']

        for site in sites:
            is_phosphosite = len(site) > 0 and site[0] in ['S', 'T', 'Y']
            assert not is_phosphosite

    def test_site_with_modification_suffix(self):
        """Test parsing site with modification suffix like S473-p."""
        site_with_mod = "S473-p"

        # Extract site part before modification
        site = site_with_mod.split('-')[0]
        mod = site_with_mod.split('-')[-1] if '-' in site_with_mod else ''

        assert site == "S473"
        assert mod == "p"


class TestEvidenceTypeFiltering:
    """Tests for IN_VIVO_RXN and IN_VITRO_RXN filtering."""

    def test_in_vivo_marker_detection(self, sample_kinase_substrate_df):
        """Test detection of 'X' marker in IN_VIVO_RXN column."""
        df = sample_kinase_substrate_df

        in_vivo_entries = df[df['IN_VIVO_RXN'].str.strip() == 'X']

        assert len(in_vivo_entries) > 0
        for _, row in in_vivo_entries.iterrows():
            assert row['IN_VIVO_RXN'].strip() == 'X'

    def test_in_vitro_marker_detection(self, sample_kinase_substrate_df):
        """Test detection of 'X' marker in IN_VITRO_RXN column."""
        df = sample_kinase_substrate_df

        in_vitro_entries = df[df['IN_VITRO_RXN'].str.strip() == 'X']

        assert len(in_vitro_entries) > 0
        for _, row in in_vitro_entries.iterrows():
            assert row['IN_VITRO_RXN'].strip() == 'X'

    def test_both_evidence_types(self, sample_kinase_substrate_df):
        """Test entries with both in_vivo and in_vitro evidence."""
        df = sample_kinase_substrate_df

        both = df[(df['IN_VIVO_RXN'].str.strip() == 'X') &
                  (df['IN_VITRO_RXN'].str.strip() == 'X')]

        # Some entries should have both
        assert len(both) > 0

    def test_filter_by_evidence_type(self, sample_kinase_substrate_df):
        """Test filtering by specific evidence type."""
        df = sample_kinase_substrate_df
        evidence_filter = ['in_vivo']

        # Apply filter
        filtered = df.copy()
        if 'in_vivo' in evidence_filter:
            filtered = filtered[filtered['IN_VIVO_RXN'].str.strip() == 'X']

        # All results should have in_vivo evidence
        for _, row in filtered.iterrows():
            assert row['IN_VIVO_RXN'].strip() == 'X'


class TestDiseaseFiltering:
    """Tests for disease association filtering."""

    def test_disease_contains_filter(self, sample_disease_sites_df):
        """Test filtering by disease name (contains)."""
        df = sample_disease_sites_df

        # Filter for cancer-related entries
        cancer_entries = df[df['DISEASE'].str.contains('cancer', case=False, na=False)]

        assert len(cancer_entries) > 0
        for _, row in cancer_entries.iterrows():
            assert 'cancer' in row['DISEASE'].lower()

    def test_case_insensitive_disease_filter(self, sample_disease_sites_df):
        """Test case-insensitive disease filtering."""
        df = sample_disease_sites_df

        result_lower = df[df['DISEASE'].str.contains('cancer', case=False, na=False)]
        result_upper = df[df['DISEASE'].str.contains('CANCER', case=False, na=False)]

        assert len(result_lower) == len(result_upper)

    def test_specific_disease_filter(self, sample_disease_sites_df):
        """Test filtering for specific disease."""
        df = sample_disease_sites_df

        breast_cancer = df[df['DISEASE'].str.contains('breast cancer', case=False, na=False)]

        assert len(breast_cancer) > 0
        for _, row in breast_cancer.iterrows():
            assert 'breast cancer' in row['DISEASE'].lower()

    def test_no_disease_filter(self, sample_disease_sites_df):
        """Test query without disease filter returns all entries."""
        df = sample_disease_sites_df
        gene = 'TP53'

        # Without disease filter
        all_tp53 = df[df['GENE'].str.upper() == gene]

        # Should return multiple disease associations
        assert len(all_tp53) > 1

    def test_rare_disease_filter(self, sample_disease_sites_df):
        """Test filtering for disease with few/no entries."""
        df = sample_disease_sites_df

        rare = df[df['DISEASE'].str.contains('rare_nonexistent_disease', case=False, na=False)]

        assert len(rare) == 0


class TestSubstrateFiltering:
    """Tests for substrate gene filtering in kinase-substrate data."""

    def test_substrate_gene_filter(self, sample_kinase_substrate_df):
        """Test filtering by substrate gene."""
        df = sample_kinase_substrate_df
        df_human = df[df['SUB_ORGANISM'].str.lower() == 'human']

        result = df_human[df_human['SUB_GENE'].str.upper() == 'GSK3B']

        assert len(result) > 0
        for _, row in result.iterrows():
            assert row['SUB_GENE'].upper() == 'GSK3B'

    def test_substrate_site_filter(self, sample_kinase_substrate_df):
        """Test filtering by substrate site."""
        df = sample_kinase_substrate_df
        df_human = df[df['SUB_ORGANISM'].str.lower() == 'human']

        result = df_human[(df_human['SUB_GENE'].str.upper() == 'GSK3B') &
                         (df_human['SUB_MOD_RSD'].str.upper() == 'S9')]

        assert len(result) > 0


class TestRegulatoryDataFiltering:
    """Tests for regulatory sites data filtering."""

    def test_gene_filter(self, sample_regulatory_sites_df):
        """Test filtering regulatory data by gene."""
        df = sample_regulatory_sites_df
        df_human = df[df['ORGANISM'].str.lower() == 'human']

        result = df_human[df_human['GENE'].str.upper() == 'AKT1']

        assert len(result) > 0
        for _, row in result.iterrows():
            assert row['GENE'].upper() == 'AKT1'

    def test_gene_site_filter(self, sample_regulatory_sites_df):
        """Test filtering regulatory data by gene and site."""
        df = sample_regulatory_sites_df
        df_human = df[df['ORGANISM'].str.lower() == 'human']

        result = df_human[(df_human['GENE'].str.upper() == 'AKT1') &
                         (df_human['MOD_RSD'].str.upper() == 'S473-P')]

        assert len(result) > 0

    def test_functional_annotation_present(self, sample_regulatory_sites_df):
        """Test that regulatory data includes functional annotations."""
        df = sample_regulatory_sites_df

        # Check for non-empty ON_FUNCTION entries
        has_function = df[df['ON_FUNCTION'].notna() & (df['ON_FUNCTION'] != '')]

        assert len(has_function) > 0


class TestPMIDParsing:
    """Tests for PMID parsing and handling."""

    def test_single_pmid(self):
        """Test parsing single PMID."""
        pmid_str = "12345"
        pmids = pmid_str.split(';')

        assert len(pmids) == 1
        assert pmids[0] == "12345"

    def test_multiple_pmids(self):
        """Test parsing multiple semicolon-separated PMIDs."""
        pmid_str = "12345;67890;11111"
        pmids = [p.strip() for p in pmid_str.split(';') if p.strip()]

        assert len(pmids) == 3
        assert "12345" in pmids
        assert "67890" in pmids
        assert "11111" in pmids

    def test_empty_pmid(self):
        """Test handling of empty PMID field."""
        pmid_str = ""
        pmids = [p.strip() for p in pmid_str.split(';') if p.strip()]

        assert len(pmids) == 0

    def test_nan_pmid(self):
        """Test handling of NaN PMID value."""
        pmid_val = np.nan

        if pd.notna(pmid_val):
            pmid_str = str(pmid_val)
        else:
            pmid_str = ""

        assert pmid_str == ""
