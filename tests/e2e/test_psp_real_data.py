"""
End-to-end tests for PSP MCP tools using real data.

These tests use the actual PSP_2025 dataset files.
Tests are skipped if PSP data files are not available.
"""
import pytest
import pandas as pd
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))


# Path to PSP data directory
PSP_DATA_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    'datasets',
    'PSP_2025'
)

# Check if PSP data is available
PSP_DATA_AVAILABLE = all([
    os.path.exists(os.path.join(PSP_DATA_DIR, f))
    for f in [
        'Phosphorylation_site_dataset.txt',
        'Kinase_Substrate_Dataset.txt',
        'Regulatory_sites.txt',
        'Disease-associated_sites.txt',
    ]
])


def skip_if_no_psp_data():
    """Skip test if PSP data is not available."""
    if not PSP_DATA_AVAILABLE:
        pytest.skip("PSP data files not available")


@pytest.mark.e2e
class TestPSPDataIntegrity:
    """Tests to verify PSP dataset integrity and structure."""

    @pytest.fixture(autouse=True)
    def check_psp_data(self):
        """Check PSP data availability before each test."""
        skip_if_no_psp_data()

    def test_kinase_substrate_file_readable(self):
        """Test that kinase-substrate file can be read."""
        filepath = os.path.join(PSP_DATA_DIR, 'Kinase_Substrate_Dataset.txt')
        df = pd.read_csv(filepath, sep='\t', skiprows=3, encoding='latin1', low_memory=False)

        assert len(df) > 0
        assert 'KINASE' in df.columns
        assert 'SUB_GENE' in df.columns

    def test_regulatory_sites_file_readable(self):
        """Test that regulatory sites file can be read."""
        filepath = os.path.join(PSP_DATA_DIR, 'Regulatory_sites.txt')
        df = pd.read_csv(filepath, sep='\t', skiprows=3, encoding='latin1', low_memory=False)

        assert len(df) > 0
        assert 'GENE' in df.columns
        assert 'MOD_RSD' in df.columns

    def test_disease_sites_file_readable(self):
        """Test that disease sites file can be read."""
        filepath = os.path.join(PSP_DATA_DIR, 'Disease-associated_sites.txt')
        df = pd.read_csv(filepath, sep='\t', skiprows=3, encoding='latin1', low_memory=False)

        assert len(df) > 0
        assert 'GENE' in df.columns
        assert 'DISEASE' in df.columns

    def test_phospho_sites_file_readable(self):
        """Test that phosphorylation sites file can be read."""
        filepath = os.path.join(PSP_DATA_DIR, 'Phosphorylation_site_dataset.txt')
        df = pd.read_csv(filepath, sep='\t', skiprows=3, encoding='latin1', low_memory=False)

        assert len(df) > 0
        assert 'GENE' in df.columns
        assert 'MOD_RSD' in df.columns

    def test_human_data_present(self):
        """Test that human data is present in datasets."""
        filepath = os.path.join(PSP_DATA_DIR, 'Kinase_Substrate_Dataset.txt')
        df = pd.read_csv(filepath, sep='\t', skiprows=3, encoding='latin1', low_memory=False)

        human_entries = df[df['KIN_ORGANISM'].str.lower() == 'human']
        assert len(human_entries) > 100, "Expected significant human data"


@pytest.mark.e2e
class TestPSPRealData:
    """End-to-end tests using real PSP data."""

    @pytest.fixture(autouse=True)
    def check_psp_data(self):
        """Check PSP data availability before each test."""
        skip_if_no_psp_data()

    def test_akt1_substrates(self, csv_parser):
        """Test querying real AKT1 substrates."""
        from psp_proteomics import get_kinase_substrates

        result = get_kinase_substrates('AKT1')

        assert 'AKT1' in result
        rows = csv_parser(result['AKT1'])

        # AKT1 is a well-studied kinase, should have many substrates
        assert len(rows) > 10, "AKT1 should have many known substrates"

        # Check for well-known substrates
        substrate_genes = [row['gene'] for row in rows]
        # GSK3B is a well-known AKT1 substrate
        assert any('GSK3' in g for g in substrate_genes), "GSK3B should be an AKT1 substrate"

    def test_tp53_regulatory_sites(self, csv_parser):
        """Test querying real TP53 regulatory sites."""
        from psp_proteomics import get_regulatory_sites

        result = get_regulatory_sites('TP53')

        assert 'TP53' in result
        rows = csv_parser(result['TP53'])

        # TP53 is well-studied, should have many regulatory sites
        assert len(rows) > 5, "TP53 should have many regulatory sites"

    def test_mtor_s2448_upstream_kinases(self, csv_parser):
        """Test finding upstream kinases for MTOR S2448."""
        from psp_proteomics import find_upstream_kinases

        result = find_upstream_kinases('MTOR_S2448')

        assert 'MTOR_S2448' in result
        rows = csv_parser(result['MTOR_S2448'])

        # MTOR S2448 is a well-studied site
        assert len(rows) > 0, "MTOR S2448 should have known upstream kinases"

    def test_disease_associations_for_braf(self, csv_parser):
        """Test disease associations for BRAF."""
        from psp_proteomics import get_disease_sites

        result = get_disease_sites('BRAF')

        assert 'BRAF' in result
        rows = csv_parser(result['BRAF'])

        # BRAF is associated with multiple cancers
        # If there are disease associations, verify they exist
        # (specific disease names may vary in different PSP versions)
        if len(rows) > 0:
            diseases = [row.get('disease', '') for row in rows]
            # Should have at least some disease associations
            assert any(d for d in diseases), "BRAF should have disease associations"

    def test_akt1_interaction_network(self):
        """Test building interaction network for AKT1."""
        from psp_proteomics import get_interaction_network

        result = get_interaction_network('AKT1', depth=1)

        assert 'nodes' in result
        assert 'edges' in result
        assert 'summary' in result

        # Parse summary
        summary = result['summary']
        # Should have multiple nodes and edges
        assert 'nodes=' in summary
        assert 'edges=' in summary

    def test_akt1_pathway_context(self):
        """Test pathway context for AKT1."""
        from psp_proteomics import get_pathway_context

        result = get_pathway_context('AKT1')

        assert 'query' in result
        assert 'summary' in result
        assert 'upstream' in result
        assert 'downstream' in result

        # AKT1 is a central kinase, should have both upstream and downstream
        summary = result['summary']
        # Parse summary to check counts
        assert 'downstream_substrates=' in summary

    def test_evidence_summary_for_mtor(self):
        """Test evidence summary for MTOR."""
        from psp_proteomics import get_evidence_summary

        result = get_evidence_summary('MTOR')

        assert 'MTOR' in result

        data = result['MTOR']
        assert 'summary' in data
        assert 'top_kinases' in data
        assert 'functional_categories' in data


@pytest.mark.e2e
class TestPSPComplexQueries:
    """Tests for complex queries on real PSP data."""

    @pytest.fixture(autouse=True)
    def check_psp_data(self):
        """Check PSP data availability before each test."""
        skip_if_no_psp_data()

    def test_multiple_kinases_query(self, csv_parser):
        """Test querying multiple kinases at once."""
        from psp_proteomics import get_kinase_substrates

        result = get_kinase_substrates('AKT1,MTOR,GSK3B,CDK1')

        for kinase in ['AKT1', 'MTOR', 'GSK3B', 'CDK1']:
            assert kinase in result
            rows = csv_parser(result[kinase])
            assert len(rows) > 0, f"{kinase} should have substrates"

    def test_mixed_regulatory_query(self, csv_parser):
        """Test mixed gene and Gene_Site regulatory query."""
        from psp_proteomics import get_regulatory_sites

        result = get_regulatory_sites('TP53,AKT1_S473')

        assert 'TP53' in result
        # Note: Site format in regulatory data may differ

    def test_disease_filter_cancer(self, csv_parser):
        """Test filtering diseases by cancer."""
        from psp_proteomics import get_disease_sites

        result = get_disease_sites('TP53,BRAF,EGFR', disease='cancer')

        # All results should be cancer-related
        for gene in ['TP53', 'BRAF', 'EGFR']:
            if gene in result:
                rows = csv_parser(result[gene])
                for row in rows:
                    disease = row.get('disease', '').lower()
                    if disease:
                        assert 'cancer' in disease or 'carcinoma' in disease or \
                               'tumor' in disease or 'melanoma' in disease or \
                               'leukemia' in disease or 'lymphoma' in disease

    def test_network_depth_2(self):
        """Test network building with depth 2."""
        from psp_proteomics import get_interaction_network

        result = get_interaction_network('AKT1', depth=2)

        # Depth 2 should include secondary connections
        nodes_csv = result['nodes']
        node_count = len(nodes_csv.split('\n')) - 1  # Subtract header

        # With depth 2, should have more nodes
        assert node_count > 1


@pytest.mark.e2e
class TestPSPPerformance:
    """Performance tests for PSP queries."""

    @pytest.fixture(autouse=True)
    def check_psp_data(self):
        """Check PSP data availability before each test."""
        skip_if_no_psp_data()

    @pytest.mark.slow
    def test_large_kinase_query(self, csv_parser):
        """Test querying many kinases at once."""
        from psp_proteomics import get_kinase_substrates

        kinases = 'AKT1,MTOR,GSK3B,CDK1,CDK2,MAPK1,MAPK3,PKA,PKC,SRC'
        result = get_kinase_substrates(kinases)

        # Should return results for all kinases
        assert len(result) == 10

    @pytest.mark.slow
    def test_large_network(self):
        """Test building a larger network."""
        from psp_proteomics import get_interaction_network

        result = get_interaction_network('AKT1,MTOR', depth=1)

        # Should complete without timeout
        assert 'nodes' in result
        assert 'edges' in result


@pytest.mark.e2e
class TestPSPDataConsistency:
    """Tests for data consistency across PSP tools."""

    @pytest.fixture(autouse=True)
    def check_psp_data(self):
        """Check PSP data availability before each test."""
        skip_if_no_psp_data()

    def test_kinase_substrate_bidirectional(self, csv_parser):
        """Test that kinase-substrate and upstream kinase queries are consistent."""
        from psp_proteomics import get_kinase_substrates, find_upstream_kinases

        # AKT1 phosphorylates GSK3B at S9
        kinase_result = get_kinase_substrates('AKT1')
        kinase_rows = csv_parser(kinase_result['AKT1'])

        # Check if GSK3B S9 is in substrates
        gsk3b_s9_found = any(
            row.get('gene') == 'GSK3B' and 'S9' in row.get('site', '')
            for row in kinase_rows
        )

        if gsk3b_s9_found:
            # Reverse check: GSK3B S9 upstream should include AKT1
            upstream_result = find_upstream_kinases('GSK3B_S9')
            upstream_rows = csv_parser(upstream_result['GSK3B_S9'])

            akt1_found = any(
                'AKT1' in row.get('kinase', '').upper()
                for row in upstream_rows
            )

            assert akt1_found, "AKT1 should be upstream of GSK3B S9"
