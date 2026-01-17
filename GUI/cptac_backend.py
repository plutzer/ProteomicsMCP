"""
CPTAC Data Backend

Preprocesses and caches tumor vs normal statistics for all CPTAC cancer types at startup.
Provides instant query access to preprocessed data for the GUI.

Author: Claude Code
Date: 2026-01-11
"""

import pandas as pd
import numpy as np
from typing import List, Dict, Tuple, Optional
from scipy import stats
from scipy.stats import false_discovery_control
import sys
import os
from datetime import datetime

# Import helper functions (using local data loader)
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from local_data_loader import LocalDataLoader
from cptac_proteomics import (
    get_deduplicated_phospho,
    get_normalized_phospho,
    get_deduplicated_proteomics
)

print("[BACKEND_INIT] Importing CPTAC backend modules...")


class CancerDataBackend:
    """
    Preprocesses and caches tumor vs normal statistics for all CPTAC cancer types.

    This class loads all cancer data at initialization and precomputes tumor vs normal
    statistics for every phosphosite and protein. This one-time preprocessing cost
    (5-10 minutes) enables instant queries (<100ms) during app usage.

    Attributes:
        test_mode (bool): If True, only load PDAC for fast testing (~1 min)
        cancer_objects (dict): Maps cancer name to cptac cancer object
        phospho_data (dict): Maps cancer -> {'raw': DataFrame, 'normalized': DataFrame}
        protein_data (dict): Maps cancer -> DataFrame
        phospho_choices (dict): Maps cancer -> list of "GENE_SITE" strings
        protein_choices (dict): Maps cancer -> list of gene names
    """

    def __init__(self, test_mode: bool = False):
        """
        Initialize backend and preprocess all cancer data.

        Args:
            test_mode: If True, only load PDAC for fast testing (~1 min startup)
                      If False, load all 9 cancers (~5-10 min startup)
        """
        print(f"\n{'='*70}")
        print("[BACKEND] Initializing CPTAC Data Backend")
        print(f"[BACKEND] Test mode: {test_mode}")
        print(f"[BACKEND] Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"{'='*70}\n")

        self.test_mode = test_mode
        self.phospho_data = {}
        self.protein_data = {}
        self.phospho_choices = {}
        self.protein_choices = {}

        # Initialize local data loader
        self.data_loader = LocalDataLoader()

        # Define cancer types to load
        if test_mode:
            cancer_names = ['pdac']
            print("[BACKEND] TEST MODE: Loading only PDAC for fast development")
        else:
            cancer_names = list(self.data_loader.CANCER_MAP.keys())
            print(f"[BACKEND] PRODUCTION MODE: Loading all {len(cancer_names)} cancers")

        # Store cancer names
        print("\n[BACKEND] Step 1/3: Initializing data loader...")
        self.cancer_names = cancer_names

        print(f"  [OK] Data loader initialized for {len(cancer_names)} cancers\n")

        # Preprocess each cancer
        print("[BACKEND] Step 2/3: Preprocessing phosphoproteomics data...")
        for i, cancer_name in enumerate(cancer_names, 1):
            print(f"\n  [{i}/{len(cancer_names)}] Processing {cancer_name.upper()} phospho data...")
            self.phospho_data[cancer_name] = self._preprocess_phospho(cancer_name)

            # Build choice lists
            raw_df = self.phospho_data[cancer_name]['raw']
            self.phospho_choices[cancer_name] = [
                f"{gene}_{site}"
                for gene, site in raw_df.index
            ]
            print(f"      [OK] {len(self.phospho_choices[cancer_name])} phosphosites available")

        print("\n[BACKEND] Step 3/3: Preprocessing proteomics data...")
        for i, cancer_name in enumerate(cancer_names, 1):
            print(f"\n  [{i}/{len(cancer_names)}] Processing {cancer_name.upper()} protein data...")
            self.protein_data[cancer_name] = self._preprocess_protein(cancer_name)

            # Build choice lists
            self.protein_choices[cancer_name] = self.protein_data[cancer_name].index.tolist()
            print(f"      [OK] {len(self.protein_choices[cancer_name])} proteins available")

        print(f"\n{'='*70}")
        print("[BACKEND] Initialization Complete!")
        print(f"[BACKEND] End time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"[BACKEND] Ready for instant queries (<100ms)")
        print(f"{'='*70}\n")

    def _compute_tumor_vs_normal_stats(
        self,
        data_df: pd.DataFrame,
        idx: Tuple
    ) -> Optional[Dict]:
        """
        Compute tumor vs normal statistics for a single phosphosite/protein.

        Args:
            data_df: DataFrame with samples as columns, features as rows
            idx: Index tuple (gene, site, peptide, db_id) for phospho or gene for protein

        Returns:
            Dict with statistics or None if insufficient data
        """
        # Get all column names
        all_cols = data_df.columns

        # Identify normal samples (have '.N' in column name)
        normal_cols = [col for col in all_cols if '.N' in col]

        # Find paired tumor samples
        paired_tumor = []
        paired_normal = []

        for normal_col in normal_cols:
            # Create corresponding tumor sample name
            tumor_col = normal_col.replace('.N', '')

            if tumor_col in all_cols:
                tumor_val = data_df.loc[idx, tumor_col]
                normal_val = data_df.loc[idx, normal_col]

                # Only include pairs where both values are not NaN
                if pd.notna(tumor_val) and pd.notna(normal_val):
                    paired_tumor.append(tumor_val)
                    paired_normal.append(normal_val)

        # Return None if insufficient paired samples
        if len(paired_tumor) == 0:
            return None

        # Convert to arrays
        paired_tumor = np.array(paired_tumor)
        paired_normal = np.array(paired_normal)

        # Calculate statistics
        mean_tumor = np.mean(paired_tumor)
        mean_normal = np.mean(paired_normal)
        log2_fc = mean_tumor - mean_normal

        # Perform paired t-test if we have more than 1 pair
        if len(paired_tumor) > 1:
            t_stat, p_value = stats.ttest_rel(paired_tumor, paired_normal)
        else:
            p_value = np.nan

        # Build result dict (structure depends on phospho vs protein)
        if len(idx) == 4:  # Phospho: (gene, site, peptide, db_id)
            result = {
                'gene': idx[0],
                'site': idx[1],
                'peptide': idx[2],
                'database_id': idx[3],
                'log2_fold_change': float(log2_fc),
                'p_value': float(p_value) if pd.notna(p_value) else None,
                'n_pairs': len(paired_tumor),
                'mean_tumor': float(mean_tumor),
                'mean_normal': float(mean_normal)
            }
        else:  # Protein: single gene name
            result = {
                'gene': idx if isinstance(idx, str) else idx[0],
                'log2_fold_change': float(log2_fc),
                'p_value': float(p_value) if pd.notna(p_value) else None,
                'n_pairs': len(paired_tumor),
                'mean_tumor': float(mean_tumor),
                'mean_normal': float(mean_normal)
            }

        return result

    def _preprocess_phospho(
        self,
        cancer_name: str
    ) -> Dict[str, pd.DataFrame]:
        """
        Preprocess all phosphoproteomics data for one cancer type.

        Computes tumor vs normal statistics for ~40k phosphosites and applies
        FDR correction across all sites.

        Args:
            cancer_name: Name of cancer (e.g., 'brca', 'luad')

        Returns:
            Dict with keys 'raw' and 'normalized', each containing a DataFrame
            with MultiIndex (gene, site) and columns for all statistics
        """
        start_time = datetime.now()

        # Load both raw and normalized phospho data using helper functions
        print(f"      Loading raw phospho data...")
        raw_phospho = get_deduplicated_phospho(cancer_name)

        print(f"      Loading normalized phospho data...")
        normalized_phospho = get_normalized_phospho(cancer_name)

        # Process both datasets
        results = {}
        for data_type, phospho_data in [('raw', raw_phospho), ('normalized', normalized_phospho)]:
            print(f"      Processing {data_type} data ({len(phospho_data)} sites)...")

            # Compute statistics for each phosphosite
            site_results = []
            for idx in phospho_data.index:
                stats_dict = self._compute_tumor_vs_normal_stats(phospho_data, idx)

                # Only keep sites with valid statistics
                if stats_dict and stats_dict['p_value'] is not None:
                    site_results.append(stats_dict)

            print(f"      Found {len(site_results)} sites with valid statistics")

            if len(site_results) == 0:
                # Return empty DataFrame with correct structure
                results[data_type] = pd.DataFrame(columns=[
                    'peptide', 'database_id', 'log2_fold_change', 'p_value',
                    'p_value_adjusted', 'n_pairs', 'mean_tumor', 'mean_normal'
                ])
                results[data_type].index = pd.MultiIndex.from_tuples([], names=['gene', 'site'])
                continue

            # Convert to DataFrame
            df = pd.DataFrame(site_results)

            # Apply FDR correction
            print(f"      Applying FDR correction...")
            valid_p_mask = df['p_value'].notna()
            adjusted_p = np.full(len(df), np.nan)

            if valid_p_mask.sum() > 0:
                adjusted_p[valid_p_mask] = false_discovery_control(
                    df.loc[valid_p_mask, 'p_value'].values,
                    method='bh'
                )

            df['p_value_adjusted'] = adjusted_p

            # Set MultiIndex (gene, site) for fast lookup
            df = df.set_index(['gene', 'site'])

            results[data_type] = df
            print(f"      [OK] {data_type.capitalize()} data ready: {len(df)} sites")

        elapsed = (datetime.now() - start_time).total_seconds()
        print(f"      Completed in {elapsed:.1f} seconds")

        return results

    def _preprocess_protein(
        self,
        cancer_name: str
    ) -> pd.DataFrame:
        """
        Preprocess all proteomics data for one cancer type.

        Computes tumor vs normal statistics for ~12k proteins and applies
        FDR correction.

        Args:
            cancer_name: Name of cancer

        Returns:
            DataFrame with gene as index and columns for all statistics
        """
        start_time = datetime.now()

        # Load proteomics data using helper function
        print(f"      Loading proteomics data...")
        protein_data = get_deduplicated_proteomics(cancer_name)

        print(f"      Processing {len(protein_data)} proteins...")

        # Compute statistics for each protein
        protein_results = []
        for idx in protein_data.index:
            stats_dict = self._compute_tumor_vs_normal_stats(protein_data, idx)

            # Only keep proteins with valid statistics
            if stats_dict and stats_dict['p_value'] is not None:
                protein_results.append(stats_dict)

        print(f"      Found {len(protein_results)} proteins with valid statistics")

        if len(protein_results) == 0:
            # Return empty DataFrame
            df = pd.DataFrame(columns=[
                'log2_fold_change', 'p_value', 'p_value_adjusted',
                'n_pairs', 'mean_tumor', 'mean_normal'
            ])
            df.index.name = 'gene'
            return df

        # Convert to DataFrame
        df = pd.DataFrame(protein_results)

        # Apply FDR correction
        print(f"      Applying FDR correction...")
        valid_p_mask = df['p_value'].notna()
        adjusted_p = np.full(len(df), np.nan)

        if valid_p_mask.sum() > 0:
            adjusted_p[valid_p_mask] = false_discovery_control(
                df.loc[valid_p_mask, 'p_value'].values,
                method='bh'
            )

        df['p_value_adjusted'] = adjusted_p

        # Set gene as index
        df = df.set_index('gene')

        elapsed = (datetime.now() - start_time).total_seconds()
        print(f"      Completed in {elapsed:.1f} seconds")

        return df

    # ==================== Public Query Methods ====================

    def get_cancer_types(self) -> List[str]:
        """Return list of available cancer types."""
        return self.cancer_names

    def get_phospho_choices(self, cancer: str) -> List[str]:
        """
        Get list of available phosphosites for selectize input.

        Args:
            cancer: Cancer type (e.g., 'brca', 'luad')

        Returns:
            List of phosphosite labels ["GENE_SITE", ...]
        """
        if cancer not in self.phospho_choices:
            raise ValueError(f"Cancer '{cancer}' not loaded. Available: {self.get_cancer_types()}")

        return self.phospho_choices[cancer]

    def get_protein_choices(self, cancer: str) -> List[str]:
        """
        Get list of available proteins for selectize input.

        Args:
            cancer: Cancer type

        Returns:
            List of gene names ["GENE1", "GENE2", ...]
        """
        if cancer not in self.protein_choices:
            raise ValueError(f"Cancer '{cancer}' not loaded. Available: {self.get_cancer_types()}")

        return self.protein_choices[cancer]

    def get_phospho_background(
        self,
        cancer: str,
        normalized: bool = True
    ) -> pd.DataFrame:
        """
        Get full preprocessed phospho data for volcano plot background.

        Args:
            cancer: Cancer type
            normalized: If True, use protein-normalized data

        Returns:
            DataFrame with all phosphosites (for background layer)
        """
        if cancer not in self.phospho_data:
            raise ValueError(f"Cancer '{cancer}' not loaded")

        data_type = 'normalized' if normalized else 'raw'
        df = self.phospho_data[cancer][data_type].copy()

        # Reset index to have gene and site as columns for easier plotting
        df = df.reset_index()

        return df

    def get_protein_background(self, cancer: str) -> pd.DataFrame:
        """
        Get full preprocessed protein data for volcano plot background.

        Args:
            cancer: Cancer type

        Returns:
            DataFrame with all proteins (for background layer)
        """
        if cancer not in self.protein_data:
            raise ValueError(f"Cancer '{cancer}' not loaded")

        df = self.protein_data[cancer].copy()

        # Reset index to have gene as column
        df = df.reset_index()

        return df

    def query_phospho(
        self,
        cancer: str,
        sites: List[str],
        normalized: bool = True
    ) -> pd.DataFrame:
        """
        Query specific phosphosites.

        Args:
            cancer: Cancer type
            sites: List of site labels ["GENE_SITE", ...]
            normalized: If True, use protein-normalized data

        Returns:
            DataFrame with selected phosphosites (MultiIndex preserved)
        """
        if cancer not in self.phospho_data:
            raise ValueError(f"Cancer '{cancer}' not loaded")

        data_type = 'normalized' if normalized else 'raw'
        df = self.phospho_data[cancer][data_type]

        # Parse site labels into (gene, site) tuples
        query_indices = []
        for site_label in sites:
            if '_' in site_label:
                parts = site_label.split('_', 1)  # Split on first underscore only
                gene = parts[0]
                site = parts[1]
                query_indices.append((gene, site))

        # Filter to requested sites
        # Use .loc with list of tuples for MultiIndex lookup
        found_indices = [idx for idx in query_indices if idx in df.index]

        if not found_indices:
            # Return empty DataFrame with correct structure
            return df.iloc[0:0].copy()

        result_df = df.loc[found_indices].copy()

        return result_df

    def query_protein(
        self,
        cancer: str,
        genes: List[str]
    ) -> pd.DataFrame:
        """
        Query specific proteins.

        Args:
            cancer: Cancer type
            genes: List of gene names ["GENE1", "GENE2", ...]

        Returns:
            DataFrame with selected proteins
        """
        if cancer not in self.protein_data:
            raise ValueError(f"Cancer '{cancer}' not loaded")

        df = self.protein_data[cancer]

        # Filter to requested genes
        found_genes = [gene for gene in genes if gene in df.index]

        if not found_genes:
            # Return empty DataFrame with correct structure
            return df.iloc[0:0].copy()

        result_df = df.loc[found_genes].copy()

        return result_df


def test_backend():
    """Test function for backend development."""
    print("\n" + "="*70)
    print("TESTING BACKEND")
    print("="*70 + "\n")

    # Test with PDAC only
    backend = CancerDataBackend(test_mode=True)

    # Test phospho queries
    print("\n[TEST] Testing phospho queries...")
    cancer = 'pdac'
    choices = backend.get_phospho_choices(cancer)
    print(f"  Available phosphosites: {len(choices)}")
    print(f"  First 5: {choices[:5]}")

    # Query a few sites
    test_sites = choices[:3] if len(choices) >= 3 else choices
    print(f"\n  Querying sites: {test_sites}")
    result = backend.query_phospho(cancer, test_sites, normalized=True)
    print(f"  Result shape: {result.shape}")
    print(f"  Columns: {result.columns.tolist()}")
    print(f"\n  Sample data:")
    print(result.head())

    # Test protein queries
    print("\n[TEST] Testing protein queries...")
    choices = backend.get_protein_choices(cancer)
    print(f"  Available proteins: {len(choices)}")
    print(f"  First 5: {choices[:5]}")

    test_genes = choices[:3] if len(choices) >= 3 else choices
    print(f"\n  Querying genes: {test_genes}")
    result = backend.query_protein(cancer, test_genes)
    print(f"  Result shape: {result.shape}")
    print(f"  Columns: {result.columns.tolist()}")
    print(f"\n  Sample data:")
    print(result.head())

    print("\n" + "="*70)
    print("TEST COMPLETE")
    print("="*70 + "\n")


if __name__ == "__main__":
    # Run test when module is executed directly
    test_backend()
