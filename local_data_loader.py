"""
Local Data Loader for CPTAC Proteomics Data

Loads CPTAC data directly from local files instead of using the cptac library.
This avoids Zenodo download issues in Docker environments.

Author: Claude Code
Date: 2026-01-13
"""

import os
import gzip
import logging
import pandas as pd
import numpy as np
from typing import Dict, Optional, Tuple

logging.basicConfig(level=logging.INFO)


class LocalDataLoader:
    """
    Loads CPTAC proteomics data from local files.

    Data is expected in the structure:
        data/bcm-{cancer}/{CANCER}_phospho_site_abundance_log2_reference_intensity_normalized_{Tumor|Normal}.txt
        data/bcm-{cancer}/{CANCER}_proteomics_gene_abundance_log2_reference_intensity_normalized_{Tumor|Normal}.txt.gz

    Attributes:
        data_dir: Path to the data directory
        ensp_to_gene: Mapping from ENSEMBL protein ID to gene symbol
        ensg_to_gene: Mapping from ENSEMBL gene ID to gene symbol
    """

    # Map cancer names to folder/file naming conventions
    CANCER_MAP = {
        'brca': ('bcm-brca', 'BRCA'),
        'coad': ('bcm-coad', 'COAD'),
        'ccrcc': ('bcm-ccrcc', 'CCRCC'),
        'gbm': ('bcm-gbm', 'GBM'),
        'hnscc': ('bcm-hnscc', 'HNSCC'),
        'lscc': ('bcm-lscc', 'LSCC'),
        'luad': ('bcm-luad', 'LUAD'),
        'ovarian': ('bcm-ov', 'OV'),
        'pdac': ('bcm-pdac', 'PDAC'),
    }

    def __init__(self, data_dir: str = None):
        """
        Initialize the loader with gene mappings.

        Args:
            data_dir: Path to data directory. Defaults to ./data relative to this file.
        """
        if data_dir is None:
            data_dir = os.path.join(os.path.dirname(__file__), 'data')

        self.data_dir = data_dir
        logging.info(f"LocalDataLoader initialized with data_dir: {data_dir}")

        # Load gene mappings
        self.ensp_to_gene = {}
        self.ensg_to_gene = {}
        self._load_gene_mappings()

    def _load_gene_mappings(self):
        """Load gene symbol mappings from cptac_genes.csv and gencode annotations."""
        # Load ENSP -> Gene mapping from cptac_genes.csv
        genes_file = os.path.join(self.data_dir, 'cptac_genes.csv')
        if os.path.exists(genes_file):
            genes_df = pd.read_csv(genes_file)
            self.ensp_to_gene = dict(zip(genes_df['Database_ID'], genes_df['Gene_Name']))
            logging.info(f"Loaded {len(self.ensp_to_gene)} ENSP->Gene mappings")
        else:
            logging.warning(f"Gene mapping file not found: {genes_file}")

        # Load ENSG -> Gene mapping from gencode annotation (first available)
        # Format: transcript, protein, gene, gene_name, coding, MANE, SwissProt
        for cancer_folder, _ in self.CANCER_MAP.values():
            gencode_file = os.path.join(self.data_dir, cancer_folder, 'gencode.v34.basic.annotation-mapping.txt.gz')
            if os.path.exists(gencode_file):
                try:
                    with gzip.open(gencode_file, 'rt') as f:
                        header = f.readline()  # Skip header
                        for line in f:
                            parts = line.strip().split('\t')
                            if len(parts) >= 4:
                                ensg_id = parts[2]    # gene column (ENSG...)
                                gene_name = parts[3]  # gene_name column
                                self.ensg_to_gene[ensg_id] = gene_name
                    logging.info(f"Loaded {len(self.ensg_to_gene)} ENSG->Gene mappings from {gencode_file}")
                    break
                except Exception as e:
                    logging.warning(f"Failed to load gencode mapping: {e}")

    def _parse_phospho_index(self, idx: str) -> Tuple[str, str, str, str]:
        """
        Parse phospho row index into (Gene_Name, Site, Peptide, Database_ID).

        Input format: 'ENSG...|ENSP...|Site|Peptide|1'
        Example: 'ENSG00000048028.11|ENSP00000003302.4|S1053|PPTIRPNSPYDLCSR|1'

        Returns:
            Tuple of (Gene_Name, Site, Peptide, Database_ID)
        """
        parts = idx.split('|')
        if len(parts) < 4:
            return (idx, '', '', idx)

        ensg_id = parts[0]
        ensp_id = parts[1]
        site = parts[2]
        peptide = parts[3]

        # Try ENSP mapping first, then ENSG
        gene_name = self.ensp_to_gene.get(ensp_id)
        if gene_name is None:
            gene_name = self.ensg_to_gene.get(ensg_id)
        if gene_name is None:
            # Fallback: use ENSP ID without version
            gene_name = ensp_id.split('.')[0] if '.' in ensp_id else ensp_id

        return (gene_name, site, peptide, ensp_id)

    def _parse_proteomics_index(self, idx: str) -> str:
        """
        Parse proteomics row index into Gene_Name.

        Input format: 'ENSG00000000003.15'

        Returns:
            Gene symbol
        """
        gene_name = self.ensg_to_gene.get(idx)
        if gene_name is None:
            # Try without version number
            base_id = idx.split('.')[0] if '.' in idx else idx
            gene_name = self.ensg_to_gene.get(base_id, idx)
        return gene_name

    def _load_file(self, filepath: str) -> Optional[pd.DataFrame]:
        """Load a data file (handles both .txt and .txt.gz)."""
        if not os.path.exists(filepath):
            # Try gzipped version
            if not filepath.endswith('.gz'):
                filepath_gz = filepath + '.gz'
                if os.path.exists(filepath_gz):
                    filepath = filepath_gz
                else:
                    return None
            else:
                return None

        try:
            if filepath.endswith('.gz'):
                df = pd.read_csv(filepath, sep='\t', index_col=0, compression='gzip')
            else:
                df = pd.read_csv(filepath, sep='\t', index_col=0)
            return df
        except Exception as e:
            logging.error(f"Error loading {filepath}: {e}")
            return None

    def get_phosphoproteomics(self, cancer: str, source: str = 'bcm') -> pd.DataFrame:
        """
        Load phosphoproteomics data for a cancer type.

        Args:
            cancer: Cancer name (e.g., 'brca', 'luad')
            source: Data source (default 'bcm', currently only bcm supported)

        Returns:
            DataFrame with samples as rows, phosphosites as columns (transposed format).
            Columns have MultiIndex: (Name, Site, Peptide, Database_ID)
        """
        if cancer not in self.CANCER_MAP:
            raise ValueError(f"Unknown cancer type: {cancer}. Available: {list(self.CANCER_MAP.keys())}")

        folder, prefix = self.CANCER_MAP[cancer]
        base_path = os.path.join(self.data_dir, folder)

        # Load tumor data
        tumor_file = os.path.join(base_path, f"{prefix}_phospho_site_abundance_log2_reference_intensity_normalized_Tumor.txt")
        tumor_df = self._load_file(tumor_file)

        if tumor_df is None:
            raise FileNotFoundError(f"Tumor phospho file not found for {cancer}")

        logging.info(f"Loaded tumor phospho data: {tumor_df.shape}")

        # Load normal data if available
        normal_file = os.path.join(base_path, f"{prefix}_phospho_site_abundance_log2_reference_intensity_normalized_Normal.txt")
        normal_df = self._load_file(normal_file)

        if normal_df is not None:
            logging.info(f"Loaded normal phospho data: {normal_df.shape}")
            # Add .N suffix to normal sample columns
            normal_df.columns = [f"{col}.N" for col in normal_df.columns]
            # Merge tumor and normal
            df = pd.concat([tumor_df, normal_df], axis=1)
        else:
            logging.info(f"No normal phospho data available for {cancer}")
            df = tumor_df

        # Parse index into MultiIndex
        parsed_indices = [self._parse_phospho_index(idx) for idx in df.index]
        df.index = pd.MultiIndex.from_tuples(parsed_indices, names=['Name', 'Site', 'Peptide', 'Database_ID'])

        # Transpose to match cptac library format: samples as columns
        df = df.T

        logging.info(f"Final phospho DataFrame: {df.shape} (samples x phosphosites)")
        return df

    def get_proteomics(self, cancer: str, source: str = 'bcm') -> pd.DataFrame:
        """
        Load proteomics data for a cancer type.

        Args:
            cancer: Cancer name (e.g., 'brca', 'luad')
            source: Data source (default 'bcm', currently only bcm supported)

        Returns:
            DataFrame with samples as rows, proteins as columns (transposed format).
            Columns are gene names.
        """
        if cancer not in self.CANCER_MAP:
            raise ValueError(f"Unknown cancer type: {cancer}. Available: {list(self.CANCER_MAP.keys())}")

        folder, prefix = self.CANCER_MAP[cancer]
        base_path = os.path.join(self.data_dir, folder)

        # Load tumor data (gzipped)
        tumor_file = os.path.join(base_path, f"{prefix}_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt.gz")
        tumor_df = self._load_file(tumor_file)

        if tumor_df is None:
            raise FileNotFoundError(f"Tumor proteomics file not found for {cancer}")

        logging.info(f"Loaded tumor proteomics data: {tumor_df.shape}")

        # Load normal data if available
        normal_file = os.path.join(base_path, f"{prefix}_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt.gz")
        normal_df = self._load_file(normal_file)

        if normal_df is not None:
            logging.info(f"Loaded normal proteomics data: {normal_df.shape}")
            # Add .N suffix to normal sample columns
            normal_df.columns = [f"{col}.N" for col in normal_df.columns]
            # Merge tumor and normal
            df = pd.concat([tumor_df, normal_df], axis=1)
        else:
            logging.info(f"No normal proteomics data available for {cancer}")
            df = tumor_df

        # Map ENSG IDs to gene names
        new_index = [self._parse_proteomics_index(idx) for idx in df.index]
        df.index = new_index

        # Transpose to match cptac library format: samples as columns
        df = df.T

        logging.info(f"Final proteomics DataFrame: {df.shape} (samples x proteins)")
        return df

    def get_cancer_types(self) -> list:
        """Return list of available cancer types."""
        return list(self.CANCER_MAP.keys())


# Module-level loader instance (lazy initialization)
_loader = None


def get_loader(data_dir: str = None) -> LocalDataLoader:
    """Get or create the module-level loader instance."""
    global _loader
    if _loader is None:
        _loader = LocalDataLoader(data_dir)
    return _loader


def test():
    """Test the LocalDataLoader."""
    print("Testing LocalDataLoader...")
    print("=" * 60)

    loader = LocalDataLoader()

    print(f"\nAvailable cancer types: {loader.get_cancer_types()}")
    print(f"ENSP mappings loaded: {len(loader.ensp_to_gene)}")
    print(f"ENSG mappings loaded: {len(loader.ensg_to_gene)}")

    # Test with PDAC (has both tumor and normal)
    print("\n" + "=" * 60)
    print("Testing PDAC phosphoproteomics...")
    print("=" * 60)
    try:
        phospho = loader.get_phosphoproteomics('pdac')
        print(f"Shape: {phospho.shape}")
        print(f"Index type: {type(phospho.index)}")
        print(f"Columns type: {type(phospho.columns)}")
        print(f"\nSample columns (first 5): {list(phospho.index[:5])}")
        print(f"Normal samples: {sum('.N' in str(col) for col in phospho.index)}")
        print(f"\nSample phosphosite columns:")
        for col in list(phospho.columns[:3]):
            print(f"  {col}")
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()

    print("\n" + "=" * 60)
    print("Testing PDAC proteomics...")
    print("=" * 60)
    try:
        prot = loader.get_proteomics('pdac')
        print(f"Shape: {prot.shape}")
        print(f"\nSample columns (first 5): {list(prot.index[:5])}")
        print(f"Normal samples: {sum('.N' in str(col) for col in prot.index)}")
        print(f"\nSample gene columns (first 10): {list(prot.columns[:10])}")
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()

    # Test with BRCA (tumor only)
    print("\n" + "=" * 60)
    print("Testing BRCA phosphoproteomics (tumor only)...")
    print("=" * 60)
    try:
        phospho = loader.get_phosphoproteomics('brca')
        print(f"Shape: {phospho.shape}")
        print(f"Normal samples: {sum('.N' in str(col) for col in phospho.index)}")
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()

    print("\n" + "=" * 60)
    print("Test complete!")
    print("=" * 60)


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == 'test':
        test()
    else:
        print("Usage: python local_data_loader.py test")
