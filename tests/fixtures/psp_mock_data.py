"""
PSP Mock Data Generators

This module provides functions to generate mock PhosphoSitePlus-style DataFrames
for testing purposes.
"""
import pandas as pd
import numpy as np
from typing import List, Optional


def create_kinase_substrate_df(
    kinases: Optional[List[str]] = None,
    substrates_per_kinase: int = 3,
    organism: str = 'human',
    include_mouse: bool = True,
) -> pd.DataFrame:
    """
    Create a mock kinase-substrate DataFrame matching PSP format.

    Column names match actual PSP Kinase_Substrate_Dataset.txt format.

    Parameters:
    -----------
    kinases : list of str, optional
        Kinase names. Default: ['AKT1', 'MTOR', 'GSK3B', 'CDK1', 'MAPK1']
    substrates_per_kinase : int
        Number of substrates per kinase. Default: 3
    organism : str
        Default organism. Default: 'human'
    include_mouse : bool
        Include some mouse entries for filtering tests. Default: True

    Returns:
    --------
    pd.DataFrame
        Mock kinase-substrate DataFrame
    """
    if kinases is None:
        kinases = ['AKT1', 'MTOR', 'GSK3B', 'CDK1', 'MAPK1']

    # Predefined substrate mappings for realistic data
    substrate_map = {
        'AKT1': [
            ('GSK3B', 'S9', 'RARTSSFAEpSGKGAPM', 'X', '', '12345'),
            ('TSC2', 'S939', 'RKRRVSGGpSPARNTAR', 'X', 'X', '12346'),
            ('FOXO1', 'T24', 'PRRRAAApTPAAAAFAR', '', 'X', '12347'),
        ],
        'MTOR': [
            ('RPS6KB1', 'T389', 'TFLGFTYpTATEPKDGS', 'X', 'X', '23456'),
            ('EIF4EBP1', 'T37', 'STTPGGTpTRPRATES', 'X', '', '23457'),
            ('ULK1', 'S757', 'NLKKGLSpSGLSRGET', '', 'X', '23458'),
        ],
        'GSK3B': [
            ('CTNNB1', 'S33', 'DRKAAVSpSSLNS', 'X', 'X', '34567'),
            ('SNAI1', 'S96', 'PSGKQSpSSPGSSVR', 'X', '', '34568'),
            ('MYC', 'T58', 'PLLKKTpSFFPLSSS', 'X', 'X', '34569'),
        ],
        'CDK1': [
            ('CDKN1A', 'T145', 'RRLIFSpSKRKPEPR', 'X', '', '45678'),
            ('LMNA', 'S22', 'GGSGAQpSPQRRGPL', 'X', 'X', '45679'),
            ('RB1', 'S807', 'PLKSPLpSPRRGTPL', '', 'X', '45680'),
        ],
        'MAPK1': [
            ('RSK1', 'T573', 'KITPPDpTPSPEFAA', 'X', 'X', '56789'),
            ('ELK1', 'S383', 'QKGKPRpSPAPDNSS', 'X', '', '56790'),
            ('STAT3', 'S727', 'RPMSPEpSNNSAPLT', '', 'X', '56791'),
        ],
    }

    rows = []
    for kinase in kinases:
        substrates = substrate_map.get(kinase, [
            ('SUBSTRATE1', 'S100', 'AAAAAApSAAAAAAA', 'X', '', '99999'),
            ('SUBSTRATE2', 'T200', 'AAAAAApTAAAAAAA', '', 'X', '99998'),
            ('SUBSTRATE3', 'Y300', 'AAAAAApYAAAAAAA', 'X', 'X', '99997'),
        ])[:substrates_per_kinase]

        for sub_gene, site, sequence, in_vivo, in_vitro, pmid in substrates:
            rows.append({
                'KINASE': kinase,
                'KIN_ACC_ID': f'P{hash(kinase) % 100000:05d}',
                'KIN_GENE_ID': kinase,
                'KIN_ORGANISM': organism,
                'SUBSTRATE': f'{sub_gene} ({organism})',
                'SUB_GENE': sub_gene,
                'SUB_ACC_ID': f'P{hash(sub_gene) % 100000:05d}',
                'SUB_GENE_ID': sub_gene,
                'SUB_ORGANISM': organism,
                'SUB_MOD_RSD': site,
                'SITE_+/-7_AA': sequence,
                'IN_VIVO_RXN': in_vivo,
                'IN_VITRO_RXN': in_vitro,
                'PMIDs': pmid,
            })

    # Add mouse entries if requested
    if include_mouse:
        rows.append({
            'KINASE': 'Akt1',
            'KIN_ACC_ID': 'Q00001',
            'KIN_GENE_ID': 'Akt1',
            'KIN_ORGANISM': 'mouse',
            'SUBSTRATE': 'Gsk3b (mouse)',
            'SUB_GENE': 'Gsk3b',
            'SUB_ACC_ID': 'Q00002',
            'SUB_GENE_ID': 'Gsk3b',
            'SUB_ORGANISM': 'mouse',
            'SUB_MOD_RSD': 'S9',
            'SITE_+/-7_AA': 'RARTSSFAEpSGKGAPM',
            'IN_VIVO_RXN': 'X',
            'IN_VITRO_RXN': '',
            'PMIDs': '88888',
        })

    return pd.DataFrame(rows)


def create_regulatory_sites_df(
    genes: Optional[List[str]] = None,
    organism: str = 'human',
) -> pd.DataFrame:
    """
    Create a mock regulatory sites DataFrame matching PSP format.

    Parameters:
    -----------
    genes : list of str, optional
        Gene names. Default: ['AKT1', 'TP53', 'MTOR', 'GSK3B']
    organism : str
        Organism name. Default: 'human'

    Returns:
    --------
    pd.DataFrame
        Mock regulatory sites DataFrame
    """
    if genes is None:
        genes = ['AKT1', 'TP53', 'MTOR', 'GSK3B']

    # Predefined regulatory information
    regulatory_map = {
        'AKT1': [
            ('S473-p', 'enzymatic activity, induced', 'cell survival; apoptosis', 'MTOR', '11111;11112', 50, 200),
            ('T308-p', 'enzymatic activity, induced', 'cell growth', 'PDPK1', '11113', 30, 150),
        ],
        'TP53': [
            ('S15-p', 'protein stabilization', 'apoptosis; DNA damage response', '', '22221;22222;22223', 100, 300),
            ('S392-p', 'DNA binding, induced', 'transcription', '', '22224', 20, 50),
            ('S46-p', 'apoptosis, induced', 'cell death', '', '22225;22226', 40, 100),
        ],
        'MTOR': [
            ('S2448-p', 'enzymatic activity, induced', 'cell growth; autophagy', 'AKT1', '33331;33332', 80, 250),
            ('S2481-p', 'enzymatic activity', 'protein synthesis', '', '33333', 15, 80),
        ],
        'GSK3B': [
            ('S9-p', 'enzymatic activity, inhibited', 'glycogen metabolism; cell survival', 'AKT1', '44441', 60, 180),
            ('Y216-p', 'enzymatic activity, induced', 'cell proliferation', '', '44442;44443', 25, 90),
        ],
    }

    rows = []
    for gene in genes:
        sites = regulatory_map.get(gene, [
            ('S100-p', 'unknown', 'unknown', '', '99999', 1, 5),
        ])

        for site, on_function, on_process, on_interact, pmids, lt_lit, ms_lit in sites:
            rows.append({
                'GENE': gene,
                'PROTEIN': f'{gene} ({organism})',
                'ACC_ID': f'P{hash(gene) % 100000:05d}',
                'ORGANISM': organism,
                'MOD_RSD': site,
                'DOMAIN': 'Kinase domain' if 'kinase' in gene.lower() else '',
                'ON_FUNCTION': on_function,
                'ON_PROCESS': on_process,
                'ON_PROT_INTERACT': on_interact,
                'PMIDs': pmids,
                'LT_LIT': lt_lit,
                'MS_LIT': ms_lit,
                'NOTES': '',
            })

    return pd.DataFrame(rows)


def create_disease_sites_df(
    genes: Optional[List[str]] = None,
    organism: str = 'human',
) -> pd.DataFrame:
    """
    Create a mock disease-associated sites DataFrame matching PSP format.

    Parameters:
    -----------
    genes : list of str, optional
        Gene names. Default: ['TP53', 'CTNNB1', 'BRAF', 'EGFR']
    organism : str
        Organism name. Default: 'human'

    Returns:
    --------
    pd.DataFrame
        Mock disease-associated sites DataFrame
    """
    if genes is None:
        genes = ['TP53', 'CTNNB1', 'BRAF', 'EGFR']

    # Predefined disease associations
    disease_map = {
        'TP53': [
            ('S15-p', 'breast cancer', 'increased', '55551', 30, 80, 'DNA damage marker'),
            ('S392-p', 'colorectal cancer', 'altered', '55552', 15, 40, ''),
            ('S46-p', "Alzheimer's disease", 'increased', '55553', 10, 20, 'Neurodegeneration'),
        ],
        'CTNNB1': [
            ('S33-p', 'colorectal cancer', 'decreased', '66661', 50, 120, 'Wnt pathway'),
            ('S45-p', 'hepatocellular carcinoma', 'decreased', '66662', 25, 60, ''),
        ],
        'BRAF': [
            ('V600-p', 'melanoma', 'altered', '77771;77772', 100, 300, 'Driver mutation'),
            ('S446-p', 'thyroid cancer', 'increased', '77773', 20, 50, ''),
        ],
        'EGFR': [
            ('Y1068-p', 'lung cancer', 'increased', '88881', 80, 200, 'Receptor activation'),
            ('Y1173-p', 'glioblastoma', 'increased', '88882', 40, 100, ''),
        ],
    }

    rows = []
    for gene in genes:
        sites = disease_map.get(gene, [
            ('S100-p', 'unknown disease', 'unknown', '99999', 1, 5, ''),
        ])

        for site, disease, alteration, pmids, lt_lit, ms_lit, notes in sites:
            rows.append({
                'GENE': gene,
                'PROTEIN': f'{gene} ({organism})',
                'ACC_ID': f'P{hash(gene) % 100000:05d}',
                'ORGANISM': organism,
                'MOD_RSD': site,
                'DISEASE': disease,
                'ALTERATION': alteration,
                'PMIDs': pmids,
                'LT_LIT': lt_lit,
                'MS_LIT': ms_lit,
                'NOTES': notes,
            })

    return pd.DataFrame(rows)


def create_phospho_sites_df(
    genes: Optional[List[str]] = None,
    organism: str = 'human',
) -> pd.DataFrame:
    """
    Create a mock phosphorylation sites DataFrame matching PSP format.

    This is the main phosphorylation site dataset (Phosphorylation_site_dataset.txt).

    Parameters:
    -----------
    genes : list of str, optional
        Gene names. Default: ['AKT1', 'TP53', 'MTOR', 'GSK3B', 'MAPK1']
    organism : str
        Organism name. Default: 'human'

    Returns:
    --------
    pd.DataFrame
        Mock phosphorylation sites DataFrame
    """
    if genes is None:
        genes = ['AKT1', 'TP53', 'MTOR', 'GSK3B', 'MAPK1']

    # Predefined phosphorylation sites
    site_map = {
        'AKT1': ['S473', 'T308', 'S124', 'T450'],
        'TP53': ['S15', 'S392', 'S46', 'S20', 'T18'],
        'MTOR': ['S2448', 'S2481', 'T2446', 'S1261'],
        'GSK3B': ['S9', 'Y216', 'S389', 'T390'],
        'MAPK1': ['T185', 'Y187', 'T202', 'Y204'],
    }

    rows = []
    for gene in genes:
        sites = site_map.get(gene, ['S100', 'T200', 'Y300'])

        for site in sites:
            mod_type = site[0]  # S, T, or Y
            position = site[1:]

            rows.append({
                'GENE': gene,
                'PROTEIN': f'{gene} ({organism})',
                'ACC_ID': f'P{hash(gene) % 100000:05d}',
                'HU_CHR_LOC': f'chr{hash(gene) % 22 + 1}',
                'ORGANISM': organism,
                'MOD_RSD': f'{site}-p',
                'SITE_GRP_ID': f'{hash((gene, site)) % 1000000}',
                'MW_kD': round(np.random.uniform(30, 150), 1),
                'DOMAIN': '',
                'SITE_+/-7_AA': f'AAAA{mod_type}AAAA',
                'LT_LIT': np.random.randint(1, 100),
                'MS_LIT': np.random.randint(5, 500),
                'MS_CST': np.random.randint(0, 50),
                'CST_CAT#': f'#{np.random.randint(1000, 9999)}',
            })

    return pd.DataFrame(rows)


def create_empty_psp_dataframes():
    """
    Create a set of empty PSP DataFrames with correct column structure.

    Returns:
    --------
    tuple
        (kinase_substrate_df, regulatory_df, disease_df, phospho_df)
    """
    kinase_substrate_cols = [
        'KINASE', 'KIN_ACC_ID', 'KIN_GENE_ID', 'KIN_ORGANISM',
        'SUBSTRATE', 'SUB_GENE', 'SUB_ACC_ID', 'SUB_GENE_ID', 'SUB_ORGANISM',
        'SUB_MOD_RSD', 'SITE_+/-7_AA', 'IN_VIVO_RXN', 'IN_VITRO_RXN', 'PMIDs'
    ]

    regulatory_cols = [
        'GENE', 'PROTEIN', 'ACC_ID', 'ORGANISM', 'MOD_RSD', 'DOMAIN',
        'ON_FUNCTION', 'ON_PROCESS', 'ON_PROT_INTERACT', 'PMIDs', 'LT_LIT', 'MS_LIT', 'NOTES'
    ]

    disease_cols = [
        'GENE', 'PROTEIN', 'ACC_ID', 'ORGANISM', 'MOD_RSD', 'DISEASE',
        'ALTERATION', 'PMIDs', 'LT_LIT', 'MS_LIT', 'NOTES'
    ]

    phospho_cols = [
        'GENE', 'PROTEIN', 'ACC_ID', 'HU_CHR_LOC', 'ORGANISM', 'MOD_RSD',
        'SITE_GRP_ID', 'MW_kD', 'DOMAIN', 'SITE_+/-7_AA', 'LT_LIT', 'MS_LIT', 'MS_CST', 'CST_CAT#'
    ]

    return (
        pd.DataFrame(columns=kinase_substrate_cols),
        pd.DataFrame(columns=regulatory_cols),
        pd.DataFrame(columns=disease_cols),
        pd.DataFrame(columns=phospho_cols),
    )
