import cptac
import pandas as pd
import numpy as np
from typing import List, Dict, Tuple, Optional, Union

# Correct logging for mcp
import logging
from typing import Any
from mcp.server.fastmcp import FastMCP

# Initialize FastMCP server
mcp = FastMCP("cptac_query")

cancers = {'brca': cptac.Brca(),
            'coad': cptac.Coad(),
            'hnscc': cptac.Hnscc(),
            'luad': cptac.Luad(),
            'ovarian': cptac.Ov(),
            'ccrcc': cptac.Ccrcc(),
            'gbm': cptac.Gbm(),
            'lscc': cptac.Lscc(),
            'pdac': cptac.Pdac(),
}


def get_deduplicated_phospho(cancer):
    """
    Get deduplicated phosphoproteomics data in transposed format.

    Parameters:
    -----------
    cancer : cptac cancer object
        CPTAC cancer dataset object

    Returns:
    --------
    pd.DataFrame
        Transposed, deduplicated phosphoproteomics data with samples as columns
    """
    logging.info("Loading and deduplicating phosphoproteomics data...")
    phospho = cancer.get_phosphoproteomics('bcm').T
    phospho.columns = phospho.columns.values.tolist()
    phospho = phospho.groupby(level=[0,1,2,3]).agg('mean').replace(-np.inf, np.nan).replace(np.inf, np.nan)
    logging.info(f"Phospho data shape: {phospho.shape} (phosphosites x samples)")
    return phospho


def get_normalized_phospho(cancer):
    """
    Get protein-normalized phosphoproteomics data.

    Parameters:
    -----------
    cancer : cptac cancer object
        CPTAC cancer dataset object

    Returns:
    --------
    pd.DataFrame
        Protein-normalized phosphoproteomics data
    """
    logging.info("Calculating protein-normalized phosphoproteomics...")
    phospho = get_deduplicated_phospho(cancer)
    whole_cell = cancer.get_proteomics('bcm').T
    whole_cell.columns = whole_cell.columns.values.tolist()

    normalized_phospho = phospho.copy()
    samples_normalized = 0
    samples_skipped = 0

    for colname in normalized_phospho.columns:
        if colname in whole_cell.columns:
            wc_col = whole_cell[colname]
            phospho_col = phospho[colname]
            merge = pd.merge(phospho_col, wc_col, left_index=True, right_index=True, suffixes=('_phospho', '_wc'))
            merge[colname] = merge[colname+'_phospho'] - merge[colname+'_wc']
            normalized_phospho[colname] = merge[colname]
            samples_normalized += 1
        else:
            normalized_phospho = normalized_phospho.drop(columns=[colname])
            samples_skipped += 1

    logging.info(f"Normalization complete: {samples_normalized} samples normalized, {samples_skipped} samples removed")
    logging.info(f"Normalized data shape: {normalized_phospho.shape}")
    return normalized_phospho

@mcp.tool()
def get_cancer_types():
    return cptac.get_cancer_options()

@mcp.tool()
def phospho_tumor_vs_normal(cancer, query, normalized):
    """
    Query phosphoproteomics data for specific sites and return raw values.

    Parameters:
    -----------
    cancer : cptac cancer name (str)
        Supported cancers: 'brca', 'coad', 'hnscc', 'luad', 'ovarian', 'ccrcc', 'gbm', 'lscc', 'pdac'
    query : str
        Comma-separated list of phosphosites in the format 'Gene_Site' (e.g., 'AKT1_S473,TP53_S15')
        or gene names (e.g., 'AKT1') to get all sites for that gene
    normalized : bool
        If True, use protein-normalized phospho data. If False, use raw phospho data (default: True)

    Returns:
    --------
    dict
        Dictionary containing raw phosphoproteomics data for queried sites
    """

    cancer_obj = cancers.get(cancer)
    if cancer_obj is None:
        raise ValueError(f"Unsupported cancer type: {cancer}. Supported types are: {list(cancers.keys())}")

    # Get phospho data using helper functions
    if normalized == 'true':
        phospho = get_normalized_phospho(cancer_obj)
    elif normalized == 'false':
        phospho = get_deduplicated_phospho(cancer_obj)
    else:
        raise ValueError("Parameter 'normalized' must be 'true' or 'false'")

    # Parse query string into list of phosphosites
    query_sites = [site.strip() for site in query.split(',')]
    logging.info(f"Querying {len(query_sites)} sites: {query_sites}")

    # Extract gene and site information from query
    query_results = []
    for query_site in query_sites:
        if '_' in query_site:
            parts = query_site.split('_')
            gene = parts[0]
            site = '_'.join(parts[1:])  # Handle cases like GENE_S123_S456

            # Search for matching phosphosites in the index
            # The phospho dataframe has multi-index rows: (Name, Site, Peptide, Database_ID)
            matching_rows = [idx for idx in phospho.index
                           if idx[0] == gene and idx[1] == site]

            if matching_rows:
                logging.info(f"Found {len(matching_rows)} matches for {query_site}")
                for idx in matching_rows:
                    query_results.append(idx)
            else:
                logging.warning(f"No matches found for {query_site}")
        else:
            # If no underscore, treat as gene name and get all sites for that gene
            matching_rows = [idx for idx in phospho.index if idx[0] == query_site]
            if matching_rows:
                logging.info(f"Found {len(matching_rows)} phosphosites for gene {query_site}")
                query_results.extend(matching_rows)
            else:
                logging.warning(f"No phosphosites found for gene {query_site}")

    # Return the filtered data for query sites
    if query_results:
        result_df = phospho.loc[query_results]
        logging.info(f"Returning data for {len(query_results)} phosphosites across {result_df.shape[1]} samples")
        return result_df.to_dict()
    else:
        return {"error": f"No matching phosphosites found for query: {query}"}



def test():
    print('Testing for debugging purposes')
    result = phospho_tumor_vs_normal('brca', 'A2M_S710,AAAS_S2526')
    print(f"\nResult type: {type(result)}")
    if isinstance(result, dict):
        if 'error' in result:
            print(f"Error: {result['error']}")
        else:
            print(f"Number of columns returned: {len(result)}")
            print(f"Columns: {list(result.keys())[:3]}...")  # Show first 3 columns

def main():
    logging.info("Starting CPTAC Query Tool")
    mcp.run(transport='stdio')

if __name__ == "__main__":
    main()
