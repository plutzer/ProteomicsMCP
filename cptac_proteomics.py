import cptac
import pandas as pd
import numpy as np
from typing import List, Dict, Tuple, Optional, Union
from scipy import stats
from scipy.stats import false_discovery_control


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

################################
######## MCP Tools #############
################################

@mcp.tool()
def get_cancer_types():
    """Get list of available CPTAC cancer types."""
    return cptac.get_cancer_options()

@mcp.tool()
def phospho_tumor_vs_normal(cancer, query, normalized="true"):
    """
    Query phosphoproteomics data for specific sites and calculate tumor vs normal statistics.

    Parameters:
    -----------
    cancer : cptac cancer name (str)
        Supported cancers: 'brca', 'coad', 'hnscc', 'luad', 'ovarian', 'ccrcc', 'gbm', 'lscc', 'pdac'
    query : str
        Comma-separated list of phosphosites in the format 'Gene_Site' (e.g., 'AKT1_S473,TP53_S15')
        or gene names (e.g., 'AKT1') to get all sites for that gene
    normalized : str
        If 'true', use protein-normalized phospho data. If 'false', use raw phospho data (default: 'true')

    Returns:
    --------
    dict
        Dictionary containing for each phosphosite:
        - gene: Gene name
        - site: Phosphorylation site
        - log2_fold_change: Paired log2 fold change (tumor vs normal)
        - p_value: P-value from paired t-test
        - p_value_adjusted: FDR-adjusted p-value (Benjamini-Hochberg correction across all queried sites)
        - n_pairs: Number of paired samples used in analysis
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

    # Return error if no matches found
    if not query_results:
        return {"error": f"No matching phosphosites found for query: {query}"}

    result_df = phospho.loc[query_results]
    logging.info(f"Obtained data for {len(query_results)} phosphosites across {result_df.shape[1]} samples")

    # Calculate paired tumor vs normal statistics
    results = []
    for idx in query_results:
        gene = idx[0]
        site = idx[1]
        peptide = idx[2]
        db_id = idx[3]

        # Get tumor and normal sample columns
        # Normal samples have '.N' in column name
        all_cols = result_df.columns
        normal_cols = [col for col in all_cols if '.N' in col]

        # For each normal sample, find corresponding tumor sample
        paired_tumor = []
        paired_normal = []

        for normal_col in normal_cols:
            # Create tumor sample name by removing '.N'
            tumor_col = normal_col.replace('.N', '')

            if tumor_col in all_cols:
                tumor_val = result_df.loc[idx, tumor_col]
                normal_val = result_df.loc[idx, normal_col]

                # Only include pairs where both values are not NaN
                if pd.notna(tumor_val) and pd.notna(normal_val):
                    paired_tumor.append(tumor_val)
                    paired_normal.append(normal_val)

        # Calculate statistics if we have paired samples
        if len(paired_tumor) > 0:
            paired_tumor = np.array(paired_tumor)
            paired_normal = np.array(paired_normal)

            # Calculate log2 fold change (mean of tumor - mean of normal)
            mean_tumor = np.mean(paired_tumor)
            mean_normal = np.mean(paired_normal)
            log2_fc = mean_tumor - mean_normal

            # Perform paired t-test
            if len(paired_tumor) > 1:
                t_stat, p_value = stats.ttest_rel(paired_tumor, paired_normal)
            else:
                p_value = np.nan

            results.append({
                'gene': gene,
                'site': site,
                'peptide': peptide,
                'database_id': db_id,
                'log2_fold_change': float(log2_fc),
                'p_value': float(p_value) if pd.notna(p_value) else None,
                'n_pairs': len(paired_tumor),
                'mean_tumor': float(mean_tumor),
                'mean_normal': float(mean_normal)
            })
        else:
            logging.warning(f"No paired samples found for {gene}_{site}")
            results.append({
                'gene': gene,
                'site': site,
                'peptide': peptide,
                'database_id': db_id,
                'log2_fold_change': None,
                'p_value': None,
                'n_pairs': 0,
                'mean_tumor': None,
                'mean_normal': None,
                'note': 'No paired samples available'
            })

    logging.info(f"Calculated statistics for {len(results)} phosphosites")

    # Apply multiple testing correction (FDR) to p-values
    p_values = [r['p_value'] for r in results]
    valid_p_indices = [i for i, p in enumerate(p_values) if p is not None and pd.notna(p)]

    if len(valid_p_indices) > 0:
        valid_p_values = [p_values[i] for i in valid_p_indices]

        # Apply Benjamini-Hochberg FDR correction
        adjusted_p_values = false_discovery_control(valid_p_values, method='bh')

        # Update results with adjusted p-values
        for i, idx in enumerate(valid_p_indices):
            results[idx]['p_value_adjusted'] = float(adjusted_p_values[i])

        # Add None for sites without valid p-values
        for i, r in enumerate(results):
            if i not in valid_p_indices:
                r['p_value_adjusted'] = None

        logging.info(f"Applied FDR correction to {len(valid_p_indices)} p-values")
    else:
        # No valid p-values to adjust
        for r in results:
            r['p_value_adjusted'] = None
        logging.warning("No valid p-values to adjust")

    return {'results': results}


@mcp.tool()
def protein_tumor_vs_normal(cancer, query):
    """
    Query proteomics data for specific proteins and calculate tumor vs normal statistics.

    Parameters:
    -----------
    cancer : cptac cancer name (str)
        Supported cancers: 'brca', 'coad', 'hnscc', 'luad', 'ovarian', 'ccrcc', 'gbm', 'lscc', 'pdac'
    query : str
        Comma-separated list of gene names (e.g., 'AKT1,TP53,EGFR')

    Returns:
    --------
    dict
        Dictionary containing for each protein:
        - gene: Gene name
        - log2_fold_change: Paired log2 fold change (tumor vs normal)
        - p_value: P-value from paired t-test
        - p_value_adjusted: FDR-adjusted p-value (Benjamini-Hochberg correction across all queried proteins)
        - n_pairs: Number of paired samples used in analysis
    """

    cancer_obj = cancers.get(cancer)
    if cancer_obj is None:
        raise ValueError(f"Unsupported cancer type: {cancer}. Supported types are: {list(cancers.keys())}")

    # Get proteomics data
    logging.info("Loading proteomics data...")
    proteomics = cancer_obj.get_proteomics('bcm').T
    proteomics.columns = proteomics.columns.values.tolist()
    logging.info(f"Proteomics data shape: {proteomics.shape} (proteins x samples)")

    # Check if index is MultiIndex and deduplicate if needed
    if isinstance(proteomics.index, pd.MultiIndex):
        # Deduplicate by taking mean across all index levels
        proteomics = proteomics.groupby(level=list(range(proteomics.index.nlevels))).agg('mean').replace(-np.inf, np.nan).replace(np.inf, np.nan)
        logging.info(f"Deduplicated proteomics data shape: {proteomics.shape}")

    # Parse query string into list of genes
    query_genes = [gene.strip() for gene in query.split(',')]
    logging.info(f"Querying {len(query_genes)} genes: {query_genes}")

    # Find matching proteins in the data
    query_results = []
    for gene in query_genes:
        if isinstance(proteomics.index, pd.MultiIndex):
            # MultiIndex - search in first level (gene name)
            matching_indices = [idx for idx in proteomics.index if idx[0] == gene]
            if matching_indices:
                query_results.extend(matching_indices)
                logging.info(f"Found {len(matching_indices)} entries for gene {gene}")
            else:
                logging.warning(f"No protein data found for gene {gene}")
        else:
            # Simple index - direct lookup
            if gene in proteomics.index:
                query_results.append(gene)
                logging.info(f"Found protein for gene {gene}")
            else:
                logging.warning(f"No protein data found for gene {gene}")

    # Return error if no matches found
    if not query_results:
        return {"error": f"No matching proteins found for query: {query}"}

    result_df = proteomics.loc[query_results]
    logging.info(f"Obtained data for {len(query_results)} proteins across {result_df.shape[1]} samples")

    # Calculate paired tumor vs normal statistics
    results = []
    for idx in query_results:
        # Extract gene name (handle both MultiIndex and simple index)
        if isinstance(proteomics.index, pd.MultiIndex):
            gene = idx[0]
        else:
            gene = idx
        # Get tumor and normal sample columns
        # Normal samples have '.N' in column name
        all_cols = result_df.columns
        normal_cols = [col for col in all_cols if '.N' in col]

        # For each normal sample, find corresponding tumor sample
        paired_tumor = []
        paired_normal = []

        for normal_col in normal_cols:
            # Create tumor sample name by removing '.N'
            tumor_col = normal_col.replace('.N', '')

            if tumor_col in all_cols:
                tumor_val = result_df.loc[idx, tumor_col]
                normal_val = result_df.loc[idx, normal_col]

                # Only include pairs where both values are not NaN
                if pd.notna(tumor_val) and pd.notna(normal_val):
                    paired_tumor.append(tumor_val)
                    paired_normal.append(normal_val)

        # Calculate statistics if we have paired samples
        if len(paired_tumor) > 0:
            paired_tumor = np.array(paired_tumor)
            paired_normal = np.array(paired_normal)

            # Calculate log2 fold change (mean of tumor - mean of normal)
            mean_tumor = np.mean(paired_tumor)
            mean_normal = np.mean(paired_normal)
            log2_fc = mean_tumor - mean_normal

            # Perform paired t-test
            if len(paired_tumor) > 1:
                t_stat, p_value = stats.ttest_rel(paired_tumor, paired_normal)
            else:
                p_value = np.nan

            results.append({
                'gene': gene,
                'log2_fold_change': float(log2_fc),
                'p_value': float(p_value) if pd.notna(p_value) else None,
                'n_pairs': len(paired_tumor),
                'mean_tumor': float(mean_tumor),
                'mean_normal': float(mean_normal)
            })
        else:
            logging.warning(f"No paired samples found for {gene}")
            results.append({
                'gene': gene,
                'log2_fold_change': None,
                'p_value': None,
                'n_pairs': 0,
                'mean_tumor': None,
                'mean_normal': None,
                'note': 'No paired samples available'
            })

    logging.info(f"Calculated statistics for {len(results)} proteins")

    # Apply multiple testing correction (FDR) to p-values
    p_values = [r['p_value'] for r in results]
    valid_p_indices = [i for i, p in enumerate(p_values) if p is not None and pd.notna(p)]

    if len(valid_p_indices) > 0:
        valid_p_values = [p_values[i] for i in valid_p_indices]

        # Apply Benjamini-Hochberg FDR correction
        adjusted_p_values = false_discovery_control(valid_p_values, method='bh')

        # Update results with adjusted p-values
        for i, idx in enumerate(valid_p_indices):
            results[idx]['p_value_adjusted'] = float(adjusted_p_values[i])

        # Add None for proteins without valid p-values
        for i, r in enumerate(results):
            if i not in valid_p_indices:
                r['p_value_adjusted'] = None

        logging.info(f"Applied FDR correction to {len(valid_p_indices)} p-values")
    else:
        # No valid p-values to adjust
        for r in results:
            r['p_value_adjusted'] = None
        logging.warning("No valid p-values to adjust")

    return {'results': results}


def test():
    print('Testing protein_tumor_vs_normal...')
    logging.basicConfig(level=logging.INFO)
    try:
        result = protein_tumor_vs_normal('ovarian', 'AKT1,TP53')
        print(f"\nResult type: {type(result)}")
        if isinstance(result, dict):
            if 'error' in result:
                print(f"Error: {result['error']}")
            else:
                print(f"Number of results: {len(result.get('results', []))}")
                for r in result.get('results', []):
                    print(f"  {r['gene']}: log2FC={r['log2_fold_change']:.3f}, p={r['p_value']:.3e}, n_pairs={r['n_pairs']}")
    except Exception as e:
        print(f"Exception occurred: {e}")
        import traceback
        traceback.print_exc()

def main():
    logging.info("Starting CPTAC Query Tool")
    mcp.run(transport='stdio')

if __name__ == "__main__":
    main()
