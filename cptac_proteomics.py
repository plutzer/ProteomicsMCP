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


def get_deduplicated_proteomics(cancer):
    """
    Get deduplicated proteomics data in transposed format.

    Parameters:
    -----------
    cancer : cptac cancer object
        CPTAC cancer dataset object

    Returns:
    --------
    pd.DataFrame
        Transposed, deduplicated proteomics data with samples as columns
    """
    logging.info("Loading and deduplicating proteomics data...")
    proteomics = cancer.get_proteomics('bcm').T
    proteomics.columns = proteomics.columns.values.tolist()

    # Check if index is MultiIndex and deduplicate if needed
    if isinstance(proteomics.index, pd.MultiIndex):
        proteomics = proteomics.groupby(level=list(range(proteomics.index.nlevels))).agg('mean').replace(-np.inf, np.nan).replace(np.inf, np.nan)
        logging.info(f"Deduplicated proteomics data")

    logging.info(f"Proteomics data shape: {proteomics.shape} (proteins x samples)")
    return proteomics

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

    # Convert to tabular format for better LLM retrieval accuracy
    columns = ['gene', 'site', 'peptide', 'database_id', 'log2_fold_change',
               'p_value', 'p_value_adjusted', 'n_pairs', 'mean_tumor', 'mean_normal']
    data = []
    for r in results:
        row = [
            r.get('gene'),
            r.get('site'),
            r.get('peptide'),
            r.get('database_id'),
            r.get('log2_fold_change'),
            r.get('p_value'),
            r.get('p_value_adjusted'),
            r.get('n_pairs'),
            r.get('mean_tumor'),
            r.get('mean_normal')
        ]
        data.append(row)

    return {'columns': columns, 'data': data}


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

    # Convert to tabular format for better LLM retrieval accuracy
    columns = ['gene', 'log2_fold_change', 'p_value', 'p_value_adjusted',
               'n_pairs', 'mean_tumor', 'mean_normal']
    data = []
    for r in results:
        row = [
            r.get('gene'),
            r.get('log2_fold_change'),
            r.get('p_value'),
            r.get('p_value_adjusted'),
            r.get('n_pairs'),
            r.get('mean_tumor'),
            r.get('mean_normal')
        ]
        data.append(row)

    return {'columns': columns, 'data': data}


@mcp.tool()
def correlation_analysis(cancer, query, data_type="phospho", normalized="true"):
    """
    Analyze correlations between items in a query set across tumor samples.

    Parameters:
    -----------
    cancer : cptac cancer name (str)
        Supported cancers: 'brca', 'coad', 'hnscc', 'luad', 'ovarian', 'ccrcc', 'gbm', 'lscc', 'pdac'
    query : str
        Comma-separated list of items to correlate.
        For phospho data: 'Gene_Site' format (e.g., 'AKT1_S473,TP53_S15') or gene names
        For protein data: gene names (e.g., 'AKT1,TP53')
        To force protein-only search: append '_protein' to gene name (e.g., 'TSC2_protein')
        Can mix phospho and protein queries if data_type='both'
        Examples: 'MTOR_S2448,TSC2_protein' queries MTOR S2448 phosphosite and TSC2 protein
    data_type : str
        Type of data to analyze: 'phospho', 'proteomics', or 'both' (default: 'phospho')
    normalized : str
        For phospho data, if 'true' use protein-normalized data, if 'false' use raw data (default: 'true')

    Returns:
    --------
    dict
        Dictionary containing:
        - correlation_matrix: Pairwise Pearson correlations between all queried items
        - p_value_matrix: Adjusted p-values (FDR-corrected) for each correlation
        - item_labels: Labels for each item in the matrices
        - n_samples: Number of samples used in correlation analysis
        OR if matrix would be too large (>50 items):
        - available_sites: List of available sites that can be queried
        - message: Instructions to refine query
    """

    cancer_obj = cancers.get(cancer)
    if cancer_obj is None:
        raise ValueError(f"Unsupported cancer type: {cancer}. Supported types are: {list(cancers.keys())}")

    # Validate data_type parameter
    if data_type not in ['phospho', 'proteomics', 'both']:
        raise ValueError("Parameter 'data_type' must be 'phospho', 'proteomics', or 'both'")

    # Get data based on data_type
    phospho_data = None
    protein_data = None

    if data_type in ['phospho', 'both']:
        if normalized == 'true':
            phospho_data = get_normalized_phospho(cancer_obj)
        elif normalized == 'false':
            phospho_data = get_deduplicated_phospho(cancer_obj)
        else:
            raise ValueError("Parameter 'normalized' must be 'true' or 'false'")

    if data_type in ['proteomics', 'both']:
        protein_data = get_deduplicated_proteomics(cancer_obj)

    # Parse query string
    query_items = [item.strip() for item in query.split(',')]
    logging.info(f"Querying {len(query_items)} items: {query_items}")

    # Build data matrix: rows are items, columns are samples
    item_data = {}  # Dict mapping item label to Series of values across samples
    item_labels = []  # List of item labels in order
    MAX_ITEMS = 25  # Maximum items before returning available sites instead

    for query_item in query_items:
        found = False
        force_protein = False

        # Check if user explicitly requested protein data (ends with _protein)
        if query_item.endswith('_protein'):
            force_protein = True
            query_item = query_item[:-8]  # Remove '_protein' suffix
            logging.info(f"Forcing protein-only search for {query_item}")

        # Determine if this looks like a phosphosite query (Gene_Site format)
        is_phosphosite_query = False
        if '_' in query_item and not force_protein:
            parts = query_item.split('_')
            if len(parts) >= 2:
                potential_site = parts[1]
                # Check if it looks like a phosphosite (starts with S/T/Y and has numbers)
                if len(potential_site) > 0 and potential_site[0] in ['S', 'T', 'Y']:
                    is_phosphosite_query = True

        # Try to find in phospho data first (if available and not forced to protein)
        if phospho_data is not None and not force_protein:
            if is_phosphosite_query:
                # Parse as Gene_Site
                parts = query_item.split('_')
                gene = parts[0]
                site = '_'.join(parts[1:])

                # Search for matching phosphosites
                matching_rows = [idx for idx in phospho_data.index
                               if idx[0] == gene and idx[1] == site]

                if matching_rows:
                    # Use first match if multiple
                    idx = matching_rows[0]
                    item_label = f"{idx[0]}_{idx[1]}"
                    item_labels.append(item_label)
                    item_data[item_label] = phospho_data.loc[idx]
                    found = True
                    logging.info(f"Found phosphosite {item_label}")
            else:
                # Try as gene name in phospho data - get ALL sites for the gene
                matching_rows = [idx for idx in phospho_data.index if idx[0] == query_item]
                if matching_rows:
                    # Check if adding all these sites would exceed limit
                    if len(item_labels) + len(matching_rows) > MAX_ITEMS:
                        # Return available sites instead
                        available_sites = [f"{idx[0]}_{idx[1]}" for idx in matching_rows]
                        logging.warning(f"Query would result in {len(item_labels) + len(matching_rows)} items (max {MAX_ITEMS}). Returning available sites for {query_item}.")
                        return {
                            "message": f"Too many phosphorylation sites ({len(matching_rows)}) for gene {query_item}. Query would exceed maximum of {MAX_ITEMS} total items. Please specify individual sites from the list below.",
                            "gene": query_item,
                            "available_sites": available_sites,
                            "n_sites": len(matching_rows),
                            "current_items": len(item_labels),
                            "max_items": MAX_ITEMS
                        }

                    # Add all matching phosphosites for this gene
                    for idx in matching_rows:
                        item_label = f"{idx[0]}_{idx[1]}"
                        item_labels.append(item_label)
                        item_data[item_label] = phospho_data.loc[idx]
                    found = True
                    logging.info(f"Found {len(matching_rows)} phosphosites for gene {query_item}")

        # If not found in phospho, try protein data (or if forced to protein)
        if not found and protein_data is not None:
            gene = query_item.split('_')[0] if '_' in query_item and not force_protein else query_item

            if isinstance(protein_data.index, pd.MultiIndex):
                matching_indices = [idx for idx in protein_data.index if idx[0] == gene]
                if matching_indices:
                    idx = matching_indices[0]
                    item_label = f"{idx[0]}_protein"
                    item_labels.append(item_label)
                    item_data[item_label] = protein_data.loc[idx]
                    found = True
                    logging.info(f"Found protein {item_label}")
            else:
                if gene in protein_data.index:
                    item_label = f"{gene}_protein"
                    item_labels.append(item_label)
                    item_data[item_label] = protein_data.loc[gene]
                    found = True
                    logging.info(f"Found protein {item_label}")

        if not found:
            logging.warning(f"No data found for query item: {query_item}")

    # Return error if no items found
    if len(item_data) == 0:
        return {"error": f"No matching data found for query: {query}"}

    # Create DataFrame with items as rows and samples as columns
    data_matrix = pd.DataFrame(item_data).T
    logging.info(f"Data matrix shape: {data_matrix.shape} (items x samples)")

    # Filter to tumor samples only (exclude .N samples)
    tumor_cols = [col for col in data_matrix.columns if '.N' not in col]
    data_matrix = data_matrix[tumor_cols]
    logging.info(f"Using {len(tumor_cols)} tumor samples for correlation analysis")

    # Calculate pairwise correlations
    n_items = len(item_labels)
    corr_matrix = np.zeros((n_items, n_items))
    pval_matrix = np.zeros((n_items, n_items))

    for i in range(n_items):
        for j in range(n_items):
            if i == j:
                corr_matrix[i, j] = 1.0
                pval_matrix[i, j] = 0.0
            else:
                # Get data for both items
                x = data_matrix.iloc[i].values
                y = data_matrix.iloc[j].values

                # Remove NaN pairs
                mask = ~(np.isnan(x) | np.isnan(y))
                x_clean = x[mask]
                y_clean = y[mask]

                if len(x_clean) > 2:
                    # Calculate Pearson correlation
                    corr, pval = stats.pearsonr(x_clean, y_clean)
                    corr_matrix[i, j] = corr
                    pval_matrix[i, j] = pval
                else:
                    corr_matrix[i, j] = np.nan
                    pval_matrix[i, j] = np.nan
                    logging.warning(f"Insufficient data for correlation between {item_labels[i]} and {item_labels[j]}")

    # Apply FDR correction to p-values
    # Extract upper triangle p-values (excluding diagonal)
    upper_tri_indices = np.triu_indices(n_items, k=1)
    upper_tri_pvals = pval_matrix[upper_tri_indices]

    # Filter out NaN p-values
    valid_mask = ~np.isnan(upper_tri_pvals)
    valid_pvals = upper_tri_pvals[valid_mask]

    if len(valid_pvals) > 0:
        # Apply FDR correction
        adjusted_pvals = false_discovery_control(valid_pvals, method='bh')

        # Create adjusted p-value matrix
        adjusted_pval_matrix = np.zeros((n_items, n_items))
        adjusted_pval_matrix[:] = np.nan

        # Fill in adjusted p-values
        valid_idx = 0
        for idx_pos, (i, j) in enumerate(zip(upper_tri_indices[0], upper_tri_indices[1])):
            if valid_mask[idx_pos]:
                adjusted_pval_matrix[i, j] = adjusted_pvals[valid_idx]
                adjusted_pval_matrix[j, i] = adjusted_pvals[valid_idx]  # Make symmetric
                valid_idx += 1
            else:
                adjusted_pval_matrix[i, j] = np.nan
                adjusted_pval_matrix[j, i] = np.nan

        # Diagonal should be 0
        for i in range(n_items):
            adjusted_pval_matrix[i, i] = 0.0

        logging.info(f"Applied FDR correction to {len(valid_pvals)} p-values")
    else:
        adjusted_pval_matrix = pval_matrix.copy()
        logging.warning("No valid p-values to adjust")

    # Convert to tabular format for better LLM retrieval accuracy
    # Use lower triangular matrix only (since correlation matrix is symmetric)
    correlations = []
    for i in range(n_items):
        for j in range(i + 1, n_items):  # Only upper triangle (excluding diagonal)
            corr_val = corr_matrix[i, j]
            pval_val = adjusted_pval_matrix[i, j]
            correlations.append([
                item_labels[i],
                item_labels[j],
                float(corr_val) if not np.isnan(corr_val) else None,
                float(pval_val) if not np.isnan(pval_val) else None
            ])

    return {
        'columns': ['item1', 'item2', 'correlation', 'p_value_adjusted'],
        'data': correlations,
        'n_samples': len(tumor_cols)
    }


def test():
    print('Testing correlation_analysis...')
    logging.basicConfig(level=logging.INFO)

    # Test 1: Phospho correlation analysis
    print('\n' + '='*60)
    print('Test 1: Phospho correlation analysis (COAD)')
    print('='*60)
    try:
        result = correlation_analysis('coad', 'CTNNB1_S675,GSK3B_S9,AKT1_S473', data_type='phospho', normalized='true')
        if isinstance(result, dict):
            if 'error' in result:
                print(f"Error: {result['error']}")
            else:
                print(f"\nItems analyzed: {result['item_labels']}")
                print(f"Number of samples: {result['n_samples']}")
                print(f"\nCorrelation Matrix:")
                for item_i in result['item_labels']:
                    row = []
                    for item_j in result['item_labels']:
                        val = result['correlation_matrix'][item_i][item_j]
                        if val is not None:
                            row.append(f"{val:7.3f}")
                        else:
                            row.append("    NaN")
                    print(f"  {item_i:20s}: {' '.join(row)}")

                print(f"\nAdjusted P-value Matrix:")
                for item_i in result['item_labels']:
                    row = []
                    for item_j in result['item_labels']:
                        val = result['p_value_adjusted_matrix'][item_i][item_j]
                        if val is not None:
                            if val == 0.0:
                                row.append("  0.000")
                            else:
                                row.append(f"{val:7.3e}")
                        else:
                            row.append("    NaN")
                    print(f"  {item_i:20s}: {' '.join(row)}")
    except Exception as e:
        print(f"Exception occurred: {e}")
        import traceback
        traceback.print_exc()

    # Test 2: Protein correlation analysis
    print('\n' + '='*60)
    print('Test 2: Protein correlation analysis (LUAD)')
    print('='*60)
    try:
        result = correlation_analysis('luad', 'AKT1,TP53,EGFR,KRAS', data_type='proteomics')
        if isinstance(result, dict):
            if 'error' in result:
                print(f"Error: {result['error']}")
            else:
                print(f"\nItems analyzed: {result['item_labels']}")
                print(f"Number of samples: {result['n_samples']}")
                print(f"\nCorrelation Matrix:")
                for item_i in result['item_labels']:
                    row = []
                    for item_j in result['item_labels']:
                        val = result['correlation_matrix'][item_i][item_j]
                        if val is not None:
                            row.append(f"{val:7.3f}")
                        else:
                            row.append("    NaN")
                    print(f"  {item_i:20s}: {' '.join(row)}")

                # Only print significant correlations
                print(f"\nSignificant correlations (adjusted p < 0.05):")
                for item_i in result['item_labels']:
                    for item_j in result['item_labels']:
                        if item_i < item_j:  # Only print upper triangle
                            corr = result['correlation_matrix'][item_i][item_j]
                            pval = result['p_value_adjusted_matrix'][item_i][item_j]
                            if pval is not None and pval < 0.05:
                                print(f"  {item_i} vs {item_j}: r={corr:.3f}, p_adj={pval:.3e}")
    except Exception as e:
        print(f"Exception occurred: {e}")
        import traceback
        traceback.print_exc()

    # Test 3: TSC2 protein vs MTOR phosphosites (Mixed phospho/protein)
    print('\n' + '='*60)
    print('Test 3: MTOR phosphosites vs TSC2 protein correlation (LUAD)')
    print('='*60)
    try:
        result = correlation_analysis('luad', 'MTOR_S2448,MTOR_S2478,MTOR_S2481,MTOR_S1821,TSC2_protein',
                                     data_type='both', normalized='true')
        if isinstance(result, dict):
            if 'error' in result:
                print(f"Error: {result['error']}")
            else:
                print(f"\nItems analyzed: {result['item_labels']}")
                print(f"Number of samples: {result['n_samples']}")
                print(f"\nCorrelation Matrix:")
                for item_i in result['item_labels']:
                    row = []
                    for item_j in result['item_labels']:
                        val = result['correlation_matrix'][item_i][item_j]
                        if val is not None:
                            row.append(f"{val:7.3f}")
                        else:
                            row.append("    NaN")
                    print(f"  {item_i:20s}: {' '.join(row)}")

                # Highlight TSC2 correlations
                print(f"\nMTOR phosphosite correlations with TSC2 protein:")
                tsc2_label = [l for l in result['item_labels'] if 'TSC2' in l][0]
                for item in result['item_labels']:
                    if item != tsc2_label and 'MTOR' in item:
                        corr = result['correlation_matrix'][item][tsc2_label]
                        pval = result['p_value_adjusted_matrix'][item][tsc2_label]
                        sig = '***' if pval and pval < 0.001 else '**' if pval and pval < 0.01 else '*' if pval and pval < 0.05 else ''
                        if corr is not None and pval is not None:
                            direction = 'NEGATIVE' if corr < 0 else 'positive'
                            print(f"  {item:15s}: r={corr:7.3f}, p_adj={pval:8.4f}  ({direction}) {sig}")
    except Exception as e:
        print(f"Exception occurred: {e}")
        import traceback
        traceback.print_exc()

def main():
    logging.info("Starting CPTAC Query Tool")
    mcp.run(transport='stdio')

if __name__ == "__main__":
    main()
