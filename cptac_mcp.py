"""
CPTAC MCP Server - MCP interface for CPTAC proteomics data

This module provides MCP tools for querying CPTAC proteomics datasets.
Data loading and preprocessing is handled by cptac_backend.py.
"""

import logging
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import false_discovery_control
from mcp.server.fastmcp import FastMCP

from cptac_backend import CancerManager

# Initialize FastMCP server
mcp = FastMCP("cptac_query")

# Cohort storage - persists within MCP server session
_cohorts: dict = {}


################################
######## Cohort Helpers ########
################################

def _parse_filters(filters: str) -> list:
    """
    Parse filter string into list of (column, operator, value) tuples.

    Supported formats:
    - 'column = value' (equality)
    - 'column in value1,value2,value3' (in-list)
    - 'column >= value', 'column <= value', 'column > value', 'column < value' (numeric)
    - Multiple filters separated by semicolons

    Returns list of tuples: [(column, operator, value), ...]
    """
    parsed = []

    # Split by semicolon for multiple filters
    filter_parts = [f.strip() for f in filters.split(';')]

    for part in filter_parts:
        if not part:
            continue

        # Try to match different operators (order matters - longer operators first)
        matched = False

        # Check for 'in' operator (must be surrounded by spaces)
        if ' in ' in part:
            col, values_str = part.split(' in ', 1)
            col = col.strip()
            # Split values by comma
            values = [v.strip() for v in values_str.split(',')]
            parsed.append((col, 'in', values))
            matched = True
        elif ' >= ' in part:
            col, val = part.split(' >= ', 1)
            parsed.append((col.strip(), '>=', float(val.strip())))
            matched = True
        elif ' <= ' in part:
            col, val = part.split(' <= ', 1)
            parsed.append((col.strip(), '<=', float(val.strip())))
            matched = True
        elif ' > ' in part:
            col, val = part.split(' > ', 1)
            parsed.append((col.strip(), '>', float(val.strip())))
            matched = True
        elif ' < ' in part:
            col, val = part.split(' < ', 1)
            parsed.append((col.strip(), '<', float(val.strip())))
            matched = True
        elif ' = ' in part:
            col, val = part.split(' = ', 1)
            parsed.append((col.strip(), '=', val.strip()))
            matched = True

        if not matched:
            raise ValueError(
                f"Could not parse filter: '{part}'. "
                f"Supported formats: 'column = value', 'column in val1,val2', "
                f"'column >= value', 'column <= value', 'column > value', 'column < value'"
            )

    return parsed


def _apply_filters(clinical_df: pd.DataFrame, parsed_filters: list) -> pd.DataFrame:
    """
    Apply parsed filters to clinical dataframe, return filtered df.

    Parameters:
    - clinical_df: DataFrame with clinical data (patients as rows)
    - parsed_filters: List of (column, operator, value) tuples

    Returns filtered DataFrame (subset of rows matching all filters).
    """
    mask = pd.Series([True] * len(clinical_df), index=clinical_df.index)

    for col, op, val in parsed_filters:
        if col not in clinical_df.columns:
            # Check for duplicate column handling (col_1, col_2, etc.)
            found = False
            for actual_col in clinical_df.columns:
                if actual_col == col or actual_col.startswith(f"{col}_"):
                    col = actual_col
                    found = True
                    break
            if not found:
                raise ValueError(f"Column '{col}' not found in clinical data. Available columns: {list(clinical_df.columns)[:20]}...")

        col_data = clinical_df[col]

        if op == '=':
            mask = mask & (col_data.astype(str) == str(val))
        elif op == 'in':
            # Convert both to strings for comparison
            col_str = col_data.astype(str)
            val_str = [str(v) for v in val]
            mask = mask & col_str.isin(val_str)
        elif op == '>=':
            mask = mask & (pd.to_numeric(col_data, errors='coerce') >= val)
        elif op == '<=':
            mask = mask & (pd.to_numeric(col_data, errors='coerce') <= val)
        elif op == '>':
            mask = mask & (pd.to_numeric(col_data, errors='coerce') > val)
        elif op == '<':
            mask = mask & (pd.to_numeric(col_data, errors='coerce') < val)

    return clinical_df[mask]


def _get_cohort_patient_ids(cohort_name: str) -> tuple:
    """
    Get cancer type and patient IDs for a cohort.

    Returns: (cancer_type: str, patient_ids: list[str])
    Raises ValueError if cohort not found.
    """
    if cohort_name not in _cohorts:
        available = list(_cohorts.keys()) if _cohorts else "none"
        raise ValueError(f"Cohort '{cohort_name}' not found. Available cohorts: {available}")

    cohort = _cohorts[cohort_name]
    return cohort['cancer'], cohort['patient_ids']


def _filter_columns_by_cohort(data_columns: list, patient_ids: list) -> list:
    """
    Filter data columns to only include samples from cohort patients.

    Handles both tumor (e.g., 'C3L-00004') and normal (e.g., 'C3L-00004.N') samples.
    Returns columns where the patient ID (without .N suffix) is in the cohort.
    """
    cohort_set = set(patient_ids)
    filtered = []
    for col in data_columns:
        # Get patient ID by removing .N suffix if present
        patient_id = col.replace('.N', '')
        if patient_id in cohort_set:
            filtered.append(col)
    return filtered


################################
######## MCP Tools #############
################################

@mcp.tool()
def get_cancer_types():
    """Get list of available CPTAC cancer types."""
    return CancerManager.available_cancers()


################################
######## Cohort Tools ##########
################################

@mcp.tool()
def create_cohort(name: str, cancer: str, filters: str) -> dict:
    """
    Create a new patient cohort by filtering clinical data.

    Parameters:
    -----------
    name : str
        Unique cohort identifier (e.g., 'pdac_early_stage', 'brca_female_over50')
    cancer : str
        Cancer type. Supported: 'brca', 'coad', 'hnscc', 'luad', 'ovarian', 'ccrcc', 'gbm', 'lscc', 'pdac'
    filters : str
        Filter expression to select patients. Supported formats:
        - Equality: 'column = value' (e.g., 'sex = Female')
        - In-list: 'column in value1,value2,value3' (e.g., 'tumor_stage_pathological in Stage I,Stage II')
        - Numeric: 'column >= value', 'column <= value', 'column > value', 'column < value'
        - Multiple filters: separate with semicolons (AND logic)
        Example: 'sex = Female; age >= 50'

    Returns:
    --------
    dict
        - name: Cohort name
        - cancer: Cancer type
        - filters: Original filter string
        - n_patients: Number of patients in cohort
        - patient_ids: List of patient IDs (first 20 shown if >20)
    """
    global _cohorts

    # Validate cancer type
    try:
        cancer_obj = CancerManager.get(cancer)
    except ValueError:
        return {"error": f"Unsupported cancer type: {cancer}. Supported types are: {CancerManager.available_cancers()}"}

    # Check if cohort name already exists
    if name in _cohorts:
        return {"error": f"Cohort '{name}' already exists. Use delete_cohort first or choose a different name."}

    # Get clinical data
    clinical = cancer_obj.clinical_data
    if clinical.empty:
        return {"error": f"No clinical data available for {cancer}"}

    # Parse and apply filters
    try:
        parsed_filters = _parse_filters(filters)
    except ValueError as e:
        return {"error": str(e)}

    try:
        filtered_clinical = _apply_filters(clinical, parsed_filters)
    except ValueError as e:
        return {"error": str(e)}

    if len(filtered_clinical) == 0:
        return {
            "error": "No patients match the specified filters",
            "filters": filters,
            "n_total_patients": len(clinical)
        }

    # Extract patient IDs from index
    patient_ids = filtered_clinical.index.tolist()

    # Store cohort
    _cohorts[name] = {
        'name': name,
        'cancer': cancer,
        'filters': filters,
        'patient_ids': patient_ids,
        'n_patients': len(patient_ids)
    }

    logging.info(f"Created cohort '{name}' with {len(patient_ids)} patients from {cancer}")

    # Return result (limit patient_ids display for readability)
    return {
        'name': name,
        'cancer': cancer,
        'filters': filters,
        'n_patients': len(patient_ids),
        'patient_ids': patient_ids[:20] if len(patient_ids) > 20 else patient_ids,
        'patient_ids_truncated': len(patient_ids) > 20
    }


@mcp.tool()
def list_cohorts() -> dict:
    """
    List all created cohorts.

    Returns:
    --------
    dict
        - cohorts: List of cohort summaries (name, cancer, n_patients, filters)
        - n_cohorts: Total number of cohorts
    """
    if not _cohorts:
        return {
            'cohorts': [],
            'n_cohorts': 0,
            'message': 'No cohorts have been created. Use create_cohort to create one.'
        }

    cohort_list = []
    for name, cohort in _cohorts.items():
        cohort_list.append({
            'name': cohort['name'],
            'cancer': cohort['cancer'],
            'n_patients': cohort['n_patients'],
            'filters': cohort['filters']
        })

    return {
        'cohorts': cohort_list,
        'n_cohorts': len(cohort_list)
    }


@mcp.tool()
def get_cohort(name: str) -> dict:
    """
    Get detailed information about a specific cohort.

    Parameters:
    -----------
    name : str
        Name of the cohort to retrieve

    Returns:
    --------
    dict
        - name: Cohort name
        - cancer: Cancer type
        - filters: Original filter string
        - n_patients: Number of patients
        - patient_ids: Full list of patient IDs
    """
    if name not in _cohorts:
        available = list(_cohorts.keys()) if _cohorts else []
        return {
            "error": f"Cohort '{name}' not found",
            "available_cohorts": available
        }

    cohort = _cohorts[name]

    return {
        'name': cohort['name'],
        'cancer': cohort['cancer'],
        'filters': cohort['filters'],
        'n_patients': cohort['n_patients'],
        'patient_ids': cohort['patient_ids']
    }


@mcp.tool()
def delete_cohort(name: str) -> dict:
    """
    Delete a cohort.

    Parameters:
    -----------
    name : str
        Name of the cohort to delete

    Returns:
    --------
    dict
        - deleted: Name of deleted cohort
        - remaining_cohorts: Number of remaining cohorts
    """
    global _cohorts

    if name not in _cohorts:
        available = list(_cohorts.keys()) if _cohorts else []
        return {
            "error": f"Cohort '{name}' not found",
            "available_cohorts": available
        }

    del _cohorts[name]
    logging.info(f"Deleted cohort '{name}'")

    return {
        'deleted': name,
        'remaining_cohorts': len(_cohorts)
    }


################################
######## Analysis Tools ########
################################

@mcp.tool()
def phospho_tumor_vs_normal(cancer, query, normalized="true", cohort=None):
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
    cohort : str, optional
        Name of a previously created cohort to filter analysis to specific patients.
        Must match the cancer type. Use create_cohort to create cohorts.

    Returns:
    --------
    dict
        Dictionary containing:
        - data: Dictionary mapping gene names to CSV strings
                Format: {"GENE1": "site,peptide,database_id,log2_fold_change,p_value,p_value_adjusted,n_pairs,mean_tumor,mean_normal\nS473,...\n..."}
                Each gene's data is a CSV string with phosphosites as rows
        - cohort: Name of cohort used (if specified)
        - cohort_n_patients: Number of patients in cohort (if specified)
    """
    # Get Cancer object from manager
    try:
        cancer_obj = CancerManager.get(cancer)
    except ValueError as e:
        raise ValueError(f"Unsupported cancer type: {cancer}. Supported types are: {CancerManager.available_cancers()}")

    # Validate cohort if specified
    cohort_patient_ids = None
    if cohort is not None:
        try:
            cohort_cancer, cohort_patient_ids = _get_cohort_patient_ids(cohort)
        except ValueError as e:
            return {"error": str(e)}

        if cohort_cancer != cancer:
            return {"error": f"Cohort '{cohort}' is for cancer '{cohort_cancer}', but you requested cancer '{cancer}'"}

        logging.info(f"Using cohort '{cohort}' with {len(cohort_patient_ids)} patients")

    # Get phospho data using cached properties
    if normalized == 'true':
        phospho = cancer_obj.phospho_normalized
    elif normalized == 'false':
        phospho = cancer_obj.phospho_deduplicated
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
        all_cols = list(result_df.columns)

        # Filter to cohort patients if specified
        if cohort_patient_ids is not None:
            all_cols = _filter_columns_by_cohort(all_cols, cohort_patient_ids)

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

    # Convert to CSV format grouped by gene
    # Group results by gene
    gene_groups = {}
    for r in results:
        gene = r.get('gene')
        if gene not in gene_groups:
            gene_groups[gene] = []
        gene_groups[gene].append(r)

    # Create CSV for each gene
    gene_data = {}
    for gene, gene_results in gene_groups.items():
        csv_rows = ["site,peptide,database_id,log2_fold_change,p_value,p_value_adjusted,n_pairs,mean_tumor,mean_normal"]

        for r in gene_results:
            site = r.get('site', '')
            peptide = r.get('peptide', '')
            db_id = r.get('database_id', '')
            log2_fc = f"{r.get('log2_fold_change'):.3f}" if r.get('log2_fold_change') is not None else ""
            pval = f"{r.get('p_value'):.4f}" if r.get('p_value') is not None else ""
            pval_adj = f"{r.get('p_value_adjusted'):.4f}" if r.get('p_value_adjusted') is not None else ""
            n_pairs = r.get('n_pairs', 0)
            mean_tumor = f"{r.get('mean_tumor'):.3f}" if r.get('mean_tumor') is not None else ""
            mean_normal = f"{r.get('mean_normal'):.3f}" if r.get('mean_normal') is not None else ""

            csv_rows.append(f"{site},{peptide},{db_id},{log2_fc},{pval},{pval_adj},{n_pairs},{mean_tumor},{mean_normal}")

        gene_data[gene] = "\n".join(csv_rows)

    result = {'data': gene_data}
    if cohort is not None:
        result['cohort'] = cohort
        result['cohort_n_patients'] = len(cohort_patient_ids)
    return result


@mcp.tool()
def protein_tumor_vs_normal(cancer, query, cohort=None):
    """
    Query proteomics data for specific proteins and calculate tumor vs normal statistics.

    Parameters:
    -----------
    cancer : cptac cancer name (str)
        Supported cancers: 'brca', 'coad', 'hnscc', 'luad', 'ovarian', 'ccrcc', 'gbm', 'lscc', 'pdac'
    query : str
        Comma-separated list of gene names (e.g., 'AKT1,TP53,EGFR')
    cohort : str, optional
        Name of a previously created cohort to filter analysis to specific patients.
        Must match the cancer type. Use create_cohort to create cohorts.

    Returns:
    --------
    dict
        Dictionary containing:
        - data: CSV-formatted string with header row
                Format: "gene,log2_fold_change,p_value,p_value_adjusted,n_pairs,mean_tumor,mean_normal"
                Each row contains statistics for one protein
        - cohort: Name of cohort used (if specified)
        - cohort_n_patients: Number of patients in cohort (if specified)
    """
    # Get Cancer object from manager
    try:
        cancer_obj = CancerManager.get(cancer)
    except ValueError as e:
        raise ValueError(f"Unsupported cancer type: {cancer}. Supported types are: {CancerManager.available_cancers()}")

    # Validate cohort if specified
    cohort_patient_ids = None
    if cohort is not None:
        try:
            cohort_cancer, cohort_patient_ids = _get_cohort_patient_ids(cohort)
        except ValueError as e:
            return {"error": str(e)}

        if cohort_cancer != cancer:
            return {"error": f"Cohort '{cohort}' is for cancer '{cohort_cancer}', but you requested cancer '{cancer}'"}

        logging.info(f"Using cohort '{cohort}' with {len(cohort_patient_ids)} patients")

    # Get proteomics data using cached property
    proteomics = cancer_obj.proteomics_deduplicated
    logging.info(f"Proteomics data shape: {proteomics.shape} (proteins x samples)")

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
        all_cols = list(result_df.columns)

        # Filter to cohort patients if specified
        if cohort_patient_ids is not None:
            all_cols = _filter_columns_by_cohort(all_cols, cohort_patient_ids)

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

    # Convert to CSV format
    csv_rows = ["gene,log2_fold_change,p_value,p_value_adjusted,n_pairs,mean_tumor,mean_normal"]

    for r in results:
        gene = r.get('gene')
        log2_fc = f"{r.get('log2_fold_change'):.3f}" if r.get('log2_fold_change') is not None else ""
        pval = f"{r.get('p_value'):.4f}" if r.get('p_value') is not None else ""
        pval_adj = f"{r.get('p_value_adjusted'):.4f}" if r.get('p_value_adjusted') is not None else ""
        n_pairs = r.get('n_pairs', 0)
        mean_tumor = f"{r.get('mean_tumor'):.3f}" if r.get('mean_tumor') is not None else ""
        mean_normal = f"{r.get('mean_normal'):.3f}" if r.get('mean_normal') is not None else ""

        csv_rows.append(f"{gene},{log2_fc},{pval},{pval_adj},{n_pairs},{mean_tumor},{mean_normal}")

    result = {'data': "\n".join(csv_rows)}
    if cohort is not None:
        result['cohort'] = cohort
        result['cohort_n_patients'] = len(cohort_patient_ids)
    return result


@mcp.tool()
def get_clinical_data(cancer: str, columns: str = 'summary') -> dict:
    """
    Get clinical metadata for a CPTAC cancer cohort.

    Parameters:
    -----------
    cancer : str
        CPTAC cancer name. Supported: 'brca', 'coad', 'hnscc', 'luad',
        'ovarian', 'ccrcc', 'gbm', 'lscc', 'pdac'
    columns : str
        Column selection mode (default: 'summary'):
        - 'summary': Returns metadata about each column (type, stats/categories)
        - 'col1,col2,col3': Returns patient-by-patient CSV data for only those columns

    Returns:
    --------
    dict
        If columns='summary':
            - n_samples: Number of patients/samples
            - columns: Dict mapping column names to metadata (type, stats or categories)
        If specific columns requested:
            - data: CSV-formatted string with only the requested columns
            - n_samples: Number of patients/samples
            - columns: List of returned column names
    """
    # Get Cancer object from manager
    try:
        cancer_obj = CancerManager.get(cancer)
    except ValueError as e:
        raise ValueError(f"Unsupported cancer type: {cancer}. Supported types are: {CancerManager.available_cancers()}")

    # Access clinical data via cached property
    clinical = cancer_obj.clinical_data
    logging.info(f"Clinical data shape: {clinical.shape} (samples x attributes)")

    # Handle empty clinical data
    if clinical.empty:
        return {
            'n_samples': 0,
            'columns': {},
            'error': f'No clinical data available for {cancer}'
        }

    # Reset index to include patient IDs as a column
    clinical_reset = clinical.reset_index()

    # Get the index column name (typically 'Patient_ID' or similar)
    index_col_name = clinical_reset.columns[0]

    if columns == 'summary':
        # Summary mode: return column metadata
        column_info = {}

        # Handle duplicate column names by iterating with iloc
        seen_cols = {}
        for i, col in enumerate(clinical_reset.columns):
            # Track duplicate column names with suffixes
            if col in seen_cols:
                seen_cols[col] += 1
                display_col = f"{col}_{seen_cols[col]}"
            else:
                seen_cols[col] = 0
                display_col = col

            # Use iloc to get the specific column by position (handles duplicates)
            col_data = clinical_reset.iloc[:, i]
            non_null_count = int(col_data.notna().sum())

            if pd.api.types.is_numeric_dtype(col_data):
                # Numeric column: compute stats
                column_info[display_col] = {
                    'type': 'numeric',
                    'mean': round(float(col_data.mean()), 2) if non_null_count > 0 else None,
                    'min': float(col_data.min()) if non_null_count > 0 else None,
                    'max': float(col_data.max()) if non_null_count > 0 else None,
                    'non_null': non_null_count
                }
            else:
                # Categorical column: list unique values
                unique_values = col_data.dropna().unique().tolist()
                # Convert to strings and limit display if too many
                unique_str = [str(v) for v in unique_values]
                column_info[display_col] = {
                    'type': 'categorical',
                    'categories': unique_str if len(unique_str) <= 20 else unique_str[:20] + [f'... ({len(unique_str) - 20} more)'],
                    'n_categories': len(unique_values),
                    'non_null': non_null_count
                }

        logging.info(f"Returning summary for {len(column_info)} columns")

        return {
            'n_samples': len(clinical),
            'columns': column_info
        }
    else:
        # Specific columns mode: return CSV with requested columns
        requested_cols = [col.strip() for col in columns.split(',')]

        # Build mapping from display names (with suffixes for duplicates) to column indices
        col_name_to_indices = {}
        seen_cols = {}
        for i, col in enumerate(clinical_reset.columns):
            if col in seen_cols:
                seen_cols[col] += 1
                display_col = f"{col}_{seen_cols[col]}"
            else:
                seen_cols[col] = 0
                display_col = col
            col_name_to_indices[display_col] = i

        # Always include the patient ID column first (index 0)
        indices_to_return = [0] if index_col_name not in requested_cols else []
        cols_to_return = [index_col_name] if index_col_name not in requested_cols else []

        # Validate and add requested columns
        missing_cols = []
        for col in requested_cols:
            if col in col_name_to_indices:
                idx = col_name_to_indices[col]
                if idx not in indices_to_return:
                    indices_to_return.append(idx)
                    cols_to_return.append(col)
            else:
                missing_cols.append(col)

        if missing_cols:
            logging.warning(f"Requested columns not found: {missing_cols}")

        if len(cols_to_return) == 1 and cols_to_return[0] == index_col_name:
            # No valid columns requested (only patient ID would be returned)
            return {
                'error': f"None of the requested columns were found: {requested_cols}",
                'available_columns': list(col_name_to_indices.keys())
            }

        # Subset the dataframe by column indices and convert to CSV
        subset_df = clinical_reset.iloc[:, indices_to_return]
        # Rename columns to use display names (handles duplicates)
        subset_df.columns = cols_to_return
        csv_data = subset_df.to_csv(index=False, quoting=1)

        logging.info(f"Returning {len(cols_to_return)} columns for {len(clinical)} samples")

        result = {
            'data': csv_data,
            'n_samples': len(clinical),
            'columns': cols_to_return
        }

        if missing_cols:
            result['missing_columns'] = missing_cols

        return result


@mcp.tool()
def correlation_analysis(cancer, query, data_type="phospho", normalized="true", cohort=None):
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
    cohort : str, optional
        Name of a previously created cohort to filter analysis to specific patients.
        Must match the cancer type. Use create_cohort to create cohorts.

    Returns:
    --------
    dict
        Dictionary containing:
        - correlation_matrix: CSV-formatted correlation matrix with row/column labels
                             Format: ",item1,item2,...\nitem1,1.000,0.456,...\nitem2,0.456,1.000,..."
        - p_value_matrix: CSV-formatted p-value matrix (FDR-corrected)
                         Format: ",item1,item2,...\nitem1,0.000,0.001,...\nitem2,0.001,0.000,..."
        - n_samples: Number of tumor samples used in correlation analysis
        - cohort: Name of cohort used (if specified)
        - cohort_n_patients: Number of patients in cohort (if specified)
        OR if matrix would be too large (>25 items):
        - available_sites: List of available sites that can be queried
        - message: Instructions to refine query
    """
    # Get Cancer object from manager
    try:
        cancer_obj = CancerManager.get(cancer)
    except ValueError as e:
        raise ValueError(f"Unsupported cancer type: {cancer}. Supported types are: {CancerManager.available_cancers()}")

    # Validate cohort if specified
    cohort_patient_ids = None
    if cohort is not None:
        try:
            cohort_cancer, cohort_patient_ids = _get_cohort_patient_ids(cohort)
        except ValueError as e:
            return {"error": str(e)}

        if cohort_cancer != cancer:
            return {"error": f"Cohort '{cohort}' is for cancer '{cohort_cancer}', but you requested cancer '{cancer}'"}

        logging.info(f"Using cohort '{cohort}' with {len(cohort_patient_ids)} patients")

    # Validate data_type parameter
    if data_type not in ['phospho', 'proteomics', 'both']:
        raise ValueError("Parameter 'data_type' must be 'phospho', 'proteomics', or 'both'")

    # Get data based on data_type using cached properties
    phospho_data = None
    protein_data = None

    if data_type in ['phospho', 'both']:
        if normalized == 'true':
            phospho_data = cancer_obj.phospho_normalized
        elif normalized == 'false':
            phospho_data = cancer_obj.phospho_deduplicated
        else:
            raise ValueError("Parameter 'normalized' must be 'true' or 'false'")

    if data_type in ['proteomics', 'both']:
        protein_data = cancer_obj.proteomics_deduplicated

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
                    # Only add if not already present (avoid duplicates)
                    if item_label not in item_data:
                        item_labels.append(item_label)
                        item_data[item_label] = phospho_data.loc[idx]
                        found = True
                        logging.info(f"Found phosphosite {item_label}")
                    else:
                        found = True
                        logging.info(f"Phosphosite {item_label} already in query, skipping duplicate")
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
                    added_count = 0
                    for idx in matching_rows:
                        item_label = f"{idx[0]}_{idx[1]}"
                        # Only add if not already present (avoid duplicates)
                        if item_label not in item_data:
                            item_labels.append(item_label)
                            item_data[item_label] = phospho_data.loc[idx]
                            added_count += 1
                    found = True
                    logging.info(f"Found {added_count} phosphosites for gene {query_item}")

        # If not found in phospho, try protein data (or if forced to protein)
        if not found and protein_data is not None:
            gene = query_item.split('_')[0] if '_' in query_item and not force_protein else query_item

            if isinstance(protein_data.index, pd.MultiIndex):
                matching_indices = [idx for idx in protein_data.index if idx[0] == gene]
                if matching_indices:
                    idx = matching_indices[0]
                    item_label = f"{idx[0]}_protein"
                    # Only add if not already present (avoid duplicates)
                    if item_label not in item_data:
                        item_labels.append(item_label)
                        item_data[item_label] = protein_data.loc[idx]
                        found = True
                        logging.info(f"Found protein {item_label}")
                    else:
                        found = True
                        logging.info(f"Protein {item_label} already in query, skipping duplicate")
            else:
                if gene in protein_data.index:
                    item_label = f"{gene}_protein"
                    # Only add if not already present (avoid duplicates)
                    if item_label not in item_data:
                        item_labels.append(item_label)
                        item_data[item_label] = protein_data.loc[gene]
                        found = True
                        logging.info(f"Found protein {item_label}")
                    else:
                        found = True
                        logging.info(f"Protein {item_label} already in query, skipping duplicate")

        if not found:
            logging.warning(f"No data found for query item: {query_item}")

    # Return error if no items found
    if len(item_data) == 0:
        return {"error": f"No matching data found for query: {query}"}

    # Create DataFrame with items as rows and samples as columns
    data_matrix = pd.DataFrame(item_data).T
    logging.info(f"Data matrix shape: {data_matrix.shape} (items x samples)")
    logging.info(f"Number of item_labels: {len(item_labels)}")
    logging.info(f"Number of rows in data_matrix: {len(data_matrix)}")

    # Reset index to ensure integer indexing works properly
    data_matrix.index = range(len(data_matrix))

    # Filter to tumor samples only (exclude .N samples)
    tumor_cols = [col for col in data_matrix.columns if '.N' not in col]

    # Filter to cohort patients if specified
    if cohort_patient_ids is not None:
        tumor_cols = _filter_columns_by_cohort(tumor_cols, cohort_patient_ids)

    data_matrix = data_matrix[tumor_cols]
    logging.info(f"Using {len(tumor_cols)} tumor samples for correlation analysis")
    logging.info(f"Final data_matrix shape: {data_matrix.shape}")

    # Calculate pairwise correlations
    # Use data_matrix length to ensure we don't go out of bounds
    n_items = len(data_matrix)
    logging.info(f"n_items for correlation matrix: {n_items}")

    # Verify item_labels matches data_matrix rows
    if len(item_labels) != n_items:
        logging.error(f"Mismatch: {len(item_labels)} item_labels but {n_items} rows in data_matrix")
        return {"error": f"Internal error: item count mismatch ({len(item_labels)} labels vs {n_items} data rows)"}
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

    # Convert to CSV format
    # Create two matrices: correlation and p-value
    # Header row with item labels
    corr_csv_rows = [",".join([""] + item_labels)]  # Empty cell for top-left, then column headers
    pval_csv_rows = [",".join([""] + item_labels)]

    for i in range(n_items):
        # Row label + correlation values
        corr_row = [item_labels[i]]
        pval_row = [item_labels[i]]

        for j in range(n_items):
            corr_val = corr_matrix[i, j]
            pval_val = adjusted_pval_matrix[i, j]

            if not np.isnan(corr_val):
                corr_row.append(f"{corr_val:.3f}")
            else:
                corr_row.append("")

            if not np.isnan(pval_val):
                pval_row.append(f"{pval_val:.4f}")
            else:
                pval_row.append("")

        corr_csv_rows.append(",".join(corr_row))
        pval_csv_rows.append(",".join(pval_row))

    result = {
        'correlation_matrix': "\n".join(corr_csv_rows),
        'p_value_matrix': "\n".join(pval_csv_rows),
        'n_samples': len(tumor_cols)
    }
    if cohort is not None:
        result['cohort'] = cohort
        result['cohort_n_patients'] = len(cohort_patient_ids)
    return result


################################
######## Test Function #########
################################

def test():
    """Test CPTAC MCP tools"""
    logging.basicConfig(level=logging.INFO)

    # Test 1a: get_clinical_data - summary mode (default)
    print('\n' + '='*60)
    print('Test 1a: get_clinical_data - summary mode (BRCA)')
    print('='*60)
    try:
        result = get_clinical_data(cancer='brca')

        if isinstance(result, dict):
            if 'error' in result and result.get('n_samples', 0) == 0:
                print(f"Error: {result['error']}")
            else:
                print(f"\nSuccess!")
                print(f"Number of samples: {result.get('n_samples', 'N/A')}")
                columns_info = result.get('columns', {})
                print(f"Number of columns: {len(columns_info)}")

                # Show a few column summaries
                print(f"\nColumn summaries (first 5):")
                for i, (col_name, col_info) in enumerate(columns_info.items()):
                    if i >= 5:
                        break
                    if col_info['type'] == 'numeric':
                        print(f"  {col_name}: numeric, mean={col_info.get('mean')}, range=[{col_info.get('min')}, {col_info.get('max')}], n={col_info.get('non_null')}")
                    else:
                        cats = col_info.get('categories', [])
                        cats_str = str(cats[:3]) + '...' if len(cats) > 3 else str(cats)
                        print(f"  {col_name}: categorical, {col_info.get('n_categories')} categories: {cats_str}, n={col_info.get('non_null')}")

                print(f"\nTest 1a PASSED!")
    except Exception as e:
        print(f"Exception occurred: {e}")
        import traceback
        traceback.print_exc()

    # Test 1b: get_clinical_data - specific columns mode
    print('\n' + '='*60)
    print('Test 1b: get_clinical_data - specific columns (BRCA)')
    print('='*60)
    try:
        # First get summary to find available columns
        summary = get_clinical_data(cancer='brca')
        available_cols = list(summary.get('columns', {}).keys())

        # Pick a few columns to request (skip the first which is typically patient ID)
        test_cols = available_cols[1:4] if len(available_cols) > 3 else available_cols[1:]
        test_cols_str = ','.join(test_cols)
        print(f"Requesting columns: {test_cols_str}")

        result = get_clinical_data(cancer='brca', columns=test_cols_str)

        if isinstance(result, dict):
            if 'error' in result:
                print(f"Error: {result['error']}")
            else:
                print(f"\nSuccess!")
                print(f"Number of samples: {result.get('n_samples', 'N/A')}")
                print(f"Returned columns: {result.get('columns', [])}")

                # Show first few rows of data
                if 'data' in result and result['data']:
                    data_lines = result['data'].split('\n')
                    print(f"\nClinical Data (first 5 rows):")
                    for line in data_lines[:6]:
                        print(f"  {line[:150]}..." if len(line) > 150 else f"  {line}")

                print(f"\nTest 1b PASSED!")
    except Exception as e:
        print(f"Exception occurred: {e}")
        import traceback
        traceback.print_exc()

    # Test 2: correlation_analysis
    print('\n' + '='*60)
    print('Test 2: Correlation analysis with mixed phospho sites (COAD)')
    print('='*60)
    try:
        result = correlation_analysis(
            cancer='coad',
            query='ELL_S442,ELL_S309,POLR2A_S2,POLR2A_S5,CDK1_T14,CDK1_T161,CDK2_T160,TP53_S15,TP53_S392,SMAD2_S245,SMAD2_S250,SMAD3_S423,STAT3_S727,E2F1_S364',
            data_type='both'
        )

        if isinstance(result, dict):
            if 'error' in result:
                print(f"Error: {result['error']}")
            elif 'message' in result:
                print(f"Message: {result['message']}")
                if 'available_sites' in result:
                    print(f"Available sites: {result['available_sites'][:5]}...")
            else:
                print(f"\nSuccess!")
                print(f"Number of samples: {result.get('n_samples', 'N/A')}")

                # Parse correlation matrix to show sample
                if 'correlation_matrix' in result:
                    corr_lines = result['correlation_matrix'].split('\n')
                    print(f"\nCorrelation Matrix (first 3 rows):")
                    for line in corr_lines[:4]:
                        print(f"  {line}")

                print(f"\nTest 2 PASSED!")
    except Exception as e:
        print(f"Exception occurred: {e}")
        import traceback
        traceback.print_exc()

    # Test 3: Cohort creation and management
    print('\n' + '='*60)
    print('Test 3: Cohort creation and management (BRCA)')
    print('='*60)
    try:
        # Clear any existing cohorts
        _cohorts.clear()

        # Test 3a: Create cohort with equality filter
        print("\n3a: Creating cohort with equality filter (sex = Female)...")
        result = create_cohort(name='brca_female', cancer='brca', filters='sex = Female')
        if 'error' in result:
            print(f"  Error: {result['error']}")
        else:
            print(f"  Created cohort 'brca_female' with {result['n_patients']} patients")
            print(f"  Sample patient IDs: {result['patient_ids'][:3]}...")

        # Test 3b: Create cohort with numeric filter
        print("\n3b: Creating cohort with numeric filter (Age_at_Diagnosis >= 60)...")
        result = create_cohort(name='brca_older', cancer='brca', filters='Age_at_Diagnosis >= 60')
        if 'error' in result:
            print(f"  Error: {result['error']}")
        else:
            print(f"  Created cohort 'brca_older' with {result['n_patients']} patients")

        # Test 3c: Create cohort with in-list filter
        print("\n3c: Creating cohort with in-list filter...")
        result = create_cohort(name='brca_stages', cancer='brca', filters='tumor_stage_pathological in Stage I,Stage II,Stage IA,Stage IB,Stage IIA,Stage IIB')
        if 'error' in result:
            print(f"  Error: {result['error']}")
        else:
            print(f"  Created cohort 'brca_stages' with {result['n_patients']} patients")

        # Test 3d: Create cohort with multiple filters
        print("\n3d: Creating cohort with multiple filters (sex = Female; Age_at_Diagnosis >= 50)...")
        result = create_cohort(name='brca_female_older', cancer='brca', filters='sex = Female; Age_at_Diagnosis >= 50')
        if 'error' in result:
            print(f"  Error: {result['error']}")
        else:
            print(f"  Created cohort 'brca_female_older' with {result['n_patients']} patients")

        # Test 3e: List cohorts
        print("\n3e: Listing all cohorts...")
        result = list_cohorts()
        print(f"  Total cohorts: {result['n_cohorts']}")
        for c in result['cohorts']:
            print(f"    - {c['name']}: {c['n_patients']} patients ({c['cancer']})")

        # Test 3f: Get specific cohort
        print("\n3f: Getting cohort details for 'brca_female'...")
        result = get_cohort(name='brca_female')
        if 'error' in result:
            print(f"  Error: {result['error']}")
        else:
            print(f"  Cohort: {result['name']}")
            print(f"  Patients: {result['n_patients']}")
            print(f"  Filters: {result['filters']}")

        # Test 3g: Delete cohort
        print("\n3g: Deleting cohort 'brca_stages'...")
        result = delete_cohort(name='brca_stages')
        if 'error' in result:
            print(f"  Error: {result['error']}")
        else:
            print(f"  Deleted: {result['deleted']}")
            print(f"  Remaining cohorts: {result['remaining_cohorts']}")

        print(f"\nTest 3 PASSED!")
    except Exception as e:
        print(f"Exception occurred: {e}")
        import traceback
        traceback.print_exc()

    # Test 4: Analysis with cohorts
    print('\n' + '='*60)
    print('Test 4: Analysis with cohort filtering (BRCA)')
    print('='*60)
    try:
        # Ensure we have a cohort
        if 'brca_female' not in _cohorts:
            create_cohort(name='brca_female', cancer='brca', filters='sex = Female')

        # Test 4a: phospho_tumor_vs_normal with cohort
        print("\n4a: phospho_tumor_vs_normal with cohort...")
        result = phospho_tumor_vs_normal(cancer='brca', query='AKT1_S473', cohort='brca_female')
        if 'error' in result:
            print(f"  Error: {result['error']}")
        else:
            print(f"  Cohort used: {result.get('cohort', 'N/A')}")
            print(f"  Cohort patients: {result.get('cohort_n_patients', 'N/A')}")
            if 'data' in result:
                for gene, csv_data in result['data'].items():
                    lines = csv_data.split('\n')
                    print(f"  {gene}: {len(lines)-1} phosphosites")
                    if len(lines) > 1:
                        print(f"    First row: {lines[1][:80]}...")

        # Test 4b: protein_tumor_vs_normal with cohort
        print("\n4b: protein_tumor_vs_normal with cohort...")
        result = protein_tumor_vs_normal(cancer='brca', query='AKT1,TP53', cohort='brca_female')
        if 'error' in result:
            print(f"  Error: {result['error']}")
        else:
            print(f"  Cohort used: {result.get('cohort', 'N/A')}")
            print(f"  Cohort patients: {result.get('cohort_n_patients', 'N/A')}")
            if 'data' in result:
                lines = result['data'].split('\n')
                print(f"  Results: {len(lines)-1} proteins")

        # Test 4c: correlation_analysis with cohort
        print("\n4c: correlation_analysis with cohort...")
        result = correlation_analysis(cancer='brca', query='AKT1_S473,TP53_S15', cohort='brca_female', data_type='phospho')
        if 'error' in result:
            print(f"  Error: {result['error']}")
        elif 'message' in result:
            print(f"  Message: {result['message']}")
        else:
            print(f"  Cohort used: {result.get('cohort', 'N/A')}")
            print(f"  Cohort patients: {result.get('cohort_n_patients', 'N/A')}")
            print(f"  Samples used: {result.get('n_samples', 'N/A')}")

        # Test 4d: Error case - wrong cancer type for cohort
        print("\n4d: Testing error case (cohort cancer mismatch)...")
        result = phospho_tumor_vs_normal(cancer='coad', query='AKT1', cohort='brca_female')
        if 'error' in result:
            print(f"  Expected error: {result['error']}")
        else:
            print(f"  ERROR: Should have returned an error!")

        print(f"\nTest 4 PASSED!")
    except Exception as e:
        print(f"Exception occurred: {e}")
        import traceback
        traceback.print_exc()

    # Clean up cohorts at the end
    _cohorts.clear()
    print('\n' + '='*60)
    print('All tests completed!')
    print('='*60)


################################
######## Main Entry Point ######
################################

def main():
    logging.info("Starting CPTAC Query Tool")
    mcp.run(transport='stdio')


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == 'test':
        test()
    else:
        main()
