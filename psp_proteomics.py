import pandas as pd
import numpy as np
from typing import Dict, Optional
import os
import logging
from mcp.server.fastmcp import FastMCP

# Initialize FastMCP server
mcp = FastMCP("psp_query")

# Path to PSP datasets
PSP_DIR = os.path.join(os.path.dirname(__file__), 'datasets', 'PSP_2025')

# Load PSP datasets
logging.info("Loading PhosphoSitePlus datasets...")

# Load phosphorylation site dataset
phospho_file = os.path.join(PSP_DIR, 'Phosphorylation_site_dataset.txt')
kinase_substrate_file = os.path.join(PSP_DIR, 'Kinase_Substrate_Dataset.txt')
regulatory_file = os.path.join(PSP_DIR, 'Regulatory_sites.txt')
disease_file = os.path.join(PSP_DIR, 'Disease-associated_sites.txt')

# Global variables to hold loaded data
phospho_data = None
kinase_substrate_data = None
regulatory_data = None
disease_data = None


def load_psp_data():
    """Load all PSP datasets into memory."""
    global phospho_data, kinase_substrate_data, regulatory_data, disease_data

    try:
        # Load phosphorylation site dataset
        phospho_data = pd.read_csv(phospho_file, sep='\t', skiprows=3, encoding='latin1', low_memory=False)
        logging.info(f"Loaded phosphorylation dataset: {phospho_data.shape[0]} rows")

        # Load kinase-substrate dataset
        kinase_substrate_data = pd.read_csv(kinase_substrate_file, sep='\t', skiprows=3, encoding='latin1', low_memory=False)
        logging.info(f"Loaded kinase-substrate dataset: {kinase_substrate_data.shape[0]} rows")

        # Load regulatory sites dataset
        regulatory_data = pd.read_csv(regulatory_file, sep='\t', skiprows=3, encoding='latin1', low_memory=False)
        logging.info(f"Loaded regulatory sites dataset: {regulatory_data.shape[0]} rows")

        # Load disease-associated sites dataset
        disease_data = pd.read_csv(disease_file, sep='\t', skiprows=3, encoding='latin1', low_memory=False)
        logging.info(f"Loaded disease-associated sites dataset: {disease_data.shape[0]} rows")

    except Exception as e:
        logging.error(f"Error loading PSP datasets: {e}")
        raise


# Load data on module import
load_psp_data()


################################
######## MCP Tools #############
################################

@mcp.tool()
def get_kinase_substrates(query: str, organism: str = 'human') -> Dict:
    """
    Get known substrates for specified kinases with phosphorylation sites.

    Parameters:
    -----------
    query : str
        Comma-separated kinase names (e.g., 'AKT1,MTOR,GSK3B')
    organism : str
        Organism name (default: 'human')

    Returns:
    --------
    dict
        Dictionary mapping each kinase to CSV string with substrate information
        Format: {"KINASE1": "substrate,gene,site,site_sequence,in_vivo,in_vitro,pmids\n..."}
    """
    if kinase_substrate_data is None:
        return {"error": "Kinase-substrate data not loaded"}

    # Parse query
    kinases = [k.strip() for k in query.split(',')]
    logging.info(f"Querying {len(kinases)} kinases: {kinases}")

    # Filter by organism
    organism_lower = organism.lower()
    df = kinase_substrate_data[kinase_substrate_data['KIN_ORGANISM'].str.lower() == organism_lower]

    result = {}
    for kinase in kinases:
        # Find matching kinase entries
        kinase_df = df[df['KINASE'].str.upper() == kinase.upper()]

        if len(kinase_df) == 0:
            logging.warning(f"No substrates found for kinase {kinase}")
            result[kinase] = "substrate,gene,site,site_sequence,in_vivo,in_vitro,pmids"
            continue

        # Build CSV
        csv_rows = ["substrate,gene,site,site_sequence,in_vivo,in_vitro,pmids"]

        for _, row in kinase_df.iterrows():
            substrate = str(row.get('SUBSTRATE', '')).replace(',', ';')
            gene = str(row.get('SUB_GENE', ''))
            site = str(row.get('SUB_MOD_RSD', ''))
            sequence = str(row.get('SITE_+/-7_AA', ''))
            in_vivo = 'X' if str(row.get('IN_VIVO_RXN', '')).strip() == 'X' else ''
            in_vitro = 'X' if str(row.get('IN_VITRO_RXN', '')).strip() == 'X' else ''
            pmids = str(row.get('PMIDs', '') if pd.notna(row.get('PMIDs')) else '')

            csv_rows.append(f"{substrate},{gene},{site},{sequence},{in_vivo},{in_vitro},{pmids}")

        result[kinase] = "\n".join(csv_rows)
        logging.info(f"Found {len(kinase_df)} substrates for {kinase}")

    return result


@mcp.tool()
def get_regulatory_sites(query: str, organism: str = 'human') -> Dict:
    """
    Get functional information for phosphorylation sites.

    Parameters:
    -----------
    query : str
        Mixed - gene names OR 'Gene_Site' format (e.g., 'TP53,AKT1_S473,EGFR')
    organism : str
        Organism name (default: 'human')

    Returns:
    --------
    dict
        Dictionary mapping queries to CSV strings with regulatory information
        Format: {"TP53": "site,modification,on_function,on_process,on_interaction,pmids,lt_lit,ms_lit\n..."}
    """
    if regulatory_data is None:
        return {"error": "Regulatory sites data not loaded"}

    # Parse query
    queries = [q.strip() for q in query.split(',')]
    logging.info(f"Querying {len(queries)} items: {queries}")

    # Filter by organism
    organism_lower = organism.lower()
    df = regulatory_data[regulatory_data['ORGANISM'].str.lower() == organism_lower]

    result = {}
    for query_item in queries:
        if '_' in query_item:
            # Parse as Gene_Site
            parts = query_item.split('_')
            gene = parts[0]
            site = '_'.join(parts[1:])

            # Find matching entries
            matches = df[(df['GENE'].str.upper() == gene.upper()) &
                        (df['MOD_RSD'].str.upper() == site.upper())]
        else:
            # Parse as gene name - get all sites
            gene = query_item
            matches = df[df['GENE'].str.upper() == gene.upper()]

        if len(matches) == 0:
            logging.warning(f"No regulatory sites found for {query_item}")
            result[query_item] = "site,modification,on_function,on_process,on_interaction,pmids,lt_lit,ms_lit"
            continue

        # Build CSV
        csv_rows = ["site,modification,on_function,on_process,on_interaction,pmids,lt_lit,ms_lit"]

        for _, row in matches.iterrows():
            site = str(row.get('MOD_RSD', ''))
            mod_type = site.split('-')[-1] if '-' in site else ''
            on_function = str(row.get('ON_FUNCTION', '') if pd.notna(row.get('ON_FUNCTION')) else '').replace(',', ';')
            on_process = str(row.get('ON_PROCESS', '') if pd.notna(row.get('ON_PROCESS')) else '').replace(',', ';')
            on_interact = str(row.get('ON_PROT_INTERACT', '') if pd.notna(row.get('ON_PROT_INTERACT')) else '').replace(',', ';')
            pmids = str(row.get('PMIDs', '') if pd.notna(row.get('PMIDs')) else '').replace(',', ';')
            lt_lit = str(row.get('LT_LIT', '') if pd.notna(row.get('LT_LIT')) else '')
            ms_lit = str(row.get('MS_LIT', '') if pd.notna(row.get('MS_LIT')) else '')

            csv_rows.append(f"{site},{mod_type},{on_function},{on_process},{on_interact},{pmids},{lt_lit},{ms_lit}")

        result[query_item] = "\n".join(csv_rows)
        logging.info(f"Found {len(matches)} regulatory sites for {query_item}")

    return result


@mcp.tool()
def get_disease_sites(query: str, disease: Optional[str] = None, organism: str = 'human') -> Dict:
    """
    Get disease associations for phosphorylation sites.

    Parameters:
    -----------
    query : str
        Mixed - gene names OR 'Gene_Site' format (e.g., 'CTNNB1,BRAF_V600,TP53')
    disease : str, optional
        Filter by disease name (e.g., 'breast cancer', "Alzheimer's disease")
    organism : str
        Organism name (default: 'human')

    Returns:
    --------
    dict
        Dictionary mapping queries to CSV strings with disease association information
        Format: {"CTNNB1": "site,disease,alteration,pmids,lt_lit,ms_lit,notes\n..."}
    """
    if disease_data is None:
        return {"error": "Disease-associated sites data not loaded"}

    # Parse query
    queries = [q.strip() for q in query.split(',')]
    logging.info(f"Querying {len(queries)} items: {queries}")

    # Filter by organism
    organism_lower = organism.lower()
    df = disease_data[disease_data['ORGANISM'].str.lower() == organism_lower]

    # Filter by disease if specified
    if disease:
        df = df[df['DISEASE'].str.contains(disease, case=False, na=False)]
        logging.info(f"Filtering by disease: {disease}")

    result = {}
    for query_item in queries:
        if '_' in query_item:
            # Parse as Gene_Site
            parts = query_item.split('_')
            gene = parts[0]
            site = '_'.join(parts[1:])

            # Find matching entries (site might be in MOD_RSD column)
            matches = df[(df['GENE'].str.upper() == gene.upper()) &
                        (df['MOD_RSD'].str.contains(site, case=False, na=False))]
        else:
            # Parse as gene name - get all sites
            gene = query_item
            matches = df[df['GENE'].str.upper() == gene.upper()]

        if len(matches) == 0:
            logging.warning(f"No disease-associated sites found for {query_item}")
            result[query_item] = "site,disease,alteration,pmids,lt_lit,ms_lit,notes"
            continue

        # Build CSV
        csv_rows = ["site,disease,alteration,pmids,lt_lit,ms_lit,notes"]

        for _, row in matches.iterrows():
            site = str(row.get('MOD_RSD', ''))
            disease_name = str(row.get('DISEASE', '')).replace(',', ';')
            alteration = str(row.get('ALTERATION', '')).replace(',', ';')
            pmids = str(row.get('PMIDs', '') if pd.notna(row.get('PMIDs')) else '').replace(',', ';')
            lt_lit = str(row.get('LT_LIT', '') if pd.notna(row.get('LT_LIT')) else '')
            ms_lit = str(row.get('MS_LIT', '') if pd.notna(row.get('MS_LIT')) else '')
            notes = str(row.get('NOTES', '') if pd.notna(row.get('NOTES')) else '').replace(',', ';')

            csv_rows.append(f"{site},{disease_name},{alteration},{pmids},{lt_lit},{ms_lit},{notes}")

        result[query_item] = "\n".join(csv_rows)
        logging.info(f"Found {len(matches)} disease-associated sites for {query_item}")

    return result


@mcp.tool()
def find_upstream_kinases(query: str, organism: str = 'human') -> Dict:
    """
    Find kinases that phosphorylate specified substrate sites.

    Parameters:
    -----------
    query : str
        Mixed - 'Gene_Site' OR gene names (e.g., 'AKT1_S473,GSK3B_S9,MAPK1')
    organism : str
        Organism name (default: 'human')

    Returns:
    --------
    dict
        Dictionary mapping queries to CSV strings with upstream kinase information
        Format for specific sites: {"AKT1_S473": "kinase,kinase_gene,site_sequence,in_vivo,in_vitro,pmids\n..."}
        Format for genes: {"MAPK1": "site,kinase,kinase_gene,site_sequence,in_vivo,in_vitro,pmids\n..."}
    """
    if kinase_substrate_data is None:
        return {"error": "Kinase-substrate data not loaded"}

    # Parse query
    queries = [q.strip() for q in query.split(',')]
    logging.info(f"Querying {len(queries)} items: {queries}")

    # Filter by organism
    organism_lower = organism.lower()
    df = kinase_substrate_data[kinase_substrate_data['SUB_ORGANISM'].str.lower() == organism_lower]

    result = {}
    for query_item in queries:
        if '_' in query_item and len(query_item.split('_')[1]) > 0 and query_item.split('_')[1][0] in ['S', 'T', 'Y']:
            # Parse as Gene_Site
            parts = query_item.split('_')
            gene = parts[0]
            site = '_'.join(parts[1:])

            # Find matching entries
            matches = df[(df['SUB_GENE'].str.upper() == gene.upper()) &
                        (df['SUB_MOD_RSD'].str.upper() == site.upper())]

            if len(matches) == 0:
                logging.warning(f"No upstream kinases found for {query_item}")
                result[query_item] = "kinase,kinase_gene,site_sequence,in_vivo,in_vitro,pmids"
                continue

            # Build CSV
            csv_rows = ["kinase,kinase_gene,site_sequence,in_vivo,in_vitro,pmids"]

            for _, row in matches.iterrows():
                kinase = str(row.get('KINASE', ''))
                kinase_gene = str(row.get('KIN_GENE_ID', '') if pd.notna(row.get('KIN_GENE_ID')) else '')
                sequence = str(row.get('SITE_+/-7_AA', ''))
                in_vivo = 'X' if str(row.get('IN_VIVO_RXN', '')).strip() == 'X' else ''
                in_vitro = 'X' if str(row.get('IN_VITRO_RXN', '')).strip() == 'X' else ''
                pmids = str(row.get('PMIDs', '') if pd.notna(row.get('PMIDs')) else '')

                csv_rows.append(f"{kinase},{kinase_gene},{sequence},{in_vivo},{in_vitro},{pmids}")

            result[query_item] = "\n".join(csv_rows)
            logging.info(f"Found {len(matches)} upstream kinases for {query_item}")

        else:
            # Parse as gene name - get all sites
            gene = query_item
            matches = df[df['SUB_GENE'].str.upper() == gene.upper()]

            if len(matches) == 0:
                logging.warning(f"No upstream kinases found for gene {gene}")
                result[query_item] = "site,kinase,kinase_gene,site_sequence,in_vivo,in_vitro,pmids"
                continue

            # Build CSV with site information
            csv_rows = ["site,kinase,kinase_gene,site_sequence,in_vivo,in_vitro,pmids"]

            for _, row in matches.iterrows():
                site = str(row.get('SUB_MOD_RSD', ''))
                kinase = str(row.get('KINASE', ''))
                kinase_gene = str(row.get('KIN_GENE_ID', '') if pd.notna(row.get('KIN_GENE_ID')) else '')
                sequence = str(row.get('SITE_+/-7_AA', ''))
                in_vivo = 'X' if str(row.get('IN_VIVO_RXN', '')).strip() == 'X' else ''
                in_vitro = 'X' if str(row.get('IN_VITRO_RXN', '')).strip() == 'X' else ''
                pmids = str(row.get('PMIDs', '') if pd.notna(row.get('PMIDs')) else '')

                csv_rows.append(f"{site},{kinase},{kinase_gene},{sequence},{in_vivo},{in_vitro},{pmids}")

            result[query_item] = "\n".join(csv_rows)
            logging.info(f"Found {len(matches)} kinase-site pairs for gene {gene}")

    return result


@mcp.tool()
def get_interaction_network(query: str, depth: int = 1, evidence_types: Optional[str] = None) -> Dict:
    """
    Get interaction network for specified proteins/sites.

    Parameters:
    -----------
    query : str
        Mixed - gene names OR 'Gene_Site' format (e.g., 'MTOR_S2448,AKT1')
    depth : int
        Network expansion depth (1=direct, 2=includes secondary connections)
    evidence_types : str, optional
        Comma-separated evidence types to filter (e.g., 'in_vivo,in_vitro')

    Returns:
    --------
    dict
        Dictionary with 'nodes', 'edges', and 'summary' CSV strings
    """
    if kinase_substrate_data is None:
        return {"error": "Kinase-substrate data not loaded"}

    # Parse query
    queries = [q.strip() for q in query.split(',')]
    logging.info(f"Building network for {len(queries)} items: {queries}")

    # Parse evidence types filter
    evidence_filter = None
    if evidence_types:
        evidence_filter = [e.strip().lower() for e in evidence_types.split(',')]

    nodes = []
    edges = []
    node_set = set()

    def add_node(node_id, node_type, gene, site=''):
        """Add a node if not already present."""
        if node_id not in node_set:
            nodes.append({
                'node_id': node_id,
                'type': node_type,
                'gene': gene,
                'site': site
            })
            node_set.add(node_id)

    def add_edges_for_query(query_item, current_depth):
        """Add edges for a query item."""
        if current_depth > depth:
            return

        # Determine if this is a site or gene query
        if '_' in query_item and len(query_item.split('_')[1]) > 0 and query_item.split('_')[1][0] in ['S', 'T', 'Y']:
            # Site query - find upstream kinases
            parts = query_item.split('_')
            gene = parts[0]
            site = '_'.join(parts[1:])

            add_node(query_item, 'phosphosite', gene, site)

            # Find upstream kinases
            matches = kinase_substrate_data[
                (kinase_substrate_data['SUB_GENE'].str.upper() == gene.upper()) &
                (kinase_substrate_data['SUB_MOD_RSD'].str.upper() == site.upper())
            ]

            for _, row in matches.iterrows():
                # Check evidence filter
                if evidence_filter:
                    has_evidence = False
                    if 'in_vivo' in evidence_filter and str(row.get('IN_VIVO_RXN', '')).strip() == 'X':
                        has_evidence = True
                    if 'in_vitro' in evidence_filter and str(row.get('IN_VITRO_RXN', '')).strip() == 'X':
                        has_evidence = True
                    if not has_evidence:
                        continue

                kinase = str(row.get('KINASE', ''))
                kinase_id = f"{kinase}"

                add_node(kinase_id, 'kinase', kinase)

                # Build evidence list
                evidence = []
                if str(row.get('IN_VIVO_RXN', '')).strip() == 'X':
                    evidence.append('in_vivo')
                if str(row.get('IN_VITRO_RXN', '')).strip() == 'X':
                    evidence.append('in_vitro')

                edges.append({
                    'source': kinase_id,
                    'target': query_item,
                    'type': 'phosphorylates',
                    'effect': '',
                    'evidence': ';'.join(evidence),
                    'pmids': str(row.get('PMIDs', '') if pd.notna(row.get('PMIDs')) else ''),
                    'sequence': str(row.get('SITE_+/-7_AA', ''))
                })

                # Recurse for depth > 1
                if current_depth < depth:
                    add_edges_for_query(kinase_id, current_depth + 1)
        else:
            # Gene/kinase query - find downstream substrates
            gene = query_item
            add_node(query_item, 'kinase', gene)

            # Find downstream substrates
            matches = kinase_substrate_data[kinase_substrate_data['KINASE'].str.upper() == gene.upper()]

            for _, row in matches.iterrows():
                # Check evidence filter
                if evidence_filter:
                    has_evidence = False
                    if 'in_vivo' in evidence_filter and str(row.get('IN_VIVO_RXN', '')).strip() == 'X':
                        has_evidence = True
                    if 'in_vitro' in evidence_filter and str(row.get('IN_VITRO_RXN', '')).strip() == 'X':
                        has_evidence = True
                    if not has_evidence:
                        continue

                substrate_gene = str(row.get('SUB_GENE', ''))
                substrate_site = str(row.get('SUB_MOD_RSD', ''))
                substrate_id = f"{substrate_gene}_{substrate_site}"

                add_node(substrate_id, 'phosphosite', substrate_gene, substrate_site)

                # Build evidence list
                evidence = []
                if str(row.get('IN_VIVO_RXN', '')).strip() == 'X':
                    evidence.append('in_vivo')
                if str(row.get('IN_VITRO_RXN', '')).strip() == 'X':
                    evidence.append('in_vitro')

                edges.append({
                    'source': query_item,
                    'target': substrate_id,
                    'type': 'phosphorylates',
                    'effect': '',
                    'evidence': ';'.join(evidence),
                    'pmids': str(row.get('PMIDs', '') if pd.notna(row.get('PMIDs')) else ''),
                    'sequence': str(row.get('SITE_+/-7_AA', ''))
                })

                # Recurse for depth > 1
                if current_depth < depth:
                    add_edges_for_query(substrate_id, current_depth + 1)

    # Build network
    for query_item in queries:
        add_edges_for_query(query_item, 1)

    # Convert to CSV format
    nodes_csv = ["node_id,type,gene,site"]
    for node in nodes:
        nodes_csv.append(f"{node['node_id']},{node['type']},{node['gene']},{node['site']}")

    edges_csv = ["source,target,interaction_type,effect,evidence,pmids,sequence"]
    for edge in edges:
        edges_csv.append(f"{edge['source']},{edge['target']},{edge['type']},{edge['effect']},{edge['evidence']},{edge['pmids']},{edge['sequence']}")

    # Count evidence types
    evidence_counts = {}
    total_pmids = set()
    for edge in edges:
        for ev in edge['evidence'].split(';'):
            if ev:
                evidence_counts[ev] = evidence_counts.get(ev, 0) + 1
        if edge['pmids']:
            total_pmids.update(edge['pmids'].split(';'))

    summary_parts = [f"nodes={len(nodes)}", f"edges={len(edges)}"]
    for ev, count in evidence_counts.items():
        summary_parts.append(f"evidence_{ev}={count}")
    summary_parts.append(f"total_pmids={len(total_pmids)}")

    logging.info(f"Network built: {len(nodes)} nodes, {len(edges)} edges")

    return {
        'nodes': '\n'.join(nodes_csv),
        'edges': '\n'.join(edges_csv),
        'summary': ','.join(summary_parts)
    }


@mcp.tool()
def get_pathway_context(query: str, include_downstream: bool = True, include_upstream: bool = True) -> Dict:
    """
    Get pathway context for a single protein or phosphosite.

    Parameters:
    -----------
    query : str
        Single gene OR 'Gene_Site' (e.g., 'AKT1_S473' or 'MTOR')
    include_downstream : bool
        Include downstream substrates (default: True)
    include_upstream : bool
        Include upstream kinases (default: True)

    Returns:
    --------
    dict
        Dictionary with nested CSV structure containing upstream, downstream, functional_outcomes, and summary
    """
    if kinase_substrate_data is None or regulatory_data is None:
        return {"error": "Required datasets not loaded"}

    query = query.strip()
    logging.info(f"Getting pathway context for {query}")

    result = {'query': query}

    # Determine if this is a site or gene query
    is_site = '_' in query and len(query.split('_')[1]) > 0 and query.split('_')[1][0] in ['S', 'T', 'Y']

    if is_site:
        parts = query.split('_')
        gene = parts[0]
        site = '_'.join(parts[1:])
    else:
        gene = query
        site = None

    # Get upstream kinases
    if include_upstream:
        if is_site:
            upstream_matches = kinase_substrate_data[
                (kinase_substrate_data['SUB_GENE'].str.upper() == gene.upper()) &
                (kinase_substrate_data['SUB_MOD_RSD'].str.upper() == site.upper())
            ]
        else:
            # For a kinase, we look for what phosphorylates it
            upstream_matches = kinase_substrate_data[
                kinase_substrate_data['SUB_GENE'].str.upper() == gene.upper()
            ]

        upstream_csv = ["kinase,kinase_gene,evidence,pmids"]
        for _, row in upstream_matches.iterrows():
            kinase = str(row.get('KINASE', ''))
            kinase_gene = str(row.get('KINASE', ''))  # Using KINASE as gene name
            evidence = []
            if str(row.get('IN_VIVO_RXN', '')).strip() == 'X':
                evidence.append('in_vivo')
            if str(row.get('IN_VITRO_RXN', '')).strip() == 'X':
                evidence.append('in_vitro')
            pmids = str(row.get('PMIDs', '') if pd.notna(row.get('PMIDs')) else '')

            upstream_csv.append(f"{kinase},{kinase_gene},{';'.join(evidence)},{pmids}")

        result['upstream'] = '\n'.join(upstream_csv)

    # Get downstream substrates
    if include_downstream:
        downstream_matches = kinase_substrate_data[
            kinase_substrate_data['KINASE'].str.upper() == gene.upper()
        ]

        downstream_csv = ["substrate,gene,site,effect,evidence,pmids"]
        for _, row in downstream_matches.iterrows():
            substrate = str(row.get('SUBSTRATE', ''))
            sub_gene = str(row.get('SUB_GENE', ''))
            sub_site = str(row.get('SUB_MOD_RSD', ''))
            effect = ''  # PSP doesn't have effect in kinase-substrate table
            evidence = []
            if str(row.get('IN_VIVO_RXN', '')).strip() == 'X':
                evidence.append('in_vivo')
            if str(row.get('IN_VITRO_RXN', '')).strip() == 'X':
                evidence.append('in_vitro')
            pmids = str(row.get('PMIDs', '') if pd.notna(row.get('PMIDs')) else '')

            downstream_csv.append(f"{substrate},{sub_gene},{sub_site},{effect},{';'.join(evidence)},{pmids}")

        result['downstream'] = '\n'.join(downstream_csv)

    # Get functional outcomes from regulatory sites
    if is_site:
        regulatory_matches = regulatory_data[
            (regulatory_data['GENE'].str.upper() == gene.upper()) &
            (regulatory_data['MOD_RSD'].str.upper() == site.upper())
        ]
    else:
        regulatory_matches = regulatory_data[
            regulatory_data['GENE'].str.upper() == gene.upper()
        ]

    functional_csv = ["category,description,pmid_count"]

    # Aggregate functional categories
    function_counts = {}
    process_counts = {}

    for _, row in regulatory_matches.iterrows():
        on_function = str(row.get('ON_FUNCTION', '') if pd.notna(row.get('ON_FUNCTION')) else '')
        on_process = str(row.get('ON_PROCESS', '') if pd.notna(row.get('ON_PROCESS')) else '')
        pmids = str(row.get('PMIDs', '') if pd.notna(row.get('PMIDs')) else '')
        pmid_count = len([p for p in pmids.split(';') if p.strip()]) if pmids else 0

        if on_function:
            for func in on_function.split(';'):
                func = func.strip()
                if func:
                    if func not in function_counts:
                        function_counts[func] = 0
                    function_counts[func] += pmid_count

        if on_process:
            for proc in on_process.split(';'):
                proc = proc.strip()
                if proc:
                    if proc not in process_counts:
                        process_counts[proc] = 0
                    process_counts[proc] += pmid_count

    # Add to CSV
    for func, count in sorted(function_counts.items(), key=lambda x: x[1], reverse=True):
        functional_csv.append(f"function,{func},{count}")
    for proc, count in sorted(process_counts.items(), key=lambda x: x[1], reverse=True):
        functional_csv.append(f"process,{proc},{count}")

    result['functional_outcomes'] = '\n'.join(functional_csv)

    # Generate summary
    upstream_count = len(upstream_matches) if include_upstream else 0
    downstream_count = len(downstream_matches) if include_downstream else 0
    functional_count = len(function_counts) + len(process_counts)

    # Count total PMIDs
    all_pmids = set()
    if include_upstream:
        for _, row in upstream_matches.iterrows():
            pmids = str(row.get('PMIDs', '') if pd.notna(row.get('PMIDs')) else '')
            if pmids:
                all_pmids.update(pmids.split(';'))
    if include_downstream:
        for _, row in downstream_matches.iterrows():
            pmids = str(row.get('PMIDs', '') if pd.notna(row.get('PMIDs')) else '')
            if pmids:
                all_pmids.update(pmids.split(';'))
    for _, row in regulatory_matches.iterrows():
        pmids = str(row.get('PMIDs', '') if pd.notna(row.get('PMIDs')) else '')
        if pmids:
            all_pmids.update(pmids.split(';'))

    result['summary'] = f"upstream_kinases={upstream_count},downstream_substrates={downstream_count},functional_outcomes={functional_count},total_pmids={len(all_pmids)}"

    logging.info(f"Pathway context: {upstream_count} upstream, {downstream_count} downstream, {functional_count} functions")

    return result


@mcp.tool()
def get_evidence_summary(query: str, group_by: str = 'site') -> Dict:
    """
    Get aggregated evidence summary for a single protein or phosphosite.

    Parameters:
    -----------
    query : str
        Single gene OR 'Gene_Site' (e.g., 'MTOR_S2448' or 'TP53')
    group_by : str
        Grouping method: 'site', 'kinase', 'disease', or 'function' (default: 'site')

    Returns:
    --------
    dict
        Dictionary with summary statistics and evidence aggregations
    """
    if phospho_data is None or kinase_substrate_data is None or regulatory_data is None or disease_data is None:
        return {"error": "Required datasets not loaded"}

    query = query.strip()
    logging.info(f"Getting evidence summary for {query}")

    # Determine if this is a site or gene query
    is_site = '_' in query and len(query.split('_')[1]) > 0 and query.split('_')[1][0] in ['S', 'T', 'Y']

    if is_site:
        parts = query.split('_')
        gene = parts[0]
        site = '_'.join(parts[1:])
    else:
        gene = query
        site = None

    result = {query: {}}

    # Get all relevant data
    if is_site:
        phospho_matches = phospho_data[
            (phospho_data['GENE'].str.upper() == gene.upper()) &
            (phospho_data['MOD_RSD'].str.upper() == site.upper())
        ]
        kinase_matches = kinase_substrate_data[
            (kinase_substrate_data['SUB_GENE'].str.upper() == gene.upper()) &
            (kinase_substrate_data['SUB_MOD_RSD'].str.upper() == site.upper())
        ]
        regulatory_matches = regulatory_data[
            (regulatory_data['GENE'].str.upper() == gene.upper()) &
            (regulatory_data['MOD_RSD'].str.upper() == site.upper())
        ]
        disease_matches = disease_data[
            (disease_data['GENE'].str.upper() == gene.upper()) &
            (disease_data['MOD_RSD'].str.contains(site, case=False, na=False))
        ]
    else:
        phospho_matches = phospho_data[phospho_data['GENE'].str.upper() == gene.upper()]
        kinase_matches = kinase_substrate_data[kinase_substrate_data['SUB_GENE'].str.upper() == gene.upper()]
        regulatory_matches = regulatory_data[regulatory_data['GENE'].str.upper() == gene.upper()]
        disease_matches = disease_data[disease_data['GENE'].str.upper() == gene.upper()]

    # Calculate evidence strengths
    total_pubs = set()
    in_vivo_count = 0
    in_vitro_count = 0
    ms_count = 0

    # From phospho data
    for _, row in phospho_matches.iterrows():
        lt_lit = row.get('LT_LIT', 0)
        ms_lit = row.get('MS_LIT', 0)
        if pd.notna(lt_lit):
            total_pubs.add(f"phospho_{row.name}")
        if pd.notna(ms_lit) and ms_lit > 0:
            ms_count += int(ms_lit)

    # From kinase-substrate data
    for _, row in kinase_matches.iterrows():
        if str(row.get('IN_VIVO_RXN', '')).strip() == 'X':
            in_vivo_count += 1
        if str(row.get('IN_VITRO_RXN', '')).strip() == 'X':
            in_vitro_count += 1

    # From regulatory data
    for _, row in regulatory_matches.iterrows():
        pmids = str(row.get('PMIDs', '') if pd.notna(row.get('PMIDs')) else '')
        if pmids:
            total_pubs.update(pmids.split(';'))
        ms_lit = row.get('MS_LIT', 0)
        if pd.notna(ms_lit) and ms_lit > 0:
            ms_count += int(ms_lit)

    # Summary statistics
    summary = f"total_publications={len(total_pubs)},in_vivo={in_vivo_count},in_vitro={in_vitro_count},mass_spec={ms_count}"
    result[query]['summary'] = summary

    # Top kinases
    kinase_counts = {}
    for _, row in kinase_matches.iterrows():
        kinase = str(row.get('KINASE', ''))
        if kinase not in kinase_counts:
            kinase_counts[kinase] = {
                'count': 0,
                'in_vivo': False,
                'in_vitro': False,
                'pmids': set()
            }
        kinase_counts[kinase]['count'] += 1
        if str(row.get('IN_VIVO_RXN', '')).strip() == 'X':
            kinase_counts[kinase]['in_vivo'] = True
        if str(row.get('IN_VITRO_RXN', '')).strip() == 'X':
            kinase_counts[kinase]['in_vitro'] = True
        pmids = str(row.get('PMIDs', '') if pd.notna(row.get('PMIDs')) else '')
        if pmids:
            kinase_counts[kinase]['pmids'].update(pmids.split(';'))

    top_kinases_csv = ["kinase,evidence_count,in_vivo,in_vitro,pmid_count"]
    for kinase, data in sorted(kinase_counts.items(), key=lambda x: x[1]['count'], reverse=True)[:10]:
        top_kinases_csv.append(
            f"{kinase},{data['count']},{'yes' if data['in_vivo'] else 'no'},{'yes' if data['in_vitro'] else 'no'},{len(data['pmids'])}"
        )
    result[query]['top_kinases'] = '\n'.join(top_kinases_csv)

    # Functional categories
    func_counts = {}
    for _, row in regulatory_matches.iterrows():
        on_function = str(row.get('ON_FUNCTION', '') if pd.notna(row.get('ON_FUNCTION')) else '')
        if on_function:
            for func in on_function.split(';'):
                func = func.strip()
                if func:
                    func_counts[func] = func_counts.get(func, 0) + 1

    total_func = sum(func_counts.values())
    func_csv = ["category,citation_count,percent"]
    for func, count in sorted(func_counts.items(), key=lambda x: x[1], reverse=True)[:10]:
        percent = (count / total_func * 100) if total_func > 0 else 0
        func_csv.append(f"{func},{count},{percent:.1f}")
    result[query]['functional_categories'] = '\n'.join(func_csv)

    # Disease associations
    disease_counts = {}
    for _, row in disease_matches.iterrows():
        disease = str(row.get('DISEASE', ''))
        alteration = str(row.get('ALTERATION', ''))
        if disease not in disease_counts:
            disease_counts[disease] = {'count': 0, 'alteration': alteration}
        disease_counts[disease]['count'] += 1

    disease_csv = ["disease,citation_count,alteration_type"]
    for disease, data in sorted(disease_counts.items(), key=lambda x: x[1]['count'], reverse=True)[:10]:
        disease_csv.append(f"{disease},{data['count']},{data['alteration']}")
    result[query]['disease_associations'] = '\n'.join(disease_csv)

    # Key PMIDs (top cited)
    all_pmids_list = list(total_pubs)[:10]
    result[query]['key_pmids'] = ','.join(all_pmids_list)

    logging.info(f"Evidence summary: {len(total_pubs)} publications, {len(kinase_counts)} kinases")

    return result


def main():
    logging.basicConfig(level=logging.INFO)
    logging.info("Starting PSP Query Tool")
    mcp.run(transport='stdio')


if __name__ == "__main__":
    main()
