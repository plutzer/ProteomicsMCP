"""
CPTAC Proteomics Data Explorer

Interactive Shiny web application for exploring CPTAC proteomics and phosphoproteomics data.
Provides three main query interfaces:
1. Phospho Tumor vs Normal - Query phosphoproteomics data
2. Protein Tumor vs Normal - Query proteomics data
3. Correlation Analysis - Analyze correlations between proteins/phosphosites

Author: Claude Code
Date: 2026-01-11
"""

from shiny import App, Inputs, Outputs, Session, reactive, render, ui
from shinywidgets import output_widget, render_plotly
import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
import numpy as np
import os
import sys
import tempfile
import datetime
from io import StringIO

# Add parent directory to path to import modules
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from cptac_proteomics import (
    correlation_analysis  # Still needed for correlation tab
)
from cptac_backend import CancerDataBackend

print("[INIT] CPTAC Proteomics Explorer starting...")
print("[INIT] Initializing backend...")
print("[INIT] This will take 5-10 minutes for all cancers (or 1 min in test mode)")
print("[INIT] Loading and preprocessing all tumor vs normal statistics...")

# Initialize backend - set test_mode=True for faster development
backend = CancerDataBackend(test_mode=True)

print("[INIT] Backend ready! Queries will now be instant.")

# ==================== Helper Functions ====================

def get_cancer_choices():
    """Get list of available cancer types"""
    return ["brca", "coad", "hnscc", "luad", "ovarian", "ccrcc", "gbm", "lscc", "pdac"]


def parse_phospho_results(result_dict):
    """
    Parse phospho_tumor_vs_normal result into single DataFrame

    Args:
        result_dict: Dict with 'data' key containing gene -> CSV string mapping

    Returns:
        pd.DataFrame with all phosphosites combined
    """
    print("[PARSE_PHOSPHO] Parsing phosphoproteomics results...")
    all_dfs = []

    for gene, csv_string in result_dict['data'].items():
        df = pd.read_csv(StringIO(csv_string))
        df['gene'] = gene
        all_dfs.append(df)
        print(f"[PARSE_PHOSPHO] Parsed {len(df)} sites for gene {gene}")

    combined = pd.concat(all_dfs, ignore_index=True)

    # Reorder columns for better display
    cols = ['gene', 'site', 'peptide', 'database_id', 'log2_fold_change',
            'p_value', 'p_value_adjusted', 'n_pairs', 'mean_tumor', 'mean_normal']

    print(f"[PARSE_PHOSPHO] Combined {len(combined)} total phosphosites")
    return combined[cols]


# ==================== UI Definition ====================

app_ui = ui.page_navbar(
    # ===== Tab 1: Phospho Tumor vs Normal =====
    ui.nav_panel(
        "Phospho Tumor vs Normal",
        ui.layout_columns(
            # Left column: Input controls
            ui.card(
                ui.card_header("Query Parameters"),
                ui.input_select(
                    "phospho_cancer",
                    "Cancer Type",
                    choices=get_cancer_choices(),
                    selected="pdac"
                ),
                ui.input_selectize(
                    "phospho_query",
                    "Select Phosphosites",
                    choices=[],  # Populated dynamically when cancer changes
                    multiple=True,
                    options={
                        'placeholder': 'Type to search (e.g., AKT1_S473)...',
                        'maxItems': 50,
                        'plugins': ['remove_button']
                    }
                ),
                ui.input_switch(
                    "phospho_normalized",
                    "Protein-Normalized",
                    value=True
                ),
                ui.input_action_button(
                    "run_phospho_query",
                    "Run Query",
                    class_="btn-primary"
                ),
                ui.download_button(
                    "download_phospho_csv",
                    "Download CSV"
                )
            ),
            col_widths=[4, 8]
        ),
        # Visualization rows
        ui.layout_columns(
            ui.card(
                ui.card_header("Volcano Plot"),
                output_widget("phospho_volcano")
            ),
            ui.card(
                ui.card_header("Tumor vs Normal Intensities"),
                output_widget("phospho_barplot")
            ),
            col_widths=[6, 6]
        ),
        ui.card(
            ui.card_header("Results Table"),
            ui.output_data_frame("phospho_table")
        )
    ),

    # ===== Tab 2: Protein Tumor vs Normal =====
    ui.nav_panel(
        "Protein Tumor vs Normal",
        ui.layout_columns(
            # Left column: Input controls
            ui.card(
                ui.card_header("Query Parameters"),
                ui.input_select(
                    "protein_cancer",
                    "Cancer Type",
                    choices=get_cancer_choices(),
                    selected="pdac"
                ),
                ui.input_selectize(
                    "protein_query",
                    "Select Proteins",
                    choices=[],  # Populated dynamically when cancer changes
                    multiple=True,
                    options={
                        'placeholder': 'Type to search (e.g., AKT1)...',
                        'maxItems': 50,
                        'plugins': ['remove_button']
                    }
                ),
                ui.input_action_button(
                    "run_protein_query",
                    "Run Query",
                    class_="btn-primary"
                ),
                ui.download_button(
                    "download_protein_csv",
                    "Download CSV"
                )
            ),
            col_widths=[4, 8]
        ),
        # Visualizations
        ui.layout_columns(
            ui.card(
                ui.card_header("Volcano Plot"),
                output_widget("protein_volcano")
            ),
            ui.card(
                ui.card_header("Bar Chart"),
                output_widget("protein_barplot")
            ),
            col_widths=[6, 6]
        ),
        ui.card(
            ui.card_header("Results Table"),
            ui.output_data_frame("protein_table")
        )
    ),

    # ===== Tab 3: Correlation Analysis =====
    ui.nav_panel(
        "Correlation Analysis",
        ui.layout_columns(
            # Left column: Input controls
            ui.card(
                ui.card_header("Query Parameters"),
                ui.input_select(
                    "corr_cancer",
                    "Cancer Type",
                    choices=get_cancer_choices(),
                    selected="pdac"
                ),
                ui.input_text_area(
                    "corr_query",
                    "Items to Correlate (comma-separated)",
                    placeholder="Examples:\nAKT1_S473,TP53_S15 (phosphosites)\nAKT1,TP53 (proteins)\nMTOR_S2448,TSC2_protein (mixed)",
                    rows=5
                ),
                ui.input_radio_buttons(
                    "corr_data_type",
                    "Data Type",
                    choices={
                        "phospho": "Phosphoproteomics",
                        "proteomics": "Proteomics",
                        "both": "Both (Mixed)"
                    },
                    selected="phospho"
                ),
                ui.input_switch(
                    "corr_normalized",
                    "Protein-Normalized (for phospho)",
                    value=True
                ),
                ui.input_action_button(
                    "run_corr_query",
                    "Run Query",
                    class_="btn-primary"
                ),
                ui.download_button(
                    "download_corr_csv",
                    "Download Correlation Matrix"
                )
            ),
            col_widths=[4, 8]
        ),
        # Heatmaps
        ui.layout_columns(
            ui.card(
                ui.card_header("Correlation Heatmap"),
                output_widget("corr_heatmap")
            ),
            ui.card(
                ui.card_header("P-value Heatmap (FDR-corrected)"),
                output_widget("pvalue_heatmap")
            ),
            col_widths=[6, 6]
        ),
        ui.card(
            ui.card_header("Correlation Matrix"),
            ui.output_data_frame("corr_table")
        )
    ),

    title="CPTAC Proteomics Explorer",
    id="main_nav"
)


# ==================== Server Logic ====================

def server(input: Inputs, output: Outputs, session: Session):
    print("[SERVER] CPTAC Proteomics Explorer server started")

    # ===== Reactive State =====
    phospho_results = reactive.Value(pd.DataFrame())
    protein_results = reactive.Value(pd.DataFrame())
    corr_results = reactive.Value({})

    # ===== Dynamic Choice Updaters =====

    @reactive.effect
    @reactive.event(input.phospho_cancer)
    def update_phospho_choices():
        """Update phosphosite choices when cancer type changes"""
        cancer = input.phospho_cancer()
        print(f"[PHOSPHO_CHOICES] Loading choices for {cancer}")

        choices = backend.get_phospho_choices(cancer)
        print(f"[PHOSPHO_CHOICES] {len(choices)} phosphosites available")

        ui.update_selectize(
            "phospho_query",
            choices=choices,
            selected=[]  # Clear selection when cancer changes
        )

    @reactive.effect
    @reactive.event(input.protein_cancer)
    def update_protein_choices():
        """Update protein choices when cancer type changes"""
        cancer = input.protein_cancer()
        print(f"[PROTEIN_CHOICES] Loading choices for {cancer}")

        choices = backend.get_protein_choices(cancer)
        print(f"[PROTEIN_CHOICES] {len(choices)} proteins available")

        ui.update_selectize(
            "protein_query",
            choices=choices,
            selected=[]  # Clear selection when cancer changes
        )

    # ========================================
    # PHOSPHO TAB - Server Logic
    # ========================================

    @reactive.effect
    @reactive.event(input.run_phospho_query)
    def execute_phospho_query():
        """Execute phosphoproteomics query using backend"""
        # Get inputs
        cancer = input.phospho_cancer()
        selected_sites = input.phospho_query()  # Returns list from selectize
        normalized = input.phospho_normalized()

        print(f"[PHOSPHO] Cancer: {cancer}, Selected: {len(selected_sites)} sites, Normalized: {normalized}")

        # Validate input
        if not selected_sites:
            ui.notification_show("Please select phosphosites", type="error")
            print("[PHOSPHO ERROR] No sites selected")
            return

        try:
            # Query backend (instant!)
            df = backend.query_phospho(cancer, selected_sites, normalized)
            phospho_results.set(df)

            print(f"[PHOSPHO SUCCESS] Retrieved {len(df)} phosphosites")
            ui.notification_show(f"Retrieved {len(df)} phosphosites", type="message", duration=2)

        except Exception as e:
            ui.notification_show(f"Error: {str(e)}", type="error")
            print(f"[PHOSPHO EXCEPTION] {str(e)}")
            import traceback
            traceback.print_exc()

    @render_plotly
    def phospho_volcano():
        """Render layered volcano plot: background + query"""
        query_df = phospho_results.get()
        cancer = input.phospho_cancer()
        normalized = input.phospho_normalized()

        # Get background data (all sites)
        bg_df = backend.get_phospho_background(cancer, normalized)

        print(f"[PHOSPHO_VOLCANO] Background: {len(bg_df)} sites, Query: {len(query_df)} sites")

        fig = go.Figure()

        # Layer 1: Background (grey, non-interactive)
        if not bg_df.empty:
            fig.add_trace(go.Scatter(
                x=bg_df['log2_fold_change'],
                y=-np.log10(bg_df['p_value']),
                mode='markers',
                marker=dict(color='lightgrey', size=4, opacity=0.3),
                hoverinfo='skip',
                showlegend=True,
                name=f'All sites (n={len(bg_df)})'
            ))

        # Layer 2: Query points (red, interactive)
        if not query_df.empty:
            # Create labels from MultiIndex
            labels = [f"{g}_{s}" for g, s in query_df.index]
            query_df = query_df.reset_index()

            hover_text = (
                query_df['gene'] + '_' + query_df['site'] +
                '<br>log2FC: ' + query_df['log2_fold_change'].round(3).astype(str) +
                '<br>p-value: ' + query_df['p_value'].apply(lambda x: f'{x:.4e}') +
                '<br>FDR: ' + query_df['p_value_adjusted'].apply(lambda x: f'{x:.4e}')
            )

            fig.add_trace(go.Scatter(
                x=query_df['log2_fold_change'],
                y=-np.log10(query_df['p_value']),
                mode='markers+text',
                marker=dict(color='red', size=8, opacity=1.0),
                text=labels,
                textposition='top center',
                textfont=dict(size=9),
                hovertext=hover_text,
                hoverinfo='text',
                showlegend=True,
                name=f'Selected (n={len(query_df)})'
            ))

        # Add p=0.05 threshold line
        fig.add_hline(y=-np.log10(0.05), line_dash="dash", line_color="gray",
                     annotation_text="p=0.05", annotation_position="right")

        fig.update_layout(
            title='Volcano Plot: Phosphoproteomics',
            xaxis_title='Log2 Fold Change (Tumor - Normal)',
            yaxis_title='-log10(p-value)',
            hovermode='closest'
        )

        return fig

    @render_plotly
    def phospho_barplot():
        """Render grouped bar plot: tumor vs normal intensities"""
        df = phospho_results.get()

        if df.empty:
            return go.Figure().update_layout(title="No data")

        # Create labels from MultiIndex
        labels = [f"{g}_{s}" for g, s in df.index]
        plot_df = df.reset_index()

        print(f"[PHOSPHO_BAR] Rendering {len(plot_df)} sites")

        fig = go.Figure()

        # Tumor bars (red)
        fig.add_trace(go.Bar(
            x=labels,
            y=plot_df['mean_tumor'],
            name='Tumor',
            marker_color='#e74c3c',
            text=plot_df['mean_tumor'].round(2),
            textposition='outside'
        ))

        # Normal bars (blue)
        fig.add_trace(go.Bar(
            x=labels,
            y=plot_df['mean_normal'],
            name='Normal',
            marker_color='#3498db',
            text=plot_df['mean_normal'].round(2),
            textposition='outside'
        ))

        fig.update_layout(
            barmode='group',
            title='Tumor vs Normal Intensities',
            xaxis_title='Phosphosite',
            yaxis_title='Log2 Intensity (Mean)',
            xaxis={'tickangle': -45}
        )

        return fig

    @render.data_frame
    def phospho_table():
        """Render phospho results table"""
        df = phospho_results.get()

        if df.empty:
            return pd.DataFrame()

        # Format numeric columns
        display_df = df.copy()

        for col in ['log2_fold_change', 'mean_tumor', 'mean_normal']:
            if col in display_df.columns:
                display_df[col] = display_df[col].apply(
                    lambda x: f'{x:.3f}' if pd.notna(x) else ''
                )

        for col in ['p_value', 'p_value_adjusted']:
            if col in display_df.columns:
                display_df[col] = display_df[col].apply(
                    lambda x: f'{x:.4e}' if pd.notna(x) else ''
                )

        print(f"[PHOSPHO_TABLE] Displaying {len(display_df)} rows")
        return render.DataGrid(display_df, width="100%", height="400px")

    @render.download(filename=lambda: f"phospho_results_{input.phospho_cancer()}_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.csv")
    def download_phospho_csv():
        """Download phospho results as CSV"""
        df = phospho_results.get()

        if df.empty:
            print("[PHOSPHO_DOWNLOAD] No data to download")
            return None

        path = os.path.join(tempfile.gettempdir(), "phospho_results.csv")
        df.to_csv(path, index=False)
        print(f"[PHOSPHO_DOWNLOAD] Saved {len(df)} rows to {path}")
        return path

    # ========================================
    # PROTEIN TAB - Server Logic
    # ========================================

    @reactive.effect
    @reactive.event(input.run_protein_query)
    def execute_protein_query():
        """Execute proteomics query using backend"""
        cancer = input.protein_cancer()
        selected_genes = input.protein_query()  # Returns list from selectize

        print(f"[PROTEIN] Cancer: {cancer}, Selected: {len(selected_genes)} genes")

        # Validate input
        if not selected_genes:
            ui.notification_show("Please select proteins", type="error")
            return

        try:
            # Query backend (instant!)
            df = backend.query_protein(cancer, selected_genes)
            protein_results.set(df)

            print(f"[PROTEIN SUCCESS] Retrieved {len(df)} proteins")
            ui.notification_show(f"Retrieved {len(df)} proteins", type="message")

        except Exception as e:
            ui.notification_show(f"Error: {str(e)}", type="error")
            print(f"[PROTEIN EXCEPTION] {str(e)}")
            import traceback
            traceback.print_exc()

    @render_plotly
    def protein_volcano():
        """Render layered volcano plot: background + query"""
        query_df = protein_results.get()
        cancer = input.protein_cancer()

        # Get background data (all proteins)
        bg_df = backend.get_protein_background(cancer)

        print(f"[PROTEIN_VOLCANO] Background: {len(bg_df)} proteins, Query: {len(query_df)} proteins")

        fig = go.Figure()

        # Layer 1: Background (grey, non-interactive)
        if not bg_df.empty:
            fig.add_trace(go.Scatter(
                x=bg_df['log2_fold_change'],
                y=-np.log10(bg_df['p_value']),
                mode='markers',
                marker=dict(color='lightgrey', size=4, opacity=0.3),
                hoverinfo='skip',
                showlegend=True,
                name=f'All proteins (n={len(bg_df)})'
            ))

        # Layer 2: Query points (red, interactive)
        if not query_df.empty:
            query_df = query_df.reset_index()

            hover_text = (
                query_df['gene'] +
                '<br>log2FC: ' + query_df['log2_fold_change'].round(3).astype(str) +
                '<br>p-value: ' + query_df['p_value'].apply(lambda x: f'{x:.4e}') +
                '<br>FDR: ' + query_df['p_value_adjusted'].apply(lambda x: f'{x:.4e}')
            )

            fig.add_trace(go.Scatter(
                x=query_df['log2_fold_change'],
                y=-np.log10(query_df['p_value']),
                mode='markers+text',
                marker=dict(color='red', size=8, opacity=1.0),
                text=query_df['gene'],
                textposition='top center',
                hovertext=hover_text,
                hoverinfo='text',
                showlegend=True,
                name=f'Selected (n={len(query_df)})'
            ))

        # Add p=0.05 threshold line
        fig.add_hline(y=-np.log10(0.05), line_dash="dash", line_color="gray",
                     annotation_text="p=0.05", annotation_position="right")

        fig.update_layout(
            title='Volcano Plot: Proteomics',
            xaxis_title='Log2 Fold Change (Tumor - Normal)',
            yaxis_title='-log10(p-value)',
            hovermode='closest'
        )

        return fig

    @render_plotly
    def protein_barplot():
        """Render grouped bar plot: tumor vs normal intensities"""
        df = protein_results.get()

        if df.empty:
            return go.Figure().update_layout(title="No data")

        # Get gene names from index
        genes = df.index.tolist()
        plot_df = df.reset_index()

        print(f"[PROTEIN_BAR] Rendering {len(plot_df)} proteins")

        fig = go.Figure()

        # Tumor bars (red)
        fig.add_trace(go.Bar(
            x=genes,
            y=plot_df['mean_tumor'],
            name='Tumor',
            marker_color='#e74c3c',
            text=plot_df['mean_tumor'].round(2),
            textposition='outside'
        ))

        # Normal bars (blue)
        fig.add_trace(go.Bar(
            x=genes,
            y=plot_df['mean_normal'],
            name='Normal',
            marker_color='#3498db',
            text=plot_df['mean_normal'].round(2),
            textposition='outside'
        ))

        fig.update_layout(
            barmode='group',
            title='Tumor vs Normal Intensities',
            xaxis_title='Gene',
            yaxis_title='Log2 Intensity (Mean)',
            xaxis={'tickangle': -45}
        )

        return fig

    @render.data_frame
    def protein_table():
        """Render protein results table"""
        df = protein_results.get()

        if df.empty:
            return pd.DataFrame()

        # Format numeric columns
        display_df = df.copy()

        for col in ['log2_fold_change', 'mean_tumor', 'mean_normal']:
            if col in display_df.columns:
                display_df[col] = display_df[col].apply(
                    lambda x: f'{x:.3f}' if pd.notna(x) else ''
                )

        for col in ['p_value', 'p_value_adjusted']:
            if col in display_df.columns:
                display_df[col] = display_df[col].apply(
                    lambda x: f'{x:.4e}' if pd.notna(x) else ''
                )

        print(f"[PROTEIN_TABLE] Displaying {len(display_df)} rows")
        return render.DataGrid(display_df, width="100%", height="400px")

    @render.download(filename=lambda: f"protein_results_{input.protein_cancer()}_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.csv")
    def download_protein_csv():
        """Download protein results as CSV"""
        df = protein_results.get()

        if df.empty:
            print("[PROTEIN_DOWNLOAD] No data to download")
            return None

        path = os.path.join(tempfile.gettempdir(), "protein_results.csv")
        df.to_csv(path, index=False)
        print(f"[PROTEIN_DOWNLOAD] Saved {len(df)} rows to {path}")
        return path

    # ========================================
    # CORRELATION TAB - Server Logic
    # ========================================

    @reactive.effect
    @reactive.event(input.run_corr_query)
    def execute_corr_query():
        """Execute correlation analysis query"""
        with ui.Progress(min=0, max=100) as progress:
            progress.set(message="Running correlation analysis...", value=10)

            # Get inputs
            cancer = input.corr_cancer()
            query = input.corr_query().strip()
            data_type = input.corr_data_type()
            normalized = 'true' if input.corr_normalized() else 'false'

            print(f"[CORRELATION] Cancer: {cancer}, Query: {query}, Type: {data_type}, Normalized: {normalized}")

            # Validate input
            if not query:
                ui.notification_show("Please enter items to correlate", type="error")
                print("[CORRELATION ERROR] Empty query")
                return

            progress.set(message="Calling MCP tool...", value=30)

            try:
                # Call MCP tool
                result = correlation_analysis(cancer, query, data_type, normalized)

                progress.set(message="Parsing results...", value=70)

                # Check for errors
                if 'error' in result:
                    ui.notification_show(result['error'], type="error")
                    print(f"[CORRELATION ERROR] {result['error']}")
                    return

                # Check for "too many sites" warning
                if 'message' in result:
                    msg = result['message']
                    if 'available_sites' in result:
                        msg += f"\n\nShowing first 10 of {result['n_sites']} sites:\n"
                        msg += ", ".join(result['available_sites'][:10])
                        if result['n_sites'] > 10:
                            msg += f"\n... and {result['n_sites'] - 10} more"
                    ui.notification_show(msg, type="warning", duration=10)
                    print(f"[CORRELATION WARNING] {msg}")
                    return

                # Store results
                corr_results.set(result)

                progress.set(value=100)
                print(f"[CORRELATION SUCCESS] Matrix for {result['n_samples']} samples")
                ui.notification_show(
                    f"Successfully computed correlation matrix ({result['n_samples']} samples)",
                    type="message",
                    duration=3
                )

            except Exception as e:
                ui.notification_show(f"Error: {str(e)}", type="error")
                print(f"[CORRELATION EXCEPTION] {str(e)}")
                import traceback
                traceback.print_exc()

    @render_plotly
    def corr_heatmap():
        """Render correlation heatmap"""
        result = corr_results.get()

        if not result or 'correlation_matrix' not in result:
            return go.Figure().update_layout(
                title="No data - run a query first",
                xaxis_title="Items",
                yaxis_title="Items"
            )

        print(f"[CORR_HEATMAP] Rendering correlation matrix")

        # Parse CSV matrix
        corr_df = pd.read_csv(StringIO(result['correlation_matrix']), index_col=0)

        # Create heatmap
        fig = go.Figure(data=go.Heatmap(
            z=corr_df.values,
            x=corr_df.columns.tolist(),
            y=corr_df.index.tolist(),
            colorscale='RdBu_r',
            zmid=0,
            zmin=-1,
            zmax=1,
            text=corr_df.values,
            texttemplate='%{text:.2f}',
            textfont={"size": 10},
            colorbar=dict(title="Correlation<br>Coefficient")
        ))

        fig.update_layout(
            title=f"Correlation Matrix (n={result['n_samples']} tumor samples)",
            xaxis={'side': 'bottom', 'tickangle': -45},
            yaxis={'autorange': 'reversed'},
            height=500
        )

        return fig

    @render_plotly
    def pvalue_heatmap():
        """Render p-value heatmap"""
        result = corr_results.get()

        if not result or 'p_value_matrix' not in result:
            return go.Figure().update_layout(
                title="No data",
                xaxis_title="Items",
                yaxis_title="Items"
            )

        print(f"[PVALUE_HEATMAP] Rendering p-value matrix")

        # Parse CSV matrix
        pval_df = pd.read_csv(StringIO(result['p_value_matrix']), index_col=0)

        # Convert to -log10 for better visualization (avoid log(0))
        log_pval = -np.log10(pval_df.values + 1e-300)

        # Create heatmap
        fig = go.Figure(data=go.Heatmap(
            z=log_pval,
            x=pval_df.columns.tolist(),
            y=pval_df.index.tolist(),
            colorscale='Reds',
            text=pval_df.values,
            texttemplate='%{text:.3e}',
            textfont={"size": 9},
            colorbar=dict(title="-log10<br>(p-value)")
        ))

        fig.update_layout(
            title="P-value Matrix (FDR-corrected)",
            xaxis={'side': 'bottom', 'tickangle': -45},
            yaxis={'autorange': 'reversed'},
            height=500
        )

        return fig

    @render.data_frame
    def corr_table():
        """Render correlation matrix as table"""
        result = corr_results.get()

        if not result or 'correlation_matrix' not in result:
            return pd.DataFrame()

        # Parse CSV matrix
        corr_df = pd.read_csv(StringIO(result['correlation_matrix']), index_col=0)

        # Format to 3 decimal places
        display_df = corr_df.applymap(lambda x: f'{x:.3f}' if pd.notna(x) else '')

        print(f"[CORR_TABLE] Displaying {len(corr_df)}x{len(corr_df.columns)} matrix")
        return render.DataGrid(display_df, width="100%", height="400px")

    @render.download(filename=lambda: f"correlation_matrix_{input.corr_cancer()}_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.csv")
    def download_corr_csv():
        """Download correlation matrix as CSV"""
        result = corr_results.get()

        if not result or 'correlation_matrix' not in result:
            print("[CORR_DOWNLOAD] No data to download")
            return None

        # Parse matrix
        corr_df = pd.read_csv(StringIO(result['correlation_matrix']), index_col=0)

        path = os.path.join(tempfile.gettempdir(), "correlation_matrix.csv")
        corr_df.to_csv(path)
        print(f"[CORR_DOWNLOAD] Saved {len(corr_df)}x{len(corr_df.columns)} matrix to {path}")
        return path


# ==================== App Creation ====================

app = App(app_ui, server)

if __name__ == "__main__":
    print("[MAIN] Starting CPTAC Proteomics Explorer...")
    print("[MAIN] App will be available at http://localhost:3838")
    app.run(host="0.0.0.0", port=3838)
