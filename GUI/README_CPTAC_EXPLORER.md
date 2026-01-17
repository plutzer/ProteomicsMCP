# CPTAC Proteomics Explorer - Usage Guide

## Quick Start

### Installation Requirements
Ensure you have the required Python packages:
```bash
pip install shiny shinywidgets plotly pandas numpy scipy cptac
```

### Running the App
```bash
cd GUI
python cptac_explorer_app.py
```

The app will be available at **http://localhost:3838**

## Application Features

### Tab 1: Phospho Tumor vs Normal
Query phosphoproteomics data and visualize differential expression between tumor and normal samples.

**Inputs:**
- **Cancer Type**: Select from 9 CPTAC cancer types (brca, coad, hnscc, luad, ovarian, ccrcc, gbm, lscc, pdac)
- **Query**: Enter phosphosites or gene names (comma-separated)
  - Specific sites: `AKT1_S473,TP53_S15,MTOR_S2448`
  - All sites for genes: `AKT1,TP53,MTOR`
- **Protein-Normalized**: Toggle to use protein-normalized phospho data (recommended: ON)

**Outputs:**
- **Volcano Plot**: Shows log2 fold change vs -log10(p-value), colored by significance (FDR < 0.05)
- **Bar Chart**: Top 15 phosphosites by absolute fold change
- **Box Plot**: Tumor vs normal distribution comparison (top 10 sites)
- **Results Table**: Complete data with statistics
- **CSV Download**: Export all results

### Tab 2: Protein Tumor vs Normal
Query proteomics data for differential protein expression.

**Inputs:**
- **Cancer Type**: Select cancer type
- **Query**: Enter gene names (comma-separated)
  - Example: `AKT1,TP53,EGFR,MTOR,KRAS`

**Outputs:**
- **Volcano Plot**: Log2 fold change vs -log10(p-value) with gene labels
- **Bar Chart**: All queried proteins sorted by fold change
- **Results Table**: Statistics for each protein
- **CSV Download**: Export results

### Tab 3: Correlation Analysis
Analyze correlations between proteins and/or phosphosites across tumor samples.

**Inputs:**
- **Cancer Type**: Select cancer type
- **Query**: Enter items to correlate (comma-separated)
  - Phosphosites: `AKT1_S473,TP53_S15,MTOR_S2448`
  - Proteins: `AKT1,TP53,EGFR` (use with data_type='proteomics')
  - Mixed: `MTOR_S2448,TSC2_protein,AKT1` (use with data_type='both')
- **Data Type**:
  - Phosphoproteomics: Search phospho data only
  - Proteomics: Search protein data only
  - Both (Mixed): Allow mixed phospho/protein queries (use '_protein' suffix to force protein)
- **Protein-Normalized**: For phospho data normalization

**Outputs:**
- **Correlation Heatmap**: Interactive heatmap showing Pearson correlations
- **P-value Heatmap**: FDR-corrected p-values (-log10 transformed)
- **Correlation Matrix Table**: Numerical correlation values
- **CSV Download**: Export correlation matrix

**Special Notes:**
- If querying a gene with many phosphosites, the app will return a list of available sites
- Correlation analysis uses only tumor samples (normal samples excluded)
- Maximum 25 items per correlation query

## Example Queries

### Phospho Tab Examples
1. **PI3K/AKT pathway sites:**
   ```
   AKT1_S473,AKT1_T308,GSK3B_S9,TSC2_S939,MTOR_S2448
   ```

2. **All phosphosites for specific genes:**
   ```
   TP53,EGFR,MAPK1
   ```

### Protein Tab Examples
1. **Core signaling proteins:**
   ```
   AKT1,MTOR,TP53,EGFR,KRAS,PIK3CA
   ```

2. **Cell cycle proteins:**
   ```
   CDK1,CDK2,CCND1,CCNE1,RB1,E2F1
   ```

### Correlation Tab Examples
1. **Phosphosite co-regulation (data_type='phospho'):**
   ```
   MTOR_S2448,MTOR_S2481,RPS6KB1_T389,EIF4EBP1_T37
   ```

2. **Protein-protein correlation (data_type='proteomics'):**
   ```
   MTOR,TSC2,RHEB,RPTOR
   ```

3. **Mixed analysis (data_type='both'):**
   ```
   MTOR_S2448,TSC2_protein,AKT1_S473,PTEN_protein
   ```

## Data Interpretation

### Statistical Metrics
- **log2_fold_change**: Log2 of (mean tumor / mean normal)
  - Positive = upregulated in tumor
  - Negative = downregulated in tumor
- **p_value**: Paired t-test p-value
- **p_value_adjusted**: FDR-corrected p-value (Benjamini-Hochberg)
- **n_pairs**: Number of paired tumor-normal samples
- **mean_tumor/mean_normal**: Average log2 expression values

### Significance Thresholds
- **FDR < 0.05**: Commonly used threshold for significance (highlighted in red in volcano plots)
- **|log2FC| > 1**: Indicates 2-fold change (often biologically meaningful)

## Debugging

The app includes extensive print statements for debugging:
- `[INIT]` - Initialization messages
- `[SERVER]` - Server startup
- `[PHOSPHO]`, `[PROTEIN]`, `[CORRELATION]` - Query execution
- `[PARSE_*]` - Data parsing
- `[*_VOLCANO]`, `[*_BAR]`, `[*_HEATMAP]` - Visualization rendering
- `[*_ERROR]`, `[*_EXCEPTION]` - Error messages

Check the terminal/console output for detailed debugging information.

## Data Sources

All data is sourced from:
- **CPTAC (Clinical Proteomic Tumor Analysis Consortium)**
- **Baylor College of Medicine (BCM)** proteomics pipeline
- Covers 9 cancer types with paired tumor-normal samples

## Tips

1. **Start small**: Test with 2-3 items first before querying many genes
2. **Use normalized data**: Protein-normalized phospho data accounts for total protein abundance changes
3. **Check sample sizes**: Higher n_pairs = more reliable statistics
4. **Explore correlations**: Correlation analysis reveals co-regulation patterns
5. **Export data**: Download CSV files for further analysis in R, Python, or Excel

## Troubleshooting

**Issue: Query returns no results**
- Check gene/site spelling (case-sensitive)
- Try querying just the gene name first to see available sites
- Not all genes have phospho data in all cancer types

**Issue: "Too many sites" warning in correlation tab**
- Query specific sites instead of gene names
- The app shows available sites for that gene

**Issue: Slow performance**
- First query per cancer type loads data from CPTAC (can take 30-60 seconds)
- Subsequent queries for the same cancer type are faster
- Correlation analysis with many items (>15) can be slow

**Issue: App won't start**
- Check that all dependencies are installed
- Verify `cptac_proteomics.py` is in the parent directory
- Check that port 3838 is not already in use

## Contact

For issues or questions about the CPTAC data explorer, please check:
- CPTAC documentation: https://github.com/PayneLab/cptac
- Shiny for Python: https://shiny.posit.co/py/
