# ProteomicsMCP

MCP Server for querying CPTAC proteomics data through the Model Context Protocol.

## Installation

1. **Install Python dependencies:**

```bash
pip install cptac pandas numpy scipy mcp
```

2. **Clone or download this repository to your local machine.**

3. **Configure Claude Desktop to use this MCP server:**

Edit your Claude Desktop configuration file:
- **Windows**: `%APPDATA%\Claude\claude_desktop_config.json`
- **macOS**: `~/Library/Application Support/Claude/claude_desktop_config.json`
- **Linux**: `~/.config/Claude/claude_desktop_config.json`

Add the following configuration (replace paths with your actual absolute paths):

```json
{
  "mcpServers": {
    "cptac": {
      "command": "/absolute/path/to/python",
      "args": ["/absolute/path/to/ProteomicsMCP/cptac_proteomics.py"],
      "cwd": "/absolute/path/to/ProteomicsMCP"
    }
  }
}
```

**Important:** You must use absolute paths for both the Python executable and the script.

### Finding Your Python Path

**Using a virtual environment (recommended):**
```bash
# Create and activate a virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install cptac pandas numpy scipy mcp

# Find the Python executable path
which python  # On Windows: where python
```

**Using conda:**
```bash
conda create -n proteomics python=3.11
conda activate proteomics
pip install cptac pandas numpy scipy mcp
which python  # On Windows: where python
```

**Example configurations:**

Windows with conda:
```json
{
  "mcpServers": {
    "cptac": {
      "command": "C:\\Users\\YourName\\anaconda3\\envs\\proteomics\\python.exe",
      "args": ["C:\\Users\\YourName\\Repos\\ProteomicsMCP\\cptac_proteomics.py"],
      "cwd": "C:\\Users\\YourName\\Repos\\ProteomicsMCP"
    }
  }
}
```

macOS/Linux with venv:
```json
{
  "mcpServers": {
    "cptac": {
      "command": "/home/yourname/ProteomicsMCP/venv/bin/python",
      "args": ["/home/yourname/ProteomicsMCP/cptac_proteomics.py"],
      "cwd": "/home/yourname/ProteomicsMCP"
    }
  }
}
```

4. **Restart Claude Desktop** for the changes to take effect.

## Usage

Once configured, the CPTAC MCP server provides the following tools in Claude Desktop:

- **`get_cancer_types()`**: List available cancer datasets
- **`phospho_tumor_vs_normal(cancer, query, normalized)`**: Query phosphoproteomics data with tumor vs normal statistics
  - Returns log2 fold change, p-value, and FDR-adjusted p-value for each phosphosite
  - Performs paired t-tests between tumor and normal samples
- **`protein_tumor_vs_normal(cancer, query)`**: Query whole-cell proteomics data with tumor vs normal statistics
  - Returns log2 fold change, p-value, and FDR-adjusted p-value for each protein
  - Performs paired t-tests between tumor and normal samples

### Example Queries

Ask Claude:
- "What cancer types are available in CPTAC?"
- "Get phosphoproteomics data for AKT1_S473 in breast cancer"
- "Show me all phosphosites for TP53 in lung adenocarcinoma, normalized"
- "Compare AKT1 and EGFR protein levels between tumor and normal in ovarian cancer"
- "What are the most significantly changed phosphosites for MAPK1 in colorectal cancer?"

## Supported Cancer Types

- `brca` - Breast Cancer
- `coad` - Colon Adenocarcinoma
- `hnscc` - Head and Neck Squamous Cell Carcinoma
- `luad` - Lung Adenocarcinoma
- `ovarian` - Ovarian Cancer
- `ccrcc` - Clear Cell Renal Cell Carcinoma
- `gbm` - Glioblastoma
- `lscc` - Lung Squamous Cell Carcinoma
- `pdac` - Pancreatic Ductal Adenocarcinoma