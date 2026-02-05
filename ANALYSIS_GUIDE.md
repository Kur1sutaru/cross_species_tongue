# GSEA and KEGG Pathway Analysis Guide
## For Conserved Markers Across Cell Subtypes

## Overview
This guide provides complete instructions and scripts for performing Gene Set Enrichment Analysis (GSEA) and KEGG pathway analysis on your conserved markers data.

## Your Data Structure
Based on your screenshot, you have:
- **9 cell subtypes** (endothelial cell populations)
- **Conserved markers** across 4 species (rat, marmoset, mouse, human)
- Each marker has:
  - p-values and adjusted p-values
  - Average log2 fold changes
  - Percentage of cells expressing the marker

## Cell Subtypes in Your Data:
1. Metabolically-active Capillary ECs (56 genes)
2. Capillary ECs (301 genes)
3. ECM-remodeling fibroblasts (184 genes)
4. Smooth muscle cell (262 genes)
5. Lymphatic endothelial (149 genes)
6. Angiogenic ECs (195 genes)
7. Pericytes-to-SMC transition (115 genes)
8. ECM-producing fibroblasts (303 genes)
9. Inflammatory fibroblasts (186 genes)

## Analysis Approach

### 1. **Gene Set Enrichment Analysis (GSEA)**
- Uses **ranked gene lists** based on log2FC
- Identifies pathways enriched at top/bottom of list
- No arbitrary cutoff needed
- More powerful for detecting subtle changes

### 2. **Over-Representation Analysis (ORA)**
- Uses **top significant genes** (e.g., log2FC > 0.5, p < 0.05)
- Tests if genes are enriched in specific pathways
- Requires binary classification (significant vs not)

### 3. **KEGG Pathway Analysis**
- Focuses on metabolic and signaling pathways
- Well-curated, high-quality annotations
- Organism-specific pathways

### 4. **Gene Ontology (GO) Analysis**
- Biological Process (BP): biological objectives
- Molecular Function (MF): molecular activities
- Cellular Component (CC): cellular locations

## Installation Instructions

### Option 1: R-based Analysis (Recommended)

```R
# Install required packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "clusterProfiler",
    "org.Hs.eg.db",  # Human
    "org.Mm.eg.db",  # Mouse  
    "org.Rn.eg.db",  # Rat
    "enrichplot",
    "DOSE"
))

install.packages(c(
    "ggplot2",
    "dplyr",
    "openxlsx",
    "tidyverse"
))
```

### Option 2: Python-based Analysis

```bash
pip install gseapy pandas numpy matplotlib seaborn openpyxl rpy2 pyreadr
```

## Running the Analysis

### Using R (gsea_kegg_analysis.R):

```R
# Load the script
source("gsea_kegg_analysis.R")

# The script will:
# 1. Read your RDS file
# 2. Process each subtype
# 3. Perform GSEA and ORA
# 4. Generate plots
# 5. Export to Excel
```

### Using Python (gsea_kegg_analysis_python.py):

```python
# Install gseapy
pip install gseapy

# Run the script
python gsea_kegg_analysis_python.py
```

## Output Files

### 1. **Excel File**: `enrichment_results_all_subtypes.xlsx`
Each sheet contains:
- Pathway/Term name
- NES (Normalized Enrichment Score)
- p-value and FDR q-value
- Leading edge genes
- Gene set size

### 2. **PDF Plots** (in `enrichment_plots/` directory):
For each subtype:
- `{Subtype}_GSEA_GO.pdf` - GO term enrichment
- `{Subtype}_GSEA_KEGG.pdf` - KEGG pathway enrichment
- `{Subtype}_ORA_GO.pdf` - GO over-representation
- `{Subtype}_ORA_KEGG.pdf` - KEGG over-representation

### 3. **RDS File**: `enrichment_results.rds`
R object for further analysis

### 4. **Summary**: `enrichment_summary.csv`
Overview of all results

## Interpreting Results

### GSEA Results:
- **NES (Normalized Enrichment Score)**:
  - Positive: pathway enriched in upregulated genes
  - Negative: pathway enriched in downregulated genes
  - |NES| > 1.5: strong enrichment
  
- **FDR q-value < 0.05**: statistically significant
- **FDR q-value < 0.25**: potentially interesting

### ORA Results:
- **Combined Score**: enrichment metric (higher = more significant)
- **Adjusted P-value < 0.05**: statistically significant
- **Overlap**: number of your genes in the pathway

## Example Interpretation

For "Metabolically-active Capillary ECs":
- Top markers: FLT1, PECAM1, FCGRT, TCF4
- Expected pathways:
  - Angiogenesis
  - Vascular development
  - Endothelial cell migration
  - VEGF signaling
  - Blood vessel morphogenesis

## Customization Options

### Change Species:
```R
# In gsea_kegg_analysis.R, line 162:
enrichment_results[[subtype]] <- perform_enrichment(
    conserved_markers[[subtype]], 
    subtype,
    species = "mouse"  # Change to "rat", "marm", or "human"
)
```

### Adjust Significance Thresholds:
```R
# In the gseGO/gseKEGG functions:
pvalueCutoff = 0.05,  # Adjust this value
minGSSize = 10,       # Minimum gene set size
maxGSSize = 500       # Maximum gene set size
```

### Select Specific Pathways:
```R
# Use specific gene sets
gene_sets <- c("hsa04010", "hsa04060", "hsa04514")  # Specific KEGG pathways
```

## Troubleshooting

### Issue: "No significant pathways found"
**Solutions:**
1. Lower the p-value cutoff: `pvalueCutoff = 0.1`
2. Use more genes in ORA
3. Check if gene symbols are correct
4. Try different species databases

### Issue: "Cannot convert gene symbols"
**Solutions:**
1. Check gene symbol format (should be uppercase for human)
2. Update annotation databases:
   ```R
   BiocManager::install("org.Hs.eg.db", force = TRUE)
   ```

### Issue: "Too few genes for analysis"
**Solutions:**
1. Use less stringent filtering
2. Combine related subtypes
3. Use at least 50-100 genes for meaningful results

## Advanced Options

### Compare Multiple Subtypes:
```R
# Compare enrichment across subtypes
library(clusterProfiler)

# Compare two subtypes
comp_result <- compareCluster(
    geneCluster = list(
        Subtype1 = genes1,
        Subtype2 = genes2
    ),
    fun = "enrichKEGG",
    organism = "hsa"
)

dotplot(comp_result)
```

### Custom Gene Sets:
```R
# Use custom GMT file
custom_pathways <- read.gmt("custom_pathways.gmt")

gsea_custom <- GSEA(
    gene_list_entrez,
    TERM2GENE = custom_pathways
)
```

## Citation

If you use these analyses in your publication, please cite:

1. **clusterProfiler**:
   Yu et al. (2012). clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS: A Journal of Integrative Biology, 16(5), 284-287.

2. **KEGG**:
   Kanehisa et al. (2021). KEGG: integrating viruses and cellular organisms. Nucleic Acids Research, 49(D1), D545-D551.

3. **Gene Ontology**:
   Gene Ontology Consortium (2021). The Gene Ontology resource: enriching a GOld mine. Nucleic Acids Research, 49(D1), D325-D334.

## Next Steps

1. Upload your RDS file: `conserved_markers_all_4_species.rds`
2. Run the analysis script
3. Examine the Excel file for detailed results
4. Review the plots in `enrichment_plots/`
5. Interpret the pathways in the context of your biology

## Questions?

Common questions:
- **Q: Which species should I use?**
  A: Use human for best pathway annotations. Mouse is also well-annotated.

- **Q: Should I use GSEA or ORA?**
  A: Use both! GSEA is more powerful, ORA is more interpretable.

- **Q: How many pathways should I report?**
  A: Focus on top 10-20 most significant pathways per subtype.

- **Q: What if results don't make biological sense?**
  A: Check your gene list, filtering criteria, and consider cell type-specific biology.

## Contact

For issues with the scripts or analysis questions, provide:
1. Your RDS file structure
2. Error messages
3. Number of genes per subtype
4. Expected biological pathways

Good luck with your analysis!
