# Cross-Species Tongue Single-Cell Atlas

A fully reproducible pipeline for **cross-species single-cell RNA-seq integration and conservation analysis** of tongue tissue across four mammalian species: **Human, Mouse, Rat, and Marmoset** (~135,736 cells total).

This workflow covers ortholog mapping, Harmony integration, cell type annotation, conservation scoring, functional enrichment, and publication-quality figure generation.

OSF repository available at <https://doi.org/10.17605/OSF.IO/4M98Q>

---

## 📁 Repository Structure

```
cross_species_tongue/
├── README.md                    # This file — full workflow documentation
├── ANALYSIS_GUIDE.md            # Quick-reference guide and parameter notes
├── manuscript_writing_plan.md   # Manuscript outline and figure plan
├── ortologs_convert.R           # Step 1–2: Ortholog filtering and gene symbol conversion
├── harmony.R                    # Step 3–5: Seurat merging, Harmony integration, clustering
└── gsea_kegg_analysis.R         # Step 6: GSEA and KEGG pathway enrichment
```

---

## 🧬 Biological Context

The tongue is a complex organ innervated by multiple cranial nerves and populated by diverse cell types — including epithelial/keratinocyte layers, immune cells, fibroblasts, endothelial cells, smooth and skeletal muscle, and Schwann cells. This cross-species atlas aims to:

- Identify **conserved and species-specific transcriptomic signatures** across mammalian tongues
- Characterize **cell type diversity** with a unified annotation framework
- Provide a comparative resource for studying **pain-relevant cell types** and tongue-associated pathologies (e.g., oral squamous cell carcinoma)

**Species included:**

| Species | Abbreviation | Color |
|---------|-------------|-------|
| *Homo sapiens* | Human | `#E41A1C` |
| *Mus musculus* | Mouse | `#4DAF4A` |
| *Callithrix jacchus* | Marmoset | `#377EB8` / `#FF7F00` |
| *Rattus norvegicus* | Rat | `#984EA3` |

**Cell type ordering (canonical):** Epithelial → Immune → Endothelial → Fibroblast → Muscle → Schwann

---

## ⚙️ Software Requirements

| Package | Version | Purpose |
|---------|---------|---------|
| R | ≥ 4.4 | Core runtime |
| Seurat | v5 | scRNA-seq object handling |
| harmony | latest | Cross-species batch correction |
| clusterProfiler | latest | GSEA / KEGG enrichment |
| dplyr, tidyr | latest | Data wrangling |
| ggplot2, patchwork | latest | Visualization |

> **Note:** This pipeline deliberately avoids `biomaRt`, `AnnotationHub`, and `orthogene` to ensure compatibility across R installations. All ortholog mappings use **local TSV tables from Ensembl BioMart**.

---

## 🔬 Workflow Overview

```
[1] Retrieve ortholog tables from Ensembl BioMart (web interface)
        ↓
[2] Filter to strict 1:1 orthologs → save as TSV
        ↓
[3] Convert each species' Seurat object to Human gene symbols
        ↓
[4] Merge all species into a single Seurat object
        ↓
[5] Harmony integration → Clustering → UMAP
        ↓
[6] Cell type annotation (marker-based, iterative DEG refinement)
        ↓
[7] Conservation analysis (FindConservedMarkers + CV-based scoring)
        ↓
[8] GSEA / KEGG pathway enrichment
        ↓
[9] Publication-quality figure generation
```

---

## Step 1 — Retrieve 1:1 Ortholog Tables from Ensembl BioMart

**Script:** *(manual step — web interface)*

1. Navigate to [https://www.ensembl.org/biomart/martview](https://www.ensembl.org/biomart/martview)
2. For each non-human species, select:
   - **Database:** Ensembl Genes (latest release)
   - **Dataset:**
     - Mouse → *Mus musculus* genes (GRCm39)
     - Rat → *Rattus norvegicus* genes (mRatBN7.2)
     - Marmoset → *Callithrix jacchus* genes (ASM275486v1)
3. *(Optional)* Filter to protein-coding genes: **Filters → GENE → Gene type → protein_coding**
4. Select the following **Attributes**:
   - GENE: Gene stable ID, Gene name
   - HOMOLOGS → Human orthologs: Human gene stable ID, Human gene name, Human orthology type, Human orthology confidence
5. Export as **TSV** and save as:
   - `mouse_human.tsv`
   - `rat_human.tsv`
   - `marmoset_human.tsv`

---

## Step 2 — Filter to Strict 1:1 Orthologs

**Script:** `ortologs_convert.R` (first section)

Raw BioMart tables contain one-to-many and many-to-one ortholog relationships. This step retains only **ortholog_one2one** pairs to ensure unambiguous gene mapping.

```r
df <- read.table("mouse_human.tsv", sep = "\t", header = TRUE)

# Keep only strict 1:1 orthologs
df_1to1 <- df[df$Human.orthology.type == "ortholog_one2one", ]

# Select and rename columns
df_1to1 <- df_1to1[, c("Gene.name", "Human.gene.name")]
colnames(df_1to1) <- c("gene_from", "gene_human")

# Remove rows with empty human gene names
df_1to1 <- df_1to1[df_1to1$gene_human != "", ]

write.table(df_1to1, "mouse_human_1to1.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
```

Repeat for rat and marmoset. Output files:
- `mouse_human_1to1.tsv`
- `rat_human_1to1.tsv`
- `marmoset_human_1to1.tsv`

---

## Step 3 — Convert Gene Symbols to Human Orthologs

**Script:** `ortologs_convert.R` (second section)

Each species' Seurat object has genes named in its own gene symbol space. This step renames features to their human ortholog symbols, retaining only 1:1 mappable genes.

```r
convert_to_human <- function(seurat_obj, ortholog_table, species_name) {
  genes_in_obj  <- rownames(seurat_obj)
  ortho_filt    <- ortholog_table[ortholog_table$gene_from %in% genes_in_obj, ]
  
  # Remove duplicated human gene targets (keep only unique mappings)
  ortho_filt    <- ortho_filt[!duplicated(ortho_filt$gene_human), ]
  ortho_filt    <- ortho_filt[!duplicated(ortho_filt$gene_from), ]
  
  seurat_sub    <- seurat_obj[ortho_filt$gene_from, ]
  rownames(seurat_sub[["RNA"]]) <- ortho_filt$gene_human
  seurat_sub$species <- species_name
  return(seurat_sub)
}

mouse_converted    <- convert_to_human(mouse_seurat,    mouse_ortho,    "Mouse")
rat_converted      <- convert_to_human(rat_seurat,      rat_ortho,      "Rat")
marmoset_converted <- convert_to_human(marmoset_seurat, marmoset_ortho, "Marmoset")
```

> Each converted object retains only genes with a valid 1:1 human ortholog, substantially reducing the feature space. A `species` metadata column is added to each object.

---

## Step 4 — Merge All Species

**Script:** `harmony.R` (merging section)

All four species objects (Human + three converted) are merged into a single Seurat v5 object. No integration is applied at this stage — the merge is purely additive.

```r
merged <- merge(
  human_seurat,
  y = list(mouse_converted, rat_converted, marmoset_converted),
  add.cell.ids = c("Human", "Mouse", "Rat", "Marmoset"),
  project = "CrossSpeciesTongue"
)

# Standard preprocessing on merged object
merged <- NormalizeData(merged)
merged <- FindVariableFeatures(merged, nfeatures = 3000)
merged <- ScaleData(merged)
merged <- RunPCA(merged, npcs = 50)
```

---

## Step 5 — Harmony Integration, Clustering, and UMAP

**Script:** `harmony.R` (integration section)

Harmony corrects for batch effects across species using the `species` metadata variable. After integration, standard Seurat clustering and UMAP dimensionality reduction are performed.

```r
library(harmony)

merged <- RunHarmony(merged, group.by.vars = "species", dims.use = 1:30)

merged <- FindNeighbors(merged, reduction = "harmony", dims = 1:30)
merged <- FindClusters(merged, resolution = 0.5)   # Adjust resolution as needed
merged <- RunUMAP(merged, reduction = "harmony", dims = 1:30)

DimPlot(merged, group.by = "species", cols = c(
  "Human"    = "#E41A1C",
  "Mouse"    = "#4DAF4A",
  "Marmoset" = "#377EB8",
  "Rat"      = "#984EA3"
))
```

**Key parameters to tune:**
- `dims.use` in `RunHarmony()`: typically 1:30 or 1:40
- `resolution` in `FindClusters()`: lower → fewer clusters; higher → more granular

---

## Step 6 — Cell Type Annotation

**Script:** `harmony.R` (annotation section)

Clusters are annotated using canonical markers for each major cell class. Annotation is performed iteratively: initial broad annotation → subclustering → DEG-based refinement.

### Major cell type markers

| Cell Type | Key Markers |
|-----------|-------------|
| Epithelial / Keratinocyte | `KRT5`, `KRT14`, `KRT13`, `TP63`, `EPCAM` |
| Immune (broad) | `PTPRC` (CD45) |
| — Macrophage / LAM | `MRC1`, `TREM2`, `SPP1`, `C1QA` |
| — T cell | `CD3E`, `CD3D` |
| — γδT17 | `TRDC`, `TRGC1`, `IL17A` |
| — NK / Cytotoxic ILC | `GNLY`, `NKG7`, `GZMB` |
| — B cell | `MS4A1`, `CD79A` |
| — Mast cell | `KIT`, `TPSAB1` |
| Endothelial | `PECAM1`, `CDH5`, `VWF` |
| Fibroblast | `DCN`, `LUM`, `COL1A1` |
| Smooth Muscle | `ACTA2`, `MYH11` |
| Skeletal Muscle | `MYH2`, `TNNT3`, `TTN` |
| Schwann | `MPZ`, `MBP`, `S100B` |

### Annotation workflow

```r
# Visualize marker expression
FeaturePlot(merged, features = c("KRT5", "PTPRC", "PECAM1", "DCN"))
DotPlot(merged, features = marker_list, group.by = "seurat_clusters") + RotatedAxis()

# Assign identities
cluster_annotations <- c(
  "0" = "Keratinocyte_Basal",
  "1" = "Fibroblast",
  # ... etc
)
merged$celltype <- plyr::mapvalues(
  merged$seurat_clusters,
  from = names(cluster_annotations),
  to   = cluster_annotations
)
```

> Immune subtypes (13–15 subtypes) require additional subclustering. See `harmony.R` for the full iterative DEG-based refinement loop.

---

## Step 7 — Conservation Analysis

**Script:** `harmony.R` (conservation section)

This module identifies genes that are conserved markers across species for each cell type, assigns a **stability score** using the coefficient of variation (CV), and classifies genes into conservation quartiles.

### `FindConservedMarkers` per cell type

```r
Idents(merged) <- "celltype"
DefaultAssay(merged) <- "RNA"

conserved_markers <- FindConservedMarkers(
  merged,
  ident.1       = "Keratinocyte_Basal",
  grouping.var  = "species",
  only.pos      = TRUE,
  logfc.threshold = 0.25
)
```

### CV-based gene stability classification

For each gene, the average log2FC across species is computed, and the CV (SD / mean) is used to classify stability:

```r
# Compute per-species avg_logFC columns, then:
logfc_cols <- grep("_avg_log2FC$", colnames(conserved_markers), value = TRUE)
fc_matrix  <- conserved_markers[, logfc_cols]

conserved_markers$mean_logFC <- rowMeans(fc_matrix, na.rm = TRUE)
conserved_markers$sd_logFC   <- apply(fc_matrix, 1, sd, na.rm = TRUE)
conserved_markers$cv         <- conserved_markers$sd_logFC / conserved_markers$mean_logFC

# Percentile-based quartile classification
q_breaks <- quantile(conserved_markers$cv, probs = c(0, 0.25, 0.50, 0.75, 1.0), na.rm = TRUE)
conserved_markers$conservation_class <- cut(
  conserved_markers$cv,
  breaks = q_breaks,
  labels = c("Highly Conserved", "Moderately Conserved", "Variable", "Highly Variable"),
  include.lowest = TRUE
)
```

**Interpretation:** Genes in Q1 (lowest CV) are the most consistently expressed across all four species → **Highly Conserved**.

---



---



---

## 🗂 Data Availability

Raw scRNA-seq data for each species were processed independently using standard Seurat v5 pipelines (QC filtering, normalization, dimensional reduction) prior to integration. Processed Seurat RDS objects are available upon request.

---

## 📖 Citation

If you use this pipeline, please cite:

> Villalba Silva GC, *et al.* (in preparation). A cross-species single-cell atlas of the mammalian tongue reveals conserved and divergent transcriptomic programs. *[Journal TBD]*

---

## 🤝 Contact

**Cristal Villalba, Msc PhD**  
Postdoctoral Fellow — Center for Pain Therapeutics and Addiction Research (CPTAR)  
UT Health San Antonio, School of Dentistry  
📧 [villalbasilva@uthscsa.edu]

---

## 📄 License

This repository is provided for reproducibility and transparency in scientific research. Code is available under the MIT License.
