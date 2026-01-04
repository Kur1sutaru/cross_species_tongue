# cross_species_tongue

# Cross-Species 1:1 Ortholog Mapping and Harmony Integration

This repository provides a fully reproducible workflow for performing **cross-species single-cell RNA-seq integration** using **strict 1:1 orthologs** and **Harmony**.  
The workflow supports:

- Mouse → Human  
- Rat → Human  
- Marmoset → Human  
- Human (reference)

The pipeline uses **local ortholog tables** exported from **Ensembl BioMart**, ensuring compatibility with any R installation (including R 4.4) without requiring Bioconductor packages such as `biomaRt`, `AnnotationHub`, or `orthogene`.

---

## 🧬 Overview

### Workflow summary

1. Retrieve ortholog tables from Ensembl BioMart (web interface)  
2. Filter to **strict 1:1 orthologs**  
3. Save as TSV files  
4. Validate tables using the included script  
5. Convert each species’ Seurat object to human gene symbols  
6. Merge all species  
7. Integrate using **Harmony**  
8. Perform clustering, UMAP, and QC summaries  

---

## 📥 1. Retrieve 1:1 Ortholog Tables from Ensembl BioMart

### Step-by-step instructions

1. Open BioMart:  
   https://www.ensembl.org/biomart/martview

2. Select:
   - **Database:** Ensembl Genes (latest)
   - **Dataset:**  
     - *Mus musculus genes (GRCm39)*  
     - *Rattus norvegicus genes (mRatBN7.2)*  
     - *Callithrix jacchus genes (ASM275486v1)*  

3. (Optional) Filter to protein-coding genes:  
   **Filters → GENE → Gene type → protein_coding**

4. Select attributes:

   **GENE**
   - Gene stable ID  
   - Gene name  

   **HOMOLOGS → Human orthologs**
   - Human gene stable ID  
   - Human gene name  
   - Human orthology type  
   - Human orthology confidence  

5. Click **Results** → choose **TSV** → **Go**

6. Save files as:
   - `mouse_human.tsv`
   - `rat_human.tsv`
   - `marmoset_human.tsv`

---

## 🧪 2. Clean and Filter to Strict 1:1 Orthologs

Example in R:

```r
df <- read.table("mouse_human.tsv", sep="\t", header=TRUE)

df_1to1 <- df[df$Human.orthology.type == "ortholog_one2one", ]

df_1to1 <- df_1to1[, c("Gene.name", "Human.gene.name")]
colnames(df_1to1) <- c("gene_from", "gene_human")

write.table(df_1to1, "mouse_human_1to1.tsv", sep="\t", quote=FALSE, row.names=FALSE)
